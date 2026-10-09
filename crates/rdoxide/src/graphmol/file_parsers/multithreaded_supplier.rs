//! Molecule suppliers that parse on several threads (upstream
//! `FileParsers/Multithreaded{MolSupplier,SDMolSupplier,SmilesMolSupplier}`).
//!
//! A reader thread cuts the input into records, `num_writer_threads` worker
//! threads turn records into molecules, and the consumer takes the results
//! from a bounded queue. Upstream hands out molecules in whatever order the
//! workers finish; here they are re-ordered by record id, so the sequence of
//! molecules is the same as the file's whatever the thread count.
//!
//! A record that cannot be parsed yields `None` (its text is still
//! available through [`MultithreadedMolSupplier::last_item_text`]).

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, sync_channel};
use std::sync::{Arc, Mutex, OnceLock};
use std::thread::JoinHandle;

use super::mol_file_parser::{MolFileParserParams, mol_from_mol_data_stream};
use super::sd_mol_supplier::{process_mol_property_list, strip};
use super::smiles_mol_supplier::{SmilesMolSupplierParams, tokenize};
use super::utils::{LineStream, lossy};
use crate::graphmol::romol::RWMol;
use crate::graphmol::smiles_parse::mol_from_smiles;
use crate::rdgeneral::threads::num_threads_to_use;
use crate::rdgeneral::{RDProps, common_properties as cp};
use crate::{RDKitError, Result, rd_log_warning};

/// `MultithreadedMolSupplier::Parameters`
#[derive(Debug, Clone)]
pub struct SupplierParams {
    /// Worker threads; zero or negative means all hardware threads minus
    /// `|num_writer_threads|` (at least one).
    pub num_writer_threads: i32,
    /// Capacity of the queue of records waiting to be parsed.
    pub size_input_queue: usize,
    /// Capacity of the queue of parsed molecules waiting to be taken.
    pub size_output_queue: usize,
}

impl Default for SupplierParams {
    fn default() -> Self {
        SupplierParams {
            num_writer_threads: 1,
            size_input_queue: 5,
            size_output_queue: 5,
        }
    }
}

/// One record of the input.
#[derive(Debug, Clone)]
pub struct Record {
    /// 1-based position of the record in the input.
    pub id: u32,
    /// Line number reported for the record.
    pub line: u32,
    pub text: Vec<u8>,
}

/// Cuts the input into records on the reader thread.
trait RecordSource: Send {
    /// Passes the records to `emit` until the input ends or `emit` returns
    /// `false`.
    fn read_records(self: Box<Self>, emit: &mut dyn FnMut(Record) -> bool);
}

/// Turns a record into a molecule on a worker thread; `None` when the record
/// cannot be parsed.
trait RecordParser: Send + Sync {
    fn parse(&self, record: &Record) -> Option<RWMol>;
}

/// What a worker hands to the consumer.
type Parsed = (u32, Option<RWMol>, Vec<u8>);

/// The running pipeline.
struct Pipeline {
    results: Option<Receiver<Parsed>>,
    stop: Arc<AtomicBool>,
    threads: Vec<JoinHandle<()>>,
}

impl Pipeline {
    fn start(
        source: Box<dyn RecordSource>,
        parser: Arc<dyn RecordParser>,
        params: &SupplierParams,
    ) -> Self {
        let stop = Arc::new(AtomicBool::new(false));
        let (record_tx, record_rx) = sync_channel::<Record>(params.size_input_queue.max(1));
        let (result_tx, result_rx) = sync_channel::<Parsed>(params.size_output_queue.max(1));
        let record_rx = Arc::new(Mutex::new(record_rx));

        let mut threads = Vec::new();
        let reader_stop = Arc::clone(&stop);
        threads.push(std::thread::spawn(move || {
            source.read_records(&mut |record| {
                !reader_stop.load(Ordering::Relaxed) && record_tx.send(record).is_ok()
            });
        }));
        for _ in 0..num_threads_to_use(params.num_writer_threads) {
            let (records, results, parser) = (
                Arc::clone(&record_rx),
                result_tx.clone(),
                Arc::clone(&parser),
            );
            let worker_stop = Arc::clone(&stop);
            threads.push(std::thread::spawn(move || {
                loop {
                    let next = records.lock().unwrap_or_else(|e| e.into_inner()).recv();
                    let Ok(record) = next else { break };
                    if worker_stop.load(Ordering::Relaxed) {
                        break;
                    }
                    let mol = parser.parse(&record);
                    if results.send((record.id, mol, record.text)).is_err() {
                        break;
                    }
                }
            }));
        }
        Pipeline {
            results: Some(result_rx),
            stop,
            threads,
        }
    }

    /// Stops the threads and waits for them.
    fn shut_down(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        // dropping the receiver unblocks workers waiting to deliver; their
        // exit drops the record receiver, which unblocks the reader
        self.results = None;
        for thread in self.threads.drain(..) {
            let _ = thread.join();
        }
    }
}

impl Drop for Pipeline {
    fn drop(&mut self) {
        self.shut_down();
    }
}

/// A supplier whose molecules are parsed concurrently.
pub struct MultithreadedMolSupplier {
    /// The pipeline's ingredients until the first molecule is requested.
    pending: Option<(Box<dyn RecordSource>, Arc<dyn RecordParser>)>,
    pipeline: Option<Pipeline>,
    params: SupplierParams,
    /// Parsed records that arrived ahead of their turn.
    buffered: BTreeMap<u32, (Option<RWMol>, Vec<u8>)>,
    next_id: u32,
    last_record_id: u32,
    last_item_text: Vec<u8>,
    /// Set by the SD supplier: whether property lists are applied.
    process_property_lists: Option<Arc<AtomicBool>>,
    closed: bool,
}

impl std::fmt::Debug for MultithreadedMolSupplier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MultithreadedMolSupplier")
            .field("params", &self.params)
            .field("next_id", &self.next_id)
            .finish_non_exhaustive()
    }
}

impl MultithreadedMolSupplier {
    fn new(
        source: Box<dyn RecordSource>,
        parser: Arc<dyn RecordParser>,
        params: SupplierParams,
    ) -> Self {
        MultithreadedMolSupplier {
            pending: Some((source, parser)),
            pipeline: None,
            params,
            buffered: BTreeMap::new(),
            next_id: 1,
            last_record_id: 0,
            last_item_text: Vec::new(),
            process_property_lists: None,
            closed: false,
        }
    }

    /// Supplies the molecules of an SD file.
    pub fn from_sd_file(
        file_name: &str,
        parse_params: MolFileParserParams,
        params: SupplierParams,
    ) -> Result<Self> {
        let data = read_input(file_name)?;
        let process = Arc::new(AtomicBool::new(true));
        let parser = SdParser {
            params: parse_params,
            process_property_lists: Arc::clone(&process),
        };
        let mut supplier = Self::new(Box::new(SdSource { data }), Arc::new(parser), params);
        supplier.process_property_lists = Some(process);
        Ok(supplier)
    }

    /// Supplies the molecules of a delimited SMILES file.
    pub fn from_smiles_file(
        file_name: &str,
        parse_params: SmilesMolSupplierParams,
        params: SupplierParams,
    ) -> Result<Self> {
        let data = read_input(file_name)?;
        let columns = Arc::new(OnceLock::new());
        let source = SmilesSource {
            data,
            delimiter: parse_params.delimiter.clone(),
            title_line: parse_params.title_line,
            columns: Arc::clone(&columns),
        };
        let parser = SmilesParser {
            params: parse_params,
            columns,
        };
        Ok(Self::new(Box::new(source), Arc::new(parser), params))
    }

    /// `getProcessPropertyLists` (SD suppliers only; `true` otherwise).
    pub fn process_property_lists(&self) -> bool {
        self.process_property_lists
            .as_ref()
            .is_none_or(|f| f.load(Ordering::Relaxed))
    }

    /// `setProcessPropertyLists`; ignored by suppliers that have none.
    pub fn set_process_property_lists(&self, value: bool) {
        if let Some(flag) = &self.process_property_lists {
            flag.store(value, Ordering::Relaxed);
        }
    }

    /// The queue of parsed records (starting the threads on first use);
    /// `None` once the supplier is closed.
    fn results(&mut self) -> Option<&Receiver<Parsed>> {
        if self.closed {
            return None;
        }
        if let Some((source, parser)) = self.pending.take() {
            self.pipeline = Some(Pipeline::start(source, parser, &self.params));
        }
        self.pipeline.as_ref().and_then(|p| p.results.as_ref())
    }

    /// Whether every record has been handed out (waits for the next record
    /// to be parsed if necessary). As upstream, this is `false` until the
    /// first molecule has been requested, even for an empty input.
    pub fn at_end(&mut self) -> bool {
        if self.pending.is_some() || !self.buffered.is_empty() {
            return false;
        }
        match self.results().map(Receiver::recv) {
            Some(Ok((id, mol, text))) => {
                self.buffered.insert(id, (mol, text));
                false
            }
            _ => true,
        }
    }

    /// The next molecule in file order: `None` at the end of the input and
    /// `Some(None)` for a record that could not be parsed.
    pub fn next_mol(&mut self) -> Option<Option<RWMol>> {
        loop {
            if let Some((mol, text)) = self.buffered.remove(&self.next_id) {
                self.last_record_id = self.next_id;
                self.last_item_text = text;
                self.next_id += 1;
                return Some(mol);
            }
            match self.results().map(Receiver::recv) {
                Some(Ok((id, mol, text))) => {
                    self.buffered.insert(id, (mol, text));
                }
                // no more results: nothing is left once the buffer is empty
                _ if self.buffered.is_empty() => return None,
                // (cannot happen: records are numbered without gaps)
                _ => self.next_id = *self.buffered.keys().next().expect("not empty"),
            }
        }
    }

    /// The id (1-based position) of the record returned last.
    pub fn last_record_id(&self) -> u32 {
        self.last_record_id
    }

    /// The text of the record returned last.
    pub fn last_item_text(&self) -> String {
        lossy(&self.last_item_text)
    }

    /// Stops the threads; the supplier yields nothing afterwards.
    pub fn close(&mut self) {
        self.closed = true;
        self.pending = None;
        self.buffered.clear();
        if let Some(mut pipeline) = self.pipeline.take() {
            pipeline.shut_down();
        }
    }
}

fn read_input(file_name: &str) -> Result<Vec<u8>> {
    std::fs::read(file_name).map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))
}

// ------------------------------------------------------------------ SD files

struct SdSource {
    data: Vec<u8>,
}

fn starts_with(line: &[u8], prefix: &[u8]) -> bool {
    line.starts_with(prefix)
}

fn is_blank(line: &[u8]) -> bool {
    line.iter()
        .all(|c| matches!(c, b' ' | b'\t' | b'\r' | b'\n'))
}

impl RecordSource for SdSource {
    /// `MultithreadedSDMolSupplier::extractNextRecord`: lines up to the
    /// next `$$$$` that follows a blank line or `M  END`.
    fn read_records(self: Box<Self>, emit: &mut dyn FnMut(Record) -> bool) {
        let mut stream = LineStream::new(&self.data);
        let (mut line_no, mut id) = (0u32, 1u32);
        while !stream.eof() {
            let first_line = line_no;
            let (mut prev, mut current): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
            let mut record = Vec::new();
            while !stream.eof()
                && !stream.fail()
                && ((!is_blank(&prev) && !starts_with(&prev, b"M  END"))
                    || !starts_with(&current, b"$$$$"))
            {
                std::mem::swap(&mut prev, &mut current);
                current = stream.getline_raw();
                record.extend_from_slice(&current);
                record.push(b'\n');
                line_no += 1;
            }
            // ignore trailing new lines
            if record.iter().all(|&c| c == b'\n' || c == b'\r') {
                return;
            }
            let item = Record {
                id,
                line: first_line,
                text: record,
            };
            id += 1;
            if !emit(item) {
                return;
            }
        }
    }
}

struct SdParser {
    params: MolFileParserParams,
    process_property_lists: Arc<AtomicBool>,
}

impl SdParser {
    /// `MultithreadedSDMolSupplier::readMolProps`: the data fields after the
    /// molecule block. It is its own implementation upstream, not the one of
    /// the sequential supplier: for example, a data item without a label
    /// swallows the rest of its record and fails.
    fn read_mol_props(&self, stream: &mut LineStream<'_>, mol: &mut RWMol) -> Result<()> {
        let process_lists = self.process_property_lists.load(Ordering::Relaxed);
        let mut has_prop = false;
        let mut warning_issued = false;
        let mut label = String::new();
        let mut line = stream.getline_raw();
        while !stream.eof() && !stream.fail() && !line.starts_with(b"$$$$") {
            let text = strip(&line).to_vec();
            if !text.is_empty() {
                if text[0] == b'>' {
                    has_prop = true;
                    warning_issued = false;
                    let rest = &text[1..];
                    let open = rest.iter().position(|&c| c == b'<');
                    let close = rest.iter().rposition(|&c| c == b'>');
                    match (open, close) {
                        (Some(open), Some(close)) if close != open + 1 => {
                            // (`substr` with a wrapped length reads to the end)
                            label = lossy(rest.get(open + 1..close).unwrap_or(&rest[open + 1..]));
                            self.read_data_item(stream, mol, &label, process_lists);
                        }
                        _ => {
                            // no usable label: skip the item. The loop of
                            // upstream never refreshes its blank-line test,
                            // so a non-blank first line reads to the end.
                            let first = stream.getline_raw();
                            if !strip(&first).is_empty() {
                                loop {
                                    stream.getline_raw();
                                    if stream.eof() {
                                        return Err(RDKitError::FileParse(
                                            "End of data field name not found".into(),
                                        ));
                                    }
                                }
                            }
                        }
                    }
                } else if self.params.strict_parsing {
                    return Err(RDKitError::FileParse(
                        "Problems encountered parsing data fields".into(),
                    ));
                } else if !warning_issued {
                    if has_prop {
                        rd_log_warning!(
                            "Property <{label}> will be truncated after the first blank line\n"
                        );
                    } else {
                        rd_log_warning!(
                            "Spurious data before the first property will be ignored\n"
                        );
                    }
                    warning_issued = true;
                }
            }
            line = stream.getline_raw();
        }
        Ok(())
    }

    /// The lines of one data item (up to the next blank line) as a property.
    fn read_data_item(
        &self,
        stream: &mut LineStream<'_>,
        mol: &mut RWMol,
        label: &str,
        process_lists: bool,
    ) {
        let mut line = stream.getline_raw();
        let mut value: Vec<u8> = Vec::new();
        let mut n_lines = 0;
        while !strip(&line).is_empty() || matches!(line.first(), Some(b' ' | b'\t')) {
            n_lines += 1;
            if n_lines > 1 {
                value.push(b'\n');
            }
            if line.last() == Some(&b'\r') {
                line.pop();
            }
            value.extend_from_slice(&line);
            line = stream.getline_raw();
        }
        mol.props.set_val(label, lossy(&value));
        if process_lists {
            process_mol_property_list(mol, label);
        }
    }

    fn parse_record(&self, record: &Record) -> Result<Option<RWMol>> {
        let mut stream = LineStream::new(&record.text);
        let mut line = record.line;
        let Some(mut mol) = mol_from_mol_data_stream(&mut stream, &mut line, &self.params)? else {
            return Ok(None);
        };
        self.read_mol_props(&mut stream, &mut mol)?;
        Ok(Some(mol))
    }
}

impl RecordParser for SdParser {
    /// `processMoleculeRecord`: the molecule block and its data fields.
    fn parse(&self, record: &Record) -> Option<RWMol> {
        self.parse_record(record).ok().flatten()
    }
}

// -------------------------------------------------------------- SMILES files

struct SmilesSource {
    data: Vec<u8>,
    delimiter: String,
    title_line: bool,
    /// Filled with the column names of the title line before the first
    /// record is handed out.
    columns: Arc<OnceLock<Vec<String>>>,
}

fn is_empty_or_comment(line: &[u8]) -> bool {
    line.is_empty() || line[0] == b'#' || is_blank(line)
}

impl SmilesSource {
    /// `processTitleLine`: the names of the columns.
    fn title_columns(&self, stream: &mut LineStream<'_>) -> Vec<String> {
        let mut line = Vec::new();
        while !stream.eof() && !stream.fail() && is_empty_or_comment(&line) {
            line = stream.get_line();
        }
        tokenize(&lossy(&line), &self.delimiter)
            .into_iter()
            .map(|t| lossy(strip(t.as_bytes())))
            .collect()
    }
}

impl RecordSource for SmilesSource {
    /// `MultithreadedSmilesMolSupplier::extractNextRecord`: the next line
    /// that is not blank or a comment. Like upstream, the line that ends the
    /// input without a molecule is a record too.
    fn read_records(self: Box<Self>, emit: &mut dyn FnMut(Record) -> bool) {
        let mut stream = LineStream::new(&self.data);
        // the reader never counts lines: upstream hands the record's line
        // number to the parser as its initial value of -1
        let line_number = u32::MAX;
        let mut first = true;
        let mut id = 1u32;
        while !stream.eof() {
            if first {
                let columns = if self.title_line {
                    self.title_columns(&mut stream)
                } else {
                    Vec::new()
                };
                let _ = self.columns.set(columns);
                first = false;
            }
            let mut line = Vec::new();
            while !stream.eof() && !stream.fail() && is_empty_or_comment(&line) {
                line = stream.get_line();
            }
            let item = Record {
                id,
                line: line_number,
                text: line,
            };
            id += 1;
            if !emit(item) {
                return;
            }
        }
    }
}

struct SmilesParser {
    params: SmilesMolSupplierParams,
    columns: Arc<OnceLock<Vec<String>>>,
}

impl RecordParser for SmilesParser {
    /// `processMoleculeRecord`
    fn parse(&self, record: &Record) -> Option<RWMol> {
        let line = lossy(&record.text);
        let recs: Vec<String> = tokenize(&line, &self.params.delimiter)
            .into_iter()
            .map(|t| lossy(strip(t.as_bytes())))
            .collect();
        let smiles_column = usize::try_from(self.params.smiles_column).unwrap_or(usize::MAX);
        let smiles = recs.get(smiles_column)?;
        let mut res = mol_from_smiles(smiles, &self.params.parse_parameters)
            .ok()
            .flatten()?;

        let name_column = self.params.name_column;
        if name_column == -1 {
            res.set_prop(cp::_Name, record.line.to_string(), false)
                .ok()?;
        } else {
            match usize::try_from(name_column) {
                Ok(c) if c < recs.len() => {
                    res.set_prop(cp::_Name, recs[c].clone(), false).ok()?;
                }
                _ => rd_log_warning!("WARNING: no name column found on line {}", record.line),
            }
        }

        let columns = self.columns.get().map_or(&[][..], Vec::as_slice);
        for (col, value) in recs.iter().enumerate() {
            if col == smiles_column || i64::try_from(col).ok() == Some(i64::from(name_column)) {
                continue;
            }
            let name = match columns.get(col) {
                Some(n) if !n.is_empty() => n.clone(),
                _ => format!("Column_{col}"),
            };
            res.set_prop(&name, value.clone(), false).ok()?;
        }
        Some(res)
    }
}
