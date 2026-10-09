//! `v2::FileParsers::SmilesMolSupplier`: a lazy supplier over a delimited
//! text file with one SMILES per line.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/SmilesMolSupplier.cpp`
//! (2026.09.1). The input is held in memory and positions are byte offsets
//! into it, with `std::istream`'s eof/fail state tracked by [`LineStream`]
//! so that the stream-position bookkeeping matches upstream line by line.
//!
//! Copyright (C) 2002-2024 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::sd_mol_supplier::strip;
use super::utils::{LineStream, lossy};
use crate::graphmol::romol::RWMol;
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::rdgeneral::{RDProps, common_properties as cp};
use crate::{RDKitError, Result, rd_log_error, rd_log_warning};

/// `SmilesMolSupplierParams`
#[derive(Debug, Clone)]
pub struct SmilesMolSupplierParams {
    /// every character of the string separates columns (not the string itself)
    pub delimiter: String,
    pub smiles_column: i32,
    /// `-1`: no name column, molecules are named after their line number
    pub name_column: i32,
    pub title_line: bool,
    pub parse_parameters: SmilesParserParams,
}

impl Default for SmilesMolSupplierParams {
    fn default() -> Self {
        SmilesMolSupplierParams {
            delimiter: " \t".to_string(),
            smiles_column: 0,
            name_column: 1,
            title_line: true,
            parse_parameters: SmilesParserParams {
                sanitize: true,
                allow_cx_smiles: false,
                strict_cx_smiles: true,
                parse_name: false,
                remove_hs: true,
                skip_cleanup: false,
                debug_parse: false,
                replacements: Default::default(),
            },
        }
    }
}

/// `boost::tokenizer` with `char_separator(delimiter, "", keep_empty_tokens)`:
/// a split at every delimiter character that keeps empty fields.
pub(super) fn tokenize<'a>(line: &'a str, delimiter: &str) -> Vec<&'a str> {
    // (an empty input has no tokens at all)
    if line.is_empty() {
        return Vec::new();
    }
    line.split(|c| delimiter.contains(c)).collect()
}

fn no_stream() -> RDKitError {
    RDKitError::precondition("no stream")
}

/// `std::istream::tellg`: -1 (`None`) unless the stream is good — a stream
/// that only has eofbit set fails the sentry, which sets failbit.
fn tellg(s: &mut LineStream<'_>) -> Option<usize> {
    if s.eof() && !s.fail() {
        *s = LineStream::with_state(s.data(), (s.position(), true, true));
    }
    s.tellg()
}

/// C++11 `seekg(pos)`: eofbit is cleared first; a failed stream does not move.
fn seekg(s: &mut LineStream<'_>, pos: usize) {
    let (cur, _, fail) = s.state();
    *s = if fail {
        LineStream::with_state(s.data(), (cur, false, true))
    } else {
        LineStream::with_state(s.data(), (pos.min(s.data().len()), false, false))
    };
}

/// `seekg(pos)` with `pos` as returned by `tellg()` (-1 puts the stream in the fail state).
fn seekg_opt(s: &mut LineStream<'_>, pos: Option<usize>) {
    match pos {
        Some(p) => seekg(s, p),
        None => *s = LineStream::with_state(s.data(), (s.position(), false, true)),
    }
}

/// `clear()`
fn clear(s: &mut LineStream<'_>) {
    s.seek(s.position());
}

#[derive(Debug, Clone)]
pub struct SmilesMolSupplier {
    data: Option<Vec<u8>>,
    /// `(pos, eof, fail)` of the input stream
    state: (usize, bool, bool),
    params: SmilesMolSupplierParams,
    df_end: bool,
    /// total number of molecules, -1 while unknown
    d_len: i64,
    /// the molecule `next()` returns, -1 before the first read
    d_next: i64,
    d_line: i64,
    d_molpos: Vec<Option<usize>>,
    d_line_nums: Vec<i64>,
    /// column names of the title line
    d_props: Vec<String>,
}

impl Default for SmilesMolSupplier {
    /// The default constructor: no input.
    fn default() -> Self {
        SmilesMolSupplier {
            data: None,
            state: (0, false, false),
            params: SmilesMolSupplierParams::default(),
            df_end: false,
            d_len: -1,
            d_next: -1,
            d_line: -1,
            d_molpos: Vec::new(),
            d_line_nums: Vec::new(),
            d_props: Vec::new(),
        }
    }
}

impl SmilesMolSupplier {
    /// `SmilesMolSupplier(fileName, params)`
    pub fn from_file(file_name: &str, params: SmilesMolSupplierParams) -> Result<Self> {
        let data = std::fs::read(file_name)
            .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))?;
        if data.is_empty() {
            return Err(RDKitError::BadFile(format!(
                "Invalid input file {file_name}"
            )));
        }
        let mut s = SmilesMolSupplier::default();
        s.set_data_bytes(data, params);
        Ok(s)
    }

    /// `setData(text, params)`
    pub fn set_data(&mut self, text: &[u8], params: SmilesMolSupplierParams) {
        self.set_data_bytes(text.to_vec(), params);
    }

    fn set_data_bytes(&mut self, data: Vec<u8>, params: SmilesMolSupplierParams) {
        // init()
        *self = SmilesMolSupplier {
            data: Some(data),
            params,
            ..Default::default()
        };
        self.check_for_end();
    }

    /// `close()`
    pub fn close(&mut self) {
        self.data = None;
    }

    fn with_stream<T>(&mut self, f: impl FnOnce(&mut Self, &mut LineStream<'_>) -> T) -> Result<T> {
        let data = self.data.take().ok_or_else(no_stream)?;
        let mut s = LineStream::with_state(&data, self.state);
        let r = f(self, &mut s);
        self.state = s.state();
        self.data = Some(data);
        Ok(r)
    }

    /// `checkForEnd`: makes sure there is a line to read.
    fn check_for_end(&mut self) {
        let _ = self.with_stream(|this, s| {
            if this.skip_comments(s).is_some() {
                this.d_line = -1;
                seekg(s, 0);
                this.df_end = false;
            }
        });
    }

    /// `atEnd()`
    pub fn at_end(&self) -> Result<bool> {
        self.data.as_ref().ok_or_else(no_stream)?;
        Ok(self.df_end)
    }

    /// `reset()`
    pub fn reset(&mut self) -> Result<()> {
        self.with_stream(|this, s| {
            clear(s);
            this.df_end = false;
            match this.d_molpos.first() {
                Some(&first) => {
                    seekg_opt(s, first);
                    this.d_next = 0;
                    this.d_line = 0;
                }
                None => {
                    seekg(s, 0);
                    this.d_next = -1;
                    this.d_line = -1;
                }
            }
        })
    }

    /// `processLine`: a molecule from one line of the file, `None` (after
    /// logging) when the line cannot be turned into one.
    fn process_line(&self, line: &str) -> Option<RWMol> {
        let recs: Vec<String> = tokenize(line, &self.params.delimiter)
            .into_iter()
            .map(|t| lossy(strip(t.as_bytes())))
            .collect();
        let smiles_column = usize::try_from(self.params.smiles_column).unwrap_or(usize::MAX);
        if recs.len() <= smiles_column {
            rd_log_error!("ERROR: Could not process molecule on line {}", self.d_line);
            return None;
        }
        let mut res = match mol_from_smiles(&recs[smiles_column], &self.params.parse_parameters) {
            Ok(Some(m)) => m,
            Ok(None) => {
                rd_log_error!("ERROR: Smiles parse error on line {}", self.d_line);
                rd_log_error!(
                    "ERROR: Cannot create molecule from : '{}'",
                    recs[smiles_column]
                );
                return None;
            }
            Err(e @ RDKitError::Sanitize(_)) => {
                rd_log_error!("ERROR: Could not sanitize molecule on line {}", self.d_line);
                rd_log_error!("ERROR: {e}");
                return None;
            }
            Err(_) => {
                rd_log_error!("ERROR: Could not process molecule on line {}", self.d_line);
                return None;
            }
        };

        // the name (if there's a name column)
        if self.params.name_column == -1 {
            // no name column: default to the line number
            res.set_prop(cp::_Name, self.d_line.to_string(), false)
                .ok()?;
        } else {
            match usize::try_from(self.params.name_column) {
                Ok(c) if c < recs.len() => {
                    res.set_prop(cp::_Name, recs[c].clone(), false).ok()?;
                }
                _ => rd_log_warning!("WARNING: no name column found on line {}", self.d_line),
            }
        }

        // the remaining columns become properties
        for (col, value) in recs.iter().enumerate() {
            if col == smiles_column
                || i64::try_from(col).ok() == Some(i64::from(self.params.name_column))
            {
                continue;
            }
            let mut name = self.d_props.get(col).cloned().unwrap_or_default();
            if name.is_empty() {
                name = format!("Column_{col}");
            }
            res.set_prop(&name, value.clone(), false).ok()?;
        }
        Some(res)
    }

    /// `nextLine`: the next line of the stream. Sets `df_end` when EOF is hit
    /// without reading anything and counts the lines read.
    fn next_line(&mut self, s: &mut LineStream<'_>) -> Vec<u8> {
        if self.df_end {
            return Vec::new();
        }
        let line = s.get_line();
        if line.is_empty() {
            // got an empty string, check to see if we hit EOF
            if s.eof() || s.fail() {
                self.df_end = true;
            }
        } else if s.eof() {
            // we got some data before hitting EOF: clear the flag
            clear(s);
        }
        self.d_line += 1;
        line
    }

    /// `skipComments`: the position of the beginning of the next
    /// non-comment line (`None` when there is none).
    fn skip_comments(&mut self, s: &mut LineStream<'_>) -> Option<usize> {
        if self.df_end {
            return None;
        }
        let mut prev = tellg(s);
        let mut line = self.next_line(s);
        let skippable = |l: &[u8]| l.first() == Some(&b'#') || strip(l).is_empty();
        if !self.df_end {
            // if we didn't immediately hit EOF, loop until we get a valid line
            while skippable(&line) {
                prev = tellg(s);
                line = self.next_line(s);
                if self.df_end {
                    break;
                }
            }
        }
        // if we hit EOF without getting a proper line, there is none
        if line.is_empty() || skippable(&line) {
            return None;
        }
        prev
    }

    /// `processTitleLine`: reads the column names.
    fn process_title_line(&mut self, s: &mut LineStream<'_>) {
        if let Some(pos) = self.skip_comments(s) {
            seekg(s, pos);
            let line = lossy(&s.get_line());
            for token in tokenize(&line, &self.params.delimiter) {
                self.d_props.push(lossy(strip(token.as_bytes())));
            }
        }
    }

    fn ran_out_of_lines(idx: usize) -> RDKitError {
        RDKitError::FileParse(format!(
            "ERROR: Index error (idx = {idx}): ran out of lines\n"
        ))
    }

    /// `moveTo`: positions the stream at the beginning of entry `idx`.
    fn move_to_in(&mut self, s: &mut LineStream<'_>, idx: usize) -> Result<()> {
        // the easy situations (boundary conditions) first
        if self.d_len > -1 && idx as i64 >= self.d_len {
            self.df_end = true;
            return Err(Self::ran_out_of_lines(idx));
        }
        // earlier calls to next() may have put the stream into a bad state
        clear(s);

        // case 1: we have already read the entry
        if self.d_molpos.len() > idx {
            self.df_end = false;
            seekg_opt(s, self.d_molpos[idx]);
            self.d_next = idx as i64;
            self.d_line = self.d_line_nums[idx];
            return Ok(());
        }

        // case 2: move forward until we've gone far enough
        if self.d_molpos.is_empty() {
            // just starting out: process the title line
            seekg(s, 0);
            if self.params.title_line {
                self.process_title_line(s);
            }
        } else {
            // move to the last position we've seen and read that line
            seekg_opt(s, *self.d_molpos.last().expect("not empty"));
            s.get_line();
        }

        while self.d_molpos.len() <= idx {
            let Some(next_pos) = self.skip_comments(s) else {
                return Err(Self::ran_out_of_lines(idx));
            };
            self.d_molpos.push(Some(next_pos));
            self.d_line_nums.push(self.d_line);
            if self.d_molpos.len() == idx + 1 && self.df_end {
                // we could read the point we were looking for but not the
                // next one: indicate that we've reached EOF
                clear(s);
                seekg(s, s.data().len());
                self.d_len = self.d_molpos.len() as i64;
                break;
            }
        }
        seekg_opt(s, self.d_molpos[idx]);
        self.d_next = idx as i64;
        Ok(())
    }

    /// `moveTo(idx)`
    pub fn move_to(&mut self, idx: usize) -> Result<()> {
        self.with_stream(|this, s| this.move_to_in(s, idx))?
    }

    fn next_in(&mut self, s: &mut LineStream<'_>) -> Result<Option<RWMol>> {
        if self.d_next < 0 {
            self.d_next = 0;
        }
        // this throws an exception if it fails
        self.move_to_in(s, self.d_next as usize)?;
        let next = self.d_next as usize;
        if self.d_molpos.len() <= next {
            return Err(RDKitError::invariant("bad index length"));
        }
        // set the stream to the relevant position, grab the line and process it
        clear(s);
        seekg_opt(s, self.d_molpos[next]);
        self.d_line = self.d_line_nums[next];
        let line = lossy(&s.get_line());
        let res = self.process_line(&line);
        // if we don't already know the length of the supplier, check if we
        // can read another line
        if self.d_len < 0 && self.skip_comments(s).is_none() {
            self.d_len = self.d_molpos.len() as i64;
        }
        // make sure the line number is correct
        if let Some(&n) = self.d_line_nums.get(next) {
            self.d_line = n;
        }
        self.d_next += 1;
        // if we just hit the last one, simulate EOF
        if self.d_len > 0 && self.d_next == self.d_len {
            self.df_end = true;
        }
        Ok(res)
    }

    /// `next()`: the next molecule, `None` when its line cannot be parsed.
    pub fn next_mol(&mut self) -> Result<Option<RWMol>> {
        self.with_stream(|this, s| this.next_in(s))?
    }

    /// `operator[](idx)`
    pub fn get(&mut self, idx: usize) -> Result<Option<RWMol>> {
        self.with_stream(|this, s| {
            this.move_to_in(s, idx)?;
            this.next_in(s)
        })?
    }

    /// `getItemText(idx)`: the line of entry `idx`.
    pub fn get_item_text(&mut self, idx: usize) -> Result<String> {
        self.with_stream(|this, s| {
            // upstream keeps these in an unsigned int, so -1 comes back as 2^32-1
            let holder = this.d_next as u32;
            let end_holder = this.df_end;
            // this throws the relevant exception if we go too far
            this.move_to_in(s, idx)?;
            let res = lossy(&s.get_line());
            this.d_next = i64::from(holder);
            this.df_end = end_holder;
            Ok(res)
        })?
    }

    /// `length()`: the number of molecules (reads the whole input once).
    pub fn length(&mut self) -> Result<usize> {
        self.with_stream(|this, s| {
            if this.d_len > 0 {
                return this.d_len as usize;
            }
            let old_pos = tellg(s);
            if let Some(&last) = this.d_molpos.last() {
                // we've already read some molecules: go to the last one and
                // read it in to initialize our location
                seekg_opt(s, last);
                // skip that line and then continue
                this.skip_comments(s);
            } else if this.params.title_line {
                this.process_title_line(s);
            }
            while let Some(pos) = this.skip_comments(s) {
                this.d_molpos.push(Some(pos));
                this.d_line_nums.push(this.d_line);
            }
            // remember to set the stream to its original position
            seekg_opt(s, old_pos);
            this.d_len = this.d_molpos.len() as i64;
            this.d_len as usize
        })
    }
}

#[cfg(test)]
mod tokenize_tests {
    use super::tokenize;

    #[test]
    fn splits_at_every_delimiter_character_and_keeps_empty_fields() {
        assert_eq!(tokenize("a,,b,", ","), ["a", "", "b", ""]);
        assert_eq!(tokenize(",", ","), ["", ""]);
        assert_eq!(tokenize("a b\tc", " \t"), ["a", "b", "c"]);
    }

    #[test]
    fn empty_input_has_no_tokens() {
        assert!(tokenize("", ",").is_empty());
    }
}
