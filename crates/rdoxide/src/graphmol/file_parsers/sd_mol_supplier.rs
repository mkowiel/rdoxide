//! SD file suppliers (`v2::FileParsers::ForwardSDMolSupplier` /
//! `SDMolSupplier`) and SD property lists (`processMolPropertyList`).
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp`,
//! `SDMolSupplier.cpp` and `FileParserUtils.h` (2026.09.1). The input is held
//! in memory; stream positions are byte offsets.
//!
//! Copyright (C) 2002-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::mol_file_parser::{MolFileParserParams, mol_from_mol_data_stream};
use super::utils::{LineStream, lexical_cast_int, lossy};
use crate::graphmol::romol::RWMol;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_error, rd_log_warning};

/// `FileParserUtils::strip` with the default " \t\r\n"
fn strip(s: &[u8]) -> &[u8] {
    let ws = |c: &u8| matches!(c, b' ' | b'\t' | b'\r' | b'\n');
    let start = s.iter().position(|c| !ws(c)).unwrap_or(s.len());
    let end = s.iter().rposition(|c| !ws(c)).map_or(start, |e| e + 1);
    &s[start..end]
}

fn is_dollars(t: &[u8]) -> bool {
    t.starts_with(b"$$$$")
}

/// `boost::split(tokens, s, is_any_of(" \t\n"), token_compress_on)`
fn split_compress(s: &str) -> Vec<&str> {
    let mut res = Vec::new();
    let mut cur = 0;
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() {
        if matches!(b[i], b' ' | b'\t' | b'\n') {
            res.push(&s[cur..i]);
            while i < b.len() && matches!(b[i], b' ' | b'\t' | b'\n') {
                i += 1;
            }
            cur = i;
        } else {
            i += 1;
        }
    }
    res.push(&s[cur..]);
    res
}

#[derive(Clone, Copy)]
enum ListType {
    Str,
    Int,
    Double,
    Bool,
}

fn parse_list_value(t: ListType, v: &str) -> Option<RDValue> {
    match t {
        ListType::Str => Some(RDValue::String(v.to_string())),
        ListType::Int => lexical_cast_int(v.as_bytes()).ok().map(RDValue::Int),
        ListType::Double => {
            let ok = !v.is_empty()
                && v.bytes()
                    .all(|c| c.is_ascii_digit() || b"+-.eE".contains(&c));
            if ok {
                v.parse::<f64>().ok().map(RDValue::Double)
            } else {
                None
            }
        }
        ListType::Bool => match v {
            "1" => Some(RDValue::Bool(true)),
            "0" => Some(RDValue::Bool(false)),
            _ => None,
        },
    }
}

/// `FileParserUtils::applyMolListProp`
fn apply_mol_list_prop(
    mol: &mut RWMol,
    pn: &str,
    prefix: &str,
    missing: &str,
    atoms: bool,
    t: ListType,
) {
    let itempn = &pn[prefix.len()..];
    let Some(RDValue::String(sv)) = mol.get_prop(pn).cloned() else {
        return;
    };
    let tokens = split_compress(&sv);
    let n_items = if atoms {
        mol.num_atoms()
    } else {
        mol.num_bonds()
    } as usize;
    let mut mv = missing.to_string();
    let mut first = 0;
    if tokens.len() == n_items + 1
        && tokens[0].starts_with('[')
        && tokens[0].ends_with(']')
        && !tokens[0].is_empty()
    {
        mv = tokens[0][1..tokens[0].len().max(2) - 1].to_string();
        first = 1;
    }
    if mv.is_empty() {
        rd_log_warning!("Missing value marker for property {pn} is empty.");
    }
    if tokens.len() - first != n_items {
        rd_log_warning!(
            "Property list {pn} has incompatible size, {} elements found; expecting {n_items}. Ignoring it.",
            tokens.len()
        );
        return;
    }
    for (i, tok) in tokens.iter().enumerate().skip(first) {
        if *tok == mv {
            continue;
        }
        let item = (i - first) as u32;
        match parse_list_value(t, tok) {
            Some(v) => {
                if atoms {
                    mol.atom_mut(item).props.set_val(itempn, v);
                } else {
                    mol.bond_mut(item).props.set_val(itempn, v);
                }
            }
            None => rd_log_warning!(
                "Value {tok} for property {pn} of item {item} can not be parsed. Ignoring it."
            ),
        }
    }
}

/// `FileParserUtils::processMolPropertyList(mol, pn, "n/a")`
pub fn process_mol_property_list(mol: &mut RWMol, pn: &str) {
    for (base, atoms) in [("atom.", true), ("bond.", false)] {
        if pn.starts_with(base) && pn.len() > base.len() {
            for (kind, t) in [
                ("prop.", ListType::Str),
                ("iprop.", ListType::Int),
                ("dprop.", ListType::Double),
                ("bprop.", ListType::Bool),
            ] {
                let prefix = format!("{base}{kind}");
                if pn.starts_with(&prefix) && pn.len() > prefix.len() {
                    apply_mol_list_prop(mol, pn, &prefix, "n/a", atoms, t);
                    return;
                }
            }
            return;
        }
    }
}

/// `FileParserUtils::processMolPropertyLists`
pub fn process_mol_property_lists(mol: &mut RWMol) {
    for pn in mol.get_prop_list(true, true) {
        process_mol_property_list(mol, &pn);
    }
}

/// `v2::FileParsers::SDMolSupplier` (which is also a
/// `ForwardSDMolSupplier`) over an in-memory buffer.
#[derive(Debug, Clone)]
pub struct SDMolSupplier {
    data: Option<Vec<u8>>,
    state: (usize, bool, bool),
    pub params: MolFileParserParams,
    pub process_property_lists: bool,
    df_end: bool,
    eof_hit_on_read: bool,
    d_line: i64,
    d_len: i64,
    d_last: usize,
    d_molpos: Vec<usize>,
}

fn no_stream() -> RDKitError {
    RDKitError::precondition("no stream")
}

impl Default for SDMolSupplier {
    /// The default constructor: no input.
    fn default() -> Self {
        SDMolSupplier {
            data: None,
            state: (0, false, false),
            params: MolFileParserParams::default(),
            process_property_lists: true,
            df_end: false,
            eof_hit_on_read: false,
            d_line: 0,
            d_len: -1,
            d_last: 0,
            d_molpos: Vec::new(),
        }
    }
}

impl SDMolSupplier {
    /// `SDMolSupplier(fileName, params)`
    pub fn from_file(file_name: &str, params: MolFileParserParams) -> Result<Self> {
        let data = std::fs::read(file_name)
            .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))?;
        if data.is_empty() {
            return Err(RDKitError::BadFile(format!(
                "Invalid input file {file_name}"
            )));
        }
        let mut s = SDMolSupplier {
            params,
            ..Default::default()
        };
        s.set_data_bytes(data);
        Ok(s)
    }

    /// `setData(text, params)`
    pub fn set_data(&mut self, text: &[u8], params: MolFileParserParams) {
        self.params = params;
        self.set_data_bytes(text.to_vec());
    }

    fn set_data_bytes(&mut self, data: Vec<u8>) {
        // init()
        self.df_end = false;
        self.d_line = 0;
        self.process_property_lists = true;
        self.d_len = -1;
        self.d_last = 0;
        self.d_molpos.clear();
        self.data = Some(data);
        self.state = (0, false, false);
        self.d_molpos.push(0);
        self.check_for_end();
        if self.df_end {
            self.d_len = 0;
        }
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

    /// `SDMolSupplier::checkForEnd`
    fn check_for_end(&mut self) {
        let _ = self.with_stream(|this, s| this.check_for_end_in(s));
    }

    fn check_for_end_in(&mut self, s: &mut LineStream<'_>) {
        if s.eof() {
            self.df_end = true;
            self.d_len = self.d_molpos.len() as i64;
            return;
        }
        let mut nempty = 0;
        for _ in 0..4 {
            let t = s.get_line();
            if s.eof() {
                self.df_end = true;
                self.d_len = self.d_molpos.len() as i64;
                return;
            }
            if t.iter().all(|c| matches!(c, b' ' | b'\t' | b'\r' | b'\n')) {
                nempty += 1;
            }
        }
        if nempty == 4 {
            self.df_end = true;
            self.d_len = self.d_molpos.len() as i64;
        }
    }

    /// `atEnd()`
    pub fn at_end(&self) -> Result<bool> {
        self.data.as_ref().ok_or_else(no_stream)?;
        Ok(self.df_end)
    }

    /// `getEOFHitOnRead()`
    pub fn eof_hit_on_read(&self) -> bool {
        self.eof_hit_on_read
    }

    /// `reset()`
    pub fn reset(&mut self) -> Result<()> {
        self.data.as_ref().ok_or_else(no_stream)?;
        self.state = (0, false, false);
        self.df_end = false;
        self.d_last = 0;
        self.d_line = 0;
        Ok(())
    }

    /// `ForwardSDMolSupplier::readMolProps`
    fn read_mol_props(&mut self, s: &mut LineStream<'_>, mol: &mut RWMol) -> Result<()> {
        self.d_line += 1;
        let mut has_prop = false;
        let mut warning_issued = false;
        let mut dlabel = String::new();
        let mut t = s.getline_raw();
        while !s.eof() && !s.fail() && (t.is_empty() || !is_dollars(&t)) {
            let st = strip(&t).to_vec();
            if !st.is_empty() {
                if st[0] == b'>' {
                    has_prop = true;
                    warning_issued = false;
                    let rest = &st[1..];
                    let sl = rest.iter().position(|&c| c == b'<');
                    let se = rest.iter().rposition(|&c| c == b'>');
                    match (sl, se) {
                        (Some(sl), Some(se)) if se != sl + 1 => {
                            // upstream: substr(sl + 1, se - sl - 1) (unsigned
                            // wrap-around when '>' precedes '<')
                            dlabel = if se > sl {
                                lossy(&rest[sl + 1..se])
                            } else {
                                lossy(&rest[sl + 1..])
                            };
                            self.d_line += 1;
                            let mut t2 = s.getline_raw();
                            let mut prop: Vec<u8> = Vec::new();
                            let mut stmp = strip(&t2).to_vec();
                            let mut nplines = 0;
                            while !stmp.is_empty()
                                || (!t2.is_empty() && matches!(t2[0], b' ' | b'\t'))
                            {
                                nplines += 1;
                                if nplines > 1 {
                                    prop.push(b'\n');
                                }
                                if !t2.is_empty() {
                                    if t2.last() == Some(&b'\r') {
                                        t2.pop();
                                    }
                                    prop.extend_from_slice(&t2);
                                }
                                self.d_line += 1;
                                t2 = s.getline_raw();
                                stmp = strip(&t2).to_vec();
                            }
                            mol.props.set_val(&dlabel, lossy(&prop));
                            if self.process_property_lists {
                                process_mol_property_list(mol, &dlabel);
                            }
                        }
                        _ => {
                            self.d_line += 1;
                            let t2 = s.getline_raw();
                            // upstream never updates `stmp` in this loop: a
                            // non-blank line reads on to the end of the input
                            let blank = strip(&t2).is_empty();
                            if !blank {
                                loop {
                                    self.d_line += 1;
                                    s.getline_raw();
                                    if s.eof() {
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
                            "Property <{dlabel}> will be truncated after the first blank line"
                        );
                    } else {
                        rd_log_warning!("Spurious data before the first property will be ignored");
                    }
                    warning_issued = true;
                }
            }
            self.d_line += 1;
            t = s.getline_raw();
        }
        Ok(())
    }

    /// Skips to the line after the next `$$$$` (the recovery loop of the
    /// exception handlers in `_next`).
    fn skip_to_next_record(&mut self, s: &mut LineStream<'_>, set_eof_flag: bool) {
        self.d_line += 1;
        let mut t = s.getline_raw();
        if set_eof_flag && s.eof() {
            self.eof_hit_on_read = true;
        }
        while !s.eof() && !s.fail() && (t.is_empty() || !is_dollars(&t)) {
            self.d_line += 1;
            t = s.getline_raw();
        }
    }

    /// `ForwardSDMolSupplier::_next`
    fn next_in(&mut self, s: &mut LineStream<'_>) -> Option<RWMol> {
        if s.eof() {
            self.df_end = true;
            return None;
        }
        self.eof_hit_on_read = false;
        let mut line = self.d_line as u32;
        let mut res: Option<RWMol> = None;
        let params = self.params.clone();
        let outcome: Result<()> = (|| {
            res = mol_from_mol_data_stream(s, &mut line, &params)?;
            if res.is_none() && s.eof() && (line as i64 - self.d_line) < 2 {
                self.eof_hit_on_read = true;
            }
            self.d_line = line as i64;
            if let Some(m) = res.as_mut() {
                self.read_mol_props(s, m)?;
            } else if !s.eof() {
                let mut t = s.getline_raw();
                self.d_line += 1;
                while !s.eof() && !s.fail() {
                    // tempStr.at(0) on an empty line throws std::out_of_range
                    if t.is_empty() {
                        return Err(RDKitError::Index("basic_string::at".into()));
                    }
                    if is_dollars(&t) {
                        break;
                    }
                    t = s.getline_raw();
                    self.d_line += 1;
                }
            }
            Ok(())
        })();
        if let Err(e) = outcome {
            if self.d_line < line as i64 {
                self.d_line = line as i64;
            }
            match e {
                RDKitError::FileParse(m) => {
                    rd_log_error!("ERROR: {m}");
                    rd_log_error!("ERROR: moving to the beginning of the next molecule\n");
                    self.skip_to_next_record(s, false);
                }
                RDKitError::Sanitize(se) => {
                    rd_log_error!(
                        "ERROR: Could not sanitize molecule ending on line {}",
                        self.d_line
                    );
                    rd_log_error!("ERROR: {}\n", RDKitError::Sanitize(se));
                    self.skip_to_next_record(s, true);
                }
                _ => {
                    if s.eof() {
                        self.eof_hit_on_read = true;
                    }
                    rd_log_error!("Unexpected error hit on line {}", self.d_line);
                    rd_log_error!("ERROR: moving to the beginning of the next molecule\n");
                    self.skip_to_next_record(s, true);
                }
            }
        }
        if s.eof() {
            self.df_end = true;
        }
        res
    }

    /// `SDMolSupplier::next`
    pub fn next_mol(&mut self) -> Result<Option<RWMol>> {
        self.data.as_ref().ok_or_else(no_stream)?;
        if self.df_end && self.d_last as i64 >= self.d_len {
            return Err(RDKitError::FileParse("EOF hit.".into()));
        }
        self.with_stream(|this, s| {
            s.seekg(this.d_molpos[this.d_last]);
            if s.eof() {
                this.df_end = true;
                this.d_len = this.d_molpos.len() as i64;
                return None;
            }
            let res = this.next_in(s);
            this.d_last += 1;
            let pos_hold = s.tellg();
            this.check_for_end_in(s);
            if !this.df_end && this.d_last >= this.d_molpos.len() {
                // tellg() of a failed stream is -1; then df_end is set above
                this.d_molpos.push(pos_hold.unwrap_or(usize::MAX));
            }
            res
        })
    }

    /// The data a `$$$$` scan sees from `pos` (`buildIndexTo`): positions
    /// after each `$$$$` line that starts a further record.
    fn build_index_to(&mut self, target_idx: usize) -> Result<()> {
        let data = self.data.take().ok_or_else(no_stream)?;
        let start = *self.d_molpos.last().expect("molpos");
        self.d_last = self.d_molpos.len() - 1;
        let mut ptr = start;
        // the search starts one byte into an (initially '\n'-filled) overlap
        let preceded_by_nl = |i: usize| i == start || (i > 0 && data[i - 1] == b'\n');
        while ptr + 4 <= data.len() {
            let Some(off) = data[ptr..].windows(4).position(|w| w == b"$$$$") else {
                break;
            };
            let m = ptr + off;
            if preceded_by_nl(m) {
                let mut nl = m + 4;
                while nl < data.len() && data[nl] != b'\n' {
                    nl += 1;
                }
                if nl < data.len() {
                    nl += 1;
                }
                let at_true_eof = nl >= data.len();
                if !at_true_eof {
                    // peekCheckForEnd
                    let mut empty = 0;
                    let mut decided = false;
                    for &c in &data[nl..] {
                        if !c.is_ascii_whitespace() && c != 0x0b {
                            decided = true;
                            break;
                        }
                        if c == b'\n' {
                            empty += 1;
                            if empty >= 4 {
                                self.df_end = true;
                                self.d_len = self.d_molpos.len() as i64;
                                decided = true;
                                break;
                            }
                        }
                    }
                    if !decided {
                        let mut s = LineStream::with_state(&data, (nl, false, false));
                        self.check_for_end_in(&mut s);
                    }
                    if !self.df_end {
                        self.d_molpos.push(nl);
                        self.d_last += 1;
                        if self.d_last == target_idx {
                            break;
                        }
                    }
                }
            }
            ptr = m + 4;
        }
        self.data = Some(data);
        Ok(())
    }

    /// `moveTo(idx)`
    pub fn move_to(&mut self, idx: usize) -> Result<()> {
        self.data.as_ref().ok_or_else(no_stream)?;
        self.state.1 = false;
        self.state.2 = false;
        if idx < self.d_molpos.len() {
            self.state.0 = self.d_molpos[idx];
            self.d_last = idx;
            return Ok(());
        }
        self.build_index_to(idx)?;
        if idx < self.d_molpos.len() {
            self.state = (self.d_molpos[idx], false, false);
            self.d_last = idx;
            Ok(())
        } else {
            self.df_end = true;
            self.d_len = self.d_molpos.len() as i64;
            Err(RDKitError::FileParse(format!(
                "ERROR: Index error (idx = {idx}) :  we do not have enough mol blocks"
            )))
        }
    }

    /// `operator[](idx)`
    pub fn get(&mut self, idx: usize) -> Result<Option<RWMol>> {
        self.move_to(idx)?;
        self.next_mol()
    }

    /// `length()`
    pub fn length(&mut self) -> Result<usize> {
        self.data.as_ref().ok_or_else(no_stream)?;
        if self.d_len > 0 || (self.df_end && self.d_len == 0) {
            return Ok(self.d_len as usize);
        }
        let old_last = self.d_last;
        self.build_index_to(u32::MAX as usize)?;
        self.d_len = self.d_molpos.len() as i64;
        self.d_last = old_last;
        self.state = (self.d_molpos[self.d_last], false, false);
        self.df_end = false;
        Ok(self.d_len as usize)
    }

    /// `getItemText(idx)`
    pub fn get_item_text(&mut self, idx: usize) -> Result<Vec<u8>> {
        self.data.as_ref().ok_or_else(no_stream)?;
        let holder = self.d_last;
        self.move_to(idx)?;
        let beg = self.d_molpos[idx];
        let end = match self.move_to(idx + 1) {
            Ok(()) => self.d_molpos[idx + 1],
            Err(RDKitError::FileParse(_)) => {
                let len = self.data.as_ref().map_or(0, Vec::len);
                self.state = (len, false, false);
                len
            }
            Err(e) => return Err(e),
        };
        self.d_last = holder;
        let data = self.data.as_ref().ok_or_else(no_stream)?;
        let end = end.min(data.len());
        self.state = (end, false, false);
        Ok(data[beg.min(end)..end].to_vec())
    }

    /// `setStreamIndices(locs)`
    pub fn set_stream_indices(&mut self, locs: Vec<usize>) -> Result<()> {
        self.d_molpos = locs;
        self.reset()?;
        self.d_len = self.d_molpos.len() as i64;
        Ok(())
    }
}
