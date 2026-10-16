//! `v2::FileParsers::TDTMolSupplier`: a lazy supplier over a Daylight TDT
//! file, whose records start with a `$SMI<...>` line and end with `|`.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/TDTMolSupplier.cpp`
//! (2026.09.1). Like the SMILES supplier the input is held in memory with
//! `std::istream`'s eof/fail state tracked by [`LineStream`], so that the
//! stream-position bookkeeping matches upstream line by line.
//!
//! Copyright (C) 2010 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::smiles_mol_supplier::{clear, no_stream, seekg, seekg_opt, tellg};
use super::utils::{LineStream, lossy};
use crate::geometry::Point3D;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::romol::RWMol;
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::rdgeneral::{RDProps, common_properties as cp};
use crate::{RDKitError, Result, rd_log_error};

/// `TDTMolSupplierParams`
#[derive(Debug, Clone)]
pub struct TdtMolSupplierParams {
    /// the property whose value is also the molecule's name
    pub name_record: String,
    /// the id of the 2D conformer made from the `2D` record, -1 to ignore it
    pub conf_id_2d: i32,
    /// the id of the 3D conformer made from the `3D` record, -1 to ignore it
    pub conf_id_3d: i32,
    pub parse_parameters: SmilesParserParams,
}

impl Default for TdtMolSupplierParams {
    fn default() -> Self {
        TdtMolSupplierParams {
            name_record: String::new(),
            conf_id_2d: -1,
            conf_id_3d: -1,
            parse_parameters: SmilesParserParams::default(),
        }
    }
}

/// `TDTParseUtils::ParseNumberList`: the numbers of a comma-separated list
/// that ends with `;>`, which may continue on the following lines.
fn parse_number_list(
    mut line: String,
    s: Option<&mut LineStream<'_>>,
    res: &mut Vec<f64>,
) -> Result<()> {
    let mut s = s;
    let mut found_end = false;
    loop {
        for tok in line.split(',') {
            let at_end = tok.contains(";>");
            let number = tok.trim_matches(|c| " \r\n\t;>".contains(c));
            if !number.is_empty() && !at_end {
                res.push(
                    number
                        .parse::<f64>()
                        .map_err(|_| RDKitError::Value("bad lexical cast".into()))?,
                );
            } else if at_end {
                // that's it, we're done:
                found_end = true;
                break;
            }
        }
        match s.as_deref_mut() {
            Some(st) if !found_end && !st.eof() => line = lossy(&st.getline_raw()),
            _ => break,
        }
    }
    if !found_end {
        return Err(RDKitError::FileParse(
            "no end tag found for numeric list".into(),
        ));
    }
    Ok(())
}

#[derive(Debug, Clone, Default)]
pub struct TdtMolSupplier {
    data: Option<Vec<u8>>,
    /// `(pos, eof, fail)` of the input stream
    state: (usize, bool, bool),
    params: TdtMolSupplierParams,
    df_end: bool,
    /// total number of molecules, -1 while unknown
    d_len: i64,
    /// the record `next()` reads
    d_last: usize,
    d_line: i64,
    d_molpos: Vec<Option<usize>>,
}

/// `std::string::find(c)` for a char sequence in a byte string
fn find(line: &str, c: char) -> Option<usize> {
    line.find(c)
}

impl TdtMolSupplier {
    fn blank(data: Option<Vec<u8>>, params: TdtMolSupplierParams) -> Self {
        TdtMolSupplier {
            data,
            state: (0, false, false),
            params,
            df_end: false,
            d_len: -1,
            d_last: 0,
            d_line: 0,
            d_molpos: Vec::new(),
        }
    }

    /// `TDTMolSupplier(fileName, params)`
    pub fn from_file(file_name: &str, params: TdtMolSupplierParams) -> Result<Self> {
        let data = std::fs::read(file_name)
            .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))?;
        if data.is_empty() {
            return Err(RDKitError::BadFile(format!(
                "Invalid input file {file_name}"
            )));
        }
        Ok(Self::from_bytes(data, params))
    }

    /// The supplier over `text` (`setData` on a new supplier).
    pub fn from_bytes(data: Vec<u8>, params: TdtMolSupplierParams) -> Self {
        let mut s = Self::blank(Some(data), params);
        s.start();
        s
    }

    /// The default constructor: no input.
    pub fn new() -> Self {
        Self::blank(None, TdtMolSupplierParams::default())
    }

    /// `setData(text, params)`
    pub fn set_data(&mut self, text: &[u8], params: TdtMolSupplierParams) {
        *self = Self::from_bytes(text.to_vec(), params);
    }

    fn start(&mut self) {
        let _ = self.with_stream(|this, s| {
            this.advance_to_next_record(s);
            let pos = tellg(s);
            this.d_molpos.push(pos);
            this.check_for_end(s);
        });
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

    /// `advanceToNextRecord`: positions the stream at the next `$SMI<` line.
    fn advance_to_next_record(&mut self, s: &mut LineStream<'_>) -> bool {
        let mut pos;
        loop {
            if s.eof() {
                return false;
            }
            pos = tellg(s);
            let line = lossy(&s.getline_raw());
            if line.starts_with("$SMI<") {
                break;
            }
        }
        clear(s);
        seekg_opt(s, pos);
        true
    }

    fn check_for_end(&mut self, s: &mut LineStream<'_>) {
        if s.eof() {
            self.df_end = true;
            // the -1 here is because by the time we get here we've already
            // pushed on the position of the next line
            self.d_len = self.d_molpos.len() as i64 - 1;
            return;
        }
        // we are not at the end of file, but check for blank lines:
        let temp = lossy(&s.getline_raw());
        if temp.trim_start_matches([' ', '\t', '\r', '\n']).is_empty() {
            self.df_end = true;
            self.d_len = self.d_molpos.len() as i64 - 1;
        }
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
            seekg(s, 0);
            this.df_end = false;
            this.d_last = 0;
            this.d_line = 0;
        })
    }

    fn parse_mol(&mut self, s: &mut LineStream<'_>, in_line: String) -> Result<Option<RWMol>> {
        let start_p = find(&in_line, '<');
        let end_p = in_line.rfind('>');
        // (a missing '>' gives std::string::npos, so the rest of the line)
        let from = start_p.map_or(0, |p| p + 1);
        let to = end_p.map_or(in_line.len(), |e| e.max(from));
        let smiles = &in_line[from.min(in_line.len())..to.min(in_line.len())];
        let res = mol_from_smiles(smiles, &self.params.parse_parameters)?;
        let Some(mut res) = res else {
            return Ok(None);
        };
        if res.num_atoms() == 0 {
            return Ok(Some(res));
        }
        // Process the properties:
        self.d_line += 1;
        let mut in_line = lossy(&s.getline_raw());
        while !s.eof() && !s.fail() && !in_line.starts_with('|') {
            let end_p = find(&in_line, '<');
            let prop_name = end_p
                .map_or(in_line.as_str(), |e| &in_line[..e])
                .trim_matches([' ', '\t'])
                .to_string();
            // (npos + 1 wraps to 0)
            let start_p = end_p.map_or(0, |e| e + 1);
            let rest = || in_line[start_p.min(in_line.len())..].to_string();
            if prop_name == cp::TWOD && self.params.conf_id_2d >= 0 {
                let mut coords = Vec::new();
                parse_number_list(rest(), Some(s), &mut coords)?;
                let mut conf = Conformer::new(res.num_atoms() as usize);
                conf.set_id(self.params.conf_id_2d as u32);
                conf.set_3d(false);
                for at in 0..res.num_atoms() as usize {
                    if 2 * at + 1 < coords.len() {
                        conf.set_atom_pos(
                            at as u32,
                            Point3D::new(coords[2 * at], coords[2 * at + 1], 0.0),
                        )?;
                    }
                }
                res.add_conformer(conf, false)?;
            } else if prop_name == "3D" && self.params.conf_id_3d >= 0 {
                let mut coords = Vec::new();
                parse_number_list(rest(), Some(s), &mut coords)?;
                let mut conf = Conformer::new(res.num_atoms() as usize);
                conf.set_id(self.params.conf_id_3d as u32);
                conf.set_3d(true);
                for at in 0..res.num_atoms() as usize {
                    if 3 * at + 2 < coords.len() {
                        conf.set_atom_pos(
                            at as u32,
                            Point3D::new(coords[3 * at], coords[3 * at + 1], coords[3 * at + 2]),
                        )?;
                    }
                }
                res.add_conformer(conf, false)?;
            } else {
                let Some(end_p) = in_line.rfind('>') else {
                    return Err(RDKitError::FileParse(format!(
                        "no end tag found for property{prop_name}"
                    )));
                };
                let s0 = start_p.min(in_line.len());
                let prop_val = in_line[s0..end_p.max(s0)].to_string();
                res.set_prop(&prop_name, prop_val.clone(), false)?;
                if prop_name == self.params.name_record {
                    res.set_prop(cp::_Name, prop_val, false)?;
                }
            }
            in_line = lossy(&s.getline_raw());
        }
        Ok(Some(res))
    }

    fn next_in(&mut self, s: &mut LineStream<'_>) -> Result<Option<RWMol>> {
        // set the stream to the appropriate position
        seekg_opt(s, self.d_molpos[self.d_last]);
        // finally if we reached the end of the file set end to be true
        if s.eof() {
            self.df_end = true;
            self.d_len = self.d_molpos.len() as i64;
            return Ok(None);
        }
        // start by finding the $SMI element (we're assuming that this starts
        // the block)
        self.d_line += 1;
        let mut temp = lossy(&s.getline_raw());
        while !temp.starts_with("$SMI<") && !s.eof() && !s.fail() {
            self.d_line += 1;
            temp = lossy(&s.getline_raw());
        }
        let mut res = None;
        if temp.starts_with("$SMI<") {
            match self.parse_mol(s, temp) {
                Ok(m) => res = m,
                Err(se @ RDKitError::Sanitize(_)) => {
                    // We couldn't sanitize a molecule we got - write out an
                    // error message and move to the end of the record
                    rd_log_error!(
                        "ERROR: Could not sanitize molecule ending on line {}",
                        self.d_line
                    );
                    rd_log_error!("ERROR: {se}\n");
                    let mut t = String::new();
                    while !s.eof() && !s.fail() && !t.starts_with('|') {
                        self.d_line += 1;
                        t = lossy(&s.getline_raw());
                    }
                }
                Err(e) => return Err(e),
            }
        }
        self.d_last += 1;
        if self.d_last >= self.d_molpos.len() {
            self.d_molpos.push(tellg(s));
        }
        self.check_for_end(s);
        Ok(res)
    }

    fn move_to_in(&mut self, s: &mut LineStream<'_>, idx: usize) -> Result<()> {
        // seekg() is called for all idx values and earlier calls to next()
        // may have put the stream into a bad state
        clear(s);
        // move until we hit the desired idx
        if idx < self.d_molpos.len() {
            seekg_opt(s, self.d_molpos[idx]);
            self.d_last = idx;
        } else {
            self.d_last = self.d_molpos.len() - 1;
            seekg_opt(s, *self.d_molpos.last().expect("a position"));
            while self.d_last < idx && !s.eof() && !s.fail() {
                self.d_line += 1;
                let temp = lossy(&s.getline_raw());
                if temp.starts_with('|') {
                    self.d_molpos.push(tellg(s));
                    self.d_last += 1;
                }
            }
            // if we reached end of file without reaching "idx" we have an
            // index error
            if s.eof() {
                self.d_len = self.d_molpos.len() as i64;
                return Err(RDKitError::FileParse(format!(
                    "ERROR: Index error (idx = {idx}) :  we do no have enough molecule blocks"
                )));
            }
        }
        Ok(())
    }

    /// `moveTo(idx)`
    pub fn move_to(&mut self, idx: usize) -> Result<()> {
        self.with_stream(|this, s| this.move_to_in(s, idx))?
    }

    /// `next()`: the next molecule, `None` when its record cannot be parsed.
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

    /// `getItemText(idx)`: the text of record `idx`.
    pub fn get_item_text(&mut self, idx: usize) -> Result<String> {
        self.with_stream(|this, s| {
            let holder = this.d_last;
            this.move_to_in(s, idx)?;
            let beg = this.d_molpos[idx];
            let end_holder = this.df_end;
            let end = match this.move_to_in(s, idx + 1) {
                Ok(()) => this.d_molpos[idx + 1],
                Err(RDKitError::FileParse(_)) => {
                    clear(s);
                    Some(s.data().len())
                }
                Err(e) => return Err(e),
            };
            this.d_last = holder;
            this.df_end = end_holder;
            let (b, e) = (beg.unwrap_or(0), end.unwrap_or(s.data().len()));
            Ok(lossy(&s.data()[b.min(e)..e]))
        })?
    }

    /// `length()`: the number of records (reads the whole input once).
    pub fn length(&mut self) -> Result<usize> {
        self.with_stream(|this, s| {
            if this.d_len > 0 {
                return this.d_len as usize;
            }
            this.d_len = this.d_molpos.len() as i64;
            seekg_opt(s, *this.d_molpos.last().expect("a position"));
            let _ = s.getline_raw();
            while this.advance_to_next_record(s) {
                this.d_molpos.push(tellg(s));
                this.d_len += 1;
                let _ = s.getline_raw();
            }
            // now remember to set the stream to the last position we want to
            // read
            clear(s);
            seekg_opt(s, this.d_molpos[this.d_last]);
            this.d_len as usize
        })
    }
}
