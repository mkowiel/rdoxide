//! `SynthonSpace.cpp`: the collection of synthon sets and the synthon pool.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use super::params::{CONNECTOR_SYMBOLS, MAX_CONNECTOR_NUM};
use super::stream::{Reader, Writer};
use super::synthon::Synthon;
use super::synthon_set::SynthonSet;
use crate::{RDKitError, Result};

const VERSION_MAJOR: i32 = 3;
const VERSION_MINOR: i32 = 2;
const ENDIAN_ID: i32 = 0xa100f;

#[derive(Default)]
pub struct SynthonSpace {
    file_name: String,
    /// The synthons, sorted by SMILES.
    pub(crate) pool: Vec<Synthon>,
    /// The synthon sets, sorted by reaction id.
    pub(crate) reactions: Vec<SynthonSet>,
    /// For each pool entry, the indices of the reactions it is in.
    synthon_reactions: Vec<Vec<usize>>,
    fp_type: String,
    num_conformers: u32,
    file_version: u32,
    num_products: u64,
    max_num_synthons: usize,
    has_ring_former: bool,
}

/// `std::regex_token_iterator` with submatch -1 for a pattern that is a run
/// of the characters accepted by `is_sep` (`collapse`: `+` quantified).
fn split_line(s: &str, is_sep: impl Fn(char) -> bool, collapse: bool) -> Vec<String> {
    let mut parts = Vec::new();
    let mut cur = String::new();
    let mut in_sep = false;
    let mut chars = s.chars().peekable();
    // a leading separator gives an empty first token
    while let Some(c) = chars.next() {
        if is_sep(c) {
            if !(collapse && in_sep) {
                parts.push(std::mem::take(&mut cur));
            }
            in_sep = true;
            let _ = chars.peek();
        } else {
            cur.push(c);
            in_sep = false;
        }
    }
    if !cur.is_empty() {
        parts.push(cur);
    }
    parts
}

fn fix_connectors(smiles: &mut String) {
    for (i, sym) in CONNECTOR_SYMBOLS.iter().enumerate().take(MAX_CONNECTOR_NUM) {
        *smiles = smiles.replace(sym, &format!("[{}*]", i + 1));
    }
}

fn deduce_format(line: &str) -> i32 {
    let first_line_opts: [&[&str]; 3] = [
        &["SMILES", "synton_id", "synton#", "reaction_id"],
        &["SMILES", "synton_id", "synton#", "reaction_id", "release"],
        &["SMILES", "synton_id", "synton_role", "reaction_id"],
    ];
    let matches = |parts: &[String], opt: &[&str]| {
        parts.len() == opt.len() && parts.iter().zip(opt).all(|(a, b)| a == b)
    };
    let parts = split_line(line, |c| c == '\t', false);
    for (i, o) in first_line_opts.iter().enumerate() {
        if matches(&parts, o) {
            return i as i32 + 6;
        }
    }
    let parts = split_line(line, char::is_whitespace, true);
    for (i, o) in first_line_opts.iter().enumerate() {
        if matches(&parts, o) {
            return i as i32;
        }
    }
    let parts = split_line(line, |c| c == ',', true);
    for (i, o) in first_line_opts.iter().enumerate() {
        if matches(&parts, o) {
            return (i + first_line_opts.len()) as i32;
        }
    }
    -1
}

fn get_synthon_num(format: i32, synthon: &str) -> Result<i32> {
    let bad = |_| RDKitError::Runtime("stoi".into());
    Ok(match format {
        0 | 1 | 3 | 4 | 6 | 7 => synthon.trim().parse::<i32>().map_err(bad)?,
        2 | 5 | 8 => synthon
            .get(7..)
            .ok_or_else(|| RDKitError::Runtime("basic_string::substr".into()))?
            .parse::<i32>()
            .map_err(bad)?,
        _ => i32::MAX,
    })
}

pub fn formatted_integer_string(value: i64) -> String {
    let s = value.unsigned_abs().to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(' ');
        }
        out.push(c);
    }
    if value < 0 {
        out.insert(0, '-');
    }
    out
}

impl SynthonSpace {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn num_reactions(&self) -> usize {
        self.reactions.len()
    }
    pub fn reaction_names(&self) -> Vec<String> {
        self.reactions.iter().map(|r| r.id().to_string()).collect()
    }
    pub fn reaction(&self, name: &str) -> Result<&SynthonSet> {
        self.reactions
            .binary_search_by(|r| r.id().cmp(name))
            .map(|i| &self.reactions[i])
            .map_err(|_| {
                RDKitError::Runtime(format!("Could not find synthon set for reaction {name}"))
            })
    }
    pub fn num_products(&self) -> u64 {
        self.num_products
    }
    pub fn num_synthons(&self) -> usize {
        self.pool.len()
    }
    pub fn input_file_name(&self) -> &str {
        &self.file_name
    }
    pub fn has_ring_former(&self) -> bool {
        self.has_ring_former
    }
    pub fn max_num_synthons(&self) -> usize {
        self.max_num_synthons
    }
    pub fn max_num_connectors(&self) -> u32 {
        self.reactions
            .iter()
            .map(|r| r.connectors().count_ones())
            .max()
            .unwrap_or(0)
    }
    pub fn has_fingerprints(&self) -> bool {
        !self.fp_type.is_empty()
    }
    pub fn fp_type(&self) -> &str {
        &self.fp_type
    }
    pub fn num_conformers(&self) -> u32 {
        self.num_conformers
    }
    pub fn pool(&self) -> &[Synthon] {
        &self.pool
    }
    pub fn pattern_fp_size(&self) -> Result<u32> {
        use crate::datastructs::bit_vect::BitVect;
        if self.reactions.is_empty() {
            return Err(RDKitError::precondition("No synthon sets available."));
        }
        for r in &self.reactions {
            if let Some(s) = r.synthons().first().and_then(|s| s.first()) {
                return Ok(self.pool[s.1].patt_fp().num_bits());
            }
        }
        Err(RDKitError::Runtime(
            "Could not find pattern fingerprint for any synthon.".into(),
        ))
    }
    pub fn fp_size(&self) -> Result<u32> {
        use crate::datastructs::bit_vect::BitVect;
        if self.reactions.is_empty() {
            return Err(RDKitError::precondition("No synthon sets available."));
        }
        for r in &self.reactions {
            if let Some(s) = r.synthons().first().and_then(|s| s.first()) {
                return self.pool[s.1]
                    .fp()
                    .map(|f| f.num_bits())
                    .ok_or_else(|| RDKitError::Runtime("Could not find fingerprint.".into()));
            }
        }
        Err(RDKitError::Runtime(
            "Could not find fingerprint for any synthon.".into(),
        ))
    }
    pub fn num_synthons_with_shapes(&self) -> u64 {
        self.pool.iter().filter(|s| s.has_shapes()).count() as u64
    }

    fn pool_index(&self, smiles: &str) -> Option<usize> {
        self.pool.binary_search_by(|s| s.smiles().cmp(smiles)).ok()
    }

    /// `readTextFile`
    pub fn read_text_file(&mut self, in_filename: &str) -> Result<()> {
        self.file_name = in_filename.to_string();
        let data = std::fs::read(in_filename)
            .map_err(|_| RDKitError::Runtime(format!("Couldn't open file {in_filename}")))?;
        self.read_stream(&String::from_utf8_lossy(&data))
    }

    /// `readStream`
    pub fn read_stream(&mut self, text: &str) -> Result<()> {
        let mut format = -1;
        // (reaction, synthon set number, synthon id, SMILES)
        let mut entries: Vec<(String, i32, String, String)> = Vec::new();
        let mut lines: Vec<&str> = text.split('\n').collect();
        if lines.last().is_some_and(|l| l.is_empty()) {
            lines.pop();
        }
        for (this_line, raw) in (1..).zip(lines) {
            let next_line = raw.strip_suffix('\r').unwrap_or(raw);
            if next_line.is_empty() || next_line.starts_with('#') {
                continue;
            }
            if format == -1 {
                format = deduce_format(next_line);
                if format == -1 {
                    return Err(RDKitError::Runtime(format!(
                        "Bad format for SynthonSpace file {}",
                        self.file_name
                    )));
                }
                continue;
            }
            // (upstream leaves format 3 without a splitter)
            let mut synthon = if format < 3 {
                split_line(next_line, char::is_whitespace, true)
            } else if format > 3 && format < 6 {
                split_line(next_line, |c| c == ',', true)
            } else if format > 5 {
                split_line(next_line, |c| c == '\t', true)
            } else {
                Vec::new()
            };
            if synthon.len() < 4 {
                return Err(RDKitError::Runtime(format!(
                    "Bad format for SynthonSpace file {} on line {}",
                    self.file_name, this_line
                )));
            }
            fix_connectors(&mut synthon[0]);
            let num = get_synthon_num(format, &synthon[2])?;
            entries.push((
                synthon[3].clone(),
                num,
                synthon[1].clone(),
                synthon[0].clone(),
            ));
        }

        // The synthon pool, sorted by SMILES.
        let smiles: BTreeSet<&str> = entries.iter().map(|e| e.3.as_str()).collect();
        self.pool.clear();
        for smi in smiles {
            self.pool.push(Synthon::new(smi)?);
        }
        let mut reactions: BTreeMap<String, SynthonSet> = BTreeMap::new();
        for (rxn, num, id, smi) in &entries {
            let set = reactions
                .entry(rxn.clone())
                .or_insert_with(|| SynthonSet::new(rxn));
            let idx = self.pool_index(smi).expect("synthon in pool");
            set.add_synthon(*num as usize, idx, id);
        }
        self.reactions = reactions.into_values().collect();

        // Do some final processing.
        self.num_products = 0;
        self.max_num_synthons = 0;
        for r in 0..self.reactions.len() {
            let rxn = &mut self.reactions[r];
            rxn.remove_empty_synthon_sets();
            rxn.make_synthon_search_mols(&mut self.pool)?;
            rxn.build_connector_regions(&self.pool)?;
            rxn.assign_connectors_used(&self.pool);
            rxn.assess_ring_formers();
            if rxn.num_ring_formers() > 0 {
                self.has_ring_former = true;
            }
            self.num_products += rxn.num_products();
            self.max_num_synthons = self.max_num_synthons.max(rxn.synthons().len());
            rxn.initialize_search_orders(&self.pool);
        }
        self.fill_synthon_reactions();
        Ok(())
    }

    fn fill_synthon_reactions(&mut self) {
        self.synthon_reactions = vec![Vec::new(); self.pool.len()];
        for (r, rxn) in self.reactions.iter().enumerate() {
            let n_sets = rxn.synthons().len() as u32;
            for set in rxn.synthons() {
                for (_, idx) in set {
                    self.synthon_reactions[*idx].push(r);
                    self.pool[*idx].update_max_synthon_set_size(n_sets);
                }
            }
        }
    }

    /// `writeDBFile`
    pub fn write_db_file(&self, out_filename: &str) -> Result<()> {
        let mut w = Writer::new();
        w.i32(ENDIAN_ID);
        w.i32(VERSION_MAJOR);
        w.i32(VERSION_MINOR);
        w.bool(self.has_fingerprints());
        if self.has_fingerprints() {
            w.string(self.fp_type.as_bytes());
        }
        w.u32(self.num_conformers);
        w.u64(self.pool.len() as u64);
        w.u64(self.reactions.len() as u64);
        w.u64(self.num_products);
        let header_len = w.buf.len() + 8 * (self.pool.len() + self.reactions.len());

        let mut body = Writer::new();
        let mut synthon_pos = Vec::with_capacity(self.pool.len());
        for s in &self.pool {
            synthon_pos.push((header_len + body.buf.len()) as u64);
            s.write_to_db_stream(&mut body)?;
        }
        let mut reaction_pos = Vec::with_capacity(self.reactions.len());
        for r in &self.reactions {
            reaction_pos.push((header_len + body.buf.len()) as u64);
            r.write_to_db_stream(&self.pool, &mut body)?;
        }
        for p in synthon_pos.iter().chain(reaction_pos.iter()) {
            w.u64(*p);
        }
        w.buf.extend_from_slice(&body.buf);
        std::fs::write(out_filename, &w.buf)
            .map_err(|e| RDKitError::Runtime(format!("Couldn't write file {out_filename}: {e}")))
    }

    /// `readDBFile`
    pub fn read_db_file(&mut self, in_filename: &str, _num_threads: i32) -> Result<()> {
        self.file_name = in_filename.to_string();
        let data = std::fs::read(in_filename)
            .map_err(|_| RDKitError::Runtime(format!("Couldn't open file {in_filename}")))?;
        let mut r = Reader::new(&data);
        let endian_test = r.i32()?;
        if endian_test != ENDIAN_ID {
            return Err(RDKitError::Runtime(format!(
                "Endianness mismatch in SynthonSpace file {in_filename}"
            )));
        }
        let major = r.i32()?;
        let minor = r.i32()?;
        if major > VERSION_MAJOR || (major == VERSION_MAJOR && minor > VERSION_MINOR) {
            crate::rd_log_warning!(
                "Deserializing from a version number ({major}.{minor})that is higher than our version ({VERSION_MAJOR}.{VERSION_MINOR}).\nThis probably won't work."
            );
        }
        if major > 1000 || minor > 100 {
            return Err(RDKitError::Runtime("unreasonable version numbers".into()));
        }
        let version = (1000 * major + minor * 10) as u32;
        self.file_version = version;
        if version < 3000 {
            return Err(RDKitError::Runtime(
                "This binary file version is no longer supported.  Please re-build with a recent version of the RDKit."
                    .into(),
            ));
        }
        if r.bool()? {
            self.fp_type = r.string()?;
        }
        if version > 3010 {
            self.num_conformers = r.u32()?;
        }
        let num_synthons = r.u64()? as usize;
        let num_reactions = r.u64()? as usize;
        self.num_products = r.u64()?;
        let mut synthon_pos = Vec::with_capacity(num_synthons + 1);
        for _ in 0..num_synthons {
            synthon_pos.push(r.u64()? as usize);
        }
        let mut reaction_pos = Vec::with_capacity(num_reactions + 1);
        for _ in 0..num_reactions {
            reaction_pos.push(r.u64()? as usize);
        }
        if num_reactions == 0 {
            return Err(RDKitError::Runtime(format!(
                "Database file {in_filename} appears corrupted."
            )));
        }
        synthon_pos.push(reaction_pos[0]);
        reaction_pos.push(data.len());

        self.pool.clear();
        for i in 0..num_synthons {
            let mut sr = Reader::new(&data[synthon_pos[i]..synthon_pos[i + 1]]);
            self.pool
                .push(Synthon::read_from_db_stream(&mut sr, version)?);
        }
        if !self.pool.windows(2).all(|w| w[0].smiles() <= w[1].smiles()) {
            self.pool.sort_by(|a, b| a.smiles().cmp(b.smiles()));
        }
        self.reactions.clear();
        for i in 0..num_reactions {
            let mut rr = Reader::new(&data[reaction_pos[i]..reaction_pos[i + 1]]);
            let pool = &self.pool;
            let lookup = |smi: &str| pool.binary_search_by(|s| s.smiles().cmp(smi)).ok();
            self.reactions.push(SynthonSet::read_from_db_stream(
                &mut rr,
                &lookup,
                in_filename,
                version,
            )?);
        }
        if !self.reactions.windows(2).all(|w| w[0].id() <= w[1].id()) {
            self.reactions.sort_by(|a, b| a.id().cmp(b.id()));
        }
        self.has_ring_former = false;
        self.max_num_synthons = 0;
        for rxn in &mut self.reactions {
            rxn.assess_ring_formers();
            if rxn.num_ring_formers() > 0 {
                self.has_ring_former = true;
            }
            self.max_num_synthons = self.max_num_synthons.max(rxn.synthons().len());
        }
        self.fill_synthon_reactions();
        crate::rd_log_info!(
            "Number of synthons : {} of which {} have shapes.",
            self.pool.len(),
            self.num_synthons_with_shapes()
        );
        Ok(())
    }

    /// `summarise`
    pub fn summarise(&self) -> String {
        let mut os = String::new();
        let _ = writeln!(os, "Read from file {}", self.file_name);
        let _ = writeln!(os, "Number of reactions : {}", self.reactions.len());
        let mut synth_counts = [0usize; MAX_CONNECTOR_NUM + 1];
        for rxn in &self.reactions {
            let _ = writeln!(os, "Reaction name {}", rxn.id());
            for (i, s) in rxn.synthons().iter().enumerate() {
                let _ = writeln!(os, "  Synthon set {i} has {} synthons", s.len());
            }
            if rxn.synthons().len() <= MAX_CONNECTOR_NUM {
                synth_counts[rxn.synthons().len()] += 1;
            }
        }
        let _ = writeln!(
            os,
            "Approximate number of molecules in SynthonSpace : {}",
            formatted_integer_string(self.num_products as i64)
        );
        let _ = writeln!(os, "Number of unique synthons : {}", self.pool.len());
        for (i, c) in synth_counts.iter().enumerate().take(MAX_CONNECTOR_NUM) {
            if *c > 0 {
                let _ = writeln!(os, "Number of {i} molecule reactions : {c}");
            }
        }
        os
    }

    /// `writeEnumeratedFile`
    pub fn write_enumerated_file(&self, out_filename: &str) -> Result<()> {
        crate::rd_log_warning!(
            "Writing the enumerated file may take quite some time and result in a large file."
        );
        let mut out = String::new();
        for r in &self.reactions {
            r.enumerate_to_stream(&self.pool, &mut out);
        }
        std::fs::write(out_filename, out)
            .map_err(|e| RDKitError::Runtime(format!("Couldn't write file {out_filename}: {e}")))
    }
}

/// `convertTextToDBFile` (without fingerprints or shapes)
pub fn convert_text_to_db_file(in_filename: &str, out_filename: &str) -> Result<()> {
    let mut space = SynthonSpace::new();
    space.read_text_file(in_filename)?;
    space.write_db_file(out_filename)
}
