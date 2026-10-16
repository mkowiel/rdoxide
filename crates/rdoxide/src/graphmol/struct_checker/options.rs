//! Options and rule data of the structure checker.
//!
//! Ported from RDKit `Code/GraphMol/StructChecker/StructCheckerOptions.cpp`
//! (2026.09.1).
//!
//! Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;
use std::sync::Arc;

use super::aa_data::{DEFAULT_ACIDIC_ATOMS, DEFAULT_AUGMENTED_ATOM_TRANSFORMS, DEFAULT_GOOD_ATOMS};
use super::{
    AABondType, AATopology, ANY_CHARGE, AugmentedAtom, IncEntry, Ligand, PathEntry, RadicalType,
};
use crate::graphmol::file_parsers::mol_file_parser::{
    MolFileParserParams, mol_from_mol_data_stream,
};
use crate::graphmol::file_parsers::sd_mol_supplier::SDMolSupplier;
use crate::graphmol::file_parsers::utils::LineStream;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::rdgeneral::json::Json;
use crate::{RDKitError, Result};

/// The options of the structure checker together with its rule data.
#[derive(Debug, Clone)]
pub struct StructCheckerOptions {
    pub acidity_limit: f64,
    pub remove_minor_fragments: bool,
    pub desired_charge: i32,
    pub check_collisions: bool,
    pub collision_limit_percent: i32,
    pub max_mol_size: u32,
    pub convert_s_text: bool,
    pub squeeze_identifiers: bool,
    pub strip_zeros: bool,
    pub check_stereo: bool,
    pub convert_atom_texts: bool,
    pub groups_to_s_groups: bool,
    pub verbose: bool,

    pub augmented_atom_pairs: Vec<(AugmentedAtom, AugmentedAtom)>,
    pub acidic_atoms: Vec<AugmentedAtom>,
    pub good_atoms: Vec<AugmentedAtom>,
    pub patterns: Vec<Arc<ROMol>>,
    pub rotate_patterns: Vec<Arc<ROMol>>,
    pub stereo_patterns: Vec<Arc<ROMol>>,
    pub from_tautomer: Vec<Arc<ROMol>>,
    pub to_tautomer: Vec<Arc<ROMol>>,

    /// `elneg_table[0].value`
    pub elneg0: f64,
    /// Atomic number -> electronegativity
    pub elneg_table: BTreeMap<u32, f64>,
    pub atom_acidity: Vec<IncEntry>,
    pub charge_inc_table: Vec<IncEntry>,
    /// The coefficients of the function that stretches the preliminary pKa
    /// values to the predictions: `pKa = 7 + (pKa'-7)*beta + ((pKa'-7)*alpha)^3`.
    pub alpha: f64,
    pub beta: f64,
    pub alpha_path_table: Vec<PathEntry>,
    pub beta_path_table: Vec<PathEntry>,
}

impl Default for StructCheckerOptions {
    /// The default options with the default augmented atom data.
    fn default() -> Self {
        let parse = |s: &str| string_to_augmented_atom(s).expect("INTERNAL ERROR in default data");
        let good = DEFAULT_GOOD_ATOMS.iter().map(|s| parse(s)).collect();
        let acidic = DEFAULT_ACIDIC_ATOMS.iter().map(|s| parse(s)).collect();
        let pairs = DEFAULT_AUGMENTED_ATOM_TRANSFORMS
            .iter()
            .map(|(from, to)| (parse(from), parse(to)))
            .collect();
        StructCheckerOptions {
            acidity_limit: 0.0,
            remove_minor_fragments: false,
            desired_charge: 0,
            check_collisions: false,
            collision_limit_percent: 15,
            max_mol_size: 255,
            convert_s_text: false,
            squeeze_identifiers: false,
            strip_zeros: false,
            check_stereo: false,
            convert_atom_texts: false,
            groups_to_s_groups: false,
            verbose: false,
            augmented_atom_pairs: pairs,
            acidic_atoms: acidic,
            good_atoms: good,
            patterns: Vec::new(),
            rotate_patterns: Vec::new(),
            stereo_patterns: Vec::new(),
            from_tautomer: Vec::new(),
            to_tautomer: Vec::new(),
            elneg0: 0.0,
            elneg_table: BTreeMap::new(),
            atom_acidity: Vec::new(),
            charge_inc_table: Vec::new(),
            alpha: 0.0,
            beta: 0.0,
            alpha_path_table: Vec::new(),
            beta_path_table: Vec::new(),
        }
    }
}

/// How a `boost::property_tree` value converts: JSON scalars are strings.
fn ptree_string(v: &Json) -> Option<String> {
    Some(match v {
        Json::Null => "null".into(),
        Json::Bool(b) => b.to_string(),
        Json::Int(i) => i.to_string(),
        Json::UInt(u) => u.to_string(),
        Json::Double(d) => crate::rdgeneral::strings::format_g(*d, 17),
        Json::String(s) => s.clone(),
        Json::Array(_) | Json::Object(_) => return None,
    })
}

fn ptree_f64(o: &crate::rdgeneral::json::JsonObject, key: &str, default: f64) -> Result<f64> {
    let Some(v) = o.get(key) else {
        return Ok(default);
    };
    ptree_string(v)
        .and_then(|s| s.trim().parse::<f64>().ok())
        .ok_or_else(|| bad_data(key))
}

fn ptree_i32(o: &crate::rdgeneral::json::JsonObject, key: &str, default: i32) -> Result<i32> {
    let Some(v) = o.get(key) else {
        return Ok(default);
    };
    ptree_string(v)
        .and_then(|s| s.trim().parse::<i32>().ok())
        .ok_or_else(|| bad_data(key))
}

fn ptree_bool(o: &crate::rdgeneral::json::JsonObject, key: &str, default: bool) -> Result<bool> {
    let Some(v) = o.get(key) else {
        return Ok(default);
    };
    match ptree_string(v).as_deref().map(str::trim) {
        Some("true") | Some("1") => Ok(true),
        Some("false") | Some("0") => Ok(false),
        _ => Err(bad_data(key)),
    }
}

fn bad_data(key: &str) -> RDKitError {
    RDKitError::Runtime(format!("conversion of data to type failed: {key}"))
}

/// `parseOptionsJSON`: resets `op` to the defaults and reads the scalar
/// options from the JSON object; `false` (after logging) on any error.
pub fn parse_options_json(json: &str, op: &mut StructCheckerOptions) -> bool {
    if json.is_empty() {
        return false;
    }
    let parsed = Json::parse(json).and_then(|root| {
        let mut o = StructCheckerOptions::default();
        if let Json::Object(obj) = &root {
            o.acidity_limit = ptree_f64(obj, "AcidityLimit", o.acidity_limit)?;
            o.remove_minor_fragments =
                ptree_bool(obj, "RemoveMinorFragments", o.remove_minor_fragments)?;
            o.desired_charge = ptree_i32(obj, "DesiredCharge", o.desired_charge)?;
            o.check_collisions = ptree_bool(obj, "CheckCollisions", o.check_collisions)?;
            o.collision_limit_percent =
                ptree_i32(obj, "CollisionLimitPercent", o.collision_limit_percent)?;
            o.max_mol_size = ptree_i32(obj, "MaxMolSize", o.max_mol_size as i32)? as u32;
            o.convert_s_text = ptree_bool(obj, "ConvertSText", o.convert_s_text)?;
            o.squeeze_identifiers = ptree_bool(obj, "SqueezeIdentifiers", o.squeeze_identifiers)?;
            o.strip_zeros = ptree_bool(obj, "StripZeros", o.strip_zeros)?;
            o.check_stereo = ptree_bool(obj, "CheckStereo", o.check_stereo)?;
            o.convert_atom_texts = ptree_bool(obj, "ConvertAtomTexts", o.convert_atom_texts)?;
            o.groups_to_s_groups = ptree_bool(obj, "GroupsToSGroups", o.groups_to_s_groups)?;
            o.verbose = ptree_bool(obj, "Verbose", o.verbose)?;
        }
        Ok(o)
    });
    match parsed {
        Ok(o) => {
            *op = o;
            true
        }
        Err(e) => {
            crate::rd_log_error!("JSON StructureCheck Options:{e}\n");
            false
        }
    }
}

/// `loadOptionsFromFiles`: empty file names are skipped.
pub fn load_options_from_files(
    op: &mut StructCheckerOptions,
    augmented_atom_translations_file: &str,
    pattern_file: &str,
    rotate_pattern_file: &str,
    stereo_pattern_file: &str,
    tautomer_file: &str,
) -> Result<bool> {
    let mut res = true;
    if !augmented_atom_translations_file.is_empty() {
        res &= op.load_augmented_atom_translations(augmented_atom_translations_file)?;
    }
    if !pattern_file.is_empty() {
        res &= op.load_patterns(pattern_file);
    }
    if !rotate_pattern_file.is_empty() {
        res &= op.load_rotate_patterns(rotate_pattern_file);
    }
    if !stereo_pattern_file.is_empty() {
        res &= op.load_stereo_patterns(stereo_pattern_file);
    }
    if !tautomer_file.is_empty() {
        res &= op.load_tautomer_data(tautomer_file);
    }
    Ok(res)
}

//=====================================================================
// Augmented atom strings

/// Bond symbols, ordered like `AABondType`.
const BOND_TO_STRING: [&str; 9] = ["?", "-", "=", "#", "~", "-=", "-~", "=~", "*"];

fn bond_type_from_index(i: usize) -> AABondType {
    match i {
        0 => AABondType::None,
        1 => AABondType::Single,
        2 => AABondType::Double,
        3 => AABondType::Triple,
        4 => AABondType::Aromatic,
        5 => AABondType::SingleDouble,
        6 => AABondType::SingleAromatic,
        7 => AABondType::DoubleAromatic,
        _ => AABondType::Any,
    }
}

/// A cursor over a NUL terminated C string.
struct Cursor<'a> {
    s: &'a [u8],
    p: usize,
}

impl Cursor<'_> {
    fn at(&self, off: usize) -> u8 {
        self.s.get(self.p + off).copied().unwrap_or(0)
    }
    fn rest(&self) -> String {
        String::from_utf8_lossy(self.s.get(self.p..).unwrap_or(&[])).into_owned()
    }
    /// The atom symbol: letters, digits, ',' and '#'.
    fn symbol(&mut self) -> String {
        let mut i = 0;
        while self.at(i).is_ascii_alphanumeric() || self.at(i) == b',' || self.at(i) == b'#' {
            i += 1;
        }
        let sym = String::from_utf8_lossy(&self.s[self.p..self.p + i]).into_owned();
        self.p += i;
        sym
    }
    /// A charge descriptor; `None` after logging a syntax error.
    fn charge(&mut self) -> Option<i32> {
        if self.at(0) == b'+' && self.at(1) == b'-' {
            self.p += 2;
            Some(ANY_CHARGE)
        } else if self.at(0) == b'-' || self.at(0) == b'+' {
            let sign = if self.at(0) == b'-' { -1 } else { 1 };
            self.p += 1;
            if (b'1'..=b'7').contains(&self.at(0)) {
                let c = sign * i32::from(self.at(0) - b'0');
                self.p += 1;
                Some(c)
            } else {
                self.p -= 1;
                crate::rd_log_error!("syntax error '{}'\n", self.rest());
                None
            }
        } else {
            Some(0)
        }
    }
    fn radical(&mut self) -> RadicalType {
        // a digit after the symbol is a degree specification
        let r = match self.at(0) {
            b'|' => RadicalType::Singlet,
            b'.' => RadicalType::Doublet,
            b':' => RadicalType::Triplet,
            _ => return RadicalType::None,
        };
        if self.at(1).is_ascii_digit() {
            return RadicalType::None;
        }
        self.p += 1;
        r
    }
}

/// `StringToAugmentedAtom`. The syntax is
/// `['@'|'!@'] <Atom Symbol> [<Charge>] {'(' <Bond Symbol> <Atom Symbol> [<Charge>] ')'}`;
/// '@' means the central atom is in a ring, '!@' that it is not.
/// Syntax errors are logged and give `None`.
pub fn string_to_augmented_atom(text: &str) -> Option<AugmentedAtom> {
    let mut aa = AugmentedAtom {
        short_name: text.to_string(),
        ..Default::default()
    };
    let mut c = Cursor {
        s: text.as_bytes(),
        p: 0,
    };
    if c.at(0) == b'@' {
        aa.topology = AATopology::Ring;
        c.p += 1;
    } else if c.at(0) == b'!' && c.at(1) == b'@' {
        aa.topology = AATopology::Chain;
        c.p += 2;
    }

    if !c.at(0).is_ascii_alphabetic() {
        crate::rd_log_error!("atom symbol is required '{}'\n", c.rest());
        return None;
    }
    aa.atom_symbol = c.symbol();
    aa.charge = c.charge()?;
    aa.radical = c.radical();

    while c.at(0) == b'(' {
        let mut lig = Ligand::default();
        c.p += 1;
        // the bond type
        let n = if c.at(1).is_ascii_alphabetic() { 1 } else { 2 };
        let found = BOND_TO_STRING.iter().position(|b| {
            let lit: Vec<u8> = b.bytes().chain(std::iter::once(0)).collect();
            (0..n).all(|k| c.at(k) == lit.get(k).copied().unwrap_or(0))
        });
        let Some(i) = found else {
            crate::rd_log_error!(
                "syntax error '{}'\nthere must be a bond type symbol.\n",
                c.rest()
            );
            return None;
        };
        lig.bond_type = bond_type_from_index(i);
        c.p += BOND_TO_STRING[i].len();
        if !c.at(0).is_ascii_alphabetic() {
            crate::rd_log_error!(
                "syntax error '{}'\nthere must be an atom symbol.\n",
                c.rest()
            );
            return None;
        }
        lig.atom_symbol = c.symbol();
        lig.charge = c.charge()?;
        lig.radical = c.radical();
        // substitution count descriptor
        if c.at(0) == b':' && c.at(1).is_ascii_digit() {
            lig.substitution_count = u32::from(c.at(1) - b'0');
            c.p += 2;
        } else {
            lig.substitution_count = 0;
        }
        if c.at(0) == b')' {
            c.p += 1;
        } else {
            crate::rd_log_error!("unclosed ( '{}'\n", c.rest());
            return None;
        }
        aa.ligands.push(lig);
    }

    if c.at(0) != 0 {
        crate::rd_log_error!("unexpected symbols after scan '{}'\n", c.rest());
        return None;
    }
    Some(aa)
}

/// A C `char >= ' '` (the C `char` is signed, so bytes >= 128 are not).
fn printable(b: u8) -> bool {
    (b as i8) >= b' ' as i8
}

/// The text between the first two '"' of `line` starting at `from`:
/// `(text, index after the closing quote)`.
fn quoted(line: &[u8], from: usize) -> Option<(String, usize)> {
    let mut s = from;
    while s < line.len() && printable(line[s]) && line[s] != b'"' {
        s += 1;
    }
    if s >= line.len() || line[s] != b'"' {
        return None;
    }
    s += 1;
    let mut len = 0;
    while s + len < line.len() && printable(line[s + len]) && line[s + len] != b'"' {
        len += 1;
    }
    Some((
        String::from_utf8_lossy(&line[s..s + len]).into_owned(),
        s + len + 1,
    ))
}

/// The lines of a file like `fgets` gives them (at most 1023 bytes each).
fn fgets_lines(data: &[u8]) -> Vec<&[u8]> {
    let mut lines = Vec::new();
    let mut pos = 0;
    while pos < data.len() {
        let mut end = pos;
        while end < data.len() && end - pos < 1023 {
            end += 1;
            if data[end - 1] == b'\n' {
                break;
            }
        }
        lines.push(&data[pos..end]);
        pos = end;
    }
    lines
}

/// `%d` of `sscanf`/`fscanf`: skips whitespace, then an optionally signed
/// integer. Returns the value and the number of bytes consumed.
fn scan_int(s: &[u8]) -> Option<(i64, usize)> {
    let mut i = 0;
    while i < s.len() && s[i].is_ascii_whitespace() {
        i += 1;
    }
    let start = i;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') {
        i += 1;
    }
    let digits = i;
    while i < s.len() && s[i].is_ascii_digit() {
        i += 1;
    }
    if i == digits {
        return None;
    }
    std::str::from_utf8(&s[start..i])
        .ok()?
        .parse()
        .ok()
        .map(|v| (v, i))
}

fn read_file(path: &str) -> Option<Vec<u8>> {
    match std::fs::read(path) {
        Ok(d) => Some(d),
        Err(e) => {
            crate::rd_log_error!("could not open file '{path}'\n{e}\n");
            None
        }
    }
}

/// `ReadAugmentedAtoms`: loader of the `*.chk` and `*.aci` files. Only the
/// part of a line between the two '"' is used.
fn read_augmented_atoms(path: &str) -> Option<Vec<AugmentedAtom>> {
    let data = read_file(path)?;
    let lines = fgets_lines(&data);
    // skip the comment block
    let mut it = lines.into_iter().skip_while(|l| l.first() == Some(&b'#'));
    let n = it
        .next()
        .and_then(scan_int)
        .map_or(0, |(v, _)| v as u32 as usize);
    let mut atoms = vec![AugmentedAtom::default(); n];
    for (i, atom) in atoms.iter_mut().enumerate() {
        let Some(line) = it.next() else { break };
        let parsed = quoted(line, 0).and_then(|(s, _)| string_to_augmented_atom(&s));
        match parsed {
            Some(aa) => *atom = aa,
            None => {
                crate::rd_log_error!(
                    "unsuccessful AA parsing of {} {}\n",
                    i + 1,
                    String::from_utf8_lossy(line)
                );
                return None;
            }
        }
    }
    Some(atoms)
}

/// `ReadAAPairs`: loader of the `*.trn` files; the pairs are the search and
/// the target pattern.
fn read_aa_pairs(path: &str) -> Result<Option<Vec<(AugmentedAtom, AugmentedAtom)>>> {
    let Some(data) = read_file(path) else {
        return Ok(None);
    };
    let Some((n, used)) = scan_int(&data) else {
        return Err(RDKitError::precondition("Failed Assertion: num_scan == 1"));
    };
    let n = n as u32 as usize;
    // the rest of the first line is the version string
    let rest = &data[used..];
    let first_line_end = rest
        .iter()
        .take(79)
        .position(|&b| b == b'\n')
        .map_or(rest.len().min(79), |p| p + 1);
    let lines = fgets_lines(&rest[first_line_end..]);
    let mut pairs = vec![(AugmentedAtom::default(), AugmentedAtom::default()); n];
    for (i, pair) in pairs.iter_mut().enumerate() {
        let Some(line) = lines.get(i) else { continue };
        let Some((first, next)) = quoted(line, 0) else {
            continue;
        };
        let fail = || {
            crate::rd_log_error!(
                "unsuccessful translation of {}\n",
                String::from_utf8_lossy(line)
            );
        };
        let Some(a) = string_to_augmented_atom(&first) else {
            fail();
            return Ok(None);
        };
        pair.0 = a;
        if let Some((second, _)) = quoted(line, next) {
            let Some(b) = string_to_augmented_atom(&second) else {
                fail();
                return Ok(None);
            };
            pair.1 = b;
        }
    }
    Ok(Some(pairs))
}

//=====================================================================
// Molecule data

fn load_sdf(path: &str) -> Option<Vec<Arc<ROMol>>> {
    let run = || -> Result<Vec<Arc<ROMol>>> {
        let mut suppl = SDMolSupplier::from_file(path, MolFileParserParams::default())?;
        let mut mols = Vec::new();
        while !suppl.at_end()? {
            if let Some(m) = suppl.next_mol()? {
                mols.push(Arc::new(m));
            }
        }
        Ok(mols)
    };
    match run() {
        Ok(m) => Some(m),
        Err(e) => {
            crate::rd_log_error!("{e}");
            None
        }
    }
}

fn parse_smarts(smarts: &[String]) -> Result<Vec<Arc<ROMol>>> {
    let mut mols = Vec::new();
    for patt in smarts {
        if let Some(m) = mol_from_smarts(patt, &SmartsParserParams::default())? {
            mols.push(Arc::new(m));
        }
    }
    Ok(mols)
}

/// An RD-file with reactions: the pairs of the first and second molecule
/// after each `$RXN`.
type TautomerRules = (Vec<Arc<ROMol>>, Vec<Arc<ROMol>>);

fn load_rdf(path: &str) -> Result<Option<TautomerRules>> {
    let Some(data) = read_file(path) else {
        return Ok(None);
    };
    let mut s = LineStream::new(&data);
    let params = MolFileParserParams::default();
    let (mut from, mut to) = (Vec::new(), Vec::new());
    let mut line = 0u32;
    let is = |l: &[u8], tag: &[u8]| l.len() >= 4 && l.starts_with(tag) && l.len() == tag.len();
    let corrupted = |line: u32| {
        crate::rd_log_error!("RD-file corrupted. {path} line {line}\n");
    };
    loop {
        let str_ = s.get_line();
        line += 1;
        if is(&str_, b"$RXN") {
            while !s.eof() {
                let str_ = s.get_line();
                line += 1;
                if is(&str_, b"$MOL") {
                    let Some(mol_from) = mol_from_mol_data_stream(&mut s, &mut line, &params)?
                    else {
                        corrupted(line);
                        return Ok(None);
                    };
                    let mol_from = Arc::new(mol_from);
                    while !s.eof() {
                        let str_ = s.get_line();
                        line += 1;
                        if is(&str_, b"$MOL") {
                            let Some(mol_to) =
                                mol_from_mol_data_stream(&mut s, &mut line, &params)?
                            else {
                                corrupted(line);
                                return Ok(None);
                            };
                            from.push(mol_from.clone());
                            to.push(Arc::new(mol_to));
                        }
                    }
                }
            }
        }
        if s.eof() {
            break;
        }
    }
    Ok(Some((from, to)))
}

impl StructCheckerOptions {
    pub fn new() -> Self {
        Self::default()
    }

    /// `clear`: back to the defaults.
    pub fn clear(&mut self) {
        *self = StructCheckerOptions::default();
    }

    pub fn load_augmented_atom_translations(&mut self, path: &str) -> Result<bool> {
        self.augmented_atom_pairs.clear();
        if path.is_empty() {
            return Ok(false);
        }
        match read_aa_pairs(path)? {
            Some(p) => {
                self.augmented_atom_pairs = p;
                Ok(true)
            }
            None => Ok(false),
        }
    }

    pub fn set_augmented_atom_translations(&mut self, pairs: Vec<(AugmentedAtom, AugmentedAtom)>) {
        self.augmented_atom_pairs = pairs;
    }

    pub fn load_acidic_augmented_atoms(&mut self, path: &str) -> bool {
        self.acidic_atoms.clear();
        if path.is_empty() {
            return false;
        }
        match read_augmented_atoms(path) {
            Some(a) => {
                self.acidic_atoms = a;
                true
            }
            None => false,
        }
    }

    pub fn set_acidic_augmented_atoms(&mut self, atoms: Vec<AugmentedAtom>) {
        self.acidic_atoms = atoms;
    }

    pub fn load_good_augmented_atoms(&mut self, path: &str) -> bool {
        self.good_atoms.clear();
        if path.is_empty() {
            return false;
        }
        match read_augmented_atoms(path) {
            Some(a) => {
                self.good_atoms = a;
                true
            }
            None => false,
        }
    }

    pub fn set_good_augmented_atoms(&mut self, atoms: Vec<AugmentedAtom>) {
        self.good_atoms = atoms;
    }

    /// Loads the clean patterns from an SD file.
    pub fn load_patterns(&mut self, path: &str) -> bool {
        self.patterns.clear();
        if path.is_empty() {
            return false;
        }
        match load_sdf(path) {
            Some(m) => {
                self.patterns = m;
                true
            }
            None => false,
        }
    }

    pub fn parse_patterns(&mut self, smarts: &[String]) -> Result<()> {
        self.patterns = parse_smarts(smarts)?;
        Ok(())
    }

    pub fn set_patterns(&mut self, p: Vec<Arc<ROMol>>) {
        self.patterns = p;
    }

    /// Loads the rotate patterns from an SD file.
    pub fn load_rotate_patterns(&mut self, path: &str) -> bool {
        self.rotate_patterns.clear();
        if path.is_empty() {
            return false;
        }
        match load_sdf(path) {
            Some(m) => {
                self.rotate_patterns = m;
                true
            }
            None => false,
        }
    }

    pub fn parse_rotate_patterns(&mut self, smarts: &[String]) -> Result<()> {
        self.rotate_patterns = parse_smarts(smarts)?;
        Ok(())
    }

    pub fn set_rotate_patterns(&mut self, p: Vec<Arc<ROMol>>) {
        self.rotate_patterns = p;
    }

    /// Loads the stereo patterns from an SD file.
    pub fn load_stereo_patterns(&mut self, path: &str) -> bool {
        self.stereo_patterns.clear();
        if path.is_empty() {
            return false;
        }
        match load_sdf(path) {
            Some(m) => {
                self.stereo_patterns = m;
                true
            }
            None => false,
        }
    }

    pub fn parse_stereo_patterns(&mut self, smarts: &[String]) -> Result<()> {
        self.stereo_patterns = parse_smarts(smarts)?;
        Ok(())
    }

    pub fn set_stereo_patterns(&mut self, p: Vec<Arc<ROMol>>) {
        self.stereo_patterns = p;
    }

    /// Loads the tautomer rules from an SD file (pairs of molecules) or an
    /// RD file (reactions).
    pub fn load_tautomer_data(&mut self, path: &str) -> bool {
        self.to_tautomer.clear();
        self.from_tautomer.clear();
        if path.is_empty() {
            return false;
        }
        let ext = path.rsplit_once('.').map_or(path, |(_, e)| e);
        let mut run = || -> Result<bool> {
            if ext == "sdf" || ext == "SDF" {
                let mut suppl = SDMolSupplier::from_file(path, MolFileParserParams::default())?;
                while !suppl.at_end()? {
                    let m1 = suppl.next_mol()?;
                    if suppl.at_end()? {
                        break;
                    }
                    let m2 = suppl.next_mol()?;
                    if let (Some(m1), Some(m2)) = (m1, m2) {
                        self.from_tautomer.push(Arc::new(m1));
                        self.to_tautomer.push(Arc::new(m2));
                    }
                }
                Ok(true)
            } else if ext == "rdf" || ext == "RDF" {
                match load_rdf(path)? {
                    Some((from, to)) => {
                        self.from_tautomer = from;
                        self.to_tautomer = to;
                        Ok(true)
                    }
                    None => Ok(false),
                }
            } else {
                crate::rd_log_error!("unsupported file type. {path}\n");
                Ok(false)
            }
        };
        match run() {
            Ok(b) => b,
            Err(e) => {
                crate::rd_log_error!("{e}\n");
                false
            }
        }
    }

    pub fn parse_tautomer_data(&mut self, from: &[String], to: &[String]) -> Result<()> {
        self.from_tautomer = parse_smarts(from)?;
        self.to_tautomer = parse_smarts(to)?;
        Ok(())
    }

    pub fn set_tautomer_data(&mut self, from: Vec<Arc<ROMol>>, to: Vec<Arc<ROMol>>) {
        self.from_tautomer = from;
        self.to_tautomer = to;
    }
}
