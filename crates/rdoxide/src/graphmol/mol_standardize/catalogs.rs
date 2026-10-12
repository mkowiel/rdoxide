//! The MolStandardize catalogs: acid/base pairs, fragment patterns,
//! normalization transforms and tautomer transforms.
//!
//! Ported from RDKit `Code/GraphMol/MolStandardize/{AcidBaseCatalog,
//! FragmentCatalog,TransformCatalog,TautomerCatalog}/*CatalogParams.cpp` and
//! `*CatalogUtils.cpp`. Upstream wraps each list in an `RDCatalog` that the
//! standardizers only ever use through its params object, so only the params
//! (the parsed entries) are ported.
//!
//! The default tables are compiled in upstream (`*.in` initializers); they
//! are extracted to `data/molstandardize/*.tsv` by `tools/sync_upstream.py`.

use std::sync::OnceLock;

use crate::error::{RDKitError, Result};
use crate::graphmol::bond::BondType;
use crate::graphmol::chem_reactions::ChemicalReaction;
use crate::graphmol::chem_reactions::daylight::{ReactionSmartsParserParams, reaction_from_smarts};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::rd_log_warning;
use crate::rdgeneral::common_properties;
use crate::rdgeneral::dict::RDProps;

const DEFAULT_ACID_BASE_PAIRS: &str =
    include_str!("../../../data/molstandardize/acid_base_pairs.tsv");
const DEFAULT_FRAGMENTS: &str = include_str!("../../../data/molstandardize/fragment_patterns.tsv");
const DEFAULT_NORMALIZATIONS: &str =
    include_str!("../../../data/molstandardize/normalizations.tsv");
const DEFAULT_TAUTOMER_TRANSFORMS: &str =
    include_str!("../../../data/molstandardize/tautomer_transforms.tsv");
const DEFAULT_TAUTOMER_TRANSFORMS_V1: &str =
    include_str!("../../../data/molstandardize/tautomer_transforms_v1.tsv");

/// The entries of an extracted default table (exactly `N` tab-separated
/// fields per line).
fn table<const N: usize>(text: &str) -> Vec<[String; N]> {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let fields: Vec<String> = l.split('\t').map(str::to_string).collect();
            fields
                .try_into()
                .expect("malformed MolStandardize data table")
        })
        .collect()
}

/// `defaults::defaultAcidBasePairs`: (name, acid SMARTS, base SMARTS).
pub fn default_acid_base_data() -> Vec<(String, String, String)> {
    table::<3>(DEFAULT_ACID_BASE_PAIRS)
        .into_iter()
        .map(|[a, b, c]| (a, b, c))
        .collect()
}

/// `defaults::defaultFragments`: (name, SMARTS).
pub fn default_fragment_data() -> Vec<(String, String)> {
    table::<2>(DEFAULT_FRAGMENTS)
        .into_iter()
        .map(|[a, b]| (a, b))
        .collect()
}

/// `defaults::defaultNormalizations`: (name, reaction SMARTS).
pub fn default_normalization_data() -> Vec<(String, String)> {
    table::<2>(DEFAULT_NORMALIZATIONS)
        .into_iter()
        .map(|[a, b]| (a, b))
        .collect()
}

/// (name, SMARTS, bonds, charges) — `TautomerTransformDefs`.
pub type TautomerTransformDefs = Vec<(String, String, String, String)>;

fn tautomer_defs(text: &str) -> TautomerTransformDefs {
    table::<4>(text)
        .into_iter()
        .map(|[a, b, c, d]| (a, b, c, d))
        .collect()
}

/// `defaults::defaultTautomerTransforms`
pub fn default_tautomer_transform_data() -> TautomerTransformDefs {
    tautomer_defs(DEFAULT_TAUTOMER_TRANSFORMS)
}

/// `defaults::defaultTautomerTransformsv1`
pub fn default_tautomer_transform_data_v1() -> TautomerTransformDefs {
    tautomer_defs(DEFAULT_TAUTOMER_TRANSFORMS_V1)
}

// ------------------------------------------------------------------ input

/// Reads a catalog file the way upstream's `std::ifstream` does.
fn read_file(file_name: &str) -> Result<String> {
    std::fs::read(file_name)
        .map(|b| String::from_utf8_lossy(&b).into_owned())
        .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))
}

/// The lines `while (!eof && !fail) getline(inLine, 512, '\n')` produces:
/// a line longer than 511 bytes is cut there and ends the reading.
fn stream_lines(text: &str) -> Vec<String> {
    const MAX_LINE_LEN: usize = 512;
    let mut out = Vec::new();
    let mut rest = text.as_bytes();
    loop {
        let end = rest.iter().position(|&c| c == b'\n');
        let line = &rest[..end.unwrap_or(rest.len())];
        if line.len() > MAX_LINE_LEN - 1 {
            out.push(String::from_utf8_lossy(&line[..MAX_LINE_LEN - 1]).into_owned());
            break;
        }
        out.push(String::from_utf8_lossy(line).into_owned());
        match end {
            Some(e) => rest = &rest[e + 1..],
            None => break,
        }
    }
    out
}

/// `boost::tokenizer` with `char_separator("\t")`: empty tokens are dropped.
fn tab_tokens(s: &str) -> Vec<String> {
    s.split('\t')
        .filter(|t| !t.is_empty())
        .map(str::to_string)
        .collect()
}

fn erase_spaces(s: &str) -> String {
    s.replace(' ', "")
}

/// `boost::trim` (the classic locale's whitespace).
fn boost_trim(s: &str) -> &str {
    s.trim_matches(|c: char| matches!(c, ' ' | '\t' | '\n' | '\r' | '\u{b}' | '\u{c}'))
}

fn is_blank_or_comment(s: &str) -> bool {
    s.is_empty() || s.starts_with("//")
}

fn smarts_to_mol(smarts: &str) -> Result<Option<ROMol>> {
    mol_from_smarts(smarts, &SmartsParserParams::default())
}

fn set_name<T: RDProps>(obj: &mut T, name: &str) {
    obj.dict_mut()
        .set_val(common_properties::_Name, name.to_string());
}

// --------------------------------------------------------------- acid/base

/// `AcidBaseCatalogParams`: acid/base pairs, strongest acid first.
#[derive(Debug, Clone, Default)]
pub struct AcidBaseCatalogParams {
    pairs: Vec<(ROMol, ROMol)>,
}

fn acid_base_pair(name: &str, acid_smarts: &str, base_smarts: &str) -> Result<(ROMol, ROMol)> {
    let mut acid = smarts_to_mol(acid_smarts)?
        .ok_or_else(|| RDKitError::Value(format!("Failed parsing acid SMARTS: {acid_smarts}")))?;
    let mut base = smarts_to_mol(base_smarts)?
        .ok_or_else(|| RDKitError::Value(format!("Failed parsing base SMARTS: {base_smarts}")))?;
    set_name(&mut acid, name);
    set_name(&mut base, name);
    Ok((acid, base))
}

impl AcidBaseCatalogParams {
    /// `AcidBaseCatalogParams(acidBaseFile)`: the defaults for an empty name.
    pub fn from_file(file_name: &str) -> Result<Self> {
        if file_name.is_empty() {
            return Self::from_data(&default_acid_base_data());
        }
        Self::from_text(&read_file(file_name)?)
    }

    /// `AcidBaseCatalogParams(std::istream &)`: one tab-separated
    /// `name acid base` triple per line; `//` comments.
    pub fn from_text(text: &str) -> Result<Self> {
        let mut pairs = Vec::new();
        for line in stream_lines(text) {
            let line = boost_trim(&line);
            if is_blank_or_comment(line) {
                continue;
            }
            let tokens = tab_tokens(line);
            if tokens.len() < 3 {
                return Err(RDKitError::Value(format!("bad acid/base line: {line}")));
            }
            pairs.push(acid_base_pair(
                &erase_spaces(&tokens[0]),
                &erase_spaces(&tokens[1]),
                &erase_spaces(&tokens[2]),
            )?);
        }
        Ok(AcidBaseCatalogParams { pairs })
    }

    /// `AcidBaseCatalogParams(data)`
    pub fn from_data(data: &[(String, String, String)]) -> Result<Self> {
        let pairs = data
            .iter()
            .map(|(n, a, b)| acid_base_pair(n, a, b))
            .collect::<Result<_>>()?;
        Ok(AcidBaseCatalogParams { pairs })
    }

    /// The (cached) default pairs.
    pub fn defaults() -> &'static AcidBaseCatalogParams {
        static DEFAULTS: OnceLock<AcidBaseCatalogParams> = OnceLock::new();
        DEFAULTS.get_or_init(|| {
            Self::from_data(&default_acid_base_data()).expect("default acid/base pairs")
        })
    }

    /// `getPairs()`: (acid, base) query molecules.
    pub fn pairs(&self) -> &[(ROMol, ROMol)] {
        &self.pairs
    }
}

// ---------------------------------------------------------------- fragments

/// `FragmentCatalogParams`: the fragment patterns `FragmentRemover` uses.
#[derive(Debug, Clone, Default)]
pub struct FragmentCatalogParams {
    func_groups: Vec<ROMol>,
}

fn fragment_mol(name: &str, smarts: &str) -> Result<ROMol> {
    let mut mol = smarts_to_mol(smarts)?
        .ok_or_else(|| RDKitError::Value(format!("Failed parsing fragment SMARTS: {smarts}")))?;
    set_name(&mut mol, name);
    mol.dict_mut()
        .set_val(common_properties::_fragSMARTS, smarts.to_string());
    Ok(mol)
}

impl FragmentCatalogParams {
    /// `FragmentCatalogParams(fgroupFile)`: the defaults for an empty name.
    pub fn from_file(file_name: &str) -> Result<Self> {
        if file_name.is_empty() {
            return Self::from_data(&default_fragment_data());
        }
        Self::from_text(&read_file(file_name)?)
    }

    /// `FragmentCatalogParams(std::istream &)`: tab-separated `name SMARTS`
    /// lines; `//` comments.
    pub fn from_text(text: &str) -> Result<Self> {
        let mut func_groups = Vec::new();
        for line in stream_lines(text) {
            let line = boost_trim(&line);
            if is_blank_or_comment(line) {
                continue;
            }
            let tokens = tab_tokens(line);
            if tokens.len() < 2 {
                return Err(RDKitError::Value(format!(
                    "no SMARTS found in input: {line}"
                )));
            }
            func_groups.push(fragment_mol(&tokens[0], &erase_spaces(&tokens[1]))?);
        }
        Ok(FragmentCatalogParams { func_groups })
    }

    /// `FragmentCatalogParams(data)`
    pub fn from_data(data: &[(String, String)]) -> Result<Self> {
        let func_groups = data
            .iter()
            .map(|(n, s)| fragment_mol(n, s))
            .collect::<Result<_>>()?;
        Ok(FragmentCatalogParams { func_groups })
    }

    /// The (cached) default fragments.
    pub fn defaults() -> &'static FragmentCatalogParams {
        static DEFAULTS: OnceLock<FragmentCatalogParams> = OnceLock::new();
        DEFAULTS.get_or_init(|| {
            Self::from_data(&default_fragment_data()).expect("default fragment patterns")
        })
    }

    /// `getFuncGroups()`
    pub fn func_groups(&self) -> &[ROMol] {
        &self.func_groups
    }
}

// ----------------------------------------------------------- normalizations

/// `TransformCatalogParams`: the normalization reactions. The reactant
/// matchers are initialized on construction (upstream does it lazily, in
/// `initializeTransforms`, before the first use).
#[derive(Debug, Clone, Default)]
pub struct TransformCatalogParams {
    transformations: Vec<ChemicalReaction>,
}

fn transformation(name: &str, smarts: &str) -> Result<ChemicalReaction> {
    let mut rxn = reaction_from_smarts(smarts, &ReactionSmartsParserParams::new())?;
    set_name(&mut rxn, name);
    rxn.init_reactant_matchers(false);
    Ok(rxn)
}

impl TransformCatalogParams {
    /// `TransformCatalogParams(transformFile)`: the defaults for an empty name.
    pub fn from_file(file_name: &str) -> Result<Self> {
        if file_name.is_empty() {
            return Self::from_data(&default_normalization_data());
        }
        Self::from_text(&read_file(file_name)?)
    }

    /// `TransformCatalogParams(std::istream &)`: tab-separated
    /// `name reaction-SMARTS` lines; `//` comments.
    pub fn from_text(text: &str) -> Result<Self> {
        let mut transformations = Vec::new();
        for line in stream_lines(text) {
            if is_blank_or_comment(&line) {
                continue;
            }
            let tokens = tab_tokens(&line);
            if tokens.len() < 2 {
                return Err(RDKitError::Value(format!(
                    "bad transformation line: {line}"
                )));
            }
            transformations.push(transformation(
                &erase_spaces(&tokens[0]),
                &erase_spaces(&tokens[1]),
            )?);
        }
        Ok(TransformCatalogParams { transformations })
    }

    /// `TransformCatalogParams(data)`
    pub fn from_data(data: &[(String, String)]) -> Result<Self> {
        let transformations = data
            .iter()
            .map(|(n, s)| transformation(n, s))
            .collect::<Result<_>>()?;
        Ok(TransformCatalogParams { transformations })
    }

    /// The (cached) default normalizations.
    pub fn defaults() -> &'static TransformCatalogParams {
        static DEFAULTS: OnceLock<TransformCatalogParams> = OnceLock::new();
        DEFAULTS.get_or_init(|| {
            Self::from_data(&default_normalization_data()).expect("default normalizations")
        })
    }

    /// `getTransformations()`
    pub fn transformations(&self) -> &[ChemicalReaction] {
        &self.transformations
    }
}

// ----------------------------------------------------------------- tautomers

/// `TautomerTransform`: a SMARTS whose matched atoms are re-bonded with
/// `bond_types` (default: single/double alternation flipped) and recharged
/// with `charges`.
#[derive(Debug, Clone)]
pub struct TautomerTransform {
    pub mol: ROMol,
    pub bond_types: Vec<BondType>,
    pub charges: Vec<i32>,
}

/// `stringToBondType`: `-`, `=`, `#`, `:`; other characters are ignored.
pub fn string_to_bond_type(bond_str: &str) -> Vec<BondType> {
    bond_str
        .chars()
        .filter_map(|c| match c {
            '-' => Some(BondType::Single),
            '=' => Some(BondType::Double),
            '#' => Some(BondType::Triple),
            ':' => Some(BondType::Aromatic),
            _ => None,
        })
        .collect()
}

/// `stringToCharge`: `+`, `0`, `-`.
pub fn string_to_charge(charge_str: &str) -> Result<Vec<i32>> {
    charge_str
        .chars()
        .map(|c| match c {
            '+' => Ok(1),
            '0' => Ok(0),
            '-' => Ok(-1),
            _ => Err(RDKitError::Value("Charge symbol not recognised.".into())),
        })
        .collect()
}

fn tautomer_transform(
    name: &str,
    smarts: &str,
    bond_str: &str,
    charge_str: &str,
) -> Result<TautomerTransform> {
    let bond_types = string_to_bond_type(bond_str);
    let charges = string_to_charge(charge_str)?;
    let mut mol = smarts_to_mol(smarts)?
        .ok_or_else(|| RDKitError::Value(format!("cannot parse tautomer SMARTS: {smarts}")))?;
    set_name(&mut mol, name);
    Ok(TautomerTransform {
        mol,
        bond_types,
        charges,
    })
}

/// `TautomerCatalogParams`: the tautomer transforms.
#[derive(Debug, Clone, Default)]
pub struct TautomerCatalogParams {
    transforms: Vec<TautomerTransform>,
}

impl TautomerCatalogParams {
    /// `TautomerCatalogParams(tautomerFile)`: the defaults for an empty name.
    pub fn from_file(file_name: &str) -> Result<Self> {
        if file_name.is_empty() {
            return Self::from_data(&default_tautomer_transform_data());
        }
        Self::from_text(&read_file(file_name)?)
    }

    /// `readTautomers(std::istream &)`: tab-separated
    /// `name SMARTS [bonds [charges]]` lines; `//` comments.
    pub fn from_text(text: &str) -> Result<Self> {
        let mut transforms = Vec::new();
        for line in stream_lines(text) {
            if is_blank_or_comment(&line) {
                continue;
            }
            let tokens = tab_tokens(&line);
            if tokens.len() < 2 {
                rd_log_warning!("Invalid line: {line}");
                continue;
            }
            // upstream leaves every field empty for more than four tokens
            let field = |i: usize| {
                if tokens.len() <= 4 {
                    tokens.get(i).map_or(String::new(), |t| erase_spaces(t))
                } else {
                    String::new()
                }
            };
            transforms.push(tautomer_transform(
                &field(0),
                &field(1),
                &field(2),
                &field(3),
            )?);
        }
        Ok(TautomerCatalogParams { transforms })
    }

    /// `TautomerCatalogParams(data)`
    pub fn from_data(data: &[(String, String, String, String)]) -> Result<Self> {
        let transforms = data
            .iter()
            .map(|(n, s, b, c)| tautomer_transform(n, s, b, c))
            .collect::<Result<_>>()?;
        Ok(TautomerCatalogParams { transforms })
    }

    /// The (cached) default transforms.
    pub fn defaults() -> &'static TautomerCatalogParams {
        static DEFAULTS: OnceLock<TautomerCatalogParams> = OnceLock::new();
        DEFAULTS.get_or_init(|| {
            Self::from_data(&default_tautomer_transform_data())
                .expect("default tautomer transforms")
        })
    }

    /// The (cached) transforms of `getV1TautomerEnumerator`.
    pub fn defaults_v1() -> &'static TautomerCatalogParams {
        static DEFAULTS: OnceLock<TautomerCatalogParams> = OnceLock::new();
        DEFAULTS.get_or_init(|| {
            Self::from_data(&default_tautomer_transform_data_v1())
                .expect("default v1 tautomer transforms")
        })
    }

    /// `getTransforms()`
    pub fn transforms(&self) -> &[TautomerTransform] {
        &self.transforms
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stream_lines_like_getline() {
        assert_eq!(stream_lines("a\nb\n"), vec!["a", "b", ""]);
        assert_eq!(stream_lines("a"), vec!["a"]);
        let long = format!("{}\nnext", "x".repeat(600));
        let lines = stream_lines(&long);
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0].len(), 511);
    }

    #[test]
    fn defaults_parse() {
        assert_eq!(AcidBaseCatalogParams::defaults().pairs().len(), 33);
        assert_eq!(FragmentCatalogParams::defaults().func_groups().len(), 61);
        assert_eq!(
            TransformCatalogParams::defaults().transformations().len(),
            22
        );
        assert_eq!(TautomerCatalogParams::defaults().transforms().len(), 37);
        assert_eq!(TautomerCatalogParams::defaults_v1().transforms().len(), 36);
    }
}
