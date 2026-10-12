//! The standardization `Pipeline`: parse a mol block, validate it,
//! standardize it and make its parent, logging what happened.
//!
//! Ported from RDKit `Code/GraphMol/MolStandardize/Pipeline.{h,cpp}`.

use std::collections::BTreeMap;

use super::charge::Uncharger;
use super::fragment::LargestFragmentChooser;
use super::metal::{MetalDisconnector, MetalDisconnectorOptions};
use super::normalize::Normalizer;
use super::validate::{
    FeaturesValidation, Layout2DValidation, ValidationMethod, squared_median_bond_length,
};
use crate::error::{RDKitError, Result};
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::file_parsers::mol_file_parser::{MolFileParserParams, mol_from_mol_block};
use crate::graphmol::file_parsers::mol_file_writer::{
    MolWriterParams, mol_to_v2k_mol_block, mol_to_v3k_mol_block,
};
use crate::graphmol::mol_ops::sanitize_flags::{
    SANITIZE_ALL, SANITIZE_CLEANUP, SANITIZE_CLEANUP_ORGANOMETALLICS, SANITIZE_FINDRADICALS,
    SANITIZE_SYMMRINGS,
};
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::RWMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::graphmol::wedge_bonds::reapply_mol_block_wedging;

/// The default normalizations of the pipeline.
const DEFAULT_NORMALIZER_DATA: &str = concat!(
    "// Name\tSMIRKS\n",
    "Nitro to N+(O-)=O\t[N,P,As,Sb;X3:1](=[O,S,Se,Te:2])=[O,S,Se,Te:3]>>[*+1:1]([*-1:2])=[*:3]\n",
    "Sulfone to S(=O)(=O)\t[S+2:1]([O-:2])([O-:3])>>[S+0:1](=[O-0:2])(=[O-0:3])\n",
    "Pyridine oxide to n+O-\t[nH0+0:1]=[OH0+0:2]>>[n+:1][O-:2]\n",
    "Azide to N=N+=N-\t[*:1][N:2]=[N:3]#[N:4]>>[*:1][N:2]=[N+:3]=[N-:4]\n",
    "Diazo/azo to =N+=N-\t[*:1]=[N:2]#[N:3]>>[*:1]=[N+:2]=[N-:3]\n",
    "[SH](=O)(=O) to S(=O)O\t[c,C,N,O,F,Cl,Br,I:1][SH+0:2](=[O:3])=[O:4]>>[*:1][*:2]([*:3])=[*:4]\n",
    "Phosphate to P(O-)=O\t[O-:1][P+;D4:2][O,S,Se,Te;-1:3]>>[O+0:1]=[P+0;D5:2][*-1:3]\n",
    "Generalized phosphate to P(X-)=Y\t[S,Se,Te;-1:1][P+;D4:2][S,Se,Te;-1:3]>>[*+0:1]=[P+0;D5:2][*-1:3]\n",
    "C/S+N to C/S=N+\t[C,S&!$([S+]-[O-]);X3+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]\n",
    "P+N to P=N+\t[P;X4+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]\n",
    "Recombine 1,3-separated charges\t[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[N,P,As,Sb,O,S,Se,Te;+1:3]>>[*-0:1]=[*:2]-[*+0:3]\n",
    "Recombine 1,3-separated charges\t[n,o,p,s;-1:1]:[a:2]=[N,O,P,S;+1:3]>>[*-0:1]:[*:2]-[*+0:3]\n",
    "Recombine 1,3-separated charges\t[N,O,P,S;-1:1]-[a+0:2]:[n,o,p,s;+1:3]>>[*-0:1]=[*:2]:[*+0:3]\n",
    "Recombine 1,5-separated charges\t[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[A:3]-[A:4]=[N,P,As,Sb,O,S,Se,Te;+1:5]>>[*-0:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\n",
    "Recombine 1,5-separated charges\t[n,o,p,s;-1:1]:[a:2]:[a:3]:[c:4]=[N,O,P,S;+1:5]>>[*-0:1]:[*:2]:[*:3]:[c:4]-[*+0:5]\n",
    "Recombine 1,5-separated charges\t[N,O,P,S;-1:1]-[c:2]:[a:3]:[a:4]:[n,o,p,s;+1:5]>>[*-0:1]=[c:2]:[*:3]:[*:4]:[*+0:5]\n",
    "Normalize 1,3 conjugated cation\t[N;+0!H0:1]@-[A:2]=[N!$(*~[N,O,P,S;-1]),O;+1H0:3]>>[*+1:1]=[*:2]-[*+0:3]\n",
    "Normalize 1,5 conjugated cation\t[N;+0!H0:1]@-[A:2]=[A:3]@-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\n",
    "Normalize 1,3 conjugated cation\t[N,O!$(*N);+0!H0:1]-[A:2]=[N!$(*~[N,O,P,S;-1]),O;+1H0:3]>>[*+1:1]=[*:2]-[*+0:3]\n",
    "Normalize 1,3 conjugated cation\t[n;+0!H0:1]:[c:2]=[N!$(*~[N,O,P,S;-1]),O;+1H0:3]>>[*+1:1]:[*:2]-[*+0:3]\n",
    "Normalize 1,5 conjugated cation\t[N;+0!H0:1]@-[A:2]=[A:3]-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\n",
    "Normalize 1,5 conjugated cation\t[N,O!$(*N);+0!H0:1]-[A:2]=[A:3]@-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\n",
    "Normalize 1,5 conjugated cation\t[N,O!$(*N);+0!H0:1]-[A:2]=[A:3]-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\n",
    "Normalize 1,5 conjugated cation\t[n;+0!H0:1]:[a:2]:[a:3]:[c:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[n+1:1]:[*:2]:[*:3]:[*:4]-[*+0:5]\n",
    "Charge normalization\t[F,Cl,Br,I,At;-1:1]=[O:2]>>[*-0:1][O-:2]\n",
    "Charge recombination\t[N,P,As,Sb;-1:1]=[C+;v3:2]>>[*+0:1]#[C+0:2]\n",
);

/// `PipelineOptions`
#[derive(Debug, Clone, PartialEq)]
pub struct PipelineOptions {
    pub strict_parsing: bool,
    pub report_all_failures: bool,
    pub allow_empty_molecules: bool,
    pub allow_enhanced_stereo: bool,
    pub allow_aromatic_bond_type: bool,
    pub allow_dative_bond_type: bool,
    pub is_2d_zero_threshold: f64,
    pub atom_clash_limit: f64,
    pub min_median_bond_length: f64,
    pub bond_length_limit: f64,
    pub allow_long_bonds_in_rings: bool,
    pub allow_atom_bond_clash_exemption: bool,
    pub metal_nof: String,
    pub metal_non: String,
    pub normalizer_data: String,
    pub normalizer_max_restarts: u32,
    pub scaled_median_bond_length: f64,
    pub output_v2000: bool,
}

impl Default for PipelineOptions {
    fn default() -> Self {
        PipelineOptions {
            strict_parsing: false,
            report_all_failures: true,
            allow_empty_molecules: false,
            allow_enhanced_stereo: false,
            allow_aromatic_bond_type: false,
            allow_dative_bond_type: false,
            is_2d_zero_threshold: 1e-3,
            atom_clash_limit: 0.03,
            min_median_bond_length: 1e-3,
            bond_length_limit: 100.0,
            allow_long_bonds_in_rings: true,
            allow_atom_bond_clash_exemption: true,
            metal_nof: "[Li,Na,K,Rb,Cs,Fr]~[#7,#8,F]".into(),
            metal_non: String::new(),
            normalizer_data: DEFAULT_NORMALIZER_DATA.into(),
            normalizer_max_restarts: 200,
            scaled_median_bond_length: 1.0,
            output_v2000: false,
        }
    }
}

/// `PipelineStatus` flags.
pub mod status {
    pub const NO_EVENT: u32 = 0;
    pub const INPUT_ERROR: u32 = 1 << 0;
    pub const PREPARE_FOR_VALIDATION_ERROR: u32 = 1 << 1;
    pub const FEATURES_VALIDATION_ERROR: u32 = 1 << 2;
    pub const BASIC_VALIDATION_ERROR: u32 = 1 << 3;
    pub const IS2D_VALIDATION_ERROR: u32 = 1 << 4;
    pub const LAYOUT2D_VALIDATION_ERROR: u32 = 1 << 5;
    pub const STEREO_VALIDATION_ERROR: u32 = 1 << 6;
    pub const VALIDATION_ERROR: u32 = FEATURES_VALIDATION_ERROR
        | BASIC_VALIDATION_ERROR
        | IS2D_VALIDATION_ERROR
        | LAYOUT2D_VALIDATION_ERROR
        | STEREO_VALIDATION_ERROR;
    pub const PREPARE_FOR_STANDARDIZATION_ERROR: u32 = 1 << 7;
    pub const METAL_STANDARDIZATION_ERROR: u32 = 1 << 8;
    pub const NORMALIZER_STANDARDIZATION_ERROR: u32 = 1 << 9;
    pub const FRAGMENT_STANDARDIZATION_ERROR: u32 = 1 << 10;
    pub const CHARGE_STANDARDIZATION_ERROR: u32 = 1 << 11;
    pub const STANDARDIZATION_ERROR: u32 = METAL_STANDARDIZATION_ERROR
        | NORMALIZER_STANDARDIZATION_ERROR
        | FRAGMENT_STANDARDIZATION_ERROR
        | CHARGE_STANDARDIZATION_ERROR;
    pub const OUTPUT_ERROR: u32 = 1 << 12;
    pub const PIPELINE_ERROR: u32 = INPUT_ERROR
        | PREPARE_FOR_VALIDATION_ERROR
        | VALIDATION_ERROR
        | PREPARE_FOR_STANDARDIZATION_ERROR
        | STANDARDIZATION_ERROR
        | OUTPUT_ERROR;
    pub const METALS_DISCONNECTED: u32 = 1 << 23;
    pub const NORMALIZATION_APPLIED: u32 = 1 << 24;
    pub const FRAGMENTS_REMOVED: u32 = 1 << 25;
    pub const PROTONATION_CHANGED: u32 = 1 << 26;
    pub const STRUCTURE_MODIFICATION: u32 =
        METALS_DISCONNECTED | NORMALIZATION_APPLIED | FRAGMENTS_REMOVED | PROTONATION_CHANGED;
}
use status::*;

/// `PipelineStage`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PipelineStage {
    NotStarted = 0,
    ParsingInput,
    PrepareForValidation,
    Validation,
    PrepareForStandardization,
    Standardization,
    ReapplyWedging,
    Cleanup2D,
    MakeParent,
    SerializingOutput,
    Completed,
}

/// `PipelineLogEntry`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PipelineLogEntry {
    pub status: u32,
    pub detail: String,
}

/// `PipelineResult`
#[derive(Debug, Clone, Default)]
pub struct PipelineResult {
    pub status: u32,
    pub stage: u32,
    pub log: Vec<PipelineLogEntry>,
    pub input_mol_data: String,
    pub output_mol_data: String,
    pub parent_mol_data: String,
}

impl PipelineResult {
    /// `append`
    pub fn append(&mut self, new_status: u32, info: &str) {
        self.status |= new_status;
        self.log.push(PipelineLogEntry {
            status: new_status,
            detail: info.to_string(),
        });
    }

    fn failed(&self, options: &PipelineOptions) -> bool {
        (self.status & PIPELINE_ERROR) != NO_EVENT && !options.report_all_failures
    }
}

fn smiles(mol: &mut RWMol) -> Result<String> {
    mol_to_smiles(mol, &SmilesWriteParams::default())
}

fn is_sanitize_error(e: &RDKitError) -> bool {
    matches!(e, RDKitError::Sanitize(_))
}

/// A pipeline step: returns `None` when no molecule is left.
pub type Operation = fn(RWMol, &mut PipelineResult, &PipelineOptions) -> Result<Option<RWMol>>;

/// `Pipeline`
#[derive(Debug, Clone, Default)]
pub struct Pipeline {
    pub options: PipelineOptions,
}

impl Pipeline {
    pub fn new(options: PipelineOptions) -> Self {
        Pipeline { options }
    }

    /// `run`: standardizes `molblock`.
    pub fn run(&self, molblock: &str) -> Result<PipelineResult> {
        let options = &self.options;
        let mut result = PipelineResult {
            input_mol_data: molblock.to_string(),
            ..Default::default()
        };
        result.stage = PipelineStage::ParsingInput as u32;
        let Some(mol) = parse(molblock, &mut result, options)? else {
            return Ok(result);
        };
        if result.failed(options) {
            return Ok(result);
        }
        // sanitization and validation run on a copy: the input is kept
        let validation: [(PipelineStage, Operation); 2] = [
            (PipelineStage::PrepareForValidation, prepare_for_validation),
            (PipelineStage::Validation, validate),
        ];
        let mut mol_copy = mol.clone();
        for (stage, op) in validation {
            result.stage = stage as u32;
            match op(mol_copy, &mut result, options)? {
                Some(m) if !result.failed(options) => mol_copy = m,
                _ => return Ok(result),
            }
        }
        let standardization: [(PipelineStage, Operation); 4] = [
            (
                PipelineStage::PrepareForStandardization,
                prepare_for_standardization,
            ),
            (PipelineStage::Standardization, standardize),
            (PipelineStage::ReapplyWedging, reapply_wedging),
            (PipelineStage::Cleanup2D, cleanup_2d),
        ];
        let mut mol = mol;
        for (stage, op) in standardization {
            result.stage = stage as u32;
            match op(mol, &mut result, options)? {
                Some(m) if !result.failed(options) => mol = m,
                _ => return Ok(result),
            }
        }
        result.stage = PipelineStage::MakeParent as u32;
        let Some((output, parent)) = make_parent(mol, &mut result, options)? else {
            return Ok(result);
        };
        if result.failed(options) {
            return Ok(result);
        }
        result.stage = PipelineStage::SerializingOutput as u32;
        serialize(&output, &parent, &mut result, options);
        if result.failed(options) {
            return Ok(result);
        }
        result.stage = PipelineStage::Completed as u32;
        Ok(result)
    }
}

/// `Operations::parse`
pub fn parse(
    molblock: &str,
    result: &mut PipelineResult,
    options: &PipelineOptions,
) -> Result<Option<RWMol>> {
    let params = MolFileParserParams {
        // no sanitization (and no H removal) at this stage
        sanitize: false,
        remove_hs: false,
        strict_parsing: options.strict_parsing,
        ..Default::default()
    };
    let mol = match mol_from_mol_block(molblock, &params) {
        Ok(m) => m,
        Err(RDKitError::FileParse(msg)) => {
            result.append(INPUT_ERROR, &msg);
            None
        }
        Err(e) => return Err(e),
    };
    if mol.is_none() {
        result.append(
            INPUT_ERROR,
            "Could not instantiate a valid molecule from input",
        );
    }
    Ok(mol)
}

/// `Operations::serialize`
pub fn serialize(
    output: &RWMol,
    parent: &RWMol,
    result: &mut PipelineResult,
    options: &PipelineOptions,
) {
    let params = MolWriterParams::default();
    // the output block is kept when writing the parent fails
    let write = |result: &mut PipelineResult| -> Result<()> {
        if options.output_v2000 {
            result.output_mol_data = mol_to_v2k_mol_block(output, &params, -1)?;
            result.parent_mol_data = mol_to_v2k_mol_block(parent, &params, -1)?;
        } else {
            result.output_mol_data = mol_to_v3k_mol_block(output, &params, -1)?;
            result.parent_mol_data = mol_to_v3k_mol_block(parent, &params, -1)?;
        }
        Ok(())
    };
    match write(result) {
        Ok(()) => {}
        Err(RDKitError::Value(e)) if options.output_v2000 => result.append(
            OUTPUT_ERROR,
            &format!("Can't write molecule to V2000 output format: {e}"),
        ),
        Err(e) => result.append(
            OUTPUT_ERROR,
            &format!("Can't write molecule to output format: {e}"),
        ),
    }
}

/// `Operations::prepareForValidation`: minimal sanitization, stereo
/// assignment and the original wedging.
pub fn prepare_for_validation(
    mut mol: RWMol,
    result: &mut PipelineResult,
    _options: &PipelineOptions,
) -> Result<Option<RWMol>> {
    const OPS: u32 = SANITIZE_CLEANUP
        | SANITIZE_SYMMRINGS
        | SANITIZE_CLEANUP_ORGANOMETALLICS
        | SANITIZE_FINDRADICALS;
    let res = sanitize_mol(&mut mol, OPS)
        .map_err(|(_, e)| e)
        .and_then(|()| assign_stereochemistry(&mut mol, true, true, true))
        .and_then(|()| reapply_mol_block_wedging(&mut mol, true, false));
    match res {
        Ok(()) => {}
        Err(e) if is_sanitize_error(&e) => result.append(
            PREPARE_FOR_VALIDATION_ERROR,
            "An error occurred while preparing the molecule for validation.",
        ),
        Err(e) => return Err(e),
    }
    Ok(Some(mol))
}

/// The regex `^(ERROR|INFO): \[.+\] `: the validation message prefix (`.+`
/// is greedy, so the prefix ends at the last `] ` of the first line).
pub fn remove_error_prefix(message: &str) -> String {
    for head in ["ERROR: [", "INFO: ["] {
        let Some(rest) = message.strip_prefix(head) else {
            continue;
        };
        let line = rest.split('\n').next().unwrap_or("");
        if let Some(pos) = line.rfind("] ") {
            if pos >= 1 {
                return rest[pos + 2..].to_string();
            }
        }
    }
    message.to_string()
}

/// `Operations::validate`
pub fn validate(
    mol: RWMol,
    result: &mut PipelineResult,
    options: &PipelineOptions,
) -> Result<Option<RWMol>> {
    let validations: [(ValidationMethod, u32); 7] = [
        // undesired features (e.g., query atoms/bonds)
        (
            ValidationMethod::Features(FeaturesValidation {
                allow_enhanced_stereo: options.allow_enhanced_stereo,
                allow_aromatic_bond_type: options.allow_aromatic_bond_type,
                allow_dative_bond_type: options.allow_dative_bond_type,
                ..Default::default()
            }),
            FEATURES_VALIDATION_ERROR,
        ),
        // number of atoms and valences
        (
            ValidationMethod::RDKit {
                allow_empty_molecules: options.allow_empty_molecules,
            },
            BASIC_VALIDATION_ERROR,
        ),
        (ValidationMethod::DisallowedRadical, BASIC_VALIDATION_ERROR),
        (
            ValidationMethod::Isotope { strict: true },
            BASIC_VALIDATION_ERROR,
        ),
        (
            ValidationMethod::Is2D {
                threshold: options.is_2d_zero_threshold,
            },
            IS2D_VALIDATION_ERROR,
        ),
        (
            ValidationMethod::Layout2D(Layout2DValidation {
                clash_limit: options.atom_clash_limit,
                bond_length_limit: options.bond_length_limit,
                allow_long_bonds_in_rings: options.allow_long_bonds_in_rings,
                allow_atom_bond_clash_exemption: options.allow_atom_bond_clash_exemption,
                min_median_bond_length: options.min_median_bond_length,
            }),
            LAYOUT2D_VALIDATION_ERROR,
        ),
        (ValidationMethod::Stereo, STEREO_VALIDATION_ERROR),
    ];
    for (v, st) in &validations {
        let errors = v.validate(&mol, options.report_all_failures)?;
        for e in &errors {
            result.append(*st, &remove_error_prefix(e));
        }
        if !errors.is_empty() && !options.report_all_failures {
            break;
        }
    }
    Ok(Some(mol))
}

/// `Operations::prepareForStandardization`
pub fn prepare_for_standardization(
    mut mol: RWMol,
    result: &mut PipelineResult,
    _options: &PipelineOptions,
) -> Result<Option<RWMol>> {
    match sanitize_mol(&mut mol, SANITIZE_ALL) {
        Ok(()) => {}
        Err((_, e)) if is_sanitize_error(&e) => result.append(
            PREPARE_FOR_STANDARDIZATION_ERROR,
            "An error occurred while preparing the molecule for standardization.",
        ),
        Err((_, e)) => return Err(e),
    }
    Ok(Some(mol))
}

fn smarts(text: &str) -> Result<RWMol> {
    mol_from_smarts(text, &SmartsParserParams::default())?
        .ok_or_else(|| RDKitError::Value(format!("Failed parsing SMARTS: {text}")))
}

/// `Operations::standardize`: metals, functional groups, largest fragment.
pub fn standardize(
    mut mol: RWMol,
    result: &mut PipelineResult,
    options: &PipelineOptions,
) -> Result<Option<RWMol>> {
    let mut reference = smiles(&mut mol)?;

    // bonding to metals (any error is caught upstream)
    let metals = (|| -> Result<()> {
        let mut md = MetalDisconnector::new(MetalDisconnectorOptions::default());
        md.set_metal_nof(smarts(&options.metal_nof)?);
        md.set_metal_non(smarts(&options.metal_non)?);
        md.disconnect_in_place(&mut mol)
    })();
    if metals.is_err() {
        result.append(
            METAL_STANDARDIZATION_ERROR,
            "An error occurred while processing the bonding of metal species.",
        );
        return Ok(Some(mol));
    }
    let s = smiles(&mut mol)?;
    if s != reference {
        result.append(
            METALS_DISCONNECTED,
            "One or more metal atoms were disconnected.",
        );
    }
    reference = s;

    // functional groups; normalize() and not normalizeInPlace(), which may
    // leave an ill-formed molecule when a sanitization fails (GitHub #7189)
    let normalized = (|| -> Result<RWMol> {
        let normalizer = if options.normalizer_data.is_empty() {
            Normalizer::default()
        } else {
            Normalizer::from_text(&options.normalizer_data, options.normalizer_max_restarts)?
        };
        let mut m = normalizer.normalize(&mol)?;
        m.update_property_cache(false)?;
        Ok(m)
    })();
    match normalized {
        Ok(m) => mol = m,
        Err(_) => {
            result.append(
                NORMALIZER_STANDARDIZATION_ERROR,
                "An error occurred while normalizing the representation of some functional groups",
            );
            return Ok(Some(mol));
        }
    }
    let s = smiles(&mut mol)?;
    if s != reference {
        result.append(
            NORMALIZATION_APPLIED,
            "The representation of some functional groups was adjusted.",
        );
    }
    reference = s;

    // keep the largest fragment
    if LargestFragmentChooser::default()
        .choose_in_place(&mut mol)
        .is_err()
    {
        result.append(
            FRAGMENT_STANDARDIZATION_ERROR,
            "An error occurred while removing the disconnected fragments",
        );
        return Ok(Some(mol));
    }
    let s = smiles(&mut mol)?;
    if s != reference {
        result.append(
            FRAGMENTS_REMOVED,
            "One or more disconnected fragments (e.g., counterions) were removed.",
        );
    }
    // the stereochemistry is assigned once the graph is final
    assign_stereochemistry(&mut mol, true, true, true)?;
    Ok(Some(mol))
}

/// `Operations::reapplyWedging`: the input wedging, without wavy bonds.
pub fn reapply_wedging(
    mut mol: RWMol,
    result: &mut PipelineResult,
    _options: &PipelineOptions,
) -> Result<Option<RWMol>> {
    let old: BTreeMap<u32, BondStereo> = mol
        .bonds()
        .iter()
        .enumerate()
        .map(|(i, b)| (i as u32, b.stereo()))
        .collect();
    // 1) the original wedging of the input mol block
    reapply_mol_block_wedging(&mut mol, true, false)?;
    // 2) double bonds of "either" type: restore STEREOANY on double bonds
    //    with an adjacent wavy bond that are now STEREONONE
    for i in 0..mol.num_bonds() {
        let bond = mol.bond(i);
        if bond.bond_type() != BondType::Double {
            continue;
        }
        let old_stereo = old.get(&i).copied().unwrap_or(BondStereo::StereoNone);
        let new_stereo = bond.stereo();
        let mut has_adjacent_wavy = false;
        for atom in [bond.begin_atom_idx(), bond.end_atom_idx()] {
            for adj in mol.atom_bonds(atom) {
                if adj.idx() != i && adj.bond_dir() == BondDir::Unknown {
                    has_adjacent_wavy = true;
                }
            }
        }
        if has_adjacent_wavy
            && old_stereo == BondStereo::StereoAny
            && new_stereo == BondStereo::StereoNone
        {
            mol.bond_mut(i).set_stereo(BondStereo::StereoAny)?;
            result.append(
                NORMALIZATION_APPLIED,
                &format!(
                    "Double bond {i} was assigned an undefined/unknown stereochemical configuration"
                ),
            );
        }
    }
    // 3) no "unknown" bond directions
    for i in 0..mol.num_bonds() {
        if mol.bond(i).bond_dir() != BondDir::Unknown {
            continue;
        }
        mol.bond_mut(i).set_bond_dir(BondDir::None);
        result.append(
            NORMALIZATION_APPLIED,
            &format!("The \"wavy\" style of bond {i} was removed"),
        );
    }
    Ok(Some(mol))
}

/// `Operations::cleanup2D`: scales the coordinates to the configured median
/// bond length and zeroes the z coordinates.
pub fn cleanup_2d(
    mut mol: RWMol,
    _result: &mut PipelineResult,
    options: &PipelineOptions,
) -> Result<Option<RWMol>> {
    if options.scaled_median_bond_length > 0.0 && mol.num_conformers() > 0 {
        let median = squared_median_bond_length(&mol, mol.conformer(-1)?)?.sqrt();
        if median > options.min_median_bond_length {
            let scale = options.scaled_median_bond_length / median;
            let conf = mol.conformer_mut(-1)?;
            for p in conf.positions_mut().iter_mut() {
                p.x *= scale;
                p.y *= scale;
                p.z = 0.0;
            }
        }
    }
    Ok(Some(mol))
}

fn replace_dative_bonds(mol: &mut RWMol) -> Result<()> {
    let mut modified = false;
    for i in 0..mol.num_bonds() {
        let b = mol.bond(i);
        if b.bond_type() != BondType::Dative {
            continue;
        }
        let (donor, acceptor) = (b.begin_atom_idx(), b.end_atom_idx());
        let d = mol.atom_mut(donor);
        d.set_formal_charge(d.formal_charge() + 1);
        let a = mol.atom_mut(acceptor);
        a.set_formal_charge(a.formal_charge() - 1);
        mol.bond_mut(i).set_bond_type(BondType::Single);
        modified = true;
    }
    if modified {
        mol.update_property_cache(false)?;
    }
    Ok(())
}

fn remove_hs_at_protonated_sites(mol: &mut RWMol) -> Result<()> {
    let mut protons = vec![false; mol.num_atoms() as usize];
    for i in 0..mol.num_atoms() {
        if mol.atom(i).atomic_num() != 1 || mol.atom_degree(i) != 1 {
            continue;
        }
        if mol
            .atom_neighbors(i)
            .any(|n| mol.atom(n).formal_charge() > 0)
        {
            protons[i as usize] = true;
        }
    }
    if protons.iter().any(|&p| p) {
        for idx in (0..mol.num_atoms()).rev() {
            if !protons[idx as usize] {
                continue;
            }
            let nbr = mol.atom_neighbors(idx).next();
            if let Some(nbr) = nbr {
                let n = mol.atom_mut(nbr);
                n.set_num_explicit_hs(n.num_explicit_hs() + 1);
            }
            mol.remove_atom(idx)?;
        }
        mol.update_property_cache(false)?;
    }
    Ok(())
}

/// `Operations::makeParent`: (output, parent); the output is replaced by
/// the parent when its protonation state is unsuitable.
pub fn make_parent(
    mut mol: RWMol,
    result: &mut PipelineResult,
    _options: &PipelineOptions,
) -> Result<Option<(RWMol, RWMol)>> {
    let reference = smiles(&mut mol)?;
    let mut parent = mol.clone();
    // the overall charge status (Hs next to positive sites are removed first:
    // the Uncharger would not find those sites)
    let uncharged = remove_hs_at_protonated_sites(&mut parent)
        .and_then(|()| Uncharger::new(false, true, true).uncharge_in_place(&mut parent));
    if uncharged.is_err() {
        result.append(
            CHARGE_STANDARDIZATION_ERROR,
            "An error occurred while normalizing the compound's charge status",
        );
        return Ok(None);
    }
    let parent_charge: i32 = parent.atoms().iter().map(|a| a.formal_charge()).sum();
    let mol_charge: i32 = mol.atoms().iter().map(|a| a.formal_charge()).sum();
    // mol is accepted if neutral or partially/fully balancing the parent's
    // non-neutralizable charges
    // upstream then shares one molecule between the output and the parent
    let use_parent = (mol_charge > 0 && mol_charge > parent_charge)
        || (mol_charge < 0 && mol_charge < parent_charge);
    if use_parent {
        mol = parent.clone();
    }
    if smiles(&mut mol)? != reference {
        result.append(PROTONATION_CHANGED, "The protonation state was adjusted.");
    }
    replace_dative_bonds(&mut parent)?;
    if use_parent {
        mol = parent.clone();
    }
    Ok(Some((mol, parent)))
}
