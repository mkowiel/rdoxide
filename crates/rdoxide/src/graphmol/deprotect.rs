//! Port of RDKit `Code/GraphMol/Deprotect`: removal of protecting groups,
//! described as reactions that leave behind the group being protected.

use std::sync::{Arc, OnceLock};

use crate::graphmol::chem_reactions::ChemicalReaction;
use crate::graphmol::chem_reactions::daylight::{ReactionSmartsParserParams, reaction_from_smarts};
use crate::graphmol::chem_reactions::runner::{run_reactant, run_reactant_in_place};
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::log::LogStateSetter;
use crate::{RDKitError, Result, rd_log_error};

/// A single deprotection: the reaction that removes the protecting group and
/// information about it.
#[derive(Debug, Clone)]
pub struct DeprotectData {
    /// functional group being protected (amine, alcohol, ...)
    pub deprotection_class: String,
    /// the reaction SMARTS that removes the protecting group
    pub reaction_smarts: String,
    /// common abbreviation of the protecting group (Boc, Fmoc, ...)
    pub abbreviation: String,
    pub full_name: String,
    pub example: String,
    rxn: Arc<ChemicalReaction>,
}

impl DeprotectData {
    pub fn new(
        deprotection_class: &str,
        reaction_smarts: &str,
        abbreviation: &str,
        full_name: &str,
        example: &str,
    ) -> Result<Self> {
        let mut rxn = reaction_from_smarts(reaction_smarts, &ReactionSmartsParserParams::new())?;
        if rxn.num_product_templates() != 1 {
            rd_log_error!("Deprotection reactions must have exactly one product\n");
        }
        rxn.init_reactant_matchers(false);
        Ok(DeprotectData {
            deprotection_class: deprotection_class.into(),
            reaction_smarts: reaction_smarts.into(),
            abbreviation: abbreviation.into(),
            full_name: full_name.into(),
            example: example.into(),
            rxn: Arc::new(rxn),
        })
    }

    pub fn is_valid(&self) -> bool {
        self.rxn.num_product_templates() == 1
    }

    pub fn reaction(&self) -> &ChemicalReaction {
        &self.rxn
    }
}

impl PartialEq for DeprotectData {
    fn eq(&self, other: &Self) -> bool {
        self.deprotection_class == other.deprotection_class
            && self.full_name == other.full_name
            && self.abbreviation == other.abbreviation
            && self.reaction_smarts == other.reaction_smarts
            && self.is_valid() == other.is_valid()
    }
}

/// (class, reaction SMARTS, abbreviation, full name, example)
const DEPROTECTIONS: &[(&str, &str, &str, &str, &str)] = &[
    (
        "alcohol",
        "CC(C)([Si](C)(C)[O;H0:1])C>>[O;H1:1]",
        "TBDMS",
        "tert-butyldimethylsilyl",
        "CC(C)(C)[Si](C)(C)OC>>OC",
    ),
    (
        "alcohol",
        "COc1ccc(C([O;H0:1])(c2ccccc2)c2ccccc2)cc1>>[O;H1:1]",
        "MMT",
        "methoxytrityl",
        "COc1ccc(C(ON)(c2ccccc2)c2ccccc2)cc1>>NO",
    ),
    (
        "alcohol",
        "O1C([O;H0:1])[C;H2][C;H2][C;H2][C;H2]1>>[O;H1:1]",
        "THP",
        "tetrahydropyranyl",
        "COC1CCCCO1>>CO",
    ),
    (
        "alcohol",
        "[C;H3][O;X2&R0][C;H2&R0][O;H0:1]>>[O;H1:1]",
        "MOM",
        "methoxymethyl_ether",
        "COCOC>>CO",
    ),
    (
        "alcohol",
        "[C;R0][O;R0][C;R0][C;R0][O;R0][C;R0][O;X2:1]>>[O;H1:1]",
        "MEM",
        "beta-Methoxyethoxymethyl_ether",
        "C(C)(C)OCOCCOC>>C(C)(C)O",
    ),
    (
        "alcohol",
        "[O;!$(*C(=O)):1][CH2]c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[O;H1:1]",
        "Bn",
        "benzyl",
        "NOCc1ccccc1>>ON",
    ),
    (
        "alcohol",
        "[O;H0&X2:1][Si]([C;H3])([C;H3])([C;H3])>>[O;H1:1]",
        "TMS",
        "trimethylsilyl",
        "NO[Si](C)(C)C>>NO",
    ),
    (
        "alcohol",
        "[O;H0:1]C(=O)c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[O;H1:1]",
        "Bz",
        "benzoyl",
        "COC(=O)c1ccccc1>>OC",
    ),
    (
        "alcohol",
        "[O;H0:1][C;R0](=O)[C;R0]([C;H3])([C;H3])[C;H3]>>[O;H1:1]",
        "Piv",
        "pivaloyl",
        "CC(C)(C)C(=O)ON>>NO",
    ),
    (
        "alcohol",
        "[O;H0:1][Si](C(C)C)(C(C)C)(C(C)C)>>[O;H1:1]",
        "TIPS",
        "triisopropylsilyl",
        "NO[Si](C(C)C)(C(C)C)(C(C)C)>>ON",
    ),
    (
        "alcohol",
        "[O;R0:1][C;R0](=O)[C;H3]>>[O:1]",
        "Ac",
        "acetyl",
        "NOC(=O)C>>NO",
    ),
    (
        "alcohol",
        "[c;H1]1[c;H1]c(O[C;H3])[c;H1][c;H1]c1[C;H2][O;D2&R0:1]>>[O;H1:1]",
        "PMB",
        "para-methoxybenzyl_ether",
        "c1cc(OC)ccc1CON>>NO",
    ),
    (
        "alcohol",
        "c1ccc(C([NX3;H0,H1,H2:1])(c2ccccc2)c2ccccc2)cc1>>[N:1]",
        "Tr",
        "triphenylmethyl",
        "CN(C)C(c1ccccc1)(c1ccccc1)c1ccccc1>>CNC",
    ),
    (
        "amine",
        "O=C1[N;H0:1]C(c2[c;H1][c;H1][c;H1][c;H1]c21)=O>>[N:1]",
        "Phth",
        "phthaloyl",
        "O=C1N(C(=O)c2ccccc12)C1=CC=CC=C1>>Nc1ccccc1",
    ),
    (
        "amine",
        "[#7:1]COCC[Si]([C;R0])([C;R0])[C;R0]>>[#7:1]",
        "SEM",
        "trimethylsilylethoxymethyl",
        "C[Si](C)(C)CCOCNO>>NO",
    ),
    (
        "amine",
        "[C;H3;R0][O;R0]c1[c;H1][c;H1][c]([NX3;H0,H1;!$(NC=O):1])[c;H1][c;H1]1>>[N:1]",
        "PMP",
        "para-methoxyphenyl",
        "COc1ccc(NO)cc1>>NO",
    ),
    (
        "amine",
        "[C;H3]c1[c;H1][c;H1]c(S(=O)(=O)[NX3;H0,H1;!$(NC=O):1])[c;H1][c;H1]1>>[N:1]",
        "Ts",
        "tosyl",
        "Cc1ccc(cc1)S(NO)(=O)=O>>NO",
    ),
    (
        "amine",
        "[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]",
        "Boc",
        "tert-butyloxycarbonyl",
        "N(C(=O)OC(C)(C)C)Cc1ccccc1NC(=O)OC(C)(C)C>>NCc1ccccc1N",
    ),
    (
        "amine",
        "[N;H0,H1:1]C(=O)C(F)(F)F>>[N:1]",
        "TFA",
        "trifluoroacetyl",
        "ONC(=O)C(F)(F)F>>NO",
    ),
    (
        "amine",
        "[NX3;H0,H1:1]C(=O)c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[N:1]",
        "Bz",
        "benzoyl",
        "ONC(=O)c1ccccc1>>NO",
    ),
    (
        "amine",
        "[NX3;H0,H1:1][#6](=O)-[#8]-[#6]-[#6]-1-c2ccccc2-c2ccccc-12>>[N:1]",
        "Fmoc",
        "9-fluorenylmethyloxycarbonyl",
        "CN(C)C(=O)OCC1C2=CC=CC=C2C2=C1C=CC=C2>>CNC",
    ),
    (
        "amine",
        "[NX3;H0,H1;!$(NC=O):1][C;H2]c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[N:1]",
        "Bn",
        "benzylamine",
        "CN(C)Cc1ccccc1>>CNC",
    ),
    (
        "amine",
        "[NX3;H0,H1:1][C;R0](=O)[C;H3]>>[N:1]",
        "Ac",
        "acetamide",
        "N(C)(C)C(=O)C>>CNC",
    ),
    (
        "amine",
        "[NX3;H0,H1:1][C;R0](=O)[O;R0][C;R0]c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[N:1]",
        "Cbz",
        "carbobenzyloxy",
        "N(C)(C)C(=O)OCc1ccccc1>>CNC",
    ),
    (
        "carbonyl",
        "[#6,#1:1][C:2]([#6,#1:3])([#8:4][#6])([#8][#6])>>[*:1][C:2]([*:3])[#8:4]",
        "Acetyl/Ketal",
        "Acetal/Ketal",
        "O=C(OC)C1CCC(OC2)(OC2)CC1>>COC(=O)C1CCC(O)CC1",
    ),
];

/// `getDeprotections`: the built-in deprotections
pub fn get_deprotections() -> &'static [DeprotectData] {
    static DATA: OnceLock<Vec<DeprotectData>> = OnceLock::new();
    DATA.get_or_init(|| {
        DEPROTECTIONS
            .iter()
            .map(|&(class, smarts, abbreviation, name, example)| {
                DeprotectData::new(class, smarts, abbreviation, name, example)
                    .expect("built-in deprotection reaction")
            })
            .collect()
    })
}

const MAX_DEPROTECTIONS: usize = 100;

fn annotate(mol: &mut ROMol, used: Vec<String>) -> Result<()> {
    let count = used.len() as i32;
    mol.set_prop("DEPROTECTIONS", used, false)?;
    mol.set_prop("DEPROTECTION_COUNT", count, false)
}

/// Whether the molecule can be sanitized (the logs are off meanwhile).
fn sanitizes(mol: &mut RWMol) -> bool {
    let _blocker = LogStateSetter::new();
    !matches!(
        sanitize_mol(mol, SANITIZE_ALL),
        Err((_, RDKitError::Sanitize(_)))
    )
}

/// `deprotect`: the deprotected form of the molecule, annotated with the
/// deprotections used (`DEPROTECTIONS`) and their number
/// (`DEPROTECTION_COUNT`).
pub fn deprotect(mol: &ROMol, deprotections: &[DeprotectData]) -> Result<ROMol> {
    let mut used: Vec<String> = Vec::new();
    let mut m = mol.clone();
    let mut something_happened = true;
    while something_happened {
        something_happened = false;
        for d in deprotections.iter().filter(|d| d.is_valid()) {
            for mut prods in run_reactant(d.reaction(), &m, 0, 1)? {
                // like upstream, the product replaces the molecule even if it
                // then turns out not to sanitize
                m = prods.swap_remove(0);
                if !sanitizes(&mut m) {
                    continue;
                }
                used.push(d.abbreviation.clone());
                if used.len() >= MAX_DEPROTECTIONS {
                    rd_log_error!("Too many deprotections, halting...\n");
                } else {
                    something_happened = true;
                }
                break;
            }
        }
    }
    annotate(&mut m, used)?;
    Ok(m)
}

/// `deprotectInPlace`: whether the molecule was changed.
pub fn deprotect_in_place(mol: &mut RWMol, deprotections: &[DeprotectData]) -> Result<bool> {
    let mut used: Vec<String> = Vec::new();
    let mut modified = false;
    let mut something_happened = true;
    while something_happened {
        something_happened = false;
        for d in deprotections.iter().filter(|d| d.is_valid()) {
            if run_reactant_in_place(d.reaction(), mol, true)? {
                if !sanitizes(mol) {
                    continue;
                }
                modified = true;
                used.push(d.abbreviation.clone());
                if used.len() >= MAX_DEPROTECTIONS {
                    rd_log_error!("Too many deprotections, halting...\n");
                } else {
                    something_happened = true;
                }
            }
        }
    }
    annotate(mol, used)?;
    Ok(modified)
}
