//! MMFF atom typing and parameter assignment for a molecule (upstream
//! `MMFFMolProperties` of `ForceFieldHelpers/MMFF/AtomTyper.cpp`).

use super::heavy::HeavyAtomTyper;
use super::rings::RingMembershipSize;
use crate::error::Result;
use crate::graphmol::aromaticity::set_mmff_aromaticity;
use crate::graphmol::mol_ops::add_hs::needs_hs;
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::mol_ops::sanitize_flags::*;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rd_log_warning;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

/// How the electrostatic term depends on the distance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DielectricModel {
    Constant = 1,
    Distance = 2,
}

/// The MMFF properties of one atom.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct MmffAtomProperties {
    pub atom_type: u8,
    pub formal_charge: f64,
    pub partial_charge: f64,
}

/// `sanitizeMMFFMol`: sanitizes without perceiving aromaticity. Returns the
/// operation that failed (`SANITIZE_NONE` on success).
pub fn sanitize_mmff_mol(mol: &mut RWMol) -> u32 {
    let ops = SANITIZE_CLEANUP
        | SANITIZE_PROPERTIES
        | SANITIZE_SYMMRINGS
        | SANITIZE_KEKULIZE
        | SANITIZE_FINDRADICALS
        | SANITIZE_SETCONJUGATION
        | SANITIZE_SETHYBRIDIZATION
        | SANITIZE_CLEANUPCHIRALITY
        | SANITIZE_ADJUSTHS;
    match sanitize_mol(mol, ops) {
        Ok(()) => {
            if !mol.has_prop(cp::_MMFFSanitized) {
                let _ = mol.set_prop(cp::_MMFFSanitized, 1i32, true);
            }
            SANITIZE_NONE
        }
        Err((flag, _)) => flag,
    }
}

/// The atom types, charges and settings of an MMFF calculation.
#[derive(Debug, Clone)]
pub struct MmffMolProperties {
    pub(super) valid: bool,
    pub(super) mmffs: bool,
    pub(super) bond_term: bool,
    pub(super) angle_term: bool,
    pub(super) stretch_bend_term: bool,
    pub(super) oop_term: bool,
    pub(super) torsion_term: bool,
    pub(super) vdw_term: bool,
    pub(super) ele_term: bool,
    pub(super) diel_const: f64,
    pub(super) diel_model: DielectricModel,
    pub(super) atoms: Vec<MmffAtomProperties>,
}

impl MmffMolProperties {
    /// Assigns the atom types and charges; the molecule is kekulized if it
    /// was not MMFF-sanitized and gets MMFF aromaticity. When an atom type
    /// is missing the properties are invalid and the charges stay 0.
    pub fn new(mol: &mut ROMol, mmff_variant: &str) -> Result<Self> {
        let mut props = MmffMolProperties {
            valid: true,
            mmffs: mmff_variant == "MMFF94s",
            bond_term: true,
            angle_term: true,
            stretch_bend_term: true,
            oop_term: true,
            torsion_term: true,
            vdw_term: true,
            ele_term: true,
            diel_const: 1.0,
            diel_model: DielectricModel::Constant,
            atoms: vec![MmffAtomProperties::default(); mol.num_atoms() as usize],
        };
        if needs_hs(mol)? {
            rd_log_warning!("Molecule does not have explicit Hs. Consider calling AddHs()");
        }
        if !mol.has_prop(cp::_MMFFSanitized) {
            if (0..mol.num_atoms()).any(|i| mol.atom(i).is_aromatic()) {
                kekulize(mol, true, true, 100)?;
            }
            mol.set_prop(cp::_MMFFSanitized, 1i32, true)?;
        }
        set_mmff_aromaticity(mol)?;
        let rm = RingMembershipSize::new(mol)?;
        let heavy = HeavyAtomTyper::new(mol, &rm)?;
        for i in 0..mol.num_atoms() {
            if mol.atom(i).atomic_num() != 1 {
                let atom_type = heavy.atom_type(i)?;
                props.set_atom_type(i, atom_type);
            }
        }
        for i in 0..mol.num_atoms() {
            if mol.atom(i).atomic_num() == 1 {
                let atom_type = props.hydrogen_type(mol, i);
                props.set_atom_type(i, atom_type);
            }
        }
        Ok(props)
    }

    fn set_atom_type(&mut self, idx: u32, atom_type: u32) {
        self.atoms[idx as usize].atom_type = atom_type as u8;
        if atom_type == 0 {
            self.valid = false;
        }
    }

    pub fn is_valid(&self) -> bool {
        self.valid
    }

    pub fn atom_type(&self, idx: u32) -> u8 {
        self.atoms[idx as usize].atom_type
    }

    pub fn formal_charge(&self, idx: u32) -> f64 {
        self.atoms[idx as usize].formal_charge
    }

    pub fn partial_charge(&self, idx: u32) -> f64 {
        self.atoms[idx as usize].partial_charge
    }

    /// The type of a hydrogen atom, which follows from its (typed) neighbor.
    fn hydrogen_type(&self, mol: &ROMol, idx: u32) -> u32 {
        let mut atom_type = 0;
        for nbr in mol.atom_neighbors(idx) {
            atom_type = match mol.atom(nbr).atomic_num() {
                // HC, HSI
                6 | 14 => 5,
                7 => match self.atom_type(nbr) {
                    // HNR, H3N, HPYL, HNOX
                    8 | 39 | 62 | 67 | 68 => 23,
                    // HNN+, HGD+, HIM+ and other cationic nitrogens
                    34 | 54 | 55 | 56 | 58 | 81 => 36,
                    // HN=N, HN=C
                    9 => 27,
                    // HNCC, HNCN, HNCO, HNCS, HNNC, HNNN, HNSO, HNC%, HSP2
                    _ => 28,
                },
                8 => self.hydroxyl_hydrogen_type(mol, nbr),
                // HP, HS, HS=N
                15 | 16 => 71,
                _ => 0,
            };
        }
        atom_type
    }

    /// Hydrogens on oxygen: acidic hydrogens of carboxylic and phosphoric
    /// acids, enolic or phenolic hydrogens and hydrogens on oxygens bonded to
    /// sulfur are told apart from alcohols and generic hydroxyl hydrogens.
    fn hydroxyl_hydrogen_type(&self, mol: &ROMol, oxygen: u32) -> u32 {
        use crate::graphmol::bond::BondType;
        match self.atom_type(oxygen) {
            49 => return 50, // HO+
            51 => return 52, // HO=+
            70 => return 31, // HOH
            6 => {}
            _ => return 21, // HO
        }
        let (mut is_hocc_or_hocn, mut is_hoco, mut is_hop, mut is_hos) =
            (false, false, false, false);
        for nbr2 in mol.atom_neighbors(oxygen) {
            match mol.atom(nbr2).atomic_num() {
                6 => {
                    for nbr3 in mol.atom_neighbors(nbr2) {
                        if nbr3 == oxygen {
                            continue;
                        }
                        let bt = mol
                            .bond_idx_between(nbr2, nbr3)
                            .map(|b| mol.bond(b).bond_type());
                        let z3 = mol.atom(nbr3).atomic_num();
                        if matches!(z3, 6 | 7)
                            && matches!(bt, Some(BondType::Double | BondType::Aromatic))
                        {
                            is_hocc_or_hocn = true;
                        }
                        if z3 == 8 && bt == Some(BondType::Double) {
                            is_hoco = true;
                        }
                    }
                }
                15 => is_hop = true,
                16 => is_hos = true,
                _ => {}
            }
        }
        if is_hoco || is_hop {
            24 // HOCO
        } else if is_hocc_or_hocn {
            29 // HOCC, HOCN
        } else if is_hos {
            33 // HOS
        } else {
            21 // HO, HOR
        }
    }
}
