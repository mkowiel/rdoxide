//! MMFF atom typing and parameter assignment for a molecule (upstream
//! `MMFFMolProperties` of `ForceFieldHelpers/MMFF/AtomTyper.cpp`).

use super::heavy::HeavyAtomTyper;
use super::rings::{
    RingMembershipSize, is_angle_in_ring_of_size_3_or_4, is_torsion_in_ring_of_size_4_or_5,
};
use crate::error::{RDKitError, Result};
use crate::forcefield::mmff::params::{MmffParams, Prop};
use crate::graphmol::aromaticity::set_mmff_aromaticity;
use crate::graphmol::bond::BondType;
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::mol_ops::add_hs::needs_hs;
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
        if props.valid {
            props.compute_charges(mol)?;
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

    pub(super) fn params(&self) -> &'static MmffParams {
        MmffParams::get(self.mmffs)
    }

    pub(super) fn prop(&self, idx: u32) -> &'static Prop {
        self.params()
            .prop(u32::from(self.atom_type(idx)))
            .expect("MMFF atom types always have properties")
    }

    pub(super) fn require_valid(&self) -> Result<()> {
        if self.valid {
            Ok(())
        } else {
            Err(RDKitError::precondition(
                "missing atom types - invalid force-field",
            ))
        }
    }

    /// The MMFF bond type: 1 for a single bond between two atoms that both
    /// have the `sbmb` flag or are both aromatic, else 0.
    pub fn bond_type(&self, mol: &ROMol, bond: u32) -> Result<u32> {
        self.require_valid()?;
        let bond = mol.bond(bond);
        let (p1, p2) = (
            self.prop(bond.begin_atom_idx()),
            self.prop(bond.end_atom_idx()),
        );
        Ok(u32::from(
            bond.bond_type() == BondType::Single
                && ((p1.sbmb != 0 && p2.sbmb != 0) || (p1.arom != 0 && p2.arom != 0)),
        ))
    }

    fn bond_type_between(&self, mol: &ROMol, a: u32, b: u32) -> Result<u32> {
        match mol.bond_idx_between(a, b) {
            Some(bond) => self.bond_type(mol, bond),
            None => Err(RDKitError::precondition("invalid bond")),
        }
    }

    /// The MMFF angle type of the angle idx1-idx2-idx3:
    /// 0 normal, 1/2 one/two bonds of type 1, 3/4 in a 3-/4-membered ring,
    /// 5/6 and 7/8 in a 3-/4-membered ring with a bond type sum of 1/2.
    pub fn angle_type(&self, mol: &ROMol, idx1: u32, idx2: u32, idx3: u32) -> Result<u32> {
        self.require_valid()?;
        let bond_type_sum =
            self.bond_type_between(mol, idx1, idx2)? + self.bond_type_between(mol, idx2, idx3)?;
        let mut angle_type = bond_type_sum;
        let size = is_angle_in_ring_of_size_3_or_4(mol, idx1, idx2, idx3);
        if size != 0 {
            angle_type = size;
            if bond_type_sum != 0 {
                angle_type += bond_type_sum + size - 2;
            }
        }
        Ok(angle_type)
    }

    /// The principal and the secondary torsion type of idx1-idx2-idx3-idx4;
    /// the secondary one is used when no parameters exist for the principal.
    pub fn torsion_type(
        &self,
        mol: &ROMol,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
    ) -> Result<(u32, u32)> {
        self.require_valid()?;
        let bond_jk = mol
            .bond_idx_between(idx2, idx3)
            .ok_or_else(|| RDKitError::precondition("invalid bond"))?;
        let bond_type_ij = self.bond_type_between(mol, idx1, idx2)?;
        let bond_type_jk = self.bond_type(mol, bond_jk)?;
        let bond_type_kl = self.bond_type_between(mol, idx3, idx4)?;
        let mut torsion_type = bond_type_jk;
        let mut second_torsion_type = 0;
        // MMFF.IV asks for less, but CYGUAN01 only passes with this condition
        if bond_type_jk == 0
            && mol.bond(bond_jk).bond_type() == BondType::Single
            && (bond_type_ij == 1 || bond_type_kl == 1)
        {
            torsion_type = 2;
        }
        let size = is_torsion_in_ring_of_size_4_or_5(mol, idx1, idx2, idx3, idx4);
        // the check for bonds I-K and J-L keeps two fused 3-membered rings
        // from getting type 4
        if size == 4
            && !(mol.bond_idx_between(idx1, idx3).is_some()
                || mol.bond_idx_between(idx2, idx4).is_some())
        {
            second_torsion_type = torsion_type;
            torsion_type = 4;
        } else if size == 5
            && [idx1, idx2, idx3, idx4]
                .iter()
                .any(|&i| self.atom_type(i) == 1)
        {
            second_torsion_type = torsion_type;
            torsion_type = 5;
        }
        Ok((torsion_type, second_torsion_type))
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

/// The MMFF stretch-bend type from the angle type and the types of the two
/// bonds of the angle.
pub fn stretch_bend_type(angle_type: u32, bond_type1: u32, bond_type2: u32) -> u32 {
    let first_or_equal = bond_type1 != 0 || bond_type1 == bond_type2;
    match angle_type {
        1 => {
            if first_or_equal {
                1
            } else {
                2
            }
        }
        2 => 3,
        4 => 4,
        3 => 5,
        5 => {
            if first_or_equal {
                6
            } else {
                7
            }
        }
        6 => 8,
        7 => {
            if first_or_equal {
                9
            } else {
                10
            }
        }
        8 => 11,
        _ => 0,
    }
}
