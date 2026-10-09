//! Stereoisomer enumeration.
//!
//! Ported from RDKit `Code/GraphMol/EnumerateStereoisomers/{EnumerateStereoisomers,Flippers}.{h,cpp}`
//! (2026.09.1).
//!
//! Copyright (C) 2025 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! The molecule's potential stereo elements are found with
//! [`find_potential_stereo`]; every one of them gets a [`Flipper`] that sets
//! one of its two configurations. Isomers are drawn at random (without
//! repetition) from the 2^N possible flipper settings, with the same engine
//! (`std::mt19937` and a libstdc++ `bernoulli_distribution`) as upstream so
//! seeded runs give the same sequence.

use std::collections::HashSet;

use super::atom::ChiralType;
use super::bond::{BondDir, BondStereo};
use super::chirality::bond_dirs::set_double_bond_neighbor_directions;
use super::chirality::find_stereo::{
    NOATOM, StereoInfo, StereoSpecified, StereoType, find_potential_stereo,
};
use super::chirality::legacy::assign_stereochemistry;
use super::romol::ROMol;
use super::smiles_parse::cx_write::{RestoreBondDirOption, cx_fields, mol_to_cx_smiles};
use super::smiles_parse::write::SmilesWriteParams;
use super::stereo_group::StereoGroupType;
use super::wedge_bonds::{BondWedgingParameters, wedge_mol_bonds};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::random::Mt19937;
use crate::{RDKitError, Result};

/// `StereoEnumerationOptions`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoEnumerationOptions {
    /// try to embed every isomer; not supported yet (needs distance geometry)
    pub try_embedding: bool,
    /// leave specified stereo elements alone unless they are in a relative
    /// stereo group
    pub only_unassigned: bool,
    /// only vary the stereo groups of the molecule
    pub only_stereo_groups: bool,
    /// only return isomers with different canonical CXSMILES
    pub unique: bool,
    /// the maximum number of isomers to return, 0 for no limit
    pub max_isomers: u64,
    /// seed of the random number generator, -1 for a non-deterministic one
    pub random_seed: i32,
}

impl Default for StereoEnumerationOptions {
    fn default() -> Self {
        StereoEnumerationOptions {
            try_embedding: false,
            only_unassigned: true,
            only_stereo_groups: false,
            unique: true,
            max_isomers: 0,
            random_seed: -1,
        }
    }
}

/// Sets one of the two configurations of a stereo element (`details::Flipper`).
#[derive(Debug, Clone)]
enum Flipper {
    Atom(u32),
    Bond(u32),
    /// the atoms of a stereo group with their original parities
    StereoGroup(Vec<(u32, ChiralType)>),
    Atropisomer {
        bond: u32,
        ctrl: [u32; 4],
    },
}

impl Flipper {
    fn flip(&self, mol: &mut ROMol, flag: bool) -> Result<()> {
        match self {
            Flipper::Atom(atom) => mol.atom_mut(*atom).set_chiral_tag(if flag {
                ChiralType::TetrahedralCW
            } else {
                ChiralType::TetrahedralCCW
            }),
            Flipper::Bond(bond) => mol.bond_mut(*bond).set_stereo(if flag {
                BondStereo::StereoCis
            } else {
                BondStereo::StereoTrans
            })?,
            Flipper::StereoGroup(parities) => {
                for &(atom, original) in parities {
                    let tag = match (flag, original) {
                        (true, _) => Some(original),
                        (false, ChiralType::TetrahedralCW) => Some(ChiralType::TetrahedralCCW),
                        (false, ChiralType::TetrahedralCCW) => Some(ChiralType::TetrahedralCW),
                        (false, _) => None,
                    };
                    if let Some(tag) = tag {
                        mol.atom_mut(atom).set_chiral_tag(tag);
                    }
                }
            }
            Flipper::Atropisomer { bond, ctrl } => flip_atropisomer(mol, *bond, ctrl, flag)?,
        }
        Ok(())
    }
}

fn flip_atropisomer(mol: &mut ROMol, bond: u32, ctrl: &[u32; 4], flag: bool) -> Result<()> {
    let (begin, end) = (
        mol.bond(bond).begin_atom_idx(),
        mol.bond(bond).end_atom_idx(),
    );
    let find = |atom: u32, other: u32| -> Result<u32> {
        mol.get_bond_between_atoms(atom, other)?
            .map(|b| b.idx())
            .ok_or_else(|| RDKitError::Value("atropisomer controlling atom is not bonded".into()))
    };
    let (b1, b2) = (find(begin, ctrl[0])?, find(begin, ctrl[1])?);
    let (b3, b4) = (find(end, ctrl[2])?, find(end, ctrl[3])?);
    mol.bond_mut(b3).set_bond_dir(BondDir::None);
    mol.bond_mut(b4).set_bond_dir(BondDir::None);
    let (stereo, dir) = if flag {
        (BondStereo::StereoAtropCCW, BondDir::BeginWedge)
    } else {
        (BondStereo::StereoAtropCW, BondDir::BeginDash)
    };
    mol.bond_mut(bond).set_stereo(stereo)?;
    let (marked, cleared) = if ctrl[0] < ctrl[1] {
        (b1, b2)
    } else {
        (b2, b1)
    };
    mol.bond_mut(marked).set_bond_dir(dir);
    mol.bond_mut(cleared).set_bond_dir(BondDir::None);
    Ok(())
}

/// Builds the flipper of a double bond; `None` if it cannot be flipped
/// (no usable stereo atoms).
fn bond_flipper(mol: &mut ROMol, si: &StereoInfo) -> Option<Flipper> {
    let bond = si.centered_on;
    if mol.bond(bond).stereo_atoms().is_empty() {
        let (a, b) = (si.controlling_atoms[0], si.controlling_atoms[2]);
        if a == NOATOM || b == NOATOM {
            return None;
        }
        *mol.bond_mut(bond).stereo_atoms_mut() = vec![a as i32, b as i32];
    }
    Some(Flipper::Bond(bond))
}

/// `std::bernoulli_distribution(0.5)` as implemented by libstdc++:
/// `generate_canonical<double, 53>` (two 32-bit draws) `< 0.5`.
fn random_bool(rng: &mut Mt19937) -> bool {
    let low = f64::from(rng.next_u32());
    let high = f64::from(rng.next_u32());
    let mut value = (low + high * 4_294_967_296.0) / 18_446_744_073_709_551_616.0;
    if value >= 1.0 {
        value = 1.0 - f64::EPSILON / 2.0;
    }
    value < 0.5
}

/// `StereoisomerEnumerator`: acts like a generator of stereoisomers.
#[derive(Debug)]
pub struct StereoisomerEnumerator {
    mol: ROMol,
    options: StereoEnumerationOptions,
    num_returned: u64,
    num_to_return: u64,
    total_possible: u64,
    generated_isomers: HashSet<String>,
    rng: Mt19937,
    flippers: Vec<Flipper>,
    seen: HashSet<Vec<bool>>,
}

impl StereoisomerEnumerator {
    pub fn new(mol: &ROMol, options: StereoEnumerationOptions) -> Result<Self> {
        if options.try_embedding {
            return Err(RDKitError::Value(
                "tryEmbedding is not supported yet: distance geometry is not ported to rdoxide"
                    .into(),
            ));
        }
        let mut mol = mol.clone();
        if mol.num_conformers() > 0 {
            wedge_mol_bonds(&mut mol, None, &BondWedgingParameters::default())?;
        }
        let flippers = Self::build_flippers(&mut mol, &options)?;
        for idx in 0..mol.num_atoms() {
            mol.atom_mut(idx).clear_prop(cp::_CIPCode);
        }
        for idx in 0..mol.num_bonds() {
            if matches!(
                mol.bond(idx).bond_dir(),
                BondDir::EitherDouble | BondDir::Unknown
            ) {
                mol.bond_mut(idx).set_bond_dir(BondDir::None);
            }
        }
        let total_possible = if flippers.is_empty() {
            1
        } else {
            mol.set_prop(cp::_MolFileChiralFlag, "1", false)?;
            u32::try_from(flippers.len())
                .ok()
                .and_then(|n| 1u64.checked_shl(n))
                .unwrap_or(u64::MAX)
        };
        let num_to_return = match options.max_isomers {
            0 => total_possible,
            max => total_possible.min(max),
        };
        let seed = if options.random_seed == -1 {
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.subsec_nanos() ^ (d.as_secs() as u32))
        } else {
            options.random_seed as u32
        };
        Ok(StereoisomerEnumerator {
            mol,
            options,
            num_returned: 0,
            num_to_return,
            total_possible,
            generated_isomers: HashSet::new(),
            rng: Mt19937::new(seed),
            flippers,
            seen: HashSet::new(),
        })
    }

    /// `getStereoisomerCount`
    pub fn stereoisomer_count(&self) -> u64 {
        self.total_possible
    }

    fn build_flippers(mol: &mut ROMol, options: &StereoEnumerationOptions) -> Result<Vec<Flipper>> {
        let mut flippers = Vec::new();
        for si in find_potential_stereo(mol, true, true)? {
            if options.only_unassigned
                && !matches!(
                    si.specified,
                    StereoSpecified::Unknown | StereoSpecified::Unspecified
                )
            {
                continue;
            }
            match si.stereo_type {
                StereoType::AtomTetrahedral => flippers.push(Flipper::Atom(si.centered_on)),
                StereoType::BondDouble => flippers.extend(bond_flipper(mol, &si)),
                StereoType::BondAtropisomer => flippers.push(Flipper::Atropisomer {
                    bond: si.centered_on,
                    ctrl: [
                        si.controlling_atoms[0],
                        si.controlling_atoms[1],
                        si.controlling_atoms[2],
                        si.controlling_atoms[3],
                    ],
                }),
                _ => {}
            }
        }
        if options.only_unassigned {
            // otherwise these would be counted twice
            for group in mol.stereo_groups() {
                if group.group_type() != StereoGroupType::Absolute {
                    flippers.push(Flipper::StereoGroup(
                        group
                            .atoms()
                            .iter()
                            .map(|&a| (a, mol.atom(a).chiral_tag()))
                            .collect(),
                    ));
                }
            }
        }
        Ok(flippers)
    }

    /// Another stereoisomer, or `None` when done.
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Result<Option<ROMol>> {
        if self.num_returned == self.num_to_return {
            return Ok(None);
        }
        if self.flippers.is_empty() {
            self.num_returned += 1;
            return Ok(Some(self.mol.clone()));
        }
        let isomer = self.generate_random_isomer()?;
        self.num_returned += 1;
        Ok(isomer)
    }

    fn generate_random_isomer(&mut self) -> Result<Option<ROMol>> {
        let n = self.flippers.len();
        while (self.seen.len() as u64) < self.total_possible {
            let config: Vec<bool> = (0..n).map(|_| random_bool(&mut self.rng)).collect();
            if !self.seen.insert(config.clone()) {
                continue;
            }
            for (flipper, &flag) in self.flippers.iter().zip(&config) {
                flipper.flip(&mut self.mol, flag)?;
            }
            let mut isomer = self.mol.clone();
            // the stereo groups are not needed any more
            if !isomer.stereo_groups().is_empty() {
                isomer.set_stereo_groups(Vec::new())?;
            }
            set_double_bond_neighbor_directions(&mut isomer, None)?;
            isomer.clear_computed_props(false);
            let atrop_stereo: Vec<(u32, BondStereo)> = (0..isomer.num_bonds())
                .map(|i| (i, isomer.bond(i).stereo()))
                .filter(|(_, s)| {
                    matches!(s, BondStereo::StereoAtropCW | BondStereo::StereoAtropCCW)
                })
                .collect();
            assign_stereochemistry(&mut isomer, true, true, true)?;
            // the flipper updates bond directions but not an existing 2D
            // conformer: keep the configuration it selected
            for (bond, stereo) in atrop_stereo {
                isomer.bond_mut(bond).set_stereo(stereo)?;
            }
            if self.options.unique {
                let smiles = mol_to_cx_smiles(
                    &isomer,
                    &SmilesWriteParams::default(),
                    cx_fields::ALL_BUT_COORDS,
                    RestoreBondDirOption::Clear,
                )?;
                if !self.generated_isomers.insert(smiles) {
                    continue;
                }
            }
            return Ok(Some(isomer));
        }
        Ok(None)
    }
}
