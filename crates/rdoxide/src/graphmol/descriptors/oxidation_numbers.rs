//! Oxidation numbers (states) of the atoms of a molecule from Pauling
//! electronegativities (upstream `GraphMol/Descriptors/OxidationNumbers.cpp`,
//! Copyright (c) 2023 David Cosgrove, CozChemIx Limited, BSD 3-Clause).
//!
//! The algorithm is based on the `linchemin` project (syngenta/linchemin),
//! which is MIT licensed (see LICENSES/linchemin-MIT.txt).

use crate::graphmol::bond::BondType;
use crate::graphmol::kekulize::{DEFAULT_MAX_BACKTRACKS, kekulize};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::dict::RDProps;
use crate::{RDKitError, Result};

/// Pauling electronegativity by atomic number; elements without an entry
/// count as 0.
const PAULING_EN: &[(u32, f32)] = &[
    (1, 2.2),
    (3, 0.98),
    (4, 1.57),
    (5, 2.04),
    (6, 2.55),
    (7, 3.04),
    (8, 3.44),
    (9, 3.98),
    (11, 0.93),
    (12, 1.31),
    (13, 1.61),
    (14, 1.9),
    (15, 2.19),
    (16, 2.58),
    (17, 3.16),
    (19, 0.82),
    (20, 1.0),
    (21, 1.36),
    (22, 1.54),
    (23, 1.63),
    (24, 1.66),
    (25, 1.55),
    (26, 1.83),
    (27, 1.88),
    (28, 1.91),
    (29, 1.9),
    (30, 1.65),
    (31, 1.81),
    (32, 2.01),
    (33, 2.18),
    (34, 2.55),
    (35, 2.96),
    (36, 3.0),
    (37, 0.82),
    (38, 0.95),
    (39, 1.22),
    (40, 1.33),
    (41, 1.6),
    (42, 2.16),
    (43, 1.9),
    (44, 2.2),
    (45, 2.28),
    (46, 2.2),
    (47, 1.93),
    (48, 1.69),
    (49, 1.78),
    (50, 1.96),
    (51, 2.05),
    (52, 2.1),
    (53, 2.66),
    (54, 2.6),
    (55, 0.79),
    (56, 0.89),
    (57, 1.1),
    (58, 1.12),
    (59, 1.13),
    (60, 1.14),
    (62, 1.17),
    (64, 1.2),
    (66, 1.22),
    (67, 1.23),
    (68, 1.24),
    (69, 1.25),
    (71, 1.27),
    (72, 1.3),
    (73, 1.5),
    (74, 2.36),
    (75, 1.9),
    (76, 2.2),
    (77, 2.2),
    (78, 2.28),
    (79, 2.54),
    (80, 2.0),
    (81, 1.62),
    (82, 2.33),
    (83, 2.02),
    (84, 2.0),
    (85, 2.2),
    (88, 0.9),
    (89, 1.1),
    (90, 1.3),
    (91, 1.5),
    (92, 1.38),
    (93, 1.36),
    (94, 1.28),
    (95, 1.3),
    (96, 1.3),
    (97, 1.3),
    (98, 1.3),
    (99, 1.3),
    (100, 1.3),
    (101, 1.3),
    (102, 1.39),
];

fn electronegativity(atomic_num: u32) -> f32 {
    PAULING_EN
        .binary_search_by_key(&atomic_num, |&(z, _)| z)
        .map_or(0.0, |i| PAULING_EN[i].1)
}

/// -1 when the atom is the more electronegative partner, +1 when the less.
fn sign(en_diff: f32) -> i32 {
    if en_diff > 0.0 {
        -1
    } else if en_diff < 0.0 {
        1
    } else {
        0
    }
}

fn oxidation_number_by_en(mol: &ROMol, idx: u32) -> Result<i32> {
    let atom = mol.atom(idx);
    let own_en = electronegativity(atom.atomic_num());
    let mut ox = 0i32;
    for bond in mol.atom_bonds(idx) {
        if matches!(
            bond.bond_type(),
            BondType::Dative
                | BondType::DativeOne
                | BondType::DativeL
                | BondType::DativeR
                | BondType::Zero
        ) {
            continue;
        }
        let other = mol.atom(bond.other_atom_idx(idx)?);
        if other.atomic_num() > 1 {
            let en_diff = own_en - electronegativity(other.atomic_num());
            let order = bond.bond_type_as_double()?;
            // the molecule must be kekulized: no aromatic 1.5 bonds
            if order > 1.0 && order < 2.0 {
                return Err(RDKitError::Value(
                    "Molecule appears not to be Kekulized, oxidation number calculation fails."
                        .into(),
                ));
            }
            // truncating conversion, as in the upstream `int += double`
            ox = (f64::from(ox) + order * f64::from(sign(en_diff))) as i32;
        }
    }
    ox += sign(own_en - electronegativity(1)) * mol.total_num_hs(idx, false)? as i32;
    Ok(ox + atom.formal_charge())
}

/// `calcOxidationNumbers`: sets the int atom property `OxidationNumber` on
/// every atom (the calculation runs on a kekulized copy).
pub fn calc_oxidation_numbers(mol: &mut ROMol) -> Result<()> {
    let mut kekulized = mol.clone();
    kekulize(&mut kekulized, true, true, DEFAULT_MAX_BACKTRACKS)?;
    for idx in 0..mol.num_atoms() {
        let ox = oxidation_number_by_en(&kekulized, idx)?;
        mol.atom_mut(idx).set_prop(cp::OxidationNumber, ox, false)?;
    }
    Ok(())
}
