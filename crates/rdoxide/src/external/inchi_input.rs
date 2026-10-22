//! The input side of RDKit's `MolToInchi` (`External/INCHI-API/inchi.cpp`):
//! the "reverse clean-up" of the molecule and its conversion to the InChI
//! API structures `inchi_Atom` / `inchi_Stereo0D` that `GetINCHI` takes.
//!
//! Copyright (c) 2011-2025 Novartis Institutes for BioMedical Research Inc.
//! and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-ChemReactions-Novartis-BSD-3-Clause.txt, the same Novartis
//! license text).
//!
//! One difference from upstream: when a stereocenter has fewer than three
//! explicit neighbors, the unused neighbor slots of its `inchi_Stereo0D` are
//! uninitialised memory in C++; they are zero here.

use crate::Result;
use crate::external::inchi::inchi_dll::{
    INCHI_BOND_STEREO_DOUBLE_EITHER, INCHI_BOND_STEREO_SINGLE_1DOWN,
    INCHI_BOND_STEREO_SINGLE_1EITHER, INCHI_BOND_STEREO_SINGLE_1UP, ISOTOPIC_SHIFT_FLAG, InchiAtom,
};
use crate::external::inchi::inpdef::MAXVAL;
use crate::external::inchi::readinch::{
    INCHI_PARITY_EVEN, INCHI_PARITY_ODD, INCHI_PARITY_UNKNOWN, INCHI_STEREO_TYPE_DOUBLE_BOND,
    INCHI_STEREO_TYPE_TETRAHEDRAL, InchiStereo0D, NO_ATOM,
};
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};

/// `inchi_Input`: the atoms, the 0D stereo and the (fixed) option string.
#[derive(Debug, Clone, PartialEq)]
pub struct InchiInput {
    pub atoms: Vec<InchiAtom>,
    pub stereo_0d: Vec<InchiStereo0D>,
    pub options: Option<String>,
}

/// `fixOptionSymbol`: the InChI option prefix of this platform is `-`, so a
/// `/` in the option string is replaced.
pub fn fix_option_symbol(options: &str) -> String {
    options.replace('/', "-")
}

/// `rCleanUp`: "reverse" clean up, undoing the perchlorate-like clean-up done
/// by RDKit so the molecule can be used with the InChI SDK.
pub fn r_clean_up(mol: &mut RWMol) -> Result<()> {
    let Some(q) = mol_from_smiles("[O-][Cl+3]([O-])([O-])O", &SmilesParserParams::default())?
    else {
        return Ok(());
    };
    let matches = substruct_match(mol, &q, &SubstructMatchParameters::default())?;
    for m in matches {
        let mut map = [0u32; 5];
        for &(q_idx, m_idx) in &m {
            map[q_idx as usize] = m_idx as u32;
        }
        // check charges
        if mol.atom(map[1]).formal_charge() != 3 {
            return Ok(());
        }
        let mut uncharged_found: i32 = -1;
        for i in 0..5 {
            if i == 1 {
                continue;
            }
            if mol.atom(map[i]).formal_charge() == 0 {
                if uncharged_found != -1 {
                    return Ok(()); // too many uncharged oxygen
                }
                uncharged_found = i as i32;
            }
        }
        // flip bonds and remove charges
        for i in 0..5 {
            if i == 1 || i as i32 == uncharged_found {
                continue;
            }
            let b = bond_idx(mol, map[1], map[i])?;
            if uncharged_found == -1 && i == 0 {
                mol.bond_mut(b).set_bond_type(BondType::Single);
                mol.atom_mut(map[i]).set_formal_charge(-1);
                continue;
            }
            mol.bond_mut(b).set_bond_type(BondType::Double);
            mol.atom_mut(map[i]).set_formal_charge(0);
        }
        mol.atom_mut(map[1]).set_formal_charge(0);
    }
    Ok(())
}

fn bond_idx(mol: &ROMol, a: u32, b: u32) -> Result<u32> {
    let idx = mol
        .get_bond_between_atoms(a, b)?
        .map(|bd| bd.idx())
        .ok_or_else(|| crate::RDKitError::Value("no bond between the atoms".into()))?;
    Ok(idx)
}

/// The atoms and bonds half of `MolToInchi`: returns `None` where upstream
/// returns an empty InChI before calling the library (an atom with more
/// than `MAXVAL` bonds).
pub fn mol_to_inchi_input(mol: &ROMol, options: Option<&str>) -> Result<Option<InchiInput>> {
    let mut m = mol.clone();

    // assign stereochem
    if mol.needs_update_property_cache() {
        m.update_property_cache(false)?;
    }
    // kekulize
    kekulize(&mut m, false, false, 100)?;
    // "reverse" cleanup: undo some clean up done by RDKit
    r_clean_up(&mut m)?;

    let n_atoms = m.num_atoms();
    let n_bonds = m.num_bonds();
    let mut atoms = vec![InchiAtom::default(); n_atoms as usize];
    let mut stereo_0d: Vec<InchiStereo0D> = Vec::new();
    let table = PeriodicTable::get();

    for i in 0..n_atoms {
        let iu = i as usize;
        // coordinates
        if m.num_conformers() != 0 {
            let p = m.conformers()[0].atom_pos(i)?;
            atoms[iu].x = p.x;
            atoms[iu].y = p.y;
            atoms[iu].z = p.z;
        }

        // element name
        let atomic_number = m.atom(i).atomic_num();
        let symbol = table.get_element_symbol(atomic_number)?;
        for (d, s) in atoms[iu].elname.iter_mut().zip(symbol.bytes()) {
            *d = s;
        }

        // isotopes
        let isotope = m.atom(i).isotope();
        if isotope != 0 {
            let weight = (table.get_atomic_weight(atomic_number)? + 0.5) as i32;
            atoms[iu].isotopic_mass = (ISOTOPIC_SHIFT_FLAG + isotope as i32 - weight) as i16;
        } else {
            atoms[iu].isotopic_mass = 0;
        }

        // charge
        atoms[iu].charge = m.atom(i).formal_charge() as i8;

        // number of implicit H
        let n_hs: i32 = match atomic_number {
            6 | 7 | 8 | 9 | 17 | 35 | 53 => -1,
            _ => m.total_num_hs(i, false)? as i32,
        };
        atoms[iu].num_iso_h = [n_hs as i8, 0, 0, 0];

        // radical
        atoms[iu].radical = 0;
        if m.atom(i).num_radical_electrons() != 0 {
            // the number of implicit Hs together with the atom identity and
            // charge is sufficient to describe the radical
            atoms[iu].num_iso_h[0] = m.total_num_hs(i, false)? as i8;
        }

        // convert tetrahedral chirality info to Stereo0D
        let tag = m.atom(i).chiral_tag();
        if tag == ChiralType::TetrahedralCCW || tag == ChiralType::TetrahedralCW {
            m.calc_implicit_valence(i, true)?;
            let tval = m.total_degree(i)?;
            if !(3..=4).contains(&tval) {
                crate::rd_log_warning!(
                    "tetrahedral chirality on atom with <3 or >4 neighbors will be ignored.\n"
                );
                continue;
            }
            let mut st = InchiStereo0D {
                neighbor: [0; 4],
                central_atom: i as i16,
                type_: INCHI_STEREO_TYPE_TETRAHEDRAL,
                parity: 0,
            };
            let mut nid = 0usize;
            for nb in m.atom_neighbors(i) {
                if nid < 4 {
                    st.neighbor[nid] = nb as i16;
                }
                nid += 1;
            }
            if nid == 2 {
                // C++ leaves the last two slots uninitialised; the official
                // build accepts the descriptor (with the parity unchanged, as
                // found by trying all placeholders against the library)
                let (n0, n1) = (st.neighbor[0], st.neighbor[1]);
                st.neighbor = [n0, n1, i as i16, n0];
            }
            if nid == 3 {
                for k in (1..=3).rev() {
                    st.neighbor[k] = st.neighbor[k - 1];
                }
                st.neighbor[0] = i as i16;
            }
            // the tag is CW or CCW here, so the descriptor is always used
            let even = if m.atom_degree(i) == 4 {
                tag == ChiralType::TetrahedralCW
            } else {
                tag == ChiralType::TetrahedralCCW
            };
            st.parity = if even {
                INCHI_PARITY_EVEN
            } else {
                INCHI_PARITY_ODD
            } as i8;
            stereo_0d.push(st);
        }
    }

    // read bond info
    for i in 0..n_bonds {
        let (mut a1, mut a2, bond_type, bond_dir, stereo, stereo_atoms) = {
            let b = m.bond(i);
            (
                b.begin_atom_idx(),
                b.end_atom_idx(),
                b.bond_type(),
                b.bond_dir(),
                b.stereo(),
                b.stereo_atoms().to_vec(),
            )
        };
        let mut modifier: i8 = 1;
        // update only for the atom having smaller index
        if a1 > a2 {
            std::mem::swap(&mut a1, &mut a2);
            modifier = -1;
        }

        // neighbor
        let idx = atoms[a1 as usize].num_bonds as usize;
        // the InChI code has a max number of neighbors allowed
        if idx >= MAXVAL {
            crate::rd_log_error!(
                " atom {a1} has too many bonds: {idx}. The InChI library supports at most {MAXVAL}\n"
            );
            return Ok(None);
        }
        atoms[a1 as usize].neighbor[idx] = a2 as i16;

        // bond type
        let mut bt = bond_type;
        if (bt as i32) > BondType::Triple as i32 {
            crate::rd_log_warning!(
                "bond type above 3 ({}) is treated as unspecified!\n",
                bt as i32
            );
            bt = BondType::Unspecified;
        }
        atoms[a1 as usize].bond_type[idx] = bt as i32 as i8;

        // stereo
        atoms[a1 as usize].bond_stereo[idx] = match bond_dir {
            BondDir::BeginWedge => modifier * INCHI_BOND_STEREO_SINGLE_1UP,
            BondDir::BeginDash => modifier * INCHI_BOND_STEREO_SINGLE_1DOWN,
            BondDir::EitherDouble => INCHI_BOND_STEREO_DOUBLE_EITHER,
            BondDir::Unknown => modifier * INCHI_BOND_STEREO_SINGLE_1EITHER,
            _ => 0,
        };

        // double bond stereochemistry (a single bond in a big ring gets E/Z as well)
        if (stereo as i32) > BondStereo::StereoAny as i32 && stereo_atoms.len() >= 2 {
            let mut st = InchiStereo0D {
                neighbor: [
                    stereo_atoms[0] as i16,
                    a1 as i16,
                    a2 as i16,
                    stereo_atoms[1] as i16,
                ],
                central_atom: NO_ATOM,
                type_: INCHI_STEREO_TYPE_DOUBLE_BOND,
                parity: if matches!(stereo, BondStereo::StereoZ | BondStereo::StereoCis) {
                    INCHI_PARITY_ODD
                } else {
                    INCHI_PARITY_EVEN
                } as i8,
            };
            if m.get_bond_between_atoms(st.neighbor[0] as u32, st.neighbor[1] as u32)?
                .is_none()
            {
                st.neighbor.swap(0, 3);
            }
            stereo_0d.push(st);
        } else if stereo == BondStereo::StereoAny {
            // Collapse coordinates so InChI cannot infer stereo from geometry,
            // and send a proper stereo0D with UNKNOWN parity so that -SUU
            // produces the correct unknown annotation. The stereo atoms may
            // be cleared for STEREOANY, so the neighbors are found from the
            // bonds.
            let (x, y, z) = (
                atoms[a2 as usize].x,
                atoms[a2 as usize].y,
                atoms[a2 as usize].z,
            );
            atoms[a1 as usize].x = x;
            atoms[a1 as usize].y = y;
            atoms[a1 as usize].z = z;
            let left = m
                .atom_bonds(a1)
                .filter_map(|nb| nb.other_atom_idx(a1).ok())
                .find(|&o| o != a2);
            let right = m
                .atom_bonds(a2)
                .filter_map(|nb| nb.other_atom_idx(a2).ok())
                .find(|&o| o != a1);
            if let (Some(l), Some(r)) = (left, right) {
                let mut st = InchiStereo0D {
                    neighbor: [l as i16, a1 as i16, a2 as i16, r as i16],
                    central_atom: NO_ATOM,
                    type_: INCHI_STEREO_TYPE_DOUBLE_BOND,
                    parity: INCHI_PARITY_UNKNOWN as i8,
                };
                if m.get_bond_between_atoms(st.neighbor[0] as u32, st.neighbor[1] as u32)?
                    .is_none()
                {
                    st.neighbor.swap(0, 3);
                }
                stereo_0d.push(st);
            }
        }

        // number of bonds
        atoms[a1 as usize].num_bonds += 1;
    }

    Ok(Some(InchiInput {
        atoms,
        stereo_0d,
        options: options.map(fix_option_symbol),
    }))
}
