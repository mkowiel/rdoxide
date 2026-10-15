//! Construction of MMFF force fields (upstream
//! `ForceFieldHelpers/MMFF/Builder.cpp`).
//!
//! The verbose per-term listings of upstream (`MMFFVerbosity`) are not
//! produced.

use std::sync::OnceLock;

use super::props::{DielectricModel, MmffMolProperties};
use crate::error::{RDKitError, Result};
use crate::forcefield::ForceField;
use crate::forcefield::mmff::terms::{
    AngleBendContrib, BondStretchContrib, DielModel, EleContrib, NonbondedContrib, OopBendContrib,
    StretchBendContrib, TorsionAngleContrib, VdWContrib,
};
use crate::graphmol::atom::HybridizationType;
use crate::graphmol::matrices::get_distance_mat;
use crate::graphmol::mol_ops::frags::get_mol_frags_mapping;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};

/// The smarts of the bonds around which torsions are added.
pub const DEFAULT_TORSION_BOND_SMARTS: &str = "[!$(*#*)&!D1]~[!$(*#*)&!D1]";

const RELATION_1_2: u8 = 0;
const RELATION_1_3: u8 = 1;
const RELATION_1_4: u8 = 2;
const RELATION_1_X: u8 = 3;

/// The default nonbonded cutoff of the force field construction.
pub const DEFAULT_NONBONDED_THRESH: f64 = 100.0;

fn require_valid(props: &MmffMolProperties) -> Result<()> {
    if props.is_valid() {
        Ok(())
    } else {
        Err(RDKitError::precondition(
            "missing atom types - invalid force-field",
        ))
    }
}

fn diel_model(model: DielectricModel) -> DielModel {
    match model {
        DielectricModel::Constant => DielModel::Constant,
        DielectricModel::Distance => DielModel::Distance,
    }
}

fn cell_pos(n: usize, i: usize, j: usize) -> usize {
    let (i, j) = if j < i { (j, i) } else { (i, j) };
    // i * (N - 1) - (i - 1) * i / 2 + j, with the signed arithmetic of upstream
    let (ni, nj) = (i as isize, j as isize);
    (ni * (n as isize - 1) - (ni - 1) * ni / 2 + nj) as usize
}

/// The relation of all atom pairs: directly connected (0), connected through
/// one atom (1), in a 1-4 relationship (2) or farther apart (3).
pub fn build_neighbor_matrix(mol: &ROMol) -> Vec<u8> {
    let n = mol.num_atoms() as usize;
    let mut res = vec![RELATION_1_X; n * (n + 1) / 2 + 1];
    let dmat = get_distance_mat(mol);
    for i in 0..n {
        res[cell_pos(n, i, i)] = RELATION_1_X;
        for j in i + 1..n {
            let d = dmat[i * n + j];
            res[cell_pos(n, i, j)] = if d == 1.0 {
                RELATION_1_2
            } else if d == 2.0 {
                RELATION_1_3
            } else if d == 3.0 {
                RELATION_1_4
            } else {
                RELATION_1_X
            };
        }
    }
    res
}

fn push_if_any<C: crate::forcefield::ForceFieldContrib + 'static>(
    field: &mut ForceField,
    contrib: C,
    any: bool,
) {
    if any {
        field.push_contrib(Box::new(contrib));
    }
}

pub fn add_bonds(mol: &ROMol, props: &MmffMolProperties, field: &mut ForceField) -> Result<()> {
    require_valid(props)?;
    let mut contrib = BondStretchContrib::new();
    for bond in mol.bonds() {
        let (idx1, idx2) = (bond.begin_atom_idx(), bond.end_atom_idx());
        if let Some(p) = props.bond_stretch_params(mol, idx1, idx2)? {
            contrib.add_term(idx1 as usize, idx2 as usize, &p.params);
        }
    }
    let any = !contrib.is_empty();
    push_if_any(field, contrib, any);
    Ok(())
}

pub fn add_angles(mol: &ROMol, props: &MmffMolProperties, field: &mut ForceField) -> Result<()> {
    require_valid(props)?;
    let mut contrib = AngleBendContrib::new();
    for j in 0..mol.num_atoms() {
        if mol.atom_degree(j) == 1 {
            continue;
        }
        let central = props.prop(j);
        let nbrs: Vec<u32> = mol.atom_neighbors(j).collect();
        for (a, &i) in nbrs.iter().enumerate() {
            for &k in &nbrs[a + 1..] {
                if let Some(p) = props.angle_bend_params(mol, i, j, k)? {
                    contrib.add_term([i as usize, j as usize, k as usize], &p.params, central);
                }
            }
        }
    }
    let any = !contrib.is_empty();
    push_if_any(field, contrib, any);
    Ok(())
}

pub fn add_stretch_bend(
    mol: &ROMol,
    props: &MmffMolProperties,
    field: &mut ForceField,
) -> Result<()> {
    require_valid(props)?;
    let mut contrib = StretchBendContrib::new();
    for j in 0..mol.num_atoms() {
        if mol.atom_degree(j) == 1 || props.prop(j).linh != 0 {
            continue;
        }
        let nbrs: Vec<u32> = mol.atom_neighbors(j).collect();
        for (a, &i) in nbrs.iter().enumerate() {
            for &k in &nbrs[a + 1..] {
                if let Some(p) = props.stretch_bend_params(mol, i, j, k)? {
                    contrib.add_term(
                        [i as usize, j as usize, k as usize],
                        &p.params,
                        &p.angle,
                        &p.bonds[0],
                        &p.bonds[1],
                    );
                }
            }
        }
    }
    let any = !contrib.is_empty();
    push_if_any(field, contrib, any);
    Ok(())
}

pub fn add_oop(mol: &ROMol, props: &MmffMolProperties, field: &mut ForceField) -> Result<()> {
    require_valid(props)?;
    let mut contrib = OopBendContrib::new();
    for center in 0..mol.num_atoms() {
        if mol.atom_degree(center) != 3 {
            continue;
        }
        let nbrs: Vec<u32> = mol.atom_neighbors(center).collect();
        // the neighbors are i, k and l around the center j
        let idx = [nbrs[0], center, nbrs[1], nbrs[2]];
        // without parameters the term is left out (SURDOX02)
        let Some(p) = props.oop_bend_params(mol, idx[0], idx[1], idx[2], idx[3]) else {
            continue;
        };
        for n in [[0, 1, 2, 3], [0, 1, 3, 2], [2, 1, 3, 0]] {
            contrib.add_term(n.map(|k| idx[k] as usize), &p);
        }
    }
    let any = !contrib.is_empty();
    push_if_any(field, contrib, any);
    Ok(())
}

pub fn add_torsions(
    mol: &ROMol,
    props: &MmffMolProperties,
    field: &mut ForceField,
    torsion_bond_smarts: &str,
) -> Result<()> {
    require_valid(props)?;
    static DEFAULT_QUERY: OnceLock<ROMol> = OnceLock::new();
    let owned;
    let query = if torsion_bond_smarts == DEFAULT_TORSION_BOND_SMARTS {
        DEFAULT_QUERY.get_or_init(|| {
            mol_from_smarts(DEFAULT_TORSION_BOND_SMARTS, &SmartsParserParams::default())
                .ok()
                .flatten()
                .expect("the default torsion bond smarts is valid")
        })
    } else {
        owned = mol_from_smarts(torsion_bond_smarts, &SmartsParserParams::default())?
            .ok_or_else(|| RDKitError::Value("bad torsion bond SMARTS".into()))?;
        &owned
    };
    let matches = substruct_match(mol, query, &SubstructMatchParameters::default())?;
    let mut contrib = TorsionAngleContrib::new();
    let is_sp2_or_sp3 = |idx: u32| {
        matches!(
            mol.atom(idx).hybridization(),
            HybridizationType::SP2 | HybridizationType::SP3
        )
    };
    for m in matches {
        if m.len() != 2 {
            return Err(RDKitError::invariant(
                "torsion bond query must have two atoms",
            ));
        }
        let (idx2, idx3) = (m[0].1 as u32, m[1].1 as u32);
        let bond = mol
            .get_bond_between_atoms(idx2, idx3)?
            .ok_or_else(|| RDKitError::invariant("torsion bond not found"))?;
        if !(is_sp2_or_sp3(idx2) && is_sp2_or_sp3(idx3)) {
            continue;
        }
        for t_bond1 in mol.atom_bonds(idx2) {
            if t_bond1.idx() == bond.idx() {
                continue;
            }
            let idx1 = t_bond1.other_atom_idx(idx2)?;
            for t_bond2 in mol.atom_bonds(idx3) {
                if t_bond2.idx() == bond.idx() || t_bond2.idx() == t_bond1.idx() {
                    continue;
                }
                let idx4 = t_bond2.other_atom_idx(idx3)?;
                // not a three-membered ring
                if idx4 == idx1 {
                    continue;
                }
                if let Some(p) = props.torsion_params(mol, idx1, idx2, idx3, idx4)? {
                    contrib.add_term(
                        [idx1 as usize, idx2 as usize, idx3 as usize, idx4 as usize],
                        &p.params,
                    );
                }
            }
        }
    }
    let any = !contrib.is_empty();
    push_if_any(field, contrib, any);
    Ok(())
}

fn frag_mapping(mol: &ROMol, ignore_interfrag: bool) -> Vec<i32> {
    if ignore_interfrag {
        get_mol_frags_mapping(mol).1
    } else {
        Vec::new()
    }
}

/// Van der Waals terms between atoms that are at least 1-4 related and
/// closer than `non_bonded_thresh`.
pub fn add_vdw(
    mol: &ROMol,
    conf_id: i32,
    props: &MmffMolProperties,
    field: &mut ForceField,
    neighbor_matrix: &[u8],
    non_bonded_thresh: f64,
    ignore_interfrag: bool,
) -> Result<()> {
    require_valid(props)?;
    let frags = frag_mapping(mol, ignore_interfrag);
    let n = mol.num_atoms() as usize;
    let conf = mol.conformer(conf_id)?;
    let mut contrib = VdWContrib::new();
    for i in 0..n {
        for j in i + 1..n {
            if ignore_interfrag && frags[i] != frags[j] {
                continue;
            }
            if neighbor_matrix[cell_pos(n, i, j)] < RELATION_1_4 {
                continue;
            }
            let d = (conf.atom_pos(i as u32)? - conf.atom_pos(j as u32)?).length();
            if d > non_bonded_thresh {
                continue;
            }
            if let Some(c) = props.vdw_params(i as u32, j as u32) {
                contrib.add_term(i, j, &c);
            }
        }
    }
    let any = !contrib.is_empty();
    push_if_any(field, contrib, any);
    Ok(())
}

fn charge_term(props: &MmffMolProperties, i: u32, j: u32) -> Option<f64> {
    let (qi, qj) = (props.partial_charge(i), props.partial_charge(j));
    if qi.abs() < 1.0e-10 || qj.abs() < 1.0e-10 {
        None
    } else {
        Some(qi * qj / props.dielectric_constant())
    }
}

/// Electrostatic terms between charged atoms that are at least 1-4 related
/// and closer than `non_bonded_thresh`.
pub fn add_ele(
    mol: &ROMol,
    conf_id: i32,
    props: &MmffMolProperties,
    field: &mut ForceField,
    neighbor_matrix: &[u8],
    non_bonded_thresh: f64,
    ignore_interfrag: bool,
) -> Result<()> {
    require_valid(props)?;
    let frags = frag_mapping(mol, ignore_interfrag);
    let n = mol.num_atoms() as usize;
    let conf = mol.conformer(conf_id)?;
    let model = diel_model(props.dielectric_model());
    let mut contrib = EleContrib::new();
    for i in 0..n {
        for j in i + 1..n {
            if ignore_interfrag && frags[i] != frags[j] {
                continue;
            }
            let cell = neighbor_matrix[cell_pos(n, i, j)];
            if cell < RELATION_1_4 {
                continue;
            }
            let Some(q) = charge_term(props, i as u32, j as u32) else {
                continue;
            };
            let d = (conf.atom_pos(i as u32)? - conf.atom_pos(j as u32)?).length();
            if d > non_bonded_thresh {
                continue;
            }
            contrib.add_term(i, j, q, model, cell == RELATION_1_4);
        }
    }
    let any = !contrib.is_empty();
    push_if_any(field, contrib, any);
    Ok(())
}

/// The van der Waals and electrostatic terms in one contribution, honoring
/// the term switches of the properties.
pub fn add_nonbonded(
    mol: &ROMol,
    conf_id: i32,
    props: &MmffMolProperties,
    field: &mut ForceField,
    neighbor_matrix: &[u8],
    non_bonded_thresh: f64,
    ignore_interfrag: bool,
) -> Result<()> {
    require_valid(props)?;
    let frags = frag_mapping(mol, ignore_interfrag);
    let n = mol.num_atoms() as usize;
    let conf = mol.conformer(conf_id)?;
    let model = diel_model(props.dielectric_model());
    let mut contrib = NonbondedContrib::new();
    for i in 0..n {
        for j in i + 1..n {
            if ignore_interfrag && frags[i] != frags[j] {
                continue;
            }
            let cell = neighbor_matrix[cell_pos(n, i, j)];
            if cell < RELATION_1_4 {
                continue;
            }
            let d = (conf.atom_pos(i as u32)? - conf.atom_pos(j as u32)?).length();
            if d > non_bonded_thresh {
                continue;
            }
            let vdw = if props.vdw_term() {
                props.vdw_params(i as u32, j as u32)
            } else {
                None
            };
            let ele = if props.ele_term() {
                charge_term(props, i as u32, j as u32).map(|q| (q, model, cell == RELATION_1_4))
            } else {
                None
            };
            contrib.add_term(i, j, vdw.as_ref(), ele);
        }
    }
    let any = !contrib.is_empty();
    push_if_any(field, contrib, any);
    Ok(())
}

/// `constructForceField(mol, mmffMolProperties, nonBondedThresh, confId,
/// ignoreInterfragInteractions)`
pub fn construct_force_field_with_props(
    mol: &ROMol,
    props: &MmffMolProperties,
    non_bonded_thresh: f64,
    conf_id: i32,
    ignore_interfrag: bool,
) -> Result<ForceField> {
    require_valid(props)?;
    let mut field = ForceField::new(3);
    let conf = mol.conformer(conf_id)?;
    for i in 0..mol.num_atoms() {
        let p = conf.atom_pos(i)?;
        field.push_point(&[p.x, p.y, p.z]);
    }
    if props.bond_term() {
        add_bonds(mol, props, &mut field)?;
    }
    if props.angle_term() {
        add_angles(mol, props, &mut field)?;
    }
    if props.stretch_bend_term() {
        add_stretch_bend(mol, props, &mut field)?;
    }
    if props.oop_term() {
        add_oop(mol, props, &mut field)?;
    }
    if props.torsion_term() {
        add_torsions(mol, props, &mut field, DEFAULT_TORSION_BOND_SMARTS)?;
    }
    if props.vdw_term() || props.ele_term() {
        let neighbor_matrix = build_neighbor_matrix(mol);
        add_nonbonded(
            mol,
            conf_id,
            props,
            &mut field,
            &neighbor_matrix,
            non_bonded_thresh,
            ignore_interfrag,
        )?;
    }
    Ok(field)
}

/// `constructForceField(mol, nonBondedThresh, confId,
/// ignoreInterfragInteractions)`: types the atoms of `mol` (which gets MMFF
/// aromaticity) with the default MMFF94 settings.
pub fn construct_force_field(
    mol: &mut ROMol,
    non_bonded_thresh: f64,
    conf_id: i32,
    ignore_interfrag: bool,
) -> Result<ForceField> {
    let props = MmffMolProperties::new(mol, "MMFF94")?;
    require_valid(&props)?;
    construct_force_field_with_props(mol, &props, non_bonded_thresh, conf_id, ignore_interfrag)
}
