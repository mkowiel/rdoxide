//! UFF atom typing and force-field construction (upstream
//! `ForceFieldHelpers/UFF/AtomTyper.cpp` and `Builder.cpp`).

use std::sync::OnceLock;

use crate::error::{RDKitError, Result};
use crate::forcefield::ForceField;
use crate::forcefield::uff::{
    AngleBendContrib, AtomicParams, BondStretchContrib, InversionContrib, ParamCollection, RAD2DEG,
    TorsionAngleContrib, TorsionHybridization, VdwContrib, calc_angle_force_constant,
    calc_bond_force_constant, calc_bond_rest_length,
    calc_inversion_coefficients_and_force_constant, calc_nonbonded_depth, calc_nonbonded_minimum,
    equation17, is_in_group6,
};
use crate::graphmol::atom::HybridizationType;
use crate::graphmol::conjug_hybrid::atom_has_conjugated_bond;
use crate::graphmol::mol_ops::add_hs::needs_hs;
use crate::graphmol::mol_ops::frags::get_mol_frags_mapping;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::{rd_log_error, rd_log_warning};

/// The parameters of every atom of a molecule (`None`: no parameters).
pub type AtomicParamVect = Vec<Option<AtomicParams>>;

/// `UFFBond`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UffBond {
    pub kb: f64,
    pub r0: f64,
}

/// `UFFAngle`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UffAngle {
    pub ka: f64,
    pub theta0: f64,
}

/// The smarts of the bonds around which torsions are added.
pub const DEFAULT_TORSION_BOND_SMARTS: &str = "[!$(*#*)&!D1]~[!$(*#*)&!D1]";

fn unrecognized_charge(idx: u32) {
    rd_log_error!("UFFTYPER: Unrecognized charge state for atom: {idx}");
}

/// `Tools::addAtomChargeFlags`: appends the charge flag of the atom type.
pub fn add_atom_charge_flags(
    mol: &ROMol,
    idx: u32,
    key: &mut String,
    tolerate_charge_mismatch: bool,
) -> Result<()> {
    let atom = mol.atom(idx);
    let total_valence = atom.total_valence()? as i32;
    let fc = atom.formal_charge();
    let num = atom.atomic_num();
    // elements that only come with one charge in the default UFF parameters
    let fixed_charge = match num {
        29 | 47 => Some(1),
        4 | 20 | 25 | 26 | 28 | 46 | 78 => Some(2),
        21 | 24 | 27 | 79 | 89 | 96..=103 => Some(3),
        2 | 18 | 22 | 36 | 54 | 90..=95 => Some(4),
        23 | 41 | 43 | 73 => Some(5),
        42 => Some(6),
        _ => None,
    };
    if let Some(n) = fixed_charge {
        if total_valence == n || fc == n || tolerate_charge_mismatch {
            key.push_str(&format!("+{n}"));
        } else {
            unrecognized_charge(idx);
        }
    }
    // elements whose charge flag follows the total valence; the default
    // applies when the mismatch is tolerated
    let by_valence: Option<(&[(i32, &str)], &str)> = match num {
        12 | 30 | 34 | 48 | 52 | 80 | 84 => Some((&[(2, "+2")], "+2")),
        15 => Some((&[(3, "+3"), (5, "+5")], "+5")),
        16 if atom.hybridization() != HybridizationType::SP2 => {
            Some((&[(2, "+2"), (4, "+4"), (6, "+6")], "+6"))
        }
        31 | 33 | 49 | 51 | 81 | 82 | 83 => Some((&[(3, "+3")], "+3")),
        _ => None,
    };
    if let Some((table, default)) = by_valence {
        match table.iter().find(|(v, _)| *v == total_valence) {
            Some((_, flag)) => key.push_str(flag),
            None => {
                if tolerate_charge_mismatch {
                    key.push_str(default);
                }
                unrecognized_charge(idx);
            }
        }
    }
    match num {
        13 if total_valence != 3 => unrecognized_charge(idx),
        14 if total_valence != 4 => unrecognized_charge(idx),
        75 => {
            if tolerate_charge_mismatch {
                if key == "Re6" {
                    *key = "Re6+5".into();
                } else if key == "Re3" {
                    *key = "Re3+7".into();
                }
            }
            unrecognized_charge(idx);
        }
        _ => {}
    }
    // lanthanides
    if (57..=71).contains(&num) {
        if total_valence == 6 {
            key.push_str("+3");
        } else {
            if tolerate_charge_mismatch {
                key.push_str("+3");
            }
            unrecognized_charge(idx);
        }
    }
    Ok(())
}

/// `Tools::getAtomLabel`: the UFF atom type label of an atom.
pub fn get_atom_label(mol: &ROMol, idx: u32) -> Result<String> {
    let atom = mol.atom(idx);
    let num = atom.atomic_num();
    let mut key = atom.symbol();
    if key.len() == 1 {
        key.push('_');
    }
    let table = PeriodicTable::get();
    // hybridization is not added to alkali metals or halogens
    if num != 0 {
        let outer = table.get_n_outer_elecs(num)?;
        if table.get_default_valence(num)? == -1 || (outer != 1 && outer != 7) {
            let hyb = atom.hybridization();
            match num {
                12 | 13 | 14 | 15 | 50 | 51 | 52 | 81 | 82 | 83 | 84 => {
                    key.push('3');
                    if hyb != HybridizationType::SP3 {
                        rd_log_warning!(
                            "UFFTYPER: Warning: hybridization set to SP3 for atom {idx}"
                        );
                    }
                }
                80 => {
                    key.push('1');
                    if hyb != HybridizationType::SP {
                        rd_log_warning!(
                            "UFFTYPER: Warning: hybridization set to SP for atom {idx}"
                        );
                    }
                }
                _ => match hyb {
                    HybridizationType::S => {}
                    HybridizationType::SP => key.push('1'),
                    HybridizationType::SP2 => {
                        if (atom.is_aromatic() || atom_has_conjugated_bond(mol, idx))
                            && matches!(num, 6 | 7 | 8 | 16)
                        {
                            key.push('R');
                        } else {
                            key.push('2');
                        }
                    }
                    HybridizationType::SP3 => key.push('3'),
                    HybridizationType::SP2D => key.push('4'),
                    HybridizationType::SP3D => key.push('5'),
                    HybridizationType::SP3D2 => key.push('6'),
                    _ => rd_log_error!("UFFTYPER: Unrecognized hybridization for atom: {idx}"),
                },
            }
        }
    }
    add_atom_charge_flags(mol, idx, &mut key, true)?;
    Ok(key)
}

fn params_for(mol: &ROMol, idx: u32) -> Result<Option<AtomicParams>> {
    let key = get_atom_label(mol, idx)?;
    Ok(ParamCollection::default_params().get(&key).copied())
}

/// `getAtomTypes`: the parameters of every atom and whether all were found.
pub fn get_atom_types(mol: &ROMol) -> Result<(AtomicParamVect, bool)> {
    let mut found_all = true;
    let mut params = Vec::with_capacity(mol.num_atoms() as usize);
    for i in 0..mol.num_atoms() {
        let key = get_atom_label(mol, i)?;
        let p = ParamCollection::default_params().get(&key).copied();
        if p.is_none() {
            found_all = false;
            rd_log_error!("UFFTYPER: Unrecognized atom type: {key} ({i})");
        }
        params.push(p);
    }
    Ok((params, found_all))
}

/// `getUFFBondStretchParams`
pub fn get_uff_bond_stretch_params(mol: &ROMol, idx1: u32, idx2: u32) -> Result<Option<UffBond>> {
    let Some(bond) = mol.get_bond_between_atoms(idx1, idx2)? else {
        return Ok(None);
    };
    let (Some(p1), Some(p2)) = (params_for(mol, idx1)?, params_for(mol, idx2)?) else {
        return Ok(None);
    };
    let r0 = calc_bond_rest_length(bond.bond_type_as_double()?, &p1, &p2);
    Ok(Some(UffBond {
        kb: calc_bond_force_constant(r0, &p1, &p2),
        r0,
    }))
}

/// `getUFFAngleBendParams`
pub fn get_uff_angle_bend_params(
    mol: &ROMol,
    idx1: u32,
    idx2: u32,
    idx3: u32,
) -> Result<Option<UffAngle>> {
    let (Some(b1), Some(b2)) = (
        mol.get_bond_between_atoms(idx1, idx2)?,
        mol.get_bond_between_atoms(idx2, idx3)?,
    ) else {
        return Ok(None);
    };
    let (Some(p1), Some(p2), Some(p3)) = (
        params_for(mol, idx1)?,
        params_for(mol, idx2)?,
        params_for(mol, idx3)?,
    ) else {
        return Ok(None);
    };
    let ka = calc_angle_force_constant(
        p2.theta0,
        b1.bond_type_as_double()?,
        b2.bond_type_as_double()?,
        &p1,
        &p2,
        &p3,
    );
    Ok(Some(UffAngle {
        ka,
        theta0: RAD2DEG * p2.theta0,
    }))
}

/// `getUFFTorsionParams`: the barrier V
pub fn get_uff_torsion_params(
    mol: &ROMol,
    idx1: u32,
    idx2: u32,
    idx3: u32,
    idx4: u32,
) -> Result<Option<f64>> {
    let idx = [idx1, idx2, idx3, idx4];
    for w in idx.windows(2) {
        if mol.get_bond_between_atoms(w[0], w[1])?.is_none() {
            return Ok(None);
        }
    }
    let bond = mol
        .get_bond_between_atoms(idx2, idx3)?
        .expect("checked above");
    let mut has_sp2 = false;
    for &i in &[idx1, idx4] {
        if mol.atom(i).hybridization() == HybridizationType::SP2 {
            has_sp2 = true;
        }
    }
    let (Some(p2), Some(p3)) = (params_for(mol, idx2)?, params_for(mol, idx3)?) else {
        return Ok(None);
    };
    let (n2, n3) = (mol.atom(idx2).atomic_num(), mol.atom(idx3).atomic_num());
    let (h2, h3) = (
        mol.atom(idx2).hybridization(),
        mol.atom(idx3).hybridization(),
    );
    use HybridizationType::{SP2, SP3};
    if !((h2 == SP2 || h2 == SP3) && (h3 == SP2 || h3 == SP3)) {
        return Ok(None);
    }
    let bond_order = bond.bond_type_as_double()?;
    let single = (bond_order * 10.0) as i32 == 10;
    let v = if h2 == SP3 && h3 == SP3 {
        let mut vv = (p2.v1 * p3.v1).sqrt();
        if single && is_in_group6(n2) && is_in_group6(n3) {
            let (mut v2, mut v3) = (6.8f64, 6.8f64);
            if n2 == 8 {
                v2 = 2.0;
            }
            if n3 == 8 {
                v3 = 2.0;
            }
            vv = (v2 * v3).sqrt();
        }
        vv
    } else if h2 == SP2 && h3 == SP2 {
        equation17(bond_order, &p2, &p3)
    } else {
        let mut vv = 1.0;
        if single {
            if (h2 == SP3 && is_in_group6(n2) && !is_in_group6(n3))
                || (h3 == SP3 && is_in_group6(n3) && !is_in_group6(n2))
            {
                vv = equation17(bond_order, &p2, &p3);
            } else if has_sp2 {
                vv = 2.0;
            }
        }
        vv
    };
    Ok(Some(v))
}

/// `getUFFInversionParams`: the force constant K
pub fn get_uff_inversion_params(
    mol: &ROMol,
    idx1: u32,
    idx2: u32,
    idx3: u32,
    idx4: u32,
) -> Result<Option<f64>> {
    if mol.get_bond_between_atoms(idx1, idx2)?.is_none()
        || mol.get_bond_between_atoms(idx2, idx3)?.is_none()
        || mol.get_bond_between_atoms(idx2, idx4)?.is_none()
    {
        return Ok(None);
    }
    let center = mol.atom(idx2);
    let num = center.atomic_num();
    // the central atom must be carbon, nitrogen, oxygen, phosphorous,
    // arsenic, antimony or bismuth with three neighbors
    if !matches!(num, 6 | 7 | 8 | 15 | 33 | 51 | 83) || mol.atom_degree(idx2) != 3 {
        return Ok(None);
    }
    if (6..=8).contains(&num) && center.hybridization() != HybridizationType::SP2 {
        return Ok(None);
    }
    let mut bound_to_sp2_o = false;
    for &i in &[idx1, idx3, idx4] {
        let a = mol.atom(i);
        if a.atomic_num() == 8 && a.hybridization() == HybridizationType::SP2 {
            bound_to_sp2_o = true;
        }
    }
    let (k, ..) = calc_inversion_coefficients_and_force_constant(num, bound_to_sp2_o && num == 6);
    Ok(Some(k))
}

/// `getUFFVdWParams`: (x_ij, D_ij)
pub fn get_uff_vdw_params(mol: &ROMol, idx1: u32, idx2: u32) -> Result<Option<(f64, f64)>> {
    let (Some(p1), Some(p2)) = (params_for(mol, idx1)?, params_for(mol, idx2)?) else {
        return Ok(None);
    };
    Ok(Some((
        calc_nonbonded_minimum(&p1, &p2),
        calc_nonbonded_depth(&p1, &p2),
    )))
}

// ------------------------------------------------------------------ builder

const RELATION_1_2: u8 = 0;
const RELATION_1_3: u8 = 1;
const RELATION_1_4: u8 = 2;
const RELATION_1_X: u8 = 3;

/// upper triangle position of the pair (i, j)
fn cell_pos(n_atoms: usize, i: usize, j: usize) -> usize {
    let (i, j) = if j < i { (j, i) } else { (i, j) };
    // i * (n - 1) + i * (1 - i) / 2 + j, evaluated with signed integers
    ((i as i64) * (n_atoms as i64 - 1) + (i as i64) * (1 - i as i64) / 2 + j as i64) as usize
}

/// 1-2, 1-3 or further relations of all atom pairs (`buildNeighborMatrix`).
pub fn build_neighbor_matrix(mol: &ROMol) -> Vec<u8> {
    let n = mol.num_atoms() as usize;
    let mut res = vec![RELATION_1_X; n * (n + 1) / 2 + 1];
    for bi in mol.bonds() {
        let (ib, ie) = (bi.begin_atom_idx() as usize, bi.end_atom_idx() as usize);
        res[cell_pos(n, ib, ie)] = RELATION_1_2;
        for bj in mol.bonds() {
            if bj.idx() == bi.idx() {
                continue;
            }
            let (jb, je) = (bj.begin_atom_idx() as usize, bj.end_atom_idx() as usize);
            let pair = if ib == jb {
                Some((ie, je))
            } else if ib == je {
                Some((ie, jb))
            } else if ie == jb {
                Some((ib, je))
            } else if ie == je {
                Some((ib, jb))
            } else {
                None
            };
            if let Some((a, b)) = pair {
                res[cell_pos(n, a, b)] = RELATION_1_3;
            }
        }
    }
    res
}

fn add_bonds(mol: &ROMol, params: &AtomicParamVect, field: &mut ForceField) -> Result<()> {
    for bond in mol.bonds() {
        let (i1, i2) = (bond.begin_atom_idx() as usize, bond.end_atom_idx() as usize);
        if let (Some(p1), Some(p2)) = (&params[i1], &params[i2]) {
            field.push_contrib(Box::new(BondStretchContrib::new(
                i1,
                i2,
                bond.bond_type_as_double()?,
                p1,
                p2,
            )));
        }
    }
    Ok(())
}

fn add_angles(mol: &ROMol, params: &AtomicParamVect, field: &mut ForceField) -> Result<()> {
    let rings = mol.ring_info();
    for j in 0..mol.num_atoms() {
        let Some(pj) = &params[j as usize] else {
            continue;
        };
        if mol.atom_degree(j) == 1 {
            continue;
        }
        let atom_j = mol.atom(j);
        let nbrs: Vec<u32> = mol.atom_neighbors(j).collect();
        for (a, &i) in nbrs.iter().enumerate() {
            let Some(pi) = &params[i as usize] else {
                continue;
            };
            for &k in &nbrs[a + 1..] {
                let Some(pk) = &params[k as usize] else {
                    continue;
                };
                // skip trigonal bipyramids
                if atom_j.hybridization() == HybridizationType::SP3D && mol.atom_degree(j) == 5 {
                    continue;
                }
                let b1 = mol
                    .get_bond_between_atoms(i, j)?
                    .expect("neighbors are bonded");
                let b2 = mol
                    .get_bond_between_atoms(k, j)?
                    .expect("neighbors are bonded");
                let order = match atom_j.hybridization() {
                    HybridizationType::SP => 1,
                    HybridizationType::SP2 => {
                        let mut order = 3;
                        // a hack for decent geometries of 3- and 4-membered
                        // rings with sp2 atoms: an angle between a ring
                        // substituent and a ring edge (30, 40), or a ring
                        // angle (35, 45)
                        if rings.is_atom_in_ring_of_size(j, 3)? {
                            let (ri, rk) = (
                                rings.is_atom_in_ring_of_size(i, 3)?,
                                rings.is_atom_in_ring_of_size(k, 3)?,
                            );
                            if ri != rk {
                                order = 30;
                            } else if ri && rk {
                                order = 35;
                            }
                        } else if rings.is_atom_in_ring_of_size(j, 4)? {
                            let (ri, rk) = (
                                rings.is_atom_in_ring_of_size(i, 4)?,
                                rings.is_atom_in_ring_of_size(k, 4)?,
                            );
                            if ri != rk {
                                order = 40;
                            } else if ri && rk {
                                order = 45;
                            }
                        }
                        order
                    }
                    HybridizationType::SP3D2 => 4,
                    _ => 0,
                };
                field.push_contrib(Box::new(AngleBendContrib::new(
                    i as usize,
                    j as usize,
                    k as usize,
                    b1.bond_type_as_double()?,
                    b2.bond_type_as_double()?,
                    pi,
                    pj,
                    pk,
                    order,
                )));
            }
        }
    }
    Ok(())
}

/// The angles around a trigonal-bipyramidal atom.
fn add_trigonal_bipyramid_angles(
    mol: &ROMol,
    atom_idx: u32,
    conf_id: i32,
    params: &AtomicParamVect,
    field: &mut ForceField,
) -> Result<()> {
    let conf = mol.conformer(conf_id)?;
    // the axial bonds are the two with the most negative dot product
    let bonds: Vec<&crate::graphmol::bond::Bond> = mol.atom_bonds(atom_idx).collect();
    let center = conf.atom_pos(atom_idx)?;
    let mut most_neg = 100.0;
    let (mut ax1, mut ax2) = (None, None);
    for b1 in &bonds {
        let v1 = center.direction_vector(&conf.atom_pos(b1.other_atom_idx(atom_idx)?)?)?;
        for b2 in &bonds {
            if b2.idx() > b1.idx() {
                let v2 = center.direction_vector(&conf.atom_pos(b2.other_atom_idx(atom_idx)?)?)?;
                let dot = v1.dot_product(&v2);
                if dot < most_neg {
                    most_neg = dot;
                    ax1 = Some(*b1);
                    ax2 = Some(*b2);
                }
            }
        }
    }
    let (Some(ax1), Some(ax2)) = (ax1, ax2) else {
        return Err(RDKitError::invariant("axial bond not found"));
    };
    let eq: Vec<&crate::graphmol::bond::Bond> = bonds
        .iter()
        .copied()
        .filter(|b| b.idx() != ax1.idx() && b.idx() != ax2.idx())
        .collect();
    if eq.len() < 3 {
        return Err(RDKitError::invariant("equatorial bond not found"));
    }
    let (eq1, eq2, eq3) = (eq[0], eq[1], eq[2]);
    let pa = params[atom_idx as usize].as_ref();
    let mut add = |b1: &crate::graphmol::bond::Bond,
                   b2: &crate::graphmol::bond::Bond,
                   order: u32|
     -> Result<()> {
        let i = b1.other_atom_idx(atom_idx)? as usize;
        let j = b2.other_atom_idx(atom_idx)? as usize;
        if let (Some(pi), Some(pj), Some(pa)) = (&params[i], &params[j], pa) {
            field.push_contrib(Box::new(AngleBendContrib::new(
                i,
                atom_idx as usize,
                j,
                b1.bond_type_as_double()?,
                b2.bond_type_as_double()?,
                pi,
                pa,
                pj,
                order,
            )));
        }
        Ok(())
    };
    // axial-axial
    add(ax1, ax2, 2)?;
    // equatorial-equatorial
    add(eq1, eq2, 3)?;
    add(eq1, eq3, 3)?;
    add(eq2, eq3, 3)?;
    // axial-equatorial
    for ax in [ax1, ax2] {
        for e in [eq1, eq2, eq3] {
            add(ax, e, 0)?;
        }
    }
    Ok(())
}

fn add_nonbonded(
    mol: &ROMol,
    conf_id: i32,
    params: &AtomicParamVect,
    field: &mut ForceField,
    neighbor_matrix: &[u8],
    vdw_thresh: f64,
    ignore_interfrag: bool,
) -> Result<()> {
    let frag_mapping = if ignore_interfrag {
        get_mol_frags_mapping(mol).1
    } else {
        Vec::new()
    };
    let n = mol.num_atoms() as usize;
    let conf = mol.conformer(conf_id)?;
    for i in 0..n {
        let Some(pi) = &params[i] else { continue };
        for j in i + 1..n {
            let Some(pj) = &params[j] else { continue };
            if ignore_interfrag && frag_mapping[i] != frag_mapping[j] {
                continue;
            }
            if neighbor_matrix[cell_pos(n, i, j)] >= RELATION_1_4 {
                let d = (conf.atom_pos(i as u32)? - conf.atom_pos(j as u32)?).length();
                if d < vdw_thresh * calc_nonbonded_minimum(pi, pj) {
                    field.push_contrib(Box::new(VdwContrib::new(i, j, pi, pj, 10.0)));
                }
            }
        }
    }
    Ok(())
}

fn torsion_hybridization(h: HybridizationType) -> Option<TorsionHybridization> {
    match h {
        HybridizationType::SP2 => Some(TorsionHybridization::Sp2),
        HybridizationType::SP3 => Some(TorsionHybridization::Sp3),
        _ => None,
    }
}

fn add_torsions(
    mol: &ROMol,
    params: &AtomicParamVect,
    field: &mut ForceField,
    torsion_bond_smarts: &str,
) -> Result<()> {
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
    for m in matches {
        if m.len() != 2 {
            return Err(RDKitError::invariant(
                "torsion bond query must have two atoms",
            ));
        }
        let (idx1, idx2) = (m[0].1 as u32, m[1].1 as u32);
        let (Some(p1), Some(p2)) = (&params[idx1 as usize], &params[idx2 as usize]) else {
            continue;
        };
        let bond = mol
            .get_bond_between_atoms(idx1, idx2)?
            .ok_or_else(|| RDKitError::invariant("torsion bond not found"))?;
        let (atom1, atom2) = (mol.atom(idx1), mol.atom(idx2));
        let (Some(h1), Some(h2)) = (
            torsion_hybridization(atom1.hybridization()),
            torsion_hybridization(atom2.hybridization()),
        ) else {
            continue;
        };
        let first = field.contribs().len();
        for t_bond1 in mol.atom_bonds(idx1) {
            if t_bond1.idx() == bond.idx() {
                continue;
            }
            let b_idx = t_bond1.other_atom_idx(idx1)?;
            for t_bond2 in mol.atom_bonds(idx2) {
                if t_bond2.idx() == bond.idx() || t_bond2.idx() == t_bond1.idx() {
                    continue;
                }
                let e_idx = t_bond2.other_atom_idx(idx2)?;
                // not a three-membered ring
                if e_idx == b_idx {
                    continue;
                }
                // set when either end atom is sp2
                let has_sp2 = mol.atom(b_idx).hybridization() == HybridizationType::SP2
                    || mol.atom(e_idx).hybridization() == HybridizationType::SP2;
                field.push_contrib(Box::new(TorsionAngleContrib::new(
                    [b_idx as usize, idx1 as usize, idx2 as usize, e_idx as usize],
                    bond.bond_type_as_double()?,
                    atom1.atomic_num(),
                    atom2.atomic_num(),
                    h1,
                    h2,
                    p1,
                    p2,
                    has_sp2,
                )));
            }
        }
        // divide the force constants about this bond by their number
        let count = (field.contribs().len() - first) as u32;
        field.scale_torsions_from(first, count);
    }
    Ok(())
}

fn add_inversions(mol: &ROMol, field: &mut ForceField) -> Result<()> {
    for center in 0..mol.num_atoms() {
        let atom = mol.atom(center);
        let num = atom.atomic_num();
        // carbon, nitrogen, oxygen, phosphorous, arsenic, antimony or
        // bismuth with three neighbors
        if !matches!(num, 6 | 7 | 8 | 15 | 33 | 51 | 83) || mol.atom_degree(center) != 3 {
            continue;
        }
        // carbon, nitrogen and oxygen must be sp2
        if (6..=8).contains(&num) && atom.hybridization() != HybridizationType::SP2 {
            continue;
        }
        let nbrs: Vec<u32> = mol.atom_neighbors(center).collect();
        // set if the sp2 carbon is bound to sp2 oxygen
        let bound_to_sp2_o = num == 6
            && nbrs.iter().any(|&n| {
                let a = mol.atom(n);
                a.atomic_num() == 8 && a.hybridization() == HybridizationType::SP2
            });
        let idx = [nbrs[0], center, nbrs[1], nbrs[2]];
        for order in [[0, 1, 2, 3], [0, 1, 3, 2], [2, 1, 3, 0]] {
            field.push_contrib(Box::new(InversionContrib::new(
                [
                    idx[order[0]] as usize,
                    idx[order[1]] as usize,
                    idx[order[2]] as usize,
                    idx[order[3]] as usize,
                ],
                num,
                bound_to_sp2_o,
                1.0,
            )));
        }
    }
    Ok(())
}

/// `constructForceField(mol, params, vdwThresh, confId,
/// ignoreInterfragInteractions)`
pub fn construct_force_field_with_params(
    mol: &ROMol,
    params: &AtomicParamVect,
    vdw_thresh: f64,
    conf_id: i32,
    ignore_interfrag_interactions: bool,
) -> Result<ForceField> {
    if params.len() != mol.num_atoms() as usize {
        return Err(RDKitError::precondition("bad parameters"));
    }
    if needs_hs(mol)? {
        rd_log_warning!("Molecule does not have explicit Hs. Consider calling AddHs()");
    }
    let mut field = ForceField::new(3);
    // the atomic positions
    let conf = mol.conformer(conf_id)?;
    for i in 0..mol.num_atoms() {
        let p = conf.atom_pos(i)?;
        field.push_point(&[p.x, p.y, p.z]);
    }
    add_bonds(mol, params, &mut field)?;
    add_angles(mol, params, &mut field)?;
    // trigonal bipyramids
    for i in 0..mol.num_atoms() {
        let atom = mol.atom(i);
        if atom.hybridization() == HybridizationType::SP3D && mol.atom_degree(i) == 5 {
            add_trigonal_bipyramid_angles(mol, i, conf_id, params, &mut field)?;
        }
    }
    let neighbor_matrix = build_neighbor_matrix(mol);
    add_nonbonded(
        mol,
        conf_id,
        params,
        &mut field,
        &neighbor_matrix,
        vdw_thresh,
        ignore_interfrag_interactions,
    )?;
    add_torsions(mol, params, &mut field, DEFAULT_TORSION_BOND_SMARTS)?;
    add_inversions(mol, &mut field)?;
    Ok(field)
}

/// `constructForceField(mol, vdwThresh, confId, ignoreInterfragInteractions)`
pub fn construct_force_field(
    mol: &ROMol,
    vdw_thresh: f64,
    conf_id: i32,
    ignore_interfrag_interactions: bool,
) -> Result<ForceField> {
    let (params, _) = get_atom_types(mol)?;
    construct_force_field_with_params(
        mol,
        &params,
        vdw_thresh,
        conf_id,
        ignore_interfrag_interactions,
    )
}
