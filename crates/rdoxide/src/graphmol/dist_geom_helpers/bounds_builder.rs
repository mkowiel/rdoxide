//! Distance bounds from the molecular topology (upstream
//! `BoundsMatrixBuilder.cpp`).
//!
//! The 1-2 bounds come from force-field bond lengths, the 1-3 bounds from bond
//! angles, the 1-4 bounds from the allowed torsion angles, and optionally 1-5
//! bounds for paths with rigid 1-4 fragments. The remaining lower bounds are
//! van der Waals contacts.

use std::collections::{BTreeMap, HashSet};
use std::f64::consts::PI;

use super::details::{
    BitSet, Bounds, Path14Configuration, Path14Vect, SymmGrid, TorsionType, TorsionValue, Type14,
    merge_bounds, unified_id2, unified_id3,
};
use super::embed_params::{EmbedFF, EmbedParameters, InitialEmbeddingMode};
use super::internal_coordinates::InternalCoordinates;
use crate::distgeom::BoundsMatrix;
use crate::error::{RDKitError, Result};
use crate::forcefield::uff::calc_bond_rest_length;
use crate::geometry::utils::{
    compute_13_dist, compute_14_dist_3d, compute_14_dist_cis, compute_14_dist_trans,
};
use crate::graphmol::atom::{Atom, HybridizationType};
use crate::graphmol::bond::{Bond, BondStereo, BondType};
use crate::graphmol::chirality::nontetrahedral::{
    get_ideal_angle_between_ligands, has_non_tetrahedral_stereo,
};
use crate::graphmol::forcefield_helpers::mmff::MmffMolProperties;
use crate::graphmol::forcefield_helpers::uff::{AtomicParamVect, get_atom_types};
use crate::graphmol::matrices::get_distance_mat;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;
use crate::rd_log_warning;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;

const DIST12_DELTA: f64 = 0.01;
const ANGLE_DELTA: f64 = 0.035;
const GEN_DIST_TOL: f64 = 0.06;
const DIST15_TOL: f64 = 0.08;
const VDW_SCALE_15: f64 = 0.7;
const H_BOND_LENGTH: f64 = 1.8;
const MAX_UPPER: f64 = 1000.0;
const MIN_MACROCYCLE_RING_SIZE: usize = 9;
const FAILED12: f64 = -1.0;

/// How far a 1-4 path is known to be cis or trans, and what is needed to
/// classify it.
#[derive(Debug, Clone, Copy, Default)]
struct Optional14Info {
    force_trans_amides: bool,
    use_macrocycle_14_config: bool,
    ring_size: usize,
    prefer_trans: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum DistType {
    Dist12,
    Dist13,
    Dist14,
}

/// What has been computed about the molecule so far.
struct ComputedData {
    bond_lengths: Vec<f64>,
    /// the central atom of a pair of adjacent bonds (-1: not adjacent)
    bond_adj: SymmGrid<i32>,
    bond_angles: SymmGrid<f64>,
    paths14: Path14Vect,
    cis_paths: HashSet<usize>,
    trans_paths: HashSet<usize>,
    set15_atoms: BitSet,
    visited12: BitSet,
    visited13: BitSet,
    visited14: BitSet,
}

impl ComputedData {
    fn new(num_atoms: usize, num_bonds: usize) -> Self {
        ComputedData {
            bond_lengths: vec![0.0; num_bonds],
            bond_adj: SymmGrid::new(num_bonds, -1),
            bond_angles: SymmGrid::new(num_bonds, -1.0),
            paths14: Vec::new(),
            cis_paths: HashSet::new(),
            trans_paths: HashSet::new(),
            set15_atoms: BitSet::new(num_atoms * num_atoms),
            visited12: BitSet::new(num_atoms * num_atoms),
            visited13: BitSet::new(num_atoms * num_atoms),
            visited14: BitSet::new(num_atoms * num_atoms),
        }
    }

    /// Whether the pair was handled as a distance of at most `max_dist_type`.
    fn visited_bound(&self, pid: usize, max_dist_type: DistType) -> bool {
        (max_dist_type >= DistType::Dist12 && self.visited12.get(pid))
            || (max_dist_type >= DistType::Dist13 && self.visited13.get(pid))
            || (max_dist_type >= DistType::Dist14 && self.visited14.get(pid))
    }
}

/// `_checkAndSetBounds`: tightens the bounds of (i, j) pessimistically, or
/// with `set_if_better` replaces them by the overlap of old and new bounds.
pub fn check_and_set_bounds(
    i: u32,
    j: u32,
    lb: f64,
    ub: f64,
    mmat: &mut BoundsMatrix,
    set_if_better: bool,
) -> Result<()> {
    let (i, j) = (i as usize, j as usize);
    let clb = mmat.lower_bound(i, j);
    let cub = mmat.upper_bound(i, j);

    if ub <= lb {
        return Err(RDKitError::invariant(
            "upper bound not greater than lower bound",
        ));
    }
    if !(lb > DIST12_DELTA || clb > DIST12_DELTA) {
        return Err(RDKitError::invariant("bad lower bound"));
    }

    if set_if_better {
        let mut nlb = clb.max(lb);
        let mut nub = cub.min(ub);
        if nub <= nlb {
            // the ranges do not overlap: be conservative
            nlb = clb.min(lb);
            nub = cub.max(ub);
        }
        mmat.set_lower_bound(i, j, nlb)?;
        mmat.set_upper_bound(i, j, nub)?;
    } else {
        if clb <= DIST12_DELTA {
            mmat.set_lower_bound(i, j, lb)?;
        } else if lb < clb && lb > DIST12_DELTA {
            mmat.set_lower_bound(i, j, lb)?;
        }
        if cub >= MAX_UPPER {
            mmat.set_upper_bound(i, j, ub)?;
        } else if ub > cub && ub < MAX_UPPER {
            mmat.set_upper_bound(i, j, ub)?;
        }
    }
    Ok(())
}

/// Larger heteroatoms in conjugated 5-rings need a bit of extra flex.
fn squish_bond(mol: &ROMol, bond: &Bond) -> bool {
    bond.is_conjugated()
        && (mol.atom(bond.begin_atom_idx()).atomic_num() > 10
            || mol.atom(bond.end_atom_idx()).atomic_num() > 10)
        && mol
            .ring_info()
            .is_bond_in_ring_of_size(bond.idx(), 5)
            .unwrap_or(false)
}

/// Atom parameters of the force field that supplies the 1-2 and 1-3 values.
enum Parametrization {
    Uff(AtomicParamVect),
    Mmff(Box<MmffMolProperties>),
}

fn parametrize_uff(mol: &ROMol) -> Result<(Parametrization, bool)> {
    let (params, found_all) = get_atom_types(mol)?;
    Ok((Parametrization::Uff(params), found_all))
}

/// Works on a copy of the molecule: typing kekulizes it, and setting bounds
/// must not change the molecule that is embedded.
fn parametrize_mmff(mol: &ROMol) -> Result<(Parametrization, bool)> {
    let mut copy = mol.clone();
    let props = MmffMolProperties::new(&mut copy, "MMFF94")?;
    let valid = props.is_valid();
    Ok((Parametrization::Mmff(Box::new(props)), valid))
}

/// The ideal bond length of the bond i-j, or `FAILED12`.
fn calc12(mol: &ROMol, bond: &Bond, params: &Parametrization, i: u32, j: u32) -> Result<f64> {
    match params {
        Parametrization::Uff(atom_params) => {
            let order = bond.bond_type_as_double()?;
            match (atom_params[i as usize], atom_params[j as usize]) {
                (Some(pi), Some(pj)) if order > 0.0 => Ok(calc_bond_rest_length(order, &pi, &pj)),
                _ => Ok(FAILED12),
            }
        }
        Parametrization::Mmff(props) => Ok(props
            .bond_stretch_params(mol, i, j)?
            .map_or(FAILED12, |p| p.params.r0)),
    }
}

const FORCE_FIELD_FALLBACK_WARNING: &str = "Molecule was not parametrized succesfully with the selected force field. Falling back to the Universal Force Field!";

/// Sets the 1-2 bounds from bond lengths of the force field (mostly), widened
/// by a small tolerance. Returns whether all bonds were parametrized
/// by the force field, and its parameters.
fn set12_bounds(
    mol: &ROMol,
    mmat: &mut BoundsMatrix,
    accum: &mut ComputedData,
    parametrize: fn(&ROMol) -> Result<(Parametrization, bool)>,
    is_fallback: bool,
) -> Result<(bool, Parametrization)> {
    let na = mol.num_atoms();
    if mmat.num_rows() != na as usize {
        return Err(RDKitError::invariant("Wrong size metric matrix"));
    }
    if accum.bond_lengths.len() < mol.num_bonds() as usize {
        return Err(RDKitError::invariant("Wrong size accumData"));
    }
    let (atom_params, found_all) = parametrize(mol)?;
    if !is_fallback && !found_all {
        rd_log_warning!("{FORCE_FIELD_FALLBACK_WARNING}");
        return Ok((found_all, atom_params));
    }
    let mut squish_atoms = vec![false; na as usize];
    if mol.ring_info().is_initialized() {
        for bond in mol.bonds() {
            if squish_bond(mol, bond) {
                squish_atoms[bond.begin_atom_idx() as usize] = true;
                squish_atoms[bond.end_atom_idx() as usize] = true;
            }
        }
    }

    let mut success = true;
    for bond in mol.bonds() {
        let beg = bond.begin_atom_idx();
        let end = bond.end_atom_idx();
        let bl = calc12(mol, bond, &atom_params, beg, end)?;
        if bl > FAILED12 {
            let extra_squish = if squish_atoms[beg as usize] && squish_atoms[end as usize] {
                0.2 // empirical
            } else {
                0.0
            };
            accum.bond_lengths[bond.idx() as usize] = bl;
            mmat.set_upper_bound(beg as usize, end as usize, bl + extra_squish + DIST12_DELTA)?;
            mmat.set_lower_bound(beg as usize, end as usize, bl - extra_squish - DIST12_DELTA)?;
        } else if is_fallback {
            // no force-field value: estimate from the periodic table
            let pt = PeriodicTable::get();
            let z1 = mol.atom(beg).atomic_num();
            let z2 = mol.atom(end).atomic_num();
            let mut bl = pt.get_rcovalent(z1)? + pt.get_rcovalent(z2)?;
            let (upper_scale, lower_scale);
            let bt = bond.bond_type();
            if bt > BondType::Aromatic || bt < BondType::Single {
                // weird bond types: use the mean van der Waals radius and
                // allow a lot more flex
                bl = (pt.get_rvdw(z1)? + pt.get_rvdw(z2)?) / 2.0;
                upper_scale = 1.5;
                lower_scale = 0.75;
            } else {
                // Pauling's formula for the bond order (from the UFF code)
                const PAULING_LAMBDA: f64 = 0.1332;
                bl -= PAULING_LAMBDA * bond.bond_type_as_double()?.ln() * bl;
                upper_scale = 1.1;
                lower_scale = 0.9;
            }
            accum.bond_lengths[bond.idx() as usize] = bl;
            mmat.set_upper_bound(beg as usize, end as usize, upper_scale * bl)?;
            mmat.set_lower_bound(beg as usize, end as usize, lower_scale * bl)?;
            success = false;
        } else {
            // the force field failed for this bond: redo with the fallback
            rd_log_warning!("{FORCE_FIELD_FALLBACK_WARNING}");
            return Ok((false, atom_params));
        }
        accum.visited12.set(unified_id2(beg, end, na));
    }
    Ok((success, atom_params))
}

fn is_h_bond_acceptor(atom: &Atom) -> bool {
    atom.atomic_num() == 7 || atom.atomic_num() == 8
}

fn is_h_in_h_bond_donor(mol: &ROMol, idx: u32) -> bool {
    mol.atom(idx).atomic_num() == 1
        && mol
            .atom_neighbors(idx)
            .any(|n| matches!(mol.atom(n).atomic_num(), 7 | 8))
}

/// Lower bounds from van der Waals radii for pairs without other bounds.
///
/// Pairs 4 or 5 bonds apart get scaled radii so the atoms can get closer; a
/// hydrogen of a donor and an acceptor get the length of a hydrogen bond.
/// (Like upstream, the flags of atom 0 are never set.)
fn set_lower_bound_vdw(mol: &ROMol, mmat: &mut BoundsMatrix, dmat: &[f64]) -> Result<()> {
    let npt = mmat.num_rows();
    if npt != mol.num_atoms() as usize {
        return Err(RDKitError::precondition("Wrong size metric matrix"));
    }
    let pt = PeriodicTable::get();
    let mut h_in_h_bond_donors = vec![false; npt];
    let mut h_bond_acceptors = vec![false; npt];
    for i in 1..npt {
        let atom_i = mol.atom(i as u32);
        let vw1 = pt.get_rvdw(atom_i.atomic_num())?;
        if is_h_in_h_bond_donor(mol, i as u32) {
            h_in_h_bond_donors[i] = true;
        }
        if is_h_bond_acceptor(atom_i) {
            h_bond_acceptors[i] = true;
        }
        for j in 0..i {
            let vw2 = pt.get_rvdw(mol.atom(j as u32).atomic_num())?;
            if mmat.lower_bound(i, j) < DIST12_DELTA {
                let d = dmat[i * npt + j];
                let lower = if (h_in_h_bond_donors[i] && h_bond_acceptors[j])
                    || (h_bond_acceptors[i] && h_in_h_bond_donors[j])
                {
                    H_BOND_LENGTH
                } else if d == 4.0 {
                    VDW_SCALE_15 * (vw1 + vw2)
                } else if d == 5.0 {
                    (VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * (vw1 + vw2)
                } else {
                    vw1 + vw2
                };
                mmat.set_lower_bound(i, j, lower)?;
            }
        }
    }
    Ok(())
}

fn is_larger_sp2_atom(mol: &ROMol, idx: u32) -> bool {
    let atom = mol.atom(idx);
    atom.atomic_num() > 13
        && atom.hybridization() == HybridizationType::SP2
        && mol.ring_info().num_atom_rings(idx).unwrap_or(0) > 0
}

/// The angle tolerance is doubled for each larger atom outside the first row
/// of the periodic table.
fn scale_tolerance_for_large_sp2(mol: &ROMol, a1: u32, a2: u32, a3: u32) -> f64 {
    let mut tolerance = ANGLE_DELTA;
    for idx in [a1, a2, a3] {
        if is_larger_sp2_atom(mol, idx) {
            tolerance *= 2.0;
        }
    }
    tolerance
}

fn get_angle_range(angle: f64, tolerance: f64) -> (f64, f64) {
    if angle + tolerance >= PI {
        (PI - tolerance, PI)
    } else if angle - tolerance <= 0.0 {
        (0.0, tolerance)
    } else {
        (angle - tolerance, angle + tolerance)
    }
}

fn set13_bounds_helper(
    aid1: u32,
    aid: u32,
    aid3: u32,
    angle: f64,
    mmat: &mut BoundsMatrix,
    mol: &ROMol,
) -> Result<()> {
    let tolerance = scale_tolerance_for_large_sp2(mol, aid1, aid, aid3);
    let (angle_l, angle_u) = get_angle_range(angle, tolerance);
    let (a1, a, a3) = (aid1 as usize, aid as usize, aid3 as usize);
    let du = compute_13_dist(mmat.upper_bound(a1, a), mmat.upper_bound(a, a3), angle_u);
    let dl = compute_13_dist(mmat.lower_bound(a1, a), mmat.lower_bound(a, a3), angle_l);
    check_and_set_bounds(aid1, aid3, dl, du, mmat, false)
}

/// The angle at an atom in a ring of the given size.
///
/// This assumes that all angles in a ring are equal, which is not always true
/// (e.g. aromatic rings with heteroatoms like s1cncc1).
fn ring_angle(atom: &Atom, ring_size: usize) -> f64 {
    let hyb = atom.hybridization();
    if (hyb == HybridizationType::SP2 && ring_size <= 8) || ring_size == 3 || ring_size == 4 {
        PI * (1.0 - 2.0 / ring_size as f64)
    } else if hyb == HybridizationType::SP3 {
        if ring_size == 5 {
            104.0 * PI / 180.0
        } else {
            109.5 * PI / 180.0
        }
    } else if hyb == HybridizationType::SP3D {
        105.0 * PI / 180.0
    } else if hyb == HybridizationType::SP3D2 {
        90.0 * PI / 180.0
    } else {
        120.0 * PI / 180.0
    }
}

fn check_13_sizes(mol: &ROMol, mmat: &BoundsMatrix, accum: &ComputedData) -> Result<()> {
    if mmat.num_rows() != mol.num_atoms() as usize {
        return Err(RDKitError::invariant("Wrong size metric matrix"));
    }
    if accum.bond_angles.num_rows() != mol.num_bonds() as usize {
        return Err(RDKitError::invariant("Wrong size bond angle matrix"));
    }
    if accum.bond_adj.num_rows() != mol.num_bonds() as usize {
        return Err(RDKitError::invariant("Wrong size bond adjacency matrix"));
    }
    Ok(())
}

/// Sets the 1-3 bounds from bond angles and lengths.
///
/// Most special cases come from rings, so the atoms of all simple rings are
/// handled first (smallest rings first), remembering how much of the angle
/// around each atom has been used. Then the remaining pairs of bonds at a ring
/// atom (non-ring neighbours, or atoms of different fused rings) get the rest
/// of the angle, and finally the pairs at non-ring atoms are handled by
/// hybridization.
fn set13_bounds(mol: &ROMol, mmat: &mut BoundsMatrix, accum: &mut ComputedData) -> Result<()> {
    check_13_sizes(mol, mmat, accum)?;
    let npt = mmat.num_rows();
    let na = mol.num_atoms();
    let nb = mol.num_bonds();
    let rinfo = mol.ring_info();

    let mut atom_rings: Vec<Vec<i32>> = rinfo.atom_rings().to_vec();
    libcxx_sort_by(&mut atom_rings, |a, b| a.len() < b.len());

    let mut visited = vec![0u32; npt];
    let mut angle_taken = vec![0.0f64; npt];
    let mut done_paths = BitSet::new(nb as usize * nb as usize);

    // first the atoms in rings
    for ring in &atom_rings {
        let r_size = ring.len();
        let mut aid1 = ring[r_size - 1] as u32;
        for i in 0..r_size {
            let aid2 = ring[i] as u32;
            let aid3 = (if i == r_size - 1 {
                ring[0]
            } else {
                ring[i + 1]
            }) as u32;
            let bid1 = mol
                .bond_idx_between(aid1, aid2)
                .ok_or_else(|| RDKitError::invariant("no bond found"))?;
            let bid2 = mol
                .bond_idx_between(aid2, aid3)
                .ok_or_else(|| RDKitError::invariant("no bond found"))?;
            let bond_pair_id = unified_id2(bid1, bid2, nb);
            if !done_paths.get(bond_pair_id) {
                // bridged systems cover the same ring paths several times,
                // which must not add to the angle taken at the central atom
                // (Issue 215)
                let angle = ring_angle(mol.atom(aid2), r_size);
                let pid = unified_id2(aid1, aid3, na);
                if !accum.visited_bound(pid, DistType::Dist12) {
                    set13_bounds_helper(aid1, aid2, aid3, angle, mmat, mol)?;
                    accum.visited13.set(pid);
                }
                accum.bond_angles.set(bid1, bid2, angle);
                accum.bond_adj.set(bid1, bid2, aid2 as i32);
                visited[aid2 as usize] += 1;
                angle_taken[aid2 as usize] += angle;
                done_paths.set(bond_pair_id);
            }
            aid1 = aid2;
        }
    }

    // then the remaining pairs of bonds
    for aid2 in 0..npt as u32 {
        let atom = mol.atom(aid2);
        let deg = mol.atom_degree(aid2);
        let n13 = (deg * deg.wrapping_sub(1)) / 2;
        if n13 == visited[aid2 as usize] {
            continue;
        }
        let ahyb = atom.hybridization();
        let bonds: Vec<u32> = mol.atom_bond_indices(aid2).collect();
        if visited[aid2 as usize] >= 1 {
            // atoms already seen in rings: set the 1-3 distances for
            //  1) non-ring atoms with a ring atom in between,
            //  2) a non-ring atom and a ring atom with a ring atom in between,
            //  3) ring atoms of different rings of a fused system
            for (i1, &bid1) in bonds.iter().enumerate() {
                let aid1 = mol.bond(bid1).other_atom_idx(aid2)?;
                for &bid2 in &bonds[..i1] {
                    let aid3 = mol.bond(bid2).other_atom_idx(aid2)?;
                    if accum.bond_angles.get(bid1, bid2) >= 0.0 {
                        continue;
                    }
                    let angle = if ahyb == HybridizationType::SP2 {
                        // planar: divide the remaining angle among the
                        // remaining 1-3 configurations (there should be one)
                        (2.0 * PI - angle_taken[aid2 as usize])
                            / f64::from(n13.wrapping_sub(visited[aid2 as usize]))
                    } else if ahyb == HybridizationType::SP3 {
                        // the tetrahedral angle, but 3- and 4-membered rings
                        // are not close to it
                        if rinfo.is_atom_in_ring_of_size(aid2, 3)? {
                            116.0 * PI / 180.0
                        } else if rinfo.is_atom_in_ring_of_size(aid2, 4)? {
                            112.0 * PI / 180.0
                        } else {
                            109.5 * PI / 180.0
                        }
                    } else if has_non_tetrahedral_stereo(mol, aid2) {
                        get_ideal_angle_between_ligands(mol, aid2, aid1, aid3) * PI / 180.0
                    } else if deg == 5 {
                        105.0 * PI / 180.0
                    } else if deg == 6 {
                        135.0 * PI / 180.0
                    } else {
                        120.0 * PI / 180.0 // FIX: probably not the best default
                    };

                    let pid = unified_id2(aid1, aid3, na);
                    if !accum.visited_bound(pid, DistType::Dist12) {
                        set13_bounds_helper(aid1, aid2, aid3, angle, mmat, mol)?;
                        accum.visited13.set(pid);
                    }
                    accum.bond_angles.set(bid1, bid2, angle);
                    accum.bond_adj.set(bid1, bid2, aid2 as i32);
                    angle_taken[aid2 as usize] += angle;
                    visited[aid2 as usize] += 1;
                }
            }
        } else if visited[aid2 as usize] == 0 {
            // non-ring atoms: angles from the hybridization
            let non_tetrahedral = has_non_tetrahedral_stereo(mol, aid2);
            for (i1, &bid1) in bonds.iter().enumerate() {
                let aid1 = mol.bond(bid1).other_atom_idx(aid2)?;
                for &bid2 in &bonds[..i1] {
                    let aid3 = mol.bond(bid2).other_atom_idx(aid2)?;
                    let angle = if non_tetrahedral {
                        get_ideal_angle_between_ligands(mol, aid2, aid1, aid3) * PI / 180.0
                    } else if ahyb == HybridizationType::SP {
                        PI
                    } else if ahyb == HybridizationType::SP2 {
                        2.0 * PI / 3.0
                    } else if ahyb == HybridizationType::SP3 {
                        109.5 * PI / 180.0
                    } else if ahyb == HybridizationType::SP3D {
                        // FIX: this and the next default should probably be
                        // special cased, particularly when the stereochemistry
                        // is known
                        105.0 * PI / 180.0
                    } else if ahyb == HybridizationType::SP3D2 {
                        135.0 * PI / 180.0
                    } else {
                        120.0 * PI / 180.0
                    };
                    let pid = unified_id2(aid1, aid3, na);
                    if !accum.visited_bound(pid, DistType::Dist12) {
                        if deg <= 4
                            || (non_tetrahedral && atom.props.has_val(cp::_chiralPermutation))
                        {
                            set13_bounds_helper(aid1, aid2, aid3, angle, mmat, mol)?;
                        } else {
                            // just use 180 as the max angle and an arbitrary
                            // min angle
                            let dmax = accum.bond_lengths[bid1 as usize]
                                + accum.bond_lengths[bid2 as usize];
                            check_and_set_bounds(aid1, aid3, 1.0, dmax * 1.2, mmat, false)?;
                        }
                        accum.visited13.set(pid);
                    }
                    accum.bond_angles.set(bid1, bid2, angle);
                    accum.bond_adj.set(bid1, bid2, aid2 as i32);
                    angle_taken[aid2 as usize] += angle;
                    visited[aid2 as usize] += 1;
                }
            }
        }
    }
    Ok(())
}

/// The atoms (aid1, aid2, aid3) of two bonds sharing an atom, with aid2 the
/// shared one; `None` if the bonds are not adjacent.
fn adjacent_bond_atoms(bondi: &Bond, bondj: &Bond) -> Option<(u32, u32, u32)> {
    let (a11, a12) = (bondi.begin_atom_idx(), bondi.end_atom_idx());
    let (a21, a22) = (bondj.begin_atom_idx(), bondj.end_atom_idx());
    if a12 == a21 {
        Some((a11, a12, a22))
    } else if a12 == a22 {
        Some((a11, a12, a21))
    } else if a11 == a21 {
        Some((a12, a11, a22))
    } else if a11 == a22 {
        Some((a12, a11, a21))
    } else {
        None
    }
}

/// Sets the 1-3 bounds from the MMFF equilibrium angles.
fn set13_bounds_mmff(
    mol: &ROMol,
    mmat: &mut BoundsMatrix,
    accum: &mut ComputedData,
    mmff: &MmffMolProperties,
) -> Result<()> {
    check_13_sizes(mol, mmat, accum)?;
    let na = mol.num_atoms();
    for bondi in mol.bonds() {
        let bid1 = bondi.idx();
        for bondj in &mol.bonds()[bid1 as usize + 1..] {
            let bid2 = bondj.idx();
            let Some((aid1, aid2, aid3)) = adjacent_bond_atoms(bondi, bondj) else {
                continue;
            };
            let params = mmff.angle_bend_params(mol, aid1, aid2, aid3)?;
            if params.is_none() {
                rd_log_warning!(
                    "Bounds matrix builder: Invalid MMFF angle parameters for ({aid1}, {aid2}, {aid3})"
                );
            }
            // theta0 is in degrees
            let angle = params.map_or(0.0, |p| p.params.theta0) * PI / 180.0;
            let pid = unified_id2(aid1, aid3, na);
            if !accum.visited_bound(pid, DistType::Dist12) {
                set13_bounds_helper(aid1, aid2, aid3, angle, mmat, mol)?;
                accum.visited13.set(pid);
            }
            accum.bond_angles.set(bid1, bid2, angle);
            accum.bond_adj.set(bid1, bid2, aid2 as i32);
        }
    }
    Ok(())
}

/// The stereo of a double bond seen from the atoms aid1 and aid4.
fn get_atom_stereo(bond: &Bond, aid1: u32, aid4: u32) -> BondStereo {
    let mut stype = bond.stereo();
    let st_atoms = bond.stereo_atoms();
    if stype > BondStereo::StereoAny && st_atoms.len() >= 2 {
        if (st_atoms[0] as u32 != aid1) ^ (st_atoms[1] as u32 != aid4) {
            stype = match stype {
                BondStereo::StereoZ => BondStereo::StereoE,
                BondStereo::StereoE => BondStereo::StereoZ,
                BondStereo::StereoCis => BondStereo::StereoTrans,
                BondStereo::StereoTrans => BondStereo::StereoCis,
                other => other,
            };
        }
    }
    stype
}

fn is_cis(stype: BondStereo) -> bool {
    matches!(stype, BondStereo::StereoZ | BondStereo::StereoCis)
}

fn is_trans(stype: BondStereo) -> bool {
    matches!(stype, BondStereo::StereoE | BondStereo::StereoTrans)
}

fn both_sp2(mol: &ROMol, a2: u32, a3: u32) -> bool {
    mol.atom(a2).hybridization() == HybridizationType::SP2
        && mol.atom(a3).hybridization() == HybridizationType::SP2
}

/// The atoms of a 1-4 path.
#[derive(Clone, Copy)]
struct Atoms14 {
    a1: u32,
    a2: u32,
    a3: u32,
    a4: u32,
}

/// The bonds of a 1-4 path.
#[derive(Clone, Copy)]
struct Bonds14 {
    b1: u32,
    b2: u32,
    b3: u32,
}

fn get_in_ring_14_type(mol: &ROMol, bnd2: &Bond, at: Atoms14, ring_size: usize) -> TorsionValue {
    let stype = get_atom_stereo(bnd2, at.a1, at.a4);
    let sp2 = both_sp2(mol, at.a2, at.a3);
    // there is no reason to assume cis bonds in bigger rings (github #1240:
    // failure to embed larger aromatic rings)
    if (ring_size <= 8 && sp2 && !is_trans(stype)) || is_cis(stype) {
        TorsionValue::of(TorsionType::Cis)
    } else if is_trans(stype) {
        TorsionValue::of(TorsionType::Trans)
    } else if sp2 {
        TorsionValue::of(TorsionType::CisTrans)
    } else {
        TorsionValue::of(TorsionType::Flexible)
    }
}

fn get_two_in_same_ring_14_type(
    mol: &ROMol,
    bnd2: &Bond,
    at: Atoms14,
    prefer_trans: bool,
) -> TorsionValue {
    let stype = get_atom_stereo(bnd2, at.a1, at.a4);
    let sp2 = both_sp2(mol, at.a2, at.a3);
    if (prefer_trans && sp2 && !is_cis(stype)) || is_trans(stype) {
        // basically a flat ring with an external substituent
        TorsionValue::of(TorsionType::Trans)
    } else if is_cis(stype) {
        TorsionValue::of(TorsionType::Cis)
    } else if sp2 {
        TorsionValue::of(TorsionType::CisTrans)
    } else {
        TorsionValue::of(TorsionType::Flexible)
    }
}

/// An amide or ester seen from the atoms of the path
///
/// ```text
///        4    <- 4 is the O
///        |    <- the double bond
///    1   3
///     \ / \
///      2   5  <- 2 is an oxygen/nitrogen
/// ```
fn check_amide_ester_14(
    mol: &ROMol,
    bnd1: &Bond,
    bnd3: &Bond,
    a2: u32,
    a3: u32,
    a4: u32,
) -> Result<bool> {
    let a2_num = mol.atom(a2).atomic_num();
    let a3_num = mol.atom(a3).atomic_num();
    let a4_num = mol.atom(a4).atomic_num();
    Ok(a3_num == 6
        && bnd3.bond_type() == BondType::Double
        && (a4_num == 8 || a4_num == 7)
        && bnd1.bond_type() == BondType::Single
        && (a2_num == 8 || (a2_num == 7 && mol.total_num_hs(a2, true)? == 1)))
}

/// Like [`check_amide_ester_14`] when all three bonds are part of the
/// macrocycle ring. (This is a rewrite with more explicit logic; it samples
/// macrocycles better.)
fn check_macrocycle_all_in_same_ring_amide_ester_14(mol: &ROMol, at: Atoms14) -> bool {
    let a2_num = mol.atom(at.a2).atomic_num();
    if mol.atom(at.a3).atomic_num() != 6 {
        return false;
    }
    if !(a2_num == 7 || a2_num == 8)
        || !(mol.atom_degree(at.a2) == 3 && mol.atom_degree(at.a3) == 3)
    {
        return false;
    }
    if let Some(nbr) = mol
        .atom_neighbors(at.a2)
        .find(|&n| n != at.a1 && n != at.a3)
    {
        // a (methylated) amide
        let res = mol.atom(nbr);
        let resbnd = mol.bond(mol.bond_idx_between(at.a2, nbr).expect("neighbor bond"));
        if (res.atomic_num() != 6 && res.atomic_num() != 1)
            || resbnd.bond_type() != BondType::Single
        {
            return false;
        }
    }
    if let Some(nbr) = mol
        .atom_neighbors(at.a3)
        .find(|&n| n != at.a2 && n != at.a4)
    {
        // the carbonyl oxygen
        let res = mol.atom(nbr);
        let resbnd = mol.bond(mol.bond_idx_between(at.a3, nbr).expect("neighbor bond"));
        if res.atomic_num() != 8 || resbnd.bond_type() != BondType::Double {
            return false;
        }
    }
    true
}

fn is_carbonyl(mol: &ROMol, idx: u32) -> bool {
    let atom = mol.atom(idx);
    if atom.atomic_num() == 6 && mol.atom_degree(idx) > 2 {
        for nbr in mol.atom_neighbors(idx) {
            let z = mol.atom(nbr).atomic_num();
            if (z == 8 || z == 7)
                && mol
                    .bond_idx_between(idx, nbr)
                    .is_some_and(|b| mol.bond(b).bond_type() == BondType::Double)
            {
                return true;
            }
        }
    }
    false
}

fn check_amide_ester_15(mol: &ROMol, bnd1: &Bond, bnd3: &Bond, a2: u32, a3: u32) -> Result<bool> {
    let a2_num = mol.atom(a2).atomic_num();
    Ok(
        (a2_num == 8 || (a2_num == 7 && mol.total_num_hs(a2, true)? == 1))
            && bnd1.bond_type() == BondType::Single
            && mol.atom(a3).atomic_num() == 6
            && bnd3.bond_type() == BondType::Single
            && is_carbonyl(mol, a3),
    )
}

/// A secondary amide's hydrogen at the end of the path
fn is_secondary_amide_h(mol: &ROMol, at: Atoms14) -> Result<bool> {
    let side = |h: u32, n: u32| -> Result<bool> {
        Ok(mol.atom(h).atomic_num() == 1
            && mol.atom(n).atomic_num() == 7
            && mol.atom_degree(n) == 3
            && mol.total_num_hs(n, true)? == 1)
    };
    Ok(side(at.a1, at.a2)? || side(at.a4, at.a3)?)
}

fn get_chain_14_type(
    mol: &ROMol,
    bonds: [&Bond; 3],
    at: Atoms14,
    force_trans_amides: bool,
) -> Result<TorsionValue> {
    let [bnd1, bnd2, bnd3] = bonds;
    match bnd2.bond_type() {
        BondType::Double => {
            if bnd1.bond_type() == BondType::Double || bnd3.bond_type() == BondType::Double {
                // CC=C=C: the torsion angle is zero
                Ok(TorsionValue::of(TorsionType::Cis))
            } else if bnd2.stereo() > BondStereo::StereoAny {
                if is_cis(get_atom_stereo(bnd2, at.a1, at.a4)) {
                    Ok(TorsionValue::of(TorsionType::Cis))
                } else {
                    Ok(TorsionValue::of(TorsionType::Trans))
                }
            } else if both_sp2(mol, at.a2, at.a3) {
                Ok(TorsionValue::of(TorsionType::CisTrans))
            } else {
                Ok(TorsionValue::of(TorsionType::Flexible))
            }
        }
        BondType::Single => {
            if mol.atom(at.a2).atomic_num() == 16
                && mol.atom(at.a3).atomic_num() == 16
                && mol.atom_degree(at.a2) == 2
                && mol.atom_degree(at.a3) == 2
            {
                // *S-S*
                return Ok(TorsionValue {
                    value: Some(PI / 2.0),
                    ..TorsionValue::of(TorsionType::Custom)
                });
            }
            if check_amide_ester_14(mol, bnd1, bnd3, at.a2, at.a3, at.a4)?
                || check_amide_ester_14(mol, bnd3, bnd1, at.a3, at.a2, at.a1)?
            {
                // Here we set the distance between atoms 1 and 4, 1 and 5 is
                // handled below. Issue284: forcing amides to be trans is
                // convenient but bogus; the torsion may roam from cis to
                // trans and the force field planarizes things later.
                return Ok(if !force_trans_amides {
                    TorsionValue::of(TorsionType::CisTrans)
                } else if is_secondary_amide_h(mol, at)? {
                    // the H should be trans to the O
                    TorsionValue::forced(TorsionType::Trans)
                } else {
                    TorsionValue::forced(TorsionType::Cis)
                });
            }
            if check_amide_ester_15(mol, bnd1, bnd3, at.a2, at.a3)?
                || check_amide_ester_15(mol, bnd3, bnd1, at.a3, at.a2)?
            {
                // the 1-4 contact was set above, here is 1-5 (Issue355)
                return Ok(if !force_trans_amides {
                    TorsionValue::of(TorsionType::CisTrans)
                } else if is_secondary_amide_h(mol, at)? {
                    // the H is cis to atom 5
                    TorsionValue::forced(TorsionType::Cis)
                } else {
                    TorsionValue::forced(TorsionType::Trans)
                });
            }
            Ok(TorsionValue::of(TorsionType::Flexible))
        }
        _ => Ok(TorsionValue::of(TorsionType::Flexible)),
    }
}

fn record_14_path(
    mol: &ROMol,
    bonds: Bonds14,
    accum: &mut ComputedData,
    r_size: usize,
) -> Result<()> {
    let Bonds14 { b1, b2, b3 } = bonds;
    let nb = mol.num_bonds();
    let aid2 = accum.bond_adj.get(b1, b2) as u32;
    let aid3 = accum.bond_adj.get(b2, b3) as u32;
    let aid1 = mol.bond(b1).other_atom_idx(aid2)?;
    let aid4 = mol.bond(b2).other_atom_idx(aid3)?;

    let value = if both_sp2(mol, aid2, aid3) {
        // FIX: check for trans
        accum.cis_paths.insert(unified_id3(b1, b2, b3, nb));
        TorsionValue::of(TorsionType::Cis)
    } else {
        TorsionValue::of(TorsionType::Flexible)
    };
    accum.paths14.push(Path14Configuration {
        bid1: b1,
        bid2: b2,
        bid3: b3,
        aid1,
        aid2,
        aid3,
        aid4,
        value,
        type14: Type14::InRing,
        r_size,
    });
    Ok(())
}

/// Adapted from [`check_amide_ester_14`] without the hydrogen count of atom 2
/// (the count misses hydrogens that are present), so that methylated amides
/// are recognised; only used for macrocycles.
fn check_macrocycle_two_in_same_ring_amide_ester_14(
    mol: &ROMol,
    bnd1: &Bond,
    bnd3: &Bond,
    a1: u32,
    a2: u32,
    a3: u32,
    a4: u32,
) -> bool {
    let a2_num = mol.atom(a2).atomic_num();
    let a4_num = mol.atom(a4).atomic_num();
    mol.atom(a1).atomic_num() != 1
        && mol.atom(a3).atomic_num() == 6
        && bnd3.bond_type() == BondType::Double
        && (a4_num == 8 || a4_num == 7)
        && bnd1.bond_type() == BondType::Single
        && (a2_num == 8 || a2_num == 7)
}

fn get_macrocycle_two_in_same_ring_14_type(
    mol: &ROMol,
    bonds: [&Bond; 3],
    at: Atoms14,
    use_macrocycle_14_config: bool,
) -> TorsionValue {
    let [bnd1, bnd2, bnd3] = bonds;
    let stype = get_atom_stereo(bnd2, at.a1, at.a4);
    if check_macrocycle_two_in_same_ring_amide_ester_14(mol, bnd1, bnd3, at.a1, at.a2, at.a3, at.a4)
        || check_macrocycle_two_in_same_ring_amide_ester_14(
            mol, bnd3, bnd1, at.a4, at.a3, at.a2, at.a1,
        )
    {
        return if use_macrocycle_14_config {
            TorsionValue::forced(TorsionType::Cis)
        } else {
            TorsionValue::of(TorsionType::CisTrans)
        };
    }
    if is_cis(stype) {
        TorsionValue::of(TorsionType::Cis)
    } else if is_trans(stype) {
        TorsionValue::of(TorsionType::Trans)
    } else if both_sp2(mol, at.a2, at.a3) {
        TorsionValue::of(TorsionType::CisTrans)
    } else {
        TorsionValue::of(TorsionType::Flexible)
    }
}

fn get_macrocycle_all_in_same_ring_14_type(
    mol: &ROMol,
    bonds: [&Bond; 3],
    at: Atoms14,
    use_macrocycle_14_config: bool,
) -> Result<TorsionValue> {
    let [bnd1, bnd2, bnd3] = bonds;
    match bnd2.bond_type() {
        BondType::Double => {
            if bnd1.bond_type() == BondType::Double || bnd3.bond_type() == BondType::Double {
                // CC=C=C: the torsion angle is zero
                Ok(TorsionValue::of(TorsionType::Cis))
            } else if bnd2.stereo() > BondStereo::StereoAny {
                if is_cis(get_atom_stereo(bnd2, at.a1, at.a4)) {
                    Ok(TorsionValue::of(TorsionType::Cis))
                } else {
                    Ok(TorsionValue::of(TorsionType::Trans))
                }
            } else if both_sp2(mol, at.a2, at.a3) {
                Ok(TorsionValue::of(TorsionType::CisTrans))
            } else {
                Ok(TorsionValue::of(TorsionType::Flexible))
            }
        }
        BondType::Single => {
            if mol.atom(at.a2).atomic_num() == 16
                && mol.atom(at.a3).atomic_num() == 16
                && mol.atom_degree(at.a2) == 2
                && mol.atom_degree(at.a3) == 2
            {
                // *S-S*
                return Ok(TorsionValue {
                    value: Some(PI / 2.0),
                    ..TorsionValue::of(TorsionType::Custom)
                });
            }
            let reversed = Atoms14 {
                a1: at.a4,
                a2: at.a3,
                a3: at.a2,
                a4: at.a1,
            };
            if check_macrocycle_all_in_same_ring_amide_ester_14(mol, at)
                || check_macrocycle_all_in_same_ring_amide_ester_14(mol, reversed)
            {
                return Ok(if use_macrocycle_14_config {
                    TorsionValue {
                        extra_dist: Some(0.1),
                        ..TorsionValue::forced(TorsionType::Trans)
                    }
                } else {
                    TorsionValue::of(TorsionType::CisTrans)
                });
            }
            if check_amide_ester_15(mol, bnd1, bnd3, at.a2, at.a3)?
                || check_amide_ester_15(mol, bnd3, bnd1, at.a3, at.a2)?
            {
                return Ok(if use_macrocycle_14_config {
                    if mol.atom(at.a2).atomic_num() == 7
                        && mol.atom_degree(at.a2) == 3
                        && mol.atom(at.a1).atomic_num() == 1
                        && mol.total_num_hs(at.a2, true)? == 1
                    {
                        // a secondary amide, this is the H
                        TorsionValue::forced(TorsionType::Cis)
                    } else {
                        TorsionValue::forced(TorsionType::Trans)
                    }
                } else {
                    TorsionValue::of(TorsionType::CisTrans)
                });
            }
            Ok(TorsionValue::of(TorsionType::Flexible))
        }
        _ => Ok(TorsionValue::of(TorsionType::Flexible)),
    }
}

/// Computes the bounds of the 1-4 path through the three bonds and records
/// them in `collected`; they are merged and applied later since a pair of
/// atoms can be reached by several paths. With `collect_all` paths that are
/// not true 1-4 contacts are recorded as well, but without bounds.
#[allow(clippy::too_many_arguments)]
fn collect_14_bounds(
    mol: &ROMol,
    bonds: Bonds14,
    kind: Type14,
    accum: &mut ComputedData,
    mmat: &BoundsMatrix,
    dmat: &[f64],
    info: Optional14Info,
    collected: &mut BTreeMap<usize, Vec<Bounds>>,
    collect_all: bool,
) -> Result<()> {
    let mut set_it = true;
    let Bonds14 {
        b1: bid1,
        b2: bid2,
        b3: bid3,
    } = bonds;
    let (bnd1, bnd2, bnd3) = (mol.bond(bid1), mol.bond(bid2), mol.bond(bid3));
    let na = mol.num_atoms();
    let nb = mol.num_bonds();

    let a2 = accum.bond_adj.get(bid1, bid2) as u32;
    let a3 = accum.bond_adj.get(bid2, bid3) as u32;
    let a1 = bnd1.other_atom_idx(a2)?;
    let a4 = bnd3.other_atom_idx(a3)?;
    let at = Atoms14 { a1, a2, a3, a4 };

    let pid = unified_id2(a1, a4, na);

    // check that the bound was not set before and this actually is a 1-4
    // contact
    if accum.visited_bound(pid, DistType::Dist13)
        || dmat[a1.max(a4) as usize * na as usize + a1.min(a4) as usize] < 2.9
    {
        if collect_all {
            set_it = false;
        } else {
            return Ok(());
        }
    }

    let (u1, u2, u3) = (
        mmat.upper_bound(a1 as usize, a2 as usize),
        mmat.upper_bound(a2 as usize, a3 as usize),
        mmat.upper_bound(a3 as usize, a4 as usize),
    );
    let (l1, l2, l3) = (
        mmat.lower_bound(a1 as usize, a2 as usize),
        mmat.lower_bound(a2 as usize, a3 as usize),
        mmat.lower_bound(a3 as usize, a4 as usize),
    );

    let ba12 = accum.bond_angles.get(bid1, bid2);
    let ba23 = accum.bond_angles.get(bid2, bid3);
    if ba12 <= 0.0 || ba23 <= 0.0 {
        return Err(RDKitError::invariant("bad bond angle"));
    }

    // in fused rings this may not be a 1-4 contact (sf.net bug 2835784)
    if kind > Type14::InChain
        && (mol.bond_idx_between(a1, a3).is_some() || mol.bond_idx_between(a4, a2).is_some())
    {
        if collect_all {
            set_it = false;
        } else {
            return Ok(());
        }
    }

    let bond_refs = [bnd1, bnd2, bnd3];
    let torsion = match kind {
        Type14::InChain => get_chain_14_type(mol, bond_refs, at, info.force_trans_amides)?,
        Type14::InRing => get_in_ring_14_type(mol, bnd2, at, info.ring_size),
        Type14::MacrocycleAllInSameRing => get_macrocycle_all_in_same_ring_14_type(
            mol,
            bond_refs,
            at,
            info.use_macrocycle_14_config,
        )?,
        Type14::MacrocycleTwoInSameRing => get_macrocycle_two_in_same_ring_14_type(
            mol,
            bond_refs,
            at,
            info.use_macrocycle_14_config,
        ),
        // similar to all bonds in one ring; some fine tuning is probably
        // possible if a2 and a3 are not sp2, for now 0-180 degrees is used
        Type14::ShareRingBond | Type14::TwoInDiffRing => get_in_ring_14_type(mol, bnd2, at, 0),
        Type14::TwoInSameRing => get_two_in_same_ring_14_type(mol, bnd2, at, info.prefer_trans),
    };

    let tol12 = scale_tolerance_for_large_sp2(mol, a1, a2, a3);
    let tol23 = scale_tolerance_for_large_sp2(mol, a2, a3, a4);
    let (ba_l12, ba_u12) = get_angle_range(ba12, tol12);
    let (ba_l23, ba_u23) = get_angle_range(ba23, tol23);

    let extra = torsion.extra_dist.unwrap_or(0.0);
    let (mut dl, mut du);
    match torsion.kind {
        TorsionType::Cis => {
            dl = compute_14_dist_cis(l1, l2, l3, ba_l12, ba_l23) + extra;
            du = compute_14_dist_cis(u1, u2, u3, ba_u12, ba_u23) + extra;
            accum.cis_paths.insert(unified_id3(bid1, bid2, bid3, nb));
        }
        TorsionType::Trans => {
            dl = compute_14_dist_trans(l1, l2, l3, ba_l12, ba_l23) + extra;
            du = compute_14_dist_trans(u1, u2, u3, ba_u12, ba_u23) + extra;
            accum.trans_paths.insert(unified_id3(bid1, bid2, bid3, nb));
        }
        TorsionType::CisTrans | TorsionType::Flexible => {
            dl = compute_14_dist_cis(l1, l2, l3, ba_l12, ba_l23);
            du = compute_14_dist_trans(u1, u2, u3, ba_u12, ba_u23);
            // in highly strained situations these can get mixed up
            if du < dl {
                std::mem::swap(&mut du, &mut dl);
            }
            if (du - dl).abs() < DIST12_DELTA {
                dl -= GEN_DIST_TOL;
                du += GEN_DIST_TOL;
            }
        }
        TorsionType::Custom => {
            let value = torsion
                .value
                .ok_or_else(|| RDKitError::invariant("Missing value for custom torsion type"))?;
            dl = compute_14_dist_3d(l1, l2, l3, ba_l12, ba_l23, value);
            du = compute_14_dist_3d(u1, u2, u3, ba_u12, ba_u23, value);
        }
    }

    if (du - dl).abs() < DIST12_DELTA {
        // make sure the lower and upper bounds are not too close
        dl -= GEN_DIST_TOL;
        du += GEN_DIST_TOL;
    }

    accum.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        aid1: a1,
        aid2: a2,
        aid3: a3,
        aid4: a4,
        value: torsion,
        type14: kind,
        r_size: info.ring_size,
    });
    accum.visited14.set(pid);

    if set_it {
        collected
            .entry(pid)
            .or_default()
            .push(Bounds::new(dl, du, a1, a4));
    }
    Ok(())
}

/// Sets the 1-4 bounds from the range of allowed torsion angles.
///
/// 1-4 atoms in the same ring are handled first, then the remaining paths by
/// how their bonds relate to rings, with special cases for cis/trans double
/// bonds, amides and esters, and macrocycles.
fn set14_bounds(
    mol: &ROMol,
    mmat: &mut BoundsMatrix,
    accum: &mut ComputedData,
    dmat: &[f64],
    use_macrocycle_14_config: bool,
    force_trans_amides: bool,
    collect_all: bool,
) -> Result<()> {
    let npt = mmat.num_rows();
    if npt != mol.num_atoms() as usize {
        return Err(RDKitError::invariant("Wrong size metric matrix"));
    }
    // ids of bond triples must fit 64 bits (2.6 million bonds)
    let max_num_bonds = (u64::MAX as f64).powf(1.0 / 3.0) as usize;
    if mol.num_bonds() as usize >= max_num_bonds {
        return Err(RDKitError::Value(
            "Too many bonds in the molecule, cannot compute 1-4 bounds".into(),
        ));
    }
    let rinfo = mol.ring_info();
    let mut bond_rings: Vec<Vec<i32>> = rinfo.bond_rings().to_vec();

    // handle bigger rings first
    libcxx_sort_by(&mut bond_rings, |a, b| a.len() > b.len());

    let mut bid_is_macrocycle: HashSet<u32> = HashSet::new();
    let nb = mol.num_bonds();
    let nbu = nb as usize;
    let mut ring_bond_pairs = BitSet::new(nbu * nbu);
    let mut done_paths = BitSet::new(nbu * nbu * nbu);
    let mut cis_ring_bond_pairs = BitSet::new(nbu * nbu);
    let mut collected: BTreeMap<usize, Vec<Bounds>> = BTreeMap::new();

    // first the 1-4 atoms that belong to the same ring
    for bring in &bond_rings {
        let r_size = bring.len();
        if r_size < 3 {
            continue; // rings with less than 3 bonds are not useful
        }
        let mut bid1 = bring[r_size - 1] as u32;
        for i in 0..r_size {
            let bid2 = bring[i] as u32;
            let bid3 = bring[(i + 1) % r_size] as u32;
            let pid = unified_id2(bid1, bid2, nb);
            let id = unified_id3(bid1, bid2, bid3, nb);
            ring_bond_pairs.set(pid);
            done_paths.set(id);

            let bonds = Bonds14 {
                b1: bid1,
                b2: bid2,
                b3: bid3,
            };
            if r_size > 5 {
                if use_macrocycle_14_config && r_size >= MIN_MACROCYCLE_RING_SIZE {
                    collect_14_bounds(
                        mol,
                        bonds,
                        Type14::MacrocycleAllInSameRing,
                        accum,
                        mmat,
                        dmat,
                        Optional14Info {
                            use_macrocycle_14_config,
                            ..Default::default()
                        },
                        &mut collected,
                        collect_all,
                    )?;
                    bid_is_macrocycle.insert(bid2);
                } else {
                    collect_14_bounds(
                        mol,
                        bonds,
                        Type14::InRing,
                        accum,
                        mmat,
                        dmat,
                        Optional14Info {
                            ring_size: r_size,
                            ..Default::default()
                        },
                        &mut collected,
                        collect_all,
                    )?;
                    cis_ring_bond_pairs.set_to(pid, r_size <= 8);
                }
            } else {
                record_14_path(mol, bonds, accum, r_size)?;
                cis_ring_bond_pairs.set(pid);
            }
            bid1 = bid2;
        }
    }

    for bond in mol.bonds() {
        let bid2 = bond.idx();
        let aid2 = bond.begin_atom_idx();
        let aid3 = bond.end_atom_idx();
        for bnd1 in mol.atom_bonds(aid2) {
            let bid1 = bnd1.idx();
            if bid1 == bid2 {
                continue;
            }
            for bnd3 in mol.atom_bonds(aid3) {
                let bid3 = bnd3.idx();
                if bid3 == bid2 {
                    continue;
                }
                if done_paths.get(unified_id3(bid1, bid2, bid3, nb)) {
                    continue; // dealt with before
                }
                let bonds = Bonds14 {
                    b1: bid1,
                    b2: bid2,
                    b3: bid3,
                };
                let pid1 = unified_id2(bid1, bid2, nb);
                let pid2 = unified_id2(bid2, bid3, nb);

                let (kind, info) = if ring_bond_pairs.get(pid1) || ring_bond_pairs.get(pid2) {
                    // either (bid1, bid2) or (bid2, bid3) are in the same ring
                    // (not all three, which was dealt with before)
                    if use_macrocycle_14_config && bid_is_macrocycle.contains(&bid2) {
                        (
                            Type14::MacrocycleTwoInSameRing,
                            Optional14Info {
                                use_macrocycle_14_config,
                                ..Default::default()
                            },
                        )
                    } else {
                        (
                            Type14::TwoInSameRing,
                            Optional14Info {
                                prefer_trans: cis_ring_bond_pairs.get(pid1)
                                    || cis_ring_bond_pairs.get(pid2),
                                ..Default::default()
                            },
                        )
                    }
                } else if (rinfo.num_bond_rings(bid1)? > 0 && rinfo.num_bond_rings(bid2)? > 0)
                    || (rinfo.num_bond_rings(bid2)? > 0 && rinfo.num_bond_rings(bid3)? > 0)
                {
                    // (bid1, bid2) or (bid2, bid3) are ring bonds of
                    // different rings; the third bond is either external or
                    // in a third ring
                    (Type14::TwoInDiffRing, Optional14Info::default())
                } else if rinfo.num_bond_rings(bid2)? > 0 {
                    // the middle bond is a ring bond and the others do not
                    // belong to the same ring or are not ring bonds
                    (Type14::ShareRingBond, Optional14Info::default())
                } else {
                    // the middle bond is not in a ring
                    (
                        Type14::InChain,
                        Optional14Info {
                            force_trans_amides,
                            ..Default::default()
                        },
                    )
                };
                collect_14_bounds(
                    mol,
                    bonds,
                    kind,
                    accum,
                    mmat,
                    dmat,
                    info,
                    &mut collected,
                    collect_all,
                )?;
            }
        }
    }

    for bounds in collected.into_values() {
        let merged = merge_bounds(bounds)?;
        check_and_set_bounds(
            merged.aid1,
            merged.aid4,
            merged.lower,
            merged.upper,
            mmat,
            false,
        )?;
    }
    Ok(())
}

/// Copies the torsion ranges of the 1-4 paths into the internal coordinates,
/// merging them with the ranges that are already present.
fn collect_torsions(
    mol: &ROMol,
    internal_coords: &mut InternalCoordinates,
    accum: &ComputedData,
) -> Result<()> {
    for path in &accum.paths14 {
        let mut torsion_range = path.to_torsion_range();
        let pid = unified_id3(path.bid1, path.bid2, path.bid3, mol.num_bonds());
        if let Some(given) = internal_coords.torsion_range.get(&pid) {
            torsion_range = torsion_range.merge(given)?;
        }
        // like upstream's `emplace`, a range that is present stays
        internal_coords
            .torsion_range
            .entry(pid)
            .or_insert(torsion_range);
    }
    Ok(())
}

/// Copies bond lengths, angles and torsions from the computed data.
fn collect_internal_coordinates(
    mol: &ROMol,
    internal_coords: &mut InternalCoordinates,
    accum: &ComputedData,
) -> Result<()> {
    let nb = mol.num_bonds();
    internal_coords.lengths.resize(nb as usize, 0.0);
    for bondi in mol.bonds() {
        let bid = bondi.idx();
        internal_coords.lengths[bid as usize] = accum.bond_lengths[bid as usize];
        for bondj in &mol.bonds()[bid as usize + 1..] {
            if adjacent_bond_atoms(bondi, bondj).is_none() {
                continue;
            }
            internal_coords.angles.insert(
                unified_id2(bid, bondj.idx(), nb),
                accum.bond_angles.get(bid, bondj.idx()),
            );
        }
    }
    collect_torsions(mol, internal_coords, accum)
}

/// Fills the off-diagonal elements with the default bounds.
pub fn init_bounds_mat(mmat: &mut BoundsMatrix, default_min: f64, default_max: f64) -> Result<()> {
    for i in 1..mmat.num_rows() {
        for j in 0..i {
            mmat.set_upper_bound(i, j, default_max)?;
            mmat.set_lower_bound(i, j, default_min)?;
        }
    }
    Ok(())
}

/// Options of [`set_topol_bounds`].
#[derive(Debug, Clone, Copy)]
pub struct TopolBoundsOptions {
    /// ignored (kept for compatibility)
    pub scale_vdw: bool,
    pub set15_bounds: bool,
    pub set14_bounds: bool,
    pub set13_bounds: bool,
    /// force field that supplies the 1-2 and 1-3 values
    pub embed_force_field: EmbedFF,
}

impl Default for TopolBoundsOptions {
    fn default() -> Self {
        TopolBoundsOptions {
            scale_vdw: false,
            set15_bounds: true,
            set14_bounds: true,
            set13_bounds: true,
            embed_force_field: EmbedFF::Uff,
        }
    }
}

/// Optional extra results of [`set_topol_bounds`].
#[derive(Debug, Default)]
pub struct TopolBoundsOutput<'a> {
    /// receives the 1-4 paths
    pub paths14: Option<&'a mut Path14Vect>,
    /// receives the internal coordinates, if the embedding mode needs them
    pub internal_coords: Option<&'a mut InternalCoordinates>,
}

/// Sets the distance bounds from the topology: 1-2, 1-3 and 1-4 distances from
/// bond lengths, bond angles and torsion angle ranges, optionally 1-5 bounds
/// for paths that contain rigid 1-4 paths, and finally the lower bounds of the
/// remaining pairs from the van der Waals radii.
///
/// For some strained systems the bounds with 1-5 bounds fail triangle
/// smoothing; then recompute them without 1-5 bounds.
pub fn set_topol_bounds(
    mol: &ROMol,
    mmat: &mut BoundsMatrix,
    params: &EmbedParameters,
    options: &TopolBoundsOptions,
    output: TopolBoundsOutput<'_>,
) -> Result<()> {
    let nb = mol.num_bonds();
    let na = mol.num_atoms();
    if na == 0 {
        return Err(RDKitError::Value("molecule has no atoms".into()));
    }
    let max_num_bonds = (u64::MAX as f64).powf(1.0 / 3.0) as usize;
    if nb as usize >= max_num_bonds {
        return Err(RDKitError::Value(
            "Too many bonds in the molecule, cannot compute 1-4 bounds".into(),
        ));
    }

    let mut accum = ComputedData::new(na as usize, nb as usize);
    let dist_matrix = get_distance_mat(mol);
    let mut mmff_params: Option<Box<MmffMolProperties>> = None;

    let mut parametrized = false;
    if options.embed_force_field == EmbedFF::Mmff {
        let (success, parametrization) =
            set12_bounds(mol, mmat, &mut accum, parametrize_mmff, false)?;
        if success {
            if let Parametrization::Mmff(props) = parametrization {
                mmff_params = Some(props);
            }
            parametrized = true;
        }
    }
    if !parametrized {
        set12_bounds(mol, mmat, &mut accum, parametrize_uff, true)?;
    }

    let internal_coordinate_mode =
        params.initial_embedding_mode == InitialEmbeddingMode::InternalCoordinateEmbedding;
    if options.set13_bounds {
        match &mmff_params {
            Some(props) => set13_bounds_mmff(mol, mmat, &mut accum, props)?,
            None => set13_bounds(mol, mmat, &mut accum)?,
        }
        if options.set14_bounds {
            set14_bounds(
                mol,
                mmat,
                &mut accum,
                &dist_matrix,
                params.use_macrocycle_14_config,
                params.force_trans_amides,
                internal_coordinate_mode,
            )?;
            if options.set15_bounds {
                set15_bounds(mol, mmat, &mut accum, &dist_matrix)?;
            }
        }
    }
    set_lower_bound_vdw(mol, mmat, &dist_matrix)?;

    if internal_coordinate_mode {
        if let Some(coords) = output.internal_coords {
            collect_internal_coordinates(mol, coords, &accum)?;
        }
    }
    if let Some(paths14) = output.paths14 {
        *paths14 = accum.paths14;
    }
    Ok(())
}

/// A bond as its pair of atoms, and an angle as (atom1, atom2, atom3, flag)
/// with the flag set for triple bonds and consecutive double bonds.
pub type BondList = Vec<(i32, i32)>;
pub type AngleList = Vec<[i32; 4]>;

/// The bonds and the angles of the molecule, for the experimental torsion
/// preferences.
pub fn collect_bonds_and_angles(mol: &ROMol) -> (BondList, AngleList) {
    let mut bonds = Vec::with_capacity(mol.num_bonds() as usize);
    let mut angles = Vec::new();
    for bondi in mol.bonds() {
        bonds.push((bondi.begin_atom_idx() as i32, bondi.end_atom_idx() as i32));
        for bondj in &mol.bonds()[bondi.idx() as usize + 1..] {
            let Some((a1, a2, a3)) = adjacent_bond_atoms(bondi, bondj) else {
                continue;
            };
            let mut flag = 0;
            if bondi.bond_type() == BondType::Triple || bondj.bond_type() == BondType::Triple {
                flag = 1;
            } else if bondi.bond_type() == BondType::Double
                && bondj.bond_type() == BondType::Double
                && mol.atom_degree(a2) == 2
            {
                // consecutive double bonds
                flag = 1;
            }
            angles.push([a1 as i32, a2 as i32, a3 as i32, flag]);
        }
    }
    (bonds, angles)
}

/// [`set_topol_bounds`] that also returns the bonds and angles.
pub fn set_topol_bounds_with_bonds_and_angles(
    mol: &ROMol,
    mmat: &mut BoundsMatrix,
    params: &EmbedParameters,
    options: &TopolBoundsOptions,
    output: TopolBoundsOutput<'_>,
) -> Result<(BondList, AngleList)> {
    set_topol_bounds(mol, mmat, params, options, output)?;
    Ok(collect_bonds_and_angles(mol))
}

// the 1-5 helper distances: the first four atoms are cis or trans and the
// fifth atom is placed cis or trans to them

/// Cos of the angle 143, clamped to the valid range.
fn clamp_cos(cval: f64) -> f64 {
    cval.clamp(-1.0, 1.0)
}

/// 1-5 distance with atoms 1-4 cis and atom 5 on the same side:
///
/// ```text
///         5
///          \
///     1     4
///      \   /
///       2-3
/// ```
fn compute15_dists_cis_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = clamp_cos((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14);
    let ang145 = ang34 - cval.acos();
    compute_13_dist(d14, d4, ang145)
}

/// 1-5 distance with atoms 1-4 cis and atom 5 on the other side:
///
/// ```text
///  1     4-5
///   \   /
///    2-3
/// ```
fn compute15_dists_cis_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = clamp_cos((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14);
    let ang145 = ang34 + cval.acos();
    compute_13_dist(d14, d4, ang145)
}

/// 1-5 distance with atoms 1-4 trans and atom 5 on the far side:
///
/// ```text
///  1
///   \
///    2-3
///       \
///        4-5
/// ```
fn compute15_dists_trans_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = clamp_cos((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14);
    let ang145 = ang34 + cval.acos();
    compute_13_dist(d14, d4, ang145)
}

/// 1-5 distance with atoms 1-4 trans and atom 5 on the near side:
///
/// ```text
///                    1
///                     \
///                      2-3
///                         \
///                          4
///                         /
///                        5
/// ```
fn compute15_dists_trans_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = clamp_cos((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14);
    let ang145 = ang34 - cval.acos();
    compute_13_dist(d14, d4, ang145)
}

/// The 1-5 bounds for the atoms at the ends of the 1-4 path through the three
/// bonds and any fourth bond at its last atom.
fn set15_bounds_helper(
    mol: &ROMol,
    bonds: Bonds14,
    kind: TorsionType,
    accum: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
) -> Result<()> {
    let Bonds14 {
        b1: bid1,
        b2: bid2,
        b3: bid3,
    } = bonds;
    let nb = mol.num_bonds();
    let na = mol.num_atoms();

    let aid2 = accum.bond_adj.get(bid1, bid2) as u32;
    let aid1 = mol.bond(bid1).other_atom_idx(aid2)?;
    let aid3 = accum.bond_adj.get(bid2, bid3) as u32;
    let aid4 = mol.bond(bid3).other_atom_idx(aid3)?;
    let d1 = accum.bond_lengths[bid1 as usize];
    let d2 = accum.bond_lengths[bid2 as usize];
    let d3 = accum.bond_lengths[bid3 as usize];
    let ang12 = accum.bond_angles.get(bid1, bid2);
    let ang23 = accum.bond_angles.get(bid2, bid3);

    for i in 0..nb {
        if accum.bond_adj.get(bid3, i) != aid4 as i32 {
            continue;
        }
        let aid5 = mol.bond(i).other_atom_idx(aid4)?;
        // make sure we did not come back to the first atom of the path, which
        // is possible with 4-membered rings (Issue 244)
        let pid = unified_id2(aid1, aid5, na);
        if accum.visited_bound(pid, DistType::Dist14) {
            return Ok(());
        }
        // check that this actually is a 1-5 contact
        if dmat[aid1.max(aid5) as usize * mmat.num_rows() + aid1.min(aid5) as usize] < 3.9 {
            continue;
        }
        if aid1 == aid5 {
            continue;
        }
        if !(mmat.lower_bound(aid1 as usize, aid5 as usize) < DIST12_DELTA
            || accum.set15_atoms.get(pid))
        {
            continue;
        }
        let d4 = accum.bond_lengths[i as usize];
        let ang34 = accum.bond_angles.get(bid3, i);
        let path_id = unified_id3(bid2, bid3, i, nb);
        let is_cis_path = accum.cis_paths.contains(&path_id);
        let is_trans_path = accum.trans_paths.contains(&path_id);
        let forward = |f: fn(f64, f64, f64, f64, f64, f64, f64) -> f64| {
            f(d1, d2, d3, d4, ang12, ang23, ang34)
        };
        let backward = |f: fn(f64, f64, f64, f64, f64, f64, f64) -> f64| {
            f(d4, d3, d2, d1, ang34, ang23, ang12)
        };

        let (mut dl, mut du);
        match kind {
            TorsionType::Cis => {
                if is_cis_path {
                    dl = forward(compute15_dists_cis_cis);
                    du = dl + DIST15_TOL;
                    dl -= DIST15_TOL;
                } else if is_trans_path {
                    dl = forward(compute15_dists_cis_trans);
                    du = dl + DIST15_TOL;
                    dl -= DIST15_TOL;
                } else {
                    dl = forward(compute15_dists_cis_cis) - DIST15_TOL;
                    du = forward(compute15_dists_cis_trans) + DIST15_TOL;
                }
            }
            TorsionType::Trans => {
                if is_cis_path {
                    dl = forward(compute15_dists_trans_cis);
                    du = dl + DIST15_TOL;
                    dl -= DIST15_TOL;
                } else if is_trans_path {
                    dl = forward(compute15_dists_trans_trans);
                    du = dl + DIST15_TOL;
                    dl -= DIST15_TOL;
                } else {
                    dl = forward(compute15_dists_trans_cis) - DIST15_TOL;
                    du = forward(compute15_dists_trans_trans) + DIST15_TOL;
                }
            }
            _ => {
                if is_cis_path {
                    dl = backward(compute15_dists_cis_cis) - DIST15_TOL;
                    du = backward(compute15_dists_cis_trans) + DIST15_TOL;
                } else if is_trans_path {
                    dl = backward(compute15_dists_trans_cis) - DIST15_TOL;
                    du = backward(compute15_dists_trans_trans) + DIST15_TOL;
                } else {
                    let pt = PeriodicTable::get();
                    let vw1 = pt.get_rvdw(mol.atom(aid1).atomic_num())?;
                    let vw5 = pt.get_rvdw(mol.atom(aid5).atomic_num())?;
                    dl = VDW_SCALE_15 * (vw1 + vw5);
                    du = -1.0; // no upper bound
                }
            }
        }
        if du < 0.0 {
            du = MAX_UPPER;
        }
        check_and_set_bounds(aid1, aid5, dl, du, mmat, false)?;
        accum.set15_atoms.set(pid);
    }
    Ok(())
}

/// Sets the 1-5 bounds: for every recorded 1-4 path, in both directions.
fn set15_bounds(
    mol: &ROMol,
    mmat: &mut BoundsMatrix,
    accum: &mut ComputedData,
    dmat: &[f64],
) -> Result<()> {
    for k in 0..accum.paths14.len() {
        let path = accum.paths14[k];
        let kind = path.value.kind;
        set15_bounds_helper(
            mol,
            Bonds14 {
                b1: path.bid1,
                b2: path.bid2,
                b3: path.bid3,
            },
            kind,
            accum,
            mmat,
            dmat,
        )?;
        // going the other way: reverse the 1-4 path
        set15_bounds_helper(
            mol,
            Bonds14 {
                b1: path.bid3,
                b2: path.bid2,
                b3: path.bid1,
            },
            kind,
            accum,
            mmat,
            dmat,
        )?;
    }
    Ok(())
}
