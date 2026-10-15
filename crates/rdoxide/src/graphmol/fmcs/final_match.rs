//! The checks applied to every complete match of a seed against a target
//! (upstream `FinalMatchCheckFunction` and friends in `FMCS/FMCS.cpp`): ring
//! fusion, chirality and the user-defined checker.

use super::graph::Graph;
use super::params::MCSParameters;
use crate::graphmol::atom::ChiralType;
use crate::graphmol::ring_info::RingInfo;
use crate::graphmol::{BondStereo, BondType, ROMol};
use crate::rdgeneral::utils::count_swaps_to_interconvert;
use crate::{RDKitError, Result};

/// `FinalMatchCheckFunction`: `c1[i]` / `c2[i]` are the query and target graph
/// vertices of the i-th matched pair.
pub fn final_match_check(
    c1: &[u32],
    c2: &[u32],
    mol1: &ROMol,
    query: &Graph,
    mol2: &ROMol,
    target: &Graph,
    p: &MCSParameters,
) -> Result<bool> {
    let bcp = &p.bond_compare_parameters;
    if (bcp.match_fused_rings || bcp.match_fused_rings_strict)
        && !ring_fusion_check(c1, c2, mol1, query, mol2, target, p)?
    {
        return Ok(false);
    }
    if p.atom_compare_parameters.match_chiral_tag
        && !final_chirality_check(c1, c2, mol1, query, mol2, target)?
    {
        return Ok(false);
    }
    if let Some(user) = &p.final_match_checker {
        return user(c1, c2, mol1, query, mol2, target, p);
    }
    Ok(true)
}

#[derive(Default, Clone)]
struct BondCount {
    is_mcs_ring_bond: Vec<bool>,
    is_mcs_ring_bond_non_fused: Vec<bool>,
    non_fused_count_pass1: usize,
    is_mcs_ring_bond_fused: Vec<bool>,
    fused_count_pass1: usize,
}

fn count(bits: &[bool]) -> usize {
    bits.iter().filter(|&&b| b).count()
}

struct RingBondCountVect<'a> {
    mol: &'a ROMol,
    ring_info: &'a RingInfo,
    ring_bond_counts: Vec<BondCount>,
    is_mcs_bond: Vec<bool>,
}

impl<'a> RingBondCountVect<'a> {
    fn new(mol: &'a ROMol) -> Result<Self> {
        let ring_info = mol.ring_info();
        let n_bonds = mol.num_bonds() as usize;
        let template = BondCount {
            is_mcs_ring_bond: vec![false; n_bonds],
            is_mcs_ring_bond_non_fused: vec![false; n_bonds],
            is_mcs_ring_bond_fused: vec![false; n_bonds],
            ..Default::default()
        };
        Ok(RingBondCountVect {
            mol,
            ring_info,
            ring_bond_counts: vec![template; ring_info.num_rings()? as usize],
            is_mcs_bond: vec![false; n_bonds],
        })
    }

    fn set_mcs_bond_bits_pass1(
        &mut self,
        begin: u32,
        end: u32,
        c: &[u32],
        graph: &Graph,
    ) -> Result<()> {
        let bond =
            self.mol
                .bond_idx_between(graph.atom(c[begin as usize]), graph.atom(c[end as usize]))
                .ok_or_else(|| RDKitError::invariant("MCS bond not found"))? as usize;
        self.is_mcs_bond[bond] = true;
        let n_rings = self.ring_info.num_bond_rings(bond as u32)?;
        if n_rings == 0 {
            return Ok(());
        }
        let members = self.ring_info.bond_members(bond as u32)?;
        if n_rings == 1 {
            let ring = &mut self.ring_bond_counts[members[0] as usize];
            ring.is_mcs_ring_bond[bond] = true;
            ring.is_mcs_ring_bond_non_fused[bond] = true;
        } else {
            for &ring_idx in members {
                let ring = &mut self.ring_bond_counts[ring_idx as usize];
                ring.is_mcs_ring_bond[bond] = true;
                ring.is_mcs_ring_bond_fused[bond] = true;
            }
        }
        Ok(())
    }

    fn set_mcs_bond_bits_pass2(&mut self) -> Result<()> {
        for ring in &mut self.ring_bond_counts {
            ring.non_fused_count_pass1 = count(&ring.is_mcs_ring_bond_non_fused);
            ring.fused_count_pass1 = count(&ring.is_mcs_ring_bond_fused);
        }
        for bi in 0..self.mol.num_bonds() as usize {
            if !self.is_mcs_bond[bi] {
                continue;
            }
            let mut fused_bond_ring_idx: Option<usize> = None;
            let mut fused_bond_count = 0;
            for (ring_idx, ring) in self.ring_bond_counts.iter_mut().enumerate() {
                if !ring.is_mcs_ring_bond_fused[bi] {
                    continue;
                }
                if ring.non_fused_count_pass1 == 0
                    && ring.fused_count_pass1 < self.ring_info.bond_rings()[ring_idx].len()
                {
                    ring.is_mcs_ring_bond_fused[bi] = false;
                } else {
                    fused_bond_count += 1;
                    fused_bond_ring_idx = Some(ring_idx);
                }
            }
            if fused_bond_count == 1 {
                let ring = &mut self.ring_bond_counts[fused_bond_ring_idx.expect("one ring")];
                ring.is_mcs_ring_bond_non_fused[bi] = true;
                ring.is_mcs_ring_bond_fused[bi] = false;
            }
        }
        Ok(())
    }

    fn is_ring_fusion_honored(&self) -> Result<bool> {
        for (ring_idx, ring) in self.ring_bond_counts.iter().enumerate() {
            let bond_ring = &self.ring_info.bond_rings()[ring_idx];
            let num_ring_bonds_in_mcs = count(&ring.is_mcs_ring_bond);
            if num_ring_bonds_in_mcs == 0 || num_ring_bonds_in_mcs == bond_ring.len() {
                continue;
            }
            let num_non_fused = count(&ring.is_mcs_ring_bond_non_fused);
            let num_fused = count(&ring.is_mcs_ring_bond_fused);
            let mut num_missing_fused = 0;
            for &bi in bond_ring {
                if self.ring_info.num_bond_rings(bi as u32)? > 1
                    && !ring.is_mcs_ring_bond_non_fused[bi as usize]
                    && !ring.is_mcs_ring_bond_fused[bi as usize]
                {
                    num_missing_fused += 1;
                }
            }
            if num_missing_fused + num_fused + num_non_fused == bond_ring.len() {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

/// `ringFusionCheck`: in permissive mode one of the two molecules may fail
/// the match, in strict mode none.
fn ring_fusion_check(
    c1: &[u32],
    c2: &[u32],
    mol1: &ROMol,
    query: &Graph,
    mol2: &ROMol,
    target: &Graph,
    p: &MCSParameters,
) -> Result<bool> {
    if target.num_edges() < query.num_edges() {
        return Ok(true);
    }
    let mut mol1_counts = RingBondCountVect::new(mol1)?;
    let mut mol2_counts = RingBondCountVect::new(mol2)?;
    for &(begin, end, _) in query.edges() {
        mol1_counts.set_mcs_bond_bits_pass1(begin, end, c1, query)?;
        mol2_counts.set_mcs_bond_bits_pass1(begin, end, c2, target)?;
    }
    mol1_counts.set_mcs_bond_bits_pass2()?;
    mol2_counts.set_mcs_bond_bits_pass2()?;
    let mol1_honored = mol1_counts.is_ring_fusion_honored()?;
    let mol2_honored = mol2_counts.is_ring_fusion_honored()?;
    Ok(if p.bond_compare_parameters.match_fused_rings_strict {
        mol1_honored && mol2_honored
    } else {
        mol1_honored || mol2_honored
    })
}

fn is_tetrahedral(tag: ChiralType) -> bool {
    matches!(tag, ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW)
}

/// `FinalChiralityCheckFunction`
pub fn final_chirality_check(
    c1: &[u32],
    c2: &[u32],
    mol1: &ROMol,
    query: &Graph,
    mol2: &ROMol,
    target: &Graph,
) -> Result<bool> {
    let qna = query.num_vertices();
    for i in 0..qna {
        let a1_idx = query.atom(c1[i]);
        let a1 = mol1.atom(a1_idx);
        let a2_idx = target.atom(c2[i]);
        let a2 = mol2.atom(a2_idx);
        let (ac1, ac2) = (a1.chiral_tag(), a2.chiral_tag());

        // #688: doesn't deal with "explicit" Hs properly
        if mol1.atom_degree(a1_idx) < 3 || !is_tetrahedral(ac1) {
            continue; // skip non chiral center QUERY atoms
        }
        if !is_tetrahedral(ac2) {
            return Ok(false);
        }
        let a1_degree = query.out_degree(c1[i]);
        let a2_degree = mol2.atom_degree(a2_idx);
        if a1_degree > a2_degree {
            if a1_degree == 1 && mol1.atom_degree(a1_idx) == a2_degree {
                continue; // continue to grow the seed
            }
            return Ok(false);
        }

        let mut q_order: Vec<i32> = Vec::new();
        for j in 0..qna {
            if q_order.len() == a1_degree as usize {
                break;
            }
            if let Some(b) = mol1.bond_idx_between(a1_idx, query.atom(c1[j])) {
                q_order.push(b as i32);
            }
        }
        let qmo_order: Vec<i32> = mol1
            .atom_bond_indices(a1_idx)
            .map(|b| b as i32)
            .filter(|b| q_order.contains(b))
            .collect();
        let q_perm_count = count_swaps_to_interconvert(&qmo_order, &q_order)?;

        let mut m_order: Vec<i32> = Vec::new();
        for j in 0..qna {
            if m_order.len() == a2_degree as usize {
                break;
            }
            if let Some(b) = mol2.bond_idx_between(a2_idx, target.atom(c2[j])) {
                m_order.push(b as i32);
            }
        }
        while m_order.len() < a2_degree as usize {
            m_order.push(-1);
        }
        let mo_order: Vec<i32> = mol2
            .atom_bond_indices(a2_idx)
            .map(|b| {
                if m_order.contains(&(b as i32)) {
                    b as i32
                } else {
                    -1
                }
            })
            .collect();
        let m_perm_count = count_swaps_to_interconvert(&mo_order, &m_order)?;

        if (q_perm_count % 2 == m_perm_count % 2 && a1.chiral_tag() != a2.chiral_tag())
            || (q_perm_count % 2 != m_perm_count % 2 && a1.chiral_tag() == a2.chiral_tag())
        {
            return Ok(false);
        }
    }

    // double bond stereo
    let mut q_map = std::collections::BTreeMap::new();
    for j in 0..qna {
        q_map.insert(query.atom(c1[j]), j);
    }
    for &(_, _, bond_idx) in query.edges() {
        let q_bnd = mol1.bond(bond_idx);
        if q_bnd.bond_type() != BondType::Double || q_bnd.stereo() <= BondStereo::StereoAny {
            continue;
        }
        let q_stereo = q_bnd.stereo_atoms();
        if q_stereo.len() != 2 {
            continue; // MUST check it in the seed, not in the molecule
        }
        // like std::map::operator[], atoms outside the seed map to position 0
        let t_atom =
            |query_atom: u32| target.atom(c2[q_map.get(&query_atom).copied().unwrap_or(0)]);
        let m_bond_idx = mol2
            .bond_idx_between(t_atom(q_bnd.begin_atom_idx()), t_atom(q_bnd.end_atom_idx()))
            .ok_or_else(|| RDKitError::invariant("Matching bond not found"))?;
        let m_bnd = mol2.bond(m_bond_idx);
        if m_bnd.bond_type() != BondType::Double || m_bnd.stereo() <= BondStereo::StereoAny {
            continue;
        }
        let m_stereo = m_bnd.stereo_atoms();
        if m_stereo.len() != 2 {
            continue;
        }
        let (mut end1_matches, mut end2_matches) = (0, 0);
        let hit = |q: i32, m: i32| t_atom(q as u32) as i32 == m;
        if t_atom(q_bnd.begin_atom_idx()) == m_bnd.begin_atom_idx() {
            if hit(q_stereo[0], m_stereo[0]) {
                end1_matches = 1;
            }
            if hit(q_stereo[1], m_stereo[1]) {
                end2_matches = 1;
            }
        } else {
            if hit(q_stereo[0], m_stereo[1]) {
                end1_matches = 1;
            }
            if hit(q_stereo[1], m_stereo[0]) {
                end2_matches = 1;
            }
        }
        let one_end = end1_matches + end2_matches == 1;
        if m_bnd.stereo() == q_bnd.stereo() && one_end {
            return Ok(false);
        }
        if m_bnd.stereo() != q_bnd.stereo() && !one_end {
            return Ok(false);
        }
    }
    Ok(true)
}
