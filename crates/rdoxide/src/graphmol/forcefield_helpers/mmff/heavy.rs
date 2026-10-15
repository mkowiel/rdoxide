//! MMFF atom types of heavy atoms (upstream
//! `MMFFMolProperties::setMMFFHeavyAtomType`).
//!
//! Upstream assigns the type inside one long `switch`; each element is a
//! method here that returns the type (0: none), with early returns where the
//! original `break`s out of the switch.

use super::rings::{RingMembershipSize, is_atom_n_oxide};
use crate::error::Result;
use crate::graphmol::atom::ValenceType;
use crate::graphmol::bond::BondType;
use crate::graphmol::romol::ROMol;

pub(super) struct HeavyAtomTyper<'a> {
    mol: &'a ROMol,
    rm: &'a RingMembershipSize,
    /// total degree of every atom
    degree: Vec<u32>,
    /// explicit valence plus implicit hydrogens of every atom
    bond_order: Vec<u32>,
}

impl<'a> HeavyAtomTyper<'a> {
    pub(super) fn new(mol: &'a ROMol, rm: &'a RingMembershipSize) -> Result<Self> {
        let n = mol.num_atoms();
        let mut degree = Vec::with_capacity(n as usize);
        let mut bond_order = Vec::with_capacity(n as usize);
        for i in 0..n {
            degree.push(mol.total_degree(i)?);
            let atom = mol.atom(i);
            bond_order.push(atom.valence(ValenceType::Explicit)? + atom.num_implicit_hs()?);
        }
        Ok(HeavyAtomTyper {
            mol,
            rm,
            degree,
            bond_order,
        })
    }

    fn z(&self, idx: u32) -> u32 {
        self.mol.atom(idx).atomic_num()
    }

    fn nbrs(&self, idx: u32) -> impl Iterator<Item = u32> + 'a {
        self.mol.atom_neighbors(idx)
    }

    fn bond_type(&self, a: u32, b: u32) -> BondType {
        self.mol
            .bond_idx_between(a, b)
            .map_or(BondType::Unspecified, |i| self.mol.bond(i).bond_type())
    }

    fn in_ring_of_size(&self, idx: u32, size: u32) -> bool {
        self.mol
            .ring_info()
            .is_atom_in_ring_of_size(idx, size)
            .unwrap_or(false)
    }

    fn n_oxide(&self, idx: u32) -> Result<bool> {
        is_atom_n_oxide(self.mol, idx)
    }

    /// The type of a heavy atom (0 if there is none).
    pub(super) fn atom_type(&self, idx: u32) -> Result<u32> {
        let mut atom_type = 0;
        if self.mol.atom(idx).is_aromatic() {
            if self.rm.is_atom_in_aromatic_ring_of_size(idx, 5) {
                atom_type = self.aromatic_5_ring(idx)?;
            }
            if atom_type == 0 && self.rm.is_atom_in_aromatic_ring_of_size(idx, 6) {
                atom_type = self.aromatic_6_ring(idx)?;
            }
        }
        if atom_type == 0 {
            atom_type = self.aliphatic(idx)?;
        }
        Ok(atom_type)
    }

    /// Whether `nbr` (in a 5-membered aromatic ring shared with the ipso atom)
    /// is an oxygen, sulfur or non-N-oxide trivalent nitrogen.
    fn is_ring_heteroatom(&self, ipso: u32, nbr: u32) -> Result<bool> {
        Ok(self.rm.are_atoms_in_same_ring_of_size(5, &[ipso, nbr])
            && (self.z(nbr) == 8
                || self.z(nbr) == 16
                || (self.z(nbr) == 7 && self.degree[nbr as usize] == 3 && !self.n_oxide(nbr)?)))
    }

    fn aromatic_5_ring(&self, idx: u32) -> Result<u32> {
        let mut alpha_het = Vec::new();
        let mut beta_het = Vec::new();
        let mut is_alpha_os = false;
        let mut is_beta_os = false;
        let mut alpha_or_beta_in_same_ring = false;
        if self.z(idx) == 6 || self.z(idx) == 7 {
            for nbr in self.nbrs(idx) {
                if !self.rm.is_atom_in_aromatic_ring_of_size(nbr, 5) {
                    continue;
                }
                if self.is_ring_heteroatom(idx, nbr)? {
                    alpha_het.push(nbr);
                }
                for nbr2 in self.nbrs(nbr) {
                    if nbr2 == idx || !self.rm.is_atom_in_aromatic_ring_of_size(nbr2, 5) {
                        continue;
                    }
                    if self.is_ring_heteroatom(idx, nbr2)? {
                        beta_het.push(nbr2);
                    }
                }
            }
            is_alpha_os = alpha_het.iter().any(|&a| matches!(self.z(a), 8 | 16));
            is_beta_os = beta_het.iter().any(|&b| matches!(self.z(b), 8 | 16));
            if !alpha_het.is_empty() && !beta_het.is_empty() {
                alpha_or_beta_in_same_ring = alpha_het.iter().any(|&a| {
                    beta_het
                        .iter()
                        .any(|&b| self.rm.are_atoms_in_same_ring_of_size(5, &[a, b]))
                });
            }
        }
        let (n_alpha, n_beta) = (alpha_het.len(), beta_het.len());
        Ok(match self.z(idx) {
            6 => {
                // no beta heteroatoms: imidazolium carbon (CIM+) needs at
                // least two trivalent nitrogens, one of them in a 5-ring and
                // none in a 6-ring, one of them protonated
                if beta_het.is_empty() {
                    let (mut n_n, mut n_charged, mut in5, mut in6) = (0, 0, 0, 0);
                    for nbr in self.nbrs(idx) {
                        if self.z(nbr) == 7 && self.degree[nbr as usize] == 3 {
                            n_n += 1;
                            if self.mol.atom(nbr).formal_charge() > 0 && !self.n_oxide(nbr)? {
                                n_charged += 1;
                            }
                            if self.rm.is_atom_in_aromatic_ring_of_size(nbr, 5) {
                                in5 += 1;
                            }
                            if self.rm.is_atom_in_aromatic_ring_of_size(nbr, 6) {
                                in6 += 1;
                            }
                        }
                    }
                    if ((n_n == 2 && in5 > 0) || (n_n == 3 && in5 == 2))
                        && n_charged > 0
                        && in6 == 0
                    {
                        return Ok(80);
                    }
                }
                // neither alpha nor beta heteroatoms, or both but in
                // different rings
                if n_alpha == n_beta {
                    let mut surrounded_by_benzene_c = true;
                    let mut surrounded_by_arom = true;
                    for nbr in self.nbrs(idx) {
                        if self.z(nbr) != 6 || !self.in_ring_of_size(nbr, 6) {
                            surrounded_by_benzene_c = false;
                        }
                        if self.rm.are_atoms_in_same_ring_of_size(5, &[idx, nbr])
                            && !self.mol.atom(nbr).is_aromatic()
                        {
                            surrounded_by_arom = false;
                        }
                    }
                    if (n_alpha == 0
                        && n_beta == 0
                        && !surrounded_by_benzene_c
                        && surrounded_by_arom)
                        || (n_alpha > 0
                            && n_beta > 0
                            && (!alpha_or_beta_in_same_ring || (!is_alpha_os && !is_beta_os)))
                    {
                        return Ok(78); // C5
                    }
                }
                if n_alpha > 0 && (n_beta == 0 || is_alpha_os) {
                    return Ok(63); // C5A
                }
                if n_beta > 0 && (n_alpha == 0 || is_beta_os) {
                    return Ok(64); // C5B
                }
                0
            }
            7 => {
                if self.n_oxide(idx)? {
                    return Ok(82); // N5AX, N5BX, N5OX
                }
                if n_alpha == 0 && n_beta == 0 {
                    // pyrrole nitrogen (NPYL) or 5-ring aromatic anion (N5M)
                    return Ok(if self.degree[idx as usize] == 3 {
                        39
                    } else {
                        76
                    });
                }
                if self.degree[idx as usize] == 3 && n_alpha != n_beta {
                    return Ok(81); // NIM+, N5A+, N5B+, N5+
                }
                if n_alpha > 0 && (n_beta == 0 || is_alpha_os) {
                    return Ok(65); // N5A
                }
                if n_beta > 0 && (n_alpha == 0 || is_beta_os) {
                    return Ok(66); // N5B
                }
                if n_alpha > 0 && n_beta > 0 {
                    return Ok(79); // N5
                }
                0
            }
            8 => 59,  // OFUR
            16 => 44, // STHI
            _ => 0,
        })
    }

    fn aromatic_6_ring(&self, idx: u32) -> Result<u32> {
        Ok(match self.z(idx) {
            6 => 37, // CB
            7 => {
                if self.n_oxide(idx)? {
                    69 // NPOX
                } else if self.degree[idx as usize] == 3 {
                    58 // NPD+
                } else {
                    38 // NPYD
                }
            }
            _ => 0,
        })
    }

    fn aliphatic(&self, idx: u32) -> Result<u32> {
        match self.z(idx) {
            6 => Ok(self.aliphatic_carbon(idx)),
            7 => self.aliphatic_nitrogen(idx),
            8 => Ok(self.oxygen(idx)),
            15 => Ok(self.phosphorus(idx)),
            16 => Ok(self.sulfur(idx)),
            17 => Ok(self.chlorine(idx)),
            z => Ok(self.simple_element(idx, z)),
        }
    }

    fn aliphatic_carbon(&self, idx: u32) -> u32 {
        match self.degree[idx as usize] {
            4 => {
                if self.in_ring_of_size(idx, 3) {
                    22 // CR3R
                } else if self.in_ring_of_size(idx, 4) {
                    20 // CR4R
                } else {
                    1 // CR
                }
            }
            3 => {
                let (mut n_n2, mut n_n3, mut n_o, mut n_s) = (0, 0, 0, 0);
                let mut double_bonded_element = 0;
                for nbr in self.nbrs(idx) {
                    let bt = self.bond_type(nbr, idx);
                    if bt == BondType::Double {
                        double_bonded_element = self.z(nbr);
                    }
                    if self.degree[nbr as usize] == 1 {
                        match self.z(nbr) {
                            8 => n_o += 1,
                            16 => n_s += 1,
                            _ => {}
                        }
                    } else if self.z(nbr) == 7 {
                        if self.degree[nbr as usize] == 3 {
                            n_n3 += 1;
                        } else if self.degree[nbr as usize] == 2 && bt == BondType::Double {
                            n_n2 += 1;
                        }
                    }
                }
                if n_n3 >= 2 && n_n2 == 0 && double_bonded_element == 7 {
                    57 // CNN+, CGD+
                } else if n_o == 2 || n_s == 2 {
                    41 // CO2M, CS2M
                } else if self.in_ring_of_size(idx, 4) && double_bonded_element == 6 {
                    30 // CR4E
                } else if matches!(double_bonded_element, 7 | 8 | 15 | 16) {
                    3 // C=N, C=O, C=S, ...
                } else {
                    2 // C=C, CSP2
                }
            }
            2 => 4,  // CSP, =C=
            1 => 60, // C%-
            _ => 0,
        }
    }

    /// Number of neighbors of `idx` that are terminal oxygens.
    fn terminal_oxygens(&self, idx: u32) -> u32 {
        self.nbrs(idx)
            .filter(|&n| self.z(n) == 8 && self.degree[n as usize] == 1)
            .count() as u32
    }

    fn aliphatic_nitrogen(&self, idx: u32) -> Result<u32> {
        let degree = self.degree[idx as usize];
        let bond_order = self.bond_order[idx as usize];
        let mut n_term_o_bonded_to_n = 0;
        let mut is_nso2_or_nso3_or_ncn = false;
        for nbr in self.nbrs(idx) {
            if self.z(nbr) == 8 && self.degree[nbr as usize] == 1 {
                n_term_o_bonded_to_n += 1;
            }
            // a sulfonamide nitrogen has a P or S neighbor with two or more
            // terminal oxygens
            if bond_order >= 3 && matches!(self.z(nbr), 15 | 16) && !is_nso2_or_nso3_or_ncn {
                is_nso2_or_nso3_or_ncn = self.terminal_oxygens(nbr) >= 2;
            }
        }
        if degree == 4 {
            return Ok(if self.n_oxide(idx)? { 68 } else { 34 });
        }
        if degree == 3 {
            if bond_order >= 4 {
                let mut double_bonded_cn = false;
                for nbr in self.nbrs(idx) {
                    if self.bond_type(nbr, idx) == BondType::Double {
                        double_bonded_cn = matches!(self.z(nbr), 6 | 7);
                        if self.z(nbr) == 6 {
                            for nbr2 in self.nbrs(nbr) {
                                if !double_bonded_cn {
                                    break;
                                }
                                if nbr2 == idx {
                                    continue;
                                }
                                double_bonded_cn =
                                    !(self.z(nbr2) == 7 && self.degree[nbr2 as usize] == 3);
                            }
                        }
                    }
                }
                if n_term_o_bonded_to_n == 1 {
                    return Ok(67); // N2OX
                }
                if n_term_o_bonded_to_n >= 2 {
                    return Ok(45); // NO2, NO3
                }
                if double_bonded_cn {
                    return Ok(54); // N+=C, N+=N
                }
            }
            if bond_order >= 3 {
                if let Some(t) = self.trivalent_nitrogen(idx, &mut is_nso2_or_nso3_or_ncn)? {
                    return Ok(t);
                }
            }
        }
        if degree == 2 {
            if bond_order == 4 {
                let isonitrile = self
                    .nbrs(idx)
                    .any(|n| self.bond_type(idx, n) == BondType::Triple);
                return Ok(if isonitrile { 61 } else { 53 }); // NR%, =N=
            }
            if bond_order == 3 {
                let (mut is_nitroso, mut is_imine_or_azo) = (false, false);
                for nbr in self.nbrs(idx) {
                    if self.bond_type(idx, nbr) == BondType::Double {
                        is_nitroso = self.z(nbr) == 8 && n_term_o_bonded_to_n == 1;
                        is_imine_or_azo = matches!(self.z(nbr), 6 | 7);
                    }
                }
                if is_nitroso && !is_imine_or_azo {
                    return Ok(46); // N=O
                }
                if is_imine_or_azo {
                    return Ok(9); // N=C, N=N
                }
            }
            if bond_order >= 2 {
                let is_nso = self
                    .nbrs(idx)
                    .any(|n| self.z(n) == 16 && self.terminal_oxygens(n) == 1);
                if is_nso {
                    return Ok(48); // NSO
                }
                if !is_nso2_or_nso3_or_ncn {
                    return Ok(62); // NM
                }
            }
        }
        // sulfonamide (3 neighbors) or cyano (2 neighbors)
        if is_nso2_or_nso3_or_ncn {
            return Ok(43);
        }
        if degree == 1 {
            let (mut is_nsp, mut is_nazt) = (false, false);
            for nbr in self.nbrs(idx) {
                if is_nsp || is_nazt {
                    break;
                }
                is_nsp = self.bond_type(idx, nbr) == BondType::Triple;
                if self.z(nbr) == 7 && self.degree[nbr as usize] == 2 {
                    is_nazt = self.nbrs(nbr).any(|n2| {
                        (self.z(n2) == 7 && self.degree[n2 as usize] == 2)
                            || (self.z(n2) == 6 && self.degree[n2 as usize] == 3)
                    });
                }
            }
            if is_nsp {
                return Ok(42); // NSP
            }
            if is_nazt {
                return Ok(47); // NAZT
            }
        }
        Ok(8) // NR
    }

    /// The nitrogen rules for three neighbors and a total bond order of at
    /// least 3 (amide, enamine, amidinium, guanidinium, ...).
    #[allow(unused_assignments)] // the carbon counters keep the last neighbor's values
    fn trivalent_nitrogen(
        &self,
        idx: u32,
        is_nso2_or_nso3_or_ncn: &mut bool,
    ) -> Result<Option<u32>> {
        let mut is_nco_or_ncs = false;
        let mut is_ncn_plus = false;
        let mut is_ngd_plus = false;
        let mut is_nnn_or_nnc = false;
        let mut is_nbr_c = false;
        let mut is_nbr_benzene_c = false;
        let mut element_double_bonded_to_c = 0;
        let mut element_triple_bonded_to_c = 0;
        // the carbon counters keep the values of the last carbon neighbor
        let (mut n_n2_bonded_to_c, mut n_n3_bonded_to_c) = (0, 0);
        let (mut n_o_bonded_to_c, mut n_s_bonded_to_c) = (0, 0);
        for nbr in self.nbrs(idx) {
            if self.z(nbr) == 6 {
                is_nbr_c = true;
                if self.mol.atom(nbr).is_aromatic() && self.in_ring_of_size(nbr, 6) {
                    is_nbr_benzene_c = true;
                }
                n_n2_bonded_to_c = 0;
                n_n3_bonded_to_c = 0;
                n_o_bonded_to_c = 0;
                n_s_bonded_to_c = 0;
                let (mut n_formal_charge, mut n_in_aromatic_6_ring) = (0, 0);
                for nbr2 in self.nbrs(nbr) {
                    let bt = self.bond_type(nbr, nbr2);
                    let z2 = self.z(nbr2);
                    if bt == BondType::Double && matches!(z2, 8 | 16) {
                        is_nco_or_ncs = true;
                    }
                    // an aromatic bond also counts for carbon, or for a
                    // nitrogen that is in a single ring
                    if bt == BondType::Double
                        || (bt == BondType::Aromatic
                            && (z2 == 6
                                || (z2 == 7 && self.mol.ring_info().num_atom_rings(nbr2)? == 1)))
                    {
                        element_double_bonded_to_c = z2;
                    }
                    if bt == BondType::Triple {
                        element_triple_bonded_to_c = z2;
                    }
                    if z2 == 7 && self.degree[nbr2 as usize] == 3 {
                        if self.mol.atom(nbr2).formal_charge() == 1 {
                            n_formal_charge += 1;
                        }
                        if self.rm.is_atom_in_aromatic_ring_of_size(nbr, 6) {
                            n_in_aromatic_6_ring += 1;
                        }
                        // fewer than two oxygens: neither nitro nor nitrate
                        let n_o_bonded_to_n3 =
                            self.nbrs(nbr2).filter(|&n3| self.z(n3) == 8).count();
                        if n_o_bonded_to_n3 < 2 {
                            n_n3_bonded_to_c += 1;
                        }
                    }
                    if z2 == 7
                        && self.degree[nbr2 as usize] == 2
                        && matches!(bt, BondType::Double | BondType::Aromatic)
                    {
                        n_n2_bonded_to_c += 1;
                    }
                    if self.mol.atom(nbr2).is_aromatic() {
                        if z2 == 8 {
                            n_o_bonded_to_c += 1;
                        }
                        if z2 == 16 {
                            n_s_bonded_to_c += 1;
                        }
                    }
                }
                if element_double_bonded_to_c == 7 {
                    // amidinium nitrogen (>N-C=N+<)
                    if n_n3_bonded_to_c == 2
                        && n_n2_bonded_to_c == 0
                        && n_formal_charge > 0
                        && n_in_aromatic_6_ring == 0
                        && self.degree[nbr as usize] < 4
                    {
                        is_ncn_plus = true;
                    }
                    // guanidinium nitrogen
                    if n_n3_bonded_to_c == 3 {
                        is_ngd_plus = true;
                    }
                }
            }
            if self.z(nbr) == 7 {
                for nbr2 in self.nbrs(nbr) {
                    if self.bond_type(nbr, nbr2) != BondType::Double {
                        continue;
                    }
                    if self.z(nbr2) == 6 {
                        let (mut n_n, mut n_o, mut n_s) = (0, 0, 0);
                        for nbr3 in self.nbrs(nbr2) {
                            if nbr3 == nbr {
                                continue;
                            }
                            match self.z(nbr3) {
                                7 => n_n += 1,
                                8 => n_o += 1,
                                16 => n_s += 1,
                                _ => {}
                            }
                        }
                        // N=N-C
                        if n_o == 0 && n_s == 0 && n_n == 0 && !is_nbr_benzene_c {
                            is_nnn_or_nnc = true;
                        }
                    }
                    // N=N-N
                    if self.z(nbr2) == 7 && !is_nbr_benzene_c {
                        is_nnn_or_nnc = true;
                    }
                }
            }
        }
        if is_nbr_c {
            // N-C%N
            if element_triple_bonded_to_c == 7 {
                *is_nso2_or_nso3_or_ncn = true;
            }
            if is_ncn_plus {
                return Ok(Some(55)); // NCN+
            }
            if is_ngd_plus {
                return Ok(Some(56)); // NGD+
            }
            // aniline nitrogen, or a nitrogen on a carbon that is double
            // bonded to C, N or P or triple bonded to C
            if !is_nco_or_ncs
                && !*is_nso2_or_nso3_or_ncn
                && ((n_o_bonded_to_c == 0 && n_s_bonded_to_c == 0 && is_nbr_benzene_c)
                    || matches!(element_double_bonded_to_c, 6 | 7 | 15)
                    || element_triple_bonded_to_c == 6)
            {
                return Ok(Some(40)); // NC=C, NC=N, NC=P, NC%C
            }
        }
        if !*is_nso2_or_nso3_or_ncn && (is_nco_or_ncs || is_nnn_or_nnc) {
            return Ok(Some(10)); // NC=O, NC=S, NN=C, NN=N
        }
        Ok(None)
    }

    fn oxygen(&self, idx: u32) -> u32 {
        let degree = self.degree[idx as usize];
        if degree == 3 {
            return 49; // O+
        }
        let atom = self.mol.atom(idx);
        if degree == 2 {
            if self.bond_order[idx as usize] == 3 {
                return 51; // O=+
            }
            let n_h = self.nbrs(idx).filter(|&n| self.z(n) == 1).count() as u32;
            if n_h + atom.num_implicit_hs().unwrap_or(0) == 2 {
                return 70; // OH2
            }
            return 6; // divalent oxygen
        }
        if self.mol.atom_degree(idx) <= 1 {
            return self.terminal_oxygen(idx);
        }
        0
    }

    fn terminal_oxygen(&self, idx: u32) -> u32 {
        let atom = self.mol.atom(idx);
        let (mut n_n, mut n_o, mut n_s) = (0, 0, 0);
        let mut is_oxide_o_bonded_to_h =
            atom.num_explicit_hs() + atom.num_implicit_hs().unwrap_or(0) > 0;
        let (mut is_carboxylate_o, mut is_carbonyl_o, mut is_oxide_o_bonded_to_c) =
            (false, false, false);
        let (mut is_nitroso_o, mut is_oxide_o_bonded_to_n, mut is_n_oxide_o, mut is_nitro_o) =
            (false, false, false, false);
        let (mut is_thio_sulfinate_o, mut is_sulfate_o, mut is_sulfoxide_o) = (false, false, false);
        let mut is_phosphate_or_perchlorate_o = false;
        for nbr in self.nbrs(idx) {
            if is_oxide_o_bonded_to_c
                || is_oxide_o_bonded_to_n
                || is_oxide_o_bonded_to_h
                || is_carboxylate_o
                || is_nitro_o
                || is_n_oxide_o
                || is_thio_sulfinate_o
                || is_sulfate_o
                || is_phosphate_or_perchlorate_o
                || is_carbonyl_o
                || is_nitroso_o
                || is_sulfoxide_o
            {
                break;
            }
            let bt = self.bond_type(idx, nbr);
            let z = self.z(nbr);
            if matches!(z, 6 | 7 | 16) {
                for nbr2 in self.nbrs(nbr) {
                    let (z2, d2) = (self.z(nbr2), self.degree[nbr2 as usize]);
                    if z2 == 7 && d2 == 2 {
                        n_n += 1;
                    }
                    if z2 == 8 && d2 == 1 {
                        n_o += 1;
                    }
                    if z2 == 16 && d2 == 1 {
                        n_s += 1;
                    }
                }
            }
            is_oxide_o_bonded_to_h = z == 1;
            if z == 6 {
                is_carboxylate_o = n_o == 2;
                is_carbonyl_o = bt == BondType::Double;
                is_oxide_o_bonded_to_c = bt == BondType::Single && n_o == 1;
            }
            if z == 7 {
                is_nitroso_o = bt == BondType::Double;
                if bt == BondType::Single && n_o == 1 {
                    let bo = self.bond_order[nbr as usize];
                    is_oxide_o_bonded_to_n = self.degree[nbr as usize] == 2 || bo == 3;
                    is_n_oxide_o = bo == 4;
                }
                is_nitro_o = n_o >= 2;
            }
            if z == 16 {
                is_thio_sulfinate_o = n_s == 1;
                is_sulfate_o = bt == BondType::Single || (bt == BondType::Double && n_o + n_n > 1);
                is_sulfoxide_o = bt == BondType::Double && n_o + n_n == 1;
            }
            is_phosphate_or_perchlorate_o = matches!(z, 15 | 17);
        }
        if is_oxide_o_bonded_to_c || is_oxide_o_bonded_to_n || is_oxide_o_bonded_to_h {
            return 35; // OM, OM2
        }
        if is_carboxylate_o
            || is_nitro_o
            || is_n_oxide_o
            || is_thio_sulfinate_o
            || is_sulfate_o
            || is_phosphate_or_perchlorate_o
        {
            return 32; // O2CM, ONX, O2N, O2S, OP, ...
        }
        if is_carbonyl_o || is_nitroso_o || is_sulfoxide_o {
            return 7; // O=C, O=N, O=S
        }
        0
    }

    fn phosphorus(&self, idx: u32) -> u32 {
        match self.degree[idx as usize] {
            4 => 25, // PO4, PO3, PO2, PO, PTET
            3 => 26, // P
            2 => 75, // -P=C
            _ => 0,
        }
    }

    fn sulfur(&self, idx: u32) -> u32 {
        let degree = self.degree[idx as usize];
        if degree == 3 || degree == 4 {
            let (mut n_o_or_n, mut n_s) = (0, 0);
            let mut is_c_double_bonded_to_s = false;
            for nbr in self.nbrs(idx) {
                if self.z(nbr) == 6 && self.bond_type(idx, nbr) == BondType::Double {
                    is_c_double_bonded_to_s = true;
                }
                let nbr_degree = self.mol.atom_degree(nbr);
                if (nbr_degree == 1 && self.z(nbr) == 8)
                    || (self.degree[nbr as usize] == 2 && self.z(nbr) == 7)
                {
                    n_o_or_n += 1;
                }
                if nbr_degree == 1 && self.z(nbr) == 16 {
                    n_s += 1;
                }
            }
            if (degree == 3 && n_o_or_n == 2 && is_c_double_bonded_to_s) || degree == 4 {
                return 18; // =SO2
            }
            if (n_o_or_n > 0 && n_s > 0) || (n_o_or_n == 2 && !is_c_double_bonded_to_s) {
                return 73; // SSOM
            }
            return 17; // S=O, >S=N
        }
        if degree == 2 {
            let is_o_double_bonded_to_s = self
                .nbrs(idx)
                .any(|n| self.z(n) == 8 && self.bond_type(idx, n) == BondType::Double);
            return if is_o_double_bonded_to_s { 74 } else { 15 }; // =S=O, S
        }
        if self.mol.atom_degree(idx) == 1 {
            let mut n_term_s_bonded_to_nbr = 0;
            let mut is_c_double_bonded_to_s = false;
            for nbr in self.nbrs(idx) {
                n_term_s_bonded_to_nbr += self
                    .nbrs(nbr)
                    .filter(|&n2| self.z(n2) == 16 && self.degree[n2 as usize] == 1)
                    .count();
                if self.z(nbr) == 6 && self.bond_type(idx, nbr) == BondType::Double {
                    is_c_double_bonded_to_s = true;
                }
            }
            // a thioketone, not a dithiocarboxylate
            return if is_c_double_bonded_to_s && n_term_s_bonded_to_nbr != 2 {
                16 // S=C
            } else {
                72 // S-P, SM, SSMO
            };
        }
        0
    }

    fn chlorine(&self, idx: u32) -> u32 {
        let degree = self.degree[idx as usize];
        if degree == 4 && self.nbrs(idx).filter(|&n| self.z(n) == 8).count() == 4 {
            return 77; // CLO4
        }
        if degree == 1 {
            return 12; // Cl
        }
        if self.mol.atom_degree(idx) == 0 {
            return 90; // Cl-
        }
        0
    }

    /// Elements whose type depends on their degree (and charge) only.
    fn simple_element(&self, idx: u32, z: u32) -> u32 {
        let degree = self.mol.atom_degree(idx);
        let charge = self.mol.atom(idx).formal_charge();
        match (z, degree) {
            (3, 0) => 92,  // LI+
            (9, 1) => 11,  // F
            (9, 0) => 89,  // F-
            (11, 0) => 93, // NA+
            (12, 0) => 99, // MG+2
            (14, _) => 19, // SI
            (19, 0) => 94, // K+
            (20, 0) => 96, // CA+2
            (26, 0) if charge == 2 => 87,
            (26, 0) if charge == 3 => 88,
            (29, 0) if charge == 1 => 97,
            (29, 0) if charge == 2 => 98,
            (30, 0) => 95, // ZN+2
            (35, 1) => 13, // Br
            (35, 0) => 91, // BR-
            (53, 1) => 14, // I
            _ => 0,
        }
    }
}
