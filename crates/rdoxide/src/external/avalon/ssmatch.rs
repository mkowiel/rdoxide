//! Port of the substructure matcher of Avalon Toolkit's `ssmatch.c`: `SSMatch`, `SSMatchExt`, `FGMatch`,
//! `RSSMatch` and their helpers (the fingerprint code of that file lives in `fingerprint.rs`).
//!
//! A match is the list of molecule atom indices (0-based) the query atoms `0..n` map to. The C code keeps
//! the matches in a linked list into which new matches are prepended; the functions here return the
//! matches in that list order (the last match found first).

use super::perceive::{
    AROMATIC, CHAIN, DOUBLE, Neighbourhood, RING, SINGLE, TRIPLE, perceive_aromatic_bonds,
    perceive_dy_aromaticity, ring_state, setup_neighbourhood,
};
use super::reaccs::{Atom, Bond, MAXATOMS, MAXBONDS, Molecule, NONE};
use super::utilities::atom_symbol_match;

pub const NOT_ASSIGNED: i32 = -1;
pub const BTF_SINGLE: i32 = 0x01;
pub const BTF_DOUBLE: i32 = 0x02;
pub const BTF_TRIPLE: i32 = 0x04;
pub const BTF_AROMATIC: i32 = 0x08;
pub const BTF_DYAROMATIC: i32 = 0x10;

const SINGLE_DOUBLE: i32 = 5;
const SINGLE_AROMATIC: i32 = 6;
const DOUBLE_AROMATIC: i32 = 7;
const ANY_BOND: i32 = 8;
const ZERO_COUNT: i32 = 1;
const SUB_ZERO: i32 = -1;
const SUB_AS_IS: i32 = -2;
const SUB_MORE: i32 = 6;

/// `SetupNeighbourhood(mp, nbp, nlimit)`: only bonds between atoms with an index below `nlimit` count
fn setup_neighbourhood_limited(mp: &Molecule, nlimit: usize) -> Option<Vec<Neighbourhood>> {
    let mut nb = vec![Neighbourhood::default(); mp.atoms.len()];
    for (i, b) in mp.bonds.iter().enumerate() {
        let (a0, a1) = ((b.atoms[0] - 1) as usize, (b.atoms[1] - 1) as usize);
        if a0 >= nlimit || a1 >= nlimit {
            continue;
        }
        nb[a0].atoms.push(a1);
        nb[a0].bonds.push(i);
        if nb[a0].atoms.len() > 20 {
            return None;
        }
        nb[a1].atoms.push(a0);
        nb[a1].bonds.push(i);
        if nb[a1].atoms.len() > 20 {
            return None;
        }
    }
    Some(nb)
}

/// `AtomMatch`: does the query atom `ssap` (query atom number `atom`, 1-based) match `map`?
fn atom_match(map: &Atom, ssap: &Atom, mol: &Molecule, atom: usize) -> bool {
    if map.rsize_flags & ssap.rsize_flags != ssap.rsize_flags {
        return false;
    }
    let props = (ssap.charge == NONE || map.charge == ssap.charge)
        && (ssap.radical == NONE || map.radical == ssap.radical)
        && (ssap.query_h_count == NONE
            || (ssap.query_h_count == ZERO_COUNT && map.query_h_count == ZERO_COUNT)
            || ssap.query_h_count <= map.query_h_count);
    if ssap.atom_symbol == "L" {
        for sl in &mol.symbol_lists {
            if sl.atom as usize == atom {
                if !props {
                    return false;
                }
                let r = atom_symbol_match(&map.atom_symbol, &sl.string);
                return if sl.logic == 0 { !r } else { r };
            }
        }
        return false;
    }
    props && atom_symbol_match(&map.atom_symbol, &ssap.atom_symbol)
}

/// `BondMatch`: does the query bond `ssbp` match `mbp`?
fn bond_match(mbp: &Bond, ssbp: &Bond) -> bool {
    if ssbp.topography != NONE && ssbp.topography != mbp.topography {
        return false;
    }
    if mbp.rsize_flags & ssbp.rsize_flags != ssbp.rsize_flags {
        return false;
    }
    match ssbp.bond_type {
        SINGLE | DOUBLE | TRIPLE => ssbp.bond_type == mbp.bond_type,
        AROMATIC => {
            ssbp.bond_type == mbp.bond_type || ssbp.bond_type_flags & mbp.bond_type_flags != 0
        }
        SINGLE_DOUBLE => mbp.bond_type == SINGLE || mbp.bond_type == DOUBLE,
        DOUBLE_AROMATIC => mbp.bond_type == AROMATIC || mbp.bond_type == DOUBLE,
        ANY_BOND => true,
        SINGLE_AROMATIC => mbp.bond_type == AROMATIC || mbp.bond_type == SINGLE,
        _ => false,
    }
}

/// the substitution description test of `RecSSMatch` and `Refine`
fn sub_desc_ok(ap1: &Atom, n_lig: usize, ss_lig: usize) -> bool {
    match ap1.sub_desc {
        NONE => true,
        SUB_AS_IS => n_lig == ss_lig,
        SUB_ZERO => n_lig == 0,
        SUB_MORE => n_lig >= SUB_MORE as usize,
        d => n_lig as i32 == d,
    }
}

/// `RecSSMatch`: extends the partial matching `amatch[0..n]`; matches are appended to `found`
/// (in order of discovery; the C list has them reversed)
#[allow(clippy::too_many_arguments)]
fn rec_ss_match(
    mp: &Molecule,
    ssp: &Molecule,
    amatch: &[i32],
    bmatch: &[i32],
    n: usize,
    nlimit: usize,
    single_match: bool,
    nbp: &[Neighbourhood],
    ssnbp: &[Neighbourhood],
    found: &mut Vec<Vec<usize>>,
    signature: Option<&[u64]>,
) {
    if n == nlimit {
        // the bond match must be compatible with either pure standard or pure Daylight aromaticity
        let mut dy_mismatch = false;
        let mut dy_required_match = false;
        for (j, &bm) in bmatch.iter().enumerate().take(ssp.bonds.len()) {
            if bm >= 0 {
                let sbp = &ssp.bonds[j];
                let mbp = &mp.bonds[bm as usize];
                if sbp.bond_type < AROMATIC && sbp.bond_type == mbp.bond_type {
                    continue;
                }
                if sbp.bond_type_flags & BTF_AROMATIC != 0 {
                    if mbp.bond_type_flags & BTF_DYAROMATIC != 0 {
                        dy_required_match = true;
                    }
                } else if mbp.bond_type_flags & BTF_DYAROMATIC != 0 {
                    dy_mismatch = true;
                }
            }
        }
        if dy_required_match && dy_mismatch {
            return;
        }
        found.push(amatch[..n].iter().map(|&a| a as usize).collect());
        return;
    }
    let nm = mp.atoms.len();
    let mut cand = vec![false; nm];
    match signature {
        Some(sig) => {
            for i in 0..nm {
                if sig[i] & (1u64 << (n % 31)) != 0 {
                    cand[i] = true;
                }
            }
        }
        None => cand.iter_mut().for_each(|c| *c = true),
    }
    for &a in &amatch[..n] {
        cand[a as usize] = false;
    }
    // a candidate for an atom linked to the mapped part must neighbour the map of that ligand
    let mut candidates: Vec<usize> = Vec::new();
    for &l in &ssnbp[n].atoms {
        if l < n {
            let lig = amatch[l] as usize;
            for &a in &nbp[lig].atoms {
                if cand[a] {
                    candidates.push(a);
                }
            }
        }
    }
    if !candidates.is_empty() {
        cand.iter_mut().for_each(|c| *c = false);
        for &c in &candidates {
            cand[c] = true;
        }
    }
    for i in 0..nm {
        if !(cand[i]
            && nbp[i].atoms.len() >= ssnbp[n].atoms.len()
            && atom_match(&mp.atoms[i], &ssp.atoms[n], ssp, n + 1))
        {
            continue;
        }
        if !sub_desc_ok(&ssp.atoms[n], nbp[i].atoms.len(), ssnbp[n].atoms.len()) {
            continue;
        }
        let mut new_atom = amatch.to_vec();
        new_atom[n] = i as i32;
        let mut new_bond = bmatch.to_vec();
        let mut ok = true;
        for j in 0..ssnbp[n].atoms.len() {
            let mut false_match = false;
            let sl = ssnbp[n].atoms[j];
            if new_atom[sl] != NOT_ASSIGNED {
                if new_bond[ssnbp[n].bonds[j]] == NOT_ASSIGNED {
                    let mut k = 0;
                    while k < nbp[i].atoms.len() {
                        if nbp[i].atoms[k] as i32 == new_atom[sl] {
                            if bond_match(&mp.bonds[nbp[i].bonds[k]], &ssp.bonds[ssnbp[n].bonds[j]])
                            {
                                new_bond[ssnbp[n].bonds[j]] = nbp[i].bonds[k] as i32;
                                break;
                            }
                            false_match = true;
                        }
                        k += 1;
                    }
                    if k == nbp[i].atoms.len() {
                        false_match = true;
                    }
                }
            }
            if false_match {
                ok = false;
                break;
            }
        }
        if ok {
            rec_ss_match(
                mp,
                ssp,
                &new_atom,
                &new_bond,
                n + 1,
                nlimit,
                single_match,
                nbp,
                ssnbp,
                found,
                signature,
            );
        }
        if !found.is_empty() && single_match {
            break;
        }
    }
}

/// `GoodFGEnvironment`
fn good_fg_environment(
    m: &[usize],
    mp: &Molecule,
    nbp: &[Neighbourhood],
    ssp: &Molecule,
    ssnbp: &[Neighbourhood],
    nlimit: usize,
) -> bool {
    for sindex in 0..nlimit {
        let mindex = m[sindex];
        for i in 0..nbp[mindex].atoms.len() {
            let mligand = nbp[mindex].atoms[i];
            let mut j = 0;
            while j < ssnbp[sindex].atoms.len() {
                let sligand = ssnbp[sindex].atoms[j];
                if sligand < nlimit && m[sligand] == mligand {
                    break;
                }
                if sligand < nlimit {
                    j += 1;
                    continue;
                }
                if atom_match(&mp.atoms[mligand], &ssp.atoms[sligand], ssp, sligand + 1)
                    && bond_match(
                        &mp.bonds[nbp[mindex].bonds[i]],
                        &ssp.bonds[ssnbp[sindex].bonds[j]],
                    )
                {
                    break;
                }
                j += 1;
            }
            if j == ssnbp[sindex].atoms.len() {
                return false;
            }
        }
    }
    true
}

/// `FGMatch`: matches of the first `nlimit` atoms of `ssp`, whose neighbourhood is constrained by the other atoms
pub fn fg_match(mp: &Molecule, ssp: &Molecule, nlimit: usize) -> Vec<Vec<usize>> {
    if mp.atoms.len() > MAXATOMS as usize || ssp.atoms.len() > MAXATOMS as usize {
        return Vec::new();
    }
    let (Some(nbp), Some(ssnbp)) = (
        setup_neighbourhood(mp),
        setup_neighbourhood_limited(ssp, nlimit),
    ) else {
        return Vec::new();
    };
    let mut amatch = vec![NOT_ASSIGNED; ssp.atoms.len()];
    let bmatch = vec![NOT_ASSIGNED; ssp.bonds.len()];
    let mut found = Vec::new();
    for i in 0..mp.atoms.len() {
        if nbp[i].atoms.len() >= ssnbp[0].atoms.len()
            && atom_match(&mp.atoms[i], &ssp.atoms[0], ssp, 1)
        {
            amatch[0] = i as i32;
            rec_ss_match(
                mp, ssp, &amatch, &bmatch, 1, nlimit, false, &nbp, &ssnbp, &mut found, None,
            );
        }
    }
    let Some(full) = setup_neighbourhood(ssp) else {
        return Vec::new();
    };
    // pruning walks the list from its head and prepends the good ones
    let mut kept: Vec<Vec<usize>> = Vec::new();
    for m in found.iter().rev() {
        if good_fg_environment(m, mp, &nbp, ssp, &full, nlimit) {
            kept.push(m.clone());
        }
    }
    // `kept` is in walk order; the resulting list has the last kept match first
    kept.reverse();
    kept
}

/// `RSSMatch`: matches in which query atoms without mapping have no other neighbours than mapped ones
pub fn rss_match(mp: &Molecule, ssp: &Molecule) -> Vec<Vec<usize>> {
    if mp.atoms.len() > MAXATOMS as usize || ssp.atoms.len() > MAXATOMS as usize {
        return Vec::new();
    }
    let (Some(nbp), Some(ssnbp)) = (setup_neighbourhood(mp), setup_neighbourhood(ssp)) else {
        return Vec::new();
    };
    let mut amatch = vec![NOT_ASSIGNED; ssp.atoms.len()];
    let bmatch = vec![NOT_ASSIGNED; ssp.bonds.len()];
    let mut found = Vec::new();
    for i in 0..mp.atoms.len() {
        if nbp[i].atoms.len() >= ssnbp[0].atoms.len()
            && atom_match(&mp.atoms[i], &ssp.atoms[0], ssp, 1)
        {
            amatch[0] = i as i32;
            rec_ss_match(
                mp,
                ssp,
                &amatch,
                &bmatch,
                1,
                ssp.atoms.len(),
                false,
                &nbp,
                &ssnbp,
                &mut found,
                None,
            );
        }
    }
    let mut kept: Vec<Vec<usize>> = Vec::new();
    for m in found.iter().rev() {
        let bad = (0..ssp.atoms.len())
            .any(|i| ssp.atoms[i].mapping == NONE && ssnbp[i].atoms.len() != nbp[m[i]].atoms.len());
        if !bad {
            kept.push(m.clone());
        }
    }
    kept.reverse();
    kept
}

/// `Refine`: candidate atoms for query atom 0 by constraint refinement; fills `is_candidate` and
/// (if given) the match signature; returns the number of candidates
#[allow(clippy::too_many_arguments)]
pub fn refine(
    mp: &Molecule,
    nbp: &[Neighbourhood],
    aa_mask: Option<&[u64]>,
    ssp: &Molecule,
    ssnbp: &[Neighbourhood],
    saa_mask: Option<&[u64]>,
    is_candidate: &mut [bool],
    mut signature: Option<&mut [u64]>,
) -> usize {
    let (nm, ns) = (mp.atoms.len(), ssp.atoms.len());
    if nm > MAXATOMS as usize || ns == 0 || ns > MAXATOMS as usize {
        return 0;
    }
    // compatibility of the bond type counts
    let mut sscounts = [0usize; AROMATIC as usize + 1];
    let mut mcounts = [0usize; AROMATIC as usize + 1];
    for b in &ssp.bonds {
        if b.bond_type <= AROMATIC {
            sscounts[b.bond_type as usize] += 1;
        }
    }
    for b in &mp.bonds {
        if b.bond_type <= AROMATIC {
            mcounts[b.bond_type as usize] += 1;
        }
    }
    for i in 0..AROMATIC as usize {
        if sscounts[i] > mcounts[i] {
            return 0;
        }
    }
    let mut atom_cand = vec![vec![false; nm]; ns];
    let mut ncand = vec![0i32; ns];
    let mut bond_cand = vec![vec![false; mp.bonds.len()]; ssp.bonds.len()];
    let mut refine_ok = false;
    'violation: {
        for i in 0..ns {
            let ap1 = &ssp.atoms[i];
            let sym = ap1.atom_symbol.as_bytes();
            let quick = sym.len() == 1
                && sym[0] != b'A'
                && sym[0] != b'L'
                && sym[0] != b'Q'
                && sym[0] != b'R'
                && ap1.sub_desc == NONE
                && ap1.query_h_count == NONE;
            if quick {
                for j in 0..nm {
                    let ap2 = &mp.atoms[j];
                    let s2 = ap2.atom_symbol.as_bytes();
                    if sym[0] == s2.first().copied().unwrap_or(0)
                        && s2.len() == 1
                        && ap2.rsize_flags & ap1.rsize_flags == ap1.rsize_flags
                        && nbp[j].atoms.len() >= ssnbp[i].atoms.len()
                    {
                        atom_cand[i][j] = true;
                        ncand[i] += 1;
                    }
                }
                if let (Some(saa), Some(aa)) = (saa_mask, aa_mask) {
                    for j in 0..nm {
                        if saa[i] & aa[j] != saa[i] {
                            atom_cand[i][j] = false;
                        }
                    }
                }
            } else {
                for j in 0..nm {
                    let ap2 = &mp.atoms[j];
                    if nbp[j].atoms.len() >= ssnbp[i].atoms.len()
                        && ap2.rsize_flags & ap1.rsize_flags == ap1.rsize_flags
                        && atom_match(ap2, ap1, ssp, i + 1)
                    {
                        let mut c = sub_desc_ok(ap1, nbp[j].atoms.len(), ssnbp[i].atoms.len());
                        if let (Some(saa), Some(aa)) = (saa_mask, aa_mask) {
                            if c && saa[i] & aa[j] != saa[i] {
                                c = false;
                            }
                        }
                        atom_cand[i][j] = c;
                        if c {
                            ncand[i] += 1;
                        }
                    }
                }
            }
            if ncand[i] == 0 {
                break 'violation;
            }
        }
        // possible bond matches
        for i in 0..ssp.bonds.len() {
            let bp1 = &ssp.bonds[i];
            let mut nbond_cand = 0;
            if bp1.bond_type < AROMATIC && bp1.rsize_flags == 0 && bp1.topography == NONE {
                for (j, bp2) in mp.bonds.iter().enumerate() {
                    if bp2.bond_type == bp1.bond_type {
                        bond_cand[i][j] = true;
                        nbond_cand += 1;
                    }
                }
            } else {
                for (j, bp2) in mp.bonds.iter().enumerate() {
                    if (bp1.bond_type >= AROMATIC || bp2.bond_type == bp1.bond_type)
                        && bp2.rsize_flags & bp1.rsize_flags == bp1.rsize_flags
                        && bond_match(bp2, bp1)
                    {
                        bond_cand[i][j] = true;
                        nbond_cand += 1;
                    }
                }
            }
            if nbond_cand <= 0 {
                let (ai1, ai2) = ((bp1.atoms[0] - 1) as usize, (bp1.atoms[1] - 1) as usize);
                for j in 0..nm {
                    atom_cand[ai1][j] = false;
                    atom_cand[ai2][j] = false;
                }
                break 'violation;
            }
        }
        // refinement
        'loop_back: loop {
            let mut changed;
            loop {
                changed = false;
                for i in 0..ssp.bonds.len() {
                    let bp1 = &ssp.bonds[i];
                    let mut tmp1 = vec![false; nm];
                    let mut tmp2 = vec![false; nm];
                    let (ai1, ai2) = ((bp1.atoms[0] - 1) as usize, (bp1.atoms[1] - 1) as usize);
                    for (j, bp2) in mp.bonds.iter().enumerate() {
                        if bond_cand[i][j] {
                            let (m0, m1) =
                                ((bp2.atoms[0] - 1) as usize, (bp2.atoms[1] - 1) as usize);
                            let mut ok = false;
                            if atom_cand[ai1][m0] && atom_cand[ai2][m1] {
                                tmp1[m0] = true;
                                tmp2[m1] = true;
                                ok = true;
                            }
                            if atom_cand[ai1][m1] && atom_cand[ai2][m0] {
                                tmp1[m1] = true;
                                tmp2[m0] = true;
                                ok = true;
                            }
                            if !ok {
                                changed = true;
                                bond_cand[i][j] = false;
                            }
                        }
                    }
                    for j in 0..nm {
                        if atom_cand[ai1][j] && !tmp1[j] {
                            atom_cand[ai1][j] = false;
                            ncand[ai1] -= 1;
                            if ncand[ai1] == 0 {
                                break 'violation;
                            }
                            changed = true;
                        }
                    }
                    for j in 0..nm {
                        if atom_cand[ai2][j] && !tmp2[j] {
                            atom_cand[ai2][j] = false;
                            ncand[ai2] -= 1;
                            if ncand[ai2] == 0 {
                                break 'violation;
                            }
                            changed = true;
                        }
                    }
                }
                if ncand.contains(&0) {
                    break 'violation;
                }
                if !changed {
                    break;
                }
            }
            for i in 0..ns {
                if ncand[i] == 1 {
                    let nmatch = (0..nm).find(|&j| atom_cand[i][j]).unwrap_or(nm);
                    for j in 0..ns {
                        if j != i && nmatch < nm && atom_cand[j][nmatch] {
                            ncand[j] -= 1;
                            if ncand[j] == 0 {
                                break 'violation;
                            }
                            atom_cand[j][nmatch] = false;
                            changed = true;
                        }
                    }
                }
            }
            if changed {
                continue 'loop_back;
            }
            refine_ok = true;
            break;
        }
    }
    let mut result = 0;
    for i in 0..nm {
        if let Some(sig) = signature.as_deref_mut() {
            sig[i] = 0;
            for j in 0..ns {
                if refine_ok && atom_cand[j][i] {
                    sig[i] |= 1u64 << (j % 31);
                }
            }
        }
        if refine_ok && atom_cand[0][i] {
            result += 1;
            is_candidate[i] = true;
        } else {
            is_candidate[i] = false;
        }
    }
    if ncand.contains(&0) {
        result = 0;
    }
    if result != 0 {
        for i in 0..ns {
            let mut tmp = vec![false; nm];
            for &l in &ssnbp[i].atoms {
                for k in 0..nm {
                    if atom_cand[l][k] {
                        tmp[k] = true;
                    }
                }
            }
            if tmp.iter().filter(|&&t| t).count() < ssnbp[i].atoms.len() {
                result = 0;
                break;
            }
        }
    }
    result
}

/// `PerceiveAromaticity` of `ssmatch.c`: the standard and Daylight aromaticity go into `bond_type_flags`
pub fn perceive_aromaticity(mp: &mut Molecule, nbp: &[Neighbourhood]) {
    perceive_aromatic_bonds(mp);
    let h_counts: Vec<i32> = mp.atoms.iter().map(|a| a.query_h_count).collect();
    let bond_types: Vec<i32> = mp.bonds.iter().map(|b| b.bond_type).collect();
    for b in &mut mp.bonds {
        match b.bond_type {
            SINGLE => b.bond_type_flags |= BTF_SINGLE,
            DOUBLE => b.bond_type_flags |= BTF_DOUBLE,
            TRIPLE => b.bond_type_flags |= BTF_TRIPLE,
            AROMATIC => b.bond_type_flags |= BTF_AROMATIC,
            SINGLE_DOUBLE => b.bond_type_flags |= BTF_SINGLE | BTF_DOUBLE,
            SINGLE_AROMATIC => b.bond_type_flags |= BTF_SINGLE | BTF_AROMATIC,
            DOUBLE_AROMATIC => b.bond_type_flags |= BTF_AROMATIC | BTF_DOUBLE,
            _ => {}
        }
    }
    for (a, &h) in mp.atoms.iter_mut().zip(&h_counts) {
        a.query_h_count = h;
    }
    for (b, &t) in mp.bonds.iter_mut().zip(&bond_types) {
        b.bond_type = t;
    }
    perceive_dy_aromaticity(mp, nbp);
    for (b, &t) in mp.bonds.iter_mut().zip(&bond_types) {
        if b.bond_type == AROMATIC && b.bond_type_flags & BTF_AROMATIC == 0 {
            b.bond_type_flags |= BTF_DYAROMATIC;
        }
        if b.bond_type == AROMATIC {
            b.bond_type_flags |= BTF_AROMATIC;
        }
        b.bond_type = t;
    }
    for (a, &h) in mp.atoms.iter_mut().zip(&h_counts) {
        a.query_h_count = h;
    }
}

/// `SSMatchExt`: matches of `ssp` in `mp`. As in C, the aromaticity flags of both molecules (and the
/// topography of `ssp`) are left in the structures.
pub fn ss_match_ext(
    mp: &mut Molecule,
    aa_mask: Option<&[u64]>,
    ssp: &mut Molecule,
    saa_mask: Option<&[u64]>,
    single_match: bool,
) -> Vec<Vec<usize>> {
    if mp.atoms.len() > MAXATOMS as usize
        || mp.bonds.len() > MAXBONDS as usize
        || ssp.atoms.len() > MAXATOMS as usize
    {
        return Vec::new();
    }
    let (Some(nbp), Some(ssnbp)) = (setup_neighbourhood(mp), setup_neighbourhood(ssp)) else {
        return Vec::new();
    };
    let hcounts: Vec<i32> = mp.atoms.iter().map(|a| a.query_h_count).collect();
    let bond_orders: Vec<i32> = mp.bonds.iter().map(|b| b.bond_type).collect();
    let topographies: Vec<i32> = mp.bonds.iter().map(|b| b.topography).collect();
    perceive_aromaticity(mp, &nbp);
    let (_, bond_status) = ring_state(mp);
    for (b, &s) in mp.bonds.iter_mut().zip(&bond_status) {
        b.topography = if s > 0 { RING } else { CHAIN };
    }
    perceive_aromaticity(ssp, &ssnbp);
    let (_, bond_status) = ring_state(ssp);
    for (b, &s) in ssp.bonds.iter_mut().zip(&bond_status) {
        if s > 0 {
            b.topography = RING;
        }
    }
    let restore = |mp: &mut Molecule| {
        for (j, b) in mp.bonds.iter_mut().enumerate() {
            b.bond_type = bond_orders[j];
            b.topography = topographies[j];
        }
        for (a, &h) in mp.atoms.iter_mut().zip(&hcounts) {
            a.query_h_count = h;
        }
    };
    let mut signature = vec![0u64; mp.atoms.len() + 1];
    let mut is_candidate = vec![false; mp.atoms.len()];
    let ncand = refine(
        mp,
        &nbp,
        aa_mask,
        ssp,
        &ssnbp,
        saa_mask,
        &mut is_candidate,
        Some(&mut signature),
    );
    let mut found: Vec<Vec<usize>> = Vec::new();
    if ncand == 0 {
        restore(mp);
        return found;
    }
    let mut amatch = vec![NOT_ASSIGNED; ssp.atoms.len()];
    let bmatch = vec![NOT_ASSIGNED; ssp.bonds.len()];
    for i in 0..mp.atoms.len() {
        if is_candidate[i]
            && nbp[i].atoms.len() >= ssnbp[0].atoms.len()
            && atom_match(&mp.atoms[i], &ssp.atoms[0], ssp, 1)
        {
            amatch[0] = i as i32;
            rec_ss_match(
                mp,
                ssp,
                &amatch,
                &bmatch,
                1,
                ssp.atoms.len(),
                single_match,
                &nbp,
                &ssnbp,
                &mut found,
                Some(&signature),
            );
            if !found.is_empty() && single_match {
                break;
            }
        }
    }
    restore(mp);
    found.reverse();
    found
}

/// `SSMatch`
pub fn ss_match(mp: &mut Molecule, ssp: &mut Molecule, single_match: bool) -> Vec<Vec<usize>> {
    ss_match_ext(mp, None, ssp, None, single_match)
}

/// `PruneDuplicateAtomSets`: drops matches that hit the atom set of an already kept match
pub fn prune_duplicate_atom_sets(matches: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let mut result: Vec<Vec<usize>> = Vec::new();
    for m in matches {
        if !result.iter().any(|kept| m.iter().all(|a| kept.contains(a))) {
            result.push(m);
        }
    }
    result.reverse();
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::reaccsio::mol_str_to_mol;

    fn pm(tag: &str, m: &[Vec<usize>]) -> String {
        let mut s = tag.to_string();
        for x in m {
            s += &format!(
                " [{}]",
                x.iter()
                    .map(|a| a.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            );
        }
        s + "\n"
    }

    fn pflags(tag: &str, m: &Molecule) -> String {
        let mut s = tag.to_string();
        for b in &m.bonds {
            s += &format!(" {}/{}/{}", b.bond_type, b.bond_type_flags, b.topography);
        }
        s + "\n"
    }

    /// the substructure matcher must find what the Avalon C code finds, in the same order
    #[test]
    fn ssmatch_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_ssmatch.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let want = v["output"].as_str().unwrap();
            let (Some(mp), Some(ssp)) = (
                mol_str_to_mol(v["mol"].as_str().unwrap()),
                mol_str_to_mol(v["query"].as_str().unwrap()),
            ) else {
                assert_eq!(want, "NULL\n");
                continue;
            };
            let mut out = String::new();
            let (mut m, mut s) = (mp.clone(), ssp.clone());
            out += &pm("ss", &ss_match(&mut m, &mut s, false));
            out += &pflags("mflags", &m);
            out += &pflags("sflags", &s);
            let (mut m, mut s) = (mp.clone(), ssp.clone());
            out += &pm("single", &ss_match(&mut m, &mut s, true));
            let (mut m, mut s) = (mp.clone(), ssp.clone());
            out += &pm(
                "prune",
                &prune_duplicate_atom_sets(ss_match(&mut m, &mut s, false)),
            );
            out += &pm("rss", &rss_match(&mp, &ssp));
            if ssp.atoms.len() > 1 {
                out += &pm("fg", &fg_match(&mp, &ssp, ssp.atoms.len() - 1));
                out += &pm("fg2", &fg_match(&mp, &ssp, ssp.atoms.len().div_ceil(2)));
            }
            let (m, s) = (mp.clone(), ssp.clone());
            let nbp = setup_neighbourhood(&m).unwrap();
            let snbp = setup_neighbourhood(&s).unwrap();
            let mut cand = vec![false; m.atoms.len()];
            let mut sig = vec![0u64; m.atoms.len() + 1];
            let r = refine(&m, &nbp, None, &s, &snbp, None, &mut cand, Some(&mut sig));
            out += &format!("refine {r}:");
            for i in 0..m.atoms.len() {
                out += &format!(" {}/{}", i32::from(cand[i]), sig[i]);
            }
            out += "\n";
            if out != want {
                bad.push(format!("{}\nwant:\n{want}\ngot:\n{out}", v["name"]));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }
}
