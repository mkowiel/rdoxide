//! Path and subgraph enumeration (upstream `GraphMol/Subgraphs/{Subgraphs,SubgraphUtils}.{h,cpp}`).

use std::collections::{BTreeMap, HashMap};

use crate::error::{RDKitError, Result};
use crate::graphmol::fingerprints::{hash_combine, hash_vect};
use crate::graphmol::matrices::get_distance_mat;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;

/// `PATH_TYPE`
pub type Path = Vec<i32>;
/// `INT_PATH_LIST_MAP`
pub type PathMap = BTreeMap<i32, Vec<Path>>;

/// `Subgraphs::getNbrsList`: for each (non-H) bond, its neighboring bonds.
fn get_nbrs_list(mol: &ROMol, use_hs: bool) -> Result<BTreeMap<i32, Vec<i32>>> {
    let mut nbrs: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
    for i in 0..mol.num_atoms() {
        if !use_hs && mol.atom(i).atomic_num() == 1 {
            continue;
        }
        let bonds: Vec<_> = mol.atom_bonds(i).collect();
        for b1 in &bonds {
            if !use_hs && mol.atom(b1.other_atom_idx(i)?).atomic_num() == 1 {
                continue;
            }
            let bid1 = b1.idx() as i32;
            nbrs.entry(bid1).or_default();
            for b2 in &bonds {
                let bid2 = b2.idx() as i32;
                if bid1 != bid2 && (use_hs || mol.atom(b2.other_atom_idx(i)?).atomic_num() != 1) {
                    nbrs.get_mut(&bid1).expect("inserted above").push(bid2);
                }
            }
        }
    }
    Ok(nbrs)
}

/// `Subgraphs::recurseWalkRange` (`recurseWalk` is the `lower == upper` case
/// collected into a single list).
fn recurse_walk_range(
    nbrs: &BTreeMap<i32, Vec<i32>>,
    spath: &Path,
    cands: &mut Vec<i32>,
    lower: u32,
    upper: u32,
    mut forbidden: Vec<bool>,
    res: &mut PathMap,
) {
    let nsize = spath.len() as u32;
    if nsize >= lower && nsize <= upper {
        res.entry(nsize as i32).or_default().push(spath.clone());
    }
    if nsize >= upper {
        return;
    }
    while let Some(next) = cands.pop() {
        if !forbidden[next as usize] {
            forbidden[next as usize] = true;
            let mut tstack = cands.clone();
            if let Some(nl) = nbrs.get(&next) {
                for &bid in nl {
                    if !forbidden[bid as usize] {
                        tstack.push(bid);
                    }
                }
            }
            let mut tpath = spath.clone();
            tpath.push(next);
            recurse_walk_range(
                nbrs,
                &tpath,
                &mut tstack,
                lower,
                upper,
                forbidden.clone(),
                res,
            );
        }
    }
}

fn forbidden_bonds(mol: &ROMol, ignore_atoms: Option<&[bool]>) -> Result<Vec<bool>> {
    let mut forbidden = vec![false; mol.num_bonds() as usize];
    if let Some(ia) = ignore_atoms {
        if ia.len() != mol.num_atoms() as usize {
            return Err(RDKitError::precondition("bad ignoreAtoms size"));
        }
        for i in 0..mol.num_bonds() {
            let b = mol.bond(i);
            if ia[b.begin_atom_idx() as usize] || ia[b.end_atom_idx() as usize] {
                forbidden[i as usize] = true;
            }
        }
    }
    Ok(forbidden)
}

fn subgraphs_impl(
    mol: &ROMol,
    lower: u32,
    upper: u32,
    use_hs: bool,
    rooted_at_atom: i32,
    ignore_atoms: Option<&[bool]>,
) -> Result<PathMap> {
    let mut forbidden = forbidden_bonds(mol, ignore_atoms)?;
    let nbrs = get_nbrs_list(mol, use_hs)?;
    let mut res = PathMap::new();
    for (&i, nl) in &nbrs {
        if forbidden[i as usize] {
            continue;
        }
        let b = mol.bond(i as u32);
        if rooted_at_atom >= 0
            && b.begin_atom_idx() != rooted_at_atom as u32
            && b.end_atom_idx() != rooted_at_atom as u32
        {
            continue;
        }
        forbidden[i as usize] = true;
        let mut cands = nl.clone();
        recurse_walk_range(
            &nbrs,
            &vec![i],
            &mut cands,
            lower,
            upper,
            forbidden.clone(),
            &mut res,
        );
    }
    Ok(res)
}

/// `findAllSubgraphsOfLengthN`
pub fn find_all_subgraphs_of_length_n(
    mol: &ROMol,
    target_len: u32,
    use_hs: bool,
    rooted_at_atom: i32,
    ignore_atoms: Option<&[bool]>,
) -> Result<Vec<Path>> {
    let mut m = subgraphs_impl(
        mol,
        target_len,
        target_len,
        use_hs,
        rooted_at_atom,
        ignore_atoms,
    )?;
    Ok(m.remove(&(target_len as i32)).unwrap_or_default())
}

/// `findAllSubgraphsOfLengthsMtoN`
pub fn find_all_subgraphs_of_lengths_m_to_n(
    mol: &ROMol,
    lower: u32,
    upper: u32,
    use_hs: bool,
    rooted_at_atom: i32,
    ignore_atoms: Option<&[bool]>,
) -> Result<PathMap> {
    if lower > upper {
        return Err(RDKitError::precondition(""));
    }
    let mut res = subgraphs_impl(mol, lower, upper, use_hs, rooted_at_atom, ignore_atoms)?;
    for l in lower..=upper {
        res.entry(l as i32).or_default();
    }
    Ok(res)
}

/// `findUniqueSubgraphsOfLengthN`
pub fn find_unique_subgraphs_of_length_n(
    mol: &ROMol,
    target_len: u32,
    use_hs: bool,
    use_bo: bool,
    rooted_at_atom: i32,
    ignore_atoms: Option<&[bool]>,
) -> Result<Vec<Path>> {
    let all =
        find_all_subgraphs_of_length_n(mol, target_len, use_hs, rooted_at_atom, ignore_atoms)?;
    uniquify_paths(mol, &all, use_bo)
}

/// `Subgraphs::extendPaths`
fn extend_paths(
    adj: &[bool],
    dim: usize,
    paths: &[Path],
    allow_ring_closures: i32,
    dist_mat: Option<&[f64]>,
    ignore_atoms: Option<&[bool]>,
) -> Vec<Path> {
    let mut res = Vec::new();
    for path in paths {
        let end_idx = *path.last().expect("paths are never empty") as usize;
        for other in 0..dim {
            if ignore_atoms.is_some_and(|ia| ia[other]) {
                continue;
            }
            if !adj[end_idx * dim + other] {
                continue;
            }
            if let Some(dm) = dist_mat
                && dm[path[0] as usize * dim + other] - (path.len() as f64) < -0.001
            {
                continue;
            }
            if !path.contains(&(other as i32)) {
                let mut p = path.clone();
                p.push(other as i32);
                res.push(p);
            } else if allow_ring_closures > 2
                && path.len() as i32 == allow_ring_closures - 1
                && path[path.len() - 2] != other as i32
            {
                let mut p = path.clone();
                p.push(other as i32);
                res.push(p);
            }
        }
    }
    res
}

/// `Subgraphs::pathFinderHelper`
fn path_finder_helper(
    adj: &[bool],
    dim: usize,
    min_len: u32,
    max_len: u32,
    rooted_at_atom: i32,
    dist_mat: Option<&[f64]>,
    ignore_atoms: Option<&[bool]>,
) -> Result<PathMap> {
    if min_len > max_len {
        return Err(RDKitError::precondition("bad lengths provided"));
    }
    let mut res = PathMap::new();
    let mut paths: Vec<Path> = if rooted_at_atom < 0 {
        (0..dim)
            .filter(|&i| !ignore_atoms.is_some_and(|ia| ia[i]))
            .map(|i| vec![i as i32])
            .collect()
    } else if (rooted_at_atom as usize) < dim
        && !ignore_atoms.is_some_and(|ia| ia[rooted_at_atom as usize])
    {
        vec![vec![rooted_at_atom]]
    } else {
        return Ok(res);
    };
    for length in 1..max_len {
        if length >= min_len {
            res.insert(length as i32, paths.clone());
        }
        paths = extend_paths(adj, dim, &paths, max_len as i32, dist_mat, ignore_atoms);
    }
    res.insert(max_len as i32, paths);
    Ok(res)
}

/// `findAllPathsOfLengthsMtoN`
#[allow(clippy::too_many_arguments)]
pub fn find_all_paths_of_lengths_m_to_n(
    mol: &ROMol,
    mut lower: u32,
    mut upper: u32,
    use_bonds: bool,
    use_hs: bool,
    rooted_at_atom: i32,
    only_shortest_paths: bool,
    ignore_atoms: Option<&[bool]>,
) -> Result<PathMap> {
    if lower > upper {
        return Err(RDKitError::precondition(""));
    }
    let dim = mol.num_atoms() as usize;
    if ignore_atoms.is_some_and(|ia| ia.len() != dim) {
        return Err(RDKitError::precondition("bad ignoreAtoms size"));
    }
    let dist_mat = only_shortest_paths.then(|| get_distance_mat(mol));
    let mut adj = vec![false; dim * dim];
    match &dist_mat {
        None => {
            for i in 0..mol.num_bonds() {
                let b = mol.bond(i);
                let (s, e) = (b.begin_atom_idx() as usize, b.end_atom_idx() as usize);
                if use_hs
                    || (mol.atom(s as u32).atomic_num() != 1
                        && mol.atom(e as u32).atomic_num() != 1)
                {
                    adj[s * dim + e] = true;
                    adj[e * dim + s] = true;
                }
            }
        }
        Some(dm) => {
            for i in 0..dim {
                for j in i + 1..dim {
                    if (dm[i * dim + j] - 1.0).abs() < 1e-4 {
                        adj[i * dim + j] = true;
                        adj[j * dim + i] = true;
                    }
                }
            }
        }
    }
    if use_bonds {
        lower += 1;
        upper += 1;
    }
    let mut atom_paths = path_finder_helper(
        &adj,
        dim,
        lower,
        upper,
        rooted_at_atom,
        dist_mat.as_deref(),
        ignore_atoms,
    )?;
    let mut res = PathMap::new();
    if !use_bonds && lower >= 1 {
        res.insert(1, atom_paths.entry(1).or_default().clone());
    }
    if use_bonds || upper > 1 {
        for i in lower..=upper {
            if i <= 1 {
                continue;
            }
            let mut invars: Vec<Vec<bool>> = Vec::new();
            let paths = atom_paths.entry(i as i32).or_default().clone();
            for resi in &paths {
                let mut invar = vec![false; mol.num_bonds() as usize];
                let mut loc = Vec::with_capacity(i as usize);
                for j in 0..(i - 1) as usize {
                    let bond = mol
                        .get_bond_between_atoms(resi[j] as u32, resi[j + 1] as u32)?
                        .ok_or_else(|| RDKitError::invariant("no bond between path atoms"))?;
                    loc.push(bond.idx() as i32);
                    invar[bond.idx() as usize] = true;
                }
                if !invars.contains(&invar) {
                    invars.push(invar);
                    if use_bonds {
                        res.entry(i as i32 - 1).or_default().push(loc);
                    } else {
                        res.entry(i as i32).or_default().push(resi.clone());
                    }
                }
            }
        }
    }
    Ok(res)
}

/// `findAllPathsOfLengthN`
pub fn find_all_paths_of_length_n(
    mol: &ROMol,
    target_len: u32,
    use_bonds: bool,
    use_hs: bool,
    rooted_at_atom: i32,
    only_shortest_paths: bool,
    ignore_atoms: Option<&[bool]>,
) -> Result<Vec<Path>> {
    let mut m = find_all_paths_of_lengths_m_to_n(
        mol,
        target_len,
        target_len,
        use_bonds,
        use_hs,
        rooted_at_atom,
        only_shortest_paths,
        ignore_atoms,
    )?;
    Ok(m.remove(&(target_len as i32)).unwrap_or_default())
}

/// `findAtomEnvironmentOfRadiusN`; also returns the atom -> distance map.
pub fn find_atom_environment_of_radius_n(
    mol: &ROMol,
    radius: u32,
    rooted_at_atom: u32,
    use_hs: bool,
    enforce_size: bool,
) -> Result<(Path, HashMap<u32, u32>)> {
    if rooted_at_atom >= mol.num_atoms() {
        return Err(RDKitError::Value("bad atom index".into()));
    }
    let mut atom_map = HashMap::new();
    atom_map.insert(rooted_at_atom, 0);
    let mut res = Path::new();
    if radius == 0 {
        return Ok((res, atom_map));
    }
    let mut nbr_stack: Vec<(u32, u32)> = Vec::new();
    for b in mol.atom_bonds(rooted_at_atom) {
        if use_hs || mol.atom(b.other_atom_idx(rooted_at_atom)?).atomic_num() != 1 {
            nbr_stack.push((rooted_at_atom, b.idx()));
        }
    }
    let mut bonds_in = vec![false; mol.num_bonds() as usize];
    let mut i = 0;
    while i < radius {
        if nbr_stack.is_empty() {
            break;
        }
        let mut next_layer = Vec::new();
        for (start, bidx) in nbr_stack.drain(..) {
            if bonds_in[bidx as usize] {
                continue;
            }
            bonds_in[bidx as usize] = true;
            res.push(bidx as i32);
            let o = mol.bond(bidx).other_atom_idx(start)?;
            let e = atom_map.entry(o).or_insert(i + 1);
            *e = (*e).min(i + 1);
            if i < radius - 1 {
                for b in mol.atom_bonds(o) {
                    if !bonds_in[b.idx() as usize]
                        && (use_hs || mol.atom(b.other_atom_idx(o)?).atomic_num() != 1)
                    {
                        next_layer.push((o, b.idx()));
                    }
                }
            }
        }
        nbr_stack = next_layer;
        i += 1;
    }
    if i != radius && enforce_size {
        res.clear();
        atom_map.clear();
    }
    Ok((res, atom_map))
}

/// `Subgraphs::calcPathDiscriminators` -> (pathInvar, path length, atom count)
pub fn calc_path_discriminators(
    mol: &ROMol,
    path: &[i32],
    use_bo: bool,
) -> Result<(u32, usize, usize)> {
    let pt = PeriodicTable::get();
    let mut atoms_used = vec![-1i32; mol.num_atoms() as usize];
    let mut atoms: Vec<u32> = Vec::new();
    let mut path_degrees: Vec<u32> = Vec::new();
    for &b in path {
        let bond = mol.bond(b as u32);
        for a in [bond.begin_atom_idx(), bond.end_atom_idx()] {
            if atoms_used[a as usize] < 0 {
                atoms_used[a as usize] = atoms.len() as i32;
                atoms.push(a);
                path_degrees.push(1);
            } else {
                path_degrees[atoms_used[a as usize] as usize] += 1;
            }
        }
    }
    let n_atoms = atoms.len();
    let mut invars = Vec::with_capacity(n_atoms);
    for (k, &a) in atoms.iter().enumerate() {
        let atom = mol.atom(a);
        let mut invar = atom.atomic_num();
        hash_combine(&mut invar, path_degrees[k]);
        hash_combine(&mut invar, atom.formal_charge() as u32);
        let delta_mass = (atom.mass() - pt.get_atomic_weight(atom.atomic_num())?) as i32;
        hash_combine(&mut invar, delta_mass as u32);
        if atom.is_aromatic() {
            hash_combine(&mut invar, 1);
        }
        invars.push(invar);
    }
    let n_cycles = path.len() / 2 + 1;
    for _ in 0..n_cycles {
        let mut loc: Vec<Vec<u32>> = vec![Vec::new(); n_atoms];
        for &b in path {
            let bond = mol.bond(b as u32);
            let bi = atoms_used[bond.begin_atom_idx() as usize] as usize;
            let ei = atoms_used[bond.end_atom_idx() as usize] as usize;
            let mut v1 = invars[bi];
            let mut v2 = invars[ei];
            if use_bo {
                hash_combine(&mut v1, bond.bond_type() as u32);
                hash_combine(&mut v2, bond.bond_type() as u32);
            }
            loc[bi].push(v2);
            loc[ei].push(v1);
        }
        for (i, l) in loc.iter_mut().enumerate() {
            l.sort_unstable();
            invars[i] = hash_vect(l);
        }
    }
    invars.sort_unstable();
    Ok((hash_vect(&invars), path.len(), n_atoms))
}

/// `Subgraphs::uniquifyPaths`
pub fn uniquify_paths(mol: &ROMol, paths: &[Path], use_bo: bool) -> Result<Vec<Path>> {
    let mut seen = Vec::new();
    let mut res = Vec::new();
    for p in paths {
        let d = calc_path_discriminators(mol, p, use_bo)?;
        if !seen.contains(&d) {
            seen.push(d);
            res.push(p.clone());
        }
    }
    Ok(res)
}
