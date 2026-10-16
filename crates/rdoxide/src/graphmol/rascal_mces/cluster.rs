//! `RascalCluster.cpp` and `RascalButinaCluster.cpp`: clustering with the
//! RASCAL MCES similarity.

use std::collections::{BTreeSet, VecDeque};
use std::sync::Arc;

use super::mces::rascal_mces;
use super::options::{RascalClusterOptions, RascalOptions};
use super::result::RascalResult;
use crate::graphmol::mol_ops::frags::get_mol_frags_as_mols;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::rdgeneral::threads::par_map;
use crate::{RDKitError, Result};

/// An edge of the proximity graph.
#[derive(Debug, Clone, Default)]
pub struct ClusNode {
    pub res: Option<Arc<RascalResult>>,
    pub sim: f64,
    pub mol1_num: u32,
    pub mol2_num: u32,
}

fn calc_mol_mol_similarity(
    i: usize,
    j: usize,
    mols: &[Arc<ROMol>],
    opts: &RascalOptions,
    clus_opts: &RascalClusterOptions,
) -> Result<ClusNode> {
    let mut cn = ClusNode {
        mol1_num: i as u32,
        mol2_num: j as u32,
        ..Default::default()
    };
    let mut res = rascal_mces(&mols[i], &mols[j], opts)?;
    // (if there is no result, the tier 1 and 2 similarities were above the
    // threshold but no MCES was found)
    if let Some(first) = res.first_mut() {
        if !first.bond_matches().is_empty() {
            first.trim_small_frags(3)?;
            first.largest_frags_only(clus_opts.max_num_frags)?;
            cn.sim = first.similarity();
            if cn.sim >= opts.similarity_threshold {
                cn.res = Some(Arc::new(first.clone()));
            }
        }
    }
    Ok(cn)
}

/// The graph of the molecules, with an edge for every pair that has an MCES
/// similar enough.
pub fn build_proximity_graph(
    mols: &[Arc<ROMol>],
    clus_opts: &RascalClusterOptions,
) -> Result<Vec<Vec<ClusNode>>> {
    if mols.len() < 2 {
        return Ok(Vec::new());
    }
    let n = mols.len();
    let opts = RascalOptions {
        similarity_threshold: clus_opts.similarity_cutoff,
        ..Default::default()
    };
    let to_do: Vec<(usize, usize)> = (0..n - 1)
        .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
        .collect();
    let sims = par_map(&to_do, clus_opts.num_threads, |&(i, j)| {
        calc_mol_mol_similarity(i, j, mols, &opts, clus_opts)
    });
    let mut prox_graph = vec![vec![ClusNode::default(); n]; n];
    for cn in sims {
        let cn = cn?;
        let (i, j) = (cn.mol1_num as usize, cn.mol2_num as usize);
        prox_graph[i][j] = cn.clone();
        prox_graph[j][i] = cn;
    }
    Ok(prox_graph)
}

/// The connected subgraphs of the proximity graph.
fn disconnect_proximity_graphs(prox_graph: &[Vec<ClusNode>]) -> Vec<Vec<u32>> {
    let mut sub_graphs = Vec::new();
    let mut done = vec![false; prox_graph.len()];
    while let Some(start) = done.iter().position(|d| !d) {
        let mut nodes: Vec<u32> = Vec::new();
        let mut to_do: VecDeque<u32> = VecDeque::from([start as u32]);
        while let Some(next) = to_do.pop_front() {
            if !done[next as usize] {
                nodes.push(next);
            }
            done[next as usize] = true;
            for i in 0..prox_graph.len() {
                if !done[i] && prox_graph[next as usize][i].res.is_some() {
                    to_do.push_back(i as u32);
                    nodes.push(i as u32);
                    done[i] = true;
                }
            }
        }
        nodes.sort_unstable();
        sub_graphs.push(nodes);
    }
    sub_graphs
}

fn g_ij(mol: &ROMol, a: f64, b: f64, p: u32) -> Result<f64> {
    let (frags, _) = get_mol_frags_as_mols(mol, false, true)?;
    let mut num_big_frags = frags.iter().filter(|mf| mf.num_bonds() > p).count() as i32;
    if num_big_frags == 0 {
        num_big_frags = frags.len() as i32;
    }
    let mut g = f64::from(mol.num_atoms());
    g += b * (1.0 - a * f64::from(num_big_frags - 1)) * f64::from(mol.num_bonds());
    Ok(g)
}

fn make_sub_clusters(
    nbors: &[ClusNode],
    clus_opts: &RascalClusterOptions,
) -> Result<Vec<Vec<u32>>> {
    let mut sub_clusters: Vec<Vec<u32>> = Vec::new();
    let mut tmp_nbors: Vec<Option<&ClusNode>> = nbors.iter().map(Some).collect();
    let mcess = |n: &ClusNode| -> Result<Arc<ROMol>> {
        n.res
            .as_ref()
            .and_then(|r| r.mces_mol().transpose())
            .transpose()?
            .ok_or_else(|| RDKitError::Runtime("an edge has no MCES molecule".into()))
    };

    while let Some(front) = tmp_nbors.first().copied().flatten() {
        sub_clusters.push(vec![front.mol1_num, front.mol2_num]);
        let m1 = mcess(front)?;
        let g_12 = g_ij(&m1, clus_opts.a, clus_opts.b, clus_opts.min_frag_size)?;
        for i in 1..tmp_nbors.len() {
            let Some(nb) = tmp_nbors[i] else { continue };
            let m2 = mcess(nb)?;
            let g_13 = g_ij(&m2, clus_opts.a, clus_opts.b, clus_opts.min_frag_size)?;

            let results = rascal_mces(&m1, &m2, &RascalOptions::default())?;
            let Some(res) = results.first() else { continue };
            if res.bond_matches().is_empty() {
                continue;
            }
            let Some(mces) = res.mces_mol()? else {
                continue;
            };
            let g_12_13 = g_ij(&mces, clus_opts.a, clus_opts.b, clus_opts.min_frag_size)?;
            let sim = g_12_13 / g_12.min(g_13);
            if sim > clus_opts.min_intra_cluster_sim {
                let sc = sub_clusters.last_mut().expect("just pushed");
                sc.push(nb.mol2_num);
                sc.push(nb.mol1_num);
                tmp_nbors[i] = None;
            }
        }
        tmp_nbors[0] = None;
        tmp_nbors.retain(Option::is_some);
        let sc = sub_clusters.last_mut().expect("just pushed");
        sc.sort_unstable();
        sc.dedup();
    }
    Ok(sub_clusters)
}

fn form_initial_clusters(
    sub_graph: &[u32],
    prox_graph: &[Vec<ClusNode>],
    clus_opts: &RascalClusterOptions,
) -> Result<Vec<Vec<u32>>> {
    let mut clusters: Vec<Vec<u32>> = Vec::new();
    if sub_graph.len() < 2 {
        return Ok(clusters);
    }
    for &i in sub_graph {
        let mut nbors: Vec<ClusNode> = sub_graph
            .iter()
            .filter(|&&j| prox_graph[i as usize][j as usize].res.is_some())
            .map(|&j| prox_graph[i as usize][j as usize].clone())
            .collect();
        libcxx_sort_by(&mut nbors, |c1, c2| c1.sim > c2.sim);
        if !nbors.is_empty() {
            clusters.extend(make_sub_clusters(&nbors, clus_opts)?);
        }
    }
    libcxx_sort_by(&mut clusters, |c1, c2| {
        if c1.len() == c2.len() {
            c1[0] < c2[0]
        } else {
            c1.len() > c2.len()
        }
    });
    clusters.dedup();
    Ok(clusters)
}

fn merge_clusters(clusters: &[Vec<u32>], clus_opts: &RascalClusterOptions) -> Vec<Vec<u32>> {
    let mut out: Vec<Vec<u32>> = clusters.to_vec();
    if out.len() < 2 {
        return out;
    }
    let mut i = 0;
    while i + 1 < out.len() {
        for j in i + 1..out.len() {
            let in_common = {
                let a: BTreeSet<_> = out[i].iter().collect();
                out[j].iter().filter(|x| a.contains(x)).count()
            };
            let s = in_common as f64 / (out[i].len() as f64).min(out[j].len() as f64);
            if s > clus_opts.cluster_merge_sim {
                let moved = std::mem::take(&mut out[j]);
                out[i].extend(moved);
                out[i].sort_unstable();
                out[i].dedup();
            }
        }
        out.retain(|c| !c.is_empty());
        i += 1;
    }
    out
}

fn sort_cluster_members_by_mean_sim(prox_graph: &[Vec<ClusNode>], clusters: &mut [Vec<u32>]) {
    for clus in clusters.iter_mut() {
        let mut clus_sims: Vec<(u32, f64)> = Vec::new();
        for i in 0..clus.len() {
            let mut tot = 0.0;
            for j in 0..clus.len() {
                if i != j {
                    tot += prox_graph[clus[i] as usize][clus[j] as usize].sim;
                }
            }
            clus_sims.push((clus[i], tot / (clus.len() as f64 - 1.0)));
        }
        libcxx_sort_by(&mut clus_sims, |p1, p2| p1.1 > p2.1);
        for (c, s) in clus.iter_mut().zip(&clus_sims) {
            *c = s.0;
        }
    }
}

fn make_clusters(
    sub_graphs: &[Vec<u32>],
    prox_graph: &[Vec<ClusNode>],
    clus_opts: &RascalClusterOptions,
) -> Result<Vec<Vec<u32>>> {
    let mut clusters = Vec::new();
    for sg in sub_graphs {
        let these = form_initial_clusters(sg, prox_graph, clus_opts)?;
        clusters.extend(merge_clusters(&these, clus_opts));
    }
    libcxx_sort_by(&mut clusters, |c1, c2| c1.len() > c2.len());
    Ok(clusters)
}

fn collect_singletons(prox_graph: &[Vec<ClusNode>]) -> Vec<u32> {
    (0..prox_graph.len())
        .filter(|&i| prox_graph[i].iter().all(|cn| cn.res.is_none()))
        .map(|i| i as u32)
        .collect()
}

/// Fuzzy clustering with the RASCAL MCES similarity (Raymond et al., J.
/// Chem. Inf. Comput. Sci. 2002). The last cluster holds the molecules that
/// didn't fit into another cluster (the singletons).
pub fn rascal_cluster(
    mols: &[Arc<ROMol>],
    clus_opts: &RascalClusterOptions,
) -> Result<Vec<Vec<u32>>> {
    let prox_graph = build_proximity_graph(mols, clus_opts)?;
    let sub_graphs = disconnect_proximity_graphs(&prox_graph);
    let mut clusters = make_clusters(&sub_graphs, &prox_graph, clus_opts)?;
    let singletons = collect_singletons(&prox_graph);
    clusters.push(singletons);
    sort_cluster_members_by_mean_sim(&prox_graph, &mut clusters);
    Ok(clusters)
}

/// The molecules in the order of the cluster centres: each list is a molecule
/// followed by its neighbours by decreasing similarity.
fn build_nbor_lists(prox_graph: &[Vec<ClusNode>]) -> Vec<Vec<u32>> {
    let mut nbor_lists = Vec::new();
    for (i, row) in prox_graph.iter().enumerate() {
        let mut tmp: Vec<(u32, f64)> = Vec::new();
        for cn in row {
            if cn.res.is_some() {
                if i as u32 == cn.mol1_num {
                    tmp.push((cn.mol2_num, cn.sim));
                } else {
                    tmp.push((cn.mol1_num, cn.sim));
                }
            }
        }
        libcxx_sort_by(&mut tmp, |p1, p2| p1.1 > p2.1);
        let mut list = vec![i as u32];
        list.extend(tmp.iter().map(|p| p.0));
        nbor_lists.push(list);
    }
    sort_nbor_lists(&mut nbor_lists);
    nbor_lists
}

fn sort_nbor_lists(nbor_lists: &mut [Vec<u32>]) {
    libcxx_sort_by(nbor_lists, |nl1, nl2| {
        if nl1.len() == nl2.len() {
            nl1 > nl2
        } else {
            nl1.len() > nl2.len()
        }
    });
}

fn form_clusters(mut nbor_lists: Vec<Vec<u32>>) -> Vec<Vec<u32>> {
    let mut clusters = Vec::new();
    while !nbor_lists.is_empty() {
        clusters.push(nbor_lists[0].clone());
        let in_nbor_list: BTreeSet<u32> = nbor_lists[0].iter().copied().collect();
        nbor_lists[0].clear();
        for nl in nbor_lists.iter_mut() {
            nl.retain(|n| !in_nbor_list.contains(n));
        }
        nbor_lists.retain(|nl| !nl.is_empty());
        sort_nbor_lists(&mut nbor_lists);
    }
    clusters
}

/// Butina clustering (J. Chem. Inf. Comput. Sci. 39, 747 (1999)) with the
/// RASCAL MCES similarity.
pub fn rascal_butina_cluster(
    mols: &[Arc<ROMol>],
    clus_opts: &RascalClusterOptions,
) -> Result<Vec<Vec<u32>>> {
    let prox_graph = build_proximity_graph(mols, clus_opts)?;
    Ok(form_clusters(build_nbor_lists(&prox_graph)))
}
