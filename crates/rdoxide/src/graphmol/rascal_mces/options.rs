//! `RascalOptions.h` and `RascalClusterOptions.h`.

/// Options controlling an MCES determination.
#[derive(Debug, Clone, PartialEq)]
pub struct RascalOptions {
    /// If calculated below this, no MCES will be evaluated.
    pub similarity_threshold: f64,
    /// If true, partial aromatic rings won't be returned.
    pub complete_aromatic_rings: bool,
    /// If true, ring bonds won't match non-ring bonds.
    pub ring_matches_ring_only: bool,
    /// If true, only complete rings present in both input molecules' ring
    /// info will be returned. Implies `complete_aromatic_rings` and
    /// `ring_matches_ring_only`.
    pub complete_smallest_rings: bool,
    /// If true, atoms will only match atoms if they have the same number of
    /// explicit connections. E.g. the central atom of `C(C)(C)` won't match
    /// either atom in `CC`.
    pub exact_connections_match: bool,
    /// If true, only return a single fragment for the MCES (the largest of
    /// those of the largest MCES). May not produce the largest possible
    /// single fragment that the molecules have in common.
    pub single_largest_frag: bool,
    /// Minimum number of atoms in any fragment; -1 means no minimum.
    pub min_frag_size: i32,
    /// Biggest through-bond distance that bonds can match; -1 means no limit.
    pub max_frag_separation: i32,
    /// If true, all MCESs are returned, in order of diminishing score. This
    /// is likely to result in higher run times.
    pub all_best_mcess: bool,
    /// Caps the number of equivalent solutions, to keep pathological cases
    /// from exhausting the memory.
    pub max_best_mcess: u32,
    /// Maximum run time in seconds; -1 means no maximum.
    pub timeout: i32,
    /// Might make the code run a bit faster in some circumstances, but on
    /// average it is very marginal.
    pub do_equiv_bond_pruning: bool,
    /// If true, when the similarity thresholds aren't matched still return a
    /// result with the tier 1 and tier 2 similarities filled in.
    pub return_empty_mces: bool,
    /// Too many matching bond (vertex) pairs can cause it to run out of
    /// memory.
    pub max_bond_match_pairs: u32,
    /// SMARTS strings defining atoms that should be considered equivalent,
    /// e.g. `[F,Cl,Br,I]`. A space-separated list allows more than 1 class.
    pub equivalent_atoms: String,
    /// If true, all bonds are treated the same, irrespective of order.
    pub ignore_bond_orders: bool,
    /// If true, atoms are matched just on atomic number; if false, aromatic
    /// and aliphatic atoms are different.
    pub ignore_atom_aromaticity: bool,
    /// Normally the minimum clique size is specified via the similarity
    /// threshold. If this is > 0, it overrides it. It is the minimum number of
    /// BONDS in the MCES.
    pub min_clique_size: u32,
}

impl Default for RascalOptions {
    fn default() -> Self {
        RascalOptions {
            similarity_threshold: 0.7,
            complete_aromatic_rings: true,
            ring_matches_ring_only: false,
            complete_smallest_rings: false,
            exact_connections_match: false,
            single_largest_frag: false,
            min_frag_size: -1,
            max_frag_separation: -1,
            all_best_mcess: false,
            max_best_mcess: 10000,
            timeout: 60,
            do_equiv_bond_pruning: false,
            return_empty_mces: false,
            max_bond_match_pairs: 1000,
            equivalent_atoms: String::new(),
            ignore_bond_orders: false,
            ignore_atom_aromaticity: true,
            min_clique_size: 0,
        }
    }
}

/// Options for the clustering.
#[derive(Debug, Clone, PartialEq)]
pub struct RascalClusterOptions {
    /// Similarity cutoff for clustering. Initial clusters will have molecule
    /// pairs of at least this similarity.
    pub similarity_cutoff: f64,
    /// Penalty score for each unconnected component in the MCES.
    pub a: f64,
    /// Weight of matched bonds over matched atoms.
    pub b: f64,
    /// Minimum number of atoms in a fragment for it to be included in the
    /// MCES (also p in the paper).
    pub min_frag_size: u32,
    /// Two pairs of molecules are included in the same cluster if the
    /// similarity between their MCESs is greater than this (S_a in the
    /// paper).
    pub min_intra_cluster_sim: f64,
    /// Two clusters are merged if the fraction of molecules they have in
    /// common is greater than this (S_b in the paper).
    pub cluster_merge_sim: f64,
    /// The maximum number of fragments in any MCES.
    pub max_num_frags: u32,
    /// The number of threads to use. If > 0, that number; if <= 0, the number
    /// of hardware threads plus this number.
    pub num_threads: i32,
}

impl Default for RascalClusterOptions {
    fn default() -> Self {
        RascalClusterOptions {
            similarity_cutoff: 0.7,
            a: 0.05,
            b: 2.0,
            min_frag_size: 3,
            min_intra_cluster_sim: 0.9,
            cluster_merge_sim: 0.6,
            max_num_frags: 2,
            num_threads: -1,
        }
    }
}
