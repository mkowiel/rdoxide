//! Topological distance bounds of small molecules.

use rdoxide::distgeom::{BoundsMatrix, triangle_smooth_bounds};
use rdoxide::graphmol::dist_geom_helpers::{
    EmbedFF, EmbedParameters, InitialEmbeddingMode, InternalCoordinates, Path14Vect,
    TopolBoundsOptions, TopolBoundsOutput, TorsionType, collect_bonds_and_angles, init_bounds_mat,
    set_topol_bounds,
};
use rdoxide::graphmol::mol_ops::add_hs::{AddHsParameters, add_hs};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::{ROMol, RWMol};

fn mol_with_hs(smiles: &str) -> RWMol {
    let mut mol = mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    add_hs(&mut mol, &AddHsParameters::default(), None).unwrap();
    mol
}

fn bounds(mol: &ROMol, params: &EmbedParameters, options: &TopolBoundsOptions) -> BoundsMatrix {
    let mut mmat = BoundsMatrix::new(mol.num_atoms() as usize);
    init_bounds_mat(&mut mmat, 0.0, 1000.0).unwrap();
    set_topol_bounds(
        mol,
        &mut mmat,
        params,
        options,
        TopolBoundsOutput::default(),
    )
    .unwrap();
    mmat
}

fn default_bounds(smiles: &str) -> (RWMol, BoundsMatrix) {
    let mol = mol_with_hs(smiles);
    let mmat = bounds(
        &mol,
        &EmbedParameters::default(),
        &TopolBoundsOptions::default(),
    );
    (mol, mmat)
}

#[test]
fn bounds_are_ordered_and_smoothable() {
    for smiles in [
        "CCO",
        "CC(=O)NC",
        "c1ccccc1C=CC(=O)O",
        "C1CC2CCC1C2",
        "c1ccsc1",
        "C#CC=C=CC",
        "OCC[S][S]CCO",
    ] {
        let (mol, mut mmat) = default_bounds(smiles);
        for i in 1..mol.num_atoms() as usize {
            for j in 0..i {
                assert!(
                    mmat.lower_bound(i, j) <= mmat.upper_bound(i, j),
                    "{smiles}: bounds ({i}, {j})"
                );
                assert!(mmat.lower_bound(i, j) > 0.0, "{smiles}: lower ({i}, {j})");
            }
        }
        assert!(triangle_smooth_bounds(&mut mmat, 0.0), "{smiles}");
    }
}

#[test]
fn ethane_bounds_follow_bond_lengths_and_angles() {
    let (mol, mmat) = default_bounds("CC");
    // atoms: C0 C1 H2 H3 H4 H5 H6 H7
    assert_eq!(mol.num_atoms(), 8);
    // 1-2: UFF C_3-C_3 and C_3-H bond lengths +- 0.01
    let (l, u) = (mmat.lower_bound(0, 1), mmat.upper_bound(0, 1));
    assert!((u - l - 0.02).abs() < 1e-12 && (1.45..1.6).contains(&l));
    let (l, u) = (mmat.lower_bound(0, 2), mmat.upper_bound(0, 2));
    assert!((u - l - 0.02).abs() < 1e-12 && (1.0..1.2).contains(&l));
    // 1-3 H-C-C is about 2.1 A, H-C-H about 1.8 A
    let hcc = (mmat.lower_bound(1, 2) + mmat.upper_bound(1, 2)) / 2.0;
    assert!((2.0..2.3).contains(&hcc), "{hcc}");
    let hch = (mmat.lower_bound(2, 3) + mmat.upper_bound(2, 3)) / 2.0;
    assert!((1.7..1.95).contains(&hch), "{hch}");
    // 1-4 H-C-C-H spans the cis (about 2.5) to the trans (about 3.1) distance
    let (l, u) = (mmat.lower_bound(2, 5), mmat.upper_bound(2, 5));
    assert!(l < 2.6 && u > 3.0 && u < 3.2, "({l}, {u})");
}

#[test]
fn benzene_para_carbons_are_fixed_by_the_ring() {
    let (_, mmat) = default_bounds("c1ccccc1");
    // ring carbons 0 and 3 are para: 2 * 1.4 A
    let (l, u) = (mmat.lower_bound(0, 3), mmat.upper_bound(0, 3));
    assert!(l > 2.6 && u < 2.9 && u - l < 0.25, "({l}, {u})");
    // 1-5 across the ring is fixed as well (meta, a 1-3 distance of 2.42)
    let (l, u) = (mmat.lower_bound(0, 2), mmat.upper_bound(0, 2));
    assert!(l > 2.3 && u < 2.55, "({l}, {u})");
}

#[test]
fn trans_amides_narrow_the_oxygen_hydrogen_bounds() {
    // CC(=O)NC: C0 C1 O2 N3 C4, then H on N is the first hydrogen after the
    // methyl hydrogens
    let mol = mol_with_hs("CC(=O)NC");
    let h_on_n = mol
        .atom_neighbors(3)
        .find(|&n| mol.atom(n).atomic_num() == 1)
        .unwrap() as usize;
    let forced = bounds(
        &mol,
        &EmbedParameters::default(),
        &TopolBoundsOptions::default(),
    );
    let free = bounds(
        &mol,
        &EmbedParameters {
            force_trans_amides: false,
            ..Default::default()
        },
        &TopolBoundsOptions::default(),
    );
    let width = |m: &BoundsMatrix| m.upper_bound(2, h_on_n) - m.lower_bound(2, h_on_n);
    assert!(width(&forced) < 0.2, "{}", width(&forced));
    assert!(width(&free) > 0.5, "{}", width(&free));
    // the H is trans to the O: far apart
    let lower = forced.lower_bound(2, h_on_n);
    assert!(
        lower > 3.0 && lower > free.lower_bound(2, h_on_n),
        "{lower}"
    );
}

#[test]
fn paths_and_internal_coordinates_are_collected() {
    let mol = mol_with_hs("C/C=C/C");
    let params = EmbedParameters {
        initial_embedding_mode: InitialEmbeddingMode::InternalCoordinateEmbedding,
        ..Default::default()
    };
    let mut mmat = BoundsMatrix::new(mol.num_atoms() as usize);
    init_bounds_mat(&mut mmat, 0.0, 1000.0).unwrap();
    let mut paths = Path14Vect::new();
    let mut coords = InternalCoordinates::new(mol.num_bonds() as usize);
    set_topol_bounds(
        &mol,
        &mut mmat,
        &params,
        &TopolBoundsOptions::default(),
        TopolBoundsOutput {
            paths14: Some(&mut paths),
            internal_coords: Some(&mut coords),
        },
    )
    .unwrap();

    assert_eq!(coords.lengths.len(), mol.num_bonds() as usize);
    assert!(coords.lengths.iter().all(|&l| l > 0.9 && l < 1.6));
    let (_, angles) = collect_bonds_and_angles(&mol);
    assert_eq!(coords.angles.len(), angles.len());
    assert!(!paths.is_empty());
    assert_eq!(coords.torsion_range.len(), paths.len());
    // the trans double bond: C0-C1=C2-C3
    let trans = paths
        .iter()
        .find(|p| (p.aid1, p.aid2, p.aid3, p.aid4) == (0, 1, 2, 3))
        .expect("path over the double bond");
    assert_eq!(trans.value.kind, TorsionType::Trans);
}

#[test]
fn mmff_bounds_use_mmff_parameters() {
    let mol = mol_with_hs("CCO");
    let uff = bounds(
        &mol,
        &EmbedParameters::default(),
        &TopolBoundsOptions::default(),
    );
    let mmff = bounds(
        &mol,
        &EmbedParameters::default(),
        &TopolBoundsOptions {
            embed_force_field: EmbedFF::Mmff,
            ..Default::default()
        },
    );
    // MMFF C-O is 1.418 A
    let mid = (mmff.lower_bound(1, 2) + mmff.upper_bound(1, 2)) / 2.0;
    assert!((mid - 1.418).abs() < 0.01, "{mid}");
    assert_ne!(uff.lower_bound(1, 2), mmff.lower_bound(1, 2));
    // an unparametrizable molecule falls back to UFF without failing
    let weird = mol_with_hs("[Zn]CC");
    let _ = bounds(
        &weird,
        &EmbedParameters::default(),
        &TopolBoundsOptions {
            embed_force_field: EmbedFF::Mmff,
            ..Default::default()
        },
    );
}

#[test]
fn empty_molecule_is_rejected() {
    let mol = ROMol::new();
    let mut mmat = BoundsMatrix::new(0);
    let r = set_topol_bounds(
        &mol,
        &mut mmat,
        &EmbedParameters::default(),
        &TopolBoundsOptions::default(),
        TopolBoundsOutput::default(),
    );
    assert!(r.is_err());
}
