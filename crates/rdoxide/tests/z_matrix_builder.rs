//! Z-matrices built from the internal coordinates of a molecule.

use rdoxide::distgeom::{BoundsMatrix, ZMatrix, compute_z_matrix_coords};
use rdoxide::geometry::Point3D;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::dist_geom_helpers::{
    EmbedParameters, InitialEmbeddingMode, InternalCoordinates, TopolBoundsOptions,
    TopolBoundsOutput, correct_chiral_centers, init_bounds_mat, set_molecule_dfs,
    set_molecule_dfs_from, set_topol_bounds,
};
use rdoxide::graphmol::mol_ops::add_hs::{AddHsParameters, add_hs};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::random::MinStdRand;

fn prepare(smiles: &str) -> (RWMol, BoundsMatrix, InternalCoordinates) {
    let mut mol = mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    add_hs(&mut mol, &AddHsParameters::default(), None).unwrap();
    let mut mmat = BoundsMatrix::new(mol.num_atoms() as usize);
    init_bounds_mat(&mut mmat, 0.0, 1000.0).unwrap();
    let mut coords = InternalCoordinates::new(mol.num_bonds() as usize);
    let params = EmbedParameters {
        initial_embedding_mode: InitialEmbeddingMode::InternalCoordinateEmbedding,
        ..Default::default()
    };
    set_topol_bounds(
        &mol,
        &mut mmat,
        &params,
        &TopolBoundsOptions::default(),
        TopolBoundsOutput {
            paths14: None,
            internal_coords: Some(&mut coords),
        },
    )
    .unwrap();
    (mol, mmat, coords)
}

fn place(mol: &RWMol, zmat: &ZMatrix, seed: u32) -> Vec<Point3D> {
    let mut positions = vec![Point3D::new(0.0, 0.0, 0.0); mol.num_atoms() as usize];
    let mut rng = MinStdRand::new(seed);
    compute_z_matrix_coords(zmat, &mut positions, &mut rng).unwrap();
    positions
}

#[test]
fn z_matrix_reproduces_bond_lengths_and_angles() {
    for smiles in [
        "CC(=O)NC",
        "C1CCCCC1",
        "c1ccccc1C",
        "N[C@@H](C)C(=O)O",
        "C#CC=C",
    ] {
        let (mol, mmat, coords) = prepare(smiles);
        let n = mol.num_atoms() as usize;
        let mut zmat = ZMatrix::new(n);
        set_molecule_dfs(&mol, &mut zmat, &coords).unwrap();
        assert_eq!(zmat.len(), n, "{smiles}");
        correct_chiral_centers(&mol, &mut zmat).unwrap();
        let pos = place(&mol, &zmat, 7);
        // each atom is at its bond length from the atom it is placed from
        for row in zmat.rows().skip(1) {
            let r = row.atom_idx;
            let b = row.internal.bond_ref.unwrap();
            let d = (pos[r] - pos[b]).length();
            assert!(
                d > mmat.lower_bound(r, b) - 0.02 && d < mmat.upper_bound(r, b) + 0.02,
                "{smiles}: ({r}, {b}) {d}"
            );
        }
        // and at its bond angle from the angle reference
        for row in zmat.rows().skip(2) {
            let (r, b, a) = (
                row.atom_idx,
                row.internal.bond_ref.unwrap(),
                row.internal.angle_ref.unwrap(),
            );
            let (v1, v2) = (pos[r] - pos[b], pos[a] - pos[b]);
            let angle = (v1.dot_product(&v2) / (v1.length() * v2.length())).acos();
            assert!(
                (angle - row.internal.angle.unwrap()).abs() < 1e-3,
                "{smiles}: angle at {b}"
            );
        }
    }
}

#[test]
fn chirality_is_reproduced() {
    // L-alanine: the two enantiomers must come out mirrored
    let signed_volume = |smiles: &str| {
        let (mol, _, coords) = prepare(smiles);
        let n = mol.num_atoms() as usize;
        let mut zmat = ZMatrix::new(n);
        set_molecule_dfs(&mol, &mut zmat, &coords).unwrap();
        correct_chiral_centers(&mol, &mut zmat).unwrap();
        let pos = place(&mol, &zmat, 11);
        // the centre is atom 1; its neighbours N0, C2, C3 (carboxyl)
        let c = pos[1];
        let (a, b, d) = (pos[0] - c, pos[2] - c, pos[3] - c);
        a.dot_product(&b.cross_product(&d))
    };
    let l = signed_volume("N[C@@H](C)C(=O)O");
    let d = signed_volume("N[C@H](C)C(=O)O");
    assert!(l * d < 0.0, "{l} {d}");
}

#[test]
fn single_atoms_and_pairs() {
    let (mol, _, coords) = prepare("[Na+]");
    let mut zmat = ZMatrix::new(1);
    set_molecule_dfs(&mol, &mut zmat, &coords).unwrap();
    assert_eq!(zmat.len(), 1);

    let (mol, _, coords) = prepare("[H][H]");
    let mut zmat = ZMatrix::new(2);
    set_molecule_dfs_from(&mol, &mut zmat, &coords, 1, 0).unwrap();
    assert_eq!(zmat.len(), 2);
    assert_eq!(zmat.row(0).unwrap().atom_idx, 1);
}
