//! Random flip sampling in the depictor and `compute_2d_coords_mimic_dist_mat`;
//! the reference structures are the ones of RDKit 2026.09.1's `testDepictor.py`.

use std::path::Path;

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::depictor::{
    Compute2DCoordParameters, compute_2d_coords, compute_2d_coords_mimic_dist_mat,
};
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_file};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::random::MinStdRand;

const DATA: &str = "../../upstream/rdkit/Code/GraphMol/Depictor/test_data";

fn have_data() -> bool {
    Path::new(DATA).is_dir()
}

fn load(name: &str) -> RWMol {
    mol_from_mol_file(&format!("{DATA}/{name}"), &MolFileParserParams::default())
        .unwrap()
        .unwrap()
}

fn coords(mol: &RWMol) -> Vec<(f64, f64)> {
    mol.conformer(-1)
        .unwrap()
        .positions()
        .iter()
        .map(|p| (p.x, p.y))
        .collect()
}

fn assert_same_coords(mol: &RWMol, reference: &RWMol) {
    let (a, b) = (coords(mol), coords(reference));
    assert_eq!(a.len(), b.len());
    for (i, (p, q)) in a.iter().zip(&b).enumerate() {
        assert!(
            (p.0 - q.0).abs() < 1e-3 && (p.1 - q.1).abs() < 1e-3,
            "atom {i}: {p:?} != {q:?}"
        );
    }
}

fn dist_mat(mol: &RWMol) -> Vec<f64> {
    let pos = mol.conformer(-1).unwrap().positions().to_vec();
    let n = pos.len();
    let mut res = vec![0.0; n * (n - 1) / 2];
    for i in 1..n {
        for j in 0..i {
            res[i * (i - 1) / 2 + j] = (pos[j] - pos[i]).length();
        }
    }
    res
}

#[test]
fn sampling_spread_matches_upstream() {
    if !have_data() {
        return;
    }
    let mut mol = load("7UPJ_xtal.mol");
    let params = Compute2DCoordParameters {
        canon_orient: false,
        n_flips_per_sample: 3,
        n_samples: 100,
        sample_seed: 100,
        permute_deg4_nodes: true,
        ..Default::default()
    };
    compute_2d_coords(&mut mol, &params).unwrap();
    assert_same_coords(&mol, &load("7UPJ_spread.mol"));
}

#[test]
fn mimic_distmat_matches_upstream() {
    if !have_data() {
        return;
    }
    let mut mol = load("7UPJ_xtal.mol");
    let dmat = dist_mat(&mol);
    compute_2d_coords_mimic_dist_mat(&mut mol, &dmat, false, true, 0.001, 3, 100, 100, true)
        .unwrap();
    assert_same_coords(&mol, &load("7UPJ_mimic3D_1.mol"));
    compute_2d_coords_mimic_dist_mat(&mut mol, &dmat, false, true, 0.003, 3, 100, 100, true)
        .unwrap();
    assert_same_coords(&mol, &load("7UPJ_mimic3D_2.mol"));
}

#[test]
fn seeded_sampling_is_reproducible() {
    let smi = "CCCCC(CC)C(=O)NCC1CCCCC1";
    let depict = || {
        let mut mol = mol_from_smiles(smi, &SmilesParserParams::default())
            .unwrap()
            .unwrap();
        let params = Compute2DCoordParameters {
            n_flips_per_sample: 2,
            n_samples: 20,
            sample_seed: 7,
            ..Default::default()
        };
        compute_2d_coords(&mut mol, &params).unwrap();
        coords(&mol)
    };
    assert_eq!(depict(), depict());
}

#[test]
fn minstd_uniform_int_stays_in_range() {
    let mut rng = MinStdRand::new(42);
    assert_eq!(rng.uniform_int(5, 5), 5);
    let mut seen = [false; 7];
    for _ in 0..500 {
        let v = rng.uniform_int(0, 6);
        seen[v as usize] = true;
    }
    assert!(seen.iter().all(|&s| s));
}
