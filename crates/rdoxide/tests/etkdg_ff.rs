//! The ETKDG force fields refine a distance-geometry embedding.

use rdoxide::distgeom::{
    BoundsMatrix, ChiralSet, add_torsion_terms, compute_random_coords, construct_3d_force_field,
    construct_3d_improper_force_field, construct_all_in_one_force_field, construct_force_field,
    construct_plain_3d_force_field, triangle_smooth_bounds,
};
use rdoxide::forcefield::ForceField;
use rdoxide::graphmol::dist_geom_helpers::{
    EmbedParameters, TopolBoundsOptions, TopolBoundsOutput, collect_bonds_and_angles,
    init_bounds_mat, set_topol_bounds,
};
use rdoxide::graphmol::forcefield_helpers::crystal_ff::{
    CrystalFFDetails, TorsionOptions, get_experimental_torsions, populate_ref_table,
};
use rdoxide::graphmol::matrices::get_distance_mat;
use rdoxide::graphmol::mol_ops::add_hs::{AddHsParameters, add_hs};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::{ROMol, RWMol};
use rdoxide::rdgeneral::random::MinStdRand;

struct Setup {
    mol: RWMol,
    mmat: BoundsMatrix,
    details: CrystalFFDetails,
}

fn setup(smiles: &str, version: u32) -> Setup {
    let mut mol = mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    add_hs(&mut mol, &AddHsParameters::default(), None).unwrap();
    let n = mol.num_atoms() as usize;
    let mut mmat = BoundsMatrix::new(n);
    init_bounds_mat(&mut mmat, 0.0, 1000.0).unwrap();
    let mut paths = Vec::new();
    set_topol_bounds(
        &mol,
        &mut mmat,
        &EmbedParameters::default(),
        &TopolBoundsOptions::default(),
        TopolBoundsOutput {
            paths14: Some(&mut paths),
            internal_coords: None,
        },
    )
    .unwrap();
    assert!(triangle_smooth_bounds(&mut mmat, 0.0));
    let (bonds, angles) = collect_bonds_and_angles(&mol);
    let mut details = CrystalFFDetails {
        bonds,
        angles,
        atom_nums: (0..n as u32)
            .map(|i| mol.atom(i).atomic_num() as i32)
            .collect(),
        bounds_mat_force_scaling: 1.0,
        dist_mat: get_distance_mat(&mol),
        path14_configs: paths,
        ..Default::default()
    };
    let options = TorsionOptions {
        use_exp_torsions: true,
        use_basic_knowledge: true,
        version,
        ..Default::default()
    };
    get_experimental_torsions(&mol, &mut details, &options).unwrap();
    if version == 4 {
        populate_ref_table(&mut details).unwrap();
    }
    Setup { mol, mmat, details }
}

/// Random coordinates relaxed into the bounds, the starting point for ETKDG.
fn dg_coordinates(mmat: &BoundsMatrix, seed: u32) -> Vec<f64> {
    let n = mmat.num_rows();
    let mut rng = MinStdRand::new(seed);
    let mut pos = vec![0.0; 3 * n];
    compute_random_coords(&mut pos, 10.0, &mut rng).unwrap();
    let mut field = construct_force_field(
        mmat,
        &pos,
        3,
        &Vec::<ChiralSet>::new(),
        1.0,
        0.0,
        0.0,
        None,
        5.0,
        None,
    )
    .unwrap();
    field.minimize(400, 1e-3, 1e-6).unwrap();
    field.positions().to_vec()
}

fn dist(pos: &[f64], i: usize, j: usize) -> f64 {
    (0..3)
        .map(|k| (pos[3 * i + k] - pos[3 * j + k]).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn refine(field: &mut ForceField) -> Vec<f64> {
    for _ in 0..10 {
        if field.minimize(400, 1e-3, 1e-6).unwrap() == 0 {
            break;
        }
    }
    field.positions().to_vec()
}

fn check_geometry(mol: &ROMol, mmat: &BoundsMatrix, pos: &[f64]) {
    for b in mol.bonds() {
        let (i, j) = (b.begin_atom_idx() as usize, b.end_atom_idx() as usize);
        let d = dist(pos, i, j);
        assert!(
            d > mmat.lower_bound(i, j) - 0.05 && d < mmat.upper_bound(i, j) + 0.05,
            "bond {i}-{j}: {d}"
        );
    }
}

#[test]
fn etkdg_field_refines_cosine_torsions() {
    let s = setup("CC(=O)NC", 2);
    let start = dg_coordinates(&s.mmat, 42);
    let mut field = construct_3d_force_field(&s.mmat, &start, &s.details, None).unwrap();
    let before = field.calc_energy();
    let pos = refine(&mut field);
    assert!(field.calc_energy() < before);
    check_geometry(&s.mol, &s.mmat, &pos);
    // the forced amide: the carbonyl O (2) and the H on N (3) are far apart
    let h_on_n = s
        .mol
        .atom_neighbors(3)
        .find(|&n| s.mol.atom(n).atomic_num() == 1)
        .unwrap() as usize;
    assert!(dist(&pos, 2, h_on_n) > 3.0, "{}", dist(&pos, 2, h_on_n));
}

#[test]
fn gaussian_field_and_plain_field_build() {
    let s = setup("CC(=O)NC", 4);
    let start = dg_coordinates(&s.mmat, 7);
    let mut field = construct_3d_force_field(&s.mmat, &start, &s.details, None).unwrap();
    let before = field.calc_energy();
    let pos = refine(&mut field);
    assert!(field.calc_energy() <= before);
    check_geometry(&s.mol, &s.mmat, &pos);

    let plain = construct_plain_3d_force_field(&s.mmat, &start, &s.details).unwrap();
    assert!(plain.calc_energy().is_finite());
    assert!(plain.contribs().len() < field.contribs().len());
}

#[test]
fn improper_and_all_in_one_fields() {
    let s = setup("CC(=O)NC", 2);
    let start = dg_coordinates(&s.mmat, 3);
    let mut improper = construct_3d_improper_force_field(
        &s.mmat,
        &start,
        &s.details.improper_atoms,
        &s.details.angles,
    )
    .unwrap();
    let before = improper.calc_energy();
    refine(&mut improper);
    assert!(improper.calc_energy() <= before);

    let cpci = std::collections::BTreeMap::from([((0u32, 3u32), 0.5)]);
    let mut aio = construct_all_in_one_force_field(
        &s.mmat,
        &start,
        3,
        &s.details,
        &[],
        Some(&cpci),
        None,
        None,
    )
    .unwrap();
    let n_before = aio.contribs().len();
    add_torsion_terms(&mut aio, &s.details, true, true).unwrap();
    assert!(aio.contribs().len() > n_before);
    assert!(aio.calc_energy().is_finite());
}

#[test]
fn mismatched_inputs_are_rejected() {
    let s = setup("CC", 2);
    assert!(construct_3d_force_field(&s.mmat, &[0.0; 3], &s.details, None).is_err());
}
