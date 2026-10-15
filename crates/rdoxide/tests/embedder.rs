//! Distance-geometry embedding of small molecules.

use rdoxide::geometry::Point3D;
use rdoxide::graphmol::dist_geom_helpers::{
    EmbedParameters, InitialEmbeddingMode, embed_molecule, embed_multiple_confs,
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

fn params(base: EmbedParameters, seed: i32) -> EmbedParameters {
    EmbedParameters {
        random_seed: seed,
        ..base
    }
}

fn positions(mol: &ROMol, conf: i32) -> Vec<Point3D> {
    mol.conformer(conf).unwrap().positions().to_vec()
}

fn bond_lengths_are_sane(mol: &ROMol, pos: &[Point3D]) {
    for b in mol.bonds() {
        let d = (pos[b.begin_atom_idx() as usize] - pos[b.end_atom_idx() as usize]).length();
        assert!(d > 0.9 && d < 1.8, "bond {}: {d}", b.idx());
    }
}

/// signed volume of the neighbors of an atom
fn volume(pos: &[Point3D], center: usize, a: usize, b: usize, c: usize) -> f64 {
    let (v1, v2, v3) = (
        pos[a] - pos[center],
        pos[b] - pos[center],
        pos[c] - pos[center],
    );
    v1.dot_product(&v2.cross_product(&v3))
}

#[test]
fn every_parameter_set_embeds_a_simple_molecule() {
    let presets = [
        ("dg", EmbedParameters::dg()),
        ("kdg", EmbedParameters::kdg()),
        ("etdg", EmbedParameters::etdg()),
        ("etkdg", EmbedParameters::etkdg()),
        ("etkdg_v2", EmbedParameters::etkdg_v2()),
        ("etkdg_v3", EmbedParameters::etkdg_v3()),
        ("etkdg_v4", EmbedParameters::etkdg_v4()),
    ];
    for (name, preset) in presets {
        let mut mol = mol_with_hs("CC(=O)NCC(O)C");
        let id = embed_molecule(&mut mol, &mut params(preset, 42)).unwrap();
        assert_eq!(id, 0, "{name}");
        bond_lengths_are_sane(&mol, &positions(&mol, id));
    }
}

#[test]
fn embedding_is_reproducible_and_independent_of_the_thread_count() {
    let run = |threads: i32| {
        let mut mol = mol_with_hs("CCC(C)NC(=O)c1ccccc1");
        let mut p = params(EmbedParameters::etkdg_v3(), 7);
        p.num_threads = threads;
        let ids = embed_multiple_confs(&mut mol, 6, &mut p).unwrap();
        assert_eq!(ids, (0..6).collect::<Vec<_>>());
        ids.iter().map(|&i| positions(&mol, i)).collect::<Vec<_>>()
    };
    let single = run(1);
    assert_eq!(single, run(1));
    assert_eq!(single, run(3));
    // different seeds give different conformers
    assert_ne!(single[0], single[1]);
}

#[test]
fn enantiomers_embed_as_mirror_images() {
    let sign = |smiles: &str| {
        let mut mol = mol_with_hs(smiles);
        let id = embed_molecule(&mut mol, &mut params(EmbedParameters::etkdg_v3(), 11)).unwrap();
        assert_eq!(id, 0);
        let pos = positions(&mol, id);
        // the centre is atom 1; neighbors N0, C2, C3
        volume(&pos, 1, 0, 2, 3)
    };
    let l = sign("N[C@@H](C)C(=O)O");
    let d = sign("N[C@H](C)C(=O)O");
    assert!(l * d < 0.0, "{l} {d}");
    // and with the all-in-one embedding
    let sign_v4 = |smiles: &str| {
        let mut mol = mol_with_hs(smiles);
        let id = embed_molecule(&mut mol, &mut params(EmbedParameters::etkdg_v4(), 11)).unwrap();
        assert_eq!(id, 0);
        volume(&positions(&mol, id), 1, 0, 2, 3)
    };
    assert!(sign_v4("N[C@@H](C)C(=O)O") * sign_v4("N[C@H](C)C(=O)O") < 0.0);
}

#[test]
fn double_bond_stereo_is_respected() {
    let dist_h = |smiles: &str| {
        let mut mol = mol_with_hs(smiles);
        let id = embed_molecule(&mut mol, &mut params(EmbedParameters::etkdg_v3(), 5)).unwrap();
        assert_eq!(id, 0);
        let pos = positions(&mol, id);
        // the terminal methyl carbons C0 and C3
        (pos[0] - pos[3]).length()
    };
    let trans = dist_h("C/C=C/C");
    let cis = dist_h("C/C=C\\C");
    assert!(trans > cis + 0.5, "{trans} {cis}");
}

#[test]
fn fragments_are_embedded_separately() {
    let mut mol = mol_with_hs("CC.O");
    let ids =
        embed_multiple_confs(&mut mol, 2, &mut params(EmbedParameters::etkdg_v3(), 3)).unwrap();
    assert_eq!(ids.len(), 2);
    let pos = positions(&mol, 0);
    bond_lengths_are_sane(&mol, &pos);
}

#[test]
fn pruning_removes_similar_conformers() {
    let mut mol = mol_with_hs("CCO");
    let mut p = params(EmbedParameters::etkdg_v3(), 1);
    p.prune_rms_thresh = 5.0;
    let ids = embed_multiple_confs(&mut mol, 8, &mut p).unwrap();
    // everything is within 5 A of the first conformer
    assert_eq!(ids, vec![0]);
}

#[test]
fn internal_coordinate_and_random_coordinate_modes_work() {
    for mode in [
        InitialEmbeddingMode::InternalCoordinateEmbedding,
        InitialEmbeddingMode::RandomCoordinateEmbedding,
    ] {
        let mut mol = mol_with_hs("CC(C)C(=O)O");
        let mut p = params(EmbedParameters::etkdg_v3(), 21);
        p.initial_embedding_mode = mode;
        p.max_iterations = 200;
        let id = embed_molecule(&mut mol, &mut p).unwrap();
        assert_eq!(id, 0, "{mode:?}");
        bond_lengths_are_sane(&mol, &positions(&mol, id));
    }
}

#[test]
fn failures_can_be_tracked_and_bad_input_is_rejected() {
    let mut mol = mol_with_hs("CCO");
    let mut p = params(EmbedParameters::dg(), 4);
    p.track_failures = true;
    embed_molecule(&mut mol, &mut p).unwrap();
    assert_eq!(p.failures.len(), 15);

    let mut p = EmbedParameters {
        et_version: 3,
        ..EmbedParameters::dg()
    };
    assert!(embed_molecule(&mut mol_with_hs("CC"), &mut p).is_err());
    assert!(embed_molecule(&mut ROMol::new(), &mut EmbedParameters::dg()).is_err());
}

#[test]
fn clear_confs_false_keeps_existing_conformers() {
    let mut mol = mol_with_hs("CCO");
    embed_molecule(&mut mol, &mut params(EmbedParameters::dg(), 2)).unwrap();
    let mut p = params(EmbedParameters::dg(), 3);
    p.clear_confs = false;
    let id = embed_molecule(&mut mol, &mut p).unwrap();
    assert_eq!(id, 1);
    assert_eq!(mol.num_conformers(), 2);
}

#[test]
fn rings_macrocycles_and_larger_molecules_embed() {
    let cases = [
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "C1CC2CCC1C2",
        "C1CCCCCCCCCCC1",
        "O=C1NC(=O)c2ccccc12",
        "C[C@H]1CC[C@@H](C)CC1",
        "CN1CCC[C@H]1c1cccnc1",
    ];
    for smiles in cases {
        for preset in [EmbedParameters::etkdg_v3(), EmbedParameters::etkdg_v4()] {
            let mut mol = mol_with_hs(smiles);
            let mut p = params(preset, 17);
            p.max_iterations = 500;
            let id = embed_molecule(&mut mol, &mut p).unwrap();
            assert_eq!(id, 0, "{smiles} v{}", p.et_version);
            bond_lengths_are_sane(&mol, &positions(&mol, id));
        }
    }
}
