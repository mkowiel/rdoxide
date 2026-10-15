//! Experimental torsion preferences and knowledge terms of small molecules.

use rdoxide::distgeom::BoundsMatrix;
use rdoxide::graphmol::dist_geom_helpers::{
    EmbedParameters, TopolBoundsOptions, TopolBoundsOutput, init_bounds_mat, set_topol_bounds,
};
use rdoxide::graphmol::forcefield_helpers::crystal_ff::{
    AMIDE_TRANS_IDX, CrystalFFDetails, EtkdgForceConsts, TorsionAngleRef, TorsionOptions,
    TorsionParamKind, TorsionParams, get_experimental_torsions,
    get_experimental_torsions_with_bonds, populate_ref_table,
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

/// The details with the 1-4 paths of the molecule, as the embedder has them.
fn details_for(mol: &ROMol) -> CrystalFFDetails {
    let mut mmat = BoundsMatrix::new(mol.num_atoms() as usize);
    init_bounds_mat(&mut mmat, 0.0, 1000.0).unwrap();
    let mut paths = Vec::new();
    set_topol_bounds(
        mol,
        &mut mmat,
        &EmbedParameters::default(),
        &TopolBoundsOptions::default(),
        TopolBoundsOutput {
            paths14: Some(&mut paths),
            internal_coords: None,
        },
    )
    .unwrap();
    CrystalFFDetails {
        path14_configs: paths,
        ..Default::default()
    }
}

#[test]
fn experimental_torsions_match_the_rotatable_bonds() {
    // ethyl benzoate: the ester and aryl-C(=O) torsions have preferences
    let mol = mol_with_hs("CCOC(=O)c1ccccc1");
    let mut details = CrystalFFDetails::default();
    let options = TorsionOptions {
        use_exp_torsions: true,
        ..Default::default()
    };
    let bonds = get_experimental_torsions_with_bonds(&mol, &mut details, &options).unwrap();
    assert_eq!(details.torsion_param_kind, TorsionParamKind::Cosine);
    assert!(!details.exp_torsion_atoms.is_empty());
    assert_eq!(
        details.exp_torsion_atoms.len(),
        details.exp_torsion_angles.len()
    );
    assert_eq!(details.exp_torsion_atoms.len(), details.torsion_idx.len());
    assert_eq!(bonds.len(), details.exp_torsion_atoms.len());
    // every central bond is used once, and exists
    let mut central: Vec<_> = bonds.iter().map(|b| b.0).collect();
    central.sort_unstable();
    central.dedup();
    assert_eq!(central.len(), bonds.len());
    for (bid, atoms, angle) in &bonds {
        let bond = mol.bond(*bid);
        let ends = [bond.begin_atom_idx(), bond.end_atom_idx()];
        assert!(ends.contains(&atoms[1]) && ends.contains(&atoms[2]));
        assert!(matches!(angle, TorsionAngleRef::Cosine(_)));
    }
    for params in &details.exp_torsion_angles {
        let TorsionParams::Cosine { signs, v } = params else {
            panic!("cosine parameters expected");
        };
        assert_eq!((signs.len(), v.len()), (6, 6));
    }
}

#[test]
fn the_scaling_multiplies_the_cosine_constants() {
    let mol = mol_with_hs("CCOC(=O)c1ccccc1");
    let options = TorsionOptions {
        use_exp_torsions: true,
        ..Default::default()
    };
    let mut plain = CrystalFFDetails::default();
    get_experimental_torsions(&mol, &mut plain, &options).unwrap();
    let mut scaled = CrystalFFDetails {
        force_consts: EtkdgForceConsts {
            et_term_scaling: 0.5,
            ..EtkdgForceConsts::SEQ_COSINE
        },
        ..Default::default()
    };
    get_experimental_torsions(&mol, &mut scaled, &options).unwrap();
    let (TorsionParams::Cosine { v: a, .. }, TorsionParams::Cosine { v: b, .. }) =
        (&plain.exp_torsion_angles[0], &scaled.exp_torsion_angles[0])
    else {
        panic!("cosine parameters expected");
    };
    for (x, y) in a.iter().zip(b) {
        assert!((x * 0.5 - y).abs() < 1e-12);
    }
}

#[test]
fn basic_knowledge_adds_amides_inversions_and_flat_rings() {
    let mol = mol_with_hs("CC(=O)Nc1ccccc1");
    let mut details = details_for(&mol);
    let options = TorsionOptions {
        use_basic_knowledge: true,
        ..Default::default()
    };
    get_experimental_torsions(&mol, &mut details, &options).unwrap();
    // the benzene ring has six flat torsions, the forced amide adds more
    let arom = details
        .torsion_idx
        .iter()
        .filter(|&&i| i == rdoxide::graphmol::forcefield_helpers::crystal_ff::AROM_TORSION_IDX)
        .count();
    assert_eq!(arom, 6);
    assert!(details.torsion_idx.contains(&AMIDE_TRANS_IDX) || details.torsion_idx.contains(&2001));
    // inversion centres: the carbonyl C, the amide N and the 6 ring carbons
    // that have three heavy neighbours or hydrogens
    assert!(details.improper_atoms.len() >= 8);
    for imp in &details.improper_atoms {
        assert_eq!(mol.atom_degree(imp.atoms[1]), 3);
    }
    assert!(details.improper_atoms.iter().any(|i| i.bound_to_sp2_o));
}

#[test]
fn gaussian_torsions_need_version_4() {
    let mol = mol_with_hs("CCOC(=O)c1ccccc1");
    let mut details = CrystalFFDetails::default();
    let bad = TorsionOptions {
        use_exp_torsions: true,
        version: 3,
        ..Default::default()
    };
    assert!(get_experimental_torsions(&mol, &mut details, &bad).is_err());

    let options = TorsionOptions {
        use_exp_torsions: true,
        use_basic_knowledge: true,
        version: 4,
        ..Default::default()
    };
    let mut details = details_for(&mol);
    get_experimental_torsions(&mol, &mut details, &options).unwrap();
    assert_eq!(details.torsion_param_kind, TorsionParamKind::Gaussian);
    populate_ref_table(&mut details).unwrap();
    assert!(!details.phi_to_energy.is_empty());
    assert_eq!(details.phi_to_energy.len(), details.phi_to_grad.len());
    // the indices are ranks into the tables
    assert!(
        details
            .torsion_idx
            .iter()
            .all(|&i| i < details.phi_to_energy.len())
    );
    assert!(details.phi_to_energy.iter().all(|t| t.len() == 180));
}

#[test]
fn empty_molecule_is_rejected() {
    let mut details = CrystalFFDetails::default();
    assert!(
        get_experimental_torsions(&ROMol::new(), &mut details, &TorsionOptions::default()).is_err()
    );
}
