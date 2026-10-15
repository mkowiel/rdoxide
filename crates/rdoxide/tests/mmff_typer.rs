//! MMFF atom types and charges of small molecules (the types follow the
//! MMFF94 symbolic type definitions).

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::forcefield_helpers::mmff::MmffMolProperties;
use rdoxide::graphmol::mol_ops::add_hs::{AddHsParameters, add_hs};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn typed(smiles: &str) -> (RWMol, MmffMolProperties) {
    let mut mol = mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    add_hs(&mut mol, &AddHsParameters::default(), None).unwrap();
    let props = MmffMolProperties::new(&mut mol, "MMFF94").unwrap();
    (mol, props)
}

fn types(smiles: &str) -> Vec<u8> {
    let (mol, props) = typed(smiles);
    assert!(props.is_valid(), "{smiles}");
    (0..mol.num_atoms()).map(|i| props.atom_type(i)).collect()
}

#[test]
fn alcohol_and_water() {
    assert_eq!(types("CCO"), [1, 1, 6, 5, 5, 5, 5, 5, 21]);
    assert_eq!(types("O"), [70, 31, 31]);
}

#[test]
fn carbonyl_compounds() {
    // CR, C=O, O=C, OC=O; HC x3, HOCO
    assert_eq!(types("CC(=O)O"), [1, 3, 7, 6, 5, 5, 5, 24]);
    // acetate: CO2M, O2CM x2
    assert_eq!(types("CC(=O)[O-]"), [1, 41, 32, 32, 5, 5, 5]);
}

#[test]
fn aromatic_rings() {
    assert_eq!(
        types("c1ccccc1"),
        [37, 37, 37, 37, 37, 37, 5, 5, 5, 5, 5, 5]
    );
    // NPYD
    assert_eq!(types("c1ccncc1")[3], 38);
    // pyrrole: C5B x2, C5A x2 around NPYL, and HPYL on the nitrogen
    let t = types("c1cc[nH]c1");
    assert_eq!(&t[..5], &[64, 64, 63, 39, 63][..]);
    assert_eq!(t[8], 23);
}

#[test]
fn charges_sum_to_the_net_charge() {
    for (smiles, net) in [
        ("CCO", 0.0),
        ("CC(=O)[O-]", -1.0),
        ("C[NH3+]", 1.0),
        ("c1ccc(cc1)C(=O)NC", 0.0),
        ("CS(=O)(=O)N", 0.0),
    ] {
        let (mol, props) = typed(smiles);
        assert!(props.is_valid(), "{smiles}");
        let sum: f64 = (0..mol.num_atoms()).map(|i| props.partial_charge(i)).sum();
        assert!((sum - net).abs() < 1e-9, "{smiles}: {sum}");
    }
}

#[test]
fn formal_charges() {
    let (_, props) = typed("CC(=O)[O-]");
    assert_eq!(props.formal_charge(2), -0.5);
    assert_eq!(props.formal_charge(3), -0.5);
    let (_, props) = typed("C[NH3+]");
    assert_eq!(props.atom_type(1), 34);
    assert_eq!(props.formal_charge(1), 1.0);
}

#[test]
fn parameters_of_ethanol() {
    let (mol, props) = typed("CCO");
    let b = props.bond_stretch_params(&mol, 0, 1).unwrap().unwrap();
    assert_eq!((b.bond_type, b.params.kb, b.params.r0), (0, 4.258, 1.508));
    let a = props.angle_bend_params(&mol, 0, 1, 2).unwrap().unwrap();
    assert_eq!(a.angle_type, 0);
    assert!(a.params.ka > 0.0);
    assert!(props.stretch_bend_params(&mol, 0, 1, 2).unwrap().is_some());
    assert!(props.torsion_params(&mol, 3, 0, 1, 2).unwrap().is_some());
    // out-of-plane terms need a planar center
    assert!(props.oop_bend_params(&mol, 0, 1, 2, 3).is_none());
}

mod force_field {
    use std::path::PathBuf;

    use rdoxide::graphmol::file_parsers::mol_file_parser::{
        MolFileParserParams, mol_from_mol_block,
    };
    use rdoxide::graphmol::forcefield_helpers::mmff::construct_force_field;
    use rdoxide::graphmol::mol_ops::{sanitize_flags::SANITIZE_ALL, sanitize_mol};

    fn minimized_energy_drop(name: &str, max_its: u32) -> (i32, f64, f64) {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/data/molfiles")
            .join(name);
        let block = std::fs::read_to_string(path).unwrap();
        let params = MolFileParserParams {
            remove_hs: false,
            ..MolFileParserParams::default()
        };
        let mut mol = mol_from_mol_block(&block, &params).unwrap().unwrap();
        sanitize_mol(&mut mol, SANITIZE_ALL).unwrap();
        let mut field = construct_force_field(&mut mol, 100.0, -1, true).unwrap();
        let before = field.calc_energy();
        let needs_more = field.minimize(max_its, 1e-6, 1e-3).unwrap();
        let after = field.calc_energy();
        (needs_more, before, after)
    }

    #[test]
    fn minimization_converges() {
        // the molecules of upstream's testMMFFHelpers
        for (name, its) in [
            ("mmff_small1.mol", 200),
            ("mmff_small2.mol", 100),
            ("mmff_benzene.mol", 200),
            ("mmff_toluene.mol", 200),
            ("mmff_complex1.mol", 200),
        ] {
            let (needs_more, before, after) = minimized_energy_drop(name, its);
            assert_eq!(needs_more, 0, "{name}");
            assert!(after <= before + 1e-9, "{name}: {before} -> {after}");
        }
    }
}
