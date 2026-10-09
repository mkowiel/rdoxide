//! Bond wedging; expectations from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::bond::BondDir;
use rdoxide::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use rdoxide::graphmol::file_parsers::mol_file_writer::{MolWriterParams, mol_to_mol_block};
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_block};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::wedge_bonds::{
    BondWedgingParameters, clear_mol_block_wedging_info, invert_mol_block_wedging_info,
    reapply_mol_block_wedging, wedge_bond, wedge_mol_bonds,
};

fn depicted(smiles: &str) -> RWMol {
    let mut mol = mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    let params = Compute2DCoordParameters {
        canon_orient: true,
        ..Default::default()
    };
    compute_2d_coords(&mut mol, &params).unwrap();
    mol
}

/// `(begin, end, direction)` of the bonds carrying a wedge.
fn wedged(mol: &RWMol) -> Vec<(u32, u32, BondDir)> {
    mol.bonds()
        .iter()
        .filter(|b| b.bond_dir() != BondDir::None)
        .map(|b| (b.begin_atom_idx(), b.end_atom_idx(), b.bond_dir()))
        .collect()
}

#[test]
fn wedges_one_bond_per_center() {
    use BondDir::{BeginDash, BeginWedge};
    for (smiles, expected) in [
        ("C[C@H](F)Cl", vec![(1, 0, BeginWedge)]),
        ("N[C@@H](C)C(=O)O", vec![(1, 2, BeginDash)]),
        (
            "C[C@H]1CC[C@@H](C)CC1",
            vec![(1, 0, BeginDash), (4, 5, BeginDash)],
        ),
        (
            "C1CC[C@@H]2CCCC[C@H]2C1",
            vec![(3, 2, BeginWedge), (8, 7, BeginDash)],
        ),
    ] {
        let mut mol = depicted(smiles);
        wedge_mol_bonds(&mut mol, None, &BondWedgingParameters::default()).unwrap();
        assert_eq!(wedged(&mol), expected, "{smiles}");
    }
}

#[test]
fn can_wedge_a_second_bond_at_degree_four_centers() {
    let two = BondWedgingParameters {
        wedge_two_bonds_if_possible: true,
    };
    let mut mol = depicted("F[C@](Cl)(Br)I");
    wedge_mol_bonds(&mut mol, None, &BondWedgingParameters::default()).unwrap();
    assert_eq!(wedged(&mol), vec![(1, 0, BondDir::BeginDash)]);

    let mut mol = depicted("F[C@](Cl)(Br)I");
    wedge_mol_bonds(&mut mol, None, &two).unwrap();
    assert_eq!(
        wedged(&mol),
        vec![(1, 0, BondDir::BeginDash), (1, 4, BondDir::BeginWedge)]
    );
}

#[test]
fn wedging_needs_a_conformer() {
    let mut mol = mol_from_smiles("C[C@H](F)Cl", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    let err = wedge_mol_bonds(&mut mol, None, &BondWedgingParameters::default())
        .unwrap_err()
        .to_string();
    assert!(err.contains("no conformer available"), "{err}");
}

#[test]
fn wedge_bond_sets_a_direction_from_the_given_atom() {
    let mut mol = depicted("C[C@H](F)Cl");
    let conf = mol.conformer(-1).unwrap().clone();
    wedge_bond(&mut mol, 0, 1, &conf).unwrap();
    assert_eq!(wedged(&mol), vec![(0, 1, BondDir::BeginWedge)]);
    // not from a chiral atom
    let err = wedge_bond(&mut mol, 0, 0, &conf).unwrap_err().to_string();
    assert!(err.contains("Test Assert"), "{err}");
    // non-single bonds are left alone
    let mut alkene = depicted("C/C=C/C");
    let before = wedged(&alkene);
    let conf = alkene.conformer(-1).unwrap().clone();
    wedge_bond(&mut alkene, 1, 1, &conf).unwrap();
    assert_eq!(wedged(&alkene), before);
}

#[test]
fn mol_block_wedging_round_trips() {
    let mut mol = depicted("C[C@H](F)Cl");
    let conf = mol.conformer(-1).unwrap().clone();
    wedge_bond(&mut mol, 0, 1, &conf).unwrap();
    wedge_bond(&mut mol, 2, 1, &conf).unwrap();
    let block = mol_to_mol_block(&mol, &MolWriterParams::default(), -1).unwrap();
    let parse = || {
        mol_from_mol_block(&block, &MolFileParserParams::default())
            .unwrap()
            .unwrap()
    };

    let mut all = parse();
    reapply_mol_block_wedging(&mut all, true, false).unwrap();
    assert_eq!(
        wedged(&all),
        vec![(0, 1, BondDir::BeginWedge), (1, 3, BondDir::BeginWedge)]
    );
    let mut verified = parse();
    reapply_mol_block_wedging(&mut verified, true, true).unwrap();
    assert_eq!(wedged(&verified), vec![(1, 3, BondDir::BeginWedge)]);

    // inverting swaps wedges and dashes, clearing forgets them
    let mut inverted = parse();
    invert_mol_block_wedging_info(&mut inverted);
    reapply_mol_block_wedging(&mut inverted, true, false).unwrap();
    assert_eq!(
        wedged(&inverted),
        vec![(0, 1, BondDir::BeginDash), (1, 3, BondDir::BeginDash)]
    );
    let mut cleared = parse();
    clear_mol_block_wedging_info(&mut cleared);
    reapply_mol_block_wedging(&mut cleared, true, false).unwrap();
    assert!(wedged(&cleared).is_empty());
}
