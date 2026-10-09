//! Reaction SMARTS/SMILES parsing, running and writing; expectations from RDKit 2026.03.6.

use rdoxide::graphmol::chem_reactions::daylight::{
    ReactionSmartsParserParams, reaction_from_smarts, reaction_from_smiles,
};
use rdoxide::graphmol::chem_reactions::mdl::{
    reaction_from_rxn_block, reaction_from_rxn_mol, reaction_to_rxn_block, reaction_to_rxn_mol,
    reaction_to_v3k_rxn_block,
};
use rdoxide::graphmol::chem_reactions::runner::run_reactants;
use rdoxide::graphmol::chem_reactions::writer::{reaction_to_smarts, reaction_to_smiles};
use rdoxide::graphmol::chem_reactions::{ChemicalReaction, get_reacting_atoms};
use rdoxide::graphmol::file_parsers::MolFileParserParams;
use rdoxide::graphmol::romol::ROMol;
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

const AMIDE: &str = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]";

fn mol(smiles: &str) -> ROMol {
    mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn reaction(smarts: &str) -> ChemicalReaction {
    let mut rxn = reaction_from_smarts(smarts, &ReactionSmartsParserParams::new()).unwrap();
    rxn.init_reactant_matchers(true);
    rxn
}

fn canonical_smiles(sets: Vec<Vec<ROMol>>) -> Vec<Vec<String>> {
    sets.into_iter()
        .map(|set| {
            set.into_iter()
                .map(|mut m| {
                    rdoxide::graphmol::mol_ops::sanitize_mol(
                        &mut m,
                        rdoxide::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL,
                    )
                    .unwrap();
                    mol_to_smiles(&mut m, &SmilesWriteParams::default()).unwrap()
                })
                .collect()
        })
        .collect()
}

#[test]
fn amide_formation() {
    let rxn = reaction(AMIDE);
    let products = run_reactants(&rxn, &[mol("CC(=O)O"), mol("NC")], 1000, 1).unwrap();
    assert_eq!(canonical_smiles(products), [["CNC(C)=O"]]);
}

#[test]
fn products_do_not_depend_on_thread_count() {
    let rxn = reaction(AMIDE);
    let reactants = [mol("OC(=O)c1ccccc1C(=O)O"), mol("NCC")];
    let serial = canonical_smiles(run_reactants(&rxn, &reactants, 1000, 1).unwrap());
    let parallel = canonical_smiles(run_reactants(&rxn, &reactants, 1000, 4).unwrap());
    assert_eq!(serial.len(), 2);
    assert_eq!(serial, parallel);
}

#[test]
fn chirality_inversion() {
    let rxn = reaction("[C@:1]([F:2])([Cl:3])([Br:4])>>[C@@:1]([F:2])([Cl:3])([Br:4])");
    let products = run_reactants(&rxn, &[mol("F[C@](Cl)(Br)C")], 1000, 1).unwrap();
    assert_eq!(canonical_smiles(products), [["C[C@](F)(Cl)Br"]]);
}

#[test]
fn reactant_count_is_checked() {
    let rxn = reaction("[C:1][OH]>>[C:1]Cl");
    let err = run_reactants(&rxn, &[mol("CCO"), mol("C")], 1000, 1).unwrap_err();
    assert!(err.to_string().contains("Number of reactants"));
}

#[test]
fn running_requires_initialization() {
    let rxn = reaction_from_smarts(AMIDE, &ReactionSmartsParserParams::new()).unwrap();
    assert!(!rxn.is_initialized());
    assert!(run_reactants(&rxn, &[mol("CC(=O)O"), mol("NC")], 1000, 1).is_err());
}

#[test]
fn smarts_round_trip() {
    let rxn = reaction(AMIDE);
    assert_eq!(
        reaction_to_smarts(&rxn, &SmilesWriteParams::default()).unwrap(),
        AMIDE
    );
    let rxn = reaction("[C:1][OH]>>[C:1]Cl");
    assert_eq!(
        reaction_to_smarts(&rxn, &SmilesWriteParams::default()).unwrap(),
        "[C:1][O&H1]>>[C:1]Cl"
    );
}

#[test]
fn smiles_reaction_and_agents() {
    let rxn =
        reaction_from_smiles("CC(=O)O.N>[Pd]>CC(=O)N", &ReactionSmartsParserParams::new()).unwrap();
    assert_eq!(
        (
            rxn.num_reactant_templates(),
            rxn.num_agent_templates(),
            rxn.num_product_templates()
        ),
        (2, 1, 1)
    );
    assert_eq!(
        reaction_to_smiles(&rxn, &SmilesWriteParams::default()).unwrap(),
        "CC(=O)O.N>[Pd]>CC(N)=O"
    );
}

#[test]
fn malformed_reactions_are_rejected() {
    let p = ReactionSmartsParserParams::new();
    assert!(reaction_from_smarts("CC>", &p).is_err());
    assert!(reaction_from_smarts("C>C>C>C", &p).is_err());
}

#[test]
fn reacting_atoms() {
    let rxn = reaction(AMIDE);
    assert_eq!(
        get_reacting_atoms(&rxn, false).unwrap(),
        [vec![0, 2], vec![0]]
    );
}

fn round_trip_smarts(rxn: &ChemicalReaction) -> String {
    reaction_to_smarts(rxn, &SmilesWriteParams::default()).unwrap()
}

#[test]
fn rxn_blocks_round_trip() {
    let rxn = reaction("[C:1](=[O:2])O.[N:3]>[Pd]>[C:1](=[O:2])[N:3]");
    let expected = "[#6:1](=[#8:2])[#8].[#7:3]>[Pd]>[#6:1](=[#8:2])[#7:3]";
    let params = MolFileParserParams {
        sanitize: false,
        remove_hs: false,
        ..MolFileParserParams::default()
    };
    for block in [
        reaction_to_rxn_block(&rxn, true, false).unwrap(),
        reaction_to_rxn_block(&rxn, false, false).unwrap(),
        reaction_to_v3k_rxn_block(&rxn, true).unwrap(),
    ] {
        let back = reaction_from_rxn_block(&block, &params).unwrap();
        assert!(!back.implicit_properties_flag());
        // agents sit with the reactants unless written separately
        if block.contains("BEGIN AGENT") || block.starts_with("$RXN\n\n      RDKit\n\n  2  1  1") {
            assert_eq!(round_trip_smarts(&back), expected);
        } else {
            assert_eq!(back.num_reactant_templates(), 3);
        }
    }
}

#[test]
fn malformed_rxn_blocks_are_rejected() {
    let params = MolFileParserParams::default();
    for block in ["", "$RXN\n", "x", "$RXN\n\n  RDKit\n\n  1  1\n"] {
        assert!(
            reaction_from_rxn_block(block, &params).is_err(),
            "{block:?}"
        );
    }
}

#[test]
fn rxn_mol_round_trip() {
    let rxn = reaction("[C:1](=[O:2])O.[N:3]>[Pd]>[C:1](=[O:2])[N:3]");
    let back = reaction_from_rxn_mol(&reaction_to_rxn_mol(&rxn).unwrap()).unwrap();
    assert_eq!(
        (
            back.num_reactant_templates(),
            back.num_agent_templates(),
            back.num_product_templates()
        ),
        (2, 1, 1)
    );
    assert_eq!(round_trip_smarts(&back), round_trip_smarts(&rxn));
}

#[test]
fn pickles_round_trip() {
    use rdoxide::graphmol::chem_reactions::pickler::{
        pickle_reaction_default, reaction_from_pickle,
    };
    let mut rxn = reaction("[C:1](=[O:2])O.[N:3]>[Pd]>[C:1](=[O:2])[N:3]");
    rxn.substruct_params_mut().use_chirality = true;
    let back = reaction_from_pickle(&pickle_reaction_default(&rxn).unwrap()).unwrap();
    assert_eq!(round_trip_smarts(&back), round_trip_smarts(&rxn));
    assert!(back.is_initialized() && back.implicit_properties_flag());
    assert!(back.substruct_params().use_chirality);
    assert_eq!(back.num_agent_templates(), 1);
    assert!(reaction_from_pickle(b"not a pickle").is_err());
}

#[test]
fn sanitize_maps_rgroup_dummies() {
    use rdoxide::graphmol::chem_reactions::sanitize::{
        default_rxn_adjust_params, flags::SANITIZE_ALL, sanitize_rxn,
    };
    let mut rxn =
        reaction_from_smarts("[*:1]C>>[*:1]CC", &ReactionSmartsParserParams::new()).unwrap();
    sanitize_rxn(&mut rxn, SANITIZE_ALL, &default_rxn_adjust_params()).unwrap();
    assert_eq!(rxn.num_reactant_templates(), 1);
}

#[test]
fn preprocess_adds_functional_group_queries() {
    use rdoxide::graphmol::chem_reactions::preprocess::preprocess_reaction_with_default_groups;
    use rdoxide::rdgeneral::{RDProps, common_properties::molFileValue};
    let mut rxn = reaction_from_smarts(
        "[C:1]Br.[N:2]>>[C:1][N:2]",
        &ReactionSmartsParserParams::new(),
    )
    .unwrap();
    rxn.reactants_mut()[0]
        .atom_mut(1)
        .dict_mut()
        .set_val(molFileValue, "Halogen.Aromatic");
    let result = preprocess_reaction_with_default_groups(&mut rxn, molFileValue).unwrap();
    assert!(result.valid);
    assert_eq!(
        result.reactant_labels,
        [vec![(1, "halogen.aromatic".to_string())], vec![]]
    );
    rxn.reactants_mut()[0]
        .atom_mut(1)
        .dict_mut()
        .set_val(molFileValue, "boromicacid");
    let err = preprocess_reaction_with_default_groups(&mut rxn, molFileValue).unwrap_err();
    assert!(matches!(err, rdoxide::RDKitError::Key(k) if k == "boromicacid"));
}

#[test]
fn cx_extensions_are_numbered_over_all_templates() {
    use rdoxide::graphmol::chem_reactions::writer::reaction_to_cx_smarts;
    use rdoxide::graphmol::smiles_parse::cx_write::cx_fields;
    // expectations from RDKit 2026.03.6
    let text = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3] |$_R1;;;;;;$,atomProp:1.foo.bar|";
    let rxn = reaction_from_smarts(text, &ReactionSmartsParserParams::new()).unwrap();
    use rdoxide::rdgeneral::RDProps;
    let label = |t: usize, a: u32| rxn.reactants()[t].atom(a).has_prop("atomLabel");
    assert!(label(0, 0));
    assert!(!label(1, 0));
    let written =
        reaction_to_cx_smarts(&rxn, &SmilesWriteParams::default(), cx_fields::ALL).unwrap();
    assert_eq!(written, text);
}

#[test]
fn reaction_fingerprints_match_reference() {
    use rdoxide::datastructs::BitVect;
    use rdoxide::graphmol::chem_reactions::fingerprints::{
        ReactionFingerprintParams, difference_fingerprint, structural_fingerprint,
    };
    // expectations from RDKit 2026.03.6
    let rxn = reaction_from_smarts(AMIDE, &ReactionSmartsParserParams::new()).unwrap();
    let diff =
        difference_fingerprint(&rxn, &ReactionFingerprintParams::default_difference()).unwrap();
    let elements: Vec<(i128, i32)> = diff
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(elements, [(367, 10), (399, -10), (497, 10), (1731, -10)]);
    let structural =
        structural_fingerprint(&rxn, &ReactionFingerprintParams::default_structural()).unwrap();
    assert_eq!(structural.num_bits(), 4096);
    assert_eq!(structural.on_bits()[..5], [71, 287, 292, 293, 319]);
}
