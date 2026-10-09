//! Reaction SMARTS/SMILES parsing, running and writing; expectations from RDKit 2026.03.6.

use rdoxide::graphmol::chem_reactions::daylight::{
    ReactionSmartsParserParams, reaction_from_smarts, reaction_from_smiles,
};
use rdoxide::graphmol::chem_reactions::runner::run_reactants;
use rdoxide::graphmol::chem_reactions::writer::{reaction_to_smarts, reaction_to_smiles};
use rdoxide::graphmol::chem_reactions::{ChemicalReaction, get_reacting_atoms};
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
