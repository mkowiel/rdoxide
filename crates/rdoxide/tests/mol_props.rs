//! Molecular weight and formula tests; expected values from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::mol_props::{get_avg_mol_wt, get_exact_mol_wt, get_mol_formula};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn smi(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}

#[test]
fn weights() {
    let m = smi("CC(=O)Nc1ccc(O)cc1");
    close(get_avg_mol_wt(&m, false).unwrap(), 151.165);
    close(get_avg_mol_wt(&m, true).unwrap(), 142.093);
    close(get_exact_mol_wt(&m, false).unwrap(), 151.063328528);
    close(get_exact_mol_wt(&m, true).unwrap(), 141.99290324);
}

#[test]
fn charge_changes_the_exact_weight_by_the_electron_mass() {
    close(
        get_exact_mol_wt(&smi("OC(=O)[O-]"), false).unwrap(),
        60.99311747191,
    );
}

#[test]
fn formulas() {
    let formula = |s: &str, iso: bool, abbrev: bool| get_mol_formula(&smi(s), iso, abbrev).unwrap();
    assert_eq!(formula("C[13CH2]O", false, true), "C2H6O");
    assert_eq!(formula("C[13CH2]O", true, true), "C[13C]H6O");
    assert_eq!(formula("[2H]OC", true, true), "CH3DO");
    assert_eq!(formula("[2H]OC", true, false), "CH3[2H]O");
    assert_eq!(formula("[NH4+].[Cl-]", false, true), "H4ClN");
    assert_eq!(formula("OC(=O)[O-]", false, true), "CHO3-");
    assert_eq!(formula("[Fe+3]", false, true), "Fe+3");
    assert_eq!(formula("BrCCl", false, true), "CH2BrCl");
}
