"""Molecular weights, formula and atom counts (expectations from RDKit 2026.09.1)."""

import pytest

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D

# smiles: (MW, MW heavy, exact MW, exact MW heavy, formula, formula with isotopes,
#          formula with isotopes and no D/T abbreviation, heavy atoms, all atoms)
CASES = {
    "CC(=O)Nc1ccc(O)cc1": (151.165, 142.093, 151.063328528, 141.99290324, "C8H9NO2", "C8H9NO2",
                           "C8H9NO2", 11, 20),
    "C[13CH2]O": (47.06135484, 41.01335484, 47.045219652, 40.99826946, "C2H6O", "C[13C]H6O",
                  "C[13C]H6O", 3, 9),
    "[NH4+].[Cl-]": (53.492, 49.46, 53.003226808, 48.97192668, "H4ClN", "H4ClN", "H4ClN", 2, 6),
    "[2H]OC": (33.048101778, 28.01, 33.032491494, 27.99491462, "CH4O", "CH3DO", "CH3[2H]O", 2, 6),
    "OC(=O)[O-]": (61.016, 60.008, 60.99311747191, 59.98529243991, "CHO3-", "CHO3-", "CHO3-", 4, 5),
    "[Fe+3]": (55.845, 55.845, 55.93329176027, 55.93329176027, "Fe+3", "Fe+3", "Fe+3", 1, 1),
    "*C": (15.035, 12.011, 15.023475096, 12.0, "CH3*", "CH3*", "CH3*", 1, 5),
}


@pytest.mark.parametrize("smiles", CASES)
def test_weights_formula_and_counts(smiles):
    amw, amw_heavy, exact, exact_heavy, formula, formula_iso, formula_iso_full, heavy, atoms = CASES[
        smiles]
    m = Chem.MolFromSmiles(smiles)
    assert D._CalcMolWt(m) == pytest.approx(amw, abs=1e-9)
    assert D._CalcMolWt(m, onlyHeavy=True) == pytest.approx(amw_heavy, abs=1e-9)
    assert D.CalcExactMolWt(m) == pytest.approx(exact, abs=1e-9)
    assert D.CalcExactMolWt(m, True) == pytest.approx(exact_heavy, abs=1e-9)
    assert D.CalcMolFormula(m) == formula
    assert D.CalcMolFormula(m, True) == formula_iso
    assert D.CalcMolFormula(m, separateIsotopes=True, abbreviateHIsotopes=False) == formula_iso_full
    assert D.CalcNumHeavyAtoms(m) == heavy
    assert D.CalcNumAtoms(m) == atoms


def test_hill_order_puts_carbon_and_hydrogen_first():
    assert D.CalcMolFormula(Chem.MolFromSmiles("BrCCl")) == "CH2BrCl"
    assert D.CalcMolFormula(Chem.MolFromSmiles("OS(=O)(=O)O")) == "H2O4S"


def test_bool_arguments_take_truthy_values():
    m = Chem.MolFromSmiles("CCO")
    assert D._CalcMolWt(m, 1) == D._CalcMolWt(m, True)
    assert D.CalcExactMolWt(m, onlyHeavy=0) == D.CalcExactMolWt(m, onlyHeavy=False)


def test_versions():
    assert D._CalcMolWt_version == "1.0.0"
    assert D._CalcExactMolWt_version == "1.0.0"
    assert D._CalcMolFormula_version == "1.3.0"
    assert D._CalcNumHeavyAtoms_version == "1.0.0"
    assert D._CalcNumAtoms_version == "1.0.0"
