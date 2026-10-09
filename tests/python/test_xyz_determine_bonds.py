"""XYZ I/O and rdDetermineBonds; expected values come from reference RDKit."""

import pytest
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds as DB

ETHANOL = """9
ethanol
C 0.0 0.0 0.0
C 1.52 0.0 0.0
O 2.0 1.35 0.0
H -0.4 1.0 0.0
H -0.4 -0.5 0.9
H -0.4 -0.5 -0.9
H 1.9 -0.5 0.9
H 1.9 -0.5 -0.9
H 2.95 1.35 0.0
"""

ACETIC = """8
acetic acid
C 0.0 0.0 0.0
C 1.5 0.0 0.0
O 2.1 1.05 0.0
O 2.1 -1.2 0.0
H -0.4 1.0 0.0
H -0.4 -0.5 0.9
H -0.4 -0.5 -0.9
H 3.05 -1.2 0.0
"""


def test_xyz_roundtrip():
    m = Chem.MolFromXYZBlock(ETHANOL)
    assert m.GetNumAtoms() == 9 and m.GetNumBonds() == 0
    assert m.GetConformer().GetAtomPosition(2).x == pytest.approx(2.0)
    out = Chem.MolToXYZBlock(m)
    lines = out.splitlines()
    assert lines[0] == "9" and lines[2].split()[0] == "C"
    assert Chem.MolFromXYZBlock(out).GetNumAtoms() == 9


def test_xyz_errors():
    assert Chem.MolFromXYZBlock("3\nx\nC 0 0 0\n") is None
    assert Chem.MolFromXYZFile("/no/such/file.xyz") is None
    assert Chem.MolToXYZBlock(Chem.MolFromSmiles("CC")) == ""


def test_determine_connectivity_and_orders():
    m = Chem.MolFromXYZBlock(ACETIC)
    DB.DetermineConnectivity(m)
    assert m.GetNumBonds() == 7
    DB.DetermineBondOrders(m, charge=0)
    assert Chem.MolToSmiles(Chem.RemoveHs(m)) == "CC(=O)O"


def test_determine_bonds_one_step():
    m = Chem.MolFromXYZBlock(ETHANOL)
    DB.DetermineBonds(m, charge=0)
    assert Chem.MolToSmiles(Chem.RemoveHs(m)) == "CCO"


def test_charged_fragment():
    m = Chem.MolFromXYZBlock(
        "5\nammonium\nN 0 0 0\nH 1.0 0 0\nH -0.33 0.94 0\nH -0.33 -0.47 0.82\nH -0.33 -0.47 -0.82\n"
    )
    DB.DetermineBonds(m, charge=1)
    assert Chem.GetFormalCharge(m) == 1
    assert m.GetAtomWithIdx(0).GetFormalCharge() == 1
