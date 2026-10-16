"""Gasteiger partial charges (expectations from RDKit 2026.09.1)."""

import math

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdMolDescriptors as D
from rdoxide.Chem import rdPartialCharges as P

# smiles: (atom charges, charges of the implicit hydrogens on each atom)
CASES = {
    "CC(=O)O": ([0.033767593407, 0.299684553158, -0.252820433497, -0.481433475172],
                [0.104545081232, 0.0, 0.0, 0.296256680872]),
    "c1ccccc1N": ([-0.039173741204, -0.060221846524, -0.062118053027, -0.060221846524,
                   -0.039173741204, 0.031347592591, -0.398728829487],
                  [0.06444875297, 0.06234520427, 0.062271059658, 0.06234520427, 0.06444875297, 0.0,
                   0.31243149124]),
    "CS(=O)(=O)C": ([0.038502576555, 0.144104177445, -0.229414459142, -0.229414459142,
                     0.038502576555], [0.118859793864, 0.0, 0.0, 0.0, 0.118859793864]),
    "[NH4+]": ([-0.368859027281], [1.368859027281]),
    # the formal charge is shared by the conjugated nitrogens
    "NC(=[NH2+])c1ccccc1": ([-0.287040272236, 0.270171173916, -0.287040272236, 0.062728292976,
                             -0.046778139557, -0.061366487955, -0.062215913427, -0.061366487955,
                             -0.046778139557],
                            [0.603199700436, 0.0, 0.603199700436, 0.0, 0.063213486412,
                             0.062295265238, 0.062269341858, 0.062295265238, 0.063213486412]),
    "C[Si](C)(C)C": ([-0.06974707149, 0.041113738196, -0.06974707149, -0.06974707149,
                      -0.06974707149], [0.059468636941, 0.0, 0.059468636941, 0.059468636941,
                                        0.059468636941]),
}


@pytest.mark.parametrize("smiles", CASES)
def test_charges(smiles):
    charges, h_charges = CASES[smiles]
    m = Chem.MolFromSmiles(smiles)
    P.ComputeGasteigerCharges(m)
    assert [a.GetDoubleProp("_GasteigerCharge") for a in m.GetAtoms()] == pytest.approx(
        charges, abs=1e-11)
    assert [a.GetDoubleProp("_GasteigerHCharge") for a in m.GetAtoms()] == pytest.approx(
        h_charges, abs=1e-11)


def test_charges_are_computed_properties():
    m = Chem.MolFromSmiles("CC(=O)O")
    atom = m.GetAtomWithIdx(0)
    assert not atom.HasProp("_GasteigerCharge")
    P.ComputeGasteigerCharges(m)
    assert atom.GetProp("_GasteigerCharge") == "0.033767593406775237"
    assert "_GasteigerCharge" in atom.GetPropsAsDict(True, True)
    assert "_GasteigerCharge" not in atom.GetPropsAsDict(True, False)


def test_iteration_count():
    m = Chem.MolFromSmiles("CC(=O)O")
    P.ComputeGasteigerCharges(m, 0)
    assert [a.GetDoubleProp("_GasteigerCharge") for a in m.GetAtoms()] == [0.0] * 4
    P.ComputeGasteigerCharges(m, nIter=3)
    assert [a.GetDoubleProp("_GasteigerCharge") for a in m.GetAtoms()] == pytest.approx(
        [0.026903408433, 0.316035593437, -0.240780891572, -0.47604212939], abs=1e-11)


def test_elements_without_parameters():
    m = Chem.MolFromSmiles("O=[Se]=O")
    P.ComputeGasteigerCharges(m)  # zero parameters remove the atom from the iteration
    assert all(math.isnan(a.GetDoubleProp("_GasteigerCharge")) for a in m.GetAtoms())
    with pytest.raises(ValueError, match="Element: Se Mode: sp2"):
        P.ComputeGasteigerCharges(m, 12, True)
    with pytest.raises(ValueError, match="Element: Se Mode: sp2"):
        P.ComputeGasteigerCharges(m, throwOnParamFailure=1)


def test_peoe_vsa():
    bins = D.PEOE_VSA_(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))
    assert [round(x, 10) for x in bins] == [
        10.4233159988, 4.7945371841, 0.0, 0.0, 0.0, 0.0, 24.2654682738, 12.6111234744, 0.0,
        5.7495118333, 0.0, 5.9071797294, 0.0, 0.0
    ]
