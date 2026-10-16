"""MQNs and 2D autocorrelation descriptors (expectations from RDKit 2026.09.1)."""

import pytest

from rdoxide import Chem
from rdoxide.Chem import Descriptors
from rdoxide.Chem import rdMolDescriptors as D

MQNS = {
    "CC(=O)Nc1ccc(O)cc1": [8, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 11, 4, 1, 0, 3, 3, 0, 1, 5, 3, 2, 2,
                           0, 0, 3, 1, 1, 0, 4, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    "C12CC3CC(CC(C3)C1)C2": [10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 12, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 0, 6, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 10, 12],
    "[O-][N+](=O)c1ccccc1": [6, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 9, 2, 1, 0, 3, 3, 0, 1, 6, 3, 0, 0,
                             1, 1, 2, 0, 1, 0, 5, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
}


@pytest.mark.parametrize("smiles", MQNS)
def test_mqns(smiles):
    assert D.MQNs_(Chem.MolFromSmiles(smiles)) == MQNS[smiles]


def test_mqns_need_ring_info_and_ignore_force():
    m = Chem.MolFromSmiles("CCO", sanitize=False)
    with pytest.raises(RuntimeError):
        D.MQNs_(m)
    m = Chem.MolFromSmiles("CCO")
    assert D.MQNs_(m, force=True) == D.MQNs_(m)


def test_autocorr2d_layout():
    values = D.CalcAUTOCORR2D(Chem.MolFromSmiles("CCO"))
    assert len(values) == 192
    # weighting 0 (mass): plain, centered, Moran and Geary families over lags 1-8
    assert values[0:3] == [1.204, 0.847, 0.0]
    assert values[48:50] == [0.037, 0.025]
    assert values[96:98] == [-0.25, -1.0]
    assert values[144:146] == [0.75, 1.5]
    assert D._CalcAUTOCORR2D_version == "1.0.0"


def test_autocorr2d_custom_property():
    m = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    for atom, v in zip(m.GetAtoms(), range(10)):
        atom.SetDoubleProp("p", v + 1.5)
    values = D.CalcAUTOCORR2D(m, "p")
    assert values == [5.95, 6.205, 5.775, 5.159, 5.043, 4.208, 3.882, 0.0,
                      78.514, 78.44, 59.748, 57.81, 55.215, 30.039, 24.087, 0.0,
                      0.317, 0.301, -0.138, -0.032, -0.962, -1.047, -1.259, 0.0,
                      0.285, -0.036, 0.033, -0.481, -0.78, -0.721, -0.784, 0.0]
    # a missing property weighs every atom 1
    assert D.CalcAUTOCORR2D(Chem.MolFromSmiles("CCO"), "missing")[:3] == [1.099, 0.693, 0.0]


def test_autocorr2d_bad_property_and_descriptor_list():
    m = Chem.MolFromSmiles("CCO")
    m.GetAtomWithIdx(0).SetIntProp("q", 3)
    with pytest.raises(RuntimeError, match="bad any cast"):
        D.CalcAUTOCORR2D(m, "q")
    Descriptors.setupAUTOCorrDescriptors()
    names = [name for name, _ in Descriptors._descList]
    assert "AUTOCORR2D_1" in names and "AUTOCORR2D_192" in names
