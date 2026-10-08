"""Chi, Hall-Kier alpha and Kappa indices (expectations from RDKit 2026.09.1)."""

import pytest

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D

NAMES = ("Chi0v", "Chi1v", "Chi2v", "Chi3v", "Chi4v", "Chi0n", "Chi1n", "Chi2n", "Chi3n", "Chi4n",
         "HallKierAlpha", "Kappa1", "Kappa2", "Kappa3")

CASES = {
    "CC(=O)Nc1ccc(O)cc1": [
        6.164862962722325, 3.2490981480278296, 2.228060242189872, 1.1876733148984595,
        0.7480265991539318, 6.164862962722325, 3.2490981480278296, 2.228060242189872,
        1.1876733148984595, 0.7480265991539318, -1.5099999999999998, 7.595374077976819,
        3.0531381005525, 2.292303972882069
    ],
    "c1ccccc1": [
        3.464101615137755, 2.0000000000000004, 1.1547005383792521, 0.6666666666666671,
        0.38490017945975075, 3.464101615137755, 2.0000000000000004, 1.1547005383792521,
        0.6666666666666671, 0.38490017945975075, -0.78, 3.4115708812260532, 1.6057694396735218,
        0.5823992601400448
    ],
    "C1CC1": [
        2.1213203435596424, 1.4999999999999998, 1.060660171779821, 0.3535533905932737, 0.0,
        2.1213203435596424, 1.4999999999999998, 1.060660171779821, 0.3535533905932737, 0.0, 0.0,
        1.3333333333333333, 0.2222222222222222, 0.0
    ],
    # beyond the first row the valence and outer-electron variants differ
    "O=C(O)c1ccccc1Br": [
        6.628823974846256, 3.4870617182335026, 2.63933574163861, 1.7957666577305706,
        1.089023902310419, 5.042827435731551, 2.6940634486761503, 1.7849988424640895,
        1.104265451624193, 0.6346499871085173, -0.8300000000000001, 7.279051254089422,
        2.8358235933154585, 1.3594078806149144
    ],
    "CCP(C)(C)=O": [
        5.456995858150284, 4.886794209741975, 5.670434130282488, 2.2846649307217692, 0.0,
        4.5625686671503685, 2.100335924038357, 1.8901447100941628, 0.7615549769072563, 0.0,
        0.22999999999999998, 6.23, 1.7902176271069716, 4.230000000000002
    ],
}


@pytest.mark.parametrize("smiles", CASES)
def test_indices(smiles):
    m = Chem.MolFromSmiles(smiles)
    values = [getattr(D, "Calc" + name)(m) for name in NAMES]
    assert values == pytest.approx(CASES[smiles], rel=1e-12, abs=1e-12)


def test_general_order_functions():
    m = Chem.MolFromSmiles("CCO")
    # order 0 of the general functions counts atom paths of length one, unlike Chi0v
    assert D.CalcChiNv(m, 0) == D.CalcChiNn(m, 0) == 3.0
    assert D.CalcChi0v(m) == pytest.approx(2.1543203766865053, rel=1e-12)
    assert D.CalcChiNn(m, 1) == pytest.approx(1.0233345472033855, rel=1e-12)
    assert D.CalcChiNv(m, 5) == 0.0
    with pytest.raises(OverflowError):
        D.CalcChiNv(m, -1)


def test_hall_kier_alpha_contributions():
    m = Chem.MolFromSmiles("*C(=O)Nc1ccccc1")
    contribs = [0.0] * m.GetNumAtoms()
    assert D.CalcHallKierAlpha(m, contribs) == pytest.approx(-1.31)
    # the dummy atom does not take part
    assert contribs == pytest.approx([0.0, -0.13, -0.2, -0.2, -0.13, -0.13, -0.13, -0.13, -0.13, -0.13])
    with pytest.raises(ValueError, match="length of atomContribs list"):
        D.CalcHallKierAlpha(m, [0.0])
    with pytest.raises(TypeError, match="Expecting an object of type list"):
        D.CalcHallKierAlpha(m, "x")


def test_deltas_are_cached_on_the_molecule():
    m = Chem.MolFromSmiles("CCO")
    assert not m.HasProp("_connectivityHKDeltas")
    D.CalcChi0v(m)
    assert m.HasProp("_connectivityHKDeltas")


def test_versions():
    assert D._CalcChi0v_version == D._CalcChiNn_version == D._CalcHallKierAlpha_version == "1.2.0"
    assert D._CalcKappa1_version == D._CalcKappa3_version == "1.1.0"
