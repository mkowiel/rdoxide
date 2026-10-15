from rdkit import Chem
from rdkit.Chem import rdFMCS


def _mols(*smiles):
    return [Chem.MolFromSmiles(s) for s in smiles]


def test_find_mcs_defaults():
    res = rdFMCS.FindMCS(_mols("c1ccccc1O", "c1ccccc1N"))
    assert (res.numAtoms, res.numBonds) == (6, 6)
    assert res.smartsString == "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
    assert not res.canceled
    assert res.queryMol.GetNumAtoms() == 6


def test_parameters_object_shares_compare_parameters():
    ps = rdFMCS.MCSParameters()
    ps.AtomCompareParameters.RingMatchesRingOnly = True
    assert ps.AtomCompareParameters.RingMatchesRingOnly
    ps.Threshold = 0.5
    assert ps.Threshold == 0.5
    assert ps.AtomTyper == rdFMCS.AtomCompare.CompareElements


def test_custom_atom_compare_and_progress():
    class AnyAtom(rdFMCS.MCSAtomCompare):
        def __call__(self, p, mol1, atom1, mol2, atom2):
            return True

    ps = rdFMCS.MCSParameters()
    ps.AtomTyper = AnyAtom()
    res = rdFMCS.FindMCS(_mols("CCO", "CCN"), ps)
    assert res.numAtoms == 3

    class Stop(rdFMCS.MCSProgress):
        def __call__(self, stat, params):
            return False

    ps = rdFMCS.MCSParameters()
    ps.ProgressCallback = Stop()
    assert rdFMCS.FindMCS(_mols("c1ccccc1CCC", "c1ccccc1CCO"), ps).canceled


def test_callback_exceptions_propagate():
    class Boom(rdFMCS.MCSAtomCompare):
        def __call__(self, p, mol1, atom1, mol2, atom2):
            raise KeyError("boom")

    ps = rdFMCS.MCSParameters()
    ps.AtomTyper = Boom()
    try:
        rdFMCS.FindMCS(_mols("CC", "CC"), ps)
    except KeyError as e:
        assert "boom" in str(e)
    else:
        raise AssertionError("no exception")
