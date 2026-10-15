"""ResonanceMolSupplier; the expectations come from RDKit's own rough_test.py."""

import pytest
from rdkit import Chem


def total_charge(mol):
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


AMIDINE = "NC(=[NH2+])c1ccc(cc1)C(=O)[O-]"


def test_no_conjugated_group():
    suppl = Chem.ResonanceMolSupplier(Chem.MolFromSmiles("CC"))
    assert suppl.GetNumConjGrps() == 0
    assert len(suppl) == 1


def test_counts_for_the_flags():
    mol = Chem.MolFromSmiles(AMIDINE)
    suppl = Chem.ResonanceMolSupplier(mol)
    assert not suppl.GetIsEnumerated()
    assert len(suppl) == 4
    assert suppl.GetIsEnumerated()
    assert len(Chem.ResonanceMolSupplier(mol, Chem.KEKULE_ALL)) == 8
    flags = Chem.ALLOW_INCOMPLETE_OCTETS | Chem.UNCONSTRAINED_CATIONS | Chem.UNCONSTRAINED_ANIONS
    suppl = Chem.ResonanceMolSupplier(mol, flags)
    assert len(suppl) == 32
    assert all(total_charge(suppl[i]) == total_charge(mol) for i in range(len(suppl)))
    assert len(Chem.ResonanceMolSupplier(mol, flags, 10)) == 10


def test_iteration_and_indexing():
    mol = Chem.MolFromSmiles("CC(=O)[O-]")
    suppl = Chem.ResonanceMolSupplier(mol)
    first = next(suppl)
    assert not suppl.atEnd()
    next(suppl)
    assert suppl.atEnd()
    with pytest.raises(StopIteration):
        next(suppl)
    suppl.reset()
    assert Chem.MolToSmiles(next(suppl)) == Chem.MolToSmiles(first)
    assert Chem.MolToSmiles(suppl[-1]) == "CC(=O)[O-]"
    with pytest.raises(IndexError):
        suppl[2]
    assert len(list(suppl)) == 2


def test_alternating_kekule_structures():
    suppl = Chem.ResonanceMolSupplier(Chem.MolFromSmiles(AMIDINE), Chem.KEKULE_ALL)
    kinds = {suppl[i].GetBondBetweenAtoms(3, 4).GetBondType() for i in range(2)}
    assert len(kinds) == 2


def test_conjugated_group_indices():
    mol = Chem.MolFromSmiles("CC(=O)[O-]")
    suppl = Chem.ResonanceMolSupplier(mol)
    assert suppl.GetNumConjGrps() == 1
    assert [suppl.GetAtomConjGrpIdx(i) for i in range(4)] == [-1, 0, 0, 0]
    assert [suppl.GetBondConjGrpIdx(i) for i in range(3)] == [-1, 0, 0]


def test_substructure_matches():
    suppl = Chem.ResonanceMolSupplier(Chem.MolFromSmiles("CC(=O)[O-]"))
    query = Chem.MolFromSmarts("C(=O)[O-]")
    assert suppl.GetSubstructMatches(query) == ((1, 2, 3), (1, 3, 2))
    assert suppl.GetSubstructMatches(query, uniquify=True) == ((1, 2, 3),)
    assert suppl.GetSubstructMatch(query) in [(1, 2, 3), (1, 3, 2)]
    dmap = Chem.ResonanceMolSupplier(Chem.MolFromSmiles("C(C)Nc1cc[nH+]cc1"))
    matches = dmap.GetSubstructMatches(Chem.MolFromSmarts("[#7+]"), False, False, False)
    assert sorted(m[0] for m in matches) == [2, 6]


def test_threads_give_the_same_matches():
    mol = Chem.MolFromSmiles("NC(=[NH2+])N.CC(=O)[O-]")
    suppl = Chem.ResonanceMolSupplier(mol)
    suppl.SetNumThreads(0)
    query = Chem.MolFromSmarts("NC(=[NH2+])N")
    assert suppl.GetSubstructMatches(query, numThreads=0) == suppl.GetSubstructMatches(query)


class StopAfter(Chem.ResonanceMolSupplierCallback):
    def __init__(self, limit):
        super().__init__()
        self.limit = limit

    def __call__(self):
        assert self.GetNumConjGrps() == 1
        return self.GetNumStructures(0) < self.limit


MOL = Chem.MolFromSmiles(
    "ClC1=NC(NC2=CC=CC3=C2C(=O)C2=CC=CC=C2C3=O)=NC(NC2=CC=CC3=C2C(=O)C2=CC=CC=C2C3=O)=N1"
)


def test_progress_callback_cancels():
    suppl = Chem.ResonanceMolSupplier(MOL, Chem.ALLOW_CHARGE_SEPARATION)
    assert suppl.GetProgressCallback() is None
    callback = StopAfter(12)
    suppl.SetProgressCallback(callback)
    assert suppl.GetProgressCallback() is callback
    suppl.SetProgressCallback(None)
    assert suppl.GetProgressCallback() is None
    suppl.SetProgressCallback(StopAfter(12))
    assert len(suppl) == 12
    assert suppl.WasCanceled()


def test_big_molecule_counts():
    assert len(Chem.ResonanceMolSupplier(MOL)) == 1
    assert len(Chem.ResonanceMolSupplier(MOL, Chem.KEKULE_ALL)) == 32
    suppl = Chem.ResonanceMolSupplier(MOL, Chem.ALLOW_CHARGE_SEPARATION, 10)
    assert len(suppl) == 10
    assert not suppl.WasCanceled()


def test_callback_validation():
    class NotACallback(Chem.ResonanceMolSupplier):
        def __call__(self):
            return True

    class NoCall(Chem.ResonanceMolSupplierCallback):
        pass

    suppl = Chem.ResonanceMolSupplier(MOL, Chem.ALLOW_CHARGE_SEPARATION)
    with pytest.raises(TypeError):
        suppl.SetProgressCallback(NotACallback())
    with pytest.raises(TypeError):
        suppl.SetProgressCallback(object())
    with pytest.raises(AttributeError):
        suppl.SetProgressCallback(NoCall())


def test_callback_exceptions_propagate():
    class Raises(Chem.ResonanceMolSupplierCallback):
        def __call__(self):
            raise ValueError("boom")

    suppl = Chem.ResonanceMolSupplier(MOL, Chem.ALLOW_CHARGE_SEPARATION)
    suppl.SetProgressCallback(Raises())
    with pytest.raises(ValueError, match="boom"):
        suppl.Enumerate()
