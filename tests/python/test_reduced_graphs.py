"""rdReducedGraphs (expectations from the reference RDKit)."""

import pytest
from rdoxide import Chem
from rdoxide.Chem import rdReducedGraphs as rg


def sparse(fp):
    return {i: round(float(v), 6) for i, v in enumerate(fp) if v}


def test_fingerprints():
    fp = rg.GetErGFingerprint(Chem.MolFromSmiles("OCCc1ccccc1"))
    assert len(fp) == 315
    assert sparse(fp) == {77: 0.3, 78: 1.0, 79: 0.3, 152: 0.3, 153: 1.0, 154: 0.3}
    short = rg.GetErGFingerprint(Chem.MolFromSmiles("OCCc1ccccc1"), fuzzIncrement=0.1, minPath=2, maxPath=6)
    assert len(short) == 105
    assert sparse(short)[27] == 1.0


def test_reduced_graph_and_its_fingerprint():
    mol = Chem.MolFromSmiles("OCCc1ccccc1")
    graph = rg.GenerateMolExtendedReducedGraph(mol)
    graph.UpdatePropertyCache(False)
    assert Chem.MolToSmiles(graph) == "*cCCO"
    types = [tuple(a.GetPropsAsDict().get("_ErGAtomTypes")) for a in graph.GetAtoms()]
    assert types == [(0, 1), (), (), (), (5,)]
    assert list(rg.GenerateErGFingerprintForReducedGraph(graph)) == pytest.approx(
        list(rg.GetErGFingerprint(mol))
    )


def test_errors():
    mol = Chem.MolFromSmiles("CCO")
    with pytest.raises(ValueError, match="atom types"):
        rg.GetErGFingerprint(mol, atomTypes=[1])
    with pytest.raises(RuntimeError):
        rg.GetErGFingerprint(mol, minPath=3, maxPath=3)
