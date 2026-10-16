"""Bulk fingerprint functions of rdFingerprintGenerator."""

import pytest
from rdoxide import Chem
from rdoxide.Chem import rdFingerprintGenerator as rfg

SMILES = ["CCC", "OCCCCC", "c1ccccc1O", "CC(=O)Nc1ccc(O)cc1", "C1CCC1C(=O)O"] * 8
GENERATORS = {
    rfg.MorganFP: lambda: rfg.GetMorganGenerator(2),
    rfg.RDKitFP: rfg.GetRDKitFPGenerator,
    rfg.AtomPairFP: rfg.GetAtomPairGenerator,
    rfg.TopologicalTorsionFP: rfg.GetTopologicalTorsionGenerator,
}


@pytest.mark.parametrize("fp_type", list(GENERATORS))
@pytest.mark.parametrize("threads", [1, 3, 0])
def test_bulk_matches_single(fp_type, threads):
    mols = [Chem.MolFromSmiles(s) for s in SMILES]
    gen = GENERATORS[fp_type]()
    got = rfg.GetFPs(mols, fp_type, threads)
    assert len(got) == len(mols)
    assert all(g == gen.GetFingerprint(m) for g, m in zip(got, mols))
    got = rfg.GetSparseCountFPs(mols, fp_type, threads)
    assert all(g == gen.GetSparseCountFingerprint(m) for g, m in zip(got, mols))


def test_default_is_morgan_and_empty_input():
    mols = [Chem.MolFromSmiles("CCO")]
    assert rfg.GetCountFPs(mols)[0] == rfg.GetMorganGenerator(2).GetCountFingerprint(mols[0])
    assert rfg.GetSparseFPs([]) == []
