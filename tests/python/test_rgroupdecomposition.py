import pytest

from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition as rgd
from rdkit.Chem.rdRGroupDecomposition import (RGroupCoreAlignment, RGroupDecompose,
                                              RGroupDecomposition, RGroupDecompositionParameters,
                                              RGroupLabelling, RGroupLabels, RGroupMatching,
                                              RGroupScore, RelabelMappedDummies)

SYMMETRY = ["c1(Cl)ccccc1", "c1c(Cl)cccc1", "c1cccc(Cl)c1", "c1cc(Cl)ccc1", "c1ccc(Cl)cc1"]


def test_enums_are_exported():
    assert RGroupLabels.IsotopeLabels == 1
    assert rgd.AutoDetect == RGroupLabels.AutoDetect == 0xFF
    assert RGroupMatching.GA == 0x10
    assert RGroupCoreAlignment.MCS == 1
    assert getattr(RGroupCoreAlignment, "None") == RGroupCoreAlignment.NoAlignment
    assert RGroupScore.FingerprintVariance == 4
    assert RGroupLabelling.AtomMap | RGroupLabelling.MDLRGroup == 5


def test_parameters():
    params = RGroupDecompositionParameters()
    assert params.labels == RGroupLabels.AutoDetect
    assert params.matchingStrategy == RGroupMatching.GreedyChunks
    assert params.rgroupLabelling == RGroupLabelling.AtomMap | RGroupLabelling.MDLRGroup
    assert params.chunkSize == 5
    assert params.timeout == -1.0
    params.chunkSize = 3
    params.onlyMatchAtRGroups = True
    assert params.chunkSize == 3 and params.onlyMatchAtRGroups
    with pytest.raises(AttributeError):
        params.notAnAttribute = 1
    # the substructure parameters are a live object
    assert params.substructMatchParams.useChirality
    params.substructMatchParams.useChirality = False
    assert not params.substructMatchParams.useChirality
    with pytest.raises(AttributeError):
        params.substructMatchParams = None


def test_rows_and_columns():
    core = Chem.MolFromSmiles("c1ccccc1")
    mols = [Chem.MolFromSmiles(s) for s in SYMMETRY] + [Chem.MolFromSmiles("CC")]
    rows, unmatched = RGroupDecompose([core], mols, asSmiles=True)
    assert unmatched == [5]
    assert rows == [{"Core": "c1ccc([*:1])cc1", "R1": "Cl[*:1]"}] * 5
    cols, unmatched = RGroupDecompose(core, mols, asSmiles=True, asRows=False)
    assert unmatched == [5]
    assert cols == {"Core": ["c1ccc([*:1])cc1"] * 5, "R1": ["Cl[*:1]"] * 5}


def test_class_interface():
    core = Chem.MolFromSmiles("c1ccccc1")
    decomp = RGroupDecomposition(core)
    for i, smi in enumerate(SYMMETRY):
        assert decomp.Add(Chem.MolFromSmiles(smi)) == i
    assert decomp.Add(Chem.MolFromSmiles("CC")) == -1
    ok, score = decomp.ProcessAndScore()
    assert ok and score > 0
    assert decomp.GetRGroupLabels() == ["Core", "R1"]
    rows = decomp.GetRGroupsAsRows()
    assert all(isinstance(m, Chem.Mol) for m in rows[0].values())
    assert decomp.GetMatchingCoreIdx(Chem.MolFromSmiles("c1ccccc1O")) == 0
    matches = []
    assert decomp.GetMatchingCoreIdx(Chem.MolFromSmiles("c1ccccc1O"), matches) == 0
    # all 12 symmetry-related mappings of the ring are reported
    assert len(matches) == 12 and all(len(m) == 6 for m in matches)


def test_none_molecule_is_rejected():
    core = Chem.MolFromSmiles("c1ccccc1")
    with pytest.raises(ValueError):
        RGroupDecompose([core], [None])
    with pytest.raises(ValueError):
        RGroupDecomposition([None])


def test_fingerprint_variance_and_labels():
    core = Chem.MolFromSmiles("c1ccccc1")
    params = RGroupDecompositionParameters()
    params.scoreMethod = RGroupScore.FingerprintVariance
    params.rgroupLabelling = RGroupLabelling.Isotope
    rows, _ = RGroupDecompose([core], [Chem.MolFromSmiles(s) for s in SYMMETRY], asSmiles=True,
                              options=params)
    assert rows[0] == {"Core": "[1*]c1ccccc1", "R1": "[1*]Cl"}


def test_unsupported_features_are_reported():
    core = Chem.MolFromSmiles("c1ccccc1")
    params = RGroupDecompositionParameters()
    params.matchingStrategy = RGroupMatching.GA
    decomp = RGroupDecomposition(core, params)
    assert decomp.Add(Chem.MolFromSmiles(SYMMETRY[0])) == 0
    assert not decomp.Process()
    params = RGroupDecompositionParameters()
    params.doEnumeration = True
    with pytest.raises(RuntimeError):
        RGroupDecomposition(core, params)


def test_relabel_mapped_dummies():
    mol = Chem.MolFromSmiles("[*:1]c1ccccc1[2*]")
    RelabelMappedDummies(mol)
    assert Chem.MolToCXSmiles(mol).startswith("*c1ccccc1*")
    assert {a.GetPropsAsDict().get("dummyLabel") for a in mol.GetAtoms()} >= {"R1", "R2"}
