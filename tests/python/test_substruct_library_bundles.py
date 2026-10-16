"""SubstructLibrary searches with MolBundle and ExtendedQueryMol queries (the latter is
upstream's test_using_xqms, which upstream skips without Boost serialization)."""
import pytest
from rdkit import Chem
from rdkit.Chem import rdGeneralizedSubstruct, rdSubstructLibrary

SMIS = ["COCC=O", "COOCC=O", "COOOCC=O", "COOOOCC=O"]


@pytest.mark.parametrize(
    "holder_cls", [rdSubstructLibrary.CachedSmilesMolHolder,
                   rdSubstructLibrary.CachedTrustedSmilesMolHolder])
@pytest.mark.parametrize("with_fps", [False, True])
@pytest.mark.parametrize("smiles,expected", [
    ("COCC", [0]),
    ("COC=CO", [0]),
    ("COCC |LN:1:1.3|", [0, 1, 2]),
    ("COC=CO |LN:1:1.3|", [0, 1, 2]),
    ("CNC=CO |LN:1:1.3|", []),
])
def test_using_xqms(holder_cls, with_fps, smiles, expected):
    holder = holder_cls()
    for smi in SMIS:
        holder.AddSmiles(smi)
    lib = rdSubstructLibrary.SubstructLibrary(holder)
    if with_fps:
        rdSubstructLibrary.AddPatterns(lib, rdSubstructLibrary.TautomerPatternHolder())
    xqm = rdGeneralizedSubstruct.CreateExtendedQueryMol(Chem.MolFromSmiles(smiles))
    assert list(lib.GetMatches(xqm)) == expected
    assert lib.HasMatch(xqm) == bool(expected)
    assert lib.CountMatches(xqm) == len(expected)
    for n in (1, 2, 4):
        assert list(lib.GetMatches(xqm, numThreads=n)) == expected


def test_bundle_max_results_dedup_and_threads():
    ssl = rdSubstructLibrary.SubstructLibrary(rdSubstructLibrary.MolHolder(),
                                              rdSubstructLibrary.PatternHolder())
    for smi in ('CCOC', 'CCNC', 'COOCOO', 'CCNC', 'CCCC', 'CCOCC'):
        ssl.AddMol(Chem.MolFromSmiles(smi))
    bndl = Chem.MolBundle()
    for smi in ('COC', 'CCC', 'CC'):  # 'CC' re-matches most molecules found before
        bndl.AddMol(Chem.MolFromSmiles(smi))
    allm = list(ssl.GetMatches(bndl))
    assert len(allm) == len(set(allm))
    assert sorted(allm) == [0, 1, 3, 4, 5]
    assert ssl.CountMatches(bndl) == 5
    assert ssl.HasMatch(bndl)
    assert len(ssl.GetMatches(bndl, maxResults=2)) == 2
    for n in (1, 3):
        assert sorted(ssl.GetMatches(bndl, numThreads=n)) == [0, 1, 3, 4, 5]
