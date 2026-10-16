"""rdSynthonSpaceSearch: API surface not covered by the upstream tests."""
import os
from pathlib import Path

import pytest
from rdoxide import Chem
from rdoxide.Chem import rdSynthonSpaceSearch

DATA = Path(__file__).resolve().parents[2] / "upstream" / "rdkit" / "Code" / "GraphMol" / "SynthonSpaceSearch" / "data"
pytestmark = pytest.mark.skipif(not DATA.exists(), reason="upstream test data not checked out")


def test_text_round_trip(tmp_path):
    space = rdSynthonSpaceSearch.SynthonSpace()
    space.ReadTextFile(str(DATA / "amide_space.txt"))
    assert space.GetNumReactions() == 1
    assert space.GetNumProducts() == 12
    assert space.GetNumSynthons() == 7
    out = tmp_path / "amide.spc"
    space.WriteDBFile(str(out))
    space2 = rdSynthonSpaceSearch.SynthonSpace()
    space2.ReadDBFile(str(out))
    assert space2.GetNumProducts() == 12
    q = Chem.MolFromSmarts("c1ccccc1C(=O)N1CCCC1")
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.maxHits = -1
    r1 = space.SubstructureSearch(q, params=params)
    r2 = space2.SubstructureSearch(q, params=params)
    n1 = sorted(m.GetProp("_Name") for m in r1.GetHitMolecules())
    n2 = sorted(m.GetProp("_Name") for m in r2.GetHitMolecules())
    assert n1 and n1 == n2


def test_params_are_checked():
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    with pytest.raises(AttributeError):
        params.noSuchOption = 1
    params.maxHits = 5
    assert params.maxHits == 5
    assert params.similarityCutoff == 0.5


def test_enumerated_file(tmp_path):
    space = rdSynthonSpaceSearch.SynthonSpace()
    space.ReadTextFile(str(DATA / "amide_space.txt"))
    out = tmp_path / "enum.smi"
    space.WriteEnumeratedFile(str(out))
    assert len(out.read_text().splitlines()) == 12


def test_formatted_integer_string():
    assert rdSynthonSpaceSearch.FormattedIntegerString(1234567) == "1 234 567"
    assert rdSynthonSpaceSearch.FormattedIntegerString(12) == "12"


def test_rascal_and_fingerprint_search():
    from rdoxide.Chem import rdFingerprintGenerator, rdRascalMCES

    space = rdSynthonSpaceSearch.SynthonSpace()
    space.ReadTextFile(str(DATA / "amide_space.txt"))
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.maxHits = -1
    q = Chem.MolFromSmiles("c1ccccc1C(=O)N1CCCC1")
    res = space.RascalSearch(q, rdRascalMCES.RascalOptions(), params)
    assert res.GetHitMolecules()
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    params.similarityCutoff = 0.3
    res = space.FingerprintSearch(q, fpgen, params)
    assert all(float(m.GetProp("Similarity")) >= 0.3 for m in res.GetHitMolecules())
    assert space.GetSynthonFingerprintType() == fpgen.GetInfoString()


SHAPE_QUERY = (
    "O=C(c1ccccc1)N1CCCC1 |(0.0443291,-1.81486,-1.76886;0.0506321,-0.858174,-0.921491;"
    "1.37975,-0.430412,-0.483603;2.18964,-1.35506,0.144714;3.47088,-1.00454,0.585539;"
    "3.93803,0.297573,0.388032;3.1267,1.22739,-0.242406;1.85597,0.849751,-0.670531;"
    "-1.14837,-0.261434,-0.446583;-1.26073,0.836916,0.520219;-2.73583,1.04666,0.696614;"
    "-3.34033,-0.283345,0.290893;-2.46516,-0.679843,-0.874401)|")


def test_shape_search():
    space = rdSynthonSpaceSearch.SynthonSpace()
    space.ReadDBFile(str(DATA / "amide_space_shapes.spc"))
    assert space.GetNumSynthonsWithShapes() == 7
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.similarityCutoff = 0.8
    params.fragSimilarityAdjuster = 0.2
    params.approxSimilarityAdjuster = 0.2
    params.numConformers = 100
    params.timeOut = 0
    params.randomSeed = 0xdac
    params.bestHit = True
    query = Chem.MolFromSmiles(SHAPE_QUERY)
    hits = space.ShapeSearch(query, params).GetHitMolecules()
    assert sorted(m.GetProp("_Name") for m in hits) == [
        "1-1;2-1;amide-1", "1-2;2-1;amide-1", "1-3;2-1;amide-1"]
    assert all(float(m.GetProp("Similarity")) > 0.9 for m in hits)


def test_build_shapes_and_convert(tmp_path):
    params = rdSynthonSpaceSearch.ShapeBuildParams()
    params.numConfs = 10
    params.randomSeed = 0xdac
    out = tmp_path / "amide_shapes.spc"
    rdSynthonSpaceSearch.ConvertTextToDBFile(str(DATA / "amide_space.txt"), str(out), py_shapeParams=params)
    space = rdSynthonSpaceSearch.SynthonSpace()
    space.ReadDBFile(str(out))
    assert space.GetNumSynthonsWithShapes() == 7
    # a second build is a no-op
    space.BuildSynthonShapes(params)
    assert space.GetNumSynthonsWithShapes() == 7
