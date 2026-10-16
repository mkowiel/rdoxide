"""Fingerprint generator JSON; expected strings come from RDKit 2026.03.6."""

import json

import pytest
from rdoxide import Chem
from rdoxide.Chem import rdFingerprintGenerator as rfg

GENERATORS = {
    "atompair": lambda: rfg.GetAtomPairGenerator(),
    "atompair_opts": lambda: rfg.GetAtomPairGenerator(
        minDistance=2, maxDistance=9, includeChirality=True, fpSize=512
    ),
    "morgan": lambda: rfg.GetMorganGenerator(radius=3, countSimulation=True, fpSize=512),
    "morgan_feature": lambda: rfg.GetMorganGenerator(
        atomInvariantsGenerator=rfg.GetMorganFeatureAtomInvGen()
    ),
    "morgan_bond": lambda: rfg.GetMorganGenerator(
        bondInvariantsGenerator=rfg.GetMorganBondInvGen(useBondTypes=False, useChirality=True)
    ),
    "rdkit": lambda: rfg.GetRDKitFPGenerator(),
    "rdkit_opts": lambda: rfg.GetRDKitFPGenerator(
        minPath=2, maxPath=5, useHs=False, branchedPaths=False, numBitsPerFeature=3
    ),
    "torsion": lambda: rfg.GetTopologicalTorsionGenerator(),
    "torsion_opts": lambda: rfg.GetTopologicalTorsionGenerator(
        torsionAtomCount=5, countSimulation=False, includeChirality=True
    ),
}

MOLS = [Chem.MolFromSmiles(s) for s in ("CCC", "c1ccccc1O", "C[C@H](F)C(=O)N")]


def test_default_atompair_json_text():
    assert rfg.GetAtomPairGenerator().ToJSON() == (
        '{"name":"FingerprintGenerator","fingerprintArguments":{"type":"AtomPairArguments",'
        '"use2D":"true","minDistance":"1","maxDistance":"30","countSimulation":"true",'
        '"fpSize":"2048","numBitsPerFeature":"1","includeChirality":"false",'
        '"countBounds":["1","2","4","8"]},"atomEnvironmentGenerator":'
        '{"type":"AtomPairEnvGenerator"},"atomInvariantsGenerator":'
        '{"type":"AtomPairAtomInvGenerator","includeChirality":"false",'
        '"topologicalTorsionCorrection":"false"}}'
    )


@pytest.mark.parametrize("name", GENERATORS)
def test_round_trip_preserves_fingerprints(name):
    g1 = GENERATORS[name]()
    text = g1.ToJSON()
    json.loads(text)
    g2 = rfg.FingerprintGeneratorFromJSON(text)
    assert g2.ToJSON() == text
    for m in MOLS:
        assert g1.GetFingerprint(m) == g2.GetFingerprint(m)
        assert g1.GetSparseCountFingerprint(m) == g2.GetSparseCountFingerprint(m)


def test_missing_values_use_constructor_defaults():
    g = rfg.FingerprintGeneratorFromJSON(
        '{"fingerprintArguments":{"type":"MorganArguments","radius":"x","countBounds":[1,2],'
        '"countSimulation":true},"atomEnvironmentGenerator":{"type":"MorganEnvGenerator"}}'
    )
    assert g.ToJSON() == (
        '{"name":"FingerprintGenerator","fingerprintArguments":{"type":"MorganArguments",'
        '"onlyNonzeroInvariants":"false","radius":"3","countSimulation":"true","fpSize":"2048",'
        '"numBitsPerFeature":"1","includeChirality":"false","countBounds":["1","2"]},'
        '"atomEnvironmentGenerator":{"type":"MorganEnvGenerator"}}'
    )


@pytest.mark.parametrize(
    "text,exc",
    [
        ("", RuntimeError),
        ("{", RuntimeError),
        ('{"fingerprintArguments":{"type":"Foo"}}', ValueError),
        (
            '{"fingerprintArguments":{"type":"MorganArguments"},'
            '"atomEnvironmentGenerator":{"type":"Bar"}}',
            ValueError,
        ),
    ],
)
def test_bad_json(text, exc):
    with pytest.raises(exc):
        rfg.FingerprintGeneratorFromJSON(text)
