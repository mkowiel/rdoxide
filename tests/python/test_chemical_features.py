"""rdChemicalFeatures / rdMolChemicalFeatures (expectations from RDKit 2026.09.1)."""

import pickle

import pytest

from rdoxide import Chem
from rdoxide.Chem import ChemicalFeatures
from rdoxide.Geometry import Point3D

FDEF = """
# donors and acceptors
AtomType NDonor [N&!H0&v3,N&!H0&+1&v4,n&H1&+0]
AtomType Acceptor [O,N&H0]

DefineFeature HDonor [{NDonor};!$(NC=O)]
  Family Donor
  Weights 1.0
EndFeature

DefineFeature Carbonyl [O]=[C]\\
   -[N,O]
  Family Group
  Weights 1,2,1
EndFeature

DefineFeature HAcceptor [{Acceptor}]
  Family Acceptor
EndFeature
"""

# a planar fragment: O=C(N)OC with 3D coordinates
MOLBLOCK = """
  RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    1.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000   -1.0000    0.5000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    2.0000    1.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  1  3  1  0
  1  4  1  0
  3  5  1  0
M  END
"""


@pytest.fixture(scope="module")
def factory():
    return ChemicalFeatures.BuildFeatureFactoryFromString(FDEF)


def test_factory_definitions(factory):
    assert factory.GetNumFeatureDefs() == 3
    assert factory.GetFeatureFamilies() == ("Donor", "Group", "Acceptor")
    defs = factory.GetFeatureDefs()
    assert list(defs) == ["Donor.HDonor", "Group.Carbonyl", "Acceptor.HAcceptor"]
    assert defs["Group.Carbonyl"] == "[O]=[C]-[N,O]"
    assert defs["Acceptor.HAcceptor"].startswith("[$([O,N&H0])]")


def test_features_of_a_molecule(factory):
    mol = Chem.MolFromSmiles("NC(=O)OC")
    feats = factory.GetFeaturesForMol(mol)
    summary = [(f.GetFamily(), f.GetType(), f.GetAtomIds(), f.GetId()) for f in feats]
    # the amide nitrogen is excluded from the donors by the definition
    assert summary == [
        ("Group", "Carbonyl", (2, 1, 0), 1),
        ("Group", "Carbonyl", (2, 1, 3), 2),
        ("Acceptor", "HAcceptor", (2,), 3),
        ("Acceptor", "HAcceptor", (3,), 4),
    ]
    assert factory.GetNumMolFeatures(mol) == len(feats)
    assert factory.GetNumMolFeatures(mol, includeOnly="Group") == 2
    assert len(factory.GetFeaturesForMol(mol, includeOnly="Nothing")) == 0
    assert feats[0].GetMol().GetNumAtoms() == 5
    assert feats[0].GetFactory().GetNumFeatureDefs() == 3


def test_feature_positions_use_the_weights(factory):
    mol = Chem.MolFromMolBlock(MOLBLOCK)
    feats = {f.GetType(): f for f in factory.GetFeaturesForMol(mol, includeOnly="Group")}
    carbonyl = feats["Carbonyl"]
    atoms = carbonyl.GetAtomIds()
    conf = mol.GetConformer()
    weights = [0.25, 0.5, 0.25]
    expected = [sum(w * list(conf.GetAtomPosition(a))[k] for w, a in zip(weights, atoms)) for k in range(3)]
    assert list(carbonyl.GetPos()) == pytest.approx(expected)
    assert list(carbonyl.GetPos(0)) == pytest.approx(expected)
    assert carbonyl.GetActiveConformer() == -1

    # positions are cached until the cache is cleared
    conf.SetAtomPosition(atoms[0], Point3D(9, 9, 9))
    assert list(carbonyl.GetPos()) == pytest.approx(expected)
    carbonyl.ClearCache()
    assert list(carbonyl.GetPos()) != pytest.approx(expected)

    single = factory.GetFeaturesForMol(mol, includeOnly="Acceptor")[0]
    assert list(single.GetPos()) == list(conf.GetAtomPosition(single.GetAtomIds()[0]))


def test_position_errors(factory):
    flat = Chem.MolFromSmiles("NC(=O)OC")
    feat = factory.GetFeaturesForMol(flat)[0]
    with pytest.raises(RuntimeError, match="molecule has no conformers"):
        feat.GetPos()
    mol = Chem.MolFromMolBlock(MOLBLOCK)
    feat = factory.GetFeaturesForMol(mol)[0]
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        feat.GetPos(5)


def test_get_mol_feature(factory):
    mol = Chem.MolFromSmiles("NC(=O)OC")
    n = factory.GetNumMolFeatures(mol)
    assert factory.GetMolFeature(mol, n - 1).GetId() == n
    with pytest.raises(IndexError):
        factory.GetMolFeature(mol, n)
    with pytest.raises(IndexError):
        factory.GetMolFeature(mol, -1)
    feat = factory.GetMolFeature(mol, 0, confId=3)
    assert feat.GetActiveConformer() == 3


def test_get_atom_match(factory):
    mol = Chem.MolFromSmiles("NC(=O)OC")
    feats = factory.GetFeaturesForMol(mol, includeOnly="Group")
    assert ChemicalFeatures.GetAtomMatch(feats[:1]) == [list(feats[0].GetAtomIds())]
    assert ChemicalFeatures.GetAtomMatch(feats) == []  # they share atoms


@pytest.mark.parametrize(
    "text, message",
    [
        ("Bogus x\n", "line 1:\nbad or missing keyword"),
        ("AtomType A\n", "bad AtomType line, missing definition"),
        ("AtomType A [N\n", "invalid SMARTS in AtomType ({A}): [$([N)]"),
        ("DefineFeature a [N]\nFamily F\n", "could not find EndFeature line for feature: a"),
        ("DefineFeature a [N]\nEndFeature\n", "did not find Family definition for feature: a"),
        ("DefineFeature a [N]\nFamily F\nWeights q\nEndFeature", "bad weight value (q) for feature: a"),
        ("DefineFeature a [N]\nFamily F\nWeights 1,2\nEndFeature\n", "len(feature weight vector)"),
    ],
)
def test_definition_errors(text, message):
    with pytest.raises(ValueError) as err:
        ChemicalFeatures.BuildFeatureFactoryFromString(text)
    assert message in str(err.value)


def test_missing_file():
    with pytest.raises(OSError, match="could not be opened"):
        ChemicalFeatures.BuildFeatureFactory("/nonexistent/features.fdef")


def test_free_chemical_feature():
    feat = ChemicalFeatures.FreeChemicalFeature("Donor", "Typ", Point3D(1, 2, 3), 5)
    assert (feat.GetId(), feat.GetFamily(), feat.GetType()) == (5, "Donor", "Typ")
    assert list(feat.GetPos()) == [1, 2, 3]
    other = pickle.loads(pickle.dumps(feat))
    assert (other.GetId(), other.GetFamily(), other.GetType(), list(other.GetPos())) == (5, "Donor", "Typ", [1, 2, 3])

    partial = ChemicalFeatures.FreeChemicalFeature("Acc", Point3D(0, 1, 0))
    assert (partial.GetId(), partial.GetType()) == (-1, "")
    blank = ChemicalFeatures.FreeChemicalFeature()
    blank.SetId(3)
    blank.SetFamily("x")
    blank.SetType("y")
    blank.SetPos(Point3D(1, 1, 1))
    assert (blank.GetId(), blank.GetFamily(), blank.GetType(), list(blank.GetPos())) == (3, "x", "y", [1, 1, 1])
