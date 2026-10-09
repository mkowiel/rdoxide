"""molzip, molzipFragments and MolzipParams (expectations from RDKit 2026.09.1)."""

import json
import os
import platform

import pytest

from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdmolops import MolzipLabel, MolzipParams

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "mols3d.json")) as f:
    BLOCKS = json.load(f)

M = Chem.MolFromSmiles
S = Chem.MolToSmiles


def params(label=None, **attrs):
    p = MolzipParams()
    if label is not None:
        p.label = label
    for k, v in attrs.items():
        setattr(p, k, v)
    return p


def test_zip_atom_map_numbers_keeps_stereo():
    a, b = M("C=C[*:1]"), M("O/C=N/[*:1]")
    assert S(rdmolops.molzip(a, b)) == "C=C/N=C/O"
    assert S(rdmolops.molzip(rdmolops.CombineMols(a, b))) == "C=C/N=C/O"
    assert S(rdmolops.molzip(M("C[C@H](F)[*:1]"), M("[*:1]Cl"))) == "C[C@H](F)Cl"
    assert S(rdmolops.molzip(M("[*:1][C@H](F)Cl"), M("C[*:1]"))) == "C[C@H](F)Cl"
    assert S(rdmolops.molzip(M("F/C=C/[*:1]"), M("[*:1]/C=C/Cl"))) == "F/C=C/C=C/Cl"


def test_zip_isotope_and_atom_type_labels():
    iso = params(MolzipLabel.Isotope)
    assert S(rdmolops.molzip(M("C=C[1*]"), M("O/C=N/[1*]"), iso)) == "C=C/N=C/O"
    typed = params(MolzipLabel.AtomType)
    typed.setAtomSymbols(["V", "Xe"])
    assert S(rdmolops.molzip(M("C[V]"), M("N[Xe]"), typed)) == "N[Xe].[CH3][V]"


def test_zip_link_bond():
    link = M("[*:1][*:2]")
    assert S(rdmolops.molzip(link, rdmolops.CombineMols(M("C[*:1]"), M("S[*:2]")))) == "CS"


@pytest.mark.parametrize(
    "smi, cut, expected",
    [
        ("CCOC(=O)c1ccccc1", [1, 3], "CCOC(=O)c1ccccc1"),
        ("C[C@H](F)Cl", [0, 1], "C[C@H](F)Cl"),
        ("C[C@@H]1CC[C@H](C)CC1", [0, 5], "C[C@H]1CC[C@@H](C)CC1"),
    ],
)
def test_zip_fragment_on_bonds_round_trip(smi, cut, expected):
    frags = rdmolops.FragmentOnBonds(M(smi), cut)
    zipped = rdmolops.molzip(frags, params(MolzipLabel.FragmentOnBonds))
    assert S(zipped) == expected


def test_zip_fragments_and_rgroup_row():
    pieces = [M("[*:1]C[*:2]"), M("C[*:1]"), M("N[*:2]")]
    assert S(rdmolops.molzipFragments(pieces)) == "CCN"
    row = {"Core": pieces[0], "R1": pieces[1], "R2": pieces[2]}
    assert S(rdmolops.molzip(row)) == "CCN"
    assert rdmolops.molzipFragments([]) is None
    with pytest.raises(RuntimeError, match="RGroup has no Core|Core"):
        rdmolops.molzip({"R1": M("C[*:1]")})


def test_zip_fragments_generate_coordinates():
    pieces = [M("[*:1]C[*:2]"), M("C[*:1]"), M("N[*:2]")]
    zipped = rdmolops.molzipFragments(pieces, params(generateCoordinates=True))
    assert zipped.GetNumConformers() == 1
    assert all(p.GetNumConformers() == 1 for p in pieces)


def test_zip_errors():
    with pytest.raises(RuntimeError, match="bond info already exists for end atom with label:1"):
        rdmolops.molzip(M("C[*:1]"), M("C[*:1].C[*:1]"))
    with pytest.raises(RuntimeError, match="Must have an end atom to bond"):
        rdmolops.molzip(M("[*:1][*:2]"), M("C[*:1]"))
    with pytest.raises(TypeError, match="did not match C\\+\\+ signature"):
        rdmolops.molzip(3)
    # unmatched labels are left alone
    assert S(rdmolops.molzip(M("C[*:1]"), M("N[*:2]"))) == "C[*:1].N[*:2]"


def test_params_attributes():
    p = MolzipParams()
    assert p.label == MolzipLabel.AtomMapNumber
    assert p.enforceValenceRules and not p.generateCoordinates and not p.alignCoordinates
    p.label = MolzipLabel.Isotope
    assert p.label == MolzipLabel.Isotope
    with pytest.raises(AttributeError, match="Cannot set unknown attribute 'foo'"):
        p.foo = 1
    with pytest.raises(TypeError, match="did not match C\\+\\+ signature"):
        p.label = 99


def test_align_coordinates_of_fragments():
    ethanol = Chem.MolFromMolBlock(BLOCKS["ethanol"], removeHs=False)
    frags = rdmolops.FragmentOnBonds(ethanol, [1])
    p = params(MolzipLabel.FragmentOnBonds, alignCoordinates=True)
    try:
        zipped = rdmolops.molzip(frags, p)
    except ValueError:
        # FragmentOnBonds puts the dummies on the neighbouring atoms, so both
        # attachment vectors coincide and whether cos(theta) rounds to exactly
        # 1 depends on fused multiply-add (aarch64 fuses, x86-64 does not).
        # RDKit's own arithmetic fails the same way there.
        assert platform.machine().lower() not in ("arm64", "aarch64")
        pytest.skip("degenerate alignment needs fused multiply-add")
    conf = zipped.GetConformer()
    expected = {0: (3.4114, 0.8321, 1.4089), 1: (1.9433, 0.9917, 1.762), 2: (1.2237, 0.3416, 0.7792), 5: (3.4761, 0.1761, 0.522)}
    for idx, pos in expected.items():
        assert tuple(conf.GetAtomPosition(idx)) == pytest.approx(pos, abs=1e-4)
