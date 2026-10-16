"""rdkit.Chem.rdEnumerateStereoisomers; expected values come from reference RDKit."""

import pytest
from rdkit import Chem
from rdkit.Chem import rdEnumerateStereoisomers as E


def _all(mol, **opts):
    options = E.StereoEnumerationOptions()
    for k, v in opts.items():
        setattr(options, k, v)
    en = E.StereoisomerEnumerator(mol, options, False)
    out = []
    while (iso := en.next()) is not None:
        out.append(Chem.MolToCXSmiles(iso))
    return en, out


def test_basic_count():
    en, isos = _all(Chem.MolFromSmiles("CC(F)=CC(Cl)C"))
    assert en.GetStereoisomerCount() == 4
    assert len({Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in isos}) == 4


def test_seeded_sequence_is_reproducible():
    mol = Chem.MolFromSmiles("FC(Cl)C(F)C(Cl)C(F)Cl")
    first = _all(mol, randomSeed=42, maxIsomers=3)[1]
    assert len(first) == 3
    assert first == _all(mol, randomSeed=42, maxIsomers=3)[1]


def test_many_centers_and_chiral_flag():
    mol = Chem.MolFromSmiles("CC(F)=CC(Cl)C" * 101)
    opts = E.StereoEnumerationOptions()
    opts.maxIsomers = 13
    en = E.StereoisomerEnumerator(mol, opts, False)
    count = 0
    while (iso := en.next()) is not None:
        assert iso.GetProp("_MolFileChiralFlag") == "1"
        count += 1
    assert count == 13


def test_unknown_option_rejected():
    with pytest.raises(AttributeError):
        E.StereoEnumerationOptions().bogus = 1


def test_try_embedding():
    # upstream catch_tests "Embedding"
    mol = Chem.MolFromSmiles("BrC=CC1OC(C2)(F)C2(Cl)C1")
    opts = E.StereoEnumerationOptions()
    opts.tryEmbedding = True
    enum = E.StereoisomerEnumerator(mol, opts)
    got = set()
    while (iso := enum.next()) is not None:
        got.add(Chem.MolToSmiles(iso))
        assert iso.GetNumConformers() == 1
    assert got == {
        "F[C@@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2",
        "F[C@@]12C[C@]1(Cl)C[C@@H](/C=C\\Br)O2",
        "F[C@@]12C[C@]1(Cl)C[C@H](/C=C/Br)O2",
        "F[C@@]12C[C@]1(Cl)C[C@H](/C=C\\Br)O2",
        "F[C@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2",
        "F[C@]12C[C@@]1(Cl)C[C@@H](/C=C\\Br)O2",
        "F[C@]12C[C@@]1(Cl)C[C@H](/C=C/Br)O2",
        "F[C@]12C[C@@]1(Cl)C[C@H](/C=C\\Br)O2",
    }
