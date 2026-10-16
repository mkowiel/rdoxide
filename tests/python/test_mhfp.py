"""rdoxide.Chem.rdMHFPFingerprint; expected values come from reference RDKit."""

import pytest
from rdoxide import Chem
from rdoxide.Chem import rdMHFPFingerprint as M


@pytest.fixture(scope="module")
def encoder():
    return M.MHFPEncoder(16, 42)


def test_from_array(encoder):
    assert list(encoder.FromArray([1, 2, 3])) == [
        734796314, 152805157, 1046592797, 839703249, 2584932064, 1099380393,
        2221926165, 283729444, 1482179631, 287907250, 330695626, 2713252445,
        1872093025, 915334311, 747481211, 720244703,
    ]


def test_encode_smiles(encoder):
    assert list(encoder.EncodeSmiles("CCO")) == [
        998603975, 113829445, 99459723, 323656167, 2238542606, 634808226,
        563983127, 697535188, 1029691159, 886737314, 147326604, 1461843357,
        1051908113, 102127419, 169658043, 248371099,
    ]


def test_shingling_is_order_independent_set(encoder):
    assert sorted(encoder.CreateShinglingFromSmiles("CCO")) == ["C(C)O", "CC", "CCO", "OC", "OCC"]


def test_secfp(encoder):
    assert list(encoder.EncodeSECFPSmiles("CCO", length=64).GetOnBits()) == [4, 7, 12, 15, 35]
    fp = encoder.EncodeSECFPMol(Chem.MolFromSmiles("CCO"), length=64)
    assert list(fp.GetOnBits()) == [4, 7, 12, 15, 35]


def test_distance(encoder):
    assert M.MHFPEncoder.Distance(encoder.EncodeSmiles("CCO"), encoder.EncodeSmiles("CCN")) == 0.875


def test_bulk_matches_single(encoder):
    smiles = ["CCO", "c1ccccc1O", "CC(=O)N"]
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    single = [list(encoder.EncodeSmiles(s)) for s in smiles]
    assert [list(v) for v in encoder.EncodeSmilesBulk(smiles)] == single
    assert [list(v) for v in encoder.EncodeMolsBulk(mols)] == single
    fps = encoder.EncodeSECFPMolsBulk(mols, length=128)
    assert [list(f.GetOnBits()) for f in fps] == [
        list(encoder.EncodeSECFPSmiles(s, length=128).GetOnBits()) for s in smiles
    ]


def test_bad_smiles_raises(encoder):
    with pytest.raises(RuntimeError):
        encoder.EncodeSmiles("not")
