"""rdMolAlign (expectations from RDKit 2026.09.1; the 1oir pair is Code/GraphMol/MolAlign/test_data)."""

import json
import math
import os

import numpy as np
import pytest

from rdoxide import Chem
from rdoxide.Chem import rdMolAlign as A
from rdoxide.Chem import rdMolTransforms as T
from rdoxide.Geometry import Point3D

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
with open(os.path.join(DATA, "align_1oir.json")) as f:
    PAIR = json.load(f)
with open(os.path.join(DATA, "mols3d.json")) as f:
    BLOCKS = json.load(f)


def pair():
    return Chem.MolFromMolBlock(PAIR["probe"]), Chem.MolFromMolBlock(PAIR["ref"])


def positions(mol, cid=-1):
    return mol.GetConformer(cid).GetPositions()


def moved(mol, angle=1.1, shift=(3.0, -2.0, 5.0)):
    """A copy of `mol` rotated about z and translated."""
    m = Chem.Mol(mol)
    c, s = math.cos(angle), math.sin(angle)
    t = np.array([[c, -s, 0, shift[0]], [s, c, 0, shift[1]], [0, 0, 1, shift[2]], [0, 0, 0, 1]])
    T.TransformConformer(m.GetConformer(), t)
    return m


def test_alignment_transform_matches_upstream():
    prb, ref = pair()
    rmsd, trans = A.GetAlignmentTransform(prb, ref)
    assert rmsd == pytest.approx(0.65786, abs=1e-5)
    assert trans.shape == (4, 4)
    assert trans[0] == pytest.approx([0.90854, -0.28099, 0.30918, 0.84544], abs=1e-5)
    assert trans[3] == pytest.approx([0, 0, 0, 1])
    rmsd, trans = A.GetAlignmentTransform(prb, ref, reflect=True)
    assert rmsd == pytest.approx(0.99808, abs=1e-5)


def test_atom_map_and_weights():
    prb, ref = pair()
    amap = [(i, i) for i in range(10)]
    assert A.GetAlignmentTransform(prb, ref, atomMap=amap)[0] == pytest.approx(1.78771, abs=1e-5)
    weights = [float(i) for i in range(1, 11)]
    assert A.GetAlignmentTransform(prb, ref, atomMap=amap, weights=weights)[0] == pytest.approx(3.22936, abs=1e-5)
    with pytest.raises(ValueError, match="Incorrect number of weights specified"):
        A.GetAlignmentTransform(prb, ref, atomMap=amap, weights=[1.0] * 9)
    with pytest.raises(ValueError, match="Incorrect format for an atomMap"):
        A.GetAlignmentTransform(prb, ref, atomMap=[(1, 2, 3)])


def test_align_mol_moves_the_probe_onto_the_reference():
    ref = Chem.MolFromMolBlock(BLOCKS["aspirin"], removeHs=False)
    prb = moved(ref)
    rmsd = A.AlignMol(prb, ref)
    assert rmsd == pytest.approx(0.0, abs=1e-6)
    assert np.abs(positions(prb) - positions(ref)).max() < 1e-6

    prb, ref = pair()
    assert A.AlignMol(prb, ref) == pytest.approx(0.65786, abs=1e-5)
    assert positions(prb)[0] == pytest.approx([-1.6176, 30.4574, 10.765], abs=1e-4)


def test_reflect_aligns_a_mirror_image():
    ref = Chem.MolFromMolBlock(BLOCKS["chiral"], removeHs=False)
    prb = Chem.Mol(ref)
    pos = positions(prb)
    pos[:, 0] *= -1
    conf = prb.GetConformer()
    for i, p in enumerate(pos):
        conf.SetAtomPosition(i, Point3D(*p))
    plain = A.GetAlignmentTransform(prb, ref)[0]
    mirrored = A.GetAlignmentTransform(prb, ref, reflect=True)[0]
    assert plain > 0.1
    assert mirrored == pytest.approx(0.0, abs=1e-6)


def test_best_alignment_transform_and_symmetry():
    prb, ref = pair()
    rmsd, trans, match = A.GetBestAlignmentTransform(prb, ref)
    assert rmsd == pytest.approx(0.65786, abs=1e-5)
    assert len(match) == prb.GetNumAtoms() and match[0] == (0, 9)
    assert trans.shape == (4, 4)
    # threads do not change the answer
    assert A.GetBestAlignmentTransform(prb, ref, numThreads=4)[2] == match
    params = A.BestAlignmentParams()
    params.ignoreHs = True
    assert A.GetBestAlignmentTransform(prb, ref, params)[0] == pytest.approx(0.65786, abs=1e-5)


def test_symmetric_atom_order_is_found():
    """A relabelled copy is aligned perfectly only when the match search permutes atoms."""
    ref = Chem.MolFromMolBlock(BLOCKS["benzene"], removeHs=False)
    order = list(reversed(range(ref.GetNumAtoms())))
    prb = Chem.RenumberAtoms(moved(ref), order)
    assert A.GetBestRMS(prb, ref) == pytest.approx(0.0, abs=1e-6)
    # the index-wise alignment cannot be perfect for the reversed numbering
    prb2 = Chem.RenumberAtoms(moved(ref), [3, 4, 5, 0, 1, 2, 6, 7, 8, 9, 10, 11][: ref.GetNumAtoms()])
    assert A.CalcRMS(prb2, ref, map=[[(i, i) for i in range(ref.GetNumAtoms())]]) > 0.5


def test_best_rms_is_in_place_aligned_and_threads_agree():
    prb, ref = pair()
    rms = A.GetBestRMS(Chem.Mol(prb), ref)
    assert rms == pytest.approx(0.65786, abs=1e-5)
    assert A.GetBestRMS(Chem.Mol(prb), ref, numThreads=3) == rms
    params = A.BestAlignmentParams()
    params.numThreads = 0
    moved_prb = Chem.Mol(prb)
    assert A.GetBestRMS(moved_prb, ref, params) == pytest.approx(rms)
    assert positions(moved_prb)[0] == pytest.approx([-1.6176, 30.4574, 10.765], abs=1e-4)
    # a molecule against itself needs no work
    assert A.GetBestRMS(Chem.Mol(ref), ref) == pytest.approx(0.0, abs=1e-6)
    same = Chem.Mol(ref)
    assert A.GetBestRMS(same, same) == pytest.approx(0.0, abs=1e-6)


def test_calc_rms_does_not_align():
    prb, ref = pair()
    before = positions(prb).copy()
    assert A.CalcRMS(prb, ref) == pytest.approx(30.64688, abs=1e-5)
    amap = [(i, i) for i in range(10)]
    assert A.CalcRMS(prb, ref, map=[amap]) == pytest.approx(29.99032, abs=1e-5)
    assert np.array_equal(before, positions(prb))


def test_weights_must_match_the_map():
    prb, ref = pair()
    amap = [(i, i) for i in range(10)]
    with pytest.raises(ValueError, match="Incorrect number of weights specified"):
        A.GetBestAlignmentTransform(prb, ref, map=[amap], weights=[1.0] * 3)
    params = A.BestAlignmentParams()
    params.map = [amap]
    with pytest.raises(ValueError, match="do not match number of atoms"):
        params.weights = [1.0] * 3
    params.weights = [1.0] * 10
    assert params.weights == (1.0,) * 10 and params.map == (amap,)
    assert (params.maxMatches, params.numThreads) == (1000000, 1)
    with pytest.raises(AttributeError):
        params.bogus = 1


def test_no_match_error():
    prb, ref = pair()
    other = Chem.MolFromMolBlock(BLOCKS["ethanol"], removeHs=False)
    with pytest.raises(RuntimeError, match="No sub-structure match found between the reference and probe mol"):
        A.GetBestRMS(other, ref)


def multi_conformer():
    prb, ref = pair()
    mol = Chem.Mol(ref)
    for conf in prb.GetConformers():
        mol.AddConformer(conf, assignId=True)
    mol.AddConformer(Chem.Conformer(ref.GetConformer()), assignId=True)
    return mol, ref


def test_all_conformer_best_rms():
    mol, _ = multi_conformer()
    rms = A.GetAllConformerBestRMS(mol)
    assert isinstance(rms, tuple) and len(rms) == 3
    assert rms[0] == pytest.approx(4.52572, abs=1e-5)
    assert rms[1] == pytest.approx(0.0, abs=1e-6)  # conformers 2 and 0 are identical
    assert rms[2] == pytest.approx(4.52572, abs=1e-5)
    assert A.GetAllConformerBestRMS(mol, 3) == rms
    params = A.BestAlignmentParams()
    params.numThreads = 2
    assert A.GetAllConformerBestRMS(mol, params) == rms


@pytest.mark.skipif(not hasattr(A, "GetAllConformerBestRMSToRef"), reason="needs a newer reference RDKit")
def test_all_conformer_best_rms_to_ref():
    mol, ref = multi_conformer()
    rms = A.GetAllConformerBestRMSToRef(mol, ref)
    assert len(rms) == 3
    # conformers 0 and 2 are copies of the reference's conformer
    assert rms[0] == pytest.approx(0.0, abs=1e-6)
    assert rms[1] == pytest.approx(4.52572, abs=1e-5)
    assert rms[2] == pytest.approx(0.0, abs=1e-6)
    assert A.GetAllConformerBestRMSToRef(mol, mol, A.BestAlignmentParams())[0] == pytest.approx(0.0, abs=1e-6)


def test_align_mol_conformers():
    mol, _ = multi_conformer()
    rms = []
    A.AlignMolConformers(mol, RMSlist=rms)
    assert len(rms) == 2
    assert rms[1] == pytest.approx(0.0, abs=1e-6)
    ref_pos = positions(mol, 0)
    assert np.abs(positions(mol, 2) - ref_pos).max() < 1e-6

    mol, _ = multi_conformer()
    rms = []
    A.AlignMolConformers(mol, atomIds=[0, 1, 2, 3, 4], confIds=[0, 2, 1], RMSlist=rms)
    assert len(rms) == 2
    with pytest.raises(ValueError):
        A.AlignMolConformers(mol, confIds=[0, 9])


def test_random_transform_is_reproducible_and_rigid():
    ref = Chem.MolFromMolBlock(BLOCKS["aspirin"], removeHs=False)
    a, b = Chem.Mol(ref), Chem.Mol(ref)
    A.RandomTransform(a, seed=42)
    A.RandomTransform(b, seed=42)
    assert np.allclose(positions(a), positions(b))
    assert not np.allclose(positions(a), positions(ref))
    # a rigid motion keeps the pairwise distances and is undone by an alignment
    da = np.linalg.norm(positions(a)[:, None] - positions(a)[None], axis=-1)
    dr = np.linalg.norm(positions(ref)[:, None] - positions(ref)[None], axis=-1)
    assert np.allclose(da, dr, atol=1e-9)
    assert A.AlignMol(a, ref) == pytest.approx(0.0, abs=1e-6)
