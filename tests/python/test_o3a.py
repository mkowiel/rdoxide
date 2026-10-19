import pytest

from rdoxide import Chem
from rdoxide.Chem import AllChem, rdMolAlign


def _confs(smiles, n=3, seed=7):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMultipleConfs(m, n, randomSeed=seed)
    AllChem.MMFFOptimizeMoleculeConfs(m)
    return m


def test_o3a_aligns_conformers_of_the_same_molecule():
    ref = _confs("CC(=O)Nc1ccc(O)cc1", 1)
    prb = _confs("CC(=O)Nc1ccc(O)cc1", 3, seed=11)
    o3a = rdMolAlign.GetO3A(prb, ref, prbCid=0, refCid=0)
    assert o3a.Score() > 0
    matches = o3a.Matches()
    assert len(matches) >= 3 and len(matches) == len(o3a.Weights())
    rmsd, mat = o3a.Trans()
    assert mat.shape == (4, 4)
    assert o3a.Align() == pytest.approx(rmsd, abs=1e-6)


def test_crippen_o3a_and_probe_conformer_lists_agree_with_single_calls():
    ref = _confs("c1ccccc1CCN", 1)
    prb = _confs("c1ccccc1CCO", 3)
    single = [rdMolAlign.GetCrippenO3A(prb, ref, prbCid=cid).Score() for cid in range(3)]
    many = [o.Score() for o in rdMolAlign.GetCrippenO3AForProbeConfs(prb, ref, numThreads=2)]
    assert many == pytest.approx(single)
    mmff = rdMolAlign.GetO3AForProbeConfs(prb, ref)
    assert len(mmff) == 3


def test_constraints_are_validated():
    ref = _confs("CCO", 1)
    prb = _confs("CCO", 1)
    with pytest.raises(ValueError, match="out of range"):
        rdMolAlign.GetO3A(prb, ref, constraintMap=[(99, 0)])
    with pytest.raises(ValueError, match="heavy"):
        rdMolAlign.GetO3A(prb, ref, constraintMap=[(3, 0)])
    with pytest.raises(ValueError, match="number of weights"):
        rdMolAlign.GetO3A(prb, ref, constraintMap=[(0, 0)], constraintWeights=[1.0, 2.0])
