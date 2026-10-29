"""rdEHTTools against results recorded from the official RDKit wheel (tools/make_yaehmop_corpus.py)."""
import json
from pathlib import Path

import numpy as np
import pytest

from rdoxide import Chem
from rdoxide.Chem import rdEHTTools

CASES = json.loads((Path(__file__).resolve().parents[1] / "data" / "yaehmop_eht.json").read_text())


@pytest.mark.parametrize("case", CASES, ids=[c["smiles"] for c in CASES])
def test_matches_wheel(case):
    mol = Chem.MolFromMolBlock(case["molblock"], removeHs=False, sanitize=False)
    ok, res = rdEHTTools.RunMol(mol, keepOverlapAndHamiltonianMatrices=True)
    assert ok
    assert res.numOrbitals == case["numOrbitals"]
    assert res.numElectrons == case["numElectrons"]
    assert res.fermiEnergy == pytest.approx(case["fermiEnergy"], abs=1e-9)
    assert res.totalEnergy == pytest.approx(case["totalEnergy"], abs=1e-8)
    np.testing.assert_allclose(res.GetAtomicCharges(), case["charges"], atol=1e-9)
    np.testing.assert_allclose(res.GetOrbitalEnergies(), case["energies"], atol=1e-9)
    np.testing.assert_allclose(res.GetReducedOverlapPopulationMatrix(), case["rop"], atol=1e-9)
    np.testing.assert_allclose(res.GetReducedChargeMatrix().ravel(), case["rchg"], atol=1e-8)
    assert res.GetHamiltonian().sum() == pytest.approx(case["hamiltonian_sum"], abs=1e-8)
    assert res.GetOverlapMatrix().sum() == pytest.approx(case["overlap_sum"], abs=1e-8)


def test_matrices_need_the_flag():
    mol = Chem.MolFromSmiles("O")
    mol = Chem.AddHs(mol)
    from rdoxide.Chem import AllChem

    AllChem.EmbedMolecule(mol, randomSeed=1)
    ok, res = rdEHTTools.RunMol(mol)
    assert ok
    with pytest.raises(ValueError):
        res.GetHamiltonian()
    with pytest.raises(ValueError):
        res.GetOverlapMatrix()
