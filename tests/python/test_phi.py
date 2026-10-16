"""CalcPhi (Hall-Kier based Phi, ConnectivityDescriptors) with expectations from RDKit 2026.09.1."""

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdMolDescriptors as D


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CCCC", 3.0),
        ("c1ccccc1", 0.9130327103921437),
        ("CC(=O)O", 0.7589395748987852),
        ("OCCN", 2.8616),
        ("c1ccncc1", 0.8653513194512445),
        ("CS(=O)(=O)C", 0.9612633641691973),
        ("FC(F)(F)Cl", 1.0637709972299172),
        ("C1CCCCC1", 1.5432098765432098),
        ("CCBr", 2.8768),
        # heteroatoms without an alpha entry use the periodic-table fallback
        ("B(O)(O)C", 1.3149794850034497),
        ("CP(C)C", 1.9066112973760918),
        ("C[Si](C)(C)C", 1.1780339967422768),
        ("CCO[Se]C", 4.909045248777196),
        ("[H]C([H])([H])[H]", 0.0),
    ],
)
def test_calc_phi(smiles, expected):
    assert D.CalcPhi(Chem.MolFromSmiles(smiles)) == pytest.approx(expected, rel=1e-12)
