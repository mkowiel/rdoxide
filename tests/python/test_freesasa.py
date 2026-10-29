import os

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdFreeSASA

PDB = os.path.join(os.path.dirname(__file__), "..", "..", "upstream", "rdkit", "External", "FreeSASA", "test_data", "1d3z.pdb")

pytestmark = pytest.mark.skipif(not os.path.exists(PDB), reason="upstream checkout not present")


def test_protor_lee_richards_total_and_classes():
    m = Chem.MolFromPDBFile(PDB)
    radii = rdFreeSASA.classifyAtoms(m)
    assert len(radii) == m.GetNumAtoms()
    total = rdFreeSASA.CalcSASA(m, radii)
    assert total == pytest.approx(m.GetDoubleProp("SASA"))
    polar = rdFreeSASA.CalcSASA(m, radii, query=rdFreeSASA.MakeFreeSasaPolarAtomQuery())
    apolar = rdFreeSASA.CalcSASA(m, radii, query=rdFreeSASA.MakeFreeSasaAPolarAtomQuery())
    assert polar + apolar == pytest.approx(total)
    assert polar > 0 and apolar > 0


def test_shrake_rupley_close_to_lee_richards():
    m = Chem.MolFromPDBFile(PDB)
    radii = rdFreeSASA.classifyAtoms(m)
    lr = rdFreeSASA.CalcSASA(m, radii)
    opts = rdFreeSASA.SASAOpts(rdFreeSASA.ShrakeRupley, rdFreeSASA.Protor)
    sr = rdFreeSASA.CalcSASA(m, radii, opts=opts)
    assert abs(sr - lr) / lr < 0.05


# totals of the official wheel (all algorithm/classifier pairs, probe radius 0.5), equal to the last digit
WHEEL = {
    (0, 0): 6867.607195284994, (0, 1): 6819.27260178962, (0, 2): 6518.776597089304,
    (1, 0): 6897.546244392468, (1, 1): 6865.9410808814955, (1, 2): 6526.8739369295645,
}


@pytest.mark.parametrize("alg,cls", sorted(WHEEL))
def test_totals_equal_the_wheel(alg, cls):
    m = Chem.MolFromPDBFile(PDB)
    opts = rdFreeSASA.SASAOpts(alg, cls, 0.5)
    radii = rdFreeSASA.classifyAtoms(m, opts)
    total = rdFreeSASA.CalcSASA(m, radii, opts=opts)
    assert total == pytest.approx(WHEEL[(alg, cls)], rel=1e-12)
    assert total == pytest.approx(m.GetDoubleProp("SASA"))
