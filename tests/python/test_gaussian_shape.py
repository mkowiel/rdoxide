"""rdGaussianShape: API surface not covered by the upstream tests."""
import pytest
from rdoxide import Chem
from rdoxide.Chem import AllChem, rdGaussianShape


def _mol(smiles):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=42)
    return m


def test_enums_and_options():
    assert rdGaussianShape.ROTATE_180_WIGGLE == rdGaussianShape.StartMode.ROTATE_180_WIGGLE
    opts = rdGaussianShape.ShapeOverlayOptions()
    assert opts.startMode == rdGaussianShape.StartMode.A_LA_PUBCHEM
    opts.optimMode = rdGaussianShape.OptimMode.SHAPE_ONLY
    assert opts.optimMode == rdGaussianShape.OptimMode.SHAPE_ONLY
    with pytest.raises(AttributeError):
        opts.noSuchOption = 1


def test_self_score_is_one():
    m = _mol("CC(=O)Oc1ccccc1C(=O)O")
    tpl = rdGaussianShape.ScoreMol(m, m)
    assert tpl[0] == pytest.approx(1.0, abs=1e-6)


def test_align_moves_fit():
    ref = _mol("CC(=O)Oc1ccccc1C(=O)O")
    fit = Chem.Mol(ref)
    conf = fit.GetConformer()
    p = conf.GetAtomPosition(0)
    conf.SetAtomPosition(0, p.__class__(p.x + 3, p.y, p.z))
    before = fit.GetConformer().GetPositions().copy()
    tpl = rdGaussianShape.AlignMol(ref, fit)
    assert tpl[1] > 0.9
    assert (fit.GetConformer().GetPositions() != before).any()


def test_best_similarity_nothing_above_threshold():
    a = rdGaussianShape.ShapeInput(_mol("CCCCCCCC"), -1, rdGaussianShape.ShapeInputOptions(),
                                   rdGaussianShape.ShapeOverlayOptions())
    b = rdGaussianShape.ShapeInput(_mol("c1ccccc1"), -1, rdGaussianShape.ShapeInputOptions(),
                                   rdGaussianShape.ShapeOverlayOptions())
    sim, *_ = a.BestSimilarity(b, 2.0)
    assert sim == (-1.0, -1.0, -1.0)
