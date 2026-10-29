from rdoxide import Geometry
from rdoxide.Chem import MolFromSmiles, rdCoordGen, rdDepictor


def _bond_lengths(m):
    c = m.GetConformer()
    return [(c.GetAtomPosition(b.GetBeginAtomIdx()) - c.GetAtomPosition(b.GetEndAtomIdx())).Length() for b in m.GetBonds()]


def test_add_coords_gives_unit_bonds():
    m = MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    rdCoordGen.AddCoords(m)
    assert m.GetNumConformers() == 1
    assert not m.GetConformer().Is3D()
    assert all(abs(x - 1.0) < 0.05 for x in _bond_lengths(m))


def test_macrocycle():
    m = MolFromSmiles("C1COCCOCCOCCOCCOCCO1")
    rdCoordGen.AddCoords(m)
    assert all(abs(x - 1.0) < 0.05 for x in _bond_lengths(m))


def test_coord_map():
    m = MolFromSmiles("CCCCC")
    p = rdCoordGen.CoordGenParams()
    p.SetCoordMap({0: Geometry.Point2D(0, 0), 1: Geometry.Point2D(1.5, 0)})
    rdCoordGen.AddCoords(m, p)
    c = m.GetConformer()
    assert abs(c.GetAtomPosition(1).x - c.GetAtomPosition(0).x) > 1.0


def test_prefer_coordgen_context():
    assert rdDepictor.IsCoordGenSupportAvailable()
    assert not rdDepictor.GetPreferCoordGen()
    with rdDepictor.UsingCoordGen(True):
        assert rdDepictor.GetPreferCoordGen()
        m = MolFromSmiles("CCO")
        rdDepictor.Compute2DCoords(m)
    assert not rdDepictor.GetPreferCoordGen()
