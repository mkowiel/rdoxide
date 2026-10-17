"""The pure-Python packages copied from RDKit import and work on top of the native modules."""

import numpy
from rdoxide import Chem, Geometry
from rdoxide.Chem import ChemicalFeatures
from rdoxide.Chem.FeatMaps import FeatMapPoint
from rdoxide.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdoxide.ML.Data import Stats
from rdoxide.VLib.NodeLib import SmilesSupply  # noqa: F401  (imports the node library)


def test_free_chemical_feature_can_be_subclassed():
  pt = FeatMapPoint.FeatMapPoint()
  pt.initFromFeat(ChemicalFeatures.FreeChemicalFeature("Aromatic", "Foo", Geometry.Point3D(1, 2, 3)))
  assert pt.GetFamily() == "Aromatic"
  assert list(pt.GetPos()) == [1.0, 2.0, 3.0]
  assert pt.featDirs == []


def test_pharm2d_signature_of_a_molecule():
  fp = Generate.Gen2DFingerprint(Chem.MolFromSmiles("OCC=CC(=O)O"), Gobbi_Pharm2D.factory)
  assert fp.GetNumOnBits() > 0


def test_ml_stats():
  assert abs(Stats.StandardizeMatrix(numpy.array([1.0, 2.0, 3.0, 4.0])).sum()) < 1e-9
