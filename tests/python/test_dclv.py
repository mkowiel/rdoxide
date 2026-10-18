"""DoubleCubicLatticeVolume of a single atom has closed forms."""

import math

import pytest
from rdoxide import Chem, Geometry
from rdoxide.Chem import rdMolDescriptors


def single_atom():
  mol = Chem.RWMol()
  mol.AddAtom(Chem.Atom(18))
  conf = Chem.Conformer(1)
  conf.SetAtomPosition(0, Geometry.Point3D(0, 0, 0))
  mol.AddConformer(conf)
  return mol


def test_isolated_sphere():
  radius, probe = 2.0, 1.4
  dclv = rdMolDescriptors.DoubleCubicLatticeVolume(single_atom(), [radius], probeRadius=probe)
  assert dclv.GetSurfaceArea() == pytest.approx(4 * math.pi * (radius + probe)**2, rel=1e-4)
  assert dclv.GetVolume() == pytest.approx(4 / 3 * math.pi * (radius + probe)**3, rel=1e-4)
  assert dclv.GetVDWVolume() == pytest.approx(4 / 3 * math.pi * radius**3, rel=1e-4)
  assert dclv.GetPackingDensity() == pytest.approx((radius / (radius + probe))**3, rel=1e-4)
  assert len(dclv.GetSurfacePoints()[0]) == 320
  assert len(dclv.GetSurfacePoints(allPoints=True)[0]) == 320


def test_arguments_and_errors():
  mol = single_atom()
  by_position = rdMolDescriptors.DoubleCubicLatticeVolume(mol, [2.0], False, True, 1.2)
  by_name = rdMolDescriptors.DoubleCubicLatticeVolume(mol, radii=[2.0], probeRadius=1.2)
  assert by_position.GetSurfaceArea() == by_name.GetSurfaceArea()
  default_radius = rdMolDescriptors.DoubleCubicLatticeVolume(mol, isProtein=False)
  assert default_radius.GetSurfaceArea() > 0
  with pytest.raises(ValueError):
    rdMolDescriptors.DoubleCubicLatticeVolume(mol, [1.0, 2.0])
  with pytest.raises(ValueError):
    by_name.GetPartialSurfaceArea([])
