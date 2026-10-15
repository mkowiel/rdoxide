import pytest

from rdkit import Chem
from rdkit.Chem import ChemicalForceFields, rdtrajectory
from rdkit.Chem.rdtrajectory import Snapshot, Trajectory


def test_snapshot_points_need_a_trajectory():
    s = Snapshot([1.0, 2.0, 3.0], 4.5)
    assert s.GetEnergy() == 4.5
    s.SetEnergy(1.0)
    assert s.GetEnergy() == 1.0
    with pytest.raises(Exception):
        s.GetPoint3D(0)


def test_trajectory_editing():
    traj = Trajectory(3, 2)
    assert traj.Dimension() == 3 and traj.NumPoints() == 2
    assert traj.AddSnapshot(Snapshot([0, 0, 0, 1, 1, 1], 1.0)) == 0
    assert traj.AddSnapshot(Snapshot([2, 2, 2, 3, 3, 3], 2.0)) == 1
    assert traj.InsertSnapshot(1, Snapshot([9, 9, 9, 8, 8, 8], 3.0)) == 1
    assert len(traj) == 3
    p = traj.GetSnapshot(1).GetPoint3D(1)
    assert (p.x, p.y, p.z) == (8.0, 8.0, 8.0)
    traj.RemoveSnapshot(1)
    assert [traj.GetSnapshot(i).GetEnergy() for i in range(len(traj))] == [1.0, 2.0]
    copy = Trajectory(traj)
    traj.Clear()
    assert len(traj) == 0 and len(copy) == 2
    with pytest.raises(Exception):
        traj.GetSnapshot(0)


def test_minimize_trajectory_to_conformers():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, randomSeed=42)
    field = ChemicalForceFields.UFFGetMoleculeForceField(mol)
    res, snapshots = field.MinimizeTrajectory(1, 200, 1e-3)
    assert res in (0, 1) and len(snapshots) >= 1
    traj = Trajectory(3, mol.GetNumAtoms(), snapshots)
    n = mol.GetNumConformers()
    assert traj.AddConformersToMol(mol) == len(snapshots)
    assert mol.GetNumConformers() == n + len(snapshots)
    energies = [s.GetEnergy() for s in snapshots]
    assert energies[-1] <= energies[0] + 1e-6
