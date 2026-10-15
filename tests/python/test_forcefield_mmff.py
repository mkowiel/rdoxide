"""MMFF bindings (rdkit.Chem.rdForceFieldHelpers, rdkit.ForceField)."""

import pytest
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers as ffh


def mol_with_conf(smiles="CC(=O)O"):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    conf = Chem.Conformer(mol.GetNumAtoms())
    # a crude zig-zag so that no two atoms coincide
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (1.2 * i, 0.7 * (i % 3), 0.5 * (i % 2)))
    mol.AddConformer(conf)
    return mol


def test_properties():
    mol = mol_with_conf()
    props = ffh.MMFFGetMoleculeProperties(mol)
    assert props is not None
    assert [props.GetMMFFAtomType(i) for i in range(4)] == [1, 3, 7, 6]
    assert props.GetMMFFVariant() == "MMFF94"
    props.SetMMFFVariant("MMFF94s")
    assert props.GetMMFFVariant() == "MMFF94s"
    with pytest.raises(Exception):
        props.SetMMFFVariant("MMFF95")
    total = sum(props.GetMMFFPartialCharge(i) for i in range(mol.GetNumAtoms()))
    assert total == pytest.approx(0.0, abs=1e-9)
    bond = props.GetMMFFBondStretchParams(mol, 0, 1)
    assert bond[0] == 0 and bond[1] > 0 and bond[2] > 0
    assert props.GetMMFFOopBendParams(mol, 0, 1, 2, 3) is not None
    assert len(props.GetMMFFVdWParams(0, 2)) == 4


def test_invalid_molecule_has_no_properties():
    mol = Chem.AddHs(Chem.MolFromSmiles("[Pt]"))
    assert ffh.MMFFGetMoleculeProperties(mol) is None
    assert not ffh.MMFFHasAllMoleculeParams(mol)


def test_force_field_and_optimization():
    mol = mol_with_conf()
    props = ffh.MMFFGetMoleculeProperties(mol)
    ff = ffh.MMFFGetMoleculeForceField(mol, props)
    e0 = ff.CalcEnergy()
    assert ff.Minimize(maxIts=500) == 0
    assert ff.CalcEnergy() < e0
    assert ffh.MMFFOptimizeMolecule(mol, maxIters=500) == 0
    ff2 = ffh.MMFFGetMoleculeForceField(mol, props)
    props.SetMMFFEleTerm(False)
    ff3 = ffh.MMFFGetMoleculeForceField(mol, props)
    assert ff3.CalcEnergy() != ff2.CalcEnergy()


def test_conformers_independent_of_threads():
    mol = mol_with_conf()
    for shift in (0.1, 0.2):
        conf = Chem.Conformer(mol.GetConformer())
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x + shift * (i % 2), p.y, p.z - shift))
        mol.AddConformer(conf, assignId=True)
    one = ffh.MMFFOptimizeMoleculeConfs(Chem.Mol(mol), numThreads=1, maxIters=500)
    many = ffh.MMFFOptimizeMoleculeConfs(Chem.Mol(mol), numThreads=0, maxIters=500)
    assert one == many
    assert len(one) == 3


def test_optimize_with_a_generic_force_field():
    mol = mol_with_conf()
    props = ffh.MMFFGetMoleculeProperties(mol)
    ff = ffh.MMFFGetMoleculeForceField(mol, props)
    assert ffh.OptimizeMolecule(ff, maxIters=500) == 0
    empty = ffh.CreateEmptyForceFieldForMol(mol)
    assert empty.NumPoints() == mol.GetNumAtoms()
