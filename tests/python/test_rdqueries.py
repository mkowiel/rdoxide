"""rdqueries, QueryAtom/QueryBond; verified against RDKit 2026.09.1."""

import pytest
from rdoxide import Chem
from rdoxide.Chem import rdqueries as q


def matching(mol, query_atom):
    return [a.GetIdx() for a in mol.GetAtomsMatchingQuery(query_atom)]


def test_factories_return_query_atoms():
    qa = q.AtomNumEqualsQueryAtom(6)
    assert type(qa).__name__ == "QueryAtom"
    assert isinstance(qa, Chem.Atom) and isinstance(qa, Chem.QueryAtom)
    assert qa.HasQuery()
    assert qa.GetSmarts() == "[#6]"
    assert qa.DescribeQuery() == "AtomAtomicNum 6 = val\n"
    assert q.AtomNumEqualsQueryAtom(6, negate=True).GetSmarts() == "[!#6]"


def test_matching_atoms():
    mol = Chem.MolFromSmiles("CC(=O)Nc1ccccc1")
    assert matching(mol, q.AtomNumEqualsQueryAtom(7)) == [3]
    assert matching(mol, q.IsAromaticQueryAtom()) == [4, 5, 6, 7, 8, 9]
    assert matching(mol, q.IsAromaticQueryAtom(negate=True)) == [0, 1, 2, 3]
    assert matching(mol, q.IsInRingQueryAtom()) == [4, 5, 6, 7, 8, 9]
    assert matching(mol, q.HCountEqualsQueryAtom(3)) == [0]


def test_less_and_greater_compare_the_atom_with_the_value():
    mol = Chem.MolFromSmiles("CCO")
    assert matching(mol, q.AtomNumLessQueryAtom(7)) == [0, 1]
    assert matching(mol, q.AtomNumLessQueryAtom(6)) == []
    assert matching(mol, q.AtomNumGreaterQueryAtom(6)) == [2]
    assert matching(mol, q.AtomNumLessQueryAtom(7, negate=True)) == [2]
    assert q.AtomNumLessQueryAtom(7).DescribeQuery() == "AtomNumLess 7 > \n"


def test_property_queries():
    mol = Chem.MolFromSmiles("CCO")
    mol.GetAtomWithIdx(0).SetIntProp("i", 5)
    mol.GetAtomWithIdx(1).SetDoubleProp("d", 2.5)
    mol.GetAtomWithIdx(2).SetProp("s", "x")
    assert matching(mol, q.HasPropQueryAtom("i")) == [0]
    assert matching(mol, q.HasPropQueryAtom("i", True)) == [1, 2]
    assert matching(mol, q.HasIntPropWithValueQueryAtom("i", 4, False, 1)) == [0]
    assert matching(mol, q.HasDoublePropWithValueQueryAtom("d", 2.4, False, 0.2)) == [1]
    assert matching(mol, q.HasStringPropWithValueQueryAtom("s", "x")) == [2]
    bond_query = q.HasPropQueryBond("p")
    assert type(bond_query).__name__ == "QueryBond"


def test_expand_and_set_query():
    qa = q.AtomNumEqualsQueryAtom(6)
    qa.ExpandQuery(q.IsAromaticQueryAtom())
    assert qa.GetSmarts() == "[#6&a]"
    qa.ExpandQuery(q.HCountEqualsQueryAtom(1), Chem.CompositeQueryType.COMPOSITE_OR)
    assert qa.DescribeQuery().startswith("AtomOr")
    qa.SetQuery(q.AtomNumEqualsQueryAtom(7))
    assert qa.GetSmarts() == "[#7]"
    with pytest.raises(Exception) as err:
        qa.ExpandQuery(Chem.Atom(6))
    assert type(err.value).__name__ == "ArgumentError"


def test_atoms_of_query_molecules_are_query_atoms():
    pattern = Chem.MolFromSmarts("[#6;R]-,=[#7]")
    assert [type(a).__name__ for a in pattern.GetAtoms()] == ["QueryAtom"] * 2
    assert type(pattern.GetBondWithIdx(0)).__name__ == "QueryBond"
    assert type(Chem.MolFromSmiles("CC").GetAtomWithIdx(0)).__name__ == "Atom"
    assert type(Chem.AtomFromSmarts("[C;R]")).__name__ == "QueryAtom"
    assert type(Chem.BondFromSmarts("-,=")).__name__ == "QueryBond"
    with pytest.raises(RuntimeError):
        Chem.QueryAtom(6)


def test_query_atoms_match_other_atoms():
    ring_carbon = Chem.AtomFromSmarts("[C;R]")
    assert ring_carbon.Match(Chem.MolFromSmiles("C1CC1").GetAtomWithIdx(0))
    assert not ring_carbon.Match(Chem.MolFromSmiles("CC").GetAtomWithIdx(0))
    assert q.AtomNumEqualsQueryAtom(6).Match(Chem.Atom(6))
    double_or_single = Chem.BondFromSmarts("-,=")
    assert double_or_single.Match(Chem.MolFromSmiles("C=N").GetBondWithIdx(0))
    assert not double_or_single.Match(Chem.MolFromSmiles("C#N").GetBondWithIdx(0))


def test_recursive_query_atoms_in_matching():
    mol = Chem.MolFromSmiles("CCO")
    assert matching(mol, Chem.AtomFromSmarts("[$(C-O)]")) == [1]


def test_replace_atom_with_query_atom():
    mol = Chem.RWMol(Chem.MolFromSmiles("CCO"))
    replaced = q.ReplaceAtomWithQueryAtom(mol, mol.GetAtomWithIdx(1))
    assert type(replaced).__name__ == "QueryAtom"
    assert replaced.GetIdx() == 1
    replaced.ExpandQuery(q.IsInRingQueryAtom(True))
    assert Chem.MolToSmarts(mol) == "[#6]-[#6&!R]-[#8]"
    # nothing changes when there already is a query
    assert q.ReplaceAtomWithQueryAtom(mol, replaced).GetIdx() == 1


def test_bad_arguments_raise_argument_error():
    for call in (
        lambda: q.AtomNumEqualsQueryAtom("x"),
        lambda: q.AtomNumEqualsQueryAtom(),
        lambda: q.IsAromaticQueryAtom(1, 2),
        lambda: q.AtomNumEqualsQueryAtom(6, val=6),
        lambda: q.HasPropQueryAtom(3),
    ):
        with pytest.raises(Exception) as err:
            call()
        assert type(err.value).__name__ == "ArgumentError"
