import pytest

from rdoxide import Chem
from rdoxide.Chem import rdChemReactions as rc


def test_initialize_uses_only_the_sizes_and_validates_molecules():
    rxn = rc.ChemicalReaction()
    groups = [[Chem.MolFromSmiles("C")] * 3, [Chem.MolFromSmiles("N")] * 2]
    s = rc.CartesianProductStrategy()
    s.Initialize(rxn, groups)
    assert s.GetNumPermutations() == 6
    assert [tuple(next(s)) for _ in range(2)] == [(0, 0), (1, 0)]
    with pytest.raises(ValueError, match="non molecule"):
        s.Initialize(rxn, [[1]])
