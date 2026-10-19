import pytest
from rdoxide import Chem


def test_set_expand_query_none_raises_value_error():
    m = Chem.MolFromSmarts("CC")
    for obj in (m.GetAtomWithIdx(0), m.GetBondWithIdx(0)):
        with pytest.raises(ValueError):
            obj.SetQuery(None)
        with pytest.raises(ValueError):
            obj.ExpandQuery(None)
