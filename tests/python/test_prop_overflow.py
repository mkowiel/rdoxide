import pytest
from rdoxide import Chem


def test_integer_overflow_messages():
    m = Chem.MolFromSmiles("CC")
    for ob in (m, m.GetAtomWithIdx(0), m.GetBondWithIdx(0)):
        ob.SetIntProp("foo", -1)
        with pytest.raises(ValueError, match="negative overflow"):
            ob.GetUnsignedProp("foo")
        ob.SetUnsignedProp("foo", 4294967295)
        assert ob.GetUnsignedProp("foo") == 4294967295
        with pytest.raises(ValueError, match="positive overflow"):
            ob.GetIntProp("foo")
        ob.SetDoubleProp("foo", 2.0)
        with pytest.raises(ValueError, match="bad any cast"):
            ob.GetIntProp("foo")


def test_none_property_name_is_an_argument_error():
    m = Chem.Mol()
    for name, extra in [("HasProp", ()), ("GetProp", ()), ("SetProp", ("",)), ("SetIntProp", (1,))]:
        with pytest.raises(TypeError, match="Python argument types"):
            getattr(m, name)(None, *extra)
