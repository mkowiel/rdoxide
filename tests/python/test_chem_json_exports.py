from rdoxide import Chem


def test_json_functions_are_exported_from_chem():
    for name in ("MolToJSON", "MolsToJSON", "JSONToMols", "JSONParseParameters", "JSONWriteParameters"):
        assert hasattr(Chem, name), name
    mols = Chem.JSONToMols(Chem.MolToJSON(Chem.MolFromSmiles("CCO")))
    assert [Chem.MolToSmiles(m) for m in mols] == ["CCO"]
