import pytest

from rdoxide import Chem
from rdoxide.Chem import AllChem

ETHANOL = """<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
<CDXML BondLength="14.4">
<page id="1">
<fragment id="2">
<n id="3" p="100 100" Element="6"/>
<n id="4" p="114.4 100" Element="6"/>
<n id="5" p="128.8 100" Element="8"/>
<b id="6" B="3" E="4"/>
<b id="7" B="4" E="5"/>
</fragment>
</page>
</CDXML>"""


def test_molecules_from_a_cdxml_block():
    mols = Chem.MolsFromCDXML(ETHANOL)
    assert len(mols) == 1
    assert Chem.MolToSmiles(mols[0]) == "CCO"
    conf = mols[0].GetConformer()
    # coordinates are rescaled to the RDKit bond length
    assert conf.GetAtomPosition(0).Distance(conf.GetAtomPosition(1)) == pytest.approx(1.5)


def test_params_and_options():
    params = Chem.CDXMLParserParams()
    assert params.sanitize and params.removeHs and not params.parseQueries
    params.sanitize = False
    assert len(Chem.MolsFromCDXML(ETHANOL, params)) == 1
    assert len(Chem.MolsFromCDXML(ETHANOL, False, False)) == 1
    with pytest.raises(AttributeError):
        params.notAnAttribute = 1
    assert Chem.HasChemDrawCDXSupport() is False
    with pytest.raises(RuntimeError):
        Chem.MolsFromCDXML(ETHANOL, Chem.CDXMLParserParams(True, True, Chem.CDXMLFormat.CDX))


def test_bad_input():
    with pytest.raises(RuntimeError):
        Chem.MolsFromCDXML("<CDXML")
    with pytest.raises(OSError):
        Chem.MolsFromCDXMLFile("/does/not/exist.cdxml")
    assert AllChem.ReactionsFromCDXMLBlock("") == ()
    assert AllChem.ReactionsFromCDXMLFile("/does/not/exist.cdxml") == ()
    with pytest.raises(RuntimeError):
        Chem.MolToCDXMLBlock(Chem.MolFromSmiles("C"))
