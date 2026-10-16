"""rdoxide.Chem.rdAbbreviations; expected values come from reference RDKit."""

import pytest
from rdoxide import Chem
from rdoxide.Chem import rdAbbreviations as A


def test_default_tables():
    abbrevs = A.GetDefaultAbbreviations()
    assert len(abbrevs) == 37 and abbrevs[0].label == "CO2Et"
    assert abbrevs[0].displayLabelW == "EtO<sub>2</sub>C"
    assert abbrevs[0].includesXBonds
    assert [a.label for a in A.GetDefaultLinkers()][:2] == ["PEG6", "PEG5"]


def test_condense_labels_atoms_and_tracks_mapping():
    m = Chem.MolFromSmiles("CCOC(=O)c1ccc(C(F)(F)F)cc1")
    res = A.CondenseMolAbbreviations(m, A.GetDefaultAbbreviations(), maxCoverage=0.8)
    labels = [a.GetProp("atomLabel") for a in res.GetAtoms() if a.HasProp("atomLabel")]
    assert labels == ["CO2Et", "CF3"]
    assert res.GetNumAtoms() == 8
    assert res.GetProp("_origAtomMapping") == "[3,5,6,7,8,9,13,14]"


def test_label_then_condense_substance_groups():
    m = Chem.MolFromSmiles("c1ccccc1C(F)(F)F")
    labelled = A.LabelMolAbbreviations(m, A.GetDefaultAbbreviations(), maxCoverage=0.8)
    groups = Chem.GetMolSubstanceGroups(labelled)
    assert [g.GetProp("LABEL") for g in groups] == ["CF3"]
    condensed = A.CondenseAbbreviationSubstanceGroups(labelled)
    assert condensed.GetNumAtoms() == 7


def test_parse_abbreviations_and_linkers():
    defs = A.ParseAbbreviations("FOO C(F)(F)F Foo FooW\nBAR OC\n")
    assert [(d.label, d.displayLabel, d.displayLabelW) for d in defs] == [
        ("FOO", "Foo", "FooW"),
        ("BAR", "", ""),
    ]
    assert [d.label for d in A.ParseLinkers("L1 *CC*\nL2 *OCC*")] == ["L1", "L2"]


def test_max_coverage_skips_large_abbreviations():
    m = Chem.MolFromSmiles("CC(F)(F)F")
    res = A.CondenseMolAbbreviations(m, A.GetDefaultAbbreviations(), maxCoverage=0.4)
    assert res.GetNumAtoms() == m.GetNumAtoms()


def test_missing_substance_groups_is_a_no_op():
    m = Chem.MolFromSmiles("CCO")
    assert A.CondenseAbbreviationSubstanceGroups(m).GetNumAtoms() == 3
