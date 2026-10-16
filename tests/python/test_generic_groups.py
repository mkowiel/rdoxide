"""Generic group matching; cases from upstream GenericGroups/generic_tests.cpp."""
import pytest
from rdoxide import Chem

CASES = [
    (['Alkyl', 'ALK'], [('O=CCC', 1), ('O=CC=C', 0), ('O=CCC=C', 0), ('O=CCO', 0), ('O=CC', 1), ('O=CC[H]', 1), ('O=C[H]', 0), ('O=CC.CC(C=O)C', 2), ('O=CCC1CC1', 0)]),
    (['AlkylH', 'ALH'], [('O=CCC', 1), ('O=CC=C', 0), ('O=CCC=C', 0), ('O=CCO', 0), ('O=CC', 1), ('O=CC[H]', 1), ('O=C[H]', 1), ('O=CC.CC(C=O)C', 2), ('O=CCC1CC1', 0)]),
    (['Carboacyclic', 'ABC'], [('O=CC1CC1', 0), ('O=CCCCC1CCC1', 0), ('O=CCCCC', 1), ('O=CC([H])([H])([H])', 1), ('O=CCCCO', 0), ('O=C[H]', 0)]),
    (['CarboacyclicH', 'ABH'], [('O=CC1CC1', 0), ('O=CCCCC1CCC1', 0), ('O=CCCCC', 1), ('O=CC([H])([H])([H])', 1), ('O=CCCCO', 0), ('O=C[H]', 1)]),
    (['Acyclic', 'ACY'], [('O=CC1CC1', 0), ('O=CCCCC1CCC1', 0), ('O=CCCCC', 1), ('O=CCCCO', 1), ('O=C[H]', 0)]),
    (['AcyclicH', 'ACH'], [('O=CC1CC1', 0), ('O=CCCCC1CCC1', 0), ('O=CCCCC', 1), ('O=CCCCO', 1), ('O=C[H]', 1)]),
    (['Alkenyl', 'AEL'], [('O=CCC', 0), ('O=CC=C', 1), ('O=CCC=C', 1), ('O=CCC#C', 0), ('O=CC=N', 0), ('O=CC=C[H]', 1), ('O=CCC1CC=C1', 0), ('O=C[H]', 0)]),
    (['AlkenylH', 'AEH'], [('O=CCC', 0), ('O=CC=C', 1), ('O=CCC=C', 1), ('O=CCC#C', 0), ('O=CC=N', 0), ('O=CC=C[H]', 1), ('O=CCC1CC=C1', 0), ('O=C[H]', 1)]),
    (['Alkynyl', 'AYL'], [('O=CCC', 0), ('O=CC#C', 1), ('O=CCC#C', 1), ('O=CCC=C', 0), ('O=CC#N', 0), ('O=CC#C[H]', 1), ('O=CCC1CC#C1', 0), ('O=C([H])[H]', 0)]),
    (['AlkynylH', 'AYH'], [('O=CCC', 0), ('O=CC#C', 1), ('O=CCC#C', 1), ('O=CCC=C', 0), ('O=CC#N', 0), ('O=CC#C[H]', 1), ('O=CCC1CC#C1', 0), ('O=C([H])[H]', 2)]),
    (['Heteroacyclic', 'AHC'], [('O=CC1OC1', 0), ('O=CCCCC1CCC1', 0), ('O=CCCCC', 0), ('O=CCCCO', 1), ('O=C[H]', 0)]),
    (['HeteroacyclicH', 'AHH'], [('O=CC1OC1', 0), ('O=CCCCC1CCC1', 0), ('O=CCCCC', 0), ('O=CCCCO', 1), ('O=C[H]', 1)]),
    (['Alkoxy', 'AOX'], [('O=CCC', 0), ('O=COC', 1), ('O=COC=C', 0), ('O=COC=C', 0), ('O=COCO', 0), ('O=CO', 0), ('O=CO[H]', 0), ('O=COC1CC1', 0), ('O=C([H])[H]', 0)]),
    (['AlkoxyH', 'AOH'], [('O=CCC', 0), ('O=COC', 1), ('O=COC=C', 0), ('O=COC=C', 0), ('O=COCO', 0), ('O=CO', 0), ('O=CO[H]', 0), ('O=COC1CC1', 0), ('O=C([H])[H]', 2)]),
    (['Carboaryl', 'ARY'], [('O=CC1CCCCC1', 0), ('O=Cc1ccccc1', 1), ('O=CC1=CC=CC2=C1CCCC2', 0), ('O=Cc1cccnc1', 0), ('O=CC1=CC=CC2=C1CCNC2', 0), ('O=C([H])[H]', 0)]),
    (['CarboarylH', 'ARH'], [('O=CC1CCCCC1', 0), ('O=Cc1ccccc1', 1), ('O=CC1=CC=CC2=C1CCCC2', 0), ('O=Cc1cccnc1', 0), ('O=CC1=CC=CC2=C1CCNC2', 0), ('O=C([H])[H]', 2)]),
    (['Carbocycloalkyl', 'CAL'], [('O=CCC', 0), ('O=CC1CC1', 1), ('O=CC1C(O)C1', 1), ('O=CCC1CC1', 0), ('O=CC1C=C1', 0), ('O=CC1OC1', 0), ('O=CC1CC(=O)C1', 1), ('O=CC1CC2CCC1CC2', 1), ('O=CC1CC2CCC1NC2', 0), ('O=CC1CCCC2CC=CCC12', 0), ('O=C([H])[H]', 0)]),
    (['CarbocycloalkylH', 'CAH'], [('O=CCC', 0), ('O=CC1CC1', 1), ('O=CC1C(O)C1', 1), ('O=CCC1CC1', 0), ('O=CC1C=C1', 0), ('O=CC1OC1', 0), ('O=CC1CC(=O)C1', 1), ('O=CC1CC2CCC1CC2', 1), ('O=CC1CC2CCC1NC2', 0), ('O=CC1CCCC2CC=CCC12', 0), ('O=C([H])[H]', 2)]),
    (['Carbocyclic', 'CBC'], [('O=CC1CCCCC1', 1), ('O=Cc1ccccc1', 1), ('O=Cc1ccc(O)cc1', 1), ('O=CC1=CC=CC2=C1CCCC2', 1), ('O=Cc1cccnc1', 0), ('O=CC1=CC=CC2=C1CCNC2', 0), ('O=C([H])[H]', 0)]),
    (['CarbocyclicH', 'CBH'], [('O=CC1CCCCC1', 1), ('O=Cc1ccccc1', 1), ('O=Cc1ccc(O)cc1', 1), ('O=CC1=CC=CC2=C1CCCC2', 1), ('O=Cc1cccnc1', 0), ('O=CC1=CC=CC2=C1CCNC2', 0), ('O=C([H])[H]', 2)]),
    (['Carbocycloalkenyl', 'CEL'], [('O=CCC', 0), ('O=Cc1ccccc1', 1), ('O=CC1CC=CCC1', 1), ('O=CC1CC2CCC1CC2', 0), ('O=CC1CC2CCC1C=C2', 0), ('O=CC1CC2C=CC1NC2', 0), ('O=CC1CCCC2=C1CCCC2', 1), ('O=CC1CC=CC2C1C=CC1CCC=CC21', 1), ('O=CC1CC=CC2C1C=CC1CNC=CC21', 0), ('O=C([H])[H]', 0)]),
    (['CarbocycloalkenylH', 'CEH'], [('O=CCC', 0), ('O=Cc1ccccc1', 1), ('O=CC1CC=CCC1', 1), ('O=CC1CC2CCC1CC2', 0), ('O=CC1CC2CCC1C=C2', 0), ('O=CC1CC2C=CC1NC2', 0), ('O=CC1CCCC2=C1CCCC2', 1), ('O=CC1CC=CC2C1C=CC1CCC=CC21', 1), ('O=CC1CC=CC2C1C=CC1CNC=CC21', 0), ('O=C([H])[H]', 2)]),
    (['Heterocyclic', 'CHC'], [('O=CC1OC1', 1), ('O=CC1CC1', 0), ('O=CCC1COC1', 0), ('O=CC1CCC2C1CNC2', 1), ('O=Cc1ccccc1', 0), ('O=Cc1cccnc1', 1), ('O=C([H])[H]', 0)]),
    (['HeterocyclicH', 'CHH'], [('O=CC1OC1', 1), ('O=CC1CC1', 0), ('O=CCC1COC1', 0), ('O=CC1CCC2C1CNC2', 1), ('O=Cc1ccccc1', 0), ('O=Cc1cccnc1', 1), ('O=C([H])[H]', 2)]),
    (['NoCarbonRing', 'CXX'], [('O=CN1NC1', 0), ('O=CN1NN1', 1), ('O=CCN1NNN1', 0), ('O=CN1NNNN2N1N2', 1), ('O=CN1NNNN2N1C2', 0), ('O=CNNN', 0), ('O=C1NNN1', 0), ('O=C([H])[H]', 0)]),
    (['NoCarbonRingH', 'CXH'], [('O=CN1NC1', 0), ('O=CN1NN1', 1), ('O=CCN1NNN1', 0), ('O=CN1NNNN2N1N2', 1), ('O=CN1NNNN2N1C2', 0), ('O=CNNN', 0), ('O=C1NNN1', 0), ('O=C([H])[H]', 2)]),
    (['Cyclic', 'CYC'], [('O=CC1CCC1', 1), ('O=C1CCCC1', 0), ('O=CCC1CCC1', 0), ('O=CC1CNC1', 1), ('O=Cc1ccccc1', 1), ('O=CC1CC2C1C2', 1), ('O=CC1CC2CCC1CC2', 1), ('O=C([H])[H]', 0)]),
    (['CyclicH', 'CYH'], [('O=CC1CCC1', 1), ('O=C1CCCC1', 0), ('O=CCC1CCC1', 0), ('O=CC1CNC1', 1), ('O=Cc1ccccc1', 1), ('O=CC1CC2C1C2', 1), ('O=CC1CC2CCC1CC2', 1), ('O=C([H])[H]', 2)]),
    (['Heteroaryl', 'HAR'], [('O=CC1OC1', 0), ('O=CCC1COC1', 0), ('O=CC1CCC2C1CNC2', 0), ('O=Cc1cccc2c1ccnc2', 1), ('O=Cc1ccccc1', 0), ('O=Cc1cccnc1', 1), ('O=C([H])[H]', 0)]),
    (['HeteroarylH', 'HAH'], [('O=CC1OC1', 0), ('O=CCC1COC1', 0), ('O=CC1CCC2C1CNC2', 0), ('O=Cc1cccc2c1ccnc2', 1), ('O=Cc1ccccc1', 0), ('O=Cc1cccnc1', 1), ('O=C([H])[H]', 2)]),
    (['Group', 'G'], [('O=CC1CCC1', 1), ('O=CCCCC', 1), ('O=CCC1CCC1', 1), ('O=CC1CNC1', 1), ('O=Cc1ccccc1', 1), ('O=CC1CC2C1C2', 1), ('O=CC1CC2CCC1CC2', 1), ('O=C([H])[H]', 0)]),
    (['GroupH', 'GH'], [('O=CC1CCC1', 1), ('O=CCCCC', 1), ('O=CCC1CCC1', 1), ('O=CC1CNC1', 1), ('O=Cc1ccccc1', 1), ('O=CC1CC2C1C2', 1), ('O=CC1CC2CCC1CC2', 1), ('O=C([H])[H]', 2)]),
    (['Group*', 'G*'], [('O=CCCCC', 0), ('O=CCC1CCC1', 1), ('O=CC1CNC1', 1), ('O=Cc1ccccc1', 1), ('O=CC1CC2C1C2', 1), ('O=CC1CC2CCC1CC2', 1), ('O=C([H])[H]', 0)]),
    (['GroupH*', 'GH*'], [('O=CC1CCC1', 1), ('O=CCCCC', 0), ('O=CCC1CCC1', 1), ('O=CC1CNC1', 1), ('O=Cc1ccccc1', 1), ('O=CC1CC2C1C2', 1), ('O=CC1CC2CCC1CC2', 1), ('O=C([H])[H]', 2)]),
]


def _query(label, via):
    if via == "cxsmiles":
        q = Chem.MolFromSmarts(f"*C=O |${label};;$|")
    else:
        q = Chem.MolFromSmarts("*C=O |$;;$|")
        q.GetAtomWithIdx(0).SetProp("atomLabel", label)
    Chem.SetGenericQueriesFromProperties(q, True, True)
    return q


def _target(smi):
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    return Chem.MolFromSmiles(smi, ps)


@pytest.mark.parametrize("via", ["cxsmiles", "atomLabel"])
@pytest.mark.parametrize("labels,tests", CASES)
def test_generic_groups(labels, tests, via):
    ps = Chem.SubstructMatchParameters()
    ps.useGenericMatchers = True
    for label in labels:
        q = _query(label, via)
        assert q.GetAtomWithIdx(0).HasProp("_QueryAtomGenericLabel"), label
        for smi, count in tests:
            assert len(_target(smi).GetSubstructMatches(q, ps)) == count, (label, smi)


def test_generic_label_roundtrip_and_sgroups():
    m = Chem.MolFromSmiles("CO* |$;foo_p;ARY_p$|")
    Chem.SetGenericQueriesFromProperties(m)
    a = m.GetAtomWithIdx(2)
    assert a.GetProp("_QueryAtomGenericLabel") == "ARY"
    assert not a.HasProp("atomLabel")
    assert m.GetAtomWithIdx(1).HasProp("atomLabel")
    Chem.ConvertGenericQueriesToSubstanceGroups(m)
    assert not a.HasProp("_QueryAtomGenericLabel")
    sgs = Chem.GetMolSubstanceGroups(m)
    assert [(s.GetProp("TYPE"), s.GetProp("LABEL")) for s in sgs] == [("SUP", "ARY")]
    Chem.SetGenericQueriesFromProperties(m, useAtomLabels=False)
    assert a.GetProp("_QueryAtomGenericLabel") == "ARY"
    assert len(Chem.GetMolSubstanceGroups(m)) == 0


def test_adjust_with_generic_groups():
    m = Chem.MolFromSmiles("O=C* |$;;ALK_p$|")
    q = Chem.AdjustQueryPropertiesWithGenericGroups(m)
    ps = Chem.SubstructMatchParameters()
    ps.useGenericMatchers = True
    assert q.GetAtomWithIdx(2).GetProp("_QueryAtomGenericLabel") == "ALK"
    assert len(Chem.MolFromSmiles("O=CCC").GetSubstructMatches(q, ps)) == 1
    assert len(Chem.MolFromSmiles("O=CC=C").GetSubstructMatches(q, ps)) == 0
