"""rdChemReactions (expectations checked against RDKit 2026.03.6 and the port)."""

import pytest

from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions

AMIDE = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]"


def smiles_sets(products):
  res = []
  for ps in products:
    row = []
    for m in ps:
      Chem.SanitizeMol(m)
      row.append(Chem.MolToSmiles(m))
    res.append(tuple(row))
  return res


def test_parse_and_counts():
  rxn = rdChemReactions.ReactionFromSmarts(AMIDE)
  assert (rxn.GetNumReactantTemplates(), rxn.GetNumProductTemplates(),
          rxn.GetNumAgentTemplates()) == (2, 1, 0)
  assert rxn.Validate(True) == (0, 0)
  assert rdChemReactions.ReactionToSmarts(rxn) == AMIDE
  assert rxn._getImplicitPropertiesFlag()


def test_run_reactants():
  rxn = rdChemReactions.ReactionFromSmarts(AMIDE)
  ps = rxn.RunReactants((Chem.MolFromSmiles("CC(=O)O"), Chem.MolFromSmiles("NC")))
  assert smiles_sets(ps) == [("CNC(C)=O", )]
  assert rxn.IsInitialized()


def test_run_reactants_thread_count_is_irrelevant():
  rxn = rdChemReactions.ReactionFromSmarts(AMIDE)
  reactants = (Chem.MolFromSmiles("OC(=O)c1ccccc1C(=O)O"), Chem.MolFromSmiles("NCC"))
  try:
    serial = smiles_sets(rxn.RunReactants(reactants, numThreads=1))
    parallel = smiles_sets(rxn.RunReactants(reactants, numThreads=4))
  except Exception as e:  # the reference RDKit has no numThreads argument
    pytest.skip(f"numThreads is an rdoxide extension: {e}")
  assert serial == parallel == [("CCNC(=O)c1ccccc1C(=O)O", )] * 2


def test_chirality_inversion():
  rxn = rdChemReactions.ReactionFromSmarts(
    "[C@:1]([F:2])([Cl:3])([Br:4])>>[C@@:1]([F:2])([Cl:3])([Br:4])")
  assert rdChemReactions.ReactionToSmarts(rxn) == "[F:2][C@:1]([Cl:3])[Br:4]>>[F:2][C@@:1]([Cl:3])[Br:4]"
  ps = rxn.RunReactants((Chem.MolFromSmiles("F[C@](Cl)(Br)C"), ))
  assert smiles_sets(ps) == [("C[C@](F)(Cl)Br", )]


def test_bad_reactant_count():
  rxn = rdChemReactions.ReactionFromSmarts("[C:1][OH]>>[C:1]Cl")
  with pytest.raises(ValueError, match="Number of reactants provided"):
    rxn.RunReactants((Chem.MolFromSmiles("CCO"), Chem.MolFromSmiles("C")))
  with pytest.raises(ValueError, match="None reactants"):
    rxn.RunReactants((None, ))


def test_run_reactant_and_in_place():
  rxn = rdChemReactions.ReactionFromSmarts(AMIDE)
  acid = Chem.MolFromSmiles("CC(=O)O")
  assert [Chem.MolToSmiles(p[0]) for p in rxn.RunReactant(acid, 0)] == ["CC(N)=O"]
  rxn = rdChemReactions.ReactionFromSmarts("[C:1][OH:2]>>[C:1][N:2]")
  m = Chem.MolFromSmiles("CCO")
  assert rxn.RunReactantInPlace(m)
  assert Chem.MolToSmiles(m) == "CCN"
  m = Chem.MolFromSmiles("CCC")
  assert not rxn.RunReactantInPlace(m)
  assert Chem.MolToSmiles(m) == "CCC"


def test_molecule_roles_and_reacting_atoms():
  rxn = rdChemReactions.ReactionFromSmarts(AMIDE)
  rxn.Initialize()
  assert rxn.IsMoleculeReactant(Chem.MolFromSmiles("CC(=O)O"))
  assert rxn.IsMoleculeProduct(Chem.MolFromSmiles("CC(=O)NC"))
  assert not rxn.IsMoleculeAgent(Chem.MolFromSmiles("C"))
  assert rxn.GetReactingAtoms() == ((0, 2), (0, ))


def test_smiles_reactions_and_agents():
  rxn = rdChemReactions.ReactionFromSmiles("CC(=O)O.N>[Pd]>CC(=O)N")
  assert (rxn.GetNumReactantTemplates(), rxn.GetNumAgentTemplates(),
          rxn.GetNumProductTemplates()) == (2, 1, 1)
  assert rdChemReactions.ReactionToSmiles(rxn) == "CC(=O)O.N>[Pd]>CC(N)=O"
  assert Chem.MolToSmiles(rxn.GetAgentTemplate(0)) == "[Pd]"
  with pytest.raises(ValueError, match="template index too high"):
    rxn.GetAgentTemplate(1)


@pytest.mark.parametrize("text", ["CC>", "C>C>C>C"])
def test_malformed_reactions(text):
  with pytest.raises(ValueError, match="ChemicalReactionParserException"):
    rdChemReactions.ReactionFromSmarts(text)


def test_unmapped_templates_and_agents():
  rxn = rdChemReactions.ReactionFromSmarts("[C:1].O>>[C:1]")
  rxn.RemoveUnmappedReactantTemplates()
  assert (rxn.GetNumReactantTemplates(), rxn.GetNumAgentTemplates()) == (1, 1)
  removed = []
  rxn.RemoveAgentTemplates(removed)
  assert len(removed) == 1 and rxn.GetNumAgentTemplates() == 0


def test_properties():
  rxn = rdChemReactions.ReactionFromSmarts("[C:1]>>[C:1]")
  rxn.SetProp("a", "b")
  rxn.SetIntProp("i", 3)
  assert rxn.GetPropsAsDict() == {"a": "b", "i": 3}
  assert rxn.HasProp("a")
  assert rxn.GetIntProp("i") == 3


def test_helpers():
  rxn = rdChemReactions.ReactionFromSmarts(AMIDE)
  assert rdChemReactions.HasReactionAtomMapping(rxn)
  assert rdChemReactions.IsReactionTemplateMoleculeAgent(Chem.MolFromSmarts("[C:1]O"), 0.5) is False
  product = rxn.RunReactants((Chem.MolFromSmiles("CC(=O)O"), Chem.MolFromSmiles("NC")))[0][0]
  side = rdChemReactions.ReduceProductToSideChains(product)
  assert Chem.MolToSmiles(side) == "C[*:1].C[*:3]"
  rdChemReactions.RemoveMappingNumbersFromReactions(rxn)
  assert rdChemReactions.ReactionToSmarts(rxn) == "C(=O)O.N>>C(=O)N"
  mol = rdChemReactions.ReactionToMolecule(rxn)
  assert sorted(a.GetIntProp("molRxnRole") for a in mol.GetAtoms()) == [1, 1, 1, 1, 2, 2, 2]


def test_substruct_params_are_live():
  rxn = rdChemReactions.ReactionFromSmarts("[C@:1]([F:2])([Cl:3])Br>>[C:1]([F:2])([Cl:3])Br")
  rxn.GetSubstructParams().useChirality = True
  rxn.Initialize()
  assert not rxn.IsMoleculeReactant(Chem.MolFromSmiles("F[C@@H](Cl)Br"))
  assert rxn.GetSubstructParams().useChirality


def test_allchem_exports():
  assert AllChem.ReactionFromSmarts is rdChemReactions.ReactionFromSmarts


SMARTS_WITH_AGENT = "[C:1](=[O:2])O.[N:3]>[Pd]>[C:1](=[O:2])[N:3]"
ROUND_TRIP = "[#6:1](=[#8:2])[#8].[#7:3]>[Pd]>[#6:1](=[#8:2])[#7:3]"


def test_rxn_block_v2000():
  rxn = rdChemReactions.ReactionFromSmarts(SMARTS_WITH_AGENT)
  block = rdChemReactions.ReactionToRxnBlock(rxn)
  assert block.startswith("$RXN\n\n      RDKit\n\n  3  1\n$MOL\n")
  assert block.count("$MOL") == 4
  back = rdChemReactions.ReactionFromRxnBlock(rdChemReactions.ReactionToRxnBlock(rxn, True))
  assert (back.GetNumReactantTemplates(), back.GetNumAgentTemplates(),
          back.GetNumProductTemplates()) == (2, 1, 1)
  assert rdChemReactions.ReactionToSmarts(back) == ROUND_TRIP
  assert not back._getImplicitPropertiesFlag()


def test_rxn_block_v3000():
  rxn = rdChemReactions.ReactionFromSmarts(SMARTS_WITH_AGENT)
  block = rdChemReactions.ReactionToV3KRxnBlock(rxn, True)
  assert block.startswith("$RXN V3000\n\n      RDKit\n\nM  V30 COUNTS 2 1 1\nM  V30 BEGIN REACTANT\n")
  assert "M  V30 BEGIN AGENT" in block
  assert rdChemReactions.ReactionToRxnBlock(rxn, False, True).startswith("$RXN V3000")
  back = rdChemReactions.ReactionFromRxnBlock(block)
  assert rdChemReactions.ReactionToSmarts(back) == ROUND_TRIP


@pytest.mark.parametrize("block", ["", "$RXN\n", "x", "$RXN\n\n  RDKit\n\n  1  1\n", "$RXN\n\n  RDKit\n\n  a  1\n$MOL\n"])
def test_bad_rxn_blocks(block):
  with pytest.raises(ValueError, match="ChemicalReactionParserException|ChemicalParserException"):
    rdChemReactions.ReactionFromRxnBlock(block)


def test_rxn_file(tmp_path):
  rxn = rdChemReactions.ReactionFromSmarts(AMIDE)
  path = tmp_path / "amide.rxn"
  path.write_text(rdChemReactions.ReactionToRxnBlock(rxn))
  back = rdChemReactions.ReactionFromRxnFile(str(path))
  assert rdChemReactions.ReactionToSmarts(back) == "[#6:1](=[#8:2])[#8].[#7:3]>>[#6:1](=[#8:2])[#7:3]"
  assert rdChemReactions.ReactionFromRxnFile(str(tmp_path / "missing.rxn")) is None


def test_reaction_molecule_round_trip():
  rxn = rdChemReactions.ReactionFromSmarts(SMARTS_WITH_AGENT)
  mol = rdChemReactions.ReactionToMolecule(rxn)
  back = rdChemReactions.ReactionFromMolecule(mol)
  assert (back.GetNumReactantTemplates(), back.GetNumAgentTemplates(),
          back.GetNumProductTemplates()) == (2, 1, 1)
  assert rdChemReactions.ReactionToSmarts(back) == rdChemReactions.ReactionToSmarts(rxn)


def _bromide_reaction():
  rxn = rdChemReactions.ReactionFromSmarts("[C:1][Br].[N:2]>>[C:1][N:2]")
  rxn.Initialize()
  return rxn


def test_preprocess_reaction_with_default_groups():
  rxn = _bromide_reaction()
  r1 = rxn.GetReactantTemplate(0)
  r1.GetAtomWithIdx(1).SetProp("molFileValue", "Halogen.Aromatic")
  assert r1.HasSubstructMatch(r1)
  assert Chem.MolFromSmiles("CBr").HasSubstructMatch(r1)
  warnings, errors, n_reactants, n_products, labels = rdChemReactions.PreprocessReaction(rxn)
  assert (warnings, errors, n_reactants, n_products) == (0, 0, 2, 1)
  assert labels == (((1, "halogen.aromatic"), ), ())
  # the handed-out template sees the new query, as upstream's reference does
  assert not Chem.MolFromSmiles("CBr").HasSubstructMatch(r1)
  assert not Chem.MolFromSmiles("c1ccccc1Br").HasSubstructMatch(r1)


def test_preprocess_reaction_unknown_group():
  rxn = _bromide_reaction()
  rxn.GetReactantTemplate(0).GetAtomWithIdx(1).SetProp("molFileValue", "boromicacid")
  with pytest.raises(KeyError, match="boromicacid"):
    rdChemReactions.PreprocessReaction(rxn)


def test_add_recursive_queries_to_reaction():
  rxn = _bromide_reaction()
  rxn.GetReactantTemplate(0).GetAtomWithIdx(1).SetProp("query", "ar")
  queries = {"ar": Chem.MolFromSmarts("[Br]-a")}
  assert rxn.AddRecursiveQueriesToReaction(queries, "query") is None
  q = rxn.GetReactantTemplate(0)
  assert Chem.MolFromSmiles("c1ccccc1Br").HasSubstructMatch(Chem.MolFromSmarts("[Br]-a"))
  assert not Chem.MolFromSmiles("CBr").HasSubstructMatch(q)
  rxn = _bromide_reaction()
  rxn.GetReactantTemplate(0).GetAtomWithIdx(1).SetProp("query", "ar")
  labels = rxn.AddRecursiveQueriesToReaction(queries, "query", True)
  assert labels == (((1, "ar"), ), ())


def test_sanitize_flags_and_rxn_params():
  assert rdChemReactions.SANITIZE_NONE == 0
  assert int(rdChemReactions.SanitizeFlags.SANITIZE_MERGEHS) == 8
  assert not rdChemReactions.GetDefaultAdjustParams().adjustDegree
  assert rdChemReactions.MatchOnlyAtRgroupsAdjustParams().adjustDegree


def test_sanitize_rxn():
  rxn = rdChemReactions.ReactionFromRxnBlock(
    rdChemReactions.ReactionToRxnBlock(rdChemReactions.ReactionFromSmarts(AMIDE)))
  assert int(rdChemReactions.SanitizeRxn(rxn, catchErrors=True)) == 0
  rdChemReactions.SanitizeRxnAsMols(rxn, Chem.SanitizeFlags.SANITIZE_NONE)


def test_pickle_round_trip():
  import pickle
  rxn = rdChemReactions.ReactionFromSmarts(SMARTS_WITH_AGENT)
  rxn.Initialize()
  rxn.GetSubstructParams().useChirality = True
  copy = pickle.loads(pickle.dumps(rxn))
  assert rdChemReactions.ReactionToSmarts(copy) == rdChemReactions.ReactionToSmarts(rxn)
  assert copy.IsInitialized() and copy._getImplicitPropertiesFlag()
  assert copy.GetSubstructParams().useChirality
  assert (copy.GetNumReactantTemplates(), copy.GetNumAgentTemplates(),
          copy.GetNumProductTemplates()) == (2, 1, 1)
  from_binary = rdChemReactions.ChemicalReaction(rxn.ToBinary())
  assert rdChemReactions.ReactionToSmarts(from_binary) == rdChemReactions.ReactionToSmarts(rxn)
  ps = copy.RunReactants((Chem.MolFromSmiles("CC(=O)O"), Chem.MolFromSmiles("NC")))
  assert smiles_sets(ps) == [("CNC(C)=O", )]


def test_handed_out_templates_are_live():
  rxn = rdChemReactions.ReactionFromSmarts("[C:1][O:2]>>[C:1][N:2]")
  template = rxn.GetReactantTemplate(0)
  template.GetAtomWithIdx(0).SetProp("edited", "yes")
  rxn.Initialize()
  assert rxn.GetReactantTemplate(0).GetAtomWithIdx(0).GetProp("edited") == "yes"
  assert rxn.GetReactants()[0].GetAtomWithIdx(0).GetProp("edited") == "yes"


def test_cxsmarts_reaction():
  # CXExtensionsTests.test_cxsmarts_reaction
  text = ("[CH3:1][CH:2]([CH3:3])[*:4].[OH:5][CH2:6][*:7]>O=C=O>"
          "[CH3:1][CH:2]([CH3:3])[CH2:6][OH:5] |$;;;_AP1;;;_AP1;;;;;;;;$|")
  rxn = rdChemReactions.ReactionFromSmarts(text)
  flags = Chem.CXSmilesFields.CX_ALL ^ Chem.CXSmilesFields.CX_ATOM_PROPS
  assert rdChemReactions.ReactionToCXSmarts(rxn, Chem.SmilesWriteParams(), flags) == (
    "[C&H3:1][C&H1:2]([C&H3:3])[*:4].[O&H1:5][C&H2:6][*:7]>O=C=O>"
    "[C&H3:1][C&H1:2]([C&H3:3])[C&H2:6][O&H1:5] |$;;;_AP1;;;_AP1;;;;;;;;$|")


def test_cx_extensions_span_all_templates():
  # labels and properties are numbered over reactants, agents, then products
  rxn = rdChemReactions.ReactionFromSmarts(
    "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3] |$_R1;;;;;;$,atomProp:1.foo.bar|")
  assert rxn.GetReactantTemplate(0).GetAtomWithIdx(0).GetProp("atomLabel") == "_R1"
  assert rxn.GetReactantTemplate(0).GetAtomWithIdx(1).GetProp("foo") == "bar"
  assert not rxn.GetReactantTemplate(1).GetAtomWithIdx(0).HasProp("atomLabel")
  assert rdChemReactions.ReactionToCXSmarts(rxn) == (
    "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3] |$_R1;;;;;;$,atomProp:1.foo.bar|")
  assert rdChemReactions.ReactionToCXSmiles(rxn) == (
    "O[C:1]=[O:2].[NH3:3]>>[C:1](=[O:2])[N:3] |$;_R1;;;;;$,atomProp:2.foo.bar|")


def test_cxsmiles_reaction_second_template():
  rxn = rdChemReactions.ReactionFromSmiles("CC.CO>>CCO |$;;;a;;;b$|")
  assert rxn.GetReactantTemplate(1).GetAtomWithIdx(1).GetProp("atomLabel") == "a"
  assert rxn.GetProductTemplate(0).GetAtomWithIdx(2).GetProp("atomLabel") == "b"
