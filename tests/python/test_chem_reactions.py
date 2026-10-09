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
