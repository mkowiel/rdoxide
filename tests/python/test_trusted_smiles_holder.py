"""Molecules of a trusted-SMILES holder are not sanitized: their rings are not perceived."""

import pytest
from rdoxide import Chem
from rdoxide.Chem import rdSubstructLibrary


def test_ring_queries_need_the_library():
  holder = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
  lib = rdSubstructLibrary.SubstructLibrary(holder)
  lib.AddMol(Chem.MolFromSmiles("C1CC1"))
  ring_atom, ring_bond = Chem.MolFromSmarts("[C&R1]"), Chem.MolFromSmarts("C@C")
  for query in (ring_atom, ring_bond):
    with pytest.raises(RuntimeError):
      holder.GetMol(0).HasSubstructMatch(query)
    assert len(lib.GetMatches(query)) == 1
    assert lib.CountMatches(query) == 1
