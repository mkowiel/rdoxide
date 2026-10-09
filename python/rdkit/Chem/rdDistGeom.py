#
# rdoxide: PROVISIONAL placeholder for upstream's native rdkit.Chem.rdDistGeom.
#
# Distance geometry (embedding, ETKDG, bounds matrices) is not ported yet (see
# docs/STATUS.md). The module exists so that pure-Python modules which import
# it at the top level (EnumerateStereoisomers) can be used; everything that
# needs an embedding raises NotImplementedError. It is replaced by the native
# module when distance geometry is ported.
__all__ = ['EmbedMolecule']


def EmbedMolecule(mol, *args, **kwargs):
  """Placeholder: distance geometry embedding is not ported yet."""
  raise NotImplementedError("rdkit.Chem.rdDistGeom.EmbedMolecule is not ported to rdoxide yet")
