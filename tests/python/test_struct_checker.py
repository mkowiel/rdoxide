"""rdoxide.Chem.rdStructChecker (beyond the upstream rough_test.py)."""
import unittest

from rdoxide import Chem
from rdoxide.Chem import rdStructChecker


class TestStructChecker(unittest.TestCase):

  def testFlagStrings(self):
    flags = rdStructChecker.StructureFlags
    s = rdStructChecker.StructChecker.StructureFlagsToString(flags.STEREO_ERROR | flags.TRANSFORMED)
    self.assertEqual(s, "STEREO_ERROR,TRANSFORMED")
    self.assertEqual(rdStructChecker.StructChecker.StringToStructureFlags(s),
                     flags.STEREO_ERROR | flags.TRANSFORMED)

  def testDefaults(self):
    chk = rdStructChecker.StructChecker()
    m = Chem.MolFromSmiles("COC(=O)C")
    self.assertEqual(chk.CheckMolStructure(m), rdStructChecker.StructureFlags.NO_CHANGE)

  def testOptions(self):
    opts = rdStructChecker.StructCheckerOptions()
    self.assertEqual(opts.MaxMolSize, 255)
    self.assertEqual(opts.CollisionLimitPercent, 15)
    opts.MaxMolSize = 3
    chk = rdStructChecker.StructChecker(opts)
    m = Chem.MolFromSmiles("CCCCC")
    self.assertEqual(chk.CheckMolStructure(m), rdStructChecker.StructureFlags.SIZE_CHECK_FAILED)
    self.assertFalse(opts.LoadGoodAugmentedAtoms("/does/not/exist.chk"))


if __name__ == "__main__":
  unittest.main()
