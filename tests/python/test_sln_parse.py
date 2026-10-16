"""rdoxide.Chem.rdSLNParse (beyond the upstream testSLN.py)."""
import unittest

from rdoxide import Chem
from rdoxide.Chem import rdSLNParse


class TestSLNParse(unittest.TestCase):

  def testMols(self):
    m = rdSLNParse.MolFromSLN('CH3CH(CH3)OCH3')
    self.assertEqual(Chem.MolToSmiles(m), 'COC(C)C')
    m = rdSLNParse.MolFromSLN('CH3ZCH3{Z:O}')
    self.assertEqual(Chem.MolToSmiles(m), 'COC')

  def testProps(self):
    m = rdSLNParse.MolFromSLN('CH4<name=methane;blah;foo="1">')
    self.assertEqual(m.GetProp('_Name'), 'methane')
    self.assertTrue(m.HasProp('blah'))
    self.assertEqual(m.GetProp('foo'), '1')

  def testFailures(self):
    self.assertIsNone(rdSLNParse.MolFromSLN('CH4[I=13]'))
    self.assertIsNone(rdSLNParse.MolFromSLN('C1CC'))

  def testQueries(self):
    patt = rdSLNParse.MolFromQuerySLN('Any[charge=+1;!HC=2]')
    m = Chem.MolFromSmiles('C[CH2+](CC)[NH+]', sanitize=False)
    self.assertEqual(m.GetSubstructMatches(patt), ((4, ), ))
    patt = rdSLNParse.MolFromQuerySLN('C[IS=O]')
    self.assertEqual(patt.GetNumAtoms(), 1)


if __name__ == "__main__":
  unittest.main()
