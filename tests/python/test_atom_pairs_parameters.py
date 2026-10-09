"""rdMolDescriptors.AtomPairsParameters and the decoding helpers built on it."""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.AtomPairs import Pairs, Utils


def test_constants_match_the_reference():
    p = rdMolDescriptors.AtomPairsParameters
    assert (p.version, p.numTypeBits, p.numPiBits, p.numBranchBits) == ("1.1.0", 4, 2, 3)
    assert (p.numChiralBits, p.codeSize, p.numPathBits, p.numAtomPairFingerprintBits) == (2, 9, 5, 23)
    assert list(p.atomTypes) == [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53, 0]


def test_atom_code_round_trip():
    mol = Chem.MolFromSmiles("CC(=O)N")
    code = rdMolDescriptors.GetAtomPairAtomCode(mol.GetAtomWithIdx(3))
    assert Utils.ExplainAtomCode(code) == ("N", 1, 0)


def test_pair_fingerprint_explanation():
    fp = Pairs.GetAtomPairFingerprint(Chem.MolFromSmiles("CCO"))
    assert len(fp.GetNonzeroElements()) == 3
