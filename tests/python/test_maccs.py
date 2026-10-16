"""MACCS keys (expectations from RDKit 2026.09.1)."""

from rdoxide import Chem
from rdoxide.Chem import MACCSkeys, rdMolDescriptors


def test_maccs():
    m = Chem.MolFromSmiles("OC(=O)c1ccccc1N.Cl")
    fp = MACCSkeys.GenMACCSKeys(m)
    assert fp.GetNumBits() == 167
    assert list(fp.GetOnBits()) == [53, 84, 97, 103, 123, 131, 133, 134, 135, 139, 144, 150, 151, 154, 156, 157, 158, 159, 161, 162, 163, 164, 165, 166]
    assert rdMolDescriptors.GetMACCSKeysFingerprint(m) == fp
    assert MACCSkeys.GenMACCSKeys(Chem.Mol()).GetNumOnBits() == 0
