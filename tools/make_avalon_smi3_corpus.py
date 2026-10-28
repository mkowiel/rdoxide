#!/usr/bin/env python3
"""Builds tests/data/avalon_smi2mol.jsonl: SMILES and the molecule SMIToMOL of the Avalon C code
returns for flags 0, TRUE_DB_STEREO and EXPECT_SMARTS (driver: .work/avalon_drv/smi3.c, no layout).

    .venv-ref/bin/python tools/make_avalon_smi3_corpus.py DRIVER > tests/data/avalon_smi2mol.jsonl
"""
import json
import subprocess
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

SMILES = [
    "CCO", "c1ccccc1O", "C[N+](C)(C)C", "[O-]C(=O)C", "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "C/C=C/C", "C/C=C\\C", "[13CH3]C",
    "[2H]C([2H])Cl", "[CH3]", "[Na+].[Cl-]", "C1CC1", "O=S(=O)(N)c1ccccc1", "[Fe+3]", "CC(C)(C)c1ccc(O)cc1", "C#CC", "[NH4+]",
    "C[Si](C)(C)C", "FC(F)(F)C(=O)O", "[14C]", "[18F]C", "OC[C@@H](O)[C@H](O)CO", "[CH2]C", "c1ccc2ccccc2c1", "[Se]", "[3H]O",
    "[N+](=O)[O-]", "C12C3C4C1C5C2C3C45", "CC12CCC3C(CCc4cc(O)ccc34)C1CCC2O", "c1cc[nH]c1", "c1ccoc1", "c1ccsc1", "c1ccncc1",
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "O=C1C=CC(=O)C=C1", "F/C=C/C=C/Cl", "F/C=C\\C=C/Cl", "N[C@@H](C)C(=O)N[C@H](C(=O)O)Cc1c[nH]c2ccccc12",
    "C[C@H]1CC[C@@H](C)CC1", "CC(=O)Oc1ccccc1C(=O)O", "C=C=C", "CC=C=CC", "C1CCC1C1CCC1", "c1ccc(cc1)-c1ccccc1", "[Cl-].[Cl-].[Mg+2]",
    "CC[S+](C)C", "C[P](C)(C)=O", "CS(=O)C", "O=C1CCCN1", "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O", "[C@H](F)(Cl)Br",
    "C[C@](F)(Cl)Br", "N#Cc1ccccc1", "[H][H]", "O=C(c1ccccc1)c1ccccc1", "CCN(CC)CC", "NC(=N)N", "C1CC2CCC1C2", "CN1CCC[C@H]1c1cccnc1",
    "[O-][n+]1ccccc1", "c1ccc2[nH]ccc2c1", "C(/C=C/c1ccccc1)=O", "[nH]1cccc1", "O=s1cccc1", "[O-][N+](=O)c1ccccc1", "C[N+]([O-])=O",
    "[1*]CC", "[*]CC(=O)O", "[*:3]C", "[R]C", "C[Mg]Br", "F[C@](Cl)(Br)I", "C[C@@H]1CC[C@H](C)CC1", "CC1=CC=CC=C1", "C1=CC=CC=C1",
    "N[C@@H](CS)C(=O)O", "C1CC1.C1CC1", "C%11CC%11", "C=1CC=1", "C-1CC1", "C#1CC#1", "[C@@H]1(F)CCCC1", "[C@H]1(F)CCCC1Cl",
]
EXTRA = [
    "{abc}C", "C{x}", "{}C", "{ab", "C.1CC.1", "C(C)(C)(C)(C)C", "C(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)",
    "C(F)(Cl)(Br)(I)", "[C@](F)(Cl)(Br)I", "[C@@](F)(Cl)(Br)(I)N", "F/C=C/F", "F\\C=C\\F", "F/C=C(/F)Cl", "C/C=C(/C)\\C", "FC=CF",
    "F/C=C/C=C/F", "c1ccccc1/C=C/c1ccccc1", "C/C=C/C=C\\C", "[CH-]", "[CH+]", "[C-]", "O[SH]", "[PH4+]", "[OH3+]", "[nH+]1ccccc1", "C=[N+]=[N-]",
    "OP(=O)(O)O", "S(=O)(=O)(O)O", "Cl[Cl]", "[Xe]", "[Uuo]", "[U]", "CC(=O)[O-].[Na+]", "n1ccccc1", "c1ccccc1c1ccccc1", "c1ccccc1-c1ccccc1",
    "c1cc2ccc3cccc4ccc(c1)c2c34", "C1=CC=C[CH]1", "[CH]1=CC=CC=C1", "[cH-]1cccc1", "[nH]1cccc1.C", "Cl.C", "C.", "C..C", "C(C)", "CC)", "(C)C",
]


def main():
    smis = []
    for smi in SMILES + EXTRA:
        smis.append(smi)
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        smis.append(Chem.MolToSmiles(m))
        smis.append(Chem.MolToSmiles(m, kekuleSmiles=True) if not Chem.MolToSmiles(m, kekuleSmiles=True) is None else smi)
        for _ in range(3):
            smis.append(Chem.MolToSmiles(m, doRandom=True))
    seen, uniq = set(), []
    for s in smis:
        if s not in seen and "\n" not in s:
            seen.add(s)
            uniq.append(s)
    out = subprocess.run([sys.argv[1]], input=("\n".join(uniq) + "\n").encode(), capture_output=True, timeout=300).stdout.decode()
    parts = out.split("@@\n")
    assert len(parts) == 3 * len(uniq) + 1, (len(parts), len(uniq))
    for i, s in enumerate(uniq):
        print(json.dumps({"input": s, "outputs": parts[3 * i : 3 * i + 3]}))


if __name__ == "__main__":
    main()
