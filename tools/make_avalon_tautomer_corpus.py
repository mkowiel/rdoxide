#!/usr/bin/env python3
"""Builds tests/data/avalon_tautomer.jsonl and the rule / profile files in tests/data/avalon_patterns/: molecules
sent through Avalon's struchk with tautomer rules (-tt) and option profiles (-f) and what the C code answers (driver:
.work/avalon_drv/sc.c). `{P}` in the options stands for the directory with the files.

    .venv-ref/bin/python tools/make_avalon_tautomer_corpus.py DRIVER PATTERNDIR > tests/data/avalon_tautomer.jsonl
"""
import json
import os
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

RULES = {
    "rule_amide.rdf": ("[C:1](=[O:2])[N:3]", "[C:1]([O:2])=[N:3]"),
    "rule_nitro.rdf": ("[O-:1][N+:2]=[O:3]", "[O:1]=[N+:2][O-:3]"),
    "rule_enol.rdf": ("[C:1]([O:2])=[C:3]", "[C:1](=[O:2])[C:3]"),
}
PROFILES = {
    "profile_a.txt": "# a profile\ntm\ncc   # check collisions\ncl 3\ncn \"999\"\ncs\n\n",
    "profile_b.txt": "cn 4\nda\n# end\n",
    "profile_c.txt": "cs\nzz 1\ncn 9\n",
    "profile_d.txt": "tm\n\n\t  cs  \n# only comments and blanks\n",
}
OPTS = [
    "-tt {P}/rule_amide.rdf", "-tt {P}/rule_nitro.rdf", "-tt {P}/rule_enol.rdf", "-f {P}/profile_a.txt",
    "-f {P}/profile_b.txt|-tm", "-f {P}/profile_c.txt", "-f {P}/profile_d.txt", "-cs|-f {P}/profile_a.txt|-cn 8",
    "-tm|-tt {P}/rule_amide.rdf|-tt {P}/rule_nitro.rdf", "-tt {P}/rule_enol.rdf|-tt {P}/rule_amide.rdf",
]
SMILES = [
    "CC(=O)N", "CC(=O)NC", "CC(=O)N(C)C", "O=C(N)c1ccccc1", "[CH3][C:1](=[O:2])[NH2:3]", "C[C:1](=[O:2])[NH:3]C",
    "[O-:1][N+:2](=[O:3])c1ccccc1", "C[N+](=O)[O-]", "[O-][N+](=O)c1ccc(cc1)C(=O)N", "[CH3:5][C:1](=[O:2])[NH2:3]",
    "OC=C", "O[C:1]=[C:3]C", "C/C=C/O", "CC(O)=CC", "[O:2][C:1]=[C:3]", "O=C(N)C(=O)N", "NC(=O)CC(=O)N",
    "CC(=O)N.CC(=O)N", "CC(C)C(=O)N", "c1ccccc1", "CCO", "N#CC", "CC(=O)O", "[O-:1][N+:2](=[O:3])C", "[O-:1][N+:2](=[O:3])[O-]",
    "O=C1NC(=O)C=C1", "O=C1[NH:3][C:1](=[O:2])C=C1", "CC(=O)[N-:3]C",
]


def rule_file(path, frm, to):
    def block(smarts):
        m = Chem.MolFromSmiles(smarts, sanitize=False)
        m.UpdatePropertyCache(strict=False)
        AllChem.Compute2DCoords(m)
        return Chem.MolToMolBlock(m, kekulize=False)

    with open(path, "w") as f:
        f.write("$RDFILE 1\n$DATM 01/01/26 00:00\n$RFMT $RIREG 1\n$RXN\n\n  rule\n\n  1  1\n$MOL\n" + block(frm) + "$MOL\n" + block(to) + "$$$$\n")


def main():
    driver, pdir = sys.argv[1], sys.argv[2]
    os.makedirs(pdir, exist_ok=True)
    for name, (a, b) in RULES.items():
        rule_file(os.path.join(pdir, name), a, b)
    for name, text in PROFILES.items():
        open(os.path.join(pdir, name), "w").write(text)
    inputs = [("SMILES:" + s, s) for s in SMILES]
    for s in SMILES[::2]:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        inputs.append((Chem.MolToMolBlock(m), "2d " + s))
    recs = []
    for k, (body, name) in enumerate(inputs):
        for j in (k % len(OPTS), (k + 3) % len(OPTS)):
            recs.append((OPTS[j], body, name))
    for o, b, n in recs:
        inp = o.replace("{P}", pdir) + "\n" + b + ("\n" if not b.endswith("\n") else "")
        out = subprocess.run([driver], input=inp.encode(), capture_output=True, timeout=600).stdout.decode()
        assert out.startswith("@@\n"), (o, n)
        print(json.dumps({"options": o, "input": b, "name": n, "output": out[3:]}))


if __name__ == "__main__":
    main()
