#!/usr/bin/env python3
"""Builds tests/data/avalon_fixcharges.jsonl: charged molecules sent through Avalon's struchk with the charge table
(`-tc`, `-tq`, `-tl`) and what the C code (RechargeMolecule) answers: result bits, messages and molecule
(driver: .work/avalon_drv/sc.c, the same as for avalon_struchk.jsonl). `{D}` stands for the directory with checkfgs.aci.

    python3 tools/make_avalon_fixcharges_corpus.py DRIVER DATADIR > tests/data/avalon_fixcharges.jsonl
"""
import json
import subprocess
import sys

OPTS = [
    "-tc {D}/checkfgs.aci|-tq 0|-tl 24",
    "-tc {D}/checkfgs.aci|-tq -1|-tl 24",
    "-tc {D}/checkfgs.aci|-tq 1|-tl 24",
    "-tc {D}/checkfgs.aci|-tq 0|-tl 8",
    "-tl 24|-tc {D}/checkfgs.aci|-tq -2",
]

SMILES = [
    "[O-]C(=O)C", "[O-]C(=O)CC(=O)[O-]", "[O-]C(=O)CCC(=O)[O-]", "[O-]C(=O)CCCC(=O)[O-]", "[O-]C(=O)c1ccccc1C(=O)[O-]",
    "[O-]C(=O)c1ccc(cc1)C(=O)[O-]", "[O-]c1ccccc1", "[O-]c1ccc(cc1)C(=O)[O-]", "[O-]c1ccc(cc1)[N+](=O)[O-]",
    "[S-]CC(=O)[O-]", "[S-]c1ccccc1", "[O-]S(=O)(=O)CC(=O)[O-]", "[O-]S(=O)(=O)c1ccccc1", "[O-]P(=O)([O-])OCC",
    "[O-]P(=O)([O-])[O-]", "[O-]P(=O)([O-])OP(=O)([O-])[O-]", "CC([O-])=O.[Na+]", "[NH3+]CC(=O)[O-]", "[NH3+]CC[O-]",
    "[NH3+]CCC(=O)[O-]", "C[NH+](C)CC(=O)[O-]", "[NH3+]C(CC(=O)[O-])C(=O)[O-]", "[NH3+]C(CCC(=O)[O-])C(=O)[O-]",
    "[NH3+]CCCC[C@H]([NH3+])C(=O)[O-]", "[NH3+][C@@H](Cc1ccc([O-])cc1)C(=O)[O-]", "[NH3+][C@@H](CS)C(=O)[O-]",
    "[NH3+][C@@H](C[S-])C(=O)[O-]", "[O-]C(=O)[C@H](C)C[C@@H](C)C(=O)[O-]", "C[C@H](C(=O)[O-])[C@@H](C)C(=O)[O-]",
    "[O-]C(=O)[C@H](O)[C@@H](O)C(=O)[O-]", "[O-]C(=O)C(O)C(O)C(=O)[O-]", "[O-]C(=O)CC(O)(CC(=O)[O-])C(=O)[O-]",
    "[O-]C(=O)CC(CC(=O)[O-])C(=O)[O-]", "[N-]S(=O)(=O)c1ccccc1", "C[N-]S(=O)(=O)C", "c1nnn[n-]1", "c1ccc2[n-]nnc2c1",
    "[n-]1cccc1", "c1cc[n-]c1C(=O)[O-]", "[O-]C(=O)c1ccncc1", "[O-]C(=O)c1cccnc1.[Cl-]", "CC(=O)[CH-]C(C)=O",
    "N#C[CH-]C#N", "[O-]N=O", "[O-][N+](=O)[O-]", "[O-]C(=O)O", "[O-]C(=O)[O-]", "OC(=O)C(=O)[O-]", "[O-]C(=O)C(=O)[O-]",
    "[O-]C(=O)C=CC(=O)[O-]", "[O-]C(=O)/C=C/C(=O)[O-]", "[O-]C(=O)/C=C\\C(=O)[O-]", "[O-]C(=O)C#CC(=O)[O-]",
    "C[N+](C)(C)CC(=O)[O-]", "C[N+](C)(C)CC[O-]", "[O-]C(=O)C[N+](C)(C)C.[O-]C(=O)C", "CC(=O)[O-].CC(=O)[O-]",
    "[O-]C(=O)C.[O-]C(=O)CC", "[O-]c1ccccc1[O-]", "[O-]c1cccc([O-])c1", "[O-]c1ccc([O-])cc1", "[O-]C(=O)c1ccccc1[O-]",
    "[O-]C(=O)CCCCCCCC(=O)[O-]", "[O-]CC[O-]", "[O-]CCC[O-]", "OCC[O-]", "[O-]C(F)(F)F", "CC[S-]", "CS(=O)(=O)[N-]C(=O)C",
    "O=C([O-])c1ccc2ccccc2c1", "[O-]C(=O)C1CCCCC1C(=O)[O-]", "[O-]C(=O)[C@H]1CCCC[C@@H]1C(=O)[O-]",
    "[O-]C(=O)[C@H]1CCCC[C@H]1C(=O)[O-]", "[O-]C(=O)C(Cl)Cl", "[O-]C(=O)C(F)(F)F", "[O-]C(=O)c1ccccc1.[O-]c1ccccc1",
    "[Fe+2].[O-]C(=O)C.[O-]C(=O)C", "[Cu+2].[O-]C(=O)C.[O-]C(=O)C", "[Na+].[Na+].[O-]C(=O)CC(=O)[O-]",
    "N=C(N)NCCC[C@H]([NH3+])C(=O)[O-]", "[NH3+]CC(=O)NCC(=O)[O-]", "[NH3+]CC(=O)N[C@@H](CC(=O)[O-])C(=O)[O-]",
    "C(=O)([O-])C[N+]1=CC=CC=C1", "[O-]P(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O",
    "[O-]P(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O",
    "[O-]C(=O)C[C@H](N)C(=O)[O-]", "[O-]C(=O)CS", "[O-]C(=O)CO", "[O-]C(=O)CN", "[O-]C(=O)C[N+](=O)[O-]",
    "[O-]S(=O)(=O)O", "[O-]S(=O)(=O)[O-]", "O=S([O-])[O-]", "[O-][Cl+3]([O-])([O-])[O-]", "[O-]B(O)O",
    "[CH2-]C(=O)C", "C[CH-]C#N", "[C-]#N.[Na+]", "O=C1[N-]C(=O)c2ccccc12", "O=C1[N-]C(=O)NC(=O)C1",
    "[O-]C(=O)[C@H](N)[C@@H](N)C(=O)[O-]", "[O-]C(=O)[C@H](Cl)[C@H](Cl)C(=O)[O-]", "[O-]C(=O)[C@@H](C)[C@@H](C)C(=O)[O-]",
    "[O-]C(=O)[C@H](C)[C@H](C(=O)[O-])[C@@H](C)C(=O)[O-]", "C[C@H]([O-])[C@@H](C)[O-]", "C[C@H](C(=O)[O-])N(C)[C@H](C)C(=O)[O-]",
    "[O-]C(=O)C(C)(C)C(=O)[O-]", "[O-]C(=O)C1(CC1)C(=O)[O-]", "[O-]C(=O)C=C(C(=O)[O-])C(=O)[O-]",
    "O=C1CC(=O)[N-]C(=O)N1", "CC1=CC(=O)[N-]C(=O)N1", "c1ccc2c(c1)[nH]c1ccccc12.[O-]C(=O)C", "CC(C)C(=O)[O-].CC(C)C(=O)[O-]",
]


def main():
    driver, datadir = sys.argv[1], sys.argv[2]
    recs = []
    for k, s in enumerate(SMILES):
        for j in (k % len(OPTS), (k + 2) % len(OPTS)):
            recs.append((OPTS[j], "SMILES:" + s, s))
    inp = "$$$$\n".join(o.replace("{D}", datadir) + "\n" + b + "\n" for o, b, _ in recs)
    out = subprocess.run([driver], input=inp.encode(), capture_output=True, timeout=1800).stdout.decode()
    parts = out.split("@@\n")[1:]
    assert len(parts) == len(recs), (len(parts), len(recs))
    for (o, b, n), p in zip(recs, parts):
        print(json.dumps({"options": o, "input": b, "name": n, "output": p}))


if __name__ == "__main__":
    main()
