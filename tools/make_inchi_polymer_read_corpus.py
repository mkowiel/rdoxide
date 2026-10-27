#!/usr/bin/env python3
"""Builds tests/data/inchi_polymers_read.jsonl: the polymer InChI strings (`InChI=1B/.../z...`)
produced from the polymer molfile corpora and what the official RDKit 2026.09.1 makes of them
(`MolFromInchi`, no sanitization, as a molblock; null when it fails, with the log).

    .venv-ref/bin/python tools/make_inchi_polymer_read_corpus.py > tests/data/inchi_polymers_read.jsonl
"""
import json
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdinchi

DATA = Path(__file__).resolve().parents[1] / "tests" / "data"


def convert(s):
    RDLogger.DisableLog("rdApp.*")
    m, ret, msg, log = None, None, None, None
    r = rdinchi.InchiToMol(s, False, False)
    m, ret, msg, log = r[0], r[1], r[2], r[3]
    return {"inchi": s, "molblock": Chem.MolToMolBlock(m) if m is not None else None,
            "ret": ret, "message": msg, "log": log}


def main():
    if len(sys.argv) > 1:
        print(json.dumps(convert(sys.stdin.read())))
        return
    seen = []
    for name in ("inchi_polymers.jsonl", "inchi_polymers_zz.jsonl"):
        for line in (DATA / name).read_text().splitlines():
            for res in json.loads(line)["results"].values():
                s = res[0]
                if s.startswith("InChI=1B") and s not in seen:
                    seen.append(s)
    for s in seen:
        p = subprocess.run([sys.executable, __file__, "child"], input=s, capture_output=True, text=True)
        if p.returncode != 0:  # the reference library crashes on some input: skip
            continue
        print(p.stdout.strip())


if __name__ == "__main__":
    main()
