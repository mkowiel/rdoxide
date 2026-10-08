#!/usr/bin/env python3
"""Differential testing of rdoxide against the official RDKit wheel.

    python tools/differential/compare.py MODE [--limit N] [FILES...]

Runs tools/differential/dump.py under the reference interpreter
(official `rdkit` of the ported version, e.g. a conda-forge env:
`micromamba create -p REF -c conda-forge python=3.13 rdkit=2026.09.1`,
selected with --ref-python or RDOXIDE_REF_PYTHON) and the rdoxide interpreter (.venv) and
reports molecules whose dumps differ. Default corpus: SMILES files from the
upstream checkout. Exit code 1 if anything differs.
"""
import argparse
import os
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DUMP = Path(__file__).resolve().parent / "dump.py"
DEFAULT_CORPUS = [
    "upstream/rdkit/Data/NCI/first_5K.smi",
    "upstream/rdkit/Code/GraphMol/test_data/canonSmiles.long.smi",
]


def run(python, mode, path):
    out = subprocess.run([python, "-I", str(DUMP), mode, str(path)], capture_output=True, text=True, check=True)
    return [json.loads(line) for line in out.stdout.splitlines()]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["parse", "rings", "sanitize", "ranks", "molfromsmiles", "smiles", "smarts", "substruct", "pickle", "picklesmarts", "molblock", "sdf", "molwrite", "depict", "sdwrite", "morgan", "atompair", "subgraphs", "torsion", "rdkitfp", "rdkfingerprint", "moldescfp", "maccs", "crippen", "addhs", "addhs3d", "lipinski"])
    ap.add_argument("files", nargs="*")
    ap.add_argument("--limit", type=int, default=20, help="number of differences to print")
    ap.add_argument("--ref-python", default=os.environ.get("RDOXIDE_REF_PYTHON", str(ROOT / ".venv-ref/bin/python")),
                    help="interpreter with the official RDKit of the ported version (env RDOXIDE_REF_PYTHON)")
    ap.add_argument("--test-python", default=str(ROOT / ".venv/bin/python"))
    args = ap.parse_args()
    files = [Path(f) for f in (args.files or [ROOT / f for f in DEFAULT_CORPUS])]
    total = ndiff = 0
    for f in files:
        ref = run(args.ref_python, args.mode, f)
        test = run(args.test_python, args.mode, f)
        if len(ref) != len(test):
            total += 1
            ndiff += 1
            print(f"DIFF {f}: {len(ref)} records vs {len(test)}")
            continue
        for (smi, r), (smi2, t) in zip(ref, test):
            total += 1
            if smi != smi2:
                r, t = [smi, r], [smi2, t]
            if r != t:
                ndiff += 1
                if ndiff <= args.limit:
                    print(f"DIFF {f.name if args.mode in ('sdf', 'sdwrite') else ''} {smi}\n  ref:  {json.dumps(r)[:400]}\n  test: {json.dumps(t)[:400]}")
    print(f"{args.mode}: {total - ndiff}/{total} identical")
    return 1 if ndiff else 0


if __name__ == "__main__":
    sys.exit(main())
