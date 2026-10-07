#!/usr/bin/env python3
"""Differential testing of rdoxide against the official RDKit wheel.

    python tools/differential/compare.py MODE [--limit N] [FILES...]

Runs tools/differential/dump.py under the reference interpreter
(.venv-ref, official `rdkit`) and the rdoxide interpreter (.venv) and
reports molecules whose dumps differ. Default corpus: SMILES files from the
upstream checkout. Exit code 1 if anything differs.
"""
import argparse
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
    ap.add_argument("mode", choices=["parse", "sanitize", "smiles"])
    ap.add_argument("files", nargs="*")
    ap.add_argument("--limit", type=int, default=20, help="number of differences to print")
    ap.add_argument("--ref-python", default=str(ROOT / ".venv-ref/bin/python"))
    ap.add_argument("--test-python", default=str(ROOT / ".venv/bin/python"))
    args = ap.parse_args()
    files = [Path(f) for f in (args.files or [ROOT / f for f in DEFAULT_CORPUS])]
    total = ndiff = 0
    for f in files:
        ref = run(args.ref_python, args.mode, f)
        test = run(args.test_python, args.mode, f)
        for (smi, r), (smi2, t) in zip(ref, test):
            assert smi == smi2
            total += 1
            if r != t:
                ndiff += 1
                if ndiff <= args.limit:
                    print(f"DIFF {smi}\n  ref:  {json.dumps(r)[:400]}\n  test: {json.dumps(t)[:400]}")
    print(f"{args.mode}: {total - ndiff}/{total} identical")
    return 1 if ndiff else 0


if __name__ == "__main__":
    sys.exit(main())
