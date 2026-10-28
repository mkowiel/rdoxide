#!/usr/bin/env python3
"""Builds tests/data/avalon_checkmol.jsonl: what the official rdkit.Avalon.pyAvalonTools answers for CheckMolecule,
CheckMoleculeString and GetCheckMolLog with the in-memory MolKey style option set.

    .venv-ref/bin/python tools/make_avalon_checkmol_corpus.py > tests/data/avalon_checkmol.jsonl
"""
import json
import os
import sys

from rdkit import Chem, RDConfig, RDLogger
from rdkit.Avalon import pyAvalonTools as p
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
sys.path.insert(0, os.path.dirname(__file__))
from make_avalon_struchk_corpus import inputs  # noqa: E402

DATA = os.path.join(RDConfig.RDDataDir, "struchk", "")
INIT = f"""-tm
-ta {DATA}checkfgs.trn
-tm
-or
-ca {DATA}checkfgs.chk
-cc
-cl 3
-cs
-cn 999"""


def main():
    assert p.InitializeCheckMol(INIT) == 0
    for body, name in inputs():
        is_smiles = body.startswith("SMILES:")
        data = body[7:].strip() if is_smiles else body
        rec = {"name": name, "input": data, "isSmiles": is_smiles}
        try:
            err, mol = p.CheckMolecule(data, is_smiles)
            rec["errs"] = int(err)
            rec["molblock"] = Chem.MolToMolBlock(mol) if mol is not None else None
        except Exception as e:  # noqa: BLE001
            rec["exception"] = type(e).__name__
        rec["log"] = p.GetCheckMolLog()
        try:
            e2, s = p.CheckMoleculeString(data, is_smiles)
            rec["str_errs"] = int(e2)
            rec["str"] = s
        except Exception as e:  # noqa: BLE001
            rec["str_exception"] = type(e).__name__
        print(json.dumps(rec))
    p.CloseCheckMolFiles()


if __name__ == "__main__":
    main()
