import json
import os
import pathlib

import pytest

from rdoxide import Chem, RDConfig
from rdoxide.Avalon import pyAvalonTools as p

DATA = pathlib.Path(__file__).parent.parent / "data" / "avalon_checkmol.jsonl"
CASES = [json.loads(line) for line in DATA.read_text().splitlines()]
CONF = os.path.join(RDConfig.RDDataDir, "struchk", "")
INIT = f"""-tm
-ta {CONF}checkfgs.trn
-tm
-or
-ca {CONF}checkfgs.chk
-cc
-cl 3
-cs
-cn 999"""


@pytest.fixture(scope="module", autouse=True)
def _init():
    assert p.InitializeCheckMol(INIT) == 0
    yield
    p.CloseCheckMolFiles()


def test_all_cases_match_reference():
    bad = []
    for rec in CASES:
        if "exception" in rec:
            with pytest.raises(Exception):
                p.CheckMolecule(rec["input"], rec["isSmiles"])
            continue
        err, mol = p.CheckMolecule(rec["input"], rec["isSmiles"])
        got = {"errs": int(err), "molblock": Chem.MolToMolBlock(mol) if mol is not None else None,
               "log": p.GetCheckMolLog()}
        e2, s = p.CheckMoleculeString(rec["input"], rec["isSmiles"])
        got.update(str_errs=int(e2), str=s)
        want = {k: rec[k] for k in ("errs", "molblock", "log", "str_errs", "str")}
        if got != want:
            bad.append((rec["name"], [k for k in want if want[k] != got[k]]))
    assert not bad, f"{len(bad)} of {len(CASES)} differ, e.g. {bad[:4]}"


def test_enums_and_helpers():
    assert p.StruChkResult.success == 0
    assert p.StruChkFlag.bad_molecule | p.StruChkFlag.stereo_error == 33
    assert p.StruChkResult.bad_set & p.StruChkFlag.atom_clash
    err, mol = p.CheckMolecule(Chem.MolFromSmiles("CCO"))
    assert err == 0 and mol.GetNumAtoms() == 3
    assert p.CheckMolecule("not a smiles", True) == (1, None)
