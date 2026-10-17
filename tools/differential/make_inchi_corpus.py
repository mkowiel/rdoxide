"""Generate tests/data/inchi_corpus.jsonl with the official RDKit (InChI 1.07.3).

usage (reference interpreter, see compare.py):
    $RDOXIDE_REF_PYTHON -I tools/differential/make_inchi_corpus.py > tests/data/inchi_corpus.jsonl

One JSON object per line:
  {"kind": "mol", "smiles", "options", "inchi", "retcode", "message", "log", "aux", "key"}
      rdinchi.MolToInchi(mol, options) + InchiToInchiKey(inchi)
  {"kind": "molblock", "name", "molblock", "options", "inchi", "retcode", "message", "log", "aux", "key"}
      rdinchi.MolBlockToInchi(molblock, options)
  {"kind": "key", "inchi", "key"}
      rdinchi.InchiToInchiKey(inchi) for hand-written (also invalid) InChIs

Every InChI of the mol/molblock records is also a key test case.
"""
import json
import sys
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdinchi

RDLogger.DisableLog("rdApp.*")
assert rdBase.rdkitVersion == "2026.09.1", rdBase.rdkitVersion
assert rdinchi.GetInchiVersion() == "1.07.3", rdinchi.GetInchiVersion()

ROOT = Path(__file__).resolve().parents[2]
UPSTREAM = ROOT / "upstream" / "rdkit"
N_NCI = 150

# Hand-picked molecules exercising the InChI layers (charges, protons, isotopes,
# stereo, metals, radicals, multi-component, tautomers) and edge cases.
EDGE_SMILES = [
    "C", "[H][H]", "[H+]", "[H-]", "[2H]C([2H])([2H])[2H]", "[13CH4]", "O", "[OH-]", "[OH3+]",
    "CC(=O)[O-]", "CC(=O)O", "C[N+](C)(C)C", "[Na+].[Cl-]", "[Na]Cl", "OC(=O)CC(O)(CC(=O)O)C(=O)O",
    "OC(=O)CC(O)(CC(=O)[O-])C(=O)[O-].[Na+].[Na+]", "[O-]C(=O)CC(=O)[O-].[O-]C(=O)CC(=O)[O-]",
    "[O-]S(=O)(=O)[O-]", "OP(=O)(O)O", "[O-]P(=O)([O-])[O-]",
    "C(C(=O)[O-])C(=O)[O-].C(C(=O)[O-])C(=O)[O-].C(C(=O)[O-])C(=O)[O-].C(C(=O)[O-])C(=O)[O-]"
    ".C(C(=O)[O-])C(=O)[O-].C(C(=O)[O-])C(=O)[O-].C(C(=O)[O-])C(=O)[O-]",
    "[NH4+].[NH4+].[NH4+].[NH4+].[NH4+].[NH4+].[NH4+].[NH4+].[NH4+].[NH4+].[NH4+].[NH4+].[NH4+]",
    "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "CC(N)C(=O)O", "C/C=C/C", "C/C=C\\C", "CC=CC",
    "C[C@H]1CC[C@@H](C)CC1", "C1CC[C@H]2CCCC[C@@H]2C1", "CC=C=CC", "C[C@H](O)[C@@H](O)C",
    "O=C1NC=CC=C1", "Oc1ccccn1", "CC(=O)CC(C)=O", "N=C(N)N", "c1ccccc1", "C1=CC=CC=C1",
    "c1ccc2ccccc2c1", "c1cc[nH]c1", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "[CH3]", "[CH2]", "C[O]",
    "[Fe]", "[Fe+2]", "Cl[Pt](Cl)([NH3])[NH3]", "[Cu+2].[O-]S(=O)(=O)[O-]", "C[Mg]Br",
    "[Li]C", "CC[Sn](CC)(CC)CC", "[U]", "[Xe]", "[He]", "[Og]", "B(F)(F)F", "F[B-](F)(F)F",
    "C1CCCCCCCC/C=C/1", "C1CCCCCCC/C=C/CCCCCCC1", "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](O)CC[C@]34C)"
    "[C@@H]1CC[C@@H]2O", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "CC(C)(C)c1ccc(O)cc1",
    "[2H]C(Cl)(Cl)Cl", "[3H]O[3H]", "[18F]c1ccccc1", "[C-]#[O+]", "[N-]=[N+]=N", "N#N",
    "C1=CC=CC=CC=C1", "c1ccc(cc1)-c1ccccc1", "OC1=CC=CC1=O", "*C", "",
]

OPTION_SETS = ["", "/FixedH /RecMet", "/SUU /SLUUD /SNon", "/KET /15T"]

KEY_CASES = [
    "", "foo", "InChI=", "InChI=1", "InChI=1S", "InChI=1S/", "InChI=2S/CH4/h1H4", "InChI=1X/CH4/h1H4",
    "InChI=1S/CH4/h1H4", "InChI=1/CH4/h1H4", "InChI=1B/CH4/h1H4", "InChI=1S//", "InChI=1//",
    "InChI=1S/?", "InChI=1S/%", "InChI=1S/CH4/h1H4 trailing text", "InChI=1S/CH4/h1H4\n",
    "InChI=1S/CH4/h1H4\tAuxInfo=1/1/N:1", "InChI=1S/CH4/h1H4/", "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)/p-1",
    "InChI=1S/H3N/h1H3/p+1", "InChI=1S/H3N/h1H3/p+12", "InChI=1S/H3N/h1H3/p+13", "InChI=1S/H3N/h1H3/p-12",
    "InChI=1S/H3N/h1H3/p-13", "InChI=1S/H3N/h1H3/p", "InChI=1S/H3N/h1H3/p+0", "InChI=1S/H3N/h1H3/px",
    "InChI=1S/H3N/h1H3/p+99999999999999999999", "InChI=1S/CH4/h1H4/f/h", "InChI=1/CH4/h1H4/f/h",
    "InChI=1S/CH4/h1H4/r", "InChI=1/CH4/h1H4/rCH4/h1H4", "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+",
    "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1",
    "InChI=1S/CH4/i1+1", "InChI=1S/ClH/h1H/q", "InChI=1S/x", "InChI=1S/c", "xInChI=1S/CH4/h1H4",
    "InChI=1S/CH4/h1H4/p+1/q+1", "InChI=1S/" + "C" * 300, "InChI=1/" + "C" * 120 + "/b" + "1-2+," * 60,
]


def mol_record(smi, options):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    inchi, retcode, message, log, aux = rdinchi.MolToInchi(mol, options)
    return {"kind": "mol", "smiles": smi, "options": options, "inchi": inchi, "retcode": retcode,
            "message": message, "log": log, "aux": aux, "key": rdinchi.InchiToInchiKey(inchi)}


def main():
    out = sys.stdout
    nci = []
    with open(UPSTREAM / "Data" / "NCI" / "first_5K.smi") as fh:
        for line in fh:
            parts = line.split()
            if parts and len(nci) < N_NCI:
                nci.append(parts[0])
    stereo = [l.split()[0] for l in open(ROOT / "tests" / "data" / "stereo_edge_cases.smi") if l.strip()]
    seen = set()
    for i, smi in enumerate(EDGE_SMILES + stereo + nci):
        # option variants for the hand-picked molecules and every 25th corpus molecule
        opts = OPTION_SETS if i < len(EDGE_SMILES) or i % 25 == 0 else [""]
        for options in opts:
            if (smi, options) in seen:
                continue
            seen.add((smi, options))
            rec = mol_record(smi, options)
            if rec is not None:
                out.write(json.dumps(rec) + "\n")
    for path in sorted((UPSTREAM / "External" / "INCHI-API" / "test_data").glob("*.mol")):
        blk = path.read_text()
        for options in ("", "/FixedH"):
            inchi, retcode, message, log, aux = rdinchi.MolBlockToInchi(blk, options)
            out.write(json.dumps({"kind": "molblock", "name": path.name, "molblock": blk, "options": options,
                                  "inchi": inchi, "retcode": retcode, "message": message, "log": log,
                                  "aux": aux, "key": rdinchi.InchiToInchiKey(inchi)}) + "\n")
    for inchi in KEY_CASES:
        out.write(json.dumps({"kind": "key", "inchi": inchi, "key": rdinchi.InchiToInchiKey(inchi)}) + "\n")


if __name__ == "__main__":
    main()
