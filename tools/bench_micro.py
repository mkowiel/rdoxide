#!/usr/bin/env python3
"""Microbenchmarks of common cheminformatics operations (ms per 1000 molecules).

    python tools/bench_micro.py rdoxide
    python tools/bench_micro.py rdkit      # original RDKit, other venv
Prints one JSON object {operation: milliseconds per pass over the molecule set}.
"""

import json
import sys
import time
from pathlib import Path

lib = sys.argv[1] if len(sys.argv) > 1 else "rdoxide"
if lib == "rdkit":
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem, Descriptors, rdFingerprintGenerator, rdMolDescriptors
else:
    from rdoxide import Chem, DataStructs, RDLogger
    from rdoxide.Chem import AllChem, Descriptors, rdFingerprintGenerator, rdMolDescriptors
RDLogger.DisableLog("rdApp.*")

DATA = Path(__file__).resolve().parent.parent / "upstream" / "rdkit" / "Docs" / "Book" / "data"
mols = []
for name in ("solubility.train.sdf", "bzr.sdf", "actives_5ht3.sdf"):
    mols += [m for m in Chem.SDMolSupplier(str(DATA / name)) if m is not None]
mols = mols[:1000]
smiles = [Chem.MolToSmiles(m) for m in mols]
molblocks = [Chem.MolToMolBlock(m) for m in mols]
unsanitized = [Chem.MolFromSmiles(s, sanitize=False) for s in smiles]
noh = [Chem.RemoveHs(m) for m in mols]
embed_mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in smiles[:30]]


def each(f, seq=None):
    seq = mols if seq is None else seq

    def run():
        for x in seq:
            f(x)

    return run


morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
rdk = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)
ap = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)
tt = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)
fps = [morgan.GetFingerprint(m) for m in mols]


def mol_from_smiles_each():
    for s in smiles:
        Chem.MolFromSmiles(s)


def bulk_tanimoto():
    for f in fps[:100]:
        DataStructs.BulkTanimotoSimilarity(f, fps)


def embed():
    for m in embed_mols:
        AllChem.EmbedMolecule(Chem.Mol(m), randomSeed=42)


def mmff():
    for m in embed_mols[:10]:
        m = Chem.Mol(m)
        AllChem.EmbedMolecule(m, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(m)


BENCH = {
    "MolFromSmiles": mol_from_smiles_each,
    "MolFromSmiles sanitize=False": lambda: [Chem.MolFromSmiles(s, sanitize=False) for s in smiles],
    "MolToSmiles": each(Chem.MolToSmiles),
    "MolToSmiles isomeric=False": each(lambda m: Chem.MolToSmiles(m, isomericSmiles=False)),
    "MolFromMolBlock": lambda: [Chem.MolFromMolBlock(b) for b in molblocks],
    "MolToMolBlock": each(Chem.MolToMolBlock),
    "SanitizeMol": lambda: [Chem.SanitizeMol(Chem.Mol(m)) for m in unsanitized],
    "Mol copy": each(Chem.Mol),
    "AddHs": each(Chem.AddHs),
    "RemoveHs": each(Chem.RemoveHs, [Chem.AddHs(m) for m in mols[:300]]),
    "Kekulize": each(lambda m: Chem.Kekulize(Chem.Mol(m), clearAromaticFlags=True)),
    "AssignStereochemistry": each(lambda m: Chem.AssignStereochemistry(Chem.Mol(m), cleanIt=True, force=True)),
    "GetSymmSSSR": each(Chem.GetSymmSSSR),
    "FindMolChiralCenters": each(lambda m: Chem.FindMolChiralCenters(m, includeUnassigned=True)),
    "Morgan r2 2048": each(morgan.GetFingerprint),
    "Morgan count": each(morgan.GetCountFingerprint),
    "RDKit FP": each(rdk.GetFingerprint),
    "AtomPair": each(ap.GetFingerprint),
    "TopologicalTorsion": each(tt.GetFingerprint),
    "MACCS": each(lambda m: rdMolDescriptors.GetMACCSKeysFingerprint(m)),
    "BulkTanimoto 100x1000": bulk_tanimoto,
    "MolWt": each(Descriptors.MolWt),
    "MolLogP (Crippen)": each(Descriptors.MolLogP),
    "TPSA": each(rdMolDescriptors.CalcTPSA),
    "all Descriptors": each(lambda m: Descriptors.CalcMolDescriptors(m), mols[:100]),
    "Gasteiger": each(lambda m: AllChem.ComputeGasteigerCharges(Chem.Mol(m))),
    "GetDistanceMatrix": each(Chem.GetDistanceMatrix),
    "Embed 30 mols": embed,
    "Embed+MMFF 10 mols": mmff,
}
out = {}
for label, fn in BENCH.items():
    # embedding is slow: one run; everything else is best of three
    repeats = 1 if fn is embed or fn is mmff else 3
    best = float("inf")
    for _ in range(repeats):
        t = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t)
    out[label] = best * 1e3
print(json.dumps(out, indent=1))
