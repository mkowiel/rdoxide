"""Dump molecule features as JSON lines. Runs unchanged under the official
RDKit and under rdoxide; tools/differential/compare.py diffs the outputs.

usage: python -I dump.py MODE INPUT.smi > out.jsonl
  MODE: parse (sanitize=False) | rings | sanitize (SanitizeMol) | molfromsmiles | smiles (canonical SMILES)
"""
import json
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def atom_record(a):
    rec = [a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsotope(), a.GetNumExplicitHs(),
           a.GetNoImplicit(), a.GetIsAromatic(), int(a.GetChiralTag()), a.GetNumRadicalElectrons(),
           sorted(a.GetPropNames())]
    return rec


STEREO_PROPS = ("_CIPCode", "_ChiralityPossible", "_CIPRank", "_ringStereoAtoms", "_UnknownStereo",
                "_StereochemDone")


def private_props(obj):
    d = obj.GetPropsAsDict(includePrivate=True, includeComputed=True)
    return {k: (list(v) if not isinstance(v, (int, float, str)) else v)
            for k, v in sorted(d.items()) if k in STEREO_PROPS}


def bond_record(b):
    return [b.GetBeginAtomIdx(), b.GetEndAtomIdx(), int(b.GetBondType()), int(b.GetBondDir()),
            int(b.GetStereo()), list(b.GetStereoAtoms()), b.GetIsAromatic()]


def dump(mode, smi):
    if mode == "parse":
        p = Chem.SmilesParserParams()
        p.sanitize = False
        p.removeHs = False
        m = Chem.MolFromSmiles(smi, p)
        if m is None:
            return None
        return {"atoms": [atom_record(a) for a in m.GetAtoms()],
                "bonds": [bond_record(b) for b in m.GetBonds()]}
    if mode == "rings":
        p = Chem.SmilesParserParams()
        p.sanitize = False
        p.removeHs = False
        m = Chem.MolFromSmiles(smi, p)
        if m is None:
            return None
        sssr = [list(r) for r in Chem.GetSSSR(m)]
        symm = [list(r) for r in Chem.GetSymmSSSR(Chem.Mol(m))]
        alg = getattr(Chem, "SymmetrizeSSSRAlgorithm", None)
        legacy = None
        if alg is not None:
            legacy = [list(r) for r in Chem.GetSymmSSSR(Chem.Mol(m), algorithm=alg.LEGACY)]
        m2 = Chem.Mol(m)
        Chem.FastFindRings(m2)
        fast = [list(r) for r in m2.GetRingInfo().AtomRings()]
        return {"sssr": sssr, "symm": symm, "legacy": legacy, "fast": fast,
                "bonds": [list(r) for r in m.GetRingInfo().BondRings()]}
    if mode == "sanitize":
        p = Chem.SmilesParserParams()
        p.sanitize = False
        p.removeHs = False
        m = Chem.MolFromSmiles(smi, p)
        if m is None:
            return None
        try:
            Chem.SanitizeMol(m)
        except Exception as e:  # noqa: BLE001
            return {"error": type(e).__name__, "msg": str(e)}
        recs = []
        for a in m.GetAtoms():
            recs.append(atom_record(a) + [a.GetTotalNumHs(), int(a.GetHybridization()),
                                          a.GetTotalValence()])
        return {"atoms": recs,
                "bonds": [bond_record(b) + [b.GetIsConjugated(), private_props(b)] for b in m.GetBonds()],
                "rings": [list(r) for r in m.GetRingInfo().AtomRings()],
                "mol": private_props(m)}
    if mode == "ranks":
        p = Chem.SmilesParserParams()
        p.sanitize = False
        p.removeHs = False
        m = Chem.MolFromSmiles(smi, p)
        if m is None:
            return None
        try:
            Chem.SanitizeMol(m)
        except Exception as e:  # noqa: BLE001
            return {"error": type(e).__name__}
        n = m.GetNumAtoms()
        res = {"ties": list(Chem.CanonicalRankAtoms(m, breakTies=False)),
               "ranks": list(Chem.CanonicalRankAtoms(m)),
               "nochiral": list(Chem.CanonicalRankAtoms(m, includeChirality=False))}
        if n > 2:
            half = list(range(n // 2))
            bonds = [b.GetIdx() for b in m.GetBonds()
                     if b.GetBeginAtomIdx() in half and b.GetEndAtomIdx() in half]
            res["frag"] = list(Chem.CanonicalRankAtomsInFragment(m, half, bonds))
        m2 = Chem.Mol(m)
        Chem.Kekulize(m2)
        res["kekule"] = [int(b.GetBondType()) for b in m2.GetBonds()]
        return res
    if mode == "molfromsmiles":
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        recs = []
        for a in m.GetAtoms():
            recs.append(atom_record(a) + [a.GetTotalNumHs(), int(a.GetHybridization()),
                                          a.IsInRing(), private_props(a)])
        return {"atoms": recs,
                "bonds": [bond_record(b) + [b.GetIsConjugated(), private_props(b)] for b in m.GetBonds()],
                "rings": [list(r) for r in m.GetRingInfo().AtomRings()],
                "mol": private_props(m)}
    if mode == "smiles":
        m = Chem.MolFromSmiles(smi)
        return None if m is None else Chem.MolToSmiles(m)
    raise SystemExit(f"unknown mode {mode}")


def main():
    mode, path = sys.argv[1], sys.argv[2]
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if not parts or parts[0].lower() == "smiles":
                continue
            smi = parts[0]
            try:
                res = dump(mode, smi)
            except Exception as e:  # noqa: BLE001
                res = {"exception": type(e).__name__}
            print(json.dumps([smi, res]))


if __name__ == "__main__":
    main()
