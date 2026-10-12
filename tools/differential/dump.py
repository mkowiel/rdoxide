"""Dump molecule features as JSON lines. Runs unchanged under the official
RDKit and under rdoxide; tools/differential/compare.py diffs the outputs.

usage: python -I dump.py MODE INPUT.smi > out.jsonl
  MODE: parse (sanitize=False) | rings | sanitize (SanitizeMol) | molfromsmiles | smiles (SMILES writer) | smarts (MolFromSmarts)
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


def molblock_record(blk, **kw):
    m = Chem.MolFromMolBlock(blk, **kw)
    if m is None:
        return None
    rec = {"smiles": Chem.MolToSmiles(m),
           "atoms": [atom_record(a) + [private_props(a)] for a in m.GetAtoms()],
           "bonds": [bond_record(b) + [private_props(b)] for b in m.GetBonds()],
           "mol": private_props(m)}
    try:
        rec["pickle"] = m.ToBinary(Chem.PropertyPickleOptions.AllProps).hex()
    except Exception as e:  # noqa: BLE001
        rec["pickle"] = type(e).__name__
    return rec


def sdwrite_records(path):
    """SDWriter output for the molecules of an SD file."""
    import io
    try:
        mols = [m for m in Chem.SDMolSupplier(path) if m is not None][:300]
    except Exception as e:  # noqa: BLE001
        yield ["open", type(e).__name__]
        return
    for key, setup in (("default", lambda w: None), ("v3000", lambda w: w.SetForceV3000(True)),
                       ("nokekule", lambda w: w.SetKekulize(False)),
                       ("props", lambda w: w.SetProps(["NAME", "ID", "AMW", "missing"]))):
        sio = io.StringIO()
        w = Chem.SDWriter(sio)
        setup(w)
        out = []
        for m in mols:
            try:
                w.write(m)
            except Exception as e:  # noqa: BLE001
                out.append(type(e).__name__)
        w.flush()
        yield [key, [sio.getvalue(), out, w.NumMols()]]
    def text(m):
        try:
            return Chem.SDWriter.GetText(m)
        except Exception as e:  # noqa: BLE001
            return type(e).__name__
    yield ["gettext", [text(m) for m in mols[:20]]]


def sdf_records(path):
    """Everything an SDMolSupplier reports for one SD file."""
    def mrec(m):
        if m is None:
            return None
        props = m.GetPropsAsDict(includePrivate=False, includeComputed=False)
        return {"smiles": Chem.MolToSmiles(m), "props": sorted((k, str(v)) for k, v in props.items()),
                "atomprops": [sorted((k, str(v)) for k, v in a.GetPropsAsDict(False, False).items())
                              for a in m.GetAtoms()],
                "bondprops": [sorted((k, str(v)) for k, v in b.GetPropsAsDict(False, False).items())
                              for b in m.GetBonds()]}

    def guard(f):
        try:
            return f()
        except Exception as e:  # noqa: BLE001
            return {"exception": type(e).__name__}

    for key, kw in (("default", {}), ("keepHs", {"removeHs": False}),
                    ("lax", {"strictParsing": False}), ("nosanitize", {"sanitize": False})):
        try:
            suppl = Chem.SDMolSupplier(path, **kw)
        except Exception as e:  # noqa: BLE001
            yield [f"{key} open", type(e).__name__]
            continue
        mols = []
        while True:
            try:
                m = next(suppl)
            except StopIteration:
                break
            except Exception as e:  # noqa: BLE001
                mols.append({"exception": type(e).__name__})
                break
            mols.append(mrec(m))
            if len(mols) > 2000:
                break
        yield [f"{key} iterate", mols]
        if key != "default":
            continue
        suppl = Chem.SDMolSupplier(path, **kw)
        n = guard(lambda: len(suppl))
        yield ["len", n]
        idxs = [0, 1, 3, -1, 7, 2]
        yield ["getitem", [guard(lambda i=i: mrec(suppl[i])) for i in idxs]]
        yield ["text", [guard(lambda i=i: suppl.GetItemText(i)) for i in (0, 2, 1)]]
        yield ["after", [guard(lambda: mrec(next(suppl))), guard(lambda: suppl.atEnd())]]
        suppl2 = Chem.SDMolSupplier()
        suppl2.SetData(open(path, encoding="utf-8", errors="replace").read(), **kw)
        yield ["setdata", [mrec(m) for m in suppl2]]


def molwrite_record(blk):
    res = {}
    for key, kw in (("default", {}), ("keepHs", {"removeHs": False}), ("nosanitize", {"sanitize": False})):
        try:
            m = Chem.MolFromMolBlock(blk, **kw)
        except Exception as e:  # noqa: BLE001
            res[key] = type(e).__name__
            continue
        if m is None:
            res[key] = None
            continue
        out = {}
        for wkey, f in (("block", lambda: Chem.MolToMolBlock(m)),
                        ("nokekule", lambda: Chem.MolToMolBlock(m, kekulize=False)),
                        ("v3k", lambda: Chem.MolToV3KMolBlock(m)),
                        ("v2k", lambda: Chem.MolToV2KMolBlock(m)),
                        ("nostereo", lambda: Chem.MolToMolBlock(m, includeStereo=False)),
                        ("roundtrip", lambda: Chem.MolToSmiles(Chem.MolFromMolBlock(Chem.MolToMolBlock(m))))):
            try:
                out[wkey] = f()
            except Exception as e:  # noqa: BLE001
                out[wkey] = type(e).__name__
        res[key] = out
    return res


MORGAN_CONFIGS = (
    {}, {"radius": 2, "fpSize": 1024}, {"includeChirality": True}, {"countSimulation": True},
    {"useBondTypes": False, "onlyNonzeroInvariants": True},
    {"includeRedundantEnvironments": True, "radius": 1},
)


def morgan_record(smi):
    """Morgan generator output (all fingerprint kinds + AdditionalOutput)."""
    from rdkit.Chem import rdFingerprintGenerator as G
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    res = {}
    for i, kw in enumerate(MORGAN_CONFIGS):
        g = G.GetMorganGenerator(**kw)
        res[f"bits{i}"] = list(g.GetFingerprint(m).GetOnBits())
    g = G.GetMorganGenerator()
    res["sparsecount"] = sorted(g.GetSparseCountFingerprint(m).GetNonzeroElements().items())
    res["count"] = sorted(g.GetCountFingerprint(m).GetNonzeroElements().items())
    res["sparse"] = list(G.GetMorganGenerator(countSimulation=True).GetSparseFingerprint(m).GetOnBits())
    res["fromatoms"] = list(g.GetFingerprint(m, fromAtoms=[0]).GetOnBits())
    ao = G.AdditionalOutput()
    ao.AllocateBitInfoMap()
    ao.AllocateAtomCounts()
    ao.AllocateAtomToBits()
    ao.AllocateAtomsPerBit()
    G.GetMorganGenerator(radius=2, fpSize=1024).GetFingerprint(m, additionalOutput=ao)
    res["ao"] = [sorted(ao.GetBitInfoMap().items()), ao.GetAtomCounts(), ao.GetAtomToBits(),
                 sorted(ao.GetAtomsPerBit().items())]
    return res


ATOMPAIR_CONFIGS = (
    {}, {"countSimulation": False}, {"includeChirality": True}, {"minDistance": 2, "maxDistance": 5, "fpSize": 512},
    {"atomInvariantsGenerator": "chiral"},
)


def atompair_record(smi):
    """Atom-pair generator output (all fingerprint kinds + AdditionalOutput)."""
    from rdkit.Chem import rdFingerprintGenerator as G
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    res = {}
    for i, kw in enumerate(ATOMPAIR_CONFIGS):
        kw = dict(kw)
        if kw.get("atomInvariantsGenerator") == "chiral":
            kw["atomInvariantsGenerator"] = G.GetAtomPairAtomInvGen(includeChirality=True)
        g = G.GetAtomPairGenerator(**kw)
        res[f"bits{i}"] = list(g.GetFingerprint(m).GetOnBits())
        res[f"sparsecount{i}"] = sorted(g.GetSparseCountFingerprint(m).GetNonzeroElements().items())
    g = G.GetAtomPairGenerator()
    res["count"] = sorted(g.GetCountFingerprint(m).GetNonzeroElements().items())
    res["sparse"] = list(g.GetSparseFingerprint(m).GetOnBits())
    res["fromatoms"] = list(g.GetFingerprint(m, fromAtoms=[0], ignoreAtoms=[1]).GetOnBits())
    ao = G.AdditionalOutput()
    ao.AllocateBitInfoMap()
    ao.AllocateAtomCounts()
    ao.AllocateAtomToBits()
    ao.AllocateAtomsPerBit()
    G.GetAtomPairGenerator(countSimulation=False).GetSparseCountFingerprint(m, additionalOutput=ao)
    res["ao"] = [sorted(ao.GetBitInfoMap().items()), ao.GetAtomCounts(), ao.GetAtomToBits(),
                 sorted(ao.GetAtomsPerBit().items())]
    return res


def subgraphs_record(smi):
    """Subgraph / path enumeration results."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    lst = lambda r: [list(p) for p in r]  # noqa: E731
    res = {}
    for n in (1, 2, 3, 4):
        res[f"sub{n}"] = lst(Chem.FindAllSubgraphsOfLengthN(m, n))
        res[f"path{n}"] = lst(Chem.FindAllPathsOfLengthN(m, n))
        res[f"apath{n}"] = lst(Chem.FindAllPathsOfLengthN(m, n + 1, useBonds=False))
        res[f"spath{n}"] = lst(Chem.FindAllPathsOfLengthN(m, n, onlyShortestPaths=True))
    res["uniq3"] = lst(Chem.FindUniqueSubgraphsOfLengthN(m, 3))
    res["uniq3nobo"] = lst(Chem.FindUniqueSubgraphsOfLengthN(m, 3, useBO=False))
    res["mton"] = [lst(x) for x in Chem.FindAllSubgraphsOfLengthMToN(m, 2, 3)]
    res["rooted"] = lst(Chem.FindAllSubgraphsOfLengthN(m, 3, rootedAtAtom=0))
    res["rootedpath"] = lst(Chem.FindAllPathsOfLengthN(m, 3, rootedAtAtom=0))
    env = []
    for r in (1, 2, 3):
        d = {}
        env.append([list(Chem.FindAtomEnvironmentOfRadiusN(m, r, 0, atomMap=d)), sorted(d.items()),
                    list(Chem.FindAtomEnvironmentOfRadiusN(m, r + 5, 0, enforceSize=False))])
    res["env"] = env
    return res


def torsion_record(smi):
    """Topological-torsion generator output (all fingerprint kinds + AdditionalOutput)."""
    from rdkit.Chem import rdFingerprintGenerator as G
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    res = {}
    configs = ({}, {"countSimulation": False}, {"includeChirality": True}, {"torsionAtomCount": 3, "fpSize": 512},
               {"torsionAtomCount": 5}, {"atomInvariantsGenerator": "ap"})
    for i, kw in enumerate(configs):
        kw = dict(kw)
        if kw.get("atomInvariantsGenerator") == "ap":
            kw["atomInvariantsGenerator"] = G.GetAtomPairAtomInvGen()
        g = G.GetTopologicalTorsionGenerator(**kw)
        res[f"bits{i}"] = list(g.GetFingerprint(m).GetOnBits())
        res[f"sparsecount{i}"] = sorted(g.GetSparseCountFingerprint(m).GetNonzeroElements().items())
    g = G.GetTopologicalTorsionGenerator()
    res["count"] = sorted(g.GetCountFingerprint(m).GetNonzeroElements().items())
    res["fromatoms"] = list(g.GetFingerprint(m, fromAtoms=[0], ignoreAtoms=[2]).GetOnBits())
    g.GetOptions().onlyShortestPaths = True
    res["shortest"] = sorted(g.GetSparseCountFingerprint(m).GetNonzeroElements().items())
    ao = G.AdditionalOutput()
    ao.AllocateBitPaths()
    ao.AllocateAtomCounts()
    ao.AllocateAtomToBits()
    ao.AllocateAtomsPerBit()
    G.GetTopologicalTorsionGenerator(countSimulation=False).GetSparseCountFingerprint(m, additionalOutput=ao)
    res["ao"] = [sorted(ao.GetBitPaths().items()), ao.GetAtomCounts(), ao.GetAtomToBits(),
                 sorted(ao.GetAtomsPerBit().items())]
    return res


def rdkitfp_record(smi):
    """RDKit FP generator output (all fingerprint kinds + AdditionalOutput)."""
    from rdkit.Chem import rdFingerprintGenerator as G
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    res = {}
    configs = ({}, {"numBitsPerFeature": 1}, {"numBitsPerFeature": 4, "fpSize": 1024}, {"branchedPaths": False},
               {"useBondOrder": False, "maxPath": 5}, {"minPath": 2, "maxPath": 4, "countSimulation": True})
    for i, kw in enumerate(configs):
        g = G.GetRDKitFPGenerator(**kw)
        res[f"bits{i}"] = list(g.GetFingerprint(m).GetOnBits())
    g = G.GetRDKitFPGenerator(maxPath=5)
    res["sparsecount"] = sorted(g.GetSparseCountFingerprint(m).GetNonzeroElements().items())
    res["count"] = sorted(g.GetCountFingerprint(m).GetNonzeroElements().items())
    res["fromatoms"] = list(g.GetFingerprint(m, fromAtoms=[0, 1], ignoreAtoms=[2]).GetOnBits())
    ao = G.AdditionalOutput()
    ao.AllocateBitPaths()
    ao.AllocateAtomCounts()
    ao.AllocateAtomToBits()
    ao.AllocateAtomsPerBit()
    G.GetRDKitFPGenerator(maxPath=4, fpSize=512).GetFingerprint(m, additionalOutput=ao)
    res["ao"] = [sorted(ao.GetBitPaths().items()), ao.GetAtomCounts(), ao.GetAtomToBits(),
                 sorted(ao.GetAtomsPerBit().items())]
    return res


def patternfp_record(smi):
    """PatternFingerprint of the SMILES-parsed and the SMARTS-parsed molecule."""
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    rec = {}
    for kind, mol in (("smiles", Chem.MolFromSmiles(smi)), ("smarts", Chem.MolFromSmarts(smi))):
        if mol is None:
            continue
        rec[kind] = {
            f"{size}{'t' if taut else ''}": list(Chem.PatternFingerprint(mol, size, tautomerFingerprints=taut).GetOnBits())
            for size, taut in ((2048, False), (2048, True), (512, False))
        }
    return rec or None


def cx_record(m):
    """CXSMILES of a molecule with several field selections."""
    if m is None:
        return None
    from rdkit.Chem import rdmolfiles
    F = rdmolfiles.CXSmilesFields
    res = {}
    for key, flags in (("all", F.CX_ALL), ("nocoords", F.CX_ALL_BUT_COORDS), ("none", F.CX_NONE),
                       ("stereo", F.CX_ENHANCEDSTEREO | F.CX_BOND_CFG), ("props", F.CX_ATOM_PROPS | F.CX_ATOM_LABELS)):
        p = Chem.SmilesWriteParams()
        try:
            res[key] = Chem.MolToCXSmiles(m, p, int(flags))
        except Exception as e:  # noqa: BLE001
            res[key] = type(e).__name__
    try:
        res["nonisomeric"] = Chem.MolToCXSmiles(m, isomericSmiles=False)
    except Exception as e:  # noqa: BLE001
        res["nonisomeric"] = type(e).__name__
    return res


def molhash_record(smi):
    """Every MolHash function of a molecule (the CXSMILES variants are not ported)."""
    from rdkit.Chem import rdMolHash
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return {f.name: rdMolHash.MolHash(m, f) for f in sorted(rdMolHash.HashFunction.values.values(), key=int)}


def erg_record(smi):
    """Reduced graph and ErG fingerprints of a molecule."""
    from rdkit import RDLogger
    from rdkit.Chem import rdReducedGraphs as rg
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    mrg = rg.GenerateMolExtendedReducedGraph(m)
    mrg.UpdatePropertyCache(False)
    rec = {
        "smiles": Chem.MolToSmiles(mrg),
        # (sorted: the order of the ring centroids follows the order of the rings,
        # which differs between RDKit releases)
        "types": sorted(list(a.GetPropsAsDict().get("_ErGAtomTypes", ())) for a in mrg.GetAtoms()),
        "fp": [round(float(x), 9) for x in rg.GetErGFingerprint(m)],
        "fp-nofuzz": [round(float(x), 9) for x in rg.GetErGFingerprint(m, fuzzIncrement=0)],
        "fp-short": [round(float(x), 9) for x in rg.GetErGFingerprint(m, minPath=2, maxPath=8)],
        "fp-graph": [round(float(x), 9) for x in rg.GenerateErGFingerprintForReducedGraph(mrg)],
    }
    return rec


LAYERED_VARIANTS = (
    {},
    {"layerFlags": 0x07},
    {"layerFlags": 0x38, "maxPath": 5, "fpSize": 512},
    {"layerFlags": 0x2A, "branchedPaths": False, "fromAtoms": [0, 1]},  # unrooted unbranched crashes upstream
    {"minPath": 2, "maxPath": 4, "fromAtoms": [0]},
)


def layeredfp_record(smi):
    """LayeredFingerprint (bits + atomCounts) of the SMILES- and SMARTS-parsed molecule."""
    from rdkit import DataStructs, RDLogger
    RDLogger.DisableLog("rdApp.*")
    rec = {}
    for kind, mol in (("smiles", Chem.MolFromSmiles(smi)), ("smarts", Chem.MolFromSmarts(smi))):
        if mol is None or mol.GetNumAtoms() > 40:
            continue
        out = []
        for kw in LAYERED_VARIANTS:
            counts = [0] * mol.GetNumAtoms()
            fp = Chem.LayeredFingerprint(mol, atomCounts=counts, **kw)
            # setOnlyBits: drop the first on-bit of the result and compute again
            mask = DataStructs.ExplicitBitVect(fp.GetNumBits())
            for bit in list(fp.GetOnBits())[1:]:
                mask.SetBit(bit)
            masked_counts = [0] * mol.GetNumAtoms()
            masked = Chem.LayeredFingerprint(mol, setOnlyBits=mask, atomCounts=masked_counts, **kw)
            out.append([list(masked.GetOnBits()), masked_counts])
            out.append([list(fp.GetOnBits()), counts])
        rec[kind] = out
    return rec or None


def rdkfingerprint_record(smi):
    """Legacy Chem.RDKFingerprint / UnfoldedRDKFingerprintCountBased."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    res = {}
    res["default"] = list(Chem.RDKFingerprint(m).GetOnBits())
    res["density"] = list(Chem.RDKFingerprint(m, maxPath=5, tgtDensity=0.3, minSize=64).GetOnBits())
    res["linear"] = list(Chem.RDKFingerprint(m, maxPath=6, nBitsPerHash=1, branchedPaths=False, fpSize=1024).GetOnBits())
    ab, bi = [], {}
    res["out"] = list(Chem.RDKFingerprint(m, maxPath=4, fromAtoms=[0], atomBits=ab, bitInfo=bi).GetOnBits())
    res["ab"] = ab
    res["bi"] = sorted(bi.items())
    ab, bi = [], {}
    u = Chem.UnfoldedRDKFingerprintCountBased(m, maxPath=5, atomBits=ab, bitInfo=bi)
    res["unfolded"] = sorted(u.GetNonzeroElements().items())
    res["uab"] = ab
    res["ubi"] = sorted(bi.items())
    return res


def moldescfp_record(smi):
    """Legacy rdMolDescriptors fingerprints (Morgan / atom pairs / torsions)."""
    from rdkit.Chem import rdMolDescriptors as D
    from rdkit.Chem import rdFingerprintGenerator as G
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    nz = lambda v: sorted(v.GetNonzeroElements().items())  # noqa: E731
    res = {}
    bi = {}
    res["morgan2"] = nz(D.GetMorganFingerprint(m, 2, bitInfo=bi))
    res["morgan2bi"] = sorted(bi.items())
    res["morgan2nc"] = nz(D.GetMorganFingerprint(m, 2, useCounts=False, useChirality=True))
    res["morgan3feat"] = nz(D.GetMorganFingerprint(m, 3, useFeatures=True))
    res["hmorgan"] = nz(D.GetHashedMorganFingerprint(m, 2, nBits=1024, useBondTypes=False))
    res["morganbv"] = list(D.GetMorganFingerprintAsBitVect(m, 2, nBits=512, fromAtoms=[0]).GetOnBits())
    res["conninv"] = D.GetConnectivityInvariants(m)
    res["featinv"] = D.GetFeatureInvariants(m)
    res["featgen"] = list(G.GetMorganGenerator(atomInvariantsGenerator=G.GetMorganFeatureAtomInvGen()).GetFingerprint(m).GetOnBits())
    res["ap"] = nz(D.GetAtomPairFingerprint(m))
    res["apchiral"] = nz(D.GetAtomPairFingerprint(m, includeChirality=True, maxLength=10))
    res["hap"] = nz(D.GetHashedAtomPairFingerprint(m, nBits=512))
    res["hapbv"] = list(D.GetHashedAtomPairFingerprintAsBitVect(m, ignoreAtoms=[0]).GetOnBits())
    res["hapbv3"] = list(D.GetHashedAtomPairFingerprintAsBitVect(m, nBitsPerEntry=3).GetOnBits())
    res["tt"] = nz(D.GetTopologicalTorsionFingerprint(m))
    res["ttchiral"] = nz(D.GetTopologicalTorsionFingerprint(m, includeChirality=True, fromAtoms=[0, 1]))
    res["htt"] = nz(D.GetHashedTopologicalTorsionFingerprint(m, nBits=512))
    res["httbv"] = list(D.GetHashedTopologicalTorsionFingerprintAsBitVect(m).GetOnBits())
    res["apcode"] = [D.GetAtomPairAtomCode(a) for a in m.GetAtoms()]
    return res


def props_record(obj):
    """All properties (private and computed), or the exception name."""
    try:
        d = obj.GetPropsAsDict(includePrivate=True, includeComputed=True)
    except Exception as e:  # noqa: BLE001
        return type(e).__name__
    return {k: d[k] for k in sorted(d)}


def mol_state_record(m):
    rec = {"smiles": Chem.MolToSmiles(m),
           "atoms": [atom_record(a) + [props_record(a)] for a in m.GetAtoms()],
           "bonds": [bond_record(b) + [props_record(b)] for b in m.GetBonds()],
           "mol": props_record(m)}
    if m.GetNumConformers():
        conf = m.GetConformer()
        rec["coords"] = [[round(p.x, 6), round(p.y, 6), round(p.z, 6)]
                         for p in (conf.GetAtomPosition(i) for i in range(m.GetNumAtoms()))]
    return rec


def addhs_records(m):
    """AddHs with the default, explicit-only, coordinate and atom-restricted options."""
    n = m.GetNumAtoms()
    only = [i for i in range(n) if i % 2 == 0]
    cases = (
        ("default", lambda mm: Chem.AddHs(mm)),
        ("explicitOnly", lambda mm: Chem.AddHs(mm, explicitOnly=True)),
        ("addCoords", lambda mm: Chem.AddHs(mm, addCoords=True)),
        ("onlyOnAtoms", lambda mm: Chem.AddHs(mm, onlyOnAtoms=only)),
        ("params", lambda mm: Chem.AddHs(mm, Chem.AddHsParameters())),
    )
    res = {}
    for key, fn in cases:
        try:
            res[key] = mol_state_record(fn(Chem.Mol(m)))
        except Exception as e:  # noqa: BLE001
            res[key] = type(e).__name__
    return res


def crippen_record(smi):
    from rdkit.Chem import rdMolDescriptors as D
    from rdkit.Chem import Crippen
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None

    def attempt(fn):
        try:
            return fn()
        except Exception as e:  # noqa: BLE001
            return type(e).__name__

    def contribs(mm, force=False):
        n = mm.GetNumAtoms()
        types, labels = [0] * n, [""] * n
        c = D._CalcCrippenContribs(mm, force, types, labels)
        return [[list(t) for t in c], types, labels]

    res = {}
    res["calc"] = attempt(lambda: list(D.CalcCrippenDescriptors(Chem.Mol(m))))
    res["calcNoH"] = attempt(lambda: list(D.CalcCrippenDescriptors(Chem.Mol(m), False)))
    res["contribs"] = attempt(lambda: contribs(Chem.Mol(m)))
    res["contribsNoForce"] = attempt(lambda: contribs(Chem.Mol(m), True))
    res["MolLogP"] = attempt(lambda: Crippen.MolLogP(Chem.Mol(m)))
    res["MolMR"] = attempt(lambda: Crippen.MolMR(Chem.Mol(m)))

    def cached():
        mm = Chem.Mol(m)
        D.CalcCrippenDescriptors(mm, False)
        return contribs(mm)
    res["cachedContribs"] = attempt(cached)

    def cached_calc():
        mm = Chem.Mol(m)
        contribs(mm)
        return list(D.CalcCrippenDescriptors(mm, False))
    res["cachedCalc"] = attempt(cached_calc)

    def forced():
        mm = Chem.Mol(m)
        D.CalcCrippenDescriptors(mm)
        return list(D.CalcCrippenDescriptors(mm, True, True))
    res["forced"] = attempt(forced)

    mh = Chem.AddHs(m)
    res["explicitH"] = attempt(lambda: list(D.CalcCrippenDescriptors(Chem.Mol(mh), False)))
    res["explicitHContribs"] = attempt(lambda: contribs(Chem.Mol(mh)))
    res["badTypes"] = attempt(lambda: D._CalcCrippenContribs(Chem.Mol(m), False, [0] * (m.GetNumAtoms() + 1)))
    return res


def molprops_record(smi):
    """MolWt/ExactMolWt/formula/atom counts (MolProps.cpp), also with explicit Hs, isotopes and charges."""
    from rdkit.Chem import rdMolDescriptors as D
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None

    def attempt(fn):
        try:
            return fn()
        except Exception as e:  # noqa: BLE001
            return type(e).__name__

    def record(mol):
        return {
            "amw": attempt(lambda: D._CalcMolWt(mol)),
            "amwHeavy": attempt(lambda: D._CalcMolWt(mol, True)),
            "exact": attempt(lambda: D.CalcExactMolWt(mol)),
            "exactHeavy": attempt(lambda: D.CalcExactMolWt(mol, onlyHeavy=True)),
            "formula": attempt(lambda: D.CalcMolFormula(mol)),
            "formulaIso": attempt(lambda: D.CalcMolFormula(mol, True)),
            "formulaIsoNoAbbrev": attempt(lambda: D.CalcMolFormula(mol, True, False)),
            "heavy": attempt(lambda: D.CalcNumHeavyAtoms(mol)),
            "atoms": attempt(lambda: D.CalcNumAtoms(mol)),
        }

    def isotopes(mol):
        mm = Chem.RWMol(mol)
        for a in mm.GetAtoms():
            if a.GetAtomicNum() == 1:
                a.SetIsotope(2 if a.GetIdx() % 2 else 3)
            elif a.GetIdx() == 0:
                a.SetIsotope(13)
        return mm

    def charged(mol):
        mm = Chem.RWMol(mol)
        if mm.GetNumAtoms():
            mm.GetAtomWithIdx(0).SetFormalCharge(mm.GetAtomWithIdx(0).GetFormalCharge() + 1)
            mm.UpdatePropertyCache(strict=False)
        return mm

    withHs = Chem.AddHs(m)
    return {
        "plain": record(m),
        "withHs": record(withHs),
        "isotopes": record(isotopes(withHs)),
        "charged": record(charged(m)),
    }


def molsurf_record(smi):
    """Labute ASA, TPSA and the VSA binnings (MolSurf.cpp), each on a fresh copy (they cache on the mol)."""
    from rdkit.Chem import rdMolDescriptors as D
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None

    def attempt(fn):
        try:
            return fn()
        except Exception as e:  # noqa: BLE001
            return type(e).__name__

    def fresh(fn):
        return attempt(lambda: fn(Chem.Mol(m)))

    def labute_contribs(mol):
        contribs, h = D._CalcLabuteASAContribs(mol)
        return [list(contribs), h]

    def charged_prop(mol):
        for a in mol.GetAtoms():
            a.SetDoubleProp("p", (a.GetIdx() % 7) / 3.0)
        return D.CustomProp_VSA_(mol, "p", [0.2, 0.9, 1.5])

    return {
        "labute": fresh(lambda x: D.CalcLabuteASA(x)),
        "labuteNoHs": fresh(lambda x: D.CalcLabuteASA(x, includeHs=False)),
        "labuteContribs": fresh(labute_contribs),
        "tpsa": fresh(lambda x: D.CalcTPSA(x)),
        "tpsaSP": fresh(lambda x: D.CalcTPSA(x, includeSandP=True)),
        "tpsaContribs": fresh(lambda x: list(D._CalcTPSAContribs(x))),
        "tpsaContribsSP": fresh(lambda x: list(D._CalcTPSAContribs(x, includeSandP=True))),
        "slogpVsa": fresh(lambda x: D.SlogP_VSA_(x)),
        "smrVsa": fresh(lambda x: D.SMR_VSA_(x)),
        "slogpVsaBins": fresh(lambda x: D.SlogP_VSA_(x, [0.5, -0.5, 0.0])),
        "customVsa": fresh(charged_prop),
        "tpsaWithHs": attempt(lambda: D.CalcTPSA(Chem.AddHs(m), includeSandP=True)),
    }


def gasteiger_record(smi):
    """Gasteiger charges (rdPartialCharges) with several iteration counts, and PEOE_VSA_."""
    from rdkit.Chem import rdMolDescriptors as D
    from rdkit.Chem import rdPartialCharges as P
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None

    def attempt(fn):
        try:
            return fn()
        except Exception as e:  # noqa: BLE001
            return type(e).__name__

    def charges(mol, *args, **kwargs):
        P.ComputeGasteigerCharges(mol, *args, **kwargs)
        return [[a.GetProp("_GasteigerCharge"), a.GetProp("_GasteigerHCharge"),
                 a.GetDoubleProp("_GasteigerCharge"), a.GetDoubleProp("_GasteigerHCharge")]
                for a in mol.GetAtoms()]

    def charged(mol):
        mm = Chem.RWMol(mol)
        if mm.GetNumAtoms():
            mm.GetAtomWithIdx(0).SetFormalCharge(mm.GetAtomWithIdx(0).GetFormalCharge() + 1)
            mm.UpdatePropertyCache(strict=False)
        return mm

    return {
        "default": attempt(lambda: charges(Chem.Mol(m))),
        "iter0": attempt(lambda: charges(Chem.Mol(m), 0)),
        "iter3": attempt(lambda: charges(Chem.Mol(m), 3)),
        "iter30": attempt(lambda: charges(Chem.Mol(m), nIter=30)),
        "throw": attempt(lambda: charges(Chem.Mol(m), 12, True)),
        "withHs": attempt(lambda: charges(Chem.AddHs(m))),
        "charged": attempt(lambda: charges(charged(m))),
        "peoeVsa": attempt(lambda: D.PEOE_VSA_(Chem.Mol(m))),
        "peoeVsaBins": attempt(lambda: D.PEOE_VSA_(Chem.Mol(m), [-0.1, 0.1])),
    }


def connectivity_record(smi):
    """Chi/Kappa/Hall-Kier/Phi (ConnectivityDescriptors.cpp), on the molecule and on its explicit-H form."""
    from rdkit.Chem import rdMolDescriptors as D
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None

    def attempt(fn):
        try:
            return fn()
        except Exception as e:  # noqa: BLE001
            return type(e).__name__

    def record(mol):
        res = {}
        for name in ("Chi0v", "Chi1v", "Chi2v", "Chi3v", "Chi4v", "Chi0n", "Chi1n", "Chi2n", "Chi3n", "Chi4n",
                     "HallKierAlpha", "Kappa1", "Kappa2", "Kappa3", "Phi"):
            res[name] = attempt(lambda: getattr(D, "Calc" + name)(Chem.Mol(mol)))
        res["ChiNv5"] = attempt(lambda: D.CalcChiNv(Chem.Mol(mol), 5))
        res["ChiNn1"] = attempt(lambda: D.CalcChiNn(Chem.Mol(mol), 1, True))
        res["ChiNnForce"] = attempt(lambda: D.CalcChiNn(Chem.Mol(mol), 3, force=True))

        def contribs():
            lst = [0.0] * mol.GetNumAtoms()
            total = D.CalcHallKierAlpha(Chem.Mol(mol), lst)
            return [total, lst]
        res["alphaContribs"] = attempt(contribs)
        return res

    return {"plain": record(m), "withHs": record(Chem.AddHs(m))}


def matrices_record(smi):
    """Distance/adjacency matrices (Matrices.cpp) as digests of the raw bytes, plus the cache property names."""
    import hashlib
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None

    def digest(fn):
        try:
            arr = fn()
            return [str(arr.dtype), list(arr.shape), hashlib.md5(arr.tobytes()).hexdigest()]
        except Exception as e:  # noqa: BLE001
            return type(e).__name__

    def record(mol):
        res = {}
        for useBO in (False, True):
            for wts in (False, True):
                res["dist%d%d" % (useBO, wts)] = digest(
                    lambda: Chem.GetDistanceMatrix(Chem.Mol(mol), useBO=useBO, useAtomWts=wts))
            res["adj%d" % useBO] = digest(lambda: Chem.GetAdjacencyMatrix(Chem.Mol(mol), useBO=useBO))
        res["adjEmpty"] = digest(lambda: Chem.GetAdjacencyMatrix(Chem.Mol(mol), emptyVal=-1))
        cached = Chem.Mol(mol)
        Chem.GetDistanceMatrix(cached)
        Chem.GetAdjacencyMatrix(cached, useBO=True)
        res["props"] = sorted(cached.GetPropNames(includePrivate=True, includeComputed=True))
        return res

    return {"plain": record(m), "withHs": record(Chem.AddHs(m))}


def transforms_record(smi):
    """DeleteSubstructs (ChemTransforms.cpp) and FindPotentialStereoBonds (Chirality.cpp)."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None

    def attempt(fn):
        try:
            return fn()
        except Exception as e:  # noqa: BLE001
            return type(e).__name__

    res = {}
    for name, sma in (("hetero", "[!#6;!#1]"), ("carbonyl", "C=O"), ("halogen", "[F,Cl,Br,I]"), ("ring", "[R]"),
                      ("chiral", "[C@,C@@]"), ("methyl", "[CH3]"), ("arom", "a")):
        q = Chem.MolFromSmarts(sma)
        for opts in ((False, False), (True, False), (False, True)):
            res["del-%s-%d%d" % ((name,) + tuple(map(int, opts)))] = attempt(
                lambda: Chem.MolToSmiles(Chem.DeleteSubstructs(m, q, onlyFrags=opts[0], useChirality=opts[1])))

    def potential(clean, hs=False):
        mm = Chem.AddHs(m) if hs else Chem.Mol(m)
        Chem.FindPotentialStereoBonds(mm, cleanIt=clean)
        return [[b.GetIdx(), str(b.GetStereo()), list(b.GetStereoAtoms())] for b in mm.GetBonds()
                if b.GetBondType() == Chem.BondType.DOUBLE] + [sorted(mm.GetPropNames(True, True))]
    res["potential"] = attempt(lambda: potential(False))
    res["potentialClean"] = attempt(lambda: potential(True))
    res["potentialHs"] = attempt(lambda: potential(False, hs=True))
    return res


def descriptors_record(smi):
    """Every descriptor of rdkit.Chem.Descriptors (minus SPS, which needs the new stereo perception)."""
    from rdkit.Chem import Descriptors
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    skip = {"SPS"}
    res = {}
    for name, fn in Descriptors._descList:
        if name in skip:
            continue
        try:
            res[name] = fn(Chem.Mol(m))
            if name.startswith("BCUT2D"):
                # eigenvalues come from different solvers (Eigen vs rdoxide's QL): compare to 9 decimals
                res[name] = round(res[name], 9) + 0.0
        except Exception as e:  # noqa: BLE001
            res[name] = type(e).__name__
    return res


def lipinski_record(smi):
    """rdMolDescriptors counts (Lipinski.cpp), with the rotatable-bond modes and the atom lists."""
    from rdkit.Chem import rdMolDescriptors as D
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None

    def attempt(fn):
        try:
            return fn()
        except Exception as e:  # noqa: BLE001
            return type(e).__name__

    E = D.NumRotatableBondsOptions
    res = {}
    for name in ("CalcNumLipinskiHBA", "CalcNumLipinskiHBD", "CalcNumHBD", "CalcNumHBA", "CalcNumHeteroatoms",
                 "CalcNumAmideBonds", "CalcNumRings", "CalcFractionCSP3", "CalcNumHeterocycles",
                 "CalcNumAromaticRings", "CalcNumSaturatedRings", "CalcNumAliphaticRings",
                 "CalcNumAromaticHeterocycles", "CalcNumAromaticCarbocycles", "CalcNumAliphaticHeterocycles",
                 "CalcNumAliphaticCarbocycles", "CalcNumSaturatedHeterocycles", "CalcNumSaturatedCarbocycles",
                 "CalcNumAtomStereoCenters", "CalcNumUnspecifiedAtomStereoCenters"):
        res[name] = attempt(lambda: getattr(D, name)(Chem.Mol(m)))
    res["rotDefault"] = attempt(lambda: D.CalcNumRotatableBonds(Chem.Mol(m)))
    for key, opt in (("rotNonStrict", E.NonStrict), ("rotStrict", E.Strict), ("rotLinkages", E.StrictLinkages),
                     ("rotBoolTrue", True), ("rotBoolFalse", False)):
        res[key] = attempt(lambda: D.CalcNumRotatableBonds(Chem.Mol(m), opt))

    def with_list(fn):
        atoms = []
        n = fn(Chem.Mol(m), atoms)
        return [n, atoms]
    res["spiro"] = attempt(lambda: with_list(D.CalcNumSpiroAtoms))
    res["bridgehead"] = attempt(lambda: with_list(D.CalcNumBridgeheadAtoms))

    def unassigned():
        mm = Chem.Mol(m)
        mm.ClearProp("_StereochemDone")
        return [D.CalcNumAtomStereoCenters(mm), D.CalcNumUnspecifiedAtomStereoCenters(mm)]
    res["stereoUnassigned"] = attempt(unassigned)
    return res


def standardize_record(smi):
    """rdMolStandardize results; RDOXIDE_STANDARDIZE (comma-separated keys)
    limits the comparison to the steps ported so far."""
    import os

    from rdkit.Chem.MolStandardize import rdMolStandardize as ms
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    steps = {
        "removeFragments": lambda x: ms.RemoveFragments(x),
        "fragmentRemover": lambda x: ms.FragmentRemover().remove(x),
        "fragmentRemoverSkip": lambda x: ms.FragmentRemover(skip_if_all_match=True).remove(x),
        "largestFragment": lambda x: ms.LargestFragmentChooser().choose(x),
        "largestOrganic": lambda x: ms.LargestFragmentChooser(preferOrganic=True).choose(x),
        "largestHeavy": lambda x: ms.LargestFragmentChooser(_heavy_params()).choose(x),
        "normalize": lambda x: ms.Normalize(x),
        "normalizeInPlace": lambda x: (ms.NormalizeInPlace(x), x)[1],
        "reionize": lambda x: ms.Reionize(x),
        "uncharge": lambda x: ms.Uncharger().uncharge(x),
        "unchargeNoncanon": lambda x: ms.Uncharger(canonicalOrder=False).uncharge(x),
        "disconnect": lambda x: ms.DisconnectOrganometallics(x),
        "metalDisconnector": lambda x: ms.MetalDisconnector().Disconnect(x),
        "metalDisconnectorAll": lambda x: ms.MetalDisconnector(_metal_opts()).Disconnect(x),
        "cleanup": lambda x: ms.Cleanup(x),
        "fragmentParent": lambda x: ms.FragmentParent(x),
        "chargeParent": lambda x: ms.ChargeParent(x),
        "isotopeParent": lambda x: ms.IsotopeParent(x),
        "stereoParent": lambda x: ms.StereoParent(x),
        "canonicalTautomer": lambda x: ms.CanonicalTautomer(x),
        "tautomerParent": lambda x: ms.TautomerParent(x),
        "superParent": lambda x: ms.SuperParent(x),
        "standardizeSmiles": lambda x: ms.StandardizeSmiles(Chem.MolToSmiles(x)),
        "validate": lambda x: ms.ValidateSmiles(Chem.MolToSmiles(x)),
        "validators": lambda x: _validators(ms, x),
        "tautomers": lambda x: _tautomers(ms, x),
    }
    keys = os.environ.get("RDOXIDE_STANDARDIZE")
    keys = keys.split(",") if keys else list(steps)
    res = {}
    for k in keys:
        try:
            out = steps[k](Chem.Mol(m))
        except Exception as e:  # noqa: BLE001
            res[k] = type(e).__name__
            continue
        res[k] = Chem.MolToSmiles(out) if isinstance(out, Chem.Mol) else out
    return res


def _validators(ms, m):
    from rdkit.Chem import rdDepictor
    res = {}
    vals = {
        "rdkit": ms.RDKitValidation(), "molvs": ms.MolVSValidation(), "neutral": ms.NeutralValidation(),
        "isotope": ms.IsotopeValidation(), "isotopeStrict": ms.IsotopeValidation(strict=True),
        "radical": ms.DisallowedRadicalValidation(), "features": ms.FeaturesValidation(),
        "allowed": ms.AllowedAtomsValidation([Chem.Atom(x) for x in (6, 7, 8, 1)]),
        "disallowed": ms.DisallowedAtomsValidation([Chem.Atom(x) for x in (9, 17, 35)]),
        "is2d": ms.Is2DValidation(), "layout": ms.Layout2DValidation(), "stereo": ms.StereoValidation(),
    }
    m2 = Chem.Mol(m)
    rdDepictor.Compute2DCoords(m2)
    Chem.WedgeMolBonds(m2, m2.GetConformer())
    for k, v in vals.items():
        for all_ in (False, True):
            for tag, mol in (("", m), ("2d", m2)):
                try:
                    res[f"{k}{tag}{int(all_)}"] = v.validate(mol, reportAllFailures=all_)
                except Exception as e:  # noqa: BLE001
                    res[f"{k}{tag}{int(all_)}"] = type(e).__name__
    return res


def _metal_opts():
    from rdkit.Chem.MolStandardize import rdMolStandardize as ms
    o = ms.MetalDisconnectorOptions()
    o.splitGrignards = True
    o.splitAromaticC = True
    return o


def _heavy_params():
    from rdkit.Chem.MolStandardize import rdMolStandardize as ms
    p = ms.CleanupParameters()
    p.largestFragmentChooserCountHeavyAtomsOnly = True
    return p


def _tautomers(ms, m):
    te = ms.TautomerEnumerator()
    te.SetMaxTautomers(50)
    res = te.Enumerate(m)
    return {"status": str(res.status), "smiles": list(res.smiles),
            "canonical": Chem.MolToSmiles(te.Canonicalize(m)),
            "score": te.ScoreTautomer(m)}


def dump(mode, smi):
    if mode == "standardize":
        return standardize_record(smi)
    if mode == "lipinski":
        return lipinski_record(smi)
    if mode == "molprops":
        return molprops_record(smi)
    if mode == "molsurf":
        return molsurf_record(smi)
    if mode == "gasteiger":
        return gasteiger_record(smi)
    if mode == "connectivity":
        return connectivity_record(smi)
    if mode == "matrices":
        return matrices_record(smi)
    if mode == "transforms":
        return transforms_record(smi)
    if mode == "descriptors":
        return descriptors_record(smi)
    if mode == "crippen":
        return crippen_record(smi)
    if mode == "addhs":
        m = Chem.MolFromSmiles(smi)
        return None if m is None else addhs_records(m)
    if mode == "addhs3d":
        m = Chem.MolFromMolBlock(smi)
        return None if m is None else addhs_records(m)
    if mode == "erg":
        return erg_record(smi)
    if mode == "molhash":
        return molhash_record(smi)
    if mode == "layeredfp":
        return layeredfp_record(smi)
    if mode == "patternfp":
        return patternfp_record(smi)
    if mode == "maccs":
        from rdkit.Chem import MACCSkeys
        m = Chem.MolFromSmiles(smi)
        return None if m is None else list(MACCSkeys.GenMACCSKeys(m).GetOnBits())
    if mode == "moldescfp":
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
        return moldescfp_record(smi)
    if mode == "rdkfingerprint":
        return rdkfingerprint_record(smi)
    if mode == "rdkitfp":
        return rdkitfp_record(smi)
    if mode == "torsion":
        return torsion_record(smi)
    if mode == "subgraphs":
        return subgraphs_record(smi)
    if mode == "morgan":
        return morgan_record(smi)
    if mode == "atompair":
        return atompair_record(smi)
    if mode == "depict":
        from rdkit.Chem import rdDepictor
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        res = {}
        for key, kw in (("default", {}), ("nocanon", {"canonOrient": False})):
            mm = Chem.Mol(m)
            try:
                rdDepictor.Compute2DCoords(mm, **kw)
                conf = mm.GetConformer()
                res[key] = [[round(conf.GetAtomPosition(i).x, 6), round(conf.GetAtomPosition(i).y, 6)] for i in range(mm.GetNumAtoms())]
            except Exception as e:  # noqa: BLE001
                res[key] = type(e).__name__
        try:
            res["block"] = Chem.MolToMolBlock(m)
        except Exception as e:  # noqa: BLE001
            res["block"] = type(e).__name__
        return res
    if mode == "molwrite":
        return molwrite_record(smi)
    if mode == "molblock":
        res = {}
        for key, kw in (("default", {}), ("nosanitize", {"sanitize": False}),
                        ("keepHs", {"removeHs": False}), ("lax", {"strictParsing": False})):
            try:
                res[key] = molblock_record(smi, **kw)
            except Exception as e:  # noqa: BLE001
                res[key] = {"exception": type(e).__name__}
        return res
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
    if mode == "smarts":
        m = Chem.MolFromSmarts(smi)
        if m is None:
            return None
        atoms = [atom_record(a) + [a.GetQueryType(), a.DescribeQuery(), a.GetAtomMapNum()]
                 for a in m.GetAtoms()]
        bonds = [bond_record(b) + [b.DescribeQuery()] for b in m.GetBonds()]
        out = {"atoms": atoms, "bonds": bonds}

        def run(key, f):
            try:
                out[key] = f()
            except Exception as e:  # noqa: BLE001
                out[key] = "EXC:" + type(e).__name__

        run("smarts", lambda: Chem.MolToSmarts(m))
        run("order", lambda: list(m.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"]))
        run("noniso", lambda: Chem.MolToSmarts(m, isomericSmiles=False))
        run("rooted", lambda: Chem.MolToSmarts(m, rootedAtAtom=m.GetNumAtoms() - 1))
        run("fragment", lambda: Chem.MolFragmentToSmarts(m, atomsToUse=list(range(max(1, m.GetNumAtoms() // 2)))))
        run("atom_smarts", lambda: [a.GetSmarts() for a in m.GetAtoms()])
        run("bond_smarts", lambda: [b.GetSmarts() for b in m.GetBonds()])
        run("roundtrip", lambda: Chem.MolToSmarts(Chem.MolFromSmarts(Chem.MolToSmarts(m))))
        return out
    if mode in ("pickle", "picklesmarts"):
        m = Chem.MolFromSmiles(smi) if mode == "pickle" else Chem.MolFromSmarts(smi)
        if m is None:
            return None
        P = Chem.PropertyPickleOptions
        out = {}
        for key, flags in (("default", None), ("all", P.AllProps),
                           ("all_double", P.AllProps | P.CoordsAsDouble),
                           ("noconf_mol", P.MolProps | P.NoConformers)):
            try:
                b = m.ToBinary() if flags is None else m.ToBinary(flags)
                back = Chem.Mol(b)
                out[key] = [b.hex(), back.GetNumAtoms(), back.GetNumBonds(),
                            Chem.MolToSmiles(back) if mode == "pickle" else Chem.MolToSmarts(back),
                            sorted(back.GetPropsAsDict(True, True).items(), key=str),
                            [sorted(a.GetPropsAsDict(True, True).items(), key=str) for a in back.GetAtoms()],
                            list(map(list, back.GetRingInfo().AtomRings()))]
            except Exception as e:  # noqa: BLE001
                out[key] = "EXC:" + type(e).__name__ + ":" + str(e)[:80]
        return out
    if mode == "substruct":
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        out = []
        qs = GENERIC_QUERIES + [QUERIES[(QUERY_OFFSET + k) % len(QUERIES)] for k in range(25)]
        for qsma in qs:
            q = Chem.MolFromSmarts(qsma)
            if q is None:
                continue
            rec = [qsma]
            for f in (lambda: m.GetSubstructMatches(q),
                      lambda: m.GetSubstructMatches(q, uniquify=False, maxMatches=50),
                      lambda: m.GetSubstructMatches(q, useChirality=True),
                      lambda: m.HasSubstructMatch(q)):
                try:
                    rec.append(f())
                except Exception as e:  # noqa: BLE001
                    rec.append("EXC:" + type(e).__name__)
            out.append(rec)
        # chirality: the molecule against itself and against its mirror image
        mirror = Chem.Mol(m)
        for a in mirror.GetAtoms():
            if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            elif a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        for q in (m, mirror):
            try:
                out.append([m.GetSubstructMatches(q, useChirality=True, uniquify=False, maxMatches=20),
                            m.HasSubstructMatch(q, useChirality=True)])
            except Exception as e:  # noqa: BLE001
                out.append("EXC:" + type(e).__name__)
        return out
    if mode == "smiles":
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        out = {}

        def run(key, f):
            try:
                out[key] = f()
            except Exception as e:  # noqa: BLE001
                out[key] = "EXC:" + type(e).__name__

        run("canonical", lambda: Chem.MolToSmiles(m))
        run("order", lambda: list(m.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"]))
        run("bond_order", lambda: list(m.GetPropsAsDict(True, True)["_smilesBondOutputOrder"]))
        run("noncanonical", lambda: Chem.MolToSmiles(m, canonical=False))
        run("kekule", lambda: Chem.MolToSmiles(Chem.Mol(m), kekuleSmiles=True))
        run("noniso", lambda: Chem.MolToSmiles(m, isomericSmiles=False))
        run("explicit", lambda: Chem.MolToSmiles(m, allBondsExplicit=True, allHsExplicit=True))
        run("rooted", lambda: Chem.MolToSmiles(m, rootedAtAtom=m.GetNumAtoms() - 1))
        run("random", lambda: Chem.MolToRandomSmilesVect(m, 3, randomSeed=42))
        n = m.GetNumAtoms()
        run("fragment", lambda: Chem.MolFragmentToSmiles(m, atomsToUse=list(range(max(1, n // 2)))))
        run("fragment_kekule", lambda: Chem.MolFragmentToSmiles(
            Chem.Mol(m), atomsToUse=list(range(n)), kekuleSmiles=True, canonical=False))
        run("roundtrip", lambda: Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(m))))
        run("tosmarts", lambda: Chem.MolToSmarts(m))
        run("atom_smarts", lambda: [a.GetSmarts() for a in m.GetAtoms()])
        return out
    raise SystemExit(f"unknown mode {mode}")


QUERIES = []
QUERY_OFFSET = 0
GENERIC_QUERIES = [
    "*", "C", "c", "[#6]~[#6]", "C=O", "[OH]", "[#7]", "c1ccccc1", "C1CCCCC1", "[R]", "[R2]",
    "[!R]~[R]", "*~*~*", "[#6](~[#8])~[#8]", "[$(C=O)]", "[$([#6]~[#7])]~*", "[a;r6]", "[x3]",
    "C-,=C", "[C,c]:,-[N,n]", "**(*)*", "[D3]", "[#6;!$([#6]=O)]", "c:c", "C@C", "C!@C",
    "[N,O;!H0]", "*1~*~*~1", "[$(*~[#8]),$(*~[#7])]", "[#6]=,:[#6]-[#6]",
]


def main():
    global QUERY_OFFSET
    mode, path = sys.argv[1], sys.argv[2]
    if mode == "substruct":
        import os
        qpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tests", "data",
                             "substruct_queries.smarts")
        QUERIES.extend(l.split()[0] for l in open(qpath) if l.strip())
    if mode == "sdwrite":
        for rec in sdwrite_records(path):
            print(json.dumps(rec, default=list))
        return
    if mode == "sdf":
        for rec in sdf_records(path):
            print(json.dumps(rec, default=list))
        return
    if mode == "cxsmiles":
        for line in open(path):
            line = line.strip()
            if not line or line.lower().startswith("smiles"):
                continue
            try:
                res = cx_record(Chem.MolFromSmiles(line))
            except Exception as e:  # noqa: BLE001
                res = {"exception": type(e).__name__}
            print(json.dumps([line, res], default=list))
        return
    if mode == "cxmolblock":
        txt = open(path, encoding="utf-8").read()
        for i, blk in enumerate(txt.split("$$$$\n")):
            if not blk.strip():
                continue
            try:
                res = cx_record(Chem.MolFromMolBlock(blk))
            except Exception as e:  # noqa: BLE001
                res = {"exception": type(e).__name__}
            print(json.dumps([f"#{i} {blk.splitlines()[0][:40] if blk else ''}", res], default=list))
        return
    if mode in ("molblock", "molwrite", "addhs3d"):
        txt = open(path, encoding="utf-8").read()
        for i, blk in enumerate(txt.split("$$$$\n")):
            if not blk.strip():
                continue
            try:
                res = dump(mode, blk)
            except Exception as e:  # noqa: BLE001
                res = {"exception": type(e).__name__}
            print(json.dumps([f"#{i} {blk.splitlines()[0][:40] if blk else ''}", res], default=list))
        return
    with open(path) as fh:
        for line in fh:
            QUERY_OFFSET = (QUERY_OFFSET + 7) % max(1, len(QUERIES) or 1)
            parts = line.split()
            if not parts or parts[0].lower() == "smiles":
                continue
            smi = parts[0]
            try:
                res = dump(mode, smi)
            except Exception as e:  # noqa: BLE001
                res = {"exception": type(e).__name__}
            print(json.dumps([smi, res], default=list))


if __name__ == "__main__":
    main()
