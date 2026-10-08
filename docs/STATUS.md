# Port status

Upstream: RDKit 2026.09.1. Legend: ✅ done · 🟡 partial · ⬜ not started.

| Stage | Component | Status | Notes |
|-------|-----------|--------|-------|
| 0 | Workspace, packaging (maturin, abi3-py310), CI, licenses, NOTICE | ✅ | |
| 0 | `tools/sync_upstream.py` (fetch/regen/diff/stamp), `porting/manifest.toml` | ✅ | |
| 0 | Upstream test ratchet (`tools/run_upstream_tests.py`) | ✅ | 58 upstream tests passing |
| 1 | RDGeneral: Dict/RDValue/RDProps | ✅ | |
| 1 | RDGeneral: RDLog + `rdBase` logging API | ✅ | upstream `UnitTestLogging.py` 32/32 |
| 1 | RDGeneral: RNG (`boost::minstd_rand`, glibc `rand`) | ✅ | |
| 1 | `rdBase` vector converters (`_vecti`, ...) | ⬜ | come with stage 2 return types |
| 2 | PeriodicTable (+ data extraction from upstream) | ✅ | |
| 2 | Atom, Bond, Conformer, ROMol/RWMol, RingInfo, StereoGroup (core) | ✅ | boost adjacency ordering, RWMol edits, batch edits, valence calc |
| 2 | `rdchem` Python: Atom, Bond, Mol, RWMol, EditableMol, Conformer, RingInfo, StereoGroup, enums, props API, sanitize exceptions | ✅ | `_wrapperType = 'nanobind'` (PyO3 semantics match nanobind's) |
| 2 | `rdGeometry` Point3D / Point2D | 🟡 | PointND, UniformGrid3D later |
| 2b | SubstanceGroup, MonomerInfo / PDBResidueInfo, bookmarks API in Python, QueryAtom/QueryBond | ⬜ | queries with stage 4 |
| 3 | Ring finding: RingDecomposerLib port, findSSSR, symmetrizeSSSR (RDL + legacy), fastFindRings, ring families | ✅ | differential: 10312/10312 identical to RDKit 2026.09.1 |
| 3 | sanitizeMol: cleanUp, cleanUpOrganometallics, Kekulize (canonical + input order), aromaticity (RDKit/simple/MDL/MMFF94), conjugation, hybridization, radicals, adjustHs, chirality/atropisomer cleanup | ✅ | differential `sanitize` 10312/10312 |
| 3 | Canonical atom ranking (new_canon: rankMolAtoms, rankFragmentAtoms, chiralRankMolAtoms); `CanonicalRankAtoms(InFragment)` | ✅ | differential `ranks` 10312/10312 |
| 3 | removeHs/removeAllHs (`RemoveHsParameters`), legacy assignStereochemistry (CIP ranks/codes, ring stereo, E/Z), setDoubleBondNeighborDirections | ✅ | differential `molfromsmiles` 10312/10312 + 35140 stereo-rich + 72 edge cases; SubstanceGroup checks wait for stage 2b |
| 3 | Stereo from coordinates: assignChiralTypesFromBondDirs (pseudo-3D wedges), assignChiralTypesFrom3D (incl. non-tetrahedral), detectBondStereochemistry, atropisomer perception (2D/3D/no conformer) | ✅ | differential `molblock` |
| 3 | addHs, new (non-legacy) stereo perception, CIPLabeler | ⬜ | `SetUseLegacyStereoPerception(False)` currently raises |
| 4 | SMILES parser (lexer + grammar + SmilesParseOps, non-tetrahedral stereo) | ✅ | differential: 10312/10312 identical to RDKit (sanitize=False) |
| 4 | `MolFromSmiles` sanitize/removeHs/stereo path (defaults) | ✅ | differential `molfromsmiles` identical |
| 4 | SMILES writer + canonical SMILES (Canon.cpp DFS, SmilesWrite.cpp, getMolFrags/copyMolSubset): `MolToSmiles` (both overloads), `MolFragmentToSmiles`, `MolToRandomSmilesVect`, `SmilesWriteParams` | ✅ | differential `smiles` (canonical, output order, non-canonical, kekule, non-isomeric, explicit, rooted, random, fragment, roundtrip) 10312 + 35973 identical |
| 4 | Query model (Query/*, QueryOps, QueryAtom/QueryBond expandQuery, null-query algebra), SMARTS parser (`MolFromSmarts`, `AtomFromSmarts`, `BondFromSmarts`, `SmartsParserParams`), `DescribeQuery` | ✅ | differential `smarts` (atoms, bonds, query descriptions) 1275 library SMARTS + 151 edge cases identical |
| 4 | SMARTS writer: `MolToSmarts` (both overloads), `MolFragmentToSmarts`, `Atom/Bond.GetSmarts` | ✅ | differential `smarts` (+ writer: default, non-isomeric, rooted, fragment, per-atom/bond, roundtrip) 1426/1426; `smiles` mode also checks `MolToSmarts` of plain molecules |
| 4 | `mergeQueryHs`, Python `QueryAtom`/`QueryBond` classes, `rdqueries` | ⬜ | `mergeHs=True` currently returns None |
| 4 | CXSMILES (read/write) | ⬜ | |
| 5 | Substructure matching: VF2, `SubstructMatch`/`SubstructMatchCount`, recursive queries, chirality + double-bond stereo, enhanced stereo, query-query matches, aromaticMatchesConjugated/SingleOrDouble, atom/bond properties; `Mol.HasSubstructMatch/GetSubstructMatch/GetSubstructMatches` (both overloads), `SubstructMatchParameters` | ✅ | differential `substruct` (55 queries x 4 variants per molecule + chirality self/mirror) 10312 + 2418 identical |
| 5 | `extraFinalCheck`/`extraAtomCheck`/`extraBondCheck` callbacks, generic matchers, MolBundle, ResonanceMolSupplier | ⬜ | |
| 6 | MolPickler (byte-identical pickles; reads versions >= 2 incl. queries, rings, stereo groups, conformers, props): `Mol.ToBinary`, `Mol(pkl[, flags])`, Python `pickle`, `Get/SetDefaultPickleProperties`, `PropertyPickleOptions` | ✅ | differential `pickle`/`picklesmarts` byte-identical on 10312 + 1101 + 1426 + edge cases; upstream canonSmiles.v2.pkl stream reads |
| 6 | Pickles with SubstanceGroups / monomer info, v1 pickles, ExplicitBitVect property handler | ⬜ | raise RuntimeError for now |
| 6 | MolBlock V2000 reader: atoms/bonds incl. query atoms (A/Q/X/M/R/*, generic groups, H count), query bonds (5-8, topology), M CHG/RAD/ISO/RGP/ALS/RBC/SUB/UNS/PXA, old atom lists, aliases/values, S SKP, 2D/3D stereo; `MolFromMolBlock`, `MolFromMolFile` | ✅ | differential `molblock` (4 parse variants, byte-identical pickles with all props): 2D 4990/4991, 3D 490/490, upstream test_data 309/364 (all diffs are the unported features below) |
| 6 | V3000 CTAB, SGroups, M ZBO/ZCH/HYD/MRV/APO/LIN, SDMolSupplier, MolToMolBlock, SDWriter | ⬜ | parse as failures (None) for now; next |
| 7–13 | see PLAN.md | ⬜ | |

## Python package files

| File | Kind |
|------|------|
| `rdkit/__init__.py`, `RDLogger.py`, `RDRandom.py`, `RDConfig.py` | verbatim upstream (porting/copies.toml) |
| `rdkit/RDPaths.py` | rdoxide's version of the build-generated file (package-relative paths, like the PyPI wheels) |
| `rdkit/rdBase.py`, `rdkit/Chem/rdchem.py`, `rdkit/Geometry/rdGeometry.py` | shims → `rdkit._rdoxide` |
| `rdkit/Chem/__init__.py`, `rdkit/Geometry/__init__.py` | provisional (until rdmolfiles/rdmolops/DataStructs exist) |

## Iteration log

* 2026-10-07 — Stage 0 + 1 + PeriodicTable.
* 2026-10-07 — Stage 2: core graph + rdchem bindings. Next: stage 3/4 critical path —
  SMILES parser (no sanitize) → FindRings (SSSR/symm) → Kekulize/Aromaticity/ConjugHybrid →
  sanitizeMol → SMILES writer + canonical ranking. Then rough_test.py starts importing.
* 2026-10-07 — Stage 3a–c: rings, sanitizeMol, canonical ranking.
* 2026-10-07 — Stage 3d: removeHs + legacy stereo perception; default `MolFromSmiles`
  is now identical to RDKit. Next: SMILES writer (SmilesWrite.cpp, Canon.cpp
  canonicalizeFragment, MolToSmiles/MolFragmentToSmiles) → rough_test.py imports.
* 2026-10-08 — Stage 4b: SMILES writer identical to RDKit. Next: SMARTS + query atoms/bonds
  (needed by rough_test.py imports and substructure search), then CXSMILES.
* 2026-10-08 — Stage 4c (part 1): query model + SMARTS parser identical to RDKit.
* 2026-10-08 — Stage 4c (part 2): SMARTS writer identical to RDKit. Next: substructure matching
  (stage 5: vf2, SubstructMatch, HasSubstructMatch/GetSubstructMatch(es), recursive queries),
  then mergeQueryHs + QueryAtom classes.
* 2026-10-08 — Stage 5 (core): substructure matching identical to RDKit (match order included).
  Next: what rough_test.py needs to import (Chem/__init__ pieces), then MolPickler and file
  formats (stage 6).
* 2026-10-08 — RDConfig/RDPaths added. Remaining upstream-test blockers by count: DataStructs
  (~19 tests), MolBlock/SDMolSupplier/MolToMolBlock (~22), AllChem (11), pickling (3).
  Next: stage 6 (MolPickler, then mol file reader/writer, SDMolSupplier).
* 2026-10-08 — Stage 6a: MolPickler byte-identical to RDKit; Python pickling. Ratchet 53 -> 56.
  Next: MolBlock reader/writer + SDMolSupplier (stage 6b).
* 2026-10-08 — Stage 6b (part 1): V2000 MolBlock reader + stereo from coordinates identical to
  RDKit (`molblock` differential). Ratchet 56 -> 58. Next: V3000 + SGroups, SDMolSupplier,
  MolToMolBlock.

## Differential testing

`tools/differential/compare.py MODE` runs the same dump script under the
official RDKit of the ported version and under rdoxide (`.venv`) on the
upstream SMILES corpora (NCI first_5K, canonSmiles.long) and reports
differences. The reference interpreter comes from conda-forge
(`micromamba create -p REF -c conda-forge python=3.13 rdkit=2026.09.1`) and is
selected with `--ref-python` or `RDOXIDE_REF_PYTHON`. (PyPI only has
2026.03.x, which predates the RDL-based SSSR, so it is not a valid oracle for
ring order.)

Current results: `parse`, `rings`, `sanitize`, `ranks`, `molfromsmiles`, `smiles`: 10312/10312 each;
`smarts`: 1426/1426 (RDKit's SMARTS libraries + `tests/data/smarts_edge_cases.smarts`);
`substruct`: 10312/10312 + 2418/2418 stereo-rich (queries in `tests/data/substruct_queries.smarts`).
`molblock`: corpora built with `tools/differential/make_molblock_corpus.py` under the reference
interpreter (NCI 2D with wedges/explicit Hs, embedded 3D, upstream FileParsers test_data).
`molfromsmiles` (atoms, bonds, stereo props `_CIPCode`/`_CIPRank`/`_ringStereoAtoms`/
`_ChiralityPossible`) also matches on 35140 molecules from the stereo-rich upstream
corpora (zinc, ChEMBL, idorsia enumeration, MMFF) and on `tests/data/stereo_edge_cases.smi`.
New Python test files are also run against the official wheel to validate expectations.

## Known deviations (to revisit)

* Atom/Bond objects returned from a Mol hold a reference to the Python Mol plus
  an index (upstream holds a C++ pointer); after removing atoms the index may
  point at a different atom (upstream: dangling pointer).
* `Atom.IsInRing()` needs stage 3 (`findSSSR`).
* `MolFragmentToSmiles(mol, params, atomsToUse)` with `params` passed positionally works
  in rdoxide; upstream's Boost overload resolution raises AttributeError for it (keyword
  form works in both).
* `GetPropsAsDict(includePrivate=True, includeComputed=True)` returns `__computedProps`
  as a list (upstream: a `_vectSt6string` wrapper).
* `Chem.SetAllowNontetrahedralChirality` & co. keep the value in-process;
  upstream also writes it to the process environment.
