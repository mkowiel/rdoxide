# Port status

Upstream: RDKit 2026.09.1. Legend: ✅ done · 🟡 partial · ⬜ not started.

| Stage | Component | Status | Notes |
|-------|-----------|--------|-------|
| 0 | Workspace, packaging (maturin, abi3-py310), CI, licenses, NOTICE | ✅ | |
| 0 | `tools/sync_upstream.py` (fetch/regen/diff/stamp), `porting/manifest.toml` | ✅ | |
| 0 | Upstream test ratchet (`tools/run_upstream_tests.py`) | ✅ | 34 upstream tests passing |
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
| 3 | Kekulize, aromaticity, conjugation/hybridization, cleanUp, sanitizeMol, AddHs/RemoveHs, stereo perception | ⬜ | next |
| 4 | SMILES parser (lexer + grammar + SmilesParseOps, non-tetrahedral stereo) | ✅ | differential: 10312/10312 identical to RDKit (sanitize=False) |
| 4 | `MolFromSmiles` sanitize/removeHs path | ⬜ | needs stage 3 |
| 4 | CXSMILES, SMARTS, query atoms/bonds (`~`), SMILES writer, canonicalization | ⬜ | |
| 5–13 | see PLAN.md | ⬜ | |

## Python package files

| File | Kind |
|------|------|
| `rdkit/__init__.py`, `RDLogger.py`, `RDRandom.py` | verbatim upstream (porting/copies.toml) |
| `rdkit/rdBase.py`, `rdkit/Chem/rdchem.py`, `rdkit/Geometry/rdGeometry.py` | shims → `rdkit._rdoxide` |
| `rdkit/Chem/__init__.py`, `rdkit/Geometry/__init__.py` | provisional (until rdmolfiles/rdmolops/DataStructs exist) |

## Iteration log

* 2026-10-07 — Stage 0 + 1 + PeriodicTable.
* 2026-10-07 — Stage 2: core graph + rdchem bindings. Next: stage 3/4 critical path —
  SMILES parser (no sanitize) → FindRings (SSSR/symm) → Kekulize/Aromaticity/ConjugHybrid →
  sanitizeMol → SMILES writer + canonical ranking. Then rough_test.py starts importing.

## Differential testing

`tools/differential/compare.py MODE` runs the same dump script under the
official RDKit of the ported version and under rdoxide (`.venv`) on the
upstream SMILES corpora (NCI first_5K, canonSmiles.long) and reports
differences. The reference interpreter comes from conda-forge
(`micromamba create -p REF -c conda-forge python=3.13 rdkit=2026.09.1`) and is
selected with `--ref-python` or `RDOXIDE_REF_PYTHON`. (PyPI only has
2026.03.x, which predates the RDL-based SSSR, so it is not a valid oracle for
ring order.)

Current results: `parse` 10312/10312, `rings` 10312/10312.

## Known deviations (to revisit)

* Atom/Bond objects returned from a Mol hold a reference to the Python Mol plus
  an index (upstream holds a C++ pointer); after removing atoms the index may
  point at a different atom (upstream: dangling pointer).
* `Mol(pickle)` / `ToBinary` not available until MolPickler (stage 6).
* `Atom.IsInRing()` needs stage 3 (`findSSSR`).
* `Chem.SetAllowNontetrahedralChirality` & co. keep the value in-process;
  upstream also writes it to the process environment.
