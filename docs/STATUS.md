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
| 2 | Atom, Bond, Conformer, ROMol/RWMol, RingInfo | ⬜ | next |
| 3 | Perception (rings, valence, aromaticity, kekulize, sanitize, stereo) | ⬜ | |
| 4 | SMILES/SMARTS | ⬜ | |
| 5–13 | see PLAN.md | ⬜ | |

## Python package files

| File | Kind |
|------|------|
| `rdkit/__init__.py`, `RDLogger.py`, `RDRandom.py` | verbatim upstream (porting/copies.toml) |
| `rdkit/rdBase.py`, `rdkit/Chem/rdchem.py` | shims → `rdkit._rdoxide` |
| `rdkit/Chem/__init__.py` | provisional (until rdmolfiles/rdmolops/... exist) |

## Iteration log

* 2026-10-07 — Stage 0 + 1 + PeriodicTable.
