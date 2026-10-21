# rdoxide — porting plan

rdoxide ("rusted RDKit") is a pure-Rust re-implementation of the
[RDKit](https://github.com/rdkit/rdkit) cheminformatics toolkit that ships the
**same Python API** (`import rdoxide`, `from rdoxide import Chem`, ...) so that it
can replace RDKit 1:1 in existing Python code.

* Upstream reference version: **RDKit 2026.09.1** (`Release_2026_09_1`,
  commit `fece8caa860bdf6c9c82bdb9253f22a0cf87150c`) — see `UPSTREAM.toml`.
* Supported Python: **CPython 3.10 – 3.15** (the wheel matrix RDKit publishes
  for its 2026 releases; built as a single `abi3` wheel, `abi3-py310`).
* No C/C++ dependencies: everything below the Python boundary is Rust, so there
  is no FFI overhead between components (SMILES parser → sanitizer →
  substructure matcher → fingerprints, ...).

## Architecture

```
crates/
  rdoxide/            core library (pure Rust, no Python). Module tree mirrors
                      upstream Code/ (rdgeneral, graphmol, graphmol::smiles_parse,
                      graphmol::substruct, data_structs, ...). Each ported file
                      carries an "Ported from RDKit <path>" header.
  rdoxide-python/     PyO3 bindings, mirrors upstream */Wrap/*.cpp. Built as
                      the native extension `rdoxide._rdoxide` (abi3).
python/rdoxide/        The importable `rdoxide` package. Pure-Python modules are
                      copied verbatim from upstream (BSD-3) by tools/sync_upstream.py;
                      native modules (rdoxide.Chem.rdchem, rdoxide.rdBase, ...) are
                      thin shims that re-export from rdoxide._rdoxide.
data/                 data tables extracted from upstream (periodic table,
                      isotopes, ...) by tools/, embedded with include_str!.
tests/python/         Python tests: upstream tests ported to pytest + our own.
porting/manifest.toml Map: rust module <-> upstream files (+ git blob hash at the
                      ported version). Drives upgrade diffs.
LICENSES/             Third-party licenses (RDKit BSD-3-Clause, Novartis BSD, ...).
```

### Updating to a new RDKit release

1. Bump `UPSTREAM.toml` (`tag`, `commit`).
2. `python3 tools/sync_upstream.py fetch` — clones upstream into `upstream/` (git-ignored).
3. `python3 tools/sync_upstream.py regen` — re-extracts data tables and re-copies the
   pure-Python package files / test data.
4. `python3 tools/sync_upstream.py diff` — lists, per rust module, the upstream
   files that changed since the last port (blob hash mismatch in
   `porting/manifest.toml`). Port those diffs, then `sync_upstream.py stamp` to
   record the new hashes.
5. Run the full test suite (`cargo test`, `pytest tests/python`).

## Stages

Each stage ends with: Rust unit tests ported from upstream `catch_*.cpp` /
`test*.cpp`, Python tests ported from upstream `rough_test.py` /
`UnitTest*.py` for the covered API, and green CI. Upstream behaviour is the
oracle: expected values are taken from upstream tests, and (where available)
a differential test compares against the official `rdkit` wheel.

| # | Stage | Upstream scope | Python surface |
|---|-------|----------------|----------------|
| 0 | Infrastructure | — | packaging, CI, licenses, sync tooling, plan |
| 1 | RDGeneral | RDGeneral (Dict/RDValue, exceptions, logging, versions) | `rdoxide.rdBase`, `RDLogger` |
| 2 | Core graph | PeriodicTable, Atom, Bond, Conformer, ROMol/RWMol, RingInfo, StereoGroup, SubstanceGroup, props | `Chem.rdchem` (Atom, Bond, Mol, RWMol, EditableMol, Conformer, PeriodicTable) |
| 3 | Perception | FindRings (SSSR/symmetrized), valence, Kekulize, Aromaticity, ConjugHybrid, sanitize, AddHs/RemoveHs, Chirality/FindStereo, CIP (legacy + new CIPLabeler) | `Chem.rdmolops` core |
| 4 | SMILES / SMARTS | SmilesParse (parser + writer), Canon/new_canon, QueryOps/QueryAtom/QueryBond, CXSMILES | `MolFromSmiles`, `MolToSmiles`, `MolFromSmarts`, `MolToSmarts`, `CanonSmiles` |
| 5 | Substructure | Substruct (VF2), SubstructUtils, recursive SMARTS, chirality-aware match | `HasSubstructMatch`, `GetSubstructMatch(es)` |
| 6 | File formats | MolFileParser/Writer (V2000/V3000), SDMolSupplier/SDWriter, SmilesMolSupplier, PDB, XYZ, Mol2, MolPickler (binary-compatible pickles), MolInterchange (JSON) | `Chem.rdmolfiles`, pickling |
| 7 | DataStructs | ExplicitBitVect, SparseIntVect, similarity metrics, BitOps, fingerprint pickles | `rdoxide.DataStructs` |
| 8 | Descriptors & fingerprints | Descriptors (Crippen, TPSA, Lipinski, MolWt, ...), Fingerprints (Morgan, RDKit, AtomPair, Torsion, MACCS, Pattern, Layered), FingerprintGenerator | `rdMolDescriptors`, `rdFingerprintGenerator`, `Descriptors`, `AllChem` |
| 9 | Transforms | ChemTransforms, ChemReactions, MolStandardize, Fragments, Scaffolds, MolHash, Resonance, EnumerateStereoisomers | `rdChemReactions`, `rdMolStandardize`, ... |
| 10 | Geometry & 3D | Geometry, Numerics, MolTransforms, Depictor (2D coords), DistGeom + ETKDG, ForceField (UFF, MMFF), MolAlign, ShapeHelpers | `rdDepictor`, `rdDistGeom`, `rdForceFieldHelpers`, `rdMolAlign` |
| 11 | Search & analysis | FMCS, RascalMCES, RGroupDecomposition, SimDivPickers, ML (clustering), FilterCatalog | corresponding modules |
| 12 | Drawing | MolDraw2D (SVG backend in pure Rust; Cairo-free PNG via a Rust rasterizer) | `rdMolDraw2D`, `Draw` |
| 13 | Long tail | InChI (pure-Rust port of InChI 1.07.3, MIT; see INCHI_PORT.md), Avalon, CDXML, Marvin, SLN, MolEnumerator, MMPA, Abbreviations, ... | — |

Priorities inside each stage are driven by real-world usage: the API used by
`rdoxide.Chem.__init__`, `AllChem`, `Descriptors` and `PandasTools` comes first.

### Python API compatibility rules

* Module, class, function, keyword-argument names and default values match
  upstream's Boost.Python signatures (`rdkit-stubs` is used as the reference).
* Exceptions match (`ValueError`, `RuntimeError`, `KeyError`, rdkit's
  `AtomValenceException` etc. live in `rdoxide.Chem.rdchem`).
* `__module__` of every class equals the upstream module (pickling works).
* Pickles are byte-compatible with RDKit's `MolPickler` format.
* Where a feature is not ported yet, the symbol is absent (never a silent stub),
  so `hasattr` checks behave like a build of RDKit without that feature.

### Engineering practices

* `cargo fmt`, `cargo clippy -D warnings`, `ruff` for Python; CI on Linux/macOS/Windows.
* No `unsafe` in the core crate (`#![forbid(unsafe_code)]`).
* Every ported module lists its upstream sources in `porting/manifest.toml`.
* Benchmarks (`criterion` + Python `pytest-benchmark`) against the official
  wheel for hot paths (SMILES parse/write, sanitize, substruct, fingerprints).

Progress is tracked in [STATUS.md](STATUS.md).

## Follow-ups after the InChI port (I12)

Requested by the project owner; start only once the InChI port is finished.
Upstream ratchet run of 2026-10-21 (1884/1987 passing) shows these failures:

* **Tuple vs list return types** (most of them): upstream `nbWrap/rough_test.py`
  expects Python `list`s where rdoxide returns `tuple`s or the other way round
  (`testBitVectQuery`, `testCrambin`, `testDetectChemistryProblems`,
  `testGithub4138`, `testReplaceAtomWithQueryAtom`, `testSubstructMatchAcetate`,
  `testSubstructParameters*`, `testSubstanceGroupProperties`); audit every
  wrapper's return type against upstream (nanobind returns `list` for
  `std::vector`, `tuple` for `std::tuple`/`std::pair`/explicit `make_tuple`).
* Other `rough_test.py` failures: `test12Smarts`, `test57AddRecursiveQuery`,
  `test16Props`, `testGithub4570` (`HasProp(None)`), `testGithub8877`,
  `testSetEnhancedStereoGroupOwnershipCheck`, `testSetEnhancedStereoTypeChecking`,
  `test_conformer_lifetime` (`PyBorrowError` panic).
* Missing modules/methods: `rdMolDraw2DQt`, `Chem.PandasTools`,
  `BitEnsemble.InitScoreTable`, `testMaxMin.py` (`499500.0` passed where a
  sequence of ints is expected), `testGetPropDefault.py`, `testThreads.py`
  (collected without tests).

