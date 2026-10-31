# Changelog

All notable changes of rdoxide. Versions follow the ported RDKit release
(`2026.9.1` = RDKit 2026.09.1) with a pre-release suffix while the port is
incomplete (`a0` = first alpha, PEP 440; `2026.9.1-alpha.0` in Cargo).

## 2026.9.1a0 — first alpha

First public pre-release. A pure-Rust port of RDKit 2026.09.1 with the same
Python API under the package name `rdoxide` (Python 3.10 – 3.15, one `abi3`
wheel per platform, Linux/macOS/Windows).

### Ported

* **Core and formats**: `rdchem`, sanitization, aromaticity, stereo perception
  (legacy, new, CIP labeler), canonical ranking, `rdmolops`; SMILES/SMARTS/
  CXSMILES (including substance groups, link nodes and variable attachment
  points), mol blocks V2000/V3000, SD/SMILES/TDT suppliers and writers, PDB,
  XYZ, Mol2, TPL, CML, MolJSON, Marvin, CDXML, SCSR, binary pickles compatible
  with RDKit's.
* **Search and analysis**: substructure matching, `SubstructLibrary`, MolBundle,
  tautomer and generalized queries, FMCS, RASCAL, MMPA, SynthonSpace,
  R-group decomposition (all matching strategies, including the genetic
  algorithm), descriptors (all 217 of `Descriptors._descList`), fingerprints
  (vectorised popcount kernels), pickers and clustering, filter catalogs.
* **Chemistry transforms**: reactions and enumeration, standardization,
  tautomers, scaffolds, hashing, resonance, stereoisomer enumeration.
* **3D**: depiction (templates, constrained, CoordGen), distance geometry/ETKDG,
  UFF/MMFF, alignment, shape tools, `DetermineBonds` (including extended
  Hückel connectivity), trajectories.
* **Drawing**: `rdMolDraw2D` with SVG and PNG (tiny-skia) output, `Chem.Draw`.
* **External libraries translated to Rust**: InChI 1.07.3 (including polymers),
  Avalon Toolkit, CoordGen, FreeSASA, YAeHMOP (`rdEHTTools`), pubchem-align3d
  (`rdShapeAlign`), RingDecomposerLib.

### Verification

* 2081 of 2105 tests of RDKit's own Python test-suite pass (the ratchet in
  `tests/upstream/expected_pass.txt`), alongside Rust and pytest suites and
  differential comparisons against the official RDKit wheel.
* Multicore: batch APIs use all cores with results independent of the thread
  count (`RDOXIDE_NUM_THREADS`).

### Known gaps

* Not ported: binary ChemDraw CDX (proprietary SDK), Qt drawing
  (`rdMolDraw2DQt`), PyMOL helpers, Boost-archive serialisation, InChI CRU
  folding (the official library crashes on it), the InChI time limit (`-W`).
* Differences from the official wheel: 12 of the 2105 upstream tests fail (with the long tests enabled by
  `--long`: 2087 pass) and 6 are skipped by upstream itself (long tests, missing test database); each
  one is listed with its reason in [docs/UPSTREAM_TESTS.md](docs/UPSTREAM_TESTS.md).
* The API is not frozen; symbols that are not ported are absent rather than
  stubbed.
