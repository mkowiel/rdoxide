# rdoxide

**RDKit, rusted.** A pure-Rust port of the [RDKit](https://www.rdkit.org)
cheminformatics toolkit with the *same Python API* under its own package
name, so code written for RDKit only needs its imports changed, and both
libraries can be installed side by side:

```python
from rdoxide import Chem          # instead of: from rdkit import Chem
```

rdoxide is an independent project and is not affiliated with or endorsed by
the RDKit project.

* Ported from **RDKit 2026.09.1** (see [`UPSTREAM.toml`](UPSTREAM.toml)).
* Python **3.10 – 3.15** (single `abi3` wheel).
* No C/C++ dependencies — the whole stack below Python is Rust.

> **Status: alpha (2026.9.1a0).** Most of RDKit's Python API is ported and
> checked against the official 2026.09.1 build, but expect gaps and rough
> edges. 2076 of 2101 tests of RDKit's own Python test-suite pass (see
> [docs/UPSTREAM_TESTS.md](docs/UPSTREAM_TESTS.md) for the rest and [CHANGELOG.md](CHANGELOG.md)). Symbols
> that are not ported are absent rather than stubbed.

## What is ported

Ported, with the same Python module and function names as RDKit (per-stage
detail and test evidence in [docs/STATUS.md](docs/STATUS.md)):

* **Core**: `rdchem` (atoms, bonds, molecules, conformers, rings, stereo and
  substance groups), sanitization, aromaticity, stereo perception (legacy and
  new, CIP labeler), canonical ranking, `rdmolops`, `rdqueries`.
* **Formats**: SMILES/SMARTS/CXSMILES, mol blocks (V2000/V3000), SD/SMILES/TDT
  suppliers and writers (incl. multithreaded), PDB, XYZ, Mol2, TPL, CML, MolJSON,
  binary pickles compatible with RDKit's.
* **Search**: substructure matching, `SubstructLibrary`, MolBundle, tautomer and
  generalized queries, FMCS, RASCAL, R-group decomposition, MMPA, SynthonSpace.
* **Descriptors and fingerprints**: `rdMolDescriptors`, `Descriptors`, Morgan,
  RDKit, atom-pair, torsion, MACCS, pattern, MHFP, `DataStructs`, pickers,
  clustering, filter catalogs, feature factories.
* **Chemistry transforms**: reactions (incl. enumeration), standardization,
  tautomers, scaffolds, hashing, resonance, stereoisomer enumeration.
* **3D**: depiction, distance geometry/ETKDG, UFF/MMFF, alignment, shape
  and Gaussian-shape tools, interaction fields, trajectories.

* **External libraries** (translated to Rust, each under its own license):
  InChI 1.07.3 (including polymers), Avalon Toolkit, CoordGen, FreeSASA,
  YAeHMOP (`rdEHTTools`), pubchem-align3d (`rdShapeAlign`), RingDecomposerLib.
* **Drawing**: `rdMolDraw2D` with SVG and PNG (tiny-skia) output.

Not ported yet: binary ChemDraw CDX (needs the proprietary ChemDraw SDK),
PyMOL helpers, Boost-archive serialisation, the InChI time limit, and the
remaining gaps listed as partial in [docs/STATUS.md](docs/STATUS.md).

Not planned: the C++ API and headers, the Java/C# wrappers, the PostgreSQL
cartridge, the Cairo/Qt/FreeType drawing back ends (drawing is SVG and PNG),
`Contrib/` and the demos/GUI. Symbols that are not ported are absent rather
than stubbed, so `hasattr` checks behave like an RDKit build without them.

## Installing

Wheels of each release (one `abi3` wheel per platform: Linux x86-64/aarch64,
macOS x86-64/arm64, Windows x86-64) are attached to the
[GitHub releases](https://github.com/mkowiel/rdoxide/releases); no Rust
toolchain is needed:

```bash
pip install --pre rdoxide-2026.9.1a0-*.whl
```

### From the repository

To get the latest development version, install straight from the repository. This
builds from source and needs a [Rust toolchain](https://rustup.rs) (no C/C++
compiler); `maturin` is fetched automatically as the build backend.

```bash
# uv (into the active environment)
uv pip install "git+https://github.com/mkowiel/rdoxide.git"

# uv (add to a project)
uv add "rdoxide @ git+https://github.com/mkowiel/rdoxide.git"

# pip
pip install "git+https://github.com/mkowiel/rdoxide.git"

# a specific branch, tag or commit
pip install "git+https://github.com/mkowiel/rdoxide.git@main"
```

Every CI run also uploads the built `abi3` wheel as an artifact (Actions tab),
which installs without a Rust toolchain: `pip install rdoxide-*.whl`.

```python
from rdoxide import Chem
Chem.MolToSmiles(Chem.MolFromSmiles("c1ccccc1O"))
```

## Threads

rdoxide uses all cores by default wherever the work can be split without
changing the result (conformer embedding, force-field optimisation of
conformers, bulk fingerprints and similarities, substructure libraries,
reaction runs, topological distance matrices, synthon search, ...). Parallel
results are identical to single-threaded ones.

The number of threads any operation may use is capped by the environment
variable `RDOXIDE_NUM_THREADS`:

```bash
RDOXIDE_NUM_THREADS=4 python my_script.py     # at most 4 threads
RDOXIDE_NUM_THREADS=-1 python my_script.py    # all cores but one
RDOXIDE_NUM_THREADS=1 python my_script.py     # no parallelism
```

A positive number is that many threads, `0` or `all` (the default) all
hardware threads, and a negative number all but that many. At run time use
`rdBase.SetMaxThreads(n)` and `rdBase.GetMaxThreads()`. Functions and
parameter objects that take `numThreads` keep their meaning: a positive value
asks for that many threads (never more than the cap), `0` for the cap and a
negative value for the cap minus that many. `EmbedParameters.numThreads` still
reads `1` as upstream does, but embedding uses the cap until you set it.

## Building

```bash
uv venv --python 3.13 .venv && source .venv/bin/activate
uv pip install maturin pytest numpy
maturin develop --release
pytest tests/python
cargo test --workspace
```

## Testing against RDKit's own test-suite

```bash
python3 tools/sync_upstream.py fetch          # clone upstream RDKit at the ported tag
python tools/run_upstream_tests.py -v         # run upstream Python tests against rdoxide
```

`tests/upstream/expected_pass.txt` is a ratchet: tests listed there must keep
passing; `--update` records newly passing ones.

## Updating to a new RDKit release

See "Updating to a new RDKit release" in [docs/PLAN.md](docs/PLAN.md):
bump `UPSTREAM.toml`, `sync_upstream.py fetch && regen && diff`, port the
listed changes, `stamp`.

## License

rdoxide is released under the BSD 3-Clause License; code, data and Python modules
derived from RDKit remain under RDKit's BSD 3-Clause license. Ported third-party
libraries keep their own licenses, e.g. the InChI support translates the InChI
Software 1.07.3 (Copyright IUPAC and InChI Trust, MIT license). See [NOTICE](NOTICE),
[LICENSE](LICENSE) and [LICENSES/](LICENSES/).
