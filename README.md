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

> **Status: pre-alpha.** The port proceeds in stages; see
> [docs/PLAN.md](docs/PLAN.md) and [docs/STATUS.md](docs/STATUS.md) for what
> is available. Symbols that are not ported yet are absent rather than stubbed.

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
derived from RDKit remain under RDKit's BSD 3-Clause license. See [NOTICE](NOTICE),
[LICENSE](LICENSE) and [LICENSES/](LICENSES/).
