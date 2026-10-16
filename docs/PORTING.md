# Porting guide

Conventions for translating RDKit C++ into rdoxide.

## Where code goes

| Upstream | rdoxide |
|----------|---------|
| `Code/RDGeneral/X.{h,cpp}` | `crates/rdoxide/src/rdgeneral/x.rs` |
| `Code/GraphMol/X.{h,cpp}` | `crates/rdoxide/src/graphmol/x.rs` |
| `Code/GraphMol/Sub/X.cpp` | `crates/rdoxide/src/graphmol/sub/x.rs` |
| `Code/**/Wrap/rdFoo.cpp` | `crates/rdoxide-python/src/<module>.rs` |
| `rdkit/**/*.py` | `python/rdoxide/**` (verbatim via `porting/copies.toml` when possible) |
| data embedded in C++ | `crates/rdoxide/data/` extracted by `tools/sync_upstream.py regen` |

Every new Rust module gets an entry in `porting/manifest.toml` listing its
upstream sources; run `python3 tools/sync_upstream.py stamp` afterwards.

## File header

```rust
//! Ported from RDKit `Code/GraphMol/Foo.cpp`.
//!
//! Copyright (C) <years> <upstream copyright holders>.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
```

## Fidelity rules

* Behaviour, including numeric output, string formatting, atom/bond
  ordering and error messages, must match upstream. When upstream iterates a
  `std::map`/`std::set`, use `BTreeMap`/`BTreeSet`; never iterate a `HashMap`
  where order is observable.
* Keep upstream function names (snake_cased) and structure so diffs between
  RDKit releases can be applied mechanically.
* `PRECONDITION`/`CHECK_INVARIANT` → `Err(RDKitError::precondition(..))` /
  `RDKitError::invariant(..)`; these surface as `RuntimeError` in Python.
* Random numbers: use `rdgeneral::random` generators that reproduce upstream's.
* No `unsafe`; no C/C++ dependencies.

## Tests

* Translate upstream `catch_*.cpp`/`test*.cpp` cases into `#[cfg(test)]`
  modules (or `crates/rdoxide/tests/`), with a comment naming the source test.
* Python behaviour is tested by (1) upstream Python tests run verbatim through
  `tools/run_upstream_tests.py`, and (2) focused pytest files in
  `tests/python/` ported from upstream tests.
* After a change: `cargo test --workspace`, `maturin develop`,
  `pytest tests/python`, `python tools/run_upstream_tests.py --update`.
