# Third-party licenses

| File | Applies to |
|------|-----------|
| `RDKit-BSD-3-Clause.txt` | All code translated from RDKit, verbatim copies of RDKit Python modules, RDKit-derived data files (`crates/rdoxide/data/`) and ported tests. Source: `license.txt` in the RDKit repository. |
| `RingDecomposerLib-BSD-3-Clause.txt` | Code translated from RingDecomposerLib v1.1.3_rdkit (University of Hamburg, ZBH; Andresen, Flachsenberg, Rarey), used by RDKit for ring families / SSSR: `crates/rdoxide/src/external/ring_decomposer_lib/`. |
| `LLVM-Apache-2.0-with-LLVM-exception.txt` | `crates/rdoxide/src/rdgeneral/cxx_sort.rs`, a translation of libc++'s `std::sort` (LLVM Project), used to reproduce the element order of RDKit's unstable sorts. |
| `RDKit-ChemReactions-Novartis-BSD-3-Clause.txt` | Code translated from `Code/GraphMol/ChemReactions` (Novartis). Source: `Code/GraphMol/ChemReactions/license.txt`. |

When porting a component that carries its own license (anything under
upstream `External/`, fonts in `Data/Fonts`, ...), copy its license file here,
add a row to this table and to `NOTICE`.

## Rust crate dependencies

`CRATES.md` lists every crate in `Cargo.lock` with its license. It is generated
and enforced by `tools/check_licenses.py` (run in CI): only the permissive
licenses on the script's allow-list are accepted, so adding a dependency under
any other license needs a deliberate change there. After changing
dependencies run `python tools/check_licenses.py --write` and commit the result.
