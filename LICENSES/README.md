# Third-party licenses

| File | Applies to |
|------|-----------|
| `RDKit-BSD-3-Clause.txt` | All code translated from RDKit, verbatim copies of RDKit Python modules, RDKit-derived data files (`crates/rdoxide/data/`) and ported tests. Source: `license.txt` in the RDKit repository. |
| `RingDecomposerLib-BSD-3-Clause.txt` | Code translated from RingDecomposerLib v1.1.3_rdkit (University of Hamburg, ZBH; Andresen, Flachsenberg, Rarey), used by RDKit for ring families / SSSR: `crates/rdoxide/src/external/ring_decomposer_lib/`. |
| `LLVM-Apache-2.0-with-LLVM-exception.txt` | `crates/rdoxide/src/rdgeneral/cxx_sort.rs`, a translation of libc++'s `std::sort` (LLVM Project), used to reproduce the element order of RDKit's unstable sorts. |
| `RDKit-ChemReactions-Novartis-BSD-3-Clause.txt` | Code translated from `Code/GraphMol/ChemReactions` and `Code/GraphMol/Descriptors/USRDescriptor.cpp` (Novartis). Source: `Code/GraphMol/ChemReactions/license.txt` and the USRDescriptor.cpp header (same text). |
| `RDKit-Descriptors3D-Godin-BSD-3-Clause.txt` | Code translated from the descriptor files by Guillaume Godin (Institute of Cancer Research): `AUTOCORR2D/3D`, `MolData3Ddescriptors`, `Data3Ddescriptors` (`WHIM`, `RDF`, `MORSE`, `EEM`, `CoulombMat`, `GETAWAY`): `crates/rdoxide/src/graphmol/descriptors/data_3d.rs`, `autocorr2d.rs`, `autocorr3d.rs`, `rdf.rs`, `morse.rs`, `whim.rs`, `eem.rs`, `coulomb_mat.rs`, `getaway.rs`. Source: the file headers in `Code/GraphMol/Descriptors/`. |
| `Eigen-MPL-2.0.txt` | `crates/rdoxide/src/numerics/eigen_3x3.rs`, `full_piv_lu.rs` and `jacobi_svd.rs`, translations of Eigen 5.0.1 (`SelfAdjointEigenSolver`, `Tridiagonalization`, `JacobiRotation`, `FullPivLU` with its triangular solves and matrix-vector kernel, `JacobiSVD`; Copyright Gael Guennebaud, Benoit Jacob, Jitse Niesen), which RDKit uses for principal axes, EEM charges and WHIM. Mozilla Public License 2.0: the file stays under MPL-2.0 (file-level copyleft) and names its sources. |
| `RDKit-PBF-InstituteOfCancerResearch-BSD-3-Clause.txt` | Code translated from `Code/GraphMol/Descriptors/PBF.cpp` (plane of best fit; Nicholas Firth, Institute of Cancer Research): `crates/rdoxide/src/graphmol/descriptors/pbf.rs`. Source: the file header. |

When porting a component that carries its own license (anything under
upstream `External/`, fonts in `Data/Fonts`, ...), copy its license file here,
add a row to this table and to `NOTICE`.

## Rust crate dependencies

`CRATES.md` lists every crate in `Cargo.lock` with its license. It is generated
and enforced by `tools/check_licenses.py` (run in CI): only the permissive
licenses on the script's allow-list are accepted, so adding a dependency under
any other license needs a deliberate change there. After changing
dependencies run `python tools/check_licenses.py --write` and commit the result.
