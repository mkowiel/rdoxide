# Third-party licenses

| File | Applies to |
|------|-----------|
| `RDKit-BSD-3-Clause.txt` | All code translated from RDKit, verbatim copies of RDKit Python modules, RDKit-derived data files (`crates/rdoxide/data/`) and ported tests. Source: `license.txt` in the RDKit repository. |
| `RingDecomposerLib-BSD-3-Clause.txt` | Code translated from RingDecomposerLib v1.1.3_rdkit (University of Hamburg, ZBH; Andresen, Flachsenberg, Rarey), used by RDKit for ring families / SSSR: `crates/rdoxide/src/external/ring_decomposer_lib/`. |
| `LLVM-Apache-2.0-with-LLVM-exception.txt` | `crates/rdoxide/src/rdgeneral/cxx_sort.rs`, a translation of libc++'s `std::sort` (LLVM Project), used to reproduce the element order of RDKit's unstable sorts. |
| `RDKit-ChemReactions-Novartis-BSD-3-Clause.txt` | Code translated from `Code/GraphMol/ChemReactions`, `Code/GraphMol/Descriptors/USRDescriptor.cpp` and RDKit's InChI wrapper `External/INCHI-API` (`inchi.cpp`, `Wrap/pyInchi.cpp`, `python/inchi.py`; Novartis): `crates/rdoxide/src/external/inchi_api.rs`, `crates/rdoxide-python/src/rdinchi.rs`, `python/rdoxide/Chem/inchi.py`. Source: `Code/GraphMol/ChemReactions/license.txt` and the USRDescriptor.cpp / inchi.cpp headers (same text). |
| `InChI-MIT.txt` | Code translated from the InChI Software library v1.07.3 (IUPAC and InChI Trust; the version RDKit downloads in `External/INCHI-API`): `crates/rdoxide/src/external/inchi/` (except `sha256.rs`, which is original rdoxide code). MIT. Source: the C file headers of the v1.07.3 release `INCHI-1-SRC.zip` and the repository `LICENSE` at tag v1.07.3 (see `docs/INCHI_PORT.md`). |
| `Avalon-BSD-3-Clause-Novartis.txt` | Code translated from the Avalon Toolkit 2.0.5-pre.3 (Novartis Pharma AG, 2001-2011; the version RDKit downloads in `External/AvalonTools`): `crates/rdoxide/src/external/avalon/` (except `scan.rs`, which only emulates C `sscanf`/`printf`), including the rule files `data/checkfgs.chk` and `data/checkfgs.trn` copied from the toolkit's `struchk` resources (also shipped as `python/rdoxide/share/RDKit/Data/struchk/`, RDKit's `Data/struchk`). BSD 3-clause. See `docs/AVALON_PORT.md`. |
| `RDKit-Descriptors3D-Godin-BSD-3-Clause.txt` | Code translated from the descriptor files by Guillaume Godin (Institute of Cancer Research): `AUTOCORR2D/3D`, `MolData3Ddescriptors`, `Data3Ddescriptors` (`WHIM`, `RDF`, `MORSE`, `EEM`, `CoulombMat`, `GETAWAY`): `crates/rdoxide/src/graphmol/descriptors/data_3d.rs`, `autocorr2d.rs`, `autocorr3d.rs`, `rdf.rs`, `morse.rs`, `whim.rs`, `eem.rs`, `coulomb_mat.rs`, `getaway.rs`. Source: the file headers in `Code/GraphMol/Descriptors/`. |
| `Eigen-MPL-2.0.txt` | `crates/rdoxide/src/numerics/eigen_3x3.rs`, `full_piv_lu.rs` and `jacobi_svd.rs`, translations of Eigen 5.0.1 (`SelfAdjointEigenSolver`, `Tridiagonalization`, `JacobiRotation`, `FullPivLU` with its triangular solves and matrix-vector kernel, `JacobiSVD`; Copyright Gael Guennebaud, Benoit Jacob, Jitse Niesen), which RDKit uses for principal axes, EEM charges and WHIM. Mozilla Public License 2.0: the file stays under MPL-2.0 (file-level copyleft) and names its sources. |
| `RDKit-PBF-InstituteOfCancerResearch-BSD-3-Clause.txt` | Code translated from `Code/GraphMol/Descriptors/PBF.cpp` (plane of best fit; Nicholas Firth, Institute of Cancer Research): `crates/rdoxide/src/graphmol/descriptors/pbf.rs`. Source: the file header. |
| `Boost-BSL-1.0.txt` | `crates/rdoxide/src/rdgeneral/json.rs` (`number` module): the digit bookkeeping and `dec_to_float` of Boost.JSON 1.85's `basic_parser_impl.hpp` (Copyright (c) 2019 Vinnie Falco, (c) 2020 Krystian Stasiowski), reproduced so that doubles read from MolJSON are bit-identical to what RDKit reads. Boost Software License 1.0. |
| `SciPy-BSD-3-Clause.txt` | `crates/rdoxide/src/graphmol/rascal_mces/lap.rs`, a translation of `lap_a_la_scipy.cpp` (Crouse's shortest augmenting path rectangular assignment algorithm as implemented by PM Larsen in SciPy; Copyright (c) 2001-2002 Enthought, Inc., 2003-2023 SciPy Developers), used by the RASCAL tier 2 similarity. BSD 3-Clause. |
| `Roboto-Apache-2.0.txt` | `crates/rdoxide/data/fonts/Roboto-Regular.ttf`, the Roboto Regular font (developed by Google) extracted from RDKit's `roboto_regular.h`; used to lay out and outline text in the molecule drawing code (`graphmol/mol_draw2d`). Apache License 2.0. |
| `Telex-OFL-1.1.txt` | `crates/rdoxide/data/fonts/Telex-Regular.ttf`, the Telex Regular font (Copyright 2011 The Telex Project Authors) shipped with RDKit as its default drawing font. SIL Open Font License 1.1. |

When porting a component that carries its own license (anything under
upstream `External/`, fonts in `Data/Fonts`, ...), copy its license file here,
add a row to this table and to `NOTICE`.

## Rust crate dependencies

`CRATES.md` lists every crate in `Cargo.lock` with its license. It is generated
and enforced by `tools/check_licenses.py` (run in CI): only the permissive
licenses on the script's allow-list are accepted, so adding a dependency under
any other license needs a deliberate change there. After changing
dependencies run `python tools/check_licenses.py --write` and commit the result.
