# Avalon Toolkit port: notes

Status: stage A0 done (2026-10-28): tarball downloaded to `.work/avalon/` (not committed) and its
MD5 verified; license copied to `LICENSES/Avalon-BSD-3-Clause-Novartis.txt`.

## 1. Source and license

* Source RDKit builds against: `ava-formake` tag `AvalonToolkit_2.0.5-pre.3`,
  <https://github.com/rdkit/ava-formake/archive/refs/tags/AvalonToolkit_2.0.5-pre.3.tar.gz>,
  MD5 `7a20c25a7e79f3344e0f9f49afa03351` (from `External/AvalonTools/CMakeLists.txt`).
  The `README` in RDKit's tree still talks about a patch, but the CMake file of this tag applies
  none: the `rdkit/ava-formake` tag is already RDKit's patched fork (A0 finding).
* License (`license.txt` at the top of the repository, read 2026-10-28): a BSD-style
  3-clause license, **Copyright 2001-2011 Novartis Pharma AG**; source and binary
  redistribution allowed with or without modification if (1) source redistributions keep the
  notice, (2) binary redistributions reproduce notice, conditions and disclaimer in the
  documentation, (3) the names of the copyright holder/contributors are not used to endorse
  derived products; "AS IS" disclaimer. No copyleft, no linking exception needed. The text
  has two sloppy spots (condition 3 names "Novartis Institutes for BioMedical Research Inc.",
  the disclaimer keeps a `<COPYRIGHT HOLDER>` placeholder), so copy it verbatim and do not
  edit it.
* Compatibility: permissive, so it fits the existing mix (BSD-3-Clause AND MIT AND Apache-2.0
  WITH LLVM-exception). When code is ported:
  1. add `LICENSES/Avalon-BSD-3-Clause-Novartis.txt` with the verbatim `license.txt`,
  2. add a row to `LICENSES/README.md` naming the Rust files translated from it,
  3. put the copyright line in the header of every translated Rust file,
  4. note it in the wheel metadata like the InChI MIT license.
  Verify against the downloaded tarball (the license text there is the authority, not this note).

## 2. What RDKit uses (port only this)

Entry points (`AvalonTools.h`, `Wrap/pyAvalonTools.cpp`, Python `rdkit.Avalon.pyAvalonTools`):

| RDKit function | Avalon C functions it calls |
|---|---|
| `GetAvalonFP`, `GetAvalonCountFP`, `GetAvalonFPAsWords` (mol, SMILES or molblock) | `SMIToMOL` or `MolStr2Mol` (reader), `SetFingerprintBits`, `SetFingerprintCountsWithFocus` (`fingerprint.c`) |
| `GetCanonSmiles` (mol, SMILES, molblock) | `MolStr2Mol`, `MOLToSMI`, `CanSmiles` (`canonizer.c`, `smi2mol.c`) |
| `Generate2DCoords` | `SMIToMOL`, `RecolorMolecule`, `LayoutMolecule` (`layout.c`, `didepict.c`, `depictutil.c`) |
| `CheckMolecule`, `CheckMoleculeString`, `InitializeCheckMol`, `GetCheckMolLog`, `CloseCheckMolFiles`, enums `StruChkFlag`, `StruChkResult` | `RunStruchk`, `ClearParameters`, `InitCheckMol`, `GetMsgList`, `FreeMsgList`, `CloseOpenFiles` (`struchk.c` and its pattern/parameter data) |
| plumbing | `MolToMolStr`, `FreeMolecule`, `MyFree` (`reaccsio.c`, `utilities.c`) |

RDKit's CMake compiles these Avalon files (`avalon_clib_srcs`): `layout symboltable patclean
utilities symbol_lists stereo set perceive local graph geometry forio depictutil denormal
casutils ssmatch rtutils smi2mol didepict pattern canonizer` plus (continuing in the CMake
list) the reaccs reader/writer, `fingerprint` and the `struchk` group. The exact list has to be
read from the CMake file once the tarball is available; everything not reachable from the
functions in the table (the `programs/` tools, the Java/JNI part, GUI, converters) is out of scope.

Python surface to reproduce: `rdoxide.Avalon.pyAvalonTools` (function names above),
`rdkit/Avalon/__init__.py` is empty, `Chem.MolKey` (uses `CheckMolecule`/`GetCanonSmiles`;
needs the InChI port, which is done).

Upstream tests to turn green: `External/AvalonTools/Wrap/testAvalonTools.py` (14 tests),
`External/AvalonTools/test1.cpp` (data in `test_data/`), `UnitTestMolKey` (5 of 15 not passing
now, per STATUS row 280).

## 3. Proposed stages (one branch each, differential corpus against official RDKit for each)

| Stage | Content | Check |
|---|---|---|
| A0 | Download tarball into `.work/`, verify MD5, apply RDKit's patch, copy the license, produce the exact function-reachability list from the C sources (call graph from the entry points above) | list committed in this file |
| A1 | `reaccs` data structures, molblock reader/writer (`MolStr2Mol`, `MolToMolStr`), `utilities`, `set`, `graph`, `symboltable`/`symbol_lists`, `local` | molblock round trip vs RDKit's wrapper on the SDF test data |
| A2 | `perceive`, `ssmatch`, `denormal`, `stereo`, `smi2mol` (`SMIToMOL`, `MOLToSMI`) | `GetCanonSmiles` inputs |
| A3 | `canonizer` (`CanSmiles`) | `GetCanonSmiles` vs RDKit on a corpus |
| A4 | `fingerprint` (`SetFingerprintBits`, `...CountsWithFocus`) | `GetAvalonFP` bit vectors and counts, all `nBits`/flag options |
| A5 | `struchk` (checkmol) with its parameter/pattern tables | `CheckMolecule` result tuples and log text |
| A6 | `layout`, `didepict`, `depictutil` (`LayoutMolecule`) | `Generate2DCoords` coordinates vs RDKit (float tolerance only if C uses libm differently) |
| A7 | Python module `rdoxide.Avalon.pyAvalonTools`, `Chem.MolKey` wiring, ratchet update | `testAvalonTools.py`, `UnitTestMolKey` |

Order is chosen so the high-value, self-contained pieces (fingerprints, canonical SMILES,
check) come before the largest and riskiest one (layout). A1 is the prerequisite of all.

## 4. Risks

* The C code relies on global state, `static` buffers and macros (`MyFree`, `TypeAlloc`);
  port with explicit state and plain `Vec`s; compare against a build of the C sources made
  the way `docs/INCHI_PORT.md` describes for InChI (second oracle when the wheel is not enough).
* Memory bugs in the C code (RDKit already patches one): build the oracle with the guard-page
  allocator trick (see the InChI fold investigation) before trusting an odd reference output.
* `struchk` reads parameter files at run time in the original; they must be embedded or
  generated into the crate (as `tools/gen_inchi_eldata.py` does for InChI).
* Locale: RDKit wraps reader calls in `LocaleSwitcher`; Rust parsing must be locale independent.
* Size is unknown until A0 counts the reachable lines; the earlier estimate (50-70k lines in
  the whole toolkit) shrinks to the reachable subset.

## 5. Reachability (A0 result)

`tools/avalon_reach.py` builds the call graph of the 473 function definitions in the 26 C files
RDKit compiles (clang AST, function references included) and walks it from the entry points of
section 2. Result: **324 functions, about 23.4k of 27.3k C lines (86 %)** are reachable; the list
is `docs/avalon_reachable.txt`. So "only what RDKit uses" saves little here: almost the whole
library is used. Biggest reachable parts: `layout.c` 5.0k, `ssmatch.c` 3.8k, `smi2mol.c` 3.0k,
`utilities.c` 1.9k, `struchk.c` 1.3k, `reaccsio.c` 1.2k, `patclean.c` 1.1k, `canonizer.c` 1.0k,
`fixcharges.c` 0.9k, `stereo.c` 0.9k. Not reachable: `didepict.c` (0.95k), `depictutil.c`,
`rtutils.c`, `shortcut.c`, `sketch.c`, `jni`, `programs/*` except `struchk.c`.
The oracle is the official wheel (it contains `rdkit.Avalon.pyAvalonTools`); a C build is
possible as a second oracle from the same tarball.

## Status (updated)

Done and checked against the Avalon C code (oracle driver built from the tarball) or against the
official wheel: A1 (MOL reader/writer), A2a-d (sets, rings, aromaticity, denormalization, `MOLToSMI`),
A4 (all fingerprint classes) and the first part of A7: `rdoxide.Avalon.pyAvalonTools` with
`GetAvalonFP`, `GetAvalonCountFP` and `GetAvalonFPAsWords` for molecules and MOL blocks
(`tests/python/test_avalon_fp.py`, 780 cases vs the official wheel). Findings worth keeping:

* the Python fingerprint is two passes (the second one with Daylight aromaticity, OR-ed in unless
  `isQuery`); the count fingerprint ignores `bitFlags` (always `USE_ALL_FEATURES`) and adds a second
  accumulating pass in query mode;
* `GetAvalonFPAsWords` packs signed `char`s, so bytes >= 0x80 sign-extend;
* still open: SMILES input (`SMIToMOL` needs layout), `GetCanonSmiles`, `CheckMolecule`,
  `Generate2DCoords`, `PerceiveMarkush`, `MolKey` wiring.

### A5 (structure check), progress

* A5a `stereo.rs`: `AtomParity`, `CheckStereo`, `NoParityDefined`, `DubiousStereochemistry`,
  `RemoveDubiousStereochemistry`, `FixDubious3DMolecule`, `AtomClash`, `CisTransPerception`,
  `PermutationParity` (432 cases identical to the C code, `tests/data/avalon_stereo.jsonl`). Messages
  go to a `Vec<String>` instead of the global list. Not ported: `IsStereoMatch` (needs `ssmatch.c`)
  and `AllCentersRefined` (reads one element past its sort buffer in C).
