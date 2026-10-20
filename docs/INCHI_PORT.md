# InChI port

RDKit exposes InChI through `rdkit.Chem.rdinchi` (native) and `rdkit.Chem.inchi`
(Python, star-imported into `rdkit.Chem`): `MolToInchi`, `MolToInchiAndAuxInfo`,
`MolBlockToInchi(AndAuxInfo)`, `MolToInchiKey`, `InchiToInchiKey`, `MolFromInchi`,
`MolFromInchiAndAuxInfo`, `GetInchiVersion`, `InchiReadWriteError`,
`INCHI_AVAILABLE`, plus free-form option strings (`"/FixedH /RecMet"`, ...) that
are passed through to the InChI library. This document covers the license audit
and the plan for a faithful pure-Rust port of the InChI library that RDKit
2026.09.1 builds against.

## 1. License

### What RDKit 2026.09.1 builds against

* `Code/cmake/Modules/FindInchi.cmake` (RDKit tag `Release_2026_09_1`, commit
  `fece8caa860bdf6c9c82bdb9253f22a0cf87150c`) downloads, when no system InChI
  is found,
  `https://github.com/IUPAC-InChI/InChI/releases/download/v1.07.3/INCHI-1-SRC.zip`
  with MD5 `b761f37a7d3a3ff8574daf6906993759`.
  (`External/INCHI-API/download-inchi.sh` and `README` still mention 1.04/1.05;
  they are stale and not used by the CMake build.)
* `External/INCHI-API/CMakeLists.txt` compiles 45 files of
  `INCHI_BASE/src` and 12 of `INCHI_API/libinchi/src` (incl. `ixa/`) into
  `libInchi`, then RDKit's own wrapper `inchi.cpp` into `RDInchiLib`.
* Confirmed against the official wheel (conda-forge `rdkit=2026.09.1`):
  `rdinchi.GetInchiVersion() == "1.07.3"`.
* InChI tag `v1.07.3` = commit `0b3e941d29289f3e5024c9ecfd45186285319420`
  in <https://github.com/IUPAC-InChI/InChI>. I downloaded the release zip,
  checked the MD5 (matches the pin), and read the license text in it.

### License of InChI 1.07.3 (primary sources)

* **Every C source and header that RDKit compiles carries the MIT license**,
  checked file by file in the release zip (all `INCHI_BASE/src/*.{c,h}`,
  `INCHI_API/libinchi/src/**/*.{c,h}`), header:
  "International Chemical Identifier (InChI) / Version 1 / Software version
  1.07 / April 30, 2024 / MIT License / Copyright (c) 2024 IUPAC and InChI
  Trust", followed by the standard MIT permission text and the note that the
  library was developed under IUPAC, originally at NIST, with additions by
  IUPAC, the InChI Trust and external contributors.
* The repository `LICENSE` at tag v1.07.3
  (<https://github.com/IUPAC-InChI/InChI/blob/v1.07.3/LICENSE>) is the same
  MIT text with "Copyright (c) 2024 InChI Project"; the README carries an
  MIT badge.
* **Exception, `INCHI_BASE/src/sha2.{c,h}`** (SHA-256; Copyright (C) 2010 Paul
  Bakker, originally Christophe Devine): "free software ... as a part of the
  International Chemical Identifier (InChI) software, under the terms of the
  MIT license ... Note that this licensing is only valid in combination with
  the InChI software package, for all other cases contact Paul Bakker."
  rdoxide does not translate this file. Its SHA-256 is original code written
  from FIPS 180-4 (`external/inchi/sha256.rs`, BSD-3-Clause like the rest of
  rdoxide's own code). Hashes are fully specified, so InChIKeys are identical.
* `INCHI_API/tbb/LICENSE` (Apache-2.0) belongs to an optional Windows
  threading helper that RDKit does not compile; nothing is ported from it.
* History, for context: up to v1.06 (Dec 2020) the sources were under the
  "IUPAC/InChI-Trust Licence No.1.0" (an LGPL-2.1-derived copyleft licence;
  see the headers at tag `v1.06`, e.g. `INCHI-1-SRC/INCHI_BASE/src/ichican2.c`).
  v1.07 (2024) relicensed to MIT. **The version RDKit 2026.09.1 uses, 1.07.3,
  is MIT.** Later 1.07.4/1.07.5 tags exist but are not what RDKit pins, so
  they are not the port target (their InChIs can differ).
* The file `External-contributors` that the headers mention is not in the
  release zip. It names contributors, not license terms (see open questions).

### RDKit wrapper

`External/INCHI-API/inchi.{h,cpp}`, `Wrap/pyInchi.cpp` and `python/inchi.py`:
BSD 3-Clause, "Copyright (c) 2011-2025 Novartis Institutes for BioMedical
Research Inc. and other RDkit contributors". This is the same Novartis text as
`LICENSES/RDKit-ChemReactions-Novartis-BSD-3-Clause.txt`, which now also lists
these files. It adds no new obligation beyond keeping the notice.

### Project policy and verdict

* Project policy (`tools/check_licenses.py`, `LICENSES/README.md`, `NOTICE`,
  `docs/PLAN.md`): crate dependencies must be on a permissive allow-list
  (MIT is on it). Translated third-party code keeps its license, which is
  recorded in `LICENSES/` + `NOTICE`, with source named in each file header.
  No C/C++ dependencies and no `unsafe`. So **linking libinchi is not an
  option; a Rust translation is required.**
* **Verdict: compatible.** A faithful Rust translation of the InChI 1.07.3 C
  sources is a derivative work under MIT, which allows modification and
  redistribution under any terms as long as the notice is kept. MIT is
  already part of rdoxide's combined license expression (linchemin, xyz2mol).
  No new crate is needed, so `check_licenses.py` is unaffected.
* **Obligations** (all done in this stage):
  1. `LICENSES/InChI-MIT.txt`: the full MIT text with both copyright lines,
     the IUPAC/NIST provenance note, the source URL/MD5/tag, and the
     `sha2.c` exception.
  2. Row in `LICENSES/README.md` and entry in the `NOTICE` license index
     (MIT now also "files derived from ... InChI").
  3. Header in every translated file: source file(s) + "Copyright (c) 2024
     IUPAC and InChI Trust. Licensed under the MIT license (see
     LICENSES/InChI-MIT.txt)."
  4. Do not translate `sha2.c`/`sha2.h`.
  5. The `UPSTREAM.toml` `[external.InChI]` pin (url, md5, license), so
     `tools/sync_upstream.py fetch` puts the exact sources in
     `upstream/external/InChI`. `fetch` now handles `.zip` archives.
  6. Credit in README is via `NOTICE`/`LICENSES/` (existing convention). The
     user-facing name "InChI" is a trademark of the InChI Trust. rdoxide
     states that it implements the InChI 1.07.3 algorithm and claims no
     endorsement.

## 2. Architecture

```
crates/rdoxide/src/external/
  inchi_api.rs          RDKit External/INCHI-API/inchi.{h,cpp} (BSD-3 Novartis):
                        ExtraInchiReturnValues, mol_to_inchi, mol_block_to_inchi,
                        inchi_to_mol, inchi_to_inchi_key, mol_to_inchi_key,
                        get_inchi_version
  inchi/                InChI Software 1.07.3 (MIT), one Rust module per C file,
    mod.rs              constants from inchi_api.h / bcf_s.h, InchiRet codes
    ikey_dll.rs         ikey_dll.c (+ extract_inchi_substring of util.c)   DONE
    ikey_base26.rs      ikey_base26.c                                       DONE
    sha256.rs           original FIPS 180-4 implementation (not sha2.c)     DONE
    ...                 see the map below
crates/rdoxide-python/src/rdinchi.rs   Wrap/pyInchi.cpp -> rdoxide.Chem.rdinchi
python/rdoxide/Chem/rdinchi.py         shim to the native module
python/rdoxide/Chem/inchi.py           provisional (see below), then verbatim
                                       External/INCHI-API/python/inchi.py
```

The modules keep the C file names (`ichican2.rs`, `ichi_bns.rs`, ...) and C
function names (snake_cased), as the RingDecomposerLib port does. That way,
when RDKit moves to a newer InChI, the upstream diff can be applied file by
file.

Error model: `inchi_api` returns `RDKitError::NotImplemented` (Python:
`NotImplementedError`) for the not-yet-ported conversions. The library-level
functions return the InChI return codes (`InchiRet`, `InchiKeyError`) like the
C API. The C library's own messages go into `ExtraInchiReturnValues.message` /
`log` verbatim.

Python, current state: `rdinchi` has the full upstream surface.
`InchiToInchiKey` and `GetInchiVersion` work; the rest raise
`NotImplementedError`. `Chem/inchi.py` is provisional and matches RDKit built
without InChI (`noinchi.py`): `INCHI_AVAILABLE = False` and no
`MolToInchi`/`MolFromInchi`, so upstream tests skip. It also exports
`InchiToInchiKey`, `GetInchiVersion` and `InchiReadWriteError`. It becomes
the verbatim upstream file (via `porting/copies.toml`) once both directions
work.

## 3. Upstream source map and port order

Line counts are raw C lines (comments included) of InChI 1.07.3 as compiled
by RDKit. Expect about 0.6–0.8 Rust lines per C line. The order follows the
data flow of `GetINCHI` (structure → InChI), then `GetStructFromINCHI`
(InChI → structure).

| # | Upstream C sources | Lines | Role | Rust module(s) |
|---|--------------------|------:|------|----------------|
| 0 | `ikey_dll.c`, `ikey_base26.c`, `sha2.c` (not ported, replaced) | 2.5k | InChIKey from InChI string | `ikey_dll`, `ikey_base26`, `sha256` ✅ |
| 1 | headers `mode.h`, `incomdef.h`, `ichi.h`, `inpdef.h`, `extr_ct.h`, `ichicomn.h`, `ichicant.h`, `ichitaut.h`, `ichi_bns.h`, `ichirvrs.h`, `inchi_api.h` (structs/constants); `util.c`, `ichisort.c`, `ichiring.c`, `ichiqueu.c`, `ichierr.c`, `ichi_io.c` | 7.2k (+12.6k headers) | data model (`inp_ATOM`, `sp_ATOM`, `INChI`, `INChI_Aux`, `T_GROUP_INFO`, `BN_STRUCT`), element table and valences, InChI's own `inchi_qsort`/insertion sorts (deterministic, no libc dependence), ring/BFS queues, error strings, in-memory output streams | `mode`, `types`, `util`, `ichisort`, `ichiring`, `ichiqueu`, `ichierr`, `ichi_io` |
| 2 | `inchi_dll.c` (input side: `ExtractOneStructure`, `SetAtomProperties`, ...), `inchi_dll_a.c` (`parse_options_string`), `ichiparm.c` (option parsing), RDKit `inchi.cpp` `MolToInchi` input building (atoms, bonds, `inchi_Stereo0D` from RDKit stereo, `fixOptionSymbol`) | 8.9k | `inchi_Input` from an RDKit molecule; options | `inchi_dll`, `ichiparm`, `inchi_api` |
| 3 | `runichi.c`, `runichi2.c`, `runichi3.c` (`PreprocessOneStructure`, unusual valences, polymer units), `strutil.c` (components, metal disconnection / `RecMet`, salts) | 16.2k | preprocessing and the per-structure driver `ProcessOneStructureEx` | `runichi*`, `strutil` |
| 4 | `ichinorm.c`, `ichi_bns.c`, `ichitaut.c`, `ichiisot.c` | 27.5k | normalization: ring systems, aromatic bonds, balanced-network search (charges, radicals, mobile H), tautomeric groups, isotopic keys | `ichinorm`, `ichi_bns`, `ichitaut`, `ichiisot` |
| 5 | `ichister.c`, `ichicans.c` | 7.3k | stereo perception (0D parities, 2D/3D geometry) and stereo canonicalization | `ichister`, `ichicans` |
| 6 | `ichican2.c`, `ichicano.c`, `ichimap1.c`, `ichimap2.c`, `ichimap4.c` | 16.5k | canonical numbering (`CanonGraph`, connection tables, symmetry/stereo mapping) | `ichican2`, `ichicano`, `ichimap*` |
| 7 | `ichimake.c`, `ichimak2.c`, `inchi_dll_a2.c`, `runichi4.c` | 12.7k | `INChI`/`INChI_Aux` construction, Hill formula, component sorting, composite structures | `ichimake`, `ichimak2`, `inchi_dll_a2`, `runichi4` |
| 8 | `ichiprt1.c`, `ichiprt2.c`, `ichiprt3.c` | 12.3k | InChI string + AuxInfo serialization, warnings/messages | `ichiprt*` |
| 9 | `mol_fmt1.c`, `mol_fmt2.c`, `mol_fmt3.c`, `mol2atom.c`, `readinch.c`, `inchi_dll_b.c` | 11.0k | InChI's own MOL V2000/V3000 reader (`MakeINCHIFromMolfileText`, used by `MolBlockToInchi`) | `mol_fmt*`, `mol2atom`, `readinch`, `inchi_dll_b` |
| 10 | `ichiread.c` | 12.5k | InChI string parser (`InChILine2Data`) | `ichiread` |
| 11 | `ichirvr1.c`–`ichirvr7.c`, `ichilnct.c`, `inchi_dll.c` (`GetStructFromINCHI`) | 31.1k | structure restoration from InChI | `ichirvr*`, `ichilnct` |
| 12 | RDKit `inchi.cpp` `InchiToMol` + `cleanUp`/`_Valence*CleanUp*`/`rCleanUp` | 0.6k | RDKit molecule from restored structure | `inchi_api` |
| — | not ported: `ixa/*.c` (IXA API, compiled by RDKit but unused by `inchi.cpp`), `mol_fmt4.c` (SD writer, executable only), `inchi_gui.c`, `bcf_s.c`, `INCHI_EXE`, `sha2.c` | ~9k | | |

Total to port: about 165k C lines, so roughly 100–130k Rust lines. This is the
largest single subsystem after MolDraw2D.

## 4. Stages (PR-sized)

Stage numbers in `docs/STATUS.md` are assigned when each PR lands. The
InChI-local labels are:

| Stage | Content | Test gate |
|-------|---------|-----------|
| I0 ✅ | license audit, skeleton, InChIKey (`ikey_dll`, `ikey_base26`, own SHA-256), `rdinchi` bindings, provisional `Chem.inchi`, fixture corpus, `inchi`/`inchikey` differential modes | `crates/rdoxide/tests/inchi.rs`, `tests/python/test_inchi.py`, differential `inchikey` 588/588 (+20558/20558 on the full default corpus) |
| I1 | data model + utilities (row 1). Done so far: `ichisort` (`inchi_qsort`, the insertion sorts, the neighbor-list comparisons) `ichierr` and the structure-independent part of `util` (element table and helpers); `inpdef` (`inp_ATOM`) and `ichiring` (queue, smallest ring through a bond, 3-membered ring test; `tools/inchi_ring_reference.py`); the `inp_ATOM` helpers of `util` (`nBondsValToMetal`, `num_of_H`, `nNoMetal*`, ion groups, ilists; `tools/inchi_atoms_reference.py`); the string helpers of `util` (`normalize_string`, `lrtrim`, `read_upto_delim`, `extract_inchi_substring`, `MakeRemovedProtonsString`, ...; `tools/inchi_strings_reference.py`), which completes `util`; `ichi_io` (string/file streams, line readers, print buffers, `Output_RecordInfo`; `tools/inchi_io_reference.py`; formatting is done by the caller, so `GetMaxPrintfLength` is not translated); the remaining plain structs (`sp_ATOM`, `INChI`, `T_GROUP_INFO`, `BN_STRUCT`) have no behavior to test and are defined together with the first code that uses them (I2 onwards), so I1 is complete, `ichi_io`, the structs. (`ichiqueu.c`, listed with them, is the tautomeric path search; it moves to the `ichitaut` stage I4) | unit tests of element table/valences; sorts against the output of the C code (`tools/inchi_sort_reference.py`, `tests/inchi_utils.rs`, 372 cases incl. the order of ties; `ErrMsg`/`AddErrorMessage` against `tools/inchi_err_reference.py`) |
| I2 | input from RDKit mol + options (row 2). Done so far: `ichiparm::set_common_options_by_parg` (`tools/inchi_parm_reference.py`, 568 cases) and `parse_options_string` (`tools/inchi_optstr_reference.py`, 600 cases); `ReadCommandLineParms` as built into the library (`ichiparm::read_command_line_parms`, `INPUT_PARMS`, no GUI/developer options; `tools/inchi_readparms_reference.py`, 700 cases; the C `errno` that the timeout option tests but never resets is modelled); next `GetINCHIEx` option use (`enforce_std_format`),  `ReadCommandLineParms` pieces, the `inchi_Input` building | dump of `inchi_Input` vs a reference trace (C sources built with a print hook, `tools/reference_trace`) |
| I3 | preprocessing driver (row 3) | trace of preprocessed `inp_ATOM` arrays |
| I4a/b | `ichinorm` + `ichitaut`/`ichiisot`; then `ichi_bns` alone (12.5k lines) | traces of normalized structures, t-groups |
| I5 | stereo (row 5) | trace of parities |
| I6 | canonicalization (row 6) | trace of canonical ranks |
| I7 | INChI/Aux construction (row 7) | trace of `INChI` structs |
| I8 | serialization (row 8): **first end-to-end `MolToInchi`**, `MolToInchiKey`, `MolToInchiAndAuxInfo` | `mol_to_inchi_matches_rdkit` exact on the fixture (InChI, AuxInfo, return code, message); differential `inchi` on 10312 molecules |
| I9 | `MolBlockToInchi` (row 9) | fixture `molblock` records, upstream `test_data/*.mol` |
| I10 | InChI parser (row 10) | parse/round-trip of fixture InChIs |
| I11a/b | structure restoration (row 11, 31k lines, split rvr1–3 / rvr4–7) | traces of restored atoms |
| I12 | RDKit `InchiToMol`, `MolFromInchi(AndAuxInfo)`; switch `Chem/inchi.py` to verbatim upstream (`INCHI_AVAILABLE = True`) | upstream `rdkit/Chem/UnitTestInchi.py` and `testMolStandardize.py` InChI test enter the ratchet; differential `molfrominchi` |

Fidelity notes for later stages:
* InChI uses its own `inchi_qsort` and insertion sorts, so the order is
  deterministic and defined by the C code. Port them literally; don't
  substitute Rust's `sort`.
* Floating-point use is limited to stereo geometry (`ichister.c`) and
  coordinate output (`ichiprt1.c` `WriteOrigCoord`, `%f`-style formatting).
  Reproduce the C formatting exactly.
* Global state: the C library keeps some static state (e.g. timeouts, the
  `bLibInchiSemaphore`). In Rust it becomes explicit context structs, so
  that calls are reentrant. RDKit serializes calls with a mutex, which we do
  not need.

## 5. Open questions

* `External-contributors` is referenced by every header but missing from
  the v1.07.3 zip. Ask the InChI Trust or check later tags. It names people
  and does not change the license terms.
* Should rdoxide also offer the newer InChI 1.07.4/1.07.5 behind an option?
  Not planned: RDKit parity requires 1.07.3.
* `INCHI_AVAILABLE` stays `False` until both `MolToInchi` and `MolFromInchi`
  work, because the upstream tests guard on it for both directions. An
  alternative is to flip it after I8 and accept failing reverse tests in the
  ratchet.
