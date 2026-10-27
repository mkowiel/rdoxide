# InChI port

RDKit exposes InChI through `rdkit.Chem.rdinchi` (native) and `rdkit.Chem.inchi`
(Python, star-imported into `rdkit.Chem`): `MolToInchi`, `MolToInchiAndAuxInfo`,
`MolBlockToInchi(AndAuxInfo)`, `MolToInchiKey`, `InchiToInchiKey`, `MolFromInchi`,
`MolFromInchiAndAuxInfo`, `GetInchiVersion`, `InchiReadWriteError`,
`INCHI_AVAILABLE`, plus free-form option strings (`"/FixedH /RecMet"`, ...) that
are passed through to the InChI library. This document covers the license audit
and the plan for a faithful pure-Rust port of the InChI library that RDKit
2026.09.1 builds against.

**Status (2026-10-22): the port is complete for everything RDKit uses**, apart
from the polymer code paths. Generation (`MolToInchi`, `MolBlockToInchi`,
InChIKey) and restoration (`MolFromInchi`, `MolFromInchiAndAuxInfo`) work and
agree with the official RDKit 2026.09.1 and with the C library on every corpus
tried; section 4 lists what was verified how and section 6 what is missing.

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
  inchi_input.rs        inchi.cpp MolToInchi input side (rCleanUp, inchi_Atom and
                        inchi_Stereo0D from a molecule, fixOptionSymbol)
  inchi_to_mol.rs       inchi.cpp InchiToMol (assignBondDirs, findAlternatingBonds,
                        the cleanUp family)
  inchi/                InChI Software 1.07.3 (MIT), one Rust module per C file
                        (about 84k Rust lines for 165k C lines); see the map below
crates/rdoxide-python/src/rdinchi.rs   Wrap/pyInchi.cpp -> rdoxide.Chem.rdinchi
python/rdoxide/Chem/rdinchi.py         shim to the native module
python/rdoxide/Chem/inchi.py           verbatim External/INCHI-API/python/inchi.py
```

The modules keep the C file names (`ichican2.rs`, `ichi_bns.rs`, ...) and C
function names (snake_cased), as the RingDecomposerLib port does. That way,
when RDKit moves to a newer InChI, the upstream diff can be applied file by
file. Exceptions: `ichimap1/2` live in `ichimap.rs`,
`mol_fmt1-3` in `molfmt*.rs`, `ichicans` in `ichister.rs`, `ichiisot` in
`ichican2.rs`, `ichiqueu` in `ichiring.rs` and `ichitaut.rs`, `GetStructFromINCHI`
in `inchi_dll_c.rs`, and the C allocators, `Free*` functions and GUI/display code
have no Rust counterpart.

Error model: the library-level functions return the InChI return codes
(`InchiRet`, `InchiKeyError`) like the C API; the C library's own messages go
into `ExtraInchiReturnValues.message` / `log` verbatim. Global C state
(timeouts, the `bLibInchiSemaphore`) became explicit context structs, so calls
are reentrant.

Python: `rdinchi` has the full upstream surface; `Chem/inchi.py` is the
upstream file (copied by `tools/sync_upstream.py regen`), `INCHI_AVAILABLE` is
`True`, and `MolFromInchiAndAuxInfo` parses the AuxInfo in Python as upstream
does.

## 3. Upstream source map

Line counts are raw C lines (comments included) of InChI 1.07.3 as compiled
by RDKit. The order of the rows follows the data flow of `GetINCHI` (structure
→ InChI), then `GetStructFromINCHI` (InChI → structure). Status: ✅ everything
RDKit can reach is ported; 🟡 ported with the gaps named in the last column;
⬜ not ported.

| # | Upstream C sources | Lines | Role | Rust module(s) | Status | Not ported |
|---|--------------------|------:|------|----------------|--------|------------|
| 0 | `ikey_dll.c`, `ikey_base26.c`, `sha2.c` | 2.5k | InChIKey from the InChI string | `ikey_dll`, `ikey_base26`, `sha256` (own FIPS 180-4 code instead of `sha2.c`) | ✅ | `sha2.c` self-test and HMAC |
| 1 | headers `mode.h`, `incomdef.h`, `ichi.h`, `inpdef.h`, `extr_ct.h`, `ichicomn.h`, `ichicant.h`, `ichitaut.h`, `ichi_bns.h`, `ichirvrs.h`, `inchi_api.h`; `util.c`, `ichisort.c`, `ichiring.c`, `ichiqueu.c`, `ichierr.c`, `ichi_io.c` | 7.2k (+12.6k headers) | data model (`inp_ATOM`, `sp_ATOM`, `INChI`, ...), utilities, sorting, queues, I/O streams | `ichi`, `inpdef`, `incomdef`, `util`, `el_data`, `ichisort`, `ichiring`, `ichierr`, `ichi_io` | ✅ | allocators, Windows text window output, `_inchi_trace` |
| 2 | `inchi_dll.c` (input side), `inchi_dll_a.c` (`parse_options_string`), `ichiparm.c`, RDKit `inchi.cpp` (`MolToInchi` input) | 8.9k | `inchi_Input` from a molecule; option parsing | `inchi_dll`, `ichiparm`, `inchi_input` | ✅ | staged `INCHIGEN_*`/`STDINCHIGEN_*` API, `GetINCHIfromINCHI`, command-line help |
| 3 | `runichi.c`, `runichi2.c`, `runichi3.c`, `strutil.c` | 16.2k | preprocessing (unusual valences, salts, metal disconnection, components) and the per-structure driver `ProcessOneStructureEx` | `runichi`, `runichi2`, `runichi3`, `runichi4`, `strutil` | 🟡 | polymer units (`OAD_Polymer*`, CRU folding, frame shift), AuxInfo input (`InchiToInchi_Input`), display/debug code |
| 4 | `ichinorm.c`, `ichi_bns.c`, `ichitaut.c`, `ichiisot.c` | 27.5k | normalization: ring systems, aromatic bonds, balanced-network search (charges, radicals, mobile H), tautomeric groups (incl. keto-enol and 1,5), isotopic keys | `ichinorm`, `ichi_bns`, `ichitaut`, `ichican2` (isotopic keys) | 🟡 | the ring-chain/underivatization code of `ichinorm.c` (engineering mode), `TempFix_NH_NH_Bonds` and a few `#if`-ed helpers |
| 5 | `ichister.c`, `ichicans.c` | 7.3k | stereo perception (0D parities, 2D/3D geometry) and stereo canonicalization | `ichister` | ✅ | |
| 6 | `ichican2.c`, `ichicano.c`, `ichimap1.c`, `ichimap2.c`, `ichimap4.c` | 16.5k | canonical numbering (`CanonGraph`, connection tables, symmetry/stereo mapping) | `ichican2`, `ichicano`, `ichimap`, `ichimap4` | ✅ | time limit (`InchiClock`, `bInchiTimeIsOver`), user-quit callbacks, `CMODE_NOEQ_STEREO`, CRC32 debug code |
| 7 | `ichimake.c`, `ichimak2.c`, `inchi_dll_a2.c`, `runichi4.c` | 12.7k | `INChI`/`INChI_Aux` construction, Hill formula, component sorting, composite structures | `ichimake`, `ichimak2`, `runichi4` | ✅ | stepwise normalization/canonicalization API of `inchi_dll_a2.c` |
| 8 | `ichiprt1.c`, `ichiprt2.c`, `ichiprt3.c` | 12.3k | InChI string and AuxInfo serialization, warnings/messages | `ichiprt1`, `ichiprt2`, `ichiprt3` | 🟡 | polymer layers (`OutputINCHI_PolymerLayer*`, `OutputAUXINFO_PolymerInfo`) |
| 9 | `mol_fmt1.c`, `mol_fmt2.c`, `mol_fmt3.c`, `mol2atom.c`, `readinch.c`, `inchi_dll_b.c` | 11.0k | InChI's own MOL V2000/V3000 reader (`MakeINCHIFromMolfileText`, used by `MolBlockToInchi`) | `molfmt`, `molfmt3`, `mol2atom`, `readinch`, `inchi_dll_b` | ✅ | `mol_fmt4.c` writer half (S-group lists, SD output; executable only) |
| 10 | `ichiread.c` | 12.5k | InChI string parser (`InChILine2Data`) | `ichiread` | 🟡 | AuxInfo segment parsers (`ParseAuxSegment*`, `ReadInChICoord`'s AuxInfo line) and `OutputInChIAsRequested` |
| 11 | `ichirvr1.c`–`ichirvr7.c`, `ichilnct.c`, `inchi_dll.c` (`GetStructFromINCHI`) | 31.1k | structure restoration from InChI | `ichirvrs`, `ichirvr1`–`ichirvr7`, `inchi_dll_c` | ✅ | display code, `Get_*_Input_FromAuxInfo` |
| 12 | RDKit `inchi.cpp` `InchiToMol` + `cleanUp`/`_Valence*CleanUp*`/`rCleanUp` | 0.6k | RDKit molecule from the restored structure | `inchi_api`, `inchi_to_mol`, `inchi_input` | ✅ | |
| — | `ixa/*.c` (IXA API, compiled by RDKit but unused by `inchi.cpp`), `inchi_gui.c`, `bcf_s.c`, `INCHI_EXE` | ~9k | | | ⬜ | by design: not reachable from RDKit |

Total: about 165k C lines became about 84k Rust lines (the C code repeats
itself across the `CanonGraph01..12` size variants and spends lines on
allocation and display).

## 4. Stages (PR-sized)

Stage numbers in `docs/STATUS.md` are assigned when each PR lands (the rows
125–277 and 281 there). The InChI-local labels are:

| Stage | Content | Status | Test gate |
|-------|---------|--------|-----------|
| I0 | license audit, skeleton, InChIKey (`ikey_dll`, `ikey_base26`, own SHA-256), `rdinchi` bindings, fixture corpus, `inchi`/`inchikey` differential modes | ✅ | `tests/inchi.rs`, `tests/python/test_inchi.py`; differential `inchikey` 588/588 and 20558/20558 |
| I1a–e | data model and utilities: `ichisort`, `ichierr`, `util`, `inpdef`, `ichiring`, `ichi_io` | ✅ | 400–600 random inputs per function vs the C code (`tools/inchi_*_reference.py`) |
| I2a–h | options (`ichiparm`: `set_common_options_by_parg`, `parse_options_string`, `ReadCommandLineParms`, standard-format enforcement), input conversion (`SetAtomProperties`, `Extract0DParities`, `InpAtom0DToInchiAtom`), RDKit's `MolToInchi` input side | ✅ | 568–1500 cases per function vs the C code; AuxInfo of the reference corpus |
| I3a–h | preprocessing (`strutil`): drawing corrections, ion pairs, odd things, bond removal, ammonium and metal salts, metal disconnection, terminal H, components, `PreprocessOneStructure` | ✅ | 1500 random atom arrays per step vs the C code |
| I4a–z | normalization: ring systems, tautomer endpoints, charge groups, the balanced network structure and search (`ichi_bns`, 12.5k lines), proton add/remove, alternating bonds, salt charge groups | ✅ | random graphs/structures per function vs the C code |
| I5a–l | tautomer group rules (1,3; keto-enol; 1,5; 4-pyridinol; pyrazole; tropolone) and stereo perception (`ichister`) | ✅ | 4000–8000 random molecules per function vs the C code |
| I6a–u | t-group sorting, stereo canonicalization (`ichicans`), rank refinement and stereo mapping (`ichimap*`), canonical numbering (`ichican2`, `ichicano`: `CanonGraph`, `GetBaseCanonRanking`, `Canon_INChI3`) | ✅ | 3000–6000 random graphs/tables per function vs the C code |
| I7a–l | `INChI`/`INChI_Aux` construction (`ichimake`, `ichimak2`): `Create_INChI`, comparators, `CompareReversedINChI*`, `CreateOneComponentINChI`, `CreateOneStructureINChI` | ✅ | 3000–4000 random structures vs the C code |
| I8a–c | serialization (`ichiprt2`, `ichiprt3`; `ichiprt1` came with I9) | ✅ | 6000–9000 inputs per layer vs the C code |
| I9a–h | InChI's Molfile reader (`molfmt`, `mol2atom`), the molfile-text → InChI pipeline, `GetINCHI`/`GetStdINCHI`, `MolToInchi`/`MolBlockToInchi`/`MolToInchiKey` wired to RDKit | ✅ | 2200–3900 (fuzzed) molfiles and 3000 `inchi_Input` structures vs the C code; 539-record official-RDKit corpus |
| I10a–i | InChI parser (`ichiread`): `/q`, `/p`, formula, connections, mobile H, stereo, isotopic layers, `ReadInChILine`, `InChILine2Data` | 🟡 | 9800–12700 reader cases vs the C code; AuxInfo segments not ported |
| I11a–z | structure restoration, first half (`ichirvr1`–`ichirvr5`): restoration types, BNS construction for the restore, charge/H fixing functions, `MakeOneInChIOutOfStrFromINChI`, `FixFixedHRestoredStructure`, `FixMobileHRestoredStructure` | ✅ | I11a–j: 41745–135000 random cases per group vs the C code (`tests/inchi_rvr*.rs`); I11k–z end to end only |
| I12a–q | restoration, second half (`ichirvr6`, `ichirvr7`: stereo fixing, `InChI2Atom`, comparisons, `ConvertInChI2Struct`), `GetStructFromINCHI` (`inchi_dll_c`), RDKit `InchiToMol`, `inchi.py` switched to upstream, `MolKey` | ✅ | restored structures identical to C for 1576 InChI; `InchiToMol` identical to the official RDKit on 39584 molecules; upstream `UnitTestInchi` 18/18 |

## 5. Fidelity notes

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


## 6. What is missing or incomplete

* **Polymers**: `MolBlockToInchi` of polymer molfiles agrees with RDKit on the
  oracle corpora (`tests/data/inchi_polymers*.jsonl`, 238 molfiles x all polymer
  option sets): validation, the `/z` layer, `/Z:` and `/U/` AuxInfo, 1.05/1.05+
  cyclization, backbone bonds and CRU frame shift; polymer InChI are read back
  (`tests/data/inchi_polymers_read.jsonl`, `InchiToMol`). Not ported: CRU
  folding (`analyze_CRU_folding`). The official library crashes whenever the
  fold analysis runs on a unit it could fold, because of two memory bugs in
  InChI 1.07.3 itself: `OAD_StructureEdits_Apply` sizes the new atom array as
  `nat0 - del_atom->used` although `del_atom` may hold duplicates (heap overflow),
  and `analyze_CRU_folding` passes the NULL `bkbonds` of a unit that was not
  pre-cyclized (HH/EU units, or units without two caps) to
  `OAD_CollectBackboneBonds` (NULL dereference). rdoxide skips the fold step and
  returns the unfolded result, which equals the C code with folding switched off
  (checked on 198 Zz-capped molfiles against a build of the C sources). RDKit's `MolToInchi` never
  passes polymer data, so only `MolBlockToInchi` can differ.
* **AuxInfo input**: `ParseAuxSegment*`, `ReadInChICoord` (AuxInfo line) and
  `Get_inchi_Input_FromAuxInfo` are not ported; RDKit does not call them
  (`MolFromInchiAndAuxInfo` parses the AuxInfo string in Python).
* **Time limits and user-quit callbacks** (`-W`, `InchiClock`) are accepted but
  never enforced; `CMODE_NOEQ_STEREO` is not supported.
* **Ring-chain/underivatization code** of `ichinorm.c` (engineering mode) and
  the display/debug code are not ported.
* **Restoration of InChI with `/KET`/`/15T` layers** fails, as it does in the C
  library (no reversal of these tautomer rules).
* **Verification gaps**: the restoration functions of `ichirvr4`–`ichirvr7`
  (stages 255–277, e.g. `FixRestoredStructureStereo`, `FixMobileHRestoredStructure`,
  `RunBnsRestore1`) are compared with the C code only through the end-to-end
  tests (1576 InChI; the corpus round trips); a hydrogen lost in
  `RunBnsRestore1` for 6 of these 1576 InChI was found that way (stage 281).
  A per-function harness for them was not built.
* **Not planned**: the IXA API, `GetINCHIfromINCHI`, the staged `INCHIGEN_*`
  API, `inchi_gui`, the InChI executable, newer InChI versions (RDKit parity
  requires 1.07.3).
* **License follow-up**: `External-contributors`, referenced by every InChI
  header, is missing from the v1.07.3 zip; it names people and does not change
  the license terms (ask the InChI Trust or check later tags).
