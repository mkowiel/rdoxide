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
* A5b `pattern.rs` + `aacheck.rs`: augmented atoms (`StringToAugmentedAtom`, `ReadAugmentedAtoms`, `ReadAAPairs`),
  `AAMatch`, `TransformAA`, `AAFix`, `TransformAugmentedAtoms`, `AlternativeAATransformation`,
  `ConvertMoleculeToString`, `AAToString` and `CheckAtoms`, tested with the toolkit's real `checkfgs.chk`/`.trn`
  rule files (180 molecules, 915 pattern strings identical to the C code).
* A3a `canonizer.rs`: `CanonicalNumbering` (property refinement, parity/cis-trans splits, tie breaking, stereo clean-up);
  192 molfiles identical to the C code (`tests/data/avalon_canon.jsonl`). `CanSmiles` itself still needs `SMIToMOL`.
* A3b `canonizer.rs`/`utilities.rs`: `StripMolecule` (+ `NewParity`), `FloodFillFragments`, `SplitIonicBonds`,
  `ForceOctet`, `RemoveInorganicFragments`, `ScrambleMolecule` (151 molfiles identical to the C code,
  `tests/data/avalon_strip.jsonl`). `StripMolecule` rebuilds symbol lists and property lines by prepending, so
  both come out reversed; the corpus only reaches the `M  RGP` branch (the C reader does not turn `M  ALS`
  into symbol lists).
* A6a `smi2mol.rs`: `SmilesToMDLAtom`, `ComputeSmilesValenceAndHCount`, `QuickCheck` (2525 SMILES strings and
  351 molfiles identical to the C code). Still missing for `SMIToMOL`: the main parser (ring closures, branches,
  stereo, coordinates section) and the layout it calls.
* A6b `smi2mol.rs`: `SMIToMOL` without layout (flags 0, `TRUE_DB_STEREO`, `EXPECT_SMARTS`; `DO_LAYOUT` returns an
  error until `layout.c` is ported): 1152 results (384 SMILES x 3 flag sets) identical to the C code
  (`tests/data/avalon_smi2mol.jsonl`).
* A6c `geometry.rs`: all of `geometry.c` (`TransformPoints`, `PointSetMatchTransformation`, `NextGraphNeighbour`,
  `HasCusp`, `PointsIntoRing`, `NextSubstituentPoint`, `TransformGraphSubstituent`), bit-identical to the C code on
  2027 cases (`tests/data/avalon_geometry.jsonl`; the C driver must be built with `-ffp-contract=off`). C quirks are
  kept. `CanSmiles` needs `LayoutMolecule` (`layout.c`, ~5400 lines), which comes next in slices.
* L1 `layout.rs` (layout.c lines 74-735): `RandomCoordinates`, `RecolorMolecule`, `FloodColor`/`FloodClearColor`,
  `MakeBondTrans`, `InvertFragmentColor`, ring node table, `LayoutRingSegment`, `LayoutRings`,
  `HideMetalComplexProblems`, `MergeConnectedFragmentColors` (keeps the `if (color = col2)` typo); 72 molfiles
  identical to the C code including the single precision steps (`tests/data/avalon_layout1.jsonl`).
* L2 `layout.rs` (layout.c lines 736-1176): `GetColoredGraph`, `ColorStrain`, `ImproveBondByFlip`,
  `ImproveBondByStretch`, `ImproveMoleculeByBondFlip`, `ImproveSPAtoms`, `CountColor`; 94 molfiles identical to the C
  code (`tests/data/avalon_layout2.jsonl`). `ColorStrain` adds two single precision products to a double, as C does.
* L3 `layout.rs` (layout.c lines 1177-1727): `LayoutLargeLactams`, `FlipRingCarbonyls`, `SproutRingSubstituents`,
  `LinkRingFragments`, `LinkRemainingFragments`, `BranchQuality`, `ColorSize`; 126 molfiles identical to the C code
  (`tests/data/avalon_layout3.jsonl`; the lactam functions are driven with a synthetic ring size array because
  `RingState` only counts rings per bond).
* L4 `layout.rs` (layout.c lines 1728-2390): `LayoutChainAtoms`, `LayoutRubberFragments`, `BestBondColor`, `Diameter`,
  `OrientFragment`, `LayoutFragments`; 156 molfiles identical to the C code (`tests/data/avalon_layout4.jsonl`).
  Watch for `int * float` products (`200*ap->x`), which C evaluates in single precision.

### L5 — layout.c stereo part

Ported `LayoutAtomStereo`, `ClearFlipFlags`, `ClearDBStereoInSmallRings`, `LayoutBondStereo` and `MakeLandscape`.
Verified against the C code on 155 SMILES (`tests/data/avalon_layout5.jsonl`, driver `lay5.c`, atoms with parity column).

### L6 — layout.c collisions and bridgeheads

Ported `BondLengthSQR`, `AtomBondClash`, `CountCollisions`, `BondLinksAtoms`, `StretchBond`, `ImproveCollisions`,
`ScaleByFixedFragments`, `FixLinearBridgeHeads` and `FixInwardHydrogens`. Float coordinate differences are formed in
single precision as in C (`fd`, and `SQR` of float differences in `ScaleByFixedFragments`). Checked on 90 SMILES
(`tests/data/avalon_layout6.jsonl`, driver `lay6.c`), including randomly placed colliding layouts that exercise the
stretch and random-nudge paths. `FixInwardHydrogens` never fires on this corpus (the test only proves it a no-op there).

### L7 — layout.c fragment bond flips

Ported `TopBottomSeparation`, `AttachmentCenterDistance`, `FloodInvertColor`, `ChangeAtomColors`,
`ImproveFragmentByBondFlip` and `ImproveFragmentByAllBondFlips`. Checked on 56 SMILES with a synthetic coloring
(ring atoms = fragment, others = top/bottom attachments) in `tests/data/avalon_layout7.jsonl` (driver `lay7.c`).
The flip-and-undo trials reproduce C's float round-off since the undo is a second reflection, not a restore.

### L8 — layout.c shortcut chains

Ported `IsShortcutLikeAtom`, `PerceiveShortcutClasses` (without the stderr diagnostics) and `LayoutShortcutChains`.
69 molfiles with `R` + atext atoms and three `max_line_size` values (0, 8, 14) in `tests/data/avalon_layout8.jsonl`
(driver `lay8.c`, records separated by `$$$$`); line splitting changes the classes in 5 of them.

### L9 — layout.c fragment info and alignment

Ported `CollectFragmentInfo` (fragments as chains: a shortcut fragment followed by the atom fragments linked to it),
`ColorScore`, `CollectAttachments`, `TfmPoint`, `ScoreTransformation`, `FindBestTransformation`, `AlignFragment`,
`FragmentRenderedHeight` and `FragmentHeight`. Checked on 69 molfiles in `tests/data/avalon_layout9.jsonl`
(driver `lay9.c`). C quirks kept: the x collision term multiplies by `ap->x - ap->x` (zero), the y term does not;
float differences are squared in single precision.

### L10 — layout.c: AlignShortcutsAndFragments and LayoutMolecule

`layout.c` is complete: `AlignShortcutsAndFragments` and `LayoutMolecule` (`layout_molecule`) are ported, and 221 inputs
(SMILES read with `SMIToMOL` plus shortcut molfiles, `tests/data/avalon_layout10.jsonl`, driver `lay10.c`) give
identical coordinates, wedges and bond marks. Note for C's global `rand()`: `CombineRings` does `srand(1)` (unless the
ring list is empty), so the stream `ImproveCollisions` draws from is the state left by the last ring perception
(`LayoutRings`), not the one after `RandomCoordinates`. `combine_rings_rng` / `layout_rings_rng` carry that state.

### L11 — SMIToMOL with DO_LAYOUT

`smi_to_mol` no longer rejects `DO_LAYOUT`: it recolors, applies an optional coordinate string after `|` (`KEEP_POSITION`
atoms, `DONT_FLIP_BOND` between two of them), calls `layout_molecule` and restores the atom ordinals as colors.
159 inputs incl. coordinate strings with missing/garbled fields in `tests/data/avalon_layout11.jsonl` (driver `smil.c`).

### A3 — CanSmiles

`can_smiles` (`CanSmilesStep`, `UnifyFragmentDBStereo`, `makeComponentSet`) and `remove_standard_components` are ported
in `canonizer.rs`; the scramble seed is the C file-level `seed` (starts at 13, `set_seed`). 627 cases (SMILES x flag
combinations incl. `MAIN_PART|COMPONENT_SET`, `FORCE_OCTET`, `SCRAMBLE`, `TO_MONOMER`, `DY_AROMATICITY`, random SMILES of
the same molecules, and salt removal) in `tests/data/avalon_cansmiles.jsonl` (driver `cs.c`) are identical to the C code.

### A7b — pyAvalonTools: GetCanonSmiles, Generate2DCoords, data-string fingerprints

`rdoxide.Avalon.pyAvalonTools` gained `GetCanonSmiles(mol|data, ...)`, `Generate2DCoords(mol|data, ...)`, the SMILES/molblock
data flavours of `GetAvalonFP`, `GetAvalonCountFP`, `GetAvalonFPAsWords`, and `avalonSSSBits` / `avalonSimilarityBits`
(`External/AvalonTools/nbWrap/pyAvalonTools.cpp`). Only the structure checker (`CheckMolecule`, ...) is still missing.
Quirk kept: the data flavour of `GetAvalonCountFP` honours `bitFlags`, the molecule flavour does not (an upstream overload
mix-up). `tests/data/avalon_pyapi.jsonl` (1248 cases) is made by `tools/make_avalon_pyapi_corpus.py`.
Platform note: the layout uses the C library's `rand()`; the official macOS wheel gives other coordinates than a glibc
build (verified: the C oracle with macOS `rand` reproduces the wheel exactly, the glibc emulation matches Linux). The
corpus therefore takes `Generate2DCoords(data, ...)` from the C oracle with the glibc generator and, for
`Generate2DCoords(mol)`, only compares success/failure and the conformer id with the wheel.

### A5a — ssmatch.c: substructure matcher

`ssmatch.rs` ports the matching code of `ssmatch.c` (the fingerprint part of that file is in `fingerprint.rs`):
`AtomMatch`, `BondMatch`, `RecSSMatch`, `Refine`, `PerceiveAromaticity` (flags `BTF_*` in `bond_type_flags`),
`SSMatchExt`/`SSMatch`, `FGMatch`, `RSSMatch`, `PruneDuplicateAtomSets`. Matches are returned in the order of the C
linked list (last found first). As in C, `SSMatchExt` leaves the aromaticity flags in both molecules and the topography in
the query. 380 molecule/query pairs (queries with `A`/`Q` atoms, bond types 5-8, mapping numbers) in
`tests/data/avalon_ssmatch.jsonl` (driver `ss.c`) are identical to the C code. Next in the struchk chain: `fixcharges.c`,
`patclean.c`, `shortcut.c`, the `struchk.c` driver itself.

### A5b — struchk.c: InitCheckMol and RunStruchk (library options)

`struchk.rs` ports the library entry points of `programs/struchk.c` for the options RDKit's MolKey uses (`-ta -ca -cc -cl
-cn -cs -tm -da -dg -ds -or -l`): the option tokenizer of `InitCheckMol`, `ClearParameters` (the collision limit and the
open log file survive it, as in C), `RunStruchk`, `ConvertAtomAliases`, `ConvertGroupsToSGroups`, `StripSmallFragments`
(incl. `SplitMolecule`) and the global message list. Quirks kept: with `-l` the messages are flushed to the log file at
the end of every check (newest first), so `GetCheckMolLog` is then empty; a failing `-ta`/`-ca` file only prints an
error. The charge (`-tc`), tautomer (`-tt`), pattern (`-pc -pr -ps`) and profile (`-f`) options need `fixcharges.c`,
`patclean.c` and the RD file reader; `init_check_mol` counts them as option errors for now. 306 option set x input
combinations (SMILES, 2D/3D/shrunk/either-bond/alias molfiles) in `tests/data/avalon_struchk.jsonl` (driver `sc.c`)
are identical to the C code. `fixcharges.c` was not needed for this step after all (MolKey does not use `-tc`).

### A7c — pyAvalonTools: structure checker, count fingerprints as sparse vectors

`CheckMolecule`, `CheckMoleculeString`, `InitializeCheckMol`, `CloseCheckMolFiles`, `GetCheckMolLog` and the
`StruChkFlag` / `StruChkResult` enums are exposed (the enums as `enum.IntEnum`); `RDConfig.RDDataDir/struchk` now ships the
three struchk files. `GetAvalonCountFP` returns a `UIntSparseIntVect` like upstream (it returned a list before).
`tests/data/avalon_checkmol.jsonl` (153 inputs, from the official wheel with the MolKey option set) is compared in
`tests/python/test_avalon_checkmol.py`. The upstream tests `External/AvalonTools/Wrap/testAvalonTools.py` (12) and
`Chem/MolKey/UnitTestMolKey.py` (all 16, 6 more than before) are now in the ratchet.

### A5c — fixcharges.c: charge fixing for struchk (`-tc`, `-tq`, `-tl`)

`fixcharges.rs` ports `InitializeChargeDataTables` (the `.aci` acidity table, shipped as `checkfgs.aci`),
`SetpKaValues`, `MarkMostAcidicAtoms`, `RefineAcidicAtoms`, `AllCentersRefined`, `RechargeMolecule` and
`SetAcidityLimit`; `struchk.rs` uses them for `-tc <file>`, `-tq <charge>` and `-tl <limit>` (the C code exits the process
on a malformed table, here the option is counted as an error). Numerics as in C: `float` atom values updated in double,
the cubic pKa transformation in double. Quirks kept: the acidity limit is global and survives `ClearParameters`; `cip`
keeps pointing into the charge table between the increment steps (a bond-order increment without a matched atom uses
0); molecules with symbol lists or property lines are reported as not rechargeable. The logging for parameter fitting
(`SetChargeLog`) is not ported. `tests/data/avalon_fixcharges.jsonl` (226 inputs x option sets, driver `sc.c`) is
identical to the C code. Still open for struchk: `-tt`, `-pc/-pr/-ps` (`patclean.c`), `-f`.

### A5d — patclean.c: template cleaning, rotation and forced stereo (`-pc`, `-pr`, `-ps`)

`patclean.rs` ports `ForceStereoTemplate`, `TemplateClean`, `TemplateRotate` and their helpers (`ClosestMatch`,
`AlignMoleculeWithSubstructure`, `AdjustToSubstructureTemplate`, `ColorSubstituents`, `TransformColoredSubstituent`,
`ScaleByTemplate`, `AddStereoHydrogen`); `struchk.rs` reads the template files (`read_pattern_molecules`) and applies
them after the charge fixing in the C order (stereo, clean, rotate). The templates are modified by the substructure
matcher as in C and kept between molecules. Quirks kept: the template reader skips every second record (the buffer is
already at `$$$$`, `GetBuffer` then eats the next molecule); `STEREO_TRANSFORMED` is never set (C compares the result with
15); `FlipMolecule` also swaps the wedges; `SQR()` of a `float` difference and the sum of two of them are `float`
operations; `ring_size` is read after `free` for the `either` bond fix (first four entries count as non-zero, no corpus
case reaches it). The C code crashes when a later option set frees templates of an earlier one, so the corpus is built with
one process per record. `tests/data/avalon_patclean.jsonl` (196 inputs x option sets, 12 template files in
`tests/data/avalon_patterns/`, driver `sc.c`) is identical to the C code. Still open for struchk: `-tt`, `-f`.

### A5e — struchk: tautomer rules (`-tt`) and option profiles (`-f`)

`-tt <rd file>` reads the source and target pattern of a `$RXN` record (`read_tautomer_rules`; the C loop ends after the
first rule, so only one rule is ever used) and `ApplyTautomer` rewrites bond orders, charges and radicals of the first
match, at most three times per rule, after the fragment removal. Quirks kept: the bond and charge changes are driven by the
atom mapping numbers of the *molecule* (the colours set from the pattern mapping are restored without being used), so a rule
only changes molecules that carry atom maps, while the message "has been tautomerized with rule" is added for every match.
`-f <profile>` makes the following options come from the file (`<option> [<parameter>] [<comment>]`, `#` comments, quoted
parameters) before the rest of the command line, as `GetOption` does in the library path; an undefined option counts as an
error. A profile option without its parameter and a missing rule file crash or exit the C code and are errors here.
`tests/data/avalon_tautomer.jsonl` (84 inputs x option sets, rule and profile files in `tests/data/avalon_patterns/`,
driver `sc.c`) is identical to the C code. All struchk options of the library are ported now.
