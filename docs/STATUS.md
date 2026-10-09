# Port status

Upstream: RDKit 2026.09.1. Legend: ✅ done · 🟡 partial · ⬜ not started.

| Stage | Component | Status | Notes |
|-------|-----------|--------|-------|
| 0 | Workspace, packaging (maturin, abi3-py310), CI, licenses, NOTICE | ✅ | |
| 0 | `tools/sync_upstream.py` (fetch/regen/diff/stamp), `porting/manifest.toml` | ✅ | |
| 0 | Multi-core: `rdgeneral::threads` (`getNumThreadsToUse`, order-preserving `par_map` on scoped threads, no dependency); `rdFingerprintGenerator.Get*Fingerprints(mols, numThreads)` run in parallel with the GIL released | ✅ | results independent of the thread count; 20k molecules RDKit-FP 2.4s -> 0.4s on 14 cores. Other batch APIs adopt `par_map` as they are ported |
| 0 | Dependency license audit (`tools/check_licenses.py`, `LICENSES/CRATES.md`, CI step) | ✅ | permissive-only allow-list; 53 crates |
| 0 | Upstream test ratchet (`tools/run_upstream_tests.py`) | ✅ | 298 upstream tests passing |
| 1 | RDGeneral: Dict/RDValue/RDProps | ✅ | |
| 1 | RDGeneral: RDLog + `rdBase` logging API | ✅ | upstream `UnitTestLogging.py` 32/32 |
| 1 | RDGeneral: RNG (`boost::minstd_rand`, glibc `rand`) | ✅ | |
| 1 | `rdBase` vector types `_vecti`, `_vectj`, `_vectd`, list of int vectors (mutable sequences without value equality, like Boost.Python's indexing suites) | 🟡 | used by DataStructs and Subgraphs; string/nested-vector variants later |
| 2 | PeriodicTable (+ data extraction from upstream) | ✅ | |
| 2 | Atom, Bond, Conformer, ROMol/RWMol, RingInfo, StereoGroup (core) | ✅ | boost adjacency ordering, RWMol edits, batch edits, valence calc |
| 2 | `rdchem` Python: Atom, Bond, Mol, RWMol, EditableMol, Conformer, RingInfo, StereoGroup, enums, props API, sanitize exceptions | ✅ | `_wrapperType = 'nanobind'` (PyO3 semantics match nanobind's) |
| 2 | `rdGeometry` Point3D / Point2D | 🟡 | PointND, UniformGrid3D later |
| 2b | SubstanceGroup, MonomerInfo / PDBResidueInfo, bookmarks API in Python, QueryAtom/QueryBond | 🟡 | SubstanceGroup in stage 13e; MonomerInfo / PDBResidueInfo outstanding |
| 3 | Ring finding: RingDecomposerLib port, findSSSR, symmetrizeSSSR (RDL + legacy), fastFindRings, ring families | ✅ | differential: 10312/10312 identical to RDKit 2026.09.1 |
| 3 | sanitizeMol: cleanUp, cleanUpOrganometallics, Kekulize (canonical + input order), aromaticity (RDKit/simple/MDL/MMFF94), conjugation, hybridization, radicals, adjustHs, chirality/atropisomer cleanup | ✅ | differential `sanitize` 10312/10312 |
| 3 | Canonical atom ranking (new_canon: rankMolAtoms, rankFragmentAtoms, chiralRankMolAtoms); `CanonicalRankAtoms(InFragment)` | ✅ | differential `ranks` 10312/10312 |
| 3 | removeHs/removeAllHs (`RemoveHsParameters`), legacy assignStereochemistry (CIP ranks/codes, ring stereo, E/Z), setDoubleBondNeighborDirections | ✅ | differential `molfromsmiles` 10312/10312 + 35140 stereo-rich + 72 edge cases; SubstanceGroup checks wait for stage 2b |
| 3 | Stereo from coordinates: assignChiralTypesFromBondDirs (pseudo-3D wedges), assignChiralTypesFrom3D (incl. non-tetrahedral), detectBondStereochemistry, atropisomer perception (2D/3D/no conformer) | ✅ | differential `molblock` |
| 3 | addHs (`AddHs`, `AddHsParameters`, `onlyOnAtoms`, explicit/implicit Hs with coordinates from `Transform3D`; residue info waits for MonomerInfo) | ✅ | differential `addhs` 5000/5000, `addhs3d` 7706/7721 (the 15 are upstream molblocks with SGroups, which the reader cannot read yet) |
| 3 | new (non-legacy) stereo perception, CIPLabeler | ⬜ | `SetUseLegacyStereoPerception(False)` currently raises |
| 4 | SMILES parser (lexer + grammar + SmilesParseOps, non-tetrahedral stereo) | ✅ | differential: 10312/10312 identical to RDKit (sanitize=False) |
| 4 | `MolFromSmiles` sanitize/removeHs/stereo path (defaults) | ✅ | differential `molfromsmiles` identical |
| 4 | SMILES writer + canonical SMILES (Canon.cpp DFS, SmilesWrite.cpp, getMolFrags/copyMolSubset): `MolToSmiles` (both overloads), `MolFragmentToSmiles`, `MolToRandomSmilesVect`, `SmilesWriteParams` | ✅ | differential `smiles` (canonical, output order, non-canonical, kekule, non-isomeric, explicit, rooted, random, fragment, roundtrip) 10312 + 35973 identical |
| 4 | Query model (Query/*, QueryOps, QueryAtom/QueryBond expandQuery, null-query algebra), SMARTS parser (`MolFromSmarts`, `AtomFromSmarts`, `BondFromSmarts`, `SmartsParserParams`), `DescribeQuery` | ✅ | differential `smarts` (atoms, bonds, query descriptions) 1275 library SMARTS + 151 edge cases identical |
| 4 | SMARTS writer: `MolToSmarts` (both overloads), `MolFragmentToSmarts`, `Atom/Bond.GetSmarts` | ✅ | differential `smarts` (+ writer: default, non-isomeric, rooted, fragment, per-atom/bond, roundtrip) 1426/1426; `smiles` mode also checks `MolToSmarts` of plain molecules |
| 4 | Python `QueryAtom`/`QueryBond` classes, `rdqueries` | ✅ | stage 13c; `HasBitVectPropWithValueQueryAtom` not ported (no bit-vector property queries yet) |
| 4 | CXSMILES (read/write) | ⬜ | |
| 5 | Substructure matching: VF2, `SubstructMatch`/`SubstructMatchCount`, recursive queries, chirality + double-bond stereo, enhanced stereo, query-query matches, aromaticMatchesConjugated/SingleOrDouble, atom/bond properties; `Mol.HasSubstructMatch/GetSubstructMatch/GetSubstructMatches` (both overloads), `SubstructMatchParameters` | ✅ | differential `substruct` (55 queries x 4 variants per molecule + chirality self/mirror) 10312 + 2418 identical |
| 5 | `extraFinalCheck`/`extraAtomCheck`/`extraBondCheck` callbacks, generic matchers, MolBundle, ResonanceMolSupplier | ⬜ | |
| 6 | MolPickler (byte-identical pickles; reads versions >= 2 incl. queries, rings, stereo groups, conformers, props): `Mol.ToBinary`, `Mol(pkl[, flags])`, Python `pickle`, `Get/SetDefaultPickleProperties`, `PropertyPickleOptions` | ✅ | differential `pickle`/`picklesmarts` byte-identical on 10312 + 1101 + 1426 + edge cases; upstream canonSmiles.v2.pkl stream reads |
| 6 | Pickles with SubstanceGroups / monomer info, v1 pickles, ExplicitBitVect property handler | ⬜ | raise RuntimeError for now |
| 6 | MolBlock reader, V2000 + V3000 (atom lists, NOT lists, query props HCOUNT/UNSAT/RBCNT/SUBST, RGROUPS, enhanced stereo collections, LINKNODE, OBJ3D skipping, continuation lines): atoms/bonds incl. query atoms (A/Q/X/M/R/*, generic groups, H count), query bonds (5-8, topology), M CHG/RAD/ISO/RGP/ALS/RBC/SUB/UNS/PXA, old atom lists, aliases/values, S SKP, 2D/3D stereo; `MolFromMolBlock`, `MolFromMolFile` | ✅ | differential `molblock` (4 parse variants, byte-identical pickles with all props): 2D 4991/4991, V3000 1876/1876 (116 with stereo groups), 3D 490/490, upstream test_data 349/364 (remaining diffs need SGroups) |
| 6 | SDMolSupplier (lazy + random access, `$$$$` indexing, data fields, atom/bond property lists, error recovery; Python iteration/indexing/len/SetData/GetItemText/context manager) | ✅ | differential `sdf` (4 parse variants + len/getitem/GetItemText/SetData) on all 389 upstream SD files: identical except 43 files that need SGroups |
| 6 | MolBlock writer: `MolToMolBlock` (both overloads), `MolToV3KMolBlock`, `MolToV2KMolBlock`, `MolToMolFile`, `MolToV3KMolFile`, `MolWriterParams`; wedging (`pickBondsToWedge`, libc++ `std::sort` tie order), crossed double bonds, parity, query atoms/bonds, atom lists, R groups, enhanced stereo collections | ✅ | differential `molwrite` (5 writer variants x 3 reader variants + roundtrip): 2D 4991/4991, V3000 1876/1876, 3D 490/490, upstream test_data 327/364 (rest need SGroups) |
| 6 | SDWriter (file name, stdout, Python file objects; SetProps, kekulize, V3000, GetText) | ✅ | differential `sdwrite`: identical on all upstream SD files except those needing SGroups |
| 6 | Atropisomer wedging (`wedgeBondsFromAtropisomers`: no conformer, 2D, 3D) | ✅ | covered by `sdwrite` (all upstream atropisomer SD files identical) |
| 6 | SGroups (V2000/V3000; also needed for M ZBO/ZCH/HYD/MRV and for writing SMARTS-query / generic-group / zero-order-bond molecules), M APO/LIN, ATTCHORD templates, ForwardSDMolSupplier | ⬜ | reader returns None, writer raises RuntimeError for now |
| 7 | 2D coordinates: RDDepict::compute2DCoords (EmbeddedFrag, DepictUtils, fused rings, cis/trans and non-tetrahedral fragments, collision removal: bond/spiro flips, open angles, shorten bonds, path angle expansion, canonical orientation, coordMap, bondLength); `rdDepictor.Compute2DCoords`; `MolToMolBlock` of molecules without conformers | ✅ | differential `depict` (canonOrient on/off coordinates to 6 decimals + MolBlock): NCI 4999/4999, canonSmiles.long + stereo/pickle edge cases 5400/5400 identical to the arm64 build |
| 7 | Ring-system templates, random flip sampling (nSample/nFlipsPerSample), Compute2DCoordsMimicDistmat, generateDepictionMatching2D/3DStructure, straighten/normalizeDepiction, CoordGen | ⬜ | raise RuntimeError |
| 8 | DataStructs: ExplicitBitVect, SparseBitVect (bit ops, pickles byte-identical, base64, indexing), all similarity metrics + Bulk*/\*Neighbors, Tversky, on/off-bit projections, FoldFingerprint, AllProbeBitsMatch, text/FPS/binary conversions, ConvertToNumpyArray; `rdkit.DataStructs` Python package (verbatim upstream pure-Python modules) | ✅ | 3000 random pairs identical to the arm64 build (x86 differs in last bits through FMA, as upstream does); upstream `testBV.py` and `UnitTestcBitVect.py` pass. Lists of ints/doubles are `rdBase._vecti`/`_vectd` like upstream |
| 8 | SparseIntVect (Int/Long/UInt/ULong): fuzzy &/|, arithmetic, scalar ops (incl. upstream's stored zeros and missing `/=`), pickles, Dice/Tanimoto/Tversky (+ bounds, Bulk*), ConvertToNumpyArray | ✅ | 3000 random pairs identical to the arm64 build; upstream `testSparseIntVect.py` 10/10 incl. legacy pickles |
| 8 | DiscreteValueVect (+ `DiscreteValueType`), RealValueVect, ComputeL1Norm, ConvertToNumpyArray for all vector types | ✅ | 2000 random pairs identical (incl. error cases); upstream tests 12/12 |
| 8 | FPBReader, InitFromDaylightString | ⬜ | |
| 9 | Fingerprint generator framework (bit / sparse / count / sparse-count fingerprints, count simulation, numBitsPerFeature via boost mersenne_twister (n=4, m=2), fromAtoms, custom invariants, AdditionalOutput incl. atomsPerBit, GetOptions views, bulk + NumPy) and the Morgan generator (`GetMorganGenerator`, `GetMorganAtomInvGen`, `GetMorganBondInvGen`); `rdkit.Chem.rdFingerprintGenerator` | ✅ | differential `morgan` (6 configs + all output kinds): ZINC 10000 + ChEMBL chiral 5000 identical. SparseBitVect now iterates in signed (`std::set<int>`) order |
| 9 | Atom-pair generator (`GetAtomPairGenerator`, `GetAtomPairAtomInvGen`, `AtomPairFingerprintOptions`; 2D and 3D distances, hashed and unhashed codes), `get3DDistanceMat` | ✅ | differential `atompair` (5 configs + all output kinds): 15000/15000 identical; 3D (`use2D=False`) 490/490 |
| 9 | Subgraphs: findAllSubgraphsOfLengthN / MtoN, findUniqueSubgraphsOfLengthN (path discriminators), findAllPathsOfLengthN / MtoN (incl. onlyShortestPaths, rooted), findAtomEnvironmentOfRadiusN; rdmolops wrappers | ✅ | differential `subgraphs` 5000/5000 identical. `PathToSubmol` waits for copyMolSubset(BONDS) |
| 9 | Topological-torsion generator (`GetTopologicalTorsionGenerator`, `TopologicalTorsionFingerprintOptions`; hashed and unhashed codes, onlyShortestPaths, bitPaths) | ✅ | differential `torsion` (6 configs + all output kinds): 15000/15000 identical |
| 9 | RDKit FP generator (`GetRDKitFPGenerator`, `GetRDKitAtomInvGen`, `RDKitFingerprintOptions`; branched/linear paths, bond orders, query bonds, numBitsPerFeature via upstream's 4-word Mersenne twister) | ✅ | differential `rdkitfp` (6 configs + all output kinds): 5000/5000 identical |
| 9 | Legacy `Chem.RDKFingerprint` (tgtDensity folding, atomBits, bitInfo) and `UnfoldedRDKFingerprintCountBased` | ✅ | differential `rdkfingerprint` 1300/1300 identical |
| 9 | `rdMolDescriptors` legacy fingerprints: GetMorganFingerprint / GetHashedMorganFingerprint / GetMorganFingerprintAsBitVect (useFeatures, bitInfo), GetConnectivityInvariants, GetFeatureInvariants, GetAtomPair(Atom)Code, (Hashed) atom pairs and topological torsions (+AsBitVect), deprecation warnings; `GetMorganFeatureAtomInvGen` | ✅ | differential `moldescfp` 5000/5000 identical |
| 9 | MACCS keys (`rdMolDescriptors.GetMACCSKeysFingerprint`, verbatim `Chem/MACCSkeys.py`); `maccs.rs` is generated from upstream `MACCS.cpp` by `tools/gen_maccs.py` (then `cargo fmt`) | ✅ | differential `maccs` 15300/15300 identical |
| 12b | `UniformGrid3D`/`UniformRealValueGrid3D` (binary pickles, set operations, sphere occupancy, `.grd` export), grid utilities (Tversky/Tanimoto/protrude, centroid, terminal points); `PointND`, Point2D angle/in-place operations, `ComputeDihedralAngle`/`ComputeSignedDihedralAngle` | ✅ | upstream `testGeometry.py` 20/20; `test_grids.py` identical on RDKit 2026.09.1 |
| 12c | `SmilesWriter` (header, name and property columns, delimiter, isomeric/kekule options; file, stdout and Python file-object sinks shared with `SDWriter` via `writer_sink`; closed-writer errors identical to upstream's) | ✅ | probe transcripts (6 layouts, file, closed writer) identical to RDKit incl. error text; `test_smiles_writer.py` |
| 12d | `ForwardSDMolSupplier` (binary file objects incl. gzip, file names; shares the record engine of `SDMolSupplier`; forward-only iteration, `atEnd`, `GetEOFHitOnRead`, property-list switch, context manager, closed-supplier errors); `Mol.GetPropNames` now returns RDKit's string vector type; log-message newlines fixed | ✅ | 40-line probe transcript identical to RDKit (damaged records, missing final `$$$$`, trailing blanks, text-mode rejection); `test_forward_sd_supplier.py` |
| 12e | `MQNs_` (42 counts); `CalcAUTOCORR2D` (192 values, or 32 for a custom atom property) with the Dragon-based weight tables, intrinsic state weights (shared with the later 3D descriptors) and upstream's dummy-atom table read; new license file for Godin's descriptor code | ✅ | 62 molecules x 2 descriptors + custom-property/error cases identical to RDKit; `test_autocorr_mqn.py` |
| 12f | `rdMolTransforms` (centroid, principal axes/moments from inertia and gyration matrices with upstream's caching, canonical transforms, `CanonicalizeConformer/Mol`, `TransformConformer`, get/set bond length, angle and dihedral); `numerics::eigen_3x3` reproduces Eigen's `SelfAdjointEigenSolver<Matrix3d>` (vector signs and order) with the clang FMA pattern; invariant error texts now end with upstream's `BOOST:` line and range errors name the C++ source line | ✅ | upstream `testMolTransforms.py` 8/8; 10 molecules x flags compared at full precision (last-digit differences only on degenerate spectra); `test_moltransforms.py` |
| 12g | 3D shape descriptors on `rdMolDescriptors`: `CalcPBF`, `CalcNPR1/2`, `CalcPMI1/2/3`, `CalcRadiusOfGyration`, `CalcInertialShapeFactor`, `CalcEccentricity`, `CalcAsphericity`, `CalcSpherocityIndex` (moment caching, `Descriptors3D` verbatim module); conformer lookups raise upstream's `ValueError: Bad Conformer Id` | ✅ | 5 molecules x 11 descriptors x flag variants identical to RDKit (benzene: last-digit differences from its degenerate spectrum); `test_pmi.py` |
| 12h | `CalcAUTOCORR3D`, `CalcRDF`, `CalcMORSE` (default and `CustomAtomProperty` variants) on the shared atom weights; the 3D distance matrix is cached as `_3DDistanceMatrix_Conf<id>` like upstream | ✅ | 10 molecules x with/without Hs x custom/missing property, dummy-atom and error cases identical to RDKit; `test_3d_profiles.py`; upstream test3D.py gains 6 passes (332/508) |
| 12i | `CalcCoulombMat`, `GetUSR`, `GetUSRCAT` (default SMARTS feature classes or 1-based atom selections), `GetUSRDistributions[FromPoints]`, `GetUSRFromDistributions`, `GetUSRScore`; Boost-style `TypeError` texts via shared `expect_list`/`no_converter` helpers | ✅ | 10 molecules x with/without Hs x selections, error cases and distribution helpers identical to RDKit; `test_usr.py`; upstream testMolDescriptors gains testUSR/testUSRCAT (334/508) |
| 12j | `CalcWHIM` (114 Dragon values or 17 for a custom property), `CalcEEMcharges`; Eigen ports `FullPivLU` (+ its `solve`, panel/FMA order of the triangular solves) and `JacobiSVD` for 3x3 | ✅ | 10 molecules x with/without Hs x threshold/custom property: WHIM and EEM output bit-identical to RDKit, including the noise-dominated terms of planar molecules; `test_whim_eem.py`; upstream test3D.py gains test6WHIM (335/508) |
| 12k | `CalcGETAWAY` (273 values, or 45 for a custom property) with `numerics::dense::Matrix` (Eigen-order products) and `JacobiSVD` generalized to n x n with `U` and the vectorized column rotations | ✅ | 10 molecules x with/without Hs x precision/custom property identical to RDKit (including NaN degenerate cases and the cached matrices); `test_getaway.py`; upstream test3D.py gains test3GETAWAY (337/508) |
| 12l | `GetMolFrags` (asMols/sanitizeFrags/frags/fragsMolAtomMapping), `CombineMols`, `RenumberAtoms` (conformers, rings, ring-stereo props, stereo groups), `GetFormalCharge`, `GetShortestPath`, `RemoveStereochemistry`; shared `pyconv::python_object_to_vect`, Boost `ArgumentError`, int-vector props surface as `_vecti` | ✅ | probe diff against RDKit identical on 12 molecules (CXSMILES-dependent cases excluded: CXSMILES not yet ported); `test_mol_structure.py`, `tests/mol_structure.rs`. |
| 12m | `ReplaceSubstructs`, `ReplaceSidechains`, `ReplaceCore` (both overloads: query core, explicit match vector; ring-stereo repair, chirality fix-ups, dummy coordinates, stereo groups), `MurckoDecompose`; `chem_transforms` split into `mod.rs` + `replace_core.rs` | ✅ | 2000 molecules x 26 operations identical to RDKit, except the order of equally-labelled dummies in unlabelled `ReplaceCore` (upstream breaks ties by `Atom*` address, i.e. allocator-dependent; we use creation order); `test_chem_transforms.py`, `tests/chem_transforms.rs`. |
| 12n | `FragmentOnBonds` (dummy labels, bond types, cuts per atom, stereo/query/conformer transfer), `FragmentOnSomeBonds` (combinations built in upstream's bit-permutation order, fragmented in parallel with `par_map`), `FragmentOnBRICSBonds`; BRICS atom/bond definition tables and the generic `constructFragmenter*Types` parsers | ✅ | 2000 molecules (BRICS, bond subsets with/without dummies) and ~700 hand-picked cases incl. error texts identical to RDKit; `test_fragment_on_bonds.py`, `tests/mol_fragmenter.rs`. |
| 12o | `molzip` (atom map / isotope / atom type / atom property / `FragmentOnBonds` labels, link bonds, tetrahedral and double-bond stereo restoration, wedge-direction fusion, `alignCoordinates`), `molzipFragments` (R group decompositions, `generateCoordinates` via the 2D depictor), R group row zipping; `MolzipParams`/`MolzipLabel`; `Transform3D` matrix product; `RDKitError::{invariant,precondition}_at` for errors that name the upstream line | ✅ | non-bulk probe cases identical to RDKit incl. error texts; bulk `FragmentOnBonds` round trips only comparable up to the first ~50 molecules because upstream reads freed atoms (undefined behaviour) when a double bond is cut at both ends and then crashes, so we preserve the stereo properly there; `test_molzip.py`, `tests/mol_zip.rs`. |
| 12p | Bindings for the stereo/atom helpers the core already had: `AssignChiralTypesFromBondDirs`, `AssignStereochemistryFrom3D`, `AssignAtomChiralTagsFromStructure`, `DetectBondStereo(c)hemistry`, `SetBondStereoFromDirections`, `SetTerminalAtomCoords`, `CleanupChirality/Atropisomers/Organometallics`, `CountAtomElec`, `AtomHasConjugatedBond`; new core `needsHs` and `assignChiralTypesFromMolParity` (`AssignAtomChiralTagsFromMolParity`, `NeedsHs`); `Get/SetAtomRLabel/Alias/Value/SupplementalSmilesLabel`, `GetNumPiElectrons`; exported `CHI_*`, `EXPLICIT/IMPLICIT`, `STEREO_*`, `COMPOSITE_*`, `ResonanceFlags` and pickle-option constants in `rdchem` | ✅ | probe cases identical to RDKit incl. error texts and upstream line numbers (argument-type errors are plain `TypeError`s instead of Boost `ArgumentError`s); `test_stereo_misc.py`, `tests/stereo_misc.rs`. |
| 12q | Bond wedging API moved out of the molfile writer into `graphmol::wedge_bonds`: `WedgeMolBonds` (incl. `BondWedgingParameters.wedgeTwoBondsIfPossible`), `WedgeBond`, `ReapplyMolBlockWedging`, clear/invert of the stored mol block wedging, `RemoveNonExplicit3DChirality`; `RDKitError::test_assert_at` | ✅ | probe cases and 51 stereo molecules (both parameter settings) identical to RDKit incl. error texts; `test_wedging.py`, `tests/wedge_bonds.rs`. |
| 12r | `copyMolSubset` as a full port in `graphmol::subset` (atoms+bonds or atom/bond paths, `SubsetOptions`/`SubsetInfo`, coordinates, stereo groups, double-bond stereo remapping, query copies, sanitize), used by `getMolFrags`; `pathToSubmol`; Python `CopyMolSubset` (4 overloads), `SubsetOptions`, `SubsetInfo`, `SubsetMethod`, `UIntUIntMap`, `PathToSubmol`; Boost-style unsigned-int sequence conversion | ✅ | 13,641 random atom/bond subsets of stereo-rich NCI molecules and ~60 probe cases identical to RDKit incl. error texts; `test_copy_mol_subset.py`, `tests/subset.rs`. |
| 12s | Enhanced stereo helpers in `chirality::stereo_groups`: `AddStereoAnnotations` (label templates), `SimplifyEnhancedStereo`, `FindMesoCenters` (via the canonical ranking, incl. ring stereo), `TranslateChiralFlagToStereoGroups`, `CleanupStereoGroups` binding; duplicate `assignStereoGroupIds` removed from the molfile writer | ✅ | ~120 probe cases (groups built through the API since CXSMILES is not ported) and 70 stereo molecules identical to RDKit; `test_stereo_groups.py`, `tests/stereo_groups.rs`. |
| 13a | CXSMILES reading (`smiles_parse::cx`, wired into `MolFromSmiles`/`MolFromSmarts`): coordinates, atom labels/values/props, radicals, enhanced stereo, wedge/wiggly bonds, `ctu`/`c`/`t`, dative/hydrogen/zero bonds, `u`/`rb`/`s` query features, special labels (`Q_e`, `star_e`, …), `_CXSMILES_Data`, strict/lenient failure; SMILES `~` any-bond. Substance groups (`Sg`/`SgD`), link nodes (`LN`) and variable attachment points (`m`) are not ported: such extensions fail to parse | ✅ | ~50 probe cases × strict/sanitize modes identical to RDKit except the unported sections; `test_cxsmiles_read.py` (passes on RDKit too), `tests/cxsmiles_read.rs`. |
| 13b | `MergeQueryHs` (`mol_ops::query_hs`: `[H]` neighbors → `!H<n>` queries, mapped/isotope options, recursive queries), `HasQueryHs`, `MolFromSmarts(mergeHs=True)` | ✅ | 36 SMARTS × 6 call variants identical to RDKit; `test_merge_query_hs.py` (passes on RDKit too), `tests/merge_query_hs.rs`. |
| 13c | `Chem.QueryAtom`/`QueryBond` (atoms and bonds of query molecules and `AtomFromSmarts`/`BondFromSmarts` objects are instances, `ExpandQuery`, `SetQuery`, query-aware `Match`), `rdkit.Chem.rdqueries` (all `...Equals/Less/GreaterQueryAtom`, flag and property query factories, `ReplaceAtomWithQueryAtom`, Boost-style `ArgumentError`s), `Mol.GetAtomsMatchingQuery` (`substruct::atoms_matching_query`, recursive queries) | ✅ | ~1000 factory/molecule combinations and the error texts identical to RDKit; `test_rdqueries.py` (passes on RDKit too), `tests/atoms_matching_query.rs`. |
| 13d | `AdjustQueryProperties` / `AdjustQueryParameters` / `AdjustQueryWhichFlags` (`graphmol::adjust_query`: degree, heavy-atom degree, ring count and ring/chain queries with the which-flags, generic atoms and bonds, dummy queries, aromatization, stereo-care bonds, MDL five-ring aromaticity, conjugated five-rings, single bonds next to aromatic atoms). The JSON parameter format and `AdjustQueryPropertiesWithGenericGroups` are not ported | ✅ | 21 SMILES + 6 SMARTS + 3 CTAB molecules × ~70 parameter sets identical to RDKit (including the errors); `test_adjust_query.py` (passes on RDKit too), `tests/adjust_query.rs`. |
| 13e | `SubstanceGroup` (`graphmol::substance_group`: atoms, parent atoms, bonds, brackets, CSTATEs, attach points, properties, type checks; kept consistent when atoms or bonds are removed, molecules are inserted, copied, split into fragments or subset, and by `removeHs`) and its Python API (`GetMolSubstanceGroups`, `GetMolSubstanceGroupWithIdx`, `CreateMolSubstanceGroup`, `CreateMolDataSubstanceGroup`, `AddMolSubstanceGroup`, `ClearMolSubstanceGroups`, `SubstanceGroup`, `SubstanceGroupCState`, `SubstanceGroupAttach`, `SubstanceGroup_VECT`). Pickling molecules with groups, molfile and CXSMILES input/output of groups follow in 13f–13h; `GetOwningMol()` returns the same object where RDKit returns a new wrapper | ✅ | ~110 probe lines (API, errors, 40 molecule operations) identical to RDKit; `test_substance_groups.py` (passes on RDKit too), `tests/substance_groups.rs`. |
| 13f | SubstanceGroups in mol blocks: V2000 `M  STY/SST/SLB/SCN/SDS/SAL/SBL/SPA/SMT/SDI/SBV/SDT/SDD/SCD/SED/SPL/SNC/SAP/SCL/SBT` lines and the V3000 `BEGIN SGROUP` block (`file_parsers::sgroup_parsing`: strict and lenient parsing, DEFAULT labels, quoted values, malformed SAP repair, `MolFileUnhandledFeatureException` skipping the molecule), and the structure-changing groups of `processSGroups` (MarvinSketch coordinate bonds, implicit Hs and SMARTS queries, RDKit ZBO/ZCH/HYD). The `M  ZBO/ZCH/HYD/MRV/APO/LIN` lines themselves and SGroup writing (13g) follow | ✅ | 162 SGroup-bearing mol files of the upstream tree × strict/lenient and ~60 hand-made good and malformed blocks identical to RDKit; `test_molfile_sgroups.py` (passes on RDKit too), `tests/molfile_sgroups.rs`. |
| 13g | SubstanceGroup output in mol blocks (`file_parsers::sgroup_writing`: V2000 `M  S..` lines and the V3000 `BEGIN SGROUP` block with 78-column wrapping, quoting, BRKTYP checks and the 200 character data limit), SMARTS-query atoms and zero-order bonds written as DAT groups (`moveAdditionalPropertiesToSGroups`), `M  ZBO/ZCH/HYD` output and input. Generic group queries (`_QueryAtomGenericLabel`) are still rejected | ✅ | 102 written blocks (V2000 and V3000 of 51 molecules, including every SGroup-bearing mol file of the upstream tree) identical to RDKit; `test_molfile_sgroup_write.py` (passes on RDKit too), `tests/molfile_sgroups.rs`. |
| 13h | SubstanceGroup pickling (`MolPickler` BEGINSGROUP: properties, atoms, parent atoms, bonds, brackets, cstates, attach points; byte-identical pickles, narrow and wide index format). CXSMILES `Sg`/`LN`/`m` sections still follow | ✅ | `ToBinary` of SRU/SUP/DAT/300-atom GEN groups md5-identical to RDKit and round-tripping through `Chem.Mol`/`pickle`; `test_pickle_sgroups.py` (passes on RDKit too), `tests/substance_groups.rs`. |
| 10a | Wildman-Crippen logP/MR: `CalcCrippenDescriptors` (includeHs, cached sums), `_CalcCrippenContribs` (atom types and labels, cache semantics), `rdkit.Chem.Crippen` (verbatim `MolLogP`/`MolMR`) | ✅ | differential `crippen` 15000/15000 (ZINC 10k + chiral 5k); upstream `UnitTestCrippen` 6/6 |
| 10a | `Transform3D` (identity, rotations, point transform) for addHs coordinates; the rest of Transform3D waits for its consumers | 🟡 | covered by `addhs3d` |
| 10b | Lipinski counts: HBA/HBD (Lipinski and SMARTS variants), rotatable bonds (all four definitions; `Default` is Strict as in the released build), heteroatoms, amide bonds, ring families, FractionCSP3, spiro/bridgehead atoms, stereo-center counts | ✅ | differential `lipinski` 15300/15300 identical (NCI 300 + ZINC 10k + chiral 5k); upstream bool/enum argument handling reproduced |
| 10b | `CalcPhi` (Hall-Kier alpha, kappa helpers of ConnectivityDescriptors) and verbatim `Chem/Lipinski.py` | ✅ | identical to RDKit on 14 molecules; upstream `UnitTestLipinski` 2/3 (`testLipinskiLong` is skipped, and fails on the reference build too: its regression data predates the strict default) |
| 10e | `_CalcMolWt`, `CalcExactMolWt`, `CalcMolFormula` (Hill order, isotopes, D/T), `CalcNumHeavyAtoms`, `CalcNumAtoms` (`MolOps::getAvgMolWt`/`getExactMolWt`/`getMolFormula`) | ✅ | differential `molprops` 15312/15312 identical (also with explicit Hs, isotopes, charges); bool arguments take truthy values |
| 10f | `CalcLabuteASA`, `_CalcLabuteASAContribs`, `CalcTPSA` (+ `includeSandP`), `_CalcTPSAContribs`, `SlogP_VSA_`, `SMR_VSA_`, `CustomProp_VSA_` (`MolSurf.cpp`; `std::upper_bound` search reproduced for user bins; per-mol caches like upstream) | ✅ | differential `molsurf` 15312/15312 bit-identical (NCI, ZINC, chiral set) |
| 10g | Gasteiger charges (`rdPartialCharges.ComputeGasteigerCharges`, parameter tables, conjugated charge sharing), `PEOE_VSA_`, verbatim `Chem/MolSurf.py` | ✅ | differential `gasteiger` 15312/15312 identical (string and double properties); upstream `UnitTestSurf` 9/13 + 4 long tests that also pass with `doLong` (one of them fails identically on the reference: bad SMILES in its data), `testPartialCharges` 6/7 |
| 10h | Connectivity descriptors: `CalcChi{0..4}{v,n}`, `CalcChiNv`/`CalcChiNn`, `CalcHallKierAlpha` (with `atomContribs`), `CalcKappa1/2/3` (+ `CalcPhi` moved to its own binding module); the HK deltas / nVals caches follow upstream (including its unsigned `nOuter - nHs`) | ✅ | differential `connectivity` 15312/15312 identical (plain and explicit-H forms) |
| 11a | `Descriptors` module prerequisites: property registry (`rdMolDescriptors.Properties`, `PythonPropertyFunctor`, `PropertyRangeQuery`, `MakePropertyRangeQuery`; 43 built-in functors in upstream order, Python-defined functors via `RegisterProperty`), `rdkit.ML.InfoTheory.rdInfoTheory` (`InfoEntropy`, `InfoGain`, `ChiSquare`), `_vect<string>` type, `FragmentDescriptors.csv` data file; verbatim `Descriptors.py`, `Fragments.py`, `GraphDescriptors.py`, `Graphs.py`, `QED.py`, `SpacialScore.py`, `Descriptors3D.py`, `EState/*`, `ChemUtils` | 🟡 | `from rdkit.Chem import Descriptors` works with 209 of 217 descriptors listed; 178 of them evaluate (rest need `GetDistanceMatrix`, `DeleteSubstructs`, `FindPotentialStereoBonds`); `BCUT2D_*` (8) missing; `InfoBitRanker`/`BitCorrMatGenerator` not ported |
| 11b | `Chem.GetDistanceMatrix` (bond orders, atom weights; Floyd-Warshall with upstream's operation order), `GetAdjacencyMatrix`, `Get3DDistanceMatrix` with upstream's computed-property caching (`DistanceMatrix`, `_Paths`, `force`, `prefix`); `DeleteSubstructs`; `FindPotentialStereoBonds`; `Get/SetUseLegacyStereoPerception`, `Get/SetAllowNontetrahedralChirality`; `_vect<string>` for string-vector props; pure-Python tail of `Chem/__init__.py` (`FindMolChiralCenters`, `CanonSmiles`, ...) | ✅ | differential `matrices`, `transforms` and `descriptors` (all 200+ `Descriptors._descList` values on 10312 molecules) identical. `SPS` waits for the new stereo perception (`FindPotentialStereo`), `BCUT2D_*` for BCUT2D |
| 11c | BCUT2D (`rdMolDescriptors.BCUT2D` with all overloads, `_BCUT2D_version`, `Descriptors.BCUT2D_*`); `numerics::symmetric_eigenvalues` (Householder + implicit QL replaces Eigen's solver; agrees to ~1e-14, not bit-identical); `removeAllHs` in core | ✅ | differential `descriptors` (BCUT values to 9 decimals) 10312/10312; `Descriptors._descList` now has all 217 entries of RDKit |
| 12a | `SmilesMolSupplier` (lazy + random access over delimited text: title line, name/SMILES columns, delimiter set, property columns, comment/blank lines, CRLF, `std::istream` state emulation incl. upstream's quirks), `SmilesMolSupplierFromText`; Python iteration/indexing/len/SetData/GetItemText/context manager | ✅ | 22 scenario transcripts (basic, no title, tabs, name column -1, columns swapped, unsanitized, CRLF, no final newline, empty, comment-only, extra columns, short lines, files, bad/empty files, `GetItemText` before `next`) identical to RDKit incl. error messages |
| 10a | rest of `rdMolDescriptors` (MQNs, BCUT, AUTOCORR2D, ...) and the `Descriptors` module (needs EState, QED, SpacialScore, GraphDescriptors) | ⬜ | next |
| 9 | LayeredFingerprint, PatternFingerprint, JSON (de)serialization, GetFPs bulk helpers | ⬜ | |
| 7–13 | see PLAN.md | ⬜ | |

## Python package files

| File | Kind |
|------|------|
| `rdkit/__init__.py`, `RDLogger.py`, `RDRandom.py`, `RDConfig.py` | verbatim upstream (porting/copies.toml) |
| `rdkit/RDPaths.py` | rdoxide's version of the build-generated file (package-relative paths, like the PyPI wheels) |
| `rdkit/rdBase.py`, `rdkit/Chem/rdchem.py`, `rdkit/Geometry/rdGeometry.py` | shims → `rdkit._rdoxide` |
| `rdkit/DataStructs/{__init__,TopNContainer,BitUtils,VectCollection,LazySignature,BitEnsemble}.py` | verbatim upstream |
| `rdkit/DataStructs/cDataStructs.py` | shim → `rdkit._rdoxide` |
| `rdkit/Chem/MACCSkeys.py` | verbatim upstream |
| `rdkit/Chem/rdFingerprintGenerator.py`, `rdkit/Chem/rdMolDescriptors.py` | shims → `rdkit._rdoxide` |
| `rdkit/Chem/__init__.py`, `rdkit/Geometry/__init__.py` | provisional (until rdmolfiles/rdmolops/DataStructs exist) |

## Iteration log

* 2026-10-07 — Stage 0 + 1 + PeriodicTable.
* 2026-10-07 — Stage 2: core graph + rdchem bindings. Next: stage 3/4 critical path —
  SMILES parser (no sanitize) → FindRings (SSSR/symm) → Kekulize/Aromaticity/ConjugHybrid →
  sanitizeMol → SMILES writer + canonical ranking. Then rough_test.py starts importing.
* 2026-10-07 — Stage 3a–c: rings, sanitizeMol, canonical ranking.
* 2026-10-07 — Stage 3d: removeHs + legacy stereo perception; default `MolFromSmiles`
  is now identical to RDKit. Next: SMILES writer (SmilesWrite.cpp, Canon.cpp
  canonicalizeFragment, MolToSmiles/MolFragmentToSmiles) → rough_test.py imports.
* 2026-10-08 — Stage 4b: SMILES writer identical to RDKit. Next: SMARTS + query atoms/bonds
  (needed by rough_test.py imports and substructure search), then CXSMILES.
* 2026-10-08 — Stage 4c (part 1): query model + SMARTS parser identical to RDKit.
* 2026-10-08 — Stage 4c (part 2): SMARTS writer identical to RDKit. Next: substructure matching
  (stage 5: vf2, SubstructMatch, HasSubstructMatch/GetSubstructMatch(es), recursive queries),
  then mergeQueryHs + QueryAtom classes.
* 2026-10-08 — Stage 5 (core): substructure matching identical to RDKit (match order included).
  Next: what rough_test.py needs to import (Chem/__init__ pieces), then MolPickler and file
  formats (stage 6).
* 2026-10-08 — RDConfig/RDPaths added. Remaining upstream-test blockers by count: DataStructs
  (~19 tests), MolBlock/SDMolSupplier/MolToMolBlock (~22), AllChem (11), pickling (3).
  Next: stage 6 (MolPickler, then mol file reader/writer, SDMolSupplier).
* 2026-10-08 — Stage 6a: MolPickler byte-identical to RDKit; Python pickling. Ratchet 53 -> 56.
  Next: MolBlock reader/writer + SDMolSupplier (stage 6b).
* 2026-10-08 — Stage 6b (part 1): V2000 MolBlock reader + stereo from coordinates identical to
  RDKit (`molblock` differential), then the V3000 CTAB. Ratchet 56 -> 58. Next: SDMolSupplier,
  MolToMolBlock, SGroups.
* 2026-10-08 — Stage 6c: SDMolSupplier identical to RDKit (`sdf` differential). Ratchet 58 -> 60.
  Next: MolToMolBlock/SDWriter, then SGroups.
* 2026-10-08 — Stage 6d: MolBlock writer byte-identical to RDKit (`molwrite` differential); libc++
  `std::sort` port (`rdgeneral::cxx_sort`, Apache-2.0 WITH LLVM-exception) for RDKit's unstable
  sorts. Ratchet 60 -> 64. Next: RDDepict (2D coordinates), SGroups, SDWriter.

## Differential testing

`tools/differential/compare.py MODE` runs the same dump script under the
official RDKit of the ported version and under rdoxide (`.venv`) on the
upstream SMILES corpora (NCI first_5K, canonSmiles.long) and reports
differences. The reference interpreter comes from conda-forge
(`micromamba create -p REF -c conda-forge python=3.13 rdkit=2026.09.1`) and is
selected with `--ref-python` or `RDOXIDE_REF_PYTHON`. (PyPI only has
2026.03.x, which predates the RDL-based SSSR, so it is not a valid oracle for
ring order.)

Current results: `parse`, `rings`, `sanitize`, `ranks`, `molfromsmiles`, `smiles`: 10312/10312 each;
`smarts`: 1426/1426 (RDKit's SMARTS libraries + `tests/data/smarts_edge_cases.smarts`);
`substruct`: 10312/10312 + 2418/2418 stereo-rich (queries in `tests/data/substruct_queries.smarts`).
`molblock`: corpora built with `tools/differential/make_molblock_corpus.py` under the reference
interpreter (NCI 2D with wedges/explicit Hs, embedded 3D, upstream FileParsers test_data).
`molfromsmiles` (atoms, bonds, stereo props `_CIPCode`/`_CIPRank`/`_ringStereoAtoms`/
`_ChiralityPossible`) also matches on 35140 molecules from the stereo-rich upstream
corpora (zinc, ChEMBL, idorsia enumeration, MMFF) and on `tests/data/stereo_edge_cases.smi`.
New Python test files are also run against the official wheel to validate expectations.

## Known deviations (to revisit)

* Floating-point results follow the native RDKit build of the platform: on aarch64 rdoxide fuses
  multiply-adds exactly where clang does (`rdgeneral::fp::fmuladd`), on x86-64 it does not. The
  official arm64 and x86-64 builds themselves differ in the last bits of 2D coordinates (and so in
  a few depictions); the `depict` differential uses the arm64 build. Windows' libm can flip the
  sign of rounded zeros (`-0.0000`) in MolBlocks.

* Where RDKit uses the unstable `std::sort` (e.g. picking wedge bonds), rdoxide reproduces
  libc++'s element order (the macOS/conda-forge reference). RDKit builds that use libstdc++
  (Linux GCC) can order ties differently.

* Atom/Bond objects returned from a Mol hold a reference to the Python Mol plus
  an index (upstream holds a C++ pointer); after removing atoms the index may
  point at a different atom (upstream: dangling pointer).
* `Atom.IsInRing()` needs stage 3 (`findSSSR`).
* `MolFragmentToSmiles(mol, params, atomsToUse)` with `params` passed positionally works
  in rdoxide; upstream's Boost overload resolution raises AttributeError for it (keyword
  form works in both).
* `Chem.SetAllowNontetrahedralChirality` & co. keep the value in-process;
  upstream also writes it to the process environment.
* 2026-10-08 — Stage 8a: DataStructs bit vectors + similarity. Next: SparseIntVect/DiscreteValueVect, then fingerprints (Morgan, RDKit, atom pairs) or SGroups.
* 2026-10-08 — Stage 8b: SparseIntVect. Next: DiscreteValueVect/RealValueVect, then fingerprints.
* 2026-10-08 — Stage 8c: DiscreteValueVect, RealValueVect. CI: macOS jobs commented out at the user's request. Next: fingerprints (Morgan, RDKit, atom pairs, topological torsions) or SGroups.
* 2026-10-08 — Stage 9a: fingerprint generator framework + Morgan. Next: RDKit FP / atom pairs / torsions generators.
* 2026-10-08 — Stage 9b: atom-pair generator. Next: Subgraphs (path/subgraph enumeration) → topological torsion + RDKit FP generators.
* 2026-10-08 — Stage 9c: Subgraphs + rdBase vector types. Next: topological torsion + RDKit FP generators.
* 2026-10-08 — Stage 9d: topological-torsion generator. CI blocked by GitHub billing (jobs not started); PRs #17/#18 wait. Next: RDKit FP generator.
* 2026-10-08 — Stage 9e: RDKit FP generator (CI back; #17-#19 merged, merged branches deleted). Next: legacy fingerprint functions (rdMolDescriptors / Chem.RDKFingerprint, AllChem.GetMorganFingerprint*), MACCS, pattern FP.
* 2026-10-08 — Stage 9f: legacy RDKFingerprint. Next: rdMolDescriptors legacy Morgan / atom pairs / torsions, MACCS.
* 2026-10-08 — Stage 9g: rdMolDescriptors legacy fingerprints. Next: MACCS keys, then descriptors (rdMolDescriptors / Descriptors).
* 2026-10-08 — Stage 9h: MACCS keys. Next: descriptors (rdMolDescriptors: CalcNumAtoms, formula, MW, ring counts, Lipinski counts, TPSA, logP/MR (Crippen), ...).
* 2026-10-08 — Stage 10a (part 1): Crippen logP/MR (`CalcCrippenDescriptors`, `_CalcCrippenContribs`, `rdkit.Chem.Crippen`) and addHs (with coordinates) identical to RDKit; ratchet 169. Next: Lipinski and the other descriptors (MolWt, TPSA, counts), then the rest of `rdMolDescriptors`.
* 2026-10-08 — Stage 10b: Lipinski counts identical to RDKit (rotatable-bond definitions with the released strict default, bool/enum argument handling, spiro/bridgehead lists). Next: CalcPhi and verbatim Lipinski.py, then MolWt/TPSA/Labute ASA.
* 2026-10-09 — Stage 10c: CalcPhi and verbatim Lipinski.py; Crippen/AddHs-style bool arguments accept any truthy value like Boost.Python; multi-core batch fingerprints (`numThreads`); crate license audit in CI. Ratchet 172. Next: MolWt/ExactMolWt, TPSA, Labute ASA.
* 2026-10-09 — Stage 10e: MolWt/ExactMolWt/formula/atom counts identical to RDKit (`molprops` differential). Next: TPSA, Labute ASA (MolSurf), then the `Descriptors` module.
* 2026-10-09 — Stage 10f: Labute ASA, TPSA and the SlogP/SMR/custom VSA binnings identical to RDKit (`molsurf` differential). Next: Gasteiger charges, `PEOE_VSA_`, verbatim `MolSurf.py`.
* 2026-10-09 — Stage 10g: Gasteiger charges, `PEOE_VSA_` and verbatim MolSurf.py identical to RDKit (`gasteiger` differential). Next: the rest of `rdMolDescriptors` (Chi/Kappa, MQNs, ...) and the `Descriptors` module.
* 2026-10-09 — Stage 10h: Chi/Kappa/Hall-Kier connectivity descriptors identical to RDKit (`connectivity` differential; `res += a*b` of Chi1 is fused on arm64). Next: the `Descriptors` module prerequisites (EState, GraphDescriptors, ...).
* 2026-10-09 — Stage 11a: `Descriptors` module imports (property registry, InfoTheory, Fragments data, verbatim pure-Python modules). Ratchet 231. Next: `GetDistanceMatrix`, `DeleteSubstructs`, `FindPotentialStereoBonds` (EState, GraphDescriptors, QED, SPS), then BCUT2D.
* 2026-10-09 — Stage 11b: distance/adjacency matrices, DeleteSubstructs, FindPotentialStereoBonds; every `Descriptors` value but SPS and BCUT2D identical to RDKit. Ratchet 253. Next: BCUT2D, AUTOCORR2D, new stereo perception (SPS), then the remaining `rdMolDescriptors` functions (MQNs, ...).
* 2026-10-09 — Stage 11c: BCUT2D with an own symmetric eigenvalue solver; provisional `Chem/AllChem.py` (imports what is ported). Ratchet 292. Next: SmilesMolSupplier, AUTOCORR2D/MQNs, new stereo perception (SPS), 3D descriptors (PMI, PBF, ...).
* 2026-10-09 — Stage 12a: SmilesMolSupplier. Ratchet 298. Next: SmilesWriter, ForwardSDMolSupplier, then AUTOCORR2D/MQNs and the 3D descriptors.
* 2026-10-09 — Stage 12b: uniform grids, PointND, dihedrals. Ratchet 316. Next: SmilesWriter, ForwardSDMolSupplier, AUTOCORR2D/MQNs, 3D descriptors.
* 2026-10-09 — Stage 12c: SmilesWriter, shared writer sink. Next: ForwardSDMolSupplier, AUTOCORR2D/MQNs, 3D descriptors.
* 2026-10-09 — Stage 12d: ForwardSDMolSupplier. Next: AUTOCORR2D/MQNs, 3D descriptors, FPBReader.
* 2026-10-09 — Stage 12e: MQNs, AUTOCORR2D (LICENSES/RDKit-Descriptors3D-Godin-BSD-3-Clause.txt). Next: 3D descriptors (PMI, PBF, WHIM, RDF, ...), FPBReader, new stereo perception (SPS).
* 2026-10-09 — Stage 12f: rdMolTransforms, Eigen-compatible 3x3 eigensolver (MPL-2.0, LICENSES/Eigen-MPL-2.0.txt). Next: PMI/PBF/other 3D descriptors on top of it, FPBReader, new stereo perception (SPS).
* 2026-10-09 — Stage 12g: PMI/PBF family (LICENSES/RDKit-PBF-InstituteOfCancerResearch-BSD-3-Clause.txt). Next: WHIM/RDF/MORSE/GETAWAY/AUTOCORR3D on the shared weights, FPBReader, new stereo perception (SPS).
* 2026-10-09 — Stage 12h: AUTOCORR3D/RDF/MORSE (shared `radial` pair-sum kernel). Next: WHIM/GETAWAY/CoulombMat/EEM/USR, FPBReader, new stereo perception (SPS).
* 2026-10-09 — Stage 12i: CoulombMat + USR/USRCAT (LICENSES README now lists USRDescriptor under the Novartis file). Next: WHIM/GETAWAY/EEM, FPBReader, new stereo perception (SPS).
* 2026-10-09 — Stage 12j: WHIM + EEM (new numerics: `full_piv_lu.rs`, `jacobi_svd.rs`, both MPL-2.0 translations of Eigen 5.0.1, already listed in LICENSES). Next: GETAWAY, FPBReader, new stereo perception (SPS).
* 2026-10-09 — Stage 12k: GETAWAY completes the Descriptors3D set. Next: FPBReader, new stereo perception (SPS), SubstanceGroups, rdqueries.
* 2026-10-09 — Stage 12l: mol-structure rdmolops (frags, combine, renumber, formal charge, shortest path, remove stereo). Next: CXSMILES, FPBReader, new stereo perception (SPS), SubstanceGroups, rdqueries.
* 2026-10-09 — Stage 12m: ChemTransforms (ReplaceSubstructs/Sidechains/Core, MurckoDecompose). Next: FragmentOnBonds (MolFragmenter), molzip, CXSMILES, FPBReader, new stereo perception (SPS).
* 2026-10-09 — Stage 12n: MolFragmenter (FragmentOnBonds/SomeBonds/BRICSBonds). Next: molzip, CXSMILES, FPBReader, new stereo perception (SPS).
* 2026-10-09 — Stage 12o: molzip. Next: CXSMILES, FPBReader, new stereo perception (SPS), SubstanceGroups, rdqueries.
* 2026-10-09 — Stage 12p: stereo/atom helper bindings and rdchem enum exports.
* 2026-10-09 — Stage 12q: bond wedging API (WedgeMolBonds, WedgeBond, ReapplyMolBlockWedging).
* 2026-10-09 — Stage 12r: copyMolSubset / PathToSubmol.
* 2026-10-09 — Stage 12s: enhanced stereo helpers. Next: CXSMILES, FPBReader, new stereo perception (SPS), SubstanceGroups, rdqueries.
* 2026-10-09 — Stage 13a: CXSMILES reading. Next: CXSMILES writing + SubstanceGroups, FPBReader, new stereo perception (SPS), rdqueries.
* 2026-10-09 — Stage 13b: MergeQueryHs/HasQueryHs. Next: QueryAtom/QueryBond/rdqueries, CXSMILES writing + SubstanceGroups, FPBReader, new stereo perception (SPS).
* 2026-10-09 — Stage 13c: QueryAtom/QueryBond classes and rdqueries. Next: CXSMILES writing + SubstanceGroups, FPBReader, new stereo perception (SPS), AdjustQueryProperties.
* 2026-10-09 — Stage 13d: AdjustQueryProperties. Next: CXSMILES writing + SubstanceGroups, FPBReader, new stereo perception (SPS), GenericGroups.
* 2026-10-09 — Stage 13e: SubstanceGroup core + Python API. Next: SGroups in molfiles (13f read, 13g write), pickling and CXSMILES (13h), FPBReader, new stereo perception (SPS), GenericGroups.
* 2026-10-09 — Stage 13f: SGroups in mol blocks (V2000/V3000 reading).
* 2026-10-09 — Stage 13g: SGroup writing. Next: pickling and CXSMILES (13h), M  MRV/APO/LIN lines, FPBReader, SPS, GenericGroups.
* 2026-10-09 — Stage 13h: SubstanceGroup pickling. Next: CXSMILES SGroup sections, M  MRV/APO/LIN lines, FPBReader, SPS, GenericGroups.
