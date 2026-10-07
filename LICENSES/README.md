# Third-party licenses

| File | Applies to |
|------|-----------|
| `RDKit-BSD-3-Clause.txt` | All code translated from RDKit, verbatim copies of RDKit Python modules, RDKit-derived data files (`crates/rdoxide/data/`) and ported tests. Source: `license.txt` in the RDKit repository. |
| `RingDecomposerLib-BSD-3-Clause.txt` | Code translated from RingDecomposerLib v1.1.3_rdkit (University of Hamburg, ZBH; Andresen, Flachsenberg, Rarey), used by RDKit for ring families / SSSR: `crates/rdoxide/src/external/ring_decomposer_lib/`. |
| `RDKit-ChemReactions-Novartis-BSD-3-Clause.txt` | Code translated from `Code/GraphMol/ChemReactions` (Novartis). Source: `Code/GraphMol/ChemReactions/license.txt`. |

When porting a component that carries its own license (anything under
upstream `External/`, fonts in `Data/Fonts`, ...), copy its license file here,
add a row to this table and to `NOTICE`.
