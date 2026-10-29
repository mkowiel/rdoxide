# FreeSASA port

Goal: `rdoxide.Chem.rdFreeSASA` (RDKit's `External/FreeSASA`): `CalcSASA`, `classifyAtoms`, `SASAOpts`, the polar/apolar
queries. RDKit builds FreeSASA 2.0.3 (MIT, Simon Mitternacht); `LICENSES/FreeSASA-MIT.txt`.

| step | content | verification |
|---|---|---|
| F1 `external::freesasa::{nb, sr, lr}` | cell-list neighbor lists, Shrake & Rupley, Lee & Richards, `freesasa_calc_coord` | 82 cases (20 drug-like molecules with hydrogens, a 700-atom protein; both algorithms, probes 1.4 and 0.8) equal the C library bit for bit on macOS (`tests/data/freesasa_sasa.txt`, `tools/freesasa_drv/sasa.c`, `tools/make_freesasa_corpus.py`); perturbing each of the four files fails the test |
| F2 classifiers | Protor, NACCESS, OONS radii and polar/apolar classes by residue and atom name | planned |
| F3 wrapper | `calcSASA`, `classifyAtoms` (atom/molecule properties `SASA`, `SASAClass`, `SASAClassName`), Python `rdFreeSASA`, upstream `testFreeSASA.py` | planned |

Deliberate difference: the C code reads an uninitialized neighbor for an atom without neighbors (Shrake & Rupley); here such
an atom is fully exposed. Threads are not used: the areas of the atoms do not depend on the thread count.
