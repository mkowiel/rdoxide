# pubchem_shape (`rdkit.Chem.rdShapeAlign`) port

RDKit's `External/pubchem_shape` wraps NCBI's pubchem-align3d (commit daefab3, public domain, ~2000 lines of float code;
`exp` is `expApproxDouble`, the lookup tables are unused) in `PubChemShape.cpp` (766 lines) and `Wrap/cshapealign.cpp`.
It is independent of the newer `rdGaussianShape` already ported (stage 45): the scores differ (upstream `test_rdshapealign.py`
expects 0.773/0.305 for `test1.sdf`; `rdGaussianShape.AlignMol` gives other numbers), so it needs its own port.

| step | content | oracle |
|---|---|---|
| S1 `external::pubchem_shape::functions` | `shape_functions1/2.cpp`: Gaussian alphas, overlaps and analytic gradients, quaternion/matrix algebra, steric frame (eigen decomposition), colour type sets | C++ driver built from the NCBI sources, hex floats |
| S2 `neighbor` | `shape_neighbor.cpp`: `Neighbor_Conformers` (starting points, two-phase BFGS-like optimization) | same driver |
| S3 `ShapeInput`, `PrepareConformer`, `AlignShape`, `AlignMolecule`, `ScoreShape`, `TransformConformer` | `PubChemShape.cpp` | official wheel's `rdShapeAlign` on `test1.sdf` and more pairs |
| S4 Python `rdoxide.Chem.rdShapeAlign` | `Wrap/cshapealign.cpp`, `ShapeInputOptions`, `FloatVector` | upstream `test_rdshapealign.py` in the ratchet |

Licence file to add with S1: `LICENSES/PubChemAlign3D-PublicDomain.txt` (the NCBI notice).
