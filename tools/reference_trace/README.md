# Reference traces

Small C++ programs linked against the official RDKit libraries of the ported
version, printing intermediate results that the Python API does not expose.
They are used to find the first point where rdoxide diverges.

Setup (conda-forge reference environment, see docs/STATUS.md):

```bash
micromamba install -p $REF -c conda-forge librdkit-dev=2026.09.1 libboost-headers
clang++ -std=c++20 -O2 -I$REF/include/rdkit -I$REF/include depict_trace.cpp \
  -L$REF/lib -Wl,-rpath,$REF/lib -lRDKitDepictor -lRDKitSmilesParse \
  -lRDKitGraphMol -lRDKitRDGeneral -lRDKitRDGeometryLib -o depict_trace
```

`depict_trace SMILES` prints the 2D coordinates after each stage of
`compute2DCoords`; `RDOXIDE_DEPICT_TRACE=1` makes rdoxide print the same
lines (on stderr) from `rdDepictor.Compute2DCoords(mol, canonOrient=False)`.
