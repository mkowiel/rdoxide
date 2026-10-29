#!/usr/bin/env python3
"""Record extended Hueckel results of the official RDKit wheel (rdEHTTools) for tests/data/yaehmop_eht.json.

usage: .venv-ref/bin/python tools/make_yaehmop_corpus.py"""
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdEHTTools

SMILES = """CCO c1ccccc1 CC(=O)O CCN(CC)CC c1ccncc1 C1CCCCC1 CS(=O)(=O)C CC(C)(C)Cl BrCCBr CP(C)C O=C=O N#CC
FC(F)(F)F c1ccc2ccccc2c1 C[Si](C)(C)C OB(O)O C1=CC=CS1 c1ccoc1 CC[N+](C)(C)C [O-]C(=O)C CCS CCI
O=[N+]([O-])c1ccccc1 CC(C)C(=O)NC1CC1 C#C C=C=C F[B-](F)(F)F Fc1ccccc1 O=C1CCCN1 O OO
CC(C)Cc1ccc(cc1)C(C)C(=O)O CN1C=NC2=C1C(=O)N(C(=O)N2C)C [Li]C [Mg](C)C C[Al](C)C [Zn](C)C
Cl[Ti](Cl)(Cl)Cl [Fe](C)(C)(C)C""".split()

cases = []
for smi in SMILES:
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    if AllChem.EmbedMolecule(m, randomSeed=7) != 0:
        if AllChem.EmbedMolecule(m, randomSeed=7, useRandomCoords=True) != 0:
            continue
    block = Chem.MolToMolBlock(m)
    m = Chem.MolFromMolBlock(block, removeHs=False, sanitize=False)
    ok, r = rdEHTTools.RunMol(m, keepOverlapAndHamiltonianMatrices=True)
    cases.append(
        dict(
            smiles=smi,
            molblock=block,
            numOrbitals=r.numOrbitals,
            numElectrons=r.numElectrons,
            fermiEnergy=r.fermiEnergy,
            totalEnergy=r.totalEnergy,
            charges=r.GetAtomicCharges().tolist(),
            energies=r.GetOrbitalEnergies().tolist(),
            rop=r.GetReducedOverlapPopulationMatrix().tolist(),
            rchg=r.GetReducedChargeMatrix().ravel().tolist(),
            hamiltonian_sum=float(r.GetHamiltonian().sum()),
            overlap_sum=float(r.GetOverlapMatrix().sum()),
        )
    )
Path("tests/data/yaehmop_eht.json").write_text(json.dumps(cases, separators=(",", ":")))
print(len(cases), "cases")
