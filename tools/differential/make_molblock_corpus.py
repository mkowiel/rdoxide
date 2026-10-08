"""Build Mol block corpora for the `molblock` differential mode.

Run with an interpreter that has the official RDKit (the writer is used to
produce the blocks):

    $RDOXIDE_REF_PYTHON -I tools/differential/make_molblock_corpus.py OUTDIR

Writes OUTDIR/molblocks_2d.sdf (NCI first 5K with 2D coordinates and
wedges, with and without explicit Hs), OUTDIR/molblocks_3d.sdf (embedded
conformers) and OUTDIR/molblocks_upstream.sdf (the upstream test-data Mol
files and SD records).
"""
import glob
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDepictor

RDLogger.DisableLog("rdApp.*")
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")


def smiles(path="upstream/rdkit/Data/NCI/first_5K.smi"):
    with open(os.path.join(ROOT, path)) as fh:
        for line in fh:
            parts = line.split()
            if parts:
                yield parts[0]


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "molblocks_2d.sdf"), "w") as fh:
        for i, smi in enumerate(smiles()):
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            if i % 5 == 0:
                m = Chem.AddHs(m)
            rdDepictor.Compute2DCoords(m)
            try:
                fh.write(Chem.MolToMolBlock(m, kekulize=(i % 7 != 0)) + "$$$$\n")
            except Exception:  # noqa: BLE001
                pass
    with open(os.path.join(out, "molblocks_3d.sdf"), "w") as fh:
        for i, smi in enumerate(smiles()):
            if i % 10:
                continue
            m = Chem.MolFromSmiles(smi)
            if m is None or m.GetNumAtoms() > 40:
                continue
            mh = Chem.AddHs(m)
            if AllChem.EmbedMolecule(mh, randomSeed=0xF00D) != 0:
                continue
            if i % 20 == 0:
                mh = Chem.RemoveHs(mh)
            fh.write(Chem.MolToMolBlock(mh) + "$$$$\n")
    with open(os.path.join(out, "molblocks_v3000.sdf"), "w") as fh:
        stereo = list(smiles("tests/data/stereo_edge_cases.smi"))
        for i, smi in enumerate(list(smiles()) + stereo * 3):
            if i % 3 and i < 5000:
                continue
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            if i % 2:
                m = Chem.AddHs(m)
            rdDepictor.Compute2DCoords(m)
            centers = [a.GetIdx() for a in m.GetAtoms()
                       if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]
            if centers and (i % 4 == 0 or i >= 5000):
                rw = Chem.RWMol(m)
                kind = (Chem.StereoGroupType.STEREO_OR, Chem.StereoGroupType.STEREO_AND,
                        Chem.StereoGroupType.STEREO_ABSOLUTE)[i % 3]
                rw.SetStereoGroups([Chem.CreateStereoGroup(kind, rw, centers[:2])])
                m = rw.GetMol()
            try:
                fh.write(Chem.MolToV3KMolBlock(m) + "$$$$\n")
            except Exception:  # noqa: BLE001
                pass
    data = os.path.join(ROOT, "upstream/rdkit/Code/GraphMol/FileParsers/test_data")
    with open(os.path.join(out, "molblocks_upstream.sdf"), "w") as fh:
        for f in sorted(glob.glob(os.path.join(data, "*.mol"))):
            txt = open(f, "rb").read().decode("utf-8", "replace").replace("\r\n", "\n")
            if "$$$$" in txt:
                continue
            fh.write(txt.rstrip("\n") + "\n$$$$\n")
        for f in sorted(glob.glob(os.path.join(data, "*.sdf"))):
            txt = open(f, "rb").read().decode("utf-8", "replace").replace("\r\n", "\n")
            for rec in txt.split("$$$$\n")[:50]:
                if rec.strip():
                    fh.write(rec.split("M  END")[0] + "M  END\n$$$$\n" if "M  END" in rec else rec + "$$$$\n")


if __name__ == "__main__":
    main()
