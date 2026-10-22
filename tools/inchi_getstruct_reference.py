"""Generates tests/data/inchi_getstruct.json from the InChI C library:
`GetStructFromINCHI` (atoms with bonds, charges, isotopes, hydrogen counts,
and the 0D stereo) for the InChI strings of tests/data/inchi_corpus.jsonl and
of the PubChem "hard set" of upstream's UnitTestInchi.

Compiles the whole InChI library (INCHI_BASE/src plus the libinchi API files).

    python tools/inchi_getstruct_reference.py
"""
import json
import pickle
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
B = ROOT / "upstream/external/InChI"
SRC = B / "INCHI_BASE/src"
API = B / "INCHI_API/libinchi/src"
HARD = ROOT / "upstream/rdkit/rdkit/Chem/test_data/pubchem-hard-set.inchi"

HARNESS = r"""
#include <stdio.h>
#include <string.h>
#include "inchi_api.h"
int main(void){ static char line[65536]; char opt[1]="";
 while(fgets(line,sizeof line,stdin)){ line[strcspn(line,"\n")]=0; inchi_InputINCHI in; inchi_OutputStruct out; in.szInChI=line; in.szOptions=opt;
  memset(&out,0,sizeof out);
  int r=GetStructFromINCHI(&in,&out); printf("R %d %d %d\n",r,out.num_atoms,out.num_stereo0D);
  if(r<2){ for(int i=0;i<out.num_atoms;i++){ inchi_Atom*a=&out.atom[i]; printf("A %s %d %d %d %d %d %d %d %d",a->elname,a->charge,a->radical,a->isotopic_mass,a->num_iso_H[0],a->num_iso_H[1],a->num_iso_H[2],a->num_iso_H[3],a->num_bonds);
    for(int b=0;b<a->num_bonds;b++) printf(" %d,%d,%d",a->neighbor[b],a->bond_type[b],a->bond_stereo[b]); printf("\n"); }
   for(int i=0;i<out.num_stereo0D;i++){ inchi_Stereo0D*s=&out.stereo0D[i]; printf("S %d %d %d %d %d %d %d\n",s->neighbor[0],s->neighbor[1],s->neighbor[2],s->neighbor[3],s->central_atom,s->type,s->parity);}
   FreeStructFromINCHI(&out);} printf("E\n"); }
 return 0;}
"""


def inchis():
    out = set()
    for line in open(ROOT / "tests/data/inchi_corpus.jsonl"):
        v = json.loads(line)
        if v["kind"] == "mol" and v["inchi"]:
            out.add(v["inchi"])
    with open(HARD, "rb") as fh:
        data = pickle.loads(fh.read().replace(b"\r\n", b"\n"), encoding="latin1")
    out.update(data.values())
    return sorted(out)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        files = [f for f in SRC.glob("*.c") if f.stem != "inchi_gui"]
        files += [API / f"{n}.c" for n in
                  ("inchi_dll", "inchi_dll_a", "inchi_dll_a2", "inchi_dll_b", "inchi_dll_main", "ichilnct")]
        (tmp / "gs.c").write_text(HARNESS)
        files.append(tmp / "gs.c")
        procs = []
        objs = []
        for f in files:
            o = tmp / (f.stem + ".o")
            objs.append(str(o))
            procs.append(subprocess.Popen(
                ["cc", "-O1", "-w", "-DTARGET_API_LIB", f"-I{SRC}", f"-I{API}", "-c", str(f), "-o", str(o)]))
            if len(procs) >= 8:
                for p in procs:
                    assert p.wait() == 0
                procs = []
        for p in procs:
            assert p.wait() == 0
        subprocess.run(["cc", *objs, "-lm", "-o", str(tmp / "gs")], check=True)
        ins = inchis()
        res = subprocess.run([str(tmp / "gs")], input="\n".join(ins) + "\n", capture_output=True, text=True, check=True)
    blocks = res.stdout.split("E\n")[:-1]
    assert len(blocks) == len(ins)
    cases = [{"in": i, "out": b.strip().split("\n")} for i, b in zip(ins, blocks)]
    (ROOT / "tests/data/inchi_getstruct.json").write_text(json.dumps(cases))
    print(len(cases), "cases")


main()
