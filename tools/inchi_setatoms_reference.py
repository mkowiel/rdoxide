"""Generates tests/data/inchi_setatoms.json from the InChI C code: the
conversion of `inchi_Atom` arrays into `inp_ATOM` (SetAtomProperties,
SetBondProperties, SetAtomAndBondProperties, SetNumImplicitH and the loop of
ExtractOneStructure; WriteCoord) on random, often invalid, molecules. The
function texts are cut out of inchi_dll.c and ichimak2.c. Needs a C compiler
and the fetched InChI sources."""

import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src"
DLL = ROOT / "upstream/external/InChI/INCHI_API/libinchi/src/inchi_dll.c"

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mode.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "inchi_api.h"
#include "mol_fmt.h"
@FUNCS@
static void hex(const char *s) { printf("x"); for (; s && *s; s++) printf("%02x", (unsigned char) *s); }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, dna;
        scanf("%d %d", &n, &dna);
        inchi_Atom *ati = calloc(n ? n : 1, sizeof(inchi_Atom));
        for (int i = 0; i < n; i++) {
            char h[64]; int nb, iso, rad, chg, nh[4];
            double x, y, z;
            scanf("%s %lf %lf %lf %d %d %d %d %d %d %d %d", h, &x, &y, &z, &nb, &iso, &rad, &chg, &nh[0], &nh[1], &nh[2], &nh[3]);
            int m = 0;
            for (const char *q = h + 1; q[0] && q[1]; q += 2) { unsigned v; sscanf(q, "%2x", &v); ati[i].elname[m++] = (char) v; }
            ati[i].x = x; ati[i].y = y; ati[i].z = z;
            ati[i].num_bonds = nb; ati[i].isotopic_mass = iso; ati[i].radical = rad; ati[i].charge = chg;
            for (int k = 0; k < 4; k++) ati[i].num_iso_H[k] = nh[k];
            for (int j = 0; j < nb; j++) { int a, b, s; scanf("%d %d %d", &a, &b, &s); ati[i].neighbor[j] = a; ati[i].bond_type[j] = b; ati[i].bond_stereo[j] = s; }
        }
        inp_ATOM *at = calloc(n ? n : 1, sizeof(inp_ATOM));
        MOL_COORD *szCoord = calloc(n ? n : 1, sizeof(MOL_COORD));
        char pStrErr[STR_ERR_LEN]; memset(pStrErr, 0, sizeof pStrErr);
        int err = 0, nDim = 0, nNumBonds = 0, a1, j, valence, stop = 0;
        for (a1 = 0; a1 < n && !stop; a1++) {
            SetAtomProperties(at, szCoord, ati, a1, &nDim, pStrErr, &err);
            if (err) { stop = 1; break; }
            valence = ati[a1].num_bonds;
            for (j = 0; j < valence; j++) SetBondProperties(at, ati, a1, j, n, &nNumBonds, pStrErr, &err);
            if (err) { stop = 1; break; }
        }
        if (!stop) {
            for (a1 = 0; a1 < n && !stop; a1++) {
                SetAtomAndBondProperties(at, ati, a1, dna, pStrErr, &err);
                if (err) stop = 1;
            }
            if (!stop) SetNumImplicitH(at, n);
        }
        printf("%d %d %d ", err, nDim, nNumBonds); hex(pStrErr); printf("\n");
        for (int i = 0; i < n; i++) {
            inp_ATOM *a = &at[i];
            printf("A "); hex(a->elname);
            printf(" %d %d %d %d %d %d %d %d %d %d %d %d %d %d", a->el_number, a->charge, a->radical, a->valence, a->chem_bonds_valence, a->num_H,
                   a->num_iso_H[0], a->num_iso_H[1], a->num_iso_H[2], a->iso_atw_diff, a->at_type, a->orig_at_number, (int) a->x * 0, 0);
            printf(" %.6f %.6f %.6f", a->x, a->y, a->z);
            for (int k = 0; k < a->valence && k < MAXVAL; k++) printf(" %d:%d:%d", a->neighbor[k], a->bond_type[k], a->bond_stereo[k]);
            printf(" |"); { char t[31]; memcpy(t, szCoord[i], 30); t[30] = 0; hex(t); }
            printf("\n");
        }
        free(ati); free(at); free(szCoord);
    }
    return 0;
}
"""

ELEMS = ["C", "C", "N", "O", "O", "H", "D", "T", "S", "P", "Cl", "Fe", "Na", "Xx", "CH3", "NH2", "OH-", "H2", "CH4+", "N+", "O-2", "HD", "Zn", "C2", "", "Si"]


def case(rng, clean):
    n = rng.randrange(1, 8)
    elems = ["C", "C", "N", "O", "O", "H", "D", "T", "S", "P", "Cl", "Fe", "Na", "CH3", "NH2", "OH-", "H2", "CH4+", "N+", "Zn", "Si", "HD"] if clean else ELEMS
    adj = [[] for _ in range(n)]
    atoms = []
    if clean:
        # a random graph with symmetric adjacency lists
        edges = set()
        for i in range(1, n):
            edges.add((rng.randrange(i), i))
        for _ in range(rng.randrange(0, 3)):
            a_, b_ = rng.randrange(n), rng.randrange(n)
            if a_ != b_:
                edges.add((min(a_, b_), max(a_, b_)))
        for a_, b_ in sorted(edges):
            bt = rng.choice([1, 1, 1, 2, 3, 4])
            st = rng.choice([0, 0, 0, 1, 4, 6, -1, 3])
            if rng.random() < 0.8:
                adj[a_].append((b_, bt, st))
                adj[b_].append((a_, bt, -st if st in (1, 4, 6, -1) else st))
            else:
                adj[rng.choice([a_, b_])].append((b_ if True else a_, bt, st))
    for i in range(n):
        if clean:
            nb = adj[i]
            rad = rng.choice([0, 0, 0, 0, 1, 2, 3])
            chg = rng.choice([0, 0, 0, 1, -1])
        else:
            nb = []
            for _ in range(rng.choice([0, 1, 1, 2, 2, 3, 4])):
                tgt = rng.choice(list(range(n)) * 4 + [n + 3, i])
                nb.append((tgt, rng.choice([1, 1, 1, 2, 3, 4, 4, 0, 5]),
                           rng.choice([0, 0, 0, 1, 4, 6, -1, -4, -6, 3, -3, 2, 7])))
            if rng.random() < 0.05:
                nb = [(rng.randrange(n), 1, 0) for _ in range(23)]
            rad = rng.choice([0, 0, 0, 1, 2, 3, 4, 5, 6, 7, -1])
            chg = rng.choice([0, 0, 0, 1, -1, 2, -3])
        atoms.append({
            "el": rng.choice(elems), "xyz": [rng.choice([0, 0, 1.5, -2.25, 12345678.9, -1e8, 0.0000001, 100000.5]) for _ in range(3)],
            "iso": rng.choice([0, 0, 0, 13, 2, 3, 10002, 10003, 9950, 10120, 1, 10000, 20000, 65535]),
            "rad": rad, "chg": chg,
            "nh": [rng.choice([-1, -1, 0, 0, 1, 2, 3]), rng.choice([0, 0, 1]), rng.choice([0, 0, 1]), rng.choice([0, 0, 1])],
            "nb": nb[:20]})
    return {"n": n, "dna": rng.choice([0, 0, 1]), "atoms": atoms}


def main():
    rng = random.Random(31)
    cases = [case(rng, i % 3 != 0) for i in range(1200)]
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{c['n']} {c['dna']}\n"
        for a in c["atoms"]:
            stdin += "x" + (a["el"].encode().hex() or "00")[:60] + " " + " ".join(repr(float(v)) for v in a["xyz"]) + f" {len(a['nb'])} {a['iso']} {a['rad']} {a['chg']} " + " ".join(map(str, a["nh"])) + "\n"
            stdin += " ".join(f"{t} {b} {s}" for t, b, s in a["nb"]) + "\n"
    src = DLL.read_text(errors="replace")
    funcs = []
    for start, end in (("void SetNumImplicitH( inp_ATOM* at, int num_atoms )\n{", "/****"), ("#define MIN_BOND_LENGTH", "int InpAtom0DToInchiAtom( inp_ATOM *at,\n                          int num_inp_atoms,")):
        a = src.index(start, src.index("int parse_options_string( char *cmd, const char *argv[], int maxargs )\n{") if start.startswith("#define") else 0)
        b = src.index("\n}\n", a) + 3 if end == "/****" else src.index(end, a + 5)
        funcs.append(src[a:b])
    mak = (SRC / "ichimak2.c").read_text(errors="replace")
    a = mak.index("void WriteCoord( char *str, double x )")
    wc = mak[a:mak.index("#endif", a)]
    text = "void WriteCoord( char *str, double x );\n" + "\n".join(funcs) + "\n" + wc
    # the stub of the C-level prototypes used before definition
    protos = """int SetAtomProperties( inp_ATOM *at, MOL_COORD *szCoord, inchi_Atom *ati, int a1, int *nDim, char *pStrErr, int *err );
int SetBondProperties( inp_ATOM *at, inchi_Atom *ati, int a1, int j, int nNumAtoms, int *nNumBonds, char *pStrErr, int *err );
int SetAtomAndBondProperties( inp_ATOM *at, inchi_Atom *ati, int a1, int bDoNotAddH, char *pStrErr, int *err );
void SetNumImplicitH( inp_ATOM* at, int num_atoms );
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNCS@", protos + text))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    k = 0
    data = []
    for c in cases:
        head = out[k]; k += 1
        atoms = out[k:k + c["n"]]; k += c["n"]
        data.append({"in": c, "head": head, "atoms": atoms})
    path = ROOT / "tests/data/inchi_setatoms.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
