"""Generates tests/data/inchi_mol2atom.json from the InChI C code: the Molfile
-> `inp_ATOM` conversion of mol2atom.c (`CreateOrigInpDataFromMolfile`,
`ReadMolfileToInpAtoms`, `MakeInpAtomsFromMolfileData`, `calculate_valences`,
`SetInpAtomsXYZ`, `SetExtOrigAtDataByMolfileExtInput`) on random and fuzzed
V2000/V3000 molfiles. Links the upstream InChI base library (see
inchi_markalt_reference.py).

    python tools/inchi_mol2atom_reference.py
"""
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


m1 = load("molfmt1", "tools/inchi_molfmt1_reference.py")
m2 = load("molfmt2", "tools/inchi_molfmt2_reference.py")
m3 = load("molfmt3", "tools/inchi_molfmt3_reference.py")
m4 = load("molfmt4", "tools/inchi_molfmt4_reference.py")
ma = m1.ma
SRC = ma.SRC

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mode.h"
#include "ichi.h"
#include "ichierr.h"
#include "util.h"
#include "ichi_io.h"
#include "mol_fmt.h"
#include "inpdef.h"
#include "ichicomn.h"
#include "ichimain.h"
/* mol_fmt1.c is compiled here without optimisation: it overflows `stmp[4]` in the S-group reader (undefined
   behaviour whose effect depends on the compiler); the other upstream objects come from the library */
#include "mol_fmt1.c"
int CreateCompositeNormAtom(COMP_ATOM_DATA *a, INP_ATOM_DATA2 *b, int c) { return 0; }
extern int CreateOrigInpDataFromMolfile(INCHI_IOSTREAM *inp_file, ORIG_ATOM_DATA *orig_at_data, int bMergeAllInputStructures,
    int bGetOrigCoord, int bDoNotAddH, int treat_polymers, int treat_NPZz, const char *pSdfLabel, char *pSdfValue,
    unsigned long *lSdfId, long *lMolfileNumber, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr, int bNoWarnings);
static void hex(const char *s, int n) { for (int i = 0; i < n; i++) printf("%02x", (unsigned char) s[i]); if (!n) printf("-"); }
static void hexs(const char *s) { hex(s, (int) strlen(s)); }
static void pd(double d) { if (isnan(d)) { printf("nan"); return; } unsigned long long u; memcpy(&u, &d, 8); printf("%016llx", u); }
static int unhex(const char *h, char *out) {
    int n = 0;
    out[0] = 0;
    if (h[0] == '-') return 0;
    for (; h[0] && h[1]; h += 2) { unsigned v; sscanf(h, "%2x", &v); out[n++] = (char) v; }
    out[n] = 0;
    return n;
}
static char text[1 << 17], hx[1 << 18];
static void dump_at(inp_ATOM *at, int n) {
    for (int i = 0; i < n && i < 100; i++) {
        inp_ATOM *a = &at[i];
        printf("A %d ", i); hexs(a->elname);
        printf(" %d %d %d %d %d %d %d %d %d %d %d ", a->el_number, a->orig_at_number, a->valence, a->chem_bonds_valence, a->num_H,
               a->num_iso_H[0], a->num_iso_H[1], a->num_iso_H[2], a->iso_atw_diff, a->charge, a->radical);
        pd(a->x); printf(" "); pd(a->y); printf(" "); pd(a->z);
        for (int j = 0; j < a->valence && j < MAXVAL; j++) printf(" %d:%d:%d", a->neighbor[j], a->bond_type[j], a->bond_stereo[j]);
        printf("\n");
    }
}
static void dump_lists_arr(int **l, int n, int len_at, int add) {
    printf(" %d", n);
    if (l) for (int i = 0; i < n; i++) {
        printf(" [");
        for (int k = 0; k < l[i][len_at] + add && k < 40; k++) printf("%d,", l[i][k]);
        printf("]");
    }
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int use_od, merge, bcoord, noaddh, polymers, npzz, large, label_kind, use_val, use_id, nowarn;
        printf("#CASE %d\n", cs);
        scanf("%d %d %d %d %d %d %d %d %d %d %d %s", &use_od, &merge, &bcoord, &noaddh, &polymers, &npzz, &large, &label_kind, &use_val, &use_id, &nowarn, hx);
        int n = unhex(hx, text);
        INCHI_IOSTREAM ios;
        inchi_ios_init(&ios, INCHI_IOS_TYPE_STRING, NULL);
        ios.s.pStr = (char *) malloc(n + 2);
        memcpy(ios.s.pStr, text, n + 1);
        ios.s.nUsedLength = n; ios.s.nAllocatedLength = n + 2; ios.s.nPtr = 0;
        char err[STR_ERR_LEN]; err[0] = 0;
        char value[300]; memset(value, 0, sizeof(value));
        unsigned long id = 777;
        long molnum = -5;
        const char *labels[] = { NULL, "", "CAS", "x", "MolfileName", "MolfileLine2", "MolfileComment", "MolfileIntRegNo" };
        ORIG_ATOM_DATA od; memset(&od, 0, sizeof(od));
        INCHI_MODE flags = large ? FLAG_SET_INP_LARGE_MOLS : 0;
        int e = 0;
        int na = CreateOrigInpDataFromMolfile(&ios, use_od ? &od : NULL, merge, bcoord, noaddh, polymers, npzz, labels[label_kind],
                                              use_val ? value : NULL, use_id ? &id : NULL, &molnum, &flags, &e, err, nowarn);
        printf("R %d %d %ld %lu %lu ", na, e, molnum, (unsigned long) flags, id); hexs(value); printf(" "); hexs(err); printf("\n");
        if (use_od && !e) {
            printf("O %d %d %d\n", od.num_inp_atoms, od.num_inp_bonds, od.num_dimensions);
            dump_at(od.at, od.num_inp_atoms);
            if (od.szCoord) for (int i = 0; i < od.num_inp_atoms && i < 100; i++) { printf("K %d ", i); hexs(od.szCoord[i]); printf("\n"); }
            if (od.polymer) {
                OAD_Polymer *p = od.polymer;
                printf("P %d %d %d %d %d %d %d\n", p->n, p->n_pzz, p->really_do_frame_shift, p->frame_shift_scheme, p->treat, p->is_in_reconn, p->edit_repeats);
                for (int k = 0; k < p->n; k++) {
                    OAD_PolymerUnit *u = p->units[k];
                    printf("U %d %d %d %d %d %d %d ", u->id, u->type, u->subtype, u->conn, u->label, u->na, u->nb);
                    for (int q = 0; q < 4; q++) { pd(u->xbr1[q]); printf(" "); }
                    for (int q = 0; q < 4; q++) { pd(u->xbr2[q]); printf(" "); }
                    hexs(u->smt);
                    for (int q = 0; q < u->na; q++) printf(" %d", u->alist[q]);
                    printf(" |");
                    for (int q = 0; q < 2 * u->nb; q++) printf(" %d", u->blist[q]);
                    printf("\n");
                }
            }
            if (od.v3000) {
                OAD_V3000 *v = od.v3000;
                printf("V %d %d %d %d %d %d %d %d %d %d %d %d", v->n_non_star_atoms, v->n_star_atoms, v->n_sgroups, v->n_3d_constraints,
                       v->n_collections, v->n_non_haptic_bonds, v->n_haptic_bonds, v->n_steabs, v->n_sterel, v->n_sterac,
                       v->atom_index_orig != NULL, v->atom_index_fin != NULL);
                if (v->atom_index_orig) for (int i = 0; i < 16; i++) printf(" %d", ((unsigned char *) v->atom_index_orig)[i]);
                printf("\nW"); dump_lists_arr(v->lists_haptic_bonds, v->n_haptic_bonds, 2, 3);
                printf("\nS"); dump_lists_arr(v->lists_steabs, v->n_steabs, 1, 2); dump_lists_arr(v->lists_sterel, v->n_sterel, 1, 2); dump_lists_arr(v->lists_sterac, v->n_sterac, 1, 2);
                printf("\n");
            }
        }
        inchi_ios_free_str(&ios);
    }
    return 0;
}
"""

R = random.Random(5150)
m2.SYMS = ["C", "C", "C", "N", "O", "O", "H", "Cl", "CL", "Br", "S", "P", "F", "Na", "Fe", "D", "T", "Zz", "*", "Zy", "He", "Xx", "Si", "B"]
m3.SYMS = m2.SYMS


def sdf_tail(r):
    L = []
    if r.random() < 0.7:
        L += [r.choice(["> <x>", "> <NAME>", "> <CAS_NO>"]), r.choice(["ethanol", "64-17-5", "a b"]), ""]
    if r.random() < 0.8:
        L.append("$$$$")
    return L


def mk_clean(r):
    na = r.randrange(1, 12)
    syms = ["C", "C", "C", "C", "N", "N", "O", "O", "H", "Cl", "S", "P", "F", "Na", "Fe", "Br", "Si", "B", "D", "T", "Zz", "*"]
    nb = r.randrange(0, na + 3) if na > 1 else 0
    props = []
    for _ in range(r.randrange(0, 4)):
        k = r.choice(["CHG", "RAD", "ISO", "A", "SGR"])
        a = r.randrange(1, na + 1)
        if k == "CHG":
            props.append(f"M  CHG  1 {a:3d} {r.randrange(-2, 3):3d}")
        elif k == "RAD":
            props.append(f"M  RAD  1 {a:3d} {r.randrange(0, 4):3d}")
        elif k == "ISO":
            props.append(f"M  ISO  1 {a:3d} {r.choice([13, 2, 18, 35, 1, 3, 15]):3d}")
        elif k == "A":
            props += [f"A  {a:3d}", r.choice(["CH3", "OH", "D", "T", "NH2", "Me+", "H2O", "C2H5", "XYZ"])]
        else:
            props += [f"M  STY  1   1 SRU", f"M  SAL   1  2 {a:3d} {r.randrange(1, na + 1):3d}"]
            if nb:
                props.append(f"M  SBL   1  1 {r.randrange(1, nb + 1):3d}")
            props.append("M  SMT   1 n")
    L = ["Structure #3", "  InChI 1.07  SDfile Output", "", f"{na:3d}{nb:3d}  0  0{r.choice([0, 1]):3d}  0  0  0  0  0{len(props) + 1:3d} V2000"]
    for _ in range(na):
        L.append(f"{r.uniform(-6, 6):10.4f}{r.uniform(-6, 6):10.4f}{r.choice([0, 0, 1]) * r.uniform(-3, 3):10.4f} {r.choice(syms):<3} 0{r.choice([0, 0, 0, 0, 3, 5, 4]):3d}{r.choice([0, 0, 1, 2]):3d}  0  0{r.choice([0] * 5 + [1, 3, 4, 15]):3d}  0  0  0  0  0")
    for _ in range(nb):
        a1, a2 = r.randrange(1, na + 1), r.randrange(1, na + 1)
        L.append(f"{a1:3d}{a2:3d}{r.choice([1] * 5 + [2, 2, 3, 4, 4, 4, 0, 8]):3d}{r.choice([0] * 6 + [1, 4, 6, 3, 2]):3d}  0  0")
    L += props + ["M  END"]
    return "\n".join(L) + "\n"


def mk_text(r, merge):
    def one():
        x = r.random()
        body = mk_clean(r) if x < 0.7 else (m2.mk_mol2(r) if x < 0.85 else m3.mk_mol3(r))
        nl = "\r\n" if "\r\n" in body else "\n"
        return body + nl.join(sdf_tail(r)) + nl
    t = one()
    if merge and r.random() < 0.7:
        t += one()
    return t


def main():
    cases = []
    for _ in range(2400):
        merge = R.choice([0, 0, 0, 1])
        cases.append(f"{R.choice([1] * 19 + [0])} {merge} {R.choice([0, 1])} {R.choice([0, 0, 1])} {R.choice([0, 1, 1])} {R.choice([0, 1, 2])} "
                     f"{R.choice([0, 0, 1])} {R.randrange(0, 8)} {R.choice([0, 1])} {R.choice([0, 1])} {R.choice([0, 1])} {mk_text(R, merge).encode().hex()}")
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-ffp-contract=off", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        kept, out = [], []
        for c in cases:
            try:
                r1 = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=10)
            except subprocess.TimeoutExpired:
                continue
            if r1.returncode == 0:
                again = [subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=10) for _ in range(3)]
                if any(a.returncode != 0 or a.stdout != r1.stdout for a in again):
                    continue
                kept.append(c)
                out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode(errors="replace").split("\n")[1:]
        print(len(cases) - len(kept), "cases dropped")
        cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    data = [x for x in data if all(len(l.split()[2]) <= 62 for l in x["out"] if l.startswith("K "))]
    # merging several structures leaves dangling polymer/V3000 pointers in the C code (freed after being stored)
    data = [x for x in data if not (x["in"].split()[1] == "1" and any(l[:2] in ("P ", "V ", "W ", "S ", "U ") for l in x["out"]))]
    path = ROOT / "tests/data/inchi_mol2atom.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
