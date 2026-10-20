"""Generates tests/data/inchi_metals.json from the InChI C code:
`DisconnectMetals`, `DisconnectOneLigand`, `move_explicit_Hcation`,
`GetMinDistDistribution` and `dist3D` of strutil.c on random metal complexes
with 2D coordinates. The function texts are cut out of strutil.c. Needs a C
compiler and the fetched InChI sources."""

import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src"

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "mode.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "extr_ct.h"
@FUNCS@
void inchi_swap( char *a, char *b, size_t width ) { char tmp; if (a != b) { while (width--) { tmp = *a; *a++ = *b; *b++ = tmp; } } }
static uint64_t bits(double d) { uint64_t u; memcpy(&u, &d, 8); return u; }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, check;
        scanf("%d %d", &n, &check);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, chg, rad, iso, nh, ih0, ih1, ih2, comp;
            double x, y, z;
            scanf("%d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf", &v, &cbv, &el, &chg, &rad, &iso, &nh, &ih0, &ih1, &ih2, &comp, &x, &y, &z);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].charge = chg;
            at[i].radical = rad; at[i].iso_atw_diff = iso; at[i].num_H = nh; at[i].component = comp;
            at[i].num_iso_H[0] = ih0; at[i].num_iso_H[1] = ih1; at[i].num_iso_H[2] = ih2;
            at[i].x = x; at[i].y = y; at[i].z = z; at[i].orig_at_number = i + 1;
            at[i].elname[0] = 'X';
            for (int j = 0; j < v; j++) { int q, b, s; scanf("%d %d %d", &q, &b, &s); at[i].neighbor[j] = q; at[i].bond_type[j] = b; at[i].bond_stereo[j] = s; }
        }
        ORIG_ATOM_DATA od; memset(&od, 0, sizeof(od));
        AT_NUMB oldcomp[16];
        for (int q = 0; q < 16; q++) oldcomp[q] = 1;
        od.at = at; od.num_inp_atoms = n; od.nOldCompNumber = oldcomp;
        INCHI_MODE fl = 0;
        int may = bMayDisconnectMetals(&od, check, &fl);
        int r = DisconnectMetals(&od, check, &fl);
        printf("R %d %d %d %d %d %d %d %d\n", may, od.bDisconnectCoord, r, od.num_inp_atoms, (int)fl, oldcomp[0], oldcomp[1], oldcomp[2]);
        at = od.at;
        for (int i = 0; i < od.num_inp_atoms; i++) {
            printf("A %d %d %d %d %d %d %d %d %d %d %d %d %llu %llu %llu", at[i].el_number, at[i].valence, at[i].chem_bonds_valence, at[i].charge, at[i].num_H,
                   at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2], at[i].iso_atw_diff, at[i].radical, at[i].component, at[i].orig_at_number,
                   (unsigned long long)bits(at[i].x), (unsigned long long)bits(at[i].y), (unsigned long long)bits(at[i].z));
            for (int k = 0; k < at[i].valence; k++) printf(" %d:%d:%d", at[i].neighbor[k], at[i].bond_type[k], at[i].bond_stereo[k]);
            printf("\n");
        }
    }
    return 0;
}
"""

METALS = [3, 11, 19, 12, 20, 30, 13, 26, 29, 37, 50, 82, 47]
RNG = random.Random(3)


def new(atoms, el, comp=1, chg=0, nh=0, rad=0):
    atoms.append({"el": el, "chg": chg, "rad": rad, "iso": 0, "nh": nh, "ih": [0, 0, 0], "comp": comp,
                  "xyz": [round(RNG.uniform(-4, 4), 3), round(RNG.uniform(-4, 4), 3), RNG.choice([0.0, 0.0, 0.5])], "nb": []})
    return len(atoms) - 1


def link(atoms, a, b, t=1):
    atoms[a]["nb"].append((b, t, 0))
    atoms[b]["nb"].append((a, t, 0))


def gen(rng):
    atoms = []
    comp = 1
    metals = []
    for _ in range(rng.randrange(1, 3)):
        m = new(atoms, rng.choice(METALS), comp, rng.choice([0, 0, 1, 2, -1]), rng.choice([0, 0, 0, 1, 2]))
        if rng.random() < 0.2:
            atoms[m]["ih"][rng.randrange(3)] = 1
        metals.append(m)
    # the explicit H added by DisconnectMetals need the coordinates of metals only
    for m in metals:
        for _ in range(rng.randrange(1, 4)):
            el = rng.choice([8, 8, 7, 7, 6, 17, 16, 9, 15, 35])
            lig = new(atoms, el, comp, rng.choice([0, 0, -1, 1]), rng.choice([0, 0, 1]))
            link(atoms, m, lig, rng.choice([1, 1, 1, 2, 4]))
            for _ in range(rng.randrange(0, 3)):
                x = new(atoms, rng.choice([6, 6, 1, 8, 7]), comp)
                link(atoms, lig, x, rng.choice([1, 1, 2, 4]))
            if rng.random() < 0.3 and len(metals) > 1:
                link(atoms, lig, rng.choice([q for q in metals if q != m]), 1) if all(
                    q != lig for q, _, _ in atoms[metals[0]]["nb"] + atoms[metals[-1]]["nb"]) else None
    if len(metals) > 1 and rng.random() < 0.4 and all(q != metals[1] for q, _, _ in atoms[metals[0]]["nb"]):
        link(atoms, metals[0], metals[1], 1)
    # make adjacency lists unique/symmetric
    for i, a in enumerate(atoms):
        seen = set()
        keep = []
        for e in a["nb"]:
            if e[0] not in seen and e[0] != i:
                seen.add(e[0])
                keep.append(e)
        a["nb"] = keep
    for i, a in enumerate(atoms):
        for e in list(a["nb"]):
            if all(q != i for q, _, _ in atoms[e[0]]["nb"]):
                atoms[e[0]]["nb"].append((i, e[1], 0))
    for a in atoms:
        a["cbv"] = sum(t if t < 4 else 1 for _, t, _ in a["nb"])
        if rng.random() < 0.04:
            a["rad"] = 2
    return atoms


def main():
    rng = random.Random(29)
    cases = []
    for _ in range(4000):
        cases.append({"atoms": gen(rng), "check": rng.choice([0, 1])})
    stdin = f"{len(cases)}\n"
    for c in cases:
        atoms = c["atoms"]
        stdin += f"{len(atoms)} {c['check']}\n"
        for a in atoms:
            stdin += (f"{len(a['nb'])} {a['cbv']} {a['el']} {a['chg']} {a['rad']} {a['iso']} {a['nh']} "
                      f"{a['ih'][0]} {a['ih'][1]} {a['ih'][2]} {a['comp']} {a['xyz'][0]} {a['xyz'][1]} {a['xyz'][2]} "
                      + " ".join(f"{x} {t} {s}" for x, t, s in a["nb"]) + "\n")
    text = (SRC / "strutil.c").read_text(errors="replace")
    ister = (SRC / "ichister.c").read_text(errors="replace")

    def cut(t, start):
        a = t.rindex(start)
        return t[a:t.index("\n}\n", a) + 3]

    funcs = [cut(ister, "int get_opposite_sb_atom( inp_ATOM *at,\n"),
             cut(text, "int RemoveInpAtBond( inp_ATOM *atom, int iat, int k )\n{"),
             cut(text, "int DisconnectInpAtBond( inp_ATOM *at,\n"),
             cut(text, "int bIsAmmoniumSalt( inp_ATOM *at,\n"),
             cut(text, "int bIsMetalSalt( inp_ATOM *at, int i )\n{"),
             cut(text, "int bIsMetalToDisconnect( inp_ATOM *at, int i, int bCheckMetalValence )\n{"),
             cut(text, "int bMayDisconnectMetals( ORIG_ATOM_DATA *orig_inp_data,\n"),
             cut(text, "int DisconnectMetals( ORIG_ATOM_DATA *orig_inp_data,\n"),
             cut(text, "int DisconnectOneLigand( inp_ATOM *at,\n"),
             cut(text, "double dist3D( inp_ATOM *at1, inp_ATOM *at2 )\n{"),
             cut(text, "double GetMinDistDistribution( inp_ATOM *at,\n"),
             cut(text, "int move_explicit_Hcation( inp_ATOM *at,\n")]
    protos = """int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt, int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );
int RemoveInpAtBond( inp_ATOM *atom, int iat, int k );
int DisconnectInpAtBond( inp_ATOM *at, AT_NUMB *nOldCompNumber, int iat, int neigh_ord );
int bIsAmmoniumSalt( inp_ATOM *at, int i, int *piO, int *pk, S_CHAR *num_explicit_H );
int bIsMetalSalt( inp_ATOM *at, int i );
int bIsMetalToDisconnect( inp_ATOM *at, int i, int bCheckMetalValence );
int bMayDisconnectMetals( ORIG_ATOM_DATA *orig_inp_data, int bCheckMetalValence, INCHI_MODE *bTautFlagsDone );
int DisconnectMetals( ORIG_ATOM_DATA *orig_inp_data, int bCheckMetalValence, INCHI_MODE *bTautFlagsDone );
int DisconnectOneLigand( inp_ATOM *at, AT_NUMB *nOldCompNumber, S_CHAR *bMetal, char *elnumber_Heteroat, int num_halogens, int num_atoms, int iMetal, int jLigand, INCHI_MODE *bTautFlagsDone );
double dist3D( inp_ATOM *at1, inp_ATOM *at2 );
double GetMinDistDistribution( inp_ATOM *at, int num_at, int iat, int iat_H, int bInAllComponents, double min_dist[], int num_segm );
int move_explicit_Hcation( inp_ATOM *at, int num_at, int iat, int iat_H, int bInAllComponents );
#define MIN_BOND_LENGTH   (1.0e-6)
#define MIN_COS           (1.0e-6)
#define MIN_BOND_LENGTH2  (MIN_BOND_LENGTH*MIN_BOND_LENGTH)
#define MAX_BOND_LENGTH   (1.0e30)
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNCS@", protos + "\n".join(funcs)))
        subprocess.run(["cc", "-w", "-ffp-contract=off", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c"), "-lm"], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    k = 0
    data = []
    changed = 0
    for c in cases:
        head = out[k].split()
        n = len(c["atoms"])
        r = int(head[3])
        n_after = int(head[4])
        lines = out[k:k + 1 + n_after]
        k += 1 + n_after
        changed += r > 0
        data.append({"in": c, "out": lines})
    path = ROOT / "tests/data/inchi_metals.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases,", changed, "disconnected ->", path)


main()
