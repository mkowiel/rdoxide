"""Generates tests/data/inchi_components.json from the InChI C code:
`MarkDisconnectedComponents`, `ExtractConnectedComponent`,
`SetConnectedComponentNumber`, `bNumHeterAtomHasIsotopicH`,
`bHeteroAtomMayHaveXchgIsoH` and `get_iat_number` of strutil.c on random
atom arrays. The function texts are cut out of strutil.c. Needs a C compiler
and the fetched InChI sources."""

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
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, bproc, prevnc, hasold;
        scanf("%d %d %d %d", &n, &bproc, &prevnc, &hasold);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, chg, rad, iso, nh, ih0, ih1, ih2, comp;
            scanf("%d %d %d %d %d %d %d %d %d %d %d", &v, &cbv, &el, &chg, &rad, &iso, &nh, &ih0, &ih1, &ih2, &comp);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].charge = chg;
            at[i].radical = rad; at[i].iso_atw_diff = iso; at[i].num_H = nh; at[i].component = comp;
            at[i].num_iso_H[0] = ih0; at[i].num_iso_H[1] = ih1; at[i].num_iso_H[2] = ih2;
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        ORIG_ATOM_DATA od; memset(&od, 0, sizeof(od));
        od.at = at; od.num_inp_atoms = n; od.num_components = prevnc;
        if (hasold) {
            od.nOldCompNumber = calloc(prevnc + 2, sizeof(AT_NUMB));
            for (int q = 0; q < prevnc + 1; q++) { int x; scanf("%d", &x); od.nOldCompNumber[q] = x; }
        }
        printf("B %d\n", bNumHeterAtomHasIsotopicH(at, n));
        printf("H");
        for (int i = 0; i < n; i++) printf(" %d:%d", get_iat_number(at[i].el_number), bHeteroAtomMayHaveXchgIsoH(at, i));
        printf("\n");
        int hadold = od.nOldCompNumber != NULL;
        int r = MarkDisconnectedComponents(&od, bproc);
        printf("R %d %d\n", r, od.num_components);
        printf("C");
        for (int i = 0; i < n; i++) printf(" %d", at[i].component);
        printf("\nL");
        for (int i = 0; i < r; i++) printf(" %d", od.nCurAtLen[i]);
        printf("\nO");
        int osz = (r > prevnc ? r : prevnc);
        if (r > 0) for (int i = 0; i < osz + 1; i++) printf(" %d", od.nOldCompNumber[i]);
        printf("\n");
        for (int k = 1; k <= r; k++) {
            inp_ATOM *ca = calloc(n, sizeof(inp_ATOM));
            int m = ExtractConnectedComponent(at, n, k, ca);
            printf("X %d", m);
            for (int i = 0; i < m; i++) { printf(" |%d", ca[i].orig_compt_at_numb); for (int j = 0; j < ca[i].valence; j++) printf(" %d", ca[i].neighbor[j]); }
            printf("\n");
            free(ca);
        }
        SetConnectedComponentNumber(at, n, 7);
        printf("S %d %d\n", at[0].component, at[n - 1].component);
        free(at);
    }
    return 0;
}
"""

ELS = [1, 6, 6, 7, 8, 8, 15, 16, 34, 52, 9, 17, 35, 53, 11, 26]


def gen(rng):
    n = rng.randrange(1, 13)
    nb = [[] for _ in range(n)]
    for _ in range(rng.randrange(0, n + 2)):
        a, b = rng.randrange(n), rng.randrange(n)
        if a != b and all(x != b for x, _ in nb[a]):
            t = rng.choice([1, 1, 2, 3])
            nb[a].append((b, t))
            nb[b].append((a, t))
    prevnc = rng.randrange(0, 5)
    atoms = []
    for i in range(n):
        el = rng.choice(ELS)
        chg = rng.choice([0, 0, 0, 1, -1, 2])
        rad = rng.choice([0] * 8 + [1, 2])
        iso = rng.choice([0, 0, 0, 1, 2])
        ih = [rng.choice([0, 0, 0, 1]) for _ in range(3)]
        atoms.append({"el": el, "chg": chg, "rad": rad, "iso": iso, "nh": rng.choice([0, 0, 1, 2]), "ih": ih,
                      "comp": rng.randrange(0, prevnc + 1) if rng.random() < 0.9 else 0,
                      "nb": nb[i], "cbv": sum(t for _, t in nb[i])})
    old = [rng.choice([0, 0, 1, 2, 3, 4]) for _ in range(prevnc + 1)]
    return atoms, prevnc, old


def main():
    rng = random.Random(53)
    cases = []
    for _ in range(12000):
        atoms, prevnc, old = gen(rng)
        cases.append({"atoms": atoms, "bproc": rng.choice([0, 1, 1]), "prevnc": prevnc, "old": old,
                      "hasold": rng.choice([0, 1, 1, 1])})
    stdin = f"{len(cases)}\n"
    for c in cases:
        atoms = c["atoms"]
        stdin += f"{len(atoms)} {c['bproc']} {c['prevnc']} {c['hasold']}\n"
        for a in atoms:
            stdin += (f"{len(a['nb'])} {a['cbv']} {a['el']} {a['chg']} {a['rad']} {a['iso']} {a['nh']} "
                      f"{a['ih'][0]} {a['ih'][1]} {a['ih'][2]} {a['comp']} " + " ".join(f"{x} {t}" for x, t in a["nb"]) + "\n")
        if c["hasold"]:
            stdin += " ".join(map(str, c["old"])) + "\n"
    text = (SRC / "strutil.c").read_text(errors="replace")

    def cut(start):
        a = text.rindex(start)
        return text[a:text.index("\n}\n", a) + 3]

    funcs = [cut("int get_iat_number( int el_number )\n{"),
             cut("int bHeteroAtomMayHaveXchgIsoH( inp_ATOM *atom, int iat )\n{"),
             cut("int bNumHeterAtomHasIsotopicH( inp_ATOM *atom, int num_atoms )\n{"),
             cut("int cmp_components( const void *a1, const void *a2 )\n{"),
             cut("int MarkDisconnectedComponents( ORIG_ATOM_DATA *orig_at_data,\n"),
             cut("int ExtractConnectedComponent( inp_ATOM *at,\n"),
             cut("int SetConnectedComponentNumber( inp_ATOM *at, int num_at, int component_number )\n{")]
    protos = """int get_iat_number( int el_number );
int bHeteroAtomMayHaveXchgIsoH( inp_ATOM *atom, int iat );
int bNumHeterAtomHasIsotopicH( inp_ATOM *atom, int num_atoms );
int cmp_components( const void *a1, const void *a2 );
int MarkDisconnectedComponents( ORIG_ATOM_DATA *orig_at_data, int bProcessOldCompNumbers );
int ExtractConnectedComponent( inp_ATOM *at, int num_at, int component_number, inp_ATOM *component_at );
int SetConnectedComponentNumber( inp_ATOM *at, int num_at, int component_number );
enum { IAT_H = 0, IAT_C, IAT_N, IAT_P, IAT_O, IAT_S, IAT_Se, IAT_Te, IAT_F, IAT_Cl, IAT_Br, IAT_I, IAT_MAX = 12 };
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNCS@", protos + "\n".join(funcs)))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    k = 0
    data = []
    for c in cases:
        # lines: B, H, R, C, L, O, X*r, S
        r = int(out[k + 2].split()[1])
        cnt = 6 + max(r, 0) + 1
        data.append({"in": c, "out": out[k:k + cnt]})
        k += cnt
    path = ROOT / "tests/data/inchi_components.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
