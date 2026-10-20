"""Generates tests/data/inchi_strutil.json from the InChI C code: the drawing
corrections of strutil.c (`fix_non_uniform_drawn_oxoanions`,
`fix_non_uniform_drawn_amidiniums`, `the_only_doublet_neigh`, `nFindOneOM`,
`cmp_iso_atw_diff_component_no`) on random atom arrays. The function texts are
cut out of strutil.c. Needs a C compiler and the fetched InChI sources."""

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
@FUNCS@
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n;
        scanf("%d", &n);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, chg, rad, iso, nh, comp;
            scanf("%d %d %d %d %d %d %d %d", &v, &cbv, &el, &chg, &rad, &iso, &nh, &comp);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].charge = chg;
            at[i].radical = rad; at[i].iso_atw_diff = iso; at[i].num_H = nh; at[i].component = comp;
            for (int j = 0; j < v; j++) { int x, b, s; scanf("%d %d %d", &x, &b, &s); at[i].neighbor[j] = x; at[i].bond_type[j] = b; at[i].bond_stereo[j] = s; }
        }
        /* queries */
        int nq; scanf("%d", &nq);
        for (int q = 0; q < nq; q++) {
            int atno, num, ord[MAXVAL];
            scanf("%d %d", &atno, &num);
            for (int k = 0; k < num; k++) scanf("%d", &ord[k]);
            int r = nFindOneOM(at, atno, ord, num);
            printf("Q %d", r);
            for (int k = 0; k < num && k < MAXVAL; k++) printf(" %d", ord[k]);
            printf("\n");
        }
        for (int i = 0; i < n; i++) {
            int a = -9, b2 = -9;
            int r = the_only_doublet_neigh(at, i, &a, &b2);
            printf("D %d %d %d\n", r, r >= 0 ? a : -9, r >= 0 ? b2 : -9);
        }
        for (int i = 0; i + 1 < n; i++) printf("C %d\n", cmp_iso_atw_diff_component_no(&at[i], &at[i + 1]));
        int nch1 = 0, nch2 = 0;
        fix_non_uniform_drawn_oxoanions(n, at, &nch1);
        fix_non_uniform_drawn_amidiniums(n, at, &nch2);
        printf("N %d %d\n", nch1, nch2);
        for (int i = 0; i < n; i++) {
            printf("A %d %d %d", at[i].charge, at[i].chem_bonds_valence, at[i].valence);
            for (int k = 0; k < at[i].valence; k++) printf(" %d:%d", at[i].bond_type[k], at[i].bond_stereo[k]);
            printf("\n");
        }
        free(at);
    }
    return 0;
}
"""


def gen(rng):
    n = rng.randrange(2, 10)
    atoms = [{"el": rng.choice([6, 7, 8, 8, 8, 16, 17, 35, 53, 85, 34, 52, 15, 1, 9]), "chg": rng.choice([0, 0, 0, -1, -1, 1, 1, 2]),
              "rad": rng.choice([0, 0, 0, 0, 1, 2, 3]), "iso": rng.choice([0, 0, 0, 1, 2, -1, 3]), "nh": rng.choice([0, 0, 0, 1, 2]),
              "comp": rng.randrange(1, 4), "nb": []} for _ in range(n)]
    # centers with terminal neighbors
    for i in range(n):
        if rng.random() < 0.5:
            if rng.random() < 0.6:
                atoms[i]["el"] = rng.choice([17, 35, 53, 16, 34, 52, 85])
                atoms[i]["chg"] = -1
                atoms[i]["rad"] = rng.choice([0, 0, 0, 1])
            for _ in range(rng.randrange(1, 5)):
                if len(atoms[i]["nb"]) < 6:
                    atoms.append({"el": rng.choice([8, 8, 8, 16, 34, 52, 7, 9]), "chg": rng.choice([0, 0, 0, -1]),
                                  "rad": rng.choice([0, 0, 0, 1, 2]), "iso": rng.choice([0, 0, 1, 2, -1]), "nh": 0, "comp": 1, "nb": []})
                    j = len(atoms) - 1
                    t = rng.choice([1, 2, 2, 2])
                    atoms[i]["nb"].append((j, t, 0))
                    atoms[j]["nb"].append((i, t, 0))
    if rng.random() < 0.25:
        # an amidinium cation: Y+ with N, N and a third single-bonded neighbor
        y = len(atoms)
        atoms.append({"el": rng.choice([6, 6, 16, 15]), "chg": 1, "rad": 0, "iso": 0, "nh": 0, "comp": 1, "nb": []})
        for el in (7, 7, rng.choice([6, 8, 7])):
            atoms.append({"el": el, "chg": rng.choice([0, 0, 0, 1]), "rad": 0, "iso": 0, "nh": rng.choice([0, 1, 2]), "comp": 1, "nb": []})
            j = len(atoms) - 1
            atoms[y]["nb"].append((j, 1, rng.choice([0, 0, 1])))
            atoms[j]["nb"].append((y, 1, 0))
    m = len(atoms)
    for _ in range(rng.randrange(0, 4)):
        a, b = rng.randrange(m), rng.randrange(m)
        if a != b and all(x != b for x, _, _ in atoms[a]["nb"]) and len(atoms[a]["nb"]) < 8 and len(atoms[b]["nb"]) < 8:
            t = rng.choice([1, 1, 2, 3])
            st = rng.choice([0, 0, 1, -1, 4])
            atoms[a]["nb"].append((b, t, st))
            atoms[b]["nb"].append((a, t, -st))
    return atoms


def main():
    rng = random.Random(61)
    cases = []
    for _ in range(1500):
        atoms = gen(rng)
        queries = []
        for _ in range(rng.randrange(0, 4)):
            i = rng.randrange(len(atoms))
            v = len(atoms[i]["nb"])
            if v:
                k = rng.randrange(1, v + 1)
                queries.append([i, rng.sample(range(v), k)])
        cases.append({"atoms": atoms, "q": queries})
    stdin = f"{len(cases)}\n"
    for c in cases:
        atoms = c["atoms"]
        stdin += f"{len(atoms)}\n"
        for a in atoms:
            cbv = sum(t for _, t, _ in a["nb"])
            a["cbv"] = cbv
            stdin += f"{len(a['nb'])} {cbv} {a['el']} {a['chg']} {a['rad']} {a['iso']} {a['nh']} {a['comp']} " + " ".join(f"{x} {t} {s}" for x, t, s in a["nb"]) + "\n"
        stdin += f"{len(c['q'])}\n" + "".join(f"{i} {len(o)} " + " ".join(map(str, o)) + "\n" for i, o in c["q"])
    text = (SRC / "strutil.c").read_text(errors="replace")

    def cut(start):
        a = text.index(start)
        return text[a:text.index("\n}\n", a) + 3]

    funcs = [cut("int cmp_iso_atw_diff_component_no( const void *a1, const void *a2 )\n{"),
             cut("int the_only_doublet_neigh( inp_ATOM *at,\n"),
             cut("int fix_non_uniform_drawn_oxoanions( int num_atoms,\n"),
             cut("int fix_non_uniform_drawn_amidiniums( int num_atoms,\n"),
             cut("int nFindOneOM( inp_ATOM *at, int at_no, int ord_OM[], int num_OM )\n{")]
    protos = """int cmp_iso_atw_diff_component_no( const void *a1, const void *a2 );
int the_only_doublet_neigh( inp_ATOM *at, int i1, int *ineigh1, int *ineigh2 );
int fix_non_uniform_drawn_oxoanions( int num_atoms, inp_ATOM *at, int *num_changes );
int fix_non_uniform_drawn_amidiniums( int num_atoms, inp_ATOM *at, int *num_changes );
int nFindOneOM( inp_ATOM *at, int at_no, int ord_OM[], int num_OM );
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNCS@", protos + "\n".join(funcs)))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    k = 0
    data = []
    for c in cases:
        n = len(c["atoms"])
        nq = len(c["q"])
        lines = out[k:k + nq + n + max(n - 1, 0) + 1 + n]
        k += len(lines)
        data.append({"in": c, "out": lines})
    path = ROOT / "tests/data/inchi_strutil.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
