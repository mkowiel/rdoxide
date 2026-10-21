"""Generates tests/data/inchi_atomtype.json from the InChI C code:
`GetAtomChargeType` (the charge/H type and bit mask of an atom and the
running type totals), `IsZOX` and `is_Z_atom` of ichi_bns.c on random small
structures. The function texts are cut out of ichi_bns.c. Needs a C compiler
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
#include "ichitaut.h"
#include "ichi_bns.h"
@MACROS@
@FUNCS@
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n;
        scanf("%d", &n);
        inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, chg, nh, rad;
            scanf("%d %d %d %d %d %d", &v, &cbv, &el, &chg, &nh, &rad);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].charge = chg; at[i].num_H = nh; at[i].radical = rad;
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        printf("N %d\n", c);
        int sub[4] = {0, 1, 2, -2};
        for (int k = 0; k < 4; k++) {
            int totals[ATTOT_ARRAY_LEN];
            memset(totals, 0, sizeof(totals));
            for (int i = 0; i < n; i++) {
                int mask = -12345;
                int t = GetAtomChargeType(at, i, (i & 1) ? NULL : totals, &mask, sub[k]);
                printf("%d,%d ", t, mask);
            }
            printf("|");
            for (int i = 0; i < ATTOT_ARRAY_LEN; i++) printf(" %d", totals[i]);
            printf("\n");
        }
        for (int i = 0; i < n; i++) {
            printf("z%d", is_Z_atom(at[i].el_number));
            if (at[i].valence) printf(" x%d", IsZOX(at, i, 0));
            printf("\n");
        }
        free(at);
    }
    return 0;
}
"""

ELS = [6, 6, 7, 7, 8, 8, 8, 15, 16, 34, 52, 17, 35, 53, 1, 9, 11, 26, 5]


def gen(rng):
    n = rng.randrange(1, 7)
    nb = [[] for _ in range(n)]
    for i in range(1, n):
        if rng.random() < 0.95:
            j = rng.randrange(i)
            t = rng.choice([1, 1, 1, 2, 3])
            nb[i].append((j, t))
            nb[j].append((i, t))
    atoms = []
    for i in range(n):
        v = len(nb[i])
        bs = sum(t for _, t in nb[i])
        cbv = bs + rng.choice([0, 0, 0, 1, 2]) if rng.random() < 0.1 else bs
        bl = []
        for j, t in nb[i]:
            if rng.random() < 0.08:
                t |= 0x10
            bl.append((j, t))
        atoms.append({"v": v, "cbv": cbv, "el": rng.choice(ELS), "chg": rng.choice([0] * 5 + [1, 1, -1, -1, 2, -2]),
                      "nh": rng.choice([0, 0, 1, 1, 2, 3]), "rad": rng.choice([0] * 12 + [1, 2, 3]), "nb": bl})
    return atoms


def main():
    rng = random.Random(19)
    cases = []
    stdin = ""
    for _ in range(12000):
        atoms = gen(rng)
        cases.append(atoms)
    stdin = f"{len(cases)}\n"
    for atoms in cases:
        stdin += f"{len(atoms)}\n"
        for a in atoms:
            stdin += f"{a['v']} {a['cbv']} {a['el']} {a['chg']} {a['nh']} {a['rad']} " + " ".join(f"{x} {t}" for x, t in a["nb"]) + "\n"
    text = (SRC / "ichi_bns.c").read_text(errors="replace")

    def cut(start):
        a = text.rindex(start)
        return text[a:text.index("\n}\n", a) + 3]

    funcs = [cut("int is_Z_atom( U_CHAR el_number )\n{"), cut("int IsZOX( inp_ATOM *atom, int at_x, int ord )\n{"),
             cut("int GetAtomChargeType( inp_ATOM *atom,\n")]
    a = text.index("typedef enum tagAtTypeTotals")
    b = text.index("#define ATBIT_Errors")
    macros = text[a:text.index("\n", b)]
    protos = "int is_Z_atom( U_CHAR el_number );\nint IsZOX( inp_ATOM *atom, int at_x, int ord );\nint GetAtomChargeType( inp_ATOM *atom, int at_no, int nAtTypeTotals[], int *pMask, int bSubtract );\n"
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@MACROS@", macros).replace("@FUNCS@", protos + "\n".join(funcs)))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c"), "-lm"], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_atomtype.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
