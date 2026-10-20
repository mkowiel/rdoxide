"""Generates tests/data/inchi_ringsys.json from the InChI C code:
`MarkRingSystemsInp` of ichinorm.c (blocks, cut vertices, ring systems) and
the small `UnMark*`/`subtract_DT_from_num_H` helpers on random molecular
graphs (trees, rings, fused/bridged systems, several components, invalid
start atoms). The function texts are cut out of the C sources. Needs a C
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
#include "mode.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "extr_ct.h"
#include "ichinorm.h"
@FUNCS@
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, start;
        scanf("%d %d", &n, &start);
        inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, nh, ih0, ih1, ih2;
            scanf("%d %d %d %d %d", &v, &nh, &ih0, &ih1, &ih2);
            at[i].valence = v; at[i].num_H = nh;
            at[i].num_iso_H[0] = ih0; at[i].num_iso_H[1] = ih1; at[i].num_iso_H[2] = ih2;
            for (int j = 0; j < v; j++) { int q; scanf("%d", &q); at[i].neighbor[j] = q; }
            at[i].at_type = 5; at[i].cFlags = 3; at[i].component = 4; at[i].orig_compt_at_numb = 6;
        }
        int r = MarkRingSystemsInp(at, n, start);
        printf("R %d\n", r);
        for (int i = 0; i < n; i++) printf("A %d %d %d %d\n", at[i].nBlockSystem, at[i].bCutVertex, at[i].nRingSystem, at[i].nNumAtInRingSystem);
        UnMarkRingSystemsInp(at, n);
        UnMarkOtherIndicators(at, n);
        UnMarkOneComponent(at, n);
        subtract_DT_from_num_H(n, at);
        for (int i = 0; i < n; i++) printf("B %d %d %d %d %d %d %d\n", at[i].nBlockSystem, at[i].bCutVertex, at[i].nRingSystem, at[i].at_type, at[i].cFlags, at[i].component, at[i].num_H);
        free(at);
    }
    return 0;
}
"""


def gen(rng):
    n = rng.randrange(1, 25)
    nb = [[] for _ in range(n)]

    def link(a, b):
        if a != b and b not in nb[a] and len(nb[a]) < 6 and len(nb[b]) < 6:
            nb[a].append(b)
            nb[b].append(a)

    kind = rng.choice(["tree", "rings", "random", "fused"])
    for i in range(1, n):
        if rng.random() < 0.9:
            link(i, rng.randrange(i))
    extra = {"tree": 0, "rings": rng.randrange(1, 4), "random": rng.randrange(0, 8), "fused": rng.randrange(2, 7)}[kind]
    for _ in range(extra):
        a = rng.randrange(n)
        b = (a + rng.randrange(2, 7)) % n if kind in ("rings", "fused") else rng.randrange(n)
        link(a, b)
    for l in nb:
        rng.shuffle(l)
    return nb


def main():
    rng = random.Random(8)
    cases = []
    for _ in range(15000):
        nb = gen(rng)
        n = len(nb)
        start = rng.randrange(n) if rng.random() < 0.97 else n + rng.randrange(0, 2)
        atoms = [{"nb": l, "nh": rng.choice([0, 1, 2, 3]), "ih": [rng.choice([0, 0, 1]) for _ in range(3)]} for l in nb]
        cases.append({"atoms": atoms, "start": start})
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{len(c['atoms'])} {c['start']}\n"
        for a in c["atoms"]:
            stdin += f"{len(a['nb'])} {a['nh']} {a['ih'][0]} {a['ih'][1]} {a['ih'][2]} " + " ".join(map(str, a["nb"])) + "\n"
    text = (SRC / "ichinorm.c").read_text(errors="replace")
    run3 = (SRC / "runichi3.c").read_text(errors="replace")

    def cut(t, start):
        a = t.rindex(start)
        return t[a:t.index("\n}\n", a) + 3]

    funcs = [cut(text, "int MarkRingSystemsInp( inp_ATOM *at, int num_atoms, int start )\n{"),
             cut(run3, "int UnMarkRingSystemsInp( inp_ATOM *at, int num_atoms )\n{"),
             cut(text, "int UnMarkOtherIndicators( inp_ATOM *at, int num_atoms )\n{"),
             cut(text, "int UnMarkOneComponent( inp_ATOM *at, int num_atoms )\n{"),
             cut(text, "int subtract_DT_from_num_H( int num_atoms, inp_ATOM *at )\n")]
    protos = """int MarkRingSystemsInp( inp_ATOM *at, int num_atoms, int start );
int UnMarkRingSystemsInp( inp_ATOM *at, int num_atoms );
int UnMarkOtherIndicators( inp_ATOM *at, int num_atoms );
int UnMarkOneComponent( inp_ATOM *at, int num_atoms );
int subtract_DT_from_num_H( int num_atoms, inp_ATOM *at );
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
        cnt = 1 + (n if out[k] != "R -30000" else 0) * 0
        # an invalid start gives "R -30000" followed by the untouched atoms
        cnt = 1 + n + n
        data.append({"in": c, "out": out[k:k + cnt]})
        k += cnt
    path = ROOT / "tests/data/inchi_ringsys.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
