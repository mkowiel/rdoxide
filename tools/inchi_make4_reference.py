"""Generates tests/data/inchi_make4.json from the InChI C code: `GetDfsOrder4CT`
and `mystrrev` of ichimake.c on random linear connection tables. Links the
upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_make4_reference.py
"""
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("markalt", ROOT / "tools/inchi_markalt_reference.py")
ma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ma)
SRC = ma.SRC

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
#include "ichicano.h"
#include "ichicant.h"
#include "ichinorm.h"
#include "ichimake.h"
int CreateCompositeNormAtom() { return 0; }
void mystrrev( char *p );
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, len, mode, hasH;
        scanf("%d %d %d %d", &n, &len, &mode, &hasH);
        AT_NUMB *ct = calloc(len + 1, sizeof(AT_NUMB));
        for (int i = 0; i < len; i++) { int x; scanf("%d", &x); ct[i] = x; }
        S_CHAR *h = calloc(n + 1, 1);
        for (int i = 0; i < n; i++) { int x; scanf("%d", &x); h[i] = x; }
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        AT_NUMB *o = GetDfsOrder4CT(&cg, ct, len, hasH ? h : NULL, n, mode);
        printf("#CASE %d\n", c);
        if (!o) printf("NULL\n");
        else {
            for (int k = 0; ; k += 3) {
                if (!o[k]) break;
                printf(" %d %d %d", o[k], o[k + 1], o[k + 2]);
            }
            printf("\n");
        }
        char s[16];
        int sl;
        scanf("%d", &sl);
        for (int i = 0; i < sl; i++) { int x; scanf("%d", &x); s[i] = x; }
        s[sl] = 0;
        mystrrev(s);
        printf("%s\n", s);
    }
    return 0;
}
"""


def main():
    rng = random.Random(77)
    cases = []
    for _ in range(6000):
        n = rng.randrange(1, 11)
        edges = set()
        mode = rng.choice([0, 0, 8])
        # mostly connected random graphs with rings, sometimes disconnected
        for v in range(1, n):
            if rng.random() < 0.9:
                edges.add((rng.randrange(0, v), v))
        for _e in range(rng.choice([0, 0, 1, 2, 3])):
            if n > 1:
                a, b = sorted(rng.sample(range(n), 2))
                edges.add((a, b))
        nb = {v: sorted(a for a, b in edges if b == v) for v in range(n)}
        ct = []
        for v in range(n):
            ct.append(v + 1)
            ct += [a + 1 for a in nb[v]]
        has_h = rng.choice([0, 1])
        h = [rng.choice([0, 0, 1, 2, 3, 12]) for _ in range(n)]
        sl = rng.randrange(0, 10)
        s = [rng.randrange(65, 91) for _ in range(sl)]
        cases.append(f"{n} {len(ct)} {mode} {has_h}\n" + " ".join(map(str, ct)) + "\n" + " ".join(map(str, h)) + f"\n{sl} " + " ".join(map(str, s)))
    stdin = f"{len(cases)}\n" + "\n".join(cases) + "\n"
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True)
    out = r.stdout.decode().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != "" or True][:2]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_make4.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
