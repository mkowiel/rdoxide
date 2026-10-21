"""Generates tests/data/inchi_can2e.json from the InChI C code:
`FixCanonEquivalenceInfo` (ichicano.c) and `make_iso_sort_key` (ichiisot.c),
helpers of the base canonical numbering.

    python tools/inchi_can2e_reference.py
"""
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"tools/{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ma = load("inchi_markalt_reference")
SRC = ma.SRC

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ichican2.c"
int CreateCompositeNormAtom() { return 0; }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        int n;
        scanf("%d", &n);
        AT_RANK sym[64] = {0}, cur[64] = {0}, tmp[64] = {0}, an[64] = {0};
        for (int i = 0; i < n; i++) { int x; scanf("%d", &x); sym[i] = x; }
        for (int i = 0; i < n; i++) { int x; scanf("%d", &x); cur[i] = x; }
        for (int i = 0; i < n; i++) { int x; scanf("%d", &x); an[i] = x; }
        int d[4]; scanf("%d %d %d %d", &d[0], &d[1], &d[2], &d[3]);
        int bch = -1;
        int r = FixCanonEquivalenceInfo(&cg, n, sym, cur, tmp, an, &bch);
        printf("N %d\n", c);
        printf("F %d %d |", r, bch);
        for (int i = 0; i < n; i++) printf(" %d", sym[i]);
        printf(" |"); for (int i = 0; i < n; i++) printf(" %d", cur[i]);
        printf(" |"); for (int i = 0; i < n; i++) printf(" %d", tmp[i]);
        printf(" |"); for (int i = 0; i < n; i++) printf(" %d", an[i]);
        printf("\nK %ld\n", (long) make_iso_sort_key(d[0], d[1], d[2], d[3]));
    }
    return 0;
}
"""


def gen(rng):
    n = rng.randrange(1, 14)
    # classes: sizes, symm rank = first rank of class (1-based); atoms permuted
    order = list(range(n))
    rng.shuffle(order)
    sizes = []
    left = n
    while left:
        s = rng.randrange(1, left + 1)
        sizes.append(s)
        left -= s
    sym = [0] * n
    cur = [0] * n
    pos = 0
    for s in sizes:
        for j in range(s):
            sym[order[pos + j]] = pos + 1
            cur[order[pos + j]] = pos + s
        pos += s
    if rng.random() < 0.4:
        # coarser/finer current ranks
        k = rng.randrange(1, n + 1)
        cur = [rng.randrange(1, k + 1) if rng.random() < 0.5 else v for v in cur]
    if rng.random() < 0.15:
        i = rng.randrange(n)
        sym[i] = rng.randrange(1, n + 1)
    an = list(range(n))
    rng.shuffle(an)
    d = [rng.randrange(-5, 40), rng.randrange(0, 5), rng.randrange(0, 5), rng.randrange(0, 5)]
    return {"n": n, "sym": sym, "cur": cur, "an": an, "d": d}


def fmt(c):
    return "\n".join([str(c["n"]), " ".join(map(str, c["sym"])), " ".join(map(str, c["cur"])),
                      " ".join(map(str, c["an"])), " ".join(map(str, c["d"]))]) + "\n"


def main():
    rng = random.Random(777)
    cases = [gen(rng) for _ in range(6000)]
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, timeout=900)
        assert r.returncode == 0, (r.returncode, r.stderr.decode()[:2000])
        out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_can2e.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
