"""Generates tests/data/inchi_can2d.json from the InChI C code: `CanonGraph`
of ichican2.c (the canonical numbering search), run on random and symmetric
graphs with and without the additional layers and with a previous result as
the fixed connection table (compiled with the upstream source included).

    python tools/inchi_can2d_reference.py
"""
import importlib.util
import itertools
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
static void pa(const char *t, AT_RANK *a, int n) { printf("%s", t); for (int i = 0; i < n; i++) printf(" %d", a[i]); printf("\n"); }
static void dump(const char *t, ConTable *c, int n) {
    printf("%s %d %d %d %d %d |", t, c->lenCt, c->lenPos, c->lenNumH, c->len_iso_sort_key, c->len_iso_exchg_atnos);
    for (int i = 0; i <= c->lenCt && i < c->maxlenCt; i++) printf(" %d", c->Ctbl[i]);
    printf(" |"); for (int i = 0; i < c->lenPos; i++) printf(" %d", c->nextCtblPos[i]);
    printf(" |"); for (int i = 0; i < c->lenPos; i++) printf(" %d", c->nextAtRank[i]);
    printf(" |"); if (c->NumH) for (int i = 0; i < c->lenNumH; i++) printf(" %d", c->NumH[i]);
    printf(" |"); if (c->NumHfixed) for (int i = 0; i < n; i++) printf(" %d", c->NumHfixed[i]);
    printf(" |"); if (c->iso_sort_key) for (int i = 0; i < c->len_iso_sort_key; i++) printf(" %ld", c->iso_sort_key[i]);
    printf(" |"); if (c->iso_exchg_atnos) for (int i = 0; i < c->len_iso_exchg_atnos; i++) printf(" %d", c->iso_exchg_atnos[i]);
    printf("\n");
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        int n, m, dig, large, useaux, nst, masks[4];
        scanf("%d %d %d %d %d %d", &n, &m, &dig, &large, &useaux, &nst);
        for (int i = 0; i < nst; i++) scanf("%d", &masks[i]);
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        int col[64], ed[200][2];
        for (int i = 0; i < n; i++) scanf("%d", &col[i]);
        for (int i = 0; i < m; i++) scanf("%d %d", &ed[i][0], &ed[i][1]);
        for (int i = 0; i < m; i++) {
            int a = ed[i][0], b = ed[i][1];
            at[a].neighbor[at[a].valence++] = b;
            at[b].neighbor[at[b].valence++] = a;
        }
        NUM_H nh[64], nhf[64]; AT_ISO_SORT_KEY isk[64]; AT_RANK aux[64];
        for (int i = 0; i < n; i++) { int x; scanf("%d", &x); nh[i] = x; }
        for (int i = 0; i < n; i++) { int x; scanf("%d", &x); nhf[i] = x; }
        for (int i = 0; i < n; i++) { long x; scanf("%ld", &x); isk[i] = x; }
        for (int i = 0; i < n; i++) { int x; scanf("%d", &x); aux[i] = x; }
        printf("N %d\n", c);
        ConTable *fix = NULL;
        for (int st = 0; st < nst; st++) {
            NEIGH_LIST *G = CreateNeighList(n, n, at, 0, NULL);
            Partition *pi = calloc(n + 4, sizeof(Partition));
            for (int i = 0; i < n + 4; i++) PartitionCreate(&pi[i], n + 2);
            AT_RANK *prev = calloc(n + 4, 2);
            for (int i = 0; i < n; i++) { int r = 0; for (int j = 0; j < n; j++) if (col[j] <= col[i]) r++; pi[0].Rank[i] = r; pi[0].AtNumber[i] = i; }
            int nd = 0; { int seen[64] = {0}; for (int i = 0; i < n; i++) if (!seen[pi[0].Rank[i]]++) nd++; }
            long it = 0;
            DifferentiateRanksBasic(&cg, n, G, nd, pi[0].Rank, prev, pi[0].AtNumber, &it, 1);
            CANON_DATA cd; memset(&cd, 0, sizeof(cd));
            cd.nMaxLenLinearCT = 4 * (n + m) + 8; cd.nLenCTAtOnly = n + m;
            cd.LinearCT = calloc(cd.nMaxLenLinearCT + 2, 2);
            if (masks[st] & 1) { cd.NumH = nh; cd.maxlenNumH = n + 2; }
            if (masks[st] & 2) { cd.NumHfixed = nhf; cd.maxlenNumHfixed = n + 2; }
            if (masks[st] & 4) { cd.iso_sort_key = isk; cd.maxlen_iso_sort_key = n + 2; }
            if (useaux) cd.nAuxRank = aux;
            AT_RANK sym[70] = {0}, can[70] = {0}, atn[70] = {0};
            CANON_COUNTS cc; memset(&cc, 0, sizeof(cc));
            ConTable *out = NULL;
            int ret = CanonGraph(NULL, &cg, n, n, n + 2, dig, G, pi, sym, can, atn, &cd, &cc, &fix, &out, large);
            printf("S %d %d %d %ld %ld %ld %ld %ld %.0f %ld %ld\n", st, ret, cd.nLenLinearCT, cc.lNumBreakTies, cc.lNumDecreasedCT, cc.lNumRejectedCT, cc.lNumEqualCT, cc.lNumTotCT, cc.dGroupSize, cc.lNumGenerators, cc.lNumStoredIsomorphisms);
            if (ret >= 0) {
                pa("sym", sym, n); pa("can", can, n); pa("atn", atn, n);
                printf("ct"); for (int i = 0; i < cd.nLenLinearCT; i++) printf(" %d", cd.LinearCT[i]); printf("\n");
                if (out) dump("out", out, n); else printf("out none\n");
            }
            printf("nl"); for (int i = 0; i < n; i++) { printf(" |"); for (int j = 1; j <= G[i][0]; j++) printf(" %d", G[i][j]); } printf("\n");
            if (ret >= 0 && out) fix = out; else fix = NULL;
        }
    }
    return 0;
}
"""


def edges_of(name, rng):
    if name == "cycle":
        n = rng.randrange(3, 13)
        return n, [(i, (i + 1) % n) for i in range(n)]
    if name == "path":
        n = rng.randrange(2, 12)
        return n, [(i, i + 1) for i in range(n - 1)]
    if name == "star":
        n = rng.randrange(3, 8)
        return n, [(0, i) for i in range(1, n)]
    if name == "cube":
        return 8, [(a, a ^ b) for a in range(8) for b in (1, 2, 4) if a < a ^ b]
    if name == "k33":
        return 6, [(i, 3 + j) for i in range(3) for j in range(3)]
    if name == "k44":
        return 8, [(i, 4 + j) for i in range(4) for j in range(4)]
    if name == "k4":
        return 4, [(i, j) for i in range(4) for j in range(i + 1, 4)]
    if name == "petersen":
        e = [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
        return 10, e
    if name == "prism":
        k = rng.randrange(3, 6)
        return 2 * k, [(i, (i + 1) % k) for i in range(k)] + [(k + i, k + (i + 1) % k) for i in range(k)] + [(i, k + i) for i in range(k)]
    if name == "octa":
        return 6, [(i, j) for i in range(6) for j in range(i + 1, 6) if j != (i + 3) % 6]
    if name == "circ":
        n = rng.randrange(5, 13)
        k = rng.randrange(2, n // 2 + 1)
        e = set()
        for i in range(n):
            e.add(tuple(sorted((i, (i + 1) % n))))
            e.add(tuple(sorted((i, (i + k) % n))))
        return n, sorted(e)
    if name == "union":
        a = rng.randrange(3, 6)
        k = rng.randrange(1, 3)
        e = []
        for c in range(k + 1):
            e += [(c * a + i, c * a + (i + 1) % a) for i in range(a)]
        return a * (k + 1), e
    if name == "benzene-like":
        e = [(i, (i + 1) % 6) for i in range(6)] + [(i, 6 + i) for i in range(6)]
        return 12, e
    # random
    n = rng.randrange(1, 11)
    e = set()
    deg = [0] * n
    for _ in range(rng.randrange(0, 2 * n)):
        a, b = rng.randrange(n), rng.randrange(n)
        if a != b and deg[a] < 4 and deg[b] < 4 and (min(a, b), max(a, b)) not in e:
            e.add((min(a, b), max(a, b)))
            deg[a] += 1
            deg[b] += 1
    return n, sorted(e)


NAMES = ["cycle", "path", "star", "cube", "k33", "k44", "k4", "petersen", "prism", "octa", "circ", "union", "benzene-like", "random", "random", "random"]


def gen(rng):
    name = rng.choice(NAMES)
    n, e = edges_of(name, rng)
    perm = list(range(n))
    rng.shuffle(perm)
    e = [(perm[a], perm[b]) if rng.random() < 0.5 else (perm[b], perm[a]) for a, b in e]
    rng.shuffle(e)
    kind = rng.random()
    if kind < 0.5:
        col = [0] * n
    elif kind < 0.8:
        col = [rng.randrange(0, 3) for _ in range(n)]
    else:
        deg = [0] * n
        for a, b in e:
            deg[a] += 1
            deg[b] += 1
        col = deg
    nst = rng.choice([1, 1, 2, 3])
    full = rng.randrange(8)
    masks = sorted({0 if rng.random() < 0.5 else 0, *[full & rng.randrange(8) for _ in range(nst)]})[:nst]
    masks = [0] + [m for m in masks if m][: nst - 1]
    masks = sorted(set(masks)) if len(masks) > 1 else masks
    vals = lambda lo, hi: [rng.randrange(lo, hi) for _ in range(n)]
    sym = rng.random() < 0.6
    return {
        "n": n, "edges": e, "dig": 1 if rng.random() < 0.25 else 0, "large": 1 if rng.random() < 0.01 else 0,
        "aux": 1 if rng.random() < 0.2 else 0, "masks": masks, "col": col,
        "nh": [0] * n if sym else vals(0, 3), "nhf": [0] * n if sym else vals(0, 3),
        "isk": [0] * n if sym else vals(0, 3), "auxr": vals(1, n + 1),
    }


def fmt(c):
    out = [f"{c['n']} {len(c['edges'])} {c['dig']} {c['large']} {c['aux']} {len(c['masks'])}", " ".join(map(str, c["masks"])),
           " ".join(map(str, c["col"]))]
    out += [f"{a} {b}" for a, b in c["edges"]]
    for k in ("nh", "nhf", "isk", "auxr"):
        out.append(" ".join(map(str, c[k])))
    return "\n".join(out) + "\n"


def main():
    rng = random.Random(4242)
    cases = [gen(rng) for _ in range(3000)]
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
    path = ROOT / "tests/data/inchi_can2d.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
