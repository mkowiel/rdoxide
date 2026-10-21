"""Generates tests/data/inchi_knownsc.json from the InChI C code: the stereo
center parts of ichicans.c that work on parities known in advance
(`SetKnownStereoCenterParities`, `MarkKnownEqualStereoCenterParities`,
`RemoveKnownNonStereoCenterParities`) and the whole `FillOutStereoParities`,
on random, partly mirror-symmetric molecules with real equivalence ranks.

    python tools/inchi_knownsc_reference.py
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


ks = load("inchi_knownsb_reference")
ma = ks.ma
SRC = ma.SRC

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mode.h"
#include "ichicant.h"
#include "ichicomn.h"
#include "ichister.h"
int SetKnownStereoCenterParities( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nRank, const AT_RANK *nAtomNumber );
int MarkKnownEqualStereoCenterParities( sp_ATOM *at, int num_atoms, const AT_RANK *nRank, const AT_RANK *nAtomNumber );
int RemoveKnownNonStereoCenterParities( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nRank, CANON_STAT *pCS );
int FillAllStereoDescriptors( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nAtomNumberCanon, CANON_STAT *pCS );
int FillOutStereoParities( sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nAtomNumberCanon, const AT_RANK *nRank, const AT_RANK *nAtomNumber, CANON_STAT *pCS, CANON_GLOBALS *pCG, int bIsotopic );
int CreateCompositeNormAtom() { return 0; }
static void dump(sp_ATOM *at, int n) {
    for (int i = 0; i < n; i++) {
        printf(" |%d %d %d %d", at[i].parity, at[i].stereo_atom_parity, at[i].final_parity, at[i].bHasStereoOrEquToStereo);
        for (int j = 0; j < 3; j++) printf(" %d %d %d %d", at[i].stereo_bond_neighbor[j], at[i].stereo_bond_ord[j], at[i].stereo_bond_z_prod[j], at[i].stereo_bond_parity[j]);
    }
    printf("\n");
}
static void dumpcs(CANON_STAT *s) {
    printf("L %d %d", s->nLenLinearCTStereoCarb, s->nLenLinearCTStereoDble);
    for (int i = 0; i < s->nLenLinearCTStereoCarb; i++) printf(" c%d,%d", s->LinearCTStereoCarb[i].at_num, s->LinearCTStereoCarb[i].parity);
    for (int i = 0; i < s->nLenLinearCTStereoDble; i++) printf(" d%d,%d,%d", s->LinearCTStereoDble[i].at_num1, s->LinearCTStereoDble[i].at_num2, s->LinearCTStereoDble[i].parity);
    printf("\n");
}
static void initcs(CANON_STAT *s) {
    memset(s, 0, sizeof(*s));
    s->nMaxLenLinearCTStereoCarb = 16; s->nMaxLenLinearCTStereoDble = 16;
    s->LinearCTStereoCarb = calloc(17, sizeof(AT_STEREO_CARB));
    s->LinearCTStereoDble = calloc(17, sizeof(AT_STEREO_DBLE));
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    CANON_GLOBALS cg;
    memset(&cg, 0, sizeof(cg));
    for (int c = 0; c < ncases; c++) {
        int n, iso;
        scanf("%d %d", &n, &iso);
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        AT_RANK *rank = calloc(n + 1, sizeof(AT_RANK)), *crank = calloc(n + 1, sizeof(AT_RANK));
        AT_RANK *anum = calloc(n + 1, sizeof(AT_RANK)), *acanon = calloc(n + 1, sizeof(AT_RANK));
        for (int i = 0; i < n; i++) {
            int v, p, sap, r, cr, nh, rs, an, ac;
            scanf("%d %d %d %d %d %d %d %d %d", &v, &p, &sap, &r, &cr, &nh, &rs, &an, &ac);
            at[i].valence = v; at[i].parity = p; at[i].stereo_atom_parity = sap; at[i].num_H = nh; at[i].nRingSystem = rs;
            rank[i] = r; crank[i] = cr; anum[i] = an; acanon[i] = ac;
            for (int j = 0; j < v; j++) { int q; scanf("%d", &q); at[i].neighbor[j] = q; }
            for (int j = 0; j < 3; j++) { int a, b, z, pp; scanf("%d %d %d %d", &a, &b, &z, &pp); at[i].stereo_bond_neighbor[j] = a; at[i].stereo_bond_ord[j] = b; at[i].stereo_bond_z_prod[j] = z; at[i].stereo_bond_parity[j] = pp; }
        }
        printf("N %d\n", c);
        sp_ATOM *w = malloc((n + 1) * sizeof(sp_ATOM));
        memcpy(w, at, (n + 1) * sizeof(sp_ATOM));
        printf("K %d", SetKnownStereoCenterParities(&cg, w, n, crank, rank, anum)); dump(w, n);
        printf("M %d", MarkKnownEqualStereoCenterParities(w, n, rank, anum)); dump(w, n);
        CANON_STAT s;
        initcs(&s);
        printf("F %d\n", FillAllStereoDescriptors(&cg, w, n, crank, acanon, &s));
        printf("R %d", RemoveKnownNonStereoCenterParities(&cg, w, n, crank, rank, &s)); dump(w, n);
        dumpcs(&s);
        free(s.LinearCTStereoCarb); free(s.LinearCTStereoDble);
        memcpy(w, at, (n + 1) * sizeof(sp_ATOM));
        initcs(&s);
        printf("O %d", FillOutStereoParities(w, n, crank, acanon, rank, anum, &s, &cg, iso)); dump(w, n);
        dumpcs(&s);
        free(s.LinearCTStereoCarb); free(s.LinearCTStereoDble);
        free(w); free(at); free(rank); free(crank); free(anum); free(acanon);
    }
    return 0;
}
"""


def build(rng):
    adj = []

    def new():
        adj.append([])
        return len(adj) - 1

    def link(a, b):
        adj[a].append(b)
        adj[b].append(a)

    def tree(d):
        return [tree(d + 1) for _ in range(rng.choice([0, 1, 2, 3] if d < 2 else [0, 1, 2]))] if d < 3 else []

    def grow(parent, kids):
        for sub in kids:
            c = new()
            link(parent, c)
            grow(c, sub)
    half = [tree(2) for _ in range(rng.choice([2, 3]))]
    r1 = new()
    grow(r1, half)
    r2 = new()
    grow(r2, half if rng.random() < 0.8 else [tree(2) for _ in range(rng.choice([2, 3]))])
    link(r1, r2)
    if rng.random() < 0.3:
        # a lone third center connected to one of the roots
        a = new()
        link(a, rng.choice([r1, r2]))
        for _ in range(rng.choice([2, 3])):
            c = new()
            link(a, c)
    return adj


def gen(rng):
    adj = build(rng)
    n = len(adj)
    if n > 24 or any(len(a) > 4 for a in adj):
        return gen(rng)
    atoms = []
    for i in range(n):
        v = len(adj[i])
        par = rng.choice([1, 2]) if v >= 3 and rng.random() < 0.8 else rng.choice([0, 0, 3, 4])
        sap = rng.choice([6, 6, 6, 6, 0, 1, 2, 6 | 0x40]) if par else rng.choice([0, 0, 6])
        atoms.append({"nb": list(adj[i]), "sb": [[0, 0, 0, 0] for _ in range(3)], "par": par, "sap": sap,
                      "nh": rng.choice([0, 0, 1]), "rs": i + 1})
    colors = [(len(a["nb"]), a["par"] != 0, a["sap"] != 0) for a in atoms]
    order = sorted(set(colors))
    rank = ks.wl_ranks(adj, [order.index(c) for c in colors])
    idx = sorted(range(n), key=lambda i: (rank[i], rng.random()))
    crank = [0] * n
    for pos, i in enumerate(idx):
        crank[i] = pos + 1
    anum = sorted(range(n), key=lambda i: (rank[i], i))
    acanon = sorted(range(n), key=lambda i: crank[i])
    for i, a in enumerate(atoms):
        a["rank"], a["crank"], a["an"], a["ac"] = rank[i], crank[i], anum[i], acanon[i]
    return {"atoms": atoms, "iso": rng.choice([0, 1])}


def fmt(c):
    return f"{len(c['atoms'])} {c['iso']}\n" + ks.fmt(c).split("\n", 1)[1]


def main():
    rng = random.Random(91)
    cases = [gen(rng) for _ in range(6000)]
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True)
        assert r.returncode == 0, r.stderr
        out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_knownsc.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
