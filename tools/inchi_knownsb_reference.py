"""Generates tests/data/inchi_knownsb.json from the InChI C code: the
stereo bond parts of ichicans.c that work on parities known in advance
(`SetKnownStereoBondParities`, `MarkKnownEqualStereoBondParities`,
`RemoveKnownNonStereoBondParities` with `PathsHaveIdenticalKnownParities`)
on random, partly mirror-symmetric molecules with real equivalence ranks.

    python tools/inchi_knownsb_reference.py
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
#include "ichicant.h"
#include "ichicomn.h"
#include "ichister.h"
int SetKnownStereoBondParities( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nRank, const AT_RANK *nAtomNumber );
int MarkKnownEqualStereoBondParities( sp_ATOM *at, int num_atoms, const AT_RANK *nRank, const AT_RANK *nAtomNumber );
int RemoveKnownNonStereoBondParities( sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nRank, CANON_STAT *pCS );
int FillAllStereoDescriptors( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nAtomNumberCanon, CANON_STAT *pCS );
int CreateCompositeNormAtom() { return 0; }
static void dump(sp_ATOM *at, int n) {
    for (int i = 0; i < n; i++) {
        printf(" |%d %d %d", at[i].parity, at[i].stereo_atom_parity, at[i].bHasStereoOrEquToStereo);
        for (int j = 0; j < 3; j++) printf(" %d %d %d %d", at[i].stereo_bond_neighbor[j], at[i].stereo_bond_ord[j], at[i].stereo_bond_z_prod[j], at[i].stereo_bond_parity[j]);
    }
    printf("\n");
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    CANON_GLOBALS cg;
    memset(&cg, 0, sizeof(cg));
    for (int c = 0; c < ncases; c++) {
        int n;
        scanf("%d", &n);
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
        printf("K %d", SetKnownStereoBondParities(&cg, at, n, crank, rank, anum)); dump(at, n);
        printf("M %d", MarkKnownEqualStereoBondParities(at, n, rank, anum)); dump(at, n);
        CANON_STAT s;
        memset(&s, 0, sizeof(s));
        s.nMaxLenLinearCTStereoCarb = 16; s.nMaxLenLinearCTStereoDble = 16;
        s.LinearCTStereoCarb = calloc(17, sizeof(AT_STEREO_CARB));
        s.LinearCTStereoDble = calloc(17, sizeof(AT_STEREO_DBLE));
        printf("F %d %d\n", FillAllStereoDescriptors(&cg, at, n, crank, acanon, &s), s.nLenLinearCTStereoDble);
        printf("R %d", RemoveKnownNonStereoBondParities(at, n, crank, rank, &s)); dump(at, n);
        printf("L %d", s.nLenLinearCTStereoDble);
        for (int i = 0; i < s.nLenLinearCTStereoDble; i++) printf(" d%d,%d,%d", s.LinearCTStereoDble[i].at_num1, s.LinearCTStereoDble[i].at_num2, s.LinearCTStereoDble[i].parity);
        printf("\n");
        free(s.LinearCTStereoCarb); free(s.LinearCTStereoDble);
        free(at); free(rank); free(crank); free(anum); free(acanon);
    }
    return 0;
}
"""


def half_tree(rng, depth, maxdeg):
    """list of (parent, ...) edges as nested structure: returns list of children counts"""
    def rec(d):
        k = rng.choice([0, 1, 2]) if d < depth else 0
        return [rec(d + 1) for _ in range(k)]
    return [rec(1) for _ in range(rng.choice([1, 2]))]


def build(rng):
    """returns (adjacency lists, list of stereo bonds [(a, b, chain_atoms)])"""
    adj = []
    def new():
        adj.append([])
        return len(adj) - 1
    def link(a, b):
        adj[a].append(b)
        adj[b].append(a)
    def grow(parent, kids):
        for sub in kids:
            c = new()
            link(parent, c)
            grow(c, sub)
    mirror = rng.random() < 0.7
    tree = half_tree(rng, 2, 3)
    r1 = new()
    grow(r1, tree)
    sbs = []
    chain = rng.choice([0, 0, 0, 1, 2])
    if mirror:
        r2 = new()
        grow(r2, tree if rng.random() < 0.85 else half_tree(rng, 2, 3))
    else:
        r2 = new()
        grow(r2, half_tree(rng, 2, 3))
    prev = r1
    mids = []
    for _ in range(chain):
        m = new()
        link(prev, m)
        mids.append(m)
        prev = m
    link(prev, r2)
    sbs.append((r1, r2, mids))
    # an extra random stereo bond sometimes
    if rng.random() < 0.3:
        a = rng.randrange(len(adj))
        if len(adj[a]) < 3:
            b = new()
            link(a, b)
            c = new()
            link(b, c)
            if rng.random() < 0.6:
                d = new(); link(c, d)
            sbs.append((b, c, []))
    return adj, sbs


def wl_ranks(adj, colors):
    n = len(adj)
    col = list(colors)
    for _ in range(n):
        sig = [(col[i], tuple(sorted(col[j] for j in adj[i]))) for i in range(n)]
        order = sorted(set(sig))
        new = [order.index(s) for s in sig]
        if len(set(new)) == len(set(col)):
            col = new
            break
        col = new
    # rank = number of atoms with colour <= own colour
    return [sum(1 for c in col if c <= col[i]) for i in range(n)]


def gen(rng):
    adj, sbs = build(rng)
    n = len(adj)
    if n > 20 or any(len(a) > 3 for a in adj):
        return gen(rng)
    atoms = [{"nb": list(adj[i]), "sb": [[0, 0, 0, 0] for _ in range(3)], "par": 0, "sap": rng.choice([0] * 6 + [6, 1]),
              "nh": rng.choice([0, 0, 1]), "rs": i + 1} for i in range(n)]
    for (a, b, mids) in sbs:
        first_a = mids[0] if mids else b
        first_b = mids[-1] if mids else a
        bits = len(mids) * 8
        par = rng.choice([6, 6, 6, 6, 1, 2, 3, 4, 6 | 0x40, 1 | 0x40]) | bits
        zp = rng.choice([100, -100, 100, 60, 20, -20])
        pa = rng.choice([1, 2])
        pb = rng.choice([1, 2]) if rng.random() < 0.7 else pa
        for (x, y, first, p, z) in ((a, b, first_a, pa, zp), (b, a, first_b, pb, zp * rng.choice([1, 1, -1]) if False else zp)):
            slot = next(s for s in range(3) if atoms[x]["sb"][s][0] == 0)
            atoms[x]["sb"][slot] = [y + 1, atoms[x]["nb"].index(first), z, par]
            atoms[x]["par"] = p if rng.random() < 0.9 else rng.choice([0, 3, 4])
    colors = [(len(a["nb"]), a["par"] != 0, a["sap"] != 0, tuple(sorted(s[3] & 0x3f for s in a["sb"] if s[0]))) for a in atoms]
    order = sorted(set(colors))
    rank = wl_ranks(adj, [order.index(c) for c in colors])
    idx = sorted(range(n), key=lambda i: (rank[i], rng.random()))
    crank = [0] * n
    for pos, i in enumerate(idx):
        crank[i] = pos + 1
    anum = sorted(range(n), key=lambda i: (rank[i], i))
    acanon = sorted(range(n), key=lambda i: crank[i])
    for i, a in enumerate(atoms):
        a["rank"], a["crank"], a["an"], a["ac"] = rank[i], crank[i], anum[i], acanon[i]
    return {"atoms": atoms}


def fmt(c):
    t = f"{len(c['atoms'])}\n"
    for a in c["atoms"]:
        t += f"{len(a['nb'])} {a['par']} {a['sap']} {a['rank']} {a['crank']} {a['nh']} {a['rs']} {a['an']} {a['ac']} " + " ".join(map(str, a["nb"])) + "\n"
        t += " ".join(" ".join(map(str, s)) for s in a["sb"]) + "\n"
    return t


def main():
    rng = random.Random(81)
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
    path = ROOT / "tests/data/inchi_knownsb.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
