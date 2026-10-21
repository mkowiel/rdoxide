"""Generates tests/data/inchi_ranks.json from the InChI C code: neighbor
lists and rank refinement of ichisort.c / ichimap2.c (`CreateNeighList`,
`CreateNeighListFromLinearCT`, `SortedEquInfoToRanks`,
`SetNewRanksFromNeighLists*`, `DifferentiateRanks*`, `BreakAllTies`,
`CompAtomInvariants2`), linked from the upstream library.

    python tools/inchi_ranks_reference.py
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
int SortNeighLists3( int num_atoms, AT_RANK *nRank, NEIGH_LIST *NeighList, AT_RANK *nAtomNumber );
int CreateCompositeNormAtom() { return 0; }
static NEIGH_LIST *dup(NEIGH_LIST *nl, int n) {
    NEIGH_LIST *r = calloc(n + 1, sizeof(NEIGH_LIST));
    for (int i = 0; i < n; i++) { r[i] = malloc((nl[i][0] + 1) * sizeof(AT_NUMB)); memcpy(r[i], nl[i], (nl[i][0] + 1) * sizeof(AT_NUMB)); }
    return r;
}
static void pnl(NEIGH_LIST *nl, int n) {
    for (int i = 0; i < n; i++) { printf(" ["); for (int j = 0; j <= nl[i][0]; j++) printf("%d%s", nl[i][j], j < nl[i][0] ? "," : ""); printf("]"); }
    printf("\n");
}
static void pa(const char *t, AT_RANK *a, int n) { printf("%s", t); for (int i = 0; i < n; i++) printf(" %d", a[i]); printf("\n"); }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        int n, ntg;
        scanf("%d %d", &n, &ntg);
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        AT_RANK *rank = calloc(n + 8, 2), *anum = calloc(n + 8, 2), *symm = calloc(n + 8, 2);
        for (int i = 0; i < n; i++) {
            int v, r, an, ep, sy;
            scanf("%d %d %d %d %d", &v, &r, &an, &ep, &sy);
            at[i].valence = v; at[i].endpoint = ep; rank[i] = r; anum[i] = an; symm[i] = sy;
            for (int j = 0; j < v; j++) { int q, bt; scanf("%d %d", &q, &bt); at[i].neighbor[j] = q; at[i].bond_type[j] = bt; }
        }
        T_GROUP_INFO tgi;
        memset(&tgi, 0, sizeof(tgi));
        T_GROUP *tg = calloc(ntg + 1, sizeof(T_GROUP));
        AT_NUMB *ends = calloc(n + 1, sizeof(AT_NUMB));
        int pos = 0;
        for (int g = 0; g < ntg; g++) {
            int ne; scanf("%d", &ne);
            tg[g].nNumEndpoints = ne; tg[g].nFirstEndpointAtNoPos = pos;
            for (int k = 0; k < ne; k++) { int e; scanf("%d", &e); ends[pos++] = e; }
        }
        tgi.t_group = tg; tgi.num_t_groups = ntg; tgi.nEndpointAtomNumber = ends;
        int nct; scanf("%d", &nct);
        AT_NUMB *lct = calloc(nct + 1, sizeof(AT_NUMB));
        for (int i = 0; i < nct; i++) { int x; scanf("%d", &x); lct[i] = x; }
        int ninv; scanf("%d", &ninv);
        ATOM_INVARIANT2 *inv = calloc(ninv + 1, sizeof(ATOM_INVARIANT2));
        for (int i = 0; i < ninv; i++) {
            for (int j = 0; j < AT_INV_LENGTH; j++) { int x; scanf("%d", &x); inv[i].val[j] = x; }
            long long k; int ak; scanf("%lld %d", &k, &ak); inv[i].iso_sort_key = k; inv[i].iso_aux_key = ak;
        }
        printf("N %d\n", c);
        for (int sq = 0; sq < 2; sq++) {
            NEIGH_LIST *nl = CreateNeighList(n, n, at, sq, &tgi);
            printf("C%d", sq); pnl(nl, n);
            if (ntg) { NEIGH_LIST *nt = CreateNeighList(n, n + ntg, at, sq, &tgi); printf("T%d", sq); pnl(nt, n + ntg); }
        }
        {
            NEIGH_LIST *nl = CreateNeighListFromLinearCT(lct, nct, n);
            if (nl) { printf("L"); pnl(nl, n); } else printf("L null\n");
        }
        NEIGH_LIST *base = CreateNeighList(n, n, at, 0, &tgi);
        /* equ info <-> ranks */
        {
            AT_RANK *r2 = calloc(n + 8, 2); int ch = -1;
            printf("E %d", SortedEquInfoToRanks(symm, r2, anum, n, &ch)); printf(" %d\n", ch); pa("e", r2, n);
            AT_RANK *s2 = calloc(n + 8, 2);
            printf("Q %d\n", n ? SortedRanksToEquInfo(s2, rank, anum, n) : 0); pa("q", s2, n);
            AT_RANK *inpl = calloc(n + 8, 2); memcpy(inpl, symm, n * 2);
            printf("I %d", SortedEquInfoToRanks(inpl, inpl, anum, n, NULL)); pa("i", inpl, n);
        }
        /* SortNeighLists and SetNewRanks* */
        {
            NEIGH_LIST *nl = dup(base, n);
            AT_RANK *an = calloc(n + 8, 2); memcpy(an, anum, n * 2);
            SortNeighLists2(n, rank, nl, an); printf("S2"); pnl(nl, n);
            nl = dup(base, n);
            SortNeighLists3(n, rank, nl, an); printf("S3"); pnl(nl, n);
            SortNeighListsBySymmAndCanonRank(n, nl, symm, rank); printf("SB"); pnl(nl, n);
        }
        for (int v = 0; v < 4; v++) {
            NEIGH_LIST *nl = dup(base, n);
            AT_RANK *an = calloc(n + 8, 2); memcpy(an, anum, n * 2);
            AT_RANK *nr = calloc(n + 8, 2);
            int r;
            if (v == 0) r = SetNewRanksFromNeighLists3(&cg, n, nl, rank, nr, an);
            else if (v == 1) r = SetNewRanksFromNeighLists4(&cg, n, nl, rank, nr, an, (AT_RANK) (n / 2 + 1));
            else r = SetNewRanksFromNeighLists(&cg, n, nl, rank, nr, an, v - 2, CompNeighListRanksOrd);
            printf("R%d %d", v, r); pa("", nr, n); pa("a", an, n);
        }
        for (int alt = 0; alt < 2; alt++) {
            NEIGH_LIST *nl = dup(base, n);
            AT_RANK *an = calloc(n + 8, 2); memcpy(an, anum, n * 2);
            AT_RANK *nr = calloc(n + 8, 2);
            printf("RB%d %d", alt, SetNewRanksFromNeighLists(&cg, n, nl, rank, nr, an, alt, CompNeighListRanks)); pa("", nr, n);
        }
        /* Differentiate */
        for (int v = 0; v < 6; v++) {
            NEIGH_LIST *nl = dup(base, n);
            AT_RANK *an = calloc(n + 8, 2); memcpy(an, anum, n * 2);
            AT_RANK *cur = calloc(n + 8, 2), *prev = calloc(n + 8, 2);
            memcpy(cur, rank, n * 2);
            long it = 0;
            int r;
            if (v == 0) r = DifferentiateRanks2(&cg, n, nl, 1, cur, prev, an, &it, 0);
            else if (v == 1) r = DifferentiateRanks2(&cg, n, nl, 1, cur, prev, an, &it, 1);
            else if (v == 2) r = DifferentiateRanks3(&cg, n, nl, 1, cur, prev, an, &it);
            else if (v == 3) r = DifferentiateRanks4(&cg, n, nl, 1, cur, prev, an, (AT_RANK) (n / 2 + 1), &it);
            else r = DifferentiateRanksBasic(&cg, n, nl, 1, cur, prev, an, &it, v - 4);
            printf("D%d %d %ld", v, r, it); pa("", cur, n); pa("p", prev, n); pa("a", an, n); printf("n"); pnl(nl, n);
        }
        /* BreakAllTies */
        {
            NEIGH_LIST *nl = dup(base, n);
            AT_RANK *stack[4]; stack[0] = calloc(n + 8, 2); stack[1] = calloc(n + 8, 2); stack[2] = NULL; stack[3] = NULL;
            memcpy(stack[0], rank, n * 2); memcpy(stack[1], anum, n * 2);
            AT_RANK *tmp = calloc(n + 8, 2);
            CANON_STAT cs; memset(&cs, 0, sizeof(cs));
            int r = BreakAllTies(&cg, n, n + 4, stack, nl, tmp, &cs);
            printf("B %d %ld %ld", r, cs.lNumNeighListIter, cs.lNumBreakTies); pa("", stack[2], n); pa("b", stack[3], n);
        }
        /* atom invariants */
        {
            cg.m_pAtomInvariant2ForSort = inv;
            printf("V");
            for (int i = 0; i < ninv; i++) for (int j = 0; j < ninv; j++) {
                AT_RANK a = i, b = j;
                printf(" %d,%d", CompAtomInvariants2Only(&a, &b, &cg), CompAtomInvariants2(&a, &b, &cg));
            }
            printf("\n");
        }
    }
    return 0;
}
"""


def gen(rng):
    n = rng.randrange(1, 13)
    adj = [[] for _ in range(n)]
    for i in range(1, n):
        if rng.random() < 0.9:
            p = rng.randrange(i)
            adj[i].append(p)
            adj[p].append(i)
    for _ in range(rng.randrange(0, n)):
        a, b = rng.randrange(n), rng.randrange(n)
        if a != b and b not in adj[a]:
            adj[a].append(b)
            adj[b].append(a)
    for a in adj:
        rng.shuffle(a)
    adj = [a[:12] for a in adj]
    # initial colours -> ranks (rank = number of atoms with colour <= own)
    ncol = rng.choice([1, 2, 3, n])
    col = [rng.randrange(ncol) for _ in range(n)]
    rank = [sum(1 for c in col if c <= col[i]) for i in range(n)]
    anum = sorted(range(n), key=lambda i: (rank[i], i))
    symm_col = [rng.randrange(max(1, ncol)) for _ in range(n)]
    symm = [sum(1 for c in symm_col if c <= symm_col[i]) for i in range(n)]
    anum_sym = sorted(range(n), key=lambda i: (symm_col[i], i)) if rng.random() < 0.5 else anum
    ntg = rng.choice([0, 0, 1, 2]) if n >= 2 else 0
    groups = []
    endpoint = [0] * n
    free = list(range(n))
    rng.shuffle(free)
    for g in range(ntg):
        k = min(len(free), rng.randrange(1, 4))
        eps = [free.pop() for _ in range(k)]
        if not eps:
            break
        groups.append(sorted(eps))
        for e in eps:
            endpoint[e] = g + 1
    ntg = len(groups)
    bonds = {}
    for i in range(n):
        for q in adj[i]:
            bonds[(min(i, q), max(i, q))] = bonds.get((min(i, q), max(i, q)), rng.choice([1, 1, 2, 4]))
    # linear CT (canonical numbers = a permutation) sometimes corrupted
    perm = list(range(1, n + 1))
    rng.shuffle(perm)
    lct = []
    for v in range(1, n + 1):
        i = perm.index(v)
        nbrs = sorted(perm[q] for q in adj[i] if perm[q] < v)
        lct.append(v)
        lct.extend(nbrs)
    mode = rng.random()
    if mode < 0.15 and len(lct) > 2:
        lct[rng.randrange(1, len(lct))] = rng.randrange(1, n + 3)
    elif mode < 0.2 and len(lct) > 1:
        lct = lct[:-1] or [1]
    ninv = rng.randrange(1, 5)
    inv = []
    for _ in range(ninv):
        base = [rng.randrange(3) for _ in range(10)]
        inv.append({"val": base, "k": rng.choice([0, 0, 1, -5, 10 ** 12]), "ak": rng.choice([0, 0, 1, -1])})
    for i in range(1, ninv):
        if rng.random() < 0.5:
            inv[i]["val"] = list(inv[i - 1]["val"][:7]) + [rng.randrange(3) for _ in range(3)]
    return {"n": n, "ntg": ntg, "adj": adj, "bonds": [[bonds[(min(i, q), max(i, q))] for q in adj[i]] for i in range(n)],
            "rank": rank, "anum": anum_sym if False else anum, "symm": symm, "endpoint": endpoint, "groups": groups,
            "lct": lct, "inv": inv}


def fmt(c):
    n = c["n"]
    t = f"{n} {c['ntg']}\n"
    for i in range(n):
        t += f"{len(c['adj'][i])} {c['rank'][i]} {c['anum'][i]} {c['endpoint'][i]} {c['symm'][i]} " + " ".join(f"{q} {b}" for q, b in zip(c["adj"][i], c["bonds"][i])) + "\n"
    for g in c["groups"]:
        t += f"{len(g)} " + " ".join(map(str, g)) + "\n"
    t += f"{len(c['lct'])} " + " ".join(map(str, c["lct"])) + "\n"
    t += f"{len(c['inv'])}\n"
    for v in c["inv"]:
        t += " ".join(map(str, v["val"])) + f" {v['k']} {v['ak']}\n"
    return t


def main():
    rng = random.Random(101)
    cases = [gen(rng) for _ in range(3000)]
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True)
        assert r.returncode == 0, (r.returncode, r.stderr)
        out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_ranks.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
