"""Generates tests/data/inchi_mapbreak.json from the InChI C code: the
tie breaking and symmetric path functions of ichimap2.c
(`BreakNeighborsTie`, `CheckNextSymmNeighborsAndBonds`,
`CreateCheckSymmPaths`), linked from the upstream
library.

    python tools/inchi_mapbreak_reference.py
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


ms = load("inchi_mapstereo_reference")
ma = ms.ma
SRC = ma.SRC

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mode.h"
#include "ichicant.h"
#include "ichicomn.h"
#include "ichister.h"
int CreateCompositeNormAtom() { return 0; }
int BreakNeighborsTie( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, int num_at_tg, int ib, int ia, AT_RANK *neigh_num, int in1, int in2, int mode, AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList, const AT_RANK *nSymmRank, AT_RANK *nCanonRank, NEIGH_LIST *nl1, NEIGH_LIST *nl2, long *lNumIter );
int CheckNextSymmNeighborsAndBonds( sp_ATOM *at, AT_RANK cur1, AT_RANK cur2, AT_RANK n1, AT_RANK n2, AT_RANK *nAvoidCheckAtom, AT_RANK *nVisited1, AT_RANK *nVisited2, AT_RANK *nVisitOrd1, AT_RANK *nVisitOrd2, const AT_RANK *nRank1, const AT_RANK *nRank2 );
int CreateCheckSymmPaths( sp_ATOM *at, AT_RANK prev1, AT_RANK cur1, AT_RANK prev2, AT_RANK cur2, AT_RANK *nAvoidCheckAtom, AT_RANK *nVisited1, AT_RANK *nVisited2, AT_RANK *nVisitOrd1, AT_RANK *nVisitOrd2, NEIGH_LIST *nl1, NEIGH_LIST *nl2, const AT_RANK *nRank1, const AT_RANK *nRank2, AT_RANK *nCanonRank, AT_RANK *nLength, int *bParitiesInverted, int mode );
static void pa(AT_RANK *a, int n) { for (int i = 0; i < n; i++) printf(" %d", a[i]); printf(" |"); }
static void pnl(NEIGH_LIST *nl, int n) { for (int i = 0; i < n; i++) { printf(" ["); for (int j = 0; j <= nl[i][0]; j++) printf("%d%s", nl[i][j], j < nl[i][0] ? "," : ""); printf("]"); } printf("\n"); }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        int n;
        scanf("%d", &n);
        int nmax = n + 2;
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        AT_RANK *crank = calloc(nmax + 8, 2), *r1 = calloc(nmax + 8, 2), *r2 = calloc(nmax + 8, 2);
        for (int i = 0; i < n; i++) {
            int v, p, sap, r, cr, nh, an, ac, rr2, aa2, rs;
            scanf("%d %d %d %d %d %d %d %d %d %d %d", &v, &p, &sap, &r, &cr, &nh, &an, &ac, &rr2, &aa2, &rs);
            at[i].valence = v; at[i].parity = p; at[i].stereo_atom_parity = sap; at[i].num_H = nh; at[i].nRingSystem = rs;
            r1[i] = r; r2[i] = rr2; crank[i] = cr;
            for (int j = 0; j < v; j++) { int q; scanf("%d", &q); at[i].neighbor[j] = q; }
            for (int j = 0; j < 3; j++) { int a, b, z, pp; scanf("%d %d %d %d", &a, &b, &z, &pp); at[i].stereo_bond_neighbor[j] = a; at[i].stereo_bond_ord[j] = b; at[i].stereo_bond_z_prod[j] = z; at[i].stereo_bond_parity[j] = pp; }
        }
        printf("N %d\n", c);
        int nc; scanf("%d", &nc);
        for (int k = 0; k < nc; k++) {
            int cur1, cur2, n1, n2, a0, a1, same;
            scanf("%d %d %d %d %d %d %d", &cur1, &cur2, &n1, &n2, &a0, &a1, &same);
            AT_RANK av[2] = {a0, a1};
            AT_RANK *v1 = calloc(nmax, 2), *v2 = calloc(nmax, 2), *o1 = calloc(nmax, 2), *o2 = calloc(nmax, 2);
            for (int i = 0; i < n; i++) { int x; scanf("%d", &x); v1[i] = x; }
            for (int i = 0; i < n; i++) { int x; scanf("%d", &x); v2[i] = x; }
            for (int i = 0; i < n; i++) { int x; scanf("%d", &x); o1[i] = x; }
            for (int i = 0; i < n; i++) { int x; scanf("%d", &x); o2[i] = x; }
            printf("C %d\n", CheckNextSymmNeighborsAndBonds(at, cur1, cur2, n1, n2, av, v1, v2, o1, o2, r1, same ? r1 : r2));
        }
        int np; scanf("%d", &np);
        NEIGH_LIST *nlp1 = CreateNeighList(n, n, at, 0, NULL), *nlp2 = CreateNeighList(n, n, at, 0, NULL);
        for (int k = 0; k < np; k++) {
            int p1, c1, p2, c2, a0, a1, inv, mode, same;
            scanf("%d %d %d %d %d %d %d %d %d", &p1, &c1, &p2, &c2, &a0, &a1, &inv, &mode, &same);
            AT_RANK av[2] = {a0, a1};
            AT_RANK *v1 = calloc(nmax, 2), *v2 = calloc(nmax, 2), *o1 = calloc(nmax, 2), *o2 = calloc(nmax, 2);
            AT_RANK len = 0;
            int ret = CreateCheckSymmPaths(at, p1, c1, p2, c2, av, v1, v2, o1, o2, nlp1, nlp2, r1, same ? r1 : r2, crank, &len, &inv, mode);
            printf("P %d %d %d |", ret, len, inv); pa(v1, n); pa(v2, n); pa(o1, n); pa(o2, n); printf("\n");
        }
        int nb; scanf("%d", &nb);
        for (int k = 0; k < nb; k++) {
            int ib, ia, in1, in2, mode, symm2;
            scanf("%d %d %d %d %d %d", &ib, &ia, &in1, &in2, &mode, &symm2);
            AT_RANK perm[64];
            for (int i = 0; i < n; i++) { int x; scanf("%d", &x); perm[i] = x; }
            AT_RANK nn[4] = {0, 0, 0, 0};
            for (int i = 0; i < at[ia].valence && i < 4; i++) nn[i] = at[ia].neighbor[i];
            AT_RANK *s1[3], *s2[3];
            for (int i = 0; i < 2; i++) { s1[i] = calloc(nmax, 2); s2[i] = calloc(nmax, 2); }
            s1[2] = s2[2] = NULL;
            for (int i = 0; i < n; i++) s1[1][i] = perm[i];
            AT_RANK *tmp = calloc(nmax + 8, 2);
            NEIGH_LIST *nl = CreateNeighList(n, n, at, 0, NULL), *nl1 = CreateNeighList(n, n, at, 0, NULL), *nl2 = CreateNeighList(n, n, at, 0, NULL);
            long iter = 0;
            int ret = BreakNeighborsTie(&cg, at, n, n, ib, ia, nn, in1, in2, mode, s1, s2, tmp, nl, symm2 ? r2 : r1, crank, nl1, nl2, &iter);
            printf("B %d %ld |", ret, iter); pa(s1[0], n); pa(s1[1], n); pa(s2[0], n); pa(s2[1], n); printf("\n");
            pnl(nl1, n); pnl(nl2, n);
        }
    }
    return 0;
}
"""


def gen(rng):
    c = ms.gen(rng)
    atoms = c["atoms"]
    n = len(atoms)
    c["rs"] = [rng.choice([1, 1, 1, 2]) for _ in range(n)]
    r1 = [a["rank"] for a in atoms]
    def vis():
        return [rng.choice([0, 0, rng.randrange(1, n + 1)]) for _ in range(n)]
    c["checks"] = []
    for _ in range(rng.randrange(1, 8)):
        cur1 = rng.randrange(n)
        cur2 = rng.choice([cur1, rng.randrange(n)])
        nb1 = atoms[cur1]["nb"]
        nb2 = atoms[cur2]["nb"]
        n1 = rng.choice(nb1) if nb1 else cur1
        n2 = rng.choice(nb2) if nb2 else cur2
        vv1, vv2 = vis(), vis()
        if rng.random() < 0.5:
            vv1[n1] = n2 + 1
            vv2[n2] = n1 + 1
        o1, o2 = [rng.randrange(0, 3) for _ in range(n)], None
        o2 = list(o1) if rng.random() < 0.7 else [rng.randrange(0, 3) for _ in range(n)]
        c["checks"].append([cur1, cur2, n1, n2, rng.randrange(n), rng.randrange(n), rng.randrange(2), vv1, vv2, o1, o2])
    c["paths"] = []
    for _ in range(rng.randrange(1, 6)):
        cur1 = rng.randrange(n)
        cur2 = cur1 if rng.random() < 0.5 else rng.randrange(n)
        nb1 = atoms[cur1]["nb"]
        nb2 = atoms[cur2]["nb"]
        p1 = rng.choice(nb1) if nb1 and rng.random() < 0.8 else n + 3
        p2 = rng.choice(nb2) if nb2 and rng.random() < 0.8 else n + 3
        c["paths"].append([p1, cur1, p2, cur2, rng.randrange(n), rng.randrange(n), rng.choice([-1, -1, 0, 1]), rng.randrange(4), rng.randrange(2)])
    c["breaks"] = []
    cand = [i for i in range(n) if 2 <= len(atoms[i]["nb"]) <= 4]
    for _ in range(rng.randrange(0, 6) if cand else 0):
        ia = rng.choice(cand)
        v = len(atoms[ia]["nb"])
        in1, in2 = rng.sample(range(v), 2)
        mode = rng.choice([0, 0, 0, 1, 2, 3])
        if mode == 3 and rng.random() < 0.8:
            in1 = 0
            in2 = rng.randrange(1, v)
        ib = rng.choice([n, n, rng.randrange(n)])
        perm = list(range(n))
        rng.shuffle(perm)
        c["breaks"].append([ib, ia, in1, in2, mode, rng.randrange(2)] + perm)
    return c


def fmt(c):
    n = len(c["atoms"])
    lines = ms.fmt(c).split("\n")
    out = [lines[0]]
    for i in range(n):
        toks = lines[1 + 2 * i].split(" ")
        toks.insert(10, str(c["rs"][i]))
        out.append(" ".join(toks))
        out.append(lines[2 + 2 * i])
    out.append(str(len(c["checks"])))
    for ch in c["checks"]:
        out.append(" ".join(map(str, ch[:7])))
        for k in range(7, 11):
            out.append(" ".join(map(str, ch[k])))
    out.append(str(len(c["paths"])))
    out += [" ".join(map(str, p)) for p in c["paths"]]
    out.append(str(len(c["breaks"])))
    out += [" ".join(map(str, b)) for b in c["breaks"]]
    return "\n".join(out) + "\n"


def main():
    rng = random.Random(3113)
    cases = [gen(rng) for _ in range(3000)]
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True)
        assert r.returncode == 0, (r.returncode, r.stderr.decode()[:2000])
        out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_mapbreak.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
