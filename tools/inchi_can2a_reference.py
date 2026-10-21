"""Generates tests/data/inchi_can2a.json from the InChI C code: the
the partition, cell, node set and connection table helpers of ichican2.c
(compiled with the upstream source included so that its types are available), linked from the upstream
library.

    python tools/inchi_can2a_reference.py
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
#include "ichican2.c"
int CreateCompositeNormAtom() { return 0; }
static void pa(const char *t, AT_RANK *a, int n) { printf("%s", t); for (int i = 0; i < n; i++) printf(" %d", a[i]); printf("\n"); }
static void pset(const char *t, NodeSet *s, int l) { printf("%s", t); for (int i = 0; i < s->len_set; i++) printf(" %d", s->bitword[l][i]); printf("\n"); }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        SetBitCreate(&cg);
        int n;
        scanf("%d", &n);
        int nmax = n + 2;
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        Partition p0, pB;
        PartitionCreate(&p0, nmax); PartitionCreate(&pB, nmax);
        for (int i = 0; i < n; i++) {
            int v, p, sap, r, cr, nh, an, ac, rr2, aa2, rs;
            scanf("%d %d %d %d %d %d %d %d %d %d %d", &v, &p, &sap, &r, &cr, &nh, &an, &ac, &rr2, &aa2, &rs);
            at[i].valence = v;
            p0.Rank[i] = r; p0.AtNumber[i] = an; pB.Rank[i] = rr2; pB.AtNumber[i] = aa2;
            for (int j = 0; j < v; j++) { int q; scanf("%d", &q); at[i].neighbor[j] = q; }
            for (int j = 0; j < 3; j++) { int a, b, z, pp; scanf("%d %d %d %d", &a, &b, &z, &pp); }
        }
        int marks[64], useaux, aux[64];
        for (int i = 0; i < n; i++) scanf("%d", &marks[i]);
        scanf("%d", &useaux);
        for (int i = 0; i < n; i++) scanf("%d", &aux[i]);
        int nv1, nv2; Node v1[64], v2[64];
        scanf("%d", &nv1); for (int i = 0; i < nv1; i++) { int x; scanf("%d", &x); v1[i] = x; }
        scanf("%d", &nv2); for (int i = 0; i < nv2; i++) { int x; scanf("%d", &x); v2[i] = x; }
        UnorderedPartition ua, ub;
        UnorderedPartitionCreate(&ua, nmax); UnorderedPartitionCreate(&ub, nmax);
        for (int i = 0; i < n; i++) { int x; scanf("%d", &x); ua.equ2[i] = x; }
        for (int i = 0; i < n; i++) { int x; scanf("%d", &x); ub.equ2[i] = x; }
        int njoin; int jn[16][2];
        scanf("%d", &njoin); for (int i = 0; i < njoin; i++) scanf("%d %d", &jn[i][0], &jn[i][1]);
        int ncv, cvs[8], digraph;
        scanf("%d %d", &ncv, &digraph); for (int i = 0; i < ncv; i++) scanf("%d", &cvs[i]);
        int ctmax, ctat, hasnumh, hasfix, hasiso, hasexc, clearK, ctfill;
        scanf("%d %d %d %d %d %d %d %d", &ctmax, &ctat, &hasnumh, &hasfix, &hasiso, &hasexc, &clearK, &ctfill);
        int nl[8], nnl, maxrj;
        scanf("%d %d", &nnl, &maxrj); for (int i = 0; i < nnl; i++) scanf("%d", &nl[i]);
        printf("N %d\n", c);

        printf("D %d %d\n", PartitionIsDiscrete(&p0, n), PartitionIsDiscrete(&pB, n));
        Partition *ps[2] = {&p0, &pB};
        for (int q = 0; q < 2; q++) {
            Cell W[70];
            memset(W, 0, sizeof(W));
            for (int k = 1; k < 66; k++) {
                int sz = PartitionGetFirstCell(ps[q], W, k, n);
                printf("C %d %d %d %d\n", k, sz, W[k - 1].first, W[k - 1].next);
                if (!sz) { CellMakeEmpty(W, k); printf("E %d %d %d\n", W[k - 1].first, W[k - 1].next, W[k - 1].prev); break; }
            }
        }
        {
            NodeSet Mcr, Fix;
            NodeSetCreate(&cg, &Mcr, n, 2); NodeSetCreate(&cg, &Fix, n, 2);
            PartitionGetMcrAndFixSet(&cg, &p0, &Mcr, &Fix, n, 1);
            PartitionGetMcrAndFixSet(&cg, &pB, &Mcr, &Fix, n, 2);
            pset("M1", &Mcr, 0); pset("F1", &Fix, 0); pset("M2", &Mcr, 1); pset("F2", &Fix, 1);
            printf("AN %d %d %d %d\n", AllNodesAreInSet(&Mcr, 1, &Fix, 1), AllNodesAreInSet(&Fix, 1, &Mcr, 1), AllNodesAreInSet(&Mcr, 2, &Mcr, 1), AllNodesAreInSet(&Fix, 2, &Mcr, 2));
            Partition pw;
            PartitionCreate(&pw, nmax);
            PartitionCopy(&pw, &p0, n);
            for (int i = 0; i < n; i++) if (marks[i]) pw.Rank[i] |= rank_mark_bit;
            CANON_DATA cd; memset(&cd, 0, sizeof(cd));
            AT_RANK auxr[64]; for (int i = 0; i < n; i++) auxr[i] = aux[i];
            cd.nAuxRank = useaux ? auxr : NULL;
            Cell W[70]; memset(W, 0, sizeof(W));
            for (int k = 1; k < 66; k++) {
                int sz = PartitionGetFirstCell(&p0, W, k, n);
                if (!sz) break;
                int nn = CellGetNumberOfNodes(&pw, &W[k - 1]);
                int ki = CellIntersectWithSet(&cg, &pw, &W[k - 1], &Mcr, 1);
                int nn2 = CellGetNumberOfNodes(&pw, &W[k - 1]);
                printf("W %d %d %d %d |", k, nn, ki, nn2);
                Node v = 0;
                for (int it = 0; it < 70; it++) { v = CellGetMinNode(&pw, &W[k - 1], v, &cd); printf(" %d", v); if (v == INCHI_CANON_INFINITY) break; }
                printf("\n");
            }
            pa("R", pw.Rank, n);
            Cell Wemp; CellMakeEmpty(&Wemp, 1);
            printf("WE %d %d\n", CellGetMinNode(&pw, &Wemp, 0, &cd), CellGetMinNode(&pw, &Wemp, 3, NULL));
        }
        {
            NodeSet S; NodeSetCreate(&cg, &S, n, 2);
            NodeSetFromVertices(&cg, &S, 1, v1, nv1);
            NodeSetFromVertices(&cg, &S, 2, v2, nv2);
            pset("S1", &S, 0); pset("S2", &S, 1);
            printf("SA %d %d %d\n", AllNodesAreInSet(&S, 1, &S, 2), AllNodesAreInSet(&S, 2, &S, 1), AllNodesAreInSet(&S, 1, &S, 1));
        }
        {
            Transposition g; TranspositionCreate(&g, nmax);
            PartitionGetTransposition(&p0, &pB, n, &g);
            pa("T", g.nAtNumb, n);
        }
        {
            UnorderedPartition a2, b2;
            UnorderedPartitionCreate(&a2, nmax); UnorderedPartitionCreate(&b2, nmax);
            memcpy(a2.equ2, ua.equ2, nmax * 2); memcpy(b2.equ2, ub.equ2, nmax * 2);
            printf("G"); for (int i = 0; i < n; i++) printf(" %d", nGetMcr2(a2.equ2, i)); printf("\n");
            pa("GA", a2.equ2, n);
            printf("GN"); for (int i = 1; i <= n; i++) printf(" %d", GetUnorderedPartitionMcrNode(&b2, i)); printf("\n");
            for (int i = 0; i < njoin; i++) { printf("J %d", nJoin2Mcrs2(b2.equ2, jn[i][0], jn[i][1])); printf("\n"); }
            pa("JB", b2.equ2, n);
            printf("UJ %d\n", UnorderedPartitionJoin(&ua, &b2, n));
            pa("UB", b2.equ2, n);
            UnorderedPartitionMakeDiscrete(&a2, n); pa("MD", a2.equ2, n);
        }
        printf("L %d %d\n", PartitionSatisfiesLemma_2_25(&p0, n), PartitionSatisfiesLemma_2_25(&pB, n));
        for (int t = 0; t < ncv; t++) {
            NEIGH_LIST *g = CreateNeighList(n, n, at, 0, NULL);
            Partition pp[3]; memset(pp, 0, sizeof(pp));
            PartitionCreate(&pp[0], nmax);
            PartitionCopy(&pp[0], &p0, n);
            int ret = PartitionColorVertex(&cg, g, pp, cvs[t], n, n, nmax, digraph, 3);
            printf("PV %d |", ret); for (int i = 0; i < n; i++) printf(" %d", pp[1].Rank[i]); printf(" |"); for (int i = 0; i < n; i++) printf(" %d", pp[1].AtNumber[i]); printf("\n");
        }
        {
            CANON_DATA cd; memset(&cd, 0, sizeof(cd));
            cd.nMaxLenLinearCT = ctmax; cd.nLenCTAtOnly = ctat;
            NUM_H nh[4]; AT_ISO_SORT_KEY isk[4]; S_CHAR exc[4];
            if (hasnumh) { cd.NumH = nh; cd.maxlenNumH = 5; }
            if (hasfix) { cd.NumHfixed = nh; cd.maxlenNumHfixed = 6; }
            if (hasiso) { cd.iso_sort_key = isk; cd.maxlen_iso_sort_key = 7; }
            if (hasexc) { cd.iso_exchg_atnos = exc; cd.maxlen_iso_exchg_atnos = 8; }
            ConTable ct;
            int rc = CTableCreate(&ct, n, &cd);
            printf("CT %d %d %d %d %d %d %d %d %d %d %d %d\n", rc, ct.maxlenCt, ct.maxPos, ct.maxVert, ct.maxlenNumH, ct.maxlen_iso_sort_key, ct.maxlen_iso_exchg_atnos, ct.lenCt, ct.lenPos, ct.nLenCTAtOnly, ct.NumH != NULL, ct.NumHfixed != NULL);
            for (int i = 0; i < ctmax && i < ctfill; i++) ct.Ctbl[i] = 100 + i;
            for (int i = 0; i <= n; i++) ct.nextCtblPos[i] = i < ctfill ? i : ctfill;
            ct.lenCt = ctfill < ctmax ? ctfill : ctmax; ct.lenPos = n;
            CtPartClear(&ct, clearK);
            printf("CC %d %d |", ct.lenCt, ct.lenPos); for (int i = 0; i < ctmax; i++) printf(" %d", ct.Ctbl[i]); printf("\n");
            CTableFree(&ct);
        }
        {
            AT_RANK rk[64]; for (int i = 0; i < n; i++) rk[i] = p0.Rank[i] | (marks[i] ? rank_mark_bit : 0);
            NEIGH_LIST l = calloc(10, 2);
            l[0] = nnl; for (int i = 0; i < nnl; i++) l[i + 1] = nl[i];
            insertions_sort_NeighList_AT_NUMBERS2(l, rk, maxrj);
            printf("IS"); for (int i = 0; i <= nnl; i++) printf(" %d", l[i]); printf("\n");
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
    for a in atoms:
        a["sb"] = [[0, 0, 0, 0] for _ in range(3)]
        if a["par"] and rng.random() < 0.8:
            a["sap"] = rng.choice([6, 6, 6, 3, 4])
    edges = [(i, j) for i in range(n) for j in atoms[i]["nb"] if i < j and len(atoms[i]["nb"]) == 3 and len(atoms[j]["nb"]) == 3]
    rng.shuffle(edges)
    for (i, j) in edges[: rng.randrange(0, 3)]:
        if atoms[i]["sb"][0][0] or atoms[j]["sb"][0][0]:
            continue
        pp = rng.choice([6, 6, 6, 3, 4, 1, 2])
        atoms[i]["sb"][0] = [j + 1, atoms[i]["nb"].index(j), 1, pp]
        atoms[j]["sb"][0] = [i + 1, atoms[j]["nb"].index(i), 1, pp]
        for a in (atoms[i], atoms[j]):
            a["par"] = a["par"] or rng.choice([1, 2])
    carb, dble = [], []
    for i, a in enumerate(atoms):
        if a["par"] and not a["sb"][0][0]:
            carb.append([a["crank"], a["sap"] & 7 or 1])
    for i, a in enumerate(atoms):
        if a["sb"][0][0]:
            j = a["sb"][0][0] - 1
            if i < j:
                dble.append([max(a["crank"], atoms[j]["crank"]), min(a["crank"], atoms[j]["crank"]), a["sb"][0][3] & 7 or 1])
    if rng.random() < 0.1 and carb:
        carb.pop()
    c["carb"], c["dble"] = carb, dble
    return c


def fmt(c):
    n = len(c["atoms"])
    lines = ms.fmt(c).split("\n")
    out = [lines[0]]
    for i in range(n):
        toks = lines[1 + 2 * i].split(" ")
        toks.insert(10, str(c["rs"][i]))
        out.append(" ".join(toks))
        out.append(" ".join(" ".join(map(str, s)) for s in c["atoms"][i]["sb"]))
    out.append(str(len(c["carb"])))
    out += [" ".join(map(str, x)) for x in c["carb"]]
    out.append(str(len(c["dble"])))
    out += [" ".join(map(str, x)) for x in c["dble"]]
    return "\n".join(out) + "\n"




def gen2(rng):
    c = gen(rng)
    n = len(c["atoms"])
    c["marks"] = [rng.choice([0, 0, 1]) for _ in range(n)]
    c["useaux"] = rng.randrange(2)
    c["aux"] = [rng.randrange(1, 4) for _ in range(n)]
    c["v1"] = [rng.randrange(1, n + 1) for _ in range(rng.randrange(0, 5))]
    c["v2"] = [rng.randrange(1, n + 1) for _ in range(rng.randrange(0, 5))]
    def equ():
        e = [rng.randrange(n) for _ in range(n)]
        return [min(e[i], i) if rng.random() < 0.7 else i for i in range(n)]
    c["ua"], c["ub"] = equ(), equ()
    c["joins"] = [[rng.randrange(n), rng.randrange(n)] for _ in range(rng.randrange(0, 5))]
    c["cv"] = [rng.randrange(1, n + 1) for _ in range(rng.randrange(0, 4))]
    c["digraph"] = rng.randrange(2)
    ctmax = rng.randrange(3, 12)
    c["ct"] = [ctmax, rng.randrange(0, 4), rng.randrange(2), rng.randrange(2), rng.randrange(2), rng.randrange(2), rng.randrange(0, 5), rng.randrange(0, ctmax + 2)]
    nnl = rng.randrange(0, 6)
    c["nl"] = [rng.randrange(n) for _ in range(nnl)]
    c["maxrj"] = rng.choice([0x7FFF, 3, 5, 100])
    return c


def fmt2(c):
    n = len(c["atoms"])
    lines = fmt(c).rstrip("\n").split("\n")
    out = lines[: 1 + 2 * n]
    out.append(" ".join(map(str, c["marks"])))
    out.append(str(c["useaux"]))
    out.append(" ".join(map(str, c["aux"])))
    out.append(" ".join(map(str, [len(c["v1"])] + c["v1"])))
    out.append(" ".join(map(str, [len(c["v2"])] + c["v2"])))
    out.append(" ".join(map(str, c["ua"])))
    out.append(" ".join(map(str, c["ub"])))
    out.append(" ".join(map(str, [len(c["joins"])] + [x for j in c["joins"] for x in j])))
    out.append(" ".join(map(str, [len(c["cv"]), c["digraph"]] + c["cv"])))
    out.append(" ".join(map(str, c["ct"])))
    out.append(" ".join(map(str, [len(c["nl"]), c["maxrj"]] + c["nl"])))
    return "\n".join(out) + "\n"


def main():
    rng = random.Random(7117)
    cases = [gen2(rng) for _ in range(3000)]
    stdin = f"{len(cases)}\n" + "".join(fmt2(c) for c in cases)
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
    path = ROOT / "tests/data/inchi_can2a.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
