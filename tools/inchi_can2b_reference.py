"""Generates tests/data/inchi_can2b.json from the InChI C code: the
the connection table part helpers of ichican2.c (`CtPartFill`, `CtPartCompare`,
`CtFullCompare`, `CtPartCopy`, the layer helpers, ...)
(compiled with the upstream source included so that its types are available), linked from the upstream
library.

    python tools/inchi_can2b_reference.py
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
static void build(ConTable *Ct, CANON_DATA *cd, sp_ATOM *at, int n, int *crank, int *ms, int nms) {
    CTableCreate(Ct, n, cd);
    NEIGH_LIST *G = CreateNeighList(n, n, at, 0, NULL);
    Partition p; PartitionCreate(&p, n + 2);
    for (int t = 0; t < nms; t++) {
        int m = ms[t], pos = 0;
        for (int r = 1; r <= m; r++) for (int x = 0; x < n; x++) if (crank[x] == r) p.AtNumber[pos++] = x;
        for (int x = 0; x < n; x++) if (crank[x] > m) p.AtNumber[pos++] = x;
        for (int x = 0; x < n; x++) p.Rank[x] = crank[x] <= m ? crank[x] : n;
        CtPartFill(G, cd, &p, Ct, t + 1, n, n, n + 2);
    }
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        SetBitCreate(&cg);
        int n;
        scanf("%d", &n);
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        int crank[64], crankB[64];
        for (int i = 0; i < n; i++) {
            int v, p, sap, r, cr, nh, an, ac, rr2, aa2, rs;
            scanf("%d %d %d %d %d %d %d %d %d %d %d", &v, &p, &sap, &r, &cr, &nh, &an, &ac, &rr2, &aa2, &rs);
            at[i].valence = v; crank[i] = cr;
            for (int j = 0; j < v; j++) { int q; scanf("%d", &q); at[i].neighbor[j] = q; }
            for (int j = 0; j < 3; j++) { int a, b, z, pp; scanf("%d %d %d %d", &a, &b, &z, &pp); }
        }
        for (int i = 0; i < n; i++) scanf("%d", &crankB[i]);
        int na, nb, msA[64], msB[64];
        scanf("%d", &na); for (int i = 0; i < na; i++) scanf("%d", &msA[i]);
        scanf("%d", &nb); for (int i = 0; i < nb; i++) scanf("%d", &msB[i]);
        int flags[8];
        for (int i = 0; i < 8; i++) scanf("%d", &flags[i]);
        int bOnly = flags[0], bSplit = flags[1], atOnly = flags[2];
        NUM_H nh[80], nhf[80]; AT_ISO_SORT_KEY isk[80]; S_CHAR exc[80];
        for (int i = 0; i < n + 2; i++) { int a, b, c2, d; scanf("%d %d %d %d", &a, &b, &c2, &d); nh[i] = a; nhf[i] = b; isk[i] = c2; exc[i] = d; }
        int perm[64];
        for (int i = 0; i < n; i++) scanf("%d", &perm[i]);
        CANON_DATA cd; memset(&cd, 0, sizeof(cd));
        cd.nMaxLenLinearCT = 8 * n + 8; cd.nLenLinearCT = flags[3]; cd.nLenCTAtOnly = atOnly;
        if (flags[4] & 1) { cd.NumH = nh; cd.maxlenNumH = n + 2; }
        if (flags[4] & 2) { cd.NumHfixed = nhf; cd.maxlenNumHfixed = n + 2; }
        if (flags[4] & 4) { cd.iso_sort_key = isk; cd.maxlen_iso_sort_key = n + 2; }
        if (flags[4] & 8) { cd.iso_exchg_atnos = exc; cd.maxlen_iso_exchg_atnos = n + 2; }
        printf("N %d\n", c);
        ConTable A, B;
        build(&A, &cd, at, n, crank, msA, na);
        build(&B, &cd, at, n, crankB, msB, nb);
        dump("A", &A, n); dump("B", &B, n);
        {
            int kmax = A.lenPos < B.lenPos ? A.lenPos : B.lenPos;
            S_CHAR cmp[128]; kLeast kl[MAX_LAYERS];
            memset(cmp, 0, sizeof(cmp)); memset(kl, 0, sizeof(kl));
            for (int k = 1; k <= kmax; k++) {
                int d1 = CtPartCompare(&A, &B, cmp, kl, k, bOnly, bSplit);
                printf("P1 %d %d %d |", k, d1, cmp[k - 1]);
                for (int l = 0; l < 7; l++) printf(" %d,%d", kl[l].k, kl[l].i);
                printf("\n");
                int d2 = CtPartCompare(&A, &B, NULL, NULL, k, bOnly, bSplit);
                int d3 = CtPartCompare(&B, &A, NULL, NULL, k, !bOnly, bSplit);
                printf("P2 %d %d\n", d2, d3);
            }
            printf("FC %d %d %d %d\n", CtFullCompare(&A, &B, bOnly, bSplit), CtFullCompare(&A, &B, !bOnly, !bSplit), CtFullCompare(&B, &A, 0, 0), CtFullCompare(&A, &A, 0, 0));
            printf("FL %d\n", CtFullCompareLayers(kl));
            int L, I, K;
            int g1 = CtCompareLayersGetFirstDiff(kl, 2, &L, &I, &K);
            printf("GD %d %d %d %d\n", g1, L, I, K);
            g1 = CtCompareLayersGetFirstDiff(NULL, 2, &L, &I, &K);
            printf("GN %d\n", g1);
            printf("PL %d %d %d\n", CtPartCompareLayers(kl, 5, 0), CtPartCompareLayers(kl, 0, 3), CtPartCompareLayers(NULL, 1, 1));
            UpdateCompareLayers(kl, 2);
            printf("UL"); for (int l = 0; l < 7; l++) printf(" %d,%d", kl[l].k, kl[l].i); printf("\n");
            UpdateCompareLayers(NULL, 2);
        }
        {
            ConTable C, D;
            CTableCreate(&C, n, &cd); CTableCreate(&D, n, &cd);
            for (int k = 1; k <= A.lenPos; k++) CtPartCopy(&C, &A, k);
            dump("C", &C, n);
            CtFullCopy(&D, &B);
            dump("D", &D, n);
            S_CHAR cmp[128]; memset(cmp, 1, sizeof(cmp));
            CtPartINCHI_CANON_INFINITY(&C, cmp, A.lenPos);
            printf("INF"); for (int i = 0; i <= C.lenCt; i++) printf(" %d", C.Ctbl[i]); printf(" |"); for (int i = 0; i < 8; i++) printf(" %d", cmp[i]); printf("\n");
            CtPartINCHI_CANON_INFINITY(&D, NULL, 1);
            printf("INF1 %d\n", D.Ctbl[0]);
            CtPartClear(&C, 1);
            printf("CLR %d\n", GetOneAdditionalLayer(&cd, &C));
            printf("GOL %d %d\n", GetOneAdditionalLayer(NULL, &A), GetOneAdditionalLayer(&cd, NULL));
        }
        {
            Transposition g; TranspositionCreate(&g, n + 2);
            for (int i = 0; i < n; i++) g.nAtNumb[i] = perm[i];
            NodeSet M, F; NodeSetCreate(&cg, &M, n, 2); NodeSetCreate(&cg, &F, n, 2);
            UnorderedPartition up; UnorderedPartitionCreate(&up, n + 2);
            TranspositionGetMcrAndFixSetAndUnorderedPartition(&cg, &g, &M, &F, n, 2, &up);
            printf("TG"); for (int i = 0; i < n; i++) printf(" %d", g.nAtNumb[i]);
            printf(" | M"); for (int i = 0; i < M.len_set; i++) printf(" %d", M.bitword[1][i]);
            printf(" | F"); for (int i = 0; i < F.len_set; i++) printf(" %d", F.bitword[1][i]);
            printf(" | U"); for (int i = 0; i < n; i++) printf(" %d", up.equ2[i]); printf("\n");
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
    crankA = [a["crank"] for a in c["atoms"]]
    perm = list(range(1, n + 1))
    rng.shuffle(perm)
    c["crankB"] = perm if rng.random() < 0.6 else list(crankA)
    if c["crankB"] == crankA and rng.random() < 0.5:
        # swap two ranks
        i, j = rng.sample(range(n), 2)
        c["crankB"][i], c["crankB"][j] = c["crankB"][j], c["crankB"][i]
    def levels():
        k = rng.randrange(1, min(n, 4) + 1)
        ms = sorted(rng.sample(range(1, n + 1), k))
        return ms
    c["msA"], c["msB"] = levels(), levels()
    if rng.random() < 0.5:
        c["msB"] = list(c["msA"])
    c["flags"] = [rng.randrange(2), rng.randrange(2), rng.randrange(0, 40), rng.randrange(0, 60), rng.randrange(16), 0, 0, 0]
    base = [[rng.randrange(0, 3), rng.randrange(0, 3), rng.randrange(0, 3), rng.randrange(-1, 2)] for _ in range(n + 2)]
    c["vals"] = base
    c["perm"] = list(range(n))
    rng.shuffle(c["perm"])
    return c


def fmt2(c):
    n = len(c["atoms"])
    lines = fmt(c).rstrip("\n").split("\n")
    out = lines[: 1 + 2 * n]
    out.append(" ".join(map(str, c["crankB"])))
    out.append(" ".join(map(str, [len(c["msA"])] + c["msA"])))
    out.append(" ".join(map(str, [len(c["msB"])] + c["msB"])))
    out.append(" ".join(map(str, c["flags"])))
    for v in c["vals"]:
        out.append(" ".join(map(str, v)))
    out.append(" ".join(map(str, c["perm"])))
    return "\n".join(out) + "\n"


def main():
    rng = random.Random(8118)
    cases = [gen2(rng) for _ in range(4000)]
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
    path = ROOT / "tests/data/inchi_can2b.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
