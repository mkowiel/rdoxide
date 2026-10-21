"""Generates tests/data/inchi_can2f.json from the InChI C code:
`GetBaseCanonRanking` of ichican2.c (the base canonical numbering of a
tautomeric / non-tautomeric structure with its fixed-H and isotopic variants),
compiled with the upstream source included.

    python tools/inchi_can2f_reference.py
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
d2 = load("inchi_can2d_reference")
SRC = ma.SRC

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ichican2.c"
int CreateCompositeNormAtom() { return 0; }
static void pu(const char *t, AT_RANK *a, int n) { printf("%s", t); if (!a) { printf(" -\n"); return; } for (int i = 0; i < n; i++) printf(" %d", a[i]); printf("\n"); }
static void ph(const char *t, NUM_H *a, int n) { printf("%s", t); if (!a) { printf(" -\n"); return; } for (int i = 0; i < n; i++) printf(" %d", a[i]); printf("\n"); }
static void pk(const char *t, AT_ISO_SORT_KEY *a, int n) { printf("%s", t); if (!a) { printf(" -\n"); return; } for (int i = 0; i < n; i++) printf(" %ld", (long) a[i]); printf("\n"); }
static void pc(const char *t, S_CHAR *a, int n) { printf("%s", t); if (!a) { printf(" -\n"); return; } for (int i = 0; i < n; i++) printf(" %d", a[i]); printf("\n"); }
static void read_atoms(sp_ATOM *at, int n, int *edges, int m, int withTaut) {
    for (int i = 0; i < n; i++) {
        int el, nh, aw, i1, i2, i3, cf, ep;
        scanf("%d %d %d %d %d %d %d %d", &el, &nh, &aw, &i1, &i2, &i3, &cf, &ep);
        const char *names[] = {"C", "N", "O", "S", "Cl", "H"};
        strcpy(at[i].elname, names[el]);
        at[i].num_H = nh; at[i].iso_atw_diff = aw;
        at[i].num_iso_H[0] = i1; at[i].num_iso_H[1] = i2; at[i].num_iso_H[2] = i3;
        at[i].cFlags = cf; at[i].endpoint = ep;
    }
    for (int i = 0; i < m; i++) {
        int a = edges[2 * i], b = edges[2 * i + 1];
        at[a].neighbor[at[a].valence++] = b;
        at[b].neighbor[at[b].valence++] = a;
    }
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        int n, m, ntg, nE, haveYes, haveNon, haveTgi, bfix, large;
        scanf("%d %d %d %d %d %d %d %d %d", &n, &m, &ntg, &nE, &haveYes, &haveNon, &haveTgi, &bfix, &large);
        int edges[200];
        for (int i = 0; i < 2 * m; i++) scanf("%d", &edges[i]);
        ATOM_SIZES s[2]; memset(s, 0, sizeof(s));
        for (int k = 0; k < 2; k++) {
            scanf("%d %d %d %d %d %d %d", &s[k].nLenCT, &s[k].nLenCTAtOnly, &s[k].num_isotopic_atoms, &s[k].bIgnoreIsotopic, &s[k].nLenLinearCTTautomer, &s[k].bHasIsotopicTautGroups, &s[k].nLenIsotopicEndpoints);
        }
        sp_ATOM *atY = calloc(n + 1, sizeof(sp_ATOM)), *atN = calloc(n + 1, sizeof(sp_ATOM));
        read_atoms(atY, n, edges, m, 1);
        read_atoms(atN, n, edges, m, 0);
        T_GROUP *tg = calloc(ntg + 2, sizeof(T_GROUP));
        AT_NUMB *epl = calloc(nE + 2, sizeof(AT_NUMB));
        for (int k = 0; k < ntg; k++) {
            int ne, fp; long w;
            scanf("%d %d %ld", &ne, &fp, &w);
            tg[k].nNumEndpoints = ne; tg[k].nFirstEndpointAtNoPos = fp; tg[k].nGroupNumber = k + 1;
            for (int j = 0; j < T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC; j++) { int x; scanf("%d", &x); tg[k].num[j] = x; }
            long iw; scanf("%ld", &iw); tg[k].iWeight = iw;
        }
        for (int i = 0; i < nE; i++) { int x; scanf("%d", &x); epl[i] = x; }
        int nIso, nNumIsoEnd, flagsDone, normFlags, ignIso, tgNull;
        scanf("%d %d %d %d %d %d", &nIso, &nNumIsoEnd, &flagsDone, &normFlags, &ignIso, &tgNull);
        AT_NUMB *isoList = calloc(n + 3, sizeof(AT_NUMB));
        isoList[0] = nIso;
        for (int i = 0; i < nIso; i++) { int x; scanf("%d", &x); isoList[i + 1] = x; }
        T_GROUP_INFO tgi; memset(&tgi, 0, sizeof(tgi));
        tgi.t_group = tgNull ? NULL : tg; tgi.nEndpointAtomNumber = epl; tgi.num_t_groups = ntg; tgi.max_num_t_groups = ntg;
        tgi.nNumEndpoints = nE;
        tgi.nIsotopicEndpointAtomNumber = isoList; tgi.nNumIsotopicEndpoints = nNumIsoEnd;
        tgi.bTautFlagsDone = flagsDone; tgi.tni.bNormalizationFlags = normFlags; tgi.bIgnoreIsotopic = ignIso;
        BCN bcn; memset(&bcn, 0, sizeof(bcn));
        sp_ATOM *at[2];
        at[TAUT_NON] = haveNon ? atN : NULL;
        at[TAUT_YES] = haveYes ? atY : NULL;
        int ntot = n + ntg;
        printf("N %d\n", c);
        int ret = GetBaseCanonRanking(NULL, n, ntot, at, haveTgi ? &tgi : NULL, s, &bcn, NULL, &cg, bfix, large);
        printf("R %d %d %d %d %d\n", ret, bcn.num_max, bcn.num_at_tg, bcn.num_atoms, bcn.nMaxLenRankStack);
        printf("TI %d\n", tgi.bIgnoreIsotopic);
        for (int j = 0; j < bcn.nMaxLenRankStack / 2; j++) {
            AT_RANK *a = bcn.pRankStack ? bcn.pRankStack[2 * j] : NULL, *b = bcn.pRankStack ? bcn.pRankStack[2 * j + 1] : NULL;
            char t[16]; sprintf(t, "RS%d", j); pu(t, a, bcn.num_max);
            sprintf(t, "RT%d", j); pu(t, b, bcn.num_max);
        }
        if (ret == 0) {
            for (int k = 0; k < 2; k++) {
                FTCN *f = &bcn.ftcn[k];
                printf("F%d %d %d %d %d %d %d %d %d %d\n", k, f->num_at_tg, f->num_atoms, f->nCanonFlags, f->nLenLinearCtAtOnly, f->nLenLinearCt, f->nMaxLenLinearCt, f->nLenNumH, f->nLenNumHFixH, f->len_iso_sort_keys);
                printf("NL");
                if (f->NeighList) for (int i = 0; i < f->num_at_tg; i++) { printf(" |"); for (int j = 1; j <= f->NeighList[i][0]; j++) printf(" %d", f->NeighList[i][j]); }
                else printf(" -");
                printf("\n");
                pu("LC", f->LinearCt, f->LinearCt ? f->nLenLinearCt + 1 : 0);
                pu("PR", f->PartitionCt.Rank, bcn.num_max); pu("PA", f->PartitionCt.AtNumber, bcn.num_max);
                pu("SY", f->nSymmRankCt, bcn.num_max);
                ph("NHO", f->nNumHOrig, f->nLenNumH); ph("NH", f->nNumH, f->nLenNumH);
                ph("NFO", f->nNumHOrigFixH, f->nLenNumH); ph("NF", f->nNumHFixH, f->nLenNumH);
                pu("IR", f->PartitionCtIso.Rank, bcn.num_max); pu("IA", f->PartitionCtIso.AtNumber, bcn.num_max);
                pu("IS", f->nSymmRankCtIso, bcn.num_max);
                pk("K", f->iso_sort_keys, f->len_iso_sort_keys); pk("KO", f->iso_sort_keysOrig, f->len_iso_sort_keys);
                pc("X", f->iso_exchg_atnos, f->len_iso_sort_keys); pc("XO", f->iso_exchg_atnosOrig, f->len_iso_sort_keys);
            }
        }
    }
    return 0;
}
"""

ELS = ["C", "N", "O", "S", "Cl"]
NORM = 1


def gen(rng):
    name = rng.choice(d2.NAMES)
    n, e = d2.edges_of(name, rng)
    n = max(n, 1)
    if n > 14:
        n, e = 6, [(i, (i + 1) % 6) for i in range(6)]
    perm = list(range(n))
    rng.shuffle(perm)
    e = [(perm[a], perm[b]) for a, b in e]
    m = len(e)
    mode = rng.choices(["taut", "non", "yes"], [6, 2, 2])[0]
    el = [rng.choice([0, 0, 1, 2, 3, 4]) if rng.random() < 0.7 else 0 for _ in range(n)]
    if rng.random() < 0.3:
        el = [0] * n
    hetero = [i for i in range(n) if el[i] in (1, 2, 3)]
    # t-groups
    groups = []
    if mode == "taut":
        avail = hetero[:]
        rng.shuffle(avail)
        for _ in range(rng.choice([0, 1, 1, 2, 3])):
            k = rng.randrange(1, 4)
            if len(avail) < k:
                break
            groups.append(sorted(avail[:k]))
            avail = avail[k:]
    ntg = len(groups)
    ep = [0] * n
    for g, mem in enumerate(groups):
        for a in mem:
            ep[a] = g + 1
    epl = [a + 0 for mem in groups for a in mem]
    nE = len(epl)
    first = []
    pos = 0
    for mem in groups:
        first.append(pos)
        pos += len(mem)

    def atom(i, endpoint, base=None):
        if base is not None and rng.random() < 0.8:
            a = list(base)
            a[7] = endpoint
            if rng.random() < 0.3:
                a[1] = rng.randrange(0, 3)
            return a
        return [el[i], rng.randrange(0, 4) if el[i] != 4 else 0, rng.choice([0, 0, 0, 0, 1, -1, 2]),
                rng.choice([0, 0, 0, 1, 2]), rng.choice([0, 0, 0, 1]), rng.choice([0, 0, 0, 1]),
                1 if rng.random() < 0.3 else 0, endpoint]
    atY = [atom(i, ep[i]) for i in range(n)]
    atN = [atom(i, 0, atY[i]) for i in range(n)]
    tgs = []
    for g, mem in enumerate(groups):
        num = [rng.randrange(0, 4), rng.randrange(0, 2)] + [rng.choice([0, 0, 1, 2]) for _ in range(3)]
        tgs.append([len(mem), first[g]] + num)
    isoatoms = [i for i in range(n) if ep[i] == 0 and atY[i][6] & 1]
    nIso = len(isoatoms)
    if rng.random() < 0.15:
        nIso = max(0, nIso + rng.choice([-1, 1]))
        isoatoms = (isoatoms + [0])[:nIso]
    def isoany(a):
        return a[2] != 0 or a[3] or a[4] or a[5]
    sY = [n + ntg + m + nE, n + m, rng.choice([0, 0, sum(isoany(a) for a in atY), 1, 2]), rng.choice([0, 0, 1]),
          rng.choice([0, 1, 5, n + m]) if mode == "taut" else rng.choice([0, 3]),
          rng.choice([0, 0, 1]), rng.choice([0, 1, 2, 3])]
    sN = [n + m, n + m, rng.choice([0, 0, sum(isoany(a) for a in atN), 1, 2]), rng.choice([0, 0, 1]), 0, 0, 0]
    if mode == "taut":
        sY[4] = max(sY[4], 1) if rng.random() < 0.9 else sY[4]
    elif mode == "yes":
        sY = [n + m, n + m] + sY[2:4] + [0] + sY[5:]
    c = {
        "n": n, "m": m, "ntg": ntg, "nE": nE,
        "haveYes": 0 if mode == "non" else 1,
        "haveNon": 1 if mode == "non" else (1 if mode == "taut" and rng.random() < 0.7 else (1 if mode == "yes" and rng.random() < 0.2 else 0)),
        "haveTgi": 1 if (mode == "taut" or rng.random() < 0.5) else 0,
        "bfix": rng.randrange(2), "large": 0,
        "edges": e, "sizes": sN + sY, "atY": atY, "atN": atN, "tgs": tgs, "epl": epl,
        "iso": [nIso, rng.choice([0, 1, nIso + 1, nIso + 2]), rng.choice([0, 0x4000, 0x8000, 0xC000]),
                rng.choice([0, NORM]) if mode == "taut" and ntg == 0 else (NORM if mode == "taut" and rng.random() < 0.1 else 0),
                rng.choice([0, 0, 1]), 1 if (mode == "taut" and ntg == 0 and False) else 0],
        "isolist": isoatoms,
    }
    if mode == "taut" and ntg == 0:
        c["iso"][3] = NORM
    return c


def fmt(c):
    o = [" ".join(map(str, [c["n"], c["m"], c["ntg"], c["nE"], c["haveYes"], c["haveNon"], c["haveTgi"], c["bfix"], c["large"]]))]
    o.append(" ".join(f"{a} {b}" for a, b in c["edges"]))
    o.append(" ".join(map(str, c["sizes"][:7])))
    o.append(" ".join(map(str, c["sizes"][7:])))
    o += [" ".join(map(str, a)) for a in c["atY"]]
    o += [" ".join(map(str, a)) for a in c["atN"]]
    for t in c["tgs"]:
        o.append(" ".join(map(str, [t[0], t[1], 0] + t[2:7] + [rng_iw(t)])))
    o.append(" ".join(map(str, c["epl"])))
    o.append(" ".join(map(str, c["iso"])))
    o.append(" ".join(map(str, c["isolist"])))
    return "\n".join(o) + "\n"


def rng_iw(t):
    return (t[2] * 7 + t[4] * 3 + t[5]) % 5  # iWeight


def main():
    rng = random.Random(31337)
    cases = [gen(rng) for _ in range(3000)]
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, timeout=1800)
        assert r.returncode == 0, (r.returncode, r.stderr.decode()[:2000], r.stdout.decode()[-500:])
        out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_can2f.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
