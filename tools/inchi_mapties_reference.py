"""Generates tests/data/inchi_mapties.json from the InChI C code: the
mapping functions of ichimap2.c (`NumberOfTies`, `map_an_atom2`,
`parity_of_mapped_atom2`, `parity_of_mapped_half_bond`,
`might_change_other_atom_parity`, `GetMinNewRank`), linked from the upstream
library.

    python tools/inchi_mapties_reference.py
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
int NumberOfTies( AT_RANK **pRankStack1, AT_RANK **pRankStack2, int length, int at_no1, int at_no2, AT_RANK *nNewRank, int *bAddStack, int *bMapped1 );
int parity_of_mapped_half_bond( int from_at, int to_at, int from_neigh, int to_neigh, sp_ATOM *at, EQ_NEIGH *pEN, const AT_RANK *nCanonRankFrom, const AT_RANK *nRankFrom, const AT_RANK *nRankTo );
int map_an_atom2( CANON_GLOBALS *pCG, int num_atoms, int num_max, int at_no1, int at_no2, AT_RANK *nTempRank, int nNumMappedRanks, int *pnNewNumMappedRanks, CANON_STAT *pCS, NEIGH_LIST *NeighList, AT_RANK **pRankStack1, AT_RANK **pRankStack2, int *bAddStack );
int might_change_other_atom_parity( sp_ATOM *at, int num_atoms, int at_no, AT_RANK *nRank2, AT_RANK *nRank1 );
AT_RANK GetMinNewRank( AT_RANK *nAtomRank, AT_RANK *nAtomNumb, AT_RANK nRank1 );
static void pen(EQ_NEIGH *e) { printf(" | %d %d %d %d %d %d %d %d", e->num_to, e->to_at[0], e->to_at[1], e->to_at[2], e->to_at[3], e->from_at, e->rank, e->canon_rank); }
static void pa(const char *t, AT_RANK *a, int n) { printf("%s", t); for (int i = 0; i < n; i++) printf(" %d", a[i]); }
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
        AT_RANK *crank = calloc(nmax + 8, 2), *acanon = calloc(nmax + 8, 2);
        AT_RANK *r1 = calloc(nmax + 8, 2), *a1 = calloc(nmax + 8, 2), *r2 = calloc(nmax + 8, 2), *a2 = calloc(nmax + 8, 2);
        for (int i = 0; i < n; i++) {
            int v, p, sap, r, cr, nh, an, ac, rr2, aa2, hs;
            scanf("%d %d %d %d %d %d %d %d %d %d %d", &v, &p, &sap, &r, &cr, &nh, &an, &ac, &rr2, &aa2, &hs);
            at[i].valence = v; at[i].parity = p; at[i].stereo_atom_parity = sap; at[i].num_H = nh; at[i].bHasStereoOrEquToStereo = hs;
            r1[i] = r; a1[i] = an; r2[i] = rr2; a2[i] = aa2; crank[i] = cr; acanon[i] = ac;
            for (int j = 0; j < v; j++) { int q; scanf("%d", &q); at[i].neighbor[j] = q; }
            for (int j = 0; j < 3; j++) { int a, b, z, pp; scanf("%d %d %d %d", &a, &b, &z, &pp); at[i].stereo_bond_neighbor[j] = a; at[i].stereo_bond_ord[j] = b; at[i].stereo_bond_z_prod[j] = z; at[i].stereo_bond_parity[j] = pp; }
        }
        printf("N %d\n", c);
        int np; scanf("%d", &np);
        for (int k = 0; k < np; k++) {
            int f, t, kf, kt, useen;
            scanf("%d %d %d %d %d", &f, &t, &kf, &kt, &useen);
            EQ_NEIGH e1, e2;
            memset(&e1, 0x55, sizeof(e1)); memset(&e2, 0x55, sizeof(e2));
            int ra = parity_of_mapped_atom2(&cg, f, t, at, useen ? &e1 : NULL, crank, r1, r2);
            printf("A %d", ra); if (useen) pen(&e1); printf("\n");
            int rh = parity_of_mapped_half_bond(f, t, at[f].neighbor[kf], at[t].neighbor[kt], at, useen ? &e2 : NULL, crank, r1, r2);
            printf("H %d", rh); if (useen) pen(&e2); printf("\n");
            ra = parity_of_mapped_atom2(&cg, f, f, at, NULL, crank, r1, r1);
            rh = parity_of_mapped_half_bond(f, f, at[f].neighbor[kf], at[f].neighbor[kf], at, NULL, crank, r1, r1);
            printf("S %d %d\n", ra, rh);
        }
        printf("G");
        for (int i = 0; i < n; i++) printf(" %d %d", GetMinNewRank(r1, a1, r1[i]), GetMinNewRank(r2, a2, r2[i]));
        printf("\n");
        printf("M");
        for (int i = 0; i < n; i++) printf(" %d", might_change_other_atom_parity(at, n, i, r2, r1));
        printf("\n");
        NEIGH_LIST *nl = CreateNeighList(n, n, at, 0, NULL);
        AT_RANK *st1[8], *st2[8];
        memset(st1, 0, sizeof(st1)); memset(st2, 0, sizeof(st2));
        st1[0] = r1; st1[1] = a1; st2[0] = r2; st2[1] = a2;
        for (int i = 2; i < 6; i++) st1[i] = calloc(nmax, 2);
        for (int i = 2; i < 4; i++) st2[i] = calloc(nmax, 2);
        AT_RANK *tmp = calloc(nmax + 8, 2);
        CANON_STAT cs;
        memset(&cs, 0, sizeof(cs));
        int nm; scanf("%d", &nm);
        for (int k = 0; k < nm; k++) {
            int x, y, nmr;
            scanf("%d %d %d", &x, &y, &nmr);
            AT_RANK nr; int badd, bm;
            int rt = NumberOfTies(st1, st2, nmax * 2, x, y, &nr, &badd, &bm);
            printf("T %d %d %d %d", rt, nr, badd, bm);
            int newn = -1;
            int rr = map_an_atom2(&cg, n, nmax, x, y, tmp, nmr, &newn, &cs, nl, st1, st2, &badd);
            printf("\nR %d %d %d %ld %ld", rr, newn, badd, cs.lNumNeighListIter, cs.lNumBreakTies);
            pa(" |", st1[2], n); pa(" |", st1[3], n); pa(" |", st2[2], n); pa(" |", st2[3], n); printf(" | %d %d", st1[4][0], st1[5][0]);
            printf("\n");
        }
    }
    return 0;
}
"""


def gen(rng):
    c = ms.gen(rng)
    n = len(c["atoms"])
    c["hs"] = [rng.choice([0, 1, 1]) for _ in range(n)]
    c["pairs"] = []
    for _ in range(rng.randrange(1, 10)):
        f = rng.randrange(n)
        cand = [i for i in range(n) if c["atoms"][i]["rank2"] == c["atoms"][f]["rank"]] if rng.random() < 0.8 else []
        t = rng.choice(cand) if cand else rng.randrange(n)
        c["pairs"].append([f, t, rng.randrange(4), rng.randrange(4), rng.randrange(2)])
    c["maps"] = []
    for _ in range(rng.randrange(1, 6)):
        x = rng.randrange(n)
        cand = [i for i in range(n) if c["atoms"][i]["rank2"] == c["atoms"][x]["rank"]]
        y = rng.choice(cand) if cand and rng.random() < 0.85 else rng.randrange(n)
        c["maps"].append([x, y, rng.randrange(1, n + 1)])
    return c


def fmt(c):
    n = len(c["atoms"])
    lines = ms.fmt(c).split("\n")
    out = [lines[0]]
    for i in range(n):
        toks = lines[1 + 2 * i].split(" ")
        # "v par sap rank crank nh an ac rank2 an2 nb... sb..." -> insert the flag after an2
        toks.insert(10, str(c["hs"][i]))
        out.append(" ".join(toks))
        out.append(lines[2 + 2 * i])
    out.append(str(len(c["pairs"])))
    out += [" ".join(map(str, p)) for p in c["pairs"]]
    out.append(str(len(c["maps"])))
    out += [" ".join(map(str, m)) for m in c["maps"]]
    return "\n".join(out) + "\n"


def main():
    rng = random.Random(2112)
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
    path = ROOT / "tests/data/inchi_mapties.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
