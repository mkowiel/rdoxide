"""Generates tests/data/inchi_mapstereo4.json from the InChI C code: the
`map_stereo_bonds4` and `map_stereo_atoms4` of ichimap4.c, linked from the upstream
library.

    python tools/inchi_mapstereo4_reference.py
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
int map_stereo_bonds4( struct tagINCHI_CLOCK *ic, CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, int num_at_tg, int num_max, int bAllene, const AT_RANK *nCanonRankFrom, const AT_RANK *nAtomNumberCanonFrom, AT_RANK *nCanonRankTo, const AT_RANK *nSymmRank, AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, int nNumMappedRanksInput, AT_RANK *nSymmStereo, NEIGH_LIST *NeighList, CANON_STAT *pCS, CUR_TREE *cur_tree, int nNumMappedBonds, int vABParityUnknown );
static void pa(const char *t, AT_RANK *a, int n) { printf("%s", t); for (int i = 0; i < n; i++) printf(" %d", a[i]); printf("\n"); }
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
        AT_RANK *crank = calloc(nmax + 8, 2), *r1 = calloc(nmax + 8, 2), *acanon = calloc(nmax + 8, 2), *an = calloc(nmax + 8, 2);
        for (int i = 0; i < n; i++) {
            int v, p, sap, r, cr, nh, an_, ac, rr2, aa2, rs;
            scanf("%d %d %d %d %d %d %d %d %d %d %d", &v, &p, &sap, &r, &cr, &nh, &an_, &ac, &rr2, &aa2, &rs);
            at[i].valence = v; at[i].parity = p; at[i].stereo_atom_parity = sap; at[i].num_H = nh; at[i].nRingSystem = rs;
            r1[i] = r; crank[i] = cr; acanon[i] = ac; an[i] = an_;
            for (int j = 0; j < v; j++) { int q; scanf("%d", &q); at[i].neighbor[j] = q; }
            for (int j = 0; j < 3; j++) { int a, b, z, pp; scanf("%d %d %d %d", &a, &b, &z, &pp); at[i].stereo_bond_neighbor[j] = a; at[i].stereo_bond_ord[j] = b; at[i].stereo_bond_z_prod[j] = z; at[i].stereo_bond_parity[j] = pp; }
        }
        int bFirst, bKeep, nmode, nmri;
        scanf("%d %d %d %d", &bFirst, &bKeep, &nmode, &nmri);
        printf("N %d\n", c);
        CANON_STAT cs;
        memset(&cs, 0, sizeof(cs));
        int nc; scanf("%d", &nc);
        cs.LinearCTStereoCarb = calloc(nc + 1, sizeof(AT_STEREO_CARB));
        cs.nLenLinearCTStereoCarb = nc;
        for (int i = 0; i < nc; i++) { int a, p; scanf("%d %d", &a, &p); cs.LinearCTStereoCarb[i].at_num = a; cs.LinearCTStereoCarb[i].parity = p; }
        int nd; scanf("%d", &nd);
        cs.LinearCTStereoDble = calloc(nd + 1, sizeof(AT_STEREO_DBLE));
        cs.nLenLinearCTStereoDble = nd;
        for (int i = 0; i < nd; i++) { int a, b, p; scanf("%d %d %d", &a, &b, &p); cs.LinearCTStereoDble[i].at_num1 = a; cs.LinearCTStereoDble[i].at_num2 = b; cs.LinearCTStereoDble[i].parity = p; }
        cs.nSymmRank = r1;
        cs.bFirstCT = bFirst; cs.bKeepSymmRank = bKeep; cs.nMode = nmode;
        cs.bRankUsedForStereo = calloc(nmax + 8, 1);
        cs.bAtomUsedForStereo = calloc(nmax + 8, 1);
        cs.nPrevAtomNumber = calloc(nmax + 8, 2);
        NEIGH_LIST *nl = CreateNeighList(n, n, at, 0, NULL);
        AT_RANK *s1[64], *s2[64];
        memset(s1, 0, sizeof(s1)); memset(s2, 0, sizeof(s2));
        s1[0] = calloc(nmax, 2); s1[1] = calloc(nmax, 2); s2[0] = calloc(nmax, 2); s2[1] = calloc(nmax, 2);
        for (int i = 0; i < n; i++) { s1[0][i] = r1[i]; s1[1][i] = an[i]; s2[0][i] = r1[i]; s2[1][i] = an[i]; }
        AT_RANK *tmp = calloc(nmax + 8, 2);
        AT_RANK *canonTo = calloc(nmax + 8, 2), *symmStereo = calloc(nmax + 8, 2);
        for (int i = 0; i < n; i++) symmStereo[i] = i;
        CUR_TREE ct;
        memset(&ct, 0, sizeof(ct));
        CurTreeAlloc(&ct, n);
        int ret = map_stereo_bonds4(NULL, &cg, at, n, n, nmax, 0, crank, acanon, canonTo, r1, s1, s2, tmp, nmri, symmStereo, nl, &cs, &ct, 0, 3);
        printf("R %d %d %d %ld %ld %ld %ld %ld %ld %d %d\n", ret, cs.bFirstCT, cs.bStereoIsBetter, cs.lNumBreakTies, cs.lNumNeighListIter, cs.lNumTotCT, cs.lNumDecreasedCT, cs.lNumRejectedCT, cs.lNumEqualCT, cs.nLenLinearCTStereoCarb, cs.nLenLinearCTStereoDble);
        for (int i = 0; i < cs.nLenLinearCTStereoCarb; i++) printf("c %d %d\n", cs.LinearCTStereoCarb[i].at_num, cs.LinearCTStereoCarb[i].parity);
        for (int i = 0; i < cs.nLenLinearCTStereoDble; i++) printf("d %d %d %d\n", cs.LinearCTStereoDble[i].at_num1, cs.LinearCTStereoDble[i].at_num2, cs.LinearCTStereoDble[i].parity);
        for (int i = 0; i < n; i++) {
            printf("a %d %d %d |", at[i].parity, at[i].stereo_atom_parity, at[i].final_parity);
            for (int j = 0; j < 3; j++) printf(" %d,%d,%d", at[i].stereo_bond_neighbor[j], at[i].stereo_bond_parity[j], at[i].stereo_bond_ord[j]);
            printf("\n");
        }
        { printf("U"); for (int i = 0; i < n; i++) printf(" %d,%d", cs.bRankUsedForStereo[i], cs.bAtomUsedForStereo[i]); printf("\n"); }
        pa("T", canonTo, n); pa("S", symmStereo, n); pa("P", cs.nPrevAtomNumber, n);
        { printf("K %d %d %d |", ct.cur_len, ct.max_len, ct.incr_len); for (int i = 0; i < ct.cur_len; i++) printf(" %d", ct.tree[i]); printf("\n"); }
        { int k1 = 0, k2 = 0; while (s1[k1]) k1++; while (s2[k2]) k2++; printf("Z %d %d\n", k1, k2);
          for (int k = 0; k < k1; k++) pa(" s1", s1[k], n);
          for (int k = 0; k < k2; k++) pa(" s2", s2[k], n); }
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
    c["ctl"] = [rng.randrange(2), rng.randrange(2), rng.choice([0, 0, 0x80]), rng.randrange(1, len(c["atoms"]) + 1)]
    return c


def fmt2(c):
    lines = fmt(c).rstrip("\n").split("\n")
    n = len(c["atoms"])
    ctl = " ".join(map(str, c["ctl"]))
    return "\n".join(lines[: 1 + 2 * n] + [ctl] + lines[1 + 2 * n:]) + "\n"


def main():
    rng = random.Random(6116)
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
    path = ROOT / "tests/data/inchi_mapstereo4.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
