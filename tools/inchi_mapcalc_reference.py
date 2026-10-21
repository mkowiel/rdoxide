"""Generates tests/data/inchi_mapcalc.json from the InChI C code: the
`CalculatedPathsParitiesAreIdentical` of ichimap2.c, linked from the upstream
library.

    python tools/inchi_mapcalc_reference.py
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
int CalculatedPathsParitiesAreIdentical( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_RANK *nAtomNumberCanon, AT_RANK *nAtomNumberCanon1, AT_RANK *nAtomNumberCanon2, AT_RANK *nVisited1, AT_RANK *nVisited2, AT_RANK prev_sb_neigh, AT_RANK cur, AT_RANK next1, AT_RANK next2, int nNeighMode, int bParitiesInverted, int mode, CANON_STAT *pCS, int vABParityUnknown );
static void pa(AT_RANK *a, int n) { for (int i = 0; i < n; i++) printf(" %d", a[i]); printf(" |"); }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        int n;
        scanf("%d", &n);
        int nmax = n + 2;
        sp_ATOM *at = calloc(32770, sizeof(sp_ATOM));
        AT_RANK *crank = calloc(nmax + 8, 2), *r1 = calloc(nmax + 8, 2), *acanon = calloc(nmax + 8, 2);
        for (int i = 0; i < n; i++) {
            int v, p, sap, r, cr, nh, an, ac, rr2, aa2, rs;
            scanf("%d %d %d %d %d %d %d %d %d %d %d", &v, &p, &sap, &r, &cr, &nh, &an, &ac, &rr2, &aa2, &rs);
            at[i].valence = v; at[i].parity = p; at[i].stereo_atom_parity = sap; at[i].num_H = nh; at[i].nRingSystem = rs;
            r1[i] = r; crank[i] = cr; acanon[i] = ac;
            for (int j = 0; j < v; j++) { int q; scanf("%d", &q); at[i].neighbor[j] = q; }
            for (int j = 0; j < 3; j++) { int a, b, z, pp; scanf("%d %d %d %d", &a, &b, &z, &pp); at[i].stereo_bond_neighbor[j] = a; at[i].stereo_bond_ord[j] = b; at[i].stereo_bond_z_prod[j] = z; at[i].stereo_bond_parity[j] = pp; }
        }
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
        int nq; scanf("%d", &nq);
        for (int k = 0; k < nq; k++) {
            int prev, cur, nx1, nx2, nm, inv, mode;
            scanf("%d %d %d %d %d %d %d", &prev, &cur, &nx1, &nx2, &nm, &inv, &mode);
            AT_RANK *v1 = calloc(nmax, 2), *v2 = calloc(nmax, 2), *c1 = calloc(nmax, 2), *c2 = calloc(nmax, 2);
            for (int i = 0; i < n; i++) { int x; scanf("%d", &x); v1[i] = x; }
            for (int i = 0; i < n; i++) { int x; scanf("%d", &x); v2[i] = x; }
            int ret = CalculatedPathsParitiesAreIdentical(&cg, at, n, r1, crank, acanon, c1, c2, v1, v2, prev, cur, nx1, nx2, nm, inv, mode, &cs, 3);
            printf("Q %d |", ret); pa(v1, n); pa(v2, n); pa(c1, n); pa(c2, n); printf("\n");
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
    c["carb"] = [[rng.randrange(1, n + 1), rng.randrange(1, 5)] for _ in range(rng.randrange(0, 5))]
    c["dble"] = [[rng.randrange(1, n + 1), rng.randrange(1, n + 1), rng.randrange(1, 5)] for _ in range(rng.randrange(0, 4))]
    c["calls"] = []
    for _ in range(rng.randrange(1, 8)):
        cur = rng.randrange(n)
        nb = atoms[cur]["nb"]
        if rng.random() < 0.5:
            prev = 32767
        else:
            prev = rng.randrange(n)
        nx1 = rng.choice(nb) if nb else 0
        nx2 = rng.choice(nb) if nb else 0
        nm = rng.choice([1, 1, 2, 2, 0])
        inv = rng.choice([0, 0, 1, -1])
        mode = rng.choice([0, 1, 2, 3])
        full = inv != 0 or rng.random() < 0.3
        def vis():
            if rng.random() < 0.4:
                v = list(range(1, n + 1))
            else:
                v = list(range(1, n + 1))
                rng.shuffle(v)
            return [x if (full or rng.random() < 0.7) else 0 for x in v]
        c["calls"].append([prev, cur, nx1, nx2, nm, inv, mode, vis(), vis()])
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
    out.append(str(len(c["carb"])))
    out += [" ".join(map(str, x)) for x in c["carb"]]
    out.append(str(len(c["dble"])))
    out += [" ".join(map(str, x)) for x in c["dble"]]
    out.append(str(len(c["calls"])))
    for k in c["calls"]:
        out.append(" ".join(map(str, k[:7])))
        out.append(" ".join(map(str, k[7])))
        out.append(" ".join(map(str, k[8])))
    return "\n".join(out) + "\n"


def main():
    rng = random.Random(4114)
    cases = [gen(rng) for _ in range(6000)]
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
    path = ROOT / "tests/data/inchi_mapcalc.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
