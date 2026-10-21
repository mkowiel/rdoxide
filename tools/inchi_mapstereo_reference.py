"""Generates tests/data/inchi_mapstereo.json from the InChI C code: the
stereo mapping helpers of ichimap1.c (`All_SC_Same`, `Next_SC_At_CanonRank2`,
`All_SB_Same`, `Next_SB_At_CanonRanks2`, `NextStereoParity2Test`,
`CompareLinCtStereo*`, `bUniqueAtNbrFromMappingRank`, `nGetMcr`,
`nJoin2Mcrs`, `SetUseAtomForStereo`, `CUR_TREE`), linked from the upstream
library.

    python tools/inchi_mapstereo_reference.py
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
int CreateCompositeNormAtom() { return 0; }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n;
        scanf("%d", &n);
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        AT_RANK *crank = calloc(n + 8, 2), *acanon = calloc(n + 8, 2);
        AT_RANK *r1 = calloc(n + 8, 2), *a1 = calloc(n + 8, 2), *r2 = calloc(n + 8, 2), *a2 = calloc(n + 8, 2);
        for (int i = 0; i < n; i++) {
            int v, p, sap, r, cr, nh, an, ac, rr2, aa2;
            scanf("%d %d %d %d %d %d %d %d %d %d", &v, &p, &sap, &r, &cr, &nh, &an, &ac, &rr2, &aa2);
            at[i].valence = v; at[i].parity = p; at[i].stereo_atom_parity = sap; at[i].num_H = nh;
            r1[i] = r; a1[i] = an; r2[i] = rr2; a2[i] = aa2; crank[i] = cr; acanon[i] = ac;
            for (int j = 0; j < v; j++) { int q; scanf("%d", &q); at[i].neighbor[j] = q; }
            for (int j = 0; j < 3; j++) { int a, b, z, pp; scanf("%d %d %d %d", &a, &b, &z, &pp); at[i].stereo_bond_neighbor[j] = a; at[i].stereo_bond_ord[j] = b; at[i].stereo_bond_z_prod[j] = z; at[i].stereo_bond_parity[j] = pp; }
        }
        AT_RANK *st1[2] = {r1, a1}, *st2[2] = {r2, a2};
        S_CHAR *used = calloc(n + 8, 1);
        SetUseAtomForStereo(used, at, n);
        int nmask; scanf("%d", &nmask);
        for (int i = 0; i < nmask; i++) { int k, v; scanf("%d %d", &k, &v); used[k] = v; }
        printf("N %d\n", c);
        printf("U"); for (int i = 0; i < n; i++) printf(" %d", used[i]); printf("\n");
        printf("A");
        for (int cr = 1; cr <= n; cr++) printf(" %d", All_SC_Same(cr, st1, st2, acanon, at));
        printf("\n");
        printf("P");
        for (int cr = 1; cr <= n; cr++) { AT_NUMB u = 999; int ok = bUniqueAtNbrFromMappingRank(st1, cr, &u); printf(" %d,%d", ok, ok ? u : -1); }
        printf("\n");
        {
            AT_RANK cr1 = 0, min = 0; int first = 1;
            printf("X");
            for (int it = 0; it < 12; it++) {
                int ok = Next_SC_At_CanonRank2(&cr1, &min, &first, used, st1, st2, acanon, n);
                printf(" %d,%d,%d,%d", ok, cr1, min, first);
                if (!ok) break;
                used[acanon[cr1 - 1] % (n ? n : 1)] = used[acanon[cr1 - 1] % (n ? n : 1)];
            }
            printf("\n");
        }
        printf("B");
        for (int a = 1; a <= n; a++) for (int b = 1; b <= n; b++) printf(" %d", All_SB_Same(a, b, st1, st2, acanon, at));
        printf("\n");
        for (int al = 0; al < 2; al++) {
            AT_RANK cr1 = 0, cr2 = 0, m1 = 0, m2 = 0; int first = 1;
            printf("Y%d", al);
            for (int it = 0; it < 12; it++) {
                int ok = Next_SB_At_CanonRanks2(&cr1, &cr2, &m1, &m2, &first, used, st1, st2, crank, acanon, at, n, al);
                printf(" %d,%d,%d,%d,%d,%d", ok, cr1, cr2, m1, m2, first);
                if (!ok) break;
            }
            printf("\n");
        }
        printf("Z");
        for (int sp = 0; sp < 7; sp++) for (int cp = 0; cp < 3; cp++) for (int m = 0; m < 8; m++) {
            int p = sp, q = cp;
            int r = NextStereoParity2Test(&p, &q, m & 1, (m >> 1) & 1, (m >> 2) & 1, (m + 1) & 1, (m + 3) & 1, 3);
            printf(" %d,%d,%d", r, p, q);
        }
        printf("\n");
        /* stereo CT comparisons */
        {
            int nd1, nd2, nc1, nc2, f;
            scanf("%d %d %d %d %d", &nd1, &nd2, &nc1, &nc2, &f);
            AT_STEREO_DBLE *d1 = calloc(nd1 + 1, sizeof(AT_STEREO_DBLE)), *d2 = calloc(nd2 + 1, sizeof(AT_STEREO_DBLE));
            AT_STEREO_CARB *c1 = calloc(nc1 + 1, sizeof(AT_STEREO_CARB)), *c2 = calloc(nc2 + 1, sizeof(AT_STEREO_CARB));
            for (int i = 0; i < nd1; i++) { int a, b, p; scanf("%d %d %d", &a, &b, &p); d1[i].at_num1 = a; d1[i].at_num2 = b; d1[i].parity = p; }
            for (int i = 0; i < nd2; i++) { int a, b, p; scanf("%d %d %d", &a, &b, &p); d2[i].at_num1 = a; d2[i].at_num2 = b; d2[i].parity = p; }
            for (int i = 0; i < nc1; i++) { int a, p; scanf("%d %d", &a, &p); c1[i].at_num = a; c1[i].parity = p; }
            for (int i = 0; i < nc2; i++) { int a, p; scanf("%d %d", &a, &p); c2[i].at_num = a; c2[i].parity = p; }
            AT_STEREO_DBLE *pd1 = (f & 1) ? NULL : d1, *pd2 = (f & 2) ? NULL : d2;
            AT_STEREO_CARB *pc1 = (f & 4) ? NULL : c1, *pc2 = (f & 8) ? NULL : c2;
            printf("K %d %d %d", CompareLinCtStereoDble(pd1, nd1, pd2, nd2), CompareLinCtStereoCarb(pc1, nc1, pc2, nc2), CompareLinCtStereo(pd1, nd1, pc1, nc1, pd2, nd2, pc2, nc2));
            if (nc1 > 0 && nd1 > 0) printf(" %d %d %d", CompareLinCtStereoAtomToValues(&c1[0], c2[0].at_num ? c2[0].at_num : 1, c2[0].parity), CompareLinCtStereoDoubleToValues(&d1[0], d2[0].at_num1, d2[0].at_num2, d2[0].parity), CompareLinCtStereoDoubleToValues(&d1[0], d1[0].at_num1, d1[0].at_num2, d1[0].parity));
            printf("\n");
        }
        /* equivalence classes */
        {
            AT_RANK *eq = calloc(n + 8, 2);
            for (int i = 0; i < n; i++) { int v; scanf("%d", &v); eq[i] = v; }
            int nj; scanf("%d", &nj);
            printf("E");
            for (int i = 0; i < n; i++) printf(" %d", nGetMcr(eq, i));
            printf("\n");
            for (int j = 0; j < nj; j++) { int a, b; scanf("%d %d", &a, &b); printf("J %d", nJoin2Mcrs(eq, a, b)); for (int i = 0; i < n; i++) printf(" %d", eq[i]); printf("\n"); }
        }
        /* CUR_TREE */
        {
            int nops; scanf("%d", &nops);
            CUR_TREE ct; memset(&ct, 0, sizeof(ct));
            int saved[64], nsaved = 0;
            AT_NUMB *symm = calloc(n + 8, 2);
            for (int i = 0; i < n; i++) { int v; scanf("%d", &v); symm[i] = v; }
            for (int o = 0; o < nops; o++) {
                int op, arg, arg2; scanf("%d %d %d", &op, &arg, &arg2);
                int r = 0;
                switch (op) {
                    case 0: r = CurTreeAlloc(&ct, arg > 2 ? arg : 3); nsaved = 0; break;
                    case 1: r = CurTreeAddRank(&ct, arg); break;
                    case 2: r = CurTreeAddAtom(&ct, arg % (n > 0 ? n : 1)); break;
                    case 3: { int tp = nsaved ? saved[arg % nsaved] : 0; CurTreeKeepLastAtomsOnly(&ct, tp, 1); r = tp; break; }
                    case 4: r = CurTreeRemoveIfLastAtom(&ct, arg % (n > 0 ? n : 1)); break;
                    case 5: r = CurTreeRemoveLastRank(&ct); break;
                    case 6: r = CurTreeRemoveLastRankIfNoAtoms(&ct); break;
                    case 7: if (nsaved) { r = CurTreeSetPos(&ct, saved[arg % nsaved]); } else r = -7; break;
                    case 8: r = CurTreeIsLastRank(&ct, arg); break;
                    case 9: if (ct.tree && ct.cur_len > 1 && ct.tree[ct.cur_len-1] - 1 > ct.cur_len - 1) { r = -9; break; } r = CurTreeIsLastAtomEqu(&ct, arg % (n > 0 ? n : 1), symm); break;
                    case 10: CurTreeFree(&ct); r = 0; nsaved = 0; break;
                    case 11: r = CurTreeGetPos(&ct); if (nsaved < 64) saved[nsaved++] = r; break;
                }
                if (ct.cur_len < 0) ct.cur_len = 0;
                printf("T %d %d %d %d:", op, r, ct.cur_len, ct.max_len);
                for (int i = 0; i < ct.cur_len && i < ct.max_len && ct.tree; i++) printf(" %d", ct.tree[i]);
                printf("\n");
            }
            CurTreeFree(&ct);
        }
        free(at);
    }
    return 0;
}
"""


def gen(rng):
    c = ks.gen(rng)
    atoms = c["atoms"]
    n = len(atoms)
    if rng.random() < 0.7:
        r2, a2 = [a["rank"] for a in atoms], [a["an"] for a in atoms]
    else:
        ncol = rng.choice([1, 2, 3, n])
        col = [rng.randrange(ncol) for _ in range(n)]
        r2 = [sum(1 for x in col if x <= col[i]) for i in range(n)]
        a2 = sorted(range(n), key=lambda i: (r2[i], i))
    for i, a in enumerate(atoms):
        a["rank2"], a["an2"] = r2[i], a2[i]
    mask = [[rng.randrange(n), rng.choice([0, 1, 2, 3, 8])] for _ in range(rng.randrange(0, 3))]
    nd1, nd2, nc1, nc2 = (rng.randrange(0, 4) for _ in range(4))
    base_d = [[rng.randrange(1, 5), rng.randrange(1, 5), rng.randrange(1, 4)] for _ in range(4)]
    base_c = [[rng.randrange(1, 5), rng.randrange(1, 4)] for _ in range(4)]
    d1 = [list(base_d[i]) for i in range(nd1)]
    d2 = [list(base_d[i]) if rng.random() < 0.7 else [rng.randrange(1, 5), rng.randrange(1, 5), rng.randrange(1, 4)] for i in range(nd2)]
    c1 = [list(base_c[i]) for i in range(nc1)]
    c2 = [list(base_c[i]) if rng.random() < 0.7 else [rng.randrange(1, 5), rng.randrange(1, 4)] for i in range(nc2)]
    eq = [rng.randrange(n) for _ in range(n)]
    eq = [min(eq[i], i) if rng.random() < 0.7 else i for i in range(n)]
    joins = [[rng.randrange(n), rng.randrange(n)] for _ in range(rng.randrange(0, 6))]
    ops = []
    for _ in range(rng.randrange(5, 40)):
        ops.append([rng.choice([0, 1, 1, 1, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11]), rng.randrange(0, 12), 0])
    if rng.random() < 0.8:
        ops.insert(0, [0, rng.choice([n, n * 2, 3]), 0])
    symm = [rng.randrange(1, 4) for _ in range(n)]
    return {"atoms": atoms, "mask": mask, "ct": [nd1, nd2, nc1, nc2, rng.randrange(16)], "d1": d1, "d2": d2, "c1": c1, "c2": c2,
            "eq": eq, "joins": joins, "ops": ops, "symm": symm}


def fmt(c):
    atoms = c["atoms"]
    t = f"{len(atoms)}\n"
    for a in atoms:
        t += f"{len(a['nb'])} {a['par']} {a['sap']} {a['rank']} {a['crank']} {a['nh']} {a['an']} {a['ac']} {a['rank2']} {a['an2']} " + " ".join(map(str, a["nb"])) + "\n"
        t += " ".join(" ".join(map(str, s)) for s in a["sb"]) + "\n"
    t += f"{len(c['mask'])} " + " ".join(f"{k} {v}" for k, v in c["mask"]) + "\n"
    t += " ".join(map(str, c["ct"])) + "\n"
    for key in ("d1", "d2", "c1", "c2"):
        for row in c[key]:
            t += " ".join(map(str, row)) + "\n"
    t += " ".join(map(str, c["eq"])) + "\n"
    t += f"{len(c['joins'])}\n" + "".join(f"{a} {b}\n" for a, b in c["joins"])
    t += f"{len(c['ops'])}\n" + " ".join(map(str, c["symm"])) + "\n" + "".join(" ".join(map(str, o)) + "\n" for o in c["ops"])
    return t


def main():
    rng = random.Random(111)
    cases = [gen(rng) for _ in range(3000)]
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True)
        assert r.returncode == 0, (r.returncode, r.stderr.decode()[:3000])
        out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_mapstereo.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
