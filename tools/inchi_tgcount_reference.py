"""Generates tests/data/inchi_tgcount.json from the InChI C code:
`set_tautomer_iso_sort_keys`, `CountTautomerGroups` and
`SortTautomerGroupsAndEndpoints` of ichitaut.c on random t-group lists.
Links the upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_tgcount_reference.py
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
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichi_bns.h"
#include "ichicano.h"
#include "ichicant.h"
#include "ichinorm.h"
int CreateCompositeNormAtom() { return 0; }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, g, maxg, niso, flagsdone;
        unsigned long long norm;
        scanf("%d %d %d %d %d %llu", &n, &g, &maxg, &niso, &flagsdone, &norm);
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        for (int i = 0; i < n; i++) { int e; scanf("%d", &e); at[i].endpoint = e; }
        T_GROUP_INFO tgi;
        memset(&tgi, 0, sizeof(tgi));
        tgi.max_num_t_groups = maxg;
        tgi.num_t_groups = g;
        tgi.t_group = calloc(g + 3, sizeof(T_GROUP));
        for (int i = 0; i < g; i++) {
            int gn, ne, nums[5];
            scanf("%d %d %d %d %d %d %d", &gn, &ne, nums, nums + 1, nums + 2, nums + 3, nums + 4);
            tgi.t_group[i].nGroupNumber = gn; tgi.t_group[i].nNumEndpoints = ne;
            for (int q = 0; q < 5; q++) tgi.t_group[i].num[q] = nums[q];
        }
        tgi.nNumIsotopicEndpoints = niso;
        tgi.bTautFlagsDone = flagsdone ? (TG_FLAG_FOUND_ISOTOPIC_H_DONE) : 0;
        tgi.tni.bNormalizationFlags = norm;
        printf("N %d\n", c);
        int w = set_tautomer_iso_sort_keys(&tgi);
        printf("K %d", w);
        for (int i = 0; i < g; i++) printf(" %ld", (long)tgi.t_group[i].iWeight);
        printf("\n");
        int ret = CountTautomerGroups(at, n, &tgi);
        printf("R %d %d %d\n", ret, tgi.num_t_groups, tgi.nNumEndpoints);
        for (int i = 0; i < tgi.num_t_groups; i++) {
            printf("T %d %d %d |", tgi.t_group[i].nGroupNumber, tgi.t_group[i].nNumEndpoints, tgi.t_group[i].nFirstEndpointAtNoPos);
            for (int q = 0; q < 5; q++) printf(" %d", tgi.t_group[i].num[q]);
            printf("\n");
        }
        printf("P");
        for (int i = 0; i < n; i++) printf(" %d", at[i].endpoint);
        printf("\n");
        printf("E");
        if (tgi.nEndpointAtomNumber) for (int i = 0; i < tgi.nNumEndpoints; i++) printf(" %d", tgi.nEndpointAtomNumber[i]);
        printf("\n");
        printf("W");
        if (tgi.tGroupNumber) for (int i = 0; i < tgi.num_t_groups * TGSO_TOTAL_LEN; i++) printf(" %d", tgi.tGroupNumber[i]);
        printf("\n");
        int nat = n + tgi.num_t_groups;
        AT_RANK *rank = calloc(n + g + 2, sizeof(AT_RANK));
        for (int i = 0; i < n + g + 1; i++) { int r; scanf("%d", &r); rank[i] = r; }
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        int sr = SortTautomerGroupsAndEndpoints(&cg, &tgi, n, nat, rank);
        printf("S %d\n", sr);
        printf("E");
        if (tgi.nEndpointAtomNumber) for (int i = 0; i < tgi.nNumEndpoints; i++) printf(" %d", tgi.nEndpointAtomNumber[i]);
        printf("\n");
        printf("W");
        if (tgi.tGroupNumber) for (int i = 0; i < tgi.num_t_groups * TGSO_TOTAL_LEN; i++) printf(" %d", tgi.tGroupNumber[i]);
        printf("\n");
        free(at); free(rank);
    }
    return 0;
}
"""


def fmt(c):
    t = f"{len(c['ep'])} {len(c['groups'])} {c['maxg']} {c['niso']} {c['done']} {c['norm']}\n"
    t += " ".join(map(str, c["ep"])) + "\n"
    for g in c["groups"]:
        t += " ".join(map(str, g)) + "\n"
    t += " ".join(map(str, c["rank"])) + "\n"
    return t


def main():
    rng = random.Random(5)
    cases = []
    for _ in range(3000):
        n = rng.randrange(1, 14)
        g = rng.randrange(0, 5)
        ep = [rng.randrange(1, g + 1) if g and rng.random() < 0.45 else 0 for _ in range(n)]
        if rng.random() < 0.03 and g:
            ep[rng.randrange(n)] = g + 1  # out of range
        groups = []
        for i in range(g):
            cnt = sum(1 for e in ep if e == i + 1)
            if rng.random() < 0.08:
                cnt += rng.choice([-1, 1])
            gn = i + 1
            n0 = rng.randrange(0, 4)
            n1 = n0 if rng.random() < 0.2 else rng.randrange(0, 4)
            groups.append([gn, max(cnt, 0), n0, n1, rng.randrange(0, 3), rng.randrange(0, 3), rng.randrange(0, 3)])
        if g and rng.random() < 0.15:
            groups[rng.randrange(g)][1] = 0  # a group of another component
        maxg = g + rng.choice([1, 1, 1, 2]) if rng.random() > 0.05 else 0
        rank = [rng.randrange(1, 6) for _ in range(n + g + 1)]
        cases.append({"ep": ep, "groups": groups, "maxg": maxg, "niso": rng.choice([0, 0, 1, 2]), "done": rng.choice([0, 1]),
                      "norm": rng.choice([0, 0, 0, 1, 4, 0x80]), "rank": rank})
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True)
        if r.returncode == 0:
            out = r.stdout.decode().strip().split("\n")
        else:
            kept, lines = [], []
            for c in cases:
                r1 = subprocess.run([str(exe)], input=("1\n" + fmt(c)).encode(), capture_output=True)
                if r1.returncode == 0:
                    kept.append(c)
                    lines += [f"N {len(kept) - 1}"] + r1.stdout.decode().strip().split("\n")[1:]
            print(len(cases) - len(kept), "cases dropped")
            cases = kept
            out = lines
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_tgcount.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
