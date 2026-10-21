"""Generates tests/data/inchi_can2c.json from the InChI C code: the
ranking set-up of ichican2.c (`FillOutAtomInvariant2`, `SetInitialRanks2`,
`CleanNumH`, `CleanCt`, `CleanIsoSortKeys`), compiled with the upstream
source included so that its types are available.

    python tools/inchi_can2c_reference.py
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
SRC = ma.SRC

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ichican2.c"
int CreateCompositeNormAtom() { return 0; }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        int n, ntg, nat, flags[6], hasArr, ignIso;
        scanf("%d %d %d", &n, &nat, &ntg);
        for (int i = 0; i < 6; i++) scanf("%d", &flags[i]);
        scanf("%d %d", &hasArr, &ignIso);
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        for (int i = 0; i < n; i++) {
            char el[16]; int v, ep, nh; long isk;
            scanf("%s %d %d %d %ld", el, &v, &ep, &nh, &isk);
            if (strcmp(el, "-")) strcpy(at[i].elname, el);
            at[i].valence = v; at[i].endpoint = ep; at[i].num_H = nh; at[i].iso_sort_key = isk;
        }
        T_GROUP *tg = calloc(ntg + 8, sizeof(T_GROUP));
        for (int k = 0; k < ntg; k++) {
            int ne; scanf("%d", &ne); tg[k].nNumEndpoints = ne;
            for (int j = 0; j < T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC; j++) { int x; scanf("%d", &x); tg[k].num[j] = x; }
        }
        T_GROUP_INFO tgi; memset(&tgi, 0, sizeof(tgi));
        tgi.t_group = hasArr ? tg : NULL; tgi.num_t_groups = ntg; tgi.bIgnoreIsotopic = ignIso;
        int withInfo = flags[5];
        ATOM_INVARIANT2 *inv = calloc(nat + 2, sizeof(ATOM_INVARIANT2));
        FillOutAtomInvariant2(at, n, nat, inv, flags[0], flags[1], flags[2], flags[3], flags[4], withInfo ? &tgi : NULL);
        printf("N %d\n", c);
        for (int i = 0; i < nat; i++) {
            printf("INV");
            for (int j = 0; j < AT_INV_LENGTH; j++) printf(" %d", inv[i].val[j]);
            printf(" %ld\n", (long) inv[i].iso_sort_key);
        }
        AT_RANK *rk = calloc(nat + 2, sizeof(AT_RANK)), *an = calloc(nat + 2, sizeof(AT_RANK));
        int nd = SetInitialRanks2(nat, inv, rk, an, &cg);
        printf("SR %d |", nd);
        for (int i = 0; i < nat; i++) printf(" %d", rk[i]);
        printf(" |");
        for (int i = 0; i < nat; i++) printf(" %d", an[i]);
        printf("\n");
        int m, len;
        scanf("%d", &m);
        NUM_H nh[40]; for (int i = 0; i < m; i++) { int x; scanf("%d", &x); nh[i] = x; }
        CleanNumH(nh, m);
        printf("CNH"); for (int i = 0; i < m; i++) printf(" %d", nh[i]); printf("\n");
        scanf("%d %d", &m, &len);
        AT_RANK ct[40]; for (int i = 0; i < m; i++) { int x; scanf("%d", &x); ct[i] = x; }
        int r = CleanCt(ct, len);
        printf("CCT %d", r); for (int i = 0; i < m; i++) printf(" %d", ct[i]); printf("\n");
        r = CleanCt(NULL, 0);
        printf("CCN %d\n", r);
        scanf("%d", &m);
        AT_ISO_SORT_KEY k[40]; for (int i = 0; i < m; i++) { long x; scanf("%ld", &x); k[i] = x; }
        CleanIsoSortKeys(k, m);
        printf("CIK"); for (int i = 0; i < m; i++) printf(" %ld", (long) k[i]); printf("\n");
    }
    return 0;
}
"""

ELS = ["C", "C", "C", "H", "D", "T", "N", "O", "Cl", "Br", "S", "P", "F", "Na", "Si", "-", "X", "Zn", "Se", "B"]
LONGMAX = 9223372036854775807


def gen(rng):
    n = rng.randrange(0, 12)
    ntg = rng.randrange(0, 4)
    nat = n + (rng.randrange(0, ntg + 3) if rng.random() < 0.6 else 0)
    c = {"n": n, "nat": nat, "ntg": ntg}
    c["flags"] = [rng.randrange(2) for _ in range(5)] + [1 if rng.random() < 0.85 else 0]
    c["hasArr"] = 1 if rng.random() < 0.9 else 0
    c["ign"] = rng.randrange(2)
    c["atoms"] = [
        [rng.choice(ELS), rng.randrange(0, 5), rng.randrange(0, ntg + 3), rng.randrange(0, 4),
         rng.choice([0, 0, 1, 5, 100, 4097, -3])]
        for _ in range(n)
    ]
    c["tg"] = [[rng.randrange(0, 6)] + [rng.randrange(0, 5) for _ in range(5)] for _ in range(ntg)]
    c["nh"] = [rng.choice([0x7FFE, 0x3FFE, 0x3FFF, 0x3FFE + 2, 3, 0]) for _ in range(rng.randrange(0, 10))]
    m = rng.randrange(1, 10)
    c["ct"] = [rng.choice([0, 0, 1, 7, 100]) for _ in range(m)]
    c["ctlen"] = rng.randrange(0, m)
    c["isk"] = [rng.choice([LONGMAX, 0, 5, -7, LONGMAX]) for _ in range(rng.randrange(0, 10))]
    return c


def fmt(c):
    out = [f"{c['n']} {c['nat']} {c['ntg']}", " ".join(map(str, c["flags"])), f"{c['hasArr']} {c['ign']}"]
    out += [" ".join(map(str, a)) for a in c["atoms"]]
    out += [" ".join(map(str, t)) for t in c["tg"]]
    out.append(" ".join(map(str, [len(c["nh"])] + c["nh"])))
    out.append(" ".join(map(str, [len(c["ct"]), c["ctlen"]] + c["ct"])))
    out.append(" ".join(map(str, [len(c["isk"])] + c["isk"])))
    return "\n".join(out) + "\n"


def main():
    rng = random.Random(9021)
    cases = [gen(rng) for _ in range(5000)]
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
    path = ROOT / "tests/data/inchi_can2c.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
