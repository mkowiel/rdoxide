"""Generates tests/data/inchi_make1.json from the InChI C code: `inp2spATOM`,
`set_atom_iso_sort_keys` and the Hill formula
comparisons (`CompareHillFormulas`, `CompareHillFormulasNoH`) of ichimake.c /
ichiisot.c. Links the upstream InChI base library (see
inchi_markalt_reference.py).

    python tools/inchi_make1_reference.py
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
int inp2spATOM( inp_ATOM *inp_at, int num_inp_at, sp_ATOM *at );
int CompareHillFormulas( const char *f1, const char *f2 );
int CompareHillFormulasNoH( const char *f1, const char *f2, int *num_H1, int *num_H2 );
AT_ISO_SORT_KEY make_iso_sort_key( int iso_atw_diff, int num_1H, int num_2H, int num_3H );
int set_atom_iso_sort_keys( int num_at, sp_ATOM *at, T_GROUP_INFO* t_group_info, int *bHasIsotopicInTautomerGroups );
static int sgn(int x) { return x < 0 ? -1 : x > 0; }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        char mode[4];
        scanf("%s", mode);
        printf("#CASE %d\n", c);
        if (mode[0] == 'A') {
            int n;
            scanf("%d", &n);
            inp_ATOM *ia = calloc(n + 1, sizeof(inp_ATOM));
            sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
            for (int i = 0; i < n; i++) {
                char el[16]; int v, k, t;
                scanf("%s %d", el, &v);
                strncpy(ia[i].elname, el, sizeof(ia[i].elname));
                ia[i].valence = v;
                for (int j = 0; j < v; j++) { scanf("%d %d", &k, &t); ia[i].neighbor[j] = k; ia[i].bond_type[j] = t; }
                int a[17];
                for (int q = 0; q < 17; q++) scanf("%d", a + q);
                ia[i].chem_bonds_valence = a[0]; ia[i].orig_at_number = a[1]; ia[i].orig_compt_at_numb = a[2];
                ia[i].endpoint = a[3]; ia[i].iso_atw_diff = a[4]; ia[i].num_H = a[5]; ia[i].cFlags = a[6];
                ia[i].num_iso_H[0] = a[7]; ia[i].num_iso_H[1] = a[8]; ia[i].num_iso_H[2] = a[9];
                ia[i].charge = a[10]; ia[i].radical = a[11]; ia[i].nBlockSystem = a[12]; ia[i].bCutVertex = a[13];
                ia[i].nRingSystem = a[14]; ia[i].nNumAtInRingSystem = a[15];
            }
            for (int i = 0; i < n; i++) memset(&at[i], 0x5a, sizeof(at[i]));
            inp2spATOM(ia, n, at);
            for (int i = 0; i < n; i++) {
                printf("%s %d %d |", at[i].elname, at[i].el_number, at[i].valence);
                for (int j = 0; j < at[i].valence; j++) printf(" %d %d", at[i].neighbor[j], at[i].bond_type[j]);
                printf(" | %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", at[i].chem_bonds_valence, at[i].orig_at_number,
                    at[i].orig_compt_at_numb, at[i].endpoint, at[i].iso_atw_diff, at[i].num_H, at[i].cFlags,
                    at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2], at[i].charge, at[i].radical,
                    at[i].nBlockSystem, at[i].bCutVertex, at[i].nRingSystem, at[i].nNumAtInRingSystem);
            }
        } else if (mode[0] == 'B') {
            int n, hastg, ng, hasiso, wantflag;
            scanf("%d %d %d %d %d", &n, &hastg, &ng, &hasiso, &wantflag);
            sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
            for (int i = 0; i < n; i++) {
                int a[7];
                for (int q = 0; q < 7; q++) scanf("%d", a + q);
                at[i].endpoint = a[0]; at[i].iso_atw_diff = a[1]; at[i].cFlags = a[2];
                at[i].num_iso_H[0] = a[3]; at[i].num_iso_H[1] = a[4]; at[i].num_iso_H[2] = a[5];
            }
            T_GROUP_INFO tgi; memset(&tgi, 0, sizeof(tgi));
            T_GROUP tg[2]; memset(tg, 0, sizeof(tg));
            AT_NUMB iso[2] = {1, 2};
            if (hastg) { tgi.t_group = tg; tgi.num_t_groups = ng; }
            if (hasiso) tgi.nIsotopicEndpointAtomNumber = iso;
            int flag = -7;
            int nis = set_atom_iso_sort_keys(n, at, hastg >= 0 ? &tgi : NULL, wantflag ? &flag : NULL);
            printf("%d %d\n", nis, flag);
            for (int i = 0; i < n; i++) {
                printf("%lld\n", (long long)at[i].iso_sort_key);
            }
        } else {
            char f1[64], f2[64];
            scanf("%s %s", f1, f2);
            if (f1[0] == '-') f1[0] = 0;
            if (f2[0] == '-') f2[0] = 0;
            int a = 0, b = 0, c1 = 0, c2 = 0;
            int r1 = CompareHillFormulas(f1, f2);
            int r2 = CompareHillFormulasNoH(f1, f2, &a, &b);
            printf("%d %d %d %d\n", sgn(r1), sgn(r2), a, b);
        }
    }
    return 0;
}
"""

ELS = ["C", "H", "N", "O", "S", "P", "Cl", "Br", "Na", "Fe", "Xx", "D"]


def formula(rng):
    if rng.random() < 0.05:
        return rng.choice(["-", "1C", "c2", "C2H"])
    k = rng.randrange(0, 5)
    parts = []
    for _ in range(k):
        e = rng.choice(ELS[:10])
        n = rng.choice(["", "", "2", "3", "10", "22", "100"])
        parts.append(e + n)
    return "".join(parts) or "-"


def main():
    rng = random.Random(11)
    cases = []
    for _ in range(600):
        n = rng.randrange(1, 8)
        toks = [f"A {n}"]
        for _i in range(n):
            v = rng.randrange(0, 5)
            el = rng.choice(ELS + ["Zz", "Uuo", "He"])
            t = f"{el} {v}"
            for _j in range(v):
                t += f" {rng.randrange(0, 100)} {rng.choice([1, 1, 2, 3, 4, 0x11, 0x21])}"
            t += " " + " ".join(str(x) for x in [
                rng.randrange(0, 9), rng.randrange(1, 60), rng.randrange(1, 60), rng.randrange(0, 4),
                rng.randrange(-5, 6), rng.randrange(0, 4), rng.choice([0, 1, 2]),
                rng.randrange(0, 3), rng.randrange(0, 3), rng.randrange(0, 3),
                rng.randrange(-2, 3), rng.randrange(0, 4), rng.randrange(0, 5), rng.randrange(0, 2),
                rng.randrange(0, 5), rng.randrange(0, 9), 0])
            toks.append(t)
        cases.append("\n".join(toks))
    for _ in range(1500):
        n = rng.randrange(0, 8)
        hastg = rng.choice([0, 1, 1, -1])
        ng = rng.choice([0, 1, 2])
        t = f"B {n} {hastg} {ng} {rng.choice([0, 1])} {rng.choice([0, 1])}"
        rows = []
        for _i in range(n):
            iso = rng.choice([0, 0, 0, 1, -1, 3, -4, 120, -120, 30])
            rows.append(" ".join(str(x) for x in [
                rng.choice([0, 0, 1, 2]), iso, rng.choice([0, 0, 1, 2, 3]),
                rng.choice([0, 0, 1, 2]), rng.choice([0, 0, 1, 31]), rng.choice([0, 0, 1, 2]), 0]))
        cases.append(t + "\n" + "\n".join(rows))
    for _ in range(3000):
        f1 = formula(rng)
        f2 = f1 if rng.random() < 0.2 else formula(rng)
        cases.append(f"D {f1} {f2}")
    stdin = f"{len(cases)}\n" + "\n".join(cases) + "\n"
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True)
    out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_make1.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
