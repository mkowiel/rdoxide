"""Generates tests/data/inchi_cansparity.json from the InChI C code: the
sp_ATOM parity bookkeeping of ichicans.c / ichimap2.c
(`find_atoms_with_parity`, `Set/Remove...Stereo...`, `GetStereoBondParity`,
`HalfStereoBondParity`, `GetPermutationParity`, `GetStereoCenterParity`,
`GetStereoNeighborPos`). The functions are linked from the upstream library.

    python tools/inchi_cansparity_reference.py
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
#include "ichicant.h"
#include "ichicomn.h"
#include "ichister.h"
int SetHalfStereoBondIllDefPariy( sp_ATOM *at, int jn, int k1, int new_parity );
int RemoveHalfStereoBond( sp_ATOM *at, int jn, int k1 );
int CreateCompositeNormAtom() { return 0; }
static void dump(const char *tag, sp_ATOM *at, int n) {
    printf("%s", tag);
    for (int i = 0; i < n; i++) {
        printf(" |%d %d %d", at[i].parity, at[i].stereo_atom_parity, at[i].final_parity);
        for (int j = 0; j < 3; j++) printf(" %d %d %d %d", at[i].stereo_bond_neighbor[j], at[i].stereo_bond_ord[j], at[i].stereo_bond_z_prod[j], at[i].stereo_bond_parity[j]);
    }
    printf("\n");
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    CANON_GLOBALS cg;
    memset(&cg, 0, sizeof(cg));
    for (int c = 0; c < ncases; c++) {
        int n;
        scanf("%d", &n);
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        AT_RANK *rank = calloc(n + 1, sizeof(AT_RANK)), *crank = calloc(n + 1, sizeof(AT_RANK));
        for (int i = 0; i < n; i++) {
            int v, p, sap, fp, r, cr;
            scanf("%d %d %d %d %d %d", &v, &p, &sap, &fp, &r, &cr);
            at[i].valence = v; at[i].parity = p; at[i].stereo_atom_parity = sap; at[i].final_parity = fp;
            rank[i] = r; crank[i] = cr;
            for (int j = 0; j < v; j++) { int q; scanf("%d", &q); at[i].neighbor[j] = q; }
            for (int j = 0; j < 3; j++) { int a, b, z, pp; scanf("%d %d %d %d", &a, &b, &z, &pp); at[i].stereo_bond_neighbor[j] = a; at[i].stereo_bond_ord[j] = b; at[i].stereo_bond_z_prod[j] = z; at[i].stereo_bond_parity[j] = pp; }
        }
        printf("N %d\n", c);
        S_CHAR *vis = calloc(n + 1, 1);
        for (int i = 0; i < n; i++) {
            memset(vis, 0, n + 1);
            printf("c %d %d %d\n", GetStereoCenterParity(&cg, at, i, rank), find_atoms_with_parity(at, vis, -1, i), GetPermutationParity(&cg, at + i, (AT_RANK) (i % 3), crank));
            for (int j = 0; j < 3; j++) {
                printf("h %d", HalfStereoBondParity(at, i, j, rank));
                if (at[i].stereo_bond_neighbor[j]) {
                    int m = at[i].stereo_bond_neighbor[j] - 1;
                    printf(" %d %d", GetStereoNeighborPos(at, i, m), GetStereoBondParity(at, i, m, rank));
                }
                printf("\n");
            }
            printf("p %d\n", GetStereoNeighborPos(at, i, (i + 1) % n));
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                sp_ATOM *cp = malloc((n + 1) * sizeof(sp_ATOM));
                int r;
                memcpy(cp, at, (n + 1) * sizeof(sp_ATOM));
                r = at[i].stereo_bond_neighbor[j] ? SetOneStereoBondIllDefParity(cp, i, j, 1 + (i + j) % 4) : -9;
                printf("s %d", r); dump("", cp, n);
                memcpy(cp, at, (n + 1) * sizeof(sp_ATOM));
                r = at[i].stereo_bond_neighbor[j] ? RemoveOneStereoBond(cp, i, j) : -9;
                printf("r %d", r); dump("", cp, n);
                memcpy(cp, at, (n + 1) * sizeof(sp_ATOM));
                r = SetHalfStereoBondIllDefPariy(cp, i, j, 3 + j);
                printf("sh %d", r); dump("", cp, n);
                memcpy(cp, at, (n + 1) * sizeof(sp_ATOM));
                r = RemoveHalfStereoBond(cp, i, j);
                printf("rh %d", r); dump("", cp, n);
                free(cp);
            }
            sp_ATOM *cp = malloc((n + 1) * sizeof(sp_ATOM));
            memcpy(cp, at, (n + 1) * sizeof(sp_ATOM));
            printf("rc %d", RemoveOneStereoCenter(cp, i)); dump("", cp, n);
            free(cp);
        }
        free(at); free(rank); free(crank); free(vis);
    }
    return 0;
}
"""


def gen(rng):
    n = rng.randrange(2, 8)
    nb = [[] for _ in range(n)]
    for i in range(1, n):
        p = rng.randrange(i)
        nb[i].append(p)
        nb[p].append(i)
    for _ in range(rng.randrange(0, 3)):
        a, b = rng.sample(range(n), 2)
        if b not in nb[a]:
            nb[a].append(b)
            nb[b].append(a)
    nb = [x[:6] for x in nb]
    for i in range(n):
        nb[i] = [q for q in nb[i] if i in nb[q] or True]
    atoms = []
    for i in range(n):
        rng.shuffle(nb[i])
        atoms.append({"nb": list(nb[i]), "sb": [[0, 0, 0, 0] for _ in range(3)]})
    # stereo bonds: symmetric most of the time
    for i in range(n):
        for q in nb[i]:
            if q < i and rng.random() < 0.4 and i in nb[q]:
                zp = rng.choice([100, -100, 60, -60, 20, 0])
                par = rng.choice([1, 2, 3, 4, 6, 1 | 8, 2 | 8, 6 | 16, 1 | 0x40, 2 | 0x40, 3, 4])
                for a, b in ((i, q), (q, i)):
                    slot = next((s for s in range(3) if atoms[a]["sb"][s][0] == 0), None)
                    if slot is not None:
                        atoms[a]["sb"][slot] = [b + 1, nb[a].index(b) if rng.random() < 0.95 else rng.randrange(len(nb[a])), zp if a == i else rng.choice([zp, zp, -zp]), par if rng.random() < 0.8 else rng.choice([1, 2, 6, 3])]
    for i, a in enumerate(atoms):
        has_sb = a["sb"][0][0] != 0
        a["par"] = rng.choice([1, 2, 1, 2, 3, 4, 6, 0]) if (has_sb or rng.random() < 0.5) else 0
        a["sap"] = rng.choice([0, 6, a["par"], 1 | 0x40])
        a["fp"] = rng.choice([0, 1, 2, 3])
        a["rank"] = rng.choice([0, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2])
        a["crank"] = rng.choice([0, 1, 2, 3, 4, 5, 6, 7, 8])
    return {"atoms": atoms}


def fmt(c):
    t = f"{len(c['atoms'])}\n"
    for a in c["atoms"]:
        t += f"{len(a['nb'])} {a['par']} {a['sap']} {a['fp']} {a['rank']} {a['crank']} " + " ".join(map(str, a["nb"])) + "\n"
        t += " ".join(" ".join(map(str, s)) for s in a["sb"]) + "\n"
    return t


def main():
    rng = random.Random(61)
    cases = [gen(rng) for _ in range(5000)]
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True)
        assert r.returncode == 0, r.stderr
        out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_cansparity.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
