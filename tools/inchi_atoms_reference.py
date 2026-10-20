"""Generates tests/data/inchi_atoms.json from the InChI C code: the
`inp_ATOM` helpers of util.c (bonds to metals, explicit H count, ion
groups, `nNoMetal*`) on random atom arrays. Needs a C compiler and the
fetched InChI sources."""

import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src"

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include "mode.h"
#include "inpdef.h"
#include "util.h"
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int el = 0; el < 100; el++) printf("%d,%d ", ion_el_group(el), get_endpoint_valence(el));
    printf("\n");
    for (int c = 0; c < ncases; c++) {
        int n;
        scanf("%d", &n);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int el, ch, v, cbv, nh, i1, i2, i3;
            scanf("%d %d %d %d %d %d %d %d", &el, &ch, &v, &cbv, &nh, &i1, &i2, &i3);
            at[i].el_number = el; at[i].charge = ch; at[i].valence = v;
            at[i].chem_bonds_valence = cbv; at[i].num_H = nh;
            at[i].num_iso_H[0] = i1; at[i].num_iso_H[1] = i2; at[i].num_iso_H[2] = i3;
            for (int j = 0; j < v; j++) {
                int x, b;
                scanf("%d %d", &x, &b);
                at[i].neighbor[j] = x; at[i].bond_type[j] = b;
            }
        }
        for (int i = 0; i < n; i++) {
            int n0 = at[i].valence ? at[i].neighbor[0] : -1;
            int n1 = at[i].valence > 1 ? at[i].neighbor[1] : -1;
            printf("%d %d %d %d %d %d %d %d", nBondsValToMetal(at, i), num_of_H(at, i),
                   nNoMetalNumBonds(at, i), nNoMetalBondsValence(at, i), nNoMetalNeighIndex(at, i),
                   nNoMetalOtherNeighIndex(at, i, n0), nNoMetalOtherNeighIndex2(at, i, n0, n1),
                   n0 >= 0 ? has_other_ion_neigh(at, i, n0) : -1);
            printf(" %d\n", n0 >= 0 ? has_other_ion_in_sphere_2(at, i, n0) : -1);
        }
        free(at);
    }
    return 0;
}
"""

ELEMENTS = [1, 6, 6, 7, 7, 8, 8, 15, 16, 34, 52, 11, 26, 29, 12, 17]


def case(rng):
    n = rng.randrange(2, 10)
    nb = [[] for _ in range(n)]
    for i in range(1, n):
        j = rng.randrange(i)
        b = rng.choice([1, 1, 2, 3, 4])
        nb[i].append((j, b))
        nb[j].append((i, b))
    atoms = []
    for i in range(n):
        cbv = sum(b for _, b in nb[i]) + rng.randrange(0, 2)
        atoms.append({
            "el": rng.choice(ELEMENTS), "ch": rng.choice([0, 0, 0, 1, -1, 2]),
            "cbv": cbv, "h": [rng.randrange(0, 3), rng.choice([0, 0, 1]), rng.choice([0, 0, 1]), rng.choice([0, 0, 1])],
            "nb": nb[i]})
    return atoms


def main():
    rng = random.Random(5)
    cases = [case(rng) for _ in range(500)]
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{len(c)}\n"
        for a in c:
            stdin += f"{a['el']} {a['ch']} {len(a['nb'])} {a['cbv']} " + " ".join(map(str, a["h"])) + "\n"
            stdin += " ".join(f"{j} {b}" for j, b in a["nb"]) + "\n"
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS)
        srcs = [SRC / "util.c", SRC / "ichierr.c"]
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c")] + [str(s) for s in srcs], check=True)
        out = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True, check=True).stdout.strip().split("\n")
    groups = [[int(x) for x in p.split(",")] for p in out[0].split()]
    k = 1
    for c in cases:
        for a in c:
            a["out"] = [int(x) for x in out[k].split()]
            k += 1
    path = ROOT / "tests/data/inchi_atoms.json"
    path.write_text(json.dumps({"groups": groups, "cases": cases}, separators=(",", ":")))
    print(len(cases), "cases ->", path)


main()
