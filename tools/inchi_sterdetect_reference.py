"""Generates tests/data/inchi_sterdetect.json from the InChI C code: the
stereo center / stereo bond candidate checks of ichister.c
(`bCanInpAtomBeAStereoCenter`, `bInpAtomHasRequirdNeigh`,
`can_be_a_stereo_bond_with_isotopic_H`, `bIsOxide`, ...). Like
inchi_stergeom_reference.py the harness `#include`s ichister.c.

    python tools/inchi_sterdetect_reference.py
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
#include <math.h>
#define pow(x, y) ((x) * (x))
#include "ichister.c"
int CreateCompositeNormAtom() { return 0; }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int kind;
        scanf("%d", &kind);
        printf("N %d\n", c);
        if (kind == 0) {
            int p, unk, iso, vab;
            scanf("%d %d %d %d", &p, &unk, &iso, &vab);
            printf("H %d\n", half_stereo_bond_action(p, unk, iso, vab));
            continue;
        }
        int n; long mode;
        scanf("%d %ld", &n, &mode);
        inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            char el[8]; int ch, rad, nh, ep, cv, v;
            scanf("%s %d %d %d %d %d %d", el, &ch, &rad, &nh, &ep, &cv, &v);
            strcpy(at[i].elname, el);
            at[i].el_number = (U_CHAR) get_periodic_table_number(el);
            at[i].charge = ch; at[i].radical = rad; at[i].num_H = nh; at[i].endpoint = ep;
            at[i].chem_bonds_valence = cv; at[i].valence = v;
            for (int j = 0; j < v; j++) { int q, bt; scanf("%d %d", &q, &bt); at[i].neighbor[j] = q; at[i].bond_type[j] = bt; }
        }
        for (int i = 0; i < n; i++) {
            for (int pe = 0; pe < 16; pe += 7) for (int zz = 0; zz < 2; zz++)
                printf("S %d\n", bCanInpAtomBeAStereoCenter(at, i, pe, zz));
            for (int req = 0; req < 4; req++) for (int nd = 0; nd < 3; nd++)
                printf("R %d\n", bInpAtomHasRequirdNeigh(at, i, req, nd, (req + nd) & 1));
            printf("E %d %d %d\n",
                   bCanAtomHaveAStereoBond(at[i].elname, at[i].charge, at[i].radical),
                   bCanAtomBeMiddleAllene(at[i].elname, at[i].charge, at[i].radical),
                   bCanAtomBeTerminalAllene(at[i].elname, at[i].charge, at[i].radical));
            printf("T %d\n", bIsSuitableHeteroInpAtom(at + i));
            for (int j = 0; j < at[i].valence; j++) printf("B %d\n", get_allowed_stereo_bond_type(at[i].bond_type[j]));
        }
        for (int i = 0; i < n; i++) {
            inp_ATOM *cp = malloc((n + 1) * sizeof(inp_ATOM));
            memcpy(cp, at, (n + 1) * sizeof(inp_ATOM));
            printf("O %d", bIsOxide(cp, i));
            for (int j = 0; j < cp[i].valence; j++) printf(" %d", cp[i].bond_type[j]);
            printf("\n");
            memcpy(cp, at, (n + 1) * sizeof(inp_ATOM));
            printf("C %d", can_be_a_stereo_bond_with_isotopic_H(cp, i, mode));
            for (int a = 0; a < n; a++) for (int j = 0; j < cp[a].valence; j++) printf(" %d", cp[a].bond_type[j]);
            printf("\n");
            free(cp);
        }
        free(at);
    }
    return 0;
}
"""

ELEMS = ["C", "C", "C", "C", "C", "C", "Si", "Ge", "Sn", "As", "B", "S", "Se", "N", "N", "P", "O", "O", "H", "Zz", "Zy", "F", "Cl"]


def gen_sb(rng):
    """a double (or alternating, tautomeric) bond between 2 atoms with few substituents"""
    n = rng.randrange(4, 8)
    nb = [[] for _ in range(n)]

    def edge(a, b, bt):
        nb[a].append((b, bt))
        nb[b].append((a, bt))
    edge(1, 0, rng.choice([2, 2, 2, 2, 4, 4, 8, 2 | 0x10, 1]))
    for k in range(2, n):
        a = rng.choice([0, 1] if k < 4 else range(k))
        if len(nb[a]) < 3:
            edge(a, k, rng.choice([1, 1, 1, 2, 4, 8, 9]))
    if rng.random() < 0.3:
        a, b = rng.sample(range(n), 2)
        if all(q != b for q, _ in nb[a]) and len(nb[a]) < 3 and len(nb[b]) < 3:
            edge(a, b, rng.choice([1, 2, 4]))
    atoms = []
    for i in range(n):
        el = rng.choice(["C", "C", "C", "N", "N", "Si", "O", "O", "S"])
        order = sum(2 if (bt & 0xF) == 2 else (3 if (bt & 0xF) == 3 else 1) for _, bt in nb[i])
        if rng.random() < 0.3:
            order = max(1, order - rng.choice([0, 1]))
        valence = len(nb[i])
        nh = max(0, rng.choice([3, 3, 2, 4]) - valence)
        atoms.append({"el": el, "ch": rng.choice([0] * 8 + [1]), "rad": rng.choice([0] * 12 + [1]), "nh": nh, "ep": 0,
                      "cv": order + rng.choice([0] * 6 + [1, -1]), "nb": nb[i]})
    return {"k": 1, "mode": rng.choice([0, 0, 0x100]), "atoms": atoms}


def gen(rng):
    if rng.random() < 0.4:
        return gen_sb(rng)
    if rng.random() < 0.1:
        return {"k": 0, "p": rng.randrange(-8, 9), "unk": rng.choice([0, 1]), "iso": rng.choice([0, 1]), "vab": rng.choice([3, 3, 3, 6])}
    n = rng.randrange(3, 8)
    edges = {}
    nb = [[] for _ in range(n)]
    for _ in range(rng.randrange(n - 1, 2 * n)):
        a, b = rng.sample(range(n), 2)
        if (min(a, b), max(a, b)) in edges or len(nb[a]) >= 4 or len(nb[b]) >= 4:
            continue
        bt = rng.choice([1, 1, 1, 1, 2, 2, 2, 4, 8, 9, 3, 1 | 0x10, 2 | 0x20, 4 | 0x40])
        edges[(min(a, b), max(a, b))] = bt
        nb[a].append((b, bt))
        nb[b].append((a, bt))
    atoms = []
    for i in range(n):
        el = rng.choice(ELEMS)
        order = sum((bt & 0xF) if (bt & 0xF) < 4 else (1 if (bt & 0xF) != 4 else 1) for _, bt in nb[i])
        nh = rng.choice([0, 0, 1, 2])
        if rng.random() < 0.6:
            nh = max(0, rng.choice([4, 4, 3, 3, 5, 6]) - len(nb[i]) - rng.choice([0, 0, 1]))
        cv = order if rng.random() < 0.85 else rng.choice([3, 4, 5, 6])
        atoms.append({"el": el, "ch": rng.choice([0, 0, 0, 0, 0, 1, -1]), "rad": rng.choice([0] * 8 + [1, 2]), "nh": nh,
                      "ep": rng.choice([0, 0, 0, 0, 1]), "cv": cv, "nb": nb[i]})
    return {"k": 1, "mode": rng.choice([0, 0x100]), "atoms": atoms}


def fmt(c):
    if c["k"] == 0:
        return f"0\n{c['p']} {c['unk']} {c['iso']} {c['vab']}\n"
    t = f"1\n{len(c['atoms'])} {c['mode']}\n"
    for a in c["atoms"]:
        t += f"{a['el']} {a['ch']} {a['rad']} {a['nh']} {a['ep']} {a['cv']} {len(a['nb'])} " + " ".join(f"{q} {b}" for q, b in a["nb"]) + "\n"
    return t


def main():
    rng = random.Random(21)
    cases = [gen(rng) for _ in range(6000)]
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-ffp-contract=off", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True)
        assert r.returncode == 0, r.stderr
        out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_sterdetect.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
