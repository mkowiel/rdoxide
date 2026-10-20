"""Generates tests/data/inchi_ionpairs.json from the InChI C code: `remove_ion_pairs`
and `fix_odd_things` of strutil.c (with the helpers they call) on random atom
arrays rich in charged N/O/C centres. The function texts are cut out of
strutil.c. Needs a C compiler and the fetched InChI sources."""

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
#include <string.h>
#include "mode.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
@FUNCS@
static void dump(int n, inp_ATOM *at) {
    for (int i = 0; i < n; i++) {
        printf("A %d %d %d %d", at[i].charge, at[i].chem_bonds_valence, at[i].radical, at[i].valence);
        for (int k = 0; k < at[i].valence; k++) printf(" %d:%d:%d", at[i].neighbor[k], at[i].bond_type[k], at[i].bond_stereo[k]);
        printf("\n");
    }
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, mode, fixbug, fixnu;
        scanf("%d %d %d %d", &n, &mode, &fixbug, &fixnu);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, chg, rad, iso, nh;
            scanf("%d %d %d %d %d %d %d", &v, &cbv, &el, &chg, &rad, &iso, &nh);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].charge = chg;
            at[i].radical = rad; at[i].iso_atw_diff = iso; at[i].num_H = nh;
            for (int j = 0; j < v; j++) { int x, b, s; scanf("%d %d %d", &x, &b, &s); at[i].neighbor[j] = x; at[i].bond_type[j] = b; at[i].bond_stereo[j] = s; }
        }
        int r = mode ? fix_odd_things(n, at, fixbug, fixnu) : remove_ion_pairs(n, at);
        printf("R %d\n", r);
        dump(n, at);
        free(at);
    }
    return 0;
}
"""

ELS = [6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 15, 16, 1, 17, 34]


def gen(rng):
    n = rng.randrange(2, 9)
    atoms = [{"el": rng.choice(ELS), "chg": rng.choice([0, 0, 1, 1, -1, -1, -1, 2]),
              "rad": rng.choice([0] * 8 + [1, 2, 3]), "iso": 0, "nh": rng.choice([0, 0, 0, 1, 2]), "nb": []}
             for _ in range(n)]
    # terminal oxygens / nitrogens around centres
    for i in range(n):
        if rng.random() < 0.5:
            for _ in range(rng.randrange(1, 4)):
                if len(atoms[i]["nb"]) < 6:
                    atoms.append({"el": rng.choice([8, 8, 8, 7, 6, 16]), "chg": rng.choice([0, 0, -1, -1, 1]),
                                  "rad": rng.choice([0] * 6 + [1, 2]), "iso": 0, "nh": rng.choice([0, 0, 0, 1]), "nb": []})
                    j = len(atoms) - 1
                    t = rng.choice([1, 1, 2, 2, 3])
                    atoms[i]["nb"].append((j, t, 0))
                    atoms[j]["nb"].append((i, t, 0))
    m = len(atoms)
    for _ in range(rng.randrange(1, 6)):
        a, b = rng.randrange(m), rng.randrange(m)
        if a != b and all(x != b for x, _, _ in atoms[a]["nb"]) and len(atoms[a]["nb"]) < 6 and len(atoms[b]["nb"]) < 6:
            t = rng.choice([1, 1, 1, 2, 3])
            st = rng.choice([0, 0, 1, -1, 4])
            atoms[a]["nb"].append((b, t, st))
            atoms[b]["nb"].append((a, t, -st))
    return atoms


def link(atoms, a, b, t):
    atoms[a]["nb"].append((b, t, 0))
    atoms[b]["nb"].append((a, t, 0))


def new(atoms, el, chg=0, nh=0, rad=0):
    atoms.append({"el": el, "chg": chg, "rad": rad, "iso": 0, "nh": nh, "nb": []})
    return len(atoms) - 1


def gen2(rng):
    """charge-separated motifs: A(+)-B(-), A(+)-X-B(-), A(+)(-B(-))-..."""
    atoms = []
    a = new(atoms, rng.choice([7, 7, 8, 6]), 1)
    b = new(atoms, rng.choice([8, 8, 7, 6]), -1, rng.choice([0, 0, 1]))
    shape = rng.choice(["pair", "pair", "bridge", "triple", "dpair"])
    if shape == "pair":
        link(atoms, a, b, rng.choice([1, 1, 2]))
    elif shape == "bridge":
        x = new(atoms, rng.choice([8, 7, 6]), 0, rng.choice([0, 0, 1]))
        link(atoms, a, x, rng.choice([1, 2]))
        link(atoms, x, b, rng.choice([1, 2]))
    elif shape == "triple":
        x = new(atoms, rng.choice([7, 7, 8, 6]), 0)
        b2 = new(atoms, rng.choice([6, 6, 7]), 1 if atoms[b]["chg"] < 0 else -1)
        link(atoms, a, x, rng.choice([1, 2]))
        link(atoms, x, b, rng.choice([1, 2]))
        link(atoms, x, b2, rng.choice([1, 2]))
        if rng.random() < 0.4:
            atoms[b2]["chg"] = 0
            atoms[b]["chg"] = 0
            atoms[a]["chg"] = 0
    else:
        link(atoms, a, b, rng.choice([1, 2]))
        o = new(atoms, 8, -1)
        link(atoms, a, o, 1)
    for i in range(len(atoms)):
        target = NORMAL.get((atoms[i]["el"], atoms[i]["chg"]), 4)
        while rng.random() < 0.6:
            d = target - sum(t for _, t, _ in atoms[i]["nb"])
            if d <= 0 or len(atoms[i]["nb"]) >= 5:
                break
            j = new(atoms, rng.choice([8, 8, 7, 6, 6, 16]), 0)
            t = rng.randint(1, min(d, 3))
            if NORMAL.get((atoms[j]["el"], 0), 4) < t:
                t = 1
            link(atoms, i, j, t)
    fill(rng, atoms)
    return atoms


NORMAL = {(6, 0): 4, (6, 1): 3, (6, -1): 3, (7, 0): 3, (7, 1): 4, (7, -1): 2, (8, 0): 2, (8, 1): 3, (8, -1): 1,
          (16, 0): 2, (1, 0): 1}


def fill(rng, atoms):
    """give most atoms their normal valence by implicit hydrogens"""
    for a in atoms:
        target = NORMAL.get((a["el"], a["chg"]))
        d = None if target is None else target - sum(t for _, t, _ in a["nb"])
        if d is not None and d > 0 and rng.random() < 0.85:
            a["nh"] = d if rng.random() < 0.9 else max(0, d - 1)


def main():
    rng = random.Random(77)
    cases = []
    for _ in range(12000):
        cases.append({"atoms": gen2(rng) if rng.random() < 0.8 else gen(rng), "mode": rng.choice([0, 1]), "fb": rng.choice([0, 1]), "fn": rng.choice([0, 1])})
    stdin = f"{len(cases)}\n"
    for c in cases:
        atoms = c["atoms"]
        stdin += f"{len(atoms)} {c['mode']} {c['fb']} {c['fn']}\n"
        for a in atoms:
            a["cbv"] = sum(t for _, t, _ in a["nb"])
            stdin += f"{len(a['nb'])} {a['cbv']} {a['el']} {a['chg']} {a['rad']} {a['iso']} {a['nh']} " + " ".join(f"{x} {t} {s}" for x, t, s in a["nb"]) + "\n"
    text = (SRC / "strutil.c").read_text(errors="replace")

    def cut(start):
        a = text.index(start)
        return text[a:text.index("\n}\n", a) + 3]

    funcs = [cut("int the_only_doublet_neigh( inp_ATOM *at,\n"),
             cut("int fix_non_uniform_drawn_oxoanions( int num_atoms,\n"),
             cut("int fix_non_uniform_drawn_amidiniums( int num_atoms,\n"),
             cut("int nFindOneOM( inp_ATOM *at, int at_no, int ord_OM[], int num_OM )\n{"),
             cut("int remove_ion_pairs( int num_atoms, inp_ATOM *at )\n{"),
             cut("int fix_odd_things( int num_atoms,\n")]
    protos = """int the_only_doublet_neigh( inp_ATOM *at, int i1, int *ineigh1, int *ineigh2 );
int fix_non_uniform_drawn_oxoanions( int num_atoms, inp_ATOM *at, int *num_changes );
int fix_non_uniform_drawn_amidiniums( int num_atoms, inp_ATOM *at, int *num_changes );
int nFindOneOM( inp_ATOM *at, int at_no, int ord_OM[], int num_OM );
int remove_ion_pairs( int num_atoms, inp_ATOM *at );
int fix_odd_things( int num_atoms, inp_ATOM *at, int bFixBug, int bFixNonUniformDraw );
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNCS@", protos + "\n".join(funcs)))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    k = 0
    data = []
    changed = 0
    for c in cases:
        n = len(c["atoms"])
        lines = out[k:k + n + 1]
        k += n + 1
        changed += lines[0] != "R 0"
        data.append({"in": c, "out": lines})
    path = ROOT / "tests/data/inchi_ionpairs.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases,", changed, "with changes ->", path)


main()
