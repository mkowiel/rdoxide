"""Generates tests/data/inchi_termh.json from the InChI C code:
`remove_terminal_HDT` and `add_DT_to_num_H` of strutil.c (with
`cmp_iso_atw_diff_component_no` and `insertions_sort_AT_RANK`) on random atom
arrays with explicit H, D and T atoms and stereobond parities. The function
texts are cut out of the C sources. Needs a C compiler and the fetched InChI
sources."""

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
#include "extr_ct.h"
@FUNCS@
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, fix, add;
        scanf("%d %d %d", &n, &fix, &add);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, chg, iso, nh, ih0, ih1, ih2;
            char e0[8], e1[8];
            scanf("%d %d %s %s %d %d %d %d %d %d", &v, &cbv, e0, e1, &chg, &iso, &nh, &ih0, &ih1, &ih2);
            at[i].elname[0] = e0[0] == '_' ? 0 : e0[0];
            at[i].elname[1] = e1[0] == '_' ? 0 : e1[0];
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].charge = chg;
            at[i].iso_atw_diff = iso; at[i].num_H = nh;
            at[i].num_iso_H[0] = ih0; at[i].num_iso_H[1] = ih1; at[i].num_iso_H[2] = ih2;
            for (int j = 0; j < v; j++) { int q, b, s; scanf("%d %d %d", &q, &b, &s); at[i].neighbor[j] = q; at[i].bond_type[j] = b; at[i].bond_stereo[j] = s; }
            for (int m = 0; m < MAX_NUM_STEREO_BONDS; m++) {
                int a, b, p; scanf("%d %d %d", &a, &b, &p);
                at[i].sb_ord[m] = a; at[i].sn_ord[m] = b; at[i].sb_parity[m] = p;
            }
        }
        int r = remove_terminal_HDT(n, at, fix);
        if (add) add_DT_to_num_H(r > 0 ? r : n, at);
        printf("R %d\n", r);
        for (int i = 0; i < n; i++) {
            printf("A %d %d %d %d %d %d %d %d %d %d %d", at[i].elname[0], at[i].elname[1], at[i].valence, at[i].chem_bonds_valence, at[i].charge,
                   at[i].iso_atw_diff, at[i].num_H, at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2], at[i].component);
            for (int k = 0; k < at[i].valence; k++) printf(" %d:%d:%d", at[i].neighbor[k], at[i].bond_type[k], at[i].bond_stereo[k]);
            for (int m = 0; m < MAX_NUM_STEREO_BONDS; m++) printf(" S%d:%d:%d", at[i].sb_ord[m], at[i].sn_ord[m], at[i].sb_parity[m]);
            printf("\n");
        }
        free(at);
    }
    return 0;
}
"""


def gen(rng):
    n_heavy = rng.randrange(1, 6)
    atoms = []

    def new(name, chg=0, iso=0, nh=0):
        atoms.append({"e0": name[0], "e1": name[1] if len(name) > 1 else "_", "chg": chg, "iso": iso, "nh": nh,
                      "ih": [0, 0, 0], "nb": [], "sb": []})
        return len(atoms) - 1

    def link(a, b, t=1):
        atoms[a]["nb"].append((b, t, 0))
        atoms[b]["nb"].append((a, t, 0))

    heavy = [new(rng.choice(["C", "C", "N", "O", "Cl", "P"]), rng.choice([0, 0, 0, 1, -1])) for _ in range(n_heavy)]
    for _ in range(rng.randrange(0, n_heavy + 1)):
        a, b = rng.sample(heavy, 2) if n_heavy > 1 else (None, None)
        if a is not None and all(x != b for x, _, _ in atoms[a]["nb"]):
            link(a, b, rng.choice([1, 1, 2]))
    for h in heavy:
        for _ in range(rng.choice([0, 0, 1, 2, 3])):
            name = rng.choice(["H", "H", "H", "D", "T", "H"])
            hh = new(name, rng.choice([0, 0, 0, 1, -1]), rng.choice([0, 0, 1, 2, 3]) if name == "H" else 0, 0)
            if rng.random() < 0.05:
                atoms[hh]["nh"] = 1
            link(h, hh, 1)
    if rng.random() < 0.15:
        # an isolated pair of hydrogens: H-D, H-H, ...
        atoms[:] = []
        a = new(rng.choice(["H", "D", "T"]), rng.choice([0, 0, 1]), rng.choice([0, 0, 1, 2]))
        b = new(rng.choice(["H", "D", "T"]), rng.choice([0, 0, -1]), rng.choice([0, 0, 1, 2]))
        link(a, b)
    if rng.random() < 0.1:
        atoms.append({"e0": "_", "e1": "_", "chg": 0, "iso": 0, "nh": 0, "ih": [0, 0, 0], "nb": [], "sb": []})
    # shuffle atom order so that H precede their neighbors in some cases
    perm = list(range(len(atoms)))
    rng.shuffle(perm)
    inv = {old: new_ for new_, old in enumerate(perm)}
    out = []
    for new_, old in enumerate(perm):
        a = atoms[old]
        out.append({**a, "nb": [(inv[x], t, s) for x, t, s in a["nb"]], "sb": []})
    atoms = out
    # stereobonds: pairs of double-bonded atoms with parities
    for i, a in enumerate(atoms):
        for k, (j, t, _) in enumerate(a["nb"]):
            if i < j and t == 2 and len(a["sb"]) < 2 and len(atoms[j]["sb"]) < 2 and rng.random() < 0.8:
                kj = [x for x, _, _ in atoms[j]["nb"]].index(i)
                for (x, ko) in ((i, k), (j, kj)):
                    others = [q for q in range(len(atoms[x]["nb"])) if q != ko]
                    sn = rng.choice(others) if others else -1
                    atoms[x]["sb"].append([ko, sn, rng.choice([1, 2, 1, 2, 3])])
    for a in atoms:
        a["cbv"] = sum(t for _, t, _ in a["nb"])
    return atoms


def main():
    rng = random.Random(41)
    cases = []
    for _ in range(12000):
        cases.append({"atoms": gen(rng), "fix": rng.choice([0, 1]), "add": rng.choice([0, 1])})
    stdin = f"{len(cases)}\n"
    for c in cases:
        atoms = c["atoms"]
        stdin += f"{len(atoms)} {c['fix']} {c['add']}\n"
        for a in atoms:
            sb = a["sb"][:3] + [[0, 0, 0]] * (3 - len(a["sb"][:3]))
            stdin += (f"{len(a['nb'])} {a['cbv']} {a['e0']} {a['e1']} {a['chg']} {a['iso']} {a['nh']} {a['ih'][0]} {a['ih'][1]} {a['ih'][2]} "
                      + " ".join(f"{x} {t} {s}" for x, t, s in a["nb"]) + " " + " ".join(" ".join(map(str, e)) for e in sb) + "\n")
    text = (SRC / "strutil.c").read_text(errors="replace")
    sort_text = (SRC / "ichisort.c").read_text(errors="replace")

    def cut(t, start):
        a = t.rindex(start)
        return t[a:t.index("\n}\n", a) + 3]

    funcs = [cut(sort_text, "int insertions_sort_AT_RANK( AT_RANK *base, int num )\n{"),
             cut(text, "int cmp_iso_atw_diff_component_no( const void *a1, const void *a2 )\n{"),
             cut(text, "int add_DT_to_num_H( int num_atoms, inp_ATOM *at )\n"),
             cut(text, "int remove_terminal_HDT( int num_atoms, inp_ATOM *at, int bFixTermHChrg )\n{")]
    protos = """int insertions_sort_AT_RANK( AT_RANK *base, int num );
int cmp_iso_atw_diff_component_no( const void *a1, const void *a2 );
int add_DT_to_num_H( int num_atoms, inp_ATOM *at );
int remove_terminal_HDT( int num_atoms, inp_ATOM *at, int bFixTermHChrg );
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
        changed += lines[0] != f"R {n}"
        data.append({"in": c, "out": lines})
    path = ROOT / "tests/data/inchi_termh.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases,", changed, "with terminal H ->", path)


main()
