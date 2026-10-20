"""Generates tests/data/inchi_util.json: results of the element helpers of
InChI's util.c (and the element table itself) for `crates/rdoxide/tests/inchi_utils.rs`.
Needs a C compiler and the fetched InChI sources."""

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
#include "util.h"
extern const int nElDataLen;
int el_number_in_internal_ref_table(const char *elname);
int main(void) {
    char op[8], buf[256];
    printf("LEN %d\n", nElDataLen);
    while (scanf("%7s", op) == 1) {
        if (!strcmp(op, "S")) {
            int n; char s[8]; scanf("%d", &n);
            int r = get_element_chemical_symbol(n, s);
            char t[8]; int r2 = get_element_or_pseudoelement_symbol(n, t);
            printf("%d %s %d %s\n", r, s, r2, t);
        } else if (!strcmp(op, "P")) {
            scanf("%s", buf);
            printf("%d %d %d %d\n", el_number_in_internal_ref_table(buf), get_periodic_table_number(buf), get_atomic_mass(buf), 0);
        } else if (!strcmp(op, "E")) {
            int n; scanf("%d", &n);
            printf("%d %d %d %d\n", get_atomic_mass_from_elnum(n), if_skip_add_H(n), get_el_type(n), is_el_a_metal(n));
        } else if (!strcmp(op, "V")) {
            int n, c, i; scanf("%d %d %d", &n, &c, &i);
            printf("%d\n", get_el_valence(n, c, i));
        } else if (!strcmp(op, "U")) {
            int n, c, r, bv, abv, nh, nb; scanf("%d %d %d %d %d %d %d", &n, &c, &r, &bv, &abv, &nh, &nb);
            printf("%d %d %d\n", get_unusual_el_valence(n, c, r, bv, nh, nb), needed_unusual_el_valence(n, c, r, bv, abv, nh, nb), detect_unusual_el_valence(n, c, r, bv, nh, nb));
        } else if (!strcmp(op, "C")) {
            scanf("%s", buf);
            int rad, ch;
            int r = extract_charges_and_radicals(buf, &rad, &ch);
            printf("%d %d %d [%s]\n", r, rad, ch, buf);
        } else if (!strcmp(op, "H")) {
            S_CHAR iso[3] = {0, 0, 0};
            scanf("%s", buf);
            int h = extract_H_atoms(buf, iso);
            printf("%d %d %d %d [%s]\n", h, iso[0], iso[1], iso[2], buf);
        } else if (!strcmp(op, "N")) {
            char el[8]; int inh, ch, rad, cbv, aiv, al, dn, hm; int iso[3]; S_CHAR siso[3];
            scanf("%s %d %d %d %d %d %d %d %d %d %d %d", el, &inh, &iso[0], &iso[1], &iso[2], &ch, &rad, &cbv, &aiv, &al, &dn, &hm);
            for (int i = 0; i < 3; i++) siso[i] = (S_CHAR)iso[i];
            printf("%d\n", get_num_H(el, inh, siso, ch, rad, cbv, aiv, al, dn, hm));
        }
    }
    return 0;
}
"""


def main():
    rng = random.Random(11)
    els = ["H", "D", "T", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si", "Fe", "Na", "Zz", "Zy", "Xx", "Og", "He", "Li", "Pt", "Ca"]
    lines = []
    for n in range(-2, 125):
        lines.append(f"S {n}")
    for e in els:
        lines.append(f"P {e}")
    for n in range(0, 121):
        lines.append(f"E {n}")
    for n in range(0, 121):
        for c in range(-3, 4):
            lines.append(f"V {n} {c} {rng.randrange(0, 6)}")
    for _ in range(3000):
        n = rng.choice([1, 5, 6, 7, 8, 9, 15, 16, 17, 35, 53, 14, 26, rng.randrange(1, 100)])
        lines.append("U %d %d %d %d %d %d %d" % (n, rng.randrange(-3, 4), rng.randrange(0, 4), rng.randrange(0, 8), rng.randrange(0, 8), rng.randrange(0, 5), rng.randrange(0, 6)))
    strs = ["C", "Na+", "Fe+2", "O--", "N-3", "Cl-", "C^", "C^^", "CH3.", "CH2..", "CH:", "O+-", "Fe+++2", "N+0", "X+2", "C..", "+", "-", "Ca+2.", "C^-", "H2O", "NH4+", "Mn2+", "S-2x", "O^^^", "Au+1"]
    for s in strs:
        lines.append(f"C {s}")
    hs = ["CH3", "CD3", "OH", "NH2", "H", "D", "T", "H2O", "CHD", "CH12", "CT2D", "pH4d", "Hg", "HD2T3", "CH3D", "SiH3", "Ho", "H2", "D2T", "Cl", "NaH", "CHH", "Dy", "Tc", "HeH", "C2H5"]
    for s in hs:
        lines.append(f"H {s}")
    for _ in range(3000):
        e = rng.choice(["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si", "Fe", "Na", "Zz", "Xx", "Se", "As", "Al", "Cu", "H"])
        lines.append("N %s %d %d %d %d %d %d %d %d %d %d %d" % (e, rng.randrange(0, 4), rng.randrange(0, 3), rng.randrange(0, 2), rng.randrange(0, 2), rng.randrange(-3, 4), rng.randrange(0, 4), rng.randrange(0, 8), rng.choice([0, 0, 0, 3, 4, 15]), rng.randrange(2) , rng.randrange(2), rng.randrange(2)))
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS)
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c")], check=True)
        out = subprocess.run([str(exe)], input="\n".join(lines) + "\n", capture_output=True, text=True, check=True).stdout.split("\n")
    assert out[0].startswith("LEN")
    results = [o for o in out[1:] if o != ""]
    assert len(results) == len(lines), (len(results), len(lines))
    data = {"len": int(out[0].split()[1]), "cases": [{"in": i, "out": o} for i, o in zip(lines, results)]}
    path = ROOT / "tests/data/inchi_util.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(lines), "cases ->", path)


main()
