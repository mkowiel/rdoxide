"""Generates tests/data/inchi_salts.json from the InChI C code: the ammonium
and metal salt functions of strutil.c (`bIsAmmoniumSalt`,
`DisconnectAmmoniumSalt`, `bIsMetalSalt`, `DisconnectMetalSalt`,
`DisconnectSalts`, `bIsMetalToDisconnect`, `bMayDisconnectMetals`) on random
atom arrays. The function texts are cut out of strutil.c. Needs a C compiler
and the fetched InChI sources."""

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
        int n, mode, nb0;
        scanf("%d %d %d", &n, &mode, &nb0);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, chg, rad, iso, nh, ih0, ih1, ih2;
            double x, y, z;
            scanf("%d %d %d %d %d %d %d %d %d %d %lf %lf %lf", &v, &cbv, &el, &chg, &rad, &iso, &nh, &ih0, &ih1, &ih2, &x, &y, &z);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].charge = chg;
            at[i].radical = rad; at[i].iso_atw_diff = iso; at[i].num_H = nh;
            at[i].num_iso_H[0] = ih0; at[i].num_iso_H[1] = ih1; at[i].num_iso_H[2] = ih2;
            at[i].x = x; at[i].y = y; at[i].z = z; at[i].orig_at_number = i + 1;
            for (int j = 0; j < v; j++) { int q, b, s; scanf("%d %d %d", &q, &b, &s); at[i].neighbor[j] = q; at[i].bond_type[j] = b; at[i].bond_stereo[j] = s; }
        }
        ORIG_ATOM_DATA od; memset(&od, 0, sizeof(od));
        od.at = at; od.num_inp_atoms = n; od.num_inp_bonds = nb0;
        if (mode == 0) {
            int r = DisconnectSalts(&od, 0);
            printf("R %d %d\n", r, od.num_inp_bonds);
        } else if (mode == 1) {
            int r = DisconnectSalts(&od, 1);
            printf("R %d %d\n", r, od.num_inp_bonds);
        } else {
            INCHI_MODE fl = 0;
            int r = bMayDisconnectMetals(&od, mode == 3, &fl);
            printf("R %d %d %d\n", r, od.bDisconnectCoord, (int)fl);
        }
        for (int i = 0; i < n; i++) {
            printf("A %d %d %d %d %d %d %d %d %d", at[i].valence, at[i].chem_bonds_valence, at[i].charge, at[i].num_H,
                   at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2], at[i].iso_atw_diff, at[i].radical);
            for (int k = 0; k < at[i].valence; k++) printf(" %d:%d:%d", at[i].neighbor[k], at[i].bond_type[k], at[i].bond_stereo[k]);
            printf("\n");
        }
        free(at);
    }
    return 0;
}
"""

METALS = [3, 11, 19, 12, 20, 30, 13, 26, 29, 37]


def new(atoms, el, chg=0, nh=0, rad=0, iso=0, ih=(0, 0, 0)):
    atoms.append({"el": el, "chg": chg, "rad": rad, "iso": iso, "nh": nh, "ih": list(ih), "nb": [],
                  "xyz": [rng_xyz() for _ in range(3)]})
    return len(atoms) - 1


RNG = random.Random(11)


def rng_xyz():
    return RNG.choice([0.0, 0.5, 1.0, 1.5, 2.0, -1.0, 3.25])


def link(atoms, a, b, t=1):
    atoms[a]["nb"].append((b, t, 0))
    atoms[b]["nb"].append((a, t, 0))


def gen(rng):
    atoms = []
    kind = rng.choice(["amm", "amm", "metal", "metal", "mix"])
    if kind in ("amm", "mix"):
        nitrogen = new(atoms, 7, rng.choice([0, 0, 1, 1, -1]), rng.choice([0, 0, 1, 2]))
        nh = 4
        # leave one slot for O or halogen
        want_h = rng.choice([3, 3, 3, 2, 4])
        impl = rng.choice([0, 0, 1, 2])
        atoms[nitrogen]["nh"] = min(impl, want_h)
        if rng.random() < 0.3:
            atoms[nitrogen]["ih"][rng.randrange(3)] = rng.choice([1, 1, 2])
        expl = max(0, want_h - atoms[nitrogen]["nh"] - sum(atoms[nitrogen]["ih"]))
        for _ in range(expl):
            h = new(atoms, 1, 0, iso=rng.choice([0, 0, 0, 1, 2, 3]))
            link(atoms, nitrogen, h)
        if rng.random() < 0.7:
            o = new(atoms, 8, rng.choice([0, 0, -1, -1, 1]))
            link(atoms, nitrogen, o)
            if rng.random() < 0.8:
                c = new(atoms, 6, 0, rng.choice([0, 0, 1]))
                link(atoms, o, c)
                for _ in range(rng.randrange(0, 3)):
                    x = new(atoms, rng.choice([6, 8, 1, 7]), 0)
                    link(atoms, c, x, rng.choice([1, 1, 2]))
        else:
            x = new(atoms, rng.choice([9, 17, 35, 53, 8]), rng.choice([0, 0, -1]))
            link(atoms, nitrogen, x)
    if kind in ("metal", "mix"):
        m = new(atoms, rng.choice(METALS), rng.choice([0, 0, 0, 1]), rng.choice([0, 0, 0, 1]))
        for _ in range(rng.randrange(1, 4)):
            if rng.random() < 0.75:
                o = new(atoms, 8)
                link(atoms, m, o)
                c = new(atoms, 6, 0, rng.choice([0, 0, 0, 1]))
                link(atoms, o, c)
                link(atoms, c, new(atoms, 8), 2)
                if rng.random() < 0.7:
                    link(atoms, c, new(atoms, 6), 1)
                if rng.random() < 0.15:
                    link(atoms, c, new(atoms, 1), 1)
            else:
                link(atoms, m, new(atoms, rng.choice([9, 17, 35, 53])))
    for a in atoms:
        if rng.random() < 0.05:
            a["rad"] = rng.choice([1, 2, 3])
    return atoms


def main():
    rng = random.Random(23)
    cases = []
    for _ in range(9000):
        atoms = gen(rng)
        for a in atoms:
            a["cbv"] = sum(t for _, t, _ in a["nb"])
            # fix the implicit-H tally of carbon-like atoms to a normal valence sometimes
        cases.append({"atoms": atoms, "mode": rng.choice([0, 1, 1, 2, 3]),
                      "nb": sum(len(a["nb"]) for a in atoms) // 2})
    stdin = f"{len(cases)}\n"
    for c in cases:
        atoms = c["atoms"]
        stdin += f"{len(atoms)} {c['mode']} {c['nb']}\n"
        for a in atoms:
            stdin += (f"{len(a['nb'])} {a['cbv']} {a['el']} {a['chg']} {a['rad']} {a['iso']} {a['nh']} "
                      f"{a['ih'][0]} {a['ih'][1]} {a['ih'][2]} {a['xyz'][0]} {a['xyz'][1]} {a['xyz'][2]} "
                      + " ".join(f"{x} {t} {s}" for x, t, s in a["nb"]) + "\n")
    text = (SRC / "strutil.c").read_text(errors="replace")

    def cut(start):
        a = text.rindex(start)
        return text[a:text.index("\n}\n", a) + 3]

    ister = (SRC / "ichister.c").read_text(errors="replace")
    a = ister.rindex("int get_opposite_sb_atom( inp_ATOM *at,\n")
    funcs = [ister[a:ister.index("\n}\n", a) + 3],
             cut("int RemoveInpAtBond( inp_ATOM *atom, int iat, int k )\n{"),
             cut("int bIsAmmoniumSalt( inp_ATOM *at,\n"),
             cut("int DisconnectAmmoniumSalt( inp_ATOM *at,\n"),
             cut("int bIsMetalSalt( inp_ATOM *at, int i )\n{"),
             cut("int DisconnectMetalSalt( inp_ATOM *at, int i )\n{"),
             cut("int DisconnectSalts( ORIG_ATOM_DATA *orig_inp_data, int bDisconnect )\n{"),
             cut("int bIsMetalToDisconnect( inp_ATOM *at, int i, int bCheckMetalValence )\n{"),
             cut("int bMayDisconnectMetals( ORIG_ATOM_DATA *orig_inp_data,\n")]
    protos = """int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt, int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );
int RemoveInpAtBond( inp_ATOM *atom, int iat, int k );
int bIsAmmoniumSalt( inp_ATOM *at, int i, int *piO, int *pk, S_CHAR *num_explicit_H );
int DisconnectAmmoniumSalt( inp_ATOM *at, int iN, int iO, int k, S_CHAR *num_explicit_H );
int bIsMetalSalt( inp_ATOM *at, int i );
int DisconnectMetalSalt( inp_ATOM *at, int i );
int DisconnectSalts( ORIG_ATOM_DATA *orig_inp_data, int bDisconnect );
int bIsMetalToDisconnect( inp_ATOM *at, int i, int bCheckMetalValence );
int bMayDisconnectMetals( ORIG_ATOM_DATA *orig_inp_data, int bCheckMetalValence, INCHI_MODE *bTautFlagsDone );
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
        changed += lines[0].split()[1] != "0"
        data.append({"in": c, "out": lines})
    path = ROOT / "tests/data/inchi_salts.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases,", changed, "with changes ->", path)


main()
