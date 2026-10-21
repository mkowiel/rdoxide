"""Generates tests/data/inchi_stereo.json from the InChI C code:
`set_stereo_parity` (with `set_stereo_atom_parity`,
`GetStereocenter0DParity`, `can_be_a_stereo_atom_with_isotopic_H`) of
ichister.c on random small molecules with stereocenters and stereo bonds.
The harness `#include`s ichister.c.

    python tools/inchi_stereo_reference.py
"""
import importlib.util
import json
import math
import random
import struct
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"tools/{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


sb = load("inchi_sbparity_reference")
ma = sb.ma
SRC = ma.SRC
bits = sb.bits

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define pow(x, y) ((x) * (x))
#include "ichister.c"
int CreateCompositeNormAtom() { return 0; }
static double rd(void) { unsigned long long u; double d; scanf("%llu", &u); memcpy(&d, &u, 8); return d; }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    CANON_GLOBALS cg;
    memset(&cg, 0, sizeof(cg));
    for (int c = 0; c < ncases; c++) {
        int n, nr, pe, vab, loose, zz; long mode;
        scanf("%d %d %d %d %d %d %ld", &n, &nr, &pe, &vab, &loose, &zz, &mode);
        int tot = n + nr;
        inp_ATOM *at = calloc(tot + 1, sizeof(inp_ATOM));
        for (int i = 0; i < tot; i++) {
            char el[8]; int ch, rad, nh, cf, v, orig, adiff, u0, cv, ep, cp, rs, nrs, pp, iso[3], sbp[3], sbo[3], sno[3], po[4];
            scanf("%s %d %d %d %d %d %d %d %d %d %d %d %d %d %d", el, &ch, &rad, &nh, &cf, &v, &orig, &adiff, &u0, &cv, &ep, &cp, &rs, &nrs, &pp);
            scanf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", iso, iso + 1, iso + 2, sbp, sbp + 1, sbp + 2, sbo, sbo + 1, sbo + 2, sno, sno + 1, sno + 2, po, po + 1, po + 2, po + 3);
            strcpy(at[i].elname, el);
            at[i].el_number = (U_CHAR) get_periodic_table_number(el);
            at[i].charge = ch; at[i].radical = rad; at[i].num_H = nh; at[i].cFlags = cf; at[i].valence = v;
            at[i].orig_at_number = orig; at[i].iso_atw_diff = adiff; at[i].bUsed0DParity = u0;
            at[i].chem_bonds_valence = cv; at[i].endpoint = ep; at[i].c_point = cp; at[i].nRingSystem = rs;
            at[i].nNumAtInRingSystem = nrs; at[i].p_parity = pp;
            at[i].x = rd(); at[i].y = rd(); at[i].z = rd();
            for (int j = 0; j < 3; j++) { at[i].num_iso_H[j] = iso[j]; at[i].sb_parity[j] = sbp[j]; at[i].sb_ord[j] = sbo[j]; at[i].sn_orig_at_num[j] = sno[j]; }
            for (int j = 0; j < 4; j++) at[i].p_orig_at_num[j] = po[j];
            for (int j = 0; j < v; j++) { int q, bt, bs; scanf("%d %d %d", &q, &bt, &bs); at[i].neighbor[j] = q; at[i].bond_type[j] = bt; at[i].bond_stereo[j] = bs; }
        }
        sp_ATOM *out = calloc(tot + 1, sizeof(sp_ATOM));
        int ma = -1, mb = -1;
        printf("N %d\n", c);
        int r = set_stereo_parity(&cg, at, out, n, nr, &ma, &mb, mode, pe, vab, loose, zz);
        printf("R %d %d %d\n", r, ma, mb);
        for (int i = 0; i < n; i++) {
            sp_ATOM *o = out + i;
            printf("A");
            for (int j = 0; j < 3; j++) printf(" %d %d %d %d %d %d %d %d", o->stereo_bond_neighbor[j], o->stereo_bond_neighbor2[j], o->stereo_bond_ord[j], o->stereo_bond_ord2[j], o->stereo_bond_z_prod[j], o->stereo_bond_z_prod2[j], o->stereo_bond_parity[j], o->stereo_bond_parity2[j]);
            printf(" %d %d %d %d %d %d %d %d | %d %d\n", o->parity, o->parity2, o->stereo_atom_parity, o->stereo_atom_parity2, o->bAmbiguousStereo, o->z_dir[0], o->z_dir[1], o->z_dir[2], at[i].bAmbiguousStereo, at[i].bUsed0DParity);
        }
        free(at); free(out);
    }
    return 0;
}
"""


def gen(rng):
    n = rng.randrange(3, 10)
    bonds = {}
    deg = [0] * n
    for i in range(1, n):
        for _ in range(10):
            p = rng.randrange(i)
            if deg[p] < 4:
                break
        else:
            return gen(rng)
        bonds[(p, i)] = 2 if rng.random() < 0.18 else 1
        deg[p] += 1
        deg[i] += 1
    if rng.random() < 0.3:
        a, b = sorted(rng.sample(range(n), 2))
        if (a, b) not in bonds and deg[a] < 4 and deg[b] < 4:
            bonds[(a, b)] = rng.choice([1, 1, 2, 4])
    nb = [[] for _ in range(n)]
    for (a, b), t in bonds.items():
        nb[a].append((b, t))
        nb[b].append((a, t))
    rs = sb.ring_systems(n, bonds)
    sizes = {}
    for x in rs:
        sizes[x] = sizes.get(x, 0) + 1
    flat2d = rng.random() < 0.5
    out = []
    for i in range(n):
        el = rng.choice(["C"] * 8 + ["N", "Si", "P", "S", "B", "O"])
        order = sum({1: 1, 2: 2, 4: 1}[t] for _, t in nb[i])
        tv = {"C": 4, "N": 3, "Si": 4, "P": 3, "S": 2, "B": 3, "O": 2}[el]
        ch = 0
        if el == "N" and rng.random() < 0.4:
            ch, tv = 1, 4
        if el == "P" and rng.random() < 0.3:
            tv = rng.choice([4, 5])
        if el == "S" and rng.random() < 0.4:
            tv = rng.choice([4, 6, 3])
            ch = 1 if tv == 3 else 0
        if rng.random() < 0.05:
            tv += rng.choice([-1, 1])
        nh = max(0, tv - order)
        cv = order + (1 if any(t == 4 for _, t in nb[i]) and rng.random() < 0.5 else 0)
        ang = rng.uniform(0, 6.28)
        z = 0.0 if flat2d else rng.uniform(-1, 1)
        xyz = [math.cos(i * 1.3 + rng.uniform(-0.3, 0.3)) * (1 + 0.3 * i), math.sin(i * 1.3 + rng.uniform(-0.3, 0.3)) * (1 + 0.2 * i), z]
        if rng.random() < 0.02:
            xyz = [0.0, 0.0, 0.0]
        stereo_bonds = rng.random() < 0.7
        nbl = []
        for q, t in nb[i]:
            bs = 0
            if flat2d and stereo_bonds and t == 1 and rng.random() < 0.5:
                bs = rng.choice([1, 1, 6, 6, 4, -1, -6])
            elif not flat2d and rng.random() < 0.05:
                bs = 4
            nbl.append([q, t, bs])
        iso = [0, 0, 0]
        if nh and rng.random() < 0.3:
            iso[rng.randrange(3)] = 1
            if rng.random() < 0.2:
                iso[rng.randrange(3)] += 1
        nbo = [q + 1 for q, _ in nb[i]]
        po = nbo + [i + 1] * (4 - len(nbo)) if rng.random() < 0.5 else nbo + [0] * (4 - len(nbo))
        rng.shuffle(po)
        out.append({"el": el, "ch": ch, "rad": rng.choice([0] * 20 + [1, 2]), "nh": nh, "cf": rng.choice([0] * 5 + [1]), "orig": i + 1, "adiff": 0,
                    "u0": 0, "cv": cv, "ep": rng.choice([0] * 25 + [1]), "cp": rng.choice([0] * 40 + [1]), "rs": rs[i], "nrs": sizes[rs[i]],
                    "pp": rng.choice([0, 0, 1, 2, 1, 2, 3, 4]) if rng.random() < 0.35 else 0,
                    "iso": iso, "sbp": [0] * 3, "sbo": [0] * 3, "sno": [0] * 3, "po": po, "xyz": xyz, "nb": nbl})
    if rng.random() < 0.25:
        for a in out:
            dbl = [j for j, (q, t, _) in enumerate(a["nb"]) if t == 2]
            if dbl and rng.random() < 0.8:
                a["u0"] = rng.choice([2, 2, 0])
                a["sbp"][0] = rng.choice([1, 2, 3, 4])
                a["sbo"][0] = rng.choice(dbl)
                a["sno"][0] = rng.choice([out[q]["orig"] for q, _, _ in a["nb"]] + [0])
    nr = rng.choice([0, 0, 0, 1, 2])
    for j in range(nr):
        tgt = rng.randrange(n)
        h = {"el": "H", "ch": 0, "rad": 0, "nh": 0, "cf": 0, "orig": 100 + j, "adiff": rng.choice([0, 0, 1, 2, 3]), "u0": 0, "cv": 1, "ep": 0, "cp": 0, "rs": 0, "nrs": 0, "pp": 0,
             "iso": [0] * 3, "sbp": [0] * 3, "sbo": [0] * 3, "sno": [0] * 3, "po": [0] * 4,
             "xyz": [rng.uniform(-1, 1), rng.uniform(-1, 1), rng.choice([0.0, 0.5])], "nb": [[tgt, 1, rng.choice([0, 0, 1, -1, 6, -6, 4])]]}
        out.append(h)
    for i in range(n):
        att = sum(1 for h in out[n:] if h["nb"][0][0] == i)
        out[i]["nh"] = max(out[i]["nh"], att)
        if att > 2:
            out[i]["nh"] = 3
    # at most 4 removed H on one atom is guaranteed (nr <= 2)
    return {"n": n, "nr": nr, "pe": rng.choice([0, 1, 1, 3, 7, 8, 15]), "vab": rng.choice([3, 3, 3, 6]), "loose": rng.choice([0, 1]),
            "zz": rng.choice([0, 1]), "mode": rng.choice([0, 5, 8]) << 16 | rng.choice([0, 0, 0x100]), "atoms": out}


def fmt(c):
    t = f"{c['n']} {c['nr']} {c['pe']} {c['vab']} {c['loose']} {c['zz']} {c['mode']}\n"
    for a in c["atoms"]:
        t += f"{a['el']} {a['ch']} {a['rad']} {a['nh']} {a['cf']} {len(a['nb'])} {a['orig']} {a['adiff']} {a['u0']} {a['cv']} {a['ep']} {a['cp']} {a['rs']} {a['nrs']} {a['pp']}\n"
        t += " ".join(map(str, a["iso"] + a["sbp"] + a["sbo"] + a["sno"] + a["po"])) + "\n"
        t += " ".join(str(bits(x)) for x in a["xyz"]) + " " + " ".join(f"{q} {b} {s}" for q, b, s in a["nb"]) + "\n"
    return t


def main():
    rng = random.Random(51)
    cases = [gen(rng) for _ in range(8000)]
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
    for c in cases:
        for a in c["atoms"]:
            a["xyz"] = [bits(x) for x in a["xyz"]]
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_stereo.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
