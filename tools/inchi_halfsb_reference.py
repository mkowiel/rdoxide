"""Generates tests/data/inchi_halfsb.json from the InChI C code:
`half_stereo_bond_parity` (and with it `GetHalfStereobond0DParity`) of
ichister.c. The harness `#include`s ichister.c.

    python tools/inchi_halfsb_reference.py
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
static double rd(void) { unsigned long long u; double d; scanf("%llu", &u); memcpy(&d, &u, 8); return d; }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, nr, pe, vab, usez, zi[3];
        scanf("%d %d %d %d %d %d %d %d", &n, &nr, &pe, &vab, &usez, zi, zi + 1, zi + 2);
        inp_ATOM *at = calloc(n + nr + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n + nr; i++) {
            char el[8]; int ch, rad, nh, cf, v, orig, adiff, u0, iso[3], sbp[3], sbo[3], sno[3];
            scanf("%s %d %d %d %d %d %d %d %d %d %d %d", el, &ch, &rad, &nh, &cf, &v, &orig, &adiff, &u0, iso, iso + 1, iso + 2);
            scanf("%d %d %d %d %d %d %d %d %d", sbp, sbp + 1, sbp + 2, sbo, sbo + 1, sbo + 2, sno, sno + 1, sno + 2);
            strcpy(at[i].elname, el);
            at[i].charge = ch; at[i].radical = rad; at[i].num_H = nh; at[i].cFlags = cf; at[i].valence = v;
            at[i].orig_at_number = orig; at[i].iso_atw_diff = adiff; at[i].bUsed0DParity = u0;
            at[i].x = rd(); at[i].y = rd(); at[i].z = rd();
            for (int j = 0; j < 3; j++) { at[i].num_iso_H[j] = iso[j]; at[i].sb_parity[j] = sbp[j]; at[i].sb_ord[j] = sbo[j]; at[i].sn_orig_at_num[j] = sno[j]; }
            for (int j = 0; j < v; j++) { int q, bs; scanf("%d %d", &q, &bs); at[i].neighbor[j] = q; at[i].bond_stereo[j] = bs; }
        }
        S_CHAR z[3] = {zi[0], zi[1], zi[2]};
        int r = half_stereo_bond_parity(at, 0, at + n, nr, usez ? z : NULL, pe, vab);
        printf("%d %d %d %d %d %d %d\n", r, usez ? z[0] : 0, usez ? z[1] : 0, usez ? z[2] : 0, at[0].bAmbiguousStereo, at[0].bUsed0DParity, 0);
        free(at);
    }
    return 0;
}
"""


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def gen(rng):
    nn = rng.randrange(1, 4)  # neighbors of the center
    nr = rng.choice([0, 0, 1, 1, 2, 3])
    n = 1 + nn + 2  # center, neighbors, extra atoms for the stereo-bond partner
    planar = rng.random() < 0.6
    flat = rng.random() < 0.15

    def xyz(ang=None):
        if flat:
            return [0.0, 0.0, 0.0]
        a = ang if ang is not None else rng.uniform(0, 2 * math.pi)
        r = rng.uniform(0.5, 1.6)
        z = 0.0 if planar else rng.choice([0.0, 0.5, -0.5, 1.0])
        return [r * math.cos(a), r * math.sin(a), z]

    def atom(el, nb, orig, c=None):
        return {"el": el, "ch": rng.choice([0] * 9 + [1]), "rad": rng.choice([0] * 12 + [1]),
                "nh": rng.choice([0, 0, 1, 1, 2, 3]), "cf": rng.choice([0, 0, 0, 1]), "orig": orig,
                "adiff": 0, "u0": 0, "iso": [rng.choice([0, 0, 0, 1, 2]) for _ in range(3)],
                "sbp": [0, 0, 0], "sbo": [0, 0, 0], "sno": [0, 0, 0],
                "xyz": c if c else [0.0, 0.0, 0.0], "nb": nb}

    bs = lambda: rng.choice([0] * 4 + [1, -1, 4, -4, 6, -6, 3, 2, 9])
    nbl = list(range(1, 1 + nn))
    center = atom(rng.choice(["C", "C", "C", "N", "N", "Si", "O"]), [(q, bs()) for q in nbl], 1)
    center["nh"] = rng.choice([0, 0, 1, 1, 2])
    if rng.random() < 0.4:
        center["iso"] = [0, 0, 0]
    atoms = [center]
    for q in nbl:
        atoms.append(atom("C", [(0, 0)], 1 + q, xyz()))
    # filler atoms so that the array has n atoms
    for q in range(n - 1 - nn):
        atoms.append(atom("C", [(0, 0)], 1 + nn + 1 + q, xyz()))
    atoms = atoms[:1 + nn + 1]  # keep 1 + nn + 1 main atoms
    for a in atoms[1:]:
        a["nb"] = [(0, rng.choice([0, 0, 1, -1]))]
    # a stereo bond partner among the neighbors for the 0D parities
    if rng.random() < 0.5:
        center["sbp"][0] = rng.choice([1, 2, 3, 4, 1, 2])
        center["sbo"][0] = rng.randrange(nn)
        center["sno"][0] = rng.choice([a["orig"] for a in atoms[1:]] + [0, 99])
        if rng.random() < 0.3:
            center["sbp"][1] = rng.choice([1, 2, 3, 4])
            center["sbo"][1] = rng.randrange(nn)
            center["sno"][1] = rng.choice([a["orig"] for a in atoms[1:]] + [0])
        partner = atoms[1 + center["sbo"][0]]
        partner["nb"] = [(0, 0)] + [(q, 0) for q in range(1, 1 + rng.randrange(0, 3)) if q != 1 + center["sbo"][0]]
    if rng.random() < 0.25:
        for a in atoms[1:]:
            a["u0"] = rng.choice([0, 2])
    n = len(atoms)
    for j in range(nr):
        tgt = 0 if rng.random() < 0.85 else rng.randrange(n)
        h = atom("H", [(tgt, rng.choice([0, 0, 1, -1, 4, -4, 6, -6]))], 100 + j, xyz())
        h["adiff"] = rng.choice([0, 0, 1, 2, 3, 0, 1, 4 if rng.random() < 0.03 else 2])
        h["nh"] = 0
        atoms.append(h)
    attached = sum(1 for a in atoms[n:] if a["nb"][0][0] == 0)
    center["nh"] = max(center["nh"], attached)
    return {"n": n, "nr": nr, "pe": rng.choice([0, 1, 1, 3, 7]), "vab": rng.choice([3, 3, 3, 6]), "usez": rng.choice([0, 1, 1]),
            "z": [rng.choice([0, 0, 100, -100, 30]) for _ in range(3)], "atoms": atoms}


def fmt(c):
    t = f"{c['n']} {c['nr']} {c['pe']} {c['vab']} {c['usez']} " + " ".join(map(str, c["z"])) + "\n"
    for a in c["atoms"]:
        t += f"{a['el']} {a['ch']} {a['rad']} {a['nh']} {a['cf']} {len(a['nb'])} {a['orig']} {a['adiff']} {a['u0']} " + " ".join(map(str, a["iso"])) + "\n"
        t += " ".join(map(str, a["sbp"] + a["sbo"] + a["sno"])) + "\n"
        t += " ".join(str(bits(x)) for x in a["xyz"]) + " " + " ".join(f"{q} {b}" for q, b in a["nb"]) + "\n"
    return t


def main():
    rng = random.Random(31)
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
    assert len(out) == len(cases)
    for c in cases:
        for a in c["atoms"]:
            a["xyz"] = [bits(x) for x in a["xyz"]]
    data = [{"in": c, "out": o} for c, o in zip(cases, out)]
    path = ROOT / "tests/data/inchi_halfsb.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
