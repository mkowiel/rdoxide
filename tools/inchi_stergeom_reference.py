"""Generates tests/data/inchi_stergeom.json from the InChI C code: the
geometry helpers of ichister.c (`get_z_coord`, `triple_prod*`,
`Get2DTetrahedralAmbiguity`, `are_*_in_one_plane`, `triple_prod_char`, ...).
The harness `#include`s ichister.c (the helpers are static) and links the
rest of the upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_stergeom_reference.py
"""
import importlib.util
import json
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
/* ichister.c only uses pow(x, 2.0); libm's pow is not always the correctly
   rounded x*x, which is the intended (and portable) value */
#define pow(x, y) ((x) * (x))
#include "ichister.c"
int CreateCompositeNormAtom() { return 0; }
static double rd(void) { unsigned long long u; double d; scanf("%llu", &u); memcpy(&d, &u, 8); return d; }
static void pd(const char *tag, double d) { unsigned long long u; memcpy(&u, &d, 8); printf("%s %llu\n", tag, u); }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    CANON_GLOBALS cg;
    memset(&cg, 0, sizeof(cg));
    for (int c = 0; c < ncases; c++) {
        int kind;
        scanf("%d", &kind);
        printf("N %d\n", c);
        if (kind == 0) { /* vector ops */
            double a[3], b[3], cc[3], r[3], s;
            for (int i = 0; i < 3; i++) a[i] = rd();
            for (int i = 0; i < 3; i++) b[i] = rd();
            for (int i = 0; i < 3; i++) cc[i] = rd();
            pd("len3", len3(a)); pd("len2", len2(a));
            diff3(a, b, r); pd("d0", r[0]); pd("d1", r[1]); pd("d2", r[2]);
            cross_prod3(a, b, r); pd("x0", r[0]); pd("x1", r[1]); pd("x2", r[2]);
            pd("dot", dot_prod3(a, b));
            pd("tp", triple_prod(a, b, cc, &s)); pd("sine", s);
            pd("tp0", triple_prod(a, b, cc, NULL));
        } else if (kind == 1) { /* Get2DTetrahedralAmbiguity */
            int add, fix; double ang; double co[4][3];
            scanf("%d %d", &add, &fix);
            ang = rd();
            for (int i = 0; i < 4; i++) for (int j = 0; j < 3; j++) co[i][j] = rd();
            printf("T %d\n", Get2DTetrahedralAmbiguity(&cg, co, add, fix, ang));
        } else if (kind == 2) { /* sines / planes */
            double co[4][3], central[3], ms, minsine, vmin; int add, amb;
            for (int i = 0; i < 4; i++) for (int j = 0; j < 3; j++) co[i][j] = rd();
            for (int j = 0; j < 3; j++) central[j] = rd();
            vmin = rd(); minsine = rd();
            scanf("%d", &add);
            pd("p1", triple_prod_and_min_abs_sine(co, &ms)); pd("ms1", ms);
            printf("P3 %d\n", are_3_vect_in_one_plane(co, minsine));
            printf("P4 %d\n", are_4at_in_one_plane(co, minsine));
            amb = 0;
            pd("p2", triple_prod_and_min_abs_sine2(co, central, add, &ms, &amb, vmin)); pd("ms2", ms); printf("A %d\n", amb);
            amb = 5;
            pd("p2b", triple_prod_and_min_abs_sine2(co, central, add, &ms, &amb, vmin)); pd("ms2b", ms); printf("A %d\n", amb);
            pd("p2c", triple_prod_and_min_abs_sine2(co, central, add, &ms, NULL, vmin)); pd("ms2c", ms);
            pd("p2d", triple_prod_and_min_abs_sine2(co, central, add, NULL, NULL, vmin));
        } else { /* atoms: get_z_coord and triple_prod_char */
            int n; scanf("%d", &n);
            inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
            for (int i = 0; i < n; i++) {
                int v, u0;
                at[i].x = rd(); at[i].y = rd(); at[i].z = rd();
                scanf("%d %d", &v, &u0);
                at[i].valence = v; at[i].bUsed0DParity = u0;
                for (int j = 0; j < v; j++) { int q, bs; scanf("%d %d", &q, &bs); at[i].neighbor[j] = q; at[i].bond_stereo[j] = bs; }
            }
            int pe; scanf("%d", &pe);
            for (int i = 0; i < n; i++) for (int j = 0; j < at[i].valence; j++) {
                int t; double z = get_z_coord(at, i, j, &t, pe);
                pd("z", z); printf("Z %d\n", t);
            }
            int a1, i1, a2, i2, zd1[3], zd2[3];
            scanf("%d %d %d %d %d %d %d %d %d %d", &a1, &i1, &a2, &i2, zd1, zd1 + 1, zd1 + 2, zd2, zd2 + 1, zd2 + 2);
            S_CHAR z1[3] = {zd1[0], zd1[1], zd1[2]}, z2[3] = {zd2[0], zd2[1], zd2[2]};
            printf("C %d\n", triple_prod_char(at, a1, i1, z1, a2, i2, z2));
            printf("D %d\n", dot_prodchar3(z1, z2));
            free(at);
        }
    }
    return 0;
}
"""


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def rnd(rng):
    r = rng.random()
    if r < 0.15:
        return 0.0
    if r < 0.25:
        return float(rng.randrange(-3, 4))
    if r < 0.35:
        return rng.uniform(-1e-3, 1e-3)
    return rng.uniform(-3, 3)


def vec(rng):
    return [rnd(rng) for _ in range(3)]


def gen(rng):
    k = rng.choice([0, 1, 1, 1, 2, 2, 3])
    if k == 0:
        return {"k": 0, "v": [vec(rng) for _ in range(3)]}
    if k == 1:
        co = []
        for _ in range(4):
            ang = rng.uniform(0, 6.2831853)
            z = rng.choice([0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.5])
            r = rng.uniform(0.5, 1.5)
            import math
            # sometimes duplicate / opposite directions
            co.append([r * math.cos(ang), r * math.sin(ang), z])
        if rng.random() < 0.3:
            co[1][0], co[1][1] = co[0][0] + rng.choice([0, 1e-7, -1e-7, 1e-5]), co[0][1]
        if rng.random() < 0.2:
            co[2][0], co[2][1] = -co[0][0], -co[0][1]
        return {"k": 1, "add": rng.choice([0, 1]), "fix": rng.choice([0, 1]), "ang": rng.choice([0.1, 0.001, 0.087156]), "co": co}
    if k == 2:
        co = [vec(rng) for _ in range(4)]
        if rng.random() < 0.3:
            co[3] = [0.0, 0.0, 0.0]
        return {"k": 2, "co": co, "central": vec(rng), "vmin": rng.choice([0.03, 0.001, 0.125]), "minsine": rng.choice([0.03, 0.1, 0.001]), "add": rng.choice([0, 1])}
    n = rng.randrange(3, 7)
    atoms = []
    for i in range(n):
        v = rng.randrange(1, min(n, 4))
        nbs = rng.sample([q for q in range(n) if q != i], v)
        atoms.append({"xyz": [rng.choice([0.0, rnd(rng)]) for _ in range(3)] if rng.random() < 0.7 else [0.0, 0.0, 0.0],
                      "u0": rng.choice([0, 0, 1]),
                      "nb": [(q, rng.choice([0, 0, 1, 4, 6, -1, -4, -6, 2])) for q in nbs]})
    if rng.random() < 0.5:
        for a in atoms:
            a["xyz"][2] = 0.0
    a1, a2 = rng.randrange(n), rng.randrange(n)
    return {"k": 3, "atoms": atoms, "pe": rng.choice([0, 1, -1]), "a1": a1, "i1": rng.randrange(len(atoms[a1]["nb"])),
            "a2": a2, "i2": rng.randrange(len(atoms[a2]["nb"])), "z1": [rng.randrange(-100, 101) for _ in range(3)], "z2": [rng.randrange(-100, 101) for _ in range(3)]}


def fmt(c):
    f = lambda v: " ".join(str(bits(x)) for x in v)
    if c["k"] == 0:
        return "0\n" + "\n".join(f(v) for v in c["v"]) + "\n"
    if c["k"] == 1:
        return f"1\n{c['add']} {c['fix']}\n{bits(c['ang'])}\n" + "\n".join(f(v) for v in c["co"]) + "\n"
    if c["k"] == 2:
        return "2\n" + "\n".join(f(v) for v in c["co"]) + "\n" + f(c["central"]) + f"\n{bits(c['vmin'])} {bits(c['minsine'])}\n{c['add']}\n"
    t = f"3\n{len(c['atoms'])}\n"
    for a in c["atoms"]:
        t += f(a["xyz"]) + f" {len(a['nb'])} {a['u0']} " + " ".join(f"{q} {b}" for q, b in a["nb"]) + "\n"
    t += f"{c['pe']}\n{c['a1']} {c['i1']} {c['a2']} {c['i2']} " + " ".join(map(str, c["z1"] + c["z2"])) + "\n"
    return t


def main():
    rng = random.Random(11)
    cases = [gen(rng) for _ in range(4000)]
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
        if c["k"] == 3:
            for a in c["atoms"]:
                a["xyz"] = [bits(x) for x in a["xyz"]]
        else:
            for key in ("v", "co"):
                if key in c:
                    c[key] = [[bits(x) for x in v] for v in c[key]]
            for key in ("central",):
                if key in c:
                    c[key] = [bits(x) for x in c[key]]
            for key in ("ang", "vmin", "minsine"):
                if key in c:
                    c[key] = bits(c[key])
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_stergeom.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
