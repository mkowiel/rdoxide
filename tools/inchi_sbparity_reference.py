"""Generates tests/data/inchi_sbparity.json from the InChI C code:
`set_stereo_bonds_parity` (with `save_a_stereo_bond`, `FixSb0DParities`,
`half_stereo_bond_action`, ...) of ichister.c on random small molecules.
The harness `#include`s ichister.c.

    python tools/inchi_sbparity_reference.py
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
        int n, nr, pe, vab, minring; long mode;
        scanf("%d %d %d %d %d %ld", &n, &nr, &pe, &vab, &minring, &mode);
        int tot = n + nr;
        inp_ATOM *at = calloc(tot + 1, sizeof(inp_ATOM));
        for (int i = 0; i < tot; i++) {
            char el[8]; int ch, rad, nh, cf, v, orig, adiff, u0, cv, ep, cp, rs, iso[3], sbp[3], sbo[3], sno[3];
            scanf("%s %d %d %d %d %d %d %d %d %d %d %d %d", el, &ch, &rad, &nh, &cf, &v, &orig, &adiff, &u0, &cv, &ep, &cp, &rs);
            scanf("%d %d %d %d %d %d %d %d %d %d %d %d", iso, iso + 1, iso + 2, sbp, sbp + 1, sbp + 2, sbo, sbo + 1, sbo + 2, sno, sno + 1, sno + 2);
            strcpy(at[i].elname, el);
            at[i].el_number = (U_CHAR) get_periodic_table_number(el);
            at[i].charge = ch; at[i].radical = rad; at[i].num_H = nh; at[i].cFlags = cf; at[i].valence = v;
            at[i].orig_at_number = orig; at[i].iso_atw_diff = adiff; at[i].bUsed0DParity = u0;
            at[i].chem_bonds_valence = cv; at[i].endpoint = ep; at[i].c_point = cp; at[i].nRingSystem = rs;
            at[i].x = rd(); at[i].y = rd(); at[i].z = rd();
            for (int j = 0; j < 3; j++) { at[i].num_iso_H[j] = iso[j]; at[i].sb_parity[j] = sbp[j]; at[i].sb_ord[j] = sbo[j]; at[i].sn_orig_at_num[j] = sno[j]; }
            for (int j = 0; j < v; j++) { int q, bt, bs; scanf("%d %d %d", &q, &bt, &bs); at[i].neighbor[j] = q; at[i].bond_type[j] = bt; at[i].bond_stereo[j] = bs; }
        }
        sp_ATOM *out = calloc(tot + 1, sizeof(sp_ATOM));
        QUEUE *q = QueueCreate(tot + 1, sizeof(AT_RANK));
        AT_RANK *lev = calloc(tot + 1, sizeof(AT_RANK));
        S_CHAR *src = calloc(tot + 1, sizeof(S_CHAR));
        printf("N %d\n", c);
        for (int a1 = 0; a1 < n; a1++) {
            int r = set_stereo_bonds_parity(out, at, a1, at + n, nr, mode, q, lev, src, minring, pe, vab);
            printf("R %d\n", r);
        }
        for (int i = 0; i < n; i++) {
            sp_ATOM *o = out + i;
            printf("A");
            for (int j = 0; j < 3; j++) printf(" %d %d %d %d %d %d %d %d", o->stereo_bond_neighbor[j], o->stereo_bond_neighbor2[j], o->stereo_bond_ord[j], o->stereo_bond_ord2[j], o->stereo_bond_z_prod[j], o->stereo_bond_z_prod2[j], o->stereo_bond_parity[j], o->stereo_bond_parity2[j]);
            printf(" %d %d %d %d %d %d %d | %d %d\n", o->parity, o->parity2, o->bAmbiguousStereo, o->z_dir[0], o->z_dir[1], o->z_dir[2], 0, at[i].bAmbiguousStereo, at[i].bUsed0DParity);
        }
        QueueDelete(q);
        free(at); free(out); free(lev); free(src);
    }
    return 0;
}
"""


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def build(rng):
    """random small molecule: backbone with double bonds, substituents, maybe a ring"""
    bl = rng.randrange(2, 9)
    atoms = [{"el": rng.choice(["C"] * 7 + ["N", "Si"])} for _ in range(bl)]
    bonds = {}

    def add(a, b, t):
        bonds[(min(a, b), max(a, b))] = t
    for i in range(bl - 1):
        add(i, i + 1, 2 if rng.random() < 0.4 and not (i > 0 and bonds.get((i - 1, i)) == 2 and rng.random() < 0.6) else 1)
    # substituents
    for i in range(bl):
        for _ in range(rng.randrange(0, 3)):
            atoms.append({"el": rng.choice(["C", "C", "C", "N", "O", "F", "Cl", "C"])})
            add(i, len(atoms) - 1, 1)
    if rng.random() < 0.3 and bl >= 3:
        a, b = sorted(rng.sample(range(bl), 2))
        if (a, b) not in bonds:
            add(a, b, rng.choice([1, 1, 2, 4]))
    if rng.random() < 0.15:
        for k in list(bonds):
            if bonds[k] == 2 and rng.random() < 0.5:
                bonds[k] = rng.choice([4, 8, 9])
    return atoms, bonds


def ring_systems(n, bonds):
    """2-edge-connected components (bridges removed); each atom gets an id"""
    adj = [[] for _ in range(n)]
    for (a, b) in bonds:
        adj[a].append(b)
        adj[b].append(a)
    # find bridges via DFS lowpoint
    disc, low, bridges, t = [-1] * n, [0] * n, set(), [0]

    def dfs(u, p):
        disc[u] = low[u] = t[0]
        t[0] += 1
        for v in adj[u]:
            if v == p:
                continue
            if disc[v] < 0:
                dfs(v, u)
                low[u] = min(low[u], low[v])
                if low[v] > disc[u]:
                    bridges.add((min(u, v), max(u, v)))
            else:
                low[u] = min(low[u], disc[v])
    for s in range(n):
        if disc[s] < 0:
            dfs(s, -1)
    comp = [-1] * n
    cid = 0
    for s in range(n):
        if comp[s] >= 0:
            continue
        comp[s] = cid
        st = [s]
        while st:
            u = st.pop()
            for v in adj[u]:
                if comp[v] < 0 and (min(u, v), max(u, v)) not in bridges:
                    comp[v] = cid
                    st.append(v)
        cid += 1
    return [c + 1 for c in comp]


def gen(rng):
    atoms, bonds = build(rng)
    n = len(atoms)
    if n > 12:
        return gen(rng)
    nb = [[] for _ in range(n)]
    for (a, b), t in bonds.items():
        nb[a].append((b, t))
        nb[b].append((a, t))
    if any(len(x) > 3 for x in nb):
        return gen(rng)
    rs = ring_systems(n, bonds)
    planar = rng.random() < 0.6
    out = []
    for i, a in enumerate(atoms):
        order = sum({1: 1, 2: 2, 4: 1, 8: 1, 9: 1}[t] for _, t in nb[i])
        tv = {"C": 4, "N": 3, "Si": 4, "O": 2, "F": 1, "Cl": 1}[a["el"]]
        if rng.random() < 0.1:
            tv += rng.choice([-1, 1])
        nh = max(0, tv - order)
        cv = order + (1 if any(t in (4, 8, 9) for _, t in nb[i]) and rng.random() < 0.5 else 0)
        xyz = [math.cos(i * 1.1 + rng.uniform(-0.4, 0.4)) * (1 + i * 0.2) + 0.1 * i, math.sin(i * 0.9 + rng.uniform(-0.4, 0.4)) * (1 + 0.1 * i), 0.0 if planar else rng.choice([0.0, 0.7, -0.7])]
        if rng.random() < 0.03:
            xyz = [0.0, 0.0, 0.0]
        out.append({"el": a["el"], "ch": 0, "rad": 0, "nh": nh, "cf": rng.choice([0] * 5 + [1]), "orig": i + 1, "adiff": 0,
                    "u0": 0, "cv": cv, "ep": rng.choice([0] * 20 + [1]), "cp": rng.choice([0] * 30 + [1]), "rs": rs[i],
                    "iso": [rng.choice([0, 0, 0, 0, 1]) for _ in range(3)], "sbp": [0] * 3, "sbo": [0] * 3, "sno": [0] * 3,
                    "xyz": xyz,
                    "nb": [(q, t | (rng.choice([0] * 15 + [0x10, 0x20]) if t < 8 else 0), rng.choice([0] * 5 + [1, -1, 6, -6, 4, 3])) for q, t in nb[i]]})
    # symmetric bond_stereo is not needed; 0D parities
    if rng.random() < 0.3:
        for i, a in enumerate(out):
            dbl = [j for j, (q, t, _) in enumerate(a["nb"]) if (t & 0xF) == 2]
            if dbl and rng.random() < 0.8:
                a["u0"] = rng.choice([2, 2, 2, 0])
                a["sbp"][0] = rng.choice([1, 2, 1, 2, 3, 4])
                a["sbo"][0] = rng.choice(dbl)
                a["sno"][0] = rng.choice([out[q]["orig"] for q, _, _ in a["nb"]] + [0])
    # removed explicit H
    nr = rng.choice([0, 0, 0, 1, 2])
    for j in range(nr):
        tgt = rng.randrange(n)
        h = {"el": "H", "ch": 0, "rad": 0, "nh": 0, "cf": 0, "orig": 100 + j, "adiff": rng.choice([0, 0, 1, 2]), "u0": 0, "cv": 1, "ep": 0, "cp": 0, "rs": 0,
             "iso": [0] * 3, "sbp": [0] * 3, "sbo": [0] * 3, "sno": [0] * 3,
             "xyz": [rng.uniform(-1, 1), rng.uniform(-1, 1), rng.choice([0.0, 0.5])], "nb": [(tgt, 1, rng.choice([0, 0, 1, -1, 6, -6]))]}
        out.append(h)
        out[tgt]["nh"] = max(out[tgt]["nh"], 1) + (0 if out[tgt]["nh"] else 0)
    for i in range(n):
        att = sum(1 for h in out[n:] if h["nb"][0][0] == i)
        out[i]["nh"] = max(out[i]["nh"], att)
    return {"n": n, "nr": nr, "pe": rng.choice([0, 1, 1, 3]), "vab": rng.choice([3, 3, 3, 6]), "minring": rng.choice([8, 8, 0, 5]),
            "mode": rng.choice([0, 0, 0x100]), "atoms": out}


def fmt(c):
    t = f"{c['n']} {c['nr']} {c['pe']} {c['vab']} {c['minring']} {c['mode']}\n"
    for a in c["atoms"]:
        t += f"{a['el']} {a['ch']} {a['rad']} {a['nh']} {a['cf']} {len(a['nb'])} {a['orig']} {a['adiff']} {a['u0']} {a['cv']} {a['ep']} {a['cp']} {a['rs']}\n"
        t += " ".join(map(str, a["iso"] + a["sbp"] + a["sbo"] + a["sno"])) + "\n"
        t += " ".join(str(bits(x)) for x in a["xyz"]) + " " + " ".join(f"{q} {b} {s}" for q, b, s in a["nb"]) + "\n"
    return t


def main():
    rng = random.Random(41)
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
    for c in cases:
        for a in c["atoms"]:
            a["xyz"] = [bits(x) for x in a["xyz"]]
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_sbparity.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
