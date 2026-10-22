"""Generates tests/data/inchi_getinchi.json from the InChI C code: `GetINCHI`
and `GetStdINCHI` (option parsing, `ExtractOneStructure`, the whole InChI
generation and the output assembly) on random, often invalid `inchi_Input`
structures with 0D stereo descriptors and option strings. Links the upstream
InChI base library and the libinchi API sources (see
inchi_mol2inchi_reference.py).

    python tools/inchi_getinchi_reference.py
"""
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


mi = load("mol2inchi", "tools/inchi_mol2inchi_reference.py")
sa = load("setatoms", "tools/inchi_setatoms_reference.py")
ma = mi.ma
SRC, API = mi.SRC, mi.API

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inchi_api.h"
static void hexs(const char *s) { if (!s) { printf("-"); return; } int n = (int) strlen(s); for (int i = 0; i < n; i++) printf("%02x", (unsigned char) s[i]); if (!n) printf("="); }
static int unhex(const char *h, char *out) {
    int n = 0;
    out[0] = 0;
    if (h[0] == '=') return 0;
    for (; h[0] && h[1]; h += 2) { unsigned v; sscanf(h, "%2x", &v); out[n++] = (char) v; }
    out[n] = 0;
    return n;
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int std, use_opt, n, ns;
        char hx[512], opt[512];
        printf("#CASE %d\n", cs);
        scanf("%d %d %s %d %d", &std, &use_opt, hx, &n, &ns);
        unhex(hx, opt);
        inchi_Atom *ati = calloc(n ? n : 1, sizeof(inchi_Atom));
        for (int i = 0; i < n; i++) {
            char h[64]; int nb, iso, rad, chg, nh[4];
            double x, y, z;
            scanf("%s %lf %lf %lf %d %d %d %d %d %d %d %d", h, &x, &y, &z, &nb, &iso, &rad, &chg, &nh[0], &nh[1], &nh[2], &nh[3]);
            int m = 0;
            for (const char *q = h + 1; q[0] && q[1]; q += 2) { unsigned v; sscanf(q, "%2x", &v); ati[i].elname[m++] = (char) v; }
            ati[i].x = x; ati[i].y = y; ati[i].z = z;
            ati[i].num_bonds = nb; ati[i].isotopic_mass = iso; ati[i].radical = rad; ati[i].charge = chg;
            for (int k = 0; k < 4; k++) ati[i].num_iso_H[k] = nh[k];
            for (int j = 0; j < nb; j++) { int a, b, s; scanf("%d %d %d", &a, &b, &s); ati[i].neighbor[j] = a; ati[i].bond_type[j] = b; ati[i].bond_stereo[j] = s; }
        }
        inchi_Stereo0D *st = calloc(ns ? ns : 1, sizeof(inchi_Stereo0D));
        for (int i = 0; i < ns; i++) {
            int t, p, c, a[4];
            scanf("%d %d %d %d %d %d %d", &t, &p, &c, &a[0], &a[1], &a[2], &a[3]);
            st[i].type = t; st[i].parity = p; st[i].central_atom = c;
            for (int k = 0; k < 4; k++) st[i].neighbor[k] = a[k];
        }
        inchi_Input inp;
        memset(&inp, 0, sizeof(inp));
        inp.atom = n ? ati : NULL; inp.num_atoms = n;
        inp.stereo0D = ns ? st : NULL; inp.num_stereo0D = ns;
        inp.szOptions = use_opt ? opt : NULL;
        inchi_Output out;
        memset(&out, 0, sizeof(out));
        int r = std ? GetStdINCHI(&inp, &out) : GetINCHI(&inp, &out);
        printf("R %d\n", r);
        printf("I "); hexs(out.szInChI); printf("\n");
        printf("X "); hexs(out.szAuxInfo); printf("\n");
        printf("M "); hexs(out.szMessage); printf("\n");
        printf("L "); hexs(out.szLog); printf("\n");
    }
    return 0;
}
"""

R = random.Random(20261023)
OPTS = mi.OPTS + ["-Polymers", "-FoldCRU", "-W5", "-SUU -SLUUD", "-LooseTSACheck", "-Tabbed", "-NoLabels"]


def clean_case(r):
    n = r.randrange(1, 13)
    pal = ["C"] * 10 + ["N"] * 3 + ["O"] * 3 + ["S", "P", "Cl", "F", "Br", "Na", "Fe", "H", "B", "Si", "I"]
    els = [r.choice(pal) for _ in range(n)]
    edges = {(r.randrange(i), i) for i in range(1, n) if r.random() < 0.92}
    for _ in range(r.randrange(0, 3)):
        if n > 2:
            a, b = r.sample(range(n), 2)
            edges.add((min(a, b), max(a, b)))
    adj = [[] for _ in range(n)]
    d3 = r.random() < 0.4
    for a, b in sorted(edges):
        bt = r.choice([1] * 6 + [2, 2, 3, 4])
        st = r.choice([0] * 6 + [1, 4, 6])
        adj[a].append((b, bt, st))
        adj[b].append((a, bt, -st))
    atoms = []
    for i in range(n):
        atoms.append({"el": els[i], "xyz": [r.uniform(-4, 4), r.uniform(-4, 4), r.uniform(-2, 2) if d3 else 0.0],
                      "iso": r.choice([0] * 8 + [13, 2, 18]), "rad": r.choice([0] * 10 + [1, 2, 3]), "chg": r.choice([0] * 8 + [1, -1, 2]),
                      "nh": [-1, 0, 0, 0] if r.random() < 0.9 else [r.choice([0, 1, 2]), 0, 0, 0], "nb": adj[i][:20]})
    return {"n": n, "dna": 0, "atoms": atoms}


def make_case(r):
    c = clean_case(r) if r.random() < 0.6 else sa.case(r, r.random() < 0.8)
    n = c["n"]
    stereo = []
    for _ in range(r.choice([0, 0, 0, 1, 2, 3])):
        t = r.choice([1, 2, 2, 3, 0, 5])
        p = r.choice([1, 2, 1, 2, 3, 4, 0, 9, 10, 17])
        pick = lambda: r.randrange(-1, n + 1)
        if t == 2 and n > 1:
            cen = r.randrange(n)
            nb = [x[0] for x in c["atoms"][cen]["nb"]][:4]
            nb = nb + [cen] * (4 - len(nb))
            stereo.append((t, p, cen, nb))
        elif t == 1:
            stereo.append((t, p, 32767, [pick() for _ in range(4)]))
        else:
            stereo.append((t, p, r.randrange(n), [pick() for _ in range(4)]))
    x = r.random()
    opt = None if x < 0.3 else ("" if x < 0.45 else " ".join(r.sample(OPTS, r.choice([1, 1, 2, 3]))))
    return c, stereo, opt


def fmt(c, stereo, opt, std):
    t = f"{std} {0 if opt is None else 1} {opt.encode().hex() if opt else '='} {c['n']} {len(stereo)} "
    for a in c["atoms"]:
        t += "x" + (a["el"].encode().hex() or "00")[:60] + " " + " ".join(repr(float(v)) for v in a["xyz"]) + f" {len(a['nb'])} {a['iso']} {a['rad']} {a['chg']} " + " ".join(map(str, a["nh"])) + " "
        t += " ".join(f"{a_} {b} {s}" for a_, b, s in a["nb"]) + " "
    for ty, p, cen, nb in stereo:
        t += f"{ty} {p} {cen} " + " ".join(map(str, nb)) + " "
    return " ".join(t.split())


def main():
    cases = []
    while len(cases) < 3000:
        c, stereo, opt = make_case(R)
        if opt and "-OutputSDF" in opt:
            continue
        cases.append(fmt(c, stereo, opt, R.choice([0, 0, 1])))
    objs = ma.build_lib() + mi.build_api()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-ffp-contract=off", "-DTARGET_API_LIB", f"-I{SRC}", f"-I{API}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        kept, out = [], []
        for c in cases:
            try:
                r1 = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=20)
            except subprocess.TimeoutExpired:
                continue
            if r1.returncode == 0:
                again = [subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=20) for _ in range(2)]
                if any(a.returncode != 0 or a.stdout != r1.stdout for a in again):
                    continue
                kept.append(c)
                out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode(errors="replace").split("\n")[1:]
        print(len(cases) - len(kept), "cases dropped")
        cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_getinchi.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
