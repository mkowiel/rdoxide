"""Generates tests/data/inchi_mol2inchi.json from the InChI C code: the whole
Molfile text -> InChI path of `MakeINCHIFromMolfileText` (option parsing,
`GetOneStructure`, `ProcessOneStructureEx`, the output assembly of
`produce_generation_output`) on random and fuzzed V2000/V3000 molfiles and
option strings. Links the upstream InChI base library (see
inchi_markalt_reference.py) and the API sources of libinchi.

    python tools/inchi_mol2inchi_reference.py
"""
import importlib.util
import json
import random
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


m2a = load("mol2atom", "tools/inchi_mol2atom_reference.py")
m2, m3 = m2a.m2, m2a.m3
ma = m2a.ma
SRC = ma.SRC
API = ROOT / "upstream/external/InChI/INCHI_API/libinchi/src"
APIOBJ = ma.CACHE / "api"
API_FILES = ["inchi_dll", "inchi_dll_b", "inchi_dll_a", "inchi_dll_a2", "ichilnct"]

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inchi_api.h"
static void hexs(const char *s) { if (!s) { printf("-"); return; } int n = (int) strlen(s); for (int i = 0; i < n; i++) printf("%02x", (unsigned char) s[i]); if (!n) printf("="); }
static char opt[1 << 14], mol[1 << 18], hx[1 << 19];
static int unhex(const char *h, char *out) {
    int n = 0;
    out[0] = 0;
    if (h[0] == '=') return 0;
    for (; h[0] && h[1]; h += 2) { unsigned v; sscanf(h, "%2x", &v); out[n++] = (char) v; }
    out[n] = 0;
    return n;
}
int MakeINCHIFromMolfileText(const char *moltext, char *szOptions, inchi_Output *result);
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int use_opt;
        char o[1 << 14];
        printf("#CASE %d\n", cs);
        scanf("%d %s", &use_opt, hx);
        unhex(hx, o);
        scanf("%s", hx);
        unhex(hx, mol);
        inchi_Output out;
        memset(&out, 0, sizeof(out));
        int r = MakeINCHIFromMolfileText(mol, use_opt ? o : NULL, &out);
        printf("R %d\n", r);
        printf("I "); hexs(out.szInChI); printf("\n");
        printf("X "); hexs(out.szAuxInfo); printf("\n");
        printf("M "); hexs(out.szMessage); printf("\n");
        printf("L "); hexs(out.szLog); printf("\n");
    }
    return 0;
}
"""

R = random.Random(20261022)
OPTS = ["-SNon", "-FixedH", "-RecMet", "-NoWarnings", "-SUU", "-AuxNone", "-DoNotAddH", "-SAbs", "-SRel", "-SRac", "-SUCF",
        "-KET", "-15T", "-SLUUD", "-NPZz", "-SaveOpt", "-WarnOnEmptyStructure", "-LargeMolecules", "-Foo", "-ChiralFlagON",
        "-ChiralFlagOFF", "-InChI2InChI", "-SNon -SRel", "-DoNotAddH -FixedH", "-NEWPSOFF", "-SPXYZ", "-DISABLEFIXEDH", "-Key", "-XHash1", "-ASC", "-NoWarnings -SUU"]


def pick_opts(r):
    x = r.random()
    if x < 0.15:
        return None
    if x < 0.4:
        return ""
    return " ".join(r.sample(OPTS, r.choice([1, 1, 2, 3])))


def mk_graph(r):
    """A chemically plausible random structure: tree + ring closures, wedges, charges, isotopes, fragments"""
    n = r.randrange(1, 15)
    pal = ["C"] * 10 + ["N"] * 3 + ["O"] * 3 + ["S", "P", "Cl", "F", "Br", "Na", "K", "Fe", "H", "B", "Si", "I", "Se"]
    syms = [r.choice(pal) for _ in range(n)]
    bonds = []
    for i in range(1, n):
        if r.random() < 0.9:
            bonds.append([r.randrange(0, i), i, r.choice([1] * 6 + [2, 2, 3])])
    for _ in range(r.randrange(0, 3)):
        if n > 2:
            a, b = r.sample(range(n), 2)
            bonds.append([a, b, r.choice([1, 1, 2, 4])])
    chg = {}
    for _ in range(r.randrange(0, 3)):
        chg[r.randrange(0, n)] = r.choice([-1, 1, 1, -2, 2])
    d3 = r.random() < 0.3
    L = ["mol", "  fuzz", ""]
    flag = r.choice([0, 1])
    props = []
    for a, c in chg.items():
        props.append(f"M  CHG  1 {a + 1:3d} {c:3d}")
    for _ in range(r.randrange(0, 2)):
        props.append(f"M  ISO  1 {r.randrange(1, n + 1):3d} {r.choice([13, 2, 18, 3, 35, 15]):3d}")
    if r.random() < 0.1:
        props.append(f"M  RAD  1 {r.randrange(1, n + 1):3d} {r.randrange(1, 4):3d}")
    L.append(f"{n:3d}{len(bonds):3d}  0  0{flag:3d}  0  0  0  0  0{len(props) + 1:3d} V2000")
    for s in syms:
        z = r.uniform(-2, 2) if d3 else 0.0
        L.append(f"{r.uniform(-5, 5):10.4f}{r.uniform(-5, 5):10.4f}{z:10.4f} {s:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    for a, b, o in bonds:
        st = r.choice([0] * 7 + [1, 6, 4, 3])
        L.append(f"{a + 1:3d}{b + 1:3d}{o:3d}{st:3d}  0  0")
    L += props + ["M  END"]
    return "\n".join(L) + "\n"


def mk_input(r):
    x = r.random()
    if x < 0.55:
        return mk_graph(r)
    if x < 0.75:
        return m2a.mk_clean(r) + ("$$$$\n" if r.random() < 0.5 else "")
    if x < 0.9:
        return m2.mk_mol2(r)
    return m3.mk_mol3(r)


POLY = re.compile(r"STY|SGROUP|SRU|\bMON\b|\bCOP\b")


def build_api():
    APIOBJ.mkdir(parents=True, exist_ok=True)
    jobs = []
    for n in API_FILES:
        f, o = API / f"{n}.c", APIOBJ / f"{n}.o"
        if not o.exists() or o.stat().st_mtime < f.stat().st_mtime:
            jobs.append(subprocess.Popen(["cc", "-O1", "-w", "-DTARGET_API_LIB", f"-I{SRC}", f"-I{API}", "-c", str(f), "-o", str(o)]))
    for j in jobs:
        assert j.wait() == 0
    return [str(APIOBJ / f"{n}.o") for n in API_FILES]


def main():
    cases = []
    while len(cases) < 3000:
        text = mk_input(R)
        if POLY.search(text):
            continue
        o = pick_opts(R)
        cases.append(f"{0 if o is None else 1} {o.encode().hex() if o else '='} {text.encode().hex()}")
    objs = ma.build_lib() + build_api()
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
    path = ROOT / "tests/data/inchi_mol2inchi.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
