"""Generates tests/data/inchi_cmp3.json from the InChI C code:
`CompareReversedINChI3` (ichirvr7.c) on random pairs of INChI structures.

    python tools/inchi_cmp3_reference.py
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


r1 = load("rvr1", "tools/inchi_rvr1_reference.py")
r2, ma, SRC = r1.r2, r1.ma, r1.SRC

MAIN = r1.PRE + r"""
static INChI *rd(void) {
    int has;
    scanf("%d", &has);
    if (!has) return NULL;
    INChI *p = calloc(1, sizeof(INChI));
    int n, k, t;
    scanf("%d %d", &p->nErrorCode, &n);
    p->nNumberOfAtoms = n;
    p->nAtom = calloc(n + 1, 1); p->nNum_H = calloc(n + 1, 1);
    for (int i = 0; i < n; i++) { scanf("%d", &t); p->nAtom[i] = t; }
    for (int i = 0; i < n; i++) { scanf("%d", &t); p->nNum_H[i] = t; }
    scanf("%d", &t);
    if (t) { p->nNum_H_fixed = calloc(n + 1, 1); for (int i = 0; i < n; i++) { scanf("%d", &k); p->nNum_H_fixed[i] = k; } }
    char hill[64]; scanf("%63s", hill);
    p->szHillFormula = strcmp(hill, "-") ? strdup(hill) : NULL;
    scanf("%d", &p->lenConnTable);
    p->nConnTable = calloc(p->lenConnTable + 1, sizeof(AT_NUMB));
    for (int i = 0; i < p->lenConnTable; i++) { scanf("%d", &t); p->nConnTable[i] = t; }
    scanf("%d", &p->lenTautomer);
    p->nTautomer = calloc(p->lenTautomer + 8, sizeof(AT_NUMB));
    for (int i = 0; i < p->lenTautomer; i++) { scanf("%d", &t); p->nTautomer[i] = t; }
    scanf("%d", &p->nNumberOfIsotopicAtoms);
    p->IsotopicAtom = calloc(p->nNumberOfIsotopicAtoms + 1, sizeof(p->IsotopicAtom[0]));
    for (int i = 0; i < p->nNumberOfIsotopicAtoms; i++) {
        int a, b, c, d, e; scanf("%d %d %d %d %d", &a, &b, &c, &d, &e);
        p->IsotopicAtom[i].nAtomNumber = a; p->IsotopicAtom[i].nIsoDifference = b;
        p->IsotopicAtom[i].nNum_H = c; p->IsotopicAtom[i].nNum_D = d; p->IsotopicAtom[i].nNum_T = e;
    }
    scanf("%d", &p->nTotalCharge);
    for (int w = 0; w < 2; w++) {
        scanf("%d", &has);
        if (!has) continue;
        INChI_Stereo *s = calloc(1, sizeof(INChI_Stereo));
        scanf("%d", &s->nNumberOfStereoCenters);
        s->nNumber = calloc(s->nNumberOfStereoCenters + 1, sizeof(AT_NUMB)); s->t_parity = calloc(s->nNumberOfStereoCenters + 1, 1);
        for (int i = 0; i < s->nNumberOfStereoCenters; i++) { scanf("%d %d", &t, &k); s->nNumber[i] = t; s->t_parity[i] = k; }
        scanf("%d", &s->nCompInv2Abs);
        scanf("%d", &s->nNumberOfStereoBonds);
        s->nBondAtom1 = calloc(s->nNumberOfStereoBonds + 1, sizeof(AT_NUMB)); s->nBondAtom2 = calloc(s->nNumberOfStereoBonds + 1, sizeof(AT_NUMB)); s->b_parity = calloc(s->nNumberOfStereoBonds + 1, 1);
        for (int i = 0; i < s->nNumberOfStereoBonds; i++) { int a, b, c; scanf("%d %d %d", &a, &b, &c); s->nBondAtom1[i] = a; s->nBondAtom2[i] = b; s->b_parity[i] = c; }
        if (w == 0) p->Stereo = s; else p->StereoIsotopic = s;
    }
    return p;
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        INChI *a = rd(), *b = rd();
        int err = 0;
        INCHI_MODE r = CompareReversedINChI3(a, b, NULL, NULL, &err);
        printf("%llu %d\n", (unsigned long long)r, err);
    }
    return 0;
}
"""
R = random.Random(2468)
HILLS = ["C2H6", "C2H6O", "C3H8", "CH4", "C2H4O2", "C2H5N", "C2H6O2", "C3H6O", "H2O"]


def stereo():
    if R.random() < 0.5:
        return None
    nsc = R.randrange(0, 4)
    sc = sorted(R.sample(range(1, 9), nsc))
    nsb = R.randrange(0, 3)
    sb = sorted(set((R.randrange(1, 7), R.randrange(2, 9)) for _ in range(nsb)))
    return {"sc": [(a, R.randrange(1, 5)) for a in sc], "inv": R.randrange(0, 3),
            "sb": [(a, b, R.randrange(1, 5)) for a, b in sb]}


def base():
    n = R.randrange(0, 6)
    ng = R.choice([0, 0, 1, 2])
    taut = [ng]
    for _ in range(ng):
        k = R.randrange(1, 4)
        taut += [k + 2, R.randrange(0, 3), R.randrange(0, 2)] + sorted(R.sample(range(1, 8), k))
    if ng == 0 and R.random() < 0.3:
        taut = [0]
    elif ng == 0:
        taut = []
    return {"err": 0, "n": n, "atoms": [R.choice([6, 7, 8]) for _ in range(n)],
            "nh": [R.randrange(0, 4) for _ in range(n)],
            "fh": [R.choice([0, 0, 0, 1, 2]) for _ in range(n)] if R.random() < 0.4 else None,
            "hill": R.choice(HILLS), "conn": [R.randrange(1, 9) for _ in range(R.randrange(0, 8))],
            "taut": taut,
            "iso": [(R.randrange(1, 9), R.randrange(-1, 3), R.randrange(0, 2), R.randrange(0, 2), R.randrange(0, 2)) for _ in range(R.choice([0, 0, 1, 2]))],
            "chg": R.choice([0, 0, 0, 1, -1]), "st": [stereo(), stereo()]}


def mutate(b):
    import copy
    c = copy.deepcopy(b)
    for _ in range(R.randrange(0, 3)):
        w = R.randrange(14)
        if w == 0: c["err"] = R.choice([0, 0, 1])
        elif w == 1 and c["n"]: c["atoms"][R.randrange(c["n"])] = R.choice([6, 7, 8, 16])
        elif w == 2 and c["n"]: c["nh"][R.randrange(c["n"])] += R.choice([-1, 1])
        elif w == 3: c["fh"] = None if c["fh"] else [R.choice([0, 1]) for _ in range(c["n"])]
        elif w == 4 and c["fh"]: c["fh"][R.randrange(c["n"])] += R.choice([-1, 1, 2])
        elif w == 5: c["hill"] = R.choice(HILLS)
        elif w == 6 and c["conn"]: c["conn"][R.randrange(len(c["conn"]))] += 1
        elif w == 7: c["conn"] = c["conn"] + [R.randrange(1, 9)]
        elif w == 8: c["taut"] = base()["taut"]
        elif w == 9 and len(c["taut"]) > 4: c["taut"][-1] += 1
        elif w == 10: c["iso"] = base()["iso"]
        elif w == 11: c["chg"] += 1
        elif w == 12: c["st"][R.randrange(2)] = stereo()
        elif w == 13 and c["st"][0]:
            s = c["st"][0]
            if s["sc"]: s["sc"][0] = (s["sc"][0][0], R.randrange(1, 5))
            s["inv"] = R.randrange(0, 3)
    return c


def fmt(c):
    if c is None:
        return "0"
    t = [1, c["err"], c["n"], *c["atoms"], *c["nh"]]
    t += [1, *c["fh"]] if c["fh"] is not None else [0]
    t += [c["hill"] if c["hill"] else "-", len(c["conn"]), *c["conn"], len(c["taut"]), *c["taut"], len(c["iso"])]
    for x in c["iso"]:
        t += list(x)
    t.append(c["chg"])
    for s in c["st"]:
        if s is None:
            t.append(0)
            continue
        t += [1, len(s["sc"])]
        for a, p in s["sc"]:
            t += [a, p]
        t += [s["inv"], len(s["sb"])]
        for a, b, p in s["sb"]:
            t += [a, b, p]
    return " ".join(map(str, t))


def main():
    cases = []
    for _ in range(6000):
        b = base()
        a = mutate(b) if R.random() < 0.9 else base()
        pair = [a, b]
        if R.random() < 0.03:
            pair[R.randrange(2)] = None
        cases.append(pair)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(r2.HARNESS.split("int main(void)")[0] + MAIN)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        inp = f"{len(cases)}\n" + "\n".join(fmt(a) + "\n" + fmt(b) for a, b in cases) + "\n"
        out = subprocess.run([str(exe)], input=inp.encode(), capture_output=True, check=True).stdout.decode().split("\n")
    data = [{"in": [fmt(a), fmt(b)], "out": o} for (a, b), o in zip(cases, out)]
    (ROOT / "tests/data/inchi_cmp3.json").write_text(json.dumps(data, separators=(",", ":")))
    from collections import Counter
    print(len(data), len(Counter(d["out"] for d in data)), "distinct results")


if __name__ == "__main__":
    main()
