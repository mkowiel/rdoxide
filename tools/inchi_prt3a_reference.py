"""Generates tests/data/inchi_prt3a.json from the InChI C code: the first
whole-structure substring builders of ichiprt3.c (`str_HillFormula`,
`str_HillFormula2`, `str_Connections`, `str_H_atoms`, `str_Charge2`,
`str_FixedH_atoms`) on random multi-component structures. Links the upstream
InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_prt3a_reference.py
"""
import copy
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


m2 = load("make2", "tools/inchi_make2_reference.py")
p2 = load("prt2", "tools/inchi_prt2_reference.py")
ma = m2.ma
SRC = ma.SRC

HARNESS = m2.HARNESS.split("int main(void)")[0] + r"""
#include "ichi_io.h"
#include "ichimain.h"
int str_HillFormula( INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int num_components, int bUseMulipliers );
int str_HillFormula2( INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int num_components, int bUseMulipliers );
int str_Connections( CANON_GLOBALS *pCG, INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int ATOM_MODE, int num_components, int bUseMulipliers );
int str_H_atoms( INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int ATOM_MODE, int TAUT_MODE, int num_components, int bUseMulipliers );
int str_Charge2( INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int num_components, int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers );
int str_FixedH_atoms( INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int ATOM_MODE, int num_components, int bUseMulipliers );
static void zero_stack(void) { volatile char b[32768]; memset((void*)b, 0, sizeof(b)); }
static INCHI_SORT *rd_sorts(int n) {
    INCHI_SORT *s = calloc(n + 1, sizeof(INCHI_SORT));
    for (int i = 0; i < n; i++) { s[i].pINChI[0] = rd_inchi(); s[i].pINChI[1] = rd_inchi(); s[i].ord_number = i; }
    return s;
}
int main(void) {
    int ncases;
    char cmd[32];
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        INCHI_IOS_STRING sb; inchi_strbuf_init(&sb, 16, 16);
        int ovf = 0, ret = 0;
        scanf("%s", cmd);
        int n = rd(), ot = rd(), mult = rd();
        INCHI_SORT *a = rd_sorts(n);
        INCHI_SORT *b = rd_sorts(n);
        int am = rd(), tm = rd(), second = rd(), omit = rd();
        CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
        printf("#CASE %d\n", c);
        zero_stack();
        if (!strcmp(cmd, "hill")) ret = str_HillFormula(a, &sb, &ovf, ot, n, mult);
        else if (!strcmp(cmd, "hill2")) ret = str_HillFormula2(a, b, &sb, &ovf, ot, n, mult);
        else if (!strcmp(cmd, "conn")) ret = str_Connections(&cg, a, &sb, &ovf, ot, am, n, mult);
        else if (!strcmp(cmd, "hatoms")) ret = str_H_atoms(a, &sb, &ovf, ot, am, tm, n, mult);
        else if (!strcmp(cmd, "charge")) ret = str_Charge2(a, b, &sb, &ovf, ot, n, second, omit, mult);
        else if (!strcmp(cmd, "fixedh")) ret = str_FixedH_atoms(a, &sb, &ovf, ot, am, n, mult);
        printf("%d %d [%s]\n", ret, ovf, sb.pStr ? sb.pStr : "");
    }
    return 0;
}
"""

R = random.Random(8080)
MODES = [0, 1, 8, 9, 2, 10, 11, 16, 17, 24, 25]
FORMS = ["C2H6O", "CH4", "C6H6", "ClH", "H2O", "C3H8O2", "C2H6"]


def taut(r, n):
    if r.random() < 0.55:
        return []
    groups = []
    for _ in range(r.choice([1, 1, 2])):
        eps = sorted(r.sample(range(1, n + 1), r.randrange(1, min(n, 3) + 1)))
        groups.append([len(eps) + 2, r.choice([0, 1, 2]), r.choice([0, 0, 1])] + eps)
    return [len(groups)] + [x for g in groups for x in g]


def base(r):
    c = m2.comp(r)
    n = c["n"]
    c["f"] = r.choice(FORMS)
    c["numh"] = [r.choice([0, 0, 1, 2, 3]) for _ in range(n)]
    c["fixed"] = [r.choice([0, 0, 0, 1, 2]) for _ in range(n)]
    c["charge"] = r.choice([0, 0, 0, 1, -1, 2])
    c["st"] = None
    c["sti"] = None
    c["iso"] = []
    c["isotg"] = []
    c["del"] = r.choice([0, 0, 0, 0, 1])
    _, ct = p2.gen_ct(r) if False else gen_ct_n(r, n)
    c["ct"] = ct
    c["taut"] = taut(r, n)
    return c


def gen_ct_n(r, n):
    edges = set()
    for v in range(1, n):
        if r.random() < 0.9:
            edges.add((r.randrange(0, v), v))
    nb = {v: sorted(a for a, b in edges if b == v) for v in range(n)}
    ct = []
    for v in range(n):
        ct.append(v + 1)
        ct += [a + 1 for a in nb[v]]
    return n, ct


def slot(r, c):
    """(non-taut, taut) INChI pair for one component"""
    k = r.choice(["N", "T", "B", "TN", "none", "N", "B"])
    n_ = copy.deepcopy(c)
    n_["taut"] = []
    t_ = copy.deepcopy(c)
    if not t_["taut"]:
        t_["taut"] = taut(r, 3) or [1, 3, 1, 0, 1, 2]
        t_["taut"] = [x for x in t_["taut"]]
        # endpoints must stay within the atom count
        t_["taut"] = [min(x, c["n"]) if i > 3 else x for i, x in enumerate(t_["taut"])]
    tn = copy.deepcopy(c)
    tn["taut"] = []
    if k == "N":
        return n_, None
    if k == "T":
        return None, t_
    if k == "B":
        return n_, t_
    if k == "TN":
        return None, tn
    return None, None


def sorts(r, pool, n):
    out = []
    for _ in range(n):
        c = copy.deepcopy(r.choice(pool))
        if r.random() < 0.2:
            c = m2.mutate(r, c)
            c["st"] = c["sti"] = None
            c["iso"], c["isotg"] = [], []
            n = c["n"]
            c["ct"] = gen_ct_n(r, n)[1]
            c["taut"] = taut(r, n)
            c["numh"] = [r.choice([0, 0, 1, 2, 3]) for _ in range(n)]
            c["fixed"] = [r.choice([0, 0, 0, 1, 2]) for _ in range(n)]
        out.append(slot(r, c))
    return out


def tok_sorts(ss):
    return " ".join((m2.tok(a) if a else "0") + " " + (m2.tok(b) if b else "0") for a, b in ss)


def main():
    cases = []
    for _ in range(6000):
        cmd = r_cmd = R.choice(["hill", "hill2", "conn", "hatoms", "charge", "fixedh"])
        pool = [base(R) for _ in range(R.randrange(1, 4))]
        n = R.randrange(1, 5)
        a = sorts(R, pool, n)
        b = sorts(R, pool, n)
        ot = R.choice([0, 1, 2, 3, 4])
        mult = R.choice([0, 1, 1, 1])
        am = R.choice(MODES)
        tm = R.choice([0, 2])
        second = R.choice([0, 1])
        omit = R.choice([0, 1, 1])
        cases.append(f"{cmd} {n} {ot} {mult} {tok_sorts(a)} {tok_sorts(b)} {am} {tm} {second} {omit}")
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=(f"{len(cases)}\n" + "\n".join(cases) + "\n").encode(), capture_output=True)
        if r.returncode == 0:
            out = r.stdout.decode(errors="replace").split("\n")
        else:
            kept, out = [], []
            for c in cases:
                r1 = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True)
                if r1.returncode == 0:
                    kept.append(c)
                    out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode(errors="replace").split("\n")[1:]
            print(len(cases) - len(kept), "cases dropped")
            cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases), (len(starts), len(cases))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    data = [x for x in data if all(ch.isprintable() and ord(ch) < 127 for l in x["out"] for ch in l)]
    path = ROOT / "tests/data/inchi_prt3a.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
