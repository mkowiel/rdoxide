"""Generates tests/data/inchi_prt1h.json from the InChI C code: the small
helpers of the polymer output in ichiprt1.c (`OutputINChIPlainError`,
`IsBondAtomNumsLesser`, `CountPseudoElementInFormula`,
`InternallyGetCanoNumsAndComponentNums`, `EditINCHI_HidePolymerZz`). The
harness includes ichiprt1.c itself. Links the upstream InChI base library
(see inchi_markalt_reference.py).

    python tools/inchi_prt1h_reference.py
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


p3d = load("prt3d", "tools/inchi_prt3d_reference.py")
p3c = p3d.p3c
ma = p3c.ma
SRC = ma.SRC

HARNESS = p3c.HARNESS.split("int str_AuxEqu(")[0] + r"""
#define inchi_strbuf_reset(b) (zero_stack(), inchi_strbuf_reset(b))
#include "ichiprt1.c"
static void decode(char *t) {
    char *o = t;
    for (char *p = t; *p; p++) {
        if (p[0] == '~' && p[1] == 'n') { *o++ = '\n'; p++; }
        else *o++ = *p;
    }
    *o = 0;
}
int main(void) {
    int ncases;
    char cmd[32];
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        char t[4096], u[256];
        scanf("%s", cmd);
        printf("#CASE %d\n", c);
        if (!strcmp(cmd, "err")) {
            int kind = rd(); scanf("%s", t); decode(t);
            INCHI_IOSTREAM out; inchi_ios_init(&out, INCHI_IOS_TYPE_STRING, NULL);
            int ret = OutputINChIPlainError(&out, t, kind);
            printf("%d [%s]\n", ret, out.s.pStr ? out.s.pStr : "");
        } else if (!strcmp(cmd, "bondlt")) {
            int b1[2], b2[2];
            b1[0] = rd(); b1[1] = rd(); b2[0] = rd(); b2[1] = rd();
            printf("%d\n", IsBondAtomNumsLesser(b1, b2));
        } else if (!strcmp(cmd, "count")) {
            scanf("%s", t);
            printf("%d\n", CountPseudoElementInFormula("Zz", t));
        } else if (!strcmp(cmd, "hide")) {
            int pzz = rd(), zy = rd();
            scanf("%s", t); decode(t);
            INCHI_IOSTREAM out; inchi_ios_init(&out, INCHI_IOS_TYPE_STRING, NULL);
            inchi_ios_print_nodisplay(&out, "%s", t);
            EditINCHI_HidePolymerZz(&out, pzz, zy);
            printf("%d [%s]\n", out.s.nUsedLength, out.s.pStr ? out.s.pStr : "");
        } else {
            CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
            INCHI_IOS_STRING sb; inchi_strbuf_init(&sb, 64, 16);
            INCHI_OUT_CTL io; memset(&io, 0, sizeof(io));
            int n = rd(), ot = rd();
            INCHI_SORT *a = rd_sorts(n);
            io.TAUT_MODE = rd(); io.bSecondNonTautPass = 0; io.bOmitRepetitions = rd();
            int nat = rd();
            io.bOutType = ot; io.num_components = n; io.pINChISort = a; io.pINChISort2 = a;
            int *cano = calloc(nat + 2, sizeof(int)), *comp = calloc(nat + 2, sizeof(int));
            int ret = InternallyGetCanoNumsAndComponentNums(&cg, &sb, &io, nat, cano, comp);
            printf("%d", ret);
            for (int i = 0; i <= nat; i++) printf(" %d", cano[i]);
            printf(" |");
            for (int i = 0; i <= nat; i++) printf(" %d", comp[i]);
            printf("\n");
        }
    }
    return 0;
}
"""

R = random.Random(1999)
LAYERS = ["/c1-2-3(4)5", "/h2-3H,1H3", "/q+1", "/p-1", "/b3-2+", "/t2-,3+;1?", "/m1", "/s1", "/i1D", "/f", "/r", "/z101-1-3(8-10,2-4)", "/z201-1-3(8,10-2)",
          "/c1-2", "/h1H", "/o(1,2)", "/t1-", "/m0", "/s1/i1D"]
FORMS = ["C2H6OZz2", "Zz", "CZz2", "Zz3", "C4H8Zz2.2Zz3", "Zz2.CH4", "H2Zz", "C3Zz2.Zz4Zz", "OZz2Zz5", "2Zz", "C6H6"]


def hide_text(r):
    s = ("InChI=1B/" if r.random() < 0.93 else "InChI=1S/") + r.choice(FORMS)
    for _ in range(r.randrange(0, 6)):
        s += r.choice(LAYERS)
    if r.random() < 0.5:
        s += "/z" + r.choice(["101-1-3(8-10,2-4)", "201-1-3(8,10-2)", "101-1-2(1-3,4-5)"])
    if r.random() < 0.3:
        s += "~n" + "AuxInfo=1/0/N:1,2,3/rA:3CZz.Zy/rB:;s1/rC:;"
    if r.random() < 0.6:
        s += "~n"
    if r.random() < 0.2:
        s += "/r" + r.choice(FORMS) + r.choice(LAYERS)
    return s


def main():
    cases = []
    for _ in range(7000):
        cmd = R.choice(["err", "bondlt", "count", "hide", "hide", "hide", "cano"])
        if cmd == "err":
            cases.append(f"err {R.choice([1, 2, 3, 0, 5])} {R.choice(['bad_input', 'x', 'Structure~nhas~nproblem', 'A_B'])}")
        elif cmd == "bondlt":
            cases.append("bondlt " + " ".join(str(R.randrange(0, 5)) for _ in range(4)))
        elif cmd == "count":
            f = R.choice(FORMS) + R.choice(["", "/c1-2", "/h1H", ".Zz2", ".2Zz3", ".3C2Zz2", "..Zz"])
            cases.append(f"count {f}")
        elif cmd == "hide":
            cases.append(f"hide {R.choice([0, 1, 2, 3, 5])} {R.choice([0, 0, 1, 2])} {hide_text(R)}")
        else:
            pool = [p3d.pool_item(R) for _ in range(R.randrange(1, 4))]
            n = R.randrange(1, 4)
            a = p3c.sorts(R, pool, n)
            cases.append(f"cano {n} {R.choice([0, 1, 2, 3, 4])} {a} {R.choice([0, 2])} {R.choice([0, 1])} {R.randrange(20, 26)}")
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
    path = ROOT / "tests/data/inchi_prt1h.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
