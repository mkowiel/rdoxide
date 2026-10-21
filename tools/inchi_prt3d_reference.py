"""Generates tests/data/inchi_prt3d.json from the InChI C code: the remaining
AuxInfo builders of ichiprt3.c (`str_AuxNumb`, `str_AuxInvSp3Numb`,
`str_AuxIsoNumb`, `str_AuxInvIsoSp3Numb`, `str_AuxChargeRadVal`,
`bin_AuxTautTrans` + `str_AuxTautTrans`) on random multi-component
structures. Links the upstream InChI base library (see
inchi_markalt_reference.py).

    python tools/inchi_prt3d_reference.py
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


p3c = load("prt3c", "tools/inchi_prt3c_reference.py")
p2 = p3c.p2
ma = p3c.ma
SRC = p3c.SRC

HARNESS = p3c.HARNESS.split("int str_AuxEqu(")[0] + r"""
#include <stdio.h>
int str_AuxNumb( CANON_GLOBALS *pCG, INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int second, int omit );
int str_AuxInvSp3Numb( CANON_GLOBALS *pCG, INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int second, int omit );
int str_AuxIsoNumb( CANON_GLOBALS *pCG, INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int second, int omit );
int str_AuxInvIsoSp3Numb( CANON_GLOBALS *pCG, INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int second, int omit );
int str_AuxChargeRadVal( INCHI_SORT *a, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int mult );
int bin_AuxTautTrans( INCHI_SORT *a, INCHI_SORT *b, AT_NUMB **tn, AT_NUMB **ts, int ot, int n );
int str_AuxTautTrans( CANON_GLOBALS *pCG, AT_NUMB *tn, AT_NUMB *ts, INCHI_IOS_STRING *sb, int *ovf, int TM, int n );
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
        int tm = rd(), second = rd(), omit = rd();
        for (int i = 0; i < n; i++) a[i].ord_number = rd();
        for (int i = 0; i < n; i++) b[i].ord_number = rd();
        CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
        printf("#CASE %d\n", c);
        zero_stack();
        if (!strcmp(cmd, "numb")) ret = str_AuxNumb(&cg, a, b, &sb, &ovf, ot, tm, n, second, omit);
        else if (!strcmp(cmd, "invnumb")) ret = str_AuxInvSp3Numb(&cg, a, b, &sb, &ovf, ot, tm, n, second, omit);
        else if (!strcmp(cmd, "isonumb")) ret = str_AuxIsoNumb(&cg, a, b, &sb, &ovf, ot, tm, n, second, omit);
        else if (!strcmp(cmd, "invisonumb")) ret = str_AuxInvIsoSp3Numb(&cg, a, b, &sb, &ovf, ot, tm, n, second, omit);
        else if (!strcmp(cmd, "crv")) ret = str_AuxChargeRadVal(a, &sb, &ovf, ot, tm, n, mult);
        else if (!strcmp(cmd, "trans")) {
            AT_NUMB *tn = NULL, *ts = NULL;
            bin_AuxTautTrans(a, b, &tn, &ts, ot, n);
            ret = str_AuxTautTrans(&cg, tn, ts, &sb, &ovf, tm, n);
        }
        printf("%d %d [%s]\n", ret, ovf, sb.pStr ? sb.pStr : "");
    }
    return 0;
}
"""

R = random.Random(4242)


def pool_item(r):
    c, a = p3c.mk_pool_item(r)
    for k in (5, 6, 7):
        if r.random() < 0.6:
            a["arrs"][k] = list(a["arrs"][r.choice([4, 4, 5, 6, 7])])
    return c, a


def main():
    cases = []
    for _ in range(7000):
        cmd = R.choice(["numb", "invnumb", "isonumb", "invisonumb", "crv", "trans"])
        pool = [pool_item(R) for _ in range(R.randrange(1, 4))]
        n = R.randrange(1, 5)
        a = p3c.sorts(R, pool, n)
        b = p3c.sorts(R, pool, n)
        ords = " ".join(str(R.randrange(0, n)) for _ in range(2 * n))
        cases.append(f"{cmd} {n} {R.choice([0, 1, 2, 3, 4])} {R.choice([0, 1, 1, 1])} {a} {b} "
                     f"{R.choice([0, 2])} {R.choice([0, 1])} {R.choice([0, 1, 1])} {R.choice([0, 1])} {ords}")
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
    path = ROOT / "tests/data/inchi_prt3d.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
