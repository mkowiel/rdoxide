"""Generates tests/data/inchi_prt3c.json from the InChI C code: the AuxInfo
equivalence and inversion builders of ichiprt3.c (`str_AuxEqu`,
`str_AuxIsoEqu`, `str_AuxTgroupEqu`, `str_AuxIsoTgroupEqu`, `str_AuxInvSp3`,
`str_AuxInvIsoSp3`) on random multi-component structures. Links the upstream
InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_prt3c_reference.py
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


p3 = load("prt3a", "tools/inchi_prt3a_reference.py")
p3b = load("prt3b", "tools/inchi_prt3b_reference.py")
p2 = load("prt2", "tools/inchi_prt2_reference.py")
m2 = p3.m2
ma = p3.ma
SRC = ma.SRC

_h2 = p2.HARNESS
AUX_READERS = _h2[_h2.index("static AT_NUMB *rd_at"):_h2.index("static INChI *rd_iso_inchi")]

HARNESS = p3.HARNESS.split("static INCHI_SORT *rd_sorts")[0] + AUX_READERS + r"""
static INCHI_SORT *rd_sorts(int n) {
    INCHI_SORT *s = calloc(n + 1, sizeof(INCHI_SORT));
    for (int i = 0; i < n; i++) { s[i].pINChI[0] = rd_inchi(); s[i].pINChI[1] = rd_inchi(); s[i].pINChI_Aux[0] = rd_aux(); s[i].pINChI_Aux[1] = rd_aux(); s[i].ord_number = i; }
    return s;
}
int str_AuxEqu( INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int second, int omit, int mult );
int str_AuxIsoEqu( INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int second, int omit, int mult );
int str_AuxTgroupEqu( INCHI_SORT *a, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int mult );
int str_AuxIsoTgroupEqu( INCHI_SORT *a, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int omit, int mult );
int str_AuxInvSp3( INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int second, int omit, int mult );
int str_AuxInvIsoSp3( INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int second, int omit, int mult );
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
        printf("#CASE %d\n", c);
        zero_stack();
        if (!strcmp(cmd, "equ")) ret = str_AuxEqu(a, b, &sb, &ovf, ot, tm, n, second, omit, mult);
        else if (!strcmp(cmd, "isoequ")) ret = str_AuxIsoEqu(a, b, &sb, &ovf, ot, tm, n, second, omit, mult);
        else if (!strcmp(cmd, "tg")) ret = str_AuxTgroupEqu(a, &sb, &ovf, ot, tm, n, mult);
        else if (!strcmp(cmd, "isotg")) ret = str_AuxIsoTgroupEqu(a, &sb, &ovf, ot, tm, n, omit, mult);
        else if (!strcmp(cmd, "invsp3")) ret = str_AuxInvSp3(a, b, &sb, &ovf, ot, tm, n, second, omit, mult);
        else if (!strcmp(cmd, "inviso")) ret = str_AuxInvIsoSp3(a, b, &sb, &ovf, ot, tm, n, second, omit, mult);
        printf("%d %d [%s]\n", ret, ovf, sb.pStr ? sb.pStr : "");
    }
    return 0;
}
"""

R = random.Random(777)


def mk_pool_item(r):
    c = p3b.base(r)
    n = c["n"]
    t = c["taut"][0] if c["taut"] else 0
    a = p2.mk_aux(r)
    a["n"] = n
    a["t"] = t
    a["arrs"][0] = p2.gen_equ(r, n)
    a["arrs"][1] = p2.gen_equ(r, max(t, 1))
    a["arrs"][2] = p2.gen_equ(r, n)
    a["arrs"][3] = p2.gen_equ(r, max(t, 1))
    a["present"] = [r.random() < 0.9 for _ in range(8)]
    a["iso"] = r.choice([0, 1, 1])
    return c, a


def sorts(r, pool, n):
    toks = []
    for _ in range(n):
        c, a = copy.deepcopy(r.choice(pool))
        if r.random() < 0.25:
            c["st"] = p3b.stereo(r, c["n"])
            if r.random() < 0.5:
                c["sti"] = copy.deepcopy(c["st"])
            a = p2.mut_aux(r, a)
        s0, s1 = p3.slot(r, c)
        t = (m2.tok(s0) if s0 else "0") + " " + (m2.tok(s1) if s1 else "0")
        a0 = p2.tok_aux(a) if s0 and r.random() < 0.92 else "0"
        a1 = p2.tok_aux(a) if s1 and r.random() < 0.92 else "0"
        toks.append(f"{t} {a0} {a1}")
    return " ".join(toks)


def main():
    cases = []
    for _ in range(7000):
        cmd = R.choice(["equ", "isoequ", "tg", "isotg", "invsp3", "inviso"])
        pool = [mk_pool_item(R) for _ in range(R.randrange(1, 4))]
        n = R.randrange(1, 5)
        a = sorts(R, pool, n)
        b = sorts(R, pool, n)
        cases.append(f"{cmd} {n} {R.choice([0, 1, 2, 3, 4])} {R.choice([0, 1, 1, 1])} {a} {b} "
                     f"{R.choice([0, 2])} {R.choice([0, 1])} {R.choice([0, 1, 1])}")
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
    path = ROOT / "tests/data/inchi_prt3c.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
