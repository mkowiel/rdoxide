"""Generates tests/data/inchi_prt3b.json from the InChI C code: the stereo and
isotopic whole-structure substring builders of ichiprt3.c (`str_Sp2`,
`str_Sp3`, `str_StereoAbsInv`, `str_IsoAtoms`, `str_IsoSp2`, `str_IsoSp3`,
`str_IsoStereoAbsInv`) on random multi-component structures. Links the
upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_prt3b_reference.py
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
m2 = p3.m2
ma = p3.ma
SRC = ma.SRC

HARNESS = p3.HARNESS.split("int main(void)")[0] + r"""
int str_Sp2( INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int second, int omit, int mult );
int str_Sp3( INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int bRelRac, int second, int omit, int mult );
int str_StereoAbsInv( INCHI_SORT *a, INCHI_IOS_STRING *sb, int *ovf, int ot, int n );
int str_IsoAtoms( INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int abc, int second, int omit, int mult );
int str_IsoSp2( INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int second, int omit, int mult );
int str_IsoSp3( INCHI_SORT *a, INCHI_SORT *b, INCHI_IOS_STRING *sb, int *ovf, int ot, int TM, int n, int bRelRac, int second, int omit, int mult );
int str_IsoStereoAbsInv( INCHI_SORT *a, INCHI_IOS_STRING *sb, int *ovf, int ot, int n );
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
        int tm = rd(), second = rd(), omit = rd(), abc = rd();
        printf("#CASE %d\n", c);
        zero_stack();
        if (!strcmp(cmd, "sp2")) ret = str_Sp2(a, b, &sb, &ovf, ot, tm, n, second, omit, mult);
        else if (!strcmp(cmd, "sp3")) ret = str_Sp3(a, b, &sb, &ovf, ot, tm, n, 0, second, omit, mult);
        else if (!strcmp(cmd, "abs")) ret = str_StereoAbsInv(a, &sb, &ovf, ot, n);
        else if (!strcmp(cmd, "isoatoms")) ret = str_IsoAtoms(a, b, &sb, &ovf, ot, tm, n, abc, second, omit, mult);
        else if (!strcmp(cmd, "isosp2")) ret = str_IsoSp2(a, b, &sb, &ovf, ot, tm, n, second, omit, mult);
        else if (!strcmp(cmd, "isosp3")) ret = str_IsoSp3(a, b, &sb, &ovf, ot, tm, n, 0, second, omit, mult);
        else if (!strcmp(cmd, "isoabs")) ret = str_IsoStereoAbsInv(a, &sb, &ovf, ot, n);
        printf("%d %d [%s]\n", ret, ovf, sb.pStr ? sb.pStr : "");
    }
    return 0;
}
"""

R = random.Random(9091)


def stereo(r, n):
    if r.random() < 0.3:
        return None
    ns = sorted(r.sample(range(1, max(n, 2) + 1), min(r.choice([0, 1, 2]), max(n, 2))))
    bonds = sorted({(r.randrange(1, 4), r.randrange(4, 7)) for _ in range(r.choice([0, 0, 1, 2]))})
    return {"sc": [[a, r.choice([1, 2, 3, 4])] for a in ns],
            "sb": [[a, b, r.choice([1, 2])] for a, b in bonds],
            "inv": r.choice([0, 0, 1, -1, 2])}


def base(r):
    c = p3.base(r)
    n = c["n"]
    c["st"] = stereo(r, n)
    c["sti"] = stereo(r, n) if r.random() < 0.7 else (copy.deepcopy(c["st"]) if r.random() < 0.5 else None)
    c["iso"] = [[r.randrange(1, n + 1), r.choice([0, 1, -1]), r.choice([0, 1]), r.choice([0, 1]), r.choice([0, 1])]
                for _ in range(r.choice([0, 0, 1, 2]))]
    c["isotg"] = [[r.randrange(1, 3)] + [r.choice([0, 1]) for _ in range(3)] for _ in range(r.choice([0, 0, 1]))]
    return c


def sorts(r, pool, n):
    out = []
    for _ in range(n):
        c = copy.deepcopy(r.choice(pool))
        if r.random() < 0.25:
            # a variant: change the stereo or isotopes only
            c["st"] = stereo(r, c["n"])
            if r.random() < 0.5:
                c["sti"] = copy.deepcopy(c["st"])
        out.append(p3.slot(r, c))
    return out


def main():
    cases = []
    for _ in range(7000):
        cmd = R.choice(["sp2", "sp3", "abs", "isoatoms", "isosp2", "isosp3", "isoabs"])
        pool = [base(R) for _ in range(R.randrange(1, 4))]
        n = R.randrange(1, 5)
        a = sorts(R, pool, n)
        b = sorts(R, pool, n)
        cases.append(f"{cmd} {n} {R.choice([0, 1, 2, 3, 4])} {R.choice([0, 1, 1, 1])} {p3.tok_sorts(a)} {p3.tok_sorts(b)} "
                     f"{R.choice([0, 2])} {R.choice([0, 1])} {R.choice([0, 1, 1])} {R.choice([0, 1])}")
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
    path = ROOT / "tests/data/inchi_prt3b.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
