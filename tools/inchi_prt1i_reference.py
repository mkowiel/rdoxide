"""Generates tests/data/inchi_prt1i.json from the InChI C code: the whole
InChI serialization of ichiprt1.c (`OutputINChI2`, `OutputINChI1`) on random
multi-component structures (mobile-H/fixed-H pairs, stereo, isotopes,
AuxInfo, reversibility). Links the upstream InChI base library (see
inchi_markalt_reference.py).

    python tools/inchi_prt1i_reference.py
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
p3 = p3c.p3
ma = p3c.ma
SRC = ma.SRC

HARNESS = p3c.HARNESS.split("int str_AuxEqu(")[0] + r"""
#define inchi_strbuf_reset(b) (zero_stack(), inchi_strbuf_reset(b))
#include "ichiprt1.c"
static void rd_text(char *t) { scanf("%s", t); for (char *p = t; *p; p++) if (*p == '_') *p = ' '; }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
        printf("#CASE %d\n", c);
        int n = rd();
        INCHI_SORT *sorts[2][2];
        for (int k = 0; k < 4; k++) sorts[k / 2][k % 2] = rd_sorts(n);
        int nc[2] = { rd(), rd() }, nnt[2] = { rd(), rd() }, nt[2] = { rd(), rd() };
        INPUT_PARMS ip; memset(&ip, 0, sizeof(ip));
        ip.bAbcNumbers = rd(); ip.bCtPredecessors = rd(); ip.bPolymers = rd(); ip.bLargeMolecules = rd();
        char lbl[64], val[64];
        int has_lbl = rd(), has_val = rd();
        lbl[0] = val[0] = 0;
        if (has_lbl) strcpy(lbl, "lbl");
        if (has_val) strcpy(val, "val");
        ip.pSdfLabel = has_lbl ? lbl : NULL; ip.pSdfValue = has_val ? val : NULL;
        int disc = rd(), otype = rd(), opts = rd(), saveopt = rd(), alloc = rd(), flags = rd(), bor = rd();
        int use_os = rd();
        ORIG_STRUCT os; memset(&os, 0, sizeof(os));
        char atoms[512], bonds[512], coord[512];
        if (use_os) {
            os.num_atoms = rd(); os.n_zy = rd();
            rd_text(atoms); rd_text(bonds); rd_text(coord);
            os.szAtoms = rd() ? atoms : NULL; os.szBonds = rd() ? bonds : NULL; os.szCoord = rd() ? coord : NULL;
        }
        int use_od = rd();
        ORIG_ATOM_DATA od; memset(&od, 0, sizeof(od));
        od.n_zy = rd();
        int use_log = rd();
        INCHI_IOS_STRING sb; inchi_strbuf_init(&sb, alloc, 16);
        INCHI_IOSTREAM out, lg;
        inchi_ios_init(&out, INCHI_IOS_TYPE_STRING, NULL);
        inchi_ios_init(&lg, INCHI_IOS_TYPE_STRING, NULL);
        int dis2 = rd();
        zero_stack();
        int ret;
        if (dis2) ret = OutputINChI2(&cg, &sb, sorts, bor, use_od ? &od : NULL, use_os ? &os : NULL, &ip, disc, otype, opts, nc, nnt, nt, &out, use_log ? &lg : NULL, 7, &flags, (unsigned char) saveopt);
        else ret = OutputINChI1(&cg, &sb, sorts, bor, use_od ? &od : NULL, use_os ? &os : NULL, &ip, disc, otype, opts, nc, nnt, nt, &out, use_log ? &lg : NULL, 7, &flags, (unsigned char) saveopt);
        printf("%d %d [%s] [%s]\n", ret, flags, out.s.pStr ? out.s.pStr : "", lg.s.pStr ? lg.s.pStr : "");
    }
    return 0;
}
"""

R = random.Random(404)
ATOMS = ["3CNO", "5CCHOH", "2C", "4C2O", "1H"]
BONDS = ["s1;", "s1;d2;", ";s1;", "s1;s2;t3;"]
COORD = ["0,0,0;", "1.5,2,3;4,5,6;", "0,0,1;0,0,2;0,0,3;"]


def mk_ostruct(r):
    n = r.randrange(1, 40)
    return (f"{n} {r.choice([0, 0, 1, 2])} " +
            f"{''.join(r.choice('CNOHcZyH0123456789.+-') for _ in range(r.randrange(1, 150)))}".replace(" ", "_") + " " +
            ";".join("s" + str(r.randrange(1, 9)) + "d" + str(r.randrange(1, 9)) for _ in range(r.randrange(1, 40))) + ";" + " " +
            ";".join(f"{r.randrange(-9, 9)}.{r.randrange(0, 999)},{r.randrange(-9, 9)},0" for _ in range(r.randrange(1, 30))) + ";" +
            f" {r.choice([1, 1, 1, 0])} {r.choice([1, 1, 1, 0])} {r.choice([1, 1, 1, 0])}")


def main():
    cases = []
    for _ in range(9000):
        pool = [p3d.pool_item(R) for _ in range(R.randrange(1, 4))]
        n = R.randrange(1, 5)
        sorts = " ".join(p3c.sorts(R, pool, n) for _ in range(4))
        nums = lambda: R.randrange(0, n + 1)
        opts = 0
        for bit, pr in [(1, .25), (4, .12), (8, .3), (0x40, .5), (0x80, .25), (0x200, .1), (0x400, .3), (0x2000, .4), (0x4000, .3), (0x8000, .3)]:
            if R.random() < pr:
                opts |= bit
        use_os = R.choice([0, 1, 1])
        cases.append(" ".join(map(str, [n])) + f" {sorts} {nums()} {nums()} {nums()} {nums()} {nums()} {nums()} "
                     + " ".join(map(str, [R.choice([0, 1]), R.choice([0, 1]), R.choice([0, 0, 1, 2]), R.choice([0, 0, 1]),
                                          R.choice([0, 1]), R.choice([0, 1]),
                                          R.choice([0, 1, 1]), R.choice([0, 1, 2, 3, 3, 3, 4]), opts, R.randrange(0, 64),
                                          R.choice([100, 200, 1000, 5000]), R.choice([0, 0, 4, 8, 16, 32]), R.choice([0, 1]), use_os]))
                     + (" " + mk_ostruct(R) if use_os else "")
                     + f" {R.choice([0, 1])} {R.choice([0, 0, 1, 2])} {R.choice([0, 1])} {R.choice([0, 1])}")
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
    path = ROOT / "tests/data/inchi_prt1i.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
