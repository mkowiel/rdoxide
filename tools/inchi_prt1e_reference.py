"""Generates tests/data/inchi_prt1e.json from the InChI C code: the original
structure writers of ichiprt1.c (`CleanOrigCoord`, `WriteOrigCoord`,
`WriteOrigAtoms`). The harness includes ichiprt1.c itself (they are static).
Links the upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_prt1e_reference.py
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


p3c = load("prt3c", "tools/inchi_prt3c_reference.py")
p3 = p3c.p3
ma = p3c.ma
SRC = ma.SRC

HARNESS = p3c.HARNESS.split("int str_AuxEqu(")[0] + r"""
#include "ichiprt1.c"
static void rd_coord(MOL_COORD c) {
    char t[64]; scanf("%s", t);
    memset(c, 'X', sizeof(MOL_COORD));
    int n = (int) strlen(t);
    for (int i = 0; i < n && i < (int) sizeof(MOL_COORD); i++) c[i] = t[i] == '_' ? ' ' : t[i] == '~' ? '\0' : t[i];
}
static void pr_coord(MOL_COORD c) {
    for (int i = 0; i < (int) sizeof(MOL_COORD); i++) printf("%02x", (unsigned char) c[i]);
}
int main(void) {
    int ncases;
    char cmd[32];
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
        STRUCT_DATA sd; memset(&sd, 0, sizeof(sd));
        scanf("%s", cmd);
        printf("#CASE %d\n", c);
        if (!strcmp(cmd, "clean")) {
            MOL_COORD mc; rd_coord(mc);
            int delim = rd();
            int nz = CleanOrigCoord(mc, delim);
            printf("%d ", nz); pr_coord(mc); printf("\n");
        } else if (!strcmp(cmd, "wcoord")) {
            int n = rd(), start = rd(), buf_len = rd();
            MOL_COORD *mc = calloc(n + 1, sizeof(MOL_COORD));
            for (int k = 0; k < n; k++) rd_coord(mc[k]);
            char buf[1024]; memset(buf, 'Q', sizeof(buf));
            int i = start;
            int len = WriteOrigCoord(n, mc, &i, buf, buf_len);
            printf("%d %d [%s]\n", len, i, buf);
        } else {
            int n = rd(), start = rd(), buf_len = rd();
            sd.bChiralFlag = rd();
            inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
            for (int k = 0; k < n; k++) {
                char el[16]; scanf("%s", el);
                strcpy(at[k].elname, el);
                at[k].el_number = rd();
                at[k].valence = rd();
                for (int m = 0; m < at[k].valence; m++) at[k].neighbor[m] = rd();
                for (int m = 0; m < at[k].valence; m++) at[k].bond_type[m] = rd();
                at[k].orig_at_number = rd();
                at[k].chem_bonds_valence = rd();
                at[k].num_H = rd();
                for (int m = 0; m < 3; m++) at[k].num_iso_H[m] = rd();
                at[k].iso_atw_diff = rd();
                at[k].charge = rd();
                at[k].radical = rd();
                at[k].p_parity = rd();
                for (int m = 0; m < 4; m++) at[k].p_orig_at_num[m] = rd();
            }
            char buf[1024]; memset(buf, 'Q', sizeof(buf));
            int i = start;
            int len = WriteOrigAtoms(&cg, n, at, &i, buf, buf_len, &sd);
            if (len == 0 && buf_len > 0 && buf[0] == 'Q') buf[0] = 0;
            buf[len] = 0;
            printf("%d %d [%s]\n", len, i, buf);
        }
    }
    return 0;
}
"""

R = random.Random(555)
ELS = [("H", 1), ("C", 6), ("C", 6), ("C", 6), ("N", 7), ("O", 8), ("F", 9), ("Cl", 17), ("S", 16), ("P", 15),
       ("Na", 11), ("Br", 35), ("Fe", 26), ("B", 5), ("Si", 14), ("I", 53), ("Zn", 30)]


def mk_coord(r):
    def one():
        k = r.random()
        if k < 0.1:
            return "0.0000"
        if k < 0.2:
            return "%.4f" % r.uniform(-9, 9)
        if k < 0.3:
            return r.choice(["1.5e-3", "-2.5E+02", "3.0d1", "1.20E-05", "0.0000001", "+1.500", "-0", "12e", "1e+0", "5.", ".5", "100.000", "0.00100", "007.5", "-0.50", "1.2e-007", "9.99999e9"])
        return "%.*f" % (r.randrange(0, 5), r.uniform(-99, 99))
    fs = []
    for _ in range(3):
        t = one()[:10]
        pad = r.choice([" ", " ", "~"])
        t = (" " * r.choice([0, 0, 1, 2])) + t
        t = t[:10]
        t = t + ("_" if pad == " " else "~") * (10 - len(t)) if r.random() < 0.9 else t
        fs.append(t.replace(" ", "_"))
    return "".join(fs)


def mk_atoms(r, n):
    adj = [[] for _ in range(n)]
    bt = {}
    for v in range(1, n):
        if r.random() < 0.9:
            u = r.randrange(0, v)
            adj[v].append(u)
            adj[u].append(v)
            bt[(u, v)] = bt[(v, u)] = r.choice([1, 1, 1, 2, 3, 4])
    out = []
    for v in range(n):
        el, num = r.choice(ELS)
        nb = sorted(adj[v])
        valence = len(nb)
        types = [bt[(v, u)] for u in nb]
        orig = v + 1 if r.random() < 0.95 else r.randrange(0, 8)
        cbv = sum(types) + (r.choice([0, 0, 0, 1, -1]))
        numh = r.choice([0, 0, 1, 2, 3, 4])
        iso = [r.choice([0, 0, 0, 1, 2]) for _ in range(3)]
        diff = r.choice([0, 0, 0, 1, -1, 2, -2, 5])
        ch = r.choice([0, 0, 0, 1, -1, 2, -3])
        rad = r.choice([0, 0, 0, 0, 1, 2, 3])
        parity = r.choice([0, 0, 0, 1, 2, 3, 4])
        cand = [u + 1 for u in nb] + [v + 1]
        po = [r.choice(cand) if r.random() < 0.9 else r.randrange(0, 9) for _ in range(4)]
        if r.random() < 0.6:
            po = (cand + [0] * 4)[:4]
            r.shuffle(po)
        out.append(" ".join(map(str, [el, num, valence, *nb, *types, orig, cbv, numh, *iso, diff, ch, rad, parity, *po])))
    return out


def main():
    cases = []
    for _ in range(9000):
        cmd = R.choice(["clean", "wcoord", "atoms", "atoms"])
        if cmd == "clean":
            cases.append(f"clean {mk_coord(R)} {R.choice([44, 59, 44, 32])}")
        elif cmd == "wcoord":
            n = R.randrange(0, 6)
            cases.append(f"wcoord {n} {R.randrange(0, n + 1)} {R.choice([12, 30, 60, 100, 300])} " + " ".join(mk_coord(R) for _ in range(n)))
        else:
            n = R.randrange(1, 7)
            cases.append(f"atoms {n} {R.choice([0, 0, 0, 1, 2])} {R.choice([8, 20, 40, 100, 300])} {R.choice([0, 1, 2, 3])} " + " ".join(mk_atoms(R, n)))
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
    path = ROOT / "tests/data/inchi_prt1e.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
