"""Generates tests/data/inchi_prt1f.json from the InChI C code: the bond
writer of ichiprt1.c (`WriteOrigBonds`, with stereo bonds, cumulenes and
metals). The harness includes ichiprt1.c itself (it is static). Links the
upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_prt1f_reference.py
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
ma = p3c.ma
SRC = ma.SRC

HARNESS = p3c.HARNESS.split("int str_AuxEqu(")[0] + r"""
#include "ichiprt1.c"
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
        STRUCT_DATA sd; memset(&sd, 0, sizeof(sd));
        printf("#CASE %d\n", c);
        int n = rd(), start = rd(), buf_len = rd();
        inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
        for (int k = 0; k < n; k++) {
            at[k].el_number = rd();
            at[k].valence = rd();
            for (int m = 0; m < at[k].valence; m++) at[k].neighbor[m] = rd();
            for (int m = 0; m < at[k].valence; m++) at[k].bond_type[m] = rd();
            for (int m = 0; m < at[k].valence; m++) at[k].bond_stereo[m] = rd();
            at[k].chem_bonds_valence = rd();
            for (int m = 0; m < 3; m++) at[k].sb_parity[m] = rd();
            for (int m = 0; m < 3; m++) at[k].sb_ord[m] = rd();
            for (int m = 0; m < 3; m++) at[k].sn_ord[m] = rd();
        }
        char buf[2048]; memset(buf, 'Q', sizeof(buf));
        int i = start;
        int len = WriteOrigBonds(&cg, n, at, &i, buf, buf_len, &sd);
        if (len == 0 && buf[0] == 'Q') buf[0] = 0;
        buf[len] = 0;
        printf("%d %d [%s] [%s]\n", len, i, buf, sd.pStrErrStruct);
    }
    return 0;
}
"""

R = random.Random(8128)
ELS = [6, 6, 6, 7, 8, 1, 17, 11, 26, 30, 16, 15]


def mk_atoms(r):
    """a random graph with stereo double bonds and cumulenes"""
    nodes = []  # element numbers
    bonds = {}  # frozenset -> type
    stereo = {}  # node -> list of (other, sn_neighbor, parity)

    def new(el=None):
        nodes.append(el if el is not None else r.choice(ELS))
        return len(nodes) - 1

    def bond(a, b, t):
        bonds[frozenset((a, b))] = t

    for _ in range(r.randrange(0, 3)):
        u = new()
        v = new()
        mids = [new(6) for _ in range(r.choice([0, 0, 1, 2]))]
        chain = [u] + mids + [v]
        for a, b in zip(chain, chain[1:]):
            bond(a, b, 2)
        for end in (u, v):
            for _ in range(r.choice([1, 1, 2])):
                s = new()
                bond(end, s, r.choice([1, 1, 1, 2]))
        for end, other in ((u, chain[1]), (v, chain[-2])):
            nb = [x for x in nodes_nb(bonds, end) if x != other]
            if nb:
                stereo.setdefault(end, []).append((other, r.choice(nb), r.choice([0, 1, 1, 2, 2, 3, 4, 1 | 8 | 16, 2 | 8 | 16, 1 | 16, 3 | 24, 0o71])))
    for _ in range(r.randrange(0, 5)):
        a = new()
        if len(nodes) > 1 and r.random() < 0.9:
            b = r.randrange(0, len(nodes) - 1)
            bond(a, b, r.choice([1, 1, 1, 2, 3, 4]))
    n = len(nodes)
    perm = list(range(n))
    r.shuffle(perm)  # new index of node v is perm[v]
    atoms = [None] * n
    for v in range(n):
        nb = nodes_nb(bonds, v)
        r.shuffle(nb)
        types = [bonds[frozenset((v, u))] for u in nb]
        st = []
        for u in nb:
            t = bonds[frozenset((v, u))]
            if t == 1:
                st.append(r.choice([0, 0, 0, 1, -1, 6, -6, 4, -4]))
            elif t == 2:
                st.append(r.choice([0, 0, 0, 3, -3]))
            else:
                st.append(0)
        sb = [0, 0, 0]
        so = [0, 0, 0]
        sn = [0, 0, 0]
        for m, (other, snn, par) in enumerate(stereo.get(v, [])[:3]):
            sb[m] = par
            so[m] = nb.index(other)
            sn[m] = nb.index(snn)
        cbv = sum(types) if r.random() < 0.97 else sum(types) + 1
        atoms[perm[v]] = (nodes[v], len(nb), [perm[u] for u in nb], types, st, cbv, sb, so, sn)
    return atoms


def nodes_nb(bonds, v):
    return [next(iter(k - {v})) for k in bonds if v in k]


def tok(atoms):
    out = []
    for el, val, nb, types, st, cbv, sb, so, sn in atoms:
        out.append(" ".join(map(str, [el, val, *nb, *types, *st, cbv, *sb, *so, *sn])))
    return " ".join(out)


def main():
    cases = []
    for _ in range(9000):
        atoms = mk_atoms(R)
        n = len(atoms)
        cases.append(f"{n} {R.choice([0, 0, 0, 1, 2])} {R.choice([10, 20, 40, 100, 400])} {tok(atoms)}")
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
    path = ROOT / "tests/data/inchi_prt1f.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
