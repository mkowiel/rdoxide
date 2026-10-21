"""Generates tests/data/inchi_molfmt3.json from the InChI C code: the V3000
Molfile reader (mol_fmt3.c: CTAB head, atoms with key=value pairs, bonds,
haptic bonds, S-group/3D/collection blocks, star atoms, continuation lines)
inside `MolfileReadDataLines` on random and fuzzed V3000 molfiles. Links the
upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_molfmt3_reference.py
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


m1 = load("molfmt1", "tools/inchi_molfmt1_reference.py")
m2 = load("molfmt2", "tools/inchi_molfmt2_reference.py")
ma = m1.ma
SRC = ma.SRC

# the harness of molfmt2 without the V3000 stubs (the real mol_fmt3.c is linked)
HARNESS = "\n".join(l for l in m2.HARNESS.split("\n") if not l.startswith(("int MolfileV3000", "int DeleteMolfileV3000")))
HARNESS = HARNESS.replace("static void dump_all(MOL_FMT_CTAB *ct) {", r"""static void dump_lists(NUM_LISTS *l, int extra) {
    if (!l) { printf(" -"); return; }
    printf(" %d", l->used);
    for (int i = 0; i < l->used; i++) {
        int n = extra == 3 ? l->lists[i][2] : l->lists[i][1];
        printf(" [");
        for (int k = 0; k < n + extra && k < 40; k++) printf("%d,", l->lists[i][k]);
        printf("]");
    }
}
static void dump_v3(MOL_FMT_CTAB *ct) {
    MOL_FMT_v3000 *v = ct->v3000;
    if (!v) return;
    printf("V %d %d %d %d %d %d %d %d %d %d %d", v->n_non_star_atoms, v->n_star_atoms, v->n_sgroups, v->n_3d_constraints, v->n_collections,
           v->n_non_haptic_bonds, v->n_haptic_bonds, v->n_steabs, v->n_sterel, v->n_sterac, v->atom_index_orig != NULL);
    if (v->atom_index_orig) for (int i = 0; i < v->n_non_star_atoms + v->n_star_atoms; i++) printf(" %d:%d", v->atom_index_orig[i], v->atom_index_fin[i]);
    printf("\nW"); dump_lists(v->haptic_bonds, 3);
    printf("\nS"); dump_lists(v->steabs, 2); dump_lists(v->sterel, 2); dump_lists(v->sterac, 2);
    printf("\n");
}
static void dump_all(MOL_FMT_CTAB *ct) {
    dump_v3(ct);""", 1)
HARNESS = HARNESS.replace("int main(void) {", "int main(void) {\n    /* indices into atom_index arrays are limited to the allocated n_atoms */", 1)

R = random.Random(30003)
SYMS = ["C", "N", "O", "H", "Cl", "CL", "Br", "Zz", "*", "*", "D", "T", "Fe", "S", "P", "Na"]


def num(r):
    return r.choice([f"{r.uniform(-9, 9):.4f}", f"{r.uniform(-900, 900):.3f}", "0", "1.5", "-2", "0.0001", "12345.678"])


def atom_line(r, idx):
    sym = r.choice(SYMS)
    kv = []
    for k, vals in [("CHG", [-1, 1, 2, -3, 200]), ("RAD", [1, 2, 3]), ("CFG", [1, 2, 3, 0]), ("MASS", [13, 2, 18, 35, 99999999]),
                    ("VAL", [-1, 1, 3, 4, 0]), ("HCOUNT", [1, -1]), ("STBOX", [1]), ("INVRET", [1]), ("EXACHG", [0]), ("SUBST", [2]),
                    ("UNSAT", [1]), ("RBCNT", [3]), ("ATTCHPT", [1]), ("SEQID", [5]), ("RGROUPS", ["(1 2)"]), ("ATTCHORD", ["(2 1)"]),
                    ("CLASS", ["x"]), ("FOO", [1])]:
        if r.random() < 0.12:
            kv.append(f"{k}={r.choice(vals)}")
    r.shuffle(kv)
    return f"{idx} {sym} {num(r)} {num(r)} {num(r)} {r.choice([0, 0, 0, 3])} " + " ".join(kv)


def bond_line(r, i, labels):
    a1, a2 = r.choice(labels + [0, 99]), r.choice(labels + [0, 99])
    kv = []
    if r.random() < 0.2:
        kv.append(f"CFG={r.choice([1, 2, 3, 0])}")
    if r.random() < 0.08:
        kv.append(r.choice(["TOPO=1", "RXCTR=1", "STBOX=1", "DISP=x", "ATTACH=ANY"]))
    if r.random() < 0.12:
        m = r.randrange(0, 4)
        kv.append(f"ENDPTS=({m} " + " ".join(str(r.choice(labels + [99])) for _ in range(m)) + ")" + r.choice([" ATTACH=ALL", " ATTACH=ALL", "", " ATTACH=ANY"]))
    r.shuffle(kv)
    return f"{i} {r.choice([1, 1, 2, 3, 4, 9, 0])} {a1} {a2} " + " ".join(kv)


def fuzz_content(r, s, p):
    out = []
    for ch in s:
        x = r.random()
        if x < p * 0.4:
            out.append(r.choice("0123456789 -+.e=()*A/"))
        elif x < p * 0.55:
            continue
        elif x < p * 0.7:
            out += [ch, r.choice(" 1-.a")]
        else:
            out.append(ch)
    return "".join(out)


def mk_mol3(r):
    na = r.choice([0, 1, 2, 3, 4, 6, 9])
    nb = r.choice([0, 1, 2, 3, 5]) if na > 1 else r.choice([0, 0, 1])
    labels = r.sample(range(1, 40), na) if r.random() < 0.5 else list(range(1, na + 1))
    L = ["Structure #7", "  InChI 1.07  SDfile Output 3D", "c",
         "  0  0  0     0  0            999 V3000", "M  V30 BEGIN CTAB"]
    L.append(f"M  V30 COUNTS {na} {nb} {r.choice([0, 0, 1])} 0 {r.choice([0, 1])}")
    L.append("M  V30 BEGIN ATOM")
    p = r.choice([0, 0, 0, 0, 0.01, 0.03])
    content = [atom_line(r, labels[i]) for i in range(na)]
    for c in content:
        c = fuzz_content(r, c, p)
        if r.random() < 0.1 and len(c) > 6:
            k = r.randrange(3, len(c) - 2)
            L += ["M  V30 " + c[:k] + "-", "M  V30 " + c[k:]]
        else:
            L.append("M  V30 " + c)
    L.append("M  V30 END ATOM")
    if nb:
        L.append("M  V30 BEGIN BOND")
        for i in range(nb):
            L.append("M  V30 " + fuzz_content(r, bond_line(r, i + 1, labels), p))
        L.append("M  V30 END BOND")
    if r.random() < 0.2:
        L += ["M  V30 BEGIN SGROUP", "M  V30 1 SRU 1 ATOMS=(1 1)", "M  V30 END SGROUP"]
    if r.random() < 0.1:
        L += ["M  V30 BEGIN OBJ3D", "M  V30 END OBJ3D"] if r.random() < 0.8 else ["M  V30 BEGIN OBJ3D", "M  V30 junk"]
    if r.random() < 0.1:
        L.append("M  V30 LINKNODE")
    if r.random() < 0.35 and labels:
        L.append("M  V30 BEGIN COLLECTION")
        for _ in range(r.randrange(0, 4)):
            kind = r.choice(["STEABS", "STEREL1", "STEREL2", "STERAC1", "STERAC12", "STEREL", "XXX", "STEABS"])
            m = r.randrange(0, 4)
            at = " ".join(str(r.choice(labels + [99])) for _ in range(m))
            L.append(f"M  V30 MDLV30/{kind} ATOMS=({m} {at})" if r.random() < 0.9 else f"M  V30 MDLV30/{kind} JUNK")
        L.append("M  V30 END COLLECTION")
    if r.random() < 0.95:
        L.append("M  V30 END CTAB")
    L.append("M  END")
    if r.random() < 0.1:
        L = L[:r.randrange(5, len(L))]
        # never leave an SGROUP block open (the C loop would not terminate)
        if "M  V30 BEGIN SGROUP" in L and "M  V30 END SGROUP" not in L:
            L.append("M  V30 END SGROUP")
    nl = r.choice(["\n", "\n", "\n", "\r\n"])
    return nl.join(L) + r.choice([nl, nl, ""])


def main():
    cases = []
    for _ in range(3500):
        cases.append(f"0 {R.choice([0, 1])} {R.choice([0, 1, 1])} {R.choice([0, 1])} {R.choice([0, 1])} {mk_mol3(R).encode().hex()}")
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-ffp-contract=off", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        kept, out = [], []
        for c in cases:
            try:
                r1 = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=10)
            except subprocess.TimeoutExpired:
                continue
            if r1.returncode == 0:
                kept.append(c)
                out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode(errors="replace").split("\n")[1:]
        print(len(cases) - len(kept), "cases dropped")
        cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    # a coordinate string longer than MOL_COORD overflows into the next slot (C undefined behaviour)
    data = [x for x in data if all(len(l.split()[2]) <= 62 for l in x["out"] if l.startswith("O "))]
    path = ROOT / "tests/data/inchi_molfmt3.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
