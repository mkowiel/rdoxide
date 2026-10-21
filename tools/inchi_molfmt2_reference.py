"""Generates tests/data/inchi_molfmt2.json from the InChI C code: the whole
V2000 reading path of mol_fmt1.c (`MolfileReadDataLines` with the property
block: charges, radicals, isotopes, aliases, polymer S-groups, valence and
hydrogen-isotope checks) and `MolfileTreatPseudoElementAtoms` on random and
fuzzed molfiles. Links the upstream InChI base library (see
inchi_markalt_reference.py).

    python tools/inchi_molfmt2_reference.py
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
ma = m1.ma
SRC = ma.SRC

HARNESS = re.search(r'HARNESS = r"""(.*?)"""', (ROOT / "tools/inchi_molfmt1_reference.py").read_text(), re.S).group(1)
HARNESS = HARNESS.split("static void dump_ctab")[0] + r"""
static void dump_ctab(MOL_FMT_CTAB *c) {
    printf("C %d %d %d %d %d ", c->n_atoms, c->n_bonds, c->chiral_flag, c->n_stext_entries, c->n_property_lines);
    hexs(c->version_string);
    printf(" %d\n", c->v3000 != NULL);
}
static void dump_hdr(MOL_FMT_HEADER_BLOCK *h) {
    printf("H "); hexs(h->molname); printf(" "); hexs(h->line2); printf(" "); hexs(h->user_initls); printf(" "); hexs(h->prog_name);
    printf(" %d %d %d %d %d ", h->month, h->day, h->year, h->hour, h->minute); hexs(h->dim_code);
    printf(" %d ", h->scaling_factor1); pd(h->scaling_factor2); printf(" "); pd(h->energy); printf(" %ld ", h->internal_regno); hexs(h->comment);
    printf("\n");
}
static void dump_all(MOL_FMT_CTAB *ct) {
    for (int i = 0; i < ct->n_atoms && i < 300; i++) {
        MOL_FMT_ATOM *a = &ct->atoms[i];
        printf("A %d ", i); pd(a->fx); printf(" "); pd(a->fy); printf(" "); pd(a->fz); printf(" "); hexs(a->symbol);
        printf(" %d %d %d %d %d %d\n", a->mass_difference, a->charge, a->radical, a->stereo_parity, a->valence, a->atom_aliased_flag);
        if (ct->coords) { printf("O %d ", i); hexs(ct->coords[i]); printf("\n"); }
    }
    if (ct->bonds) for (int i = 0; i < ct->n_bonds && i < 300; i++) {
        MOL_FMT_BOND *b = &ct->bonds[i];
        printf("B %d %d %d %d %d\n", i, b->atnum1, b->atnum2, b->bond_type, b->bond_stereo);
    }
    for (int g = 0; g < ct->sgroups.used; g++) {
        MOL_FMT_SGROUP *s = ct->sgroups.group[g];
        printf("G %d %d %d %d %d ", s->id, s->type, s->subtype, s->conn, s->label);
        for (int q = 0; q < 4; q++) { pd(s->xbr1[q]); printf(" "); }
        for (int q = 0; q < 4; q++) { pd(s->xbr2[q]); printf(" "); }
        hexs(s->smt); printf(" %d", s->alist.used);
        for (int q = 0; q < s->alist.used; q++) printf(" %d", s->alist.item[q]);
        printf(" %d", s->blist.used);
        for (int q = 0; q < s->blist.used; q++) printf(" %d", s->blist.item[q]);
        printf("\n");
    }
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int mode, bcoord, polymers, nowarn, pseudo;
        printf("#CASE %d\n", cs);
        scanf("%d %d %d %d %d %s", &mode, &bcoord, &polymers, &nowarn, &pseudo, hx);
        int n = unhex(hx, text);
        INCHI_IOSTREAM ios;
        inchi_ios_init(&ios, INCHI_IOS_TYPE_STRING, NULL);
        ios.s.pStr = (char *) malloc(n + 2);
        memcpy(ios.s.pStr, text, n + 1);
        ios.s.nUsedLength = n; ios.s.nAllocatedLength = n + 2; ios.s.nPtr = 0;
        char err[STR_ERR_LEN]; err[0] = 0;
        int e = 0;
        MOL_FMT_DATA *m;
        MOL_FMT_HEADER_BLOCK hb; MOL_FMT_CTAB ctb;
        if (mode == 0) {
            m = MolfileReadDataLines(&ios, NULL, NULL, bcoord, polymers, &e, err, nowarn);
            if (m) {
                dump_hdr(&m->hdr); dump_ctab(&m->ctab);
                int e2 = 0; char s2[STR_ERR_LEN]; s2[0] = 0;
                int nzz = MolfileTreatPseudoElementAtoms(&m->ctab, pseudo, &e2, s2);
                dump_all(&m->ctab);
                printf("Z %d %d ", nzz, e2); hexs(s2); printf("\n");
            }
        } else {
            MolfileReadDataLines(&ios, &hb, &ctb, bcoord, polymers, &e, err, nowarn);
            dump_hdr(&hb); dump_ctab(&ctb);
        }
        printf("R %d ", e); hexs(err); printf("\n");
        inchi_ios_free_str(&ios);
    }
    return 0;
}
"""

R = random.Random(777)
SYMS = ["C", "N", "O", "H", "Cl", "CL", "Br", "Zz", "*", "D", "T", "Fe", "S", "P", "F", "Zy", "Na", "Xx", "He"]


def sgroup_lines(r, na, nb):
    lines = []
    ids = r.sample(range(1, 12), r.randrange(1, 4))
    clean = r.random() < 0.6
    typ = lambda: r.choice(["SRU", "MON", "COP", "MOD", "CRO", "MER"] if clean else ["SRU", "MON", "XXX", " SRU", "SR"])
    ent = lambda xs: "".join(xs)
    lines.append(f"M  STY{len(ids):3d}" + "".join(f" {i:3d} {typ():<3}" for i in ids))
    if r.random() < 0.6:
        lines.append(f"M  SST{len(ids):3d}" + "".join(f" {i:3d} {r.choice((['ALT', 'RAN', 'BLK', 'BLO'] if clean else ['ALT', 'XYZ'])):<3}" for i in ids))
    if r.random() < 0.6:
        lines.append(f"M  SLB{len(ids):3d}" + "".join(f" {i:3d} {r.randrange(0, 9):3d}" for i in ids))
    if r.random() < 0.6:
        lines.append(f"M  SCN{len(ids):3d}" + "".join(f" {i:3d} {r.choice((['HT', 'HH', 'EU'] if clean else ['HT', 'XX'])):<3}" for i in ids))
    for i in ids:
        if r.random() < 0.7:
            k = r.randrange(1, 8)
            lines.append(f"M  SAL {i:3d}{k:3d}" + "".join(f" {r.randrange(1 if clean else 0, na + (1 if clean else 2)):3d}" for _ in range(k)))
        if r.random() < 0.5:
            k = r.randrange(1, 5)
            lines.append(f"M  SBL {i:3d}{k:3d}" + "".join(f" {r.randrange(1 if clean else 0, nb + (1 if clean else 2)) if nb or not clean else 1:3d}" for _ in range(k)))
        if r.random() < 0.5:
            lines.append(f"M  SDI {i:3d}  4" + "".join(f"{r.uniform(-9, 9):10.4f}" for _ in range(4)))
            if r.random() < 0.7:
                lines.append(f"M  SDI {i:3d}  4" + "".join(f"{r.uniform(-9, 9):10.4f}" for _ in range(4)))
        if r.random() < 0.5:
            lines.append(f"M  SMT {i:3d} {r.choice(['n', 'x', 'abc def', ' pad '])}")
    return lines


def prop_lines(r, na, nb):
    lines = []
    for _ in range(r.randrange(0, 6)):
        k = r.choice(["CHG", "CHG", "RAD", "ISO", "A", "V", "G", "REG", "SKP", "SGR", "SGR", "SGR", "SGR", "UNK"])
        a = lambda: r.randrange(0, na + 2)
        if k in ("CHG", "RAD", "ISO"):
            n = r.choice([1, 1, 2, 3, 8, 9, 0])
            val = {"CHG": lambda: r.randrange(-3, 4), "RAD": lambda: r.randrange(-1, 5), "ISO": lambda: r.randrange(1, 240)}[k]
            lines.append(f"M  {k}{n:3d}" + "".join(f" {a():3d} {val():3d}" for _ in range(r.choice([n, n, max(n - 1, 0)]))))
        elif k == "A":
            lines.append(f"A  {a():3d}")
            lines.append(r.choice(["R", "D", "T", "Me+", "Et--", "  Ph  ", "", "   ", "Xyzabcdef", "C^", "N:", "A B"]))
        elif k == "V":
            lines.append(f"V  {a():3d} value")
        elif k == "G":
            lines.append(f"G  {a():3d}   1")
            lines.append("grp")
        elif k == "REG":
            lines.append(f"M  REG {r.choice(['12', '  34', '5x', '', '-7', '99999999999'])}")
        elif k == "SKP":
            lines.append(f"M  SKP{r.choice([1, 2, 0]):3d}")
            lines.append("skipped M  CHG  1   1   1")
        elif k == "SGR":
            lines += sgroup_lines(r, na, nb)
        else:
            lines.append("M  UNK  1")
    return lines


def mk_mol2(r):
    na = r.choice([2, 3, 5, 8, 12, 1])
    nb = r.choice([1, 2, 4, 7]) if na > 1 else 0
    lines = [r.choice(["mol", "Structure #5", ""]),
             r.choice(["  rdkit     01012012342D", "  InChI 1.07  SDfile Output 3D", ""]), "c"]
    props = prop_lines(r, na, nb)
    vers = r.choice(["V2000"] * 5 + [""])
    nprop = len(props) + 1 + r.choice([0, 0, 0, 0, 1, -1])
    lines.append(f"{na:3d}{nb:3d}  0  0{r.choice([0, 1]):3d}{r.choice([0, 0, 0, 1]):3d}  0  0  0  0{nprop:3d} {vers}")
    for _ in range(na):
        sym = r.choice(SYMS)
        lines.append(f"{r.uniform(-9, 9):10.4f}{r.uniform(-9, 9):10.4f}{r.choice([0, 0, 1]) * r.uniform(-3, 3):10.4f} {sym:<3}"
                     f"{r.choice([0] * 12 + [1, 2, 127, -1, 5]):2d}{r.choice([0, 0, 0, 3, 5, 4, 7, 9]):3d}{r.choice([0, 1, 2, 3]):3d}"
                     f"  0  0{r.choice([0] * 20 + [1, 2, 4, 15, 21]):3d}  0  0  0  0  0")
    for _ in range(nb):
        lines.append(f"{r.randrange(1, na + 1):3d}{r.randrange(1, na + 1):3d}{r.choice([1, 1, 2, 3]):3d}{r.choice([0, 0, 1, 4, 6]):3d}  0  0")
    for _ in range(r.choice([0, 0, 0, 1])):
        lines += ["stext1", "stext2"]
    lines += props
    if r.random() < 0.9:
        lines.append("M  END")
    if r.random() < 0.3:
        lines.append("$$$$")
    p = r.choice([0, 0, 0, 0.01, 0.03])
    lines = [m1.fuzz(r, l, p) if i > 3 else l for i, l in enumerate(lines)]
    if r.random() < 0.1:
        lines = lines[:r.randrange(0, len(lines))]
    if r.random() < 0.05 and lines:
        lines.insert(r.randrange(0, len(lines)), "$$$$")
    nl = r.choice(["\n", "\n", "\n", "\r\n"])
    return nl.join(lines) + nl


def main():
    cases = []
    for _ in range(2200):
        cases.append(f"{R.choice([0, 0, 0, 1])} {R.choice([0, 1])} {R.choice([0, 1, 1])} {R.choice([0, 1])} {R.choice([0, 1])} "
                     f"{mk_mol2(R).encode().hex()}")
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-ffp-contract=off", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
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
    path = ROOT / "tests/data/inchi_molfmt2.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
