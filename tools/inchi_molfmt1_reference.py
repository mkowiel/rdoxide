"""Generates tests/data/inchi_molfmt1.json from the InChI C code: the V2000
Molfile reader of mol_fmt1.c / mol_fmt2.c (field reader, header, counts line,
atom and bond blocks, STEXT block, structure number, coordinate
normalization) on random and fuzzed molfile text. Links the upstream InChI
base library (see inchi_markalt_reference.py).

    python tools/inchi_molfmt1_reference.py
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


ma = load("markalt", "tools/inchi_markalt_reference.py")
SRC = ma.SRC

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mode.h"
#include "mol_fmt.h"
#include "ichierr.h"
#include "util.h"
#include "ichi_io.h"
/* the V3000 reader (mol_fmt3.c) is not under test: stubs keep the link small */
int CreateCompositeNormAtom() { return 0; }
int MolfileV3000Init(MOL_FMT_CTAB *c, char *e) { return 0; }
int MolfileV3000ReadCTABBeginAndCountsLine(MOL_FMT_CTAB *c, INCHI_IOSTREAM *f, char *e) { return 0; }
int MolfileV3000ReadAtomsBlock(MOL_FMT_CTAB *c, INCHI_IOSTREAM *f, int r, char *e) { return 0; }
int MolfileV3000ReadBondsBlock(MOL_FMT_CTAB *c, INCHI_IOSTREAM *f, int r, char *e) { return 0; }
int MolfileV3000ReadTailOfCTAB(MOL_FMT_CTAB *c, INCHI_IOSTREAM *f, int r, char *e) { return 0; }
int DeleteMolfileV3000Info(MOL_FMT_v3000 *v) { return 0; }
#include "mol_fmt1.c"
#include "mol_fmt2.c"

static void hex(const char *s, int n) { for (int i = 0; i < n; i++) printf("%02x", (unsigned char) s[i]); if (!n) printf("-"); }
static void hexs(const char *s) { hex(s, (int) strlen(s)); }
static void pd(double d) { if (isnan(d)) { printf("nan"); return; } unsigned long long u; memcpy(&u, &d, 8); printf("%016llx", u); }
static int unhex(const char *h, char *out) {
    int n = 0;
    out[0] = 0;
    if (h[0] == '-') return 0;
    for (; h[0] && h[1]; h += 2) { unsigned v; sscanf(h, "%2x", &v); out[n++] = (char) v; }
    out[n] = 0;
    return n;
}
static char text[1 << 17], hx[1 << 18];

static void dump_ctab(MOL_FMT_CTAB *c) {
    printf("C %d %d %d %d %d ", c->n_atoms, c->n_bonds, c->chiral_flag, c->n_stext_entries, c->n_property_lines);
    hexs(c->version_string);
    printf(" %d\n", c->v3000 != NULL);
}

int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        char kind;
        printf("#CASE %d\n", cs);
        scanf(" %c", &kind);
        if (kind == 'F') {
            char t; int len;
            scanf(" %c %d %s", &t, &len, hx);
            int n = unhex(hx, text);
            char *p = text;
            union { char s[400]; S_CHAR c; S_SHORT h; long l; double d; float f; } u;
            memset(&u, 0, sizeof(u));
            int ret = MolfileReadField(&u, len, t, &p);
            printf("%d %ld ", ret, (long) (p - text));
            switch (t) {
            case 'S': hexs(u.s); break;
            case 'C': printf("%d", u.c); break;
            case 'N': printf("%d", u.h); break;
            case 'L': printf("%ld", u.l); break;
            case 'D': pd(u.d); break;
            case 'F': pd((double) u.f); break;
            default: printf("-");
            }
            printf("\n");
            continue;
        }
        int bcoord;
        scanf("%d %s", &bcoord, hx);
        int n = unhex(hx, text);
        INCHI_IOSTREAM ios;
        inchi_ios_init(&ios, INCHI_IOS_TYPE_STRING, NULL);
        ios.s.pStr = (char *) malloc(n + 2);
        memcpy(ios.s.pStr, text, n + 1);
        ios.s.nUsedLength = n; ios.s.nAllocatedLength = n + 2; ios.s.nPtr = 0;
        char err[STR_ERR_LEN]; err[0] = 0;
        MOL_FMT_DATA *m = (MOL_FMT_DATA *) calloc(1, sizeof(MOL_FMT_DATA));
        MOL_FMT_CTAB *ct = &m->ctab;
        int retcode, prevcode, data_ended = 0;
        retcode = MolfileReadHeaderLines(&m->hdr, &ios, err);
        MOL_FMT_HEADER_BLOCK *h = &m->hdr;
        printf("H "); hexs(h->molname); printf(" "); hexs(h->line2); printf(" "); hexs(h->user_initls); printf(" "); hexs(h->prog_name);
        printf(" %d %d %d %d %d ", h->month, h->day, h->year, h->hour, h->minute); hexs(h->dim_code);
        printf(" %d ", h->scaling_factor1); pd(h->scaling_factor2); printf(" "); pd(h->energy); printf(" %ld ", h->internal_regno); hexs(h->comment);
        printf(" %d\n", retcode);
        printf("N %ld\n", MolfileExtractStrucNum(h));
        if (retcode) { retcode += 10; goto fin; }
        retcode = MolfileReadCountsLine(ct, &ios, err);
        dump_ctab(ct);
        if (retcode) { retcode += 20; goto fin; }
        if (ct->v3000) goto fin;
        {
            int na = inchi_max(ct->n_atoms, 1);
            ct->atoms = (MOL_FMT_ATOM *) calloc(na, sizeof(MOL_FMT_ATOM));
            if (bcoord) ct->coords = (MOL_COORD *) calloc(na, sizeof(MOL_COORD));
        }
        retcode = MolfileReadAtomsBlock(ct, &ios, retcode, err);
        if (retcode) { if (retcode < 0) { retcode = -retcode; data_ended = 1; } retcode += 30; }
        if (retcode < 30 && !data_ended) {
            ct->bonds = (MOL_FMT_BOND *) calloc(inchi_max(ct->n_bonds, 1), sizeof(MOL_FMT_BOND));
        }
        prevcode = retcode;
        if (!data_ended) {
            retcode = MolfileReadBondsBlock(ct, &ios, retcode, err);
            if (retcode) { if (retcode < 0) { retcode = -retcode; data_ended = 1; } retcode = prevcode ? prevcode : retcode + 40; }
        }
        prevcode = retcode;
        if (!data_ended) {
            retcode = MolfileReadSTextBlock(ct, &ios, retcode, err);
            if (retcode) retcode = prevcode ? prevcode : retcode + 50;
        }
        for (int i = 0; i < ct->n_atoms && i < 300; i++) {
            MOL_FMT_ATOM *a = &ct->atoms[i];
            printf("A %d ", i); pd(a->fx); printf(" "); pd(a->fy); printf(" "); pd(a->fz); printf(" "); hexs(a->symbol);
            printf(" %d %d %d %d %d\n", a->mass_difference, a->charge, a->radical, a->stereo_parity, a->valence);
            if (ct->coords) { printf("O %d ", i); hexs(ct->coords[i]); printf("\n"); }
        }
        if (ct->bonds) for (int i = 0; i < ct->n_bonds && i < 300; i++) {
            MOL_FMT_BOND *b = &ct->bonds[i];
            printf("B %d %d %d %d %d\n", i, b->atnum1, b->atnum2, b->bond_type, b->bond_stereo);
        }
        printf("HN %d\n", MolfileHasNoChemStruc(m));
        for (int f = 0; f < 2; f++) {
            double x0, y0, z0, xm, ym, zm, sc; int e2 = 0; char e2s[STR_ERR_LEN]; e2s[0] = 0;
            int d = MolfileGetXYZDimAndNormFactors(m, f, &x0, &y0, &z0, &xm, &ym, &zm, &sc, &e2, e2s);
            printf("X %d %d ", d, e2); pd(x0); printf(" "); pd(y0); printf(" "); pd(z0); printf(" "); pd(xm); printf(" "); pd(ym); printf(" "); pd(zm); printf(" "); pd(sc); printf(" "); hexs(e2s); printf("\n");
        }
    fin:
        printf("R %d %d ", retcode, data_ended); hexs(err); printf("\n");
        inchi_ios_free_str(&ios);
    }
    return 0;
}
"""

R = random.Random(2026)
SYMS = ["C", "N", "O", "H", "Cl", "CL", "Br", "Zz", "*", "D", "T", "Fe", "S", "P", "F", "", "A", "R#", "Xx"]


def fuzz(r, s, p):
    s = list(s)
    out = []
    for ch in s:
        x = r.random()
        if x < p * 0.4:
            out.append(r.choice("0123456789 -+.eEdx\tZ_"))
        elif x < p * 0.55:
            continue
        elif x < p * 0.7:
            out.append(ch)
            out.append(r.choice(" 1-.a"))
        else:
            out.append(ch)
    return "".join(out)


def num(r, w, d, v):
    s = f"{v:{w}.{d}f}" if d is not None else f"{v:{w}d}"
    return s


def mk_mol(r):
    na = r.choice([0, 1, 2, 3, 5, 8, 12, 20])
    nb = r.choice([0, 1, 2, 4, 7, 10]) if na > 1 else r.choice([0, 0, 1])
    lines = []
    name = r.choice(["Structure #%d" % r.randrange(1, 9999), "Structure #%d " % r.randrange(1, 99), "Structure #12x", "structure #7",
                     "mol", "", "  indented", "x" * r.randrange(1, 230), "name with trailing   "])
    lines.append(name)
    stamp = f"{r.choice(['  ', 'IN'])}{r.choice(['rdkit   ', 'InChI 1.', 'prog    '])}{r.randrange(0, 99):02d}{r.randrange(0, 99):02d}{r.randrange(0, 99):02d}{r.randrange(0, 99):02d}{r.randrange(0, 99):02d}{r.choice(['2D', '3D', '  ', '2d'])}"
    stamp += r.choice(["", " 1", "10", f"{r.random() * 3:10.5f}"]) + r.choice(["", f"{r.random():10.5f}", f"{-r.random() * 3000:12.5f}"]) + r.choice(["", "123456", "   12"])
    if r.random() < 0.15:
        stamp = "InChI 1.07  SDfile Output" if r.random() < 0.5 else "  InChI-1.07 SDfile Output 3D"
    lines.append(stamp)
    lines.append(r.choice(["", "comment", "c" * r.randrange(1, 120)]))
    cl = f"{na:3d}{nb:3d}  0  0{r.choice([0, 1, 0]):3d}{r.choice([0, 0, 0, 1, 2]):3d}  0  0  0  0{r.choice([999, 0, 3]):3d} {r.choice(['V2000', 'V2000', 'V2000', 'V3000', 'V200', '', 'v2000'])}"
    if r.random() < 0.1:
        cl = f"{r.randrange(-9, 120):3d}{r.randrange(-9, 120):3d}" + cl[6:]
    lines.append(cl)
    for _ in range(na):
        sym = r.choice(SYMS)
        a = (f"{r.uniform(-9, 9):10.4f}{r.uniform(-9, 9):10.4f}{r.choice([0, 0, 1]) * r.uniform(-3, 3):10.4f} {sym:<3}"
             f"{r.choice([0, 0, 0, 1, -1, 2, 127, -3, 4]):2d}{r.choice([0, 0, 0, 3, 5, 4, 7, 9, -2, 20]):3d}{r.choice([0, 1, 2, 3]):3d}"
             f"  0  0{r.choice([0, 0, 1, 2, 4, 15, 21]):3d}  0  0  0  0  0")
        lines.append(a)
    for _ in range(nb):
        a1, a2 = r.randrange(0, na + 2), r.randrange(0, na + 2)
        lines.append(f"{a1:3d}{a2:3d}{r.choice([1, 1, 2, 3, 4, 8, 0]):3d}{r.choice([0, 0, 1, 4, 6]):3d}  0  0")
    for _ in range(r.choice([0, 0, 0, 1, 2])):
        lines.append("x" * r.randrange(1, 30))
    lines.append("M  END")
    if r.random() < 0.3:
        lines.append("$$$$")
    # fuzz
    p = r.choice([0, 0, 0, 0.01, 0.03, 0.08])
    lines = [fuzz(r, l, p) if i > 2 else l for i, l in enumerate(lines)]
    # cut
    if r.random() < 0.15:
        lines = lines[:r.randrange(0, len(lines))]
    if r.random() < 0.07:
        k = r.randrange(0, len(lines)) if lines else 0
        lines.insert(k, "$$$$")
    nl = r.choice(["\n", "\n", "\n", "\r\n", "\r"])
    return nl.join(lines) + r.choice(["", nl, ""])


def mk_field(r):
    t = r.choice("SCNLDFJ")
    ln = r.choice([1, 2, 3, 3, 6, 8, 10, 12, 0, -1, 32, 33, 40, 200])
    nums = ["12", "-5", "  7", "7  ", " 1 2", "1x", "x1", "", "   ", "+3", "0x1F", "1e3", "1.5", "-0.0", "123456789012345678901234567890", "9223372036854775807",
            "-9223372036854775808", "99999999999", "inf", "-INF", "nan", "NaN(abc)", "infinity", "1e999", "1e-999", ".5", "5.", "1.5e", "1.5e+2", "0x1.8p1", "0x", " 12abc", "1,5", "3e-39", "-3e38", "3.5e38"]
    if t in "SJ":
        txt = r.choice(["abc", "  abc  ", "a b c d e f", "", " ", "abcdefghij", "\tx", "ZnZn", "C  "])
    else:
        txt = r.choice(nums)
    if r.random() < 0.5:
        txt = txt + r.choice(["", "abc", " 12", "  ", "99999"])
    if r.random() < 0.3:
        txt = " " * r.randrange(0, 4) + txt
    return f"F {t} {ln} {txt.encode().hex() or '-'}"


def main():
    cases = []
    for _ in range(2500):
        cases.append(mk_field(R))
    for _ in range(1400):
        cases.append(f"M {R.choice([0, 1])} {mk_mol(R).encode().hex() or '-'}")
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
    path = ROOT / "tests/data/inchi_molfmt1.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
