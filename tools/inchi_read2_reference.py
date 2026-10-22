"""Generates tests/data/inchi_read2.json from the InChI C code: the stereo
layer parsers of ichiread.c (`/b`, `/t`, `/m`, `/s`, their isotopic and
fixed-H variants) and `CopySegment`. The harness includes ichiread.c (its
parsers are static). Links the upstream InChI base library (see
inchi_markalt_reference.py).

    python tools/inchi_read2_reference.py
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
#include "ichiread.c"
static char hx[1 << 16], buf[1 << 16];
static int unhex(const char *h, char *out) {
    int n = 0;
    out[0] = 0;
    if (h[0] == '=') return 0;
    for (; h[0] && h[1]; h += 2) { unsigned v; sscanf(h, "%2x", &v); out[n++] = (char) v; }
    out[n] = 0;
    return n;
}
int CreateCompositeNormAtom(COMP_ATOM_DATA *a, INP_ATOM_DATA2 *b, int c) { return 0; }

static void dump_stereo(INChI_Stereo *s) {
    if (!s) { printf(" NULL"); return; }
    printf(" [%d %d %d %d %d%d%d%d%d", s->nNumberOfStereoCenters, s->nNumberOfStereoBonds, s->nCompInv2Abs, s->bTrivialInv,
           s->t_parity != NULL, s->nNumber != NULL, s->b_parity != NULL, s->nBondAtom1 != NULL, s->nBondAtom2 != NULL);
    if (s->t_parity && s->nNumber) { printf(" |"); for (int i = 0; i < s->nNumberOfStereoCenters; i++) printf(" %d:%d", s->nNumber[i], s->t_parity[i]); }
    if (s->b_parity && s->nBondAtom1 && s->nBondAtom2) { printf(" |"); for (int i = 0; i < s->nNumberOfStereoBonds; i++) printf(" %d-%d:%d", s->nBondAtom1[i], s->nBondAtom2[i], s->b_parity[i]); }
    printf("]");
}
static void dump_inchi(INChI *c) {
    printf("%d %d %d", c->nNumberOfAtoms, c->nNumberOfIsotopicAtoms, c->IsotopicAtom != NULL);
    for (int i = 0; i < c->nNumberOfIsotopicAtoms && c->IsotopicAtom; i++)
        printf(" (%d %d %d %d %d)", c->IsotopicAtom[i].nAtomNumber, c->IsotopicAtom[i].nIsoDifference, c->IsotopicAtom[i].nNum_H, c->IsotopicAtom[i].nNum_D, c->IsotopicAtom[i].nNum_T);
    dump_stereo(c->Stereo);
    dump_stereo(c->StereoIsotopic);
    printf("\n");
}
/* a stereo built from a compact description: nsc nsb inv triv tl bl */
static INChI_Stereo *mk_stereo(int nsc, int nsb, int inv, int triv, int tl, int bl) {
    INChI_Stereo *s = calloc(1, sizeof(*s));
    s->nNumberOfStereoCenters = nsc; s->nNumberOfStereoBonds = nsb; s->nCompInv2Abs = inv; s->bTrivialInv = triv;
    if (tl) { s->t_parity = calloc(nsc + 1, 1); s->nNumber = calloc(nsc + 1, sizeof(AT_NUMB));
              for (int i = 0; i < nsc; i++) { s->t_parity[i] = (i * 3 + 1) % 5; s->nNumber[i] = i * 7 + 2; } }
    if (bl) { s->b_parity = calloc(nsb + 1, 1); s->nBondAtom1 = calloc(nsb + 1, sizeof(AT_NUMB)); s->nBondAtom2 = calloc(nsb + 1, sizeof(AT_NUMB));
              for (int i = 0; i < nsb; i++) { s->b_parity[i] = i % 4 + 1; s->nBondAtom1[i] = i + 10; s->nBondAtom2[i] = i + 3; } }
    return s;
}
static INChI_Stereo *rd_stereo(void) {
    int pr; scanf("%d", &pr);
    if (!pr) return NULL;
    int a, b, c, d, e, f; scanf("%d %d %d %d %d %d", &a, &b, &c, &d, &e, &f);
    return mk_stereo(a, b, c, d, e, f);
}
static void rd_inchi(INChI *c) {
    int na, ni, ip; scanf("%d %d %d", &na, &ni, &ip);
    c->nNumberOfAtoms = na; c->nNumberOfIsotopicAtoms = ni;
    if (ip) { c->IsotopicAtom = calloc(ni + 1, sizeof(INChI_IsotopicAtom));
              for (int i = 0; i < ni; i++) { c->IsotopicAtom[i].nAtomNumber = i + 1; c->IsotopicAtom[i].nIsoDifference = i - 1; c->IsotopicAtom[i].nNum_H = i; c->IsotopicAtom[i].nNum_D = 2 * i; c->IsotopicAtom[i].nNum_T = 3; } }
    c->Stereo = rd_stereo();
    c->StereoIsotopic = rd_stereo();
}

int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int t;
        printf("#CASE %d\n", cs);
        scanf("%d", &t);
        if (t == 11) {
            int abc, nc[2], nseg;
            scanf("%d %d %d", &abc, &nc[0], &nc[1]);
            INChI *pp[2] = { NULL, NULL };
            for (int k = 0; k < 2; k++) if (nc[k] > 0) {
                pp[k] = calloc(nc[k] + 2, sizeof(INChI));
                for (int i = 0; i < nc[k]; i++) { int v; scanf("%d", &v); pp[k][i].nNumberOfAtoms = v; }
            }
            int sv[2][2] = { { 0, 0 }, { 0, 0 } };
            REM_PROTONS rp[2]; memset(rp, 0, sizeof(rp));
            scanf("%d", &nseg);
            for (int g = 0; g < nseg; g++) {
                int state; scanf("%d %s", &state, hx);
                unhex(hx, buf);
                int mob = (state >= IST_FIXED_H_LAYER_FORK) ? 0 : 1;
                int r;
                switch (state) {
                case IST_MOBILE_H_SP2: case IST_FIXED_H_SP2: case IST_MOBILE_H_ISO_SP2: case IST_FIXED_H_ISO_SP2:
                    r = ParseSegmentSp2(buf, mob, pp, nc, state, &abc); break;
                case IST_MOBILE_H_SP3: case IST_FIXED_H_SP3: case IST_MOBILE_H_ISO_SP3: case IST_FIXED_H_ISO_SP3:
                    r = ParseSegmentSp3(buf, mob, pp, nc, state, &abc); break;
                case IST_MOBILE_H_SP3_M: case IST_FIXED_H_SP3_M: case IST_MOBILE_H_ISO_SP3_M: case IST_FIXED_H_ISO_SP3_M:
                    r = ParseSegmentSp3m(buf, mob, pp, nc, state); break;
                case IST_MOBILE_H_ISO_ATOMS: case IST_FIXED_H_ISO_ATOMS:
                    r = ParseSegmentIsoAtoms(buf, mob, pp, nc, state, &abc); break;
                case IST_MOBILE_H_ISO_EXCH_H:
                    r = ParseSegmentIsoExchgH(buf, mob, rp, nc, state, &abc); break;
                case IST_FIXED_H_PERMUTATION: case IST_FIXED_H_ISO_PERMUTATION:
                    r = ParseSegmentPerm(buf, mob, pp, nc, state, &abc); break;
                default:
                    r = ParseSegmentSp3s(buf, mob, pp, sv, nc, state); break;
                }
                printf("R %d %d %d %d %d %d\n", state, r, abc, sv[0][0], sv[0][1], sv[1][0] * 10 + sv[1][1]);
                if (r < 0) break;
            }
            printf("P %d %d %d\n", rp[1].nNumRemovedIsotopicH[0], rp[1].nNumRemovedIsotopicH[1], rp[1].nNumRemovedIsotopicH[2]);
            for (int k = 0; k < 2; k++) for (int i = 0; i < nc[k]; i++) { printf("C%d.%d ", k, i); dump_inchi(&pp[k][i]); }
        } else if (t == 12) {
            int ty, bto, bfrom, hasfrom;
            INChI to, from;
            memset(&to, 0, sizeof(to)); memset(&from, 0, sizeof(from));
            scanf("%d %d %d %d", &ty, &bto, &bfrom, &hasfrom);
            rd_inchi(&to);
            if (hasfrom) rd_inchi(&from);
            int r = CopySegment(&to, hasfrom ? &from : NULL, ty, bto, bfrom);
            printf("R %d\nC ", r);
            dump_inchi(&to);
        }
    }
    return 0;
}
"""

R = random.Random(4242)


def hx(b):
    return b.hex() if b else "="


ST = {
    "mob": dict(b=5, t=6, m=7, s=8),
    "mobi": dict(b=12, t=13, m=14, s=15, i=10, h=11),
    "fix": dict(b=20, t=21, m=22, s=23, o=24),
    "fixi": dict(b=28, t=29, m=30, s=31, i=26, o=32),
}


def abcnum(n):
    ds = []
    while n:
        ds.append(n % 27)
        n //= 27
    ds.reverse()
    out = ""
    for i, d in enumerate(ds):
        if i == 0:
            out += chr(ord("A") + d - 1)
        else:
            out += "@" if d == 0 else chr(ord("a") + d - 1)
    return out


def gen_items(kind, n, abc):
    """items of one component with n atoms"""
    if n < 1:
        return ""
    par = "-+u?"
    if kind == "t":
        k = R.randrange(0, min(n, 6) + 1) if R.random() < 0.9 else R.randrange(1, 8)
        atoms = sorted(R.sample(range(1, n + 1), min(k, n)))
        if abc:
            return "".join(abcnum(a) + str(R.randrange(1, 5)) for a in atoms)
        return ",".join(f"{a}{R.choice(par)}" for a in atoms)
    pairs = sorted({(R.randrange(2, n + 1), 0) for _ in range(6)} if n > 1 else set())
    pairs = sorted({(a, R.randrange(1, a)) for a, _ in pairs})
    if R.random() < 0.7:
        pairs = pairs[: R.randrange(0, len(pairs) + 1)]
    if abc:
        return "".join(abcnum(a) + abcnum(b) + str(R.randrange(1, 5)) for a, b in pairs)
    return ",".join(f"{a}-{b}{R.choice(par)}" for a, b in pairs)


def gen_layer(kind, atoms, mode, abc, first_of=None):
    """atoms: list of atom counts; mode: mob/mobi/fix/fixi"""
    n = len(atoms)
    abbr = {"mob": "", "mobi": "me", "fix": "me", "fixi": "mMne"}[mode]
    groups = []
    i = 0
    while i < n:
        left = n - i
        r = R.random()
        if abbr and r < 0.35:
            cnt = R.randrange(1, left + 1)
            groups.append((R.choice(["", "", str(cnt)]) if cnt == 1 else str(cnt)) + R.choice(abbr))
            i += cnt
        elif r < 0.55 and left > 1:
            cnt = R.randrange(2, left + 1)
            groups.append(f"{cnt}*" + gen_items(kind, atoms[i], abc))
            i += cnt
        else:
            groups.append(gen_items(kind, atoms[i], abc))
            i += 1
    return kind + ";".join(groups)


def mutate(s):
    c = list(s)
    for _ in range(R.choice([0, 0, 0, 0, 1, 1, 2, 3])):
        k = R.randrange(0, len(c) + 1)
        op = R.choice(["del", "ins", "rep"])
        ch = R.choice(list("0123456789,-+u?;*mMneN.AbcZ@ "))
        if op == "del" and c:
            c.pop(min(k, len(c) - 1))
        elif op == "ins":
            c.insert(k, ch)
        elif c:
            c[min(k, len(c) - 1)] = ch
    return "".join(c)


def gen_iso_items(n, abc):
    k = R.randrange(0, min(n, 5) + 1)
    atoms = sorted(R.sample(range(1, n + 1), min(k, n)))
    out = []
    for a in atoms:
        if abc:
            iso = R.choice([0, 0, 1, -1, 2, -3, 12])
            cnt = "".join(c + (str(R.choice([2, 3, 10])) if R.random() < 0.3 else "") for c in R.sample("tdh", R.randrange(0, 4)))
            out.append(abcnum(a) + str(iso) + cnt)
        else:
            iso = R.choice(["", "", "+1", "-1", "+2", "-3", "+0", "-0", "+12"])
            cnt = "".join(c + (str(R.choice([2, 3, 10])) if R.random() < 0.3 else "") for c in R.sample("TDH", R.randrange(0, 4)))
            out.append(f"{a}{iso}{cnt}")
    return "".join(out) if abc else ",".join(out)


def gen_iso(atoms, mode, abc):
    n = len(atoms)
    abbr = "me" if mode == "fixi" else ""
    groups = []
    i = 0
    while i < n:
        left = n - i
        r = R.random()
        if abbr and r < 0.35:
            cnt = R.randrange(1, left + 1)
            groups.append((R.choice(["", str(cnt)]) if cnt == 1 else str(cnt)) + R.choice(abbr))
            i += cnt
        elif r < 0.55 and left > 1:
            cnt = R.randrange(2, left + 1)
            groups.append(f"{cnt}*" + gen_iso_items(atoms[i], abc))
            i += cnt
        else:
            groups.append(gen_iso_items(atoms[i], abc))
            i += 1
    return "i" + ";".join(groups)


def gen_exchg(abc):
    if abc:
        parts = [f"{R.choice([1, 2, 5, 12])}{c}" for c in "tdh" if R.random() < 0.5]
        return "h" + "".join(parts)
    return "h" + "".join(c + R.choice(["", "", "2", "7"]) for c in "TDH" if R.random() < 0.5)


def gen_perm(n, abc):
    comps = list(range(1, n + 1))
    R.shuffle(comps)
    cycles = []
    while len(comps) >= 2 and R.random() < 0.8:
        k = R.randrange(2, min(len(comps), 4) + 1)
        cycles.append(comps[:k])
        comps = comps[k:]
    if abc:
        return "o" + "".join("(" + "".join(abcnum(c) for c in cy) + ")" for cy in cycles)
    return "o" + "".join("(" + ",".join(map(str, cy)) + ")" for cy in cycles)


def gen_m(n, present=True):
    if not present:
        return R.choice(["", "x"])
    r = R.random()
    if r < 0.1:
        return "m"
    return "m" + "".join(R.choice("01.") for _ in range(n))


def seq_case():
    abc0 = R.choice([-1, -1, -1, 0, 1])
    n = R.choice([1, 1, 2, 3, 4, 6])
    atoms = [R.randrange(1, 14) for _ in range(n)]
    if R.random() < 0.1:
        atoms[R.randrange(n)] = R.choice([0, 1, 30000, 40000])
    ncs = [n, n]
    if R.random() < 0.15:
        ncs[R.randrange(2)] = 0
    atom_lists = [atoms if ncs[0] else [], atoms if ncs[1] else []]
    abc = {-1: R.random() < 0.5, 0: False, 1: True}[abc0]
    segs = []
    modes = []
    if R.random() < 0.9:
        modes.append(("mob", None))
    if R.random() < 0.5:
        modes.append(("mobi", None))
    if R.random() < 0.5:
        modes.append(("fix", None))
    if R.random() < 0.4:
        modes.append(("fixi", None))
    R.shuffle(modes) if R.random() < 0.1 else None
    for mode, _ in modes:
        tgt = 1 if mode.startswith("mob") else 0
        a = atom_lists[tgt] or atoms
        kinds = ["b", "t", "m", "s"]
        if mode == "mobi":
            kinds = ["i", "h"] + kinds
        elif mode == "fix":
            kinds = kinds + ["o"]
        elif mode == "fixi":
            kinds = ["i"] + kinds + ["o"]
        for kind in kinds:
            if R.random() < 0.2:
                continue
            if kind == "i":
                s = gen_iso(a, mode, abc)
            elif kind == "h":
                s = gen_exchg(abc)
            elif kind == "o":
                s = gen_perm(len(a), abc)
            elif kind in "bt":
                s = gen_layer(kind, a, mode, abc)
                if R.random() < 0.08:
                    s = kind
            elif kind == "m":
                s = gen_m(len(a), R.random() < 0.6)
            else:
                s = R.choice(["s", "s1", "s2", "s3", "s0", "s4", "s", "x", "s1"])
            segs.append((ST[mode][kind], mutate(s) if R.random() < 0.3 else s))
    if R.random() < 0.05:
        segs = segs[: R.randrange(0, len(segs) + 1)]
    flat = " ".join(f"{st} {hx(s.encode())}" for st, s in segs)
    return f"11 {abc0} {ncs[0]} {ncs[1]} " + " ".join(map(str, atom_lists[0] + atom_lists[1])) + f" {len(segs)} {flat}"


def stereo_desc():
    if R.random() < 0.25:
        return "0"
    nsc = R.choice([0, 0, 1, 2, 3, 5])
    nsb = R.choice([0, 0, 1, 2, 4])
    return f"1 {nsc} {nsb} {R.choice([0, 0, 1, -1, 9999])} {R.choice([0, 1, 1, 2])} {R.choice([0, 1, 1])} {R.choice([0, 1, 1])}"


def inchi_desc(to=False):
    na = R.choice([0, 3, 5, 12])
    ni = R.choice([0, 0, 1, 3])
    # C copies len + 1 isotopic atoms into an existing array and dumps
    # nNumberOfIsotopicAtoms of them: keep the target array consistent
    ip = int(ni > 0) if to else R.choice([0, 1, 1])
    return f"{na} {ni} {ip} {stereo_desc()} {stereo_desc()}"


def main():
    cases = [seq_case() for _ in range(26000)]
    for _ in range(4000):
        hasfrom = int(R.random() < 0.9)
        cases.append(f"12 {R.choice([0, 1, 2, 3, 4, 4])} {R.choice([0, 1])} {R.choice([-1, 0, 1])} {hasfrom} {inchi_desc(True)}" + (f" {inchi_desc()}" if hasfrom else ""))
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        kept, out = [], []
        for c in cases:
            try:
                r1 = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=1)
            except subprocess.TimeoutExpired:
                continue
            if r1.returncode != 0:
                continue
            kept.append(c)
            out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode(errors="replace").split("\n")[1:]
        print(len(cases) - len(kept), "cases dropped")
        cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_read2.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
