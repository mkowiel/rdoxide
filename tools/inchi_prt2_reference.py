"""Generates tests/data/inchi_prt2.json from the InChI C code: the string
construction and comparison helpers of ichiprt2.c (`MakeAbcNumber`,
`MakeDecNumber`, `MakeMult`, `MakeDelim`, `MakeEqStr`, `MakeCtString`,
`MakeHString`, `MakeTautString`, `MakeCRVString`, `MakeEquString`,
`MakeIsoAtomString`, `MakeIsoTautString`, `MakeIsoHString`, `MakeStereoString`,
`print_sequence_of_nums_compressing_ranges`, `Eql_INChI_Aux_Equ/Num`,
`Eql_INChI_Isotopic`, `bHasOrigInfo`, `EqlOrigInfo`, `bHasEquString`) on random
inputs. Links the upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_prt2_reference.py
"""
import copy
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("markalt", ROOT / "tools/inchi_markalt_reference.py")
ma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ma)
SRC = ma.SRC

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mode.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichi_bns.h"
#include "ichicano.h"
#include "ichicant.h"
#include "ichinorm.h"
#include "ichimake.h"
#include "ichi_io.h"
int CreateCompositeNormAtom() { return 0; }
int Eql_INChI_Aux_Equ( INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2 );
int Eql_INChI_Aux_Num( INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2 );
int Eql_INChI_Isotopic( INChI *i1, INChI *i2 );
int bHasEquString( AT_NUMB *LinearCT, int nLenCT );
static void zero_stack(void) { volatile char b[32768]; memset((void*)b, 0, sizeof(b)); }
static int rd(void) { int x; scanf("%d", &x); return x; }
static AT_NUMB *rd_at(int *n) { if (!rd()) { *n = 0; return NULL; } int l = rd(); *n = l; AT_NUMB *a = calloc(l + 8, sizeof(AT_NUMB)); for (int i = 0; i < l; i++) a[i] = rd(); return a; }
static S_CHAR *rd_sc(int *n) { if (!rd()) { *n = 0; return NULL; } int l = rd(); *n = l; S_CHAR *a = calloc(l + 8, 1); for (int i = 0; i < l; i++) a[i] = rd(); return a; }
static ORIG_INFO *rd_oi(int *n) { if (!rd()) { *n = 0; return NULL; } int l = rd(); *n = l; ORIG_INFO *a = calloc(l + 8, sizeof(ORIG_INFO)); for (int i = 0; i < l; i++) { a[i].cCharge = rd(); a[i].cRadical = rd(); a[i].cUnusualValence = rd(); } return a; }
static INChI_Aux *rd_aux(void) {
    if (!rd()) return NULL;
    INChI_Aux *a = calloc(1, sizeof(*a)); int n;
    a->nNumberOfAtoms = rd(); a->nNumberOfTGroups = rd(); a->bDeleted = rd(); a->bIsIsotopic = rd();
    a->nConstitEquNumbers = rd_at(&n); a->nConstitEquTGroupNumbers = rd_at(&n); a->nConstitEquIsotopicNumbers = rd_at(&n);
    a->nConstitEquIsotopicTGroupNumbers = rd_at(&n); a->nOrigAtNosInCanonOrd = rd_at(&n); a->nIsotopicOrigAtNosInCanonOrd = rd_at(&n);
    a->nOrigAtNosInCanonOrdInv = rd_at(&n); a->nIsotopicOrigAtNosInCanonOrdInv = rd_at(&n);
    a->OrigInfo = rd_oi(&n);
    return a;
}
static INChI *rd_iso_inchi(void) {
    if (!rd()) return NULL;
    INChI *c = calloc(1, sizeof(*c));
    c->bDeleted = rd();
    int l = c->nNumberOfIsotopicAtoms = rd();
    if (rd()) { c->IsotopicAtom = calloc(l + 4, sizeof(INChI_IsotopicAtom)); for (int i = 0; i < l; i++) { c->IsotopicAtom[i].nAtomNumber = rd(); c->IsotopicAtom[i].nIsoDifference = rd();
        c->IsotopicAtom[i].nNum_H = rd(); c->IsotopicAtom[i].nNum_D = rd(); c->IsotopicAtom[i].nNum_T = rd(); } }
    l = c->nNumberOfIsotopicTGroups = rd();
    if (rd()) { c->IsotopicTGroup = calloc(l + 4, sizeof(INChI_IsotopicTGroup)); for (int i = 0; i < l; i++) { c->IsotopicTGroup[i].nTGroupNumber = rd();
        c->IsotopicTGroup[i].nNum_H = rd(); c->IsotopicTGroup[i].nNum_D = rd(); c->IsotopicTGroup[i].nNum_T = rd(); } }
    return c;
}
int main(void) {
    int ncases;
    char cmd[32], d1[16];
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        INCHI_IOS_STRING sb; inchi_strbuf_init(&sb, 16, 16);
        int ovf = 0, ret = 0, n;
        scanf("%s", cmd);
        printf("#CASE %d\n", c);
        if (!strcmp(cmd, "abc") || !strcmp(cmd, "dec")) {
            int hasd = rd(); scanf("%s", d1); int v = rd(), bl = rd();
            char buf[64]; memset(buf, 0x7f, sizeof(buf));
            ret = !strcmp(cmd, "abc") ? MakeAbcNumber(buf, bl, hasd ? d1 : NULL, v) : MakeDecNumber(buf, bl, hasd ? d1 : NULL, v);
            printf("%d [%s]\n", ret, ret >= 0 ? buf : "");
            continue;
        } else if (!strcmp(cmd, "mult")) {
            int m = rd(); scanf("%s", d1); int mode = rd(); ovf = rd();
            ret = MakeMult(m, d1, &sb, mode, &ovf);
        } else if (!strcmp(cmd, "delim")) {
            int has = rd(); scanf("%s", d1); ovf = rd();
            ret = MakeDelim(has ? d1 : NULL, &sb, &ovf);
        } else if (!strcmp(cmd, "eqstr")) {
            int has = rd(); scanf("%s", d1); int m = rd(); ovf = rd();
            ret = MakeEqStr(has ? d1 : NULL, m, &sb, &ovf);
        } else if (!strcmp(cmd, "ct")) {
            int len = rd(); AT_NUMB *ct = calloc(len + 8, sizeof(AT_NUMB)); for (int i = 0; i < len; i++) ct[i] = rd();
            int na = rd(); int hasH = rd(); S_CHAR *h = calloc(na + 8, 1); for (int i = 0; i < na; i++) h[i] = rd();
            int ad = rd(), mode = rd();
            CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
            zero_stack();
            ret = MakeCtString(&cg, ct, len, ad, hasH ? h : NULL, na, &sb, mode, &ovf);
        } else if (!strcmp(cmd, "h")) {
            int ad = rd(); int len = rd(); S_CHAR *a = calloc(len + 8, 1); for (int i = 0; i < len; i++) a[i] = rd(); int mode = rd();
            ret = MakeHString(ad, a, len, &sb, mode, &ovf);
        } else if (!strcmp(cmd, "taut")) {
            int len = rd(); AT_NUMB *a = calloc(len + 8, sizeof(AT_NUMB)); for (int i = 0; i < len; i++) a[i] = rd(); int ad = rd(), mode = rd();
            ret = MakeTautString(a, len, ad, &sb, mode, &ovf);
        } else if (!strcmp(cmd, "crv")) {
            int len = rd(); ORIG_INFO *a = calloc(len + 8, sizeof(ORIG_INFO)); for (int i = 0; i < len; i++) { a[i].cCharge = rd(); a[i].cRadical = rd(); a[i].cUnusualValence = rd(); }
            int ad = rd(), mode = rd();
            ret = MakeCRVString(a, len, ad, &sb, mode, &ovf);
        } else if (!strcmp(cmd, "equ")) {
            int len = rd(); AT_NUMB *a = calloc(len + 8, sizeof(AT_NUMB)); for (int i = 0; i < len; i++) a[i] = rd(); int ad = rd(), mode = rd();
            ret = MakeEquString(a, len, ad, &sb, mode, &ovf);
        } else if (!strcmp(cmd, "isoatom")) {
            int len = rd(); INChI_IsotopicAtom *a = calloc(len + 8, sizeof(*a));
            for (int i = 0; i < len; i++) { a[i].nAtomNumber = rd(); a[i].nIsoDifference = rd(); a[i].nNum_H = rd(); a[i].nNum_D = rd(); a[i].nNum_T = rd(); }
            int mode = rd();
            ret = MakeIsoAtomString(a, len, &sb, mode, &ovf);
        } else if (!strcmp(cmd, "isotaut")) {
            int len = rd(); INChI_IsotopicTGroup *a = calloc(len + 8, sizeof(*a));
            for (int i = 0; i < len; i++) { a[i].nTGroupNumber = rd(); a[i].nNum_H = rd(); a[i].nNum_D = rd(); a[i].nNum_T = rd(); }
            int mode = rd();
            ret = MakeIsoTautString(a, len, &sb, mode, &ovf);
        } else if (!strcmp(cmd, "isoh")) {
            int h[3]; h[0] = rd(); h[1] = rd(); h[2] = rd(); int mode = rd();
            zero_stack();
            ret = MakeIsoHString(h, &sb, mode, &ovf);
        } else if (!strcmp(cmd, "stereo")) {
            int len = rd(); int n1, n2, n3;
            AT_NUMB *a1 = rd_at(&n1), *a2 = rd_at(&n2); S_CHAR *p = rd_sc(&n3); int ad = rd(), mode = rd();
            ret = MakeStereoString(a1, a2, p, ad, len, &sb, mode, &ovf);
        } else if (!strcmp(cmd, "seq")) {
            int len = rd(); int *a = calloc(len + 1, sizeof(int)); for (int i = 0; i < len; i++) a[i] = rd();
            print_sequence_of_nums_compressing_ranges(len, a, &sb);
        } else if (!strcmp(cmd, "eqlaux")) {
            INChI_Aux *a1 = rd_aux(), *a2 = rd_aux(); int e1 = rd(), e2 = rd();
            printf("%d %d %d\n", Eql_INChI_Aux_Equ(a1, e1, a2, e2), Eql_INChI_Aux_Num(a1, e1 & 3, a2, e2 & 3), EqlOrigInfo(a1, a2));
            printf("%d\n", a1 && a1->OrigInfo ? bHasOrigInfo(a1->OrigInfo, a1->nNumberOfAtoms) : 0);
            continue;
        } else if (!strcmp(cmd, "eqliso")) {
            INChI *i1 = rd_iso_inchi(), *i2 = rd_iso_inchi();
            printf("%d\n", Eql_INChI_Isotopic(i1, i2));
            continue;
        } else if (!strcmp(cmd, "hasequ")) {
            int len = rd(); AT_NUMB *a = calloc(len + 8, sizeof(AT_NUMB)); for (int i = 0; i < len; i++) a[i] = rd();
            printf("%d\n", bHasEquString(a, len));
            continue;
        }
        printf("%d %d [%s]\n", ret, ovf, sb.pStr ? sb.pStr : "");
    }
    return 0;
}
"""

R = random.Random(2024)
MODES = [0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 34, 42, 43, 35]


def ints(*a):
    return " ".join(str(x) for x in a)


def gen_ct(r):
    n = r.randrange(1, 9)
    edges = set()
    for v in range(1, n):
        if r.random() < 0.9:
            edges.add((r.randrange(0, v), v))
    for _ in range(r.choice([0, 0, 1, 2])):
        if n > 1:
            a, b = sorted(r.sample(range(n), 2))
            edges.add((a, b))
    nb = {v: sorted(a for a, b in edges if b == v) for v in range(n)}
    ct = []
    for v in range(n):
        ct.append(v + 1)
        ct += [a + 1 for a in nb[v]]
    return n, ct


def optarr(r, vals, p=0.85):
    return "0" if r.random() > p else f"1 {len(vals)} " + " ".join(map(str, vals))


def gen_equ(r, n):
    # canonical equivalence numbers: each class is represented by its smallest member
    cls = [r.randrange(1, n + 1) for _ in range(n)]
    out = []
    for i in range(n):
        out.append(min(j for j in range(n) if cls[j] == cls[i]) + 1)
    return out


def mk_aux(r):
    n = r.randrange(1, 7)
    t = r.randrange(0, 4)
    d = {"n": n, "t": t, "del": r.choice([0, 0, 0, 1]), "iso": r.choice([0, 1, 1]),
         "arrs": [gen_equ(r, n), gen_equ(r, max(t, 1)), gen_equ(r, n), gen_equ(r, max(t, 1)),
                  r.sample(range(1, n + 1), n), r.sample(range(1, n + 1), n), r.sample(range(1, n + 1), n), r.sample(range(1, n + 1), n)],
         "present": [r.random() < 0.9 for _ in range(8)],
         "oi": [[r.choice([0, 0, 1, -1]), r.choice([0, 0, 0, 1, 2, 3]), r.choice([0, 0, 0, 3, 4])] for _ in range(n)],
         "oip": r.random() < 0.8}
    return d


def tok_aux(a):
    if a is None:
        return "0"
    t = f"1 {a['n']} {a['t']} {a['del']} {a['iso']}"
    for arr, p in zip(a["arrs"], a["present"]):
        t += " " + (f"1 {len(arr)} " + " ".join(map(str, arr)) if p else "0")
    t += " " + (f"1 {len(a['oi'])} " + " ".join(" ".join(map(str, x)) for x in a["oi"]) if a["oip"] else "0")
    return t


def mut_aux(r, a):
    a = copy.deepcopy(a)
    k = r.random()
    if k < 0.35:
        return a
    if k < 0.5:
        a["del"] ^= 1
    elif k < 0.6:
        a["iso"] ^= 1
    elif k < 0.8:
        i = r.randrange(8)
        if a["arrs"][i]:
            a["arrs"][i][r.randrange(len(a["arrs"][i]))] = r.randrange(1, 5)
    elif k < 0.9:
        a["oi"][r.randrange(len(a["oi"]))][0] += 1
    else:
        a["present"][r.randrange(8)] ^= True
    return a


def mk_isoinchi(r):
    na = r.randrange(0, 3)
    nt = r.randrange(0, 3)
    return {"del": r.choice([0, 0, 0, 1]),
            "atoms": [[r.randrange(1, 5), r.choice([0, 1, -1, 2]), r.choice([0, 1]), r.choice([0, 1]), r.choice([0, 1])] for _ in range(na)],
            "tg": [[r.randrange(1, 4), r.choice([0, 1]), r.choice([0, 1]), r.choice([0, 1])] for _ in range(nt)],
            "pa": r.random() < 0.9, "pt": r.random() < 0.9}


def tok_isoinchi(d):
    if d is None:
        return "0"
    t = f"1 {d['del']} {len(d['atoms'])} " + ("1 " + " ".join(" ".join(map(str, a)) for a in d["atoms"]) if d["pa"] else "0")
    t += f" {len(d['tg'])} " + ("1 " + " ".join(" ".join(map(str, a)) for a in d["tg"]) if d["pt"] else "0")
    return t


def mut_iso(r, d):
    d = copy.deepcopy(d)
    k = r.random()
    if k < 0.4:
        return d
    if k < 0.6 and d["atoms"]:
        d["atoms"][0][1] += 1
    elif k < 0.75 and d["tg"]:
        d["tg"][0][1] ^= 1
    elif k < 0.85:
        d["del"] ^= 1
    else:
        d["pa"] ^= True
    return d


def gen_case(r):
    kind = r.choice(["abc", "dec", "mult", "delim", "eqstr", "ct", "ct", "ct", "h", "h", "taut", "crv", "equ", "isoatom", "isotaut",
                     "isoh", "stereo", "seq", "eqlaux", "eqlaux", "eqliso", "hasequ"])
    if kind in ("abc", "dec"):
        v = r.choice([0, 1, 5, 26, 27, 28, 100, 728, 19683, 32766, -1, -27, -400, 123456789, -2147483647])
        hd = r.randrange(2)
        return f"{kind} {hd} {r.choice([',', '-', '(', '+', ''])or ','} {v} {r.choice([2, 3, 4, 6, 8, 12, 20, 64])}"
    if kind == "mult":
        return f"mult {r.choice([0, 1, 2, 5, 30, 100])} {r.choice([',', ')', 'x'])} {r.choice([0, 2])} {r.choice([0, 0, 1])}"
    if kind == "delim":
        return f"delim {r.randrange(2)} {r.choice([',', '/', ';'])} {r.choice([0, 0, 1])}"
    if kind == "eqstr":
        return f"eqstr {r.randrange(2)} {r.choice([',', '*'])} {r.choice([1, 2, 7, 35])} {r.choice([0, 0, 1])}"
    if kind == "ct":
        n, ct = gen_ct(r)
        hh = [r.choice([0, 0, 1, 2, 3]) for _ in range(n)]
        return f"ct {len(ct)} {ints(*ct)} {n} {r.randrange(2)} {ints(*hh)} {r.randrange(2)} {r.choice(MODES)}"
    if kind == "h":
        n = r.randrange(1, 9)
        a = [r.choice([0, 0, 1, 1, 2, 3, -1, -2, 20, -5, 17]) for _ in range(n)]
        return f"h {r.randrange(2)} {n} {ints(*a)} {r.choice([0, 2, 16, 18])}"
    if kind == "taut":
        groups = []
        for _ in range(r.choice([1, 1, 2, 3])):
            eps = sorted(r.sample(range(1, 9), r.randrange(1, 5)))
            groups.append([len(eps) + 2, r.choice([0, 1, 2]), r.choice([0, 0, 1])] + eps)
        a = [len(groups)] + [x for g in groups for x in g]
        if r.random() < 0.1:
            a = [0]
        return f"taut {len(a)} {ints(*a)} {r.randrange(2)} {r.choice([0, 2])}"
    if kind == "crv":
        n = r.randrange(1, 6)
        a = [[r.choice([0, 0, 1, -1, 2]), r.choice([0, 0, 1, 2, 3]), r.choice([0, 0, 3, 4])] for _ in range(n)]
        return f"crv {n} {' '.join(' '.join(map(str, x)) for x in a)} {r.randrange(2)} {r.choice([0, 2])}"
    if kind == "equ":
        n = r.randrange(1, 8)
        a = gen_equ(r, n)
        return f"equ {n} {ints(*a)} {r.randrange(2)} {r.choice([0, 2])}"
    if kind == "isoatom":
        n = r.randrange(0, 4)
        a = [[r.randrange(1, 40), r.choice([0, 1, -1, 2, -3, 5]), r.choice([0, 1, 2]), r.choice([0, 1, 12]), r.choice([0, 1, 3])] for _ in range(n)]
        return f"isoatom {n} {' '.join(' '.join(map(str, x)) for x in a)} {r.choice([0, 2])}"
    if kind == "isotaut":
        n = r.randrange(0, 4)
        a = [[r.randrange(1, 6), r.choice([0, 1, 2]), r.choice([0, 1]), r.choice([0, 1, 4])] for _ in range(n)]
        return f"isotaut {n} {' '.join(' '.join(map(str, x)) for x in a)} {r.choice([0, 2])}"
    if kind == "isoh":
        h = [r.choice([0, 1, 2]), r.choice([0, 1, 3]), r.choice([0, 1, 5])]
        if not any(h):  # all zero prints the uninitialised szValue
            h[r.randrange(3)] = 1
        return f"isoh {ints(*h)} {r.choice([0, 2])}"
    if kind == "stereo":
        n = r.randrange(0, 6)
        a1 = [r.randrange(1, 20) for _ in range(n)]
        a2 = [r.randrange(1, 20) for _ in range(n)]
        p = [r.choice([1, 2, 3, 4, 0, 7, -1]) for _ in range(n)]
        def opt(x):
            return "0" if r.random() < 0.2 else f"1 {n} {ints(*x)}"
        parts = [opt(a1), opt(a2), opt(p)]
        if all(x == "0" for x in parts):  # nothing written: the C code prints the uninitialised szValue
            parts[r.randrange(3)] = f"1 {n} {ints(*a1)}"
        return f"stereo {n} {parts[0]} {parts[1]} {parts[2]} {r.randrange(2)} {r.choice([0, 2])}"
    if kind == "seq":
        n = r.randrange(1, 12)
        a = []
        x = r.randrange(1, 5)
        for _ in range(n):
            a.append(x)
            x += r.choice([1, 1, 1, 2, 3])
        return f"seq {n} {ints(*a)}"
    if kind == "eqlaux":
        a1 = mk_aux(r) if r.random() > 0.05 else None
        a2 = (mut_aux(r, a1) if a1 else None) if r.random() > 0.05 else None
        return f"eqlaux {tok_aux(a1)} {tok_aux(a2)} {r.choice([0, 1, 2, 3])} {r.choice([0, 1, 2, 3])}"
    if kind == "eqliso":
        d1 = mk_isoinchi(r) if r.random() > 0.05 else None
        d2 = (mut_iso(r, d1) if d1 else None) if r.random() > 0.05 else None
        return f"eqliso {tok_isoinchi(d1)} {tok_isoinchi(d2)}"
    n = r.randrange(0, 8)
    return f"hasequ {n} {ints(*gen_equ(r, n)) if n else ''}"


def main():
    cases = [gen_case(R) for _ in range(9000)]
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
                    out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode().split("\n")[1:]
            print(len(cases) - len(kept), "cases dropped")
            cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases), (len(starts), len(cases))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    # MakeCtStringNew prints an unterminated szValue when the first DFS item has no number but H
    data = [x for x in data if all(ch.isprintable() and ord(ch) < 127 for l in x["out"] for ch in l)]
    path = ROOT / "tests/data/inchi_prt2.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
