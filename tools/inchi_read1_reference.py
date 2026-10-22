"""Generates tests/data/inchi_read1.json from the InChI C code: the first
parts of the InChI string reader of ichiread.c (error names, `/q` and `/p`
layers, `ParseSegmentReadDelimitedNumbers`, `nFillOutProtonMobileH`,
`insertions_sort_AT_NUMB`, `nGetInChISegment`). The harness includes
ichiread.c (its parsers are static). Links the upstream InChI base library
(see inchi_markalt_reference.py).

    python tools/inchi_read1_reference.py
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
static void hexs(const char *s, int n) { for (int i = 0; i < n; i++) printf("%02x", (unsigned char) s[i]); if (!n) printf("="); }
int CreateCompositeNormAtom(COMP_ATOM_DATA *a, INP_ATOM_DATA2 *b, int c) { return 0; }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int t;
        printf("#CASE %d\n", cs);
        scanf("%d", &t);
        if (t == 1) {
            int mob, n0, n1; INChI *pp[2];
            scanf("%d %d %d %s", &mob, &n0, &n1, hx);
            unhex(hx, buf);
            int nn = n0 > n1 ? n0 : n1;
            pp[0] = calloc(nn + 2, sizeof(INChI)); pp[1] = calloc(nn + 2, sizeof(INChI));
            for (int i = 0; i < n1; i++) { int v; scanf("%d", &v); pp[1][i].nTotalCharge = v; }
            int nc[2] = { n0, n1 };
            int r = ParseSegmentCharge(buf, mob, pp, nc);
            printf("R %d\nQ", r);
            for (int k = 0; k < 2; k++) { printf(" |"); for (int i = 0; i < nn; i++) printf(" %d", pp[k][i].nTotalCharge); }
            printf("\n");
        } else if (t == 2) {
            int mob; REM_PROTONS rp[2];
            scanf("%d %s", &mob, hx);
            unhex(hx, buf);
            memset(rp, 0, sizeof(rp)); rp[mob].nNumRemovedProtons = -777;
            int r = ParseSegmentProtons(buf, mob, rp, NULL);
            printf("R %d %d\n", r, rp[mob].nNumRemovedProtons);
        } else if (t == 3) {
            int from, pend; int ret; char d, s; INT_ARRAY a; memset(&a, 0, sizeof(a)); IntArray_Alloc(&a, 4);
            int dd, ss;
            scanf("%s %d %d %d %d", hx, &from, &pend, &dd, &ss);
            int n = unhex(hx, buf);
            d = dd; s = ss;
            const char *p = ParseSegmentReadDelimitedNumbers(buf + from, buf + pend, &a, d, s, &ret);
            printf("R %d %ld", ret, p ? (long) (p - buf) : -1L);
            for (int i = 0; i < a.used; i++) printf(" %d", a.item[i]);
            printf("\n");
        } else if (t == 4) {
            int used; char tok[64], tokh[128];
            scanf("%s %d %s", hx, &used, tokh);
            int n = unhex(hx, buf);
            unhex(tokh, tok);
            INCHI_IOSTREAM ios; inchi_ios_init(&ios, INCHI_IOS_TYPE_STRING, NULL);
            ios.s.pStr = buf; ios.s.nUsedLength = used > n ? n : used; ios.s.nPtr = 0;
            SEGM_LINE line; memset(&line, 0, sizeof(line));
            for (int k = 0; k < 40; k++) {
                int c = nGetInChISegment(&ios, &line, tok);
                printf("S %d %d ", c, line.len); hexs(line.str ? line.str : "", line.len); printf("\n");
                if (c <= 0) break;
            }
        } else if (t == 5) {
            int n; AT_NUMB v[40];
            scanf("%d", &n);
            for (int i = 0; i < n; i++) { int x; scanf("%d", &x); v[i] = x; }
            int r = insertions_sort_AT_NUMB(v, n);
            printf("R %d", r);
            for (int i = 0; i < n; i++) printf(" %d", v[i]);
            printf("\n");
        } else if (t == 6) {
            int st; char m[256];
            scanf("%d", &st);
            getInchiStateReadErr(st, m);
            printf("M %s\n%s\n", m, getInchiErrName(st));
        } else if (t == 8) {
            char h1[4096], h2[4096], b1[2048], b2[2048];
            int state, first;
            scanf("%s %s %d %d", h1, h2, &state, &first);
            unhex(h1, b1); unhex(h2, b2);
            INChI *pp[2] = { NULL, NULL }; int pn[2] = { 0, 0 };
            for (int k = 0; k < 2; k++) {
                int mob = (k == 0) == (first == 1) ? 1 : 0;
                int na = state;
                int r = ParseSegmentFormula(k == 0 ? b1 : b2, mob, pp, pn, &na);
                printf("R%d %d %d %d %d\n", mob, r, pn[mob], na, pn[1 - mob]);
                if (r > 0 || r == 0) for (int i = 0; i < pn[mob] && pp[mob]; i++) {
                    INChI *c = &pp[mob][i];
                    printf("C %d %s %d ", i, c->szHillFormula ? c->szHillFormula : "NULL", c->nNumberOfAtoms);
                    if (c->nAtom) for (int j = 0; j <= c->nNumberOfAtoms; j++) printf("%d,", c->nAtom[j]); else printf("NULL");
                    printf(" %d %d\n", c->nTotalCharge, c->bDeleted);
                }
                if (r < 0 && r != RI_ERR_SYNTAX && r != RI_ERR_PROGR) break;
            }
        } else if (t == 7) {
            INChI c; memset(&c, 0, sizeof(c));
            int r = nFillOutProtonMobileH(&c);
            printf("R %d %d %s %d %d %d %d %d %d %d\n", r, c.bDeleted, c.szHillFormula, c.nNumberOfAtoms, c.nAtom[0], c.nTotalCharge, c.nConnTable[0], c.lenConnTable, c.nTautomer[0], c.lenTautomer);
        }
    }
    return 0;
}
"""

R = random.Random(777)


def hx(b):
    return b.hex() if b else "="


def charge_str():
    parts = []
    for _ in range(R.randrange(0, 6)):
        k = R.random()
        if k < 0.3:
            parts.append(R.choice(["+", "-"]) + str(R.randrange(0, 5)))
        elif k < 0.5:
            parts.append(f"{R.randrange(1, 4)}*{R.choice('+-')}{R.randrange(0, 4)}")
        elif k < 0.7:
            parts.append(f"{R.choice(['', '2', '3', '0'])}{R.choice('mnMNe')}")
        elif k < 0.8:
            parts.append("")
        elif k < 0.9:
            parts.append(R.choice(["+300", "-999", "x", "++1", "+", "2*", "*+1", " 2", "-0", "+0"]))
        else:
            parts.append(R.choice(["+", "-"]) + str(R.randrange(0, 300)))
    s = ";".join(parts)
    if R.random() < 0.05:
        s = ""
    return "q" + s


def main():
    cases = []
    for _ in range(4000):
        mob = R.choice([0, 1])
        n1 = R.randrange(0, 6)
        n0 = R.choice([n1, n1, R.randrange(0, 6)])
        init = [R.choice([0, 0, 1, -1, 2, 5, 9999]) for _ in range(n1)]
        cases.append(f"1 {mob} {n0} {n1} {hx(charge_str().encode())} " + " ".join(map(str, init)))
    for _ in range(1500):
        s = "p" + R.choice(["+", "-", "", "x"]) + R.choice(["1", "2", "0", "12", "", "100000", "7;", "3;4"]) + R.choice(["", ";", "a"])
        if R.random() < 0.1:
            s = R.choice(["", "q+1", "pp", "p+"])
        cases.append(f"2 {R.choice([0, 1])} {hx(s.encode())}")
    for _ in range(1500):
        nums = ",".join(str(R.choice([R.randrange(0, 20), R.randrange(0, 70000), -3])) for _ in range(R.randrange(0, 8)))
        s = R.choice(["-", "c", "/", " "]) + nums + R.choice([";", "/", ",", "x", ",/", ";", "/", ";"]) + R.choice(["", "tail"])
        b = s.encode()
        cases.append(f"3 {hx(b)} 0 {R.choice([len(b), len(b) - 1, len(b) + 1, 3, 0]) if b else 0} {ord(',')} {ord(R.choice([';', '/', '-']))}")
    for _ in range(2500):
        parts = []
        for _ in range(R.randrange(0, 6)):
            parts.append(R.choice(["InChI=1S", "C2H6O", "c1-2-3", "h3H,2H2", "", "q+1", "é", "ab\tcd", "x\r\ny"]))
        b = "/".join(parts).encode()
        if R.random() < 0.1:
            b += b"\x00" + b"zz"
        tok = R.choice(["/\n\r\t\\", "/", "\n", ""]).encode()
        cases.append(f"4 {hx(b)} {R.choice([len(b), len(b), len(b) + 5, max(len(b) - 3, 0)])} {hx(tok)}")
    for _ in range(600):
        n = R.randrange(0, 12)
        cases.append(f"5 {n} " + " ".join(str(R.randrange(0, 40)) for _ in range(n)))
    for st in list(range(-3, 40)) + [100, 105, 136, 137, 99, -1, -2, -3, -4, -9]:
        cases.append(f"6 {st}")
    els = ["C", "H", "O", "N", "Cl", "Br", "Na", "Fe", "Xx", "c", "Zz", "S", "P", "D", "Ca", "Uuo"]

    def formula():
        parts = []
        for _ in range(R.randrange(0, 5)):
            f = "".join(R.choice(els) + R.choice(["", "", "2", "3", "10", "0", "40000"]) for _ in range(R.randrange(0, 5)))
            parts.append(R.choice(["", "", "2", "3", "0"]) + f)
        z = ".".join(parts)
        if R.random() < 0.2:
            z += R.choice(["/q+1", "q+1", "/c1-2", "p+1", "x", "/h", ".", "(", "-"])
        return z

    for _ in range(3000):
        a, b = formula(), formula()
        if R.random() < 0.5:
            b = a if R.random() < 0.5 else b
        cases.append(f"8 {hx(a.encode())} {hx(b.encode())} {R.choice([0, 0, 17, 1, 0])} {R.choice([0, 1])}")
    cases.append("7")
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
                continue  # the C code loops forever
            if r1.returncode != 0:
                continue
            kept.append(c)
            out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode(errors="replace").split("\n")[1:]
        print(len(cases) - len(kept), "cases dropped")
        cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_read1.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
