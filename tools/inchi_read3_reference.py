"""Generates tests/data/inchi_read3.json from the InChI C code: the segment
state machine `ReadInChILine` of ichiread.c over whole InChI strings (the
strings of the official RDKit corpus, mutated). The harness includes
ichiread.c (its parsers are static). Links the upstream InChI base library
(see inchi_markalt_reference.py).

    python tools/inchi_read3_reference.py
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


r2 = load("read2", "tools/inchi_read2_reference.py")
ma = r2.ma
SRC = ma.SRC

MAIN = r"""
static void dump_full(INChI *c) {
    printf("%d %d %d %lu %s |", c->nNumberOfAtoms, c->bDeleted, c->nTotalCharge, (unsigned long) c->nFlags, c->szHillFormula ? c->szHillFormula : "NULL");
    if (c->nAtom) for (int j = 0; j <= c->nNumberOfAtoms; j++) printf("%d,", c->nAtom[j]); else printf("NULL");
    printf(" | ");
    if (c->nNum_H) for (int j = 0; j <= c->nNumberOfAtoms; j++) printf("%d,", c->nNum_H[j]); else printf("NULL");
    printf(" | ");
    if (c->nNum_H_fixed) for (int j = 0; j <= c->nNumberOfAtoms; j++) printf("%d,", c->nNum_H_fixed[j]); else printf("NULL");
    printf(" | %d |", c->lenTautomer);
    if (c->nTautomer) for (int j = 0; j <= c->lenTautomer; j++) printf(" %d", c->nTautomer[j]); else printf(" NULL");
    printf(" | %d |", c->lenConnTable);
    if (c->nConnTable) for (int j = 0; j <= c->lenConnTable; j++) printf(" %d", c->nConnTable[j]); else printf(" NULL");
    printf(" | ");
    dump_inchi(c);
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int t;
        printf("#CASE %d\n", cs);
        scanf("%d", &t);
        if (t == 13) {
            scanf("%s", hx);
            int n = unhex(hx, buf);
            INCHI_IOSTREAM ios; inchi_ios_init(&ios, INCHI_IOS_TYPE_STRING, NULL);
            ios.s.pStr = buf; ios.s.nUsedLength = n; ios.s.nPtr = 0;
            SEGM_LINE line; memset(&line, 0, sizeof(line));
            char *pstr = NULL; int state = -77;
            INChI *pp[2][2] = { { NULL, NULL }, { NULL, NULL } };
            int nc[2][2]; memset(nc, 0, sizeof(nc));
            REM_PROTONS rp[2][2]; memset(rp, 0, sizeof(rp));
            int s[2][2][2]; memset(s, 0, sizeof(s));
            int stdf = -5, hso = -5; unsigned char bits = 0;
            OAD_Polymer *poly = NULL; OAD_V3000 *v3 = NULL;
            int r = ReadInChILine(&ios, &line, &pstr, &state, pp, nc, rp, s, &stdf, &hso, &bits, 1, &poly, &v3);
            printf("R %d %d %d %d %d %d\n", r, state, stdf, hso, bits, line.len);
            printf("N %d %d %d %d\n", nc[0][0], nc[0][1], nc[1][0], nc[1][1]);
            printf("S %d %d %d %d %d %d %d %d\n", s[0][0][0], s[0][0][1], s[0][1][0], s[0][1][1], s[1][0][0], s[1][0][1], s[1][1][0], s[1][1][1]);
            for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++)
                printf("P%d%d %d %d %d %d\n", a, b, rp[a][b].nNumRemovedProtons, rp[a][b].nNumRemovedIsotopicH[0], rp[a][b].nNumRemovedIsotopicH[1], rp[a][b].nNumRemovedIsotopicH[2]);
            for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) for (int i = 0; i < nc[a][b] && pp[a][b]; i++) {
                printf("C%d%d.%d ", a, b, i); dump_full(&pp[a][b][i]);
            }
        }
    }
    return 0;
}
"""
HARNESS = r2.HARNESS.split("int main(void)")[0] + MAIN

R = random.Random(99)
hx = r2.hx


def mutate_inchi(s):
    c = list(s)
    alpha = list("0123456789/,-+;*()=.mnMNefhitbsoqpzABCDEHNOPSc") + ["\n", "\r", "\\"]
    alpha = [a for a in alpha if a != "z"]
    for _ in range(R.choice([0, 0, 0, 1, 1, 2, 4])):
        k = R.randrange(0, len(c) + 1)
        op = R.choice(["del", "ins", "rep", "swap"])
        ch = R.choice(alpha)
        if op == "del" and c:
            c.pop(min(k, len(c) - 1))
        elif op == "ins":
            c.insert(k, ch)
        elif op == "swap" and len(c) > 1:
            j = min(k, len(c) - 2)
            c[j], c[j + 1] = c[j + 1], c[j]
        elif c:
            c[min(k, len(c) - 1)] = ch
    return "".join(c)


def seg_mutate(s):
    """segment-level mutations (drop, repeat, reorder layers)"""
    ly = s.split("/")
    r = R.random()
    if r < 0.15 and len(ly) > 2:
        ly.pop(R.randrange(1, len(ly)))
    elif r < 0.25 and len(ly) > 2:
        i = R.randrange(1, len(ly))
        ly.insert(i, ly[i])
    elif r < 0.32 and len(ly) > 3:
        i = R.randrange(1, len(ly) - 1)
        ly[i], ly[i + 1] = ly[i + 1], ly[i]
    elif r < 0.37:
        ly = ly[: R.randrange(1, len(ly) + 1)]
    return "/".join(ly)


def main():
    corpus = [json.loads(l) for l in open(ROOT / "tests/data/inchi_corpus.jsonl")]
    inchis = sorted({c["inchi"] for c in corpus if c.get("inchi", "").startswith("InChI=1") and "/z" not in c["inchi"]})
    cases = []
    for _ in range(26000):
        ic = R.choice(inchis)
        if R.random() < 0.3:
            ic = ic + R.choice(["\n", "\r\n", "\n" + R.choice(inchis), "\tfoo", "\\AAA", "\\BCD"])
        if R.random() < 0.05:
            ic = "Structure #1\n" + ic
        if R.random() < 0.04:
            ic = "junk line\n" + ic
        if R.random() < 0.5:
            ic = seg_mutate(ic)
        if R.random() < 0.5:
            ic = mutate_inchi(ic)
        cases.append(f"13 {hx(ic.encode())}")
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        kept, out = [], []
        for c in cases:
            try:
                r1 = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=2)
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
    path = ROOT / "tests/data/inchi_read3.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
