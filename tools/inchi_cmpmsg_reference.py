"""Generates tests/data/inchi_cmpmsg.json from the InChI C code:
`AddOneMsg`/`FillOutCompareMessage` (ichirvr7.c) on random difference bits.

    python tools/inchi_cmpmsg_reference.py
"""
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src"
text = (SRC / "ichirvr7.c").read_text().splitlines()
start = next(i for i, l in enumerate(text) if l.startswith("/* message group names */"))
end = next(i for i, l in enumerate(text) if l.startswith("#endif") and i > start)
body = "\n".join(text[start:end])
PRE = r"""
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "inchicmp.h"
typedef unsigned long long INCHI_MODE;
#define TAUT_NON 0
#define TAUT_YES 1
#define TAUT_NUM 2
#define ICHICONST const
"""
MAIN = r"""
int main(void) {
    int n;
    scanf("%d", &n);
    for (int c = 0; c < n; c++) {
        unsigned long long b0, b1;
        int len;
        char hex[900];
        scanf("%llu %llu %d %899s", &b0, &b1, &len, hex);
        char msg[1024];
        int nn = 0;
        if (strcmp(hex, "-")) for (; hex[nn * 2]; nn++) { unsigned v; sscanf(hex + nn * 2, "%2x", &v); msg[nn] = (char)v; }
        msg[nn] = 0;
        INCHI_MODE bits[2] = { b0, b1 };
        int r = FillOutCompareMessage(msg, len, bits);
        printf("%d|%s\n", r, msg);
    }
    return 0;
}
"""
R = random.Random(4242)
ALL = [1 << i for i in range(30)]
cases = []
for _ in range(3000):
    def bits():
        if R.random() < 0.2:
            return 0
        return sum(b for b in ALL if R.random() < R.choice([0.05, 0.15, 0.4]))
    init = R.choice(["-", "-", "-", "Abc", " Problems/mismatches: Mobile-H(x)", "x" * R.randrange(1, 200)])
    cases.append([bits(), bits(), R.choice([40, 120, 250, 400, 1000]), init])
with tempfile.TemporaryDirectory() as d:
    c = Path(d) / "t.c"
    c.write_text(PRE + body + MAIN)
    exe = Path(d) / "t"
    subprocess.run(["cc", "-w", "-I", str(SRC), "-o", str(exe), str(c)], check=True)
    inp = f"{len(cases)}\n" + "\n".join(f"{a} {b} {l} {s.encode().hex() if s != '-' else '-'}" for a, b, l, s in cases)
    out = subprocess.run([str(exe)], input=inp, capture_output=True, text=True, check=True).stdout.splitlines()
res = []
for (a, b, l, s), o in zip(cases, out):
    r, m = o.split("|", 1)
    res.append({"in": [a, b, l, s], "out": [int(r), m]})
(ROOT / "tests/data/inchi_cmpmsg.json").write_text(json.dumps(res))
print(len(res), sum(1 for x in res if x["out"][0]))
