"""Generates tests/data/inchi_err.json from the InChI C code: `ErrMsg` for a
range of codes and `AddErrorMessage` for random message sequences (the
buffer is 256 bytes, so overflow is covered). Needs a C compiler and the
fetched InChI sources."""

import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src"

HARNESS = r"""
#include <stdio.h>
#include <string.h>
#include "mode.h"
#include "ichierr.h"
int main(void) {
    int n;
    char line[512];
    if (scanf("%d\n", &n) != 1) return 1;
    for (int i = 0; i < n; i++) {
        int code;
        scanf("%d\n", &code);
        printf("E%s\n", ErrMsg(code));
    }
    int ncases;
    scanf("%d\n", &ncases);
    for (int c = 0; c < ncases; c++) {
        int nm;
        char buf[STR_ERR_LEN] = "";
        scanf("%d\n", &nm);
        for (int m = 0; m < nm; m++) {
            fgets(line, sizeof line, stdin);
            line[strcspn(line, "\n")] = 0;
            int r = AddErrorMessage(buf, line);
            printf("A%d|%s\n", r, buf);
        }
        printf("-\n");
    }
    return 0;
}
"""

WORDS = ["Error", "Warning:", "Cannot", "stereo", "x", "Out of RAM", "Time limit exceeded", "a;", "b:", "message: more", "ab", "Error:", "long message " + "z" * 40]


def main():
    rng = random.Random(7)
    codes = [0, 1, -1, 5, -30000, -30001, -30019, -30020, -30100, -9988, -9986, -9985, -9999, -9990, -29999] + list(range(-30022, -29995))
    cases = []
    for _ in range(300):
        cases.append([rng.choice(WORDS) + ("" if rng.random() < 0.6 else " " + rng.choice(WORDS)) for _ in range(rng.randrange(1, 14))])
    stdin = f"{len(codes)}\n" + "\n".join(map(str, codes)) + f"\n{len(cases)}\n"
    for c in cases:
        stdin += f"{len(c)}\n" + "\n".join(c) + "\n"
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS)
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True, check=True).stdout.split("\n")
    errs = [l[1:] for l in out if l.startswith("E")]
    adds, cur = [], []
    for l in out:
        if l.startswith("A"):
            r, buf = l[1:].split("|", 1)
            cur.append([int(r), buf])
        elif l == "-":
            adds.append(cur)
            cur = []
    data = {
        "codes": codes,
        "messages": errs,
        "cases": [{"in": c, "out": a} for c, a in zip(cases, adds)],
    }
    path = ROOT / "tests/data/inchi_err.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(codes), "codes,", len(cases), "sequences ->", path)


main()
