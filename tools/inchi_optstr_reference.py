"""Generates tests/data/inchi_optstr.json from the InChI C code:
`parse_options_string` (inchi_dll.c) on random option strings with quotes,
backslashes, spaces and tabs. The function text is cut out of inchi_dll.c.
Needs a C compiler and the fetched InChI sources."""

import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_API/libinchi/src/inchi_dll.c"

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
@FUNC@
static void hex(const char *s) { printf("x"); for (; *s; s++) printf("%02x", (unsigned char) *s); }
int main(void) {
    int n;
    scanf("%d\n", &n);
    for (int c = 0; c < n; c++) {
        char line[512], buf[512];
        int maxargs;
        const char *argv[64];
        scanf("%d", &maxargs); getchar();
        fgets(line, sizeof line, stdin);
        line[strcspn(line, "\n")] = 0;
        for (char *q = line; *q; q++) { if (*q == '\x01') *q = '\t'; }
        strcpy(buf, line);
        int r = parse_options_string(buf, argv, maxargs);
        printf("%d", r);
        for (int i = 0; i < r; i++) { printf(" "); hex(argv[i]); }
        printf("\n");
    }
    return 0;
}
"""

PIECES = ["/Key", " ", " ", "  ", "\t", '"', '"', '""', "\\", "\\\\", '\\"', "SNON", "a b", "x", "/FixedH", "'", "=", "FrameShift:Move Stars"]


def main():
    rng = random.Random(4)
    cases = []
    for _ in range(600):
        s = "".join(rng.choice(PIECES) for _ in range(rng.randrange(0, 12)))
        cases.append((rng.choice([2, 3, 4, 6, 32]), s))
    src = SRC.read_text(errors="replace")
    a = src.index("int parse_options_string( char *cmd, const char *argv[], int maxargs )\n{")
    fn = src[a:src.index("/****", a)]
    stdin = f"{len(cases)}\n" + "".join(f"{m} {s.replace(chr(9), chr(1))}\n" for m, s in cases)
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNC@", fn))
        subprocess.run(["cc", "-w", "-o", str(exe), str(Path(d) / "h.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    data = []
    for (m, s), o in zip(cases, out):
        f = o.split(" ")
        data.append({"max": m, "s": "x" + s.encode().hex(), "n": int(f[0]), "argv": f[1:]})
    path = ROOT / "tests/data/inchi_optstr.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
