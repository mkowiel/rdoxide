"""Generates tests/data/inchi_sort.json: inputs and outputs of `inchi_qsort`
and `insertions_sort` of the InChI C sources (upstream/external/InChI), for
`crates/rdoxide/tests/inchi_utils.rs`. Needs a C compiler and the fetched
InChI sources (`tools/sync_upstream.py fetch`).

Elements are (key, id) pairs compared by key only, so the order of equal keys
shows how the unstable sorts treat ties.
"""

import json
import random
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src/ichisort.c"

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct { int key; int id; } EL;
void inchi_swap(char *a, char *b, size_t width);
static int comp(const void *a, const void *b, void *p) {
    return ((const EL *)a)->key - ((const EL *)b)->key;
}
%s
int main(void) {
    int ncases;
    if (scanf("%%d", &ncases) != 1) return 1;
    for (int c = 0; c < ncases; c++) {
        int n, mode;
        scanf("%%d %%d", &n, &mode);
        EL *v = malloc(sizeof(EL) * (n ? n : 1));
        for (int i = 0; i < n; i++) scanf("%%d", &v[i].key), v[i].id = i;
        int trans = 0;
        if (mode == 0) inchi_qsort(NULL, v, n, sizeof(EL), comp);
        else trans = insertions_sort(NULL, v, n, sizeof(EL), comp);
        printf("%%d", trans);
        for (int i = 0; i < n; i++) printf(" %%d", v[i].id);
        printf("\n");
        free(v);
    }
    return 0;
}
"""


def extract(text, start, end):
    a = text.index(start)
    b = text.index(end, a)
    return text[a:b]


def main():
    text = SRC.read_text()
    # the sort functions, without the CANON_GLOBALS-dependent rest
    body = "#define STKSIZ (8*sizeof(void*) - 2)\n"
    body += extract(text, "void inchi_qsort(", "/****************************************************************************\n Sort by insertions\n")
    body += extract(text, "int insertions_sort( void *pCG,", "/****************************************************************************\n Sort by insertions\n****************************************************************************/\nint insertions_sort_AT_NUMBERS")
    rng = random.Random(20261020)
    cases = []
    for n in list(range(0, 40)) + [64, 100, 257, 1000]:
        for modulus in (2, 3, 7, 1000000):
            for mode in (0, 1):
                if mode == 1 and n > 300:
                    continue
                cases.append((mode, [rng.randrange(modulus) for _ in range(n)]))
    # sorted, reversed and constant runs
    for n in (10, 50, 200):
        for keys in (list(range(n)), list(range(n, 0, -1)), [5] * n, [i // 3 for i in range(n)]):
            cases.append((0, keys))
            cases.append((1, keys))
    with tempfile.TemporaryDirectory() as d:
        c = Path(d) / "h.c"
        c.write_text(HARNESS % body)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-O1", "-w", "-o", str(exe), str(c)], check=True)
        stdin = f"{len(cases)}\n" + "".join(
            f"{len(k)} {m}\n{' '.join(map(str, k))}\n" for m, k in cases
        )
        out = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True, check=True).stdout
    results = [list(map(int, line.split())) for line in out.strip("\n").split("\n")]
    data = [
        {"mode": "qsort" if m == 0 else "insertions", "keys": k, "swaps": r[0], "order": r[1:]}
        for (m, k), r in zip(cases, results)
    ]
    path = ROOT / "tests/data/inchi_sort.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(f"{len(data)} cases -> {path} ({path.stat().st_size} bytes)")


if __name__ == "__main__":
    sys.exit(main())
