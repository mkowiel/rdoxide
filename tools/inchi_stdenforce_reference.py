"""Generates tests/data/inchi_stdenforce.json from the InChI C code: the
"Ensure standardness" block of GetINCHIEx (inchi_dll.c) on random option
states. The block is cut out of inchi_dll.c. Needs a C compiler and the
fetched InChI sources."""

import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src"
DLL = ROOT / "upstream/external/InChI/INCHI_API/libinchi/src/inchi_dll.c"

HARNESS = r"""
#include <stdio.h>
#include <string.h>
#include "mode.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
static void enforce(INPUT_PARMS *ip) {
    int enforce_std_format = 1;
@BLOCK@
}
int main(void) {
    int n;
    scanf("%d", &n);
    for (int c = 0; c < n; c++) {
        INPUT_PARMS ip; memset(&ip, 0, sizeof ip);
        unsigned long long nm, tf; int oo;
        scanf("%llu %llu %d", &nm, &tf, &oo);
        ip.nMode = nm; ip.bTautFlags = tf; ip.bINChIOutputOptions = oo;
        enforce(&ip);
        printf("%llu %llu %d\n", (unsigned long long) ip.nMode, (unsigned long long) ip.bTautFlags, ip.bINChIOutputOptions);
    }
    return 0;
}
"""


def main():
    rng = random.Random(23)
    src = DLL.read_text(errors="replace")
    a = src.index("    /* Ensure standardness */")
    b = src.index("PrintInputParms( log_file, ip );", a)
    block = src[a:b]
    block = block[: block.rindex("    /* */")]
    cases = [(rng.randrange(0, 1 << 20), rng.randrange(0, 1 << 29), rng.randrange(0, 1 << 16)) for _ in range(600)]
    stdin = f"{len(cases)}\n" + "".join(f"{a} {b} {c}\n" for a, b, c in cases)
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@BLOCK@", block))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    data = [{"in": list(c), "out": [int(x) for x in o.split()]} for c, o in zip(cases, out)]
    path = ROOT / "tests/data/inchi_stdenforce.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
