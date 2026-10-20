"""Generates tests/data/inchi_readparms.json from the InChI C code:
`ReadCommandLineParms` (ichiparm.c) as built into the library, for random
argument lists made of real and invented options. The text of the function
(with `set_common_options_by_parg`) is cut out of ichiparm.c. Needs a C
compiler and the fetched InChI sources."""

import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src"

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include "mode.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "ichimain.h"
#include "extr_ct.h"
#include "ichicomp.h"
#include "ichi_io.h"
const char gsMissing[] = "is missing";
const char gsEmpty[] = "";
const char gsSpace[] = " ";
const char gsEqual[] = "=";
@FUNC@
static void hex(const char *s) { printf("x"); for (; s && *s; s++) printf("%02x", (unsigned char) *s); }
int main(void) {
    int n;
    scanf("%d", &n);
    for (int c = 0; c < n; c++) {
        int rel, argc;
        scanf("%d %d", &rel, &argc);
        const char *argv[40];
        static char store[40][128];
        argv[0] = "";
        for (int i = 1; i < argc; i++) {
            char h[300];
            scanf("%s", h);
            int m = 0;
            for (const char *q = h + 1; q[0] && q[1]; q += 2) { unsigned v; sscanf(q, "%2x", &v); store[i][m++] = (char) v; }
            store[i][m] = 0;
            argv[i] = store[i];
        }
        INPUT_PARMS ip; errno = 0;
        char sdfval[128] = "";
        unsigned long ul = 0;
        INCHI_IOSTREAM log;
        inchi_ios_init(&log, INCHI_IOS_TYPE_STRING, NULL);
        int r = ReadCommandLineParms(argc, argv, &ip, sdfval, &ul, rel, &log);
        printf("%d", r);
        printf(" %d %d %d %d", ip.nInputType, ip.num_paths, ip.bAbcNumbers, ip.bCtPredecessors);
        printf(" %llu %llu %llu", (unsigned long long) ip.nMode, (unsigned long long) ip.bTautFlags, (unsigned long long) ip.bTautFlagsDone);
        printf(" %d %d %d %d %d %d %d", ip.bINChIOutputOptions, ip.bINChIOutputOptions2, ip.bDisplayEachComponentINChI, ip.bDisplay, ip.bDisplayIfRestoreWarnings, ip.bMergeAllInputStructures, ip.bSaveWarningStructsAsProblem);
        printf(" %d %d %d %d %d %d %d %d", ip.bSaveAllGoodStructsAsProblem, ip.bGetSdfileId, ip.bGetMolfileNumber, ip.bCompareComponents, ip.bDisplayCompositeResults, ip.bDoNotAddH, ip.bNoStructLabels, ip.bFixNonUniformDraw);
        printf(" %d %d %d %d %d %d %d %d", ip.bCalcInChIHash, ip.bChiralFlag, ip.bAllowEmptyStructure, ip.bLargeMolecules, ip.bLooseTSACheck, ip.bPolymers, ip.bFoldPolymerSRU, ip.bFrameShiftScheme);
        printf(" %d %d %d %d %d", ip.bNPZz, ip.bStereoAtZz, ip.bMergeHash, ip.bNoWarnings, ip.bHideInChI);
        printf(" %d %d %d %d %d", ip.bReadInChIOptions, ip.bUnderivatize, ip.bRing2Chain, ip.bIgnoreUnchanged, ip.bFilterSS);
        printf(" %ld %ld %ld %ld %ld %lu %ld", ip.first_struct_number, ip.last_struct_number, ip.msec_MaxTime, ip.msec_LeftTime, ip.ulDisplTime, ip.lSdfId, ip.lMolfileNumber);
        printf(" %d %d ", ip.pSdfLabel != NULL, ip.pSdfValue != NULL);
        hex(ip.szSdfDataHeader);
        for (int i = 0; i < MAX_NUM_PATHS; i++) { printf(" "); if (ip.path[i]) hex(ip.path[i]); else printf("N"); }
        printf(" ");
        hex(log.s.pStr ? log.s.pStr : "");
        printf("\n");
    }
    return 0;
}
"""

OPTS = ["-Key", "-SNON", "-FixedH", "-RecMet", "-KET", "-15T", "-SUU", "-SLUUD", "-SRel", "-SRac", "-SUCF", "-ChiralFlagON", "-ChiralFlagOFF",
        "-AuxNone", "-NoLabels", "-DoNotAddH", "-LooseTSACheck", "-NEWPSOFF", "-OutputSDF", "-SdfAtomsDT", "-WarnOnEmptyStructure",
        "-InChI2InChI", "-InChI2Struct", "-DoDrvReport", "-STDIO", "-SDF:name", "-SDF:", "-EQU", "-D", "-F12", "-Fx", "-XHash1", "-XHash2",
        "-Polymers", "-Polymers105", "-LargeMolecules", "-NPZz", "-NoWarnings", "-MergeHash", "-HideInChI", "-FoldSRU", "-FrameShift:MoveStars",
        "-NoEdits", "-SATZZ", "-SAVEOPT", "-PT_22_00", "-W10", "-W0.5", "-W1e2", "-Wabc", "-W5x", "-W-3", "-W1e400", "-Winf", "-WM500", "-WM", "-WM0", "-WM-4", "-WM12x",
        "-Bogus", "-.ext", "-.a", "-.b", "-.c", "-.d", "-Start:3", "-End:9", "-Record:4", "-INPAUX", "-FNUDOFF", "-PGO", "-FilterSS",
        "mol1.mol", "dir/in.sdf", "NUL", "x", "-", "a.b", "dir.d/file", "/abs/f.mol"]


def main():
    rng = random.Random(17)
    cases = []
    for _ in range(700):
        args = [rng.choice(OPTS) for _ in range(rng.randrange(0, 7))]
        # an extension option after a path makes the C code free a string literal
        args.sort(key=lambda a: 0 if a.startswith('-.') else 1)
        # "NUL" as the first path is used after being freed by the C code
        pos = [i for i, a in enumerate(args) if not a.startswith('-') or a == '-']
        if pos and args[pos[0]] == 'NUL':
            args[pos[0]] = 'x'
        cases.append((rng.choice([0, 1]), args))
    src = (SRC / "ichiparm.c").read_text(errors="replace")
    a = src.index("#define VER103_DEFAULT_MODE")
    b = src.index("Print input parameters")
    b = src.rindex("/****", 0, b)
    fn = src[a:b]
    stdin = f"{len(cases)}\n"
    for rel, args in cases:
        stdin += f"{rel} {len(args) + 1} " + " ".join("x" + x.encode().hex() for x in args) + "\n"
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNC@", fn))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c"), str(SRC / "ichi_io.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    data = [{"rel": rel, "args": args, "out": o} for (rel, args), o in zip(cases, out)]
    path = ROOT / "tests/data/inchi_readparms.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
