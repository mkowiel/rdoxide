"""Generates tests/data/inchi_parm.json from the InChI C code:
`set_common_options_by_parg` (ichiparm.c) for every option name, in several
letter cases and with numeric arguments, from random initial states. The
function text is cut out of ichiparm.c (the rest of that file is the command
line program). Needs a C compiler and the fetched InChI sources."""

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
@FUNC@
int main(void) {
    int n;
    scanf("%d", &n);
    for (int c = 0; c < n; c++) {
        char arg[256];
        int dev;
        unsigned long long vm, nm;
        int v[34];
        scanf("%255s %d %llu %llu", arg, &dev, &vm, &nm);
        for (int i = 0; i < 34; i++) scanf("%d", &v[i]);
        INPUT_PARMS ip; memset(&ip, 0, sizeof ip);
        INCHI_MODE bVer1 = vm; int mode = nm;
        int o1 = v[0], o2 = v[1], std = v[2], key = v[3], x1 = v[4], x2 = v[5], sp3 = v[6], fb2 = v[7], ph = v[8], ars = v[9], nol = v[10], pe = v[11], dna = v[12], fcf = v[13], rec = v[14], ket = v[15], t15 = v[16], p6 = v[17], p13 = v[18], p16 = v[19], p18 = v[20], p22 = v[21], p39 = v[22], loose = v[23], lm = v[24], pol = v[25], fold = v[26], fss = v[27], saz = v[28], npzz = v[29], nowarn = v[30], mh = v[31], hide = v[32];
        ip.nInputType = v[33];
        int got = set_common_options_by_parg(arg, dev, &ip, &bVer1, &mode, &o1, &o2, &std, &key, &x1, &x2, &sp3, &fb2, &ph, &ars, &nol, &pe, &dna, &fcf, &rec, &ket, &t15, &p6, &p13, &p16, &p18, &p22, &p39, &loose, &lm, &pol, &fold, &fss, &saz, &npzz, &nowarn, &mh, &hide);
        printf("%d %llu %d", got, (unsigned long long) bVer1, mode);
        int out[33] = {o1, o2, std, key, x1, x2, sp3, fb2, ph, ars, nol, pe, dna, fcf, rec, ket, t15, p6, p13, p16, p18, p22, p39, loose, lm, pol, fold, fss, saz, npzz, nowarn, mh, hide};
        for (int i = 0; i < 33; i++) printf(" %d", out[i]);
        printf(" %d %ld %ld %d %d %d\n", ip.nInputType, ip.first_struct_number, ip.last_struct_number, ip.bSaveAllGoodStructsAsProblem, ip.bFilterSS, ip.bFixNonUniformDraw);
    }
    return 0;
}
"""

NAMES = ["INPAUX", "START:12", "END:-5", "RECORD:7", "RECORD:x", "NOLABELS", "SAVEOPT", "AUXNONE", "MISMATCHISERROR", "OUTERRINCHI", "Key", "XHash1", "XHash2", "SNON", "NEWPSOFF", "DONOTADDH", "LooseTSACheck", "SREL", "SRAC", "SUCF", "ChiralFlagON", "ChiralFlagOFF", "SUU", "SLUUD", "FIXEDH", "RECMET", "KET", "15T", "PT_22_00", "PT_16_00", "PT_06_00", "PT_39_00", "PT_13_00", "PT_18_00", "LargeMolecules", "Polymers", "Polymers105", "NPZz", "NoWarnings", "MergeHash", "NoInChI", "HideInChI", "FoldCRU", "FoldSRU", "FrameShift:", "FrameShift:None", "FrameShift: Cyclize ", "FrameShift:MoveStars", "FrameShift:MoveBrackets", "FrameShift:bogus", "NoFrameShift", "NoEdits", "SATZZ", "PGO", "FilterSS", "InvFilterSS", "FNUDOFF", "FixSp3bugOFF", "FBOFF", "FB2OFF", "SPXYZOFF", "SASXYZOFF", "Polymers105+", "Unknown", "", "start:3", "end:", "KEY", "snon", "Polymers1", "S"]


def main():
    rng = random.Random(8)
    cases = []
    for name in NAMES:
        for dev in (0, 1):
            for _ in range(4):
                cases.append({"arg": name, "dev": dev, "vm": rng.randrange(0, 1 << 20), "nm": rng.randrange(0, 1 << 20),
                              "v": [rng.choice([0, 1, 0x8000, 0xff, 4, 8, 7]) for _ in range(33)] + [rng.choice([0, 0, 4, 2])]})
    # the C code reads the argument with %s, so empty / blank arguments are given as "-"
    stdin = f"{len(cases)}\n"
    for c in cases:
        arg = c["arg"] if c["arg"] else "~EMPTY~"
        stdin += f"{arg.replace(' ', chr(1))} {c['dev']} {c['vm']} {c['nm']} " + " ".join(map(str, c["v"])) + "\n"
    src = (SRC / "ichiparm.c").read_text(errors="replace")
    a = src.index("int set_common_options_by_parg(const char* pArg,\n    int  developer_options,", src.index("Returns 1 if pArg recognized"))
    fn = src[a:src.index("Parse command line and set processing options")]
    fn = fn[: fn.rindex("/****")]
    # the harness maps the placeholder bytes back inside the C function's input
    harness = HARNESS.replace("@FUNC@", fn).replace(
        'scanf("%255s %d %llu %llu", arg, &dev, &vm, &nm);',
        'scanf("%255s %d %llu %llu", arg, &dev, &vm, &nm); if (!strcmp(arg, "~EMPTY~")) arg[0] = 0; for (char *q = arg; *q; q++) if (*q == 1) *q = 32;')
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(harness)
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    for c, o in zip(cases, out):
        c["out"] = [int(x) for x in o.split()]
    path = ROOT / "tests/data/inchi_parm.json"
    path.write_text(json.dumps(cases, separators=(",", ":")))
    print(len(cases), "cases ->", path)


main()
