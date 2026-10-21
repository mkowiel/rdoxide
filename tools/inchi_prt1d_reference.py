"""Generates tests/data/inchi_prt1d.json from the InChI C code: the AuxInfo
output layers of ichiprt1.c (`OutputAUXINFO_HeaderAndNormalization_type`,
`..._OriginalNumbersAndEquivalenceClasses`, `..._TautomericGroupsEquivalence`,
`..._Stereo`, `..._IsotopicInfo`, `..._ChargesRadicalsAndUnusualValences`).
The harness includes ichiprt1.c itself (the layer functions are static).
Links the upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_prt1d_reference.py
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


p3d = load("prt3d", "tools/inchi_prt3d_reference.py")
p3c = p3d.p3c
p3 = p3c.p3
ma = p3c.ma
SRC = ma.SRC

ARRAYS = ["bAtomEqu", "bTautEqu", "bInvStereo", "bInvStereoOrigNumb", "bIsotopicOrigNumb", "bIsotopicAtomEqu",
          "bIsotopicTautEqu", "bInvIsotopicStereo", "bInvIsotopicStereoOrigNumb", "bIgn_UU_Sp3_Iso",
          "bIgn_UU_Sp2_Iso", "bChargesRadVal"]

HARNESS = p3c.HARNESS.split("int str_AuxEqu(")[0] + r"""
/* each str_* call of a layer function starts with a zeroed stack (some of the
   C buffers are read uninitialised for inconsistent inputs) */
#define inchi_strbuf_reset(b) (zero_stack(), inchi_strbuf_reset(b))
#include "ichiprt1.c"
int main(void) {
    int ncases;
    char cmd[32];
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        INCHI_IOS_STRING sb; inchi_strbuf_init(&sb, 16, 16);
        INCHI_IOSTREAM out; inchi_ios_init(&out, INCHI_IOS_TYPE_STRING, NULL);
        INCHI_OUT_CTL io; memset(&io, 0, sizeof(io));
        CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
        int ret = 0, flags = 0;
        scanf("%s", cmd);
        int n = rd(), ot = rd(), mult = rd();
        INCHI_SORT *a = rd_sorts(n);
        INCHI_SORT *b = rd_sorts(n);
        io.TAUT_MODE = rd(); io.nTag = rd(); io.bFhTag = rd(); io.bPlainTextTags = rd();
        io.bSecondNonTautPass = rd(); io.bOmitRepetitions = rd(); io.iCurTautMode = rd();
        io.bTautomeric = rd(); io.bTautomericOutputAllowed = rd(); io.bIsotopic = rd();
        flags = rd();
        int nc2[2]; nc2[0] = rd(); nc2[1] = rd();
        int outopts = rd(), rec = rd();
        io.bOverflow = rd();
""" + "\n".join(f"        for (int i = 0; i < 2; i++) io.{a}[i] = rd();" for a in ARRAYS) + r"""
        io.bUseMulipliers = mult; io.bOutType = ot; io.num_components = n;
        io.pINChISort = a; io.pINChISort2 = b; io.pSortPrintINChIFlags = &flags;
        char *pLF, *pTAB;
        set_line_separators( outopts, &pLF, &pTAB );
        printf("#CASE %d\n", c);
        zero_stack();
        if (!strcmp(cmd, "head")) ret = OutputAUXINFO_HeaderAndNormalization_type( &cg, &out, &sb, outopts, &rec, nc2, &io, pLF, pTAB );
        else if (!strcmp(cmd, "orig")) ret = OutputAUXINFO_OriginalNumbersAndEquivalenceClasses( &cg, &out, &sb, nc2, &io, pLF, pTAB );
        else if (!strcmp(cmd, "tg")) ret = OutputAUXINFO_TautomericGroupsEquivalence( &cg, &out, &sb, &io );
        else if (!strcmp(cmd, "stereo")) ret = OutputAUXINFO_Stereo( &cg, &out, &sb, &io, pLF, pTAB );
        else if (!strcmp(cmd, "iso")) ret = OutputAUXINFO_IsotopicInfo( &cg, &out, &sb, &rec, &io, pLF, pTAB );
        else ret = OutputAUXINFO_ChargesRadicalsAndUnusualValences( &cg, &out, &sb, &io, pLF, pTAB );
        printf("%d %d %d %d %d %d %d [%s] [%s] [%s] [%s] [%s]\n", ret, io.bOverflow, io.tot_len, io.bTag1, io.bTag2, io.bTag3, io.bAlways,
               io.szTag1, io.szTag2, io.szTag3, out.s.pStr ? out.s.pStr : "", sb.pStr ? sb.pStr : "");
    }
    return 0;
}
"""

R = random.Random(99)


def main():
    cases = []
    for _ in range(8000):
        cmd = R.choice(["head", "orig", "tg", "stereo", "iso", "iso", "crv"])
        pool = [p3d.pool_item(R) for _ in range(R.randrange(1, 4))]
        n = R.randrange(1, 5)
        a = p3c.sorts(R, pool, n)
        b = p3c.sorts(R, pool, n)
        ot = R.choice([0, 1, 2, 3, 4])
        cases.append(" ".join(map(str, [cmd, n, ot, R.choice([0, 1, 1])])) + f" {a} {b} " + " ".join(map(str, [
            R.choice([0, 2]), R.choice([0, 1, 2, 3, 4]), R.choice([0, 0, 1, 2]), R.choice([0, 1, 1]),
            R.choice([0, 1]), R.choice([0, 1, 1]), R.choice([0, 1]),
            R.choice([0, 1, 1]), R.choice([0, 1, 1]), R.choice([0, 1, 1]),
            R.choice([0, 0, 16, 32]),
            R.choice([0, 1, 2]), R.choice([0, 1, 2]),
            R.choice([0, 0x200, 0x80, 0x280]), R.choice([0, 1, 1, 2]),
            R.choice([0, 0, 0, 1]),
            *[R.choice([0, 1, 1]) for _ in range(2 * len(ARRAYS))]])))
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
                    out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode(errors="replace").split("\n")[1:]
            print(len(cases) - len(kept), "cases dropped")
            cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases), (len(starts), len(cases))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    data = [x for x in data if all(ch.isprintable() and ord(ch) < 127 for l in x["out"] for ch in l)]
    path = ROOT / "tests/data/inchi_prt1d.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
