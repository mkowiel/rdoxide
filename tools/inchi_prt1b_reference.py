"""Generates tests/data/inchi_prt1b.json from the InChI C code: the stereo
and isotopic output layers of ichiprt1.c (`OutputINCHI_StereoLayer`,
`OutputINCHI_IsotopicLayer`). The harness includes ichiprt1.c itself (the
layer functions are static). Links the upstream InChI base library (see
inchi_markalt_reference.py).

    python tools/inchi_prt1b_reference.py
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


p3c = load("prt3c", "tools/inchi_prt3c_reference.py")
p3 = p3c.p3
ma = p3c.ma
SRC = ma.SRC

HARNESS = p3c.HARNESS.split("int str_AuxEqu(")[0] + r"""
#include "ichiprt1.c"
int bin_AuxTautTrans( INCHI_SORT *a, INCHI_SORT *b, AT_NUMB **tn, AT_NUMB **ts, int ot, int n );
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
        for (int i = 0; i < n; i++) a[i].ord_number = rd();
        for (int i = 0; i < n; i++) b[i].ord_number = rd();
        io.TAUT_MODE = rd(); io.nTag = rd(); io.bFhTag = rd(); io.bPlainTextTags = rd();
        io.bSecondNonTautPass = rd(); io.bOmitRepetitions = rd(); io.iCurTautMode = rd(); io.nCurINChISegment = rd();
        io.bAbcNumbers = rd(); io.bOutputType = rd();
        for (int i = 0; i < 2; i++) io.bRelativeStereo[i] = rd();
        for (int i = 0; i < 2; i++) io.bRacemicStereo[i] = rd();
        for (int i = 0; i < 2; i++) io.bIsotopicRelativeStereo[i] = rd();
        for (int i = 0; i < 2; i++) io.bIsotopicRacemicStereo[i] = rd();
        for (int i = 0; i < 3; i++) io.num_iso_H[i] = rd();
        int rec = rd();
        io.bOverflow = rd();
        int cs = io.nCurINChISegment;
        io.sDifSegs[cs][DIFS_b_SBONDS] = rd();
        io.sDifSegs[cs][DIFS_t_SATOMS] = rd();
        io.sDifSegs[cs][DIFS_m_SP3INV] = rd();
        io.sDifSegs[cs][DIFS_s_STYPE] = rd();
        io.sDifSegs[cs][DIFS_i_IATOMS] = rd();
        io.sDifSegs[cs][DIFS_h_H_ATOMS] = rd();
        io.sDifSegs[DIFL_F][DIFS_o_TRANSP] = rd();
        io.bUseMulipliers = mult; io.bOutType = ot; io.num_components = n;
        io.pINChISort = a; io.pINChISort2 = b; io.pSortPrintINChIFlags = &flags;
        char *pLF, *pTAB;
        set_line_separators( rd() ? 0x80 : 0, &pLF, &pTAB );
        printf("#CASE %d\n", c);
        zero_stack();
        if (!strcmp(cmd, "iso") && io.bOutType == OUT_NONTAUT && io.bOutputType == OUT_TN && io.bSecondNonTautPass &&
            INCHI_SEGM_FILL == INChI_SegmentAction( io.sDifSegs[DIFL_F][DIFS_o_TRANSP] )) {
            int has = 0;
            for (int i = 0; i < n; i++) {
                INCHI_SORT *sa = a + i, *sb2 = b + i; int ii = GET_II(ot, sa), ii2 = GET_II(OUT_T1, sb2);
                INChI *p = ii >= 0 ? a[i].pINChI[ii] : NULL, *t = ii2 >= 0 ? b[i].pINChI[ii2] : NULL;
                if (p && p->nNumberOfAtoms > 0 && t && t->nNumberOfAtoms > 0 && a[i].ord_number != b[i].ord_number) has = 1;
            }
            if (!has) { printf("UNDEF\n"); continue; }
        }
        if (!strcmp(cmd, "stereo")) ret = OutputINCHI_StereoLayer( &cg, &out, &sb, &io, pLF, pTAB );
        else ret = OutputINCHI_IsotopicLayer( &cg, &out, &sb, &rec, &io, pLF, pTAB );
        printf("%d %d %d %d %d %d %d %d %d %d %d %d %d [%s] [%s] [%s] [%s] [%s]\n", ret, io.bOverflow, io.tot_len, io.tot_len2,
               io.bTag1, io.bTag2, io.bTag3, io.bAlways, io.bNonTautNonIsoIdentifierNotEmpty, io.bNonTautIsoIdentifierNotEmpty,
               io.nSegmAction, io.bRelRac, flags, io.szTag1, io.szTag2, io.szTag3,
               out.s.pStr ? out.s.pStr : "", sb.pStr ? sb.pStr : "");
    }
    return 0;
}
"""

R = random.Random(2718)


def main():
    cases = []
    for _ in range(8000):
        cmd = R.choice(["stereo", "iso"])
        pool = [p3c.mk_pool_item(R) for _ in range(R.randrange(1, 4))]
        n = R.randrange(1, 5)
        a = p3c.sorts(R, pool, n)
        b = p3c.sorts(R, pool, n)
        ords = " ".join(str(R.randrange(0, n)) for _ in range(2 * n))
        ot = R.choice([0, 1, 2, 3, 4, 4])
        seg = lambda: R.choice([0, 1, 2, 3, 4, 5, 7, 8, 12, 15, 16, 17, 31, 63, 64, 127, 255, R.randrange(256)])
        rb = lambda: R.choice([0, 0, 1])
        cases.append(" ".join(map(str, [cmd, n, ot, R.choice([0, 1, 1])])) + f" {a} {b} {ords} " + " ".join(map(str, [
            R.choice([0, 2]), R.choice([0, 1, 2, 3, 4]), R.choice([0, 0, 1, 2]), R.choice([0, 1, 1]),
            R.choice([0, 1]), R.choice([0, 1, 1]), R.choice([0, 1]), R.choice([0, 0, 1]),
            R.choice([0, 1]), R.choice([0, 3, 3, 1]),
            rb(), rb(), rb(), rb(), rb(), rb(), rb(), rb(),
            R.randrange(0, 3), R.randrange(0, 3), R.randrange(0, 3),
            R.choice([0, 1]), R.choice([0, 0, 0, 1]),
            seg(), seg(), seg(), seg(), seg(), seg(), seg(), R.choice([0, 1])])))
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
    data = [x for x in data if all(ch.isprintable() and ord(ch) < 127 for l in x["out"] for ch in l) and x["out"] != ["UNDEF"]]
    path = ROOT / "tests/data/inchi_prt1b.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
