"""Generates tests/data/inchi_prt1a.json from the InChI C code: the tag
helpers and the first output layers of ichiprt1.c (`szGetTag`,
`str_LineEnd`, `GetSaveOptLetters`, `set_line_separators`,
`MergeZzInHillFormula`, `OutputINCHI_VersionAndKind`,
`OutputINCHI_MainLayerFormula/Connections/Hydrogens`,
`OutputINCHI_ChargeAndRemovedAddedProtonsLayers`). The harness includes
ichiprt1.c itself (its layer functions are static). Links the upstream
InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_prt1a_reference.py
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


p3 = load("prt3a", "tools/inchi_prt3a_reference.py")
m2 = p3.m2
ma = p3.ma
SRC = ma.SRC

HARNESS = p3.HARNESS.split("int main(void)")[0] + r"""
#include "ichiprt1.c"
static void rd_str(char *s) { scanf("%s", s); }
int main(void) {
    int ncases;
    char cmd[32];
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        INCHI_IOS_STRING sb; inchi_strbuf_init(&sb, 16, 16);
        INCHI_IOSTREAM out; inchi_ios_init(&out, INCHI_IOS_TYPE_STRING, NULL);
        INCHI_OUT_CTL io; memset(&io, 0, sizeof(io));
        CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
        int ret = 0;
        scanf("%s", cmd);
        int n = rd(), ot = rd(), mult = rd();
        INCHI_SORT *a = rd_sorts(n);
        INCHI_SORT *b = rd_sorts(n);
        io.ATOM_MODE = rd(); io.TAUT_MODE = rd(); io.nTag = rd(); io.bFhTag = rd(); io.bPlainTextTags = rd();
        io.bSecondNonTautPass = rd(); io.bOmitRepetitions = rd(); io.iCurTautMode = rd(); io.nCurINChISegment = rd();
        io.nNumRemovedProtons = rd(); io.n_pzz = rd(); io.n_zy = rd();
        int nc2[2]; nc2[0] = rd(); nc2[1] = rd();
        int rec = rd(), outopts = rd(), beta = rd();
        io.bOverflow = rd();
        io.sDifSegs[io.nCurINChISegment][DIFS_h_H_ATOMS] = rd();
        io.sDifSegs[io.nCurINChISegment][DIFS_q_CHARGE] = rd();
        io.sDifSegs[io.nCurINChISegment][DIFS_p_PROTONS] = rd();
        io.bTag1 = rd();
        io.bUseMulipliers = mult; io.bOutType = ot; io.num_components = n;
        io.pINChISort = a; io.pINChISort2 = b;
        char *pLF, *pTAB;
        set_line_separators( outopts, &pLF, &pTAB );
        printf("#CASE %d\n", c);
        zero_stack();
        if (!strcmp(cmd, "tag")) {
            int always = -7; char tag[256]; tag[0] = 0;
            szGetTag( io.bTag1 & 0x40000 ? AuxLbl : IdentLbl, io.nTag, io.bTag1 & 0x3ffff, tag, &always );
            printf("%s %d\n", tag, always);
        } else if (!strcmp(cmd, "opt")) {
            char l1, l2; GetSaveOptLetters( (unsigned char)io.bTag1, &l1, &l2 );
            printf("%c%c %d %d %d\n", l1, l2, pLF[0], pTAB[0], outopts);
        } else if (!strcmp(cmd, "lineend")) {
            inchi_strbuf_printf(&sb, "%s", io.sDifSegs[0][0] ? "abc" : "");
            if (io.bTag1 & 1) inchi_strbuf_printf(&sb, "x,y");
            char tg[16]; strcpy(tg, io.bTag1 & 2 ? "/c" : io.bTag1 & 4 ? "/qq" : "");
            ret = str_LineEnd( tg, &io.bOverflow, &sb, io.nCurINChISegment - 2, io.bPlainTextTags );
            printf("%d %d [%s] %d\n", ret, io.bOverflow, sb.pStr ? sb.pStr : "", sb.nUsedLength);
        } else if (!strcmp(cmd, "ver")) {
            ret = OutputINCHI_VersionAndKind( &out, &sb, outopts, beta, pLF, pTAB );
            printf("%d [%s] [%s]\n", ret, out.s.pStr ? out.s.pStr : "", sb.pStr ? sb.pStr : "");
        } else {
            if (!strcmp(cmd, "formula")) ret = OutputINCHI_MainLayerFormula( &cg, &out, &sb, nc2, &rec, &io, pLF, pTAB );
            else if (!strcmp(cmd, "conn")) ret = OutputINCHI_MainLayerConnections( &cg, &out, &sb, nc2, &rec, &io, pLF, pTAB );
            else if (!strcmp(cmd, "hyd")) ret = OutputINCHI_MainLayerHydrogens( &cg, &out, &sb, nc2, &rec, &io, pLF, pTAB );
            else if (!strcmp(cmd, "charge")) ret = OutputINCHI_ChargeAndRemovedAddedProtonsLayers( &cg, &out, &sb, &io, pLF, pTAB );
            printf("%d %d %d %d %d %d %d %d [%s] [%s] [%s]\n", ret, io.bOverflow, io.tot_len, io.tot_len2, io.bTag1, io.bAlways,
                   io.bNonTautNonIsoIdentifierNotEmpty, io.nSegmAction, io.szTag1,
                   out.s.pStr ? out.s.pStr : "", sb.pStr ? sb.pStr : "");
        }
    }
    return 0;
}
"""

R = random.Random(31337)
FORMS = ["C2H6O", "CH4", "C6H6", "ClH", "H2O", "Zz", "Zz2", "C2H6OZz", "CZz2", "Zz3Zz4", "H2Zz2", "OZz2Zz5", "C3H8O2"]


def main():
    cases = []
    for _ in range(7000):
        cmd = R.choice(["formula", "conn", "hyd", "charge", "tag", "opt", "lineend", "ver", "formula", "charge"])
        pool = [p3.base(R) for _ in range(R.randrange(1, 4))]
        for c in pool:
            c["f"] = R.choice(FORMS)
        n = R.randrange(1, 5)
        a = p3.sorts(R, pool, n)
        b = p3.sorts(R, pool, n)
        bt = R.choice([0, R.randrange(0, 1 << 19), R.randrange(0, 1 << 19), 1 << R.randrange(0, 19)])
        if cmd == "tag" and R.random() < 0.3:
            bt |= 1 << 18  # aux table selector
        if cmd == "tag":
            bt = (bt & 0x3ffff) | (R.choice([0, 0x40000]))
            # keep the number of selected comment bits small (64-byte buffer)
            bits = [i for i in range(18) if bt >> i & 1]
            for i in bits[3:]:
                bt &= ~(1 << i)
        if cmd == "opt":
            bt = R.randrange(0, 64)
        segs = [R.choice([0, 1, 2, 3, 4, 5, 7, 8, 12, 15, 16, 17, 31, 63, 64, 127, 255, R.randrange(256)]) for _ in range(3)]
        cases.append(" ".join(map(str, [
            cmd, n, R.choice([0, 1, 2, 3, 4]), R.choice([0, 1, 1])])) + f" {p3.tok_sorts(a)} {p3.tok_sorts(b)} " +
            " ".join(map(str, [
                R.choice(p3.MODES), R.choice([0, 2]), R.choice([0, 1, 2, 3, 4]), R.choice([0, 0, 1, 2]),
                R.choice([0, 1, 1]), R.choice([0, 1]), R.choice([0, 1, 1]), R.choice([0, 1]), R.choice([0, 0, 1]),
                R.choice([0, 1, -1, 2, -3, 0]), R.choice([0, 1, 1, 2]), R.choice([0, 1, 1, 2]),
                R.choice([0, 1]), R.choice([0, 1]), R.choice([0, 1]),
                R.choice([0, 0x40, 0x80, 0xC0, 0x4000, 0x40C0, 0x400, 0x4040]), R.choice([0, 0, 1]),
                R.choice([0, 0, 0, 1]), *segs, bt])))
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
    path = ROOT / "tests/data/inchi_prt1a.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
