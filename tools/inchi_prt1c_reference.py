"""Generates tests/data/inchi_prt1c.json from the InChI C code: the fixed-H
output layer of ichiprt1.c (`OutputINCHI_FixedHLayerWithSublayers`). The
harness includes ichiprt1.c itself (the layer function is static). Links the
upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_prt1c_reference.py
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


p1b = load("prt1b", "tools/inchi_prt1b_reference.py")
p3c = p1b.p3c
p3 = p3c.p3
ma = p3c.ma
SRC = ma.SRC

HARNESS = p3c.HARNESS.split("int str_AuxEqu(")[0] + r"""
#include "ichiprt1.c"
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        INCHI_IOS_STRING sb; inchi_strbuf_init(&sb, 16, 16);
        INCHI_IOSTREAM out; inchi_ios_init(&out, INCHI_IOS_TYPE_STRING, NULL);
        INCHI_OUT_CTL io; memset(&io, 0, sizeof(io));
        CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
        int flags = 0, repeat = -1;
        int n = rd(), ot = rd(), mult = rd();
        INCHI_SORT *a = rd_sorts(n);
        INCHI_SORT *b = rd_sorts(n);
        INCHI_SORT *t0 = rd_sorts(n);
        INCHI_SORT *t1 = rd_sorts(n);
        io.ATOM_MODE = rd(); io.TAUT_MODE = rd(); io.nTag = rd(); io.bFhTag = rd(); io.bPlainTextTags = rd();
        io.bSecondNonTautPass = rd(); io.bOmitRepetitions = rd(); io.iCurTautMode = rd(); io.nCurINChISegment = rd();
        io.bTautomeric = rd(); io.bNonTautomeric = rd(); io.bNonTautIsIdenticalToTaut = rd(); io.bOutputType = rd();
        io.num_comp[0] = rd(); io.num_comp[1] = rd();
        io.n_pzz = rd(); io.n_zy = rd();
        int rec = rd();
        io.bOverflow = rd();
        flags = rd();
        io.bNonTautNonIsoIdentifierNotEmpty = rd(); io.bNonTautIsoIdentifierNotEmpty = rd();
        io.sDifSegs[DIFL_F][DIFS_f_FORMULA] = rd();
        io.sDifSegs[DIFL_F][DIFS_h_H_ATOMS] = rd();
        io.bUseMulipliers = mult; io.bOutType = ot; io.num_components = n;
        io.pINChISort = a; io.pINChISort2 = b; io.pSortPrintINChIFlags = &flags;
        io.pINChISortTautAndNonTaut = (INCHI_SORT **) calloc(2, sizeof(INCHI_SORT *));
        io.pINChISortTautAndNonTaut[0] = t0; io.pINChISortTautAndNonTaut[1] = t1;
        char *pLF, *pTAB;
        set_line_separators( rd() ? 0x80 : 0, &pLF, &pTAB );
        printf("#CASE %d\n", c);
        zero_stack();
        int ret = OutputINCHI_FixedHLayerWithSublayers( &cg, &out, &sb, &rec, &io, pLF, pTAB, &repeat );
        int which = io.pINChISort == a ? 0 : io.pINChISort == t0 ? 1 : io.pINChISort == t1 ? 2 : 3;
        printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d [%s] [%s] [%s]\n", ret, repeat, io.bOutType, io.iCurTautMode, io.bSecondNonTautPass,
               io.nCurINChISegment, io.num_components, io.bFhTag, io.bOverflow, io.tot_len, io.tot_len2, io.bTag1, io.bAlways,
               io.bNonTautNonIsoIdentifierNotEmpty * 100 + io.bNonTautIsoIdentifierNotEmpty, io.nSegmAction * 1000 + flags, which,
               io.szTag1, out.s.pStr ? out.s.pStr : "", sb.pStr ? sb.pStr : "");
    }
    return 0;
}
"""

R = random.Random(1618)
FORMS = ["C2H6O", "CH4", "C6H6", "ClH", "H2O", "Zz", "Zz2", "C2H6OZz", "Zz3Zz4", "OZz2Zz5", "C3H8O2"]


def main():
    cases = []
    for _ in range(8000):
        pool = [p3c.mk_pool_item(R) for _ in range(R.randrange(1, 4))]
        for c, _a in pool:
            c["f"] = R.choice(FORMS)
        n = R.randrange(1, 5)
        ss = [p3c.sorts(R, pool, n) for _ in range(4)]
        seg = lambda: R.choice([0, 1, 2, 3, 4, 5, 7, 8, 12, 15, 16, 17, 31, 63, 64, 127, 255, R.randrange(256)])
        ot = R.choice([0, 1, 2, 3, 3, 3, 4, 4])
        cases.append(" ".join(map(str, [n, ot, R.choice([0, 1, 1])])) + " " + " ".join(ss) + " " + " ".join(map(str, [
            R.choice(p3.MODES), R.choice([0, 2]), R.choice([0, 1, 2, 3, 4]), R.choice([0, 0, 1, 2]), R.choice([0, 1, 1]),
            R.choice([0, 1]), R.choice([0, 1, 1]), R.choice([0, 1]), R.choice([0, 1, 2]),
            R.choice([0, 1]), R.choice([0, 1]), R.choice([0, 0, 1]), R.choice([0, 3, 3, 1]),
            R.randrange(0, n + 1), R.randrange(0, n + 1),
            R.choice([0, 1, 2]), R.choice([0, 1, 2]),
            R.choice([0, 1]), R.choice([0, 0, 0, 1]),
            R.choice([0, 1, 2, 4, 8, 16, 32]),
            R.choice([0, 0, 1, 2]), R.choice([0, 0, 1, 2]),
            seg(), seg(), R.choice([0, 1])])))
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
    path = ROOT / "tests/data/inchi_prt1c.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
