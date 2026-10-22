"""Generates tests/data/inchi_read5.json from the InChI C code: the
preparation of the read InChI for the structure restoration of ichiread.c
and ichirvr7.c: `SetProtonsAndXchgIsoH`, `RemoveFixHInChIIdentical2MobH`,
`DetectInpInchiCreationOptions`, `MarkDisconectedIdenticalToReconnected`
after `InChILine2Data`. The harness includes ichiread.c and ichirvr7.c.
Links the upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_read5_reference.py
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


r3 = load("read3", "tools/inchi_read3_reference.py")
r2 = r3.r2
ma = r3.ma
SRC = ma.SRC

PRE = r3.MAIN.split("int main(void)")[0]
MAIN = PRE + r"""
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int t;
        printf("#CASE %d\n", cs);
        scanf("%d", &t);
        if (t == 15) {
            int mode, reqnt, reqst, reqpe; scanf("%d %d %d %d %s", &mode, &reqnt, &reqst, &reqpe, hx);
            int n = unhex(hx, buf);
            INCHI_IOSTREAM ios; inchi_ios_init(&ios, INCHI_IOS_TYPE_STRING, NULL);
            ios.s.pStr = buf; ios.s.nUsedLength = n; ios.s.nPtr = 0;
            SEGM_LINE line; memset(&line, 0, sizeof(line));
            char *pstr = NULL; int state = -77, nerr = -9;
            InpInChI one; memset(&one, 0, sizeof(one));
            int stdf = -5, hso = -5; unsigned char bits = 0;
            int r = InChILine2Data(&ios, &line, &pstr, &state, &nerr, one.pInpInChI, one.nNumComponents, one.nNumProtons, one.s, 1, 1, mode, &stdf, &hso, &bits, &one.polymer, &one.v3000);
            printf("R %d %d %d\n", r, state, nerr);
            int ok = (r == RI_ERR_EOL || r == RI_ERR_EOF) && !nerr && (one.nNumComponents[0][0] + one.nNumComponents[0][1]);
            if (ok) {
                int numc[2]; MODE_PIXH mp[2];
                int r2 = SetProtonsAndXchgIsoH(1, 1, reqpe ? reqnt : reqnt, reqnt, reqst, numc, mp, &one);
                printf("SP %d %d %d %d %d\n", r2, numc[0], numc[1], mp[0], mp[1]);
                if (r2 >= 0) {
                    RemoveFixHInChIIdentical2MobH(&one);
                    int hr, hm, hf, nm, nt;
                    int rd = DetectInpInchiCreationOptions(&one, &hr, &hm, &hf, &nm, &nt);
                    printf("D %d %d %d %d %d %d\n", rd, hr, hm, hf, nm, nt);
                    int nmk = MarkDisconectedIdenticalToReconnected(&one);
                    printf("M %d\n", nmk);
                }
            }
            printf("N %d %d %d %d\n", one.nNumComponents[0][0], one.nNumComponents[0][1], one.nNumComponents[1][0], one.nNumComponents[1][1]);
            for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) {
                REM_PROTONS *p = &one.nNumProtons[a][b];
                printf("P%d%d %d %d %d %d |", a, b, p->nNumRemovedProtons, p->nNumRemovedIsotopicH[0], p->nNumRemovedIsotopicH[1], p->nNumRemovedIsotopicH[2]);
                if (p->pNumProtons) { int n2 = one.nNumComponents[a][0] > one.nNumComponents[a][1] ? one.nNumComponents[a][0] : one.nNumComponents[a][1];
                    for (int q = 0; q < n2; q++) printf(" %d,%d,%d,%d", p->pNumProtons[q].nNumRemovedProtons, p->pNumProtons[q].nNumRemovedIsotopicH[0], p->pNumProtons[q].nNumRemovedIsotopicH[1], p->pNumProtons[q].nNumRemovedIsotopicH[2]); } else printf(" NULL");
                printf("\n");
            }
            for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) for (int i = 0; i < one.nNumComponents[a][b] && one.pInpInChI[a][b]; i++) {
                printf("C%d%d.%d ", a, b, i); dump_full(&one.pInpInChI[a][b][i]); printf("L%d%d.%d %d\n", a, b, i, one.pInpInChI[a][b][i].nLink);
            }
        }
    }
    return 0;
}
"""
HARNESS = r2.HARNESS.split("int main(void)")[0] + MAIN

R = random.Random(1357)
hx = r2.hx


def main():
    corpus = [json.loads(l) for l in open(ROOT / "tests/data/inchi_corpus.jsonl")]
    inchis = sorted({c["inchi"] for c in corpus if c.get("inchi", "").startswith("InChI=1") and "/z" not in c["inchi"] and "AuxInfo" not in c["inchi"]})
    r3.R = R
    cases = []
    for _ in range(30000):
        ic = R.choice(inchis)
        if R.random() < 0.4:
            ic = r3.seg_mutate(ic)
        if R.random() < 0.2:
            ic = r3.mutate_inchi(ic)
        ic += R.choice(["", "", "\n"])
        mode = R.choice([0, 0, 0, 0x200, 0x400])
        reqnt = R.choice([1, 1, 1, 0])
        reqst = R.choice([1, 1, 1, 0])
        cases.append(f"15 {mode} {reqnt} {reqst} 0 {hx(ic.encode())}")
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        kept, out = [], []
        for c in cases:
            try:
                r1 = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=2)
            except subprocess.TimeoutExpired:
                continue
            if r1.returncode != 0:
                continue
            kept.append(c)
            out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode(errors="replace").split("\n")[1:]
        print(len(cases) - len(kept), "cases dropped")
        cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_read5.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
