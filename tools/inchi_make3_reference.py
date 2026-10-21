"""Generates tests/data/inchi_make3.json from the InChI C code: `Eql_INChI_Stereo`,
`CompINChILayers`, `MarkUnusedAndEmptyLayers` and `INChI_SegmentAction` on
random pairs of related INChI structures (see inchi_make2_reference.py).

    python tools/inchi_make3_reference.py
"""
import copy
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("make2", ROOT / "tools/inchi_make2_reference.py")
m2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m2)
ma, SRC = m2.ma, m2.SRC

head, _ = m2.HARNESS.split("int main(void) {")
MAIN = r"""
int Eql_INChI_Stereo( INChI_Stereo *s1, int eql1, INChI_Stereo *s2, int eql2, int bRelRac );
int CompINChILayers( const INCHI_SORT *p1, const INCHI_SORT *p2, char sDifSegs[][DIFS_LENGTH], int bFixTranspChargeBug );
int MarkUnusedAndEmptyLayers( char sDifSegs[][DIFS_LENGTH] );
int INChI_SegmentAction( char cDifSegs );
int main(void) {
    int ncases = rd();
    for (int k = 0; k < ncases; k++) {
        INCHI_SORT p[2];
        memset(p, 0, sizeof(p));
        int ord1 = rd(), ord2 = rd(), fix = rd();
        for (int q = 0; q < 2; q++) { p[q].pINChI[0] = rd_inchi(); p[q].pINChI[1] = rd_inchi(); }
        p[0].ord_number = ord1; p[1].ord_number = ord2;
        printf("#CASE %d\n", k);
        INChI *a = p[0].pINChI[1] ? p[0].pINChI[1] : p[0].pINChI[0];
        INChI *b = p[1].pINChI[1] ? p[1].pINChI[1] : p[1].pINChI[0];
        if (a && b) {
            int e2[4] = { 8, 2, 4, 1 };
            for (int x = 0; x < 4; x++) for (int y = 0; y < 4; y++)
                printf("%d", Eql_INChI_Stereo(a->Stereo, e2[x], b->StereoIsotopic, e2[y], 0));
            printf("\n");
        } else printf("-\n");
        char s[DIFL_LENGTH][DIFS_LENGTH];
        memset(s, 0, sizeof(s));
        int ret = CompINChILayers(&p[0], &p[1], s, fix);
        printf("%d", ret);
        for (int l = 0; l < DIFL_LENGTH; l++) for (int i = 0; i < DIFS_LENGTH; i++) printf(" %d", s[l][i]);
        printf("\n");
        MarkUnusedAndEmptyLayers(s);
        for (int l = 0; l < DIFL_LENGTH; l++) for (int i = 0; i < DIFS_LENGTH; i++) printf(" %d", s[l][i]);
        printf("\n");
        for (int v = 0; v < 16; v++) printf(" %d", INChI_SegmentAction(v));
        printf("\n");
    }
    return 0;
}
"""
HARNESS = (head + MAIN).replace("c->nNum_H_fixed = calloc(n + 1, 1)", "c->nNum_H_fixed = calloc(64, 1)")


def main():
    rng = random.Random(33)
    cases = []
    for _ in range(10000):
        base = m2.comp(rng)
        if rng.random() < 0.5:
            base["sti"] = copy.deepcopy(base["st"]) if rng.random() < 0.7 else base["sti"]
        if rng.random() < 0.3:
            base["fixed"] = [rng.choice([0, 0, 1]) for _ in range(base["n"])]
        slots = []
        for _p in range(2):
            has_t = rng.random() < 0.6
            has_n = rng.random() < 0.7
            slots.append(m2.mutate(rng, base) if has_t else None)
            slots.append(m2.mutate(rng, base) if has_n else None)
        pre = f"{rng.randrange(0, 3)} {rng.randrange(0, 3)} {rng.choice([0, 1, 1])} "
        cases.append(pre + " ".join(m2.tok(s) for s in slots))
    stdin = f"{len(cases)}\n" + "\n".join(cases) + "\n"
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True)
    out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_make3.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
