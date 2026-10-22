"""Generates tests/data/inchi_read4.json from the InChI C code:
`InChILine2Data` (ReadInChILine plus the restoration of the omitted
layers, the balance checks) over single InChI strings of the official
RDKit corpus, mutated and with layers dropped. The harness includes
ichiread.c. Links the upstream InChI base library (see
inchi_markalt_reference.py).

    python tools/inchi_read4_reference.py
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

MAIN = r3.MAIN.replace("if (t == 13) {\n            scanf(\"%s\", hx);", "if (t == 14) {\n            int mode; scanf(\"%d %s\", &mode, hx);").replace(
    "int r = ReadInChILine(&ios, &line, &pstr, &state, pp, nc, rp, s, &stdf, &hso, &bits, 1, &poly, &v3);\n            printf(\"R %d %d %d %d %d %d\\n\", r, state, stdf, hso, bits, line.len);",
    "int nerr = -9;\n            int r = InChILine2Data(&ios, &line, &pstr, &state, &nerr, pp, nc, s_dummy_placeholder, 1, 1, mode, &stdf, &hso, &bits, &poly, &v3);\n            printf(\"R %d %d %d %d %d %d %d\\n\", r, state, stdf, hso, bits, line.len, nerr);",
)
MAIN = MAIN.replace("s_dummy_placeholder, 1, 1, mode", "rp, s, 1, 1, mode")
# the protons: print the per-component lists too
MAIN = MAIN.replace(
    "printf(\"P%d%d %d %d %d %d\\n\", a, b, rp[a][b].nNumRemovedProtons, rp[a][b].nNumRemovedIsotopicH[0], rp[a][b].nNumRemovedIsotopicH[1], rp[a][b].nNumRemovedIsotopicH[2]);",
    "{ printf(\"P%d%d %d %d %d %d |\", a, b, rp[a][b].nNumRemovedProtons, rp[a][b].nNumRemovedIsotopicH[0], rp[a][b].nNumRemovedIsotopicH[1], rp[a][b].nNumRemovedIsotopicH[2]);\n"
    "                  if (rp[a][b].pNumProtons) { int n2 = nc[a][0] > nc[a][1] ? nc[a][0] : nc[a][1]; for (int q = 0; q < n2; q++) printf(\" %d,%d,%d,%d\", rp[a][b].pNumProtons[q].nNumRemovedProtons, rp[a][b].pNumProtons[q].nNumRemovedIsotopicH[0], rp[a][b].pNumProtons[q].nNumRemovedIsotopicH[1], rp[a][b].pNumProtons[q].nNumRemovedIsotopicH[2]); } else printf(\" NULL\");\n"
    "                  printf(\"\\n\"); }",
)
HARNESS = r2.HARNESS.split("int main(void)")[0] + MAIN

R = random.Random(2468)
hx = r2.hx


def main():
    corpus = [json.loads(l) for l in open(ROOT / "tests/data/inchi_corpus.jsonl")]
    inchis = sorted({c["inchi"] for c in corpus if c.get("inchi", "").startswith("InChI=1") and "/z" not in c["inchi"] and "AuxInfo" not in c["inchi"]})
    r3.R = R
    cases = []
    for _ in range(30000):
        ic = R.choice(inchis)
        if R.random() < 0.5:
            ic = r3.seg_mutate(ic)
        if R.random() < 0.3:
            ic = r3.seg_mutate(ic)
        if R.random() < 0.3:
            ic = r3.mutate_inchi(ic)
        ic += R.choice(["", "", "\n", "\r\n", "\t"])
        mode = R.choice([0, 0, 0, 0x200, 0x400])
        cases.append(f"14 {mode} {hx(ic.encode())}")
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
    path = ROOT / "tests/data/inchi_read4.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
