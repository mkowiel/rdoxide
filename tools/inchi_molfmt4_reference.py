"""Generates tests/data/inchi_molfmt4.json from the InChI C code: `ReadMolfile`
(the Molfile reader plus `SDFileSkipExtraData`, `SDFileIdentifyLabel`,
`SDFileExtractCASNo` and the pseudo-element treatment) on V2000/V3000 molfiles
followed by random SD-file data items. Links the upstream InChI base library
(see inchi_markalt_reference.py).

    python tools/inchi_molfmt4_reference.py
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


m1 = load("molfmt1", "tools/inchi_molfmt1_reference.py")
m2 = load("molfmt2", "tools/inchi_molfmt2_reference.py")
m3 = load("molfmt3", "tools/inchi_molfmt3_reference.py")
ma = m1.ma
SRC = ma.SRC

HARNESS = m3.HARNESS.split("int main(void) {")[0] + r"""
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int bcoord, polymers, npzz, nowarn, use_name, use_id, use_val, label_kind;
        printf("#CASE %d\n", cs);
        scanf("%d %d %d %d %d %d %d %d %s", &bcoord, &polymers, &npzz, &nowarn, &use_name, &use_id, &use_val, &label_kind, hx);
        int n = unhex(hx, text);
        INCHI_IOSTREAM ios;
        inchi_ios_init(&ios, INCHI_IOS_TYPE_STRING, NULL);
        ios.s.pStr = (char *) malloc(n + 2);
        memcpy(ios.s.pStr, text, n + 1);
        ios.s.nUsedLength = n; ios.s.nAllocatedLength = n + 2; ios.s.nPtr = 0;
        char err[STR_ERR_LEN]; err[0] = 0;
        char name[100]; memset(name, 0, sizeof(name));
        char value[300]; memset(value, 0, sizeof(value));
        unsigned long id = 777;
        const char *labels[] = { NULL, "", "CAS", "CAS_NO", "MyLabel", "NAME", "COMMENT", "MolfileName", "MolfileLine2", "MolfileComment", "MolfileIntRegNo", "x" };
        int e = 0;
        MOL_FMT_DATA *m = ReadMolfile(&ios, NULL, NULL, bcoord, polymers, npzz, use_name ? name : NULL, use_name ? 64 : 0, use_id ? &id : NULL,
                                      labels[label_kind], use_val ? value : NULL, &e, err, nowarn);
        if (m) {
            dump_hdr(&m->hdr); dump_ctab(&m->ctab);
            dump_all(&m->ctab);
        }
        printf("N "); hexs(name); printf(" %lu ", id); hexs(value); printf("\n");
        printf("R %d ", e); hexs(err); printf("\n");
        inchi_ios_free_str(&ios);
    }
    return 0;
}
"""

R = random.Random(4444)


def sdf_tail(r):
    L = []
    if r.random() < 0.3:
        L.append("M  END")
    for _ in range(r.randrange(0, 5)):
        hdr = r.choice(["> <NAME>", ">  <NAME> 1", "> <CAS_NO>", "> <CAS>", "> <MyLabel>", "> <COMMENT>", "> <other>", ">  <MyLabel>  (1)", "> <CAS RN>", ">", "> <x y>", "> <NAME >", "> <x>", "> <x>", "> <x>"])
        L.append(hdr)
        for _ in range(r.randrange(0, 3)):
            L.append(r.choice(["ethanol", "64-17-5", "  spaced   out  ", "CAS 50-00-0 x", "-12-3", "", "9" * 25, "x" * r.randrange(1, 260), "a\tb"]))
        if r.random() < 0.8:
            L.append("")
    if r.random() < 0.15:
        L.append(r.choice(["junk line", "$$$ ", "$$$$ x"]))
        L.append("")
    if r.random() < 0.9:
        L.append("$$$$")
    if r.random() < 0.2:
        L += ["molB", "second record"]
    return L


def main():
    cases = []
    for _ in range(2300):
        body = m2.mk_mol2(R) if R.random() < 0.65 else m3.mk_mol3(R)
        nl = "\r\n" if "\r\n" in body else "\n"
        text = body + nl.join(sdf_tail(R)) + nl
        cases.append(f"{R.choice([0, 1])} {R.choice([0, 1, 1])} {R.choice([0, 1, 2])} {R.choice([0, 1])} {R.choice([0, 1])} {R.choice([0, 1])} "
                     f"{R.choice([0, 1])} {R.randrange(0, 12)} {text.encode().hex()}")
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-ffp-contract=off", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        kept, out = [], []
        for c in cases:
            try:
                r1 = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=10)
            except subprocess.TimeoutExpired:
                continue
            if r1.returncode == 0:
                # the C code reads freed memory in some inputs: keep only reproducible outputs
                again = [subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=10) for _ in range(4)]
                if any(a.returncode != 0 or a.stdout != r1.stdout for a in again):
                    continue
                kept.append(c)
                out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode(errors="replace").split("\n")[1:]
        print(len(cases) - len(kept), "cases dropped")
        cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    data = [x for x in data if all(len(l.split()[2]) <= 62 for l in x["out"] if l.startswith("O "))]
    path = ROOT / "tests/data/inchi_molfmt4.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
