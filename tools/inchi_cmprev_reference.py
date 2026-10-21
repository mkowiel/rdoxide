"""Generates tests/data/inchi_cmprev.json from the InChI C code:
`CompareReversedINChI`, `CompareReversedStereoINChI` and
`CompareReversedStereoINChI2` (with the `ICR` lists) of ichimake.c on random
pairs of related INChI structures. Links the upstream InChI base library (see
inchi_markalt_reference.py).

    python tools/inchi_cmprev_reference.py
"""
import copy
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


m2 = load("make2", "tools/inchi_make2_reference.py")
ma = m2.ma
SRC = ma.SRC

HARNESS = m2.HARNESS.split("int main(void)")[0] + r"""
int CompareReversedINChI( INChI *i1, INChI *i2, INChI_Aux *a1, INChI_Aux *a2 );
int CompareReversedStereoINChI( INChI_Stereo *s1, INChI_Stereo *s2 );
int CompareReversedStereoINChI2( INChI_Stereo *s1, INChI_Stereo *s2, ICR *picr );
static INChI_Aux *rd_aux(void) {
    if (!rd()) return NULL;
    INChI_Aux *a = calloc(1, sizeof(*a));
    a->nNumRemovedProtons = rd();
    for (int i = 0; i < 3; i++) a->nNumRemovedIsotopicH[i] = rd();
    return a;
}
static void dump(const char *nm, AT_NUMB *a, int n) {
    printf(" %s", nm);
    for (int i = 0; i < n; i++) printf(",%d", a[i]);
}
int main(void) {
    int ncases = rd();
    for (int k = 0; k < ncases; k++) {
        INChI *c1 = rd_inchi(), *c2 = rd_inchi();
        INChI_Aux *a1 = rd_aux(), *a2 = rd_aux();
        printf("#CASE %d\n", k);
        printf("%d\n", CompareReversedINChI(c1, c2, a1, a2));
        INChI_Stereo *s1 = c1 ? c1->Stereo : NULL, *s2 = c2 ? c2->Stereo : NULL;
        INChI_Stereo *t1 = c1 ? c1->StereoIsotopic : NULL, *t2 = c2 ? c2->StereoIsotopic : NULL;
        printf("%d %d %d %d\n", CompareReversedStereoINChI(s1, s2), CompareReversedStereoINChI(t1, t2),
               CompareReversedStereoINChI(s1, t2), CompareReversedStereoINChI(t1, s2));
        ICR icr;
        memset(&icr, 0, sizeof(icr));
        int r1 = CompareReversedStereoINChI2(s1, s2, &icr);
        int r2 = CompareReversedStereoINChI2(t1, t2, &icr);
        printf("%d %d", r1, r2);
        dump("sc1", icr.sc_in1_only, icr.num_sc_in1_only);
        dump("sc2", icr.sc_in2_only, icr.num_sc_in2_only);
        dump("sb1", icr.sb_in1_only, icr.num_sb_in1_only);
        dump("sb2", icr.sb_in2_only, icr.num_sb_in2_only);
        dump("sbu1", icr.sb_undef_in1_only, icr.num_sb_undef_in1_only);
        dump("sbu2", icr.sb_undef_in2_only, icr.num_sb_undef_in2_only);
        dump("scu1", icr.sc_undef_in1_only, icr.num_sc_undef_in1_only);
        dump("scu2", icr.sc_undef_in2_only, icr.num_sc_undef_in2_only);
        printf("\n");
    }
    return 0;
}
"""


def stereo(rng):
    if rng.random() < 0.35:
        return None
    ns = sorted(rng.sample(range(1, 9), rng.choice([0, 0, 1, 2, 3, 4])))
    bonds = sorted({(rng.randrange(1, 7), rng.randrange(7, 12)) for _ in range(rng.choice([0, 0, 1, 2, 3]))})
    return {
        "sc": [[n, rng.choice([1, 2, 3, 4, 4])] for n in ns],
        "sb": [[a, b, rng.choice([1, 2, 4])] for a, b in bonds],
        "inv": rng.choice([0, 0, 1, -1, 2]),
    }


def mut_stereo(rng, s):
    if s is None:
        return stereo(rng) if rng.random() < 0.5 else None
    s = copy.deepcopy(s)
    r = rng.random()
    if r < 0.2:
        return s
    if r < 0.6 and s["sc"]:
        i = rng.randrange(len(s["sc"]))
        if rng.random() < 0.5:
            del s["sc"][i]
        else:
            s["sc"][i][1] = rng.choice([1, 2, 3, 4])
    elif r < 0.8 and s["sb"]:
        i = rng.randrange(len(s["sb"]))
        if rng.random() < 0.5:
            del s["sb"][i]
        else:
            s["sb"][i][2] = rng.choice([1, 2, 4])
    elif r < 0.9:
        n = rng.randrange(1, 9)
        if all(x[0] != n for x in s["sc"]):
            s["sc"] = sorted(s["sc"] + [[n, rng.choice([1, 2, 4])]])
    else:
        s["inv"] = rng.choice([0, 1, -1, 2])
    return s


def comp(rng):
    c = m2.comp(rng)
    c["st"] = stereo(rng)
    c["sti"] = stereo(rng)
    c["del"] = rng.choice([0, 0, 0, 0, 0, 1])
    c["err"] = 0
    return c


def mutate(rng, c):
    c = copy.deepcopy(c)
    if rng.random() < 0.4:
        c["st"] = mut_stereo(rng, c["st"])
    if rng.random() < 0.5:
        c["sti"] = mut_stereo(rng, c["sti"])
    if rng.random() < 0.15:
        c["sti"] = None
        c["st"] = mut_stereo(rng, c["st"])
    if rng.random() < 0.3:
        c = m2.mutate(rng, c)
        c["err"] = 0
    return c


def aux(rng):
    if rng.random() < 0.3:
        return None
    return [rng.choice([0, 0, 1, 2]), rng.choice([0, 0, 1]), rng.choice([0, 0, 1]), rng.choice([0, 0, 1])]


def tok_aux(a):
    return "0" if a is None else "1 " + " ".join(map(str, a))


def main():
    rng = random.Random(4242)
    cases = []
    for _ in range(4000):
        c1 = comp(rng)
        c2 = mutate(rng, c1)
        # the compared fixed-H lists need numh/fixed of the same length; same n in a mutation
        a1 = aux(rng)
        a2 = a1 if rng.random() < 0.6 else aux(rng)
        miss = rng.random()
        t1 = None if miss < 0.04 else m2.tok(c1)
        t2 = None if 0.04 <= miss < 0.08 else m2.tok(c2)
        cases.append(" ".join([t1 or "0", t2 or "0", tok_aux(a1), tok_aux(a2)]))
    stdin = f"{len(cases)}\n" + "\n".join(cases) + "\n"
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True)
    out = r.stdout.decode().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_cmprev.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
