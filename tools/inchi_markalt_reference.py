"""Generates tests/data/inchi_markalt.json from the InChI C code:
`mark_alt_bonds_and_taut_groups` of ichi_bns.c, the main normalization
driver (protonation, charge/tautomeric/salt groups, alternating bonds), on
random molecules. Unlike the other harnesses it links the whole upstream
InChI library (compiled with `cc` into a cache directory).

    python tools/inchi_markalt_reference.py
"""
import importlib.util
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src"
CACHE = Path(tempfile.gettempdir()) / "rdoxide_inchi_lib"

spec = importlib.util.spec_from_file_location("marktaut", ROOT / "tools/inchi_marktaut_reference.py")
mt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mt)

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mode.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichi_bns.h"
#include "ichicano.h"
#include "ichicant.h"
#include "ichinorm.h"

int CreateCompositeNormAtom() { return 0; } /* referenced, never called */
static const char *elsym(int el) {
    switch (el) { case 1: return "H"; case 6: return "C"; case 7: return "N"; case 8: return "O"; case 15: return "P"; case 16: return "S"; case 17: return "Cl"; case 35: return "Br"; default: return "Xx"; }
}
static void dumpat(const char *tag, inp_ATOM *at, int n) {
    printf("%s", tag);
    for (int i = 0; i < n; i++) {
        printf(" %d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d:", at[i].charge, at[i].num_H, at[i].radical, at[i].chem_bonds_valence, at[i].endpoint, at[i].c_point, at[i].at_type, at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2], at[i].cFlags);
        for (int j = 0; j < at[i].valence; j++) printf("%d.%d.", at[i].neighbor[j], at[i].bond_type[j]);
    }
    printf("\n");
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, nrem, use_tgi;
        unsigned long long flags, done;
        scanf("%d %d %llu %llu %d", &n, &nrem, &flags, &done, &use_tgi);
        inp_ATOM *at = calloc(n + 8, sizeof(inp_ATOM));
        inp_ATOM *out = calloc(n + 8, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, chg, nh, i0, i1, i2, orig, ring, cv, rsys, att;
            scanf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d", &v, &cbv, &el, &rad, &chg, &nh, &i0, &i1, &i2, &orig, &ring, &cv, &rsys, &att);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; strcpy(at[i].elname, elsym(el)); at[i].radical = rad;
            at[i].charge = chg; at[i].num_H = nh; at[i].num_iso_H[0] = i0; at[i].num_iso_H[1] = i1; at[i].num_iso_H[2] = i2; at[i].orig_at_number = orig;
            at[i].nNumAtInRingSystem = ring; at[i].bCutVertex = cv; at[i].nRingSystem = rsys; at[i].at_type = att;
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        cg.m_num_bit = 16;
        cg.m_bBit = calloc(16, sizeof(bitWord));
        for (int i = 0; i < 16; i++) cg.m_bBit[i] = 1 << i;
        T_GROUP_INFO tgi;
        memset(&tgi, 0, sizeof(tgi));
        tgi.bTautFlags = flags;
        tgi.bTautFlagsDone = done;
        tgi.tni.nNumRemovedExplicitH = nrem;
        INCHI_MODE inpflags = flags, inpdone = done;
        printf("N %d\n", c);
        int ret = mark_alt_bonds_and_taut_groups(NULL, &cg, at, out, n, NULL, use_tgi ? &tgi : NULL, &inpflags, &inpdone, 0, NULL);
        printf("R %d\n", ret);
        dumpat("A", at, n);
        dumpat("O", out, n + nrem > n + 8 ? n + 8 : n + nrem);
        printf("F %llu %llu\n", (unsigned long long)(use_tgi ? tgi.bTautFlagsDone : inpdone), (unsigned long long)(use_tgi ? tgi.bTautFlags : inpflags));
        if (use_tgi) {
            printf("I %d %d %d %d %d %d %d %d %d |", tgi.tni.nNumRemovedExplicitH, tgi.tni.nNumRemovedProtons, tgi.tni.nNumRemovedProtonsIsotopic[0], tgi.tni.nNumRemovedProtonsIsotopic[1], tgi.tni.nNumRemovedProtonsIsotopic[2], (int)tgi.tni.bNormalizationFlags, tgi.num_iso_H[0], tgi.num_iso_H[1], tgi.num_iso_H[2]);
            printf(" %d |", tgi.nNumIsotopicEndpoints);
            if (tgi.nIsotopicEndpointAtomNumber) for (int i = 0; i <= tgi.nIsotopicEndpointAtomNumber[0]; i++) printf(" %d", tgi.nIsotopicEndpointAtomNumber[i]);
            printf("\n");
            printf("G %d\n", tgi.num_t_groups);
            for (int i = 0; i < tgi.num_t_groups; i++) {
                printf("T %d %d |", tgi.t_group[i].nGroupNumber, tgi.t_group[i].nNumEndpoints);
                for (int q = 0; q < 5; q++) printf(" %d", tgi.t_group[i].num[q]);
                printf(" |");
                for (int q = 0; q < 6; q++) printf(" %d", tgi.t_group[i].num_DA[q]);
                printf("\n");
            }
        }
        free(at); free(out);
    }
    return 0;
}
"""


def build_lib():
    CACHE.mkdir(exist_ok=True)
    jobs = []
    for f in sorted(SRC.glob("*.c")):
        if f.name == "inchi_gui.c":
            continue
        o = CACHE / (f.stem + ".o")
        if not o.exists() or o.stat().st_mtime < f.stat().st_mtime:
            jobs.append(subprocess.Popen(["cc", "-O1", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-c", str(f), "-o", str(o)]))
    for j in jobs:
        assert j.wait() == 0
    lib = CACHE / "libinchi_base.a"
    objs = sorted(str(p) for p in CACHE.glob("*.o"))
    if not lib.exists() or any(Path(o).stat().st_mtime > lib.stat().st_mtime for o in objs):
        lib.unlink(missing_ok=True) if hasattr(lib, "unlink") else None
        subprocess.run(["ar", "rcs", str(lib), *objs], check=True)
    return [str(lib)]


FLAGS = [0x1, 0x4, 0x8, 0x10, 0x20, 0x40, 0x400, 0x800, 0x1000, 0x80000, 0x100000, 0x800000, 0x1000000, 0x2000000, 0x4000000, 0x8000000, 0x10000000]


def fmt(c):
    t = f"{len(c['atoms'])} {c['nrem']} {c['flags']} {c['done']} {c['use_tgi']}\n"
    for a in c["atoms"]:
        t += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['chg']} {a['nh']} {a['iso'][0]} {a['iso'][1]} {a['iso'][2]} {a['orig']} {a['ring']} {a['cv']} {a['rs']} 0\n"
        t += " ".join(f"{x} {t2}" for x, t2 in a["nb"]) + "\n"
    return t


def main():
    rng = random.Random(77)
    cases = []
    for _ in range(1500):
        atoms = mt.gen(rng)
        for i, a in enumerate(atoms):
            a["iso"] = [rng.choice([0, 0, 0, 1]) for _ in range(3)]
            a["orig"] = i + 1
            a["ring"], a["cv"] = (a["mk"][0], a["mk"][1]) if a.get("mk") else (rng.choice([0, 1, 5, 6]), 0)
            a["rs"] = a["mk"][2] if a.get("mk") else i % 2
            if rng.random() < 0.02:
                a["rad"] = rng.choice([1, 2, 3])
        r = rng.random()
        if r < 0.12:
            # a pair of adjacent radical centers
            cand = [(i, j) for i, a in enumerate(atoms) for j, t in a["nb"] if i < j and a["el"] in (6, 7, 8) and atoms[j]["el"] in (6, 7, 8)]
            if cand:
                i, j = rng.choice(cand)
                atoms[i]["rad"] = atoms[j]["rad"] = 2
                for k in (i, j):
                    atoms[k]["nh"] = max(0, atoms[k]["nh"] - 1)
        elif r < 0.18:
            # an aromatic-like radical: valence 2, both bonds alternating
            cand = [i for i, a in enumerate(atoms) if a["v"] == 2 and a["el"] == 6]
            if cand:
                i = rng.choice(cand)
                atoms[i]["rad"] = 2
                for j, t in atoms[i]["nb"]:
                    pass
                atoms[i]["nb"] = [(j, 4) for j, t in atoms[i]["nb"]]
                for j, t in atoms[i]["nb"]:
                    atoms[j]["nb"] = [(q, 4 if q == i else t2) for q, t2 in atoms[j]["nb"]]
        flags = 0
        for f in FLAGS:
            if rng.random() < 0.45:
                flags |= f
        if rng.random() < 0.8:
            flags |= 1
        done = rng.choice([0, 0, 0x2000, 0x4000, 0x8000, 0xc000, 0x400])
        cases.append({"atoms": atoms, "nrem": rng.choice([0, 0, 0, 1, 2]), "flags": flags, "done": done,
                      "use_tgi": 0 if rng.random() < 0.08 else 1})
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True)
        if r.returncode == 0:
            out = r.stdout.decode().strip().split("\n")
        else:
            kept, lines = [], []
            for c in cases:
                r1 = subprocess.run([str(exe)], input=("1\n" + fmt(c)).encode(), capture_output=True)
                if r1.returncode == 0:
                    kept.append(c)
                    lines += [f"N {len(kept) - 1}"] + r1.stdout.decode().strip().split("\n")[1:]
            print(len(cases) - len(kept), "cases dropped")
            cases = kept
            out = lines
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_markalt.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
