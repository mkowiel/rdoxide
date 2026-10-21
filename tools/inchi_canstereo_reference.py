"""Generates tests/data/inchi_canstereo.json from the InChI C code: the
CANON_STAT based parts of ichicans.c (`UnmarkNonStereo`,
`FillAllStereoDescriptors`/`FillSingleStereoDescriptors`, `InvertStereo`,
`SwitchAtomStereoAndIsotopicStereo`), linked from the upstream library.

    python tools/inchi_canstereo_reference.py
"""
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"tools/{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


cs = load("inchi_cansparity_reference")
ma = cs.ma
SRC = ma.SRC

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mode.h"
#include "ichicant.h"
#include "ichicomn.h"
#include "ichister.h"
int UnmarkNonStereo( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nRank, const AT_RANK *nAtomNumber, int bIsotopic );
int FillAllStereoDescriptors( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nAtomNumberCanon, CANON_STAT *pCS );
void SwitchAtomStereoAndIsotopicStereo( sp_ATOM *at, int num_atoms, int *bSwitched );
int CreateCompositeNormAtom() { return 0; }
static void dump(sp_ATOM *at, int n) {
    for (int i = 0; i < n; i++) {
        printf(" |%d %d %d %d %d %d %d %d", at[i].parity, at[i].stereo_atom_parity, at[i].final_parity, at[i].bHasStereoOrEquToStereo, at[i].parity2, at[i].stereo_atom_parity2, at[i].final_parity2, at[i].bHasStereoOrEquToStereo2);
        for (int j = 0; j < 3; j++) printf(" %d %d %d %d %d %d %d %d", at[i].stereo_bond_neighbor[j], at[i].stereo_bond_ord[j], at[i].stereo_bond_z_prod[j], at[i].stereo_bond_parity[j], at[i].stereo_bond_neighbor2[j], at[i].stereo_bond_ord2[j], at[i].stereo_bond_z_prod2[j], at[i].stereo_bond_parity2[j]);
    }
    printf("\n");
}
static void dumpcs(CANON_STAT *s) {
    printf("%d %d", s->nLenLinearCTStereoCarb, s->nLenLinearCTStereoDble);
    for (int i = 0; i < s->nLenLinearCTStereoCarb; i++) printf(" c%d,%d", s->LinearCTStereoCarb[i].at_num, s->LinearCTStereoCarb[i].parity);
    for (int i = 0; i < s->nLenLinearCTStereoDble; i++) printf(" d%d,%d,%d", s->LinearCTStereoDble[i].at_num1, s->LinearCTStereoDble[i].at_num2, s->LinearCTStereoDble[i].parity);
    printf("\n");
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    CANON_GLOBALS cg;
    memset(&cg, 0, sizeof(cg));
    for (int c = 0; c < ncases; c++) {
        int n, iso, maxc, maxd, nostereo;
        scanf("%d %d %d %d %d", &n, &iso, &maxc, &maxd, &nostereo);
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        AT_RANK *rank = calloc(n + 1, sizeof(AT_RANK)), *crank = calloc(n + 1, sizeof(AT_RANK));
        AT_RANK *anum = calloc(n + 1, sizeof(AT_RANK)), *acanon = calloc(n + 1, sizeof(AT_RANK));
        for (int i = 0; i < n; i++) {
            int v, p, sap, fp, r, cr, ep, nh, i0, i1, i2, rs, an, ac;
            scanf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d", &v, &p, &sap, &fp, &r, &cr, &ep, &nh, &i0, &i1, &i2, &rs, &an, &ac);
            at[i].valence = v; at[i].parity = p; at[i].stereo_atom_parity = sap; at[i].final_parity = fp;
            at[i].endpoint = ep; at[i].num_H = nh; at[i].num_iso_H[0] = i0; at[i].num_iso_H[1] = i1; at[i].num_iso_H[2] = i2; at[i].nRingSystem = rs;
            at[i].bHasStereoOrEquToStereo = 1;
            rank[i] = r; crank[i] = cr; anum[i] = an; acanon[i] = ac;
            for (int j = 0; j < v; j++) { int q; scanf("%d", &q); at[i].neighbor[j] = q; }
            for (int j = 0; j < 3; j++) { int a, b, z, pp; scanf("%d %d %d %d", &a, &b, &z, &pp); at[i].stereo_bond_neighbor[j] = a; at[i].stereo_bond_ord[j] = b; at[i].stereo_bond_z_prod[j] = z; at[i].stereo_bond_parity[j] = pp; at[i].stereo_bond_neighbor2[j] = pp ? a : 0; at[i].stereo_bond_parity2[j] = pp; }
            at[i].parity2 = p; at[i].stereo_atom_parity2 = sap; at[i].final_parity2 = fp; at[i].bHasStereoOrEquToStereo2 = 0;
        }
        printf("N %d\n", c);
        sp_ATOM *w = malloc((n + 1) * sizeof(sp_ATOM));
        memcpy(w, at, (n + 1) * sizeof(sp_ATOM));
        printf("U %d", UnmarkNonStereo(&cg, w, n, rank, anum, iso)); dump(w, n);
        CANON_STAT s;
        memset(&s, 0, sizeof(s));
        s.nMaxLenLinearCTStereoCarb = maxc; s.nMaxLenLinearCTStereoDble = maxd;
        if (!nostereo) {
            s.LinearCTStereoCarb = calloc(maxc + 1, sizeof(AT_STEREO_CARB));
            s.LinearCTStereoDble = calloc(maxd + 1, sizeof(AT_STEREO_DBLE));
        }
        printf("F %d ", FillAllStereoDescriptors(&cg, w, n, crank, acanon, &s)); dumpcs(&s);
        if (!nostereo) {
            for (int inv = 0; inv < 2; inv++) {
                sp_ATOM *w2 = malloc((n + 1) * sizeof(sp_ATOM));
                memcpy(w2, w, (n + 1) * sizeof(sp_ATOM));
                AT_STEREO_CARB *sc = malloc((maxc + 1) * sizeof(AT_STEREO_CARB));
                AT_STEREO_DBLE *sd = malloc((maxd + 1) * sizeof(AT_STEREO_DBLE));
                memcpy(sc, s.LinearCTStereoCarb, (maxc + 1) * sizeof(AT_STEREO_CARB));
                memcpy(sd, s.LinearCTStereoDble, (maxd + 1) * sizeof(AT_STEREO_DBLE));
                CANON_STAT s2 = s;
                s2.LinearCTStereoCarb = sc; s2.LinearCTStereoDble = sd;
                AT_RANK *an2 = calloc(n + 1, sizeof(AT_RANK));
                printf("I %d", InvertStereo(w2, n, crank, an2, &s2, inv));
                for (int i = 0; i < n; i++) printf(" %d", an2[i]);
                dump(w2, n); dumpcs(&s2);
                free(w2); free(sc); free(sd); free(an2);
            }
            free(s.LinearCTStereoCarb); free(s.LinearCTStereoDble);
        }
        int sw = 0;
        SwitchAtomStereoAndIsotopicStereo(w, n, &sw);
        printf("S %d", sw); dump(w, n);
        free(w); free(at); free(rank); free(crank); free(anum); free(acanon);
    }
    return 0;
}
"""


def gen(rng):
    c = cs.gen(rng)
    atoms = c["atoms"]
    n = len(atoms)
    # equivalence ranks with ties, atoms sorted by rank
    rank = [rng.choice([1, 2, 3, 4, 5, 6, 7, 8][:max(2, n)]) if rng.random() < 0.6 else i + 1 for i in range(n)]
    anum = sorted(range(n), key=lambda i: (rank[i], i))
    perm = list(range(1, n + 1))
    rng.shuffle(perm)
    acanon = sorted(range(n), key=lambda i: perm[i])
    for i, a in enumerate(atoms):
        a["rank"] = rank[i]
        a["crank"] = perm[i]
        a["ep"] = rng.choice([0] * 8 + [1])
        a["nh"] = rng.choice([0, 0, 1, 2, 3])
        a["iso"] = [rng.choice([0, 0, 1, 2]) for _ in range(3)]
        a["rs"] = rng.choice([0, 1, 1, 2, 3])
        a["an"] = anum[i]
        a["ac"] = acanon[i]
    return {"atoms": atoms, "iso": rng.choice([0, 1]), "maxc": rng.choice([0, 1, 2, 8]), "maxd": rng.choice([0, 1, 2, 8]),
            "nostereo": rng.choice([0, 0, 0, 1])}


def fmt(c):
    t = f"{len(c['atoms'])} {c['iso']} {c['maxc']} {c['maxd']} {c['nostereo']}\n"
    for a in c["atoms"]:
        t += f"{len(a['nb'])} {a['par']} {a['sap']} {a['fp']} {a['rank']} {a['crank']} {a['ep']} {a['nh']} {a['iso'][0]} {a['iso'][1]} {a['iso'][2]} {a['rs']} {a['an']} {a['ac']} " + " ".join(map(str, a["nb"])) + "\n"
        t += " ".join(" ".join(map(str, s)) for s in a["sb"]) + "\n"
    return t


def main():
    rng = random.Random(71)
    cases = [gen(rng) for _ in range(5000)]
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True)
        assert r.returncode == 0, r.stderr
        out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_canstereo.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
