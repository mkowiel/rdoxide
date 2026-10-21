"""Generates tests/data/inchi_make2.json from the InChI C code: the component
comparison functions of ichimake.c (`CompareTautNonIsoPartOfINChI`,
`CompareInchiStereo`, `GetSp3RelRacAbs`, `CompINChITautVsNonTaut`,
`CompINChI2`, `CompINChINonTaut2`, `CompINChITaut2`) on random pairs of
related INChI structures. Links the upstream InChI base library (see
inchi_markalt_reference.py).

    python tools/inchi_make2_reference.py
"""
import copy
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("markalt", ROOT / "tools/inchi_markalt_reference.py")
ma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ma)
SRC = ma.SRC

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
int CreateCompositeNormAtom() { return 0; }
int CompareTautNonIsoPartOfINChI( const INChI *i1, const INChI *i2 );
int CompareInchiStereo( INChI_Stereo *Stereo1, INCHI_MODE nFlags1, INChI_Stereo *Stereo2, INCHI_MODE nFlags2 );
int GetSp3RelRacAbs( const INChI *pINChI, INChI_Stereo *Stereo );
int CompINChITautVsNonTaut( const INCHI_SORT *p1, const INCHI_SORT *p2, int bCompareIsotopic );
int CompINChI2( const INCHI_SORT *p1, const INCHI_SORT *p2, int bTaut, int bCompareIsotopic );
int CompINChINonTaut2( const void *p1, const void *p2 );
int CompINChITaut2( const void *p1, const void *p2 );
static int sgn(int x) { return x < 0 ? -1 : x > 0; }
static int rd(void) { int x; scanf("%d", &x); return x; }
static INChI_Stereo *rd_stereo(void) {
    if (!rd()) return NULL;
    INChI_Stereo *s = calloc(1, sizeof(*s));
    int n = rd();
    s->nNumberOfStereoCenters = n;
    s->nNumber = calloc(n + 1, sizeof(AT_NUMB)); s->t_parity = calloc(n + 1, 1);
    for (int i = 0; i < n; i++) { s->nNumber[i] = rd(); s->t_parity[i] = rd(); }
    n = rd();
    s->nNumberOfStereoBonds = n;
    s->nBondAtom1 = calloc(n + 1, sizeof(AT_NUMB)); s->nBondAtom2 = calloc(n + 1, sizeof(AT_NUMB)); s->b_parity = calloc(n + 1, 1);
    for (int i = 0; i < n; i++) { s->nBondAtom1[i] = rd(); s->nBondAtom2[i] = rd(); s->b_parity[i] = rd(); }
    s->nCompInv2Abs = rd();
    return s;
}
static INChI *rd_inchi(void) {
    if (!rd()) return NULL;
    INChI *c = calloc(1, sizeof(*c));
    char f[64];
    c->nFlags = rd(); c->nTotalCharge = rd();
    int n = c->nNumberOfAtoms = rd();
    scanf("%s", f); c->szHillFormula = strdup(f);
    c->nAtom = calloc(n + 1, 1);
    for (int i = 0; i < n; i++) c->nAtom[i] = rd();
    int l = c->lenConnTable = rd();
    c->nConnTable = calloc(l + 1, sizeof(AT_NUMB));
    for (int i = 0; i < l; i++) c->nConnTable[i] = rd();
    l = c->lenTautomer = rd();
    c->nTautomer = calloc(l + 1, sizeof(AT_NUMB));
    for (int i = 0; i < l; i++) c->nTautomer[i] = rd();
    c->nNum_H = calloc(n + 1, 1);
    for (int i = 0; i < n; i++) c->nNum_H[i] = rd();
    if (rd()) { c->nNum_H_fixed = calloc(n + 1, 1); for (int i = 0; i < n; i++) c->nNum_H_fixed[i] = rd(); }
    l = c->nNumberOfIsotopicAtoms = rd();
    c->IsotopicAtom = calloc(l + 1, sizeof(INChI_IsotopicAtom));
    for (int i = 0; i < l; i++) { c->IsotopicAtom[i].nAtomNumber = rd(); c->IsotopicAtom[i].nIsoDifference = rd();
        c->IsotopicAtom[i].nNum_H = rd(); c->IsotopicAtom[i].nNum_D = rd(); c->IsotopicAtom[i].nNum_T = rd(); }
    l = c->nNumberOfIsotopicTGroups = rd();
    c->IsotopicTGroup = calloc(l + 1, sizeof(INChI_IsotopicTGroup));
    for (int i = 0; i < l; i++) { c->IsotopicTGroup[i].nTGroupNumber = rd(); c->IsotopicTGroup[i].nNum_H = rd();
        c->IsotopicTGroup[i].nNum_D = rd(); c->IsotopicTGroup[i].nNum_T = rd(); }
    c->Stereo = rd_stereo();
    c->StereoIsotopic = rd_stereo();
    c->bDeleted = rd();
    return c;
}
int main(void) {
    int ncases = rd();
    for (int k = 0; k < ncases; k++) {
        INCHI_SORT p[2];
        memset(p, 0, sizeof(p));
        for (int q = 0; q < 2; q++) {
            p[q].pINChI[0] = rd_inchi(); p[q].pINChI[1] = rd_inchi(); p[q].ord_number = q;
        }
        printf("#CASE %d\n", k);
        INChI *a = p[0].pINChI[1] ? p[0].pINChI[1] : p[0].pINChI[0];
        INChI *b = p[1].pINChI[1] ? p[1].pINChI[1] : p[1].pINChI[0];
        if (a && b) printf("%d\n", sgn(CompareTautNonIsoPartOfINChI(a, b)));
        else printf("0\n");
        if (a) printf("%d %d\n", GetSp3RelRacAbs(a, a->Stereo), GetSp3RelRacAbs(a, a->StereoIsotopic));
        else printf("0 0\n");
        if (a && b) printf("%d %d\n", sgn(CompareInchiStereo(a->Stereo, a->nFlags, b->Stereo, b->nFlags)),
                           sgn(CompareInchiStereo(a->StereoIsotopic, a->nFlags, b->Stereo, b->nFlags)));
        else printf("0 0\n");
        printf("%d %d\n", sgn(CompINChITautVsNonTaut(&p[0], &p[1], 1)), sgn(CompINChITautVsNonTaut(&p[0], &p[1], 0)));
        printf("%d %d %d %d\n", sgn(CompINChI2(&p[0], &p[1], 0, 1)), sgn(CompINChI2(&p[0], &p[1], 1, 1)),
               sgn(CompINChI2(&p[0], &p[1], 0, 0)), sgn(CompINChI2(&p[0], &p[1], 1, 0)));
        printf("%d %d\n", sgn(CompINChINonTaut2(&p[0], &p[1])), sgn(CompINChITaut2(&p[0], &p[1])));
    }
    return 0;
}
"""

FORMS = ["C2H6O", "C2H6", "CH4", "C6H6", "ClH", "C3H8O2", "C2", "H2O", "Br", "C10H22N", "CH2O"]


def stereo(rng):
    if rng.random() < 0.4:
        return None
    return {
        "sc": [[rng.randrange(1, 6), rng.choice([1, 2, 3, 4])] for _ in range(rng.choice([0, 0, 1, 2, 3]))],
        "sb": [[rng.randrange(1, 6), rng.randrange(1, 6), rng.choice([1, 2])] for _ in range(rng.choice([0, 0, 1, 2]))],
        "inv": rng.choice([0, 0, 1, -1, 2]),
    }


def comp(rng):
    n = rng.randrange(1, 5)
    return {
        "flags": rng.choice([0, 0, 0, 2, 4, 1]),
        "charge": rng.choice([0, 0, 0, 1, -1, 2]),
        "n": n,
        "f": rng.choice(FORMS),
        "atoms": [rng.choice([6, 6, 7, 8]) for _ in range(n)],
        "ct": [rng.randrange(0, 5) for _ in range(rng.choice([0, 1, 3, 5]))],
        "taut": rng.choice([[], [], [1, 2, 3], [2, 1, 1]]),
        "numh": [rng.choice([0, 0, 1, 2, 3]) for _ in range(n)],
        "fixed": rng.choice([None, None, [rng.choice([0, 0, 1]) for _ in range(n)]]),
        "iso": [[rng.randrange(1, n + 1), rng.choice([0, 1, -1]), rng.choice([0, 1]), rng.choice([0, 1]), rng.choice([0, 1])]
                for _ in range(rng.choice([0, 0, 1, 2]))],
        "isotg": [[rng.randrange(1, 3)] + [rng.choice([0, 1]) for _ in range(3)] for _ in range(rng.choice([0, 0, 1]))],
        "st": stereo(rng),
        "sti": stereo(rng),
        "del": rng.choice([0, 0, 0, 0, 1]),
    }


def mutate(rng, c):
    c = copy.deepcopy(c)
    r = rng.random()
    if r < 0.3:
        return c
    k = rng.choice(["flags", "charge", "f", "ct", "taut", "numh", "fixed", "iso", "isotg", "st", "sti", "del", "atoms"])
    new = comp(rng)
    if k in ("f",):
        c["f"] = new["f"]
    elif k == "ct":
        c["ct"] = new["ct"]
    elif k == "numh":
        c["numh"] = [rng.choice([0, 1, 2, 3]) for _ in range(c["n"])]
    elif k == "fixed":
        c["fixed"] = None if c["fixed"] else [rng.choice([0, 1]) for _ in range(c["n"])]
    elif k == "atoms":
        c["atoms"] = [rng.choice([6, 7, 8]) for _ in range(c["n"])]
    else:
        c[k] = new[k]
    return c


def tok_stereo(s):
    if s is None:
        return "0"
    t = f"1 {len(s['sc'])} " + " ".join(f"{a} {b}" for a, b in s["sc"])
    t += f" {len(s['sb'])} " + " ".join(f"{a} {b} {c}" for a, b, c in s["sb"])
    return t + f" {s['inv']}"


def tok(c):
    if c is None:
        return "0"
    t = f"1 {c['flags']} {c['charge']} {c['n']} {c['f']} " + " ".join(map(str, c["atoms"]))
    t += f" {len(c['ct'])} " + " ".join(map(str, c["ct"]))
    t += f" {len(c['taut'])} " + " ".join(map(str, c["taut"]))
    t += " " + " ".join(map(str, c["numh"]))
    t += " 0" if c["fixed"] is None else " 1 " + " ".join(map(str, c["fixed"]))
    t += f" {len(c['iso'])} " + " ".join(" ".join(map(str, x)) for x in c["iso"])
    t += f" {len(c['isotg'])} " + " ".join(" ".join(map(str, x)) for x in c["isotg"])
    t += " " + tok_stereo(c["st"]) + " " + tok_stereo(c["sti"]) + f" {c['del']}"
    return t


def main():
    rng = random.Random(21)
    cases = []
    for _ in range(8000):
        base = comp(rng)
        slots = []
        for _p in range(2):
            has_t = rng.random() < 0.6
            has_n = rng.random() < 0.7
            slots.append(mutate(rng, base) if has_t else None)
            slots.append(mutate(rng, base) if has_n else None)
        cases.append(" ".join(tok(s) for s in slots))
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
    path = ROOT / "tests/data/inchi_make2.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
