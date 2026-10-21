"""Generates tests/data/inchi_mak2a.json from the InChI C code: the Hill formula
functions, `Copy2StereoBondOrAllene`, `CopyLinearCTStereoToINChIStereo`,
`MarkAmbiguousStereo` and `UnmarkAllUndefinedUnknownStereo` of ichimak2.c.
Links the upstream InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_mak2a_reference.py
"""
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
char *AllocateAndFillHillFormula( INChI *pINChI );
int GetHillFormulaCounts( U_CHAR *nAtom, S_CHAR *nNum_H, int num_atoms, AT_NUMB *nTautomer, int lenTautomer,
                          int *pnum_C, int *pnum_H, int *pnLen, int *pnNumNonHAtoms );
int CopyLinearCTStereoToINChIStereo( INChI_Stereo *Stereo, AT_STEREO_CARB *LinearCTStereoCarb, int nLenLinearCTStereoCarb,
    AT_STEREO_DBLE *LinearCTStereoDble, int nLenLinearCTStereoDble, AT_NUMB *pCanonOrd, AT_RANK *pCanonRank, sp_ATOM *at,
    int bIsotopic, AT_STEREO_CARB *LinearCTStereoCarbInv, AT_STEREO_DBLE *LinearCTStereoDbleInv, AT_NUMB *pCanonOrdInv, AT_RANK *pCanonRankInv );
int MarkAmbiguousStereo( sp_ATOM *at, inp_ATOM *norm_at, int bIsotopic, AT_NUMB *pCanonOrd, AT_STEREO_CARB *LinearCTStereoCarb,
    int nLenLinearCTStereoCarb, AT_STEREO_DBLE *LinearCTStereoDble, int nLenLinearCTStereoDble );
INCHI_MODE UnmarkAllUndefinedUnknownStereo( INChI_Stereo *Stereo, INCHI_MODE nUserMode );
static int rd(void) { int x; scanf("%d", &x); return x; }
#define CAP 48
static void pr_stereo(INChI_Stereo *s) {
    printf("S %d %d %d %d\n", s->nNumberOfStereoCenters, s->nNumberOfStereoBonds, s->nCompInv2Abs, s->bTrivialInv);
    printf("C");
    for (int i = 0; i < s->nNumberOfStereoCenters; i++) printf(" %d %d %d %d", s->nNumber[i], s->t_parity[i], s->nNumberInv[i], s->t_parityInv[i]);
    printf("\nB");
    for (int i = 0; i < s->nNumberOfStereoBonds; i++) printf(" %d %d %d", s->nBondAtom1[i], s->nBondAtom2[i], s->b_parity[i]);
    printf("\n");
}
int main(void) {
    int ncases = rd();
    for (int k = 0; k < ncases; k++) {
        int mode;
        scanf("%d", &mode);
        printf("#CASE %d\n", k);
        if (mode == 0) {
            INChI c; memset(&c, 0, sizeof(c));
            int n = rd();
            c.nNumberOfAtoms = n;
            c.nAtom = calloc(n + 1, 1); c.nNum_H = calloc(n + 1, 1);
            for (int i = 0; i < n; i++) c.nAtom[i] = rd();
            for (int i = 0; i < n; i++) c.nNum_H[i] = rd();
            int lt = rd();
            c.lenTautomer = lt;
            c.nTautomer = calloc(lt + 2, sizeof(AT_NUMB));
            for (int i = 0; i < lt; i++) c.nTautomer[i] = rd();
            int nc, nh, len, nn;
            int r = GetHillFormulaCounts(c.nAtom, c.nNum_H, n, lt ? c.nTautomer : NULL, lt, &nc, &nh, &len, &nn);
            if (r) printf("E\n"); else printf("%d %d %d %d\n", nc, nh, len, nn);
            char *f = AllocateAndFillHillFormula(&c);
            printf("%s\n", f ? f : "NULL");
        } else {
            int n = rd();
            sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
            inp_ATOM *na = calloc(n + 1, sizeof(inp_ATOM));
            for (int i = 0; i < n; i++) {
                at[i].valence = rd();
                for (int j = 0; j < 3; j++) at[i].neighbor[j] = rd();
                at[i].stereo_bond_parity[0] = rd(); at[i].stereo_bond_parity2[0] = rd();
                at[i].stereo_bond_neighbor[1] = rd(); at[i].stereo_bond_neighbor2[1] = rd();
                at[i].stereo_bond_ord[0] = rd(); at[i].stereo_bond_ord2[0] = rd();
                at[i].bAmbiguousStereo = rd(); na[i].bAmbiguousStereo = at[i].bAmbiguousStereo;
            }
            int iso = rd(), hasord = rd(), hasordinv = rd();
            AT_NUMB *ord = calloc(n + 1, sizeof(AT_NUMB)), *ordi = calloc(n + 1, sizeof(AT_NUMB));
            AT_RANK *rank = calloc(n + 1, sizeof(AT_RANK)), *ranki = calloc(n + 1, sizeof(AT_RANK));
            for (int i = 0; i < n; i++) ord[i] = rd();
            for (int i = 0; i < n; i++) rank[i] = rd();
            for (int i = 0; i < n; i++) ordi[i] = rd();
            for (int i = 0; i < n; i++) ranki[i] = rd();
            int nca = rd(), ndb = rd();
            AT_STEREO_CARB *ca = calloc(nca + 1, sizeof(*ca)), *cai = calloc(nca + 1, sizeof(*cai));
            AT_STEREO_DBLE *db = calloc(ndb + 1, sizeof(*db)), *dbi = calloc(ndb + 1, sizeof(*dbi));
            for (int i = 0; i < nca; i++) { ca[i].at_num = rd(); ca[i].parity = rd(); }
            for (int i = 0; i < nca; i++) { cai[i].at_num = rd(); cai[i].parity = rd(); }
            for (int i = 0; i < ndb; i++) { db[i].at_num1 = rd(); db[i].at_num2 = rd(); db[i].parity = rd(); }
            for (int i = 0; i < ndb; i++) { dbi[i].at_num1 = rd(); dbi[i].at_num2 = rd(); dbi[i].parity = rd(); }
            int umode = rd();
            INChI_Stereo s; memset(&s, 0, sizeof(s));
            s.nNumber = calloc(CAP, sizeof(AT_NUMB)); s.t_parity = calloc(CAP, 1);
            s.nNumberInv = calloc(CAP, sizeof(AT_NUMB)); s.t_parityInv = calloc(CAP, 1);
            s.nBondAtom1 = calloc(CAP, sizeof(AT_NUMB)); s.nBondAtom2 = calloc(CAP, sizeof(AT_NUMB)); s.b_parity = calloc(CAP, 1);
            int r = CopyLinearCTStereoToINChIStereo(&s, ca, nca, db, ndb, hasord ? ord : NULL, hasord ? rank : NULL, at, iso,
                                                    cai, dbi, hasordinv ? ordi : NULL, hasordinv ? ranki : NULL);
            printf("%d\n", r);
            pr_stereo(&s);
            int m = MarkAmbiguousStereo(at, na, iso, hasord ? ord : NULL, ca, nca, db, ndb);
            printf("%d\n", m);
            for (int i = 0; i < n; i++) printf(" %d/%d", at[i].bAmbiguousStereo, na[i].bAmbiguousStereo);
            printf("\n");
            long long u = UnmarkAllUndefinedUnknownStereo(&s, umode);
            printf("%lld\n", u);
            pr_stereo(&s);
        }
    }
    return 0;
}
"""

HILL = [6, 7, 8, 9, 15, 16, 17, 35, 53, 11, 26]  # C N O F P S Cl Br I Na Fe


def sym_key(z):
    names = {6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S", 17: "Cl", 35: "Br", 53: "I", 11: "Na", 26: "Fe", 1: "H"}
    return names.get(z, "?")


def hill_case(rng):
    n = rng.randrange(0, 9)
    atoms = [rng.choice(HILL + [1, 1]) for _ in range(n)]
    r = rng.random()
    if r < 0.8:
        atoms.sort(key=lambda z: (0 if z == 6 else 2 if z == 1 else 1, sym_key(z)))
    elif r < 0.9:
        rng.shuffle(atoms)
    if rng.random() < 0.05 and atoms:
        atoms[rng.randrange(len(atoms))] = rng.choice([0, 200, 119])
    nh = [rng.choice([0, 0, 1, 2, 3]) for _ in range(len(atoms))]
    taut = []
    if rng.random() < 0.4:
        g = rng.randrange(1, 4)
        taut = [g]
        for _ in range(g):
            L = rng.randrange(2, 5)
            taut += [L, rng.randrange(0, 4)] + [rng.randrange(1, 9) for _ in range(L - 1)]
    return f"0 {len(atoms)} " + " ".join(map(str, atoms)) + " " + " ".join(map(str, nh)) + f" {len(taut)} " + " ".join(map(str, taut))


def stereo_case(rng):
    n = rng.randrange(3, 9)
    t = [f"1 {n}"]
    for _ in range(n):
        val = rng.randrange(1, 4)
        nb = [rng.randrange(0, n) for _ in range(3)]
        p1 = rng.choice([0, 8 * rng.randrange(0, 3) + rng.choice([1, 2, 3, 4, 6])])
        p2 = rng.choice([0, 8 * rng.randrange(0, 3) + rng.choice([1, 2, 3, 4, 6])])
        t.append(" ".join(map(str, [val, *nb, p1, p2, rng.choice([0, 0, 0, 1]), rng.choice([0, 0, 0, 1]),
                                     rng.randrange(0, val), rng.randrange(0, val), rng.choice([0, 0, 1, 2])])))
    t.append(f"{rng.choice([0, 1])} {rng.choice([1, 1, 0])} {rng.choice([1, 1, 0])}")
    for _ in range(2):
        perm = list(range(n))
        rng.shuffle(perm)
        t.append(" ".join(map(str, perm)))
        rank = [rng.randrange(1, n + 1) for _ in range(n)]
        t.append(" ".join(map(str, rank)))
    nca = rng.randrange(0, 4)
    nd = rng.randrange(0, 4)
    t.append(f"{nca} {nd}")
    pa = lambda: f"{rng.randrange(1, n + 1)} {rng.choice([1, 2, 2, 3, 4, 5])}"
    pd = lambda: f"{rng.randrange(1, n + 1)} {rng.randrange(1, n + 1)} {rng.choice([1, 2, 2, 3, 4, 5])}"
    base_c = [pa() for _ in range(nca)]
    base_d = [pd() for _ in range(nd)]
    inv_c = [c if rng.random() < 0.9 else pa() for c in base_c]
    inv_d = [d if rng.random() < 0.95 else pd() for d in base_d]
    t += [" ".join(base_c), " ".join(inv_c), " ".join(base_d), " ".join(inv_d)]
    t.append(str(rng.choice([0, 0x800, 0x1000, 0x1800])))
    return "\n".join(t)


def main():
    rng = random.Random(55)
    cases = [hill_case(rng) for _ in range(3000)] + [stereo_case(rng) for _ in range(6000)]
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
    path = ROOT / "tests/data/inchi_mak2a.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
