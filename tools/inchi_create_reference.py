"""Generates tests/data/inchi_create.json from the InChI C code: `Create_INChI`
of ichimake.c (normalization, tautomer groups, stereo, canonicalization and
filling of INChI/INChI_Aux for one component) on random small molecules with
3D coordinates, isotopes, charges and implicit/isotopic H. Links the upstream
InChI base library (see inchi_markalt_reference.py).

    python tools/inchi_create_reference.py
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
#include "ichimake.h"
#include "ichitime.h"
int CreateCompositeNormAtom() { return 0; }
int Create_INChI( CANON_GLOBALS *pCG, INCHI_CLOCK *ic, INPUT_PARMS *ip, INChI **ppINChI, INChI_Aux **ppINChI_Aux,
                  ORIG_ATOM_DATA *orig_inp_data, inp_ATOM *inp_at, INP_ATOM_DATA *out_norm_data[2], int num_inp_at,
                  INCHI_MODE nUserMode, INCHI_MODE *pbTautFlags, INCHI_MODE *pbTautFlagsDone, struct tagInchiTime *ulMaxTime,
                  T_GROUP_INFO *ti_out, char *pStrErrStruct );
static const char *elsym(int el) {
    switch (el) { case 6: return "C"; case 7: return "N"; case 8: return "O"; case 15: return "P"; case 16: return "S"; case 17: return "Cl"; default: return "Xx"; }
}
static void pr_stereo(const char *tag, INChI_Stereo *s) {
    if (!s) { printf("%s -\n", tag); return; }
    printf("%s %d %d %d %d |", tag, s->nNumberOfStereoCenters, s->nNumberOfStereoBonds, s->nCompInv2Abs, s->bTrivialInv);
    for (int i = 0; i < s->nNumberOfStereoCenters; i++) printf(" %d,%d,%d,%d", s->nNumber[i], s->t_parity[i], s->nNumberInv[i], s->t_parityInv[i]);
    printf(" |");
    for (int i = 0; i < s->nNumberOfStereoBonds; i++) printf(" %d,%d,%d", s->nBondAtom1[i], s->nBondAtom2[i], s->b_parity[i]);
    printf("\n");
}
static void pr_arr(const char *tag, AT_NUMB *a, int n) {
    printf("%s", tag);
    if (!a) { printf(" -\n"); return; }
    for (int i = 0; i < n; i++) printf(" %d", a[i]);
    printf("\n");
}
int main(void) {
    int ncases;
    setvbuf(stdout, NULL, _IONBF, 0);
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, which;
        unsigned long long mode, flags, done;
        int loose;
        scanf("%d %d %llu %llu %llu %d", &n, &which, &mode, &flags, &done, &loose);
        inp_ATOM *at = calloc(n + 8, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, chg, nh, i0, i1, i2, aw;
            double x, y, z;
            scanf("%d %d %d %d %d %d %d %d %d %d %lf %lf %lf", &v, &cbv, &el, &rad, &chg, &nh, &i0, &i1, &i2, &aw, &x, &y, &z);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; strcpy(at[i].elname, elsym(el));
            at[i].radical = rad; at[i].charge = chg; at[i].num_H = nh;
            at[i].num_iso_H[0] = i0; at[i].num_iso_H[1] = i1; at[i].num_iso_H[2] = i2; at[i].iso_atw_diff = aw;
            at[i].orig_at_number = i + 1; at[i].x = x; at[i].y = y; at[i].z = z;
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        cg.m_num_bit = 16;
        cg.m_bBit = calloc(16, sizeof(bitWord));
        for (int i = 0; i < 16; i++) cg.m_bBit[i] = 1 << i;
        INCHI_CLOCK ic; memset(&ic, 0, sizeof(ic));
        INPUT_PARMS ip; memset(&ip, 0, sizeof(ip));
        ip.bLooseTSACheck = loose;
        INP_ATOM_DATA nd[2];
        memset(nd, 0, sizeof(nd));
        INP_ATOM_DATA *pnd[2] = { &nd[0], &nd[1] };
        if (which & 1) { nd[0].at = calloc(n + 8, sizeof(inp_ATOM)); nd[0].num_at = n; }
        if (which & 2) { nd[1].at = calloc(n + 8, sizeof(inp_ATOM)); nd[1].num_at = n; if (which & 4) nd[1].at_fixed_bonds = calloc(n + 8, sizeof(inp_ATOM)); }
        INChI *pi[2] = {NULL, NULL};
        INChI_Aux *pa[2] = {NULL, NULL};
        int nb, ni;
        for (int k = 0; k < 2; k++) {
            if (!(which & (1 << k))) continue;
            pi[k] = Alloc_INChI(at, n, &nb, &ni, mode);
            pa[k] = Alloc_INChI_Aux(n, ni, mode, 0);
        }
        T_GROUP_INFO ti; memset(&ti, 0, sizeof(ti));
        char se[STR_ERR_LEN]; se[0] = 0;
        INCHI_MODE fl = flags, dn = done;
        printf("#CASE %d\n", c);
        int ret = Create_INChI(&cg, &ic, &ip, pi, pa, NULL, at, pnd, n, mode, &fl, &dn, NULL, &ti, se);
        printf("R %d %llu %llu %d %d\n", ret, (unsigned long long)fl, (unsigned long long)dn, ti.num_t_groups, ti.nNumEndpoints);
        printf("E %s\n", se);
        for (int k = 0; k < 2; k++) {
            INP_ATOM_DATA *d = &nd[k];
            printf("D%d %d %d %d %d %d %d %d %d %d %d %d %llu %llu %llu\n", k, d->num_at, d->num_removed_H, d->bDeleted, d->bTautomeric, d->bTautPreprocessed,
                   d->nNumRemovedProtons, d->nNumRemovedProtonsIsotopic[0], d->nNumRemovedProtonsIsotopic[1], d->nNumRemovedProtonsIsotopic[2],
                   d->num_iso_H[0], d->num_iso_H[1] * 100 + d->num_iso_H[2], (unsigned long long)d->bTautFlags, (unsigned long long)d->bTautFlagsDone, (unsigned long long)d->bNormalizationFlags);
            if (d->at) {
                printf("DA%d", k);
                for (int i = 0; i < n; i++) {
                    inp_ATOM *a = &d->at[i];
                    printf(" %d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d:", a->charge, a->num_H, a->radical, a->chem_bonds_valence, a->endpoint, a->c_point, a->at_type, a->num_iso_H[0], a->num_iso_H[1], a->num_iso_H[2], a->valence, a->orig_at_number);
                    for (int j = 0; j < a->valence; j++) printf("%d.%d.", a->neighbor[j], a->bond_type[j]);
                }
                printf("\n");
            }
            INChI *p = pi[k]; INChI_Aux *a = pa[k];
            if (!p) continue;
            printf("I%d %d %lu %d %d %s\n", k, p->nErrorCode, (unsigned long)p->nFlags, p->nTotalCharge, p->nNumberOfAtoms, p->szHillFormula ? p->szHillFormula : "NULL");
            printf("IA%d", k); for (int i = 0; i < p->nNumberOfAtoms; i++) printf(" %d", p->nAtom[i]); printf("\n");
            printf("IC%d %d", k, p->lenConnTable); for (int i = 0; i < p->lenConnTable; i++) printf(" %d", p->nConnTable[i]); printf("\n");
            printf("IT%d %d", k, p->lenTautomer); for (int i = 0; i < p->lenTautomer; i++) printf(" %d", p->nTautomer[i]); printf("\n");
            printf("IH%d", k); for (int i = 0; i < p->nNumberOfAtoms; i++) printf(" %d/%d", p->nNum_H[i], p->nNum_H_fixed ? p->nNum_H_fixed[i] : -1); printf("\n");
            printf("II%d %d %d", k, p->nNumberOfIsotopicAtoms, p->nNumberOfIsotopicTGroups);
            for (int i = 0; i < p->nNumberOfIsotopicAtoms; i++) printf(" %d,%d,%d,%d,%d", p->IsotopicAtom[i].nAtomNumber, p->IsotopicAtom[i].nIsoDifference, p->IsotopicAtom[i].nNum_H, p->IsotopicAtom[i].nNum_D, p->IsotopicAtom[i].nNum_T);
            printf(" |");
            for (int i = 0; i < p->nNumberOfIsotopicTGroups; i++) printf(" %d,%d,%d,%d", p->IsotopicTGroup[i].nTGroupNumber, p->IsotopicTGroup[i].nNum_H, p->IsotopicTGroup[i].nNum_D, p->IsotopicTGroup[i].nNum_T);
            printf("\n");
            pr_arr("IL", p->nPossibleLocationsOfIsotopicH, p->nPossibleLocationsOfIsotopicH ? p->nPossibleLocationsOfIsotopicH[0] : 0);
            pr_stereo("IS", p->Stereo);
            pr_stereo("IZ", p->StereoIsotopic);
            printf("X%d %d %d %d %d %d\n", k, a->nErrorCode, a->nNumberOfAtoms, a->nNumberOfTGroups, a->bIsTautomeric, a->nCanonFlags);
            pr_arr("XA", a->nOrigAtNosInCanonOrd, n); pr_arr("XB", a->nOrigAtNosInCanonOrdInv, n);
            pr_arr("XC", a->nIsotopicOrigAtNosInCanonOrd, n); pr_arr("XD", a->nIsotopicOrigAtNosInCanonOrdInv, n);
            pr_arr("XE", a->nConstitEquNumbers, n); pr_arr("XF", a->nConstitEquIsotopicNumbers, n);
            pr_arr("XG", a->nConstitEquTGroupNumbers, a->nNumberOfTGroups); pr_arr("XH", a->nConstitEquIsotopicTGroupNumbers, a->nNumberOfTGroups);
            printf("XO"); if (a->OrigInfo) for (int i = 0; i < n; i++) printf(" %d,%d,%d", a->OrigInfo[i].cCharge, a->OrigInfo[i].cRadical, a->OrigInfo[i].cUnusualValence); printf("\n");
            printf("XP %d %d %d %d\n", a->nNumRemovedProtons, a->nNumRemovedIsotopicH[0], a->nNumRemovedIsotopicH[1], a->nNumRemovedIsotopicH[2]);
        }
    }
    return 0;
}
"""

MAXV = {6: 4, 7: 3, 8: 2, 15: 3, 16: 2, 17: 1}


def gen(rng):
    n = rng.randrange(1, 10)
    els = [rng.choice([6, 6, 6, 7, 7, 8, 8, 8, 16, 15, 17]) for _ in range(n)]
    edges = []
    for v in range(1, n):
        if rng.random() < 0.96:
            edges.append((rng.randrange(0, v), v))
    for _ in range(rng.choice([0, 0, 1, 1, 2])):
        if n > 2:
            a, b = rng.sample(range(n), 2)
            if (a, b) not in edges and (b, a) not in edges:
                edges.append((a, b))
    cap = [MAXV[e] for e in els]
    bonds = []
    nb = [[] for _ in range(n)]
    for a, b in edges:
        if cap[a] < 1 or cap[b] < 1:
            continue
        order = rng.choices([1, 2, 3], [70, 26, 4])[0]
        order = min(order, cap[a], cap[b])
        cap[a] -= order
        cap[b] -= order
        nb[a].append((b, order))
        nb[b].append((a, order))
    atoms = []
    for i in range(n):
        cbv = sum(t for _, t in nb[i])
        el = els[i]
        chg = 0
        nh = cap[i]
        r = rng.random()
        if r < 0.1 and el == 7:
            chg, nh = 1, nh + 1
        elif r < 0.2 and el == 8 and nh >= 1:
            chg, nh = -1, nh - 1
        elif r < 0.25 and nh >= 1 and el in (6, 7):
            nh -= 1
        iso = [0, 0, 0]
        if nh > 0 and rng.random() < 0.15:
            for k in range(3):
                if nh > 0 and rng.random() < 0.5:
                    iso[k] = 1
                    nh -= 1
        aw = rng.choice([0, 0, 0, 0, 0, 1, -1, 2]) if rng.random() < 0.2 else 0
        xyz = [round(rng.uniform(-2, 2), 3) for _ in range(3)]
        atoms.append({"v": len(nb[i]), "cbv": cbv, "el": el, "rad": 0, "chg": chg, "nh": nh, "iso": iso, "aw": aw,
                      "xyz": xyz, "nb": nb[i]})
    return atoms


FLAGS = [0x1, 0x4, 0x10, 0x40, 0x400, 0x800, 0x1000, 0x80000, 0x100000, 0x800000, 0x1000000, 0x2000000,
         0x2000, 0x8000, 0x10000, 0x40000, 0x200000, 0x400000]
MODES = [0x1 | 0x2 | 0x4 | 0x8 | 0x10 | 0x20, 0x1 | 0x2 | 0x8 | 0x10, 0x1 | 0x2 | 0x4 | 0x10 | 0x20, 0x2 | 0x8 | 0x10,
         0x1 | 0x8 | 0x10, 0x1 | 0x2 | 0x8, 0x1 | 0x2 | 0x4 | 0x8 | 0x10 | 0x20 | 0x200, 0x1 | 0x2 | 0x4 | 0x8 | 0x10 | 0x20 | 0x400,
         0x1 | 0x2 | 0x4 | 0x8 | 0x10 | 0x20 | 0x800 | 0x1000]


def fmt(c):
    t = f"{len(c['atoms'])} {c['which']} {c['mode']} {c['flags']} {c['done']} {c['loose']}\n"
    for a in c["atoms"]:
        t += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['chg']} {a['nh']} {a['iso'][0]} {a['iso'][1]} {a['iso'][2]} {a['aw']} " \
             f"{a['xyz'][0]!r} {a['xyz'][1]!r} {a['xyz'][2]!r}\n"
        t += " ".join(f"{x} {t2}" for x, t2 in a["nb"]) + "\n"
    return t


def main():
    rng = random.Random(4242)
    cases = []
    for _ in range(3000):
        atoms = gen(rng)
        flags = 0
        for f in FLAGS:
            if rng.random() < 0.5:
                flags |= f
        flags |= 0x1
        which = rng.choice([3, 3, 3, 7, 1, 2, 6])
        cases.append({"atoms": atoms, "which": which, "mode": rng.choice(MODES), "flags": flags,
                      "done": rng.choice([0, 0, 0, 0x400]), "loose": rng.randrange(2)})
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=(f"{len(cases)}\n" + "".join(fmt(c) for c in cases)).encode(), capture_output=True)
        if r.returncode == 0:
            out = r.stdout.decode().strip().split("\n")
        else:
            kept, lines = [], []
            for c in cases:
                r1 = subprocess.run([str(exe)], input=("1\n" + fmt(c)).encode(), capture_output=True)
                if r1.returncode == 0:
                    kept.append(c)
                    lines += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode().strip().split("\n")[1:]
            print(len(cases) - len(kept), "cases dropped")
            cases = kept
            out = lines
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": fmt(c), "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_create.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
