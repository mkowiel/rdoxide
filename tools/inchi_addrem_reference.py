"""Generates tests/data/inchi_addrem.json from the InChI C code: the
adding and removing of protons that makes the normalization of a
reconstructed structure give the original one (`AddRemoveProtonsRestr`,
`AddRemoveIsoProtonsRestr`, `bHasChargedNeighbor`, `mark_at_type`,
`bHeteroAtomMayHaveXchgIsoH`) on random acids, bases and salts. The function
texts are cut out of ichi_bns.c and strutil.c. Needs a C compiler and the
fetched InChI sources."""

import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src"

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
#define LOG_NO_ARGS(m)
#define LOG_MULT_ARGS(f, ...)
#define RI_ERR_PROGR (-3)
#define FLAG_PROTON_AC_HARD_ADDED 0x0020
typedef enum tagIonAtomType { IAT_H = 0, IAT_C, IAT_N, IAT_P, IAT_O, IAT_S, IAT_Se, IAT_Te, IAT_F, IAT_Cl, IAT_Br, IAT_I, IAT_MAX = 12 } ION_ATOM_TYPE;
@MACROS@
@FUNCS@
static void dumpatoms(inp_ATOM *at, int total, T_GROUP_INFO *t, int *totals) {
    for (int i = 0; i < total; i++) {
        printf("A %d %d %d %d %d,%d,%d t%d s%d,%d,%d o%d,%d,%d n%d,%d,%d g%d,%d,%d r%d nb%d e%d w%d p%d v%d\n", at[i].charge, at[i].num_H, at[i].radical, at[i].chem_bonds_valence, at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2], at[i].at_type,
               at[i].sb_parity[0], at[i].sb_parity[1], at[i].sb_parity[2], at[i].sb_ord[0], at[i].sb_ord[1], at[i].sb_ord[2],
               at[i].sn_ord[0], at[i].sn_ord[1], at[i].sn_ord[2], at[i].sn_orig_at_num[0], at[i].sn_orig_at_num[1], at[i].sn_orig_at_num[2], at[i].orig_at_number, at[i].neighbor[0], at[i].endpoint, at[i].iso_atw_diff, at[i].p_parity, at[i].valence);
        printf("B");
        for (int j = 0; j < at[i].valence; j++) printf(" %d.%d", at[i].neighbor[j], at[i].bond_type[j]);
        printf("\n");
    }
    printf("I %d %d %d,%d,%d\n", t->tni.nNumRemovedExplicitH, t->tni.nNumRemovedProtons, t->tni.nNumRemovedProtonsIsotopic[0], t->tni.nNumRemovedProtonsIsotopic[1], t->tni.nNumRemovedProtonsIsotopic[2]);
    printf("T");
    for (int i = 0; i < ATTOT_ARRAY_LEN; i++) printf(" %d", totals[i]);
    printf("\n");
}
static void queries(inp_ATOM *at, int n) {
    for (int i = 0; i < n; i++) {
        int sub = 0;
        int r = bIsHardRemHCandidate(at, i, &sub);
        printf("Q %d %d %d %d %d\n", r, sub, bHasAcidicHydrogen(at, i), bHasOtherExchangableH(at, i), bHasAcidicMinus(at, i));
    }
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, nH, nops;
        scanf("%d %d %d", &n, &nH, &nops);
        inp_ATOM *at = calloc(n + nH + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, chg, nh, i0, i1, i2, orig, ep, pp, iw, sp[3], so[3], sn[3], sg[3];
            scanf("%d %d %d %d %d %d %d %d %d %d %d %d %d", &v, &cbv, &el, &rad, &chg, &nh, &i0, &i1, &i2, &orig, &ep, &pp, &iw);
            scanf("%d %d %d %d %d %d %d %d %d %d %d %d", sp, sp + 1, sp + 2, so, so + 1, so + 2, sn, sn + 1, sn + 2, sg, sg + 1, sg + 2);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].radical = rad;
            at[i].charge = chg; at[i].num_H = nh; at[i].num_iso_H[0] = i0; at[i].num_iso_H[1] = i1; at[i].num_iso_H[2] = i2; at[i].orig_at_number = orig; at[i].endpoint = ep; at[i].p_parity = pp; at[i].iso_atw_diff = iw;
            for (int k = 0; k < 3; k++) { at[i].sb_parity[k] = sp[k]; at[i].sb_ord[k] = so[k]; at[i].sn_ord[k] = sn[k]; at[i].sn_orig_at_num[k] = sg[k]; }
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        for (int i = 0; i < nH; i++) {
            int nb0, orig, iw;
            scanf("%d %d %d", &nb0, &orig, &iw);
            at[n + i].neighbor[0] = nb0; at[n + i].orig_at_number = orig; at[n + i].el_number = 1; at[n + i].valence = 1; at[n + i].iso_atw_diff = iw;
        }
        T_GROUP_INFO tgi;
        memset(&tgi, 0, sizeof(tgi));
        tgi.tni.nNumRemovedExplicitH = nH;
        tgi.tni.nNumRemovedProtons = 1;
        int totals[ATTOT_ARRAY_LEN];
        memset(totals, 0, sizeof(totals));
        BN_AATG aatg;
        memset(&aatg, 0, sizeof(aatg));
        aatg.nAtTypeTotals = totals;
        aatg.t_group_info = &tgi;
        for (int i = 0; i < n; i++) { int mask; GetAtomChargeType(at, i, totals, &mask, 0); }
        printf("N %d\n", c);
        queries(at, n);
        dumpatoms(at, n + nH, &tgi, totals);
        for (int k = 0; k < nops; k++) {
            int op, a1, a2, a3, a4, a5, a6;
            scanf("%d %d %d %d %d %d %d", &op, &a1, &a2, &a3, &a4, &a5, &a6);
            if (op == 0) {
                int np = a1;
                int r = AddRemoveProtonsRestr(at, n, &np, a2, a3 ? FLAG_PROTON_AC_HARD_ADDED : 0, a4, 0, 0);
                printf("R 0 %d %d\n", r, np);
            } else {
                NUM_H np[3] = {a1, a2, a3};
                int r = AddRemoveIsoProtonsRestr(at, n, np, a4);
                printf("R 1 %d %d %d %d\n", r, np[0], np[1], np[2]);
            }
            queries(at, n);
            dumpatoms(at, n + nH, &tgi, totals);
        }
        free(at);
    }
    return 0;
}
"""

# functional groups: (element, charge, H count) atoms, the first one binds to
# the skeleton; bonds as (i, j, order)
GROUPS = [
    ([(6, 0, 0), (8, 0, 0), (8, 0, 1)], [(0, 1, 2), (0, 2, 1)]),           # COOH
    ([(6, 0, 0), (8, 0, 0), (8, -1, 0)], [(0, 1, 2), (0, 2, 1)]),          # COO-
    ([(16, 0, 0), (8, 0, 0), (8, 0, 0), (8, 0, 1)], [(0, 1, 2), (0, 2, 2), (0, 3, 1)]),
    ([(16, 0, 0), (8, 0, 0), (8, 0, 0), (8, -1, 0)], [(0, 1, 2), (0, 2, 2), (0, 3, 1)]),
    ([(15, 0, 0), (8, 0, 0), (8, 0, 1), (8, 0, 1)], [(0, 1, 2), (0, 2, 1), (0, 3, 1)]),
    ([(15, 0, 0), (8, 0, 0), (8, -1, 0), (8, 0, 1)], [(0, 1, 2), (0, 2, 1), (0, 3, 1)]),
    ([(7, 1, 3)], []),                                                      # NH3+
    ([(7, 0, 2)], []),                                                      # NH2
    ([(7, 0, 1), (6, 0, 0), (8, 0, 0)], [(0, 1, 1), (1, 2, 2)]),            # amide
    ([(7, 1, 0), (8, 0, 0), (8, -1, 0)], [(0, 1, 2), (0, 2, 1)]),           # nitro
    ([(8, 0, 1)], []),                                                      # OH
    ([(8, -1, 0)], []),                                                     # O-
    ([(8, 1, 1)], []),                                                      # OH+
    ([(7, 0, 0), (6, 0, 1)], [(0, 1, 2)]),                                  # imine
    ([(7, -1, 0)], []),                                                     # N-
    ([(7, 1, 1), (6, 0, 1)], [(0, 1, 2)]),                                  # iminium
    ([(7, 0, 0), (8, 0, 1)], [(0, 1, 1)]),                                  # N-OH (hydroxylamine-like)
    ([(7, 0, 0), (8, -1, 0)], [(0, 1, 1)]),
    ([(7, 0, 0), (8, 0, 0)], [(0, 1, 2)]),                                  # nitroso
    ([(16, 0, 1)], []),                                                     # SH
    ([(16, -1, 0)], []),                                                    # S-
    ([(6, 0, 0), (16, 0, 0), (8, 0, 1)], [(0, 1, 2), (0, 2, 1)]),           # C(=S)OH
    ([(6, 0, 0), (16, -1, 0), (8, 0, 0)], [(0, 1, 1), (0, 2, 2)]),
]
LONE = [(17, -1, 0), (17, 0, 1), (35, -1, 0), (7, 1, 4), (8, 0, 2), (8, 1, 3), (1, 1, 0)]


def gen(rng):
    atoms = []  # [el, chg, nh, neighbors]

    def add(el, chg, nh):
        atoms.append([el, chg, nh, []])
        return len(atoms) - 1

    def bond(a, b, t):
        atoms[a][3].append((b, t))
        atoms[b][3].append((a, t))

    skel = [add(6, 0, 0)]
    for _ in range(rng.randrange(0, 3)):
        c = add(6, 0, 0)
        bond(rng.choice(skel), c, 1)
        skel.append(c)
    for _ in range(rng.randrange(1, 4)):
        els, bonds = rng.choice(GROUPS)
        ids = [add(*e) for e in els]
        for i, j, t in bonds:
            bond(ids[i], ids[j], t)
        bond(rng.choice(skel), ids[0], 1)
    if rng.random() < 0.3:
        add(*rng.choice(LONE))
    out = []
    for el, chg, nh, nb in atoms:
        cbv = sum(t for _, t in nb)
        if el == 6 and nb and chg == 0:
            nh = max(0, 4 - cbv) if nh == 0 else nh
        out.append({"v": len(nb), "cbv": cbv, "el": el, "rad": 0, "chg": chg, "nh": nh, "nb": nb})
    return out


def main():
    rng = random.Random(21)
    cases = []
    for _ in range(900):
        atoms = gen(rng)
        n = len(atoms)
        nh = rng.randrange(0, 4)
        origs = [i + 1 for i in range(n)] + [100 + j for j in range(nh)] + [0]
        for i, a in enumerate(atoms):
            a["iso"] = [rng.choice([0, 0, 1]) for _ in range(3)]
            a["orig"] = i + 1
            k = rng.randrange(0, 4)
            a["sbp"] = [rng.choice([1, 2, 3, 4]) if j < k else 0 for j in range(3)]
            v = max(a["v"], 1)
            a["sbo"] = [rng.randrange(v) for _ in range(3)]
            a["sno"] = [rng.randrange(-1, 2) for _ in range(3)]
            a["sng"] = [rng.choice(origs) for _ in range(3)]
        hs = [(rng.randrange(n), 100 + j, rng.choice([0, 0, 1, 2])) for j in range(nh)]
        for i, a in enumerate(atoms):
            a["ep"] = rng.choice([0, 0, 0, 1, 1, 2])
            a["pp"] = rng.choice([0] * 12 + [1])
            a["iw"] = rng.choice([0] * 6 + [1])
        # explicit H atoms (n..) are neighbors at the start of the list of the parent
        for h, (parent, _, _) in enumerate(hs):
            a = atoms[parent]
            a["nb"] = [(n + h, 1)] + list(a["nb"])
            a["v"] = len(a["nb"])
        ops = []
        for _ in range(rng.randrange(1, 5)):
            if rng.random() < 0.6:
                ops.append((0, rng.randrange(-3, 4), rng.randrange(-2, 3), rng.choice([0, 1]), rng.choice([0, 1, 1, 2]), 0, 0))
            else:
                ops.append((1, rng.choice([0, 1, 2, 3, -1]), rng.choice([0, 1, 2]), rng.choice([0, 1, 2]), rng.choice([0, 1]), 0, 0))
        cases.append({"atoms": atoms, "hs": hs, "ops": ops})

    def chunk(cs):
        t = f"{len(cs)}\n"
        for c in cs:
            t += f"{len(c['atoms'])} {len(c['hs'])} {len(c['ops'])}\n"
            for a in c["atoms"]:
                t += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['chg']} {a['nh']} {a['iso'][0]} {a['iso'][1]} {a['iso'][2]} {a['orig']} {a['ep']} {a['pp']} {a['iw']}\n"
                t += " ".join(map(str, a["sbp"] + a["sbo"] + a["sno"] + a["sng"])) + "\n"
                t += " ".join(f"{x} {t2}" for x, t2 in a["nb"]) + "\n"
            for h in c["hs"]:
                t += f"{h[0]} {h[1]} {h[2]}\n"
            for o in c["ops"]:
                t += " ".join(map(str, o)) + "\n"
        return t

    text = (SRC / "ichi_bns.c").read_text(errors="replace")

    def cut(start, src=text):
        a = src.rindex(start)
        return src[a:src.index("\n}\n", a) + 3]

    names = [
        "int GetAtomChargeType( inp_ATOM *atom,\n", "int IsZOX( inp_ATOM *atom, int at_x, int ord )\n{", "int is_Z_atom( U_CHAR el_number )\n{",
        "int AddOrRemoveExplOrImplH( int nDelta,\n",
        "int SimpleRemoveHplusNPO( inp_ATOM *at,\n", "int mark_at_type( inp_ATOM *atom, int num_atoms, int nAtTypeTotals[] )\n{", "int bHasChargedNeighbor( inp_ATOM *at, int iat )\n{", "int AddRemoveProtonsRestr( inp_ATOM *at,\n", "int AddRemoveIsoProtonsRestr( inp_ATOM *at,\n", "int bIsAtomTypeHard( inp_ATOM *at,\n", "int bIsHDonorAccAtomType( inp_ATOM *at, int endpoint, int *cSubType )\n{",
        "int bIsNegAtomType( inp_ATOM *at, int endpoint, int *cSubType )\n{", "int bIsHardRemHCandidate( inp_ATOM *at, int i, int *cSubType )\n{",
        "int SimpleRemoveAcidicProtons( inp_ATOM *at,\n", "int bHasAcidicHydrogen( inp_ATOM *at, int i )\n{", "int bHasOtherExchangableH( inp_ATOM *at, int i )\n{",
        "int SimpleAddAcidicProtons( inp_ATOM *at,\n", "int bHasAcidicMinus( inp_ATOM *at, int i )\n{"]
    funcs = []
    for arr in ("int ArTypMask[] =", "int AaTypMask[] ="):
        a = text.rindex(arr)
        funcs.append(text[a:text.index("};", a) + 3])
    funcs += [cut(n) for n in names]
    sut = (SRC / "strutil.c").read_text(errors="replace")
    funcs.append(cut("int get_iat_number( int el_number )\n{", sut))
    funcs.append(cut("int bHeteroAtomMayHaveXchgIsoH( inp_ATOM *atom, int iat )\n{", sut))
    st = (SRC / "ichister.c").read_text(errors="replace")
    funcs.append(cut("int get_opposite_sb_atom(", st))
    taut = (SRC / "ichitaut.c").read_text(errors="replace")
    funcs.append(cut("int is_centerpoint_elem( U_CHAR el_number )\n{", taut))
    ea = text.index("typedef enum tagAtTypeTotals")
    eb = text.index("#define AA_HARD_TYP_CO")
    macros = text[ea:text.index("\n", eb)]
    protos = """int mark_at_type( inp_ATOM *atom, int num_atoms, int nAtTypeTotals[] );
int bHasChargedNeighbor( inp_ATOM *at, int iat );
int AddRemoveProtonsRestr( inp_ATOM *at, int num_atoms, int *num_protons_to_add, int nNumProtAddedByRestr, INCHI_MODE bNormalizationFlags, int num_tg, int nChargeRevrs, int nChargeInChI );
int AddRemoveIsoProtonsRestr( inp_ATOM *at, int num_atoms, NUM_H num_protons_to_add[], int num_tg );
int get_iat_number( int el_number );
int bHeteroAtomMayHaveXchgIsoH( inp_ATOM *atom, int iat );
int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt, int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );
int AddOrRemoveExplOrImplH( int nDelta, inp_ATOM *at, int num_atoms, AT_NUMB at_no, T_GROUP_INFO *t_group_info );
int GetAtomChargeType( inp_ATOM *atom, int at_no, int nAtTypeTotals[], int *pMask, int bSubtract );
int IsZOX( inp_ATOM *atom, int at_x, int ord );
int is_Z_atom( U_CHAR el_number );
int SimpleRemoveHplusNPO( inp_ATOM *at, int num_atoms, int nAtTypeTotals[], T_GROUP_INFO *t_group_info );
int bIsHardRemHCandidate( inp_ATOM *at, int i, int *cSubType );
int SimpleRemoveAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2remove );
int SimpleAddAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2add );
int bHasAcidicHydrogen( inp_ATOM *at, int i );
int bHasOtherExchangableH( inp_ATOM *at, int i );
int bHasAcidicMinus( inp_ATOM *at, int i );
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        code = HARNESS.replace("@MACROS@", macros).replace("@FUNCS@", protos + "\n".join(funcs))
        (Path(d) / "h.c").write_text(code)
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c"), "-lm"], check=True)
        r = subprocess.run([str(exe)], input=chunk(cases).encode(), capture_output=True)
        if r.returncode == 0:
            out = r.stdout.decode().strip().split("\n")
        else:
            kept, lines = [], []
            for c in cases:
                r1 = subprocess.run([str(exe)], input=chunk([c]).encode(), capture_output=True)
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
    for op in (0, 1):
        print("op", op, "changed:", sum(int(l.split()[2]) for d in data for l in d["out"] if l.startswith(f"R {op} ")))
    path = ROOT / "tests/data/inchi_addrem.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
