"""Generates tests/data/inchi_cano2.json from the InChI C code: `AllocateCS`
and `Canon_INChI` / `Canon_INChI3` of ichicano.c (the final canonical
numbering with the isotopic and stereo optimization) run on the result of
`GetBaseCanonRanking`, followed by `FillOutINChI` of ichimak2.c (the "FO"
lines), compiled with the upstream source included.

    python tools/inchi_cano2_reference.py
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


ma = load("inchi_markalt_reference")
f2 = load("inchi_can2f_reference")
SRC = ma.SRC

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ichicano.c"
int FillOutINChI( INChI *pINChI, INChI_Aux *pINChI_Aux, int num_atoms, int num_at_tg, int num_removed_H, sp_ATOM *at, inp_ATOM *norm_at,
                  CANON_STAT *pCS, CANON_GLOBALS *pCG, int bTautomeric, INCHI_MODE nUserMode, char *pStrErrStruct, int bNoWarnings );
int CreateCompositeNormAtom() { return 0; }
static void pu(const char *t, AT_RANK *a, int n) { printf("%s", t); if (!a) { printf(" -\n"); return; } for (int i = 0; i < n; i++) printf(" %d", a[i]); printf("\n"); }
static void read_atoms(sp_ATOM *at, int n, int *edges, int m) {
    for (int i = 0; i < n; i++) {
        int el, nh, aw, i1, i2, i3, cf, ep;
        scanf("%d %d %d %d %d %d %d %d", &el, &nh, &aw, &i1, &i2, &i3, &cf, &ep);
        const char *names[] = {"C", "N", "O", "S", "Cl", "H"};
        strcpy(at[i].elname, names[el]);
        at[i].num_H = nh; at[i].iso_atw_diff = aw;
        at[i].num_iso_H[0] = i1; at[i].num_iso_H[1] = i2; at[i].num_iso_H[2] = i3;
        at[i].cFlags = cf; at[i].endpoint = ep;
        at[i].iso_sort_key = make_iso_sort_key(aw, i1, i2, i3);
    }
    for (int i = 0; i < m; i++) {
        int a = edges[2 * i], b = edges[2 * i + 1];
        at[a].neighbor[at[a].valence++] = b;
        at[b].neighbor[at[b].valence++] = a;
    }
}
static void read_stereo(sp_ATOM *at, int n, int st[][27]) {
    for (int i = 0; i < n; i++) {
        at[i].parity = st[i][0]; at[i].parity2 = st[i][1]; at[i].stereo_atom_parity = st[i][2]; at[i].stereo_atom_parity2 = st[i][3];
        at[i].nRingSystem = st[i][4];
        for (int j = 0; j < 3; j++) {
            at[i].stereo_bond_neighbor[j] = st[i][5 + 4 * j]; at[i].stereo_bond_ord[j] = st[i][6 + 4 * j];
            at[i].stereo_bond_z_prod[j] = st[i][7 + 4 * j]; at[i].stereo_bond_parity[j] = st[i][8 + 4 * j];
        }
    }
}
static void read_stereo2(sp_ATOM *at, int n, int st[][27]) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < 3; j++) {
            /* the isotopic sb arrays follow the plain ones in the row; with 3*4+5=17 used, 10 left: only
               the first two are stored (see the generator): neighbor, ord, z, parity of j=0 and j=1, 2 spare */
        }
}
static const char *slot(void *p, void *b, void *inv, void *iso, void *isoinv) {
    if (p == b) return "B"; if (p == inv) return "I"; if (p == iso) return "O"; if (p == isoinv) return "P"; return "-";
}
int main(void) {
    int ncases; setvbuf(stdout, NULL, _IONBF, 0);
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        INCHI_CLOCK ic; memset(&ic, 0, sizeof(ic));
        int n, m, ntg, nE, haveYes, haveNon, haveTgi, bfix, large;
        scanf("%d %d %d %d %d %d %d %d %d", &n, &m, &ntg, &nE, &haveYes, &haveNon, &haveTgi, &bfix, &large);
        int edges[200];
        for (int i = 0; i < 2 * m; i++) scanf("%d", &edges[i]);
        ATOM_SIZES s[2]; memset(s, 0, sizeof(s));
        int sv[2][7];
        for (int k = 0; k < 2; k++) for (int j = 0; j < 7; j++) scanf("%d", &sv[k][j]);
        sp_ATOM *atY = calloc(n + 1, sizeof(sp_ATOM)), *atN = calloc(n + 1, sizeof(sp_ATOM));
        read_atoms(atY, n, edges, m);
        read_atoms(atN, n, edges, m);
        T_GROUP *tg = calloc(ntg + 2, sizeof(T_GROUP));
        AT_NUMB *epl = calloc(nE + 2, sizeof(AT_NUMB));
        for (int k = 0; k < ntg; k++) {
            int ne, fp; long w;
            scanf("%d %d %ld", &ne, &fp, &w);
            tg[k].nNumEndpoints = ne; tg[k].nFirstEndpointAtNoPos = fp; tg[k].nGroupNumber = k + 1;
            for (int j = 0; j < T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC; j++) { int x; scanf("%d", &x); tg[k].num[j] = x; }
            long iw; scanf("%ld", &iw); tg[k].iWeight = iw;
        }
        for (int i = 0; i < nE; i++) { int x; scanf("%d", &x); epl[i] = x; }
        int nIso, nNumIsoEnd, flagsDone, normFlags, ignIso, tgNull;
        scanf("%d %d %d %d %d %d", &nIso, &nNumIsoEnd, &flagsDone, &normFlags, &ignIso, &tgNull);
        AT_NUMB *isoList = calloc(n + 3, sizeof(AT_NUMB));
        isoList[0] = nIso;
        for (int i = 0; i < nIso; i++) { int x; scanf("%d", &x); isoList[i + 1] = x; }
        /* stereo: plain row [17] + isotopic row [17] per atom */
        static int st1[64][27], st2[64][27];
        for (int i = 0; i < n; i++) for (int j = 0; j < 17; j++) scanf("%d", &st1[i][j]);
        for (int i = 0; i < n; i++) for (int j = 0; j < 17; j++) scanf("%d", &st2[i][j]);
        int modeBits, uu;
        scanf("%d %d", &modeBits, &uu);
        /* the plain stereo goes to parity.., the isotopic one to parity2.. */
        for (int a = 0; a < 2; a++) {
            sp_ATOM *at = a ? atN : atY;
            for (int i = 0; i < n; i++) {
                at[i].parity = st1[i][0]; at[i].stereo_atom_parity = st1[i][2]; at[i].nRingSystem = st1[i][4];
                at[i].parity2 = st2[i][0]; at[i].stereo_atom_parity2 = st2[i][2];
                for (int j = 0; j < 3; j++) {
                    at[i].stereo_bond_neighbor[j] = st1[i][5 + 4 * j]; at[i].stereo_bond_ord[j] = st1[i][6 + 4 * j];
                    at[i].stereo_bond_z_prod[j] = st1[i][7 + 4 * j]; at[i].stereo_bond_parity[j] = st1[i][8 + 4 * j];
                    at[i].stereo_bond_neighbor2[j] = st2[i][5 + 4 * j]; at[i].stereo_bond_ord2[j] = st2[i][6 + 4 * j];
                    at[i].stereo_bond_z_prod2[j] = st2[i][7 + 4 * j]; at[i].stereo_bond_parity2[j] = st2[i][8 + 4 * j];
                }
            }
        }
        T_GROUP_INFO tgi; memset(&tgi, 0, sizeof(tgi));
        tgi.t_group = tgNull ? NULL : tg; tgi.nEndpointAtomNumber = epl; tgi.num_t_groups = ntg; tgi.max_num_t_groups = ntg;
        tgi.nNumEndpoints = nE;
        tgi.nIsotopicEndpointAtomNumber = isoList; tgi.nNumIsotopicEndpoints = nNumIsoEnd;
        tgi.bTautFlagsDone = flagsDone; tgi.tni.bNormalizationFlags = normFlags; tgi.bIgnoreIsotopic = ignIso;
        tgi.tGroupNumber = calloc(4 * ntg + 4, sizeof(AT_NUMB));
        /* sizes from GetCanonLengths */
        sp_ATOM *ats[2]; ats[TAUT_NON] = haveNon ? atN : NULL; ats[TAUT_YES] = haveYes ? atY : NULL;
        for (int k = 0; k < 2; k++) {
            s[k].nLenLinearCTTautomer = sv[k][4];
            if (ats[k]) GetCanonLengths(n, ats[k], &s[k], haveTgi ? &tgi : NULL);
            s[k].num_isotopic_atoms = sv[k][2]; s[k].bIgnoreIsotopic = sv[k][3];
            s[k].bHasIsotopicTautGroups = sv[k][5]; s[k].nLenIsotopicEndpoints = sv[k][6] > s[k].nLenIsotopicEndpoints ? sv[k][6] : s[k].nLenIsotopicEndpoints;
        }
        BCN bcn; memset(&bcn, 0, sizeof(bcn));
        int ntot = n + ntg;
        printf("N %d\n", c);
        int ret = GetBaseCanonRanking(&ic, n, ntot, ats, haveTgi ? &tgi : NULL, s, &bcn, NULL, &cg, bfix, large);
        printf("R %d\n", ret);
        if (ret != 0) continue;
        for (int k = 0; k < 2; k++) {
            FTCN *f = &bcn.ftcn[k];
            if (!f->LinearCt) continue;
            sp_ATOM *at = k == TAUT_YES ? atY : atN;
            if (!(k == TAUT_YES ? haveYes : haveNon) && f->LinearCt) at = k == TAUT_YES ? atN : atY;
            CANON_STAT cs; memset(&cs, 0, sizeof(cs));
            int nat = f->num_at_tg, nA = f->num_atoms;
            int mode = CMODE_CT | modeBits;
            if (nat > nA) mode |= CMODE_TAUT;
            int nds = s[k].nLenLinearCTStereoDble, ncs = s[k].nLenLinearCTStereoCarb;
            int ids, ics;
            {
                ATOM_SIZES si; memset(&si, 0, sizeof(si));
                int sw = 0;
                SwitchAtomStereoAndIsotopicStereo(at, nA, &sw);
                GetCanonLengths(nA, at, &si, NULL);
                SwitchAtomStereoAndIsotopicStereo(at, nA, &sw);
                ids = si.nLenLinearCTStereoDble; ics = si.nLenLinearCTStereoCarb;
            }
            int lenIso = 0;
            for (int q = 0; q < nA; q++) { int kk = !at[q].endpoint && !(at[q].cFlags & AT_FLAG_ISO_H_POINT) && (at[q].num_iso_H[0] || at[q].num_iso_H[1] || at[q].num_iso_H[2]); lenIso += (at[q].iso_atw_diff || kk); }
            int lenIT = 0; for (int q = 0; q < ntg; q++) lenIT += tg[q].iWeight != 0;
            if (haveTgi && tgi.nNumIsotopicEndpoints) lenIT = 0;
            AllocateCS(&cs, nA, nat, f->nLenLinearCt, f->nLenLinearCtAtOnly, nds, ids, ncs, ics,
                       3 * ntg + nE + 1, lenIT, lenIso, mode, &bcn);
            cs.pBCN = &bcn; cs.t_group_info = haveTgi ? &tgi : NULL;
            cs.bIgnoreIsotopic = 0;
            void *pb_d = cs.LinearCTStereoDble, *pi_d = cs.LinearCTStereoDbleInv, *po_d = cs.LinearCTIsotopicStereoDble, *pp_d = cs.LinearCTIsotopicStereoDbleInv;
            void *pb_c = cs.LinearCTStereoCarb, *pi_c = cs.LinearCTStereoCarbInv, *po_c = cs.LinearCTIsotopicStereoCarb, *pp_c = cs.LinearCTIsotopicStereoCarbInv;
            int r2 = Canon_INChI(&ic, nA, nat, at, &cs, &cg, mode | (uu ? REQ_MODE_DIFF_UU_STEREO : 0), k);
            printf("C%d %d %d %d %d %d %d\n", k, r2, mode, cs.nLenLinearCT, cs.nLenLinearCTAtOnly, cs.nCanonFlags, cs.bIgnoreIsotopic);
            printf("SL %s%s%s%s\n", slot(cs.LinearCTStereoDble, pb_d, pi_d, po_d, pp_d), slot(cs.LinearCTStereoCarb, pb_c, pi_c, po_c, pp_c),
                   slot(cs.LinearCTStereoDbleInv, pb_d, pi_d, po_d, pp_d), slot(cs.LinearCTStereoCarbInv, pb_c, pi_c, po_c, pp_c));
            pu("LC", cs.LinearCT, cs.LinearCT ? cs.nMaxLenLinearCT : 0);
            pu("L2", cs.LinearCT2, cs.LinearCT2 ? cs.nMaxLenLinearCT : 0);
            printf("LE %d %d %d %d %d %d %d %d\n", cs.nLenLinearCT2, cs.nLenLinearCTAtOnly2, cs.nLenLinearCTTautomer, cs.nLenLinearCTIsotopic, cs.nLenLinearCTIsotopicTautomer, cs.bCmpStereo, cs.bCmpIsotopicStereo, cs.nMode);
            printf("NH"); if (cs.nNum_H) for (int i = 0; i < nA; i++) printf(" %d", cs.nNum_H[i]); printf("\n");
            printf("NF"); if (cs.nNum_H_fixed) for (int i = 0; i < nA; i++) printf(" %d", cs.nNum_H_fixed[i]); printf("\n");
            printf("NX"); if (cs.nExchgIsoH) for (int i = 0; i < nA; i++) printf(" %d", cs.nExchgIsoH[i]); printf("\n");
            pu("LT", cs.LinearCTTautomer, cs.nMaxLenLinearCTTautomer);
            printf("IT"); for (int i = 0; i < cs.nMaxLenLinearCTIsotopicTautomer; i++) printf(" %d:%d,%d,%d", cs.LinearCTIsotopicTautomer[i].tgroup_num, cs.LinearCTIsotopicTautomer[i].num[0], cs.LinearCTIsotopicTautomer[i].num[1], cs.LinearCTIsotopicTautomer[i].num[2]); printf("\n");
            printf("IA"); for (int i = 0; i < cs.nMaxLenLinearCTIsotopic; i++) printf(" %d,%d,%d,%d,%d", cs.LinearCTIsotopic[i].at_num, cs.LinearCTIsotopic[i].num_1H, cs.LinearCTIsotopic[i].num_D, cs.LinearCTIsotopic[i].num_T, cs.LinearCTIsotopic[i].iso_atw_diff); printf("\n");
            pu("SR", cs.nSymmRank, cs.nSymmRank ? nat : 0);
            pu("CO", cs.nCanonOrd, cs.nCanonOrd ? nat : 0);
            pu("SI", cs.nSymmRankIsotopic, cs.nSymmRankIsotopic ? nat : 0);
            pu("CI", cs.nCanonOrdIsotopic, cs.nCanonOrdIsotopic ? nat : 0);
            pu("ST", cs.nCanonOrdStereo, cs.nCanonOrdStereo ? nat : 0);
            pu("SV", cs.nCanonOrdStereoInv, cs.nCanonOrdStereoInv ? nat : 0);
            pu("IS", cs.nCanonOrdIsotopicStereo, cs.nCanonOrdIsotopicStereo ? nat : 0);
            pu("IV", cs.nCanonOrdIsotopicStereoInv, cs.nCanonOrdIsotopicStereoInv ? nat : 0);
            int ntgs = haveTgi ? tgi.num_t_groups : 0;
            pu("TO", cs.nCanonOrdTaut, cs.nCanonOrdTaut ? ntgs : 0);
            pu("TS", cs.nSymmRankTaut, cs.nSymmRankTaut ? ntgs : 0);
            pu("TI", cs.nCanonOrdIsotopicTaut, cs.nCanonOrdIsotopicTaut ? ntgs : 0);
            pu("TJ", cs.nSymmRankIsotopicTaut, cs.nSymmRankIsotopicTaut ? ntgs : 0);
            pu("TT", cs.nCanonOrdStereoTaut, cs.nCanonOrdStereoTaut ? ntgs : 0);
            pu("TU", cs.nCanonOrdIsotopicStereoTaut, cs.nCanonOrdIsotopicStereoTaut ? ntgs : 0);
            printf("CL %d %d %d %d %d %d %d %d\n", cs.nLenCanonOrd, cs.nLenCanonOrdTaut, cs.nLenCanonOrdIsotopic, cs.nLenCanonOrdIsotopicTaut, cs.nLenCanonOrdStereo, cs.nLenCanonOrdStereoTaut, cs.nLenCanonOrdIsotopicStereo, cs.nLenCanonOrdIsotopicStereoTaut);
            printf("SS %d %d %d %d %d %d %d %d\n", cs.nLenLinearCTStereoDble, cs.nLenLinearCTStereoCarb, cs.nLenLinearCTStereoDbleInv, cs.nLenLinearCTStereoCarbInv, cs.nLenLinearCTIsotopicStereoDble, cs.nLenLinearCTIsotopicStereoCarb, cs.nLenLinearCTIsotopicStereoDbleInv, cs.nLenLinearCTIsotopicStereoCarbInv);
            printf("SM %d %d %d %d\n", cs.nMaxLenLinearCTStereoDble, cs.nMaxLenLinearCTStereoCarb, cs.nMaxLenLinearCTIsotopicStereoDble, cs.nMaxLenLinearCTIsotopicStereoCarb);
            AT_STEREO_DBLE *dd[4] = {pb_d, pi_d, po_d, pp_d}; AT_STEREO_CARB *cc[4] = {pb_c, pi_c, po_c, pp_c};
            int md[4] = {cs.nMaxLenLinearCTStereoDble, cs.nMaxLenLinearCTStereoDble, cs.nMaxLenLinearCTIsotopicStereoDble, cs.nMaxLenLinearCTIsotopicStereoDble};
            int mc[4] = {cs.nMaxLenLinearCTStereoCarb, cs.nMaxLenLinearCTStereoCarb, cs.nMaxLenLinearCTIsotopicStereoCarb, cs.nMaxLenLinearCTIsotopicStereoCarb};
            int ad[4] = {nds, nds, ids, ids}, ac[4] = {ncs, ncs, ics, ics};
            for (int q = 0; q < 4; q++) {
                printf("D%d", q); if (dd[q]) for (int i = 0; i < ad[q]; i++) printf(" %d,%d,%d", dd[q][i].at_num1, dd[q][i].at_num2, dd[q][i].parity); printf("\n");
                printf("K%d", q); if (cc[q]) for (int i = 0; i < ac[q]; i++) printf(" %d,%d", cc[q][i].at_num, cc[q][i].parity); printf("\n");
            }
            for (int i = 0; i < n; i++) {
                printf("A %d %d %d %d %d %d %d %d |", at[i].parity, at[i].parity2, at[i].final_parity, at[i].final_parity2, at[i].stereo_atom_parity, at[i].stereo_atom_parity2, at[i].bHasStereoOrEquToStereo, at[i].bHasStereoOrEquToStereo2);
                for (int j = 0; j < 3; j++) printf(" %d,%d,%d,%d", at[i].stereo_bond_neighbor[j], at[i].stereo_bond_parity[j], at[i].stereo_bond_ord[j], at[i].stereo_bond_z_prod[j]);
                printf(" |");
                for (int j = 0; j < 3; j++) printf(" %d,%d,%d,%d", at[i].stereo_bond_neighbor2[j], at[i].stereo_bond_parity2[j], at[i].stereo_bond_ord2[j], at[i].stereo_bond_z_prod2[j]);
                printf("\n");
            }
            printf("TG"); if (haveTgi) { for (int i = 0; i < 4 * ntgs; i++) printf(" %d", tgi.tGroupNumber[i]); printf(" |"); for (int i = 0; i < nE; i++) printf(" %d", epl[i]); } printf("\n");
            printf("TI %d\n", haveTgi ? tgi.bIgnoreIsotopic : -1);
            printf("PN %d\n", f->NeighList != NULL);
            printf("RS %d\n", bcn.nMaxLenRankStack);
            for (int j = 0; j < bcn.nMaxLenRankStack / 2; j++) {
                AT_RANK *a = bcn.pRankStack ? bcn.pRankStack[2 * j] : NULL, *b = bcn.pRankStack ? bcn.pRankStack[2 * j + 1] : NULL;
                char t[16]; sprintf(t, "Q%d", j); pu(t, a, a ? bcn.num_max : 0);
                sprintf(t, "W%d", j); pu(t, b, b ? bcn.num_max : 0);
            }

            if (r2 > 0) {
                int umodes[6] = {0, 0x200, 0x400, 0x800, 0x1000, 0x1800};
                INCHI_MODE um = umodes[c % 6];
                INChI pi; memset(&pi, 0, sizeof(pi));
                INChI_Aux pa; memset(&pa, 0, sizeof(pa));
                int cap = 4 * (nat + nE) + 40;
                pi.nAtom = calloc(nat + 2, 1); pi.nConnTable = calloc(f->nLenLinearCt + nat + 10, sizeof(AT_NUMB));
                pi.nTautomer = calloc(cap, sizeof(AT_NUMB)); pi.nNum_H = calloc(nat + 2, 1); pi.nNum_H_fixed = calloc(nat + 2, 1);
                pi.IsotopicAtom = calloc(nA + 2, sizeof(INChI_IsotopicAtom)); pi.IsotopicTGroup = calloc(ntg + 2, sizeof(INChI_IsotopicTGroup));
                pi.nPossibleLocationsOfIsotopicH = calloc(nA + 3, sizeof(AT_NUMB));
                INChI_Stereo *sts[2];
                for (int q = 0; q < 2; q++) {
                    INChI_Stereo *s = calloc(1, sizeof(*s));
                    s->nNumber = calloc(64, sizeof(AT_NUMB)); s->t_parity = calloc(64, 1);
                    s->nNumberInv = calloc(64, sizeof(AT_NUMB)); s->t_parityInv = calloc(64, 1);
                    s->nBondAtom1 = calloc(64, sizeof(AT_NUMB)); s->nBondAtom2 = calloc(64, sizeof(AT_NUMB)); s->b_parity = calloc(64, 1);
                    sts[q] = s;
                }
                pi.Stereo = sts[0]; pi.StereoIsotopic = sts[1];
                int al = 2 * nA + 8;
                pa.nOrigAtNosInCanonOrd = calloc(al, sizeof(AT_NUMB)); pa.nIsotopicOrigAtNosInCanonOrd = calloc(al, sizeof(AT_NUMB));
                pa.nOrigAtNosInCanonOrdInv = calloc(al, sizeof(AT_NUMB)); pa.nIsotopicOrigAtNosInCanonOrdInv = calloc(al, sizeof(AT_NUMB));
                pa.nConstitEquNumbers = calloc(al, sizeof(AT_NUMB)); pa.nConstitEquIsotopicNumbers = calloc(al, sizeof(AT_NUMB));
                pa.nConstitEquTGroupNumbers = calloc(al, sizeof(AT_NUMB)); pa.nConstitEquIsotopicTGroupNumbers = calloc(al, sizeof(AT_NUMB));
                pa.OrigInfo = calloc(al, sizeof(ORIG_INFO));
                if (haveTgi) { pa.bNormalizationFlags = tgi.tni.bNormalizationFlags; pa.bTautFlags = tgi.bTautFlags; pa.bTautFlagsDone = tgi.bTautFlagsDone; }
                inp_ATOM *na = calloc(nA + 2, sizeof(inp_ATOM));
                for (int i = 0; i < nA; i++) {
                    at[i].el_number = get_periodic_table_number(at[i].elname);
                    at[i].charge = (i * 7 + 3) % 4 - 1;
                    at[i].radical = i % 5 == 4 ? 1 : (i % 7 == 3 ? 2 : 0);
                    at[i].orig_at_number = 101 + i;
                    at[i].chem_bonds_valence = at[i].valence;
                    na[i].valence = at[i].valence; na[i].num_H = at[i].num_H; na[i].charge = at[i].charge; na[i].radical = at[i].radical;
                    na[i].el_number = at[i].el_number; na[i].chem_bonds_valence = at[i].chem_bonds_valence;
                }
                char se[STR_ERR_LEN]; se[0] = 0;
                int fr = FillOutINChI(&pi, &pa, nA, nat, 0, at, na, &cs, &cg, k, um, se, 0);
                printf("FO %d %d %d %lu %d %d %s\n", fr, pi.nErrorCode, pa.nErrorCode, (unsigned long)pi.nFlags, pi.nTotalCharge, pi.nNumberOfAtoms, pi.szHillFormula ? pi.szHillFormula : "NULL");
                printf("FA"); for (int i = 0; i < nA; i++) printf(" %d", pi.nAtom[i]); printf("\n");
                printf("FC %d", pi.lenConnTable); for (int i = 0; i < pi.lenConnTable; i++) printf(" %d", pi.nConnTable[i]); printf("\n");
                printf("FT %d", pi.lenTautomer); for (int i = 0; i < pi.lenTautomer; i++) printf(" %d", pi.nTautomer[i]); printf("\n");
                printf("FH"); for (int i = 0; i < nA; i++) printf(" %d/%d", pi.nNum_H[i], pi.nNum_H_fixed[i]); printf("\n");
                printf("FI %d %d", pi.nNumberOfIsotopicAtoms, pi.nNumberOfIsotopicTGroups);
                for (int i = 0; i < pi.nNumberOfIsotopicAtoms; i++) printf(" %d,%d,%d,%d,%d", pi.IsotopicAtom[i].nAtomNumber, pi.IsotopicAtom[i].nIsoDifference, pi.IsotopicAtom[i].nNum_H, pi.IsotopicAtom[i].nNum_D, pi.IsotopicAtom[i].nNum_T);
                printf(" |");
                for (int i = 0; i < pi.nNumberOfIsotopicTGroups; i++) printf(" %d,%d,%d,%d", pi.IsotopicTGroup[i].nTGroupNumber, pi.IsotopicTGroup[i].nNum_H, pi.IsotopicTGroup[i].nNum_D, pi.IsotopicTGroup[i].nNum_T);
                printf("\n");
                printf("FL"); for (int i = 0; i < nA + 3; i++) printf(" %d", pi.nPossibleLocationsOfIsotopicH[i]); printf("\n");
                for (int q = 0; q < 2; q++) {
                    INChI_Stereo *s = sts[q];
                    printf("FS%d %d %d %d %d |", q, s->nNumberOfStereoCenters, s->nNumberOfStereoBonds, s->nCompInv2Abs, s->bTrivialInv);
                    for (int i = 0; i < 64; i++) printf(" %d,%d,%d,%d", s->nNumber[i], s->t_parity[i], s->nNumberInv[i], s->t_parityInv[i]);
                    printf(" |");
                    for (int i = 0; i < 64; i++) printf(" %d,%d,%d", s->nBondAtom1[i], s->nBondAtom2[i], s->b_parity[i]);
                    printf("\n");
                }
                printf("FX %d %d %d %d %d", pa.nNumberOfAtoms, pa.nNumberOfTGroups, pa.nNumRemovedProtons, pa.nNumRemovedIsotopicH[0], pa.nNumRemovedIsotopicH[1]);
                printf(" %d\n", pa.nNumRemovedIsotopicH[2]);
                AT_NUMB *aa[8] = {pa.nOrigAtNosInCanonOrd, pa.nIsotopicOrigAtNosInCanonOrd, pa.nOrigAtNosInCanonOrdInv, pa.nIsotopicOrigAtNosInCanonOrdInv,
                                  pa.nConstitEquNumbers, pa.nConstitEquIsotopicNumbers, pa.nConstitEquTGroupNumbers, pa.nConstitEquIsotopicTGroupNumbers};
                for (int q = 0; q < 8; q++) { printf("FQ%d", q); for (int i = 0; i < al; i++) printf(" %d", aa[q][i]); printf("\n"); }
                printf("FO2"); for (int i = 0; i < nA; i++) printf(" %d,%d,%d", pa.OrigInfo[i].cCharge, pa.OrigInfo[i].cRadical, pa.OrigInfo[i].cUnusualValence); printf("\n");
                printf("FE %s\n", se);
                printf("FB"); for (int i = 0; i < nA; i++) printf(" %d/%d", at[i].bAmbiguousStereo, na[i].bAmbiguousStereo); printf("\n");
                pu("FY", cs.nCanonOrdStereo, cs.nCanonOrdStereo ? nat : 0);
                pu("FZ", cs.nCanonOrdIsotopicStereo, cs.nCanonOrdIsotopicStereo ? nat : 0);
                printf("FG"); if (haveTgi) { for (int i = 0; i < 4 * ntgs; i++) printf(" %d", tgi.tGroupNumber[i]); printf(" |"); for (int i = 0; i < nE; i++) printf(" %d", epl[i]); } printf("\n");
            }
        }
    }
    return 0;
}
"""


def stereo_rows(rng, n, edges, valence):
    """per atom [par, par2(unused), sap, sap2, rs, 3 x (nbr+1, ord, z, parity)] (17 ints)"""
    nb = [[] for _ in range(n)]
    for a, b in edges:
        nb[a].append(b)
        nb[b].append(a)
    rows = [[0] * 17 for _ in range(n)]
    used = [False] * n
    cand = list(range(n))
    rng.shuffle(cand)
    for i in cand:
        if used[i]:
            continue
        if len(nb[i]) >= 3 and rng.random() < 0.5:
            rows[i][0] = rng.choice([1, 2, 1, 2, 3, 4])
            rows[i][2] = rng.choice([6, 6, 6, 3, 4])
            used[i] = True
    be = [(a, b) for a, b in edges if len(nb[a]) in (2, 3) and len(nb[b]) in (2, 3)]
    rng.shuffle(be)
    for a, b in be[: rng.randrange(0, 3)]:
        if used[a] or used[b]:
            continue
        par = rng.choice([6, 6, 1, 2, 3, 4])
        z = rng.choice([1, -1])
        rows[a][5:9] = [b + 1, nb[a].index(b), z, par]
        rows[b][5:9] = [a + 1, nb[b].index(a), z, par]
        rows[a][0] = rows[a][0] or rng.choice([1, 2])
        rows[b][0] = rows[b][0] or rng.choice([1, 2])
        used[a] = used[b] = True
    for i in range(n):
        rows[i][4] = rng.randrange(1, 3)
    return rows


def gen(rng):
    c = f2.gen(rng)
    n = c["n"]
    e = c["edges"]
    plain = stereo_rows(rng, n, e, None)
    iso = stereo_rows(rng, n, e, None) if rng.random() < 0.8 else [r[:] for r in plain]
    if rng.random() < 0.4:
        for i in range(n):
            if plain[i][0] and rng.random() < 0.7:
                iso[i] = plain[i][:]
    c["st1"], c["st2"] = plain, iso
    c["modeBits"] = rng.choice([0, 0x08, 0x02, 0x0A, 0x1A, 0x1A, 0x1A, 0x08 | 0x10, 0x02 | 0x08 | 0x10])
    c["uu"] = rng.randrange(2)
    return c


def fmt(c):
    o = f2.fmt(c).rstrip("\n").split("\n")
    # f2.fmt wrote sizes as two lines of 7 ints, which the harness reads the same way
    o += [" ".join(map(str, r)) for r in c["st1"]]
    o += [" ".join(map(str, r)) for r in c["st2"]]
    o.append(f"{c['modeBits']} {c['uu']}")
    return "\n".join(o) + "\n"


def main():
    rng = random.Random(777555)
    cases = [gen(rng) for _ in range(2500)]
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        dropped = 0
        while True:
            stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
            r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, timeout=1800)
            out = r.stdout.decode().strip().split("\n")
            if r.returncode == 0:
                break
            # the C code aborts (heap corruption) on a few random inputs: drop the case
            last = max(i for i, l in enumerate(out) if l.startswith("N "))
            idx = int(out[last].split()[1])
            del cases[idx]
            dropped += 1
    print("dropped", dropped, "cases on which the C code aborts")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_cano2.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
