"""Generates tests/data/inchi_bnchargeh.json from the InChI C code: how the
alternating paths found by the balanced network search change the charges
and H counts of the structure (`MarkAtomsAtTautGroups`,
`EliminatePlusMinusChargeAmbiguity`, `SubtractOrChangeAtHChargeBNS`,
`AddChangedAtHChargeBNS`, `AddOrRemoveExplOrImplH`) on networks with fictitious
group vertices built from random molecules. The function texts are cut out of
ichi_bns.c. Needs a C compiler and the fetched InChI sources."""

import json
import random
import re
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
#define IS_FORBIDDEN(EDGE_FORBIDDEN, PBNS)     (EDGE_FORBIDDEN)
#define prim(v) (Vertex)((v)^1)
#define SwitchEdge_Vert1(u) SwitchEdge[u][0]
#define SwitchEdge_Vert2(u) Get2ndEdgeVertex( pBNS, SwitchEdge[u] )
#define SwitchEdge_IEdge(u) SwitchEdge[u][1]
#define CHECK_TACN  1
#define bInchiTimeIsOver(a, b) 0
@MACROS@
@FUNCS@
static void dump(BN_STRUCT *p) {
    printf("C %d %d\n", p->num_vertices, p->num_edges);
    for (int i = 0; i < p->num_vertices; i++) {
        BNS_VERTEX *v = &p->vert[i];
        printf("V %d %d %d %d %d %d %d %d |", v->type, v->st_edge.cap, v->st_edge.cap0, v->st_edge.flow, v->st_edge.flow0, v->st_edge.pass, v->num_adj_edges, v->max_adj_edges);
        for (int k = 0; k < v->num_adj_edges; k++) printf(" %d", v->iedge[k]);
        printf("\n");
    }
    for (int i = 0; i < p->num_edges; i++) {
        BNS_EDGE *e = &p->edge[i];
        printf("E %d %d %d %d %d %d %d %d %d\n", e->cap, e->cap0, e->flow, e->flow0, e->pass, e->neighbor1, e->neighbor12, e->neigh_ord[0], e->neigh_ord[1]);
    }
}
static unsigned lcg;
static int rnd(int m) { lcg = (lcg * 1103515245u + 12345u) & 0x7fffffffu; return (int)((lcg >> 8) % (unsigned)m); }
static void dumpatoms(inp_ATOM *at, int total, T_GROUP_INFO *t) {
    for (int i = 0; i < total; i++) {
        printf("A %d %d %d %d %d,%d,%d t%d s%d,%d,%d o%d,%d,%d n%d,%d,%d g%d,%d,%d r%d nb%d\n", at[i].charge, at[i].num_H, at[i].radical, at[i].chem_bonds_valence, at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2], at[i].at_type,
               at[i].sb_parity[0], at[i].sb_parity[1], at[i].sb_parity[2], at[i].sb_ord[0], at[i].sb_ord[1], at[i].sb_ord[2],
               at[i].sn_ord[0], at[i].sn_ord[1], at[i].sn_ord[2], at[i].sn_orig_at_num[0], at[i].sn_orig_at_num[1], at[i].sn_orig_at_num[2], at[i].orig_at_number, at[i].neighbor[0]);
    }
    printf("I %d %d %d,%d,%d\n", t->tni.nNumRemovedExplicitH, t->tni.nNumRemovedProtons, t->tni.nNumRemovedProtonsIsotopic[0], t->tni.nNumRemovedProtonsIsotopic[1], t->tni.nNumRemovedProtonsIsotopic[2]);
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, nH, maxadd_a, maxadd_e, maxaltp, bcf, nops, e1, e2, mk, mkseed, walk, seed, npaths;
        scanf("%d %d %d %d %d %d %d %d %d %d %d %d %d", &n, &nH, &maxadd_a, &maxadd_e, &maxaltp, &bcf, &nops, &e1, &e2, &mkseed, &walk, &seed, &npaths);
        inp_ATOM *at = calloc(n + nH + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, vt, chg, nh, i0, i1, i2, orig, sp[3], so[3], sn[3], sg[3];
            scanf("%d %d %d %d %d %d %d %d %d %d %d", &v, &cbv, &el, &rad, &vt, &chg, &nh, &i0, &i1, &i2, &orig);
            scanf("%d %d %d %d %d %d %d %d %d %d %d %d", sp, sp + 1, sp + 2, so, so + 1, so + 2, sn, sn + 1, sn + 2, sg, sg + 1, sg + 2);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].radical = rad; at[i].c_point = vt;
            at[i].charge = chg; at[i].num_H = nh; at[i].num_iso_H[0] = i0; at[i].num_iso_H[1] = i1; at[i].num_iso_H[2] = i2; at[i].orig_at_number = orig;
            for (int k = 0; k < 3; k++) { at[i].sb_parity[k] = sp[k]; at[i].sb_ord[k] = so[k]; at[i].sn_ord[k] = sn[k]; at[i].sn_orig_at_num[k] = sg[k]; }
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        for (int i = 0; i < nH; i++) {
            int nb0, orig;
            scanf("%d %d", &nb0, &orig);
            at[n + i].neighbor[0] = nb0; at[n + i].orig_at_number = orig; at[n + i].el_number = 1; at[n + i].valence = 1;
        }
        int nchg = 0;
        BN_STRUCT *pBNS = AllocateAndInitBnStruct(at, n, maxadd_a, maxadd_e, maxaltp, &nchg);
        for (int i = 0; i < n; i++) pBNS->vert[i].type |= (at[i].c_point & 0x1c) << 1 | (at[i].c_point & 3) << 9;
        printf("N %d\n", c);
        for (int k = 0; k < nops; k++) {
            int op, a, b, x, y, z;
            scanf("%d %d %d %d %d %d", &op, &a, &b, &x, &y, &z);
            int dots = 5;
            int nv = pBNS->num_vertices;
            if (op == 10) {
                if (a >= nv) continue;
                bAddNewVertex(pBNS, a, x, y, z, &dots);
            } else if (op == 11) {
                if (a < pBNS->num_atoms || a >= nv) continue;
                pBNS->vert[a].type = x;
            } else {
                if (a >= nv || b >= nv || a == b) continue;
                AddNewEdge(pBNS->vert + a, pBNS->vert + b, pBNS, x, y);
            }
        }
        BN_DATA *pBD = AllocateAndInitBnData(pBNS->max_vertices);
        int ret = 0;
        if (!walk) {
            ret = RunBalancedNetworkSearch(pBNS, pBD, bcf);
        } else {
            lcg = seed;
            pBNS->num_altp = 0;
            for (int pass = 0; pass < npaths && pass < pBNS->max_altp; pass++) {
                BNS_ALT_PATH *ap = pBNS->altp[pass];
                int v = -1;
                for (int t = 0; t < 20; t++) { int c2 = rnd(pBNS->num_vertices); if (pBNS->vert[c2].num_adj_edges > 0) { v = c2; break; } }
                if (v < 0) continue;
                int L = 1 + rnd(10);
                if (L > pBNS->max_vertices) L = pBNS->max_vertices;
                int start = v;
                int delta = rnd(2) ? 1 : -1;
                int k = 0;
                for (; k < L; k++) {
                    BNS_VERTEX *pv = &pBNS->vert[v];
                    if (pv->num_adj_edges == 0) break;
                    int o = rnd(pv->num_adj_edges);
                    BNS_EDGE *e = &pBNS->edge[pv->iedge[o]];
                    int v2 = e->neighbor12 ^ v;
                    int o2 = (v2 == e->neighbor1) ? e->neigh_ord[0] : e->neigh_ord[1];
                    ap[iALTP_NEIGHBOR + k].ineigh[0] = o;
                    ap[iALTP_NEIGHBOR + k].ineigh[1] = o2;
                    v = v2;
                }
                if (k == 0) continue;
                ap[iALTP_START_ATOM].number = start;
                ap[iALTP_END_ATOM].number = v;
                ap[iALTP_PATH_LEN].number = k;
                ap[iALTP_FLOW].flow[0] = delta;
                pBNS->num_altp++;
            }
        }
        printf("X %d %d %d\n", ret, pBNS->num_altp, pBNS->num_vertices);
        for (int a = 0; a < pBNS->num_altp && a < pBNS->max_altp; a++) {
            BNS_ALT_PATH *ap = pBNS->altp[a];
            printf("P %d %d %d %d", ap[iALTP_START_ATOM].number, ap[iALTP_END_ATOM].number, ap[iALTP_PATH_LEN].number, ap[iALTP_FLOW].flow[0]);
            for (int k = 0; k < ap[iALTP_PATH_LEN].number; k++) printf(" %d,%d", ALTP_THIS_ATOM_NEIGHBOR(ap, k), ALTP_NEXT_ATOM_NEIGHBOR(ap, k));
            printf("\n");
        }
        BN_AATG aatg;
        memset(&aatg, 0, sizeof(aatg));
        aatg.nMarkedAtom = calloc(2 * pBNS->max_vertices + 8, 1);
        for (int i = 0; i < n; i++) aatg.nMarkedAtom[i] = (mkseed >> (i % 8)) & 1;
        int r1 = MarkAtomsAtTautGroups(pBNS, n, &aatg, e1, e2);
        printf("M %d %d |", r1, aatg.nNumFound);
        for (int i = 0; i < pBNS->num_vertices; i++) printf(" %d", aatg.nMarkedAtom[i]);
        printf("\n");
        int r2 = EliminatePlusMinusChargeAmbiguity(pBNS, n);
        printf("E %d\n", r2);
        printf("V");
        for (int i = 0; i < pBNS->num_vertices; i++) printf(" %d:%d:%d", pBNS->vert[i].type, pBNS->vert[i].st_edge.cap, pBNS->vert[i].st_edge.flow);
        printf("\nF");
        for (int i = 0; i < pBNS->num_edges; i++) printf(" %d", pBNS->edge[i].flow);
        printf("\n");
        T_GROUP_INFO tgi;
        memset(&tgi, 0, sizeof(tgi));
        tgi.tni.nNumRemovedExplicitH = nH;
        tgi.tni.nNumRemovedProtons = 1;
        S_CHAR *mark2 = calloc(n + 1, 1);
        int totals[ATTOT_ARRAY_LEN];
        memset(totals, 0, sizeof(totals));
        int r3 = SubtractOrChangeAtHChargeBNS(pBNS, at, n, totals, mark2, NULL, 1);
        printf("S %d |", r3);
        for (int i = 0; i < n; i++) printf(" %d", mark2[i]);
        printf(" |");
        for (int i = 0; i < ATTOT_ARRAY_LEN; i++) printf(" %d", totals[i]);
        printf("\n");
        int r4 = AddChangedAtHChargeBNS(at, n, totals, mark2);
        printf("D %d |", r4);
        for (int i = 0; i < ATTOT_ARRAY_LEN; i++) printf(" %d", totals[i]);
        printf("\n");
        int r5 = SubtractOrChangeAtHChargeBNS(pBNS, at, n, NULL, NULL, &tgi, 0);
        printf("C %d\n", r5);
        dumpatoms(at, n + nH, &tgi);
        for (int k = 0; k < 3; k++) {
            int delta, no;
            scanf("%d %d", &delta, &no);
            int r6 = AddOrRemoveExplOrImplH(delta, at, n, no, &tgi);
            printf("H %d\n", r6);
            dumpatoms(at, n + nH, &tgi);
        }
        DeAllocateBnData(pBD);
        DeAllocateBnStruct(pBNS);
        free(at);
    }
    return 0;
}
"""

ELS = [6, 6, 6, 7, 7, 8, 8, 16, 15, 17, 1]


def gen(rng):
    n = rng.randrange(2, 9)
    nb = [[] for _ in range(n)]
    for i in range(1, n):
        if rng.random() < 0.95:
            j = rng.randrange(i)
            t = rng.choice([1, 1, 2, 3, 4, 1, 2, 0x11, 9, 8])
            nb[i].append((j, t))
            nb[j].append((i, t))
    for _ in range(rng.randrange(0, 3)):
        a, b = rng.randrange(n), rng.randrange(n)
        if a != b and all(x != b for x, _ in nb[a]):
            t = rng.choice([1, 2, 4])
            nb[a].append((b, t))
            nb[b].append((a, t))
    atoms = []
    for i in range(n):
        v = len(nb[i])
        cbv = sum(min(t & 15, 3) for _, t in nb[i]) + rng.choice([0, 0, 1, 1, 2])
        atoms.append({"v": v, "cbv": cbv, "el": rng.choice(ELS), "rad": rng.choice([0] * 6 + [2, 1, 3]),
                      "vt": rng.randrange(0, 32), "nb": nb[i]})
    return atoms


def main():
    rng = random.Random(20)
    cases = []
    for _ in range(2300):
        atoms = gen(rng)
        n = len(atoms)
        nh = rng.randrange(0, 4)
        ma, me = rng.randrange(1, 5), rng.randrange(1, 6)
        nv = n + ma
        origs = [i + 1 for i in range(n)] + [100 + j for j in range(nh)] + [0]
        for i, a in enumerate(atoms):
            a["chg"] = rng.choice([0, 0, 0, 1, -1])
            a["nh"] = rng.choice([0, 1, 2, 3])
            a["iso"] = [rng.choice([0, 0, 1]) for _ in range(3)]
            a["orig"] = i + 1
            k = rng.randrange(0, 4)
            a["sbp"] = [rng.choice([1, 2, 3, 4]) if j < k else 0 for j in range(3)]
            v = max(a["v"], 1)
            a["sbo"] = [rng.randrange(v) for _ in range(3)]
            a["sno"] = [rng.randrange(-1, 2) for _ in range(3)]
            a["sng"] = [rng.choice(origs) for _ in range(3)]
        hs = [(rng.randrange(n), 100 + j) for j in range(nh)]
        ops = []
        for _ in range(rng.randrange(2, 12)):
            op = rng.choice([10, 10, 11, 11, 12, 12])
            if op == 10:
                ops.append((10, rng.randrange(n), 0, rng.randrange(0, 3), rng.randrange(0, 2), rng.randrange(1, 5)))
            elif op == 11:
                ops.append((11, n + rng.randrange(ma), 0, rng.choice([4, 0x10, 0x110, 0x10, 0x40, 4]), 0, 0))
            else:
                ops.append((12, rng.randrange(nv), rng.randrange(nv), rng.randrange(0, 3), rng.randrange(0, 2), 0))
        cases.append({"atoms": atoms, "hs": hs, "ma": ma, "me": me, "mp": rng.randrange(1, 5), "bcf": rng.choice([1, 1, 1, 5, 13]),
                      "walk": rng.random() < 0.9, "seed": rng.randrange(1 << 30), "np": rng.randrange(1, 4), "ops": ops, "e1": rng.randrange(-1, nv), "e2": rng.randrange(-1, nv), "mk": rng.randrange(256),
                      "hc": [(rng.randrange(-3, 4), rng.randrange(n)) for _ in range(3)]})
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{len(c['atoms'])} {len(c['hs'])} {c['ma']} {c['me']} {c['mp']} {c['bcf']} {len(c['ops'])} {c['e1']} {c['e2']} {c['mk']} {int(c['walk'])} {c['seed']} {c['np']}\n"
        for a in c["atoms"]:
            stdin += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['vt']} {a['chg']} {a['nh']} {a['iso'][0]} {a['iso'][1]} {a['iso'][2]} {a['orig']}\n"
            stdin += " ".join(map(str, a["sbp"] + a["sbo"] + a["sno"] + a["sng"])) + "\n"
            stdin += " ".join(f"{x} {t}" for x, t in a["nb"]) + "\n"
        for h in c["hs"]:
            stdin += f"{h[0]} {h[1]}\n"
        for o in c["ops"]:
            stdin += " ".join(map(str, o)) + "\n"
        for d, no in c["hc"]:
            stdin += f"{d} {no}\n"
    text = (SRC / "ichi_bns.c").read_text(errors="replace")

    def cut(start):
        a = text.rindex(start)
        return text[a:text.index("\n}\n", a) + 3]

    names = [
        "BN_STRUCT* AllocateAndInitBnStruct( inp_ATOM *at,\n", "BN_STRUCT* DeAllocateBnStruct( BN_STRUCT *pBNS )\n{",
        "int ReInitBnStructAltPaths( BN_STRUCT *pBNS )\n{",
        "void ClearAllBnDataVertices( Vertex *v, Vertex value, int size )\n{", "void ClearAllBnDataEdges( Edge *e, Vertex value, int size )\n{",
        "BN_DATA *DeAllocateBnData( BN_DATA *pBD )\n{", "BN_DATA *AllocateAndInitBnData( int max_num_vertices )\n{",
        "int ReInitBnData( BN_DATA *pBD )\n{", "int GetVertexDegree( BN_STRUCT* pBNS, Vertex v )\n{",
        "Vertex Get2ndEdgeVertex( BN_STRUCT* pBNS, Edge uv )\n{", "Vertex GetVertexNeighbor( BN_STRUCT* pBNS,\n",
        "int GetEdgePointer( BN_STRUCT* pBNS,\n", "int AugmentEdge( BN_STRUCT* pBNS,\n", "int rescap_mark( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv )\n{",
        "Vertex GetPrevVertex( BN_STRUCT* pBNS, Vertex y, Edge *SwitchEdge, EdgeIndex *iuv )\n{",
        "int bIgnoreVertexNonTACN_atom( BN_STRUCT* pBNS, Vertex u, Vertex v )\n{", "int bIgnoreVertexNonTACN_group( BN_STRUCT* pBNS,\n",
        "int rescap( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv )\n{",
        "int bIsBnsEndpoint( BN_STRUCT *pBNS, int v )\n{", "int bRadChangesAtomType( BN_STRUCT *pBNS,\n",
        "int RegisterRadEndpoint( BN_STRUCT *pBNS, BN_DATA *pBD, Vertex u )\n{",
        "Vertex FindBase( Vertex u, Vertex *BasePtr )\n{", "int FindPathToVertex_s( Vertex x,\n", "Vertex MakeBlossom( BN_STRUCT* pBNS,\n",
        "int PullFlow( BN_STRUCT *pBNS,\n", "int FindPathCap( BN_STRUCT* pBNS, Edge *SwitchEdge, Vertex x, Vertex y, int delta )\n{",
        "int BalancedNetworkSearch( BN_STRUCT* pBNS, BN_DATA *pBD, int bChangeFlow )\n",
        "int RunBalancedNetworkSearch( BN_STRUCT *pBNS, BN_DATA *pBD, int bChangeFlow )\n{",
        "int bAddNewVertex( BN_STRUCT *pBNS,\n", "int AddNewEdge( BNS_VERTEX *p1,\n", "BNS_IEDGE GetEdgeToGroupVertex( BN_STRUCT *pBNS, Vertex v1, AT_NUMB type )\n{",
        "Vertex GetGroupVertex( BN_STRUCT *pBNS, Vertex v1, AT_NUMB type )\n{", "int bAddStCapToAVertex( BN_STRUCT *pBNS,\n",
        "int bRestoreFlowAfterCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd )\n{", "int bSetFlowToCheckOneBond( BN_STRUCT *pBNS,\n",
        "int bRestoreBnsAfterCheckAltPath( BN_STRUCT *pBNS,\n", "int bSetBnsToCheckAltPath( BN_STRUCT *pBNS,\n",
        "int MarkAtomsAtTautGroups( BN_STRUCT *pBNS,\n", "int EliminatePlusMinusChargeAmbiguity( BN_STRUCT *pBNS, int num_atoms )\n{",
        "int AddChangedAtHChargeBNS( inp_ATOM *at,\n", "int AddOrRemoveExplOrImplH( int nDelta,\n", "int SubtractOrChangeAtHChargeBNS( BN_STRUCT *pBNS,\n",
        "int GetAtomChargeType( inp_ATOM *atom,\n", "int IsZOX( inp_ATOM *atom, int at_x, int ord )\n{", "int is_Z_atom( U_CHAR el_number )\n{"]
    st = (SRC / "ichister.c").read_text(errors="replace")
    a2 = st.rindex("int get_opposite_sb_atom(")
    opp = st[a2:st.index("\n}\n", a2) + 3]
    funcs = [cut(n) for n in names] + [opp]
    taut = (SRC / "ichitaut.c").read_text(errors="replace")
    a = taut.rindex("int is_centerpoint_elem( U_CHAR el_number )\n{")
    funcs.append(taut[a:taut.index("\n}\n", a) + 3])
    a = text.index("#define TREE_NOT_IN_M")
    b = text.index("#define TREE_MARK")
    macros = text[a:text.index("\n", b)]
    ea = text.index("typedef enum tagAtTypeTotals")
    eb = text.index("#define ATBIT_Errors")
    macros += "\n" + text[ea:text.index("\n", eb)] + "\n#define MAX_ALT_AATG_ARRAY_LEN 127\n#define AATG_MARK_IN_PATH 1\n"
    alt = re.search(r"typedef struct tagAltPathChanges.*?\} ALT_PATH_CHANGES;", text, re.S).group(0)
    protos = """int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt, int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );
int MarkAtomsAtTautGroups( BN_STRUCT *pBNS, int num_atoms, BN_AATG *pAATG, int nEnd1, int nEnd2 );
int EliminatePlusMinusChargeAmbiguity( BN_STRUCT *pBNS, int num_atoms );
int AddChangedAtHChargeBNS( inp_ATOM *at, int num_atoms, int nAtTypeTotals[], S_CHAR *mark );
int AddOrRemoveExplOrImplH( int nDelta, inp_ATOM *at, int num_atoms, AT_NUMB at_no, T_GROUP_INFO *t_group_info );
int SubtractOrChangeAtHChargeBNS( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int nAtTypeTotals[], S_CHAR *mark, T_GROUP_INFO *t_group_info, int bSubtract );
int GetAtomChargeType( inp_ATOM *atom, int at_no, int nAtTypeTotals[], int *pMask, int bSubtract );
int IsZOX( inp_ATOM *atom, int at_x, int ord );
int is_Z_atom( U_CHAR el_number );
int bRestoreFlowAfterCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd );
int bSetFlowToCheckOneBond( BN_STRUCT *pBNS, int iedge, int flow, BNS_FLOW_CHANGES *fcd );
int bRestoreBnsAfterCheckAltPath( BN_STRUCT *pBNS, ALT_PATH_CHANGES *apc, int bChangeFlow );
int bSetBnsToCheckAltPath( BN_STRUCT *pBNS, int nVertDoubleBond, int nVertSingleBond, AT_NUMB type, int path_type, ALT_PATH_CHANGES *apc, BNS_FLOW_CHANGES *fcd, int *nDots );
int bAddNewVertex( BN_STRUCT *pBNS, int nVertDoubleBond, int nCap, int nFlow, int nMaxAdjEdges, int *nDots );
int AddNewEdge( BNS_VERTEX *p1, BNS_VERTEX *p2, BN_STRUCT *pBNS, int nEdgeCap, int nEdgeFlow );
BNS_IEDGE GetEdgeToGroupVertex( BN_STRUCT *pBNS, Vertex v1, AT_NUMB type );
Vertex GetGroupVertex( BN_STRUCT *pBNS, Vertex v1, AT_NUMB type );
int bAddStCapToAVertex( BN_STRUCT *pBNS, Vertex v1, Vertex v2, VertexFlow *nOldCapVertSingleBond, int *nDots, int bAdjacentDonors );
int rescap( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv );
BN_STRUCT* DeAllocateBnStruct( BN_STRUCT *pBNS );
BN_DATA *DeAllocateBnData( BN_DATA *pBD );
int ReInitBnStructAltPaths( BN_STRUCT *pBNS );
int GetVertexDegree( BN_STRUCT* pBNS, Vertex v );
Vertex GetVertexNeighbor( BN_STRUCT* pBNS, Vertex v, int neigh, EdgeIndex *iedge );
Vertex Get2ndEdgeVertex( BN_STRUCT* pBNS, Edge uv );
int GetEdgePointer( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv, BNS_EDGE **uv, S_CHAR *s_or_t );
int AugmentEdge( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv, int delta, S_CHAR bReverse, int bChangeFlow );
Vertex GetPrevVertex( BN_STRUCT* pBNS, Vertex y, Edge *SwitchEdge, EdgeIndex *iuv );
int bRadChangesAtomType( BN_STRUCT *pBNS, BN_DATA *pBD, Vertex v, Vertex v_1, Vertex v_2 );
Vertex FindBase( Vertex u, Vertex *BasePtr );
int FindPathToVertex_s( Vertex x, Edge *SwitchEdge, Vertex *BasePtr, Vertex *Path, int MaxPathLen );
int PullFlow( BN_STRUCT *pBNS, Edge *SwitchEdge, Vertex x, Vertex y, int delta, S_CHAR bReverse, int bChangeFlow );
int FindPathCap( BN_STRUCT* pBNS, Edge *SwitchEdge, Vertex x, Vertex y, int delta );
int RegisterRadEndpoint( BN_STRUCT *pBNS, BN_DATA *pBD, Vertex u );
int BalancedNetworkSearch( BN_STRUCT* pBNS, BN_DATA *pBD, int bChangeFlow );
Vertex MakeBlossom( BN_STRUCT* pBNS, Vertex *ScanQ, int *pQSize, Vertex *Pu, Vertex *Pv, int max_len_Pu_Pv, Edge *SwitchEdge, Vertex *BasePtr, Vertex u, Vertex v, EdgeIndex iuv, Vertex b_u, Vertex b_v, S_CHAR *Tree );
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        code = HARNESS.replace("@MACROS@", alt + "\n" + macros + "\n#define BNS_CHK_ALTP_NO_ALTPATH 0\n#define BNS_CHK_ALTP_SAME_TGROUP 1\n#define BNS_CHK_ALTP_SAME_VERTEX 2\n#define BNS_CHK_ALTP_SET_SUCCESS 4\n").replace("@FUNCS@", protos + "\n".join(funcs))
        (Path(d) / "h.c").write_text(code)
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c"), "-lm"], check=True)
        # a few random cases hit undefined behaviour of the C code (reads/writes
        # past an array); they are dropped
        def chunk(c):
            t = f"1\n{len(c['atoms'])} {len(c['hs'])} {c['ma']} {c['me']} {c['mp']} {c['bcf']} {len(c['ops'])} {c['e1']} {c['e2']} {c['mk']} {int(c['walk'])} {c['seed']} {c['np']}\n"
            for a in c["atoms"]:
                t += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['vt']} {a['chg']} {a['nh']} {a['iso'][0]} {a['iso'][1]} {a['iso'][2]} {a['orig']}\n"
                t += " ".join(map(str, a["sbp"] + a["sbo"] + a["sno"] + a["sng"])) + "\n"
                t += " ".join(f"{x} {t2}" for x, t2 in a["nb"]) + "\n"
            for h in c["hs"]:
                t += f"{h[0]} {h[1]}\n"
            for o in c["ops"]:
                t += " ".join(map(str, o)) + "\n"
            for d, no in c["hc"]:
                t += f"{d} {no}\n"
            return t
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True)
        if r.returncode == 0:
            out = r.stdout.decode().strip().split("\n")
        else:
            kept, lines = [], []
            for c in cases:
                r1 = subprocess.run([str(exe)], input=chunk(c).encode(), capture_output=True)
                if r1.returncode == 0:
                    kept.append(c)
                    lines += r1.stdout.decode().strip().split("\n")[0:0] + [f"N {len(kept) - 1}"] + r1.stdout.decode().strip().split("\n")[1:]
            print(len(cases) - len(kept), "cases dropped")
            cases = kept
            out = lines
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_bnchargeh.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
