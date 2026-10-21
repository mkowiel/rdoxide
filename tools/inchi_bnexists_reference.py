"""Generates tests/data/inchi_bnexists.json from the InChI C code:
`bExistsAltPath` and `bExistsAnyAltPath` of ichi_bns.c (whether H or a charge
can move between two atoms along an alternating path) on networks with
fictitious group vertices built from random molecules. The function texts are
cut out of the InChI sources. Needs a C compiler and the fetched InChI
sources."""

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
#include "ichicano.h"
#include "ichicant.h"
#define LOG_NO_ARGS(m)
#define LOG_MULT_ARGS(f, ...)
#define IS_FORBIDDEN(EDGE_FORBIDDEN, PBNS)     (EDGE_FORBIDDEN)
#define prim(v) (Vertex)((v)^1)
#define SwitchEdge_Vert1(u) SwitchEdge[u][0]
#define SwitchEdge_Vert2(u) Get2ndEdgeVertex( pBNS, SwitchEdge[u] )
#define SwitchEdge_IEdge(u) SwitchEdge[u][1]
#define CHECK_TACN  1
#define bInchiTimeIsOver(a, b) 0
#define MAX_NUM_RAD 256
#define BNS_MAX_NUM_FLOW_CHANGES (1+2*MAX_BOND_EDGE_CAP)
int SetBitCreate( CANON_GLOBALS *pCG ) { return 0; }
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
static void dumpstate(BN_STRUCT *p, inp_ATOM *at, int n) {
    printf("C %d %d %d %d\n", p->num_vertices, p->num_edges, p->tot_st_cap, p->tot_st_flow);
    printf("V");
    for (int i = 0; i < p->num_vertices; i++) printf(" %d:%d:%d:%d:%d:%d", p->vert[i].type, p->vert[i].st_edge.cap, p->vert[i].st_edge.cap0, p->vert[i].st_edge.flow, p->vert[i].st_edge.flow0, p->vert[i].num_adj_edges);
    printf("\nF");
    for (int i = 0; i < p->num_edges; i++) printf(" %d:%d:%d:%d", p->edge[i].cap, p->edge[i].cap0, p->edge[i].flow, p->edge[i].flow0);
    printf("\nA");
    for (int i = 0; i < n; i++) {
        printf(" %d,%d,%d,%d,%d,%d,%d,%d,%d:", at[i].charge, at[i].num_H, at[i].radical, at[i].chem_bonds_valence, at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2], at[i].at_type, at[i].endpoint);
        for (int j = 0; j < at[i].valence; j++) printf("%d.", at[i].bond_type[j]);
    }
    printf("\n");
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, nH, maxadd_a, maxadd_e, maxaltp, nops, ncalls, tflags;
        scanf("%d %d %d %d %d %d %d %d", &n, &nH, &maxadd_a, &maxadd_e, &maxaltp, &nops, &ncalls, &tflags);
        inp_ATOM *at = calloc(n + nH + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, vt, chg, nh, i0, i1, i2, orig, ep, att, ring, sp[3], so[3], sn[3], sg[3];
            scanf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d", &v, &cbv, &el, &rad, &vt, &chg, &nh, &i0, &i1, &i2, &orig, &ep, &att, &ring);
            scanf("%d %d %d %d %d %d %d %d %d %d %d %d", sp, sp + 1, sp + 2, so, so + 1, so + 2, sn, sn + 1, sn + 2, sg, sg + 1, sg + 2);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].radical = rad; at[i].c_point = vt;
            at[i].charge = chg; at[i].num_H = nh; at[i].num_iso_H[0] = i0; at[i].num_iso_H[1] = i1; at[i].num_iso_H[2] = i2; at[i].orig_at_number = orig;
            at[i].endpoint = ep; at[i].at_type = att; at[i].nNumAtInRingSystem = ring; at[i].nRingSystem = i % 2;
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
        INCHI_MODE tfl = tflags;
        pBNS->pbTautFlags = &tfl;
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
        pBD->RadEndpoints = calloc(2 * pBNS->max_vertices + 8, sizeof(Vertex));
        pBD->RadEdges = calloc(2 * pBNS->max_vertices + 8, sizeof(EdgeIndex));
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        cg.m_num_bit = 16;
        cg.m_bBit = calloc(16, sizeof(bitWord));
        for (int i = 0; i < 16; i++) cg.m_bBit[i] = 1 << i;
        dumpstate(pBNS, at, n);
        for (int k = 0; k < ncalls; k++) {
            int vd, vs, pt, any, am;
            scanf("%d %d %d %d %d", &vd, &vs, &pt, &any, &am);
            int ret;
            BN_AATG aatg;
            memset(&aatg, 0, sizeof(aatg));
            T_GROUP_INFO tgi;
            memset(&tgi, 0, sizeof(tgi));
            tgi.tni.nNumRemovedExplicitH = nH;
            int totals[ATTOT_ARRAY_LEN];
            memset(totals, 0, sizeof(totals));
            if (am > 0) aatg.nMarkedAtom = calloc(2 * pBNS->max_vertices + 8, 1);
            if (am > 1) { aatg.nAtTypeTotals = totals; aatg.t_group_info = &tgi; }
            if (any) ret = bExistsAnyAltPath(&cg, pBNS, pBD, at, n, vd, vs, pt);
            else ret = bExistsAltPath(&cg, pBNS, pBD, am ? &aatg : NULL, at, n, vd, vs, pt);
            printf("R %d %d %d\n", ret, aatg.nNumFound, tgi.tni.nNumRemovedProtons);
            if (am > 0) { printf("K"); for (int i = 0; i < n; i++) printf(" %d", aatg.nMarkedAtom[i]); printf("\n"); }
            if (am > 1) { printf("T"); for (int i = 0; i < ATTOT_ARRAY_LEN; i++) printf(" %d", totals[i]); printf("\n"); }
            dumpstate(pBNS, at, n);
            free(aatg.nMarkedAtom);
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


ELS2 = [6, 6, 7, 7, 8, 8, 8, 15, 16, 34, 52, 33, 51, 17, 35, 53, 1, 9]
MODES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]


def main():
    rng = random.Random(21)
    cases = []
    for _ in range(3000):
        atoms = gen(rng)
        n = len(atoms)
        nh = rng.randrange(0, 3)
        ma, me = rng.randrange(1, 5), rng.randrange(1, 6)
        nv = n + ma
        for i, a in enumerate(atoms):
            a["el"] = rng.choice(ELS2)
            a["chg"] = rng.choice([0, 0, 0, 1, -1])
            a["nh"] = rng.choice([0, 1, 1, 2])
            a["iso"] = [rng.choice([0, 0, 0, 1]) for _ in range(3)]
            a["orig"] = i + 1
            a["ep"] = rng.choice([0, 0, 0, 1, 2])
            a["att"] = rng.choice([0, 0, 1, 1, 2])
            a["ring"] = rng.choice([0, 1, 5, 6])
            a["sbp"] = [0, 0, 0]
            a["sbo"] = [0, 0, 0]
            a["sno"] = [0, 0, 0]
            a["sng"] = [0, 0, 0]
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
        calls = [(rng.randrange(n), rng.randrange(n), rng.choice(MODES), int(rng.random() < 0.2), rng.choice([0, 0, 1, 2]))
                 for _ in range(rng.randrange(1, 5))]
        cases.append({"atoms": atoms, "hs": hs, "ma": ma, "me": me, "mp": rng.randrange(1, 5), "tf": rng.choice([0, 0x10, 0x10]),
                      "ops": ops, "calls": calls})
    def fmt(c):
        t = f"{len(c['atoms'])} {len(c['hs'])} {c['ma']} {c['me']} {c['mp']} {len(c['ops'])} {len(c['calls'])} {c['tf']}\n"
        for a in c["atoms"]:
            t += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['vt']} {a['chg']} {a['nh']} {a['iso'][0]} {a['iso'][1]} {a['iso'][2]} {a['orig']} {a['ep']} {a['att']} {a['ring']}\n"
            t += " ".join(map(str, a["sbp"] + a["sbo"] + a["sno"] + a["sng"])) + "\n"
            t += " ".join(f"{x} {t2}" for x, t2 in a["nb"]) + "\n"
        for h in c["hs"]:
            t += f"{h[0]} {h[1]}\n"
        for o in c["ops"]:
            t += " ".join(map(str, o)) + "\n"
        for cl in c["calls"]:
            t += " ".join(map(str, cl)) + "\n"
        return t

    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
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
        "int GetAtomChargeType( inp_ATOM *atom,\n", "int IsZOX( inp_ATOM *atom, int at_x, int ord )\n{", "int is_Z_atom( U_CHAR el_number )\n{",
        "int RemoveRadEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at )\n{\n    BNS_EDGE   *e;", "int RestoreRadicalsOnly( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at )\n{\n    BNS_EDGE   *e;",
        "int SetRadEndpoints2( CANON_GLOBALS *pCG,\n", "int cmp_rad_endpoints( const void *a1, const void *a2 )\n{",
        "int RestoreEdgeFlow( BNS_EDGE *edge, int delta, int bChangeFlow )\n{", "int SetAtomBondType( BNS_EDGE *edge,\n",
        "int SetAtomRadAndChemValFromVertexCapFlow( BN_STRUCT *pBNS,\n", "int SetBondsFromBnStructFlow( BN_STRUCT *pBNS,\n",
        "int RestoreBnStructFlow( BN_STRUCT *pBNS, int bChangeFlow )\n{", "int bSetBondsAfterCheckOneBond( BN_STRUCT *pBNS,\n",
        "int bExistsAltPath( CANON_GLOBALS *pCG,\n", "int bExistsAnyAltPath( CANON_GLOBALS *pCG,\n"]
    ican = (SRC / "ichican2.c").read_text(errors="replace")
    def cut2(start):
        a = ican.rindex(start)
        return ican[a:ican.index("\n}\n", a) + 3]
    nodefuncs = [cut2(x) for x in ["int NodeSetCreate( struct tagCANON_GLOBALS *pCG,\n", "void NodeSetFree( struct tagCANON_GLOBALS *pCG, NodeSet *pSet )\n{",
        "void NodeSetFromRadEndpoints( CANON_GLOBALS *pCG,\n", "void RemoveFromNodeSet( CANON_GLOBALS *pCG,\n", "int DoNodeSetsIntersect( NodeSet *cur_nodes, int k1, int k2 )\n{",
        "int IsNodeSetEmpty( NodeSet *cur_nodes, int k )\n{", "void AddNodeSet2ToNodeSet1( NodeSet *cur_nodes, int k1, int k2 )\n{", "int AddNodesToRadEndpoints( CANON_GLOBALS *pCG,\n"]]
    tt = (SRC / "ichitaut.c").read_text(errors="replace")
    def cut3(start):
        a = tt.rindex(start)
        return tt[a:tt.index("\n}\n", a) + 3]
    ctype = tt[tt.index("typedef struct tagChargeType"):tt.index("#define NUM_C_TYPES")] + "#define NUM_C_TYPES  (int)(sizeof( CType )/sizeof(CType[0]))\n"
    endfuncs = [ctype, cut3("int is_centerpoint_elem_KET( U_CHAR el_number )\n{"), cut3("int is_centerpoint_elem_strict( U_CHAR el_number )\n{"),
        cut3("int bCanBeACPoint( inp_ATOM *at,\n"), cut3("int GetChargeType( inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype )\n{")]
    for v in ["", "_PT_22_00", "_PT_16_00", "_PT_06_00", "_PT_39_00", "_PT_13_00", "_PT_18_00"]:
        sig = "int nGetEndpointInfo%s(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)" % v if v else "int nGetEndpointInfo( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif )"
        endfuncs.append(cut3(sig))
    endfuncs.append(cut3("int nGetEndpointInfo_KET( inp_ATOM *atom,\n"))
    st = (SRC / "ichister.c").read_text(errors="replace")
    a2 = st.rindex("int get_opposite_sb_atom(")
    opp = st[a2:st.index("\n}\n", a2) + 3]
    funcs = endfuncs[:1] + [cut(n) for n in names] + [opp] + nodefuncs + endfuncs[1:]
    taut = (SRC / "ichitaut.c").read_text(errors="replace")
    a = taut.rindex("int is_centerpoint_elem( U_CHAR el_number )\n{")
    funcs.append(taut[a:taut.index("\n}\n", a) + 3])
    a = text.index("#define TREE_NOT_IN_M")
    b = text.index("#define TREE_MARK")
    macros = text[a:text.index("\n", b)]
    ea = text.index("typedef enum tagAtTypeTotals")
    eb = text.index("#define ATBIT_Errors")
    macros += "\n#define BNS_CHK_ALTP_NO_ALTPATH 0\n#define BNS_CHK_ALTP_SAME_TGROUP 1\n#define BNS_CHK_ALTP_SAME_VERTEX 2\n#define BNS_CHK_ALTP_SET_SUCCESS 4\n" + text[ea:text.index("\n", eb)] + "\n#define MAX_ALT_AATG_ARRAY_LEN 127\n#define AATG_MARK_IN_PATH 1\n"
    alt = re.search(r"typedef struct tagAltPathChanges.*?\} ALT_PATH_CHANGES;", text, re.S).group(0)
    protos = """int bExistsAltPath( CANON_GLOBALS *pCG, BN_STRUCT *pBNS, BN_DATA *pBD, BN_AATG *pAATG, inp_ATOM *at, int num_atoms, int nVertDoubleBond, int nVertSingleBond, int path_type );
int bExistsAnyAltPath( CANON_GLOBALS *pCG, BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms, int nVert2, int nVert1, int path_type );
int RemoveRadEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at );
int RestoreRadicalsOnly( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at );
int SetRadEndpoints2( CANON_GLOBALS *pCG, BN_STRUCT *pBNS, BN_DATA *pBD, BRS_MODE bRadSrchMode );
int cmp_rad_endpoints( const void *a1, const void *a2 );
int RestoreEdgeFlow( BNS_EDGE *edge, int delta, int bChangeFlow );
int SetAtomBondType( BNS_EDGE *edge, U_CHAR *bond_type12, U_CHAR *bond_type21, int delta, int bChangeFlow );
int SetAtomRadAndChemValFromVertexCapFlow( BN_STRUCT *pBNS, inp_ATOM *atom, int v1 );
int SetBondsFromBnStructFlow( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bChangeFlow0 );
int RestoreBnStructFlow( BN_STRUCT *pBNS, int bChangeFlow );
int bSetBondsAfterCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd, int nTestFlow, inp_ATOM *at, int num_atoms, int bChangeFlow );
int BalancedNetworkSearch( BN_STRUCT* pBNS, BN_DATA *pBD, int bChangeFlow );
int RegisterRadEndpoint( BN_STRUCT *pBNS, BN_DATA *pBD, Vertex u );
int bRadChangesAtomType( BN_STRUCT *pBNS, BN_DATA *pBD, Vertex v, Vertex v_1, Vertex v_2 );
Vertex FindBase( Vertex u, Vertex *BasePtr );
int FindPathToVertex_s( Vertex x, Edge *SwitchEdge, Vertex *BasePtr, Vertex *Path, int MaxPathLen );
int PullFlow( BN_STRUCT *pBNS, Edge *SwitchEdge, Vertex x, Vertex y, int delta, S_CHAR bReverse, int bChangeFlow );
int FindPathCap( BN_STRUCT* pBNS, Edge *SwitchEdge, Vertex x, Vertex y, int delta );
Vertex MakeBlossom( BN_STRUCT* pBNS, Vertex *ScanQ, int *pQSize, Vertex *Pu, Vertex *Pv, int max_len_Pu_Pv, Edge *SwitchEdge, Vertex *BasePtr, Vertex u, Vertex v, EdgeIndex iuv, Vertex b_u, Vertex b_v, S_CHAR *Tree );
int nGetEndpointInfo_PT_22_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_PT_16_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_PT_06_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_PT_39_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_PT_13_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_PT_18_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_KET( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif );
int is_centerpoint_elem_KET( U_CHAR el_number );
int is_centerpoint_elem_strict( U_CHAR el_number );
int GetChargeType( inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype );
int bCanBeACPoint( inp_ATOM *at, S_CHAR cCharge, S_CHAR cChangeValence, S_CHAR neutral_bonds_valence, S_CHAR neutral_valence, S_CHAR nEndpointValence, S_CHAR *cChargeSubtype );
int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt, int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );
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
            return "1\n" + fmt(c)
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
    path = ROOT / "tests/data/inchi_bnexists.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
