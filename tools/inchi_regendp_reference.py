"""Generates tests/data/inchi_regendp.json from the InChI C code:
`RegisterEndPoints` (creating, extending and merging t-groups and renumbering
the endpoints, with or without rebuilding the balanced network) and
`SetTautomericBonds` of ichitaut.c on random molecules with random t-groups
and endpoint lists. The function texts are cut out of the InChI sources.
Needs a C compiler and the fetched InChI sources."""

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
#define FIX_CPOINT_BOND_CAP 1
#define RESET_EDGE_FORBIDDEN_MASK 1
#define MAX_STACK_ARRAY_LEN 127
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
static void dumpstate(BN_STRUCT *p, inp_ATOM *at, int n, T_GROUP_INFO *t, C_GROUP_INFO *cg) {
    printf("C %d %d %d %d %d %d\n", p->num_vertices, p->num_edges, p->num_t_groups, p->num_c_groups, p->num_added_atoms, p->num_added_edges);
    printf("V");
    for (int i = 0; i < p->num_vertices; i++) {
        BNS_VERTEX *v = &p->vert[i];
        printf(" %d:%d:%d:%d:%d:%d:%d:%d|", v->type, v->st_edge.cap, v->st_edge.cap0, v->st_edge.flow, v->st_edge.flow0, v->num_adj_edges, v->max_adj_edges, (int)(v->iedge - p->iedge));
        for (int k = 0; k < v->num_adj_edges; k++) printf("%d.", v->iedge[k]);
    }
    printf("\nE");
    for (int i = 0; i < p->num_edges; i++) {
        BNS_EDGE *e = &p->edge[i];
        printf(" %d:%d:%d:%d:%d:%d:%d:%d:%d:%d", e->cap, e->cap0, e->flow, e->flow0, e->pass, e->forbidden, e->neighbor1, e->neighbor12, e->neigh_ord[0], e->neigh_ord[1]);
    }
    printf("\nG");
    for (int i = 0; i < t->num_t_groups; i++) printf(" %d", t->t_group[i].nGroupNumber);
    printf(" |");
    for (int i = 0; i < cg->num_c_groups; i++) printf(" %d", cg->c_group[i].nGroupNumber);
    printf("\nA");
    for (int i = 0; i < n; i++) printf(" %d,%d", at[i].endpoint, at[i].c_point);
    printf("\n");
}
#define NNUM 5
#define NDA 6
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, maxadd_a, maxadd_e, maxaltp, tfl, tflags, ntg, ncg, usebns, maxt, nep, nbp;
        scanf("%d %d %d %d %d %d %d %d %d %d %d %d", &n, &maxadd_a, &maxadd_e, &maxaltp, &tfl, &tflags, &ntg, &ncg, &usebns, &maxt, &nep, &nbp);
        inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, vt, chg, nh, ep, cp;
            scanf("%d %d %d %d %d %d %d %d %d", &v, &cbv, &el, &rad, &vt, &chg, &nh, &ep, &cp);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].radical = rad;
            at[i].charge = chg; at[i].num_H = nh; at[i].endpoint = ep; at[i].c_point = cp;
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        T_GROUP_INFO tgi;
        memset(&tgi, 0, sizeof(tgi));
        tgi.t_group = calloc(maxt + 2, sizeof(T_GROUP));
        tgi.num_t_groups = ntg;
        tgi.max_num_t_groups = maxt;
        tgi.bTautFlags = tflags;
        for (int i = 0; i < ntg; i++) {
            int g, ne;
            scanf("%d %d", &g, &ne);
            tgi.t_group[i].nGroupNumber = g; tgi.t_group[i].nNumEndpoints = ne;
            for (int k = 0; k < NNUM; k++) { int x; scanf("%d", &x); tgi.t_group[i].num[k] = x; }
            for (int k = 0; k < NDA; k++) { int x; scanf("%d", &x); tgi.t_group[i].num_DA[k] = x; }
        }
        C_GROUP_INFO cgi;
        memset(&cgi, 0, sizeof(cgi));
        cgi.c_group = calloc(ncg + 1, sizeof(C_GROUP));
        cgi.num_c_groups = ncg;
        for (int i = 0; i < ncg; i++) { int g, ne; scanf("%d %d", &g, &ne); cgi.c_group[i].nGroupNumber = g; cgi.c_group[i].num_CPoints = ne; }
        T_ENDPOINT *ep = calloc(nep + 1, sizeof(T_ENDPOINT));
        for (int i = 0; i < nep; i++) {
            int a1, g, e;
            scanf("%d %d %d", &a1, &g, &e);
            ep[i].nAtomNumber = a1; ep[i].nGroupNumber = g; ep[i].nEquNumber = e;
            for (int k = 0; k < NNUM; k++) { int x; scanf("%d", &x); ep[i].num[k] = x; }
            for (int k = 0; k < NDA; k++) { int x; scanf("%d", &x); ep[i].num_DA[k] = x; }
        }
        T_BONDPOS *bp = calloc(nbp + 1, sizeof(T_BONDPOS));
        for (int i = 0; i < nbp; i++) { int a1, b1; scanf("%d %d", &a1, &b1); bp[i].nAtomNumber = a1; bp[i].neighbor_index = b1; }
        int nchg = 0;
        BN_STRUCT *pBNS = AllocateAndInitBnStruct(at, n, maxadd_a, maxadd_e, maxaltp, &nchg);
        INCHI_MODE tf2 = tfl;
        pBNS->pbTautFlags = &tf2;
        printf("N %d\n", c);
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        for (int rep = 0; rep < 2; rep++) {
            int ret = RegisterEndPoints(&cg, &tgi, ep, nep, at, n, ncg ? &cgi : NULL, usebns ? pBNS : NULL);
            printf("R %d %d\n", ret, tgi.num_t_groups);
            for (int i = 0; i < tgi.num_t_groups; i++) {
                printf("T %d %d |", tgi.t_group[i].nGroupNumber, tgi.t_group[i].nNumEndpoints);
                for (int k = 0; k < NNUM; k++) printf(" %d", tgi.t_group[i].num[k]);
                printf(" |");
                for (int k = 0; k < NDA; k++) printf(" %d", tgi.t_group[i].num_DA[k]);
                printf("\n");
            }
            printf("I");
            if (tgi.tGroupNumber) for (int i = 0; i < 2 * maxt; i++) printf(" %d", tgi.tGroupNumber[i]);
            printf("\nQ");
            for (int i = 0; i < nep; i++) printf(" %d,%d,%d", ep[i].nAtomNumber, ep[i].nGroupNumber, ep[i].nEquNumber);
            printf("\n");
            dumpstate(pBNS, at, n, &tgi, &cgi);
            /* the endpoint list as processed by the first call is the input of the second */
            for (int i = 0; i < nep; i++) ep[i].nGroupNumber = at[ep[i].nAtomNumber].endpoint;
        }
        printf("S %d\n", SetTautomericBonds(at, nbp, bp));
        for (int i = 0; i < n; i++) {
            printf("B");
            for (int j = 0; j < at[i].valence; j++) printf(" %d", at[i].bond_type[j]);
            printf("\n");
        }
        DeAllocateBnStruct(pBNS);
        free(tgi.t_group);
        free(cgi.c_group);
        free(ep);
        free(bp);
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
FLAGS = [0x08, 0x08 | 0x80000, 0x08 | 0x80000 | 0x00800000 | 0x01000000 | 0x02000000 | 0x04000000 | 0x08000000 | 0x10000000, 0]


def main():
    rng = random.Random(23)
    cases = []
    for _ in range(3000):
        atoms = gen(rng)
        n = len(atoms)
        ntg, ncg = rng.randrange(0, 4), rng.randrange(0, 3)
        for a in atoms:
            a["el"] = rng.choice(ELS2)
            a["chg"] = rng.choice([0, 0, 0, 1, -1])
            a["nh"] = rng.choice([0, 1, 1, 2])
            a["ep"] = rng.randrange(1, ntg + 1) if ntg and rng.random() < 0.5 else 0
            a["cp"] = rng.randrange(1, ncg + 1) if ncg and rng.random() < 0.5 else 0
        tgl = []
        for g in range(1, ntg + 1):
            ne = sum(1 for a in atoms if a["ep"] == g)
            tgl.append((g, ne, [rng.randrange(0, 4) for _ in range(5)], [rng.randrange(0, 3) for _ in range(6)]))
        cgl = [(g, sum(1 for a in atoms if a["cp"] == g)) for g in range(1, ncg + 1)]
        maxt = ntg + rng.choice([0, 1, 1, 2, 3])
        nep = rng.randrange(1, min(n, 5) + 1)
        atn = rng.sample(range(n), nep)
        mode = rng.randrange(3)
        eps = []
        for at_no in atn:
            g = atoms[at_no]["ep"]
            if mode == 0:
                e = 0
            elif mode == 1:
                e = rng.choice([g, g, 0 if False else rng.randrange(1, ntg + 2), n + 1 + rng.randrange(0, 3)]) if g or True else g
                e = e or n + 1
            else:
                e = rng.choice([g, rng.randrange(1, ntg + 2), n + 1]) or 1
            eps.append((at_no, g, e, [rng.randrange(0, 2) for _ in range(5)], [rng.randrange(0, 2) for _ in range(6)]))
        bps = []
        for _ in range(rng.randrange(0, 4)):
            i = rng.randrange(n)
            if atoms[i]["v"]:
                bps.append((i, rng.randrange(atoms[i]["v"])))
        cases.append({"atoms": atoms, "ma": rng.randrange(0, 3), "me": rng.randrange(0, 3), "mp": 1,
                      "tfl": rng.choice(FLAGS), "tf": rng.choice(FLAGS), "tg": tgl, "cg": cgl, "usebns": rng.randrange(2),
                      "maxt": maxt, "ep": eps, "bp": bps})

    def fmt(c):
        t = f"{len(c['atoms'])} {c['ma']} {c['me']} {c['mp']} {c['tfl']} {c['tf']} {len(c['tg'])} {len(c['cg'])} {c['usebns']} {c['maxt']} {len(c['ep'])} {len(c['bp'])}\n"
        for a in c["atoms"]:
            t += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['vt']} {a['chg']} {a['nh']} {a['ep']} {a['cp']} " + " ".join(f"{x} {t2}" for x, t2 in a["nb"]) + "\n"
        for g, ne, nums, da in c["tg"]:
            t += f"{g} {ne} " + " ".join(map(str, nums + da)) + "\n"
        for g, ne in c["cg"]:
            t += f"{g} {ne}\n"
        for at_no, g, e, nums, da in c["ep"]:
            t += f"{at_no} {g} {e} " + " ".join(map(str, nums + da)) + "\n"
        for i, j in c["bp"]:
            t += f"{i} {j}\n"
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
        "int ReInitBnStruct( BN_STRUCT *pBNS,\n", "int ReInitBnStructAddGroups( CANON_GLOBALS *pCG,\n", "int AddTGroups2BnStruct( CANON_GLOBALS *pCG,\n",
        "int AddCGroups2BnStruct( CANON_GLOBALS *pCG,\n", "int CompTGroupNumber( const void *tg1, const void *tg2, void *p )\n{", "int CompCGroupNumber( const void *cg1, const void *cg2, void *p )\n{"]
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
    endfuncs = [ctype, cut3("int RegisterEndPoints( CANON_GLOBALS *pCG,\n"), cut3("int SetTautomericBonds( inp_ATOM *at,\n"), cut3("int is_centerpoint_elem_KET( U_CHAR el_number )\n{"), cut3("int is_centerpoint_elem_strict( U_CHAR el_number )\n{"),
        cut3("int bCanBeACPoint( inp_ATOM *at,\n"), cut3("int GetChargeType( inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype )\n{")]
    for v in ["", "_PT_22_00", "_PT_16_00", "_PT_06_00", "_PT_39_00", "_PT_13_00", "_PT_18_00"]:
        sig = "int nGetEndpointInfo%s(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)" % v if v else "int nGetEndpointInfo( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif )"
        endfuncs.append(cut3(sig))
    endfuncs.append(cut3("int nGetEndpointInfo_KET( inp_ATOM *atom,\n"))
    st = (SRC / "ichister.c").read_text(errors="replace")
    a2 = st.rindex("int get_opposite_sb_atom(")
    opp = st[a2:st.index("\n}\n", a2) + 3]
    isrt = (SRC / "ichisort.c").read_text(errors="replace")
    a3 = isrt.rindex("int insertions_sort( void *pCG,")
    swp = isrt.rindex("void inchi_swap(")
    funcs = endfuncs[:1] + [isrt[swp:isrt.index("\n}\n", swp) + 3], isrt[a3:isrt.index("\n}\n", a3) + 3]] + [cut(n) for n in names] + [opp] + nodefuncs + endfuncs[1:]
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
    protos = """int RegisterEndPoints( CANON_GLOBALS *pCG, T_GROUP_INFO *t_group_info, T_ENDPOINT *EndPoint, int nNumEndPoints, inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi, struct BalancedNetworkStructure *pBNS );
int SetTautomericBonds( inp_ATOM *at, int nNumBondPos, T_BONDPOS *BondPos );
int ReInitBnStruct( BN_STRUCT *pBNS, inp_ATOM *at, int num_at, int bRemoveGroupsFromAtoms );
int ReInitBnStructAddGroups( CANON_GLOBALS *pCG, BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, T_GROUP_INFO *tgi, C_GROUP_INFO *cgi );
int AddTGroups2BnStruct( CANON_GLOBALS *pCG, BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, T_GROUP_INFO *tgi );
int AddCGroups2BnStruct( CANON_GLOBALS *pCG, BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi );
int insertions_sort( void *pCG, void *base, size_t num, size_t width, int( *compare )( const void *, const void *, void * ) );
void inchi_swap( char *a, char *b, size_t width );
int CompTGroupNumber( const void *tg1, const void *tg2, void *p );
int CompCGroupNumber( const void *cg1, const void *cg2, void *p );
int bExistsAltPath( CANON_GLOBALS *pCG, BN_STRUCT *pBNS, BN_DATA *pBD, BN_AATG *pAATG, inp_ATOM *at, int num_atoms, int nVertDoubleBond, int nVertSingleBond, int path_type );
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
    path = ROOT / "tests/data/inchi_regendp.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
