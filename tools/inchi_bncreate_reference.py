"""Generates tests/data/inchi_bncreate.json from the InChI C code: adding one
t-group or c-group to a balanced network by atom type and mask and removing the
last group (`CreateCGroupInBnStruct`, `CreateTGroupInBnStruct`,
`RemoveLastGroupFromBnStruct`, `SetInitCapFlowToCurrent`) in random sequences
on networks built from random molecules. The function texts are cut out of
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
    printf("\nA");
    for (int i = 0; i < n; i++) printf(" %d,%d", at[i].endpoint, at[i].c_point);
    printf("\n");
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, maxadd_a, maxadd_e, maxaltp, tacn, nops;
        scanf("%d %d %d %d %d %d", &n, &maxadd_a, &maxadd_e, &maxaltp, &tacn, &nops);
        inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, vt, chg, nh, ep, cp;
            scanf("%d %d %d %d %d %d %d %d %d", &v, &cbv, &el, &rad, &vt, &chg, &nh, &ep, &cp);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].radical = rad;
            at[i].charge = chg; at[i].num_H = nh; at[i].endpoint = ep; at[i].c_point = cp;
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        int nchg = 0;
        BN_STRUCT *pBNS = AllocateAndInitBnStruct(at, n, maxadd_a, maxadd_e, maxaltp, &nchg);
        pBNS->type_TACN = tacn;
        printf("N %d\n", c);
        dumpstate(pBNS, at, n);
        for (int k = 0; k < nops; k++) {
            int op, ty, mask, chg;
            scanf("%d %d %d %d", &op, &ty, &mask, &chg);
            int r = 0;
            if (op == 0) r = CreateCGroupInBnStruct(at, n, pBNS, ty, mask, chg);
            else if (op == 1) r = CreateTGroupInBnStruct(at, n, pBNS, ty, mask);
            else if (op == 2) r = RemoveLastGroupFromBnStruct(at, n, ty < 0 ? (pBNS->num_vertices > pBNS->num_atoms ? pBNS->num_vertices - 1 : 1000) : ty, pBNS);
            else r = SetInitCapFlowToCurrent(pBNS);
            printf("R %d\n", r);
            dumpstate(pBNS, at, n);
        }
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

A = {n: 1 << i for i, n in enumerate("NP_PLUS NP_PROTON NP_H N_MINUS NP ON COH CSH ZOH OOH ZOOH NOH N_OH CO ZO NO N_O CO_MINUS CS_MINUS ZO_MINUS OO_MINUS ZOO_MINUS NO_MINUS N_O_MINUS O_MINUS OH_PLUS O_PLUS PROTON HAL_ANION HAL_ACID".split())}
MSK_NP = A["NP_PLUS"] | A["NP_PROTON"] | A["NP_H"] | A["N_MINUS"] | A["NP"]
MSK_OS = sum(A[k] for k in "COH CSH ZOH OOH ZOOH NOH N_OH CO ZO NO N_O CO_MINUS CS_MINUS ZO_MINUS OO_MINUS ZOO_MINUS NO_MINUS N_O_MINUS".split())
MSK_H = sum(A[k] for k in "NP_PROTON NP_H COH CSH ZOH OOH ZOOH NOH N_OH".split())
TYPEMASKS = [(0x40 | 0x80 | 0x800, A["NP_PROTON"] | A["OH_PLUS"]), (0x40 | 0x80, MSK_NP), (0x1 | 0x2 | 0x4 | 0x8 | 0x10 | 0x200, MSK_OS),
             (0x1 | 0x2 | 0x4 | 0x8 | 0x10 | 0x4000, A["COH"] | A["CSH"] | A["OOH"] | A["ZOOH"] | A["NOH"] | A["HAL_ACID"]), (0x40 | 0x80 | 0x1 | 0x2 | 0x4 | 0x8 | 0x10 | 0x200, MSK_H),
             (0x100 | 0x400, A["O_MINUS"]), (0x1 | 0x2 | 0x4 | 0x8 | 0x10 | 0x200 | 0x100 | 0x400, MSK_OS | A["O_MINUS"])]


def main():
    rng = random.Random(23)
    cases = []
    for _ in range(2500):
        atoms = gen(rng)
        n = len(atoms)
        for a in atoms:
            a["el"] = rng.choice([6, 7, 8, 8, 7, 16, 15])
            a["chg"] = rng.choice([0, 0, 0, 1, -1])
            if rng.random() < 0.8:
                val = {6: 4, 7: 3, 8: 2, 16: 2, 15: 3}[a["el"]] + a["chg"] * (1 if a["el"] in (7, 15) else -1 if a["el"] in (8, 16) else 0)
                a["nh"] = max(0, val - a["cbv"])
            else:
                a["nh"] = rng.choice([0, 1, 1, 2, 3])
            a["ep"] = 0
            a["cp"] = 0
        ops = []
        for _ in range(rng.randrange(1, 7)):
            op = rng.choice([0, 1, 1, 0, 2, 2, 3])
            ty, mask = rng.choice(TYPEMASKS)
            if rng.random() < 0.3:
                ty |= 1 << rng.randrange(0, 15)
                mask |= 1 << rng.randrange(0, 30)
            if op == 2:
                ty = -1 if rng.random() < 0.85 else rng.randrange(0, n + 4)
            ops.append((op, ty, mask, rng.choice([1, -1, 0])))
        cases.append({"atoms": atoms, "ma": rng.randrange(0, 4), "me": rng.randrange(0, 3), "mp": 1,
                      "tacn": rng.choice([0, 0, 0x200]), "ops": ops})

    def fmt(c):
        t = f"{len(c['atoms'])} {c['ma']} {c['me']} {c['mp']} {c['tacn']} {len(c['ops'])}\n"
        for a in c["atoms"]:
            t += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['vt']} {a['chg']} {a['nh']} {a['ep']} {a['cp']} " + " ".join(f"{x} {t2}" for x, t2 in a["nb"]) + "\n"
        for o in c["ops"]:
            t += " ".join(map(str, o)) + "\n"
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
        "int CreateCGroupInBnStruct( inp_ATOM *at,\n", "int CreateTGroupInBnStruct( inp_ATOM *at,\n", "int RemoveLastGroupFromBnStruct( inp_ATOM *at,\n", "int SetInitCapFlowToCurrent( BN_STRUCT *pBNS )\n{",
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
    endfuncs = [ctype, cut3("int is_centerpoint_elem_KET( U_CHAR el_number )\n{"), cut3("int is_centerpoint_elem_strict( U_CHAR el_number )\n{"),
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
    macros += "\n#define KNOWN_ACIDIC_TYPE (ATT_ACIDIC_CO | ATT_ACIDIC_S | ATT_OO | ATT_ZOO | ATT_NO)\n#define BNS_CHK_ALTP_NO_ALTPATH 0\n#define BNS_CHK_ALTP_SAME_TGROUP 1\n#define BNS_CHK_ALTP_SAME_VERTEX 2\n#define BNS_CHK_ALTP_SET_SUCCESS 4\n" + text[ea:text.index("\n", eb)] + "\n#define MAX_ALT_AATG_ARRAY_LEN 127\n#define AATG_MARK_IN_PATH 1\n"
    alt = re.search(r"typedef struct tagAltPathChanges.*?\} ALT_PATH_CHANGES;", text, re.S).group(0)
    protos = """int CreateCGroupInBnStruct( inp_ATOM *at, int num_atoms, BN_STRUCT *pBNS, int nType, int nMask, int nCharge );
int CreateTGroupInBnStruct( inp_ATOM *at, int num_atoms, BN_STRUCT *pBNS, int nType, int nMask );
int RemoveLastGroupFromBnStruct( inp_ATOM *at, int num_atoms, int tg, BN_STRUCT *pBNS );
int SetInitCapFlowToCurrent( BN_STRUCT *pBNS );
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
    path = ROOT / "tests/data/inchi_bncreate.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
