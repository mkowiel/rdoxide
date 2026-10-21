"""Generates tests/data/inchi_bnbonds.json from the InChI C code: how the flow
of a balanced network is transferred to the bonds of the structure
(`SetBondsFromBnStructFlow`, `SetAtomBondType`, `RestoreBnStructFlow`,
`RestoreEdgeFlow`, `SetAtomRadAndChemValFromVertexCapFlow`,
`bSetBondsAfterCheckOneBond`) after a balanced network search on networks built
from random molecules. The function texts are cut out of ichi_bns.c. Needs a C
compiler and the fetched InChI sources."""

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
    for (int i = 0; i < p->num_vertices; i++) {
        BNS_VERTEX *v = &p->vert[i];
        printf("V %d %d %d %d %d %d\n", v->type, v->st_edge.cap, v->st_edge.cap0, v->st_edge.flow, v->st_edge.flow0, v->st_edge.pass);
    }
    for (int i = 0; i < p->num_edges; i++) {
        BNS_EDGE *e = &p->edge[i];
        printf("E %d %d %d %d %d\n", e->cap, e->cap0, e->flow, e->flow0, e->pass);
    }
    printf("P %d %d %d\n", p->num_altp, p->bChangeFlow, p->bNotASimplePath);
    for (int a = 0; a < p->num_altp && a < p->max_altp; a++) {
        BNS_ALT_PATH *ap = p->altp[a];
        printf("A %d %d %d %d %d", ap[iALTP_START_ATOM].number, ap[iALTP_END_ATOM].number, ap[iALTP_PATH_LEN].number, ap[iALTP_FLOW].flow[0], ap[iALTP_FLOW].flow[1]);
        for (int k = 0; k < ap[iALTP_PATH_LEN].number; k++) printf(" %d,%d", ALTP_THIS_ATOM_NEIGHBOR(ap, k), ALTP_NEXT_ATOM_NEIGHBOR(ap, k));
        printf("\n");
    }
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, maxadd_a, maxadd_e, maxaltp, bcf, mode, bcf2, ntest, eidx;
        scanf("%d %d %d %d %d %d %d %d %d", &n, &maxadd_a, &maxadd_e, &maxaltp, &bcf, &mode, &bcf2, &ntest, &eidx);
        inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, vt;
            scanf("%d %d %d %d %d", &v, &cbv, &el, &rad, &vt);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].radical = rad; at[i].c_point = vt;
            at[i].nRingSystem = (i * 7 + c) % 3;
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        int nchg = 0;
        BN_STRUCT *pBNS = AllocateAndInitBnStruct(at, n, maxadd_a, maxadd_e, maxaltp, &nchg);
        printf("N %d\n", c);
        BN_DATA *pBD = AllocateAndInitBnData(pBNS->max_vertices);
        BNS_FLOW_CHANGES fcd[64];
        for (int i = 0; i < 64; i++) fcd[i].iedge = NO_VERTEX;
        if (mode == 2 && pBNS->num_edges > 0) printf("f %d\n", bSetFlowToCheckOneBond(pBNS, eidx % pBNS->num_edges, ntest % 3, fcd));
        int ret = RunBalancedNetworkSearch(pBNS, pBD, bcf);
        printf("X %d\n", ret);
        if (mode == 0) {
            printf("r %d\n", SetBondsFromBnStructFlow(pBNS, at, n, bcf2));
        } else if (mode == 1) {
            printf("r %d\n", RestoreBnStructFlow(pBNS, bcf2));
        } else {
            int r1 = bSetBondsAfterCheckOneBond(pBNS, fcd, ntest - 1, at, n, bcf2);
            int r2 = SetBondsFromBnStructFlow(pBNS, at, n, bcf2);
            int r3 = RestoreBnStructFlow(pBNS, bcf2 & 3);
            bRestoreFlowAfterCheckOneBond(pBNS, fcd);
            printf("r %d %d %d\n", r1, r2, r3);
        }
        for (int i = 0; i < n; i++) {
            printf("T %d %d %d |", at[i].valence, at[i].chem_bonds_valence, at[i].radical);
            for (int j = 0; j < at[i].valence; j++) printf(" %d", at[i].bond_type[j]);
            printf("\n");
        }
        dump(pBNS);
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
            t = rng.choice([1, 1, 2, 3, 4, 1, 2, 0x11, 9, 8, 5, 6, 7, 0x21, 0x31, 0x41, 0x51, 0x14, 0x59])
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
    rng = random.Random(17)
    cases = []
    for _ in range(3000):
        atoms = gen(rng)
        cases.append({"atoms": atoms, "ma": rng.randrange(0, 3), "me": rng.randrange(0, 3), "mp": rng.randrange(1, 5),
                      "bcf": rng.choice([1, 1, 5, 13, 29, 61, 21, 17, 9, 3, 0]), "mode": rng.randrange(3),
                      "bcf2": rng.choice([1, 5, 13, 29, 45, 61, 21, 17, 9, 3, 7, 0, 12, 36, 60, 8, 40]),
                      "nt": rng.randrange(0, 3), "ei": rng.randrange(0, 1000)})
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{len(c['atoms'])} {c['ma']} {c['me']} {c['mp']} {c['bcf']} {c['mode']} {c['bcf2']} {c['nt']} {c['ei']}\n"
        for a in c["atoms"]:
            stdin += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['vt']} " + " ".join(f"{x} {t}" for x, t in a["nb"]) + "\n"
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
        "int RestoreEdgeFlow( BNS_EDGE *edge, int delta, int bChangeFlow )\n{", "int SetAtomBondType( BNS_EDGE *edge,\n",
        "int SetAtomRadAndChemValFromVertexCapFlow( BN_STRUCT *pBNS,\n", "int SetBondsFromBnStructFlow( BN_STRUCT *pBNS,\n",
        "int RestoreBnStructFlow( BN_STRUCT *pBNS, int bChangeFlow )\n{", "int bSetBondsAfterCheckOneBond( BN_STRUCT *pBNS,\n",
        "int bRestoreFlowAfterCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd )\n{", "int bSetFlowToCheckOneBond( BN_STRUCT *pBNS,\n"]
    funcs = [cut(n) for n in names]
    taut = (SRC / "ichitaut.c").read_text(errors="replace")
    a = taut.rindex("int is_centerpoint_elem( U_CHAR el_number )\n{")
    funcs.append(taut[a:taut.index("\n}\n", a) + 3])
    a = text.index("#define TREE_NOT_IN_M")
    b = text.index("#define TREE_MARK")
    macros = text[a:text.index("\n", b)]
    alt = re.search(r"typedef struct tagAltPathChanges.*?\} ALT_PATH_CHANGES;", text, re.S).group(0)
    protos = """int RestoreEdgeFlow( BNS_EDGE *edge, int delta, int bChangeFlow );
int SetAtomBondType( BNS_EDGE *edge, U_CHAR *bond_type12, U_CHAR *bond_type21, int delta, int bChangeFlow );
int SetAtomRadAndChemValFromVertexCapFlow( BN_STRUCT *pBNS, inp_ATOM *atom, int v1 );
int SetBondsFromBnStructFlow( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bChangeFlow0 );
int RestoreBnStructFlow( BN_STRUCT *pBNS, int bChangeFlow );
int bSetBondsAfterCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd, int nTestFlow, inp_ATOM *at, int num_atoms, int bChangeFlow );
int bRestoreFlowAfterCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd );
int bSetFlowToCheckOneBond( BN_STRUCT *pBNS, int iedge, int flow, BNS_FLOW_CHANGES *fcd );
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
        code = HARNESS.replace("@MACROS@", alt + "\n" + macros).replace("@FUNCS@", protos + "\n".join(funcs))
        (Path(d) / "h.c").write_text(code)
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c"), "-lm"], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_bnbonds.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
