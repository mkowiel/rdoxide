"""Generates tests/data/inchi_bnvertex.json from the InChI C code: the helpers
of ichi_bns.c that grow a balanced network (`bAddNewVertex`, `AddNewEdge`,
`GetEdgeToGroupVertex`, `GetGroupVertex`, `bAddStCapToAVertex`), applied in
random sequences to networks built from random molecules. The function texts
are cut out of ichi_bns.c. Needs a C compiler and the fetched InChI sources."""

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
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, maxadd_a, maxadd_e, maxaltp, nops;
        scanf("%d %d %d %d %d", &n, &maxadd_a, &maxadd_e, &maxaltp, &nops);
        inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, vt;
            scanf("%d %d %d %d %d", &v, &cbv, &el, &rad, &vt);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].radical = rad; at[i].c_point = vt;
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        int nchg = 0;
        BN_STRUCT *pBNS = AllocateAndInitBnStruct(at, n, maxadd_a, maxadd_e, maxaltp, &nchg);
        for (int i = 0; i < n; i++) pBNS->vert[i].type |= (at[i].c_point & 0x1c) << 1 | (at[i].c_point & 3) << 9;
        printf("N %d\n", c);
        for (int k = 0; k < nops; k++) {
            int op, a, b, x, y, z;
            scanf("%d %d %d %d %d %d", &op, &a, &b, &x, &y, &z);
            int dots = 5;
            if ((op == 0 && a >= pBNS->num_vertices) || ((op == 1 || op == 4) && (a >= pBNS->num_vertices || b >= pBNS->num_vertices))) { printf("skip\n"); continue; }
            if (op == 0) {
                int r = bAddNewVertex(pBNS, a, x, y, z, &dots);
                printf("o0 %d %d\n", r, dots);
            } else if (op == 1) {
                int r = AddNewEdge(pBNS->vert + a, pBNS->vert + b, pBNS, x, y);
                printf("o1 %d\n", r);
            } else if (op == 2) {
                printf("o2 %d\n", GetEdgeToGroupVertex(pBNS, a, x));
            } else if (op == 3) {
                printf("o3 %d\n", GetGroupVertex(pBNS, a, x));
            } else {
                VertexFlow buf[64];
                int r = bAddStCapToAVertex(pBNS, a, b, buf, &dots, x);
                printf("o4 %d %d:", r, dots);
                for (int i = 0; i < r; i++) printf(" %d", buf[i]);
                printf("\n");
            }
            dump(pBNS);
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


def main():
    rng = random.Random(16)
    cases = []
    for _ in range(3000):
        atoms = gen(rng)
        n = len(atoms)
        ma, me = rng.randrange(0, 5), rng.randrange(0, 6)
        nv = n + ma
        ops = []
        for _ in range(rng.randrange(1, 7)):
            op = rng.choice([0, 0, 1, 1, 2, 3, 4, 4])
            if op == 1:
                a = rng.randrange(nv + 1)
                b = rng.randrange(nv + 1)
                if a == b:
                    b = (a + 1) % (nv + 1)
                ops.append((1, a, b, rng.randrange(0, 4), rng.randrange(0, 3), 0))
            else:
                ops.append((op, rng.randrange(nv + 1), rng.randrange(nv + 1), rng.choice([0, 1, 2, 3, 0x40, 4, 8, 2, 0x10]) if op in (2, 3) else rng.randrange(0, 4), rng.randrange(0, 3), rng.randrange(0, 4)))
        cases.append({"atoms": atoms, "ma": ma, "me": me, "mp": 1, "ops": ops})
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{len(c['atoms'])} {c['ma']} {c['me']} {c['mp']} {len(c['ops'])}\n"
        for a in c["atoms"]:
            stdin += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['vt']} " + " ".join(f"{x} {t}" for x, t in a["nb"]) + "\n"
        for o in c["ops"]:
            stdin += " ".join(map(str, o)) + "\n"
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
        "Vertex GetGroupVertex( BN_STRUCT *pBNS, Vertex v1, AT_NUMB type )\n{", "int bAddStCapToAVertex( BN_STRUCT *pBNS,\n"]
    funcs = [cut(n) for n in names]
    taut = (SRC / "ichitaut.c").read_text(errors="replace")
    a = taut.rindex("int is_centerpoint_elem( U_CHAR el_number )\n{")
    funcs.append(taut[a:taut.index("\n}\n", a) + 3])
    a = text.index("#define TREE_NOT_IN_M")
    b = text.index("#define TREE_MARK")
    macros = text[a:text.index("\n", b)]
    alt = re.search(r"typedef struct tagAltPathChanges.*?\} ALT_PATH_CHANGES;", text, re.S).group(0)
    protos = """int bAddNewVertex( BN_STRUCT *pBNS, int nVertDoubleBond, int nCap, int nFlow, int nMaxAdjEdges, int *nDots );
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
        code = HARNESS.replace("@MACROS@", alt + "\n" + macros).replace("@FUNCS@", protos + "\n".join(funcs))
        (Path(d) / "h.c").write_text(code)
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c"), "-lm"], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_bnvertex.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
