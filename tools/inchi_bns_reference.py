"""Generates tests/data/inchi_bns.json from the InChI C code: the balanced
network structure of ichi_bns.c (`AllocateAndInitBnStruct`, `ReInitBnStruct`,
`AllocateAndInitBnData`, `ReInitBnData`, `GetVertexDegree`, `Get2ndEdgeVertex`,
`GetVertexNeighbor`, `GetEdgePointer`, `AugmentEdge`, `rescap`,
`rescap_mark`, `GetPrevVertex`, `bIgnoreVertexNonTACN_atom/group`) on
networks built from random molecules. The function texts are cut out of
ichi_bns.c. Needs a C compiler and the fetched InChI sources."""

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
#define IS_FORBIDDEN(EDGE_FORBIDDEN, PBNS)     (EDGE_FORBIDDEN)
#define TREE_NOT_IN_M  0
#define prim(v) (Vertex)((v)^1)
#define SwitchEdge_Vert1(u) SwitchEdge[u][0]
#define SwitchEdge_Vert2(u) Get2ndEdgeVertex( pBNS, SwitchEdge[u] )
#define SwitchEdge_IEdge(u) SwitchEdge[u][1]
#define CHECK_TACN  1
@FUNCS@
static void dump(BN_STRUCT *p) {
    printf("H %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", p->num_atoms, p->num_added_atoms, p->nMaxAddAtoms, p->num_c_groups, p->num_t_groups,
           p->num_vertices, p->num_bonds, p->num_edges, p->num_added_edges, p->nMaxAddEdges, p->max_vertices, p->max_edges, p->max_iedges,
           p->tot_st_cap, p->tot_st_flow, p->bNotASimplePath);
    for (int i = 0; i < p->num_vertices; i++) {
        BNS_VERTEX *v = &p->vert[i];
        printf("V %d %d %d %d %d %d %d %d %d %d\n", v->type, v->num_adj_edges, v->max_adj_edges, (int)(v->iedge - p->iedge), v->st_edge.cap, v->st_edge.cap0, v->st_edge.flow, v->st_edge.flow0, v->st_edge.pass, 0);
    }
    for (int i = 0; i < p->num_edges; i++) {
        BNS_EDGE *e = &p->edge[i];
        printf("E %d %d %d %d %d %d %d %d %d %d\n", e->neighbor1, e->neighbor12, e->neigh_ord[0], e->neigh_ord[1], e->cap, e->cap0, e->flow, e->flow0, e->pass, e->forbidden);
    }
    printf("I");
    for (int i = 0; i < p->num_atoms; i++) for (int k = 0; k < p->vert[i].num_adj_edges; k++) printf(" %d", p->vert[i].iedge[k]);
    printf("\n");
    if (p->alt_path) {
        printf("P %d %d %d", p->alt_path[iALTP_START_ATOM].number, p->alt_path[iALTP_END_ATOM].number, p->alt_path[iALTP_PATH_LEN].number);
        printf(" %d %d", p->alt_path[iALTP_FLOW].flow[0], p->alt_path[iALTP_FLOW].flow[1]);
        for (int k = 0; k < p->alt_path[iALTP_PATH_LEN].number; k++) printf(" %d,%d", ALTP_THIS_ATOM_NEIGHBOR(p->alt_path, k), ALTP_NEXT_ATOM_NEIGHBOR(p->alt_path, k));
        printf("\n");
    }
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, maxadd_a, maxadd_e, maxaltp, ttacn, tt, tcn;
        scanf("%d %d %d %d %d %d %d", &n, &maxadd_a, &maxadd_e, &maxaltp, &ttacn, &tt, &tcn);
        inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, vt;
            scanf("%d %d %d %d %d", &v, &cbv, &el, &rad, &vt);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].radical = rad; at[i].c_point = vt;
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        int nchg = 0;
        BN_STRUCT *pBNS = AllocateAndInitBnStruct(at, n, maxadd_a, maxadd_e, maxaltp, &nchg);
        printf("N %d\n", nchg);
        for (int i = 0; i < n; i++) { printf("B %d", at[i].valence); for (int j = 0; j < at[i].valence; j++) printf(" %d", at[i].bond_type[j]); printf("\n"); }
        dump(pBNS);
        pBNS->type_TACN = ttacn; pBNS->type_T = tt; pBNS->type_CN = tcn;
        for (int i = 0; i < n; i++) pBNS->vert[i].type |= (at[i].c_point & 0x1c) << 1 | (at[i].c_point & 3) << 9;
        pBNS->alt_path = pBNS->altp[0];
        /* switch-edge table filled from the atoms */
        static Edge SwitchEdge[4096];
        for (int i = 0; i < 4096; i++) { SwitchEdge[i][0] = NO_VERTEX; SwitchEdge[i][1] = 0; }
        int cnt = 0;
        for (int y = 3; y < 2 * pBNS->num_vertices + 2 && pBNS->num_edges > 0; y++) {
            if ((y * 7 + n) % 3 == 0) { SwitchEdge[y][0] = 2 + (y * 5 + 1) % (y - 2); SwitchEdge[y][1] = (y * 3) % pBNS->num_edges; }
        }
        for (int u = 0; u < 2 * pBNS->num_vertices + 2; u++) {
            int deg = GetVertexDegree(pBNS, u);
            printf("D %d %d\n", u, deg);
            for (int k = 0; k < deg; k++) {
                EdgeIndex ie = 12345;
                Vertex w = GetVertexNeighbor(pBNS, u, k, &ie);
                printf("W %d %d %d", k, w, ie);
                if (w != NO_VERTEX) {
                    int rc = rescap(pBNS, u, w, ie);
                    printf(" r%d a%d", rc, bIgnoreVertexNonTACN_atom(pBNS, u, w));
                    if (u > 1 && w > 1) {
                        if (SwitchEdge[u][0] != NO_VERTEX) {
                            printf(" g%d", bIgnoreVertexNonTACN_group(pBNS, u, w, SwitchEdge));
                            EdgeIndex iuv = -77;
                            Vertex pv = GetPrevVertex(pBNS, u, SwitchEdge, &iuv);
                            printf(" p%d,%d", pv, iuv);
                        }
                    }
                    int rm = rescap_mark(pBNS, u, w, ie);
                    printf(" m%d", rm);
                    int delta = (cnt % 4) - 1;
                    int rv = (cnt >> 2) & 1;
                    int bcf = (cnt >> 3) & 1;
                    int ar = AugmentEdge(pBNS, u, w, ie, delta, rv, bcf);
                    printf(" A%d", ar);
                    cnt++;
                }
                printf("\n");
            }
        }
        printf("S %d\n", pBNS->bNotASimplePath);
        dump(pBNS);
        pBNS->alt_path = NULL;
        printf("R %d\n", ReInitBnStruct(pBNS, at, n, 1));
        dump(pBNS);
        BN_DATA *pBD = AllocateAndInitBnData(n + 1);
        pBD->QSize = 2;
        pBD->ScanQ[0] = 2; pBD->ScanQ[1] = 5; pBD->ScanQ[2] = 3;
        pBD->BasePtr[2] = pBD->BasePtr[3] = 7; pBD->Tree[5] = pBD->Tree[4] = 2; pBD->SwitchEdge[2][0] = 9;
        printf("Z %d %d %d %d %d\n", pBD->max_num_vertices, pBD->max_len_Pu_Pv, ReInitBnData(pBD), pBD->QSize, pBD->BasePtr[2]);
        printf("T %d %d %d %d\n", pBD->Tree[4], pBD->Tree[5], pBD->SwitchEdge[2][0], pBD->BasePtr[3]);
        DeAllocateBnData(pBD);
        DeAllocateBnStruct(pBNS);
        free(at);
    }
    return 0;
}
"""

ELS = [6, 6, 6, 7, 7, 8, 8, 16, 15, 17, 1]


def gen(rng):
    n = rng.randrange(1, 8)
    nb = [[] for _ in range(n)]
    for i in range(1, n):
        if rng.random() < 0.9:
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
        cbv = sum(min(t & 15, 3) for _, t in nb[i]) + rng.choice([0, 0, 1])
        atoms.append({"v": v, "cbv": cbv, "el": rng.choice(ELS), "rad": rng.choice([0] * 6 + [2, 1, 3]),
                      "vt": rng.randrange(0, 32), "nb": nb[i]})
    return atoms


def main():
    rng = random.Random(14)
    cases = []
    for _ in range(3000):
        atoms = gen(rng)
        cases.append({"atoms": atoms, "ma": rng.randrange(0, 4), "me": rng.randrange(0, 3), "mp": rng.randrange(1, 4),
                      "tacn": rng.choice([0, 0x200]), "t": rng.choice([0, 4]), "cn": rng.choice([0, 0x110])})
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{len(c['atoms'])} {c['ma']} {c['me']} {c['mp']} {c['tacn']} {c['t']} {c['cn']}\n"
        for a in c["atoms"]:
            stdin += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['vt']} " + " ".join(f"{x} {t}" for x, t in a["nb"]) + "\n"
    text = (SRC / "ichi_bns.c").read_text(errors="replace")

    def cut(start):
        a = text.rindex(start)
        return text[a:text.index("\n}\n", a) + 3]

    names = [
        "BN_STRUCT* AllocateAndInitBnStruct( inp_ATOM *at,\n", "BN_STRUCT* DeAllocateBnStruct( BN_STRUCT *pBNS )\n{",
        "int ReInitBnStructAltPaths( BN_STRUCT *pBNS )\n{", "int ReInitBnStruct( BN_STRUCT *pBNS,\n",
        "void ClearAllBnDataVertices( Vertex *v, Vertex value, int size )\n{", "void ClearAllBnDataEdges( Edge *e, Vertex value, int size )\n{",
        "BN_DATA *DeAllocateBnData( BN_DATA *pBD )\n{", "BN_DATA *AllocateAndInitBnData( int max_num_vertices )\n{",
        "int ReInitBnData( BN_DATA *pBD )\n{", "int GetVertexDegree( BN_STRUCT* pBNS, Vertex v )\n{",
        "Vertex Get2ndEdgeVertex( BN_STRUCT* pBNS, Edge uv )\n{", "Vertex GetVertexNeighbor( BN_STRUCT* pBNS,\n",
        "int GetEdgePointer( BN_STRUCT* pBNS,\n", "int AugmentEdge( BN_STRUCT* pBNS,\n", "int rescap_mark( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv )\n{",
        "Vertex GetPrevVertex( BN_STRUCT* pBNS, Vertex y, Edge *SwitchEdge, EdgeIndex *iuv )\n{",
        "int bIgnoreVertexNonTACN_atom( BN_STRUCT* pBNS, Vertex u, Vertex v )\n{", "int bIgnoreVertexNonTACN_group( BN_STRUCT* pBNS,\n",
        "int rescap( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv )\n{"]
    funcs = [cut(n) for n in names]
    taut = (SRC / "ichitaut.c").read_text(errors="replace")
    a = taut.rindex("int is_centerpoint_elem( U_CHAR el_number )\n{")
    funcs.append(taut[a:taut.index("\n}\n", a) + 3])
    protos = """int rescap( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv );
BN_STRUCT* DeAllocateBnStruct( BN_STRUCT *pBNS );
BN_DATA *DeAllocateBnData( BN_DATA *pBD );
int ReInitBnStructAltPaths( BN_STRUCT *pBNS );
int GetVertexDegree( BN_STRUCT* pBNS, Vertex v );
Vertex GetVertexNeighbor( BN_STRUCT* pBNS, Vertex v, int neigh, EdgeIndex *iedge );
Vertex Get2ndEdgeVertex( BN_STRUCT* pBNS, Edge uv );
int GetEdgePointer( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv, BNS_EDGE **uv, S_CHAR *s_or_t );
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNCS@", protos + "\n".join(funcs)))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_bns.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
