"""Generates tests/data/inchi_bnblocks.json from the InChI C code: the blocks
of the graph of alternating bonds and the non-stereo bond marking
(`ReInitBnStructForAltBns`, `MarkRingSystemsAltBns`, `MarkNonStereoAltBns`)
on networks built from random molecules with alternating and marked bonds.
The function texts are cut out of ichi_bns.c. Needs a C compiler and the
fetched InChI sources."""

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
static int ob[256][20];
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
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; ob[i][j] = b; }
        }
        int nchg = 0;
        BN_STRUCT *pBNS = AllocateAndInitBnStruct(at, n, maxadd_a, maxadd_e, maxaltp, &nchg);
        printf("N %d\n", c);
        BN_DATA *pBD = AllocateAndInitBnData(pBNS->max_vertices);
        /* the allocation changes alternating bonds: restore the input bonds */
        for (int i = 0; i < n; i++) for (int j = 0; j < at[i].valence; j++) at[i].bond_type[j] = ob[i][j];
        for (int i = 0; i < n; i++) at[i].endpoint = at[i].c_point;
        int bu = bcf2 & 1;
        int r0 = ReInitBnStructForAltBns(pBNS, at, n, bu);
        printf("a %d\n", r0);
        for (int i = 0; i < pBNS->num_edges; i++) printf("e %d %d %d %d\n", pBNS->edge[i].nBondTypeInpAltBns, pBNS->edge[i].nBlockNumberAltBns, pBNS->edge[i].nNumAtInBlockAltBns, pBNS->edge[i].forbidden);
        int r1 = MarkRingSystemsAltBns(pBNS, bu);
        printf("b %d\n", r1);
        for (int i = 0; i < pBNS->num_edges; i++) printf("e %d %d %d\n", pBNS->edge[i].nBondTypeInpAltBns, pBNS->edge[i].nBlockNumberAltBns, pBNS->edge[i].nNumAtInBlockAltBns);
        int r2 = MarkNonStereoAltBns(pBNS, at, n, bu);
        printf("c %d\n", r2);
        for (int i = 0; i < n; i++) {
            printf("s");
            for (int j = 0; j < at[i].valence; j++) printf(" %d.%d", at[i].bond_type[j], at[i].bond_stereo[j]);
            printf("\n");
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


def gen_real(rng):
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


MARKS = [0x10, 0x20, 0x30, 0x40, 0x50]


def gen(rng):
    atoms = gen_real(rng)
    for _ring in range(rng.choice([0, 1, 1, 2])):
        # a Kekule ring (benzene, pyridine, ...) attached to the first atom
        m = len(atoms)
        size = rng.choice([5, 6, 6, 6])
        for k in range(size):
            el, chg, nh = (7, 0, 0) if rng.random() < 0.2 else (6, 0, 1)
            atoms.append({"v": 0, "cbv": 0, "el": el, "rad": 0, "chg": chg, "nh": nh, "nb": []})
        ring_alt = rng.random() < 0.6
        for k in range(size):
            j = (k + 1) % size
            t = 2 if k % 2 == 0 and not (size == 5 and k == size - 1) else 1
            if ring_alt:
                t = 4
            atoms[m + k]["nb"].append((m + j, t))
            atoms[m + j]["nb"].append((m + k, t))
        atoms[m]["nb"].append((0, 1))
        atoms[0]["nb"].append((m, 1))
        for k in range(m, m + size):
            a = atoms[k]
            a["v"] = len(a["nb"])
            a["cbv"] = sum(t for _, t in a["nb"])
            if a["el"] == 6:
                a["nh"] = max(0, 4 - a["cbv"])
            else:
                a["nh"] = max(0, 3 - a["cbv"])
        atoms[0]["v"] = len(atoms[0]["nb"])
        atoms[0]["cbv"] = sum(t for _, t in atoms[0]["nb"])
        if atoms[0]["el"] == 6:
            atoms[0]["nh"] = max(0, 4 - atoms[0]["cbv"])
    for i, a in enumerate(atoms):
        a["vt"] = 1 if rng.random() < 0.08 else 0
        if rng.random() < 0.3:
            a["rad"] = rng.choice([2, 1, 3])
    if rng.random() < 0.7:
        for i, a in enumerate(atoms):
            for k, (j, t) in enumerate(a["nb"]):
                if j > i and rng.random() < 0.3:
                    nt = rng.choice([t | rng.choice(MARKS), 4, 9, 4 | rng.choice(MARKS)])
                    a["nb"][k] = (j, nt)
                    atoms[j]["nb"] = [(x, nt if x == i else y) for x, y in atoms[j]["nb"]]
    return atoms


def main():
    rng = random.Random(17)
    cases = []
    for _ in range(1500):
        atoms = gen(rng)
        cases.append({"atoms": atoms, "ma": rng.randrange(0, 3), "me": rng.randrange(0, 3), "mp": rng.randrange(1, 5),
                      "bcf": 1, "mode": rng.randrange(2),
                      "bcf2": rng.choice([1, 5, 13, 29, 45, 61, 21, 17, 9, 3, 7, 13, 29]),
                      "nt": rng.choice([0, 0, 0, 1, 2, 3]), "ei": rng.randrange(0, 1000)})
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
        "int bRestoreFlowAfterCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd )\n{", "int bSetFlowToCheckOneBond( BN_STRUCT *pBNS,\n", "int ReInitBnStruct( BN_STRUCT *pBNS,\n", "int ReInitBnStructForAltBns( BN_STRUCT *pBNS,\n", "int MarkRingSystemsAltBns( BN_STRUCT* pBNS, int bUnknAltAsNoStereo )\n{", "int MarkNonStereoAltBns( BN_STRUCT *pBNS,\n",
        "int nMaxFlow2Check( BN_STRUCT *pBNS, int iedge )\n{", "int nCurFlow2Check( BN_STRUCT *pBNS, int iedge )\n{", "int nMinFlow2Check( BN_STRUCT *pBNS, int iedge )\n{",
        "int bNeedToTestTheFlow( int bond_type,\n", "int nBondsValenceInpAt( const inp_ATOM *at,\n", "int BnsAdjustFlowBondsRad( BN_STRUCT *pBNS,\n", "int BnsTestAndMarkAltBonds( BN_STRUCT *pBNS,\n"]
    funcs = [cut(n) for n in names]
    taut = (SRC / "ichitaut.c").read_text(errors="replace")
    a = taut.rindex("int is_centerpoint_elem( U_CHAR el_number )\n{")
    funcs.append(taut[a:taut.index("\n}\n", a) + 3])
    a = text.index("#define TREE_NOT_IN_M")
    b = text.index("#define TREE_MARK")
    macros = text[a:text.index("\n", b)]
    ba = text.index("#define BT_ALTERN_BOND")
    bb = text.index("#define valenceAltBns")
    macros += "\n#define BNS_MARK_ONLY_BLOCKS 1\n#define RESET_EDGE_FORBIDDEN_MASK 1\n" + text[ba:text.index("\n", bb)] + "\n"
    alt = re.search(r"typedef struct tagAltPathChanges.*?\} ALT_PATH_CHANGES;", text, re.S).group(0)
    protos = """int ReInitBnStruct( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bRemoveGroupsFromAtoms );
int ReInitBnStructForAltBns( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bUnknAltAsNoStereo );
int MarkRingSystemsAltBns( BN_STRUCT* pBNS, int bUnknAltAsNoStereo );
int MarkNonStereoAltBns( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bUnknAltAsNoStereo );
int nMaxFlow2Check( BN_STRUCT *pBNS, int iedge );
int nCurFlow2Check( BN_STRUCT *pBNS, int iedge );
int nMinFlow2Check( BN_STRUCT *pBNS, int iedge );
int bNeedToTestTheFlow( int bond_type, int nTestFlow, int bTestForNonStereoBond );
int nBondsValenceInpAt( const inp_ATOM *at, int *nNumAltBonds, int *nNumWrongBonds );
int BnsAdjustFlowBondsRad( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms );
int BnsTestAndMarkAltBonds( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms, BNS_FLOW_CHANGES *fcd, int bChangeFlow, int nBondTypeToTest );
int RestoreEdgeFlow( BNS_EDGE *edge, int delta, int bChangeFlow );
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
    path = ROOT / "tests/data/inchi_bnblocks.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
