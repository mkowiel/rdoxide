"""Generates tests/data/inchi_isoh.json from the InChI C code:
`MakeIsotopicHGroup` of ichitaut.c (acceptor and donor pairs of "salt"
candidates connected by alternating paths, merged into one t-group) on
networks built from random acids, bases and salts with t-groups and charge
groups. The function texts are cut out of the InChI sources. The field "tf"
packs the network taut flags (low 16 bits), the number of candidates (bits
24-31) and the t-group flag choices (bits 16-23). Needs a C compiler and the
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
#define FIX_CPOINT_BOND_CAP 1
#define RESET_EDGE_FORBIDDEN_MASK 1
#define MAX_STACK_ARRAY_LEN 127
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
static const char *elsym(int el) {
    switch (el) { case 1: return "H"; case 6: return "C"; case 7: return "N"; case 8: return "O"; case 15: return "P"; case 16: return "S"; case 17: return "Cl"; case 35: return "Br"; default: return "Xx"; }
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
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; strcpy(at[i].elname, elsym(el)); at[i].radical = rad; at[i].c_point = vt;
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
        INCHI_MODE tfl = tflags & 0xffff;
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
        BnsAdjustFlowBondsRad(pBNS, pBD, at, n);
        dumpstate(pBNS, at, n);
        C_GROUP_INFO cgi;
        memset(&cgi, 0, sizeof(cgi));
        int maxcp = 0, maxep = 0;
        for (int i = 0; i < n; i++) { if (at[i].c_point > maxcp) maxcp = at[i].c_point; if (at[i].endpoint > maxep) maxep = at[i].endpoint; }
        cgi.num_c_groups = maxcp;
        cgi.max_num_c_groups = maxcp;
        cgi.c_group = calloc(maxcp + 1, sizeof(C_GROUP));
        for (int g = 1; g <= maxcp; g++) {
            cgi.c_group[g - 1].nGroupNumber = g;
            for (int i = 0; i < n; i++) if (at[i].c_point == g) { cgi.c_group[g - 1].num_CPoints++; cgi.c_group[g - 1].num[0] += (at[i].charge == 1); }
        }
        int fch = (tflags >> 16) & 0xff;
        T_GROUP_INFO tgi;
        memset(&tgi, 0, sizeof(tgi));
        tgi.max_num_t_groups = maxep + 3;
        tgi.t_group = calloc(maxep + 5, sizeof(T_GROUP));
        tgi.num_t_groups = maxep;
        for (int g = 1; g <= maxep; g++) {
            tgi.t_group[g - 1].nGroupNumber = g;
            for (int i = 0; i < n; i++) if (at[i].endpoint == g) { tgi.t_group[g - 1].nNumEndpoints++; tgi.t_group[g - 1].num[0] += (at[i].num_H > 0); tgi.t_group[g - 1].num[1] += (at[i].charge == -1); }
        }
        tgi.tGroupNumber = calloc(maxep + 5, sizeof(AT_NUMB));
        for (int g = 1; g <= maxep; g++) tgi.tGroupNumber[g] = (fch >= 60 && g == 1) ? 0 : g;
        tgi.bTautFlags = (fch & 1) ? 0x20 : 0;
        tgi.bTautFlagsDone = (fch & 2) ? 0x2000 : 0;
        { static const int nf[4] = {0, 0x2, 0x10, 0x20}; tgi.tni.bNormalizationFlags = nf[(fch >> 2) & 3]; }
        S_GROUP_INFO sgi;
        memset(&sgi, 0, sizeof(sgi));
        sgi.max_num_candidates = (tflags >> 24) & 0xff;
        sgi.s_candidate = calloc(sgi.max_num_candidates + 2, sizeof(S_CANDIDATE));
        for (int k = 0; k < ncalls; k++) {
            int a1, g, nbi;
            scanf("%d %d %d", &a1, &g, &nbi);
            int ret = MakeIsotopicHGroup(at, n, &sgi, &tgi);
            printf("R %d %d %d\n", ret, sgi.num_candidates, tgi.num_t_groups);
            for (int i = 0; i < sgi.max_num_candidates + 1; i++) printf("K %d %d %d %d\n", sgi.s_candidate[i].atnumber, sgi.s_candidate[i].type, sgi.s_candidate[i].subtype, sgi.s_candidate[i].endpoint);
            for (int i = 0; i < tgi.num_t_groups; i++) {
                printf("T %d %d |", tgi.t_group[i].nGroupNumber, tgi.t_group[i].nNumEndpoints);
                for (int q = 0; q < 5; q++) printf(" %d", tgi.t_group[i].num[q]);
                printf(" |");
                for (int q = 0; q < 6; q++) printf(" %d", tgi.t_group[i].num_DA[q]);
                printf("\n");
            }
            printf("F %d\n", (int)tgi.bTautFlagsDone);
            printf("I %d %d %d %d |", tgi.num_iso_H[0], tgi.num_iso_H[1], tgi.num_iso_H[2], tgi.nNumIsotopicEndpoints);
            if (tgi.nIsotopicEndpointAtomNumber) for (int i = 0; i <= tgi.nIsotopicEndpointAtomNumber[0]; i++) printf(" %d", tgi.nIsotopicEndpointAtomNumber[i]);
            printf(" |");
            for (int i = 0; i < n; i++) printf(" %d", at[i].cFlags);
            printf("\n");
            printf("P");
            for (int i = 0; i < n; i++) printf(" %d,%d", at[i].endpoint, at[i].c_point);
            printf("\n");
            dumpstate(pBNS, at, n);
        }
        free(cgi.c_group);
        free(tgi.t_group);
        free(sgi.s_candidate);
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
NGROUPS = [
    ([(6, 0, 0), (7, 0, 1), (7, 0, 2)], [(0, 1, 2), (0, 2, 1)]),            # amidine
    ([(6, 0, 0), (7, 1, 2), (7, 0, 2)], [(0, 1, 2), (0, 2, 1)]),            # amidinium
    ([(6, 0, 0), (7, 1, 2), (7, 0, 2), (7, 0, 2)], [(0, 1, 2), (0, 2, 1), (0, 3, 1)]),
    ([(6, 0, 0), (7, 0, 1), (7, 0, 1)], [(0, 1, 2), (0, 2, 1)]),
    ([(6, 0, 0), (7, 1, 1), (7, 0, 1)], [(0, 1, 2), (0, 2, 1)]),
    ([(7, 1, 0), (8, -1, 0)], [(0, 1, 1)]),
    ([(7, 1, 0), (6, 0, 0), (6, 0, 0)], [(0, 1, 2), (0, 2, 1)]),
    ([(7, 0, 0), (6, 0, 1)], [(0, 1, 2)]),
    ([(7, 1, 1), (6, 0, 1)], [(0, 1, 2)]),
    ([(7, -1, 0), (6, 0, 0)], [(0, 1, 1)]),
]
SGROUPS = []
for h2, h7 in [(1, 0), (0, 1), (1, 1), (0, 0)]:
    chg2 = -1 if not h2 else 0
    chg7 = -1 if not h7 else 0
    SGROUPS.append(([(6, 0, 0), (8, 0, 0), (8, chg2, h2), (6, 0, 0), (6, 0, 0), (6, 0, 0), (8, 0, 0), (8, chg7, h7)],
                    [(0, 1, 2), (0, 2, 1), (0, 3, 1), (3, 4, 2), (4, 5, 1), (5, 6, 2), (5, 7, 1)]))
    SGROUPS.append(([(6, 0, 0), (8, 0, 0), (8, chg2, h2), (6, 0, 0), (7, 0, 1), (6, 0, 0), (8, 0, 0), (8, chg7, h7)],
                    [(0, 1, 2), (0, 2, 1), (0, 3, 1), (3, 4, 2), (4, 5, 1), (5, 6, 2), (5, 7, 1)]))
    SGROUPS.append(([(6, 0, 0), (8, 0, 0), (8, chg2, h2), (6, 0, 0), (7, 0, 0), (6, 0, 0), (7, 0, 1), (8, chg7, h7)],
                    [(0, 1, 2), (0, 2, 1), (0, 3, 1), (3, 4, 2), (4, 5, 1), (5, 6, 2), (5, 7, 1)]))
for h1, h7 in [(1, 1), (0, 1), (1, 0), (0, 0)]:
    SGROUPS.append(([(6, 0, 0), (8, -1 if not h1 else 0, h1), (6, 0, 1), (6, 0, 1), (6, 0, 1), (8, -1 if not h7 else 0, h7)],
                    [(0, 1, 1), (0, 2, 2), (2, 3, 1), (3, 4, 2), (4, 5, 1)]))
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
        els, bonds = rng.choice(GROUPS + 3 * SGROUPS)
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


MODES = []


def main():
    rng = random.Random(21)
    cases = []
    for _ in range(1200):
        atoms = gen(rng)
        n = len(atoms)
        nh = rng.randrange(0, 3)
        ma, me = rng.randrange(6, 12), rng.randrange(20, 40)
        for i, a in enumerate(atoms):
            a["iso"] = [rng.choice([0, 0, 0, 1]) for _ in range(3)]
            a["orig"] = i + 1
            a["vt"] = 0
            a["ep"] = 0
            a["att"] = 0
            a["ring"] = rng.choice([0, 1, 5, 6])
            a["sbp"] = [0, 0, 0]
            a["sbo"] = [0, 0, 0]
            a["sno"] = [0, 0, 0]
            a["sng"] = [0, 0, 0]
        hs = [(rng.randrange(n), 100 + j) for j in range(nh)]
        ntg = rng.randrange(0, 4)
        for a in atoms:
            a["ep"] = rng.randrange(1, ntg + 1) if ntg and rng.random() < 0.08 else 0
            a["vt"] = rng.randrange(1, 3) if rng.random() < 0.25 else 0
        calls = [(0, 0, 0)] * rng.randrange(1, 3)
        fch = rng.randrange(0, 64)
        mcand = rng.choice([2, 3, n, n + 2, 2 * n])
        cases.append({"atoms": atoms, "hs": hs, "ma": ma, "me": me, "mp": rng.randrange(4, 9), "tf": rng.choice([0x10, 0x10, 0x1010]) | (fch << 16) | (min(mcand, 255) << 24),
                      "ops": [], "calls": calls})
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
        "int bExistsAltPath( CANON_GLOBALS *pCG,\n", "int bExistsAnyAltPath( CANON_GLOBALS *pCG,\n",
        "int CreateCGroupInBnStruct( inp_ATOM *at,\n", "int CreateTGroupInBnStruct( inp_ATOM *at,\n", "int RemoveLastGroupFromBnStruct( inp_ATOM *at,\n",
        "int SimpleRemoveHplusNPO( inp_ATOM *at,\n", "int bIsAtomTypeHard( inp_ATOM *at,\n", "int SimpleRemoveAcidicProtons( inp_ATOM *at,\n",
        "int SimpleAddAcidicProtons( inp_ATOM *at,\n", "int HardRemoveAcidicProtons( CANON_GLOBALS *pCG,\n", "int HardAddAcidicProtons( CANON_GLOBALS *pCG,\n",
        "int HardRemoveHplusNP( CANON_GLOBALS *pCG,\n", "int BnsAdjustFlowBondsRad( BN_STRUCT *pBNS,\n", "int ReInitBnStruct( BN_STRUCT *pBNS,\n", "int AddTGroups2BnStruct( CANON_GLOBALS *pCG,\n", "int AddCGroups2BnStruct( CANON_GLOBALS *pCG,\n", "int CompTGroupNumber( const void *tg1, const void *tg2, void *p )\n{", "int CompCGroupNumber( const void *cg1, const void *cg2, void *p )\n{", "int bIsAtomTypeHard( inp_ATOM *at,\n", "int bIsHDonorAccAtomType( inp_ATOM *at, int endpoint, int *cSubType )\n{", "int bIsNegAtomType( inp_ATOM *at, int endpoint, int *cSubType )\n{", "int bIsHardRemHCandidate( inp_ATOM *at, int i, int *cSubType )\n{", "int nBondsValenceInpAt( const inp_ATOM *at,\n"]
    ican = (SRC / "ichican2.c").read_text(errors="replace")
    def cut2(start):
        a = ican.rindex(start)
        return ican[a:ican.index("\n}\n", a) + 3]
    nodefuncs = [cut2(x) for x in ["int NodeSetCreate( struct tagCANON_GLOBALS *pCG,\n", "void NodeSetFree( struct tagCANON_GLOBALS *pCG, NodeSet *pSet )\n{",
        "void NodeSetFromRadEndpoints( CANON_GLOBALS *pCG,\n", "void RemoveFromNodeSet( CANON_GLOBALS *pCG,\n", "int DoNodeSetsIntersect( NodeSet *cur_nodes, int k1, int k2 )\n{",
        "int IsNodeSetEmpty( NodeSet *cur_nodes, int k )\n{", "void AddNodeSet2ToNodeSet1( NodeSet *cur_nodes, int k1, int k2 )\n{", "int AddNodesToRadEndpoints( CANON_GLOBALS *pCG,\n"]]
    tt = (SRC / "ichitaut.c").read_text(errors="replace")
    def cut_first(start):
        a = tt.index(start)
        return tt[a:tt.index("\n}\n", a) + 3]

    def cut3(start):
        a = tt.rindex(start)
        return tt[a:tt.index("\n}\n", a) + 3]
    ctype = tt[tt.index("typedef struct tagChargeType"):tt.index("#define NUM_C_TYPES")] + "#define NUM_C_TYPES  (int)(sizeof( CType )/sizeof(CType[0]))\n"
    endfuncs = [ctype, cut3("int CmpCCandidates( const void *a1, const void *a2 )\n{"), cut3("int RegisterCPoints( C_GROUP *c_group,\n"), cut3("int MakeIsotopicHGroup( inp_ATOM *at,\n"), cut3("int GetSaltChargeType( inp_ATOM *at,\n"), cut3("int GetOtherSaltChargeType( inp_ATOM *at,\n"), cut3("int GetOtherSaltType( inp_ATOM *at, int at_no, int *s_subtype )\n{"), cut3("int bDoNotMergeNonTautAtom( inp_ATOM *at, int at_no )\n{"), cut3("int comp_candidates( const void *a1, const void *a2 )\n{"), cut3("int RegisterEndPoints( CANON_GLOBALS *pCG,\n"), cut3("int AddEndPoint( T_ENDPOINT *pEndPoint, inp_ATOM *at, int iat )\n{"), cut3("int AddAtom2num( AT_RANK num[], inp_ATOM *atom, int at_no, int bSubtract )\n{"), cut3("void AddAtom2DA( AT_RANK num_DA[], inp_ATOM *atom, int at_no, int bSubtract )\n{"), cut3("int is_centerpoint_elem_KET( U_CHAR el_number )\n{"), cut3("int is_centerpoint_elem_strict( U_CHAR el_number )\n{"),
        cut3("int bCanBeACPoint( inp_ATOM *at,\n"), cut3("int GetChargeType( inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype )\n{")]
    for v in ["", "_PT_22_00", "_PT_16_00", "_PT_06_00", "_PT_39_00", "_PT_13_00", "_PT_18_00"]:
        sig = "int nGetEndpointInfo%s(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)" % v if v else "int nGetEndpointInfo( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif )"
        endfuncs.append(cut3(sig))
    endfuncs.append(cut3("int nGetEndpointInfo_KET( inp_ATOM *atom,\n"))
    for sig in ("int bHasAcidicHydrogen( inp_ATOM *at, int i )\n{", "int bHasAcidicMinus( inp_ATOM *at, int i )\n{", "int bHasOtherExchangableH( inp_ATOM *at, int i )\n{"):
        k = text.rindex(sig)
        endfuncs.append(text[k:text.index("\n}\n", k) + 3])
    st = (SRC / "ichister.c").read_text(errors="replace")
    a2 = st.rindex("int get_opposite_sb_atom(")
    opp = st[a2:st.index("\n}\n", a2) + 3]
    arrs = []
    for arr in ("int ArTypMask[] =", "int AaTypMask[] ="):
        i0 = text.rindex(arr)
        arrs.append(text[i0:text.index("};", i0) + 3])
    isrt = (SRC / "ichisort.c").read_text(errors="replace")
    a3 = isrt.rindex("int insertions_sort( void *pCG,")
    swp = isrt.rindex("void inchi_swap(")
    names = list(dict.fromkeys(names))
    funcs = endfuncs[:1] + [isrt[swp:isrt.index("\n}\n", swp) + 3], isrt[a3:isrt.index("\n}\n", a3) + 3]] + [cut(n) for n in names] + [opp] + nodefuncs + endfuncs[1:]
    # the hard proton functions need the arrays before the functions
    funcs = funcs[:1] + arrs + funcs[1:]
    taut = (SRC / "ichitaut.c").read_text(errors="replace")
    a = taut.rindex("int is_centerpoint_elem( U_CHAR el_number )\n{")
    funcs.append(taut[a:taut.index("\n}\n", a) + 3])
    a = text.index("#define TREE_NOT_IN_M")
    b = text.index("#define TREE_MARK")
    macros = text[a:text.index("\n", b)]
    pa = text.index("#define PR_SIMPLE_MSK")
    macros += "\n" + text[pa:text.index("\n", text.index("#define AA_HARD_TYP_H"))] + "\n"
    ea = text.index("typedef enum tagAtTypeTotals")
    eb = text.index("#define ATBIT_Errors")
    macros += "\n#define BNS_CHK_ALTP_NO_ALTPATH 0\n#define BNS_CHK_ALTP_SAME_TGROUP 1\n#define BNS_CHK_ALTP_SAME_VERTEX 2\n#define BNS_CHK_ALTP_SET_SUCCESS 4\n" + text[ea:text.index("\n", eb)] + "\n#define MAX_ALT_AATG_ARRAY_LEN 127\n#define AATG_MARK_IN_PATH 1\n"
    alt = re.search(r"typedef struct tagAltPathChanges.*?\} ALT_PATH_CHANGES;", text, re.S).group(0)
    macros += "\n#define MAX_LOCAL_TGNUM 0\n#ifndef FLAG_FORCE_SALT_TAUT\n#endif\n"
    protos = """int MakeIsotopicHGroup( inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info, T_GROUP_INFO *t_group_info );
int bHasOtherExchangableH( inp_ATOM *at, int i );
int GetSaltChargeType( inp_ATOM *at, int at_no, T_GROUP_INFO *t_group_info, int *s_subtype );
int GetOtherSaltChargeType( inp_ATOM *at, int at_no, T_GROUP_INFO *t_group_info, int *s_subtype, int bAccept_O );
int GetOtherSaltType( inp_ATOM *at, int at_no, int *s_subtype );
int bDoNotMergeNonTautAtom( inp_ATOM *at, int at_no );
int comp_candidates( const void *a1, const void *a2 );
int RegisterEndPoints( CANON_GLOBALS *pCG, T_GROUP_INFO *t_group_info, T_ENDPOINT *EndPoint, int nNumEndPoints, inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi, struct BalancedNetworkStructure *pBNS );
int AddEndPoint( T_ENDPOINT *pEndPoint, inp_ATOM *at, int iat );
int AddAtom2num( AT_RANK num[], inp_ATOM *atom, int at_no, int bSubtract );
void AddAtom2DA( AT_RANK num_DA[], inp_ATOM *atom, int at_no, int bSubtract );
int ReInitBnStruct( BN_STRUCT *pBNS, inp_ATOM *at, int num_at, int bRemoveGroupsFromAtoms );
int AddTGroups2BnStruct( CANON_GLOBALS *pCG, BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, T_GROUP_INFO *tgi );
int AddCGroups2BnStruct( CANON_GLOBALS *pCG, BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi );
int insertions_sort( void *pCG, void *base, size_t num, size_t width, int( *compare )( const void *, const void *, void * ) );
void inchi_swap( char *a, char *b, size_t width );
int CompTGroupNumber( const void *tg1, const void *tg2, void *p );
int CompCGroupNumber( const void *cg1, const void *cg2, void *p );
int bIsHardRemHCandidate( inp_ATOM *at, int i, int *cSubType );
int GetNeutralRepsIfNeeded( AT_NUMB *pri, AT_NUMB *prj, inp_ATOM *at, int num_atoms, T_ENDPOINT *EndPoint, int nNumEndPoints, C_GROUP_INFO *cgi );
int FindAccessibleEndPoints( CANON_GLOBALS *pCG, T_ENDPOINT *EndPoint, int *nNumEndPoints, T_BONDPOS *BondPos, int *nNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi, int taut_mode );
int BnsAdjustFlowBondsRad( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms );
int nBondsValenceInpAt( const inp_ATOM *at, int *nNumAltBonds, int *nNumWrongBonds );
int CmpCCandidates( const void *a1, const void *a2 );
int RegisterCPoints( C_GROUP *c_group, int *pnum_c, int max_num_c, T_GROUP_INFO *t_group_info, int point1, int point2, int ctype, inp_ATOM *at, int num_atoms );
int MarkChargeGroups( struct tagCANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, C_GROUP_INFO *c_group_info, T_GROUP_INFO *t_group_info, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD );
int CreateCGroupInBnStruct( inp_ATOM *at, int num_atoms, BN_STRUCT *pBNS, int nType, int nMask, int nCharge );
int CreateTGroupInBnStruct( inp_ATOM *at, int num_atoms, BN_STRUCT *pBNS, int nType, int nMask );
int RemoveLastGroupFromBnStruct( inp_ATOM *at, int num_atoms, int tg, BN_STRUCT *pBNS );
int SimpleRemoveHplusNPO( inp_ATOM *at, int num_atoms, int nAtTypeTotals[], T_GROUP_INFO *t_group_info );
int SimpleRemoveAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2remove );
int SimpleAddAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2add );
int HardRemoveAcidicProtons( CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2remove, int *nNumCanceledCharges, BN_STRUCT *pBNS, BN_DATA *pBD );
int HardAddAcidicProtons( CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2add, int *nNumCanceledCharges, BN_STRUCT *pBNS, BN_DATA *pBD );
int HardRemoveHplusNP( CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, int bCancelChargesAlways, int *nNumCanceledCharges, BN_AATG *pAATG, BN_STRUCT *pBNS, BN_DATA *pBD );
int mark_at_type( inp_ATOM *atom, int num_atoms, int nAtTypeTotals[] );
int RemoveNPProtonsAndAcidCharges( CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, BN_AATG *pAATG, BN_STRUCT *pBNS, BN_DATA *pBD );
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
    path = ROOT / "tests/data/inchi_isoh.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
