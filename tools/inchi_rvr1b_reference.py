"""Generates tests/data/inchi_rvr1b.json from the InChI C code: small
helpers of ichirvr1.c -- `BondFlowMaxcapMinorder`, `AtomStcapStflow`,
`nTautEndpointEdgeCap`, `get_pVA_atom_type`, `comp_cc_cand`,
`RegisterTCGroup` and the EDGE_LIST functions.

    python tools/inchi_rvr1b_reference.py
"""
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


r1 = load("rvr1", "tools/inchi_rvr1_reference.py")
r2, ma, SRC = r1.r2, r1.ma, r1.SRC

MAIN = r1.PRE + r"""
static int cmpcc(const void *a, const void *b) { return comp_cc_cand(a, b); }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int t;
        printf("#CASE %d\n", cs);
        scanf("%d", &t);
        if (t == 30) {
            SRM s; memset(&s, 0, sizeof(s));
            scanf("%d %d %d %d %d %d %d %d %d", &s.bMetalAddFlower, &s.nMetalMinBondOrder, &s.nMetalInitEdgeFlow, &s.nMetalInitBondOrder,
                  &s.nMetal2EndpointMinBondOrder, &s.nMetal2EndpointInitBondOrder, &s.nMetal2EndpointInitEdgeFlow, &s.nMetalMaxCharge_D, &s.bFixStereoBonds);
            int n; scanf("%d", &n);
            inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
            VAL_AT *va = calloc(n + 1, sizeof(VAL_AT));
            for (int i = 0; i < n; i++) {
                int el, nh, cbv, e, ch, rad, pp, nsb, nb, v[8];
                scanf("%d %d %d %d %d %d %d %d", &el, &nh, &cbv, &e, &ch, &rad, &pp, &nsb);
                at[i].el_number = el; at[i].num_H = nh; at[i].chem_bonds_valence = cbv; at[i].endpoint = e;
                at[i].charge = ch; at[i].radical = rad; at[i].p_parity = pp;
                for (int k = 0; k < nsb; k++) { int p, o; scanf("%d %d", &p, &o); at[i].sb_parity[k] = p; at[i].sb_ord[k] = o; }
                scanf("%d", &nb); at[i].valence = nb;
                for (int k = 0; k < nb; k++) { int nn, bt; scanf("%d %d", &nn, &bt); at[i].neighbor[k] = nn; at[i].bond_type[k] = bt; }
                for (int k = 0; k < 8; k++) scanf("%d", &v[k]);
                va[i].cMetal = v[0]; va[i].cInitValenceToMetal = v[1]; va[i].cInitOrigValenceToMetal = v[2]; va[i].cInitFreeValences = v[3];
                va[i].cnListIndex = v[4]; va[i].cNumValenceElectrons = v[5]; va[i].cPeriodicRowNumber = v[6]; va[i].cPeriodicNumber = v[7];
            }
            for (int i = 0; i < n; i++) {
                for (int k = 0; k < at[i].valence; k++) {
                    int mc, mo, nf; int fl = BondFlowMaxcapMinorder(at, va, &s, i, k, &mc, &mo, &nf);
                    printf("B%d.%d %d %d %d %d\n", i, k, fl, mc, mo, nf);
                }
                int stcap, stflow; EdgeFlow mcap, mflow;
                int r = AtomStcapStflow(at, va, &s, i, &stcap, &stflow, &mcap, &mflow);
                printf("S%d %d %d %d %d %d\n", i, r, stcap, stflow, mcap, mflow);
                printf("T%d %d\n", i, nTautEndpointEdgeCap(at, va, i));
                printf("P%d %d %d %d\n", i, get_pVA_atom_type(va, at, i, 1), get_pVA_atom_type(va, at, i, 2), get_pVA_atom_type(va, at, i, 0));
            }
        } else if (t == 31) {
            int n; scanf("%d", &n);
            CC_CAND *c = calloc(n, sizeof(CC_CAND));
            for (int i = 0; i < n; i++) {
                int v[8]; for (int k = 0; k < 8; k++) scanf("%d", &v[k]);
                c[i].iat = v[0]; c[i].num_bonds = v[1]; c[i].chem_valence = v[2]; c[i].cMetal = v[3]; c[i].cNumBondsToMetal = v[4];
                c[i].cNumValenceElectrons = v[5]; c[i].cPeriodicRowNumber = v[6]; c[i].cNumChargeStates = 0; c[i].el_number = v[7];
            }
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) printf("%d ", comp_cc_cand(&c[i], &c[j]));
            printf("\n");
        } else if (t == 32) {
            ALL_TC_GROUPS g; memset(&g, 0, sizeof(g));
            int m; scanf("%d", &m);
            for (int i = 0; i < m; i++) {
                int v[7]; for (int k = 0; k < 7; k++) scanf("%d", &v[k]);
                int r = RegisterTCGroup(&g, v[0], v[1], v[2], v[3], v[4], v[5], v[6]);
                printf("r%d %d %d %d\n", i, r, g.num_tc_groups, g.max_tc_groups);
            }
            for (int i = 0; i < g.num_tc_groups; i++) printf("g%d %d %d %d %d %d %d %d\n", i, g.pTCG[i].type, g.pTCG[i].ord_num, g.pTCG[i].num_edges, g.pTCG[i].st_cap, g.pTCG[i].st_flow, g.pTCG[i].edges_cap, g.pTCG[i].edges_flow);
        } else if (t == 34 || t == 40 || t == 41) {
            int xa = 0, xe = 0, coa = 0, fsel = 0;
            if (t >= 40) scanf("%d %d %d", &xa, &xe, &coa);
            if (t == 41) scanf("%d", &fsel);
            SRM s; memset(&s, 0, sizeof(s));
            scanf("%d %d %d %d %d %d %d %d %d", &s.bMetalAddFlower, &s.nMetalMinBondOrder, &s.nMetalInitEdgeFlow, &s.nMetalInitBondOrder,
                  &s.nMetal2EndpointMinBondOrder, &s.nMetal2EndpointInitBondOrder, &s.nMetal2EndpointInitEdgeFlow, &s.nMetalMaxCharge_D, &s.bFixStereoBonds);
            s.nMetalFlowerParam_D = s.nMetalMaxCharge_D ? 16 : 3;
            int n, ng, tch; scanf("%d %d %d", &n, &ng, &tch);
            T_GROUP_INFO ti; memset(&ti, 0, sizeof(ti));
            ti.num_t_groups = ng;
            ti.t_group = calloc(ng + 1, sizeof(T_GROUP));
            for (int g = 0; g < ng; g++) {
                int a, b, ne; scanf("%d %d %d", &a, &b, &ne);
                ti.t_group[g].num[0] = a; ti.t_group[g].num[1] = b; ti.t_group[g].nGroupNumber = g + 1; ti.t_group[g].nNumEndpoints = ne;
            }
            inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
            VAL_AT *va = calloc(n + 1, sizeof(VAL_AT));
            for (int i = 0; i < n; i++) {
                int el, nh, cbv, e, nb, nsb, v[7];
                scanf("%d %d %d %d %d", &el, &nh, &cbv, &e, &nb);
                at[i].el_number = el; at[i].num_H = nh; at[i].chem_bonds_valence = cbv; at[i].endpoint = e; at[i].valence = nb;
                for (int k = 0; k < nb; k++) { int nn, bt; scanf("%d %d", &nn, &bt); at[i].neighbor[k] = nn; at[i].bond_type[k] = bt; }
                for (int k = 0; k < 7; k++) scanf("%d", &v[k]);
                va[i].cMetal = v[0]; va[i].cNumBondsToMetal = v[1]; va[i].cInitValenceToMetal = v[2]; va[i].cInitOrigValenceToMetal = v[3];
                va[i].cInitFreeValences = v[4]; va[i].cnListIndex = v[5]; va[i].cPeriodicNumber = el; if (t == 41) { int row; va[i].cNumValenceElectrons = get_sp_element_type(el, &row); va[i].cPeriodicRowNumber = row; va[i].cNumValenceElectrons--; }
            }
            ALL_TC_GROUPS g; memset(&g, 0, sizeof(g));
            for (int k = 0; k < NUM_TCGROUP_TYPES; k++) g.nGroup[k] = -1;
            g.total_charge = tch; g.charge_on_atoms = coa;
            int r = nCountBnsSizes(at, n, 0, 0, &ti, va, &s, &g);
            printf("C %d | %d %d %d %d %d %d %d %d %d %d | %d %d\n", r, g.num_atoms, g.num_bonds, g.num_tgroups, g.num_tgroup_edges, g.tgroup_charge,
                   g.nVertices, g.nEdges, g.num_metal_atoms, g.num_metal_bonds, g.total_electrons, g.total_electrons_metals, g.num_tc_groups);
            for (int i = 0; i < g.num_tc_groups; i++) printf("g%d %d %d %d %d %d %d %d | %d %d %d %d\n", i, g.pTCG[i].type, g.pTCG[i].ord_num, g.pTCG[i].num_edges, g.pTCG[i].st_cap, g.pTCG[i].st_flow, g.pTCG[i].edges_cap, g.pTCG[i].edges_flow,
                   g.pTCG[i].tg_num_H, g.pTCG[i].tg_num_Minus, 0, 0);
            if (r >= 0) {
                int r2 = nAddSuperCGroups(&g);
                printf("U %d %d %d %d |", r2, g.nVertices, g.nEdges, g.num_tc_groups);
                for (int k = 0; k < NUM_TCGROUP_TYPES; k++) printf(" %d", g.nGroup[k]);
                printf("\n");
                for (int i = 0; i < g.num_tc_groups; i++) printf("u%d %d %d %d\n", i, g.pTCG[i].type, g.pTCG[i].ord_num, g.pTCG[i].num_edges);
                if (t >= 40 && r2 >= 0) {
                    StrFromINChI st; memset(&st, 0, sizeof(st)); st.at = at; st.num_atoms = n; st.pSrm = &s;
                    int nch = -1;
                    BN_STRUCT *b = AllocateAndInitTCGBnStruct(&st, va, &g, xa, xe, 1, &nch);
                    if (!b) { printf("NULL\n"); }
                    else {
                        int r3 = AddTGroups2TCGBnStruct(b, &st, va, &g, xe);
                        printf("T %d\n", r3);
                        int r4 = r3 < 0 ? -1 : AddCGroups2TCGBnStruct(b, &st, va, &g, xe);
                        printf("Q %d %d\n", r4, nch);
                        printf("N %d %d %d %d %d %d %d %d %d %d %d\n", b->num_atoms, b->num_added_atoms, b->num_c_groups, b->num_t_groups, b->num_vertices, b->num_bonds, b->num_edges,
                               b->max_vertices, b->max_edges, b->max_iedges, b->tot_st_cap * 1000 + b->tot_st_flow);
                        for (int i = 0; i < b->num_vertices; i++) {
                            BNS_VERTEX *v = b->vert + i;
                            printf("V%d %d %d %d %d %d %d %d %d %d |", i, v->st_edge.cap, v->st_edge.cap0, v->st_edge.flow, v->st_edge.flow0, v->type, v->num_adj_edges, v->max_adj_edges, (int)(v->iedge - b->iedge), 0);
                            for (int k = 0; k < v->num_adj_edges; k++) printf(" %d", v->iedge[k]);
                            printf("\n");
                        }
                        for (int i = 0; i < b->num_edges; i++) {
                            BNS_EDGE *e = b->edge + i;
                            printf("E%d %d %d %d %d %d %d %d %d %d %d\n", i, e->neighbor1, e->neighbor12, e->neigh_ord[0], e->neigh_ord[1], e->cap, e->cap0, e->flow, e->flow0, e->pass, e->forbidden);
                        }
                        for (int i = 0; i < n; i++) printf("W%d %d %d %d %d | %d %d\n", i, va[i].nCMinusGroupEdge, va[i].nCPlusGroupEdge, va[i].nMetalGroupEdge, va[i].nTautGroupEdge, at[i].bond_type[0], at[i].valence ? at[i].bond_type[at[i].valence - 1] : 0);
                        for (int i = 0; i < g.num_tc_groups; i++) printf("X%d %d %d %d %d %d %d %d %d\n", i, g.pTCG[i].nVertexNumber, g.pTCG[i].nForwardEdge, g.pTCG[i].nBackwardEdge, g.pTCG[i].st_cap, g.pTCG[i].st_flow, g.pTCG[i].edges_cap, g.pTCG[i].edges_flow, g.pTCG[i].num_edges);
                        printf("Y %d %d %d %d\n", g.nEdgePlus, g.nEdgeMinus, g.nEdge4charge, g.added_charge);
                        if (r4 >= 0 && t == 41) {
                            for (int j = 0; j < at[0].valence && j < 4; j++) if (fsel >> j & 1) b->edge[b->vert[0].iedge[j]].forbidden |= 4;
                            int nc = AdjustTgroupsToForbiddenEdges2(b, at, va, n, 4);
                            printf("A %d %d %d\n", nc, b->tot_st_cap, b->tot_st_flow);
                            for (int i = 0; i < b->num_edges; i++) printf("p%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                            for (int i = 0; i < b->num_vertices; i++) printf("q%d %d\n", i, b->vert[i].st_edge.flow);
                        }
                        if (r4 >= 0) {
                            for (int i = 0; i < b->num_edges; i++) printf("U%d %d\n", i, GetChargeFlowerUpperEdge(b, va, i));
                            EDGE_LIST el; memset(&el, 0, sizeof(el));
                            for (int i = 0; i < b->num_edges; i += 2) AddToEdgeList(&el, i, 4);
                            SetForbiddenEdgeMask(b, &el, 2);
                            printf("F"); for (int i = 0; i < b->num_edges; i++) printf(" %d", b->edge[i].forbidden); printf("\n");
                            RemoveForbiddenEdgeMask(b, &el, 1);
                            printf("F"); for (int i = 0; i < b->num_edges; i++) printf(" %d", b->edge[i].forbidden); printf("\n");
                            RemoveForbiddenBondFlowBits(b, 2);
                            printf("F"); for (int i = 0; i < b->num_edges; i++) printf(" %d", b->edge[i].forbidden); printf("\n");
                            BN_DATA *bd = AllocateAndInitBnData(b->max_vertices);
                            Vertex vf, vl; int pl, dh = -7, dc = -7, nv = -7;
                            int rt = RunBnsTestOnce(b, bd, va, &vf, &vl, &pl, &dh, &dc, &nv);
                            printf("R %d %d %d %d %d %d %d\n", rt, vf, vl, pl, dh, dc, nv);
                            int rr = RunBnsRestoreOnce(b, bd, va, &g);
                            printf("S %d %d %d\n", rr, b->tot_st_cap, b->tot_st_flow);
                            for (int i = 0; i < b->num_edges; i++) printf("e%d %d %d\n", i, b->edge[i].flow, b->edge[i].cap);
                            for (int i = 0; i < b->num_vertices; i++) printf("v%d %d %d\n", i, b->vert[i].st_edge.flow, b->vert[i].st_edge.cap);
                            for (int zz = 0; zz < 2; zz++) {
                                inp_ATOM *ac = calloc(n + 1, sizeof(inp_ATOM));
                                memcpy(ac, at, n * sizeof(inp_ATOM));
                                for (int i = 0; i < g.num_tc_groups && i < b->num_t_groups; i++) { g.pTCG[i].tg_RestoreFlags = (n + i) & 1; g.pTCG[i].tg_set_Minus = (n * 3 + i) % (n + 2); }
                                StrFromINChI s2; memset(&s2, 0, sizeof(s2)); s2.at = ac; s2.num_atoms = n; s2.pSrm = &s;
                                int rc = CopyBnsToAtom(&s2, b, va, &g, zz);
                                printf("K%d %d\n", zz, rc);
                                for (int i = 0; i < n; i++) {
                                    printf("k%d %d %d %d %d %d %d |", i, ac[i].charge, ac[i].radical, ac[i].num_H, ac[i].endpoint, ac[i].chem_bonds_valence, ac[i].valence);
                                    for (int q = 0; q < ac[i].valence; q++) printf(" %d", ac[i].bond_type[q]);
                                    printf("\n");
                                }
                            }
                            {
                                inp_ATOM *a2 = calloc(n + 1, sizeof(inp_ATOM));
                                StrFromINChI s3; memset(&s3, 0, sizeof(s3)); s3.at = at; s3.num_atoms = n; s3.pSrm = &s;
                                int rm = MoveRadToAtomsAddCharges(b, bd, &s3, at, a2, va, &g, 4);
                                printf("M %d %d %d %d\n", rm, b->tot_st_cap, b->tot_st_flow, b->edge_forbidden_mask);
                                for (int i = 0; i < b->num_vertices; i++) printf("m%d %d %d\n", i, b->vert[i].st_edge.cap, b->vert[i].st_edge.flow);
                                for (int i = 0; i < b->num_edges; i++) printf("f%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                for (int i = 0; i < n; i++) printf("c%d %d %d %d %d\n", i, va[i].cInitCharge, a2[i].valence, a2[i].chem_bonds_valence, a2[i].charge);
                                int rra = RearrangePlusMinusEdgesFlow(b, bd, va, &g, 4);
                                printf("RA %d %d %d\n", rra, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("g%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                int nr = 0, td = 0;
                                int riz = IncrementZeroOrderBondsToHeteroat(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                printf("IZ %d %d %d %d %d\n", riz, nr, td, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("h%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                BFS_Q bq; memset(&bq, 0, sizeof(bq)); AllocBfsQueue(&bq, n, 0); s3.pbfsq = &bq;
                                nr = 0; td = 0;
                                int rmp = MovePlusFromS2DiaminoCarbon(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                printf("MP %d %d %d %d\n", rmp, nr, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("i%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                if ((n & 1) && b->num_edges > 0) b->edge[0].forbidden = 8;
                                int rec = EliminateChargeSeparationOnHeteroatoms(b, bd, &s3, at, a2, va, &g, &nr, &td, 4, (n & 1) ? 8 : 0);
                                printf("EC %d %d %d %d %d\n", rec, nr, td, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("j%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                nr = 0; td = 0;
                                int rcy = RestoreCyanoGroup(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                printf("CY %d %d %d %d %d\n", rcy, nr, td, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("k%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                int ric = RestoreIsoCyanoGroup(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                printf("IC %d %d %d %d %d\n", ric, nr, td, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("l%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                int rfm = FixMetal_Nminus_Ominus(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                printf("FM %d %d %d %d %d\n", rfm, nr, td, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("m%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                int rnn = RestoreNNNgroup(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                printf("NN %d %d %d %d %d\n", rnn, nr, td, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("o%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                int rp = EliminateNitrogen5Val3Bonds(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                printf("E5 %d %d %d %d %d\n", rp, nr, td, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("q%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                int rr = Convert_SIV_to_SVI(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                printf("SV %d %d %d %d %d\n", rr, nr, td, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("s%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                int rt = PlusFromDB_N_DB_O_to_Metal(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                printf("PM %d %d %d %d %d\n", rt, nr, td, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("u%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                int rv = MoveMobileHToAvoidFixedBonds(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                printf("MM %d %d %d %d %d\n", rv, nr, td, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("w%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                int rx = RemoveRadFromMobileHEndpoint(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                printf("RR %d %d %d %d %d\n", rx, nr, td, b->tot_st_cap, b->tot_st_flow);
                                for (int i = 0; i < b->num_edges; i++) printf("y%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                {
                                    int rmc = MoveChargeToMakeCenerpoints(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                    long long chk = 0;
                                    for (int i = 0; i < b->num_edges; i++) chk = (chk + (long long)(b->edge[i].flow * 31 + b->edge[i].forbidden) * (i + 1)) % 1000000007LL;
                                    printf("MC %d %d %d %d %d %lld\n", rmc, nr, td, b->tot_st_cap, b->tot_st_flow, chk);
                                }
                                {
                                    T_GROUP tgx; memset(&tgx, 0, sizeof(tgx));
                                    AT_NUMB *epl = calloc(n + 1, sizeof(AT_NUMB)), *epf = calloc(n + 1, sizeof(AT_NUMB));
                                    int ne = 0;
                                    for (int i = 0; i < n; i++) { if (at[i].endpoint) epl[ne++] = i; epf[i] = at[i].endpoint; }
                                    tgx.nNumEndpoints = ne;
                                    s3.ti.t_group = &tgx; s3.ti.num_t_groups = ne ? 1 : 0; s3.ti.nNumEndpoints = ne;
                                    s3.ti.nEndpointAtomNumber = epl; s3.endpoint = epf; s3.iMobileH = 0;
                                    int rfh = RemoveRadFromMobileHEndpointFixH(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                    printf("FH %d %d %d %d %d %d\n", rfh, nr, td, b->tot_st_cap, b->tot_st_flow, ne);
                                    for (int i = 0; i < b->num_edges; i++) printf("V%d %d %d\n", i, b->edge[i].flow, b->edge[i].forbidden);
                                    {
                                        EDGE_LIST nl; AllocEdgeList(&nl, EDGE_LIST_CLEAR);
                                        int rn = ForbidNintrogenPlus2BondsInSmallRings(b, a2, n, va, 6, &g, &nl, 4);
                                        long long chk = 0; for (int q = 0; q < b->num_edges; q++) chk = (chk + (long long)(b->edge[q].flow * 31 + b->edge[q].forbidden) * (q + 1)) % 1000000007LL;
                                        printf("FN %d %d %lld\n", rn, nl.num_edges, chk);
                                        RemoveForbiddenEdgeMask(b, &nl, 4);
                                        AllocEdgeList(&nl, EDGE_LIST_FREE);
                                        inp_ATOM *atf = calloc(n + 1, sizeof(inp_ATOM)), *atn = calloc(n + 1, sizeof(inp_ATOM));
                                        for (int pass = 0; pass < 3; pass++) {
                                            memcpy(atf, a2, n * sizeof(inp_ATOM)); memcpy(atn, a2, n * sizeof(inp_ATOM));
                                            int rr = 0;
                                            if (pass == 0) {
                                                for (int i = 0; i < n; i++) {
                                                    if (i % 2 == 0 && atf[i].num_H > 0) atf[i].num_H--;
                                                    for (int j = 0; j < atf[i].valence; j++) if ((i + j) & 1) {
                                                        int bt = atf[i].bond_type[j], m = bt & BOND_TYPE_MASK;
                                                        if (m == 1) atf[i].bond_type[j] = bt + 1; else if (m == 2 || m == 3) atf[i].bond_type[j] = bt - 1;
                                                    }
                                                }
                                                rr = FixLessHydrogenInFormula(b, bd, &s3, at, a2, atf, va, &g, &nr, &td, 4);
                                            } else if (pass == 1) {
                                                for (int i = 0; i < n; i++) if (a2[i].charge == -1) { atf[i].charge = 0; atf[i].num_H = a2[i].num_H + 1; }
                                                s3.bMobileH = n & 1;
                                                rr = FixMoreHydrogenInFormula(b, bd, &s3, at, a2, atf, va, &g, &nr, &td, 4);
                                            } else {
                                                ICR icr; memset(&icr, 0, sizeof(icr));
                                                for (int i = 0; i < n && i < 32; i++) icr.endp_in1_only[icr.num_endp_in1_only++] = i + 1;
                                                for (int i = 0; i < n; i++) atn[i].endpoint = 1;
                                                rr = FixRemoveExtraTautEndpoints(b, bd, &s3, at, a2, atf, atn, va, &g, &icr, &nr, &td, 4);
                                            }
                                            chk = 0; for (int q = 0; q < b->num_edges; q++) chk = (chk + (long long)(b->edge[q].flow * 31 + b->edge[q].forbidden) * (q + 1)) % 1000000007LL;
                                            int cf = 0; for (int i = 0; i < n; i++) cf += a2[i].cFlags;
                                            printf("FX%d %d %d %d %d %d %lld\n", pass, rr, nr, td, b->tot_st_flow, cf, chk);
                                        }
                                        s3.charge = (n % 3) - 1; s3.iMobileH = ((n >> 1) & 1) ? 0 : 1;
                                        s3.fixed_H = calloc(n + 1, 1); for (int i = 0; i < n; i++) s3.fixed_H[i] = (i % 7 == 0);
                                        int rmr = MoveChargeToRemoveCenerpoints(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                        chk = 0; for (int q = 0; q < b->num_edges; q++) chk = (chk + (long long)(b->edge[q].flow * 31 + b->edge[q].forbidden) * (q + 1)) % 1000000007LL;
                                        int ci = 0; for (int i = 0; i < n; i++) ci += va[i].cInitCharge * (i + 1);
                                        printf("MR %d %d %d %d %d %d %lld\n", rmr, nr, td, b->tot_st_cap, b->tot_st_flow, ci, chk);
                                        int rms2 = MakeSingleBondsMetal2ChargedHeteroat(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                        chk = 0; for (int q = 0; q < b->num_edges; q++) chk = (chk + (long long)(b->edge[q].flow * 31 + b->edge[q].forbidden) * (q + 1)) % 1000000007LL;
                                        printf("MS %d %d %d %d %d %lld\n", rms2, nr, td, b->tot_st_cap, b->tot_st_flow, chk);
                                        s3.iInchiRec = 1;
                                        int rsb = SaltBondsToCoordBonds(b, bd, &s3, at, a2, va, &g, &nr, &td, 4);
                                        chk = 0; for (int q = 0; q < b->num_edges; q++) chk = (chk + (long long)(b->edge[q].flow * 31 + b->edge[q].forbidden + b->edge[q].cap * 7) * (q + 1)) % 1000000007LL;
                                        ci = 0; for (int i = 0; i < n; i++) ci += va[i].cInitCharge * (i + 1);
                                        printf("SB %d %d %d %d %d %d %d %lld\n", rsb, nr, td, b->tot_st_cap, b->tot_st_flow, ci, s.bMetalAddFlower, chk);
                                        s3.iInchiRec = 0;
                                        free(s3.fixed_H); s3.fixed_H = NULL; s3.iMobileH = 0; s3.charge = 0;
                                        free(atf); free(atn);
                                    }
                                    s3.endpoint = NULL; s3.ti.t_group = NULL; free(epl); free(epf);
                                }
                            }
                            for (int i = 0; i < n; i++) at[i].bond_type[0] = 0;
                            IncrZeroBondsAndClearEndpts(at, n, 5);
                            for (int i = 0; i < n; i++) printf("z%d %d %d %d %d\n", i, at[i].endpoint, at[i].component, at[i].chem_bonds_valence, at[i].bond_type[0]);
                        }
                    }
                }
            }
        } else if (t == 50 || t == 51) {
            int n, nd; scanf("%d %d", &n, &nd);
            inp_ATOM *at = calloc(n + nd + 1, sizeof(inp_ATOM));
            for (int i = 0; i < n + nd; i++) {
                int orig, nh, nb, cbv, iso[3], iad, nsb;
                scanf("%d %d %d", &orig, &nh, &nb);
                at[i].orig_at_number = orig; at[i].num_H = nh; at[i].valence = nb;
                for (int k = 0; k < nb; k++) { int nn, bt, bs; scanf("%d %d %d", &nn, &bt, &bs); at[i].neighbor[k] = nn; at[i].bond_type[k] = bt; at[i].bond_stereo[k] = bs; }
                scanf("%d %d %d %d %d %d", &cbv, &iso[0], &iso[1], &iso[2], &iad, &nsb);
                at[i].chem_bonds_valence = cbv; at[i].num_iso_H[0] = iso[0]; at[i].num_iso_H[1] = iso[1]; at[i].num_iso_H[2] = iso[2]; at[i].iso_atw_diff = iad;
                for (int k = 0; k < nsb; k++) { int a, b, c, d; scanf("%d %d %d %d", &a, &b, &c, &d); at[i].sb_parity[k] = a; at[i].sb_ord[k] = b; at[i].sn_ord[k] = c; at[i].sn_orig_at_num[k] = d; }
            }
            int r = t == 50 ? ConnectDisconnectedH(at, n, nd) : DisconnectedConnectedH(at, n, nd);
            printf("H %d\n", r);
            for (int i = 0; i < n + nd; i++) {
                printf("a%d %d %d %d %d %d %d %d |", i, at[i].valence, at[i].chem_bonds_valence, at[i].num_H, at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2], at[i].iso_atw_diff);
                for (int k = 0; k < MAXVAL; k++) printf(" %d,%d,%d", at[i].neighbor[k], at[i].bond_type[k], at[i].bond_stereo[k]);
                printf(" |");
                for (int k = 0; k < 3; k++) printf(" %d,%d,%d,%d", at[i].sb_parity[k], at[i].sb_ord[k], at[i].sn_ord[k], at[i].sn_orig_at_num[k]);
                printf("\n");
            }
        } else if (t == 33) {
            EDGE_LIST e; memset(&e, 0, sizeof(e));
            int m; scanf("%d", &m);
            for (int i = 0; i < m; i++) {
                int op, a, b; scanf("%d %d %d", &op, &a, &b);
                int r = 0;
                if (op == 0) r = AddToEdgeList(&e, a, b);
                else if (op == 1) r = RemoveFromEdgeListByIndex(&e, a);
                else if (op == 2) r = FindInEdgeList(&e, a);
                else if (op == 3) r = RemoveFromEdgeListByValue(&e, a);
                else if (op == 4) r = AllocEdgeList(&e, a);
                printf("o%d %d | %d %d |", i, r, e.num_alloc, e.num_edges);
                for (int k = 0; k < e.num_alloc; k++) printf(" %d", e.pnEdges[k]);
                printf("\n");
            }
        }
    }
    return 0;
}
"""
HARNESS = r2.HARNESS.split("int main(void)")[0] + MAIN
R = random.Random(97531)
ELS = [1, 5, 6, 6, 6, 7, 7, 8, 8, 9, 14, 15, 16, 16, 17, 34, 35, 53, 11, 12, 26, 29, 13, 3, 33, 52, 20, 30, 50, 82, 79, 92]


def atom_case():
    n = R.randrange(1, 6)
    adj = [[] for _ in range(n)]
    edges = {}
    for i in range(1, n):
        edges[(R.randrange(0, i), i)] = R.choice([1, 1, 2, 3, 4, 5])
    for (i, j), b in edges.items():
        adj[i].append((j, b))
        adj[j].append((i, b))
    srm = [R.choice([0, 1, 1]), R.choice([0, 1]), R.choice([0, 1]), R.choice([0, 1, 2]), R.choice([0, 1]), R.choice([0, 1, 2]),
           R.choice([0, 1]), R.choice([0, 16]), R.choice([0, 1])]
    parts = ["30 " + " ".join(map(str, srm)) + f" {n}"]
    for i in range(n):
        el = R.choice(ELS)
        nsb = R.choice([0, 0, 1, 2]) if adj[i] else 0
        sb = []
        for _ in range(nsb):
            sb += [R.choice([1, 2]), R.randrange(0, len(adj[i]))]
        va = [R.choice([0, 0, 1]), R.choice([0, 1, 2]), R.choice([0, 1, 2, 3]), R.choice([0, 0, 1, 2]), R.choice([0, 0, 1, 3, 18]),
              R.choice([0, 3, 4, 5, 6, 7]), R.choice([0, 1, 2, 3, 5]), el]
        parts.append(f"{el} {R.choice([0, 1, 2])} {R.randrange(0, 5)} {R.choice([0, 0, 1])} {R.choice([-1, 0, 0, 1])} {R.choice([0, 0, 2])} "
                     f"{R.choice([0, 0, 1])} {nsb} " + " ".join(map(str, sb)) + f" {len(adj[i])} " + " ".join(f"{j} {b}" for j, b in adj[i]) + " " + " ".join(map(str, va)))
    return " ".join(parts)


def cc_case():
    n = R.randrange(1, 6)
    return f"31 {n} " + " ".join(" ".join(str(R.choice(c)) for c in
                                         ([0, 1, 5], [0, 1, 2, 3], [0, 1, 2, 3], [0, 1], [0, 1, 2], [0, 4, 5, 6], [0, 1, 2, 3], [6, 7, 8])) for _ in range(n))


def tcg_case():
    m = R.randrange(1, 40)
    return f"32 {m} " + " ".join(" ".join(str(R.randrange(a, b)) for a, b in ((0, 4), (0, 40), (-2, 3), (-2, 3), (-2, 3), (-2, 3), (0, 3))) for _ in range(m))


def el_case():
    m = R.randrange(1, 25)
    ops = []
    for _ in range(m):
        op = R.choice([0, 0, 0, 0, 1, 2, 3, 4])
        a = R.randrange(0, 6) if op != 4 else R.choice([-2, -1, 0, 3, 5, 8])
        b = R.choice([-1, 0, 1, 2, 4]) if op == 0 else 0
        ops.append(f"{op} {a} {b}")
    return f"33 {m} " + " ".join(ops)


def bns_case(rad=False):
    n = R.randrange(1, 6)
    ng = 0 if rad else R.choice([0, 0, 1, 2, 3])
    adj = [[] for _ in range(n)]
    edges = {}
    for i in range(1, n):
        edges[(R.randrange(0, i), i)] = R.choice([1, 1, 2, 3])
    for (i, j), b in edges.items():
        adj[i].append((j, b))
        adj[j].append((i, b))
    srm = [R.choice([0, 1, 1]), R.choice([0, 1]), R.choice([0, 1]), R.choice([0, 1, 2]), R.choice([0, 1]), R.choice([0, 1, 2]),
           R.choice([0, 1]), R.choice([0, 16]), 0]
    eps = [0 if rad else R.choice([0, 0, 0] + list(range(1, ng + 1))) for _ in range(n)]
    cnt = [sum(1 for e in eps if e == g) for g in range(1, ng + 1)]
    if R.random() < 0.7:
        tgs = [(R.randrange(0, 4), R.randrange(0, 3), cnt[g]) for g in range(ng)]
    else:
        tgs = [(R.randrange(0, 4), R.randrange(0, 3), R.randrange(0, 4)) for _ in range(ng)]
    parts = ["34 " + " ".join(map(str, srm)) + f" {n} {ng} {R.choice([-1, 0, 0, 1, 2])}"]
    for a, b, ne in tgs:
        parts.append(f"{a} {b} {ne}")
    for i in range(n):
        el = R.choice(ELS)
        va = [R.choice([0, 0, 1]), R.choice([0, 1, 2]), R.choice([0, 1]), R.choice([0, 1, 2]), R.choice([0, 0, 1, 2]),
              R.choice([0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18])]
        if rad:
            va = [0, 0, 0, 0, R.choice([0, 1, 2, 3]), R.choice([0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])]
        parts.append(f"{el} {R.choice([0, 1, 2])} {R.randrange(0, 5)} {eps[i]} {len(adj[i])} " + " ".join(f"{j} {b}" for j, b in adj[i])
                     + " " + " ".join(map(str, va)) + " 0")
    return " ".join(parts)


def pipe_case(rad=False):
    c = bns_case(rad).split(" ", 1)
    return "40 " + f"{R.choice([0, 0, 2])} {R.choice([0, 1, 2])} {R.choice([0, 0, 1, -1])} " + c[1]


def h_case():
    n = R.randrange(1, 4)
    adj = [[] for _ in range(n)]
    for i in range(1, n):
        j = R.randrange(0, i)
        b = R.choice([1, 1, 2])
        adj[i].append((j, b))
        adj[j].append((i, b))
    groups = []  # (parent, [iso diffs])
    parents = list(range(n))
    if R.random() < 0.2:
        R.shuffle(parents)
    for k in parents:
        if R.random() < 0.7:
            e = R.randrange(1, 4)
            diffs = [0] * R.randrange(0, e + 1)
            diffs += sorted(R.choice([1, 2, 3]) for _ in range(e - len(diffs))) if R.random() < 0.8 else [R.choice([1, 2, 3, 4])] * (e - len(diffs))
            groups.append((k, diffs))
    nd = sum(len(d) for _, d in groups)
    conn = R.random() < 0.5
    hidx = {}
    idx = n
    for k, diffs in groups:
        hidx[k] = list(range(idx, idx + len(diffs)))
        idx += len(diffs)
    atoms = []
    for i in range(n):
        iso = [0, 0, 0]
        e = 0
        if i in hidx:
            e = len(hidx[i])
            for gk, diffs in groups:
                if gk == i:
                    for d in diffs:
                        if 1 <= d <= 3 and R.random() < 0.93:
                            iso[d - 1] += 1
        nh = e + R.randrange(0, 3) + sum(iso) if R.random() < 0.9 else R.randrange(0, 3)
        neigh = [(j, b, R.choice([0, 0, 1, 6])) for j, b in adj[i]]
        if conn:
            neigh = [(h, 1, 0) for h in hidx.get(i, [])] + neigh
        nsb = R.choice([0, 0, 1, 2])
        sb = []
        for _ in range(nsb):
            sn = R.choice([-1, -1, 0, 1])
            sb += [R.choice([1, 2]), R.randrange(0, 2), sn, R.choice([100 + h for h in hidx.get(i, [0])] + [1])]
        cbv = sum(b for _, b, _ in neigh) + (0 if conn else R.choice([0, 0, 1]))
        atoms.append(f"{i + 1} {nh} {len(neigh)} " + " ".join(f"{a} {b} {c}" for a, b, c in neigh) + f" {cbv} {iso[0]} {iso[1]} {iso[2]} 0 {nsb} " + " ".join(map(str, sb)))
    for k, diffs in groups:
        for q, d in enumerate(diffs):
            atoms.append(f"{100 + hidx[k][q]} 0 1 {k} 1 0 1 0 0 0 {d} 0")
    return f"{51 if conn else 50} {n} {nd} " + " ".join(atoms)


def star_case():
    k = R.choice([3, 3, 4])
    n = 1 + k
    ng = R.choice([1, 2])
    center = R.choice([16, 16, 16, 7, 15, 6])
    nels = [R.choice([7, 7, 8, 8, 16]) for _ in range(k)]
    bonds = [R.choice([1, 1, 2]) for _ in range(k)]
    eps = [R.choice([0, 1, 1, 2][:2 + ng]) for _ in range(k)]
    cnt = [sum(1 for e in eps if e == g) for g in range(1, ng + 1)]
    tgs = [(R.randrange(0, 4), R.randrange(0, 2), cnt[g]) for g in range(ng)]
    srm = [1, 0, 1, 1, 1, 1, 0, 16, 0]
    cbv_c = sum(bonds)
    parts = ["41 " + f"{R.choice([0, 2])} {R.choice([0, 1])} 0 {R.randrange(0, 16)} " + " ".join(map(str, srm)) + f" {n} {ng} 0"]
    for a, b, ne in tgs:
        parts.append(f"{a} {b} {ne}")
    parts.append(f"{center} 0 {cbv_c} 0 {k} " + " ".join(f"{j + 1} {bonds[j]}" for j in range(k)) + " 0 0 0 0 0 0 0")
    for j in range(k):
        nh = R.choice([0, 1, 1, 2]) if bonds[j] == 1 else R.choice([0, 1])
        parts.append(f"{nels[j]} {nh} {bonds[j]} {eps[j]} 1 0 {bonds[j]} 0 0 0 0 0 0 0")
    return " ".join(parts)


def tmpl_case():
    """centre with 3 mobile-H neighbours as in the Table 5 cases 2, 7 and 4"""
    kind = R.choice(["5.2", "5.7", "5.4"])
    if kind == "5.2":
        center, nb, forb = 7, [(7, 1, 2), (8, 2, 0), (8, 2, 0)], 0
    elif kind == "5.7":
        center, nb, forb = 16, [(7, 1, 2), (8, 2, 0), (8, 1, 1)], 0
    else:
        center, nb, forb = 7, [(8, 2, 0), (7, 1, 2), (7, 2, 1)], 0
    order = list(range(3))
    R.shuffle(order)
    nb = [nb[o] for o in order]
    forb = order.index(forb)
    fsel = 1 << forb
    if R.random() < 0.2:
        fsel |= 1 << R.randrange(3)
    k = 3
    n = 4
    bonds = [b for _, b, _ in nb]
    srm = [1, 0, 1, 1, 1, 1, 0, 16, 0]
    parts = ["41 " + f"{R.choice([0, 2])} {R.choice([0, 1])} 0 {fsel} " + " ".join(map(str, srm)) + f" {n} 1 0"]
    parts.append(f"{R.randrange(0, 4)} {R.randrange(0, 2)} 3")
    parts.append(f"{center} 0 {sum(bonds)} 0 {k} " + " ".join(f"{j + 1} {bonds[j]}" for j in range(k)) + " 0 0 0 0 0 0 0")
    for j, (el, b, nh) in enumerate(nb):
        parts.append(f"{el} {nh} {b} 1 1 0 {b} 0 0 0 0 0 0 0")
    return " ".join(parts)


def cyano_case():
    # N(-)=C=X with the N charge structure "MN" (cnList index 14)
    x = R.choice([8, 8, 16, 7])
    srm = [1, 0, 1, 1, 1, 1, 0, 16, 0]
    parts = ["40 " + f"{R.choice([0, 2])} {R.choice([0, 1])} {R.choice([0, -1])} " + " ".join(map(str, srm)) + f" 3 0 {R.choice([-1, 0])}"]
    parts.append(f"7 0 2 0 1 1 2 0 0 0 0 0 {R.choice([14, 14, 4, 0])} 0")
    parts.append("6 0 4 0 2 0 2 2 2 0 0 0 0 0 0 0")
    parts.append(f"{x} 0 2 0 1 1 2 0 0 0 0 0 0 0")
    return " ".join(parts)


def main():
    cases = []
    for _ in range(3000):
        cases += [atom_case(), cc_case(), tcg_case(), el_case(), bns_case(), bns_case(), pipe_case(), pipe_case(), pipe_case(), h_case(), h_case(), star_case(), tmpl_case(), tmpl_case(), cyano_case(), cyano_case()]
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        kept, out = [], []
        for c in cases:
            try:
                r1_ = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=2)
            except subprocess.TimeoutExpired:
                continue
            if r1_.returncode != 0:
                continue
            kept.append(c)
            out += [f"#CASE {len(kept) - 1}"] + r1_.stdout.decode(errors="replace").split("\n")[1:]
        print(len(cases) - len(kept), "cases dropped")
        cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_rvr1b.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
