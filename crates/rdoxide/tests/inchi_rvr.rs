//! Differential tests of the InChI -> structure restoration (`ichirvr*.rs`)
//! against the C library, with fixtures from `tools/inchi_rvr*_reference.py`.

use rdoxide::external::inchi::ichi::InChI;
use rdoxide::external::inchi::ichi_bns::{BnStruct, allocate_and_init_bn_data};
use rdoxide::external::inchi::ichi_io::{
    INCHI_IOS_TYPE_STRING, InchiIostream, IosFile, inchi_ios_init,
};
use rdoxide::external::inchi::ichicano::TGSO_SYMM_IORDER;
use rdoxide::external::inchi::ichiread::*;
use rdoxide::external::inchi::ichirvr1::*;
use rdoxide::external::inchi::ichirvrs::*;
use rdoxide::external::inchi::ichitaut::{TGroup, TGroupInfo};
use rdoxide::external::inchi::inpdef::InpAtom;
use rdoxide::external::inchi::util::is_el_a_metal;

fn unhex(h: &str) -> Vec<u8> {
    if h == "=" {
        return Vec::new();
    }
    (0..h.len() / 2)
        .map(|i| u8::from_str_radix(&h[2 * i..2 * i + 2], 16).unwrap())
        .collect()
}

fn dump_va(v: &ValAt) -> String {
    format!(
        "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
        v.c_do_not_add_h,
        v.c_metal,
        v.c_num_bonds_to_metal,
        v.c_init_flow_to_metal,
        v.c_init_valence_to_metal,
        v.c_init_orig_valence_to_metal,
        v.c_max_flow_to_metal,
        v.c_init_free_valences,
        v.c_init_charge,
        v.c_num_valence_electrons,
        v.c_periodic_row_number,
        v.c_min_ring_size,
        v.c_periodic_number,
        v.c_cn_list_index,
        v.n_c_minus_group_edge,
        v.n_c_plus_group_edge,
        v.n_metal_group_edge,
        v.n_taut_group_edge
    )
}

fn dump_ti(ti: &TGroupInfo) -> String {
    let mut s = format!(
        "TI {} {} {} | ",
        ti.num_t_groups, ti.n_num_endpoints, ti.max_num_t_groups
    );
    if let Some(tg) = &ti.t_group {
        for g in tg.iter().take(ti.num_t_groups) {
            s.push_str(&format!(
                "({} {} {} {} {}) ",
                g.num[0],
                g.num[1],
                g.n_group_number,
                g.n_num_endpoints,
                g.n_first_endpoint_at_no_pos
            ));
        }
    }
    s.push_str("| ");
    for k in 0..ti.n_num_endpoints.max(0) as usize {
        if k < ti.n_endpoint_atom_number.len() {
            s.push_str(&format!("{},", ti.n_endpoint_atom_number[k]));
        }
    }
    s.push_str("| ");
    if ti.t_group_number.is_empty() {
        s.push_str("NULL");
    } else {
        for i in 0..ti.num_t_groups {
            s.push_str(&format!(
                "{},{} ",
                ti.t_group_number[i],
                ti.t_group_number[TGSO_SYMM_IORDER * ti.num_t_groups + i]
            ));
        }
    }
    s
}

#[test]
fn atom_restore_info_and_tgroups_match_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_rvr1.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let t: Vec<&str> = inp.split(' ').collect();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let p = |k: usize| t[k].parse::<i32>().unwrap();
        let mut got: Vec<String> = Vec::new();
        match t[0] {
            "20" => {
                let mob = p(1);
                let n = p(2) as usize;
                let mut k = 3;
                let mut at = vec![InpAtom::default(); n + 1];
                let mut ep = vec![0u16; n + 1];
                for i in 0..n {
                    at[i].el_number = p(k) as u8;
                    at[i].num_h = p(k + 1) as i8;
                    at[i].chem_bonds_valence = p(k + 2) as i8;
                    at[i].endpoint = p(k + 3) as u16;
                    ep[i] = p(k + 3) as u16;
                    let nb = p(k + 4) as usize;
                    at[i].valence = nb as i8;
                    k += 5;
                    for q in 0..nb {
                        at[i].neighbor[q] = p(k) as u16;
                        at[i].bond_type[q] = p(k + 1) as u8;
                        k += 2;
                    }
                }
                let mut va = vec![ValAt::default(); n + 1];
                let srm = set_up_srm();
                for i in 0..n {
                    let mut row = 0;
                    va[i].c_num_valence_electrons =
                        get_sp_element_type(i32::from(at[i].el_number), &mut row) as i8;
                    va[i].c_periodic_row_number = row as i8;
                    va[i].c_periodic_number = at[i].el_number;
                    va[i].c_num_valence_electrons -= 1;
                    if is_el_a_metal(i32::from(at[i].el_number)) {
                        va[i].c_metal = 1;
                    }
                }
                for i in 0..n {
                    va[i].c_init_free_valences = 0;
                    let r = get_atom_restore_info(&mut at, i, &mut va, &srm, mob, Some(&ep));
                    got.push(format!(
                        "A{i} {r} {} | {}",
                        at[i].chem_bonds_valence,
                        dump_va(&va[i])
                    ));
                }
            }
            "21" => {
                let mode = p(1);
                let bytes = unhex(t[2]);
                let mut ios = InchiIostream::default();
                inchi_ios_init(&mut ios, INCHI_IOS_TYPE_STRING, IosFile::None);
                ios.s.n_used_length = bytes.len() as i32;
                ios.s.p_str = Some(bytes);
                ios.s.n_ptr = 0;
                let mut line = SegmLine::default();
                let mut pstr = None;
                let (mut state, mut nerr) = (-77, -9);
                let mut d = InpInChI::default();
                let (mut stdf, mut hso, mut bits) = (-5, -5, 0u8);
                let r = inchi_line2data(
                    &mut ios, &mut line, &mut pstr, &mut state, &mut nerr, &mut d, true, 0,
                    &mut stdf, &mut hso, &mut bits,
                );
                let ok = (r == RI_ERR_EOL || r == RI_ERR_EOF)
                    && nerr == 0
                    && d.n_num_components[0][0] + d.n_num_components[0][1] != 0;
                got.push(format!("R {}", i32::from(ok)));
                if ok {
                    for a in 0..2 {
                        for b in 0..2 {
                            for i in 0..d.n_num_components[a][b].max(0) as usize {
                                let Some(c) = d.p_inp_inchi[a][b].get(i) else {
                                    break;
                                };
                                if c.n_number_of_atoms == 0 {
                                    continue;
                                }
                                let na = c.n_number_of_atoms as usize;
                                let mut at = vec![InpAtom::default(); na + 2];
                                let mut ep = vec![0u16; na + 2];
                                let mut ti = TGroupInfo::default();
                                for pass in 0..2 {
                                    let rr = get_tgroup_info_from_inchi(
                                        &mut ti,
                                        if mode & 1 != 0 {
                                            Some(&mut at[..])
                                        } else {
                                            None
                                        },
                                        if mode & 2 != 0 {
                                            Some(&mut ep[..])
                                        } else {
                                            None
                                        },
                                        Some(c),
                                    );
                                    got.push(format!("G{a}{b}.{i}.{pass} {rr} {}", dump_ti(&ti)));
                                    let mut e = String::from("E ");
                                    for k in 0..na {
                                        e.push_str(&format!("{},", at[k].endpoint));
                                    }
                                    e.push_str(" | ");
                                    for k in 0..na {
                                        e.push_str(&format!("{},", ep[k]));
                                    }
                                    got.push(e);
                                }
                                let mut ep2 = None;
                                let rr = fill_outp_struct_endpoint_from_inchi(c, &mut ep2);
                                let mut f = format!("F {rr} ");
                                for k in 0..na {
                                    f.push_str(&format!("{},", ep2.as_ref().unwrap()[k]));
                                }
                                got.push(f);
                            }
                        }
                    }
                }
                let _: Option<&InChI> = None;
            }
            "22" => {
                let mut s = String::new();
                for e in 0..130 {
                    let mut row = 0;
                    let ty = get_sp_element_type(e, &mut row);
                    s.push_str(&format!("{ty},{row} "));
                }
                got.push(s);
            }
            _ => unreachable!(),
        }
        assert_eq!(got, want, "case {ci}: {inp}");
    }
}

#[test]
fn restoration_helpers_match_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_rvr1b.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let t: Vec<i32> = inp.split_whitespace().map(|x| x.parse().unwrap()).collect();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut got: Vec<String> = Vec::new();
        let mut k = 1;
        let mut next = || {
            k += 1;
            t[k - 1]
        };
        match t[0] {
            30 => {
                let s = Srm {
                    b_metal_add_flower: next(),
                    n_metal_min_bond_order: next(),
                    n_metal_init_edge_flow: next(),
                    n_metal_init_bond_order: next(),
                    n_metal2_endpoint_min_bond_order: next(),
                    n_metal2_endpoint_init_bond_order: next(),
                    n_metal2_endpoint_init_edge_flow: next(),
                    n_metal_max_charge_d: next(),
                    b_fix_stereo_bonds: next(),
                    ..Default::default()
                };
                let n = next() as usize;
                let mut at = vec![InpAtom::default(); n + 1];
                let mut va = vec![ValAt::default(); n + 1];
                for i in 0..n {
                    at[i].el_number = next() as u8;
                    at[i].num_h = next() as i8;
                    at[i].chem_bonds_valence = next() as i8;
                    at[i].endpoint = next() as u16;
                    at[i].charge = next() as i8;
                    at[i].radical = next() as i8;
                    at[i].p_parity = next() as i8;
                    let nsb = next() as usize;
                    for q in 0..nsb {
                        at[i].sb_parity[q] = next() as i8;
                        at[i].sb_ord[q] = next() as i8;
                    }
                    let nb = next() as usize;
                    at[i].valence = nb as i8;
                    for q in 0..nb {
                        at[i].neighbor[q] = next() as u16;
                        at[i].bond_type[q] = next() as u8;
                    }
                    va[i].c_metal = next() as i8;
                    va[i].c_init_valence_to_metal = next() as i8;
                    va[i].c_init_orig_valence_to_metal = next() as i8;
                    va[i].c_init_free_valences = next() as i8;
                    va[i].c_cn_list_index = next() as i8;
                    va[i].c_num_valence_electrons = next() as i8;
                    va[i].c_periodic_row_number = next() as i8;
                    va[i].c_periodic_number = next() as u8;
                }
                for i in 0..n {
                    for q in 0..at[i].valence as usize {
                        let (fl, mc, mo, nf) = bond_flow_maxcap_minorder(&at, &va, &s, i, q);
                        got.push(format!("B{i}.{q} {fl} {mc} {mo} {}", i32::from(nf)));
                    }
                    let (r, stcap, stflow, mcap, mflow) = atom_stcap_stflow(&at, &va, &s, i);
                    got.push(format!("S{i} {r} {stcap} {stflow} {mcap} {mflow}"));
                    got.push(format!("T{i} {}", n_taut_endpoint_edge_cap(&at, &va, i)));
                    got.push(format!(
                        "P{i} {} {} {}",
                        get_pva_atom_type(&va, &at, i, 1),
                        get_pva_atom_type(&va, &at, i, 2),
                        get_pva_atom_type(&va, &at, i, 0)
                    ));
                }
            }
            31 => {
                let n = next() as usize;
                let mut cc = Vec::new();
                for _ in 0..n {
                    cc.push(CcCand {
                        iat: next() as i16,
                        num_bonds: next() as i8,
                        chem_valence: next() as i8,
                        c_metal: next() as i8,
                        c_num_bonds_to_metal: next() as i8,
                        c_num_valence_electrons: next() as i8,
                        c_periodic_row_number: next() as i8,
                        c_num_charge_states: 0,
                        el_number: next() as u8,
                    });
                }
                let mut s = String::new();
                for a in &cc {
                    for b in &cc {
                        s.push_str(&format!("{} ", comp_cc_cand(a, b)));
                    }
                }
                got.push(s);
            }
            32 => {
                let mut g = AllTcGroups::default();
                let m = next();
                for i in 0..m {
                    let v: Vec<i32> = (0..7).map(|_| next()).collect();
                    let r = register_tc_group(&mut g, v[0], v[1], v[2], v[3], v[4], v[5], v[6]);
                    got.push(format!("r{i} {r} {} {}", g.num_tc_groups, g.max_tc_groups));
                }
                for i in 0..g.num_tc_groups as usize {
                    let x = &g.p_tcg[i];
                    got.push(format!(
                        "g{i} {} {} {} {} {} {} {}",
                        x.type_,
                        x.ord_num,
                        x.num_edges,
                        x.st_cap,
                        x.st_flow,
                        x.edges_cap,
                        x.edges_flow
                    ));
                }
            }
            34 | 40 | 41 => {
                let (xa, xe, coa) = if t[0] >= 40 {
                    (next(), next(), next())
                } else {
                    (0, 0, 0)
                };
                let fsel = if t[0] == 41 { next() } else { 0 };
                let s = Srm {
                    b_metal_add_flower: next(),
                    n_metal_min_bond_order: next(),
                    n_metal_init_edge_flow: next(),
                    n_metal_init_bond_order: next(),
                    n_metal2_endpoint_min_bond_order: next(),
                    n_metal2_endpoint_init_bond_order: next(),
                    n_metal2_endpoint_init_edge_flow: next(),
                    n_metal_max_charge_d: next(),
                    b_fix_stereo_bonds: next(),
                    ..Default::default()
                };
                let s = Srm {
                    n_metal_flower_param_d: if s.n_metal_max_charge_d != 0 { 16 } else { 3 },
                    ..s
                };
                let n = next() as usize;
                let ng = next() as usize;
                let tch = next();
                let mut ti = TGroupInfo {
                    num_t_groups: ng,
                    ..Default::default()
                };
                let mut tgs = vec![TGroup::default(); ng + 1];
                for (g, tg) in tgs.iter_mut().take(ng).enumerate() {
                    tg.num[0] = next() as u16;
                    tg.num[1] = next() as u16;
                    tg.n_group_number = g as u16 + 1;
                    tg.n_num_endpoints = next() as u16;
                }
                ti.t_group = Some(tgs);
                let mut at = vec![InpAtom::default(); n + 1];
                let mut va = vec![ValAt::default(); n + 1];
                for i in 0..n {
                    at[i].el_number = next() as u8;
                    at[i].num_h = next() as i8;
                    at[i].chem_bonds_valence = next() as i8;
                    at[i].endpoint = next() as u16;
                    let nb = next() as usize;
                    at[i].valence = nb as i8;
                    for q in 0..nb {
                        at[i].neighbor[q] = next() as u16;
                        at[i].bond_type[q] = next() as u8;
                    }
                    va[i].c_metal = next() as i8;
                    va[i].c_num_bonds_to_metal = next() as i8;
                    va[i].c_init_valence_to_metal = next() as i8;
                    va[i].c_init_orig_valence_to_metal = next() as i8;
                    va[i].c_init_free_valences = next() as i8;
                    va[i].c_cn_list_index = next() as i8;
                    va[i].c_periodic_number = at[i].el_number;
                    if t[0] == 41 {
                        let mut row = 0;
                        va[i].c_num_valence_electrons =
                            get_sp_element_type(i32::from(at[i].el_number), &mut row) as i8;
                        va[i].c_periodic_row_number = row as i8;
                        va[i].c_num_valence_electrons -= 1;
                    }
                    next();
                }
                let mut g = AllTcGroups {
                    n_group: [-1; NUM_TCGROUP_TYPES],
                    total_charge: tch,
                    charge_on_atoms: coa,
                    ..Default::default()
                };
                let r = n_count_bns_sizes(&at, n, &ti, &va, &s, &mut g);
                got.push(format!(
                    "C {r} | {} {} {} {} {} {} {} {} {} {} | {} {}",
                    g.num_atoms,
                    g.num_bonds,
                    g.num_tgroups,
                    g.num_tgroup_edges,
                    g.tgroup_charge,
                    g.n_vertices,
                    g.n_edges,
                    g.num_metal_atoms,
                    g.num_metal_bonds,
                    g.total_electrons,
                    g.total_electrons_metals,
                    g.num_tc_groups
                ));
                for i in 0..g.num_tc_groups as usize {
                    let x = &g.p_tcg[i];
                    got.push(format!(
                        "g{i} {} {} {} {} {} {} {} | {} {} 0 0",
                        x.type_,
                        x.ord_num,
                        x.num_edges,
                        x.st_cap,
                        x.st_flow,
                        x.edges_cap,
                        x.edges_flow,
                        x.tg_num_h,
                        x.tg_num_minus
                    ));
                }
                if r >= 0 {
                    let r2 = n_add_super_c_groups(&mut g);
                    let mut u = format!(
                        "U {r2} {} {} {} |",
                        g.n_vertices, g.n_edges, g.num_tc_groups
                    );
                    for k in 0..NUM_TCGROUP_TYPES {
                        u.push_str(&format!(" {}", g.n_group[k]));
                    }
                    got.push(u);
                    for i in 0..g.num_tc_groups as usize {
                        let x = &g.p_tcg[i];
                        got.push(format!("u{i} {} {} {}", x.type_, x.ord_num, x.num_edges));
                    }
                    if t[0] >= 40 && r2 >= 0 {
                        match allocate_and_init_tcg_bn_struct(&mut at, n, &s, &va, &g, xa, xe, 1) {
                            None => got.push("NULL".into()),
                            Some((mut b, nch)) => {
                                let r3 = add_tgroups_to_tcg_bn_struct(
                                    &mut b, &at, n, &mut va, &mut g, xe,
                                );
                                got.push(format!("T {r3}"));
                                let r4 = if r3 < 0 {
                                    0
                                } else {
                                    add_cgroups_to_tcg_bn_struct(
                                        &mut b, &at, n, &s, &mut va, &mut g, xe,
                                    )
                                };
                                got.push(format!("Q {r4} {nch}"));
                                got.push(format!(
                                    "N {} {} {} {} {} {} {} {} {} {} {}",
                                    b.num_atoms,
                                    b.num_added_atoms,
                                    b.num_c_groups,
                                    b.num_t_groups,
                                    b.num_vertices,
                                    b.num_bonds,
                                    b.num_edges,
                                    b.max_vertices,
                                    b.max_edges,
                                    b.max_iedges,
                                    b.tot_st_cap * 1000 + b.tot_st_flow
                                ));
                                for i in 0..b.num_vertices as usize {
                                    let v = &b.vert[i];
                                    let mut l = format!(
                                        "V{i} {} {} {} {} {} {} {} {} 0 |",
                                        v.st_edge.cap,
                                        v.st_edge.cap0,
                                        v.st_edge.flow,
                                        v.st_edge.flow0,
                                        v.type_,
                                        v.num_adj_edges,
                                        v.max_adj_edges,
                                        v.iedge
                                    );
                                    for q in 0..v.num_adj_edges as usize {
                                        l.push_str(&format!(" {}", b.iedge[v.iedge + q]));
                                    }
                                    got.push(l);
                                }
                                for i in 0..b.num_edges as usize {
                                    let e = &b.edge[i];
                                    got.push(format!(
                                        "E{i} {} {} {} {} {} {} {} {} {} {}",
                                        e.neighbor1,
                                        e.neighbor12,
                                        e.neigh_ord[0],
                                        e.neigh_ord[1],
                                        e.cap,
                                        e.cap0,
                                        e.flow,
                                        e.flow0,
                                        e.pass,
                                        e.forbidden
                                    ));
                                }
                                for i in 0..n {
                                    got.push(format!(
                                        "W{i} {} {} {} {} | {} {}",
                                        va[i].n_c_minus_group_edge,
                                        va[i].n_c_plus_group_edge,
                                        va[i].n_metal_group_edge,
                                        va[i].n_taut_group_edge,
                                        at[i].bond_type[0],
                                        if at[i].valence > 0 {
                                            at[i].bond_type[at[i].valence as usize - 1]
                                        } else {
                                            0
                                        }
                                    ));
                                }
                                for i in 0..g.num_tc_groups as usize {
                                    let x = &g.p_tcg[i];
                                    got.push(format!(
                                        "X{i} {} {} {} {} {} {} {} {}",
                                        x.n_vertex_number,
                                        x.n_forward_edge,
                                        x.n_backward_edge,
                                        x.st_cap,
                                        x.st_flow,
                                        x.edges_cap,
                                        x.edges_flow,
                                        x.num_edges
                                    ));
                                }
                                got.push(format!(
                                    "Y {} {} {} {}",
                                    g.n_edge_plus, g.n_edge_minus, g.n_edge4charge, g.added_charge
                                ));
                                if r4 >= 0 && t[0] == 41 {
                                    for j in 0..(at[0].valence as usize).min(4) {
                                        if fsel >> j & 1 != 0 {
                                            let ie = b.iedge[b.vert[0].iedge + j] as usize;
                                            b.edge[ie].forbidden |= 4;
                                        }
                                    }
                                    let nc = rdoxide::external::inchi::ichirvr2::adjust_tgroups_to_forbidden_edges2(&mut b, &at, &va, n, 4);
                                    got.push(format!("A {nc} {} {}", b.tot_st_cap, b.tot_st_flow));
                                    for i in 0..b.num_edges as usize {
                                        got.push(format!(
                                            "p{i} {} {}",
                                            b.edge[i].flow, b.edge[i].forbidden
                                        ));
                                    }
                                    for i in 0..b.num_vertices as usize {
                                        got.push(format!("q{i} {}", b.vert[i].st_edge.flow));
                                    }
                                }
                                if r4 >= 0 {
                                    for i in 0..b.num_edges {
                                        got.push(format!(
                                            "U{i} {}",
                                            get_charge_flower_upper_edge(&b, i)
                                        ));
                                    }
                                    let mut el = EdgeList::default();
                                    for i in (0..b.num_edges).step_by(2) {
                                        add_to_edge_list(&mut el, i, 4);
                                    }
                                    let forb = |b: &BnStruct| {
                                        let mut l = String::from("F");
                                        for i in 0..b.num_edges as usize {
                                            l.push_str(&format!(" {}", b.edge[i].forbidden));
                                        }
                                        l
                                    };
                                    set_forbidden_edge_mask(&mut b, &el, 2);
                                    got.push(forb(&b));
                                    remove_forbidden_edge_mask(&mut b, &el, 1);
                                    got.push(forb(&b));
                                    remove_forbidden_bond_flow_bits(&mut b, 2);
                                    got.push(forb(&b));
                                    let mut bd = allocate_and_init_bn_data(b.max_vertices);
                                    let mut res = BnsTestResult {
                                        delta_h: -7,
                                        delta_charge: -7,
                                        num_visited_atoms: -7,
                                        ..Default::default()
                                    };
                                    let rt = run_bns_test_once(&mut b, &mut bd, &va, &mut res);
                                    got.push(format!(
                                        "R {rt} {} {} {} {} {} {}",
                                        res.v_first,
                                        res.v_last,
                                        res.path_len,
                                        res.delta_h,
                                        res.delta_charge,
                                        res.num_visited_atoms
                                    ));
                                    let rr = run_bns_restore_once(&mut b, &mut bd);
                                    got.push(format!("S {rr} {} {}", b.tot_st_cap, b.tot_st_flow));
                                    for i in 0..b.num_edges as usize {
                                        got.push(format!(
                                            "e{i} {} {}",
                                            b.edge[i].flow, b.edge[i].cap
                                        ));
                                    }
                                    for i in 0..b.num_vertices as usize {
                                        got.push(format!(
                                            "v{i} {} {}",
                                            b.vert[i].st_edge.flow, b.vert[i].st_edge.cap
                                        ));
                                    }
                                    for zz in 0..2 {
                                        let mut ps =
                                            rdoxide::external::inchi::ichirvr2::StrFromInChI {
                                                at: at.clone(),
                                                num_atoms: n as i32,
                                                srm: s,
                                                ..Default::default()
                                            };
                                        for i in 0..(g.num_tc_groups.min(b.num_t_groups)) as usize {
                                            g.p_tcg[i].tg_restore_flags = ((n + i) & 1) as i16;
                                            g.p_tcg[i].tg_set_minus =
                                                ((n * 3 + i) % (n + 2)) as i16;
                                        }
                                        let rc =
                                            rdoxide::external::inchi::ichirvr2::copy_bns_to_atom(
                                                &mut ps,
                                                &b,
                                                &va,
                                                &g,
                                                zz == 1,
                                            );
                                        got.push(format!("K{zz} {rc}"));
                                        for i in 0..n {
                                            let x = &ps.at[i];
                                            let mut l = format!(
                                                "k{i} {} {} {} {} {} {} |",
                                                x.charge,
                                                x.radical,
                                                x.num_h,
                                                x.endpoint,
                                                x.chem_bonds_valence,
                                                x.valence
                                            );
                                            for q in 0..x.valence as usize {
                                                l.push_str(&format!(" {}", x.bond_type[q]));
                                            }
                                            got.push(l);
                                        }
                                    }
                                    {
                                        use rdoxide::external::inchi::ichirvr2::*;
                                        let mut ps = StrFromInChI {
                                            at: at.clone(),
                                            num_atoms: n as i32,
                                            srm: s,
                                            ..Default::default()
                                        };
                                        let mut a2: Vec<_> = Vec::new();
                                        let rm = move_rad_to_atoms_add_charges(
                                            &mut b, &mut bd, &mut ps, &mut a2, &mut va, &g, 4,
                                        );
                                        got.push(format!(
                                            "M {rm} {} {} {}",
                                            b.tot_st_cap, b.tot_st_flow, b.edge_forbidden_mask
                                        ));
                                        for i in 0..b.num_vertices as usize {
                                            got.push(format!(
                                                "m{i} {} {}",
                                                b.vert[i].st_edge.cap, b.vert[i].st_edge.flow
                                            ));
                                        }
                                        for i in 0..b.num_edges as usize {
                                            got.push(format!(
                                                "f{i} {} {}",
                                                b.edge[i].flow, b.edge[i].forbidden
                                            ));
                                        }
                                        for i in 0..n {
                                            let (v, cb, ch) = if a2.is_empty() {
                                                (0, 0, 0)
                                            } else {
                                                (
                                                    a2[i].valence,
                                                    a2[i].chem_bonds_valence,
                                                    a2[i].charge,
                                                )
                                            };
                                            got.push(format!(
                                                "c{i} {} {} {} {}",
                                                va[i].c_init_charge, v, cb, ch
                                            ));
                                        }
                                    }
                                    for a in at.iter_mut().take(n) {
                                        a.bond_type[0] = 0;
                                    }
                                    incr_zero_bonds_and_clear_endpts(&mut at, n, 5);
                                    for i in 0..n {
                                        got.push(format!(
                                            "z{i} {} {} {} {}",
                                            at[i].endpoint,
                                            at[i].component,
                                            at[i].chem_bonds_valence,
                                            at[i].bond_type[0]
                                        ));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            50 | 51 => {
                let n = next() as usize;
                let nd = next() as usize;
                let mut at = vec![InpAtom::default(); n + nd + 1];
                for a in at.iter_mut().take(n + nd) {
                    a.orig_at_number = next() as u16;
                    a.num_h = next() as i8;
                    let nb = next() as usize;
                    a.valence = nb as i8;
                    for q in 0..nb {
                        a.neighbor[q] = next() as u16;
                        a.bond_type[q] = next() as u8;
                        a.bond_stereo[q] = next() as i8;
                    }
                    a.chem_bonds_valence = next() as i8;
                    for q in 0..3 {
                        a.num_iso_h[q] = next() as i8;
                    }
                    a.iso_atw_diff = next() as i8;
                    let nsb = next() as usize;
                    for q in 0..nsb {
                        a.sb_parity[q] = next() as i8;
                        a.sb_ord[q] = next() as i8;
                        a.sn_ord[q] = next() as i8;
                        a.sn_orig_at_num[q] = next() as u16;
                    }
                }
                let r = if t[0] == 50 {
                    connect_disconnected_h(&mut at, n, nd)
                } else {
                    disconnected_connected_h(&mut at, n, nd)
                };
                got.push(format!("H {r}"));
                for (i, a) in at.iter().enumerate().take(n + nd) {
                    let mut l = format!(
                        "a{i} {} {} {} {} {} {} {} |",
                        a.valence,
                        a.chem_bonds_valence,
                        a.num_h,
                        a.num_iso_h[0],
                        a.num_iso_h[1],
                        a.num_iso_h[2],
                        a.iso_atw_diff
                    );
                    for q in 0..20 {
                        l.push_str(&format!(
                            " {},{},{}",
                            a.neighbor[q], a.bond_type[q], a.bond_stereo[q]
                        ));
                    }
                    l.push_str(" |");
                    for q in 0..3 {
                        l.push_str(&format!(
                            " {},{},{},{}",
                            a.sb_parity[q], a.sb_ord[q], a.sn_ord[q], a.sn_orig_at_num[q]
                        ));
                    }
                    got.push(l);
                }
            }
            33 => {
                let mut e = EdgeList::default();
                let m = next();
                for i in 0..m {
                    let (op, a, b) = (next(), next(), next());
                    let r = match op {
                        0 => add_to_edge_list(&mut e, a, b),
                        1 => remove_from_edge_list_by_index(&mut e, a),
                        2 => find_in_edge_list(&e, a),
                        3 => remove_from_edge_list_by_value(&mut e, a),
                        _ => alloc_edge_list(&mut e, a),
                    };
                    let mut s = format!("o{i} {r} | {} {} |", e.num_alloc, e.num_edges);
                    for q in 0..e.num_alloc as usize {
                        s.push_str(&format!(" {}", e.pn_edges[q]));
                    }
                    got.push(s);
                }
            }
            _ => unreachable!(),
        }
        if got != want {
            let d = got
                .iter()
                .zip(&want)
                .position(|(a, b)| a != *b)
                .unwrap_or(got.len().min(want.len()));
            let more: Vec<String> = got
                .iter()
                .zip(&want)
                .enumerate()
                .filter(|(_, (a, b))| a != b)
                .take(8)
                .map(|(i, (a, b))| format!("{i}: {a}  <>  {b}"))
                .collect();
            panic!(
                "case {ci}: {inp}\nline {d}: got {:?}\n want {:?}\n{}",
                got.get(d),
                want.get(d),
                more.join("\n")
            );
        }
    }
}

#[test]
fn restore_atom_connections_matches_c() {
    use rdoxide::external::inchi::ichirvr2::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_rvr2.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let t: Vec<&str> = inp.split(' ').collect();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mode: i32 = t[1].parse().unwrap();
        let bytes = unhex(t[2]);
        let mut ios = InchiIostream::default();
        inchi_ios_init(&mut ios, INCHI_IOS_TYPE_STRING, IosFile::None);
        ios.s.n_used_length = bytes.len() as i32;
        ios.s.p_str = Some(bytes);
        ios.s.n_ptr = 0;
        let mut line = SegmLine::default();
        let mut pstr = None;
        let (mut state, mut nerr) = (-77, -9);
        let mut d = InpInChI::default();
        let (mut stdf, mut hso, mut bits) = (-5, -5, 0u8);
        let r = inchi_line2data(
            &mut ios, &mut line, &mut pstr, &mut state, &mut nerr, &mut d, true, 0, &mut stdf,
            &mut hso, &mut bits,
        );
        let ok = (r == RI_ERR_EOL || r == RI_ERR_EOF)
            && nerr == 0
            && d.n_num_components[0][0] + d.n_num_components[0][1] != 0;
        let mut got = vec![format!("R {}", i32::from(ok))];
        if ok {
            for a in 0..2 {
                for b in 0..2 {
                    for i in 0..d.n_num_components[a][b].max(0) as usize {
                        let Some(comp) = d.p_inp_inchi[a][b].get(i) else {
                            break;
                        };
                        if comp.n_number_of_atoms == 0 {
                            continue;
                        }
                        let mob = if mode & 1 != 0
                            && b == 0
                            && (i as i32) < d.n_num_components[a][1]
                            && d.p_inp_inchi[a][1]
                                .get(i)
                                .is_some_and(|m| m.n_number_of_atoms != 0)
                        {
                            d.p_inp_inchi[a][1].get(i)
                        } else {
                            None
                        };
                        let mut s = StrFromInChI::default();
                        let rr = restore_atom_connections_set_stereo(
                            &mut s,
                            i as i32,
                            if mode & 2 != 0 { 5 } else { 0 },
                            comp,
                            mob,
                        );
                        got.push(format!("G{a}{b}.{i} {rr}"));
                        if rr > 0 {
                            let r2 = set_stereo_bond_types_from_0d_stereo(&mut s, comp);
                            got.push(format!("B {r2}"));
                            let tot = (s.num_atoms + s.num_deleted_h) as usize;
                            got.push(format!(
                                "N {} {} {} {} {}",
                                s.num_atoms,
                                s.num_deleted_h,
                                s.charge,
                                s.b_isotopic,
                                s.ti.num_t_groups
                            ));
                            for (k, x) in s.at.iter().enumerate().take(tot) {
                                let mut l = format!(
                                    "a{k} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} |",
                                    x.el_number,
                                    x.elname_str(),
                                    x.orig_at_number,
                                    x.orig_compt_at_numb,
                                    x.component,
                                    x.num_h,
                                    x.valence,
                                    x.chem_bonds_valence,
                                    x.at_type,
                                    x.endpoint,
                                    x.iso_atw_diff,
                                    x.num_iso_h[0],
                                    x.num_iso_h[1],
                                    x.num_iso_h[2],
                                    x.p_parity
                                );
                                for q in 0..x.valence as usize {
                                    l.push_str(&format!(" {},{}", x.neighbor[q], x.bond_type[q]));
                                }
                                l.push_str(" |");
                                for q in 0..4 {
                                    l.push_str(&format!(" {}", x.p_orig_at_num[q]));
                                }
                                l.push_str(" |");
                                for q in 0..3 {
                                    l.push_str(&format!(
                                        " {},{},{},{}",
                                        x.sb_parity[q],
                                        x.sb_ord[q],
                                        x.sn_ord[q],
                                        x.sn_orig_at_num[q]
                                    ));
                                }
                                got.push(l);
                            }
                            if let Some(st) = &s.st {
                                for (k, x) in st.iter().enumerate().take(s.num_atoms as usize) {
                                    let mut l = format!("s{k} {} |", x.p_parity);
                                    for q in 0..4 {
                                        l.push_str(&format!(" {}", x.p_orig_at_num[q]));
                                    }
                                    l.push_str(" |");
                                    for q in 0..3 {
                                        l.push_str(&format!(
                                            " {},{},{},{}",
                                            x.sb_parity[q],
                                            x.sb_ord[q],
                                            x.sn_ord[q],
                                            x.sn_orig_at_num[q]
                                        ));
                                    }
                                    got.push(l);
                                }
                            }
                        }
                    }
                }
            }
        }
        if got != want {
            let dd = got
                .iter()
                .zip(&want)
                .position(|(a, b)| a != *b)
                .unwrap_or(got.len().min(want.len()));
            panic!(
                "case {ci}: {inp}\nline {dd}: got {:?}\n want {:?}",
                got.get(dd),
                want.get(dd)
            );
        }
    }
}
