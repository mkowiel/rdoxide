//! Differential tests of the InChI -> structure restoration (`ichirvr*.rs`)
//! against the C library, with fixtures from `tools/inchi_rvr*_reference.py`.

use rdoxide::external::inchi::ichi::InChI;
use rdoxide::external::inchi::ichi_io::{
    INCHI_IOS_TYPE_STRING, InchiIostream, IosFile, inchi_ios_init,
};
use rdoxide::external::inchi::ichicano::TGSO_SYMM_IORDER;
use rdoxide::external::inchi::ichiread::*;
use rdoxide::external::inchi::ichirvr1::*;
use rdoxide::external::inchi::ichirvrs::*;
use rdoxide::external::inchi::ichitaut::TGroupInfo;
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
