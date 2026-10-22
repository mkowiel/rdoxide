//! Differential tests of the InChI string reader (`ichiread.rs`) against the
//! C library, with fixtures from `tools/inchi_read1_reference.py`.

use rdoxide::external::inchi::ichi::{InChI, InChIIsotopicAtom, InChIStereo};
use rdoxide::external::inchi::ichi_io::{INCHI_IOS_TYPE_STRING, InchiIostream};
use rdoxide::external::inchi::ichiread::*;

fn unhex(h: &str) -> Vec<u8> {
    if h == "=" {
        return Vec::new();
    }
    (0..h.len() / 2)
        .map(|i| u8::from_str_radix(&h[2 * i..2 * i + 2], 16).unwrap())
        .collect()
}

fn hex(s: &[u8]) -> String {
    if s.is_empty() {
        "=".into()
    } else {
        s.iter().map(|c| format!("{c:02x}")).collect()
    }
}

#[test]
fn read_helpers_match_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_read1.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let t: Vec<&str> = inp.split(' ').collect();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut got: Vec<String> = Vec::new();
        let p = |k: usize| t[k].parse::<i32>().unwrap();
        match t[0] {
            "1" => {
                let (mob, n0, n1) = (p(1) as usize, p(2), p(3));
                let s = unhex(t[4]);
                let nn = n0.max(n1) as usize;
                let mut pp: [Vec<InChI>; 2] = [
                    vec![InChI::default(); nn + 2],
                    vec![InChI::default(); nn + 2],
                ];
                for i in 0..n1 as usize {
                    pp[1][i].n_total_charge = p(5 + i);
                }
                let r = parse_segment_charge(&s, mob, &mut pp, [n0, n1]);
                got.push(format!("R {r}"));
                let mut q = String::from("Q");
                for k in 0..2 {
                    q.push_str(" |");
                    for i in 0..nn {
                        q.push_str(&format!(" {}", pp[k][i].n_total_charge));
                    }
                }
                got.push(q);
            }
            "2" => {
                let mob = p(1) as usize;
                let s = unhex(t[2]);
                let mut rp = [RemProtons::default(), RemProtons::default()];
                rp[mob].n_num_removed_protons = -777;
                let r = parse_segment_protons(&s, mob, &mut rp);
                got.push(format!("R {r} {}", rp[mob].n_num_removed_protons));
            }
            "3" => {
                let s = unhex(t[1]);
                let (from, pend) = (p(2) as usize, p(3) as usize);
                let mut a = Vec::new();
                let (pos, ret) = parse_segment_read_delimited_numbers(
                    &s,
                    from,
                    pend,
                    &mut a,
                    p(4) as u8,
                    p(5) as u8,
                );
                let mut l = format!("R {ret} {}", pos.map_or(-1, |x| x as i64));
                for v in a {
                    l.push_str(&format!(" {v}"));
                }
                got.push(l);
            }
            "4" => {
                let s = unhex(t[1]);
                let used = p(2);
                let tok = unhex(t[3]);
                let n = s.len() as i32;
                let mut ios = InchiIostream {
                    type_: INCHI_IOS_TYPE_STRING,
                    ..Default::default()
                };
                ios.s.p_str = Some(s);
                ios.s.n_used_length = used.min(n);
                let mut line = SegmLine::default();
                for _ in 0..40 {
                    let c = n_get_inchi_segment(&mut ios, &mut line, &tok);
                    got.push(format!("S {c} {} {}", line.len, hex(&line.str[..line.len])));
                    if c <= 0 {
                        break;
                    }
                }
            }
            "5" => {
                let n = p(1) as usize;
                let mut v: Vec<u16> = (0..n).map(|i| p(2 + i) as u16).collect();
                let r = insertions_sort_at_numb(&mut v, n);
                let mut l = format!("R {r}");
                for x in v {
                    l.push_str(&format!(" {x}"));
                }
                got.push(l);
            }
            "6" => {
                let st = p(1);
                got.push(format!("M {}", get_inchi_state_read_err(st)));
                got.push(get_inchi_err_name(st).to_string());
            }
            "8" => {
                let (a, b) = (unhex(t[1]), unhex(t[2]));
                let (state, first) = (p(3), p(4));
                let mut pp: [Vec<InChI>; 2] = [Vec::new(), Vec::new()];
                let mut pn = [0i32; 2];
                for k in 0..2 {
                    let mob = if (k == 0) == (first == 1) { 1 } else { 0 };
                    let mut na = state;
                    let r = parse_segment_formula(
                        if k == 0 { &a } else { &b },
                        mob,
                        &mut pp,
                        &mut pn,
                        &mut na,
                    );
                    got.push(format!("R{mob} {r} {} {na} {}", pn[mob], pn[1 - mob]));
                    if r >= 0 {
                        for i in 0..pn[mob] as usize {
                            if i >= pp[mob].len() {
                                break;
                            }
                            let c = &pp[mob][i];
                            let atoms: String = (0..=c.n_number_of_atoms as usize)
                                .map(|j| format!("{},", c.n_atom.get(j).copied().unwrap_or(0)))
                                .collect();
                            got.push(format!(
                                "C {i} {} {} {atoms} {} {}",
                                c.sz_hill_formula.as_deref().unwrap_or("NULL"),
                                c.n_number_of_atoms,
                                c.n_total_charge,
                                c.b_deleted
                            ));
                        }
                    }
                    if r < 0 && r != RI_ERR_SYNTAX && r != RI_ERR_PROGR {
                        break;
                    }
                }
            }
            "9" => {
                let (a, b) = (unhex(t[1]), unhex(t[2]));
                let mut abc = p(3);
                let nullfirst = p(4) == 1;
                let mut pp: [Vec<InChI>; 2] = [Vec::new(), Vec::new()];
                let mut pn = [0i32; 2];
                let mut na = 0;
                let mut r = 0;
                if !nullfirst {
                    r = parse_segment_formula(&a, 1, &mut pp, &mut pn, &mut na);
                    got.push(format!("F {r} {}", pn[1]));
                }
                if nullfirst || r >= 0 {
                    let mut nb = -5;
                    let mut v = std::mem::take(&mut pp[1]);
                    r = parse_segment_connections(&b, &mut v, &mut pn[1], &mut abc, &mut nb);
                    got.push(format!("R {r} {} {abc} {nb}", pn[1]));
                    for i in 0..pn[1] as usize {
                        let Some(c) = v.get(i) else { break };
                        let atoms: String = if c.n_atom.is_empty() {
                            "NULL".to_string()
                        } else {
                            (0..=c.n_number_of_atoms as usize)
                                .map(|j| format!("{},", c.n_atom.get(j).copied().unwrap_or(0)))
                                .collect()
                        };
                        let ct: String = if c.n_conn_table.is_empty() {
                            " NULL".to_string()
                        } else {
                            (0..=c.len_conn_table as usize)
                                .map(|j| {
                                    format!(" {}", c.n_conn_table.get(j).copied().unwrap_or(0))
                                })
                                .collect()
                        };
                        got.push(format!(
                            "C {i} {} {atoms} {} |{ct}",
                            c.n_number_of_atoms, c.len_conn_table
                        ));
                    }
                }
            }
            "10" => {
                let mut abc = p(1);
                let hasfixed = p(2) == 1;
                let b: Vec<Vec<u8>> = (3..8).map(|k| unhex(t[k])).collect();
                let mut pp: [Vec<InChI>; 2] = [Vec::new(), Vec::new()];
                let mut pn = [0i32; 2];
                let mut na = IST_MOBILE_H_FORMULA;
                let mut nb = 0;
                let mut r = parse_segment_formula(&b[0], 1, &mut pp, &mut pn, &mut na);
                got.push(format!("F1 {r} {}", pn[1]));
                if r >= 0 {
                    let mut v = std::mem::take(&mut pp[1]);
                    r = parse_segment_connections(&b[1], &mut v, &mut pn[1], &mut abc, &mut nb);
                    pp[1] = v;
                    got.push(format!("C1 {r} {} {abc}", pn[1]));
                }
                if r >= 0 {
                    r = parse_segment_mobile_h(&b[2], 1, &mut pp, pn, &mut abc);
                    got.push(format!("H1 {r} {abc}"));
                }
                if r >= 0 && hasfixed {
                    na = IST_FIXED_H_FORMULA;
                    r = parse_segment_formula(&b[3], 0, &mut pp, &mut pn, &mut na);
                    got.push(format!("F2 {r} {}", pn[0]));
                    if r >= 0 {
                        r = parse_segment_mobile_h(&b[4], 0, &mut pp, pn, &mut abc);
                        got.push(format!("H2 {r} {abc}"));
                    }
                }
                for m in [1usize, 0] {
                    for i in 0..pn[m].max(0) as usize {
                        let Some(c) = pp[m].get(i) else { break };
                        let nat = c.n_number_of_atoms.max(0) as usize;
                        let list = |v: &[i8]| -> String {
                            if v.is_empty() {
                                "NULL".to_string()
                            } else {
                                (0..=nat)
                                    .map(|j| format!("{},", v.get(j).copied().unwrap_or(0)))
                                    .collect()
                            }
                        };
                        let atoms = if c.n_atom.is_empty() {
                            "NULL".to_string()
                        } else {
                            (0..=nat)
                                .map(|j| format!("{},", c.n_atom.get(j).copied().unwrap_or(0)))
                                .collect()
                        };
                        let fixed = c.n_num_h_fixed.as_deref().map_or("NULL".to_string(), list);
                        let taut = if c.n_tautomer.is_empty() {
                            " NULL".to_string()
                        } else {
                            (0..=c.len_tautomer as usize)
                                .map(|j| format!(" {}", c.n_tautomer.get(j).copied().unwrap_or(0)))
                                .collect()
                        };
                        let ct = if c.n_conn_table.is_empty() {
                            " NULL".to_string()
                        } else {
                            (0..=c.len_conn_table as usize)
                                .map(|j| {
                                    format!(" {}", c.n_conn_table.get(j).copied().unwrap_or(0))
                                })
                                .collect()
                        };
                        got.push(format!(
                            "D{m}.{i} {} {} |{atoms} | {} | {fixed} | {} |{taut} | {} |{ct}",
                            c.n_number_of_atoms,
                            c.sz_hill_formula.as_deref().unwrap_or("NULL"),
                            list(&c.n_num_h),
                            c.len_tautomer,
                            c.len_conn_table
                        ));
                    }
                }
            }
            "7" => {
                let mut c = InChI::default();
                let r = n_fill_out_proton_mobile_h(&mut c);
                got.push(
                    format!(
                        "R {r} {} {} {} {} {} {} {} {}",
                        c.b_deleted,
                        c.sz_hill_formula.as_deref().unwrap(),
                        c.n_number_of_atoms,
                        c.n_atom[0],
                        c.n_total_charge,
                        c.n_conn_table[0],
                        c.len_conn_table,
                        c.n_tautomer[0],
                    ) + &format!(" {}", c.len_tautomer),
                );
            }
            _ => unreachable!(),
        }
        assert_eq!(got, want, "case {ci}: {inp}");
    }
}

fn dump_stereo(s: &Option<InChIStereo>) -> String {
    let Some(s) = s else {
        return " NULL".into();
    };
    let mut o = format!(
        " [{} {} {} {} {}{}{}{}{}",
        s.n_number_of_stereo_centers,
        s.n_number_of_stereo_bonds,
        s.n_comp_inv2_abs,
        s.b_trivial_inv,
        u8::from(!s.t_parity.is_empty()),
        u8::from(!s.n_number.is_empty()),
        u8::from(!s.b_parity.is_empty()),
        u8::from(!s.n_bond_atom1.is_empty()),
        u8::from(!s.n_bond_atom2.is_empty()),
    );
    if !s.t_parity.is_empty() && !s.n_number.is_empty() {
        o.push_str(" |");
        for i in 0..s.n_number_of_stereo_centers.max(0) as usize {
            o.push_str(&format!(" {}:{}", s.n_number[i], s.t_parity[i]));
        }
    }
    if !s.b_parity.is_empty() && !s.n_bond_atom1.is_empty() && !s.n_bond_atom2.is_empty() {
        o.push_str(" |");
        for i in 0..s.n_number_of_stereo_bonds.max(0) as usize {
            o.push_str(&format!(
                " {}-{}:{}",
                s.n_bond_atom1[i], s.n_bond_atom2[i], s.b_parity[i]
            ));
        }
    }
    o.push(']');
    o
}

fn dump_inchi(c: &InChI) -> String {
    let mut o = format!(
        "{} {} {}",
        c.n_number_of_atoms,
        c.n_number_of_isotopic_atoms,
        u8::from(!c.isotopic_atom.is_empty())
    );
    if !c.isotopic_atom.is_empty() {
        for a in c
            .isotopic_atom
            .iter()
            .take(c.n_number_of_isotopic_atoms.max(0) as usize)
        {
            o.push_str(&format!(
                " ({} {} {} {} {})",
                a.n_atom_number, a.n_iso_difference, a.n_num_h, a.n_num_d, a.n_num_t
            ));
        }
    }
    o.push_str(&dump_stereo(&c.stereo));
    o.push_str(&dump_stereo(&c.stereo_isotopic));
    o
}

fn mk_stereo(t: &[&str], k: &mut usize) -> Option<InChIStereo> {
    let n = |k: &mut usize| {
        *k += 1;
        t[*k - 1].parse::<i32>().unwrap()
    };
    if n(k) == 0 {
        return None;
    }
    let (nsc, nsb, inv, triv, tl, bl) = (n(k), n(k), n(k), n(k), n(k), n(k));
    let mut s = InChIStereo {
        n_number_of_stereo_centers: nsc,
        n_number_of_stereo_bonds: nsb,
        n_comp_inv2_abs: inv,
        b_trivial_inv: triv,
        ..Default::default()
    };
    if tl != 0 {
        s.t_parity = (0..=nsc)
            .map(|i| if i < nsc { ((i * 3 + 1) % 5) as i8 } else { 0 })
            .collect();
        s.n_number = (0..=nsc)
            .map(|i| if i < nsc { (i * 7 + 2) as u16 } else { 0 })
            .collect();
    }
    if bl != 0 {
        s.b_parity = (0..=nsb)
            .map(|i| if i < nsb { (i % 4 + 1) as i8 } else { 0 })
            .collect();
        s.n_bond_atom1 = (0..=nsb)
            .map(|i| if i < nsb { (i + 10) as u16 } else { 0 })
            .collect();
        s.n_bond_atom2 = (0..=nsb)
            .map(|i| if i < nsb { (i + 3) as u16 } else { 0 })
            .collect();
    }
    Some(s)
}

fn mk_inchi(t: &[&str], k: &mut usize) -> InChI {
    let mut n = || {
        *k += 1;
        t[*k - 1].parse::<i32>().unwrap()
    };
    let (na, ni, ip) = (n(), n(), n());
    let mut c = InChI {
        n_number_of_atoms: na,
        n_number_of_isotopic_atoms: ni,
        ..Default::default()
    };
    if ip != 0 {
        c.isotopic_atom = (0..=ni)
            .map(|i| {
                if i < ni {
                    InChIIsotopicAtom {
                        n_atom_number: (i + 1) as u16,
                        n_iso_difference: (i - 1) as i16,
                        n_num_h: i as i16,
                        n_num_d: (2 * i) as i16,
                        n_num_t: 3,
                    }
                } else {
                    InChIIsotopicAtom::default()
                }
            })
            .collect();
    }
    c.stereo = mk_stereo(t, k);
    c.stereo_isotopic = mk_stereo(t, k);
    c
}

#[test]
fn stereo_layers_match_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_read2.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let t: Vec<&str> = inp.split(' ').collect();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut got: Vec<String> = Vec::new();
        let p = |k: usize| t[k].parse::<i32>().unwrap();
        match t[0] {
            "11" => {
                let mut abc = p(1);
                let nc = [p(2), p(3)];
                let mut k = 4;
                let mut pp: [Vec<InChI>; 2] = [Vec::new(), Vec::new()];
                for m in 0..2 {
                    if nc[m] > 0 {
                        pp[m] = vec![InChI::default(); nc[m] as usize + 2];
                        for i in 0..nc[m] as usize {
                            pp[m][i].n_number_of_atoms = p(k);
                            k += 1;
                        }
                    }
                }
                let mut sv = [[0i32; 2]; 2];
                let mut rp = [RemProtons::default(), RemProtons::default()];
                let nseg = p(k);
                k += 1;
                for _ in 0..nseg {
                    let state = p(k);
                    let s = unhex(t[k + 1]);
                    k += 2;
                    let mob = usize::from(state < IST_FIXED_H_LAYER_FORK);
                    let r = match state {
                        IST_MOBILE_H_SP2 | IST_FIXED_H_SP2 | IST_MOBILE_H_ISO_SP2
                        | IST_FIXED_H_ISO_SP2 => {
                            parse_segment_sp2(&s, mob, &mut pp, nc, state, &mut abc)
                        }
                        IST_MOBILE_H_SP3 | IST_FIXED_H_SP3 | IST_MOBILE_H_ISO_SP3
                        | IST_FIXED_H_ISO_SP3 => {
                            parse_segment_sp3(&s, mob, &mut pp, nc, state, &mut abc)
                        }
                        IST_MOBILE_H_SP3_M
                        | IST_FIXED_H_SP3_M
                        | IST_MOBILE_H_ISO_SP3_M
                        | IST_FIXED_H_ISO_SP3_M => parse_segment_sp3m(&s, mob, &mut pp, nc, state),
                        IST_MOBILE_H_ISO_ATOMS | IST_FIXED_H_ISO_ATOMS => {
                            parse_segment_iso_atoms(&s, mob, &mut pp, nc, state, &mut abc)
                        }
                        IST_MOBILE_H_ISO_EXCH_H => {
                            parse_segment_iso_exchg_h(&s, mob, &mut rp, state, &mut abc)
                        }
                        IST_FIXED_H_PERMUTATION | IST_FIXED_H_ISO_PERMUTATION => {
                            parse_segment_perm(&s, mob, &mut pp, nc, state, &mut abc)
                        }
                        _ => parse_segment_sp3s(&s, mob, &mut pp, &mut sv, nc, state),
                    };
                    got.push(format!(
                        "R {state} {r} {abc} {} {} {}",
                        sv[0][0],
                        sv[0][1],
                        sv[1][0] * 10 + sv[1][1]
                    ));
                    if r < 0 {
                        break;
                    }
                }
                let h = rp[1].n_num_removed_isotopic_h;
                got.push(format!("P {} {} {}", h[0], h[1], h[2]));
                for m in 0..2 {
                    for i in 0..nc[m] as usize {
                        got.push(format!("C{m}.{i} {}", dump_inchi(&pp[m][i])));
                    }
                }
            }
            "12" => {
                let (ty, bto, bfrom, hasfrom) = (p(1), p(2), p(3), p(4) != 0);
                let mut k = 5;
                let mut to = mk_inchi(&t, &mut k);
                let from = hasfrom.then(|| mk_inchi(&t, &mut k));
                let r = copy_segment(&mut to, from.as_ref(), ty, bto, bfrom);
                got.push(format!("R {r}"));
                got.push(format!("C {}", dump_inchi(&to)));
            }
            _ => unreachable!(),
        }
        assert_eq!(got, want, "case {ci}: {inp}");
    }
}
