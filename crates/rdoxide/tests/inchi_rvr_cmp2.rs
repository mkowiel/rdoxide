#![allow(clippy::field_reassign_with_default)]
//! Differential test of `FillOutCMP2FHINCHI` / `FillOutCMP2MHINCHI`
//! (`ichirvr4.rs`) against the C library on random `StrFromINChI` data, with
//! the fixture from `tools/inchi_rvr_cmp2_reference.py`.

use rdoxide::external::inchi::ichi::InChI;
use rdoxide::external::inchi::ichimake::InpAtomData;
use rdoxide::external::inchi::ichirvr2::StrFromInChI;
use rdoxide::external::inchi::ichirvr4::{
    Cmp2FhInchi, Cmp2MhInchi, fill_out_cmp2fh_inchi, fill_out_cmp2mh_inchi,
};
use rdoxide::external::inchi::ichirvrs::{AllTcGroups, ValAt};
use rdoxide::external::inchi::ichitaut::{TGroup, TGroupInfo};
use rdoxide::external::inchi::inpdef::InpAtom;

struct Rd<'a>(std::slice::Iter<'a, i32>);
impl Rd<'_> {
    fn n(&mut self) -> i32 {
        *self.0.next().unwrap()
    }
}

fn read_ti(r: &mut Rd) -> TGroupInfo {
    let nt = r.n() as usize;
    let mut groups = Vec::new();
    for _ in 0..nt {
        let mut g = TGroup::default();
        g.num[0] = r.n() as u16;
        g.num[1] = r.n() as u16;
        g.n_num_endpoints = r.n() as u16;
        g.n_first_endpoint_at_no_pos = r.n() as u16;
        groups.push(g);
    }
    let ne = r.n() as usize;
    let eps: Vec<u16> = (0..ne).map(|_| r.n() as u16).collect();
    TGroupInfo {
        t_group: Some(groups),
        num_t_groups: nt,
        n_num_endpoints: ne as i32,
        n_endpoint_atom_number: eps,
        ..Default::default()
    }
}

fn read_inchi(r: &mut Rd, with_fixed: bool) -> Option<InChI> {
    if r.n() == 0 {
        return None;
    }
    let mut c = InChI {
        n_number_of_atoms: r.n(),
        b_deleted: r.n(),
        ..Default::default()
    };
    if r.n() != 0 {
        let n = r.n();
        c.n_num_h = (0..n).map(|_| r.n() as i8).collect();
    }
    if with_fixed && r.n() != 0 {
        let n = r.n();
        c.n_num_h_fixed = Some((0..n).map(|_| r.n() as i8).collect());
    }
    c.n_total_charge = r.n();
    Some(c)
}

#[test]
fn cmp2_matches_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_rvr_cmp2.json")).unwrap();
    let (mut n_diff, mut n_err) = (0, 0);
    for case in data.as_array().unwrap() {
        let input = case["in"].as_str().unwrap();
        let want: Vec<&str> = case["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_str().unwrap())
            .collect();
        let vals: Vec<i32> = input
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut r = Rd(vals.iter());
        let mode = r.n();
        let n = r.n() as usize;
        let mut ps = StrFromInChI {
            num_atoms: n as i32,
            ..Default::default()
        };
        ps.ti = read_ti(&mut r);
        ps.one_ti = read_ti(&mut r);
        ps.one_ti.tni.n_num_removed_protons = r.n();
        ps.n_num_removed_protons_mobh_inchi = r.n();
        ps.endpoint = Some((0..n).map(|_| r.n() as u16).collect());
        ps.fixed_h = Some((0..n).map(|_| r.n() as i8).collect());
        let mut va = vec![ValAt::default(); n + 1];
        for v in va.iter_mut().take(n) {
            v.c_num_valence_electrons = r.n() as i8;
            v.c_periodic_row_number = r.n() as i8;
            v.c_metal = r.n() as i8;
        }
        let mut at2 = vec![InpAtom::default(); n + 1];
        for a in at2.iter_mut().take(n) {
            a.num_h = r.n() as i8;
            a.charge = r.n() as i8;
            a.valence = r.n() as _;
            a.chem_bonds_valence = r.n() as _;
            a.endpoint = r.n() as u16;
        }
        let tcg = AllTcGroups {
            tgroup_charge: r.n(),
            ..Default::default()
        };
        let pin0 = read_inchi(&mut r, false);
        let pin1 = read_inchi(&mut r, false);
        ps.p_one_inchi[0] = read_inchi(&mut r, true);
        ps.p_one_inchi[1] = read_inchi(&mut r, true);
        for k in 0..2 {
            if r.n() != 0 {
                let mut d = InpAtomData::default();
                d.at = vec![InpAtom::default(); n + 1];
                for a in d.at.iter_mut().take(n) {
                    a.endpoint = r.n() as u16;
                    a.charge = r.n() as i8;
                }
                ps.p_one_norm_data[k] = Some(d);
            }
        }
        let atno2canon: Vec<u16> = (0..n).map(|_| r.n() as u16).collect();
        let mut canon2atno = vec![0u16; n + 1];
        for (i, &c) in atno2canon.iter().enumerate() {
            canon2atno[c as usize] = i as u16;
        }
        ps.n_atno2canon[0] = atno2canon;
        ps.n_canon2atno[0] = canon2atno;
        let pin = [pin0.as_ref(), pin1.as_ref()];
        let mut got = Vec::new();
        if mode == 0 {
            if ps.p_one_inchi[0].is_none() {
                panic!("fixture case needs pOneINChI[0]");
            }
            let mut c = Cmp2FhInchi::default();
            let ret = fill_out_cmp2fh_inchi(&ps, &at2, &va, pin, &mut c);
            got.push(format!(
                "{ret} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                c.c2at.len(),
                c.n_num_rem_h_inchi,
                c.n_num_rem_h_revrs,
                c.n_num_tg_inchi,
                c.n_num_tg_revrs,
                c.n_num_endp_inchi,
                c.n_num_endp_revrs,
                c.n_num_tg_diff_minus,
                c.n_num_tg_diff_h,
                c.n_num_tg_m_inchi,
                c.n_num_tg_h_inchi,
                c.n_num_tg_m_revrs,
                c.n_num_tg_h_revrs,
                c.n_charge_fix_h_inchi,
                c.n_charge_mob_h_inchi,
                c.n_charge_fix_h_revrs,
                c.n_charge_mob_h_revrs,
                c.n_charge_fix_h_revrs_non_metal,
                c.n_charge_mob_h_revrs_non_metal,
                i32::from(c.b_fixed_h_layer_exists_revrs),
                i32::from(c.b_has_difference),
                c.n_num_diff_mob_h
            ));
            for a in &c.c2at {
                got.push(format!(
                    "a {} {} {} {} {} {} {} {} {} {} {}",
                    a.endpt_inchi,
                    a.endpt_revrs,
                    a.atom_number,
                    a.n_val_electr,
                    a.n_period_num,
                    a.n_fix_h_inchi,
                    a.n_fix_h_revrs,
                    a.n_mob_h_inchi,
                    a.n_mob_h_revrs,
                    a.n_num_h_revrs,
                    a.n_at_charge_revrs
                ));
            }
            n_diff += i32::from(c.b_has_difference);
            n_err += i32::from(ret != 0);
        } else {
            let mut c = Cmp2MhInchi::default();
            let ret = fill_out_cmp2mh_inchi(&ps, &tcg, &at2, &va, pin, &mut c);
            let j = |a: &[i16; 10]| {
                a.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(" ")
            };
            got.push(format!(
                "{ret} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                c.c2at.len(),
                c.n_num_rem_h_inchi,
                c.n_num_rem_h_revrs,
                c.n_num_tg_inchi,
                c.n_num_tg_revrs,
                c.n_num_endp_inchi,
                c.n_num_endp_revrs,
                c.n_num_tg_diff_minus,
                c.n_num_tg_diff_h,
                c.n_num_tg_m_inchi,
                c.n_num_tg_h_inchi,
                c.n_num_tg_o_inchi,
                c.n_num_tg_n_inchi,
                c.n_num_tg_m_revrs,
                c.n_num_tg_h_revrs,
                c.n_num_tg_o_revrs,
                c.n_num_tg_n_revrs,
                j(&c.revrs),
                j(&c.inchi),
                c.n_charge_mob_h_inchi,
                c.n_charge_mob_h_revrs,
                c.n_charge_mob_h_revrs_non_metal,
                i32::from(c.b_fixed_h_layer_exists_revrs),
                i32::from(c.b_has_difference)
            ));
            let last = got.last_mut().unwrap();
            last.push_str(&format!(" {}", c.n_num_diff_mob_h));
            for a in &c.c2at {
                got.push(format!(
                    "a {} {} {} {} {} {} {} {} {}",
                    a.endpt_inchi,
                    a.endpt_revrs,
                    a.atom_number,
                    a.n_val_electr,
                    a.n_period_num,
                    a.n_mob_h_inchi,
                    a.n_mob_h_revrs,
                    a.n_num_h_revrs,
                    a.n_at_charge_revrs
                ));
            }
            n_diff += i32::from(c.b_has_difference);
            n_err += i32::from(ret != 0);
        }
        assert_eq!(got, want, "case: {input}");
    }
    eprintln!("cases with difference: {n_diff}, with error: {n_err}");
}
