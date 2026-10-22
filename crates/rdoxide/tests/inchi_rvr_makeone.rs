//! Differential test of `MakeOneInChIOutOfStrFromINChI` (`ichirvr1.rs`)
//! against the C library on random small molecules used as restored
//! structures, with the fixture from `tools/inchi_rvr_makeone_reference.py`.

use rdoxide::external::inchi::ichi::InChI;
use rdoxide::external::inchi::ichiparm::InputParms;
use rdoxide::external::inchi::ichirvr1::{
    make_inchi_out_of_str_from_inchi2, make_one_inchi_out_of_str_from_inchi,
};
use rdoxide::external::inchi::ichirvr2::StrFromInChI;
use rdoxide::external::inchi::inpdef::InpAtom;
use rdoxide::external::inchi::runichi3::StructData;

fn dump_inchi(k: usize, p: &Option<InChI>) -> Vec<String> {
    let Some(p) = p else {
        return vec![format!("I{k} -")];
    };
    let na = p.n_number_of_atoms as usize;
    let join = |v: &[i32]| v.iter().map(|x| format!(" {x}")).collect::<String>();
    vec![
        format!(
            "I{k} {} {} {} {} {}",
            p.n_error_code,
            p.n_flags,
            p.n_total_charge,
            p.n_number_of_atoms,
            p.sz_hill_formula.as_deref().unwrap_or("NULL")
        ),
        format!(
            "IA{k}{}",
            join(
                &p.n_atom[..na]
                    .iter()
                    .map(|&x| i32::from(x))
                    .collect::<Vec<_>>()
            )
        ),
        format!(
            "IC{k} {}{}",
            p.len_conn_table,
            join(
                &p.n_conn_table[..p.len_conn_table as usize]
                    .iter()
                    .map(|&x| i32::from(x))
                    .collect::<Vec<_>>()
            )
        ),
        format!(
            "IT{k} {}{}",
            p.len_tautomer,
            join(
                &p.n_tautomer[..p.len_tautomer as usize]
                    .iter()
                    .map(|&x| i32::from(x))
                    .collect::<Vec<_>>()
            )
        ),
        format!(
            "IH{k}{}",
            (0..na)
                .map(|i| format!(
                    " {}/{}",
                    p.n_num_h[i],
                    p.n_num_h_fixed.as_ref().map_or(-1, |f| i32::from(f[i]))
                ))
                .collect::<String>()
        ),
    ]
}

#[test]
fn make_one_inchi_matches_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_rvr_makeone.json")).unwrap();
    let mut n_ok = 0;
    let (mut n_whole, mut n_diff) = (0, 0);
    for case in data.as_array().unwrap() {
        let input = case["in"].as_str().unwrap();
        let want: Vec<&str> = case["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_str().unwrap())
            .collect();
        let mut t = input.split_whitespace();
        let int = |t: &mut std::str::SplitWhitespace| -> i64 { t.next().unwrap().parse().unwrap() };
        let n = int(&mut t) as usize;
        let bmob = int(&mut t) as i32;
        let imob = int(&mut t) as i32;
        let mode = int(&mut t) as u64;
        let flags = int(&mut t) as u64;
        let loose = int(&mut t) as i32;
        let mut at2 = vec![InpAtom::default(); n + 8];
        for (i, a) in at2.iter_mut().take(n).enumerate() {
            let v = int(&mut t) as usize;
            a.valence = v as i8;
            a.chem_bonds_valence = int(&mut t) as i8;
            let el = int(&mut t);
            a.el_number = el as u8;
            let name = match el {
                6 => "C",
                7 => "N",
                8 => "O",
                15 => "P",
                16 => "S",
                17 => "Cl",
                _ => "Xx",
            };
            a.elname[..name.len()].copy_from_slice(name.as_bytes());
            a.radical = int(&mut t) as i8;
            a.charge = int(&mut t) as i8;
            a.num_h = int(&mut t) as i8;
            a.num_iso_h = [int(&mut t) as i8, int(&mut t) as i8, int(&mut t) as i8];
            a.iso_atw_diff = int(&mut t) as i8;
            a.x = t.next().unwrap().parse().unwrap();
            a.y = t.next().unwrap().parse().unwrap();
            a.z = t.next().unwrap().parse().unwrap();
            a.orig_at_number = i as u16 + 1;
            for j in 0..v {
                a.neighbor[j] = int(&mut t) as u16;
                a.bond_type[j] = int(&mut t) as u8;
            }
        }
        let ip = InputParms {
            b_loose_tsa_check: loose,
            n_mode: mode,
            b_taut_flags: flags,
            ..Default::default()
        };
        let mut ps = StrFromInChI {
            num_atoms: n as i32,
            b_mobile_h: bmob,
            i_mobile_h: imob,
            ..Default::default()
        };
        let mut at3 = Vec::new();
        let ret = make_one_inchi_out_of_str_from_inchi(&ip, &mut ps, &at2, &mut at3, 0);
        let mut got = vec![format!(
            "R {ret} {} {} {} {}",
            ps.n_charge_revrs,
            ps.n_one_inchi_b_mobile_h,
            ps.n_num_removed_protons,
            ps.one_ti.num_t_groups
        )];
        let mut l = format!("TI {}", ps.one_ti.n_num_endpoints);
        if let Some(g) = &ps.one_ti.t_group {
            for x in g.iter().take(ps.one_ti.num_t_groups) {
                l += &format!(
                    " ({} {} {} {} {})",
                    x.num[0],
                    x.num[1],
                    x.n_group_number,
                    x.n_num_endpoints,
                    x.n_first_endpoint_at_no_pos
                );
            }
        }
        l += " |";
        for e in ps
            .one_ti
            .n_endpoint_atom_number
            .iter()
            .take(ps.one_ti.n_num_endpoints as usize)
        {
            l += &format!(" {e}");
        }
        got.push(l);
        for k in 0..2 {
            got.extend(dump_inchi(k, &ps.p_one_inchi[k]));
            match &ps.p_one_norm_data[k] {
                None => got.push(format!("D{k} -")),
                Some(d) => {
                    let mut l = format!(
                        "D{k} {} {} {} {}",
                        d.num_at, d.num_removed_h, d.b_tautomeric, d.b_taut_preprocessed
                    );
                    for a in d.at.iter().take(d.num_at as usize) {
                        l += &format!(
                            " {},{},{},{},{}",
                            a.charge, a.num_h, a.chem_bonds_valence, a.endpoint, a.valence
                        );
                    }
                    got.push(l);
                }
            }
        }
        assert_eq!(got, want, "case: {input}");
        n_ok += i32::from(ret == 0);
        // the whole-structure InChI of the same restored structure must agree
        // with the one-component InChI
        if ret == 0 && ps.n_one_inchi_b_mobile_h == 1 {
            let mut ps2 = StrFromInChI {
                num_atoms: n as i32,
                b_mobile_h: bmob,
                i_mobile_h: imob,
                ..Default::default()
            };
            let mut at2c = at2.clone();
            let mut sd = StructData::default();
            let r = make_inchi_out_of_str_from_inchi2(&ip, &mut sd, &mut ps2, &mut at2c, 0, 1);
            if r >= 0
                && ps2.rev_inchi.num_components[0] == 1
                && ps2.rev_inchi_component_exists(0, 1, 0)
            {
                let rev = ps2.rev_inchi.p_inchi[0][0][1].as_ref().unwrap();
                let one = ps.p_one_inchi[0].as_ref().unwrap();
                n_whole += 1;
                if rev.sz_hill_formula != one.sz_hill_formula
                    || rev.n_total_charge != one.n_total_charge
                {
                    n_diff += 1;
                }
            }
        }
    }
    eprintln!("cases with result: {n_ok}, whole-structure InChI: {n_whole}, differ: {n_diff}");
    assert!(n_whole > 0);
    assert_eq!(n_diff, 0);
}
