//! Differential tests of the InChI string reader (`ichiread.rs`) against the
//! C library, with fixtures from `tools/inchi_read1_reference.py`.

use rdoxide::external::inchi::ichi::InChI;
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
