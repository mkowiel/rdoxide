//! `external::inchi::ichisort` against the output of the InChI C code
//! (`tools/inchi_sort_reference.py`): the order of equal elements after the
//! unstable `inchi_qsort`, and the swaps of `insertions_sort`.

use rdoxide::external::inchi::ichisort::{
    iisort, inchi_qsort, insertions_sort, insertions_sort_at_rank,
    insertions_sort_neigh_list_at_numbers,
};

#[test]
fn sorts_match_the_c_code() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_sort.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() > 300);
    for case in cases {
        let keys: Vec<i32> = case["keys"]
            .as_array()
            .unwrap()
            .iter()
            .map(|k| k.as_i64().unwrap() as i32)
            .collect();
        let expected: Vec<usize> = case["order"]
            .as_array()
            .unwrap()
            .iter()
            .map(|k| k.as_u64().unwrap() as usize)
            .collect();
        let mut v: Vec<(i32, usize)> = keys.iter().copied().zip(0..).collect();
        let swaps = match case["mode"].as_str().unwrap() {
            "qsort" => {
                inchi_qsort(&mut v, |a, b| a.0 - b.0);
                0
            }
            _ => insertions_sort(&mut v, |a, b| a.0 - b.0),
        };
        let order: Vec<usize> = v.iter().map(|e| e.1).collect();
        assert_eq!(order, expected, "keys {keys:?}");
        assert_eq!(
            swaps as i64,
            case["swaps"].as_i64().unwrap(),
            "keys {keys:?}"
        );
    }
}

#[test]
fn small_sorts() {
    let mut a = [3u16, 1, 2];
    assert_eq!(insertions_sort_at_rank(&mut a), 2);
    assert_eq!(a, [1, 2, 3]);
    // NEIGH_LIST layout: the count, then the neighbors
    let ranks = [0u16, 30, 10, 20];
    let mut nl = [3u16, 1, 2, 3];
    insertions_sort_neigh_list_at_numbers(&mut nl, &ranks);
    assert_eq!(nl, [3, 2, 3, 1]);
    let mut l = [5, 2, 9, 1];
    iisort(&mut l);
    assert_eq!(l, [1, 2, 5, 9]);
}

#[test]
fn error_messages_match_the_c_code() {
    use rdoxide::external::inchi::ichierr::{add_error_message, err_msg};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_err.json")).unwrap();
    for (code, msg) in data["codes"]
        .as_array()
        .unwrap()
        .iter()
        .zip(data["messages"].as_array().unwrap())
    {
        let code = code.as_i64().unwrap() as i32;
        assert_eq!(err_msg(code), msg.as_str().unwrap(), "code {code}");
    }
    let cases = data["cases"].as_array().unwrap();
    assert!(cases.len() >= 300);
    for case in cases {
        let mut buf = String::new();
        let outs = case["out"].as_array().unwrap();
        for (m, out) in case["in"].as_array().unwrap().iter().zip(outs) {
            let r = add_error_message(&mut buf, m.as_str().unwrap());
            assert_eq!(i64::from(r), out[0].as_i64().unwrap(), "{case}");
            assert_eq!(buf, out[1].as_str().unwrap(), "{case}");
        }
    }
}

#[test]
fn element_helpers_match_the_c_code() {
    use rdoxide::external::inchi::incomdef::NUM_H_ISOTOPES;
    use rdoxide::external::inchi::util::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_util.json")).unwrap();
    let cases = data["cases"].as_array().unwrap();
    assert!(cases.len() > 7000);
    let ints = |s: &str| -> Vec<i32> { s.split(' ').map(|v| v.parse().unwrap()).collect() };
    for case in cases {
        let (input, expected) = (case["in"].as_str().unwrap(), case["out"].as_str().unwrap());
        let (op, rest) = input.split_once(' ').unwrap();
        let got = match op {
            "S" => {
                let n: i32 = rest.parse().unwrap();
                let (r, s) = get_element_chemical_symbol(n);
                let (r2, t) = get_element_or_pseudoelement_symbol(n);
                format!("{r} {s} {r2} {t}")
            }
            "P" => format!(
                "{} {} {} 0",
                el_number_in_internal_ref_table(rest),
                get_periodic_table_number(rest),
                get_atomic_mass(rest)
            ),
            "E" => {
                let n: i32 = rest.parse().unwrap();
                format!(
                    "{} {} {} {}",
                    get_atomic_mass_from_elnum(n),
                    i32::from(if_skip_add_h(n)),
                    get_el_type(n),
                    i32::from(is_el_a_metal(n))
                )
            }
            "V" => {
                let a = ints(rest);
                get_el_valence(a[0], a[1], a[2] as usize).to_string()
            }
            "U" => {
                let a = ints(rest);
                format!(
                    "{} {} {}",
                    get_unusual_el_valence(a[0], a[1], a[2], a[3], a[5], a[6]),
                    needed_unusual_el_valence(a[0], a[1], a[2], a[3], a[4], a[5], a[6]),
                    detect_unusual_el_valence(a[0], a[1], a[2], a[3], a[5], a[6])
                )
            }
            "C" => {
                let mut s = rest.as_bytes().to_vec();
                let (rad, ch) = extract_charges_and_radicals(&mut s);
                format!(
                    "{} {rad} {ch} [{}]",
                    i32::from(rad != 0 || ch != 0),
                    String::from_utf8(s).unwrap()
                )
            }
            "H" => {
                let mut s = rest.as_bytes().to_vec();
                let mut iso = [0i8; NUM_H_ISOTOPES];
                let h = extract_h_atoms(&mut s, &mut iso);
                format!(
                    "{h} {} {} {} [{}]",
                    iso[0],
                    iso[1],
                    iso[2],
                    String::from_utf8(s).unwrap()
                )
            }
            "N" => {
                let (el, nums) = rest.split_once(' ').unwrap();
                let a = ints(nums);
                let iso = [a[1] as i8, a[2] as i8, a[3] as i8];
                get_num_h(
                    el,
                    a[0],
                    Some(&iso),
                    a[4],
                    a[5],
                    a[6],
                    a[7],
                    a[8] != 0,
                    a[9] != 0,
                    a[10] != 0,
                )
                .to_string()
            }
            _ => panic!("{input}"),
        };
        assert_eq!(got, expected, "{input}");
    }
}

#[test]
fn ring_search_matches_the_c_code() {
    use rdoxide::external::inchi::ichiring::{
        Queue, is_atom_in_3memb_ring, is_bond_in_nmax_memb_ring,
    };
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_ring.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() >= 400);
    let ints = |v: &serde_json::Value| -> Vec<i64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect()
    };
    for case in cases {
        let adj = case["adj"].as_array().unwrap();
        let n = adj.len();
        let mut atoms = vec![InpAtom::default(); n];
        for (i, a) in adj.iter().enumerate() {
            let nb = ints(a);
            atoms[i].valence = nb.len() as i8;
            for (j, x) in nb.iter().enumerate() {
                atoms[i].neighbor[j] = *x as u16;
            }
            atoms[i].n_ring_system = case["rs"][i].as_u64().unwrap() as u16;
            atoms[i].n_num_at_in_ring_system = case["nar"][i].as_u64().unwrap() as u16;
        }
        let mut q = Queue::create(n as i32 + 1).unwrap();
        let mut lev = vec![0u16; n];
        let mut src = vec![0i8; n];
        let expected = ints(&case["ring"]);
        for (k, query) in case["q"].as_array().unwrap().iter().enumerate() {
            let p = ints(query);
            let r = is_bond_in_nmax_memb_ring(
                &atoms,
                p[0] as usize,
                p[1] as usize,
                &mut q,
                &mut lev,
                &mut src,
                p[2] as u16,
            );
            assert_eq!(i64::from(r), expected[k], "query {p:?} in {adj:?}");
        }
        assert!(lev.iter().all(|&l| l == 0) && src.iter().all(|&s| s == 0));
        let tri = ints(&case["tri"]);
        for i in 0..n {
            assert_eq!(
                is_atom_in_3memb_ring(&atoms, i) as i64,
                tri[i],
                "atom {i} in {adj:?}"
            );
        }
    }
}

#[test]
fn atom_helpers_match_the_c_code() {
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::util::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_atoms.json")).unwrap();
    for (el, g) in data["groups"].as_array().unwrap().iter().enumerate() {
        assert_eq!(
            ion_el_group(el as u8) as i64,
            g[0].as_i64().unwrap(),
            "el {el}"
        );
        assert_eq!(
            get_endpoint_valence(el as u8) as i64,
            g[1].as_i64().unwrap(),
            "el {el}"
        );
    }
    let cases = data["cases"].as_array().unwrap();
    assert!(cases.len() >= 500);
    for case in cases {
        let atoms_json = case.as_array().unwrap();
        let mut at = vec![InpAtom::default(); atoms_json.len()];
        for (i, a) in atoms_json.iter().enumerate() {
            at[i].el_number = a["el"].as_u64().unwrap() as u8;
            at[i].charge = a["ch"].as_i64().unwrap() as i8;
            at[i].chem_bonds_valence = a["cbv"].as_i64().unwrap() as i8;
            at[i].num_h = a["h"][0].as_i64().unwrap() as i8;
            for k in 0..3 {
                at[i].num_iso_h[k] = a["h"][k + 1].as_i64().unwrap() as i8;
            }
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            for (j, e) in nb.iter().enumerate() {
                at[i].neighbor[j] = e[0].as_u64().unwrap() as u16;
                at[i].bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
        }
        for (i, a) in atoms_json.iter().enumerate() {
            let exp: Vec<i64> = a["out"]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap())
                .collect();
            let n0 = if at[i].valence > 0 {
                i32::from(at[i].neighbor[0])
            } else {
                -1
            };
            let n1 = if at[i].valence > 1 {
                i32::from(at[i].neighbor[1])
            } else {
                -1
            };
            let got = [
                i64::from(n_bonds_val_to_metal(&at, i)),
                i64::from(num_of_h(&at, i)),
                i64::from(n_no_metal_num_bonds(&at, i)),
                i64::from(n_no_metal_bonds_valence(&at, i)),
                i64::from(n_no_metal_neigh_index(&at, i)),
                i64::from(n_no_metal_other_neigh_index(&at, i, n0)),
                i64::from(n_no_metal_other_neigh_index2(&at, i, n0, n1)),
                if n0 >= 0 {
                    i64::from(has_other_ion_neigh(&at, i, n0 as usize))
                } else {
                    -1
                },
                if n0 >= 0 {
                    i64::from(has_other_ion_in_sphere_2(&mut at, i, n0 as usize))
                } else {
                    -1
                },
            ];
            assert_eq!(got.to_vec(), exp, "atom {i}: {a}");
        }
    }
}

fn unhex(s: &str) -> Vec<u8> {
    let h = &s[1..];
    (0..h.len() / 2)
        .map(|i| u8::from_str_radix(&h[2 * i..2 * i + 2], 16).unwrap())
        .collect()
}

fn hex(b: &[u8]) -> String {
    let mut s = String::from("x");
    for c in b {
        s += &format!("{c:02x}");
    }
    s
}

#[test]
fn string_helpers_match_the_c_code() {
    use rdoxide::external::inchi::util::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_strings.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() >= 600);
    for c in cases {
        let s1 = unhex(c["s1"].as_str().unwrap());
        let s2 = unhex(c["s2"].as_str().unwrap());
        let slen = c["slen"].as_u64().unwrap() as usize;
        let maxlen = c["maxlen"].as_i64().unwrap() as i32;
        let out: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let ctx = format!("{:?}", String::from_utf8_lossy(&s1));
        let mut b = s1.clone();
        let len = normalize_string(&mut b);
        assert_eq!(format!("{} {len}", hex(&b)), out[0], "normalize {ctx}");
        let mut b = s1.clone();
        let n = dotify_non_printable_chars(&mut b);
        assert_eq!(format!("{} {n}", hex(&b)), out[1], "dotify {ctx}");
        let mut b = s1.clone();
        let len = lrtrim(&mut b);
        assert_eq!(format!("{} {len}", hex(&b)), out[2], "lrtrim {ctx}");
        let mut b = s1.clone();
        remove_trailing_spaces(&mut b);
        assert_eq!(hex(&b), out[3], "trailing {ctx}");
        let mut b = s1.clone();
        remove_one_lf(&mut b);
        assert_eq!(hex(&b), out[4], "lf {ctx}");
        let mut pos = Some(0);
        let mut field = Vec::new();
        let n = read_upto_delim(&s1, &mut pos, &mut field, maxlen, b", \t");
        let fld = if n > 0 { hex(&field) } else { "x".to_string() };
        let p = pos.map_or(-1, |p| p as i64);
        assert_eq!(format!("{n} {fld} {p}"), out[5], "read_upto_delim {ctx}");
        let st = String::from_utf8_lossy(&s1).to_string();
        let ascii = s1.is_ascii() && s2.is_ascii();
        if ascii {
            let o =
                extract_inchi_substring(&st, slen).map_or("N".to_string(), |s| hex(s.as_bytes()));
            assert_eq!(o, out[6], "extract_inchi {ctx}");
            let o =
                extract_auxinfo_substring(&st, slen).map_or("N".to_string(), |s| hex(s.as_bytes()));
            assert_eq!(o, out[7], "extract_aux {ctx}");
            let a = inchi_stricmp(&st, &String::from_utf8_lossy(&s2)).signum();
            let m = inchi_memicmp(&s1, &s2, s1.len().min(s2.len())).signum();
            assert_eq!(format!("{a} {m}"), out[8], "icmp {ctx}");
        }
        let mut t = vec![b'z'; 64];
        t[63] = 0;
        mystrncpy(&mut t, &s1, maxlen as usize);
        assert_eq!(hex(&t), out[9], "mystrncpy {ctx}");
        let bi = c["bi"].as_i64().unwrap();
        let arr = |k: &str| -> [i32; 3] {
            let a = c[k].as_array().unwrap();
            [
                a[0].as_i64().unwrap() as i32,
                a[1].as_i64().unwrap() as i32,
                a[2].as_i64().unwrap() as i32,
            ]
        };
        let (e, r) = (arr("e"), arr("r"));
        let (s, niso) = make_removed_protons_string(
            c["nrp"].as_i64().unwrap() as i32,
            if bi & 2 != 0 { None } else { Some(&e) },
            if bi & 4 != 0 { None } else { Some(&r) },
            bi & 1 != 0,
        );
        assert_eq!(
            format!("{} {}", hex(s.as_bytes()), niso.unwrap_or(0)),
            out[10],
            "protons"
        );
    }
}

fn upto_nul(b: &[u8]) -> &[u8] {
    &b[..b.iter().position(|&c| c == 0).unwrap_or(b.len())]
}

#[test]
fn io_streams_match_the_c_code() {
    use rdoxide::external::inchi::ichi_io::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_io.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() >= 600);
    let ios_state = |s: &InchiIostream| {
        format!(
            " [{} {} {}]",
            s.s.n_allocated_length, s.s.n_used_length, s.s.n_ptr
        )
    };
    let ios_text = |s: &InchiIostream| match &s.s.p_str {
        Some(p) => hex(upto_nul(p)),
        None => hex(b""),
    };
    for c in cases {
        let expected = c["out"].as_str().unwrap();
        match c["kind"].as_str().unwrap() {
            "Q" => {
                let text = unhex(c["text"].as_str().unwrap());
                let mut ios = InchiIostream::default();
                inchi_ios_init(&mut ios, INCHI_IOS_TYPE_STRING, IosFile::None);
                inchi_ios_print(Some(&mut ios), &String::from_utf8(text.clone()).unwrap());
                ios.s.n_ptr = 0;
                let mut sb = InchiIosString::default();
                inchi_strbuf_init(&mut sb, 8, 4);
                let mut got = String::new();
                let line_of =
                    |r: &IosLine| format!("{} {} {}", r.length, r.too_long as i32, hex(&r.line));
                for op in c["ops"].as_array().unwrap() {
                    let o = op[0].as_str().unwrap();
                    let a = op[1].as_i64().unwrap() as i32;
                    let b = op[2].as_i64().unwrap();
                    match o {
                        "g" => got += &line_of(&inchi_ios_gets(a, &mut ios)),
                        "t" => got += &line_of(&inchi_ios_gets_tab(a, &mut ios)),
                        "1" => got += &line_of(&inchi_ios_gets_tab1(a, &mut ios)),
                        "s" => {
                            got += &inchi_sgets(a, &mut ios)
                                .map_or("N".to_string(), |l| format!("S {}", hex(&l)));
                        }
                        "f" => {
                            got += &inchi_fgets_lf(a, &mut ios)
                                .map_or("N".to_string(), |l| format!("S {}", hex(&l)));
                        }
                        "c" => {
                            got += &inchi_ios_str_getc(&mut ios)
                                .map_or(-1, i32::from)
                                .to_string();
                        }
                        _ => {
                            let r = inchi_strbuf_addline(&mut sb, &mut ios, a != 0, b != 0);
                            got += &format!(
                                "{r} {} {}",
                                hex(upto_nul(sb.p_str.as_ref().unwrap())),
                                sb.n_allocated_length
                            );
                        }
                    }
                    got += &format!(" #{};", ios.s.n_ptr);
                }
                assert_eq!(got, expected, "text {:?}", String::from_utf8_lossy(&text));
            }
            "P" => {
                let text = unhex(c["text"].as_str().unwrap());
                let mut o = InchiIostream::default();
                inchi_ios_init(&mut o, INCHI_IOS_TYPE_STRING, IosFile::None);
                let mut sb = InchiIosString::default();
                inchi_strbuf_init(&mut sb, 16, 8);
                for (k, chunk) in text.chunks(7).enumerate() {
                    let ch = String::from_utf8(chunk.to_vec()).unwrap();
                    inchi_ios_print(Some(&mut o), &ch);
                    inchi_strbuf_printf(&mut sb, &ch);
                    if k == 2 {
                        inchi_strbuf_printf_from(&mut sb, 3, "ab");
                    }
                }
                let mut got = format!(
                    "{} {}{}",
                    o.s.n_allocated_length,
                    ios_text(&o),
                    ios_state(&o)
                );
                got += &format!(
                    " | {} {} {} {}",
                    sb.n_allocated_length,
                    sb.n_used_length,
                    sb.n_ptr,
                    hex(upto_nul(sb.p_str.as_ref().unwrap()))
                );
                got += &format!(" {}", inchi_strbuf_update(&mut sb, 40));
                inchi_strbuf_reset(&mut sb);
                got += &format!(
                    " {} {} {}",
                    sb.n_allocated_length, sb.n_used_length, sb.n_ptr
                );
                inchi_ios_flush_not_displayed(Some(&mut o));
                got += &format!(" | {}{}", ios_text(&o), ios_state(&o));
                assert_eq!(got, expected, "print {:?}", String::from_utf8_lossy(&text));
            }
            _ => {
                let a = c["args"].as_array().unwrap();
                let (lab, val) = (a[0].as_str().unwrap(), a[1].as_str().unwrap());
                let mut o = InchiIostream::default();
                inchi_ios_init(&mut o, INCHI_IOS_TYPE_STRING, IosFile::None);
                output_record_info(
                    &mut o,
                    a[2].as_i64().unwrap() as i32,
                    a[3].as_i64().unwrap() != 0,
                    (!lab.is_empty()).then_some(lab),
                    (!val.is_empty()).then_some(val),
                    a[4].as_u64().unwrap(),
                    "\n",
                    "\t",
                );
                assert_eq!(
                    format!("{}{}", ios_text(&o), ios_state(&o)),
                    expected,
                    "record {a:?}"
                );
            }
        }
    }
}

#[test]
fn common_options_match_the_c_code() {
    use rdoxide::external::inchi::ichiparm::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_parm.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() >= 500);
    for c in cases {
        let ints = |v: &serde_json::Value| -> Vec<i64> {
            v.as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap())
                .collect()
        };
        let v = ints(&c["v"]);
        let mut o = CommonOptions {
            b_ver1_default_mode: c["vm"].as_u64().unwrap(),
            n_mode: c["nm"].as_u64().unwrap(),
            b_inchi_output_options: v[0] as i32,
            b_inchi_output_options2: v[1] as i32,
            b_std_format: v[2] as i32,
            b_hash_key: v[3] as i32,
            b_hash_xtra1: v[4] as i32,
            b_hash_xtra2: v[5] as i32,
            b_fix_sp3_bug: v[6] as i32,
            b_fix_fb2: v[7] as i32,
            b_add_phosphine_stereo: v[8] as i32,
            b_add_arsine_stereo: v[9] as i32,
            b_no_struct_labels: v[10] as i32,
            b_pointed_edge_stereo: v[11] as i32,
            b_do_not_add_h: v[12] as i32,
            b_forced_chiral_flag: v[13] as i32,
            b_reconnect_coord: v[14] as i32,
            b_keto_enol_taut: v[15] as i32,
            b15_taut_non_ring: v[16] as i32,
            b_pt_06_00_taut: v[17] as i32,
            b_pt_13_00_taut: v[18] as i32,
            b_pt_16_00_taut: v[19] as i32,
            b_pt_18_00_taut: v[20] as i32,
            b_pt_22_00_taut: v[21] as i32,
            b_pt_39_00_taut: v[22] as i32,
            b_loose_tsa_check: v[23] as i32,
            b_large_molecules: v[24] as i32,
            b_polymers: v[25] as i32,
            b_fold_polymer_sru: v[26] as i32,
            b_frame_shift_scheme: v[27] as i32,
            b_stereo_at_zz: v[28] as i32,
            b_np_zz: v[29] as i32,
            b_no_warnings: v[30] as i32,
            b_merge_hash: v[31] as i32,
            b_hide_inchi: v[32] as i32,
        };
        let mut ip = InputParms {
            n_input_type: v[33] as i32,
            ..Default::default()
        };
        let got = set_common_options_by_parg(
            c["arg"].as_str().unwrap(),
            c["dev"].as_i64().unwrap() != 0,
            &mut ip,
            &mut o,
        );
        let out = ints(&c["out"]);
        let mut mine = vec![
            i64::from(got),
            o.b_ver1_default_mode as i64,
            o.n_mode as i64,
        ];
        mine.extend(
            [
                o.b_inchi_output_options,
                o.b_inchi_output_options2,
                o.b_std_format,
                o.b_hash_key,
                o.b_hash_xtra1,
                o.b_hash_xtra2,
                o.b_fix_sp3_bug,
                o.b_fix_fb2,
                o.b_add_phosphine_stereo,
                o.b_add_arsine_stereo,
                o.b_no_struct_labels,
                o.b_pointed_edge_stereo,
                o.b_do_not_add_h,
                o.b_forced_chiral_flag,
                o.b_reconnect_coord,
                o.b_keto_enol_taut,
                o.b15_taut_non_ring,
                o.b_pt_06_00_taut,
                o.b_pt_13_00_taut,
                o.b_pt_16_00_taut,
                o.b_pt_18_00_taut,
                o.b_pt_22_00_taut,
                o.b_pt_39_00_taut,
                o.b_loose_tsa_check,
                o.b_large_molecules,
                o.b_polymers,
                o.b_fold_polymer_sru,
                o.b_frame_shift_scheme,
                o.b_stereo_at_zz,
                o.b_np_zz,
                o.b_no_warnings,
                o.b_merge_hash,
                o.b_hide_inchi,
            ]
            .map(i64::from),
        );
        mine.extend([
            i64::from(ip.n_input_type),
            ip.first_struct_number,
            ip.last_struct_number,
            i64::from(ip.b_save_all_good_structs_as_problem),
            i64::from(ip.b_filter_ss),
            i64::from(ip.b_fix_non_uniform_draw),
        ]);
        assert_eq!(mine, out, "arg {:?} dev {}", c["arg"], c["dev"]);
    }
}

#[test]
fn options_string_splitting_matches_the_c_code() {
    use rdoxide::external::inchi::ichiparm::parse_options_string;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_optstr.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() >= 600);
    for c in cases {
        let s = unhex(c["s"].as_str().unwrap());
        let argv = parse_options_string(&s, c["max"].as_u64().unwrap() as usize);
        assert_eq!(
            argv.len() as u64,
            c["n"].as_u64().unwrap(),
            "{:?}",
            String::from_utf8_lossy(&s)
        );
        let got: Vec<String> = argv.iter().map(|a| hex(a)).collect();
        let want: Vec<&str> = c["argv"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{:?}", String::from_utf8_lossy(&s));
    }
}

#[test]
fn read_command_line_parms_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_io::*;
    use rdoxide::external::inchi::ichiparm::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_readparms.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() >= 700);
    for c in cases {
        let args: Vec<String> = c["args"]
            .as_array()
            .unwrap()
            .iter()
            .map(|a| a.as_str().unwrap().to_string())
            .collect();
        let mut argv: Vec<&str> = vec![""];
        argv.extend(args.iter().map(String::as_str));
        let mut log = InchiIostream::default();
        inchi_ios_init(&mut log, INCHI_IOS_TYPE_STRING, IosFile::None);
        let mut ip = InputParms::default();
        let r = read_command_line_parms(&argv, &mut ip, c["rel"].as_i64().unwrap() != 0, &mut log);
        let ints: Vec<String> = [
            i64::from(r),
            i64::from(ip.n_input_type),
            ip.num_paths as i64,
            i64::from(ip.b_abc_numbers),
            i64::from(ip.b_ct_predecessors),
            ip.n_mode as i64,
            ip.b_taut_flags as i64,
            ip.b_taut_flags_done as i64,
            i64::from(ip.b_inchi_output_options),
            i64::from(ip.b_inchi_output_options2),
            i64::from(ip.b_display_each_component_inchi),
            i64::from(ip.b_display),
            i64::from(ip.b_display_if_restore_warnings),
            i64::from(ip.b_merge_all_input_structures),
            i64::from(ip.b_save_warning_structs_as_problem),
            i64::from(ip.b_save_all_good_structs_as_problem),
            i64::from(ip.b_get_sdfile_id),
            i64::from(ip.b_get_molfile_number),
            i64::from(ip.b_compare_components),
            i64::from(ip.b_display_composite_results),
            i64::from(ip.b_do_not_add_h),
            i64::from(ip.b_no_struct_labels),
            i64::from(ip.b_fix_non_uniform_draw),
            i64::from(ip.b_calc_inchi_hash),
            i64::from(ip.b_chiral_flag),
            i64::from(ip.b_allow_empty_structure),
            i64::from(ip.b_large_molecules),
            i64::from(ip.b_loose_tsa_check),
            i64::from(ip.b_polymers),
            i64::from(ip.b_fold_polymer_sru),
            i64::from(ip.b_frame_shift_scheme),
            i64::from(ip.b_np_zz),
            i64::from(ip.b_stereo_at_zz),
            i64::from(ip.b_merge_hash),
            i64::from(ip.b_no_warnings),
            i64::from(ip.b_hide_inchi),
            i64::from(ip.b_read_inchi_options),
            i64::from(ip.b_underivatize),
            i64::from(ip.b_ring2_chain),
            i64::from(ip.b_ignore_unchanged),
            i64::from(ip.b_filter_ss),
            ip.first_struct_number,
            ip.last_struct_number,
            ip.msec_max_time,
            ip.msec_left_time,
            ip.ul_displ_time,
            ip.l_sdf_id as i64,
            ip.l_molfile_number,
            i64::from(ip.p_sdf_label),
            i64::from(ip.p_sdf_value),
        ]
        .iter()
        .map(|v| v.to_string())
        .collect();
        let mut got = ints.join(" ");
        got += &format!(" {}", hex(ip.sz_sdf_data_header.as_bytes()));
        for p in &ip.path {
            got += &p
                .as_ref()
                .map_or(" N".to_string(), |p| format!(" {}", hex(p.as_bytes())));
        }
        let logtext = log
            .s
            .p_str
            .as_ref()
            .map_or(Vec::new(), |p| upto_nul(p).to_vec());
        got += &format!(" {}", hex(&logtext));
        assert_eq!(
            got,
            c["out"].as_str().unwrap(),
            "args {:?} rel {}",
            c["args"],
            c["rel"]
        );
    }
}

#[test]
fn std_format_enforcement_matches_the_c_code() {
    use rdoxide::external::inchi::ichiparm::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_stdenforce.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() >= 600);
    for c in cases {
        let i: Vec<u64> = c["in"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap())
            .collect();
        let want: Vec<u64> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap())
            .collect();
        let mut ip = InputParms {
            n_mode: i[0],
            b_taut_flags: i[1],
            b_inchi_output_options: i[2] as i32,
            ..Default::default()
        };
        enforce_std_format(&mut ip);
        let got = vec![ip.n_mode, ip.b_taut_flags, ip.b_inchi_output_options as u64];
        assert_eq!(got, want, "input {i:?}");
    }
}

#[test]
fn atom_input_conversion_matches_the_c_code() {
    use rdoxide::external::inchi::inchi_dll::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_setatoms.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() >= 1200);
    for c in cases {
        let inp = &c["in"];
        let atoms = inp["atoms"].as_array().unwrap();
        let mut ati = Vec::new();
        for a in atoms {
            let mut x = InchiAtom::default();
            for (d, s) in x
                .elname
                .iter_mut()
                .zip(a["el"].as_str().unwrap().bytes().take(5))
            {
                *d = s;
            }
            let xyz = a["xyz"].as_array().unwrap();
            x.x = xyz[0].as_f64().unwrap();
            x.y = xyz[1].as_f64().unwrap();
            x.z = xyz[2].as_f64().unwrap();
            x.isotopic_mass = a["iso"].as_u64().unwrap() as i16;
            x.radical = a["rad"].as_i64().unwrap() as i8;
            x.charge = a["chg"].as_i64().unwrap() as i8;
            for (k, v) in a["nh"].as_array().unwrap().iter().enumerate() {
                x.num_iso_h[k] = v.as_i64().unwrap() as i8;
            }
            let nb = a["nb"].as_array().unwrap();
            x.num_bonds = nb.len() as i16;
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_i64().unwrap() as i16;
                x.bond_type[j] = e[1].as_i64().unwrap() as i8;
                x.bond_stereo[j] = e[2].as_i64().unwrap() as i8;
            }
            ati.push(x);
        }
        let mut msg = String::new();
        let mut err = 0;
        let (at, n_dim, n_bonds, coords) =
            extract_atoms_and_bonds(&ati, inp["dna"].as_i64().unwrap() != 0, &mut msg, &mut err);
        let head = format!("{err} {n_dim} {n_bonds} {}", hex(msg.as_bytes()));
        assert_eq!(head, c["head"].as_str().unwrap(), "case {inp}");
        for (i, line) in c["atoms"].as_array().unwrap().iter().enumerate() {
            let a = &at[i];
            let mut s = format!("A {}", hex(upto_nul(&a.elname)));
            s += &format!(
                " {} {} {} {} {} {} {} {} {} {} {} {} 0 0",
                a.el_number,
                a.charge,
                a.radical,
                a.valence,
                a.chem_bonds_valence,
                a.num_h,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.iso_atw_diff,
                a.at_type,
                a.orig_at_number
            );
            s += &format!(" {:.6} {:.6} {:.6}", a.x, a.y, a.z);
            for k in 0..a.valence as usize {
                s += &format!(" {}:{}:{}", a.neighbor[k], a.bond_type[k], a.bond_stereo[k]);
            }
            s += &format!(" |{}", hex(upto_nul(&coords[i])));
            assert_eq!(s, line.as_str().unwrap(), "atom {i} of {inp}");
        }
    }
}

#[test]
fn zero_d_parities_match_the_c_code() {
    use rdoxide::external::inchi::inchi_dll::inp_atom_0d_to_inchi_atom;
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::readinch::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_0d.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() >= 1500);
    for c in cases {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            let ex = a["x"].as_array().unwrap();
            x.charge = ex[0].as_i64().unwrap() as i8;
            x.radical = ex[1].as_i64().unwrap() as i8;
            x.iso_atw_diff = ex[2].as_i64().unwrap() as i8;
            x.num_h = ex[3].as_i64().unwrap() as i8;
            for k in 0..3 {
                x.num_iso_h[k] = ex[4 + k].as_i64().unwrap() as i8;
            }
            x.set_elname(ex[7].as_str().unwrap());
            at.push(x);
        }
        let st: Vec<InchiStereo0D> = inp["desc"]
            .as_array()
            .unwrap()
            .iter()
            .map(|d| {
                let v: Vec<i64> = d
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_i64().unwrap())
                    .collect();
                InchiStereo0D {
                    neighbor: [v[0] as i16, v[1] as i16, v[2] as i16, v[3] as i16],
                    central_atom: v[4] as i16,
                    type_: v[5] as i8,
                    parity: v[6] as i8,
                }
            })
            .collect();
        let mut msg = String::new();
        let mut err = 0;
        extract_0d_parities(
            &mut at,
            &st,
            &mut msg,
            &mut err,
            inp["unk"].as_i64().unwrap() as i32,
        );
        assert_eq!(
            format!("{err} {}", hex(msg.as_bytes())),
            c["head"].as_str().unwrap(),
            "case {inp}"
        );
        for (i, line) in c["atoms"].as_array().unwrap().iter().enumerate() {
            let a = &at[i];
            let mut s = format!(
                "A {} {} {} | {} {} {} | {} {} {} | {} {} {} | {} | {} {} {} {} |",
                a.sb_parity[0],
                a.sb_parity[1],
                a.sb_parity[2],
                a.sb_ord[0],
                a.sb_ord[1],
                a.sb_ord[2],
                a.sn_ord[0],
                a.sn_ord[1],
                a.sn_ord[2],
                a.sn_orig_at_num[0],
                a.sn_orig_at_num[1],
                a.sn_orig_at_num[2],
                a.p_parity,
                a.p_orig_at_num[0],
                a.p_orig_at_num[1],
                a.p_orig_at_num[2],
                a.p_orig_at_num[3]
            );
            for k in 0..a.valence as usize {
                s += &format!(" {}", a.bond_stereo[k]);
            }
            assert_eq!(s, line.as_str().unwrap(), "atom {i} of {inp}");
        }
        // and back to the API structures
        let (ret, atoms, stereo) = inp_atom_0d_to_inchi_atom(&at);
        let back = c["back"].as_array().unwrap();
        assert_eq!(
            format!("R {ret} {} {}", atoms.len(), stereo.len()),
            back[0].as_str().unwrap(),
            "back of {inp}"
        );
        for (i, a) in atoms.iter().enumerate() {
            let mut s = format!(
                "B {} {} {} {} {} {} {} {} 0 {}",
                a.num_bonds,
                a.charge,
                a.radical,
                a.isotopic_mass,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.num_iso_h[3],
                hex(upto_nul(&a.elname))
            );
            for k in 0..a.num_bonds as usize {
                s += &format!(" {}:{}", a.neighbor[k], a.bond_type[k]);
            }
            assert_eq!(s, back[1 + i].as_str().unwrap(), "back atom {i} of {inp}");
        }
        for (i, t) in stereo.iter().enumerate() {
            let s = format!(
                "T {} {} {} {} {} {} {}",
                t.type_,
                t.central_atom,
                t.parity,
                t.neighbor[0],
                t.neighbor[1],
                t.neighbor[2],
                t.neighbor[3]
            );
            assert_eq!(
                s,
                back[1 + atoms.len() + i].as_str().unwrap(),
                "back stereo {i} of {inp}"
            );
        }
    }
}

#[test]
fn drawing_corrections_match_the_c_code() {
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::strutil::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_strutil.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() >= 1500);
    for c in cases {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let nb = a["nb"].as_array().unwrap();
            let mut x = InpAtom {
                valence: nb.len() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                charge: a["chg"].as_i64().unwrap() as i8,
                radical: a["rad"].as_i64().unwrap() as i8,
                iso_atw_diff: a["iso"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                component: a["comp"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
                x.bond_stereo[j] = e[2].as_i64().unwrap() as i8;
            }
            at.push(x);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut got: Vec<String> = Vec::new();
        for q in inp["q"].as_array().unwrap() {
            let atno = q[0].as_u64().unwrap() as usize;
            let mut ord: Vec<i32> = q[1]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap() as i32)
                .collect();
            let n = ord.len();
            let r = n_find_one_om(&at, atno, &mut ord, n);
            let mut s = format!("Q {r}");
            for o in &ord {
                s += &format!(" {o}");
            }
            got.push(s);
        }
        for i in 0..at.len() {
            got.push(match the_only_doublet_neigh(&at, i) {
                Some((n, a, b)) => format!("D {n} {a} {b}"),
                None => "D -1 -9 -9".to_string(),
            });
        }
        for i in 0..at.len().saturating_sub(1) {
            got.push(format!(
                "C {}",
                cmp_iso_atw_diff_component_no(&at[i], &at[i + 1])
            ));
        }
        let (mut n1, mut n2) = (0, 0);
        fix_non_uniform_drawn_oxoanions(&mut at, &mut n1);
        fix_non_uniform_drawn_amidiniums(&mut at, &mut n2);
        got.push(format!("N {n1} {n2}"));
        for a in &at {
            let mut s = format!("A {} {} {}", a.charge, a.chem_bonds_valence, a.valence);
            for k in 0..a.valence as usize {
                s += &format!(" {}:{}", a.bond_type[k], a.bond_stereo[k]);
            }
            got.push(s);
        }
        assert_eq!(got, want, "case {inp}");
    }
}

#[test]
fn ion_pairs_and_odd_things_match_the_c_code() {
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::strutil::{fix_odd_things, remove_ion_pairs};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_ionpairs.json")).unwrap();
    let cases = data.as_array().unwrap();
    assert!(cases.len() >= 12000);
    let mut changed = 0;
    for c in cases {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let nb = a["nb"].as_array().unwrap();
            let mut x = InpAtom {
                valence: nb.len() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                charge: a["chg"].as_i64().unwrap() as i8,
                radical: a["rad"].as_i64().unwrap() as i8,
                iso_atw_diff: a["iso"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                ..Default::default()
            };
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
                x.bond_stereo[j] = e[2].as_i64().unwrap() as i8;
            }
            at.push(x);
        }
        let r = if inp["mode"] == 1 {
            fix_odd_things(&mut at, inp["fb"] == 1, inp["fn"] == 1)
        } else {
            remove_ion_pairs(&mut at)
        };
        let mut got = vec![format!("R {r}")];
        for a in &at {
            let mut s = format!(
                "A {} {} {} {}",
                a.charge, a.chem_bonds_valence, a.radical, a.valence
            );
            for k in 0..a.valence as usize {
                s += &format!(" {}:{}:{}", a.neighbor[k], a.bond_type[k], a.bond_stereo[k]);
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
        changed += usize::from(r != 0);
    }
    assert!(changed > 1000);
}

#[test]
fn bond_removal_matches_the_c_code() {
    use rdoxide::external::inchi::inpdef::{InpAtom, MAX_NUM_STEREO_BONDS};
    use rdoxide::external::inchi::strutil::{disconnect_inp_at_bond, remove_inp_at_bond};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bondremove.json")).unwrap();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for (i, a) in inp["atoms"].as_array().unwrap().iter().enumerate() {
            let nb = a["nb"].as_array().unwrap();
            let mut x = InpAtom {
                valence: nb.len() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                p_parity: a["pp"].as_i64().unwrap() as i8,
                component: a["cmp"].as_u64().unwrap() as u16,
                orig_at_number: i as u16 + 1,
                ..Default::default()
            };
            for (m, v) in a["po"].as_array().unwrap().iter().enumerate() {
                x.p_orig_at_num[m] = v.as_u64().unwrap() as u16;
            }
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
                x.bond_stereo[j] = e[2].as_i64().unwrap() as i8;
            }
            for (m, e) in a["sb"].as_array().unwrap().iter().enumerate() {
                x.sb_ord[m] = e[0].as_i64().unwrap() as i8;
                x.sn_ord[m] = e[1].as_i64().unwrap() as i8;
                x.sb_parity[m] = e[2].as_i64().unwrap() as i8;
                x.sn_orig_at_num[m] = e[3].as_u64().unwrap() as u16;
            }
            at.push(x);
        }
        let mut comp = [1u16; 16];
        let mut got: Vec<String> = Vec::new();
        for r in inp["rem"].as_array().unwrap() {
            let (mode, iat, k) = (
                r[0].as_i64().unwrap(),
                r[1].as_u64().unwrap() as usize,
                r[2].as_u64().unwrap() as usize,
            );
            let ret = if mode == 1 {
                i32::from(disconnect_inp_at_bond(&mut at, Some(&mut comp), iat, k))
            } else {
                remove_inp_at_bond(&mut at, iat, k)
            };
            got.push(format!("R {ret} {} {} {}", comp[0], comp[1], comp[2]));
        }
        for a in &at {
            let mut s = format!("A {} {} {}", a.valence, a.chem_bonds_valence, a.p_parity);
            for m in 0..4 {
                s += &format!(" {}", a.p_orig_at_num[m]);
            }
            for k in 0..a.valence as usize {
                s += &format!(" {}:{}:{}", a.neighbor[k], a.bond_type[k], a.bond_stereo[k]);
            }
            for m in 0..MAX_NUM_STEREO_BONDS {
                s += &format!(
                    " S{}:{}:{}:{}",
                    a.sb_ord[m], a.sn_ord[m], a.sb_parity[m], a.sn_orig_at_num[m]
                );
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn salt_disconnection_matches_the_c_code() {
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::strutil::{disconnect_salts, may_disconnect_metals};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_salts.json")).unwrap();
    let mut changed = 0;
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for (i, a) in inp["atoms"].as_array().unwrap().iter().enumerate() {
            let nb = a["nb"].as_array().unwrap();
            let mut x = InpAtom {
                valence: nb.len() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                charge: a["chg"].as_i64().unwrap() as i8,
                radical: a["rad"].as_i64().unwrap() as i8,
                iso_atw_diff: a["iso"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: i as u16 + 1,
                x: a["xyz"][0].as_f64().unwrap(),
                y: a["xyz"][1].as_f64().unwrap(),
                z: a["xyz"][2].as_f64().unwrap(),
                ..Default::default()
            };
            for m in 0..3 {
                x.num_iso_h[m] = a["ih"][m].as_i64().unwrap() as i8;
            }
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
                x.bond_stereo[j] = e[2].as_i64().unwrap() as i8;
            }
            at.push(x);
        }
        let mode = inp["mode"].as_i64().unwrap();
        let mut got = Vec::new();
        if mode <= 1 {
            let mut nbonds = inp["nb"].as_i64().unwrap() as i32;
            let r = disconnect_salts(&mut at, &mut nbonds, mode == 1);
            got.push(format!("R {r} {nbonds}"));
            changed += usize::from(r != 0);
        } else {
            let (r, bdc) = may_disconnect_metals(&at);
            got.push(format!("R {r} {bdc} 0"));
        }
        for a in &at {
            let mut s = format!(
                "A {} {} {} {} {} {} {} {} {}",
                a.valence,
                a.chem_bonds_valence,
                a.charge,
                a.num_h,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.iso_atw_diff,
                a.radical
            );
            for k in 0..a.valence as usize {
                s += &format!(" {}:{}:{}", a.neighbor[k], a.bond_type[k], a.bond_stereo[k]);
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(changed > 500);
}

#[test]
fn metal_disconnection_matches_the_c_code() {
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::strutil::{
        OrigAtomData, disconnect_metals, may_disconnect_metals,
    };
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_metals.json")).unwrap();
    let mut coord_mismatch = 0;
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for (i, a) in inp["atoms"].as_array().unwrap().iter().enumerate() {
            let nb = a["nb"].as_array().unwrap();
            let mut x = InpAtom {
                valence: nb.len() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                charge: a["chg"].as_i64().unwrap() as i8,
                radical: a["rad"].as_i64().unwrap() as i8,
                iso_atw_diff: a["iso"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                component: a["comp"].as_u64().unwrap() as u16,
                orig_at_number: i as u16 + 1,
                x: a["xyz"][0].as_f64().unwrap(),
                y: a["xyz"][1].as_f64().unwrap(),
                z: a["xyz"][2].as_f64().unwrap(),
                ..Default::default()
            };
            for m in 0..3 {
                x.num_iso_h[m] = a["ih"][m].as_i64().unwrap() as i8;
            }
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
                x.bond_stereo[j] = e[2].as_i64().unwrap() as i8;
            }
            at.push(x);
        }
        let check = inp["check"] == 1;
        let (may, bdc) = may_disconnect_metals(&at);
        let mut od = OrigAtomData {
            num_inp_atoms: at.len(),
            at,
            b_disconnect_coord: bdc,
            n_old_comp_number: Some(vec![1u16; 16]),
            ..Default::default()
        };
        let mut flags = 0u64;
        let r = disconnect_metals(&mut od, check, Some(&mut flags));
        let old = od.n_old_comp_number.as_ref().unwrap();
        let mut got = vec![format!(
            "R {may} {bdc} {r} {} {flags} {} {} {}",
            od.num_inp_atoms, old[0], old[1], old[2]
        )];
        for a in &od.at[..od.num_inp_atoms] {
            let mut s = format!(
                "A {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                a.el_number,
                a.valence,
                a.chem_bonds_valence,
                a.charge,
                a.num_h,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.iso_atw_diff,
                a.radical,
                a.component,
                a.orig_at_number,
                a.x.to_bits(),
                a.y.to_bits(),
                a.z.to_bits()
            );
            for k in 0..a.valence as usize {
                s += &format!(" {}:{}:{}", a.neighbor[k], a.bond_type[k], a.bond_stereo[k]);
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        // the coordinates come from sin/cos/atan2 of the platform libm: compare
        // them to 1e-12 (a few boundary cases may differ), everything else exactly
        assert_eq!(got.len(), want.len(), "{inp}");
        for (g, w) in got.iter().zip(&want) {
            let (gt, wt): (Vec<&str>, Vec<&str>) = (g.split(' ').collect(), w.split(' ').collect());
            assert_eq!(gt.len(), wt.len(), "{inp}");
            let coord = g.starts_with("A ");
            for (q, (x, y)) in gt.iter().zip(&wt).enumerate() {
                if coord && (13..16).contains(&q) {
                    let fx = f64::from_bits(x.parse().unwrap());
                    let fy = f64::from_bits(y.parse().unwrap());
                    // (libm differences can flip the sector choice of an H: counted)
                    coord_mismatch += usize::from((fx - fy).abs() > 1e-12 * fy.abs().max(1.0));
                } else {
                    assert_eq!(x, y, "{inp}");
                }
            }
        }
    }
    assert!(
        coord_mismatch <= 20,
        "{coord_mismatch} coordinate mismatches"
    );
}

#[test]
fn terminal_hydrogen_removal_matches_the_c_code() {
    use rdoxide::external::inchi::inpdef::{InpAtom, MAX_NUM_STEREO_BONDS};
    use rdoxide::external::inchi::strutil::{add_dt_to_num_h, remove_terminal_hdt};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_termh.json")).unwrap();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let nb = a["nb"].as_array().unwrap();
            let ch = |k: &str| match a[k].as_str().unwrap() {
                "_" => 0u8,
                s => s.as_bytes()[0],
            };
            let mut x = InpAtom {
                valence: nb.len() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                charge: a["chg"].as_i64().unwrap() as i8,
                iso_atw_diff: a["iso"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                ..Default::default()
            };
            x.elname[0] = ch("e0");
            x.elname[1] = ch("e1");
            for m in 0..3 {
                x.num_iso_h[m] = a["ih"][m].as_i64().unwrap() as i8;
            }
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
                x.bond_stereo[j] = e[2].as_i64().unwrap() as i8;
            }
            for (m, e) in a["sb"].as_array().unwrap().iter().take(3).enumerate() {
                x.sb_ord[m] = e[0].as_i64().unwrap() as i8;
                x.sn_ord[m] = e[1].as_i64().unwrap() as i8;
                x.sb_parity[m] = e[2].as_i64().unwrap() as i8;
            }
            at.push(x);
        }
        let n = at.len();
        let r = remove_terminal_hdt(&mut at, inp["fix"] == 1);
        if inp["add"] == 1 {
            let m = if r > 0 { r as usize } else { n };
            add_dt_to_num_h(&mut at[..m]);
        }
        let mut got = vec![format!("R {r}")];
        for a in &at {
            let mut s = format!(
                "A {} {} {} {} {} {} {} {} {} {} {}",
                a.elname[0],
                a.elname[1],
                a.valence,
                a.chem_bonds_valence,
                a.charge,
                a.iso_atw_diff,
                a.num_h,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.component
            );
            for k in 0..a.valence as usize {
                s += &format!(" {}:{}:{}", a.neighbor[k], a.bond_type[k], a.bond_stereo[k]);
            }
            for m in 0..MAX_NUM_STEREO_BONDS {
                s += &format!(" S{}:{}:{}", a.sb_ord[m], a.sn_ord[m], a.sb_parity[m]);
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn component_marking_matches_the_c_code() {
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::strutil::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_components.json")).unwrap();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let nb = a["nb"].as_array().unwrap();
            let mut x = InpAtom {
                valence: nb.len() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                charge: a["chg"].as_i64().unwrap() as i8,
                radical: a["rad"].as_i64().unwrap() as i8,
                iso_atw_diff: a["iso"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                component: a["comp"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for m in 0..3 {
                x.num_iso_h[m] = a["ih"][m].as_i64().unwrap() as i8;
            }
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let mut got: Vec<String> = Vec::new();
        got.push(format!("B {}", num_hetero_atom_has_isotopic_h(&at)));
        let mut h = String::from("H");
        for i in 0..n {
            h += &format!(
                " {}:{}",
                get_iat_number(at[i].el_number),
                hetero_atom_may_have_xchg_iso_h(&at, i)
            );
        }
        got.push(h);
        let old: Vec<u16> = inp["old"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as u16)
            .collect();
        let mut od = OrigAtomData {
            at,
            num_inp_atoms: n,
            num_components: inp["prevnc"].as_i64().unwrap() as i32,
            n_old_comp_number: (inp["hasold"] == 1).then_some(old),
            ..Default::default()
        };
        let r = mark_disconnected_components(&mut od, inp["bproc"] == 1);
        got.push(format!("R {r} {}", od.num_components));
        let mut s = String::from("C");
        for a in &od.at {
            s += &format!(" {}", a.component);
        }
        got.push(s);
        let mut s = String::from("L");
        let mut o = String::from("O");
        if r > 0 {
            for v in &od.n_cur_at_len.as_ref().unwrap()[..r as usize] {
                s += &format!(" {v}");
            }
            for v in od.n_old_comp_number.as_ref().unwrap() {
                o += &format!(" {v}");
            }
        }
        got.push(s);
        got.push(o);
        for k in 1..=r {
            let mut ca = Vec::new();
            let m = extract_connected_component(&od.at, k as u16, &mut ca);
            let mut s = format!("X {m}");
            for a in &ca {
                s += &format!(" |{}", a.orig_compt_at_numb);
                for j in 0..a.valence as usize {
                    s += &format!(" {}", a.neighbor[j]);
                }
            }
            got.push(s);
        }
        set_connected_component_number(&mut od.at, 7);
        got.push(format!(
            "S {} {}",
            od.at[0].component,
            od.at[n - 1].component
        ));
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

fn dump_orig(tag: &str, od: &rdoxide::external::inchi::strutil::OrigAtomData) -> Vec<String> {
    use rdoxide::external::inchi::inpdef::MAX_NUM_STEREO_BONDS;
    let mut out = vec![format!(
        "{tag} {} {} {} {}",
        od.num_inp_atoms, od.b_disconnect_salts, od.b_disconnect_coord, od.num_components
    )];
    if od.at.is_empty() {
        out.push("-".into());
        return out;
    }
    let nc = od.num_components.max(0) as usize;
    let list = |name: &str, v: &Option<Vec<u16>>| {
        let mut s = name.to_string();
        if let Some(v) = v {
            for x in &v[..nc] {
                s += &format!(" {x}");
            }
        }
        s
    };
    out.push(list("L", &od.n_cur_at_len));
    out.push(list("O", &od.n_old_comp_number));
    for a in &od.at[..od.num_inp_atoms] {
        let mut s = format!(
            "A {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
            a.el_number,
            a.valence,
            a.chem_bonds_valence,
            a.charge,
            a.num_h,
            a.num_iso_h[0],
            a.num_iso_h[1],
            a.num_iso_h[2],
            a.iso_atw_diff,
            a.radical,
            a.component,
            a.orig_at_number,
            a.x.to_bits(),
            a.y.to_bits(),
            a.z.to_bits()
        );
        for k in 0..a.valence as usize {
            s += &format!(" {}:{}:{}", a.neighbor[k], a.bond_type[k], a.bond_stereo[k]);
        }
        for m in 0..MAX_NUM_STEREO_BONDS {
            s += &format!(
                " S{}:{}:{}:{}",
                a.sb_ord[m], a.sn_ord[m], a.sb_parity[m], a.sn_orig_at_num[m]
            );
        }
        out.push(s);
    }
    out
}

#[test]
fn preprocess_one_structure_matches_the_c_code() {
    use rdoxide::external::inchi::ichiparm::InputParms;
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::runichi3::{StructData, preprocess_one_structure};
    use rdoxide::external::inchi::strutil::OrigAtomData;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_preprocess.json")).unwrap();
    let mut coord_mismatch = 0;
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for (i, a) in inp["atoms"].as_array().unwrap().iter().enumerate() {
            let nb = a["nb"].as_array().unwrap();
            let mut x = InpAtom {
                valence: nb.len() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                charge: a["chg"].as_i64().unwrap() as i8,
                radical: a["rad"].as_i64().unwrap() as i8,
                iso_atw_diff: a["iso"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: i as u16 + 1,
                x: a["xyz"][0].as_f64().unwrap(),
                y: a["xyz"][1].as_f64().unwrap(),
                z: a["xyz"][2].as_f64().unwrap(),
                ..Default::default()
            };
            for (k, b) in a["sym"].as_str().unwrap().bytes().enumerate() {
                x.elname[k] = b;
            }
            for m in 0..3 {
                x.num_iso_h[m] = a["ih"][m].as_i64().unwrap() as i8;
            }
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
                x.bond_stereo[j] = e[2].as_i64().unwrap() as i8;
            }
            at.push(x);
        }
        let ip = InputParms {
            b_taut_flags: inp["flags"].as_u64().unwrap(),
            b_fix_non_uniform_draw: inp["fnu"].as_i64().unwrap() as i32,
            b_no_warnings: inp["nowarn"].as_i64().unwrap() as i32,
            ..Default::default()
        };
        let mut orig = OrigAtomData {
            num_inp_atoms: at.len(),
            num_inp_bonds: inp["nbonds"].as_i64().unwrap() as i32,
            at,
            ..Default::default()
        };
        let mut prep = [OrigAtomData::default(), OrigAtomData::default()];
        let mut sd = StructData::default();
        let r = preprocess_one_structure(&mut sd, &ip, &mut orig, &mut prep);
        let mut got = vec![format!(
            "R {r} {} {} [{}] {} {} {} {}",
            sd.n_error_type,
            sd.n_struct_read_error,
            sd.str_err_struct,
            sd.taut_flags[0],
            sd.taut_flags[1],
            sd.taut_flags_done[0],
            sd.taut_flags_done[1]
        )];
        got.extend(dump_orig("O", &orig));
        got.extend(dump_orig("P", &prep[0]));
        got.extend(dump_orig("Q", &prep[1]));
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        // coordinates of added hydrogens depend on the platform libm: counted
        assert_eq!(got.len(), want.len(), "{inp}");
        for (g, w) in got.iter().zip(&want) {
            let (gt, wt): (Vec<&str>, Vec<&str>) = (g.split(' ').collect(), w.split(' ').collect());
            assert_eq!(gt.len(), wt.len(), "{inp}: {g} vs {w}");
            let coord = g.starts_with("A ");
            for (q, (x, y)) in gt.iter().zip(&wt).enumerate() {
                if coord && (13..16).contains(&q) {
                    let fx = f64::from_bits(x.parse().unwrap());
                    let fy = f64::from_bits(y.parse().unwrap());
                    coord_mismatch += usize::from((fx - fy).abs() > 1e-12 * fy.abs().max(1.0));
                } else {
                    assert_eq!(x, y, "{inp}: {g} vs {w}");
                }
            }
        }
    }
    assert!(
        coord_mismatch <= 10,
        "{coord_mismatch} coordinate mismatches"
    );
}

#[test]
fn ring_systems_match_the_c_code() {
    use rdoxide::external::inchi::ichinorm::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_ringsys.json")).unwrap();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let nb = a["nb"].as_array().unwrap();
            let mut x = InpAtom {
                valence: nb.len() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                at_type: 5,
                c_flags: 3,
                component: 4,
                orig_compt_at_numb: 6,
                ..Default::default()
            };
            for m in 0..3 {
                x.num_iso_h[m] = a["ih"][m].as_i64().unwrap() as i8;
            }
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e.as_u64().unwrap() as u16;
            }
            at.push(x);
        }
        let r = mark_ring_systems_inp(&mut at, inp["start"].as_u64().unwrap() as usize);
        let mut got = vec![format!("R {r}")];
        for a in &at {
            got.push(format!(
                "A {} {} {} {}",
                a.n_block_system, a.b_cut_vertex, a.n_ring_system, a.n_num_at_in_ring_system
            ));
        }
        unmark_ring_systems_inp(&mut at);
        unmark_other_indicators(&mut at);
        unmark_one_component(&mut at);
        subtract_dt_from_num_h(&mut at);
        for a in &at {
            got.push(format!(
                "B {} {} {} {} {} {} {}",
                a.n_block_system,
                a.b_cut_vertex,
                a.n_ring_system,
                a.at_type,
                a.c_flags,
                a.component,
                a.num_h
            ));
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn tautomer_endpoints_match_the_c_code() {
    type InfoFn = fn(
        &[rdoxide::external::inchi::inpdef::InpAtom],
        usize,
        &mut rdoxide::external::inchi::ichitaut::EndpointInfo,
    ) -> i32;
    use rdoxide::external::inchi::ichitaut::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::util::get_endpoint_valence_ket;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_endpoint.json")).unwrap();
    let pe = |e: &EndpointInfo, r: i32| {
        if r > 0 {
            format!(
                " {r}:{},{},{},{},{},{},0",
                e.c_moveable_charge,
                e.c_neutral_bonds_valence,
                e.c_mobile,
                e.c_donor,
                e.c_acceptor,
                e.c_keto_enol_code
            )
        } else {
            " 0".to_string()
        }
    };
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let nb = a["nb"].as_array().unwrap();
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                charge: a["chg"].as_i64().unwrap() as i8,
                radical: a["rad"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                c_point: a["cp"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                at_type: a["att"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (k, b) in a["sym"].as_str().unwrap().bytes().enumerate() {
                x.elname[k] = b;
            }
            for m in 0..3 {
                x.num_iso_h[m] = a["ih"][m].as_i64().unwrap() as i8;
            }
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let mut got: Vec<String> = Vec::new();
        for i in 0..at.len() {
            let mut sub = 99i8;
            let ct = get_charge_type(&at, i, &mut sub);
            let el = at[i].el_number;
            let mut s = format!(
                "G {ct} {sub} {} {} {} {}",
                i32::from(is_centerpoint_elem(el)),
                i32::from(is_centerpoint_elem_ket(el)),
                i32::from(is_centerpoint_elem_strict(el)),
                get_endpoint_valence_ket(el)
            );
            let fs: [InfoFn; 8] = [
                n_get_endpoint_info,
                n_get_endpoint_info_pt_22_00,
                n_get_endpoint_info_pt_16_00,
                n_get_endpoint_info_pt_06_00,
                n_get_endpoint_info_pt_39_00,
                n_get_endpoint_info_pt_13_00,
                n_get_endpoint_info_pt_18_00,
                n_get_endpoint_info_ket,
            ];
            for f in fs {
                let mut e = EndpointInfo::default();
                let r = f(&at, i, &mut e);
                s += &pe(&e, r);
            }
            got.push(s);
            let mut ep = TEndpoint::default();
            add_end_point(&mut ep, &at, i);
            let mut s = format!(
                "E {} {} {}",
                ep.n_atom_number, ep.n_equ_number, ep.n_group_number
            );
            for v in ep.num {
                s += &format!(" {v}");
            }
            for v in ep.num_da {
                s += &format!(" {v}");
            }
            let mut num: [u16; 5] = [9, 8, 7, 6, 5];
            let mut da: [u16; 6] = [9, 8, 7, 6, 5, 4];
            for mode in 0..3 {
                let m = add_atom2num(&mut num, &at[i], mode);
                add_atom2da(&mut da, &at[i], mode);
                s += &format!(" |{m}");
                for v in num {
                    s += &format!(" {v}");
                }
                for v in da {
                    s += &format!(" {v}");
                }
                if mode == 1 {
                    for (k, v) in num.iter_mut().enumerate() {
                        *v = 9 - k as u16;
                    }
                    for (k, v) in da.iter_mut().enumerate() {
                        *v = 9 - k as u16;
                    }
                }
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn charge_groups_and_salt_types_match_the_c_code() {
    use rdoxide::external::inchi::ichitaut::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_cgroups.json")).unwrap();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let nb = a["nb"].as_array().unwrap();
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                charge: a["chg"].as_i64().unwrap() as i8,
                radical: a["rad"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                c_point: a["cp"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (k, b) in a["sym"].as_str().unwrap().bytes().enumerate() {
                x.elname[k] = b;
            }
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let hastg = inp["hastg"].as_i64().unwrap();
        let tgs: Vec<TGroup> = inp["tgs"]
            .as_array()
            .unwrap()
            .iter()
            .map(|g| {
                let mut t = TGroup {
                    n_group_number: g[0].as_u64().unwrap() as u16,
                    ..Default::default()
                };
                t.num[0] = g[1].as_u64().unwrap() as u16;
                t.num[1] = g[2].as_u64().unwrap() as u16;
                t
            })
            .collect();
        let tgi = TGroupInfo {
            num_t_groups: tgs.len(),
            t_group: (hastg != 0).then_some(tgs),
            ..Default::default()
        };
        let tg_opt = (hastg != 2).then_some(&tgi);
        let mut got: Vec<String> = Vec::new();
        for i in 0..at.len() {
            let (mut s1, mut s2, mut s3, mut s4) = (0, 0, 0, 0);
            let r1 = get_salt_charge_type(&at, i, tg_opt, &mut s1);
            let r2 = get_other_salt_charge_type(&at, i, tg_opt, &mut s2, false);
            let r3 = get_other_salt_charge_type(&at, i, tg_opt, &mut s3, true);
            let r4 = get_other_salt_type(&at, i, &mut s4);
            got.push(format!(
                "S {r1} {s1} {r2} {s2} {r3} {s3} {r4} {s4} {}",
                i32::from(do_not_merge_non_taut_atom(&at, i))
            ));
        }
        let maxc = inp["maxc"].as_u64().unwrap() as usize;
        let mut cg = vec![CGroup::default(); maxc];
        let mut numc = 0usize;
        for r in inp["regs"].as_array().unwrap() {
            let ret = register_c_points(
                &mut cg,
                &mut numc,
                r[0].as_u64().unwrap() as usize,
                r[1].as_u64().unwrap() as usize,
                r[2].as_u64().unwrap() as u8,
                &mut at,
            );
            let mut s = format!("C {ret} {numc} |");
            for g in &cg[..numc] {
                s += &format!(
                    " {},{},{},{},{}",
                    g.num[0], g.num[1], g.num_c_points, g.n_group_number, g.c_group_type
                );
            }
            s += " |";
            for a in &at {
                s += &format!(" {}", a.c_point);
            }
            got.push(s);
        }
        let mut cc: Vec<CCandidate> = inp["cands"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| CCandidate {
                atnumber: x[0].as_u64().unwrap() as u16,
                type_: x[1].as_i64().unwrap() as i8,
                subtype: x[2].as_i64().unwrap() as i8,
            })
            .collect();
        cc.sort_by(|a, b| cmp_c_candidates(a, b).cmp(&0));
        let mut s = String::from("Q");
        for x in &cc {
            s += &format!(" {},{},{}", x.atnumber, x.type_, x.subtype);
        }
        got.push(s);
        let mut sc: Vec<SCandidate> = inp["scands"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| SCandidate {
                atnumber: x[0].as_u64().unwrap() as u16,
                type_: x[1].as_i64().unwrap() as i8,
                subtype: x[2].as_i64().unwrap() as i8,
                endpoint: x[3].as_u64().unwrap() as u16,
            })
            .collect();
        sc.sort_by(|a, b| comp_candidates(a, b).cmp(&0));
        let mut s = String::from("T");
        for x in &sc {
            s += &format!(" {},{},{},{}", x.atnumber, x.type_, x.subtype, x.endpoint);
        }
        got.push(s);
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

fn dump_bns(p: &rdoxide::external::inchi::ichi_bns::BnStruct, out: &mut Vec<String>) {
    use rdoxide::external::inchi::ichi_bns::*;
    out.push(format!(
        "H {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
        p.num_atoms,
        p.num_added_atoms,
        p.n_max_add_atoms,
        p.num_c_groups,
        p.num_t_groups,
        p.num_vertices,
        p.num_bonds,
        p.num_edges,
        p.num_added_edges,
        p.n_max_add_edges,
        p.max_vertices,
        p.max_edges,
        p.max_iedges,
        p.tot_st_cap,
        p.tot_st_flow,
        p.b_not_a_simple_path
    ));
    for i in 0..p.num_vertices as usize {
        let v = &p.vert[i];
        out.push(format!(
            "V {} {} {} {} {} {} {} {} {} 0",
            v.type_,
            v.num_adj_edges,
            v.max_adj_edges,
            v.iedge,
            v.st_edge.cap,
            v.st_edge.cap0,
            v.st_edge.flow,
            v.st_edge.flow0,
            v.st_edge.pass
        ));
    }
    for i in 0..p.num_edges as usize {
        let e = &p.edge[i];
        out.push(format!(
            "E {} {} {} {} {} {} {} {} {} {}",
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
    let mut s = String::from("I");
    for i in 0..p.num_atoms as usize {
        for k in 0..p.vert[i].num_adj_edges as usize {
            s += &format!(" {}", p.iedge[p.vert[i].iedge + k]);
        }
    }
    out.push(s);
    if let Some(k) = p.alt_path {
        let ap = &p.altp[k];
        let len = ap[I_ALTP_PATH_LEN].number();
        let mut s = format!(
            "P {} {} {} {} {}",
            ap[I_ALTP_START_ATOM].number(),
            ap[I_ALTP_END_ATOM].number(),
            len,
            ap[I_ALTP_FLOW].flow[0],
            ap[I_ALTP_FLOW].flow[1]
        );
        for j in 0..len as usize {
            let e = &ap[I_ALTP_NEIGHBOR + j];
            s += &format!(" {},{}", e.ineigh(0), e.ineigh(1));
        }
        out.push(s);
    }
}

#[test]
fn balanced_network_structure_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bns.json")).unwrap();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let nb = a["nb"].as_array().unwrap();
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (j, e) in nb.iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, nchg) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        let mut got: Vec<String> = vec![format!("N {nchg}")];
        for a in &at {
            let mut s = format!("B {}", a.valence);
            for j in 0..a.valence as usize {
                s += &format!(" {}", a.bond_type[j]);
            }
            got.push(s);
        }
        dump_bns(&p, &mut got);
        p.type_tacn = inp["tacn"].as_u64().unwrap() as u16;
        p.type_t = inp["t"].as_u64().unwrap() as u16;
        p.type_cn = inp["cn"].as_u64().unwrap() as u16;
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        p.alt_path = Some(0);
        let mut switch_edge: Vec<Edge> = vec![[NO_VERTEX, 0]; 4096];
        let mut y = 3;
        while y < 2 * p.num_vertices + 2 && p.num_edges > 0 {
            if (y * 7 + n as i32) % 3 == 0 {
                switch_edge[y as usize][0] = 2 + (y * 5 + 1) % (y - 2);
                switch_edge[y as usize][1] = (y * 3) % p.num_edges;
            }
            y += 1;
        }
        let mut cnt = 0;
        for u in 0..2 * p.num_vertices + 2 {
            let deg = get_vertex_degree(&p, u);
            got.push(format!("D {u} {deg}"));
            for k in 0..deg {
                let mut ie: EdgeIndex = 12345;
                let w = get_vertex_neighbor(&p, u, k, &mut ie);
                let mut s = format!("W {k} {w} {ie}");
                if w != NO_VERTEX {
                    let rc = rescap(&p, u, w, ie);
                    s += &format!(
                        " r{rc} a{}",
                        i32::from(ignore_vertex_non_tacn_atom(&p, u, w))
                    );
                    if u > 1 && w > 1 && switch_edge[u as usize][0] != NO_VERTEX {
                        s += &format!(
                            " g{}",
                            i32::from(ignore_vertex_non_tacn_group(&p, u, w, &switch_edge))
                        );
                        let mut iuv: EdgeIndex = -77;
                        let pv = get_prev_vertex(&p, u, &switch_edge, &mut iuv);
                        s += &format!(" p{pv},{iuv}");
                    }
                    let rm = rescap_mark(&mut p, u, w, ie);
                    s += &format!(" m{rm}");
                    let delta = (cnt % 4) - 1;
                    let rv = ((cnt >> 2) & 1) == 1;
                    let bcf = (cnt >> 3) & 1;
                    let ar = augment_edge(&mut p, u, w, ie, delta, rv, bcf);
                    s += &format!(" A{ar}");
                    cnt += 1;
                }
                got.push(s);
            }
        }
        got.push(format!("S {}", p.b_not_a_simple_path));
        dump_bns(&p, &mut got);
        p.alt_path = None;
        let r = re_init_bn_struct(&mut p, &mut at, n, true);
        got.push(format!("R {r}"));
        dump_bns(&p, &mut got);
        let mut bd = allocate_and_init_bn_data(n as i32 + 1);
        bd.q_size = 2;
        bd.scan_q[0] = 2;
        bd.scan_q[1] = 5;
        bd.scan_q[2] = 3;
        bd.base_ptr[2] = 7;
        bd.base_ptr[3] = 7;
        bd.tree[5] = 2;
        bd.tree[4] = 2;
        bd.switch_edge[2][0] = 9;
        let (mn, ml) = (bd.max_num_vertices, bd.max_len_pu_pv);
        let rr = re_init_bn_data(&mut bd);
        got.push(format!("Z {mn} {ml} {rr} {} {}", bd.q_size, bd.base_ptr[2]));
        got.push(format!(
            "T {} {} {} {}",
            bd.tree[4], bd.tree[5], bd.switch_edge[2][0], bd.base_ptr[3]
        ));
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        if got.len() != want.len() {
            panic!("length {} vs {}: {inp}", got.len(), want.len());
        }
        for (i, (g, w)) in got.iter().zip(&want).enumerate() {
            assert_eq!(g, w, "line {i}: {inp}");
        }
    }
}

#[test]
fn balanced_network_search_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnsearch.json")).unwrap();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        p.type_tacn = inp["tacn"].as_u64().unwrap() as u16;
        p.type_t = inp["t"].as_u64().unwrap() as u16;
        p.type_cn = inp["cn"].as_u64().unwrap() as u16;
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        let mut got: Vec<String> = Vec::new();
        let mut s = String::new();
        for i in 0..n as i32 {
            s += &format!("e{}", i32::from(is_bns_endpoint(&p, i)));
        }
        got.push(s);
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        bd.rad_endpoints = vec![0; 2 * p.max_vertices as usize + 8];
        bd.b_rad_srch_mode = inp["rm"].as_i64().unwrap() as i32;
        let ret = run_balanced_network_search(&mut p, &mut bd, inp["bcf"].as_i64().unwrap() as i32);
        got.push(format!("X {ret} {}", bd.n_num_rad_endpoints));
        let mut s = String::new();
        for i in 0..bd.n_num_rad_endpoints as usize {
            s += &format!(" {}", bd.rad_endpoints[i]);
        }
        got.push(s);
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            got.push(format!(
                "V {} {} {} {} {} {}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, e.pass
            ));
        }
        for e in p.edge.iter().take(p.num_edges as usize) {
            got.push(format!(
                "E {} {} {} {} {}",
                e.cap, e.cap0, e.flow, e.flow0, e.pass
            ));
        }
        got.push(format!(
            "P {} {} {}",
            p.num_altp, p.b_change_flow, p.b_not_a_simple_path
        ));
        for a in 0..p.num_altp.min(p.max_altp) as usize {
            let ap = &p.altp[a];
            let len = ap[I_ALTP_PATH_LEN].number();
            let mut s = format!(
                "A {} {} {} {} {}",
                ap[I_ALTP_START_ATOM].number(),
                ap[I_ALTP_END_ATOM].number(),
                len,
                ap[I_ALTP_FLOW].flow[0],
                ap[I_ALTP_FLOW].flow[1]
            );
            for j in 0..len as usize {
                let e = &ap[I_ALTP_NEIGHBOR + j];
                s += &format!(" {},{}", e.ineigh(0), e.ineigh(1));
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn balanced_network_growth_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnvertex.json")).unwrap();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        let mut got: Vec<String> = Vec::new();
        for o in inp["ops"].as_array().unwrap() {
            let o: Vec<i32> = o
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap() as i32)
                .collect();
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            if (op == 0 && a >= p.num_vertices)
                || ((op == 1 || op == 4) && (a >= p.num_vertices || b >= p.num_vertices))
            {
                got.push("skip".into());
                continue;
            }
            match op {
                0 => {
                    let r = add_new_vertex(&mut p, a, x, y, z, &mut dots);
                    got.push(format!("o0 {r} {dots}"));
                }
                1 => got.push(format!("o1 {}", add_new_edge(&mut p, a, b, x, y))),
                2 => got.push(format!("o2 {}", get_edge_to_group_vertex(&p, a, x as u16))),
                3 => got.push(format!("o3 {}", get_group_vertex(&p, a, x as u16))),
                _ => {
                    let mut buf = [0i32; 64];
                    let r = add_st_cap_to_a_vertex(&mut p, a, b, &mut buf, &mut dots, x != 0);
                    let mut s = format!("o4 {r} {dots}:");
                    for v in &buf[..r as usize] {
                        s += &format!(" {v}");
                    }
                    got.push(s);
                }
            }
            got.push(format!("C {} {}", p.num_vertices, p.num_edges));
            for v in p.vert.iter().take(p.num_vertices as usize) {
                let e = &v.st_edge;
                let mut s = format!(
                    "V {} {} {} {} {} {} {} {} |",
                    v.type_,
                    e.cap,
                    e.cap0,
                    e.flow,
                    e.flow0,
                    e.pass,
                    v.num_adj_edges,
                    v.max_adj_edges
                );
                for k in 0..v.num_adj_edges as usize {
                    s += &format!(" {}", p.iedge[v.iedge + k]);
                }
                got.push(s);
            }
            for e in p.edge.iter().take(p.num_edges as usize) {
                got.push(format!(
                    "E {} {} {} {} {} {} {} {} {}",
                    e.cap,
                    e.cap0,
                    e.flow,
                    e.flow0,
                    e.pass,
                    e.neighbor1,
                    e.neighbor12,
                    e.neigh_ord[0],
                    e.neigh_ord[1]
                ));
            }
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn alternating_path_check_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnaltpath.json")).unwrap();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        let mut got: Vec<String> = Vec::new();
        let mut fcd = vec![BnsFlowChanges::default(); 64];
        let mut apc = AltPathChanges::default();
        p.pb_taut_flags = inp["tf"].as_u64().unwrap() as u32;
        for o in inp["ops"].as_array().unwrap() {
            let o: Vec<i32> = o
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap() as i32)
                .collect();
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            let nv = p.num_vertices;
            match op {
                0 => {
                    if p.num_edges == 0 {
                        got.push("skip".into());
                        continue;
                    }
                    let ie = a % p.num_edges;
                    let r = set_flow_to_check_one_bond(&mut p, ie, x % 3, &mut fcd);
                    got.push(format!("o0 {r}"));
                }
                1 => {
                    got.push(format!(
                        "o1 {}",
                        restore_flow_after_check_one_bond(&mut p, &fcd)
                    ));
                    fcd[0].iedge = NO_VERTEX;
                }
                2 => {
                    if a >= nv || b >= nv {
                        got.push("skip".into());
                        continue;
                    }
                    let r = set_bns_to_check_alt_path(
                        &mut p, a, b, x as u16, y, &mut apc, &mut fcd, &mut dots,
                    );
                    got.push(format!("o2 {r} {dots}"));
                }
                3 => {
                    got.push(format!(
                        "o3 {}",
                        restore_bns_after_check_alt_path(&mut p, &apc, x)
                    ));
                    apc = AltPathChanges::default();
                }
                10 => {
                    if a >= nv {
                        got.push("skip".into());
                        continue;
                    }
                    let r = add_new_vertex(&mut p, a, x, y, z, &mut dots);
                    got.push(format!("o10 {r} {dots}"));
                }
                11 => {
                    if a < p.num_atoms || a >= nv {
                        got.push("skip".into());
                        continue;
                    }
                    p.vert[a as usize].type_ = x as u16;
                    got.push("o11".into());
                }
                _ => {
                    if a >= nv || b >= nv || a == b {
                        got.push("skip".into());
                        continue;
                    }
                    got.push(format!("o12 {}", add_new_edge(&mut p, a, b, x, y)));
                }
            }
            got.push(format!("C {} {}", p.num_vertices, p.num_edges));
            for v in p.vert.iter().take(p.num_vertices as usize) {
                let e = &v.st_edge;
                let mut s = format!(
                    "V {} {} {} {} {} {} {} {} |",
                    v.type_,
                    e.cap,
                    e.cap0,
                    e.flow,
                    e.flow0,
                    e.pass,
                    v.num_adj_edges,
                    v.max_adj_edges
                );
                for k in 0..v.num_adj_edges as usize {
                    s += &format!(" {}", p.iedge[v.iedge + k]);
                }
                got.push(s);
            }
            for e in p.edge.iter().take(p.num_edges as usize) {
                got.push(format!(
                    "E {} {} {} {} {} {} {} {} {}",
                    e.cap,
                    e.cap0,
                    e.flow,
                    e.flow0,
                    e.pass,
                    e.neighbor1,
                    e.neighbor12,
                    e.neigh_ord[0],
                    e.neigh_ord[1]
                ));
            }
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn flow_to_bond_transfer_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnbonds.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for (i, a) in inp["atoms"].as_array().unwrap().iter().enumerate() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                n_ring_system: ((i * 7 + ci) % 3) as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        let mode = inp["mode"].as_i64().unwrap();
        let (bcf, bcf2) = (
            inp["bcf"].as_i64().unwrap() as i32,
            inp["bcf2"].as_i64().unwrap() as i32,
        );
        let nt = inp["nt"].as_i64().unwrap() as i32;
        let mut got: Vec<String> = Vec::new();
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let mut fcd = vec![BnsFlowChanges::default(); 64];
        if mode == 2 && p.num_edges > 0 {
            let ie = inp["ei"].as_i64().unwrap() as i32 % p.num_edges;
            got.push(format!(
                "f {}",
                set_flow_to_check_one_bond(&mut p, ie, nt % 3, &mut fcd)
            ));
        }
        let ret = run_balanced_network_search(&mut p, &mut bd, bcf);
        got.push(format!("X {ret}"));
        match mode {
            0 => got.push(format!(
                "r {}",
                set_bonds_from_bn_struct_flow(&mut p, &mut at, n as i32, bcf2)
            )),
            1 => got.push(format!("r {}", restore_bn_struct_flow(&mut p, bcf2))),
            _ => {
                let r1 =
                    set_bonds_after_check_one_bond(&mut p, &fcd, nt - 1, &mut at, n as i32, bcf2);
                let r2 = set_bonds_from_bn_struct_flow(&mut p, &mut at, n as i32, bcf2);
                let r3 = restore_bn_struct_flow(&mut p, bcf2 & 3);
                restore_flow_after_check_one_bond(&mut p, &fcd);
                got.push(format!("r {r1} {r2} {r3}"));
            }
        }
        for a in &at {
            let mut s = format!("T {} {} {} |", a.valence, a.chem_bonds_valence, a.radical);
            for j in 0..a.valence as usize {
                s += &format!(" {}", a.bond_type[j]);
            }
            got.push(s);
        }
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            got.push(format!(
                "V {} {} {} {} {} {}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, e.pass
            ));
        }
        for e in p.edge.iter().take(p.num_edges as usize) {
            got.push(format!(
                "E {} {} {} {} {}",
                e.cap, e.cap0, e.flow, e.flow0, e.pass
            ));
        }
        got.push(format!(
            "P {} {} {}",
            p.num_altp, p.b_change_flow, p.b_not_a_simple_path
        ));
        for a in 0..p.num_altp.min(p.max_altp) as usize {
            let ap = &p.altp[a];
            let len = ap[I_ALTP_PATH_LEN].number();
            let mut s = format!(
                "A {} {} {} {} {}",
                ap[I_ALTP_START_ATOM].number(),
                ap[I_ALTP_END_ATOM].number(),
                len,
                ap[I_ALTP_FLOW].flow[0],
                ap[I_ALTP_FLOW].flow[1]
            );
            for j in 0..len as usize {
                let e = &ap[I_ALTP_NEIGHBOR + j];
                s += &format!(" {},{}", e.ineigh(0), e.ineigh(1));
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn radical_endpoints_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn dump_all(p: &BnStruct, bd: &BnData, at: &[InpAtom], got: &mut Vec<String>) {
        got.push(format!(
            "C {} {} {} {}",
            p.num_vertices, p.num_edges, p.tot_st_cap, p.tot_st_flow
        ));
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            let mut s = format!(
                "V {} {} {} {} {} {} {} {} |",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, e.pass, v.num_adj_edges, v.max_adj_edges
            );
            for k in 0..v.num_adj_edges as usize {
                s += &format!(" {}", p.iedge[v.iedge + k]);
            }
            got.push(s);
        }
        for e in p.edge.iter().take(p.num_edges as usize) {
            got.push(format!(
                "E {} {} {} {} {} {} {} {} {}",
                e.cap,
                e.cap0,
                e.flow,
                e.flow0,
                e.pass,
                e.neighbor1,
                e.neighbor12,
                e.neigh_ord[0],
                e.neigh_ord[1]
            ));
        }
        let mut s = format!(
            "R {} {} {} |",
            bd.n_num_rad_endpoints, bd.n_num_rad_edges, bd.n_num_radicals
        );
        for i in 0..bd.n_num_rad_endpoints as usize {
            s += &format!(" {}", bd.rad_endpoints[i]);
        }
        s += " |";
        for i in 0..bd.n_num_rad_edges as usize {
            s += &format!(" {}", bd.rad_edges[i]);
        }
        got.push(s);
        let mut s = String::from("Q");
        for a in at {
            s += &format!(" {}", a.radical);
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnrad.json")).unwrap();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let len = 2 * p.max_vertices as usize + 8;
        bd.rad_endpoints = vec![0; len];
        bd.rad_edges = vec![0; len];
        let rm = inp["rm"].as_i64().unwrap() as i32;
        let ret = if inp["which"].as_i64().unwrap() == 1 {
            set_rad_endpoints(&mut p, &mut bd, rm)
        } else {
            set_rad_endpoints2(&mut p, &mut bd, rm)
        };
        let mut got = vec![format!("X {ret}")];
        dump_all(&p, &bd, &at, &mut got);
        let with_at = inp["withat"].as_i64().unwrap() == 1;
        if inp["rf"].as_i64().unwrap() == 1 {
            let r = restore_radicals_only(&p, &bd, if with_at { Some(&mut at) } else { None });
            got.push(format!("S {r}"));
            dump_all(&p, &bd, &at, &mut got);
        }
        let r = remove_rad_endpoints(&mut p, &mut bd, if with_at { Some(&mut at) } else { None });
        got.push(format!("Y {r}"));
        dump_all(&p, &bd, &at, &mut got);
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn atom_charge_types_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_atomtype.json")).unwrap();
    let mut types = std::collections::HashSet::new();
    for c in data.as_array().unwrap() {
        let mut at: Vec<InpAtom> = Vec::new();
        for a in c["in"].as_array().unwrap() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                radical: a["rad"].as_i64().unwrap() as i8,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let mut got: Vec<String> = Vec::new();
        for sub in [0, 1, 2, -2] {
            let mut totals = [0i32; attot::ARRAY_LEN];
            let mut s = String::new();
            for i in 0..n {
                let (t, m) = if i & 1 == 1 {
                    get_atom_charge_type(&at, i, None, sub)
                } else {
                    get_atom_charge_type(&at, i, Some(&mut totals), sub)
                };
                types.insert(t);
                s += &format!("{},{} ", t, m.unwrap_or(-12345));
            }
            s += "|";
            for t in totals {
                s += &format!(" {t}");
            }
            got.push(s);
        }
        for i in 0..n {
            let mut s = format!("z{}", i32::from(is_z_atom(at[i].el_number)));
            if at[i].valence != 0 {
                s += &format!(" x{}", is_zox(&at, i, 0));
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{}", c["in"]);
    }
    assert!(types.len() > 12, "only {} types seen", types.len());
}

#[test]
fn flow_changes_of_charges_and_hydrogens_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::TGroupInfo;
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_atoms(at: &[InpAtom], t: &TGroupInfo, got: &mut Vec<String>) {
        for a in at {
            got.push(format!(
                "A {} {} {} {} {},{},{} t{} s{},{},{} o{},{},{} n{},{},{} g{},{},{} r{} nb{}",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.sb_parity[0],
                a.sb_parity[1],
                a.sb_parity[2],
                a.sb_ord[0],
                a.sb_ord[1],
                a.sb_ord[2],
                a.sn_ord[0],
                a.sn_ord[1],
                a.sn_ord[2],
                a.sn_orig_at_num[0],
                a.sn_orig_at_num[1],
                a.sn_orig_at_num[2],
                a.orig_at_number,
                a.neighbor[0]
            ));
        }
        got.push(format!(
            "I {} {} {},{},{}",
            t.tni.n_num_removed_explicit_h,
            t.tni.n_num_removed_protons,
            t.tni.n_num_removed_protons_isotopic[0],
            t.tni.n_num_removed_protons_isotopic[1],
            t.tni.n_num_removed_protons_isotopic[2]
        ));
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnchargeh.json")).unwrap();
    let (mut n_marked, mut n_changed) = (0, 0);
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let (sbp, sbo, sno, sng) = (
                ints(&a["sbp"]),
                ints(&a["sbo"]),
                ints(&a["sno"]),
                ints(&a["sng"]),
            );
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
                x.sb_parity[k] = sbp[k] as i8;
                x.sb_ord[k] = sbo[k] as i8;
                x.sn_ord[k] = sno[k] as i8;
                x.sn_orig_at_num[k] = sng[k] as u16;
            }
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            let nv = p.num_vertices;
            if op == 10 {
                if a >= nv {
                    continue;
                }
                add_new_vertex(&mut p, a, x, y, z, &mut dots);
            } else if op == 11 {
                if a < p.num_atoms || a >= nv {
                    continue;
                }
                p.vert[a as usize].type_ = x as u16;
            } else {
                if a >= nv || b >= nv || a == b {
                    continue;
                }
                add_new_edge(&mut p, a, b, x, y);
            }
        }
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let mut ret = 0;
        if !inp["walk"].as_bool().unwrap() {
            ret = run_balanced_network_search(&mut p, &mut bd, inp["bcf"].as_i64().unwrap() as i32);
        } else {
            // a random walk along the edges stands for an alternating path
            let mut lcg = inp["seed"].as_u64().unwrap() as u32;
            let mut rnd = |m: u32| -> i32 {
                lcg = lcg.wrapping_mul(1103515245).wrapping_add(12345) & 0x7fff_ffff;
                ((lcg >> 8) % m) as i32
            };
            p.num_altp = 0;
            let np = inp["np"].as_i64().unwrap() as i32;
            for pass in 0..np.min(p.max_altp) {
                let mut v = -1;
                for _ in 0..20 {
                    let c2 = rnd(p.num_vertices as u32);
                    if p.vert[c2 as usize].num_adj_edges > 0 {
                        v = c2;
                        break;
                    }
                }
                if v < 0 {
                    continue;
                }
                let l = (1 + rnd(10)).min(p.max_vertices);
                let start = v;
                let delta = if rnd(2) != 0 { 1 } else { -1 };
                let mut k = 0usize;
                while (k as i32) < l {
                    let pv = p.vert[v as usize];
                    if pv.num_adj_edges == 0 {
                        break;
                    }
                    let o = rnd(u32::from(pv.num_adj_edges));
                    let e = p.edge[p.iedge[pv.iedge + o as usize] as usize];
                    let v2 = i32::from(e.neighbor12) ^ v;
                    let o2 = if v2 == i32::from(e.neighbor1) {
                        e.neigh_ord[0]
                    } else {
                        e.neigh_ord[1]
                    };
                    p.altp[pass as usize][I_ALTP_NEIGHBOR + k].set_ineigh(0, o as u16);
                    p.altp[pass as usize][I_ALTP_NEIGHBOR + k].set_ineigh(1, o2);
                    v = v2;
                    k += 1;
                }
                if k == 0 {
                    continue;
                }
                let ap = &mut p.altp[pass as usize];
                ap[I_ALTP_START_ATOM].set_number(start);
                ap[I_ALTP_END_ATOM].set_number(v);
                ap[I_ALTP_PATH_LEN].set_number(k as i32);
                ap[I_ALTP_FLOW].flow[0] = delta;
                p.num_altp += 1;
            }
        }
        let mut got = vec![format!("X {ret} {} {}", p.num_altp, p.num_vertices)];
        for a in 0..p.num_altp.min(p.max_altp) as usize {
            let ap = &p.altp[a];
            let len = ap[I_ALTP_PATH_LEN].number();
            let mut s = format!(
                "P {} {} {} {}",
                ap[I_ALTP_START_ATOM].number(),
                ap[I_ALTP_END_ATOM].number(),
                len,
                ap[I_ALTP_FLOW].flow[0]
            );
            for j in 0..len as usize {
                let e = &ap[I_ALTP_NEIGHBOR + j];
                s += &format!(" {},{}", e.ineigh(0), e.ineigh(1));
            }
            got.push(s);
        }
        let mkseed = inp["mk"].as_i64().unwrap();
        let mut aatg = BnAatg {
            n_num_found: 0,
            n_at_type_totals: None,
            t_group_info: None,
            n_marked_atom: (0..2 * p.max_vertices as usize + 8)
                .map(|i| {
                    if i < n {
                        ((mkseed >> (i % 8)) & 1) as i8
                    } else {
                        0
                    }
                })
                .collect(),
        };
        let r1 = mark_atoms_at_taut_groups(
            &mut p,
            n as i32,
            &mut aatg,
            inp["e1"].as_i64().unwrap() as i32,
            inp["e2"].as_i64().unwrap() as i32,
        );
        let mut s = format!("M {r1} {} |", aatg.n_num_found);
        for m in &aatg.n_marked_atom[..p.num_vertices as usize] {
            s += &format!(" {m}");
        }
        got.push(s);
        n_marked += i32::from(aatg.n_num_found > 0);
        got.push(format!(
            "E {}",
            eliminate_plus_minus_charge_ambiguity(&mut p, n as i32)
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            s += &format!(" {}:{}:{}", v.type_, v.st_edge.cap, v.st_edge.flow);
        }
        got.push(s);
        let mut s = String::from("F");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(" {}", e.flow);
        }
        got.push(s);
        let mut tgi = TGroupInfo::default();
        tgi.tni.n_num_removed_explicit_h = hs.len() as i32;
        tgi.tni.n_num_removed_protons = 1;
        let mut mark2 = vec![0i8; n + 1];
        let mut totals = [0i32; attot::ARRAY_LEN];
        let r3 = subtract_or_change_at_h_charge_bns(
            &mut p,
            &mut at,
            n as i32,
            Some(&mut totals),
            Some(&mut mark2),
            None,
            true,
        );
        let mut s = format!("S {r3} |");
        for m in &mark2[..n] {
            s += &format!(" {m}");
        }
        s += " |";
        for t in totals {
            s += &format!(" {t}");
        }
        got.push(s);
        let r4 = add_changed_at_h_charge_bns(&mut at, n as i32, Some(&mut totals), &mut mark2);
        let mut s = format!("D {r4} |");
        for t in totals {
            s += &format!(" {t}");
        }
        got.push(s);
        let r5 = subtract_or_change_at_h_charge_bns(
            &mut p,
            &mut at,
            n as i32,
            None,
            None,
            Some(&mut tgi),
            false,
        );
        got.push(format!("C {r5}"));
        n_changed += i32::from(r5 > 0);
        dump_atoms(&at, &tgi, &mut got);
        for hc in inp["hc"].as_array().unwrap() {
            let hc = ints(hc);
            let r6 =
                add_or_remove_expl_or_impl_h(hc[0], &mut at, n as i32, hc[1] as usize, &mut tgi);
            got.push(format!("H {r6}"));
            dump_atoms(&at, &tgi, &mut got);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(n_marked > 0 && n_changed > 100, "{n_marked} {n_changed}");
}

#[test]
fn alternating_path_existence_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::TGroupInfo;
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_state(p: &BnStruct, at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        got.push(format!(
            "C {} {} {} {}",
            p.num_vertices, p.num_edges, p.tot_st_cap, p.tot_st_flow
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges
            );
        }
        got.push(s);
        let mut s = String::from("F");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(" {}:{}:{}:{}", e.cap, e.cap0, e.flow, e.flow0);
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at.iter().take(n) {
            s += &format!(
                " {},{},{},{},{},{},{},{},{}:",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.endpoint
            );
            for j in 0..a.valence as usize {
                s += &format!("{}.", a.bond_type[j]);
            }
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnexists.json")).unwrap();
    let mut seen = std::collections::HashMap::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                at_type: a["att"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
            }
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        for (i, a) in at.iter_mut().enumerate() {
            a.n_ring_system = (i % 2) as u16;
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        p.pb_taut_flags = inp["tf"].as_u64().unwrap() as u32;
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            let nv = p.num_vertices;
            if op == 10 {
                if a >= nv {
                    continue;
                }
                add_new_vertex(&mut p, a, x, y, z, &mut dots);
            } else if op == 11 {
                if a < p.num_atoms || a >= nv {
                    continue;
                }
                p.vert[a as usize].type_ = x as u16;
            } else {
                if a >= nv || b >= nv || a == b {
                    continue;
                }
                add_new_edge(&mut p, a, b, x, y);
            }
        }
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let len = 2 * p.max_vertices as usize + 8;
        bd.rad_endpoints = vec![0; len];
        bd.rad_edges = vec![0; len];
        let mut got: Vec<String> = Vec::new();
        dump_state(&p, &at, n, &mut got);
        for cl in inp["calls"].as_array().unwrap() {
            let cl = ints(cl);
            let (vd, vs, pt, any, am) = (cl[0], cl[1], cl[2], cl[3], cl[4]);
            let mut aatg = BnAatg::default();
            let mut tgi = TGroupInfo::default();
            tgi.tni.n_num_removed_explicit_h = hs.len() as i32;
            if am > 0 {
                aatg.n_marked_atom = vec![0; len];
            }
            if am > 1 {
                aatg.n_at_type_totals = Some(vec![0; attot::ARRAY_LEN]);
                aatg.t_group_info = Some(tgi.clone());
            }
            let ret = if any != 0 {
                exists_any_alt_path(&mut p, &mut bd, &mut at, n as i32, vd, vs, pt)
            } else {
                exists_alt_path(
                    &mut p,
                    &mut bd,
                    if am != 0 { Some(&mut aatg) } else { None },
                    &mut at,
                    n as i32,
                    vd,
                    vs,
                    pt,
                )
            };
            *seen.entry(ret).or_insert(0) += 1;
            let protons = aatg
                .t_group_info
                .as_ref()
                .map_or(tgi.tni.n_num_removed_protons, |t| {
                    t.tni.n_num_removed_protons
                });
            got.push(format!("R {ret} {} {protons}", aatg.n_num_found));
            if am > 0 {
                let mut s = String::from("K");
                for m in &aatg.n_marked_atom[..n] {
                    s += &format!(" {m}");
                }
                got.push(s);
            }
            if am > 1 {
                let mut s = String::from("T");
                for t in aatg.n_at_type_totals.as_ref().unwrap() {
                    s += &format!(" {t}");
                }
                got.push(s);
            }
            dump_state(&p, &at, n, &mut got);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(seen.keys().any(|&k| k > 0) && seen.len() > 3, "{seen:?}");
}

#[test]
fn group_vertices_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::{CGroup, CGroupInfo, TGroup, TGroupInfo};
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_state(
        p: &BnStruct,
        at: &[InpAtom],
        t: &TGroupInfo,
        cg: &CGroupInfo,
        got: &mut Vec<String>,
    ) {
        got.push(format!(
            "C {} {} {} {} {} {}",
            p.num_vertices,
            p.num_edges,
            p.num_t_groups,
            p.num_c_groups,
            p.num_added_atoms,
            p.num_added_edges
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}:{}:{}|",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges, v.max_adj_edges, v.iedge
            );
            for k in 0..v.num_adj_edges as usize {
                s += &format!("{}.", p.iedge[v.iedge + k]);
            }
        }
        got.push(s);
        let mut s = String::from("E");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(
                " {}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
                e.cap,
                e.cap0,
                e.flow,
                e.flow0,
                e.pass,
                e.forbidden,
                e.neighbor1,
                e.neighbor12,
                e.neigh_ord[0],
                e.neigh_ord[1]
            );
        }
        got.push(s);
        let mut s = String::from("G");
        for g in t.t_group.iter().flatten().take(t.num_t_groups) {
            s += &format!(" {}", g.n_group_number);
        }
        s += " |";
        for g in cg.c_group.iter().flatten().take(cg.num_c_groups) {
            s += &format!(" {}", g.n_group_number);
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at {
            s += &format!(" {},{}", a.endpoint, a.c_point);
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bngroups.json")).unwrap();
    let mut codes = std::collections::HashSet::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                c_point: a["cp"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        p.pb_taut_flags = inp["tfl"].as_u64().unwrap() as u32;
        let mut tgi = TGroupInfo {
            t_group: Some(
                inp["tg"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|g| {
                        let g = ints(g);
                        TGroup {
                            n_group_number: g[0] as u16,
                            n_num_endpoints: g[1] as u16,
                            ..Default::default()
                        }
                    })
                    .collect(),
            ),
            num_t_groups: inp["tg"].as_array().unwrap().len(),
            b_taut_flags: inp["tf"].as_u64().unwrap(),
            ..Default::default()
        };
        let mut cgi = CGroupInfo {
            c_group: Some(
                inp["cg"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|g| {
                        let g = ints(g);
                        CGroup {
                            n_group_number: g[0] as u16,
                            num_c_points: g[1] as u16,
                            ..Default::default()
                        }
                    })
                    .collect(),
            ),
            num_c_groups: inp["cg"].as_array().unwrap().len(),
            ..Default::default()
        };
        let mode = inp["mode"].as_i64().unwrap();
        let mut got: Vec<String> = Vec::new();
        for _ in 0..2 {
            let ret = if mode == 0 {
                let has_t = tgi.num_t_groups > 0;
                let has_c = cgi.num_c_groups > 0;
                re_init_bn_struct_add_groups(
                    &mut p,
                    &mut at,
                    n,
                    if has_t { Some(&mut tgi) } else { None },
                    if has_c { Some(&mut cgi) } else { None },
                )
            } else {
                let r0 = re_init_bn_struct(&mut p, &mut at, n, false);
                let r = if mode == 1 {
                    add_t_groups_to_bn_struct(&mut p, &at, n as i32, Some(&mut tgi))
                } else {
                    add_c_groups_to_bn_struct(&mut p, &at, n as i32, Some(&mut cgi))
                };
                r + 100000 * r0
            };
            codes.insert(ret);
            got.push(format!("R {ret}"));
            dump_state(&p, &at, &tgi, &cgi, &mut got);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(codes.len() >= 3, "{codes:?}");
}

#[test]
fn single_group_creation_and_removal_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn dump_state(p: &BnStruct, at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        got.push(format!(
            "C {} {} {} {} {} {}",
            p.num_vertices,
            p.num_edges,
            p.num_t_groups,
            p.num_c_groups,
            p.num_added_atoms,
            p.num_added_edges
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}:{}:{}|",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges, v.max_adj_edges, v.iedge
            );
            for k in 0..v.num_adj_edges as usize {
                s += &format!("{}.", p.iedge[v.iedge + k]);
            }
        }
        got.push(s);
        let mut s = String::from("E");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(
                " {}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
                e.cap,
                e.cap0,
                e.flow,
                e.flow0,
                e.pass,
                e.forbidden,
                e.neighbor1,
                e.neighbor12,
                e.neigh_ord[0],
                e.neigh_ord[1]
            );
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at.iter().take(n) {
            s += &format!(" {},{}", a.endpoint, a.c_point);
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bncreate.json")).unwrap();
    let (mut created, mut removed) = (0, 0);
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                c_point: a["cp"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        p.type_tacn = inp["tacn"].as_u64().unwrap() as u16;
        let mut got: Vec<String> = Vec::new();
        dump_state(&p, &at, n, &mut got);
        for o in inp["ops"].as_array().unwrap() {
            let o: Vec<i32> = o
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap() as i32)
                .collect();
            let (op, ty, mask, chg) = (o[0], o[1], o[2], o[3]);
            let r = match op {
                0 => create_c_group_in_bn_struct(&at, n as i32, &mut p, ty, mask, chg),
                1 => create_t_group_in_bn_struct(&at, n as i32, &mut p, ty, mask),
                2 => {
                    let tg = if ty < 0 {
                        if p.num_vertices > p.num_atoms {
                            p.num_vertices - 1
                        } else {
                            1000
                        }
                    } else {
                        ty
                    };
                    remove_last_group_from_bn_struct(&mut at, n as i32, tg, &mut p)
                }
                _ => set_init_cap_flow_to_current(&mut p),
            };
            if op < 2 && r > 0 {
                created += 1;
            }
            if op == 2 && r == 0 {
                removed += 1;
            }
            got.push(format!("R {r}"));
            dump_state(&p, &at, n, &mut got);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(created > 100 && removed > 50, "{created} {removed}");
}

#[test]
fn simple_proton_moves_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::TGroupInfo;
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump(at: &[InpAtom], t: &TGroupInfo, totals: &[i32], got: &mut Vec<String>) {
        for a in at {
            got.push(format!(
                "A {} {} {} {} {},{},{} t{} s{},{},{} o{},{},{} n{},{},{} g{},{},{} r{} nb{}",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.sb_parity[0],
                a.sb_parity[1],
                a.sb_parity[2],
                a.sb_ord[0],
                a.sb_ord[1],
                a.sb_ord[2],
                a.sn_ord[0],
                a.sn_ord[1],
                a.sn_ord[2],
                a.sn_orig_at_num[0],
                a.sn_orig_at_num[1],
                a.sn_orig_at_num[2],
                a.orig_at_number,
                a.neighbor[0]
            ));
        }
        got.push(format!(
            "I {} {} {},{},{}",
            t.tni.n_num_removed_explicit_h,
            t.tni.n_num_removed_protons,
            t.tni.n_num_removed_protons_isotopic[0],
            t.tni.n_num_removed_protons_isotopic[1],
            t.tni.n_num_removed_protons_isotopic[2]
        ));
        let tot: Vec<String> = totals.iter().map(|x| format!(" {x}")).collect();
        got.push(format!("T{}", tot.concat()));
    }
    fn queries(at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        for i in 0..n {
            let mut sub = 0;
            let r = b_is_hard_rem_h_candidate(at, i, &mut sub);
            got.push(format!(
                "Q {r} {sub} {} {} {}",
                b_has_acidic_hydrogen(at, i),
                b_has_other_exchangable_h(at, i),
                b_has_acidic_minus(at, i)
            ));
        }
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnsimple.json")).unwrap();
    let mut moved = [0; 3];
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let (sbp, sbo, sno, sng) = (
                ints(&a["sbp"]),
                ints(&a["sbo"]),
                ints(&a["sno"]),
                ints(&a["sng"]),
            );
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
                x.sb_parity[k] = sbp[k] as i8;
                x.sb_ord[k] = sbo[k] as i8;
                x.sn_ord[k] = sno[k] as i8;
                x.sn_orig_at_num[k] = sng[k] as u16;
            }
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let mut tgi = TGroupInfo::default();
        tgi.tni.n_num_removed_explicit_h = hs.len() as i32;
        tgi.tni.n_num_removed_protons = 1;
        let mut totals = vec![0i32; attot::ARRAY_LEN];
        for i in 0..n {
            get_atom_charge_type(&at, i, Some(&mut totals), 0);
        }
        let mut aatg = BnAatg {
            n_at_type_totals: Some(totals.clone()),
            t_group_info: Some(tgi.clone()),
            ..Default::default()
        };
        let mut got: Vec<String> = Vec::new();
        queries(&at, n, &mut got);
        dump(&at, &tgi, &totals, &mut got);
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let (op, a) = (o[0], o[1]);
            let r = match op {
                0 => {
                    let t = if a != 0 { Some(&mut totals[..]) } else { None };
                    simple_remove_hplus_npo(&mut at, n as i32, t, &mut tgi)
                }
                1 | 2 => {
                    // the C code works on one structure for both the totals
                    // and the group info
                    aatg.n_at_type_totals = Some(totals.clone());
                    aatg.t_group_info = Some(tgi.clone());
                    let r = if op == 1 {
                        simple_remove_acidic_protons(&mut at, n as i32, &mut aatg, a)
                    } else {
                        simple_add_acidic_protons(&mut at, n as i32, &mut aatg, a)
                    };
                    totals = aatg.n_at_type_totals.clone().unwrap();
                    tgi = aatg.t_group_info.clone().unwrap();
                    r
                }
                _ => unreachable!(),
            };
            moved[op as usize] += r;
            got.push(format!("R {op} {r}"));
            queries(&at, n, &mut got);
            dump(&at, &tgi, &totals, &mut got);
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(moved.iter().all(|&m| m > 30), "{moved:?}");
}

#[test]
fn hard_proton_moves_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::TGroupInfo;
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_state(p: &BnStruct, at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        got.push(format!(
            "C {} {} {} {}",
            p.num_vertices, p.num_edges, p.tot_st_cap, p.tot_st_flow
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges
            );
        }
        got.push(s);
        let mut s = String::from("F");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(" {}:{}:{}:{}", e.cap, e.cap0, e.flow, e.flow0);
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at.iter().take(n) {
            s += &format!(
                " {},{},{},{},{},{},{},{},{}:",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.endpoint
            );
            for j in 0..a.valence as usize {
                s += &format!("{}.", a.bond_type[j]);
            }
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnhard.json")).unwrap();
    let mut seen = std::collections::HashMap::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                at_type: a["att"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
            }
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        for (i, a) in at.iter_mut().enumerate() {
            a.n_ring_system = (i % 2) as u16;
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        p.pb_taut_flags = inp["tf"].as_u64().unwrap() as u32;
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            let nv = p.num_vertices;
            if op == 10 {
                if a >= nv {
                    continue;
                }
                add_new_vertex(&mut p, a, x, y, z, &mut dots);
            } else if op == 11 {
                if a < p.num_atoms || a >= nv {
                    continue;
                }
                p.vert[a as usize].type_ = x as u16;
            } else {
                if a >= nv || b >= nv || a == b {
                    continue;
                }
                add_new_edge(&mut p, a, b, x, y);
            }
        }
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let len = 2 * p.max_vertices as usize + 8;
        bd.rad_endpoints = vec![0; len];
        bd.rad_edges = vec![0; len];
        let mut got: Vec<String> = Vec::new();
        dump_state(&p, &at, n, &mut got);
        let mut tgi = TGroupInfo::default();
        tgi.tni.n_num_removed_explicit_h = hs.len() as i32;
        tgi.tni.n_num_removed_protons = 1;
        tgi.b_taut_flags = u64::from(p.pb_taut_flags);
        let mut totals = vec![0; attot::ARRAY_LEN];
        let mut aatg = BnAatg {
            n_marked_atom: vec![0; len],
            ..Default::default()
        };
        got.push(format!(
            "M {}",
            mark_at_type(&mut at, n as i32, Some(&mut totals))
        ));
        aatg.n_at_type_totals = Some(totals);
        aatg.t_group_info = Some(tgi);
        for cl in inp["calls"].as_array().unwrap() {
            let cl = ints(cl);
            let (f, arg) = (cl[0], cl[1]);
            let mut canceled = 0;
            let ni = n as i32;
            let ret = match f {
                0 => hard_remove_acidic_protons(
                    &mut at,
                    ni,
                    &mut aatg,
                    arg,
                    Some(&mut canceled),
                    &mut p,
                    &mut bd,
                ),
                1 => hard_add_acidic_protons(
                    &mut at,
                    ni,
                    &mut aatg,
                    arg,
                    Some(&mut canceled),
                    &mut p,
                    &mut bd,
                ),
                2 => hard_remove_hplus_np(
                    &mut at,
                    ni,
                    arg != 0,
                    Some(&mut canceled),
                    &mut aatg,
                    &mut p,
                    &mut bd,
                ),
                _ => remove_np_protons_and_acid_charges(&mut at, ni, &mut aatg, &mut p, &mut bd),
            };
            *seen.entry(ret).or_insert(0) += 1;
            let t = aatg.t_group_info.as_ref().unwrap();
            got.push(format!(
                "R {ret} {canceled} {} {} {}",
                t.tni.b_normalization_flags,
                t.tni.n_num_removed_protons,
                t.tni.n_num_removed_explicit_h
            ));
            let mut s = String::from("T");
            for x in aatg.n_at_type_totals.as_ref().unwrap() {
                s += &format!(" {x}");
            }
            got.push(s);
            dump_state(&p, &at, n, &mut got);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(seen.keys().any(|&k| k > 0) && seen.len() > 3, "{seen:?}");
}

#[test]
fn alternating_bond_marking_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnaltmark.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for (i, a) in inp["atoms"].as_array().unwrap().iter().enumerate() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                n_ring_system: ((i * 7 + ci) % 3) as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        let mode = inp["mode"].as_i64().unwrap();
        let (bcf, bcf2) = (
            inp["bcf"].as_i64().unwrap() as i32,
            inp["bcf2"].as_i64().unwrap() as i32,
        );
        let nt = inp["nt"].as_i64().unwrap() as i32;
        let mut got: Vec<String> = Vec::new();
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let mut fcd = vec![BnsFlowChanges::default(); 64];
        let ni = n as i32;
        let _ = bcf;
        got.push(format!(
            "a {}",
            bns_adjust_flow_bonds_rad(&mut p, &mut bd, &mut at, ni)
        ));
        if mode == 1 {
            for i in 0..p.num_edges {
                got.push(format!(
                    "m {} {} {}",
                    n_min_flow2_check(&p, i),
                    n_max_flow2_check(&p, i),
                    n_cur_flow2_check(&p, i)
                ));
            }
            let mut s = String::from("b");
            for a in &at {
                for j in 0..a.valence as usize {
                    for f in 0..3 {
                        let bt = i32::from(a.bond_type[j]);
                        s += &format!(
                            " {}{}",
                            i32::from(b_need_to_test_the_flow(bt, f, false)),
                            i32::from(b_need_to_test_the_flow(bt, f, true))
                        );
                    }
                }
            }
            got.push(s);
            let r1 = bns_test_and_mark_alt_bonds(&mut p, &mut bd, &mut at, ni, &mut fcd, bcf2, nt);
            got.push(format!("r {r1}"));
        }
        for a in &at {
            let (bv, na, nw) = n_bonds_valence_inp_at(a);
            got.push(format!("v {bv} {na} {nw}"));
        }
        for a in &at {
            let mut s = format!("T {} {} {} |", a.valence, a.chem_bonds_valence, a.radical);
            for j in 0..a.valence as usize {
                s += &format!(" {}", a.bond_type[j]);
            }
            got.push(s);
        }
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            got.push(format!(
                "V {} {} {} {} {} {}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, e.pass
            ));
        }
        for e in p.edge.iter().take(p.num_edges as usize) {
            got.push(format!(
                "E {} {} {} {} {}",
                e.cap, e.cap0, e.flow, e.flow0, e.pass
            ));
        }
        got.push(format!(
            "P {} {} {}",
            p.num_altp, p.b_change_flow, p.b_not_a_simple_path
        ));
        for a in 0..p.num_altp.min(p.max_altp) as usize {
            let ap = &p.altp[a];
            let len = ap[I_ALTP_PATH_LEN].number();
            let mut s = format!(
                "A {} {} {} {} {}",
                ap[I_ALTP_START_ATOM].number(),
                ap[I_ALTP_END_ATOM].number(),
                len,
                ap[I_ALTP_FLOW].flow[0],
                ap[I_ALTP_FLOW].flow[1]
            );
            for j in 0..len as usize {
                let e = &ap[I_ALTP_NEIGHBOR + j];
                s += &format!(" {},{}", e.ineigh(0), e.ineigh(1));
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn forbidden_edges_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnforbid.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for (i, a) in inp["atoms"].as_array().unwrap().iter().enumerate() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                n_ring_system: ((i * 7 + ci) % 3) as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        let mode = inp["mode"].as_i64().unwrap();
        let (bcf, bcf2) = (
            inp["bcf"].as_i64().unwrap() as i32,
            inp["bcf2"].as_i64().unwrap() as i32,
        );
        let nt = inp["nt"].as_i64().unwrap() as i32;
        let mut got: Vec<String> = Vec::new();
        let ni = n as i32;
        let _ = (bcf, nt);
        let fm = if bcf2 & 1 != 0 { 1 } else { 2 };
        let r0 = set_forbidden_edges(&mut p, &at, ni, fm);
        got.push(format!("f {r0} {}", p.edge_forbidden_mask));
        let r1 = if mode != 0 {
            fix_special_bonds(&mut p, &at, ni, 4)
        } else {
            0
        };
        got.push(format!("g {r1} {}", p.edge_forbidden_mask));
        for e in p.edge.iter().take(p.num_edges as usize) {
            got.push(format!("e {}", e.forbidden));
        }
        for a in &at {
            let mut s = format!("T {} {} {} |", a.valence, a.chem_bonds_valence, a.radical);
            for j in 0..a.valence as usize {
                s += &format!(" {}", a.bond_type[j]);
            }
            got.push(s);
        }
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            got.push(format!(
                "V {} {} {} {} {} {}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, e.pass
            ));
        }
        for e in p.edge.iter().take(p.num_edges as usize) {
            got.push(format!(
                "E {} {} {} {} {}",
                e.cap, e.cap0, e.flow, e.flow0, e.pass
            ));
        }
        got.push(format!(
            "P {} {} {}",
            p.num_altp, p.b_change_flow, p.b_not_a_simple_path
        ));
        for a in 0..p.num_altp.min(p.max_altp) as usize {
            let ap = &p.altp[a];
            let len = ap[I_ALTP_PATH_LEN].number();
            let mut s = format!(
                "A {} {} {} {} {}",
                ap[I_ALTP_START_ATOM].number(),
                ap[I_ALTP_END_ATOM].number(),
                len,
                ap[I_ALTP_FLOW].flow[0],
                ap[I_ALTP_FLOW].flow[1]
            );
            for j in 0..len as usize {
                let e = &ap[I_ALTP_NEIGHBOR + j];
                s += &format!(" {},{}", e.ineigh(0), e.ineigh(1));
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn alternating_bond_blocks_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_bnblocks.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for (i, a) in inp["atoms"].as_array().unwrap().iter().enumerate() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                n_ring_system: ((i * 7 + ci) % 3) as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        let mode = inp["mode"].as_i64().unwrap();
        let (bcf, bcf2) = (
            inp["bcf"].as_i64().unwrap() as i32,
            inp["bcf2"].as_i64().unwrap() as i32,
        );
        let nt = inp["nt"].as_i64().unwrap() as i32;
        let mut got: Vec<String> = Vec::new();
        let ni = n as i32;
        let _ = (bcf, nt, mode);
        // the allocation changes the alternating bonds: restore the input
        for (i, a) in inp["atoms"].as_array().unwrap().iter().enumerate() {
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                at[i].bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at[i].endpoint = at[i].c_point;
        }
        let bu = bcf2 & 1 != 0;
        got.push(format!(
            "a {}",
            re_init_bn_struct_for_alt_bns(&mut p, &mut at, ni, bu)
        ));
        for e in p.edge.iter().take(p.num_edges as usize) {
            got.push(format!("e {} {} {} {}", e.pass, e.flow, e.cap, e.forbidden));
        }
        got.push(format!("b {}", mark_ring_systems_alt_bns(&mut p)));
        for e in p.edge.iter().take(p.num_edges as usize) {
            got.push(format!("e {} {} {}", e.pass, e.flow, e.cap));
        }
        got.push(format!(
            "c {}",
            mark_non_stereo_alt_bns(&p, &mut at, ni, bu)
        ));
        for a in &at {
            let mut s = String::from("s");
            for j in 0..a.valence as usize {
                s += &format!(" {}.{}", a.bond_type[j], a.bond_stereo[j]);
            }
            got.push(s);
        }
        for a in &at {
            let mut s = format!("T {} {} {} |", a.valence, a.chem_bonds_valence, a.radical);
            for j in 0..a.valence as usize {
                s += &format!(" {}", a.bond_type[j]);
            }
            got.push(s);
        }
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            got.push(format!(
                "V {} {} {} {} {} {}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, e.pass
            ));
        }
        for e in p.edge.iter().take(p.num_edges as usize) {
            got.push(format!(
                "E {} {} {} {} {}",
                e.cap, e.cap0, e.flow, e.flow0, e.pass
            ));
        }
        got.push(format!(
            "P {} {} {}",
            p.num_altp, p.b_change_flow, p.b_not_a_simple_path
        ));
        for a in 0..p.num_altp.min(p.max_altp) as usize {
            let ap = &p.altp[a];
            let len = ap[I_ALTP_PATH_LEN].number();
            let mut s = format!(
                "A {} {} {} {} {}",
                ap[I_ALTP_START_ATOM].number(),
                ap[I_ALTP_END_ATOM].number(),
                len,
                ap[I_ALTP_FLOW].flow[0],
                ap[I_ALTP_FLOW].flow[1]
            );
            for j in 0..len as usize {
                let e = &ap[I_ALTP_NEIGHBOR + j];
                s += &format!(" {},{}", e.ineigh(0), e.ineigh(1));
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn proton_restoration_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::TGroupInfo;
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump(at: &[InpAtom], t: &TGroupInfo, totals: &[i32], got: &mut Vec<String>) {
        for a in at {
            got.push(format!(
                "A {} {} {} {} {},{},{} t{} s{},{},{} o{},{},{} n{},{},{} g{},{},{} r{} nb{} e{} w{} p{} v{}",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.sb_parity[0],
                a.sb_parity[1],
                a.sb_parity[2],
                a.sb_ord[0],
                a.sb_ord[1],
                a.sb_ord[2],
                a.sn_ord[0],
                a.sn_ord[1],
                a.sn_ord[2],
                a.sn_orig_at_num[0],
                a.sn_orig_at_num[1],
                a.sn_orig_at_num[2],
                a.orig_at_number,
                a.neighbor[0],
                a.endpoint,
                a.iso_atw_diff,
                a.p_parity,
                a.valence
            ));
            let mut s = String::from("B");
            for j in 0..a.valence as usize {
                s += &format!(" {}.{}", a.neighbor[j], a.bond_type[j]);
            }
            got.push(s);
        }
        got.push(format!(
            "I {} {} {},{},{}",
            t.tni.n_num_removed_explicit_h,
            t.tni.n_num_removed_protons,
            t.tni.n_num_removed_protons_isotopic[0],
            t.tni.n_num_removed_protons_isotopic[1],
            t.tni.n_num_removed_protons_isotopic[2]
        ));
        let tot: Vec<String> = totals.iter().map(|x| format!(" {x}")).collect();
        got.push(format!("T{}", tot.concat()));
    }
    fn queries(at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        for i in 0..n {
            let mut sub = 0;
            let r = b_is_hard_rem_h_candidate(at, i, &mut sub);
            got.push(format!(
                "Q {r} {sub} {} {} {}",
                b_has_acidic_hydrogen(at, i),
                b_has_other_exchangable_h(at, i),
                b_has_acidic_minus(at, i)
            ));
        }
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_addrem.json")).unwrap();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let (sbp, sbo, sno, sng) = (
                ints(&a["sbp"]),
                ints(&a["sbo"]),
                ints(&a["sno"]),
                ints(&a["sng"]),
            );
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                p_parity: a["pp"].as_i64().unwrap() as i8,
                iso_atw_diff: a["iw"].as_i64().unwrap() as i8,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
                x.sb_parity[k] = sbp[k] as i8;
                x.sb_ord[k] = sbo[k] as i8;
                x.sn_ord[k] = sno[k] as i8;
                x.sn_orig_at_num[k] = sng[k] as u16;
            }
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                iso_atw_diff: h[2].as_i64().unwrap() as i8,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let mut tgi = TGroupInfo::default();
        tgi.tni.n_num_removed_explicit_h = hs.len() as i32;
        tgi.tni.n_num_removed_protons = 1;
        let mut totals = vec![0i32; attot::ARRAY_LEN];
        for i in 0..n {
            get_atom_charge_type(&at, i, Some(&mut totals), 0);
        }
        let mut got: Vec<String> = Vec::new();
        queries(&at, n, &mut got);
        dump(&at, &tgi, &totals, &mut got);
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let ni = n as i32;
            if o[0] == 0 {
                let mut np = o[1];
                let flags = if o[3] != 0 {
                    FLAG_PROTON_AC_HARD_ADDED
                } else {
                    0
                };
                let r = add_remove_protons_restr(&mut at, ni, &mut np, o[2], flags, o[4], 0, 0);
                got.push(format!("R 0 {r} {np}"));
            } else {
                let mut np = [o[1] as i16, o[2] as i16, o[3] as i16];
                let r = add_remove_iso_protons_restr(&mut at, ni, &mut np, o[4]);
                got.push(format!("R 1 {r} {} {} {}", np[0], np[1], np[2]));
            }
            queries(&at, n, &mut got);
            dump(&at, &tgi, &totals, &mut got);
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
}

#[test]
fn charge_group_marking_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::{CCandidate, CGroup, CGroupInfo, mark_charge_groups};
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_state(p: &BnStruct, at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        got.push(format!(
            "C {} {} {} {}",
            p.num_vertices, p.num_edges, p.tot_st_cap, p.tot_st_flow
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges
            );
        }
        got.push(s);
        let mut s = String::from("F");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(" {}:{}:{}:{}", e.cap, e.cap0, e.flow, e.flow0);
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at.iter().take(n) {
            s += &format!(
                " {},{},{},{},{},{},{},{},{}:",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.endpoint
            );
            for j in 0..a.valence as usize {
                s += &format!("{}.", a.bond_type[j]);
            }
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_markcg.json")).unwrap();
    let mut seen = std::collections::HashMap::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                at_type: a["att"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
            }
            x.set_elname(match x.el_number {
                1 => "H",
                6 => "C",
                7 => "N",
                8 => "O",
                15 => "P",
                16 => "S",
                17 => "Cl",
                35 => "Br",
                _ => "Xx",
            });
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        for (i, a) in at.iter_mut().enumerate() {
            a.n_ring_system = (i % 2) as u16;
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        let tf_packed = inp["tf"].as_u64().unwrap();
        p.pb_taut_flags = (tf_packed & 0xffff) as u32;
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            let nv = p.num_vertices;
            if op == 10 {
                if a >= nv {
                    continue;
                }
                add_new_vertex(&mut p, a, x, y, z, &mut dots);
            } else if op == 11 {
                if a < p.num_atoms || a >= nv {
                    continue;
                }
                p.vert[a as usize].type_ = x as u16;
            } else {
                if a >= nv || b >= nv || a == b {
                    continue;
                }
                add_new_edge(&mut p, a, b, x, y);
            }
        }
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let len = 2 * p.max_vertices as usize + 8;
        bd.rad_endpoints = vec![0; len];
        bd.rad_edges = vec![0; len];
        let mut got: Vec<String> = Vec::new();
        dump_state(&p, &at, n, &mut got);
        let ni = n as i32;
        bns_adjust_flow_bonds_rad(&mut p, &mut bd, &mut at, ni);
        dump_state(&p, &at, n, &mut got);
        let mut cgi = CGroupInfo {
            max_num_c_groups: ((tf_packed >> 16) & 0xff) as usize,
            max_num_candidates: ((tf_packed >> 24) & 0xff) as i32,
            ..Default::default()
        };
        cgi.c_group = Some(vec![CGroup::default(); cgi.max_num_c_groups + 1]);
        cgi.c_candidate = vec![CCandidate::default(); cgi.max_num_candidates as usize + 1];
        for _ in inp["calls"].as_array().unwrap() {
            let ret = mark_charge_groups(&mut at, ni, &mut cgi, &mut p, &mut bd);
            *seen.entry(ret).or_insert(0) += 1;
            got.push(format!(
                "R {ret} {} {}",
                cgi.num_c_groups, cgi.num_candidates
            ));
            for g in &cgi.c_group.as_ref().unwrap()[..cgi.num_c_groups] {
                got.push(format!(
                    "G {} {} {} {} {}",
                    g.num[0], g.num[1], g.num_c_points, g.n_group_number, g.c_group_type
                ));
            }
            for c in &cgi.c_candidate[..cgi.max_num_candidates as usize] {
                got.push(format!("K {} {} {}", c.atnumber, c.type_, c.subtype));
            }
            let mut s = String::from("P");
            for a in at.iter().take(n) {
                s += &format!(" {}", a.c_point);
            }
            got.push(s);
            dump_state(&p, &at, n, &mut got);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(seen.keys().any(|&k| k > 0) && seen.len() > 3, "{seen:?}");
}

#[test]
fn end_point_registration_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::{
        CGroup, CGroupInfo, TBondPos, TEndpoint, TGroup, TGroupInfo, register_end_points,
        set_tautomeric_bonds,
    };
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_state(
        p: &BnStruct,
        at: &[InpAtom],
        t: &TGroupInfo,
        cg: &CGroupInfo,
        got: &mut Vec<String>,
    ) {
        got.push(format!(
            "C {} {} {} {} {} {}",
            p.num_vertices,
            p.num_edges,
            p.num_t_groups,
            p.num_c_groups,
            p.num_added_atoms,
            p.num_added_edges
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}:{}:{}|",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges, v.max_adj_edges, v.iedge
            );
            for k in 0..v.num_adj_edges as usize {
                s += &format!("{}.", p.iedge[v.iedge + k]);
            }
        }
        got.push(s);
        let mut s = String::from("E");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(
                " {}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
                e.cap,
                e.cap0,
                e.flow,
                e.flow0,
                e.pass,
                e.forbidden,
                e.neighbor1,
                e.neighbor12,
                e.neigh_ord[0],
                e.neigh_ord[1]
            );
        }
        got.push(s);
        let mut s = String::from("G");
        for g in t.t_group.iter().flatten().take(t.num_t_groups) {
            s += &format!(" {}", g.n_group_number);
        }
        s += " |";
        for g in cg.c_group.iter().flatten().take(cg.num_c_groups) {
            s += &format!(" {}", g.n_group_number);
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at {
            s += &format!(" {},{}", a.endpoint, a.c_point);
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_regendp.json")).unwrap();
    let mut codes = std::collections::HashSet::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                c_point: a["cp"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        p.pb_taut_flags = inp["tfl"].as_u64().unwrap() as u32;
        let tg_in = inp["tg"].as_array().unwrap();
        let max_t = inp["maxt"].as_u64().unwrap() as usize;
        let mut t_group = vec![TGroup::default(); max_t + 2];
        for (i, g) in tg_in.iter().enumerate() {
            let g = g.as_array().unwrap();
            let nums = ints(&g[2]);
            let das = ints(&g[3]);
            t_group[i].n_group_number = g[0].as_u64().unwrap() as u16;
            t_group[i].n_num_endpoints = g[1].as_u64().unwrap() as u16;
            for k in 0..5 {
                t_group[i].num[k] = nums[k] as u16;
            }
            for k in 0..6 {
                t_group[i].num_da[k] = das[k] as u16;
            }
        }
        let mut tgi = TGroupInfo {
            t_group: Some(t_group),
            num_t_groups: tg_in.len(),
            max_num_t_groups: max_t,
            b_taut_flags: inp["tf"].as_u64().unwrap(),
            ..Default::default()
        };
        let mut cgi = CGroupInfo {
            c_group: Some(
                inp["cg"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|g| {
                        let g = ints(g);
                        CGroup {
                            n_group_number: g[0] as u16,
                            num_c_points: g[1] as u16,
                            ..Default::default()
                        }
                    })
                    .collect(),
            ),
            num_c_groups: inp["cg"].as_array().unwrap().len(),
            ..Default::default()
        };
        let use_bns = inp["usebns"].as_i64().unwrap() != 0;
        let mut eps: Vec<TEndpoint> = inp["ep"]
            .as_array()
            .unwrap()
            .iter()
            .map(|e| {
                let e = e.as_array().unwrap();
                let nums = ints(&e[3]);
                let das = ints(&e[4]);
                let mut x = TEndpoint {
                    n_atom_number: e[0].as_u64().unwrap() as u16,
                    n_group_number: e[1].as_u64().unwrap() as u16,
                    n_equ_number: e[2].as_u64().unwrap() as u16,
                    ..Default::default()
                };
                for k in 0..5 {
                    x.num[k] = nums[k] as u16;
                }
                for k in 0..6 {
                    x.num_da[k] = das[k] as u16;
                }
                x
            })
            .collect();
        let mut got: Vec<String> = Vec::new();
        for _ in 0..2 {
            let has_c = cgi.num_c_groups > 0;
            let ret = register_end_points(
                &mut tgi,
                &mut eps,
                &mut at,
                n as i32,
                if has_c { Some(&mut cgi) } else { None },
                if use_bns { Some(&mut p) } else { None },
            );
            codes.insert(ret);
            got.push(format!("R {ret} {}", tgi.num_t_groups));
            for g in &tgi.t_group.as_ref().unwrap()[..tgi.num_t_groups] {
                let nums: Vec<String> = g.num.iter().map(|x| format!(" {x}")).collect();
                let das: Vec<String> = g.num_da.iter().map(|x| format!(" {x}")).collect();
                got.push(format!(
                    "T {} {} |{} |{}",
                    g.n_group_number,
                    g.n_num_endpoints,
                    nums.concat(),
                    das.concat()
                ));
            }
            let mut s = String::from("I");
            for v in &tgi.t_group_number {
                s += &format!(" {v}");
            }
            got.push(s);
            let mut s = String::from("Q");
            for e in &eps {
                s += &format!(
                    " {},{},{}",
                    e.n_atom_number, e.n_group_number, e.n_equ_number
                );
            }
            got.push(s);
            dump_state(&p, &at, &tgi, &cgi, &mut got);
            for e in eps.iter_mut() {
                e.n_group_number = at[e.n_atom_number as usize].endpoint;
            }
        }
        let bp: Vec<TBondPos> = inp["bp"]
            .as_array()
            .unwrap()
            .iter()
            .map(|b| TBondPos {
                n_atom_number: b[0].as_u64().unwrap() as u16,
                neighbor_index: b[1].as_u64().unwrap() as u16,
            })
            .collect();
        got.push(format!("S {}", set_tautomeric_bonds(&mut at, &bp)));
        for a in at.iter().take(n) {
            let mut s = String::from("B");
            for j in 0..a.valence as usize {
                s += &format!(" {}", a.bond_type[j]);
            }
            got.push(s);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(codes.len() >= 3, "{codes:?}");
}

#[test]
fn accessible_end_points_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::{
        CGroup, CGroupInfo, TBondPos, TEndpoint, find_accessible_end_points,
        get_neutral_reps_if_needed,
    };
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_state(p: &BnStruct, at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        got.push(format!(
            "C {} {} {} {}",
            p.num_vertices, p.num_edges, p.tot_st_cap, p.tot_st_flow
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges
            );
        }
        got.push(s);
        let mut s = String::from("F");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(" {}:{}:{}:{}", e.cap, e.cap0, e.flow, e.flow0);
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at.iter().take(n) {
            s += &format!(
                " {},{},{},{},{},{},{},{},{}:",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.endpoint
            );
            for j in 0..a.valence as usize {
                s += &format!("{}.", a.bond_type[j]);
            }
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_findacc.json")).unwrap();
    let mut seen = std::collections::HashMap::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                at_type: a["att"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
            }
            x.set_elname(match x.el_number {
                1 => "H",
                6 => "C",
                7 => "N",
                8 => "O",
                15 => "P",
                16 => "S",
                17 => "Cl",
                35 => "Br",
                _ => "Xx",
            });
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        for (i, a) in at.iter_mut().enumerate() {
            a.n_ring_system = (i % 2) as u16;
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        let tf_packed = inp["tf"].as_u64().unwrap();
        p.pb_taut_flags = (tf_packed & 0xffff) as u32;
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            let nv = p.num_vertices;
            if op == 10 {
                if a >= nv {
                    continue;
                }
                add_new_vertex(&mut p, a, x, y, z, &mut dots);
            } else if op == 11 {
                if a < p.num_atoms || a >= nv {
                    continue;
                }
                p.vert[a as usize].type_ = x as u16;
            } else {
                if a >= nv || b >= nv || a == b {
                    continue;
                }
                add_new_edge(&mut p, a, b, x, y);
            }
        }
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let len = 2 * p.max_vertices as usize + 8;
        bd.rad_endpoints = vec![0; len];
        bd.rad_edges = vec![0; len];
        let mut got: Vec<String> = Vec::new();
        dump_state(&p, &at, n, &mut got);
        let ni = n as i32;
        bns_adjust_flow_bonds_rad(&mut p, &mut bd, &mut at, ni);
        dump_state(&p, &at, n, &mut got);
        let max_cp = at.iter().take(n).map(|a| a.c_point).max().unwrap_or(0) as usize;
        let mut cgi = CGroupInfo {
            num_c_groups: max_cp,
            max_num_c_groups: max_cp,
            ..Default::default()
        };
        let mut groups = vec![CGroup::default(); max_cp + 1];
        for g in 1..=max_cp {
            groups[g - 1].n_group_number = g as u16;
            for a in at.iter().take(n) {
                if a.c_point as usize == g {
                    groups[g - 1].num_c_points += 1;
                    groups[g - 1].num[0] += u16::from(a.charge == 1);
                }
            }
        }
        cgi.c_group = Some(groups);
        let calls = inp["calls"].as_array().unwrap();
        let mut eps: Vec<TEndpoint> = Vec::new();
        let mut bps: Vec<TBondPos> = Vec::new();
        for cl in calls {
            let cl = ints(cl);
            eps.push(TEndpoint {
                n_atom_number: cl[0] as u16,
                n_group_number: cl[1] as u16,
                ..Default::default()
            });
            bps.push(TBondPos {
                n_atom_number: cl[0] as u16,
                neighbor_index: cl[2] as u16,
            });
        }
        let (mut n_end, mut n_bond) = (eps.len(), bps.len());
        let cgi_ref = if max_cp > 0 { Some(&cgi) } else { None };
        let ret = find_accessible_end_points(
            &mut eps,
            &mut n_end,
            &mut bps,
            &mut n_bond,
            &mut p,
            &mut bd,
            &mut at,
            ni,
            cgi_ref,
            ((tf_packed >> 16) & 0xff) as i32,
        );
        *seen.entry(ret).or_insert(0) += 1;
        got.push(format!("R {ret} {n_end} {n_bond}"));
        for i in 0..n_end {
            got.push(format!(
                "E {} {} {} {}",
                eps[i].n_atom_number,
                eps[i].n_group_number,
                eps[i].n_equ_number,
                bps[i].neighbor_index
            ));
        }
        {
            let mut ri = if eps.is_empty() {
                0
            } else {
                eps[0].n_atom_number
            };
            let mut rj = if eps.len() > 1 {
                eps[1].n_atom_number
            } else {
                0
            };
            get_neutral_reps_if_needed(&mut ri, &mut rj, &at, ni, &eps[..n_end], cgi_ref);
            got.push(format!("G {ri} {rj}"));
        }
        dump_state(&p, &at, n, &mut got);
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(seen.keys().any(|&k| k > 0) && seen.len() > 3, "{seen:?}");
}

#[test]
fn salt_charge_groups2_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::{
        CGroup, CGroupInfo, SCandidate, SGroupInfo, TGroup, TGroupInfo, mark_salt_charge_groups2,
    };
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_state(p: &BnStruct, at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        got.push(format!(
            "C {} {} {} {}",
            p.num_vertices, p.num_edges, p.tot_st_cap, p.tot_st_flow
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges
            );
        }
        got.push(s);
        let mut s = String::from("F");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(" {}:{}:{}:{}", e.cap, e.cap0, e.flow, e.flow0);
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at.iter().take(n) {
            s += &format!(
                " {},{},{},{},{},{},{},{},{}:",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.endpoint
            );
            for j in 0..a.valence as usize {
                s += &format!("{}.", a.bond_type[j]);
            }
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_saltcg2.json")).unwrap();
    let mut seen = std::collections::HashMap::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                at_type: a["att"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
            }
            x.set_elname(match x.el_number {
                1 => "H",
                6 => "C",
                7 => "N",
                8 => "O",
                15 => "P",
                16 => "S",
                17 => "Cl",
                35 => "Br",
                _ => "Xx",
            });
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        for (i, a) in at.iter_mut().enumerate() {
            a.n_ring_system = (i % 2) as u16;
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        let tf_packed = inp["tf"].as_u64().unwrap();
        p.pb_taut_flags = (tf_packed & 0xffff) as u32;
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            let nv = p.num_vertices;
            if op == 10 {
                if a >= nv {
                    continue;
                }
                add_new_vertex(&mut p, a, x, y, z, &mut dots);
            } else if op == 11 {
                if a < p.num_atoms || a >= nv {
                    continue;
                }
                p.vert[a as usize].type_ = x as u16;
            } else {
                if a >= nv || b >= nv || a == b {
                    continue;
                }
                add_new_edge(&mut p, a, b, x, y);
            }
        }
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let len = 2 * p.max_vertices as usize + 8;
        bd.rad_endpoints = vec![0; len];
        bd.rad_edges = vec![0; len];
        let mut got: Vec<String> = Vec::new();
        dump_state(&p, &at, n, &mut got);
        let ni = n as i32;
        bns_adjust_flow_bonds_rad(&mut p, &mut bd, &mut at, ni);
        dump_state(&p, &at, n, &mut got);
        let max_cp = at.iter().take(n).map(|a| a.c_point).max().unwrap_or(0) as usize;
        let max_ep = at.iter().take(n).map(|a| a.endpoint).max().unwrap_or(0) as usize;
        let mut cgi = CGroupInfo {
            num_c_groups: max_cp,
            max_num_c_groups: max_cp,
            ..Default::default()
        };
        let mut groups = vec![CGroup::default(); max_cp + 1];
        for g in 1..=max_cp {
            groups[g - 1].n_group_number = g as u16;
            for a in at.iter().take(n) {
                if a.c_point as usize == g {
                    groups[g - 1].num_c_points += 1;
                    groups[g - 1].num[0] += u16::from(a.charge == 1);
                }
            }
        }
        cgi.c_group = Some(groups);
        let fch = (tf_packed >> 16) & 0xff;
        let mut t_group = vec![TGroup::default(); max_ep + 5];
        for g in 1..=max_ep {
            t_group[g - 1].n_group_number = g as u16;
            for a in at.iter().take(n) {
                if a.endpoint as usize == g {
                    t_group[g - 1].n_num_endpoints += 1;
                    t_group[g - 1].num[0] += u16::from(a.num_h > 0);
                    t_group[g - 1].num[1] += u16::from(a.charge == -1);
                }
            }
        }
        let mut tgi = TGroupInfo {
            t_group: Some(t_group),
            num_t_groups: max_ep,
            max_num_t_groups: max_ep + 3,
            b_taut_flags: if fch & 1 != 0 { 0x20 } else { 0 },
            b_taut_flags_done: if fch & 2 != 0 { 0x2000 } else { 0 },
            ..Default::default()
        };
        tgi.tni.b_normalization_flags = [0u64, 0x2, 0x10, 0x20][((fch >> 2) & 3) as usize];
        let max_cand = ((tf_packed >> 24) & 0xff) as i32;
        let mut sgi = SGroupInfo {
            s_candidate: vec![SCandidate::default(); max_cand as usize + 2],
            max_num_candidates: max_cand,
            ..Default::default()
        };
        let cgi_has = max_cp > 0;
        for _ in inp["calls"].as_array().unwrap() {
            let ret = mark_salt_charge_groups2(
                &mut at,
                ni,
                &mut sgi,
                &mut tgi,
                if cgi_has { Some(&mut cgi) } else { None },
                &mut p,
                &mut bd,
            );
            *seen.entry(ret).or_insert(0) += 1;
            got.push(format!(
                "R {ret} {} {}",
                sgi.num_candidates, tgi.num_t_groups
            ));
            for c in &sgi.s_candidate[..max_cand as usize + 1] {
                got.push(format!(
                    "K {} {} {} {}",
                    c.atnumber, c.type_, c.subtype, c.endpoint
                ));
            }
            for g in &tgi.t_group.as_ref().unwrap()[..tgi.num_t_groups] {
                let nums: Vec<String> = g.num.iter().map(|x| format!(" {x}")).collect();
                let das: Vec<String> = g.num_da.iter().map(|x| format!(" {x}")).collect();
                got.push(format!(
                    "T {} {} |{} |{}",
                    g.n_group_number,
                    g.n_num_endpoints,
                    nums.concat(),
                    das.concat()
                ));
            }
            got.push(format!("F {}", tgi.b_taut_flags_done));
            let mut s = String::from("P");
            for a in at.iter().take(n) {
                s += &format!(" {},{}", a.endpoint, a.c_point);
            }
            got.push(s);
            dump_state(&p, &at, n, &mut got);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(seen.keys().any(|&k| k > 0) && seen.len() > 3, "{seen:?}");
}

#[test]
fn salt_charge_groups_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::{
        CGroup, CGroupInfo, SCandidate, SGroupInfo, TGroup, TGroupInfo, mark_salt_charge_groups,
    };
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_state(p: &BnStruct, at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        got.push(format!(
            "C {} {} {} {}",
            p.num_vertices, p.num_edges, p.tot_st_cap, p.tot_st_flow
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges
            );
        }
        got.push(s);
        let mut s = String::from("F");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(" {}:{}:{}:{}", e.cap, e.cap0, e.flow, e.flow0);
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at.iter().take(n) {
            s += &format!(
                " {},{},{},{},{},{},{},{},{}:",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.endpoint
            );
            for j in 0..a.valence as usize {
                s += &format!("{}.", a.bond_type[j]);
            }
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_saltcg.json")).unwrap();
    let mut seen = std::collections::HashMap::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                at_type: a["att"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
            }
            x.set_elname(match x.el_number {
                1 => "H",
                6 => "C",
                7 => "N",
                8 => "O",
                15 => "P",
                16 => "S",
                17 => "Cl",
                35 => "Br",
                _ => "Xx",
            });
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        for (i, a) in at.iter_mut().enumerate() {
            a.n_ring_system = (i % 2) as u16;
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        let tf_packed = inp["tf"].as_u64().unwrap();
        p.pb_taut_flags = (tf_packed & 0xffff) as u32;
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            let nv = p.num_vertices;
            if op == 10 {
                if a >= nv {
                    continue;
                }
                add_new_vertex(&mut p, a, x, y, z, &mut dots);
            } else if op == 11 {
                if a < p.num_atoms || a >= nv {
                    continue;
                }
                p.vert[a as usize].type_ = x as u16;
            } else {
                if a >= nv || b >= nv || a == b {
                    continue;
                }
                add_new_edge(&mut p, a, b, x, y);
            }
        }
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let len = 2 * p.max_vertices as usize + 8;
        bd.rad_endpoints = vec![0; len];
        bd.rad_edges = vec![0; len];
        let mut got: Vec<String> = Vec::new();
        dump_state(&p, &at, n, &mut got);
        let ni = n as i32;
        bns_adjust_flow_bonds_rad(&mut p, &mut bd, &mut at, ni);
        dump_state(&p, &at, n, &mut got);
        let max_cp = at.iter().take(n).map(|a| a.c_point).max().unwrap_or(0) as usize;
        let max_ep = at.iter().take(n).map(|a| a.endpoint).max().unwrap_or(0) as usize;
        let mut cgi = CGroupInfo {
            num_c_groups: max_cp,
            max_num_c_groups: max_cp,
            ..Default::default()
        };
        let mut groups = vec![CGroup::default(); max_cp + 1];
        for g in 1..=max_cp {
            groups[g - 1].n_group_number = g as u16;
            for a in at.iter().take(n) {
                if a.c_point as usize == g {
                    groups[g - 1].num_c_points += 1;
                    groups[g - 1].num[0] += u16::from(a.charge == 1);
                }
            }
        }
        cgi.c_group = Some(groups);
        let fch = (tf_packed >> 16) & 0xff;
        let mut t_group = vec![TGroup::default(); max_ep + 5];
        for g in 1..=max_ep {
            t_group[g - 1].n_group_number = g as u16;
            for a in at.iter().take(n) {
                if a.endpoint as usize == g {
                    t_group[g - 1].n_num_endpoints += 1;
                    t_group[g - 1].num[0] += u16::from(a.num_h > 0);
                    t_group[g - 1].num[1] += u16::from(a.charge == -1);
                }
            }
        }
        let mut tgi = TGroupInfo {
            t_group: Some(t_group),
            num_t_groups: max_ep,
            max_num_t_groups: max_ep + 3,
            b_taut_flags: if fch & 1 != 0 { 0x20 } else { 0 },
            b_taut_flags_done: if fch & 2 != 0 { 0x2000 } else { 0 },
            ..Default::default()
        };
        tgi.tni.b_normalization_flags = [0u64, 0x2, 0x10, 0x20][((fch >> 2) & 3) as usize];
        let max_cand = ((tf_packed >> 24) & 0xff) as i32;
        let mut sgi = SGroupInfo {
            s_candidate: vec![SCandidate::default(); max_cand as usize + 2],
            max_num_candidates: max_cand,
            ..Default::default()
        };
        let cgi_has = max_cp > 0;
        for _ in inp["calls"].as_array().unwrap() {
            let ret = mark_salt_charge_groups(
                &mut at,
                ni,
                &mut sgi,
                &mut tgi,
                if cgi_has { Some(&mut cgi) } else { None },
                &mut p,
                &mut bd,
            );
            *seen.entry(ret).or_insert(0) += 1;
            got.push(format!(
                "R {ret} {} {}",
                sgi.num_candidates, tgi.num_t_groups
            ));
            for c in &sgi.s_candidate[..max_cand as usize + 1] {
                got.push(format!(
                    "K {} {} {} {}",
                    c.atnumber, c.type_, c.subtype, c.endpoint
                ));
            }
            for g in &tgi.t_group.as_ref().unwrap()[..tgi.num_t_groups] {
                let nums: Vec<String> = g.num.iter().map(|x| format!(" {x}")).collect();
                let das: Vec<String> = g.num_da.iter().map(|x| format!(" {x}")).collect();
                got.push(format!(
                    "T {} {} |{} |{}",
                    g.n_group_number,
                    g.n_num_endpoints,
                    nums.concat(),
                    das.concat()
                ));
            }
            got.push(format!("F {}", tgi.b_taut_flags_done));
            let mut s = String::from("P");
            for a in at.iter().take(n) {
                s += &format!(" {},{}", a.endpoint, a.c_point);
            }
            got.push(s);
            dump_state(&p, &at, n, &mut got);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(seen.keys().any(|&k| k > 0) && seen.len() > 3, "{seen:?}");
}

#[test]
fn salt_merge_groups_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::{
        CGroup, CGroupInfo, SCandidate, SGroupInfo, TGroup, TGroupInfo, merge_salt_taut_groups,
    };
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_state(p: &BnStruct, at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        got.push(format!(
            "C {} {} {} {}",
            p.num_vertices, p.num_edges, p.tot_st_cap, p.tot_st_flow
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges
            );
        }
        got.push(s);
        let mut s = String::from("F");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(" {}:{}:{}:{}", e.cap, e.cap0, e.flow, e.flow0);
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at.iter().take(n) {
            s += &format!(
                " {},{},{},{},{},{},{},{},{}:",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.endpoint
            );
            for j in 0..a.valence as usize {
                s += &format!("{}.", a.bond_type[j]);
            }
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_saltmerge.json")).unwrap();
    let mut seen = std::collections::HashMap::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                at_type: a["att"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
            }
            x.set_elname(match x.el_number {
                1 => "H",
                6 => "C",
                7 => "N",
                8 => "O",
                15 => "P",
                16 => "S",
                17 => "Cl",
                35 => "Br",
                _ => "Xx",
            });
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        for (i, a) in at.iter_mut().enumerate() {
            a.n_ring_system = (i % 2) as u16;
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        let tf_packed = inp["tf"].as_u64().unwrap();
        p.pb_taut_flags = (tf_packed & 0xffff) as u32;
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            let nv = p.num_vertices;
            if op == 10 {
                if a >= nv {
                    continue;
                }
                add_new_vertex(&mut p, a, x, y, z, &mut dots);
            } else if op == 11 {
                if a < p.num_atoms || a >= nv {
                    continue;
                }
                p.vert[a as usize].type_ = x as u16;
            } else {
                if a >= nv || b >= nv || a == b {
                    continue;
                }
                add_new_edge(&mut p, a, b, x, y);
            }
        }
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let len = 2 * p.max_vertices as usize + 8;
        bd.rad_endpoints = vec![0; len];
        bd.rad_edges = vec![0; len];
        let mut got: Vec<String> = Vec::new();
        dump_state(&p, &at, n, &mut got);
        let ni = n as i32;
        bns_adjust_flow_bonds_rad(&mut p, &mut bd, &mut at, ni);
        dump_state(&p, &at, n, &mut got);
        let max_cp = at.iter().take(n).map(|a| a.c_point).max().unwrap_or(0) as usize;
        let max_ep = at.iter().take(n).map(|a| a.endpoint).max().unwrap_or(0) as usize;
        let mut cgi = CGroupInfo {
            num_c_groups: max_cp,
            max_num_c_groups: max_cp,
            ..Default::default()
        };
        let mut groups = vec![CGroup::default(); max_cp + 1];
        for g in 1..=max_cp {
            groups[g - 1].n_group_number = g as u16;
            for a in at.iter().take(n) {
                if a.c_point as usize == g {
                    groups[g - 1].num_c_points += 1;
                    groups[g - 1].num[0] += u16::from(a.charge == 1);
                }
            }
        }
        cgi.c_group = Some(groups);
        let fch = (tf_packed >> 16) & 0xff;
        let mut t_group = vec![TGroup::default(); max_ep + 5];
        for g in 1..=max_ep {
            t_group[g - 1].n_group_number = g as u16;
            for a in at.iter().take(n) {
                if a.endpoint as usize == g {
                    t_group[g - 1].n_num_endpoints += 1;
                    t_group[g - 1].num[0] += u16::from(a.num_h > 0);
                    t_group[g - 1].num[1] += u16::from(a.charge == -1);
                }
            }
        }
        let mut tgi = TGroupInfo {
            t_group: Some(t_group),
            num_t_groups: max_ep,
            max_num_t_groups: max_ep + 3,
            b_taut_flags: if fch & 1 != 0 { 0x20 } else { 0 },
            b_taut_flags_done: if fch & 2 != 0 { 0x2000 } else { 0 },
            ..Default::default()
        };
        tgi.tni.b_normalization_flags = [0u64, 0x2, 0x10, 0x20][((fch >> 2) & 3) as usize];
        let max_cand = ((tf_packed >> 24) & 0xff) as i32;
        let mut sgi = SGroupInfo {
            s_candidate: vec![SCandidate::default(); max_cand as usize + 2],
            max_num_candidates: max_cand,
            ..Default::default()
        };
        let cgi_has = max_cp > 0;
        for _ in inp["calls"].as_array().unwrap() {
            let ret = merge_salt_taut_groups(
                &mut at,
                ni,
                &mut sgi,
                &mut tgi,
                if cgi_has { Some(&mut cgi) } else { None },
                &mut p,
            );
            *seen.entry(ret).or_insert(0) += 1;
            got.push(format!(
                "R {ret} {} {}",
                sgi.num_candidates, tgi.num_t_groups
            ));
            for c in &sgi.s_candidate[..max_cand as usize + 1] {
                got.push(format!(
                    "K {} {} {} {}",
                    c.atnumber, c.type_, c.subtype, c.endpoint
                ));
            }
            for g in &tgi.t_group.as_ref().unwrap()[..tgi.num_t_groups] {
                let nums: Vec<String> = g.num.iter().map(|x| format!(" {x}")).collect();
                let das: Vec<String> = g.num_da.iter().map(|x| format!(" {x}")).collect();
                got.push(format!(
                    "T {} {} |{} |{}",
                    g.n_group_number,
                    g.n_num_endpoints,
                    nums.concat(),
                    das.concat()
                ));
            }
            got.push(format!("F {}", tgi.b_taut_flags_done));
            let mut s = String::from("P");
            for a in at.iter().take(n) {
                s += &format!(" {},{}", a.endpoint, a.c_point);
            }
            got.push(s);
            dump_state(&p, &at, n, &mut got);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(seen.keys().any(|&k| k > 0) && seen.len() > 3, "{seen:?}");
}

#[test]
fn isotopic_h_group_matches_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::{
        CGroup, CGroupInfo, SCandidate, SGroupInfo, TGroup, TGroupInfo, make_isotopic_h_group,
    };
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_state(p: &BnStruct, at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        got.push(format!(
            "C {} {} {} {}",
            p.num_vertices, p.num_edges, p.tot_st_cap, p.tot_st_flow
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges
            );
        }
        got.push(s);
        let mut s = String::from("F");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(" {}:{}:{}:{}", e.cap, e.cap0, e.flow, e.flow0);
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at.iter().take(n) {
            s += &format!(
                " {},{},{},{},{},{},{},{},{}:",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.endpoint
            );
            for j in 0..a.valence as usize {
                s += &format!("{}.", a.bond_type[j]);
            }
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_isoh.json")).unwrap();
    let mut seen = std::collections::HashMap::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                at_type: a["att"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
            }
            x.set_elname(match x.el_number {
                1 => "H",
                6 => "C",
                7 => "N",
                8 => "O",
                15 => "P",
                16 => "S",
                17 => "Cl",
                35 => "Br",
                _ => "Xx",
            });
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        for (i, a) in at.iter_mut().enumerate() {
            a.n_ring_system = (i % 2) as u16;
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        let tf_packed = inp["tf"].as_u64().unwrap();
        p.pb_taut_flags = (tf_packed & 0xffff) as u32;
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            let nv = p.num_vertices;
            if op == 10 {
                if a >= nv {
                    continue;
                }
                add_new_vertex(&mut p, a, x, y, z, &mut dots);
            } else if op == 11 {
                if a < p.num_atoms || a >= nv {
                    continue;
                }
                p.vert[a as usize].type_ = x as u16;
            } else {
                if a >= nv || b >= nv || a == b {
                    continue;
                }
                add_new_edge(&mut p, a, b, x, y);
            }
        }
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let len = 2 * p.max_vertices as usize + 8;
        bd.rad_endpoints = vec![0; len];
        bd.rad_edges = vec![0; len];
        let mut got: Vec<String> = Vec::new();
        dump_state(&p, &at, n, &mut got);
        let ni = n as i32;
        bns_adjust_flow_bonds_rad(&mut p, &mut bd, &mut at, ni);
        dump_state(&p, &at, n, &mut got);
        let max_cp = at.iter().take(n).map(|a| a.c_point).max().unwrap_or(0) as usize;
        let max_ep = at.iter().take(n).map(|a| a.endpoint).max().unwrap_or(0) as usize;
        let mut cgi = CGroupInfo {
            num_c_groups: max_cp,
            max_num_c_groups: max_cp,
            ..Default::default()
        };
        let mut groups = vec![CGroup::default(); max_cp + 1];
        for g in 1..=max_cp {
            groups[g - 1].n_group_number = g as u16;
            for a in at.iter().take(n) {
                if a.c_point as usize == g {
                    groups[g - 1].num_c_points += 1;
                    groups[g - 1].num[0] += u16::from(a.charge == 1);
                }
            }
        }
        cgi.c_group = Some(groups);
        let fch = (tf_packed >> 16) & 0xff;
        let mut t_group = vec![TGroup::default(); max_ep + 5];
        for g in 1..=max_ep {
            t_group[g - 1].n_group_number = g as u16;
            for a in at.iter().take(n) {
                if a.endpoint as usize == g {
                    t_group[g - 1].n_num_endpoints += 1;
                    t_group[g - 1].num[0] += u16::from(a.num_h > 0);
                    t_group[g - 1].num[1] += u16::from(a.charge == -1);
                }
            }
        }
        let mut tgi = TGroupInfo {
            t_group_number: (0..max_ep + 5)
                .map(|g| {
                    if g >= 1 && g <= max_ep && !(fch >= 60 && g == 1) {
                        g as u16
                    } else {
                        0
                    }
                })
                .collect(),
            t_group: Some(t_group),
            num_t_groups: max_ep,
            max_num_t_groups: max_ep + 3,
            b_taut_flags: if fch & 1 != 0 { 0x20 } else { 0 },
            b_taut_flags_done: if fch & 2 != 0 { 0x2000 } else { 0 },
            ..Default::default()
        };
        tgi.tni.b_normalization_flags = [0u64, 0x2, 0x10, 0x20][((fch >> 2) & 3) as usize];
        let max_cand = ((tf_packed >> 24) & 0xff) as i32;
        let mut sgi = SGroupInfo {
            s_candidate: vec![SCandidate::default(); max_cand as usize + 2],
            max_num_candidates: max_cand,
            ..Default::default()
        };
        for _ in inp["calls"].as_array().unwrap() {
            let ret = make_isotopic_h_group(&mut at, ni, &mut sgi, &mut tgi);
            *seen.entry(ret).or_insert(0) += 1;
            got.push(format!(
                "R {ret} {} {}",
                sgi.num_candidates, tgi.num_t_groups
            ));
            for c in &sgi.s_candidate[..max_cand as usize + 1] {
                got.push(format!(
                    "K {} {} {} {}",
                    c.atnumber, c.type_, c.subtype, c.endpoint
                ));
            }
            for g in &tgi.t_group.as_ref().unwrap()[..tgi.num_t_groups] {
                let nums: Vec<String> = g.num.iter().map(|x| format!(" {x}")).collect();
                let das: Vec<String> = g.num_da.iter().map(|x| format!(" {x}")).collect();
                got.push(format!(
                    "T {} {} |{} |{}",
                    g.n_group_number,
                    g.n_num_endpoints,
                    nums.concat(),
                    das.concat()
                ));
            }
            got.push(format!("F {}", tgi.b_taut_flags_done));
            got.push(format!(
                "I {} {} {} {} |{} |{}",
                tgi.num_iso_h[0],
                tgi.num_iso_h[1],
                tgi.num_iso_h[2],
                tgi.n_num_isotopic_endpoints,
                tgi.n_isotopic_endpoint_atom_number
                    .iter()
                    .map(|x| format!(" {x}"))
                    .collect::<String>(),
                at.iter()
                    .take(n)
                    .map(|a| format!(" {}", a.c_flags))
                    .collect::<String>()
            ));
            let mut s = String::from("P");
            for a in at.iter().take(n) {
                s += &format!(" {},{}", a.endpoint, a.c_point);
            }
            got.push(s);
            dump_state(&p, &at, n, &mut got);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(seen.keys().any(|&k| k > 0) && seen.len() > 3, "{seen:?}");
}

#[test]
fn tautomer_groups_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::*;
    use rdoxide::external::inchi::ichitaut::{
        CGroup, CGroupInfo, TGroup, TGroupInfo, mark_tautomer_groups,
    };
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn ints(v: &serde_json::Value) -> Vec<i32> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect()
    }
    fn dump_state(p: &BnStruct, at: &[InpAtom], n: usize, got: &mut Vec<String>) {
        got.push(format!(
            "C {} {} {} {}",
            p.num_vertices, p.num_edges, p.tot_st_cap, p.tot_st_flow
        ));
        let mut s = String::from("V");
        for v in p.vert.iter().take(p.num_vertices as usize) {
            let e = &v.st_edge;
            s += &format!(
                " {}:{}:{}:{}:{}:{}",
                v.type_, e.cap, e.cap0, e.flow, e.flow0, v.num_adj_edges
            );
        }
        got.push(s);
        let mut s = String::from("F");
        for e in p.edge.iter().take(p.num_edges as usize) {
            s += &format!(" {}:{}:{}:{}", e.cap, e.cap0, e.flow, e.flow0);
        }
        got.push(s);
        let mut s = String::from("A");
        for a in at.iter().take(n) {
            s += &format!(
                " {},{},{},{},{},{},{},{},{}:",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.at_type,
                a.endpoint
            );
            for j in 0..a.valence as usize {
                s += &format!("{}.", a.bond_type[j]);
            }
        }
        got.push(s);
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_marktaut.json")).unwrap();
    let mut seen = std::collections::HashMap::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let iso = ints(&a["iso"]);
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                c_point: a["vt"].as_u64().unwrap() as u16,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                endpoint: a["ep"].as_u64().unwrap() as u16,
                at_type: a["att"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                b_cut_vertex: a["cv"].as_i64().unwrap_or(0) as i8,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = iso[k] as i8;
            }
            x.set_elname(match x.el_number {
                1 => "H",
                6 => "C",
                7 => "N",
                8 => "O",
                15 => "P",
                16 => "S",
                17 => "Cl",
                35 => "Br",
                _ => "Xx",
            });
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        for (i, a) in at.iter_mut().enumerate() {
            a.n_ring_system = inp["atoms"][i]["rs"].as_u64().unwrap_or((i % 2) as u64) as u16;
        }
        let n = at.len();
        let hs = inp["hs"].as_array().unwrap();
        for h in hs {
            let mut x = InpAtom {
                el_number: 1,
                valence: 1,
                orig_at_number: h[1].as_u64().unwrap() as u16,
                ..Default::default()
            };
            x.neighbor[0] = h[0].as_u64().unwrap() as u16;
            at.push(x);
        }
        let (mut p, _) = allocate_and_init_bn_struct(
            &mut at,
            n,
            inp["ma"].as_i64().unwrap() as i32,
            inp["me"].as_i64().unwrap() as i32,
            inp["mp"].as_i64().unwrap() as i32,
        );
        for i in 0..n {
            let vt = u32::from(at[i].c_point);
            p.vert[i].type_ |= (((vt & 0x1c) << 1) | ((vt & 3) << 9)) as u16;
        }
        let tf_packed = inp["tf"].as_u64().unwrap();
        p.pb_taut_flags = (tf_packed & 0xffff) as u32;
        for o in inp["ops"].as_array().unwrap() {
            let o = ints(o);
            let (op, a, b, x, y, z) = (o[0], o[1], o[2], o[3], o[4], o[5]);
            let mut dots = 5;
            let nv = p.num_vertices;
            if op == 10 {
                if a >= nv {
                    continue;
                }
                add_new_vertex(&mut p, a, x, y, z, &mut dots);
            } else if op == 11 {
                if a < p.num_atoms || a >= nv {
                    continue;
                }
                p.vert[a as usize].type_ = x as u16;
            } else {
                if a >= nv || b >= nv || a == b {
                    continue;
                }
                add_new_edge(&mut p, a, b, x, y);
            }
        }
        let mut bd = allocate_and_init_bn_data(p.max_vertices);
        let len = 2 * p.max_vertices as usize + 8;
        bd.rad_endpoints = vec![0; len];
        bd.rad_edges = vec![0; len];
        let mut got: Vec<String> = Vec::new();
        dump_state(&p, &at, n, &mut got);
        let ni = n as i32;
        bns_adjust_flow_bonds_rad(&mut p, &mut bd, &mut at, ni);
        dump_state(&p, &at, n, &mut got);
        let max_cp = at.iter().take(n).map(|a| a.c_point).max().unwrap_or(0) as usize;
        let max_ep = at.iter().take(n).map(|a| a.endpoint).max().unwrap_or(0) as usize;
        let mut cgi = CGroupInfo {
            num_c_groups: max_cp,
            max_num_c_groups: max_cp,
            ..Default::default()
        };
        let mut groups = vec![CGroup::default(); max_cp + 1];
        for g in 1..=max_cp {
            groups[g - 1].n_group_number = g as u16;
            for a in at.iter().take(n) {
                if a.c_point as usize == g {
                    groups[g - 1].num_c_points += 1;
                    groups[g - 1].num[0] += u16::from(a.charge == 1);
                }
            }
        }
        cgi.c_group = Some(groups);
        let fch = (tf_packed >> 16) & 0xff;
        let fresh = max_ep == 0 && fch & 32 != 0;
        let mut t_group = vec![TGroup::default(); max_ep + 5];
        for g in 1..=max_ep {
            t_group[g - 1].n_group_number = g as u16;
            for a in at.iter().take(n) {
                if a.endpoint as usize == g {
                    t_group[g - 1].n_num_endpoints += 1;
                    t_group[g - 1].num[0] += u16::from(a.num_h > 0);
                    t_group[g - 1].num[1] += u16::from(a.charge == -1);
                }
            }
        }
        let mut tgi = TGroupInfo {
            t_group: if fresh { None } else { Some(t_group) },
            num_t_groups: if fresh { 0 } else { max_ep },
            max_num_t_groups: if fresh { 0 } else { max_ep + 3 },
            b_taut_flags: inp["tg"].as_u64().unwrap(),
            b_taut_flags_done: if fch & 2 != 0 { 0x2000 } else { 0 },
            ..Default::default()
        };
        tgi.tni.b_normalization_flags = [0u64, 0x2, 0x10, 0x20][((fch >> 2) & 3) as usize];
        let cgi_has = max_cp > 0;
        for _ in inp["calls"].as_array().unwrap() {
            let ret = mark_tautomer_groups(
                &mut at,
                ni,
                &mut tgi,
                if cgi_has { Some(&mut cgi) } else { None },
                &mut p,
                &mut bd,
            );
            *seen.entry(ret).or_insert(0) += 1;
            got.push(format!("R {ret} {}", tgi.num_t_groups));
            for g in tgi.t_group.iter().flatten().take(tgi.num_t_groups) {
                let nums: Vec<String> = g.num.iter().map(|x| format!(" {x}")).collect();
                let das: Vec<String> = g.num_da.iter().map(|x| format!(" {x}")).collect();
                got.push(format!(
                    "T {} {} |{} |{}",
                    g.n_group_number,
                    g.n_num_endpoints,
                    nums.concat(),
                    das.concat()
                ));
            }
            got.push(format!("F {}", tgi.b_taut_flags_done));
            let mut s = String::from("P");
            for a in at.iter().take(n) {
                s += &format!(" {},{}", a.endpoint, a.c_point);
            }
            got.push(s);
            dump_state(&p, &at, n, &mut got);
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(seen.keys().any(|&k| k > 0) && seen.len() > 3, "{seen:?}");
}

#[test]
fn main_normalization_matches_the_c_code() {
    use rdoxide::external::inchi::ichitaut::{TGroupInfo, mark_alt_bonds_and_taut_groups};
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn dumpat(tag: &str, at: &[InpAtom], n: usize) -> String {
        let mut s = String::from(tag);
        for a in at.iter().take(n) {
            s += &format!(
                " {},{},{},{},{},{},{},{},{},{},{}:",
                a.charge,
                a.num_h,
                a.radical,
                a.chem_bonds_valence,
                a.endpoint,
                a.c_point,
                a.at_type,
                a.num_iso_h[0],
                a.num_iso_h[1],
                a.num_iso_h[2],
                a.c_flags
            );
            for j in 0..a.valence as usize {
                s += &format!("{}.{}.", a.neighbor[j], a.bond_type[j]);
            }
        }
        s
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_markalt.json")).unwrap();
    let mut seen = std::collections::BTreeMap::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom {
                valence: a["v"].as_i64().unwrap() as i8,
                chem_bonds_valence: a["cbv"].as_i64().unwrap() as i8,
                el_number: a["el"].as_u64().unwrap() as u8,
                radical: a["rad"].as_i64().unwrap() as i8,
                charge: a["chg"].as_i64().unwrap() as i8,
                num_h: a["nh"].as_i64().unwrap() as i8,
                orig_at_number: a["orig"].as_u64().unwrap() as u16,
                n_num_at_in_ring_system: a["ring"].as_u64().unwrap() as u16,
                b_cut_vertex: a["cv"].as_i64().unwrap() as i8,
                n_ring_system: a["rs"].as_u64().unwrap() as u16,
                ..Default::default()
            };
            for k in 0..3 {
                x.num_iso_h[k] = a["iso"][k].as_i64().unwrap() as i8;
            }
            x.set_elname(match x.el_number {
                1 => "H",
                6 => "C",
                7 => "N",
                8 => "O",
                15 => "P",
                16 => "S",
                17 => "Cl",
                35 => "Br",
                _ => "Xx",
            });
            for (j, e) in a["nb"].as_array().unwrap().iter().enumerate() {
                x.neighbor[j] = e[0].as_u64().unwrap() as u16;
                x.bond_type[j] = e[1].as_u64().unwrap() as u8;
            }
            at.push(x);
        }
        let n = at.len();
        let nrem = inp["nrem"].as_i64().unwrap() as i32;
        at.resize_with(n + 8, InpAtom::default);
        let mut out: Vec<InpAtom> = vec![InpAtom::default(); n + 8];
        let flags = inp["flags"].as_u64().unwrap();
        let done = inp["done"].as_u64().unwrap();
        let use_tgi = inp["use_tgi"].as_i64().unwrap() != 0;
        let mut tgi = TGroupInfo {
            b_taut_flags: flags,
            b_taut_flags_done: done,
            ..Default::default()
        };
        tgi.tni.n_num_removed_explicit_h = nrem;
        let mut inp_done = done;
        let ret = mark_alt_bonds_and_taut_groups(
            &mut at,
            Some(&mut out),
            n,
            if use_tgi { Some(&mut tgi) } else { None },
            flags,
            &mut inp_done,
            &[],
        );
        *seen.entry(ret.min(100)).or_insert(0) += 1;
        let mut got = vec![format!("R {ret}"), dumpat("A", &at, n)];
        got.push(dumpat("O", &out, (n + nrem as usize).min(n + 8)));
        got.push(format!(
            "F {} {}",
            if use_tgi {
                tgi.b_taut_flags_done
            } else {
                inp_done
            },
            if use_tgi { tgi.b_taut_flags } else { flags }
        ));
        if use_tgi {
            got.push(format!(
                "I {} {} {} {} {} {} {} {} {} | {} |{}",
                tgi.tni.n_num_removed_explicit_h,
                tgi.tni.n_num_removed_protons,
                tgi.tni.n_num_removed_protons_isotopic[0],
                tgi.tni.n_num_removed_protons_isotopic[1],
                tgi.tni.n_num_removed_protons_isotopic[2],
                tgi.tni.b_normalization_flags,
                tgi.num_iso_h[0],
                tgi.num_iso_h[1],
                tgi.num_iso_h[2],
                tgi.n_num_isotopic_endpoints,
                tgi.n_isotopic_endpoint_atom_number
                    .iter()
                    .map(|x| format!(" {x}"))
                    .collect::<String>()
            ));
            got.push(format!("G {}", tgi.num_t_groups));
            for g in tgi.t_group.iter().flatten().take(tgi.num_t_groups) {
                let nums: String = g.num.iter().map(|x| format!(" {x}")).collect();
                let das: String = g.num_da.iter().map(|x| format!(" {x}")).collect();
                got.push(format!(
                    "T {} {} |{} |{}",
                    g.n_group_number, g.n_num_endpoints, nums, das
                ));
            }
        }
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(seen.len() > 6, "{seen:?}");
}

#[test]
fn tautomer_group_counting_matches_the_c_code() {
    use rdoxide::external::inchi::ichitaut::{
        TGroup, TGroupInfo, count_tautomer_groups, set_tautomer_iso_sort_keys,
        sort_tautomer_groups_and_endpoints,
    };
    use rdoxide::external::inchi::runichi3::TG_FLAG_FOUND_ISOTOPIC_H_DONE;
    fn join(v: &[u16]) -> String {
        v.iter().map(|x| format!(" {x}")).collect()
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_tgcount.json")).unwrap();
    let mut seen = std::collections::BTreeMap::new();
    for c in data.as_array().unwrap() {
        let inp = &c["in"];
        let mut ep: Vec<u16> = inp["ep"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as u16)
            .collect();
        let n = ep.len();
        let groups = inp["groups"].as_array().unwrap();
        let g = groups.len();
        let mut t_group = vec![TGroup::default(); g + 3];
        for (i, gr) in groups.iter().enumerate() {
            let v: Vec<u64> = gr
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap())
                .collect();
            t_group[i].n_group_number = v[0] as u16;
            t_group[i].n_num_endpoints = v[1] as u16;
            for q in 0..5 {
                t_group[i].num[q] = v[2 + q] as u16;
            }
        }
        let mut tgi = TGroupInfo {
            t_group: Some(t_group),
            num_t_groups: g,
            max_num_t_groups: inp["maxg"].as_u64().unwrap() as usize,
            n_num_isotopic_endpoints: inp["niso"].as_i64().unwrap() as i32,
            b_taut_flags_done: if inp["done"].as_i64().unwrap() != 0 {
                TG_FLAG_FOUND_ISOTOPIC_H_DONE
            } else {
                0
            },
            ..Default::default()
        };
        tgi.tni.b_normalization_flags = inp["norm"].as_u64().unwrap();
        let mut got = Vec::new();
        let w = set_tautomer_iso_sort_keys(&mut tgi);
        let mut s = format!("K {w}");
        for t in &tgi.t_group.as_ref().unwrap()[..g] {
            s += &format!(" {}", t.i_weight);
        }
        got.push(s);
        let ret = count_tautomer_groups(&mut ep, n, &mut tgi);
        *seen.entry(ret.clamp(-1, 1)).or_insert(0) += 1;
        got.push(format!(
            "R {ret} {} {}",
            tgi.num_t_groups, tgi.n_num_endpoints
        ));
        for t in &tgi.t_group.as_ref().unwrap()[..tgi.num_t_groups] {
            let nums: String = t.num.iter().map(|x| format!(" {x}")).collect();
            got.push(format!(
                "T {} {} {} |{}",
                t.n_group_number, t.n_num_endpoints, t.n_first_endpoint_at_no_pos, nums
            ));
        }
        got.push(format!("P{}", join(&ep)));
        got.push(format!(
            "E{}",
            join(
                &tgi.n_endpoint_atom_number[..tgi
                    .n_endpoint_atom_number
                    .len()
                    .min(tgi.n_num_endpoints.max(0) as usize)]
            )
        ));
        got.push(format!("W{}", join(&tgi.t_group_number)));
        let rank: Vec<u16> = inp["rank"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as u16)
            .collect();
        let nat = n + tgi.num_t_groups;
        let sr = sort_tautomer_groups_and_endpoints(&mut tgi, n, nat, &rank);
        got.push(format!("S {sr}"));
        got.push(format!(
            "E{}",
            join(
                &tgi.n_endpoint_atom_number[..tgi
                    .n_endpoint_atom_number
                    .len()
                    .min(tgi.n_num_endpoints.max(0) as usize)]
            )
        ));
        got.push(format!("W{}", join(&tgi.t_group_number)));
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(got, want, "{inp}");
    }
    assert!(seen.len() == 3, "{seen:?}");
}

#[test]
fn stereo_geometry_helpers_match_the_c_code() {
    use rdoxide::external::inchi::ichister::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn f(v: &serde_json::Value) -> f64 {
        f64::from_bits(v.as_u64().unwrap())
    }
    fn v3(v: &serde_json::Value) -> [f64; 3] {
        let a = v.as_array().unwrap();
        [f(&a[0]), f(&a[1]), f(&a[2])]
    }
    fn vs(v: &serde_json::Value) -> Vec<[f64; 3]> {
        v.as_array().unwrap().iter().map(v3).collect()
    }
    fn b(x: f64) -> String {
        x.to_bits().to_string()
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_stergeom.json")).unwrap();
    for (n, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let mut got: Vec<String> = Vec::new();
        match inp["k"].as_u64().unwrap() {
            0 => {
                let v = vs(&inp["v"]);
                let (a, bb, cc) = (v[0], v[1], v[2]);
                got.push(format!("len3 {}", b(len3(&a))));
                got.push(format!("len2 {}", b(len2(&[a[0], a[1]]))));
                for (i, x) in diff3(&a, &bb).iter().enumerate() {
                    got.push(format!("d{i} {}", b(*x)));
                }
                for (i, x) in cross_prod3(&a, &bb).iter().enumerate() {
                    got.push(format!("x{i} {}", b(*x)));
                }
                got.push(format!("dot {}", b(dot_prod3(&a, &bb))));
                let mut s = 0.0;
                got.push(format!("tp {}", b(triple_prod(&a, &bb, &cc, Some(&mut s)))));
                got.push(format!("sine {}", b(s)));
                got.push(format!("tp0 {}", b(triple_prod(&a, &bb, &cc, None))));
            }
            1 => {
                let co = vs(&inp["co"]);
                let r = get_2d_tetrahedral_ambiguity(
                    &co,
                    inp["add"].as_i64().unwrap() != 0,
                    inp["fix"].as_i64().unwrap() != 0,
                    f(&inp["ang"]),
                );
                got.push(format!("T {r}"));
            }
            2 => {
                let co = vs(&inp["co"]);
                let central = v3(&inp["central"]);
                let vmin = f(&inp["vmin"]);
                let minsine = f(&inp["minsine"]);
                let add = inp["add"].as_i64().unwrap() != 0;
                let mut ms = 0.0;
                got.push(format!(
                    "p1 {}",
                    b(triple_prod_and_min_abs_sine(&co, Some(&mut ms)))
                ));
                got.push(format!("ms1 {}", b(ms)));
                got.push(format!(
                    "P3 {}",
                    i32::from(are_3_vect_in_one_plane(&co, minsine))
                ));
                got.push(format!(
                    "P4 {}",
                    i32::from(are_4at_in_one_plane(&co, minsine))
                ));
                let mut amb = 0;
                got.push(format!(
                    "p2 {}",
                    b(triple_prod_and_min_abs_sine2(
                        &co,
                        &central,
                        add,
                        Some(&mut ms),
                        Some(&mut amb),
                        vmin
                    ))
                ));
                got.push(format!("ms2 {}", b(ms)));
                got.push(format!("A {amb}"));
                amb = 5;
                got.push(format!(
                    "p2b {}",
                    b(triple_prod_and_min_abs_sine2(
                        &co,
                        &central,
                        add,
                        Some(&mut ms),
                        Some(&mut amb),
                        vmin
                    ))
                ));
                got.push(format!("ms2b {}", b(ms)));
                got.push(format!("A {amb}"));
                got.push(format!(
                    "p2c {}",
                    b(triple_prod_and_min_abs_sine2(
                        &co,
                        &central,
                        add,
                        Some(&mut ms),
                        None,
                        vmin
                    ))
                ));
                got.push(format!("ms2c {}", b(ms)));
                got.push(format!(
                    "p2d {}",
                    b(triple_prod_and_min_abs_sine2(
                        &co, &central, add, None, None, vmin
                    ))
                ));
            }
            _ => {
                let atoms = inp["atoms"].as_array().unwrap();
                let mut at: Vec<InpAtom> = Vec::new();
                for a in atoms {
                    let mut x = InpAtom::default();
                    let p = v3(&a["xyz"]);
                    x.x = p[0];
                    x.y = p[1];
                    x.z = p[2];
                    x.b_used_0d_parity = a["u0"].as_i64().unwrap() as i8;
                    let nb = a["nb"].as_array().unwrap();
                    x.valence = nb.len() as i8;
                    for (j, q) in nb.iter().enumerate() {
                        x.neighbor[j] = q[0].as_u64().unwrap() as u16;
                        x.bond_stereo[j] = q[1].as_i64().unwrap() as i8;
                    }
                    at.push(x);
                }
                let pe = inp["pe"].as_i64().unwrap() as i32;
                for i in 0..at.len() {
                    for j in 0..at[i].valence as usize {
                        let (z, t) = get_z_coord(&at, i, j, pe);
                        got.push(format!("z {}", b(z)));
                        got.push(format!("Z {t}"));
                    }
                }
                let z = |k: &str| -> [i8; 3] {
                    let a = inp[k].as_array().unwrap();
                    [
                        a[0].as_i64().unwrap() as i8,
                        a[1].as_i64().unwrap() as i8,
                        a[2].as_i64().unwrap() as i8,
                    ]
                };
                let (z1, z2) = (z("z1"), z("z2"));
                got.push(format!(
                    "C {}",
                    triple_prod_char(
                        &at,
                        inp["a1"].as_u64().unwrap() as usize,
                        inp["i1"].as_u64().unwrap() as usize,
                        &z1,
                        inp["a2"].as_u64().unwrap() as usize,
                        inp["i2"].as_u64().unwrap() as usize,
                        &z2
                    )
                ));
                got.push(format!("D {}", dot_prodchar3(&z1, &z2)));
            }
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {n}: {inp}");
    }
}

#[test]
fn stereo_candidate_checks_match_the_c_code() {
    use rdoxide::external::inchi::ichister::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::util::get_periodic_table_number;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_sterdetect.json")).unwrap();
    for (n, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let mut got: Vec<String> = Vec::new();
        if inp["k"].as_u64().unwrap() == 0 {
            let g = |k: &str| inp[k].as_i64().unwrap() as i32;
            got.push(format!(
                "H {}",
                half_stereo_bond_action(g("p"), g("unk") != 0, g("iso") != 0, g("vab"))
            ));
        } else {
            let mut at: Vec<InpAtom> = Vec::new();
            for a in inp["atoms"].as_array().unwrap() {
                let mut x = InpAtom::default();
                let el = a["el"].as_str().unwrap();
                x.set_elname(el);
                x.el_number = get_periodic_table_number(el) as u8;
                x.charge = a["ch"].as_i64().unwrap() as i8;
                x.radical = a["rad"].as_i64().unwrap() as i8;
                x.num_h = a["nh"].as_i64().unwrap() as i8;
                x.endpoint = a["ep"].as_u64().unwrap() as u16;
                x.chem_bonds_valence = a["cv"].as_i64().unwrap() as i8;
                let nb = a["nb"].as_array().unwrap();
                x.valence = nb.len() as i8;
                for (j, q) in nb.iter().enumerate() {
                    x.neighbor[j] = q[0].as_u64().unwrap() as u16;
                    x.bond_type[j] = q[1].as_u64().unwrap() as u8;
                }
                at.push(x);
            }
            let mode = inp["mode"].as_u64().unwrap();
            let nat = at.len();
            for i in 0..nat {
                for pe in [0, 7, 14] {
                    for zz in 0..2 {
                        got.push(format!(
                            "S {}",
                            b_can_inp_atom_be_a_stereo_center(&at, i, pe, zz)
                        ));
                    }
                }
                for req in 0..4 {
                    for nd in 0..3 {
                        got.push(format!(
                            "R {}",
                            i32::from(b_inp_atom_has_requird_neigh(
                                &at,
                                i,
                                req,
                                nd,
                                (req + nd) & 1
                            ))
                        ));
                    }
                }
                let a = &at[i];
                let (e, ch, r) = (a.elname_str(), a.charge, a.radical);
                got.push(format!(
                    "E {} {} {}",
                    i32::from(b_can_atom_have_a_stereo_bond(e, ch, r)),
                    i32::from(b_can_atom_be_middle_allene(e, ch, r)),
                    i32::from(b_can_atom_be_terminal_allene(e, ch, r))
                ));
                got.push(format!("T {}", b_is_suitable_hetero_inp_atom(a)));
                for j in 0..a.valence as usize {
                    got.push(format!(
                        "B {}",
                        get_allowed_stereo_bond_type(i32::from(a.bond_type[j]))
                    ));
                }
            }
            for i in 0..nat {
                let mut cp = at.clone();
                let mut s = format!("O {}", i32::from(b_is_oxide(&mut cp, i)));
                for j in 0..cp[i].valence as usize {
                    s += &format!(" {}", cp[i].bond_type[j]);
                }
                got.push(s);
                let mut cp = at.clone();
                let mut s = format!(
                    "C {}",
                    can_be_a_stereo_bond_with_isotopic_h(&mut cp, i, mode)
                );
                for a in &cp {
                    for j in 0..a.valence as usize {
                        s += &format!(" {}", a.bond_type[j]);
                    }
                }
                got.push(s);
            }
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {n}: {inp}");
    }
}

#[test]
fn half_stereo_bond_parity_matches_the_c_code() {
    use rdoxide::external::inchi::ichister::half_stereo_bond_parity;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_halfsb.json")).unwrap();
    for (n, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let ints = |v: &serde_json::Value| -> Vec<i64> {
            v.as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap())
                .collect()
        };
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom::default();
            x.set_elname(a["el"].as_str().unwrap());
            x.charge = a["ch"].as_i64().unwrap() as i8;
            x.radical = a["rad"].as_i64().unwrap() as i8;
            x.num_h = a["nh"].as_i64().unwrap() as i8;
            x.c_flags = a["cf"].as_i64().unwrap() as i8;
            x.orig_at_number = a["orig"].as_u64().unwrap() as u16;
            x.iso_atw_diff = a["adiff"].as_i64().unwrap() as i8;
            x.b_used_0d_parity = a["u0"].as_i64().unwrap() as i8;
            let xyz = a["xyz"].as_array().unwrap();
            x.x = f64::from_bits(xyz[0].as_u64().unwrap());
            x.y = f64::from_bits(xyz[1].as_u64().unwrap());
            x.z = f64::from_bits(xyz[2].as_u64().unwrap());
            for (j, v) in ints(&a["iso"]).iter().enumerate() {
                x.num_iso_h[j] = *v as i8;
            }
            for (j, v) in ints(&a["sbp"]).iter().enumerate() {
                x.sb_parity[j] = *v as i8;
            }
            for (j, v) in ints(&a["sbo"]).iter().enumerate() {
                x.sb_ord[j] = *v as i8;
            }
            for (j, v) in ints(&a["sno"]).iter().enumerate() {
                x.sn_orig_at_num[j] = *v as u16;
            }
            let nb = a["nb"].as_array().unwrap();
            x.valence = nb.len() as i8;
            for (j, q) in nb.iter().enumerate() {
                x.neighbor[j] = q[0].as_u64().unwrap() as u16;
                x.bond_stereo[j] = q[1].as_i64().unwrap() as i8;
            }
            at.push(x);
        }
        let nmain = inp["n"].as_u64().unwrap() as usize;
        let nr = inp["nr"].as_u64().unwrap() as usize;
        let z = ints(&inp["z"]);
        let mut zd = [z[0] as i8, z[1] as i8, z[2] as i8];
        let usez = inp["usez"].as_i64().unwrap() != 0;
        let r = half_stereo_bond_parity(
            &mut at,
            0,
            nmain,
            nr,
            if usez { Some(&mut zd) } else { None },
            inp["pe"].as_i64().unwrap() as i32,
            inp["vab"].as_i64().unwrap() as i32,
        );
        let got = format!(
            "{} {} {} {} {} {} 0",
            r,
            if usez { zd[0] } else { 0 },
            if usez { zd[1] } else { 0 },
            if usez { zd[2] } else { 0 },
            at[0].b_ambiguous_stereo,
            at[0].b_used_0d_parity
        );
        assert_eq!(got, c["out"].as_str().unwrap(), "case {n}");
    }
}

#[test]
fn set_stereo_bonds_parity_matches_the_c_code() {
    use rdoxide::external::inchi::ichiring::Queue;
    use rdoxide::external::inchi::ichister::{SpAtom, set_stereo_bonds_parity};
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::util::get_periodic_table_number;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_sbparity.json")).unwrap();
    for (n, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let ints = |v: &serde_json::Value| -> Vec<i64> {
            v.as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap())
                .collect()
        };
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom::default();
            let el = a["el"].as_str().unwrap();
            x.set_elname(el);
            x.el_number = get_periodic_table_number(el) as u8;
            x.charge = a["ch"].as_i64().unwrap() as i8;
            x.radical = a["rad"].as_i64().unwrap() as i8;
            x.num_h = a["nh"].as_i64().unwrap() as i8;
            x.c_flags = a["cf"].as_i64().unwrap() as i8;
            x.orig_at_number = a["orig"].as_u64().unwrap() as u16;
            x.iso_atw_diff = a["adiff"].as_i64().unwrap() as i8;
            x.b_used_0d_parity = a["u0"].as_i64().unwrap() as i8;
            x.chem_bonds_valence = a["cv"].as_i64().unwrap() as i8;
            x.endpoint = a["ep"].as_u64().unwrap() as u16;
            x.c_point = a["cp"].as_u64().unwrap() as u16;
            x.n_ring_system = a["rs"].as_u64().unwrap() as u16;
            let xyz = a["xyz"].as_array().unwrap();
            x.x = f64::from_bits(xyz[0].as_u64().unwrap());
            x.y = f64::from_bits(xyz[1].as_u64().unwrap());
            x.z = f64::from_bits(xyz[2].as_u64().unwrap());
            for (j, v) in ints(&a["iso"]).iter().enumerate() {
                x.num_iso_h[j] = *v as i8;
            }
            for (j, v) in ints(&a["sbp"]).iter().enumerate() {
                x.sb_parity[j] = *v as i8;
            }
            for (j, v) in ints(&a["sbo"]).iter().enumerate() {
                x.sb_ord[j] = *v as i8;
            }
            for (j, v) in ints(&a["sno"]).iter().enumerate() {
                x.sn_orig_at_num[j] = *v as u16;
            }
            let nb = a["nb"].as_array().unwrap();
            x.valence = nb.len() as i8;
            for (j, q) in nb.iter().enumerate() {
                x.neighbor[j] = q[0].as_u64().unwrap() as u16;
                x.bond_type[j] = q[1].as_u64().unwrap() as u8;
                x.bond_stereo[j] = q[2].as_i64().unwrap() as i8;
            }
            at.push(x);
        }
        let nmain = inp["n"].as_u64().unwrap() as usize;
        let nr = inp["nr"].as_u64().unwrap() as usize;
        let tot = nmain + nr;
        let mut out = vec![SpAtom::default(); tot + 1];
        let mut q = Queue::create(tot as i32 + 1).unwrap();
        let mut lev = vec![0u16; tot + 1];
        let mut src = vec![0i8; tot + 1];
        let mut got: Vec<String> = Vec::new();
        for a1 in 0..nmain {
            let r = set_stereo_bonds_parity(
                &mut out,
                &mut at,
                a1,
                nmain,
                nr,
                inp["mode"].as_u64().unwrap(),
                &mut q,
                &mut lev,
                &mut src,
                inp["minring"].as_u64().unwrap() as u16,
                inp["pe"].as_i64().unwrap() as i32,
                inp["vab"].as_i64().unwrap() as i32,
            );
            got.push(format!("R {r}"));
        }
        for i in 0..nmain {
            let o = &out[i];
            let mut s = String::from("A");
            for j in 0..3 {
                s += &format!(
                    " {} {} {} {} {} {} {} {}",
                    o.stereo_bond_neighbor[j],
                    o.stereo_bond_neighbor2[j],
                    o.stereo_bond_ord[j],
                    o.stereo_bond_ord2[j],
                    o.stereo_bond_z_prod[j],
                    o.stereo_bond_z_prod2[j],
                    o.stereo_bond_parity[j],
                    o.stereo_bond_parity2[j]
                );
            }
            s += &format!(
                " {} {} {} {} {} {} 0 | {} {}",
                o.parity,
                o.parity2,
                o.b_ambiguous_stereo,
                o.z_dir[0],
                o.z_dir[1],
                o.z_dir[2],
                at[i].b_ambiguous_stereo,
                at[i].b_used_0d_parity
            );
            got.push(s);
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {n}");
    }
}

#[test]
fn set_stereo_parity_matches_the_c_code() {
    use rdoxide::external::inchi::ichister::{SpAtom, set_stereo_parity};
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::util::get_periodic_table_number;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_stereo.json")).unwrap();
    for (n, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let ints = |v: &serde_json::Value| -> Vec<i64> {
            v.as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap())
                .collect()
        };
        let mut at: Vec<InpAtom> = Vec::new();
        for a in inp["atoms"].as_array().unwrap() {
            let mut x = InpAtom::default();
            let el = a["el"].as_str().unwrap();
            x.set_elname(el);
            x.el_number = get_periodic_table_number(el) as u8;
            x.charge = a["ch"].as_i64().unwrap() as i8;
            x.radical = a["rad"].as_i64().unwrap() as i8;
            x.num_h = a["nh"].as_i64().unwrap() as i8;
            x.c_flags = a["cf"].as_i64().unwrap() as i8;
            x.orig_at_number = a["orig"].as_u64().unwrap() as u16;
            x.iso_atw_diff = a["adiff"].as_i64().unwrap() as i8;
            x.b_used_0d_parity = a["u0"].as_i64().unwrap() as i8;
            x.chem_bonds_valence = a["cv"].as_i64().unwrap() as i8;
            x.endpoint = a["ep"].as_u64().unwrap() as u16;
            x.c_point = a["cp"].as_u64().unwrap() as u16;
            x.n_ring_system = a["rs"].as_u64().unwrap() as u16;
            x.n_num_at_in_ring_system = a["nrs"].as_u64().unwrap() as u16;
            x.p_parity = a["pp"].as_i64().unwrap() as i8;
            let xyz = a["xyz"].as_array().unwrap();
            x.x = f64::from_bits(xyz[0].as_u64().unwrap());
            x.y = f64::from_bits(xyz[1].as_u64().unwrap());
            x.z = f64::from_bits(xyz[2].as_u64().unwrap());
            for (j, v) in ints(&a["iso"]).iter().enumerate() {
                x.num_iso_h[j] = *v as i8;
            }
            for (j, v) in ints(&a["sbp"]).iter().enumerate() {
                x.sb_parity[j] = *v as i8;
            }
            for (j, v) in ints(&a["sbo"]).iter().enumerate() {
                x.sb_ord[j] = *v as i8;
            }
            for (j, v) in ints(&a["sno"]).iter().enumerate() {
                x.sn_orig_at_num[j] = *v as u16;
            }
            for (j, v) in ints(&a["po"]).iter().enumerate() {
                x.p_orig_at_num[j] = *v as u16;
            }
            let nb = a["nb"].as_array().unwrap();
            x.valence = nb.len() as i8;
            for (j, q) in nb.iter().enumerate() {
                x.neighbor[j] = q[0].as_u64().unwrap() as u16;
                x.bond_type[j] = q[1].as_u64().unwrap() as u8;
                x.bond_stereo[j] = q[2].as_i64().unwrap() as i8;
            }
            at.push(x);
        }
        let nmain = inp["n"].as_u64().unwrap() as usize;
        let nr = inp["nr"].as_u64().unwrap() as usize;
        let mut out = vec![SpAtom::default(); nmain + nr + 1];
        let (mut ma, mut mb) = (-1, -1);
        let r = set_stereo_parity(
            &mut at,
            &mut out,
            nmain,
            nr,
            Some(&mut ma),
            Some(&mut mb),
            inp["mode"].as_u64().unwrap(),
            inp["pe"].as_i64().unwrap() as i32,
            inp["vab"].as_i64().unwrap() as i32,
            inp["loose"].as_i64().unwrap() != 0,
            inp["zz"].as_i64().unwrap() as i32,
        );
        let mut got: Vec<String> = vec![format!("R {r} {ma} {mb}")];
        for i in 0..nmain {
            let o = &out[i];
            let mut s = String::from("A");
            for j in 0..3 {
                s += &format!(
                    " {} {} {} {} {} {} {} {}",
                    o.stereo_bond_neighbor[j],
                    o.stereo_bond_neighbor2[j],
                    o.stereo_bond_ord[j],
                    o.stereo_bond_ord2[j],
                    o.stereo_bond_z_prod[j],
                    o.stereo_bond_z_prod2[j],
                    o.stereo_bond_parity[j],
                    o.stereo_bond_parity2[j]
                );
            }
            s += &format!(
                " {} {} {} {} {} {} {} {} | {} {}",
                o.parity,
                o.parity2,
                o.stereo_atom_parity,
                o.stereo_atom_parity2,
                o.b_ambiguous_stereo,
                o.z_dir[0],
                o.z_dir[1],
                o.z_dir[2],
                at[i].b_ambiguous_stereo,
                at[i].b_used_0d_parity
            );
            got.push(s);
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {n}");
    }
}

#[test]
fn sp_atom_parity_bookkeeping_matches_the_c_code() {
    use rdoxide::external::inchi::ichister::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_cansparity.json")).unwrap();
    fn dump(at: &[SpAtom]) -> String {
        let mut s = String::new();
        for a in at {
            s += &format!(" |{} {} {}", a.parity, a.stereo_atom_parity, a.final_parity);
            for j in 0..3 {
                s += &format!(
                    " {} {} {} {}",
                    a.stereo_bond_neighbor[j],
                    a.stereo_bond_ord[j],
                    a.stereo_bond_z_prod[j],
                    a.stereo_bond_parity[j]
                );
            }
        }
        s
    }
    for (n, c) in data.as_array().unwrap().iter().enumerate() {
        let atoms = c["in"]["atoms"].as_array().unwrap();
        let num = atoms.len();
        let mut at = vec![SpAtom::default(); num + 1];
        let mut rank = vec![0u16; num + 1];
        let mut crank = vec![0u16; num + 1];
        for (i, a) in atoms.iter().enumerate() {
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            at[i].parity = a["par"].as_i64().unwrap() as i8;
            at[i].stereo_atom_parity = a["sap"].as_i64().unwrap() as i8;
            at[i].final_parity = a["fp"].as_i64().unwrap() as i8;
            rank[i] = a["rank"].as_u64().unwrap() as u16;
            crank[i] = a["crank"].as_u64().unwrap() as u16;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
            for (j, s) in a["sb"].as_array().unwrap().iter().enumerate() {
                let s: Vec<i64> = s
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_i64().unwrap())
                    .collect();
                at[i].stereo_bond_neighbor[j] = s[0] as u16;
                at[i].stereo_bond_ord[j] = s[1] as i8;
                at[i].stereo_bond_z_prod[j] = s[2] as i8;
                at[i].stereo_bond_parity[j] = s[3] as i8;
            }
        }
        let mut got: Vec<String> = Vec::new();
        for i in 0..num {
            let mut vis = vec![0i8; num + 1];
            got.push(format!(
                "c {} {} {}",
                get_stereo_center_parity(&at, i, &rank),
                i32::from(find_atoms_with_parity(&at, &mut vis, -1, i)),
                get_permutation_parity(&at[i], (i % 3) as u16, &crank)
            ));
            for j in 0..3 {
                let mut s = format!(
                    "h {}",
                    half_stereo_bond_parity_ranks(&at, i, j as i32, &rank)
                );
                if at[i].stereo_bond_neighbor[j] != 0 {
                    let m = at[i].stereo_bond_neighbor[j] as usize - 1;
                    s += &format!(
                        " {} {}",
                        get_stereo_neighbor_pos(&at, i, m),
                        get_stereo_bond_parity(&at, i, m, &rank)
                    );
                }
                got.push(s);
            }
            got.push(format!(
                "p {}",
                get_stereo_neighbor_pos(&at, i, (i + 1) % num)
            ));
        }
        for i in 0..num {
            for j in 0..3 {
                let has = at[i].stereo_bond_neighbor[j] != 0;
                let mut cp = at.clone();
                let r = if has {
                    i32::from(set_one_stereo_bond_ill_def_parity(
                        &mut cp,
                        i,
                        j,
                        1 + ((i + j) % 4) as i32,
                    ))
                } else {
                    -9
                };
                got.push(format!("s {r}{}", dump(&cp[..num])));
                let mut cp = at.clone();
                let r = if has {
                    i32::from(remove_one_stereo_bond(&mut cp, i, j))
                } else {
                    -9
                };
                got.push(format!("r {r}{}", dump(&cp[..num])));
                let mut cp = at.clone();
                let r = i32::from(set_half_stereo_bond_ill_def_pariy(
                    &mut cp,
                    i,
                    j,
                    3 + j as i32,
                ));
                got.push(format!("sh {r}{}", dump(&cp[..num])));
                let mut cp = at.clone();
                let r = i32::from(remove_half_stereo_bond(&mut cp, i, j));
                got.push(format!("rh {r}{}", dump(&cp[..num])));
            }
            let mut cp = at.clone();
            let r = i32::from(remove_one_stereo_center(&mut cp, i));
            got.push(format!("rc {r}{}", dump(&cp[..num])));
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {n}");
    }
}

#[test]
fn canonical_stereo_descriptors_match_the_c_code() {
    use rdoxide::external::inchi::ichister::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_canstereo.json")).unwrap();
    fn dump(at: &[SpAtom]) -> String {
        let mut s = String::new();
        for a in at {
            s += &format!(
                " |{} {} {} {} {} {} {} {}",
                a.parity,
                a.stereo_atom_parity,
                a.final_parity,
                a.b_has_stereo_or_equ_to_stereo,
                a.parity2,
                a.stereo_atom_parity2,
                a.final_parity2,
                a.b_has_stereo_or_equ_to_stereo2
            );
            for j in 0..3 {
                s += &format!(
                    " {} {} {} {} {} {} {} {}",
                    a.stereo_bond_neighbor[j],
                    a.stereo_bond_ord[j],
                    a.stereo_bond_z_prod[j],
                    a.stereo_bond_parity[j],
                    a.stereo_bond_neighbor2[j],
                    a.stereo_bond_ord2[j],
                    a.stereo_bond_z_prod2[j],
                    a.stereo_bond_parity2[j]
                );
            }
        }
        s
    }
    fn dumpcs(s: &CanonStatStereo) -> String {
        let mut t = format!("{} {}", s.n_len_carb, s.n_len_dble);
        for i in 0..s.n_len_carb as usize {
            let c = s.carb.as_ref().unwrap()[i];
            t += &format!(" c{},{}", c.at_num, c.parity);
        }
        for i in 0..s.n_len_dble as usize {
            let d = s.dble.as_ref().unwrap()[i];
            t += &format!(" d{},{},{}", d.at_num1, d.at_num2, d.parity);
        }
        t
    }
    for (n, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let atoms = inp["atoms"].as_array().unwrap();
        let num = atoms.len();
        let mut at = vec![SpAtom::default(); num + 1];
        let (mut rank, mut crank, mut anum, mut acanon) = (
            vec![0u16; num + 1],
            vec![0u16; num + 1],
            vec![0u16; num + 1],
            vec![0u16; num + 1],
        );
        for (i, a) in atoms.iter().enumerate() {
            let g = |k: &str| a[k].as_i64().unwrap();
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            at[i].parity = g("par") as i8;
            at[i].stereo_atom_parity = g("sap") as i8;
            at[i].final_parity = g("fp") as i8;
            at[i].endpoint = g("ep") as u16;
            at[i].num_h = g("nh") as i8;
            for j in 0..3 {
                at[i].num_iso_h[j] = a["iso"][j].as_i64().unwrap() as i8;
            }
            at[i].n_ring_system = g("rs") as u16;
            at[i].b_has_stereo_or_equ_to_stereo = 1;
            rank[i] = g("rank") as u16;
            crank[i] = g("crank") as u16;
            anum[i] = g("an") as u16;
            acanon[i] = g("ac") as u16;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
            for (j, s) in a["sb"].as_array().unwrap().iter().enumerate() {
                let s: Vec<i64> = s
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_i64().unwrap())
                    .collect();
                at[i].stereo_bond_neighbor[j] = s[0] as u16;
                at[i].stereo_bond_ord[j] = s[1] as i8;
                at[i].stereo_bond_z_prod[j] = s[2] as i8;
                at[i].stereo_bond_parity[j] = s[3] as i8;
                at[i].stereo_bond_neighbor2[j] = if s[3] != 0 { s[0] as u16 } else { 0 };
                at[i].stereo_bond_parity2[j] = s[3] as i8;
            }
            at[i].parity2 = at[i].parity;
            at[i].stereo_atom_parity2 = at[i].stereo_atom_parity;
            at[i].final_parity2 = at[i].final_parity;
        }
        let iso = inp["iso"].as_i64().unwrap() != 0;
        let (maxc, maxd) = (
            inp["maxc"].as_i64().unwrap() as i32,
            inp["maxd"].as_i64().unwrap() as i32,
        );
        let nostereo = inp["nostereo"].as_i64().unwrap() != 0;
        let mut got: Vec<String> = Vec::new();
        let mut w = at.clone();
        got.push(format!(
            "U {}{}",
            unmark_non_stereo(&mut w, num, &rank, &anum, iso),
            dump(&w[..num])
        ));
        let mut s = CanonStatStereo {
            n_max_len_carb: maxc,
            n_max_len_dble: maxd,
            ..Default::default()
        };
        if !nostereo {
            s.carb = Some(vec![AtStereoCarb::default(); maxc as usize + 1]);
            s.dble = Some(vec![AtStereoDble::default(); maxd as usize + 1]);
        }
        let r = fill_all_stereo_descriptors(&w, num, &crank, &acanon, &mut s);
        got.push(format!("F {r} {}", dumpcs(&s)));
        if !nostereo {
            for inv in [false, true] {
                let mut w2 = w.clone();
                let mut s2 = s.clone();
                let mut an2 = vec![0u16; num + 1];
                let r = invert_stereo(&mut w2, num, &crank, &mut an2, &mut s2, inv);
                let mut t = format!("I {r}");
                for x in an2.iter().take(num) {
                    t += &format!(" {x}");
                }
                got.push(format!("{t}{}", dump(&w2[..num])));
                got.push(dumpcs(&s2));
            }
        }
        let mut sw = false;
        switch_atom_stereo_and_isotopic_stereo(&mut w, num, &mut sw);
        got.push(format!("S {}{}", i32::from(sw), dump(&w[..num])));
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {n}");
    }
}

#[test]
fn known_stereo_bond_parities_match_the_c_code() {
    use rdoxide::external::inchi::ichister::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_knownsb.json")).unwrap();
    fn dump(at: &[SpAtom]) -> String {
        let mut s = String::new();
        for a in at {
            s += &format!(
                " |{} {} {}",
                a.parity, a.stereo_atom_parity, a.b_has_stereo_or_equ_to_stereo
            );
            for j in 0..3 {
                s += &format!(
                    " {} {} {} {}",
                    a.stereo_bond_neighbor[j],
                    a.stereo_bond_ord[j],
                    a.stereo_bond_z_prod[j],
                    a.stereo_bond_parity[j]
                );
            }
        }
        s
    }
    for (n, c) in data.as_array().unwrap().iter().enumerate() {
        let atoms = c["in"]["atoms"].as_array().unwrap();
        let num = atoms.len();
        let mut at = vec![SpAtom::default(); num + 1];
        let (mut rank, mut crank, mut anum, mut acanon) = (
            vec![0u16; num + 1],
            vec![0u16; num + 1],
            vec![0u16; num + 1],
            vec![0u16; num + 1],
        );
        for (i, a) in atoms.iter().enumerate() {
            let g = |k: &str| a[k].as_i64().unwrap();
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            at[i].parity = g("par") as i8;
            at[i].stereo_atom_parity = g("sap") as i8;
            at[i].num_h = g("nh") as i8;
            at[i].n_ring_system = g("rs") as u16;
            rank[i] = g("rank") as u16;
            crank[i] = g("crank") as u16;
            anum[i] = g("an") as u16;
            acanon[i] = g("ac") as u16;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
            for (j, s) in a["sb"].as_array().unwrap().iter().enumerate() {
                let s: Vec<i64> = s
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_i64().unwrap())
                    .collect();
                at[i].stereo_bond_neighbor[j] = s[0] as u16;
                at[i].stereo_bond_ord[j] = s[1] as i8;
                at[i].stereo_bond_z_prod[j] = s[2] as i8;
                at[i].stereo_bond_parity[j] = s[3] as i8;
            }
        }
        let mut got: Vec<String> = Vec::new();
        got.push(format!(
            "K {}{}",
            set_known_stereo_bond_parities(&mut at, num, &crank, &rank, &anum),
            dump(&at[..num])
        ));
        got.push(format!(
            "M {}{}",
            mark_known_equal_stereo_bond_parities(&mut at, num, &rank, &anum),
            dump(&at[..num])
        ));
        let mut s = CanonStatStereo {
            n_max_len_carb: 16,
            n_max_len_dble: 16,
            carb: Some(vec![AtStereoCarb::default(); 17]),
            dble: Some(vec![AtStereoDble::default(); 17]),
            ..Default::default()
        };
        let f = fill_all_stereo_descriptors(&at, num, &crank, &acanon, &mut s);
        got.push(format!("F {f} {}", s.n_len_dble));
        got.push(format!(
            "R {}{}",
            remove_known_non_stereo_bond_parities(&mut at, num, &crank, &rank, &mut s),
            dump(&at[..num])
        ));
        let mut l = format!("L {}", s.n_len_dble);
        for i in 0..s.n_len_dble as usize {
            let d = s.dble.as_ref().unwrap()[i];
            l += &format!(" d{},{},{}", d.at_num1, d.at_num2, d.parity);
        }
        got.push(l);
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {n}");
    }
}

#[test]
fn known_stereo_center_parities_match_the_c_code() {
    use rdoxide::external::inchi::ichister::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_knownsc.json")).unwrap();
    fn dump(at: &[SpAtom]) -> String {
        let mut s = String::new();
        for a in at {
            s += &format!(
                " |{} {} {} {}",
                a.parity, a.stereo_atom_parity, a.final_parity, a.b_has_stereo_or_equ_to_stereo
            );
            for j in 0..3 {
                s += &format!(
                    " {} {} {} {}",
                    a.stereo_bond_neighbor[j],
                    a.stereo_bond_ord[j],
                    a.stereo_bond_z_prod[j],
                    a.stereo_bond_parity[j]
                );
            }
        }
        s
    }
    fn dumpcs(s: &CanonStatStereo) -> String {
        let mut t = format!("L {} {}", s.n_len_carb, s.n_len_dble);
        for i in 0..s.n_len_carb as usize {
            let c = s.carb.as_ref().unwrap()[i];
            t += &format!(" c{},{}", c.at_num, c.parity);
        }
        for i in 0..s.n_len_dble as usize {
            let d = s.dble.as_ref().unwrap()[i];
            t += &format!(" d{},{},{}", d.at_num1, d.at_num2, d.parity);
        }
        t
    }
    fn newcs() -> CanonStatStereo {
        CanonStatStereo {
            n_max_len_carb: 16,
            n_max_len_dble: 16,
            carb: Some(vec![AtStereoCarb::default(); 17]),
            dble: Some(vec![AtStereoDble::default(); 17]),
            ..Default::default()
        }
    }
    for (n, c) in data.as_array().unwrap().iter().enumerate() {
        let atoms = c["in"]["atoms"].as_array().unwrap();
        let num = atoms.len();
        let mut at = vec![SpAtom::default(); num + 1];
        let (mut rank, mut crank, mut anum, mut acanon) = (
            vec![0u16; num + 1],
            vec![0u16; num + 1],
            vec![0u16; num + 1],
            vec![0u16; num + 1],
        );
        for (i, a) in atoms.iter().enumerate() {
            let g = |k: &str| a[k].as_i64().unwrap();
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            at[i].parity = g("par") as i8;
            at[i].stereo_atom_parity = g("sap") as i8;
            at[i].num_h = g("nh") as i8;
            at[i].n_ring_system = g("rs") as u16;
            rank[i] = g("rank") as u16;
            crank[i] = g("crank") as u16;
            anum[i] = g("an") as u16;
            acanon[i] = g("ac") as u16;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
            for (j, s) in a["sb"].as_array().unwrap().iter().enumerate() {
                let s: Vec<i64> = s
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_i64().unwrap())
                    .collect();
                at[i].stereo_bond_neighbor[j] = s[0] as u16;
                at[i].stereo_bond_ord[j] = s[1] as i8;
                at[i].stereo_bond_z_prod[j] = s[2] as i8;
                at[i].stereo_bond_parity[j] = s[3] as i8;
            }
        }
        let iso = c["in"]["iso"].as_i64().unwrap() != 0;
        let mut got: Vec<String> = Vec::new();
        let mut w = at.clone();
        got.push(format!(
            "K {}{}",
            set_known_stereo_center_parities(&mut w, num, &crank, &rank, &anum),
            dump(&w[..num])
        ));
        got.push(format!(
            "M {}{}",
            mark_known_equal_stereo_center_parities(&mut w, num, &rank, &anum),
            dump(&w[..num])
        ));
        let mut s = newcs();
        got.push(format!(
            "F {}",
            fill_all_stereo_descriptors(&w, num, &crank, &acanon, &mut s)
        ));
        got.push(format!(
            "R {}{}",
            remove_known_non_stereo_center_parities(&mut w, num, &crank, &rank, &mut s),
            dump(&w[..num])
        ));
        got.push(dumpcs(&s));
        let mut w = at.clone();
        let mut s = newcs();
        got.push(format!(
            "O {}{}",
            fill_out_stereo_parities(&mut w, num, &crank, &acanon, &rank, &anum, &mut s, iso),
            dump(&w[..num])
        ));
        got.push(dumpcs(&s));
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {n}");
    }
}

#[test]
fn rank_refinement_matches_the_c_code() {
    use rdoxide::external::inchi::ichimap::*;
    use rdoxide::external::inchi::ichister::SpAtom;
    use rdoxide::external::inchi::ichitaut::{TGroup, TGroupInfo};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_ranks.json")).unwrap();
    fn pnl(nl: &NeighList, n: usize) -> String {
        let mut s = String::new();
        for l in nl.iter().take(n) {
            s += " [";
            s += &l
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",");
            s += "]";
        }
        s
    }
    fn pa(t: &str, a: &[u16], n: usize) -> String {
        let mut s = t.to_string();
        for x in a.iter().take(n) {
            s += &format!(" {x}");
        }
        s
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let n = inp["n"].as_u64().unwrap() as usize;
        let ntg = inp["ntg"].as_u64().unwrap() as usize;
        let u16s = |v: &serde_json::Value| -> Vec<u16> {
            v.as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap() as u16)
                .collect()
        };
        let mut rank = u16s(&inp["rank"]);
        let mut anum = u16s(&inp["anum"]);
        let mut symm = u16s(&inp["symm"]);
        rank.resize(n + 8, 0);
        anum.resize(n + 8, 0);
        symm.resize(n + 8, 0);
        let endpoint = u16s(&inp["endpoint"]);
        let mut at = vec![SpAtom::default(); n + 1];
        for i in 0..n {
            let adj = u16s(&inp["adj"][i]);
            let bonds = u16s(&inp["bonds"][i]);
            at[i].valence = adj.len() as i8;
            at[i].endpoint = endpoint[i];
            for (j, q) in adj.iter().enumerate() {
                at[i].neighbor[j] = *q;
                at[i].bond_type[j] = bonds[j] as u8;
            }
        }
        let mut tgi = TGroupInfo::default();
        let mut groups = Vec::new();
        let mut ends: Vec<u16> = Vec::new();
        for g in inp["groups"].as_array().unwrap() {
            let e = u16s(g);
            let t = TGroup {
                n_num_endpoints: e.len() as u16,
                n_first_endpoint_at_no_pos: ends.len() as u16,
                ..Default::default()
            };
            ends.extend_from_slice(&e);
            groups.push(t);
        }
        tgi.t_group = Some(groups);
        tgi.num_t_groups = ntg;
        tgi.n_endpoint_atom_number = ends;
        let lct = u16s(&inp["lct"]);
        let invs: Vec<AtomInvariant2> = inp["inv"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| {
                let mut a = AtomInvariant2::default();
                for (j, x) in u16s(&v["val"]).iter().enumerate() {
                    a.val[j] = *x;
                }
                a.iso_sort_key = v["k"].as_i64().unwrap();
                a.iso_aux_key = v["ak"].as_i64().unwrap() as i8;
                a
            })
            .collect();
        let mut got: Vec<String> = Vec::new();
        for sq in 0..2 {
            let nl = create_neigh_list(n, n, &at, sq == 1, Some(&tgi));
            got.push(format!("C{sq}{}", pnl(&nl, n)));
            if ntg != 0 {
                let nt = create_neigh_list(n, n + ntg, &at, sq == 1, Some(&tgi));
                got.push(format!("T{sq}{}", pnl(&nt, n + ntg)));
            }
        }
        match create_neigh_list_from_linear_ct(&lct, lct.len(), n) {
            Some(nl) => got.push(format!("L{}", pnl(&nl, n))),
            None => got.push("L null".to_string()),
        }
        let base = create_neigh_list(n, n, &at, false, Some(&tgi));
        {
            let mut r2 = vec![0u16; n + 8];
            let (nd, ch) = sorted_equ_info_to_ranks(&symm, &mut r2, &anum, n);
            got.push(format!("E {nd} {}", i32::from(ch)));
            got.push(pa("e", &r2, n));
            let mut s2 = vec![0u16; n + 8];
            got.push(format!(
                "Q {}",
                sorted_ranks_to_equ_info(&mut s2, &rank, &anum, n)
            ));
            got.push(pa("q", &s2, n));
            let mut inpl = symm.clone();
            let (nd, _) = sorted_equ_info_to_ranks_in_place(&mut inpl, &anum, n);
            got.push(format!("I {nd}{}", pa("i", &inpl, n)));
        }
        {
            let mut nl = base.clone();
            sort_neigh_lists2(n, &rank, &mut nl, &anum);
            got.push(format!("S2{}", pnl(&nl, n)));
            let mut nl = base.clone();
            sort_neigh_lists3(n, &rank, &mut nl, &anum);
            got.push(format!("S3{}", pnl(&nl, n)));
            sort_neigh_lists_by_symm_and_canon_rank(n, &mut nl, &symm, &rank);
            got.push(format!("SB{}", pnl(&nl, n)));
        }
        for v in 0..4 {
            let nl = base.clone();
            let mut an = anum.clone();
            let mut nr = vec![0u16; n + 8];
            let r = match v {
                0 => set_new_ranks_from_neigh_lists3(n, &nl, &rank, &mut nr, &mut an),
                1 => set_new_ranks_from_neigh_lists4(
                    n,
                    &nl,
                    &rank,
                    &mut nr,
                    &mut an,
                    (n / 2 + 1) as u16,
                ),
                _ => set_new_ranks_from_neigh_lists(
                    n,
                    &nl,
                    &rank,
                    &mut nr,
                    &mut an,
                    v - 2,
                    RankSortCmp::NeighListRanksOrd,
                ),
            };
            got.push(format!("R{v} {r}{}", pa("", &nr, n)));
            got.push(pa("a", &an, n));
        }
        for alt in 0..2 {
            let nl = base.clone();
            let mut an = anum.clone();
            let mut nr = vec![0u16; n + 8];
            let r = set_new_ranks_from_neigh_lists(
                n,
                &nl,
                &rank,
                &mut nr,
                &mut an,
                alt,
                RankSortCmp::NeighListRanks,
            );
            got.push(format!("RB{alt} {r}{}", pa("", &nr, n)));
        }
        for v in 0..6 {
            let mut nl = base.clone();
            let mut an = anum.clone();
            let mut cur = rank.clone();
            let mut prev = vec![0u16; n + 8];
            let mut it = 0i64;
            let r = match v {
                0 => differentiate_ranks2(n, &mut nl, 1, &mut cur, &mut prev, &mut an, &mut it, 0),
                1 => differentiate_ranks2(n, &mut nl, 1, &mut cur, &mut prev, &mut an, &mut it, 1),
                2 => differentiate_ranks3(n, &mut nl, 1, &mut cur, &mut prev, &mut an, &mut it),
                3 => differentiate_ranks4(
                    n,
                    &mut nl,
                    1,
                    &mut cur,
                    &mut prev,
                    &mut an,
                    (n / 2 + 1) as u16,
                    &mut it,
                ),
                _ => differentiate_ranks_basic(
                    n,
                    &mut nl,
                    1,
                    &mut cur,
                    &mut prev,
                    &mut an,
                    &mut it,
                    v - 4,
                ),
            };
            got.push(format!("D{v} {r} {it}{}", pa("", &cur, n)));
            got.push(pa("p", &prev, n));
            got.push(pa("a", &an, n));
            got.push(format!("n{}", pnl(&nl, n)));
        }
        {
            let mut nl = base.clone();
            let mut stack = vec![rank[..n].to_vec(), anum[..n].to_vec()];
            let mut tmp = vec![0u16; n + 8];
            let mut cs = CanonCounters::default();
            let r = break_all_ties(n, n + 4, &mut stack, &mut nl, &mut tmp, &mut cs);
            got.push(format!(
                "B {r} {} {}{}",
                cs.l_num_neigh_list_iter,
                cs.l_num_break_ties,
                pa("", &stack[2], n)
            ));
            got.push(pa("b", &stack[3], n));
        }
        {
            let mut s = String::from("V");
            for i in 0..invs.len() {
                for j in 0..invs.len() {
                    s += &format!(
                        " {},{}",
                        comp_atom_invariants2_only(&invs[i], &invs[j]),
                        comp_atom_invariants2(&invs, i as u16, j as u16)
                    );
                }
            }
            got.push(s);
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn stereo_mapping_helpers_match_the_c_code() {
    use rdoxide::external::inchi::ichimap::*;
    use rdoxide::external::inchi::ichister::{AtStereoCarb, AtStereoDble, SpAtom};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_mapstereo.json")).unwrap();
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let atoms = inp["atoms"].as_array().unwrap();
        let n = atoms.len();
        let mut at = vec![SpAtom::default(); n + 1];
        let (mut crank, mut acanon) = (vec![0u16; n + 8], vec![0u16; n + 8]);
        let (mut r1, mut a1, mut r2, mut a2) = (
            vec![0u16; n + 8],
            vec![0u16; n + 8],
            vec![0u16; n + 8],
            vec![0u16; n + 8],
        );
        for (i, a) in atoms.iter().enumerate() {
            let g = |k: &str| a[k].as_i64().unwrap();
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            at[i].parity = g("par") as i8;
            at[i].stereo_atom_parity = g("sap") as i8;
            at[i].num_h = g("nh") as i8;
            r1[i] = g("rank") as u16;
            a1[i] = g("an") as u16;
            r2[i] = g("rank2") as u16;
            a2[i] = g("an2") as u16;
            crank[i] = g("crank") as u16;
            acanon[i] = g("ac") as u16;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
            for (j, s) in a["sb"].as_array().unwrap().iter().enumerate() {
                let s: Vec<i64> = s
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_i64().unwrap())
                    .collect();
                at[i].stereo_bond_neighbor[j] = s[0] as u16;
                at[i].stereo_bond_ord[j] = s[1] as i8;
                at[i].stereo_bond_z_prod[j] = s[2] as i8;
                at[i].stereo_bond_parity[j] = s[3] as i8;
            }
        }
        let st1 = vec![r1.clone(), a1.clone()];
        let st2 = vec![r2.clone(), a2.clone()];
        let mut used = vec![0i8; n + 8];
        set_use_atom_for_stereo(&mut used, &at, n);
        for m in inp["mask"].as_array().unwrap() {
            used[m[0].as_u64().unwrap() as usize] = m[1].as_i64().unwrap() as i8;
        }
        let ni = n as i32;
        let mut got: Vec<String> = Vec::new();
        let join = |v: Vec<String>| v.join(" ");
        got.push(format!(
            "U {}",
            join(used[..n].iter().map(|x| x.to_string()).collect())
        ));
        got.push(format!(
            "A {}",
            join(
                (1..=n as u16)
                    .map(|cr| all_sc_same(cr, &st1, &st2, &acanon, &at).to_string())
                    .collect()
            )
        ));
        got.push(format!(
            "P {}",
            join(
                (1..=n as u16)
                    .map(|cr| match b_unique_at_nbr_from_mapping_rank(&st1, cr) {
                        Some(u) => format!("1,{u}"),
                        None => "0,-1".to_string(),
                    })
                    .collect()
            )
        ));
        {
            let (mut cr1, mut min, mut first) = (0u16, 0u16, true);
            let mut s = String::from("X");
            for _ in 0..12 {
                let ok = next_sc_at_canon_rank2(
                    &mut cr1, &mut min, &mut first, &used, &st1, &st2, &acanon, ni,
                );
                s += &format!(" {},{},{},{}", i32::from(ok), cr1, min, i32::from(first));
                if !ok {
                    break;
                }
            }
            got.push(s);
        }
        {
            let mut v = Vec::new();
            for a in 1..=n as u16 {
                for b in 1..=n as u16 {
                    v.push(all_sb_same(a, b, &st1, &st2, &acanon, &at).to_string());
                }
            }
            got.push(format!("B {}", join(v)));
        }
        for al in 0..2 {
            let (mut cr1, mut cr2, mut m1, mut m2, mut first) = (0u16, 0u16, 0u16, 0u16, true);
            let mut s = format!("Y{al}");
            for _ in 0..12 {
                let ok = next_sb_at_canon_ranks2(
                    &mut cr1,
                    &mut cr2,
                    &mut m1,
                    &mut m2,
                    &mut first,
                    &used,
                    &st1,
                    &st2,
                    &crank,
                    &acanon,
                    &at,
                    ni,
                    al == 1,
                );
                s += &format!(
                    " {},{},{},{},{},{}",
                    i32::from(ok),
                    cr1,
                    cr2,
                    m1,
                    m2,
                    i32::from(first)
                );
                if !ok {
                    break;
                }
            }
            got.push(s);
        }
        {
            let mut v = Vec::new();
            for sp in 0..7 {
                for cp in 0..3 {
                    for m in 0..8 {
                        let (mut p, mut q) = (sp, cp);
                        let r = next_stereo_parity2_test(
                            &mut p,
                            &mut q,
                            m & 1,
                            (m >> 1) & 1,
                            (m >> 2) & 1,
                            (m + 1) & 1,
                            (m + 3) & 1,
                            3,
                        );
                        v.push(format!("{r},{p},{q}"));
                    }
                }
            }
            got.push(format!("Z {}", join(v)));
        }
        {
            let ct: Vec<i64> = inp["ct"]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap())
                .collect();
            let rows = |k: &str| -> Vec<Vec<i64>> {
                inp[k]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|r| {
                        r.as_array()
                            .unwrap()
                            .iter()
                            .map(|x| x.as_i64().unwrap())
                            .collect()
                    })
                    .collect()
            };
            let dbl = |k: &str| -> Vec<AtStereoDble> {
                rows(k)
                    .iter()
                    .map(|r| AtStereoDble {
                        at_num1: r[0] as u16,
                        at_num2: r[1] as u16,
                        parity: r[2] as u8,
                    })
                    .collect()
            };
            let carb = |k: &str| -> Vec<AtStereoCarb> {
                rows(k)
                    .iter()
                    .map(|r| AtStereoCarb {
                        at_num: r[0] as u16,
                        parity: r[1] as u8,
                    })
                    .collect()
            };
            let (d1, d2, c1, c2) = (dbl("d1"), dbl("d2"), carb("c1"), carb("c2"));
            let f = ct[4];
            let pd1 = if f & 1 != 0 { None } else { Some(&d1[..]) };
            let pd2 = if f & 2 != 0 { None } else { Some(&d2[..]) };
            let pc1 = if f & 4 != 0 { None } else { Some(&c1[..]) };
            let pc2 = if f & 8 != 0 { None } else { Some(&c2[..]) };
            let (nd1, nd2, nc1, nc2) = (ct[0] as i32, ct[1] as i32, ct[2] as i32, ct[3] as i32);
            let mut s = format!(
                "K {} {} {}",
                compare_lin_ct_stereo_dble(pd1, nd1, pd2, nd2),
                compare_lin_ct_stereo_carb(pc1, nc1, pc2, nc2),
                compare_lin_ct_stereo(pd1, nd1, pc1, nc1, pd2, nd2, pc2, nc2)
            );
            if nc1 > 0 && nd1 > 0 {
                // (the C harness reads c2[0]/d2[0], which are zero-filled when the list is empty)
                let c20 = c2.first().copied().unwrap_or_default();
                let d20 = d2.first().copied().unwrap_or_default();
                s += &format!(
                    " {} {} {}",
                    compare_lin_ct_stereo_atom_to_values(
                        &c1[0],
                        if c20.at_num != 0 { c20.at_num } else { 1 },
                        c20.parity
                    ),
                    compare_lin_ct_stereo_double_to_values(
                        &d1[0],
                        d20.at_num1,
                        d20.at_num2,
                        d20.parity
                    ),
                    compare_lin_ct_stereo_double_to_values(
                        &d1[0],
                        d1[0].at_num1,
                        d1[0].at_num2,
                        d1[0].parity
                    )
                );
            }
            got.push(s);
        }
        {
            let mut eq: Vec<u16> = inp["eq"]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap() as u16)
                .collect();
            eq.resize(n + 8, 0);
            let mut s = String::from("E");
            for i in 0..n {
                s += &format!(" {}", n_get_mcr(&mut eq, i as u16));
            }
            got.push(s);
            for j in inp["joins"].as_array().unwrap() {
                let mut s = format!(
                    "J {}",
                    n_join2_mcrs(
                        &mut eq,
                        j[0].as_u64().unwrap() as u16,
                        j[1].as_u64().unwrap() as u16
                    )
                );
                for x in eq.iter().take(n) {
                    s += &format!(" {x}");
                }
                got.push(s);
            }
        }
        {
            let symm: Vec<u16> = inp["symm"]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap() as u16)
                .collect();
            let mut symm = symm;
            symm.resize(n + 8, 0);
            let mut ct = CurTree::default();
            let mut saved: Vec<i32> = Vec::new();
            let nn = if n > 0 { n as i32 } else { 1 };
            for o in inp["ops"].as_array().unwrap() {
                let (op, arg) = (o[0].as_i64().unwrap() as i32, o[1].as_i64().unwrap() as i32);
                let r = match op {
                    0 => {
                        saved.clear();
                        ct.alloc(if arg > 2 { arg } else { 3 })
                    }
                    1 => ct.add_rank(arg as u16),
                    2 => ct.add_atom(arg % nn),
                    3 => {
                        let tp = if saved.is_empty() {
                            0
                        } else {
                            saved[arg as usize % saved.len()]
                        };
                        ct.keep_last_atoms_only(tp, 1);
                        tp
                    }
                    4 => ct.remove_if_last_atom(arg % nn),
                    5 => ct.remove_last_rank(),
                    6 => ct.remove_last_rank_if_no_atoms(),
                    7 => {
                        if saved.is_empty() {
                            -7
                        } else {
                            ct.set_pos(saved[arg as usize % saved.len()])
                        }
                    }
                    8 => ct.is_last_rank(arg as u16),
                    9 => {
                        if !ct.tree.is_empty()
                            && ct.cur_len > 1
                            && i32::from(ct.tree[ct.cur_len as usize - 1]) - 1 > ct.cur_len - 1
                        {
                            -9
                        } else {
                            ct.is_last_atom_equ(arg % nn, &symm)
                        }
                    }
                    10 => {
                        ct.free();
                        saved.clear();
                        0
                    }
                    _ => {
                        let p = ct.get_pos();
                        if saved.len() < 64 {
                            saved.push(p);
                        }
                        p
                    }
                };
                if ct.cur_len < 0 {
                    ct.cur_len = 0;
                }
                let mut s = format!("T {op} {r} {} {}:", ct.cur_len, ct.max_len);
                if !ct.tree.is_empty() {
                    for i in 0..ct.cur_len.min(ct.max_len).max(0) as usize {
                        s += &format!(" {}", ct.tree[i]);
                    }
                }
                got.push(s);
            }
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn mapping_of_ties_and_parities_matches_the_c_code() {
    use rdoxide::external::inchi::ichimap::*;
    use rdoxide::external::inchi::ichister::SpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_mapties.json")).unwrap();
    fn pen(e: &EqNeigh) -> String {
        format!(
            " | {} {} {} {} {} {} {} {}",
            e.num_to,
            e.to_at[0],
            e.to_at[1],
            e.to_at[2],
            e.to_at[3],
            e.from_at,
            e.rank,
            e.canon_rank
        )
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let atoms = inp["atoms"].as_array().unwrap();
        let n = atoms.len();
        let nmax = n + 2;
        let mut at = vec![SpAtom::default(); n + 1];
        let (mut crank, mut r1, mut a1, mut r2, mut a2) = (
            vec![0u16; nmax + 8],
            vec![0u16; nmax],
            vec![0u16; nmax],
            vec![0u16; nmax],
            vec![0u16; nmax],
        );
        for (i, a) in atoms.iter().enumerate() {
            let g = |k: &str| a[k].as_i64().unwrap();
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            at[i].parity = g("par") as i8;
            at[i].stereo_atom_parity = g("sap") as i8;
            at[i].num_h = g("nh") as i8;
            at[i].b_has_stereo_or_equ_to_stereo = inp["hs"][i].as_i64().unwrap() as i8;
            r1[i] = g("rank") as u16;
            a1[i] = g("an") as u16;
            r2[i] = g("rank2") as u16;
            a2[i] = g("an2") as u16;
            crank[i] = g("crank") as u16;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
            for (j, s) in a["sb"].as_array().unwrap().iter().enumerate() {
                let s: Vec<i64> = s
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_i64().unwrap())
                    .collect();
                at[i].stereo_bond_neighbor[j] = s[0] as u16;
                at[i].stereo_bond_ord[j] = s[1] as i8;
                at[i].stereo_bond_z_prod[j] = s[2] as i8;
                at[i].stereo_bond_parity[j] = s[3] as i8;
            }
        }
        let mut got: Vec<String> = Vec::new();
        for p in inp["pairs"].as_array().unwrap() {
            let p: Vec<usize> = p
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap() as usize)
                .collect();
            let (f, t, kf, kt, useen) = (p[0], p[1], p[2], p[3], p[4] != 0);
            let mut e1 = EqNeigh::default();
            let mut e2 = EqNeigh::default();
            let ra = parity_of_mapped_atom2(f, t, &at, useen.then_some(&mut e1), &crank, &r1, &r2);
            got.push(format!(
                "A {ra}{}",
                if useen { pen(&e1) } else { String::new() }
            ));
            let rh = parity_of_mapped_half_bond(
                f,
                t,
                at[f].neighbor[kf] as usize,
                at[t].neighbor[kt] as usize,
                &at,
                useen.then_some(&mut e2),
                &crank,
                &r1,
                &r2,
            );
            got.push(format!(
                "H {rh}{}",
                if useen { pen(&e2) } else { String::new() }
            ));
            let ra = parity_of_mapped_atom2(f, f, &at, None, &crank, &r1, &r1);
            let rh = parity_of_mapped_half_bond(
                f,
                f,
                at[f].neighbor[kf] as usize,
                at[f].neighbor[kf] as usize,
                &at,
                None,
                &crank,
                &r1,
                &r1,
            );
            got.push(format!("S {ra} {rh}"));
        }
        let mut s = String::from("G");
        for i in 0..n {
            s += &format!(
                " {} {}",
                get_min_new_rank(&r1, &a1, r1[i]),
                get_min_new_rank(&r2, &a2, r2[i])
            );
        }
        got.push(s);
        let mut s = String::from("M");
        for i in 0..n {
            s += &format!(" {}", might_change_other_atom_parity(&at, n, i, &r2, &r1));
        }
        got.push(s);
        let mut nl = create_neigh_list(n, n, &at, false, None);
        let mut st1 = vec![r1.clone(), a1.clone()];
        let mut st2 = vec![r2.clone(), a2.clone()];
        for _ in 2..6 {
            st1.push(vec![0u16; nmax]);
        }
        for _ in 2..4 {
            st2.push(vec![0u16; nmax]);
        }
        let mut tmp = vec![0u16; nmax + 8];
        let mut cs = CanonCounters::default();
        for m in inp["maps"].as_array().unwrap() {
            let m: Vec<i32> = m
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap() as i32)
                .collect();
            let (x, y, nmr) = (m[0] as usize, m[1] as usize, m[2]);
            let (mut nr, mut badd, mut bm) = (0u16, 0, 0);
            let rt = number_of_ties(&mut st1, &mut st2, nmax, x, y, &mut nr, &mut badd, &mut bm);
            got.push(format!("T {rt} {nr} {badd} {bm}"));
            let mut newn = -1;
            let rr = map_an_atom2(
                n, nmax, x, y, &mut tmp, nmr, &mut newn, &mut cs, &mut nl, &mut st1, &mut st2,
                &mut badd,
            );
            let mut s = format!(
                "R {rr} {newn} {badd} {} {}",
                cs.l_num_neigh_list_iter, cs.l_num_break_ties
            );
            for v in [&st1[2], &st1[3], &st2[2], &st2[3]] {
                s += " |";
                for x in &v[..n] {
                    s += &format!(" {x}");
                }
            }
            s += &format!(" | {} {}", st1[4][0], st1[5][0]);
            got.push(s);
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn neighbor_tie_breaking_and_symmetric_paths_match_the_c_code() {
    use rdoxide::external::inchi::ichimap::*;
    use rdoxide::external::inchi::ichister::SpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_mapbreak.json")).unwrap();
    fn ints(v: &serde_json::Value) -> Vec<i64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect()
    }
    fn pa(a: &[u16], n: usize) -> String {
        let mut s = String::new();
        for x in &a[..n] {
            s += &format!(" {x}");
        }
        s + " |"
    }
    fn pnl(nl: &NeighList, n: usize) -> String {
        let mut s = String::new();
        for l in nl.iter().take(n) {
            s += &format!(
                " [{}]",
                l.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            );
        }
        s
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let atoms = inp["atoms"].as_array().unwrap();
        let n = atoms.len();
        let nmax = n + 2;
        let mut at = vec![SpAtom::default(); n + 1];
        let (mut crank, mut r1, mut r2) = (
            vec![0u16; nmax + 8],
            vec![0u16; nmax + 8],
            vec![0u16; nmax + 8],
        );
        for (i, a) in atoms.iter().enumerate() {
            let g = |k: &str| a[k].as_i64().unwrap();
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            at[i].parity = g("par") as i8;
            at[i].stereo_atom_parity = g("sap") as i8;
            at[i].num_h = g("nh") as i8;
            at[i].n_ring_system = inp["rs"][i].as_u64().unwrap() as u16;
            r1[i] = g("rank") as u16;
            r2[i] = g("rank2") as u16;
            crank[i] = g("crank") as u16;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
            for (j, s) in a["sb"].as_array().unwrap().iter().enumerate() {
                let s = ints(s);
                at[i].stereo_bond_neighbor[j] = s[0] as u16;
                at[i].stereo_bond_ord[j] = s[1] as i8;
                at[i].stereo_bond_z_prod[j] = s[2] as i8;
                at[i].stereo_bond_parity[j] = s[3] as i8;
            }
        }
        let mut got: Vec<String> = Vec::new();
        let v16 = |v: &serde_json::Value| -> Vec<u16> {
            let mut x: Vec<u16> = ints(v).iter().map(|x| *x as u16).collect();
            x.resize(nmax, 0);
            x
        };
        for ch in inp["checks"].as_array().unwrap() {
            let h = ch.as_array().unwrap();
            let i = |k: usize| h[k].as_i64().unwrap() as u16;
            let (v1, v2, o1, o2) = (v16(&h[7]), v16(&h[8]), v16(&h[9]), v16(&h[10]));
            let rr2 = if h[6].as_i64().unwrap() != 0 {
                &r1
            } else {
                &r2
            };
            let r = check_next_symm_neighbors_and_bonds(
                &at,
                i(0),
                i(1),
                i(2),
                i(3),
                &[i(4), i(5)],
                &v1,
                &v2,
                &o1,
                &o2,
                &r1,
                rr2,
            );
            got.push(format!("C {r}"));
        }
        let nlp = create_neigh_list(n, n, &at, false, None);
        for p in inp["paths"].as_array().unwrap() {
            let h = ints(p);
            let (mut v1, mut v2, mut o1, mut o2) = (
                vec![0u16; nmax],
                vec![0u16; nmax],
                vec![0u16; nmax],
                vec![0u16; nmax],
            );
            let mut len = 0u16;
            let mut inv = h[6] as i32;
            let rr2 = if h[8] != 0 { &r1 } else { &r2 };
            let ret = create_check_symm_paths(
                &at,
                h[0] as u16,
                h[1] as u16,
                h[2] as u16,
                h[3] as u16,
                &[h[4] as u16, h[5] as u16],
                &mut v1,
                &mut v2,
                &mut o1,
                &mut o2,
                &nlp,
                &nlp,
                &r1,
                rr2,
                &mut len,
                &mut inv,
                h[7] as i32,
            );
            got.push(format!(
                "P {ret} {len} {inv} |{}{}{}{}",
                pa(&v1, n),
                pa(&v2, n),
                pa(&o1, n),
                pa(&o2, n)
            ));
        }
        for b in inp["breaks"].as_array().unwrap() {
            let h = ints(b);
            let (ib, ia, in1, in2, mode, symm2) = (
                h[0] as usize,
                h[1] as usize,
                h[2] as usize,
                h[3] as usize,
                h[4] as i32,
                h[5],
            );
            let mut nn = [0u16; 4];
            let nv = (at[ia].valence as usize).min(4);
            nn[..nv].copy_from_slice(&at[ia].neighbor[..nv]);
            let mut s1 = vec![vec![0u16; nmax], vec![0u16; nmax]];
            let mut s2 = vec![vec![0u16; nmax], vec![0u16; nmax]];
            for (d, x) in s1[1].iter_mut().zip(&h[6..6 + n]) {
                *d = *x as u16;
            }
            let mut tmp = vec![0u16; nmax + 8];
            let mut nl = create_neigh_list(n, n, &at, false, None);
            let mut nl1 = create_neigh_list(n, n, &at, false, None);
            let mut nl2 = create_neigh_list(n, n, &at, false, None);
            let mut iter = 0i64;
            let sr = if symm2 != 0 { &r2 } else { &r1 };
            let ret = break_neighbors_tie(
                &at, n, n, ib, ia, &nn, in1, in2, mode, &mut s1, &mut s2, &mut tmp, &mut nl, sr,
                &crank, &mut nl1, &mut nl2, &mut iter,
            );
            got.push(format!(
                "B {ret} {iter} |{}{}{}{}",
                pa(&s1[0], n),
                pa(&s1[1], n),
                pa(&s2[0], n),
                pa(&s2[1], n)
            ));
            got.push(pnl(&nl1, n));
            got.push(pnl(&nl2, n));
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn calculated_paths_parity_comparison_matches_the_c_code() {
    use rdoxide::external::inchi::ichimap::*;
    use rdoxide::external::inchi::ichister::{AtStereoCarb, AtStereoDble, SpAtom};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_mapcalc.json")).unwrap();
    fn ints(v: &serde_json::Value) -> Vec<i64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect()
    }
    fn pa(a: &[u16], n: usize) -> String {
        let mut s = String::new();
        for x in &a[..n] {
            s += &format!(" {x}");
        }
        s + " |"
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let atoms = inp["atoms"].as_array().unwrap();
        let n = atoms.len();
        let nmax = n + 2;
        let mut at = vec![SpAtom::default(); n + 1];
        let (mut crank, mut r1, mut acanon) = (
            vec![0u16; nmax + 8],
            vec![0u16; nmax + 8],
            vec![0u16; nmax + 8],
        );
        for (i, a) in atoms.iter().enumerate() {
            let g = |k: &str| a[k].as_i64().unwrap();
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            at[i].parity = g("par") as i8;
            at[i].stereo_atom_parity = g("sap") as i8;
            at[i].num_h = g("nh") as i8;
            at[i].n_ring_system = inp["rs"][i].as_u64().unwrap() as u16;
            r1[i] = g("rank") as u16;
            crank[i] = g("crank") as u16;
            acanon[i] = g("ac") as u16;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
            for (j, s) in a["sb"].as_array().unwrap().iter().enumerate() {
                let s = ints(s);
                at[i].stereo_bond_neighbor[j] = s[0] as u16;
                at[i].stereo_bond_ord[j] = s[1] as i8;
                at[i].stereo_bond_z_prod[j] = s[2] as i8;
                at[i].stereo_bond_parity[j] = s[3] as i8;
            }
        }
        let carb: Vec<AtStereoCarb> = inp["carb"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| {
                let x = ints(x);
                AtStereoCarb {
                    at_num: x[0] as u16,
                    parity: x[1] as u8,
                }
            })
            .collect();
        let dble: Vec<AtStereoDble> = inp["dble"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| {
                let x = ints(x);
                AtStereoDble {
                    at_num1: x[0] as u16,
                    at_num2: x[1] as u16,
                    parity: x[2] as u8,
                }
            })
            .collect();
        let mut got: Vec<String> = Vec::new();
        for call in inp["calls"].as_array().unwrap() {
            let k = call.as_array().unwrap();
            let i = |j: usize| k[j].as_i64().unwrap();
            let mut v1: Vec<u16> = ints(&k[7]).iter().map(|x| *x as u16).collect();
            let mut v2: Vec<u16> = ints(&k[8]).iter().map(|x| *x as u16).collect();
            v1.resize(nmax, 0);
            v2.resize(nmax, 0);
            let (mut c1, mut c2) = (vec![0u16; nmax], vec![0u16; nmax]);
            let ret = calculated_paths_parities_are_identical(
                &at,
                n,
                &r1,
                &crank,
                &acanon,
                &mut c1,
                &mut c2,
                &mut v1,
                &mut v2,
                i(0) as u16,
                i(1) as u16,
                i(2) as u16,
                i(3) as u16,
                i(4) as i32,
                i(5) as i32,
                i(6) as i32,
                &carb,
                &dble,
                3,
            );
            got.push(format!(
                "Q {ret} |{}{}{}{}",
                pa(&v1, n),
                pa(&v2, n),
                pa(&c1, n),
                pa(&c2, n)
            ));
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn removal_of_calculated_non_stereo_parities_matches_the_c_code() {
    use rdoxide::external::inchi::ichimap::*;
    use rdoxide::external::inchi::ichister::{AtStereoCarb, AtStereoDble, SpAtom};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_mapremove.json")).unwrap();
    fn ints(v: &serde_json::Value) -> Vec<i64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect()
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let atoms = inp["atoms"].as_array().unwrap();
        let n = atoms.len();
        let nmax = n + 2;
        let mut at = vec![SpAtom::default(); n + 1];
        let (mut crank, mut r1, mut acanon) = (
            vec![0u16; nmax + 8],
            vec![0u16; nmax + 8],
            vec![0u16; nmax + 8],
        );
        for (i, a) in atoms.iter().enumerate() {
            let g = |k: &str| a[k].as_i64().unwrap();
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            at[i].parity = g("par") as i8;
            at[i].stereo_atom_parity = g("sap") as i8;
            at[i].num_h = g("nh") as i8;
            at[i].n_ring_system = inp["rs"][i].as_u64().unwrap() as u16;
            r1[i] = g("rank") as u16;
            crank[i] = g("crank") as u16;
            acanon[i] = g("ac") as u16;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
            for (j, s) in a["sb"].as_array().unwrap().iter().enumerate() {
                let s = ints(s);
                at[i].stereo_bond_neighbor[j] = s[0] as u16;
                at[i].stereo_bond_ord[j] = s[1] as i8;
                at[i].stereo_bond_z_prod[j] = s[2] as i8;
                at[i].stereo_bond_parity[j] = s[3] as i8;
            }
        }
        let mut cs = RemoveStereoCs {
            linear_ct_stereo_carb: inp["carb"]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| {
                    let x = ints(x);
                    AtStereoCarb {
                        at_num: x[0] as u16,
                        parity: x[1] as u8,
                    }
                })
                .collect(),
            linear_ct_stereo_dble: inp["dble"]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| {
                    let x = ints(x);
                    AtStereoDble {
                        at_num1: x[0] as u16,
                        at_num2: x[1] as u16,
                        parity: x[2] as u8,
                    }
                })
                .collect(),
            n_symm_rank: r1.clone(),
            l_num_neigh_list_iter: 0,
        };
        let mut nl = create_neigh_list(n, n, &at, false, None);
        let mut s1 = vec![vec![0u16; nmax], vec![0u16; nmax]];
        let mut s2 = vec![vec![0u16; nmax], vec![0u16; nmax]];
        for i in 0..n {
            s1[1][i] = i as u16;
            s2[1][i] = i as u16;
        }
        let mut tmp = vec![0u16; nmax + 8];
        let ret = remove_calculated_non_stereo(
            &mut at, n, n, &mut s1, &mut s2, &mut tmp, &mut nl, &r1, &crank, &acanon, &mut cs, 3,
        );
        let mut got = vec![format!(
            "R {ret} {} {} {}",
            cs.l_num_neigh_list_iter,
            cs.linear_ct_stereo_carb.len(),
            cs.linear_ct_stereo_dble.len()
        )];
        for a in &at[..n] {
            let mut s = format!(
                "a {} {} {} |",
                a.parity, a.stereo_atom_parity, a.final_parity
            );
            for j in 0..3 {
                s += &format!(
                    " {},{},{}",
                    a.stereo_bond_neighbor[j], a.stereo_bond_parity[j], a.stereo_bond_ord[j]
                );
            }
            got.push(s);
        }
        for x in &cs.linear_ct_stereo_carb {
            got.push(format!("c {} {}", x.at_num, x.parity));
        }
        for x in &cs.linear_ct_stereo_dble {
            got.push(format!("d {} {} {}", x.at_num1, x.at_num2, x.parity));
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn stereo_canonical_mapping_matches_the_c_code() {
    use rdoxide::external::inchi::ichimap::*;
    use rdoxide::external::inchi::ichimap4::*;
    use rdoxide::external::inchi::ichister::{AtStereoCarb, AtStereoDble, SpAtom};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_mapstereo4.json")).unwrap();
    fn ints(v: &serde_json::Value) -> Vec<i64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect()
    }
    fn pa(t: &str, a: &[u16], n: usize) -> String {
        let mut s = t.to_string();
        for x in &a[..n] {
            s += &format!(" {x}");
        }
        s
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let atoms = inp["atoms"].as_array().unwrap();
        let n = atoms.len();
        let nmax = n + 2;
        let mut at = vec![SpAtom::default(); n + 1];
        let (mut crank, mut r1, mut acanon, mut an) = (
            vec![0u16; nmax + 8],
            vec![0u16; nmax + 8],
            vec![0u16; nmax + 8],
            vec![0u16; nmax + 8],
        );
        for (i, a) in atoms.iter().enumerate() {
            let g = |k: &str| a[k].as_i64().unwrap();
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            at[i].parity = g("par") as i8;
            at[i].stereo_atom_parity = g("sap") as i8;
            at[i].num_h = g("nh") as i8;
            at[i].n_ring_system = inp["rs"][i].as_u64().unwrap() as u16;
            r1[i] = g("rank") as u16;
            crank[i] = g("crank") as u16;
            acanon[i] = g("ac") as u16;
            an[i] = g("an") as u16;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
            for (j, s) in a["sb"].as_array().unwrap().iter().enumerate() {
                let s = ints(s);
                at[i].stereo_bond_neighbor[j] = s[0] as u16;
                at[i].stereo_bond_ord[j] = s[1] as i8;
                at[i].stereo_bond_z_prod[j] = s[2] as i8;
                at[i].stereo_bond_parity[j] = s[3] as i8;
            }
        }
        let ctl = ints(&inp["ctl"]);
        let mut cs = StereoMapCs {
            rm: RemoveStereoCs {
                linear_ct_stereo_carb: inp["carb"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| {
                        let x = ints(x);
                        AtStereoCarb {
                            at_num: x[0] as u16,
                            parity: x[1] as u8,
                        }
                    })
                    .collect(),
                linear_ct_stereo_dble: inp["dble"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| {
                        let x = ints(x);
                        AtStereoDble {
                            at_num1: x[0] as u16,
                            at_num2: x[1] as u16,
                            parity: x[2] as u8,
                        }
                    })
                    .collect(),
                n_symm_rank: r1.clone(),
                l_num_neigh_list_iter: 0,
            },
            b_first_ct: ctl[0] as i32,
            b_keep_symm_rank: ctl[1] as i32,
            n_mode: ctl[2] as i32,
            b_rank_used_for_stereo: vec![0; nmax + 8],
            b_atom_used_for_stereo: vec![0; nmax + 8],
            n_prev_atom_number: vec![0; nmax + 8],
            ..Default::default()
        };
        let mut nl = create_neigh_list(n, n, &at, false, None);
        let mut s1 = vec![vec![0u16; nmax], vec![0u16; nmax]];
        let mut s2 = vec![vec![0u16; nmax], vec![0u16; nmax]];
        s1[0][..n].copy_from_slice(&r1[..n]);
        s1[1][..n].copy_from_slice(&an[..n]);
        s2[0][..n].copy_from_slice(&r1[..n]);
        s2[1][..n].copy_from_slice(&an[..n]);
        let mut tmp = vec![0u16; nmax + 8];
        let mut canon_to = vec![0u16; nmax + 8];
        let mut symm_stereo: Vec<u16> = (0..nmax as u16 + 8).collect();
        for x in symm_stereo.iter_mut().skip(n) {
            *x = 0;
        }
        let mut ct = CurTree::default();
        ct.alloc(n as i32);
        let ret = map_stereo_bonds4(
            &mut at,
            n,
            n,
            nmax,
            false,
            &crank,
            &acanon,
            &mut canon_to,
            &r1,
            &mut s1,
            &mut s2,
            &mut tmp,
            ctl[3] as i32,
            &mut symm_stereo,
            &mut nl,
            &mut cs,
            &mut ct,
            0,
            3,
        );
        let mut got = vec![format!(
            "R {ret} {} {} {} {} {} {} {} {} {} {}",
            cs.b_first_ct,
            cs.b_stereo_is_better,
            cs.counters.l_num_break_ties,
            cs.counters.l_num_neigh_list_iter,
            cs.l_num_tot_ct,
            cs.l_num_decreased_ct,
            cs.l_num_rejected_ct,
            cs.l_num_equal_ct,
            cs.rm.linear_ct_stereo_carb.len(),
            cs.rm.linear_ct_stereo_dble.len()
        )];
        for x in &cs.rm.linear_ct_stereo_carb {
            got.push(format!("c {} {}", x.at_num, x.parity));
        }
        for x in &cs.rm.linear_ct_stereo_dble {
            got.push(format!("d {} {} {}", x.at_num1, x.at_num2, x.parity));
        }
        for a in &at[..n] {
            let mut s = format!(
                "a {} {} {} |",
                a.parity, a.stereo_atom_parity, a.final_parity
            );
            for j in 0..3 {
                s += &format!(
                    " {},{},{}",
                    a.stereo_bond_neighbor[j], a.stereo_bond_parity[j], a.stereo_bond_ord[j]
                );
            }
            got.push(s);
        }
        let mut u = String::from("U");
        for i in 0..n {
            u += &format!(
                " {},{}",
                cs.b_rank_used_for_stereo[i], cs.b_atom_used_for_stereo[i]
            );
        }
        got.push(u);
        got.push(pa("T", &canon_to, n));
        got.push(pa("S", &symm_stereo, n));
        got.push(pa("P", &cs.n_prev_atom_number, n));
        let mut k = format!("K {} {} {} |", ct.cur_len, ct.max_len, ct.incr_len);
        for i in 0..ct.cur_len.max(0) as usize {
            k += &format!(" {}", ct.tree[i]);
        }
        got.push(k);
        got.push(format!("Z {} {}", s1.len(), s2.len()));
        for v in &s1 {
            got.push(pa(" s1", v, n));
        }
        for v in &s2 {
            got.push(pa(" s2", v, n));
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn canonical_search_helpers_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::node_set_create;
    use rdoxide::external::inchi::ichican2::*;
    use rdoxide::external::inchi::ichimap::create_neigh_list;
    use rdoxide::external::inchi::ichister::SpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_can2a.json")).unwrap();
    fn ints(v: &serde_json::Value) -> Vec<i64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect()
    }
    fn pa(t: &str, a: &[u16], n: usize) -> String {
        let mut s = t.to_string();
        for x in &a[..n] {
            s += &format!(" {x}");
        }
        s
    }
    fn pset(t: &str, s: &rdoxide::external::inchi::ichi_bns::NodeSet, l: usize) -> String {
        let mut o = t.to_string();
        for i in 0..s.len_set {
            o += &format!(" {}", s.bitword[l][i]);
        }
        o
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let atoms = inp["atoms"].as_array().unwrap();
        let n = atoms.len();
        let nmax = n + 2;
        let mut at = vec![SpAtom::default(); n + 1];
        let mut p0 = Partition {
            rank: vec![0; nmax],
            at_number: vec![0; nmax],
        };
        let mut pb = Partition {
            rank: vec![0; nmax],
            at_number: vec![0; nmax],
        };
        for (i, a) in atoms.iter().enumerate() {
            let g = |k: &str| a[k].as_i64().unwrap();
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            p0.rank[i] = g("rank") as u16;
            p0.at_number[i] = g("an") as u16;
            pb.rank[i] = g("rank2") as u16;
            pb.at_number[i] = g("an2") as u16;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
        }
        let marks = ints(&inp["marks"]);
        let useaux = inp["useaux"].as_i64().unwrap() != 0;
        let aux: Vec<u16> = ints(&inp["aux"]).iter().map(|x| *x as u16).collect();
        let v1: Vec<u16> = ints(&inp["v1"]).iter().map(|x| *x as u16).collect();
        let v2: Vec<u16> = ints(&inp["v2"]).iter().map(|x| *x as u16).collect();
        let mut ua = UnorderedPartition {
            equ2: vec![0; nmax],
        };
        let mut ub = UnorderedPartition {
            equ2: vec![0; nmax],
        };
        for (i, x) in ints(&inp["ua"]).iter().enumerate() {
            ua.equ2[i] = *x as u16;
        }
        for (i, x) in ints(&inp["ub"]).iter().enumerate() {
            ub.equ2[i] = *x as u16;
        }
        let joins: Vec<Vec<i64>> = inp["joins"].as_array().unwrap().iter().map(ints).collect();
        let cvs = ints(&inp["cv"]);
        let digraph = inp["digraph"].as_i64().unwrap() != 0;
        let ct = ints(&inp["ct"]);
        let nl = ints(&inp["nl"]);
        let maxrj = inp["maxrj"].as_i64().unwrap() as u16;

        let mut got: Vec<String> = Vec::new();
        got.push(format!(
            "D {} {}",
            i32::from(partition_is_discrete(&p0, n)),
            i32::from(partition_is_discrete(&pb, n))
        ));
        for p in [&p0, &pb] {
            let mut w = vec![Cell::default(); 70];
            for k in 1..66 {
                let sz = partition_get_first_cell(p, &mut w, k, n);
                got.push(format!("C {k} {sz} {} {}", w[k - 1].first, w[k - 1].next));
                if sz == 0 {
                    cell_make_empty(&mut w, k);
                    got.push(format!(
                        "E {} {} {}",
                        w[k - 1].first,
                        w[k - 1].next,
                        w[k - 1].prev
                    ));
                    break;
                }
            }
        }
        {
            let mut mcr = node_set_create(n, 2);
            let mut fix = node_set_create(n, 2);
            partition_get_mcr_and_fix_set(&p0, &mut mcr, &mut fix, n, 1);
            partition_get_mcr_and_fix_set(&pb, &mut mcr, &mut fix, n, 2);
            got.push(pset("M1", &mcr, 0));
            got.push(pset("F1", &fix, 0));
            got.push(pset("M2", &mcr, 1));
            got.push(pset("F2", &fix, 1));
            got.push(format!(
                "AN {} {} {} {}",
                i32::from(all_nodes_are_in_set(&mcr, 1, &fix, 1)),
                i32::from(all_nodes_are_in_set(&fix, 1, &mcr, 1)),
                i32::from(all_nodes_are_in_set(&mcr, 2, &mcr, 1)),
                i32::from(all_nodes_are_in_set(&fix, 2, &mcr, 2))
            ));
            let mut pw = Partition {
                rank: vec![0; nmax],
                at_number: vec![0; nmax],
            };
            partition_copy(&mut pw, &p0, n);
            for i in 0..n {
                if marks[i] != 0 {
                    pw.rank[i] |= RANK_MARK_BIT;
                }
            }
            let auxr = if useaux { Some(&aux[..]) } else { None };
            let mut w = vec![Cell::default(); 70];
            for k in 1..66 {
                let sz = partition_get_first_cell(&p0, &mut w, k, n);
                if sz == 0 {
                    break;
                }
                let nn = cell_get_number_of_nodes(&pw, &w[k - 1]);
                let ki = cell_intersect_with_set(&mut pw, &w[k - 1], &mcr, 1);
                let nn2 = cell_get_number_of_nodes(&pw, &w[k - 1]);
                let mut s = format!("W {k} {nn} {ki} {nn2} |");
                let mut v = 0u16;
                for _ in 0..70 {
                    v = cell_get_min_node(&pw, &w[k - 1], v, auxr);
                    s += &format!(" {v}");
                    if v == INCHI_CANON_INFINITY {
                        break;
                    }
                }
                got.push(s);
            }
            got.push(pa("R", &pw.rank, n));
            let mut wemp = vec![Cell::default(); 1];
            cell_make_empty(&mut wemp, 1);
            got.push(format!(
                "WE {} {}",
                cell_get_min_node(&pw, &wemp[0], 0, auxr),
                cell_get_min_node(&pw, &wemp[0], 3, None)
            ));
        }
        {
            let mut s = node_set_create(n, 2);
            node_set_from_vertices(&mut s, 1, &v1, v1.len());
            node_set_from_vertices(&mut s, 2, &v2, v2.len());
            got.push(pset("S1", &s, 0));
            got.push(pset("S2", &s, 1));
            got.push(format!(
                "SA {} {} {}",
                i32::from(all_nodes_are_in_set(&s, 1, &s, 2)),
                i32::from(all_nodes_are_in_set(&s, 2, &s, 1)),
                i32::from(all_nodes_are_in_set(&s, 1, &s, 1))
            ));
        }
        {
            let mut g = Transposition {
                n_at_numb: vec![0; nmax],
            };
            partition_get_transposition(&p0, &pb, n, &mut g);
            got.push(pa("T", &g.n_at_numb, n));
        }
        {
            let mut a2 = ua.clone();
            let mut b2 = ub.clone();
            let mut s = String::from("G");
            for i in 0..n {
                s += &format!(" {}", n_get_mcr2(&mut a2.equ2, i as u16));
            }
            got.push(s);
            got.push(pa("GA", &a2.equ2, n));
            let mut s = String::from("GN");
            for i in 1..=n {
                s += &format!(" {}", get_unordered_partition_mcr_node(&mut b2, i as u16));
            }
            got.push(s);
            for j in &joins {
                got.push(format!(
                    "J {}",
                    n_join2_mcrs2(&mut b2.equ2, j[0] as u16, j[1] as u16)
                ));
            }
            got.push(pa("JB", &b2.equ2, n));
            got.push(format!("UJ {}", unordered_partition_join(&ua, &mut b2, n)));
            got.push(pa("UB", &b2.equ2, n));
            unordered_partition_make_discrete(&mut a2, n);
            got.push(pa("MD", &a2.equ2, n));
        }
        got.push(format!(
            "L {} {}",
            i32::from(partition_satisfies_lemma_2_25(&p0, n)),
            i32::from(partition_satisfies_lemma_2_25(&pb, n))
        ));
        for &cv in &cvs {
            let mut g = create_neigh_list(n, n, &at, false, None);
            let mut pp = vec![
                Partition::default(),
                Partition::default(),
                Partition::default(),
            ];
            pp[0] = Partition {
                rank: vec![0; nmax],
                at_number: vec![0; nmax],
            };
            partition_copy(&mut pp[0], &p0, n);
            let ret = partition_color_vertex(&mut g, &mut pp, cv as u16, n, n, nmax, digraph, 3);
            let mut s = format!("PV {ret} |");
            for i in 0..n {
                s += &format!(" {}", pp[1].rank[i]);
            }
            s += " |";
            for i in 0..n {
                s += &format!(" {}", pp[1].at_number[i]);
            }
            got.push(s);
        }
        {
            let (ctmax, ctat, hasnumh, hasfix, hasiso, hasexc, clear_k, ctfill) = (
                ct[0] as i32,
                ct[1] as i32,
                ct[2],
                ct[3],
                ct[4],
                ct[5],
                ct[6] as usize,
                ct[7] as i32,
            );
            let mut t = c_table_create(
                n,
                ctmax,
                ctat,
                (hasnumh != 0).then_some(5),
                (hasfix != 0).then_some(6),
                (hasiso != 0).then_some(7),
                (hasexc != 0).then_some(8),
            );
            got.push(format!(
                "CT 1 {} {} {} {} {} {} {} {} {} {} {}",
                t.maxlen_ct,
                t.max_pos,
                t.max_vert,
                t.maxlen_num_h,
                t.maxlen_iso_sort_key,
                t.maxlen_iso_exchg_atnos,
                t.len_ct,
                t.len_pos,
                t.n_len_ct_at_only,
                i32::from(!t.num_h.is_empty()),
                i32::from(!t.num_h_fixed.is_empty())
            ));
            for i in 0..ctmax.min(ctfill) {
                t.ctbl[i as usize] = 100 + i as u16;
            }
            for i in 0..=n {
                t.next_ctbl_pos[i] = if (i as i32) < ctfill {
                    i as u16
                } else {
                    ctfill as u16
                };
            }
            t.len_ct = ctfill.min(ctmax);
            t.len_pos = n as i32;
            ct_part_clear(&mut t, clear_k);
            let mut s = format!("CC {} {} |", t.len_ct, t.len_pos);
            for i in 0..ctmax as usize {
                s += &format!(" {}", t.ctbl[i]);
            }
            got.push(s);
        }
        {
            let rk: Vec<u16> = (0..n)
                .map(|i| p0.rank[i] | if marks[i] != 0 { RANK_MARK_BIT } else { 0 })
                .collect();
            let mut l = vec![0u16; 10];
            l[0] = nl.len() as u16;
            for (i, x) in nl.iter().enumerate() {
                l[i + 1] = *x as u16;
            }
            insertions_sort_neigh_list_at_numbers2(&mut l, &rk, maxrj);
            let mut s = String::from("IS");
            for x in &l[..=nl.len()] {
                s += &format!(" {x}");
            }
            got.push(s);
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn connection_table_parts_match_the_c_code() {
    use rdoxide::external::inchi::ichi_bns::node_set_create;
    use rdoxide::external::inchi::ichican2::*;
    use rdoxide::external::inchi::ichimap::create_neigh_list;
    use rdoxide::external::inchi::ichister::SpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_can2b.json")).unwrap();
    fn ints(v: &serde_json::Value) -> Vec<i64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect()
    }
    fn dump(t: &str, c: &ConTable, n: usize) -> String {
        let mut s = format!(
            "{t} {} {} {} {} {} |",
            c.len_ct, c.len_pos, c.len_num_h, c.len_iso_sort_key, c.len_iso_exchg_atnos
        );
        let mut i = 0;
        while i <= c.len_ct && i < c.maxlen_ct {
            s += &format!(" {}", c.ctbl[i as usize]);
            i += 1;
        }
        s += " |";
        for i in 0..c.len_pos as usize {
            s += &format!(" {}", c.next_ctbl_pos[i]);
        }
        s += " |";
        for i in 0..c.len_pos as usize {
            s += &format!(" {}", c.next_at_rank[i]);
        }
        s += " |";
        if !c.num_h.is_empty() {
            for i in 0..c.len_num_h as usize {
                s += &format!(" {}", c.num_h[i]);
            }
        }
        s += " |";
        if !c.num_h_fixed.is_empty() {
            for i in 0..n {
                s += &format!(" {}", c.num_h_fixed[i]);
            }
        }
        s += " |";
        if !c.iso_sort_key.is_empty() {
            for i in 0..c.len_iso_sort_key as usize {
                s += &format!(" {}", c.iso_sort_key[i]);
            }
        }
        s += " |";
        if !c.iso_exchg_atnos.is_empty() {
            for i in 0..c.len_iso_exchg_atnos as usize {
                s += &format!(" {}", c.iso_exchg_atnos[i]);
            }
        }
        s
    }
    fn layers(kl: &[KLeast]) -> String {
        kl.iter()
            .take(7)
            .map(|e| format!(" {},{}", e.k, e.i))
            .collect()
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let atoms = inp["atoms"].as_array().unwrap();
        let n = atoms.len();
        let mut at = vec![SpAtom::default(); n + 1];
        let mut crank = vec![0usize; n];
        for (i, a) in atoms.iter().enumerate() {
            let nb = a["nb"].as_array().unwrap();
            at[i].valence = nb.len() as i8;
            crank[i] = a["crank"].as_u64().unwrap() as usize;
            for (j, q) in nb.iter().enumerate() {
                at[i].neighbor[j] = q.as_u64().unwrap() as u16;
            }
        }
        let crank_b: Vec<usize> = ints(&inp["crankB"]).iter().map(|x| *x as usize).collect();
        let ms_a: Vec<usize> = ints(&inp["msA"]).iter().map(|x| *x as usize).collect();
        let ms_b: Vec<usize> = ints(&inp["msB"]).iter().map(|x| *x as usize).collect();
        let fl = ints(&inp["flags"]);
        let (b_only, b_split, at_only) = (fl[0] != 0, fl[1] != 0, fl[2] as i32);
        let vals: Vec<Vec<i64>> = inp["vals"].as_array().unwrap().iter().map(ints).collect();
        let perm = ints(&inp["perm"]);
        let cd = CanonData {
            n_canon_flags: 0,
            linear_ct: Vec::new(),
            n_max_len_linear_ct: 8 * n as i32 + 8,
            n_len_linear_ct: fl[3] as i32,
            n_len_ct_at_only: at_only,
            num_h: (fl[4] & 1 != 0).then(|| vals.iter().map(|v| v[0] as i16).collect()),
            maxlen_num_h: n as i32 + 2,
            num_h_fixed: (fl[4] & 2 != 0).then(|| vals.iter().map(|v| v[1] as i16).collect()),
            maxlen_num_h_fixed: n as i32 + 2,
            iso_sort_key: (fl[4] & 4 != 0).then(|| vals.iter().map(|v| v[2]).collect()),
            maxlen_iso_sort_key: n as i32 + 2,
            iso_exchg_atnos: (fl[4] & 8 != 0).then(|| vals.iter().map(|v| v[3] as i8).collect()),
            maxlen_iso_exchg_atnos: n as i32 + 2,
            n_aux_rank: None,
        };
        let create = || {
            c_table_create(
                n,
                cd.n_max_len_linear_ct,
                cd.n_len_ct_at_only,
                cd.num_h.is_some().then_some(cd.maxlen_num_h),
                cd.num_h_fixed.is_some().then_some(cd.maxlen_num_h_fixed),
                cd.iso_sort_key.is_some().then_some(cd.maxlen_iso_sort_key),
                cd.iso_exchg_atnos
                    .is_some()
                    .then_some(cd.maxlen_iso_exchg_atnos),
            )
        };
        let build = |crank: &[usize], ms: &[usize]| -> ConTable {
            let mut ct = create();
            let mut g = create_neigh_list(n, n, &at, false, None);
            let mut p = Partition {
                rank: vec![0; n + 2],
                at_number: vec![0; n + 2],
            };
            for (t, &m) in ms.iter().enumerate() {
                let mut pos = 0;
                for r in 1..=m {
                    for x in 0..n {
                        if crank[x] == r {
                            p.at_number[pos] = x as u16;
                            pos += 1;
                        }
                    }
                }
                for x in 0..n {
                    if crank[x] > m {
                        p.at_number[pos] = x as u16;
                        pos += 1;
                    }
                }
                for x in 0..n {
                    p.rank[x] = if crank[x] <= m {
                        crank[x] as u16
                    } else {
                        n as u16
                    };
                }
                ct_part_fill(&mut g, &cd, &p, &mut ct, t + 1, n, n, n + 2);
            }
            ct
        };
        let a = build(&crank, &ms_a);
        let b = build(&crank_b, &ms_b);
        let mut got = vec![dump("A", &a, n), dump("B", &b, n)];
        {
            let kmax = a.len_pos.min(b.len_pos) as usize;
            let mut cmp = vec![0i8; 128];
            let mut kl = vec![KLeast::default(); MAX_LAYERS];
            for k in 1..=kmax {
                let d1 = ct_part_compare(&a, &b, Some(&mut cmp), Some(&mut kl), k, b_only, b_split);
                got.push(format!("P1 {k} {d1} {} |{}", cmp[k - 1], layers(&kl)));
                let d2 = ct_part_compare(&a, &b, None, None, k, b_only, b_split);
                let d3 = ct_part_compare(&b, &a, None, None, k, !b_only, b_split);
                got.push(format!("P2 {d2} {d3}"));
            }
            got.push(format!(
                "FC {} {} {} {}",
                ct_full_compare(&a, &b, b_only, b_split),
                ct_full_compare(&a, &b, !b_only, !b_split),
                ct_full_compare(&b, &a, false, false),
                ct_full_compare(&a, &a, false, false)
            ));
            got.push(format!("FL {}", ct_full_compare_layers(&kl)));
            let (mut l, mut i, mut k) = (0, 0, 0);
            let g1 = ct_compare_layers_get_first_diff(Some(&kl), 2, &mut l, &mut i, &mut k);
            got.push(format!("GD {g1} {l} {i} {k}"));
            let g1 = ct_compare_layers_get_first_diff(None, 2, &mut l, &mut i, &mut k);
            got.push(format!("GN {g1}"));
            got.push(format!(
                "PL {} {} {}",
                ct_part_compare_layers(Some(&kl), 5, 0),
                ct_part_compare_layers(Some(&kl), 0, 3),
                ct_part_compare_layers(None, 1, 1)
            ));
            update_compare_layers(Some(&mut kl), 2);
            got.push(format!("UL{}", layers(&kl)));
        }
        {
            let mut cc = create();
            let mut dd = create();
            for k in 1..=a.len_pos as usize {
                ct_part_copy(&mut cc, &a, k);
            }
            got.push(dump("C", &cc, n));
            ct_full_copy(&mut dd, &b);
            got.push(dump("D", &dd, n));
            let mut cmp = vec![1i8; 128];
            ct_part_inchi_canon_infinity(&mut cc, Some(&mut cmp), a.len_pos as usize);
            let mut s = String::from("INF");
            for i in 0..=cc.len_ct as usize {
                s += &format!(" {}", cc.ctbl[i]);
            }
            s += " |";
            for x in &cmp[..8] {
                s += &format!(" {x}");
            }
            got.push(s);
            ct_part_inchi_canon_infinity(&mut dd, None, 1);
            got.push(format!("INF1 {}", dd.ctbl[0]));
            ct_part_clear(&mut cc, 1);
            got.push(format!(
                "CLR {}",
                get_one_additional_layer(Some(&cd), Some(&cc))
            ));
            got.push(format!(
                "GOL {} {}",
                get_one_additional_layer(None, Some(&a)),
                get_one_additional_layer(Some(&cd), None)
            ));
        }
        {
            let mut g = Transposition {
                n_at_numb: vec![0; n + 2],
            };
            for i in 0..n {
                g.n_at_numb[i] = perm[i] as u16;
            }
            let mut m = node_set_create(n, 2);
            let mut f = node_set_create(n, 2);
            let mut up = UnorderedPartition {
                equ2: vec![0; n + 2],
            };
            transposition_get_mcr_and_fix_set_and_unordered_partition(
                &mut g, &mut m, &mut f, n, 2, &mut up,
            );
            let mut s = String::from("TG");
            for i in 0..n {
                s += &format!(" {}", g.n_at_numb[i]);
            }
            s += " | M";
            for x in &m.bitword[1] {
                s += &format!(" {x}");
            }
            s += " | F";
            for x in &f.bitword[1] {
                s += &format!(" {x}");
            }
            s += " | U";
            for i in 0..n {
                s += &format!(" {}", up.equ2[i]);
            }
            got.push(s);
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn canonical_ranking_setup_matches_the_c_code() {
    use rdoxide::external::inchi::ichican2::*;
    use rdoxide::external::inchi::ichimap::AT_INV_LENGTH;
    use rdoxide::external::inchi::ichister::SpAtom;
    use rdoxide::external::inchi::ichitaut::{TGroup, TGroupInfo};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_can2c.json")).unwrap();
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let n = inp["n"].as_u64().unwrap() as usize;
        let nat = inp["nat"].as_u64().unwrap() as usize;
        let ntg = inp["ntg"].as_u64().unwrap() as usize;
        let f: Vec<bool> = inp["flags"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() != 0)
            .collect();
        let mut at = vec![SpAtom::default(); n + 1];
        for (i, a) in inp["atoms"].as_array().unwrap().iter().enumerate() {
            let el = a[0].as_str().unwrap();
            if el != "-" {
                at[i].elname[..el.len()].copy_from_slice(el.as_bytes());
            }
            at[i].valence = a[1].as_i64().unwrap() as i8;
            at[i].endpoint = a[2].as_i64().unwrap() as u16;
            at[i].num_h = a[3].as_i64().unwrap() as i8;
            at[i].iso_sort_key = a[4].as_i64().unwrap();
        }
        let mut tgs = vec![TGroup::default(); ntg + 8];
        for (k, t) in inp["tg"].as_array().unwrap().iter().enumerate() {
            tgs[k].n_num_endpoints = t[0].as_u64().unwrap() as u16;
            for j in 0..5 {
                tgs[k].num[j] = t[1 + j].as_u64().unwrap() as u16;
            }
        }
        let tgi = TGroupInfo {
            t_group: (inp["hasArr"].as_i64().unwrap() != 0).then_some(tgs),
            num_t_groups: ntg,
            b_ignore_isotopic: inp["ign"].as_i64().unwrap() as _,
            ..Default::default()
        };
        let mut inv = vec![Default::default(); nat + 2];
        fill_out_atom_invariant2(
            &at,
            n,
            nat,
            &mut inv,
            f[0],
            f[1],
            f[2],
            f[3],
            f[4],
            f[5].then_some(&tgi),
        );
        let mut got: Vec<String> = Vec::new();
        for v in &inv[..nat] {
            let mut s = String::from("INV");
            for j in 0..AT_INV_LENGTH {
                s += &format!(" {}", v.val[j]);
            }
            s += &format!(" {}", v.iso_sort_key);
            got.push(s);
        }
        let mut rk = vec![0u16; nat + 2];
        let mut an = vec![0u16; nat + 2];
        let nd = set_initial_ranks2(nat, &inv, &mut rk, &mut an);
        let join = |v: &[u16]| v.iter().map(|x| format!(" {x}")).collect::<String>();
        got.push(format!(
            "SR {nd} |{} |{}",
            join(&rk[..nat]),
            join(&an[..nat])
        ));
        let mut nh: Vec<i16> = inp["nh"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i16)
            .collect();
        let l = nh.len();
        clean_num_h(&mut nh, l);
        got.push(format!(
            "CNH{}",
            nh.iter().map(|x| format!(" {x}")).collect::<String>()
        ));
        let mut ct: Vec<u16> = inp["ct"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_u64().unwrap() as u16)
            .collect();
        let r = clean_ct(&mut ct, inp["ctlen"].as_u64().unwrap() as usize);
        got.push(format!("CCT {r}{}", join(&ct)));
        got.push("CCN 0".to_string());
        let mut k: Vec<i64> = inp["isk"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect();
        let l = k.len();
        clean_iso_sort_keys(&mut k, l);
        got.push(format!(
            "CIK{}",
            k.iter().map(|x| format!(" {x}")).collect::<String>()
        ));
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn canonical_graph_search_matches_the_c_code() {
    use rdoxide::external::inchi::ichican2::*;
    use rdoxide::external::inchi::ichimap::{create_neigh_list, differentiate_ranks_basic};
    use rdoxide::external::inchi::ichister::SpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_can2d.json")).unwrap();
    fn ints(v: &serde_json::Value) -> Vec<i64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect()
    }
    fn join<T: std::fmt::Display>(t: &str, v: &[T]) -> String {
        let mut s = t.to_string();
        for x in v {
            s += &format!(" {x}");
        }
        s
    }
    fn dump(t: &str, c: &ConTable, n: usize) -> String {
        let mut s = format!(
            "{t} {} {} {} {} {} |",
            c.len_ct, c.len_pos, c.len_num_h, c.len_iso_sort_key, c.len_iso_exchg_atnos
        );
        let mut i = 0;
        while i <= c.len_ct && i < c.maxlen_ct {
            s += &format!(" {}", c.ctbl[i as usize]);
            i += 1;
        }
        s += " |";
        for i in 0..c.len_pos as usize {
            s += &format!(" {}", c.next_ctbl_pos[i]);
        }
        s += " |";
        for i in 0..c.len_pos as usize {
            s += &format!(" {}", c.next_at_rank[i]);
        }
        s += " |";
        if !c.num_h.is_empty() {
            for i in 0..c.len_num_h as usize {
                s += &format!(" {}", c.num_h[i]);
            }
        }
        s += " |";
        if !c.num_h_fixed.is_empty() {
            for i in 0..n {
                s += &format!(" {}", c.num_h_fixed[i]);
            }
        }
        s += " |";
        if !c.iso_sort_key.is_empty() {
            for i in 0..c.len_iso_sort_key as usize {
                s += &format!(" {}", c.iso_sort_key[i]);
            }
        }
        s += " |";
        if !c.iso_exchg_atnos.is_empty() {
            for i in 0..c.len_iso_exchg_atnos as usize {
                s += &format!(" {}", c.iso_exchg_atnos[i]);
            }
        }
        s
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let n = inp["n"].as_u64().unwrap() as usize;
        let dig = inp["dig"].as_i64().unwrap() != 0;
        let large = inp["large"].as_i64().unwrap() != 0;
        let use_aux = inp["aux"].as_i64().unwrap() != 0;
        let masks = ints(&inp["masks"]);
        let col = ints(&inp["col"]);
        let edges: Vec<(usize, usize)> = inp["edges"]
            .as_array()
            .unwrap()
            .iter()
            .map(|e| {
                (
                    e[0].as_u64().unwrap() as usize,
                    e[1].as_u64().unwrap() as usize,
                )
            })
            .collect();
        let m = edges.len();
        let mut at = vec![SpAtom::default(); n + 1];
        for &(a, b) in &edges {
            let va = at[a].valence as usize;
            at[a].neighbor[va] = b as u16;
            at[a].valence += 1;
            let vb = at[b].valence as usize;
            at[b].neighbor[vb] = a as u16;
            at[b].valence += 1;
        }
        let nh: Vec<i16> = ints(&inp["nh"]).iter().map(|&x| x as i16).collect();
        let nhf: Vec<i16> = ints(&inp["nhf"]).iter().map(|&x| x as i16).collect();
        let isk = ints(&inp["isk"]);
        let aux: Vec<u16> = ints(&inp["auxr"]).iter().map(|&x| x as u16).collect();
        let mut got: Vec<String> = Vec::new();
        let mut fix: Option<ConTable> = None;
        for (st, &mask) in masks.iter().enumerate() {
            let mut g = create_neigh_list(n, n, &at, false, None);
            let mut pi: Vec<Partition> = (0..n + 4)
                .map(|_| Partition {
                    rank: vec![0; n + 2],
                    at_number: vec![0; n + 2],
                })
                .collect();
            let mut prev = vec![0u16; n + 4];
            for i in 0..n {
                pi[0].rank[i] = col.iter().filter(|&&cj| cj <= col[i]).count() as u16;
                pi[0].at_number[i] = i as u16;
            }
            let mut seen = std::collections::HashSet::new();
            let nd = (0..n).filter(|&i| seen.insert(pi[0].rank[i])).count() as i32;
            let mut it = 0i64;
            {
                let p0 = &mut pi[0];
                differentiate_ranks_basic(
                    n,
                    &mut g,
                    nd,
                    &mut p0.rank,
                    &mut prev,
                    &mut p0.at_number,
                    &mut it,
                    1,
                );
            }
            let mut cd = CanonData {
                n_max_len_linear_ct: 4 * (n + m) as i32 + 8,
                n_len_ct_at_only: (n + m) as i32,
                ..Default::default()
            };
            cd.linear_ct = vec![0; cd.n_max_len_linear_ct as usize + 2];
            if mask & 1 != 0 {
                cd.num_h = Some(nh.clone());
                cd.maxlen_num_h = n as i32 + 2;
            }
            if mask & 2 != 0 {
                cd.num_h_fixed = Some(nhf.clone());
                cd.maxlen_num_h_fixed = n as i32 + 2;
            }
            if mask & 4 != 0 {
                cd.iso_sort_key = Some(isk.clone());
                cd.maxlen_iso_sort_key = n as i32 + 2;
            }
            if use_aux {
                cd.n_aux_rank = Some(aux.clone());
            }
            let mut sym = vec![0u16; 70];
            let mut can = vec![0u16; 70];
            let mut atn = vec![0u16; 70];
            let mut cc = CanonCounts::default();
            let mut out: Option<ConTable> = None;
            let ret = canon_graph(
                &mut || false,
                n,
                n,
                n + 2,
                dig,
                &mut g,
                &mut pi,
                &mut sym,
                &mut can,
                &mut atn,
                &mut cd,
                &mut cc,
                fix.as_ref(),
                Some(&mut out),
                large,
            );
            got.push(format!(
                "S {st} {ret} {} {} {} {} {} {} {:.0} {} {}",
                cd.n_len_linear_ct,
                cc.l_num_break_ties,
                cc.l_num_decreased_ct,
                cc.l_num_rejected_ct,
                cc.l_num_equal_ct,
                cc.l_num_tot_ct,
                cc.d_group_size,
                cc.l_num_generators,
                cc.l_num_stored_isomorphisms
            ));
            if ret >= 0 {
                got.push(join("sym", &sym[..n]));
                got.push(join("can", &can[..n]));
                got.push(join("atn", &atn[..n]));
                got.push(join("ct", &cd.linear_ct[..cd.n_len_linear_ct as usize]));
                match &out {
                    Some(o) => got.push(dump("out", o, n)),
                    None => got.push("out none".to_string()),
                }
            }
            let mut s = String::from("nl");
            for l in g.iter().take(n) {
                s += " |";
                for j in 1..=l[0] as usize {
                    s += &format!(" {}", l[j]);
                }
            }
            got.push(s);
            fix = if ret >= 0 { out } else { None };
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn base_canonical_numbering_helpers_match_the_c_code() {
    use rdoxide::external::inchi::ichican2::{fix_canon_equivalence_info, make_iso_sort_key};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_can2e.json")).unwrap();
    fn u16s(v: &serde_json::Value) -> Vec<u16> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_u64().unwrap() as u16)
            .collect()
    }
    fn join(v: &[u16]) -> String {
        v.iter().map(|x| format!(" {x}")).collect()
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let n = inp["n"].as_u64().unwrap() as usize;
        let mut sym = u16s(&inp["sym"]);
        let mut cur = u16s(&inp["cur"]);
        let mut an = u16s(&inp["an"]);
        let mut tmp = vec![0u16; 64];
        let (r, bch) = fix_canon_equivalence_info(n, &mut sym, &mut cur, &mut tmp, &mut an);
        let d: Vec<i32> = inp["d"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap() as i32)
            .collect();
        let got = vec![
            format!(
                "F {r} {bch} |{} |{} |{} |{}",
                join(&sym),
                join(&cur),
                join(&tmp[..n]),
                join(&an)
            ),
            format!("K {}", make_iso_sort_key(d[0], d[1], d[2], d[3])),
        ];
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn base_canonical_ranking_matches_the_c_code() {
    use rdoxide::external::inchi::ichican2::*;
    use rdoxide::external::inchi::ichister::SpAtom;
    use rdoxide::external::inchi::ichitaut::{TGroup, TGroupInfo};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_can2f.json")).unwrap();
    fn ints(v: &serde_json::Value) -> Vec<i64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect()
    }
    fn line<T: std::fmt::Display>(t: &str, v: &[T], n: usize) -> String {
        let mut s = t.to_string();
        if v.is_empty() {
            return s + " -";
        }
        for x in v.iter().take(n) {
            s += &format!(" {x}");
        }
        s
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let n = inp["n"].as_u64().unwrap() as usize;
        let ntg = inp["ntg"].as_u64().unwrap() as usize;
        let have_yes = inp["haveYes"].as_i64().unwrap() != 0;
        let have_non = inp["haveNon"].as_i64().unwrap() != 0;
        let have_tgi = inp["haveTgi"].as_i64().unwrap() != 0;
        let edges: Vec<(usize, usize)> = inp["edges"]
            .as_array()
            .unwrap()
            .iter()
            .map(|e| {
                (
                    e[0].as_u64().unwrap() as usize,
                    e[1].as_u64().unwrap() as usize,
                )
            })
            .collect();
        let build = |key: &str| {
            let mut at = vec![SpAtom::default(); n + 1];
            for (i, a) in inp[key].as_array().unwrap().iter().enumerate() {
                let a = ints(a);
                let names = ["C", "N", "O", "S", "Cl", "H"];
                let nm = names[a[0] as usize].as_bytes();
                at[i].elname[..nm.len()].copy_from_slice(nm);
                at[i].num_h = a[1] as i8;
                at[i].iso_atw_diff = a[2] as i8;
                at[i].num_iso_h = [a[3] as i8, a[4] as i8, a[5] as i8];
                at[i].c_flags = a[6] as i8;
                at[i].endpoint = a[7] as u16;
            }
            for &(a, b) in &edges {
                let va = at[a].valence as usize;
                at[a].neighbor[va] = b as u16;
                at[a].valence += 1;
                let vb = at[b].valence as usize;
                at[b].neighbor[vb] = a as u16;
                at[b].valence += 1;
            }
            at
        };
        let at_y = build("atY");
        let at_n = build("atN");
        let sz = ints(&inp["sizes"]);
        let mut s = [AtomSizes::default(); 2];
        for k in 0..2 {
            let v = &sz[7 * k..7 * k + 7];
            s[k] = AtomSizes {
                n_len_ct: v[0] as i32,
                n_len_ct_at_only: v[1] as i32,
                num_isotopic_atoms: v[2] as i32,
                b_ignore_isotopic: v[3] as i32,
                n_len_linear_ct_tautomer: v[4] as i32,
                b_has_isotopic_taut_groups: v[5] as i32,
                n_len_isotopic_endpoints: v[6] as i32,
                ..Default::default()
            };
        }
        let iso = ints(&inp["iso"]);
        let iso_list = ints(&inp["isolist"]);
        let mut tg = vec![TGroup::default(); ntg + 2];
        for (k, t) in inp["tgs"].as_array().unwrap().iter().enumerate() {
            let t = ints(t);
            tg[k].n_num_endpoints = t[0] as u16;
            tg[k].n_first_endpoint_at_no_pos = t[1] as u16;
            tg[k].n_group_number = k as u16 + 1;
            for j in 0..5 {
                tg[k].num[j] = t[2 + j] as u16;
            }
            tg[k].i_weight = (t[2] * 7 + t[4] * 3 + t[5]) % 5;
        }
        let mut tgi = TGroupInfo {
            t_group: Some(tg),
            n_endpoint_atom_number: ints(&inp["epl"]).iter().map(|&x| x as u16).collect(),
            num_t_groups: ntg,
            max_num_t_groups: ntg,
            n_num_endpoints: inp["nE"].as_i64().unwrap() as i32,
            n_num_isotopic_endpoints: iso[1] as i32,
            b_taut_flags_done: iso[2] as u64,
            b_ignore_isotopic: iso[4] as i32,
            ..Default::default()
        };
        tgi.tni.b_normalization_flags = iso[3] as u64;
        tgi.n_isotopic_endpoint_atom_number = std::iter::once(iso[0] as u16)
            .chain(iso_list.iter().map(|&x| x as u16))
            .collect();
        let mut bcn = Bcn::default();
        let at_refs: [Option<&[SpAtom]>; 2] =
            [have_non.then_some(&at_n[..]), have_yes.then_some(&at_y[..])];
        let ret = get_base_canon_ranking(
            &mut || false,
            n,
            n + ntg,
            at_refs,
            have_tgi.then_some(&mut tgi),
            &s,
            &mut bcn,
            inp["bfix"].as_i64().unwrap() != 0,
            false,
        );
        let mut got: Vec<String> = Vec::new();
        got.push(format!(
            "R {ret} {} {} {} {}",
            bcn.num_max, bcn.num_at_tg, bcn.num_atoms, bcn.n_max_len_rank_stack
        ));
        got.push(format!("TI {}", tgi.b_ignore_isotopic));
        for j in 0..(bcn.n_max_len_rank_stack / 2) as usize {
            let p = bcn.p_rank_stack.get(j).cloned().unwrap_or_default();
            got.push(line(&format!("RS{j}"), &p.rank, bcn.num_max));
            got.push(line(&format!("RT{j}"), &p.at_number, bcn.num_max));
        }
        if ret == 0 {
            for k in 0..2 {
                let f = &bcn.ftcn[k];
                got.push(format!(
                    "F{k} {} {} {} {} {} {} {} {} {}",
                    f.num_at_tg,
                    f.num_atoms,
                    f.n_canon_flags,
                    f.n_len_linear_ct_at_only,
                    f.n_len_linear_ct,
                    f.n_max_len_linear_ct,
                    f.n_len_num_h,
                    f.n_len_num_h_fix_h,
                    f.len_iso_sort_keys
                ));
                let mut s = String::from("NL");
                match &f.neigh_list {
                    Some(nl) => {
                        for l in nl.iter().take(f.num_at_tg) {
                            s += " |";
                            for j in 1..=l[0] as usize {
                                s += &format!(" {}", l[j]);
                            }
                        }
                    }
                    None => s += " -",
                }
                got.push(s);
                let lc_n = if f.linear_ct.is_empty() {
                    0
                } else {
                    f.n_len_linear_ct as usize + 1
                };
                got.push(line("LC", &f.linear_ct, lc_n));
                got.push(line("PR", &f.partition_ct.rank, bcn.num_max));
                got.push(line("PA", &f.partition_ct.at_number, bcn.num_max));
                got.push(line("SY", &f.n_symm_rank_ct, bcn.num_max));
                let l = f.n_len_num_h as usize;
                got.push(line("NHO", &f.n_num_h_orig, l));
                got.push(line("NH", &f.n_num_h, l));
                got.push(line("NFO", &f.n_num_h_orig_fix_h, l));
                got.push(line("NF", &f.n_num_h_fix_h, l));
                got.push(line("IR", &f.partition_ct_iso.rank, bcn.num_max));
                got.push(line("IA", &f.partition_ct_iso.at_number, bcn.num_max));
                got.push(line("IS", &f.n_symm_rank_ct_iso, bcn.num_max));
                let l = f.len_iso_sort_keys as usize;
                got.push(line("K", &f.iso_sort_keys, l));
                got.push(line("KO", &f.iso_sort_keys_orig, l));
                got.push(line("X", &f.iso_exchg_atnos, l));
                got.push(line("XO", &f.iso_exchg_atnos_orig, l));
            }
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn canonicalization_linear_ct_helpers_match_the_c_code() {
    use rdoxide::external::inchi::ichican2::AtomSizes;
    use rdoxide::external::inchi::ichicano::*;
    use rdoxide::external::inchi::ichister::SpAtom;
    use rdoxide::external::inchi::ichitaut::{TGroup, TGroupInfo};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_cano1.json")).unwrap();
    fn ints(v: &serde_json::Value) -> Vec<i64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect()
    }
    fn join<T: std::fmt::Display>(v: &[T]) -> String {
        v.iter().map(|x| format!(" {x}")).collect()
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let n = inp["n"].as_u64().unwrap() as usize;
        let ntg = inp["ntg"].as_u64().unwrap() as usize;
        let n_e = inp["nE"].as_u64().unwrap() as usize;
        let ntot = n + ntg;
        let mut at = vec![SpAtom::default(); n + 1];
        for (i, a) in inp["atoms"].as_array().unwrap().iter().enumerate() {
            let a = ints(a);
            at[i].iso_sort_key = a[0];
            at[i].parity = a[1] as i8;
            for j in 0..3 {
                at[i].stereo_bond_neighbor[j] = a[2 + j] as u16;
            }
            at[i].endpoint = a[5] as u16;
            at[i].num_iso_h = [a[6] as i8, a[7] as i8, a[8] as i8];
            at[i].iso_atw_diff = a[9] as i8;
            at[i].c_flags = a[10] as i8;
        }
        for e in inp["edges"].as_array().unwrap() {
            let (a, b) = (
                e[0].as_u64().unwrap() as usize,
                e[1].as_u64().unwrap() as usize,
            );
            let va = at[a].valence as usize;
            at[a].neighbor[va] = b as u16;
            at[a].valence += 1;
            let vb = at[b].valence as usize;
            at[b].neighbor[vb] = a as u16;
            at[b].valence += 1;
        }
        let mut tg = vec![TGroup::default(); ntg + 2];
        for (k, t) in inp["tgs"].as_array().unwrap().iter().enumerate() {
            let t = ints(t);
            tg[k].n_num_endpoints = t[0] as u16;
            tg[k].n_first_endpoint_at_no_pos = t[1] as u16;
            for j in 0..5 {
                tg[k].num[j] = t[2 + j] as u16;
            }
            tg[k].i_weight = t[7];
        }
        let sv = ints(&inp["sv"]);
        let mut tgi = TGroupInfo {
            t_group: (ntg > 0).then_some(tg),
            n_endpoint_atom_number: ints(&inp["epl"]).iter().map(|&x| x as u16).collect(),
            num_t_groups: ntg,
            n_num_endpoints: n_e as i32,
            n_num_isotopic_endpoints: inp["nIsoEnd"].as_i64().unwrap() as i32,
            t_group_number: vec![0; 4 * ntg + 4],
            ..Default::default()
        };
        let has_tgi = inp["hasTgi"].as_i64().unwrap() != 0;
        let mut got: Vec<String> = Vec::new();
        {
            let mut s = AtomSizes {
                n_len_bonds: sv[0] as i32,
                n_len_ct_at_only: sv[1] as i32,
                n_len_ct: sv[2] as i32,
                n_len_isotopic: sv[3] as i32,
                n_len_linear_ct_stereo_dble: sv[4] as i32,
                n_len_linear_ct_stereo_carb: sv[5] as i32,
                n_len_isotopic_endpoints: sv[6] as i32,
                n_len_linear_ct_tautomer: sv[7] as i32,
                ..Default::default()
            };
            let r = get_canon_lengths(n, &at, &mut s, has_tgi.then_some(&tgi));
            got.push(format!(
                "L {r} {} {} {} {} {} {} {}",
                s.n_len_bonds,
                s.n_len_ct_at_only,
                s.n_len_ct,
                s.n_len_isotopic,
                s.n_len_linear_ct_stereo_dble,
                s.n_len_linear_ct_stereo_carb,
                s.n_len_isotopic_endpoints
            ));
        }
        let mut anum: Vec<Vec<u16>> = Vec::new();
        let mut symm: Vec<Vec<u16>> = Vec::new();
        let mut rnk: Vec<Vec<u16>> = Vec::new();
        for r in inp["ranks"].as_array().unwrap() {
            let a: Vec<u16> = ints(&r[0]).iter().map(|&x| x as u16).collect();
            let sy: Vec<u16> = ints(&r[1]).iter().map(|&x| x as u16).collect();
            let mut rk = vec![0u16; ntot + 2];
            for (i, &x) in a.iter().enumerate() {
                rk[x as usize] = i as u16 + 1;
            }
            anum.push(a);
            symm.push(sy);
            rnk.push(rk);
        }
        let l = ints(&inp["lens"]);
        let (nmax, plen, nmax_iso_at, plen_iso_at, nmax_t, plen_t, nmax_it, plen_it, b_iso) = (
            l[0] as i32,
            l[1] as i32,
            l[2] as i32,
            l[3] as i32,
            l[4] as i32,
            l[5] as i32,
            l[6] as i32,
            l[7] as i32,
            l[8] != 0,
        );
        {
            let mut iso = vec![AtIsotopic::default(); nmax_iso_at as usize + 2];
            let mut len = plen_iso_at;
            let r = fill_isotopic_at_linear_ct(
                n,
                &at,
                &anum[0],
                if nmax_iso_at != 0 { &mut iso } else { &mut [] },
                nmax_iso_at,
                &mut len,
            );
            let mut s = format!("I {r} {len} |");
            for e in iso.iter().take(nmax_iso_at as usize) {
                s += &format!(
                    " {},{},{},{},{}",
                    e.at_num, e.num_1h, e.num_d, e.num_t, e.iso_atw_diff
                );
            }
            got.push(s);
        }
        if has_tgi {
            let mut lt = vec![0u16; nmax_t as usize + 2];
            let mut lit = vec![AtIsoTgroup::default(); nmax_it as usize + 2];
            let mut len_t = plen_t;
            let mut len_it = plen_it;
            let r = fill_taut_linear_ct2(
                n,
                ntot,
                b_iso,
                &rnk[0],
                &anum[0],
                &symm[0],
                &anum[1],
                &symm[1],
                &mut lt,
                nmax_t,
                &mut len_t,
                &mut lit,
                nmax_it,
                &mut len_it,
                Some(&mut tgi),
            );
            let mut s = format!("T {r} {len_t} {len_it} |{} |", join(&lt[..nmax_t as usize]));
            for e in lit.iter().take(len_it.max(0).min(nmax_it) as usize) {
                s += &format!(" {}:{},{},{}", e.tgroup_num, e.num[0], e.num[1], e.num[2]);
            }
            s += &format!(" |{} |", join(&tgi.t_group_number[..4 * ntg]));
            s += &join(&tgi.n_endpoint_atom_number[..n_e]);
            got.push(s);
        }
        {
            let mut cs = CanonStat {
                n_max_len_linear_ct: nmax,
                linear_ct: vec![0; nmax as usize + 2],
                t_group_info: has_tgi.then(|| tgi.clone()),
                n_len_linear_ct: plen,
                ..Default::default()
            };
            let nn = if has_tgi { ntot } else { n };
            for v in 0..3 {
                let r = update_full_linear_ct(n, nn, &at, &rnk[v], &anum[v], &mut cs, v == 0);
                got.push(format!(
                    "U{v} {r} {} {} |{}",
                    cs.n_len_linear_ct,
                    cs.n_len_linear_ct_at_only,
                    join(&cs.linear_ct[..nmax as usize])
                ));
            }
            let mut s = String::from("UE");
            if let Some(t) = &cs.t_group_info {
                s += &join(&t.n_endpoint_atom_number[..n_e]);
            } else {
                s += &join(&tgi.n_endpoint_atom_number[..n_e]);
            }
            got.push(s);
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn final_canonical_numbering_matches_the_c_code() {
    use rdoxide::external::inchi::ichican2::*;
    use rdoxide::external::inchi::ichicano::*;
    use rdoxide::external::inchi::ichister::{SpAtom, switch_atom_stereo_and_isotopic_stereo};
    use rdoxide::external::inchi::ichitaut::{TGroup, TGroupInfo};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_cano2.json")).unwrap();
    fn ints(v: &serde_json::Value) -> Vec<i64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_i64().unwrap())
            .collect()
    }
    fn pu(t: &str, v: &[u16], n: usize) -> String {
        let mut s = t.to_string();
        if v.is_empty() {
            return s + " -";
        }
        for x in v.iter().take(n) {
            s += &format!(" {x}");
        }
        s
    }
    fn slot_ptr(cs: &CanonStat, dble: bool, sl: Slot) -> usize {
        let empty = if dble {
            cs.dble_ref(sl).is_empty()
        } else {
            cs.carb_ref(sl).is_empty()
        };
        if empty { 0 } else { sl as usize + 1 }
    }
    fn slot_name(cs: &CanonStat, dble: bool, sl: Slot) -> char {
        let p = slot_ptr(cs, dble, sl);
        let pb = slot_ptr(cs, dble, Slot::Base);
        let pi = slot_ptr(cs, dble, Slot::Inv);
        let po = slot_ptr(cs, dble, Slot::Iso);
        let pp = slot_ptr(cs, dble, Slot::IsoInv);
        if p == pb {
            'B'
        } else if p == pi {
            'I'
        } else if p == po {
            'O'
        } else if p == pp {
            'P'
        } else {
            '-'
        }
    }
    for (case, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = &c["in"];
        let n = inp["n"].as_u64().unwrap() as usize;
        let ntg = inp["ntg"].as_u64().unwrap() as usize;
        let n_e = inp["nE"].as_u64().unwrap() as usize;
        let have_yes = inp["haveYes"].as_i64().unwrap() != 0;
        let have_non = inp["haveNon"].as_i64().unwrap() != 0;
        let have_tgi = inp["haveTgi"].as_i64().unwrap() != 0;
        let edges: Vec<(usize, usize)> = inp["edges"]
            .as_array()
            .unwrap()
            .iter()
            .map(|e| {
                (
                    e[0].as_u64().unwrap() as usize,
                    e[1].as_u64().unwrap() as usize,
                )
            })
            .collect();
        let st1: Vec<Vec<i64>> = inp["st1"].as_array().unwrap().iter().map(ints).collect();
        let st2: Vec<Vec<i64>> = inp["st2"].as_array().unwrap().iter().map(ints).collect();
        let build = |key: &str| {
            let mut at = vec![SpAtom::default(); n + 1];
            for (i, a) in inp[key].as_array().unwrap().iter().enumerate() {
                let a = ints(a);
                let names = ["C", "N", "O", "S", "Cl", "H"];
                let nm = names[a[0] as usize].as_bytes();
                at[i].elname[..nm.len()].copy_from_slice(nm);
                at[i].num_h = a[1] as i8;
                at[i].iso_atw_diff = a[2] as i8;
                at[i].num_iso_h = [a[3] as i8, a[4] as i8, a[5] as i8];
                at[i].c_flags = a[6] as i8;
                at[i].endpoint = a[7] as u16;
                at[i].iso_sort_key =
                    make_iso_sort_key(a[2] as i32, a[3] as i32, a[4] as i32, a[5] as i32);
            }
            for &(a, b) in &edges {
                let va = at[a].valence as usize;
                at[a].neighbor[va] = b as u16;
                at[a].valence += 1;
                let vb = at[b].valence as usize;
                at[b].neighbor[vb] = a as u16;
                at[b].valence += 1;
            }
            for i in 0..n {
                at[i].parity = st1[i][0] as i8;
                at[i].stereo_atom_parity = st1[i][2] as i8;
                at[i].n_ring_system = st1[i][4] as u16;
                at[i].parity2 = st2[i][0] as i8;
                at[i].stereo_atom_parity2 = st2[i][2] as i8;
                for j in 0..3 {
                    at[i].stereo_bond_neighbor[j] = st1[i][5 + 4 * j] as u16;
                    at[i].stereo_bond_ord[j] = st1[i][6 + 4 * j] as i8;
                    at[i].stereo_bond_z_prod[j] = st1[i][7 + 4 * j] as i8;
                    at[i].stereo_bond_parity[j] = st1[i][8 + 4 * j] as i8;
                    at[i].stereo_bond_neighbor2[j] = st2[i][5 + 4 * j] as u16;
                    at[i].stereo_bond_ord2[j] = st2[i][6 + 4 * j] as i8;
                    at[i].stereo_bond_z_prod2[j] = st2[i][7 + 4 * j] as i8;
                    at[i].stereo_bond_parity2[j] = st2[i][8 + 4 * j] as i8;
                }
            }
            at
        };
        let mut at_y = build("atY");
        let mut at_n = build("atN");
        let sz = ints(&inp["sizes"]);
        let mut tg = vec![TGroup::default(); ntg + 2];
        for (k, t) in inp["tgs"].as_array().unwrap().iter().enumerate() {
            let t = ints(t);
            tg[k].n_num_endpoints = t[0] as u16;
            tg[k].n_first_endpoint_at_no_pos = t[1] as u16;
            tg[k].n_group_number = k as u16 + 1;
            for j in 0..5 {
                tg[k].num[j] = t[2 + j] as u16;
            }
            tg[k].i_weight = (t[2] * 7 + t[4] * 3 + t[5]) % 5;
        }
        let iso = ints(&inp["iso"]);
        let iso_list = ints(&inp["isolist"]);
        let mut tgi = TGroupInfo {
            t_group: (iso[5] == 0).then_some(tg),
            n_endpoint_atom_number: ints(&inp["epl"]).iter().map(|&x| x as u16).collect(),
            num_t_groups: ntg,
            max_num_t_groups: ntg,
            n_num_endpoints: n_e as i32,
            n_num_isotopic_endpoints: iso[1] as i32,
            b_taut_flags_done: iso[2] as u64,
            b_ignore_isotopic: iso[4] as i32,
            t_group_number: vec![0; 4 * ntg + 4],
            ..Default::default()
        };
        tgi.tni.b_normalization_flags = iso[3] as u64;
        tgi.n_isotopic_endpoint_atom_number = std::iter::once(iso[0] as u16)
            .chain(iso_list.iter().map(|&x| x as u16))
            .collect();
        let mode_bits = inp["modeBits"].as_i64().unwrap();
        let uu = inp["uu"].as_i64().unwrap() != 0;
        let present = [have_non, have_yes];
        let mut s = [AtomSizes::default(); 2];
        for k in 0..2 {
            let v = &sz[7 * k..7 * k + 7];
            s[k].n_len_linear_ct_tautomer = v[4] as i32;
            if present[k] {
                let a = if k == TAUT_YES { &at_y } else { &at_n };
                get_canon_lengths(n, a, &mut s[k], have_tgi.then_some(&tgi));
            }
            s[k].num_isotopic_atoms = v[2] as i32;
            s[k].b_ignore_isotopic = v[3] as i32;
            s[k].b_has_isotopic_taut_groups = v[5] as i32;
            s[k].n_len_isotopic_endpoints = (v[6] as i32).max(s[k].n_len_isotopic_endpoints);
        }
        let mut bcn = Bcn::default();
        let mut got: Vec<String> = Vec::new();
        let at_refs: [Option<&[SpAtom]>; 2] =
            [have_non.then_some(&at_n[..]), have_yes.then_some(&at_y[..])];
        let ret = get_base_canon_ranking(
            &mut || false,
            n,
            n + ntg,
            at_refs,
            have_tgi.then_some(&mut tgi),
            &s,
            &mut bcn,
            inp["bfix"].as_i64().unwrap() != 0,
            false,
        );
        got.push(format!("R {ret}"));
        if ret == 0 {
            for k in 0..2 {
                if bcn.ftcn[k].linear_ct.is_empty() {
                    continue;
                }
                let (nat, na) = (bcn.ftcn[k].num_at_tg, bcn.ftcn[k].num_atoms);
                let (n_len_ct, n_len_ct_ao) = (
                    bcn.ftcn[k].n_len_linear_ct,
                    bcn.ftcn[k].n_len_linear_ct_at_only,
                );
                let k_present = if k == TAUT_YES { have_yes } else { have_non };
                let use_y = if !k_present {
                    k != TAUT_YES
                } else {
                    k == TAUT_YES
                };
                let at: &mut Vec<SpAtom> = if use_y { &mut at_y } else { &mut at_n };
                let mut cs = CanonStat::default();
                let mut mode = CMODE_CT | mode_bits;
                if nat > na {
                    mode |= CMODE_TAUT;
                }
                let (nds, ncs) = (
                    s[k].n_len_linear_ct_stereo_dble,
                    s[k].n_len_linear_ct_stereo_carb,
                );
                let (ids, ics);
                {
                    let mut si = AtomSizes::default();
                    let mut sw = false;
                    switch_atom_stereo_and_isotopic_stereo(at, na, &mut sw);
                    get_canon_lengths(na, at, &mut si, None);
                    switch_atom_stereo_and_isotopic_stereo(at, na, &mut sw);
                    ids = si.n_len_linear_ct_stereo_dble;
                    ics = si.n_len_linear_ct_stereo_carb;
                }
                let len_iso = at[..na]
                    .iter()
                    .filter(|a| {
                        let k = a.endpoint == 0
                            && a.c_flags & 1 == 0
                            && (a.num_iso_h[0] != 0 || a.num_iso_h[1] != 0 || a.num_iso_h[2] != 0);
                        a.iso_atw_diff != 0 || k
                    })
                    .count() as i32;
                allocate_cs(
                    &mut cs,
                    na,
                    nat,
                    n_len_ct,
                    n_len_ct_ao,
                    nds,
                    ids,
                    ncs,
                    ics,
                    (3 * ntg + n_e + 1) as i32,
                    if have_tgi && tgi.n_num_isotopic_endpoints != 0 {
                        0
                    } else {
                        tgi.t_group.as_ref().map_or(0, |t| {
                            t[..ntg].iter().filter(|g| g.i_weight != 0).count() as i32
                        })
                    },
                    len_iso,
                    mode,
                    true,
                );
                cs.p_bcn = Some(Box::new(std::mem::take(&mut bcn)));
                cs.t_group_info = have_tgi.then(|| std::mem::take(&mut tgi));
                cs.b_ignore_isotopic = 0;
                let r2 = canon_inchi(
                    &mut || false,
                    na,
                    nat,
                    at,
                    &mut cs,
                    mode | if uu { REQ_MODE_DIFF_UU_STEREO } else { 0 },
                    k,
                );
                bcn = *cs.p_bcn.take().unwrap();
                if let Some(t) = cs.t_group_info.take() {
                    tgi = t;
                }
                got.push(format!(
                    "C{k} {r2} {mode} {} {} {} {}",
                    cs.n_len_linear_ct,
                    cs.n_len_linear_ct_at_only,
                    cs.n_canon_flags,
                    cs.b_ignore_isotopic
                ));
                got.push(format!(
                    "SL {}{}{}{}",
                    slot_name(&cs, true, cs.stereo_slots.dble),
                    slot_name(&cs, false, cs.stereo_slots.carb),
                    slot_name(&cs, true, cs.stereo_slots.dble_inv),
                    slot_name(&cs, false, cs.stereo_slots.carb_inv)
                ));
                got.push(pu("LC", &cs.linear_ct, cs.n_max_len_linear_ct as usize));
                got.push(pu("L2", &cs.linear_ct2, cs.n_max_len_linear_ct as usize));
                got.push(format!(
                    "LE {} {} {} {} {} {} {} {}",
                    cs.n_len_linear_ct2,
                    cs.n_len_linear_ct_at_only2,
                    cs.n_len_linear_ct_tautomer,
                    cs.n_len_linear_ct_isotopic,
                    cs.n_len_linear_ct_isotopic_tautomer,
                    cs.b_cmp_stereo,
                    cs.b_cmp_isotopic_stereo,
                    cs.n_mode
                ));
                let mut l = String::from("NH");
                for x in cs.n_num_h.iter().take(na) {
                    l += &format!(" {x}");
                }
                got.push(l);
                let mut l = String::from("NF");
                for x in cs.n_num_h_fixed.iter().take(na) {
                    l += &format!(" {x}");
                }
                got.push(l);
                let mut l = String::from("NX");
                for x in cs.n_exchg_iso_h.iter().take(na) {
                    l += &format!(" {x}");
                }
                got.push(l);
                got.push(pu(
                    "LT",
                    &cs.linear_ct_tautomer,
                    cs.n_max_len_linear_ct_tautomer as usize,
                ));
                let mut l = String::from("IT");
                for e in cs
                    .linear_ct_isotopic_tautomer
                    .iter()
                    .take(cs.n_max_len_linear_ct_isotopic_tautomer as usize)
                {
                    l += &format!(" {}:{},{},{}", e.tgroup_num, e.num[0], e.num[1], e.num[2]);
                }
                got.push(l);
                let mut l = String::from("IA");
                for e in cs
                    .linear_ct_isotopic
                    .iter()
                    .take(cs.n_max_len_linear_ct_isotopic as usize)
                {
                    l += &format!(
                        " {},{},{},{},{}",
                        e.at_num, e.num_1h, e.num_d, e.num_t, e.iso_atw_diff
                    );
                }
                got.push(l);
                got.push(pu("SR", &cs.n_symm_rank, nat));
                got.push(pu("CO", &cs.n_canon_ord, nat));
                got.push(pu("SI", &cs.n_symm_rank_isotopic, nat));
                got.push(pu("CI", &cs.n_canon_ord_isotopic, nat));
                got.push(pu("ST", &cs.n_canon_ord_stereo, nat));
                got.push(pu("SV", &cs.n_canon_ord_stereo_inv, nat));
                got.push(pu("IS", &cs.n_canon_ord_isotopic_stereo, nat));
                got.push(pu("IV", &cs.n_canon_ord_isotopic_stereo_inv, nat));
                let ntgs = if have_tgi { ntg } else { 0 };
                got.push(pu("TO", &cs.n_canon_ord_taut, ntgs));
                got.push(pu("TS", &cs.n_symm_rank_taut, ntgs));
                got.push(pu("TI", &cs.n_canon_ord_isotopic_taut, ntgs));
                got.push(pu("TJ", &cs.n_symm_rank_isotopic_taut, ntgs));
                got.push(pu("TT", &cs.n_canon_ord_stereo_taut, ntgs));
                got.push(pu("TU", &cs.n_canon_ord_isotopic_stereo_taut, ntgs));
                got.push(format!(
                    "CL {} {} {} {} {} {} {} {}",
                    cs.n_len_canon_ord,
                    cs.n_len_canon_ord_taut,
                    cs.n_len_canon_ord_isotopic,
                    cs.n_len_canon_ord_isotopic_taut,
                    cs.n_len_canon_ord_stereo,
                    cs.n_len_canon_ord_stereo_taut,
                    cs.n_len_canon_ord_isotopic_stereo,
                    cs.n_len_canon_ord_isotopic_stereo_taut
                ));
                got.push(format!(
                    "SS {} {} {} {} {} {} {} {}",
                    cs.n_len_linear_ct_stereo_dble,
                    cs.n_len_linear_ct_stereo_carb,
                    cs.n_len_linear_ct_stereo_dble_inv,
                    cs.n_len_linear_ct_stereo_carb_inv,
                    cs.n_len_linear_ct_isotopic_stereo_dble,
                    cs.n_len_linear_ct_isotopic_stereo_carb,
                    cs.n_len_linear_ct_isotopic_stereo_dble_inv,
                    cs.n_len_linear_ct_isotopic_stereo_carb_inv
                ));
                got.push(format!(
                    "SM {} {} {} {}",
                    cs.n_max_len_linear_ct_stereo_dble,
                    cs.n_max_len_linear_ct_stereo_carb,
                    cs.n_max_len_linear_ct_isotopic_stereo_dble,
                    cs.n_max_len_linear_ct_isotopic_stereo_carb
                ));
                let slots = [Slot::Base, Slot::Inv, Slot::Iso, Slot::IsoInv];
                let (ad, ac) = ([nds, nds, ids, ids], [ncs, ncs, ics, ics]);
                for q in 0..4 {
                    let mut l = format!("D{q}");
                    let b = cs.dble_ref(slots[q]);
                    if !b.is_empty() {
                        for e in b.iter().take(ad[q] as usize) {
                            l += &format!(" {},{},{}", e.at_num1, e.at_num2, e.parity);
                        }
                    }
                    got.push(l);
                    let mut l = format!("K{q}");
                    let b = cs.carb_ref(slots[q]);
                    if !b.is_empty() {
                        for e in b.iter().take(ac[q] as usize) {
                            l += &format!(" {},{}", e.at_num, e.parity);
                        }
                    }
                    got.push(l);
                }
                for a in at.iter().take(n) {
                    let mut l = format!(
                        "A {} {} {} {} {} {} {} {} |",
                        a.parity,
                        a.parity2,
                        a.final_parity,
                        a.final_parity2,
                        a.stereo_atom_parity,
                        a.stereo_atom_parity2,
                        a.b_has_stereo_or_equ_to_stereo,
                        a.b_has_stereo_or_equ_to_stereo2
                    );
                    for j in 0..3 {
                        l += &format!(
                            " {},{},{},{}",
                            a.stereo_bond_neighbor[j],
                            a.stereo_bond_parity[j],
                            a.stereo_bond_ord[j],
                            a.stereo_bond_z_prod[j]
                        );
                    }
                    l += " |";
                    for j in 0..3 {
                        l += &format!(
                            " {},{},{},{}",
                            a.stereo_bond_neighbor2[j],
                            a.stereo_bond_parity2[j],
                            a.stereo_bond_ord2[j],
                            a.stereo_bond_z_prod2[j]
                        );
                    }
                    got.push(l);
                }
                let mut l = String::from("TG");
                if have_tgi {
                    for x in tgi.t_group_number.iter().take(4 * ntgs) {
                        l += &format!(" {x}");
                    }
                    l += " |";
                    for x in tgi.n_endpoint_atom_number.iter().take(n_e) {
                        l += &format!(" {x}");
                    }
                }
                got.push(l);
                got.push(format!(
                    "TI {}",
                    if have_tgi { tgi.b_ignore_isotopic } else { -1 }
                ));
                got.push(format!(
                    "PN {}",
                    i32::from(bcn.ftcn[k].neigh_list.is_some())
                ));
                got.push(format!("RS {}", bcn.n_max_len_rank_stack));
                for j in 0..(bcn.n_max_len_rank_stack / 2) as usize {
                    let p = bcn.p_rank_stack.get(j).cloned().unwrap_or_default();
                    let nm = bcn.num_max;
                    got.push(pu(
                        &format!("Q{j}"),
                        &p.rank,
                        if p.rank.is_empty() { 0 } else { nm },
                    ));
                    got.push(pu(
                        &format!("W{j}"),
                        &p.at_number,
                        if p.at_number.is_empty() { 0 } else { nm },
                    ));
                }
                if r2 > 0 {
                    use rdoxide::external::inchi::ichi::{
                        InChI, InChIAux, InChIIsotopicAtom, InChIIsotopicTGroup, InChIStereo,
                        OrigInfo,
                    };
                    use rdoxide::external::inchi::ichimak2::fill_out_inchi;
                    use rdoxide::external::inchi::inpdef::InpAtom;
                    use rdoxide::external::inchi::util::get_periodic_table_number;
                    let umodes = [0u64, 0x200, 0x400, 0x800, 0x1000, 0x1800];
                    let um = umodes[case % 6];
                    let cap = 4 * (nat + n_e) + 40;
                    let mk_stereo = || InChIStereo {
                        n_number: vec![0; 64],
                        t_parity: vec![0; 64],
                        n_number_inv: vec![0; 64],
                        t_parity_inv: vec![0; 64],
                        n_bond_atom1: vec![0; 64],
                        n_bond_atom2: vec![0; 64],
                        b_parity: vec![0; 64],
                        ..Default::default()
                    };
                    let mut pi = InChI {
                        n_atom: vec![0; nat + 2],
                        n_conn_table: vec![0; n_len_ct as usize + nat + 10],
                        n_tautomer: vec![0; cap],
                        n_num_h: vec![0; nat + 2],
                        n_num_h_fixed: Some(vec![0; nat + 2]),
                        isotopic_atom: vec![InChIIsotopicAtom::default(); na + 2],
                        isotopic_t_group: vec![InChIIsotopicTGroup::default(); ntg + 2],
                        n_possible_locations_of_isotopic_h: vec![0; na + 3],
                        stereo: Some(mk_stereo()),
                        stereo_isotopic: Some(mk_stereo()),
                        ..Default::default()
                    };
                    let al = 2 * na + 8;
                    let mut pa = InChIAux {
                        n_orig_at_nos_in_canon_ord: vec![0; al],
                        n_isotopic_orig_at_nos_in_canon_ord: vec![0; al],
                        n_orig_at_nos_in_canon_ord_inv: vec![0; al],
                        n_isotopic_orig_at_nos_in_canon_ord_inv: vec![0; al],
                        n_constit_equ_numbers: vec![0; al],
                        n_constit_equ_isotopic_numbers: vec![0; al],
                        n_constit_equ_t_group_numbers: vec![0; al],
                        n_constit_equ_isotopic_t_group_numbers: vec![0; al],
                        orig_info: vec![OrigInfo::default(); al],
                        ..Default::default()
                    };
                    if have_tgi {
                        pa.b_normalization_flags = tgi.tni.b_normalization_flags;
                        pa.b_taut_flags = tgi.b_taut_flags;
                        pa.b_taut_flags_done = tgi.b_taut_flags_done;
                    }
                    let mut nrm = vec![InpAtom::default(); na + 2];
                    for i in 0..na {
                        let el = at[i].elname.iter().take_while(|&&b| b != 0).count();
                        at[i].el_number = get_periodic_table_number(
                            std::str::from_utf8(&at[i].elname[..el]).unwrap(),
                        ) as u8;
                        at[i].charge = ((i * 7 + 3) % 4) as i8 - 1;
                        at[i].radical = if i % 5 == 4 {
                            1
                        } else if i % 7 == 3 {
                            2
                        } else {
                            0
                        };
                        at[i].orig_at_number = 101 + i as u16;
                        at[i].chem_bonds_valence = at[i].valence;
                        nrm[i].valence = at[i].valence;
                        nrm[i].num_h = at[i].num_h;
                        nrm[i].charge = at[i].charge;
                        nrm[i].radical = at[i].radical;
                        nrm[i].el_number = at[i].el_number;
                        nrm[i].chem_bonds_valence = at[i].chem_bonds_valence;
                    }
                    let mut se = String::new();
                    let fr = fill_out_inchi(
                        &mut pi,
                        &mut pa,
                        na,
                        nat,
                        0,
                        at,
                        &mut nrm,
                        &mut cs,
                        have_tgi.then_some(&mut tgi),
                        k != 0,
                        um,
                        &mut se,
                        false,
                    );
                    got.push(format!(
                        "FO {fr} {} {} {} {} {} {}",
                        pi.n_error_code,
                        pa.n_error_code,
                        pi.n_flags,
                        pi.n_total_charge,
                        pi.n_number_of_atoms,
                        pi.sz_hill_formula.as_deref().unwrap_or("NULL")
                    ));
                    got.push(format!(
                        "FA{}",
                        pi.n_atom[..na]
                            .iter()
                            .map(|x| format!(" {x}"))
                            .collect::<String>()
                    ));
                    got.push(format!(
                        "FC {}{}",
                        pi.len_conn_table,
                        pi.n_conn_table[..pi.len_conn_table as usize]
                            .iter()
                            .map(|x| format!(" {x}"))
                            .collect::<String>()
                    ));
                    got.push(format!(
                        "FT {}{}",
                        pi.len_tautomer,
                        pi.n_tautomer[..pi.len_tautomer as usize]
                            .iter()
                            .map(|x| format!(" {x}"))
                            .collect::<String>()
                    ));
                    let fixed = pi.n_num_h_fixed.as_ref().unwrap();
                    got.push(format!(
                        "FH{}",
                        (0..na)
                            .map(|i| format!(" {}/{}", pi.n_num_h[i], fixed[i]))
                            .collect::<String>()
                    ));
                    let mut l = format!(
                        "FI {} {}",
                        pi.n_number_of_isotopic_atoms, pi.n_number_of_isotopic_t_groups
                    );
                    for e in pi
                        .isotopic_atom
                        .iter()
                        .take(pi.n_number_of_isotopic_atoms.max(0) as usize)
                    {
                        l += &format!(
                            " {},{},{},{},{}",
                            e.n_atom_number, e.n_iso_difference, e.n_num_h, e.n_num_d, e.n_num_t
                        );
                    }
                    l += " |";
                    for e in pi
                        .isotopic_t_group
                        .iter()
                        .take(pi.n_number_of_isotopic_t_groups.max(0) as usize)
                    {
                        l += &format!(
                            " {},{},{},{}",
                            e.n_t_group_number, e.n_num_h, e.n_num_d, e.n_num_t
                        );
                    }
                    got.push(l);
                    got.push(format!(
                        "FL{}",
                        pi.n_possible_locations_of_isotopic_h
                            .iter()
                            .map(|x| format!(" {x}"))
                            .collect::<String>()
                    ));
                    for (q, st) in [&pi.stereo, &pi.stereo_isotopic].iter().enumerate() {
                        let s = st.as_ref().unwrap();
                        let mut l = format!(
                            "FS{q} {} {} {} {} |",
                            s.n_number_of_stereo_centers,
                            s.n_number_of_stereo_bonds,
                            s.n_comp_inv2_abs,
                            s.b_trivial_inv
                        );
                        for i in 0..64 {
                            l += &format!(
                                " {},{},{},{}",
                                s.n_number[i], s.t_parity[i], s.n_number_inv[i], s.t_parity_inv[i]
                            );
                        }
                        l += " |";
                        for i in 0..64 {
                            l += &format!(
                                " {},{},{}",
                                s.n_bond_atom1[i], s.n_bond_atom2[i], s.b_parity[i]
                            );
                        }
                        got.push(l);
                    }
                    got.push(format!(
                        "FX {} {} {} {} {} {}",
                        pa.n_number_of_atoms,
                        pa.n_number_of_t_groups,
                        pa.n_num_removed_protons,
                        pa.n_num_removed_isotopic_h[0],
                        pa.n_num_removed_isotopic_h[1],
                        pa.n_num_removed_isotopic_h[2]
                    ));
                    let aa = [
                        &pa.n_orig_at_nos_in_canon_ord,
                        &pa.n_isotopic_orig_at_nos_in_canon_ord,
                        &pa.n_orig_at_nos_in_canon_ord_inv,
                        &pa.n_isotopic_orig_at_nos_in_canon_ord_inv,
                        &pa.n_constit_equ_numbers,
                        &pa.n_constit_equ_isotopic_numbers,
                        &pa.n_constit_equ_t_group_numbers,
                        &pa.n_constit_equ_isotopic_t_group_numbers,
                    ];
                    for (q, a) in aa.iter().enumerate() {
                        got.push(format!(
                            "FQ{q}{}",
                            a.iter().map(|x| format!(" {x}")).collect::<String>()
                        ));
                    }
                    got.push(format!(
                        "FO2{}",
                        (0..na)
                            .map(|i| format!(
                                " {},{},{}",
                                pa.orig_info[i].c_charge,
                                pa.orig_info[i].c_radical,
                                pa.orig_info[i].c_unusual_valence
                            ))
                            .collect::<String>()
                    ));
                    got.push(format!("FE {se}"));
                    got.push(format!(
                        "FB{}",
                        (0..na)
                            .map(|i| format!(
                                " {}/{}",
                                at[i].b_ambiguous_stereo, nrm[i].b_ambiguous_stereo
                            ))
                            .collect::<String>()
                    ));
                    got.push(pu("FY", &cs.n_canon_ord_stereo, nat));
                    got.push(pu("FZ", &cs.n_canon_ord_isotopic_stereo, nat));
                    let mut l = String::from("FG");
                    if have_tgi {
                        for x in tgi.t_group_number.iter().take(4 * ntgs) {
                            l += &format!(" {x}");
                        }
                        l += " |";
                        for x in tgi.n_endpoint_atom_number.iter().take(n_e) {
                            l += &format!(" {x}");
                        }
                    }
                    got.push(l);
                    cs.t_group_info = have_tgi.then(|| std::mem::take(&mut tgi));
                    let kr = rdoxide::external::inchi::ichimake::check_canon_numbering_correctness(
                        na, nat, at, &mut cs,
                    );
                    if let Some(t) = cs.t_group_info.take() {
                        tgi = t;
                    }
                    got.push(format!("FK {kr}"));
                    got.push(pu("FM", &cs.linear_ct, cs.n_len_linear_ct.max(0) as usize));
                    got.push(pu(
                        "FN",
                        &cs.linear_ct2,
                        cs.n_len_linear_ct2.max(0) as usize,
                    ));
                    got.push(format!(
                        "FP {} {} {}",
                        cs.n_len_linear_ct, cs.n_len_linear_ct2, cs.n_len_linear_ct_at_only2
                    ));
                }
            }
        }
        let want: Vec<String> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        assert_eq!(got, want, "case {case}");
    }
}

#[test]
fn atom_copy_isotopic_keys_and_hill_formulas_match_the_c_code() {
    use rdoxide::external::inchi::ichican2::set_atom_iso_sort_keys;
    use rdoxide::external::inchi::ichimake::{
        compare_hill_formulas, compare_hill_formulas_no_h, inp2sp_atom,
    };
    use rdoxide::external::inchi::ichister::SpAtom;
    use rdoxide::external::inchi::ichitaut::{TGroup, TGroupInfo};
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn sgn(x: i32) -> i32 {
        x.signum()
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_make1.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut lines = inp.lines();
        let head: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
        let mut got: Vec<String> = Vec::new();
        match head[0] {
            "A" => {
                let n: usize = head[1].parse().unwrap();
                let mut ia = vec![InpAtom::default(); n];
                for a in ia.iter_mut() {
                    let t: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
                    let el = t[0].as_bytes();
                    a.elname = [0; 6];
                    a.elname[..el.len().min(6)].copy_from_slice(&el[..el.len().min(6)]);
                    let v: usize = t[1].parse().unwrap();
                    a.valence = v as i8;
                    for j in 0..v {
                        a.neighbor[j] = t[2 + 2 * j].parse().unwrap();
                        a.bond_type[j] = t[3 + 2 * j].parse().unwrap();
                    }
                    let x: Vec<i32> = t[2 + 2 * v..].iter().map(|s| s.parse().unwrap()).collect();
                    a.chem_bonds_valence = x[0] as i8;
                    a.orig_at_number = x[1] as u16;
                    a.orig_compt_at_numb = x[2] as u16;
                    a.endpoint = x[3] as u16;
                    a.iso_atw_diff = x[4] as i8;
                    a.num_h = x[5] as i8;
                    a.c_flags = x[6] as i8;
                    a.num_iso_h = [x[7] as i8, x[8] as i8, x[9] as i8];
                    a.charge = x[10] as i8;
                    a.radical = x[11] as i8;
                    a.n_block_system = x[12] as u16;
                    a.b_cut_vertex = x[13] as i8;
                    a.n_ring_system = x[14] as u16;
                    a.n_num_at_in_ring_system = x[15] as u16;
                }
                let mut at = vec![SpAtom::default(); n];
                inp2sp_atom(&ia, &mut at);
                for a in &at {
                    let el = a.elname.iter().take_while(|&&b| b != 0).count();
                    let mut s = format!(
                        "{} {} {} |",
                        std::str::from_utf8(&a.elname[..el]).unwrap(),
                        a.el_number,
                        a.valence
                    );
                    for j in 0..a.valence as usize {
                        s += &format!(" {} {}", a.neighbor[j], a.bond_type[j]);
                    }
                    s += &format!(
                        " | {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                        a.chem_bonds_valence,
                        a.orig_at_number,
                        a.orig_compt_at_numb,
                        a.endpoint,
                        a.iso_atw_diff,
                        a.num_h,
                        a.c_flags,
                        a.num_iso_h[0],
                        a.num_iso_h[1],
                        a.num_iso_h[2],
                        a.charge,
                        a.radical,
                        a.n_block_system,
                        a.b_cut_vertex,
                        a.n_ring_system,
                        a.n_num_at_in_ring_system
                    );
                    got.push(s);
                }
            }
            "B" => {
                let n: usize = head[1].parse().unwrap();
                let hastg: i32 = head[2].parse().unwrap();
                let ng: usize = head[3].parse().unwrap();
                let hasiso: i32 = head[4].parse().unwrap();
                let wantflag: i32 = head[5].parse().unwrap();
                let mut at = vec![SpAtom::default(); n];
                for a in at.iter_mut() {
                    let x: Vec<i32> = lines
                        .next()
                        .unwrap()
                        .split_whitespace()
                        .map(|s| s.parse().unwrap())
                        .collect();
                    a.endpoint = x[0] as u16;
                    a.iso_atw_diff = x[1] as i8;
                    a.c_flags = x[2] as i8;
                    a.num_iso_h = [x[3] as i8, x[4] as i8, x[5] as i8];
                }
                let mut tgi = TGroupInfo::default();
                if hastg != 0 {
                    tgi.t_group = Some(vec![TGroup::default(); 2]);
                    tgi.num_t_groups = ng;
                }
                if hasiso != 0 {
                    tgi.n_isotopic_endpoint_atom_number = vec![1, 2];
                }
                let mut flag = -7;
                let nis = set_atom_iso_sort_keys(
                    &mut at,
                    if hastg >= 0 { Some(&tgi) } else { None },
                    if wantflag != 0 { Some(&mut flag) } else { None },
                );
                got.push(format!("{nis} {flag}"));
                for a in &at {
                    got.push(format!("{}", a.iso_sort_key));
                }
            }
            _ => {
                let t: Vec<&str> = inp.split_whitespace().collect();
                let f = |s: &str| {
                    if s == "-" {
                        String::new()
                    } else {
                        s.to_string()
                    }
                };
                let (f1, f2) = (f(t[1]), f(t[2]));
                let (mut a, mut b) = (0, 0);
                let r1 = compare_hill_formulas(&f1, &f2);
                let r2 = compare_hill_formulas_no_h(&f1, &f2, &mut a, &mut b);
                got.push(format!("{} {} {a} {b}", sgn(r1), sgn(r2)));
            }
        }
        assert_eq!(got, want, "case {ci}: {inp}");
    }
}

use rdoxide::external::inchi::ichi::{InChI, InChIIsotopicAtom, InChIIsotopicTGroup, InChIStereo};

struct InchiTok<'a>(std::str::SplitWhitespace<'a>);
#[allow(clippy::field_reassign_with_default)]
impl InchiTok<'_> {
    fn i(&mut self) -> i64 {
        self.0.next().unwrap().parse().unwrap()
    }
    fn s(&mut self) -> String {
        self.0.next().unwrap().to_string()
    }
    fn stereo(&mut self) -> Option<InChIStereo> {
        if self.i() == 0 {
            return None;
        }
        let mut s = InChIStereo::default();
        let n = self.i();
        s.n_number_of_stereo_centers = n as i32;
        for _ in 0..n {
            s.n_number.push(self.i() as u16);
            s.t_parity.push(self.i() as i8);
        }
        let n = self.i();
        s.n_number_of_stereo_bonds = n as i32;
        for _ in 0..n {
            s.n_bond_atom1.push(self.i() as u16);
            s.n_bond_atom2.push(self.i() as u16);
            s.b_parity.push(self.i() as i8);
        }
        s.n_comp_inv2_abs = self.i() as i32;
        Some(s)
    }
    fn inchi(&mut self) -> Option<InChI> {
        if self.i() == 0 {
            return None;
        }
        let mut c = InChI::default();
        c.n_flags = self.i() as u64;
        c.n_total_charge = self.i() as i32;
        let n = self.i();
        c.n_number_of_atoms = n as i32;
        c.sz_hill_formula = Some(self.s());
        c.n_atom = (0..n).map(|_| self.i() as u8).collect();
        let l = self.i();
        c.len_conn_table = l as i32;
        c.n_conn_table = (0..l).map(|_| self.i() as u16).collect();
        let l = self.i();
        c.len_tautomer = l as i32;
        c.n_tautomer = (0..l).map(|_| self.i() as u16).collect();
        c.n_num_h = (0..n).map(|_| self.i() as i8).collect();
        if self.i() != 0 {
            c.n_num_h_fixed = Some((0..n).map(|_| self.i() as i8).collect());
        }
        let l = self.i();
        c.n_number_of_isotopic_atoms = l as i32;
        for _ in 0..l {
            c.isotopic_atom.push(InChIIsotopicAtom {
                n_atom_number: self.i() as u16,
                n_iso_difference: self.i() as i16,
                n_num_h: self.i() as i16,
                n_num_d: self.i() as i16,
                n_num_t: self.i() as i16,
            });
        }
        let l = self.i();
        c.n_number_of_isotopic_t_groups = l as i32;
        for _ in 0..l {
            c.isotopic_t_group.push(InChIIsotopicTGroup {
                n_t_group_number: self.i() as u16,
                n_num_h: self.i() as u16,
                n_num_d: self.i() as u16,
                n_num_t: self.i() as u16,
            });
        }
        c.stereo = self.stereo();
        c.stereo_isotopic = self.stereo();
        c.b_deleted = self.i() as i32;
        Some(c)
    }
}

#[test]
fn inchi_component_comparison_matches_the_c_code() {
    use rdoxide::external::inchi::ichi::{
        InChISort, comp_inchi_non_taut2, comp_inchi_taut_vs_non_taut, comp_inchi_taut2,
        comp_inchi2, compare_inchi_stereo, compare_taut_non_iso_part_of_inchi, get_sp3_rel_rac_abs,
    };
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_make2.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(inp.split_whitespace());
        let mut p = [InChISort::default(), InChISort::default()];
        for (q, ps) in p.iter_mut().enumerate() {
            ps.p_inchi = [t.inchi(), t.inchi()];
            ps.ord_number = q as i16;
        }
        let sg = |x: i32| x.signum();
        let a = p[0].p_inchi[1].as_ref().or(p[0].p_inchi[0].as_ref());
        let b = p[1].p_inchi[1].as_ref().or(p[1].p_inchi[0].as_ref());
        let mut got = Vec::new();
        got.push(match (a, b) {
            (Some(a), Some(b)) => format!("{}", sg(compare_taut_non_iso_part_of_inchi(a, b))),
            _ => "0".to_string(),
        });
        got.push(match a {
            Some(a) => format!(
                "{} {}",
                get_sp3_rel_rac_abs(Some(a), a.stereo.as_ref()),
                get_sp3_rel_rac_abs(Some(a), a.stereo_isotopic.as_ref())
            ),
            None => "0 0".to_string(),
        });
        got.push(match (a, b) {
            (Some(a), Some(b)) => format!(
                "{} {}",
                sg(compare_inchi_stereo(
                    a.stereo.as_ref(),
                    a.n_flags,
                    b.stereo.as_ref(),
                    b.n_flags
                )),
                sg(compare_inchi_stereo(
                    a.stereo_isotopic.as_ref(),
                    a.n_flags,
                    b.stereo.as_ref(),
                    b.n_flags
                ))
            ),
            _ => "0 0".to_string(),
        });
        got.push(format!(
            "{} {}",
            sg(comp_inchi_taut_vs_non_taut(&p[0], &p[1], true)),
            sg(comp_inchi_taut_vs_non_taut(&p[0], &p[1], false))
        ));
        got.push(format!(
            "{} {} {} {}",
            sg(comp_inchi2(&p[0], &p[1], 0, true)),
            sg(comp_inchi2(&p[0], &p[1], 1, true)),
            sg(comp_inchi2(&p[0], &p[1], 0, false)),
            sg(comp_inchi2(&p[0], &p[1], 1, false))
        ));
        got.push(format!(
            "{} {}",
            sg(comp_inchi_non_taut2(&p[0], &p[1])),
            sg(comp_inchi_taut2(&p[0], &p[1]))
        ));
        assert_eq!(got, want, "case {ci}: {inp}");
    }
}

#[test]
fn inchi_layer_difference_marking_matches_the_c_code() {
    use rdoxide::external::inchi::ichi::{
        DIFL_LENGTH, DIFS_LENGTH, DifSegs, InChISort, comp_inchi_layers, eql_inchi_stereo,
        inchi_segment_action, mark_unused_and_empty_layers,
    };
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_make3.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(inp.split_whitespace());
        let (ord1, ord2, fix) = (t.i(), t.i(), t.i());
        let mut p = [InChISort::default(), InChISort::default()];
        for ps in p.iter_mut() {
            ps.p_inchi = [t.inchi(), t.inchi()];
        }
        p[0].ord_number = ord1 as i16;
        p[1].ord_number = ord2 as i16;
        let a = p[0].p_inchi[1].as_ref().or(p[0].p_inchi[0].as_ref());
        let b = p[1].p_inchi[1].as_ref().or(p[1].p_inchi[0].as_ref());
        let mut got = Vec::new();
        got.push(match (a, b) {
            (Some(a), Some(b)) => {
                let e = [8, 2, 4, 1];
                let mut s = String::new();
                for x in e {
                    for y in e {
                        s += &i32::from(eql_inchi_stereo(
                            a.stereo.as_ref(),
                            x,
                            b.stereo_isotopic.as_ref(),
                            y,
                        ))
                        .to_string();
                    }
                }
                s
            }
            _ => "-".to_string(),
        });
        let mut s: DifSegs = [[0; DIFS_LENGTH]; DIFL_LENGTH];
        let ret = comp_inchi_layers(&p[0], &p[1], &mut s, fix as i32);
        let flat = |s: &DifSegs| -> String {
            s.iter()
                .flatten()
                .map(|v| format!(" {v}"))
                .collect::<String>()
        };
        got.push(format!("{ret}{}", flat(&s)));
        mark_unused_and_empty_layers(&mut s);
        got.push(flat(&s));
        got.push(
            (0..16)
                .map(|v| format!(" {}", inchi_segment_action(v)))
                .collect(),
        );
        assert_eq!(got, want, "case {ci}: {inp}");
    }
}

#[test]
fn dfs_order_of_the_connection_table_matches_the_c_code() {
    use rdoxide::external::inchi::ichimake::{get_dfs_order_4ct, mystrrev};
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_make4.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(inp.split_whitespace());
        let n = t.i() as usize;
        let len = t.i() as usize;
        let mode = t.i() as i32;
        let has_h = t.i() != 0;
        let ct: Vec<u16> = (0..len).map(|_| t.i() as u16).collect();
        let h: Vec<i8> = (0..n).map(|_| t.i() as i8).collect();
        let sl = t.i() as usize;
        let mut s: Vec<u8> = (0..sl).map(|_| t.i() as u8).collect();
        let first = match get_dfs_order_4ct(&ct, len, has_h.then_some(&h[..]), n, mode) {
            None => "NULL".to_string(),
            Some(o) => o
                .chunks(3)
                .take_while(|t| t[0] != 0)
                .map(|t| format!(" {} {} {}", t[0], t[1], t[2]))
                .collect(),
        };
        mystrrev(&mut s);
        let got = vec![first, String::from_utf8(s).unwrap()];
        assert_eq!(got, want, "case {ci}: {inp}");
    }
}

#[test]
fn hill_formula_and_stereo_copy_match_the_c_code() {
    use rdoxide::external::inchi::ichi::{InChI, InChIStereo};
    use rdoxide::external::inchi::ichimak2::{
        allocate_and_fill_hill_formula, copy_linear_ct_stereo_to_inchi_stereo,
        get_hill_formula_counts, mark_ambiguous_stereo, unmark_all_undefined_unknown_stereo,
    };
    use rdoxide::external::inchi::ichister::{AtStereoCarb, AtStereoDble, SpAtom};
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn pr_stereo(s: &InChIStereo) -> Vec<String> {
        let mut v = vec![format!(
            "S {} {} {} {}",
            s.n_number_of_stereo_centers,
            s.n_number_of_stereo_bonds,
            s.n_comp_inv2_abs,
            s.b_trivial_inv
        )];
        let mut c = "C".to_string();
        for i in 0..s.n_number_of_stereo_centers as usize {
            c += &format!(
                " {} {} {} {}",
                s.n_number[i], s.t_parity[i], s.n_number_inv[i], s.t_parity_inv[i]
            );
        }
        v.push(c);
        let mut b = "B".to_string();
        for i in 0..s.n_number_of_stereo_bonds as usize {
            b += &format!(
                " {} {} {}",
                s.n_bond_atom1[i], s.n_bond_atom2[i], s.b_parity[i]
            );
        }
        v.push(b);
        v
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_mak2a.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(inp.split_whitespace());
        let mode = t.i();
        let mut got: Vec<String> = Vec::new();
        if mode == 0 {
            let n = t.i() as usize;
            let n_atom: Vec<u8> = (0..n).map(|_| t.i() as u8).collect();
            let n_num_h: Vec<i8> = (0..n).map(|_| t.i() as i8).collect();
            let lt = t.i() as usize;
            let taut: Vec<u16> = (0..lt).map(|_| t.i() as u16).collect();
            match get_hill_formula_counts(&n_atom, &n_num_h, n, &taut, lt as i32) {
                None => got.push("E".to_string()),
                Some(h) => got.push(format!(
                    "{} {} {} {}",
                    h.num_c, h.num_h, h.len, h.num_non_h_atoms
                )),
            }
            let p = InChI {
                n_number_of_atoms: n as i32,
                n_atom,
                n_num_h,
                len_tautomer: lt as i32,
                n_tautomer: taut,
                ..Default::default()
            };
            got.push(allocate_and_fill_hill_formula(&p).unwrap_or_else(|| "NULL".to_string()));
        } else {
            let n = t.i() as usize;
            let mut at = vec![SpAtom::default(); n];
            let mut na = vec![InpAtom::default(); n];
            for i in 0..n {
                at[i].valence = t.i() as i8;
                for j in 0..3 {
                    at[i].neighbor[j] = t.i() as u16;
                }
                at[i].stereo_bond_parity[0] = t.i() as i8;
                at[i].stereo_bond_parity2[0] = t.i() as i8;
                at[i].stereo_bond_neighbor[1] = t.i() as u16;
                at[i].stereo_bond_neighbor2[1] = t.i() as u16;
                at[i].stereo_bond_ord[0] = t.i() as i8;
                at[i].stereo_bond_ord2[0] = t.i() as i8;
                at[i].b_ambiguous_stereo = t.i() as i8;
                na[i].b_ambiguous_stereo = at[i].b_ambiguous_stereo;
            }
            let iso = t.i() != 0;
            let hasord = t.i() != 0;
            let hasordinv = t.i() != 0;
            let v = |t: &mut InchiTok| -> Vec<u16> { (0..n).map(|_| t.i() as u16).collect() };
            let ord = v(&mut t);
            let rank = v(&mut t);
            let ordi = v(&mut t);
            let ranki = v(&mut t);
            let nca = t.i() as usize;
            let ndb = t.i() as usize;
            let rc = |t: &mut InchiTok| -> Vec<AtStereoCarb> {
                (0..nca)
                    .map(|_| AtStereoCarb {
                        at_num: t.i() as u16,
                        parity: t.i() as u8,
                    })
                    .collect()
            };
            let ca = rc(&mut t);
            let cai = rc(&mut t);
            let rd = |t: &mut InchiTok| -> Vec<AtStereoDble> {
                (0..ndb)
                    .map(|_| AtStereoDble {
                        at_num1: t.i() as u16,
                        at_num2: t.i() as u16,
                        parity: t.i() as u8,
                    })
                    .collect()
            };
            let db = rd(&mut t);
            let dbi = rd(&mut t);
            let umode = t.i() as u64;
            let cap = 48;
            let mut s = InChIStereo {
                n_number: vec![0; cap],
                t_parity: vec![0; cap],
                n_number_inv: vec![0; cap],
                t_parity_inv: vec![0; cap],
                n_bond_atom1: vec![0; cap],
                n_bond_atom2: vec![0; cap],
                b_parity: vec![0; cap],
                ..Default::default()
            };
            let r = copy_linear_ct_stereo_to_inchi_stereo(
                Some(&mut s),
                Some(&ca),
                nca,
                &db,
                ndb,
                hasord.then_some(&ord[..]),
                hasord.then_some(&rank[..]),
                &at,
                iso,
                Some(&cai),
                &dbi,
                hasordinv.then_some(&ordi[..]),
                hasordinv.then_some(&ranki[..]),
            );
            got.push(format!("{r}"));
            got.extend(pr_stereo(&s));
            let m =
                mark_ambiguous_stereo(&mut at, &mut na, iso, hasord.then_some(&ord[..]), &ca, &db);
            got.push(format!("{m}"));
            got.push(
                (0..n)
                    .map(|i| format!(" {}/{}", at[i].b_ambiguous_stereo, na[i].b_ambiguous_stereo))
                    .collect(),
            );
            let u = unmark_all_undefined_unknown_stereo(Some(&mut s), umode);
            got.push(format!("{u}"));
            got.extend(pr_stereo(&s));
        }
        assert_eq!(got, want, "case {ci}: {inp}");
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn create_inchi_matches_the_c_code() {
    use rdoxide::external::inchi::ichi::{
        InChI, InChIAux, InChIIsotopicAtom, InChIIsotopicTGroup, InChIStereo, OrigInfo,
    };
    use rdoxide::external::inchi::ichimake::{InpAtomData, create_inchi};
    use rdoxide::external::inchi::ichiparm::{InputParms, REQ_MODE_ISO};
    use rdoxide::external::inchi::ichitaut::TGroupInfo;
    use rdoxide::external::inchi::inpdef::InpAtom;
    fn alloc_stereo(n: usize, nb: usize) -> InChIStereo {
        InChIStereo {
            n_number: vec![0; n],
            t_parity: vec![0; n],
            n_number_inv: vec![0; n],
            t_parity_inv: vec![0; n],
            n_bond_atom1: vec![0; nb],
            n_bond_atom2: vec![0; nb],
            b_parity: vec![0; nb],
            ..Default::default()
        }
    }
    fn pr_stereo(tag: &str, s: &Option<InChIStereo>) -> String {
        let Some(s) = s else {
            return format!("{tag} -");
        };
        let mut l = format!(
            "{tag} {} {} {} {} |",
            s.n_number_of_stereo_centers,
            s.n_number_of_stereo_bonds,
            s.n_comp_inv2_abs,
            s.b_trivial_inv
        );
        for i in 0..s.n_number_of_stereo_centers as usize {
            l += &format!(
                " {},{},{},{}",
                s.n_number[i], s.t_parity[i], s.n_number_inv[i], s.t_parity_inv[i]
            );
        }
        l += " |";
        for i in 0..s.n_number_of_stereo_bonds as usize {
            l += &format!(
                " {},{},{}",
                s.n_bond_atom1[i], s.n_bond_atom2[i], s.b_parity[i]
            );
        }
        l
    }
    fn pr_arr(tag: &str, a: &[u16], n: usize) -> String {
        if a.is_empty() {
            return format!("{tag} -");
        }
        format!(
            "{tag}{}",
            a.iter()
                .take(n)
                .map(|x| format!(" {x}"))
                .collect::<String>()
        )
    }
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_create.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(inp.split_whitespace());
        let n = t.i() as usize;
        let which = t.i();
        let mode = t.i() as u64;
        let flags = t.i() as u64;
        let done = t.i() as u64;
        let loose = t.i();
        let mut at = vec![InpAtom::default(); n];
        for a in at.iter_mut() {
            let v = t.i() as usize;
            a.valence = v as i8;
            a.chem_bonds_valence = t.i() as i8;
            let el = t.i();
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
            a.radical = t.i() as i8;
            a.charge = t.i() as i8;
            a.num_h = t.i() as i8;
            a.num_iso_h = [t.i() as i8, t.i() as i8, t.i() as i8];
            a.iso_atw_diff = t.i() as i8;
            a.x = t.0.next().unwrap().parse().unwrap();
            a.y = t.0.next().unwrap().parse().unwrap();
            a.z = t.0.next().unwrap().parse().unwrap();
            for j in 0..v {
                a.neighbor[j] = t.i() as u16;
                a.bond_type[j] = t.i() as u8;
            }
        }
        for (i, a) in at.iter_mut().enumerate() {
            a.orig_at_number = i as u16 + 1;
        }
        let mut nd = [InpAtomData::default(), InpAtomData::default()];
        if which & 1 != 0 {
            nd[0].at = vec![InpAtom::default(); n + 8];
            nd[0].num_at = n as i32;
        }
        if which & 2 != 0 {
            nd[1].at = vec![InpAtom::default(); n + 8];
            nd[1].num_at = n as i32;
            if which & 4 != 0 {
                nd[1].at_fixed_bonds = Some(vec![InpAtom::default(); n + 8]);
            }
        }
        let nb: usize = at.iter().map(|a| a.valence as usize).sum::<usize>() / 2;
        let niso = at
            .iter()
            .filter(|a| a.iso_atw_diff != 0 || a.num_iso_h.iter().any(|&x| x != 0))
            .count();
        let b_iso = mode & REQ_MODE_ISO != 0;
        let mut pi: [Option<InChI>; 2] = [None, None];
        let mut pa: [Option<InChIAux>; 2] = [None, None];
        for k in 0..2 {
            if which & (1 << k) == 0 {
                continue;
            }
            let mut p = InChI {
                n_atom: vec![0; n],
                n_conn_table: vec![0; n + nb],
                n_tautomer: vec![0; 5 * n / 2 + 1],
                n_num_h: vec![0; n],
                n_num_h_fixed: Some(vec![0; n]),
                stereo: Some(alloc_stereo(n, nb)),
                ..Default::default()
            };
            let al = n + n / 2;
            let mut a = InChIAux {
                n_orig_at_nos_in_canon_ord: vec![0; al],
                n_orig_at_nos_in_canon_ord_inv: vec![0; al],
                n_constit_equ_numbers: vec![0; al],
                orig_info: vec![OrigInfo::default(); n],
                ..Default::default()
            };
            if n > 1 {
                a.n_constit_equ_t_group_numbers = vec![0; n / 2 + 1];
            }
            if b_iso {
                if niso > 0 {
                    p.isotopic_atom = vec![InChIIsotopicAtom::default(); niso];
                    p.isotopic_t_group = vec![InChIIsotopicTGroup::default(); niso];
                }
                p.n_possible_locations_of_isotopic_h = vec![0; n + 1];
                p.stereo_isotopic = Some(alloc_stereo(n, nb));
                a.n_isotopic_orig_at_nos_in_canon_ord = vec![0; al];
                a.n_isotopic_orig_at_nos_in_canon_ord_inv = vec![0; al];
                a.n_constit_equ_isotopic_numbers = vec![0; al];
                a.n_constit_equ_isotopic_t_group_numbers = vec![0; n / 2 + 1];
            }
            pi[k] = Some(p);
            pa[k] = Some(a);
        }
        let ip = InputParms {
            b_loose_tsa_check: loose as i32,
            ..Default::default()
        };
        let mut ti = TGroupInfo::default();
        let mut se = String::new();
        let (mut fl, mut dn) = (flags, done);
        let ret = create_inchi(
            &mut pi,
            &mut pa,
            &at,
            &mut nd,
            n,
            mode,
            &mut fl,
            &mut dn,
            Some(&mut ti),
            &mut se,
            &ip,
        );
        let mut got = vec![
            format!(
                "R {ret} {fl} {dn} {} {}",
                ti.num_t_groups, ti.n_num_endpoints
            ),
            format!("E {se}"),
        ];
        for k in 0..2 {
            let d = &nd[k];
            got.push(format!(
                "D{k} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                d.num_at,
                d.num_removed_h,
                d.b_deleted,
                d.b_tautomeric,
                d.b_taut_preprocessed,
                d.n_num_removed_protons,
                d.n_num_removed_protons_isotopic[0],
                d.n_num_removed_protons_isotopic[1],
                d.n_num_removed_protons_isotopic[2],
                d.num_iso_h[0],
                i32::from(d.num_iso_h[1]) * 100 + i32::from(d.num_iso_h[2]),
                d.b_taut_flags,
                d.b_taut_flags_done,
                d.b_normalization_flags
            ));
            if !d.at.is_empty() {
                let mut l = format!("DA{k}");
                for a in &d.at[..n] {
                    l += &format!(
                        " {},{},{},{},{},{},{},{},{},{},{},{}:",
                        a.charge,
                        a.num_h,
                        a.radical,
                        a.chem_bonds_valence,
                        a.endpoint,
                        a.c_point,
                        a.at_type,
                        a.num_iso_h[0],
                        a.num_iso_h[1],
                        a.num_iso_h[2],
                        a.valence,
                        a.orig_at_number
                    );
                    for j in 0..a.valence as usize {
                        l += &format!("{}.{}.", a.neighbor[j], a.bond_type[j]);
                    }
                }
                got.push(l);
            }
            let (Some(p), Some(a)) = (&pi[k], &pa[k]) else {
                continue;
            };
            let na = p.n_number_of_atoms as usize;
            got.push(format!(
                "I{k} {} {} {} {} {}",
                p.n_error_code,
                p.n_flags,
                p.n_total_charge,
                p.n_number_of_atoms,
                p.sz_hill_formula.as_deref().unwrap_or("NULL")
            ));
            got.push(format!(
                "IA{k}{}",
                p.n_atom[..na]
                    .iter()
                    .map(|x| format!(" {x}"))
                    .collect::<String>()
            ));
            got.push(format!(
                "IC{k} {}{}",
                p.len_conn_table,
                p.n_conn_table[..p.len_conn_table as usize]
                    .iter()
                    .map(|x| format!(" {x}"))
                    .collect::<String>()
            ));
            got.push(format!(
                "IT{k} {}{}",
                p.len_tautomer,
                p.n_tautomer[..p.len_tautomer as usize]
                    .iter()
                    .map(|x| format!(" {x}"))
                    .collect::<String>()
            ));
            let fixed = p.n_num_h_fixed.as_ref();
            got.push(format!(
                "IH{k}{}",
                (0..na)
                    .map(|i| format!(
                        " {}/{}",
                        p.n_num_h[i],
                        fixed.map_or(-1, |f| i32::from(f[i]))
                    ))
                    .collect::<String>()
            ));
            let mut l = format!(
                "II{k} {} {}",
                p.n_number_of_isotopic_atoms, p.n_number_of_isotopic_t_groups
            );
            for e in p
                .isotopic_atom
                .iter()
                .take(p.n_number_of_isotopic_atoms.max(0) as usize)
            {
                l += &format!(
                    " {},{},{},{},{}",
                    e.n_atom_number, e.n_iso_difference, e.n_num_h, e.n_num_d, e.n_num_t
                );
            }
            l += " |";
            for e in p
                .isotopic_t_group
                .iter()
                .take(p.n_number_of_isotopic_t_groups.max(0) as usize)
            {
                l += &format!(
                    " {},{},{},{}",
                    e.n_t_group_number, e.n_num_h, e.n_num_d, e.n_num_t
                );
            }
            got.push(l);
            let loc_n = p
                .n_possible_locations_of_isotopic_h
                .first()
                .map_or(0, |&x| x as usize);
            got.push(pr_arr("IL", &p.n_possible_locations_of_isotopic_h, loc_n));
            got.push(pr_stereo("IS", &p.stereo));
            got.push(pr_stereo("IZ", &p.stereo_isotopic));
            got.push(format!(
                "X{k} {} {} {} {} {}",
                a.n_error_code,
                a.n_number_of_atoms,
                a.n_number_of_t_groups,
                a.b_is_tautomeric,
                a.n_canon_flags
            ));
            got.push(pr_arr("XA", &a.n_orig_at_nos_in_canon_ord, n));
            got.push(pr_arr("XB", &a.n_orig_at_nos_in_canon_ord_inv, n));
            got.push(pr_arr("XC", &a.n_isotopic_orig_at_nos_in_canon_ord, n));
            got.push(pr_arr("XD", &a.n_isotopic_orig_at_nos_in_canon_ord_inv, n));
            got.push(pr_arr("XE", &a.n_constit_equ_numbers, n));
            got.push(pr_arr("XF", &a.n_constit_equ_isotopic_numbers, n));
            got.push(pr_arr(
                "XG",
                &a.n_constit_equ_t_group_numbers,
                a.n_number_of_t_groups.max(0) as usize,
            ));
            got.push(pr_arr(
                "XH",
                &a.n_constit_equ_isotopic_t_group_numbers,
                a.n_number_of_t_groups.max(0) as usize,
            ));
            got.push(format!(
                "XO{}",
                a.orig_info
                    .iter()
                    .take(n)
                    .map(|o| format!(" {},{},{}", o.c_charge, o.c_radical, o.c_unusual_valence))
                    .collect::<String>()
            ));
            got.push(format!(
                "XP {} {} {} {}",
                a.n_num_removed_protons,
                a.n_num_removed_isotopic_h[0],
                a.n_num_removed_isotopic_h[1],
                a.n_num_removed_isotopic_h[2]
            ));
        }
        assert_eq!(got, want, "case {ci}");
    }
}

#[test]
fn reversed_inchi_comparison_matches_the_c_code() {
    use rdoxide::external::inchi::ichi::InChIAux;
    use rdoxide::external::inchi::ichimake::{
        Icr, compare_icr, compare_reversed_inchi, compare_reversed_inchi2,
        compare_reversed_stereo_inchi, compare_reversed_stereo_inchi2,
    };
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_cmprev.json")).unwrap();
    let aux = |t: &mut InchiTok| {
        (t.i() != 0).then(|| InChIAux {
            n_num_removed_protons: t.i() as i16,
            n_num_removed_isotopic_h: [t.i() as i16, t.i() as i16, t.i() as i16],
            ..InChIAux::default()
        })
    };
    let dump = |nm: &str, v: &[u16]| {
        let mut s = format!(" {nm}");
        for x in v {
            s += &format!(",{x}");
        }
        s
    };
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let (c1, c2) = (t.inchi(), t.inchi());
        let (a1, a2) = (aux(&mut t), aux(&mut t));
        let r = compare_reversed_inchi(c1.as_ref(), c2.as_ref(), a1.as_ref(), a2.as_ref());
        assert_eq!(r.to_string(), want[0], "case {ci}");
        let st = |c: &Option<InChI>| c.as_ref().and_then(|c| c.stereo.clone());
        let si = |c: &Option<InChI>| c.as_ref().and_then(|c| c.stereo_isotopic.clone());
        let (s1, s2, t1, t2) = (st(&c1), st(&c2), si(&c1), si(&c2));
        let l = format!(
            "{} {} {} {}",
            compare_reversed_stereo_inchi(s1.as_ref(), s2.as_ref()),
            compare_reversed_stereo_inchi(t1.as_ref(), t2.as_ref()),
            compare_reversed_stereo_inchi(s1.as_ref(), t2.as_ref()),
            compare_reversed_stereo_inchi(t1.as_ref(), s2.as_ref())
        );
        assert_eq!(l, want[1], "case {ci}");
        let mut icr = Icr::default();
        let r1 = compare_reversed_stereo_inchi2(s1.as_ref(), s2.as_ref(), &mut icr);
        let r2 = compare_reversed_stereo_inchi2(t1.as_ref(), t2.as_ref(), &mut icr);
        let l = format!(
            "{r1} {r2}{}{}{}{}{}{}{}{}",
            dump("sc1", &icr.sc_in1_only),
            dump("sc2", &icr.sc_in2_only),
            dump("sb1", &icr.sb_in1_only),
            dump("sb2", &icr.sb_in2_only),
            dump("sbu1", &icr.sb_undef_in1_only),
            dump("sbu2", &icr.sb_undef_in2_only),
            dump("scu1", &icr.sc_undef_in1_only),
            dump("scu2", &icr.sc_undef_in2_only),
        );
        assert_eq!(l, want[2], "case {ci}");

        let mut q = Icr::default();
        let f = compare_reversed_inchi2(c1.as_ref(), c2.as_ref(), a1.as_ref(), a2.as_ref(), &mut q);
        let dn = |nm: &str, v: &[(u16, i8)]| {
            let mut s = format!(" {nm}");
            for (a, b) in v {
                s += &format!(",{a}:{b}");
            }
            s
        };
        let l = format!(
            "{f} {} {} {} {} {} {}{}{}{}{}{}{}{}{}{}{}{}{}{}",
            q.tot_num_h1,
            q.tot_num_h2,
            q.num_taut_h1,
            q.num_taut_h2,
            q.num_taut_m1,
            q.num_taut_m2,
            dump("e1", &q.endp_in1_only),
            dump("e2", &q.endp_in2_only),
            dn("dp", &q.diff_pos_h),
            dn("f1", &q.fixed_h1_more),
            dn("f2", &q.fixed_h2_more),
            dump("sc1", &q.sc_in1_only),
            dump("sc2", &q.sc_in2_only),
            dump("sb1", &q.sb_in1_only),
            dump("sb2", &q.sb_in2_only),
            dump("sbu1", &q.sb_undef_in1_only),
            dump("sbu2", &q.sb_undef_in2_only),
            dump("scu1", &q.sc_undef_in1_only),
            dump("scu2", &q.sc_undef_in2_only),
        );
        assert_eq!(l, want[3], "case {ci}");
        let q2 = Icr {
            flags: t.i() as u64,
            ..Icr::default()
        };
        let (r, x1, x2) = compare_icr(&q, &q2, t.i() as u64);
        assert_eq!(format!("{r} {x1} {x2}"), want[4], "case {ci}");
    }
}

#[test]
fn create_one_component_inchi_matches_the_c_code() {
    use rdoxide::external::inchi::ichi::{InChI, InChIAux};
    use rdoxide::external::inchi::ichimake::InpAtomData;
    use rdoxide::external::inchi::ichiparm::InputParms;
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::runichi::{NormCanonFlags, create_one_component_inchi};
    use rdoxide::external::inchi::runichi3::StructData;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_onecomp.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let sd_done = t.i() as u64;
        let nowarn = t.i() as i32;
        let ii = t.i() as usize;
        let outopt = t.i() as i32;
        let _ = t.i();
        let n = t.i() as usize;
        let _which = t.i();
        let mode = t.i() as u64;
        let flags = t.i() as u64;
        let done = t.i() as u64;
        let loose = t.i();
        let mut at = vec![InpAtom::default(); n];
        for (i, a) in at.iter_mut().enumerate() {
            let v = t.i() as usize;
            a.valence = v as i8;
            a.chem_bonds_valence = t.i() as i8;
            let el = t.i();
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
            a.radical = t.i() as i8;
            a.charge = t.i() as i8;
            a.num_h = t.i() as i8;
            a.num_iso_h = [t.i() as i8, t.i() as i8, t.i() as i8];
            a.iso_atw_diff = t.i() as i8;
            a.x = t.0.next().unwrap().parse().unwrap();
            a.y = t.0.next().unwrap().parse().unwrap();
            a.z = t.0.next().unwrap().parse().unwrap();
            for j in 0..v {
                a.neighbor[j] = t.i() as u16;
                a.bond_type[j] = t.i() as u8;
            }
            a.orig_at_number = i as u16 + 1;
        }
        let ip = InputParms {
            b_loose_tsa_check: loose as i32,
            n_mode: mode,
            b_taut_flags: flags,
            b_taut_flags_done: done,
            b_no_warnings: nowarn,
            b_inchi_output_options: outopt,
            ..Default::default()
        };
        let mut sd = StructData::default();
        sd.taut_flags_done[0] = sd_done;
        let mut cur = InpAtomData {
            at,
            num_at: n as i32,
            ..Default::default()
        };
        let mut nd = [InpAtomData::default(), InpAtomData::default()];
        let mut pi: [Option<InChI>; 2] = [None, None];
        let mut pa: [Option<InChIAux>; 2] = [None, None];
        let mut nc = NormCanonFlags::default();
        let ret = create_one_component_inchi(
            &mut sd, &ip, &mut cur, &mut pi, &mut pa, ii, 0, &mut nd, &mut nc,
        );
        let mut got = vec![
            format!(
                "R {ret} {} {} {} {} {}",
                sd.n_error_code,
                sd.num_non_taut[ii],
                sd.num_taut[ii],
                cur.num_bonds,
                cur.num_isotopic
            ),
            format!("E {}", sd.str_err_struct),
            format!(
                "C {}",
                if n > 0 {
                    i64::from(cur.at[0].component)
                } else {
                    -1
                }
            ),
        ];
        for k in 0..2 {
            got.push(format!(
                "B{k} {} {} {} {}",
                nd[k].b_exists,
                nd[k].b_has_isotopic_layer,
                nd[k].num_at,
                i32::from(!nd[k].at.is_empty())
            ));
            got.push(format!(
                "F{k} {} {} {}",
                nc.b_taut_flags[ii][k],
                nc.b_taut_flags_done[ii][k],
                nc.b_normalization_flags[ii][k]
            ));
            let Some(p) = &pi[k] else {
                got.push(format!("I{k} -"));
                continue;
            };
            got.push(format!(
                "I{k} {} {} {} {} {} {} {} {}",
                p.n_error_code,
                p.n_number_of_atoms,
                p.n_total_charge,
                p.sz_hill_formula.as_deref().unwrap_or("NULL"),
                p.len_conn_table,
                p.len_tautomer,
                p.n_number_of_isotopic_atoms,
                p.n_number_of_isotopic_t_groups
            ));
            let a = pa[k].as_ref();
            got.push(format!(
                "A{k} {} {} {} {}",
                a.map_or(-1, |a| a.n_number_of_atoms),
                a.map_or(-1, |a| i32::from(a.n_num_removed_protons)),
                a.map_or(-1, |a| a
                    .n_num_removed_isotopic_h
                    .iter()
                    .map(|&x| i32::from(x))
                    .sum::<i32>()),
                a.map_or(-1, |a| a.n_canon_flags)
            ));
        }
        assert_eq!(got, want, "case {ci}");
    }
}

#[test]
fn create_one_structure_inchi_matches_the_c_code() {
    use rdoxide::external::inchi::ichiparm::InputParms;
    use rdoxide::external::inchi::inpdef::InpAtom;
    use rdoxide::external::inchi::runichi::{
        NormCanonFlags, StructResults, create_one_structure_inchi,
    };
    use rdoxide::external::inchi::runichi3::StructData;
    use rdoxide::external::inchi::strutil::OrigAtomData;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_onestruct.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let n = t.i() as usize;
        let nbonds = t.i() as i32;
        let tf = t.i() as u64;
        let fnu = t.i() as i32;
        let nowarn = t.i() as i32;
        let mode = t.i() as u64;
        let cmpc = t.i() as i32;
        let mut at = vec![InpAtom::default(); n];
        for (i, a) in at.iter_mut().enumerate() {
            let v = t.i() as usize;
            a.valence = v as i8;
            a.chem_bonds_valence = t.i() as i8;
            a.el_number = t.i() as u8;
            a.charge = t.i() as i8;
            a.radical = t.i() as i8;
            a.iso_atw_diff = t.i() as i8;
            a.num_h = t.i() as i8;
            a.num_iso_h = [t.i() as i8, t.i() as i8, t.i() as i8];
            a.x = t.0.next().unwrap().parse().unwrap();
            a.y = t.0.next().unwrap().parse().unwrap();
            a.z = t.0.next().unwrap().parse().unwrap();
            for (k, b) in t.s().bytes().enumerate() {
                a.elname[k] = b;
            }
            for j in 0..v {
                a.neighbor[j] = t.i() as u16;
                a.bond_type[j] = t.i() as u8;
                a.bond_stereo[j] = t.i() as i8;
            }
            a.orig_at_number = i as u16 + 1;
        }
        let ip = InputParms {
            n_mode: mode,
            b_taut_flags: tf,
            b_fix_non_uniform_draw: fnu,
            b_no_warnings: nowarn,
            b_compare_components: cmpc,
            ..Default::default()
        };
        let mut sd = StructData::default();
        let mut orig = OrigAtomData {
            num_inp_atoms: n,
            num_inp_bonds: nbonds,
            at,
            ..Default::default()
        };
        let mut prep = [OrigAtomData::default(), OrigAtomData::default()];
        let mut res = StructResults::default();
        let mut nc = NormCanonFlags::default();
        let mut got: Vec<String> = Vec::new();
        for ii in 0..2 {
            let mut log = rdoxide::external::inchi::ichi_io::InchiIostream::default();
            let ret = create_one_structure_inchi(
                &mut sd, &ip, &mut log, 1, ii, &mut orig, &mut prep, &mut res, &mut nc,
            );
            got.push(format!(
                "S{ii} {ret} {} {} {} [{}] {} {} {} {} {}",
                sd.n_error_code,
                sd.n_error_type,
                sd.n_struct_read_error,
                sd.str_err_struct,
                sd.num_non_taut[ii],
                sd.num_taut[ii],
                sd.num_components[ii],
                prep[0].num_components,
                prep[1].num_components
            ));
            for j in 0..sd.num_components[ii] as usize {
                for k in 0..2 {
                    let Some(p) = &res.p_inchi[ii][j][k] else {
                        got.push(format!("I{ii}.{j}.{k} -"));
                        continue;
                    };
                    let a = &res.p_aux[ii][j][k];
                    let mut l = format!(
                        "I{ii}.{j}.{k} {} {} {} {} {} {} {} {} |",
                        p.n_error_code,
                        p.n_number_of_atoms,
                        p.n_total_charge,
                        p.sz_hill_formula.as_deref().unwrap_or("NULL"),
                        p.len_conn_table,
                        p.len_tautomer,
                        p.n_number_of_isotopic_atoms,
                        p.n_number_of_isotopic_t_groups
                    );
                    for q in 0..p.len_conn_table as usize {
                        l += &format!(" {}", p.n_conn_table[q]);
                    }
                    l += " |";
                    for q in 0..p.len_tautomer as usize {
                        l += &format!(" {}", p.n_tautomer[q]);
                    }
                    l += " |";
                    for q in 0..p.n_number_of_atoms as usize {
                        l += &format!(" {}", p.n_num_h[q]);
                    }
                    l += &format!(
                        " | {} {}",
                        a.as_ref().map_or(-1, |a| a.n_number_of_atoms),
                        a.as_ref()
                            .map_or(-1, |a| i32::from(a.n_num_removed_protons))
                    );
                    got.push(l);
                }
            }
        }
        let mut f = String::from("F");
        for ii in 0..2 {
            for k in 0..2 {
                f += &format!(
                    " {}/{}/{}",
                    nc.b_taut_flags[ii][k],
                    nc.b_taut_flags_done[ii][k],
                    nc.b_normalization_flags[ii][k]
                );
            }
        }
        got.push(f);
        assert_eq!(got, want, "case {ci}");
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn ichiprt2_strings_match_the_c_code() {
    use rdoxide::external::inchi::ichi::{
        InChI, InChIAux, InChIIsotopicAtom, InChIIsotopicTGroup, OrigInfo,
    };
    use rdoxide::external::inchi::ichi_io::{InchiIosString, inchi_strbuf_init};
    use rdoxide::external::inchi::ichiprt2::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt2.json")).unwrap();
    fn opt_u16(t: &mut InchiTok) -> Vec<u16> {
        if t.i() == 0 {
            return Vec::new();
        }
        let l = t.i();
        (0..l).map(|_| t.i() as u16).collect()
    }
    fn aux(t: &mut InchiTok) -> Option<InChIAux> {
        if t.i() == 0 {
            return None;
        }
        let mut a = InChIAux {
            n_number_of_atoms: t.i() as i32,
            n_number_of_t_groups: t.i() as i32,
            b_deleted: t.i() as i32,
            b_is_isotopic: t.i() as i32,
            ..Default::default()
        };
        a.n_constit_equ_numbers = opt_u16(t);
        a.n_constit_equ_t_group_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_t_group_numbers = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        if t.i() != 0 {
            let l = t.i();
            a.orig_info = (0..l)
                .map(|_| OrigInfo {
                    c_charge: t.i() as i8,
                    c_radical: t.i() as i8,
                    c_unusual_valence: t.i() as i8,
                })
                .collect();
        }
        Some(a)
    }
    fn iso_inchi(t: &mut InchiTok) -> Option<InChI> {
        if t.i() == 0 {
            return None;
        }
        let mut c = InChI {
            b_deleted: t.i() as i32,
            ..Default::default()
        };
        let l = t.i();
        c.n_number_of_isotopic_atoms = l as i32;
        if t.i() != 0 {
            for _ in 0..l {
                c.isotopic_atom.push(InChIIsotopicAtom {
                    n_atom_number: t.i() as u16,
                    n_iso_difference: t.i() as i16,
                    n_num_h: t.i() as i16,
                    n_num_d: t.i() as i16,
                    n_num_t: t.i() as i16,
                });
            }
        }
        let l = t.i();
        c.n_number_of_isotopic_t_groups = l as i32;
        if t.i() != 0 {
            for _ in 0..l {
                c.isotopic_t_group.push(InChIIsotopicTGroup {
                    n_t_group_number: t.i() as u16,
                    n_num_h: t.i() as u16,
                    n_num_d: t.i() as u16,
                    n_num_t: t.i() as u16,
                });
            }
        }
        Some(c)
    }
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let cmd = t.s();
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, 16, 16);
        let mut ovf = 0i32;
        let mut ret = 0i32;
        let text = |sb: &InchiIosString| {
            let p = sb.p_str.as_ref().unwrap();
            let n = p.iter().position(|&c| c == 0).unwrap_or(p.len());
            String::from_utf8_lossy(&p[..n]).into_owned()
        };
        let got: Vec<String> = match cmd.as_str() {
            "abc" | "dec" => {
                let hd = t.i();
                let d = t.s();
                let v = t.i() as i32;
                let bl = t.i() as i32;
                let delim = (hd != 0).then_some(d.as_str());
                let mut buf = [0x7fu8; 64];
                let r = if cmd == "abc" {
                    make_abc_number(&mut buf, bl, delim, v)
                } else {
                    make_dec_number(&mut buf, bl, delim, v)
                };
                let s = if r >= 0 {
                    let n = buf.iter().position(|&c| c == 0).unwrap_or(64);
                    String::from_utf8_lossy(&buf[..n]).into_owned()
                } else {
                    String::new()
                };
                vec![format!("{r} [{s}]")]
            }
            "eqlaux" => {
                let (a1, a2) = (aux(&mut t), aux(&mut t));
                let (e1, e2) = (t.i() as i32, t.i() as i32);
                let oi = a1
                    .as_ref()
                    .filter(|a| !a.orig_info.is_empty())
                    .is_some_and(|a| has_orig_info(&a.orig_info, a.n_number_of_atoms));
                vec![
                    format!(
                        "{} {} {}",
                        i32::from(eql_inchi_aux_equ(a1.as_ref(), e1, a2.as_ref(), e2)),
                        i32::from(eql_inchi_aux_num(a1.as_ref(), e1 & 3, a2.as_ref(), e2 & 3)),
                        i32::from(eql_orig_info(a1.as_ref(), a2.as_ref()))
                    ),
                    format!("{}", i32::from(oi)),
                ]
            }
            "eqliso" => {
                let (i1, i2) = (iso_inchi(&mut t), iso_inchi(&mut t));
                vec![format!(
                    "{}",
                    i32::from(eql_inchi_isotopic(i1.as_ref(), i2.as_ref()))
                )]
            }
            "hasequ" => {
                let l = t.i();
                let a: Vec<u16> = (0..l).map(|_| t.i() as u16).collect();
                vec![format!("{}", i32::from(has_equ_string(&a, l as i32)))]
            }
            _ => {
                match cmd.as_str() {
                    "mult" => {
                        let m = t.i() as i32;
                        let d = t.s();
                        let mode = t.i() as i32;
                        ovf = t.i() as i32;
                        ret = make_mult(m, &d, &mut sb, mode, &mut ovf);
                    }
                    "delim" => {
                        let has = t.i();
                        let d = t.s();
                        ovf = t.i() as i32;
                        ret = make_delim((has != 0).then_some(d.as_str()), &mut sb, &mut ovf);
                    }
                    "eqstr" => {
                        let has = t.i();
                        let d = t.s();
                        let m = t.i() as i32;
                        ovf = t.i() as i32;
                        ret = make_eq_str((has != 0).then_some(d.as_str()), m, &mut sb, &mut ovf);
                    }
                    "ct" => {
                        let l = t.i();
                        let ct: Vec<u16> = (0..l).map(|_| t.i() as u16).collect();
                        let na = t.i() as usize;
                        let has_h = t.i();
                        let h: Vec<i8> = (0..na).map(|_| t.i() as i8).collect();
                        let ad = t.i() != 0;
                        let mode = t.i() as i32;
                        ret = make_ct_string(
                            &ct,
                            l as i32,
                            ad,
                            (has_h != 0).then_some(h.as_slice()),
                            na,
                            &mut sb,
                            mode,
                            &mut ovf,
                        );
                    }
                    "h" => {
                        let ad = t.i() != 0;
                        let l = t.i();
                        let a: Vec<i8> = (0..l).map(|_| t.i() as i8).collect();
                        let mode = t.i() as i32;
                        ret = make_h_string(ad, &a, l as i32, &mut sb, mode, &mut ovf);
                    }
                    "taut" => {
                        let l = t.i();
                        let a: Vec<u16> = (0..l).map(|_| t.i() as u16).collect();
                        let ad = t.i() != 0;
                        let mode = t.i() as i32;
                        ret = make_taut_string(&a, l as i32, ad, &mut sb, mode, &mut ovf);
                    }
                    "crv" => {
                        let l = t.i();
                        let a: Vec<OrigInfo> = (0..l)
                            .map(|_| OrigInfo {
                                c_charge: t.i() as i8,
                                c_radical: t.i() as i8,
                                c_unusual_valence: t.i() as i8,
                            })
                            .collect();
                        let ad = t.i() != 0;
                        let mode = t.i() as i32;
                        ret = make_crv_string(&a, l as i32, ad, &mut sb, mode, &mut ovf);
                    }
                    "equ" => {
                        let l = t.i();
                        let a: Vec<u16> = (0..l).map(|_| t.i() as u16).collect();
                        let ad = t.i() != 0;
                        let mode = t.i() as i32;
                        ret = make_equ_string(&a, l as i32, ad, &mut sb, mode, &mut ovf);
                    }
                    "isoatom" => {
                        let l = t.i();
                        let a: Vec<InChIIsotopicAtom> = (0..l)
                            .map(|_| InChIIsotopicAtom {
                                n_atom_number: t.i() as u16,
                                n_iso_difference: t.i() as i16,
                                n_num_h: t.i() as i16,
                                n_num_d: t.i() as i16,
                                n_num_t: t.i() as i16,
                            })
                            .collect();
                        let mode = t.i() as i32;
                        ret = make_iso_atom_string(&a, l as i32, &mut sb, mode, &mut ovf);
                    }
                    "isotaut" => {
                        let l = t.i();
                        let a: Vec<InChIIsotopicTGroup> = (0..l)
                            .map(|_| InChIIsotopicTGroup {
                                n_t_group_number: t.i() as u16,
                                n_num_h: t.i() as u16,
                                n_num_d: t.i() as u16,
                                n_num_t: t.i() as u16,
                            })
                            .collect();
                        let mode = t.i() as i32;
                        ret = make_iso_taut_string(&a, l as i32, &mut sb, mode, &mut ovf);
                    }
                    "isoh" => {
                        let h = [t.i() as i32, t.i() as i32, t.i() as i32];
                        let mode = t.i() as i32;
                        ret = make_iso_h_string(&h, &mut sb, mode, &mut ovf);
                    }
                    "stereo" => {
                        let l = t.i() as i32;
                        let a1 = opt_u16(&mut t);
                        let a2 = opt_u16(&mut t);
                        let p: Vec<i8> = opt_u16(&mut t).into_iter().map(|x| x as i8).collect();
                        let _ad = t.i();
                        let mode = t.i() as i32;
                        ret = make_stereo_string(
                            (!a1.is_empty()).then_some(a1.as_slice()),
                            (!a2.is_empty()).then_some(a2.as_slice()),
                            (!p.is_empty()).then_some(p.as_slice()),
                            l,
                            &mut sb,
                            mode,
                            &mut ovf,
                        );
                    }
                    "seq" => {
                        let l = t.i();
                        let a: Vec<i32> = (0..l).map(|_| t.i() as i32).collect();
                        print_sequence_of_nums_compressing_ranges(&a, &mut sb);
                        ret = 0;
                    }
                    other => panic!("unknown command {other}"),
                }
                vec![format!("{ret} {ovf} [{}]", text(&sb))]
            }
        };
        let _ = &mut ret;
        assert_eq!(got, want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt3_first_substrings_match_the_c_code() {
    use rdoxide::external::inchi::ichi::InChISort;
    use rdoxide::external::inchi::ichi_io::{InchiIosString, inchi_strbuf_init};
    use rdoxide::external::inchi::ichiprt3::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt3a.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let cmd = t.s();
        let n = t.i() as usize;
        let ot = t.i() as i32;
        let mult = t.i() != 0;
        let sorts = |t: &mut InchiTok| -> Vec<InChISort> {
            (0..n)
                .map(|i| InChISort {
                    p_inchi: [t.inchi(), t.inchi()],
                    ord_number: i as i16,
                    ..Default::default()
                })
                .collect()
        };
        let a = sorts(&mut t);
        let b = sorts(&mut t);
        let am = t.i() as i32;
        let tm = t.i() as i32;
        let second = t.i() != 0;
        let omit = t.i() != 0;
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, 16, 16);
        let mut ovf = 0;
        let ret = match cmd.as_str() {
            "hill" => str_hill_formula(&a, &mut sb, &mut ovf, ot, n, mult),
            "hill2" => str_hill_formula2(&a, &b, &mut sb, &mut ovf, ot, n, mult),
            "conn" => str_connections(&a, &mut sb, &mut ovf, ot, am, n, mult),
            "hatoms" => str_h_atoms(&a, &mut sb, &mut ovf, ot, am, tm, n, mult),
            "charge" => str_charge2(&a, &b, &mut sb, &mut ovf, ot, n, second, omit, mult),
            "fixedh" => str_fixed_h_atoms(&a, &mut sb, &mut ovf, ot, am, n, mult),
            other => panic!("{other}"),
        };
        let p = sb.p_str.as_ref().unwrap();
        let len = p.iter().position(|&c| c == 0).unwrap_or(p.len());
        let got = format!("{ret} {ovf} [{}]", String::from_utf8_lossy(&p[..len]));
        assert_eq!(vec![got], want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt3_stereo_and_isotopic_substrings_match_the_c_code() {
    use rdoxide::external::inchi::ichi::InChISort;
    use rdoxide::external::inchi::ichi_io::{InchiIosString, inchi_strbuf_init};
    use rdoxide::external::inchi::ichiprt3::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt3b.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let cmd = t.s();
        let n = t.i() as usize;
        let ot = t.i() as i32;
        let mult = t.i() != 0;
        let sorts = |t: &mut InchiTok| -> Vec<InChISort> {
            (0..n)
                .map(|i| InChISort {
                    p_inchi: [t.inchi(), t.inchi()],
                    ord_number: i as i16,
                    ..Default::default()
                })
                .collect()
        };
        let a = sorts(&mut t);
        let b = sorts(&mut t);
        let tm = t.i() as i32;
        let second = t.i() != 0;
        let omit = t.i() != 0;
        let abc = t.i() != 0;
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, 16, 16);
        let mut ovf = 0;
        let ret = match cmd.as_str() {
            "sp2" => str_sp2(&a, &b, &mut sb, &mut ovf, ot, tm, n, second, omit, mult),
            "sp3" => str_sp3(&a, &b, &mut sb, &mut ovf, ot, tm, n, second, omit, mult),
            "abs" => str_stereo_abs_inv(&a, &mut sb, &mut ovf, ot, n),
            "isoatoms" => str_iso_atoms(
                &a, &b, &mut sb, &mut ovf, ot, tm, n, abc, second, omit, mult,
            ),
            "isosp2" => str_iso_sp2(&a, &b, &mut sb, &mut ovf, ot, tm, n, second, omit, mult),
            "isosp3" => str_iso_sp3(&a, &b, &mut sb, &mut ovf, ot, tm, n, second, omit, mult),
            "isoabs" => str_iso_stereo_abs_inv(&a, &mut sb, &mut ovf, ot, n),
            other => panic!("{other}"),
        };
        let p = sb.p_str.as_ref().unwrap();
        let len = p.iter().position(|&c| c == 0).unwrap_or(p.len());
        let got = format!("{ret} {ovf} [{}]", String::from_utf8_lossy(&p[..len]));
        assert_eq!(vec![got], want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt3_aux_equivalence_and_inversion_match_the_c_code() {
    use rdoxide::external::inchi::ichi::{InChIAux, InChISort, OrigInfo};
    use rdoxide::external::inchi::ichi_io::{InchiIosString, inchi_strbuf_init};
    use rdoxide::external::inchi::ichiprt3::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt3c.json")).unwrap();
    fn opt_u16(t: &mut InchiTok) -> Vec<u16> {
        if t.i() == 0 {
            return Vec::new();
        }
        let l = t.i();
        (0..l).map(|_| t.i() as u16).collect()
    }
    fn aux(t: &mut InchiTok) -> Option<InChIAux> {
        if t.i() == 0 {
            return None;
        }
        let mut a = InChIAux {
            n_number_of_atoms: t.i() as i32,
            n_number_of_t_groups: t.i() as i32,
            b_deleted: t.i() as i32,
            b_is_isotopic: t.i() as i32,
            ..Default::default()
        };
        a.n_constit_equ_numbers = opt_u16(t);
        a.n_constit_equ_t_group_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_t_group_numbers = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        if t.i() != 0 {
            let l = t.i();
            a.orig_info = (0..l)
                .map(|_| OrigInfo {
                    c_charge: t.i() as i8,
                    c_radical: t.i() as i8,
                    c_unusual_valence: t.i() as i8,
                })
                .collect();
        }
        Some(a)
    }
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let cmd = t.s();
        let n = t.i() as usize;
        let ot = t.i() as i32;
        let mult = t.i() != 0;
        let sorts = |t: &mut InchiTok| -> Vec<InChISort> {
            (0..n)
                .map(|i| InChISort {
                    p_inchi: [t.inchi(), t.inchi()],
                    p_inchi_aux: [aux(t), aux(t)],
                    ord_number: i as i16,
                })
                .collect()
        };
        let a = sorts(&mut t);
        let b = sorts(&mut t);
        let tm = t.i() as i32;
        let second = t.i() != 0;
        let omit = t.i() != 0;
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, 16, 16);
        let mut ovf = 0;
        let ret = match cmd.as_str() {
            "equ" => str_aux_equ(&a, &b, &mut sb, &mut ovf, ot, tm, n, second, omit, mult),
            "isoequ" => str_aux_iso_equ(&a, &b, &mut sb, &mut ovf, ot, tm, n, second, omit, mult),
            "tg" => str_aux_tgroup_equ(&a, &mut sb, &mut ovf, ot, tm, n, mult),
            "isotg" => str_aux_iso_tgroup_equ(&a, &mut sb, &mut ovf, ot, tm, n, omit, mult),
            "invsp3" => str_aux_inv_sp3(&a, &b, &mut sb, &mut ovf, ot, tm, n, second, omit, mult),
            "inviso" => {
                str_aux_inv_iso_sp3(&a, &b, &mut sb, &mut ovf, ot, tm, n, second, omit, mult)
            }
            other => panic!("{other}"),
        };
        let p = sb.p_str.as_ref().unwrap();
        let len = p.iter().position(|&c| c == 0).unwrap_or(p.len());
        let got = format!("{ret} {ovf} [{}]", String::from_utf8_lossy(&p[..len]));
        assert_eq!(vec![got], want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt3_aux_numbering_charge_and_transposition_match_the_c_code() {
    use rdoxide::external::inchi::ichi::{InChIAux, InChISort, OrigInfo};
    use rdoxide::external::inchi::ichi_io::{InchiIosString, inchi_strbuf_init};
    use rdoxide::external::inchi::ichiprt3::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt3d.json")).unwrap();
    fn opt_u16(t: &mut InchiTok) -> Vec<u16> {
        if t.i() == 0 {
            return Vec::new();
        }
        let l = t.i();
        (0..l).map(|_| t.i() as u16).collect()
    }
    fn aux(t: &mut InchiTok) -> Option<InChIAux> {
        if t.i() == 0 {
            return None;
        }
        let mut a = InChIAux {
            n_number_of_atoms: t.i() as i32,
            n_number_of_t_groups: t.i() as i32,
            b_deleted: t.i() as i32,
            b_is_isotopic: t.i() as i32,
            ..Default::default()
        };
        a.n_constit_equ_numbers = opt_u16(t);
        a.n_constit_equ_t_group_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_t_group_numbers = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        if t.i() != 0 {
            let l = t.i();
            a.orig_info = (0..l)
                .map(|_| OrigInfo {
                    c_charge: t.i() as i8,
                    c_radical: t.i() as i8,
                    c_unusual_valence: t.i() as i8,
                })
                .collect();
        }
        Some(a)
    }
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let cmd = t.s();
        let n = t.i() as usize;
        let ot = t.i() as i32;
        let mult = t.i() != 0;
        let sorts = |t: &mut InchiTok| -> Vec<InChISort> {
            (0..n)
                .map(|i| InChISort {
                    p_inchi: [t.inchi(), t.inchi()],
                    p_inchi_aux: [aux(t), aux(t)],
                    ord_number: i as i16,
                })
                .collect()
        };
        let mut a = sorts(&mut t);
        let mut b = sorts(&mut t);
        let tm = t.i() as i32;
        let second = t.i() != 0;
        let omit = t.i() != 0;
        for s in a.iter_mut() {
            s.ord_number = t.i() as i16;
        }
        for s in b.iter_mut() {
            s.ord_number = t.i() as i16;
        }
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, 16, 16);
        let mut ovf = 0;
        let ret = match cmd.as_str() {
            "numb" => str_aux_numb(&a, &mut sb, &mut ovf, ot, tm, n, second, omit),
            "invnumb" => str_aux_inv_sp3_numb(&a, &mut sb, &mut ovf, ot, tm, n, second, omit),
            "isonumb" => str_aux_iso_numb(&a, &mut sb, &mut ovf, ot, tm, n, second, omit),
            "invisonumb" => {
                str_aux_inv_iso_sp3_numb(&a, &mut sb, &mut ovf, ot, tm, n, second, omit)
            }
            "crv" => str_aux_charge_rad_val(&a, &mut sb, &mut ovf, ot, tm, n, mult),
            "trans" => {
                let tr = bin_aux_taut_trans(&a, &b, ot, n);
                str_aux_taut_trans(tr, &mut sb, &mut ovf, tm, n)
            }
            other => panic!("{other}"),
        };
        let p = sb.p_str.as_ref().unwrap();
        let len = p.iter().position(|&c| c == 0).unwrap_or(p.len());
        let got = format!("{ret} {ovf} [{}]", String::from_utf8_lossy(&p[..len]));
        assert_eq!(vec![got], want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt1_tags_and_first_layers_match_the_c_code() {
    use rdoxide::external::inchi::ichi::{
        DIFS_H_H_ATOMS, DIFS_P_PROTONS, DIFS_Q_CHARGE, InChISort,
    };
    use rdoxide::external::inchi::ichi_io::{
        INCHI_IOS_TYPE_STRING, InchiIosString, InchiIostream, IosFile, inchi_ios_init,
        inchi_strbuf_init, inchi_strbuf_printf,
    };
    use rdoxide::external::inchi::ichiprt1::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt1a.json")).unwrap();
    let text = |b: &InchiIosString| strbuf_str(b);
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let cmd = t.s();
        let n = t.i() as usize;
        let ot = t.i() as i32;
        let mult = t.i() as i32;
        let sorts = |t: &mut InchiTok| -> Vec<InChISort> {
            (0..n)
                .map(|i| InChISort {
                    p_inchi: [t.inchi(), t.inchi()],
                    ord_number: i as i16,
                    ..Default::default()
                })
                .collect()
        };
        let a = sorts(&mut t);
        let b = sorts(&mut t);
        let mut io = InchiOutCtl {
            atom_mode: t.i() as i32,
            taut_mode: t.i() as i32,
            n_tag: t.i() as i32,
            b_fh_tag: t.i() as i32,
            b_plain_text_tags: t.i() as i32,
            b_second_non_taut_pass: t.i() as i32,
            b_omit_repetitions: t.i() as i32,
            i_cur_taut_mode: t.i() as usize,
            n_cur_inchi_segment: t.i() as usize,
            n_num_removed_protons: t.i() as i32,
            n_pzz: t.i() as i32,
            n_zy: t.i() as i32,
            ..Default::default()
        };
        let nc2 = [t.i() as i32, t.i() as i32];
        let rec = t.i() as i32;
        let outopts = t.i() as i32;
        let beta = t.i() != 0;
        io.b_overflow = t.i() as i32;
        let seg = io.n_cur_inchi_segment;
        io.s_dif_segs[seg][DIFS_H_H_ATOMS] = t.i() as u8;
        io.s_dif_segs[seg][DIFS_Q_CHARGE] = t.i() as u8;
        io.s_dif_segs[seg][DIFS_P_PROTONS] = t.i() as u8;
        io.b_tag1 = t.i() as i32;
        io.b_use_multipliers = mult;
        io.b_out_type = ot;
        io.num_components = n;
        io.p_inchi_sort = a;
        io.p_inchi_sort2 = b;
        let (lf, tab) = set_line_separators(outopts);
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, 16, 16);
        let mut out = InchiIostream::default();
        inchi_ios_init(&mut out, INCHI_IOS_TYPE_STRING, IosFile::None);
        let got = match cmd.as_str() {
            "tag" => {
                let mut always = -7;
                let (tags, bits): (&[InchiTag], i32) = if io.b_tag1 & 0x40000 != 0 {
                    (&AUX_LBL, io.b_tag1 & 0x3ffff)
                } else {
                    (&IDENT_LBL, io.b_tag1 & 0x3ffff)
                };
                let tag = sz_get_tag(tags, io.n_tag, bits, &mut always);
                format!("{tag} {always}")
            }
            "opt" => {
                let (l1, l2) = get_save_opt_letters(io.b_tag1 as u8);
                format!(
                    "{l1}{l2} {} {} {outopts}",
                    lf.bytes().next().unwrap_or(0),
                    tab.bytes().next().unwrap_or(0)
                )
            }
            "lineend" => {
                if io.b_tag1 & 1 != 0 {
                    inchi_strbuf_printf(&mut sb, "x,y");
                }
                let tg = if io.b_tag1 & 2 != 0 {
                    "/c"
                } else if io.b_tag1 & 4 != 0 {
                    "/qq"
                } else {
                    ""
                };
                let ret = str_line_end(
                    tg,
                    &mut io.b_overflow,
                    &mut sb,
                    io.n_cur_inchi_segment as i32 - 2,
                    io.b_plain_text_tags,
                );
                format!(
                    "{ret} {} [{}] {}",
                    io.b_overflow,
                    text(&sb),
                    sb.n_used_length
                )
            }
            "ver" => {
                let ret = output_inchi_version_and_kind(&mut out, &mut sb, outopts, beta, lf);
                format!("{ret} [{}] [{}]", text(&out.s), text(&sb))
            }
            other => {
                let ret = match other {
                    "formula" => {
                        output_inchi_main_layer_formula(&mut out, &mut sb, nc2, rec, &mut io, lf)
                    }
                    "conn" => output_inchi_main_layer_connections(&mut out, &mut sb, &mut io, lf),
                    "hyd" => output_inchi_main_layer_hydrogens(&mut out, &mut sb, &mut io, lf),
                    "charge" => output_inchi_charge_and_removed_added_protons_layers(
                        &mut out, &mut sb, &mut io, lf,
                    ),
                    _ => panic!("{other}"),
                };
                format!(
                    "{ret} {} {} {} {} {} {} {} [{}] [{}] [{}]",
                    io.b_overflow,
                    io.tot_len,
                    io.tot_len2,
                    io.b_tag1,
                    io.b_always,
                    io.b_non_taut_non_iso_identifier_not_empty,
                    io.n_segm_action,
                    io.sz_tag1,
                    text(&out.s),
                    text(&sb)
                )
            }
        };
        let got: Vec<&str> = got.split('\n').filter(|l| !l.is_empty()).collect();
        assert_eq!(got, want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt1_stereo_and_isotopic_layers_match_the_c_code() {
    use rdoxide::external::inchi::ichi::{
        DIFL_F, DIFS_B_SBONDS, DIFS_H_H_ATOMS, DIFS_I_IATOMS, DIFS_M_SP3INV, DIFS_O_TRANSP,
        DIFS_S_STYPE, DIFS_T_SATOMS, InChIAux, InChISort, OrigInfo,
    };
    use rdoxide::external::inchi::ichi_io::{
        INCHI_IOS_TYPE_STRING, InchiIosString, InchiIostream, IosFile, inchi_ios_init,
        inchi_strbuf_init,
    };
    use rdoxide::external::inchi::ichiprt1::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt1b.json")).unwrap();
    fn opt_u16(t: &mut InchiTok) -> Vec<u16> {
        if t.i() == 0 {
            return Vec::new();
        }
        let l = t.i();
        (0..l).map(|_| t.i() as u16).collect()
    }
    fn aux(t: &mut InchiTok) -> Option<InChIAux> {
        if t.i() == 0 {
            return None;
        }
        let mut a = InChIAux {
            n_number_of_atoms: t.i() as i32,
            n_number_of_t_groups: t.i() as i32,
            b_deleted: t.i() as i32,
            b_is_isotopic: t.i() as i32,
            ..Default::default()
        };
        a.n_constit_equ_numbers = opt_u16(t);
        a.n_constit_equ_t_group_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_t_group_numbers = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        if t.i() != 0 {
            let l = t.i();
            a.orig_info = (0..l)
                .map(|_| OrigInfo {
                    c_charge: t.i() as i8,
                    c_radical: t.i() as i8,
                    c_unusual_valence: t.i() as i8,
                })
                .collect();
        }
        Some(a)
    }
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let cmd = t.s();
        let n = t.i() as usize;
        let ot = t.i() as i32;
        let mult = t.i() as i32;
        let sorts = |t: &mut InchiTok| -> Vec<InChISort> {
            (0..n)
                .map(|i| InChISort {
                    p_inchi: [t.inchi(), t.inchi()],
                    p_inchi_aux: [aux(t), aux(t)],
                    ord_number: i as i16,
                })
                .collect()
        };
        let mut a = sorts(&mut t);
        let mut b = sorts(&mut t);
        for s in a.iter_mut() {
            s.ord_number = t.i() as i16;
        }
        for s in b.iter_mut() {
            s.ord_number = t.i() as i16;
        }
        let mut io = InchiOutCtl {
            taut_mode: t.i() as i32,
            n_tag: t.i() as i32,
            b_fh_tag: t.i() as i32,
            b_plain_text_tags: t.i() as i32,
            b_second_non_taut_pass: t.i() as i32,
            b_omit_repetitions: t.i() as i32,
            i_cur_taut_mode: t.i() as usize,
            n_cur_inchi_segment: t.i() as usize,
            b_abc_numbers: t.i() as i32,
            b_output_type: t.i() as i32,
            ..Default::default()
        };
        for k in 0..2 {
            io.b_relative_stereo[k] = t.i() as i32;
        }
        for k in 0..2 {
            io.b_racemic_stereo[k] = t.i() as i32;
        }
        for k in 0..2 {
            io.b_isotopic_relative_stereo[k] = t.i() as i32;
        }
        for k in 0..2 {
            io.b_isotopic_racemic_stereo[k] = t.i() as i32;
        }
        for k in 0..3 {
            io.num_iso_h[k] = t.i() as i32;
        }
        let rec = t.i() as i32;
        io.b_overflow = t.i() as i32;
        let cs = io.n_cur_inchi_segment;
        io.s_dif_segs[cs][DIFS_B_SBONDS] = t.i() as u8;
        io.s_dif_segs[cs][DIFS_T_SATOMS] = t.i() as u8;
        io.s_dif_segs[cs][DIFS_M_SP3INV] = t.i() as u8;
        io.s_dif_segs[cs][DIFS_S_STYPE] = t.i() as u8;
        io.s_dif_segs[cs][DIFS_I_IATOMS] = t.i() as u8;
        io.s_dif_segs[cs][DIFS_H_H_ATOMS] = t.i() as u8;
        io.s_dif_segs[DIFL_F][DIFS_O_TRANSP] = t.i() as u8;
        io.b_use_multipliers = mult;
        io.b_out_type = ot;
        io.num_components = n;
        io.p_inchi_sort = a;
        io.p_inchi_sort2 = b;
        let (lf, _) = set_line_separators(if t.i() != 0 { 0x80 } else { 0 });
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, 16, 16);
        let mut out = InchiIostream::default();
        inchi_ios_init(&mut out, INCHI_IOS_TYPE_STRING, IosFile::None);
        let ret = if cmd == "stereo" {
            output_inchi_stereo_layer(&mut out, &mut sb, &mut io, lf)
        } else {
            output_inchi_isotopic_layer(&mut out, &mut sb, rec, &mut io, lf)
        };
        let got = format!(
            "{ret} {} {} {} {} {} {} {} {} {} {} {} {} [{}] [{}] [{}] [{}] [{}]",
            io.b_overflow,
            io.tot_len,
            io.tot_len2,
            io.b_tag1,
            io.b_tag2,
            io.b_tag3,
            io.b_always,
            io.b_non_taut_non_iso_identifier_not_empty,
            io.b_non_taut_iso_identifier_not_empty,
            io.n_segm_action,
            io.b_rel_rac,
            io.sort_print_inchi_flags,
            io.sz_tag1,
            io.sz_tag2,
            io.sz_tag3,
            strbuf_str(&out.s),
            strbuf_str(&sb)
        );
        let got: Vec<&str> = got.split('\n').filter(|l| !l.is_empty()).collect();
        assert_eq!(got, want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt1_fixed_h_layer_matches_the_c_code() {
    use rdoxide::external::inchi::ichi::{
        DIFL_F, DIFS_F_FORMULA, DIFS_H_H_ATOMS, InChIAux, InChISort, OrigInfo,
    };
    use rdoxide::external::inchi::ichi_io::{
        INCHI_IOS_TYPE_STRING, InchiIosString, InchiIostream, IosFile, inchi_ios_init,
        inchi_strbuf_init,
    };
    use rdoxide::external::inchi::ichiprt1::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt1c.json")).unwrap();
    fn opt_u16(t: &mut InchiTok) -> Vec<u16> {
        if t.i() == 0 {
            return Vec::new();
        }
        let l = t.i();
        (0..l).map(|_| t.i() as u16).collect()
    }
    fn aux(t: &mut InchiTok) -> Option<InChIAux> {
        if t.i() == 0 {
            return None;
        }
        let mut a = InChIAux {
            n_number_of_atoms: t.i() as i32,
            n_number_of_t_groups: t.i() as i32,
            b_deleted: t.i() as i32,
            b_is_isotopic: t.i() as i32,
            ..Default::default()
        };
        a.n_constit_equ_numbers = opt_u16(t);
        a.n_constit_equ_t_group_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_t_group_numbers = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        if t.i() != 0 {
            let l = t.i();
            a.orig_info = (0..l)
                .map(|_| OrigInfo {
                    c_charge: t.i() as i8,
                    c_radical: t.i() as i8,
                    c_unusual_valence: t.i() as i8,
                })
                .collect();
        }
        Some(a)
    }
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let n = t.i() as usize;
        let ot = t.i() as i32;
        let mult = t.i() as i32;
        let sorts = |t: &mut InchiTok| -> Vec<InChISort> {
            (0..n)
                .map(|i| InChISort {
                    p_inchi: [t.inchi(), t.inchi()],
                    p_inchi_aux: [aux(t), aux(t)],
                    ord_number: i as i16,
                })
                .collect()
        };
        let a = sorts(&mut t);
        let b = sorts(&mut t);
        let t0 = sorts(&mut t);
        let t1 = sorts(&mut t);
        let mut io = InchiOutCtl {
            atom_mode: t.i() as i32,
            taut_mode: t.i() as i32,
            n_tag: t.i() as i32,
            b_fh_tag: t.i() as i32,
            b_plain_text_tags: t.i() as i32,
            b_second_non_taut_pass: t.i() as i32,
            b_omit_repetitions: t.i() as i32,
            i_cur_taut_mode: t.i() as usize,
            n_cur_inchi_segment: t.i() as usize,
            b_tautomeric: t.i() as i32,
            b_non_tautomeric: t.i() as i32,
            b_non_taut_is_identical_to_taut: t.i() as i32,
            b_output_type: t.i() as i32,
            ..Default::default()
        };
        io.num_comp = [t.i() as i32, t.i() as i32];
        io.n_pzz = t.i() as i32;
        io.n_zy = t.i() as i32;
        let rec = t.i() as i32;
        io.b_overflow = t.i() as i32;
        io.sort_print_inchi_flags = t.i() as i32;
        io.b_non_taut_non_iso_identifier_not_empty = t.i() as i32;
        io.b_non_taut_iso_identifier_not_empty = t.i() as i32;
        io.s_dif_segs[DIFL_F][DIFS_F_FORMULA] = t.i() as u8;
        io.s_dif_segs[DIFL_F][DIFS_H_H_ATOMS] = t.i() as u8;
        io.b_use_multipliers = mult;
        io.b_out_type = ot;
        io.num_components = n;
        io.p_inchi_sort = a.clone();
        io.p_inchi_sort2 = b;
        io.p_inchi_sort_taut_and_non_taut = [t0.clone(), t1.clone()];
        let (lf, _) = set_line_separators(if t.i() != 0 { 0x80 } else { 0 });
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, 16, 16);
        let mut out = InchiIostream::default();
        inchi_ios_init(&mut out, INCHI_IOS_TYPE_STRING, IosFile::None);
        let mut repeat = -1;
        let switched = io.b_out_type == 3
            && io.b_non_taut_is_identical_to_taut == 0
            && rdoxide::external::inchi::ichi::inchi_segment_action(
                io.s_dif_segs[DIFL_F][DIFS_F_FORMULA],
            ) != 0;
        let restored = !switched
            && io.b_out_type == 4
            && io.b_output_type == 3
            && io.b_second_non_taut_pass != 0;
        let ret = output_inchi_fixed_h_layer_with_sublayers(
            &mut out,
            &mut sb,
            rec,
            &mut io,
            lf,
            &mut repeat,
        );
        // the arrays may be equal by value: derive the pointer from the branch
        let which = if switched {
            1
        } else if restored {
            2
        } else {
            0
        };
        let got = format!(
            "{ret} {repeat} {} {} {} {} {} {} {} {} {} {} {} {} {} {which} [{}] [{}] [{}]",
            io.b_out_type,
            io.i_cur_taut_mode,
            io.b_second_non_taut_pass,
            io.n_cur_inchi_segment,
            io.num_components,
            io.b_fh_tag,
            io.b_overflow,
            io.tot_len,
            io.tot_len2,
            io.b_tag1,
            io.b_always,
            io.b_non_taut_non_iso_identifier_not_empty * 100
                + io.b_non_taut_iso_identifier_not_empty,
            io.n_segm_action * 1000 + io.sort_print_inchi_flags,
            io.sz_tag1,
            strbuf_str(&out.s),
            strbuf_str(&sb)
        );
        let got: Vec<&str> = got.split('\n').filter(|l| !l.is_empty()).collect();
        assert_eq!(got, want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt1_auxinfo_layers_match_the_c_code() {
    use rdoxide::external::inchi::ichi::{InChIAux, InChISort, OrigInfo};
    use rdoxide::external::inchi::ichi_io::{
        INCHI_IOS_TYPE_STRING, InchiIosString, InchiIostream, IosFile, inchi_ios_init,
        inchi_strbuf_init,
    };
    use rdoxide::external::inchi::ichiprt1::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt1d.json")).unwrap();
    fn opt_u16(t: &mut InchiTok) -> Vec<u16> {
        if t.i() == 0 {
            return Vec::new();
        }
        let l = t.i();
        (0..l).map(|_| t.i() as u16).collect()
    }
    fn aux(t: &mut InchiTok) -> Option<InChIAux> {
        if t.i() == 0 {
            return None;
        }
        let mut a = InChIAux {
            n_number_of_atoms: t.i() as i32,
            n_number_of_t_groups: t.i() as i32,
            b_deleted: t.i() as i32,
            b_is_isotopic: t.i() as i32,
            ..Default::default()
        };
        a.n_constit_equ_numbers = opt_u16(t);
        a.n_constit_equ_t_group_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_t_group_numbers = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        if t.i() != 0 {
            let l = t.i();
            a.orig_info = (0..l)
                .map(|_| OrigInfo {
                    c_charge: t.i() as i8,
                    c_radical: t.i() as i8,
                    c_unusual_valence: t.i() as i8,
                })
                .collect();
        }
        Some(a)
    }
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let cmd = t.s();
        let n = t.i() as usize;
        let ot = t.i() as i32;
        let mult = t.i() as i32;
        let sorts = |t: &mut InchiTok| -> Vec<InChISort> {
            (0..n)
                .map(|i| InChISort {
                    p_inchi: [t.inchi(), t.inchi()],
                    p_inchi_aux: [aux(t), aux(t)],
                    ord_number: i as i16,
                })
                .collect()
        };
        let a = sorts(&mut t);
        let b = sorts(&mut t);
        let mut io = InchiOutCtl {
            taut_mode: t.i() as i32,
            n_tag: t.i() as i32,
            b_fh_tag: t.i() as i32,
            b_plain_text_tags: t.i() as i32,
            b_second_non_taut_pass: t.i() as i32,
            b_omit_repetitions: t.i() as i32,
            i_cur_taut_mode: t.i() as usize,
            b_tautomeric: t.i() as i32,
            b_tautomeric_output_allowed: t.i() as i32,
            b_isotopic: t.i() as i32,
            sort_print_inchi_flags: t.i() as i32,
            ..Default::default()
        };
        let nc2 = [t.i() as i32, t.i() as i32];
        let outopts = t.i() as i32;
        let rec = t.i() as i32;
        io.b_overflow = t.i() as i32;
        for arr in [
            &mut io.b_atom_equ,
            &mut io.b_taut_equ,
            &mut io.b_inv_stereo,
            &mut io.b_inv_stereo_orig_numb,
            &mut io.b_isotopic_orig_numb,
            &mut io.b_isotopic_atom_equ,
            &mut io.b_isotopic_taut_equ,
            &mut io.b_inv_isotopic_stereo,
            &mut io.b_inv_isotopic_stereo_orig_numb,
            &mut io.b_ign_uu_sp3_iso,
            &mut io.b_ign_uu_sp2_iso,
            &mut io.b_charges_rad_val,
        ] {
            for k in 0..2 {
                arr[k] = t.i() as i32;
            }
        }
        io.b_use_multipliers = mult;
        io.b_out_type = ot;
        io.num_components = n;
        io.p_inchi_sort = a;
        io.p_inchi_sort2 = b;
        let (lf, _) = set_line_separators(outopts);
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, 16, 16);
        let mut out = InchiIostream::default();
        inchi_ios_init(&mut out, INCHI_IOS_TYPE_STRING, IosFile::None);
        let ret = match cmd.as_str() {
            "head" => output_auxinfo_header_and_normalization_type(
                &mut out, &mut sb, outopts, rec, nc2, &mut io, lf,
            ),
            "orig" => output_auxinfo_original_numbers_and_equivalence_classes(
                &mut out, &mut sb, nc2, &mut io, lf,
            ),
            "tg" => output_auxinfo_tautomeric_groups_equivalence(&mut out, &mut sb, &mut io),
            "stereo" => output_auxinfo_stereo(&mut out, &mut sb, &mut io, lf),
            "iso" => output_auxinfo_isotopic_info(&mut out, &mut sb, rec, &mut io, lf),
            _ => {
                output_auxinfo_charges_radicals_and_unusual_valences(&mut out, &mut sb, &mut io, lf)
            }
        };
        let got = format!(
            "{ret} {} {} {} {} {} {} [{}] [{}] [{}] [{}] [{}]",
            io.b_overflow,
            io.tot_len,
            io.b_tag1,
            io.b_tag2,
            io.b_tag3,
            io.b_always,
            io.sz_tag1,
            io.sz_tag2,
            io.sz_tag3,
            strbuf_str(&out.s),
            strbuf_str(&sb)
        );
        let got: Vec<&str> = got.split('\n').filter(|l| !l.is_empty()).collect();
        assert_eq!(got, want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt1_original_structure_writers_match_the_c_code() {
    use rdoxide::external::inchi::ichiprt1::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt1e.json")).unwrap();
    fn coord(t: &mut InchiTok) -> MolCoord {
        let s = t.s();
        let mut c = [b'X'; MOL_COORD_LEN];
        for (d, b) in c.iter_mut().zip(s.bytes()) {
            *d = match b {
                b'_' => b' ',
                b'~' => 0,
                x => x,
            };
        }
        c
    }
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let cmd = t.s();
        let got = match cmd.as_str() {
            "clean" => {
                let mut mc = coord(&mut t);
                let delim = t.i() as u8;
                let nz = clean_orig_coord(&mut mc, delim);
                let hex: String = mc.iter().map(|b| format!("{b:02x}")).collect();
                format!("{nz} {hex}")
            }
            "wcoord" => {
                let n = t.i() as usize;
                let mut i = t.i() as usize;
                let buf_len = t.i() as usize;
                let mc: Vec<MolCoord> = (0..n).map(|_| coord(&mut t)).collect();
                let s = write_orig_coord(n, &mc, &mut i, buf_len);
                format!("{} {i} [{s}]", s.len())
            }
            _ => {
                let n = t.i() as usize;
                let mut i = t.i() as usize;
                let buf_len = t.i() as usize;
                let chiral = t.i() as i32;
                let mut at = vec![InpAtom::default(); n];
                for a in at.iter_mut() {
                    let el = t.s();
                    a.set_elname(&el);
                    a.el_number = t.i() as u8;
                    a.valence = t.i() as i8;
                    for m in 0..a.valence as usize {
                        a.neighbor[m] = t.i() as u16;
                    }
                    for m in 0..a.valence as usize {
                        a.bond_type[m] = t.i() as u8;
                    }
                    a.orig_at_number = t.i() as u16;
                    a.chem_bonds_valence = t.i() as i8;
                    a.num_h = t.i() as i8;
                    for m in 0..3 {
                        a.num_iso_h[m] = t.i() as i8;
                    }
                    a.iso_atw_diff = t.i() as i8;
                    a.charge = t.i() as i8;
                    a.radical = t.i() as i8;
                    a.p_parity = t.i() as i8;
                    for m in 0..4 {
                        a.p_orig_at_num[m] = t.i() as u16;
                    }
                }
                let s = write_orig_atoms(n, &at, &mut i, buf_len, chiral);
                format!("{} {i} [{s}]", s.len())
            }
        };
        assert_eq!(vec![got], want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt1_bond_writer_matches_the_c_code() {
    use rdoxide::external::inchi::ichiprt1::write_orig_bonds;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt1f.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let n = t.i() as usize;
        let mut i = t.i() as usize;
        let buf_len = t.i() as usize;
        let mut at = vec![InpAtom::default(); n];
        for a in at.iter_mut() {
            a.el_number = t.i() as u8;
            a.valence = t.i() as i8;
            for m in 0..a.valence as usize {
                a.neighbor[m] = t.i() as u16;
            }
            for m in 0..a.valence as usize {
                a.bond_type[m] = t.i() as u8;
            }
            for m in 0..a.valence as usize {
                a.bond_stereo[m] = t.i() as i8;
            }
            a.chem_bonds_valence = t.i() as i8;
            for m in 0..3 {
                a.sb_parity[m] = t.i() as i8;
            }
            for m in 0..3 {
                a.sb_ord[m] = t.i() as i8;
            }
            for m in 0..3 {
                a.sn_ord[m] = t.i() as i8;
            }
        }
        let mut err = String::new();
        let s = write_orig_bonds(n, &at, &mut i, buf_len, &mut err);
        let got = format!("{} {i} [{s}] [{err}]", s.len());
        assert_eq!(vec![got], want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt1_orig_struct_and_reversibility_match_the_c_code() {
    use rdoxide::external::inchi::ichi_io::{
        INCHI_IOS_TYPE_STRING, InchiIosString, InchiIostream, IosFile, inchi_ios_init,
        inchi_strbuf_init,
    };
    use rdoxide::external::inchi::ichiprt1::*;
    use rdoxide::external::inchi::inpdef::InpAtom;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt1g.json")).unwrap();
    fn coord(t: &mut InchiTok) -> MolCoord {
        let s = t.s();
        let mut c = [b'X'; MOL_COORD_LEN];
        for (d, b) in c.iter_mut().zip(s.bytes()) {
            *d = match b {
                b'_' => b' ',
                b'~' => 0,
                x => x,
            };
        }
        c
    }
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let n = t.i() as usize;
        let has_coord = t.i() != 0;
        let alloc = t.i() as i32;
        let plain = t.i() as i32;
        let second = t.i() as i32;
        let ntag = t.i() as i32;
        let fh = t.i() as i32;
        let lfopt = t.i() != 0;
        let chiral = t.i() as i32;
        let n_zy = t.i() as i32;
        let mut at = vec![InpAtom::default(); n];
        let mut coords = vec![[0u8; MOL_COORD_LEN]; n];
        for (k, a) in at.iter_mut().enumerate() {
            let el = t.s();
            a.set_elname(&el);
            a.el_number = t.i() as u8;
            a.valence = t.i() as i8;
            for m in 0..a.valence as usize {
                a.neighbor[m] = t.i() as u16;
            }
            for m in 0..a.valence as usize {
                a.bond_type[m] = t.i() as u8;
            }
            for m in 0..a.valence as usize {
                a.bond_stereo[m] = t.i() as i8;
            }
            a.chem_bonds_valence = t.i() as i8;
            for m in 0..3 {
                a.sb_parity[m] = t.i() as i8;
            }
            for m in 0..3 {
                a.sb_ord[m] = t.i() as i8;
            }
            for m in 0..3 {
                a.sn_ord[m] = t.i() as i8;
            }
            a.orig_at_number = t.i() as u16;
            a.num_h = t.i() as i8;
            a.charge = t.i() as i8;
            a.radical = t.i() as i8;
            for m in 0..3 {
                a.num_iso_h[m] = t.i() as i8;
            }
            a.iso_atw_diff = t.i() as i8;
            for m in 0..3 {
                let mc = coord(&mut t);
                if m == 0 {
                    coords[k] = mc;
                }
            }
        }
        let mut err = String::new();
        let fill = orig_struct_fill_out(
            n,
            n_zy,
            has_coord.then_some(&coords[..]),
            &at,
            chiral,
            &mut err,
        );
        let (ret, os) = match &fill {
            Ok(o) => (0, Some(o.clone())),
            Err(e) => (*e, None),
        };
        let o = os.clone().unwrap_or_default();
        let f = |s: &Option<String>| s.clone().unwrap_or_else(|| "-".to_string());
        let mut lines = vec![format!(
            "{ret} {} {} [{}] [{}] [{}] [{err}]",
            o.num_atoms,
            o.n_zy,
            f(&o.sz_atoms),
            f(&o.sz_bonds),
            f(&o.sz_coord)
        )];
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, alloc, 16);
        let mut out = InchiIostream::default();
        inchi_ios_init(&mut out, INCHI_IOS_TYPE_STRING, IosFile::None);
        let mut io = InchiOutCtl {
            b_plain_text_tags: plain,
            b_second_non_taut_pass: second,
            n_tag: ntag,
            b_fh_tag: fh,
            ..Default::default()
        };
        let (lf, _) = set_line_separators(if lfopt { 0x80 } else { 0 });
        let r = output_auxinfo_reversibility_info(&mut out, &mut sb, os.as_ref(), &mut io, lf);
        lines.push(format!(
            "{r} [{}] [{}] {}",
            io.sz_tag2,
            strbuf_str(&out.s),
            sb.n_used_length
        ));
        let got: Vec<String> = lines
            .join("\n")
            .split('\n')
            .filter(|l| !l.is_empty())
            .map(String::from)
            .collect();
        assert_eq!(got, want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt1_polymer_output_helpers_match_the_c_code() {
    use rdoxide::external::inchi::ichi::{InChIAux, InChISort, OrigInfo};
    use rdoxide::external::inchi::ichi_io::{
        INCHI_IOS_TYPE_STRING, InchiIosString, InchiIostream, IosFile, inchi_ios_init,
        inchi_ios_print_nodisplay, inchi_strbuf_init,
    };
    use rdoxide::external::inchi::ichiprt1::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt1h.json")).unwrap();
    fn opt_u16(t: &mut InchiTok) -> Vec<u16> {
        if t.i() == 0 {
            return Vec::new();
        }
        let l = t.i();
        (0..l).map(|_| t.i() as u16).collect()
    }
    fn aux(t: &mut InchiTok) -> Option<InChIAux> {
        if t.i() == 0 {
            return None;
        }
        let mut a = InChIAux {
            n_number_of_atoms: t.i() as i32,
            n_number_of_t_groups: t.i() as i32,
            b_deleted: t.i() as i32,
            b_is_isotopic: t.i() as i32,
            ..Default::default()
        };
        a.n_constit_equ_numbers = opt_u16(t);
        a.n_constit_equ_t_group_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_t_group_numbers = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        if t.i() != 0 {
            let l = t.i();
            a.orig_info = (0..l)
                .map(|_| OrigInfo {
                    c_charge: t.i() as i8,
                    c_radical: t.i() as i8,
                    c_unusual_valence: t.i() as i8,
                })
                .collect();
        }
        Some(a)
    }
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let cmd = t.s();
        let got = match cmd.as_str() {
            "err" => {
                let kind = t.i() as i32;
                let text = t.s().replace("~n", "\n");
                let mut out = InchiIostream::default();
                inchi_ios_init(&mut out, INCHI_IOS_TYPE_STRING, IosFile::None);
                let ret = output_inchi_plain_error(&mut out, &text, kind);
                format!("{ret} [{}]", strbuf_str(&out.s))
            }
            "bondlt" => {
                let b1 = [t.i() as i32, t.i() as i32];
                let b2 = [t.i() as i32, t.i() as i32];
                format!("{}", is_bond_atom_nums_lesser(&b1, &b2))
            }
            "count" => {
                let f = t.s();
                format!("{}", count_pseudo_element_in_formula("Zz", f.as_bytes()))
            }
            "hide" => {
                let pzz = t.i() as i32;
                let zy = t.i() as i32;
                let text = t.s().replace("~n", "\n");
                let mut out = InchiIostream::default();
                inchi_ios_init(&mut out, INCHI_IOS_TYPE_STRING, IosFile::None);
                inchi_ios_print_nodisplay(Some(&mut out), &text);
                edit_inchi_hide_polymer_zz(&mut out, pzz, zy);
                format!("{} [{}]", out.s.n_used_length, strbuf_str(&out.s))
            }
            _ => {
                let n = t.i() as usize;
                let ot = t.i() as i32;
                let sorts: Vec<InChISort> = (0..n)
                    .map(|i| InChISort {
                        p_inchi: [t.inchi(), t.inchi()],
                        p_inchi_aux: [aux(&mut t), aux(&mut t)],
                        ord_number: i as i16,
                    })
                    .collect();
                let mut io = InchiOutCtl {
                    taut_mode: t.i() as i32,
                    b_omit_repetitions: t.i() as i32,
                    b_out_type: ot,
                    num_components: n,
                    ..Default::default()
                };
                let nat = t.i() as usize;
                io.p_inchi_sort = sorts.clone();
                io.p_inchi_sort2 = sorts;
                let mut sb = InchiIosString::default();
                inchi_strbuf_init(&mut sb, 64, 16);
                let mut cano = vec![0i32; nat + 2];
                let mut comp = vec![0i32; nat + 2];
                let ret = internally_get_cano_nums_and_component_nums(
                    &mut sb,
                    &mut io,
                    nat,
                    &mut cano[..nat + 1],
                    &mut comp[..nat + 1],
                );
                let j = |v: &[i32]| {
                    v[..=nat]
                        .iter()
                        .map(|x| format!(" {x}"))
                        .collect::<String>()
                };
                format!("{ret}{} |{}", j(&cano), j(&comp))
            }
        };
        let got: Vec<&str> = got.split('\n').filter(|l| !l.is_empty()).collect();
        assert_eq!(got, want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn ichiprt1_whole_inchi_serialization_matches_the_c_code() {
    use rdoxide::external::inchi::ichi::{InChIAux, InChISort, OrigInfo};
    use rdoxide::external::inchi::ichi_io::{
        INCHI_IOS_TYPE_STRING, InchiIosString, InchiIostream, IosFile, inchi_ios_init,
        inchi_strbuf_init,
    };
    use rdoxide::external::inchi::ichiprt1::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt1i.json")).unwrap();
    fn opt_u16(t: &mut InchiTok) -> Vec<u16> {
        if t.i() == 0 {
            return Vec::new();
        }
        let l = t.i();
        (0..l).map(|_| t.i() as u16).collect()
    }
    fn aux(t: &mut InchiTok) -> Option<InChIAux> {
        if t.i() == 0 {
            return None;
        }
        let mut a = InChIAux {
            n_number_of_atoms: t.i() as i32,
            n_number_of_t_groups: t.i() as i32,
            b_deleted: t.i() as i32,
            b_is_isotopic: t.i() as i32,
            ..Default::default()
        };
        a.n_constit_equ_numbers = opt_u16(t);
        a.n_constit_equ_t_group_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_t_group_numbers = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        if t.i() != 0 {
            let l = t.i();
            a.orig_info = (0..l)
                .map(|_| OrigInfo {
                    c_charge: t.i() as i8,
                    c_radical: t.i() as i8,
                    c_unusual_valence: t.i() as i8,
                })
                .collect();
        }
        Some(a)
    }
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let n = t.i() as usize;
        let mut sorts: [[Vec<InChISort>; 2]; 2] = Default::default();
        for k in 0..4 {
            sorts[k / 2][k % 2] = (0..n)
                .map(|i| InChISort {
                    p_inchi: [t.inchi(), t.inchi()],
                    p_inchi_aux: [aux(&mut t), aux(&mut t)],
                    ord_number: i as i16,
                })
                .collect();
        }
        let nc = [t.i() as i32, t.i() as i32];
        let _nnt = [t.i() as i32, t.i() as i32];
        let _nt = [t.i() as i32, t.i() as i32];
        let mut ip = OutputInput {
            abc_numbers: t.i() != 0,
            ct_predecessors: t.i() != 0,
            polymers: t.i() as i32,
            large_molecules: t.i() != 0,
            ..Default::default()
        };
        if t.i() != 0 {
            ip.sdf_label = "lbl".to_string();
        }
        if t.i() != 0 {
            ip.sdf_value = "val".to_string();
        }
        let disc = t.i() != 0;
        let otype = t.i() as i32;
        let opts = t.i() as i32;
        let saveopt = t.i() as u8;
        let alloc = t.i() as i32;
        let mut flags = t.i() as i32;
        let bor = t.i() as i32;
        let use_os = t.i() != 0;
        let mut os = OrigStruct::default();
        if use_os {
            os.num_atoms = t.i() as i32;
            os.n_zy = t.i() as i32;
            let a = t.s().replace('_', " ");
            let b = t.s().replace('_', " ");
            let c2 = t.s().replace('_', " ");
            os.sz_atoms = (t.i() != 0).then_some(a);
            os.sz_bonds = (t.i() != 0).then_some(b);
            os.sz_coord = (t.i() != 0).then_some(c2);
        }
        let use_od = t.i() != 0;
        let od_zy = t.i() as i32;
        let use_log = t.i() != 0;
        let dis2 = t.i() != 0;
        ip.orig_present = use_od;
        ip.n_zy = od_zy;
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, alloc, 16);
        let mut out = InchiIostream::default();
        inchi_ios_init(&mut out, INCHI_IOS_TYPE_STRING, IosFile::None);
        let mut lg = InchiIostream::default();
        inchi_ios_init(&mut lg, INCHI_IOS_TYPE_STRING, IosFile::None);
        let osr = use_os.then_some(&os);
        let ret = if dis2 {
            output_inchi2(
                &mut sb,
                &sorts,
                bor,
                osr,
                &ip,
                disc,
                otype,
                opts,
                nc,
                &mut out,
                use_log.then_some(&mut lg),
                7,
                &mut flags,
                saveopt,
            )
        } else {
            output_inchi1(
                &mut sb,
                &sorts,
                bor,
                osr,
                &ip,
                disc,
                otype,
                opts,
                nc,
                &mut out,
                use_log.then_some(&mut lg),
                7,
                &mut flags,
                saveopt,
            )
        };
        let got = format!(
            "{ret} {flags} [{}] [{}]",
            strbuf_str(&out.s),
            strbuf_str(&lg.s)
        );
        let got: Vec<&str> = got.split('\n').filter(|l| !l.is_empty()).collect();
        assert_eq!(got, want, "case {ci}: {}", c["in"]);
    }
}

#[test]
fn runichi4_sort_and_print_inchi_matches_the_c_code() {
    use rdoxide::external::inchi::ichi::{InChIAux, OrigInfo};
    use rdoxide::external::inchi::ichi_io::{
        INCHI_IOS_TYPE_STRING, InchiIosString, InchiIostream, IosFile, inchi_ios_init,
        inchi_strbuf_init,
    };
    use rdoxide::external::inchi::ichiprt1::{OrigStruct, strbuf_str};
    use rdoxide::external::inchi::runichi4::*;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_prt1j.json")).unwrap();
    fn opt_u16(t: &mut InchiTok) -> Vec<u16> {
        if t.i() == 0 {
            return Vec::new();
        }
        let l = t.i();
        (0..l).map(|_| t.i() as u16).collect()
    }
    fn aux(t: &mut InchiTok) -> Option<InChIAux> {
        if t.i() == 0 {
            return None;
        }
        let mut a = InChIAux {
            n_number_of_atoms: t.i() as i32,
            n_number_of_t_groups: t.i() as i32,
            b_deleted: t.i() as i32,
            b_is_isotopic: t.i() as i32,
            ..Default::default()
        };
        a.n_constit_equ_numbers = opt_u16(t);
        a.n_constit_equ_t_group_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_numbers = opt_u16(t);
        a.n_constit_equ_isotopic_t_group_numbers = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord = opt_u16(t);
        a.n_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        a.n_isotopic_orig_at_nos_in_canon_ord_inv = opt_u16(t);
        if t.i() != 0 {
            let l = t.i();
            a.orig_info = (0..l)
                .map(|_| OrigInfo {
                    c_charge: t.i() as i8,
                    c_radical: t.i() as i8,
                    c_unusual_valence: t.i() as i8,
                })
                .collect();
        }
        Some(a)
    }
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mut t = InchiTok(c["in"].as_str().unwrap().split_whitespace());
        let n = t.i() as usize;
        let mut pinchi: [PInChI2; 2] = Default::default();
        let mut pinchi_aux: [PInChIAux2; 2] = Default::default();
        for j in 0..2 {
            for _ in 0..n {
                pinchi[j].push([t.inchi(), t.inchi()]);
                pinchi_aux[j].push([aux(&mut t), aux(&mut t)]);
            }
        }
        let nc = [t.i() as i32, t.i() as i32];
        let _nnt = [t.i() as i32, t.i() as i32];
        let _nt = [t.i() as i32, t.i() as i32];
        let mut tf = [t.i() as u64, t.i() as u64];
        let mut tfd = [t.i() as u64, t.i() as u64];
        let mut nctf = [[0u64; 2]; 2];
        let mut nctfd = [[0u64; 2]; 2];
        for i in 0..2 {
            for k in 0..2 {
                nctf[i][k] = t.i() as u64;
            }
        }
        for i in 0..2 {
            for k in 0..2 {
                nctfd[i][k] = t.i() as u64;
            }
        }
        let mut ip = SortPrintParams {
            inchi_output_options: t.i() as i32,
            compare_components: t.i() as i32,
            ..Default::default()
        };
        ip.input.abc_numbers = t.i() != 0;
        ip.input.ct_predecessors = t.i() != 0;
        ip.input.polymers = t.i() as i32;
        ip.input.large_molecules = t.i() != 0;
        if t.i() != 0 {
            ip.input.sdf_label = "lbl".to_string();
        }
        if t.i() != 0 {
            ip.input.sdf_value = "val".to_string();
        }
        let saveopt = t.i() as u8;
        let alloc = t.i() as i32;
        let mut flags = t.i() as i32;
        let use_os = t.i() != 0;
        let mut os = OrigStruct::default();
        if use_os {
            os.num_atoms = t.i() as i32;
            os.n_zy = t.i() as i32;
            let a = t.s().replace('_', " ");
            let b = t.s().replace('_', " ");
            let c2 = t.s().replace('_', " ");
            os.sz_atoms = (t.i() != 0).then_some(a);
            os.sz_bonds = (t.i() != 0).then_some(b);
            os.sz_coord = (t.i() != 0).then_some(c2);
        }
        let use_od = t.i() != 0;
        let od_zy = t.i() as i32;
        let use_log = t.i() != 0;
        ip.input.orig_present = use_od;
        ip.input.n_zy = od_zy;
        let mut sb = InchiIosString::default();
        inchi_strbuf_init(&mut sb, alloc, 16);
        let mut out = InchiIostream::default();
        inchi_ios_init(&mut out, INCHI_IOS_TYPE_STRING, IosFile::None);
        let mut lg = InchiIostream::default();
        inchi_ios_init(&mut out, INCHI_IOS_TYPE_STRING, IosFile::None);
        inchi_ios_init(&mut lg, INCHI_IOS_TYPE_STRING, IosFile::None);
        let chiral = is_struct_chiral(&pinchi, nc);
        let ret = sort_and_print_inchi(
            &mut out,
            &mut sb,
            use_log.then_some(&mut lg),
            &ip,
            use_os.then_some(&os),
            nc,
            &mut tf,
            &mut tfd,
            nctf,
            nctfd,
            7,
            &mut pinchi,
            &mut pinchi_aux,
            &mut flags,
            saveopt,
        );
        let mut got = format!(
            "{ret} {flags} {} {} {} {} [{}] [{}]\n",
            tf[0],
            tf[1],
            tfd[0],
            tfd[1],
            strbuf_str(&out.s),
            strbuf_str(&lg.s)
        );
        got.push_str(&format!("X{}\n", i32::from(chiral)));
        for j in 0..2 {
            for i in 0..nc[j] as usize {
                got.push_str(&format!("C{j}.{i}"));
                for k in 0..2 {
                    match &pinchi[j][i][k] {
                        Some(p) => got.push_str(&format!(
                            " {}:{}:{}:{}:{}",
                            p.sz_hill_formula.clone().unwrap_or_default(),
                            p.n_number_of_atoms,
                            p.b_deleted,
                            p.len_tautomer,
                            p.n_number_of_isotopic_atoms
                        )),
                        None => got.push_str(" -"),
                    }
                    match &pinchi_aux[j][i][k] {
                        Some(q) => got.push_str(&format!(
                            "/{}:{}:{}",
                            q.n_number_of_atoms, q.b_deleted, q.n_number_of_t_groups
                        )),
                        None => got.push_str("/-"),
                    }
                }
                got.push('\n');
            }
        }
        let got: Vec<&str> = got.split('\n').filter(|l| !l.is_empty()).collect();
        assert_eq!(got, want, "case {ci}: {}", c["in"]);
    }
}
