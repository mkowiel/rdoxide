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
