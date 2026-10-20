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
