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
