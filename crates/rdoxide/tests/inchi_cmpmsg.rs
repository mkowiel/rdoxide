//! Differential test of `FillOutCompareMessage` (`ichirvr7.rs`) against the C
//! code, with the fixture from `tools/inchi_cmpmsg_reference.py`.

use rdoxide::external::inchi::ichirvr7::fill_out_compare_message;

#[test]
fn compare_message_matches_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_cmpmsg.json")).unwrap();
    let cases = data.as_array().unwrap();
    let mut bad = 0;
    for c in cases {
        let i = c["in"].as_array().unwrap();
        let bits = [i[0].as_u64().unwrap(), i[1].as_u64().unwrap()];
        let len = i[2].as_i64().unwrap() as i32;
        let init = i[3].as_str().unwrap();
        // the C input uses TAUT_NON = 0, TAUT_YES = 1 in the same order
        let mut msg = if init == "-" {
            String::new()
        } else {
            init.to_string()
        };
        let r = fill_out_compare_message(&mut msg, len, bits);
        let o = c["out"].as_array().unwrap();
        if r != o[0].as_i64().unwrap() as i32 || msg != o[1].as_str().unwrap() {
            bad += 1;
            if bad < 4 {
                eprintln!("{i:?}\n got {r} {msg:?}\n want {o:?}");
            }
        }
    }
    assert_eq!(bad, 0, "{bad} of {} cases differ", cases.len());
}
