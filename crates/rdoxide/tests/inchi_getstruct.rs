//! `GetStructFromINCHI` (InChI -> atoms and 0D stereo): the structure of a
//! standard InChI is converted back to an InChI with `GetStdINCHI`, which
//! must reproduce the input string.

use rdoxide::external::inchi::inchi_dll::InchiInput;
use rdoxide::external::inchi::inchi_dll::get_std_inchi;
use rdoxide::external::inchi::inchi_dll_c::{
    INCHI_INVALID_LAYOUT, INCHI_INVALID_PREFIX, INCHI_INVALID_VERSION, INCHI_RET_EOF,
    INCHI_RET_ERROR, INCHI_RET_OKAY, INCHI_VALID_NON_STANDARD, INCHI_VALID_STANDARD, check_inchi,
    get_struct_from_inchi, get_struct_from_std_inchi,
};

const CASES: &[&str] = &[
    "InChI=1S/CH4/h1H4",
    "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
    "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H",
    "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)",
    "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)",
    "InChI=1S/H2O/h1H2",
    "InChI=1S/C2H4/c1-2/h1-2H2",
    "InChI=1S/CH5N/c1-2/h2H2,1H3",
    "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i1D",
    "InChI=1S/C4H8O2/c1-3-4(5)6-2/h3H2,1-2H3",
    "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1",
    "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2",
    "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
    "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
    "InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)",
    "InChI=1S/C5H5N5/c6-4-3-5(9-1-7-3)10-2-8-4/h1-2H,(H3,6,7,8,9,10)",
    "InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)",
    "InChI=1S/ClH/h1H",
    "InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3",
    "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+",
    "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3-",
    "InChI=1S/C4H10O/c1-3-4(2)5/h4-5H,3H2,1-2H3/t4-/m1/s1",
    "InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3",
    "InChI=1S/C6H5NO2/c8-7(9)6-4-2-1-3-5-6/h1-5H",
    "InChI=1S/C6H7N/c7-6-4-2-1-3-5-6/h1-5H,7H2",
    "InChI=1S/C3H6O/c1-3(2)4/h1-2H3",
    "InChI=1S/CH3O/c1-2/h1H3/q-1",
    "InChI=1S/O4S/c1-5(2,3)4/q-2",
    "InChI=1S/C2H6OS/c1-4(2)3/h1-2H3",
    "InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)",
];

#[test]
fn inchi_to_structure_and_back() {
    let mut bad = Vec::new();
    for &s in CASES {
        let (ret, out) = get_struct_from_std_inchi(Some(s), None);
        if ret != INCHI_RET_OKAY || out.atom.is_empty() {
            bad.push(format!("{s}: ret {ret}, log {:?}", out.sz_log));
            continue;
        }
        let inp = InchiInput {
            atom: out.atom,
            stereo_0d: out.stereo_0d,
            options: None,
        };
        let (r2, o2) = get_std_inchi(&inp);
        if r2 > 1 || o2.sz_inchi.as_deref() != Some(s) {
            bad.push(format!(
                "{s}: back {r2} {:?} {:?}",
                o2.sz_inchi, o2.sz_message
            ));
        }
    }
    for b in &bad {
        eprintln!("BAD {b}");
    }
    assert!(bad.is_empty());
}

#[test]
fn check_inchi_codes() {
    assert_eq!(
        check_inchi("InChI=1S/CH4/h1H4", false),
        INCHI_VALID_STANDARD
    );
    assert_eq!(
        check_inchi("InChI=1/CH4/h1H4", false),
        INCHI_VALID_NON_STANDARD
    );
    assert_eq!(check_inchi("CH4", false), INCHI_INVALID_PREFIX);
    assert_eq!(
        check_inchi("InChI=2S/CH4/h1H4", false),
        INCHI_INVALID_VERSION
    );
    assert_eq!(check_inchi("InChI=1SCH4", false), INCHI_INVALID_LAYOUT);
    assert_eq!(
        check_inchi("InChI=1S/CH4/H1H4", false),
        INCHI_INVALID_LAYOUT
    );
}

#[test]
fn bad_input_is_an_error() {
    let (r, o) = get_struct_from_inchi(Some("not an inchi"), None);
    // as in C, the invalid input ends in the default branch of the result
    // translation: no structure -> inchi_Ret_EOF
    assert_eq!(r, INCHI_RET_EOF);
    assert!(o.atom.is_empty());
    let (r, _) = get_struct_from_std_inchi(Some("InChI=1/CH4/h1H4"), None);
    assert_eq!(r, INCHI_RET_ERROR);
}
