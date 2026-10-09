//! `from_daylight_string`; expected bits from the reference RDKit.

use rdoxide::datastructs::{BitVect, ExplicitBitVect, SparseBitVect, from_daylight_string};

#[test]
fn decodes_bits_into_both_vector_kinds() {
    let expected = [5, 12, 15, 16, 18, 21, 22, 23];
    let mut e = ExplicitBitVect::new(32, false);
    from_daylight_string(&mut e, "+.ab3").unwrap();
    assert_eq!(e.on_bits(), expected);
    let mut s = SparseBitVect::new(32);
    from_daylight_string(&mut s, "+.ab3").unwrap();
    assert_eq!(s.on_bits(), expected);
}

#[test]
fn padding_marker_trims_bytes_and_previous_bits_are_cleared() {
    let mut e = ExplicitBitVect::new(32, false);
    e.set_bit(31).unwrap();
    from_daylight_string(&mut e, "+.ab2\n").unwrap();
    assert_eq!(e.on_bits(), [5, 12, 15]);
}

#[test]
fn malformed_strings_are_errors() {
    let mut e = ExplicitBitVect::new(64, false);
    for bad in ["", "abcd", "abcd4"] {
        assert!(from_daylight_string(&mut e, bad).is_err(), "{bad:?}");
    }
}
