//! Port of RDKit `Code/DataStructs`.

pub mod bit_ops;
pub mod bit_vect;
pub mod bitmap;
pub mod fpb_reader;
pub mod multi_fpb_reader;
pub mod sparse_int_vect;
pub mod value_vects;

pub use bit_vect::{BitVect, ExplicitBitVect, SparseBitVect};

use crate::{RDKitError, Result};

/// `convertToExplicit`
pub fn convert_to_explicit(sbv: &SparseBitVect) -> ExplicitBitVect {
    let mut ebv = ExplicitBitVect::new(sbv.num_bits(), false);
    for b in sbv.on_bits() {
        let _ = ebv.set_bit(b as u32);
    }
    ebv
}

/// `FromBitString`
pub fn from_bit_string<T: BitVect>(bv: &mut T, s: &[u8]) -> Result<()> {
    if s.len() as u64 > bv.num_bits() as u64 {
        return Err(RDKitError::precondition("bad bitvect length"));
    }
    *bv = T::new_sized(bv.num_bits());
    for i in 0..bv.num_bits() {
        if s.get(i as usize) == Some(&b'1') {
            bv.set_bit(i)?;
        }
    }
    Ok(())
}

/// The 6-bit code of a character of Daylight's ASCII fingerprint encoding
/// (`.`, `+`, `0-9`, `A-Z`, `a-z`); `None` for any other character.
fn daylight_char_code(c: u8) -> Option<u8> {
    match c {
        b'.' => Some(0),
        b'+' => Some(1),
        b'0'..=b'9' => Some(c - b'0' + 2),
        b'A'..=b'Z' => Some(c - b'A' + 12),
        b'a'..=b'z' => Some(c - b'a' + 38),
        _ => None,
    }
}

/// `a2b`: decodes four characters of the Daylight encoding into three bytes.
/// Like upstream, a character outside the alphabet repeats the code of the
/// preceding character (0 for the first).
fn daylight_decode_group(ascii: &[u8]) -> [u8; 3] {
    let mut code = 0u8;
    let mut bits = 0u32;
    for &c in ascii {
        code = daylight_char_code(c).unwrap_or(code);
        bits = bits << 6 | u32::from(code);
    }
    [(bits >> 16) as u8, (bits >> 8) as u8, bits as u8]
}

/// `FromDaylightString`: fills `bv` from the ASCII (Daylight) encoding of a
/// fingerprint, the text Daylight tools put in the FP field of a TDT.
pub fn from_daylight_string<T: BitVect>(bv: &mut T, s: &str) -> Result<()> {
    *bv = T::new_sized(bv.num_bits());
    let s = s.as_bytes();
    let mut length = s.len();
    if s.last() == Some(&b'\n') {
        length -= 1;
    }
    // 4 characters encode 3 bytes; one extra character marks the padding
    let Some(data_len) = length.checked_sub(1) else {
        return Err(RDKitError::Runtime(
            "bad daylight fingerprint string".into(),
        ));
    };
    let mut n_bits = 3 * data_len / 4 * 8;
    match s[data_len] {
        b'1' => n_bits = n_bits.saturating_sub(16),
        b'2' => n_bits = n_bits.saturating_sub(8),
        b'3' => {}
        _ => {
            return Err(RDKitError::Runtime(
                "bad daylight fingerprint string".into(),
            ));
        }
    }
    let mut done = 0;
    for group in s[..data_len].chunks(4) {
        // a trailing partial group reads the characters beyond it, as upstream does
        let mut padded = [0u8; 4];
        padded[..group.len()].copy_from_slice(group);
        for byte in daylight_decode_group(&padded) {
            for k in (0..8).rev() {
                if done >= n_bits {
                    return Ok(());
                }
                if byte >> k & 1 == 1 {
                    bv.set_bit(done as u32)?;
                }
                done += 1;
            }
        }
    }
    Ok(())
}
