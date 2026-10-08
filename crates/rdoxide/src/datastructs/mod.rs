//! Port of RDKit `Code/DataStructs`.

pub mod bit_ops;
pub mod bit_vect;
pub mod sparse_int_vect;

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
