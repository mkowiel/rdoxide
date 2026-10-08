//! Bit-vector comparisons and similarity metrics (`BitOps`).
//!
//! Ported from RDKit `Code/DataStructs/BitOps.{h,cpp}` (2026.09.1).
//! Arithmetic follows the reference builds (see [`crate::rdgeneral::fp`]).
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::bit_vect::{BitVect, ExplicitBitVect, SparseBitVect};
use crate::rdgeneral::fp::fmuladd;
use crate::{RDKitError, Result};

/// Bit vectors with a fast common-on-bits count.
pub trait BitOpsVect: BitVect {
    /// `NumOnBitsInCommon`
    fn num_on_bits_in_common(&self, o: &Self) -> u32;
    /// `NumBitsInCommon` without the length check
    fn num_bits_in_common_raw(&self, o: &Self) -> i32 {
        self.num_bits() as i32 - self.xor(o).num_on_bits() as i32
    }
}

impl BitOpsVect for ExplicitBitVect {
    fn num_on_bits_in_common(&self, o: &Self) -> u32 {
        self.num_on_in_common(o)
    }
    fn num_bits_in_common_raw(&self, o: &Self) -> i32 {
        self.num_bits() as i32 - self.num_xor(o) as i32
    }
}

impl BitOpsVect for SparseBitVect {
    fn num_on_bits_in_common(&self, o: &Self) -> u32 {
        self.and(o).num_on_bits()
    }
}

fn same_len<T: BitVect>(a: &T, b: &T) -> Result<()> {
    if a.num_bits() != b.num_bits() {
        return Err(RDKitError::Value("BitVects must be same length".into()));
    }
    Ok(())
}

/// The similarity metrics of `BitOps.h`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Metric {
    Tanimoto,
    Cosine,
    Kulczynski,
    Dice,
    Sokal,
    McConnaughey,
    Asymmetric,
    BraunBlanquet,
    Russel,
    RogotGoldberg,
    OnBit,
    AllBit,
}

/// Similarity between two equally long vectors.
pub fn similarity<T: BitOpsVect>(m: Metric, a: &T, b: &T) -> Result<f64> {
    same_len(a, b)?;
    let common = a.num_on_bits_in_common(b);
    let (yu, zu) = (a.num_on_bits(), b.num_on_bits());
    let x = common as f64;
    let (y, z) = (yu as f64, zu as f64);
    Ok(match m {
        Metric::Tanimoto => {
            let total = yu + zu;
            if total == 0 {
                0.0
            } else {
                common as f64 / (total - common) as f64
            }
        }
        Metric::Cosine => {
            if y * z > 0.0 {
                x / (y * z).sqrt()
            } else {
                0.0
            }
        }
        Metric::Kulczynski => {
            if y * z > 0.0 {
                x * (y + z) / (2.0 * y * z)
            } else {
                0.0
            }
        }
        Metric::Dice => {
            if y + z > 0.0 {
                2.0 * x / (y + z)
            } else {
                0.0
            }
        }
        Metric::Sokal => {
            if yu == 0 || zu == 0 {
                0.0
            } else {
                // (2.*y + 2.*z) - 3.*x
                x / fmuladd(-3.0, x, fmuladd(2.0, y, 2.0 * z))
            }
        }
        Metric::McConnaughey => {
            if y * z > 0.0 {
                fmuladd(x, y + z, -(y * z)) / (y * z)
            } else {
                0.0
            }
        }
        Metric::Asymmetric => {
            let min = if z < y { z } else { y };
            if min > 0.0 { x / min } else { 0.0 }
        }
        Metric::BraunBlanquet => {
            let max = if y < z { z } else { y };
            if max > 0.0 { x / max } else { 0.0 }
        }
        Metric::Russel => x / a.num_bits() as f64,
        Metric::RogotGoldberg => {
            if yu == 0 || zu == 0 {
                return Ok(0.0);
            }
            let l = a.num_bits() as f64;
            let d = l - y - z + x;
            let denom1 = y + z;
            let denom2 = fmuladd(2.0, l, -y) - z;
            if x == l || d == l {
                1.0
            } else if denom1 == 0.0 || denom2 == 0.0 {
                0.0
            } else {
                x / (y + z) + d / (fmuladd(2.0, l, -y) - z)
            }
        }
        Metric::OnBit => {
            let denom = a.or(b).num_on_bits() as f64;
            if denom > 0.0 { x / denom } else { 0.0 }
        }
        Metric::AllBit => a.num_bits_in_common_raw(b) as f64 / a.num_bits() as f64,
    })
}

/// `TverskySimilarity`
pub fn tversky<T: BitOpsVect>(a: &T, b: &T, ca: f64, cb: f64) -> Result<f64> {
    if !(0.0..=1.0).contains(&ca) {
        return Err(RDKitError::range("a", ca, 1));
    }
    if !(0.0..=1.0).contains(&cb) {
        return Err(RDKitError::range("b", cb, 1));
    }
    same_len(a, b)?;
    let x = a.num_on_bits_in_common(b) as f64;
    let (yu, zu) = (a.num_on_bits(), b.num_on_bits());
    if yu == 0 || zu == 0 {
        return Ok(0.0);
    }
    let (y, z) = (yu as f64, zu as f64);
    // (a*y + b*z) + (1-a-b)*x
    let denom = fmuladd(1.0 - ca - cb, x, fmuladd(ca, y, cb * z));
    Ok(if denom == 0.0 { 1.0 } else { x / denom })
}

/// `NumBitsInCommon`
pub fn num_bits_in_common<T: BitOpsVect>(a: &T, b: &T) -> Result<i32> {
    same_len(a, b)?;
    Ok(a.num_bits_in_common_raw(b))
}

/// `OnBitsInCommon`
pub fn on_bits_in_common<T: BitOpsVect>(a: &T, b: &T) -> Result<Vec<i32>> {
    same_len(a, b)?;
    Ok(a.and(b).on_bits())
}

/// `OffBitsInCommon`
pub fn off_bits_in_common<T: BitOpsVect>(a: &T, b: &T) -> Result<Vec<i32>> {
    same_len(a, b)?;
    Ok(a.or(b).not().on_bits())
}

/// `OnBitProjSimilarity`
pub fn on_bit_proj_similarity<T: BitOpsVect>(a: &T, b: &T) -> Result<[f64; 2]> {
    same_len(a, b)?;
    let num = a.num_on_bits_in_common(b) as f64;
    let mut res = [0.0; 2];
    if num != 0.0 {
        res[0] = num / a.num_on_bits() as f64;
        res[1] = num / b.num_on_bits() as f64;
    }
    Ok(res)
}

/// `OffBitProjSimilarity`
pub fn off_bit_proj_similarity<T: BitOpsVect>(a: &T, b: &T) -> Result<[f64; 2]> {
    same_len(a, b)?;
    let num = a.or(b).num_off_bits() as f64;
    let mut res = [0.0; 2];
    if num != 0.0 {
        res[0] = num / a.num_off_bits() as f64;
        res[1] = num / b.num_off_bits() as f64;
    }
    Ok(res)
}

/// `AllProbeBitsMatch(probe, ref)`
pub fn all_probe_bits_match<T: BitVect>(probe: &T, reference: &T) -> Result<bool> {
    for i in 0..probe.num_bits() {
        if probe.get_bit(i)? && !reference.get_bit(i)? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// `FoldFingerprint`
pub fn fold_fingerprint<T: BitVect>(bv: &T, factor: u32) -> Result<T> {
    if factor == 0 || factor >= bv.num_bits() {
        return Err(RDKitError::Value("invalid fold factor".into()));
    }
    let res_size = bv.num_bits() / factor;
    let mut res = T::new_sized(res_size);
    for b in bv.on_bits() {
        res.set_bit(b as u32 % res_size)?;
    }
    Ok(res)
}

/// `SimilarityWrapper`: folds the longer vector first.
pub fn similarity_wrapper<T: BitOpsVect>(
    a: &T,
    b: &T,
    metric: impl Fn(&T, &T) -> Result<f64>,
    return_distance: bool,
) -> Result<f64> {
    let mut res = if a.num_bits() > b.num_bits() {
        metric(&fold_fingerprint(a, a.num_bits() / b.num_bits())?, b)?
    } else if b.num_bits() > a.num_bits() {
        metric(a, &fold_fingerprint(b, b.num_bits() / a.num_bits())?)?
    } else {
        metric(a, b)?
    };
    if return_distance {
        res = 1.0 - res;
    }
    Ok(res)
}

/// `BitVectToText`
pub fn bit_vect_to_text<T: BitVect>(bv: &T) -> String {
    let mut res = vec![b'0'; bv.num_bits() as usize];
    for b in bv.on_bits() {
        res[b as usize] = b'1';
    }
    String::from_utf8(res).expect("ascii")
}

const BIN2HEX: &[u8; 16] = b"0123456789abcdef";

/// `BitVectToFPSText`
pub fn bit_vect_to_fps_text<T: BitVect>(bv: &T) -> String {
    let bytes = bit_vect_to_binary_text(bv);
    let mut res = Vec::with_capacity(bytes.len() * 2);
    for c in bytes {
        res.push(BIN2HEX[(c >> 4) as usize]);
        res.push(BIN2HEX[(c & 15) as usize]);
    }
    String::from_utf8(res).expect("ascii")
}

/// `BitVectToBinaryText`
pub fn bit_vect_to_binary_text<T: BitVect>(bv: &T) -> Vec<u8> {
    let n = bv.num_bits();
    let mut res = vec![0u8; (n as usize).div_ceil(8)];
    for b in bv.on_bits() {
        res[(b / 8) as usize] |= 1 << (b % 8);
    }
    res
}

/// `UpdateBitVectFromFPSText`
pub fn update_bit_vect_from_fps_text<T: BitVect>(bv: &mut T, fps: &[u8]) -> Result<()> {
    if (fps.len() as u64) * 4 < bv.num_bits() as u64 || fps.len() % 2 != 0 {
        return Err(RDKitError::precondition("bad FPS length"));
    }
    let mut bit_idx = 0u32;
    let mut i = 0;
    while i < fps.len() && bit_idx < bv.num_bits() {
        // strtol(.., 16) on the two characters: invalid digits end the number
        let mut c: u16 = 0;
        for &ch in &fps[i..i + 2] {
            match (ch as char).to_digit(16) {
                Some(d) => c = c * 16 + d as u16,
                None => break,
            }
        }
        for bit in 0..8 {
            if bit_idx >= bv.num_bits() {
                break;
            }
            if c & (1 << bit) != 0 {
                bv.set_bit(bit_idx)?;
            }
            bit_idx += 1;
        }
        i += 2;
    }
    Ok(())
}

/// `UpdateBitVectFromBinaryText`
pub fn update_bit_vect_from_binary_text<T: BitVect>(bv: &mut T, data: &[u8]) -> Result<()> {
    if (data.len() as u64) * 8 < bv.num_bits() as u64 {
        return Err(RDKitError::precondition("bad FPS length"));
    }
    let mut bit_idx = 0u32;
    for &c in data {
        if bit_idx >= bv.num_bits() {
            break;
        }
        // `unsigned short c = fps[i]` sign-extends a (signed) char
        let c = c as i8 as i16 as u16;
        for bit in 0..8 {
            if bit_idx >= bv.num_bits() {
                break;
            }
            if c & (1 << bit) != 0 {
                bv.set_bit(bit_idx)?;
            }
            bit_idx += 1;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ebv(n: u32, bits: &[u32]) -> ExplicitBitVect {
        let mut v = ExplicitBitVect::new(n, false);
        for &b in bits {
            v.set_bit(b).unwrap();
        }
        v
    }

    #[test]
    fn metrics() {
        let a = ebv(16, &[0, 1, 2, 3]);
        let b = ebv(16, &[2, 3, 4]);
        assert_eq!(similarity(Metric::Tanimoto, &a, &b).unwrap(), 2.0 / 5.0);
        assert_eq!(similarity(Metric::Dice, &a, &b).unwrap(), 4.0 / 7.0);
        assert_eq!(tversky(&a, &b, 1.0, 1.0).unwrap(), 2.0 / 5.0);
        assert_eq!(num_bits_in_common(&a, &b).unwrap(), 13);
        assert!(similarity(Metric::Tanimoto, &a, &ebv(8, &[])).is_err());
        assert_eq!(bit_vect_to_fps_text(&a), "0f00");
        let folded = fold_fingerprint(&ebv(16, &[1, 9]), 2).unwrap();
        assert_eq!(folded.on_bits(), vec![1]);
    }
}
