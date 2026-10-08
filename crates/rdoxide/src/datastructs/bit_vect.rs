//! Bit vectors (`ExplicitBitVect`, `SparseBitVect`), their binary pickle
//! format and base64 encoding.
//!
//! Ported from RDKit `Code/DataStructs/{BitVect,ExplicitBitVect,
//! SparseBitVect,base64}.cpp` and the packed-integer helpers of
//! `Code/RDGeneral/StreamOps.h` (2026.09.1).
//!
//! Copyright (C) 2002-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeSet;

use crate::{RDKitError, Result};

/// `ci_BITVECT_VERSION`
const BITVECT_VERSION: i32 = 0x0020;

fn index_error(which: u32) -> RDKitError {
    RDKitError::Index(format!("Index Error: {which}"))
}

/// Operations shared by both kinds of bit vectors (`BitVect`).
pub trait BitVect: Clone {
    fn new_sized(size: u32) -> Self;
    fn num_bits(&self) -> u32;
    fn num_on_bits(&self) -> u32;
    fn num_off_bits(&self) -> u32 {
        self.num_bits() - self.num_on_bits()
    }
    fn get_bit(&self, which: u32) -> Result<bool>;
    /// Returns the previous state of the bit.
    fn set_bit(&mut self, which: u32) -> Result<bool>;
    /// Returns the previous state of the bit.
    fn unset_bit(&mut self, which: u32) -> Result<bool>;
    fn on_bits(&self) -> Vec<i32>;
    fn and(&self, o: &Self) -> Self;
    fn or(&self, o: &Self) -> Self;
    fn xor(&self, o: &Self) -> Self;
    fn not(&self) -> Self;
    /// `toString()`: the binary pickle.
    fn to_binary(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(&(-BITVECT_VERSION).to_le_bytes());
        out.extend_from_slice(&(self.num_bits() as i32).to_le_bytes());
        out.extend_from_slice(&(self.num_on_bits() as i32).to_le_bytes());
        let mut prev: i64 = -1;
        for b in self.on_bits() {
            append_packed_int(&mut out, (b as i64 - prev - 1) as u32);
            prev = b as i64;
        }
        append_packed_int(&mut out, (self.num_bits() as i64 - prev - 1) as u32);
        out
    }
}

/// `appendPackedIntToStream`
pub fn append_packed_int(out: &mut Vec<u8>, num: u32) {
    let mut res = num;
    let (val, nbytes) = 'packed: {
        if res < (1 << 7) {
            break 'packed (res << 1, 1);
        }
        res -= 1 << 7;
        if res < (1 << 14) {
            break 'packed ((res << 2) | 1, 2);
        }
        res -= 1 << 14;
        if res < (1 << 21) {
            break 'packed ((res << 3) | 3, 3);
        }
        res -= 1 << 21;
        // upstream: CHECK_INVARIANT(res < 1 << 29)
        ((res << 3) | 7, 4)
    };
    let mut v = val;
    for _ in 0..nbytes {
        out.push((v & 255) as u8);
        v >>= 8;
    }
}

/// A byte cursor with `std::stringstream` read semantics.
struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl Reader<'_> {
    /// `streamRead`: a short read leaves the value untouched (0 here)
    fn i32(&mut self) -> i32 {
        let mut b = [0u8; 4];
        let n = (self.data.len() - self.pos.min(self.data.len())).min(4);
        b[..n].copy_from_slice(&self.data[self.pos..self.pos + n]);
        self.pos += 4;
        if n < 4 {
            return 0;
        }
        i32::from_le_bytes(b)
    }
    fn u16(&mut self) -> u16 {
        let mut b = [0u8; 2];
        let n = (self.data.len() - self.pos.min(self.data.len())).min(2);
        b[..n].copy_from_slice(&self.data[self.pos..self.pos + n]);
        self.pos += 2;
        if n < 2 {
            return 0;
        }
        u16::from_le_bytes(b)
    }
    fn byte(&mut self) -> Result<u32> {
        let b = self
            .data
            .get(self.pos)
            .copied()
            .ok_or_else(|| RDKitError::Runtime("failed to read from stream".into()))?;
        self.pos += 1;
        Ok(b as u32)
    }
    /// `readPackedIntFromStream`
    fn packed_int(&mut self) -> Result<u32> {
        let mut val = self.byte()?;
        let (shift, offset);
        if val & 1 == 0 {
            shift = 1;
            offset = 0;
        } else if val & 3 == 1 {
            val |= self.byte()? << 8;
            shift = 2;
            offset = 1 << 7;
        } else if val & 7 == 3 {
            val |= self.byte()? << 8;
            val |= self.byte()? << 16;
            shift = 3;
            offset = (1 << 7) + (1 << 14);
        } else {
            val |= self.byte()? << 8;
            val |= self.byte()? << 16;
            val |= self.byte()? << 24;
            shift = 3;
            offset = (1 << 7) + (1 << 14) + (1 << 21);
        }
        Ok((val >> shift).wrapping_add(offset))
    }
}

/// `BitVect::initFromText`
pub fn bit_vect_from_binary<T: BitVect>(data: &[u8], is_base64: bool) -> Result<T> {
    let decoded;
    let data = if is_base64 {
        decoded = base64_decode(data);
        &decoded[..]
    } else {
        data
    };
    let mut r = Reader { data, pos: 0 };
    let mut size = r.i32();
    let format = if size < 0 {
        let f = match -size {
            16 => 1,
            32 => 2,
            _ => return Err(RDKitError::Value("bad version in BitVect pickle".into())),
        };
        size = r.i32();
        f
    } else {
        return Err(RDKitError::Value("invalid BitVect pickle".into()));
    };
    let n_on = r.i32() as u32;
    let mut bv = T::new_sized(size as u32);
    if format == 1 && size >= u16::MAX as i32 {
        for _ in 0..n_on {
            let b = r.i32() as u32;
            bv.set_bit(b)?;
        }
    } else if format == 1 {
        for _ in 0..n_on {
            let b = r.u16() as u32;
            bv.set_bit(b)?;
        }
    } else {
        let mut curr: u32 = 0;
        for _ in 0..n_on {
            curr = curr.wrapping_add(r.packed_int()?);
            bv.set_bit(curr)?;
            curr = curr.wrapping_add(1);
        }
    }
    Ok(bv)
}

const B64: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/// `Base64Encode`
pub fn base64_encode(input: &[u8]) -> String {
    let mut res = Vec::with_capacity(input.len() * 4 / 3 + 4);
    let mut i = 0;
    while i < input.len() {
        res.push(B64[(input[i] >> 2) as usize]);
        if i + 1 < input.len() {
            res.push(B64[(((input[i] & 3) << 4) | (input[i + 1] >> 4)) as usize]);
            if i + 2 < input.len() {
                res.push(B64[(((input[i + 1] & 0xF) << 2) | (input[i + 2] >> 6)) as usize]);
                res.push(B64[(input[i + 2] & 0x3F) as usize]);
            } else {
                res.push(B64[((input[i + 1] & 0xF) << 2) as usize]);
                res.push(b'=');
            }
        } else {
            res.push(B64[((input[i] & 3) << 4) as usize]);
            res.push(b'=');
            res.push(b'=');
        }
        i += 3;
    }
    String::from_utf8(res).expect("ascii")
}

/// `Base64Decode` (stops at a NUL like `strlen`)
pub fn base64_decode(input: &[u8]) -> Vec<u8> {
    let input = match input.iter().position(|&c| c == 0) {
        Some(p) => &input[..p],
        None => input,
    };
    let tr = |c: u8| -> Option<u8> {
        match c {
            b'A'..=b'Z' => Some(c - b'A'),
            b'a'..=b'z' => Some(c - b'a' + 26),
            b'0'..=b'9' => Some(c - b'0' + 52),
            b'+' => Some(62),
            b'/' => Some(63),
            _ => None,
        }
    };
    let mut res = Vec::with_capacity(3 * input.len() / 4);
    let mut block = [0u8; 4];
    let mut n = 0;
    for &c in input {
        if let Some(v) = tr(c) {
            block[n] = v;
            n += 1;
            if n == 4 {
                res.push((block[0] << 2) | (block[1] >> 4));
                res.push((block[1] << 4) | (block[2] >> 2));
                res.push((block[2] << 6) | block[3]);
                n = 0;
            }
        }
    }
    if n > 1 {
        res.push((block[0] << 2) | (block[1] >> 4));
        if n > 2 {
            res.push((block[1] << 4) | (block[2] >> 2));
        }
    }
    res
}

/// `ExplicitBitVect` (a `boost::dynamic_bitset`)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct ExplicitBitVect {
    size: u32,
    words: Vec<u64>,
    num_on: u32,
}

impl ExplicitBitVect {
    pub fn new(size: u32, bits_set: bool) -> Self {
        let mut v = ExplicitBitVect {
            size,
            words: vec![0; (size as usize).div_ceil(64)],
            num_on: 0,
        };
        if bits_set {
            for w in v.words.iter_mut() {
                *w = u64::MAX;
            }
            v.trim();
            v.num_on = size;
        }
        v
    }

    fn trim(&mut self) {
        let r = self.size % 64;
        if r != 0
            && let Some(last) = self.words.last_mut()
        {
            *last &= (1u64 << r) - 1;
        }
    }

    fn count(&self) -> u32 {
        self.words.iter().map(|w| w.count_ones()).sum()
    }

    fn get(&self, i: u32) -> bool {
        self.words[(i / 64) as usize] >> (i % 64) & 1 == 1
    }

    fn zip(&self, o: &Self, f: impl Fn(u64, u64) -> u64) -> Self {
        let mut words: Vec<u64> = self
            .words
            .iter()
            .enumerate()
            .map(|(i, &a)| f(a, o.words.get(i).copied().unwrap_or(0)))
            .collect();
        words.truncate(self.words.len());
        let mut r = ExplicitBitVect {
            size: self.size,
            words,
            num_on: 0,
        };
        r.trim();
        r.num_on = r.count();
        r
    }

    /// `operator+=`: appends `o`'s bits.
    pub fn append(&mut self, o: &Self) {
        let orig = self.size;
        self.size += o.size;
        self.words.resize((self.size as usize).div_ceil(64), 0);
        for i in 0..o.size {
            if o.get(i) {
                let j = orig + i;
                self.words[(j / 64) as usize] |= 1u64 << (j % 64);
            }
        }
        self.num_on = self.count();
    }

    /// `AllProbeBitsMatch(ebv, ebv)`: `is_subset_of`
    pub fn is_subset_of(&self, o: &Self) -> bool {
        self.size == o.size && self.words.iter().zip(&o.words).all(|(a, b)| a & !b == 0)
    }

    /// popcount of `self & o` over the shared words
    pub fn num_on_in_common(&self, o: &Self) -> u32 {
        self.words
            .iter()
            .zip(&o.words)
            .map(|(a, b)| (a & b).count_ones())
            .sum()
    }

    /// popcount of `self ^ o`
    pub fn num_xor(&self, o: &Self) -> u32 {
        self.zip(o, |a, b| a ^ b).num_on
    }

    /// Raw little-endian bytes (`BitVectToBinaryText`).
    pub fn to_bytes(&self) -> Vec<u8> {
        let n = (self.size as usize).div_ceil(8);
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            out.push((self.words[i / 8] >> ((i % 8) * 8)) as u8);
        }
        out
    }
}

impl BitVect for ExplicitBitVect {
    fn new_sized(size: u32) -> Self {
        ExplicitBitVect::new(size, false)
    }
    fn num_bits(&self) -> u32 {
        self.size
    }
    fn num_on_bits(&self) -> u32 {
        self.num_on
    }
    fn get_bit(&self, which: u32) -> Result<bool> {
        if which >= self.size {
            return Err(index_error(which));
        }
        Ok(self.get(which))
    }
    fn set_bit(&mut self, which: u32) -> Result<bool> {
        if which >= self.size {
            return Err(index_error(which));
        }
        if self.get(which) {
            return Ok(true);
        }
        self.words[(which / 64) as usize] |= 1u64 << (which % 64);
        self.num_on += 1;
        Ok(false)
    }
    fn unset_bit(&mut self, which: u32) -> Result<bool> {
        if which >= self.size {
            return Err(index_error(which));
        }
        if !self.get(which) {
            return Ok(false);
        }
        self.words[(which / 64) as usize] &= !(1u64 << (which % 64));
        self.num_on -= 1;
        Ok(true)
    }
    fn on_bits(&self) -> Vec<i32> {
        let mut v = Vec::with_capacity(self.num_on as usize);
        for (wi, &w) in self.words.iter().enumerate() {
            let mut w = w;
            while w != 0 {
                let b = w.trailing_zeros();
                v.push((wi as u32 * 64 + b) as i32);
                w &= w - 1;
            }
        }
        v
    }
    fn and(&self, o: &Self) -> Self {
        self.zip(o, |a, b| a & b)
    }
    fn or(&self, o: &Self) -> Self {
        self.zip(o, |a, b| a | b)
    }
    fn xor(&self, o: &Self) -> Self {
        self.zip(o, |a, b| a ^ b)
    }
    fn not(&self) -> Self {
        let mut r = ExplicitBitVect {
            size: self.size,
            words: self.words.iter().map(|w| !w).collect(),
            num_on: 0,
        };
        r.trim();
        r.num_on = r.count();
        r
    }
}

/// `SparseBitVect` (a `std::set<int>`)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct SparseBitVect {
    size: u32,
    /// `std::set<int>`: iteration is in signed order
    bits: BTreeSet<i32>,
}

impl SparseBitVect {
    pub fn new(size: u32) -> Self {
        SparseBitVect {
            size,
            bits: BTreeSet::new(),
        }
    }
}

impl BitVect for SparseBitVect {
    fn new_sized(size: u32) -> Self {
        SparseBitVect::new(size)
    }
    fn num_bits(&self) -> u32 {
        self.size
    }
    fn num_on_bits(&self) -> u32 {
        self.bits.len() as u32
    }
    fn get_bit(&self, which: u32) -> Result<bool> {
        if which >= self.size {
            return Err(index_error(which));
        }
        Ok(self.bits.contains(&(which as i32)))
    }
    fn set_bit(&mut self, which: u32) -> Result<bool> {
        if which >= self.size {
            return Err(index_error(which));
        }
        Ok(!self.bits.insert(which as i32))
    }
    fn unset_bit(&mut self, which: u32) -> Result<bool> {
        if which >= self.size {
            return Err(index_error(which));
        }
        Ok(self.bits.remove(&(which as i32)))
    }
    fn on_bits(&self) -> Vec<i32> {
        self.bits.iter().copied().collect()
    }
    fn and(&self, o: &Self) -> Self {
        SparseBitVect {
            size: self.size,
            bits: self.bits.intersection(&o.bits).copied().collect(),
        }
    }
    fn or(&self, o: &Self) -> Self {
        SparseBitVect {
            size: self.size,
            bits: self.bits.union(&o.bits).copied().collect(),
        }
    }
    fn xor(&self, o: &Self) -> Self {
        SparseBitVect {
            size: self.size,
            bits: self.bits.symmetric_difference(&o.bits).copied().collect(),
        }
    }
    fn not(&self) -> Self {
        SparseBitVect {
            size: self.size,
            bits: (0..self.size)
                .map(|b| b as i32)
                .filter(|b| !self.bits.contains(b))
                .collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packed_ints_roundtrip() {
        for n in [
            0u32,
            1,
            127,
            128,
            16511,
            16512,
            2113663,
            2113664,
            100_000_000,
        ] {
            let mut v = Vec::new();
            append_packed_int(&mut v, n);
            let mut r = Reader { data: &v, pos: 0 };
            assert_eq!(r.packed_int().unwrap(), n);
        }
    }

    #[test]
    fn pickle_roundtrip() {
        let mut bv = ExplicitBitVect::new(2048, false);
        for b in [0, 5, 63, 64, 1000, 2047] {
            bv.set_bit(b).unwrap();
        }
        let p = bv.to_binary();
        let back: ExplicitBitVect = bit_vect_from_binary(&p, false).unwrap();
        assert_eq!(back, bv);
        let b64 = base64_encode(&p);
        let back: ExplicitBitVect = bit_vect_from_binary(b64.as_bytes(), true).unwrap();
        assert_eq!(back, bv);
        let mut sbv = SparseBitVect::new(1 << 20);
        sbv.set_bit(123456).unwrap();
        let back: SparseBitVect = bit_vect_from_binary(&sbv.to_binary(), false).unwrap();
        assert_eq!(back, sbv);
    }
}
