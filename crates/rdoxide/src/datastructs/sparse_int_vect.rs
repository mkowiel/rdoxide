//! `SparseIntVect<IndexType>` (upstream `DataStructs/SparseIntVect.h`).
//!
//! The C++ class is a template over the index type; here a single struct
//! carries the index kind and stores indices as `i128`, which holds every
//! value of all four instantiations with the same ordering.

use std::collections::BTreeMap;

use crate::error::{RDKitError, Result};
use crate::rdgeneral::fp::fmuladd;

const VERSION: u32 = 0x0001;

/// The `IndexType` template argument.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IndexKind {
    I32,
    I64,
    U32,
    U64,
}

impl IndexKind {
    pub fn size(self) -> u32 {
        match self {
            IndexKind::I32 | IndexKind::U32 => 4,
            IndexKind::I64 | IndexKind::U64 => 8,
        }
    }
    pub fn max(self) -> i128 {
        match self {
            IndexKind::I32 => i32::MAX as i128,
            IndexKind::I64 => i64::MAX as i128,
            IndexKind::U32 => u32::MAX as i128,
            IndexKind::U64 => u64::MAX as i128,
        }
    }
    pub fn min(self) -> i128 {
        match self {
            IndexKind::I32 => i32::MIN as i128,
            IndexKind::I64 => i64::MIN as i128,
            IndexKind::U32 | IndexKind::U64 => 0,
        }
    }
    /// Wraps a value read from a pickle into the index type (C++ integral conversion).
    fn wrap(self, v: u64) -> i128 {
        match self {
            IndexKind::I32 => v as u32 as i32 as i128,
            IndexKind::I64 => v as i64 as i128,
            IndexKind::U32 => v as u32 as i128,
            IndexKind::U64 => v as i128,
        }
    }
    fn write(self, out: &mut Vec<u8>, v: i128) {
        match self.size() {
            4 => out.extend_from_slice(&(v as u32).to_le_bytes()),
            _ => out.extend_from_slice(&(v as u64).to_le_bytes()),
        }
    }
}

/// Values may be stored as zeros after scalar arithmetic, exactly like upstream.
#[derive(Clone, Debug)]
pub struct SparseIntVect {
    kind: IndexKind,
    length: i128,
    data: BTreeMap<i128, i32>,
}

impl PartialEq for SparseIntVect {
    fn eq(&self, o: &Self) -> bool {
        self.length == o.length && self.data == o.data
    }
}

fn size_mismatch() -> RDKitError {
    RDKitError::Value("SparseIntVect size mismatch".into())
}

impl SparseIntVect {
    pub fn new(kind: IndexKind, length: i128) -> Self {
        SparseIntVect {
            kind,
            length,
            data: BTreeMap::new(),
        }
    }
    pub fn kind(&self) -> IndexKind {
        self.kind
    }
    pub fn length(&self) -> i128 {
        self.length
    }
    pub fn nonzero_elements(&self) -> &BTreeMap<i128, i32> {
        &self.data
    }

    fn check_index(&self, idx: i128) -> Result<()> {
        if idx < 0 || idx > self.length || (idx == self.length && self.length < self.kind.max()) {
            return Err(RDKitError::Index(format!("Index Error: {}", idx as i32)));
        }
        Ok(())
    }
    pub fn get_val(&self, idx: i128) -> Result<i32> {
        self.check_index(idx)?;
        Ok(self.data.get(&idx).copied().unwrap_or(0))
    }
    pub fn set_val(&mut self, idx: i128, val: i32) -> Result<()> {
        self.check_index(idx)?;
        if val != 0 {
            self.data.insert(idx, val);
        } else {
            self.data.remove(&idx);
        }
        Ok(())
    }
    /// `getTotalVal` (int arithmetic)
    pub fn total_val(&self, do_abs: bool) -> i32 {
        self.data.values().fold(0i32, |s, &v| {
            s.wrapping_add(if do_abs { v.wrapping_abs() } else { v })
        })
    }
    /// `updateFromSequence`
    pub fn update_from_sequence(&mut self, seq: &[i128]) -> Result<()> {
        for &i in seq {
            let v = self.get_val(i)?;
            self.set_val(i, v.wrapping_add(1))?;
        }
        Ok(())
    }

    fn same_len(&self, o: &Self) -> Result<()> {
        if self.length != o.length {
            Err(size_mismatch())
        } else {
            Ok(())
        }
    }
    /// `operator&=`: keeps keys present in both with the smaller value
    pub fn and_assign(&mut self, o: &Self) -> Result<()> {
        self.same_len(o)?;
        self.data.retain(|k, v| match o.data.get(k) {
            Some(&ov) => {
                if ov < *v {
                    *v = ov;
                }
                true
            }
            None => false,
        });
        Ok(())
    }
    /// `operator|=`: union keeping the larger value
    pub fn or_assign(&mut self, o: &Self) -> Result<()> {
        self.same_len(o)?;
        for (&k, &ov) in &o.data {
            match self.data.get_mut(&k) {
                Some(v) => {
                    if ov > *v {
                        *v = ov;
                    }
                }
                None => {
                    self.data.insert(k, ov);
                }
            }
        }
        Ok(())
    }
    fn add_sub(&mut self, o: &Self, sign: bool) -> Result<()> {
        self.same_len(o)?;
        for (&k, &ov) in &o.data {
            match self.data.get_mut(&k) {
                Some(v) => {
                    *v = if sign {
                        v.wrapping_add(ov)
                    } else {
                        v.wrapping_sub(ov)
                    };
                    if *v == 0 {
                        self.data.remove(&k);
                    }
                }
                None => {
                    self.data
                        .insert(k, if sign { ov } else { ov.wrapping_neg() });
                }
            }
        }
        Ok(())
    }
    pub fn add_assign(&mut self, o: &Self) -> Result<()> {
        self.add_sub(o, true)
    }
    pub fn sub_assign(&mut self, o: &Self) -> Result<()> {
        self.add_sub(o, false)
    }
    /// Scalar `+= -= *= /=`: applied to the stored entries only (zeros are kept).
    pub fn scalar_op(&mut self, op: char, s: i32) -> Result<()> {
        if op == '/' && s == 0 {
            return Err(RDKitError::Value("division by zero".into()));
        }
        for v in self.data.values_mut() {
            *v = match op {
                '+' => v.wrapping_add(s),
                '-' => v.wrapping_sub(s),
                '*' => v.wrapping_mul(s),
                _ => v.wrapping_div(s),
            };
        }
        Ok(())
    }

    /// `toString`
    pub fn to_binary(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(&VERSION.to_le_bytes());
        out.extend_from_slice(&self.kind.size().to_le_bytes());
        self.kind.write(&mut out, self.length);
        self.kind.write(&mut out, self.data.len() as i128);
        for (&k, &v) in &self.data {
            self.kind.write(&mut out, k);
            out.extend_from_slice(&v.to_le_bytes());
        }
        out
    }

    /// `initFromText`
    pub fn from_binary(kind: IndexKind, pkl: &[u8]) -> Result<Self> {
        let mut pos = 0usize;
        let mut read = |n: usize| -> Result<u64> {
            let b = pkl
                .get(pos..pos + n)
                .ok_or_else(|| RDKitError::Value("unreadable format".into()))?;
            pos += n;
            let mut buf = [0u8; 8];
            buf[..n].copy_from_slice(b);
            Ok(u64::from_le_bytes(buf))
        };
        let vers = read(4)? as u32;
        if vers != 0x0001 {
            return Err(RDKitError::Value(
                "bad version in SparseIntVect pickle".into(),
            ));
        }
        let t = read(4)? as u32;
        if t > kind.size() {
            return Err(RDKitError::Value(
                "IndexType cannot accommodate index size in SparseIntVect pickle".into(),
            ));
        }
        let n = match t {
            1 | 4 | 8 => t as usize,
            _ => return Err(RDKitError::Value("unreadable format".into())),
        };
        let mut res = SparseIntVect::new(kind, kind.wrap(read(n)?));
        let n_entries = read(n)?;
        for _ in 0..n_entries {
            let k = kind.wrap(read(n)?);
            let v = read(4)? as u32 as i32;
            res.data.insert(k, v);
        }
        Ok(res)
    }
}

/// `calcVectParams`: (sum |v1|, sum |v2|, sum over common keys of min(|a|,|b|))
fn calc_vect_params(v1: &SparseIntVect, v2: &SparseIntVect) -> Result<(f64, f64, f64)> {
    v1.same_len(v2)?;
    let mut v1_sum = 0.0;
    let mut v2_sum = 0.0;
    let mut and_sum = 0.0;
    let mut it1 = v1.data.iter().peekable();
    let mut it2 = v2.data.iter().peekable();
    if let Some(&(_, &v)) = it1.peek() {
        v1_sum += v.wrapping_abs() as f64;
    }
    if let Some(&(_, &v)) = it2.peek() {
        v2_sum += v.wrapping_abs() as f64;
    }
    while let Some(&(&k1, &a)) = it1.peek() {
        while let Some(&(&k2, _)) = it2.peek() {
            if k2 >= k1 {
                break;
            }
            it2.next();
            if let Some(&(_, &v)) = it2.peek() {
                v2_sum += v.wrapping_abs() as f64;
            }
        }
        match it2.peek() {
            Some(&(&k2, &b)) => {
                if k2 == k1 {
                    and_sum += if b.wrapping_abs() < a.wrapping_abs() {
                        b.wrapping_abs()
                    } else {
                        a.wrapping_abs()
                    } as f64;
                    it2.next();
                    if let Some(&(_, &v)) = it2.peek() {
                        v2_sum += v.wrapping_abs() as f64;
                    }
                }
                it1.next();
                if let Some(&(_, &v)) = it1.peek() {
                    v1_sum += v.wrapping_abs() as f64;
                }
            }
            None => break,
        }
    }
    if it1.next().is_some() {
        for (_, &v) in it1 {
            v1_sum += v.wrapping_abs() as f64;
        }
    }
    if it2.next().is_some() {
        for (_, &v) in it2 {
            v2_sum += v.wrapping_abs() as f64;
        }
    }
    Ok((v1_sum, v2_sum, and_sum))
}

/// `DiceSimilarity`
pub fn dice_similarity(
    v1: &SparseIntVect,
    v2: &SparseIntVect,
    return_distance: bool,
    bounds: f64,
) -> Result<f64> {
    v1.same_len(v2)?;
    if !return_distance && bounds > 0.0 {
        let s1 = v1.total_val(true) as f64;
        let s2 = v2.total_val(true) as f64;
        let denom = s1 + s2;
        if denom.abs() < 1e-6 {
            return Ok(0.0);
        }
        let min_v = if s1 < s2 { s1 } else { s2 };
        if 2. * min_v / denom < bounds {
            return Ok(0.0);
        }
    }
    let (s1, s2, numer) = calc_vect_params(v1, v2)?;
    let denom = s1 + s2;
    let mut sim = if denom.abs() < 1e-6 {
        0.0
    } else {
        2. * numer / denom
    };
    if return_distance {
        sim = 1. - sim;
    }
    Ok(sim)
}

/// `TverskySimilarity`
pub fn tversky_similarity(
    v1: &SparseIntVect,
    v2: &SparseIntVect,
    a: f64,
    b: f64,
    return_distance: bool,
) -> Result<f64> {
    let (s1, s2, and_sum) = calc_vect_params(v1, v2)?;
    // a*v1Sum + b*v2Sum + (1-a-b)*andSum
    let denom = fmuladd(1.0 - a - b, and_sum, fmuladd(a, s1, b * s2));
    let mut sim = if denom.abs() < 1e-6 {
        0.0
    } else {
        and_sum / denom
    };
    if return_distance {
        sim = 1. - sim;
    }
    Ok(sim)
}

/// `TanimotoSimilarity`
pub fn tanimoto_similarity(
    v1: &SparseIntVect,
    v2: &SparseIntVect,
    return_distance: bool,
) -> Result<f64> {
    tversky_similarity(v1, v2, 1.0, 1.0, return_distance)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dice_and_pickle() {
        let mut v1 = SparseIntVect::new(IndexKind::I32, 5);
        let mut v2 = SparseIntVect::new(IndexKind::I32, 5);
        for (i, v) in [(0, 2), (2, 1), (3, 4), (4, 6)] {
            v1.set_val(i, v).unwrap();
        }
        for (i, v) in [(1, 2), (2, 3), (3, 4), (4, 4)] {
            v2.set_val(i, v).unwrap();
        }
        assert_eq!(dice_similarity(&v1, &v2, false, 0.0).unwrap(), 18.0 / 26.0);
        let back = SparseIntVect::from_binary(IndexKind::I32, &v1.to_binary()).unwrap();
        assert_eq!(back, v1);
        assert!(v1.get_val(5).is_err());
    }
}
