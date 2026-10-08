//! `DiscreteValueVect` and `RealValueVect` (upstream
//! `DataStructs/{DiscreteValueVect,RealValueVect}.{h,cpp}`).

use crate::error::{RDKitError, Result};

const BITS_PER_INT: u32 = 32;

fn index_error(i: u32) -> RDKitError {
    RDKitError::Index(format!("Index Error: {}", i as i32))
}

fn read_u32(pkl: &[u8], pos: &mut usize) -> u32 {
    // a short read leaves the (zero-initialized) value untouched
    let mut b = [0u8; 4];
    if let Some(s) = pkl.get(*pos..*pos + 4) {
        b.copy_from_slice(s);
    }
    *pos += 4;
    u32::from_le_bytes(b)
}

/// A vector of small unsigned values packed into 32-bit words.
#[derive(Clone, Debug, PartialEq)]
pub struct DiscreteValueVect {
    /// `DiscreteValueType`: 0 = ONEBITVALUE .. 4 = SIXTEENBITVALUE
    val_type: u32,
    bits_per_val: u32,
    vals_per_int: u32,
    length: u32,
    mask: u32,
    data: Vec<u32>,
}

impl DiscreteValueVect {
    pub fn new(val_type: u32, length: u32) -> Self {
        let bits_per_val = 1u32 << val_type;
        let vals_per_int = BITS_PER_INT / bits_per_val;
        let num_ints = length.div_ceil(vals_per_int);
        DiscreteValueVect {
            val_type,
            bits_per_val,
            vals_per_int,
            length,
            mask: ((1u64 << bits_per_val) - 1) as u32,
            data: vec![0; num_ints as usize],
        }
    }
    pub fn length(&self) -> u32 {
        self.length
    }
    pub fn value_type(&self) -> u32 {
        self.val_type
    }
    pub fn get_val(&self, i: u32) -> Result<u32> {
        if i >= self.length {
            return Err(index_error(i));
        }
        let shift = self.bits_per_val * (i % self.vals_per_int);
        Ok((self.data[(i / self.vals_per_int) as usize] >> shift) & self.mask)
    }
    pub fn set_val(&mut self, i: u32, val: u32) -> Result<()> {
        if i >= self.length {
            return Err(index_error(i));
        }
        if val & self.mask != val {
            return Err(RDKitError::Value("Value out of range".into()));
        }
        let shift = self.bits_per_val * (i % self.vals_per_int);
        let id = (i / self.vals_per_int) as usize;
        let mask = !(self.mask << shift);
        self.data[id] = (self.data[id] & mask) | (val << shift);
        Ok(())
    }
    /// `getTotalVal` (sums every slot of the storage, unsigned arithmetic)
    pub fn total_val(&self) -> u32 {
        let mut res = 0u32;
        for &w in &self.data {
            for j in 0..self.vals_per_int {
                res =
                    res.wrapping_add(w.checked_shr(j * self.bits_per_val).unwrap_or(0) & self.mask);
            }
        }
        res
    }
    fn check_len(&self, o: &Self) -> Result<()> {
        if o.length != self.length {
            return Err(RDKitError::precondition("length mismatch"));
        }
        Ok(())
    }
    fn combine(&self, o: &Self, typ: u32, f: impl Fn(u32, u32) -> u32) -> Result<Self> {
        self.check_len(o)?;
        let mut ans = DiscreteValueVect::new(typ, self.length);
        for i in 0..self.length {
            ans.set_val(i, f(self.get_val(i)?, o.get_val(i)?))?;
        }
        Ok(ans)
    }
    /// `operator&`: element-wise minimum, narrower value type
    pub fn and(&self, o: &Self) -> Result<Self> {
        self.combine(o, self.val_type.min(o.val_type), u32::min)
    }
    /// `operator|`: element-wise maximum, wider value type
    pub fn or(&self, o: &Self) -> Result<Self> {
        self.combine(o, self.val_type.max(o.val_type), u32::max)
    }
    /// `operator+=`: values are truncated at the maximum
    pub fn add_assign(&mut self, o: &Self) -> Result<()> {
        self.check_len(o)?;
        let max_val = self.mask;
        for i in 0..self.length {
            let v = self.get_val(i)?.wrapping_add(o.get_val(i)?);
            self.set_val(i, v.min(max_val))?;
        }
        Ok(())
    }
    /// `operator-=`: would-be negative values become zero
    pub fn sub_assign(&mut self, o: &Self) -> Result<()> {
        self.check_len(o)?;
        for i in 0..self.length {
            let (a, b) = (self.get_val(i)?, o.get_val(i)?);
            self.set_val(i, a.saturating_sub(b))?;
        }
        Ok(())
    }

    /// `toString`
    pub fn to_binary(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(24 + 4 * self.data.len());
        for v in [
            (-1i32) as u32,
            self.val_type,
            self.bits_per_val,
            self.mask,
            self.length,
            self.data.len() as u32,
        ] {
            out.extend_from_slice(&v.to_le_bytes());
        }
        for w in &self.data {
            out.extend_from_slice(&w.to_le_bytes());
        }
        out
    }
    /// `initFromText`
    pub fn from_binary(pkl: &[u8]) -> Result<Self> {
        let mut pos = 0;
        let vers = (read_u32(pkl, &mut pos) as i32).wrapping_neg();
        if vers != 0x1 {
            return Err(RDKitError::Value(
                "bad version in DiscreteValueVect pickle".into(),
            ));
        }
        let val_type = read_u32(pkl, &mut pos);
        let bits_per_val = read_u32(pkl, &mut pos);
        if bits_per_val == 0 {
            return Err(RDKitError::Value("bad DiscreteValueVect pickle".into()));
        }
        let mask = read_u32(pkl, &mut pos);
        let length = read_u32(pkl, &mut pos);
        let num_ints = read_u32(pkl, &mut pos);
        let data = (0..num_ints).map(|_| read_u32(pkl, &mut pos)).collect();
        Ok(DiscreteValueVect {
            val_type,
            bits_per_val,
            vals_per_int: BITS_PER_INT / bits_per_val,
            length,
            mask,
            data,
        })
    }
}

/// `computeL1Norm` for DiscreteValueVects
pub fn compute_l1_norm_discrete(v1: &DiscreteValueVect, v2: &DiscreteValueVect) -> Result<u32> {
    if v1.length != v2.length {
        return Err(RDKitError::Value(
            "Comparing vectors of different lengths".into(),
        ));
    }
    if v1.val_type != v2.val_type {
        return Err(RDKitError::Value(
            "Comparing vector of different value types".into(),
        ));
    }
    // upstream walks the raw storage byte- (or short-) wise; summing the
    // differences of every slot is the same
    let mut res = 0u32;
    for (&a, &b) in v1.data.iter().zip(&v2.data) {
        for j in 0..v1.vals_per_int {
            let s = j * v1.bits_per_val;
            let x = a.checked_shr(s).unwrap_or(0) & v1.mask;
            let y = b.checked_shr(s).unwrap_or(0) & v1.mask;
            res = res.wrapping_add(x.abs_diff(y));
        }
    }
    Ok(res)
}

/// A vector of doubles.
#[derive(Clone, Debug, PartialEq)]
pub struct RealValueVect {
    data: Vec<f64>,
}

impl RealValueVect {
    pub fn new(length: u32) -> Self {
        RealValueVect {
            data: vec![0.0; length as usize],
        }
    }
    pub fn length(&self) -> u32 {
        self.data.len() as u32
    }
    pub fn data(&self) -> &[f64] {
        &self.data
    }
    pub fn get_val(&self, i: u32) -> Result<f64> {
        self.data
            .get(i as usize)
            .copied()
            .ok_or_else(|| index_error(i))
    }
    pub fn set_val(&mut self, i: u32, v: f64) -> Result<()> {
        *self
            .data
            .get_mut(i as usize)
            .ok_or_else(|| index_error(i))? = v;
        Ok(())
    }
    pub fn total_val(&self) -> f64 {
        self.data.iter().fold(0.0, |s, &v| s + v)
    }
    /// `applyBinaryOp`
    pub fn apply(&mut self, o: &Self, f: impl Fn(f64, f64) -> f64) -> Result<()> {
        if o.data.len() != self.data.len() {
            return Err(RDKitError::precondition("length mismatch"));
        }
        for (a, &b) in self.data.iter_mut().zip(&o.data) {
            *a = f(*a, b);
        }
        Ok(())
    }
    pub fn to_binary(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(8 + 8 * self.data.len());
        out.extend_from_slice(&(-1i32).to_le_bytes());
        out.extend_from_slice(&(self.data.len() as u32).to_le_bytes());
        for v in &self.data {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out
    }
    pub fn from_binary(pkl: &[u8]) -> Result<Self> {
        let mut pos = 0;
        let vers = (read_u32(pkl, &mut pos) as i32).wrapping_neg();
        if vers != 0x1 {
            return Err(RDKitError::Value(
                "bad version in RealValueVect pickle".into(),
            ));
        }
        let n = read_u32(pkl, &mut pos) as usize;
        let mut data = vec![0.0; n];
        for (i, d) in data.iter_mut().enumerate() {
            if let Some(s) = pkl.get(pos + 8 * i..pos + 8 * i + 8) {
                *d = f64::from_le_bytes(s.try_into().expect("8 bytes"));
            }
        }
        Ok(RealValueVect { data })
    }
}

/// `std::min<double>(v1, v2)` as used by `operator&=` (returns v1 unless v2 < v1)
pub fn cxx_min(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}
/// `std::max<double>(v1, v2)` (returns v1 unless v1 < v2)
pub fn cxx_max(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}

/// `computeL1Norm` for RealValueVects
pub fn compute_l1_norm_real(v1: &RealValueVect, v2: &RealValueVect) -> Result<f64> {
    if v1.data.len() != v2.data.len() {
        return Err(RDKitError::Value(
            "Comparing vectors of different lengths".into(),
        ));
    }
    Ok(v1
        .data
        .iter()
        .zip(&v2.data)
        .fold(0.0, |s, (a, b)| s + (a - b).abs()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn discrete_basics() {
        let mut v = DiscreteValueVect::new(2, 30);
        for i in 0..30 {
            v.set_val(i, i % 16).unwrap();
        }
        assert_eq!(v.get_val(17).unwrap(), 1);
        assert!(v.set_val(0, 16).is_err());
        let w = DiscreteValueVect::from_binary(&v.to_binary()).unwrap();
        assert_eq!(w, v);
        assert_eq!(
            compute_l1_norm_discrete(&v, &DiscreteValueVect::new(2, 30)).unwrap(),
            v.total_val()
        );
    }

    #[test]
    fn real_basics() {
        let mut v = RealValueVect::new(3);
        v.set_val(1, 2.5).unwrap();
        assert_eq!(RealValueVect::from_binary(&v.to_binary()).unwrap(), v);
        assert_eq!(v.total_val(), 2.5);
    }
}
