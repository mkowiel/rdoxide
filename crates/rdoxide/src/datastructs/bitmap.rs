//! Similarity measures on fingerprints stored as byte arrays (upstream
//! `CalcBitmap*` of `DataStructs/BitOps.cpp`).
//!
//! The first `n_bytes` of both slices are compared, and both must hold at
//! least that many bytes.

fn popcount(bytes: &[u8]) -> u32 {
    bytes.iter().map(|b| b.count_ones()).sum()
}

fn zip_popcount(a: &[u8], b: &[u8], op: impl Fn(u8, u8) -> u8) -> u32 {
    a.iter().zip(b).map(|(&x, &y)| op(x, y).count_ones()).sum()
}

/// `CalcBitmapPopcount`
pub fn bitmap_popcount(fp: &[u8]) -> u32 {
    popcount(fp)
}

/// `CalcBitmapNumBitsInCommon`
pub fn bitmap_num_bits_in_common(a: &[u8], b: &[u8]) -> u32 {
    zip_popcount(a, b, |x, y| x & y)
}

/// `CalcBitmapTanimoto`
pub fn bitmap_tanimoto(a: &[u8], b: &[u8]) -> f64 {
    let union = zip_popcount(a, b, |x, y| x | y);
    if union == 0 {
        return 0.0;
    }
    f64::from(bitmap_num_bits_in_common(a, b)) / f64::from(union)
}

/// `CalcBitmapTversky`
pub fn bitmap_tversky(a: &[u8], b: &[u8], ca: f64, cb: f64) -> f64 {
    let common = f64::from(bitmap_num_bits_in_common(a, b));
    let denom =
        ca * f64::from(popcount(a)) + cb * f64::from(popcount(b)) + (1.0 - ca - cb) * common;
    if denom == 0.0 { 0.0 } else { common / denom }
}

/// `CalcBitmapAllProbeBitsMatch`: whether every bit of `probe` is set in `reference`.
pub fn bitmap_all_probe_bits_match(probe: &[u8], reference: &[u8]) -> bool {
    probe.iter().zip(reference).all(|(&p, &r)| p & r == p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn similarities() {
        let a = [0b1111_0000u8, 0b0000_0001];
        let b = [0b1100_0000u8, 0b0000_0011];
        assert_eq!(bitmap_popcount(&a), 5);
        assert_eq!(bitmap_num_bits_in_common(&a, &b), 3);
        assert!((bitmap_tanimoto(&a, &b) - 3.0 / 6.0).abs() < 1e-12);
        assert!((bitmap_tversky(&a, &b, 1.0, 1.0) - 3.0 / 6.0).abs() < 1e-12);
        assert_eq!(bitmap_tanimoto(&[0], &[0]), 0.0);
        assert!(bitmap_all_probe_bits_match(&[0b0100_0000], &b));
        assert!(!bitmap_all_probe_bits_match(&a, &b));
    }
}
