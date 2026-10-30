//! Similarity measures on fingerprints stored as byte arrays (upstream
//! `CalcBitmap*` of `DataStructs/BitOps.cpp`).
//!
//! The first `n_bytes` of both slices are compared, and both must hold at
//! least that many bytes.

/// Eight bytes at a time as a little-endian word, so that the popcount of a fingerprint runs on whole words
/// (and vectorizes) instead of on single bytes.
fn word(c: &[u8]) -> u64 {
    u64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]])
}

fn popcount(bytes: &[u8]) -> u32 {
    let blocks = bytes.chunks_exact(32);
    let rest = blocks.remainder();
    let mut acc = [0u32; 4];
    for blk in blocks {
        for (k, a) in acc.iter_mut().enumerate() {
            *a += word(&blk[8 * k..]).count_ones();
        }
    }
    let words = rest.chunks_exact(8);
    let tail = words.remainder();
    acc.iter().sum::<u32>()
        + words.map(|c| word(c).count_ones()).sum::<u32>()
        + tail.iter().map(|b| b.count_ones()).sum::<u32>()
}

fn zip_popcount(a: &[u8], b: &[u8], op: impl Fn(u64, u64) -> u64 + Copy) -> u32 {
    let n = a.len().min(b.len());
    let (a, b) = (&a[..n], &b[..n]);
    let (ba, bb) = (a.chunks_exact(32), b.chunks_exact(32));
    let (ra, rb) = (ba.remainder(), bb.remainder());
    let mut acc = [0u32; 4];
    for (x, y) in ba.zip(bb) {
        for (k, s) in acc.iter_mut().enumerate() {
            *s += op(word(&x[8 * k..]), word(&y[8 * k..])).count_ones();
        }
    }
    let (wa, wb) = (ra.chunks_exact(8), rb.chunks_exact(8));
    let (ta, tb) = (wa.remainder(), wb.remainder());
    let mid: u32 = wa
        .zip(wb)
        .map(|(x, y)| op(word(x), word(y)).count_ones())
        .sum();
    let tail: u32 = ta
        .iter()
        .zip(tb)
        .map(|(&x, &y)| op(u64::from(x), u64::from(y)).count_ones())
        .sum();
    acc.iter().sum::<u32>() + mid + tail
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

#[cfg(test)]
mod word_tests {
    use super::*;

    #[test]
    fn word_kernels_match_bytewise_counts() {
        let mut seed = 12345u32;
        let mut next = || {
            seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
            (seed >> 24) as u8
        };
        for len in [0usize, 1, 7, 8, 9, 31, 32, 33, 100, 257] {
            let a: Vec<u8> = (0..len).map(|_| next()).collect();
            let b: Vec<u8> = (0..len).map(|_| next()).collect();
            let pc: u32 = a.iter().map(|x| x.count_ones()).sum();
            let and: u32 = a.iter().zip(&b).map(|(x, y)| (x & y).count_ones()).sum();
            let or: u32 = a.iter().zip(&b).map(|(x, y)| (x | y).count_ones()).sum();
            assert_eq!(bitmap_popcount(&a), pc, "len {len}");
            assert_eq!(bitmap_num_bits_in_common(&a, &b), and, "len {len}");
            assert_eq!(zip_popcount(&a, &b, |x, y| x | y), or, "len {len}");
        }
    }
}
