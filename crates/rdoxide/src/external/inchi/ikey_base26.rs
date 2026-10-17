//! Ported from InChI 1.07.3 `INCHI_BASE/src/ikey_base26.c` (InChIKey:
//! base-26 encoding of hash bits).
//!
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).
//!
//! The InChI build RDKit uses leaves `FIX_BASE26_ENC_BUG` undefined, so the
//! (historical, now normative) bit order of the `#ifndef` branches is ported.
//! Upstream tabulates the 16384 triplets as a literal table; it is generated
//! here from the rule the source comment states (all `AAA`..`ZZZ` in order,
//! without those starting with `E` and without `TAA`..`TTV`), which yields
//! the identical table.

const fn build_t26() -> [[u8; 3]; 16384] {
    let mut t = [[0u8; 3]; 16384];
    let mut n = 0;
    let mut a = 0u8;
    while a < 26 {
        let mut b = 0u8;
        while b < 26 {
            let mut c = 0u8;
            while c < 26 {
                // 'E'-started triplets and "TAA".."TTV" are intentionally omitted
                let skip = a == b'E' - b'A'
                    || (a == b'T' - b'A'
                        && (b < b'T' - b'A' || (b == b'T' - b'A' && c <= b'V' - b'A')));
                if !skip {
                    t[n] = [b'A' + a, b'A' + b, b'A' + c];
                    n += 1;
                }
                c += 1;
            }
            b += 1;
        }
        a += 1;
    }
    assert!(n == 16384);
    t
}

/// `t26`: 14-bit value -> triplet of uppercase letters.
static T26: [[u8; 3]; 16384] = build_t26();

/// `d26`: 9-bit value -> doublet ("AA".."ZZ" in order; only 512 are used).
fn d26(h: u32) -> [u8; 2] {
    [b'A' + (h / 26) as u8, b'A' + (h % 26) as u8]
}

/// `isbase26`
pub fn is_base26(c: u8) -> bool {
    c.is_ascii_uppercase()
}

/// `base26_triplet_1`: bits 0..13.
pub fn base26_triplet_1(a: &[u8; 32]) -> [u8; 3] {
    let b0 = a[0] as u32;
    let b1 = (a[1] & 0x3f) as u32;
    let h = b0 | b1 << 8;
    T26[h as usize]
}

/// `base26_triplet_2`: bits 14..27.
pub fn base26_triplet_2(a: &[u8; 32]) -> [u8; 3] {
    let b0 = (a[1] & 0xc0) as u32;
    let b1 = a[2] as u32;
    let b2 = (a[3] & 0x0f) as u32;
    let h = (b0 | b1 << 8 | b2 << 16) >> 6;
    T26[h as usize]
}

/// `base26_triplet_3`: bits 28..41.
pub fn base26_triplet_3(a: &[u8; 32]) -> [u8; 3] {
    let b0 = (a[3] & 0xf0) as u32;
    let b1 = a[4] as u32;
    let b2 = (a[5] & 0x03) as u32;
    let h = (b0 | b1 << 8 | b2 << 16) >> 4;
    T26[h as usize]
}

/// `base26_triplet_4`: bits 42..55.
pub fn base26_triplet_4(a: &[u8; 32]) -> [u8; 3] {
    let b0 = (a[5] & 0xfc) as u32;
    let b1 = a[6] as u32;
    let h = (b0 | b1 << 8) >> 2;
    T26[h as usize]
}

/// `base26_dublet_for_bits_28_to_36`
pub fn base26_dublet_for_bits_28_to_36(a: &[u8; 32]) -> [u8; 2] {
    let b0 = (a[3] & 0xf0) as u32;
    let b1 = (a[4] & 0x1f) as u32;
    let h = (b0 | b1 << 8) >> 4;
    d26(h)
}

/// `base26_dublet_for_bits_56_to_64`
pub fn base26_dublet_for_bits_56_to_64(a: &[u8; 32]) -> [u8; 2] {
    let b0 = a[7] as u32;
    let b1 = (a[8] & 0x01) as u32;
    let h = b0 | b1 << 8;
    d26(h)
}

fn hex_tail(first: u8, rest: &[u8]) -> String {
    let mut s = format!("{first:02x}");
    for b in rest {
        s.push_str(&format!("{b:02x}"));
    }
    s
}

/// `get_xtra_hash_major_hex`: the 191 hash bits not used by the major block.
pub fn get_xtra_hash_major_hex(a: &[u8; 32]) -> String {
    hex_tail(a[8] & 0xfe, &a[9..])
}

/// `get_xtra_hash_minor_hex`: the 219 hash bits not used by the minor block.
pub fn get_xtra_hash_minor_hex(a: &[u8; 32]) -> String {
    hex_tail(a[4] & 0xe0, &a[5..])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triplet_table_landmarks() {
        assert_eq!(&T26[0], b"AAA");
        assert_eq!(&T26[4 * 676 - 1], b"DZZ");
        assert_eq!(&T26[4 * 676], b"FAA");
        // 18 first letters (A..S without E) precede "TTW"
        assert_eq!(&T26[18 * 676 - 1], b"SZZ");
        assert_eq!(&T26[18 * 676], b"TTW");
        assert_eq!(&T26[16383], b"ZZZ");
    }

    #[test]
    fn doublets() {
        assert_eq!(d26(0), *b"AA");
        assert_eq!(d26(511), *b"TR");
    }
}
