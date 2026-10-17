//! Ported from InChI 1.07.3 `INCHI_BASE/src/ikey_dll.c` (InChIKey: hash of
//! an InChI string) and `extract_inchi_substring` of `INCHI_BASE/src/util.c`.
//!
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).
//!
//! The SHA-256 used here is rdoxide's own ([`super::sha256`]), not a port of
//! InChI's `sha2.c`; SHA-256 is fully specified, so the digests are identical.
//!
//! The C code works on NUL-terminated buffers and in a few places copies the
//! terminating NUL as part of a block; the byte-level bookkeeping below keeps
//! the upstream indices and clips such copies at the end of the string, which
//! is what the NUL does in C.

use super::ikey_base26::*;
use super::sha256::sha256;
use super::{INCHI_STRING_PREFIX, LEN_INCHI_STRING_PREFIX};

/// `INCHIKEY_*` return codes of `GetINCHIKeyFromINCHI` (`inchi_api.h`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum InchiKeyError {
    /// `INCHIKEY_UNKNOWN_ERROR`
    UnknownError = 1,
    /// `INCHIKEY_EMPTY_INPUT`
    EmptyInput = 2,
    /// `INCHIKEY_INVALID_INCHI_PREFIX`
    InvalidInchiPrefix = 3,
    /// `INCHIKEY_NOT_ENOUGH_MEMORY`
    NotEnoughMemory = 4,
    /// `INCHIKEY_INVALID_INCHI`
    InvalidInchi = 20,
    /// `INCHIKEY_INVALID_STD_INCHI`
    InvalidStdInchi = 21,
}

/// Return codes of `CheckINCHIKey` (`inchi_api.h`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum InchiKeyCheck {
    /// `INCHIKEY_VALID_STANDARD`
    ValidStandard = 0,
    /// `INCHIKEY_VALID_NON_STANDARD`
    ValidNonStandard = -1,
    /// `INCHIKEY_INVALID_LENGTH`
    InvalidLength = 1,
    /// `INCHIKEY_INVALID_LAYOUT`
    InvalidLayout = 2,
    /// `INCHIKEY_INVALID_VERSION`
    InvalidVersion = 3,
}

/// Result of [`get_inchikey_from_inchi`]: the key and, when requested, the
/// hexadecimal hash extensions (`szXtra1`, `szXtra2`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InchiKey {
    pub key: String,
    pub xtra1: Option<String>,
    pub xtra2: Option<String>,
}

/// The bytes of a C string: everything before the first NUL.
fn c_bytes(s: &str) -> &[u8] {
    let b = s.as_bytes();
    match b.iter().position(|&c| c == 0) {
        Some(n) => &b[..n],
        None => b,
    }
}

/// `str[i]` of a NUL-terminated C string (0 past the end).
fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

/// `s[from..from + n]` clipped at the terminating NUL.
fn c_slice(s: &[u8], from: usize, n: usize) -> &[u8] {
    let from = from.min(s.len());
    &s[from..(from + n).min(s.len())]
}

/// `extract_inchi_substring` (util.c): the InChI starting at the first
/// `"InChI="`, up to the first character that cannot occur in an InChI
/// (treated as whitespace ending the string). `None` without `"InChI="`.
pub fn extract_inchi_substring(s: &[u8]) -> Option<&[u8]> {
    let p = s
        .windows(LEN_INCHI_STRING_PREFIX)
        .position(|w| w == INCHI_STRING_PREFIX)?;
    let s = &s[p..];
    let n = s
        .iter()
        .position(|&c| !(c.is_ascii_alphanumeric() || b"()*+,-./;=?@".contains(&c)))
        .unwrap_or(s.len());
    Some(&s[..n])
}

/// C `strtol(s, NULL, 10)` as stored into an `int` (two's-complement
/// truncation of the 64-bit `long`, saturated on overflow like `strtol`).
fn strtol_to_int(s: &[u8]) -> i32 {
    let mut i = 0;
    while i < s.len() && (s[i] == b' ' || (b'\t'..=b'\r').contains(&s[i])) {
        i += 1;
    }
    let mut neg = false;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') {
        neg = s[i] == b'-';
        i += 1;
    }
    let mut v: i64 = 0;
    let mut overflow = false;
    while i < s.len() && s[i].is_ascii_digit() {
        let d = (s[i] - b'0') as i64;
        match v.checked_mul(10).and_then(|x| {
            if neg {
                x.checked_sub(d)
            } else {
                x.checked_add(d)
            }
        }) {
            Some(x) => v = x,
            None => overflow = true,
        }
        i += 1;
    }
    if overflow {
        v = if neg { i64::MIN } else { i64::MAX };
    }
    v as i32
}

/// `GetINCHIKeyFromINCHI`: the InChIKey of an InChI string, optionally with
/// the hash extensions (`xtra1`/`xtra2`).
pub fn get_inchikey_from_inchi(
    inchi_source: &str,
    xtra1: bool,
    xtra2: bool,
) -> Result<InchiKey, InchiKeyError> {
    const FLAG_STD: u8 = b'S'; // standard key
    const FLAG_NONSTD: u8 = b'N'; // non-standard key
    const FLAG_EXPTL: u8 = b'B'; // experimental ('beta') key
    const FLAG_VER: u8 = b'A'; // InChI v. 1
    // Protonization encoding:
    // N 0
    // O +1 P +2 Q +3 R +4 S +5 T +6 U +7 V +8 W +9 X +10 Y +11 Z +12
    // M -1 L-2 K -3 J -4 I -5 H -6 G -7 F -8 E -9 D -10 C -11 B -12
    // A < -12 or > +12
    const PPLUS: &[u8] = b"OPQRSTUVWXYZ";
    const PMINUS: &[u8] = b"MLKJIHGFEDCB";
    let mut flagproto = b'N'; // no [de]protonization, by default

    let src = c_bytes(inchi_source);
    let slen = src.len();

    // .. has valid prefix
    if slen < LEN_INCHI_STRING_PREFIX + 3 {
        return Err(InchiKeyError::InvalidInchiPrefix);
    }
    if &src[..LEN_INCHI_STRING_PREFIX] != INCHI_STRING_PREFIX {
        return Err(InchiKeyError::InvalidInchiPrefix);
    }
    // .. has InChI version 1
    if src[LEN_INCHI_STRING_PREFIX] != b'1' {
        return Err(InchiKeyError::InvalidInchiPrefix);
    }
    // .. optionally has a 'standard' flag character
    // 0 non-standard, 1 standard, -1 experimental ('beta')
    let mut is_stdinchi = 0;
    let mut pos_slash1 = LEN_INCHI_STRING_PREFIX + 1;
    if src[pos_slash1] == b'S' {
        is_stdinchi = 1;
        pos_slash1 += 1;
    } else if src[pos_slash1] == b'B' {
        is_stdinchi = -1;
        pos_slash1 += 1;
    }
    // .. has trailing slash in the right place
    if at(src, pos_slash1) != b'/' {
        return Err(InchiKeyError::InvalidInchiPrefix);
    }
    // .. the rest of source string contains at least one a..Z0.9 or slash
    let next = at(src, pos_slash1 + 1);
    if !next.is_ascii_alphanumeric() && next != b'/' && next != b'?' {
        return Err(InchiKeyError::InvalidInchi);
    }

    // Ok. Will use a local copy of the source.
    let s = extract_inchi_substring(src).ok_or(InchiKeyError::NotEnoughMemory)?;
    let slen = s.len();

    // Extract the major block
    let mut jproto = 0usize;
    let mut j = pos_slash1 + 1;
    while j < slen - 1 {
        if s[j] == b'/' {
            match s[j + 1] {
                // anything allowed from a major part
                b'c' | b'h' | b'q' => {
                    j += 1;
                    continue;
                }
                // "/p"; protons now go to to special string, not to minor hash
                b'p' => {
                    jproto = j;
                    j += 1;
                    continue;
                }
                // "/f",  "/r" : may not occur in stdInChI
                b'f' | b'r' if is_stdinchi == 1 => {
                    return Err(InchiKeyError::InvalidStdInchi);
                }
                // anything allowed from a minor part
                _ => {}
            }
            break;
        }
        j += 1;
    }
    j += 1;
    if j == slen {
        j += 1;
    } else {
        j -= 1;
    }

    let ncp = if jproto != 0 {
        jproto - pos_slash1 - 1
    } else {
        j - pos_slash1 - 1
    };
    // Trim 'InChI=1[S]/'
    let smajor = c_slice(s, pos_slash1 + 1, ncp);

    // Treat protonization
    if jproto != 0 {
        let lenproto = j - jproto;
        if lenproto < 3 {
            // empty "/p", should not occur
            return Err(InchiKeyError::InvalidInchi);
        }
        let sproto = c_slice(s, pos_slash1 + ncp + 1, lenproto);
        let nprotons = strtol_to_int(&sproto[2.min(sproto.len())..]);
        if nprotons > 0 {
            flagproto = if nprotons > 12 {
                b'A'
            } else {
                PPLUS[(nprotons - 1) as usize]
            };
        } else if nprotons < 0 {
            flagproto = if nprotons < -12 {
                b'A'
            } else {
                PMINUS[(-nprotons - 1) as usize]
            };
        } else {
            // should never occur
            return Err(InchiKeyError::InvalidStdInchi);
        }
    }

    // Extract the minor block.
    let mut sminor: Vec<u8> = if j != slen + 1 {
        s[j..].to_vec()
    } else {
        Vec::new()
    };

    // Compute and compose the InChIKey string.
    let digest_major = sha256(smajor);
    let mut key = Vec::with_capacity(27);
    key.extend_from_slice(&base26_triplet_1(&digest_major));
    key.extend_from_slice(&base26_triplet_2(&digest_major));
    key.extend_from_slice(&base26_triplet_3(&digest_major));
    key.extend_from_slice(&base26_triplet_4(&digest_major));
    key.extend_from_slice(&base26_dublet_for_bits_56_to_64(&digest_major));

    // Minor hash sub-string: short minor blocks are hashed twice over.
    let minor_len = sminor.len();
    if minor_len > 0 && minor_len < 255 {
        sminor.extend_from_within(..);
    }
    let digest_minor = sha256(&sminor);
    key.push(b'-');
    key.extend_from_slice(&base26_triplet_1(&digest_minor));
    key.extend_from_slice(&base26_triplet_2(&digest_minor));
    key.extend_from_slice(&base26_dublet_for_bits_28_to_36(&digest_minor));
    // Append a standard/non-standard flag
    key.push(match is_stdinchi {
        1 => FLAG_STD,
        -1 => FLAG_EXPTL,
        _ => FLAG_NONSTD,
    });
    // Append InChI v.1 flag, dash and protonization flag
    key.push(FLAG_VER);
    key.push(b'-');
    key.push(flagproto);

    Ok(InchiKey {
        // only ASCII uppercase letters and '-' were pushed
        key: key.into_iter().map(char::from).collect(),
        xtra1: xtra1.then(|| get_xtra_hash_major_hex(&digest_major)),
        xtra2: xtra2.then(|| get_xtra_hash_minor_hex(&digest_minor)),
    })
}

/// `GetStdINCHIKeyFromStdINCHI`
pub fn get_std_inchikey_from_std_inchi(inchi_source: &str) -> Result<String, InchiKeyError> {
    let src = c_bytes(inchi_source);
    if src.len() < LEN_INCHI_STRING_PREFIX + 3 || src[LEN_INCHI_STRING_PREFIX + 1] != b'S' {
        return Err(InchiKeyError::InvalidStdInchi);
    }
    get_inchikey_from_inchi(inchi_source, false, false).map(|k| k.key)
}

/// `CheckINCHIKey`: is the string a well-formed InChIKey?
pub fn check_inchikey(inchikey: &str) -> InchiKeyCheck {
    let k = c_bytes(inchikey);
    // Proper length is 27
    if k.len() != 27 {
        return InchiKeyCheck::InvalidLength;
    }
    // Should have dash in 14-th and 25-th position
    if k[14] != b'-' || k[25] != b'-' {
        return InchiKeyCheck::InvalidLayout;
    }
    // All other should be uppercase
    if !k[..14].iter().all(|&c| is_base26(c))
        || !k[15..25].iter().all(|&c| is_base26(c))
        || !is_base26(k[26])
    {
        return InchiKeyCheck::InvalidLayout;
    }
    // No 'E' may appear in 0,3,6,and 9 positions of the 1st block ...
    // ... and 0 and 3 pos. of the second block.
    if [0, 3, 6, 9, 15, 18].iter().any(|&j| k[j] == b'E') {
        return InchiKeyCheck::InvalidLayout;
    }
    // Check for version (only 1 allowed)
    if k[24] != b'A' {
        return InchiKeyCheck::InvalidVersion;
    }
    // Check for standard-ness
    match k[23] {
        b'S' => InchiKeyCheck::ValidStandard,
        b'N' => InchiKeyCheck::ValidNonStandard,
        _ => InchiKeyCheck::InvalidLayout,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strtol_like_c() {
        assert_eq!(strtol_to_int(b"+1"), 1);
        assert_eq!(strtol_to_int(b"-13"), -13);
        assert_eq!(strtol_to_int(b"x"), 0);
        assert_eq!(strtol_to_int(b" 7/b"), 7);
    }

    #[test]
    fn substring_stops_at_whitespace() {
        assert_eq!(
            extract_inchi_substring(b"InChI=1/Ar\r\t"),
            Some(&b"InChI=1/Ar"[..])
        );
        assert_eq!(
            extract_inchi_substring(b"x InChI=1/Ar%"),
            Some(&b"InChI=1/Ar"[..])
        );
        assert_eq!(extract_inchi_substring(b"nope"), None);
    }

    #[test]
    fn check_key() {
        assert_eq!(
            check_inchikey("QTBSBXVTEAMEQO-UHFFFAOYSA-M"),
            InchiKeyCheck::ValidStandard
        );
        assert_eq!(
            check_inchikey("QTBSBXVTEAMEQO-UHFFFAOYNA-M"),
            InchiKeyCheck::ValidNonStandard
        );
        assert_eq!(
            check_inchikey("QTBSBXVTEAMEQO"),
            InchiKeyCheck::InvalidLength
        );
        assert_eq!(
            check_inchikey("ETBSBXVTEAMEQO-UHFFFAOYSA-M"),
            InchiKeyCheck::InvalidLayout
        );
        assert_eq!(
            check_inchikey("QTBSBXVTEAMEQO-UHFFFAOYSB-M"),
            InchiKeyCheck::InvalidVersion
        );
    }

    #[test]
    fn xtra_hashes() {
        let k = get_inchikey_from_inchi("InChI=1S/CH4/h1H4", true, true).unwrap();
        assert_eq!(k.key, "VNWKTOKETHGBQD-UHFFFAOYSA-N");
        assert_eq!(k.xtra1.unwrap().len(), 48);
        assert_eq!(k.xtra2.unwrap().len(), 56);
    }
}
