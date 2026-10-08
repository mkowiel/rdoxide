//! Fixed-column field parsing helpers (`FileParserUtils`) and an
//! `std::getline`-compatible line reader.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/MolFileParser.cpp`,
//! `FileParserUtils.h` and `RDGeneral/StreamOps.h` (2026.09.1).
//!
//! Copyright (C) 2002-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::{RDKitError, Result};

/// A byte stream read line by line with `std::getline` semantics: `eof` is
/// set when a read hits the end of the data (also for a final line without
/// a newline), `fail` when nothing at all could be read.
#[derive(Debug, Clone)]
pub struct LineStream<'a> {
    data: &'a [u8],
    pos: usize,
    eof: bool,
    fail: bool,
}

impl<'a> LineStream<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        LineStream {
            data,
            pos: 0,
            eof: false,
            fail: false,
        }
    }
    pub fn eof(&self) -> bool {
        self.eof
    }
    pub fn fail(&self) -> bool {
        self.fail
    }
    pub fn position(&self) -> usize {
        self.pos
    }
    pub fn data(&self) -> &'a [u8] {
        self.data
    }
    /// Reposition (`seekg` + `clear()`).
    pub fn seek(&mut self, pos: usize) {
        self.pos = pos.min(self.data.len());
        self.eof = false;
        self.fail = false;
    }

    /// `getLine(std::istream*)`: `std::getline` plus removal of one
    /// trailing `'\r'`.
    pub fn get_line(&mut self) -> Vec<u8> {
        if self.fail || self.eof {
            self.fail = true;
            return Vec::new();
        }
        if self.pos >= self.data.len() {
            self.eof = true;
            self.fail = true;
            return Vec::new();
        }
        let rest = &self.data[self.pos..];
        let mut line = match rest.iter().position(|&c| c == b'\n') {
            Some(nl) => {
                self.pos += nl + 1;
                rest[..nl].to_vec()
            }
            None => {
                self.pos = self.data.len();
                self.eof = true;
                rest.to_vec()
            }
        };
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        line
    }
}

pub(crate) fn lossy(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

/// `std::string::substr(pos, len)`: throws `std::out_of_range` when
/// `pos > size()`.
pub fn substr(s: &[u8], pos: usize, len: usize) -> Result<&[u8]> {
    if pos > s.len() {
        return Err(RDKitError::Index(format!(
            "basic_string::substr: __pos (which is {pos}) > this->size() (which is {})",
            s.len()
        )));
    }
    Ok(&s[pos..s.len().min(pos.saturating_add(len))])
}

/// `substr(pos)` to the end.
pub fn substr_from(s: &[u8], pos: usize) -> Result<&[u8]> {
    substr(s, pos, usize::MAX)
}

fn bad_cast() -> RDKitError {
    RDKitError::Value("bad lexical cast".into())
}

/// Whether an error came from a failed `boost::lexical_cast` (as opposed to
/// e.g. an out-of-range `substr`).
pub fn is_bad_cast(e: &RDKitError) -> bool {
    matches!(e, RDKitError::Value(m) if m == "bad lexical cast")
}

fn until_nul(input: &[u8]) -> &[u8] {
    match input.iter().position(|&c| c == 0) {
        Some(p) => &input[..p],
        None => input,
    }
}

fn skip_leading_spaces(input: &[u8], accept_spaces: bool) -> Option<&[u8]> {
    let mut t = input;
    if accept_spaces {
        while t.first() == Some(&b' ') {
            t = &t[1..];
        }
        if t.is_empty() {
            return None;
        }
    }
    Some(t)
}

/// `std::from_chars` for integers: an optional '-' (signed only) then
/// digits; on failure or overflow the result stays 0.
fn from_chars(t: &[u8], signed: bool) -> i64 {
    let (neg, digits) = if signed && t.first() == Some(&b'-') {
        (true, &t[1..])
    } else {
        (false, t)
    };
    let n = digits.iter().take_while(|c| c.is_ascii_digit()).count();
    if n == 0 {
        return 0;
    }
    let s = std::str::from_utf8(&digits[..n]).unwrap_or("0");
    let lim: i128 = if signed {
        if neg {
            i32::MAX as i128 + 1
        } else {
            i32::MAX as i128
        }
    } else {
        u32::MAX as i128
    };
    match s.parse::<i128>() {
        Ok(v) if v <= lim => {
            if neg {
                -(v as i64)
            } else {
                v as i64
            }
        }
        _ => 0,
    }
}

/// `FileParserUtils::toInt(input, acceptSpaces)`
pub fn to_int(input: &[u8], accept_spaces: bool) -> Result<i32> {
    let input = until_nul(input);
    for &c in input {
        if !(c.is_ascii_digit() || (accept_spaces && c == b' ') || c == b'+' || c == b'-') {
            return Err(bad_cast());
        }
    }
    let Some(t) = skip_leading_spaces(input, accept_spaces) else {
        return Ok(0);
    };
    Ok(from_chars(t, true) as i32)
}

/// `FileParserUtils::toUnsigned(input, acceptSpaces)`
pub fn to_unsigned(input: &[u8], accept_spaces: bool) -> Result<u32> {
    let input = until_nul(input);
    for &c in input {
        if !(c.is_ascii_digit() || (accept_spaces && c == b' ') || c == b'+') {
            return Err(bad_cast());
        }
    }
    let Some(t) = skip_leading_spaces(input, accept_spaces) else {
        return Ok(0);
    };
    Ok(from_chars(t, false) as u32)
}

/// `atof`: the longest decimal floating point prefix after whitespace.
pub fn atof(s: &[u8]) -> f64 {
    let s = until_nul(s);
    let mut i = 0;
    while i < s.len() && matches!(s[i], b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c) {
        i += 1;
    }
    let start = i;
    if i < s.len() && matches!(s[i], b'+' | b'-') {
        i += 1;
    }
    let mut nd = 0;
    while i < s.len() && s[i].is_ascii_digit() {
        i += 1;
        nd += 1;
    }
    if i < s.len() && s[i] == b'.' {
        i += 1;
        while i < s.len() && s[i].is_ascii_digit() {
            i += 1;
            nd += 1;
        }
    }
    if nd == 0 {
        return 0.0;
    }
    let mut end = i;
    if i < s.len() && matches!(s[i], b'e' | b'E') {
        let mut j = i + 1;
        if j < s.len() && matches!(s[j], b'+' | b'-') {
            j += 1;
        }
        let ds = j;
        while j < s.len() && s[j].is_ascii_digit() {
            j += 1;
        }
        if j > ds {
            end = j;
        }
    }
    let txt = std::str::from_utf8(&s[start..end]).unwrap_or("0");
    let txt = txt
        .strip_suffix('.')
        .map_or_else(|| txt.to_string(), |t| format!("{t}.0"));
    let txt = if txt.starts_with('.') || txt.starts_with("-.") || txt.starts_with("+.") {
        txt.replacen('.', "0.", 1)
    } else {
        txt
    };
    txt.parse::<f64>().unwrap_or(0.0)
}

/// `FileParserUtils::toDouble(text.substr(pos, len))`. Upstream validates
/// the field but then calls `atof` on the field's start pointer, which
/// reads on to the end of the line; `line[pos..]` reproduces that.
pub fn to_double_field(line: &[u8], pos: usize, len: usize) -> Result<f64> {
    let field = until_nul(substr(line, pos, len)?);
    for &c in field {
        if !(c.is_ascii_digit() || matches!(c, b' ' | b'+' | b'-' | b',' | b'.')) {
            return Err(bad_cast());
        }
    }
    Ok(atof(&line[pos..]))
}

/// `strip(input, " ")`
pub fn strip_spaces(input: &[u8]) -> &[u8] {
    let start = input.iter().position(|&c| c != b' ');
    match start {
        None => &input[input.len()..],
        Some(s) => {
            let e = input.iter().rposition(|&c| c != b' ').unwrap_or(s) + 1;
            &input[s..e]
        }
    }
}

/// `boost::lexical_cast<unsigned int>`: an optional sign then digits
/// (a leading '-' wraps around, like boost does).
pub fn lexical_cast_unsigned(t: &[u8]) -> Result<u32> {
    let (neg, digits) = match t.first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return Err(bad_cast());
    }
    let v: u64 = std::str::from_utf8(digits)
        .ok()
        .and_then(|s| s.parse().ok())
        .filter(|&v| v <= u32::MAX as u64)
        .ok_or_else(bad_cast)?;
    Ok(if neg {
        (v as u32).wrapping_neg()
    } else {
        v as u32
    })
}

/// `boost::lexical_cast<int>`
pub fn lexical_cast_int(t: &[u8]) -> Result<i32> {
    let (neg, digits) = match t.first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return Err(bad_cast());
    }
    let v: i64 = std::str::from_utf8(digits)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(bad_cast)?;
    let v = if neg { -v } else { v };
    i32::try_from(v).map_err(|_| bad_cast())
}

/// `stripSpacesAndCast<unsigned int>(input)`
pub fn strip_spaces_and_cast_unsigned(input: &[u8]) -> Result<u32> {
    lexical_cast_unsigned(strip_spaces(input))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ints() {
        assert_eq!(to_int(b"  3", true).unwrap(), 3);
        assert_eq!(to_int(b" -1", true).unwrap(), -1);
        assert_eq!(to_int(b"   ", true).unwrap(), 0);
        assert_eq!(to_int(b" +5", true).unwrap(), 0);
        assert_eq!(to_int(b"1 2", true).unwrap(), 1);
        assert!(to_int(b" a", true).is_err());
        assert!(to_unsigned(b" -1", true).is_err());
        assert_eq!(strip_spaces_and_cast_unsigned(b" 12 ").unwrap(), 12);
        assert!(strip_spaces_and_cast_unsigned(b"   ").is_err());
    }

    #[test]
    fn doubles() {
        assert_eq!(
            to_double_field(b"   -1.2500    2.0000", 0, 10).unwrap(),
            -1.25
        );
        assert_eq!(atof(b"12345.678912345.6789"), 12345.678912345);
        assert_eq!(atof(b"  1e3x"), 1000.0);
        assert_eq!(atof(b"  .5"), 0.5);
        assert!(to_double_field(b"   1.0000x", 0, 10).is_err());
    }

    #[test]
    fn lines() {
        let mut s = LineStream::new(b"a\r\nb\n\nc");
        assert_eq!(s.get_line(), b"a");
        assert_eq!(s.get_line(), b"b");
        assert_eq!(s.get_line(), b"");
        assert!(!s.eof());
        assert_eq!(s.get_line(), b"c");
        assert!(s.eof() && !s.fail());
        assert_eq!(s.get_line(), b"");
        assert!(s.fail());
    }
}
