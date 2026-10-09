//! A minimal JSON document model standing in for `boost::json`, as used by
//! the `MolInterchange` reader/writer.
//!
//! What matters for compatibility is reproduced: objects keep insertion
//! order, integers and doubles are distinct types (the typed accessors are as
//! strict as boost's), and [`Json::serialize`] writes exactly what
//! `boost::json::serialize` does (compact, doubles as shortest round-trip
//! digits in `d.dddE±x` form, non-finite doubles as `null` / `±1e99999`).
//!
//! The `number` module translates the number reader of Boost.JSON 1.85
//! (`boost/json/basic_parser_impl.hpp`, Copyright (c) 2019 Vinnie Falco,
//! (c) 2020 Krystian Stasiowski), distributed under the Boost Software
//! License 1.0 (LICENSES/Boost-BSL-1.0.txt). The rest of the file is original.

use std::fmt::Write;

use crate::{RDKitError, Result};

/// Maximum nesting depth accepted by the parser (`boost::json` default).
const MAX_DEPTH: usize = 32;

/// A JSON object: members in insertion order.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct JsonObject(Vec<(String, Json)>);

impl JsonObject {
    pub fn new() -> Self {
        JsonObject::default()
    }

    /// `obj[key] = value`: replaces an existing member in place.
    pub fn insert(&mut self, key: &str, value: impl Into<Json>) {
        let value = value.into();
        match self.0.iter_mut().find(|(k, _)| k == key) {
            Some((_, v)) => *v = value,
            None => self.0.push((key.to_string(), value)),
        }
    }

    pub fn get(&self, key: &str) -> Option<&Json> {
        self.0.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    pub fn contains(&self, key: &str) -> bool {
        self.get(key).is_some()
    }

    /// `obj.at(key)`
    pub fn at(&self, key: &str) -> Result<&Json> {
        self.get(key).ok_or_else(out_of_range)
    }

    pub fn iter(&self) -> impl Iterator<Item = (&str, &Json)> {
        self.0.iter().map(|(k, v)| (k.as_str(), v))
    }
}

/// A JSON value.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum Json {
    #[default]
    Null,
    Bool(bool),
    Int(i64),
    UInt(u64),
    Double(f64),
    String(String),
    Array(Vec<Json>),
    Object(JsonObject),
}

fn out_of_range() -> RDKitError {
    RDKitError::Runtime("out of range".into())
}

fn wrong_type(what: &str) -> RDKitError {
    RDKitError::Runtime(format!("value is not {what}"))
}

macro_rules! json_from {
    ($($t:ty => $variant:ident as $cast:ty),* $(,)?) => {$(
        impl From<$t> for Json {
            fn from(v: $t) -> Self {
                Json::$variant(v as $cast)
            }
        }
    )*};
}
json_from!(i32 => Int as i64, i64 => Int as i64, u32 => UInt as u64, usize => UInt as u64, f64 => Double as f64);

impl From<bool> for Json {
    fn from(v: bool) -> Self {
        Json::Bool(v)
    }
}

impl From<&str> for Json {
    fn from(v: &str) -> Self {
        Json::String(v.to_string())
    }
}

impl From<String> for Json {
    fn from(v: String) -> Self {
        Json::String(v)
    }
}

impl From<Vec<Json>> for Json {
    fn from(v: Vec<Json>) -> Self {
        Json::Array(v)
    }
}

impl From<JsonObject> for Json {
    fn from(v: JsonObject) -> Self {
        Json::Object(v)
    }
}

impl Json {
    /// `value::as_object()`
    pub fn as_object(&self) -> Result<&JsonObject> {
        match self {
            Json::Object(o) => Ok(o),
            _ => Err(wrong_type("an object")),
        }
    }

    /// `value::as_array()`
    pub fn as_array(&self) -> Result<&[Json]> {
        match self {
            Json::Array(a) => Ok(a),
            _ => Err(wrong_type("an array")),
        }
    }

    /// `value::as_string()`
    pub fn as_str(&self) -> Result<&str> {
        match self {
            Json::String(s) => Ok(s),
            _ => Err(wrong_type("a string")),
        }
    }

    /// `value::as_int64()`; unsigned values that fit are accepted.
    pub fn as_i64(&self) -> Result<i64> {
        match self {
            Json::Int(i) => Ok(*i),
            Json::UInt(u) => i64::try_from(*u).map_err(|_| wrong_type("an int64")),
            _ => Err(wrong_type("an int64")),
        }
    }

    /// `value::as_double()`: only genuine doubles qualify, as in boost.
    pub fn as_f64(&self) -> Result<f64> {
        match self {
            Json::Double(d) => Ok(*d),
            _ => Err(wrong_type("a double")),
        }
    }

    /// `value::as_bool()`
    pub fn as_bool(&self) -> Result<bool> {
        match self {
            Json::Bool(b) => Ok(*b),
            _ => Err(wrong_type("a bool")),
        }
    }

    pub fn is_object(&self) -> bool {
        matches!(self, Json::Object(_))
    }

    pub fn is_array(&self) -> bool {
        matches!(self, Json::Array(_))
    }

    /// `value::at(key)` on an object.
    pub fn at(&self, key: &str) -> Result<&Json> {
        self.as_object()?.at(key)
    }

    /// `value::at(idx)` on an array.
    pub fn at_idx(&self, idx: usize) -> Result<&Json> {
        self.as_array()?.get(idx).ok_or_else(out_of_range)
    }

    /// `boost::json::parse`
    pub fn parse(text: &str) -> Result<Json> {
        let mut p = Parser {
            bytes: text.as_bytes(),
            pos: 0,
        };
        p.skip_ws();
        let v = p.value(1)?;
        p.skip_ws();
        if p.pos != p.bytes.len() {
            return Err(parse_error("extra data"));
        }
        Ok(v)
    }

    /// `boost::json::serialize`
    pub fn serialize(&self) -> String {
        let mut out = String::new();
        self.write(&mut out);
        out
    }

    fn write(&self, out: &mut String) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Json::Int(i) => {
                let _ = write!(out, "{i}");
            }
            Json::UInt(u) => {
                let _ = write!(out, "{u}");
            }
            Json::Double(d) => write_double(out, *d),
            Json::String(s) => write_string(out, s),
            Json::Array(a) => {
                out.push('[');
                for (i, v) in a.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    v.write(out);
                }
                out.push(']');
            }
            Json::Object(o) => {
                out.push('{');
                for (i, (k, v)) in o.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    write_string(out, k);
                    out.push(':');
                    v.write(out);
                }
                out.push('}');
            }
        }
    }
}

fn write_double(out: &mut String, d: f64) {
    if d.is_nan() {
        out.push_str("null");
    } else if d.is_infinite() {
        out.push_str(if d < 0.0 { "-1e99999" } else { "1e99999" });
    } else {
        // shortest round-trip digits, always in scientific notation
        let _ = write!(out, "{d:E}");
    }
}

fn write_string(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn parse_error(what: &str) -> RDKitError {
    RDKitError::Runtime(what.to_string())
}

struct Parser<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl Parser<'_> {
    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn skip_ws(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.pos += 1;
        }
    }

    /// Consumes `b`, failing with "incomplete JSON" at the end of input and
    /// "syntax error" on any other byte.
    fn expect(&mut self, b: u8) -> Result<()> {
        match self.peek() {
            Some(c) if c == b => {
                self.pos += 1;
                Ok(())
            }
            Some(_) => Err(parse_error("syntax error")),
            None => Err(parse_error("incomplete JSON")),
        }
    }

    fn literal(&mut self, word: &str, value: Json) -> Result<Json> {
        for b in word.bytes() {
            self.expect(b)?;
        }
        Ok(value)
    }

    fn value(&mut self, depth: usize) -> Result<Json> {
        match self.peek() {
            None => Err(parse_error("incomplete JSON")),
            Some(b'n') => self.literal("null", Json::Null),
            Some(b't') => self.literal("true", Json::Bool(true)),
            Some(b'f') => self.literal("false", Json::Bool(false)),
            Some(b'"') => Ok(Json::String(self.string()?)),
            Some(b'[') => self.array(depth),
            Some(b'{') => self.object(depth),
            Some(b'-' | b'0'..=b'9') => self.number(),
            Some(_) => Err(parse_error("syntax error")),
        }
    }

    fn array(&mut self, depth: usize) -> Result<Json> {
        if depth > MAX_DEPTH {
            return Err(parse_error("too deep"));
        }
        self.expect(b'[')?;
        let mut items = Vec::new();
        self.skip_ws();
        if self.peek() == Some(b']') {
            self.pos += 1;
            return Ok(Json::Array(items));
        }
        loop {
            self.skip_ws();
            items.push(self.value(depth + 1)?);
            self.skip_ws();
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b']') => {
                    self.pos += 1;
                    return Ok(Json::Array(items));
                }
                Some(_) => return Err(parse_error("syntax error")),
                None => return Err(parse_error("incomplete JSON")),
            }
        }
    }

    fn object(&mut self, depth: usize) -> Result<Json> {
        if depth > MAX_DEPTH {
            return Err(parse_error("too deep"));
        }
        self.expect(b'{')?;
        let mut obj = JsonObject::new();
        self.skip_ws();
        if self.peek() == Some(b'}') {
            self.pos += 1;
            return Ok(Json::Object(obj));
        }
        loop {
            self.skip_ws();
            if self.peek() != Some(b'"') {
                return Err(parse_error(if self.peek().is_none() {
                    "incomplete JSON"
                } else {
                    "syntax error"
                }));
            }
            let key = self.string()?;
            self.skip_ws();
            self.expect(b':')?;
            self.skip_ws();
            let v = self.value(depth + 1)?;
            // duplicate keys: the last one wins
            obj.insert(&key, v);
            self.skip_ws();
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b'}') => {
                    self.pos += 1;
                    return Ok(Json::Object(obj));
                }
                Some(_) => return Err(parse_error("syntax error")),
                None => return Err(parse_error("incomplete JSON")),
            }
        }
    }

    fn hex4(&mut self) -> Result<u32> {
        let mut v = 0u32;
        for _ in 0..4 {
            let c = self.peek().ok_or_else(|| parse_error("incomplete JSON"))?;
            let d = (c as char)
                .to_digit(16)
                .ok_or_else(|| parse_error("syntax error"))?;
            v = v * 16 + d;
            self.pos += 1;
        }
        Ok(v)
    }

    fn string(&mut self) -> Result<String> {
        self.expect(b'"')?;
        let mut out = String::new();
        loop {
            // copy the run of ordinary bytes in one go (input is valid UTF-8
            // and the delimiters are ASCII)
            let start = self.pos;
            while let Some(c) = self.peek() {
                if c == b'"' || c == b'\\' || c < 0x20 {
                    break;
                }
                self.pos += 1;
            }
            out.push_str(
                std::str::from_utf8(&self.bytes[start..self.pos])
                    .map_err(|_| parse_error("syntax error"))?,
            );
            match self.peek() {
                None => return Err(parse_error("incomplete JSON")),
                Some(b'"') => {
                    self.pos += 1;
                    return Ok(out);
                }
                Some(b'\\') => {
                    self.pos += 1;
                    let esc = self.peek().ok_or_else(|| parse_error("incomplete JSON"))?;
                    self.pos += 1;
                    match esc {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'b' => out.push('\u{8}'),
                        b'f' => out.push('\u{c}'),
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'u' => out.push(self.unicode_escape()?),
                        _ => return Err(parse_error("syntax error")),
                    }
                }
                Some(_) => return Err(parse_error("syntax error")),
            }
        }
    }

    /// The part of a `\uXXXX` escape after the `u`, joining surrogate pairs.
    fn unicode_escape(&mut self) -> Result<char> {
        let first = self.hex4()?;
        let code = if (0xD800..0xDC00).contains(&first) {
            self.expect(b'\\')?;
            self.expect(b'u')?;
            let second = self.hex4()?;
            if !(0xDC00..0xE000).contains(&second) {
                return Err(parse_error("illegal trailing surrogate"));
            }
            0x10000 + ((first - 0xD800) << 10) + (second - 0xDC00)
        } else if (0xDC00..0xE000).contains(&first) {
            return Err(parse_error("illegal leading surrogate"));
        } else {
            first
        };
        char::from_u32(code).ok_or_else(|| parse_error("syntax error"))
    }

    fn number(&mut self) -> Result<Json> {
        let (value, used) = number::parse(&self.bytes[self.pos..])?;
        self.pos += used;
        Ok(value)
    }
}

/// Number parsing as `boost::json::basic_parser` does it with its default
/// (imprecise) number mode: the mantissa and decimal exponent are collected
/// first and combined with one multiplication or division, so the last bit of
/// a double can differ from a correctly rounded conversion. The results of
/// this reader are therefore bit-identical to what RDKit reads from a file.
///
/// The exact digit bookkeeping (which digits are kept, which only shift the
/// exponent) follows boost's two code paths: the fast one for numbers with
/// at least 33 bytes of input after them, the general one otherwise.
mod number {
    use std::sync::OnceLock;

    use super::{Json, parse_error};
    use crate::Result;

    const INT_MAX: i64 = i32::MAX as i64;
    /// 2^53 - 1: the mantissa stops growing past this in the general path
    const MANT_LIMIT: u64 = 9_007_199_254_740_991;

    /// Powers of ten from `1e-308` to `1e308`, correctly rounded.
    fn pow10(exp: i64) -> f64 {
        static TABLE: OnceLock<Vec<f64>> = OnceLock::new();
        if exp > 308 {
            return f64::INFINITY;
        }
        if exp < -308 {
            return 0.0;
        }
        let table = TABLE.get_or_init(|| {
            (-308..=308)
                .map(|e| format!("1e{e}").parse().unwrap_or(0.0))
                .collect()
        });
        table[(exp + 308) as usize]
    }

    /// `detail::dec_to_float`
    fn dec_to_float(mant: u64, mut exp: i64, neg: bool) -> f64 {
        let mut x = mant as f64;
        if neg {
            x = -x;
        }
        if exp < -305 {
            x *= 1e-305;
            exp += 305;
        }
        if (-22..0).contains(&exp) {
            return x / pow10(-exp);
        }
        x * pow10(exp)
    }

    struct Reader<'a> {
        bytes: &'a [u8],
        pos: usize,
        mant: u64,
        /// decimal exponent adjustment from the digits
        bias: i64,
        exp: i64,
        exp_negative: bool,
        neg: bool,
    }

    /// What follows the digits that have been read.
    enum Next {
        /// `finish_int`
        Int,
        /// `finish_dub`
        Double,
        /// `[.eE]` after the integer part
        Point,
        Exponent,
    }

    impl Reader<'_> {
        fn peek(&self) -> Option<u8> {
            self.bytes.get(self.pos).copied()
        }

        fn digit(&self) -> Option<u64> {
            self.peek()
                .filter(u8::is_ascii_digit)
                .map(|c| u64::from(c - b'0'))
        }

        fn digit_required(&self) -> Result<()> {
            match self.peek() {
                Some(c) if c.is_ascii_digit() => Ok(()),
                Some(_) => Err(parse_error("syntax error")),
                None => Err(parse_error("incomplete JSON")),
            }
        }

        /// Up to 16 leading digits (`count_digits`).
        fn count_digits(&self) -> usize {
            self.bytes[self.pos..]
                .iter()
                .take(16)
                .take_while(|c| c.is_ascii_digit())
                .count()
        }

        fn take_digits(&mut self, n: usize) {
            for _ in 0..n {
                self.mant = 10 * self.mant + self.digit().unwrap_or(0);
                self.pos += 1;
            }
        }

        /// `do_num2`: significant digits left of the decimal point.
        fn integer_digits(&mut self) -> Result<Next> {
            loop {
                let Some(d) = self.digit() else {
                    return Ok(if self.peek().is_none() {
                        Next::Int
                    } else {
                        Next::Point
                    });
                };
                self.pos += 1;
                let overflow = if self.neg {
                    // 9223372036854775808 is INT64_MIN
                    self.mant > 922_337_203_685_477_580
                        || (self.mant == 922_337_203_685_477_580 && d > 8)
                } else {
                    self.mant > 1_844_674_407_370_955_161
                        || (self.mant == 1_844_674_407_370_955_161 && d > 5)
                };
                if overflow {
                    break;
                }
                self.mant = 10 * self.mant + d;
            }
            self.bias += 1;
            self.insignificant_integer_digits()
        }

        /// `do_num3`
        fn insignificant_integer_digits(&mut self) -> Result<Next> {
            loop {
                match self.peek() {
                    Some(c) if c.is_ascii_digit() => {
                        self.pos += 1;
                        self.bias += 1;
                    }
                    Some(b'.') => {
                        self.pos += 1;
                        return self.dropped_fraction_digits();
                    }
                    Some(b'e' | b'E') => {
                        self.pos += 1;
                        return Ok(Next::Exponent);
                    }
                    _ => return Ok(Next::Double),
                }
            }
        }

        /// `do_num4`/`do_num5`: digits right of the point that are dropped.
        fn dropped_fraction_digits(&mut self) -> Result<Next> {
            self.digit_required()?;
            self.pos += 1;
            Ok(self.skip_fraction_digits())
        }

        fn skip_fraction_digits(&mut self) -> Next {
            while self.digit().is_some() {
                self.pos += 1;
            }
            self.exponent_or_end()
        }

        fn exponent_or_end(&mut self) -> Next {
            if matches!(self.peek(), Some(b'e' | b'E')) {
                self.pos += 1;
                Next::Exponent
            } else {
                Next::Double
            }
        }

        /// `do_num8`: significant digits right of the decimal point.
        fn fraction_digits(&mut self) -> Next {
            while let Some(d) = self.digit() {
                self.pos += 1;
                if self.mant > MANT_LIMIT {
                    return self.skip_fraction_digits();
                }
                self.bias -= 1;
                self.mant = 10 * self.mant + d;
            }
            self.exponent_or_end()
        }

        /// `do_exp1` to `do_exp3`
        fn exponent(&mut self) -> Result<()> {
            match self.peek() {
                Some(b'+') => self.pos += 1,
                Some(b'-') => {
                    self.pos += 1;
                    self.exp_negative = true;
                }
                _ => {}
            }
            self.digit_required()?;
            while let Some(d) = self.digit() {
                self.pos += 1;
                self.exp = (10 * self.exp + d as i64).min(INT_MAX);
            }
            Ok(())
        }

        /// The digits after the first one, in the general path.
        fn general(&mut self) -> Result<Next> {
            match self.peek() {
                Some(b'0') => {
                    self.pos += 1;
                    self.mant = 0;
                    Ok(Next::Point)
                }
                Some(c @ b'1'..=b'9') => {
                    self.pos += 1;
                    self.mant = u64::from(c - b'0');
                    self.integer_digits()
                }
                Some(_) => Err(parse_error("syntax error")),
                None => Err(parse_error("incomplete JSON")),
            }
        }

        /// `do_num6`: after the integer part.
        fn after_integer(&mut self) -> Result<Next> {
            match self.peek() {
                Some(b'.') => {
                    self.pos += 1;
                    self.digit_required()?;
                    Ok(self.fraction_digits())
                }
                Some(b'e' | b'E') => {
                    self.pos += 1;
                    Ok(Next::Exponent)
                }
                _ => Ok(Next::Int),
            }
        }

        /// The fast path for numbers with plenty of input behind them.
        fn fast(&mut self) -> Result<Next> {
            let n1 = if self.peek() == Some(b'0') {
                self.pos += 1;
                0
            } else {
                let n1 = self.count_digits();
                if n1 == 0 {
                    return Err(parse_error("syntax error"));
                }
                self.take_digits(n1);
                if n1 == 16 {
                    return self.integer_digits();
                }
                n1
            };
            match self.peek() {
                Some(b'.') => {}
                Some(b'e' | b'E') => {
                    self.pos += 1;
                    return Ok(Next::Exponent);
                }
                _ => return Ok(Next::Int),
            }
            self.pos += 1;
            let n2 = self.count_digits();
            if n2 == 0 {
                return Err(parse_error("syntax error"));
            }
            if n1 + n2 >= 19 {
                return Ok(self.fraction_digits());
            }
            self.take_digits(n2);
            self.bias -= n2 as i64;
            match self.peek() {
                Some(b'e' | b'E') => {
                    self.pos += 1;
                    Ok(Next::Exponent)
                }
                Some(c) if c.is_ascii_digit() => Ok(self.fraction_digits()),
                _ => Ok(Next::Double),
            }
        }
    }

    /// Parses the number at the start of `bytes`; returns it with its length.
    pub(super) fn parse(bytes: &[u8]) -> Result<(Json, usize)> {
        let neg = bytes.first() == Some(&b'-');
        let start = usize::from(neg);
        let mut r = Reader {
            bytes,
            pos: start,
            mant: 0,
            bias: 0,
            exp: 0,
            exp_negative: false,
            neg,
        };
        let mut next = if bytes.len() - start >= 16 + 1 + 16 {
            r.fast()?
        } else {
            r.general()?
        };
        if let Next::Point = next {
            next = r.after_integer()?;
        }
        if let Next::Exponent = next {
            r.exponent()?;
            next = Next::Double;
        }
        let value = match next {
            Next::Int => {
                if neg {
                    Json::Int((!r.mant).wrapping_add(1) as i64)
                } else if let Ok(i) = i64::try_from(r.mant) {
                    Json::Int(i)
                } else {
                    Json::UInt(r.mant)
                }
            }
            _ => {
                let exp = r.bias + if r.exp_negative { -r.exp } else { r.exp };
                Json::Double(dec_to_float(r.mant, exp, neg))
            }
        };
        Ok((value, r.pos))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serializes_like_boost() {
        let mut o = JsonObject::new();
        o.insert("a", 1);
        o.insert("d", 1.0);
        o.insert("e", 0.1);
        o.insert("s", "x\"\n\u{1}é");
        o.insert("a", 2);
        o.insert("n", f64::NAN);
        o.insert("i", f64::INFINITY);
        o.insert("z", -0.0);
        assert_eq!(
            Json::Object(o).serialize(),
            r#"{"a":2,"d":1E0,"e":1E-1,"s":"x\"\n\u0001é","n":null,"i":1e99999,"z":-0E0}"#
        );
    }

    #[test]
    fn parses_and_distinguishes_ints_and_doubles() {
        let v = Json::parse(r#" {"a":[1,-2,3.5,1e2,18446744073709551615],"b":"é😀"} "#).unwrap();
        let a = v.at("a").unwrap().as_array().unwrap();
        assert_eq!(a[0].as_i64().unwrap(), 1);
        assert!(a[0].as_f64().is_err());
        assert_eq!(a[2].as_f64().unwrap(), 3.5);
        assert_eq!(a[3].as_f64().unwrap(), 100.0);
        assert_eq!(a[4], Json::UInt(u64::MAX));
        assert_eq!(v.at("b").unwrap().as_str().unwrap(), "é😀");
    }

    #[test]
    fn doubles_are_read_like_boost_json() {
        // boost does not round correctly: with plenty of input behind it
        // this reads as 0.117, not as the nearest double 0.11699999999999999
        let padding = " ".repeat(40);
        let v = Json::parse(&format!("[1.1699999999999999E-1{padding}]")).unwrap();
        assert_eq!(v.at_idx(0).unwrap().as_f64().unwrap(), 0.117);
        // at the end of the input the general path keeps fewer digits
        let v = Json::parse("[1.1699999999999999E-1]").unwrap();
        assert_ne!(v.at_idx(0).unwrap().as_f64().unwrap(), 0.0);
        assert_eq!(Json::parse("-0").unwrap(), Json::Int(0));
        assert_eq!(
            Json::parse("9223372036854775808").unwrap(),
            Json::UInt(1 << 63)
        );
        assert_eq!(
            Json::parse("-9223372036854775808").unwrap(),
            Json::Int(i64::MIN)
        );
        assert_eq!(
            Json::parse("1e400").unwrap().as_f64().unwrap(),
            f64::INFINITY
        );
    }

    #[test]
    fn rejects_malformed_input() {
        for bad in [
            "", "{", "[1,]", "{\"a\"}", "01", "1.", "nul", "[1] x", "\"a",
        ] {
            assert!(Json::parse(bad).is_err(), "{bad}");
        }
        let deep = "[".repeat(40) + &"]".repeat(40);
        assert!(Json::parse(&deep).is_err());
        let ok = "[".repeat(32) + &"]".repeat(32);
        assert!(Json::parse(&ok).is_ok());
    }
}
