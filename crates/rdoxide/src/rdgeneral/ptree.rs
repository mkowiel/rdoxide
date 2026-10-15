//! A small stand-in for `boost::property_tree` as RDKit uses it for JSON:
//! every value is a string, arrays are children with empty keys, and the
//! output of `write_json(..., pretty = false)` is reproduced byte for byte.
//!
//! Only what the JSON (de)serialization of RDKit objects needs is provided.

use std::fmt::Write;

use crate::{RDKitError, Result};

/// A property tree node.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct PTree {
    /// The node's own text (leaf value).
    data: String,
    /// Children in insertion order; an empty key marks an array element.
    children: Vec<(String, PTree)>,
}

impl PTree {
    pub fn new() -> Self {
        PTree::default()
    }

    /// A leaf holding `value`.
    pub fn leaf(value: impl ToString) -> Self {
        PTree {
            data: value.to_string(),
            children: Vec::new(),
        }
    }

    /// `pt.put(key, value)`: sets (or replaces) the leaf `key`.
    pub fn put(&mut self, key: &str, value: impl ToString) {
        self.add_child(key, PTree::leaf(value));
    }

    /// `pt.add_child(key, child)` with `put` semantics for an existing key.
    pub fn add_child(&mut self, key: &str, child: PTree) {
        match self.children.iter_mut().find(|(k, _)| k == key) {
            Some((_, c)) => *c = child,
            None => self.children.push((key.to_string(), child)),
        }
    }

    /// `pt.push_back(("", child))`: appends an array element.
    pub fn push(&mut self, child: PTree) {
        self.children.push((String::new(), child));
    }

    pub fn child(&self, key: &str) -> Option<&PTree> {
        self.children.iter().find(|(k, _)| k == key).map(|(_, c)| c)
    }

    /// The children (array elements or object members), in order.
    pub fn children(&self) -> impl Iterator<Item = &PTree> {
        self.children.iter().map(|(_, c)| c)
    }

    /// The members as (key, child) pairs.
    pub fn members(&self) -> impl Iterator<Item = (&str, &PTree)> {
        self.children.iter().map(|(k, c)| (k.as_str(), c))
    }

    pub fn value(&self) -> &str {
        &self.data
    }

    /// `pt.get_optional<std::string>(key)`.
    pub fn get_str(&self, key: &str) -> Option<&str> {
        self.child(key).map(PTree::value)
    }

    /// `pt.get<bool>(key, default)`: unparsable values give the default.
    pub fn get_bool(&self, key: &str, default: bool) -> bool {
        match self.get_str(key).map(str::trim) {
            Some("true" | "1") => true,
            Some("false" | "0") => false,
            _ => default,
        }
    }

    /// `pt.get<unsigned>(key, default)`: unparsable values give the default.
    pub fn get_u32(&self, key: &str, default: u32) -> u32 {
        self.get_str(key)
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(default)
    }

    /// `pt.get<int>(key, default)`: unparsable values give the default.
    pub fn get_i32(&self, key: &str, default: i32) -> i32 {
        self.get_str(key)
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(default)
    }

    /// `write_json(stream, pt, false)` (without the trailing newline).
    pub fn to_json(&self) -> String {
        let mut out = String::new();
        self.write(&mut out);
        out
    }

    fn write(&self, out: &mut String) {
        if self.children.is_empty() {
            write_string(&self.data, out);
            return;
        }
        let array = self.children.iter().all(|(k, _)| k.is_empty());
        out.push(if array { '[' } else { '{' });
        for (i, (key, child)) in self.children.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            if !array {
                write_string(key, out);
                out.push(':');
            }
            child.write(out);
        }
        out.push(if array { ']' } else { '}' });
    }

    /// `read_json`: objects and arrays become children; strings, numbers,
    /// booleans and `null` become leaf text.
    pub fn from_json(text: &str) -> Result<PTree> {
        let mut p = Parser {
            chars: text.chars().collect(),
            pos: 0,
        };
        p.skip_ws();
        let tree = p.value()?;
        p.skip_ws();
        if p.pos < p.chars.len() {
            return Err(p.error("garbage after data"));
        }
        Ok(tree)
    }
}

/// `create_escapes`: boost escapes `/` as well as the usual characters.
fn write_string(s: &str, out: &mut String) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '/' => out.push_str("\\/"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

struct Parser {
    chars: Vec<char>,
    pos: usize,
}

impl Parser {
    fn error(&self, what: &str) -> RDKitError {
        let line = 1 + self.chars[..self.pos.min(self.chars.len())]
            .iter()
            .filter(|&&c| c == '\n')
            .count();
        RDKitError::Runtime(format!("<unspecified file>({line}): {what}"))
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn skip_ws(&mut self) {
        while self.peek().is_some_and(|c| c.is_ascii_whitespace()) {
            self.pos += 1;
        }
    }

    fn eat(&mut self, c: char) -> bool {
        self.skip_ws();
        if self.peek() == Some(c) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn value(&mut self) -> Result<PTree> {
        self.skip_ws();
        match self.peek() {
            Some('{') => self.object(),
            Some('[') => self.array(),
            Some('"') => Ok(PTree::leaf(self.string()?)),
            Some(_) => self.literal(),
            None => Err(self.error("expected value")),
        }
    }

    fn object(&mut self) -> Result<PTree> {
        self.pos += 1;
        let mut tree = PTree::new();
        if self.eat('}') {
            return Ok(tree);
        }
        loop {
            self.skip_ws();
            if self.peek() != Some('"') {
                return Err(self.error("expected key string"));
            }
            let key = self.string()?;
            if !self.eat(':') {
                return Err(self.error("expected ':'"));
            }
            let child = self.value()?;
            tree.children.push((key, child));
            if self.eat(',') {
                continue;
            }
            if self.eat('}') {
                return Ok(tree);
            }
            return Err(self.error("expected ',' or '}'"));
        }
    }

    fn array(&mut self) -> Result<PTree> {
        self.pos += 1;
        let mut tree = PTree::new();
        if self.eat(']') {
            return Ok(tree);
        }
        loop {
            let child = self.value()?;
            tree.push(child);
            if self.eat(',') {
                continue;
            }
            if self.eat(']') {
                return Ok(tree);
            }
            return Err(self.error("expected ',' or ']'"));
        }
    }

    fn literal(&mut self) -> Result<PTree> {
        let start = self.pos;
        while self
            .peek()
            .is_some_and(|c| !matches!(c, ',' | '}' | ']') && !c.is_ascii_whitespace())
        {
            self.pos += 1;
        }
        let text: String = self.chars[start..self.pos].iter().collect();
        let numeric = text
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_digit() || c == '-');
        if numeric || matches!(text.as_str(), "true" | "false" | "null") {
            Ok(PTree::leaf(text))
        } else {
            self.pos = start;
            Err(self.error("expected value"))
        }
    }

    fn string(&mut self) -> Result<String> {
        self.pos += 1;
        let mut out = String::new();
        loop {
            let c = self
                .peek()
                .ok_or_else(|| self.error("unterminated string"))?;
            self.pos += 1;
            match c {
                '"' => return Ok(out),
                '\\' => {
                    let e = self.peek().ok_or_else(|| self.error("bad escape"))?;
                    self.pos += 1;
                    match e {
                        'b' => out.push('\u{8}'),
                        'f' => out.push('\u{c}'),
                        'n' => out.push('\n'),
                        'r' => out.push('\r'),
                        't' => out.push('\t'),
                        'u' => out.push(self.unicode_escape()?),
                        '"' | '\\' | '/' => out.push(e),
                        _ => return Err(self.error("invalid escape sequence")),
                    }
                }
                c => out.push(c),
            }
        }
    }

    fn hex4(&mut self) -> Result<u32> {
        let digits: String = self.chars[self.pos..].iter().take(4).collect();
        if digits.len() != 4 {
            return Err(self.error("invalid code sequence"));
        }
        self.pos += 4;
        u32::from_str_radix(&digits, 16).map_err(|_| self.error("invalid code sequence"))
    }

    fn unicode_escape(&mut self) -> Result<char> {
        let hi = self.hex4()?;
        let code = if (0xD800..0xDC00).contains(&hi) {
            if self.chars[self.pos..].starts_with(&['\\', 'u']) {
                self.pos += 2;
                let lo = self.hex4()?;
                0x10000 + ((hi - 0xD800) << 10) + (lo.wrapping_sub(0xDC00) & 0x3FF)
            } else {
                return Err(self.error("invalid code sequence"));
            }
        } else {
            hi
        };
        char::from_u32(code).ok_or_else(|| self.error("invalid code sequence"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_like_boost() {
        let mut pt = PTree::new();
        pt.put("name", "a/b\"c");
        pt.put("flag", true);
        let mut arr = PTree::new();
        arr.push(PTree::leaf(1));
        arr.push(PTree::leaf(2));
        pt.add_child("list", arr);
        pt.add_child("empty", PTree::new());
        assert_eq!(
            pt.to_json(),
            r#"{"name":"a\/b\"c","flag":"true","list":["1","2"],"empty":""}"#
        );
    }

    #[test]
    fn reads_what_it_writes_and_plain_json() {
        let pt = PTree::from_json(r#" {"a": 2, "b": [true, "x"], "c": {"d": null}} "#).unwrap();
        assert_eq!(pt.get_u32("a", 0), 2);
        assert_eq!(pt.child("b").unwrap().children().count(), 2);
        assert!(pt.child("b").unwrap().children().next().unwrap().value() == "true");
        assert_eq!(pt.child("c").unwrap().get_str("d"), Some("null"));
        let again = PTree::from_json(&pt.to_json()).unwrap();
        assert_eq!(again.to_json(), pt.to_json());
    }

    #[test]
    fn typed_getters_fall_back_to_defaults() {
        let pt = PTree::from_json(r#"{"n":"x","b":"1","u":" 7 "}"#).unwrap();
        assert_eq!(pt.get_u32("n", 9), 9);
        assert!(pt.get_bool("b", false));
        assert_eq!(pt.get_u32("u", 0), 7);
        assert_eq!(pt.get_u32("missing", 3), 3);
    }

    #[test]
    fn rejects_malformed_json() {
        assert!(PTree::from_json("").is_err());
        assert!(PTree::from_json("{").is_err());
        assert!(PTree::from_json("{} x").is_err());
        assert!(PTree::from_json(r#"{"a" 1}"#).is_err());
    }
}
