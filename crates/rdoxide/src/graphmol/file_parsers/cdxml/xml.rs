//! A small XML reader producing the tree shape the CDXML parser expects from
//! `boost::property_tree::read_xml`: elements with attributes, the
//! concatenated text of their data nodes, and their child elements in
//! document order.

use crate::{RDKitError, Result};

/// An XML element.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Node {
    pub name: String,
    pub attrs: Vec<(String, String)>,
    /// concatenated character data (entities decoded, not trimmed)
    pub data: String,
    pub children: Vec<Node>,
}

impl Node {
    pub fn attr(&self, key: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }
}

struct Reader<'a> {
    src: &'a str,
    pos: usize,
}

fn error(src: &str, pos: usize, msg: &str) -> RDKitError {
    let line = src[..pos.min(src.len())].matches('\n').count() + 1;
    RDKitError::FileParse(format!("{msg} at line: {line}"))
}

impl<'a> Reader<'a> {
    fn rest(&self) -> &'a str {
        &self.src[self.pos..]
    }

    fn starts_with(&self, s: &str) -> bool {
        self.rest().starts_with(s)
    }

    fn fail<T>(&self, msg: &str) -> Result<T> {
        Err(error(self.src, self.pos, msg))
    }

    fn skip_ws(&mut self) {
        let t = self.rest().trim_start();
        self.pos = self.src.len() - t.len();
    }

    /// Skips to just after `end`.
    fn skip_past(&mut self, end: &str, what: &str) -> Result<()> {
        match self.rest().find(end) {
            Some(i) => {
                self.pos += i + end.len();
                Ok(())
            }
            None => self.fail(&format!("unterminated {what}")),
        }
    }

    fn name(&mut self) -> Result<String> {
        let t = self.rest();
        let end = t
            .find(|c: char| c.is_whitespace() || matches!(c, '=' | '>' | '/' | '<' | '"' | '\''))
            .unwrap_or(t.len());
        if end == 0 {
            return self.fail("expected element or attribute name");
        }
        self.pos += end;
        Ok(t[..end].to_string())
    }

    fn skip_prolog_item(&mut self) -> Result<bool> {
        if self.starts_with("<?") {
            self.skip_past("?>", "processing instruction")?;
        } else if self.starts_with("<!--") {
            self.skip_past("-->", "comment")?;
        } else if self.starts_with("<!DOCTYPE") {
            // an internal subset in brackets may contain '>'
            let t = self.rest();
            let mut depth = 0i32;
            let mut end = None;
            for (i, c) in t.char_indices() {
                match c {
                    '[' => depth += 1,
                    ']' => depth -= 1,
                    '>' if depth <= 0 => {
                        end = Some(i);
                        break;
                    }
                    _ => {}
                }
            }
            match end {
                Some(i) => self.pos += i + 1,
                None => return self.fail("unterminated DOCTYPE"),
            }
        } else {
            return Ok(false);
        }
        Ok(true)
    }

    fn element(&mut self) -> Result<Node> {
        // at '<'
        self.pos += 1;
        let name = self.name()?;
        let mut node = Node {
            name,
            ..Default::default()
        };
        loop {
            self.skip_ws();
            if self.starts_with("/>") {
                self.pos += 2;
                return Ok(node);
            }
            if self.starts_with(">") {
                self.pos += 1;
                break;
            }
            let key = self.name()?;
            self.skip_ws();
            if !self.starts_with("=") {
                return self.fail("expected '=' after attribute name");
            }
            self.pos += 1;
            self.skip_ws();
            let quote = match self.rest().chars().next() {
                Some(q @ ('"' | '\'')) => q,
                _ => return self.fail("expected quoted attribute value"),
            };
            self.pos += 1;
            let Some(len) = self.rest().find(quote) else {
                return self.fail("unterminated attribute value");
            };
            let raw = &self.rest()[..len];
            let value = decode(raw).map_err(|m| error(self.src, self.pos, &m))?;
            self.pos += len + 1;
            node.attrs.push((key, value));
        }
        // content
        loop {
            let rest = self.rest();
            let Some(lt) = rest.find('<') else {
                return self.fail(&format!("unterminated element <{}>", node.name));
            };
            if lt > 0 {
                let text = decode(&rest[..lt]).map_err(|m| error(self.src, self.pos, &m))?;
                node.data.push_str(&text);
                self.pos += lt;
            }
            if self.starts_with("</") {
                self.pos += 2;
                let close = self.name()?;
                self.skip_ws();
                if close != node.name || !self.starts_with(">") {
                    return self.fail(&format!(
                        "mismatched closing tag </{close}> for <{}>",
                        node.name
                    ));
                }
                self.pos += 1;
                return Ok(node);
            } else if self.starts_with("<!--") {
                self.skip_past("-->", "comment")?;
            } else if self.starts_with("<![CDATA[") {
                self.pos += 9;
                let Some(end) = self.rest().find("]]>") else {
                    return self.fail("unterminated CDATA section");
                };
                node.data.push_str(&self.rest()[..end]);
                self.pos += end + 3;
            } else if self.starts_with("<?") {
                self.skip_past("?>", "processing instruction")?;
            } else {
                let child = self.element()?;
                node.children.push(child);
            }
        }
    }
}

/// Decodes the predefined and numeric character references.
fn decode(raw: &str) -> std::result::Result<String, String> {
    if !raw.contains('&') {
        return Ok(raw.to_string());
    }
    let mut out = String::with_capacity(raw.len());
    let mut rest = raw;
    while let Some(i) = rest.find('&') {
        out.push_str(&rest[..i]);
        let tail = &rest[i + 1..];
        let Some(semi) = tail.find(';') else {
            return Err("invalid character reference".into());
        };
        let ent = &tail[..semi];
        let ch = match ent {
            "amp" => '&',
            "lt" => '<',
            "gt" => '>',
            "quot" => '"',
            "apos" => '\'',
            _ => {
                let code = ent
                    .strip_prefix("#x")
                    .or_else(|| ent.strip_prefix("#X"))
                    .map(|h| u32::from_str_radix(h, 16))
                    .or_else(|| ent.strip_prefix('#').map(|d| d.parse::<u32>()));
                match code {
                    Some(Ok(c)) => char::from_u32(c).ok_or("invalid character reference")?,
                    _ => return Err("invalid character reference".into()),
                }
            }
        };
        out.push(ch);
        rest = &tail[semi + 1..];
    }
    out.push_str(rest);
    Ok(out)
}

/// The top-level elements of an XML document, in order.
pub fn parse(src: &str) -> Result<Vec<Node>> {
    let mut r = Reader { src, pos: 0 };
    let mut roots = Vec::new();
    loop {
        r.skip_ws();
        if r.rest().is_empty() {
            break;
        }
        if r.skip_prolog_item()? {
            continue;
        }
        if r.starts_with("<") {
            roots.push(r.element()?);
        } else {
            return r.fail("expected '<'");
        }
    }
    if roots.is_empty() {
        return Err(error(src, src.len(), "no root element"));
    }
    Ok(roots)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_attributes_text_and_nesting() {
        let doc = "<?xml version=\"1.0\"?>\n<!DOCTYPE CDXML SYSTEM \"x.dtd\" >\n<CDXML a=\"1\" b='x&amp;y'><!-- c --><n id=\"2\"><t><s>R&#49;</s></t></n><b/></CDXML>";
        let roots = parse(doc).unwrap();
        assert_eq!(roots.len(), 1);
        let root = &roots[0];
        assert_eq!(root.attr("b"), Some("x&y"));
        assert_eq!(root.children.len(), 2);
        assert_eq!(root.children[0].children[0].children[0].data, "R1");
        assert_eq!(root.children[1].name, "b");
    }

    #[test]
    fn rejects_malformed_documents() {
        assert!(parse("<a><b></a>").is_err());
        assert!(parse("").is_err());
        assert!(parse("<a x=1/>").is_err());
    }
}
