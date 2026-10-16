//! A small XML reader producing the tree `boost::property_tree::read_xml`
//! would: attributes are the children of an `<xmlattr>` node, the text of an
//! element is its data and comments are `<xmlcomment>` children.

use crate::{RDKitError, Result};

/// A node of the tree: its own text and its named children in document
/// order (names may repeat).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct XmlTree {
    pub data: String,
    pub children: Vec<(String, XmlTree)>,
}

impl XmlTree {
    /// `get_child("a.b.c")`: the first child of each step.
    pub fn get_child(&self, path: &str) -> Option<&XmlTree> {
        path.split('.').try_fold(self, |node, key| node.child(key))
    }

    /// The first child called `key`.
    pub fn child(&self, key: &str) -> Option<&XmlTree> {
        self.children.iter().find(|(k, _)| k == key).map(|(_, c)| c)
    }

    /// `get<std::string>(path, default)`: the data of the node at `path`.
    pub fn get_str(&self, path: &str, default: &str) -> String {
        self.get_child(path)
            .map_or_else(|| default.to_string(), |n| n.data.clone())
    }

    /// The value of an attribute (`<xmlattr>.name`).
    pub fn attr(&self, name: &str, default: &str) -> String {
        self.child("<xmlattr>")
            .and_then(|a| a.child(name))
            .map_or_else(|| default.to_string(), |n| n.data.clone())
    }

    /// The children called `key`, in order.
    pub fn children_named<'a>(&'a self, key: &'a str) -> impl Iterator<Item = &'a XmlTree> {
        self.children
            .iter()
            .filter(move |(k, _)| k == key)
            .map(|(_, c)| c)
    }
}

fn parse_error(msg: &str, pos: usize) -> RDKitError {
    RDKitError::Runtime(format!("<unspecified file>: {msg} at byte {pos}"))
}

struct Parser<'a> {
    s: &'a [u8],
    pos: usize,
}

fn is_name_char(c: u8) -> bool {
    c.is_ascii_alphanumeric() || matches!(c, b'_' | b':' | b'-' | b'.') || c >= 0x80
}

/// Replaces the predefined and numeric character entities.
fn decode_entities(text: &str, pos: usize) -> Result<String> {
    if !text.contains('&') {
        return Ok(text.to_string());
    }
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(i) = rest.find('&') {
        out.push_str(&rest[..i]);
        let tail = &rest[i + 1..];
        let end = tail
            .find(';')
            .ok_or_else(|| parse_error("unterminated entity", pos))?;
        let entity = &tail[..end];
        match entity {
            "amp" => out.push('&'),
            "lt" => out.push('<'),
            "gt" => out.push('>'),
            "quot" => out.push('"'),
            "apos" => out.push('\''),
            _ => {
                let code = entity
                    .strip_prefix("#x")
                    .and_then(|h| u32::from_str_radix(h, 16).ok())
                    .or_else(|| entity.strip_prefix('#').and_then(|d| d.parse().ok()))
                    .and_then(char::from_u32)
                    .ok_or_else(|| parse_error("unknown entity", pos))?;
                out.push(code);
            }
        }
        rest = &tail[end + 1..];
    }
    out.push_str(rest);
    Ok(out)
}

impl Parser<'_> {
    fn starts_with(&self, p: &str) -> bool {
        self.s[self.pos..].starts_with(p.as_bytes())
    }

    fn skip_ws(&mut self) {
        while self.pos < self.s.len() && self.s[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    fn skip_past(&mut self, end: &str) -> Result<()> {
        let hay = &self.s[self.pos..];
        let n = end.len();
        let i = hay
            .windows(n)
            .position(|w| w == end.as_bytes())
            .ok_or_else(|| parse_error(&format!("expected {end}"), self.pos))?;
        self.pos += i + n;
        Ok(())
    }

    fn name(&mut self) -> Result<String> {
        let start = self.pos;
        while self.pos < self.s.len() && is_name_char(self.s[self.pos]) {
            self.pos += 1;
        }
        if start == self.pos {
            return Err(parse_error("expected a name", self.pos));
        }
        Ok(String::from_utf8_lossy(&self.s[start..self.pos]).into_owned())
    }

    /// Parses an element whose `<` has been consumed.
    fn element(&mut self) -> Result<(String, XmlTree)> {
        let name = self.name()?;
        let mut node = XmlTree::default();
        let mut attrs = XmlTree::default();
        loop {
            self.skip_ws();
            match self.s.get(self.pos) {
                Some(b'/') => {
                    self.pos += 1;
                    if self.s.get(self.pos) != Some(&b'>') {
                        return Err(parse_error("expected >", self.pos));
                    }
                    self.pos += 1;
                    self.finish_attrs(&mut node, attrs);
                    return Ok((name, node));
                }
                Some(b'>') => {
                    self.pos += 1;
                    break;
                }
                Some(_) => {
                    let key = self.name()?;
                    self.skip_ws();
                    if self.s.get(self.pos) != Some(&b'=') {
                        return Err(parse_error("expected =", self.pos));
                    }
                    self.pos += 1;
                    self.skip_ws();
                    let quote = match self.s.get(self.pos) {
                        Some(&q) if q == b'"' || q == b'\'' => q,
                        _ => return Err(parse_error("expected a quoted value", self.pos)),
                    };
                    self.pos += 1;
                    let start = self.pos;
                    while self.pos < self.s.len() && self.s[self.pos] != quote {
                        self.pos += 1;
                    }
                    if self.pos >= self.s.len() {
                        return Err(parse_error("unterminated attribute value", start));
                    }
                    let raw = String::from_utf8_lossy(&self.s[start..self.pos]).into_owned();
                    self.pos += 1;
                    let value = decode_entities(&raw, start)?;
                    // a repeated attribute keeps the first value (as boost does)
                    if attrs.child(&key).is_none() {
                        attrs.children.push((
                            key,
                            XmlTree {
                                data: value,
                                children: Vec::new(),
                            },
                        ));
                    }
                }
                None => return Err(parse_error("unexpected end of document", self.pos)),
            }
        }
        self.finish_attrs(&mut node, attrs);
        self.content(&mut node, &name)?;
        Ok((name, node))
    }

    fn finish_attrs(&self, node: &mut XmlTree, attrs: XmlTree) {
        if !attrs.children.is_empty() {
            node.children.insert(0, ("<xmlattr>".to_string(), attrs));
        }
    }

    /// The content of an element up to and including its end tag.
    fn content(&mut self, node: &mut XmlTree, name: &str) -> Result<()> {
        loop {
            if self.pos >= self.s.len() {
                return Err(parse_error("unexpected end of document", self.pos));
            }
            if self.starts_with("</") {
                self.pos += 2;
                let end = self.name()?;
                self.skip_ws();
                if end != name || self.s.get(self.pos) != Some(&b'>') {
                    return Err(parse_error("mismatched end tag", self.pos));
                }
                self.pos += 1;
                return Ok(());
            }
            if self.starts_with("<!--") {
                self.pos += 4;
                let start = self.pos;
                self.skip_past("-->")?;
                let text = String::from_utf8_lossy(&self.s[start..self.pos - 3]).into_owned();
                node.children.push((
                    "<xmlcomment>".to_string(),
                    XmlTree {
                        data: text,
                        children: Vec::new(),
                    },
                ));
            } else if self.starts_with("<![CDATA[") {
                self.pos += 9;
                let start = self.pos;
                self.skip_past("]]>")?;
                node.data
                    .push_str(&String::from_utf8_lossy(&self.s[start..self.pos - 3]));
            } else if self.starts_with("<?") {
                self.skip_past("?>")?;
            } else if self.s[self.pos] == b'<' {
                self.pos += 1;
                let child = self.element()?;
                node.children.push(child);
            } else {
                let start = self.pos;
                while self.pos < self.s.len() && self.s[self.pos] != b'<' {
                    self.pos += 1;
                }
                let raw = String::from_utf8_lossy(&self.s[start..self.pos]);
                node.data.push_str(&decode_entities(&raw, start)?);
            }
        }
    }
}

/// `read_xml`: the document as a tree whose only child is the root element.
pub fn read_xml(text: &str) -> Result<XmlTree> {
    let mut p = Parser {
        s: text.as_bytes(),
        pos: 0,
    };
    let mut tree = XmlTree::default();
    let mut found_root = false;
    loop {
        p.skip_ws();
        if p.pos >= p.s.len() {
            break;
        }
        if p.starts_with("<?") {
            p.skip_past("?>")?;
        } else if p.starts_with("<!--") {
            p.skip_past("-->")?;
        } else if p.starts_with("<!") {
            p.skip_past(">")?;
        } else if p.s[p.pos] == b'<' {
            if found_root {
                return Err(parse_error("more than one root element", p.pos));
            }
            p.pos += 1;
            let root = p.element()?;
            tree.children.push(root);
            found_root = true;
        } else {
            return Err(parse_error(
                "unexpected text outside the root element",
                p.pos,
            ));
        }
    }
    if !found_root {
        return Err(parse_error("no root element", 0));
    }
    Ok(tree)
}

impl XmlTree {
    /// `put("<xmlattr>.name", value)`: sets an attribute, replacing a
    /// previous value.
    pub fn put_attr(&mut self, name: &str, value: impl Into<String>) {
        let value = value.into();
        let attrs = match self.children.iter().position(|(k, _)| k == "<xmlattr>") {
            Some(p) => &mut self.children[p].1,
            None => {
                self.children
                    .insert(0, ("<xmlattr>".into(), XmlTree::default()));
                &mut self.children[0].1
            }
        };
        match attrs.children.iter_mut().find(|(k, _)| k == name) {
            Some((_, v)) => v.data = value,
            None => attrs.children.push((
                name.to_string(),
                XmlTree {
                    data: value,
                    children: Vec::new(),
                },
            )),
        }
    }

    /// `add_child(name, tree)`
    pub fn add_child(&mut self, name: &str, tree: XmlTree) {
        self.children.push((name.to_string(), tree));
    }

    /// `put_child(name, tree)`: replaces the first child of that name.
    pub fn put_child(&mut self, name: &str, tree: XmlTree) {
        match self.children.iter_mut().find(|(k, _)| k == name) {
            Some((_, v)) => *v = tree,
            None => self.add_child(name, tree),
        }
    }

    /// A node holding only the given text.
    pub fn with_data(data: impl Into<String>) -> XmlTree {
        XmlTree {
            data: data.into(),
            children: Vec::new(),
        }
    }
}

fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            c => out.push(c),
        }
    }
    out
}

fn write_element(out: &mut String, name: &str, tree: &XmlTree, indent: usize, pretty: bool) {
    let pad = |out: &mut String, n: usize| {
        if pretty {
            out.push_str(&"\t".repeat(n));
        }
    };
    pad(out, indent);
    out.push('<');
    out.push_str(name);
    if let Some((_, attrs)) = tree.children.iter().find(|(k, _)| k == "<xmlattr>") {
        for (k, v) in &attrs.children {
            out.push_str(&format!(" {k}=\"{}\"", escape(&v.data)));
        }
    }
    let kids: Vec<&(String, XmlTree)> = tree
        .children
        .iter()
        .filter(|(k, _)| k != "<xmlattr>")
        .collect();
    if kids.is_empty() && tree.data.is_empty() {
        out.push_str("/>");
    } else if kids.is_empty() {
        out.push('>');
        out.push_str(&escape(&tree.data));
        out.push_str(&format!("</{name}>"));
    } else {
        out.push('>');
        if pretty {
            out.push('\n');
        }
        for (k, v) in kids {
            write_element(out, k, v, indent + 1, pretty);
        }
        pad(out, indent);
        out.push_str(&format!("</{name}>"));
    }
    if pretty {
        out.push('\n');
    }
}

/// `write_xml`: the declaration and the elements below `root`.
pub fn write_xml(root: &XmlTree, pretty: bool) -> String {
    let mut out = if pretty {
        "<?xml version=\"1.0\" encoding=\"windows-1252\"?>\n".to_string()
    } else {
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n".to_string()
    };
    for (k, v) in &root.children {
        write_element(&mut out, k, v, 0, pretty);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn attributes_text_and_paths() {
        let t = read_xml(
            "<?xml version=\"1.0\"?><!-- c --><a x='1'><b y=\"&lt;2&amp;\">hi<![CDATA[ <there>]]></b><b/><c>1</c></a>",
        )
        .unwrap();
        let a = t.get_child("a").unwrap();
        assert_eq!(a.attr("x", ""), "1");
        assert_eq!(a.get_str("b.<xmlattr>.y", ""), "<2&");
        assert_eq!(a.get_str("b", ""), "hi <there>");
        assert_eq!(a.children_named("b").count(), 2);
        assert_eq!(a.get_str("c", "no"), "1");
        assert_eq!(a.get_str("d", "no"), "no");
    }

    #[test]
    fn errors() {
        assert!(read_xml("").is_err());
        assert!(read_xml("<a><b></a>").is_err());
        assert!(read_xml("<a x=1/>").is_err());
        assert!(read_xml("<a/><b/>").is_err());
    }
}
