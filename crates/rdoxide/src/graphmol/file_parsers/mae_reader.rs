//! Reader of the Maestro (`.mae`) format: the part of Schrodinger's
//! `maeparser` that `MaeMolSupplier` uses.
//!
//! A file is a sequence of outer blocks (`f_m_ct { ... }`), each holding
//! typed properties (`b_`, `i_`, `r_`, `s_` keys, first the keys, then `:::`,
//! then the values), nested blocks and indexed blocks (`m_atom[3] { ... }`)
//! with one row of values per index. As in `maeparser`, the values of the
//! indexed blocks are only scanned when the block is read and converted when
//! they are first asked for, so malformed values show up then.
//!
//! The error messages (including the `Line N, column M:` prefix) are those of
//! `maeparser`. The grammar was recovered from its behaviour, as its source is
//! not part of RDKit.
//!
//! This is an independent implementation written from the observed behaviour
//! of `maeparser` and the public file format; it contains no `maeparser` code.

use std::collections::BTreeMap;

/// A `maeparser` exception message.
pub type MaeResult<T> = std::result::Result<T, String>;

/// Property names used by the supplier and writer (`MaeConstants.hpp`).
pub mod constants {
    pub const CT_BLOCK: &str = "f_m_ct";
    pub const ATOM_BLOCK: &str = "m_atom";
    pub const BOND_BLOCK: &str = "m_bond";
    pub const CT_TITLE: &str = "s_m_title";
    pub const CT_CHIRALITY_PROP_PREFIX: &str = "s_st_Chirality_";
    pub const CT_PSEUDOCHIRALITY_PROP_PREFIX: &str = "s_st_AtomNumChirality_";
    pub const CT_EZ_PROP_PREFIX: &str = "s_st_EZ_";
    pub const ATOM_ATOMIC_NUM: &str = "i_m_atomic_number";
    pub const ATOM_X_COORD: &str = "r_m_x_coord";
    pub const ATOM_Y_COORD: &str = "r_m_y_coord";
    pub const ATOM_Z_COORD: &str = "r_m_z_coord";
    pub const ATOM_FORMAL_CHARGE: &str = "i_m_formal_charge";
    pub const BOND_ATOM_1: &str = "i_m_from";
    pub const BOND_ATOM_2: &str = "i_m_to";
    pub const BOND_ORDER: &str = "i_m_order";
}

/// A value of a row of an indexed block, before conversion.
#[derive(Debug, Clone, PartialEq)]
struct Token {
    text: String,
    quoted: bool,
}

/// An indexed block whose rows are converted on first use.
#[derive(Debug, Clone)]
pub struct IndexedBlock {
    names: Vec<String>,
    size: usize,
    /// `size` rows of `names.len() + 1` tokens (the first is the index).
    tokens: Vec<Token>,
}

/// A column of an indexed block: `None` is the undefined value (`<>`).
pub type Column<T> = Vec<Option<T>>;

/// The converted columns of an [`IndexedBlock`].
#[derive(Debug, Default)]
pub struct ParsedIndexedBlock {
    pub size: usize,
    pub bools: BTreeMap<String, Column<bool>>,
    pub ints: BTreeMap<String, Column<i32>>,
    pub reals: BTreeMap<String, Column<f64>>,
    pub strings: BTreeMap<String, Column<String>>,
}

impl ParsedIndexedBlock {
    fn find<'a, T>(map: &'a BTreeMap<String, Column<T>>, name: &str) -> MaeResult<&'a Column<T>> {
        map.get(name)
            .ok_or_else(|| "map::at:  key not found".to_string())
    }
    pub fn int_property(&self, name: &str) -> MaeResult<&Column<i32>> {
        Self::find(&self.ints, name)
    }
    pub fn real_property(&self, name: &str) -> MaeResult<&Column<f64>> {
        Self::find(&self.reals, name)
    }
    pub fn string_property(&self, name: &str) -> Option<&Column<String>> {
        self.strings.get(name)
    }
}

impl IndexedBlock {
    pub fn size(&self) -> usize {
        self.size
    }

    /// Converts the rows (`IndexedBlock` properties are typed by the first
    /// letter of their names).
    pub fn parse(&self) -> MaeResult<ParsedIndexedBlock> {
        let width = self.names.len() + 1;
        let mut res = ParsedIndexedBlock {
            size: self.size,
            ..Default::default()
        };
        for name in &self.names {
            let n = self.size;
            match name.as_bytes()[0] {
                b'b' => {
                    res.bools.insert(name.clone(), vec![None; n]);
                }
                b'i' => {
                    res.ints.insert(name.clone(), vec![None; n]);
                }
                b'r' => {
                    res.reals.insert(name.clone(), vec![None; n]);
                }
                _ => {
                    res.strings.insert(name.clone(), vec![None; n]);
                }
            }
        }
        for row in 0..self.size {
            for (col, name) in self.names.iter().enumerate() {
                let tok = &self.tokens[row * width + col + 1];
                let undefined = !tok.quoted && tok.text == "<>";
                match name.as_bytes()[0] {
                    b'b' => {
                        if !undefined {
                            let v = match tok.text.as_str() {
                                "0" => false,
                                "1" => true,
                                _ => return Err("Bogus bool.".into()),
                            };
                            res.bools.get_mut(name).expect("column")[row] = Some(v);
                        }
                    }
                    b'i' => {
                        if !undefined {
                            res.ints.get_mut(name).expect("column")[row] =
                                Some(parse_lazy_int(&tok.text)?);
                        }
                    }
                    b'r' => {
                        if !undefined {
                            let v = tok
                                .text
                                .parse::<f64>()
                                .map_err(|_| "Bad floating point representation.".to_string())?;
                            res.reals.get_mut(name).expect("column")[row] = Some(v);
                        }
                    }
                    _ => {
                        if !undefined {
                            res.strings.get_mut(name).expect("column")[row] =
                                Some(tok.text.clone());
                        }
                    }
                }
            }
        }
        Ok(res)
    }
}

/// `[-]digits` with wrapping 32-bit arithmetic (as `maeparser` does).
fn parse_lazy_int(text: &str) -> MaeResult<i32> {
    let (neg, digits) = match text.strip_prefix('-') {
        Some(d) => (true, d),
        None => (false, text),
    };
    let mut v: i32 = 0;
    for c in digits.bytes() {
        if !c.is_ascii_digit() {
            return Err("Unexpected character in integer.".into());
        }
        v = v.wrapping_mul(10).wrapping_add(i32::from(c - b'0'));
    }
    Ok(if neg { v.wrapping_neg() } else { v })
}

/// A (possibly nested) block with its typed properties.
#[derive(Debug, Clone, Default)]
pub struct Block {
    pub name: String,
    pub bools: BTreeMap<String, bool>,
    pub ints: BTreeMap<String, i32>,
    pub reals: BTreeMap<String, f64>,
    pub strings: BTreeMap<String, String>,
    pub indexed: BTreeMap<String, IndexedBlock>,
    pub blocks: BTreeMap<String, Block>,
}

impl Block {
    pub fn indexed_block(&self, name: &str) -> MaeResult<&IndexedBlock> {
        self.indexed
            .get(name)
            .ok_or_else(|| format!("Indexed block not found: {name}"))
    }
}

/// Reads blocks one after another from the text of a file.
#[derive(Debug, Clone)]
pub struct Reader {
    data: Vec<u8>,
    pos: usize,
}

fn is_ws(c: u8) -> bool {
    c.is_ascii_whitespace()
}

impl Reader {
    pub fn new(data: Vec<u8>) -> Self {
        Reader { data, pos: 0 }
    }

    /// Back to the start of the data.
    pub fn rewind(&mut self) {
        self.pos = 0;
    }

    fn err<T>(&self, pos: usize, msg: &str) -> MaeResult<T> {
        let pos = pos.min(self.data.len());
        let line = 1 + self.data[..pos].iter().filter(|&&c| c == b'\n').count();
        let line_start = self.data[..pos]
            .iter()
            .rposition(|&c| c == b'\n')
            .map_or(0, |i| i + 1);
        Err(format!(
            "Line {line}, column {}: {msg}\n",
            pos - line_start + 1
        ))
    }

    fn peek(&self) -> Option<u8> {
        self.data.get(self.pos).copied()
    }

    /// Skips whitespace and `# ... #` comments.
    fn skip_ws(&mut self) -> MaeResult<()> {
        while let Some(c) = self.peek() {
            if is_ws(c) {
                self.pos += 1;
            } else if c == b'#' {
                match self.data[self.pos + 1..].iter().position(|&c| c == b'#') {
                    Some(n) => self.pos += n + 2,
                    None => return self.err(self.data.len(), "Unterminated comment."),
                }
            } else {
                break;
            }
        }
        Ok(())
    }

    /// The next block named `outer_name`, skipping the others; `None` at the
    /// end of the data.
    pub fn next(&mut self, outer_name: &str) -> MaeResult<Option<Block>> {
        loop {
            self.skip_ws()?;
            if self.peek().is_none() {
                return Ok(None);
            }
            let name = self.outer_block_name()?;
            self.skip_ws()?;
            if self.peek() != Some(b'{') {
                return self.err(self.pos, "Missing '{' for outer block.");
            }
            self.pos += 1;
            let mut block = self.block_body()?;
            block.name = name;
            if block.name == outer_name {
                return Ok(Some(block));
            }
        }
    }

    /// `(f|p)_<author>_<name>`, or nothing for the unnamed header block.
    fn outer_block_name(&mut self) -> MaeResult<String> {
        if self.peek() == Some(b'{') {
            return Ok(String::new());
        }
        let start = self.pos;
        let bad = "Bad format for outer block name; must be (f|p)_<author>_<name>.";
        if !matches!(self.peek(), Some(b'f' | b'p')) {
            return self.err(self.pos, bad);
        }
        self.pos += 1;
        if self.peek() != Some(b'_') {
            return self.err(self.pos, bad);
        }
        self.pos += 1;
        self.author_and_name(bad)?;
        Ok(String::from_utf8_lossy(&self.data[start..self.pos]).into_owned())
    }

    /// `<author>_<name>`: both parts nonempty, ended by whitespace, `{`
    /// or (for the name) `[`.
    fn author_and_name(&mut self, bad: &str) -> MaeResult<()> {
        let author = self.pos;
        while self
            .peek()
            .is_some_and(|c| c != b'_' && !is_ws(c) && c != b'{')
        {
            self.pos += 1;
        }
        if self.pos == author || self.peek() != Some(b'_') {
            return self.err(self.pos, bad);
        }
        self.pos += 1;
        let name = self.pos;
        while self
            .peek()
            .is_some_and(|c| !is_ws(c) && c != b'{' && c != b'[')
        {
            self.pos += 1;
        }
        if self.pos == name {
            return self.err(self.pos, bad);
        }
        Ok(())
    }

    /// The property keys of a block, up to `:::`.
    fn keys(&mut self) -> MaeResult<Vec<String>> {
        let mut keys = Vec::new();
        loop {
            self.skip_ws()?;
            if self.data[self.pos..].starts_with(b":::") {
                self.pos += 3;
                return Ok(keys);
            }
            let Some(c) = self.peek() else {
                return self.err(self.pos, "Missing property key.");
            };
            let start = self.pos;
            let bad = "Bad format for property; must be (b|i|r|s)_<author>_<name>.";
            if !matches!(c, b'b' | b'i' | b'r' | b's') {
                return self.err(self.pos, bad);
            }
            self.pos += 1;
            if self.peek() != Some(b'_') {
                return self.err(self.pos, bad);
            }
            self.pos += 1;
            // the name may contain anything but whitespace
            let author = self.pos;
            while self.peek().is_some_and(|c| c != b'_' && !is_ws(c)) {
                self.pos += 1;
            }
            if self.pos == author || self.peek() != Some(b'_') {
                return self.err(self.pos, bad);
            }
            self.pos += 1;
            let name = self.pos;
            while self.peek().is_some_and(|c| !is_ws(c)) {
                self.pos += 1;
            }
            if self.pos == name {
                return self.err(self.pos, bad);
            }
            keys.push(String::from_utf8_lossy(&self.data[start..self.pos]).into_owned());
        }
    }

    /// A quoted string (the opening quote is under the cursor): a backslash
    /// takes the next character literally.
    fn quoted(&mut self) -> MaeResult<String> {
        self.pos += 1;
        let mut out = Vec::new();
        loop {
            match self.peek() {
                None => return self.err(self.data.len(), "Unterminated quoted string at EOF."),
                Some(b'"') => {
                    self.pos += 1;
                    return Ok(String::from_utf8_lossy(&out).into_owned());
                }
                Some(b'\\') => {
                    self.pos += 1;
                    match self.peek() {
                        Some(c) => out.push(c),
                        None => {
                            return self.err(self.data.len(), "Unterminated quoted string at EOF.");
                        }
                    }
                    self.pos += 1;
                }
                Some(c) => {
                    out.push(c);
                    self.pos += 1;
                }
            }
        }
    }

    /// A token of an indexed block row: quoted, or up to whitespace.
    fn token(&mut self) -> MaeResult<Token> {
        self.skip_ws()?;
        match self.peek() {
            None => self.err(self.data.len(), "Unexpected EOF in indexed block values."),
            Some(b'"') => Ok(Token {
                text: self.quoted()?,
                quoted: true,
            }),
            Some(_) => {
                let start = self.pos;
                while self.peek().is_some_and(|c| !is_ws(c)) {
                    self.pos += 1;
                }
                Ok(Token {
                    text: String::from_utf8_lossy(&self.data[start..self.pos]).into_owned(),
                    quoted: false,
                })
            }
        }
    }

    /// The end of an eagerly parsed value: whitespace or the end of the data.
    fn value_end(&self, msg: &str) -> MaeResult<()> {
        match self.peek() {
            Some(c) if !is_ws(c) => self.err(self.pos, msg),
            _ => Ok(()),
        }
    }

    fn eager_bool(&mut self) -> MaeResult<bool> {
        let v = match self.peek() {
            Some(b'0') => false,
            Some(b'1') => true,
            _ => return self.err(self.pos, "Unexpected character for boolean value."),
        };
        self.pos += 1;
        self.value_end("Unexpected character for boolean value.")?;
        Ok(v)
    }

    fn eager_int(&mut self) -> MaeResult<i32> {
        let neg = self.peek() == Some(b'-');
        if neg {
            self.pos += 1;
        }
        let start = self.pos;
        let mut v: i32 = 0;
        while let Some(c) = self.peek().filter(u8::is_ascii_digit) {
            v = v.wrapping_mul(10).wrapping_add(i32::from(c - b'0'));
            self.pos += 1;
        }
        if self.pos == start {
            return if neg {
                self.err(self.pos, "Unexpected '-' in integer.")
            } else {
                self.err(self.pos, "Unexpected character.")
            };
        }
        self.value_end("Unexpected character.")?;
        Ok(if neg { v.wrapping_neg() } else { v })
    }

    fn eager_real(&mut self) -> MaeResult<f64> {
        let start = self.pos;
        let digits = |r: &mut Self| {
            let s = r.pos;
            while r.peek().is_some_and(|c| c.is_ascii_digit()) {
                r.pos += 1;
            }
            r.pos - s
        };
        if self.peek() == Some(b'-') {
            self.pos += 1;
        }
        let mut n = digits(self);
        if self.peek() == Some(b'.') {
            self.pos += 1;
            n += digits(self);
        }
        if n == 0 {
            return self.err(start, "Unexpected character in real number.");
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.pos += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.pos += 1;
            }
            if digits(self) == 0 {
                return self.err(self.pos, "Bad real number.");
            }
        }
        self.value_end("Unexpected character in real number.")?;
        std::str::from_utf8(&self.data[start..self.pos])
            .ok()
            .and_then(|s| s.parse().ok())
            .map_or_else(|| self.err(start, "Bad real number."), Ok)
    }

    fn eager_string(&mut self) -> MaeResult<String> {
        match self.peek() {
            None => self.err(self.pos, "Unexpected EOF."),
            Some(b'"') => self.quoted(),
            Some(_) => {
                let start = self.pos;
                while self.peek().is_some_and(|c| !is_ws(c)) {
                    self.pos += 1;
                }
                Ok(String::from_utf8_lossy(&self.data[start..self.pos]).into_owned())
            }
        }
    }

    /// The part of a block after its `{`.
    fn block_body(&mut self) -> MaeResult<Block> {
        let mut block = Block::default();
        let keys = self.keys()?;
        for key in keys {
            self.skip_ws()?;
            if self.peek().is_none() {
                return self.err(self.pos, "Unexpected EOF.");
            }
            match key.as_bytes()[0] {
                b'b' => {
                    let v = self.eager_bool()?;
                    block.bools.insert(key, v);
                }
                b'i' => {
                    let v = self.eager_int()?;
                    block.ints.insert(key, v);
                }
                b'r' => {
                    let v = self.eager_real()?;
                    block.reals.insert(key, v);
                }
                _ => {
                    let v = self.eager_string()?;
                    block.strings.insert(key, v);
                }
            }
        }
        loop {
            self.skip_ws()?;
            match self.peek() {
                None => return self.err(self.pos, "Missing '}' for block."),
                Some(b'}') => {
                    self.pos += 1;
                    return Ok(block);
                }
                Some(_) => self.sub_block(&mut block)?,
            }
        }
    }

    fn sub_block(&mut self, parent: &mut Block) -> MaeResult<()> {
        let start = self.pos;
        self.author_and_name("Bad format for block name; must be <author>_<name>.")?;
        let name = String::from_utf8_lossy(&self.data[start..self.pos]).into_owned();
        self.skip_ws()?;
        let mut size = None;
        if self.peek() == Some(b'[') {
            self.pos += 1;
            self.skip_ws()?;
            let digits = self.pos;
            let mut n: usize = 0;
            while let Some(c) = self.peek().filter(u8::is_ascii_digit) {
                n = n.wrapping_mul(10).wrapping_add(usize::from(c - b'0'));
                self.pos += 1;
            }
            if digits == self.pos {
                return self.err(self.pos, "Missing integer.");
            }
            self.skip_ws()?;
            if self.peek() != Some(b']') {
                return self.err(self.pos, "Bad block index; missing ']'.");
            }
            self.pos += 1;
            self.skip_ws()?;
            size = Some(n);
        }
        if self.peek() != Some(b'{') {
            return self.err(self.pos, "Missing '{' for block.");
        }
        self.pos += 1;
        match size {
            None => {
                let mut b = self.block_body()?;
                b.name = name.clone();
                parent.blocks.insert(name, b);
            }
            Some(size) => {
                let b = self.indexed_body(size)?;
                parent.indexed.insert(name, b);
            }
        }
        Ok(())
    }

    fn indexed_body(&mut self, size: usize) -> MaeResult<IndexedBlock> {
        let names = self.keys()?;
        let width = names.len() + 1;
        let mut tokens = Vec::new();
        for _ in 0..size.saturating_mul(width) {
            tokens.push(self.token()?);
        }
        self.skip_ws()?;
        if size > 0 || self.data[self.pos..].starts_with(b":::") {
            if !self.data[self.pos..].starts_with(b":::") {
                return self.err(self.pos, "Bad ':::' token.");
            }
            self.pos += 3;
        }
        self.skip_ws()?;
        if self.peek() != Some(b'}') {
            return self.err(self.pos, "Missing closing '}' for indexed block.");
        }
        self.pos += 1;
        Ok(IndexedBlock {
            names,
            size,
            tokens,
        })
    }
}
