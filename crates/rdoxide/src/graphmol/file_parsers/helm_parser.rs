//! HELM reader (peptides and RNA/DNA with the basic connection types).
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/SequenceParsers.cpp`
//! (`HELMToMol`, 2026.09.1).
//!
//! Copyright (C) 2015,2016 Greg Landrum and NextMove Software.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::HashMap;

use super::sequence_parser::{AMINO_ACID_CODES, Builder};
use super::sequence_tables::HELM_PEPTIDE_MONOMERS;
use crate::Result;
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::RWMol;

/// The attachment points of one monomer of a chain (`HELMMonomer`).
#[derive(Clone, Copy, Default)]
struct Monomer {
    r1: Option<u32>,
    r2: Option<u32>,
    r3: Option<u32>,
    /// the terminal OXT that a later connection replaces
    oxt: Option<u32>,
}

/// A cursor over the HELM text; reads past the end give NUL like the C string.
#[derive(Clone, Copy)]
struct Cursor<'a> {
    text: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn at(&self, offset: usize) -> u8 {
        self.text.get(self.pos + offset).copied().unwrap_or(0)
    }

    fn cur(&self) -> u8 {
        self.at(0)
    }

    fn advance(&mut self, n: usize) {
        self.pos += n;
    }

    fn starts_with(&self, s: &str) -> bool {
        self.text[self.pos.min(self.text.len())..].starts_with(s.as_bytes())
    }

    fn skip_digits(&mut self) {
        while self.cur().is_ascii_digit() {
            self.pos += 1;
        }
    }

    /// `[1-9][0-9]*`
    fn number(&mut self) -> Option<u32> {
        if !(b'1'..=b'9').contains(&self.cur()) {
            return None;
        }
        let mut n = 0u32;
        while self.cur().is_ascii_digit() {
            n = n
                .wrapping_mul(10)
                .wrapping_add(u32::from(self.cur() - b'0'));
            self.pos += 1;
        }
        Some(n)
    }

    /// `PEPTIDE[1-9][0-9]*` or `RNA[1-9][0-9]*` style id: the prefix, a
    /// non-zero digit, then more digits. Returns the id text.
    fn polymer_id(&mut self, prefix: &str) -> Option<&'a str> {
        let start = self.pos;
        if !self.starts_with(prefix) || !(b'1'..=b'9').contains(&self.at(prefix.len())) {
            return None;
        }
        self.advance(prefix.len() + 1);
        self.skip_digits();
        std::str::from_utf8(&self.text[start..self.pos]).ok()
    }
}

/// `IsHELMMonomerIDChar`
fn is_monomer_id_char(c: u8) -> bool {
    c.is_ascii_alphanumeric()
}

/// `GetHELMOneLetterCode`
fn one_letter_code(c: u8) -> Option<&'static str> {
    AMINO_ACID_CODES
        .iter()
        .find(|(u, _, _)| *u == c)
        .map(|(_, l, _)| *l)
}

/// `LookupHELMPeptideMonomer`
fn peptide_monomer(id: &str) -> Option<&'static str> {
    HELM_PEPTIDE_MONOMERS
        .iter()
        .find(|(name, _)| *name == id)
        .map(|(_, code)| *code)
}

/// `ParseHELMPeptide`: builds the chain whose text follows the `{`; returns
/// the position of the closing `}` and the monomers.
fn parse_peptide(
    b: &mut Builder,
    cur: &mut Cursor<'_>,
    chain: &str,
) -> Result<Option<Vec<Monomer>>> {
    let mut seq: Vec<Monomer> = Vec::new();
    if cur.cur() == b'}' {
        return Ok(Some(seq));
    }
    b.start_chain(chain);

    if cur.starts_with("[ac]") {
        if cur.at(4) != b'.' {
            return Ok(None);
        }
        b.info.residue_number = -2;
        let [r1, r2, r3] = b.amino_acid("ACE")?;
        seq.push(Monomer {
            r1,
            r2,
            r3,
            oxt: None,
        });
        b.info.residue_number = 0;
        cur.advance(5);
    }

    loop {
        let name = if cur.cur() == b'[' {
            cur.advance(1);
            let start = cur.pos;
            while is_monomer_id_char(cur.cur()) {
                cur.advance(1);
            }
            let id = std::str::from_utf8(&cur.text[start..cur.pos]).unwrap_or_default();
            if cur.cur() != b']' {
                return Ok(None);
            }
            peptide_monomer(id)
        } else {
            one_letter_code(cur.cur())
        };
        let Some(name) = name else {
            return Ok(None);
        };
        cur.advance(1);

        let [r1, r2, r3] = b.amino_acid(name)?;
        if let (Some(last), Some(r1)) = (seq.last_mut(), r1) {
            if let Some(prev_r2) = last.r2.take() {
                b.bond(prev_r2, r1, 1)?;
            }
        }
        seq.push(Monomer {
            r1,
            r2,
            r3,
            oxt: None,
        });

        match cur.cur() {
            b'.' => {
                if cur.starts_with(".[am]}") {
                    let Some(last) = seq.last_mut().and_then(|m| m.r2.take()) else {
                        return Ok(None);
                    };
                    b.info.residue_number += 1;
                    b.info.pdb_mut().is_hetero_atom = true;
                    b.info.residue_name = "NH2".into();
                    let n = b.atom(" N  ")?;
                    b.bond(last, n, 1)?;
                    seq.push(Monomer::default());
                    cur.advance(5);
                    return Ok(Some(seq));
                }
                cur.advance(1);
            }
            b'}' => {
                let Some(last) = seq.last_mut() else {
                    return Ok(None);
                };
                let Some(r2) = last.r2 else {
                    return Ok(None);
                };
                let oxt = b.atom(" OXT")?;
                b.bond(r2, oxt, 1)?;
                last.oxt = Some(oxt);
                return Ok(Some(seq));
            }
            _ => return Ok(None),
        }
    }
}

/// The residue code and length of a nucleotide such as `R(A)` or `[dR](T)`.
fn nucleotide(cur: &Cursor<'_>) -> Option<(&'static str, usize)> {
    let (offset, prefix) = if cur.cur() == b'R' && cur.at(1) == b'(' {
        (2, "  ")
    } else if cur.starts_with("[dR](") {
        (5, " D")
    } else {
        return None;
    };
    let base = cur.at(offset);
    if cur.at(offset + 1) != b')' {
        return None;
    }
    let code = match (prefix, base) {
        ("  ", b'A') => "  A",
        ("  ", b'C') => "  C",
        ("  ", b'G') => "  G",
        ("  ", b'T') => "  T",
        ("  ", b'U') => "  U",
        (_, b'A') => " DA",
        (_, b'C') => " DC",
        (_, b'G') => " DG",
        (_, b'T') => " DT",
        (_, b'U') => " DU",
        _ => return None,
    };
    Some((code, offset + 2))
}

/// `ParseHELMNucleic`; `Ok(false)` for malformed text.
fn parse_nucleic(b: &mut Builder, cur: &mut Cursor<'_>, chain: &str) -> Result<bool> {
    if cur.cur() == b'}' {
        return Ok(true);
    }
    b.start_chain(chain);
    let mut p_cap5 = false;
    let mut prev: Option<u32> = None;

    if cur.cur() == b'P' {
        p_cap5 = true;
        cur.advance(1);
        if cur.cur() == b'.' {
            cur.advance(1);
        }
    }

    loop {
        let Some((name, len)) = nucleotide(cur) else {
            return Ok(false);
        };
        cur.advance(len);
        let [r1, r2, _] = b.nucleic_acid(name, p_cap5)?;
        if let (Some(p), Some(r1)) = (prev, r1) {
            b.bond(p, r1, 1)?;
        }
        prev = r2;

        if cur.cur() == b'}' {
            return Ok(true);
        }
        if cur.cur() == b'.' {
            cur.advance(1);
        }
        if cur.cur() != b'P' {
            return Ok(false);
        }
        cur.advance(1);
        if cur.cur() == b'}' {
            if let Some(p) = prev {
                b.p_cap3(p)?;
            }
            return Ok(true);
        }
        p_cap5 = true;
        if cur.cur() == b'.' {
            cur.advance(1);
        }
    }
}

/// One side of a connection section entry: `PEPTIDEn,PEPTIDEm,i:Ra-j:Rb`.
struct Connection<'a> {
    id1: &'a str,
    id2: &'a str,
    res1: u32,
    res1r: u32,
    res2: u32,
    res2r: u32,
}

/// Parses `PEPTIDEn,PEPTIDEm,i:Ra-j:Rb`.
fn parse_connection<'a>(cur: &mut Cursor<'a>) -> Option<Connection<'a>> {
    let id1 = cur.polymer_id("PEPTIDE")?;
    if cur.cur() != b',' {
        return None;
    }
    cur.advance(1);
    let id2 = cur.polymer_id("PEPTIDE")?;
    if cur.cur() != b',' {
        return None;
    }
    cur.advance(1);

    let res1 = cur.number()?;
    let res1r = group(cur)?;
    if cur.cur() != b'-' {
        return None;
    }
    cur.advance(1);
    let res2 = cur.number()?;
    let res2r = group(cur)?;
    Some(Connection {
        id1,
        id2,
        res1,
        res1r,
        res2,
        res2r,
    })
}

/// `:R[1-9]`
fn group(cur: &mut Cursor<'_>) -> Option<u32> {
    if cur.cur() == b':' && cur.at(1) == b'R' && (b'1'..=b'9').contains(&cur.at(2)) {
        let g = u32::from(cur.at(2) - b'0');
        cur.advance(3);
        Some(g)
    } else {
        None
    }
}

/// The chain letter for a polymer id digit (`A` for 1, `B` for 2, ...).
fn chain_letter(digit: u8) -> String {
    ((b'A' + (digit - b'1')) as char).to_string()
}

/// Applies one connection to the chains; `Ok(false)` if it is not possible.
fn connect(
    b: &mut Builder,
    seqs: &mut HashMap<String, Vec<Monomer>>,
    removed: &mut Vec<u32>,
    c: &Connection<'_>,
) -> Result<bool> {
    if c.res1 < 1 || c.res2 < 1 {
        return Ok(false);
    }
    let in_range = |id: &str, res: u32| seqs.get(id).is_some_and(|s| res as usize <= s.len());
    if !in_range(c.id1, c.res1) || !in_range(c.id2, c.res2) {
        return Ok(false);
    }
    let (i1, i2) = (c.res1 as usize - 1, c.res2 as usize - 1);
    let m1 = seqs[c.id1][i1];
    let m2 = seqs[c.id2][i2];

    match (c.res1r, c.res2r) {
        // side chain to side chain
        (3, 3) => {
            let (Some(src), Some(dst)) = (m1.r3, m2.r3) else {
                return Ok(false);
            };
            if src == dst {
                return Ok(false);
            }
            b.bond(src, dst, 1)?;
            seqs.get_mut(c.id1).expect("checked")[i1].r3 = None;
            seqs.get_mut(c.id2).expect("checked")[i2].r3 = None;
        }
        // N terminus of the first to C terminus of the second (head to tail)
        (1, 2) | (2, 1) => {
            let (n_term, c_term, n_id, c_id, n_i, c_i) = if c.res1r == 1 {
                (m1.r1, m2.r2, c.id1, c.id2, i1, i2)
            } else {
                (m2.r1, m1.r2, c.id2, c.id1, i2, i1)
            };
            let (Some(src), Some(dst), Some(oxt)) =
                (n_term, c_term, if c.res1r == 1 { m2.oxt } else { m1.oxt })
            else {
                return Ok(false);
            };
            if src == dst {
                return Ok(false);
            }
            // upstream removes the OXT now; the removal is deferred so that
            // the atom indices held by the monomers stay valid
            removed.push(oxt);
            if c.res1r == 1 {
                b.bond(src, dst, 1)?;
            } else {
                b.bond(dst, src, 1)?;
            }
            seqs.get_mut(n_id).expect("checked")[n_i].r1 = None;
            seqs.get_mut(c_id).expect("checked")[c_i].r2 = None;
        }
        _ => return Ok(false),
    }
    Ok(true)
}

/// `ParseHELM`
fn parse_helm(b: &mut Builder, text: &[u8], removed: &mut Vec<u32>) -> Result<bool> {
    let mut cur = Cursor { text, pos: 0 };
    let mut seqs: HashMap<String, Vec<Monomer>> = HashMap::new();

    // simple polymers
    loop {
        let start = cur.pos;
        if cur.starts_with("PEPTIDE") {
            let Some(id) = cur.polymer_id("PEPTIDE") else {
                return Ok(false);
            };
            if cur.cur() != b'{' {
                return Ok(false);
            }
            cur.advance(1);
            let chain = chain_letter(text[start + 7]);
            let Some(seq) = parse_peptide(b, &mut cur, &chain)? else {
                return Ok(false);
            };
            if cur.cur() != b'}' {
                return Ok(false);
            }
            seqs.insert(id.to_string(), seq);
            cur.advance(1);
        } else if cur.starts_with("RNA") {
            if cur.polymer_id("RNA").is_none() || cur.cur() != b'{' {
                return Ok(false);
            }
            cur.advance(1);
            let chain = chain_letter(text[start + 3]);
            if !parse_nucleic(b, &mut cur, &chain)? || cur.cur() != b'}' {
                return Ok(false);
            }
            cur.advance(1);
        } else {
            return Ok(false);
        }

        match cur.cur() {
            b'$' => break,
            0 => return Ok(true),
            b'|' => cur.advance(1),
            _ => return Ok(false),
        }
    }
    cur.advance(1);

    if cur.starts_with("$$$") {
        return Ok(true);
    }

    // connections
    loop {
        let Some(c) = parse_connection(&mut cur) else {
            return Ok(false);
        };
        if !connect(b, &mut seqs, removed, &c)? {
            return Ok(false);
        }
        match cur.cur() {
            b'$' => break,
            b'|' => cur.advance(1),
            _ => return Ok(false),
        }
    }
    cur.advance(1);
    Ok(cur.cur() == b'$' && cur.at(1) == b'$')
}

/// `HELMToMol`: `Ok(None)` for text that is not valid (simple) HELM.
pub fn helm_to_mol(helm: &str, sanitize: bool) -> Result<Option<RWMol>> {
    let text = helm.as_bytes();
    let text = &text[..text.iter().position(|&c| c == 0).unwrap_or(text.len())];
    if text.starts_with(b"$$$$") {
        return Ok(Some(RWMol::new()));
    }
    let mut builder = Builder::new();
    let mut removed = Vec::new();
    if !parse_helm(&mut builder, text, &mut removed)? {
        return Ok(None);
    }
    let mut mol = builder.mol;
    // remove the replaced OXT atoms, highest index first
    removed.sort_unstable_by(|a, b| b.cmp(a));
    removed.dedup();
    for idx in removed {
        mol.remove_atom(idx)?;
    }
    if sanitize {
        sanitize_mol(&mut mol, SANITIZE_ALL).map_err(|(_, e)| e)?;
    }
    Ok(Some(mol))
}
