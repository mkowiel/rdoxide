//! Port of RDKit `Code/GraphMol/GeneralizedSubstruct` (`XQMol.{h,cpp}`):
//! extended query molecules for generalized substructure searching.
//!
//! Serialization uses an rdoxide-specific format (a tag, the molecule pickles
//! and the tautomer-query indices), not upstream's Boost serialization, so
//! binary and JSON strings are not exchangeable with RDKit.

use std::sync::Arc;

use super::adjust_query::{AdjustQueryParameters, adjust_query_properties};
use super::fingerprints::pattern::pattern_fingerprint_mol;
use super::mol_bundle::MolBundle;
use super::mol_enumerator::enumerate_default;
use super::pickle::{flags, mol_from_pickle_stream, pickle_mol};
use super::romol::ROMol;
use super::substruct::{MatchVect, SubstructMatchParameters, substruct_match};
use super::tautomer_query::TautomerQuery;
use crate::datastructs::{BitVect, ExplicitBitVect};
use crate::{RDKitError, Result};

const MAGIC: &[u8; 4] = b"RDXQ";
const JSON_FORMAT: &str = "rdoxide-xqmol";

/// The value of the string field `key` in the flat JSON object written by
/// [`ExtendedQueryMol::to_json`] (its strings never contain escapes).
fn json_string_field<'a>(text: &'a str, key: &str) -> Option<&'a str> {
    let rest = &text[text.find(&format!("\"{key}\""))? + key.len() + 2..];
    let rest = rest
        .trim_start()
        .strip_prefix(':')?
        .trim_start()
        .strip_prefix('"')?;
    Some(&rest[..rest.find('"')?])
}

/// `GeneralizedSubstruct::ExtendedQueryMol`
#[derive(Debug, Clone)]
pub enum ExtendedQueryMol {
    Mol(Box<ROMol>),
    Bundle(MolBundle),
    Tautomer(TautomerQuery),
    TautomerBundle(Vec<TautomerQuery>),
}

struct Writer(Vec<u8>);

impl Writer {
    fn u32(&mut self, v: u32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn mol(&mut self, mol: &ROMol) -> Result<()> {
        let pkl = pickle_mol(mol, flags::ALL_PROPS)?;
        self.u32(pkl.len() as u32);
        self.0.extend_from_slice(&pkl);
        Ok(())
    }
    fn indices(&mut self, v: &[u32]) {
        self.u32(v.len() as u32);
        v.iter().for_each(|&i| self.u32(i));
    }
    fn tautomer_query(&mut self, tq: &TautomerQuery) -> Result<()> {
        self.mol(tq.template_molecule())?;
        self.u32(tq.tautomers().len() as u32);
        for t in tq.tautomers() {
            self.mol(t)?;
        }
        self.indices(tq.modified_atoms());
        self.indices(tq.modified_bonds());
        Ok(())
    }
}

struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl Reader<'_> {
    fn bad() -> RDKitError {
        RDKitError::Value("invalid ExtendedQueryMol data".into())
    }
    fn u32(&mut self) -> Result<u32> {
        let b = self
            .data
            .get(self.pos..self.pos + 4)
            .ok_or_else(Self::bad)?;
        self.pos += 4;
        Ok(u32::from_le_bytes(b.try_into().expect("4 bytes")))
    }
    fn mol(&mut self) -> Result<ROMol> {
        let len = self.u32()? as usize;
        let end = self.pos.checked_add(len).ok_or_else(Self::bad)?;
        let data = self.data.get(self.pos..end).ok_or_else(Self::bad)?;
        let mut p = 0;
        let mol = mol_from_pickle_stream(data, &mut p, flags::ALL_PROPS)?;
        self.pos = end;
        Ok(mol)
    }
    fn indices(&mut self) -> Result<Vec<u32>> {
        let n = self.u32()?;
        (0..n).map(|_| self.u32()).collect()
    }
    fn tautomer_query(&mut self) -> Result<TautomerQuery> {
        let template = self.mol()?;
        let n = self.u32()?;
        let tautomers = (0..n)
            .map(|_| self.mol().map(Arc::new))
            .collect::<Result<Vec<_>>>()?;
        let atoms = self.indices()?;
        let bonds = self.indices()?;
        Ok(TautomerQuery::new(tautomers, template, atoms, bonds))
    }
}

const TQ_MAGIC: &[u8; 4] = b"RDXT";

/// `TautomerQuery::toBinary` (the pickle of a `TautomerQuery`)
pub fn tautomer_query_to_binary(tq: &TautomerQuery) -> Result<Vec<u8>> {
    let mut w = Writer(TQ_MAGIC.to_vec());
    w.tautomer_query(tq)?;
    Ok(w.0)
}

/// `TautomerQuery(pickle)`
pub fn tautomer_query_from_binary(data: &[u8]) -> Result<TautomerQuery> {
    if data.len() < 4 || &data[..4] != TQ_MAGIC {
        return Err(RDKitError::Value("invalid TautomerQuery data".into()));
    }
    Reader { data, pos: 4 }.tautomer_query()
}

fn bundle_of(mols: Vec<ROMol>) -> Result<MolBundle> {
    let mut b = MolBundle::new();
    for m in mols {
        b.add_mol(Arc::new(m))?;
    }
    Ok(b)
}

impl ExtendedQueryMol {
    /// `ExtendedQueryMol(text, isJSON)`
    pub fn from_text(text: &[u8], is_json: bool) -> Result<Self> {
        if is_json {
            Self::from_json(std::str::from_utf8(text).map_err(|_| Reader::bad())?)
        } else {
            Self::from_binary(text)
        }
    }

    /// `toBinary`
    pub fn to_binary(&self) -> Result<Vec<u8>> {
        let mut w = Writer(MAGIC.to_vec());
        match self {
            Self::Mol(m) => {
                w.0.push(0);
                w.mol(m)?;
            }
            Self::Bundle(b) => {
                w.0.push(1);
                w.u32(b.size() as u32);
                for m in b.mols() {
                    w.mol(m)?;
                }
            }
            Self::Tautomer(tq) => {
                w.0.push(2);
                w.tautomer_query(tq)?;
            }
            Self::TautomerBundle(v) => {
                w.0.push(3);
                w.u32(v.len() as u32);
                for tq in v {
                    w.tautomer_query(tq)?;
                }
            }
        }
        Ok(w.0)
    }

    /// `initFromBinary`
    pub fn from_binary(data: &[u8]) -> Result<Self> {
        if data.len() < 5 || &data[..4] != MAGIC {
            return Err(Reader::bad());
        }
        let mut r = Reader { data, pos: 5 };
        match data[4] {
            0 => Ok(Self::Mol(Box::new(r.mol()?))),
            1 => {
                let n = r.u32()?;
                let mols = (0..n).map(|_| r.mol()).collect::<Result<Vec<_>>>()?;
                Ok(Self::Bundle(bundle_of(mols)?))
            }
            2 => Ok(Self::Tautomer(r.tautomer_query()?)),
            3 => {
                let n = r.u32()?;
                let v = (0..n)
                    .map(|_| r.tautomer_query())
                    .collect::<Result<Vec<_>>>()?;
                Ok(Self::TautomerBundle(v))
            }
            _ => Err(Reader::bad()),
        }
    }

    /// `toJSON`: the binary form, hex-encoded in a small JSON document.
    pub fn to_json(&self) -> Result<String> {
        let hex: String = self
            .to_binary()?
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        Ok(format!(
            "{{\"format\": \"{JSON_FORMAT}\", \"version\": 1, \"data\": \"{hex}\"}}"
        ))
    }

    /// `initFromJSON`
    pub fn from_json(text: &str) -> Result<Self> {
        let hex = json_string_field(text, "data")
            .filter(|_| json_string_field(text, "format") == Some(JSON_FORMAT))
            .ok_or_else(Reader::bad)?;
        if hex.len() % 2 != 0 || !hex.is_ascii() {
            return Err(Reader::bad());
        }
        let data = (0..hex.len() / 2)
            .map(|i| u8::from_str_radix(&hex[2 * i..2 * i + 2], 16).map_err(|_| Reader::bad()))
            .collect::<Result<Vec<_>>>()?;
        Self::from_binary(&data)
    }

    /// `patternFingerprintQuery`
    pub fn pattern_fingerprint_query(&self, fp_size: u32) -> Result<ExplicitBitVect> {
        match self {
            Self::Mol(m) => pattern_fingerprint_mol(m, fp_size, true),
            Self::Bundle(b) => bundle_pattern_fingerprint(b, fp_size),
            Self::Tautomer(tq) => tq.pattern_fingerprint_template(fp_size),
            // upstream keeps the template fingerprint of the first query only
            Self::TautomerBundle(v) => match v.first() {
                Some(tq) => tq.pattern_fingerprint_template(fp_size),
                None => Err(RDKitError::Value("empty tautomer bundle".into())),
            },
        }
    }
}

/// `PatternFingerprintMol(MolBundle, ...)`: the bits common to all the
/// molecules of the bundle.
fn bundle_pattern_fingerprint(b: &MolBundle, fp_size: u32) -> Result<ExplicitBitVect> {
    let mut res: Option<ExplicitBitVect> = None;
    for m in b.mols() {
        let fp = pattern_fingerprint_mol(m, fp_size, true)?;
        res = Some(match res {
            None => fp,
            Some(r) => r.and(&fp),
        });
    }
    res.ok_or_else(|| RDKitError::Value("empty molecule bundle".into()))
}

/// `GeneralizedSubstruct::SubstructMatch`
pub fn substruct_match_xq(
    mol: &ROMol,
    query: &ExtendedQueryMol,
    params: &SubstructMatchParameters,
) -> Result<Vec<MatchVect>> {
    match query {
        ExtendedQueryMol::Mol(q) => substruct_match(mol, q, params),
        ExtendedQueryMol::Bundle(b) => {
            for q in b.mols() {
                let res = substruct_match(mol, q, params)?;
                if !res.is_empty() {
                    return Ok(res);
                }
            }
            Ok(Vec::new())
        }
        ExtendedQueryMol::Tautomer(tq) => Ok(tq.substruct_of(mol, params)?.0),
        ExtendedQueryMol::TautomerBundle(v) => {
            for tq in v {
                let res = tq.substruct_of(mol, params)?.0;
                if !res.is_empty() {
                    return Ok(res);
                }
            }
            Ok(Vec::new())
        }
    }
}

/// `GeneralizedSubstruct::hasSubstructMatch`
pub fn has_substruct_match_xq(
    mol: &ROMol,
    query: &ExtendedQueryMol,
    params: &SubstructMatchParameters,
) -> Result<bool> {
    let mut params = params.clone();
    params.max_matches = 1;
    Ok(!substruct_match_xq(mol, query, &params)?.is_empty())
}

/// `GeneralizedSubstruct::createExtendedQueryMol`
pub fn create_extended_query_mol(
    mol: &ROMol,
    do_enumeration: bool,
    do_tautomers: bool,
    adjust_query_props: bool,
    params: &AdjustQueryParameters,
) -> Result<ExtendedQueryMol> {
    let mut bundle = if do_enumeration {
        enumerate_default(mol, 0)?
    } else {
        MolBundle::new()
    };
    if bundle.is_empty() {
        bundle.add_mol(Arc::new(mol.clone()))?;
    }
    let adjust = |m: &ROMol| -> Result<ROMol> {
        if adjust_query_props {
            adjust_query_properties(m, params)
        } else {
            Ok(m.clone())
        }
    };
    if bundle.size() == 1 {
        let lmol = adjust(&bundle.mols()[0])?;
        if do_tautomers {
            let tq = TautomerQuery::from_mol(&lmol, "")?;
            if tq.tautomers().len() > 1 {
                return Ok(ExtendedQueryMol::Tautomer(tq));
            }
        }
        return Ok(ExtendedQueryMol::Mol(Box::new(lmol)));
    }
    let mut adjusted = MolBundle::new();
    for m in bundle.mols() {
        adjusted.add_mol(Arc::new(adjust(m)?))?;
    }
    if do_tautomers {
        let mut had_tautomers = false;
        let mut tqs = Vec::new();
        for m in adjusted.mols() {
            let tq = TautomerQuery::from_mol(m, "")?;
            had_tautomers |= tq.tautomers().len() > 1;
            tqs.push(tq);
        }
        if had_tautomers {
            return Ok(ExtendedQueryMol::TautomerBundle(tqs));
        }
    }
    Ok(ExtendedQueryMol::Bundle(adjusted))
}

/// `GeneralizedSubstruct::patternFingerprintTargetMol`
pub fn pattern_fingerprint_target_mol(mol: &ROMol, fp_size: u32) -> Result<ExplicitBitVect> {
    pattern_fingerprint_mol(mol, fp_size, true)
}
