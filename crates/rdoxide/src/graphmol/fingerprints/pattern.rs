//! Pattern fingerprint (upstream `GraphMol/Fingerprints/PatternFingerprints.cpp`).
//!
//! A substructure-screening fingerprint: every match of a fixed set of small
//! SMARTS patterns sets one bit that encodes the pattern, the matched atomic
//! numbers and the matched bond types, plus one "count" bit per match.

use std::sync::OnceLock;

use super::hash_combine;
use crate::datastructs::{BitVect, ExplicitBitVect};
use crate::error::{RDKitError, Result};
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::find_rings::fast_find_rings;
use crate::graphmol::queries::is_complex_atom_query;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};

/// `PatternFingerprintMolVersion`
pub const PATTERN_FINGERPRINT_MOL_VERSION: &str = "1.0.0";

/// The SMARTS patterns of upstream's `pqs` table, in order.
const PATTERNS: &[&str] = &[
    "[*]~[*]",
    "[*]~[*]~[*]",
    "[R]~1~[R]~[R]~1",
    "[*]~[*](~[*])~[*]",
    "[R]~1[R]~[R]~[R]~1",
    "[*]~[*]~[*](~[*])~[*]",
    "[R]~1~[R]~[R]~[R]~[R]~1",
    "[R]~1~[R]~[R]~[R]~[R]~[R]~1",
    "[R](@[R])(@[R])~[R]~[R](@[R])(@[R])",
    "[R](@[R])(@[R])~[R]@[R]~[R](@[R])(@[R])",
    "[*]~[R](@[R])@[R](@[R])~[*]",
    "[*]~[R](@[R])@[R]@[R](@[R])~[*]",
    "[*]", // single atom fragment
];

/// Upstream raises `maxMatches` so high that it is effectively unlimited
/// (github #2614).
const MAX_MATCHES: u32 = 100_000_000;

fn patterns() -> Result<&'static [ROMol]> {
    static PATS: OnceLock<Vec<ROMol>> = OnceLock::new();
    if let Some(p) = PATS.get() {
        return Ok(p);
    }
    let params = SmartsParserParams::default();
    let parsed = PATTERNS
        .iter()
        .map(|sma| {
            mol_from_smarts(sma, &params)?.ok_or_else(|| RDKitError::invariant("bad smarts"))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(PATS.get_or_init(|| parsed))
}

/// An atom that carries a query feature the fingerprint cannot encode.
fn is_pattern_query_atom(a: &Atom) -> bool {
    // `isComplexQuery` no longer treats "AtomNull" as complex, but for the
    // pattern fingerprint it is definitely a query feature
    a.query()
        .is_some_and(|q| q.description() == "AtomNull" || is_complex_atom_query(Some(q)))
}

/// `isPatternComplexQuery(const Bond *)`
fn is_pattern_complex_query_bond(b: &Bond) -> bool {
    b.query()
        .is_some_and(|q| q.negate || q.description() != "BondOrder")
}

/// `isTautomerBondQuery`; only meaningful for bonds that are complex queries.
fn is_tautomer_bond_query(b: &Bond) -> bool {
    b.query().is_some_and(|q| {
        matches!(
            q.description(),
            "SingleOrDoubleOrAromaticBond" | "SingleOrAromaticBond"
        )
    })
}

/// The per-molecule query flags consulted for every match.
struct QueryFlags {
    atom: Vec<bool>,
    bond: Vec<bool>,
    tautomer_bond: Vec<bool>,
}

impl QueryFlags {
    fn new(mol: &ROMol, tautomeric: bool) -> Self {
        let atom = mol.atoms().iter().map(is_pattern_query_atom).collect();
        let mut bond = Vec::with_capacity(mol.num_bonds() as usize);
        let mut tautomer_bond = Vec::with_capacity(mol.num_bonds() as usize);
        for b in mol.bonds() {
            let complex = is_pattern_complex_query_bond(b);
            bond.push(complex);
            tautomer_bond.push(complex && tautomeric && is_tautomer_bond_query(b));
        }
        QueryFlags {
            atom,
            bond,
            tautomer_bond,
        }
    }
}

/// The bits set by one match of pattern number `pattern_idx`.
struct MatchBits {
    /// Always set: counts the occurrences of the pattern.
    count_bit: u32,
    /// The bit encoding atoms and bond types, unless the match is a query or
    /// has a tautomer query bond.
    bond_bit: Option<u32>,
    /// The tautomer-insensitive bit (tautomeric fingerprints only).
    tautomer_bit: Option<u32>,
}

fn bond_type_hash_value(b: &Bond) -> u32 {
    if b.is_aromatic() {
        BondType::Aromatic as u32
    } else {
        b.bond_type() as u32
    }
}

/// The bits for one match, or `None` when the match touches a query atom or a
/// query bond (those matches only contribute their count bit).
fn bits_for_match(
    mol: &ROMol,
    patt: &ROMol,
    mv: &[(i32, i32)],
    pattern_idx: u32,
    flags: &QueryFlags,
    tautomeric: bool,
    count_bit: u32,
) -> Result<MatchBits> {
    let mut res = MatchBits {
        count_bit,
        bond_bit: None,
        tautomer_bit: None,
    };
    let mut bit_id = pattern_idx;
    let mut amap = vec![0u32; mv.len()];
    for &(pi, mi) in mv {
        if flags.atom[mi as usize] {
            return Ok(res);
        }
        hash_combine(&mut bit_id, mol.atom(mi as u32).atomic_num());
        amap[pi as usize] = mi as u32;
    }
    let mut tautomer_bit_id = bit_id;
    let mut tautomer_query = false;
    for pbond in patt.bonds() {
        let mbond = mol
            .get_bond_between_atoms(
                amap[pbond.begin_atom_idx() as usize],
                amap[pbond.end_atom_idx() as usize],
            )?
            .ok_or_else(|| RDKitError::invariant("matched bond missing"))?;
        let bond_idx = mbond.idx() as usize;
        if flags.bond[bond_idx] {
            if !flags.tautomer_bond[bond_idx] {
                return Ok(res);
            }
            tautomer_query = true;
        }
        if tautomeric
            && (flags.tautomer_bond[bond_idx]
                || mbond.is_aromatic()
                || matches!(
                    mbond.bond_type(),
                    BondType::Single | BondType::Double | BondType::Aromatic
                ))
        {
            hash_combine(&mut tautomer_bit_id, u32::MAX); // upstream: hash of -1
        }
        if !tautomer_query {
            hash_combine(&mut bit_id, bond_type_hash_value(mbond));
        }
    }
    if !tautomer_query {
        res.bond_bit = Some(bit_id);
    }
    if tautomeric {
        res.tautomer_bit = Some(tautomer_bit_id);
    }
    Ok(res)
}

fn update_pattern_fingerprint(
    mol: &ROMol,
    fp: &mut ExplicitBitVect,
    fp_size: u32,
    tautomeric: bool,
) -> Result<()> {
    let flags = QueryFlags::new(mol, tautomeric);
    let params = SubstructMatchParameters {
        uniquify: false,
        max_matches: MAX_MATCHES,
        ..Default::default()
    };
    for (i, patt) in patterns()?.iter().enumerate() {
        let pattern_idx = i as u32 + 1;
        let mut count_hash = pattern_idx + patt.num_atoms() + patt.num_bonds();
        for mv in substruct_match(mol, patt, &params)? {
            hash_combine(&mut count_hash, 0xBEEF);
            let bits = bits_for_match(mol, patt, &mv, pattern_idx, &flags, tautomeric, count_hash)?;
            fp.set_bit(bits.count_bit % fp_size)?;
            for bit in [bits.bond_bit, bits.tautomer_bit].into_iter().flatten() {
                fp.set_bit(bit % fp_size)?;
            }
        }
    }
    Ok(())
}

/// `PatternFingerprintMol`
///
/// `atomCounts` and `setOnlyBits` are accepted by upstream's API but have no
/// effect on the result, so they are not part of this signature.
pub fn pattern_fingerprint_mol(
    mol: &ROMol,
    fp_size: u32,
    tautomeric_fingerprint: bool,
) -> Result<ExplicitBitVect> {
    if fp_size == 0 {
        return Err(RDKitError::precondition("fpSize==0"));
    }
    let mut res = ExplicitBitVect::new(fp_size, false);
    if mol.ring_info().is_find_fast_or_better() {
        update_pattern_fingerprint(mol, &mut res, fp_size, tautomeric_fingerprint)?;
    } else {
        // upstream perceives rings on the (mutable) ring info of a const mol
        let mut local = mol.clone();
        fast_find_rings(&mut local)?;
        update_pattern_fingerprint(&local, &mut res, fp_size, tautomeric_fingerprint)?;
    }
    Ok(res)
}
