//! Binary molecule pickles (`MolPickler`), byte-compatible with RDKit, and the
//! property serialization from `RDGeneral/StreamOps.h`.
//!
//! Ported from RDKit `Code/GraphMol/MolPickler.{h,cpp}` and
//! `Code/RDGeneral/StreamOps.h` (2026.09.1). Pickles written by RDKit can be
//! read here and vice versa. Not ported yet: the version-1 format and custom (ExplicitBitVect) property handlers.
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeSet;
use std::sync::atomic::{AtomicU32, Ordering};

use crate::geometry::Point3D;
use crate::graphmol::atom::{Atom, ChiralType, HybridizationType};
use crate::graphmol::bond::{Bond, BondDir, BondStereo, BondType};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::find_rings::find_ring_families;
use crate::graphmol::monomer_info::{AtomMonomerInfo, AtomMonomerType};
use crate::graphmol::queries::{
    AtomFn, AtomQuery, BondFn, BondQuery, PropValue, Query, QueryKind, QueryMol, QueryTarget,
};
use crate::graphmol::ring_info::FindRingType;
use crate::graphmol::romol::ROMol;
use crate::graphmol::stereo_group::{StereoGroup, StereoGroupType, assign_stereo_group_ids};
use crate::graphmol::substance_group::{SubstanceGroup, add_substance_group};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{Dict, RDProps, RDValue};
use crate::{RDKitError, Result};

/// `PicklerOps::PropertyPickleOptions`
pub mod flags {
    pub const NO_PROPS: u32 = 0;
    pub const MOL_PROPS: u32 = 0x1;
    pub const ATOM_PROPS: u32 = 0x2;
    pub const BOND_PROPS: u32 = 0x4;
    pub const QUERY_ATOM_DATA: u32 = 0x2;
    pub const PRIVATE_PROPS: u32 = 0x10;
    pub const COMPUTED_PROPS: u32 = 0x20;
    pub const ALL_PROPS: u32 = 0x0000_FFFF;
    pub const COORDS_AS_DOUBLE: u32 = 0x0001_0000;
    pub const NO_CONFORMERS: u32 = 0x0002_0000;
}
use flags::*;

pub const VERSION_MAJOR: i32 = 16;
pub const VERSION_MINOR: i32 = 4;
pub const VERSION_PATCH: i32 = 0;
pub const ENDIAN_ID: i32 = 0xDEAD_BEEFu32 as i32;

/// `MolPickler::Tags`
#[allow(dead_code)]
pub(crate) mod tags {
    pub const VERSION: u8 = 0;
    pub const BEGINATOM: u8 = 1;
    pub const ENDATOM: u8 = 10;
    pub const BEGINBOND: u8 = 11;
    pub const ENDBOND: u8 = 17;
    pub const BEGINPROPS: u8 = 18;
    pub const ENDPROPS: u8 = 19;
    pub const BEGINSSSR: u8 = 20;
    pub const ENDSSSR: u8 = 21;
    pub const ENDMOL: u8 = 22;
    pub const BEGINCONFS: u8 = 23;
    pub const ATOM_MAPNUMBER: u8 = 24;
    pub const BEGINQUERY: u8 = 25;
    pub const QUERY_VALUE: u8 = 26;
    pub const QUERY_ISNEGATED: u8 = 27;
    pub const QUERY_NUMCHILDREN: u8 = 28;
    pub const QUERY_BOOL: u8 = 29;
    pub const QUERY_AND: u8 = 30;
    pub const QUERY_OR: u8 = 31;
    pub const QUERY_XOR: u8 = 32;
    pub const QUERY_EQUALS: u8 = 33;
    pub const QUERY_GREATER: u8 = 34;
    pub const QUERY_GREATEREQUAL: u8 = 35;
    pub const QUERY_LESS: u8 = 36;
    pub const QUERY_LESSEQUAL: u8 = 37;
    pub const QUERY_RANGE: u8 = 38;
    pub const QUERY_SET: u8 = 39;
    pub const QUERY_NULL: u8 = 40;
    pub const QUERY_ATOMRING: u8 = 41;
    pub const QUERY_RECURSIVE: u8 = 42;
    pub const ENDQUERY: u8 = 43;
    pub const ATOM_DUMMYLABEL: u8 = 44;
    pub const BEGIN_PDB_RESIDUE: u8 = 45;
    pub const BEGINATOMPROPS: u8 = 58;
    pub const BEGINBONDPROPS: u8 = 59;
    pub const BEGINQUERYATOMDATA: u8 = 60;
    pub const BEGINSGROUP: u8 = 61;
    pub const BEGINSTEREOGROUP: u8 = 62;
    pub const BEGINCONFPROPS: u8 = 63;
    pub const BEGINCONFS_DOUBLE: u8 = 64;
    pub const QUERY_TYPELABEL: u8 = 65;
    pub const BEGINSYMMSSSR: u8 = 66;
    pub const BEGINFASTFIND: u8 = 67;
    pub const BEGINFINDOTHERORUNKNOWN: u8 = 68;
    pub const QUERY_PROPERTY: u8 = 69;
    pub const QUERY_PROPERTY_WITH_VALUE: u8 = 70;
    pub const ATOM_PDB_RESIDUE_SERIALNUMBER: u8 = 46;
    pub const ATOM_PDB_RESIDUE_ALTLOC: u8 = 47;
    pub const ATOM_PDB_RESIDUE_RESIDUENAME: u8 = 48;
    pub const ATOM_PDB_RESIDUE_CHAINID: u8 = 49;
    pub const ATOM_PDB_RESIDUE_INSERTIONCODE: u8 = 50;
    pub const ATOM_PDB_RESIDUE_OCCUPANCY: u8 = 51;
    pub const ATOM_PDB_RESIDUE_TEMPFACTOR: u8 = 52;
    pub const ATOM_PDB_RESIDUE_ISHETEROATOM: u8 = 53;
    pub const ATOM_PDB_RESIDUE_SECONDARYSTRUCTURE: u8 = 54;
    pub const ATOM_PDB_RESIDUE_RESIDUENUMBER: u8 = 55;
    pub const ATOM_PDB_RESIDUE_SEGMENTNUMBER: u8 = 56;
    pub const END_PDB_RESIDUE: u8 = 57;
    pub const ATOM_PDB_RESIDUE_MONOMERCLASS: u8 = 71;
    pub const BEGIN_ATOM_MONOMER_INFO: u8 = 72;
    pub const ATOM_MONOMER_INFO_RESIDUENAME: u8 = 73;
    pub const ATOM_MONOMER_INFO_RESIDUENUMBER: u8 = 74;
    pub const ATOM_MONOMER_INFO_CHAINID: u8 = 75;
    pub const ATOM_MONOMER_INFO_MONOMERCLASS: u8 = 76;
    pub const END_ATOM_MONOMER_INFO: u8 = 77;
    pub const INVALID_TAG: u8 = 255;
}
use tags::*;

/// `DTags` (StreamOps.h)
mod dtags {
    pub const STRING: u8 = 0;
    pub const INT: u8 = 1;
    pub const UNSIGNED_INT: u8 = 2;
    pub const BOOL: u8 = 3;
    pub const FLOAT: u8 = 4;
    pub const DOUBLE: u8 = 5;
    pub const VEC_STRING: u8 = 6;
    pub const VEC_INT: u8 = 7;
    pub const VEC_UINT: u8 = 8;
    pub const VEC_FLOAT: u8 = 10;
    pub const VEC_DOUBLE: u8 = 11;
    pub const CUSTOM: u8 = 0xFE;
}

static DEFAULT_PROPERTIES: AtomicU32 = AtomicU32::new(NO_PROPS);

/// `MolPickler::getDefaultPickleProperties`
pub fn get_default_pickle_properties() -> u32 {
    DEFAULT_PROPERTIES.load(Ordering::Relaxed)
}

/// `MolPickler::setDefaultPickleProperties`
pub fn set_default_pickle_properties(props: u32) {
    DEFAULT_PROPERTIES.store(props, Ordering::Relaxed);
}

/// `MolPicklerException` (a `RuntimeError` in Python)
fn pickle_error(msg: impl Into<String>) -> RDKitError {
    RDKitError::Runtime(msg.into())
}

// ================================================================== streams

/// Little-endian writer (`streamWrite`).
#[derive(Default)]
struct W(Vec<u8>);

impl W {
    fn u8(&mut self, v: u8) {
        self.0.push(v);
    }
    fn i8(&mut self, v: i8) {
        self.0.push(v as u8);
    }
    fn i16(&mut self, v: i16) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn u16(&mut self, v: u16) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn i32(&mut self, v: i32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn u32(&mut self, v: u32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn u64(&mut self, v: u64) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn f32(&mut self, v: f32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn f64(&mut self, v: f64) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    /// `streamWrite(ss, std::string)`: u32 length + bytes
    fn str(&mut self, s: &str) {
        self.u32(s.len() as u32);
        self.0.extend_from_slice(s.as_bytes());
    }
    /// an index-sized value (`T` = `unsigned char` or `int32_t`)
    fn t(&mut self, wide: bool, v: i64) {
        if wide {
            self.i32(v as i32);
        } else {
            self.u8(v as u8);
        }
    }
    /// `write_sstream_to_stream`: i32 size + block
    fn block(&mut self, b: &W) {
        self.i32(b.0.len() as i32);
        self.0.extend_from_slice(&b.0);
    }
}

/// Little-endian reader (`streamRead`).
struct R<'a> {
    data: &'a [u8],
    pos: usize,
}

fn eof() -> RDKitError {
    pickle_error("Bad pickle format: unexpected End-of-File while reading")
}

impl R<'_> {
    fn bytes(&mut self, n: usize) -> Result<&[u8]> {
        if self.pos + n > self.data.len() {
            return Err(eof());
        }
        let b = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(b)
    }
    fn arr<const N: usize>(&mut self) -> Result<[u8; N]> {
        let mut a = [0u8; N];
        a.copy_from_slice(self.bytes(N)?);
        Ok(a)
    }
    fn u8(&mut self) -> Result<u8> {
        Ok(self.bytes(1)?[0])
    }
    fn i8(&mut self) -> Result<i8> {
        Ok(self.u8()? as i8)
    }
    /// `streamReadPositiveChar`
    fn pos_char(&mut self) -> Result<u8> {
        let v = self.i8()?;
        if v < 0 {
            return Err(pickle_error("invalid value in pickle"));
        }
        Ok(v as u8)
    }
    fn i16(&mut self) -> Result<i16> {
        Ok(i16::from_le_bytes(self.arr()?))
    }
    fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.arr()?))
    }
    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.arr()?))
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.arr()?))
    }
    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.arr()?))
    }
    fn f32(&mut self) -> Result<f32> {
        Ok(f32::from_le_bytes(self.arr()?))
    }
    fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(self.arr()?))
    }
    fn str(&mut self) -> Result<String> {
        let l = self.u32()? as usize;
        Ok(String::from_utf8_lossy(self.bytes(l)?).into_owned())
    }
    fn t(&mut self, wide: bool) -> Result<i64> {
        Ok(if wide {
            self.i32()? as i64
        } else {
            self.u8()? as i64
        })
    }
    /// `streamRead(ss, Tags &, version)`
    fn tag(&mut self, version: i32) -> Result<u8> {
        let v = if version < 7000 {
            let t = self.i32()?;
            if t < 0 || t >= INVALID_TAG as i32 {
                return Err(pickle_error("Invalid tag found."));
            }
            t as u8
        } else {
            self.u8()?
        };
        if v == INVALID_TAG {
            return Err(pickle_error("Invalid tag found."));
        }
        Ok(v)
    }
    fn skip(&mut self, n: i32) -> Result<()> {
        self.bytes(n.max(0) as usize).map(|_| ())
    }
}

// ================================================================== properties

fn is_serializable(v: &RDValue) -> bool {
    !matches!(
        v,
        RDValue::Int64(_) | RDValue::UInt64(_) | RDValue::BitVect(_)
    )
}

/// `streamWriteProp`
fn write_prop(w: &mut W, key: &str, v: &RDValue) -> bool {
    if !is_serializable(v) {
        return false;
    }
    w.str(key);
    match v {
        RDValue::String(s) => {
            w.u8(dtags::STRING);
            w.str(s);
        }
        RDValue::Int(i) => {
            w.u8(dtags::INT);
            w.i32(*i);
        }
        RDValue::UnsignedInt(u) => {
            w.u8(dtags::UNSIGNED_INT);
            w.u32(*u);
        }
        RDValue::Bool(b) => {
            w.u8(dtags::BOOL);
            w.u8(*b as u8);
        }
        RDValue::Float(f) => {
            w.u8(dtags::FLOAT);
            w.f32(*f);
        }
        RDValue::Double(d) => {
            w.u8(dtags::DOUBLE);
            w.f64(*d);
        }
        RDValue::VecString(v) => {
            w.u8(dtags::VEC_STRING);
            w.u64(v.len() as u64);
            v.iter().for_each(|s| w.str(s));
        }
        RDValue::VecDouble(v) => {
            w.u8(dtags::VEC_DOUBLE);
            w.u64(v.len() as u64);
            v.iter().for_each(|&x| w.f64(x));
        }
        RDValue::VecFloat(v) => {
            w.u8(dtags::VEC_FLOAT);
            w.u64(v.len() as u64);
            v.iter().for_each(|&x| w.f32(x));
        }
        RDValue::VecInt(v) => {
            w.u8(dtags::VEC_INT);
            w.u64(v.len() as u64);
            v.iter().for_each(|&x| w.i32(x));
        }
        RDValue::VecUnsignedInt(v) => {
            w.u8(dtags::VEC_UINT);
            w.u64(v.len() as u64);
            v.iter().for_each(|&x| w.u32(x));
        }
        RDValue::Int64(_) | RDValue::UInt64(_) | RDValue::BitVect(_) => unreachable!(),
    }
    true
}

/// `streamWriteProps<COUNT>(ss, props, savePrivate, saveComputed, handlers, ignore)`
fn write_props(
    w: &mut W,
    props: &impl RDProps,
    private: bool,
    computed: bool,
    ignore: &[&str],
    wide_count: bool,
) -> bool {
    let names: Vec<String> = props
        .get_prop_list(private, computed)
        .into_iter()
        .filter(|n| !ignore.contains(&n.as_str()))
        .collect();
    let dict: &Dict = props.dict();
    let items: Vec<(&str, &RDValue)> = dict
        .iter()
        .filter(|(k, v)| names.iter().any(|n| n == k) && is_serializable(v))
        .collect();
    if wide_count {
        w.u32(items.len() as u32);
    } else {
        w.u16(items.len() as u16);
    }
    if items.is_empty() {
        return false;
    }
    for (k, v) in items {
        write_prop(w, k, v);
    }
    true
}

/// `streamReadProp`
fn read_prop(r: &mut R<'_>) -> Result<Option<(String, RDValue)>> {
    let key = r.str()?;
    let t = r.u8()?;
    let vec_len = |r: &mut R<'_>| -> Result<usize> { Ok(r.u64()? as usize) };
    let v = match t {
        dtags::INT => RDValue::Int(r.i32()?),
        dtags::UNSIGNED_INT => RDValue::UnsignedInt(r.u32()?),
        dtags::BOOL => RDValue::Bool(r.u8()? != 0),
        dtags::FLOAT => RDValue::Float(r.f32()?),
        dtags::DOUBLE => RDValue::Double(r.f64()?),
        dtags::STRING => RDValue::String(r.str()?),
        dtags::VEC_STRING => {
            let n = vec_len(r)?;
            RDValue::VecString((0..n).map(|_| r.str()).collect::<Result<_>>()?)
        }
        dtags::VEC_INT => {
            let n = vec_len(r)?;
            RDValue::VecInt((0..n).map(|_| r.i32()).collect::<Result<_>>()?)
        }
        dtags::VEC_UINT => {
            let n = vec_len(r)?;
            RDValue::VecUnsignedInt((0..n).map(|_| r.u32()).collect::<Result<_>>()?)
        }
        dtags::VEC_FLOAT => {
            let n = vec_len(r)?;
            RDValue::VecFloat((0..n).map(|_| r.f32()).collect::<Result<_>>()?)
        }
        dtags::VEC_DOUBLE => {
            let n = vec_len(r)?;
            RDValue::VecDouble((0..n).map(|_| r.f64()).collect::<Result<_>>()?)
        }
        dtags::CUSTOM => {
            // TODO(DataStructs): the ExplicitBitVect property handler
            let _handler = r.str()?;
            return Ok(None);
        }
        _ => return Ok(None),
    };
    Ok(Some((key, v)))
}

/// `streamReadProps<COUNT>(ss, props, handlers, reset)`
fn read_props(r: &mut R<'_>, dict: &mut Dict, wide_count: bool, reset: bool) -> Result<()> {
    let count = if wide_count {
        r.u32()? as usize
    } else {
        r.u16()? as usize
    };
    if reset {
        dict.reset();
    }
    let mut pairs = Vec::with_capacity(count);
    for _ in 0..count {
        match read_prop(r)? {
            Some(p) => pairs.push(p),
            None => {
                return Err(RDKitError::invariant(
                    "Corrupted property serialization detected",
                ));
            }
        }
    }
    dict.extend_pairs(pairs);
    Ok(())
}

const EXPLICIT_ATOM_PROPS: [(&str, u8); 4] = [
    (cp::molStereoCare, 0x1),
    (cp::molParity, 0x2),
    (cp::molInversionFlag, 0x4),
    (cp::_ChiralityPossible, 0x8),
];
const EXPLICIT_BOND_PROPS: [(&str, u8); 5] = [
    (cp::_MolFileBondType, 0x1),
    (cp::_MolFileBondStereo, 0x2),
    (cp::_MolFileBondCfg, 0x4),
    (cp::_MolFileBondQuery, 0x8),
    (cp::molStereoCare, 0x10),
];

/// `pickleExplicitProperties<SAVEAS>`: `wide` = int16 (atoms) or int8 (bonds)
fn write_explicit_props(w: &mut W, props: &impl RDProps, table: &[(&str, u8)], wide: bool) -> bool {
    let mut bprops = 0u8;
    let mut vals = Vec::new();
    for &(name, bit) in table {
        if let Some(v) = props.get_prop(name).and_then(RDValue::as_int) {
            bprops |= bit;
            vals.push(v);
        }
    }
    w.u8(bprops);
    for v in &vals {
        if wide {
            w.i16(*v as i16);
        } else {
            w.i8(*v as i8);
        }
    }
    !vals.is_empty()
}

fn read_explicit_props(
    r: &mut R<'_>,
    dict: &mut Dict,
    table: &[(&str, u8)],
    wide: bool,
    version: i32,
) -> Result<()> {
    if version >= 14000 {
        let bprops = r.u8()?;
        for &(name, bit) in table {
            if bprops & bit != 0 {
                let v = if wide {
                    r.i16()? as i32
                } else {
                    r.i8()? as i32
                };
                dict.set_val(name, v);
            }
        }
    }
    Ok(())
}

fn atom_props_to_ignore() -> Vec<&'static str> {
    let mut v: Vec<&str> = EXPLICIT_ATOM_PROPS.iter().map(|p| p.0).collect();
    v.push(cp::molAtomMapNumber);
    v.push(cp::dummyLabel);
    v
}

// ================================================================== queries

/// The serializable part of a query node (upstream
/// `PicklerOps::QueryDetails`); shared by the binary pickler and the JSON
/// writer.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum QueryDetails {
    /// AND/OR/XOR/NULL
    Tag(u8),
    /// the equality-like comparisons and the atom-ring query
    ValTol {
        tag: u8,
        val: i32,
        tol: i32,
    },
    Range {
        lower: i32,
        upper: i32,
        tol: i32,
        ends: u8,
    },
    Set(BTreeSet<i32>),
    Property(String),
    PropertyWithValue {
        prop: String,
        val: PropValue,
        tol: f64,
    },
}

/// `PicklerOps::getQueryDetails`; `None` for recursive queries, which both
/// formats handle themselves.
pub(crate) fn query_details<F: QueryTarget>(q: &Query<F>) -> Option<QueryDetails> {
    let val_tol = |tag: u8| QueryDetails::ValTol {
        tag,
        val: q.val,
        tol: q.tol,
    };
    Some(match &q.kind {
        QueryKind::Recursive { .. } => return None,
        QueryKind::And => QueryDetails::Tag(QUERY_AND),
        QueryKind::Or => QueryDetails::Tag(QUERY_OR),
        QueryKind::XOr => QueryDetails::Tag(QUERY_XOR),
        QueryKind::Base { .. } => QueryDetails::Tag(QUERY_NULL),
        QueryKind::Equality => val_tol(QUERY_EQUALS),
        QueryKind::Greater => val_tol(QUERY_GREATER),
        QueryKind::GreaterEqual => val_tol(QUERY_GREATEREQUAL),
        QueryKind::Less => val_tol(QUERY_LESS),
        QueryKind::LessEqual => val_tol(QUERY_LESSEQUAL),
        QueryKind::AtomRing => val_tol(QUERY_ATOMRING),
        QueryKind::HasProp(name) => QueryDetails::Property(name.clone()),
        QueryKind::Range {
            lower,
            upper,
            lower_open,
            upper_open,
        } => QueryDetails::Range {
            lower: *lower,
            upper: *upper,
            tol: q.tol,
            ends: ((*lower_open as u8) << 1) | *upper_open as u8,
        },
        QueryKind::Set(set) => QueryDetails::Set(set.clone()),
        QueryKind::HasPropWithValue { prop, val, tol } => QueryDetails::PropertyWithValue {
            prop: prop.clone(),
            val: val.clone(),
            // string comparisons have no tolerance
            tol: if matches!(val, PropValue::String(_)) {
                0.0
            } else {
                *tol
            },
        },
    })
}

/// `pickleQuery`
fn pickle_query<F: QueryTarget>(w: &mut W, q: &Query<F>) -> Result<()> {
    w.str(&q.description);
    if !q.type_label.is_empty() {
        w.u8(QUERY_TYPELABEL);
        w.str(&q.type_label);
    }
    if q.negate {
        w.u8(QUERY_ISNEGATED);
    }
    match (&q.kind, query_details(q)) {
        (QueryKind::Recursive { mol, .. }, _) => {
            w.u8(QUERY_RECURSIVE);
            w.u8(QUERY_VALUE);
            let p = pickle_mol(&mol.0, get_default_pickle_properties())?;
            w.0.extend_from_slice(&p);
        }
        (_, Some(QueryDetails::Tag(tag))) => w.u8(tag),
        (_, Some(QueryDetails::ValTol { tag, val, tol })) => {
            w.u8(tag);
            w.u8(QUERY_VALUE);
            w.i32(val);
            w.i32(tol);
        }
        (
            _,
            Some(QueryDetails::Range {
                lower,
                upper,
                tol,
                ends,
            }),
        ) => {
            w.u8(QUERY_RANGE);
            w.u8(QUERY_VALUE);
            w.i32(lower);
            w.i32(upper);
            w.i32(tol);
            w.u8(ends);
        }
        (_, Some(QueryDetails::Set(set))) => {
            w.u8(QUERY_SET);
            w.u8(QUERY_VALUE);
            w.i32(set.len() as i32);
            set.iter().for_each(|&v| w.i32(v));
        }
        (_, Some(QueryDetails::Property(name))) => {
            w.u8(QUERY_PROPERTY);
            w.u8(QUERY_VALUE);
            w.str(&name);
        }
        (_, Some(QueryDetails::PropertyWithValue { prop, val, tol })) => {
            w.u8(QUERY_PROPERTY_WITH_VALUE);
            w.u8(QUERY_VALUE);
            w.f64(tol);
            let v = match val {
                PropValue::Int(i) => RDValue::Int(i),
                PropValue::UInt(u) => RDValue::UnsignedInt(u),
                PropValue::Bool(b) => RDValue::Bool(b),
                PropValue::Float(f) => RDValue::Float(f),
                PropValue::Double(d) => RDValue::Double(d),
                PropValue::String(s) => RDValue::String(s),
                PropValue::BitVect(_) => {
                    return Err(pickle_error(
                        "bit vector property queries cannot be pickled",
                    ));
                }
            };
            write_prop(w, &prop, &v);
        }
        (_, None) => unreachable!("only recursive queries have no details"),
    }
    w.u8(QUERY_NUMCHILDREN);
    w.u8(q.children.len() as u8);
    for c in &q.children {
        pickle_query(w, c)?;
    }
    Ok(())
}

/// `buildBaseQuery`
fn build_base_query<F: QueryTarget>(r: &mut R<'_>, tag: u8, version: i32) -> Result<Query<F>> {
    let expect_value = |r: &mut R<'_>| -> Result<()> {
        if r.tag(version)? != QUERY_VALUE {
            return Err(pickle_error(
                "Bad pickle format: QUERY_VALUE tag not found.",
            ));
        }
        Ok(())
    };
    let val_tol = |kind: QueryKind, r: &mut R<'_>| -> Result<Query<F>> {
        let mut q = Query::new(kind);
        expect_value(r)?;
        q.val = r.i32()?;
        q.tol = r.i32()?;
        Ok(q)
    };
    Ok(match tag {
        QUERY_AND => Query::new(QueryKind::And),
        QUERY_OR => Query::new(QueryKind::Or),
        QUERY_XOR => Query::new(QueryKind::XOr),
        QUERY_EQUALS => val_tol(QueryKind::Equality, r)?,
        QUERY_GREATER => val_tol(QueryKind::Greater, r)?,
        QUERY_GREATEREQUAL => val_tol(QueryKind::GreaterEqual, r)?,
        QUERY_LESS => val_tol(QueryKind::Less, r)?,
        QUERY_LESSEQUAL => val_tol(QueryKind::LessEqual, r)?,
        QUERY_RANGE => {
            expect_value(r)?;
            let lower = r.i32()?;
            let upper = r.i32()?;
            let tol = r.i32()?;
            let ends = r.i8()? as u8;
            let mut q = Query::new(QueryKind::Range {
                lower,
                upper,
                lower_open: ends & 0x2 != 0,
                upper_open: ends & 0x1 != 0,
            });
            q.tol = tol;
            q
        }
        QUERY_SET => {
            expect_value(r)?;
            let n = r.i32()?;
            let mut set = BTreeSet::new();
            for _ in 0..n.max(0) {
                set.insert(r.i32()?);
            }
            Query::new(QueryKind::Set(set))
        }
        QUERY_NULL => Query::new(QueryKind::Base { always_true: false }),
        QUERY_PROPERTY => {
            expect_value(r)?;
            Query::new(QueryKind::HasProp(r.str()?))
        }
        QUERY_PROPERTY_WITH_VALUE => {
            expect_value(r)?;
            let tol = r.f64()?;
            let Some((key, v)) = read_prop(r)? else {
                return Err(pickle_error("unknown query-type tag encountered"));
            };
            let val = match v {
                RDValue::Int(i) => PropValue::Int(i),
                RDValue::UnsignedInt(u) => PropValue::UInt(u),
                RDValue::Bool(b) => PropValue::Bool(b),
                RDValue::Float(f) => PropValue::Float(f),
                RDValue::Double(d) => PropValue::Double(d),
                RDValue::String(s) => PropValue::String(s),
                _ => return Err(pickle_error("unknown query-type tag encountered")),
            };
            Query::new(QueryKind::HasPropWithValue {
                prop: key,
                val,
                tol,
            })
        }
        _ => return Err(pickle_error("unknown query-type tag encountered")),
    })
}

/// `QueryOps::finalizeQueryFromDescription` (atoms)
pub(crate) fn finalize_atom_query(q: &mut AtomQuery) -> Result<()> {
    let full = q.description.clone();
    let (descr, qtype) = if let Some(d) = full.strip_prefix("range_") {
        (d, 'r')
    } else if let Some(d) = full.strip_prefix("less_") {
        (d, 'l')
    } else if let Some(d) = full.strip_prefix("greater_") {
        (d, 'g')
    } else {
        (full.as_str(), 'e')
    };
    let f = match descr {
        "AtomRingBondCount" => AtomFn::RingBondCount,
        "AtomHasRingBond" => AtomFn::HasRingBond,
        "AtomRingSize" => match qtype {
            'e' => AtomFn::RingSize(q.val),
            'r' => match q.kind {
                QueryKind::Range {
                    lower,
                    upper,
                    lower_open,
                    upper_open,
                } => AtomFn::RingSizeRange {
                    lower,
                    upper,
                    lower_open,
                    upper_open,
                },
                _ => return Err(pickle_error("bad range query type")),
            },
            'l' => AtomFn::RingSizeRange {
                lower: q.val,
                upper: -1,
                lower_open: false,
                upper_open: false,
            },
            _ => AtomFn::RingSizeRange {
                lower: -1,
                upper: q.val,
                lower_open: false,
                upper_open: false,
            },
        },
        "AtomMinRingSize" => AtomFn::MinRingSize,
        "AtomImplicitValence" => AtomFn::ImplicitValence,
        "AtomTotalValence" => AtomFn::TotalValence,
        "AtomAtomicNum" => AtomFn::Num,
        "AtomExplicitDegree" => AtomFn::ExplicitDegree,
        "AtomTotalDegree" => AtomFn::TotalDegree,
        "AtomHeavyAtomDegree" => AtomFn::HeavyAtomDegree,
        "AtomHCount" => AtomFn::HCount,
        "AtomImplicitHCount" => AtomFn::ImplicitHCount,
        "AtomHasImplicitH" => AtomFn::HasImplicitH,
        "AtomIsAromatic" => AtomFn::Aromatic,
        "AtomIsAliphatic" => AtomFn::Aliphatic,
        "AtomUnsaturated" => AtomFn::Unsaturated,
        "AtomMass" => AtomFn::Mass,
        "AtomIsotope" => AtomFn::Isotope,
        "AtomFormalCharge" => AtomFn::FormalCharge,
        "AtomNegativeFormalCharge" => AtomFn::NegativeFormalCharge,
        "AtomHybridization" => AtomFn::Hybridization,
        "AtomInRing" => AtomFn::InRing,
        "AtomInNRings" => AtomFn::InNRings,
        "AtomHasHeteroatomNeighbors" => AtomFn::HasHeteroatomNbrs,
        "AtomNumHeteroatomNeighbors" => AtomFn::NumHeteroatomNbrs,
        "AtomNonHydrogenDegree" => AtomFn::NonHydrogenDegree,
        "AtomHasAliphaticHeteroatomNeighbors" => AtomFn::HasAliphaticHeteroatomNbrs,
        "AtomNumAliphaticHeteroatomNeighbors" => AtomFn::NumAliphaticHeteroatomNbrs,
        "AtomNull" => {
            q.kind = QueryKind::Base { always_true: true };
            AtomFn::Null
        }
        "AtomType" => AtomFn::Type,
        "AtomNumRadicalElectrons" => AtomFn::NumRadicalElectrons,
        "RecursiveStructure" => AtomFn::AtIdx,
        "AtomAnd" | "AtomOr" | "AtomXor" | "HasProp" | "HasPropWithValue" => return Ok(()),
        _ => {
            return Err(RDKitError::Value(format!(
                "Do not know how to finalize query: '{descr}'"
            )));
        }
    };
    q.data_func = Some(f);
    Ok(())
}

/// `QueryOps::finalizeQueryFromDescription` (bonds)
pub(crate) fn finalize_bond_query(q: &mut BondQuery) -> Result<()> {
    let f = match q.description.as_str() {
        "BondRingSize" => {
            if !(3..=20).contains(&q.val) {
                return Err(RDKitError::range("tgt", q.val, 20));
            }
            BondFn::RingSize(q.val)
        }
        "BondMinRingSize" => BondFn::MinRingSize,
        "BondOrder" => BondFn::Order,
        "BondDir" => BondFn::Dir,
        "BondInRing" => BondFn::InRing,
        "BondInNRings" => BondFn::InNRings,
        "SingleOrAromaticBond" => BondFn::IsSingleOrAromatic,
        "SingleOrDoubleBond" => BondFn::IsSingleOrDouble,
        "DoubleOrAromaticBond" => BondFn::IsDoubleOrAromatic,
        "SingleOrDoubleOrAromaticBond" => BondFn::IsSingleOrDoubleOrAromatic,
        "BondNull" => {
            q.kind = QueryKind::Base { always_true: true };
            BondFn::Null
        }
        "BondAnd" | "BondOr" | "BondXor" | "HasProp" | "HasPropWithValue" => return Ok(()),
        d => {
            return Err(RDKitError::Value(format!(
                "Do not know how to finalize query: '{d}'"
            )));
        }
    };
    q.data_func = Some(f);
    Ok(())
}

/// `unpickleQuery(ss, Atom const *, version)`
fn unpickle_atom_query(r: &mut R<'_>, version: i32) -> Result<AtomQuery> {
    let descr = r.str()?;
    let mut tag = r.tag(version)?;
    let mut type_label = String::new();
    if tag == QUERY_TYPELABEL {
        type_label = r.str()?;
        tag = r.tag(version)?;
    }
    let mut negated = false;
    if tag == QUERY_ISNEGATED {
        negated = true;
        tag = r.tag(version)?;
    }
    let mut q: AtomQuery = match tag {
        QUERY_ATOMRING => {
            if r.tag(version)? != QUERY_VALUE {
                return Err(pickle_error(
                    "Bad pickle format: QUERY_VALUE tag not found.",
                ));
            }
            let mut q = Query::new(QueryKind::AtomRing);
            q.val = r.i32()?;
            q.tol = r.i32()?;
            q.data_func = Some(AtomFn::RingMembership);
            q
        }
        QUERY_RECURSIVE => {
            if r.tag(version)? != QUERY_VALUE {
                return Err(pickle_error(
                    "Bad pickle format: QUERY_VALUE tag not found.",
                ));
            }
            let sub = mol_from_pickle_reader(r, ALL_PROPS)?;
            Query::new(QueryKind::Recursive {
                mol: QueryMol(std::sync::Arc::new(sub)),
                serial_number: 0,
                set: BTreeSet::new(),
                initialized: false,
            })
        }
        t => build_base_query(r, t, version)?,
    };
    q.negate = negated;
    q.description = descr;
    if !type_label.is_empty() {
        q.type_label = type_label;
    }
    finalize_atom_query(&mut q)?;
    if r.tag(version)? != QUERY_NUMCHILDREN {
        return Err(pickle_error(
            "Bad pickle format: QUERY_NUMCHILDREN tag not found.",
        ));
    }
    let n = r.u8()?;
    for _ in 0..n {
        let c = unpickle_atom_query(r, version)?;
        q.children.push(c);
    }
    Ok(q)
}

/// `unpickleQuery(ss, Bond const *, version)` (the type label is read but,
/// as upstream, not applied)
fn unpickle_bond_query(r: &mut R<'_>, version: i32) -> Result<BondQuery> {
    let descr = r.str()?;
    let mut tag = r.tag(version)?;
    if tag == QUERY_TYPELABEL {
        let _ = r.str()?;
        tag = r.tag(version)?;
    }
    let mut negated = false;
    if tag == QUERY_ISNEGATED {
        negated = true;
        tag = r.tag(version)?;
    }
    let mut q: BondQuery = build_base_query(r, tag, version)?;
    q.negate = negated;
    q.description = descr;
    finalize_bond_query(&mut q)?;
    if r.tag(version)? != QUERY_NUMCHILDREN {
        return Err(pickle_error(
            "Bad pickle format: QUERY_NUMCHILDREN tag not found.",
        ));
    }
    let n = r.u8()?;
    for _ in 0..n {
        let c = unpickle_bond_query(r, version)?;
        q.children.push(c);
    }
    Ok(q)
}

// ================================================================== atoms / bonds

/// `getAtomMapNumber` (int, or a string that parses as one)
fn atom_map_number(a: &Atom) -> Option<i32> {
    match a.get_prop(cp::molAtomMapNumber)? {
        RDValue::String(s) => s.trim().parse().ok(),
        v => v.as_int(),
    }
}

/// `_pickleAtomData`
fn pickle_atom_data(w: &mut W, a: &Atom) -> i32 {
    let mut flags = 0i32;
    if a.formal_charge() != 0 {
        flags |= 1 << 1;
        w.i8(a.formal_charge() as i8);
    }
    if a.chiral_tag() as u8 != 0 {
        flags |= 1 << 2;
        w.u8(a.chiral_tag() as u8);
    }
    if a.hybridization() != HybridizationType::SP3 {
        flags |= 1 << 3;
        w.u8(a.hybridization() as u8);
    }
    if a.num_explicit_hs() != 0 {
        flags |= 1 << 4;
        w.u8(a.num_explicit_hs() as u8);
    }
    if a.explicit_valence > 0 {
        flags |= 1 << 5;
        w.i8(a.explicit_valence);
    }
    if a.implicit_valence > 0 {
        flags |= 1 << 6;
        w.i8(a.implicit_valence);
    }
    if a.num_radical_electrons() != 0 {
        flags |= 1 << 7;
        w.u8(a.num_radical_electrons() as u8);
    }
    if a.isotope() > 0 {
        flags |= 1 << 8;
        w.u32(a.isotope());
    }
    flags
}

/// `_unpickleAtomData`
fn unpickle_atom_data(r: &mut R<'_>, a: &mut Atom) -> Result<()> {
    let flags = r.i32()?;
    if flags & 1 != 0 {
        let tmp = r.f32()?;
        let iso = (tmp as f64 + a.mass() + 0.0001).floor() as i32;
        a.set_isotope(iso.max(0) as u32);
    }
    a.set_formal_charge(if flags & (1 << 1) != 0 {
        r.i8()? as i32
    } else {
        0
    });
    let chiral = if flags & (1 << 2) != 0 {
        r.pos_char()?
    } else {
        0
    };
    a.set_chiral_tag(
        ChiralType::from_u8(chiral).ok_or_else(|| pickle_error("invalid value in pickle"))?,
    );
    let hyb = if flags & (1 << 3) != 0 {
        r.pos_char()?
    } else {
        HybridizationType::SP3 as u8
    };
    a.set_hybridization(
        HybridizationType::from_u8(hyb).ok_or_else(|| pickle_error("invalid value in pickle"))?,
    );
    a.set_num_explicit_hs(if flags & (1 << 4) != 0 {
        r.i8()? as u8 as u32
    } else {
        0
    });
    a.explicit_valence = if flags & (1 << 5) != 0 { r.i8()? } else { 0 };
    a.implicit_valence = if flags & (1 << 6) != 0 { r.i8()? } else { 0 };
    a.num_radical_electrons = if flags & (1 << 7) != 0 {
        r.pos_char()?
    } else {
        0
    };
    a.isotope = 0;
    if flags & (1 << 8) != 0 {
        let iso = r.u32()?;
        a.set_isotope(iso);
    }
    Ok(())
}

/// `_pickleAtom`
fn pickle_atom(w: &mut W, a: &Atom) -> Result<()> {
    w.u8((a.atomic_num() % 256) as u8);
    let mut flags = 0u8;
    if a.is_aromatic() {
        flags |= 1 << 6;
    }
    if a.no_implicit() {
        flags |= 1 << 5;
    }
    if a.has_query() {
        flags |= 1 << 4;
    }
    let map = atom_map_number(a);
    if map.is_some() {
        flags |= 1 << 3;
    }
    if a.has_prop(cp::dummyLabel) {
        flags |= 1 << 2;
    }
    if a.monomer_info().is_some() {
        flags |= 1 << 1;
    }
    w.u8(flags);
    let mut t = W::default();
    let pf = pickle_atom_data(&mut t, a);
    w.i32(pf);
    w.0.extend_from_slice(&t.0);
    if let Some(q) = a.query() {
        w.u8(BEGINQUERY);
        pickle_query(w, q)?;
        w.u8(ENDQUERY);
    }
    if let Some(m) = map {
        w.u8(ATOM_MAPNUMBER);
        if (0..128).contains(&m) {
            w.u8(m as u8);
        } else {
            w.u8(255);
            w.i32(m);
        }
    }
    if let Some(v) = a.get_prop(cp::dummyLabel) {
        w.u8(ATOM_DUMMYLABEL);
        w.str(&v.to_rd_string());
    }
    if let Some(info) = a.monomer_info() {
        pickle_monomer_info(w, info);
    }
    Ok(())
}

/// `pickleAtomPDBResidueInfo` / `pickleAtomMonomerInfo`, including the
/// surrounding begin/end tags.
fn pickle_monomer_info(w: &mut W, info: &AtomMonomerInfo) {
    let tagged_str = |w: &mut W, tag: u8, s: &str| {
        if !s.is_empty() {
            w.u8(tag);
            w.str(s);
        }
    };
    match &info.pdb {
        Some(p) => {
            w.u8(BEGIN_PDB_RESIDUE);
            w.str(&info.name);
            w.u32(info.monomer_type as u32);
            if p.serial_number != 0 {
                w.u8(ATOM_PDB_RESIDUE_SERIALNUMBER);
                w.i32(p.serial_number as i32);
            }
            tagged_str(w, ATOM_PDB_RESIDUE_ALTLOC, &p.alt_loc);
            tagged_str(w, ATOM_PDB_RESIDUE_RESIDUENAME, &info.residue_name);
            if info.residue_number != 0 {
                w.u8(ATOM_PDB_RESIDUE_RESIDUENUMBER);
                w.i32(info.residue_number);
            }
            tagged_str(w, ATOM_PDB_RESIDUE_CHAINID, &info.chain_id);
            tagged_str(w, ATOM_PDB_RESIDUE_INSERTIONCODE, &p.insertion_code);
            if p.occupancy != 0.0 {
                w.u8(ATOM_PDB_RESIDUE_OCCUPANCY);
                w.f64(p.occupancy);
            }
            if p.temp_factor != 0.0 {
                w.u8(ATOM_PDB_RESIDUE_TEMPFACTOR);
                w.f64(p.temp_factor);
            }
            if p.is_hetero_atom {
                w.u8(ATOM_PDB_RESIDUE_ISHETEROATOM);
                w.u8(1);
            }
            if p.secondary_structure != 0 {
                w.u8(ATOM_PDB_RESIDUE_SECONDARYSTRUCTURE);
                w.u32(p.secondary_structure);
            }
            if p.segment_number != 0 {
                w.u8(ATOM_PDB_RESIDUE_SEGMENTNUMBER);
                w.u32(p.segment_number);
            }
            tagged_str(w, ATOM_PDB_RESIDUE_MONOMERCLASS, &info.monomer_class);
            w.u8(END_PDB_RESIDUE);
        }
        None => {
            w.u8(BEGIN_ATOM_MONOMER_INFO);
            w.str(&info.name);
            w.u8(info.monomer_type as u8);
            tagged_str(w, ATOM_MONOMER_INFO_RESIDUENAME, &info.residue_name);
            if info.residue_number != 0 {
                w.u8(ATOM_MONOMER_INFO_RESIDUENUMBER);
                w.i32(info.residue_number);
            }
            tagged_str(w, ATOM_MONOMER_INFO_CHAINID, &info.chain_id);
            tagged_str(w, ATOM_MONOMER_INFO_MONOMERCLASS, &info.monomer_class);
            w.u8(END_ATOM_MONOMER_INFO);
        }
    }
}

/// `unpickleAtomPDBResidueInfo` (the `BEGIN_PDB_RESIDUE` tag is already read).
fn unpickle_pdb_residue_info(r: &mut R<'_>, version: i32) -> Result<AtomMonomerInfo> {
    let name = r.str()?;
    let typ = r.u32()?;
    let monomer_type = AtomMonomerType::from_u8(typ as u8).unwrap_or_default();
    if monomer_type != AtomMonomerType::PdbResidue {
        // upstream returns the bare base class and leaves the stream at the
        // PDB end tag, which the caller does not consume either
        return Ok(AtomMonomerInfo::new(monomer_type, &name));
    }
    let mut info = AtomMonomerInfo::new_pdb(&name);
    loop {
        let tag = r.tag(version)?;
        match tag {
            ATOM_PDB_RESIDUE_SERIALNUMBER => info.pdb_mut().serial_number = r.i32()? as u32,
            ATOM_PDB_RESIDUE_ALTLOC => info.pdb_mut().alt_loc = r.str()?,
            ATOM_PDB_RESIDUE_RESIDUENAME => info.residue_name = r.str()?,
            ATOM_PDB_RESIDUE_RESIDUENUMBER => info.residue_number = r.i32()?,
            ATOM_PDB_RESIDUE_CHAINID => info.chain_id = r.str()?,
            ATOM_PDB_RESIDUE_INSERTIONCODE => info.pdb_mut().insertion_code = r.str()?,
            ATOM_PDB_RESIDUE_OCCUPANCY => info.pdb_mut().occupancy = r.f64()?,
            ATOM_PDB_RESIDUE_TEMPFACTOR => info.pdb_mut().temp_factor = r.f64()?,
            ATOM_PDB_RESIDUE_ISHETEROATOM => info.pdb_mut().is_hetero_atom = r.u8()? != 0,
            ATOM_PDB_RESIDUE_SECONDARYSTRUCTURE => info.pdb_mut().secondary_structure = r.u32()?,
            ATOM_PDB_RESIDUE_SEGMENTNUMBER => info.pdb_mut().segment_number = r.u32()?,
            ATOM_PDB_RESIDUE_MONOMERCLASS => info.monomer_class = r.str()?,
            END_PDB_RESIDUE => break,
            _ => {
                return Err(pickle_error(
                    "unrecognized tag while parsing atom peptide residue info",
                ));
            }
        }
    }
    Ok(info)
}

/// `unpickleAtomMonomerInfo` (the `BEGIN_ATOM_MONOMER_INFO` tag is already read).
fn unpickle_monomer_info(r: &mut R<'_>, version: i32) -> Result<AtomMonomerInfo> {
    let name = r.str()?;
    let typ = r.u8()?;
    let mut info = AtomMonomerInfo::new(AtomMonomerType::from_u8(typ).unwrap_or_default(), &name);
    loop {
        let tag = r.tag(version)?;
        match tag {
            ATOM_MONOMER_INFO_RESIDUENAME => info.residue_name = r.str()?,
            ATOM_MONOMER_INFO_RESIDUENUMBER => info.residue_number = r.i32()?,
            ATOM_MONOMER_INFO_CHAINID => info.chain_id = r.str()?,
            ATOM_MONOMER_INFO_MONOMERCLASS => info.monomer_class = r.str()?,
            END_ATOM_MONOMER_INFO => break,
            _ => {
                return Err(pickle_error(format!(
                    "unrecognized tag while parsing atom monomer info {tag}"
                )));
            }
        }
    }
    Ok(info)
}

/// `_addAtomFromPickle`
fn atom_from_pickle(r: &mut R<'_>, version: i32) -> Result<(Atom, Option<Point3D>)> {
    let atomic_num = r.u8()? as u32;
    let flags = r.u8()?;
    let has_query = version > 5000 && flags & (1 << 4) != 0;
    let mut a = if has_query {
        let mut a = Atom::new(0);
        if atomic_num != 0 {
            a.set_atomic_num(atomic_num as i32);
        }
        a
    } else {
        Atom::new(atomic_num)
    };
    a.set_is_aromatic(flags & (1 << 6) != 0);
    a.set_no_implicit(flags & (1 << 5) != 0);
    let has_map = version >= 6020 && flags & (1 << 3) != 0;
    let has_dummy = version >= 6020 && flags & (1 << 2) != 0;
    let has_monomer = version >= 7020 && flags & (1 << 1) != 0;
    let mut pos = None;
    if flags & (1 << 7) != 0 {
        let x = r.f32()? as f64;
        let y = r.f32()? as f64;
        let z = r.f32()? as f64;
        pos = Some(Point3D::new(x, y, z));
    }
    if version <= 5000 || !has_query {
        if version < 7000 {
            if version < 6030 {
                r.i8()?;
            } else {
                r.f32()?;
            }
            a.set_formal_charge(r.i8()? as i32);
            a.set_chiral_tag(ChiralType::from_u8(r.u8()?).unwrap_or(ChiralType::Unspecified));
            a.set_hybridization(
                HybridizationType::from_u8(r.u8()?).unwrap_or(HybridizationType::Unspecified),
            );
            a.set_num_explicit_hs(r.i8()?.max(0) as u32);
            a.explicit_valence = r.i8()?;
            a.implicit_valence = r.i8()?;
            if version > 6000 {
                a.num_radical_electrons = r.u8()?;
            }
        } else {
            unpickle_atom_data(r, &mut a)?;
        }
    } else {
        if version >= 9000 {
            unpickle_atom_data(r, &mut a)?;
        }
        if r.tag(version)? != BEGINQUERY {
            return Err(pickle_error("Bad pickle format: BEGINQUERY tag not found."));
        }
        a.query = Some(Box::new(unpickle_atom_query(r, version)?));
        if r.tag(version)? != ENDQUERY {
            return Err(pickle_error("Bad pickle format: ENDQUERY tag not found."));
        }
    }
    if version > 5000 {
        if version < 6020 {
            let save = r.pos;
            if r.tag(version)? == ATOM_MAPNUMBER {
                let v = r.i8()? as i32;
                a.props.set_val(cp::molAtomMapNumber, v);
            } else {
                r.pos = save;
            }
        } else {
            if has_map {
                if r.tag(version)? != ATOM_MAPNUMBER {
                    return Err(pickle_error(
                        "Bad pickle format: ATOM_MAPNUMBER tag not found.",
                    ));
                }
                let c = r.i8()?;
                let v = if c < 0 && version > 9000 {
                    r.i32()?
                } else {
                    c as i32
                };
                a.props.set_val(cp::molAtomMapNumber, v);
            }
            if has_dummy {
                if r.tag(version)? != ATOM_DUMMYLABEL {
                    return Err(pickle_error(
                        "Bad pickle format: ATOM_DUMMYLABEL tag not found.",
                    ));
                }
                let s = r.str()?;
                a.props.set_val(cp::dummyLabel, s);
            }
        }
    }
    if has_monomer {
        let info = match r.tag(version)? {
            BEGIN_PDB_RESIDUE => unpickle_pdb_residue_info(r, version)?,
            BEGIN_ATOM_MONOMER_INFO => unpickle_monomer_info(r, version)?,
            _ => {
                return Err(pickle_error(
                    "Bad pickle format: BEGIN_PDB_RESIDUE or BEGIN_ATOM_MONOMER_INFO tag not found.",
                ));
            }
        };
        a.set_monomer_info(Some(info));
    }
    Ok((a, pos))
}

/// `_pickleBond`
fn pickle_bond(w: &mut W, b: &Bond, wide: bool) -> Result<()> {
    w.t(wide, b.begin_atom_idx() as i64);
    w.t(wide, b.end_atom_idx() as i64);
    let mut flags = 0u8;
    if b.is_aromatic() {
        flags |= 1 << 6;
    }
    if b.is_conjugated() {
        flags |= 1 << 5;
    }
    if b.has_query() {
        flags |= 1 << 4;
    }
    if b.bond_type() != BondType::Single {
        flags |= 1 << 3;
    }
    if b.bond_dir() != BondDir::None {
        flags |= 1 << 2;
    }
    if b.stereo() != BondStereo::StereoNone {
        flags |= 1 << 1;
    }
    let endpts = match b.get_prop(cp::_MolFileBondEndPts) {
        Some(RDValue::String(s)) => Some(s.clone()),
        _ => None,
    };
    if endpts.is_some() {
        flags |= 1;
    }
    w.u8(flags);
    if b.bond_type() != BondType::Single {
        w.u8(b.bond_type() as u8);
    }
    if b.bond_dir() != BondDir::None {
        w.u8(b.bond_dir() as u8);
    }
    if b.stereo() != BondStereo::StereoNone {
        w.u8(b.stereo() as u8);
        w.u8(b.stereo_atoms().len() as u8);
        for &s in b.stereo_atoms() {
            w.t(wide, s as i64);
        }
    }
    if let Some(q) = b.query() {
        w.u8(BEGINQUERY);
        pickle_query(w, q)?;
        w.u8(ENDQUERY);
    }
    if let Some(e) = endpts.filter(|e| !e.is_empty()) {
        w.str(&e);
        let attach = match b.get_prop(cp::_MolFileBondAttach) {
            Some(v) => v.to_rd_string(),
            None => "ALL".to_string(),
        };
        w.str(&attach);
    }
    Ok(())
}

/// `_addBondFromPickle` (directMap)
fn bond_from_pickle(r: &mut R<'_>, version: i32, wide: bool) -> Result<Bond> {
    let begin = r.t(wide)?;
    let end = r.t(wide)?;
    let flags = r.u8()?;
    let has_query = flags & (1 << 4) != 0;
    let mut b = Bond::new(BondType::Unspecified);
    if version <= 5000 || (version <= 7000 && !has_query) || version > 7000 {
        b.set_is_aromatic(flags & (1 << 6) != 0);
        b.set_is_conjugated(flags & (1 << 5) != 0);
        let bad = || pickle_error("invalid value in pickle");
        if version < 7000 {
            b.set_bond_type(BondType::from_u8(r.pos_char()?).ok_or_else(bad)?);
            b.set_bond_dir(BondDir::from_u8(r.pos_char()?).ok_or_else(bad)?);
            if version > 3000 {
                let stereo = BondStereo::from_u8(r.pos_char()?).ok_or_else(bad)?;
                b.stereo = stereo;
                if stereo != BondStereo::StereoNone {
                    let n = r.i8()?;
                    for _ in 0..n.max(0) {
                        let v = r.t(wide)? as i32;
                        b.stereo_atoms.push(v);
                    }
                }
            }
        } else {
            let bt = if flags & (1 << 3) != 0 {
                BondType::from_u8(r.pos_char()?).ok_or_else(bad)?
            } else {
                BondType::Single
            };
            b.set_bond_type(bt);
            let dir = if flags & (1 << 2) != 0 {
                BondDir::from_u8(r.pos_char()?).ok_or_else(bad)?
            } else {
                BondDir::None
            };
            b.set_bond_dir(dir);
            if flags & (1 << 1) != 0 {
                let stereo = BondStereo::from_u8(r.pos_char()?).ok_or_else(bad)?;
                let n = r.i8()?;
                for _ in 0..n.max(0) {
                    let v = r.t(wide)? as i32;
                    b.stereo_atoms.push(v);
                }
                b.stereo = stereo;
            } else {
                b.stereo = BondStereo::StereoNone;
            }
            if flags & 1 != 0 {
                let e = r.str()?;
                b.props.set_val(cp::_MolFileBondEndPts, e);
                let a = r.str()?;
                b.props.set_val(cp::_MolFileBondAttach, a);
            }
        }
    }
    if version > 5000 && has_query {
        if r.tag(version)? != BEGINQUERY {
            return Err(pickle_error("Bad pickle format: BEGINQUERY tag not found."));
        }
        b.query = Some(Box::new(unpickle_bond_query(r, version)?));
        if r.tag(version)? != ENDQUERY {
            return Err(pickle_error("Bad pickle format: ENDQUERY tag not found."));
        }
    }
    b.begin = begin as u32;
    b.end = end as u32;
    Ok(b)
}

// ================================================================== molecule

/// `MolPickler::_pickle<T>`
fn pickle_body(mol: &ROMol, w: &mut W, flags: u32, wide: bool) -> Result<()> {
    w.i32(mol.num_atoms() as i32);
    w.i32(mol.num_bonds() as i32);
    w.u8(1 << 7);
    w.u8(BEGINATOM);
    for a in mol.atoms() {
        pickle_atom(w, a)?;
    }
    w.u8(BEGINBOND);
    for b in mol.bonds() {
        pickle_bond(w, b, wide)?;
    }
    let ri = mol.ring_info();
    if ri.is_initialized() {
        w.u8(match ri.ring_type() {
            FindRingType::Fast => BEGINFASTFIND,
            FindRingType::Sssr => BEGINSSSR,
            FindRingType::SymmSssr => BEGINSYMMSSSR,
            FindRingType::OtherOrUnknown => BEGINFINDOTHERORUNKNOWN,
        });
        let rings = ri.atom_rings();
        w.u32(rings.len() as u32);
        for ring in rings {
            w.t(wide, ring.len() as i64);
            for &a in ring {
                w.t(wide, a as i64);
            }
        }
    }
    if !mol.substance_groups().is_empty() {
        w.u8(BEGINSGROUP);
        w.i32(mol.substance_groups().len() as i32);
        for sg in mol.substance_groups() {
            pickle_substance_group(w, sg, wide);
        }
    }
    if !mol.stereo_groups().is_empty() {
        w.u8(BEGINSTEREOGROUP);
        let mut groups = mol.stereo_groups().to_vec();
        assign_stereo_group_ids(&mut groups);
        w.t(wide, groups.len() as i64);
        for g in &groups {
            w.t(wide, g.group_type() as u8 as i64);
            if g.group_type() != StereoGroupType::Absolute {
                w.t(wide, g.write_id as i64);
            }
            w.t(wide, g.atoms().len() as i64);
            for &a in g.atoms() {
                w.t(wide, a as i64);
            }
            w.t(wide, g.bonds().len() as i64);
            for &b in g.bonds() {
                w.t(wide, b as i64);
            }
        }
    }
    if flags & NO_CONFORMERS == 0 {
        let mut t = W::default();
        let double = flags & COORDS_AS_DOUBLE != 0;
        w.u8(if double {
            BEGINCONFS_DOUBLE
        } else {
            BEGINCONFS
        });
        t.i32(mol.num_conformers() as i32);
        for conf in mol.conformers() {
            t.u8(conf.is_3d() as u8);
            t.i32(conf.id() as i32);
            t.t(wide, conf.num_atoms() as i64);
            for p in conf.positions() {
                if double {
                    t.f64(p.x);
                    t.f64(p.y);
                    t.f64(p.z);
                } else {
                    t.f32(p.x as f32);
                    t.f32(p.y as f32);
                    t.f32(p.z as f32);
                }
            }
        }
        if flags & MOL_PROPS != 0 {
            t.u8(BEGINCONFPROPS);
            let mut t2 = W::default();
            for conf in mol.conformers() {
                write_props(
                    &mut t2,
                    conf,
                    flags & PRIVATE_PROPS != 0,
                    flags & COMPUTED_PROPS != 0,
                    &[],
                    false,
                );
            }
            t.block(&t2);
        }
        w.block(&t);
    }
    if flags & MOL_PROPS != 0 {
        let mut t = W::default();
        write_props(
            &mut t,
            mol,
            flags & PRIVATE_PROPS != 0,
            flags & COMPUTED_PROPS != 0,
            &[],
            false,
        );
        if !t.0.is_empty() {
            w.u8(BEGINPROPS);
            w.block(&t);
            w.u8(ENDPROPS);
        }
    }
    if flags & ATOM_PROPS != 0 {
        let mut t = W::default();
        let ignore = atom_props_to_ignore();
        let mut any = false;
        for a in mol.atoms() {
            any |= write_props(
                &mut t,
                a,
                flags & PRIVATE_PROPS != 0,
                flags & COMPUTED_PROPS != 0,
                &ignore,
                false,
            );
            any |= write_explicit_props(&mut t, a, &EXPLICIT_ATOM_PROPS, true);
        }
        if any {
            w.u8(BEGINATOMPROPS);
            w.block(&t);
            w.u8(ENDPROPS);
        }
    }
    if flags & BOND_PROPS != 0 {
        let mut t = W::default();
        let ignore: Vec<&str> = EXPLICIT_BOND_PROPS.iter().map(|p| p.0).collect();
        let mut any = false;
        for b in mol.bonds() {
            any |= write_props(
                &mut t,
                b,
                flags & PRIVATE_PROPS != 0,
                flags & COMPUTED_PROPS != 0,
                &ignore,
                false,
            );
            any |= write_explicit_props(&mut t, b, &EXPLICIT_BOND_PROPS, false);
        }
        if any {
            w.u8(BEGINBONDPROPS);
            w.block(&t);
            w.u8(ENDPROPS);
        }
    }
    w.u8(ENDMOL);
    Ok(())
}

/// `MolPickler::pickleMol(mol, res, propertyFlags)`
pub fn pickle_mol(mol: &ROMol, flags: u32) -> Result<Vec<u8>> {
    let mut w = W::default();
    w.i32(ENDIAN_ID);
    w.i32(VERSION as i32);
    w.i32(VERSION_MAJOR);
    w.i32(VERSION_MINOR);
    w.i32(VERSION_PATCH);
    let wide = mol.num_atoms() > 255 || mol.num_bonds() > 255;
    pickle_body(mol, &mut w, flags, wide)?;
    Ok(w.0)
}

/// `MolPickler::_depickle<T>` (the molecule is new, so `directMap`)
fn depickle_body(
    r: &mut R<'_>,
    version: i32,
    num_atoms: i32,
    num_bonds: i32,
    flags: u32,
    wide: bool,
) -> Result<ROMol> {
    let mut mol = ROMol::new();
    let mut include_coords = false;
    if version >= 3000 {
        include_coords = r.u8()? & (1 << 7) != 0;
    }
    if r.tag(version)? != BEGINATOM {
        return Err(pickle_error("Bad pickle format: BEGINATOM tag not found."));
    }
    let old_coords = (2000..3000).contains(&version) && include_coords;
    let mut old_conf = old_coords.then(|| Conformer::new(num_atoms.max(0) as usize));
    let mut have_query = false;
    for i in 0..num_atoms.max(0) {
        let (a, pos) = atom_from_pickle(r, version)?;
        if let (Some(c), Some(p)) = (old_conf.as_mut(), pos) {
            c.positions_mut()[i as usize] = p;
        }
        have_query |= a.has_query();
        mol.add_atom(a, false);
    }
    if let Some(c) = old_conf {
        mol.add_conformer(c, true)?;
    }
    if r.tag(version)? != BEGINBOND {
        return Err(pickle_error("Bad pickle format: BEGINBOND tag not found."));
    }
    for _ in 0..num_bonds.max(0) {
        let b = bond_from_pickle(r, version, wide)?;
        mol.add_bond_object(b)?;
    }
    let mut tag = r.tag(version)?;
    let ring_type = match tag {
        BEGINSSSR => Some((FindRingType::Sssr, true)),
        BEGINSYMMSSSR => Some((FindRingType::SymmSssr, true)),
        BEGINFASTFIND => Some((FindRingType::Fast, false)),
        BEGINFINDOTHERORUNKNOWN => Some((FindRingType::OtherOrUnknown, false)),
        _ => None,
    };
    if let Some((rt, families)) = ring_type {
        read_rings(r, &mut mol, version, rt, wide)?;
        tag = r.tag(version)?;
        if families {
            find_ring_families(&mut mol, false, false)?;
        }
    }
    if tag == BEGINSGROUP {
        let n = r.i32()?;
        for _ in 0..n.max(0) {
            let sgroup = substance_group_from_pickle(r, &mol, version, wide)?;
            add_substance_group(&mut mol, sgroup);
        }
        tag = r.tag(version)?;
    }
    if tag == BEGINSTEREOGROUP {
        read_stereo_groups(r, &mut mol, version, wide)?;
        tag = r.tag(version)?;
    }
    if tag == BEGINCONFS || tag == BEGINCONFS_DOUBLE {
        let double = tag == BEGINCONFS_DOUBLE;
        let mut blk = 0;
        if version >= 13000 {
            blk = r.i32()?;
        }
        if version >= 13000 && flags & NO_CONFORMERS != 0 {
            r.skip(blk)?;
            tag = r.tag(version)?;
        } else {
            let n = r.i32()?;
            let mut cids = Vec::new();
            for _ in 0..n.max(0) {
                let is_3d = if version > 4000 { r.i8()? != 0 } else { true };
                let id = r.i32()? as u32;
                let na = r.t(wide)? as usize;
                let mut c = Conformer::new(na);
                c.set_id(id);
                c.is_3d = is_3d;
                for p in c.positions_mut().iter_mut() {
                    if double {
                        *p = Point3D::new(r.f64()?, r.f64()?, r.f64()?);
                    } else {
                        let (x, y, z) = (r.f32()?, r.f32()?, r.f32()?);
                        *p = Point3D::new(x as f64, y as f64, z as f64);
                    }
                }
                mol.add_conformer(c, false)?;
                cids.push(id);
            }
            tag = r.tag(version)?;
            if tag == BEGINCONFPROPS {
                let mut blk = 0;
                if version >= 13000 {
                    blk = r.i32()?;
                }
                if version >= 13000 && flags & MOL_PROPS == 0 {
                    r.skip(blk)?;
                } else {
                    for cid in cids {
                        let conf = mol.conformer_mut(cid as i32)?;
                        read_props(r, &mut conf.props, version < 14000, true)?;
                    }
                }
                tag = r.tag(version)?;
            }
        }
    }
    while tag != ENDMOL {
        match tag {
            BEGINPROPS => {
                let mut blk = 0;
                if version >= 13000 {
                    blk = r.i32()?;
                }
                if version >= 13000 && flags & MOL_PROPS == 0 {
                    r.skip(blk)?;
                } else {
                    read_props(r, &mut mol.props, version < 14000, true)?;
                }
            }
            BEGINATOMPROPS => {
                let mut blk = 0;
                if version >= 13000 {
                    blk = r.i32()?;
                }
                if version >= 13000 && flags & ATOM_PROPS == 0 {
                    r.skip(blk)?;
                } else {
                    for i in 0..mol.num_atoms() {
                        let a = mol.atom_mut(i);
                        read_props(r, &mut a.props, version < 14000, false)?;
                        read_explicit_props(r, &mut a.props, &EXPLICIT_ATOM_PROPS, true, version)?;
                    }
                }
            }
            BEGINBONDPROPS => {
                let mut blk = 0;
                if version >= 13000 {
                    blk = r.i32()?;
                }
                if version >= 13000 && flags & BOND_PROPS == 0 {
                    r.skip(blk)?;
                } else {
                    for i in 0..mol.num_bonds() {
                        let b = mol.bond_mut(i);
                        read_props(r, &mut b.props, version < 14000, false)?;
                        read_explicit_props(r, &mut b.props, &EXPLICIT_BOND_PROPS, false, version)?;
                    }
                }
            }
            BEGINQUERYATOMDATA => {
                for i in 0..mol.num_atoms() {
                    unpickle_atom_data(r, mol.atom_mut(i))?;
                }
            }
            _ => break,
        }
        if r.tag(version)? != ENDPROPS {
            return Err(pickle_error("Bad pickle format: ENDPROPS tag not found."));
        }
        tag = r.tag(version)?;
    }
    if tag != ENDMOL {
        return Err(pickle_error("Bad pickle format: ENDMOL tag not found."));
    }
    if have_query {
        for i in 0..mol.num_atoms() {
            if mol.atom(i).has_query() {
                mol.update_atom_property_cache(i, false)?;
            }
        }
    }
    Ok(mol)
}

/// `_addRingInfoFromPickle` (directMap)
fn read_rings(
    r: &mut R<'_>,
    mol: &mut ROMol,
    version: i32,
    rt: FindRingType,
    wide: bool,
) -> Result<()> {
    mol.ring_info_mut().initialize(rt);
    let n = if version >= 13002 {
        r.u32()?
    } else {
        r.t(wide)? as u32
    };
    if n == 0 {
        return Ok(());
    }
    let (na, nb) = (mol.num_atoms() as usize, mol.num_bonds() as usize);
    mol.ring_info_mut().preallocate(na, nb);
    for _ in 0..n {
        let size = r.t(wide)?;
        if size < 0 {
            return Err(pickle_error("negative ring size"));
        }
        let size = size as usize;
        let mut atoms = Vec::with_capacity(size);
        for _ in 0..size {
            let a = r.t(wide)?;
            if a < 0 || a as usize >= na {
                return Err(pickle_error("ring-atom index out of range"));
            }
            atoms.push(a as i32);
        }
        let mut bonds = vec![0i32; size];
        if version < 7000 {
            for b in bonds.iter_mut() {
                let v = r.t(wide)?;
                if v < 0 || v as usize >= nb {
                    return Err(pickle_error("ring-bond index out of range"));
                }
                *b = v as i32;
            }
        } else if size > 0 {
            let between = |x: i32, y: i32| -> Result<i32> {
                mol.bond_idx_between(x as u32, y as u32)
                    .map(|b| b as i32)
                    .ok_or_else(|| pickle_error("ring bond not found"))
            };
            for j in 1..size {
                bonds[j - 1] = between(atoms[j - 1], atoms[j])?;
            }
            bonds[size - 1] = between(atoms[0], atoms[size - 1])?;
        }
        mol.ring_info_mut().add_ring(atoms, bonds)?;
    }
    Ok(())
}

/// `_depickleStereo`
/// `MolPickler::_pickleSubstanceGroup<T>`
fn pickle_substance_group(w: &mut W, sgroup: &SubstanceGroup, wide: bool) {
    write_props(w, sgroup, false, false, &[], false);
    for indices in [sgroup.atoms(), sgroup.parent_atoms(), sgroup.bonds()] {
        w.t(wide, indices.len() as i64);
        for &i in indices {
            w.t(wide, i as i64);
        }
    }
    w.t(wide, sgroup.brackets().len() as i64);
    for bracket in sgroup.brackets() {
        for pt in bracket {
            write_point_f32(w, pt);
        }
    }
    // the vector of a CSTATE only exists for SUP groups
    let is_sup = is_sup_group(sgroup);
    w.t(wide, sgroup.cstates().len() as i64);
    for cstate in sgroup.cstates() {
        w.t(wide, cstate.bond_idx as i64);
        if is_sup {
            write_point_f32(w, &cstate.vector);
        }
    }
    w.t(wide, sgroup.attach_points().len() as i64);
    for ap in sgroup.attach_points() {
        w.t(wide, ap.a_idx as i64);
        w.i32(ap.lv_idx);
        w.str(&ap.id);
    }
}

fn is_sup_group(sgroup: &SubstanceGroup) -> bool {
    matches!(sgroup.get_prop("TYPE"), Some(RDValue::String(t)) if t == "SUP")
}

fn write_point_f32(w: &mut W, pt: &Point3D) {
    w.f32(pt.x as f32);
    w.f32(pt.y as f32);
    w.f32(pt.z as f32);
}

fn read_point_f32(r: &mut R<'_>) -> Result<Point3D> {
    let x = r.f32()? as f64;
    let y = r.f32()? as f64;
    let z = r.f32()? as f64;
    Ok(Point3D::new(x, y, z))
}

/// `MolPickler::_getSubstanceGroupFromPickle<T>`
fn substance_group_from_pickle(
    r: &mut R<'_>,
    mol: &ROMol,
    version: i32,
    wide: bool,
) -> Result<SubstanceGroup> {
    // the TYPE property is overridden by the pickled properties
    let mut sgroup = SubstanceGroup::new("");
    read_props(r, sgroup.dict_mut(), version < 14000, true)?;

    for _ in 0..r.t(wide)? {
        sgroup.add_atom_with_idx(mol, r.t(wide)? as u32)?;
    }
    for _ in 0..r.t(wide)? {
        sgroup.add_parent_atom_with_idx(mol, r.t(wide)? as u32)?;
    }
    for _ in 0..r.t(wide)? {
        sgroup.add_bond_with_idx(mol, r.t(wide)? as u32)?;
    }
    for _ in 0..r.t(wide)? {
        let bracket = [read_point_f32(r)?, read_point_f32(r)?, read_point_f32(r)?];
        sgroup.add_bracket(bracket);
    }
    let is_sup = is_sup_group(&sgroup);
    for _ in 0..r.t(wide)? {
        let bond = r.t(wide)? as u32;
        let vector = if is_sup {
            read_point_f32(r)?
        } else {
            Point3D::default()
        };
        sgroup.add_cstate(mol, bond, vector)?;
    }
    for _ in 0..r.t(wide)? {
        let a_idx = r.t(wide)? as u32;
        let lv_idx = r.i32()?;
        let id = r.str()?;
        sgroup.add_attach_point(a_idx, lv_idx, &id);
    }
    Ok(sgroup)
}

fn read_stereo_groups(r: &mut R<'_>, mol: &mut ROMol, version: i32, wide: bool) -> Result<()> {
    let n = r.t(wide)? as u32;
    if n == 0 {
        return Ok(());
    }
    let mut groups = Vec::new();
    for _ in 0..n {
        let gt = StereoGroupType::from_u8(r.t(wide)? as u8)
            .ok_or_else(|| pickle_error("bad stereo group type"))?;
        let mut gid = 0u32;
        if version >= 14010 && gt != StereoGroupType::Absolute {
            gid = r.t(wide)? as u32;
        }
        let na = r.t(wide)? as u32;
        let mut atoms = Vec::with_capacity(na as usize);
        for _ in 0..na {
            let a = r.t(wide)? as u32;
            mol.get_atom(a)?;
            atoms.push(a);
        }
        let mut bonds = Vec::new();
        if version > 16000 {
            let nb = r.t(wide)? as u32;
            for _ in 0..nb {
                let b = r.t(wide)? as u32;
                mol.get_bond(b)?;
                bonds.push(b);
            }
        }
        groups.push(StereoGroup::new(gt, atoms, bonds, gid)?);
    }
    mol.set_stereo_groups(groups)
}

fn mol_from_pickle_reader(r: &mut R<'_>, flags: u32) -> Result<ROMol> {
    if r.i32()? != ENDIAN_ID {
        return Err(pickle_error(
            "Bad pickle format: bad endian ID or invalid file format",
        ));
    }
    if r.i32()? != VERSION as i32 {
        return Err(pickle_error("Bad pickle format: no version tag"));
    }
    let (major, minor, patch) = (r.i32()?, r.i32()?, r.i32()?);
    if major > VERSION_MAJOR || (major == VERSION_MAJOR && minor > VERSION_MINOR) {
        crate::rd_log_warning!(
            "Depickling from a version number ({major}.{minor})that is higher than our version ({VERSION_MAJOR}.{VERSION_MINOR}).\nThis probably won't work."
        );
    }
    if major > 1000 || minor > 100 || patch > 100 {
        return Err(pickle_error("unreasonable version numbers"));
    }
    let version = 1000 * major + minor * 10 + patch;
    if version == 1 {
        return Err(pickle_error("rdoxide: version-1 pickles are not supported"));
    }
    let num_atoms = r.i32()?;
    let num_bonds = r.i32()?;
    let wide = num_atoms > 255 || (version >= 16040 && num_bonds > 255);
    let mut mol = depickle_body(r, version, num_atoms, num_bonds, flags, wide)?;
    mol.clear_all_atom_bookmarks();
    mol.clear_all_bond_bookmarks();
    if version < 4000 {
        crate::graphmol::chirality::legacy::assign_stereochemistry(&mut mol, true, false, false)?;
    }
    Ok(mol)
}

/// `streamWriteProps`: the properties selected by `flags` (`PRIVATE_PROPS`,
/// `COMPUTED_PROPS`), for pickles that embed properties outside a molecule.
pub(crate) fn pickle_properties(props: &impl RDProps, flags: u32) -> Vec<u8> {
    let mut w = W::default();
    write_props(
        &mut w,
        props,
        flags & PRIVATE_PROPS != 0,
        flags & COMPUTED_PROPS != 0,
        &[],
        true,
    );
    w.0
}

/// `streamReadProps`: reads properties written by [`pickle_properties`] from
/// `data` at `*pos` into `dict` and advances `*pos`.
pub(crate) fn unpickle_properties(data: &[u8], pos: &mut usize, dict: &mut Dict) -> Result<()> {
    let mut r = R { data, pos: *pos };
    read_props(&mut r, dict, true, false)?;
    *pos = r.pos;
    Ok(())
}

/// `MolPickler::molFromPickle(pickle, mol, propertyFlags)`
pub fn mol_from_pickle(data: &[u8], flags: u32) -> Result<ROMol> {
    let mut r = R { data, pos: 0 };
    mol_from_pickle_reader(&mut r, flags)
}

/// `MolPickler::molFromPickle(std::istream &, ...)`: reads one pickle starting
/// at `*pos` and advances `*pos` past it (pickles can be concatenated).
pub fn mol_from_pickle_stream(data: &[u8], pos: &mut usize, flags: u32) -> Result<ROMol> {
    let mut r = R { data, pos: *pos };
    let m = mol_from_pickle_reader(&mut r, flags)?;
    *pos = r.pos;
    Ok(m)
}
