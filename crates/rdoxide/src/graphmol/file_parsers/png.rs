//! Molecule metadata in PNG files.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/PNGParser.{h,cpp}` (2026.09.1).
//! Molecules are stored in `tEXt`/`zTXt` chunks; the zlib streams of the
//! latter are handled by the small `inflate`/`zlib_store` helpers below, so
//! no compression library is needed.
//!
//! Copyright (C) 2020 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::graphmol::file_parsers::mol_file_parser::{MolFileParserParams, mol_from_mol_block};
use crate::graphmol::file_parsers::mol_file_writer::{MolWriterParams, mol_to_mol_block};
use crate::graphmol::pickle::flags::ALL_PROPS;
use crate::graphmol::pickle::{get_default_pickle_properties, mol_from_pickle, pickle_mol};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::cx_write::{RestoreBondDirOption, cx_fields, mol_to_cx_smiles};
use crate::graphmol::smiles_parse::write::SmilesWriteParams;
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::graphmol::wedge_bonds::reapply_mol_block_wedging;
use crate::rdgeneral::RDKIT_VERSION;
use crate::{RDKitError, Result};

/// `PNGData::smilesTag`
pub const SMILES_TAG: &str = "SMILES";
/// `PNGData::molTag`
pub const MOL_TAG: &str = "MOL";
/// `PNGData::pklTag`
pub const PKL_TAG: &str = "rdkitPKL";

/// The text chunks of a PNG: keys with their (possibly binary) values.
pub type PngMetadata = Vec<(String, Vec<u8>)>;

const PNG_HEADER: [u8; 8] = [137, 80, 78, 71, 13, 10, 26, 10];

/// `PNGMetadataParams`
#[derive(Debug, Clone)]
pub struct PngMetadataParams {
    /// include the molecule pickle
    pub include_pkl: bool,
    /// include CXSMILES for the molecule
    pub include_smiles: bool,
    /// include a mol block for the molecule
    pub include_mol: bool,
    /// properties included in the pickle
    pub property_flags: u32,
    /// parameters of the CXSMILES
    pub smiles_write_params: SmilesWriteParams,
    /// CXSMILES fields to include
    pub cx_smiles_flags: u32,
    /// what to do with bond directions in the CXSMILES
    pub restore_bond_dirs: RestoreBondDirOption,
}

impl Default for PngMetadataParams {
    fn default() -> Self {
        PngMetadataParams {
            include_pkl: true,
            include_smiles: true,
            include_mol: false,
            property_flags: get_default_pickle_properties(),
            smiles_write_params: SmilesWriteParams::default(),
            cx_smiles_flags: cx_fields::ALL,
            restore_bond_dirs: RestoreBondDirOption::Clear,
        }
    }
}

fn fpe(msg: &str) -> RDKitError {
    RDKitError::FileParse(msg.into())
}

/// `augmentTagName`
fn augment_tag_name(tag: &str) -> String {
    format!("{tag} rdkit {RDKIT_VERSION}")
}

// ------------------------------------------------------------------ zlib

/// CRC-32 of the PNG chunks (`boost::crc_32_type`).
fn crc32(data: &[u8]) -> u32 {
    let mut crc = !0u32;
    for &b in data {
        crc ^= u32::from(b);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xEDB8_8320 & 0u32.wrapping_sub(crc & 1));
        }
    }
    !crc
}

fn adler32(data: &[u8]) -> u32 {
    let (mut a, mut b) = (1u32, 0u32);
    for chunk in data.chunks(5552) {
        for &x in chunk {
            a += u32::from(x);
            b += a;
        }
        a %= 65521;
        b %= 65521;
    }
    (b << 16) | a
}

/// A zlib stream holding `data` in uncompressed ("stored") deflate blocks.
fn zlib_store(data: &[u8]) -> Vec<u8> {
    let mut out = vec![0x78, 0x01];
    let mut blocks = data.chunks(0xFFFF).peekable();
    if blocks.peek().is_none() {
        out.extend([1, 0, 0, 0xFF, 0xFF]);
    }
    while let Some(block) = blocks.next() {
        let len = block.len() as u16;
        out.push(u8::from(blocks.peek().is_none()));
        out.extend(len.to_le_bytes());
        out.extend((!len).to_le_bytes());
        out.extend(block);
    }
    out.extend(adler32(data).to_be_bytes());
    out
}

/// LSB-first bit reader of a deflate stream.
struct Bits<'a> {
    data: &'a [u8],
    pos: usize,
    bit: u32,
}

impl Bits<'_> {
    fn read(&mut self, n: u32) -> Option<u32> {
        let mut v = 0;
        for i in 0..n {
            let byte = *self.data.get(self.pos)?;
            v |= u32::from((byte >> self.bit) & 1) << i;
            self.bit += 1;
            if self.bit == 8 {
                self.bit = 0;
                self.pos += 1;
            }
        }
        Some(v)
    }
}

/// Canonical Huffman code: symbol counts per length and symbols by code.
struct Huffman {
    counts: [u16; 16],
    symbols: Vec<u16>,
}

impl Huffman {
    fn new(lengths: &[u8]) -> Huffman {
        let mut counts = [0u16; 16];
        for &l in lengths {
            counts[l as usize] += 1;
        }
        counts[0] = 0;
        let mut offsets = [0u16; 16];
        for i in 1..16 {
            offsets[i] = offsets[i - 1] + counts[i - 1];
        }
        let mut symbols = vec![0u16; lengths.len()];
        for (sym, &l) in lengths.iter().enumerate() {
            if l != 0 {
                symbols[offsets[l as usize] as usize] = sym as u16;
                offsets[l as usize] += 1;
            }
        }
        Huffman { counts, symbols }
    }

    fn decode(&self, bits: &mut Bits) -> Option<u16> {
        let (mut code, mut first, mut index) = (0i32, 0i32, 0i32);
        for len in 1..16 {
            code |= bits.read(1)? as i32;
            let count = i32::from(self.counts[len]);
            if code - count < first {
                return self.symbols.get((index + (code - first)) as usize).copied();
            }
            index += count;
            first = (first + count) << 1;
            code <<= 1;
        }
        None
    }
}

const LENGTH_BASE: [u16; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];
const LENGTH_EXTRA: [u8; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
const DIST_BASE: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
];
const DIST_EXTRA: [u8; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];

fn inflate_block(bits: &mut Bits, out: &mut Vec<u8>, lit: &Huffman, dist: &Huffman) -> Option<()> {
    loop {
        let sym = lit.decode(bits)? as usize;
        match sym {
            0..=255 => out.push(sym as u8),
            256 => return Some(()),
            _ => {
                let i = sym - 257;
                let len =
                    *LENGTH_BASE.get(i)? as usize + bits.read(u32::from(LENGTH_EXTRA[i]))? as usize;
                let d = dist.decode(bits)? as usize;
                let back =
                    *DIST_BASE.get(d)? as usize + bits.read(u32::from(DIST_EXTRA[d]))? as usize;
                if back > out.len() {
                    return None;
                }
                for _ in 0..len {
                    out.push(out[out.len() - back]);
                }
            }
        }
    }
}

/// The contents of a zlib stream, or `None` when it is malformed.
fn zlib_inflate(data: &[u8]) -> Option<Vec<u8>> {
    if data.len() < 2
        || data[0] & 0x0F != 8
        || (u32::from(data[0]) << 8 | u32::from(data[1])) % 31 != 0
    {
        return None;
    }
    let mut bits = Bits {
        data,
        pos: 2,
        bit: 0,
    };
    let mut out = Vec::new();
    loop {
        let last = bits.read(1)?;
        match bits.read(2)? {
            0 => {
                if bits.bit != 0 {
                    bits.bit = 0;
                    bits.pos += 1;
                }
                let head = data.get(bits.pos..bits.pos + 4)?;
                let len = usize::from(u16::from_le_bytes([head[0], head[1]]));
                bits.pos += 4;
                out.extend(data.get(bits.pos..bits.pos + len)?);
                bits.pos += len;
            }
            1 => {
                let mut lengths = [8u8; 288];
                lengths[144..256].fill(9);
                lengths[256..280].fill(7);
                inflate_block(
                    &mut bits,
                    &mut out,
                    &Huffman::new(&lengths),
                    &Huffman::new(&[5; 30]),
                )?;
            }
            2 => {
                let nlen = bits.read(5)? as usize + 257;
                let ndist = bits.read(5)? as usize + 1;
                let ncode = bits.read(4)? as usize + 4;
                const ORDER: [usize; 19] = [
                    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
                ];
                let mut code_lengths = [0u8; 19];
                for &o in ORDER.iter().take(ncode) {
                    code_lengths[o] = bits.read(3)? as u8;
                }
                let code = Huffman::new(&code_lengths);
                let mut lengths = Vec::with_capacity(nlen + ndist);
                while lengths.len() < nlen + ndist {
                    let sym = code.decode(&mut bits)?;
                    let (value, repeat) = match sym {
                        0..=15 => (sym as u8, 1),
                        16 => (*lengths.last()?, 3 + bits.read(2)?),
                        17 => (0, 3 + bits.read(3)?),
                        _ => (0, 11 + bits.read(7)?),
                    };
                    lengths.extend(std::iter::repeat_n(value, repeat as usize));
                }
                if lengths.len() > nlen + ndist {
                    return None;
                }
                inflate_block(
                    &mut bits,
                    &mut out,
                    &Huffman::new(&lengths[..nlen]),
                    &Huffman::new(&lengths[nlen..]),
                )?;
            }
            _ => return None,
        }
        if last == 1 {
            return Some(out);
        }
    }
}

// ------------------------------------------------------------------- PNG

fn check_header(data: &[u8]) -> Result<&[u8]> {
    data.strip_prefix(&PNG_HEADER)
        .ok_or_else(|| fpe("PNG header not recognized"))
}

/// One chunk: its type, its data and the bytes after it (CRC included).
fn next_chunk(rest: &[u8]) -> Option<([u8; 4], &[u8], &[u8])> {
    let len = u32::from_be_bytes(rest.get(..4)?.try_into().ok()?) as usize;
    let kind: [u8; 4] = rest.get(4..8)?.try_into().ok()?;
    let body = rest.get(8..8usize.checked_add(len)?)?;
    Some((kind, body, rest.get(8 + len..)?))
}

/// `PNGStreamToMetadata`: the `tEXt` and `zTXt` chunks of a PNG.
pub fn png_to_metadata(data: &[u8]) -> Result<PngMetadata> {
    let mut rest = check_header(data)?;
    let mut res = PngMetadata::new();
    loop {
        let (kind, body, after) =
            next_chunk(rest).ok_or_else(|| fpe("error when reading from PNG"))?;
        if &kind == b"IEND" {
            break;
        }
        let compressed = &kind == b"zTXt";
        if !body.is_empty() && (compressed || &kind == b"tEXt") {
            let nul = body
                .iter()
                .position(|&b| b == 0)
                .ok_or_else(|| fpe("error when reading from PNG"))?;
            let key = String::from_utf8_lossy(&body[..nul]).into_owned();
            let value = &body[nul + 1..];
            let value = if compressed {
                // the first byte is the compression method
                zlib_inflate(value.get(1..).unwrap_or_default()).unwrap_or_default()
            } else {
                value.to_vec()
            };
            if !value.is_empty() {
                res.push((key, value));
            }
        }
        // the CRC is skipped
        rest = after.get(4..).unwrap_or_default();
    }
    Ok(res)
}

/// `addMetadataToPNGStream`: the PNG with `metadata` added before `IEND`.
pub fn add_metadata_to_png(
    data: &[u8],
    metadata: &[(String, Vec<u8>)],
    compressed: bool,
) -> Result<Vec<u8>> {
    let mut rest = check_header(data)?;
    let mut res = PNG_HEADER.to_vec();
    let final_crc = loop {
        let (kind, body, after) =
            next_chunk(rest).ok_or_else(|| fpe("did not find IEND block in PNG"))?;
        if &kind == b"IEND" {
            break after
                .get(..4)
                .ok_or_else(|| fpe("did not find IEND block in PNG"))?;
        }
        // chunk header, data and CRC, unchanged
        let total = 8 + body.len() + 4;
        res.extend(
            rest.get(..total)
                .ok_or_else(|| fpe("did not find IEND block in PNG"))?,
        );
        rest = after.get(4..).unwrap_or_default();
    };
    for (key, value) in metadata {
        let mut blob = Vec::new();
        blob.extend(if compressed { b"zTXt" } else { b"tEXt" });
        blob.extend(key.as_bytes());
        blob.push(0);
        if compressed {
            blob.push(0);
            blob.extend(zlib_store(value));
        } else {
            blob.extend(value);
        }
        res.extend(((blob.len() - 4) as u32).to_be_bytes());
        res.extend(&blob);
        res.extend(crc32(&blob).to_be_bytes());
    }
    res.extend(0u32.to_be_bytes());
    res.extend(b"IEND");
    res.extend(final_crc);
    Ok(res)
}

/// `addMolToPNGStream`: the PNG with the molecule added as metadata.
pub fn add_mol_to_png(mol: &ROMol, data: &[u8], params: &PngMetadataParams) -> Result<Vec<u8>> {
    let mut metadata = PngMetadata::new();
    if params.include_pkl {
        let pkl = pickle_mol(mol, params.property_flags)?;
        metadata.push((augment_tag_name(PKL_TAG), pkl));
    }
    if params.include_smiles {
        let smi = mol_to_cx_smiles(
            mol,
            &params.smiles_write_params,
            params.cx_smiles_flags,
            params.restore_bond_dirs,
        )?;
        metadata.push((augment_tag_name(SMILES_TAG), smi.into_bytes()));
    }
    if params.include_mol {
        let mut rewedged;
        let mut mol_ref = mol;
        if params.restore_bond_dirs == RestoreBondDirOption::True {
            rewedged = mol.clone();
            reapply_mol_block_wedging(&mut rewedged, true, false)?;
            mol_ref = &rewedged;
        }
        let writer = MolWriterParams {
            include_stereo: true,
            kekulize: false,
            ..MolWriterParams::default()
        };
        let mb = mol_to_mol_block(mol_ref, &writer, -1)?;
        metadata.push((augment_tag_name(MOL_TAG), mb.into_bytes()));
    }
    add_metadata_to_png(data, &metadata, true)
}

fn mol_from_entry(
    key: &str,
    value: &[u8],
    params: &SmilesParserParams,
) -> Result<Option<Option<ROMol>>> {
    let text = || String::from_utf8_lossy(value);
    Ok(if key.starts_with(PKL_TAG) {
        Some(Some(mol_from_pickle(value, ALL_PROPS)?))
    } else if key.starts_with(SMILES_TAG) {
        Some(mol_from_smiles(&text(), params)?)
    } else if key.starts_with(MOL_TAG) {
        let p = MolFileParserParams {
            sanitize: params.sanitize,
            remove_hs: params.remove_hs,
            ..MolFileParserParams::default()
        };
        Some(mol_from_mol_block(&text(), &p)?)
    } else {
        None
    })
}

/// `PNGStreamToMol`: the molecule of the first molecule tag of the PNG.
pub fn png_to_mol(data: &[u8], params: &SmilesParserParams) -> Result<Option<ROMol>> {
    for (key, value) in png_to_metadata(data)? {
        if let Some(mol) = mol_from_entry(&key, &value, params)? {
            return Ok(mol);
        }
    }
    Err(fpe("No suitable metadata found."))
}

/// `PNGStreamToMols`: the molecules of all tags starting with `tag`.
pub fn png_to_mols(data: &[u8], tag: &str, params: &SmilesParserParams) -> Result<Vec<ROMol>> {
    let mut res = Vec::new();
    for (key, value) in png_to_metadata(data)? {
        if !key.starts_with(tag) {
            continue;
        }
        if let Some(Some(mol)) = mol_from_entry(&key, &value, params)? {
            res.push(mol);
        }
    }
    Ok(res)
}
