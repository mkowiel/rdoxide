//! Reader for chemfp's FPB fingerprint files (upstream
//! `DataStructs/FPBReader.{h,cpp}`).
//!
//! An FPB file is a magic number followed by tagged chunks: `AREN` (the
//! fingerprints), `FPID` (their identifiers), `POPC` (offsets that bucket the
//! fingerprints by population count) and a few ignored ones, ending with
//! `FEND`. The reader either loads everything into memory or, in lazy mode,
//! only remembers where the `AREN` and `FPID` data sit and seeks to it on
//! demand.

use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::Mutex;

use super::bitmap::{
    bitmap_all_probe_bits_match, bitmap_popcount, bitmap_tanimoto, bitmap_tversky,
};
use super::{BitVect, ExplicitBitVect};
use crate::error::{RDKitError, Result};

const MAGIC: &[u8; 8] = b"FPB1\r\n\0\0";
/// Fingerprints fetched per read when scanning a lazy file for neighbors.
const READ_CACHE: usize = 1000;

/// A byte source an [`FpbReader`] can read from.
pub trait FpbSource: Read + Seek + Send {}
impl<T: Read + Seek + Send> FpbSource for T {}

/// A `(similarity, fingerprint index)` search hit.
pub type Neighbor = (f64, u32);

/// Where the fingerprints live once the file has been read.
#[derive(Debug)]
enum FpData {
    Memory(Vec<u8>),
    /// Absolute file offset of the first fingerprint.
    Lazy(u64),
}

/// Where the identifiers live once the file has been read.
#[derive(Debug)]
enum IdData {
    /// The whole `FPID` chunk; `offsets` is where its offset table starts.
    Memory { chunk: Vec<u8>, offsets: usize },
    /// Absolute file offsets of the chunk and of its offset table.
    Lazy { chunk: u64, offsets: u64 },
}

/// Everything learned from the file by [`FpbReader::init`].
#[derive(Debug)]
struct Loaded {
    len: u32,
    n_bits: u32,
    /// Bytes stored per fingerprint (the fingerprint plus alignment padding).
    stored_bytes: usize,
    pop_count_offsets: Vec<u32>,
    num_4byte_ids: u32,
    num_8byte_ids: u32,
    fps: FpData,
    ids: IdData,
}

/// Reader and searcher for FPB files.
///
/// All searching methods take `&self`; in lazy mode the underlying stream is
/// serialized with a mutex, so a reader can be shared between threads.
pub struct FpbReader {
    source: Mutex<Box<dyn FpbSource>>,
    lazy: bool,
    loaded: Option<Loaded>,
}

impl std::fmt::Debug for FpbReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FpbReader")
            .field("lazy", &self.lazy)
            .field("loaded", &self.loaded)
            .finish_non_exhaustive()
    }
}

fn bad_file(msg: &str) -> RDKitError {
    RDKitError::BadFile(msg.to_string())
}

fn read_u32(r: &mut (impl Read + ?Sized)) -> std::io::Result<u32> {
    let mut b = [0; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn le_u32(bytes: &[u8], at: usize) -> Option<u32> {
    Some(u32::from_le_bytes(bytes.get(at..at + 4)?.try_into().ok()?))
}

fn truncated(_: std::io::Error) -> RDKitError {
    bad_file("EOF hit before FEND record")
}

/// A chunk header: the size of the data and its four-letter tag.
fn read_chunk_header(r: &mut impl Read) -> Result<(u64, [u8; 4])> {
    let mut size = [0; 8];
    r.read_exact(&mut size).map_err(truncated)?;
    let mut tag = [0; 4];
    r.read_exact(&mut tag).map_err(truncated)?;
    Ok((u64::from_le_bytes(size), tag))
}

fn read_chunk_data(r: &mut impl Read, size: u64) -> Result<Vec<u8>> {
    let size = usize::try_from(size).map_err(|_| bad_file("chunk too large"))?;
    let mut data = vec![0; size];
    r.read_exact(&mut data).map_err(truncated)?;
    Ok(data)
}

fn parse_pop_counts(chunk: &[u8]) -> Result<Vec<u32>> {
    if chunk.len() % 4 != 0 {
        return Err(RDKitError::Value(
            "POPC chunk size must be a multiple of 4 bytes".into(),
        ));
    }
    if chunk.len() / 4 < 9 {
        return Err(RDKitError::Value(
            "POPC must contain at least 9 offsets".into(),
        ));
    }
    Ok(chunk
        .chunks_exact(4)
        .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect())
}

/// The fixed-size header of the `AREN` chunk.
#[derive(Debug, Clone, Copy)]
struct ArenaHeader {
    n_bytes: u32,
    stored_bytes: u32,
    spacer: u8,
}

impl ArenaHeader {
    const SIZE: u64 = 9;

    fn from_bytes(b: &[u8; 9]) -> Self {
        ArenaHeader {
            n_bytes: u32::from_le_bytes([b[0], b[1], b[2], b[3]]),
            stored_bytes: u32::from_le_bytes([b[4], b[5], b[6], b[7]]),
            spacer: b[8],
        }
    }

    /// Number of fingerprints in an arena chunk of `chunk_size` bytes.
    fn len(&self, chunk_size: u64) -> Result<u32> {
        if self.stored_bytes == 0 || chunk_size < Self::SIZE + u64::from(self.spacer) {
            return Err(bad_file("Invalid AREN record"));
        }
        let n = (chunk_size - Self::SIZE - u64::from(self.spacer)) / u64::from(self.stored_bytes);
        u32::try_from(n).map_err(|_| bad_file("Invalid AREN record"))
    }
}

impl FpbReader {
    /// Opens `path` for reading; call [`FpbReader::init`] before using it.
    ///
    /// With `lazy` set only the file layout is read by `init`; fingerprints
    /// and identifiers are fetched from the file on demand.
    pub fn open(path: impl AsRef<Path>, lazy: bool) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)
            .map_err(|_| RDKitError::BadFile(format!("Bad input file {}", path.display())))?;
        Ok(Self::from_source(BufReader::new(file), lazy))
    }

    /// Reads from an arbitrary seekable stream.
    pub fn from_source(source: impl FpbSource + 'static, lazy: bool) -> Self {
        FpbReader {
            source: Mutex::new(Box::new(source)),
            lazy,
            loaded: None,
        }
    }

    fn stream(&self) -> std::sync::MutexGuard<'_, Box<dyn FpbSource>> {
        self.source.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// Reads the file layout (and, unless lazy, all fingerprints and ids).
    /// Calling it again on an initialized reader does nothing.
    pub fn init(&mut self) -> Result<()> {
        if self.loaded.is_some() {
            return Ok(());
        }
        let loaded = self.read_file()?;
        self.loaded = Some(loaded);
        Ok(())
    }

    /// Releases what [`FpbReader::init`] loaded.
    pub fn cleanup(&mut self) {
        self.loaded = None;
    }

    pub fn is_initialized(&self) -> bool {
        self.loaded.is_some()
    }

    fn read_file(&self) -> Result<Loaded> {
        let mut r = self.stream();
        let mut magic = [0; 8];
        r.read_exact(&mut magic)
            .map_err(|_| bad_file("Invalid FPB magic"))?;
        if &magic != MAGIC {
            return Err(bad_file("Invalid FPB magic"));
        }
        let mut pop_count_offsets = Vec::new();
        let mut arena: Option<(ArenaHeader, u32, FpData)> = None;
        let mut ids: Option<(u32, u32, IdData)> = None;
        loop {
            let (size, tag) = read_chunk_header(&mut *r)?;
            match (&tag, self.lazy) {
                (b"AREN", true) => arena = Some(Self::locate_arena(&mut **r, size)?),
                (b"FPID", true) => ids = Some(Self::locate_ids(&mut **r, size)?),
                _ => {
                    let chunk = read_chunk_data(&mut *r, size)?;
                    match &tag {
                        b"FEND" => break,
                        b"POPC" => pop_count_offsets = parse_pop_counts(&chunk)?,
                        b"AREN" => arena = Some(Self::load_arena(chunk)?),
                        b"FPID" => ids = Some(Self::load_ids(chunk)?),
                        b"META" | b"HASH" => {}
                        other => crate::rd_log_warning!(
                            "Unknown chunk: {} ignored.\n",
                            String::from_utf8_lossy(other)
                        ),
                    }
                }
            }
        }
        let (header, len, fps) = arena.ok_or_else(|| bad_file("No AREN record found"))?;
        let (num_4byte_ids, num_8byte_ids, ids) =
            ids.ok_or_else(|| bad_file("No FPID record found"))?;
        Ok(Loaded {
            len,
            n_bits: header.n_bytes.saturating_mul(8),
            stored_bytes: header.stored_bytes as usize,
            pop_count_offsets,
            num_4byte_ids,
            num_8byte_ids,
            fps,
            ids,
        })
    }

    fn load_arena(chunk: Vec<u8>) -> Result<(ArenaHeader, u32, FpData)> {
        let head: &[u8; 9] = chunk
            .get(..9)
            .and_then(|h| h.try_into().ok())
            .ok_or_else(|| bad_file("Invalid AREN record"))?;
        let header = ArenaHeader::from_bytes(head);
        let len = header.len(chunk.len() as u64)?;
        let start = 9 + usize::from(header.spacer);
        Ok((header, len, FpData::Memory(chunk[start..].to_vec())))
    }

    /// Lazy mode: reads the `AREN` header, records where the fingerprints
    /// start and skips past them.
    fn locate_arena(r: &mut dyn FpbSource, size: u64) -> Result<(ArenaHeader, u32, FpData)> {
        let mut head = [0; 9];
        r.read_exact(&mut head).map_err(truncated)?;
        let header = ArenaHeader::from_bytes(&head);
        let len = header.len(size)?;
        r.seek(SeekFrom::Current(i64::from(header.spacer)))
            .map_err(truncated)?;
        let start = r.stream_position().map_err(truncated)?;
        r.seek(SeekFrom::Current(
            i64::from(header.stored_bytes) * i64::from(len),
        ))
        .map_err(truncated)?;
        Ok((header, len, FpData::Lazy(start)))
    }

    fn load_ids(chunk: Vec<u8>) -> Result<(u32, u32, IdData)> {
        let (n4, n8) = match (le_u32(&chunk, 0), le_u32(&chunk, 4)) {
            (Some(n4), Some(n8)) => (n4, n8),
            _ => return Err(bad_file("Invalid FPID record")),
        };
        let table = (u64::from(n4) + 1) * 4 + u64::from(n8) * 8;
        let offsets = (chunk.len() as u64)
            .checked_sub(table)
            .ok_or_else(|| bad_file("Invalid FPID record"))? as usize;
        Ok((n4, n8, IdData::Memory { chunk, offsets }))
    }

    /// Lazy mode: records where the `FPID` chunk and its offset table are.
    fn locate_ids(r: &mut dyn FpbSource, size: u64) -> Result<(u32, u32, IdData)> {
        let start = r.stream_position().map_err(truncated)?;
        let n4 = read_u32(r).map_err(truncated)?;
        let n8 = read_u32(r).map_err(truncated)?;
        let table = (u64::from(n4) + 1) * 4 + u64::from(n8) * 8;
        let offsets = (start + size)
            .checked_sub(table)
            .ok_or_else(|| bad_file("Invalid FPID record"))?;
        r.seek(SeekFrom::Start(start + size)).map_err(truncated)?;
        Ok((
            n4,
            n8,
            IdData::Lazy {
                chunk: start,
                offsets,
            },
        ))
    }

    fn loaded(&self) -> Result<&Loaded> {
        self.loaded
            .as_ref()
            .ok_or_else(|| RDKitError::precondition("not initialized"))
    }

    /// Number of fingerprints.
    pub fn len(&self) -> Result<u32> {
        Ok(self.loaded()?.len)
    }

    pub fn is_empty(&self) -> Result<bool> {
        Ok(self.len()? == 0)
    }

    /// Number of bits in each fingerprint.
    pub fn n_bits(&self) -> Result<u32> {
        Ok(self.loaded()?.n_bits)
    }

    /// Bytes a query fingerprint must provide.
    pub fn stored_bytes(&self) -> Result<usize> {
        Ok(self.loaded()?.stored_bytes)
    }

    /// `count` consecutive stored fingerprints starting at `which`.
    fn fp_bytes(&self, which: u32, count: u32) -> Result<std::borrow::Cow<'_, [u8]>> {
        let l = self.loaded()?;
        if u64::from(which) + u64::from(count) > u64::from(l.len) {
            return Err(RDKitError::Value("bad index".into()));
        }
        let offset = which as usize * l.stored_bytes;
        let size = count as usize * l.stored_bytes;
        match &l.fps {
            FpData::Memory(data) => data
                .get(offset..offset + size)
                .map(std::borrow::Cow::Borrowed)
                .ok_or_else(|| bad_file("truncated fingerprint data")),
            FpData::Lazy(start) => {
                let mut buf = vec![0; size];
                let mut r = self.stream();
                r.seek(SeekFrom::Start(start + offset as u64))
                    .and_then(|_| r.read_exact(&mut buf))
                    .map_err(|_| bad_file("truncated fingerprint data"))?;
                Ok(std::borrow::Cow::Owned(buf))
            }
        }
    }

    /// The raw stored bytes of fingerprint `idx`.
    pub fn get_bytes(&self, idx: u32) -> Result<Vec<u8>> {
        Ok(self.fp_bytes(idx, 1)?.into_owned())
    }

    /// Fingerprint `idx` as a bit vector.
    pub fn get_fp(&self, idx: u32) -> Result<ExplicitBitVect> {
        let l = self.loaded()?;
        let bytes = self.fp_bytes(idx, 1)?;
        let mut bv = ExplicitBitVect::new(l.n_bits, false);
        for bit in 0..l.n_bits {
            if bytes[(bit / 8) as usize] >> (bit % 8) & 1 == 1 {
                bv.set_bit(bit)?;
            }
        }
        Ok(bv)
    }

    /// Identifier of fingerprint `idx`.
    pub fn get_id(&self, idx: u32) -> Result<String> {
        let l = self.loaded()?;
        if u64::from(idx) >= u64::from(l.num_4byte_ids) + u64::from(l.num_8byte_ids) {
            return Err(RDKitError::Value("bad index".into()));
        }
        // Offsets of the identifier within the chunk: where it starts and
        // where the next one starts. Upstream's addressing of the 8 byte
        // entries is reproduced as is.
        let i = u64::from(idx);
        let n4 = u64::from(l.num_4byte_ids);
        let (start_at, start_width, end_at, end_width) = if i < n4 {
            (i * 4, 4, (i + 1) * 4, 4)
        } else if i == n4 {
            (i * 4, 4, (i + 1) * 4, 8)
        } else {
            (n4 * 4 + i * 8, 8, n4 * 4 + (i + 1) * 8, 8)
        };
        let table_entry = |at: u64, width: usize| -> Result<u64> {
            let mut buf = [0u8; 8];
            match &l.ids {
                IdData::Memory { chunk, offsets } => {
                    let from = *offsets + at as usize;
                    let src = chunk
                        .get(from..from + width)
                        .ok_or_else(|| bad_file("truncated FPID record"))?;
                    buf[..width].copy_from_slice(src);
                }
                IdData::Lazy { offsets, .. } => {
                    let mut r = self.stream();
                    r.seek(SeekFrom::Start(offsets + at))
                        .and_then(|_| r.read_exact(&mut buf[..width]))
                        .map_err(|_| bad_file("truncated FPID record"))?;
                }
            }
            Ok(u64::from_le_bytes(buf))
        };
        let start = table_entry(start_at, start_width)?;
        let end = table_entry(end_at, end_width)?;
        let len = end.wrapping_sub(start) as usize;
        let raw = match &l.ids {
            IdData::Memory { chunk, .. } => chunk
                .get(start as usize..start as usize + len)
                .map(<[u8]>::to_vec)
                .ok_or_else(|| bad_file("truncated FPID record"))?,
            IdData::Lazy { chunk, .. } => {
                let mut buf = vec![0; len];
                let mut r = self.stream();
                r.seek(SeekFrom::Start(chunk + start))
                    .and_then(|_| r.read_exact(&mut buf))
                    .map_err(|_| bad_file("truncated FPID record"))?;
                buf
            }
        };
        // identifiers are NUL terminated
        let text = raw.split(|&b| b == 0).next().unwrap_or_default();
        Ok(String::from_utf8_lossy(text).into_owned())
    }

    /// Fingerprint `idx` and its identifier.
    pub fn get(&self, idx: u32) -> Result<(ExplicitBitVect, String)> {
        Ok((self.get_fp(idx)?, self.get_id(idx)?))
    }

    /// Range of fingerprint indices whose population count lies in
    /// `min_count..=max_count`; `(0, 0)` when the file has no `POPC` chunk.
    pub fn get_fp_ids_in_count_range(&self, min_count: u32, max_count: u32) -> Result<(u32, u32)> {
        let l = self.loaded()?;
        if max_count > l.n_bits {
            return Err(RDKitError::range("maxCount", max_count, l.n_bits + 1));
        }
        if max_count < min_count {
            return Err(RDKitError::precondition("max < min"));
        }
        if l.has_pop_counts() {
            Ok((
                l.pop_count_offsets[min_count as usize],
                l.pop_count_offsets[max_count as usize + 1],
            ))
        } else {
            Ok((0, 0))
        }
    }

    fn check_query<'q>(&self, query: &'q [u8]) -> Result<&'q [u8]> {
        let stored = self.loaded()?.stored_bytes;
        query.get(..stored).ok_or_else(|| {
            RDKitError::Value(format!(
                "the query must provide at least {stored} bytes, got {}",
                query.len()
            ))
        })
    }

    fn check_index(&self, idx: u32) -> Result<()> {
        if idx >= self.loaded()?.len {
            return Err(RDKitError::Value("bad index".into()));
        }
        Ok(())
    }

    /// Tanimoto similarity of fingerprint `idx` to `query`.
    pub fn get_tanimoto(&self, idx: u32, query: &[u8]) -> Result<f64> {
        let query = self.check_query(query)?;
        self.check_index(idx)?;
        Ok(bitmap_tanimoto(&self.fp_bytes(idx, 1)?, query))
    }

    /// Tversky similarity of fingerprint `idx` to `query`.
    pub fn get_tversky(&self, idx: u32, query: &[u8], ca: f64, cb: f64) -> Result<f64> {
        let query = self.check_query(query)?;
        self.check_index(idx)?;
        Ok(bitmap_tversky(&self.fp_bytes(idx, 1)?, query, ca, cb))
    }

    /// Scans fingerprints `start..end` in blocks, keeping those for which
    /// `visit` returns a value.
    fn scan<T>(
        &self,
        start: u32,
        end: u32,
        mut visit: impl FnMut(u32, &[u8]) -> Option<T>,
    ) -> Result<Vec<T>> {
        let stored = self.loaded()?.stored_bytes;
        let mut hits = Vec::new();
        let mut i = start;
        while i < end {
            let count = (end - i).min(READ_CACHE as u32);
            let block = self.fp_bytes(i, count)?;
            hits.extend(
                block
                    .chunks_exact(stored)
                    .zip(i..)
                    .filter_map(|(fp, idx)| visit(idx, fp)),
            );
            i += count;
        }
        Ok(hits)
    }

    /// Neighbors with Tanimoto similarity of at least `threshold`, most
    /// similar first (ties by index).
    ///
    /// With `use_popcount_screen` the population counts stored in the file
    /// bound the fingerprints that need to be compared (Swamidass & Baldi,
    /// J. Chem. Inf. Model. 47, 302 (2007), equation 24).
    pub fn get_tanimoto_neighbors(
        &self,
        query: &[u8],
        threshold: f64,
        use_popcount_screen: bool,
    ) -> Result<Vec<Neighbor>> {
        let query = self.check_query(query)?;
        check_threshold(threshold)?;
        let l = self.loaded()?;
        let probe = f64::from(bitmap_popcount(query));
        let (start, end) = l.scan_bounds(use_popcount_screen, || {
            let min = (threshold * probe).floor() as u32;
            let max = (threshold > 1e-6).then(|| (probe / threshold).ceil() as u32);
            (min, max)
        });
        let mut hits = self.scan(start, end, |idx, fp| {
            let sim = bitmap_tanimoto(fp, query);
            (sim >= threshold).then_some((sim, idx))
        })?;
        sort_by_similarity(&mut hits);
        Ok(hits)
    }

    /// Neighbors with Tversky similarity of at least `threshold`, most
    /// similar first (equation 25 of the paper cited at
    /// [`FpbReader::get_tanimoto_neighbors`] bounds the screen).
    pub fn get_tversky_neighbors(
        &self,
        query: &[u8],
        ca: f64,
        cb: f64,
        threshold: f64,
        use_popcount_screen: bool,
    ) -> Result<Vec<Neighbor>> {
        let query = self.check_query(query)?;
        check_threshold(threshold)?;
        let l = self.loaded()?;
        let probe = f64::from(bitmap_popcount(query));
        let (start, end) = l.scan_bounds(use_popcount_screen, || {
            let min =
                ((threshold * probe * ca) / (1.0 - threshold + threshold * ca)).floor() as u32;
            let max = (threshold * cb > 1e-6).then(|| {
                (probe * (1.0 - threshold + threshold * cb) / (threshold * cb)).ceil() as u32
            });
            (min, max)
        });
        let mut hits = self.scan(start, end, |idx, fp| {
            let sim = bitmap_tversky(fp, query, ca, cb);
            (sim >= threshold).then_some((sim, idx))
        })?;
        sort_by_similarity(&mut hits);
        Ok(hits)
    }

    /// Indices of the fingerprints that have every bit of `query` set.
    pub fn get_containing_neighbors(&self, query: &[u8]) -> Result<Vec<u32>> {
        let query = self.check_query(query)?;
        let l = self.loaded()?;
        let start = if l.has_pop_counts() {
            l.pop_count_offsets
                .get(bitmap_popcount(query) as usize)
                .copied()
                .unwrap_or(l.len)
        } else {
            0
        };
        let mut hits = self.scan(start, l.len, |idx, fp| {
            bitmap_all_probe_bits_match(query, fp).then_some(idx)
        })?;
        hits.sort_unstable();
        Ok(hits)
    }
}

impl Loaded {
    fn has_pop_counts(&self) -> bool {
        self.pop_count_offsets.len() == self.n_bits as usize + 2
    }

    /// The range of fingerprints to compare. `bounds` yields the lowest and
    /// highest population counts a hit can have (`None`: no upper limit).
    fn scan_bounds(
        &self,
        use_popcount_screen: bool,
        bounds: impl FnOnce() -> (u32, Option<u32>),
    ) -> (u32, u32) {
        if !(use_popcount_screen && self.has_pop_counts()) {
            return (0, self.len);
        }
        let (min, max) = bounds();
        let max = max.unwrap_or(self.stored_bytes as u32);
        // the offsets table has one entry per possible population count plus
        // an end marker; clamp counts that fall beyond it
        let last = self.pop_count_offsets.len() - 1;
        let at = |count: usize| self.pop_count_offsets[count.min(last)];
        (at(min as usize), at(max as usize + 1))
    }
}

fn check_threshold(threshold: f64) -> Result<()> {
    if (-1e-6..=1.0 + 1e-6).contains(&threshold) {
        Ok(())
    } else {
        Err(RDKitError::Invariant {
            kind: "Range Error",
            message: format!("threshold\n\tFailed Expression: -1e-06 <= {threshold} <= 1.000001"),
        })
    }
}

/// Most similar first. Upstream sorts with an unstable `std::sort` on the
/// similarity alone; ties are ordered by index here to be deterministic.
fn sort_by_similarity(hits: &mut [Neighbor]) {
    hits.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
}
