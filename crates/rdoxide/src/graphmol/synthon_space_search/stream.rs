//! Little-endian binary stream helpers matching upstream's `streamWrite` and
//! `streamRead` (`RDGeneral/StreamOps.h`).

use crate::graphmol::pickle::{flags, mol_from_pickle_stream, pickle_mol};
use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result};

#[derive(Default)]
pub struct Writer {
    pub buf: Vec<u8>,
}

impl Writer {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn bool(&mut self, v: bool) {
        self.buf.push(u8::from(v));
    }
    pub fn i32(&mut self, v: i32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }
    pub fn u32(&mut self, v: u32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }
    pub fn u64(&mut self, v: u64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }
    pub fn f64(&mut self, v: f64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }
    /// a `std::string`: a 32-bit length and the bytes
    pub fn string(&mut self, s: &[u8]) {
        self.u32(s.len() as u32);
        self.buf.extend_from_slice(s);
    }
    pub fn mol(&mut self, m: &ROMol) -> Result<()> {
        self.buf
            .extend_from_slice(&pickle_mol(m, flags::ALL_PROPS)?);
        Ok(())
    }
    pub fn bits(&mut self, bits: u32, size: usize) {
        self.u64(size as u64);
        for i in 0..size {
            self.bool(bits & (1 << i) != 0);
        }
    }
}

pub struct Reader<'a> {
    pub data: &'a [u8],
    pub pos: usize,
}

fn eof() -> RDKitError {
    RDKitError::Value("Unexpected end of SynthonSpace data".into())
}

impl<'a> Reader<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Reader { data, pos: 0 }
    }
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self.pos.checked_add(n).ok_or_else(eof)?;
        let s = self.data.get(self.pos..end).ok_or_else(eof)?;
        self.pos = end;
        Ok(s)
    }
    pub fn bool(&mut self) -> Result<bool> {
        Ok(self.take(1)?[0] != 0)
    }
    pub fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    pub fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    pub fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    pub fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    pub fn bytes(&mut self) -> Result<Vec<u8>> {
        let n = self.u32()? as usize;
        Ok(self.take(n)?.to_vec())
    }
    pub fn string(&mut self) -> Result<String> {
        Ok(String::from_utf8_lossy(&self.bytes()?).into_owned())
    }
    pub fn mol(&mut self) -> Result<ROMol> {
        mol_from_pickle_stream(self.data, &mut self.pos, flags::ALL_PROPS)
    }
    pub fn bits(&mut self) -> Result<u32> {
        let size = self.u64()? as usize;
        let mut res = 0u32;
        for i in 0..size {
            if self.bool()? && i < 32 {
                res |= 1 << i;
            }
        }
        Ok(res)
    }
}
