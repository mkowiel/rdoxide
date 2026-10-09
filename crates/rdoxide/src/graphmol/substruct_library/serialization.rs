//! Serialization of a library. Upstream pickles through a boost text archive;
//! this port uses a compact binary format of its own, so pickles are not
//! exchangeable with RDKit's.

use std::sync::Arc;

use super::{KeyFromPropHolder, MolHolder, PatternHolder, SubstructLibrary, shared};
use crate::datastructs::bit_vect::{base64_decode, base64_encode, bit_vect_from_binary};
use crate::datastructs::{BitVect, ExplicitBitVect};
use crate::error::{RDKitError, Result};
use crate::graphmol::pickle::{get_default_pickle_properties, mol_from_pickle, pickle_mol};

const MAGIC: &[u8; 8] = b"RDOXSSL1";

#[derive(Default)]
struct Writer(Vec<u8>);

impl Writer {
    fn u8(&mut self, v: u8) {
        self.0.push(v);
    }
    fn u32(&mut self, v: u32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn bytes(&mut self, v: &[u8]) {
        self.u32(v.len() as u32);
        self.0.extend_from_slice(v);
    }
    fn strings(&mut self, v: &[String]) {
        self.u32(v.len() as u32);
        v.iter().for_each(|s| self.bytes(s.as_bytes()));
    }
}

struct Reader<'a> {
    data: &'a [u8],
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.data.len() < n {
            return Err(RDKitError::Runtime(
                "truncated SubstructLibrary data".into(),
            ));
        }
        let (head, tail) = self.data.split_at(n);
        self.data = tail;
        Ok(head)
    }
    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(
            self.take(4)?.try_into().expect("4 bytes"),
        ))
    }
    fn bytes(&mut self) -> Result<&'a [u8]> {
        let n = self.u32()? as usize;
        self.take(n)
    }
    fn strings(&mut self) -> Result<Vec<String>> {
        (0..self.u32()?)
            .map(|_| Ok(String::from_utf8_lossy(self.bytes()?).into_owned()))
            .collect()
    }
}

impl SubstructLibrary {
    /// Serializes the library with its holders and search order.
    pub fn to_bytes(&self) -> Result<Vec<u8>> {
        let mut w = Writer::default();
        w.0.extend_from_slice(MAGIC);
        match &*super::read(&self.mols) {
            MolHolder::Mols(v) => {
                w.u8(0);
                w.u32(v.len() as u32);
                for m in v {
                    w.bytes(&pickle_mol(m, get_default_pickle_properties())?);
                }
            }
            MolHolder::Binary(v) => {
                w.u8(1);
                w.u32(v.len() as u32);
                v.iter().for_each(|p| w.bytes(p));
            }
            MolHolder::Smiles(v) => {
                w.u8(2);
                w.strings(v);
            }
            MolHolder::TrustedSmiles(v) => {
                w.u8(3);
                w.strings(v);
            }
        }
        match &self.fps {
            Some(fps) => {
                let fps = super::read(fps);
                w.u8(1 + u8::from(fps.is_tautomeric()));
                w.u32(fps.num_bits());
                w.u32(fps.len());
                fps.fingerprints()
                    .iter()
                    .for_each(|fp| w.bytes(&fp.to_binary()));
            }
            None => w.u8(0),
        }
        match &self.keys {
            Some(keys) => {
                let keys = super::read(keys);
                w.u8(1);
                w.bytes(keys.prop_name().as_bytes());
                w.strings(keys.keys());
            }
            None => w.u8(0),
        }
        w.u32(self.search_order.len() as u32);
        self.search_order.iter().for_each(|&i| w.u32(i));
        Ok(w.0)
    }

    /// [`to_bytes`](Self::to_bytes) as base64 text, which survives text streams.
    pub fn to_text(&self) -> Result<String> {
        Ok(base64_encode(&self.to_bytes()?))
    }

    /// The inverse of [`to_text`](Self::to_text).
    pub fn from_text(text: &str) -> Result<Self> {
        Self::from_bytes(&base64_decode(text.trim().as_bytes()))
    }

    /// The inverse of [`to_bytes`](Self::to_bytes).
    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        let mut r = Reader { data };
        if r.take(MAGIC.len())? != MAGIC {
            return Err(RDKitError::Runtime("not a SubstructLibrary pickle".into()));
        }
        let mols = match r.u8()? {
            0 => MolHolder::Mols(
                (0..r.u32()?)
                    .map(|_| {
                        mol_from_pickle(r.bytes()?, get_default_pickle_properties()).map(Arc::new)
                    })
                    .collect::<Result<_>>()?,
            ),
            1 => MolHolder::Binary(
                (0..r.u32()?)
                    .map(|_| r.bytes().map(<[u8]>::to_vec))
                    .collect::<Result<_>>()?,
            ),
            2 => MolHolder::Smiles(r.strings()?),
            3 => MolHolder::TrustedSmiles(r.strings()?),
            k => return Err(RDKitError::Runtime(format!("bad molecule holder kind {k}"))),
        };
        let fps = match r.u8()? {
            0 => None,
            kind => {
                let num_bits = r.u32()?;
                let mut holder = if kind == 2 {
                    PatternHolder::new_tautomeric(num_bits)
                } else {
                    PatternHolder::new(num_bits)
                };
                for _ in 0..r.u32()? {
                    holder.add_fingerprint(bit_vect_from_binary::<ExplicitBitVect>(
                        r.bytes()?,
                        false,
                    )?);
                }
                Some(shared(holder))
            }
        };
        let keys = match r.u8()? {
            0 => None,
            _ => {
                let mut holder = KeyFromPropHolder::new(&String::from_utf8_lossy(r.bytes()?));
                *holder.keys_mut() = r.strings()?;
                Some(shared(holder))
            }
        };
        let order = (0..r.u32()?).map(|_| r.u32()).collect::<Result<Vec<_>>>()?;
        let mut lib = SubstructLibrary::new(shared(mols), fps, keys);
        lib.set_search_order(order)?;
        Ok(lib)
    }
}
