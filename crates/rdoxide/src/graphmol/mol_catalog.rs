//! Port of RDKit `Code/GraphMol/MolCatalog`: catalogs of molecules.

use crate::catalogs::{CatalogEntry, CatalogParams, HierarchCatalog, StreamReader, StreamWriter};
use crate::graphmol::pickle::{get_default_pickle_properties, mol_from_pickle_stream, pickle_mol};
use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result};

/// `MolCatalogParams`: nothing to store.
#[derive(Debug, Clone)]
pub struct MolCatalogParams {
    type_str: String,
}

impl Default for MolCatalogParams {
    fn default() -> Self {
        MolCatalogParams {
            type_str: "MolCatalog Parameters".into(),
        }
    }
}

impl CatalogParams for MolCatalogParams {
    fn type_str(&self) -> &str {
        &self.type_str
    }
    fn to_stream(&self, _w: &mut StreamWriter) {}
    fn init_from_stream(&mut self, _r: &mut StreamReader) -> Result<()> {
        Ok(())
    }
}

/// `MolCatalogEntry`
#[derive(Debug, Clone)]
pub struct MolCatalogEntry {
    mol: Option<ROMol>,
    bit_id: i32,
    order: u32,
    description: String,
}

impl Default for MolCatalogEntry {
    fn default() -> Self {
        MolCatalogEntry {
            mol: None,
            bit_id: -1,
            order: 0,
            description: String::new(),
        }
    }
}

impl MolCatalogEntry {
    pub fn new(mol: ROMol) -> Self {
        MolCatalogEntry {
            mol: Some(mol),
            ..Default::default()
        }
    }

    pub fn mol(&self) -> Option<&ROMol> {
        self.mol.as_ref()
    }

    pub fn set_mol(&mut self, mol: ROMol) {
        self.mol = Some(mol);
    }

    pub fn set_description(&mut self, val: &str) {
        self.description = val.to_string();
    }

    pub fn set_order(&mut self, order: u32) {
        self.order = order;
    }

    /// `Serialize`
    pub fn serialize(&self) -> Result<Vec<u8>> {
        let mut w = StreamWriter::default();
        self.to_stream(&mut w)?;
        Ok(w.0)
    }

    /// `MolCatalogEntry(pickle)`
    pub fn from_pickle(data: &[u8]) -> Result<Self> {
        Self::from_stream(&mut StreamReader::new(data))
    }
}

impl CatalogEntry for MolCatalogEntry {
    type Order = i32;

    fn bit_id(&self) -> i32 {
        self.bit_id
    }
    fn set_bit_id(&mut self, id: i32) {
        self.bit_id = id;
    }
    fn description(&self) -> String {
        self.description.clone()
    }
    fn order(&self) -> i32 {
        self.order as i32
    }

    fn to_stream(&self, w: &mut StreamWriter) -> Result<()> {
        let mol = self
            .mol
            .as_ref()
            .ok_or_else(|| RDKitError::precondition("bad mol"))?;
        w.bytes(&pickle_mol(mol, get_default_pickle_properties())?);
        w.i32(self.bit_id);
        w.i32(self.order as i32);
        w.string(&self.description);
        Ok(())
    }

    fn from_stream(r: &mut StreamReader) -> Result<Self> {
        let mol = mol_from_pickle_stream(r.data, &mut r.pos, get_default_pickle_properties())?;
        let bit_id = r.i32()?;
        let order = r.i32()? as u32;
        let description = r.string()?;
        Ok(MolCatalogEntry {
            mol: Some(mol),
            bit_id,
            order,
            description,
        })
    }
}

/// `MolCatalog`
pub type MolCatalog = HierarchCatalog<MolCatalogEntry, MolCatalogParams>;
