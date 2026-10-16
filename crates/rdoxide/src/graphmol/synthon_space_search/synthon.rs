//! `Synthon.cpp`: a building block of a SynthonSpace.

use super::details;
use super::shape_input::SynthonShapeInput;
use super::stream::{Reader, Writer};
use crate::datastructs::bit_vect::{BitVect, ExplicitBitVect, bit_vect_from_binary};
use crate::graphmol::chirality::find_stereo::find_potential_stereo;
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::fingerprints::pattern::pattern_fingerprint_mol;
use crate::graphmol::mol_ops::frags::get_mol_frags_as_mols;
use crate::graphmol::mol_props::get_avg_mol_wt;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::rdgeneral::RDProps;
use crate::{RDKitError, Result};

pub const PATT_FP_NUM_BITS: u32 = 1024;

#[derive(Clone)]
pub struct Synthon {
    smiles: String,
    // Two copies of the molecule: the original, used for building products,
    // and the search molecule whose atoms and bonds have been fiddled with to
    // match the product (the aliphatic precursor to aromatic product issue).
    orig_mol: ROMol,
    search_mol: Option<ROMol>,
    patt_fp: Option<ExplicitBitVect>,
    fp: Option<ExplicitBitVect>,
    shapes: Option<SynthonShapeInput>,
    conn_regions: Vec<ROMol>,
    num_dummies: u32,
    num_heavy_atoms: u32,
    num_chiral_atoms: u32,
    num_chiral_atoms_exc_dummies: u32,
    mol_wt: f64,
    max_synthon_set_size: u32,
}

impl Synthon {
    pub fn new(smi: &str) -> Result<Self> {
        let params = SmilesParserParams {
            sanitize: false,
            ..Default::default()
        };
        let orig_mol = mol_from_smiles(smi, &params)?
            .ok_or_else(|| RDKitError::Value(format!("Unparsable synthon SMILES {smi}")))?;
        let mut s = Synthon::blank(smi.to_string(), orig_mol);
        s.calc_properties()?;
        Ok(s)
    }

    fn blank(smiles: String, orig_mol: ROMol) -> Self {
        Synthon {
            smiles,
            orig_mol,
            search_mol: None,
            patt_fp: None,
            fp: None,
            shapes: None,
            conn_regions: Vec::new(),
            num_dummies: 0,
            num_heavy_atoms: 0,
            num_chiral_atoms: 0,
            num_chiral_atoms_exc_dummies: 0,
            mol_wt: 0.0,
            max_synthon_set_size: 0,
        }
    }

    pub fn smiles(&self) -> &str {
        &self.smiles
    }
    pub fn orig_mol(&self) -> &ROMol {
        &self.orig_mol
    }
    pub fn search_mol(&self) -> &ROMol {
        self.search_mol.as_ref().expect("synthon has a search mol")
    }
    pub fn patt_fp(&self) -> &ExplicitBitVect {
        self.patt_fp.as_ref().expect("synthon has a pattern FP")
    }
    pub fn fp(&self) -> Option<&ExplicitBitVect> {
        self.fp.as_ref()
    }
    pub fn conn_regions(&self) -> &[ROMol] {
        &self.conn_regions
    }
    pub fn num_dummies(&self) -> u32 {
        self.num_dummies
    }
    pub fn num_heavy_atoms(&self) -> u32 {
        self.num_heavy_atoms
    }
    pub fn num_chiral_atoms(&self) -> u32 {
        self.num_chiral_atoms
    }
    pub fn num_chiral_atoms_exc_dummies(&self) -> u32 {
        self.num_chiral_atoms_exc_dummies
    }
    pub fn mol_wt(&self) -> f64 {
        self.mol_wt
    }
    pub fn max_synthon_set_size(&self) -> u32 {
        self.max_synthon_set_size
    }
    pub fn update_max_synthon_set_size(&mut self, new_val: u32) {
        if new_val > self.max_synthon_set_size {
            self.max_synthon_set_size = new_val;
        }
    }
    pub fn set_fp(&mut self, fp: ExplicitBitVect) {
        self.fp = Some(fp);
    }
    pub fn has_shapes(&self) -> bool {
        self.shapes.is_some()
    }
    pub fn shapes(&self) -> Option<&SynthonShapeInput> {
        self.shapes.as_ref()
    }
    pub fn set_shapes(&mut self, shapes: SynthonShapeInput) {
        self.shapes = Some(shapes);
    }
    pub fn clear_shapes(&mut self) {
        self.shapes = None;
    }

    pub fn set_search_mol(&mut self, mut mol: ROMol) -> Result<()> {
        // There are probably extraneous props on the atoms and bonds
        for i in 0..mol.num_atoms() {
            let a = mol.atom_mut(i);
            a.clear_prop("molNum");
            a.clear_prop("idx");
        }
        for i in 0..mol.num_bonds() {
            let b = mol.bond_mut(i);
            b.clear_prop("molNum");
            b.clear_prop("idx");
        }
        self.search_mol = Some(mol);
        self.finish_initialization()
    }

    fn finish_initialization(&mut self) -> Result<()> {
        let sm = self.search_mol.as_ref().expect("search mol");
        self.patt_fp = Some(pattern_fingerprint_mol(sm, PATT_FP_NUM_BITS, false)?);
        self.conn_regions.clear();
        if let Some(cr) = details::build_conn_region(sm)? {
            let (frags, _) = get_mol_frags_as_mols(&cr, false, true)?;
            self.conn_regions = frags;
        }
        Ok(())
    }

    fn calc_properties(&mut self) -> Result<()> {
        self.num_dummies = 0;
        self.num_heavy_atoms = 0;
        assign_stereochemistry(&mut self.orig_mol, false, false, false)?;
        for atom in self.orig_mol.atoms() {
            if atom.atomic_num() == 0 {
                self.num_dummies += 1;
            } else if atom.atomic_num() > 1 {
                self.num_heavy_atoms += 1;
            }
        }
        find_potential_stereo(&mut self.orig_mol, true, true)?;
        let (n, exc) = details::count_chiral_atoms(&mut self.orig_mol)?;
        self.num_chiral_atoms = n;
        self.num_chiral_atoms_exc_dummies = exc;
        self.mol_wt = get_avg_mol_wt(&self.orig_mol, false)?;
        Ok(())
    }

    pub fn write_to_db_stream(&self, w: &mut Writer) -> Result<()> {
        w.string(self.smiles.as_bytes());
        w.mol(&self.orig_mol)?;
        w.mol(self.search_mol())?;
        w.string(&self.patt_fp().to_binary());
        w.u64(self.conn_regions.len() as u64);
        for cr in &self.conn_regions {
            w.mol(cr)?;
        }
        if let Some(fp) = &self.fp {
            w.bool(true);
            w.string(&fp.to_binary());
        } else {
            w.bool(false);
        }
        w.u32(self.num_dummies);
        w.u32(self.num_heavy_atoms);
        w.u32(self.num_chiral_atoms);
        w.f64(self.mol_wt);
        if let Some(sh) = &self.shapes {
            w.bool(true);
            w.string(sh.to_text().as_bytes());
        } else {
            w.bool(false);
        }
        Ok(())
    }

    pub fn read_from_db_stream(r: &mut Reader<'_>, version: u32) -> Result<Self> {
        let smiles = r.string()?;
        let orig_mol = r.mol()?;
        let mut s = Synthon::blank(smiles, orig_mol);
        s.search_mol = Some(r.mol()?);
        s.patt_fp = Some(bit_vect_from_binary(&r.bytes()?, false)?);
        let n = r.u64()?;
        for _ in 0..n {
            s.conn_regions.push(r.mol()?);
        }
        if r.bool()? {
            s.fp = Some(bit_vect_from_binary(&r.bytes()?, false)?);
        }
        if version > 3000 {
            s.num_dummies = r.u32()?;
            s.num_heavy_atoms = r.u32()?;
            s.num_chiral_atoms = r.u32()?;
            s.mol_wt = r.f64()?;
        } else {
            s.calc_properties()?;
        }
        if version > 3010 && r.bool()? {
            s.shapes = Some(SynthonShapeInput::from_text(&r.string()?)?);
        }
        Ok(s)
    }
}

impl std::fmt::Debug for Synthon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Synthon({})", self.smiles)
    }
}

// keep the trait in scope for `.num_bits()`/`.num_on_bits()` users
#[allow(unused_imports)]
use BitVect as _;
