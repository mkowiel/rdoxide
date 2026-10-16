//! `SynthonSet.cpp`: the synthons of one reaction.

use std::collections::BTreeSet;
use std::fmt::Write as _;

use super::details::{self, ConnectorBits, Stepper};
use super::params::MAX_CONNECTOR_NUM;
use super::stream::{Reader, Writer};
use super::synthon::{PATT_FP_NUM_BITS, Synthon};
use crate::datastructs::bit_vect::{BitVect, ExplicitBitVect, bit_vect_from_binary};
use crate::geometry::Point3D;
use crate::graphmol::chem_transforms::combine_mols;
use crate::graphmol::chem_transforms::mol_fragmenter::fragment_on_bonds;
use crate::graphmol::chem_transforms::mol_zip::{MolzipLabel, MolzipParams, molzip_single};
use crate::graphmol::find_rings::SymmetrizeSSSRAlgorithm;
use crate::graphmol::find_rings::symmetrize_sssr;
use crate::graphmol::fingerprints::pattern::pattern_fingerprint_mol;
use crate::graphmol::mol_ops::frags::get_mol_frags_as_mols;
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// A synthon and the position of a sample molecule built from it.
pub struct SampleMolRec {
    pub synthon_idx: usize,
    pub synthon_nums: Vec<usize>,
    pub synthon_set_num: usize,
    pub num_atoms: u32,
}

#[derive(Default)]
pub struct SynthonSet {
    id: String,
    /// (synthon id, index into the SynthonSpace's synthon pool)
    synthons: Vec<Vec<(String, usize)>>,
    substructure_search_orders: Vec<Vec<usize>>,
    fingerprint_search_orders: Vec<Vec<usize>>,
    rascal_search_orders: Vec<Vec<usize>>,
    connectors: ConnectorBits,
    synth_conn_patts: Vec<ConnectorBits>,
    connector_regions: Vec<ROMol>,
    conn_reg_smis: Vec<String>,
    conn_reg_fps: Vec<ExplicitBitVect>,
    add_fp: Option<ExplicitBitVect>,
    subtract_fp: Option<ExplicitBitVect>,
    num_connectors: Vec<i32>,
    num_ring_formers: u32,
}

fn smiles_string(mol: &ROMol) -> String {
    details::smiles_of(mol)
}

impl SynthonSet {
    pub fn new(id: &str) -> Self {
        SynthonSet {
            id: id.to_string(),
            ..Default::default()
        }
    }

    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn synthons(&self) -> &[Vec<(String, usize)>] {
        &self.synthons
    }
    pub fn connectors(&self) -> ConnectorBits {
        self.connectors
    }
    pub fn synthon_connector_patterns(&self) -> &[ConnectorBits] {
        &self.synth_conn_patts
    }
    pub fn connector_regions(&self) -> &[ROMol] {
        &self.connector_regions
    }
    pub fn connector_region_smiles(&self) -> &[String] {
        &self.conn_reg_smis
    }
    pub fn conn_reg_fps(&self) -> &[ExplicitBitVect] {
        &self.conn_reg_fps
    }
    pub fn add_fp(&self) -> Option<&ExplicitBitVect> {
        self.add_fp.as_ref()
    }
    pub fn subtract_fp(&self) -> Option<&ExplicitBitVect> {
        self.subtract_fp.as_ref()
    }
    pub fn num_connectors(&self) -> &[i32] {
        &self.num_connectors
    }
    pub fn num_ring_formers(&self) -> u32 {
        self.num_ring_formers
    }
    pub fn has_add_and_subtract_fps(&self) -> bool {
        self.add_fp.is_some()
    }

    fn ordered(order: &[Vec<usize>], n_sets: usize, set_num: usize, n: usize) -> Option<usize> {
        if set_num >= n_sets || n >= order.get(set_num).map_or(0, |o| o.len()) {
            return None;
        }
        Some(order[set_num][n])
    }

    /// The number (in input order) of the `n`th synthon of set `set_num` in
    /// the substructure search order.
    pub fn substructure_ordered_synthon_num(&self, set_num: usize, n: usize) -> Option<usize> {
        Self::ordered(
            &self.substructure_search_orders,
            self.synthons.len(),
            set_num,
            n,
        )
    }
    pub fn fingerprint_ordered_synthon_num(&self, set_num: usize, n: usize) -> Option<usize> {
        Self::ordered(
            &self.fingerprint_search_orders,
            self.synthons.len(),
            set_num,
            n,
        )
    }
    pub fn rascal_ordered_synthon_num(&self, set_num: usize, n: usize) -> Option<usize> {
        Self::ordered(&self.rascal_search_orders, self.synthons.len(), set_num, n)
    }

    pub fn num_products(&self) -> u64 {
        self.synthons.iter().map(|r| r.len() as u64).product()
    }

    pub fn has_fingerprints(&self, pool: &[Synthon]) -> bool {
        self.synthons
            .first()
            .and_then(|s| s.first())
            .is_some_and(|(_, idx)| pool[*idx].fp().is_some())
    }

    pub fn add_synthon(&mut self, set_num: usize, synthon_idx: usize, synthon_id: &str) {
        if set_num >= self.synthons.len() {
            self.synthons.resize(set_num + 1, Vec::new());
        }
        self.synthons[set_num].push((synthon_id.to_string(), synthon_idx));
    }

    pub fn remove_empty_synthon_sets(&mut self) {
        self.synthons.retain(|r| !r.is_empty());
    }

    fn copy_synthon(&self, pool: &[Synthon], set_num: usize, idx: usize) -> Result<RWMol> {
        let mut m = pool[self.synthons[set_num][idx].1].orig_mol().clone();
        for i in 0..m.num_atoms() {
            let a = m.atom_mut(i);
            a.set_prop("molNum", set_num as u32, false)?;
            a.set_prop("idx", i, false)?;
        }
        for i in 0..m.num_bonds() {
            let b = m.bond_mut(i);
            b.set_prop("molNum", set_num as u32, false)?;
            b.set_prop("idx", i, false)?;
        }
        Ok(m)
    }

    fn copy_synthons(&self, pool: &[Synthon]) -> Result<Vec<Vec<RWMol>>> {
        (0..self.synthons.len())
            .map(|i| {
                (0..self.synthons[i].len())
                    .map(|j| self.copy_synthon(pool, i, j))
                    .collect()
            })
            .collect()
    }

    /// Make a molecule from the given synthon numbers (indices into
    /// `synthons`), built from the copies of the synthon molecules.
    pub fn build_molecule(&self, pool: &[Synthon], synthon_nums: &[usize]) -> Result<ROMol> {
        let mut comb_mol = ROMol::new();
        let mut mol_name = String::new();
        for (i, &sn) in synthon_nums.iter().enumerate() {
            let mut synth_mol = self.copy_synthon(pool, i, sn)?;
            for idx in 0..synth_mol.num_atoms() {
                let (is_conn, iso) = {
                    let a = synth_mol.atom(idx);
                    (a.atomic_num() == 0 && a.isotope() != 0, a.isotope())
                };
                if is_conn {
                    let nbrs: Vec<u32> = synth_mol.atom_neighbors(idx).collect();
                    for nbr in nbrs {
                        add_join_info(&mut synth_mol, nbr, iso)?;
                        add_join_info(&mut synth_mol, idx, iso)?;
                    }
                }
            }
            comb_mol = combine_mols(&comb_mol, &synth_mol, Point3D::default())?;
            if i > 0 {
                mol_name.push(';');
            }
            mol_name.push_str(&self.synthons[i][sn].0);
        }
        let mz = MolzipParams {
            label: MolzipLabel::Isotope,
            ..Default::default()
        };
        comb_mol.set_prop(cp::_Name, mol_name, false)?;
        let res = (|| -> Result<ROMol> {
            let mut zip = molzip_single(&comb_mol, &mz)?
                .ok_or_else(|| RDKitError::Value("molzip failed".into()))?;
            sanitize_mol(&mut zip, SANITIZE_ALL).map_err(|(_, e)| e)?;
            Ok(zip)
        })();
        res.map_err(|e| {
            let mut msg = format!(
                "Error:: in reaction {} :: building molecule from synthons :",
                self.id
            );
            for (i, &sn) in synth_nums_iter(synthon_nums) {
                let sep = if i > 0 { " and " } else { " " };
                let (id, idx) = &self.synthons[i][sn];
                let _ = write!(msg, "{sep}{id} ({})", pool[*idx].smiles());
            }
            msg.push_str(&format!("\n{e}\n"));
            RDKitError::Runtime(msg)
        })
    }

    fn build_sample_molecules(
        &self,
        pool: &[Synthon],
        synthon_mols: &[Vec<RWMol>],
        long_vec_num: usize,
    ) -> Result<Vec<ROMol>> {
        let mut sample_molecules = Vec::with_capacity(self.synthons[long_vec_num].len());
        // Find a small synthon in each of the non-longVecNum sets.
        let mut ex_synthons = vec![0usize; synthon_mols.len()];
        for i in 0..synthon_mols.len() {
            if i == long_vec_num {
                continue;
            }
            let mut min_atoms = synthon_mols[i][0].num_atoms();
            for j in 1..synthon_mols[i].len() {
                if synthon_mols[i][j].num_atoms() < min_atoms {
                    min_atoms = synthon_mols[i][j].num_atoms();
                    ex_synthons[i] = j;
                }
            }
        }
        for i in 0..self.synthons[long_vec_num].len() {
            ex_synthons[long_vec_num] = i;
            sample_molecules.push(self.build_molecule(pool, &ex_synthons)?);
        }
        Ok(sample_molecules)
    }

    /// Each synthon is built into a product and then split on the newly
    /// formed bonds; the fragment from the original synthon is its search mol.
    pub fn make_synthon_search_mols(&self, pool: &mut [Synthon]) -> Result<()> {
        let synthon_mol_copies = self.copy_synthons(pool)?;
        let dummy_labels: Vec<(u32, u32)> =
            (1..=MAX_CONNECTOR_NUM as u32).map(|i| (i, i)).collect();
        for synth_set_num in 0..self.synthons.len() {
            let sample_mols =
                self.build_sample_molecules(pool, &synthon_mol_copies, synth_set_num)?;
            for (j, sample) in sample_mols.iter().enumerate() {
                let split_bonds: Vec<u32> = sample
                    .bonds()
                    .iter()
                    .filter(|b| !b.has_prop("molNum"))
                    .map(|b| b.idx())
                    .collect();
                let frag_mol =
                    fragment_on_bonds(sample, &split_bonds, true, Some(&dummy_labels), None, None)?;
                let (mut mol_frags, _) = get_mol_frags_as_mols(&frag_mol, false, true)?;
                let frag_we_want = find_mol_num_frag(&mol_frags, synth_set_num)
                    .ok_or_else(|| RDKitError::Runtime("no fragment found for synthon".into()))?;
                let mut frag = mol_frags.swap_remove(frag_we_want);
                if !frag.ring_info().is_initialized() {
                    symmetrize_sssr(
                        &mut frag,
                        SymmetrizeSSSRAlgorithm::Default,
                        true,
                        false,
                        false,
                    )?;
                }
                let idx = self.synthons[synth_set_num][j].1;
                pool[idx].set_search_mol(frag)?;
            }
        }
        Ok(())
    }

    pub fn build_connector_regions(&mut self, pool: &[Synthon]) -> Result<()> {
        let mut smis = BTreeSet::new();
        for rset in &self.synthons {
            for (_, idx) in rset {
                for cr in pool[*idx].conn_regions() {
                    let smi = smiles_string(cr);
                    if smis.insert(smi.clone()) {
                        self.connector_regions.push(cr.clone());
                        self.conn_reg_smis.push(smi);
                    }
                }
            }
        }
        if !self.connector_regions.is_empty() {
            for cr in &self.connector_regions {
                self.conn_reg_fps
                    .push(pattern_fingerprint_mol(cr, PATT_FP_NUM_BITS, false)?);
            }
        }
        // all synthons in a set should have the same number of connections
        for set in &self.synthons {
            self.num_connectors
                .push(details::count_connections(pool[set[0].1].orig_mol()));
        }
        Ok(())
    }

    /// Finds "[1*]", "[1*:1]", "[2*]", "[2*:2]" etc. in the synthon SMILES.
    pub fn assign_connectors_used(&mut self, pool: &[Synthon]) {
        self.connectors = 0;
        self.synth_conn_patts.clear();
        for set in &self.synthons {
            let mut patt: ConnectorBits = 0;
            let smi = pool[set[0].1].smiles();
            for i in 1..=MAX_CONNECTOR_NUM {
                if smi.contains(&format!("[{i}*]")) || smi.contains(&format!("[{i}*:{i}]")) {
                    self.connectors |= 1 << i;
                    patt |= 1 << i;
                }
            }
            self.synth_conn_patts.push(patt);
        }
    }

    pub fn assess_ring_formers(&mut self) {
        self.num_ring_formers = 0;
        for i in 1..self.synthons.len() {
            for j in 0..i {
                if (self.synth_conn_patts[i] & self.synth_conn_patts[j]).count_ones() > 1 {
                    self.num_ring_formers += 1;
                }
            }
        }
    }

    pub fn build_product_name(&self, synth_nums: &[usize]) -> String {
        let names: Vec<&str> = synth_nums
            .iter()
            .enumerate()
            .map(|(i, &n)| self.synthons[i][n].0.as_str())
            .collect();
        details::build_product_name(&self.id, &names)
    }

    pub fn build_product(&self, pool: &[Synthon], synth_nums: &[usize]) -> Option<ROMol> {
        let synths: Vec<&ROMol> = synth_nums
            .iter()
            .enumerate()
            .map(|(i, &n)| pool[self.synthons[i][n].1].orig_mol())
            .collect();
        details::build_product(&synths)
    }

    pub fn initialize_search_orders(&mut self, pool: &[Synthon]) {
        if self.substructure_search_orders.is_empty() {
            self.substructure_search_orders = self.order_synthons_for_search(pool, |a, b| {
                a.search_mol().num_atoms() < b.search_mol().num_atoms()
            });
        }
        if self.rascal_search_orders.is_empty() {
            self.rascal_search_orders = self.order_synthons_for_search(pool, |a, b| {
                a.search_mol().num_atoms() + a.search_mol().num_bonds()
                    < b.search_mol().num_atoms() + b.search_mol().num_bonds()
            });
        }
        if self.has_fingerprints(pool) && self.fingerprint_search_orders.is_empty() {
            self.fingerprint_search_orders = self.order_synthons_for_search(pool, |a, b| {
                a.fp().map_or(0, |f| f.num_on_bits()) < b.fp().map_or(0, |f| f.num_on_bits())
            });
        }
    }

    fn order_synthons_for_search(
        &self,
        pool: &[Synthon],
        cmp: impl Fn(&Synthon, &Synthon) -> bool,
    ) -> Vec<Vec<usize>> {
        self.synthons
            .iter()
            .map(|set| {
                let mut v: Vec<(usize, usize)> = set
                    .iter()
                    .enumerate()
                    .map(|(j, (_, idx))| (*idx, j))
                    .collect();
                crate::rdgeneral::cxx_sort::libcxx_sort_by(&mut v, |a, b| {
                    cmp(&pool[a.0], &pool[b.0])
                });
                v.into_iter().map(|(_, j)| j).collect()
            })
            .collect()
    }

    /// `makeSampleMolecule`
    pub fn make_sample_molecule(
        &self,
        pool: &[Synthon],
        synthon_idx: usize,
    ) -> Option<SampleMolRec> {
        let mut found = None;
        'outer: for (i, set) in self.synthons.iter().enumerate() {
            for (j, (_, idx)) in set.iter().enumerate() {
                if *idx == synthon_idx {
                    found = Some((i, j));
                    break 'outer;
                }
            }
        }
        let (set_num, j) = found?;
        let mut ex = vec![0usize; self.synthons.len()];
        ex[set_num] = j;
        if !self.substructure_search_orders.is_empty() {
            for i in 0..self.substructure_search_orders.len() {
                if i != set_num {
                    ex[i] = self.substructure_search_orders[i][0];
                }
            }
        }
        let num_atoms = (0..self.synthons.len())
            .map(|i| pool[self.synthons[i][ex[i]].1].num_heavy_atoms())
            .sum();
        Some(SampleMolRec {
            synthon_idx,
            synthon_nums: ex,
            synthon_set_num: set_num,
            num_atoms,
        })
    }

    pub fn enumerate_to_stream(&self, pool: &[Synthon], out: &mut String) {
        let sizes: Vec<usize> = self.synthons.iter().map(|s| s.len()).collect();
        let mut stepper = Stepper::new(&sizes);
        while !stepper.at_end() {
            if let Some(prod) = self.build_product(pool, &stepper.curr_state) {
                let name = self.build_product_name(&stepper.curr_state);
                let _ = writeln!(out, "{} {}", smiles_string(&prod), name);
            }
            stepper.step();
        }
    }

    pub fn write_to_db_stream(&self, pool: &[Synthon], w: &mut Writer) -> Result<()> {
        w.string(self.id.as_bytes());
        w.u64(self.connector_regions.len() as u64);
        for cr in &self.connector_regions {
            w.mol(cr)?;
        }
        w.u64(self.conn_reg_fps.len() as u64);
        for fp in &self.conn_reg_fps {
            w.string(&fp.to_binary());
        }
        w.bits(self.connectors, MAX_CONNECTOR_NUM + 1);
        w.u64(self.synth_conn_patts.len() as u64);
        for &scp in &self.synth_conn_patts {
            w.bits(scp, MAX_CONNECTOR_NUM + 1);
        }
        w.u64(self.synthons.len() as u64);
        for rs in &self.synthons {
            w.u64(rs.len() as u64);
            for (id, idx) in rs {
                w.string(id.as_bytes());
                w.string(pool[*idx].smiles().as_bytes());
            }
        }
        if let (Some(a), Some(s)) = (&self.add_fp, &self.subtract_fp) {
            w.bool(true);
            w.string(&a.to_binary());
            w.string(&s.to_binary());
        } else {
            w.bool(false);
        }
        w.u64(self.num_connectors.len() as u64);
        for &c in &self.num_connectors {
            w.i32(c);
        }
        for orders in [
            &self.substructure_search_orders,
            &self.fingerprint_search_orders,
            &self.rascal_search_orders,
        ] {
            w.u64(orders.len() as u64);
            for o in orders {
                w.u64(o.len() as u64);
                for &v in o {
                    w.u64(v as u64);
                }
            }
        }
        Ok(())
    }

    /// Reads a SynthonSet. `lookup` finds a synthon in the pool by SMILES.
    pub fn read_from_db_stream(
        r: &mut Reader<'_>,
        lookup: &dyn Fn(&str) -> Option<usize>,
        file_name: &str,
        version: u32,
    ) -> Result<Self> {
        if version < 3000 {
            return Err(RDKitError::precondition(
                "Binary database version no longer supported.",
            ));
        }
        let mut s = SynthonSet {
            id: r.string()?,
            ..Default::default()
        };
        let n = r.u64()?;
        for _ in 0..n {
            let cr = r.mol()?;
            s.conn_reg_smis.push(smiles_string(&cr));
            s.connector_regions.push(cr);
        }
        let n = r.u64()?;
        for _ in 0..n {
            s.conn_reg_fps
                .push(bit_vect_from_binary(&r.bytes()?, false)?);
        }
        s.connectors = r.bits()?;
        let n = r.u64()?;
        for _ in 0..n {
            s.synth_conn_patts.push(r.bits()?);
        }
        let n = r.u64()?;
        for _ in 0..n {
            let num_r = r.u64()?;
            let mut set = Vec::new();
            for _ in 0..num_r {
                let name = r.string()?;
                let smiles = r.string()?;
                let Some(idx) = lookup(&smiles) else {
                    println!("smiles {smiles} not in pool");
                    return Err(RDKitError::Runtime(format!(
                        "Database file {file_name} appears corrupted."
                    )));
                };
                set.push((name, idx));
            }
            s.synthons.push(set);
        }
        if r.bool()? {
            s.add_fp = Some(bit_vect_from_binary(&r.bytes()?, false)?);
            s.subtract_fp = Some(bit_vect_from_binary(&r.bytes()?, false)?);
        }
        let n = r.u64()?;
        for _ in 0..n {
            s.num_connectors.push(r.i32()?);
        }
        if version >= 3020 {
            for which in 0..3 {
                let n = r.u64()? as usize;
                let mut orders = Vec::with_capacity(n);
                for _ in 0..n {
                    let t = r.u64()? as usize;
                    let mut o = Vec::with_capacity(t);
                    for _ in 0..t {
                        o.push(r.u64()? as usize);
                    }
                    orders.push(o);
                }
                match which {
                    0 => s.substructure_search_orders = orders,
                    1 => s.fingerprint_search_orders = orders,
                    _ => s.rascal_search_orders = orders,
                }
            }
        }
        Ok(s)
    }
}

fn synth_nums_iter(nums: &[usize]) -> impl Iterator<Item = (usize, &usize)> {
    nums.iter().enumerate()
}

fn add_join_info(mol: &mut ROMol, atom_idx: u32, isotope: u32) -> Result<()> {
    let atom = mol.atom_mut(atom_idx);
    let joins = match atom.get_prop("HoldsJoin") {
        Some(RDValue::String(s)) => format!("{s}_{isotope}"),
        _ => isotope.to_string(),
    };
    atom.set_prop("HoldsJoin", joins, false)
}

fn find_mol_num_frag(mol_frags: &[ROMol], mol_num: usize) -> Option<usize> {
    for (i, frag) in mol_frags.iter().enumerate() {
        for atom in frag.atoms() {
            if atom
                .get_prop("molNum")
                .and_then(|v| v.as_uint())
                .is_some_and(|v| v as usize == mol_num)
            {
                return Some(i);
            }
        }
    }
    None
}
