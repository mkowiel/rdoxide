//! Port of RDKit `Code/GraphMol/FragCatalog`: catalogs of the fragments
//! (connected subgraphs, with functional groups marked) found in a set of
//! molecules, and fingerprints made from them.

use std::collections::{BTreeMap, HashMap};
use std::sync::{Arc, OnceLock};

use crate::catalogs::{CatalogEntry, CatalogParams, HierarchCatalog, StreamReader, StreamWriter};
use crate::datastructs::{BitVect, ExplicitBitVect};
use crate::graphmol::atom::set_supplemental_smiles_label;
use crate::graphmol::fingerprints::hash_vect;
use crate::graphmol::pickle::{get_default_pickle_properties, mol_from_pickle_stream, pickle_mol};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::smiles_parse::write::mol_to_smiles;
use crate::graphmol::subgraphs::{
    PathMap, bond_list_from_atom_list, calc_path_discriminators,
    find_all_subgraphs_of_lengths_m_to_n, path_to_submol,
};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::utils::{
    compute_int_vect_primes_product, feq, first_thousand_primes, intersect,
};
use crate::{RDKitError, Result};

/// `getSmarts`: the functional group on a line of a group file (`None` for
/// blank and comment lines).
fn func_group_from_line(line: &str) -> Result<Option<ROMol>> {
    let line = line.trim();
    if line.is_empty() || line.starts_with("//") {
        return Ok(None);
    }
    // boost's char_separator drops empty tokens
    let mut tokens = line.split('\t').filter(|t| !t.is_empty());
    let name = tokens.next().unwrap_or_default().replace(' ', "");
    let smarts = tokens
        .next()
        .ok_or_else(|| RDKitError::invariant(line))?
        .replace(' ', "");
    let mut mol = mol_from_smarts(&smarts, &SmartsParserParams::default())?
        .ok_or_else(|| RDKitError::invariant(smarts.clone()))?;
    mol.set_prop(cp::_Name, name, false)?;
    mol.set_prop(cp::_fragSMARTS, smarts, false)?;
    Ok(Some(mol))
}

/// `readFuncGroups(file)`: the functional groups of a tab-separated file.
pub fn read_func_groups_file(file_name: &str) -> Result<Vec<Arc<ROMol>>> {
    let text = std::fs::read_to_string(file_name)
        .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))?;
    read_func_groups(&text, None)
}

/// `readFuncGroups(stream, nToRead)`
pub fn read_func_groups(text: &str, n_to_read: Option<usize>) -> Result<Vec<Arc<ROMol>>> {
    let mut groups = Vec::new();
    for line in text.lines() {
        if n_to_read.is_some_and(|n| groups.len() >= n) {
            break;
        }
        if let Some(mol) = func_group_from_line(line)? {
            groups.push(Arc::new(mol));
        }
    }
    Ok(groups)
}

/// C++ stream formatting of a double (`%g` with six significant digits).
fn format_double(v: f64) -> String {
    if v == 0.0 {
        return "0".into();
    }
    let sci = format!("{v:.5e}");
    let (mantissa, exp) = sci.split_once('e').expect("exponent");
    let exp: i32 = exp.parse().expect("exponent");
    let trim = |s: &str| {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s.to_string()
        }
    };
    if !(-5..6).contains(&exp) {
        format!(
            "{}e{}{:02}",
            trim(mantissa),
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        )
    } else {
        trim(&format!("{:.*}", (5 - exp).max(0) as usize, v))
    }
}

/// `FragCatParams`
#[derive(Debug, Clone)]
pub struct FragCatParams {
    lower_frag_len: u32,
    upper_frag_len: u32,
    tolerance: f64,
    func_groups: Vec<Arc<ROMol>>,
}

impl Default for FragCatParams {
    fn default() -> Self {
        FragCatParams {
            lower_frag_len: 0,
            upper_frag_len: 0,
            tolerance: 1e-8,
            func_groups: Vec::new(),
        }
    }
}

impl FragCatParams {
    /// `FragCatParams(lLen, uLen, fgroupFile, tol)`
    pub fn new(lower: u32, upper: u32, fgroup_file: &str, tol: f64) -> Result<Self> {
        if lower > upper {
            return Err(RDKitError::invariant(
                "The upper length for fragments must be >= lower length",
            ));
        }
        Ok(FragCatParams {
            lower_frag_len: lower,
            upper_frag_len: upper,
            tolerance: tol,
            func_groups: read_func_groups_file(fgroup_file)?,
        })
    }

    pub fn lower_frag_length(&self) -> u32 {
        self.lower_frag_len
    }
    pub fn upper_frag_length(&self) -> u32 {
        self.upper_frag_len
    }
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }
    pub fn num_func_groups(&self) -> u32 {
        self.func_groups.len() as u32
    }
    pub fn func_groups(&self) -> &[Arc<ROMol>] {
        &self.func_groups
    }
    pub fn func_group(&self, fid: u32) -> Result<&ROMol> {
        self.func_groups
            .get(fid as usize)
            .map(|m| &**m)
            .ok_or_else(|| RDKitError::Index(fid.to_string()))
    }

    /// `Serialize`: the text form
    pub fn serialize(&self) -> Result<String> {
        let mut out = format!(
            "{} {} {}\n{}\n",
            self.lower_frag_len,
            self.upper_frag_len,
            format_double(self.tolerance),
            self.func_groups.len()
        );
        for g in &self.func_groups {
            let prop = |key: &str| {
                g.get_prop(key)
                    .map(|v| v.to_rd_string())
                    .ok_or_else(|| RDKitError::Key(key.to_string()))
            };
            out.push_str(&format!(
                "{}\t{}\n",
                prop(cp::_Name)?,
                prop(cp::_fragSMARTS)?
            ));
        }
        Ok(out)
    }

    /// `initFromString`
    pub fn from_text(text: &str) -> Result<Self> {
        let mut params = FragCatParams::default();
        params.parse(text)?;
        Ok(params)
    }

    /// Reads the text form from the start of `text`; returns the bytes used.
    fn parse(&mut self, text: &str) -> Result<usize> {
        let bad = || RDKitError::Value("Bad pickle format: bad fragment catalog parameters".into());
        let mut pos = 0;
        let mut token = || -> Result<&str> {
            let rest = &text[pos..];
            let start = rest.find(|c: char| !c.is_whitespace()).ok_or_else(bad)?;
            let len = rest[start..]
                .find(char::is_whitespace)
                .unwrap_or(rest.len() - start);
            pos += start + len;
            Ok(&rest[start..start + len])
        };
        self.lower_frag_len = token()?.parse().map_err(|_| bad())?;
        self.upper_frag_len = token()?.parse().map_err(|_| bad())?;
        self.tolerance = token()?.parse().map_err(|_| bad())?;
        let n: usize = token()?.parse().map_err(|_| bad())?;
        // the rest of the line holding the count, then one line per group
        let mut groups = Vec::new();
        let mut rest = &text[pos..];
        while groups.len() < n {
            let Some(nl) = rest.find('\n') else {
                if let Some(g) = func_group_from_line(rest)? {
                    groups.push(Arc::new(g));
                }
                pos = text.len();
                break;
            };
            if let Some(g) = func_group_from_line(&rest[..nl])? {
                groups.push(Arc::new(g));
            }
            pos += nl + 1;
            rest = &rest[nl + 1..];
        }
        self.func_groups = groups;
        Ok(pos)
    }
}

impl CatalogParams for FragCatParams {
    fn type_str(&self) -> &str {
        "Fragment Catalog Parameters"
    }
    fn to_stream(&self, w: &mut StreamWriter) {
        w.bytes(self.serialize().unwrap_or_default().as_bytes());
    }
    fn init_from_stream(&mut self, r: &mut StreamReader) -> Result<()> {
        let text = String::from_utf8_lossy(&r.data[r.pos..]).into_owned();
        let used = self.parse(&text)?;
        // the text was decoded lossily; count the bytes of the consumed part
        r.pos += text[..used].len();
        Ok(())
    }
}

/// `MatchVectType` of (atom, functional group id) pairs
pub type AtomFuncGroups = Vec<(i32, i32)>;

/// `INT_INT_VECT_MAP`: atom index -> functional group ids
pub type AtomToFuncGroups = BTreeMap<u32, Vec<i32>>;

/// `findFuncGroupsOnMol`: (atom, functional group) pairs, and the bonds the
/// functional groups cover.
pub fn find_func_groups_on_mol(
    mol: &ROMol,
    params: &FragCatParams,
) -> Result<(AtomFuncGroups, Vec<i32>)> {
    let mut fg_bonds: Vec<i32> = Vec::new();
    let mut aid_fgrps = Vec::new();
    let ps = SubstructMatchParameters::default();
    // the functional groups are assumed to come in decreasing order of size
    for (fid, fgrp) in params.func_groups().iter().enumerate() {
        for m in substruct_match(mol, fgrp, &ps)? {
            // the first atom of the group is the connection atom
            let aid = m.first().map_or(0, |p| p.1);
            let maids: Vec<i32> = m.iter().map(|p| p.1).collect();
            let mut all_done = true;
            for bi in bond_list_from_atom_list(mol, &maids) {
                if !fg_bonds.contains(&bi) {
                    all_done = false;
                    fg_bonds.push(bi);
                }
            }
            if !all_done {
                aid_fgrps.push((aid, fid as i32));
            }
        }
    }
    Ok((aid_fgrps, fg_bonds))
}

/// `prepareMol`: the part of the molecule without the functional groups, and
/// the functional groups on it.
pub fn prepare_mol(mol: &ROMol, params: &FragCatParams) -> Result<(ROMol, AtomFuncGroups)> {
    let (aid_to_fid, fg_bonds) = find_func_groups_on_mol(mol, params)?;
    let core_bonds: Vec<i32> = (0..mol.num_bonds() as i32)
        .filter(|b| !fg_bonds.contains(b))
        .collect();
    let (core, idx_map) = path_to_submol(mol, &core_bonds, false)?;
    let mut a_to_f = Vec::new();
    for (aid, fid) in aid_to_fid {
        a_to_f.push((idx_map.get(&(aid as u32)).copied().unwrap_or(0) as i32, fid));
    }
    Ok((core, a_to_f))
}

/// `FragCatalogEntry`
#[derive(Debug, Clone)]
pub struct FragCatalogEntry {
    mol: ROMol,
    bit_id: i32,
    description: String,
    order: u32,
    a_to_f: AtomToFuncGroups,
    discrims: OnceLock<(u32, usize, usize)>,
}

impl FragCatalogEntry {
    /// The entry for the fragment `path` of `omol`, with the functional
    /// groups `aid_to_fid` found on `omol`.
    pub fn new(omol: &ROMol, path: &[i32], aid_to_fid: &[(i32, i32)]) -> Result<Self> {
        let (mol, idx_map) = path_to_submol(omol, path, false)?;
        let mut a_to_f = AtomToFuncGroups::new();
        for &(old, fid) in aid_to_fid {
            if let Some(&new) = idx_map.get(&(old as u32)) {
                a_to_f.entry(new).or_default().push(fid);
            }
        }
        Ok(FragCatalogEntry {
            mol,
            bit_id: -1,
            description: String::new(),
            order: path.len() as u32,
            a_to_f,
            discrims: OnceLock::new(),
        })
    }

    pub fn mol(&self) -> &ROMol {
        &self.mol
    }

    pub fn func_group_map(&self) -> &AtomToFuncGroups {
        &self.a_to_f
    }

    /// `setDescription(params)`: labels the functional groups on the
    /// fragment and describes it with its SMILES.
    pub fn set_description(&mut self, params: &FragCatParams) -> Result<()> {
        for (&at_idx, groups) in &self.a_to_f {
            let name = |fid: i32| -> Result<String> {
                let g = params.func_group(fid as u32)?;
                Ok(g.get_prop(cp::_Name)
                    .map(|v| v.to_rd_string())
                    .unwrap_or_default())
            };
            let mut label = String::new();
            let (last, init) = groups.split_last().expect("at least one group");
            for &fid in init {
                label.push_str(&format!("(<{}>)", name(fid)?));
            }
            label.push_str(&format!("<{}>", name(*last)?));
            set_supplemental_smiles_label(self.mol.atom_mut(at_idx), &label);
        }
        self.description = mol_to_smiles(&mut self.mol, &Default::default())?;
        Ok(())
    }

    /// `getDiscrims`
    pub fn discrims(&self) -> Result<(u32, usize, usize)> {
        if let Some(d) = self.discrims.get() {
            return Ok(*d);
        }
        let path: Vec<i32> = (0..self.mol.num_bonds() as i32).collect();
        let invars: Vec<u32> = (0..self.mol.num_atoms())
            .map(|aid| {
                self.a_to_f.get(&aid).map_or(0, |groups| {
                    let mut g: Vec<u32> = groups.iter().map(|&x| x as u32).collect();
                    g.sort_unstable();
                    hash_vect(&g)
                })
            })
            .collect();
        let d = calc_path_discriminators(&self.mol, &path, true, Some(&invars))?;
        Ok(*self.discrims.get_or_init(|| d))
    }

    /// `match`
    pub fn matches(&self, other: &FragCatalogEntry, tol: f64) -> Result<bool> {
        if self.order != other.order {
            return Ok(false);
        }
        if other.a_to_f.len() != self.a_to_f.len() {
            return Ok(false);
        }
        let found_all = self
            .a_to_f
            .values()
            .all(|t| other.a_to_f.values().any(|o| o == t));
        if !found_all {
            return Ok(false);
        }
        let (t, o) = (self.discrims()?, other.discrims()?);
        Ok(feq(f64::from(t.0), f64::from(o.0), tol)
            && feq(t.1 as f64, o.1 as f64, tol)
            && feq(t.2 as f64, o.2 as f64, tol))
    }

    /// `Serialize`
    pub fn serialize(&self) -> Result<Vec<u8>> {
        let mut w = StreamWriter::default();
        self.to_stream(&mut w)?;
        Ok(w.0)
    }
}

impl CatalogEntry for FragCatalogEntry {
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
        self.mol.num_bonds() as i32
    }

    fn to_stream(&self, w: &mut StreamWriter) -> Result<()> {
        w.bytes(&pickle_mol(&self.mol, get_default_pickle_properties())?);
        w.i32(self.bit_id);
        w.string(&self.description);
        w.i32(self.order as i32);
        w.i32(self.a_to_f.len() as i32);
        for (&k, v) in &self.a_to_f {
            w.i32(k as i32);
            w.i32(v.len() as i32);
            v.iter().for_each(|&x| w.i32(x));
        }
        Ok(())
    }

    fn from_stream(r: &mut StreamReader) -> Result<Self> {
        let mol = mol_from_pickle_stream(r.data, &mut r.pos, get_default_pickle_properties())?;
        let bit_id = r.i32()?;
        let description = r.string()?;
        let order = r.i32()? as u32;
        let n = r.i32()?;
        let mut a_to_f = AtomToFuncGroups::new();
        for _ in 0..n {
            let key = r.i32()? as u32;
            let size = r.i32()?;
            let vals = (0..size).map(|_| r.i32()).collect::<Result<Vec<_>>>()?;
            a_to_f.insert(key, vals);
        }
        Ok(FragCatalogEntry {
            mol,
            bit_id,
            description,
            order,
            a_to_f,
            discrims: OnceLock::new(),
        })
    }
}

/// `FragCatalog`
pub type FragCatalog = HierarchCatalog<FragCatalogEntry, FragCatParams>;

/// The key a path's invariant is looked up by (`DOUBLE_INT_MAP`).
type InvarMap = HashMap<u64, i32>;

fn params_of(fcat: &FragCatalog) -> Result<FragCatParams> {
    fcat.catalog_params()
        .cloned()
        .ok_or_else(|| RDKitError::precondition("NULL parameter object"))
}

fn down_list(fcat: &FragCatalog, id: i32) -> Result<Vec<i32>> {
    Ok(fcat
        .down_entry_list(id as u32)?
        .iter()
        .map(|&x| x as i32)
        .collect())
}

/// The order k-1 catalog entries below `path`'s order k-1 subpaths and the
/// intersection of their down entries; subpaths never seen are skipped.
/// `stop_on_unmatched` (the fingerprint code): a subpath without an entry
/// (-1) empties the list.
fn subpath_candidates(
    fcat: &FragCatalog,
    path: &[i32],
    invar: f64,
    mapkm1: &InvarMap,
    stop_on_unmatched: bool,
    sub_invars: &mut Vec<f64>,
) -> Result<Vec<i32>> {
    let primes = first_thousand_primes();
    let mut candidates: Vec<i32> = Vec::new();
    let mut scnt = 0;
    for &b in path {
        let sinvar = invar / f64::from(primes[b as usize % primes.len()]);
        let Some(&ent_id) = mapkm1.get(&sinvar.to_bits()) else {
            continue;
        };
        if stop_on_unmatched && ent_id == -1 {
            return Ok(Vec::new());
        }
        sub_invars.push(sinvar);
        let down = down_list(fcat, ent_id)?;
        candidates = if scnt == 0 {
            down
        } else {
            intersect(&down, &candidates)
        };
        scnt += 1;
    }
    Ok(candidates)
}

/// `FragCatGenerator`
#[derive(Debug, Clone, Copy, Default)]
pub struct FragCatGenerator;

impl FragCatGenerator {
    /// `addFragsFromMol`: adds the fragments of `mol` to `fcat`; returns the
    /// number of new entries.
    pub fn add_frags_from_mol(&self, mol: &ROMol, fcat: &mut FragCatalog) -> Result<u32> {
        let fparams = params_of(fcat)?;
        let (l_len, u_len) = (fparams.lower_frag_length(), fparams.upper_frag_length());
        if l_len > u_len {
            return Err(RDKitError::invariant(""));
        }
        let (core, aid_to_fid) = prepare_mol(mol, &fparams)?;
        let all_paths = find_all_subgraphs_of_lengths_m_to_n(&core, 1, u_len, false, -1, None)?;
        let mut mapkm1 = InvarMap::new();
        let n1 = add_order1_paths(&all_paths, &core, fcat, &fparams, &mut mapkm1, &aid_to_fid)?;
        let nrem =
            add_higher_order_paths(&all_paths, &core, fcat, &fparams, &mut mapkm1, &aid_to_fid)?;
        Ok(n1 + nrem)
    }
}

fn add_order1_paths(
    all_paths: &PathMap,
    mol: &ROMol,
    fcat: &mut FragCatalog,
    fparams: &FragCatParams,
    mapkm1: &mut InvarMap,
    aid_to_fid: &[(i32, i32)],
) -> Result<u32> {
    let (l_len, u_len) = (fparams.lower_frag_length(), fparams.upper_frag_length());
    let tol = fparams.tolerance();
    let mut n01 = 0;
    for path in all_paths.get(&1).map_or(&[][..], Vec::as_slice) {
        let mut nent = FragCatalogEntry::new(mol, path, aid_to_fid)?;
        let invar = compute_int_vect_primes_product(path);
        let o1_entries: Vec<u32> = fcat
            .entries_of_order(1)
            .map(<[u32]>::to_vec)
            .unwrap_or_default();
        let mut found = false;
        for eti in o1_entries {
            if nent.matches(fcat.entry_with_idx(eti)?, tol)? {
                found = true;
                mapkm1.insert(invar.to_bits(), eti as i32);
                break;
            }
        }
        if !found {
            // order 1 fragments get a bit when they are of interest
            let update_sig_l = l_len <= 1 && u_len >= 1;
            nent.set_description(fparams)?;
            let vid = fcat.add_entry(nent, update_sig_l);
            mapkm1.insert(invar.to_bits(), vid as i32);
            n01 += 1;
        }
    }
    Ok(n01)
}

fn add_higher_order_paths(
    all_paths: &PathMap,
    mol: &ROMol,
    fcat: &mut FragCatalog,
    fparams: &FragCatParams,
    mapkm1: &mut InvarMap,
    aid_to_fid: &[(i32, i32)],
) -> Result<u32> {
    let (l_len, u_len) = (fparams.lower_frag_length(), fparams.upper_frag_length());
    let tol = fparams.tolerance();
    let mut nrem = 0;
    for (&order, paths) in all_paths {
        if order < 2 {
            continue;
        }
        let mut mapk = InvarMap::new();
        for path in paths {
            let mut nent = FragCatalogEntry::new(mol, path, aid_to_fid)?;
            nent.set_description(fparams)?;
            let invar = compute_int_vect_primes_product(path);
            let mut sub_invars = Vec::new();
            let candidates = subpath_candidates(fcat, path, invar, mapkm1, false, &mut sub_invars)?;
            let mut found = None;
            for iti in candidates {
                if nent.matches(fcat.entry_with_idx(iti as u32)?, tol)? {
                    found = Some(iti);
                    break;
                }
            }
            if let Some(m) = found {
                mapk.insert(invar.to_bits(), m);
            } else {
                let ordr = nent.order() as u32;
                let update_sig_l = ordr >= l_len && ordr <= u_len;
                let vid = fcat.add_entry(nent, update_sig_l);
                mapk.insert(invar.to_bits(), vid as i32);
                nrem += 1;
                for s in sub_invars {
                    let ent_id = mapkm1[&s.to_bits()];
                    fcat.add_edge(ent_id as u32, vid)?;
                }
            }
        }
        *mapkm1 = mapk;
    }
    Ok(nrem)
}

/// `FragFPGenerator`
#[derive(Debug, Clone, Copy, Default)]
pub struct FragFPGenerator;

impl FragFPGenerator {
    /// `getFPForMol`
    pub fn fp_for_mol(&self, mol: &ROMol, fcat: &FragCatalog) -> Result<ExplicitBitVect> {
        let mut fp = ExplicitBitVect::new(fcat.fp_length(), false);
        let fparams = params_of(fcat)?;
        let (core, aid_to_fid) = prepare_mol(mol, &fparams)?;
        self.compute_fp(&core, fcat, &fparams, &aid_to_fid, &mut fp)?;
        Ok(fp)
    }

    fn compute_fp(
        &self,
        mol: &ROMol,
        fcat: &FragCatalog,
        fparams: &FragCatParams,
        aid_to_fid: &[(i32, i32)],
        fp: &mut ExplicitBitVect,
    ) -> Result<()> {
        let tol = fparams.tolerance();
        let all_paths = find_all_subgraphs_of_lengths_m_to_n(
            mol,
            1,
            fparams.upper_frag_length(),
            false,
            -1,
            None,
        )?;
        let mut mapkm1 = InvarMap::new();
        let o1_entries = fcat.entries_of_order(1)?;
        for path in all_paths.get(&1).map_or(&[][..], Vec::as_slice) {
            let mut nent = FragCatalogEntry::new(mol, path, aid_to_fid)?;
            nent.set_description(fparams)?;
            let invar = compute_int_vect_primes_product(path);
            // -1 marks "no match in the catalog"
            mapkm1.insert(invar.to_bits(), -1);
            for &eti in o1_entries {
                let entry = fcat.entry_with_idx(eti)?;
                if nent.matches(entry, tol)? {
                    if entry.bit_id() >= 0 {
                        fp.set_bit(entry.bit_id() as u32)?;
                    }
                    mapkm1.insert(invar.to_bits(), eti as i32);
                    break;
                }
            }
        }
        for (&order, paths) in &all_paths {
            if order < 2 {
                continue;
            }
            let mut mapk = InvarMap::new();
            for path in paths {
                let invar = compute_int_vect_primes_product(path);
                mapk.insert(invar.to_bits(), -1);
                let mut nent = FragCatalogEntry::new(mol, path, aid_to_fid)?;
                nent.set_description(fparams)?;
                let candidates =
                    subpath_candidates(fcat, path, invar, &mapkm1, true, &mut Vec::new())?;
                for iti in candidates {
                    let entry = fcat.entry_with_idx(iti as u32)?;
                    if nent.matches(entry, tol)? {
                        mapk.insert(invar.to_bits(), iti);
                        if entry.bit_id() >= 0 {
                            fp.set_bit(entry.bit_id() as u32)?;
                        }
                        break;
                    }
                }
            }
            mapkm1 = mapk;
        }
        Ok(())
    }
}
