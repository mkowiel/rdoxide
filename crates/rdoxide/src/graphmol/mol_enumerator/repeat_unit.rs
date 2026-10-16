//! Repeat unit (SRU) enumeration (upstream `MolEnumerator/RepeatUnit.cpp`):
//! the atoms of a SRU substance group are repeated, `HT` (head to tail) or
//! `HH` (head to head) connected. The SRU labels give the numbers of repeats.
//!
//! A work in progress upstream too; SRUs may not overlap.

use std::collections::BTreeMap;

use super::MolEnumeratorOp;
use crate::graphmol::bond::BondType;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::{RDProps, RDValue, common_properties as cp};
use crate::{RDKitError, Result, rd_log_warning};

const POLYMARKER: &str = "_polymeratom";
const HEADMARKER: &str = "_headatom";
const HEADMARKER_FRAME: &str = "_headatom_frame";
const TAILMARKER: &str = "_tailatom";
const TAILMARKER_FRAME: &str = "_tailatom_frame";
const HEADHEADMARKER: &str = "_headhead";
const LADDER_OFFSET: u32 = 100_000;

/// The repeat count when the SRU has no usable label
pub const DEFAULT_REPEAT_COUNT: usize = 4;

/// `RepeatUnitOp`
#[derive(Debug, Clone)]
pub struct RepeatUnitOp {
    /// the largest number of repeat counts to enumerate
    pub max_num_rounds: usize,
    mol: Option<ROMol>,
    frame: Option<ROMol>,
    repeats: Vec<ROMol>,
    count_at_each_point: Vec<usize>,
    min_repeat_counts: Vec<usize>,
}

impl Default for RepeatUnitOp {
    fn default() -> Self {
        RepeatUnitOp {
            max_num_rounds: usize::MAX,
            mol: None,
            frame: None,
            repeats: Vec::new(),
            count_at_each_point: Vec::new(),
            min_repeat_counts: Vec::new(),
        }
    }
}

fn uint_prop(mol: &ROMol, atom: u32, name: &str) -> Option<u32> {
    mol.atom(atom).get_prop(name).and_then(RDValue::as_uint)
}

fn vec_prop(mol: &ROMol, atom: u32, name: &str) -> Option<Vec<u32>> {
    match mol.atom(atom).get_prop(name) {
        Some(RDValue::VecUnsignedInt(v)) => Some(v.clone()),
        _ => None,
    }
}

/// `tagAtoms`
fn tag_atoms(
    mol: &mut ROMol,
    bond: u32,
    sgatoms: &[bool],
    index: u32,
    marker: &str,
    frame_marker: &str,
    connect: &str,
) -> Result<()> {
    let (begin, end) = (
        mol.bond(bond).begin_atom_idx(),
        mol.bond(bond).end_atom_idx(),
    );
    let (inside, frame_atom) = if sgatoms[begin as usize] {
        (begin, end)
    } else if sgatoms[end as usize] {
        (end, begin)
    } else {
        return Err(RDKitError::Value(
            "neither atom in an SRU bond is in the polymer".into(),
        ));
    };
    mol.atom_mut(inside).set_prop(marker, index, false)?;
    if connect == "HH" {
        mol.atom_mut(inside).set_prop(HEADHEADMARKER, 1u32, false)?;
    }
    let mut vs = vec_prop(mol, frame_atom, frame_marker).unwrap_or_default();
    vs.push(index);
    mol.atom_mut(frame_atom).set_prop(frame_marker, vs, false)
}

/// Digits at the start of the text: the number and the rest
/// (`std::from_chars`).
fn leading_number(text: &str) -> Option<(usize, &str)> {
    let digits = text.bytes().take_while(u8::is_ascii_digit).count();
    if digits == 0 {
        return None;
    }
    let n = text[..digits].parse::<usize>().ok()?;
    Some((n, &text[digits..]))
}

/// `parse_repeat_counts`: (minimum, maximum + 1) repeats from an SRU label
fn parse_repeat_counts(label: &str) -> (usize, usize) {
    const DEFAULT: (usize, usize) = (0, DEFAULT_REPEAT_COUNT);
    if label.is_empty() {
        return DEFAULT;
    }
    let Some((min_repeats, rest)) = leading_number(label) else {
        return DEFAULT;
    };
    if rest.is_empty() {
        return (0, min_repeats + 1);
    }
    let Some(rest) = rest.strip_prefix('-') else {
        return DEFAULT;
    };
    match leading_number(rest) {
        Some((max_repeats, "")) if max_repeats >= min_repeats => (min_repeats, max_repeats + 1),
        _ => DEFAULT,
    }
}

/// What the enumeration needs of a SRU substance group
struct Sru {
    atoms: Vec<u32>,
    bonds: Vec<u32>,
    connect: String,
    index: u32,
    xbcorr: Vec<u32>,
    min_repeats: usize,
    max_repeats: usize,
}

impl RepeatUnitOp {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_mol(mol: &ROMol) -> Result<Self> {
        let mut op = Self::new();
        op.init_from_mol(mol)?;
        Ok(op)
    }
}

/// `flipHeadHeadGroups`
fn flip_head_head_groups(orig_atom_count: u32, mol: &mut ROMol) -> Result<()> {
    for aidx in orig_atom_count..mol.num_atoms() {
        if !mol.atom(aidx).has_prop(HEADHEADMARKER) {
            continue;
        }
        if let Some(v) = uint_prop(mol, aidx, HEADMARKER) {
            let atom = mol.atom_mut(aidx);
            atom.set_prop(TAILMARKER, v, false)?;
            atom.clear_prop(HEADMARKER);
        } else if let Some(v) = uint_prop(mol, aidx, TAILMARKER) {
            let atom = mol.atom_mut(aidx);
            atom.set_prop(HEADMARKER, v, false)?;
            atom.clear_prop(TAILMARKER);
        }
        if let Some(v) = vec_prop(mol, aidx, HEADMARKER_FRAME) {
            let atom = mol.atom_mut(aidx);
            atom.set_prop(TAILMARKER_FRAME, v, false)?;
            atom.clear_prop(HEADMARKER_FRAME);
        } else if let Some(v) = vec_prop(mol, aidx, TAILMARKER_FRAME) {
            let atom = mol.atom_mut(aidx);
            atom.set_prop(HEADMARKER_FRAME, v, false)?;
            atom.clear_prop(TAILMARKER_FRAME);
        }
    }
    Ok(())
}

/// Removes the neighbors of an atom that are not part of the polymer.
fn remove_non_polymer_neighbors(mol: &mut ROMol, atom: u32) -> Result<()> {
    let nbrs: Vec<u32> = mol.atom_neighbors(atom).collect();
    for nbr in nbrs {
        if !mol.atom(nbr).has_prop(POLYMARKER) {
            mol.remove_atom(nbr)?;
        }
    }
    Ok(())
}

/// `connectRepeatAtomsAndRemoveExtras`
fn connect_repeat_atoms_and_remove_extras(orig_atom_count: u32, mol: &mut ROMol) -> Result<()> {
    mol.begin_batch_edit()?;
    for aidx1 in 0..orig_atom_count {
        let Some(tail_idx) = uint_prop(mol, aidx1, TAILMARKER) else {
            continue;
        };
        let mut connected = false;
        for aidx2 in orig_atom_count..mol.num_atoms() {
            if uint_prop(mol, aidx2, HEADMARKER) == Some(tail_idx) {
                connected = true;
                mol.atom_mut(aidx2).clear_prop(HEADMARKER);
                remove_non_polymer_neighbors(mol, aidx2)?;
                mol.add_bond(aidx1, aidx2, BondType::Single)?;
                break;
            }
        }
        if !connected {
            return Err(RDKitError::Value("no head found for tail".into()));
        }
        mol.atom_mut(aidx1).clear_prop(TAILMARKER);
        remove_non_polymer_neighbors(mol, aidx1)?;
    }
    mol.commit_batch_edit()
}

/// `connectRepeatToFrame`
fn connect_repeat_to_frame(
    n_orig_atoms: u32,
    mol: &mut ROMol,
    head_map: &mut BTreeMap<u32, u32>,
    tail_map: &mut BTreeMap<u32, u32>,
) -> Result<()> {
    for aidx in n_orig_atoms..mol.num_atoms() {
        if let Some(val) = uint_prop(mol, aidx, HEADMARKER)
            && let Some(&tail) = tail_map.get(&val)
        {
            mol.add_bond(aidx, tail, BondType::Single)?;
            mol.atom_mut(aidx).clear_prop(HEADMARKER);
            tail_map.remove(&val);
        }
        if let Some(val) = uint_prop(mol, aidx, TAILMARKER)
            && let Some(&head) = head_map.get(&val)
        {
            mol.add_bond(aidx, head, BondType::Single)?;
            mol.atom_mut(aidx).clear_prop(TAILMARKER);
            head_map.remove(&val);
        }
        if let Some(vals) = vec_prop(mol, aidx, HEADMARKER_FRAME) {
            for val in vals {
                head_map.insert(val, aidx);
            }
            mol.atom_mut(aidx).clear_prop(HEADMARKER_FRAME);
        }
        if let Some(vals) = vec_prop(mol, aidx, TAILMARKER_FRAME) {
            for val in vals {
                tail_map.insert(val, aidx);
            }
            mol.atom_mut(aidx).clear_prop(TAILMARKER_FRAME);
        }
    }
    Ok(())
}

/// `constructHeadAndTailMaps`
fn construct_head_and_tail_maps(
    mol: &ROMol,
    head_map: &mut BTreeMap<u32, u32>,
    tail_map: &mut BTreeMap<u32, u32>,
) -> Result<()> {
    for aidx in 0..mol.num_atoms() {
        if let Some(vals) = vec_prop(mol, aidx, HEADMARKER_FRAME) {
            for val in vals {
                if head_map.insert(val, aidx).is_some() {
                    return Err(RDKitError::Value(
                        "SRU group present with multiple head atoms with the same index".into(),
                    ));
                }
            }
        }
        if let Some(vals) = vec_prop(mol, aidx, TAILMARKER_FRAME) {
            for val in vals {
                if tail_map.insert(val, aidx).is_some() {
                    return Err(RDKitError::Value(
                        "SRU group present with multiple tail atoms with the same index".into(),
                    ));
                }
            }
        }
    }
    Ok(())
}

impl MolEnumeratorOp for RepeatUnitOp {
    fn init_from_mol(&mut self, mol: &ROMol) -> Result<()> {
        let mut mol = mol.clone();
        self.repeats.clear();
        self.count_at_each_point.clear();
        let mut atoms_in_srus = vec![false; mol.num_atoms() as usize];
        let mut srus: Vec<Sru> = Vec::new();
        for sg in mol.substance_groups() {
            let is_sru = sg.get_prop("TYPE").map(RDValue::to_rd_string).as_deref() == Some("SRU");
            if !is_sru {
                continue;
            }
            let connect = sg
                .get_prop("CONNECT")
                .map(RDValue::to_rd_string)
                .unwrap_or_default();
            if !connect.is_empty() && connect != "HT" && connect != "HH" {
                rd_log_warning!("can only enumerate SRUs with CONNECT=HT or CONNECT=HH\n");
                continue;
            }
            let label = sg
                .get_prop("LABEL")
                .map(RDValue::to_rd_string)
                .unwrap_or_default();
            let (min_repeats, max_repeats) = parse_repeat_counts(&label);
            for &aidx in sg.atoms() {
                if atoms_in_srus[aidx as usize] {
                    return Err(RDKitError::Value(
                        "cannot enumerate overlapping SRU groups".into(),
                    ));
                }
                atoms_in_srus[aidx as usize] = true;
            }
            srus.push(Sru {
                atoms: sg.atoms().to_vec(),
                bonds: sg.bonds().to_vec(),
                connect,
                index: sg
                    .get_prop("index")
                    .and_then(RDValue::as_uint)
                    .ok_or_else(|| RDKitError::Key("index".into()))?,
                xbcorr: match sg.get_prop("XBCORR") {
                    Some(RDValue::VecUnsignedInt(v)) => v.clone(),
                    _ => Vec::new(),
                },
                min_repeats,
                max_repeats,
            });
        }
        // tagging comes after the scan of all groups, which does not depend on it
        for sru in &srus {
            self.min_repeat_counts.push(sru.min_repeats);
            self.count_at_each_point
                .push((sru.max_repeats - sru.min_repeats).min(self.max_num_rounds));
            let mut sgatoms = vec![false; mol.num_atoms() as usize];
            for &aidx in &sru.atoms {
                sgatoms[aidx as usize] = true;
            }
            let (index, connect) = (sru.index, sru.connect.as_str());
            match sru.bonds.as_slice() {
                [b0, b1] => {
                    tag_atoms(
                        &mut mol,
                        *b0,
                        &sgatoms,
                        index,
                        HEADMARKER,
                        TAILMARKER_FRAME,
                        connect,
                    )?;
                    tag_atoms(
                        &mut mol,
                        *b1,
                        &sgatoms,
                        index,
                        TAILMARKER,
                        HEADMARKER_FRAME,
                        connect,
                    )?;
                }
                [b0, b1, b2, b3] => {
                    let xbcorr = if sru.xbcorr.len() == 4 {
                        sru.xbcorr.clone()
                    } else {
                        vec![*b0, *b2, *b1, *b3]
                    };
                    tag_atoms(
                        &mut mol,
                        xbcorr[0],
                        &sgatoms,
                        index,
                        HEADMARKER,
                        TAILMARKER_FRAME,
                        connect,
                    )?;
                    tag_atoms(
                        &mut mol,
                        xbcorr[2],
                        &sgatoms,
                        index + LADDER_OFFSET,
                        HEADMARKER,
                        TAILMARKER_FRAME,
                        connect,
                    )?;
                    tag_atoms(
                        &mut mol,
                        xbcorr[1],
                        &sgatoms,
                        index,
                        TAILMARKER,
                        HEADMARKER_FRAME,
                        connect,
                    )?;
                    tag_atoms(
                        &mut mol,
                        xbcorr[3],
                        &sgatoms,
                        index + LADDER_OFFSET,
                        TAILMARKER,
                        HEADMARKER_FRAME,
                        connect,
                    )?;
                }
                _ => {
                    return Err(RDKitError::Value(
                        "can only handle SRUs with two or four bonds".into(),
                    ));
                }
            }
        }
        if !srus.is_empty() && mol.has_prop(cp::molFileLinkNodes) {
            return Err(RDKitError::Value(
                "cannot enumerate molecules which include both SRUs and LINKNODEs".into(),
            ));
        }
        let mut frame = mol.clone();
        frame.begin_batch_edit()?;
        for sru in &srus {
            let mut repeat = mol.clone();
            let mut sgatoms = vec![false; mol.num_atoms() as usize];
            for &aidx in &sru.atoms {
                sgatoms[aidx as usize] = true;
                frame.remove_atom(aidx)?;
            }
            repeat.begin_batch_edit()?;
            for aidx in 0..repeat.num_atoms() {
                if !sgatoms[aidx as usize] {
                    repeat.remove_atom(aidx)?;
                } else {
                    repeat.atom_mut(aidx).set_prop(POLYMARKER, 1u32, false)?;
                }
            }
            repeat.commit_batch_edit()?;
            self.repeats.push(repeat);
        }
        frame.commit_batch_edit()?;
        self.mol = Some(mol);
        self.frame = Some(frame);
        Ok(())
    }

    fn variation_counts(&self) -> Vec<usize> {
        self.count_at_each_point.clone()
    }

    fn apply(&self, which: &[usize]) -> Result<ROMol> {
        if self.mol.is_none() {
            return Err(RDKitError::precondition("no molecule"));
        }
        let Some(frame) = &self.frame else {
            return Err(RDKitError::precondition("not initialized"));
        };
        if which.len() != self.count_at_each_point.len() {
            return Err(RDKitError::Value(
                "bad element choice in enumeration".into(),
            ));
        }
        for (w, count) in which.iter().zip(&self.count_at_each_point) {
            if w >= count {
                return Err(RDKitError::Value("bad element value in enumeration".into()));
            }
        }
        let mut res = ROMol::new();
        res.insert_mol(frame)?;
        let mut head_map = BTreeMap::new();
        let mut tail_map = BTreeMap::new();
        construct_head_and_tail_maps(&res, &mut head_map, &mut tail_map)?;
        for (i, &w) in which.iter().enumerate() {
            let mut filling = ROMol::new();
            let offset = self.min_repeat_counts[i];
            for iter in 0..offset + w {
                let orig_atom_count = filling.num_atoms();
                filling.insert_mol(&self.repeats[i])?;
                if iter % 2 == 1 {
                    flip_head_head_groups(orig_atom_count, &mut filling)?;
                }
                if iter > 0 {
                    connect_repeat_atoms_and_remove_extras(orig_atom_count, &mut filling)?;
                }
            }
            let n_orig_atoms = res.num_atoms();
            res.insert_mol(&filling)?;
            connect_repeat_to_frame(n_orig_atoms, &mut res, &mut head_map, &mut tail_map)?;
        }
        for (&key, &head_atom) in &head_map {
            if let Some(&tail_atom) = tail_map.get(&key) {
                res.add_bond(head_atom, tail_atom, BondType::Single)?;
                res.atom_mut(head_atom).clear_prop(HEADMARKER_FRAME);
                res.atom_mut(tail_atom).clear_prop(TAILMARKER_FRAME);
            }
        }
        Ok(res)
    }

    fn copy(&self) -> Box<dyn MolEnumeratorOp> {
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::parse_repeat_counts;

    #[test]
    fn repeat_counts_from_labels() {
        assert_eq!(parse_repeat_counts(""), (0, 4));
        assert_eq!(parse_repeat_counts("3"), (0, 4));
        assert_eq!(parse_repeat_counts("2-5"), (2, 6));
        assert_eq!(parse_repeat_counts("5-2"), (0, 4));
        assert_eq!(parse_repeat_counts("n"), (0, 4));
        assert_eq!(parse_repeat_counts("2-x"), (0, 4));
    }
}
