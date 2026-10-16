//! Position variation bonds: `ENDPTS`/`ATTACH=ANY` bonds from a dummy atom to
//! a substituent that may be attached to any of the end points (upstream
//! `MolEnumerator/PositionVariation.cpp`).

use super::{IDX_PROP_NAME, MolEnumeratorOp, preserve_orig_indices};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::{RDProps, common_properties as cp};
use crate::{RDKitError, Result, rd_log_warning};

/// `PositionVariationOp`
#[derive(Debug, Clone, Default)]
pub struct PositionVariationOp {
    mol: Option<ROMol>,
    /// (atom to attach, the atoms it may be attached to)
    variation_points: Vec<(u32, Vec<u32>)>,
    dummies_at_each_point: Vec<u32>,
}

/// `ParseV3000Array<unsigned>(text, maxV, strict = false)`: problems are
/// warnings, and the values read so far are returned.
fn parse_v3000_array(text: &str, max_count: u32) -> Vec<u32> {
    let mut rest = text;
    if !rest.starts_with('(') {
        rd_log_warning!("WARNING: first character of V3000 array is not '('\n");
    } else {
        rest = &rest[1..];
    }
    let mut tokens = rest.split_whitespace();
    let Some(count) = tokens.next().and_then(|t| t.parse::<u32>().ok()) else {
        rd_log_warning!("invalid count value in V3000 array\n");
        return Vec::new();
    };
    if count > max_count {
        rd_log_warning!("invalid count value\n");
        return Vec::new();
    }
    let mut values = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let Some(value) = tokens
            .next()
            .map(|t| t.trim_end_matches(')'))
            .and_then(|t| t.parse::<u32>().ok())
        else {
            rd_log_warning!("invalid value in V3000 array\n");
            return values;
        };
        values.push(value);
    }
    values
}

impl PositionVariationOp {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_mol(mol: &ROMol) -> Result<Self> {
        let mut op = Self::new();
        op.init_from_mol(mol)?;
        Ok(op)
    }
}

impl MolEnumeratorOp for PositionVariationOp {
    fn init_from_mol(&mut self, mol: &ROMol) -> Result<()> {
        let mut mol = mol.clone();
        self.variation_points.clear();
        if !mol.has_prop(IDX_PROP_NAME) {
            preserve_orig_indices(&mut mol)?;
        }
        for bond_idx in 0..mol.num_bonds() {
            let bond = mol.bond(bond_idx);
            let (Some(endpts), Some(attach)) = (
                bond.get_prop(cp::_MolFileBondEndPts)
                    .map(|v| v.to_rd_string()),
                bond.get_prop(cp::_MolFileBondAttach)
                    .map(|v| v.to_rd_string()),
            ) else {
                continue;
            };
            if attach != "ANY" {
                continue;
            }
            let (begin, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
            let mut atom = begin;
            if mol.atom(atom).atomic_num() == 0 {
                atom = end;
                if mol.atom(atom).atomic_num() == 0
                    && mol
                        .atom(atom)
                        .query()
                        .is_some_and(|q| q.description() == "AtomNull")
                    && mol
                        .atom(begin)
                        .query()
                        .is_some_and(|q| q.description() != "AtomNull")
                {
                    atom = begin;
                }
            }
            self.dummies_at_each_point.push(bond.other_atom_idx(atom)?);
            let mut oats = parse_v3000_array(&endpts, mol.num_atoms());
            for oat in &mut oats {
                if *oat == 0 || *oat > mol.num_atoms() {
                    return Err(RDKitError::Value("Bad variation point index".into()));
                }
                *oat -= 1;
                let attach_atom = mol.atom_mut(*oat);
                if attach_atom.is_aromatic() && attach_atom.atomic_num() != 6 {
                    attach_atom.set_num_explicit_hs(0);
                }
            }
            self.variation_points.push((atom, oats));
        }
        self.mol = Some(mol);
        Ok(())
    }

    fn variation_counts(&self) -> Vec<usize> {
        self.variation_points.iter().map(|p| p.1.len()).collect()
    }

    fn apply(&self, which: &[usize]) -> Result<ROMol> {
        let Some(mol) = &self.mol else {
            return Err(RDKitError::precondition("no molecule"));
        };
        if which.len() != self.variation_points.len() {
            return Err(RDKitError::Value(
                "bad element choice in enumeration".into(),
            ));
        }
        for (w, point) in which.iter().zip(&self.variation_points) {
            if *w >= point.1.len() {
                return Err(RDKitError::Value("bad element value in enumeration".into()));
            }
        }
        let mut res = mol.clone();
        for (w, (beg, ends)) in which.iter().zip(&self.variation_points) {
            let end = ends[*w];
            if res.bond_idx_between(*beg, end).is_some() {
                continue;
            }
            res.add_bond(*beg, end, crate::graphmol::bond::BondType::Single)?;
        }
        res.begin_batch_edit()?;
        for &idx in &self.dummies_at_each_point {
            res.remove_atom(idx)?;
        }
        res.commit_batch_edit()?;
        Ok(res)
    }

    fn copy(&self) -> Box<dyn MolEnumeratorOp> {
        // like upstream's copy constructor, the dummies are not copied
        Box::new(PositionVariationOp {
            mol: self.mol.clone(),
            variation_points: self.variation_points.clone(),
            dummies_at_each_point: Vec::new(),
        })
    }
}
