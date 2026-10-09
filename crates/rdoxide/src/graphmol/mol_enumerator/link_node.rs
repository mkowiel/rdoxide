//! Link nodes: `LINKNODE` records of mol files / the `LN:` CXSMILES
//! extension, stored on the molecule as the `molFileLinkNodes` property
//! (upstream `MolEnumerator/LinkNode.h`).

use crate::error::{RDKitError, Result};
use crate::graphmol::romol::ROMol;
use crate::rd_log_warning;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};

/// `MolEnumerator::LinkNode`: the group of atoms between two bonds that may
/// repeat `min_rep` to `max_rep` times.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LinkNode {
    pub min_rep: u32,
    pub max_rep: u32,
    pub n_bonds: u32,
    /// `(atom, outer atom)` of each of the two bonds, zero based
    pub bond_atoms: Vec<(u32, u32)>,
}

/// Why a record was rejected: an error when `strict`, otherwise a warning
/// and the record is skipped.
fn reject(strict: bool, err: RDKitError) -> Result<()> {
    if strict {
        Err(err)
    } else {
        rd_log_warning!("{err}");
        Ok(())
    }
}

/// `utils::getMolLinkNodes`: the link nodes described by the
/// `molFileLinkNodes` property (`|`-separated records of space-separated
/// numbers).
pub fn get_mol_link_nodes(mol: &ROMol, strict: bool) -> Result<Vec<LinkNode>> {
    let Some(text) = mol
        .get_prop(cp::molFileLinkNodes)
        .map(RDValue::to_rd_string)
    else {
        return Ok(Vec::new());
    };
    let mut res = Vec::new();
    for record in text.split('|').filter(|r| !r.is_empty()) {
        let data: std::result::Result<Vec<u32>, _> = record
            .split(' ')
            .filter(|t| !t.is_empty())
            .map(str::parse::<u32>)
            .collect();
        let Ok(mut data) = data else {
            reject(
                strict,
                RDKitError::Value(format!(
                    "Cannot convert values in LINKNODE '{record}' to unsigned ints"
                )),
            )?;
            continue;
        };
        // data[2] is the number of bonds, followed by an atom pair for each
        if data.len() < 5 || data.len() < 3 + 2 * data[2] as usize {
            reject(
                strict,
                RDKitError::Value(format!("not enough values in LINKNODE '{record}'")),
            )?;
            continue;
        }
        let (min_rep, max_rep) = (data[0], data[1]);
        if min_rep == 0 || max_rep < min_rep {
            reject(
                strict,
                RDKitError::Value(format!("bad counts in LINKNODE '{record}'")),
            )?;
            continue;
        }
        let n_bonds = data[2];
        if n_bonds != 2 {
            reject(
                strict,
                RDKitError::Runtime("only link nodes with 2 bonds are currently supported".into()),
            )?;
            continue;
        }
        // both bonds must start from the same atom
        if data[3] != data[5] {
            reject(
                strict,
                RDKitError::Value(format!(
                    "bonds don't start at the same atom for LINKNODE '{record}'"
                )),
            )?;
            continue;
        }
        for v in &mut data[3..=6] {
            *v = v.wrapping_sub(1);
        }
        let bonded = |a: u32, b: u32| {
            a < mol.num_atoms() && b < mol.num_atoms() && mol.bond_idx_between(a, b).is_some()
        };
        if !bonded(data[4], data[3]) || !bonded(data[6], data[5]) {
            reject(
                strict,
                RDKitError::Value(format!(
                    "bond not found between atoms in LINKNODE '{record}'"
                )),
            )?;
            continue;
        }
        res.push(LinkNode {
            min_rep,
            max_rep,
            n_bonds,
            bond_atoms: vec![(data[3], data[4]), (data[5], data[6])],
        });
    }
    Ok(res)
}
