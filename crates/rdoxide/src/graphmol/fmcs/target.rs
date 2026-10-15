//! Targets of the search (upstream `FMCS/Target.h`, `FMCS/MatchTable.h`) and
//! the custom substructure match between graphs (upstream
//! `FMCS/SubstructMatchCustom.{h,cpp}`).

use super::final_match::final_match_check;
use super::graph::Graph;
use super::params::MCSParameters;
use super::seed::MatchV;
use crate::Result;
use crate::graphmol::ROMol;
use crate::graphmol::substruct::vf2_all;

/// `MatchTable`: rows are the atoms (bonds) of the query molecule.
#[derive(Debug, Clone, Default)]
pub struct MatchTable {
    x_size: usize,
    y_size: usize,
    data: Vec<bool>,
}

impl MatchTable {
    pub fn resize(&mut self, rows: usize, cols: usize) {
        self.data.resize(rows * cols, false);
        self.x_size = cols;
        self.y_size = rows;
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.x_size = 0;
        self.y_size = 0;
    }

    pub fn set(&mut self, row: usize, col: usize, value: bool) {
        self.data[row * self.x_size + col] = value;
    }

    pub fn at(&self, row: usize, col: usize) -> bool {
        self.data[row * self.x_size + col]
    }
}

/// `Target`: a molecule the MCS is searched in.
#[derive(Debug)]
pub struct Target<'a> {
    pub molecule: &'a ROMol,
    pub topology: Graph,
    pub atom_match_table: MatchTable,
    pub bond_match_table: MatchTable,
}

/// How a match is checked once it is complete: with the final match checker
/// of the parameters (ring fusion, chirality, user function), or not at all.
#[derive(Clone, Copy)]
pub struct FinalCheck<'p> {
    pub params: &'p MCSParameters,
}

/// `SubstructMatchCustomTable`: the first match of `query` (a seed) in
/// `target`, where atoms and bonds match as the tables say.
pub fn substruct_match_custom_table(
    target: &Graph,
    target_mol: &ROMol,
    query: &Graph,
    query_src: &ROMol,
    atom_table: &MatchTable,
    bond_table: &MatchTable,
    check: Option<FinalCheck<'_>>,
) -> Result<Option<MatchV>> {
    if query.num_vertices() > target.num_vertices() || query.num_edges() > target.num_edges() {
        return Ok(None);
    }
    let mut results = vf2_all(
        query,
        target,
        |i, j| Ok(atom_table.at(query.atom(i) as usize, target.atom(j) as usize)),
        |i, j| Ok(bond_table.at(query.bond(i) as usize, target.bond(j) as usize)),
        |c1, c2| {
            if c1.is_empty() {
                return Ok(false); // invalid index: match failed
            }
            match check {
                Some(ck) => {
                    final_match_check(c1, c2, query_src, query, target_mol, target, ck.params)
                }
                None => Ok(true),
            }
        },
        1,
        false,
    )?;
    Ok(results.pop())
}
