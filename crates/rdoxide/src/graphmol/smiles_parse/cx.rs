//! Reading of CXSMILES extensions (`C[C@H](F)Cl |o1:1,$_R1;;;$|`).
//!
//! Supported: coordinates, atom labels and values, atom properties, radicals,
//! enhanced stereo (stereo groups), wedged and wiggly bonds, cis/trans and
//! unspecified double bonds, dative/hydrogen/zero-order bonds, and the query
//! features (unsaturation, ring bond count, substitution count, the special
//! atom labels), substance groups (data `SgD:`, polymer `Sg:` and the
//! hierarchy `SgH:`) and link nodes. Variable attachment points are not
//! ported yet; extensions that use them fail to parse (the warning names
//! the feature).
//!
//! Ported from RDKit `Code/GraphMol/SmilesParse/CXSmilesOps.cpp` (2026.09.1).
//!
//! Copyright (C) 2016-2026 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::geometry::Point3D;
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::chirality::bond_dirs::set_stereo_for_bond;
use crate::graphmol::conformer::{Conformer, has_non_zero_z_coords};
use crate::graphmol::queries::{
    AtomFn, CompositeQueryType, Query, QueryKind, make_a_atom_query, make_ah_atom_query,
    make_atom_non_hydrogen_degree_query, make_atom_null_query, make_atom_ring_bond_count_query,
    make_atom_unsaturated_query, make_m_atom_query, make_mh_atom_query, make_q_atom_query,
    make_qh_atom_query, make_x_atom_query, make_xh_atom_query,
};
use crate::graphmol::query_ops::{atom_expand_query, replace_atom_with_query_atom};
use crate::graphmol::romol::RWMol;
use crate::graphmol::stereo_group::{StereoGroup, StereoGroupType};
use crate::graphmol::substance_group::{SubstanceGroup, add_substance_group};
use crate::rd_log_warning;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};

/// Marks "any value" in the ring bond and substitution queries (`*`); the
/// query is completed once the molecule's rings are known.
pub(crate) const ANY_VALUE: u32 = 0xDEADBEEF;

const PSEUDOATOMS_P: [&str; 2] = ["Pol_p", "Mod_p"];

/// `sgroupTypemap`: CXSMILES polymer group codes and their molfile types.
pub(crate) const SGROUP_TYPEMAP: [(&str, &str); 15] = [
    ("n", "SRU"),
    ("mon", "MON"),
    ("mer", "MER"),
    ("co", "COP"),
    ("xl", "CRO"),
    ("mod", "MOD"),
    ("mix", "MIX"),
    ("f", "FOR"),
    ("any", "ANY"),
    ("gen", "GEN"),
    ("c", "COM"),
    ("grf", "GRA"),
    ("alt", "COP"),
    ("ran", "COP"),
    ("blk", "COP"),
];

/// Property holding the position of a group in the CXSMILES (`cxsmilesindex`).
const CXSMILES_INDEX: &str = "_cxsmilesindex";
const HEAD_CROSSINGS: &str = "_headCrossings";
const TAIL_CROSSINGS: &str = "_tailCrossings";
const LABELS_PROCESSED: &str = "_cxsmilesLabelsProcessed";

/// `SmilesParseOps::detail::_needsDetectBondStereo`
pub(crate) const NEEDS_DETECT_BOND_STEREO: &str = "_needsDetectBondStereo";
/// `SmilesParseOps::detail::_needsDetectAtomStereo`
pub(crate) const NEEDS_DETECT_ATOM_STEREO: &str = "_needsDetectAtomStereo";

/// The CXSMILES extension could not be parsed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxParseError;

/// `std::atof`: the longest numeric prefix of `s`, or 0.
fn atof(s: &str) -> f64 {
    let s = s.trim_start();
    let b = s.as_bytes();
    let mut end = 0;
    let digits = |from: usize| b[from..].iter().take_while(|c| c.is_ascii_digit()).count();
    if matches!(b.first(), Some(b'+' | b'-')) {
        end += 1;
    }
    end += digits(end);
    if b.get(end) == Some(&b'.') {
        end += 1 + digits(end + 1);
    }
    if matches!(b.get(end), Some(b'e' | b'E')) {
        let mut exp = end + 1;
        if matches!(b.get(exp), Some(b'+' | b'-')) {
            exp += 1;
        }
        let n = digits(exp);
        if n > 0 {
            end = exp + n;
        }
    }
    s[..end].parse().unwrap_or(0.0)
}

/// A cursor over the extension text. Like the C++ iterators it may sit at
/// the end, where it reads as a NUL byte.
struct CxParser<'a> {
    text: &'a [u8],
    first: usize,
    mol: &'a mut RWMol,
    /// stereo groups seen so far, as `10 * id + type` (Github #6050)
    sg_tracker: Vec<u32>,
}

impl CxParser<'_> {
    fn last(&self) -> usize {
        self.text.len()
    }

    fn ch(&self) -> u8 {
        self.ch_at(0)
    }

    fn ch_at(&self, offset: usize) -> u8 {
        self.text.get(self.first + offset).copied().unwrap_or(0)
    }

    /// `first >= last`
    fn at_end(&self) -> bool {
        self.first >= self.last()
    }

    fn valid_atom(&self, idx: u32) -> bool {
        idx < self.mol.num_atoms()
    }

    fn valid_bond(&self, idx: u32) -> bool {
        idx < self.mol.num_bonds()
    }

    fn read_int(&mut self) -> Option<u32> {
        let start = self.first;
        while self.ch().is_ascii_digit() {
            self.first += 1;
        }
        if start == self.first {
            return None;
        }
        let digits = std::str::from_utf8(&self.text[start..self.first]).ok()?;
        Some(digits.parse().unwrap_or(u32::MAX))
    }

    /// `read_text_to`: the text up to one of `delims`, with `&#NN;` escapes
    /// resolved.
    fn read_text_to(&mut self, delims: &[u8]) -> Option<String> {
        let mut res = Vec::new();
        while self.first < self.last() && !delims.contains(&self.ch()) {
            if self.ch() == b'&' && self.last() - self.first > 2 && self.ch_at(1) == b'#' {
                let mut next = self.first + 2;
                while next < self.last() && self.text[next].is_ascii_digit() {
                    next += 1;
                }
                if next >= self.last() || self.text[next] != b';' {
                    return None; // quoted block not terminated with ';'
                }
                if next > self.first + 2 {
                    let code = std::str::from_utf8(&self.text[self.first + 2..next]).ok()?;
                    res.push(code.parse::<u32>().unwrap_or(0) as u8);
                }
                self.first = next + 1;
            } else {
                res.push(self.ch());
                self.first += 1;
            }
        }
        Some(String::from_utf8_lossy(&res).into_owned())
    }

    /// `get_bond_with_smiles_idx`: ring-closure bonds are appended after the
    /// ordinary bonds and carry their position in the SMILES.
    fn bond_with_smiles_idx(&self, idx: u32) -> Option<u32> {
        let mut earlier_ring_bonds = 0;
        for bond in self.mol.bonds() {
            if let Some(smiles_idx) = bond
                .props
                .get("_cxsmilesBondIdx")
                .and_then(RDValue::as_uint)
            {
                if smiles_idx == idx {
                    return Some(bond.idx());
                }
                earlier_ring_bonds += u32::from(smiles_idx < idx);
            }
        }
        let bond_idx = idx - earlier_ring_bonds;
        (bond_idx < self.mol.num_bonds()).then_some(bond_idx)
    }

    fn parse_atom_values(&mut self) -> bool {
        if self.at_end() || self.ch() != b':' {
            return false;
        }
        self.first += 1;
        let mut atom_idx = 0;
        while self.first <= self.last() && self.ch() != b'$' {
            let Some(token) = self.read_text_to(b";$") else {
                return false;
            };
            if !token.is_empty() && self.valid_atom(atom_idx) {
                self.mol
                    .atom_mut(atom_idx)
                    .props
                    .set_val(cp::molFileValue, token);
            }
            atom_idx += 1;
            if self.first <= self.last() && self.ch() != b'$' {
                self.first += 1;
            }
        }
        self.expect_dollar()
    }

    fn expect_dollar(&mut self) -> bool {
        if self.at_end() || self.ch() != b'$' {
            return false;
        }
        self.first += 1;
        true
    }

    fn parse_atom_props(&mut self) -> bool {
        if self.at_end() {
            return false;
        }
        while self.first <= self.last() && self.ch() != b'|' && self.ch() != b',' {
            if let Some(atom_idx) = self.read_int() {
                if self.at_end() || self.ch() != b'.' {
                    return false;
                }
                self.first += 1;
                let Some(name) = self.read_text_to(b".") else {
                    return false;
                };
                if !name.is_empty() {
                    if self.at_end() || self.ch() != b'.' {
                        return false;
                    }
                    self.first += 1;
                    let Some(value) = self.read_text_to(b":|,") else {
                        return false;
                    };
                    if self.valid_atom(atom_idx) && !value.is_empty() {
                        self.mol.atom_mut(atom_idx).props.set_val(&name, value);
                    }
                }
            }
            if self.first <= self.last() && self.ch() != b'|' && self.ch() != b',' {
                self.first += 1;
            }
        }
        if self.first > self.last() || (self.ch() != b'|' && self.ch() != b',') {
            return false;
        }
        if self.ch() != b'|' {
            self.first += 1;
        }
        true
    }

    fn parse_atom_labels(&mut self) -> bool {
        if self.at_end() || self.ch() != b'$' {
            return false;
        }
        self.first += 1;
        let mut atom_idx = 0;
        while self.first <= self.last() && self.ch() != b'$' {
            let Some(token) = self.read_text_to(b";$") else {
                return false;
            };
            if !token.is_empty() && self.valid_atom(atom_idx) {
                self.mol
                    .atom_mut(atom_idx)
                    .props
                    .set_val(cp::atomLabel, token);
            }
            atom_idx += 1;
            if self.first <= self.last() && self.ch() != b'$' {
                self.first += 1;
            }
        }
        self.expect_dollar()
    }

    fn parse_coords(&mut self, conf_idx: u32) -> bool {
        if self.at_end() || self.ch() != b'(' {
            return false;
        }
        let mut conf = Conformer::new(self.mol.num_atoms() as usize);
        conf.set_id(conf_idx);
        if self.mol.add_conformer(conf, false).is_err() {
            return false;
        }
        self.first += 1;
        let mut atom_idx = 0;
        let mut is_3d = false;
        let mut positions = Vec::new();
        while self.first <= self.last() && self.ch() != b')' {
            let Some(token) = self.read_text_to(b";)") else {
                return false;
            };
            if self.valid_atom(atom_idx) {
                let mut pt = Point3D::default();
                let coords: Vec<&str> = token.split(',').collect();
                if !token.is_empty() {
                    if let Some(x) = coords.first().filter(|c| !c.is_empty()) {
                        pt.x = atof(x);
                    }
                    if let Some(y) = coords.get(1).filter(|c| !c.is_empty()) {
                        pt.y = atof(y);
                    }
                    if let Some(z) = coords.get(2).filter(|c| !c.is_empty()) {
                        pt.z = atof(z);
                        is_3d = true;
                    }
                }
                positions.push((atom_idx, pt));
            }
            atom_idx += 1;
            if self.first <= self.last() && self.ch() != b')' {
                self.first += 1;
            }
        }
        let conf = self
            .mol
            .conformers_mut()
            .last_mut()
            .expect("conformer was just added");
        for (idx, pt) in positions {
            conf.positions[idx as usize] = pt;
        }
        // make sure that the conformer really is 3D
        let really_3d = is_3d && has_non_zero_z_coords(conf);
        conf.set_3d(really_3d);
        if self.at_end() || self.ch() != b')' {
            return false;
        }
        self.first += 1;
        true
    }

    fn parse_coordinate_bonds(&mut self, bond_type: BondType) -> bool {
        if self.at_end() || !matches!(self.ch(), b'C' | b'H') {
            return false;
        }
        self.first += 1;
        if self.at_end() || self.ch() != b':' {
            return false;
        }
        self.first += 1;
        while self.ch().is_ascii_digit() {
            let Some((atom_idx, bond_idx)) = self.read_int_pair() else {
                return false;
            };
            if self.valid_atom(atom_idx) && self.valid_bond(bond_idx) {
                let Some(bond) = self.bond_with_smiles_idx(bond_idx) else {
                    rd_log_warning!("BOND NOT FOUND! {bond_idx} involving atom {atom_idx}");
                    return false;
                };
                let b = self.mol.bond_mut(bond);
                if b.begin_atom_idx() != atom_idx && b.end_atom_idx() != atom_idx {
                    rd_log_warning!("BOND NOT FOUND! {bond_idx} involving atom {atom_idx}");
                    return false;
                }
                b.set_bond_type(bond_type);
                if b.begin_atom_idx() != atom_idx {
                    let begin = b.begin_atom_idx();
                    b.set_begin_atom_idx_unchecked(atom_idx);
                    b.set_end_atom_idx_unchecked(begin);
                }
            }
            if self.first < self.last() && self.ch() == b',' {
                self.first += 1;
            }
        }
        true
    }

    /// `read_int_pair(first, last, n1, n2, '.')`
    fn read_int_pair(&mut self) -> Option<(u32, u32)> {
        let n1 = self.read_int()?;
        if self.at_end() || self.ch() != b'.' {
            return None;
        }
        self.first += 1;
        Some((n1, self.read_int()?))
    }

    fn parse_zero_bonds(&mut self) -> bool {
        // these look like: C1CCCCC~CCCC1 |Z:5|
        if self.at_end() || self.ch() != b'Z' {
            return false;
        }
        self.first += 1;
        if self.at_end() || self.ch() != b':' {
            return false;
        }
        self.first += 1;
        while self.first < self.last() && self.ch().is_ascii_digit() {
            let Some(bond_idx) = self.read_int() else {
                return false;
            };
            if self.valid_bond(bond_idx) {
                let Some(bond) = self.bond_with_smiles_idx(bond_idx) else {
                    rd_log_warning!("bond {bond_idx} not found, cannot mark as zero order bond.");
                    return false;
                };
                self.mol.bond_mut(bond).set_bond_type(BondType::Zero);
            }
            if self.first < self.last() && self.ch() == b',' {
                self.first += 1;
            }
        }
        true
    }

    /// Adds `query` to the (query-ified) atom with the given index.
    fn expand_atom_query(&mut self, idx: u32, query: crate::graphmol::queries::AtomQuery) -> bool {
        if !self.mol.atom(idx).has_query() && replace_atom_with_query_atom(self.mol, idx).is_err() {
            return false;
        }
        atom_expand_query(self.mol.atom_mut(idx), query, CompositeQueryType::And, true).is_ok()
    }

    fn parse_unsaturation(&mut self) -> bool {
        if self.first + 1 >= self.last() || self.ch() != b'u' {
            return false;
        }
        self.first += 1;
        if self.at_end() || self.ch() != b':' {
            return false;
        }
        self.first += 1;
        while self.first < self.last() && self.ch().is_ascii_digit() {
            let Some(idx) = self.read_int() else {
                return false;
            };
            if self.valid_atom(idx) && !self.expand_atom_query(idx, make_atom_unsaturated_query()) {
                return false;
            }
            if self.first < self.last() && self.ch() == b',' {
                self.first += 1;
            }
        }
        true
    }

    /// `n1:n2` where `n2` is a number or `*`; returns `(n1, n2)`.
    fn read_atom_value_pair(&mut self) -> Option<(u32, u32)> {
        let n1 = self.read_int()?;
        // at least two more characters must be readable
        if self.first + 1 >= self.last() || self.ch() != b':' {
            return None;
        }
        self.first += 1;
        if self.ch() == b'*' {
            self.first += 1;
            if self.valid_atom(n1) {
                self.mol.props.set_val(cp::_NeedsQueryScan, 1i32);
            }
            Some((n1, ANY_VALUE))
        } else {
            Some((n1, self.read_int()?))
        }
    }

    fn parse_ring_bonds(&mut self) -> bool {
        if self.first >= self.last()
            || self.ch() != b'r'
            || self.first + 1 >= self.last()
            || self.ch_at(1) != b'b'
            || self.first + 2 >= self.last()
            || self.ch_at(2) != b':'
        {
            return false;
        }
        self.first += 3;
        while self.first < self.last() && self.ch().is_ascii_digit() {
            let Some((n1, n2)) = self.read_atom_value_pair() else {
                return false;
            };
            let mut greater_than = false;
            if n2 != ANY_VALUE {
                match n2 {
                    0 | 2 | 3 => {}
                    4 => greater_than = true,
                    _ => {
                        rd_log_warning!("unrecognized rb value: {n2}");
                        return false;
                    }
                }
            }
            if self.valid_atom(n1) {
                let query = if greater_than {
                    // "4" means "4 or more": the query is <= with the roles reversed
                    Query::with(
                        QueryKind::LessEqual,
                        n2 as i32,
                        AtomFn::RingBondCount,
                        "AtomRingBondCount",
                    )
                } else {
                    make_atom_ring_bond_count_query(n2 as i32)
                };
                if !self.expand_atom_query(n1, query) {
                    return false;
                }
            }
            if self.first < self.last() && self.ch() == b',' {
                self.first += 1;
            }
        }
        true
    }

    fn parse_substitution(&mut self) -> bool {
        if self.first >= self.last()
            || self.ch() != b's'
            || self.first + 1 >= self.last()
            || self.ch_at(1) != b':'
        {
            return false;
        }
        self.first += 2;
        while self.first < self.last() && self.ch().is_ascii_digit() {
            let Some((n1, n2)) = self.read_atom_value_pair() else {
                return false;
            };
            if self.valid_atom(n1)
                && !self.expand_atom_query(n1, make_atom_non_hydrogen_degree_query(n2 as i32))
            {
                return false;
            }
            if self.first < self.last() && self.ch() == b',' {
                self.first += 1;
            }
        }
        true
    }

    fn parse_wedged_bonds(&mut self) -> bool {
        // these look like: CC(O)Cl |w:1.0|; also wD and wU for down and up
        // wedges. The stereochemistry is not set from this, but the bond
        // properties are, in case client code wants the information.
        if self.first >= self.last() || self.ch() != b'w' || self.first + 1 >= self.last() {
            return false;
        }
        self.first += 1;
        let (state, cfg) = match self.ch() {
            b':' => (BondDir::Unknown, 2u32),
            b'U' => {
                self.first += 1;
                (BondDir::BeginWedge, 1)
            }
            b'D' => {
                self.first += 1;
                (BondDir::BeginDash, 3)
            }
            _ => (BondDir::None, 0),
        };
        if state == BondDir::None
            || self.first >= self.last()
            || self.first + 1 >= self.last()
            || self.ch() != b':'
        {
            return false;
        }
        self.first += 1;
        while self.first < self.last() && self.ch().is_ascii_digit() {
            let Some(atom_idx) = self.read_int() else {
                return false;
            };
            if self.first < self.last() && self.ch() == b'.' {
                self.first += 1;
            } else {
                rd_log_warning!("improperly formatted w block");
                return false;
            }
            let Some(bond_idx) = self.read_int() else {
                return false;
            };
            if self.valid_atom(atom_idx)
                && self.valid_bond(bond_idx)
                && !self.wedge_bond(atom_idx, bond_idx, state, cfg)
            {
                return false;
            }
            if self.first < self.last() && self.ch() == b',' {
                self.first += 1;
            }
        }
        true
    }

    /// One entry of a `w` block: wedges bond `bond_idx` starting at `atom_idx`.
    fn wedge_bond(&mut self, atom_idx: u32, bond_idx: u32, state: BondDir, cfg: u32) -> bool {
        let Some(bond) = self.bond_with_smiles_idx(bond_idx) else {
            rd_log_warning!(
                "bond {bond_idx} not found, wedge from atom {atom_idx} cannot be applied."
            );
            return false;
        };
        if self.mol.bond(bond).has_prop(cp::_MolFileBondCfg) {
            rd_log_warning!("w block attempts to set wedging on bond {bond} more than once.");
            return false;
        }
        // the atom has to be the start atom of the bond
        let b = self.mol.bond_mut(bond);
        if atom_idx != b.begin_atom_idx() {
            if atom_idx != b.end_atom_idx() {
                rd_log_warning!(
                    "atom {atom_idx} is not associated with bond {bond_idx}({}-{}) in w block",
                    b.begin_atom_idx(),
                    b.end_atom_idx()
                );
                return false;
            }
            let end = b.begin_atom_idx();
            b.set_begin_atom_idx_unchecked(atom_idx);
            b.set_end_atom_idx_unchecked(end);
        }
        b.props.set_val(cp::_MolFileBondCfg, cfg);
        b.set_bond_dir(state);
        let can_have_direction = crate::graphmol::canon::can_have_direction(b.bond_type());
        if cfg == 2 && can_have_direction {
            self.mol
                .atom_mut(atom_idx)
                .set_chiral_tag(ChiralType::Unspecified);
            let _ = self.mol.set_prop(NEEDS_DETECT_BOND_STEREO, 1i32, false);
        }
        if (cfg == 1 || cfg == 3) && can_have_direction {
            let _ = self.mol.set_prop(NEEDS_DETECT_ATOM_STEREO, 1i32, false);
        }
        true
    }

    fn parse_doublebond_stereo(&mut self, stereo: BondStereo) -> bool {
        // these look like: C1CCCC/C=C/CCC1 |ctu:5|; also c and t for cis or trans
        while self.first < self.last() && self.ch() != b':' {
            self.first += 1;
        }
        if self.at_end() || self.ch() != b':' {
            return false;
        }
        self.first += 1;
        while self.first < self.last() && self.ch().is_ascii_digit() {
            let Some(bond_idx) = self.read_int() else {
                return false;
            };
            if self.valid_bond(bond_idx) {
                let Some(bond) = self.bond_with_smiles_idx(bond_idx) else {
                    rd_log_warning!(
                        "bond {bond_idx} not found, cannot mark as stereo double bond."
                    );
                    return false;
                };
                if set_stereo_for_bond(self.mol, bond, stereo, true).is_err() {
                    return false;
                }
            }
            if self.first < self.last() && self.ch() == b',' {
                self.first += 1;
            }
        }
        true
    }

    fn process_radical_section(&mut self, num_radical_electrons: u32) -> bool {
        if self.at_end() {
            return false;
        }
        self.first += 1;
        if self.at_end() || self.ch() != b':' {
            return false;
        }
        self.first += 1;
        let Some(idx) = self.read_int() else {
            return false;
        };
        self.set_radicals(idx, num_radical_electrons);
        while self.first < self.last() && self.ch() == b',' {
            self.first += 1;
            if self.first < self.last() && !self.ch().is_ascii_digit() {
                return true;
            }
            let Some(idx) = self.read_int() else {
                return false;
            };
            self.set_radicals(idx, num_radical_electrons);
        }
        self.first < self.last()
    }

    fn set_radicals(&mut self, idx: u32, count: u32) {
        if self.valid_atom(idx) {
            self.mol.atom_mut(idx).set_num_radical_electrons(count);
        }
    }

    fn parse_radicals(&mut self) -> bool {
        if self.at_end() || self.ch() != b'^' {
            return false;
        }
        while self.ch() == b'^' {
            self.first += 1;
            if self.at_end() || !(b'1'..=b'7').contains(&self.ch()) {
                return false;
            }
            let count = match self.ch() {
                b'1' => 1,
                b'2'..=b'4' => 2,
                _ => 3,
            };
            if !self.process_radical_section(count) {
                return false;
            }
        }
        true
    }

    fn parse_enhanced_stereo(&mut self) -> bool {
        let group_type = match self.ch() {
            b'o' => StereoGroupType::Or,
            b'&' => StereoGroupType::And,
            _ => StereoGroupType::Absolute,
        };
        self.first += 1;
        // OR and AND groups carry a group number
        let mut group_id = 0;
        if group_type != StereoGroupType::Absolute {
            group_id = self.read_int().unwrap_or(0);
        }
        if self.at_end() || self.ch() != b':' {
            return false;
        }
        self.first += 1;
        let mut atoms: Vec<u32> = Vec::new();
        while self.first <= self.last() && self.ch().is_ascii_digit() {
            let Some(atom_idx) = self.read_int() else {
                return false;
            };
            if self.valid_atom(atom_idx) {
                if atoms.contains(&atom_idx) {
                    rd_log_warning!(
                        "Atom {atom_idx} appears more than once in stereo group specification!"
                    );
                    return false;
                }
                atoms.push(atom_idx);
            }
            if self.first < self.last() && self.ch() == b',' {
                self.first += 1;
            }
        }
        if atoms.is_empty() {
            return true;
        }
        // groups with the same type and id are one group (Github #6050)
        let hash = 10 * group_id + group_type as u32;
        let mut groups = self.mol.stereo_groups().to_vec();
        if groups.len() != self.sg_tracker.len() {
            return false;
        }
        let group = match self.sg_tracker.iter().position(|&h| h == hash) {
            Some(index) => {
                let mut all_atoms = groups[index].atoms().to_vec();
                all_atoms.extend(atoms);
                let merged =
                    StereoGroup::new(groups[index].group_type(), all_atoms, vec![], group_id);
                (index, merged)
            }
            None => {
                self.sg_tracker.push(hash);
                (
                    groups.len(),
                    StereoGroup::new(group_type, atoms, vec![], group_id),
                )
            }
        };
        let (index, Ok(group)) = group else {
            return false;
        };
        if index == groups.len() {
            groups.push(group);
        } else {
            groups[index] = group;
        }
        self.mol.set_stereo_groups(groups).is_ok()
    }

    /// `read_int_list(first, last, res, sep)`: never fails.
    fn read_int_list(&mut self, sep: u8) -> Vec<u32> {
        let mut res = Vec::new();
        loop {
            let start = self.first;
            while self.first <= self.last() && self.ch().is_ascii_digit() {
                self.first += 1;
            }
            if self.first > start {
                let digits = std::str::from_utf8(&self.text[start..self.first]).unwrap_or("0");
                res.push(digits.parse().unwrap_or(0));
            }
            if self.at_end() || self.ch() != sep {
                break;
            }
            self.first += 1;
        }
        res
    }

    /// `parse_data_sgroup_attr`: one `:`-terminated attribute of a data group.
    fn parse_data_sgroup_attr(
        &mut self,
        sgroup: &mut SubstanceGroup,
        keep: bool,
        field: &str,
        is_array: bool,
    ) -> bool {
        if self.at_end() {
            return false;
        }
        if self.ch() != b'|' {
            let Some(data) = self.read_text_to(b":") else {
                return false;
            };
            self.first += 1;
            if !data.is_empty() && keep {
                let value = if is_array {
                    RDValue::VecString(vec![data])
                } else {
                    RDValue::String(data)
                };
                sgroup.dict_mut().set_val(field, value);
            }
        }
        true
    }

    /// `|SgD:atoms:field:data:queryop:unit:tag:(coords)|`
    fn parse_data_sgroup(&mut self, n_sgroups: u32) -> bool {
        if !self.text[self.first..].starts_with(b"SgD:") {
            return false;
        }
        self.first += 4;
        let atoms = self.read_int_list(b',');
        let mut sgroup = SubstanceGroup::new("DAT");
        sgroup.dict_mut().set_val(CXSMILES_INDEX, n_sgroups);
        let mut keep = false;
        for idx in atoms {
            if self.valid_atom(idx) {
                keep = true;
                if sgroup.add_atom_with_idx(self.mol, idx).is_err() {
                    return false;
                }
            }
        }
        self.first += 1;
        if !self.parse_data_sgroup_attr(&mut sgroup, keep, "FIELDNAME", false) {
            return false;
        }
        if keep {
            sgroup
                .dict_mut()
                .set_val("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0");
        }
        for (field, is_array) in [
            ("DATAFIELDS", true),
            ("QUERYOP", false),
            ("FIELDINFO", false),
            ("FIELDTAG", false),
        ] {
            if !self.parse_data_sgroup_attr(&mut sgroup, keep, field, is_array) {
                return false;
            }
        }
        if !self.at_end() && self.ch() == b'(' {
            let Some(coords) = self.read_text_to(b")") else {
                return false;
            };
            self.first += 1;
            if keep {
                sgroup.dict_mut().set_val("COORDS", coords);
            }
        }
        if keep {
            // the label processing can destroy sgroup info, so do it now
            process_cx_smiles_labels(self.mol);
            let index = self.mol.substance_groups().len() as u32 + 1;
            sgroup.dict_mut().set_val("index", index);
            add_substance_group(self.mol, sgroup);
        }
        true
    }

    /// `|SgH:parent:child.child,parent:child|`
    fn parse_sgroup_hierarchy(&mut self) -> bool {
        if !self.text[self.first..].starts_with(b"SgH:") {
            return false;
        }
        self.first += 4;
        let find = |mol: &RWMol, id: u32| {
            mol.substance_groups()
                .iter()
                .position(|sg| sg.get_prop(CXSMILES_INDEX).and_then(RDValue::as_uint) == Some(id))
        };
        loop {
            let Some(mut parent_id) = self.read_int() else {
                return false;
            };
            let parent = find(self.mol, parent_id);
            if let Some(p) = parent
                && let Some(index) = self.mol.substance_groups()[p]
                    .get_prop("index")
                    .and_then(RDValue::as_uint)
            {
                parent_id = index;
            }
            if self.first > self.last() || self.ch() != b':' {
                return false;
            }
            self.first += 1;
            let children = self.read_int_list(b'.');
            if parent.is_some() {
                for child in children {
                    if child as usize >= self.mol.substance_groups().len() {
                        // "child id references non-existent SGroup"
                        return false;
                    }
                    if let Some(c) = find(self.mol, child) {
                        self.mol.substance_groups_mut()[c]
                            .dict_mut()
                            .set_val("PARENT", parent_id);
                    }
                }
            }
            if self.first <= self.last() && self.ch() == b',' {
                self.first += 1;
            } else {
                return true;
            }
        }
    }

    /// `|Sg:type:atoms:subscript:superscript:head crossings:tail crossings:|`;
    /// empty trailing fields may be missing.
    fn parse_polymer_sgroup(&mut self, n_sgroups: u32) -> bool {
        if !self.text[self.first..].starts_with(b"Sg:") {
            return false;
        }
        self.first += 3;
        let Some(type_code) = self.read_text_to(b":") else {
            return false;
        };
        self.first += 1;
        let Some((_, sgroup_type)) = SGROUP_TYPEMAP.iter().find(|(code, _)| *code == type_code)
        else {
            return false;
        };
        let mut sgroup = SubstanceGroup::new(sgroup_type);
        sgroup.dict_mut().set_val(CXSMILES_INDEX, n_sgroups);
        match type_code.as_str() {
            "alt" => sgroup.dict_mut().set_val("SUBTYPE", "ALT"),
            "ran" => sgroup.dict_mut().set_val("SUBTYPE", "RAN"),
            "blk" => sgroup.dict_mut().set_val("SUBTYPE", "BLO"),
            _ => {}
        }
        let mut keep = false;
        for idx in self.read_int_list(b',') {
            if self.valid_atom(idx) {
                if sgroup.add_atom_with_idx(self.mol, idx).is_err() {
                    return false;
                }
                keep = true;
            }
        }
        if self.first <= self.last() && self.ch() == b':' {
            self.first += 1;
            let Some(subscript) = self.read_text_to(b":|") else {
                return false;
            };
            if keep && !subscript.is_empty() {
                sgroup.dict_mut().set_val("LABEL", subscript);
            }
            if self.first <= self.last() && self.ch() == b':' {
                self.first += 1;
                let Some(superscript) = self.read_text_to(b":|,") else {
                    return false;
                };
                if keep && !superscript.is_empty() {
                    sgroup.dict_mut().set_val("CONNECT", superscript);
                }
                if self.first <= self.last() && self.ch() == b':' {
                    self.first += 1;
                    let mut head = self.read_int_list(b',');
                    if keep && !head.is_empty() {
                        keep = head.iter().all(|&b| self.valid_atom(b));
                        head.retain(|&b| self.valid_atom(b));
                        let _ =
                            sgroup.set_prop(HEAD_CROSSINGS, RDValue::VecUnsignedInt(head), true);
                    }
                    let mut tail = Vec::new();
                    if self.first <= self.last() && self.ch() == b':' {
                        self.first += 1;
                        tail = self.read_int_list(b',');
                    }
                    if keep && !tail.is_empty() {
                        keep = tail.iter().all(|&b| self.valid_atom(b));
                        tail.retain(|&b| self.valid_atom(b));
                        let _ =
                            sgroup.set_prop(TAIL_CROSSINGS, RDValue::VecUnsignedInt(tail), true);
                    }
                }
            }
        }
        if keep {
            process_cx_smiles_labels(self.mol);
            if finalize_polymer_sgroup(self.mol, &mut sgroup).is_err() {
                return false;
            }
            let index = self.mol.substance_groups().len() as u32 + 1;
            sgroup.dict_mut().set_val("index", index);
            add_substance_group(self.mol, sgroup);
        }
        true
    }

    /// `|LN:atom:min.max[.outer1.outer2],...|`: stored as the
    /// `molFileLinkNodes` property (1-based atoms, "min max 2 a n1 a n2").
    /// The outer atoms may be left out when the atom has two neighbors.
    fn parse_linknodes(&mut self) -> bool {
        if !self.text[self.first..].starts_with(b"LN:") {
            return false;
        }
        self.first += 3;
        let mut accum = String::new();
        while self.ch().is_ascii_digit() {
            let Some(atom) = self.read_int() else {
                return false;
            };
            if self.first + 1 >= self.last() || self.ch() != b':' {
                return false;
            }
            self.first += 1;
            let Some(start_reps) = self.read_int() else {
                return false;
            };
            if self.first + 1 >= self.last() || self.ch() != b'.' {
                return false;
            }
            self.first += 1;
            let Some(end_reps) = self.read_int() else {
                return false;
            };
            let (mut idx1, mut idx2) = (0, 0);
            if self.first < self.last() && self.ch() == b'.' {
                self.first += 1;
                let Some(i1) = self.read_int() else {
                    return false;
                };
                self.first += 1;
                let Some(i2) = self.read_int() else {
                    return false;
                };
                (idx1, idx2) = (i1, i2);
            } else if self.valid_atom(atom) && self.mol.atom_degree(atom) == 2 {
                let mut nbrs = self.mol.atom_neighbors(atom);
                idx1 = nbrs.next().unwrap_or(0);
                idx2 = nbrs.next().unwrap_or(0);
            } else if self.valid_atom(atom) {
                return false;
            }
            if self.first < self.last() && self.ch() == b',' {
                self.first += 1;
            }
            if self.valid_atom(atom) {
                if !accum.is_empty() {
                    accum.push('|');
                }
                accum += &format!(
                    "{start_reps} {end_reps} 2 {} {} {} {}",
                    atom + 1,
                    idx1.wrapping_add(1),
                    atom + 1,
                    idx2.wrapping_add(1)
                );
            }
        }
        if !accum.is_empty() {
            self.mol.dict_mut().set_val(cp::molFileLinkNodes, accum);
        }
        true
    }

    /// An extension that this port does not handle yet.
    fn unsupported(&self, what: &str) -> bool {
        rd_log_warning!("CXSMILES {what} are not supported by rdoxide yet");
        false
    }

    fn parse_it(&mut self) -> bool {
        if self.at_end() || self.ch() != b'|' {
            return false;
        }
        self.first += 1;
        let mut conf_index = 0;
        let mut n_sgroups = 0;
        while self.first < self.last() && self.ch() != b'|' {
            let length = self.last() - self.first;
            let ok = match self.ch() {
                b'(' => {
                    conf_index += 1;
                    self.parse_coords(conf_index - 1)
                }
                b'$' => {
                    if length > 4 && self.text[self.first + 1..].starts_with(b"_AV:") {
                        self.first += 4;
                        self.parse_atom_values()
                    } else {
                        self.parse_atom_labels()
                    }
                }
                _ if length > 9 && self.text[self.first..].starts_with(b"atomProp:") => {
                    self.first += 9;
                    self.parse_atom_props()
                }
                b'C' => self.parse_coordinate_bonds(BondType::Dative),
                b'H' => self.parse_coordinate_bonds(BondType::Hydrogen),
                b'Z' => self.parse_zero_bonds(),
                b'^' => self.parse_radicals(),
                b'a' | b'o' => self.parse_enhanced_stereo(),
                b'&' if self.first + 1 < self.last() && self.ch_at(1) != b'#' => {
                    self.parse_enhanced_stereo()
                }
                b'r' if self.first + 1 < self.last() && self.ch_at(1) == b'b' => {
                    self.parse_ring_bonds()
                }
                b'L' if self.first + 1 < self.last() && self.ch_at(1) == b'N' => {
                    self.parse_linknodes()
                }
                b'S' if self.first + 2 < self.last()
                    && self.ch_at(1) == b'g'
                    && self.ch_at(2) == b'D' =>
                {
                    n_sgroups += 1;
                    self.parse_data_sgroup(n_sgroups - 1)
                }
                b'S' if self.first + 2 < self.last()
                    && self.ch_at(1) == b'g'
                    && self.ch_at(2) == b'H' =>
                {
                    self.parse_sgroup_hierarchy()
                }
                b'S' if self.first + 1 < self.last() && self.ch_at(1) == b'g' => {
                    n_sgroups += 1;
                    self.parse_polymer_sgroup(n_sgroups - 1)
                }
                b'u' => self.parse_unsaturation(),
                b's' => self.parse_substitution(),
                b'm' => self.unsupported("variable attachment points"),
                b'w' => self.parse_wedged_bonds(),
                b'c' if self.first + 2 < self.last()
                    && self.ch_at(1) == b't'
                    && self.ch_at(2) == b'u' =>
                {
                    self.parse_doublebond_stereo(BondStereo::StereoAny)
                }
                b'c' => self.parse_doublebond_stereo(BondStereo::StereoCis),
                b't' => self.parse_doublebond_stereo(BondStereo::StereoTrans),
                _ => {
                    self.first += 1;
                    true
                }
            };
            if !ok {
                return false;
            }
        }
        if self.at_end() || self.ch() != b'|' {
            return false;
        }
        self.first += 1; // step past the last '|'
        true
    }
}

/// `setupUnmarkedPolymerSGroup`: no information about the bonds in and out of
/// the group is present, so they are found from the first and last atom.
fn setup_unmarked_polymer_sgroup(
    mol: &RWMol,
    sgroup: &SubstanceGroup,
    head: &mut Vec<u32>,
    tail: &mut Vec<u32>,
) -> Result<(), ()> {
    let atoms = sgroup.atoms();
    let (Some(&first), Some(&last)) = (atoms.first(), atoms.last()) else {
        return Err(()); // "no atoms in polymer sgroup"
    };
    for nbr in mol.atom_neighbors(first).collect::<Vec<_>>() {
        if atoms.contains(&nbr) {
            continue;
        }
        let bond = mol.bond_idx_between(first, nbr).ok_or(())?;
        // usually a head crossing; with a single atom in the group and a head
        // crossing already known, the bond may be a tail crossing
        if atoms.len() > 1 || head.is_empty() {
            head.push(bond);
        } else if atoms.len() == 1 {
            if tail.is_empty() {
                tail.push(bond);
            } else {
                rd_log_warning!(
                    " single atom polymer Sgroup has more than two bonds to external atoms. Ignoring all bonds after the first two."
                );
            }
        }
    }
    if atoms.len() > 1 {
        for nbr in mol.atom_neighbors(last).collect::<Vec<_>>() {
            if !atoms.contains(&nbr) {
                tail.push(mol.bond_idx_between(last, nbr).ok_or(())?);
            }
        }
    }
    Ok(())
}

/// `finalizePolymerSGroup`: sets up the connectivity and crossing bonds.
fn finalize_polymer_sgroup(mol: &RWMol, sgroup: &mut SubstanceGroup) -> Result<(), ()> {
    let mut flipped = false;
    let mut connect = sgroup
        .get_prop("CONNECT")
        .map_or_else(|| "EU".to_string(), RDValue::to_rd_string);
    if connect.contains(",f") {
        flipped = true;
        connect = connect.replace(",f", "");
    }
    connect = match connect.as_str() {
        "hh" => "HH".into(),
        "ht" => "HT".into(),
        "eu" => "EU".into(),
        other => {
            rd_log_warning!("unrecognized CXSMILES CONNECT value: '{other}'. Assuming 'eu'");
            "EU".into()
        }
    };
    sgroup.dict_mut().set_val("CONNECT", connect);

    let take = |sgroup: &mut SubstanceGroup, key: &str| match sgroup.get_prop(key) {
        Some(RDValue::VecUnsignedInt(v)) => {
            let v = v.clone();
            sgroup.clear_prop(key);
            v
        }
        _ => Vec::new(),
    };
    let mut head = take(sgroup, HEAD_CROSSINGS);
    let mut tail = take(sgroup, TAIL_CROSSINGS);
    if head.is_empty() && tail.is_empty() {
        setup_unmarked_polymer_sgroup(mol, sgroup, &mut head, &mut tail)?;
    }
    if head.is_empty() && tail.is_empty() {
        return Ok(()); // we tried
    }
    for &b in head.iter().chain(&tail) {
        sgroup.add_bond_with_idx(mol, b).map_err(|_| ())?;
    }
    let mut xbcorr = Vec::new();
    for i in 0..head.len().min(tail.len()) {
        let tail_idx = if flipped {
            tail[tail.len() - i - 1]
        } else {
            tail[i]
        };
        xbcorr.push(head[i]);
        xbcorr.push(tail_idx);
    }
    sgroup
        .dict_mut()
        .set_val("XBHEAD", RDValue::VecUnsignedInt(head));
    sgroup
        .dict_mut()
        .set_val("XBCORR", RDValue::VecUnsignedInt(xbcorr));
    Ok(())
}

/// `SmilesParseOps::parseCXExtensions`: applies the extension block at the
/// start of `text` to `mol`. `pos` ends up after the block (or wherever the
/// parser stopped on failure).
pub fn parse_cx_extensions(
    mol: &mut RWMol,
    text: &str,
    pos: &mut usize,
) -> Result<(), CxParseError> {
    if text.is_empty() {
        return Ok(());
    }
    if !text.starts_with('|') {
        return Err(CxParseError);
    }
    let mut parser = CxParser {
        text: text.as_bytes(),
        first: 0,
        mol,
        sg_tracker: Vec::new(),
    };
    let ok = parser.parse_it();
    *pos = parser.first;
    if !ok {
        return Err(CxParseError);
    }
    process_cx_smiles_labels(mol);
    mol.clear_prop(LABELS_PROCESSED);
    Ok(())
}

/// `addquery`: replaces atom `idx` by a query atom carrying `query`.
fn add_query(mol: &mut RWMol, idx: u32, query: crate::graphmol::queries::AtomQuery, symbol: &str) {
    let mut atom = Atom::new(0);
    atom.set_query(Some(query));
    atom.set_no_implicit(true);
    if mol.replace_atom(idx, &atom, true).is_ok() && !symbol.is_empty() {
        mol.atom_mut(idx).props.set_val(cp::atomLabel, symbol);
    }
}

/// `processCXSmilesLabels`: the special atom labels become queries.
fn process_cx_smiles_labels(mol: &mut RWMol) {
    if mol.has_prop(LABELS_PROCESSED) {
        return;
    }
    for idx in 0..mol.num_atoms() {
        let label = mol
            .atom(idx)
            .props
            .get(cp::atomLabel)
            .map(RDValue::to_rd_string);
        if let Some(symbol) = label {
            mol.atom_mut(idx).clear_prop(cp::dummyLabel);
            match symbol.as_str() {
                // "unspecified end groups" for polymers in Marvin
                "star_e" => add_query(mol, idx, make_atom_null_query(), &symbol),
                "Q_e" => add_query(mol, idx, make_q_atom_query(), &symbol),
                "QH_p" => add_query(mol, idx, make_qh_atom_query(), &symbol),
                // "any atom, including H" (so a plain `*` is "A" below)
                "AH_p" => add_query(mol, idx, make_ah_atom_query(), &symbol),
                "X_p" => add_query(mol, idx, make_x_atom_query(), &symbol),
                "XH_p" => add_query(mol, idx, make_xh_atom_query(), &symbol),
                "M_p" => add_query(mol, idx, make_m_atom_query(), &symbol),
                "MH_p" => add_query(mol, idx, make_mh_atom_query(), &symbol),
                s if PSEUDOATOMS_P.contains(&s) => {
                    let atom = mol.atom_mut(idx);
                    atom.props
                        .set_val(cp::dummyLabel, s[..s.len() - 2].to_string());
                    atom.clear_prop(cp::atomLabel);
                }
                _ => {}
            }
        } else {
            let atom = mol.atom(idx);
            if atom.atomic_num() == 0
                && !atom.has_query()
                && atom.isotope() == 0
                && atom.symbol() == "*"
            {
                add_query(mol, idx, make_a_atom_query(), "");
            }
        }
    }
    let _ = mol.set_prop(LABELS_PROCESSED, 1i32, true);
}
