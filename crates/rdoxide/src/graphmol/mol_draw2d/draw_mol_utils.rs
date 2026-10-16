//! Helpers of `DrawMol`: atom symbols and orientations, colours, highlight
//! parameters and small geometry functions (the free functions and the
//! symbol/orientation methods of `DrawMol.cpp`, RDKit 2026.09.1).
//!
//! Copyright (C) 2020-2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;
use std::sync::OnceLock;

use super::details::do_lines_intersect;
use super::helpers::{DrawColour, MolDrawOptions, OrientType, TextAlignType};
use super::string_rect::StringRect;
use crate::RDKitError;
use crate::Result;
use crate::geometry::{Point2D, Point3D, Transform3D};
use crate::graphmol::atom::Atom;
use crate::graphmol::mol_transforms::{
    compute_centroid, transform_conformer, transform_mol_substance_groups,
};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::queries::{
    get_atom_list_query_vals, is_atom_list_query, is_complex_atom_query,
};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::{RDProps, common_properties as cp};

/// `centerMolForDrawing`: moves the centroid of the conformer (and the
/// substance group geometry) to the origin.
pub fn center_mol_for_drawing(mol: &mut RWMol, conf_id: i32) -> Result<()> {
    let centroid = {
        let conf = mol.conformer(conf_id)?;
        compute_centroid(mol, conf, false, None)?
    };
    let mut tf = Transform3D::identity();
    tf.set_translation(&Point3D::new(-centroid.x, -centroid.y, -centroid.z));
    transform_conformer(mol.conformer_mut(conf_id)?, &tf);
    transform_mol_substance_groups(mol, &tf);
    Ok(())
}

/// `isLinearAtom`: a degree 2 atom whose two bonds have the same type and are
/// (almost) collinear.
pub fn is_linear_atom(mol: &ROMol, idx: u32, at_cds: &[Point2D]) -> bool {
    if mol.atom_degree(idx) != 2 {
        return false;
    }
    let at1 = at_cds[idx as usize];
    let mut bond_vecs = [Point2D::default(); 2];
    let mut bts = [crate::graphmol::bond::BondType::Unspecified; 2];
    for (i, nbr) in mol.atom_neighbors(idx).enumerate().take(2) {
        // a zero-length vector is an error that can be ignored
        if let Ok(v) = at1.direction_vector(&at_cds[nbr as usize]) {
            bond_vecs[i] = v;
        }
        bts[i] = mol
            .get_bond_between_atoms(idx, nbr)
            .ok()
            .flatten()
            .map_or(bts[i], |b| b.bond_type());
    }
    bts[0] == bts[1] && bond_vecs[0].dot_product(&bond_vecs[1]) < -0.95
}

/// `getColourByAtomicNum`: palette colour of an element, falling back to the
/// colour for -1 (undefined) and then to black.
pub fn get_colour_by_atomic_num(mut atomic_num: i32, opts: &MolDrawOptions) -> DrawColour {
    let mut res = DrawColour::new(0.0, 0.0, 0.0);
    if atomic_num == 1 && opts.no_atom_labels {
        atomic_num = 201;
    }
    if let Some(c) = opts.atom_colour_palette.get(&atomic_num) {
        res = *c;
    } else if atomic_num != -1 {
        if let Some(c) = opts.atom_colour_palette.get(&-1) {
            res = *c;
        }
    }
    res
}

/// `getHighlightBondColour`.
pub fn get_highlight_bond_colour(
    mol: &ROMol,
    bond_idx: u32,
    opts: &MolDrawOptions,
    highlight_bonds: &[i32],
    highlight_bond_map: &BTreeMap<i32, DrawColour>,
    highlight_atom_map: &BTreeMap<i32, DrawColour>,
) -> DrawColour {
    let mut col = DrawColour::new(0.0, 0.0, 0.0);
    if highlight_bonds.contains(&(bond_idx as i32)) {
        col = opts.highlight_colour;
        if let Some(c) = highlight_bond_map.get(&(bond_idx as i32)) {
            col = *c;
        } else {
            // not given explicitly: use the begin atom colour if set, else
            // the end atom colour, else the default
            let bond = mol.bond(bond_idx);
            if let Some(c) = highlight_atom_map.get(&(bond.begin_atom_idx() as i32)) {
                col = *c;
            } else if let Some(c) = highlight_atom_map.get(&(bond.end_atom_idx() as i32)) {
                col = *c;
            }
        }
    }
    col
}

/// `getHighlightBondWidth`.
pub fn get_highlight_bond_width(
    opts: &MolDrawOptions,
    bond_idx: i32,
    multipliers: Option<&BTreeMap<i32, i32>>,
) -> f64 {
    let mut bwm = opts.highlight_bond_width_multiplier;
    // without filled highlights the lines need to be narrower
    if !opts.fill_highlights {
        bwm = (bwm / 2).max(1);
    }
    if let Some(m) = multipliers.and_then(|m| m.get(&bond_idx)) {
        bwm = *m;
    }
    opts.bond_line_width * f64::from(bwm)
}

/// `calcPerpendicular`: the normalised perpendicular to the vector between
/// two points.
pub fn calc_perpendicular(cds1: Point2D, cds2: Point2D) -> Point2D {
    let bv = [cds1.x - cds2.x, cds1.y - cds2.y];
    let perp = [-bv[1], bv[0]];
    let len = (perp[0] * perp[0] + perp[1] * perp[1]).sqrt();
    Point2D::new(perp[0] / len, perp[1] / len)
}

/// `calcInnerPerpendicular`: like [`calc_perpendicular`] but pointing inside
/// the angle made between (1, 2) and (2, 3).
pub fn calc_inner_perpendicular(cds1: Point2D, cds2: Point2D, cds3: Point2D) -> Point2D {
    let mut perp = calc_perpendicular(cds1, cds2);
    let v1 = [cds1.x - cds2.x, cds1.y - cds2.y];
    let v2 = [cds2.x - cds3.x, cds2.y - cds3.y];
    let obv = [v1[0] - v2[0], v1[1] - v2[1]];
    // opposite directions: reverse
    if obv[0] * perp.x + obv[1] * perp.y < 0.0 {
        perp.x *= -1.0;
        perp.y *= -1.0;
    }
    perp
}

/// `adjustBondEndForString`: moves `move_end` so a bond to `end2` doesn't run
/// into the rectangles of a label drawn at `move_end`.
pub fn adjust_bond_end_for_string(
    end2: Point2D,
    padding: f64,
    rects: &[StringRect],
    move_end: &mut Point2D,
) {
    let label_pos = *move_end;
    for r in rects {
        let mut r = *r;
        r.trans += label_pos;
        let [tl, tr, br, bl] = r.calc_corners(padding);
        // a wide label such as C:7 can be crossed on more than one side
        for (a, b) in [(tl, tr), (tr, br), (br, bl), (bl, tl)] {
            if let Some(ip) = do_lines_intersect(*move_end, end2, a, b) {
                *move_end = ip;
            }
        }
    }
}

/// `(xmin, xmax, ymin, ymax)`
pub type Extremes = (f64, f64, f64, f64);

/// `findRectExtremes`.
pub fn find_rect_extremes(rect: &StringRect, align: TextAlignType, ext: &mut Extremes) {
    let mut this_xmax = rect.trans.x;
    let mut this_xmin = rect.trans.x;
    let mut this_ymax = rect.trans.y;
    let mut this_ymin = rect.trans.y;
    match align {
        TextAlignType::Start => this_xmax += rect.width,
        TextAlignType::End => this_xmin -= rect.width,
        TextAlignType::Middle => {
            this_xmax += rect.width / 2.0;
            this_xmin -= rect.width / 2.0;
        }
    }
    this_ymax += rect.height / 2.0;
    this_ymin -= rect.height / 2.0;
    ext.1 = ext.1.max(this_xmax);
    ext.0 = ext.0.min(this_xmin);
    ext.3 = ext.3.max(this_ymax);
    ext.2 = ext.2.min(this_ymin);
}

/// `findRadicalExtremes`.
pub fn find_radical_extremes(radicals: &[(StringRect, OrientType, i32)], ext: &mut Extremes) {
    for rad in radicals {
        find_rect_extremes(&rad.0, TextAlignType::Middle, ext);
    }
}

/// `getBondHighlightsForAtoms`: the bonds between pairs of highlighted atoms.
pub fn get_bond_highlights_for_atoms(mol: &ROMol, highlight_atoms: &[i32]) -> Vec<i32> {
    let mut res = Vec::new();
    for (i, &a) in highlight_atoms.iter().enumerate() {
        for &b in &highlight_atoms[i + 1..] {
            if let Ok(Some(bond)) = mol.get_bond_between_atoms(a as u32, b as u32) {
                res.push(bond.idx() as i32);
            }
        }
    }
    res
}

/// `areBondsTrans`: true if the vector at2->at1 points roughly the opposite
/// way to at3->at4.
pub fn are_bonds_trans(at1: Point2D, at2: Point2D, at3: Point2D, at4: Point2D) -> bool {
    (at1 - at2).dot_product(&(at4 - at3)) < 0.0
}

/// `areBondsParallel`: true if at2->at1 is collinear with at3->at4 within
/// `tol` of a dot product of 1.
pub fn are_bonds_parallel(
    at1: Point2D,
    at2: Point2D,
    at3: Point2D,
    at4: Point2D,
    tol: f64,
) -> Result<bool> {
    let v21 = at1.direction_vector(&at2)?;
    let v34 = at4.direction_vector(&at3)?;
    Ok((1.0 - v21.dot_product(&v34).abs()).abs() < tol)
}

/// `otherNeighbor`: the `nbor_num`'th neighbour of `first` that isn't
/// `second`.
pub fn other_neighbor(mol: &ROMol, first: u32, second: u32, nbor_num: usize) -> Option<u32> {
    mol.atom_neighbors(first)
        .filter(|&n| n != second)
        .nth(nbor_num)
}

/// `getAtomListText`: `[C,N,O]`, with a leading `!` when negated.
pub fn get_atom_list_text(atom: &Atom) -> Result<String> {
    let q = atom
        .query()
        .ok_or_else(|| RDKitError::precondition("no query"))?;
    if !is_atom_list_query(Some(q), atom.atomic_num()) {
        return Err(RDKitError::precondition("query is not an atom list"));
    }
    let mut res = String::new();
    if q.negate {
        res.push('!');
    }
    res.push('[');
    let mut vals = Vec::new();
    get_atom_list_query_vals(q, &mut vals)?;
    let pt = PeriodicTable::get();
    let syms: Result<Vec<&str>> = vals
        .iter()
        .map(|&v| pt.get_element_symbol(v as u32))
        .collect();
    res.push_str(&syms?.join(","));
    res.push(']');
    Ok(res)
}

fn complex_query_symbol_map() -> &'static BTreeMap<&'static str, &'static str> {
    static MAP: OnceLock<BTreeMap<&str, &str>> = OnceLock::new();
    MAP.get_or_init(|| {
        BTreeMap::from([
            ("![H]", "A"),
            ("![C,H]", "Q"),
            ("![C]", "QH"),
            ("[F,Cl,Br,I,At]", "X"),
            ("[F,Cl,Br,I,At,H]", "XH"),
            (
                "![He,B,C,N,O,F,Ne,Si,P,S,Cl,Ar,As,Se,Br,Kr,Te,I,Xe,At,Rn,H]",
                "M",
            ),
            (
                "![He,B,C,N,O,F,Ne,Si,P,S,Cl,Ar,As,Se,Br,Kr,Te,I,Xe,At,Rn]",
                "MH",
            ),
        ])
    })
}

fn complex_query_symbol_set() -> &'static BTreeSet<&'static str> {
    static SET: OnceLock<BTreeSet<&str>> = OnceLock::new();
    SET.get_or_init(|| complex_query_symbol_map().values().copied().collect())
}

/// `getComplexQueryAtomEquivalent`: A, Q, X, M (optionally with H) for the
/// atom list queries they stand for; anything else is returned unchanged.
pub fn get_complex_query_atom_equivalent(query: &str) -> String {
    complex_query_symbol_map()
        .get(query)
        .map_or_else(|| query.to_string(), |s| (*s).to_string())
}

fn has_symbol_query_type(atom: &Atom) -> bool {
    complex_query_symbol_set().contains(atom.query_type())
}

fn atom_map_number(atom: &Atom) -> Option<String> {
    atom.get_prop(cp::molAtomMapNumber)
        .map(|v| v.to_rd_string())
}

/// `DrawMol::getAtomSymbol`: the label of an atom, with XML-like markup for
/// super- and subscripts (as in `MolDrawing.py`). Empty when no label is
/// drawn.
pub fn get_atom_symbol(
    mol: &ROMol,
    idx: u32,
    at_cds: &[Point2D],
    orientation: OrientType,
    opts: &MolDrawOptions,
) -> Result<String> {
    if opts.no_atom_labels {
        return Ok(String::new());
    }
    let atom = mol.atom(idx);
    let mut symbol;
    let mut literal_symbol = true;
    let iso = atom.isotope();
    let text_prop = |key: &str| {
        atom.get_prop(key)
            .map(|v| v.to_rd_string())
            .unwrap_or_default()
    };
    if let Some(l) = opts.atom_labels.get(&(idx as i32)) {
        // specified labels trump everything
        symbol = l.clone();
    } else if atom.has_prop(cp::_displayLabel) || atom.has_prop(cp::_displayLabelW) {
        // if only one of _displayLabel and _displayLabelW is set use that; if
        // both are set use _displayLabelW for orientation W
        let mut lbl = text_prop(cp::_displayLabel);
        let lblw = text_prop(cp::_displayLabelW);
        if lbl.is_empty() {
            lbl = lblw.clone();
        }
        symbol = if orientation == OrientType::W && !lblw.is_empty() {
            lblw
        } else {
            lbl
        };
    } else if atom.has_prop(cp::atomLabel) {
        symbol = text_prop(cp::atomLabel);
    } else if opts.dummies_are_attachments && atom.atomic_num() == 0 && mol.atom_degree(idx) == 1 {
        symbol = String::new();
        literal_symbol = false;
    } else if opts.use_complex_query_atom_symbols && has_symbol_query_type(atom) {
        symbol = atom.query_type().to_string();
    } else if is_atom_list_query(atom.query(), atom.atomic_num()) {
        symbol = get_atom_list_text(atom)?;
        if opts.use_complex_query_atom_symbols {
            symbol = get_complex_query_atom_equivalent(&symbol);
        }
        if !opts.brackets_around_atom_lists && symbol.starts_with('[') {
            symbol = symbol[1..symbol.len() - 1].to_string();
        }
    } else if is_complex_atom_query(atom.query()) {
        symbol = "?".to_string();
        if let Some(n) = atom_map_number(atom) {
            symbol.push(':');
            symbol.push_str(&n);
        }
    } else if opts.atom_label_deuterium_tritium && atom.atomic_num() == 1 && (iso == 2 || iso == 3)
    {
        // the isotope is part of the symbol: no isotope label follows
        symbol = if iso == 2 { "D" } else { "T" }.to_string();
    } else {
        literal_symbol = false;
        let mut pre_text: Vec<String> = Vec::new();
        let mut post_text: Vec<String> = Vec::new();
        // the atom map is the first thing after the symbol
        if let Some(n) = atom_map_number(atom) {
            post_text.push(format!(":{n}"));
        }
        if atom.formal_charge() != 0 {
            // the charge always comes after the symbol, as a superscript
            let ichg = atom.formal_charge();
            let mut sgn = if ichg > 0 { "+" } else { "-" }.to_string();
            if ichg.abs() > 1 {
                sgn = format!("{}{sgn}", ichg.abs());
            }
            post_text.push(format!("<sup>{sgn}</sup>"));
        }
        let total_hs = mol.total_num_hs(idx, false)?;
        let mut num_h = if atom.atomic_num() == 6 && mol.atom_degree(idx) > 0 {
            0
        } else {
            total_hs
        };
        symbol = String::new();
        if opts.explicit_methyl && atom.atomic_num() == 6 && mol.atom_degree(idx) == 1 {
            symbol.push_str(&atom.symbol());
            num_h = total_hs;
        }
        if num_h > 0 && !atom.has_query() {
            // the H text comes after the atomic symbol
            let mut h = "H".to_string();
            if num_h > 1 {
                h.push_str(&format!("<sub>{num_h}</sub>"));
            }
            post_text.push(h);
        }
        if iso != 0
            && ((opts.isotope_labels && atom.atomic_num() != 0)
                || (opts.dummy_isotope_labels && atom.atomic_num() == 0))
        {
            // the isotope always comes before the symbol
            pre_text.push(format!("<sup>{iso}</sup>"));
        }
        symbol = pre_text.concat();
        // allenes need a C, but so does any atom with degree 2 and both
        // bonds in a line
        if is_linear_atom(mol, idx, at_cds)
            || atom.atomic_num() != 6
            || mol.atom_degree(idx) == 0
            || !pre_text.is_empty()
            || !post_text.is_empty()
        {
            symbol.push_str(&atom.symbol());
        }
        symbol.push_str(&post_text.concat());
    }
    if literal_symbol && !symbol.is_empty() {
        symbol = format!("<lit>{symbol}</lit>");
    }
    Ok(symbol)
}

/// `DrawMol::getAtomOrientation`: which side of the atom its label extends
/// to. Assumes `at_cds` already have y inverted into the draw frame.
pub fn get_atom_orientation(mol: &ROMol, idx: u32, at_cds: &[Point2D]) -> OrientType {
    // anything steeper than 70 degrees is vertical
    let vert_slope = (70.0 * PI / 180.0).tan();
    let atom = mol.atom(idx);
    let at1 = at_cds[idx as usize];
    let mut nbr_sum = Point2D::new(0.0, 0.0);
    for bond in mol.atom_bonds(idx) {
        if let Ok(other) = bond.other_atom_idx(idx) {
            nbr_sum += at_cds[other as usize] - at1;
        }
    }
    let degree = mol.atom_degree(idx);
    let mut orient;
    if degree > 0 {
        let mut islope = 1000.0;
        if nbr_sum.x.abs() > 1.0e-4 {
            islope = nbr_sum.y / nbr_sum.x;
        }
        if islope.abs() <= vert_slope {
            orient = if nbr_sum.x > 0.0 {
                OrientType::W
            } else {
                OrientType::E
            };
        } else {
            orient = if nbr_sum.y > 0.0 {
                OrientType::S
            } else {
                OrientType::N
            };
        }
        // atoms of single degree should be W or E, never N or S
        if orient == OrientType::N || orient == OrientType::S {
            if degree == 1 {
                if islope.abs() > vert_slope {
                    orient = OrientType::E;
                } else if nbr_sum.x > 0.0 {
                    orient = OrientType::W;
                } else {
                    orient = OrientType::E;
                }
            } else if degree == 3 {
                // a bond pointing down with S orientation (or up with N)
                // would put the H on the bond
                for bond in mol.atom_bonds(idx) {
                    let Ok(other) = bond.other_atom_idx(idx) else {
                        continue;
                    };
                    let bond_vec = at_cds[other as usize] - at1;
                    if bond_vec.x.abs() < 1.0e-16 {
                        orient = if bond_vec.y > 0.0 {
                            OrientType::S
                        } else {
                            OrientType::N
                        };
                        break;
                    }
                    let ang = (bond_vec.y / bond_vec.x).atan() * 180.0 / PI;
                    if ang > 80.0 && ang < 100.0 && orient == OrientType::S {
                        break;
                    } else if ang < -80.0 && ang > -100.0 && orient == OrientType::N {
                        break;
                    }
                }
            }
        }
    } else {
        // degree zero atoms from the last three periods list the Hs first
        const HS_LISTED_FIRST: [u32; 10] = [8, 9, 16, 17, 34, 35, 52, 53, 84, 85];
        orient = if HS_LISTED_FIRST.contains(&atom.atomic_num()) {
            OrientType::W
        } else {
            OrientType::E
        };
    }
    orient
}

/// `DrawMol::getAtomSymbolAndOrientation`.
pub fn get_atom_symbol_and_orientation(
    mol: &ROMol,
    idx: u32,
    at_cds: &[Point2D],
    opts: &MolDrawOptions,
) -> Result<(String, OrientType)> {
    let orient = get_atom_orientation(mol, idx, at_cds);
    let symbol = get_atom_symbol(mol, idx, at_cds, orient, opts)?;
    Ok((symbol, orient))
}
