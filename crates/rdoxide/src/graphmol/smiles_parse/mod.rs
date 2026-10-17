//! Port of RDKit `Code/GraphMol/SmilesParse` (SMILES reading).
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

pub mod canonicalize_stereo_groups;
pub mod cx;
pub mod cx_write;
pub mod lexer;
pub mod ops;
pub mod parser;
pub mod smarts;
pub mod smarts_write;
pub mod write;

use std::collections::BTreeMap;

use crate::Result;
use crate::graphmol::romol::{CI_RIGHTMOST_ATOM, RWMol};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use ops::SmilesParseError;
use parser::SmilesParser;

/// `SmilesParserParams` (v2 API).
#[derive(Debug, Clone)]
pub struct SmilesParserParams {
    pub sanitize: bool,
    pub allow_cx_smiles: bool,
    pub strict_cx_smiles: bool,
    pub parse_name: bool,
    pub remove_hs: bool,
    pub skip_cleanup: bool,
    pub debug_parse: bool,
    /// SMILES "macros", applied repeatedly until none match
    pub replacements: BTreeMap<String, String>,
}

impl Default for SmilesParserParams {
    fn default() -> Self {
        SmilesParserParams {
            sanitize: true,
            allow_cx_smiles: true,
            strict_cx_smiles: true,
            parse_name: true,
            remove_hs: true,
            skip_cleanup: false,
            debug_parse: false,
            replacements: BTreeMap::new(),
        }
    }
}

/// flex's `setup_smiles_string` trimming: bytes <= 32 at both ends.
fn trim_ctrl(s: &str) -> (usize, &str) {
    let b = s.as_bytes();
    let start = b.iter().position(|&c| c > 32).unwrap_or(b.len());
    let end = b.iter().rposition(|&c| c > 32).map_or(start, |e| e + 1);
    (start, &s[start..end.max(start)])
}

/// `preprocessSmiles`: split off name / CXSMILES part and apply replacements.
pub(crate) fn preprocess(smiles: &str, params: &SmilesParserParams) -> (String, String, String) {
    let mut lsmiles = String::new();
    let mut name = String::new();
    let mut cx_part = String::new();
    let sidx = smiles.find([' ', '\t']);
    if let Some(sidx) = sidx.filter(|&i| i != 0) {
        if params.parse_name && !params.allow_cx_smiles {
            lsmiles = smiles[..sidx].to_string();
            name = smiles[sidx..].trim().to_string();
        } else if params.allow_cx_smiles {
            lsmiles = smiles[..sidx].to_string();
            cx_part = smiles[sidx..].trim().to_string();
        }
    }
    if lsmiles.is_empty() {
        lsmiles = smiles.to_string();
    }
    if !params.replacements.is_empty() {
        loop {
            let mut again = false;
            for (k, v) in &params.replacements {
                if lsmiles.contains(k.as_str()) {
                    again = true;
                    lsmiles = lsmiles.replace(k.as_str(), v);
                }
            }
            if !again {
                break;
            }
        }
    }
    (lsmiles, name, cx_part)
}

/// Parse the SMILES grammar and run the post-parse fix-ups (`toMol`).
/// Returns `None` (after logging) on a parse error.
fn to_mol(inp: &str, orig: &str) -> Option<RWMol> {
    if inp.is_empty() {
        return Some(RWMol::new());
    }
    let (_, trimmed) = trim_ctrl(inp);
    let fail = |e: &str| {
        crate::rd_log_error!("SMILES Parse Error: {e} for input: '{orig}'");
        None
    };
    let mut parser = SmilesParser::new(trimmed);
    let parsed = parser.parse_mol();
    if !parser.branch_points.is_empty() && parsed.is_ok() {
        for &(_, pos) in &parser.branch_points {
            ops::print_syntax_error_message(trimmed, "extra open parentheses", pos, "SMILES");
        }
    }
    let mut mol = match parsed {
        Ok(m) if parser.branch_points.is_empty() => m,
        _ => return fail(&format!("Failed parsing SMILES '{inp}'")),
    };
    let post = (|| -> std::result::Result<(), SmilesParseError> {
        ops::close_mol_rings(
            &mut mol,
            &mut parser.ring_atoms,
            &mut parser.ring_bonds,
            false,
        )?;
        ops::check_chirality_specifications(&mut mol, true)?;
        ops::set_unspecified_bond_types(&mut mol);
        ops::adjust_atom_chirality_flags(&mut mol);
        Ok(())
    })();
    if let Err(SmilesParseError(e)) = post {
        return fail(&e);
    }
    if mol.has_atom_bookmark(CI_RIGHTMOST_ATOM) {
        mol.clear_atom_bookmark(CI_RIGHTMOST_ATOM);
    }
    Some(mol)
}

/// `handleCXPartAndName`: applies the CXSMILES extension (if any) and takes
/// the name from what follows it. Strict parsing turns failures into errors;
/// otherwise the text consumed so far is kept as `_CXSMILES_Data`.
pub(crate) fn handle_cx_part_and_name(
    res: &mut RWMol,
    allow_cx_smiles: bool,
    strict_cx_smiles: bool,
    parse_name: bool,
    cx_part: &str,
    name: &mut String,
) -> Result<()> {
    if cx_part.is_empty() {
        return Ok(());
    }
    let mut pos = 0;
    let mut cx_failed = false;
    if allow_cx_smiles {
        if cx_part.starts_with('|') {
            if cx::parse_cx_extensions(res, cx_part, &mut pos).is_err() {
                cx_failed = true;
                if strict_cx_smiles {
                    return Err(crate::RDKitError::Value(
                        "SmilesParseException: CXSMILES extension could not be parsed".into(),
                    ));
                }
            }
            res.props.set_val(
                "_CXSMILES_Data",
                cx_part[..pos.min(cx_part.len())].to_string(),
            );
        } else if strict_cx_smiles && !parse_name && !cx_part.is_empty() {
            return Err(crate::RDKitError::Value(
                "SmilesParseException: CXSMILES extension does not start with | and parseName=false"
                    .into(),
            ));
        }
    }
    if !cx_failed && parse_name && pos < cx_part.len() {
        *name = cx_part[pos..].trim().to_string();
    }
    Ok(())
}

/// Stereo from the conformers that the CX extension may have added: wedge
/// bonds give atom parities, a 3D conformer (without a 2D one) gives them
/// from the geometry.
pub(crate) fn prepare_conformer_stereo(res: &mut RWMol) -> Result<()> {
    let confs = res.conformers();
    let conf = confs.iter().find(|c| !c.is_3d()).map(|c| c.id());
    let conf3d = confs.iter().find(|c| c.is_3d()).map(|c| c.id());
    if res.has_prop(cx::NEEDS_DETECT_ATOM_STEREO) {
        res.clear_prop(cx::NEEDS_DETECT_ATOM_STEREO);
        crate::graphmol::chirality::from_coords::assign_chiral_types_from_bond_dirs(
            res,
            conf.map_or(0, |c| c as i32),
            true,
        )?;
    }
    if let (None, Some(conf3d)) = (conf, conf3d) {
        res.update_property_cache(false)?;
        crate::graphmol::chirality::from_coords::assign_chiral_types_from_3d(
            res,
            conf3d as i32,
            true,
        )?;
    }
    Ok(())
}

/// `MolFromSmiles(smiles, params)`. `Ok(None)` mirrors upstream's null
/// result for unparsable input; sanitization failures are returned as errors
/// (the Python layer logs them and returns `None`).
pub fn mol_from_smiles(smiles: &str, params: &SmilesParserParams) -> Result<Option<RWMol>> {
    let (lsmiles, mut name, cx_part) = preprocess(smiles, params);
    let Some(mut res) = to_mol(&lsmiles, &lsmiles) else {
        return Ok(None);
    };
    handle_cx_part_and_name(
        &mut res,
        params.allow_cx_smiles,
        params.strict_cx_smiles,
        params.parse_name,
        &cx_part,
        &mut name,
    )?;
    prepare_conformer_stereo(&mut res)?;
    if params.sanitize || params.remove_hs {
        crate::graphmol::mol_ops::smiles_post_parse_sanitize(
            &mut res,
            params.sanitize,
            params.remove_hs,
        )?;
    } else {
        crate::graphmol::mol_ops::clear_single_bond_dir_flags(&mut res, true);
    }
    if res.has_prop(cp::_NeedsQueryScan) {
        res.clear_prop(cp::_NeedsQueryScan);
        if !params.sanitize {
            // the ring bond queries need rings
            crate::graphmol::find_rings::fast_find_rings(&mut res)?;
        }
        crate::graphmol::queries::complete_mol_queries(&mut res, cx::ANY_VALUE)?;
    }
    if !params.skip_cleanup {
        ops::cleanup_after_parsing(&mut res);
    }
    if !name.is_empty() {
        res.props.set_val(cp::_Name, name);
    }
    Ok(Some(res))
}

/// `AtomFromSmiles`: the atom of a SMILES atom expression such as `[13CH3]`;
/// `None` (with a logged error) if it cannot be parsed.
pub fn atom_from_smiles(smiles: &str) -> Option<crate::graphmol::atom::Atom> {
    if smiles.is_empty() {
        return None;
    }
    SmilesParser::new(smiles).parse_atom().ok()
}

/// `BondFromSmiles`: the bond of a SMILES bond symbol such as `=`.
pub fn bond_from_smiles(smiles: &str) -> Option<crate::graphmol::bond::Bond> {
    if smiles.is_empty() {
        return None;
    }
    SmilesParser::new(smiles).parse_bond().ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::{BondType, ChiralType};

    fn parse(s: &str) -> Option<RWMol> {
        let p = SmilesParserParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        };
        mol_from_smiles(s, &p).unwrap()
    }

    #[test]
    fn basics() {
        let m = parse("CC(=O)O").unwrap();
        assert_eq!(m.num_atoms(), 4);
        assert_eq!(m.num_bonds(), 3);
        assert_eq!(m.bond(1).bond_type(), BondType::Double);
        assert_eq!(
            (m.bond(2).begin_atom_idx(), m.bond(2).end_atom_idx()),
            (1, 3)
        );
        assert!(m.atoms().iter().all(|a| !a.has_prop(cp::_SmilesStart)));
    }

    #[test]
    fn rings_and_aromatic() {
        let m = parse("c1ccccc1").unwrap();
        assert_eq!(m.num_bonds(), 6);
        assert!(
            m.bonds()
                .iter()
                .all(|b| b.bond_type() == BondType::Aromatic && b.is_aromatic())
        );
        // ring closure bond is added last and starts at the closing atom
        assert_eq!(
            (m.bond(5).begin_atom_idx(), m.bond(5).end_atom_idx()),
            (5, 0)
        );
        // closures are processed in ring-number order
        let m = parse("C2CC1CC2CC1").unwrap();
        let rb: Vec<(u32, u32)> = m.bonds()[6..]
            .iter()
            .map(|b| (b.begin_atom_idx(), b.end_atom_idx()))
            .collect();
        assert_eq!(rb, vec![(6, 2), (4, 0)]);
        assert!(parse("C1CC").is_none());
        assert!(parse("C11").is_none());
        let m = parse("C%10CC%10").unwrap();
        assert_eq!(m.num_bonds(), 3);
        let m = parse("C%(123)CC%(123)").unwrap();
        assert_eq!(m.num_bonds(), 3);
    }

    #[test]
    fn bracket_atoms() {
        let m = parse("[13CH3:7][N+](C)(C)[O-]").unwrap();
        let a = m.atom(0);
        assert_eq!(
            (
                a.isotope(),
                a.num_explicit_hs(),
                a.atom_map_num(),
                a.no_implicit()
            ),
            (13, 3, 7, true)
        );
        assert_eq!(m.atom(1).formal_charge(), 1);
        assert_eq!(m.atom(4).formal_charge(), -1);
        let m = parse("[Fe+2].[Cl-].[Cl-]").unwrap();
        assert_eq!(m.atom(0).formal_charge(), 2);
        assert_eq!(m.num_bonds(), 0);
        let m = parse("[2H][H]").unwrap();
        assert_eq!(m.atom(0).isotope(), 2);
        let m = parse("[#6]").unwrap();
        assert_eq!(m.atom(0).atomic_num(), 6);
        let m = parse("*C").unwrap();
        assert_eq!(m.atom(0).symbol(), "*");
        assert!(parse("[Fe+++]").is_none());
        assert!(parse("CX").is_none());
    }

    #[test]
    fn chirality() {
        let m = parse("F[C@H](Cl)Br").unwrap();
        assert_eq!(m.atom(1).chiral_tag(), ChiralType::TetrahedralCCW);
        // expected values checked against the official RDKit wheel
        let m = parse("[C@@H](F)(Cl)Br").unwrap();
        assert_eq!(m.atom(0).chiral_tag(), ChiralType::TetrahedralCCW);
        let m = parse("F[C@]1(Cl)CC1").unwrap();
        assert_eq!(m.atom(1).chiral_tag(), ChiralType::TetrahedralCCW);
        assert_eq!(
            (m.bond(4).begin_atom_idx(), m.bond(4).end_atom_idx()),
            (4, 1)
        );
        let m = parse("F[Pt@SP1](Cl)(F)Cl").unwrap();
        assert_eq!(m.atom(1).chiral_tag(), ChiralType::SquarePlanar);
        assert!(parse("C[C@TH0]").is_none());
        let m = parse("F[C@TH2](Cl)(Br)I").unwrap();
        assert_eq!(m.atom(1).chiral_tag(), ChiralType::TetrahedralCW);
    }

    #[test]
    fn branches_and_errors() {
        assert!(parse("CC(C").is_none());
        assert!(parse("CC)C").is_none());
        assert!(parse("C()C").is_none());
        let m = parse("C(C)(C)C").unwrap();
        assert_eq!(m.atom_degree(0), 3);
        let m = parse("C->[Fe]<-N").unwrap();
        assert_eq!(m.bond(0).bond_type(), BondType::Dative);
        assert_eq!(
            (m.bond(1).begin_atom_idx(), m.bond(1).end_atom_idx()),
            (2, 1)
        );
        assert_eq!(parse("").unwrap().num_atoms(), 0);
        let p = SmilesParserParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        };
        let m = mol_from_smiles("CCO ethanol", &p).unwrap().unwrap();
        assert_eq!(m.name(), "ethanol");
    }
}
