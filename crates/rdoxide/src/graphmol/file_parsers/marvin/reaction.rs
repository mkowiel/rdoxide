//! Marvin reactions: `ReactionFromMrvBlock`, `ReactionFromMrvFile`,
//! `MrvBlockIsReaction`, `MrvFileIsReaction`, `ChemicalReactionToMrvBlock`
//! and `ChemicalReactionToMrvFile`.

use super::defs::{MarvinMol, NO_COORD, ROOT};
use super::emit::double_as_text;
use super::reader::{MrvParserParams, mol_from_marvin_mol};
use super::transform::NumberingCounts;
use super::writer::mol_to_marvin_mol;
use super::xml::{XmlTree, read_xml, write_xml};
use crate::graphmol::chem_reactions::ChemicalReaction;
use crate::graphmol::query_ops::replace_atom_with_query_atom;
use crate::{RDKitError, Result};

const ARROW_MIN_LENGTH: f64 = 1.0;
const ARROW_SPACE: f64 = 0.5;
const PLUS_SPACE: f64 = 1.0;

fn bad_path(path: &str) -> RDKitError {
    RDKitError::Runtime(format!("No such node ({path})"))
}

/// `MrvBlockIsReaction`
pub fn mrv_block_is_reaction(text: &str) -> Result<bool> {
    Ok(read_xml(text)?
        .get_child("cml.MDocument.MChemicalStruct.reaction")
        .is_some())
}

/// `MrvFileIsReaction`
pub fn mrv_file_is_reaction(path: &str) -> Result<bool> {
    let bytes =
        std::fs::read(path).map_err(|_| RDKitError::BadFile(format!("Bad input file {path}")))?;
    if bytes.is_empty() {
        return Ok(false);
    }
    mrv_block_is_reaction(&String::from_utf8_lossy(&bytes))
}

/// `ReactionFromMrvBlock`
pub fn reaction_from_mrv_block(text: &str, params: &MrvParserParams) -> Result<ChemicalReaction> {
    let tree = read_xml(text)?;
    tree.get_child("cml.MDocument")
        .ok_or_else(|| bad_path("cml.MDocument"))?;
    let rxn_tree = tree
        .get_child("cml.MDocument.MChemicalStruct.reaction")
        .ok_or_else(|| bad_path("cml.MDocument.MChemicalStruct.reaction"))?;

    let parse_list = |name: &str| -> Result<Vec<MarvinMol>> {
        let mut mols = Vec::new();
        if let Some(list) = rxn_tree.child(name) {
            for (key, child) in &list.children {
                if key != "<xmlattr>" {
                    mols.push(MarvinMol::parse(child)?);
                }
            }
        }
        Ok(mols)
    };
    let mut reactants = parse_list("reactantList")?;
    let mut agents = parse_list("agentList")?;
    let mut products = parse_list("productList")?;

    // `MarvinReaction::prepSgroupsForRDKit`: one numbering for all the molecules
    let mut counts = NumberingCounts::default();
    for mol in reactants
        .iter_mut()
        .chain(agents.iter_mut())
        .chain(products.iter_mut())
    {
        mol.prep_sgroups_for_rdkit(ROOT)?;
        mol.clean_up_numbering(&mut counts);
    }

    let mut rxn = ChemicalReaction::new();
    for m in &reactants {
        rxn.add_reactant_template(mol_from_marvin_mol(m, params)?);
    }
    for m in &agents {
        rxn.add_agent_template(mol_from_marvin_mol(m, params)?);
    }
    for m in &products {
        rxn.add_product_template(mol_from_marvin_mol(m, params)?);
    }
    // the templates are queries
    for mol in rxn.reactants_and_products_mut() {
        for i in 0..mol.num_atoms() {
            replace_atom_with_query_atom(mol, i)?;
        }
    }
    Ok(rxn)
}

/// `ReactionFromMrvFile`
pub fn reaction_from_mrv_file(path: &str, params: &MrvParserParams) -> Result<ChemicalReaction> {
    let bytes =
        std::fs::read(path).map_err(|_| RDKitError::BadFile(format!("Bad input file {path}")))?;
    reaction_from_mrv_block(&String::from_utf8_lossy(&bytes), params)
}

// ---------------------------------------------------------------- writing

/// `MarvinRectangle`
#[derive(Debug, Clone, Copy)]
struct Rect {
    ul: (f64, f64),
    lr: (f64, f64),
}

impl Rect {
    fn of_atoms(m: &MarvinMol) -> Rect {
        let atoms = &m.node(ROOT).atoms;
        if atoms.is_empty() {
            return Rect {
                ul: (0.0, 0.0),
                lr: (0.0, 0.0),
            };
        }
        let mut r = Rect {
            ul: (f64::MAX, -f64::MAX),
            lr: (-f64::MAX, f64::MAX),
        };
        for &a in atoms {
            let a = &m.atoms[a];
            r.ul.0 = r.ul.0.min(a.x2);
            r.lr.0 = r.lr.0.max(a.x2);
            r.ul.1 = r.ul.1.max(a.y2);
            r.lr.1 = r.lr.1.min(a.y2);
        }
        r
    }

    fn of_rects(rects: &[Rect]) -> Rect {
        let Some(first) = rects.first() else {
            return Rect {
                ul: (0.0, 0.0),
                lr: (0.0, 0.0),
            };
        };
        let mut r = Rect {
            ul: (f64::MAX, -f64::MAX),
            lr: (-f64::MAX, f64::MAX),
        };
        let _ = first;
        for o in rects {
            r.extend(o);
        }
        r
    }

    fn extend(&mut self, o: &Rect) {
        self.ul.0 = self.ul.0.min(o.ul.0);
        self.lr.0 = self.lr.0.max(o.lr.0);
        self.ul.1 = self.ul.1.max(o.ul.1);
        self.lr.1 = self.lr.1.min(o.lr.1);
    }

    fn center(&self) -> (f64, f64) {
        ((self.lr.0 + self.ul.0) / 2.0, (self.lr.1 + self.ul.1) / 2.0)
    }

    fn overlaps_vertically(&self, o: &Rect) -> bool {
        !(o.ul.1 < self.lr.1 || o.lr.1 > self.ul.1)
    }

    fn overlaps_horizontally(&self, o: &Rect) -> bool {
        !(o.ul.0 > self.lr.0 || o.lr.0 < self.ul.0)
    }
}

fn has_coords(m: &MarvinMol) -> bool {
    let atoms = &m.node(ROOT).atoms;
    let two_d = atoms
        .iter()
        .all(|&a| m.atoms[a].x2 != NO_COORD && m.atoms[a].y2 != NO_COORD);
    two_d
        || atoms.iter().all(|&a| {
            let a = &m.atoms[a];
            a.x3 != NO_COORD && a.y3 != NO_COORD && a.z3 != NO_COORD
        })
}

#[derive(Debug, Default)]
struct Plus {
    id: String,
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
}

#[derive(Debug, Default)]
struct Arrow {
    kind: String,
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
}

/// `GetArrowPerdendicularPosition`
fn arrow_perpendicular_position(mols: &[MarvinMol], vertical: bool) -> f64 {
    let mut rects: Vec<Rect> = Vec::new();
    for m in mols {
        let r = Rect::of_atoms(m);
        let overlap = rects.iter().any(|o| {
            if vertical {
                r.overlaps_vertically(o)
            } else {
                r.overlaps_horizontally(o)
            }
        });
        if !overlap {
            rects.push(r);
        }
    }
    if vertical {
        rects.sort_by(|a, b| b.center().1.total_cmp(&a.center().1));
    } else {
        rects.sort_by(|a, b| a.center().0.total_cmp(&b.center().0));
    }
    for w in rects.windows(2) {
        let (r1, r2) = (&w[0], &w[1]);
        if vertical {
            if r2.ul.1 - r1.lr.1 >= ARROW_SPACE {
                return (r2.ul.1 + r1.lr.1) / 2.0;
            }
        } else if r2.ul.0 - r1.lr.0 >= ARROW_SPACE {
            return (r2.ul.0 + r1.lr.0) / 2.0;
        }
    }
    if vertical {
        rects[0].lr.1 - ARROW_SPACE
    } else {
        rects[0].ul.0 - ARROW_SPACE
    }
}

/// `AddMarvinPluses`
fn add_pluses(pluses: &mut Vec<Plus>, mols: &[MarvinMol], count: &mut u32) {
    let mut rows: Vec<Vec<Rect>> = Vec::new();
    for m in mols {
        if !has_coords(m) {
            return; // no good way to make the signs
        }
        let r = Rect::of_atoms(m);
        let mut found = false;
        for row in rows.iter_mut() {
            if row.iter().any(|o| r.overlaps_vertically(o)) {
                found = true;
                row.push(r);
            }
        }
        if !found {
            rows.push(vec![r]);
        }
    }
    // (upstream sorts copies of the rows, which does nothing)
    rows.sort_by(|a, b| {
        Rect::of_rects(b)
            .center()
            .1
            .total_cmp(&Rect::of_rects(a).center().1)
    });

    for (ri, row) in rows.iter().enumerate() {
        for w in row.windows(2) {
            let (r1, r2) = (&w[0], &w[1]);
            let x = (r2.ul.0 + r1.lr.0) / 2.0;
            let y = if r2.lr.0 - r1.ul.0 >= PLUS_SPACE {
                (r2.center().1 + r1.center().1) / 2.0
            } else {
                r1.lr.1.min(r2.lr.1) - PLUS_SPACE / 2.0
            };
            *count += 1;
            pluses.push(Plus {
                id: format!("o{count}"),
                x1: x,
                y1: y,
                x2: x + PLUS_SPACE,
                y2: y + PLUS_SPACE,
            });
        }
        if ri > 0 {
            let prev = Rect::of_rects(&rows[ri - 1]);
            let curr = Rect::of_rects(row);
            let (x, y) = if prev.lr.1 - curr.ul.1 > PLUS_SPACE {
                (
                    (prev.center().0 + curr.center().0) / 2.0,
                    (prev.lr.1 + curr.ul.1) / 2.0,
                )
            } else {
                (row[0].ul.0, row[0].center().1)
            };
            *count += 1;
            pluses.push(Plus {
                id: format!("o{count}"),
                x1: x,
                y1: y,
                x2: x + PLUS_SPACE,
                y2: y + PLUS_SPACE,
            });
        }
    }
}

/// `SetArrow`
fn set_arrow(reactants: &[MarvinMol], agents: &[MarvinMol], products: &[MarvinMol]) -> Arrow {
    let mut arrow = Arrow::default();
    if reactants.is_empty() || products.is_empty() {
        return arrow;
    }
    if !reactants.iter().chain(products).all(has_coords) {
        return arrow;
    }
    arrow.kind = "DEFAULT".into();
    let combine = |mols: &[MarvinMol]| {
        let mut r = Rect::of_atoms(&mols[0]);
        for m in &mols[1..] {
            r.extend(&Rect::of_atoms(m));
        }
        r
    };
    let reactant = combine(reactants);
    let product = combine(products);
    let (rc, pc) = (reactant.center(), product.center());
    let need = ARROW_MIN_LENGTH + 2.0 * ARROW_SPACE;

    if product.ul.0 - reactant.lr.0 > need {
        arrow.x1 = reactant.lr.0 + ARROW_SPACE;
        arrow.x2 = product.ul.0 - ARROW_SPACE;
        arrow.y1 = if agents.is_empty() {
            (rc.1 + pc.1) / 2.0
        } else {
            arrow_perpendicular_position(agents, true)
        };
        arrow.y2 = arrow.y1;
    } else if reactant.lr.1 - product.ul.1 > need {
        arrow.x1 = if agents.is_empty() {
            (rc.0 + pc.0) / 2.0
        } else {
            arrow_perpendicular_position(agents, false)
        };
        arrow.x2 = arrow.x1;
        arrow.y1 = reactant.lr.1 - ARROW_SPACE;
        arrow.y2 = product.ul.1 + ARROW_SPACE;
    } else if ((rc.0 - pc.0).powi(2) + (rc.1 - pc.1).powi(2)).sqrt() > need {
        arrow.x1 = rc.0;
        arrow.x2 = pc.0;
        arrow.y1 = rc.1;
        arrow.y2 = pc.1;
    } else {
        // really no good place for the arrow
        arrow.x1 = reactant.lr.0 + ARROW_SPACE;
        arrow.x2 = reactant.lr.0 + ARROW_MIN_LENGTH + ARROW_SPACE;
        arrow.y1 = (rc.1 + pc.1) / 2.0;
        arrow.y2 = arrow.y1;
    }
    arrow
}

fn plus_tree(p: &Plus) -> XmlTree {
    let mut out = XmlTree::default();
    out.put_attr("id", p.id.clone());
    out.put_attr("toptions", "NOROT");
    out.put_attr("fontScale", "14.0");
    out.put_attr("halign", "CENTER");
    out.put_attr("valign", "CENTER");
    out.put_attr("autoSize", "true");
    let mut field = XmlTree::with_data("{D font=SansSerif,size=18,bold}+");
    field.put_attr("name", "text");
    out.put_child("Field", field);
    for (x, y) in [(p.x1, p.y1), (p.x2, p.y1), (p.x2, p.y2), (p.x1, p.y2)] {
        let mut pt = XmlTree::default();
        pt.put_attr("x", double_as_text(x, 6));
        pt.put_attr("y", double_as_text(y, 6));
        out.add_child("MPoint", pt);
    }
    out
}

/// `ChemicalReactionToMrvBlock`
pub fn reaction_to_mrv_block(rxn: &ChemicalReaction, pretty_print: bool) -> Result<String> {
    let mut counts = NumberingCounts::default();
    let mut convert = |mols: &[crate::graphmol::romol::ROMol]| -> Result<Vec<MarvinMol>> {
        mols.iter()
            .map(|m| mol_to_marvin_mol(&mut m.clone(), &mut counts, -1))
            .collect()
    };
    let reactants = convert(rxn.reactants())?;
    let agents = convert(rxn.agents())?;
    let products = convert(rxn.products())?;

    let mut pluses = Vec::new();
    let mut plus_count = 0;
    add_pluses(&mut pluses, &reactants, &mut plus_count);
    add_pluses(&mut pluses, &products, &mut plus_count);
    let arrow = set_arrow(&reactants, &agents, &products);

    let mut cml = XmlTree::default();
    cml.put_attr("xmlns", "http://www.chemaxon.com");
    cml.put_attr("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    cml.put_attr(
        "xsi:schemaLocation",
        "http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd",
    );
    let mut reaction = XmlTree::default();
    for (list_name, mols) in [
        ("reactantList", &reactants),
        ("agentList", &agents),
        ("productList", &products),
    ] {
        if mols.is_empty() {
            continue;
        }
        let mut list = XmlTree::default();
        for m in mols {
            list.add_child("molecule", m.node_tree(ROOT));
        }
        reaction.add_child(list_name, list);
    }
    let mut arrow_tree = XmlTree::default();
    arrow_tree.put_attr("type", arrow.kind.clone());
    arrow_tree.put_attr("x1", double_as_text(arrow.x1, 6));
    arrow_tree.put_attr("y1", double_as_text(arrow.y1, 6));
    arrow_tree.put_attr("x2", double_as_text(arrow.x2, 6));
    arrow_tree.put_attr("y2", double_as_text(arrow.y2, 6));
    reaction.put_child("arrow", arrow_tree);

    let mut chem = XmlTree::default();
    chem.add_child("reaction", reaction);
    let mut doc = XmlTree::default();
    doc.add_child("MChemicalStruct", chem);
    for p in &pluses {
        doc.add_child("MReactionSign", plus_tree(p));
    }
    cml.add_child("MDocument", doc);
    let mut root = XmlTree::default();
    root.add_child("cml", cml);
    Ok(write_xml(&root, pretty_print))
}

/// `ChemicalReactionToMrvFile`
pub fn reaction_to_mrv_file(rxn: &ChemicalReaction, path: &str, pretty_print: bool) -> Result<()> {
    let text = reaction_to_mrv_block(rxn, pretty_print)?;
    std::fs::write(path, text).map_err(|_| RDKitError::BadFile(format!("Bad output file {path}")))
}
