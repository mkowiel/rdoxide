//! SCSR mol files: a V3000 molecule whose atoms stand for macromolecule
//! residues (amino acids, nucleotides, ...), followed by the templates that
//! define them, and their expansion into an ordinary molecule
//! (`SCSRMolFileParser.cpp` of RDKit 2026.09.1).
//!
//! Each residue of the expanded molecule is marked with a superatom
//! substance group labelled with the name of its template.
//!
//! Differences from upstream: the `molAttachOrderTemplate` property is kept
//! as `index:name` strings (see [`super::mol_file_v3000::attach_order_template`]);
//! a template without coordinates makes the expanded molecule lack
//! coordinates, where upstream would crash; and an unhandled CTAB feature
//! inside a template is an error, where upstream leaves a null pointer that
//! it then dereferences.
//!
//! Copyright (C) 2024 Tad Hurst and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;

use super::mol_file_parser::{
    MolFileParserParams, finish_mol_processing, fpe, mol_from_mol_data_stream,
};
use super::mol_file_v3000::{attach_order_template, get_v3000_line, parse_v3000_ctab};
use super::utils::{LineStream, lossy};
use crate::geometry::Point3D;
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::stereo_group::{StereoGroup, StereoGroupType};
use crate::graphmol::substance_group::{SubstanceGroup, add_substance_group};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// `SCSRTemplateNames`: which name of a template labels its residues.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ScsrTemplateNames {
    /// the name as entered in the SCSR mol
    #[default]
    AsEntered,
    /// the first name of the template definition (for amino acids, the
    /// three-letter code)
    UseFirstName,
    /// the last name of the template definition (for amino acids, the
    /// one-letter code)
    UseSecondName,
}

/// `SCSRBaseHbondOptions`: how hydrogen bonds between bases are expanded.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ScsrBaseHbondOptions {
    /// no hydrogen bonds in the output
    Ignore,
    /// all the hydrogen bonds defined by the attachment points of the
    /// templates; there can be several per base
    #[default]
    UseSapAll,
    /// only the first hydrogen bond of each base
    UseSapOne,
    /// the standard Watson-Crick (and wobble) pairing of C, G, A, T, U and
    /// inosine, found by substructure matching; attachment points are ignored
    Auto,
}

/// `MolFromSCSRParams`
#[derive(Debug, Clone, PartialEq)]
pub struct MolFromScsrParams {
    /// keep the leaving groups of atoms that are not exo-bonded
    pub include_leaving_groups: bool,
    pub scsr_template_names: ScsrTemplateNames,
    pub scsr_base_hbond_options: ScsrBaseHbondOptions,
}

impl Default for MolFromScsrParams {
    fn default() -> Self {
        MolFromScsrParams {
            include_leaving_groups: true,
            scsr_template_names: ScsrTemplateNames::default(),
            scsr_base_hbond_options: ScsrBaseHbondOptions::default(),
        }
    }
}

/// `SCSRMol`: the macromolecule and the templates of its residues.
#[derive(Debug)]
pub struct ScsrMol {
    pub mol: RWMol,
    pub templates: Vec<RWMol>,
}

// ---------------------------------------------------------------- parsing

fn skip_spaces(line: &[u8], pos: &mut usize) {
    while line.get(*pos) == Some(&b' ') {
        *pos += 1;
    }
}

/// `getToken`: the text up to `delim` (or the end of the line), which is skipped.
fn get_token(line: &[u8], pos: &mut usize, delim: u8) -> Vec<u8> {
    skip_spaces(line, pos);
    let start = *pos;
    while *pos < line.len() && line[*pos] != delim {
        *pos += 1;
    }
    let res = line[start..*pos].to_vec();
    if line.get(*pos) == Some(&delim) {
        *pos += 1;
    }
    res
}

/// `getQuotedToken`: a `"quoted"` string with backslash escapes; empty if the
/// closing quote is missing.
fn get_quoted_token(line: &[u8], pos: &mut usize) -> Vec<u8> {
    skip_spaces(line, pos);
    let mut res = Vec::new();
    if line.get(*pos) != Some(&b'"') {
        return res;
    }
    *pos += 1;
    while *pos < line.len() && line[*pos] != b'"' {
        if line[*pos] == b'\\' {
            *pos += 1;
            if *pos >= line.len() {
                break;
            }
        }
        res.push(line[*pos]);
        *pos += 1;
    }
    if line.get(*pos) != Some(&b'"') {
        return Vec::new();
    }
    *pos += 1;
    res
}

/// `parseTemplateLine`:
/// `TEMPLATE 1 AA/Cya/Cya/ NATREPLACE=AA/A COMMENT=comment FULLNAME="full name" ...`
/// The class and names become properties of the template, as do the other
/// attributes (kept as strings).
fn parse_template_line(template: &mut RWMol, line_str: &[u8], line: u32) -> Result<()> {
    if !line_str.starts_with(b"TEMPLATE ") {
        return Err(fpe(format!("Expected \"TEMPLATE\" at line  {line}")));
    }
    let mut pos = 9;
    // the template ID
    if get_token(line_str, &mut pos, b' ').is_empty() {
        return Err(fpe(format!("Expected a Template ID at line  {line}")));
    }
    // the class and template names
    let token = get_token(line_str, &mut pos, b' ');
    if token.is_empty() {
        return Err(fpe(format!(
            "Type/Name(s) string of the form \"AA/Gly/G/\" was not found at line  {line}"
        )));
    }
    let sub_tokens: Vec<&[u8]> = token.split(|&c| c == b'/').collect();
    if sub_tokens.len() < 3 {
        return Err(fpe(format!(
            "Type/Name(s) string is not of the form \"AA/Gly/G/\" at line  {line}"
        )));
    }
    template
        .props
        .set_val(cp::molAtomClass, lossy(sub_tokens[0]));
    let names: Vec<String> = sub_tokens[1..]
        .iter()
        .filter(|s| !s.is_empty())
        .map(|s| lossy(s))
        .collect();
    template.props.set_val(cp::templateNames, names);

    // now the attributes of the form ATTRNAME=VALUE or ATTRNAME="VALUE WITH SPACES"
    loop {
        let name = get_token(line_str, &mut pos, b'=');
        if name.is_empty() {
            break;
        }
        let value = if line_str.get(pos) == Some(&b'"') {
            get_quoted_token(line_str, &mut pos)
        } else {
            get_token(line_str, &mut pos, b' ')
        };
        template.props.set_val(&lossy(&name), lossy(&value));
    }
    if pos < line_str.len() {
        return Err(fpe(format!(
            "extra characters at the end of a TEMPLATE definition line at line  {line}"
        )));
    }
    Ok(())
}

/// The string property `key` of an atom or molecule.
fn string_prop(props: &impl RDProps, key: &str) -> Option<String> {
    props.get_prop(key).map(RDValue::to_rd_string)
}

/// The names of a template.
fn template_names(template: &ROMol) -> Vec<String> {
    match template.props.get(cp::templateNames) {
        Some(RDValue::VecString(names)) => names.clone(),
        _ => Vec::new(),
    }
}

/// The superatom group of a template with the given `CLASS`.
fn find_sup_sgroup<'a>(template: &'a ROMol, class: &str) -> Option<&'a SubstanceGroup> {
    template.substance_groups().iter().find(|sg| {
        string_prop(*sg, "TYPE").as_deref() == Some("SUP")
            && string_prop(*sg, "CLASS").as_deref() == Some(class)
    })
}

/// `SCSRMolFromSCSRDataStream`: the macromolecule and its templates, checked
/// against each other.
pub fn scsr_mol_from_scsr_data_stream(
    s: &mut LineStream<'_>,
    line: &mut u32,
    params: &MolFileParserParams,
) -> Result<Option<ScsrMol>> {
    let mut chirality_possible = false;
    if s.eof() {
        return Ok(None);
    }
    let local_params = MolFileParserParams {
        parsing_scsr_mol: true,
        ..params.clone()
    };
    let mol = mol_from_mol_data_stream(s, line, &local_params)?
        .ok_or_else(|| RDKitError::precondition("bad molecule"))?;
    let mut res = ScsrMol {
        mol,
        templates: Vec::new(),
    };

    // now get all of the templates
    let mut temp = get_v3000_line(s, line)?;
    if temp != b"BEGIN TEMPLATE" {
        return Err(fpe(format!("BEGIN TEMPLATE not found at line  {line}")));
    }
    temp = get_v3000_line(s, line)?;
    while temp.starts_with(b"TEMPLATE") {
        let mut template = RWMol::new();
        parse_template_line(&mut template, &temp, *line)?;
        let complete = parse_v3000_ctab(
            s,
            line,
            &mut template,
            &mut chirality_possible,
            params.strict_parsing,
            false,
            false,
        )?;
        if !complete {
            return Err(fpe(format!(
                "Problems encountered parsing Mol data, M  END or BEGIN TEMPLATE missing around line {line}"
            )));
        }
        let template_params = MolFileParserParams {
            sanitize: false,
            remove_hs: false,
            ..params.clone()
        };
        finish_mol_processing(&mut template, chirality_possible, &template_params)?;
        res.templates.push(template);
        temp = get_v3000_line(s, line)?;
    }
    if temp != b"END TEMPLATE" {
        return Err(fpe(format!("END TEMPLATE not found at line  {line}")));
    }
    *line += 1;
    if !s.get_line().starts_with(b"M  END") {
        return Err(fpe(format!("M  END not found at line  {line}")));
    }

    check_templates(&res)?;
    Ok(Some(res))
}

/// Each atom of the macromolecule that needs a template must match one, with
/// all the attachment points it uses.
fn check_templates(scsr: &ScsrMol) -> Result<()> {
    for atom in scsr.mol.atoms() {
        let (Some(label), Some(class)) = (
            string_prop(atom, cp::dummyLabel),
            string_prop(atom, cp::molAtomClass),
        ) else {
            continue;
        };
        if label.is_empty() || class.is_empty() {
            continue;
        }
        let attach_ords = attach_order_template(atom);
        let mut found = false;
        'templates: for template in &scsr.templates {
            if string_prop(template, cp::molAtomClass).as_deref() != Some(class.as_str()) {
                continue;
            }
            for name in template_names(template) {
                if found {
                    break 'templates;
                }
                if name != label {
                    continue;
                }
                // the SUP group for the main connections (not leaving groups)
                let Some(main_sup) = find_sup_sgroup(template, &class) else {
                    break; // on to the next template
                };
                // tentative, until some attachment point of the atom is not
                // found in the template
                found = true;
                for (_, attach_label) in &attach_ords {
                    if !main_sup
                        .attach_points()
                        .iter()
                        .any(|a| &a.id == attach_label)
                    {
                        found = false;
                        break;
                    }
                }
            }
        }
        if !found {
            return Err(fpe(format!(
                "No template found for atom {} class {class} label {label}",
                atom.idx()
            )));
        }
    }
    Ok(())
}

// -------------------------------------------------------------- expansion

/// `UINT_MAX`: the template atom of an atom that is not from a template
const NOT_TEMPLATE: u32 = u32::MAX;

#[derive(Debug, Clone, Copy)]
struct HydrogenBondConnection {
    template_atom_idx: u32,
    is_donor: bool,
}

const fn hb(template_atom_idx: u32, is_donor: bool) -> HydrogenBondConnection {
    HydrogenBondConnection {
        template_atom_idx,
        is_donor,
    }
}

/// The base pairs recognised by [`ScsrBaseHbondOptions::Auto`]: SMARTS, name
/// and the hydrogen bonding atoms (indices into the SMARTS).
const HBOND_QUERIES: [(&str, &str, &[HydrogenBondConnection]); 5] = [
    (
        "[NH]1C=NC2=C1N=C([NH2])[NH]C2=O",
        "Guanine",
        &[hb(7, true), hb(8, true), hb(10, false)],
    ),
    (
        "[NH]1C=NC2=C1N=[CH]N=C2[NH2]",
        "Adenine",
        &[hb(6, true), hb(7, false), hb(9, true)],
    ),
    (
        "[NH]1C=CC([NH2])=NC1=O",
        "Cytosine",
        &[hb(7, false), hb(5, false), hb(4, true)],
    ),
    (
        "[NH]1C=CC(=O)[NH]C1=O",
        "Thyamine-Uracil",
        &[hb(7, false), hb(5, true), hb(4, false)],
    ),
    (
        "[NH]1C=NC2=C1N=C[NH]C2=O",
        "Inosine",
        &[hb(7, true), hb(9, false)],
    ),
];

struct Converter<'a> {
    templates: Vec<RWMol>,
    res: RWMol,
    params: &'a MolFromScsrParams,
    h_bond_sites_for_template: BTreeMap<u32, Vec<HydrogenBondConnection>>,
    atom_idx_to_template_idx: BTreeMap<u32, u32>,
    /// (main atom, template atom) -> new atom
    origin_atom_map: BTreeMap<(u32, u32), u32>,
    /// (main atom, attachment point name) -> template atoms
    attach_map: BTreeMap<(u32, String), Vec<u32>>,
}

impl Converter<'_> {
    fn new_atom_for_bond(&self, mol: &ROMol, atom: u32, other_atom_idx: u32) -> Result<u32> {
        let a = mol.atom(atom);
        if string_prop(a, cp::molAtomClass).is_none() {
            return self
                .origin_atom_map
                .get(&(atom, NOT_TEMPLATE))
                .copied()
                .ok_or_else(|| fpe("Error getting new atom for bond"));
        }
        // a template atom; this is not called for hydrogen bonds, so there is
        // only one atom to attach
        for (idx, label) in attach_order_template(a) {
            if idx == other_atom_idx {
                let first = self
                    .attach_map
                    .get(&(atom, label))
                    .and_then(|atoms| atoms.first())
                    .ok_or_else(|| fpe("Attachment ord not found"))?;
                return Ok(self.origin_atom_map[&(atom, *first)]);
            }
        }
        Err(fpe("Attachment ord not found"))
    }

    fn new_atoms_for_hydrogen_bond(
        &mut self,
        mol: &ROMol,
        atom: u32,
        other_atom_idx: u32,
    ) -> Result<Vec<HydrogenBondConnection>> {
        let a = mol.atom(atom);
        if string_prop(a, cp::molAtomClass).is_none() {
            return Ok(Vec::new()); // hatoms only allowed in templates
        }
        // There could be more than one atom in the template that matches the
        // atom: the hydrogen bonds to the template can result in multiple
        // hydrogen bonds in the expanded molecule.
        let template_idx = self.atom_idx_to_template_idx[&atom] as usize;
        let mut connections = Vec::new();
        for (idx, label) in attach_order_template(a) {
            if idx != other_atom_idx {
                continue;
            }
            let Some(attach_atoms) = self.attach_map.get(&(atom, label)) else {
                return Ok(Vec::new()); // error, attachment ord not found
            };
            for &template_atom_idx in attach_atoms {
                let template = &mut self.templates[template_idx];
                template.update_atom_property_cache(template_atom_idx, true)?;
                let is_donor = template.total_num_hs(template_atom_idx, false)? > 0;
                connections.push(hb(template_atom_idx, is_donor));
                if self.params.scsr_base_hbond_options == ScsrBaseHbondOptions::UseSapOne {
                    return Ok(connections);
                }
            }
            return Ok(connections);
        }
        // error, attachment ord not found
        Ok(Vec::new())
    }

    /// `copySgroupIntoResult`: adds the atoms of a template group to the
    /// result, as a new superatom group labelled `sgroup_name`.
    fn copy_sgroup_into_result(
        &mut self,
        atom_idx: u32,
        sgroup: &SubstanceGroup,
        sgroup_name: &str,
        new_sgroups: &mut Vec<SubstanceGroup>,
        new_conf: &mut Option<Conformer>,
        coord_offset: Point3D,
    ) -> Result<()> {
        // These new superatom groups are not put into the output molecule
        // yet, because the bonds have not been added to the molecule nor to
        // the group.
        let mut new_sgroup = SubstanceGroup::new("SUP");
        new_sgroup.set_prop("LABEL", sgroup_name.to_string(), false)?;

        let template = &self.templates[self.atom_idx_to_template_idx[&atom_idx] as usize];
        if let Some(conf) = new_conf {
            conf.resize(conf.num_atoms() as usize + template.num_atoms() as usize);
        }
        for &template_atom_idx in sgroup.atoms() {
            let new_idx = self
                .res
                .add_atom(template.atom(template_atom_idx).clone(), true);
            new_sgroup.add_atom_with_idx(&self.res, new_idx)?;
            self.origin_atom_map
                .insert((atom_idx, template_atom_idx), new_idx);
            if let Some(conf) = new_conf {
                let pos = template.conformer(-1)?.atom_pos(template_atom_idx)?;
                conf.set_atom_pos(new_idx, coord_offset + pos)?;
            }
        }
        new_sgroups.push(new_sgroup);
        Ok(())
    }

    fn add_bond(&mut self, bond_type: BondType, begin: u32, end: u32) -> Result<()> {
        let mut bond = Bond::new(bond_type);
        bond.set_begin_atom_idx_unchecked(begin);
        bond.set_end_atom_idx_unchecked(end);
        self.res.add_bond_object(bond)?;
        Ok(())
    }

    fn add_bonds(
        &mut self,
        bond_type: BondType,
        main_begin: u32,
        main_end: u32,
        begin_atoms: &[HydrogenBondConnection],
        end_atoms: &[HydrogenBondConnection],
    ) -> Result<()> {
        for (b, e) in begin_atoms.iter().zip(end_atoms) {
            let begin = self.origin_atom_map[&(main_begin, b.template_atom_idx)];
            let end = self.origin_atom_map[&(main_end, e.template_atom_idx)];
            self.add_bond(bond_type, begin, end)?;
        }
        Ok(())
    }

    fn process_bond_in_main_mol(&mut self, mol: &ROMol, bond_idx: u32) -> Result<()> {
        let bond = mol.bond(bond_idx);
        let bond_type = bond.bond_type();
        let (main_begin, main_end) = (bond.begin_atom_idx(), bond.end_atom_idx());
        if bond_type != BondType::Hydrogen {
            let new_begin = self.new_atom_for_bond(mol, main_begin, main_end)?;
            let new_end = self.new_atom_for_bond(mol, main_end, main_begin)?;
            return self.add_bond(bond_type, new_begin, new_end);
        }

        // A hydrogen bond. How these are processed is controlled by
        // `ScsrBaseHbondOptions`; see its documentation.
        //
        // With all the attachment points of the templates (`UseSapAll`) they
        // must be defined in the right order, from the atom nearest the Al
        // connection on. If there are 3 sites on each side and they are
        // complementary (the donors match acceptors and vice versa), the
        // bonds are added.
        //
        // If not, the sites are checked for a wobble pair. The four generally
        // accepted wobble pairs are I-C, I-U, I-A and G-U. "I" (inosine) has
        // only two sites, which must be DA, so the two atoms form the two
        // hydrogen bonds; the other side (C, U or A) is AAD (C), ADA (U) or
        // DAD (A): for C and A the second and third atoms are used (AD), for
        // U the first two (AD). G and U both have three sites, which are not
        // complementary: the second and third sites of G (DA) and the first
        // two of U (AD) are used.
        //
        // In any other case just one bond is added between the first atoms of
        // both sides, even if they are not complementary: it shows that the
        // sides are hydrogen bonded somehow and keeps the pairing straight.
        let mut begin_connections;
        let mut end_connections;
        match self.params.scsr_base_hbond_options {
            ScsrBaseHbondOptions::Ignore => return Ok(()),
            ScsrBaseHbondOptions::UseSapAll | ScsrBaseHbondOptions::UseSapOne => {
                begin_connections = self.new_atoms_for_hydrogen_bond(mol, main_begin, main_end)?;
                if begin_connections.is_empty() {
                    return Err(fpe("No Attach Point for bond"));
                }
                end_connections = self.new_atoms_for_hydrogen_bond(mol, main_end, main_begin)?;
                if end_connections.is_empty() {
                    return Err(fpe("No Attach Point for bond"));
                }
            }
            ScsrBaseHbondOptions::Auto => {
                // the sites already found for each base template
                let sites = |atom: u32, this: &Self| {
                    this.atom_idx_to_template_idx
                        .get(&atom)
                        .and_then(|t| this.h_bond_sites_for_template.get(t))
                        .cloned()
                        .unwrap_or_default()
                };
                begin_connections = sites(main_begin, self);
                end_connections = sites(main_end, self);
            }
        }
        let (mut main_begin, mut main_end) = (main_begin, main_end);
        if begin_connections.is_empty() || end_connections.is_empty() {
            return Ok(()); // no hydrogen bond sites found, so skip the bond
        }
        if begin_connections.len() == 3 && end_connections.len() == 3 {
            // are the two donor flag lists complementary?
            if (0..3).all(|i| begin_connections[i].is_donor != end_connections[i].is_donor) {
                return self.add_bonds(
                    bond_type,
                    main_begin,
                    main_end,
                    &begin_connections,
                    &end_connections,
                );
            }
            // a G-U type pair in a wobble configuration: first see if one has
            // the configuration DDA (like G)
            let is_dda =
                |c: &[HydrogenBondConnection]| c[0].is_donor && c[1].is_donor && !c[2].is_donor;
            let mut found_dda = false;
            if is_dda(&begin_connections) {
                found_dda = true;
            } else if is_dda(&end_connections) {
                found_dda = true;
                std::mem::swap(&mut begin_connections, &mut end_connections);
                std::mem::swap(&mut main_begin, &mut main_end);
            }
            if found_dda {
                let wobble_begin = [begin_connections[1], begin_connections[2]];
                // does the other side have the configuration ADA (like U)?
                if !end_connections[0].is_donor
                    && end_connections[1].is_donor
                    && !end_connections[2].is_donor
                {
                    // use the first two atoms
                    let wobble_end = [end_connections[0], end_connections[1]];
                    return self.add_bonds(
                        bond_type,
                        main_begin,
                        main_end,
                        &wobble_begin,
                        &wobble_end,
                    );
                }
            }
        } else if begin_connections.len() * end_connections.len() == 6 {
            // one has 2 sites and the other 3; make the first the one with 2
            if end_connections.len() == 2 {
                std::mem::swap(&mut begin_connections, &mut end_connections);
                std::mem::swap(&mut main_begin, &mut main_end);
            }
            if begin_connections[0].is_donor && !begin_connections[1].is_donor {
                // like I (DA)
                if !end_connections[1].is_donor && end_connections[2].is_donor {
                    // the 2nd and 3rd are AD, like A (DAD) or C (AAD)
                    let wobble_end = [end_connections[1], end_connections[2]];
                    return self.add_bonds(
                        bond_type,
                        main_begin,
                        main_end,
                        &begin_connections,
                        &wobble_end,
                    );
                } else if !end_connections[0].is_donor && end_connections[1].is_donor {
                    // like U (ADA)
                    let wobble_end = [end_connections[0], end_connections[1]];
                    return self.add_bonds(
                        bond_type,
                        main_begin,
                        main_end,
                        &begin_connections,
                        &wobble_end,
                    );
                }
            }
        }
        // no wobble pair: just add the bond between the first atoms on both sides
        let begin = self.origin_atom_map[&(main_begin, begin_connections[0].template_atom_idx)];
        let end = self.origin_atom_map[&(main_end, end_connections[0].template_atom_idx)];
        self.add_bond(bond_type, begin, end)
    }

    /// The hydrogen bonding sites of the base templates, found by
    /// substructure matching (`Auto`).
    fn find_base_hydrogen_bond_sites(&mut self) -> Result<()> {
        let queries = HBOND_QUERIES
            .iter()
            .map(|(smarts, _, sites)| {
                let query = mol_from_smarts(smarts, &SmartsParserParams::default())?
                    .ok_or_else(|| RDKitError::precondition("bad hydrogen bond query"))?;
                Ok((query, *sites))
            })
            .collect::<Result<Vec<_>>>()?;
        for (template_idx, template) in self.templates.iter_mut().enumerate() {
            template.update_property_cache(false)?;
            if string_prop(template, cp::molAtomClass).as_deref() != Some("BASE") {
                continue;
            }
            let found = self
                .h_bond_sites_for_template
                .entry(template_idx as u32)
                .or_default();
            for (query, sites) in &queries {
                let matches =
                    substruct_match(template, query, &SubstructMatchParameters::default())?;
                let Some(first) = matches.first() else {
                    continue;
                };
                for site in *sites {
                    found.push(hb(
                        first[site.template_atom_idx as usize].1 as u32,
                        site.is_donor,
                    ));
                }
                break; // only take the first match
            }
        }
        Ok(())
    }
}

/// The index of the first bond of `mol` between two atoms.
fn bond_between(mol: &ROMol, a: u32, b: u32) -> Result<Option<u32>> {
    Ok(mol.get_bond_between_atoms(a, b)?.map(Bond::idx))
}

/// `MolFromSCSRMolConverter::convert`: expands the residues of the
/// macromolecule into the atoms and bonds of their templates.
pub fn mol_from_scsr_mol(
    scsr: ScsrMol,
    mol_file_params: &MolFileParserParams,
    params: &MolFromScsrParams,
) -> Result<RWMol> {
    let ScsrMol { mol, templates } = scsr;
    let mut cv = Converter {
        templates,
        res: RWMol::new(),
        params,
        h_bond_sites_for_template: BTreeMap::new(),
        atom_idx_to_template_idx: BTreeMap::new(),
        origin_atom_map: BTreeMap::new(),
        attach_map: BTreeMap::new(),
    };

    // First get some information from the templates to be used when creating
    // the coordinates of the new atoms. This is a dirty approach that simply
    // expands the original macro atom coordinates to be big enough to hold
    // any expanded macro atom. No attempt is made to make this look nice, or
    // to avoid overlaps.
    let mut template_centroids: Vec<Point3D> = Vec::new();
    let mut max_size = 0.0f64;
    let mut main_conf: Option<&Conformer> = None;
    let mut new_conf: Option<Conformer> = None;
    if mol.num_conformers() != 0 {
        main_conf = Some(mol.conformer(0)?);
        let mut conf = Conformer::new(mol.num_atoms() as usize);
        conf.set_3d(mol.conformer(0)?.is_3d());
        new_conf = Some(conf);
        for template in &cv.templates {
            if template.num_conformers() == 0 {
                main_conf = None;
                new_conf = None;
                break;
            }
            let template_conf = template.conformer(0)?;
            let mut sum = Point3D::default();
            let mut max_coord = template_conf.atom_pos(0)?;
            let mut min_coord = max_coord;
            for atom_idx in 0..template.num_atoms() {
                let c = template_conf.atom_pos(atom_idx)?;
                sum += c;
                max_coord = Point3D::new(
                    max_coord.x.max(c.x),
                    max_coord.y.max(c.y),
                    max_coord.z.max(c.z),
                );
                min_coord = Point3D::new(
                    min_coord.x.min(c.x),
                    min_coord.y.min(c.y),
                    min_coord.z.min(c.z),
                );
            }
            template_centroids.push(sum / f64::from(template.num_atoms()));
            for extent in [
                max_coord.x - min_coord.x,
                max_coord.y - min_coord.y,
                max_coord.z - min_coord.z,
            ] {
                if extent > max_size {
                    max_size = extent;
                }
            }
        }
    }

    // if we are perceiving the H-bonding locations for base templates, do
    // that here
    if params.scsr_base_hbond_options == ScsrBaseHbondOptions::Auto {
        cv.find_base_hydrogen_bond_sites()?;
    }

    // for each atom in the main molecule, expand it to its full form
    let mut new_stereo_groups: Vec<StereoGroup> = Vec::new();
    let mut absolute_atoms: Vec<u32> = Vec::new();
    let mut absolute_bonds: Vec<u32> = Vec::new();
    let mut new_sgroups: Vec<SubstanceGroup> = Vec::new();

    for atom_idx in 0..mol.num_atoms() {
        let atom = mol.atom(atom_idx);
        let dummy_label = string_prop(atom, cp::dummyLabel).unwrap_or_default();
        let atom_class = string_prop(atom, cp::molAtomClass).unwrap_or_default();

        if dummy_label.is_empty() || atom_class.is_empty() {
            // not a template atom: just copy it
            let new_idx = cv.res.add_atom(atom.clone(), true);
            cv.origin_atom_map.insert((atom_idx, NOT_TEMPLATE), new_idx);
            if let (Some(conf), Some(new_conf)) = (main_conf, new_conf.as_mut()) {
                new_conf.set_atom_pos(new_idx, conf.atom_pos(atom_idx)? * max_size)?;
            }
            continue;
        }

        // a macro atom: expand it
        let seq_id = atom
            .get_prop(cp::molAtomSeqId)
            .and_then(RDValue::as_int)
            .unwrap_or(0) as u32;
        let seq_name = string_prop(atom, cp::molAtomSeqName).unwrap_or_default();

        // find the template that matches the class and label
        let mut found: Option<(usize, String)> = None;
        for (template_idx, template) in cv.templates.iter().enumerate() {
            let names = template_names(template);
            if string_prop(template, cp::molAtomClass).as_deref() == Some(atom_class.as_str())
                && names.contains(&dummy_label)
            {
                let name = match params.scsr_template_names {
                    ScsrTemplateNames::UseFirstName => names[0].clone(),
                    ScsrTemplateNames::UseSecondName => names[names.len() - 1].clone(),
                    ScsrTemplateNames::AsEntered => dummy_label.clone(),
                };
                found = Some((template_idx, name));
                break;
            }
        }
        let Some((template_idx, template_name)) = found else {
            return Err(fpe(format!(
                "No template found for atom {atom_idx} class {atom_class} label {dummy_label}"
            )));
        };
        cv.atom_idx_to_template_idx
            .insert(atom_idx, template_idx as u32);
        let attach_ords = attach_order_template(atom);

        // first find the group that is the base for the template of this atom
        let Some(sgroup) = find_sup_sgroup(&cv.templates[template_idx], &atom_class).cloned()
        else {
            return Err(fpe("No SUP sgroup found for atom"));
        };

        // add the atoms of the main template to the new molecule
        let mut sgroup_name = template_name;
        if seq_id != 0 {
            sgroup_name += &format!("_{seq_id}");
        }
        if !seq_name.is_empty() {
            sgroup_name += &format!("_{seq_name}");
        }
        let coord_offset = match main_conf {
            Some(conf) => conf.atom_pos(atom_idx)? * max_size - template_centroids[template_idx],
            None => Point3D::default(),
        };
        cv.copy_sgroup_into_result(
            atom_idx,
            &sgroup,
            &sgroup_name,
            &mut new_sgroups,
            &mut new_conf,
            coord_offset,
        )?;

        // find the attachment points used by this atom and record them
        for (_, label) in &attach_ords {
            let attach_atoms: Vec<u32> = sgroup
                .attach_points()
                .iter()
                .filter(|a| &a.id == label)
                .map(|a| a.a_idx)
                .collect();
            if attach_atoms.is_empty() {
                return Err(fpe("Attachment point not found"));
            }
            cv.attach_map
                .insert((atom_idx, label.clone()), attach_atoms);
        }

        // If we are including the atoms of leaving groups, go through the
        // attachment points of the main group. If the attachment point is not
        // used by the atom, find the group for that point and add its atoms
        // to the molecule.
        if params.include_leaving_groups {
            for attach_point in sgroup.attach_points() {
                if cv
                    .attach_map
                    .contains_key(&(atom_idx, attach_point.id.clone()))
                {
                    continue;
                }
                let lg_sgroup = cv.templates[template_idx]
                    .substance_groups()
                    .iter()
                    .find(|lg| {
                        string_prop(*lg, "TYPE").as_deref() == Some("SUP")
                            && string_prop(*lg, "CLASS").as_deref() == Some("LGRP")
                            && lg.atoms().iter().any(|&a| a as i32 == attach_point.lv_idx)
                    })
                    .cloned();
                let Some(lg_sgroup) = lg_sgroup else {
                    return Err(fpe(format!(
                        "Leaving group SGroup {} not found",
                        attach_point.id
                    )));
                };
                let mut lg_name = dummy_label.clone();
                if seq_id != 0 {
                    lg_name += &format!("_{seq_id}");
                }
                if !seq_name.is_empty() {
                    lg_name += &format!("_{seq_name}");
                }
                lg_name += &format!("_{}", attach_point.id);
                cv.copy_sgroup_into_result(
                    atom_idx,
                    &lg_sgroup,
                    &lg_name,
                    &mut new_sgroups,
                    &mut new_conf,
                    coord_offset,
                )?;
            }
        }

        // Copy the bonds of the template into the new molecule if they are
        // between atoms in the new molecule. Bonds to atoms in leaving groups
        // that "left" are not copied.
        let template = &cv.templates[template_idx];
        for bond in template.bonds() {
            let (Some(&new_begin), Some(&new_end)) = (
                cv.origin_atom_map.get(&(atom_idx, bond.begin_atom_idx())),
                cv.origin_atom_map.get(&(atom_idx, bond.end_atom_idx())),
            ) else {
                continue; // the bond is not in the new molecule
            };
            let mut new_bond = Bond::new(bond.bond_type());
            new_bond.set_begin_atom_idx_unchecked(new_begin);
            new_bond.set_end_atom_idx_unchecked(new_end);
            new_bond.update_props(bond.dict(), false);
            cv.res.add_bond_object(new_bond)?;
        }

        // The stereo groups of the template. Absolute groups are added to the
        // list of absolute atoms and bonds, so that one absolute group can
        // be added later.
        for sg in template.stereo_groups() {
            let group_atoms: Vec<u32> = sg
                .atoms()
                .iter()
                .filter_map(|&a| cv.origin_atom_map.get(&(atom_idx, a)).copied())
                .collect();
            let mut group_bonds = Vec::new();
            for &b in sg.bonds() {
                let tb = template.bond(b);
                if let (Some(&nb), Some(&ne)) = (
                    cv.origin_atom_map.get(&(atom_idx, tb.begin_atom_idx())),
                    cv.origin_atom_map.get(&(atom_idx, tb.end_atom_idx())),
                ) && let Some(new_bond) = bond_between(&cv.res, nb, ne)?
                {
                    group_bonds.push(new_bond);
                }
            }
            if sg.group_type() == StereoGroupType::Absolute {
                absolute_atoms.extend(group_atoms);
                absolute_bonds.extend(group_bonds);
            } else {
                new_stereo_groups.push(StereoGroup::new(
                    sg.group_type(),
                    group_atoms,
                    group_bonds,
                    0,
                )?);
            }
        }
    }

    if cv.res.num_atoms() == 0 {
        return Ok(cv.res);
    }

    if let Some(mut conf) = new_conf {
        conf.resize(cv.res.num_atoms() as usize);
        cv.res.add_conformer(conf, false)?;
    }

    // now deal with the bonds of the original molecule
    for bond_idx in 0..mol.num_bonds() {
        cv.process_bond_in_main_mol(&mol, bond_idx)?;
    }

    // copy any attributes of the main molecule
    for prop in mol.get_prop_list(false, false) {
        if let Some(value) = string_prop(&mol, &prop) {
            cv.res.props.set_val(&prop, value);
        }
    }

    // copy the groups of the main molecule for atoms that are not in a template
    for sg in mol.substance_groups() {
        if !sg.is_valid() {
            continue;
        }
        let mut new_atoms = Vec::new();
        for &atom in sg.atoms() {
            match cv.origin_atom_map.get(&(atom, NOT_TEMPLATE)) {
                Some(&new_idx) => new_atoms.push(new_idx),
                None => {
                    // some atoms were in templates and others were not:
                    // this group cannot be added
                    new_atoms.clear();
                    break;
                }
            }
        }
        if !new_atoms.is_empty() {
            let mut new_sg = SubstanceGroup::new("SUP");
            new_sg.set_atoms(&cv.res, new_atoms)?;
            new_sgroups.push(new_sg);
        }
    }

    // Now that we have all the groups from the templates and from the atoms
    // that are not from templates, and all the bonds, find the crossing bonds
    // of each group and add them.
    for bond_idx in 0..cv.res.num_bonds() {
        let (begin, end) = {
            let b = cv.res.bond(bond_idx);
            (b.begin_atom_idx(), b.end_atom_idx())
        };
        for sg in &mut new_sgroups {
            // exactly one atom of the bond is in the group
            if sg.atoms().contains(&begin) != sg.atoms().contains(&end) {
                sg.add_bond_with_idx(&cv.res, bond_idx)?;
            }
        }
    }
    for sg in new_sgroups {
        add_substance_group(&mut cv.res, sg);
    }

    // The stereo groups of the main molecule, for atoms that are not
    // template references. As upstream, the groups are only collected inside
    // the loop over their bonds.
    for sg in mol.stereo_groups() {
        let group_atoms: Vec<u32> = sg
            .atoms()
            .iter()
            .filter_map(|&a| cv.origin_atom_map.get(&(a, NOT_TEMPLATE)).copied())
            .collect();
        let mut group_bonds: Vec<u32> = Vec::new();
        for &b in sg.bonds() {
            let mb = mol.bond(b);
            if let (Some(&nb), Some(&ne)) = (
                cv.origin_atom_map.get(&(mb.begin_atom_idx(), NOT_TEMPLATE)),
                cv.origin_atom_map.get(&(mb.end_atom_idx(), NOT_TEMPLATE)),
            ) && let Some(new_bond) = bond_between(&cv.res, nb, ne)?
            {
                group_bonds.push(new_bond);
            }
            if sg.group_type() == StereoGroupType::Absolute {
                absolute_atoms.extend(group_atoms.iter().copied());
                absolute_bonds.extend(group_bonds.iter().copied());
            } else {
                new_stereo_groups.push(StereoGroup::new(
                    sg.group_type(),
                    group_atoms.clone(),
                    group_bonds.clone(),
                    0,
                )?);
            }
        }
    }

    // an absolute group with any absolute atoms or bonds of the main
    // molecule or the template instances
    if !absolute_atoms.is_empty() || !absolute_bonds.is_empty() {
        new_stereo_groups.push(StereoGroup::new(
            StereoGroupType::Absolute,
            absolute_atoms,
            absolute_bonds,
            0,
        )?);
    }
    if !new_stereo_groups.is_empty() {
        cv.res.set_stereo_groups(new_stereo_groups)?;
    }

    finish_mol_processing(&mut cv.res, false, mol_file_params)?;
    Ok(cv.res)
}

/// `MolFromSCSRDataStream`
pub fn mol_from_scsr_data_stream(
    s: &mut LineStream<'_>,
    line: &mut u32,
    mol_file_params: &MolFileParserParams,
    params: &MolFromScsrParams,
) -> Result<Option<RWMol>> {
    match scsr_mol_from_scsr_data_stream(s, line, mol_file_params)? {
        Some(scsr) => mol_from_scsr_mol(scsr, mol_file_params, params).map(Some),
        None => Ok(None),
    }
}

/// `MolFromSCSRBlock`
pub fn mol_from_scsr_block(
    mol_block: &str,
    mol_file_params: &MolFileParserParams,
    params: &MolFromScsrParams,
) -> Result<Option<RWMol>> {
    let mut s = LineStream::new(mol_block.as_bytes());
    let mut line = 0;
    mol_from_scsr_data_stream(&mut s, &mut line, mol_file_params, params)
}

/// `MolFromSCSRFile`
pub fn mol_from_scsr_file(
    file_name: &str,
    mol_file_params: &MolFileParserParams,
    params: &MolFromScsrParams,
) -> Result<Option<RWMol>> {
    let data = std::fs::read(file_name)
        .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))?;
    let mut s = LineStream::new(&data);
    let mut line = 0;
    mol_from_scsr_data_stream(&mut s, &mut line, mol_file_params, params)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokens_of_a_template_line() {
        let mut template = RWMol::new();
        parse_template_line(
            &mut template,
            b"TEMPLATE 1 AA/Cya/Cya/ NATREPLACE=AA/A COMMENT=\"a \\\"quoted\\\" comment\" UNIQUEID=u1",
            1,
        )
        .unwrap();
        assert_eq!(string_prop(&template, cp::molAtomClass).unwrap(), "AA");
        assert_eq!(template_names(&template), ["Cya", "Cya"]);
        assert_eq!(string_prop(&template, "NATREPLACE").unwrap(), "AA/A");
        assert_eq!(
            string_prop(&template, "COMMENT").unwrap(),
            "a \"quoted\" comment"
        );
        assert_eq!(string_prop(&template, "UNIQUEID").unwrap(), "u1");
    }

    #[test]
    fn bad_template_lines() {
        let mut t = RWMol::new();
        assert!(parse_template_line(&mut t, b"TEMPLAT 1 AA/G/G/", 1).is_err());
        assert!(parse_template_line(&mut t, b"TEMPLATE ", 1).is_err());
        assert!(parse_template_line(&mut t, b"TEMPLATE 1", 1).is_err());
        assert!(parse_template_line(&mut t, b"TEMPLATE 1 AA/G", 1).is_err());
    }
}
