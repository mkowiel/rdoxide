//! SubstanceGroup output of MDL mol blocks (V2000 `M  S..` lines and the V3000
//! `BEGIN SGROUP` block), and the groups that carry properties mol files have
//! no other place for (SMARTS queries, zero-order bonds).
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/MolSGroupWriting.{h,cpp}` and
//! the `moveAdditionalPropertiesToSGroups` part of `MolFileWriter.cpp`
//! (2026.09.1).
//!
//! Copyright (C) 2002-2018 Greg Landrum and T5 Informatics GmbH.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::mol_file_writer::is_list_query;
use crate::graphmol::bond::BondType;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts_write::get_atom_smarts;
use crate::graphmol::smiles_parse::write::SmilesWriteParams;
use crate::graphmol::substance_group::{SGroupBondType, SubstanceGroup, add_substance_group};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

// ------------------------------------------------------------ formatters

fn int_field(value: i64) -> String {
    format!(" {value:3}")
}

fn num_entries_field(value: usize) -> String {
    format!(" {value:2}")
}

fn double_field(value: f64) -> String {
    format!("{value:10.4}")
}

fn v3000_double(value: f64) -> String {
    double_field(value).trim().to_string()
}

/// `FormatV2000StringField`: cut to `size`, padded when `pad` is set.
fn string_field(value: &str, size: usize, pad: bool, separator: bool) -> String {
    let mut res = String::new();
    if separator {
        res.push(' ');
    }
    let bytes = value.as_bytes();
    if bytes.len() >= size {
        res += &String::from_utf8_lossy(&bytes[..size]);
    } else if pad {
        res += &format!("{value:<size$}");
    } else {
        res += value;
    }
    res
}

fn string_prop(sg: &SubstanceGroup, key: &str) -> Option<String> {
    sg.get_prop(key).map(RDValue::to_rd_string)
}

fn uint_prop(sg: &SubstanceGroup, key: &str) -> Option<u32> {
    sg.get_prop(key).and_then(RDValue::as_uint)
}

fn type_of(sg: &SubstanceGroup) -> String {
    string_prop(sg, "TYPE").unwrap_or_default()
}

// ----------------------------------------------------------------- V2000

/// Lines of the form `M  XXX nn` followed by up to `per_line` entries, one
/// entry for every group whose `entry` is `Some`.
fn entry_lines<I>(code: &str, per_line: usize, entries: I) -> String
where
    I: IntoIterator<Item = String>,
{
    let entries: Vec<String> = entries.into_iter().collect();
    let mut res = String::new();
    for chunk in entries.chunks(per_line) {
        res += &format!(
            "M  {code}{}{}\n",
            num_entries_field(chunk.len()),
            chunk.concat()
        );
    }
    res
}

fn group_entries<'a, F>(
    groups: &'a [SubstanceGroup],
    mut entry: F,
) -> impl Iterator<Item = String> + 'a
where
    F: FnMut(i64, &SubstanceGroup) -> Option<String> + 'a,
{
    groups
        .iter()
        .enumerate()
        .filter_map(move |(i, sg)| entry(i as i64 + 1, sg))
}

/// `BuildV2000IdxVectorDataLines`: `M  SAL` and friends, 15 indices per line.
fn idx_vector_lines(group_id: usize, code: &str, indices: &[u32]) -> String {
    let mut res = String::new();
    for chunk in indices.chunks(15) {
        let entries: String = chunk.iter().map(|&i| int_field(i64::from(i) + 1)).collect();
        res += &format!(
            "M  {code}{}{}{entries}\n",
            int_field(group_id as i64),
            num_entries_field(chunk.len())
        );
    }
    res
}

fn smt_line(idx: usize, sg: &SubstanceGroup) -> String {
    let mult = (type_of(sg) == "MUL")
        .then(|| string_prop(sg, "MULT"))
        .flatten();
    match mult.or_else(|| string_prop(sg, "LABEL")) {
        Some(v) => format!(
            "M  SMT{}{}\n",
            int_field(idx as i64),
            string_field(&v, 69, false, true)
        ),
        None => String::new(),
    }
}

fn sdi_lines(idx: usize, sg: &SubstanceGroup) -> String {
    sg.brackets()
        .iter()
        .map(|b| {
            let coords: String = b[..2]
                .iter()
                .map(|p| double_field(p.x) + &double_field(p.y))
                .collect();
            format!(
                "M  SDI{}{}{coords}\n",
                int_field(idx as i64),
                num_entries_field(4)
            )
        })
        .collect()
}

fn sbv_lines(idx: usize, sg: &SubstanceGroup) -> String {
    let is_sup = type_of(sg) == "SUP";
    sg.cstates()
        .iter()
        .map(|c| {
            let mut line = format!(
                "M  SBV{}{}",
                int_field(idx as i64),
                int_field(i64::from(c.bond_idx) + 1)
            );
            if is_sup {
                line += &(double_field(c.vector.x) + &double_field(c.vector.y));
            }
            line + "\n"
        })
        .collect()
}

fn sdt_line(idx: usize, sg: &SubstanceGroup) -> String {
    let Some(name) = string_prop(sg, "FIELDNAME") else {
        return String::new();
    };
    let mut line = format!(
        "M  SDT{}{}",
        int_field(idx as i64),
        string_field(&name, 30, true, true)
    );
    match string_prop(sg, "FIELDTYPE") {
        Some(v) => line += &string_field(&v, 2, true, false),
        None => line += " T",
    }
    for (key, width) in [("FIELDINFO", 20), ("QUERYTYPE", 2), ("QUERYOP", 15)] {
        if let Some(v) = string_prop(sg, key) {
            line += &string_field(&v, width, true, false);
        }
    }
    line + "\n"
}

fn single_string_line(code: &str, key: &str, idx: usize, sg: &SubstanceGroup) -> String {
    match string_prop(sg, key) {
        Some(v) => format!(
            "M  {code}{}{}\n",
            int_field(idx as i64),
            string_field(&v, 69, false, true)
        ),
        None => String::new(),
    }
}

/// `BuildV2000SCDSEDLines`: data fields in chunks of 69 characters.
fn scd_sed_lines(idx: usize, sg: &SubstanceGroup) -> Result<String> {
    let Some(RDValue::VecString(fields)) = sg.get_prop("DATAFIELDS") else {
        return Ok(String::new());
    };
    let mut res = String::new();
    for data in fields {
        let bytes = data.as_bytes();
        if bytes.len() > 200 {
            return Err(RDKitError::Runtime(format!(
                "Data field '{data}' in SGroup {idx} is longer than limit of 200 characters."
            )));
        }
        // upstream takes `substr(start, end)`, i.e. up to `end` characters from
        // `start`; the field formatter then cuts the chunk to 69
        let chunk = |start: usize, end: usize| {
            let from = start.min(bytes.len());
            String::from_utf8_lossy(&bytes[from..bytes.len().min(from + end)]).into_owned()
        };
        let (mut start, mut end) = (0, 69);
        while bytes.len() > end {
            res += &format!(
                "M  SCD{}{}\n",
                int_field(idx as i64),
                string_field(&chunk(start, end), 69, true, true)
            );
            start += 69;
            end += 69;
        }
        res += &format!(
            "M  SED{}{}\n",
            int_field(idx as i64),
            string_field(&chunk(start, end), 69, false, true)
        );
    }
    Ok(res)
}

fn sap_lines(idx: usize, sg: &SubstanceGroup) -> String {
    const PER_LINE: usize = 6;
    let mut res = String::new();
    for chunk in sg.attach_points().chunks(PER_LINE) {
        let entries: String = chunk
            .iter()
            .map(|sap| {
                // lvIdx == -1 becomes 0, which is right (see the spec)
                int_field(i64::from(sap.a_idx) + 1)
                    + &int_field(i64::from(sap.lv_idx) + 1)
                    + &string_field(&sap.id, 2, false, true)
            })
            .collect();
        // upstream formats the entry count of full lines as a 3 wide field
        let count = if chunk.len() == PER_LINE {
            int_field(PER_LINE as i64)
        } else {
            num_entries_field(chunk.len())
        };
        res += &format!("M  SAP{}{count}{entries}\n", int_field(idx as i64));
    }
    res
}

/// `GetMolFileSGroupInfo`: the V2000 lines of all groups of the molecule.
pub(super) fn mol_file_sgroup_info(mol: &ROMol) -> Result<String> {
    let groups = mol.substance_groups();
    let mut res = String::new();
    res += &entry_lines(
        "STY",
        8,
        group_entries(groups, |i, sg| {
            Some(int_field(i) + &string_field(&type_of(sg), 3, true, true))
        }),
    );
    res += &entry_lines(
        "SLB",
        8,
        group_entries(groups, |i, sg| {
            uint_prop(sg, "ID").map(|id| int_field(i) + &int_field(i64::from(id)))
        }),
    );
    for (key, code) in [("SUBTYPE", "SST"), ("CONNECT", "SCN")] {
        res += &entry_lines(
            code,
            8,
            group_entries(groups, |i, sg| {
                string_prop(sg, key).map(|v| int_field(i) + &string_field(&v, 3, true, true))
            }),
        );
    }
    // the SDS lines have a fixed `EXP` after the code
    let sds = group_entries(groups, |i, sg| {
        (string_prop(sg, "ESTATE").as_deref() == Some("E")).then(|| int_field(i))
    });
    res += &entry_lines("SDS EXP", 15, sds);
    for (key, code) in [("PARENT", "SPL"), ("COMPNO", "SNC")] {
        res += &entry_lines(
            code,
            8,
            group_entries(groups, |i, sg| {
                uint_prop(sg, key).map(|v| int_field(i) + &int_field(i64::from(v)))
            }),
        );
    }
    let mut sbt = Vec::new();
    for (i, sg) in groups.iter().enumerate() {
        if let Some(kind) = string_prop(sg, "BRKTYP") {
            let code = match kind.as_str() {
                "BRACKET" => 0,
                "PAREN" => 1,
                _ => {
                    return Err(RDKitError::Runtime(format!(
                        "Invalid BRKTYP value '{kind}' for SGroup {}",
                        i + 1
                    )));
                }
            };
            sbt.push(int_field(i as i64 + 1) + &int_field(code));
        }
    }
    res += &entry_lines("SBT", 8, sbt);

    for (i, sg) in groups.iter().enumerate() {
        let idx = i + 1;
        res += &idx_vector_lines(idx, "SAL", sg.atoms());
        // the SPA line must come after the SAL line
        res += &idx_vector_lines(idx, "SPA", sg.parent_atoms());
        res += &idx_vector_lines(idx, "SBL", sg.bonds());
        res += &sdi_lines(idx, sg);
        res += &smt_line(idx, sg);
        res += &sbv_lines(idx, sg);
        res += &sdt_line(idx, sg);
        res += &single_string_line("SDD", "FIELDDISP", idx, sg);
        // SCD/SED must come after SDT
        res += &scd_sed_lines(idx, sg)?;
        res += &sap_lines(idx, sg);
        res += &single_string_line("SCL", "CLASS", idx, sg);
    }
    Ok(res)
}

// ----------------------------------------------------------------- V3000

/// `BuildV3000IdxVectorDataBlock`
fn v3000_idx_block(key: &str, indices: &[u32]) -> String {
    if indices.is_empty() {
        return String::new();
    }
    let items: String = indices.iter().map(|i| format!(" {}", i + 1)).collect();
    format!(" {key}=({}{items})", indices.len())
}

/// `BuildV3000BondsBlock`: crossing bonds first, then the contained ones.
fn v3000_bonds_block(mol: &ROMol, sg: &SubstanceGroup) -> Result<String> {
    let mut crossing = Vec::new();
    let mut contained = Vec::new();
    for &b in sg.bonds() {
        match sg.bond_type(mol, b)? {
            SGroupBondType::XBond => crossing.push(b),
            SGroupBondType::CBond => contained.push(b),
        }
    }
    let mut res = v3000_idx_block("XBONDS", &crossing) + &v3000_idx_block("CBONDS", &contained);
    for key in ["XBHEAD", "XBCORR"] {
        if let Some(RDValue::VecUnsignedInt(v)) = sg.get_prop(key) {
            res += &v3000_idx_block(key, v);
        }
    }
    Ok(res)
}

/// `FormatV3000StringPropertyBlock`: values with blanks, quotes or opening
/// parentheses are quoted, quotes inside are doubled.
fn v3000_string_block(key: &str, sg: &SubstanceGroup) -> String {
    let Some(value) = string_prop(sg, key).filter(|v| !v.is_empty()) else {
        return String::new();
    };
    let quote = value.contains([' ', '"', '(']);
    let body = value.replace('"', "\"\"");
    if quote {
        format!(" {key}=\"{body}\"")
    } else {
        format!(" {key}={body}")
    }
}

fn v3000_uint_block(key: &str, sg: &SubstanceGroup) -> String {
    uint_prop(sg, key).map_or_else(String::new, |v| format!(" {key}={v}"))
}

fn v3000_bracket_block(sg: &SubstanceGroup) -> String {
    sg.brackets()
        .iter()
        .map(|b| {
            // z is 0 and the third point is 0 by format specification
            let pts: String = b[..2]
                .iter()
                .map(|p| format!(" {} {} 0", v3000_double(p.x), v3000_double(p.y)))
                .collect();
            format!(" BRKXYZ=(9{pts} 0 0 0)")
        })
        .collect()
}

fn v3000_field_data_block(sg: &SubstanceGroup) -> String {
    match sg.get_prop("DATAFIELDS") {
        Some(RDValue::VecString(fields)) => fields
            .iter()
            .map(|d| format!(" FIELDDATA=\"{d}\""))
            .collect(),
        _ => String::new(),
    }
}

fn v3000_cstate_block(sg: &SubstanceGroup) -> String {
    let is_sup = type_of(sg) == "SUP";
    sg.cstates()
        .iter()
        .map(|c| {
            let bond = c.bond_idx + 1;
            if is_sup {
                format!(
                    " CSTATE=(4 {bond} {} {} 0)",
                    v3000_double(c.vector.x),
                    v3000_double(c.vector.y)
                )
            } else {
                format!(" CSTATE=(1 {bond})")
            }
        })
        .collect()
}

fn v3000_attach_point_block(sg: &SubstanceGroup) -> String {
    sg.attach_points()
        .iter()
        .map(|sap| {
            let leaving = if sap.lv_idx != -1 && i64::from(sap.a_idx) == i64::from(sap.lv_idx) {
                "aidx".to_string()
            } else {
                // lvIdx == -1 becomes 0, which is right (see the spec)
                (i64::from(sap.lv_idx) + 1).to_string()
            };
            format!(" SAP=(3 {} {leaving} {})", sap.a_idx + 1, sap.id)
        })
        .collect()
}

/// `addBlockToSGroupString`: appends a block to the current line, wrapping
/// at 78 columns with the `-` continuation of V3000 lines.
fn add_block(block: &str, current: &mut Vec<u8>, out: &mut Vec<u8>) {
    let block = block.as_bytes();
    if block.is_empty() {
        return;
    }
    if current.len() + block.len() < 78 {
        current.extend_from_slice(block);
        return;
    }
    out.extend_from_slice(current);
    out.extend_from_slice(b" -\n");
    let length = block.len();
    let mut start = 0;
    while length - start >= 73 {
        out.extend_from_slice(b"M  V30");
        if start > 0 {
            out.push(b' ');
        }
        out.extend_from_slice(&block[start..start + 72]);
        start += 72;
        if start < length {
            out.extend_from_slice(b"-\n");
        }
    }
    if start < length {
        *current = b"M  V30".to_vec();
        if start > 0 {
            current.push(b' ');
        }
        current.extend_from_slice(&block[start..length.min(start + 73)]);
    }
}

/// `GetV3000MolFileSGroupLines`: the lines of one group, with the labels in
/// the order of the spec.
pub(super) fn v3000_sgroup_lines(mol: &ROMol, idx: usize, sg: &SubstanceGroup) -> Result<String> {
    let id = uint_prop(sg, "ID").unwrap_or(0);
    let mut current = format!("M  V30 {idx} {} {id}", type_of(sg)).into_bytes();
    let mut out = Vec::new();
    let blocks = [
        v3000_idx_block("ATOMS", sg.atoms()),
        // also writes XBHEAD and XBCORR
        v3000_bonds_block(mol, sg)?,
        v3000_idx_block("PATOMS", sg.parent_atoms()),
        v3000_string_block("SUBTYPE", sg),
        v3000_string_block("MULT", sg),
        v3000_string_block("CONNECT", sg),
        v3000_uint_block("PARENT", sg),
        v3000_uint_block("COMPNO", sg),
        v3000_string_block("LABEL", sg),
        v3000_bracket_block(sg),
        v3000_string_block("ESTATE", sg),
        v3000_cstate_block(sg),
        v3000_string_block("FIELDNAME", sg),
        v3000_string_block("FIELDINFO", sg),
        v3000_string_block("FIELDDISP", sg),
        v3000_string_block("QUERYTYPE", sg),
        v3000_string_block("QUERYOP", sg),
        v3000_field_data_block(sg),
        v3000_string_block("CLASS", sg),
        v3000_attach_point_block(sg),
        v3000_string_block("BRKTYP", sg),
        v3000_string_block("SEQID", sg),
    ];
    for block in &blocks {
        add_block(block, &mut current, &mut out);
    }
    if !current.is_empty() && current != b"M  V30" {
        out.extend_from_slice(&current);
        out.push(b'\n');
    } else {
        // there is an extra " -" at the end
        out.truncate(out.len().saturating_sub(3));
        out.push(b'\n');
    }
    Ok(String::from_utf8_lossy(&out).into_owned())
}

/// The `BEGIN SGROUP` block of a V3000 CTAB (empty without groups).
pub(super) fn v3000_sgroup_block(mol: &ROMol) -> Result<String> {
    if mol.substance_groups().is_empty() {
        return Ok(String::new());
    }
    let mut res = String::from("M  V30 BEGIN SGROUP\n");
    for (i, sg) in mol.substance_groups().iter().enumerate() {
        res += &v3000_sgroup_lines(mol, i + 1, sg)?;
    }
    res += "M  V30 END SGROUP\n";
    Ok(res)
}

// ------------------------------------------------- additional properties

fn is_redundant_query(q: &crate::graphmol::queries::AtomQuery) -> bool {
    q.description == "AtomAnd"
        && q.children.len() == 2
        && q.children[0].description == "AtomAtomicNum"
        && !q.children[0].negate
        && !q.children[1].negate
        && matches!(
            q.children[1].description.as_str(),
            "AtomIsotope" | "AtomFormalCharge"
        )
}

/// `createSMARTSQSubstanceGroups`: DAT groups hold the SMARTS of query atoms
/// whose query a mol file cannot express otherwise.
fn create_smartsq_groups(mol: &mut RWMol) -> Result<()> {
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom(idx);
        let Some(query) = atom.query() else {
            continue;
        };
        let marvin = atom.get_prop(cp::MRV_SMA).map(RDValue::to_rd_string);
        let needs_smarts = marvin.is_none()
            && !is_list_query(mol, idx)
            && query.description != "AtomNull"
            && !query.description.starts_with("AtomType")
            && !query.description.starts_with("AtomAtomicNum")
            && !is_redundant_query(query);
        let sma = match marvin {
            Some(sma) => sma,
            None if needs_smarts => get_atom_smarts(mol, idx, &SmilesWriteParams::default())?,
            None => String::new(),
        };
        if sma.is_empty() {
            continue;
        }
        let mut sg = SubstanceGroup::new("DAT");
        sg.set_prop("QUERYTYPE", "SMARTSQ", false)?;
        sg.set_prop("QUERYOP", "=", false)?;
        sg.set_prop("DATAFIELDS", vec![sma], false)?;
        sg.add_atom_with_idx(mol, idx)?;
        add_substance_group(mol, sg);
    }
    Ok(())
}

/// `createZBOSubstanceGroups`: zero-order bonds, with the hydrogen counts
/// and charges of their atoms.
fn create_zbo_groups(mol: &mut RWMol) -> Result<()> {
    let mut bonds = SubstanceGroup::new("DAT");
    bonds.set_prop("FIELDNAME", "ZBO", false)?;
    let mut affected = vec![false; mol.num_atoms() as usize];
    for bond in mol
        .bonds()
        .iter()
        .filter(|b| b.bond_type() == BondType::Zero)
    {
        bonds.add_bond_with_idx(mol, bond.idx())?;
        affected[bond.begin_atom_idx() as usize] = true;
        affected[bond.end_atom_idx() as usize] = true;
    }
    if !affected.contains(&true) {
        return Ok(());
    }
    let mut hydrogens = SubstanceGroup::new("DAT");
    hydrogens.set_prop("FIELDNAME", "HYD", false)?;
    let mut charges = SubstanceGroup::new("DAT");
    charges.set_prop("FIELDNAME", "ZCH", false)?;
    let (mut hyd_text, mut zch_text) = (Vec::new(), Vec::new());
    for idx in (0..mol.num_atoms()).filter(|&i| affected[i as usize]) {
        bonds.add_atom_with_idx(mol, idx)?;
        hydrogens.add_atom_with_idx(mol, idx)?;
        charges.add_atom_with_idx(mol, idx)?;
        hyd_text.push(mol.total_num_hs(idx, false)?.to_string());
        zch_text.push(mol.atom(idx).formal_charge().to_string());
    }
    hydrogens.set_prop("DATAFIELDS", vec![hyd_text.join(";")], false)?;
    charges.set_prop("DATAFIELDS", vec![zch_text.join(";")], false)?;
    for sg in [bonds, hydrogens, charges] {
        add_substance_group(mol, sg);
    }
    Ok(())
}

/// `moveAdditionalPropertiesToSGroups`
pub(super) fn move_additional_properties_to_sgroups(mol: &mut RWMol) -> Result<()> {
    if mol
        .atoms()
        .iter()
        .any(|a| a.has_prop(cp::_QueryAtomGenericLabel))
    {
        return Err(RDKitError::Runtime(
            "writing generic group queries as SubstanceGroups is not supported by rdoxide yet"
                .into(),
        ));
    }
    create_smartsq_groups(mol)?;
    create_zbo_groups(mol)
}
