//! SubstanceGroup parsing of MDL mol blocks: the V2000 `M  S..` lines and the
//! V3000 `BEGIN SGROUP` block.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/MolSGroupParsing.{h,cpp}`
//! (2026.09.1) and the SGroup post-processing of `MolFileParser.cpp`.
//!
//! Copyright (C) 2002-2018 Greg Landrum and T5 Informatics GmbH.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Groups are collected in a map keyed by their sequence number while the
//! block is read (V2000 spreads one group over many lines, in any order) and
//! are added to the molecule at the end. A group that cannot be read is
//! marked invalid; strict parsing turns the problem into an error, lenient
//! parsing logs a warning and drops the group.

use std::collections::BTreeMap;

use super::mol_file_parser::fpe;
use super::mol_file_v3000::get_v3000_line;
use super::utils::{LineStream, lossy, substr, to_double_field, to_int, to_unsigned};
use crate::geometry::Point3D;
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::queries::make_recursive_structure_query;
use crate::graphmol::query_ops::query_atom_from_atom;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substance_group::{SubstanceGroup, add_substance_group, checks};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_warning};

/// `IDX_TO_SGROUP_MAP`
pub(super) type SGroupMap = BTreeMap<i32, SubstanceGroup>;

/// `SGroupWarnOrThrow`: an error in strict mode, a warning otherwise.
fn warn_or_throw(strict: bool, msg: impl Into<String>) -> Result<()> {
    let msg = msg.into();
    if strict {
        return Err(fpe(msg));
    }
    rd_log_warning!("{msg}");
    Ok(())
}

/// `SGroupWarnOrThrow<MolFileUnhandledFeatureException>`
fn warn_or_throw_unhandled(strict: bool, msg: impl Into<String>) -> Result<()> {
    let msg = msg.into();
    if strict {
        return Err(RDKitError::UnhandledFeature(msg));
    }
    rd_log_warning!("{msg}");
    Ok(())
}

/// The `ok`-flag overloads of the field parsers: in strict mode the error
/// is raised, otherwise it is logged and the line is abandoned (`None`).
fn lenient<T>(r: Result<T>, strict: bool) -> Result<Option<T>> {
    match r {
        Ok(v) => Ok(Some(v)),
        Err(e) if strict => Err(e),
        Err(e) => {
            rd_log_warning!("{e}");
            Ok(None)
        }
    }
}

fn trim_right(s: &[u8]) -> &[u8] {
    let end = s
        .iter()
        .rposition(|c| !c.is_ascii_whitespace())
        .map_or(0, |p| p + 1);
    &s[..end]
}

fn type_of(sg: &SubstanceGroup) -> String {
    sg.get_prop("TYPE")
        .map(RDValue::to_rd_string)
        .unwrap_or_default()
}

fn set_index_prop(sg: &mut SubstanceGroup, key: &str, value: u32) -> Result<()> {
    sg.set_prop(key, value, false)
}

/// `FindSgIdx`
fn find_sg(map: &mut SGroupMap, idx: u32, line: u32) -> Option<&mut SubstanceGroup> {
    let found = map.get_mut(&(idx as i32));
    if found.is_none() {
        rd_log_warning!("SGroup {idx} referenced on line {line} not found.");
    }
    found
}

// ---------------------------------------------------------------- V2000

/// Reads the fixed-width fields of an `M  S..` line.
struct Fields<'a> {
    text: &'a [u8],
    line: u32,
    strict: bool,
    pos: usize,
}

impl<'a> Fields<'a> {
    fn new(text: &'a [u8], line: u32, strict: bool, pos: usize) -> Self {
        Fields {
            text,
            line,
            strict,
            pos,
        }
    }

    fn too_short(&self, what: &str) -> RDKitError {
        fpe(format!(
            "SGroup {what}line too short: '{}' on line {}",
            lossy(self.text),
            self.line
        ))
    }

    /// `ParseSGroupIntField` (throwing version); a counter is a smaller field.
    fn int_strict(&mut self, is_counter: bool) -> Result<u32> {
        self.pos += 1; // the separating space
        let len = 3 - usize::from(is_counter);
        let value = substr(self.text, self.pos, len)
            .map_err(|_| self.too_short(""))
            .and_then(|field| {
                to_int(field, true).map_err(|_| {
                    fpe(format!(
                        "Cannot convert '{}' to int on line {}",
                        lossy(field),
                        self.line
                    ))
                })
            })?;
        self.pos += len;
        Ok(value as u32)
    }

    /// `ParseSGroupIntField` (the `ok` overload)
    fn int(&mut self, is_counter: bool) -> Result<Option<u32>> {
        let r = self.int_strict(is_counter);
        lenient(r, self.strict)
    }

    /// `ParseSGroupDoubleField` (the `ok` overload)
    fn double(&mut self) -> Result<Option<f64>> {
        let r = to_double_field(self.text, self.pos, 10)
            .map_err(|e| match e {
                RDKitError::Index(_) => self.too_short(""),
                _ => fpe(format!(
                    "Cannot convert '{}' to double on line {}",
                    lossy(&self.text[self.pos.min(self.text.len())..]),
                    self.line
                )),
            })
            .inspect(|_| self.pos += 10);
        lenient(r, self.strict)
    }
}

/// The groups of a V2000 block while its lines are read.
#[derive(Default)]
pub(super) struct V2000SGroups {
    groups: SGroupMap,
    data_fields: BTreeMap<i32, Vec<String>>,
    scd_counter: u32,
    last_data_sgroup: u32,
    current_data_field: Vec<u8>,
}

/// The line prefixes [`V2000SGroups::parse_line`] understands.
pub(super) fn is_sgroup_line(line_beg: &[u8]) -> bool {
    matches!(
        line_beg,
        b"M  STY"
            | b"M  SST"
            | b"M  SLB"
            | b"M  SCN"
            | b"M  SDS"
            | b"M  SAL"
            | b"M  SBL"
            | b"M  SPA"
            | b"M  SMT"
            | b"M  SDI"
            | b"M  SBV"
            | b"M  SDT"
            | b"M  SDD"
            | b"M  SCD"
            | b"M  SED"
            | b"M  SPL"
            | b"M  SNC"
            | b"M  SAP"
            | b"M  SCL"
            | b"M  SBT"
    )
}

impl V2000SGroups {
    /// Parses one `M  S..` line (`is_sgroup_line(&text[..6])`).
    pub(super) fn parse_line(
        &mut self,
        mol: &RWMol,
        text: &[u8],
        line: u32,
        strict: bool,
    ) -> Result<()> {
        let mut f = Fields::new(text, line, strict, 6);
        match &text[..6] {
            b"M  STY" => self.sty(&mut f),
            b"M  SST" => self.per_group(&mut f, "SST ", 8, |sg, f| {
                let sub_type = lossy(substr(f.text, f.pos + 1, 3)?);
                f.pos += 1;
                if !checks::is_valid_sub_type(&sub_type) {
                    warn_or_throw(
                        f.strict,
                        format!("Unsupported SGroup subtype '{sub_type}' on line {}", f.line),
                    )?;
                    sg.set_is_valid(false);
                    return Ok(false);
                }
                sg.set_prop("SUBTYPE", sub_type, false)?;
                f.pos += 3;
                Ok(true)
            }),
            b"M  SLB" => self.per_group(&mut f, "SLB ", 8, |sg, f| {
                let Some(id) = f.int(false)? else {
                    sg.set_is_valid(false);
                    return Ok(false);
                };
                if id != 0 && !checks::is_substance_group_id_free(mol, id) {
                    warn_or_throw(
                        f.strict,
                        format!(
                            "SGroup ID '{id}' is assigned to more than one SGroup, on line {}",
                            f.line
                        ),
                    )?;
                    sg.set_is_valid(false);
                    return Ok(false);
                }
                set_index_prop(sg, "ID", id)?;
                Ok(true)
            }),
            b"M  SCN" => self.per_group(&mut f, "SCN ", 7, |sg, f| {
                let connect = lossy(substr(f.text, f.pos + 1, 2)?);
                f.pos += 1;
                if !checks::is_valid_connect_type(&connect) {
                    warn_or_throw(
                        f.strict,
                        format!(
                            "Unsupported SGroup connection type '{connect}' on line {}",
                            f.line
                        ),
                    )?;
                    sg.set_is_valid(false);
                    return Ok(false);
                }
                sg.set_prop("CONNECT", connect, false)?;
                f.pos += 3;
                Ok(true)
            }),
            b"M  SDS" => {
                let mut f = Fields::new(text, line, strict, 10);
                self.per_group(&mut f, "SDS ", 4, |sg, _| {
                    sg.set_prop("ESTATE", "E", false)?;
                    Ok(true)
                })
            }
            b"M  SPL" => self.per_group(&mut f, "SPL ", 8, |sg, f| {
                let parent = f.int_strict(false)?;
                set_index_prop(sg, "PARENT", parent)?;
                Ok(true)
            }),
            b"M  SNC" => self.per_group(&mut f, "SNC ", 8, |sg, f| {
                let Some(compno) = f.int(false)? else {
                    sg.set_is_valid(false);
                    return Ok(false);
                };
                if compno > 256 {
                    warn_or_throw(
                        f.strict,
                        format!("SGroup SNC value over 256: '{compno}' on line {}", f.line),
                    )?;
                    sg.set_is_valid(false);
                    return Ok(false);
                }
                set_index_prop(sg, "COMPNO", compno)?;
                Ok(true)
            }),
            b"M  SBT" => self.per_group(&mut f, "SBT ", 8, |sg, f| {
                let Some(bracket_type) = f.int(false)? else {
                    sg.set_is_valid(false);
                    return Ok(false);
                };
                let name = match bracket_type {
                    0 => "BRACKET",
                    1 => "PAREN",
                    _ => {
                        warn_or_throw(
                            f.strict,
                            format!("Invalid SBT value '{bracket_type}' on line {}", f.line),
                        )?;
                        sg.set_is_valid(false);
                        return Ok(false);
                    }
                };
                sg.set_prop("BRKTYP", name, false)?;
                Ok(true)
            }),
            b"M  SAL" | b"M  SBL" | b"M  SPA" => self.vector_data(mol, &mut f),
            b"M  SDI" => self.sdi(&mut f),
            b"M  SMT" => self.smt(&mut f),
            b"M  SBV" => self.sbv(mol, &mut f),
            b"M  SDT" => self.sdt(&mut f),
            b"M  SDD" => self.sdd(&mut f),
            b"M  SCD" | b"M  SED" => self.scd_sed(&mut f),
            b"M  SAP" => self.sap(mol, &mut f),
            b"M  SCL" => self.scl(&mut f),
            _ => unreachable!("not an SGroup line"),
        }
    }

    /// `M  STY`: creates the groups.
    fn sty(&mut self, f: &mut Fields<'_>) -> Result<()> {
        let Some(nent) = f.int(true)? else {
            return Ok(());
        };
        for _ in 0..nent {
            if f.text.len() < f.pos + 8 {
                return warn_or_throw(f.strict, f.too_short("STY ").to_string());
            }
            let Some(sequence_id) = f.int(false)? else {
                return Ok(());
            };
            let typ = lossy(substr(f.text, f.pos + 1, 3)?);
            if checks::is_valid_type(&typ) {
                let mut sg = SubstanceGroup::new(&typ);
                set_index_prop(&mut sg, "index", sequence_id)?;
                self.groups.entry(sequence_id as i32).or_insert(sg);
            } else {
                warn_or_throw_unhandled(f.strict, format!("S group {typ} on line {}", f.line))?;
            }
            f.pos += 4;
        }
        Ok(())
    }

    /// Lines that start with a count followed by `count` entries which each
    /// begin with a group number: `apply` handles one entry and returns
    /// whether to go on.
    fn per_group(
        &mut self,
        f: &mut Fields<'_>,
        name: &str,
        min_remaining: usize,
        mut apply: impl FnMut(&mut SubstanceGroup, &mut Fields<'_>) -> Result<bool>,
    ) -> Result<()> {
        let Some(nent) = f.int(true)? else {
            return Ok(());
        };
        for _ in 0..nent {
            if f.text.len() < f.pos + min_remaining {
                let mut msg = f.too_short(name).to_string();
                if name == "SCN " {
                    msg += &format!("\n needed: {} found: {}", f.pos + 7, f.text.len());
                }
                return warn_or_throw(f.strict, msg);
            }
            let Some(sg_idx) = f.int(false)? else {
                return Ok(());
            };
            let Some(sg) = find_sg(&mut self.groups, sg_idx, f.line) else {
                return Ok(());
            };
            if !apply(sg, f)? {
                return Ok(());
            }
        }
        Ok(())
    }

    /// The group number at the start of a line that describes one group.
    fn group_of<'s>(&'s mut self, f: &mut Fields<'_>) -> Result<Option<&'s mut SubstanceGroup>> {
        let Some(sg_idx) = f.int(false)? else {
            return Ok(None);
        };
        Ok(find_sg(&mut self.groups, sg_idx, f.line))
    }

    /// `M  SAL`, `M  SBL`, `M  SPA`
    fn vector_data(&mut self, mol: &RWMol, f: &mut Fields<'_>) -> Result<()> {
        let kind = f.text[3..6].to_vec();
        let Some(sg) = self.group_of(f)? else {
            return Ok(());
        };
        let Some(nent) = f.int(true)? else {
            sg.set_is_valid(false);
            return Ok(());
        };
        for _ in 0..nent {
            if f.text.len() < f.pos + 4 {
                warn_or_throw(f.strict, f.too_short("").to_string())?;
                sg.set_is_valid(false);
                return Ok(());
            }
            let Some(nbr) = f.int(false)? else {
                sg.set_is_valid(false);
                return Ok(());
            };
            let added = match kind.as_slice() {
                b"SAL" => sg.add_atom_with_bookmark(mol, nbr as i32),
                b"SBL" => sg.add_bond_with_bookmark(mol, nbr as i32),
                _ => sg.add_parent_atom_with_bookmark(mol, nbr as i32),
            };
            if let Err(e) = added {
                warn_or_throw(f.strict, e.to_string())?;
                sg.set_is_valid(false);
                return Ok(());
            }
        }
        Ok(())
    }

    /// `M  SDI`: a bracket.
    fn sdi(&mut self, f: &mut Fields<'_>) -> Result<()> {
        let Some(sg) = self.group_of(f)? else {
            return Ok(());
        };
        let Some(n_coords) = f.int(true)? else {
            sg.set_is_valid(false);
            return Ok(());
        };
        if n_coords != 4 {
            warn_or_throw(
                f.strict,
                format!(
                    "Unexpected number of coordinates for SDI on line {}",
                    f.line
                ),
            )?;
            sg.set_is_valid(false);
            return Ok(());
        }
        let mut bracket = [Point3D::new(0.0, 0.0, 0.0); 3];
        for point in &mut bracket[..2] {
            let (Some(x), Some(y)) = (f.double()?, f.double()?) else {
                sg.set_is_valid(false);
                return Ok(());
            };
            *point = Point3D::new(x, y, 0.0);
        }
        sg.add_bracket(bracket);
        Ok(())
    }

    /// `M  SMT`: the label (or the multiplier of a MUL group).
    fn smt(&mut self, f: &mut Fields<'_>) -> Result<()> {
        let Some(sg) = self.group_of(f)? else {
            return Ok(());
        };
        f.pos += 1;
        if f.pos >= f.text.len() {
            warn_or_throw(f.strict, f.too_short("").to_string())?;
            sg.set_is_valid(false);
            return Ok(());
        }
        let label = lossy(&f.text[f.pos..]);
        // abbreviation groups might not have seen their SCL line yet
        let key = if type_of(sg) == "MUL" {
            "MULT"
        } else {
            "LABEL"
        };
        sg.set_prop(key, label, false)
    }

    /// `M  SBV`: the vector of a crossing bond.
    fn sbv(&mut self, mol: &RWMol, f: &mut Fields<'_>) -> Result<()> {
        let Some(sg) = self.group_of(f)? else {
            return Ok(());
        };
        let Some(bond_mark) = f.int(false)? else {
            sg.set_is_valid(false);
            return Ok(());
        };
        let bond = mol.get_unique_bond_with_bookmark(bond_mark as i32)?;
        let mut vector = Point3D::new(0.0, 0.0, 0.0);
        if type_of(sg) == "SUP" {
            let (Some(x), Some(y)) = (f.double()?, f.double()?) else {
                sg.set_is_valid(false);
                return Ok(());
            };
            vector = Point3D::new(x, y, 0.0);
        }
        if let Err(e) = sg.add_cstate(mol, bond, vector) {
            warn_or_throw(f.strict, e.to_string())?;
            sg.set_is_valid(false);
        }
        Ok(())
    }

    /// `M  SDT`: the data field definition. Short lines are accepted: what
    /// could be read is used.
    fn sdt(&mut self, f: &mut Fields<'_>) -> Result<()> {
        let Some(sg) = self.group_of(f)? else {
            return Ok(());
        };
        let mut values: [Vec<u8>; 5] = Default::default();
        let widths = [30, 2, 20, 2, usize::MAX];
        f.pos += 1;
        for (value, width) in values.iter_mut().zip(widths) {
            let Ok(field) = substr(f.text, f.pos, width) else {
                break;
            };
            *value = trim_right(field).to_vec();
            f.pos = f.pos.saturating_add(width);
        }
        let keys = [
            "FIELDNAME",
            "FIELDTYPE",
            "FIELDINFO",
            "QUERYTYPE",
            "QUERYOP",
        ];
        for (key, value) in keys.into_iter().zip(values) {
            // only add entries for the properties that are not blank
            if !value.is_empty() {
                sg.set_prop(key, lossy(&value), false)?;
            }
        }
        Ok(())
    }

    /// `M  SDD`: the field display, stored as is.
    fn sdd(&mut self, f: &mut Fields<'_>) -> Result<()> {
        let Some(sg) = self.group_of(f)? else {
            return Ok(());
        };
        f.pos += 1;
        if f.pos < f.text.len() {
            sg.set_prop("FIELDDISP", lossy(&f.text[f.pos..]), false)?;
        }
        Ok(())
    }

    /// `M  SCD` / `M  SED`: the lines of a data field.
    fn scd_sed(&mut self, f: &mut Fields<'_>) -> Result<()> {
        let kind = f.text[3..6].to_vec();
        f.pos = 6;
        let Some(sg_idx) = f.int(false)? else {
            return Ok(());
        };
        let Some(sg) = find_sg(&mut self.groups, sg_idx, f.line) else {
            return Ok(());
        };
        if self.last_data_sgroup != 0 && self.last_data_sgroup != sg_idx {
            warn_or_throw(
                f.strict,
                format!(
                    "Found a Data Field not matching the SGroup of the last Data Field at line {}",
                    f.line
                ),
            )?;
            sg.set_is_valid(false);
            return Ok(());
        } else if self.last_data_sgroup == 0 && kind == b"SCD" {
            self.last_data_sgroup = sg_idx;
        } else if kind == b"SED" {
            self.last_data_sgroup = 0;
        }
        if !sg.has_prop("FIELDNAME") {
            // the docs can be read both ways, and CTABs with empty fieldnames
            // in SDT lines exist, so this is only a warning
            rd_log_warning!(
                "Found a SCD/SED line with missing/empty SDT specification at line {}",
                f.line
            );
        }
        if f.strict && kind == b"SCD" && self.scd_counter > 2 {
            return Err(fpe(format!(
                "Found too many consecutive SCD lines, (#{} at line {}) for SGroup {sg_idx}",
                self.scd_counter + 1,
                f.line
            )));
        }
        if f.pos + 1 < f.text.len() {
            f.pos += 1;
            self.current_data_field
                .extend_from_slice(substr(f.text, f.pos, 69)?);
            if kind == b"SED" {
                let data = trim_right(&self.current_data_field);
                let data = lossy(&data[..data.len().min(200)]);
                self.data_fields
                    .entry(sg_idx as i32)
                    .or_default()
                    .push(data);
                self.current_data_field.clear();
                self.scd_counter = 0;
            } else {
                self.scd_counter += 1;
            }
        }
        Ok(())
    }

    /// `M  SAP`: attachment points. A malformed line without the leaving
    /// atom and id gets `lvIdx = numAtoms`, which [`Self::finish`] repairs.
    fn sap(&mut self, mol: &RWMol, f: &mut Fields<'_>) -> Result<()> {
        let Some(sg) = self.group_of(f)? else {
            return Ok(());
        };
        let Some(nent) = f.int(true)? else {
            sg.set_is_valid(false);
            return Ok(());
        };
        for _ in 0..nent {
            let mut lv_idx = -1;
            if f.text.len() < f.pos + 11 {
                let msg = f.too_short("SAP ").to_string();
                if f.strict {
                    return Err(fpe(msg));
                }
                rd_log_warning!("{msg}");
                if f.text.len() < f.pos + 4 {
                    sg.set_is_valid(false);
                    return Ok(());
                }
                lv_idx = mol.num_atoms() as i32;
            }
            let mut id = "  ".to_string();
            let Some(a_idx_mark) = f.int(false)? else {
                sg.set_is_valid(false);
                return Ok(());
            };
            let a_idx = mol.get_atom_with_bookmark(a_idx_mark as i32)?;
            if lv_idx == -1 {
                let Some(lv_idx_mark) = f.int(false)? else {
                    sg.set_is_valid(false);
                    return Ok(());
                };
                if lv_idx_mark != 0 {
                    lv_idx = mol.get_atom_with_bookmark(lv_idx_mark as i32)? as i32;
                }
                if f.text.len() >= f.pos + 3 {
                    id = lossy(substr(f.text, f.pos + 1, 2)?);
                    f.pos += 3;
                }
            }
            sg.add_attach_point(a_idx, lv_idx, &id);
        }
        Ok(())
    }

    /// `M  SCL`: the template class.
    fn scl(&mut self, f: &mut Fields<'_>) -> Result<()> {
        let Some(sg) = self.group_of(f)? else {
            return Ok(());
        };
        if f.pos + 1 >= f.text.len() {
            warn_or_throw(f.strict, f.too_short("SCL ").to_string())?;
            sg.set_is_valid(false);
            return Ok(());
        }
        f.pos += 1;
        sg.set_prop("CLASS", lossy(&f.text[f.pos..]), false)
    }

    /// The final updates once `M  END` is reached: the data fields are
    /// attached, attachment points checked and the valid groups added to the
    /// molecule.
    pub(super) fn finish(mut self, mol: &mut RWMol, strict: bool) -> Result<()> {
        for (idx, mut sg) in std::mem::take(&mut self.groups) {
            if sg.is_valid() {
                let fields = self.data_fields.remove(&idx).unwrap_or_default();
                sg.set_prop("DATAFIELDS", fields, false)?;
                let valid = attachment_points_are_valid(mol, &mut sg, idx);
                sg.set_is_valid(valid);
            }
            if sg.is_valid() {
                add_substance_group(mol, sg);
            } else {
                let msg = format!("SGroup {idx} is invalid");
                if strict {
                    return Err(fpe(msg));
                }
                rd_log_warning!("{msg} and will be ignored");
            }
        }
        Ok(())
    }
}

/// `checkAttachmentPointsAreValid`: infers the leaving atom of attachment
/// points of malformed `SAP` lines when the group has a single bond.
fn attachment_points_are_valid(mol: &ROMol, sg: &mut SubstanceGroup, idx: i32) -> bool {
    let n_atoms = mol.num_atoms() as i32;
    let single_bond = match sg.bonds() {
        [b] => Some(mol.bond(*b)),
        _ => None,
    };
    let mut valid = true;
    let mut fixes = Vec::new();
    for (i, ap) in sg.attach_points().iter().enumerate() {
        let mut lv_idx = ap.lv_idx;
        if lv_idx == n_atoms
            && let Some(bond) = single_bond
            && (bond.begin_atom_idx() == ap.a_idx || bond.end_atom_idx() == ap.a_idx)
            && let Ok(other) = bond.other_atom_idx(ap.a_idx)
        {
            lv_idx = other as i32;
        }
        if lv_idx == n_atoms {
            rd_log_warning!("Could not infer missing lvIdx on malformed SAP line for SGroup {idx}");
            valid = false;
        }
        fixes.push((i, lv_idx));
    }
    for (i, lv_idx) in fixes {
        sg.attach_points_mut()[i].lv_idx = lv_idx;
    }
    valid
}

// ---------------------------------------------------------------- V3000

/// The parts of `std::stringstream` the V3000 SGroup parser relies on,
/// including its sticky `eof`/`fail` flags.
struct SStream {
    data: Vec<u8>,
    pos: usize,
    eof: bool,
    fail: bool,
}

impl SStream {
    fn new(data: &[u8]) -> Self {
        SStream {
            data: data.to_vec(),
            pos: 0,
            eof: false,
            fail: false,
        }
    }

    fn reset(&mut self, data: &[u8]) {
        *self = SStream::new(data);
    }

    fn is_good(&self) -> bool {
        !self.eof && !self.fail
    }

    /// The sentry of the formatted and unformatted input functions: they do
    /// nothing (but set `fail`) on a stream that is not good.
    fn sentry(&mut self) -> bool {
        if !self.is_good() {
            self.fail = true;
        }
        self.is_good()
    }

    /// `get()`: the next character; at the end both flags are set.
    fn get(&mut self) -> Option<u8> {
        if !self.sentry() {
            return None;
        }
        match self.data.get(self.pos) {
            Some(&c) => {
                self.pos += 1;
                Some(c)
            }
            None => {
                self.eof = true;
                self.fail = true;
                None
            }
        }
    }

    /// `peek()`: only the end-of-file flag is set at the end.
    fn peek(&mut self) -> Option<u8> {
        if !self.is_good() {
            return None;
        }
        let c = self.data.get(self.pos).copied();
        if c.is_none() {
            self.eof = true;
        }
        c
    }

    /// `putback(c)`: a `stringstream` replaces the character it backs up over.
    fn putback(&mut self, c: u8) {
        self.eof = false;
        if self.pos > 0 {
            self.pos -= 1;
            self.data[self.pos] = c;
        }
    }

    fn skip_ws(&mut self) {
        while self.data.get(self.pos).is_some_and(u8::is_ascii_whitespace) {
            self.pos += 1;
        }
    }

    /// Extracts the longest prefix of the next token accepted by `accept`.
    fn token(&mut self, accept: impl Fn(&[u8]) -> usize) -> Option<Vec<u8>> {
        if !self.sentry() {
            return None;
        }
        self.skip_ws();
        let rest = &self.data[self.pos..];
        let len = accept(rest);
        if len == 0 {
            self.fail = true;
            if rest.is_empty() {
                self.eof = true;
            }
            return None;
        }
        let token = rest[..len].to_vec();
        self.pos += len;
        if self.pos >= self.data.len() {
            self.eof = true;
        }
        Some(token)
    }

    /// `stream >> unsigned int`: 0 when the extraction fails.
    fn read_u32(&mut self) -> u32 {
        let digits = |s: &[u8]| {
            let sign = usize::from(matches!(s.first(), Some(b'+' | b'-')));
            let n = s[sign..].iter().take_while(|c| c.is_ascii_digit()).count();
            if n == 0 { 0 } else { sign + n }
        };
        match self.token(digits) {
            Some(t) => {
                let neg = t[0] == b'-';
                let digits = t
                    .strip_prefix(b"+")
                    .or_else(|| t.strip_prefix(b"-"))
                    .unwrap_or(&t);
                let n: u64 = lossy(digits).parse().unwrap_or(u64::MAX);
                if n > u64::from(u32::MAX) {
                    self.fail = true;
                    return u32::MAX;
                }
                if neg {
                    (n as u32).wrapping_neg()
                } else {
                    n as u32
                }
            }
            None => 0,
        }
    }

    /// `stream >> double`: 0 when the extraction fails.
    fn read_f64(&mut self) -> f64 {
        let number = |s: &[u8]| {
            let mut i = usize::from(matches!(s.first(), Some(b'+' | b'-')));
            let digits = |i: &mut usize| {
                let start = *i;
                while s.get(*i).is_some_and(u8::is_ascii_digit) {
                    *i += 1;
                }
                *i - start
            };
            let mut n = digits(&mut i);
            if s.get(i) == Some(&b'.') {
                i += 1;
                n += digits(&mut i);
            }
            if n == 0 {
                return 0;
            }
            if matches!(s.get(i), Some(b'e' | b'E')) {
                let mut j = i + 1 + usize::from(matches!(s.get(i + 1), Some(b'+' | b'-')));
                if digits(&mut j) > 0 {
                    i = j;
                }
            }
            i
        };
        self.token(number)
            .map_or(0.0, |t| lossy(&t).parse().unwrap_or(0.0))
    }

    /// `stream >> std::string`
    fn read_word(&mut self) -> String {
        let word = |s: &[u8]| s.iter().take_while(|c| !c.is_ascii_whitespace()).count();
        self.token(word).map(|t| lossy(&t)).unwrap_or_default()
    }

    /// `std::getline(stream, s, delim)`
    fn getline(&mut self, delim: u8) -> String {
        if !self.sentry() {
            return String::new();
        }
        let rest = &self.data[self.pos..];
        let len = rest.iter().position(|&c| c == delim);
        let taken = &rest[..len.unwrap_or(rest.len())];
        self.pos += taken.len();
        if len.is_some() {
            self.pos += 1;
        } else {
            self.eof = true;
            if taken.is_empty() {
                self.fail = true;
            }
        }
        lossy(taken)
    }
}

/// `ParseV3000Array`: `(count v1 v2 ...)`. Problems are reported as
/// errors, as the label parser calls it in strict mode.
fn parse_array<T>(
    s: &mut SStream,
    max_count: Option<u32>,
    read: fn(&mut SStream) -> T,
) -> Result<Vec<T>> {
    if s.get() != Some(b'(') {
        rd_log_warning!("WARNING: first character of V3000 array is not '('");
    }
    let count = s.read_u32();
    if s.fail {
        return Err(fpe("invalid count value in V3000 array"));
    }
    if max_count.is_some_and(|max| count > max) {
        return Err(fpe("invalid count value"));
    }
    let mut values = Vec::new();
    for _ in 0..count {
        let value = read(s);
        if s.fail {
            return Err(fpe("invalid value in V3000 array"));
        }
        values.push(value);
    }
    if s.get() != Some(b')') {
        rd_log_warning!("WARNING: final character of V3000 array is not ')'");
    }
    Ok(values)
}

/// `ParseV3000CStateLabel`
fn parse_cstate_label(
    mol: &RWMol,
    sg: &mut SubstanceGroup,
    s: &mut SStream,
    line: u32,
) -> Result<()> {
    if s.get() != Some(b'(') {
        return Err(fpe(format!(
            "Missing parentheses for CSTATE field on line {line}"
        )));
    }
    let count = s.read_u32();
    let bond_mark = s.read_u32();
    let is_sup = type_of(sg) == "SUP";
    if (!is_sup && count != 1) || (is_sup && count != 4) {
        return Err(fpe(format!(
            "Unexpected number of fields for CSTATE field on line {line}"
        )));
    }
    let bond = mol.get_unique_bond_with_bookmark(bond_mark as i32)?;
    let mut vector = Point3D::new(0.0, 0.0, 0.0);
    if is_sup {
        vector = Point3D::new(s.read_f64(), s.read_f64(), s.read_f64());
    }
    sg.add_cstate(mol, bond, vector)?;
    if s.get() != Some(b')') {
        rd_log_warning!("WARNING: final character of CSTATE field is not ')'");
    }
    Ok(())
}

/// `ParseV3000SAPLabel`
fn parse_sap_label(mol: &RWMol, sg: &mut SubstanceGroup, s: &mut SStream) -> Result<()> {
    if s.get() != Some(b'(') {
        return Err(fpe("Missing parentheses for SAP field"));
    }
    let _count = s.read_u32();
    let a_idx_mark = s.read_u32();
    let lv_idx_str = s.read_word();
    let mut sap_id = s.read_word();
    if s.fail || lv_idx_str.is_empty() || !sap_id.ends_with(')') {
        return Err(fpe("Invalid SAP field"));
    }
    sap_id.pop();
    let a_idx = mol.get_atom_with_bookmark(a_idx_mark as i32)?;
    let lv_idx = if lv_idx_str.eq_ignore_ascii_case("AIDX") {
        a_idx as i32
    } else {
        match to_int(lv_idx_str.as_bytes(), true)? as u32 {
            0 => -1,
            mark => mol.get_atom_with_bookmark(mark as i32)? as i32,
        }
    };
    sg.add_attach_point(a_idx, lv_idx, &sap_id);
    Ok(())
}

/// `ParseV3000StringPropLabel`: a bare, `"` quoted (`""` is a quote) or `'`
/// quoted value.
fn parse_string_prop_label(s: &mut SStream) -> String {
    let mut value = Vec::new();
    match s.peek() {
        Some(b' ') => return String::new(), // an empty value
        Some(b'"') => {
            s.get();
            while let Some(c) = s.get() {
                if c == b'"' {
                    if s.peek() != Some(b'"') {
                        break;
                    }
                    s.get(); // the second quote of a doubled one
                }
                value.push(c);
            }
            return lossy(trim_right(&value));
        }
        Some(b'\'') => {
            s.get();
            return lossy(trim_right(s.getline(b'\'').as_bytes()));
        }
        _ => {}
    }
    lossy(trim_right(s.read_word().as_bytes()))
}

/// `ParseV3000ParseLabel`: one `LABEL=value` of an SGroup line.
fn parse_label(
    label: &str,
    s: &mut SStream,
    data_fields: &mut Vec<String>,
    line: u32,
    sg: &mut SubstanceGroup,
    mol: &RWMol,
    strict: bool,
) -> Result<()> {
    let result = parse_label_value(label, s, data_fields, line, sg, mol, strict);
    if let Err(e) = result {
        warn_or_throw(strict, e.to_string())?;
        sg.set_is_valid(false);
    }
    Ok(())
}

fn parse_label_value(
    label: &str,
    s: &mut SStream,
    data_fields: &mut Vec<String>,
    line: u32,
    sg: &mut SubstanceGroup,
    mol: &RWMol,
    strict: bool,
) -> Result<()> {
    match label {
        "XBHEAD" | "XBCORR" => {
            let bonds = parse_array(s, Some(mol.num_bonds()), SStream::read_u32)?;
            let bonds: Vec<u32> = bonds.into_iter().map(|b| b.wrapping_sub(1)).collect();
            sg.set_prop(label, bonds, false)?;
        }
        "ATOMS" => {
            for mark in parse_array(s, Some(mol.num_atoms()), SStream::read_u32)? {
                sg.add_atom_with_bookmark(mol, mark as i32)?;
            }
        }
        "PATOMS" => {
            for mark in parse_array(s, Some(mol.num_atoms()), SStream::read_u32)? {
                sg.add_parent_atom_with_bookmark(mol, mark as i32)?;
            }
        }
        "CBONDS" | "XBONDS" => {
            for mark in parse_array(s, Some(mol.num_bonds()), SStream::read_u32)? {
                sg.add_bond_with_bookmark(mol, mark as i32)?;
            }
        }
        "BRKXYZ" => {
            let coords = parse_array(s, Some(9), SStream::read_f64)?;
            if coords.len() != 9 {
                return Err(fpe(format!(
                    "Unexpected number of coordinates for BRKXYZ on line {line}"
                )));
            }
            let mut bracket = [Point3D::new(0.0, 0.0, 0.0); 3];
            for (point, xyz) in bracket.iter_mut().zip(coords.chunks(3)) {
                *point = Point3D::new(xyz[0], xyz[1], xyz[2]);
            }
            sg.add_bracket(bracket);
        }
        "CSTATE" => parse_cstate_label(mol, sg, s, line)?,
        "SAP" => parse_sap_label(mol, sg, s)?,
        "PARENT" => {
            // the relationship is resolved once all the groups have been read
            if s.eof {
                return Err(fpe(format!("PARENT label not found on line {line}")));
            }
            let parent = s.read_u32();
            if s.fail {
                return Err(fpe(format!("Invalid PARENT label found on line {line}")));
            }
            set_index_prop(sg, "PARENT", parent)?;
        }
        "COMPNO" => {
            let compno = s.read_u32();
            if compno > 256 {
                return Err(fpe(format!(
                    "SGroup SNC value over 256: '{compno}' on line {line}"
                )));
            }
            set_index_prop(sg, "COMPNO", compno)?;
        }
        "FIELDDATA" => {
            let mut value = parse_string_prop_label(s);
            if strict {
                value.truncate(200);
            }
            data_fields.push(value);
        }
        _ => {
            let value = parse_string_prop_label(s);
            if label == "SUBTYPE" && !checks::is_valid_sub_type(&value) {
                return Err(fpe(format!(
                    "Unsupported SGroup subtype '{value}' on line {line}"
                )));
            } else if label == "CONNECT" && !checks::is_valid_connect_type(&value) {
                return Err(fpe(format!(
                    "Unsupported SGroup connection type '{value}' on line {line}"
                )));
            } else if label == "CLASS" && !checks::is_valid_class(&value) {
                return Err(fpe(format!(
                    "Unsupported SGroup template class '{value}' on line {line}"
                )));
            }
            // NATREPLACE is neither validated nor used
            sg.set_prop(label, value, false)?;
        }
    }
    Ok(())
}

/// Reads `LABEL=value` items from `s` until it is exhausted or `sg` becomes
/// invalid; `skip_parsed` makes the items of `parsed` (the labels already
/// given on the group's own line) be skipped, as the DEFAULT line needs.
struct LabelScan<'a> {
    line: u32,
    context: &'a str,
}

#[allow(clippy::too_many_arguments)]
fn scan_labels(
    s: &mut SStream,
    scan: &LabelScan<'_>,
    sg: &mut SubstanceGroup,
    data_fields: &mut Vec<String>,
    parsed: &mut Vec<String>,
    skip_parsed: bool,
    mol: &RWMol,
    strict: bool,
) -> Result<()> {
    while sg.is_valid() && s.is_good() {
        let Some(spacer) = s.get() else {
            continue;
        };
        if spacer != b' ' {
            let msg = format!(
                "Found character '{}' when expecting a separator (space) {}on line {}",
                spacer as char, scan.context, scan.line
            );
            warn_or_throw(strict, msg)?;
            sg.set_is_valid(false);
            continue;
        }
        let label = s.getline(b'=');
        if label.is_empty() {
            continue;
        }
        if !skip_parsed {
            parse_label(&label, s, data_fields, scan.line, sg, mol, strict)?;
            parsed.push(label);
        } else if !parsed.contains(&label) {
            parse_label(&label, s, data_fields, scan.line, sg, mol, strict)?;
        } else {
            // the group's own line wins: skip the value of the default
            match s.peek() {
                Some(b' ') => {
                    warn_or_throw(
                        strict,
                        format!("Found unexpected whitespace at DEFAULT label {label}"),
                    )?;
                    sg.set_is_valid(false);
                }
                Some(b'(') => {
                    s.getline(b')');
                    s.get();
                }
                Some(b'"') => {
                    s.get();
                    s.getline(b'"');
                }
                _ => {
                    s.getline(b' ');
                    s.putback(b' ');
                }
            }
        }
    }
    Ok(())
}

/// `ParseV3000SGroupsBlock`: reads the groups of a `BEGIN SGROUP` block,
/// adds the valid ones to the molecule and returns the last line read.
pub(super) fn parse_v3000_sgroups_block(
    s: &mut LineStream<'_>,
    line: &mut u32,
    mut n_sgroups: u32,
    mol: &mut RWMol,
    strict: bool,
) -> Result<Vec<u8>> {
    let mut default_line_num = 0;
    let mut default_string = Vec::new();
    // groups may come in any order, so they are collected before being added
    let mut groups = SGroupMap::new();

    let mut temp = get_v3000_line(s, line)?;
    if temp.starts_with(b"DEFAULT") && temp.len() > 8 {
        default_string = trim_right(&temp[7..]).to_vec();
        default_line_num = *line;
        temp = trim_right(&get_v3000_line(s, line)?).to_vec();
    }

    for _ in 0..n_sgroups {
        let mut ls = SStream::new(&temp);
        let sequence_id = ls.read_u32();
        let typ = ls.read_word();
        let external_id = ls.read_u32();

        if strict && !checks::is_valid_type(&typ) {
            return Err(RDKitError::UnhandledFeature(format!(
                "Unsupported SGroup type '{typ}' on line {line}"
            )));
        } else if !strict && n_sgroups == u32::MAX && ls.fail {
            // we did not know how many groups to expect and this does not look
            // like the start of one: the block is over
            n_sgroups = 0;
            break;
        }

        let mut sg = SubstanceGroup::new(&typ);
        let mut data_fields = Vec::new();
        set_index_prop(&mut sg, "index", sequence_id)?;
        if external_id > 0 {
            if !checks::is_substance_group_id_free(mol, external_id) {
                let msg = format!(
                    "Existing SGroup ID '{external_id}' assigned to a second SGroup on line {line}"
                );
                if strict {
                    return Err(fpe(msg));
                }
                rd_log_warning!("{msg}");
                sg.set_is_valid(false);
            }
            set_index_prop(&mut sg, "ID", external_id)?;
        }

        let mut parsed = Vec::new();
        let own = LabelScan {
            line: *line,
            context: "",
        };
        scan_labels(
            &mut ls,
            &own,
            &mut sg,
            &mut data_fields,
            &mut parsed,
            false,
            mol,
            strict,
        )?;

        ls.reset(&default_string);
        let defaults = LabelScan {
            line: default_line_num,
            context: "in DEFAULTS ",
        };
        scan_labels(
            &mut ls,
            &defaults,
            &mut sg,
            &mut data_fields,
            &mut parsed,
            true,
            mol,
            strict,
        )?;

        sg.set_prop("DATAFIELDS", data_fields, false)?;
        groups.entry(sequence_id as i32).or_insert(sg);
        temp = trim_right(&get_v3000_line(s, line)?).to_vec();
    }

    if groups.len() != n_sgroups as usize {
        warn_or_throw(
            strict,
            format!(
                "Found {} SGroups when {n_sgroups} were expected.\n",
                groups.len()
            ),
        )?;
    }
    for (idx, sg) in groups {
        if sg.is_valid() {
            add_substance_group(mol, sg);
        } else {
            rd_log_warning!("SGroup {idx} is invalid and will be ignored");
        }
    }
    Ok(temp)
}

// ------------------------------------------------- structure-changing groups

fn data_fields_of(sg: &SubstanceGroup) -> Option<Vec<String>> {
    match sg.get_prop("DATAFIELDS") {
        Some(RDValue::VecString(v)) => Some(v.clone()),
        _ => None,
    }
}

/// The `;` separated values of the data fields of ZCH and HYD groups, one
/// per atom of the group.
fn per_atom_values(sg: &SubstanceGroup, kind: &str) -> Option<Vec<Vec<i32>>> {
    let fields = data_fields_of(sg)?;
    if fields.is_empty() {
        rd_log_warning!("ignoring {kind} SGroup without data fields.");
        return None;
    }
    let mut res = Vec::new();
    for field in fields {
        let values: Vec<&str> = field.trim().split(';').collect();
        if values.len() < sg.atoms().len() {
            rd_log_warning!(
                "DATAFIELDS in {kind} SGroup is shorter than the number of atoms in the SGroup. Ignoring it."
            );
            continue;
        }
        let parsed = values[..sg.atoms().len()]
            .iter()
            .map(|v| {
                if v.is_empty() {
                    Ok(0)
                } else {
                    to_int(v.as_bytes(), true)
                }
            })
            .collect::<Result<Vec<_>>>();
        res.push(parsed.ok()?);
    }
    Some(res)
}

/// `processMrvCoordinateBond`: MarvinSketch marks coordinate bonds with a
/// DAT group that names the (unspecified) bond to turn into a dative one.
fn process_mrv_coordinate_bond(mol: &mut RWMol, sg: &SubstanceGroup) -> Result<()> {
    let Some(fields) = data_fields_of(sg) else {
        return Ok(());
    };
    let Some(first) = fields.first() else {
        rd_log_warning!("ignoring MRV_COORDINATE_BOND_TYPE SGroup without data fields.");
        return Ok(());
    };
    let bond_idx = to_unsigned(first.as_bytes(), true)?.wrapping_sub(1);
    if fields.len() > 1 {
        rd_log_warning!(
            "ignoring extra data fields in MRV_COORDINATE_BOND_TYPE SGroup for bond {bond_idx}."
        );
    }
    let Ok(old) = mol.get_bond(bond_idx) else {
        rd_log_warning!(
            "molecule does not contain a bond matching the MRV_COORDINATE_BOND_TYPE SGroup for bond {bond_idx}, ignoring."
        );
        return Ok(());
    };
    if old.bond_type() != BondType::Unspecified {
        rd_log_warning!(
            "MRV_COORDINATE_BOND_TYPE SGroup with value {bond_idx} does not reference a query bond, ignoring."
        );
        return Ok(());
    }
    mol.replace_bond(bond_idx, &Bond::new(BondType::Dative), true, true)
}

/// `processSMARTSQ`: a DAT group holding a SMARTS query for its atoms.
fn process_smartsq(mol: &mut RWMol, sg: &SubstanceGroup) -> Result<()> {
    if let Some(op) = sg.get_prop("QUERYOP").map(RDValue::to_rd_string)
        && op != "="
    {
        rd_log_warning!("unrecognized QUERYOP '{op}' for SMARTSQ. Query ignored.");
        return Ok(());
    }
    let fields = data_fields_of(sg).unwrap_or_default();
    let Some(sma) = fields.first() else {
        rd_log_warning!("empty FIELDDATA for SMARTSQ. Query ignored.");
        return Ok(());
    };
    if fields.len() > 1 {
        rd_log_warning!("multiple FIELDDATA values for SMARTSQ. Taking the first.");
    }
    if sma.is_empty() {
        rd_log_warning!("Skipping empty SMARTS value for SMARTSQ.");
        return Ok(());
    }
    for &aidx in sg.atoms() {
        let parsed = mol_from_smarts(sma, &SmartsParserParams::default())
            .ok()
            .flatten()
            .filter(|m| m.num_atoms() > 0);
        let Some(query_mol) = parsed else {
            rd_log_warning!(
                "SMARTS for SMARTSQ '{sma}' could not be parsed or has no atoms. Ignoring it."
            );
            return Ok(());
        };
        if !mol.get_atom(aidx)?.has_query() {
            let qa = query_atom_from_atom(mol.atom(aidx))?;
            mol.replace_atom(aidx, &qa, false)?;
        }
        let query = if query_mol.num_atoms() == 1 {
            query_mol.atom(0).query().cloned()
        } else {
            Some(make_recursive_structure_query(query_mol, 0))
        };
        let atom = mol.atom_mut(aidx);
        atom.set_query(query);
        atom.set_prop(cp::MRV_SMA, sma.clone(), false)?;
        atom.set_prop(cp::_MolFileAtomQuery, 1i32, false)?;
    }
    Ok(())
}

/// `processMrvImplicitH`: explicit hydrogen counts for aromatic atoms.
fn process_mrv_implicit_h(mol: &mut RWMol, sg: &SubstanceGroup) -> Result<()> {
    for field in data_fields_of(sg).unwrap_or_default() {
        if !field.starts_with("IMPL_H") {
            continue;
        }
        let val = to_int(&field.as_bytes()[6..], true)?;
        for &idx in sg.atoms() {
            if idx >= mol.num_atoms() {
                rd_log_warning!(
                    "bad atom index, {idx}, found in MRV_IMPLICIT_H SGroup. Ignoring it."
                );
            } else if mol
                .atom_bonds(idx)
                .any(|b| b.is_aromatic() || b.bond_type() == BondType::Aromatic)
            {
                mol.atom_mut(idx).set_num_explicit_hs(val as u32);
            } else {
                rd_log_warning!(
                    "MRV_IMPLICIT_H SGroup on atom without aromatic bonds, {idx}, ignored."
                );
            }
        }
    }
    Ok(())
}

/// `processSGroups`: applies (and removes) the groups that carry structure
/// information, namely MarvinSketch's coordinate bonds, implicit Hs and
/// SMARTS queries and RDKit's zero-order bond extensions (ZBO, ZCH, HYD).
pub(super) fn process_sgroups(mol: &mut RWMol) -> Result<()> {
    let mut remove = Vec::new();
    for (idx, sg) in mol.substance_groups().to_vec().iter().enumerate() {
        if type_of(sg) != "DAT" {
            continue;
        }
        let field = sg.get_prop("FIELDNAME").map(RDValue::to_rd_string);
        match field.as_deref() {
            Some("MRV_COORDINATE_BOND_TYPE") => process_mrv_coordinate_bond(mol, sg)?,
            Some("MRV_IMPLICIT_H") => process_mrv_implicit_h(mol, sg)?,
            Some("ZBO") => {
                for &bidx in sg.bonds() {
                    mol.get_bond_mut(bidx)?.set_bond_type(BondType::Zero);
                }
            }
            Some("ZCH") => {
                for values in per_atom_values(sg, "ZCH").unwrap_or_default() {
                    for (&aid, charge) in sg.atoms().iter().zip(values) {
                        mol.get_atom_mut(aid)?.set_formal_charge(charge);
                    }
                }
            }
            Some("HYD") => {
                for values in per_atom_values(sg, "HYD").unwrap_or_default() {
                    for (&aid, hs) in sg.atoms().iter().zip(values) {
                        let atom = mol.get_atom_mut(aid)?;
                        atom.set_prop("_ZBO_H", true, false)?;
                        atom.set_num_explicit_hs(hs as u32);
                    }
                }
            }
            _ => {
                let query_type = sg.get_prop("QUERYTYPE").map(RDValue::to_rd_string);
                if matches!(query_type.as_deref(), Some("SMARTSQ" | "SQ")) {
                    process_smartsq(mol, sg)?;
                } else {
                    continue;
                }
            }
        }
        remove.push(idx);
    }
    for idx in remove.into_iter().rev() {
        mol.substance_groups_mut().remove(idx);
    }
    Ok(())
}
