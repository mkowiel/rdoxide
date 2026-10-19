//! Atom and bond property lists: the individual atom or bond values of one
//! property written as a single molecule property (`atom.iprop.<name>` and
//! its relatives), as `FileParserUtils::createAtomIntPropertyList` and the
//! other `create*PropertyList` functions do.
//!
//! Ported from RDKit 2026.09.1 `Code/GraphMol/FileParsers/FileParserUtils.h`.
//! Copyright (C) 2018-2019 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::Result;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::dict::{RDProps, RDValue};

/// `FileParserUtils::DEFAULT_LINESIZE`
pub const DEFAULT_LINESIZE: usize = 190;

/// The value type of a property list.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ListKind {
    Int,
    Double,
    Bool,
    String,
}

impl ListKind {
    /// The marker in the name of the list property (`iprop`, `dprop`, ...).
    pub fn type_marker(self) -> &'static str {
        match self {
            ListKind::Int => "iprop",
            ListKind::Double => "dprop",
            ListKind::Bool => "bprop",
            ListKind::String => "prop",
        }
    }

    /// `getPropIfPresent<T>` followed by `lexical_cast<std::string>`: the text of
    /// the value if it is, or converts to, this type.
    fn text(self, value: &RDValue) -> Option<String> {
        match self {
            ListKind::Int => match value {
                RDValue::String(s) => s.trim().parse::<i32>().ok().map(|v| v.to_string()),
                v => v.as_int().map(|v| v.to_string()),
            },
            ListKind::Double => match value {
                RDValue::String(s) => s
                    .trim()
                    .parse::<f64>()
                    .ok()
                    .map(|v| RDValue::Double(v).to_rd_string()),
                v => v.as_double().map(|v| RDValue::Double(v).to_rd_string()),
            },
            ListKind::Bool => match value {
                RDValue::String(s) => match s.trim() {
                    "1" | "true" => Some("1".to_string()),
                    "0" | "false" => Some("0".to_string()),
                    _ => None,
                },
                RDValue::Bool(b) => Some(if *b { "1" } else { "0" }.to_string()),
                v => v
                    .as_int()
                    .map(|v| if v != 0 { "1" } else { "0" }.to_string()),
            },
            ListKind::String => Some(value.to_rd_string()),
        }
    }
}

/// `getPropertyList`: the values, space separated and wrapped before
/// `line_size`; a missing value shows as `missing_marker` (`n/a` by default,
/// in which case the list does not announce the marker).
pub fn property_list<'a>(
    values: impl IntoIterator<Item = Option<&'a RDValue>>,
    kind: ListKind,
    missing_value_marker: &str,
    line_size: usize,
) -> String {
    let mut res = String::new();
    let mut prop_val = String::new();
    let marker = if missing_value_marker.is_empty() {
        "n/a"
    } else {
        prop_val += &format!("[{missing_value_marker}] ");
        missing_value_marker
    };
    for value in values {
        let item = value
            .and_then(|v| kind.text(v))
            .unwrap_or_else(|| marker.to_string());
        if prop_val.len() + item.len() + 1 >= line_size {
            // remove the trailing space
            prop_val.pop();
            res += &prop_val;
            res.push('\n');
            prop_val.clear();
        }
        prop_val += &item;
        prop_val.push(' ');
    }
    if !prop_val.is_empty() {
        prop_val.pop();
        res += &prop_val;
    }
    res
}

/// `createAtom<Type>PropertyList`: sets `atom.<marker>.<prop_name>` on the molecule.
pub fn create_atom_property_list(
    mol: &mut ROMol,
    kind: ListKind,
    prop_name: &str,
    missing_value_marker: &str,
    line_size: usize,
) -> Result<()> {
    let list = property_list(
        mol.atoms().iter().map(|a| a.get_prop(prop_name)),
        kind,
        missing_value_marker,
        line_size,
    );
    mol.set_prop(
        &format!("atom.{}.{prop_name}", kind.type_marker()),
        list,
        false,
    )
}

/// `createBond<Type>PropertyList`: sets `bond.<marker>.<prop_name>` on the molecule.
pub fn create_bond_property_list(
    mol: &mut ROMol,
    kind: ListKind,
    prop_name: &str,
    missing_value_marker: &str,
    line_size: usize,
) -> Result<()> {
    let list = property_list(
        mol.bonds().iter().map(|b| b.get_prop(prop_name)),
        kind,
        missing_value_marker,
        line_size,
    );
    mol.set_prop(
        &format!("bond.{}.{prop_name}", kind.type_marker()),
        list,
        false,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_values_and_wrapping() {
        let (a, b) = (RDValue::Int(1), RDValue::Int(3));
        let vals = [Some(&a), None, Some(&b)];
        assert_eq!(property_list(vals, ListKind::Int, "", 190), "1 n/a 3");
        assert_eq!(property_list(vals, ListKind::Int, "?", 190), "[?] 1 ? 3");
        // each line stays below the line size
        assert_eq!(property_list(vals, ListKind::Int, "", 7), "1 n/a\n3");
    }

    #[test]
    fn values_of_other_types_convert() {
        let s = RDValue::String("7".into());
        let d = RDValue::Double(0.5);
        assert_eq!(
            property_list([Some(&s), Some(&d)], ListKind::Int, "", 190),
            "7 n/a"
        );
        assert_eq!(
            property_list([Some(&s), Some(&d)], ListKind::Double, "", 190),
            "7 0.5"
        );
    }
}
