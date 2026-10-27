//! `reaccs.h`, `mdl.h`, `symbol_lists.h`: the MOL file data model of the Avalon Toolkit.

/// `MAXATOMS`, `MAXBONDS`
pub const MAXATOMS: i32 = 999;
pub const MAXBONDS: i32 = 999;

/// `NONE`
pub const NONE: i32 = 0;
/// charge/radical codes of `mdl.h`
pub const RADICAL: i32 = 4;
pub const ANY_CHARGE: i32 = 8;
pub const DOUBLET: i32 = 2;
/// bond stereo codes
pub const CIS_TRANS_EITHER: i32 = 0x03;
pub const CIS_TRANS_SWAPPED: i32 = 0x08;
/// `symbol_lists.h`
pub const EXCLUSIVE: i32 = 0;
pub const INCLUSIVE: i32 = 1;

/// `struct reaccs_atom_t`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Atom {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub atom_symbol: String,
    pub mass_difference: i32,
    pub charge: i32,
    pub radical: i32,
    pub stereo_parity: i32,
    pub query_h_count: i32,
    pub query_stereo_box: i32,
    pub dummy1: i32,
    pub dummy2: i32,
    pub dummy3: i32,
    pub dummy4: i32,
    pub mapping: i32,
    pub second_stereo_parity: i32,
    pub dummy5: i32,
    pub sub_desc: i32,
    pub value: f32,
    pub rsize_flags: u32,
    pub color: i32,
    pub atext: String,
}

/// `struct reaccs_bond_t`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Bond {
    pub atoms: [i32; 2],
    pub bond_type: i32,
    pub stereo_symbol: i32,
    pub dummy: i32,
    pub rsize_flags: u32,
    pub topography: i32,
    pub reaction_mark: i32,
    pub value: f32,
    pub bond_type_flags: i32,
    pub color: i32,
}

/// `struct stext_line_t`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct STextLine {
    pub x: f32,
    pub y: f32,
    pub text: String,
}

/// `struct symbol_list_t`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SymbolList {
    /// 1-based atom number
    pub atom: i32,
    pub string: String,
    pub logic: i32,
}

/// `struct reaccs_molecule_t` (the linked lists of C are vectors in list order)
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Molecule {
    pub name: String,
    pub user_initials: String,
    pub program_name: String,
    pub date: String,
    pub time: String,
    pub dimensionality: String,
    pub scale1: i32,
    pub scale2: f32,
    pub energy: f32,
    pub registry_number: i64,
    pub comment: String,
    pub n_atom_lists: i32,
    pub dummy1: i32,
    pub chiral_flag: i32,
    pub stext_lines: Vec<STextLine>,
    pub n_props: i32,
    pub prop_lines: Vec<String>,
    pub version: String,
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub symbol_lists: Vec<SymbolList>,
    pub color: i32,
}

impl Molecule {
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn n_bonds(&self) -> usize {
        self.bonds.len()
    }
}

/// `CombineChargeRadical`: MDL's combined charge/radical code
pub fn combine_charge_radical(charge: i32, radical: i32) -> i32 {
    if radical == NONE && (-3..=3).contains(&charge) && charge != 0 {
        RADICAL - charge
    } else if charge == 0 && radical == NONE {
        NONE
    } else if charge == 0 && radical == DOUBLET {
        RADICAL
    } else {
        ANY_CHARGE
    }
}

/// `SplitChargeRadical`: `(charge, radical)`
pub fn split_charge_radical(charge_radical: i32) -> (i32, i32) {
    if charge_radical == NONE {
        (NONE, NONE)
    } else if charge_radical == RADICAL {
        (NONE, DOUBLET)
    } else {
        (RADICAL - charge_radical, NONE)
    }
}

/// `IntToString(periodic_table, id)`: the first symbol with that atomic number
pub fn symbol_of(id: i32) -> Option<&'static str> {
    super::tables::PERIODIC_TABLE
        .iter()
        .find(|(_, n)| *n == id)
        .map(|(s, _)| *s)
}

/// `StringToInt(periodic_table, s)`: the atomic number, `0` (the id of the sentinel) if unknown
pub fn number_of(symbol: &str) -> i32 {
    super::tables::PERIODIC_TABLE
        .iter()
        .find(|(s, _)| *s == symbol)
        .map_or(0, |(_, n)| *n)
}

/// `LookupMassDifference`: the MOL file mass difference for an isotopic mass
pub fn lookup_mass_difference(mass: f64, symbol: &str) -> i32 {
    for (s, _, m, _, _) in super::tables::PTABLE {
        if *s == symbol {
            let diff = mass - f64::from(*m);
            return if diff > 0.0 {
                (diff + 0.49) as i32
            } else {
                (diff - 0.49) as i32
            };
        }
    }
    if mass > 0.0 {
        (mass + 0.49) as i32
    } else {
        (mass - 0.49) as i32
    }
}
