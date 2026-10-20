//! The input atom of InChI (`INCHI_BASE/src/inpdef.h`, `inp_ATOM`) and its
//! constants. Arrays have the fixed sizes of the C structure so that the
//! translated algorithms index them the same way.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::incomdef::NUM_H_ISOTOPES;

/// `MAXVAL`: maximum number of bonds per atom
pub const MAXVAL: usize = 20;
/// `ATOM_EL_LEN`: length of the element symbol field, with the terminating zero
pub const ATOM_EL_LEN: usize = 6;
/// `MAX_NUM_STEREO_ATOM_NEIGH`
pub const MAX_NUM_STEREO_ATOM_NEIGH: usize = 4;
/// `MAX_NUM_STEREO_BONDS`
pub const MAX_NUM_STEREO_BONDS: usize = 3;
/// `MAX_ATOMS` (`ichisize.h`)
pub const MAX_ATOMS: usize = 32766;

// `inp_ATOM::at_type`
pub const ATT_NONE: u16 = 0x0000;
pub const ATT_ACIDIC_CO: u16 = 0x0001;
pub const ATT_ACIDIC_S: u16 = 0x0002;
pub const ATT_OO: u16 = 0x0004;
pub const ATT_ZOO: u16 = 0x0008;
pub const ATT_NO: u16 = 0x0010;
pub const ATT_N_O: u16 = 0x0020;
pub const ATT_ATOM_N: u16 = 0x0040;
pub const ATT_ATOM_P: u16 = 0x0080;
pub const ATT_OTHER_NEG_O: u16 = 0x0100;
/// `-Z=O` or `=Z=O`
pub const ATT_OTHER_ZO: u16 = 0x0200;
/// `OH(-)`, `O=O,S,Se,Te`
pub const ATT_OH_MINUS: u16 = 0x0400;
/// `-OH2(+)`, `=OH(+)`, `-OH(+)-`, `OH3(+)`, `=O(+)-`, etc.
pub const ATT_O_PLUS: u16 = 0x0800;
pub const ATT_PROTON: u16 = 0x1000;
pub const ATT_HAL_ANION: u16 = 0x2000;
pub const ATT_HAL_ACID: u16 = 0x4000;
/// `=N(-)` or `=P(-)` where `=` previously was triple (`FIX_NP_MINUS_BUG`)
pub const ATT_NP_MINUS_V23: u16 = 0x8000;

/// `AT_FLAG_ISO_H_POINT`: may have isotopic H
pub const AT_FLAG_ISO_H_POINT: i8 = 0x01;
pub const PERIODIC_NUMBER_H: u8 = 1;

pub const FLAG_SC_0D: i8 = 1;
pub const FLAG_SB_0D: i8 = 2;

/// mask for the disconnected metal parity if it is different
pub const SB_PARITY_FLAG: i8 = 0x38;
pub const SB_PARITY_SHFT: i32 = 3;
pub const SB_PARITY_MASK: i32 = 0x07;

/// `SB_PARITY_1`: refers to the connected structure
pub fn sb_parity_1(x: i8) -> i32 {
    i32::from(x) & SB_PARITY_MASK
}

/// `SB_PARITY_2`
pub fn sb_parity_2(x: i8) -> i32 {
    (i32::from(x) >> SB_PARITY_SHFT) & SB_PARITY_MASK
}

/// `inp_ATOM` (`tagInputAtom`)
#[derive(Debug, Clone, PartialEq)]
pub struct InpAtom {
    /// name of the chemical element (zero terminated)
    pub elname: [u8; ATOM_EL_LEN],
    /// number of the element in the Periodic Table
    pub el_number: u8,
    /// positions (from 0) of the neighbors in the atom array
    pub neighbor: [u16; MAXVAL],
    /// original atom number
    pub orig_at_number: u16,
    /// original atom number within the component, before terminal H removal
    pub orig_compt_at_numb: u16,
    /// 1=Up, 4=Either, 6=Down; this atom is at the pointing wedge, negative
    /// => on the opposite side; 3=Either double bond
    pub bond_stereo: [i8; MAXVAL],
    /// 1..4; 4="aromatic", should be discouraged on input
    pub bond_type: [u8; MAXVAL],
    /// coordination number: number of bonds = number of neighbors
    pub valence: i8,
    /// sum of bond types (type 4 needs special treatment)
    pub chem_bonds_valence: i8,
    /// number of implicit hydrogens, including D and T
    pub num_h: i8,
    /// number of implicit 1H, 2H(D), 3H(T) < 16
    pub num_iso_h: [i8; NUM_H_ISOTOPES],
    /// 0 => natural isotopic abundances; >0 => (mass) - (mass of the most
    /// abundant isotope) + 1; <0 => (mass) - (mass of the most abundant isotope)
    pub iso_atw_diff: i8,
    pub charge: i8,
    /// `RADICAL_SINGLET`, `RADICAL_DOUBLET`, or `RADICAL_TRIPLET`
    pub radical: i8,
    pub b_ambiguous_stereo: i8,
    /// `AT_FLAG_ISO_H_POINT`
    pub c_flags: i8,
    /// `ATT_NONE`, `ATT_ACIDIC_CO`, ...
    pub at_type: u16,
    /// number of the structure component > 0
    pub component: u16,
    /// id of a tautomeric group
    pub endpoint: u16,
    /// id of a positive charge group
    pub c_point: u16,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// bit=1 => stereobond; bit=2 => stereocenter
    pub b_used_0d_parity: i8,
    /// 0D tetrahedral parity
    pub p_parity: i8,
    pub p_orig_at_num: [u16; MAX_NUM_STEREO_ATOM_NEIGH],
    /// stereo bond/neighbor ordering number, starts from 0
    pub sb_ord: [i8; MAX_NUM_STEREO_BONDS],
    /// ordering number of the neighbor adjacent to the stereo bond, starts
    /// from 0; -1 means a removed explicit H
    pub sn_ord: [i8; MAX_NUM_STEREO_BONDS],
    pub sb_parity: [i8; MAX_NUM_STEREO_BONDS],
    /// original atom number of the `sn_ord` neighbors
    pub sn_orig_at_num: [u16; MAX_NUM_STEREO_BONDS],
    // `FIND_RING_SYSTEMS == 1`
    pub b_cut_vertex: i8,
    pub n_ring_system: u16,
    pub n_num_at_in_ring_system: u16,
    pub n_block_system: u16,
}

impl Default for InpAtom {
    fn default() -> Self {
        InpAtom {
            elname: [0; ATOM_EL_LEN],
            el_number: 0,
            neighbor: [0; MAXVAL],
            orig_at_number: 0,
            orig_compt_at_numb: 0,
            bond_stereo: [0; MAXVAL],
            bond_type: [0; MAXVAL],
            valence: 0,
            chem_bonds_valence: 0,
            num_h: 0,
            num_iso_h: [0; NUM_H_ISOTOPES],
            iso_atw_diff: 0,
            charge: 0,
            radical: 0,
            b_ambiguous_stereo: 0,
            c_flags: 0,
            at_type: 0,
            component: 0,
            endpoint: 0,
            c_point: 0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            b_used_0d_parity: 0,
            p_parity: 0,
            p_orig_at_num: [0; MAX_NUM_STEREO_ATOM_NEIGH],
            sb_ord: [0; MAX_NUM_STEREO_BONDS],
            sn_ord: [0; MAX_NUM_STEREO_BONDS],
            sb_parity: [0; MAX_NUM_STEREO_BONDS],
            sn_orig_at_num: [0; MAX_NUM_STEREO_BONDS],
            b_cut_vertex: 0,
            n_ring_system: 0,
            n_num_at_in_ring_system: 0,
            n_block_system: 0,
        }
    }
}

impl InpAtom {
    /// The element name as a string (up to the terminating zero).
    pub fn elname_str(&self) -> &str {
        let n = self
            .elname
            .iter()
            .position(|&c| c == 0)
            .unwrap_or(ATOM_EL_LEN);
        std::str::from_utf8(&self.elname[..n]).unwrap_or("")
    }

    /// Sets the element name (truncated to fit with the terminating zero).
    pub fn set_elname(&mut self, name: &str) {
        self.elname = [0; ATOM_EL_LEN];
        for (d, s) in self
            .elname
            .iter_mut()
            .zip(name.bytes().take(ATOM_EL_LEN - 1))
        {
            *d = s;
        }
    }

    /// `NUM_ISO_H`
    pub fn num_iso_h_total(&self) -> i32 {
        self.num_iso_h.iter().map(|v| i32::from(*v)).sum()
    }

    /// `NUMH`: the implicit hydrogens, isotopic ones included.
    pub fn numh(&self) -> i32 {
        i32::from(self.num_h) + self.num_iso_h_total()
    }
}

/// `BOND_TYPE_MASK` (`extr_ct.h`)
pub const BOND_TYPE_MASK: i32 = 0x0f;
