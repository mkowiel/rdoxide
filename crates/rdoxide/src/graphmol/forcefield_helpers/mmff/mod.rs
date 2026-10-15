//! MMFF atom typing and force-field construction (upstream
//! `ForceFieldHelpers/MMFF`).

mod builder;
mod charges;
mod empirical;
mod heavy;
mod param_getters;
mod props;
mod rings;

pub use empirical::angle_bend_empirical_rule_params;
pub use param_getters::{AngleBendParams, BondStretchParams, StretchBendParams, TorsionParams};
pub use props::{
    DielectricModel, MmffAtomProperties, MmffMolProperties, sanitize_mmff_mol, stretch_bend_type,
};

pub use builder::{
    DEFAULT_NONBONDED_THRESH, DEFAULT_TORSION_BOND_SMARTS, add_angles, add_bonds, add_ele,
    add_nonbonded, add_oop, add_stretch_bend, add_torsions, add_vdw, build_neighbor_matrix,
    construct_force_field, construct_force_field_with_props,
};
pub use rings::{
    are_atoms_in_same_aromatic_ring, are_atoms_in_same_ring_of_size,
    is_angle_in_ring_of_size_3_or_4, is_atom_in_aromatic_ring_of_size, is_atom_n_oxide,
    is_ring_aromatic, is_torsion_in_ring_of_size_4_or_5, periodic_table_row, periodic_table_row_hl,
};
