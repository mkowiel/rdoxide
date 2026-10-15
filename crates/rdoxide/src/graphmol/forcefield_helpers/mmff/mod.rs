//! MMFF atom typing and force-field construction (upstream
//! `ForceFieldHelpers/MMFF`).

mod charges;
mod empirical;
mod heavy;
mod props;
mod rings;

pub use empirical::angle_bend_empirical_rule_params;
pub use props::{
    DielectricModel, MmffAtomProperties, MmffMolProperties, sanitize_mmff_mol, stretch_bend_type,
};

pub use rings::{
    are_atoms_in_same_aromatic_ring, are_atoms_in_same_ring_of_size,
    is_angle_in_ring_of_size_3_or_4, is_atom_in_aromatic_ring_of_size, is_atom_n_oxide,
    is_ring_aromatic, is_torsion_in_ring_of_size_4_or_5, periodic_table_row, periodic_table_row_hl,
};
