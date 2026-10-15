//! Experimental torsion-angle preferences and the other knowledge terms of
//! ETKDG (upstream `GraphMol/ForceFieldHelpers/CrystalFF`).

mod params;
mod torsion_preferences;

pub use params::{
    ExpTorsionAngle, GaussianExpTorsionAngle, cosine_torsion_params, gaussian_torsion_params,
};
pub use torsion_preferences::{
    AMIDE_CIS_IDX, AMIDE_TRANS_IDX, AROM_TORSION_IDX, CrystalFFDetails, EtkdgForceConsts,
    ImproperAtoms, TorsionAngleRef, TorsionBond, TorsionOptions, TorsionParamKind, TorsionParams,
    get_experimental_torsions, get_experimental_torsions_with_bonds, populate_ref_table,
};
