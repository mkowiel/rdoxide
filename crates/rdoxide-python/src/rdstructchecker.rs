//! `rdoxide.Chem.rdStructChecker`: checking and fixing the structure of
//! molecules (upstream `GraphMol/StructChecker/Wrap/structchecker.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::struct_checker::structure_flags as sf;
use rdoxide::graphmol::struct_checker::{StructChecker, StructCheckerOptions};

use crate::enums::add_enum;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

const STRUCTURE_FLAGS: &[(&str, i64)] = &[
    ("NO_CHANGE", sf::NO_CHANGE as i64),
    ("BAD_MOLECULE", sf::BAD_MOLECULE as i64),
    (
        "ALIAS_CONVERSION_FAILED",
        sf::ALIAS_CONVERSION_FAILED as i64,
    ),
    ("STEREO_ERROR", sf::STEREO_ERROR as i64),
    ("STEREO_FORCED_BAD", sf::STEREO_FORCED_BAD as i64),
    ("ATOM_CLASH", sf::ATOM_CLASH as i64),
    ("ATOM_CHECK_FAILED", sf::ATOM_CHECK_FAILED as i64),
    ("SIZE_CHECK_FAILED", sf::SIZE_CHECK_FAILED as i64),
    ("TRANSFORMED", sf::TRANSFORMED as i64),
    ("FRAGMENTS_FOUND", sf::FRAGMENTS_FOUND as i64),
    ("EITHER_WARNING", sf::EITHER_WARNING as i64),
    ("DUBIOUS_STEREO_REMOVED", sf::DUBIOUS_STEREO_REMOVED as i64),
    ("RECHARGED", sf::RECHARGED as i64),
    ("STEREO_TRANSFORMED", sf::STEREO_TRANSFORMED as i64),
    ("TEMPLATE_TRANSFORMED", sf::TEMPLATE_TRANSFORMED as i64),
    ("TAUTOMER_TRANSFORMED", sf::TAUTOMER_TRANSFORMED as i64),
];

/// The options of the structure checker.
#[pyclass(
    name = "StructCheckerOptions",
    module = "rdoxide.Chem.rdStructChecker",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyStructCheckerOptions {
    opts: StructCheckerOptions,
}

macro_rules! option_fields {
    ($($getter:ident, $setter:ident, $name:literal, $field:ident: $ty:ty;)*) => {
        #[pymethods]
        impl PyStructCheckerOptions {
            #[new]
            fn py_new() -> Self {
                Self::default()
            }

            $(
                #[getter($name)]
                fn $getter(&self) -> $ty {
                    self.opts.$field
                }
                #[setter($name)]
                fn $setter(&mut self, v: $ty) {
                    self.opts.$field = v;
                }
            )*

            /// Load the set of good augmented atoms from the specified file path
            #[pyo3(name = "LoadGoodAugmentedAtoms", signature = (path))]
            fn load_good_augmented_atoms(&mut self, path: &str) -> bool {
                self.opts.load_good_augmented_atoms(path)
            }

            /// Load the set of acidic augmented atoms from the specified file path
            #[pyo3(name = "LoadAcidicAugmentedAtoms", signature = (path))]
            fn load_acidic_augmented_atoms(&mut self, path: &str) -> bool {
                self.opts.load_acidic_augmented_atoms(path)
            }

            /// Load the set of acidic augmented atoms from the specified file path
            #[pyo3(name = "LoadAugmentedAtomTranslations", signature = (path))]
            fn load_augmented_atom_translations(&mut self, path: &str) -> PyResult<bool> {
                self.opts.load_augmented_atom_translations(path).py()
            }
        }
    };
}

option_fields! {
    acidity_limit, set_acidity_limit, "AcidityLimit", acidity_limit: f64;
    remove_minor_fragments, set_remove_minor_fragments, "RemoveMinorFragments", remove_minor_fragments: bool;
    desired_charge, set_desired_charge, "DesiredCharge", desired_charge: i32;
    check_collisions, set_check_collisions, "CheckCollisions", check_collisions: bool;
    collision_limit_percent, set_collision_limit_percent, "CollisionLimitPercent", collision_limit_percent: i32;
    max_mol_size, set_max_mol_size, "MaxMolSize", max_mol_size: u32;
    convert_s_text, set_convert_s_text, "ConvertSText", convert_s_text: bool;
    strip_zeros, set_strip_zeros, "StripZeros", strip_zeros: bool;
    check_stereo, set_check_stereo, "CheckStereo", check_stereo: bool;
    convert_atom_texts, set_convert_atom_texts, "ConvertAtomTexts", convert_atom_texts: bool;
    groups_to_s_groups, set_groups_to_s_groups, "GroupsToSGroups", groups_to_s_groups: bool;
    verbose, set_verbose, "Verbose", verbose: bool;
}

/// Checks and fixes the structure of molecules.
#[pyclass(name = "StructChecker", module = "rdoxide.Chem.rdStructChecker")]
pub struct PyStructChecker {
    checker: StructChecker,
}

#[pymethods]
impl PyStructChecker {
    #[new]
    #[pyo3(signature = (options = None))]
    fn py_new(options: Option<PyStructCheckerOptions>) -> Self {
        PyStructChecker {
            checker: match options {
                Some(o) => StructChecker::with_options(o.opts),
                None => StructChecker::new(),
            },
        }
    }

    /// Check the structure and return a set of structure flags
    #[pyo3(name = "CheckMolStructure", signature = (mol))]
    fn check_mol_structure(&self, py: Python<'_>, mut mol: PyRefMut<'_, PyMol>) -> PyResult<u32> {
        let m = &mut mol.mol;
        py.detach(|| self.checker.check_mol_structure(m)).py()
    }

    /// Return the structure flags as a human readable string
    #[staticmethod]
    #[pyo3(name = "StructureFlagsToString", signature = (flags))]
    fn structure_flags_to_string(flags: u32) -> String {
        StructChecker::structure_flags_to_string(flags)
    }

    /// Convert a comma separated string to the appropriate structure flags
    #[staticmethod]
    #[pyo3(name = "StringToStructureFlags", signature = (str))]
    fn string_to_structure_flags(str: &str) -> u32 {
        StructChecker::string_to_structure_flags(str)
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add_enum(m, "StructureFlags", STRUCTURE_FLAGS, None)?;
    m.add_class::<PyStructCheckerOptions>()?;
    m.add_class::<PyStructChecker>()?;
    Ok(())
}
