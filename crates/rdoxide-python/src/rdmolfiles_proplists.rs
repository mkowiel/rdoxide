//! `Create{Atom,Bond}{Int,Double,Bool,String}PropertyList` of `rdoxide.Chem.rdmolfiles`
//! (upstream `FileParsers/FileParserUtils.h`).

use pyo3::prelude::*;
use rdoxide::graphmol::file_parsers::property_lists::{
    DEFAULT_LINESIZE, ListKind, create_atom_property_list, create_bond_property_list,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

macro_rules! list_function {
    ($rust:ident, $name:literal, $create:ident, $kind:expr, $what:literal) => {
        #[doc = concat!("creates a list property on the molecule from individual ", $what, " property values")]
        #[pyfunction]
        #[pyo3(name = $name, signature = (mol, propName, missingValueMarker = "", lineSize = DEFAULT_LINESIZE))]
        #[allow(non_snake_case)]
        fn $rust(
            mut mol: PyRefMut<'_, PyMol>,
            propName: &str,
            missingValueMarker: &str,
            lineSize: usize,
        ) -> PyResult<()> {
            $create(&mut mol.mol, $kind, propName, missingValueMarker, lineSize).py()
        }
    };
}

list_function!(
    atom_int,
    "CreateAtomIntPropertyList",
    create_atom_property_list,
    ListKind::Int,
    "atom"
);
list_function!(
    atom_double,
    "CreateAtomDoublePropertyList",
    create_atom_property_list,
    ListKind::Double,
    "atom"
);
list_function!(
    atom_bool,
    "CreateAtomBoolPropertyList",
    create_atom_property_list,
    ListKind::Bool,
    "atom"
);
list_function!(
    atom_string,
    "CreateAtomStringPropertyList",
    create_atom_property_list,
    ListKind::String,
    "atom"
);
list_function!(
    bond_int,
    "CreateBondIntPropertyList",
    create_bond_property_list,
    ListKind::Int,
    "bond"
);
list_function!(
    bond_double,
    "CreateBondDoublePropertyList",
    create_bond_property_list,
    ListKind::Double,
    "bond"
);
list_function!(
    bond_bool,
    "CreateBondBoolPropertyList",
    create_bond_property_list,
    ListKind::Bool,
    "bond"
);
list_function!(
    bond_string,
    "CreateBondStringPropertyList",
    create_bond_property_list,
    ListKind::String,
    "bond"
);

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(atom_int, m)?)?;
    m.add_function(wrap_pyfunction!(atom_double, m)?)?;
    m.add_function(wrap_pyfunction!(atom_bool, m)?)?;
    m.add_function(wrap_pyfunction!(atom_string, m)?)?;
    m.add_function(wrap_pyfunction!(bond_int, m)?)?;
    m.add_function(wrap_pyfunction!(bond_double, m)?)?;
    m.add_function(wrap_pyfunction!(bond_bool, m)?)?;
    m.add_function(wrap_pyfunction!(bond_string, m)?)?;
    Ok(())
}
