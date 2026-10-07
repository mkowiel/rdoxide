//! Python enum wrappers for the core rdchem enums.

use rdoxide::graphmol::stereo_group::StereoGroupType;
use rdoxide::graphmol::{
    BondDir, BondStereo, BondType, ChiralType, HybridizationType, ValenceType,
};

use crate::enums::py_enum;

py_enum!(PyBondType, BondType, "BondType");
py_enum!(PyBondDir, BondDir, "BondDir");
py_enum!(PyBondStereo, BondStereo, "BondStereo");
py_enum!(PyChiralType, ChiralType, "ChiralType");
py_enum!(PyHybridizationType, HybridizationType, "HybridizationType");
py_enum!(PyValenceType, ValenceType, "ValenceType");
py_enum!(PyStereoGroupType, StereoGroupType, "StereoGroupType");
