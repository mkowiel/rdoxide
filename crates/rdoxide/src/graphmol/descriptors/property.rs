//! Property functors and range queries (upstream `GraphMol/Descriptors/Property.{h,cpp}`).
//!
//! A [`PropertyFunctor`] is a named, versioned function from a molecule to a
//! number. [`builtin_properties`] lists the ones upstream registers by
//! default (`registerDescriptors`), in the same order. The mutable registry
//! that users add their own functors to lives in the Python bindings, as those
//! functors are Python objects.

use crate::error::{RDKitError, Result};
use crate::graphmol::descriptors::{connectivity as conn, crippen, lipinski as lip, mol_surf};
use crate::graphmol::mol_props;
use crate::graphmol::romol::ROMol;

/// `double (*d_dataFunc)(const ROMol &)`. It takes the molecule mutably because
/// several descriptors cache their intermediate results on it.
pub type DataFunc = fn(&mut ROMol) -> Result<f64>;

/// `PropertyFunctor`
#[derive(Clone, Debug)]
pub struct PropertyFunctor {
    pub name: String,
    pub version: String,
    pub data_func: Option<DataFunc>,
}

impl PropertyFunctor {
    pub fn new(name: &str, version: &str, data_func: Option<DataFunc>) -> Self {
        PropertyFunctor {
            name: name.to_string(),
            version: version.to_string(),
            data_func,
        }
    }

    /// `operator()`: the value of the property for `mol`.
    pub fn compute(&self, mol: &mut ROMol) -> Result<f64> {
        match self.data_func {
            Some(f) => f(mol),
            None => Err(RDKitError::Value(
                "PropertyFunctor has no data function".into(),
            )),
        }
    }
}

/// `PROP_RANGE_QUERY` made by `makePropertyRangeQuery`: matches molecules whose
/// property lies strictly between `lower` and `upper` (the ends are open).
#[derive(Clone, Debug)]
pub struct PropertyRangeQuery {
    lower: f64,
    upper: f64,
    data_func: Option<DataFunc>,
}

impl PropertyRangeQuery {
    pub fn new(lower: f64, upper: f64, data_func: Option<DataFunc>) -> Self {
        PropertyRangeQuery {
            lower,
            upper,
            data_func,
        }
    }

    /// `Match`
    pub fn matches(&self, mol: &mut ROMol) -> Result<bool> {
        let data_func = self
            .data_func
            .ok_or_else(|| RDKitError::precondition("no data function"))?;
        let value = data_func(mol)?;
        // `queryCmp(bound, value, tol = 0) < 0` for the lower and `> 0` for the upper bound
        Ok(self.lower - value < 0.0 && self.upper - value > 0.0)
    }
}

/// Wraps a descriptor expression as a [`DataFunc`], converting its value to a double.
macro_rules! data_func {
    (|$m:ident| $body:expr) => {
        Some((|$m: &mut ROMol| -> Result<f64> { Ok(f64::from($body?)) }) as DataFunc)
    };
}

fn builtin(name: &str, version: &str, data_func: Option<DataFunc>) -> PropertyFunctor {
    PropertyFunctor::new(name, version, data_func)
}

/// The properties `registerDescriptors` registers, in registration order.
pub fn builtin_properties() -> Vec<PropertyFunctor> {
    vec![
        builtin(
            "exactmw",
            "1.1.0",
            data_func!(|m| mol_props::get_exact_mol_wt(m, false)),
        ),
        builtin(
            "amw",
            "1.0.0",
            data_func!(|m| mol_props::get_avg_mol_wt(m, false)),
        ),
        builtin(
            "lipinskiHBA",
            lip::LIPINSKI_HBA_VERSION,
            data_func!(|m| Ok::<_, RDKitError>(lip::calc_lipinski_hba(m))),
        ),
        builtin(
            "lipinskiHBD",
            lip::LIPINSKI_HBD_VERSION,
            data_func!(|m| lip::calc_lipinski_hbd(m)),
        ),
        builtin(
            "NumRotatableBonds",
            lip::NUM_ROTATABLE_BONDS_VERSION,
            data_func!(|m| lip::calc_num_rotatable_bonds(
                m,
                lip::NumRotatableBondsOptions::Default
            )),
        ),
        builtin(
            "NumHBD",
            lip::NUM_HBD_VERSION,
            data_func!(|m| lip::calc_num_hbd(m)),
        ),
        builtin(
            "NumHBA",
            lip::NUM_HBA_VERSION,
            data_func!(|m| lip::calc_num_hba(m)),
        ),
        builtin(
            "NumHeavyAtoms",
            "1.0.0",
            data_func!(|m| Ok::<_, RDKitError>(m.num_heavy_atoms())),
        ),
        builtin("NumAtoms", "1.0.0", data_func!(|m| m.num_atoms_with_hs())),
        builtin(
            "NumHeteroatoms",
            lip::NUM_HETEROATOMS_VERSION,
            data_func!(|m| lip::calc_num_heteroatoms(m)),
        ),
        builtin(
            "NumAmideBonds",
            lip::NUM_AMIDE_BONDS_VERSION,
            data_func!(|m| lip::calc_num_amide_bonds(m)),
        ),
        builtin(
            "FractionCSP3",
            lip::FRACTION_CSP3_VERSION,
            Some(|m: &mut ROMol| lip::calc_fraction_csp3(m)),
        ),
        builtin(
            "NumRings",
            lip::NUM_RINGS_VERSION,
            data_func!(|m| lip::calc_num_rings(m)),
        ),
        builtin(
            "NumAromaticRings",
            lip::NUM_AROMATIC_RINGS_VERSION,
            data_func!(|m| Ok::<_, RDKitError>(lip::calc_num_aromatic_rings(m))),
        ),
        builtin(
            "NumAliphaticRings",
            lip::NUM_ALIPHATIC_RINGS_VERSION,
            data_func!(|m| Ok::<_, RDKitError>(lip::calc_num_aliphatic_rings(m))),
        ),
        builtin(
            "NumSaturatedRings",
            lip::NUM_SATURATED_RINGS_VERSION,
            data_func!(|m| Ok::<_, RDKitError>(lip::calc_num_saturated_rings(m))),
        ),
        builtin(
            "NumHeterocycles",
            lip::NUM_HETEROCYCLES_VERSION,
            data_func!(|m| Ok::<_, RDKitError>(lip::calc_num_heterocycles(m))),
        ),
        builtin(
            "NumAromaticHeterocycles",
            lip::NUM_AROMATIC_HETEROCYCLES_VERSION,
            data_func!(|m| Ok::<_, RDKitError>(lip::calc_num_aromatic_heterocycles(m))),
        ),
        builtin(
            "NumSaturatedHeterocycles",
            lip::NUM_SATURATED_HETEROCYCLES_VERSION,
            data_func!(|m| Ok::<_, RDKitError>(lip::calc_num_saturated_heterocycles(m))),
        ),
        builtin(
            "NumAliphaticHeterocycles",
            lip::NUM_ALIPHATIC_HETEROCYCLES_VERSION,
            data_func!(|m| Ok::<_, RDKitError>(lip::calc_num_aliphatic_heterocycles(m))),
        ),
        builtin(
            "NumSpiroAtoms",
            lip::NUM_SPIRO_ATOMS_VERSION,
            data_func!(|m| lip::spiro_atoms(m).map(|atoms| atoms.len() as u32)),
        ),
        builtin(
            "NumBridgeheadAtoms",
            lip::NUM_BRIDGEHEAD_ATOMS_VERSION,
            data_func!(|m| lip::bridgehead_atoms(m).map(|atoms| atoms.len() as u32)),
        ),
        builtin(
            "NumAtomStereoCenters",
            lip::NUM_ATOM_STEREO_CENTERS_VERSION,
            data_func!(|m| lip::num_atom_stereo_centers(m)),
        ),
        builtin(
            "NumUnspecifiedAtomStereoCenters",
            lip::NUM_UNSPECIFIED_ATOM_STEREO_CENTERS_VERSION,
            data_func!(|m| lip::num_unspecified_atom_stereo_centers(m)),
        ),
        builtin(
            "labuteASA",
            mol_surf::LABUTE_ASA_VERSION,
            Some(|m: &mut ROMol| mol_surf::calc_labute_asa(m, true, false)),
        ),
        builtin(
            "tpsa",
            mol_surf::TPSA_VERSION,
            Some(|m: &mut ROMol| mol_surf::calc_tpsa(m, false, false)),
        ),
        builtin(
            "CrippenClogP",
            crippen::CRIPPEN_VERSION,
            Some(|m: &mut ROMol| crippen::calc_clogp(m)),
        ),
        builtin(
            "CrippenMR",
            crippen::CRIPPEN_VERSION,
            Some(|m: &mut ROMol| crippen::calc_mr(m)),
        ),
        builtin(
            "chi0v",
            conn::CHI_VERSION,
            Some(|m: &mut ROMol| conn::calc_chi0v(m, false)),
        ),
        builtin(
            "chi1v",
            conn::CHI_VERSION,
            Some(|m: &mut ROMol| conn::calc_chi1v(m, false)),
        ),
        builtin(
            "chi2v",
            conn::CHI_VERSION,
            Some(|m: &mut ROMol| conn::calc_chi_nv(m, 2, false)),
        ),
        builtin(
            "chi3v",
            conn::CHI_VERSION,
            Some(|m: &mut ROMol| conn::calc_chi_nv(m, 3, false)),
        ),
        builtin(
            "chi4v",
            conn::CHI_VERSION,
            Some(|m: &mut ROMol| conn::calc_chi_nv(m, 4, false)),
        ),
        builtin(
            "chi0n",
            conn::CHI_VERSION,
            Some(|m: &mut ROMol| conn::calc_chi0n(m, false)),
        ),
        builtin(
            "chi1n",
            conn::CHI_VERSION,
            Some(|m: &mut ROMol| conn::calc_chi1n(m, false)),
        ),
        builtin(
            "chi2n",
            conn::CHI_VERSION,
            Some(|m: &mut ROMol| conn::calc_chi_nn(m, 2, false)),
        ),
        builtin(
            "chi3n",
            conn::CHI_VERSION,
            Some(|m: &mut ROMol| conn::calc_chi_nn(m, 3, false)),
        ),
        builtin(
            "chi4n",
            conn::CHI_VERSION,
            Some(|m: &mut ROMol| conn::calc_chi_nn(m, 4, false)),
        ),
        builtin(
            "hallKierAlpha",
            conn::HALL_KIER_ALPHA_VERSION,
            Some(|m: &mut ROMol| conn::calc_hall_kier_alpha(m, None)),
        ),
        builtin(
            "kappa1",
            conn::KAPPA_VERSION,
            Some(|m: &mut ROMol| conn::calc_kappa1(m)),
        ),
        builtin(
            "kappa2",
            conn::KAPPA_VERSION,
            Some(|m: &mut ROMol| conn::calc_kappa2(m)),
        ),
        builtin(
            "kappa3",
            conn::KAPPA_VERSION,
            Some(|m: &mut ROMol| conn::calc_kappa3(m)),
        ),
        builtin(
            "Phi",
            conn::PHI_VERSION,
            Some(|m: &mut ROMol| conn::calc_phi(m)),
        ),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

    fn smi(s: &str) -> ROMol {
        mol_from_smiles(s, &SmilesParserParams::default())
            .unwrap()
            .unwrap()
    }

    #[test]
    fn builtin_list_is_in_upstream_order() {
        let props = builtin_properties();
        assert_eq!(props.len(), 43);
        assert_eq!(props[0].name, "exactmw");
        assert_eq!(props[42].name, "Phi");
        assert!(props.iter().all(|p| p.data_func.is_some()));
    }

    #[test]
    fn functors_compute_and_report_missing_functions() {
        let mut m = smi("CC(=O)Nc1ccc(O)cc1");
        let props = builtin_properties();
        let by_name = |n: &str| props.iter().find(|p| p.name == n).unwrap();
        assert!((by_name("exactmw").compute(&mut m).unwrap() - 151.063328528).abs() < 1e-9);
        assert_eq!(by_name("NumHeavyAtoms").compute(&mut m).unwrap(), 11.0);
        let empty = PropertyFunctor::new("x", "1", None);
        assert!(empty.compute(&mut m).is_err());
    }

    #[test]
    fn range_queries_are_open_intervals() {
        let mut m = smi("CCO");
        let exact_mw = builtin_properties()[0].data_func;
        let value = 46.041864812;
        assert!(
            PropertyRangeQuery::new(0.0, 100.0, exact_mw)
                .matches(&mut m)
                .unwrap()
        );
        assert!(
            !PropertyRangeQuery::new(0.0, value, exact_mw)
                .matches(&mut m)
                .unwrap()
        );
        assert!(
            PropertyRangeQuery::new(0.0, 100.0, None)
                .matches(&mut m)
                .is_err()
        );
    }
}
