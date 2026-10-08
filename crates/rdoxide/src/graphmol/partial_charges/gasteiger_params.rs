//! Gasteiger partial charge parameters (upstream `GasteigerParams.{h,cpp}`).
//!
//! Upstream parses these tables from strings at startup; here they are
//! static arrays in the same order. Each row is `(element, mode, [a, b, c])`
//! for the orbital electronegativity `a + b*q + c*q^2`.

use crate::error::{RDKitError, Result};

/// Electronegativity of a hydrogen cation, used for the hydrogen atoms in the
/// charge equalization.
pub const IONXH: f64 = 20.02;
/// Factor the damping is scaled by after each iteration.
pub const DAMP_SCALE: f64 = 0.5;
/// Initial damping of the charge update.
pub const DAMP: f64 = 0.5;

type Row = (&'static str, &'static str, [f64; 3]);

/// `defaultParamData`
const DEFAULT_PARAMS: &[Row] = &[
    ("H", "*", [7.17, 6.24, -0.56]),
    ("C", "sp3", [7.98, 9.18, 1.88]),
    ("C", "sp2", [8.79, 9.32, 1.51]),
    ("C", "sp", [10.39, 9.45, 0.73]),
    ("N", "sp3", [11.54, 10.82, 1.36]),
    ("N", "sp2", [12.87, 11.15, 0.85]),
    ("N", "sp", [15.68, 11.7, -0.27]),
    ("O", "sp3", [14.18, 12.92, 1.39]),
    ("O", "sp2", [17.07, 13.79, 0.47]),
    ("F", "sp3", [14.66, 13.85, 2.31]),
    ("Cl", "sp3", [11.00, 9.69, 1.35]),
    ("Br", "sp3", [10.08, 8.47, 1.16]),
    ("I", "sp3", [9.9, 7.96, 0.96]),
    ("S", "sp3", [10.14, 9.13, 1.38]),
    ("S", "so", [10.14, 9.13, 1.38]),
    ("S", "so2", [12.00, 10.81, 1.20]),
    ("S", "sp2", [10.88, 9.49, 1.33]),
    ("P", "sp3", [8.90, 8.24, 0.96]),
    ("X", "*", [0.00, 0.00, 0.00]),
];

/// `additionalParamData` (generated with PyBabel)
const ADDITIONAL_PARAMS: &[Row] = &[
    ("P", "sp2", [9.665, 8.530, 0.735]),
    ("Si", "sp3", [7.300, 6.567, 0.657]),
    ("Si", "sp2", [7.905, 6.748, 0.443]),
    ("Si", "sp", [9.065, 7.027, -0.002]),
    ("B", "sp3", [5.980, 6.820, 1.605]),
    ("B", "sp2", [6.420, 6.807, 1.322]),
    ("Be", "sp3", [3.845, 6.755, 3.165]),
    ("Be", "sp2", [4.005, 6.725, 3.035]),
    ("Mg", "sp2", [3.565, 5.572, 2.197]),
    ("Mg", "sp3", [3.300, 5.587, 2.447]),
    ("Mg", "sp", [4.040, 5.472, 1.823]),
    ("Al", "sp3", [5.375, 4.953, 0.867]),
    ("Al", "sp2", [5.795, 5.020, 0.695]),
];

fn lookup(elem: &str, mode: &str) -> Option<[f64; 3]> {
    // later rows win, as when upstream fills its map in order
    DEFAULT_PARAMS
        .iter()
        .chain(ADDITIONAL_PARAMS)
        .rev()
        .find(|(e, m, _)| *e == elem && *m == mode)
        .map(|(_, _, p)| *p)
}

/// `GasteigerParams::getParams(elem, mode, throwOnFailure)`: the parameters of
/// an element in a bonding mode. Unknown combinations raise `ValueError` when
/// `throw_on_failure` is set and otherwise get the all-zero `X *` parameters,
/// which removes the atom from the iteration.
pub fn get_params(elem: &str, mode: &str, throw_on_failure: bool) -> Result<[f64; 3]> {
    if let Some(p) = lookup(elem, mode) {
        return Ok(p);
    }
    if throw_on_failure {
        return Err(RDKitError::Value(format!(
            "ERROR: No Gasteiger Partial Charge parameters for Element: {elem} Mode: {mode}"
        )));
    }
    lookup("X", "*").ok_or_else(|| {
        RDKitError::Value("ERROR: Default Gasteiger Partial Charge parameters are missing".into())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_and_unknown_parameters() {
        assert_eq!(get_params("C", "sp3", true).unwrap(), [7.98, 9.18, 1.88]);
        assert_eq!(
            get_params("Al", "sp2", true).unwrap(),
            [5.795, 5.020, 0.695]
        );
        assert_eq!(get_params("Xe", "sp3", false).unwrap(), [0.0, 0.0, 0.0]);
        let err = get_params("Xe", "sp3", true).unwrap_err();
        assert!(err.to_string().contains("Element: Xe Mode: sp3"));
    }
}
