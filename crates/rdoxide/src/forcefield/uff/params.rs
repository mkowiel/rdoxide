//! UFF atomic parameters (upstream `ForceField/UFF/Params.{h,cpp}`).

use std::collections::BTreeMap;
use std::sync::OnceLock;

/// The parameters of one UFF atom type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtomicParams {
    /// valence bond radius
    pub r1: f64,
    /// valence angle (radians)
    pub theta0: f64,
    /// vdW characteristic length
    pub x1: f64,
    /// vdW atomic energy
    pub d1: f64,
    /// vdW scaling term
    pub zeta: f64,
    /// effective charge
    pub z1: f64,
    /// sp3 torsional barrier parameter
    pub v1: f64,
    /// torsional contribution for sp2-sp3 bonds
    pub u1: f64,
    /// GMP electronegativity
    pub gmp_xi: f64,
    /// GMP hardness
    pub gmp_hardness: f64,
    /// GMP radius
    pub gmp_radius: f64,
}

/// scaling factor for the rBO correction
pub const LAMBDA: f64 = 0.1332;
/// bond force constant prefactor
pub const G: f64 = 332.06;
/// special case bond order for amide C-N bonds
pub const AMIDE_BOND_ORDER: f64 = 1.41;

pub const DEG2RAD: f64 = std::f64::consts::PI / 180.0;
pub const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;

pub(crate) fn is_double_zero(x: f64) -> bool {
    x < 1.0e-10 && x > -1.0e-10
}

pub(crate) fn clip_to_one(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

/// The parameters by UFF atom label (`ParamCollection`).
#[derive(Debug, Clone)]
pub struct ParamCollection {
    params: BTreeMap<String, AtomicParams>,
}

impl ParamCollection {
    /// The default UFF parameters.
    pub fn default_params() -> &'static ParamCollection {
        static DEFAULT: OnceLock<ParamCollection> = OnceLock::new();
        DEFAULT.get_or_init(|| {
            ParamCollection::parse(include_str!("../../../data/uff_params.tsv"))
                .expect("the bundled UFF parameters are valid")
        })
    }

    /// Parses tab-separated parameter data (`#` starts a comment line).
    pub fn parse(data: &str) -> Result<ParamCollection, String> {
        let mut params = BTreeMap::new();
        for line in data.lines() {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut fields = line.split('\t');
            let label = fields.next().unwrap_or_default().to_string();
            let mut num = || -> Result<f64, String> {
                fields
                    .next()
                    .ok_or_else(|| format!("too few fields in UFF parameter line '{line}'"))?
                    .trim()
                    .parse()
                    .map_err(|_| format!("bad number in UFF parameter line '{line}'"))
            };
            let r1 = num()?;
            let theta0 = num()? * std::f64::consts::PI / 180.;
            let p = AtomicParams {
                r1,
                theta0,
                x1: num()?,
                d1: num()?,
                zeta: num()?,
                z1: num()?,
                v1: num()?,
                u1: num()?,
                gmp_xi: num()?,
                gmp_hardness: num()?,
                gmp_radius: num()?,
            };
            params.insert(label, p);
        }
        Ok(ParamCollection { params })
    }

    /// The parameters of a UFF atom label like `C_3`.
    pub fn get(&self, symbol: &str) -> Option<&AtomicParams> {
        self.params.get(symbol)
    }

    pub fn len(&self) -> usize {
        self.params.len()
    }

    pub fn is_empty(&self) -> bool {
        self.params.is_empty()
    }
}
