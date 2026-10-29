//! Port of FreeSASA 2.0.3 (Simon Mitternacht, MIT; see `LICENSES/FreeSASA-MIT.txt`): solvent accessible surface areas by
//! Lee & Richards' or Shrake & Rupley's algorithm and the atom classifiers RDKit's `External/FreeSASA` uses.

pub mod lr;
pub mod nb;
pub mod sr;

/// `freesasa_algorithm`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Algorithm {
    LeeRichards,
    ShrakeRupley,
}

/// `freesasa_parameters` (the thread count is not needed: the results do not depend on it)
#[derive(Debug, Clone, Copy)]
pub struct Parameters {
    pub alg: Algorithm,
    pub probe_radius: f64,
    pub shrake_rupley_n_points: usize,
    pub lee_richards_n_slices: usize,
}

impl Default for Parameters {
    fn default() -> Self {
        Parameters {
            alg: Algorithm::LeeRichards,
            probe_radius: 1.4,
            shrake_rupley_n_points: 100,
            lee_richards_n_slices: 20,
        }
    }
}

/// `freesasa_result`
#[derive(Debug, Clone)]
pub struct SasaResult {
    pub total: f64,
    pub sasa: Vec<f64>,
}

/// `freesasa_calc_coord`: `xyz` holds x, y, z of every atom, `radii` the van der Waals radii
pub fn calc_coord(xyz: &[f64], radii: &[f64], params: &Parameters) -> SasaResult {
    assert_eq!(xyz.len(), 3 * radii.len());
    assert!(!radii.is_empty());
    let sasa = match params.alg {
        Algorithm::ShrakeRupley => sr::shrake_rupley(
            xyz,
            radii,
            params.probe_radius,
            params.shrake_rupley_n_points,
        ),
        Algorithm::LeeRichards => lr::lee_richards(
            xyz,
            radii,
            params.probe_radius,
            params.lee_richards_n_slices,
        ),
    };
    let mut total = 0.0;
    for s in &sasa {
        total += s;
    }
    SasaResult { total, sasa }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex(s: &str) -> f64 {
        f64::from_bits(u64::from_str_radix(s, 16).unwrap())
    }

    /// areas equal those of freesasa 2.0.3 (`tests/data/freesasa_sasa.txt`, `tools/freesasa_drv/sasa.c`); bit for bit
    /// where the platform's libm agrees with the one that made the corpus (macOS), within 1e-9 relative elsewhere
    #[test]
    fn sasa_matches_c() {
        let text = include_str!("../../../../../tests/data/freesasa_sasa.txt");
        let strict = cfg!(target_os = "macos");
        let close = |a: f64, b: f64| {
            if strict {
                a == b
            } else {
                (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
            }
        };
        let mut n = 0;
        for rec in text.split("$$\n").filter(|r| !r.is_empty()) {
            let mut lines = rec.lines();
            let f: Vec<&str> = lines.next().unwrap().split(' ').collect();
            let (alg, count, probe) = (f[1] == "1", f[2].parse::<usize>().unwrap(), hex(f[3]));
            let (mut xyz, mut radii) = (Vec::new(), Vec::new());
            for _ in 0..count {
                let v: Vec<f64> = lines.next().unwrap().split(' ').map(hex).collect();
                xyz.extend_from_slice(&v[..3]);
                radii.push(v[3]);
            }
            let params = Parameters {
                alg: if alg {
                    Algorithm::ShrakeRupley
                } else {
                    Algorithm::LeeRichards
                },
                probe_radius: probe,
                ..Default::default()
            };
            let got = calc_coord(&xyz, &radii, &params);
            let head: Vec<&str> = lines.next().unwrap().split(' ').collect();
            assert!(
                close(got.total, hex(head[3])),
                "total {} vs {}",
                got.total,
                hex(head[3])
            );
            for (i, l) in lines.enumerate() {
                let w = hex(l.split(' ').nth(2).unwrap());
                assert!(
                    close(got.sasa[i], w),
                    "case {n} atom {i}: {} vs {w}",
                    got.sasa[i]
                );
            }
            n += 1;
        }
        assert!(n > 60);
    }
}
