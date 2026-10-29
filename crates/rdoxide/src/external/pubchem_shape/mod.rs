//! Port of RDKit's `External/pubchem_shape` and the NCBI pubchem-align3d library it wraps (public domain; see
//! `LICENSES/PubChemAlign3D-PublicDomain.txt` and `docs/PUBCHEM_SHAPE_PORT.md`): shape overlays with Gaussian volumes and
//! colour features, optimized from up to 16 starting poses.

pub mod consts;
pub mod functions;
pub mod neighbor;
pub mod shape;

#[cfg(test)]
mod tests {
    use super::functions::*;
    use super::neighbor::*;

    fn dhex(s: &str) -> f64 {
        f64::from_bits(u64::from_str_radix(s, 16).unwrap())
    }

    fn fhex(s: &str) -> f32 {
        f32::from_bits(u32::from_str_radix(s, 16).unwrap())
    }

    struct Side {
        coord: Vec<f32>,
        alpha: Vec<f64>,
        types: Vec<u32>,
        vol: Vec<u32>,
        cmap: ColorMap,
        sov: f64,
        sof: f64,
    }

    fn side(lines: &[&str]) -> Side {
        let (mut coord, mut rad, mut types) = (Vec::new(), Vec::new(), Vec::new());
        for l in lines {
            let f: Vec<&str> = l.split(' ').collect();
            coord.extend([fhex(f[0]), fhex(f[1]), fhex(f[2])]);
            rad.push(dhex(f[3]));
            types.push(f[4].parse().unwrap());
        }
        let alpha = set_alpha(&rad);
        let vol = get_volume_atom_index_vector(&types);
        let cmap = get_color_atom_type_to_index_vector_map(&types);
        let sov = shape_overlap(&coord, &alpha, &vol, &coord, &alpha, &vol);
        let sof = feature_overlap(&coord, &alpha, &cmap, &coord, &alpha, &cmap, None);
        Side {
            coord,
            alpha,
            types,
            vol,
            cmap,
            sov,
            sof,
        }
    }

    /// volumes, steric frame classes and the optimized overlay equal those of pubchem-align3d
    /// (`tests/data/pubchem_shape_align.txt`, `tools/pcshape_drv/align.cpp`, `tools/make_pubchem_shape_corpus.py`)
    #[test]
    fn align_matches_c() {
        let text = include_str!("../../../../../tests/data/pubchem_shape_align.txt");
        let strict = cfg!(target_os = "macos");
        let close = |a: f64, b: f64| {
            if strict {
                a == b
            } else {
                (a - b).abs() <= 1e-6 * (1.0 + a.abs().max(b.abs()))
            }
        };
        let mut n = 0;
        for rec in text.split("$$\n").filter(|r| !r.is_empty()) {
            let lines: Vec<&str> = rec.lines().collect();
            let h: Vec<&str> = lines[0].split(' ').collect();
            let (nref, nfit) = (
                h[1].parse::<usize>().unwrap(),
                h[2].parse::<usize>().unwrap(),
            );
            let (color, opt) = (h[3] == "1", dhex(h[4]));
            let (pre, post) = (h[5].parse::<u32>().unwrap(), h[6].parse::<u32>().unwrap());
            let r = side(&lines[1..1 + nref]);
            let f = side(&lines[1 + nref..1 + nref + nfit]);
            let out = &lines[1 + nref + nfit..];
            let vol: Vec<f64> = out[0].split(' ').skip(1).map(dhex).collect();
            assert!(
                close(r.sov, vol[0]) && close(r.sof, vol[1]),
                "ref volumes, case {n}"
            );
            assert!(
                close(f.sov, vol[2]) && close(f.sof, vol[3]),
                "fit volumes, case {n}"
            );
            let qrat: Vec<u32> = out[1]
                .split(' ')
                .skip(1)
                .map(|v| v.parse().unwrap())
                .collect();
            assert_eq!(qrat_of(&r.coord, 0.7225), qrat[0], "case {n}");
            assert_eq!(qrat_of(&f.coord, 0.7225), qrat[1], "case {n}");

            let joint = get_joint_color_type_set(&r.types, &f.types);
            let (mut rm, mut fm) = (r.cmap.clone(), f.cmap.clone());
            restrict_color_map(&mut rm, &joint);
            restrict_color_map(&mut fm, &joint);
            let rs = ShapeData {
                coord: &r.coord,
                alpha: &r.alpha,
                volume_atoms: &r.vol,
                color_map: &rm,
                sov: r.sov,
                sof: r.sof,
            };
            let fs = ShapeData {
                coord: &f.coord,
                alpha: &f.alpha,
                volume_atoms: &f.vol,
                color_map: &fm,
                sov: f.sov,
                sof: f.sof,
            };
            let res = neighbor_conformers(&rs, &fs, color && !joint.is_empty(), pre, post, opt);
            let m: Vec<f32> = out[2].split(' ').skip(1).map(fhex).collect();
            for (i, (&a, &b)) in res.matrix.iter().zip(&m).enumerate() {
                let ok = if strict {
                    a == b
                } else {
                    (a - b).abs() <= 1e-3
                };
                assert!(ok, "case {n} matrix[{i}]: {a} vs {b}");
            }
            let sc: Vec<f64> = out[3]
                .split(' ')
                .filter(|t| t.len() == 16)
                .map(dhex)
                .collect();
            assert!(close(res.st, sc[0]), "case {n} st {} vs {}", res.st, sc[0]);
            assert!(close(res.ct, sc[1]), "case {n} ct {} vs {}", res.ct, sc[1]);
            n += 1;
        }
        assert!(n > 100);
    }
}
