//! JSON (de)serialization of [`EmbedParameters`] (upstream
//! `DistGeomHelpers/EmbedderUtils.cpp`). Like upstream's property tree, every
//! value is written as a string; the callback and CPCI are not serialized.

use super::embed_params::{EmbedFF, EmbedParameters, InitialEmbeddingMode};
use crate::error::{RDKitError, Result};
use crate::geometry::Point3D;
use crate::rdgeneral::ptree::PTree;

/// Text form of a parameter value, as a C++ stream writes it.
trait ParamValue: Sized {
    fn to_text(&self) -> String;
    fn from_text(text: &str) -> Option<Self>;
}

macro_rules! int_param {
    ($($t:ty),*) => {$(
        impl ParamValue for $t {
            fn to_text(&self) -> String {
                self.to_string()
            }
            fn from_text(text: &str) -> Option<Self> {
                text.trim().parse().ok()
            }
        }
    )*};
}
int_param!(u32, i32);

impl ParamValue for bool {
    fn to_text(&self) -> String {
        self.to_string()
    }
    fn from_text(text: &str) -> Option<Self> {
        match text.trim() {
            "true" | "1" => Some(true),
            "false" | "0" => Some(false),
            _ => None,
        }
    }
}

impl ParamValue for f64 {
    fn to_text(&self) -> String {
        format_g17(*self)
    }
    fn from_text(text: &str) -> Option<Self> {
        text.trim().parse().ok()
    }
}

impl ParamValue for EmbedFF {
    fn to_text(&self) -> String {
        match self {
            EmbedFF::Uff => "UFF",
            EmbedFF::Mmff => "MMFF",
        }
        .to_string()
    }
    fn from_text(text: &str) -> Option<Self> {
        match text.trim() {
            "UFF" => Some(EmbedFF::Uff),
            "MMFF" => Some(EmbedFF::Mmff),
            _ => None,
        }
    }
}

impl ParamValue for InitialEmbeddingMode {
    fn to_text(&self) -> String {
        match self {
            InitialEmbeddingMode::DgEmbedding => "DG_EMBEDDING",
            InitialEmbeddingMode::InternalCoordinateEmbedding => "INTERNAL_COORDINATE_EMBEDDING",
            InitialEmbeddingMode::RandomCoordinateEmbedding => "RANDOM_COORDINATE_EMBEDDING",
        }
        .to_string()
    }
    fn from_text(text: &str) -> Option<Self> {
        match text.trim() {
            "DG_EMBEDDING" => Some(InitialEmbeddingMode::DgEmbedding),
            "INTERNAL_COORDINATE_EMBEDDING" => {
                Some(InitialEmbeddingMode::InternalCoordinateEmbedding)
            }
            "RANDOM_COORDINATE_EMBEDDING" => Some(InitialEmbeddingMode::RandomCoordinateEmbedding),
            _ => None,
        }
    }
}

/// `printf("%.17g")`: the format boost's property tree uses for doubles.
fn format_g17(value: f64) -> String {
    if value == 0.0 {
        return if value.is_sign_negative() { "-0" } else { "0" }.to_string();
    }
    if !value.is_finite() {
        return if value.is_nan() {
            "nan"
        } else if value > 0.0 {
            "inf"
        } else {
            "-inf"
        }
        .to_string();
    }
    // 17 significant digits: d.dddddddddddddddde[+-]x
    let sci = format!("{value:.16e}");
    let (mantissa, exp) = sci.split_once('e').expect("exponent");
    let exp: i32 = exp.parse().expect("exponent");
    let (sign, digits) = match mantissa.strip_prefix('-') {
        Some(rest) => ("-", rest.replace('.', "")),
        None => ("", mantissa.replace('.', "")),
    };
    let digits = digits.trim_end_matches('0');
    let digits = if digits.is_empty() { "0" } else { digits };
    if !(-4..17).contains(&exp) {
        let frac = &digits[1..];
        let dot = if frac.is_empty() { "" } else { "." };
        let exp_sign = if exp < 0 { '-' } else { '+' };
        return format!(
            "{sign}{}{dot}{frac}e{exp_sign}{:02}",
            &digits[..1],
            exp.abs()
        );
    }
    if exp < 0 {
        format!("{sign}0.{}{digits}", "0".repeat((-exp - 1) as usize))
    } else {
        let int_len = exp as usize + 1;
        if digits.len() <= int_len {
            format!("{sign}{digits}{}", "0".repeat(int_len - digits.len()))
        } else {
            format!("{sign}{}.{}", &digits[..int_len], &digits[int_len..])
        }
    }
}

/// Applies `$m!(json_name, rust_field)` to every plain parameter, in
/// upstream's order.
macro_rules! embed_params_fields {
    ($m:ident) => {
        $m!(basinThresh, basin_thresh);
        $m!(boundsMatForceScaling, bounds_mat_force_scaling);
        $m!(boxSizeMult, box_size_mult);
        $m!(clearConfs, clear_confs);
        $m!(embedForceField, embed_force_field);
        $m!(embedFragmentsSeparately, embed_fragments_separately);
        $m!(enableSequentialRandomSeeds, enable_sequential_random_seeds);
        $m!(enforceChirality, enforce_chirality);
        $m!(ETversion, et_version);
        $m!(forceTransAmides, force_trans_amides);
        $m!(ignoreSmoothingFailures, ignore_smoothing_failures);
        $m!(maxIterations, max_iterations);
        $m!(numThreads, num_threads);
        $m!(numZeroFail, num_zero_fail);
        $m!(onlyHeavyAtomsForRMS, only_heavy_atoms_for_rms);
        $m!(optimizerForceTol, optimizer_force_tol);
        $m!(pruneRmsThresh, prune_rms_thresh);
        $m!(randNegEig, rand_neg_eig);
        $m!(randomSeed, random_seed);
        $m!(
            symmetrizeConjugatedTerminalGroupsForPruning,
            symmetrize_conjugated_terminal_groups_for_pruning
        );
        $m!(timeout, timeout);
        $m!(trackFailures, track_failures);
        $m!(useBasicKnowledge, use_basic_knowledge);
        $m!(useExpTorsionAnglePrefs, use_exp_torsion_angle_prefs);
        $m!(useLegacyImplementation, use_legacy_implementation);
        $m!(useMacrocycle14config, use_macrocycle_14_config);
        $m!(useMacrocycleTorsions, use_macrocycle_torsions);
        $m!(useRandomCoords, use_random_coords);
        $m!(useSmallRingTorsions, use_small_ring_torsions);
        $m!(useSymmetryForPruning, use_symmetry_for_pruning);
        $m!(verbose, verbose);
        $m!(onlyInitialEmbedding, only_initial_embedding);
    };
}

/// Serializes the parameters (and the coordinate map and bounds matrix, when
/// set) to a compact JSON string.
pub fn embed_parameters_to_json(params: &EmbedParameters) -> String {
    let mut pt = PTree::new();
    macro_rules! put {
        ($name:ident, $field:ident) => {
            pt.put(stringify!($name), params.$field.to_text());
        };
    }
    embed_params_fields!(put);
    pt.put(
        "initialEmbeddingMode",
        params.initial_embedding_mode.to_text(),
    );
    if let Some(coord_map) = &params.coord_map {
        let mut coord_pt = PTree::new();
        for (idx, p) in coord_map {
            let mut point = PTree::new();
            for v in [p.x, p.y, p.z] {
                point.push(PTree::leaf(format!("{v:.6}")));
            }
            coord_pt.add_child(&idx.to_string(), point);
        }
        pt.add_child("coordMap", coord_pt);
    }
    if let Some(mat) = &params.bounds_mat {
        let n = mat.num_rows();
        let mut matrix = PTree::new();
        for i in 0..n {
            let mut row = PTree::new();
            for j in 0..n {
                row.push(PTree::leaf(mat.val(i, j).to_text()));
            }
            matrix.push(row);
        }
        pt.add_child("boundsMatrix", matrix);
    }
    pt.to_json()
}

/// Updates `params` from a JSON string; keys that are absent (or whose
/// value cannot be converted) keep their current value.
pub fn update_embed_parameters_from_json(params: &mut EmbedParameters, json: &str) -> Result<()> {
    if json.is_empty() {
        return Ok(());
    }
    let pt = PTree::from_json(json)?;
    macro_rules! get {
        ($name:ident, $field:ident) => {
            if let Some(v) = pt
                .get_str(stringify!($name))
                .and_then(ParamValue::from_text)
            {
                params.$field = v;
            }
        };
    }
    embed_params_fields!(get);
    if pt.get_str("numThreads").is_some() {
        params.num_threads_is_default = false;
    }

    if let Some(entries) = pt.child("coordMap") {
        let mut map = std::collections::BTreeMap::new();
        for (key, value) in entries.members() {
            let idx: i32 = key
                .trim()
                .parse()
                .map_err(|_| RDKitError::Value(format!("bad coordMap index '{key}'")))?;
            let coords: Vec<Option<f32>> = value
                .children()
                .map(|c| c.value().trim().parse().ok())
                .collect();
            match coords[..] {
                [Some(x), Some(y), Some(z)] => {
                    map.insert(idx, Point3D::new(x.into(), y.into(), z.into()));
                }
                _ => {
                    return Err(RDKitError::Value(
                        "coordMap entries need three numbers".into(),
                    ));
                }
            }
        }
        params.coord_map = Some(map);
    }
    if let Some(text) = pt.get_str("initialEmbeddingMode") {
        params.initial_embedding_mode = ParamValue::from_text(text).ok_or_else(|| {
            RDKitError::Value(
                "Initial embedding option must be either 'DG_EMBEDDING', \
                 'INTERNAL_COORDINATE_EMBEDDING', or 'RANDOM_COORDINATE_EMBEDDING'"
                    .into(),
            )
        })?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g17_matches_printf() {
        assert_eq!(format_g17(5.0), "5");
        assert_eq!(format_g17(0.001), "0.001");
        assert_eq!(format_g17(-1.0), "-1");
        assert_eq!(format_g17(0.980_254_204_001_361_5), "0.98025420400136154");
        assert_eq!(format_g17(1e-5), "1.0000000000000001e-05");
        assert_eq!(format_g17(1e20), "1e+20");
        assert_eq!(format_g17(123456.0), "123456");
    }

    #[test]
    fn round_trip() {
        let mut p = EmbedParameters::etkdg_v3();
        p.random_seed = 42;
        p.initial_embedding_mode = InitialEmbeddingMode::RandomCoordinateEmbedding;
        let mut q = EmbedParameters::default();
        update_embed_parameters_from_json(&mut q, &embed_parameters_to_json(&p)).unwrap();
        assert_eq!(embed_parameters_to_json(&p), embed_parameters_to_json(&q));
        assert_eq!(q.random_seed, 42);
        let mut r = EmbedParameters::default();
        update_embed_parameters_from_json(&mut r, r#"{"randomSeed": 7, "clearConfs": false}"#)
            .unwrap();
        assert_eq!((r.random_seed, r.clear_confs), (7, false));
    }
}
