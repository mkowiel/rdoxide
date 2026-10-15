//! Scoring of R group assignments by the variance of the fingerprints of the
//! R groups (upstream `RGroupDecomposition/RGroupFingerprintScore.{h,cpp}`).

use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;

use super::data::{FINGERPRINT_SIZE, RGroupData, RGroupFingerprint, RGroupMatch};
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::{RDKitError, Result};

/// The fingerprint of a hydrogen R group, used for molecules that have no R
/// group at a label.
fn dummy_hydrogen_fingerprint() -> Result<&'static RGroupFingerprint> {
    static FP: OnceLock<RGroupFingerprint> = OnceLock::new();
    if let Some(fp) = FP.get() {
        return Ok(fp);
    }
    let mol = mol_from_smiles("*[H]", &SmilesParserParams::default())?
        .ok_or_else(|| RDKitError::invariant("could not parse *[H]"))?;
    let mut data = RGroupData::new();
    data.add(&std::sync::Arc::new(mol), &[0])?;
    let fp = data.fingerprint()?.clone();
    Ok(FP.get_or_init(|| fp))
}

/// `VarianceDataForLabel`: bit counts of the fingerprints of one label.
#[derive(Debug, Clone)]
pub struct VarianceDataForLabel {
    pub number_fingerprints: i32,
    pub bit_counts: Vec<i32>,
}

impl Default for VarianceDataForLabel {
    fn default() -> Self {
        Self::new()
    }
}

impl VarianceDataForLabel {
    pub fn new() -> Self {
        VarianceDataForLabel {
            number_fingerprints: 0,
            bit_counts: vec![0; FINGERPRINT_SIZE],
        }
    }

    pub fn add_fingerprint(&mut self, fp: &RGroupFingerprint) {
        self.number_fingerprints += 1;
        for &b in &fp.on_bits {
            self.bit_counts[b as usize] += 1;
        }
    }

    pub fn remove_fingerprint(&mut self, fp: &RGroupFingerprint) {
        self.number_fingerprints -= 1;
        for &b in &fp.on_bits {
            self.bit_counts[b as usize] -= 1;
        }
    }

    /// `variance`
    pub fn variance(&self) -> f64 {
        let n = f64::from(self.number_fingerprints);
        self.bit_counts.iter().fold(0.0, |sum, &count| {
            if count == 0 {
                return sum;
            }
            let bit_count = n / 2.0 + f64::from(count) / 2.0;
            let ss = bit_count - (bit_count * bit_count) / n;
            sum + ss / n
        })
    }
}

/// `FingerprintVarianceScoreData`
#[derive(Debug, Clone, Default)]
pub struct FingerprintVarianceScoreData {
    pub number_of_missing_user_rgroups: usize,
    pub number_of_molecules: usize,
    pub labels_to_variance_data: BTreeMap<i32, VarianceDataForLabel>,
}

impl FingerprintVarianceScoreData {
    /// `fingerprintVarianceGroupScore`
    pub fn fingerprint_variance_group_score(&self) -> Result<f64> {
        let sum: f64 = self
            .labels_to_variance_data
            .values()
            .map(VarianceDataForLabel::variance)
            .sum();
        if self.number_of_molecules == 0 {
            return Err(RDKitError::invariant("No compounds to be scored!"));
        }
        Ok(-sum.sqrt())
    }

    fn modify_variance_data(
        &mut self,
        match_number: usize,
        permutation_number: usize,
        matches: &[Vec<RGroupMatch>],
        labels: &BTreeSet<i32>,
        add: bool,
    ) -> Result<()> {
        let m = &matches[match_number][permutation_number];
        for &l in labels {
            let data = self.labels_to_variance_data.entry(l).or_default();
            let apply = |data: &mut VarianceDataForLabel, fp: &RGroupFingerprint| {
                if add {
                    data.add_fingerprint(fp);
                } else {
                    data.remove_fingerprint(fp);
                }
            };
            match m.rgroups.get(&l) {
                Some(rg) => apply(data, rg.read().expect("rgroup lock").fingerprint()?),
                None => apply(data, dummy_hydrogen_fingerprint()?),
            }
        }
        if add {
            self.number_of_missing_user_rgroups += m.number_missing_user_rgroups;
            self.number_of_molecules += 1;
        } else {
            self.number_of_missing_user_rgroups -= m.number_missing_user_rgroups;
            self.number_of_molecules -= 1;
        }
        Ok(())
    }

    pub fn add_variance_data(
        &mut self,
        match_number: usize,
        permutation_number: usize,
        matches: &[Vec<RGroupMatch>],
        labels: &BTreeSet<i32>,
    ) -> Result<()> {
        self.modify_variance_data(match_number, permutation_number, matches, labels, true)
    }

    pub fn remove_variance_data(
        &mut self,
        match_number: usize,
        permutation_number: usize,
        matches: &[Vec<RGroupMatch>],
        labels: &BTreeSet<i32>,
    ) -> Result<()> {
        self.modify_variance_data(match_number, permutation_number, matches, labels, false)
    }

    pub fn clear(&mut self) {
        *self = FingerprintVarianceScoreData::default();
    }
}

/// `fingerprintVarianceScore`: the score of `permutation`, accumulated into
/// `score_data` when given.
pub fn fingerprint_variance_score(
    permutation: &[usize],
    matches: &[Vec<RGroupMatch>],
    labels: &BTreeSet<i32>,
    score_data: Option<&mut FingerprintVarianceScoreData>,
) -> Result<f64> {
    if permutation.len() > matches.len() {
        return Err(RDKitError::invariant(""));
    }
    let mut local = FingerprintVarianceScoreData::default();
    let data = score_data.unwrap_or(&mut local);
    for &l in labels {
        let for_label = data.labels_to_variance_data.entry(l).or_default();
        for (m, &p) in permutation.iter().enumerate() {
            match matches[m][p].rgroups.get(&l) {
                Some(rg) => {
                    for_label.add_fingerprint(rg.read().expect("rgroup lock").fingerprint()?)
                }
                None => for_label.add_fingerprint(dummy_hydrogen_fingerprint()?),
            }
        }
    }
    let missing: usize = permutation
        .iter()
        .enumerate()
        .map(|(m, &p)| matches[m][p].number_missing_user_rgroups)
        .sum();
    data.number_of_missing_user_rgroups += missing;
    data.number_of_molecules += permutation.len();
    data.fingerprint_variance_group_score()
}
