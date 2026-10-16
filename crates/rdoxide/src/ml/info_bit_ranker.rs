//! Ranking the bits of labelled fingerprints by information gain or
//! chi-squared (`InfoBitRanker` of `ML/InfoTheory/InfoBitRanker.{h,cpp}`),
//! and the pairwise correlation counts of a list of bits
//! (`BitCorrMatGenerator` of `CorrMatGenerator.h`).
//!
//! Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::info_theory::{chi_square, info_entropy_gain};
use crate::datastructs::bit_vect::{BitVect, ExplicitBitVect};
use crate::rdgeneral::cxx_sort::LibcxxPriorityQueue;
use crate::rdgeneral::strings::format_g;
use crate::{RDKitError, Result};

const RANKER_FILE: &str = "Code/ML/InfoTheory/InfoBitRanker.cpp";

/// `InfoBitRanker::InfoType`: the measure of information.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InfoType {
    Entropy = 1,
    BiasEntropy = 2,
    ChiSquare = 3,
    BiasChiSquare = 4,
}

/// Ranks bits by how well they separate the classes of the fingerprints
/// whose votes are accumulated.
#[derive(Debug, Clone)]
pub struct InfoBitRanker {
    dims: u32,
    classes: u32,
    info_type: InfoType,
    /// for each class, the number of fingerprints that hit each bit
    counts: Vec<Vec<u16>>,
    /// the number of fingerprints of each class
    class_count: Vec<u16>,
    /// the ranked bits: `[bit, info, hits of class 0, hits of class 1, ...]`
    top_bits: Vec<f64>,
    top: u32,
    n_inst: u32,
    bias_list: Vec<i32>,
    mask_bits: Option<ExplicitBitVect>,
}

impl InfoBitRanker {
    pub fn new(n_bits: u32, n_classes: u32, info_type: InfoType) -> Self {
        InfoBitRanker {
            dims: n_bits,
            classes: n_classes,
            info_type,
            counts: vec![vec![0; n_bits as usize]; n_classes as usize],
            class_count: vec![0; n_classes as usize],
            top_bits: Vec::new(),
            top: 0,
            n_inst: 0,
            bias_list: Vec::new(),
            mask_bits: None,
        }
    }

    pub fn num_instances(&self) -> u32 {
        self.n_inst
    }

    pub fn num_classes(&self) -> u32 {
        self.classes
    }

    /// Sets the classes the ranking is biased towards (for the biased
    /// measures): a bit must hit a larger fraction of one of these classes
    /// than of any other class.
    pub fn set_bias_list(&mut self, class_list: &[i32]) -> Result<()> {
        if class_list.len() as u32 >= self.classes {
            return Err(RDKitError::range_at(
                "classList.size()",
                class_list.len(),
                self.classes,
                RANKER_FILE,
                33,
            ));
        }
        let mut list = class_list.to_vec();
        list.sort_unstable();
        if list.windows(2).any(|w| w[0] == w[1]) {
            return Err(RDKitError::invariant_at(
                "There are duplicates in the class bias list",
                "bi == d_biasList.end()",
                RANKER_FILE,
                39,
            ));
        }
        // finally make sure all the class ID are within range
        for &c in &list {
            if c as u32 >= self.classes {
                return Err(RDKitError::range_at(
                    "static_cast<unsigned int>(*bi)",
                    c as u32,
                    self.classes,
                    RANKER_FILE,
                    43,
                ));
            }
        }
        self.bias_list = list;
        Ok(())
    }

    /// Restricts the ranking to the given bits.
    pub fn set_mask_bits(&mut self, mask_bits: &[i32]) -> Result<()> {
        let mut mask = ExplicitBitVect::new(self.dims, false);
        for &b in mask_bits {
            mask.set_bit(b as u32)?;
        }
        self.mask_bits = Some(mask);
        Ok(())
    }

    /// Counts a fingerprint of the class `label` for each of its on bits.
    pub fn accumulate_votes<B: BitVect>(&mut self, bv: &B, label: u32) -> Result<()> {
        if label >= self.classes {
            return Err(RDKitError::range_at(
                "label",
                label,
                self.classes,
                RANKER_FILE,
                115,
            ));
        }
        if bv.num_bits() != self.dims {
            return Err(RDKitError::invariant_at(
                "Incorrect bit vector size",
                "bv.getNumBits() == d_dims",
                RANKER_FILE,
                116,
            ));
        }
        self.n_inst += 1;
        let cls = label as usize;
        self.class_count[cls] = self.class_count[cls].wrapping_add(1);
        for bit in bv.on_bits() {
            if self
                .mask_bits
                .as_ref()
                .is_none_or(|m| m.get_bit(bit as u32).unwrap_or(false))
            {
                let c = &mut self.counts[cls][bit as usize];
                *c = c.wrapping_add(1);
            }
        }
        Ok(())
    }

    /// `BiasCheckBit`: whether the bit is good enough to be scored, given
    /// the bias list. `res_mat` holds the hits of each class, then the misses.
    fn bias_check_bit(&self, res_mat: &[u16]) -> bool {
        if self.bias_list.is_empty() || self.bias_list.len() as u32 == self.classes {
            return true;
        }
        let classes = self.classes as usize;
        let mut fracs = vec![0.0f64; classes];
        let mut max_cor = 0.0f64;
        for i in 0..classes {
            fracs[i] = if self.class_count[i] > 0 {
                f64::from(res_mat[i]) / f64::from(self.class_count[i])
            } else {
                0.0
            };
            if !self.bias_list.contains(&(i as i32)) && fracs[i] > max_cor {
                max_cor = fracs[i];
            }
        }
        self.bias_list.iter().any(|&b| fracs[b as usize] >= max_cor)
    }

    fn info_of(&self, res_mat: &[u16]) -> f64 {
        let classes = self.classes as usize;
        match self.info_type {
            InfoType::Entropy => info_entropy_gain(res_mat, 2, classes),
            InfoType::BiasEntropy => {
                if self.bias_check_bit(res_mat) {
                    info_entropy_gain(res_mat, 2, classes)
                } else {
                    0.0
                }
            }
            InfoType::ChiSquare => chi_square(res_mat, 2, classes),
            InfoType::BiasChiSquare => {
                if self.bias_check_bit(res_mat) {
                    chi_square(res_mat, 2, classes)
                } else {
                    0.0
                }
            }
        }
    }

    /// Ranks the bits and returns the top `num` as rows of
    /// `[bit, info, hits of class 0, hits of class 1, ...]`, best first.
    pub fn get_top_n(&mut self, num: u32) -> Result<&[f64]> {
        if num > self.dims {
            return Err(RDKitError::Value(
                "attempt to rank more bits than present in the bit vectors".into(),
            ));
        }
        if let Some(mask) = &self.mask_bits
            && num > mask.num_on_bits()
        {
            return Err(RDKitError::invariant_at(
                "Can't rank more bits than the ensemble size",
                "num <= dp_maskBits->getNumOnBits()",
                RANKER_FILE,
                155,
            ));
        }
        let classes = self.classes as usize;
        let mut res_mat = vec![0u16; 2 * classes];
        // the entry with the smallest info is on top
        let mut top_n = LibcxxPriorityQueue::new(|a: &(f64, i32), b: &(f64, i32)| a.0 > b.0);
        for i in 0..self.dims as usize {
            if self
                .mask_bits
                .as_ref()
                .is_some_and(|m| !m.get_bit(i as u32).unwrap_or(false))
            {
                continue;
            }
            for j in 0..classes {
                res_mat[j] = self.counts[j][i];
                res_mat[classes + j] = self.class_count[j].wrapping_sub(self.counts[j][i]);
            }
            let info = self.info_of(&res_mat);
            let entry = (info, i as i32);
            if info >= 0.0 {
                if top_n.len() < num as usize {
                    top_n.push(entry);
                } else if info > top_n.top().map_or(f64::INFINITY, |t| t.0) {
                    top_n.pop();
                    top_n.push(entry);
                }
            }
        }

        self.top = num;
        let ncols = 2 + classes;
        let mut top_bits = vec![0.0f64; num as usize * ncols];
        let mask_list = match &self.mask_bits {
            Some(m) if top_n.len() < num as usize => m.on_bits(),
            _ => Vec::new(),
        };
        for i in (0..num as usize).rev() {
            let offset = i * ncols;
            let bid = match top_n.top().copied() {
                None => {
                    top_bits[offset + 1] = 0.0;
                    if self.mask_bits.is_some() {
                        mask_list[i] as usize
                    } else {
                        i
                    }
                }
                Some((info, bit)) => {
                    top_bits[offset + 1] = info;
                    top_n.pop();
                    bit as usize
                }
            };
            top_bits[offset] = bid as f64;
            for j in 0..classes {
                top_bits[offset + 2 + j] = f64::from(self.counts[j][bid]);
            }
        }
        self.top_bits = top_bits;
        Ok(&self.top_bits)
    }

    /// `writeTopBitsToStream`: the ranked bits as a table.
    pub fn top_bits_table(&self) -> String {
        let mut out = format!("{:>12}{:>12}", "Bit", "InfoContent");
        for ic in 0..self.classes {
            out += &format!("{:>10}{ic}", "class");
        }
        out.push('\n');
        let ncols = 2 + self.classes as usize;
        for i in 0..self.top as usize {
            out += &format!(
                "{:>12}{:>12}",
                self.top_bits[i * ncols] as i32,
                format_g(self.top_bits[i * ncols + 1], 5)
            );
            for ic in 0..self.classes as usize {
                out += &format!("{:>10}", self.top_bits[i * ncols + 2 + ic] as i32);
            }
            out.push('\n');
        }
        out
    }

    /// `writeTopBitsToFile`
    pub fn write_top_bits_to_file(&self, file_name: &str) -> Result<()> {
        std::fs::write(file_name, self.top_bits_table())
            .map_err(|_| RDKitError::FileParse(format!("Bad output file {file_name}")))
    }
}

/// Counts, for each pair of a list of bits, how many fingerprints have both
/// on. The result is the lower triangle of the symmetric matrix as one
/// vector.
#[derive(Debug, Clone, Default)]
pub struct BitCorrMatGenerator {
    descs: Vec<i32>,
    corr_mat: Vec<f64>,
    n_examples: i32,
}

impl BitCorrMatGenerator {
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the bits to correlate (e.g. the top ranked ones) and clears the
    /// counts.
    pub fn set_bit_id_list(&mut self, bit_id_list: &[i32]) {
        self.descs = bit_id_list.to_vec();
        let nd = self.descs.len();
        self.corr_mat = vec![0.0; nd * nd.saturating_sub(1) / 2];
    }

    pub fn num_examples(&self) -> i32 {
        self.n_examples
    }

    pub fn corr_bit_list(&self) -> &[i32] {
        &self.descs
    }

    pub fn corr_mat(&self) -> &[f64] {
        &self.corr_mat
    }

    /// For each pair of on bits (bi, bj) of the listed bits in `fp`,
    /// increases the count of the pair by 1.
    pub fn collect_votes<B: BitVect>(&mut self, fp: &B) -> Result<()> {
        let nd = self.descs.len();
        let mut on = Vec::with_capacity(nd);
        for &bi in &self.descs {
            on.push(fp.get_bit(bi as u32)?);
        }
        for i in 1..nd {
            if on[i] {
                let itab = i * (i - 1) / 2;
                for j in 0..i {
                    if on[j] {
                        self.corr_mat[itab + j] += 1.0;
                    }
                }
            }
        }
        self.n_examples += 1;
        Ok(())
    }
}
