//! `freesasa/src/classifier.c`: radii and polar/apolar classes of atoms by residue and atom name (ProtOr, NACCESS, OONS).

use super::classifier_data::{self, Residue};

/// `freesasa_atom_class`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomClass {
    Apolar = 0,
    Polar = 1,
    Unknown = 2,
}

impl AtomClass {
    /// `freesasa_classifier_class2str`
    pub fn name(self) -> &'static str {
        match self {
            AtomClass::Apolar => "Apolar",
            AtomClass::Polar => "Polar",
            AtomClass::Unknown => "Unknown",
        }
    }
}

/// the built-in classifiers
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Classifier {
    Protor,
    Naccess,
    Oons,
}

impl Classifier {
    pub fn name(self) -> &'static str {
        match self {
            Classifier::Protor => classifier_data::PROTOR_NAME,
            Classifier::Naccess => classifier_data::NACCESS_NAME,
            Classifier::Oons => classifier_data::OONS_NAME,
        }
    }

    fn residues(self) -> &'static [Residue] {
        match self {
            Classifier::Protor => classifier_data::PROTOR,
            Classifier::Naccess => classifier_data::NACCESS,
            Classifier::Oons => classifier_data::OONS,
        }
    }

    /// `find_atom`: the atom of the residue, else of the residue `ANY`
    fn find(self, res_name: &str, atom_name: &str) -> Option<(f64, u8)> {
        // upstream compares the first whitespace-delimited word of the key
        let first = |s: &str| s.split_whitespace().next().unwrap_or("").to_string();
        let (res_name, atom_name) = (first(res_name), first(atom_name));
        let residues = self.residues();
        let lookup = |r: &Residue| r.1.iter().find(|a| a.0 == atom_name).map(|a| (a.1, a.2));
        let any = || residues.iter().find(|r| r.0 == "ANY").and_then(lookup);
        match residues.iter().find(|r| r.0 == res_name) {
            None => any(),
            Some(r) => lookup(r).or_else(any),
        }
    }

    /// `freesasa_classifier_radius`: -1 for an unknown atom
    pub fn radius(self, res_name: &str, atom_name: &str) -> f64 {
        self.find(res_name, atom_name).map_or(-1.0, |a| a.0)
    }

    /// `freesasa_classifier_class`
    pub fn class(self, res_name: &str, atom_name: &str) -> AtomClass {
        match self.find(res_name, atom_name) {
            Some((_, 0)) => AtomClass::Apolar,
            Some((_, 1)) => AtomClass::Polar,
            _ => AtomClass::Unknown,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// radii and classes equal those of freesasa 2.0.3 for every atom of every residue of the three classifiers and for
    /// unknown residues, atoms and the `ANY` fallback (`tests/data/freesasa_classify.txt`, `tools/freesasa_drv/classify.c`)
    #[test]
    fn classifiers_match_c() {
        let text = include_str!("../../../../../tests/data/freesasa_classify.txt");
        let mut n = 0;
        for line in text.lines() {
            let f: Vec<&str> = line.split(' ').collect();
            let c = match f[0] {
                "P" => Classifier::Protor,
                "N" => Classifier::Naccess,
                _ => Classifier::Oons,
            };
            let (radius, class): (f64, i32) = (f[3].parse().unwrap(), f[4].parse().unwrap());
            assert_eq!(c.radius(f[1], f[2]), radius, "{line}");
            assert_eq!(c.class(f[1], f[2]) as i32, class, "{line}");
            n += 1;
        }
        assert!(n > 1000);
    }
}
