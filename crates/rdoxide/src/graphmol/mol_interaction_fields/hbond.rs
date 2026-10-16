//! Hydrogen bonding energy between a probe and a molecule (similar to the
//! GRID hydrogen bonding descriptors; J. Med. Chem. 1989, 32, 1083 and 1993,
//! 36, 140 and 148) and the hydrophilic field built from it.

use std::sync::OnceLock;

use super::{CUTOFF, Descriptor, MIN_CUTOFF_VAL};
use crate::geometry::Point3D;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::conjug_hybrid::atom_has_conjugated_bond;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::{RDKitError, Result, rd_log_error};

const M_55D: f64 = 0.959_931_088_596_881_3;
/// angle C-O-H, C-N-H
const M_70_5D: f64 = 1.230_457_122_656_002_4;
const M_110D: f64 = 1.919_862_177_193_762_5;
const PI_2: f64 = std::f64::consts::FRAC_PI_2;

const EM: [[f64; 2]; 2] = [[-8.368, -11.715], [-11.715, -16.736]];
const RM: [[f64; 2]; 2] = [[3.2, 3.0], [3.0, 2.8]];
const K1: f64 = 562.25380293;
const K2: f64 = 0.11697778;
const BOND_LENGTH: [f64; 2] = [-0.972, -1.019];

/// The target atom type of an interaction (also the index into the tables).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    N = 0,
    O = 1,
}

/// The angular dependence of an interaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Ang {
    Cos2,
    Cos20,
    Cos2Rot,
    Cos4,
    Cos4Rot,
    Cos6,
    Cos6Rot,
    CosAcc,
    NoDep,
}

impl Ang {
    /// The function of the angle `t` (and the out-of-plane and in-plane
    /// angles `t0` and `ti` for the plane dependent ones).
    fn eval(self, t: f64, t0: f64, ti: f64) -> f64 {
        let cos_pow = |t: f64, n: i32| {
            if t < PI_2 {
                let c = t.cos();
                match n {
                    2 => c * c,
                    4 => {
                        let c2 = c * c;
                        c2 * c2
                    }
                    _ => {
                        let c2 = c * c;
                        c2 * (c2 * c2)
                    }
                }
            } else {
                0.0
            }
        };
        match self {
            Ang::Cos2 => cos_pow(t, 2),
            Ang::Cos20 => {
                if ti < M_55D {
                    let c = t0.cos();
                    c * c
                } else if ti < PI_2 {
                    let c = t.cos();
                    c * c
                } else {
                    0.0
                }
            }
            Ang::Cos2Rot => cos_pow(t - M_70_5D, 2),
            Ang::Cos4 => cos_pow(t, 4),
            Ang::Cos4Rot => cos_pow(t - M_70_5D, 4),
            Ang::Cos6 => cos_pow(t, 6),
            Ang::Cos6Rot => cos_pow(t - M_70_5D, 6),
            Ang::CosAcc => {
                if ti < PI_2 {
                    t0.cos() * (0.9 + 0.1 * (2.0 * ti).sin())
                } else if ti < M_110D {
                    let c = ti.cos();
                    let mut temp = K2 - c * c;
                    temp = temp * temp * temp;
                    temp * (t0.cos() * K1)
                } else {
                    0.0
                }
            }
            Ang::NoDep => 1.0,
        }
    }
}

fn normalized(mut p: Point3D) -> Result<Point3D> {
    p.normalize()?;
    Ok(p)
}

fn angle(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let dot = dot.clamp(-1.0, 1.0);
    dot.acos()
}

fn normalize3(v: &mut [f64; 3]) {
    let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    v[0] /= l;
    v[1] /= l;
    v[2] /= l;
}

/// `D` for a donor probe (the molecule accepts), `A` for an acceptor probe.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Probe {
    /// the probe acts as the acceptor ("OH"/"NH" are donating probes in
    /// GRID's naming; upstream calls them `A` and looks for the molecule's
    /// acceptors)
    A,
    D,
}

/// Hydrogen bonding descriptor.
#[derive(Debug, Clone)]
pub struct HBond {
    da: Probe,
    cutoff: f64,
    probe: Kind,
    n_interact: usize,
    target_types: Vec<Kind>,
    pos: Vec<f64>,
    direction: Vec<f64>,
    lengths: Vec<f64>,
    plane: Vec<f64>,
    function: Vec<Ang>,
}

/// The state shared by the search functions: the geometry of the molecule
/// and the values that, like upstream's locals, are kept from one atom to
/// the next.
struct Search<'a> {
    mol: &'a ROMol,
    conf: &'a Conformer,
    pos: Point3D,
    dir: Point3D,
    hbonddir: Point3D,
    plane: Point3D,
    bond_direction: [Point3D; 12],
}

impl<'a> Search<'a> {
    fn new(mol: &'a ROMol, conf: &'a Conformer) -> Self {
        Search {
            mol,
            conf,
            pos: Point3D::default(),
            dir: Point3D::default(),
            hbonddir: Point3D::default(),
            plane: Point3D::default(),
            bond_direction: [Point3D::default(); 12],
        }
    }

    fn at(&self, i: u32) -> Result<Point3D> {
        self.conf.atom_pos(i)
    }

    fn z(&self, i: u32) -> u32 {
        self.mol.atom(i).atomic_num()
    }

    fn neighbors(&self, i: u32) -> Vec<u32> {
        self.mol.atom_neighbors(i).collect()
    }

    /// The normalized vector from the neighbour to the atom.
    fn dir_from(&self, atom_pos: Point3D, nb: u32) -> Result<Point3D> {
        normalized(atom_pos - self.at(nb)?)
    }
}

impl HBond {
    /// `HBond(mol, confId, probeAtomType, fixed, cutoff)`: the probe type is
    /// "OH", "O", "NH" or "N".
    pub fn new(
        mol: &ROMol,
        conf_id: i32,
        probe_atom_type: &str,
        fixed: bool,
        cutoff: f64,
    ) -> Result<Self> {
        let mut cutoff = cutoff * cutoff;
        if cutoff < MIN_CUTOFF_VAL * MIN_CUTOFF_VAL {
            cutoff = CUTOFF * CUTOFF;
        }
        let (da, probe) = match probe_atom_type {
            "O" => (Probe::D, Kind::O),
            "OH" => (Probe::A, Kind::O),
            "N" => (Probe::D, Kind::N),
            "NH" => (Probe::A, Kind::N),
            _ => {
                let msg = format!("Probe atom type not supported: {probe_atom_type}");
                rd_log_error!("{msg}");
                return Err(RDKitError::Value(msg));
            }
        };
        let mut res = HBond {
            da,
            cutoff,
            probe,
            n_interact: mol.num_atoms() as usize,
            target_types: Vec::new(),
            pos: Vec::new(),
            direction: Vec::new(),
            lengths: Vec::new(),
            plane: Vec::new(),
            function: Vec::new(),
        };
        let conf = mol.conformer(conf_id)?.clone();
        let mut s = Search::new(mol, &conf);
        let specials = res.find_specials(&mut s, fixed)?;
        match (da, fixed) {
            (Probe::A, true) => res.find_acceptors(&mut s, &specials)?,
            (Probe::A, false) => res.find_acceptors_unfixed(&mut s, &specials)?,
            (Probe::D, true) => res.find_donors(&mut s, &specials)?,
            (Probe::D, false) => res.find_donors_unfixed(&mut s, &specials)?,
        }
        res.n_interact = res.target_types.len();
        Ok(res)
    }

    /// The number of interactions.
    pub fn num_interactions(&self) -> usize {
        self.n_interact
    }

    fn add(&mut self, kind: Kind, f: Ang, pos: Point3D, dir: Point3D, plane: Point3D) {
        self.target_types.push(kind);
        self.function.push(f);
        self.pos.extend([pos.x, pos.y, pos.z]);
        let len = dir.length();
        self.lengths.push(len);
        self.direction
            .extend([dir.x / len, dir.y / len, dir.z / len]);
        let len = plane.length();
        if len > 1.0e-16 {
            self.plane
                .extend([plane.x / len, plane.y / len, plane.z / len]);
        } else {
            self.plane.extend([0.0, 0.0, 0.0]);
        }
    }

    fn add_dir(&mut self, kind: Kind, f: Ang, pos: Point3D, dir: Point3D) {
        self.add(kind, f, pos, dir, Point3D::default());
    }

    /// `findSpecials`: the serine side chain, treated as a whole.
    fn find_specials(&mut self, s: &mut Search, fixed: bool) -> Result<Vec<u32>> {
        static SERINE: OnceLock<ROMol> = OnceLock::new();
        let ser = SERINE.get_or_init(|| {
            mol_from_smiles("[CH2]([C@@H](C(=O))N)O", &SmilesParserParams::default())
                .ok()
                .flatten()
                .expect("the serine pattern is valid")
        });
        let params = SubstructMatchParameters {
            max_matches: 1,
            ..Default::default()
        };
        let matches = substruct_match(s.mol, ser, &params)?;
        let mut specials = Vec::new();
        for &(_, idx) in matches.first().map_or(&[][..], Vec::as_slice) {
            let idx = idx as u32;
            if s.z(idx) != 8 {
                continue;
            }
            let nbrs = s.neighbors(idx);
            s.pos = s.at(idx)?;
            let pos = s.pos;
            match self.da {
                Probe::A => {
                    let mut n = 0;
                    for &nb in &nbrs {
                        if s.z(nb) == 1 {
                            s.hbonddir = s.dir_from(pos, nb)?;
                        } else {
                            s.bond_direction[n] = s.dir_from(pos, nb)?;
                            n += 1;
                        }
                    }
                    if n == 2 {
                        s.dir = s.bond_direction[0] + s.hbonddir;
                        s.plane = s.bond_direction[0].cross_product(&s.hbonddir);
                        s.plane = s.plane.cross_product(&s.dir);
                        if fixed {
                            self.add(Kind::O, Ang::Cos20, pos, s.dir, s.plane);
                        } else {
                            self.add_dir(Kind::O, Ang::Cos4Rot, pos, s.bond_direction[0]);
                        }
                        specials.push(idx);
                    }
                }
                Probe::D => {
                    if fixed {
                        for &nb in &nbrs {
                            if s.z(nb) == 1 {
                                s.dir = s.at(nb)? - pos;
                                self.add_dir(Kind::O, Ang::Cos6, pos, s.dir);
                            }
                        }
                    } else {
                        for &nb in &nbrs {
                            if s.z(nb) != 1 {
                                s.dir = pos - s.at(nb)?;
                                self.add_dir(Kind::O, Ang::Cos6Rot, pos, s.dir);
                            }
                        }
                    }
                    specials.push(idx);
                }
            }
        }
        Ok(specials)
    }

    /// The acceptors with fixed hydrogens (`findAcceptors`).
    fn find_acceptors(&mut self, s: &mut Search, specials: &[u32]) -> Result<()> {
        for i in 0..self.n_interact as u32 {
            if specials.contains(&i) {
                continue;
            }
            let atom = s.mol.atom(i);
            let charge = atom.formal_charge();
            match atom.atomic_num() {
                7 => {
                    let nbrs = s.neighbors(i);
                    s.pos = s.at(i)?;
                    let pos = s.pos;
                    for (k, &nb) in nbrs.iter().enumerate() {
                        s.bond_direction[k] = s.dir_from(pos, nb)?;
                    }
                    let bd = s.bond_direction;
                    match nbrs.len() {
                        1 => {
                            if charge <= 0 {
                                self.add_dir(Kind::N, Ang::Cos2, pos, bd[0]);
                            }
                        }
                        2 => {
                            if charge <= 0 {
                                s.dir = bd[0] + bd[1];
                                s.plane = bd[0].cross_product(&bd[1]);
                                self.add(Kind::N, Ang::Cos2, pos, s.dir, s.plane);
                            }
                        }
                        3 => {
                            if charge <= 0 && !atom_has_conjugated_bond(s.mol, i) {
                                s.dir = bd[0] + bd[1] + bd[2];
                                self.add_dir(Kind::N, Ang::Cos2, pos, s.dir);
                            }
                        }
                        0 | 4 => {}
                        _ => rd_log_error!(
                            "HBond: Nitrogen atom bound to more than 4 neighbor atoms: Atom: {i}"
                        ),
                    }
                }
                8 => {
                    let nbrs = s.neighbors(i);
                    s.pos = s.at(i)?;
                    let pos = s.pos;
                    let mut aromatic = false;
                    for (k, &nb) in nbrs.iter().enumerate() {
                        s.bond_direction[k] = s.dir_from(pos, nb)?;
                        aromatic |= s.mol.atom(nb).is_aromatic();
                    }
                    let bd = s.bond_direction;
                    match nbrs.len() {
                        1 => {
                            self.carbonyl_like(s, i, nbrs[0], charge, bd[0])?;
                        }
                        2 => {
                            s.dir = bd[0] + bd[1];
                            s.plane = bd[0].cross_product(&bd[1]);
                            s.plane = s.plane.cross_product(&s.dir);
                            let f = if aromatic { Ang::Cos2 } else { Ang::Cos20 };
                            self.add(Kind::O, f, pos, s.dir, s.plane);
                        }
                        0 | 3 => {}
                        _ => rd_log_error!(
                            "HBond: Oxygen atom bound to more than 3 neighbor atoms: Atom: {i}"
                        ),
                    }
                }
                9 | 17 | 35 | 53 => self.halogen(s, i)?,
                _ => {}
            }
        }
        Ok(())
    }

    /// Oxygens with one neighbour: carbonyl, X=O, anions.
    fn carbonyl_like(
        &mut self,
        s: &mut Search,
        i: u32,
        nb: u32,
        charge: i32,
        bd0: Point3D,
    ) -> Result<()> {
        let pos = s.pos;
        for sec in s.neighbors(nb) {
            if sec != i {
                // the bond direction of the neighbour to its other neighbour
                s.dir = s.at(sec)? - s.at(nb)?;
                break;
            }
        }
        s.plane = bd0.cross_product(&s.dir);
        if charge == 0 {
            self.add(Kind::O, Ang::CosAcc, pos, bd0, s.plane);
        } else if charge == -1 {
            let f = if atom_has_conjugated_bond(s.mol, i) {
                Ang::CosAcc
            } else {
                Ang::Cos2
            };
            self.add(Kind::O, f, pos, bd0, s.plane);
        }
        Ok(())
    }

    fn halogen(&mut self, s: &mut Search, i: u32) -> Result<()> {
        s.pos = s.at(i)?;
        s.dir = Point3D::new(1.0, 1.0, 1.0);
        // halogens are like nitrogens; no direction needed
        self.add_dir(Kind::N, Ang::NoDep, s.pos, s.dir);
        Ok(())
    }

    /// The acceptors that may rotate (`findAcceptorsUnfixed`).
    fn find_acceptors_unfixed(&mut self, s: &mut Search, specials: &[u32]) -> Result<()> {
        let zero = Point3D::default();
        for i in 0..self.n_interact as u32 {
            if specials.contains(&i) {
                continue;
            }
            let atom = s.mol.atom(i);
            let charge = atom.formal_charge();
            match atom.atomic_num() {
                7 => {
                    let nbrs = s.neighbors(i);
                    s.pos = s.at(i)?;
                    let pos = s.pos;
                    let mut nonh = 0;
                    for &nb in &nbrs {
                        if s.z(nb) != 1 {
                            s.bond_direction[nonh] = s.dir_from(pos, nb)?;
                            nonh += 1;
                        }
                    }
                    // the directions to the hydrogens among the last `k` neighbours
                    let add_h_dirs = |s: &mut Search, k: usize, mut n: usize| -> Result<()> {
                        for &nb in &nbrs[nbrs.len() - k..] {
                            if s.z(nb) == 1 {
                                s.bond_direction[n] = s.dir_from(pos, nb)?;
                                n += 1;
                            }
                        }
                        Ok(())
                    };
                    match nbrs.len() {
                        1 => {
                            if charge <= 0 {
                                self.add_dir(Kind::N, Ang::Cos2, pos, s.bond_direction[0]);
                            }
                        }
                        2 => {
                            if charge <= 0 {
                                if nonh == 2 {
                                    let bd = s.bond_direction;
                                    s.dir = bd[0] + bd[1];
                                    s.plane = bd[0].cross_product(&bd[1]);
                                    self.add(Kind::N, Ang::Cos2, pos, s.dir, s.plane);
                                } else if nonh == 1 {
                                    add_h_dirs(s, 2, nonh)?;
                                    let bd = s.bond_direction;
                                    s.plane = bd[0].cross_product(&bd[1]);
                                    self.add(Kind::N, Ang::CosAcc, pos, bd[0], s.plane);
                                } else {
                                    self.add_dir(Kind::N, Ang::NoDep, pos, zero);
                                }
                            }
                        }
                        3 => {
                            if charge <= 0 {
                                if nonh == 0 {
                                    self.add_dir(Kind::N, Ang::NoDep, pos, zero);
                                } else if nonh == 1 {
                                    self.add_dir(Kind::N, Ang::Cos2Rot, pos, s.bond_direction[0]);
                                } else if !atom_has_conjugated_bond(s.mol, i) {
                                    add_h_dirs(s, 3, nonh)?;
                                    let bd = s.bond_direction;
                                    s.dir = bd[0] + bd[1] + bd[2];
                                    self.add_dir(Kind::N, Ang::Cos2, pos, s.dir);
                                }
                            }
                        }
                        0 | 4 => {}
                        _ => rd_log_error!(
                            "HBond: Nitrogen atom bound to more than 4 neighbor atoms: Atom: {i}"
                        ),
                    }
                }
                8 => {
                    let nbrs = s.neighbors(i);
                    s.pos = s.at(i)?;
                    let pos = s.pos;
                    let mut nonh = 0;
                    let mut aromatic = false;
                    for &nb in &nbrs {
                        if s.z(nb) != 1 {
                            s.bond_direction[nonh] = s.dir_from(pos, nb)?;
                            nonh += 1;
                        }
                        aromatic |= s.mol.atom(nb).is_aromatic();
                    }
                    match nbrs.len() {
                        1 => {
                            let bd0 = s.bond_direction[0];
                            self.carbonyl_like(s, i, nbrs[0], charge, bd0)?;
                        }
                        2 => {
                            if nonh == 0 {
                                self.add_dir(Kind::O, Ang::NoDep, pos, zero);
                            } else if nonh == 1 {
                                if aromatic {
                                    let mut n = nonh;
                                    for &nb in &nbrs {
                                        if s.z(nb) == 1 {
                                            s.bond_direction[n] = s.dir_from(pos, nb)?;
                                            n += 1;
                                        }
                                    }
                                    let bd = s.bond_direction;
                                    s.plane = bd[0].cross_product(&bd[1]);
                                    self.add(Kind::O, Ang::CosAcc, pos, bd[0], s.plane);
                                } else {
                                    self.add_dir(Kind::O, Ang::Cos2Rot, pos, s.bond_direction[0]);
                                }
                            } else {
                                let bd = s.bond_direction;
                                s.dir = bd[0] + bd[1];
                                s.plane = bd[0].cross_product(&bd[1]);
                                s.plane = s.plane.cross_product(&s.dir);
                                self.add(Kind::O, Ang::Cos20, pos, s.dir, s.plane);
                            }
                        }
                        0 | 3 => {}
                        _ => rd_log_error!(
                            "HBond: Oxygen atom bound to more than 3 neighbor atoms: Atom: {i}"
                        ),
                    }
                }
                9 | 17 | 35 | 53 => self.halogen(s, i)?,
                _ => {}
            }
        }
        Ok(())
    }

    /// The donors with fixed hydrogens (`findDonors`).
    fn find_donors(&mut self, s: &mut Search, specials: &[u32]) -> Result<()> {
        for i in 0..self.n_interact as u32 {
            if specials.contains(&i) {
                continue;
            }
            let (kind, f) = match s.z(i) {
                7 => (Kind::N, Ang::Cos2),
                8 => (Kind::O, Ang::Cos4),
                _ => continue,
            };
            s.pos = s.at(i)?;
            for nb in s.neighbors(i) {
                if s.z(nb) == 1 {
                    // not normalized: the length is needed for donors
                    s.dir = s.at(nb)? - s.pos;
                    self.add_dir(kind, f, s.pos, s.dir);
                }
            }
        }
        Ok(())
    }

    /// Swaps the hydrogen and the lone pair (phenolic oxygens and imines).
    fn swap_h(s: &mut Search) {
        s.plane = s.dir.cross_product(&s.hbonddir);
        s.plane = s.plane.cross_product(&s.dir);
        s.dir = s.plane * s.hbonddir.dot_product(&s.plane);
        s.hbonddir -= s.dir * 2.0;
    }

    /// The donors that may rotate (`findDonorsUnfixed`).
    fn find_donors_unfixed(&mut self, s: &mut Search, specials: &[u32]) -> Result<()> {
        let zero = Point3D::default();
        for i in 0..self.n_interact as u32 {
            if specials.contains(&i) {
                continue;
            }
            match s.z(i) {
                7 => {
                    let nbrs = s.neighbors(i);
                    s.pos = s.at(i)?;
                    let pos = s.pos;
                    let nonh = nbrs.iter().filter(|&&nb| s.z(nb) != 1).count();
                    if nonh == nbrs.len() {
                        continue; // no hydrogens, no donation
                    }
                    let n_h = BOND_LENGTH[Kind::N as usize];
                    match nbrs.len() {
                        2 => {
                            if nonh == 0 {
                                self.add_dir(Kind::N, Ang::NoDep, pos, zero);
                            } else {
                                for &nb in &nbrs {
                                    if s.z(nb) != 1 {
                                        s.dir = normalized(pos - s.at(nb)?)?;
                                    } else {
                                        s.hbonddir = s.at(nb)? - pos;
                                    }
                                }
                                self.add_dir(Kind::N, Ang::Cos2, pos, s.hbonddir);
                                Self::swap_h(s);
                                self.add_dir(Kind::N, Ang::Cos2, pos, s.hbonddir);
                            }
                        }
                        3 => {
                            if nonh == 2 {
                                for &nb in &nbrs {
                                    if s.z(nb) == 1 {
                                        s.hbonddir = s.at(nb)? - pos;
                                    }
                                }
                                self.add_dir(Kind::N, Ang::Cos2, pos, s.hbonddir);
                            } else if nonh == 1 {
                                for &nb in &nbrs {
                                    if s.z(nb) != 1 {
                                        s.dir = normalized(pos - s.at(nb)?)?;
                                    }
                                }
                                self.add_dir(Kind::N, Ang::Cos2Rot, pos, s.dir * (-n_h));
                            } else {
                                self.add_dir(Kind::N, Ang::NoDep, pos, zero);
                            }
                        }
                        4 => {
                            if nonh == 0 {
                                self.add_dir(Kind::N, Ang::NoDep, pos, zero);
                            } else if nonh == 1 {
                                for &nb in &nbrs {
                                    if s.z(nb) != 1 {
                                        s.dir = pos - s.at(nb)?;
                                    }
                                }
                                self.add_dir(Kind::N, Ang::Cos2Rot, pos, s.dir * (-n_h));
                            } else {
                                for &nb in &nbrs {
                                    if s.z(nb) == 1 {
                                        s.dir = pos - s.at(nb)?;
                                        self.add_dir(Kind::N, Ang::Cos2Rot, pos, s.dir);
                                    }
                                }
                            }
                        }
                        _ => rd_log_error!(
                            "HBond: Nitrogen atom bound to more than 4 neighbor atoms: Atom: {i}"
                        ),
                    }
                }
                8 => {
                    let nbrs = s.neighbors(i);
                    s.pos = s.at(i)?;
                    let pos = s.pos;
                    let nonh = nbrs.iter().filter(|&&nb| s.z(nb) != 1).count();
                    let aromatic = nbrs.iter().any(|&nb| s.mol.atom(nb).is_aromatic());
                    if nonh == nbrs.len() {
                        continue;
                    }
                    let o_h = BOND_LENGTH[Kind::O as usize];
                    match nbrs.len() {
                        1 => self.add_dir(Kind::O, Ang::NoDep, pos, zero),
                        2 => {
                            if nonh == 0 {
                                self.add_dir(Kind::O, Ang::NoDep, pos, zero);
                            } else {
                                for &nb in &nbrs {
                                    if s.z(nb) != 1 {
                                        s.dir = normalized(pos - s.at(nb)?)?;
                                    } else {
                                        s.hbonddir = s.at(nb)? - pos;
                                    }
                                }
                                if aromatic {
                                    self.add_dir(Kind::O, Ang::Cos4, pos, s.hbonddir);
                                    Self::swap_h(s);
                                    self.add_dir(Kind::O, Ang::Cos4, pos, s.hbonddir);
                                } else {
                                    self.add_dir(Kind::O, Ang::Cos4Rot, pos, s.dir * (-o_h));
                                }
                            }
                        }
                        3 => {
                            if nonh == 0 {
                                self.add_dir(Kind::O, Ang::NoDep, pos, zero);
                            } else if nonh == 1 {
                                for &nb in &nbrs {
                                    if s.z(nb) != 1 {
                                        s.dir = normalized(pos - s.at(nb)?)?;
                                        self.add_dir(Kind::O, Ang::Cos4Rot, pos, s.dir * (-o_h));
                                    }
                                }
                            } else {
                                for &nb in &nbrs {
                                    if s.z(nb) == 1 {
                                        s.dir = s.at(nb)? - pos;
                                        self.add_dir(Kind::O, Ang::Cos4, pos, s.dir);
                                    }
                                }
                            }
                        }
                        _ => rd_log_error!(
                            "HBond: Oxygen atom bound to more than 3 neighbor atoms: Atom: {i}"
                        ),
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }
}

impl Descriptor for HBond {
    fn value(&self, x: f64, y: f64, z: f64, thres: f64) -> f64 {
        let n = self.n_interact;
        if n < 1 {
            return 0.0;
        }
        let mut vtp = vec![0.0; 3 * n]; // target -> probe vectors
        let mut ene = vec![0.0; n];
        let mut min_id = 0;
        let mut min_ene = f64::MAX;
        for i in 0..n {
            let j = 3 * i;
            vtp[j] = x - self.pos[j];
            vtp[j + 1] = y - self.pos[j + 1];
            vtp[j + 2] = z - self.pos[j + 2];
            let mut dist2 = vtp[j] * vtp[j] + vtp[j + 1] * vtp[j + 1] + vtp[j + 2] * vtp[j + 2];
            if dist2 < thres {
                let dis = dist2.sqrt();
                let dist_n = [vtp[j] / dis, vtp[j + 1] / dis, vtp[j + 2] / dis];
                dist2 = dist2.max(self.cutoff);
                let dir = [
                    self.direction[j],
                    self.direction[j + 1],
                    self.direction[j + 2],
                ];
                let t = angle(dist_n, dir);
                let rm = RM[self.probe as usize][self.target_types[i] as usize];
                let mut term1 = rm * rm;
                term1 /= dist2;
                let mut term2 = term1 * term1;
                term1 *= term2;
                term2 *= term2;
                term2 *= 4.0 * term1 - 3.0 * term2;
                term2 *= EM[self.probe as usize][self.target_types[i] as usize];
                ene[i] = term2;

                let (t0, ti) = if self.da == Probe::A
                    && matches!(self.function[i], Ang::CosAcc | Ang::Cos20)
                {
                    let plane = [self.plane[j], self.plane[j + 1], self.plane[j + 2]];
                    let dot = vtp[j] * plane[0] + vtp[j + 1] * plane[1] + vtp[j + 2] * plane[2];
                    let mut in_plane = [
                        vtp[j] - plane[0] * dot,
                        vtp[j + 1] - plane[1] * dot,
                        vtp[j + 2] - plane[2] * dot,
                    ];
                    normalize3(&mut in_plane);
                    (angle(dist_n, in_plane), angle(dir, in_plane))
                } else {
                    (0.0, 1.0)
                };
                ene[i] *= self.function[i].eval(t, t0, ti);
                if ene[i] < min_ene {
                    min_ene = ene[i];
                    min_id = i;
                }
            } else {
                ene[i] = 0.0;
            }
        }
        let m = 3 * min_id;
        let mut probe_dir = [-vtp[m], -vtp[m + 1], -vtp[m + 2]];
        normalize3(&mut probe_dir);

        let mut res = 0.0;
        for i in 0..n {
            let j = 3 * i;
            let mut v = match self.da {
                Probe::A => {
                    let bl = BOND_LENGTH[self.probe as usize];
                    [
                        probe_dir[0] * bl - vtp[j],
                        probe_dir[1] * bl - vtp[j + 1],
                        probe_dir[2] * bl - vtp[j + 2],
                    ]
                }
                Probe::D => [
                    self.direction[j] * self.lengths[i] - vtp[j],
                    self.direction[j + 1] * self.lengths[i] - vtp[j + 1],
                    self.direction[j + 2] * self.lengths[i] - vtp[j + 2],
                ],
            };
            normalize3(&mut v);
            // the angle between the best hydrogen bond direction of the
            // probe and the hydrogen-acceptor vector
            let p = match self.da {
                Probe::A => angle(probe_dir, v),
                Probe::D => angle(v, probe_dir),
            };
            res += ene[i] * Ang::Cos2.eval(p, 0.0, 1.0);
        }
        res
    }
}

/// The hydrophilic potential: the more favourable of the hydrogen bond
/// interaction with a water oxygen ("O" probe) and a water hydrogen ("OH"
/// probe).
#[derive(Debug, Clone)]
pub struct Hydrophilic {
    oh: HBond,
    o: HBond,
}

impl Hydrophilic {
    pub fn new(mol: &ROMol, conf_id: i32, fixed: bool, cutoff: f64) -> Result<Self> {
        Ok(Hydrophilic {
            oh: HBond::new(mol, conf_id, "OH", fixed, cutoff)?,
            o: HBond::new(mol, conf_id, "O", fixed, cutoff)?,
        })
    }
}

impl Descriptor for Hydrophilic {
    fn value(&self, x: f64, y: f64, z: f64, thres: f64) -> f64 {
        let o = self.o.value(x, y, z, thres);
        let oh = self.oh.value(x, y, z, thres);
        o.min(oh)
    }
}
