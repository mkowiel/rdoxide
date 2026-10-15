//! Resonance structures (upstream `GraphMol/Resonance.{h,cpp}`,
//! `ResonanceMolSupplier`, Copyright (C) 2015 Paolo Tosco).
//!
//! Every conjugated group of the Kekulized molecule is enumerated on its own:
//! the multiple-bond arrangements (`ConjElectrons`) are generated depth first,
//! the non-bonded electrons are distributed, implausible charge patterns are
//! pruned and the survivors are ranked. The complete resonance structures
//! combine the ranked arrangements of the groups breadth first.
//!
//! Deliberate differences from upstream: duplicates are detected on the exact
//! fingerprint / metrics instead of their `boost::hash_range` hash, and
//! structures that compare equal in every ranking criterion keep the order
//! in which they were generated (upstream sorts them in the iteration order
//! of an `std::unordered_map`, which differs between standard libraries).

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use crate::error::{RDKitError, Result};
use crate::graphmol::bond::BondType;
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::matrices::get_distance_mat;
use crate::graphmol::mol_ops::sanitize_flags::{SANITIZE_ADJUSTHS, SANITIZE_FINDRADICALS};
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::{MatchVect, SubstructMatchParameters, substruct_match};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::threads::{num_threads_to_use, par_map};

/// `ResonanceMolSupplier::ResonanceFlags`
pub const ALLOW_INCOMPLETE_OCTETS: u32 = 1 << 0;
pub const ALLOW_CHARGE_SEPARATION: u32 = 1 << 1;
pub const KEKULE_ALL: u32 = 1 << 2;
pub const UNCONSTRAINED_CATIONS: u32 = 1 << 3;
pub const UNCONSTRAINED_ANIONS: u32 = 1 << 4;

const MAX_STRUCTS: u32 = 1_000_000;

/// The state a progress callback can inspect
/// (`ResonanceMolSupplierCallback::getNumStructures` and friends).
#[derive(Debug, Clone, Default)]
pub struct ResonanceProgress {
    pub max_structs: usize,
    /// per conjugated group: (structures so far, non-degenerate structures)
    pub groups: Vec<(usize, usize)>,
}

/// Called at each iteration of the enumeration; returning `false` stops it.
pub type ProgressCallback = Box<dyn FnMut(&ResonanceProgress) -> bool + Send + Sync>;

// ---------------------------------------------------------------- electrons

const LAST_BOND: u8 = 1 << 0;
const DEFINITIVE: u8 = 1 << 1;
const STACKED: u8 = 1 << 2;
const HAS_MULTIPLE_BOND: u8 = 1 << 3;

const HAVE_CATION_RIGHT_OF_N: u8 = 1 << 0;
const HAVE_CATION: u8 = 1 << 1;
const HAVE_ANION: u8 = 1 << 2;

const FP_BONDS: u32 = 1 << 0;
const FP_ATOMS: u32 = 1 << 1;

const NEED_CHARGE_BIT: u8 = 1;

/// 1000 * Allen electronegativity (Allen, L.C. J. Am. Chem. Soc. 1989, 111,
/// 9003-9014)
const ELECTRONEGATIVITY: [i32; 89] = [
    1000, 2300, 4160, 912, 1576, 2051, 2544, 3066, 3610, 4193, 4789, 869, 1293, 1613, 1916, 2253,
    2589, 2869, 3242, 734, 1034, 1190, 1380, 1530, 1650, 1750, 1800, 1840, 1880, 1850, 1590, 1756,
    1994, 2211, 2434, 2685, 2966, 706, 963, 1120, 1320, 1410, 1470, 1510, 1540, 1560, 1590, 1870,
    1520, 1656, 1824, 1984, 2158, 2359, 2582, 659, 881, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
    1000, 1000, 1000, 1000, 1000, 1000, 1000, 1090, 1160, 1340, 1470, 1600, 1650, 1680, 1720, 1920,
    1760, 1789, 1854, 2010, 2190, 2390, 2600, 670, 890,
];

/// The read-only data shared by all conjugated groups.
struct Ctx<'a> {
    mol: &'a ROMol,
    flags: u32,
    max_structs: usize,
    atom_grp: &'a [i32],
    bond_grp: &'a [i32],
    dist: &'a [f64],
}

#[derive(Clone)]
struct AtomElectrons {
    nb: u8,
    tv: u8,
    fc: i8,
    flags: u8,
    atomic_num: u32,
    aromatic: bool,
}

impl AtomElectrons {
    fn oe(&self) -> u32 {
        PeriodicTable::get()
            .get_n_outer_elecs(self.atomic_num)
            .unwrap_or(0) as u32
    }

    fn has_octet(&self) -> bool {
        u32::from(self.nb) + u32::from(self.tv) * 2 == 8
    }

    fn is_definitive(&self) -> bool {
        self.flags & DEFINITIVE != 0
    }

    fn is_last_bond(&self) -> bool {
        self.flags & LAST_BOND != 0
    }

    /// `8 - (2 * tv + nb)` as an unsigned integer (wrapping like upstream)
    fn needed_nb_for_octet(&self) -> u32 {
        (8 - (2 * i32::from(self.tv) + i32::from(self.nb))) as u32
    }

    fn assign_formal_charge(&mut self) {
        self.fc = (self.oe() as i32 - (i32::from(self.nb) + i32::from(self.tv))) as i8;
    }
}

#[derive(Clone)]
struct BondElectrons {
    bo: u8,
    flags: u8,
    begin: u32,
    end: u32,
}

impl BondElectrons {
    fn is_definitive(&self) -> bool {
        self.flags & DEFINITIVE != 0
    }
}

fn order_from_bond(mol: &ROMol, bi: u32) -> Result<u8> {
    let bond = mol.bond(bi);
    match bond.bond_type() {
        BondType::Single => Ok(1),
        BondType::Double => Ok(2),
        BondType::Triple => Ok(3),
        _ => Err(RDKitError::Runtime(format!(
            "Bond idx {bi} between atoms {} and {} has an invalid bond type",
            bond.begin_atom_idx(),
            bond.end_atom_idx()
        ))),
    }
}

#[derive(Clone, Default)]
struct CeMetrics {
    abs_formal_charges: i32,
    fc_same_sign_dist: i32,
    fc_opp_sign_dist: i32,
    nb_missing: i32,
    wtd_formal_charges: i32,
    sum_formal_charge_idxs: i32,
    sum_multiple_bond_idxs: i32,
    /// stands in for upstream's hash of the first five metrics; `None`
    /// until the structure is counted
    key: Option<[i32; 5]>,
}

#[derive(Default)]
struct CeStats {
    init: bool,
    have_no_cations_right_of_n: bool,
    have_no_charge_separation: bool,
    min_nb_missing: u32,
}

/// `ConjElectrons`: one arrangement of bond orders, non-bonded electrons and
/// formal charges of a conjugated group.
#[derive(Clone)]
struct ConjElectrons {
    total_electrons: u32,
    curr_electrons: u32,
    total_formal_charge: i32,
    allowed_chg_left_of_n: i32,
    flags: u8,
    metrics: CeMetrics,
    bonds: BTreeMap<u32, BondElectrons>,
    atoms: BTreeMap<u32, AtomElectrons>,
    begin_stack: Vec<u32>,
}

/// The structures of one group by their fingerprint, in insertion order.
#[derive(Default)]
struct CeMap {
    slots: Vec<Option<(Vec<u8>, ConjElectrons)>>,
    index: HashMap<Vec<u8>, usize>,
}

impl CeMap {
    fn len(&self) -> usize {
        self.index.len()
    }

    fn contains(&self, fp: &[u8]) -> bool {
        self.index.contains_key(fp)
    }

    fn insert(&mut self, fp: Vec<u8>, ce: ConjElectrons) {
        self.index.insert(fp.clone(), self.slots.len());
        self.slots.push(Some((fp, ce)));
    }

    fn get_mut(&mut self, fp: &[u8]) -> Option<&mut ConjElectrons> {
        let i = *self.index.get(fp)?;
        self.slots[i].as_mut().map(|(_, ce)| ce)
    }

    fn into_structures(self) -> Vec<ConjElectrons> {
        self.slots.into_iter().flatten().map(|(_, ce)| ce).collect()
    }
}

type CeDegCount = HashMap<[i32; 5], u32>;

/// 2 bits per bond order, 4 bits of total valence and 4 of non-bonded
/// electrons per atom
fn pack_fp(ce: &ConjElectrons, flags: u32) -> Vec<u8> {
    let mut fp = Vec::new();
    if flags & FP_ATOMS != 0 {
        for ae in ce.atoms.values() {
            fp.push((u32::from(ae.tv) | (u32::from(ae.nb) << 4)) as u8);
        }
    }
    if flags & FP_BONDS != 0 {
        let mut i = 0;
        let mut byte = 0u8;
        for be in ce.bonds.values() {
            if i != 0 && i % 4 == 0 {
                fp.push(byte);
                byte = 0;
                i = 0;
            }
            byte |= be.bo << (i * 2);
            i += 1;
        }
        if i != 0 {
            fp.push(byte);
        }
    }
    fp
}

/// number of permutations and the first binary code of `n_cand` unsatisfied
/// atoms out of `n_total`
fn num_comb_start_v(n_cand: u32, n_total: u32) -> (u32, u32) {
    let mut num_comb = 1u32;
    let mut v = 0u32;
    for i in 0..n_cand {
        num_comb = num_comb.wrapping_mul(n_total.wrapping_sub(i));
        num_comb /= i + 1;
        v |= 1u32.wrapping_shl(i);
    }
    (num_comb, v)
}

/// the next permutation (Gosper's hack)
fn update_v(v: u32) -> u32 {
    let t = (v | v.wrapping_sub(1)).wrapping_add(1);
    let low_t = t & (!t).wrapping_add(1);
    let low_v = v & (!v).wrapping_add(1);
    if low_v == 0 {
        return 0;
    }
    t | ((low_t / low_v) >> 1).wrapping_sub(1)
}

impl ConjElectrons {
    fn new(ctx: &Ctx, group_idx: u32) -> Result<Self> {
        let mol = ctx.mol;
        let mut ce = ConjElectrons {
            total_electrons: 0,
            curr_electrons: 0,
            total_formal_charge: 0,
            allowed_chg_left_of_n: 0,
            flags: 0,
            metrics: CeMetrics::default(),
            bonds: BTreeMap::new(),
            atoms: BTreeMap::new(),
            begin_stack: Vec::new(),
        };
        for ai in 0..mol.num_atoms() {
            if ctx.atom_grp[ai as usize] != -1 {
                ce.total_formal_charge += mol.atom(ai).formal_charge();
            }
        }
        ce.allowed_chg_left_of_n = ce.total_formal_charge;
        for bi in 0..mol.num_bonds() {
            if ctx.bond_grp[bi as usize] != group_idx as i32 {
                continue;
            }
            let bond = mol.bond(bi);
            ce.bonds.insert(
                bi,
                BondElectrons {
                    bo: 1,
                    flags: 0,
                    begin: bond.begin_atom_idx(),
                    end: bond.end_atom_idx(),
                },
            );
            for ai in [bond.begin_atom_idx(), bond.end_atom_idx()] {
                if !ce.atoms.contains_key(&ai) {
                    let atom = mol.atom(ai);
                    let mut tv = mol.total_degree(ai)? as u8;
                    for b in mol.atom_bonds(ai) {
                        if tv != 0 && b.bond_type() >= BondType::DativeOne {
                            tv -= 1;
                        }
                    }
                    ce.atoms.insert(
                        ai,
                        AtomElectrons {
                            nb: 0,
                            tv,
                            fc: 0,
                            flags: 0,
                            atomic_num: atom.atomic_num(),
                            aromatic: atom.is_aromatic(),
                        },
                    );
                }
            }
        }
        // total number of valence electrons in the conjugated group
        let mut total = 0u32;
        for &bi in ce.bonds.keys() {
            total = total.wrapping_add(2 * u32::from(order_from_bond(mol, bi)?));
        }
        for (&ai, ae) in &ce.atoms {
            let atom = mol.atom(ai);
            total = total
                .wrapping_add(ae.oe())
                .wrapping_sub(atom.total_valence()?)
                .wrapping_sub(atom.formal_charge() as u32);
        }
        ce.total_electrons = total;
        ce.curr_electrons = total;
        Ok(ce)
    }

    fn decr_curr_electrons(&mut self, d: u32) -> Result<()> {
        if self.curr_electrons < d {
            return Err(RDKitError::Runtime(format!(
                "d_currElectrons = {}, d = {d}",
                self.curr_electrons
            )));
        }
        self.curr_electrons -= d;
        Ok(())
    }

    fn push_to_begin_stack(&mut self, ai: u32) {
        let ae = self.atoms.get_mut(&ai).expect("atom of the group");
        if ae.flags & STACKED == 0 {
            ae.flags |= STACKED;
            self.begin_stack.push(ai);
        }
    }

    fn pop_from_begin_stack(&mut self) -> Option<u32> {
        while let Some(ai) = self.begin_stack.pop() {
            let ae = self.atoms.get_mut(&ai).expect("atom of the group");
            ae.flags &= !STACKED;
            if !ae.is_definitive() {
                return Some(ai);
            }
        }
        None
    }

    /// The user-supplied bond orders and charges.
    fn init_from_mol(&mut self, ctx: &Ctx) -> Result<()> {
        for (&ai, ae) in &mut self.atoms {
            let atom = ctx.mol.atom(ai);
            ae.tv = atom.total_valence()? as u8;
            ae.fc = atom.formal_charge() as i8;
            ae.nb = (ae.oe() as i32 - i32::from(ae.tv) - i32::from(ae.fc)) as u8;
        }
        for (&bi, be) in &mut self.bonds {
            be.bo = order_from_bond(ctx.mol, bi)?;
        }
        self.curr_electrons = 0;
        Ok(())
    }

    fn assign_non_bonded(&mut self) {
        let mut curr = self.curr_electrons;
        for ae in self.atoms.values_mut() {
            let nb = ae.needed_nb_for_octet().min(curr);
            curr -= nb;
            ae.nb = nb as u8;
        }
        self.curr_electrons = curr;
    }

    fn assign_formal_charge(&mut self) {
        for ae in self.atoms.values_mut() {
            ae.assign_formal_charge();
        }
    }

    fn set_order(&mut self, bi: u32, bo: u8) -> Result<()> {
        let (begin, end) = {
            let be = &self.bonds[&bi];
            (be.begin, be.end)
        };
        for ai in [begin, end] {
            let ae = self.atoms.get_mut(&ai).expect("atom of the group");
            ae.tv = ae.tv.wrapping_add(bo.wrapping_sub(1));
        }
        let be = self.bonds.get_mut(&bi).expect("bond of the group");
        be.flags |= DEFINITIVE;
        self.decr_curr_electrons(u32::from(bo) * 2)?;
        self.bonds.get_mut(&bi).expect("bond of the group").bo = bo;
        Ok(())
    }

    /// Is any conjugated neighbor of atom `ai` charged (and with outer
    /// electrons `oe_constraint`, if that is not zero)?
    fn is_nbr_charged(&self, ctx: &Ctx, ai: u32, bo: u32, oe_constraint: u32) -> bool {
        let grp = ctx.atom_grp[ai as usize];
        for bond in ctx.mol.atom_bonds(ai) {
            let bi_nbr = bond.idx();
            if ctx.bond_grp[bi_nbr as usize] != grp {
                continue;
            }
            let Some(be_nbr) = self.bonds.get(&bi_nbr) else {
                continue;
            };
            let Ok(ai_nbr) = bond.other_atom_idx(ai) else {
                continue;
            };
            let Some(ae_nbr) = self.atoms.get(&ai_nbr) else {
                continue;
            };
            let res = ((be_nbr.is_definitive() && !ae_nbr.has_octet())
                || (!be_nbr.is_definitive() && ae_nbr.is_definitive() && ae_nbr.oe() < 5 - bo))
                && (oe_constraint == 0 || ae_nbr.oe() == oe_constraint);
            if res {
                return true;
            }
        }
        false
    }

    /// Marks atom `ai` as having its last non-definitive bond when every
    /// other conjugated bond is definitive.
    fn all_conj_bonds_definitive_but(&mut self, ctx: &Ctx, ai: u32, bi: u32) {
        let grp = ctx.atom_grp[ai as usize];
        let mut all_definitive = true;
        for bond in ctx.mol.atom_bonds(ai) {
            if !all_definitive {
                break;
            }
            let nbi = bond.idx();
            if nbi != bi && ctx.bond_grp[nbi as usize] == grp {
                all_definitive = self.bonds[&nbi].is_definitive();
            }
        }
        if all_definitive {
            self.atoms.get_mut(&ai).expect("atom of the group").flags |= LAST_BOND;
        }
    }

    /// bit 0: the atom can accept a bond of order `bo`; bit 1: it needs a
    /// formal charge to do so
    fn can_add_bond_with_order(&self, ctx: &Ctx, ai: u32, bo: u32) -> u8 {
        let ae = &self.atoms[&ai];
        let mut can_add = u8::from(!ae.is_definitive());
        if can_add != 0 {
            can_add = u8::from(u32::from(ae.tv) <= 5 - bo);
        }
        if can_add != 0 && ae.is_last_bond() {
            // up to 2 on atoms right of carbon, not more than 1 left of N
            let right_of_c = u32::from(ae.oe() > 4);
            let fc_inc = if right_of_c != 0 {
                u32::from(!self.is_nbr_charged(ctx, ai, bo, 4))
            } else {
                let isnc = self.is_nbr_charged(ctx, ai, bo, 0);
                u32::from(
                    (!isnc && !(self.allowed_chg_left_of_n == 0 && self.total_formal_charge != 0))
                        || (isnc && bo == 3 && ae.oe() < 5),
                )
            };
            let e = ae.oe() + u32::from(ae.tv) - 1 + bo;
            can_add = u8::from(e + fc_inc + right_of_c >= 8);
            if can_add != 0 && e < 8 {
                if self.total_formal_charge != 0
                    || ctx.flags & (UNCONSTRAINED_CATIONS | UNCONSTRAINED_ANIONS) != 0
                    || right_of_c != 0
                {
                    can_add |= 1 << NEED_CHARGE_BIT;
                } else {
                    can_add = 0;
                }
            }
        }
        can_add
    }

    /// The 6-bit mask of `AtomElectrons::findAllowedBonds`: per bond order
    /// 1-3 one bit for "can accept" and one for "needs a charge".
    fn find_allowed_bonds(&mut self, ctx: &Ctx, ai: u32, bi: u32) -> u8 {
        self.all_conj_bonds_definitive_but(ctx, ai, bi);
        let mut res = 0u8;
        for i in 0..3 {
            res |= self.can_add_bond_with_order(ctx, ai, i + 1) << (i * 2);
        }
        res
    }

    /// called after all bonds to the atom have been marked as definitive
    fn finalize_atom(&mut self, ai: u32) {
        let ae = &self.atoms[&ai];
        if ae.oe() < 5 {
            let nb = ae.needed_nb_for_octet();
            if nb != 0 {
                let mut fc = (nb / 2) as i32;
                if self.allowed_chg_left_of_n != 0 {
                    if self.allowed_chg_left_of_n < 0 {
                        fc = -fc;
                    }
                    self.allowed_chg_left_of_n -= fc;
                }
            }
        }
    }

    fn check_charges_and_bond_orders(&mut self, ctx: &Ctx) -> bool {
        let mut acceptable = true;
        let mut have_incomplete_octet_right_of_c = false;
        let mut have_nitrogen_group = false;
        let mut have_pos_left_of_n = false;
        let mut have_neg_left_of_n = false;
        for ae in self.atoms.values() {
            if !acceptable {
                break;
            }
            // formal charges should be between -2 and +1
            acceptable = ae.fc < 2 && ae.fc > -3;
            if acceptable {
                if ae.fc > 0 {
                    self.flags |= HAVE_CATION;
                } else if ae.fc < 0 {
                    self.flags |= HAVE_ANION;
                }
                if ae.oe() > 4 {
                    if ae.oe() == 5 {
                        have_nitrogen_group = true;
                    }
                    if !ae.has_octet() {
                        have_incomplete_octet_right_of_c = true;
                    }
                    if ae.fc > 0 && ae.oe() > 5 {
                        self.flags |= HAVE_CATION_RIGHT_OF_N;
                    }
                } else if ae.fc < 0 {
                    have_neg_left_of_n = true;
                } else if ae.fc > 0 {
                    have_pos_left_of_n = true;
                }
                // no carbanions should coexist with incomplete octets
                acceptable = !(have_incomplete_octet_right_of_c && have_neg_left_of_n);
            }
        }
        if acceptable && have_pos_left_of_n && ctx.flags & UNCONSTRAINED_CATIONS == 0 {
            // positively charged atoms left of N with an incomplete octet
            // need a positive total charge, no anions and no nitrogens
            acceptable = self.total_formal_charge > 0
                && self.flags & HAVE_ANION == 0
                && !have_nitrogen_group;
        }
        if acceptable && have_neg_left_of_n && ctx.flags & UNCONSTRAINED_ANIONS == 0 {
            // negatively charged atoms left of N need a negative total
            // charge and no cations
            acceptable = self.total_formal_charge < 0 && self.flags & HAVE_CATION == 0;
        }
        // incomplete octets right of C are only allowed in the absence of
        // anions
        if acceptable && ctx.flags & ALLOW_INCOMPLETE_OCTETS == 0 {
            acceptable = !(have_incomplete_octet_right_of_c && self.flags & HAVE_ANION != 0);
        }
        let bond_ids: Vec<u32> = self.bonds.keys().copied().collect();
        for bi in bond_ids {
            if !acceptable {
                break;
            }
            let (bo, begin, end) = {
                let be = &self.bonds[&bi];
                (be.bo, be.begin, be.end)
            };
            for ai in [begin, end] {
                if !acceptable {
                    break;
                }
                // cumulated bonds are disallowed in aromatic rings
                if bo > 1 && self.atoms[&ai].aromatic {
                    if self.atoms[&ai].flags & HAS_MULTIPLE_BOND != 0 {
                        acceptable = false;
                    } else {
                        self.atoms.get_mut(&ai).expect("atom").flags |= HAS_MULTIPLE_BOND;
                    }
                }
                let ae = &self.atoms[&ai];
                if acceptable && ae.oe() < 5 {
                    // no carbocations on carbons bearing double or triple
                    // bonds, no double charges left of N
                    acceptable = bo == 1 || (bo > 1 && ae.fc < 1);
                }
            }
            if acceptable {
                // charged neighboring atoms left of N are not acceptable
                let (a0, a1) = (&self.atoms[&begin], &self.atoms[&end]);
                acceptable = !(a0.oe() < 5 && a0.fc != 0 && a1.oe() < 5 && a1.fc != 0);
            }
        }
        acceptable
    }

    fn compute_metrics(&mut self, ctx: &Ctx) {
        let en = &ELECTRONEGATIVITY;
        for ae in self.atoms.values() {
            self.metrics.abs_formal_charges += i32::from(ae.fc).abs();
            let an = ae.atomic_num as usize;
            self.metrics.wtd_formal_charges +=
                i32::from(ae.fc) * if an >= en.len() { 1000 } else { en[an] };
            self.metrics.nb_missing = self
                .metrics
                .nb_missing
                .wrapping_add(ae.needed_nb_for_octet() as i32);
        }
        // sums of the shortest path distances between all pairs of formal
        // charges of the same and of opposite signs
        let n = ctx.mol.num_atoms() as usize;
        let charged: Vec<(u32, i8)> = self
            .atoms
            .iter()
            .filter(|(_, ae)| ae.fc != 0)
            .map(|(&ai, ae)| (ai, ae.fc))
            .collect();
        for (i, &(a1, fc1)) in charged.iter().enumerate() {
            for &(a2, fc2) in &charged[i + 1..] {
                let dist = (ctx.dist[a1 as usize * n + a2 as usize] + 0.1) as u32 as i32;
                if i32::from(fc1) * i32::from(fc2) > 0 {
                    self.metrics.fc_same_sign_dist += dist;
                } else {
                    self.metrics.fc_opp_sign_dist += dist;
                }
            }
        }
        for (&ai, ae) in &self.atoms {
            if ae.fc != 0 {
                self.metrics.sum_formal_charge_idxs += ai as i32;
            }
        }
        for (&bi, be) in &self.bonds {
            if be.bo > 1 {
                self.metrics.sum_multiple_bond_idxs += bi as i32;
            }
        }
    }

    fn has_cation_right_of_n(&self) -> bool {
        self.flags & HAVE_CATION_RIGHT_OF_N != 0
    }

    fn has_charge_separation(&self) -> bool {
        self.flags & HAVE_CATION != 0 && self.flags & HAVE_ANION != 0
    }

    fn nb_missing(&self) -> u32 {
        self.metrics.nb_missing as u32
    }

    fn check_metrics(&self, ctx: &Ctx, stats: &mut CeStats, changed: &mut bool) -> bool {
        *changed = false;
        if !stats.init || self.nb_missing() < stats.min_nb_missing {
            *changed = stats.init;
            stats.init = true;
            stats.min_nb_missing = self.nb_missing();
        }
        if ctx.flags & ALLOW_INCOMPLETE_OCTETS != 0 || self.nb_missing() <= stats.min_nb_missing {
            if !self.has_cation_right_of_n() && !stats.have_no_cations_right_of_n {
                stats.have_no_cations_right_of_n = true;
                *changed = true;
            }
            if !self.has_charge_separation() && !stats.have_no_charge_separation {
                stats.have_no_charge_separation = true;
                *changed = true;
            }
        }
        // structures with fewer electrons than the most electron-complete
        // one are discarded unless incomplete octets are allowed
        !((ctx.flags & ALLOW_INCOMPLETE_OCTETS == 0 && self.nb_missing() > stats.min_nb_missing)
            || (ctx.flags & UNCONSTRAINED_CATIONS == 0
                && self.has_cation_right_of_n()
                && stats.have_no_cations_right_of_n)
            || (ctx.flags & ALLOW_CHARGE_SEPARATION == 0
                && self.has_charge_separation()
                && stats.have_no_charge_separation))
    }

    fn update_deg_count(&mut self, deg_count: &mut CeDegCount) {
        let m = &self.metrics;
        let key = [
            m.abs_formal_charges,
            m.fc_same_sign_dist,
            m.fc_opp_sign_dist,
            m.nb_missing,
            m.wtd_formal_charges,
        ];
        self.metrics.key = Some(key);
        *deg_count.entry(key).or_insert(0) += 1;
    }

    /// Removes the stored structures that no longer pass the metrics.
    /// Returns the keys that were removed.
    fn purge_maps(
        ctx: &Ctx,
        ce_map: &mut CeMap,
        deg_count: &mut CeDegCount,
        stats: &mut CeStats,
    ) -> Vec<Vec<u8>> {
        let mut removed = Vec::new();
        let mut changed = true;
        while changed {
            changed = false;
            for slot in 0..ce_map.slots.len() {
                let Some((fp, ce)) = &ce_map.slots[slot] else {
                    continue;
                };
                if !ce.check_metrics(ctx, stats, &mut changed) {
                    if let Some(key) = ce.metrics.key {
                        if let Some(n) = deg_count.get_mut(&key) {
                            *n -= 1;
                            if *n == 0 {
                                deg_count.remove(&key);
                            }
                        }
                    }
                    let fp = fp.clone();
                    ce_map.index.remove(&fp);
                    ce_map.slots[slot] = None;
                    removed.push(fp);
                    changed = true;
                }
                if changed {
                    break;
                }
            }
        }
        removed
    }

    /// Assigns the formal charges and, if they are acceptable, stores the
    /// structure.
    fn assign_formal_charges_and_store(
        mut self,
        ctx: &Ctx,
        ce_map: &mut CeMap,
        deg_count: &mut CeDegCount,
        stats: &mut CeStats,
        fp_flags: u32,
    ) {
        self.assign_formal_charge();
        let mut ok = self.check_charges_and_bond_orders(ctx);
        let mut changed = false;
        if ok {
            self.compute_metrics(ctx);
            ok = self.check_metrics(ctx, stats, &mut changed);
        }
        let mut key = None;
        if ok {
            let fp = pack_fp(&self, fp_flags);
            if ce_map.contains(&fp) {
                ok = false;
            } else {
                ce_map.insert(fp.clone(), self);
                key = Some(fp);
            }
        }
        if changed {
            let removed = Self::purge_maps(ctx, ce_map, deg_count, stats);
            if let Some(k) = &key {
                if removed.contains(k) {
                    ok = false;
                }
            }
        }
        if ok {
            if let Some(ce) = key.and_then(|k| ce_map.get_mut(&k)) {
                ce.update_deg_count(deg_count);
            }
        }
    }

    /// Enumerates all the permutations of non-bonded electrons.
    fn enumerate_non_bonded(
        &self,
        ctx: &Ctx,
        ce_map: &mut CeMap,
        deg_count: &mut CeDegCount,
        stats: &mut CeStats,
    ) -> Result<()> {
        // with all Kekule structures requested the bond orders take part in
        // the fingerprint too
        let fp_flags = FP_ATOMS
            | if ctx.flags & KEKULE_ALL != 0 {
                FP_BONDS
            } else {
                0
            };
        // the atoms that need non-bonded electrons to complete their octet
        let mut ai_vec = Vec::new();
        let mut nb_total = 0u32;
        for (&ai, ae) in &self.atoms {
            let nb = ae.needed_nb_for_octet();
            if nb != 0 {
                nb_total = nb_total.wrapping_add(nb);
                ai_vec.push(ai);
            }
        }
        if nb_total > self.curr_electrons {
            // some atoms will be satisfied and some won't: enumerate all
            // permutations of the unsatisfied atoms
            let missing = nb_total - self.curr_electrons;
            let num_cand = (missing - 1) / 2 + 1;
            let (num_comb, mut v) = num_comb_start_v(num_cand, ai_vec.len() as u32);
            for _ in 0..num_comb {
                let mut ce = self.clone();
                let mut vc = v;
                for &i in &ai_vec {
                    let mut e = ce.atoms[&i].needed_nb_for_octet();
                    // an atom chosen to stay octet-unsatisfied gets one
                    // electron pair less than it needs
                    if vc & 1 != 0 {
                        e = e.wrapping_sub(2);
                    }
                    ce.decr_curr_electrons(e)?;
                    ce.atoms.get_mut(&i).expect("atom").nb = e as u8;
                    vc >>= 1;
                }
                ce.assign_formal_charges_and_store(ctx, ce_map, deg_count, stats, fp_flags);
                v = update_v(v);
            }
        } else if nb_total == self.curr_electrons {
            // the assignment is univocal
            let mut ce = self.clone();
            ce.assign_non_bonded();
            ce.assign_formal_charges_and_store(ctx, ce_map, deg_count, stats, fp_flags);
        }
        // fewer electrons required than available: the bond assignment failed
        Ok(())
    }

    /// Writes the bond orders and formal charges into `mol`.
    fn assign_bonds_formal_charges_to_mol(&self, mol: &mut ROMol) -> Result<()> {
        for (&ai, ae) in &self.atoms {
            mol.atom_mut(ai).set_formal_charge(i32::from(ae.fc));
        }
        for (&bi, be) in &self.bonds {
            if !(1..=3).contains(&be.bo) {
                return Err(RDKitError::Runtime(format!(
                    "bond order for bond with index {bi} is {}; it should be between 1 and 3",
                    be.bo
                )));
            }
            mol.bond_mut(bi).set_bond_type(
                [BondType::Single, BondType::Double, BondType::Triple][be.bo as usize - 1],
            );
        }
        Ok(())
    }

    /// The ranking criteria in decreasing priority: unsatisfied octets, number
    /// of formal charges, charges weighted by electronegativity, distance
    /// between charges of the same sign (larger first) and of opposite signs
    /// (larger first), then the sums of the indices of charged atoms and of
    /// multiple bonds.
    fn rank_key(
        &self,
    ) -> (
        u32,
        u32,
        i32,
        std::cmp::Reverse<u32>,
        std::cmp::Reverse<u32>,
        u32,
        u32,
    ) {
        let m = &self.metrics;
        (
            m.nb_missing as u32,
            m.abs_formal_charges as u32,
            m.wtd_formal_charges,
            std::cmp::Reverse(m.fc_same_sign_dist as u32),
            std::cmp::Reverse(m.fc_opp_sign_dist as u32),
            m.sum_formal_charge_idxs as u32,
            m.sum_multiple_bond_idxs as u32,
        )
    }
}

/// The ranked structures of one conjugated group, grouped by degeneracy
/// (`CEVect2`).
struct CeVect2 {
    ce_vect: Vec<ConjElectrons>,
    deg_vect: Vec<u32>,
}

impl CeVect2 {
    fn new(map: CeMap) -> Self {
        let mut ce_vect = map.into_structures();
        ce_vect.sort_by_key(ConjElectrons::rank_key);
        let mut deg_vect: Vec<u32> = Vec::new();
        let mut prev = None;
        for ce in &ce_vect {
            if deg_vect.is_empty() || ce.metrics.key != prev {
                prev = ce.metrics.key;
                deg_vect.push(1);
            } else {
                *deg_vect.last_mut().expect("non-empty") += 1;
            }
        }
        CeVect2 { ce_vect, deg_vect }
    }

    fn ce_count(&self) -> usize {
        self.ce_vect.len()
    }

    fn depth(&self) -> usize {
        self.deg_vect.len()
    }

    fn ce_count_at_depth(&self, depth: usize) -> u32 {
        self.deg_vect[depth]
    }

    fn get_ce(&self, depth: usize, width: u32) -> Result<&ConjElectrons> {
        if depth >= self.deg_vect.len() {
            return Err(RDKitError::Runtime(format!(
                "depth = {depth}, d_degVect.size() = {}",
                self.deg_vect.len()
            )));
        }
        if width >= self.deg_vect[depth] {
            return Err(RDKitError::Runtime(format!(
                "width = {width}, d_degVect[{depth}] = {}",
                self.deg_vect[depth]
            )));
        }
        let i: u32 = self.deg_vect[..depth].iter().sum::<u32>() + width;
        Ok(&self.ce_vect[i as usize])
    }

    fn resize(&mut self, size: usize) {
        let keep = if size == 0 {
            0
        } else {
            self.deg_vect[..size].iter().sum::<u32>() as usize
        };
        self.ce_vect.truncate(keep);
        self.deg_vect.truncate(size);
    }

    fn idx_to_depth_width(&self, mut idx: usize) -> Result<(u32, u32)> {
        if idx >= self.ce_vect.len() {
            return Err(RDKitError::Runtime(format!(
                "idx = {idx}, d_ceVect.size() = {}",
                self.ce_vect.len()
            )));
        }
        let mut d = 0;
        while idx >= self.deg_vect[d] as usize {
            idx -= self.deg_vect[d] as usize;
            d += 1;
        }
        Ok((d as u32, idx as u32))
    }
}

// ----------------------------------------------------------------- supplier

/// `ResonanceMolSupplier`
#[allow(missing_debug_implementations)]
pub struct ResonanceMolSupplier {
    n_conj_grp: u32,
    length: u32,
    flags: u32,
    max_structs: u32,
    idx: u32,
    num_threads: u32,
    is_enumerated: bool,
    was_canceled: bool,
    ce_vect3: Vec<CeVect2>,
    mol: ROMol,
    dist: Vec<f64>,
    bond_conj_grp_idx: Vec<i32>,
    atom_conj_grp_idx: Vec<i32>,
    enum_idx: Vec<usize>,
    callback: Option<ProgressCallback>,
}

impl ResonanceMolSupplier {
    pub fn new(mol: &ROMol, flags: u32, max_structs: u32) -> Result<Self> {
        let mut flags = flags;
        if flags & UNCONSTRAINED_CATIONS != 0 {
            flags |= ALLOW_CHARGE_SEPARATION | ALLOW_INCOMPLETE_OCTETS;
        }
        if flags & UNCONSTRAINED_ANIONS != 0 {
            flags |= ALLOW_CHARGE_SEPARATION;
        }
        let max_structs = max_structs.min(MAX_STRUCTS);
        let mut mol = mol.clone();
        let dist = get_distance_mat(&mol);
        kekulize(&mut mol, false, false, 100)?;
        let mut res = ResonanceMolSupplier {
            n_conj_grp: 0,
            length: 1.min(max_structs),
            flags,
            max_structs,
            idx: 0,
            num_threads: 1,
            is_enumerated: false,
            was_canceled: false,
            ce_vect3: Vec::new(),
            mol,
            dist,
            bond_conj_grp_idx: Vec::new(),
            atom_conj_grp_idx: Vec::new(),
            enum_idx: Vec::new(),
            callback: None,
        };
        res.assign_conj_grp_idx();
        Ok(res)
    }

    /// The Kekulized molecule the supplier was created from.
    pub fn mol(&self) -> &ROMol {
        &self.mol
    }

    pub fn flags(&self) -> u32 {
        self.flags
    }

    pub fn get_num_conj_grps(&self) -> u32 {
        self.n_conj_grp
    }

    pub fn was_canceled(&self) -> bool {
        self.was_canceled
    }

    pub fn is_enumerated(&self) -> bool {
        self.is_enumerated
    }

    pub fn get_bond_conj_grp_idx(&self, bi: u32) -> Result<i32> {
        self.bond_conj_grp_idx
            .get(bi as usize)
            .copied()
            .ok_or_else(|| {
                RDKitError::Runtime(format!(
                    "d_bondConjGrpIdx.size() = {}, bi = {bi}",
                    self.bond_conj_grp_idx.len()
                ))
            })
    }

    pub fn get_atom_conj_grp_idx(&self, ai: u32) -> Result<i32> {
        self.atom_conj_grp_idx
            .get(ai as usize)
            .copied()
            .ok_or_else(|| {
                RDKitError::Runtime(format!(
                    "d_atomConjGrpIdx.size() = {}, ai = {ai}",
                    self.atom_conj_grp_idx.len()
                ))
            })
    }

    /// `setNumThreads`: 1 by default; 0 or a negative number selects
    /// relative to the hardware threads.
    pub fn set_num_threads(&mut self, num_threads: i32) {
        self.num_threads = (num_threads_to_use(num_threads) as u32).min(self.n_conj_grp);
    }

    /// Passing a callback makes the enumeration run on one thread.
    pub fn set_progress_callback(&mut self, callback: Option<ProgressCallback>) {
        self.callback = callback;
    }

    pub fn has_progress_callback(&self) -> bool {
        self.callback.is_some()
    }

    /// Each bond and atom gets the index of its conjugated group, or -1.
    fn assign_conj_grp_idx(&mut self) {
        let na = self.mol.num_atoms() as usize;
        let nb = self.mol.num_bonds() as usize;
        self.atom_conj_grp_idx = vec![-1; na];
        if nb == 0 {
            return;
        }
        self.bond_conj_grp_idx = vec![-1; nb];
        let mut n_conj_grp = 0;
        loop {
            let mut stack = Vec::new();
            for b in self.mol.bonds() {
                let i = b.idx() as usize;
                if b.is_conjugated() && self.bond_conj_grp_idx[i] == -1 {
                    stack.push(i as u32);
                    break;
                }
            }
            if stack.is_empty() {
                break;
            }
            while let Some(i) = stack.pop() {
                let bond = self.mol.bond(i);
                for atom in [bond.begin_atom_idx(), bond.end_atom_idx()] {
                    for nbr in self.mol.atom_bonds(atom) {
                        let bi_nbr = nbr.idx() as usize;
                        if nbr.is_conjugated() && self.bond_conj_grp_idx[bi_nbr] == -1 {
                            self.bond_conj_grp_idx[bi_nbr] = n_conj_grp;
                            stack.push(bi_nbr as u32);
                        }
                    }
                }
            }
            n_conj_grp += 1;
        }
        self.n_conj_grp = n_conj_grp as u32;
        for i in 0..nb {
            if self.bond_conj_grp_idx[i] != -1 {
                let bond = self.mol.bond(i as u32);
                self.atom_conj_grp_idx[bond.begin_atom_idx() as usize] = self.bond_conj_grp_idx[i];
                self.atom_conj_grp_idx[bond.end_atom_idx() as usize] = self.bond_conj_grp_idx[i];
            }
        }
    }

    /// Enumerates the resonance structures (done automatically by the first
    /// access).
    pub fn enumerate(&mut self) -> Result<()> {
        if self.is_enumerated {
            return Ok(());
        }
        let n = self.n_conj_grp as usize;
        let mut progress = ResonanceProgress {
            max_structs: self.max_structs as usize,
            groups: vec![(0, 0); n],
        };
        let ctx = Ctx {
            mol: &self.mol,
            flags: self.flags,
            max_structs: self.max_structs as usize,
            atom_grp: &self.atom_conj_grp_idx,
            bond_grp: &self.bond_conj_grp_idx,
            dist: &self.dist,
        };
        let mut canceled = false;
        let groups: Vec<u32> = (0..n as u32).collect();
        let results: Vec<Result<CeVect2>> = if self.num_threads > 1 && self.callback.is_none() {
            par_map(&groups, self.num_threads as i32, |&g| {
                build_ce_map(&ctx, g, None, &mut ResonanceProgress::default(), &mut false)
                    .map(CeVect2::new)
            })
        } else {
            let mut callback = self.callback.take();
            let res = groups
                .iter()
                .map(|&g| {
                    build_ce_map(&ctx, g, callback.as_mut(), &mut progress, &mut canceled)
                        .map(CeVect2::new)
                })
                .collect();
            self.callback = callback;
            res
        };
        let mut ce_vect3 = Vec::with_capacity(n);
        for r in results {
            ce_vect3.push(r?);
        }
        self.ce_vect3 = ce_vect3;
        self.was_canceled |= canceled;
        self.set_length();
        self.trim_ce_vect2();
        self.prep_enum_idx_vect()?;
        self.is_enumerated = true;
        Ok(())
    }

    fn set_length(&mut self) {
        for v in &self.ce_vect3 {
            if self.length >= self.max_structs {
                break;
            }
            let p = u64::from(self.length) * v.ce_count() as u64;
            self.length = if p < u64::from(self.max_structs) {
                p as u32
            } else {
                self.max_structs
            };
        }
    }

    /// If the number of resonance structures reaches the maximum, keeps only
    /// the most likely structures of each group.
    fn trim_ce_vect2(&mut self) {
        if self.length != self.max_structs {
            return;
        }
        let n = self.n_conj_grp as usize;
        let mut s = vec![0usize; n];
        let mut t = vec![0u32; n];
        let mut curr_size = 0u32;
        while curr_size < self.length {
            curr_size = 1;
            let mut g = 0;
            while curr_size < self.length && g < n {
                if s[g] < self.ce_vect3[g].depth() {
                    t[g] += self.ce_vect3[g].ce_count_at_depth(s[g]);
                    s[g] += 1;
                }
                curr_size = curr_size.wrapping_mul(t[g]);
                g += 1;
            }
        }
        for g in 0..n {
            self.ce_vect3[g].resize(s[g]);
        }
    }

    /// The (depth, width) of the structure of every group that makes up the
    /// complete structure `idx`.
    fn idx_to_ce_perm(&self, mut idx: usize) -> Result<Vec<u32>> {
        let n = self.n_conj_grp as usize;
        let mut c = vec![0u32; n * 2];
        for g in (0..n).rev() {
            let d: usize = self.ce_vect3[..g].iter().map(CeVect2::ce_count).product();
            if d != 0 && self.ce_vect3[g].ce_count() != 0 {
                let (depth, width) = self.ce_vect3[g].idx_to_depth_width(idx / d)?;
                c[g * 2] = depth;
                c[g * 2 + 1] = width;
                idx %= d;
            }
        }
        Ok(c)
    }

    /// Sorts the complete structures so that the most likely combinations
    /// come first, in a breadth-first fashion.
    fn prep_enum_idx_vect(&mut self) -> Result<()> {
        let mut perms = Vec::with_capacity(self.length as usize);
        for i in 0..self.length as usize {
            perms.push((i, self.idx_to_ce_perm(i)?));
        }
        perms.sort_by(|(_, a), (_, b)| {
            let sum = |v: &[u32]| v.iter().step_by(2).sum::<u32>();
            let max = |v: &[u32]| v.iter().step_by(2).copied().max().unwrap_or(0);
            sum(a)
                .cmp(&sum(b))
                .then_with(|| max(a).cmp(&max(b)))
                .then_with(|| a.iter().step_by(2).cmp(b.iter().step_by(2)))
                .then_with(|| a.iter().skip(1).step_by(2).cmp(b.iter().skip(1).step_by(2)))
        });
        self.enum_idx = perms.into_iter().map(|(i, _)| i).collect();
        Ok(())
    }

    /// The number of resonance structures.
    pub fn length(&mut self) -> Result<u32> {
        self.enumerate()?;
        Ok(self.length)
    }

    pub fn reset(&mut self) -> Result<()> {
        self.enumerate()?;
        self.idx = 0;
        Ok(())
    }

    pub fn at_end(&mut self) -> Result<bool> {
        self.enumerate()?;
        Ok(self.idx == self.length)
    }

    /// The next resonance structure, or `None` at the end.
    pub fn next_mol(&mut self) -> Result<Option<ROMol>> {
        self.enumerate()?;
        if self.idx == self.length {
            return Ok(None);
        }
        let i = self.idx;
        self.idx += 1;
        self.get(i).map(Some)
    }

    pub fn move_to(&mut self, idx: u32) -> Result<()> {
        self.enumerate()?;
        if idx >= self.length {
            return Err(RDKitError::Runtime(format!(
                "d_length = {}, idx = {idx}",
                self.length
            )));
        }
        self.idx = idx;
        Ok(())
    }

    /// The resonance structure number `idx`: the complete structures combine
    /// the ranked structures of the groups breadth first.
    pub fn get(&mut self, idx: u32) -> Result<ROMol> {
        self.enumerate()?;
        self.structure(idx)
    }

    /// Like [`get`](Self::get) for a supplier that is already enumerated.
    fn structure(&self, idx: u32) -> Result<ROMol> {
        if idx >= self.length {
            return Err(RDKitError::Runtime(format!(
                "d_length = {}, idx = {idx}",
                self.length
            )));
        }
        let c = self.idx_to_ce_perm(self.enum_idx[idx as usize])?;
        let mut mol = self.mol.clone();
        for g in 0..self.n_conj_grp as usize {
            self.ce_vect3[g]
                .get_ce(c[g * 2] as usize, c[g * 2 + 1])?
                .assign_bonds_formal_charges_to_mol(&mut mol)?;
        }
        // fix the number of explicit and implicit Hs
        mol.clear_computed_props(false);
        for i in 0..mol.num_atoms() {
            mol.atom_mut(i).clear_computed_props();
            let hs = mol.atom(i).num_implicit_hs()? + mol.atom(i).num_explicit_hs();
            mol.atom_mut(i).set_num_explicit_hs(hs);
            mol.update_atom_property_cache(i, true)?;
        }
        sanitize_mol(&mut mol, SANITIZE_FINDRADICALS | SANITIZE_ADJUSTHS).map_err(|(_, e)| e)?;
        Ok(mol)
    }

    /// `SubstructMatch(resMolSupplier, query, params)`: the matches of
    /// `query` in all the resonance structures, in sorted order. With
    /// `uniquify` only the match with the smallest atom order is kept for
    /// each set of atoms.
    pub fn substruct_matches(
        &mut self,
        query: &ROMol,
        params: &SubstructMatchParameters,
    ) -> Result<Vec<MatchVect>> {
        let length = self.length()?;
        let nt = (length as usize)
            .min(num_threads_to_use(params.num_threads))
            .max(1);
        let this: &Self = self;
        if nt == 1 {
            let mut matches = BTreeSet::new();
            res_substruct_chunk(this, query, params, &mut matches, 0, length)?;
            return Ok(matches.into_iter().collect());
        }
        let per_thread = f64::from(length) / nt as f64;
        let mut bounds = Vec::with_capacity(nt);
        let (mut ei, mut dc) = (0u32, 0.0f64);
        for _ in 0..nt {
            let bi = ei;
            dc += per_thread;
            ei = dc.floor() as u32;
            bounds.push((bi, ei));
        }
        let partial: Vec<Result<BTreeSet<MatchVect>>> = std::thread::scope(|s| {
            let handles: Vec<_> = bounds
                .iter()
                .map(|&(bi, ei)| {
                    s.spawn(move || {
                        let mut matches = BTreeSet::new();
                        res_substruct_chunk(this, query, params, &mut matches, bi, ei)
                            .map(|()| matches)
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().expect("substructure search panicked"))
                .collect()
        });
        let mut matches = BTreeSet::new();
        for chunk in partial {
            for m in chunk? {
                if !try_to_insert(&mut matches, m, params) {
                    break;
                }
            }
        }
        Ok(matches.into_iter().collect())
    }
}

fn atoms_of(m: &MatchVect) -> BTreeSet<i32> {
    m.iter().map(|&(_, a)| a).collect()
}

/// `tryToInsert`: false once `maxMatches` matches are collected.
fn try_to_insert(
    matches: &mut BTreeSet<MatchVect>,
    m: MatchVect,
    params: &SubstructMatchParameters,
) -> bool {
    if matches.len() == params.max_matches as usize {
        return false;
    }
    if !params.uniquify {
        matches.insert(m);
        return true;
    }
    let atoms = atoms_of(&m);
    let same = matches.iter().find(|e| atoms_of(e) == atoms).cloned();
    match same {
        Some(existing) if m < existing => {
            matches.remove(&existing);
            matches.insert(m);
        }
        Some(_) => {}
        None => {
            matches.insert(m);
        }
    }
    true
}

/// `ResSubstructMatchHelper_`: the structures `bi..ei`
fn res_substruct_chunk(
    suppl: &ResonanceMolSupplier,
    query: &ROMol,
    params: &SubstructMatchParameters,
    matches: &mut BTreeSet<MatchVect>,
    bi: u32,
    ei: u32,
) -> Result<()> {
    let mut i = bi;
    while matches.len() < params.max_matches as usize && i < ei {
        let mol = suppl.structure(i)?;
        for m in substruct_match(&mol, query, params)? {
            if !try_to_insert(matches, m, params) {
                break;
            }
        }
        i += 1;
    }
    Ok(())
}

/// Enumerates all the multiple-bond arrangements of conjugated group
/// `group_idx` (depth first), distributes the non-bonded electrons and keeps
/// the plausible structures.
fn build_ce_map(
    ctx: &Ctx,
    group_idx: u32,
    mut callback: Option<&mut ProgressCallback>,
    progress: &mut ResonanceProgress,
    canceled: &mut bool,
) -> Result<CeMap> {
    let mol = ctx.mol;
    let nb = mol.num_bonds();
    let na = mol.num_atoms();
    let fp_flags = FP_ATOMS
        | if ctx.flags & KEKULE_ALL != 0 {
            FP_BONDS
        } else {
            0
        };
    let fp_flags_tmp = if ctx.flags & KEKULE_ALL != 0 {
        FP_BONDS
    } else {
        FP_ATOMS
    };
    let mut ce_map = CeMap::default();
    let mut ce_set_tmp: HashSet<Vec<u8>> = HashSet::new();
    let mut deg_count = CeDegCount::new();
    let mut stats = CeStats::default();
    let mut ce = ConjElectrons::new(ctx, group_idx)?;
    let ce_copy = ce.clone();
    // the first structure has the user-supplied bond and charge arrangement
    ce.init_from_mol(ctx)?;
    // the result is ignored, but the call sets the HAVE_CATION_RIGHT_OF_N flag
    ce.check_charges_and_bond_orders(ctx);
    ce.compute_metrics(ctx);
    ce.update_deg_count(&mut deg_count);
    ce_map.insert(pack_fp(&ce, fp_flags), ce);
    // partial bond arrangements are kept on a stack whenever multiple paths
    // arise; each structure has a stack of atoms to continue from
    let mut ce_stack: Vec<ConjElectrons> = vec![ce_copy];
    'stack: while let Some(mut ce) = ce_stack.pop() {
        // prime the atom stack with any atom of this group
        if ce.begin_stack.is_empty() {
            let Some(first) = (0..na).find(|&a| ctx.atom_grp[a as usize] == group_idx as i32)
            else {
                continue;
            };
            ce.push_to_begin_stack(first);
        }
        let mut dead = false;
        while let Some(ai_begin) = ce.pop_from_begin_stack() {
            let mut ai = [ai_begin, na];
            let mut bi = nb;
            let mut found = false;
            for ai_nbr in mol.atom_neighbors(ai_begin) {
                if ctx.atom_grp[ai_nbr as usize] != group_idx as i32 {
                    continue;
                }
                // skip neighbors already dealt with
                match ce.atoms.get(&ai_nbr) {
                    Some(ae) if !ae.is_definitive() => {}
                    _ => continue,
                }
                let bi_nbr = mol
                    .get_bond_between_atoms(ai_begin, ai_nbr)?
                    .expect("neighbors are bonded")
                    .idx();
                // skip bonds whose order is already assigned
                match ce.bonds.get(&bi_nbr) {
                    Some(be) if !be.is_definitive() => {}
                    _ => continue,
                }
                if !found {
                    found = true;
                    bi = bi_nbr;
                    ai[1] = ai_nbr;
                } else {
                    ce.push_to_begin_stack(ai_nbr);
                }
            }
            if !found {
                continue;
            }
            let mut allowed_per_atom = [0u8; 2];
            for i in 0..2 {
                allowed_per_atom[i] = ce.find_allowed_bonds(ctx, ai[i], bi);
                if ce.atoms[&ai[i]].is_last_bond() {
                    ce.atoms.get_mut(&ai[i]).expect("atom").flags |= DEFINITIVE;
                } else {
                    ce.push_to_begin_stack(ai[i]);
                }
            }
            let allowed = allowed_per_atom[0] & allowed_per_atom[1];
            let last = [
                ce.atoms[&ai[0]].is_last_bond(),
                ce.atoms[&ai[1]].is_last_bond(),
            ];
            let oe = [ce.atoms[&ai[0]].oe(), ce.atoms[&ai[1]].oe()];
            let ce_copy = ce.clone();
            let mut any_allowed = false;
            let mut need_to_fork = false;
            for i in 0..3u32 {
                let t = i * 2;
                let order_mask = 1u8 << t;
                let chg_mask = 1u8 << (t + 1);
                let bo = i + 1;
                // enough electrons, both atoms can accept the bond type and
                // it does not require charges on atoms left of N
                if ce.curr_electrons >= bo * 2
                    && allowed & order_mask != 0
                    && !(allowed & chg_mask != 0 && (oe[0] < 5 || oe[1] < 5))
                {
                    any_allowed = true;
                    // a further feasible bond type forks the structure
                    let target: &mut ConjElectrons = if !need_to_fork {
                        need_to_fork = true;
                        &mut ce
                    } else {
                        ce_stack.push(ce_copy.clone());
                        ce_stack.last_mut().expect("just pushed")
                    };
                    target.set_order(bi, bo as u8)?;
                    for j in 0..2 {
                        if last[j] {
                            target.finalize_atom(ai[j]);
                        }
                    }
                }
            }
            // a dead end discards this structure
            if !any_allowed {
                dead = true;
                break;
            }
        }
        if dead {
            continue;
        }
        // a bond arrangement stored before is discarded
        if !ce_set_tmp.insert(pack_fp(&ce, fp_flags_tmp)) {
            continue;
        }
        ce.enumerate_non_bonded(ctx, &mut ce_map, &mut deg_count, &mut stats)?;
        if let Some(cb) = callback.as_deref_mut() {
            progress.groups[group_idx as usize] = (ce_map.len(), deg_count.len());
            if !cb(progress) {
                *canceled = true;
                break 'stack;
            }
        }
        // stop early once the number of non-degenerate structures exceeds
        // the maximum
        if deg_count.len() >= ctx.max_structs {
            break;
        }
    }
    Ok(ce_map)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gosper_enumerates_combinations() {
        let (n, mut v) = num_comb_start_v(2, 4);
        assert_eq!(n, 6);
        let mut seen = Vec::new();
        for _ in 0..n {
            seen.push(v);
            v = update_v(v);
        }
        assert_eq!(seen, [0b11, 0b101, 0b110, 0b1001, 0b1010, 0b1100]);
    }
}
