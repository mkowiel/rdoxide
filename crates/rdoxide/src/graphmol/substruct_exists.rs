//! Existence-only substructure search (`HasSubstructMatch`).
//!
//! Whether a match exists does not depend on the order the query atoms are
//! placed in, so unlike [`substruct_match`](super::substruct::substruct_match),
//! which follows upstream's VF2 node order to return the same matches in the
//! same order, this search picks its own: it starts at the query atom whose
//! element is rarest in the molecule and grows along bonds, so that molecules
//! without a match fail early. The idea (rare root first, atom tests in one
//! cheap comparison before the full query, results of atom and bond tests
//! cached) follows <https://github.com/mkowiel/roaring-substructure>.
//!
//! The atom and bond tests are the matcher's own, so the answer is the one
//! the VF2 search gives. Parameters whose meaning depends on the whole
//! mapping (chirality, final-check and generic matchers) or on callbacks and
//! query-query matching are left to the VF2 search.

use crate::Result;
use crate::graphmol::queries::{AtomFn, AtomQuery, QueryKind, get_atom_type_atomic_num};
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::{
    Ctx, SubstructMatchParameters, prepare_recursion, substruct_match,
};

/// Heaviest element tracked by the element counts; heavier ones share a slot.
const MAX_Z: usize = 255;

/// The atomic number every molecule atom matching `q` must have, if the
/// query says so.
fn required_atomic_num(q: &AtomQuery) -> Option<u32> {
    if q.negate {
        return None;
    }
    match q.kind {
        QueryKind::Equality if q.tol == 0 => match q.data_func {
            Some(AtomFn::Num) => u32::try_from(q.val).ok(),
            Some(AtomFn::Type) => u32::try_from(get_atom_type_atomic_num(q.val)).ok(),
            _ => None,
        },
        QueryKind::And => q.children.iter().find_map(required_atomic_num),
        _ => None,
    }
}

/// One placed query atom: the bond to an earlier atom that proposes its
/// candidates (`anchor`) and the bonds to other earlier atoms that must exist
/// (`closures`, a range of [`Scratch::pool`]).
struct Step {
    atom: u32,
    z: Option<u32>,
    anchor: Option<(usize, u32)>,
    closures: std::ops::Range<usize>,
}

/// Buffers reused between searches, so that a search of a small molecule
/// allocates nothing.
#[derive(Default)]
struct Scratch {
    zs: Vec<Option<u32>>,
    step_of: Vec<usize>,
    order: Vec<usize>,
    steps: Vec<Step>,
    pool: Vec<(usize, u32)>,
    back: Vec<(usize, u32)>,
    mapped: Vec<u32>,
    used: Vec<bool>,
    atom_cache: Vec<u8>,
    bond_cache: Vec<u8>,
}

thread_local! {
    static SCRATCH: std::cell::RefCell<Scratch> = std::cell::RefCell::default();
}

impl Scratch {
    /// Orders the query atoms: the rarest first, then growing along bonds.
    fn plan(&mut self, query: &ROMol, have: &[u32; MAX_Z + 1], n_mol: usize) {
        let n = query.num_atoms() as usize;
        self.zs.clear();
        self.zs.extend((0..n as u32).map(|i| {
            let a = query.atom(i);
            match a.query() {
                Some(q) => required_atomic_num(q),
                None => Some(a.atomic_num()),
            }
        }));
        let zs = &self.zs;
        let rarity = |i: usize| match zs[i] {
            Some(z) => have[(z as usize).min(MAX_Z)] as usize,
            None => n_mol,
        };
        let degree = |i: usize| query.adj[i].len();
        let (step_of, order) = (&mut self.step_of, &mut self.order);
        step_of.clear();
        step_of.resize(n, usize::MAX);
        order.clear();
        while order.len() < n {
            // a new component starts at its rarest atom
            let root = (0..n)
                .filter(|&i| step_of[i] == usize::MAX)
                .min_by_key(|&i| (rarity(i), std::cmp::Reverse(degree(i))))
                .expect("an unplaced atom");
            step_of[root] = order.len();
            order.push(root);
            // grow by the atom with the most bonds to placed atoms, then the rarest
            while let Some(next) = (0..n)
                .filter(|&i| step_of[i] == usize::MAX)
                .filter_map(|i| {
                    let placed = query.adj[i]
                        .iter()
                        .filter(|&&(m, _)| step_of[m as usize] != usize::MAX)
                        .count();
                    (placed > 0).then_some((std::cmp::Reverse(placed), rarity(i), i))
                })
                .min()
                .map(|(_, _, i)| i)
            {
                step_of[next] = order.len();
                order.push(next);
            }
        }
        self.steps.clear();
        self.pool.clear();
        for (k, &i) in self.order.iter().enumerate() {
            self.back.clear();
            self.back.extend(
                query.adj[i]
                    .iter()
                    .filter(|&&(m, _)| self.step_of[m as usize] < k)
                    .map(|&(m, b)| (self.step_of[m as usize], b)),
            );
            self.back.sort_unstable();
            let start = self.pool.len();
            self.pool.extend(self.back.iter().skip(1));
            self.steps.push(Step {
                atom: i as u32,
                z: self.zs[i],
                anchor: self.back.first().copied(),
                closures: start..self.pool.len(),
            });
        }
    }
}

struct Search<'a> {
    ctx: Ctx<'a>,
    steps: &'a [Step],
    pool: &'a [(usize, u32)],
    /// Target atom of each placed step.
    mapped: &'a mut [u32],
    used: &'a mut [bool],
    /// Results of the atom and bond tests: 0 not yet known, 1 passes, 2 fails.
    atom_cache: &'a mut [u8],
    bond_cache: &'a mut [u8],
}

impl Search<'_> {
    fn atom_ok(&mut self, qi: u32, mi: u32) -> Result<bool> {
        let (query, mol) = (self.ctx.query, self.ctx.mol);
        if query.adj[qi as usize].len() > mol.adj[mi as usize].len() {
            return Ok(false);
        }
        let slot = qi as usize * mol.num_atoms() as usize + mi as usize;
        if self.atom_cache[slot] == 0 {
            self.atom_cache[slot] = 1 + !self.ctx.atom_label(qi, mi)? as u8;
        }
        Ok(self.atom_cache[slot] == 1)
    }

    fn bond_ok(&mut self, qb: u32, mb: u32) -> Result<bool> {
        let slot = qb as usize * self.ctx.mol.num_bonds() as usize + mb as usize;
        if self.bond_cache[slot] == 0 {
            self.bond_cache[slot] = 1 + !self.ctx.bond_label(qb, mb)? as u8;
        }
        Ok(self.bond_cache[slot] == 1)
    }

    /// Whether `cand` can take step `k`, given the bonds of the step.
    fn accepts(&mut self, k: usize, cand: u32) -> Result<bool> {
        let (steps, pool) = (self.steps, self.pool);
        let step = &steps[k];
        if self.used[cand as usize] {
            return Ok(false);
        }
        if step
            .z
            .is_some_and(|z| self.ctx.mol.atom(cand).atomic_num() != z)
        {
            return Ok(false);
        }
        if !self.atom_ok(step.atom, cand)? {
            return Ok(false);
        }
        for &(s, qb) in &pool[step.closures.clone()] {
            let from = self.mapped[s];
            let Some(mb) = self.ctx.mol.bond_idx_between(from, cand) else {
                return Ok(false);
            };
            if !self.bond_ok(qb, mb)? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn place(&mut self, k: usize, cand: u32) -> Result<bool> {
        self.mapped[k] = cand;
        self.used[cand as usize] = true;
        let found = self.search(k + 1)?;
        self.used[cand as usize] = false;
        Ok(found)
    }

    fn search(&mut self, k: usize) -> Result<bool> {
        let steps = self.steps;
        let Some(step) = steps.get(k) else {
            return Ok(true);
        };
        let mol = self.ctx.mol;
        match step.anchor {
            Some((s, qb)) => {
                for &(cand, mb) in &mol.adj[self.mapped[s] as usize] {
                    if self.bond_ok(qb, mb)? && self.accepts(k, cand)? && self.place(k, cand)? {
                        return Ok(true);
                    }
                }
            }
            None => {
                for cand in 0..mol.num_atoms() {
                    if self.accepts(k, cand)? && self.place(k, cand)? {
                        return Ok(true);
                    }
                }
            }
        }
        Ok(false)
    }
}

/// Whether `query` matches somewhere in `mol`: `!substruct_match(..).is_empty()`
/// without building the match.
pub fn substruct_has_match(
    mol: &ROMol,
    query: &ROMol,
    ps: &SubstructMatchParameters,
) -> Result<bool> {
    let (mn, qn) = (mol.num_atoms() as usize, query.num_atoms() as usize);
    if mn == 0 || qn == 0 || qn > mn {
        return Ok(false);
    }
    let order_independent = !ps.use_chirality
        && !ps.use_generic_matchers
        && !ps.use_query_query_matches
        && ps.extra_final_check.is_none()
        && ps.extra_atom_check.is_none()
        && ps.extra_bond_check.is_none();
    if !order_independent {
        let mut ps = ps.clone();
        ps.max_matches = 1;
        return Ok(!substruct_match(mol, query, &ps)?.is_empty());
    }
    // every query atom needs an atom of its own with the right element
    let mut have = [0u32; MAX_Z + 1];
    for a in mol.atoms() {
        have[(a.atomic_num() as usize).min(MAX_Z)] += 1;
    }
    SCRATCH.with(|cell| match cell.try_borrow_mut() {
        Ok(mut scratch) => search_with(&mut scratch, mol, query, ps, &have),
        Err(_) => search_with(&mut Scratch::default(), mol, query, ps, &have),
    })
}

fn search_with(
    scratch: &mut Scratch,
    mol: &ROMol,
    query: &ROMol,
    ps: &SubstructMatchParameters,
    have: &[u32; MAX_Z + 1],
) -> Result<bool> {
    let (mn, qn) = (mol.num_atoms() as usize, query.num_atoms() as usize);
    scratch.plan(query, have, mn);
    let mut need = [0u32; MAX_Z + 1];
    for z in scratch.steps.iter().filter_map(|s| s.z) {
        need[(z as usize).min(MAX_Z)] += 1;
    }
    if need.iter().zip(have).any(|(n, h)| n > h) {
        return Ok(false);
    }
    let rec = prepare_recursion(mol, query, ps)?;
    let Scratch {
        steps,
        pool,
        mapped,
        used,
        atom_cache,
        bond_cache,
        ..
    } = scratch;
    mapped.clear();
    mapped.resize(qn, 0);
    used.clear();
    used.resize(mn, false);
    atom_cache.clear();
    atom_cache.resize(qn * mn, 0);
    bond_cache.clear();
    bond_cache.resize(query.num_bonds() as usize * mol.num_bonds() as usize, 0);
    Search {
        ctx: Ctx::new(query, mol, ps, &rec),
        steps,
        pool,
        mapped,
        used,
        atom_cache,
        bond_cache,
    }
    .search(0)
}
