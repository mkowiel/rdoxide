//! `ssmatch.c` (fingerprint half) and `hashcode.c`: the Avalon fingerprints
//! (`SetFingerprintBits`, `SetFingerprintCountsWithFocus`). Copyright 2001-2011 Novartis Pharma AG.

use super::perceive::{
    AROMATIC, DOUBLE, Neighbourhood, SINGLE, TRIPLE, perceive_aromatic_bonds,
    perceive_dy_aromaticity, ring_state, set_ring_size_flags, setup_neighbourhood,
};
use super::reaccs::{Molecule, NONE, number_of};
use super::utilities::{atom_symbol_match, compute_implicit_h};

/// feature selection flags of `ssmatch.h`
pub const USE_RING_PATTERN: i32 = 0x000001;
pub const USE_RING_PATH: i32 = 0x000002;
pub const USE_ATOM_SYMBOL_PATH: i32 = 0x000004;
pub const USE_ATOM_CLASS_PATH: i32 = 0x000008;
pub const USE_ATOM_COUNT: i32 = 0x000010;
pub const USE_AUGMENTED_ATOM: i32 = 0x000020;
pub const USE_HCOUNT_PATH: i32 = 0x000040;
pub const USE_HCOUNT_CLASS_PATH: i32 = 0x000080;
pub const USE_HCOUNT_PAIR: i32 = 0x000100;
pub const USE_BOND_PATH: i32 = 0x000200;
pub const USE_AUGMENTED_BOND: i32 = 0x000400;
pub const USE_RING_SIZE_COUNTS: i32 = 0x000800;
pub const USE_DEGREE_PATH: i32 = 0x001000;
pub const USE_CLASS_SPIDERS: i32 = 0x002000;
pub const USE_FEATURE_PAIRS: i32 = 0x004000;
pub const USE_ALL_FEATURES: i32 = 0x007FFF;
pub const USE_SCAFFOLD_IDS: i32 = 0x100000;
pub const USE_SCAFFOLD_COLORS: i32 = 0x200000;
pub const USE_SCAFFOLD_LINKS: i32 = 0x400000;
pub const USE_SHORTCUT_LABELS: i32 = 0x800000;
pub const USE_NON_SSS_BITS: i32 = 0xF00000;
/// `fpflags`
pub const USE_DY_AROMATICITY: i32 = 0x0001;
pub const ACCUMULATE_BITS: i32 = 0x0002;

const ATOM_COUNT_SEED: u64 = 31;
const AUGMENTED_ATOM_SEED: u64 = 37;
const HCOUNT_PAIR_SEED: u64 = 47;
const AUGMENTED_BOND_SEED: u64 = 61;
const ATOM_SYMBOL_PATH_SEED: u64 = 17;
const HCOUNT_PATH_SEED: u64 = 41;
const RING_PATH_SEED: u64 = 13;
const BOND_PATH_SEED: u64 = 53;
const HCOUNT_CLASS_PATH_SEED: u64 = 43;
const ATOM_CLASS_PATH_SEED: u64 = 23;
const SPECIAL_RING: u32 = 0xFC & !(1 << 6);
const RING_CLOSURE_SEED: i64 = 101;

const PROCESS_RING_CLOSURES: i32 = 0x0001;
const PROCESS_CHAINS: i32 = 0x0002;
const FORCED_HETERO_END: i32 = 0x0004;
const IGNORE_PATH_SYMBOL: i32 = 0x0008;
const IGNORE_TERM_SYMBOL: i32 = 0x0010;
const FORCED_RING_PATH: i32 = 0x0020;
const STOP_AT_HEAVY_ATOM: i32 = 0x0040;
const ANY_COLOR: i32 = 113;

const SUB_AS_IS: i32 = -2;
const SUB_ONE: i32 = 1;
const SUB_MORE: i32 = 6;
const ZERO_COUNT: i32 = 1;

/// `next_hash`
pub fn next_hash(mut hash: u64, data: u64) -> u64 {
    hash = hash.wrapping_add(data);
    hash = hash.wrapping_add(hash << 10);
    hash ^= hash >> 6;
    hash
}

/// `hash_position`
pub fn hash_position(mut hash: u64, nslots: usize) -> usize {
    hash = hash.wrapping_add(hash << 3);
    hash ^= hash >> 11;
    hash = hash.wrapping_add(hash << 15);
    (hash % nslots as u64) as usize
}

/// `hash_string` (the characters are signed, like `char` of the reference builds)
pub fn hash_string(s: &str) -> u64 {
    let mut hash = 1001u64;
    for b in s.bytes() {
        hash = next_hash(hash, i64::from(b as i8) as u64);
    }
    hash
}

/// `NEXT_SEED(seed, increment)` with a C `int` increment
fn ns(seed: u64, inc: i32) -> u64 {
    next_hash(seed, i64::from(inc) as u64)
}

/// `GuessHCountsFromSubstitution`
fn guess_h_counts_from_substitution(mp: &mut Molecule, nbp: &[Neighbourhood]) {
    for i in 0..mp.atoms.len() {
        if mp.atoms[i].query_h_count != NONE || mp.atoms[i].charge != 0 || mp.atoms[i].radical != 0
        {
            continue;
        }
        let (mut nsingle, mut ndouble, mut naromatic, mut ntriple, mut nother) = (0, 0, 0, 0, 0);
        let mut nexplicit = 0;
        for &aj in &nbp[i].atoms {
            if atom_symbol_match(&mp.atoms[aj].atom_symbol, "H,D,T") {
                nexplicit += 1;
            }
        }
        for &bj in &nbp[i].bonds {
            match mp.bonds[bj].bond_type {
                SINGLE => nsingle += 1,
                DOUBLE => ndouble += 1,
                AROMATIC => naromatic += 1,
                TRIPLE => ntriple += 1,
                _ => nother += 1,
            }
        }
        if nother > 0 {
            continue;
        }
        let ap = &mut mp.atoms[i];
        match ap.atom_symbol.as_str() {
            "C" => {
                if naromatic == 2 {
                    nsingle += 1;
                    ndouble += 1;
                    naromatic = 0;
                }
                if naromatic > 0 {
                    continue;
                }
                if ap.sub_desc == NONE {
                    if nsingle + 2 * ndouble + 3 * ntriple == 4 {
                        ap.sub_desc = SUB_AS_IS;
                    }
                } else if ap.sub_desc == SUB_AS_IS {
                    ap.query_h_count =
                        nexplicit + ZERO_COUNT + 4 - nsingle - 2 * ndouble - 3 * ntriple;
                }
            }
            "O" => {
                if ntriple + naromatic != 0 {
                    continue;
                }
                if ap.sub_desc == SUB_AS_IS {
                    ap.query_h_count = nexplicit + ZERO_COUNT + 2 - nsingle - 2 * ndouble;
                }
            }
            "N" => {
                if naromatic != 0 {
                    continue;
                }
                if ap.sub_desc == SUB_AS_IS {
                    ap.query_h_count = nexplicit + ZERO_COUNT + 3 - nsingle - 2 * ndouble;
                }
            }
            _ => {}
        }
    }
}

struct PathCtx<'a> {
    mp: &'a Molecule,
    nbp: &'a [Neighbourhood],
    counts: &'a mut [i32],
    touched: &'a mut [i32],
    minbonds: i32,
    maxbonds: i32,
    flags: i32,
    exclude_atom: i32,
}

impl PathCtx<'_> {
    fn add_bit(&mut self, seed: u64) {
        let p = hash_position(seed, self.counts.len());
        self.counts[p] += 1;
    }
}

/// `SetPathBitsRec`
fn set_path_bits_rec(
    c: &mut PathCtx,
    mut seed: u64,
    nbonds: i32,
    sprout_index: usize,
    last_index: i32,
) -> i32 {
    let mut result = 0;
    if nbonds > c.maxbonds {
        return result;
    }
    let flags = c.flags;
    let nb = &c.nbp[sprout_index];
    for i in 0..nb.atoms.len() {
        let ai = nb.atoms[i];
        if ai as i32 == last_index || ai as i32 + 1 == c.exclude_atom {
            continue;
        }
        let bi = nb.bonds[i];
        let ap_color = c.mp.atoms[ai].color;
        if flags & FORCED_RING_PATH != 0 && c.mp.bonds[bi].rsize_flags == 0 {
            continue;
        }
        if ap_color >= 18 && ap_color != ANY_COLOR && flags & STOP_AT_HEAVY_ATOM != 0 {
            continue;
        }
        let bcolor = c.mp.bonds[bi].color;
        if c.touched[ai] > 0 {
            if flags & PROCESS_RING_CLOSURES == 0 || c.touched[ai] > 1 || bcolor == 0 {
                continue;
            }
            let acolor = ap_color;
            if acolor == 0 && flags & IGNORE_PATH_SYMBOL == 0 {
                continue;
            }
            let old_seed = seed;
            seed = ns(seed, bcolor * 16);
            if flags & IGNORE_PATH_SYMBOL == 0 {
                seed = ns(seed, acolor);
            }
            seed = ns(seed, c.touched[ai]);
            c.add_bit(seed);
            seed = next_hash(
                seed,
                (RING_CLOSURE_SEED * i64::from(nbonds - c.touched[ai])) as u64,
            );
            c.add_bit(seed);
            result += 1;
            seed = old_seed;
        } else {
            if bcolor == 0 {
                continue;
            }
            let acolor = ap_color;
            if acolor == 0 && flags & IGNORE_PATH_SYMBOL == 0 {
                continue;
            }
            c.touched[ai] = nbonds + 1;
            let old_seed = seed;
            seed = ns(seed, bcolor * 16);
            if flags & IGNORE_PATH_SYMBOL == 0 {
                seed = ns(seed, acolor);
            }
            if nbonds >= c.minbonds
                && (acolor > 0 || flags & IGNORE_TERM_SYMBOL != 0)
                && flags & PROCESS_CHAINS != 0
                && (flags & FORCED_HETERO_END == 0
                    || flags & IGNORE_TERM_SYMBOL != 0
                    || acolor != 6)
                && acolor > 1
            {
                if flags & IGNORE_TERM_SYMBOL != 0 {
                    c.add_bit(seed);
                } else {
                    c.add_bit(ns(seed, 17 * acolor));
                }
                result += 1;
            }
            if nbonds < c.maxbonds {
                result += set_path_bits_rec(c, seed, nbonds + 1, ai, sprout_index as i32);
            }
            c.touched[ai] = 0;
            seed = old_seed;
        }
    }
    result
}

/// `SetFingerprintCountsWithFocus`: multiplicities of the hashed features in `fp_counts`; `Err` if
/// a requested feature class is not ported yet. The atom colours, ring size flags and query
/// hydrogen counts of `mp` are left changed, as in the C code.
pub fn set_fingerprint_counts_with_focus(
    mp: &mut Molecule,
    fp_counts: &mut [i32],
    which_bits: i32,
    as_query: bool,
    fpflags: i32,
    exclude_atom: i32,
) -> Result<i32, &'static str> {
    const IMPLEMENTED: i32 = USE_ATOM_COUNT
        | USE_ATOM_SYMBOL_PATH
        | USE_AUGMENTED_ATOM
        | USE_AUGMENTED_BOND
        | USE_HCOUNT_PAIR
        | USE_HCOUNT_PATH
        | USE_RING_PATH
        | USE_BOND_PATH
        | USE_HCOUNT_CLASS_PATH
        | USE_ATOM_CLASS_PATH
        | USE_SHORTCUT_LABELS;
    if which_bits & !IMPLEMENTED != 0 {
        return Err("fingerprint feature class not ported yet");
    }
    let ncounts = fp_counts.len();
    let n = mp.atoms.len();
    let mut result = 0;
    if fpflags & ACCUMULATE_BITS == 0 {
        fp_counts.fill(0);
    }
    let nbp = setup_neighbourhood(mp).ok_or("too many neighbours")?;
    let mut touched = vec![0i32; n];
    let old_bond_types: Vec<i32> = mp.bonds.iter().map(|b| b.bond_type).collect();

    let mut h_count = vec![0i32; n + 1];
    let hydrogen_bonds = |mp: &Molecule, h_count: &mut [i32]| {
        for b in &mp.bonds {
            if mp.atoms[(b.atoms[0] - 1) as usize].atom_symbol == "H" {
                h_count[b.atoms[1] as usize] += 1;
            } else if mp.atoms[(b.atoms[1] - 1) as usize].atom_symbol == "H" {
                h_count[b.atoms[0] as usize] += 1;
            }
        }
    };
    if as_query {
        hydrogen_bonds(mp, &mut h_count);
        guess_h_counts_from_substitution(mp, &nbp);
        for (i, a) in mp.atoms.iter().enumerate() {
            if a.query_h_count != NONE {
                h_count[i + 1] = a.query_h_count - ZERO_COUNT;
            }
        }
    } else {
        compute_implicit_h(mp, &mut h_count);
        hydrogen_bonds(mp, &mut h_count);
    }

    if fpflags & USE_DY_AROMATICITY != 0 {
        perceive_dy_aromaticity(mp, &nbp);
    } else {
        perceive_aromatic_bonds(mp);
    }
    let mut degree = vec![0i32; n];
    let mut cdegree = vec![0i32; n];
    let mut nspecial = vec![0i32; n];
    let mut unsaturated = vec![false; n];
    let (atom_status, bond_status) = ring_state(mp);
    let _ = &bond_status;
    set_ring_size_flags(mp, 14, &nbp);

    let mut nrare_atoms = 0;
    let mut is_rare = false;
    for i in 0..n {
        unsaturated[i] = false;
        let sym = mp.atoms[i].atom_symbol.clone();
        let mut color = number_of(&sym);
        if color <= 1 {
            color = 0;
        }
        if color > 115 {
            color = -1;
        }
        if sym == "A" {
            color = -1;
        }
        if sym == "L" {
            color = -1;
            for sl in mp.symbol_lists.iter().filter(|l| l.atom == i as i32 + 1) {
                is_rare = sl.logic != 0;
                for e in ["C", "H", "O", "N", "Cl", "F"] {
                    is_rare = is_rare && !atom_symbol_match(e, &sl.string);
                }
            }
        } else {
            is_rare = color > 0 && !atom_symbol_match(&sym, "C,H,O,N,S,P,Cl,F");
        }
        if is_rare && (exclude_atom != i as i32 + 1 || exclude_atom <= 0) {
            nrare_atoms += 1;
        }
        if which_bits & USE_SHORTCUT_LABELS != 0 && sym == "R" && !mp.atoms[i].atext.is_empty() {
            color = ((0xFFFF00 & hash_string(&mp.atoms[i].atext)) | 119) as i32;
        }
        mp.atoms[i].color = color;
    }

    let (mut ndouble, mut naromatic, mut nfusionb) = (0, 0, 0);
    for i in 0..mp.bonds.len() {
        let (a0, a1, bt) = (
            mp.bonds[i].atoms[0],
            mp.bonds[i].atoms[1],
            mp.bonds[i].bond_type,
        );
        let color = match bt {
            SINGLE => 1,
            DOUBLE => {
                if a0 != exclude_atom && a1 != exclude_atom {
                    ndouble += 1;
                }
                2
            }
            TRIPLE => 3,
            AROMATIC => {
                if a0 != exclude_atom && a1 != exclude_atom {
                    naromatic += 1;
                }
                4
            }
            _ => 0,
        };
        mp.bonds[i].color = color;
        let (i0, i1) = ((a0 - 1) as usize, (a1 - 1) as usize);
        if color > 1 {
            unsaturated[i0] = true;
            unsaturated[i1] = true;
        }
        if mp.atoms[i0].color != 0 && mp.atoms[i1].color != 0 {
            degree[i0] += 1;
            degree[i1] += 1;
        }
        if bt == DOUBLE {
            nspecial[i0] += 1;
            nspecial[i1] += 1;
        } else if bt == TRIPLE {
            nspecial[i0] += 2;
            nspecial[i1] += 2;
        }
        if bt != SINGLE {
            continue;
        }
        if mp.atoms[i0].color != 0 && mp.atoms[i1].color == 6 {
            cdegree[i0] += 1;
        }
        if mp.atoms[i1].color != 0 && mp.atoms[i0].color == 6 {
            cdegree[i1] += 1;
        }
        if (exclude_atom <= 0 || (a0 != exclude_atom && a1 != exclude_atom))
            && atom_status[i0] > 2
            && atom_status[i1] > 2
        {
            nfusionb += 1;
        }
    }
    for a in &mut mp.atoms {
        if a.color < 0 {
            a.color = 0;
        }
    }

    macro_rules! add_bit {
        ($seed:expr) => {{
            let p = hash_position($seed, ncounts);
            fp_counts[p] += 1;
        }};
    }
    macro_rules! add_bit_count {
        ($seed:expr, $count:expr) => {{
            let p = hash_position($seed, ncounts);
            fp_counts[p] += $count;
        }};
    }
    macro_rules! path_bits {
        ($seed:expr, $nbonds:expr, $minb:expr, $maxb:expr, $sprout:expr, $last:expr, $flags:expr) => {{
            let mut c = PathCtx {
                mp: &*mp,
                nbp: &nbp,
                counts: &mut *fp_counts,
                touched: &mut touched,
                minbonds: $minb,
                maxbonds: $maxb,
                flags: $flags,
                exclude_atom,
            };
            set_path_bits_rec(&mut c, $seed, $nbonds, $sprout, $last)
        }};
    }

    if which_bits & USE_ATOM_COUNT != 0 {
        const NCOUNT_HASH: usize = 128;
        // the C macro is not parenthesised, so `seed % NCOUNT_SEED_HASH` means `(seed % 128) * 128`
        const NCOUNT_SEED_HASH: usize = 128 * 128;
        let mut type_hash = [0i32; NCOUNT_HASH];
        let mut seed_hash = vec![0i32; NCOUNT_SEED_HASH];
        let (mut nringch2, mut nfusionch, mut nspiro) = (0, 0, 0);
        for i in 0..n {
            if i as i32 + 1 == exclude_atom {
                continue;
            }
            let color = mp.atoms[i].color;
            if color == 0 {
                continue;
            }
            if color == 6 {
                if atom_status[i] > 0 && h_count[i + 1] >= 2 {
                    nringch2 += 1;
                }
                if atom_status[i] > 2 && h_count[i + 1] >= 1 {
                    nfusionch += 1;
                }
                if atom_status[i] > 3 {
                    nspiro += 1;
                }
                continue;
            }
            let mut hash = 0i32;
            let mut seed = next_hash(317 * ATOM_COUNT_SEED, 507);
            seed = ns(seed, color + 17);
            seed_hash[((seed % 128) * 128) as usize] += 1;
            if (color == 7 || color == 8) && h_count[i + 1] <= 0 {
                continue;
            }
            hash = hash * 7 + color + 13;
            type_hash[hash as usize % NCOUNT_HASH] += 1;
            seed = ns(seed, color + 13);
            seed_hash[((seed % 128) * 128) as usize] += 1;
            if color != 7 && color != 8 {
                hash = hash * 7 + color + 2 * 13;
                type_hash[hash as usize % NCOUNT_HASH] += 1;
                seed = ns(seed, color + 2 * 13);
                seed_hash[((seed % 128) * 128) as usize] += 1;
            }
        }
        let mut nbits = 0;
        if ncounts <= 2048 {
            for (i, &c) in type_hash.iter().enumerate() {
                let i = i as i32;
                for (thr, mul, add) in [
                    (0, 19, 3),
                    (1, 23, 5),
                    (2, 29, 7),
                    (4, 31, 11),
                    (8, 37, 13),
                    (16, 41, 17),
                ] {
                    if c > thr {
                        add_bit_count!(ns(i64::from(i * mul) as u64, add), c);
                        nbits += 1;
                    }
                }
            }
        } else {
            for (i, &c) in seed_hash.iter().enumerate() {
                let i = i as i32;
                for (thr, mul, add) in [
                    (0, 19, 3),
                    (1, 23, 5),
                    (2, 29, 7),
                    (3, 29, 71),
                    (4, 31, 11),
                    (6, 31, 113),
                    (8, 37, 13),
                    (12, 37, 133),
                    (16, 41, 17),
                ] {
                    if c > thr {
                        add_bit_count!(ns(i64::from(i * mul) as u64, add), c);
                        nbits += 1;
                    }
                }
            }
        }
        for (cond, add, cnt) in [
            (naromatic >= 10, 341, naromatic),
            (naromatic >= 14, 441, naromatic),
            (naromatic >= 18, 541, naromatic),
            (naromatic >= 22, 641, naromatic),
            (naromatic >= 26, 741, naromatic),
            (ndouble >= 3, 371, naromatic),
            (ndouble >= 5, 471, naromatic),
            (ndouble >= 8, 571, naromatic),
            (nringch2 >= 6, 411, nringch2),
            (nringch2 >= 12, 511, nringch2),
            (nringch2 >= 22, 611, nringch2),
            (nfusionch >= 2, 421, nfusionch),
            (nfusionch >= 4, 521, nfusionch),
            (nspiro >= 1, 322, nspiro),
            (nspiro >= 1, 422, nspiro),
            (nspiro >= 2, 522, nspiro),
            (nspiro >= 2, 622, nspiro),
            (nfusionb >= 1, 425, nfusionb),
            (nfusionb >= 2, 525, nfusionb),
            (nfusionb >= 3, 625, nfusionb),
            (nfusionb >= 5, 825, nfusionb),
        ] {
            if cond {
                add_bit_count!(ns(ATOM_COUNT_SEED, add), cnt);
                nbits += 1;
            }
        }
        if nrare_atoms > 0 {
            let seed = ns(317 * ATOM_COUNT_SEED, 5);
            add_bit_count!(seed, nrare_atoms);
            add_bit_count!(ns(seed, 101), nrare_atoms);
            add_bit_count!(ns(seed, 103), nrare_atoms);
            nbits += 3;
            if nrare_atoms > 1 {
                add_bit_count!(ns(seed, 203), nrare_atoms);
                nbits += 1;
            }
            if nrare_atoms > 3 {
                add_bit_count!(ns(seed, 211), nrare_atoms);
                nbits += 1;
            }
        }
        result += nbits;
    }

    if which_bits & USE_ATOM_SYMBOL_PATH != 0 {
        let mut seed = ATOM_SYMBOL_PATH_SEED;
        for i in 0..n {
            if i as i32 + 1 == exclude_atom {
                continue;
            }
            let color = mp.atoms[i].color;
            if color <= 0 {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            seed = ns(seed, color);
            if color != 6 && color != 7 && color != 8 {
                add_bit!(seed);
                result += 1;
            }
            let sub_desc = mp.atoms[i].sub_desc;
            if color != 6
                && (!as_query
                    || sub_desc == SUB_AS_IS
                    || (sub_desc != NONE
                        && sub_desc != SUB_MORE
                        && sub_desc == degree[i] + SUB_ONE - 1))
            {
                result += path_bits!(
                    seed.wrapping_add(12347),
                    1,
                    1,
                    1,
                    i,
                    -1,
                    FORCED_HETERO_END | PROCESS_CHAINS
                );
            }
            if color >= 10
                || (color == 7 && degree[i] > 0)
                || (color == 8 && degree[i] > 1)
                || (color == 6 && degree[i] > 2 && atom_status[i] > 0)
            {
                result += path_bits!(
                    seed,
                    1,
                    1,
                    2,
                    i,
                    -1,
                    STOP_AT_HEAVY_ATOM | FORCED_HETERO_END | PROCESS_CHAINS
                );
            }
            if color > 6 {
                result += path_bits!(
                    seed.wrapping_add(i64::from(217i32.wrapping_mul(color)) as u64),
                    1,
                    3,
                    4,
                    i,
                    -1,
                    IGNORE_PATH_SYMBOL
                        | PROCESS_RING_CLOSURES
                        | STOP_AT_HEAVY_ATOM
                        | PROCESS_CHAINS
                );
                if color > 10 && color <= 18 {
                    result += path_bits!(
                        seed.wrapping_add(17),
                        1,
                        5,
                        7,
                        i,
                        -1,
                        FORCED_HETERO_END
                            | IGNORE_PATH_SYMBOL
                            | STOP_AT_HEAVY_ATOM
                            | PROCESS_CHAINS
                    );
                }
            }
            if (color == 7 || color == 8) && cdegree[i] > 2 {
                seed = ns(seed, color * 23);
                result += path_bits!(
                    seed,
                    1,
                    2,
                    6,
                    i,
                    -1,
                    STOP_AT_HEAVY_ATOM | IGNORE_TERM_SYMBOL | PROCESS_CHAINS
                );
            }
            seed = old_seed;
            if degree[i] >= 4 && color != 5 && color < 18 {
                seed = ns(seed, color);
                result += path_bits!(
                    ns(seed, 3 * 107),
                    1,
                    2,
                    4,
                    i,
                    -1,
                    STOP_AT_HEAVY_ATOM | IGNORE_TERM_SYMBOL | PROCESS_CHAINS
                );
            }
            seed = old_seed;
            if atom_status[i] >= 4 && color != 5 && color < 18 {
                seed = ns(seed, color + 55);
                result += path_bits!(
                    ns(seed, 3 * 109),
                    1,
                    2,
                    4,
                    i,
                    -1,
                    STOP_AT_HEAVY_ATOM
                        | IGNORE_TERM_SYMBOL
                        | IGNORE_PATH_SYMBOL
                        | PROCESS_RING_CLOSURES
                        | PROCESS_CHAINS
                );
            }
            for j in 0..nbp[i].atoms.len() {
                let bp = &mp.bonds[nbp[i].bonds[j]];
                let ai = nbp[i].atoms[j];
                if ai as i32 + 1 == exclude_atom || bp.color == 0 {
                    continue;
                }
                if bp.bond_type != DOUBLE && bp.bond_type != TRIPLE {
                    continue;
                }
                let bcolor = bp.color;
                seed = old_seed;
                seed = ns(seed, color);
                seed = ns(seed, bcolor * 613);
                touched[ai] = 1;
                result += path_bits!(
                    seed,
                    2,
                    5,
                    5,
                    ai,
                    i as i32,
                    STOP_AT_HEAVY_ATOM | IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES
                );
                touched[ai] = 0;
            }
            seed = old_seed;
            touched[i] = 0;
        }
    }

    if which_bits & USE_AUGMENTED_ATOM != 0 {
        let col = |mp: &Molecule, k: usize| mp.atoms[k].color;
        let bcol = |mp: &Molecule, k: usize| mp.bonds[k].color;
        for i in 0..n {
            let color = mp.atoms[i].color;
            if color <= 0 || i as i32 + 1 == exclude_atom {
                continue;
            }
            let nb = &nbp[i];
            if nspecial[i] >= 2 {
                add_bit!(ns(color as u64 * AUGMENTED_ATOM_SEED, 101));
                result += 1;
            }
            if (h_count[i + 1] > 0 || color != 6) && degree[i] >= 2 {
                for i1 in 0..nb.atoms.len() {
                    if col(mp, nb.atoms[i1]) == 0
                        || nb.atoms[i1] as i32 + 1 == exclude_atom
                        || bcol(mp, nb.bonds[i1]) == 0
                    {
                        continue;
                    }
                    for i2 in i1 + 1..nb.atoms.len() {
                        let mut seed = ns(AUGMENTED_ATOM_SEED, 97);
                        seed = ns(seed, color);
                        if col(mp, nb.atoms[i2]) == 0
                            || nb.atoms[i2] as i32 + 1 == exclude_atom
                            || bcol(mp, nb.bonds[i2]) == 0
                        {
                            continue;
                        }
                        if col(mp, nb.atoms[i1]) == 6 && col(mp, nb.atoms[i2]) == 6 {
                            continue;
                        }
                        let mut sum = 0u64;
                        sum = sum.wrapping_add(i64::from(
                            col(mp, nb.atoms[i1]) * bcol(mp, nb.bonds[i1]),
                        ) as u64);
                        sum = sum.wrapping_add(i64::from(
                            col(mp, nb.atoms[i2]) * bcol(mp, nb.bonds[i2]),
                        ) as u64);
                        let mut prod = 1u64;
                        prod = prod.wrapping_mul(i64::from(col(mp, nb.atoms[i1])) as u64) & 0xFFF;
                        prod = prod.wrapping_mul(i64::from(col(mp, nb.atoms[i2])) as u64) & 0xFFF;
                        seed = next_hash(seed, sum);
                        seed = next_hash(seed, prod);
                        add_bit!(seed);
                        result += 1;
                        let mut nmulti = 0;
                        if bcol(mp, nb.bonds[i1]) >= 2 {
                            nmulti += 1;
                        }
                        if bcol(mp, nb.bonds[i2]) >= 2 {
                            nmulti += 1;
                        }
                        let sd = mp.atoms[i].sub_desc;
                        if color != 6
                            && nmulti > 0
                            && (!as_query
                                || sd == SUB_AS_IS
                                || (sd != NONE && sd != SUB_MORE && sd == degree[i] + SUB_ONE - 1))
                        {
                            add_bit!(ns(seed, nmulti + 173 * degree[i]));
                            result += 1;
                        }
                    }
                }
            }
            if degree[i] <= 2 {
                continue;
            }
            for i1 in 0..nb.atoms.len() {
                if nb.atoms[i1] as i32 + 1 == exclude_atom {
                    continue;
                }
                for i2 in i1 + 1..nb.atoms.len() {
                    if nb.atoms[i2] as i32 + 1 == exclude_atom {
                        continue;
                    }
                    for i3 in i2 + 1..nb.atoms.len() {
                        if nb.atoms[i3] as i32 + 1 == exclude_atom {
                            continue;
                        }
                        let mut seed = AUGMENTED_ATOM_SEED;
                        seed = ns(seed, color);
                        let ids = [i1, i2, i3];
                        if ids
                            .iter()
                            .any(|&k| col(mp, nb.atoms[k]) == 0 || bcol(mp, nb.bonds[k]) == 0)
                        {
                            continue;
                        }
                        let nqtmp =
                            ids.iter().filter(|&&k| col(mp, nb.atoms[k]) != 6).count() as i32;
                        let nmulti = ids
                            .iter()
                            .map(|&k| {
                                i32::from(bcol(mp, nb.bonds[k]) == 2)
                                    + i32::from(bcol(mp, nb.bonds[k]) == 3)
                            })
                            .sum::<i32>();
                        let mut sum = 0u64;
                        let mut prod = 1u64;
                        for &k in &ids {
                            sum = sum.wrapping_add(i64::from(
                                col(mp, nb.atoms[k]) * bcol(mp, nb.bonds[k]),
                            ) as u64);
                            prod =
                                prod.wrapping_mul(i64::from(col(mp, nb.atoms[k])) as u64) & 0xFFF;
                        }
                        seed = next_hash(seed, sum);
                        seed = next_hash(seed, prod);
                        if (nqtmp > 2 || nmulti >= 2) && color == 6 {
                            add_bit!(seed);
                            add_bit!(ns(seed, 73));
                            result += 2;
                        }
                        seed = next_hash(seed, sum.wrapping_mul(prod) & 0xFFF);
                        add_bit!(seed);
                        result += 1;
                        if nmulti >= 2 || color > 6 {
                            add_bit!(ns(seed, 53));
                            result += 1;
                        }
                    }
                }
            }
        }
    }

    if which_bits & USE_AUGMENTED_BOND != 0 {
        let col = |mp: &Molecule, k: usize| i64::from(mp.atoms[k].color) as u64;
        let bcol = |mp: &Molecule, k: usize| i64::from(mp.bonds[k].color) as u64;
        for i in 0..mp.bonds.len() {
            let bp = &mp.bonds[i];
            if bp.atoms[0] == exclude_atom || bp.atoms[1] == exclude_atom {
                continue;
            }
            let (ai1, ai2) = ((bp.atoms[0] - 1) as usize, (bp.atoms[1] - 1) as usize);
            if degree[ai1] <= 2 || degree[ai2] <= 2 {
                continue;
            }
            let bond_color = bp.color;
            // sum and product over a pair of ligands of one end atom
            let pair = |mp: &Molecule, ai: usize, i1: usize, i2: usize| -> Option<(u64, u64)> {
                let nb = &nbp[ai];
                if nb.bonds[i1] == i || nb.bonds[i2] == i {
                    return None;
                }
                if nb.atoms[i1] as i32 + 1 == exclude_atom
                    || nb.atoms[i2] as i32 + 1 == exclude_atom
                {
                    return None;
                }
                if mp.atoms[nb.atoms[i1]].color == 0
                    || mp.atoms[nb.atoms[i2]].color == 0
                    || mp.bonds[nb.bonds[i1]].color == 0
                    || mp.bonds[nb.bonds[i2]].color == 0
                {
                    return None;
                }
                let mut sum = 0u64;
                sum = sum.wrapping_add(col(mp, nb.atoms[i1]).wrapping_mul(bcol(mp, nb.bonds[i1])));
                sum = sum.wrapping_add(col(mp, nb.atoms[i2]).wrapping_mul(bcol(mp, nb.bonds[i2])));
                let mut prod = 1u64;
                prod =
                    prod.wrapping_mul(col(mp, nb.atoms[i1]).wrapping_add(bcol(mp, nb.bonds[i1])));
                prod =
                    prod.wrapping_mul(col(mp, nb.atoms[i2]).wrapping_add(bcol(mp, nb.bonds[i2])));
                Some((sum & 0x0FFF, prod & 0x0FFF))
            };
            for i1 in 0..nbp[ai1].atoms.len() {
                for i2 in i1 + 1..nbp[ai1].atoms.len() {
                    let Some((sumi, prodi)) = pair(mp, ai1, i1, i2) else {
                        continue;
                    };
                    for j1 in 0..nbp[ai2].atoms.len() {
                        for j2 in j1 + 1..nbp[ai2].atoms.len() {
                            let Some((sumj, prodj)) = pair(mp, ai2, j1, j2) else {
                                continue;
                            };
                            let mut seed = AUGMENTED_BOND_SEED;
                            seed = ns(seed, bond_color);
                            seed = next_hash(seed, prodi + prodj);
                            seed = next_hash(seed, sumi.wrapping_mul(sumj));
                            add_bit!(seed);
                            result += 1;
                        }
                    }
                }
            }
        }
    }

    if which_bits & USE_HCOUNT_PAIR != 0 {
        for i in 0..mp.bonds.len() {
            let bp = &mp.bonds[i];
            let (a0, a1) = (bp.atoms[0] as usize, bp.atoms[1] as usize);
            if h_count[a0] == 0 && h_count[a1] == 0 {
                continue;
            }
            if bp.atoms[0] == exclude_atom || bp.atoms[1] == exclude_atom {
                continue;
            }
            let (c0, c1) = (mp.atoms[a0 - 1].color, mp.atoms[a1 - 1].color);
            if c0 == 6 && c1 == 6 && bp.color == 1 {
                continue;
            }
            if c0 == 0 || c1 == 0 {
                continue;
            }
            for j1 in 0..=h_count[a0] {
                for j2 in 0..=h_count[a1] {
                    if j1 + j2 == 0 {
                        continue;
                    }
                    if !unsaturated[a0 - 1] && !unsaturated[a1 - 1] && j1 * j2 == 0 && bp.color <= 1
                    {
                        continue;
                    }
                    if j1 + j2 > 3 {
                        continue;
                    }
                    let mut seed = ns(HCOUNT_PAIR_SEED, 7 * (j1 + 1) * (j2 + 1));
                    seed = ns(seed, 53 * (j1 + j2));
                    seed = ns(seed, bp.color);
                    seed = ns(seed, c0 + c1);
                    seed = ns(seed, c0 * c1);
                    add_bit!(seed);
                    result += 1;
                    if bp.color > 1 {
                        seed = ns(seed, 83);
                        add_bit!(seed);
                        result += 1;
                        if j1 + j2 == 1 && (c0 != 6 || c1 != 6) && bp.color <= 3 {
                            seed = ns(seed, 91);
                            add_bit!(seed);
                            seed = ns(seed, 97);
                            add_bit!(seed);
                            seed = ns(seed, 103);
                            add_bit!(seed);
                            result += 3;
                        }
                    }
                }
            }
        }
    }

    if which_bits & USE_HCOUNT_PATH != 0 {
        let mut seed = HCOUNT_PATH_SEED;
        for i in 0..n {
            let color = mp.atoms[i].color;
            if color <= 0 || i as i32 + 1 == exclude_atom || h_count[i + 1] == 0 {
                continue;
            }
            if color == 6 && h_count[i + 1] < 3 && degree[i] < 3 {
                continue;
            }
            touched[i] = 1;
            let mut old_seed = seed;
            seed = ns(seed, color);
            if color != 6 {
                result += path_bits!(seed, 1, 2, 5, i, -1, IGNORE_PATH_SYMBOL | PROCESS_CHAINS);
            } else {
                if degree[i] > 2 {
                    result += path_bits!(
                        ns(seed, 103),
                        1,
                        2,
                        5,
                        i,
                        -1,
                        IGNORE_PATH_SYMBOL | FORCED_HETERO_END | PROCESS_CHAINS
                    );
                }
                if h_count[i + 1] >= 3 {
                    result += path_bits!(
                        ns(seed, 1105),
                        1,
                        3,
                        6,
                        i,
                        -1,
                        IGNORE_PATH_SYMBOL | FORCED_HETERO_END | PROCESS_CHAINS
                    );
                }
            }
            touched[i] = 0;
            seed = old_seed;
            if color == 6 {
                continue;
            }
            if h_count[i + 1] > 1 {
                touched[i] = 1;
                old_seed = seed;
                seed = ns(seed, 113);
                seed = ns(seed, color);
                result += path_bits!(seed, 1, 1, 5, i, -1, IGNORE_PATH_SYMBOL | PROCESS_CHAINS);
                touched[i] = 0;
                seed = old_seed;
            }
            for j in 1..h_count[i + 1] {
                seed = ns(seed, 61 * j);
                seed = ns(seed, color);
                add_bit!(seed);
                result += 1;
            }
            seed = old_seed;
        }
    }

    // ring paths: only ring atoms and bonds take part
    for i in 0..n {
        if atom_status[i] <= 0 || i as i32 + 1 == exclude_atom {
            mp.atoms[i].color = 0;
        }
    }
    for i in 0..mp.bonds.len() {
        let (a0, a1) = (mp.bonds[i].atoms[0], mp.bonds[i].atoms[1]);
        if bond_status[i] <= 0
            && atom_status[(a0 - 1) as usize] == 0
            && atom_status[(a1 - 1) as usize] == 0
        {
            mp.bonds[i].color = 0;
        }
        if a0 == exclude_atom || a1 == exclude_atom {
            mp.bonds[i].color = 0;
        }
    }
    if which_bits & USE_RING_PATH != 0 {
        let mut seed = RING_PATH_SEED;
        for i in 0..n {
            let color = mp.atoms[i].color;
            if color <= 0 {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            seed = ns(seed, color);
            if color > 5 && color < 10 && atom_status[i] > 2 {
                seed = ns(seed, 61);
                result += path_bits!(
                    seed,
                    1,
                    3,
                    8,
                    i,
                    -1,
                    STOP_AT_HEAVY_ATOM | PROCESS_RING_CLOSURES
                );
            }
            seed = old_seed;
            touched[i] = 0;
        }
    }

    // all atoms but hydrogens in one class, bond types retained
    for i in 0..n {
        let c = number_of(&mp.atoms[i].atom_symbol);
        mp.atoms[i].color = if c == 1 || i as i32 + 1 == exclude_atom {
            0
        } else {
            ANY_COLOR
        };
    }
    for b in &mut mp.bonds {
        b.color = match b.bond_type {
            SINGLE => 1,
            DOUBLE => 2,
            TRIPLE => 3,
            AROMATIC => 4,
            _ => 0,
        };
        if b.atoms[0] == exclude_atom || b.atoms[1] == exclude_atom {
            b.color = 0;
        }
    }
    if which_bits & USE_BOND_PATH != 0 {
        let mut seed = BOND_PATH_SEED;
        for i in 0..n {
            if mp.atoms[i].color <= 0 || i as i32 + 1 == exclude_atom {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            if degree[i] > 2 && atom_status[i] > 1 {
                result += path_bits!(seed, 1, 4, 4, i, -1, PROCESS_CHAINS);
                if mp.atoms[i].rsize_flags & SPECIAL_RING != 0 {
                    result += path_bits!(
                        ns(seed, 217),
                        1,
                        5,
                        5,
                        i,
                        -1,
                        IGNORE_PATH_SYMBOL | PROCESS_CHAINS
                    );
                }
            }
            seed = old_seed;
            seed = ns(seed, 11);
            result += path_bits!(
                seed,
                1,
                4,
                6,
                i,
                -1,
                FORCED_RING_PATH | IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES
            );
            for j in 0..nbp[i].atoms.len() {
                let bp = &mp.bonds[nbp[i].bonds[j]];
                let ai = nbp[i].atoms[j];
                if ai as i32 + 1 == exclude_atom || bp.color == 0 {
                    continue;
                }
                if bp.bond_type != DOUBLE && bp.bond_type != TRIPLE {
                    continue;
                }
                if atom_status[i] <= 0 && bp.bond_type != TRIPLE {
                    continue;
                }
                seed = old_seed;
                seed = ns(seed, bp.color * 413);
                touched[ai] = 1;
                result += path_bits!(
                    seed,
                    2,
                    4,
                    5,
                    ai,
                    i as i32,
                    IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES | PROCESS_CHAINS
                );
                touched[ai] = 0;
            }
            seed = old_seed;
            touched[i] = 0;
        }
    }

    // atom type classes: carbon, hetero atoms, nothing
    for i in 0..n {
        if i as i32 + 1 == exclude_atom {
            mp.atoms[i].color = 0;
            continue;
        }
        let sym = mp.atoms[i].atom_symbol.clone();
        let mut color;
        if sym == "L" {
            color = 0;
            for sl in mp.symbol_lists.iter().filter(|l| l.atom == i as i32 + 1) {
                color = 0;
                if sl.logic == super::reaccs::EXCLUSIVE {
                    break;
                }
                for tok in sl.string.split(',').filter(|t| !t.is_empty()) {
                    color |= match tok {
                        "H" => 1,
                        "C" => 2,
                        "O" | "N" | "S" => 4,
                        _ => 8,
                    };
                }
                color = match color {
                    2 => 6,
                    4 | 8 => 8,
                    _ => 0,
                };
            }
        } else {
            color = match sym.as_str() {
                "H" | "D" | "T" | "A" => 0,
                "C" => 6,
                "N" | "O" | "S" | "Q" => 8,
                _ => {
                    let t = number_of(&sym);
                    if 1 < t && t < 115 { 8 } else { 0 }
                }
            };
        }
        if color > 115 {
            color = 0;
        }
        mp.atoms[i].color = color;
    }
    if which_bits & USE_HCOUNT_CLASS_PATH != 0 {
        let mut seed = HCOUNT_CLASS_PATH_SEED;
        for i in 0..n {
            let color = mp.atoms[i].color;
            if color <= 0 || i as i32 + 1 == exclude_atom || h_count[i + 1] == 0 {
                continue;
            }
            if color == 6 && h_count[i + 1] < 2 {
                continue;
            }
            let old_seed = seed;
            seed = ns(seed, color);
            touched[i] = 1;
            if color == 6 {
                result += path_bits!(seed, 1, 2, 4, i, -1, FORCED_HETERO_END | PROCESS_CHAINS);
            } else {
                result += path_bits!(
                    seed,
                    1,
                    2,
                    5,
                    i,
                    -1,
                    IGNORE_PATH_SYMBOL | FORCED_HETERO_END | PROCESS_CHAINS
                );
            }
            seed = old_seed;
            touched[i] = 0;
        }
    }
    if which_bits & USE_ATOM_CLASS_PATH != 0 {
        let mut seed = ns(ATOM_CLASS_PATH_SEED, 117);
        for i in 0..n {
            let color = mp.atoms[i].color;
            if color <= 0 || i as i32 + 1 == exclude_atom || (color == 6 && degree[i] < 3) {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            seed = ns(seed, color);
            if color != 6 {
                result += path_bits!(
                    seed,
                    1,
                    3,
                    4,
                    i,
                    -1,
                    FORCED_RING_PATH | PROCESS_RING_CLOSURES | PROCESS_CHAINS
                );
            }
            touched[i] = 0;
            seed = old_seed;
        }
    }
    // a single bond class
    for b in &mut mp.bonds {
        b.color = if (SINGLE..=8).contains(&b.bond_type) {
            5
        } else {
            0
        };
        if b.atoms[0] == exclude_atom || b.atoms[1] == exclude_atom {
            b.color = 0;
        }
    }
    if which_bits & USE_ATOM_CLASS_PATH != 0 {
        let mut seed = ATOM_CLASS_PATH_SEED;
        for i in 0..n {
            let color = mp.atoms[i].color;
            if color <= 0 || i as i32 + 1 == exclude_atom || color == 6 {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            seed = ns(seed, color);
            seed = ns(seed, 23 + color * 19);
            result += path_bits!(
                seed,
                1,
                3,
                9,
                i,
                -1,
                IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES
            );
            seed = old_seed;
            touched[i] = 0;
        }
        let (mut qq_count, mut qc_count) = (0i32, 0i32);
        for i in 0..mp.bonds.len() {
            let bp = &mp.bonds[i];
            if bond_status[i] == 0
                || bp.color == 0
                || bp.atoms[0] == exclude_atom
                || bp.atoms[1] == exclude_atom
            {
                continue;
            }
            let c1 = mp.atoms[(bp.atoms[0] - 1) as usize].color;
            let c2 = mp.atoms[(bp.atoms[1] - 1) as usize].color;
            if c1 == 0 || c2 == 0 || (c1 == 6 && c2 == 6) {
                continue;
            }
            let mul = if c1 != 6 && c2 != 6 {
                qq_count += 1;
                17
            } else {
                qc_count += 1;
                19
            };
            for j in 3..9 {
                if bp.rsize_flags & (1 << j) != 0 {
                    add_bit!(ns(ATOM_CLASS_PATH_SEED * mul, j * 8));
                    result += 1;
                }
            }
        }
        let mut seed = 2 * ATOM_CLASS_PATH_SEED + 3;
        let mut i = 1;
        while i <= qq_count {
            seed = ns(seed, i * 153);
            add_bit!(seed);
            result += 1;
            seed = ns(seed, 53);
            if i <= 1 {
                add_bit!(seed);
                result += 1;
            }
            i = (1.0 + f64::from(i) * 1.5) as i32;
        }
        seed = 3 * ATOM_CLASS_PATH_SEED + 5;
        for i in 1..=qc_count.min(2) {
            seed = ns(seed, i * 157);
            add_bit!(seed);
            result += 1;
        }
        let mut i = 3;
        while i <= qc_count {
            seed = ns(seed, i * 157);
            add_bit!(seed);
            result += 1;
            i = (f64::from(i) * 1.8) as i32;
        }
    }

    // ring pattern preparation: only bonds with a ring atom, atom classes C and Q
    for i in 0..mp.bonds.len() {
        let (a0, a1) = (mp.bonds[i].atoms[0], mp.bonds[i].atoms[1]);
        if a0 == exclude_atom || a1 == exclude_atom {
            continue;
        }
        mp.bonds[i].color = 5;
        let (i0, i1) = ((a0 - 1) as usize, (a1 - 1) as usize);
        if atom_status[i0] == 0 && atom_status[i1] == 0 {
            mp.bonds[i].color = 0;
        } else {
            for k in [i0, i1] {
                let ap = &mut mp.atoms[k];
                if ap.color != 0 {
                    if ap.color != 6 {
                        ap.color = 8;
                    }
                    if ap.rsize_flags == 0 {
                        ap.color = 0;
                    }
                }
            }
        }
    }

    for (b, t) in mp.bonds.iter_mut().zip(old_bond_types) {
        b.bond_type = t;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::reaccsio::mol_str_to_mol;

    #[test]
    fn hash_functions() {
        // values of the reference implementation for a few inputs
        assert_eq!(next_hash(0, 0), 0);
        assert!(hash_position(next_hash(1001, 65), 1024) < 1024);
    }

    /// the hashed feature counts must be those of the Avalon C code
    #[test]
    fn counts_match_c() {
        let text = include_str!("../../../../../tests/data/avalon_fp.jsonl");
        let mut bad = Vec::new();
        let mut n = 0;
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().expect("input");
            let which = v["which"].as_i64().expect("which") as i32;
            let query = v["query"].as_i64().expect("query") != 0;
            let fpflags = v["fpflags"].as_i64().expect("fpflags") as i32;
            let nbits = v["nbits"].as_u64().expect("nbits") as usize;
            let want = v["output"].as_str().expect("output");
            let mut m = mol_str_to_mol(input).expect("molecule");
            let mut counts = vec![0; nbits];
            let r =
                set_fingerprint_counts_with_focus(&mut m, &mut counts, which, query, fpflags, 0)
                    .expect("ok");
            let mut got = format!("r {r}\n");
            for (i, c) in counts.iter().enumerate() {
                if *c != 0 {
                    got += &format!("{i}:{c} ");
                }
            }
            got.push('\n');
            n += 1;
            if got != want {
                bad.push(format!("which {which:x} query {query} fpflags {fpflags} nbits {nbits}\n{input}\nwant {want}\ngot  {got}"));
            }
        }
        assert!(
            bad.is_empty(),
            "{} of {n} differ; first:\n{}\n{}",
            bad.len(),
            bad[0],
            bad.iter()
                .map(|b| b.lines().next().unwrap_or("").to_string())
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
}
