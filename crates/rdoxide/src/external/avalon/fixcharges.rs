//! Port of `fixcharges.c`: fixing the charges of molecules with a net charge by removing protons from the
//! most acidic atoms (`InitializeChargeDataTables`, `SetpKaValues`, `RechargeMolecule`).
//! Copyright 2010 Novartis Institutes for BioMedical Research Inc.
//!
//! The tables are global as in C. The `charge_log` support for parameter fitting is not ported.

use std::sync::Mutex;

use super::pattern::{ANY_RADICAL, AugmentedAtom, Ligand, aa_match, string_to_augmented_atom};
use super::perceive::{Neighbourhood, SINGLE, setup_neighbourhood};
use super::reaccs::{ANY_CHARGE, Molecule, NONE};
use super::reaccsio::{Lines, Status};
use super::scan::{Scan, fmt_gp};
use super::stereo::{EVEN_PARITY, ODD_PARITY, add_msg, atom_parity, cis_trans_perception};
use super::utilities::{atom_symbol_match, reset_colors};

/// `ANY_BOND`
const ANY_BOND: i32 = 8;
/// `MAXACIDIC`, `MAXCHARGEINC`, `MAXATOMACIDITY`, `MAXALPHAPATH`, `MAXBETAPATH`
const MAXACIDIC: i32 = 45;

/// `struct inc_entry_t`
#[derive(Debug, Clone, Default)]
struct IncEntry {
    atom_symbol: String,
    local_inc: f64,
    alpha_inc: f64,
    beta_inc: f64,
    mult_inc: f64,
}

/// `struct path_entry_t`
#[derive(Debug, Clone, Default)]
struct PathEntry {
    path: AugmentedAtom,
    cond: f64,
}

/// the own variables of `fixcharges.c`
#[derive(Debug, Clone, Default)]
struct Tables {
    acidic_atoms: Vec<AugmentedAtom>,
    charge_inc: Vec<IncEntry>,
    atom_acidity: Vec<IncEntry>,
    alpha_path: Vec<PathEntry>,
    beta_path: Vec<PathEntry>,
    alpha: f64,
    beta: f64,
    elneg: Vec<(String, f64)>,
    pka_version: String,
    acidity_limit: f64,
}

impl Tables {
    const fn new() -> Self {
        Tables {
            acidic_atoms: Vec::new(),
            charge_inc: Vec::new(),
            atom_acidity: Vec::new(),
            alpha_path: Vec::new(),
            beta_path: Vec::new(),
            alpha: 0.0,
            beta: 0.0,
            elneg: Vec::new(),
            pka_version: String::new(),
            acidity_limit: 24.0,
        }
    }
}

static TABLES: Mutex<Tables> = Mutex::new(Tables::new());

fn tables() -> std::sync::MutexGuard<'static, Tables> {
    TABLES.lock().unwrap_or_else(|e| e.into_inner())
}

/// `SetAcidityLimit`
pub fn set_acidity_limit(limit: f64) {
    tables().acidity_limit = limit;
}

/// `pKa_version`
pub fn pka_version() -> String {
    tables().pka_version.clone()
}

/// `SkipCommentLines`
fn skip_comment_lines(fp: &mut Lines) -> Status {
    while fp.status == Status::Normal && fp.buffer.starts_with('#') {
        fp.get_buffer();
    }
    fp.status
}

/// the augmented atom between the quotes of the current line, with `NONE` ligand charges widened to any
fn quoted_augmented_atom(fp: &Lines) -> Option<(AugmentedAtom, usize)> {
    let start = fp.buffer.find('"')? + 1;
    let len = fp.buffer[start..].find('"')?;
    let mut aa = string_to_augmented_atom(&fp.buffer[start..start + len])?;
    for l in &mut aa.ligands {
        if l.charge == NONE {
            l.charge = ANY_CHARGE;
        }
    }
    Some((aa, start + len + 1))
}

fn read_count(fp: &Lines) -> Option<i32> {
    Scan::new(fp.buffer.as_bytes()).int(0)
}

/// `ReadAcidicAtoms`
fn read_acidic_atoms(fp: &mut Lines, t: &mut Tables) -> bool {
    if skip_comment_lines(fp) != Status::Normal {
        return false;
    }
    let n = read_count(fp).unwrap_or(0);
    fp.get_buffer();
    if !(1..=MAXACIDIC).contains(&n) {
        return false;
    }
    t.acidic_atoms.clear();
    for _ in 0..n {
        if fp.status != Status::Normal || skip_comment_lines(fp) != Status::Normal {
            return false;
        }
        let Some((aa, _)) = quoted_augmented_atom(fp) else {
            return false;
        };
        t.acidic_atoms.push(aa);
        fp.get_buffer();
    }
    true
}

/// `ReadChargeIncrements` (`mult` columns present) and `ReadAtomAcidityIncrements`
fn read_increments(fp: &mut Lines, with_mult: bool) -> Option<Vec<IncEntry>> {
    if skip_comment_lines(fp) != Status::Normal {
        return None;
    }
    let n = read_count(fp).filter(|&n| n >= 0)?;
    fp.get_buffer();
    let mut out = Vec::new();
    for _ in 0..n {
        if skip_comment_lines(fp) != Status::Normal {
            break;
        }
        let mut sc = Scan::new(fp.buffer.as_bytes());
        let mut e = IncEntry {
            atom_symbol: sc.word()?,
            local_inc: sc.double(0)?,
            alpha_inc: sc.double(0)?,
            beta_inc: sc.double(0)?,
            mult_inc: 0.0,
        };
        if with_mult {
            e.mult_inc = sc.double(0)?;
        }
        out.push(e);
        fp.get_buffer();
    }
    matches!(fp.status, Status::Normal | Status::Eof).then_some(out)
}

/// `ReadAlphaPathTable` and `ReadBetaPathTable`
fn read_path_table(fp: &mut Lines, max: i32) -> Option<Vec<PathEntry>> {
    if skip_comment_lines(fp) != Status::Normal {
        return None;
    }
    let n = read_count(fp).unwrap_or(-1);
    fp.get_buffer();
    if !(0..=max).contains(&n) {
        return None;
    }
    let mut out = Vec::new();
    for _ in 0..n {
        if fp.status != Status::Normal || skip_comment_lines(fp) != Status::Normal {
            break;
        }
        let (aa, rest) = quoted_augmented_atom(fp)?;
        let cond = Scan::new(&fp.buffer.as_bytes()[rest..]).double(0)?;
        out.push(PathEntry { path: aa, cond });
        fp.get_buffer();
    }
    Some(out)
}

/// `ReadTransformation`
fn read_transformation(fp: &mut Lines) -> Option<(f64, f64)> {
    if skip_comment_lines(fp) != Status::Normal {
        return None;
    }
    let mut sc = Scan::new(fp.buffer.as_bytes());
    let r = (sc.double(0)?, sc.double(0)?);
    fp.get_buffer();
    Some(r)
}

/// `ReadElNegTable`: the number of entries read, 0 on failure
fn read_elneg_table(fp: &mut Lines, out: &mut Vec<(String, f64)>) -> usize {
    if skip_comment_lines(fp) != Status::Normal {
        return 0;
    }
    let Some(n) = read_count(fp) else {
        return 0;
    };
    fp.get_buffer();
    out.clear();
    for i in 0..n.max(0) as usize {
        if skip_comment_lines(fp) != Status::Normal {
            break;
        }
        let mut sc = Scan::new(fp.buffer.as_bytes());
        match (sc.word(), sc.double(0)) {
            (Some(s), Some(v)) => out.push((s, v)),
            _ => return i,
        }
        fp.get_buffer();
    }
    if matches!(fp.status, Status::Normal | Status::Eof) {
        n.max(0) as usize
    } else {
        0
    }
}

/// `InitializeChargeDataTables`: reads the acidity table (`.aci` file) from `text`; the C code exits the
/// process on a malformed table, here `false` is returned
pub fn initialize_charge_data_tables(text: &str) -> bool {
    let mut fp = Lines::open(text);
    let mut t = Tables::new();
    t.acidity_limit = tables().acidity_limit;
    if let Some(cp) = fp.buffer.find('_') {
        let rest = &fp.buffer[cp + 1..];
        let v = match rest.find('_') {
            Some(cph) => &rest[..cph],
            None => &rest[..rest.len().min(9)],
        };
        t.pka_version = v.to_string();
    }
    if !read_acidic_atoms(&mut fp, &mut t) {
        eprintln!("Error reading acidic atoms");
        return false;
    }
    match read_increments(&mut fp, true) {
        Some(v) => t.charge_inc = v,
        None => {
            eprintln!("Error reading charge increments");
            return false;
        }
    }
    match read_increments(&mut fp, false) {
        Some(v) => t.atom_acidity = v,
        None => {
            eprintln!("Error reading atom acidity increments");
            return false;
        }
    }
    match read_path_table(&mut fp, 10) {
        Some(v) => t.alpha_path = v,
        None => {
            eprintln!("Error reading alpha path table");
            return false;
        }
    }
    match read_path_table(&mut fp, 20) {
        Some(v) => t.beta_path = v,
        None => {
            eprintln!("Error reading beta path table");
            return false;
        }
    }
    match read_transformation(&mut fp) {
        Some((a, b)) => {
            t.alpha = a;
            t.beta = b;
        }
        None => {
            eprintln!("Error reading transformation polynom");
            return false;
        }
    }
    let n = read_elneg_table(&mut fp, &mut t.elneg);
    *tables() = t;
    n != 0
}

/// `LigandMatches`
fn ligand_matches(ap: &super::reaccs::Atom, bond_type: i32, lp: &Ligand, use_charge: bool) -> bool {
    if lp.bond_type != ANY_BOND && bond_type != lp.bond_type {
        return false;
    }
    if lp.radical != ANY_RADICAL && ap.radical != lp.radical {
        return false;
    }
    if (lp.charge != NONE || use_charge) && lp.charge != ANY_CHARGE && ap.charge != lp.charge {
        return false;
    }
    atom_symbol_match(&ap.atom_symbol, &lp.atom_symbol)
}

fn find_inc(table: &[IncEntry], symbol: &str) -> usize {
    table
        .iter()
        .position(|e| atom_symbol_match(symbol, &e.atom_symbol))
        .unwrap_or(table.len())
}

fn name10(mp: &Molecule) -> String {
    format!("{:>10}", mp.name)
}

/// `ResetValues`
fn reset_values(mp: &mut Molecule) {
    for a in &mut mp.atoms {
        a.value = 0.0;
    }
    for b in &mut mp.bonds {
        b.value = 0.0;
    }
}

/// `+=` on a `float` field: the sum is computed in double
fn add_f32(v: &mut f32, x: f64) {
    *v = (f64::from(*v) + x) as f32;
}

/// `SetpKaValues`: estimates the pKa value of the acidic atoms (`color` > 0) and stores it in `value`
fn set_pka_values(t: &Tables, mp: &mut Molecule, msgs: &mut Vec<String>) -> bool {
    let mut result = true;
    let Some(nb) = setup_neighbourhood(mp) else {
        return false;
    };
    reset_colors(mp);
    reset_values(mp);
    let name = name10(mp);
    for i in 0..mp.atoms.len() {
        if !t
            .acidic_atoms
            .iter()
            .any(|aa| aa_match(mp, i, aa, None, &nb).is_some())
        {
            continue;
        }
        mp.atoms[i].color += 1;
        mp.atoms[i].value = 0.0;
        let (ap_symbol, ap_charge) = (mp.atoms[i].atom_symbol.clone(), mp.atoms[i].charge);
        let nct = t.charge_inc.len();
        // `cip` keeps pointing into the table between the steps, as in C
        let mut cip = find_inc(&t.charge_inc, &ap_symbol);
        if cip < nct {
            add_f32(
                &mut mp.atoms[i].value,
                f64::from(ap_charge) * t.charge_inc[cip].local_inc,
            );
        }
        let ip = find_inc(&t.atom_acidity, &ap_symbol);
        if ip < t.atom_acidity.len() {
            add_f32(&mut mp.atoms[i].value, t.atom_acidity[ip].local_inc);
        }
        let nbp: &Neighbourhood = &nb[i];
        for j in 0..nbp.atoms.len() {
            let aap = mp.atoms[nbp.atoms[j]].clone();
            let abt = mp.bonds[nbp.bonds[j]].bond_type;
            let ppa = t.alpha_path.iter().position(|p| {
                atom_symbol_match(&ap_symbol, &p.path.atom_symbol)
                    && ligand_matches(&aap, abt, &p.path.ligands[0], false)
            });
            let Some(ppa) = ppa else {
                add_msg(msgs, format!("{name}: no alpha path for atom {}", i + 1));
                result = false;
                break;
            };
            if abt != SINGLE {
                let mult = t.charge_inc.get(cip).map_or(0.0, |e| e.mult_inc);
                add_f32(&mut mp.atoms[i].value, f64::from(abt - SINGLE) * mult);
            }
            if aap.charge != 0 {
                cip = find_inc(&t.charge_inc, &aap.atom_symbol);
                if cip == nct {
                    add_msg(
                        msgs,
                        format!("{name}: no alpha increment for atom {}", i + 1),
                    );
                    result = false;
                    break;
                }
                add_f32(
                    &mut mp.atoms[i].value,
                    f64::from(aap.charge) * t.charge_inc[cip].alpha_inc,
                );
            }
            let ip = find_inc(&t.atom_acidity, &aap.atom_symbol);
            if ip == t.atom_acidity.len() {
                add_msg(
                    msgs,
                    format!("{name}: no alpha increment for atom {}", i + 1),
                );
                result = false;
                break;
            }
            add_f32(
                &mut mp.atoms[i].value,
                t.alpha_path[ppa].cond * t.atom_acidity[ip].alpha_inc,
            );
            let nbph = &nb[nbp.atoms[j]];
            for k in 0..nbph.atoms.len() {
                if nbph.atoms[k] == i {
                    continue;
                }
                let bap = mp.atoms[nbph.atoms[k]].clone();
                let bbt = mp.bonds[nbph.bonds[k]].bond_type;
                let ppb = t.beta_path.iter().position(|p| {
                    atom_symbol_match(&ap_symbol, &p.path.atom_symbol)
                        && ligand_matches(&aap, abt, &p.path.ligands[0], false)
                        && ligand_matches(&bap, bbt, &p.path.ligands[1], false)
                });
                let Some(ppb) = ppb else {
                    add_msg(
                        msgs,
                        format!("{name}: no beta increment for atom {}", i + 1),
                    );
                    result = false;
                    break;
                };
                let ip = find_inc(&t.atom_acidity, &bap.atom_symbol);
                if ip == t.atom_acidity.len() {
                    add_msg(
                        msgs,
                        format!("{name}: no beta increment for atom {}", i + 1),
                    );
                    result = false;
                    break;
                }
                add_f32(
                    &mut mp.atoms[i].value,
                    t.beta_path[ppb].cond * t.atom_acidity[ip].beta_inc,
                );
                if bap.charge != 0 {
                    cip = find_inc(&t.charge_inc, &bap.atom_symbol);
                    if cip == nct {
                        add_msg(
                            msgs,
                            format!("{name}: no beta increment for atom {}", i + 1),
                        );
                        result = false;
                        break;
                    }
                    add_f32(
                        &mut mp.atoms[i].value,
                        f64::from(bap.charge) * t.charge_inc[cip].beta_inc,
                    );
                }
            }
        }
        let v = f64::from(mp.atoms[i].value);
        mp.atoms[i].value = (7.0
            + t.beta * (v - 7.0)
            + t.alpha * (v - 7.0) * t.alpha * (v - 7.0) * t.alpha * (v - 7.0))
            as f32;
    }
    result
}

/// `MarkMostAcidicAtoms`: keeps only the marks of the atoms with the lowest pKa; returns their number,
/// the minimum pKa and the gap to the next acidic atoms
fn mark_most_acidic_atoms(mp: &mut Molecule) -> (i32, f64, f64) {
    let epsilon = 0.000001;
    let mut min_pka = 1000.0;
    for a in &mp.atoms {
        if a.color != NONE && f64::from(a.value) < min_pka {
            min_pka = f64::from(a.value);
        }
    }
    let mut next_pka = 1000.0;
    let mut result = 0;
    for a in &mut mp.atoms {
        if a.color != NONE && f64::from(a.value) < min_pka + epsilon {
            result += 1;
            a.color = NONE + 1;
        } else {
            if a.color != NONE && f64::from(a.value) < next_pka {
                next_pka = f64::from(a.value);
            }
            a.color = NONE;
            a.value = 0.0;
        }
    }
    (result, min_pka, next_pka - min_pka)
}

/// `DecrementMarkedCharges`
fn decrement_marked_charges(mp: &mut Molecule) {
    for a in &mut mp.atoms {
        if a.color != NONE {
            a.charge -= 1;
        }
    }
}

/// `TotalCharge`
pub fn total_charge(mp: &Molecule) -> i32 {
    mp.atoms.iter().map(|a| a.charge).sum()
}

/// `struct atom_rank_t`
#[derive(Clone, Default)]
struct AtomRank {
    atom: i32,
    rank: i32,
    n_ligands: i32,
    elneg: f64,
    rank_sum: i32,
}

/// the insertion sort of the C code: `j` moves down while `swap(cur, prev)` holds
fn insertion_sort(a: &mut [AtomRank], swap: impl Fn(&AtomRank, &AtomRank) -> bool) {
    for i in 1..a.len() {
        for j in (1..=i).rev() {
            if swap(&a[j], &a[j - 1]) {
                a.swap(j, j - 1);
            } else {
                break;
            }
        }
    }
}

/// `AllCentersRefined`: all defined stereo centers have four different ligand numbers
fn all_centers_refined(mp: &Molecule, numbering: &[i32], msgs: &mut Vec<String>) -> bool {
    let Some(nb) = setup_neighbourhood(mp) else {
        return true;
    };
    for (i, nbp) in nb.iter().enumerate() {
        let n = nbp.atoms.len();
        let sym = mp.atoms[i].atom_symbol.as_str();
        if !(matches!(sym, "C" | "S" | "N" | "O" | "P" | "Si") && n > 2 && n <= 4) {
            continue;
        }
        let mut err = "";
        let parity = atom_parity(mp, i + 1, nbp, &mut err);
        if parity == EVEN_PARITY || parity == ODD_PARITY {
            let mut refnum: Vec<i32> = nbp.atoms.iter().map(|&a| numbering[a]).collect();
            refnum.sort_unstable();
            if refnum.windows(2).any(|w| w[0] == w[1]) {
                add_msg(
                    msgs,
                    format!("{}    unrefined parity at atom {}", name10(mp), i + 1),
                );
                return false;
            }
        }
    }
    true
}

/// `RefineAcidicAtoms`: restricts the marked atoms to the class of lowest rank; returns its size (-1 if an
/// atom symbol has no electronegativity)
fn refine_acidic_atoms(
    t: &Tables,
    mp: &mut Molecule,
    numbering: &mut [i32],
    msgs: &mut Vec<String>,
) -> i32 {
    let epsilon = 0.000001;
    let na = mp.atoms.len();
    let mut ranks = vec![AtomRank::default(); na];
    for i in 0..na {
        ranks[i].atom = i as i32 + 1;
        let Some((_, v)) = t.elneg.iter().find(|(s, _)| *s == mp.atoms[i].atom_symbol) else {
            eprintln!(
                "atom symbol '{}' not in periodic table",
                mp.atoms[i].atom_symbol
            );
            return -1;
        };
        ranks[i].elneg = *v + 3.0 * f64::from(mp.atoms[i].charge) - t.elneg[0].1;
    }
    // preliminary ranking on the electronegativity
    insertion_sort(&mut ranks, |a, b| a.elneg > b.elneg + epsilon);
    let mut j = 0;
    for i in 1..na {
        if ranks[i].elneg < ranks[i - 1].elneg - epsilon {
            j = i as i32;
        }
        ranks[i].rank = j;
    }
    insertion_sort(&mut ranks, |a, b| a.atom < b.atom);
    // unsaturation splits the ranks (the rank sum is misused)
    for r in &mut ranks {
        r.rank_sum = 0;
    }
    for b in &mp.bonds {
        if b.bond_type != SINGLE {
            ranks[(b.atoms[0] - 1) as usize].rank_sum += 1;
            ranks[(b.atoms[1] - 1) as usize].rank_sum += 1;
        }
    }
    insertion_sort(&mut ranks, |a, b| {
        a.rank < b.rank || (a.rank == b.rank && a.rank_sum > b.rank_sum)
    });
    let mut j = 0;
    for i in 1..na {
        if ranks[i].rank > ranks[i - 1].rank || ranks[i].rank_sum < ranks[i - 1].rank_sum {
            j = i as i32;
        }
        ranks[i].rank = j;
    }
    insertion_sort(&mut ranks, |a, b| a.atom < b.atom);
    for r in &mut ranks {
        r.n_ligands = 0;
    }
    for b in &mp.bonds {
        ranks[(b.atoms[0] - 1) as usize].n_ligands += 1;
        ranks[(b.atoms[1] - 1) as usize].n_ligands += 1;
    }
    // refinement on the neighbour rank sums
    let mut do_cis_trans = false;
    loop {
        for i in 0..na {
            numbering[i] = ranks[i].rank;
        }
        for r in &mut ranks {
            r.rank_sum = 0;
        }
        for b in &mp.bonds {
            let (a0, a1) = ((b.atoms[0] - 1) as usize, (b.atoms[1] - 1) as usize);
            ranks[a0].rank_sum += ranks[a1].rank;
            ranks[a1].rank_sum += ranks[a0].rank;
        }
        if do_cis_trans {
            let _ = cis_trans_perception(mp, numbering);
            for b in &mp.bonds {
                if b.color != NONE {
                    ranks[(b.atoms[0] - 1) as usize].rank_sum += b.color;
                    ranks[(b.atoms[1] - 1) as usize].rank_sum += b.color;
                }
            }
        }
        for r in &mut ranks {
            if r.n_ligands > 0 {
                r.rank_sum *= 10;
                r.rank_sum /= r.n_ligands;
            }
        }
        insertion_sort(&mut ranks, |a, b| {
            a.rank < b.rank || (a.rank == b.rank && a.rank_sum < b.rank_sum)
        });
        let mut changed = false;
        let mut j = 0;
        for i in 1..na {
            if ranks[i].rank > ranks[i - 1].rank || ranks[i].rank_sum > ranks[i - 1].rank_sum {
                j = i as i32;
            }
            changed = changed || ranks[i].rank != j;
            ranks[i].rank = j;
        }
        insertion_sort(&mut ranks, |a, b| a.atom < b.atom);
        if !changed && !do_cis_trans {
            do_cis_trans = true;
            changed = true;
        } else {
            do_cis_trans = false;
        }
        if !changed {
            break;
        }
    }
    // the coloured atoms of the smallest rank stay
    let mut min_rank = na as i32;
    for i in 0..na {
        if mp.atoms[i].color != NONE && ranks[i].rank < min_rank {
            min_rank = ranks[i].rank;
        }
    }
    let mut result = 0;
    for i in 0..na {
        if mp.atoms[i].color != NONE && ranks[i].rank == min_rank {
            result += 1;
            mp.atoms[i].value = (ranks[i].rank + 1) as f32;
        } else {
            mp.atoms[i].color = NONE;
            mp.atoms[i].value = 0.0;
        }
    }
    if result > 1 {
        for i in 0..na {
            if mp.atoms[i].color != NONE {
                add_msg(msgs, format!("atom {} in minimal rank class", i + 1));
            }
        }
    }
    result
}

/// `RechargeMolecule`: removes protons from `mp` until `desired_charge` is reached; returns success and the
/// numbers of deprotonation and refinement cycles
pub fn recharge_molecule(
    mp: &mut Molecule,
    desired_charge: i32,
    msgs: &mut Vec<String>,
) -> (bool, i32, i32) {
    let t = tables().clone();
    let (mut ndeprot, mut nrefine) = (0, 0);
    loop {
        while total_charge(mp) > desired_charge {
            set_pka_values(&t, mp, msgs);
            let (nacid, pka, gap) = mark_most_acidic_atoms(mp);
            if pka > t.acidity_limit {
                add_msg(msgs, format!("pKa_value = {}", fmt_gp(pka, 2)));
                return (false, ndeprot, nrefine);
            } else if nacid == 1 || (nacid == total_charge(mp) && gap > 8.0) {
                add_msg(msgs, format!("pKa = {}", fmt_gp(pka, 2)));
                decrement_marked_charges(mp);
            } else {
                break;
            }
            ndeprot += 1;
        }
        if total_charge(mp) > desired_charge {
            nrefine += 1;
            let mut numbering = vec![0; mp.atoms.len()];
            let mut nacid = refine_acidic_atoms(&t, mp, &mut numbering, msgs);
            if nacid == 1 {
                decrement_marked_charges(mp);
            } else if all_centers_refined(mp, &numbering, msgs) {
                for a in &mut mp.atoms {
                    if a.color != NONE && nacid != 1 {
                        nacid -= 1;
                        a.color = NONE;
                    }
                }
                decrement_marked_charges(mp);
            } else {
                add_msg(msgs, format!("{}: could not fix charges", name10(mp)));
                return (false, ndeprot, nrefine);
            }
            ndeprot += 1;
        } else {
            reset_values(mp);
            return (true, ndeprot, nrefine);
        }
    }
}
