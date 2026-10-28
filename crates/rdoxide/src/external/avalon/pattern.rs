//! Port of Avalon Toolkit's `pattern.c`: augmented atom patterns, their matching, the
//! transformations built on them and the string notation used in the STRUCHK rule files.
//! The message list of C is a `Vec<String>`; buffers are `String`s.

use super::perceive::{
    AROMATIC, CHAIN, DOUBLE, Neighbourhood, RING, SINGLE, TRIPLE, ring_state, setup_neighbourhood,
};
use super::reaccs::{ANY_CHARGE, Molecule, NONE};
use super::stereo::add_msg;
use super::utilities::{atom_symbol_match, reset_colors, sort_neighbourhood};

pub const ANY_RADICAL: i32 = 8;
pub const MAXNEIGHBOURS: usize = 20;
const SINGLET: i32 = 1;
const DOUBLET: i32 = 2;
const TRIPLET: i32 = 3;
const SINGLE_DOUBLE: i32 = 5;
const SINGLE_AROMATIC: i32 = 6;
const DOUBLE_AROMATIC: i32 = 7;
const ANY_BOND: i32 = 8;
const NOT_CENTER: i32 = -1;
const UNCHANGED: i32 = 2;
/// colour that marks molecules produced by a transformation
pub const TRANSFORMED: i32 = 1;

/// `bond_to_string`
const BOND_TO_STRING: &[(i32, &str)] = &[
    (NONE, "?"),
    (SINGLE, "-"),
    (DOUBLE, "="),
    (TRIPLE, "#"),
    (AROMATIC, "~"),
    (SINGLE_DOUBLE, "-="),
    (SINGLE_AROMATIC, "-~"),
    (DOUBLE_AROMATIC, "=~"),
    (ANY_BOND, "*"),
];

/// `charge_to_string`
const CHARGE_TO_STRING: &[(i32, &str)] = &[
    (NONE, ""),
    (7, "+7"),
    (6, "+6"),
    (5, "+5"),
    (4, "+4"),
    (3, "+3"),
    (2, "+2"),
    (1, "+1"),
    (-1, "-1"),
    (-2, "-2"),
    (-3, "-3"),
    (-4, "-4"),
    (-5, "-5"),
    (-6, "-6"),
    (-7, "-7"),
    (ANY_CHARGE, "+-"),
];

/// `radical_to_string`
const RADICAL_TO_STRING: &[(i32, &str)] =
    &[(NONE, ""), (SINGLET, "|"), (DOUBLET, "."), (TRIPLET, ":")];

/// `bondtext` of `ConvertMoleculeToString`
const BONDTEXT: [&str; 9] = ["?", "-", "=", "^", "%", ":", "?", "#", "~"];

/// `H_counttext`
const H_COUNTTEXT: [&str; 8] = ["", "H0", "H", "H", "H", "H", "H", "H"];

/// `IdToString` with `SAFEGUARD`
fn id_to_string(table: &[(i32, &'static str)], id: i32) -> &'static str {
    table.iter().find(|e| e.0 == id).map_or("(NULL)", |e| e.1)
}

/// `ligand_t`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Ligand {
    pub atom_symbol: String,
    pub charge: i32,
    pub radical: i32,
    pub bond_type: i32,
    /// substitution count, 0 means "don't care"
    pub s_count: i32,
}

/// `augmented_atom_t`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AugmentedAtom {
    pub charge: i32,
    pub radical: i32,
    pub atom_symbol: String,
    pub short_name: String,
    /// 1 means ring, 2 means chain, 0 "don't care"
    pub topography: i32,
    pub ligands: Vec<Ligand>,
}

/// `RecMatch`
fn rec_match(
    mp: &Molecule,
    matched: &mut [usize],
    level: usize,
    aap: &AugmentedAtom,
    nbp: &[Neighbourhood],
) -> bool {
    if level == aap.ligands.len() {
        return true;
    }
    let lig = &aap.ligands[level];
    let nbph = &nbp[matched[0]];
    for i in 0..nbph.atoms.len() {
        let a = &mp.atoms[nbph.atoms[i]];
        if (a.charge == lig.charge || lig.charge == ANY_CHARGE)
            && (a.radical == lig.radical || lig.radical == ANY_RADICAL)
            && mp.bonds[nbph.bonds[i]].bond_type == lig.bond_type
            && (lig.s_count == NONE || nbp[nbph.atoms[i]].atoms.len() as i32 == lig.s_count)
            && atom_symbol_match(&a.atom_symbol, &lig.atom_symbol)
        {
            matched[level + 1] = nbph.atoms[i];
            let is_new = !matched[..=level].contains(&nbph.atoms[i]);
            if is_new && rec_match(mp, matched, level + 1, aap, nbp) {
                return true;
            }
        }
    }
    false
}

/// `AAMatch`: the atom mapping (center first, then the ligands) if atom `i` matches `aap`
pub fn aa_match(
    mp: &Molecule,
    i: usize,
    aap: &AugmentedAtom,
    atom_ring_status: Option<&[i32]>,
    nbp: &[Neighbourhood],
) -> Option<Vec<usize>> {
    let a = &mp.atoms[i];
    if nbp[i].atoms.len() == aap.ligands.len()
        && (aap.charge == ANY_CHARGE || a.charge == aap.charge)
        && (aap.radical == ANY_RADICAL || a.radical == aap.radical)
        && atom_symbol_match(&a.atom_symbol, &aap.atom_symbol)
    {
        if let Some(st) = atom_ring_status {
            if aap.topography == RING && st[i] == 0 {
                return None;
            }
            if aap.topography == CHAIN && st[i] != 0 {
                return None;
            }
        }
        let mut matched = vec![0usize; MAXNEIGHBOURS + 1];
        matched[0] = i;
        if rec_match(mp, &mut matched, 0, aap, nbp) {
            return Some(matched);
        }
    }
    None
}

/// `TransformAA`
pub fn transform_aa(
    mp: &mut Molecule,
    i: usize,
    matched: &[usize],
    aapin: &AugmentedAtom,
    aapout: &AugmentedAtom,
    nbp: &[Neighbourhood],
) {
    if aapin.charge != aapout.charge {
        mp.atoms[i].charge = aapout.charge;
    }
    if aapin.radical != aapout.radical {
        mp.atoms[i].radical = aapout.radical;
    }
    if aapin.atom_symbol != aapout.atom_symbol {
        mp.atoms[i].atom_symbol = aapout.atom_symbol.clone();
    }
    for j in 0..nbp[i].atoms.len() {
        for k in 0..aapout.ligands.len() {
            if nbp[i].atoms[j] == matched[k + 1] {
                let (li, lo) = (&aapin.ligands[k], &aapout.ligands[k]);
                let at = nbp[i].atoms[j];
                if li.charge != lo.charge {
                    mp.atoms[at].charge = lo.charge;
                }
                if li.radical != lo.radical {
                    mp.atoms[at].radical = lo.radical;
                }
                if li.atom_symbol != lo.atom_symbol {
                    mp.atoms[at].atom_symbol = lo.atom_symbol.clone();
                }
                if li.bond_type != lo.bond_type {
                    mp.bonds[nbp[i].bonds[j]].bond_type = lo.bond_type;
                }
            }
        }
    }
}

/// `AAFix`: replaces all occurrences of `aap[0]` with `aap[1]`
pub fn aa_fix(
    mp: &mut Molecule,
    aap: &[AugmentedAtom; 2],
    nbp: &[Neighbourhood],
    msgs: &mut Vec<String>,
) -> bool {
    let (atom_status, _) = ring_state(mp);
    let mut matched = false;
    for i in 0..mp.atoms.len() {
        if let Some(m) = aa_match(mp, i, &aap[0], Some(&atom_status), nbp) {
            matched = true;
            transform_aa(mp, i, &m, &aap[0], &aap[1], nbp);
            let mut s = format!(
                "{} {} : '{}' -> '{}'",
                mp.name,
                i + 1,
                aap[0].short_name,
                aap[1].short_name
            );
            s.truncate(1023);
            add_msg(msgs, s);
        }
    }
    matched
}

/// `TransformAugmentedAtoms`
pub fn transform_augmented_atoms(
    mp: &mut Molecule,
    table: &[[AugmentedAtom; 2]],
    msgs: &mut Vec<String>,
) -> Result<bool, &'static str> {
    if mp.atoms.is_empty() {
        return Ok(false);
    }
    let nbp = setup_neighbourhood(mp).ok_or("too many neighbours")?;
    let mut result = false;
    for pair in table {
        if aa_fix(mp, pair, &nbp, msgs) {
            result = true;
        }
    }
    mp.bonds.retain(|b| b.bond_type != NONE);
    Ok(result)
}

/// `TraceBranch`
fn trace_branch(
    root: usize,
    nmapped: &mut i32,
    branch_bond: usize,
    branch_atom: usize,
    mp: &mut Molecule,
    nbarr: &[Neighbourhood],
) -> String {
    let a = &mp.atoms[branch_atom];
    let mut buffer = String::new();
    buffer += BONDTEXT
        .get(mp.bonds[branch_bond].bond_type as usize)
        .copied()
        .unwrap_or("?");
    buffer += &a.atom_symbol;
    buffer += H_COUNTTEXT
        .get(a.query_h_count as usize)
        .copied()
        .unwrap_or("");
    buffer += id_to_string(CHARGE_TO_STRING, a.charge);
    buffer += id_to_string(RADICAL_TO_STRING, a.radical);
    *nmapped += 1;
    mp.atoms[branch_atom].color = *nmapped;

    let nbp = &nbarr[branch_atom];
    let mut local = String::new();
    for i in 0..nbp.atoms.len() {
        if nbp.atoms[i] != root && mp.atoms[nbp.atoms[i]].color < mp.atoms[branch_atom].color {
            if !local.is_empty() {
                buffer += "(";
                buffer += &local;
                buffer += ")";
                local.clear();
            }
            let c = mp.atoms[nbp.atoms[i]].color;
            if c != NONE {
                local = format!(
                    "{}@{}",
                    BONDTEXT
                        .get(mp.bonds[nbp.bonds[i]].bond_type as usize)
                        .copied()
                        .unwrap_or("?"),
                    c
                );
            } else {
                local = trace_branch(branch_atom, nmapped, nbp.bonds[i], nbp.atoms[i], mp, nbarr);
            }
        }
    }
    buffer += &local;
    buffer
}

/// `ConvertMoleculeToString`; also leaves the traversal order in the atom colours
pub fn convert_molecule_to_string(mp: &mut Molecule) -> Result<String, &'static str> {
    let nbarr = setup_neighbourhood(mp).ok_or("too many neighbours")?;
    reset_colors(mp);
    let mut string = String::new();
    let mut nmapped = 0;
    loop {
        let mut minlig = MAXNEIGHBOURS;
        let mut root = None;
        for i in 0..mp.atoms.len() {
            if mp.atoms[i].color == NONE && nbarr[i].atoms.len() < minlig {
                minlig = nbarr[i].atoms.len();
                root = Some(i);
            }
        }
        let Some(root) = root else { break };
        if !string.is_empty() {
            string += " ";
        }
        nmapped += 1;
        mp.atoms[root].color = nmapped;
        let a = &mp.atoms[root];
        string += &a.atom_symbol;
        string += H_COUNTTEXT
            .get(a.query_h_count as usize)
            .copied()
            .unwrap_or("");
        string += id_to_string(CHARGE_TO_STRING, a.charge);
        string += id_to_string(RADICAL_TO_STRING, a.radical);
        let nbp = &nbarr[root];
        let mut local = String::new();
        for i in 0..nbp.atoms.len() {
            if mp.atoms[nbp.atoms[i]].color == NONE {
                if !local.is_empty() {
                    string += "(";
                    string += &local;
                    string += ")";
                    local.clear();
                }
                local = trace_branch(root, &mut nmapped, nbp.bonds[i], nbp.atoms[i], mp, &nbarr);
            }
        }
        string += &local;
    }
    string += "/";
    Ok(string)
}

/// `AAToString` (and the text `AAPrint` writes): the sorted neighbourhood of atom `i`
pub fn aa_to_string(mp: &Molecule, i: usize, np: &mut Neighbourhood) -> String {
    sort_neighbourhood(np, mp);
    let a = &mp.atoms[i];
    let mut s = format!(
        "{}{}{}",
        a.atom_symbol,
        id_to_string(CHARGE_TO_STRING, a.charge),
        id_to_string(RADICAL_TO_STRING, a.radical)
    );
    for j in 0..np.atoms.len() {
        let l = &mp.atoms[np.atoms[j]];
        s += &format!(
            "({}{}{}{})",
            id_to_string(BOND_TO_STRING, mp.bonds[np.bonds[j]].bond_type),
            l.atom_symbol,
            id_to_string(CHARGE_TO_STRING, l.charge),
            id_to_string(RADICAL_TO_STRING, l.radical)
        );
    }
    s
}

fn is_alnum_sym(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b',' || c == b'#'
}

/// reads a charge or radical code at the start of `s` from `table`; `None` if no entry matches
fn read_code(table: &[(i32, &str)], s: &[u8]) -> Option<(i32, usize)> {
    table
        .iter()
        .find(|e| !e.1.is_empty() && s.starts_with(e.1.as_bytes()))
        .map(|e| (e.0, e.1.len()))
}

/// `StringToAugmentedAtom`
pub fn string_to_augmented_atom(text: &str) -> Option<AugmentedAtom> {
    let b = text.as_bytes();
    let mut p = 0;
    let at = |i: usize| b.get(i).copied().unwrap_or(0);
    let mut aap = AugmentedAtom::default();
    if at(p) == b'@' {
        aap.topography = RING;
        p += 1;
    } else if at(p) == b'!' && at(p + 1) == b'@' {
        aap.topography = CHAIN;
        p += 2;
    }
    if !at(p).is_ascii_alphabetic() {
        return None;
    }
    let st = p;
    while is_alnum_sym(at(p)) {
        p += 1;
    }
    let atsym = text[st..p].to_string();
    // charge and radical codes
    let code = |p: &mut usize, charge: &mut i32, radical: &mut i32| -> Option<()> {
        if at(*p) == b'+' || at(*p) == b'-' {
            let (v, n) = read_code(CHARGE_TO_STRING, &b[*p..])?;
            *charge = v;
            *p += n;
        } else {
            *charge = NONE;
        }
        if matches!(at(*p), b'.' | b'|' | b':') && !at(*p + 1).is_ascii_digit() {
            let (v, n) = read_code(RADICAL_TO_STRING, &b[*p..])?;
            *radical = v;
            *p += n;
        } else {
            *radical = NONE;
        }
        Some(())
    };
    let (mut c, mut r) = (0, 0);
    code(&mut p, &mut c, &mut r)?;
    aap.charge = c;
    aap.radical = r;
    while at(p) == b'(' {
        p += 1;
        if aap.ligands.len() >= MAXNEIGHBOURS {
            return None;
        }
        let mut lig = Ligand::default();
        let e = BOND_TO_STRING
            .iter()
            .find(|e| b[p..].starts_with(e.1.as_bytes()))?;
        lig.bond_type = e.0;
        p += e.1.len();
        if !at(p).is_ascii_alphabetic() {
            return None;
        }
        let st = p;
        while is_alnum_sym(at(p)) {
            p += 1;
        }
        lig.atom_symbol = text[st..p].to_string();
        code(&mut p, &mut lig.charge, &mut lig.radical)?;
        if at(p) == b':' && at(p + 1).is_ascii_digit() {
            lig.s_count = i32::from(at(p + 1) - b'0');
            p += 2;
        } else {
            lig.s_count = 0;
        }
        if at(p) != b')' {
            return None;
        }
        p += 1;
        aap.ligands.push(lig);
    }
    if p != b.len() {
        return None;
    }
    aap.atom_symbol = atsym;
    aap.short_name = text.to_string();
    Some(aap)
}

/// emulation of the `fscanf`/`SearchChar` reading of the rule files
struct Cursor<'a> {
    s: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    /// `SearchChar`: consumes everything up to and including `tag`
    fn search(&mut self, tag: u8) {
        while self.pos < self.s.len() {
            let c = self.s[self.pos];
            self.pos += 1;
            if c == tag {
                break;
            }
        }
    }

    /// `%[^"]` with a maximum width
    fn until_quote(&mut self, max: usize) -> String {
        let st = self.pos;
        while self.pos < self.s.len() && self.s[self.pos] != b'"' && self.pos - st < max {
            self.pos += 1;
        }
        String::from_utf8_lossy(&self.s[st..self.pos]).into_owned()
    }

    /// `fscanf("%d")` followed by `fgets(buffer, 80, fp)`: the count and the version found
    /// between the first two underscores of the rest of the line
    fn count_and_version(&mut self) -> (i32, String) {
        while self.pos < self.s.len() && self.s[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
        let st = self.pos;
        while self.pos < self.s.len() && self.s[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        let n = std::str::from_utf8(&self.s[st..self.pos])
            .ok()
            .and_then(|t| t.parse().ok())
            .unwrap_or(0);
        let ls = self.pos;
        while self.pos < self.s.len() && self.pos - ls < 79 {
            let c = self.s[self.pos];
            self.pos += 1;
            if c == b'\n' {
                break;
            }
        }
        let line = String::from_utf8_lossy(&self.s[ls..self.pos]).into_owned();
        let mut version = String::new();
        if let Some(a) = line.find('_') {
            if let Some(l) = line[a + 1..].find('_') {
                version = line[a + 1..a + 1 + l].to_string();
            }
        }
        (n, version)
    }
}

/// `ReadAugmentedAtoms`: the atoms and the version string, `None` if a pattern is malformed
pub fn read_augmented_atoms(text: &str) -> Option<(Vec<AugmentedAtom>, String)> {
    let mut cur = Cursor {
        s: text.as_bytes(),
        pos: 0,
    };
    let (n, version) = cur.count_and_version();
    let mut atoms = Vec::new();
    let mut success = true;
    for _ in 0..n {
        cur.search(b'"');
        let s = cur.until_quote(499);
        cur.search(b'\n');
        match string_to_augmented_atom(&s) {
            Some(a) => atoms.push(a),
            None => success = false,
        }
    }
    success.then_some((atoms, version))
}

/// `ReadAAPairs`: the transformation pairs and the version string
pub fn read_aa_pairs(text: &str) -> Option<(Vec<[AugmentedAtom; 2]>, String)> {
    let mut cur = Cursor {
        s: text.as_bytes(),
        pos: 0,
    };
    let (n, version) = cur.count_and_version();
    let mut pairs = Vec::new();
    let mut success = true;
    for _ in 0..n {
        cur.search(b'"');
        let s0 = cur.until_quote(999);
        cur.search(b'"');
        cur.search(b'"');
        let s1 = cur.until_quote(999);
        cur.search(b'\n');
        match (string_to_augmented_atom(&s0), string_to_augmented_atom(&s1)) {
            (Some(a), Some(b)) => pairs.push([a, b]),
            _ => success = false,
        }
    }
    success.then_some((pairs, version))
}

/// `IsUseAtom`: atom number `n` (1 based) has a reaction center bond
fn is_use_atom(n: i32, mp: &Molecule) -> bool {
    mp.bonds.iter().any(|b| {
        b.reaction_mark != NONE
            && b.reaction_mark != NOT_CENTER
            && b.reaction_mark != UNCHANGED
            && (b.atoms[0] == n || b.atoms[1] == n)
    })
}

/// `AlternativeAATransformation`: all molecules that result from the alternative
/// transformations in `tfm_table`, prepended to `old_list` (list order, head first)
pub fn alternative_aa_transformation(
    mut mp: Molecule,
    mut next_atom: usize,
    tfm_table: &[[AugmentedAtom; 2]],
    old_list: Vec<Molecule>,
    use_old: bool,
) -> Result<Vec<Molecule>, &'static str> {
    let (atom_status, _) = ring_state(&mp);
    let nbp = setup_neighbourhood(&mp).ok_or("too many neighbours")?;
    let matches_at = |mp: &Molecule, i: usize, j: usize| {
        if is_use_atom(i as i32 + 1, mp) {
            aa_match(mp, i, &tfm_table[j][0], Some(&atom_status), &nbp)
        } else {
            None
        }
    };
    loop {
        let n = mp.atoms.len();
        let mut i = next_atom;
        let (mut nmatch, mut jmatch) = (0, 0);
        while i < n {
            nmatch = 0;
            for j in 0..tfm_table.len() {
                if matches_at(&mp, i, j).is_some() {
                    jmatch = j;
                    nmatch += 1;
                }
            }
            if nmatch > 0 {
                break;
            }
            i += 1;
        }
        let branch = |mp: &Molecule, result: &mut Vec<Molecule>| -> Result<(), &'static str> {
            for j in 0..tfm_table.len() {
                if let Some(m) = matches_at(mp, i, j) {
                    let mut mph = mp.clone();
                    transform_aa(&mut mph, i, &m, &tfm_table[j][0], &tfm_table[j][1], &nbp);
                    mph.color = TRANSFORMED;
                    let old = std::mem::take(result);
                    *result = alternative_aa_transformation(mph, i + 1, tfm_table, old, use_old)?;
                }
            }
            Ok(())
        };
        if use_old {
            if i == n {
                let mut r = vec![mp];
                r.extend(old_list);
                return Ok(r);
            }
            let mut result = old_list;
            branch(&mp, &mut result)?;
            return alternative_aa_transformation(mp, i + 1, tfm_table, result, use_old);
        } else if i == n {
            if next_atom == 0 {
                return Ok(old_list);
            }
            let mut r = vec![mp];
            r.extend(old_list);
            return Ok(r);
        } else if nmatch == 1 {
            let m = matches_at(&mp, i, jmatch).expect("match");
            transform_aa(
                &mut mp,
                i,
                &m,
                &tfm_table[jmatch][0],
                &tfm_table[jmatch][1],
                &nbp,
            );
            mp.color = TRANSFORMED;
            next_atom = i + 1;
        } else {
            let mut result = old_list;
            branch(&mp, &mut result)?;
            return Ok(result);
        }
    }
}

/// `AAPrint` writes this text and returns its length
pub fn aa_print(mp: &Molecule, i: usize, np: &mut Neighbourhood) -> String {
    aa_to_string(mp, i, np)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::aacheck::check_atoms;
    use crate::external::avalon::reaccsio::{mol_str_to_mol, mol_to_mol_str};

    const CHK: &str = include_str!("data/checkfgs.chk");
    const TRN: &str = include_str!("data/checkfgs.trn");

    fn flush(tag: &str, msgs: &mut Vec<String>) -> String {
        let mut s = format!("{tag}:\n");
        for m in msgs.drain(..) {
            s += &m;
            s += "\n";
        }
        s
    }

    #[test]
    fn augmented_atom_strings_match_c() {
        let text = include_str!("../../../../../tests/data/avalon_aaparse.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let got = match string_to_augmented_atom(input) {
                None => "FAIL".to_string(),
                Some(a) => {
                    let mut s = format!(
                        "{} {} {} {} {} {}",
                        a.ligands.len(),
                        a.charge,
                        a.radical,
                        a.atom_symbol,
                        a.short_name,
                        a.topography
                    );
                    for l in &a.ligands {
                        s += &format!(
                            " [{} {} {} {} {}]",
                            l.atom_symbol, l.charge, l.radical, l.bond_type, l.s_count
                        );
                    }
                    s
                }
            };
            if got != want {
                bad.push(format!("{input:?}\nwant {want}\ngot  {got}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    #[test]
    fn pattern_functions_match_c() {
        let (good, chk_version) = read_augmented_atoms(CHK).expect("rules");
        let (pairs, trn_version) = read_aa_pairs(TRN).expect("pairs");
        let text = include_str!("../../../../../tests/data/avalon_pattern.jsonl");
        let mut bad = Vec::new();
        let mut n = 0;
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let mp = mol_str_to_mol(input).expect("molecule");
            let mut out = format!(
                "read {} {} [{}] [{}]\n",
                good.len(),
                pairs.len(),
                chk_version,
                trn_version
            );
            let mut c = mp.clone();
            out += &format!("str={}\n", convert_molecule_to_string(&mut c).unwrap());
            let mut msgs = Vec::new();
            out += &format!(
                "check={}\n",
                i32::from(check_atoms(&mp, &good, &mut msgs).unwrap())
            );
            out += &flush("chkmsgs", &mut msgs);
            let mut c = mp.clone();
            out += &format!(
                "trans={}\n",
                i32::from(transform_augmented_atoms(&mut c, &pairs, &mut msgs).unwrap())
            );
            out += &flush("trnmsgs", &mut msgs);
            out += &mol_to_mol_str(&c);
            for use_old in [false, true] {
                let mut c = mp.clone();
                for (i, b) in c.bonds.iter_mut().enumerate() {
                    b.reaction_mark = if i % 3 == 2 { 0 } else { 1 };
                }
                let list =
                    alternative_aa_transformation(c, 0, &pairs, Vec::new(), use_old).unwrap();
                out += &format!("alt{}={}\n", i32::from(use_old), list.len());
                for h in &list {
                    out += &format!("color={}\n{}", h.color, mol_to_mol_str(h));
                }
            }
            n += 1;
            if out != want {
                bad.push(format!("{input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(
            bad.is_empty(),
            "{} of {n} differ; first:\n{}",
            bad.len(),
            bad[0]
        );
    }
}
