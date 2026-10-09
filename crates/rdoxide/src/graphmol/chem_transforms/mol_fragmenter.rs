//! Fragmentation of molecules on bonds (`fragmentOnBonds`,
//! `fragmentOnSomeBonds`, `fragmentOnBRICSBonds`) and the parsers for the
//! atom-type and bond-type definition tables they use (BRICS and friends).
//!
//! Ported from RDKit 2026.09.1 `GraphMol/ChemTransforms/MolFragmenter.cpp`.
//!
//! Copyright (C) 2013-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::sync::{Arc, OnceLock};

use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::threads::par_map;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_error, rd_log_warning};

/// `FragmenterBondType`: a bond definition, i.e. a two-atom SMARTS query plus
/// the labels of the dummies that replace the bond's atoms.
#[derive(Clone, Debug)]
pub struct FragmenterBondType {
    pub atom1_label: u32,
    pub atom2_label: u32,
    pub atom1_type: u32,
    pub atom2_type: u32,
    pub bond_type: BondType,
    pub query: Option<Arc<ROMol>>,
}

/// Property used to remember a neighboring double bond's stereo while its
/// stereo atoms are being changed.
const MOLFRAG_SAVE_STEREO: &str = "_molfragSaveStereo";
/// Property remembering the new bond order around an atom whose neighbor
/// moved (github #1734).
const NEW_BOND_ORDER: &str = "_newBondOrder";

fn smarts_to_mol(smarts: &str) -> Option<ROMol> {
    mol_from_smarts(smarts, &SmartsParserParams::default())
        .ok()
        .flatten()
}

/// `boost::split(tokens, line, is_any_of(" \t"), token_compress_on)`:
/// runs of separators count as one, but a leading or trailing separator still
/// delimits an empty token.
fn split_tokens(line: &str) -> Vec<&str> {
    let mut tokens = Vec::new();
    let mut start = 0;
    let mut previous_was_sep = false;
    for (i, c) in line.char_indices() {
        let is_sep = c == ' ' || c == '\t';
        if is_sep && !previous_was_sep {
            tokens.push(&line[start..i]);
        }
        if is_sep {
            start = i + c.len_utf8();
        }
        previous_was_sep = is_sep;
    }
    tokens.push(&line[start..]);
    tokens
}

/// The meaningful lines (with their 1-based numbers) of a definition table:
/// everything but empty lines and lines starting with `comment`.
fn definition_lines<'a>(text: &'a str, comment: &'a str) -> impl Iterator<Item = (usize, &'a str)> {
    text.split('\n')
        .map(|l| l.strip_suffix('\r').unwrap_or(l))
        .enumerate()
        .filter(move |(_, l)| !l.is_empty() && !l.starts_with(comment))
        .map(|(i, l)| (i + 1, l))
}

fn lexical_cast_uint(token: &str) -> Result<u32> {
    token.parse().map_err(|_| {
        RDKitError::Runtime(
            "bad lexical cast: source type value could not be interpreted as target".into(),
        )
    })
}

/// `constructFragmenterAtomTypes`: reads `index SMARTS` lines. With `validate`
/// (or when `environs` is given) the SMARTS must parse; the parsed
/// environments are stored in `environs`.
pub fn construct_fragmenter_atom_types(
    text: &str,
    comment: &str,
    validate: bool,
    mut environs: Option<&mut BTreeMap<u32, ROMol>>,
) -> Result<BTreeMap<u32, String>> {
    let mut defs = BTreeMap::new();
    for (line, content) in definition_lines(text, comment) {
        let tokens = split_tokens(content);
        if tokens.len() < 2 {
            rd_log_warning!("line {line} is too short\n");
            continue;
        }
        let idx = lexical_cast_uint(tokens[0])?;
        if defs.contains_key(&idx) {
            rd_log_warning!(
                "definition #{idx} encountered more than once. Using the first occurrence.\n"
            );
            continue;
        }
        if validate || environs.is_some() {
            let Some(p) = smarts_to_mol(tokens[1]) else {
                rd_log_warning!(
                    "cannot convert SMARTS {} to molecule at line {line}\n",
                    tokens[1]
                );
                continue;
            };
            if let Some(environs) = environs.as_deref_mut() {
                environs.insert(idx, p);
            }
        }
        defs.insert(idx, tokens[1].to_string());
    }
    Ok(defs)
}

/// `constructBRICSAtomTypes`: the BRICS environments (L1-L16).
pub fn construct_brics_atom_types(
    environs: Option<&mut BTreeMap<u32, ROMol>>,
) -> Result<BTreeMap<u32, String>> {
    // The L2 definitions ("N.pl3" in the original paper) have been
    // incorporated into a (almost) general purpose amine definition in L5
    // ("N.sp3" in the paper), for consistency of the fragmentations.
    const BRICS_DEFS: &str = "1 [C;D3]([#0,#6,#7,#8])(=O)
3 [O;D2]-;!@[#0,#6,#1]
5 [N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]
9 [n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]
10 [N;R;$(N(@C(=O))@[C,N,O,S])]
11 [S;D2](-;!@[#0,#6])
12 [S;D4]([#6,#0])(=O)(=O)
6 [C;D3;!R](=O)-;!@[#0,#6,#7,#8]
13 [C;$(C(-;@[C,N,O,S])-;@[N,O,S])]
14 [c;$(c(:[c,n,o,s]):[n,o,s])]
15 [C;$(C(-;@C)-;@C)]
4 [C;!D1;!$(C=*)]-;!@[#6]
7 [C;D2,D3]-[#6]
8 [C;!R;!D1;!$(C!-*)]
16 [c;$(c(:c):c)]";
    construct_fragmenter_atom_types(BRICS_DEFS, "//", true, environs)
}

/// `constructFragmenterBondTypes`: reads `atomType1 atomType2 bondSMARTS`
/// lines into bond definitions.
pub fn construct_fragmenter_bond_types(
    text: &str,
    atom_types: &BTreeMap<u32, String>,
    comment: &str,
    validate: bool,
    label_by_connector: bool,
) -> Result<Vec<FragmenterBondType>> {
    let mut defs = Vec::new();
    for (line, content) in definition_lines(text, comment) {
        let tokens = split_tokens(content);
        if tokens.len() < 3 {
            rd_log_warning!("line {line} is too short\n");
            continue;
        }
        let atom_type = |token: &str| -> Result<Option<(u32, &String)>> {
            let idx = lexical_cast_uint(token)?;
            let sma = atom_types.get(&idx);
            if sma.is_none() {
                rd_log_warning!("atom type #{idx} not recognized.\n");
            }
            Ok(sma.map(|s| (idx, s)))
        };
        let Some((idx1, sma1)) = atom_type(tokens[0])? else {
            continue;
        };
        let Some((idx2, sma2)) = atom_type(tokens[1])? else {
            continue;
        };
        let smarts = format!("[$({sma1})]{}[$({sma2})]", tokens[2]);
        let p = smarts_to_mol(&smarts);
        if validate && p.is_none() {
            rd_log_warning!("cannot convert SMARTS {smarts} to molecule at line {line}\n");
            continue;
        }
        let (atom1_label, atom2_label) = if label_by_connector {
            (idx1, idx2)
        } else {
            (idx2, idx1)
        };
        let (bond_type, query) = match p {
            Some(p) => {
                // for the purposes of replacing the bond, use just the first
                // character to set the bond type (if we recognize it)
                let bond_type = match tokens[2].chars().next() {
                    Some('-') => BondType::Single,
                    Some('=') => BondType::Double,
                    Some('#') => BondType::Triple,
                    Some(':') => BondType::Aromatic,
                    _ => p.get_bond(0)?.bond_type(),
                };
                (bond_type, Some(Arc::new(p)))
            }
            None => (BondType::Unspecified, None),
        };
        defs.push(FragmenterBondType {
            atom1_label,
            atom2_label,
            atom1_type: idx1,
            atom2_type: idx2,
            bond_type,
            query,
        });
    }
    Ok(defs)
}

/// `constructBRICSBondTypes`: the BRICS bond definitions.
pub fn construct_brics_bond_types() -> Result<Vec<FragmenterBondType>> {
    const BRICS_DEFS: &str = "// L1
1 3 -;!@
1 5 -;!@
1 10 -;!@
// L3 
3 4 -;!@
3 13 -;!@
3 14 -;!@
3 15 -;!@
3 16 -;!@
// L4
4 5 -;!@
4 11 -;!@
// L5
5 12 -;!@
5 14 -;!@
5 16 -;!@
5 13 -;!@
5 15 -;!@
// L6
6 13 -;!@
6 14 -;!@
6 15 -;!@
6 16 -;!@
// L7
7 7 =;!@
// L8
8 9 -;!@
8 10 -;!@
8 13 -;!@
8 14 -;!@
8 15 -;!@
8 16 -;!@
// L9
9 13 -;!@ // not in original paper
9 14 -;!@ // not in original paper
9 15 -;!@
9 16 -;!@
// L10
10 13 -;!@
10 14 -;!@
10 15 -;!@
10 16 -;!@
// L11
11 13 -;!@
11 14 -;!@
11 15 -;!@
11 16 -;!@
// L12
// none left
// L13
13 14 -;!@
13 15 -;!@
13 16 -;!@
// L14
14 14 -;!@ // not in original paper
14 15 -;!@
14 16 -;!@
// L15
15 16 -;!@
// L16
16 16 -;!@ // not in original paper";
    let atom_types = construct_brics_atom_types(None)?;
    construct_fragmenter_bond_types(BRICS_DEFS, &atom_types, "//", true, false)
}

/// The next integer with the same number of set bits
/// (<http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation>).
fn next_bit_combo(v: u64) -> u64 {
    let t = (v | v.wrapping_sub(1)).wrapping_add(1);
    t | ((((t & t.wrapping_neg()) / (v & v.wrapping_neg())) >> 1).wrapping_sub(1))
}

/// Dummy labels and bond types are given per cut bond.
fn check_per_bond_args(
    bond_indices: &[u32],
    dummy_labels: Option<&[(u32, u32)]>,
    bond_types: Option<&[BondType]>,
) -> Result<()> {
    if dummy_labels.is_some_and(|d| d.len() < bond_indices.len()) {
        return Err(RDKitError::precondition("bad dummyLabel vector"));
    }
    if bond_types.is_some_and(|b| b.len() != bond_indices.len()) {
        return Err(RDKitError::precondition("bad bondType vector"));
    }
    let unique: HashSet<_> = bond_indices.iter().collect();
    if unique.len() < bond_indices.len() {
        return Err(RDKitError::Value("bondIndices contains duplicates".into()));
    }
    Ok(())
}

/// `fragmentOnSomeBonds`: the fragmentations obtained by cutting every
/// combination of exactly `max_to_cut` of `bond_indices` (in the order of
/// the bit-permutation enumeration of upstream). The molecules are produced
/// in parallel. Also returns the number of cuts per atom of each result when
/// `want_cuts_per_atom` is set.
pub fn fragment_on_some_bonds(
    mol: &ROMol,
    bond_indices: &[u32],
    max_to_cut: u32,
    add_dummies: bool,
    dummy_labels: Option<&[(u32, u32)]>,
    bond_types: Option<&[BondType]>,
    want_cuts_per_atom: bool,
) -> Result<(Vec<RWMol>, Vec<Vec<u32>>)> {
    check_per_bond_args(bond_indices, dummy_labels, bond_types)?;
    if bond_indices.len() > 63 {
        return Err(RDKitError::Value(
            "currently can only fragment on up to 63 bonds".into(),
        ));
    }
    if max_to_cut == 0
        || max_to_cut as usize > bond_indices.len()
        || mol.num_atoms() == 0
        || bond_indices.is_empty()
    {
        return Ok((Vec::new(), Vec::new()));
    }

    // every combination of bonds to cut, with the matching labels and types
    struct Cut {
        bonds: Vec<u32>,
        labels: Option<Vec<(u32, u32)>>,
        types: Option<Vec<BondType>>,
    }
    let mut cuts = Vec::new();
    let stop = 1u64 << bond_indices.len();
    let mut state = (1u64 << max_to_cut) - 1;
    while state < stop {
        let chosen: Vec<usize> = (0..bond_indices.len())
            .filter(|i| state & (1u64 << i) != 0)
            .take(max_to_cut as usize)
            .collect();
        cuts.push(Cut {
            bonds: chosen.iter().map(|&i| bond_indices[i]).collect(),
            labels: dummy_labels.map(|d| chosen.iter().map(|&i| d[i]).collect()),
            types: bond_types.map(|t| chosen.iter().map(|&i| t[i]).collect()),
        });
        state = next_bit_combo(state);
    }

    let results = par_map(&cuts, 0, |cut| {
        let mut n_cuts = vec![0u32; mol.num_atoms() as usize];
        let res = fragment_on_bonds(
            mol,
            &cut.bonds,
            add_dummies,
            cut.labels.as_deref(),
            cut.types.as_deref(),
            want_cuts_per_atom.then_some(n_cuts.as_mut_slice()),
        )?;
        Ok::<_, RDKitError>((res, n_cuts))
    });
    let mut mols = Vec::with_capacity(results.len());
    let mut cuts_per_atom = Vec::new();
    for r in results {
        let (m, n_cuts) = r?;
        mols.push(m);
        if want_cuts_per_atom {
            cuts_per_atom.push(n_cuts);
        }
    }
    Ok((mols, cuts_per_atom))
}

/// `checkChiralityPostMove`: sets the chirality of `new_atom` (which keeps the
/// place of `old_atom` after `bond` was moved to a new neighbor) so that the
/// permutation of its bonds is accounted for.
fn check_chirality_post_move(
    mol: &ROMol,
    old_atom: u32,
    res: &mut RWMol,
    new_atom: u32,
    bond: u32,
) -> Result<()> {
    // We may be called more than once for an atom, so the new atom order is
    // kept in a property (github #1734).
    let mut new_order: Vec<i32> = match res.atom(new_atom).get_prop(NEW_BOND_ORDER) {
        Some(RDValue::VecInt(incoming)) => incoming
            .iter()
            .copied()
            .filter(|&b| b != bond as i32)
            .collect(),
        _ => mol
            .atom_bond_indices(old_atom)
            .filter(|&b| b != bond)
            .map(|b| b as i32)
            .collect(),
    };
    new_order.push(bond as i32);
    res.atom_mut(new_atom)
        .set_prop(NEW_BOND_ORDER, new_order.clone(), true)?;
    let n_swaps = mol.perturbation_order(old_atom, &new_order)?;
    let tag = mol.atom(old_atom).chiral_tag();
    let atom = res.atom_mut(new_atom);
    atom.set_chiral_tag(tag);
    if n_swaps % 2 != 0 {
        atom.invert_chirality();
    }
    Ok(())
}

/// The bonds next to `bond` that carry stereo atoms, with those atoms. Their
/// stereo is saved in a property and restored once the stereo atoms are
/// updated.
fn get_nbr_bond_stereo(mol: &mut RWMol, bond: u32) -> Result<Vec<(u32, Vec<i32>)>> {
    let (bgn, end) = (
        mol.bond(bond).begin_atom_idx(),
        mol.bond(bond).end_atom_idx(),
    );
    let mut res = Vec::new();
    for atom in [bgn, end] {
        let bonds: Vec<u32> = mol.atom_bond_indices(atom).collect();
        for obnd in bonds {
            if obnd != bond && !mol.bond(obnd).stereo_atoms().is_empty() {
                let stereo = mol.bond(obnd).stereo() as i32;
                mol.bond_mut(obnd)
                    .set_prop(MOLFRAG_SAVE_STEREO, stereo, false)?;
                res.push((obnd, mol.bond(obnd).stereo_atoms().to_vec()));
            }
        }
    }
    Ok(res)
}

/// Replaces the bond `bond_idx` by one carrying the query of the bond that is
/// being cut.
fn set_bond_query(
    res: &mut RWMol,
    bond_idx: u32,
    query: &crate::graphmol::queries::BondQuery,
) -> Result<()> {
    let mut qb = res.bond(bond_idx).copy_detached();
    qb.set_query(Some(query.clone()));
    res.replace_bond(bond_idx, &qb, false, true)
}

/// `fragmentOnBonds(mol, bondIndices, addDummies, dummyLabels, bondTypes,
/// nCutsPerAtom)`: a copy of `mol` with the given bonds broken. With
/// `add_dummies` each broken bond is replaced by two bonds to new dummy atoms
/// (labelled with `dummy_labels`, or with the indices of the atoms they
/// replace; `bond_types` gives the new bond types, default the broken bond's
/// type). `n_cuts_per_atom`, when given, receives the number of cuts at each
/// atom.
pub fn fragment_on_bonds(
    mol: &ROMol,
    bond_indices: &[u32],
    add_dummies: bool,
    dummy_labels: Option<&[(u32, u32)]>,
    bond_types: Option<&[BondType]>,
    mut n_cuts_per_atom: Option<&mut [u32]>,
) -> Result<RWMol> {
    check_per_bond_args(bond_indices, dummy_labels, bond_types)?;
    if n_cuts_per_atom
        .as_deref()
        .is_some_and(|n| n.len() != mol.num_atoms() as usize)
    {
        return Err(RDKitError::precondition("bad nCutsPerAtom vector"));
    }
    if let Some(n) = n_cuts_per_atom.as_deref_mut() {
        n.fill(0);
    }
    let mut res = mol.clone();
    if mol.num_atoms() == 0 || bond_indices.is_empty() {
        return Ok(res);
    }
    for &b in bond_indices {
        res.get_bond(b)?;
    }

    let mut atoms_to_update: BTreeMap<u32, u32> = BTreeMap::new();
    res.begin_batch_edit()?;
    for (i, &bond_idx) in bond_indices.iter().enumerate() {
        let bond = res.bond(bond_idx);
        let (bidx, eidx) = (bond.begin_atom_idx(), bond.end_atom_idx());
        let mut bond_type = bond.bond_type();
        let bond_dir = bond.bond_dir();
        let nbr_bond_stereo = get_nbr_bond_stereo(&mut res, bond_idx)?;
        // grab a copy of any query on the outgoing bond if it will be needed
        let out_query = if bond_types.is_none() {
            res.bond(bond_idx).query().cloned()
        } else {
            None
        };
        res.remove_bond(bidx, eidx)?;
        if let Some(n) = n_cuts_per_atom.as_deref_mut() {
            n[bidx as usize] += 1;
            n[eidx as usize] += 1;
        }
        if add_dummies {
            let (label1, label2) = dummy_labels.map_or((bidx, eidx), |d| d[i]);
            let mut at1 = Atom::new(0);
            at1.set_isotope(label1);
            let mut at2 = Atom::new(0);
            at2.set_isotope(label2);

            let idx1 = res.add_atom(at1, false);
            *atoms_to_update.entry(idx1).or_default() += 1;
            if let Some(types) = bond_types {
                bond_type = types[i];
            }
            let new_bond = res.add_bond(idx1, eidx, bond_type)? - 1;
            // the dummy replaces the original start atom, so the direction
            // will be ok as long as it's one of the states associated with
            // double bond stereo
            if matches!(bond_dir, BondDir::EndDownRight | BondDir::EndUpRight) {
                res.bond_mut(new_bond).set_bond_dir(bond_dir);
            }
            if let Some(q) = &out_query {
                set_bond_query(&mut res, new_bond, q)?;
            }

            let idx2 = res.add_atom(at2, false);
            *atoms_to_update.entry(idx2).or_default() += 1;
            let new_bond = res.add_bond(bidx, idx2, bond_type)? - 1;
            // this bond starts at the same atom, so its direction should
            // always be correct
            res.bond_mut(new_bond).set_bond_dir(bond_dir);
            if let Some(q) = &out_query {
                set_bond_query(&mut res, new_bond, q)?;
            }

            // restore the stereo of the neighboring bonds
            for (nbr_bond, mut stereo_atoms) in nbr_bond_stereo {
                for a in &mut stereo_atoms {
                    if *a == bidx as i32 {
                        *a = idx1 as i32;
                    }
                }
                for a in &mut stereo_atoms {
                    if *a == eidx as i32 {
                        *a = idx2 as i32;
                    }
                }
                let saved = res.bond(nbr_bond).get_prop(MOLFRAG_SAVE_STEREO).cloned();
                let b = res.bond_mut(nbr_bond);
                *b.stereo_atoms_mut() = stereo_atoms;
                if let Some(RDValue::Int(s)) = saved {
                    b.set_stereo(BondStereo::from_u8(s as u8).unwrap_or_default())?;
                }
                b.clear_prop(MOLFRAG_SAVE_STEREO);
            }

            // figure out if we need to change the stereo tags on the atoms
            let connecting = mol
                .bond_idx_between(bidx, eidx)
                .ok_or_else(|| RDKitError::invariant("bond to cut not found"))?;
            for atom in [bidx, eidx] {
                if matches!(
                    mol.atom(atom).chiral_tag(),
                    ChiralType::TetrahedralCCW | ChiralType::TetrahedralCW
                ) {
                    check_chirality_post_move(mol, atom, &mut res, atom, connecting)?;
                }
            }

            for conf in res.conformers_mut() {
                let (p1, p2) = (conf.atom_pos(bidx)?, conf.atom_pos(eidx)?);
                conf.set_atom_pos(idx1, p1)?;
                conf.set_atom_pos(idx2, p2)?;
            }
        }
        // keep track of the atoms so that we can adjust H counts later
        *atoms_to_update.entry(bidx).or_default() += 1;
        *atoms_to_update.entry(eidx).or_default() += 1;
    }
    res.commit_batch_edit()?;
    res.clear_computed_props(true);

    // adjust H counts on dummies and atoms where bonds were broken
    // (github issues 429, 6034)
    for (idx, count) in atoms_to_update {
        let atom = res.atom_mut(idx);
        if !add_dummies && (atom.no_implicit() || (atom.is_aromatic() && atom.atomic_num() != 6)) {
            let hs = atom.num_explicit_hs() + count;
            atom.set_num_explicit_hs(hs);
        }
        res.update_atom_property_cache(idx, false)?;
    }
    Ok(res)
}

/// `fragmentOnBonds(mol, bondPatterns, atomEnvirons, nCutsPerAtom)`: breaks
/// the bonds matched by the bond definitions, the more specific (later)
/// definitions first, labelling the dummies as the definitions say. The
/// definitions are only used when all their `atom_environs` match somewhere
/// in `mol`.
pub fn fragment_on_bond_patterns(
    mol: &ROMol,
    bond_patterns: &[FragmenterBondType],
    atom_environs: Option<&BTreeMap<u32, ROMol>>,
    n_cuts_per_atom: Option<&mut [u32]>,
) -> Result<RWMol> {
    if n_cuts_per_atom
        .as_deref()
        .is_some_and(|n| n.len() != mol.num_atoms() as usize)
    {
        return Err(RDKitError::precondition("bad nCutsPerAtom vector"));
    }
    let first_match = SubstructMatchParameters {
        max_matches: 1,
        ..SubstructMatchParameters::default()
    };
    let mut environs_match: BTreeMap<u32, bool> = BTreeMap::new();
    if let Some(environs) = atom_environs {
        for (&idx, env) in environs {
            environs_match.insert(idx, !substruct_match(mol, env, &first_match)?.is_empty());
        }
    }

    let mut bond_indices = Vec::new();
    let mut dummy_labels = Vec::new();
    let mut bond_types = Vec::new();
    let mut bonds_used = BTreeSet::new();
    // the bond definitions are organized (more or less) general -> specific,
    // so loop over them backwards
    for fbt in bond_patterns.iter().rev() {
        let Some(query) = fbt.query.as_deref() else {
            continue;
        };
        if query.num_atoms() != 2 || query.num_bonds() != 1 {
            rd_log_error!("fragmentation queries must have 2 atoms and 1 bond\n");
            continue;
        }
        if atom_environs.is_some()
            && (!environs_match
                .get(&fbt.atom1_type)
                .copied()
                .unwrap_or(false)
                || !environs_match
                    .get(&fbt.atom2_type)
                    .copied()
                    .unwrap_or(false))
        {
            continue;
        }
        for mv in substruct_match(mol, query, &SubstructMatchParameters::default())? {
            let bond = mol
                .bond_idx_between(mv[0].1 as u32, mv[1].1 as u32)
                .ok_or_else(|| RDKitError::invariant("bond between matched atoms not found"))?;
            if !bonds_used.insert(bond) {
                continue;
            }
            bond_indices.push(bond);
            dummy_labels.push(if mol.bond(bond).begin_atom_idx() == mv[0].1 as u32 {
                (fbt.atom1_label, fbt.atom2_label)
            } else {
                (fbt.atom2_label, fbt.atom1_label)
            });
            bond_types.push(fbt.bond_type);
        }
    }
    fragment_on_bonds(
        mol,
        &bond_indices,
        true,
        Some(&dummy_labels),
        Some(&bond_types),
        n_cuts_per_atom,
    )
}

/// The BRICS bond and atom-environment definitions, built on first use.
struct BricsDefinitions {
    bond_patterns: Vec<FragmenterBondType>,
    atom_environs: BTreeMap<u32, ROMol>,
}

fn brics_definitions() -> Result<&'static BricsDefinitions> {
    static DEFS: OnceLock<BricsDefinitions> = OnceLock::new();
    if let Some(d) = DEFS.get() {
        return Ok(d);
    }
    let mut atom_environs = BTreeMap::new();
    construct_brics_atom_types(Some(&mut atom_environs))?;
    let bond_patterns = construct_brics_bond_types()?;
    Ok(DEFS.get_or_init(|| BricsDefinitions {
        bond_patterns,
        atom_environs,
    }))
}

/// `fragmentOnBRICSBonds(mol)`: a copy of `mol` with all BRICS bonds broken.
pub fn fragment_on_brics_bonds(mol: &ROMol) -> Result<RWMol> {
    let defs = brics_definitions()?;
    fragment_on_bond_patterns(mol, &defs.bond_patterns, Some(&defs.atom_environs), None)
}
