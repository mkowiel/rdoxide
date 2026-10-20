//! The RDKit side of `MolToInchi`: the conversion of a molecule into the
//! `inchi_Atom` / `inchi_Stereo0D` arrays. The `rA:`/`rB:` parts of the
//! AuxInfo that the official RDKit 2026.09.1 / InChI 1.07.3 produced for the
//! corpus (`tests/data/inchi_corpus.jsonl`) record the structure that the
//! InChI library received; the elements, charges and bonds are compared with
//! them. Hand-written cases check the stereo and isotope encoding.

use std::path::PathBuf;

use rdoxide::external::inchi::inchi_dll::extract_atoms_and_bonds;
use rdoxide::external::inchi::readinch::{
    INCHI_PARITY_EVEN, INCHI_PARITY_ODD, INCHI_STEREO_TYPE_DOUBLE_BOND,
    INCHI_STEREO_TYPE_TETRAHEDRAL,
};
use rdoxide::external::inchi_input::{fix_option_symbol, mol_to_inchi_input};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use serde_json::Value;

fn input(smiles: &str) -> rdoxide::external::inchi_input::InchiInput {
    let mol = mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    mol_to_inchi_input(&mol, None).unwrap().unwrap()
}

#[test]
fn option_symbols() {
    assert_eq!(fix_option_symbol("/FixedH /RecMet"), "-FixedH -RecMet");
    assert_eq!(fix_option_symbol("-SNon"), "-SNon");
    let mol = mol_from_smiles("CC", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    let inp = mol_to_inchi_input(&mol, Some("/SUU")).unwrap().unwrap();
    assert_eq!(inp.options.as_deref(), Some("-SUU"));
}

#[test]
fn atoms_and_isotopes() {
    let i = input("[13CH3][2H]");
    assert_eq!(i.atoms.len(), 2);
    assert_eq!(i.atoms[0].elname_str(), "C");
    // mass 13 vs average 12.011 -> 12
    assert_eq!(i.atoms[0].isotopic_mass, 10001);
    assert_eq!(i.atoms[0].num_iso_h[0], -1);
    assert_eq!(i.atoms[1].elname_str(), "H");
    assert_eq!(i.atoms[1].isotopic_mass, 10001); // 2 - round(1.008) = 1
    assert_eq!(i.atoms[0].num_bonds, 1);
    assert_eq!(i.atoms[0].neighbor[0], 1);
    assert_eq!(i.atoms[0].bond_type[0], 1);
    assert_eq!(i.atoms[1].num_bonds, 0);
}

#[test]
fn charges_hydrogens_radicals() {
    let i = input("[NH4+]");
    assert_eq!(i.atoms[0].charge, 1);
    assert_eq!(i.atoms[0].num_iso_h[0], -1); // nitrogen: H added by InChI
    let i = input("[Na+]");
    assert_eq!(i.atoms[0].num_iso_h[0], 0); // other elements: explicit count
    let i = input("[CH3]");
    assert_eq!(i.atoms[0].num_iso_h[0], 3); // radical: the H count is given
}

#[test]
fn double_bond_stereo() {
    let i = input("C/C=C/C");
    assert_eq!(i.stereo_0d.len(), 1);
    let s = i.stereo_0d[0];
    assert_eq!(s.type_, INCHI_STEREO_TYPE_DOUBLE_BOND);
    assert_eq!(s.central_atom, -1);
    assert_eq!(s.neighbor, [0, 1, 2, 3]);
    assert_eq!(i32::from(s.parity), INCHI_PARITY_EVEN); // trans
    let i = input("C/C=C\\C");
    assert_eq!(i32::from(i.stereo_0d[0].parity), INCHI_PARITY_ODD); // cis
}

#[test]
fn tetrahedral_stereo() {
    let i = input("C[C@H](N)C(=O)O");
    assert_eq!(i.stereo_0d.len(), 1);
    let s = i.stereo_0d[0];
    assert_eq!(s.type_, INCHI_STEREO_TYPE_TETRAHEDRAL);
    assert_eq!(s.central_atom, 1);
    // three explicit neighbors: the center itself stands for the implicit H
    assert_eq!(s.neighbor, [1, 0, 2, 3]);
    let j = input("C[C@@H](N)C(=O)O");
    assert_ne!(i.stereo_0d[0].parity, j.stereo_0d[0].parity);
    // four explicit neighbors
    let k = input("C[C@](N)(O)F");
    assert_eq!(k.stereo_0d[0].neighbor, [0, 2, 3, 4]);
}

#[test]
fn bond_types_and_the_bond_limit() {
    let i = input("C#CC=O");
    assert_eq!(i.atoms[0].bond_type[0], 3);
    assert_eq!(i.atoms[2].bond_type[0], 2);
    // an atom with more than MAXVAL (20) bonds: no input
    let smiles = format!("C{}", "(C)".repeat(21));
    let mol = mol_from_smiles(
        &smiles,
        &SmilesParserParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap()
    .unwrap();
    assert!(mol_to_inchi_input(&mol, None).unwrap().is_none());
}

fn corpus() -> Vec<Value> {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/inchi_corpus.jsonl");
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}

/// Elements and charges of the `rA:` part of an AuxInfo
fn parse_ra(aux: &str) -> Option<Vec<(String, i32)>> {
    let ra = aux.split("/rA:").nth(1)?.split("/rB:").next()?;
    let digits = ra.chars().take_while(|c| c.is_ascii_digit()).count();
    let n: usize = ra[..digits].parse().ok()?;
    let b: Vec<char> = ra[digits..].chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if !b[i].is_ascii_uppercase() {
            return None;
        }
        let mut el = b[i].to_string();
        i += 1;
        if i < b.len()
            && b[i].is_ascii_lowercase()
            && !(b[i] == 'i' && i + 1 < b.len() && b[i + 1].is_ascii_digit() && false)
        {
            el.push(b[i]);
            i += 1;
        }
        let mut charge = 0;
        while i < b.len() && (b[i] == '+' || b[i] == '-') {
            let sign = if b[i] == '+' { 1 } else { -1 };
            i += 1;
            let d: String = b[i..].iter().take_while(|c| c.is_ascii_digit()).collect();
            i += d.len();
            charge += sign * d.parse::<i32>().unwrap_or(1);
        }
        // annotations: ".i13" isotope, ".e" / ".o" parity, ".m"/".u"/".?" ...
        while i < b.len() && b[i] == '.' {
            i += 1;
            while i < b.len() && !b[i].is_ascii_uppercase() && b[i] != '.' {
                i += 1;
            }
        }
        out.push((el, charge));
    }
    (out.len() == n).then_some(out)
}

#[test]
fn input_matches_the_structure_the_inchi_library_received() {
    let mut checked = 0;
    for r in corpus() {
        if r["kind"] != "mol" || r["options"] != "" {
            continue;
        }
        let smiles = r["smiles"].as_str().unwrap();
        let aux = r["aux"].as_str().unwrap();
        let Some(ra) = parse_ra(aux) else { continue };
        let mol = mol_from_smiles(smiles, &SmilesParserParams::default())
            .unwrap()
            .unwrap();
        let inp = mol_to_inchi_input(&mol, None).unwrap().unwrap();
        assert_eq!(inp.atoms.len(), ra.len(), "{smiles}");
        for (a, (el, charge)) in inp.atoms.iter().zip(&ra) {
            assert_eq!(&a.elname_str(), el, "{smiles}");
            assert_eq!(i32::from(a.charge), *charge, "{smiles}");
        }
        // bonds: segment k of rB lists the bonds of atom k + 2 to lower atoms
        let rb = aux
            .split("/rB:")
            .nth(1)
            .unwrap()
            .split("/rC:")
            .next()
            .unwrap();
        let mut want: Vec<(usize, usize, char)> = Vec::new();
        for (k, seg) in rb.split(';').enumerate() {
            // bonds of one atom follow each other: "s2s7" is two bonds
            let chars: Vec<char> = seg.chars().collect();
            let mut i = 0;
            while i < chars.len() {
                let ty = chars[i];
                i += 1;
                while i < chars.len() && (chars[i] == '+' || chars[i] == '-') {
                    i += 1;
                }
                let num: String = chars[i..]
                    .iter()
                    .take_while(|c| c.is_ascii_digit())
                    .collect();
                i += num.len();
                let lower: usize = num.parse::<usize>().unwrap() - 1;
                want.push((lower, k + 1, ty));
            }
        }
        let mut have: Vec<(usize, usize, char)> = Vec::new();
        for (i, a) in inp.atoms.iter().enumerate() {
            for j in 0..a.num_bonds as usize {
                let ty = match a.bond_type[j] {
                    1 => 's',
                    2 => 'd',
                    3 => 't',
                    _ => '?',
                };
                have.push((i, a.neighbor[j] as usize, ty));
            }
        }
        want.sort();
        have.sort();
        assert_eq!(have, want, "{smiles}");
        // the structure goes through the library's own input conversion
        let mut msg = String::new();
        let mut err = 0;
        let (at, _, n_bonds, _) = extract_atoms_and_bonds(&inp.atoms, false, &mut msg, &mut err);
        assert_eq!(err, 0, "{smiles}: {msg}");
        assert_eq!(n_bonds as usize, have.len(), "{smiles}");
        assert_eq!(at.len(), inp.atoms.len());
        checked += 1;
    }
    assert!(checked > 150, "only {checked} records compared");
}
