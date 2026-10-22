//! Differential test of `CompareReversedINChI3` (`ichirvr7.rs`) against the C
//! library on random pairs of INChI structures, with the fixture from
//! `tools/inchi_cmp3_reference.py`.

use rdoxide::external::inchi::ichi::{InChI, InChIIsotopicAtom, InChIStereo};
use rdoxide::external::inchi::ichirvr7::compare_reversed_inchi3;

#[allow(clippy::field_reassign_with_default)]
fn parse(line: &str) -> Option<InChI> {
    let mut it = line.split_whitespace();
    let mut n = || it.next().unwrap();
    if n() == "0" {
        return None;
    }
    let mut p = InChI::default();
    p.n_error_code = n().parse().unwrap();
    let na: usize = n().parse().unwrap();
    p.n_number_of_atoms = na as i32;
    p.n_atom = (0..na).map(|_| n().parse().unwrap()).collect();
    p.n_num_h = (0..na).map(|_| n().parse().unwrap()).collect();
    if n() == "1" {
        p.n_num_h_fixed = Some((0..na).map(|_| n().parse().unwrap()).collect());
    }
    let hill = n();
    p.sz_hill_formula = (hill != "-").then(|| hill.to_string());
    p.len_conn_table = n().parse().unwrap();
    p.n_conn_table = (0..p.len_conn_table)
        .map(|_| n().parse().unwrap())
        .collect();
    p.len_tautomer = n().parse().unwrap();
    p.n_tautomer = (0..p.len_tautomer).map(|_| n().parse().unwrap()).collect();
    p.n_number_of_isotopic_atoms = n().parse().unwrap();
    p.isotopic_atom = (0..p.n_number_of_isotopic_atoms)
        .map(|_| InChIIsotopicAtom {
            n_atom_number: n().parse().unwrap(),
            n_iso_difference: n().parse().unwrap(),
            n_num_h: n().parse().unwrap(),
            n_num_d: n().parse().unwrap(),
            n_num_t: n().parse().unwrap(),
        })
        .collect();
    p.n_total_charge = n().parse().unwrap();
    for w in 0..2 {
        if n() == "0" {
            continue;
        }
        let mut s = InChIStereo::default();
        s.n_number_of_stereo_centers = n().parse().unwrap();
        for _ in 0..s.n_number_of_stereo_centers {
            s.n_number.push(n().parse().unwrap());
            s.t_parity.push(n().parse().unwrap());
        }
        s.n_comp_inv2_abs = n().parse().unwrap();
        s.n_number_of_stereo_bonds = n().parse().unwrap();
        for _ in 0..s.n_number_of_stereo_bonds {
            s.n_bond_atom1.push(n().parse().unwrap());
            s.n_bond_atom2.push(n().parse().unwrap());
            s.b_parity.push(n().parse().unwrap());
        }
        if w == 0 {
            p.stereo = Some(s);
        } else {
            p.stereo_isotopic = Some(s);
        }
    }
    Some(p)
}

#[test]
fn compare_reversed_inchi3_matches_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_cmp3.json")).unwrap();
    let cases = data.as_array().unwrap();
    let mut bad = 0;
    for c in cases {
        let i = c["in"].as_array().unwrap();
        let (a, b) = (parse(i[0].as_str().unwrap()), parse(i[1].as_str().unwrap()));
        let got = format!(
            "{} 0",
            compare_reversed_inchi3(a.as_ref(), b.as_ref(), None, None)
        );
        if got != c["out"].as_str().unwrap() {
            bad += 1;
            if bad < 4 {
                eprintln!("{i:?}\n got {got}\n want {}", c["out"]);
            }
        }
    }
    assert_eq!(bad, 0, "{bad} of {} cases differ", cases.len());
}
