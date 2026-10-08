//! Labute ASA / TPSA / VSA tests; expected values from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::descriptors::mol_surf::{
    calc_labute_asa, calc_slogp_vsa, calc_tpsa, get_labute_atom_contribs, get_tpsa_atom_contribs,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn smi(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}

#[test]
fn tpsa_with_and_without_sulfur_and_phosphorus() {
    close(
        calc_tpsa(&mut smi("CS(=O)(=O)Nc1ccccc1"), false, false).unwrap(),
        46.17,
    );
    close(
        calc_tpsa(&mut smi("CS(=O)(=O)Nc1ccccc1"), false, true).unwrap(),
        54.55,
    );
    close(
        calc_tpsa(&mut smi("OP(=O)(O)O"), false, true).unwrap(),
        87.57,
    );
    close(calc_tpsa(&mut smi("C1CN1"), false, false).unwrap(), 21.94);
}

#[test]
fn tpsa_contributions_are_per_atom() {
    let (contribs, total) =
        get_tpsa_atom_contribs(&mut smi("CC(=O)Nc1ccc(O)cc1"), false, false).unwrap();
    assert_eq!(contribs.len(), 11);
    close(contribs[2], 17.07);
    close(contribs[3], 12.03);
    close(total, 49.33);
}

#[test]
fn labute_asa_is_cached_like_upstream() {
    let mut m = smi("CCO");
    let with_hs = calc_labute_asa(&mut m, true, false).unwrap();
    close(with_hs, 19.89842689442217);
    // the cached value wins whatever includeHs says, unless forced
    assert_eq!(calc_labute_asa(&mut m, false, false).unwrap(), with_hs);
    close(
        calc_labute_asa(&mut m, false, true).unwrap(),
        18.682762484374372,
    );
    close(
        get_labute_atom_contribs(&mut smi("CCO"), true, false)
            .unwrap()
            .h,
        1.2612803353779227,
    );
}

#[test]
fn vsa_bins_cover_all_contributions() {
    let mut m = smi("CC(=O)Nc1ccc(O)cc1");
    let bins = calc_slogp_vsa(&mut m, Some(&[0.0, 0.5]), false).unwrap();
    assert_eq!(bins.len(), 3);
    // the bins hold the atom contributions; the implicit-H term is separate
    let atoms = get_labute_atom_contribs(&mut m, true, false).unwrap().atoms;
    close(bins.iter().sum::<f64>(), atoms.iter().sum::<f64>());
}
