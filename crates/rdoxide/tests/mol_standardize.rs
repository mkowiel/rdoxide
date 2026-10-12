//! MolStandardize; inputs from upstream `testMolStandardize.py` and the C++
//! tests, expectations checked against RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::mol_standardize::charge::{Reionizer, Uncharger};
use rdoxide::graphmol::mol_standardize::fragment::{FragmentRemover, LargestFragmentChooser};
use rdoxide::graphmol::mol_standardize::metal::MetalDisconnector;
use rdoxide::graphmol::mol_standardize::normalize::Normalizer;
use rdoxide::graphmol::mol_standardize::pipeline::{Pipeline, status};
use rdoxide::graphmol::mol_standardize::tautomer::{
    TautomerEnumerator, TautomerEnumeratorStatus, default_tautomer_score_substructs,
    score_hetero_hs, score_rings, score_substructs, score_tautomer,
};
use rdoxide::graphmol::mol_standardize::validate::validate_smiles;
use rdoxide::graphmol::mol_standardize::{
    self as ms, CleanupParameters, default_cleanup_parameters, standardize_multiple_mols_in_place,
    update_cleanup_params_from_json,
};
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn mol(smiles: &str) -> RWMol {
    mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn smi(mol: &RWMol) -> String {
    mol_to_smiles(&mut mol.clone(), &SmilesWriteParams::default()).unwrap()
}

type Step = fn(&RWMol) -> RWMol;

fn check(step: Step, cases: &[(&str, &str)]) {
    for (input, expected) in cases {
        assert_eq!(smi(&step(&mol(input))), *expected, "{input}");
    }
}

fn params() -> &'static CleanupParameters {
    default_cleanup_parameters()
}

#[test]
fn cleanup() {
    check(
        |m| ms::cleanup(m, params()).unwrap(),
        &[
            ("CCC(=O)O[Na]", "CCC(=O)[O-].[Na+]"),
            ("[Na]OC(=O)c1ccccc1", "O=C([O-])c1ccccc1.[Na+]"),
            ("C[S+2]([O-])([O-])O", "CS(=O)(=O)O"),
            ("CN(=O)=O", "C[N+](=O)[O-]"),
            ("O=C([O-])c1ccccc1.[Na+]", "O=C([O-])c1ccccc1.[Na+]"),
        ],
    );
}

#[test]
fn standardize_smiles() {
    assert_eq!(
        ms::standardize_smiles("CCC(=O)O[Na]").unwrap(),
        "CCC(=O)[O-].[Na+]"
    );
    assert_eq!(
        ms::standardize_smiles("[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1").unwrap(),
        "O=C([O-])c1ccc(C[S](=O)=O)cc1.[Na+]"
    );
    assert!(ms::standardize_smiles("C1CC").is_err());
}

#[test]
fn parents() {
    check(
        |m| ms::fragment_parent(m, params(), false).unwrap(),
        &[
            ("[Na]OC(=O)c1ccccc1", "O=C([O-])c1ccccc1"),
            ("[Cl-].[NH4+].CCCCCC", "CCCCCC"),
            ("c1ccccc1.CC(=O)O", "c1ccccc1"),
        ],
    );
    check(
        |m| ms::charge_parent(m, params(), false).unwrap(),
        &[
            ("C[NH+](C)C.[Cl-]", "CN(C)C"),
            ("[O-]C(=O)c1ccccc1.[Na+]", "O=C(O)c1ccccc1"),
            ("C[N+](C)(C)CC(=O)[O-]", "C[N+](C)(C)CC(=O)[O-]"),
        ],
    );
    check(
        |m| ms::isotope_parent(m, params(), false).unwrap(),
        &[("[13CH3]C(=O)O", "CC(=O)O")],
    );
    check(
        |m| ms::stereo_parent(m, params(), false).unwrap(),
        &[("C[C@H](F)Cl", "CC(F)Cl")],
    );
    check(
        |m| ms::super_parent(m, params(), false).unwrap(),
        &[
            ("[O-]C(=O)C(C)Oc1ccccc1.[Na+]", "CC(Oc1ccccc1)C(=O)O"),
            ("O=C(O)C[13C@H](N)CC(O)=O", "NC(CC(=O)O)CC(=O)O"),
        ],
    );
}

#[test]
fn fragments() {
    check(
        |m| ms::remove_fragments(m, params()).unwrap(),
        &[
            ("CN(C)C.Cl.Cl.Br", "CN(C)C"),
            ("[Na+].[Cl-]", "[Na+]"),
            ("CCCC.Cl", "CCCC"),
        ],
    );
    let skip = {
        let mut r = FragmentRemover::default();
        r.skip_if_all_match = true;
        r
    };
    assert_eq!(
        smi(&skip.remove(&mol("[Na+].[Cl-]")).unwrap()),
        "[Cl-].[Na+]"
    );
    let custom = FragmentRemover::from_text("// comment\nsodium\t[Na+]\n", true, false).unwrap();
    assert_eq!(smi(&custom.remove(&mol("CC.[Na+]")).unwrap()), "CC");
    let chooser = LargestFragmentChooser::new(true);
    assert_eq!(
        smi(&chooser
            .choose(&mol("O=C(O)c1ccccc1.[Pt](Cl)(Cl)(Cl)(Cl)(Cl)Cl"))
            .unwrap()),
        "O=C(O)c1ccccc1"
    );
}

#[test]
fn normalize() {
    check(
        |m| ms::normalize(m, params()).unwrap(),
        &[
            ("CS(C)=O", "C[S+](C)[O-]"),
            ("CN=N#N", "CN=[N+]=[N-]"),
            ("C[n+]1ccccc1[O-]", "Cn1ccccc1=O"),
            ("CS(=O)(=O)O", "CS(=O)(=O)O"),
        ],
    );
    let mut m = mol("CS(C)=O");
    Normalizer::default().normalize_in_place(&mut m).unwrap();
    assert_eq!(smi(&m), "C[S+](C)[O-]");
}

#[test]
fn charge() {
    check(
        |m| ms::reionize(m, params()).unwrap(),
        &[
            (
                "C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O",
                "O=S(O)c1ccc(S(=O)(=O)[O-])cc1",
            ),
            ("[Na].O=C(O)c1ccccc1", "O=C([O-])c1ccccc1.[Na+]"),
        ],
    );
    check(
        |m| Uncharger::default().uncharge(m).unwrap(),
        &[
            ("C[NH+](C)C", "CN(C)C"),
            ("[O-]C(=O)C[NH3+]", "NCC(=O)O"),
            ("C[N+](C)(C)CC(=O)[O-]", "C[N+](C)(C)CC(=O)[O-]"),
            ("[O-]c1ccccc1", "Oc1ccccc1"),
        ],
    );
    let r =
        Reionizer::from_text("carboxylic acid\tC(=O)[OH]\tC(=O)[O-]\n", Some(Vec::new())).unwrap();
    assert_eq!(r.params().pairs().len(), 1);
}

#[test]
fn metals() {
    check(
        |m| ms::disconnect_organometallics(m, None).unwrap(),
        &[
            ("CCC(=O)O[Na]", "CCC(=O)O.[Na]"),
            ("[Mg](Cl)CC", "CC.Cl.[Mg]"),
            ("CN(C)C(=O)[Fe](Cl)Cl", "CN(C)C=O.Cl.Cl.[Fe]"),
        ],
    );
    let md = MetalDisconnector::default();
    assert_eq!(
        smi(&md.disconnect(&mol("CCC(=O)O[Na]")).unwrap()),
        "CCC(=O)[O-].[Na+]"
    );
}

#[test]
fn validation() {
    for (input, expected) in [
        (
            "ClCCCl.c1ccccc1O",
            vec!["INFO: [FragmentValidation] 1,2-dichloroethane is present"],
        ),
        (
            "[2H]C(Cl)(Cl)Cl",
            vec!["INFO: [IsotopeValidation] Molecule contains isotope 2H"],
        ),
        (
            "C[N+](C)(C)C",
            vec!["INFO: [NeutralValidation] Not an overall neutral system (+1)"],
        ),
        (
            "O=C([O-])c1ccccc1",
            vec!["INFO: [NeutralValidation] Not an overall neutral system (-1)"],
        ),
        ("", vec!["ERROR: [NoAtomValidation] Molecule has no atoms"]),
    ] {
        assert_eq!(validate_smiles(input).unwrap(), expected, "{input}");
    }
}

#[test]
fn tautomers() {
    check(
        |m| ms::canonical_tautomer(m, params()).unwrap(),
        &[
            ("Oc1ccccn1", "O=c1cccc[nH]1"),
            ("CC(O)=CC(C)=O", "CC(=O)CC(C)=O"),
            ("N=c1[nH]cccc1", "Nc1ccccn1"),
            ("C1(=CCCCC1)O", "O=C1CCCCC1"),
        ],
    );
    let te = TautomerEnumerator::default();
    let res = te.enumerate(&mol("Oc1ccccn1")).unwrap();
    assert_eq!(res.smiles(), ["O=C1CC=CC=N1", "O=c1cccc[nH]1", "Oc1ccccn1"]);
    assert_eq!(res.modified_atom_indices(), [0, 2, 6]);
    assert_eq!(res.modified_bond_indices(), [0, 1, 2, 3, 4, 5, 6]);
    assert_eq!(res.status, TautomerEnumeratorStatus::Completed);
    let res = te.enumerate(&mol("CC(=O)CC(=O)C")).unwrap();
    assert_eq!(
        res.smiles(),
        [
            "C=C(O)C=C(C)O",
            "C=C(O)CC(=C)O",
            "C=C(O)CC(C)=O",
            "CC(=O)C=C(C)O",
            "CC(=O)CC(C)=O"
        ]
    );
    assert_eq!(res.modified_atom_indices(), [0, 2, 3, 5, 6]);
    assert_eq!(
        TautomerEnumerator::v1()
            .enumerate(&mol("CC(=O)CC(=O)C"))
            .unwrap()
            .len(),
        5
    );
    let mut limited = TautomerEnumerator::default();
    limited.max_tautomers = 2;
    let res = limited.enumerate(&mol("CC(=O)CC(=O)C")).unwrap();
    assert_eq!(res.status, TautomerEnumeratorStatus::MaxTautomersReached);
}

#[test]
fn tautomer_scores() {
    assert_eq!(score_tautomer(&mol("Oc1ccccn1")).unwrap(), 100);
    assert_eq!(score_rings(&mol("c1ccccc1C1=CC=CC=N1")).unwrap(), 350);
    assert_eq!(score_hetero_hs(&mol("CS")).unwrap(), -1);
    assert_eq!(
        score_substructs(
            &mol("CC\\C=C(/O)[C@@H](C)C(C)=O"),
            default_tautomer_score_substructs()
        )
        .unwrap(),
        6
    );
}

#[test]
fn json_parameters() {
    let mut p = CleanupParameters::default();
    update_cleanup_params_from_json(
        &mut p,
        r#"{"maxRestarts": 5, "preferOrganic": true,
            "fragmentData": [{"name": "hydrogen", "smarts": "[H]"}, {"name": "", "smarts": "C"}],
            "tautomerTransformData": [{"name": "1,3 keto/enol", "smarts": "[CX4!H0]-[C]=[O]", "bonds": "", "charges": ""}]}"#,
    )
    .unwrap();
    assert_eq!(p.max_restarts, 5);
    assert!(p.prefer_organic);
    assert_eq!(
        p.fragment_data,
        [("hydrogen".to_string(), "[H]".to_string())]
    );
    assert_eq!(p.tautomer_transform_data.len(), 1);
    let te = TautomerEnumerator::from_params(&p).unwrap();
    assert_eq!(te.params().transforms().len(), 1);
}

#[test]
fn multithreaded_in_place() {
    let inputs = [
        "CCC(=O)O[Na]",
        "[Na]OC(=O)c1ccccc1",
        "CN(=O)=O",
        "C[NH+](C)C.[Cl-]",
    ];
    let mut mols: Vec<RWMol> = inputs.iter().map(|s| mol(s)).collect();
    let mut refs: Vec<&mut RWMol> = mols.iter_mut().collect();
    standardize_multiple_mols_in_place(&mut refs, 4, |m| ms::cleanup_in_place(m, params()))
        .unwrap();
    let expected: Vec<String> = inputs
        .iter()
        .map(|s| smi(&ms::cleanup(&mol(s), params()).unwrap()))
        .collect();
    assert_eq!(mols.iter().map(smi).collect::<Vec<_>>(), expected);
}

#[test]
fn pipeline() {
    let molblock = "
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000000 0.000000 0.000000 0
M  V30 2 C 1.299038 0.750000 0.000000 0
M  V30 3 O 2.598076 -0.000000 0.000000 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let res = Pipeline::default().run(molblock).unwrap();
    assert_eq!(res.status & status::PIPELINE_ERROR, 0);
    assert_ne!(res.status & status::PROTONATION_CHANGED, 0);
    assert!(res.output_mol_data.contains("M  V30 3 O"));
    let bad = Pipeline::default().run("garbage").unwrap();
    assert_ne!(bad.status & status::INPUT_ERROR, 0);
}
