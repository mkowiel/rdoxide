//! R-group decomposition (`graphmol::rgroup_decomposition`); cases ported
//! from upstream `testRGroupDecomp.cpp`.

use rdoxide::graphmol::ROMol;
use rdoxide::graphmol::rgroup_decomposition::{
    RGroupDecomposition, RGroupDecompositionParameters, RGroupRows, rgroup_labels, rgroup_score,
};
use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn mol(smi: &str) -> ROMol {
    mol_from_smiles(smi, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn smarts(sma: &str) -> ROMol {
    mol_from_smarts(sma, &SmartsParserParams::default())
        .unwrap()
        .unwrap()
}

fn smiles_of(m: &ROMol) -> String {
    mol_to_smiles(&mut m.clone(), &SmilesWriteParams::default()).unwrap()
}

/// `label:smiles` of a row in label order, as `CHECK_RGROUP` builds it.
fn row_string(rows: &RGroupRows, i: usize) -> String {
    rows[i]
        .iter()
        .map(|(label, m)| format!("{label}:{}", smiles_of(m)))
        .collect::<Vec<_>>()
        .join(" ")
}

fn decompose(
    cores: &[&ROMol],
    smiles: &[&str],
    params: RGroupDecompositionParameters,
) -> (RGroupDecomposition, Vec<i32>) {
    let mut decomp = RGroupDecomposition::new_multi(cores, params).unwrap();
    let results = smiles
        .iter()
        .map(|s| decomp.add(&mol(s)).unwrap())
        .collect();
    assert!(decomp.process());
    (decomp, results)
}

const SYM_DATA: [&str; 5] = [
    "c1(Cl)ccccc1",
    "c1c(Cl)cccc1",
    "c1cccc(Cl)c1",
    "c1cc(Cl)ccc1",
    "c1ccc(Cl)cc1",
];

#[test]
fn symmetry_matching() {
    for score_method in [rgroup_score::FINGERPRINT_VARIANCE, rgroup_score::MATCH] {
        let core = mol("c1ccccc1");
        let params = RGroupDecompositionParameters {
            score_method,
            ..Default::default()
        };
        let (decomp, results) = decompose(&[&core], &SYM_DATA, params);
        assert_eq!(results, [0, 1, 2, 3, 4]);
        let rows = decomp.rgroups_as_rows().unwrap();
        for i in 0..5 {
            assert_eq!(row_string(&rows, i), "Core:c1ccc([*:1])cc1 R1:Cl[*:1]");
        }
    }
}

#[test]
fn rgroup_only_matching() {
    let core = mol("c1ccccc1[1*]");
    let params = RGroupDecompositionParameters {
        labels: rgroup_labels::ISOTOPE_LABELS,
        only_match_at_rgroups: true,
        ..Default::default()
    };
    let smiles = [
        "c1(Cl)ccccc1",
        "c1c(Cl)cccc1",
        "c1cc(Cl)ccc1",
        "c1ccc(Cl)cc1",
        "c1c(Cl)cccc(I)1",
    ];
    let (decomp, results) = decompose(&[&core], &smiles, params);
    assert_eq!(results, [0, 1, 2, 3, -1]);
    let rows = decomp.rgroups_as_rows().unwrap();
    assert_eq!(rows.len(), 4);
    for i in 0..4 {
        assert_eq!(row_string(&rows, i), "Core:c1ccc([*:1])cc1 R1:Cl[*:1]");
    }
}

#[test]
fn ring_matching() {
    let core = mol("c1ccc[1*]1");
    let mut params = RGroupDecompositionParameters {
        labels: rgroup_labels::ISOTOPE_LABELS,
        ..Default::default()
    };
    assert!(RGroupDecomposition::new(&core, params.clone()).is_err());

    params.allow_non_terminal_rgroups = true;
    let (decomp, results) = decompose(&[&core], &["c1cocc1", "c1c[nH]cc1", "c1cscc1"], params);
    assert_eq!(results, [0, 1, 2]);
    let rows = decomp.rgroups_as_rows().unwrap();
    for (i, expected) in ["Core:c1ccoc1", "Core:c1cc[nH]c1", "Core:c1ccsc1"]
        .iter()
        .enumerate()
    {
        assert_eq!(row_string(&rows, i), *expected);
    }
}

#[test]
fn ring_matching_2() {
    let core = smarts("*1***[*:1]1C[*:2]");
    let params = RGroupDecompositionParameters {
        allow_non_terminal_rgroups: true,
        ..Default::default()
    };
    let (decomp, results) = decompose(
        &[&core],
        &["c1cocc1CCl", "c1c[nH]cc1CI", "c1cscc1CF"],
        params,
    );
    assert_eq!(results, [0, 1, 2]);
    let rows = decomp.rgroups_as_rows().unwrap();
    for (i, expected) in [
        "Core:c1cc(C[*:2])co1 R2:Cl[*:2]",
        "Core:c1cc(C[*:2])c[nH]1 R2:I[*:2]",
        "Core:c1cc(C[*:2])cs1 R2:F[*:2]",
    ]
    .iter()
    .enumerate()
    {
        assert_eq!(row_string(&rows, i), *expected);
    }
}

#[test]
fn ring_matching_3() {
    let core = smarts("*1***[*:1]1");
    for score_method in [rgroup_score::MATCH, rgroup_score::FINGERPRINT_VARIANCE] {
        let params = RGroupDecompositionParameters {
            score_method,
            allow_non_terminal_rgroups: true,
            ..Default::default()
        };
        let (decomp, _) = decompose(
            &[&core],
            &["c1cocc1CCl", "c1c[nH]cc1CI", "c1cscc1CF"],
            params,
        );
        let rows = decomp.rgroups_as_rows().unwrap();
        for (i, expected) in [
            "Core:c1cc([*:1])co1 R1:ClC[*:1]",
            "Core:c1cc([*:1])c[nH]1 R1:IC[*:1]",
            "Core:c1cc([*:1])cs1 R1:FC[*:1]",
        ]
        .iter()
        .enumerate()
        {
            assert_eq!(row_string(&rows, i), *expected);
        }
    }
}

#[test]
fn multi_core() {
    let cores = [smarts("C1CCNCCC1"), mol("C1CCOCCC1"), mol("C1CCSCCC1")];
    let core_refs: Vec<&ROMol> = cores.iter().collect();
    let smiles = [
        "C1CCNC(Cl)CC1",
        "C1CC(Cl)NCCC1",
        "C1CCNC(I)CC1",
        "C1CC(I)NCCC1",
        "C1CCSC(Cl)CC1",
        "C1CC(Cl)SCCC1",
        "C1CCSC(I)CC1",
        "C1CC(I)SCCC1",
        "C1CCOC(Cl)CC1",
        "C1CC(Cl)OCCC1",
        "C1CCOC(I)CC1",
        "C1CC(I)OCCC1",
    ];
    let expected = [
        "Core:C1CCNC([*:1])CC1 R1:Cl[*:1]",
        "Core:C1CCNC([*:1])CC1 R1:Cl[*:1]",
        "Core:C1CCNC([*:1])CC1 R1:I[*:1]",
        "Core:C1CCNC([*:1])CC1 R1:I[*:1]",
        "Core:C1CCSC([*:1])CC1 R1:Cl[*:1]",
        "Core:C1CCSC([*:1])CC1 R1:Cl[*:1]",
        "Core:C1CCSC([*:1])CC1 R1:I[*:1]",
        "Core:C1CCSC([*:1])CC1 R1:I[*:1]",
        "Core:C1CCOC([*:1])CC1 R1:Cl[*:1]",
        "Core:C1CCOC([*:1])CC1 R1:Cl[*:1]",
        "Core:C1CCOC([*:1])CC1 R1:I[*:1]",
        "Core:C1CCOC([*:1])CC1 R1:I[*:1]",
    ];
    let (decomp, results) = decompose(&core_refs, &smiles, Default::default());
    assert_eq!(results, (0..12).collect::<Vec<i32>>());
    let rows = decomp.rgroups_as_rows().unwrap();
    for (i, e) in expected.iter().enumerate() {
        assert_eq!(row_string(&rows, i), *e);
    }
}

#[test]
fn ga_matching_is_reported_unsupported() {
    let core = mol("c1ccccc1");
    let params = RGroupDecompositionParameters {
        matching_strategy: rdoxide::graphmol::rgroup_decomposition::rgroup_matching::GA,
        ..Default::default()
    };
    let mut decomp = RGroupDecomposition::new(&core, params).unwrap();
    assert_eq!(decomp.add(&mol(SYM_DATA[0])).unwrap(), 0);
    assert!(!decomp.process());
}
