//! Gaussian shape overlay; the cases of upstream
//! `GaussianShape/Wrap/test_rdgaussian_shape.py`.

use rdoxide::geometry::Point3D;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::file_parsers::mol_file_parser::MolFileParserParams;
use rdoxide::graphmol::file_parsers::sd_mol_supplier::SDMolSupplier;
use rdoxide::graphmol::gaussian_shape::{
    CustomFeature, OptimMode, ShapeInput, ShapeInputOptions, ShapeOverlayOptions, align_molecule,
    score_molecule, score_molecule_all_conformers, score_shape,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

const DATA: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../upstream/rdkit/Code/GraphMol/GaussianShape/test_data/test1.sdf"
);

fn mols() -> (RWMol, RWMol) {
    let mut s = SDMolSupplier::from_file(DATA, MolFileParserParams::default()).unwrap();
    (s.get(0).unwrap().unwrap(), s.get(1).unwrap().unwrap())
}

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{a} != {b}");
}

#[test]
fn defaults() {
    let (r, mut p) = mols();
    let (s, _) = align_molecule(
        &r,
        &mut p,
        &ShapeInputOptions::default(),
        &ShapeInputOptions::default(),
        &ShapeOverlayOptions::default(),
        -1,
        -1,
    )
    .unwrap();
    close(s[0], 0.497, 5e-4);
    close(s[1], 0.760, 5e-4);
    close(s[2], 0.233, 5e-4);
}

#[test]
fn no_color() {
    let (r, mut p) = mols();
    let ov = ShapeOverlayOptions {
        optim_mode: OptimMode::ShapeOnly,
        ..Default::default()
    };
    let so = ShapeInputOptions {
        use_colors: false,
        ..Default::default()
    };
    let (s, _) = align_molecule(&r, &mut p, &so, &so, &ov, -1, -1).unwrap();
    close(s[0], 0.760, 5e-4);
    close(s[1], 0.760, 5e-4);
    close(s[2], 0.0, 5e-4);
}

#[test]
fn from_shape() {
    let (r, _) = mols();
    let shp = ShapeInput::new(
        &r,
        -1,
        &ShapeInputOptions::default(),
        &ShapeOverlayOptions::default(),
    )
    .unwrap();
    close(shp.shape_volume(), 591.058, 5e-4);
    close(shp.color_volume(), 31.935, 5e-4);
    assert_eq!(shp.shape_to_mol(false, false).unwrap().num_atoms(), 13);
    assert_eq!(shp.shape_to_mol(true, false).unwrap().num_atoms(), 19);
}

#[test]
fn fixed_score() {
    let (r, p) = mols();
    let o = ShapeInputOptions::default();
    let ov = ShapeOverlayOptions::default();
    let (s, _) = score_molecule(&r, &r, &o, &o, &ov, -1, -1).unwrap();
    close(s[0], 1.0, 5e-4);
    let shp1 = ShapeInput::new(&p, -1, &o, &ov).unwrap();
    let shp2 = ShapeInput::new(&p, -1, &o, &ov).unwrap();
    let (s, vols) = score_shape(&shp1, &shp2, &ov);
    close(s[0], 1.0, 5e-4);
    close(vols.0, 751.0, 0.05);
    close(vols.1, 42.5, 0.05);
}

#[test]
fn custom_atom_radii_and_subset() {
    let (r, _) = mols();
    let ov = ShapeOverlayOptions::default();
    let opts = ShapeInputOptions {
        all_carbon_radii: false,
        atom_radii: vec![(4, 1.9)],
        ..Default::default()
    };
    let shp = ShapeInput::new(&r, -1, &opts, &ov).unwrap();
    close(shp.shape_volume(), 559.361, 5e-4);

    let opts = ShapeInputOptions {
        atom_subset: vec![4, 5, 6, 7, 8, 9],
        use_colors: false,
        all_carbon_radii: false,
        ..Default::default()
    };
    let shp = ShapeInput::new(&r, -1, &opts, &ov).unwrap();
    close(shp.shape_volume(), 259.144, 5e-4);
    assert_eq!(shp.color_volume(), 0.0);
}

#[test]
fn custom_features_are_used() {
    let (r, _) = mols();
    let opts = ShapeInputOptions {
        custom_features: vec![vec![CustomFeature {
            kind: 1,
            pos: Point3D::new(0.0, 0.0, 0.0),
            rad: 1.0,
            atoms: vec![],
        }]],
        ..Default::default()
    };
    let shp = ShapeInput::new(&r, -1, &opts, &ShapeOverlayOptions::default()).unwrap();
    assert_eq!(shp.num_features(), 1);
}

#[test]
fn tversky() {
    let (r, mut p) = mols();
    let o = ShapeInputOptions::default();
    let mut ov = ShapeOverlayOptions {
        sim_alpha: 0.95,
        sim_beta: 0.05,
        ..Default::default()
    };
    let (s, _) = align_molecule(&r, &mut p, &o, &o, &ov, -1, -1).unwrap();
    close(s[0], 0.700, 5e-4);
    close(s[1], 0.968, 5e-4);
    close(s[2], 0.433, 5e-4);
    ov.sim_alpha = 0.05;
    ov.sim_beta = 0.95;
    let (s, _) = align_molecule(&r, &mut p, &o, &o, &ov, -1, -1).unwrap();
    close(s[0], 0.557, 5e-4);
    close(s[1], 0.780, 5e-4);
    close(s[2], 0.335, 5e-4);
}

fn load_multi(file: &str) -> RWMol {
    let path = DATA.replace("test1.sdf", file);
    let text = std::fs::read_to_string(path).unwrap();
    let mut mol: Option<RWMol> = None;
    for line in text.lines().map(str::trim).filter(|l| !l.is_empty()) {
        let m = mol_from_smiles(line, &SmilesParserParams::default())
            .unwrap()
            .unwrap();
        match mol.as_mut() {
            None => mol = Some(m),
            Some(mol) => {
                mol.add_conformer(m.conformer(-1).unwrap().clone(), true)
                    .unwrap();
            }
        }
    }
    mol.unwrap()
}

#[test]
fn multiple_conformers() {
    let eso = load_multi("esomeprazole_multi.smi");
    let ran = load_multi("ranitidine_multi.smi");
    assert_eq!(eso.num_conformers(), 10);
    let ov = ShapeOverlayOptions::default();
    let o = ShapeInputOptions::default();
    let mut shapes1 = ShapeInput::new(&eso, -1, &o, &ov).unwrap();
    let shapes2 = ShapeInput::new(&ran, -1, &o, &ov).unwrap();
    assert_eq!(shapes1.num_shapes(), 10);
    assert_eq!(shapes2.num_shapes(), 10);
    let (sim, pair) = shapes1.best_similarity(&shapes2, -1.0, &ov).unwrap();
    close(sim[0], 0.449, 5e-4);
    let pair = pair.unwrap();
    assert_eq!((pair.this_shape, pair.fit_shape), (2, 6));

    let all = score_molecule_all_conformers(&eso, &ran, &o, &o, &ov).unwrap();
    assert_eq!((all.ref_conf_id, all.fit_conf_id), (8, 3));
    close(all.comb_scores[8][3], 0.449, 5e-4);
}
