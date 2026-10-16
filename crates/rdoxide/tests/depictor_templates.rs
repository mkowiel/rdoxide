//! Ring system templates in the depictor; the cases are those of RDKit
//! 2026.09.1's `testDepictor.py`.

use std::sync::Mutex;

use rdoxide::geometry::Point3D;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::depictor::templates::{
    add_ring_system_templates, load_default_ring_system_templates, set_ring_system_templates,
};
use rdoxide::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

/// The templates are a process-wide singleton: tests that change them must
/// not overlap with the ones that use them.
static TEMPLATES: Mutex<()> = Mutex::new(());

fn lock() -> std::sync::MutexGuard<'static, ()> {
    TEMPLATES.lock().unwrap_or_else(|e| e.into_inner())
}

fn smiles(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn depict_with_templates(mol: &mut RWMol) {
    let params = Compute2DCoordParameters {
        use_ring_templates: true,
        ..Default::default()
    };
    compute_2d_coords(mol, &params).unwrap();
}

fn positions(mol: &RWMol) -> Vec<Point3D> {
    mol.conformer(-1).unwrap().positions().to_vec()
}

/// Whether two layouts are the same up to a rigid motion (all the
/// interatomic distances agree).
fn same_layout(a: &[Point3D], b: &[Point3D]) -> bool {
    a.len() == b.len()
        && (0..a.len()).all(|i| {
            (0..i).all(|j| {
                let da = ((a[i].x - a[j].x).powi(2) + (a[i].y - a[j].y).powi(2)).sqrt();
                let db = ((b[i].x - b[j].x).powi(2) + (b[i].y - b[j].y).powi(2)).sqrt();
                (da - db).abs() < 0.02
            })
        })
}

const CAGE: &str = "C1CC2CCOC3OC4CCC(C1)C23OO4";
const CAGE_DEFAULT: &str = "C1CC2CCOC3OC4CCC(C1)C23OO4 |(3.53,-1.22,;3.53,0.3,;2.21,1.06,;2.21,2.59,;0.89,3.35,;-0.43,2.59,;-0.43,1.06,;-1.9,0.65,;-2.47,-0.76,;-1.71,-2.08,;-0.2,-2.29,;0.89,-1.22,;2.21,-1.99,;0.89,0.3,;0.12,-0.83,;-1.19,-1.25,)|";
const CAGE_USER: &str = "C1CC2CCOC3OC4CCC(C1)C23OO4 |(-0.5537,-3.1595,;-1.6057,-2.003,;-1.4262,-0.4072,;-2.9804,0.0271,;-3.5191,1.502,;-2.2028,2.3562,;-0.6818,1.8511,;1.0592,1.4391,;2.6123,1.8366,;3.5191,0.5341,;2.6067,-0.7521,;1.0061,-0.773,;0.7888,-2.3546,;-0.0405,0.5251,;0.4049,2.3,;1.7604,3.1594,)|";

#[test]
fn default_templates_are_used_on_request() {
    let _guard = lock();
    load_default_ring_system_templates();
    let template = smiles(CAGE_DEFAULT);
    let mut with = smiles(CAGE);
    depict_with_templates(&mut with);
    assert!(same_layout(&positions(&with), &positions(&template)));
    let mut without = smiles(CAGE);
    compute_2d_coords(&mut without, &Compute2DCoordParameters::default()).unwrap();
    assert!(!same_layout(&positions(&without), &positions(&template)));
}

#[test]
fn templates_can_be_replaced_and_restored() {
    let _guard = lock();
    let dir = std::env::temp_dir().join(format!("rdoxide_templates_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let user = dir.join("user.smi");
    std::fs::write(&user, format!("{CAGE_USER}\n")).unwrap();

    set_ring_system_templates(user.to_str().unwrap()).unwrap();
    let mut mol = smiles(CAGE);
    depict_with_templates(&mut mol);
    assert!(same_layout(
        &positions(&mol),
        &positions(&smiles(CAGE_USER))
    ));

    load_default_ring_system_templates();
    let mut mol = smiles(CAGE);
    depict_with_templates(&mut mol);
    assert!(same_layout(
        &positions(&mol),
        &positions(&smiles(CAGE_DEFAULT))
    ));

    // a template added later takes precedence over the existing ones
    add_ring_system_templates(user.to_str().unwrap()).unwrap();
    let mut mol = smiles(CAGE);
    depict_with_templates(&mut mol);
    assert!(same_layout(
        &positions(&mol),
        &positions(&smiles(CAGE_USER))
    ));
    load_default_ring_system_templates();
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn invalid_templates_are_errors() {
    let _guard = lock();
    let dir = std::env::temp_dir().join(format!("rdoxide_bad_templates_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for (i, content) in [
        "invalidsmiles",
        "C |(-0.5537,-3.1595)|",
        "C1CCCCC1",
        "c1ccccc1-c1ccccc1",
    ]
    .into_iter()
    .enumerate()
    {
        let path = dir.join(format!("bad{i}.smi"));
        std::fs::write(&path, content).unwrap();
        let err = set_ring_system_templates(path.to_str().unwrap()).unwrap_err();
        assert!(err.to_string().contains("Depict error"), "{content}: {err}");
    }
    assert!(set_ring_system_templates("/does/not/exist.smi").is_err());
    // the templates in use were not touched
    let mut mol = smiles(CAGE);
    depict_with_templates(&mut mol);
    assert!(same_layout(
        &positions(&mol),
        &positions(&smiles(CAGE_DEFAULT))
    ));
    std::fs::remove_dir_all(&dir).unwrap();
}
