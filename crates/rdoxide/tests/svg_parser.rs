//! Reading molecules from RDKit SVG metadata.

use rdoxide::graphmol::file_parsers::svg_parser::rdkit_svg_to_mol;
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};

const SVG: &str = r#"<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' xmlns='http://www.w3.org/2000/svg' xmlns:rdkit='http://www.rdkit.org/xml' width='200px' height='200px'>
<metadata>
<rdkit:mol xmlns:rdkit = "http://www.rdkit.org/xml" version="0.9">
<rdkit:atom idx="1" atom-smiles="[CH3]" x="-2.78651" y="0.295614" z="0" />
<rdkit:atom idx="2" atom-smiles="[NH]" x="-1.35482" y="0.743114" z="0" />
<rdkit:atom idx="3" atom-smiles="[C@H]" x="-0.251428" y="-0.273019" z="0" />
<rdkit:atom idx="4" atom-smiles="[Cl]" x="-0.579728" y="-1.73665" z="0" />
<rdkit:atom idx="5" atom-smiles="[C]" x="1.18027" y="0.174481" z="0" />
<rdkit:atom idx="6" atom-smiles="[O]" x="1.50857" y="1.63811" z="0" />
<rdkit:atom idx="7" atom-smiles="[OH]" x="2.28366" y="-0.841652" z="0" />
<rdkit:bond idx="1" begin-atom-idx="1" end-atom-idx="2" bond-smiles="-" />
<rdkit:bond idx="2" begin-atom-idx="2" end-atom-idx="3" bond-smiles="-" />
<rdkit:bond idx="3" begin-atom-idx="3" end-atom-idx="4" bond-smiles="-" />
<rdkit:bond idx="4" begin-atom-idx="3" end-atom-idx="5" bond-smiles="-" />
<rdkit:bond idx="5" begin-atom-idx="5" end-atom-idx="6" bond-smiles="=" />
<rdkit:bond idx="6" begin-atom-idx="5" end-atom-idx="7" bond-smiles="-" />
</rdkit:mol></metadata>
</svg>"#;

#[test]
fn reads_embedded_molecule() {
    let mut mol = rdkit_svg_to_mol(SVG, true, true).unwrap().unwrap();
    assert_eq!(mol.num_atoms(), 7);
    assert!(!mol.conformer(-1).unwrap().is_3d());
    assert_eq!(
        mol_to_smiles(&mut mol, &SmilesWriteParams::default()).unwrap(),
        "CN[C@H](Cl)C(=O)O"
    );
}

#[test]
fn svg_without_metadata_gives_none() {
    let svg = "<svg xmlns='http://www.w3.org/2000/svg'><rect/></svg>";
    assert!(rdkit_svg_to_mol(svg, true, true).unwrap().is_none());
}

#[test]
fn malformed_xml_is_an_error() {
    assert!(rdkit_svg_to_mol("bad svg", true, true).is_err());
}
