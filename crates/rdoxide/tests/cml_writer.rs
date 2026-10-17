//! CML writing.

use rdoxide::geometry::Point3D;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::atom::{Atom, ChiralType};
use rdoxide::graphmol::bond::BondType;
use rdoxide::graphmol::conformer::Conformer;
use rdoxide::graphmol::file_parsers::cml_writer::mol_to_cml_block;
use rdoxide::rdgeneral::RDProps;

/// S-lactic acid, as built by upstream's `testCMLWriter`.
fn lactic_acid() -> RWMol {
    let coords = [
        [-0.95330, 0.60416, 1.01609],
        [-1.00832, 1.68746, 0.83520],
        [-1.96274, 0.16103, 0.94471],
        [-0.57701, 0.44737, 2.04167],
        [0.0, 0.0, 0.0],
        [-0.43038, 0.18596, -1.01377],
        [0.22538, -1.36531, 0.19373],
        [1.21993, -1.33937, 0.14580],
        [1.38490, 0.73003, 0.0],
        [1.38490, 1.96795, 0.0],
        [2.35253, -0.07700, 0.0],
    ];
    let mut mol = RWMol::new();
    for z in [6, 1, 1, 1, 6, 1, 8, 1, 6, 8, 8] {
        mol.add_atom(Atom::new(z), true);
    }
    let bonds = [
        (0, 1, BondType::Single),
        (0, 2, BondType::Single),
        (0, 3, BondType::Single),
        (0, 4, BondType::Single),
        (4, 5, BondType::Single),
        (4, 6, BondType::Single),
        (4, 8, BondType::Single),
        (6, 7, BondType::Single),
        (8, 9, BondType::Double),
        (8, 10, BondType::Single),
    ];
    for (a, b, t) in bonds {
        mol.add_bond(a, b, t).unwrap();
    }
    let mut conf = Conformer::new(11);
    for (i, c) in coords.iter().enumerate() {
        conf.set_atom_pos(i as u32, Point3D::new(c[0], c[1], c[2]))
            .unwrap();
    }
    mol.add_conformer(conf, true).unwrap();
    mol.set_prop("_Name", "S-lactic acid".to_string(), false)
        .unwrap();
    mol.atom_mut(7).set_isotope(2);
    mol.atom_mut(10).set_formal_charge(-1);
    mol.atom_mut(4).set_chiral_tag(ChiralType::TetrahedralCCW);
    mol
}

const EXPECTED: &str = r#"<?xml version="1.0" encoding="utf-8"?>
<cml xmlns="http://www.xml-cml.org/schema" xmlns:convention="http://www.xml-cml.org/convention/" convention="convention:molecular">
  <molecule id="m-1" formalCharge="-1" spinMultiplicity="1">
    <name>S-lactic acid</name>
    <atomArray>
      <atom id="a0" elementType="C" formalCharge="0" hydrogenCount="3" x3="-0.953300" y3="0.604160" z3="1.016090"/>
      <atom id="a1" elementType="H" formalCharge="0" hydrogenCount="0" x3="-1.008320" y3="1.687460" z3="0.835200"/>
      <atom id="a2" elementType="H" formalCharge="0" hydrogenCount="0" x3="-1.962740" y3="0.161030" z3="0.944710"/>
      <atom id="a3" elementType="H" formalCharge="0" hydrogenCount="0" x3="-0.577010" y3="0.447370" z3="2.041670"/>
      <atom id="a4" elementType="C" formalCharge="0" hydrogenCount="1" x3="0.000000" y3="0.000000" z3="0.000000">
        <atomParity atomRefs4="a0 a5 a6 a8">1</atomParity>
      </atom>
      <atom id="a5" elementType="H" formalCharge="0" hydrogenCount="0" x3="-0.430380" y3="0.185960" z3="-1.013770"/>
      <atom id="a6" elementType="O" formalCharge="0" hydrogenCount="1" x3="0.225380" y3="-1.365310" z3="0.193730"/>
      <atom id="a7" elementType="H" formalCharge="0" hydrogenCount="0" isotopeNumber="2" x3="1.219930" y3="-1.339370" z3="0.145800"/>
      <atom id="a8" elementType="C" formalCharge="0" hydrogenCount="0" x3="1.384900" y3="0.730030" z3="0.000000"/>
      <atom id="a9" elementType="O" formalCharge="0" hydrogenCount="0" x3="1.384900" y3="1.967950" z3="0.000000"/>
      <atom id="a10" elementType="O" formalCharge="-1" hydrogenCount="0" x3="2.352530" y3="-0.077000" z3="0.000000"/>
    </atomArray>
    <bondArray>
      <bond atomRefs2="a0 a1" id="b0" order="S"/>
      <bond atomRefs2="a0 a2" id="b1" order="S"/>
      <bond atomRefs2="a0 a3" id="b2" order="S"/>
      <bond atomRefs2="a0 a4" id="b3" order="S"/>
      <bond atomRefs2="a4 a5" id="b4" order="S" bondStereo="H"/>
      <bond atomRefs2="a4 a6" id="b5" order="S"/>
      <bond atomRefs2="a4 a8" id="b6" order="S"/>
      <bond atomRefs2="a6 a7" id="b7" order="S"/>
      <bond atomRefs2="a8 a9" id="b8" order="D"/>
      <bond atomRefs2="a8 a10" id="b9" order="S"/>
    </bondArray>
  </molecule>
</cml>
"#;

#[test]
fn lactic_acid_matches_upstream() {
    let mol = lactic_acid();
    assert_eq!(mol_to_cml_block(&mol, -1, true).unwrap(), EXPECTED);
}

#[test]
fn escapes_names_and_handles_empty_molecules() {
    let mut mol = RWMol::new();
    mol.set_prop("_Name", "a<b & \"c\"".to_string(), false)
        .unwrap();
    let out = mol_to_cml_block(&mol, -1, true).unwrap();
    assert!(out.contains("<name>a&lt;b &amp; &quot;c&quot;</name>"));
    assert!(out.contains("<atomArray/>"));
    assert!(out.contains("<bondArray/>"));
}
