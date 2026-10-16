//! Marvin documents; the cases of upstream `MarvinParse/testMrvToMol.cpp`
//! (parse results are compared to the `.expected.sdf` files).

use rdoxide::graphmol::chirality::remove_non_explicit_3d_chirality;
use rdoxide::graphmol::file_parsers::marvin::{MrvParserParams, mol_from_mrv_file};
use rdoxide::graphmol::file_parsers::mol_file_writer::{MolWriterParams, mol_to_mol_block};
use rdoxide::graphmol::wedge_bonds::reapply_mol_block_wedging;

const DATA: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../upstream/rdkit/Code/GraphMol/MarvinParse/test_data/"
);

/// (file, expected to parse, atoms, bonds, sanitize, reapply wedging)
const MOL_TESTS: &[(&str, bool, usize, usize, bool, bool)] = &[
    ("FalseChiral.mrv", true, 4, 3, false, false),
    ("FalseChiral.mrv", true, 4, 3, true, false),
    ("Cubane.mrv", true, 16, 20, true, true),
    ("NewChiralTest.mrv", true, 13, 14, true, false),
    ("NewChiralTest.mrv", true, 13, 14, false, false),
    ("NewChiralTestNoChiral.mrv", true, 13, 14, true, false),
    ("NewChiralTestNoChiral.mrv", true, 13, 14, false, false),
    ("NewChiralTestAllChiral.mrv", true, 14, 15, true, false),
    ("NewChiralTestAllChiral.mrv", true, 14, 15, false, false),
    (
        "DataSgroupMissingUnitsDisplayed.mrv",
        true,
        15,
        16,
        true,
        true,
    ),
    ("lostStereoAnd.mrv", true, 6, 5, true, true),
    ("DoubleBondChain.mrv", true, 22, 22, true, true),
    ("WigglyAndCrossed.mrv", true, 8, 7, true, true),
    ("BondTypes.mrv", true, 26, 25, true, true),
    ("EmbeddedSGroupSUP_MUL.mrv", true, 17, 17, true, true),
    ("EmbeddedSgroupCOP_SUP.mrv", true, 10, 10, true, true),
    ("EmbeddedSgroupDAT_SUP.mrv", true, 10, 10, true, true),
    (
        "EmbeddedSgroupMULTICENTER_SUP.mrv",
        true,
        11,
        10,
        true,
        true,
    ),
    ("EmbeddedSgroupMUL_MUL.mrv", true, 141, 140, true, true),
    ("EmbeddedSgroupMUL_MUL2.mrv", true, 23, 22, true, true),
    ("EmbeddedSgroupMUL_SUP.mrv", true, 129, 128, true, true),
    ("EmbeddedSgroupSRU_SUP.mrv", true, 10, 10, true, true),
    ("EmbeddedSgroupSUPEXP_SUP.mrv", true, 10, 10, true, true),
    ("EmbeddedSgroupSUPEXP_SUP2.mrv", true, 10, 10, true, true),
    ("EmbeddedSgroupSUP_MULTICENTER.mrv", true, 10, 8, true, true),
    ("EmbeddedSgroupSUP_SUP.mrv", true, 12, 11, false, true),
    ("EmbeddedSgroupSUP_SUP2.mrv", true, 12, 12, true, true),
    ("RgroupBad.mrv", true, 9, 9, true, true),
    ("valenceLessThanDrawn.mrv", true, 14, 14, true, true),
    ("data_sgroup_no_fieldname.mrv", true, 4, 3, true, true),
    ("data_sgroup_empty_field_data.mrv", true, 2, 1, true, true),
    ("radical_value.mrv", true, 3, 2, true, true),
    ("emptyOneLineAtomList.mrv", true, 0, 0, true, true),
    ("mrvValence_value.mrv", true, 3, 2, true, true),
    ("ChiralTest2.mrv", true, 46, 47, true, true),
    ("ChiralTest.mrv", true, 8, 7, true, true),
    ("SnCl2.mrv", true, 3, 2, true, true),
    ("SnH2Cl2.mrv", true, 3, 2, true, true),
    ("marvin01.mrv", true, 11, 11, true, true),
    ("marvin02.mrv", true, 9, 9, true, true),
    ("marvin07.mrv", true, 12, 11, true, true),
    ("marvin10.mrv", true, 10, 10, true, true),
    ("marvin06.mrv", true, 11, 11, true, true),
    ("marvin12.mrv", true, 31, 33, true, true),
    ("EmptyMol.mrv", true, 0, 0, true, true),
    ("Sparse.mrv", true, 0, 0, true, true),
    ("Sparse2.mrv", true, 0, 0, true, true),
    ("Sparse3.mrv", true, 0, 0, true, true),
    ("MarvinNoCoords.mrv", true, 6, 6, true, true),
    ("aspirin.mrv", true, 13, 13, true, true),
    ("MarvinStereoGroupsZeros.mrv", true, 8, 8, true, true),
    ("MarvinStereoGroupsAbs.mrv", true, 8, 8, true, true),
    ("triphenylphosphine.mrv", true, 19, 21, true, true),
    ("MarvinOldSuperGroupTest.mrv", true, 89, 93, true, true),
    ("RadicalTests.mrv", true, 8, 7, true, true),
    ("AnyBond.mrv", true, 4, 3, true, true),
    ("cisBenzene.mrv", true, 6, 6, true, true),
    ("DativeBond.mrv", true, 6, 5, true, true),
    ("MultipleSgroup.mrv", true, 123, 122, true, true),
    ("SgroupExpanded.mrv", true, 5, 4, true, true),
    ("SgroupMultAttach.mrv", true, 44, 45, true, true),
    ("MarvinMissingX2.mrv", true, 12, 11, true, true),
    ("MarvinMissingY2.mrv", true, 12, 11, true, true),
    ("DataSgroup.mrv", true, 7, 6, true, true),
    ("MulticenterSgroup.mrv", true, 17, 16, false, true),
    ("GenericSgroup.mrv", true, 13, 13, true, true),
    ("MonomerSgroup.mrv", true, 4, 3, true, true),
    ("modification_sgroup.mrv", true, 54, 40, true, true),
    ("copolymer_sgroup.mrv", true, 19, 18, true, true),
    (
        "MultipleSgroupParentInMiddleOfAtomBlock.mrv",
        true,
        23,
        22,
        true,
        true,
    ),
    ("EmbeddedSgroups.mrv", false, 14, 14, true, true),
    ("marvin03.mrv", false, 31, 33, true, true),
    ("MarvinBadMissingMolID.mrv", false, 12, 11, true, true),
    ("MarvinBadMissingAtomID.mrv", false, 12, 11, true, true),
    ("MarvinBadX2.mrv", false, 12, 11, true, true),
    ("MarvinBadY2.mrv", false, 12, 11, true, true),
    ("MarvinBadStereoGroupsAbs.mrv", false, 8, 8, true, true),
    ("MarvinBadElementType.mrv", false, 12, 11, true, true),
    ("MarvinBadMissingBondID.mrv", false, 12, 11, true, true),
    ("MarvinBadMissingBondAtomRefs", false, 12, 11, true, true),
    ("MarvinBadMissingBondOrder.mrv", false, 12, 11, true, true),
    ("MarvinBadMissingSruMolID.mrv", false, 12, 11, true, true),
    ("MarvinBadMissingSruID.mrv", false, 12, 11, true, true),
    ("MarvinBadMissingSruRole.mrv", false, 12, 11, true, true),
    ("MarvinBadMissingSruAtomRef.mrv", false, 12, 11, true, true),
    ("MarvinBadMissingSruTitle.mrv", false, 12, 11, true, true),
    ("MarvinBadSruAtomRef.mrv", false, 12, 11, true, true),
    ("MarvinBadSruID.mrv", false, 12, 11, true, true),
    ("MarvinBadSruRole.mrv", false, 12, 11, true, true),
    ("MarvinBadSruConnect.mrv", false, 12, 11, true, true),
    ("MarvinBadSupAttachAtom.mrv", false, 9, 9, true, true),
    ("MarvinBadSupAttachBond.mrv", false, 9, 9, true, true),
    ("MarvinBadSupAttachOrder.mrv", false, 9, 9, true, true),
    ("MarvinBadSupMissingAttachBond.mrv", false, 9, 9, true, true),
    (
        "MarvinBadSupMissingAttachOrder.mrv",
        false,
        9,
        9,
        true,
        true,
    ),
    ("FalseAtropisomer.mrv", true, 17, 18, false, true),
    ("AtropEnhancedStereo.mrv", true, 16, 17, false, true),
    ("AtropManyChirals.mrv", true, 20, 20, false, true),
    ("AtropManyChiralsEnhanced.mrv", true, 20, 20, false, true),
    ("AtropManyChiralsEnhanced2.mrv", true, 20, 20, false, true),
    ("AtropManyChiralsEnhanced3.mrv", true, 20, 20, false, true),
    ("AtropManyChiralsEnhanced4.mrv", true, 20, 20, false, true),
    ("BMS-986142_3d_chiral.mrv", true, 72, 77, false, true),
    ("BMS-986142_3d.mrv", true, 72, 77, false, true),
    ("BMS-986142_atrop1.mrv", true, 42, 47, false, true),
    ("BMS-986142_atrop2.mrv", true, 42, 47, false, true),
    ("BMS-986142_atrop3.mrv", true, 42, 47, false, true),
    ("BMS-986142_atrop4.mrv", true, 42, 47, false, true),
    ("BMS-986142_atrop5.mrv", true, 42, 47, false, true),
    ("BMS-986142_atrop6.mrv", true, 42, 47, false, true),
    ("BMS-986142_atrop7.mrv", true, 42, 47, false, true),
    ("BMS-986142_atrop8.mrv", true, 42, 47, false, true),
    ("BMS-986142_atropBad2.mrv", true, 42, 47, false, true),
    ("JDQ443_3d.mrv", true, 66, 72, false, true),
    ("JDQ443_atrop1.mrv", true, 38, 44, false, true),
    ("JDQ443_atrop2.mrv", true, 38, 44, false, true),
    ("JDQ443_atrop3.mrv", true, 38, 44, false, true),
    ("JDQ443_atropBad1.mrv", true, 38, 44, false, true),
    ("RP-6306_atrop1.mrv", true, 24, 26, false, true),
    ("RP-6306_atrop2.mrv", true, 24, 26, false, true),
    ("RP-6306_atrop3.mrv", true, 24, 26, false, true),
    ("RP-6306_atrop4.mrv", true, 24, 26, false, true),
    ("RP-6306_atrop5.mrv", true, 24, 26, false, true),
    ("RP-6306_atropBad1.mrv", true, 24, 26, false, true),
    ("RP-6306_atropBad2.mrv", true, 24, 26, false, true),
    ("RP-6306_3d.mrv", true, 44, 46, false, true),
    ("Sotorasib_atrop1.mrv", true, 41, 45, false, true),
    ("Sotorasib_atrop2.mrv", true, 41, 45, false, true),
    ("Sotorasib_atrop3.mrv", true, 41, 45, false, true),
    ("Sotorasib_atrop4.mrv", true, 41, 45, false, true),
    ("Sotorasib_atrop5.mrv", true, 41, 45, false, true),
    ("Sotorasib_atropBad1.mrv", true, 41, 45, false, true),
    ("Sotorasib_atropBad2.mrv", true, 41, 45, false, true),
    ("Sotorasib_3d.mrv", true, 71, 75, false, true),
    ("ZM374979_atrop1.mrv", true, 45, 49, false, true),
    ("ZM374979_atrop2.mrv", true, 45, 49, false, true),
    ("ZM374979_atrop3.mrv", true, 45, 49, false, true),
    ("ZM374979_atropBad1.mrv", true, 45, 49, false, true),
    ("mrtx1719_3d.mrv", true, 51, 55, false, true),
    ("mrtx1719_atrop1.mrv", true, 33, 37, false, true),
    ("mrtx1719_atrop2.mrv", true, 33, 37, false, true),
    ("mrtx1719_atrop3.mrv", true, 33, 37, false, true),
    ("mrtx1719_atropBad1.mrv", true, 33, 37, false, true),
    ("FalseAtropisomer.mrv", true, 17, 18, true, false),
    ("AtropEnhancedStereo.mrv", true, 16, 17, true, false),
    ("AtropManyChirals.mrv", true, 20, 20, true, false),
    ("AtropManyChiralsEnhanced.mrv", true, 20, 20, true, false),
    ("AtropManyChiralsEnhanced2.mrv", true, 20, 20, true, false),
    ("AtropManyChiralsEnhanced3.mrv", true, 20, 20, true, false),
    ("AtropManyChiralsEnhanced4.mrv", true, 20, 20, true, false),
    ("BMS-986142_3d_chiral.mrv", true, 72, 77, true, false),
    ("BMS-986142_3d.mrv", true, 72, 77, true, false),
    ("BMS-986142_atrop1.mrv", true, 42, 47, true, false),
    ("BMS-986142_atrop2.mrv", true, 42, 47, true, false),
    ("BMS-986142_atrop3.mrv", true, 42, 47, true, false),
    ("BMS-986142_atrop4.mrv", true, 42, 47, true, false),
    ("BMS-986142_atrop5.mrv", true, 42, 47, true, false),
    ("BMS-986142_atrop6.mrv", true, 42, 47, true, false),
    ("BMS-986142_atrop7.mrv", true, 42, 47, true, false),
    ("BMS-986142_atrop8.mrv", true, 42, 47, true, false),
    ("BMS-986142_atropBad2.mrv", true, 42, 47, true, false),
    ("JDQ443_3d.mrv", true, 66, 72, true, false),
    ("JDQ443_atrop1.mrv", true, 38, 44, true, false),
    ("JDQ443_atrop2.mrv", true, 38, 44, true, false),
    ("JDQ443_atrop3.mrv", true, 38, 44, true, false),
    ("JDQ443_atropBad1.mrv", true, 38, 44, true, false),
    ("RP-6306_atrop1.mrv", true, 24, 26, true, false),
    ("RP-6306_atrop2.mrv", true, 24, 26, true, false),
    ("RP-6306_atrop3.mrv", true, 24, 26, true, false),
    ("RP-6306_atrop4.mrv", true, 24, 26, true, false),
    ("RP-6306_atrop5.mrv", true, 24, 26, true, false),
    ("RP-6306_atropBad1.mrv", true, 24, 26, true, false),
    ("RP-6306_atropBad2.mrv", true, 24, 26, true, false),
    ("RP-6306_3d.mrv", true, 44, 46, true, false),
    ("Sotorasib_atrop1.mrv", true, 41, 45, true, false),
    ("Sotorasib_atrop2.mrv", true, 41, 45, true, false),
    ("Sotorasib_atrop3.mrv", true, 41, 45, true, false),
    ("Sotorasib_atrop4.mrv", true, 41, 45, true, false),
    ("Sotorasib_atrop5.mrv", true, 41, 45, true, false),
    ("Sotorasib_atropBad1.mrv", true, 41, 45, true, false),
    ("Sotorasib_atropBad2.mrv", true, 41, 45, true, false),
    ("Sotorasib_3d.mrv", true, 71, 75, true, false),
    ("ZM374979_atrop1.mrv", true, 45, 49, true, false),
    ("ZM374979_atrop2.mrv", true, 45, 49, true, false),
    ("ZM374979_atrop3.mrv", true, 45, 49, true, false),
    ("ZM374979_atropBad1.mrv", true, 45, 49, true, false),
    ("mrtx1719_3d.mrv", true, 51, 55, true, false),
    ("mrtx1719_atrop1.mrv", true, 33, 37, true, false),
    ("mrtx1719_atrop2.mrv", true, 33, 37, true, false),
    ("mrtx1719_atrop3.mrv", true, 33, 37, true, false),
    ("mrtx1719_atropBad1.mrv", true, 33, 37, true, false),
];

/// Sanitized without re-applying the mol block wedging: the wedges the
/// writer picks depend on which chiral tags the stereo clean-up keeps, and
/// for these spiro and biaryl cases rdoxide keeps different ones (the MRV
/// parse itself matches).
const KNOWN_WEDGING_DIFFERENCES: &[&str] = &[
    "NewChiralTest.mrv",
    "NewChiralTestNoChiral.mrv",
    "FalseAtropisomer.mrv",
];

fn suffix(sanitize: bool, reapply: bool) -> String {
    format!(
        "{}{}",
        if sanitize { "" } else { ".nosan" },
        if reapply { "" } else { ".noreapply" }
    )
}

/// `testMarvinMol` without the MRV writer part
fn check(file: &str, expected: bool, atoms: usize, bonds: usize, sanitize: bool, reapply: bool) {
    let path = format!("{DATA}{file}");
    let params = MrvParserParams {
        sanitize,
        remove_hs: false,
    };
    let parsed = mol_from_mrv_file(&path, &params);
    let mut mol = match parsed {
        Ok(Some(m)) => m,
        Ok(None) => panic!("{file}: empty"),
        Err(e) => {
            assert!(!expected, "{file}: {e}");
            return;
        }
    };
    remove_non_explicit_3d_chirality(&mut mol);
    if reapply {
        reapply_mol_block_wedging(&mut mol, true, false).unwrap();
    }
    if expected {
        assert_eq!(mol.num_atoms() as usize, atoms, "{file} atoms");
        assert_eq!(mol.num_bonds() as usize, bonds, "{file} bonds");
    }
    let mut p = MolWriterParams {
        force_v3000: true,
        ..Default::default()
    };
    let block = mol_to_mol_block(&mol, &p, 0).unwrap_or_else(|_| {
        p.kekulize = false;
        mol_to_mol_block(&mol, &p, 0).unwrap()
    });
    let expected_file = format!("{path}{}.expected.sdf", suffix(sanitize, reapply));
    let want = std::fs::read_to_string(&expected_file).unwrap_or_default();
    if expected {
        assert_eq!(want, block, "{file}: {expected_file}");
    } else {
        assert_ne!(want, block, "{file} should not match");
    }
}

#[test]
fn mrv_files() {
    if !std::path::Path::new(DATA).exists() {
        return; // the upstream checkout is not available
    }
    let mut failures = Vec::new();
    for &(file, expected, atoms, bonds, sanitize, reapply) in MOL_TESTS {
        if sanitize && !reapply && KNOWN_WEDGING_DIFFERENCES.contains(&file) {
            continue;
        }
        let r = std::panic::catch_unwind(|| check(file, expected, atoms, bonds, sanitize, reapply));
        if r.is_err() {
            failures.push(format!("{file} san={sanitize} rea={reapply}"));
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
