//! Lipinski-style atom, bond and ring counts (upstream
//! `GraphMol/Descriptors/Lipinski.{h,cpp}`). `calcPhi` (ConnectivityDescriptors)
//! is not part of this module yet.

use std::sync::OnceLock;

use crate::graphmol::bond::BondType;
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::find_rings::find_sssr;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::utils::intersect;
use crate::{RDKitError, Result};

use crate::graphmol::atom::ChiralType;

pub const LIPINSKI_HBA_VERSION: &str = "1.0.0";
pub const LIPINSKI_HBD_VERSION: &str = "2.0.0";
pub const NUM_ROTATABLE_BONDS_VERSION: &str = "3.2.0";
pub const NUM_HBD_VERSION: &str = "2.0.1";
pub const NUM_HBA_VERSION: &str = "2.0.2";
pub const NUM_HETEROATOMS_VERSION: &str = "1.0.1";
pub const NUM_AMIDE_BONDS_VERSION: &str = "1.0.0";
pub const FRACTION_CSP3_VERSION: &str = "1.0.0";
pub const NUM_RINGS_VERSION: &str = "1.0.1";
pub const NUM_HETEROCYCLES_VERSION: &str = "1.0.0";
pub const NUM_AROMATIC_RINGS_VERSION: &str = "1.0.0";
pub const NUM_SATURATED_RINGS_VERSION: &str = "1.0.0";
pub const NUM_ALIPHATIC_RINGS_VERSION: &str = "1.0.0";
pub const NUM_AROMATIC_HETEROCYCLES_VERSION: &str = "1.0.0";
pub const NUM_AROMATIC_CARBOCYCLES_VERSION: &str = "1.0.0";
pub const NUM_SATURATED_HETEROCYCLES_VERSION: &str = "1.0.0";
pub const NUM_SATURATED_CARBOCYCLES_VERSION: &str = "1.0.0";
pub const NUM_ALIPHATIC_HETEROCYCLES_VERSION: &str = "1.0.0";
pub const NUM_ALIPHATIC_CARBOCYCLES_VERSION: &str = "1.0.0";
pub const NUM_SPIRO_ATOMS_VERSION: &str = "1.0.0";
pub const NUM_BRIDGEHEAD_ATOMS_VERSION: &str = "2.0.0";
pub const NUM_ATOM_STEREO_CENTERS_VERSION: &str = "1.0.1";
pub const NUM_UNSPECIFIED_ATOM_STEREO_CENTERS_VERSION: &str = "1.0.1";

/// `NumRotatableBondsOptions`
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NumRotatableBondsOptions {
    Default = -1,
    NonStrict = 0,
    Strict = 1,
    StrictLinkages = 2,
}

/// Counts the matches of a SMARTS query (compiled once per call site, like
/// upstream's `pattern_flyweight`). Matches are uniquified, as in `SubstructMatch`.
fn count_matches(mol: &ROMol, cell: &'static OnceLock<ROMol>, smarts: &str) -> Result<u32> {
    let query = match cell.get() {
        Some(q) => q,
        None => {
            let q = mol_from_smarts(smarts, &SmartsParserParams::default())?
                .ok_or_else(|| RDKitError::invariant("bad Lipinski SMARTS"))?;
            cell.get_or_init(|| q)
        }
    };
    Ok(substruct_match(mol, query, &SubstructMatchParameters::default())?.len() as u32)
}

/// `calcLipinskiHBA`: N and O atoms.
pub fn calc_lipinski_hba(mol: &ROMol) -> u32 {
    (0..mol.num_atoms())
        .filter(|&i| matches!(mol.atom(i).atomic_num(), 7 | 8))
        .count() as u32
}

/// `calcLipinskiHBD`: the hydrogens (explicit and implicit) on N and O atoms.
pub fn calc_lipinski_hbd(mol: &ROMol) -> Result<u32> {
    let mut res = 0;
    for i in 0..mol.num_atoms() {
        if matches!(mol.atom(i).atomic_num(), 7 | 8) {
            res += mol.total_num_hs(i, true)?;
        }
    }
    Ok(res)
}

/// `calcNumRotatableBonds`. `Default` is the Strict definition: upstream's
/// `RDK_USE_STRICT_ROTOR_DEFINITION` CMake option is ON by default and the
/// released builds keep it.
pub fn calc_num_rotatable_bonds(mol: &ROMol, strict: NumRotatableBondsOptions) -> Result<u32> {
    static NON_STRICT: OnceLock<ROMol> = OnceLock::new();
    static STRICT: OnceLock<ROMol> = OnceLock::new();
    static ROT_BONDS: OnceLock<ROMol> = OnceLock::new();
    static SYM_RINGS: OnceLock<ROMol> = OnceLock::new();
    static TERMINAL_TRIPLE: OnceLock<ROMol> = OnceLock::new();
    static RING_AMIDES: OnceLock<ROMol> = OnceLock::new();

    let strict = match strict {
        NumRotatableBondsOptions::Default => NumRotatableBondsOptions::Strict,
        s => s,
    };
    match strict {
        NumRotatableBondsOptions::NonStrict => {
            count_matches(mol, &NON_STRICT, "[!$(*#*)&!D1]-,:;!@[!$(*#*)&!D1]")
        }
        NumRotatableBondsOptions::Strict => count_matches(
            mol,
            &STRICT,
            concat!(
                "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])(",
                "[CH3])[CH3])&!$([CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=",
                "[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[!$",
                "(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])(",
                "[CH3])[CH3])&!$([CH3])]"
            ),
        ),
        _ => strict_linkages(mol, &ROT_BONDS, &SYM_RINGS, &TERMINAL_TRIPLE, &RING_AMIDES),
    }
}

/// The `StrictLinkages` branch of `calcNumRotatableBonds`.
fn strict_linkages(
    mol: &ROMol,
    rot_bonds: &'static OnceLock<ROMol>,
    sym_rings: &'static OnceLock<ROMol>,
    terminal_triple: &'static OnceLock<ROMol>,
    ring_amides: &'static OnceLock<ROMol>,
) -> Result<u32> {
    // initialize to the number of bonds matching the base pattern
    let mut res = count_matches(mol, rot_bonds, "[!$([D1&!#1])]-,:;!@[!$([D1&!#1])]")? as i64;
    if res == 0 {
        return Ok(0);
    }
    // remove symmetrical rings
    res -= count_matches(
        mol,
        sym_rings,
        "[a;r6;$(a(-,:;!@[a;r6])(a[!#1])a[!#1])]-,:;!@[a;r6;$(a(-,:;!@[a;r6])(a[!#1])a)]",
    )? as i64;
    res = res.max(0);
    // remove triple bonds
    res -= count_matches(mol, terminal_triple, "C#[#6,#7]")? as i64;
    res = res.max(0);
    // removing amides is more complex: an amide only counts once per distinct atom set
    let query = match ring_amides.get() {
        Some(q) => q,
        None => {
            let q = mol_from_smarts("[C&!R](=O)NC", &SmartsParserParams::default())?
                .ok_or_else(|| RDKitError::invariant("bad Lipinski SMARTS"))?;
            ring_amides.get_or_init(|| q)
        }
    };
    let mut atoms_seen = vec![false; mol.num_atoms() as usize];
    for iv in substruct_match(mol, query, &SubstructMatchParameters::default())? {
        let mut distinct = true;
        for (_, mol_atom) in &iv {
            let a = *mol_atom as usize;
            if atoms_seen[a] {
                distinct = false;
            }
            atoms_seen[a] = true;
        }
        if distinct && res > 0 {
            res -= 1;
        }
    }
    res = res.max(0);
    Ok(res as u32)
}

/// `calcNumRotatableBonds(mol, strict)` with a bool: `true` is `Strict`.
pub fn calc_num_rotatable_bonds_bool(mol: &ROMol, strict: bool) -> Result<u32> {
    calc_num_rotatable_bonds(
        mol,
        if strict {
            NumRotatableBondsOptions::Strict
        } else {
            NumRotatableBondsOptions::NonStrict
        },
    )
}

/// `calcNumHBD`
pub fn calc_num_hbd(mol: &ROMol) -> Result<u32> {
    static Q: OnceLock<ROMol> = OnceLock::new();
    count_matches(mol, &Q, "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]")
}

/// `calcNumHBA`
pub fn calc_num_hba(mol: &ROMol) -> Result<u32> {
    static Q: OnceLock<ROMol> = OnceLock::new();
    count_matches(
        mol,
        &Q,
        concat!(
            "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$(",
            "[N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]"
        ),
    )
}

/// `calcNumHeteroatoms`
pub fn calc_num_heteroatoms(mol: &ROMol) -> Result<u32> {
    static Q: OnceLock<ROMol> = OnceLock::new();
    count_matches(mol, &Q, "[!#6;!#1]")
}

/// `calcNumAmideBonds`
pub fn calc_num_amide_bonds(mol: &ROMol) -> Result<u32> {
    static Q: OnceLock<ROMol> = OnceLock::new();
    count_matches(mol, &Q, "C(=[O;!R])N")
}

/// `calcNumRings`
pub fn calc_num_rings(mol: &ROMol) -> Result<u32> {
    mol.ring_info().num_rings()
}

/// `calcFractionCSP3`: sp3 carbons (total degree 4) over all carbons.
pub fn calc_fraction_csp3(mol: &ROMol) -> Result<f64> {
    let mut n_csp3 = 0u32;
    let mut n_c = 0u32;
    for i in 0..mol.num_atoms() {
        if mol.atom(i).atomic_num() == 6 {
            n_c += 1;
            if mol.total_degree(i)? == 4 {
                n_csp3 += 1;
            }
        }
    }
    if n_c == 0 {
        return Ok(0.0);
    }
    Ok(f64::from(n_csp3) / f64::from(n_c))
}

fn bond_atoms_are_carbon(mol: &ROMol, bond: i32) -> bool {
    let b = mol.bond(bond as u32);
    mol.atom(b.begin_atom_idx()).atomic_num() == 6 && mol.atom(b.end_atom_idx()).atomic_num() == 6
}

/// `calcNumHeterocycles`: atom rings with at least one non-carbon atom.
pub fn calc_num_heterocycles(mol: &ROMol) -> u32 {
    let mut res = 0;
    for ring in mol.ring_info().atom_rings() {
        if ring.iter().any(|&i| mol.atom(i as u32).atomic_num() != 6) {
            res += 1;
        }
    }
    res
}

/// `calcNumAromaticRings`: rings whose bonds are all aromatic.
pub fn calc_num_aromatic_rings(mol: &ROMol) -> u32 {
    let mut res = 0;
    for ring in mol.ring_info().bond_rings() {
        if ring.iter().all(|&i| mol.bond(i as u32).is_aromatic()) {
            res += 1;
        }
    }
    res
}

/// `calcNumSaturatedRings`: rings whose bonds are all single and not aromatic.
pub fn calc_num_saturated_rings(mol: &ROMol) -> u32 {
    let mut res = 0;
    for ring in mol.ring_info().bond_rings() {
        let saturated = ring.iter().all(|&i| {
            let b = mol.bond(i as u32);
            b.bond_type() == BondType::Single && !b.is_aromatic()
        });
        if saturated {
            res += 1;
        }
    }
    res
}

/// `calcNumAliphaticRings`: rings with at least one non-aromatic bond.
pub fn calc_num_aliphatic_rings(mol: &ROMol) -> u32 {
    let mut res = 0;
    for ring in mol.ring_info().bond_rings() {
        if ring.iter().any(|&i| !mol.bond(i as u32).is_aromatic()) {
            res += 1;
        }
    }
    res
}

/// `calcNumAromaticHeterocycles`
pub fn calc_num_aromatic_heterocycles(mol: &ROMol) -> u32 {
    let mut res = 0;
    for ring in mol.ring_info().bond_rings() {
        let mut count_it = false;
        for &i in ring {
            if !mol.bond(i as u32).is_aromatic() {
                count_it = false;
                break;
            }
            if !count_it && !bond_atoms_are_carbon(mol, i) {
                count_it = true;
            }
        }
        if count_it {
            res += 1;
        }
    }
    res
}

/// `calcNumAromaticCarbocycles`
pub fn calc_num_aromatic_carbocycles(mol: &ROMol) -> u32 {
    let mut res = 0;
    for ring in mol.ring_info().bond_rings() {
        let mut count_it = true;
        for &i in ring {
            if !mol.bond(i as u32).is_aromatic() {
                count_it = false;
                break;
            }
            if !bond_atoms_are_carbon(mol, i) {
                count_it = false;
                break;
            }
        }
        if count_it {
            res += 1;
        }
    }
    res
}

/// `calcNumAliphaticHeterocycles`: rings with an aliphatic bond and a heteroatom bond.
pub fn calc_num_aliphatic_heterocycles(mol: &ROMol) -> u32 {
    let mut res = 0;
    for ring in mol.ring_info().bond_rings() {
        let mut has_aliph = false;
        let mut has_hetero = false;
        for &i in ring {
            if !mol.bond(i as u32).is_aromatic() {
                has_aliph = true;
            }
            if !has_hetero && !bond_atoms_are_carbon(mol, i) {
                has_hetero = true;
            }
        }
        if has_hetero && has_aliph {
            res += 1;
        }
    }
    res
}

/// `calcNumAliphaticCarbocycles`: rings with an aliphatic bond and no heteroatom bond.
pub fn calc_num_aliphatic_carbocycles(mol: &ROMol) -> u32 {
    let mut res = 0;
    for ring in mol.ring_info().bond_rings() {
        let mut has_aliph = false;
        let mut has_hetero = false;
        for &i in ring {
            if !mol.bond(i as u32).is_aromatic() {
                has_aliph = true;
            }
            if !bond_atoms_are_carbon(mol, i) {
                has_hetero = true;
                break;
            }
        }
        if has_aliph && !has_hetero {
            res += 1;
        }
    }
    res
}

/// `calcNumSaturatedHeterocycles`
pub fn calc_num_saturated_heterocycles(mol: &ROMol) -> u32 {
    let mut res = 0;
    for ring in mol.ring_info().bond_rings() {
        let mut count_it = false;
        for &i in ring {
            let b = mol.bond(i as u32);
            if b.bond_type() != BondType::Single || b.is_aromatic() {
                count_it = false;
                break;
            }
            if !count_it && !bond_atoms_are_carbon(mol, i) {
                count_it = true;
            }
        }
        if count_it {
            res += 1;
        }
    }
    res
}

/// `calcNumSaturatedCarbocycles`
pub fn calc_num_saturated_carbocycles(mol: &ROMol) -> u32 {
    let mut res = 0;
    for ring in mol.ring_info().bond_rings() {
        let mut count_it = true;
        for &i in ring {
            let b = mol.bond(i as u32);
            if b.bond_type() != BondType::Single || b.is_aromatic() {
                count_it = false;
                break;
            }
            if !bond_atoms_are_carbon(mol, i) {
                count_it = false;
                break;
            }
        }
        if count_it {
            res += 1;
        }
    }
    res
}

/// Makes sure the SSSR is available (`MolOps::findSSSR` when it is not).
fn ensure_sssr(mol: &mut ROMol) -> Result<()> {
    if !mol.ring_info().is_sssr_or_better() {
        find_sssr(mol, false, false)?;
    }
    Ok(())
}

/// `calcNumSpiroAtoms`: the atoms shared by exactly two rings that share one atom.
pub fn spiro_atoms(mol: &mut ROMol) -> Result<Vec<u32>> {
    ensure_sssr(mol)?;
    let rings = mol.ring_info().atom_rings().to_vec();
    let mut atoms: Vec<u32> = Vec::new();
    for i in 0..rings.len() {
        for j in i + 1..rings.len() {
            let inter = intersect(&rings[i], &rings[j]);
            if inter.len() == 1 {
                let a = inter[0] as u32;
                if !atoms.contains(&a) {
                    atoms.push(a);
                }
            }
        }
    }
    Ok(atoms)
}

/// `calcNumBridgeheadAtoms`: the atoms shared by two rings that share more than one bond.
pub fn bridgehead_atoms(mol: &mut ROMol) -> Result<Vec<u32>> {
    ensure_sssr(mol)?;
    let rings = mol.ring_info().bond_rings().to_vec();
    let n = mol.num_atoms() as usize;
    let mut atoms: Vec<u32> = Vec::new();
    for i in 0..rings.len() {
        for j in i + 1..rings.len() {
            let inter = intersect(&rings[i], &rings[j]);
            if inter.len() > 1 {
                let mut counts = vec![0u32; n];
                for &ii in &inter {
                    let b = mol.bond(ii as u32);
                    counts[b.begin_atom_idx() as usize] += 1;
                    counts[b.end_atom_idx() as usize] += 1;
                }
                for (ti, &c) in counts.iter().enumerate() {
                    if c == 1 && !atoms.contains(&(ti as u32)) {
                        atoms.push(ti as u32);
                    }
                }
            }
        }
    }
    Ok(atoms)
}

/// The molecule with stereochemistry assigned, when it is not assigned already.
fn with_stereo_assigned(mol: &ROMol) -> Result<Option<ROMol>> {
    if mol.has_prop(cp::_StereochemDone) {
        return Ok(None);
    }
    let mut tmol = mol.clone();
    assign_stereochemistry(&mut tmol, true, true, true)?;
    Ok(Some(tmol))
}

/// `numAtomStereoCenters`: atoms that can be stereocenters (specified or not).
pub fn num_atom_stereo_centers(mol: &ROMol) -> Result<u32> {
    let owned = with_stereo_assigned(mol)?;
    let m = owned.as_ref().unwrap_or(mol);
    Ok((0..m.num_atoms())
        .filter(|&i| m.atom(i).has_prop(cp::_ChiralityPossible))
        .count() as u32)
}

/// `numUnspecifiedAtomStereoCenters`: stereocenter candidates without a chiral tag.
pub fn num_unspecified_atom_stereo_centers(mol: &ROMol) -> Result<u32> {
    let owned = with_stereo_assigned(mol)?;
    let m = owned.as_ref().unwrap_or(mol);
    Ok((0..m.num_atoms())
        .filter(|&i| {
            let a = m.atom(i);
            a.has_prop(cp::_ChiralityPossible) && a.chiral_tag() == ChiralType::Unspecified
        })
        .count() as u32)
}
