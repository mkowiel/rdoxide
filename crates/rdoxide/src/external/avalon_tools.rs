//! Avalon fingerprints for RDKit molecules (upstream `External/AvalonTools/AvalonTools.cpp`).
//! The molecule is written as a MOL block and read by the Avalon reader, as upstream does.

use crate::datastructs::bit_vect::{BitVect, ExplicitBitVect};
use crate::external::avalon::canonizer::{CENTER_STEREO, DB_STEREO, can_smiles};
use crate::external::avalon::fingerprint::{
    ACCUMULATE_BITS, USE_ALL_FEATURES, USE_DY_AROMATICITY, set_fingerprint_counts_with_focus,
};
use crate::external::avalon::layout::{layout_molecule, recolor_molecule};
use crate::external::avalon::reaccs::Molecule;
use crate::external::avalon::reaccsio::{mol_str_to_mol, mol_to_mol_str};
use crate::external::avalon::smi2mol::smi_to_mol;
use crate::external::avalon::smiles::{DY_AROMATICITY, ISOMERIC_SMILES, mol_to_smi};
use crate::geometry::Point3D;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::file_parsers::mol_file_parser::{MolFileParserParams, mol_from_mol_block};
use crate::graphmol::file_parsers::mol_file_writer::{MolWriterParams, mol_to_mol_block};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// the default `bitFlags` of `GetAvalonFP` (all SSS features plus the non-SSS ones)
pub const DEFAULT_BIT_FLAGS: i32 = 15761407;

fn molblock_to_reaccs(molblock: &str) -> Result<Molecule> {
    mol_str_to_mol(molblock).ok_or_else(|| {
        RDKitError::Value("could not convert the molecule to the Avalon format".into())
    })
}

fn mol_to_reaccs(mol: &ROMol) -> Result<Molecule> {
    molblock_to_reaccs(&mol_to_mol_block(mol, &MolWriterParams::default(), -1)?)
}

/// `getFp`: the bytes of the fingerprint. Unless it is a query, a second pass with Daylight style
/// aromaticity is OR-ed in (the molecule keeps the changes of the first pass, as upstream).
fn fingerprint_bytes(
    mp: &mut Molecule,
    n_bytes: usize,
    is_query: bool,
    bit_flags: i32,
) -> Result<Vec<u8>> {
    let n_bytes = n_bytes.next_multiple_of(4);
    let mut bytes = vec![0u8; n_bytes];
    let passes: &[(bool, i32)] = if is_query {
        &[(true, 0)]
    } else {
        &[(false, 0), (false, ACCUMULATE_BITS | USE_DY_AROMATICITY)]
    };
    for &(as_query, fpflags) in passes {
        let mut counts = vec![0i32; n_bytes * 8];
        set_fingerprint_counts_with_focus(mp, &mut counts, bit_flags, as_query, fpflags, 0)
            .map_err(|e| RDKitError::Value(e.into()))?;
        for (i, &c) in counts.iter().enumerate() {
            if c > 0 {
                bytes[(i / 8) % n_bytes] |= 1 << (i % 8);
            }
        }
    }
    Ok(bytes)
}

/// `getAvalonFP`: sets the bits of `res` (which is cleared first if `reset_vect`)
pub fn get_avalon_fp_into(
    mol: &ROMol,
    res: &mut ExplicitBitVect,
    n_bits: u32,
    is_query: bool,
    reset_vect: bool,
    bit_flags: i32,
) -> Result<()> {
    let mut mp = mol_to_reaccs(mol)?;
    let bytes = fingerprint_bytes(&mut mp, (n_bits / 8) as usize, is_query, bit_flags)?;
    if reset_vect {
        *res = ExplicitBitVect::new(n_bits, false);
    }
    for (i, b) in bytes.iter().enumerate() {
        for j in 0..8 {
            if b & (1 << j) != 0 {
                res.set_bit((i * 8 + j) as u32)?;
            }
        }
    }
    Ok(())
}

/// `GetAvalonFP`
pub fn get_avalon_fp(
    mol: &ROMol,
    n_bits: u32,
    is_query: bool,
    bit_flags: i32,
) -> Result<ExplicitBitVect> {
    let mut res = ExplicitBitVect::new(n_bits, false);
    get_avalon_fp_into(mol, &mut res, n_bits, is_query, false, bit_flags)?;
    Ok(res)
}

/// `GetAvalonFP` of a MOL block
pub fn get_avalon_fp_from_molblock(
    molblock: &str,
    n_bits: u32,
    is_query: bool,
    bit_flags: i32,
) -> Result<ExplicitBitVect> {
    let mut mp = molblock_to_reaccs(molblock)?;
    let bytes = fingerprint_bytes(&mut mp, (n_bits / 8) as usize, is_query, bit_flags)?;
    let mut res = ExplicitBitVect::new(n_bits, false);
    for (i, b) in bytes.iter().enumerate() {
        for j in 0..8 {
            if b & (1 << j) != 0 {
                res.set_bit((i * 8 + j) as u32)?;
            }
        }
    }
    Ok(res)
}

/// `GetAvalonCountFP`: the feature multiplicities per bit position (`getCountFp`).
///
/// Upstream does not forward `bitFlags` to the counting code, which always uses `USE_ALL_FEATURES`
/// (32767); the argument is accepted and ignored here as well.
pub fn get_avalon_count_fp(
    mol: &ROMol,
    n_bits: u32,
    is_query: bool,
    bit_flags: i32,
) -> Result<Vec<u32>> {
    let _ = bit_flags;
    let bit_flags = USE_ALL_FEATURES;
    let mut mp = mol_to_reaccs(mol)?;
    let mut counts = vec![0i32; n_bits as usize];
    set_fingerprint_counts_with_focus(&mut mp, &mut counts, bit_flags, is_query, 0, 0)
        .map_err(|e| RDKitError::Value(e.into()))?;
    if !is_query {
        set_fingerprint_counts_with_focus(
            &mut mp,
            &mut counts,
            bit_flags,
            true,
            ACCUMULATE_BITS | USE_DY_AROMATICITY,
            0,
        )
        .map_err(|e| RDKitError::Value(e.into()))?;
    }
    Ok(counts.iter().map(|&c| c.max(0) as u32).collect())
}

/// `GetAvalonFPAsWords`: the fingerprint as 32 bit words
pub fn get_avalon_fp_as_words(
    mol: &ROMol,
    n_bits: u32,
    is_query: bool,
    bit_flags: i32,
) -> Result<Vec<u32>> {
    let mut mp = mol_to_reaccs(mol)?;
    let bytes = fingerprint_bytes(&mut mp, (n_bits / 8) as usize, is_query, bit_flags)?;
    Ok(bytes
        .chunks(4)
        .map(|c| {
            // upstream packs signed `char`s, so bytes >= 0x80 sign-extend into the higher bytes
            let sx = |i: usize| i32::from(c.get(i).copied().unwrap_or(0) as i8);
            (sx(0) | (sx(1) << 8) | (sx(2) << 16) | (sx(3) << 24)) as u32
        })
        .collect())
}

fn string_to_reaccs(data: &str, is_smiles: bool) -> Option<Molecule> {
    if is_smiles {
        smi_to_mol(data, DY_AROMATICITY).ok().flatten()
    } else {
        mol_str_to_mol(data)
    }
}

/// `getAvalonFP(data, isSmiles, ...)`: `None` if the data cannot be read
pub fn get_avalon_fp_from_data(
    data: &str,
    is_smiles: bool,
    n_bits: u32,
    is_query: bool,
    bit_flags: i32,
) -> Result<Option<ExplicitBitVect>> {
    let Some(mut mp) = string_to_reaccs(data, is_smiles) else {
        return Ok(None);
    };
    let bytes = fingerprint_bytes(&mut mp, (n_bits / 8) as usize, is_query, bit_flags)?;
    let mut res = ExplicitBitVect::new(n_bits, false);
    for (i, b) in bytes.iter().enumerate() {
        for j in 0..8 {
            if b & (1 << j) != 0 {
                res.set_bit((i * 8 + j) as u32)?;
            }
        }
    }
    Ok(Some(res))
}

/// `getAvalonFP(data, isSmiles, words, ...)`
pub fn get_avalon_fp_as_words_from_data(
    data: &str,
    is_smiles: bool,
    n_bits: u32,
    is_query: bool,
    bit_flags: i32,
) -> Result<Option<Vec<u32>>> {
    let Some(mut mp) = string_to_reaccs(data, is_smiles) else {
        return Ok(None);
    };
    let bytes = fingerprint_bytes(&mut mp, (n_bits / 8) as usize, is_query, bit_flags)?;
    Ok(Some(
        bytes
            .chunks(4)
            .map(|c| {
                let sx = |i: usize| i32::from(c.get(i).copied().unwrap_or(0) as i8);
                (sx(0) | (sx(1) << 8) | (sx(2) << 16) | (sx(3) << 24)) as u32
            })
            .collect(),
    ))
}

/// `getAvalonCountFP(data, isSmiles, ...)`; unlike the molecule flavour this one does use `bit_flags`
pub fn get_avalon_count_fp_from_data(
    data: &str,
    is_smiles: bool,
    n_bits: u32,
    is_query: bool,
    bit_flags: i32,
) -> Result<Option<Vec<u32>>> {
    let Some(mut mp) = string_to_reaccs(data, is_smiles) else {
        return Ok(None);
    };
    let mut counts = vec![0i32; n_bits as usize];
    let err = |e: &str| RDKitError::Value(e.into());
    set_fingerprint_counts_with_focus(&mut mp, &mut counts, bit_flags, is_query, 0, 0)
        .map_err(err)?;
    if !is_query {
        set_fingerprint_counts_with_focus(
            &mut mp,
            &mut counts,
            bit_flags,
            true,
            ACCUMULATE_BITS | USE_DY_AROMATICITY,
            0,
        )
        .map_err(err)?;
    }
    Ok(Some(counts.iter().map(|&c| c.max(0) as u32).collect()))
}

/// `getCanonSmiles(data, isSmiles, flags)`; an empty string if no SMILES could be made
pub fn get_canon_smiles_from_data(data: &str, is_smiles: bool, flags: i32) -> String {
    let flags = if flags == -1 {
        DB_STEREO | CENTER_STEREO
    } else {
        flags
    };
    let can = if is_smiles {
        can_smiles(data, flags)
    } else {
        string_to_reaccs(data, false).and_then(|mut mp| {
            let smiles = mol_to_smi(&mut mp, ISOMERIC_SMILES).ok().flatten()?;
            can_smiles(&smiles, flags)
        })
    };
    can.unwrap_or_default()
}

/// `getCanonSmiles(mol, flags)`
pub fn get_canon_smiles(mol: &ROMol, flags: i32) -> Result<String> {
    if mol.num_conformers() == 0 {
        let mut m = mol.clone();
        let smi = mol_to_smiles(
            &mut m,
            &SmilesWriteParams {
                do_isomeric_smiles: true,
                ..SmilesWriteParams::default()
            },
        )?;
        Ok(get_canon_smiles_from_data(&smi, true, flags))
    } else {
        let mb = mol_to_mol_block(mol, &MolWriterParams::default(), -1)?;
        Ok(get_canon_smiles_from_data(&mb, false, flags))
    }
}

fn reaccs_get_coords(mut mp: Molecule) -> Molecule {
    recolor_molecule(&mut mp);
    layout_molecule(&mp)
}

/// `set2DCoords(data, isSmiles)`: a MOL block with Avalon 2D coordinates (empty if the data cannot be read)
pub fn set_2d_coords_from_data(data: &str, is_smiles: bool) -> String {
    match string_to_reaccs(data, is_smiles) {
        Some(mp) => mol_to_mol_str(&reaccs_get_coords(mp)),
        None => String::new(),
    }
}

/// `set2DCoords(mol, clearConfs)`: adds a 2D conformer made by Avalon and returns its id
pub fn set_2d_coords(mol: &mut ROMol, clear_confs: bool) -> Result<u32> {
    let n = mol.num_atoms() as usize;
    let orig_has_hs = (0..mol.num_atoms()).any(|i| {
        let a = mol.atom(i);
        a.atomic_num() == 1 && a.isotope() == 0
    });
    let smiles = mol_to_smiles(
        mol,
        &SmilesWriteParams {
            do_isomeric_smiles: true,
            ..SmilesWriteParams::default()
        },
    )?;
    let mp = string_to_reaccs(&smiles, true)
        .ok_or_else(|| RDKitError::Value("could not build a molecule".into()))?;
    let mp2 = reaccs_get_coords(mp);
    if mp2.atoms.len() < n {
        return Err(RDKitError::Value("Avalon lost atoms".into()));
    }
    let mut matched = vec![0usize; n];
    if orig_has_hs || mp2.atoms.len() > n {
        let mb = mol_to_mol_str(&mp2);
        let params = MolFileParserParams {
            sanitize: false,
            remove_hs: false,
            ..MolFileParserParams::default()
        };
        let avmol = mol_from_mol_block(&mb, &params)?.ok_or_else(|| {
            RDKitError::Value("could not parse mol block from avalon toolkit".into())
        })?;
        let ms = substruct_match(&avmol, mol, &SubstructMatchParameters::default())?;
        let m = ms.first().ok_or_else(|| {
            RDKitError::Value("no substructure match found between avalon mol and input mol".into())
        })?;
        for &(q, t) in m {
            matched[q as usize] = t as usize;
        }
    } else {
        let order = match mol.get_prop(cp::_smilesAtomOutputOrder) {
            Some(RDValue::VecUnsignedInt(v)) => v.clone(),
            Some(RDValue::VecInt(v)) => v.iter().map(|&i| i as u32).collect(),
            _ => return Err(RDKitError::Key(cp::_smilesAtomOutputOrder.into())),
        };
        for (ai, m) in matched.iter_mut().enumerate() {
            *m = order[ai] as usize;
        }
    }
    let mut conf = Conformer::new(n);
    conf.set_3d(false);
    for i in 0..n {
        let a = &mp2.atoms[matched[i]];
        conf.set_atom_pos(i as u32, Point3D::new(f64::from(a.x), f64::from(a.y), 0.0))?;
    }
    if clear_confs {
        mol.clear_conformers();
        conf.set_id(0);
        mol.add_conformer(conf, false)
    } else {
        mol.add_conformer(conf, true)
    }
}

/// the structure checker entry points (`initCheckMol`, `getCheckMolLog`, `closeCheckMolFiles`)
pub use crate::external::avalon::struchk::{
    close_check_mol_files, get_msg_list as get_check_mol_log, init_check_mol,
};

/// `checkMolString(data, isSmiles, &mp)`: runs the structure checker; the molecule is `None` if the data
/// could not be read (the result is then `BAD_MOLECULE`)
pub fn check_mol_string_to_reaccs(data: &str, is_smiles: bool) -> (i32, Option<Molecule>) {
    use crate::external::avalon::struchk::{BAD_MOLECULE, free_msg_list, run_struchk};
    free_msg_list();
    match string_to_reaccs(data, is_smiles) {
        Some(mut mp) => {
            let errs = run_struchk(&mut mp);
            (errs, Some(mp))
        }
        None => (BAD_MOLECULE, None),
    }
}

/// `checkMolString(data, isSmiles)`: the result bits and the checked molecule as MOL block (empty if unreadable)
pub fn check_mol_string(data: &str, is_smiles: bool) -> (i32, String) {
    let (errs, mp) = check_mol_string_to_reaccs(data, is_smiles);
    (errs, mp.map(|m| mol_to_mol_str(&m)).unwrap_or_default())
}

fn reaccs_to_romol(mp: &Molecule) -> Result<Option<ROMol>> {
    let mb = mol_to_mol_str(mp);
    mol_from_mol_block(&mb, &MolFileParserParams::default())
}

/// `checkMol(errs, data, isSmiles)`
pub fn check_mol_data(data: &str, is_smiles: bool) -> Result<(i32, Option<ROMol>)> {
    let (errs, mp) = check_mol_string_to_reaccs(data, is_smiles);
    match mp {
        Some(mp) => Ok((errs, reaccs_to_romol(&mp)?)),
        None => Ok((errs, None)),
    }
}

/// `checkMol(errs, mol)`
pub fn check_mol(mol: &ROMol) -> Result<(i32, Option<ROMol>)> {
    use crate::external::avalon::struchk::{free_msg_list, run_struchk};
    free_msg_list();
    let mut mp = mol_to_reaccs(mol)?;
    let errs = run_struchk(&mut mp);
    Ok((errs, reaccs_to_romol(&mp)?))
}
