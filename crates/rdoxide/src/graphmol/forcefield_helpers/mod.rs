//! Force-field helpers (upstream `GraphMol/ForceFieldHelpers`).

pub mod mmff;
pub mod uff;

use crate::error::Result;
use crate::forcefield::ForceField;
use crate::geometry::Point3D;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::threads::par_map;

/// `ForceFieldsHelper::OptimizeMolecule`: (0 if the optimization converged,
/// 1 if more iterations are needed; the final energy)
pub fn optimize_molecule(ff: &mut ForceField, max_iters: u32) -> Result<(i32, f64)> {
    let res = ff.minimize(max_iters, 1e-4, 1e-6)?;
    Ok((res, ff.calc_energy()))
}

/// The positions of a conformer as the flat array of a force field.
pub fn flat_positions(positions: &[Point3D]) -> Vec<f64> {
    positions.iter().flat_map(|p| [p.x, p.y, p.z]).collect()
}

fn write_positions(mol: &mut ROMol, conf_id: i32, flat: &[f64]) -> Result<()> {
    let conf = mol.conformer_mut(conf_id)?;
    for (p, c) in conf.positions_mut().iter_mut().zip(flat.chunks_exact(3)) {
        *p = Point3D::new(c[0], c[1], c[2]);
    }
    Ok(())
}

/// Minimizes the force field of conformer `conf_id` and writes the optimized
/// coordinates back to the molecule.
pub fn optimize_conformer(
    mol: &mut ROMol,
    ff: &mut ForceField,
    conf_id: i32,
    max_iters: u32,
) -> Result<(i32, f64)> {
    let res = optimize_molecule(ff, max_iters)?;
    write_positions(mol, conf_id, ff.positions())?;
    Ok(res)
}

/// `ForceFieldsHelper::OptimizeMoleculeConfs`: optimizes every conformer
/// with a copy of `ff` (whose terms were built for one conformer). The
/// results are in conformer order whatever the number of threads.
pub fn optimize_molecule_confs(
    mol: &mut ROMol,
    ff: &ForceField,
    num_threads: i32,
    max_iters: u32,
) -> Result<Vec<(i32, f64)>> {
    let jobs: Vec<Vec<f64>> = mol
        .conformers()
        .iter()
        .map(|c| flat_positions(c.positions()))
        .collect();
    let results = par_map(
        &jobs,
        num_threads,
        |start| -> Result<(Vec<f64>, i32, f64)> {
            let mut field = ff.clone();
            field.set_positions(start)?;
            let (res, e) = optimize_molecule(&mut field, max_iters)?;
            Ok((field.positions().to_vec(), res, e))
        },
    );
    let mut out = Vec::with_capacity(results.len());
    for (conf, r) in mol.conformers_mut().iter_mut().zip(results) {
        let (flat, res, e) = r?;
        for (p, c) in conf.positions_mut().iter_mut().zip(flat.chunks_exact(3)) {
            *p = Point3D::new(c[0], c[1], c[2]);
        }
        out.push((res, e));
    }
    Ok(out)
}

/// `UFF::UFFOptimizeMolecule`
pub fn uff_optimize_molecule(
    mol: &mut ROMol,
    max_iters: u32,
    vdw_thresh: f64,
    conf_id: i32,
    ignore_interfrag_interactions: bool,
) -> Result<(i32, f64)> {
    let mut ff =
        uff::construct_force_field(mol, vdw_thresh, conf_id, ignore_interfrag_interactions)?;
    optimize_conformer(mol, &mut ff, conf_id, max_iters)
}

/// `UFF::UFFOptimizeMoleculeConfs`
pub fn uff_optimize_molecule_confs(
    mol: &mut ROMol,
    num_threads: i32,
    max_iters: u32,
    vdw_thresh: f64,
    ignore_interfrag_interactions: bool,
) -> Result<Vec<(i32, f64)>> {
    let ff = uff::construct_force_field(mol, vdw_thresh, -1, ignore_interfrag_interactions)?;
    optimize_molecule_confs(mol, &ff, num_threads, max_iters)
}

/// `ForceFieldsHelper::createEmptyForceFieldForMol`: a force field with the
/// points of a conformer and no terms.
pub fn create_empty_force_field_for_mol(mol: &ROMol, conf_id: i32) -> Result<ForceField> {
    let mut field = ForceField::new(3);
    for p in mol.conformer(conf_id)?.positions() {
        field.push_point(&[p.x, p.y, p.z]);
    }
    Ok(field)
}

/// `MMFF::MMFFOptimizeMolecule`: (-1 if parameters were missing, 0 if the
/// optimization converged, 1 if more iterations are required; the energy).
/// The atoms are typed with the given variant.
pub fn mmff_optimize_molecule(
    mol: &mut ROMol,
    max_iters: u32,
    mmff_variant: &str,
    non_bonded_thresh: f64,
    conf_id: i32,
    ignore_interfrag_interactions: bool,
) -> Result<(i32, f64)> {
    let props = mmff::MmffMolProperties::new(mol, mmff_variant)?;
    if !props.is_valid() {
        return Ok((-1, -1.0));
    }
    let mut ff = mmff::construct_force_field_with_props(
        mol,
        &props,
        non_bonded_thresh,
        conf_id,
        ignore_interfrag_interactions,
    )?;
    optimize_conformer(mol, &mut ff, conf_id, max_iters)
}

/// `MMFF::MMFFOptimizeMoleculeConfs`: (-1, -1) for every conformer when
/// parameters are missing.
pub fn mmff_optimize_molecule_confs(
    mol: &mut ROMol,
    num_threads: i32,
    max_iters: u32,
    mmff_variant: &str,
    non_bonded_thresh: f64,
    ignore_interfrag_interactions: bool,
) -> Result<Vec<(i32, f64)>> {
    let props = mmff::MmffMolProperties::new(mol, mmff_variant)?;
    if !props.is_valid() {
        return Ok(vec![(-1, -1.0); mol.conformers().len()]);
    }
    let ff = mmff::construct_force_field_with_props(
        mol,
        &props,
        non_bonded_thresh,
        -1,
        ignore_interfrag_interactions,
    )?;
    optimize_molecule_confs(mol, &ff, num_threads, max_iters)
}
