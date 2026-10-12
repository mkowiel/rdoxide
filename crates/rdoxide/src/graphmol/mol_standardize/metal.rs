//! The `MetalDisconnector`: breaks covalent bonds between metals and
//! non-metals and adjusts the charges of the atoms involved.
//!
//! Ported from RDKit `Code/GraphMol/MolStandardize/Metal.{h,cpp}`.

use std::collections::BTreeMap;
use std::sync::Arc;

use crate::error::{RDKitError, Result};
use crate::graphmol::bond::BondType;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rd_log_info;
use crate::rdgeneral::common_properties;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::rdgeneral::dict::RDProps;

/// `MetalDisconnectorOptions`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MetalDisconnectorOptions {
    /// split Grignard-type complexes
    pub split_grignards: bool,
    /// split metal-aromatic C bonds
    pub split_aromatic_c: bool,
    /// adjust the charges of the ligand atoms
    pub adjust_charges: bool,
    /// remove the dummy atoms that represent haptic bonds
    pub remove_haptic_dummies: bool,
}

impl Default for MetalDisconnectorOptions {
    fn default() -> Self {
        MetalDisconnectorOptions {
            split_grignards: false,
            split_aromatic_c: false,
            adjust_charges: true,
            remove_haptic_dummies: false,
        }
    }
}

impl MetalDisconnectorOptions {
    /// The options of `disconnectOrganometallics`: `{true, true, false, true}`.
    pub fn organometallics() -> Self {
        MetalDisconnectorOptions {
            split_grignards: true,
            split_aromatic_c: true,
            adjust_charges: false,
            remove_haptic_dummies: true,
        }
    }
}

fn smarts(text: &str) -> Result<ROMol> {
    mol_from_smarts(text, &SmartsParserParams::default())?
        .ok_or_else(|| RDKitError::Value(format!("Failed parsing SMARTS: {text}")))
}

/// `MetalDisconnector`
#[derive(Debug, Clone)]
pub struct MetalDisconnector {
    metal_nof: Arc<ROMol>,
    metal_non: Arc<ROMol>,
    metal_dummy: Arc<ROMol>,
    options: MetalDisconnectorOptions,
}

impl Default for MetalDisconnector {
    fn default() -> Self {
        MetalDisconnector::new(MetalDisconnectorOptions::default())
    }
}

/// A non-metal that lost bonds to metals.
#[derive(Debug, Default)]
struct NonMetal {
    cut_bonds: i32,
    bound_metal_indices: Vec<i32>,
}

impl MetalDisconnector {
    pub fn new(options: MetalDisconnectorOptions) -> Self {
        rd_log_info!("Initializing MetalDisconnector");
        let metal_nof = smarts(concat!(
            "[Li,Na,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,",
            "Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,",
            "W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[#7,#8,F]"
        ))
        .expect("metal_nof SMARTS");
        let mut metal_list = String::from(
            "Al,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,Hf,Ta,W,Re,Os,Ir,Pt,Au]~",
        );
        let mut non_metal_list = String::from("Si,P,As,Sb,S,Se,Te,Cl,Br,I,At]");
        metal_list = if options.split_grignards {
            format!("[Li,Na,Mg,K,{metal_list}")
        } else {
            format!("[{metal_list}")
        };
        non_metal_list = if options.split_aromatic_c {
            format!("[B,#6,{non_metal_list}")
        } else {
            format!("[B,C,{non_metal_list}")
        };
        let metal_non = smarts(&format!("{metal_list}{non_metal_list}")).expect("metal_non SMARTS");
        let metal_dummy = smarts(&format!("{metal_list}[*]")).expect("metal dummy SMARTS");
        MetalDisconnector {
            metal_nof: Arc::new(metal_nof),
            metal_non: Arc::new(metal_non),
            metal_dummy: Arc::new(metal_dummy),
            options,
        }
    }

    pub fn options(&self) -> &MetalDisconnectorOptions {
        &self.options
    }

    /// `getMetalNof`
    pub fn metal_nof(&self) -> &ROMol {
        &self.metal_nof
    }

    /// `getMetalNon`
    pub fn metal_non(&self) -> &ROMol {
        &self.metal_non
    }

    /// `setMetalNof`
    pub fn set_metal_nof(&mut self, mol: ROMol) {
        self.metal_nof = Arc::new(mol);
    }

    /// `setMetalNon`
    pub fn set_metal_non(&mut self, mol: ROMol) {
        self.metal_non = Arc::new(mol);
    }

    /// `disconnect(const ROMol &)`
    pub fn disconnect(&self, mol: &ROMol) -> Result<ROMol> {
        let mut res = mol.clone();
        self.disconnect_in_place(&mut res)?;
        Ok(res)
    }

    /// `disconnect(RWMol &)` / `disconnectInPlace`
    pub fn disconnect_in_place(&self, mol: &mut RWMol) -> Result<()> {
        rd_log_info!("Running MetalDisconnector");
        let mut non_metals: BTreeMap<i32, NonMetal> = BTreeMap::new();
        let mut metal_charge_excess: BTreeMap<i32, i32> = BTreeMap::new();
        for query in [&self.metal_nof, &self.metal_non] {
            let matches = substruct_match(mol, query, &SubstructMatchParameters::default())?;
            for m in matches {
                let metal_idx = m[0].1;
                metal_charge_excess.insert(metal_idx, 0);
                let non_idx = m[1].1;
                let bond = mol
                    .get_bond_between_atoms(metal_idx as u32, non_idx as u32)?
                    .ok_or_else(|| RDKitError::invariant("bond not found"))?;
                let bt = bond.bond_type();
                let order = if (BondType::DativeOne as u8..=BondType::DativeR as u8)
                    .contains(&(bt as u8))
                {
                    0
                } else {
                    bond.bond_type_as_double()? as i32
                };
                // disconnect the metal-R bond
                mol.remove_bond(metal_idx as u32, non_idx as u32)?;
                // record the cut bond order and the metal, needed later to
                // adjust the metal charge
                let value = non_metals.entry(non_idx).or_default();
                value.cut_bonds += order;
                if let Err(pos) = value.bound_metal_indices.binary_search(&metal_idx) {
                    value.bound_metal_indices.insert(pos, metal_idx);
                }
                rd_log_info!(
                    "Removed covalent bond between {} and {}",
                    mol.atom(metal_idx as u32).symbol(),
                    mol.atom(non_idx as u32).symbol()
                );
            }
        }
        if self.options.adjust_charges {
            adjust_charges(mol, &mut non_metals, &mut metal_charge_excess)?;
        }
        if self.options.remove_haptic_dummies {
            self.remove_haptic_dummies(mol)?;
        }
        Ok(())
    }

    fn remove_haptic_dummies(&self, mol: &mut RWMol) -> Result<()> {
        let matches =
            substruct_match(mol, &self.metal_dummy, &SubstructMatchParameters::default())?;
        let mut dummies_to_go = Vec::new();
        for m in matches {
            let (metal_idx, dummy_idx) = (m[0].1 as u32, m[1].1 as u32);
            let Some(bond) = mol.get_bond_between_atoms(metal_idx, dummy_idx)? else {
                continue;
            };
            if let Some(v) = bond.get_prop(common_properties::_MolFileBondEndPts) {
                let sprop = v.to_rd_string();
                if sprop.len() > 4 && sprop.starts_with('(') && sprop.ends_with(')') {
                    dummies_to_go.push(dummy_idx);
                }
            }
        }
        dummies_to_go.sort_unstable_by(|a, b| b.cmp(a));
        mol.begin_batch_edit()?;
        for a in dummies_to_go {
            mol.remove_atom(a)?;
        }
        mol.commit_batch_edit()
    }
}

fn adjust_charges(
    mol: &mut RWMol,
    non_metals: &mut BTreeMap<i32, NonMetal>,
    metal_charge_excess: &mut BTreeMap<i32, i32>,
) -> Result<()> {
    let table = PeriodicTable::get();
    for (&idx, nm) in non_metals.iter_mut() {
        let idx = idx as u32;
        // the formal charge is often wrong: work out the most appropriate one
        // from the allowed valences and the current valence/lone electrons;
        // without standard valences the original valence is assumed correct
        // and the non-metal neutral
        let a = mol.atom(idx);
        let valence_before_cut = a.total_valence()? as i32;
        let rad_before_cut = a.num_radical_electrons() as i32;
        let atomic_num = a.atomic_num();
        let mut fc_after_cut = -nm.cut_bonds;
        let mut valence_after_cut = 0;
        let mut lone_electrons = 0;
        let valens = table.get_valence_list(atomic_num)?;
        if !valens.is_empty() && valens[0] != -1 {
            for (i, &v) in valens.iter().enumerate() {
                valence_after_cut = valence_before_cut + rad_before_cut - nm.cut_bonds;
                if valence_after_cut > v {
                    if i + 1 != valens.len() {
                        continue;
                    }
                    valence_after_cut = v;
                    break;
                }
                fc_after_cut = valence_after_cut - v;
                // radicals before and a negative charge now: a carbene-like
                // system (e.g., [Me]-[C]=O), or something silly such as
                // [Me]-[S]=X or [Me]-[P](X)(X)X
                if rad_before_cut % 2 != 0 && fc_after_cut < 0 {
                    fc_after_cut += 1;
                    lone_electrons += 1;
                    // no radical doublets on N and higher
                    let rad = if atomic_num < 7 {
                        rad_before_cut + 1
                    } else {
                        0
                    };
                    mol.atom_mut(idx).set_num_radical_electrons(rad as u32);
                }
                break;
            }
        }
        let total_degree = mol.total_degree(idx)? as i32;
        // do not put a negative charge on sp2 carbon
        if fc_after_cut == -1 && valence_after_cut == total_degree + 1 {
            fc_after_cut = 0;
        }
        mol.atom_mut(idx).set_formal_charge(fc_after_cut);
        if nm.cut_bonds == 0
            || (fc_after_cut == -1 && atomic_num == 6 && valence_after_cut == total_degree + 2)
        {
            // do not take electrons from the metal if it was a dative bond
            // (e.g., [C-]#[O+] coordinated to metal)
            fc_after_cut = 0;
        }
        libcxx_sort_by(&mut nm.bound_metal_indices, |a, b| {
            metal_charge_excess[a] < metal_charge_excess[b]
        });
        fc_after_cut += lone_electrons;
        while fc_after_cut < 0 {
            for &i in &nm.bound_metal_indices {
                // a non-dative bond: the non-metal took electrons from the
                // metal(s), so take electrons from the once-bonded metal(s)
                let was = fc_after_cut;
                fc_after_cut += 1;
                if was >= 0 {
                    break;
                }
                *metal_charge_excess.entry(i).or_insert(0) += 1;
            }
        }
        mol.update_atom_property_cache(idx, true)?;
    }
    // adjust the formal charges of the metal atoms
    for (&idx, &excess) in metal_charge_excess.iter() {
        let idx = idx as u32;
        let a = mol.atom(idx);
        let current_fc = a.formal_charge();
        let valens = table.get_valence_list(a.atomic_num())?;
        // the valence lists contain at least -1
        let max_valence = valens.iter().copied().max().unwrap_or(-1);
        // don't go over the maximum real valence
        if max_valence != -1 && current_fc >= max_valence {
            continue;
        }
        let mut fc_after_cut = excess;
        if current_fc > 0 {
            // a positive original charge on the metal is trusted
            fc_after_cut += current_fc;
        }
        if !valens.is_empty() && valens[0] != -1 {
            for (i, &v) in valens.iter().enumerate() {
                // some metals (e.g. Mg and Ba) end with a catch-all -1
                if v == -1 {
                    continue;
                }
                if fc_after_cut > v {
                    if i + 1 != valens.len() && fc_after_cut >= v {
                        continue;
                    }
                    fc_after_cut = v;
                    break;
                }
            }
        }
        let a = mol.atom_mut(idx);
        if fc_after_cut > current_fc {
            a.set_formal_charge(fc_after_cut);
        }
        // no radical electrons on metals, and none added by sanitization
        a.set_num_radical_electrons(0);
        a.set_num_explicit_hs(0);
        a.set_no_implicit(false);
        mol.update_atom_property_cache(idx, true)?;
    }
    Ok(())
}
