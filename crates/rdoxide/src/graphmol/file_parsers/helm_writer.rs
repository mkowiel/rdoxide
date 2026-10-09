//! HELM writer for molecules carrying PDB residue information.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/SequenceWriters.cpp`
//! (`MolToHELM`, 2026.09.1).
//!
//! Copyright (C) 2015,2016 Greg Landrum and NextMove Software.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::sequence_tables::HELM_PEPTIDE_WRITER;
use crate::graphmol::bond::BondType;
use crate::graphmol::monomer_info::{AtomMonomerInfo, AtomMonomerType};
use crate::graphmol::romol::ROMol;

/// The most chains HELM output supports (`PEPTIDE1` ... `PEPTIDE9`).
const MAX_CHAINS: usize = 9;

fn insertion_code(i: &AtomMonomerInfo) -> &str {
    i.pdb.as_ref().map_or("", |p| p.insertion_code.as_str())
}

/// `getHELMAAMonomer`
fn peptide_monomer(residue: &str) -> Option<&'static str> {
    HELM_PEPTIDE_WRITER
        .iter()
        .find(|(res, _)| *res == residue)
        .map(|(_, mono)| *mono)
}

/// `getHELMNAMonomer`: `R(A)` for RNA and `[dR](A)` for DNA nucleotides.
fn nucleic_monomer(residue: &str) -> Option<String> {
    let [b' ', kind, base] = *residue.as_bytes() else {
        return None;
    };
    if !b"ACGTU".contains(&base) {
        return None;
    }
    match kind {
        b' ' => Some(format!("R({})", base as char)),
        b'D' => Some(format!("[dR]({})", base as char)),
        _ => None,
    }
}

/// `IsEupeptideBond`: whether the bond is an ordinary backbone link, which
/// HELM does not list in the connection section.
fn is_backbone_bond(src: &AtomMonomerInfo, dst: &AtomMonomerInfo) -> bool {
    if src.chain_id != dst.chain_id {
        return false;
    }
    let (sres, dres) = (src.residue_number, dst.residue_number);
    let (sins, dins) = (insertion_code(src), insertion_code(dst));
    // "later residue" for the first atom pair, "earlier" for the reversed one
    let later = || {
        if sres != dres {
            dres > sres
        } else {
            dins > sins
        }
    };
    let earlier = || {
        if sres != dres {
            dres < sres
        } else {
            dins < sins
        }
    };
    match (src.name.as_str(), dst.name.as_str()) {
        (" C  ", " N  ") => later(),
        (" N  ", " C  ") => earlier(),
        (" O3'", " P  ") | (" P  ", " O3'") => later(),
        _ => false,
    }
}

/// `IsSupportedHELMBond`
fn is_supported_bond(src: &AtomMonomerInfo, dst: &AtomMonomerInfo) -> bool {
    let cysteine = |i: &AtomMonomerInfo| matches!(i.residue_name.as_str(), "CYS" | "DCY");
    let disulfide = src.name == " SG " && dst.name == " SG " && cysteine(src) && cysteine(dst);
    let peptide = matches!(
        (src.name.as_str(), dst.name.as_str()),
        (" N  ", " C  ") | (" C  ", " N  ")
    );
    disulfide || peptide
}

fn same_residue(a: &AtomMonomerInfo, b: &AtomMonomerInfo) -> bool {
    a.residue_number == b.residue_number
        && a.residue_name == b.residue_name
        && a.chain_id == b.chain_id
        && insertion_code(a) == insertion_code(b)
}

/// `FindHELMAtom`: `PEPTIDEn` and the `m:Rk` position of an attachment atom.
fn find_atom(chains: &[Vec<&AtomMonomerInfo>], info: &AtomMonomerInfo) -> Option<(String, String)> {
    let group = match info.name.as_str() {
        " SG " => '3',
        " N  " => '1',
        " C  " => '2',
        _ => return None,
    };
    chains.iter().enumerate().skip(1).find_map(|(i, chain)| {
        let j = chain.iter().position(|t| same_residue(t, info))?;
        Some((format!("PEPTIDE{i}"), format!("{}:R{group}", j + 1)))
    })
}

/// `NameHELMBond`
fn name_bond(
    chains: &[Vec<&AtomMonomerInfo>],
    src: &AtomMonomerInfo,
    dst: &AtomMonomerInfo,
) -> Option<String> {
    let (id1, pos1) = find_atom(chains, src)?;
    let (id2, pos2) = find_atom(chains, dst)?;
    Some(format!("{id1},{id2},{pos1}-{pos2}"))
}

/// Where the nucleic acid writer is in `P`/base/`P` sequences (upstream's
/// `Pstate`; its state 3, "after 5'-base", is never entered and is left out).
#[derive(Clone, Copy, PartialEq, Eq)]
enum NaState {
    /// right after a 5' phosphate
    After5Phosphate,
    /// after a base
    AfterBase,
    /// after the phosphate following a base
    AfterBasePhosphate,
    /// after the terminal 3' phosphate
    Terminal,
}

/// `MolToHELM`: an empty string if the molecule cannot be expressed.
pub fn mol_to_helm(mol: &ROMol) -> String {
    write_helm(mol).unwrap_or_default()
}

/// Opening of a new polymer when the chain changes (or none is open yet).
struct Writer<'a> {
    out: String,
    chains: Vec<Vec<&'a AtomMonomerInfo>>,
    id: usize,
    peptide: bool,
    first: bool,
    chain: String,
}

impl Writer<'_> {
    /// Starts `PEPTIDEn{`/`RNAn{` (closing the previous polymer); `None` when
    /// there are too many chains.
    fn open(&mut self, kind: &str, chain: &str, extra: &str) -> Option<()> {
        if self.first {
            self.out = format!("{kind}1{{{extra}");
            self.first = false;
        } else {
            if self.id == MAX_CHAINS {
                return None;
            }
            self.id += 1;
            self.out.push_str(&format!("}}|{kind}{}{{{extra}", self.id));
        }
        self.chain = chain.to_string();
        self.peptide = kind == "PEPTIDE";
        Some(())
    }
}

fn write_helm(mol: &ROMol) -> Option<String> {
    use NaState::*;
    let mut w = Writer {
        out: String::new(),
        chains: vec![Vec::new(); MAX_CHAINS + 1],
        id: 1,
        peptide: true,
        first: true,
        chain: String::new(),
    };
    let mut state: Option<NaState> = None;

    // first pass: the monomers
    for atom in mol.atoms() {
        // HELM needs PDB residue information on every atom
        let info = atom
            .monomer_info()
            .filter(|i| i.monomer_type == AtomMonomerType::PdbResidue)?;
        let new_chain = info.chain_id != w.chain;
        match info.name.as_str() {
            " CA " => {
                let mono = peptide_monomer(&info.residue_name)?;
                if w.first || !w.peptide || new_chain {
                    w.open("PEPTIDE", &info.chain_id, "")?;
                } else {
                    w.out.push('.');
                }
                w.out.push_str(mono);
                w.chains[w.id].push(info);
            }
            " N  " if info.residue_name == "NH2" => {
                if w.first || !w.peptide {
                    return None;
                }
                w.out.push_str(".[am]");
            }
            " C  " if info.residue_name == "ACE" => {
                if w.first || new_chain {
                    w.open("PEPTIDE", &info.chain_id, "[ac]")?;
                } else {
                    return None;
                }
                w.chains[w.id].push(info);
            }
            " C1'" => {
                let mono = nucleic_monomer(&info.residue_name)?;
                if w.first || w.peptide || new_chain {
                    w.open("RNA", &info.chain_id, "")?;
                    state = Some(AfterBase);
                } else if matches!(state, Some(After5Phosphate | AfterBasePhosphate)) {
                    w.out.push('.');
                    state = Some(AfterBase);
                } else {
                    return None;
                }
                w.out.push_str(&mono);
            }
            " P  " => {
                if nucleic_monomer(&info.residue_name).is_some() {
                    if w.first || w.peptide || new_chain {
                        w.open("RNA", &info.chain_id, "P")?;
                        state = Some(After5Phosphate);
                    } else {
                        match state {
                            Some(AfterBase) => {
                                w.out.push('P');
                                state = Some(AfterBasePhosphate);
                            }
                            _ => return None,
                        }
                    }
                } else if info.residue_name == "2PO" {
                    if w.first || w.peptide || new_chain {
                        return None;
                    }
                    match state {
                        Some(AfterBase) => w.out.push('P'),
                        _ => return None,
                    }
                    state = Some(Terminal);
                }
            }
            _ => {}
        }
    }
    if w.first {
        return None;
    }
    w.out.push_str("}$");

    // second pass: the connections that are not the plain backbone
    let mut connections: Vec<String> = Vec::new();
    for bond in mol.bonds() {
        let (Some(binfo), Some(einfo)) = (
            mol.atom(bond.begin_atom_idx()).monomer_info(),
            mol.atom(bond.end_atom_idx()).monomer_info(),
        ) else {
            continue;
        };
        if same_residue(binfo, einfo) {
            continue;
        }
        if bond.bond_type() != BondType::Single {
            return None;
        }
        if is_backbone_bond(binfo, einfo) {
            continue;
        }
        if !is_supported_bond(binfo, einfo) {
            return None;
        }
        connections.push(name_bond(&w.chains, binfo, einfo)?);
    }
    w.out.push_str(&connections.join("|"));
    w.out.push_str("$$$");
    Some(w.out)
}
