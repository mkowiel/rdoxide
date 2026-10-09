//! Sequence and FASTA parsers: peptides and nucleic acids from one letter codes.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/SequenceParsers.cpp`
//! (2026.09.1); the residue definitions live in the generated
//! `sequence_tables.rs`.
//!
//! Copyright (C) 2015,2016 Greg Landrum and NextMove Software.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::sequence_tables::{AMINO_ACIDS, NUCLEIC_BASES};
use crate::Result;
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::BondType;
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::monomer_info::AtomMonomerInfo;
use crate::graphmol::romol::RWMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

/// A place holding an atom while a residue is being built: the three
/// attachment points handed back to the caller or one of the scratch atoms.
#[derive(Debug, Clone, Copy)]
pub(super) enum Slot {
    R1,
    R2,
    R3,
    A(u8),
}

/// One construction step of a residue (the straight-line statements of the
/// upstream code).
#[derive(Debug, Clone, Copy)]
pub(super) enum Op {
    /// N, CA, C, O, CB with the backbone bonds; `R1` = N, `R2` = C, `A(0)` =
    /// CB. The sign of the argument is the stereo (+1 L, -1 D).
    Backbone(i8),
    /// a new atom named by its four character PDB name
    Atom(Slot, &'static str),
    /// a bond of order 1 or 2
    Bond(Slot, Slot, u8),
    Chiral(Slot, ChiralType),
    /// set the HETATM flag of the following atoms
    Hetero(bool),
    /// `dst = src`
    Alias(Slot, Slot),
}

/// The attachment points of a built residue.
#[derive(Default)]
struct Frame {
    r: [Option<u32>; 3],
    a: [Option<u32>; 32],
}

impl Frame {
    fn slot(&mut self, s: Slot) -> &mut Option<u32> {
        match s {
            Slot::R1 => &mut self.r[0],
            Slot::R2 => &mut self.r[1],
            Slot::R3 => &mut self.r[2],
            Slot::A(i) => &mut self.a[i as usize],
        }
    }

    fn get(&mut self, s: Slot) -> u32 {
        self.slot(s).expect("residue tables use only created atoms")
    }
}

/// The running residue info of a chain: serial numbers from 1, residues from 0.
fn fresh_info(chain: &str) -> AtomMonomerInfo {
    let mut info = AtomMonomerInfo::new_pdb("");
    info.chain_id = chain.to_string();
    let pdb = info.pdb_mut();
    pdb.serial_number = 1;
    pdb.alt_loc = " ".into();
    pdb.insertion_code = " ".into();
    info
}

/// The molecule under construction with the running PDB residue info.
pub(super) struct Builder {
    pub(super) mol: RWMol,
    pub(super) info: AtomMonomerInfo,
}

impl Builder {
    pub(super) fn new() -> Self {
        Builder {
            mol: RWMol::new(),
            info: fresh_info("A"),
        }
    }

    /// Starts numbering a new chain (a fresh `AtomPDBResidueInfo`).
    pub(super) fn start_chain(&mut self, chain: &str) {
        self.info = fresh_info(chain);
    }

    fn residue_number(&self) -> i32 {
        self.info.residue_number
    }

    /// `info.setResidueNumber(resno + 1)`, `setIsHeteroAtom(false)`,
    /// `setResidueName(name)`: the start of a residue.
    pub(super) fn begin_residue(&mut self, name: &str) {
        self.info.residue_number += 1;
        self.info.pdb_mut().is_hetero_atom = false;
        self.info.residue_name = name.to_string();
    }

    /// `CreateAAAtom`: the element follows from the PDB atom name; the info
    /// is a copy of the running one, whose serial number advances.
    pub(super) fn atom(&mut self, name: &str) -> Result<u32> {
        let b = name.as_bytes();
        let num = match (b[0], b[1]) {
            (b' ', b'C') => 6,
            (b' ', b'N') => 7,
            (b' ', b'O') => 8,
            (b' ', b'P') => 15,
            (b' ', b'S') => 16,
            (b'S', b'E') => 34,
            _ => 0,
        };
        let idx = self.mol.add_atom(Atom::new(num), true);
        let mut info = self.info.clone();
        info.name = name.to_string();
        self.mol.atom_mut(idx).set_monomer_info(Some(info));
        self.info.pdb_mut().serial_number += 1;
        Ok(idx)
    }

    /// `CreateAABond`
    pub(super) fn bond(&mut self, a: u32, b: u32, order: u8) -> Result<()> {
        let kind = if order == 2 {
            BondType::Double
        } else {
            BondType::Single
        };
        self.mol.add_bond(a, b, kind)?;
        Ok(())
    }

    fn apply(&mut self, frame: &mut Frame, op: &Op) -> Result<()> {
        match *op {
            Op::Backbone(stereo) => {
                let n = self.atom(" N  ")?;
                let ca = self.atom(" CA ")?;
                let c = self.atom(" C  ")?;
                let o = self.atom(" O  ")?;
                let cb = self.atom(" CB ")?;
                self.bond(n, ca, 1)?;
                self.bond(ca, c, 1)?;
                self.bond(c, o, 2)?;
                self.bond(ca, cb, 1)?;
                match stereo.cmp(&0) {
                    std::cmp::Ordering::Greater => self
                        .mol
                        .atom_mut(ca)
                        .set_chiral_tag(ChiralType::TetrahedralCCW),
                    std::cmp::Ordering::Less => self
                        .mol
                        .atom_mut(ca)
                        .set_chiral_tag(ChiralType::TetrahedralCW),
                    std::cmp::Ordering::Equal => {}
                }
                frame.r[0] = Some(n);
                frame.r[1] = Some(c);
                frame.a[0] = Some(cb);
            }
            Op::Atom(slot, name) => {
                let idx = self.atom(name)?;
                *frame.slot(slot) = Some(idx);
            }
            Op::Bond(a, b, order) => {
                let (a, b) = (frame.get(a), frame.get(b));
                self.bond(a, b, order)?;
            }
            Op::Chiral(slot, tag) => {
                let idx = frame.get(slot);
                self.mol.atom_mut(idx).set_chiral_tag(tag);
            }
            Op::Hetero(flag) => self.info.pdb_mut().is_hetero_atom = flag,
            Op::Alias(dst, src) => {
                let idx = frame.get(src);
                *frame.slot(dst) = Some(idx);
            }
        }
        Ok(())
    }

    fn run(&mut self, frame: &mut Frame, ops: &[Op]) -> Result<()> {
        ops.iter().try_for_each(|op| self.apply(frame, op))
    }

    /// `CreateAminoAcid`: returns the N terminal, C terminal and side chain
    /// attachment atoms (`r1`, `r2`, `r3`).
    pub(super) fn amino_acid(&mut self, code: &str) -> Result<[Option<u32>; 3]> {
        self.begin_residue(code);
        let mut frame = Frame::default();
        if let Some((_, ops)) = AMINO_ACIDS.iter().find(|(name, _)| *name == code) {
            self.run(&mut frame, ops)?;
        }
        Ok(frame.r)
    }

    /// `CreateNucleicAcid`: `code` is the three character residue name
    /// (`"  A"` for RNA, `" DA"` for DNA).
    pub(super) fn nucleic_acid(&mut self, code: &str, p_cap5: bool) -> Result<[Option<u32>; 3]> {
        use Op::*;
        use Slot::*;
        let first = self.residue_number() == 0;
        self.begin_residue(code);
        let mut frame = Frame::default();

        let phosphate: &[Op] = &[
            Atom(A(1), " P  "),
            Atom(A(2), " OP1"),
            Atom(A(3), " OP2"),
            Atom(A(4), " O5'"),
            Bond(A(1), A(2), 2),
            Bond(A(1), A(3), 1),
            Bond(A(1), A(4), 1),
            Alias(R1, A(1)),
        ];
        let capped_phosphate: &[Op] = &[
            Atom(A(0), " OP3"),
            Atom(A(1), " P  "),
            Atom(A(2), " OP1"),
            Atom(A(3), " OP2"),
            Atom(A(4), " O5'"),
            Bond(A(0), A(1), 1),
            Bond(A(1), A(2), 2),
            Bond(A(1), A(3), 1),
            Bond(A(1), A(4), 1),
            Alias(R1, A(0)),
        ];
        let no_phosphate: &[Op] = &[Atom(A(4), " O5'"), Alias(R1, A(4))];
        self.run(
            &mut frame,
            if !first {
                phosphate
            } else if p_cap5 {
                capped_phosphate
            } else {
                no_phosphate
            },
        )?;

        const CW: ChiralType = ChiralType::TetrahedralCW;
        let sugar: &[Op] = &[
            Atom(A(5), " C5'"),
            Atom(A(6), " C4'"),
            Atom(A(7), " O4'"),
            Atom(A(8), " C3'"),
            Atom(A(9), " O3'"),
            Atom(A(10), " C2'"),
            Bond(A(4), A(5), 1),
            Bond(A(5), A(6), 1),
            Bond(A(6), A(7), 1),
            Bond(A(6), A(8), 1),
            Bond(A(8), A(9), 1),
            Bond(A(8), A(10), 1),
            Chiral(A(6), CW),
            Chiral(A(8), CW),
        ];
        self.run(&mut frame, sugar)?;
        let rna = code.as_bytes()[1] == b' ';
        let ribose: &[Op] = if rna {
            &[
                Atom(A(11), " O2'"),
                Atom(A(12), " C1'"),
                Bond(A(10), A(11), 1),
                Bond(A(7), A(12), 1),
                Bond(A(10), A(12), 1),
                Chiral(A(10), CW),
            ]
        } else {
            &[
                Atom(A(12), " C1'"),
                Bond(A(7), A(12), 1),
                Bond(A(10), A(12), 1),
            ]
        };
        self.run(&mut frame, ribose)?;
        frame.r[1] = frame.a[9];

        if let Some((_, ops)) = NUCLEIC_BASES
            .iter()
            .find(|(base, _)| *base == code.as_bytes()[2])
        {
            self.run(&mut frame, ops)?;
        }
        Ok(frame.r)
    }

    /// `CreatePCap3`: the 3' terminal phosphate group.
    pub(super) fn p_cap3(&mut self, prev: u32) -> Result<()> {
        self.info.residue_number += 1;
        self.info.pdb_mut().is_hetero_atom = true;
        self.info.residue_name = "2PO".into();
        let p = self.atom(" P  ")?;
        let op1 = self.atom(" OP1")?;
        let op2 = self.atom(" OP2")?;
        let op3 = self.atom(" OP3")?;
        self.bond(prev, p, 1)?;
        self.bond(p, op1, 2)?;
        self.bond(p, op2, 1)?;
        self.bond(p, op3, 1)
    }

    /// The chain break of a `.`: the next chain letter (up to `Z`) and
    /// numbering from the start.
    fn next_chain(&mut self) {
        if let Some(c) = self.info.chain_id.bytes().next().filter(|&c| c < b'Z') {
            self.info.chain_id = ((c + 1) as char).to_string();
        }
        self.info.residue_number = 0;
    }
}

/// (upper case one letter code, L amino acid, D amino acid)
pub(super) const AMINO_ACID_CODES: [(u8, &str, &str); 20] = [
    (b'A', "ALA", "DAL"),
    (b'C', "CYS", "DCY"),
    (b'D', "ASP", "DAS"),
    (b'E', "GLU", "DGL"),
    (b'F', "PHE", "DPN"),
    (b'G', "GLY", "GLY"),
    (b'H', "HIS", "DHI"),
    (b'I', "ILE", "DIL"),
    (b'K', "LYS", "DLY"),
    (b'L', "LEU", "DLE"),
    (b'M', "MET", "MED"),
    (b'N', "ASN", "DSG"),
    (b'P', "PRO", "DPR"),
    (b'Q', "GLN", "DGN"),
    (b'R', "ARG", "DAR"),
    (b'S', "SER", "DSN"),
    (b'T', "THR", "DTH"),
    (b'V', "VAL", "DVA"),
    (b'W', "TRP", "DTR"),
    (b'Y', "TYR", "DTY"),
];

/// The three letter code of a one letter amino acid code; `lower_d` selects
/// the D-amino acids for lower case letters (glycine has no D form).
fn amino_acid_code(c: u8, lower_d: bool) -> Option<&'static str> {
    let (_, l, d) = AMINO_ACID_CODES
        .iter()
        .find(|(u, _, _)| *u == c.to_ascii_uppercase())?;
    Some(if c.is_ascii_lowercase() && lower_d {
        d
    } else {
        l
    })
}

/// `AASequenceToMol`; `None` for an unknown character.
fn aa_sequence_to_mol(seq: &[u8], lower_d: bool) -> Result<Option<RWMol>> {
    let mut b = Builder::new();
    let mut prev: Option<u32> = None;
    for &c in seq {
        let mut r = [None; 3];
        match c {
            b'\n' | b'\r' | b'-' => continue,
            b' ' | b'\t' => {}
            b'.' => {
                if let Some(p) = prev.take() {
                    let oxt = b.atom(" OXT")?;
                    b.bond(p, oxt, 1)?;
                    b.next_chain();
                }
                continue;
            }
            _ => match amino_acid_code(c, lower_d) {
                Some(code) => r = b.amino_acid(code)?,
                None => return Ok(None),
            },
        }
        if let (Some(p), Some(r1)) = (prev, r[0]) {
            b.bond(p, r1, 1)?;
        }
        prev = r[1];
    }
    if let Some(p) = prev {
        let oxt = b.atom(" OXT")?;
        b.bond(p, oxt, 1)?;
    }
    Ok(Some(b.mol))
}

/// `NASequenceToMol`; `None` for an unknown character.
fn na_sequence_to_mol(seq: &[u8], dna: bool, p_cap5: bool, p_cap3: bool) -> Result<Option<RWMol>> {
    let mut b = Builder::new();
    let mut prev: Option<u32> = None;
    for &c in seq {
        let mut r = [None; 3];
        match c {
            b'\n' | b'\r' | b'-' => continue,
            b' ' | b'\t' => {}
            b'.' => {
                if let Some(p) = prev.take() {
                    if p_cap3 {
                        b.p_cap3(p)?;
                    }
                    b.next_chain();
                }
                continue;
            }
            _ => {
                let base = c.to_ascii_uppercase();
                if !b"ACGTU".contains(&base) {
                    return Ok(None);
                }
                let code = format!("{}{}", if dna { " D" } else { "  " }, base as char);
                r = b.nucleic_acid(&code, p_cap5)?;
            }
        }
        if let (Some(p), Some(r1)) = (prev, r[0]) {
            b.bond(p, r1, 1)?;
        }
        prev = r[1];
    }
    if let (Some(p), true) = (prev, p_cap3) {
        b.p_cap3(p)?;
    }
    Ok(Some(b.mol))
}

/// `SequenceToMol`: `flavor` 0/1 protein (L/D), 2-5 RNA and 6-9 DNA (the
/// bits 1 and 2 of `flavor - 2` mod 4 add the 5' and 3' phosphate caps).
/// `Ok(None)` for unknown characters or flavors.
pub fn sequence_to_mol(seq: &str, sanitize: bool, flavor: i32) -> Result<Option<RWMol>> {
    let seq = seq.as_bytes();
    // the C string ends at the first NUL
    let seq = &seq[..seq.iter().position(|&c| c == 0).unwrap_or(seq.len())];
    let mol = match flavor {
        0 | 1 => aa_sequence_to_mol(seq, flavor == 1)?,
        2..=9 => {
            let (dna, caps) = (flavor >= 6, (flavor - 2) % 4);
            na_sequence_to_mol(seq, dna, caps & 1 != 0, caps & 2 != 0)?
        }
        _ => return Ok(None),
    };
    let Some(mut mol) = mol else {
        return Ok(None);
    };
    if sanitize {
        sanitize_mol(&mut mol, SANITIZE_ALL).map_err(|(_, e)| e)?;
    }
    Ok(Some(mol))
}

/// `FASTAToMol`: a leading `>` line names the molecule.
pub fn fasta_to_mol(seq: &str, sanitize: bool, flavor: i32) -> Result<Option<RWMol>> {
    let mut body = seq;
    let mut title = "";
    if let Some(rest) = seq.strip_prefix('>') {
        let end = rest.find(['\n', '\r']).unwrap_or(rest.len());
        (title, body) = rest.split_at(end);
    }
    let Some(mut mol) = sequence_to_mol(body, sanitize, flavor)? else {
        return Ok(None);
    };
    if !title.is_empty() {
        mol.set_prop(cp::_Name, title.to_string(), false)?;
    }
    Ok(Some(mol))
}
