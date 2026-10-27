//! `smi2mol.c` (the MOL to SMILES half): `MOLToSMI` and its helpers. Copyright 2001-2011 Novartis
//! Pharma AG.

use super::graph::{combine_rings, ring_list};
use super::perceive::{
    AROMATIC, CHAIN, DOUBLE, Neighbourhood, RING, SINGLE, TRIPLE, perceive_dy_aromaticity,
    setup_neighbourhood,
};
use super::reaccs::{CIS_TRANS_EITHER, CIS_TRANS_SWAPPED, Molecule, NONE};
use super::tables::PTABLE;
use super::utilities::{
    DOWN, UP, atom_symbol_match, compute_implicit_h, get_num_property, make_hydrogens_implicit,
    reset_colors, sort_bonds_with_numbering, sort_neighbours_by_numbering,
};

/// flags of `smi2mol.h`
pub const SMARTS_PERCEPTION: i32 = 0x0002;
pub const ISOMERIC_SMILES: i32 = 0x0004;
pub const CLASSIC_AROMATICITY: i32 = 0x0020;
pub const AROMATIC_SMILES: i32 = 0x0040;
pub const TO_MONOMER: i32 = 0x1000;
pub const DY_AROMATICITY: i32 = 0x2000;
pub const USE_Z: i32 = 0x4000;

/// bond types of `mdl.h` not in `perceive`
const SINGLE_DOUBLE: i32 = 5;
const SINGLE_AROMATIC: i32 = 6;
const DOUBLE_AROMATIC: i32 = 7;
const ANY_BOND: i32 = 8;

const DY_ABOVE_FIRST: i32 = 1;
const DY_BELOW_FIRST: i32 = 2;
const DY_ABOVE_SECOND: i32 = 4;
const DY_BELOW_SECOND: i32 = 8;
const DY_ABOVE_BOTH: i32 = DY_ABOVE_FIRST | DY_ABOVE_SECOND;
const DY_BELOW_BOTH: i32 = DY_BELOW_FIRST | DY_BELOW_SECOND;
const MAXOPEN: usize = 99;
const VISITED: i32 = 1;

/// `StrAppend`: `result` + `sp` + `cp`; with no `result` and no `sp`, `cp` is taken as it is
fn str_append(result: Option<String>, cp: Option<String>, sp: Option<&str>) -> Option<String> {
    if cp.is_none() && sp.is_none() {
        return result;
    }
    match result {
        Some(mut r) => {
            r.push_str(sp.unwrap_or(""));
            r.push_str(cp.as_deref().unwrap_or(""));
            Some(r)
        }
        None if sp.is_none() => cp,
        None => Some(format!("{}{}", sp.unwrap_or(""), cp.unwrap_or_default())),
    }
}

/// `BndAppend`
fn bnd_append(
    result: Option<String>,
    cph: Option<String>,
    btype: i32,
    from_first: bool,
    topography: i32,
    bcolor: i32,
    expect_smarts: bool,
) -> Option<String> {
    let ring = topography == RING && expect_smarts;
    let chain = topography == CHAIN && expect_smarts;
    let by_topo =
        |result: Option<String>, cph: Option<String>, ring_s: &str, chain_s: &str, other: &str| {
            if ring {
                str_append(result, cph, Some(ring_s))
            } else if chain {
                str_append(result, cph, Some(chain_s))
            } else {
                str_append(result, cph, Some(other))
            }
        };
    match btype {
        NONE => str_append(result, cph, Some(".")),
        SINGLE => {
            let plain = |result: Option<String>, cph: Option<String>| {
                if ring {
                    str_append(result, cph, Some("@;-"))
                } else if chain {
                    str_append(result, cph, Some("!@;-"))
                } else if result.is_none() {
                    cph
                } else {
                    result
                }
            };
            if bcolor == 0 {
                return plain(result, cph);
            }
            let slash = if from_first {
                if bcolor & DY_ABOVE_FIRST != 0 {
                    Some("/")
                } else if bcolor & DY_BELOW_FIRST != 0 {
                    Some("\\")
                } else if bcolor & DY_ABOVE_SECOND != 0 {
                    Some("\\")
                } else if bcolor & DY_BELOW_SECOND != 0 {
                    Some("/")
                } else {
                    None
                }
            } else if bcolor & DY_ABOVE_FIRST != 0 {
                Some("\\")
            } else if bcolor & DY_BELOW_FIRST != 0 {
                Some("/")
            } else if bcolor & DY_ABOVE_SECOND != 0 {
                Some("/")
            } else if bcolor & DY_BELOW_SECOND != 0 {
                Some("\\")
            } else {
                None
            };
            match slash {
                Some(s) => str_append(result, cph, Some(s)),
                None => plain(result, cph),
            }
        }
        DOUBLE => by_topo(result, cph, "@;=", "!@;=", "="),
        TRIPLE => str_append(result, cph, Some("#")),
        AROMATIC => match result {
            None => cph,
            Some(_) => {
                if cph
                    .as_deref()
                    .is_none_or(|c| !c.starts_with('c') && !c.starts_with('n'))
                {
                    str_append(result, cph, Some(""))
                } else {
                    str_append(result, cph, Some(":"))
                }
            }
        },
        SINGLE_DOUBLE => {
            if expect_smarts {
                by_topo(result, cph, "@;=,@;-", "!@;=,!@;-", "-,=")
            } else {
                str_append(result, cph, Some(""))
            }
        }
        SINGLE_AROMATIC => {
            if expect_smarts {
                by_topo(result, cph, "@;-,@;:", "!@;-,!@;:", "-,:")
            } else {
                str_append(result, cph, Some(""))
            }
        }
        DOUBLE_AROMATIC => {
            if expect_smarts {
                by_topo(result, cph, "@;=,@;:", "!@;=,!@;:", "=,:")
            } else {
                str_append(result, cph, Some(""))
            }
        }
        ANY_BOND => {
            if expect_smarts {
                by_topo(result, cph, "@;~", "!@;~", "~")
            } else {
                str_append(result, cph, Some(""))
            }
        }
        _ => {
            if expect_smarts {
                str_append(result, cph, Some("?"))
            } else {
                str_append(result, cph, Some(""))
            }
        }
    }
}

/// `Volume`
fn volume(t: &[[f64; 3]; 4]) -> f64 {
    let (ax, ay, az) = (t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]);
    let (bx, by, bz) = (t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]);
    let (cx, cy, cz) = (t[3][0] - t[0][0], t[3][1] - t[0][1], t[3][2] - t[0][2]);
    ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)
}

/// `PerceiveSMILESAtomStereo`
fn perceive_smiles_atom_stereo(
    mp: &Molecule,
    icurr: usize,
    lig_ind: &[i32],
    nbp: &Neighbourhood,
    use_z: bool,
) -> i32 {
    let nlig = lig_ind.len();
    if nbp.atoms.len() > 4 || nbp.atoms.len() < 3 || nlig != 4 {
        return 0;
    }
    let found = nbp.bonds.iter().any(|&bi| {
        let bp = &mp.bonds[bi];
        bp.atoms[0] == icurr as i32 + 1
            && (bp.stereo_symbol == UP
                || bp.stereo_symbol == DOWN
                || (use_z && mp.atoms[(bp.atoms[1] - 1) as usize].z != 0.0))
    });
    if !found {
        return 0;
    }
    let mut tetra = [[0.0f64; 3]; 4];
    for i in 0..nlig {
        if lig_ind[i] as usize == icurr {
            tetra[i] = [
                f64::from(mp.atoms[icurr].x),
                f64::from(mp.atoms[icurr].y),
                0.0,
            ];
        } else {
            let a = &mp.atoms[lig_ind[i] as usize];
            tetra[i][0] = f64::from(a.x);
            tetra[i][1] = f64::from(a.y);
            for j in 0..nbp.atoms.len() {
                if nbp.atoms[j] as i32 == lig_ind[i] {
                    let bp = &mp.bonds[nbp.bonds[j]];
                    tetra[i][2] = if bp.atoms[0] != icurr as i32 + 1 {
                        0.0
                    } else if bp.stereo_symbol == UP {
                        1.0
                    } else if bp.stereo_symbol == DOWN {
                        -1.0
                    } else if use_z {
                        f64::from(a.z)
                    } else {
                        0.0
                    };
                    break;
                }
            }
        }
    }
    let vol = volume(&tetra);
    if vol < -0.01 {
        -1
    } else if vol > 0.01 {
        1
    } else {
        0
    }
}

/// the character `c` shifted by `n` (C char arithmetic)
fn shift_last(buf: &mut String, n: i32) {
    if let Some(c) = buf.pop() {
        buf.push(((c as i32 + n) as u8) as char);
    }
}

fn ring_digits(i: usize) -> String {
    if i >= 9 {
        format!("%{}", i + 1)
    } else {
        format!("{}", i + 1)
    }
}

/// `SmilesAtom`
fn smiles_atom(
    mp: &Molecule,
    index: usize,
    lig_ind: &[i32],
    nbp: &Neighbourhood,
    h_count: &[i32],
    expect_smarts: bool,
    aromatic_smiles: bool,
    use_z: bool,
) -> String {
    let nlig = lig_ind.len();
    let mut atom_text: Option<&str> = None;
    for (k, line) in mp.prop_lines.iter().enumerate() {
        if line.starts_with("A  ") {
            let n = super::scan::Scan::new(line.get(3..).unwrap_or("").as_bytes())
                .int(0)
                .unwrap_or(0);
            if n == index as i32 + 1 {
                atom_text = mp.prop_lines.get(k + 1).map(String::as_str);
            }
        }
    }
    let ap = &mp.atoms[index];
    let sym = ap.atom_symbol.as_str();
    let stereo = perceive_smiles_atom_stereo(mp, index, lig_ind, nbp, use_z);
    let btype = |k: usize| mp.bonds[nbp.bonds[k]].bond_type;

    if atom_symbol_match(sym, "A,B,C,N,O,P,S,F,Cl,Br,I,R")
        && ap.value == 0.0
        && ap.charge == NONE
        && ap.radical == NONE
        && ap.mass_difference == NONE
        && ap.query_h_count == NONE
        && (ap.mapping == NONE || sym == "R")
        && (ap.sub_desc == NONE || !expect_smarts)
        && atom_text.is_none()
        && !(atom_symbol_match(sym, "F,Cl,Br,I") && nlig > 1)
        && stereo == 0
    {
        let mut result = if sym.starts_with('R') || sym.starts_with('A') {
            if sym == "R" && !ap.atext.is_empty() {
                format!("{{{}}}", ap.atext)
            } else {
                "*".to_string()
            }
        } else {
            sym.to_string()
        };
        if atom_symbol_match(sym, "C,N,O,P,S,B") {
            let mut aromatic_level = 0;
            let mut exo_methyl = false;
            for i in 0..nbp.atoms.len() {
                let b = &mp.bonds[nbp.bonds[i]];
                if b.bond_type == SINGLE {
                    exo_methyl = true;
                }
                if b.bond_type == AROMATIC {
                    aromatic_level |= 1;
                } else if matches!(
                    b.bond_type,
                    SINGLE_DOUBLE | SINGLE_AROMATIC | DOUBLE_AROMATIC | ANY_BOND
                ) {
                    aromatic_level |= 2;
                }
                if b.topography == CHAIN {
                    aromatic_level |= 4;
                }
            }
            if nbp.atoms.len() != 1 || !atom_symbol_match(sym, "C") || !expect_smarts {
                exo_methyl = false;
            }
            let lower = |s: &str| -> String {
                let mut b = s.as_bytes().to_vec();
                b[0] += b'c' - b'C';
                String::from_utf8(b).expect("ascii")
            };
            if aromatic_level == 1 {
                result = if expect_smarts || aromatic_smiles {
                    lower(&result)
                } else {
                    "*".to_string()
                };
            } else if aromatic_level & 1 != 0 && (expect_smarts || aromatic_smiles) {
                result = lower(&result);
            } else if (aromatic_level == 5 || aromatic_level & 2 != 0)
                && !(expect_smarts || aromatic_smiles)
            {
                result = "*".to_string();
            } else if aromatic_level & 2 != 0 || exo_methyl {
                let c = sym.as_bytes()[0];
                result = format!("[{},{}]", c as char, (c + (b'c' - b'C')) as char);
            }
        }
        return result;
    } else if let Some(text) = atom_text {
        if expect_smarts {
            let mut r = String::new();
            if !text.starts_with('[') {
                r.push('[');
            }
            r.push_str(text);
            if !text.starts_with('[') {
                r.push(']');
            }
            return r;
        }
        return if text == "C" { "C".into() } else { "*".into() };
    } else if sym.starts_with("R#") {
        let rgnum = get_num_property(&mp.prop_lines, "M  RGP", index as i32 + 1).unwrap_or(
            if ap.mass_difference > 0 {
                ap.mass_difference
            } else {
                0
            },
        );
        return if rgnum == 1 {
            "[*+]".into()
        } else {
            format!("[*+{rgnum}]")
        };
    } else if sym.len() >= 2 && sym.starts_with('R') && (b'1'..=b'9').contains(&sym.as_bytes()[1]) {
        let rgnum = i32::from(sym.as_bytes()[1] - b'0');
        return if rgnum == 1 {
            "[*+]".into()
        } else {
            format!("[*+{rgnum}]")
        };
    } else if sym == "R" {
        let mut rgnum = ap.value as i32;
        if ap.mass_difference > 0 && rgnum == 0 {
            rgnum = ap.mass_difference;
        }
        if ap.charge > 0 && rgnum == 0 {
            rgnum = ap.charge;
        }
        return if rgnum == 1 {
            "[*+]".into()
        } else {
            format!("[*+{rgnum}]")
        };
    } else if sym == "L" {
        if !expect_smarts {
            return "*".into();
        }
        let mut buffer = String::from("[");
        if let Some(slp) = mp.symbol_lists.iter().find(|l| l.atom == index as i32 + 1) {
            let s = slp.string.as_bytes();
            let mut start_of_symbol = true;
            let mut p = 0;
            while p < s.len() {
                let c = s[p];
                if c == b' ' || c == b',' {
                    start_of_symbol = true;
                    p += 1;
                    continue;
                }
                if start_of_symbol && p != 0 {
                    buffer.push(if slp.logic != 0 { ',' } else { ';' });
                }
                if slp.logic == 0 && start_of_symbol {
                    buffer.push('!');
                }
                if start_of_symbol {
                    start_of_symbol = false;
                    let next_alpha = s.get(p + 1).is_some_and(u8::is_ascii_alphabetic);
                    let num = match c {
                        b'C' => Some("#6"),
                        b'N' => Some("#7"),
                        b'O' => Some("#8"),
                        b'S' => Some("#16"),
                        _ => None,
                    };
                    if let (Some(n), false) = (num, next_alpha) {
                        buffer.push_str(n);
                        p += 1;
                        continue;
                    }
                }
                buffer.push(c as char);
                p += 1;
            }
        } else {
            buffer.push('*');
        }
        push_sub_desc(&mut buffer, ap.sub_desc, nbp.atoms.len());
        buffer.push(']');
        return buffer;
    } else if sym == "Q" || sym == "A" {
        return if expect_smarts {
            if sym == "Q" {
                "[!#6;!#1]".into()
            } else {
                "[!#1]".into()
            }
        } else {
            "*".into()
        };
    }

    let mut buffer = String::from("[");
    if ap.mass_difference != NONE {
        match PTABLE.iter().find(|e| e.0 == sym) {
            None => buffer.push_str(&ap.mass_difference.to_string()),
            Some(e) => {
                buffer.push_str(&(((f64::from(e.2) + 0.5) as i32) + ap.mass_difference).to_string())
            }
        }
    } else if sym == "D" {
        buffer.push('2');
    } else if sym == "T" {
        buffer.push('3');
    }
    if sym == "D" || sym == "T" {
        buffer.push('H');
    } else if sym == "A" {
        buffer.push('*');
    } else {
        buffer.push_str(sym);
    }
    if atom_symbol_match(sym, "C,N,O,P,S,B") {
        let mut aromatic_level = 0;
        for i in 0..nbp.atoms.len() {
            match btype(i) {
                AROMATIC => aromatic_level |= 1,
                SINGLE_AROMATIC | DOUBLE_AROMATIC | ANY_BOND => aromatic_level |= 2,
                _ => {}
            }
        }
        if aromatic_level == 1 {
            shift_last(&mut buffer, i32::from(b'c' - b'C'));
        } else if aromatic_level & 2 != 0 {
            buffer.pop();
            let c = sym.as_bytes()[0];
            buffer.push_str(&format!("{},{}", c as char, (c + (b'c' - b'C')) as char));
        }
    }
    if stereo == 1 {
        buffer.push_str("@@");
    } else if stereo == -1 {
        buffer.push('@');
    }
    if h_count[index + 1] != 0
        && !expect_smarts
        && (ap.sub_desc == NONE || !expect_smarts)
        && ap.query_h_count == 0
    {
        if h_count[index + 1] == 1 {
            buffer.push('H');
        } else {
            buffer.push_str("H0");
            shift_last(&mut buffer, h_count[index + 1]);
        }
    } else if ap.query_h_count != 0 {
        if ap.query_h_count > NONE + 1 && h_count[index + 1] == 0 && expect_smarts {
            buffer.push(';');
            let mut i = ap.query_h_count - 1;
            loop {
                buffer.push_str("*H0,");
                let l = buffer.len();
                let mut b = buffer.clone().into_bytes();
                b[l - 2] = (b[l - 2] as i32 + i) as u8;
                buffer = String::from_utf8(b).expect("ascii");
                if i >= 4 {
                    break;
                }
                i += 1;
            }
            buffer.pop();
        } else if !expect_smarts {
            if ap.query_h_count == 2 {
                buffer.push('H');
            } else if ap.query_h_count == 1 {
                // nothing
            } else {
                buffer.push_str("H0");
                shift_last(&mut buffer, ap.query_h_count - 1);
            }
        }
    }
    if ap.charge > 0 {
        if buffer.contains(';') {
            buffer.push_str(";*");
        }
        if ap.charge == 1 {
            buffer.push('+');
        } else if ap.charge <= 9 {
            buffer.push_str("+0");
            shift_last(&mut buffer, ap.charge);
        } else {
            buffer.push_str("+00");
            let l = buffer.len();
            let mut b = buffer.into_bytes();
            b[l - 2] = (b[l - 2] as i32 + (ap.charge / 10) % 10) as u8;
            b[l - 1] = (b[l - 1] as i32 + ap.charge % 10) as u8;
            buffer = String::from_utf8(b).expect("ascii");
        }
    } else if ap.charge < 0 {
        if buffer.contains(';') {
            buffer.push_str(";*");
        }
        if ap.charge == -1 {
            buffer.push('-');
        } else if ap.charge >= -9 {
            buffer.push_str("-0");
            shift_last(&mut buffer, -ap.charge);
        } else {
            buffer.push_str("-00");
            let l = buffer.len();
            let mut b = buffer.into_bytes();
            b[l - 2] = (b[l - 2] as i32 + (-ap.charge / 10) % 10) as u8;
            b[l - 1] = (b[l - 1] as i32 + (-ap.charge) % 10) as u8;
            buffer = String::from_utf8(b).expect("ascii");
        }
    }
    if ap.mapping > 0 {
        buffer.push(':');
        buffer.push_str(&(ap.mapping - 1).to_string());
    }
    if expect_smarts {
        push_sub_desc(&mut buffer, ap.sub_desc, nbp.atoms.len());
    }
    buffer.push(']');
    buffer
}

/// the `;D<n>` substitution descriptors of query atoms
fn push_sub_desc(buffer: &mut String, sub_desc: i32, n_ligands: usize) {
    if sub_desc != NONE {
        if sub_desc > 0 {
            buffer.push_str(&format!(";D{sub_desc}"));
        } else if sub_desc == -1 {
            buffer.push_str(";D0");
        } else if sub_desc == -2 {
            buffer.push_str(&format!(";D{n_ligands}"));
        }
    }
}

struct Ctx<'a> {
    open_branches: [i32; MAXOPEN],
    open_neighbours: [i32; MAXOPEN],
    h_count: &'a [i32],
    nsmi_atoms: i32,
    expect_smarts: bool,
    aromatic_smiles: bool,
    use_z: bool,
}

/// `SmilesBranch`
fn smiles_branch(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    ctx: &mut Ctx,
    mut istart: usize,
    mut iprev: i32,
) -> Option<String> {
    let mut result: Option<String> = None;
    loop {
        let mut ligands: Option<String>;
        let mut closures: Option<String> = None;
        let mut lig_ind: Vec<i32> = Vec::new();
        if iprev >= 0 {
            lig_ind.push(iprev);
        }
        if nbp[istart].atoms.len() == 3
            && atom_symbol_match(&mp.atoms[istart].atom_symbol, "C,Si,Ge")
        {
            lig_ind.push(istart as i32);
        }
        ctx.nsmi_atoms += 1;
        mp.atoms[istart].color = ctx.nsmi_atoms;

        let mut nbranch = 0;
        let mut jbranch = 0;
        for j in 0..nbp[istart].atoms.len() {
            let aj = nbp[istart].atoms[j];
            if aj as i32 == iprev {
                continue;
            }
            let bj = nbp[istart].bonds[j];
            if mp.atoms[aj].color != NONE {
                // ring closure
                let b = &mp.bonds[bj];
                let (bt, ff, topo, col) = (
                    b.bond_type,
                    b.atoms[0] == istart as i32 + 1,
                    b.topography,
                    b.color,
                );
                if let Some(i) = (0..MAXOPEN).find(|&i| {
                    ctx.open_branches[i] == istart as i32 && ctx.open_neighbours[i] == aj as i32
                }) {
                    if lig_ind.len() < 4 {
                        lig_ind.push(aj as i32);
                    }
                    closures = bnd_append(closures, None, bt, ff, topo, col, ctx.expect_smarts);
                    closures = str_append(closures, None, Some(&ring_digits(i)));
                    ctx.open_branches[i] = -2;
                    ctx.open_neighbours[i] = -1;
                    continue;
                }
                if !(0..MAXOPEN).any(|i| ctx.open_branches[i] >= 0 || ctx.open_neighbours[i] >= 0) {
                    ctx.open_branches.fill(-1);
                    ctx.open_neighbours.fill(-1);
                }
                if let Some(i) = (0..MAXOPEN).find(|&i| ctx.open_branches[i] == -1) {
                    if lig_ind.len() < 4 {
                        lig_ind.push(aj as i32);
                    }
                    closures = bnd_append(closures, None, bt, ff, topo, col, ctx.expect_smarts);
                    closures = str_append(closures, None, Some(&ring_digits(i)));
                    ctx.open_branches[i] = aj as i32;
                    ctx.open_neighbours[i] = istart as i32;
                }
            } else {
                jbranch = j;
                nbranch += 1;
            }
        }

        let bond_of = |mp: &Molecule, k: usize| {
            let b = &mp.bonds[nbp[istart].bonds[k]];
            (
                b.bond_type,
                b.atoms[0] == istart as i32 + 1,
                b.topography,
                b.color,
            )
        };
        if nbranch == 1 {
            lig_ind.push(nbp[istart].atoms[jbranch] as i32);
            let symbol = smiles_atom(
                mp,
                istart,
                &lig_ind,
                &nbp[istart],
                ctx.h_count,
                ctx.expect_smarts,
                ctx.aromatic_smiles,
                ctx.use_z,
            );
            let (bt, ff, topo, col) = bond_of(mp, jbranch);
            ligands = bnd_append(None, None, bt, ff, topo, col, ctx.expect_smarts);
            result = str_append(result, Some(symbol), Some(""));
            result = str_append(result, closures, Some(""));
            result = str_append(result, ligands, Some(""));
            iprev = istart as i32;
            istart = nbp[istart].atoms[jbranch];
            continue;
        }

        let mut cp: Option<String> = None;
        ligands = None;
        for j in 0..nbp[istart].atoms.len() {
            let aj = nbp[istart].atoms[j];
            if aj as i32 == iprev {
                continue;
            }
            if mp.atoms[aj].color == NONE {
                let cph = smiles_branch(mp, nbp, ctx, aj, istart as i32);
                let (bt, ff, topo, col) = bond_of(mp, j);
                let cph = bnd_append(None, cph, bt, ff, topo, col, ctx.expect_smarts);
                if cp.is_some() {
                    ligands = str_append(ligands, None, Some("("));
                    ligands = str_append(ligands, cp, Some(""));
                    ligands = str_append(ligands, None, Some(")"));
                }
                cp = cph;
            }
        }
        if cp.is_some() {
            ligands = str_append(ligands, cp, Some(""));
        }
        for i in 0..MAXOPEN {
            if ctx.open_branches[i] == istart as i32 {
                if lig_ind.len() < 4 {
                    lig_ind.push(ctx.open_neighbours[i]);
                }
                closures = str_append(closures, None, Some(&ring_digits(i)));
                ctx.open_branches[i] = -2;
                ctx.open_neighbours[i] = -1;
            }
        }
        if !(0..MAXOPEN).any(|i| ctx.open_branches[i] >= 0 || ctx.open_neighbours[i] >= 0) {
            ctx.open_branches.fill(-1);
            ctx.open_neighbours.fill(-1);
        }
        for i in 0..nbp[istart].atoms.len() {
            if !lig_ind.contains(&(nbp[istart].atoms[i] as i32)) && lig_ind.len() < 4 {
                lig_ind.push(nbp[istart].atoms[i] as i32);
            }
        }
        let symbol = smiles_atom(
            mp,
            istart,
            &lig_ind,
            &nbp[istart],
            ctx.h_count,
            ctx.expect_smarts,
            ctx.aromatic_smiles,
            ctx.use_z,
        );
        result = str_append(result, Some(symbol), None);
        result = str_append(result, closures, None);
        result = str_append(result, ligands, None);
        return result;
    }
}

/// `DepthFirstVisit`
fn depth_first_visit(
    root: usize,
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    result: &mut Vec<usize>,
) {
    for i in 0..nbp[root].atoms.len() {
        let (ai, bi) = (nbp[root].atoms[i], nbp[root].bonds[i]);
        if mp.bonds[bi].color == 0 {
            mp.atoms[ai].color = VISITED;
            mp.bonds[bi].color = VISITED;
            result.push(bi);
            depth_first_visit(ai, mp, nbp, result);
        }
    }
}

/// `DepthFirstBondList`
fn depth_first_bond_list(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    numbering: Option<&[i32]>,
) -> Vec<usize> {
    let mut result = Vec::new();
    for a in &mut mp.atoms {
        a.color = 0;
    }
    for b in &mut mp.bonds {
        b.color = 0;
    }
    let n = mp.atoms.len();
    loop {
        let ibranch = match numbering {
            None => mp.atoms.iter().position(|a| a.color == 0).unwrap_or(n),
            Some(num) => {
                let mut ib = None;
                let mut best = n as i32 + 1;
                for i in 0..n {
                    if mp.atoms[i].color == 0 && num[i] < best {
                        ib = Some(i);
                        best = num[i];
                    }
                }
                ib.unwrap_or(n)
            }
        };
        if ibranch >= n {
            break;
        }
        mp.atoms[ibranch].color = VISITED;
        depth_first_visit(ibranch, mp, nbp, &mut result);
    }
    for a in &mut mp.atoms {
        a.color = 0;
    }
    for b in &mut mp.bonds {
        b.color = 0;
    }
    result
}

/// `PerceiveSmallRingBonds`
fn perceive_small_ring_bonds(mp: &mut Molecule) {
    let graph: Vec<[u32; 2]> = mp
        .bonds
        .iter()
        .map(|b| [b.atoms[0] as u32, b.atoms[1] as u32])
        .collect();
    let mut rings = ring_list(&graph);
    combine_rings(&mut rings);
    for r in &rings {
        for (i, b) in mp.bonds.iter_mut().enumerate() {
            if r.bond_set.is_member(i) && r.cardinality < 8 {
                b.topography = RING;
            }
        }
    }
}

/// `PerceiveDBStereo`: the `/` and `\` directions of the single bonds around double bonds end up in
/// the bond colours
fn perceive_db_stereo(mp: &mut Molecule, nbp: &[Neighbourhood], numbering: Option<&[i32]>) {
    if mp.bonds.is_empty() {
        return;
    }
    let topography_sav: Vec<i32> = mp.bonds.iter().map(|b| b.topography).collect();
    let df = depth_first_bond_list(mp, nbp, numbering);
    perceive_small_ring_bonds(mp);
    for &bi in &df {
        let bp = mp.bonds[bi].clone();
        if bp.bond_type != DOUBLE || bp.topography == RING || bp.stereo_symbol == CIS_TRANS_EITHER {
            continue;
        }
        let (ai1, ai2) = ((bp.atoms[0] - 1) as usize, (bp.atoms[1] - 1) as usize);
        if !(2..=3).contains(&nbp[ai1].atoms.len()) || !(2..=3).contains(&nbp[ai2].atoms.len()) {
            continue;
        }
        if !atom_symbol_match(&mp.atoms[ai1].atom_symbol, "C,N")
            || !atom_symbol_match(&mp.atoms[ai2].atom_symbol, "C,N")
        {
            continue;
        }
        let nmulti = |mp: &Molecule, ai: usize| {
            nbp[ai]
                .bonds
                .iter()
                .filter(|&&b| mp.bonds[b].bond_type != SINGLE)
                .count()
        };
        if nmulti(mp, ai1) > 1 || nmulti(mp, ai2) > 1 {
            continue;
        }
        let pos = |mp: &Molecule, k: usize| (f64::from(mp.atoms[k].x), f64::from(mp.atoms[k].y));
        let (x1, y1) = pos(mp, ai1);
        let (x2, y2) = pos(mp, ai2);
        let (mut dx1, mut dy1) = (x2 - x1, y2 - y1);
        let swap_needed = numbering.is_some_and(|n| n[ai1] > n[ai2]);
        let _ = swap_needed; // SwapFirstSecond() is the identity
        let mut collision = false;
        for pass in 0..2 {
            if pass == 1 {
                if !collision {
                    break;
                }
                collision = false;
                dx1 = -dx1;
                dy1 = -dy1;
            }
            for (ai, other_atom_no) in [(ai1, bp.atoms[1]), (ai2, bp.atoms[0])] {
                for j in 0..nbp[ai].atoms.len() {
                    if nbp[ai].atoms[j] as i32 + 1 == other_atom_no {
                        continue;
                    }
                    let (ax, ay) = pos(mp, ai);
                    let (lx, ly) = pos(mp, nbp[ai].atoms[j]);
                    let (dx2, dy2) = (lx - ax, ly - ay);
                    let above = dx1 * dy2 - dx2 * dy1 > 0.0;
                    let b = &mut mp.bonds[nbp[ai].bonds[j]];
                    let first = ai as i32 + 1 == b.atoms[0];
                    if pass == 0 {
                        b.color |= match (above, first) {
                            (true, true) => DY_ABOVE_FIRST,
                            (true, false) => DY_ABOVE_SECOND,
                            (false, true) => DY_BELOW_FIRST,
                            (false, false) => DY_BELOW_SECOND,
                        };
                    } else {
                        match (above, first) {
                            (true, true) => {
                                b.color &= !DY_BELOW_FIRST;
                                b.color |= DY_ABOVE_FIRST;
                            }
                            (true, false) => {
                                b.color &= !DY_BELOW_SECOND;
                                b.color |= DY_ABOVE_SECOND;
                            }
                            (false, true) => {
                                b.color &= !DY_ABOVE_FIRST;
                                b.color |= DY_BELOW_FIRST;
                            }
                            (false, false) => {
                                b.color &= !DY_ABOVE_SECOND;
                                b.color |= DY_BELOW_SECOND;
                            }
                        }
                    }
                    if b.color == DY_ABOVE_BOTH || b.color == DY_BELOW_BOTH {
                        collision = true;
                    }
                }
            }
        }
        if bp.stereo_symbol == CIS_TRANS_SWAPPED {
            for j in 0..nbp[ai2].atoms.len() {
                if nbp[ai2].atoms[j] as i32 + 1 == bp.atoms[0] {
                    continue;
                }
                let b = &mut mp.bonds[nbp[ai2].bonds[j]];
                b.color = match b.color {
                    DY_ABOVE_FIRST => DY_BELOW_FIRST,
                    DY_BELOW_FIRST => DY_ABOVE_FIRST,
                    DY_ABOVE_SECOND => DY_BELOW_SECOND,
                    DY_BELOW_SECOND => DY_ABOVE_SECOND,
                    c => c,
                };
            }
        }
    }
    for (b, t) in mp.bonds.iter_mut().zip(topography_sav) {
        b.topography = t;
    }
}

/// `SortNeighboursByShortcut`
fn sort_neighbours_by_shortcut(mp: &Molecule, nbp: &mut [Neighbourhood]) {
    let is_shortcut = |a: usize| mp.atoms[a].atom_symbol == "R" && !mp.atoms[a].atext.is_empty();
    for nb in nbp.iter_mut() {
        for i in 1..nb.atoms.len() {
            for j in (0..i).rev() {
                let (s1, s2) = (is_shortcut(nb.atoms[j]), is_shortcut(nb.atoms[j + 1]));
                let (t1, t2) = (
                    mp.bonds[nb.bonds[j]].bond_type,
                    mp.bonds[nb.bonds[j + 1]].bond_type,
                );
                if (s2 && !s1)
                    || (s1 && s2 && t1 < t2)
                    || (s1 && s2 && t1 == t2 && nb.atoms[j] > nb.atoms[j + 1])
                {
                    nb.bonds.swap(j, j + 1);
                    nb.atoms.swap(j, j + 1);
                } else {
                    break;
                }
            }
        }
    }
}

/// `MOLToSMIExt` without the coordinate output; `Err` for inputs this port does not cover yet
pub fn mol_to_smi_ext(
    mp: &mut Molecule,
    flags: i32,
    numbering: Option<&[i32]>,
) -> Result<Option<String>, &'static str> {
    if mp.atoms.is_empty() {
        return Ok(None);
    }
    if flags & SMARTS_PERCEPTION != 0 {
        return Err("SMARTS perception (PerceiveMarkush) is not ported yet");
    }
    if flags & DY_AROMATICITY != 0 {
        make_hydrogens_implicit(mp);
        for a in &mut mp.atoms {
            a.query_h_count = NONE;
        }
    }
    let n = mp.atoms.len();
    let mut h_count = vec![0; n + 1];
    compute_implicit_h(mp, &mut h_count);
    reset_colors(mp);
    let mut nbp = setup_neighbourhood(mp).ok_or("too many neighbours")?;
    if flags & DY_AROMATICITY != 0 {
        perceive_dy_aromaticity(mp, &nbp);
    }
    if let Some(num) = numbering {
        sort_bonds_with_numbering(mp, num);
        nbp = setup_neighbourhood(mp).ok_or("too many neighbours")?;
        sort_neighbours_by_numbering(&mut nbp, num);
    } else {
        sort_neighbours_by_shortcut(mp, &mut nbp);
    }
    perceive_db_stereo(mp, &nbp, numbering);

    let mut ctx = Ctx {
        open_branches: [-1; MAXOPEN],
        // the C code leaves this array uninitialised; zeros are what the reference builds see
        open_neighbours: [0; MAXOPEN],
        h_count: &h_count,
        nsmi_atoms: 0,
        expect_smarts: flags & SMARTS_PERCEPTION != 0,
        aromatic_smiles: flags & (AROMATIC_SMILES | DY_AROMATICITY) != 0,
        use_z: flags & USE_Z != 0,
    };
    let mut result: Option<String> = None;
    if numbering.is_none() {
        // strings of shortcuts first
        while let Some(mut ibranch) = (0..n).find(|&i| {
            mp.atoms[i].atom_symbol == "R"
                && !mp.atoms[i].atext.is_empty()
                && mp.atoms[i].color == 0
        }) {
            if result.is_some() {
                result = str_append(result, None, Some("."));
            }
            result = str_append(result, None, Some("{"));
            result = str_append(result, None, Some(&mp.atoms[ibranch].atext.clone()));
            result = str_append(result, None, Some("}"));
            ctx.nsmi_atoms += 1;
            mp.atoms[ibranch].color = ctx.nsmi_atoms;
            loop {
                let mut first_shortcut: Option<usize> = None;
                let mut shortcut_bond_type = NONE;
                let mut closures: Option<String> = None;
                for j in 0..nbp[ibranch].atoms.len() {
                    let aj = nbp[ibranch].atoms[j];
                    if mp.atoms[aj].color != NONE {
                        continue;
                    }
                    if first_shortcut.is_none()
                        && mp.atoms[aj].atom_symbol == "R"
                        && !mp.atoms[aj].atext.is_empty()
                    {
                        first_shortcut = Some(aj);
                        shortcut_bond_type = mp.bonds[nbp[ibranch].bonds[j]].bond_type;
                    } else {
                        let i = (0..MAXOPEN)
                            .find(|&i| ctx.open_branches[i] == -1)
                            .unwrap_or(MAXOPEN);
                        closures = str_append(closures, None, Some(&ring_digits(i)));
                        if i < MAXOPEN {
                            ctx.open_branches[i] = aj as i32;
                            ctx.open_neighbours[i] = ibranch as i32;
                        }
                    }
                }
                result = str_append(result, closures, Some(""));
                let Some(fs) = first_shortcut else { break };
                if shortcut_bond_type == NONE {
                    result = str_append(result, None, Some("."));
                }
                ibranch = fs;
                result = str_append(result, None, Some("{"));
                result = str_append(result, None, Some(&mp.atoms[ibranch].atext.clone()));
                result = str_append(result, None, Some("}"));
                for j in 0..MAXOPEN {
                    if ctx.open_branches[j] == -1 {
                        continue;
                    }
                    if ibranch as i32 == ctx.open_branches[j] {
                        result = str_append(result, None, Some(&format!("{}", j + 1)));
                        ctx.open_branches[j] = -1;
                        ctx.open_neighbours[j] = -1;
                    }
                }
                ctx.nsmi_atoms += 1;
                mp.atoms[ibranch].color = ctx.nsmi_atoms;
            }
        }
    }
    loop {
        let ibranch = match numbering {
            None => {
                let mut ib = n;
                let mut best_degree = 100;
                for i in 0..n {
                    if mp.atoms[i].color == 0 && (nbp[i].atoms.len() as i32) < best_degree {
                        ib = i;
                        best_degree = nbp[i].atoms.len() as i32;
                    }
                }
                ib
            }
            Some(num) => {
                let mut ib: Option<usize> = None;
                let mut best = n as i32 + 1;
                for i in 0..n {
                    if mp.atoms[i].color == 0
                        && num[i] < best
                        && mp.atoms[i].atom_symbol.starts_with('R')
                    {
                        ib = Some(i);
                        best = num[i];
                    }
                }
                if flags & TO_MONOMER != 0 {
                    best = n as i32 + 1;
                }
                for i in 0..n {
                    if mp.atoms[i].color == 0
                        && num[i] < best
                        && !mp.atoms[i].atom_symbol.starts_with('R')
                    {
                        ib = Some(i);
                        best = num[i];
                    }
                }
                ib.unwrap_or(n)
            }
        };
        if ibranch == n {
            break;
        }
        for i in 0..MAXOPEN {
            if ctx.open_branches[i] < 0 {
                ctx.open_branches[i] = -1;
            }
        }
        let cp = smiles_branch(mp, &nbp, &mut ctx, ibranch, -1);
        result = if result.is_some() {
            str_append(result, cp, Some("."))
        } else {
            str_append(result, cp, Some(""))
        };
    }
    Ok(result)
}

/// `MOLToSMI`
pub fn mol_to_smi(mp: &mut Molecule, flags: i32) -> Result<Option<String>, &'static str> {
    mol_to_smi_ext(mp, flags, None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::reaccsio::mol_str_to_mol;

    /// MOLToSMI must write what the Avalon C code writes
    #[test]
    fn mol_to_smi_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_smiles.jsonl");
        let mut bad = Vec::new();
        let mut n = 0;
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().expect("input");
            let flags = v["flags"].as_i64().expect("flags") as i32;
            let want = v["output"].as_str().expect("output").trim_end_matches('\n');
            let mut m = mol_str_to_mol(input).expect("molecule");
            let got = mol_to_smi(&mut m, flags)
                .expect("supported")
                .unwrap_or_else(|| "NULL".into());
            n += 1;
            if got != want {
                bad.push(format!("flags {flags}\n{input}\nwant {want}\ngot  {got}"));
            }
        }
        assert!(
            bad.is_empty(),
            "{} of {n} differ; first:\n{}",
            bad.len(),
            bad[0]
        );
    }
}
