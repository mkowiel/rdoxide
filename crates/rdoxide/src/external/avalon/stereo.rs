//! Port of Avalon Toolkit's `stereo.c`: stereo parity perception from wedge bonds and the
//! checks built on it. Messages go to a `Vec<String>` instead of the C message list; the
//! `stereo_error` global is returned through a `&mut &str`. `IsStereoMatch` needs the
//! substructure matcher and `AllCentersRefined` reads one element past its sort buffer in C,
//! so neither is ported yet.

use super::perceive::{DOUBLE, Neighbourhood, SINGLE, setup_neighbourhood};
use super::reaccs::{CIS_TRANS_EITHER, Molecule, NONE};
use super::scan::fmt_gp;
use super::utilities::{DOWN, EITHER, UP, atom_symbol_match};

pub const CIS: i32 = 1;
pub const TRANS: i32 = 2;
pub const ALLENE_PARITY: i32 = -2;
pub const ILLEGAL_REPRESENTATION: i32 = -1;
pub const UNDEFINED_PARITY: i32 = 0;
pub const ODD_PARITY: i32 = 1;
pub const EVEN_PARITY: i32 = 2;

pub const EITHER_BOND_FOUND: i32 = 1;
pub const STEREO_BOND_AT_NON_STEREO_ATOM: i32 = 2;
pub const ZEROED_Z_COORDINATES: i32 = 4;
pub const CONVERTED_TO_2D: i32 = 8;

#[allow(clippy::approx_constant)] // the C code's 12-digit value is part of the behaviour
const PI: f64 = 3.14159265359;
const ANGLE_EPSILON: f64 = PI / 180.0 * 5.0;
const EPS: f64 = 0.0000001;
const MAXMSG: usize = 255;

/// `INVERT_PARITY`
pub fn invert_parity(p: i32) -> i32 {
    if p == 0 { 0 } else { 3 - p }
}

/// `AddMsgToList`: messages are kept in the order they were added
pub fn add_msg(msgs: &mut Vec<String>, msg: String) {
    let mut m = msg;
    if m.len() > MAXMSG {
        let mut cut = MAXMSG;
        while !m.is_char_boundary(cut) {
            cut -= 1;
        }
        m.truncate(cut);
    }
    msgs.push(m);
}

/// `Angle` of `local.c`
pub fn angle(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    let l1 = (x1 * x1 + y1 * y1).sqrt();
    let l2 = (x2 * x2 + y2 * y2).sqrt();
    if l1 < 0.00001 || l2 < 0.00001 {
        return 0.0;
    }
    // safeguard against round-off errors (a NaN passes through, as in C)
    let cos_alpha = ((x1 * x2 + y1 * y2) / (l1 * l2)).clamp(-1.0, 1.0);
    let sin_alpha = (x1 * y2 - x2 * y1) / (l1 * l2);
    let mut result = cos_alpha.acos();
    if sin_alpha < 0.0 {
        result = 2.0 * PI - result;
    }
    result
}

#[derive(Clone, Copy, Default)]
struct NPoint {
    x: f64,
    y: f64,
    z: f64,
    number: i32,
}

/// `Volume`: signed volume of the tetrahedron
fn volume(t: &[NPoint; 4]) -> f64 {
    let (ax, ay, az) = (t[1].x - t[0].x, t[1].y - t[0].y, t[1].z - t[0].z);
    let (bx, by, bz) = (t[2].x - t[0].x, t[2].y - t[0].y, t[2].z - t[0].z);
    let (cx, cy, cz) = (t[3].x - t[0].x, t[3].y - t[0].y, t[3].z - t[0].z);
    ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)
}

#[derive(Clone, Copy, Default)]
struct StereoBond {
    x: f64,
    y: f64,
    symbol: i32,
    number: i32,
    angle: f64,
}

fn sort_tetra(t: &mut [NPoint; 4]) {
    for i in 1..4 {
        for j in (1..=i).rev() {
            if t[j].number < t[j - 1].number {
                t.swap(j, j - 1);
            } else {
                break;
            }
        }
    }
}

fn parity_of_volume(t: &[NPoint; 4]) -> i32 {
    if volume(t) > 0.0 {
        EVEN_PARITY
    } else {
        ODD_PARITY
    }
}

fn atom3_parity(l: &[StereoBond], err: &mut &'static str) -> i32 {
    let maxnum = l[..3].iter().map(|b| b.number).max().unwrap();
    let mut reference: Option<usize> = None;
    for (i, b) in l[..3].iter().enumerate() {
        if b.symbol != NONE {
            if reference.is_none() {
                reference = Some(i);
            } else {
                *err = "three attachments with more than 2 stereobonds";
                return ILLEGAL_REPRESENTATION;
            }
        }
    }
    let Some(reference) = reference else {
        return UNDEFINED_PARITY;
    };
    let (a, b) = match reference {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };
    let ang = angle(l[a].x, l[a].y, l[b].x, l[b].y);
    if ang < ANGLE_EPSILON || (PI - ang).abs() < ANGLE_EPSILON {
        *err = "three attachments: colinearity violation";
        return ILLEGAL_REPRESENTATION;
    }
    let mut t = [NPoint::default(); 4];
    t[0].number = maxnum + 1;
    for i in 0..3 {
        t[i + 1].x = l[i].x;
        t[i + 1].y = l[i].y;
        t[i + 1].z = if l[i].symbol == UP {
            1.0
        } else if l[i].symbol == DOWN {
            -1.0
        } else if l[i].symbol == NONE {
            0.0
        } else {
            *err = "three attachments: illegal bond symbol";
            return ILLEGAL_REPRESENTATION;
        };
        t[i + 1].number = l[i].number;
    }
    sort_tetra(&mut t);
    parity_of_volume(&t)
}

fn atom4_parity(l: &[StereoBond], err: &mut &'static str) -> i32 {
    let mut t = [NPoint::default(); 4];
    let (mut nup, mut ndown) = (0, 0);
    for i in 0..4 {
        t[i].x = l[i].x;
        t[i].y = l[i].y;
        t[i].z = 0.0;
        t[i].number = l[i].number;
        if l[i].symbol == UP {
            nup += 1;
            t[i].z = 1.0;
        } else if l[i].symbol == DOWN {
            ndown += 1;
            t[i].z = -1.0;
        } else if l[i].symbol != NONE {
            *err = "illegal bond symbol";
            return ILLEGAL_REPRESENTATION;
        }
    }
    if nup == 0 && ndown == 0 {
        return UNDEFINED_PARITY;
    }
    if nup > 2 || ndown > 2 {
        *err = "too many stereobonds";
        return ILLEGAL_REPRESENTATION;
    }
    if nup + ndown == 1 {
        let i = (0..4)
            .find(|&i| l[i].symbol == UP || l[i].symbol == DOWN)
            .unwrap_or(4);
        let mut nopposite = 0;
        for j in 0..4 {
            if i == j {
                continue;
            }
            if l[i].x * l[j].x + l[i].y * l[j].y < 0.0 {
                nopposite += 1;
            }
        }
        if nopposite > 2 {
            *err = "UMBRELLA: all non-stereo bonds opposite to single stereo bond";
            return ILLEGAL_REPRESENTATION;
        }
    }
    for i in 0..2 {
        if (l[i].symbol == UP && l[i + 2].symbol == DOWN)
            || (l[i].symbol == DOWN && l[i + 2].symbol == UP)
        {
            *err = "UP/DOWN opposition";
            return ILLEGAL_REPRESENTATION;
        }
    }
    for i in 0..4 {
        let n = l[(i + 1) % 4].symbol;
        if (l[i].symbol == UP && n == UP) || (l[i].symbol == DOWN && n == DOWN) {
            *err = "Adjacent like stereobonds";
            return ILLEGAL_REPRESENTATION;
        }
    }
    for i in 0..4 {
        let (b, c) = (&l[(i + 1) % 4], &l[(i + 2) % 4]);
        if l[i].symbol == NONE && b.symbol == NONE && c.symbol == NONE {
            let ang = angle(l[i].x - b.x, l[i].y - b.y, c.x - b.x, c.y - b.y);
            if ang < 185.0 * PI / 180.0 {
                *err = "colinearity or triangle rule violation";
                return ILLEGAL_REPRESENTATION;
            }
        }
    }
    sort_tetra(&mut t);
    parity_of_volume(&t)
}

/// `AtomParity`: parity of atom `iatom` (1 based) relative to its numbering; `err` receives
/// the text of the `stereo_error` global when `ILLEGAL_REPRESENTATION` is returned.
pub fn atom_parity(
    mp: &Molecule,
    iatom: usize,
    nbp: &Neighbourhood,
    err: &mut &'static str,
) -> i32 {
    let n = nbp.atoms.len();
    if !(3..=4).contains(&n) {
        return ILLEGAL_REPRESENTATION;
    }
    let center = &mp.atoms[iatom - 1];
    let sym = center.atom_symbol.as_str();
    let (mut multiple, mut allene, mut stereo) = (false, false, false);
    let mut lig = [StereoBond::default(); 4];
    for i in 0..n {
        let b = &mp.bonds[nbp.bonds[i]];
        if b.bond_type != SINGLE {
            multiple = true;
            let jatom = nbp.atoms[i] as i32 + 1;
            let ndb = mp
                .bonds
                .iter()
                .filter(|x| (x.atoms[0] == jatom || x.atoms[1] == jatom) && x.bond_type == DOUBLE)
                .count();
            if ndb == 2 {
                allene = true;
            }
        }
        let a = &mp.atoms[nbp.atoms[i]];
        lig[i].x = f64::from(a.x - center.x);
        lig[i].y = f64::from(a.y - center.y);
        lig[i].number = nbp.atoms[i] as i32 + 1;
        if b.atoms[0] == iatom as i32 {
            lig[i].symbol = b.stereo_symbol;
            if lig[i].symbol == UP || lig[i].symbol == DOWN {
                stereo = true;
            }
        } else {
            lig[i].symbol = NONE;
        }
    }
    if multiple && stereo && sym != "P" && sym != "S" {
        if allene {
            return ALLENE_PARITY;
        }
        *err = "AtomParity: Stereobond at unsaturated atom";
        return ILLEGAL_REPRESENTATION;
    } else if multiple && sym != "S" {
        return UNDEFINED_PARITY;
    }
    lig[0].angle = 0.0;
    for i in 1..n {
        lig[i].angle = angle(lig[0].x, lig[0].y, lig[i].x, lig[i].y);
    }
    for i in 2..n {
        for j in (2..=i).rev() {
            if lig[j].angle < lig[j - 1].angle {
                lig.swap(j, j - 1);
            } else {
                break;
            }
        }
    }
    if n == 3 {
        atom3_parity(&lig, err)
    } else {
        atom4_parity(&lig, err)
    }
}

fn candidate(sym: &str, n: usize) -> bool {
    matches!(sym, "C" | "S" | "N" | "O" | "P" | "Si") && n > 2 && n <= 4
}

fn nbs(mp: &Molecule) -> Result<Vec<Neighbourhood>, &'static str> {
    setup_neighbourhood(mp).ok_or("too many neighbours")
}

fn is_wedge(s: i32) -> bool {
    s == UP || s == DOWN
}

fn name10(mp: &Molecule) -> String {
    format!("{:>10}", mp.name)
}

/// `CheckStereo`: true if all potential stereocenters are either undefined or drawn according
/// to the MDL rules
pub fn check_stereo(mp: &Molecule, msgs: &mut Vec<String>) -> Result<bool, &'static str> {
    let nba = nbs(mp)?;
    let mut result = true;
    let mut center_defined = false;
    for (i, nbp) in nba.iter().enumerate() {
        if candidate(&mp.atoms[i].atom_symbol, nbp.atoms.len()) {
            let mut err = "";
            let parity = atom_parity(mp, i + 1, nbp, &mut err);
            if parity == ILLEGAL_REPRESENTATION {
                add_msg(
                    msgs,
                    format!("{}    atom {:3} : {}", name10(mp), i + 1, err),
                );
                result = false;
            } else if parity == EVEN_PARITY || parity == ODD_PARITY || parity == ALLENE_PARITY {
                center_defined = true;
            }
        } else {
            for &bj in &nbp.bonds {
                let b = &mp.bonds[bj];
                if b.atoms[0] == i as i32 + 1 && is_wedge(b.stereo_symbol) {
                    add_msg(
                        msgs,
                        format!(
                            "{}    atom {:3} : {}",
                            name10(mp),
                            i + 1,
                            "stereobond to non-stereogenic atom"
                        ),
                    );
                    result = false;
                }
            }
        }
    }
    if mp.chiral_flag != 0 && !center_defined {
        add_msg(msgs, "chiral flag set but no stereocenter defined".into());
        result = false;
    }
    Ok(result)
}

/// the C code compares `bond_type` with the stereo symbols UP (= SINGLE) and DOWN
fn is_flat_stereo_bond_type(t: i32) -> bool {
    t == UP || t == DOWN
}

/// `FixDubious3DMolecule`: turns tilted or 3D structures with wedges into 2D ones
pub fn fix_dubious_3d_molecule(
    mp: &mut Molecule,
    msgs: &mut Vec<String>,
) -> Result<i32, &'static str> {
    let non_zero_z = mp.atoms.iter().any(|a| a.z != 0.0);
    if !non_zero_z && mp.dimensionality == "2D" {
        return Ok(0);
    }
    let nstereo = mp
        .bonds
        .iter()
        .filter(|b| is_wedge(b.stereo_symbol))
        .count();
    let mut length = 0.0f64;
    for b in &mp.bonds {
        let (a0, a1) = (
            &mp.atoms[(b.atoms[0] - 1) as usize],
            &mp.atoms[(b.atoms[1] - 1) as usize],
        );
        let d: f32 = (a0.x - a1.x) * (a0.x - a1.x) + (a0.y - a1.y) * (a0.y - a1.y);
        length += f64::from(d);
    }
    length /= mp.bonds.len() as f64;
    length = length.sqrt();

    let nba = nbs(mp)?;
    let mut nflat_sp3 = 0;
    for (i, nbp) in nba.iter().enumerate() {
        let n = nbp.atoms.len();
        if n < 3 {
            continue;
        }
        let sym = mp.atoms[i].atom_symbol.as_str();
        if !matches!(sym, "C" | "N" | "P" | "S") {
            continue;
        }
        if mp
            .bonds
            .iter()
            .any(|b| is_flat_stereo_bond_type(b.bond_type) && i as i32 + 1 == b.atoms[0])
        {
            continue;
        }
        let pt = |k: usize| {
            let a = &mp.atoms[k];
            (f64::from(a.x), f64::from(a.y), f64::from(a.z))
        };
        let mut tetra = [NPoint::default(); 4];
        (tetra[0].x, tetra[0].y, tetra[0].z) = pt(i);
        let mut i1 = 0;
        while i1 < n {
            if mp.bonds[nbp.bonds[i1]].bond_type != SINGLE && sym != "P" && sym != "S" {
                break;
            }
            i1 += 1;
        }
        if i1 >= n {
            continue;
        }
        let mut stereo_triple = 0;
        i1 = 0;
        while i1 < n {
            (tetra[1].x, tetra[1].y, tetra[1].z) = pt(nbp.atoms[i1]);
            if is_flat_stereo_bond_type(mp.bonds[nbp.bonds[i1]].bond_type) {
                stereo_triple |= 1;
            }
            let mut i2 = i1 + 1;
            while i2 < n {
                (tetra[2].x, tetra[2].y, tetra[2].z) = pt(nbp.atoms[i2]);
                if is_flat_stereo_bond_type(mp.bonds[nbp.bonds[i2]].bond_type) {
                    stereo_triple |= 2;
                }
                let mut i3 = i2 + 1;
                while i3 < n {
                    (tetra[3].x, tetra[3].y, tetra[3].z) = pt(nbp.atoms[i3]);
                    if is_flat_stereo_bond_type(mp.bonds[nbp.bonds[i3]].bond_type) {
                        stereo_triple |= 4;
                    }
                    let vol = volume(&tetra).abs();
                    if stereo_triple != 0 {
                        if vol < 0.01 * length * length * length {
                            nflat_sp3 += 1;
                            break;
                        }
                        stereo_triple &= !4;
                    }
                    i3 += 1;
                }
                stereo_triple &= !2;
                if i3 >= n {
                    break;
                }
                i2 += 1;
            }
            stereo_triple &= !1;
            if i2 >= n {
                break;
            }
            i1 += 1;
        }
    }

    let mut result = 0;
    if non_zero_z && mp.dimensionality == "2D" {
        for a in &mut mp.atoms {
            a.z = 0.0;
        }
        result |= ZEROED_Z_COORDINATES;
        add_msg(msgs, "Cleared z-coordinates in 2D MOL file".into());
    } else if non_zero_z && nstereo > 0 && nflat_sp3 > 0 {
        mp.dimensionality = "2D".into();
        result |= CONVERTED_TO_2D;
        for a in &mut mp.atoms {
            a.z = 0.0;
        }
        result |= ZEROED_Z_COORDINATES;
        add_msg(msgs, "Cleared z-coordinates of tilted 2D MOL file".into());
    }
    Ok(result)
}

/// `DubiousStereochemistry`: bit set describing ill-defined stereochemistry in `mp`
pub fn dubious_stereochemistry(mp: &Molecule, msgs: &mut Vec<String>) -> Result<i32, &'static str> {
    let nba = nbs(mp)?;
    let mut result = 0;
    for b in &mp.bonds {
        if b.stereo_symbol == EITHER {
            add_msg(
                msgs,
                format!(
                    "{}    EITHER bond {}-{}",
                    name10(mp),
                    b.atoms[0],
                    b.atoms[1]
                ),
            );
            result |= EITHER_BOND_FOUND;
        }
    }
    // `ndb` is declared once for the whole function in C and only set for double bonds
    let mut ndb = 0;
    for (i, nbp) in nba.iter().enumerate() {
        let n = nbp.atoms.len();
        let mut nmulti = 0;
        let mut is_allene = false;
        for j in 0..n {
            let bt = mp.bonds[nbp.bonds[j]].bond_type;
            if bt != SINGLE {
                nmulti += 1;
                if bt == DOUBLE {
                    let jatom = nbp.atoms[j];
                    ndb = nba[jatom]
                        .bonds
                        .iter()
                        .filter(|&&k| mp.bonds[k].bond_type == DOUBLE)
                        .count();
                }
                if ndb == 2 {
                    is_allene = true;
                }
            }
        }
        let sym = mp.atoms[i].atom_symbol.as_str();
        let stereogenic = (sym == "C" && n > 2 && n <= 4 && nmulti == 0)
            || (sym == "C" && (2..=3).contains(&n) && nmulti == 1 && is_allene)
            || (sym == "S" && n > 2 && n <= 4)
            || (sym == "N" && n > 3 && n <= 4 && nmulti == 0)
            || (sym == "P" && n > 2 && n <= 4)
            || (sym == "Si" && n > 2 && n <= 4 && nmulti == 0);
        if !stereogenic {
            for &bj in &nbp.bonds {
                let b = &mp.bonds[bj];
                if b.atoms[0] == i as i32 + 1 && is_wedge(b.stereo_symbol) {
                    add_msg(
                        msgs,
                        format!(
                            "{}    atom {:3} : {}",
                            name10(mp),
                            i + 1,
                            "stereobond to non-stereogenic atom"
                        ),
                    );
                    result |= STEREO_BOND_AT_NON_STEREO_ATOM;
                }
            }
        }
    }
    Ok(result)
}

/// `RemoveDubiousStereochemistry`: removes EITHER marks and wedges at non-stereogenic atoms
pub fn remove_dubious_stereochemistry(mp: &mut Molecule) -> Result<(), &'static str> {
    let nba = nbs(mp)?;
    for b in &mut mp.bonds {
        if b.stereo_symbol == EITHER {
            b.stereo_symbol = NONE;
        }
    }
    for (i, nbp) in nba.iter().enumerate() {
        let n = nbp.atoms.len();
        let nmulti = nbp
            .bonds
            .iter()
            .filter(|&&k| mp.bonds[k].bond_type != SINGLE)
            .count();
        let sym = mp.atoms[i].atom_symbol.as_str();
        let stereogenic = (sym == "C" && n > 2 && n <= 4 && nmulti == 0)
            || (sym == "S" && n > 2 && n <= 4)
            || (sym == "N" && n > 3 && n <= 4 && nmulti == 0)
            || (sym == "P" && n > 2 && n <= 4)
            || (sym == "Si" && n > 2 && n <= 4 && nmulti == 0);
        if !stereogenic {
            for &bj in &nbp.bonds {
                let b = &mut mp.bonds[bj];
                if b.atoms[0] == i as i32 + 1 && is_wedge(b.stereo_symbol) {
                    b.stereo_symbol = NONE;
                }
            }
        }
    }
    Ok(())
}

/// Atom/atom and atom/bond collision limit as a fraction of the median bond length
#[derive(Debug, Clone, Copy)]
pub struct ClashLimit(pub f64);

impl Default for ClashLimit {
    fn default() -> Self {
        ClashLimit(0.15)
    }
}

impl ClashLimit {
    /// `SetCollisionLimit`
    pub fn from_percent(percent: i32) -> Self {
        ClashLimit(f64::from(percent) / 100.0)
    }
}

/// `AtomClash`: true if two atoms, or an atom and a bond, come too close
pub fn atom_clash(mp: &Molecule, limit: ClashLimit, msgs: &mut Vec<String>) -> bool {
    if mp.bonds.is_empty() {
        return false;
    }
    let clash_limit = limit.0;
    let mut bl: Vec<f64> = mp
        .bonds
        .iter()
        .map(|b| {
            let (a1, a2) = (
                &mp.atoms[(b.atoms[0] - 1) as usize],
                &mp.atoms[(b.atoms[1] - 1) as usize],
            );
            let d: f32 = (a1.x - a2.x) * (a1.x - a2.x) + (a1.y - a2.y) * (a1.y - a2.y);
            f64::from(d)
        })
        .collect();
    for i in 1..bl.len() {
        for j in (0..i).rev() {
            if bl[j] > bl[j + 1] {
                bl.swap(j, j + 1);
            } else {
                break;
            }
        }
    }
    let bond_square_median = bl[mp.bonds.len() / 2];

    let n = mp.atoms.len();
    for i in 0..n {
        for j in i + 1..n {
            let (a1, a2) = (&mp.atoms[i], &mp.atoms[j]);
            let d: f32 = (a1.x - a2.x) * (a1.x - a2.x) + (a1.y - a2.y) * (a1.y - a2.y);
            let dist = f64::from(d);
            if dist < clash_limit * clash_limit * bond_square_median {
                add_msg(
                    msgs,
                    format!(
                        "{}    atom {:3} : only {}% of avg. bond away from atom {}",
                        name10(mp),
                        j + 1,
                        fmt_gp(100.0 * (dist / (bond_square_median + EPS)).sqrt(), 2),
                        i + 1
                    ),
                );
                return true;
            }
        }
    }

    for b in &mp.bonds {
        for (j, ap) in mp.atoms.iter().enumerate() {
            if (b.atoms[0] - 1) as usize == j || (b.atoms[1] - 1) as usize == j {
                continue;
            }
            let ap1 = &mp.atoms[(b.atoms[0] - 1) as usize];
            let ap2 = &mp.atoms[(b.atoms[1] - 1) as usize];
            let rr = f64::from((ap.x - ap1.x) * (ap.x - ap1.x) + (ap.y - ap1.y) * (ap.y - ap1.y));
            let bb =
                f64::from((ap2.x - ap1.x) * (ap2.x - ap1.x) + (ap2.y - ap1.y) * (ap2.y - ap1.y));
            let rb = f64::from((ap.x - ap1.x) * (ap2.x - ap1.x) + (ap.y - ap1.y) * (ap2.y - ap1.y));
            if 0.0 <= rb
                && rb <= bb
                && (rr * bb - rb * rb) / (bb + EPS) < clash_limit * clash_limit * bond_square_median
            {
                add_msg(
                    msgs,
                    format!(
                        "{}    atom {:3} : only {}% {} from bond {}-{}",
                        name10(mp),
                        j + 1,
                        fmt_gp(
                            100.0 * ((rr * bb - rb * rb) / (bb * bond_square_median + EPS)).sqrt(),
                            2
                        ),
                        "of average bond length",
                        b.atoms[0],
                        b.atoms[1]
                    ),
                );
                return true;
            }
        }
    }
    false
}

/// `CisTransPerception`: sets the bond colour of defined double bonds to CIS, TRANS or NONE
/// depending on the ligands with the lowest `numbering`; returns the number of defined bonds
pub fn cis_trans_perception(mp: &mut Molecule, numbering: &[i32]) -> Result<i32, &'static str> {
    let nba = nbs(mp)?;
    let mut result = 0;
    for b in &mut mp.bonds {
        b.color = NONE;
    }
    let maxnum = numbering
        .iter()
        .take(mp.atoms.len())
        .copied()
        .fold(-1, i32::max)
        + 1;
    for i in 0..mp.bonds.len() {
        let bp = &mp.bonds[i];
        if bp.bond_type != DOUBLE || bp.stereo_symbol == CIS_TRANS_EITHER {
            continue;
        }
        let (j1, j2) = ((bp.atoms[0] - 1) as usize, (bp.atoms[1] - 1) as usize);
        if !atom_symbol_match(&mp.atoms[j1].atom_symbol, "C,N")
            || !atom_symbol_match(&mp.atoms[j2].atom_symbol, "C,N")
        {
            continue;
        }
        let (n1, n2) = (nba[j1].atoms.len(), nba[j2].atoms.len());
        if n1 <= 1 || n2 <= 1 || n1 > 3 || n2 > 3 {
            continue;
        }
        let lowest = |j: usize, other: usize| -> Option<usize> {
            let mut equal = false;
            let mut nmin = maxnum;
            let mut at = 0;
            for &k in &nba[j].atoms {
                if k != other {
                    if numbering[k] < nmin {
                        at = k;
                        nmin = numbering[k];
                    } else if numbering[k] == nmin {
                        equal = true;
                    }
                }
            }
            if equal { None } else { Some(at) }
        };
        let Some(at1) = lowest(j1, j2) else { continue };
        let Some(at2) = lowest(j2, j1) else { continue };
        let (a1, a2) = (&mp.atoms[j1], &mp.atoms[j2]);
        let x21 = f64::from(mp.atoms[at1].x - a1.x);
        let y21 = f64::from(mp.atoms[at1].y - a1.y);
        let x23 = f64::from(a2.x - a1.x);
        let y23 = f64::from(a2.y - a1.y);
        let (x32, y32) = (-x23, -y23);
        let x34 = f64::from(mp.atoms[at2].x - a2.x);
        let y34 = f64::from(mp.atoms[at2].y - a2.y);
        let sign = (x21 * y23 - x23 * y21) * (x32 * y34 - x34 * y32);
        if -0.001 < sign && sign < 0.001 {
            continue;
        }
        result += 1;
        mp.bonds[i].color = if sign > 0.0 { CIS } else { TRANS };
    }
    Ok(result)
}

/// `NoParityDefined`: true (despite the name) if `mp` has a defined even or odd stereocenter
pub fn no_parity_defined(mp: &Molecule) -> Result<bool, &'static str> {
    let nba = nbs(mp)?;
    for (i, nbp) in nba.iter().enumerate() {
        if candidate(&mp.atoms[i].atom_symbol, nbp.atoms.len()) {
            let mut err = "";
            let parity = atom_parity(mp, i + 1, nbp, &mut err);
            if parity == EVEN_PARITY || parity == ODD_PARITY {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// `PermutationParity`: parity of the number of exchanges that sort `indices`
pub fn permutation_parity(indices: &mut [i32]) -> i32 {
    let mut nswap = 0;
    for i in 1..indices.len() {
        for j in (0..i).rev() {
            if indices[j] > indices[j + 1] {
                nswap += 1;
                indices.swap(j, j + 1);
            } else {
                break;
            }
        }
    }
    if nswap % 2 == 0 {
        EVEN_PARITY
    } else {
        ODD_PARITY
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::reaccsio::{mol_str_to_mol, mol_to_mol_str};

    /// every check must report what the Avalon C code reports
    #[test]
    fn stereo_checks_match_c() {
        let text = include_str!("../../../../../tests/data/avalon_stereo.jsonl");
        let mut bad = Vec::new();
        let mut n = 0;
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().expect("input");
            let want = v["output"].as_str().expect("output");
            let mut m = mol_str_to_mol(input).expect("molecule");
            let mut msgs = Vec::new();
            let mut out = String::new();
            out += &format!(
                "check={}\n",
                i32::from(check_stereo(&m, &mut msgs).unwrap())
            );
            out += &format!("noparity={}\n", i32::from(no_parity_defined(&m).unwrap()));
            out += &format!(
                "dubious={}\n",
                dubious_stereochemistry(&m, &mut msgs).unwrap()
            );
            out += &format!(
                "clash={}\n",
                i32::from(atom_clash(&m, ClashLimit::default(), &mut msgs))
            );
            out += &format!(
                "clash40={}\n",
                i32::from(atom_clash(&m, ClashLimit::from_percent(40), &mut msgs))
            );
            let numbering: Vec<i32> = (0..m.atoms.len() as i32)
                .map(|i| if i % 3 == 0 { i / 2 } else { i })
                .collect();
            out += &format!(
                "cistrans={}\n",
                cis_trans_perception(&mut m, &numbering).unwrap()
            );
            for b in &m.bonds {
                out += &format!("c{} ", b.color);
            }
            out += "\n";
            out += &format!(
                "fix={}\n",
                fix_dubious_3d_molecule(&mut m, &mut msgs).unwrap()
            );
            remove_dubious_stereochemistry(&mut m).unwrap();
            out += "msgs:\n";
            for msg in &msgs {
                out += msg;
                out += "\n";
            }
            out += "mol:\n";
            out += &mol_to_mol_str(&m);
            n += 1;
            if out != want {
                bad.push(format!("{input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(
            bad.is_empty(),
            "{} of {n} differ; first:\n{}",
            bad.len(),
            bad[0]
        );
    }

    #[test]
    fn permutation_parities() {
        assert_eq!(permutation_parity(&mut [1, 2, 3]), EVEN_PARITY);
        assert_eq!(permutation_parity(&mut [2, 1, 3]), ODD_PARITY);
        assert_eq!(permutation_parity(&mut [3, 2, 1]), ODD_PARITY);
        assert_eq!(invert_parity(EVEN_PARITY), ODD_PARITY);
    }
}
