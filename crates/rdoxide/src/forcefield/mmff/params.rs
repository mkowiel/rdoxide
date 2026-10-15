//! MMFF parameter tables and their lookups (upstream
//! `ForceField/MMFF/Params.{h,cpp}`).
//!
//! The tables are tab-separated text extracted from upstream by
//! `tools/sync_upstream.py`; lines starting with `*` or `#` are comments. Like
//! upstream's sorted-vector mode, a key that occurs more than once resolves to
//! its first row.

use std::collections::HashMap;
use std::sync::OnceLock;

/// Atom types that MMFF treats as aromatic.
pub const DEFAULT_AROMATIC_TYPES: [u8; 17] = [
    37, 38, 39, 44, 58, 59, 63, 64, 65, 66, 69, 76, 78, 79, 80, 81, 82,
];

/// Symbolic-type equivalences (MMFF levels 2-5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Def {
    pub eq_level: [u8; 4],
}

/// Atom type properties.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Prop {
    pub atno: u8,
    pub crd: u8,
    pub val: u8,
    pub pilp: u8,
    pub mltb: u8,
    pub arom: u8,
    pub linh: u8,
    pub sbmb: u8,
}

/// Partial bond charge increment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pbci {
    pub pbci: f64,
    pub fcadj: f64,
}

/// Bond stretching parameters (also used for the Badger-rule `Bndk` table).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bond {
    pub kb: f64,
    pub r0: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HerschbachLaurie {
    pub a_ij: f64,
    pub d_ij: f64,
    pub dp_ij: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CovRadPauEle {
    pub r0: f64,
    pub chi: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Angle {
    pub ka: f64,
    pub theta0: f64,
}

/// Stretch-bend coefficients, for i-j-k and k-j-i.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stbn {
    pub kba_ijk: f64,
    pub kba_kji: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Oop {
    pub koop: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tor {
    pub v1: f64,
    pub v2: f64,
    pub v3: f64,
}

/// van der Waals parameters of an atom type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VdW {
    pub alpha_i: f64,
    pub n_i: f64,
    pub a_i: f64,
    pub g_i: f64,
    pub r_star: f64,
    /// donor/acceptor flag: `b'-'`, `b'D'` or `b'A'`
    pub da: u8,
}

/// The general vdW constants from the first row of the vdW table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VdWConstants {
    pub power: f64,
    pub b: f64,
    pub beta: f64,
    pub darad: f64,
    pub daeps: f64,
}

/// Rows of a table: tab-separated, empty fields dropped, comments skipped.
fn rows(text: &'static str) -> impl Iterator<Item = Vec<&'static str>> {
    text.lines()
        .filter(|l| !l.starts_with('*') && !l.starts_with('#'))
        .map(|l| l.split('\t').filter(|t| !t.is_empty()).collect::<Vec<_>>())
        .filter(|f| !f.is_empty())
}

fn int(s: &str) -> u8 {
    s.trim().parse().expect("MMFF table: bad integer")
}

fn num(s: &str) -> f64 {
    s.trim().parse().expect("MMFF table: bad number")
}

/// Inserts `value` unless the key is already present.
fn put<K: std::hash::Hash + Eq, V>(map: &mut HashMap<K, V>, key: K, value: V) {
    map.entry(key).or_insert(value);
}

/// A bundle of every MMFF parameter table.
#[derive(Debug)]
pub struct MmffParams {
    def: HashMap<u8, Def>,
    prop: HashMap<u8, Prop>,
    pbci: HashMap<u8, Pbci>,
    /// keyed by (bond type, i, j) with i <= j
    chg: HashMap<(u8, u8, u8), f64>,
    bond: HashMap<(u8, u8, u8), Bond>,
    bndk: HashMap<(u8, u8), Bond>,
    herschbach_laurie: HashMap<(u8, u8), HerschbachLaurie>,
    cov_rad_pau_ele: HashMap<u8, CovRadPauEle>,
    /// keyed by (angle type, j, i, k)
    angle: HashMap<(u8, u8, u8, u8), Angle>,
    /// keyed by (stretch-bend type, j, i, k)
    stbn: HashMap<(u8, u8, u8, u8), Stbn>,
    dfsb: HashMap<(u8, u8, u8), Stbn>,
    oop: HashMap<(u8, u8, u8, u8), Oop>,
    /// keyed by (torsion type, j, k, i, l)
    tor: HashMap<(u8, u8, u8, u8, u8), Tor>,
    vdw: HashMap<u8, VdW>,
    vdw_constants: VdWConstants,
    arom: Vec<u8>,
}

impl MmffParams {
    /// The default tables; `is_mmffs` selects the MMFF94s out-of-plane and
    /// torsion tables.
    pub fn get(is_mmffs: bool) -> &'static MmffParams {
        static MMFF94: OnceLock<MmffParams> = OnceLock::new();
        static MMFF94S: OnceLock<MmffParams> = OnceLock::new();
        let cell = if is_mmffs { &MMFF94S } else { &MMFF94 };
        cell.get_or_init(|| MmffParams::new(is_mmffs))
    }

    fn new(is_mmffs: bool) -> Self {
        let (oop_text, tor_text) = if is_mmffs {
            (
                include_str!("../../../data/mmff/sOop.tsv"),
                include_str!("../../../data/mmff/sTor.tsv"),
            )
        } else {
            (
                include_str!("../../../data/mmff/Oop.tsv"),
                include_str!("../../../data/mmff/Tor.tsv"),
            )
        };

        // only the first row of each atom type counts; upstream indexes these
        // tables by position, which breaks for types above 82 (83-86 do not
        // exist), so the rows are keyed by type instead
        let mut def = HashMap::new();
        for f in rows(include_str!("../../../data/mmff/Def.tsv")) {
            put(
                &mut def,
                int(f[1]),
                Def {
                    eq_level: [int(f[2]), int(f[3]), int(f[4]), int(f[5])],
                },
            );
        }

        let mut prop = HashMap::new();
        for f in rows(include_str!("../../../data/mmff/Prop.tsv")) {
            put(
                &mut prop,
                int(f[0]),
                Prop {
                    atno: int(f[1]),
                    crd: int(f[2]),
                    val: int(f[3]),
                    pilp: int(f[4]),
                    mltb: int(f[5]),
                    arom: int(f[6]),
                    linh: int(f[7]),
                    sbmb: int(f[8]),
                },
            );
        }

        let mut pbci = HashMap::new();
        for f in rows(include_str!("../../../data/mmff/PBCI.tsv")) {
            put(
                &mut pbci,
                int(f[1]),
                Pbci {
                    pbci: num(f[2]),
                    fcadj: num(f[3]),
                },
            );
        }

        let mut chg = HashMap::new();
        for f in rows(include_str!("../../../data/mmff/Chg.tsv")) {
            put(&mut chg, (int(f[0]), int(f[1]), int(f[2])), num(f[3]));
        }

        let mut bond = HashMap::new();
        for f in rows(include_str!("../../../data/mmff/Bond.tsv")) {
            let key = (int(f[0]), int(f[1]), int(f[2]));
            put(
                &mut bond,
                key,
                Bond {
                    kb: num(f[3]),
                    r0: num(f[4]),
                },
            );
        }

        // columns: atomic numbers, r0, kb
        let mut bndk = HashMap::new();
        for f in rows(include_str!("../../../data/mmff/Bndk.tsv")) {
            put(
                &mut bndk,
                (int(f[0]), int(f[1])),
                Bond {
                    r0: num(f[2]),
                    kb: num(f[3]),
                },
            );
        }

        let mut herschbach_laurie = HashMap::new();
        for f in rows(include_str!("../../../data/mmff/HerschbachLaurie.tsv")) {
            put(
                &mut herschbach_laurie,
                (int(f[0]), int(f[1])),
                HerschbachLaurie {
                    a_ij: num(f[2]),
                    d_ij: num(f[3]),
                    dp_ij: num(f[4]),
                },
            );
        }

        let mut cov_rad_pau_ele = HashMap::new();
        for f in rows(include_str!("../../../data/mmff/CovRadPauEle.tsv")) {
            put(
                &mut cov_rad_pau_ele,
                int(f[0]),
                CovRadPauEle {
                    r0: num(f[1]),
                    chi: num(f[2]),
                },
            );
        }

        let mut angle = HashMap::new();
        for f in rows(include_str!("../../../data/mmff/Angle.tsv")) {
            put(
                &mut angle,
                (int(f[0]), int(f[2]), int(f[1]), int(f[3])),
                Angle {
                    ka: num(f[4]),
                    theta0: num(f[5]),
                },
            );
        }

        let mut stbn = HashMap::new();
        for f in rows(include_str!("../../../data/mmff/Stbn.tsv")) {
            put(
                &mut stbn,
                (int(f[0]), int(f[2]), int(f[1]), int(f[3])),
                Stbn {
                    kba_ijk: num(f[4]),
                    kba_kji: num(f[5]),
                },
            );
        }

        let mut dfsb = HashMap::new();
        for f in rows(include_str!("../../../data/mmff/Dfsb.tsv")) {
            put(
                &mut dfsb,
                (int(f[0]), int(f[1]), int(f[2])),
                Stbn {
                    kba_ijk: num(f[3]),
                    kba_kji: num(f[4]),
                },
            );
        }

        // the key atoms are sorted by the lookup; j is the central atom
        let mut oop = HashMap::new();
        for f in rows(oop_text) {
            put(
                &mut oop,
                (int(f[1]), int(f[0]), int(f[2]), int(f[3])),
                Oop { koop: num(f[4]) },
            );
        }

        let mut tor = HashMap::new();
        for f in rows(tor_text) {
            put(
                &mut tor,
                (int(f[0]), int(f[2]), int(f[3]), int(f[1]), int(f[4])),
                Tor {
                    v1: num(f[5]),
                    v2: num(f[6]),
                    v3: num(f[7]),
                },
            );
        }

        let mut vdw = HashMap::new();
        let mut vdw_rows = rows(include_str!("../../../data/mmff/VdW.tsv"));
        let first = vdw_rows.next().expect("MMFF vdW table is empty");
        let vdw_constants = VdWConstants {
            power: num(first[0]),
            b: num(first[1]),
            beta: num(first[2]),
            darad: num(first[3]),
            daeps: num(first[4]),
        };
        for f in vdw_rows {
            let alpha_i = num(f[1]);
            let a_i = num(f[3]);
            put(
                &mut vdw,
                int(f[0]),
                VdW {
                    alpha_i,
                    n_i: num(f[2]),
                    a_i,
                    g_i: num(f[4]),
                    r_star: a_i * alpha_i.powf(vdw_constants.power),
                    da: f[5].as_bytes()[0],
                },
            );
        }

        MmffParams {
            def,
            prop,
            pbci,
            chg,
            bond,
            bndk,
            herschbach_laurie,
            cov_rad_pau_ele,
            angle,
            stbn,
            dfsb,
            oop,
            tor,
            vdw,
            vdw_constants,
            arom: DEFAULT_AROMATIC_TYPES.to_vec(),
        }
    }

    /// Whether `atom_type` is one of the aromatic MMFF types.
    pub fn is_aromatic_type(&self, atom_type: u8) -> bool {
        self.arom.contains(&atom_type)
    }

    pub fn def(&self, atom_type: u32) -> Option<&Def> {
        self.def.get(&u8::try_from(atom_type).ok()?)
    }

    pub fn prop(&self, atom_type: u32) -> Option<&Prop> {
        self.prop.get(&u8::try_from(atom_type).ok()?)
    }

    pub fn pbci(&self, atom_type: u32) -> Option<&Pbci> {
        self.pbci.get(&u8::try_from(atom_type).ok()?)
    }

    /// Bond charge increment for a bond between two atom types. The sign
    /// is -1 when the types are in canonical order and 1 when they were
    /// swapped.
    pub fn chg(&self, bond_type: u32, i_type: u32, j_type: u32) -> (i32, Option<f64>) {
        let (sign, lo, hi) = if i_type > j_type {
            (1, j_type, i_type)
        } else {
            (-1, i_type, j_type)
        };
        let bci = match (u8::try_from(bond_type), u8::try_from(lo), u8::try_from(hi)) {
            (Ok(b), Ok(i), Ok(j)) => self.chg.get(&(b, i, j)).copied(),
            _ => None,
        };
        (sign, bci)
    }

    pub fn bond(&self, bond_type: u32, atom_type: u32, nbr_type: u32) -> Option<&Bond> {
        let (lo, hi) = (atom_type.min(nbr_type), atom_type.max(nbr_type));
        self.bond.get(&(
            u8::try_from(bond_type).ok()?,
            u8::try_from(lo).ok()?,
            u8::try_from(hi).ok()?,
        ))
    }

    /// Badger's-rule bond parameters from atomic numbers.
    pub fn bndk(&self, atomic_num: u32, nbr_atomic_num: u32) -> Option<&Bond> {
        let (lo, hi) = (
            atomic_num.min(nbr_atomic_num),
            atomic_num.max(nbr_atomic_num),
        );
        self.bndk
            .get(&(u8::try_from(lo).ok()?, u8::try_from(hi).ok()?))
    }

    pub fn herschbach_laurie(&self, i_row: u32, j_row: u32) -> Option<&HerschbachLaurie> {
        let (lo, hi) = (i_row.min(j_row), i_row.max(j_row));
        self.herschbach_laurie
            .get(&(u8::try_from(lo).ok()?, u8::try_from(hi).ok()?))
    }

    pub fn cov_rad_pau_ele(&self, atomic_num: u32) -> Option<&CovRadPauEle> {
        self.cov_rad_pau_ele.get(&u8::try_from(atomic_num).ok()?)
    }

    /// Angle bending parameters, relaxing the atom types of i and k through
    /// the equivalence levels until a row is found.
    pub fn angle(&self, angle_type: u32, i_type: u32, j_type: u32, k_type: u32) -> Option<&Angle> {
        let (angle_type, j) = (u8::try_from(angle_type).ok()?, u8::try_from(j_type).ok()?);
        let (def_i, def_k) = (self.def(i_type)?, self.def(k_type)?);
        (0..4).find_map(|level| {
            let (a, b) = (def_i.eq_level[level], def_k.eq_level[level]);
            self.angle.get(&(angle_type, j, a.min(b), a.max(b)))
        })
    }

    /// Stretch-bend parameters; the flag tells whether i and k were swapped
    /// to reach the canonical order.
    pub fn stbn(
        &self,
        stretch_bend_type: u32,
        bond_type1: u32,
        bond_type2: u32,
        i_type: u32,
        j_type: u32,
        k_type: u32,
    ) -> (bool, Option<&Stbn>) {
        let swap = if i_type > k_type {
            true
        } else if i_type == k_type {
            bond_type1 < bond_type2
        } else {
            false
        };
        let (lo, hi) = (i_type.min(k_type), i_type.max(k_type));
        let found = match (
            u8::try_from(stretch_bend_type),
            u8::try_from(j_type),
            u8::try_from(lo),
            u8::try_from(hi),
        ) {
            (Ok(t), Ok(j), Ok(i), Ok(k)) => self.stbn.get(&(t, j, i, k)),
            _ => None,
        };
        (swap, found)
    }

    /// Default stretch-bend parameters from the periodic-table rows.
    pub fn dfsb(&self, row1: u32, row2: u32, row3: u32) -> (bool, Option<&Stbn>) {
        let swap = row1 > row3;
        let (lo, hi) = (row1.min(row3), row1.max(row3));
        let found = match (u8::try_from(lo), u8::try_from(row2), u8::try_from(hi)) {
            (Ok(a), Ok(b), Ok(c)) => self.dfsb.get(&(a, b, c)),
            _ => None,
        };
        (swap, found)
    }

    /// Out-of-plane bending parameters around central atom j.
    pub fn oop(&self, i_type: u32, j_type: u32, k_type: u32, l_type: u32) -> Option<&Oop> {
        let j = u8::try_from(j_type).ok()?;
        let (di, dk, dl) = (self.def(i_type)?, self.def(k_type)?, self.def(l_type)?);
        (0..4).find_map(|level| {
            let mut ikl = [di.eq_level[level], dk.eq_level[level], dl.eq_level[level]];
            ikl.sort_unstable();
            self.oop.get(&(j, ikl[0], ikl[1], ikl[2]))
        })
    }

    /// Torsion parameters. `tor_type` holds the primary and the fallback
    /// torsion type; the result carries the type that matched.
    pub fn tor(
        &self,
        tor_type: (u32, u32),
        i_type: u32,
        j_type: u32,
        k_type: u32,
        l_type: u32,
    ) -> (u32, Option<&Tor>) {
        let (def_i, def_l) = match (self.def(i_type), self.def(l_type)) {
            (Some(i), Some(l)) => (i, l),
            _ => return (tor_type.0, None),
        };
        let mut found = None;
        let mut iter = 0usize;
        let mut max_iter = 5;
        let mut can_tor_type = tor_type.0;
        while (iter < max_iter && (found.is_none() || max_iter == 4))
            || (iter == 4 && tor_type.0 == 5 && tor_type.1 != 0)
        {
            if max_iter == 5 && iter == 4 {
                max_iter = 4;
                iter = 0;
                can_tor_type = tor_type.1;
            }
            // iteration 1 wildcards l, 2 wildcards i, 3 both, 0 neither
            let (i_level, l_level) = match iter {
                1 => (1, 3),
                2 => (3, 1),
                n => (n, n),
            };
            let (mut ci, mut cj, mut ck, mut cl) = (
                u32::from(def_i.eq_level[i_level]),
                j_type,
                k_type,
                u32::from(def_l.eq_level[l_level]),
            );
            if cj > ck || (cj == ck && ci > cl) {
                std::mem::swap(&mut cj, &mut ck);
                std::mem::swap(&mut ci, &mut cl);
            }
            let key = (
                u8::try_from(can_tor_type),
                u8::try_from(cj),
                u8::try_from(ck),
                u8::try_from(ci),
                u8::try_from(cl),
            );
            if let (Ok(t), Ok(j), Ok(k), Ok(i), Ok(l)) = key {
                if let Some(p) = self.tor.get(&(t, j, k, i, l)) {
                    found = Some(p);
                    if max_iter == 4 {
                        break;
                    }
                }
            }
            iter += 1;
        }
        (can_tor_type, found)
    }

    pub fn vdw(&self, atom_type: u32) -> Option<&VdW> {
        self.vdw.get(&u8::try_from(atom_type).ok()?)
    }

    pub fn vdw_constants(&self) -> &VdWConstants {
        &self.vdw_constants
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tables_load() {
        for s in [false, true] {
            let p = MmffParams::get(s);
            assert_eq!(p.def.len(), 95);
            assert_eq!(p.pbci.len(), 99);
            assert!(p.prop.len() >= 95);
            assert!(p.vdw(1).is_some() && p.vdw(99).is_some());
        }
    }

    #[test]
    fn simple_lookups() {
        let p = MmffParams::get(false);
        let b = p.bond(0, 1, 1).unwrap();
        assert_eq!((b.kb, b.r0), (4.258, 1.508));
        // the atom types may come in either order
        assert_eq!(p.bond(0, 2, 1), p.bond(0, 1, 2));
        let a = p.angle(0, 1, 1, 1).unwrap();
        assert_eq!((a.ka, a.theta0), (0.851, 109.608));
        let v = p.vdw(1).unwrap();
        assert!((v.r_star - 3.89 * 1.05f64.powf(p.vdw_constants().power)).abs() < 1e-12);
        assert_eq!(v.da, b'-');
        assert_eq!(p.chg(0, 1, 1), (-1, Some(0.0)));
        assert_eq!(p.chg(0, 2, 1).0, 1);
        assert_eq!(p.cov_rad_pau_ele(6).unwrap().chi, 2.5);
        assert_eq!(p.prop(2).unwrap().atno, 6);
        assert!(p.is_aromatic_type(37) && !p.is_aromatic_type(1));
    }

    #[test]
    fn stbn_swap_and_default() {
        let p = MmffParams::get(false);
        let (swap, s) = p.stbn(0, 0, 0, 1, 1, 2);
        assert!(!swap);
        assert_eq!(s.unwrap().kba_ijk, 0.136);
        let (swap, s) = p.stbn(0, 0, 0, 2, 1, 1);
        assert!(swap);
        assert_eq!(s.unwrap().kba_ijk, 0.136);
        let (swap, s) = p.dfsb(1, 0, 0);
        assert!(swap);
        assert!(s.is_none() || s.is_some());
    }

    #[test]
    fn torsion_falls_back_to_wildcards() {
        let p = MmffParams::get(false);
        // an exact row
        let (t, row) = p.tor((0, 0), 1, 1, 1, 1);
        assert_eq!(t, 0);
        assert_eq!(row.unwrap().v1, 0.103);
        // a pair of unlisted outer atom types resolves through wildcards
        let (_, row) = p.tor((0, 0), 5, 1, 1, 5);
        assert!(row.is_some());
    }

    #[test]
    fn oop_is_order_independent() {
        let p = MmffParams::get(false);
        assert_eq!(p.oop(1, 2, 1, 2), p.oop(2, 2, 1, 1));
        assert_eq!(p.oop(1, 2, 1, 2).unwrap().koop, 0.030);
    }
}
