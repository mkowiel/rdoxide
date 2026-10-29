//! Overlap integrals between Slater-type orbitals (YAeHMOP `abfns.c`, `lovlap.c`, `mov.c` and
//! `R_overlap_mat.c`, molecular case).

use super::params::Atom;

const SQRT3: f64 = 1.732050807568877;
const SQRT6: f64 = 2.449489742783178;
const SQRT10: f64 = 3.16227766016838;
const SQRT15: f64 = 3.872983346207417;
/// 1/bohr radius
const AUI: f64 = 1.889644746;

const BP: usize = 1;
const EP: usize = 3;
const BD: usize = 4;
const ED: usize = 8;
const BF: usize = 9;
const EF: usize = 15;

fn fifmin0(a: i32, b: i32) -> i32 {
    a.min(b)
}

/// The A and B auxiliary functions of the overlap evaluation (`abfns`)
fn abfns(a: &mut [f64; 30], b: &mut [f64; 30], sk1: f64, sk2: f64, rr: f64, maxcal: i32) {
    let j = maxcal + 1;
    let rho1 = 0.5 * (sk1 + sk2) * rr;
    let rho2 = 0.5 * (sk1 - sk2) * rr;
    if rho1.abs() > 165.0 || rho2.abs() > 165.0 {
        a[..20].fill(0.0);
        b[..20].fill(0.0);
        return;
    }
    let c = (-rho1).exp();
    a[0] = c / rho1;
    for i in 2..=j as usize {
        a[i - 1] = ((i - 1) as f64 * a[i - 2] + c) / rho1;
    }
    let ix = j;
    let ir = (2.0 * rho2).abs() as i32;
    let is = fifmin0(ir + 1, 19);
    if rho2 == 0.0 {
        let mut i = 1;
        while i <= ix {
            b[i as usize - 1] = 2.0 / i as f64;
            b[i as usize] = 0.0;
            i += 2;
        }
        return;
    }
    let d = rho2.exp();
    let h = 1.0 / d;
    let mut r = d - h;
    if r.abs() < 0.1 {
        let mut ra = rho2;
        let rho22 = rho2 * rho2;
        let mut t = rho2;
        let mut i = 2;
        while i <= 50 {
            t = t * rho22 / ((i * i + i) as f64);
            ra += t;
            if t < 1.0e-30 {
                break;
            }
            i += 2;
        }
        r = ra + ra;
    }
    b[0] = r / rho2;
    let mut i = 2;
    let mut d2 = (ix - i + is) / is;
    while d2 > 0 {
        if ir != 0 {
            let il = is - 1;
            for jj in 1..=il {
                let k = i + jj - 1;
                let ku = k as usize;
                if k % 2 == 0 {
                    b[ku - 1] = -((d + h - (k - 1) as f64 * b[ku - 2]) / rho2);
                } else {
                    b[ku - 1] = (r + (k - 1) as f64 * b[ku - 2]) / rho2;
                }
            }
        }
        let inn = i + is - 1;
        if inn - ix > 0 {
            return;
        }
        let iu = inn as usize;
        if inn % 2 != 0 {
            let mut tr = 1.0;
            b[iu - 1] = 2.0 * tr / inn as f64;
            for jj in 1..=500 {
                tr = tr * rho2.powf(2.0) / ((2 * jj * (2 * jj - 1)) as f64);
                if (tr / b[iu - 1]).abs() <= 1.0e-7 {
                    break;
                }
                b[iu - 1] += 2.0 * tr / ((inn + 2 * jj) as f64);
            }
        } else {
            let mut tr = rho2;
            b[iu - 1] = -(2.0 * tr / (inn + 1) as f64);
            for jj in 1..=500 {
                tr = tr * rho2.powf(2.0) / ((2 * jj * (2 * jj + 1)) as f64);
                if (tr / b[iu - 1]).abs() <= 1.0e-7 {
                    break;
                }
                b[iu - 1] -= 2.0 * tr / ((inn + 1 + 2 * jj) as f64);
            }
        }
        d2 -= 1;
        i += is;
    }
}

const BINCOE: [[f64; 8]; 8] = [
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
    [0.0, 0.0, 1.0, 3.0, 6.0, 10.0, 15.0, 21.0],
    [0.0, 0.0, 0.0, 1.0, 4.0, 10.0, 20.0, 35.0],
    [0.0, 0.0, 0.0, 0.0, 1.0, 5.0, 15.0, 35.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 6.0, 21.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 7.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
];

/// Overlap of two Slater orbitals for angular momentum projection `m1` (`lovlap`)
#[allow(clippy::too_many_arguments)]
fn lovlap(
    a: &[f64; 30],
    b: &[f64; 30],
    sk1: f64,
    sk2: f64,
    r: f64,
    l1: i32,
    l2: i32,
    m1: i32,
    n1: i32,
    n2: i32,
) -> f64 {
    let mut fact = [0.0f64; 25];
    fact[0] = 1.0;
    for i in 1..25 {
        fact[i] = fact[i - 1] * i as f64;
    }
    let f = |k: i32| fact[k as usize];
    let m2 = m1;
    let rhoa = r * sk1;
    let rhob = r * sk2;
    let mut rhoap = rhoa.powf(n1 as f64);
    rhoap *= rhoap;
    rhoap *= rhoa;
    let mut rhoab = rhob.powf(n2 as f64);
    rhoab *= rhoab;
    rhoab *= rhob;
    let rhopo = rhoap * rhoab;
    let terma = 0.5f64.powf((l1 + l2 + 1) as f64)
        * (((l1 + l1 + 1) * (l2 + l2 + 1)) as f64 * f(l1 - m1) * f(l2 - m1)
            / (f(n1 + n1) * f(n2 + n2) * f(l1 + m1) * f(l2 + m1))
            * rhopo)
            .sqrt();
    let jend = 1 + (l1 - m1) / 2;
    let kend = 1 + (l2 - m2) / 2;
    let ieb = m1 + 1;
    let mut strad = 0.0;
    for ju in 0..jend {
        let iab = n1 - l1 + ju + ju + 1;
        let icb = l1 - m1 - ju - ju + 1;
        let con1 = f(l1 + l1 - ju - ju) / (f(l1 - m1 - ju - ju) * f(ju) * f(l1 - ju));
        for ku in 0..kend {
            let mut con12 =
                con1 * f(l2 + l2 - ku - ku) / (f(l2 - m2 - ku - ku) * f(ku) * f(l2 - ku));
            let iev = ju + ku + l2;
            if 2 * (iev / 2) != iev {
                con12 = -con12;
            }
            let ibb = n2 - l2 + ku + ku + 1;
            let idb = l2 - m2 - ku - ku + 1;
            let mut value = 0.0;
            for i6 in 1..=ieb {
                for i5 in 1..=ieb {
                    let mut value1 = BINCOE[i6 as usize - 1][ieb as usize - 1]
                        * BINCOE[i5 as usize - 1][ieb as usize - 1];
                    let iev = i5 + i6;
                    if 2 * (iev / 2) != iev {
                        value1 = -value1;
                    }
                    for i4 in 1..=idb {
                        value1 = -value1;
                        let value2 = BINCOE[i4 as usize - 1][idb as usize - 1] * value1;
                        for i3 in 1..=icb {
                            let mut value3 = BINCOE[i3 as usize - 1][icb as usize - 1] * value2;
                            for i2 in 1..=ibb {
                                value3 = -value3;
                                let value4 = BINCOE[i2 as usize - 1][ibb as usize - 1] * value3;
                                for i1 in 1..=iab {
                                    let term = value4 * BINCOE[i1 as usize - 1][iab as usize - 1];
                                    let ir =
                                        i1 + i2 + ieb + ieb - i6 - i6 - i3 + idb - i4 + icb - 1;
                                    let ip = iab - i1 + ibb - i2 + ieb + ieb - i5 - i5 + icb - i3
                                        + idb
                                        - i4
                                        + 1;
                                    value += a[ip as usize - 1] * b[ir as usize - 1] * term;
                                }
                            }
                        }
                    }
                }
            }
            strad += value * con12;
        }
    }
    strad * terma
}

/// sigma, pi, delta and phi overlaps between an orbital of atom 1 and one of atom 2 along the bond axis (`mov`)
#[allow(clippy::too_many_arguments)]
fn mov(a1: &Atom, a2: &Atom, dist: f64, q1: i32, q2: i32, l1: i32, l2: i32) -> [f64; 4] {
    let max = q1 + q2;
    let nn = l1.min(l2);
    let mut res = [0.0f64; 4];
    let mut ang = [0.0f64; 4];
    let num_zeta1 = if (l1 == 2 && a1.coeff_d2 != 0.0) || (l1 == 3 && a1.coeff_f2 != 0.0) {
        2
    } else {
        1
    };
    let num_zeta2 = if (l2 == 2 && a2.coeff_d2 != 0.0) || (l2 == 3 && a2.coeff_f2 != 0.0) {
        2
    } else {
        1
    };
    let exp_of = |a: &Atom, l: i32| match l {
        0 => a.exp_s,
        1 => a.exp_p,
        2 => a.exp_d,
        _ => a.exp_f,
    };
    let mut sk1 = exp_of(a1, l1);
    let mut sk2 = exp_of(a2, l2);
    let mut coeff_1 = if !(a1.coeff_d1 != 0.0 || a1.coeff_f1 != 0.0) || l1 <= 1 {
        1.0
    } else if l1 == 2 {
        a1.coeff_d1
    } else {
        a1.coeff_f1
    };
    let mut coeff_2 = if !(a2.coeff_d1 != 0.0 || a2.coeff_f1 != 0.0) || l2 <= 1 {
        1.0
    } else if l2 == 2 {
        a2.coeff_d1
    } else {
        a2.coeff_f1
    };
    let mut av = [0.0f64; 30];
    let mut bv = [0.0f64; 30];
    let mut accumulate = |sk1: f64, sk2: f64, c1: f64, c2: f64, res: &mut [f64; 4]| {
        abfns(&mut av, &mut bv, sk1, sk2, dist, max);
        for i in 0..=nn as usize {
            ang[i] = lovlap(&av, &bv, sk1, sk2, dist, l1, l2, i as i32, q1, q2);
        }
        for i in 0..4 {
            res[i] += c1 * c2 * ang[i];
        }
    };
    accumulate(sk1, sk2, coeff_1, coeff_2, &mut res);
    if num_zeta2 == 2 {
        if l2 == 2 {
            sk2 = a2.exp_d2;
            coeff_2 = a2.coeff_d2;
        }
        if l2 == 3 {
            sk2 = a2.exp_f2;
            coeff_2 = a2.coeff_f2;
        }
        accumulate(sk1, sk2, coeff_1, coeff_2, &mut res);
    }
    if num_zeta1 == 2 {
        if l1 == 2 {
            sk1 = a1.exp_d2;
            coeff_1 = a1.coeff_d2;
        }
        if l1 == 3 {
            sk1 = a1.exp_f2;
            coeff_1 = a1.coeff_f2;
        }
        accumulate(sk1, sk2, coeff_1, coeff_2, &mut res);
        if num_zeta2 == 2 {
            if l2 == 2 {
                sk2 = a2.exp_d;
                coeff_2 = a2.coeff_d1;
            }
            if l2 == 3 {
                sk2 = a2.exp_f;
                coeff_2 = a2.coeff_f1;
            }
            accumulate(sk1, sk2, coeff_1, coeff_2, &mut res);
        }
    }
    res
}

/// The number of orbitals of every atom and where they start (`build_orbital_lookup_table`)
pub fn orbital_lookup_table(atoms: &[Atom]) -> (Vec<i32>, usize) {
    let mut tab = Vec::with_capacity(atoms.len());
    let mut n = 0usize;
    for a in atoms {
        tab.push(n as i32);
        if a.ns != 0 {
            n += 1;
        }
        if a.np != 0 {
            n += 3;
        }
        if a.nd != 0 {
            n += 5;
        }
        if a.nf != 0 {
            n += 7;
        }
    }
    (tab, n)
}

/// The upper triangle of the overlap matrix of a molecule with unit diagonal (`R_space_overlap_matrix`
/// for a molecular calculation); atoms further apart than `rho` angstrom do not overlap.
#[allow(clippy::identity_op)]
pub fn overlap_matrix(
    atoms: &[Atom],
    pos: &[[f64; 3]],
    tab: &[i32],
    num_orbs: usize,
    rho: f64,
) -> Vec<f64> {
    let n = num_orbs;
    let mut overlap = vec![0.0f64; n * n];
    let mut p = [0.0f64; 9 + BP];
    let mut d = [0.0f64; 25 + BD];
    let mut f = [0.0f64; 49 + BF];
    for i in 0..atoms.len() {
        let i_tab = tab[i] as usize;
        for j in 0..i {
            let j_tab = tab[j] as usize;
            let dx = pos[i][0] - pos[j][0];
            let dy = pos[i][1] - pos[j][1];
            let dz = pos[i][2] - pos[j][2];
            let temp = dx * dx + dy * dy;
            let xy_dist = temp.sqrt();
            let mut tot_dist = (temp + dz * dz).sqrt();
            if !(tot_dist >= 1e-6 && tot_dist <= rho) {
                continue;
            }
            let (cos_b, sin_b, sin_a);
            if xy_dist < 1e-5 {
                cos_b = 1.0;
                sin_b = 0.0;
                sin_a = 0.0;
            } else {
                cos_b = dx / xy_dist;
                sin_b = dy / xy_dist;
                sin_a = xy_dist / tot_dist;
            }
            let cos_a = dz / tot_dist;
            let ai = &atoms[i];
            let aj = &atoms[j];
            let c2a = (cos_a * cos_a) - (sin_a * sin_a);
            let c2b = (cos_b * cos_b) - (sin_b * sin_b);
            let s2b = 2.0 * sin_b * cos_b;
            let sa3 = sin_a * sin_a * sin_a;
            let c3b = (c2b * cos_b) - (s2b * sin_b);
            let s3b = (c2b * sin_b) + (s2b * cos_b);
            let _ = (sa3, c3b, s3b, c2a);
            p[0 + BP] = sin_a * cos_b;
            p[1 + BP] = sin_a * sin_b;
            p[2 + BP] = cos_a;
            p[3 + BP] = cos_a * cos_b;
            p[4 + BP] = cos_a * sin_b;
            p[5 + BP] = -sin_a;
            p[6 + BP] = -sin_b;
            p[7 + BP] = cos_b;
            p[8 + BP] = 0.0;
            if ai.nd != 0 || aj.nd != 0 {
                d[0 + BD] = SQRT3 * 0.5 * sin_a * sin_a * c2b;
                d[1 + BD] = 1.0 - 1.5 * sin_a * sin_a;
                d[2 + BD] = SQRT3 * cos_b * sin_b * sin_a * sin_a;
                d[3 + BD] = SQRT3 * cos_a * sin_a * cos_b;
                d[4 + BD] = SQRT3 * cos_a * sin_a * sin_b;
                d[5 + BD] = cos_a * sin_a * c2b;
                d[6 + BD] = -SQRT3 * cos_a * sin_a;
                d[7 + BD] = cos_a * sin_a * s2b;
                d[8 + BD] = cos_b * c2a;
                d[9 + BD] = sin_b * c2a;
                d[10 + BD] = -sin_a * s2b;
                d[11 + BD] = 0.0;
                d[12 + BD] = sin_a * c2b;
                d[13 + BD] = -p[4 + BP];
                d[14 + BD] = p[3 + BP];
            }
            if (ai.nd != 0 || ai.nf != 0) && (aj.nd != 0 || aj.nf != 0) {
                d[15 + BD] = 0.50 * (1.0 + cos_a * cos_a) * c2b;
                d[16 + BD] = 0.50 * SQRT3 * sin_a * sin_a;
                d[17 + BD] = cos_b * sin_b * (1.0 + cos_a * cos_a);
                d[18 + BD] = -cos_a * sin_a * cos_b;
                d[19 + BD] = -cos_a * sin_a * sin_b;
                d[20 + BD] = -cos_a * s2b;
                d[21 + BD] = 0.0;
                d[22 + BD] = cos_a * c2b;
                d[23 + BD] = p[1 + BP];
                d[24 + BD] = -p[0 + BP];
            }
            if ai.nf != 0 || aj.nf != 0 {
                f[0 + BF] = 0.5 * cos_a * (5.0 * cos_a * cos_a - 3.0);
                f[1 + BF] = SQRT6 * 0.25 * cos_b * sin_a * (5.0 * cos_a * cos_a - 1.0);
                f[2 + BF] = SQRT6 * 0.25 * sin_b * sin_a * (5.0 * cos_a * cos_a - 1.0);
                f[3 + BF] = SQRT15 * sin_b * cos_b * cos_a * sin_a * sin_a;
                f[4 + BF] = SQRT15 * 0.5 * c2b * cos_a * sin_a * sin_a;
                f[5 + BF] = SQRT10 * 0.25 * c3b * sa3;
                f[6 + BF] = SQRT10 * 0.25 * s3b * sa3;
                f[7 + BF] = (-SQRT6) * 0.25 * sin_a * (5.0 * cos_a * cos_a - 1.0);
                f[8 + BF] = 0.25 * cos_b * cos_a * (15.0 * cos_a * cos_a - 11.0);
                f[9 + BF] = 0.25 * sin_b * cos_a * (15.0 * cos_a * cos_a - 11.0);
                f[10 + BF] = SQRT10 * 0.25 * s2b * sin_a * (3.0 * cos_a * cos_a - 1.0);
                f[11 + BF] = SQRT10 * 0.25 * c2b * sin_a * (3.0 * cos_a * cos_a - 1.0);
                f[12 + BF] = SQRT15 * 0.25 * c3b * cos_a * sin_a * sin_a;
                f[13 + BF] = SQRT15 * 0.25 * s3b * cos_a * sin_a * sin_a;
                f[14 + BF] = 0.0;
                f[15 + BF] = (-0.25) * sin_b * (5.0 * cos_a * cos_a - 1.0);
                f[16 + BF] = 0.25 * cos_b * (5.0 * cos_a * cos_a - 1.0);
                f[17 + BF] = SQRT10 * 0.5 * c2b * sin_a * cos_a;
                f[18 + BF] = (-SQRT10) * sin_a * cos_a * sin_b * cos_b;
                f[19 + BF] = (-SQRT15) * 0.25 * s3b * sin_a * sin_a;
                f[20 + BF] = SQRT15 * 0.25 * sin_a * sin_a * cos_b * (4.0 * cos_b * cos_b - 3.0);
            }
            if (ai.nd != 0 || ai.nf != 0) && (aj.nd != 0 || aj.nf != 0) {
                f[21 + BF] = 0.0;
                f[22 + BF] = SQRT10 * 0.5 * sin_b * cos_a * sin_a;
                f[23 + BF] = (-SQRT10) * 0.5 * cos_b * cos_a * sin_a;
                f[24 + BF] = c2b * (cos_a * cos_a - sin_a * sin_a);
                f[25 + BF] = (-s2b) * (cos_a * cos_a - sin_a * sin_a);
                f[26 + BF] = (-SQRT6) * 0.5 * s3b * sin_a * cos_a;
                f[27 + BF] = SQRT6 * 0.5 * cos_a * sin_a * cos_b * (4.0 * cos_b * cos_b - 3.0);
                f[28 + BF] = SQRT15 * 0.5 * cos_a * sin_a * sin_a;
                f[29 + BF] = SQRT10 * 0.25 * cos_b * sin_a * (1.0 - 3.0 * cos_a * cos_a);
                f[30 + BF] = SQRT10 * 0.25 * sin_a * sin_b * (1.0 - 3.0 * cos_a * cos_a);
                f[31 + BF] = sin_b * cos_b * cos_a * (3.0 * cos_a * cos_a - 1.0);
                f[32 + BF] = 0.5 * c2b * cos_a * (3.0 * cos_a * cos_a - 1.0);
                f[33 + BF] = SQRT6 * 0.25 * c3b * sin_a * (1.0 + cos_a * cos_a);
                f[34 + BF] = SQRT6 * 0.25 * s3b * sin_a * (1.0 + cos_a * cos_a);
            }
            if ai.nf != 0 && aj.nf != 0 {
                f[35 + BF] = SQRT10 * 0.25 * sin_a * (cos_a * cos_a - 1.0);
                f[36 + BF] = SQRT15 * 0.25 * cos_b * cos_a * sin_a * sin_a;
                f[37 + BF] = SQRT15 * 0.25 * sin_b * cos_a * sin_a * sin_a;
                f[38 + BF] = (-SQRT6) * 0.5 * sin_b * cos_b * sin_a * (1.0 + cos_a * cos_a);
                f[39 + BF] = (-SQRT6) * 0.25 * c2b * sin_a * (1.0 + cos_a * cos_a);
                f[40 + BF] = 0.25 * c3b * cos_a * (3.0 + cos_a * cos_a);
                f[41 + BF] = 0.25 * s3b * cos_a * (3.0 + cos_a * cos_a);
                f[42 + BF] = 0.0;
                f[43 + BF] = (-SQRT15) * 0.25 * sin_b * sin_a * sin_a;
                f[44 + BF] = SQRT15 * 0.25 * cos_b * sin_a * sin_a;
                f[45 + BF] = (-SQRT6) * 0.5 * c2b * sin_a * cos_a;
                f[46 + BF] = SQRT6 * sin_b * cos_b * sin_a * cos_a;
                f[47 + BF] = (-0.25) * s3b * (1.0 + 3.0 * cos_a * cos_a);
                f[48 + BF] = 0.25 * c3b * (1.0 + 3.0 * cos_a * cos_a);
            }
            tot_dist *= AUI;
            // < S(i) | S(j) >
            let (q1, q2) = (ai.ns, aj.ns);
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 0, 0);
                overlap[i_tab * n + j_tab] = m[0];
            }
            // < P(i) | S(j) >
            let q1 = ai.np;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 1, 0);
                let sigma = -m[0];
                for jo in BP..=EP {
                    overlap[(i_tab + jo) * n + j_tab] = p[jo] * sigma;
                }
            }
            // < P(i) | P(j) >
            let q2 = aj.np;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 1, 1);
                let (sigma, pi) = (-m[0], m[1]);
                for jo in BP..=EP {
                    for io in jo..=EP {
                        let v = p[jo] * p[io] * sigma
                            + (p[io + 3] * p[jo + 3] + p[io + 6] * p[jo + 6]) * pi;
                        overlap[(i_tab + io) * n + j_tab + jo] = v;
                        overlap[(i_tab + jo) * n + j_tab + io] = v;
                    }
                }
            }
            // < S(i) | P(j) >
            let q1 = ai.ns;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 0, 1);
                for jo in BP..=EP {
                    overlap[i_tab * n + j_tab + jo] = p[jo] * m[0];
                }
            }
            // < S(i) | D(j) >
            let q2 = aj.nd;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 0, 2);
                for jo in BD..=ED {
                    overlap[i_tab * n + j_tab + jo] = d[jo] * m[0];
                }
            }
            // < P(i) | D(j) >
            let q1 = ai.np;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 1, 2);
                let (sigma, pi) = (-m[0], m[1]);
                for io in BP..=EP {
                    for jo in BD..=ED {
                        overlap[(i_tab + io) * n + j_tab + jo] = p[io] * d[jo] * sigma
                            + (d[jo + 5] * p[io + 3] + d[jo + 10] * p[io + 6]) * pi;
                    }
                }
            }
            // < D(i) | S(j) >
            let q1 = ai.nd;
            let q2 = aj.ns;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 2, 0);
                for jo in BD..=ED {
                    overlap[(i_tab + jo) * n + j_tab] = d[jo] * m[0];
                }
            }
            // < D(i) | P(j) >
            let q2 = aj.np;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 2, 1);
                let (sigma, pi) = (m[0], -m[1]);
                for io in BP..=EP {
                    for jo in BD..=ED {
                        overlap[(i_tab + jo) * n + j_tab + io] = p[io] * d[jo] * sigma
                            + (p[io + 3] * d[jo + 5] + d[jo + 10] * p[io + 6]) * pi;
                    }
                }
            }
            // < D(i) | D(j) >
            let q2 = aj.nd;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 2, 2);
                let (sigma, pi, delta) = (m[0], -m[1], m[2]);
                for io in BD..=ED {
                    for jo in BD..=ED {
                        let v = d[io] * d[jo] * sigma
                            + (d[io + 5] * d[jo + 5] + d[jo + 10] * d[io + 10]) * pi
                            + (d[io + 15] * d[jo + 15] + d[io + 20] * d[jo + 20]) * delta;
                        overlap[(i_tab + jo) * n + (j_tab + io)] = v;
                        overlap[(i_tab + io) * n + j_tab + jo] = v;
                    }
                }
            }
            // < S(i) | F(j) >
            let q1 = ai.ns;
            let q2 = aj.nf;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 0, 3);
                for jo in BF..=EF {
                    overlap[i_tab * n + j_tab + jo] = f[jo] * m[0];
                }
            }
            // < P(i) | F(j) >
            let q1 = ai.np;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 1, 3);
                let (sigma, pi) = (-m[0], m[1]);
                for io in BP..=EP {
                    for jo in BF..=EF {
                        overlap[(i_tab + io) * n + j_tab + jo] = p[io] * f[jo] * sigma
                            + (p[io + 3] * f[jo + 7] + p[io + 6] * f[jo + 14]) * pi;
                    }
                }
            }
            // < D(i) | F(j) >
            let q1 = ai.nd;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 2, 3);
                let (sigma, pi, delta) = (m[0], -m[1], m[2]);
                for io in BD..=ED {
                    for jo in BF..=EF {
                        overlap[(i_tab + io) * n + j_tab + jo] = d[io] * f[jo] * sigma
                            + (d[io + 5] * f[jo + 7] + d[io + 10] * f[jo + 14]) * pi
                            + (d[io + 15] * f[jo + 28] + d[io + 20] * f[jo + 21]) * delta;
                    }
                }
            }
            // < F(i) | F(j) >
            let q1 = ai.nf;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 3, 3);
                let (sigma, pi, delta, phi) = (-m[0], m[1], m[2], m[3]);
                for io in BF..=EF {
                    for jo in BF..=EF {
                        let v = f[io] * f[jo] * sigma
                            + (f[io + 7] * f[jo + 7] + f[io + 14] * f[jo + 14]) * pi
                            + (f[io + 21] * f[jo + 21] + f[io + 28] * f[jo + 28]) * delta
                            + (f[io + 35] * f[jo + 35] + f[io + 42] * f[jo + 42]) * phi;
                        overlap[(i_tab + jo) * n + j_tab + io] = v;
                        overlap[(i_tab + io) * n + j_tab + jo] = v;
                    }
                }
            }
            // < F(i) | S(j) >
            let q2 = aj.ns;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 3, 0);
                let sigma = -m[0];
                for io in BF..=EF {
                    overlap[(i_tab + io) * n + j_tab] = f[io] * sigma;
                }
            }
            // < F(i) | P(j) >
            let q2 = aj.np;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 3, 1);
                let (sigma, pi) = (-m[0], m[1]);
                for io in BF..=EF {
                    for jo in BP..=EP {
                        overlap[(i_tab + io) * n + j_tab + jo] = f[io] * p[jo] * sigma
                            + (f[io + 7] * p[jo + 3] + f[io + 14] * p[jo + 6]) * pi;
                    }
                }
            }
            // < F(i) | D(j) >
            let q2 = aj.nd;
            if q1 != 0 && q2 != 0 {
                let m = mov(ai, aj, tot_dist, q1, q2, 3, 2);
                let (sigma, pi, delta) = (-m[0], m[1], -m[2]);
                for io in BF..=EF {
                    for jo in BD..=ED {
                        overlap[(i_tab + io) * n + j_tab + jo] = f[io] * d[jo] * sigma
                            + (f[io + 7] * d[jo + 5] + f[io + 14] * d[jo + 10]) * pi
                            + (f[io + 21] * d[jo + 20] + f[io + 28] * d[jo + 15]) * delta;
                    }
                }
            }
        }
    }
    // the matrix above is the lower triangle: move it to the upper one and put in the unit diagonal
    for j in 0..n {
        overlap[j * n + j] = 1.0;
        for k in 0..j {
            overlap[k * n + j] = overlap[j * n + k];
            overlap[j * n + k] = 0.0;
        }
    }
    overlap
}
