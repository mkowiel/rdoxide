//! Extended Hueckel parameters (`fill_atomic_parms` of YAeHMOP's `fileio.c`, reading the table that is
//! compiled into the library).

use super::params_data::ROWS;

pub(super) struct Row {
    pub sym: &'static str,
    pub atnum: i32,
    pub nval: i32,
    pub nzeta: i32,
    pub nquant: i32,
    pub ang: u8,
    pub hii: f64,
    pub exp1: f64,
    pub exp2: f64,
    pub c1: f64,
    pub c2: f64,
}

/// The basis set of one atom
#[derive(Debug, Clone, Default)]
pub struct Atom {
    pub at_number: i32,
    pub num_valence: i32,
    /// principal quantum numbers of the s, p, d and f functions (0: not used)
    pub ns: i32,
    pub np: i32,
    pub nd: i32,
    pub nf: i32,
    pub coul_s: f64,
    pub coul_p: f64,
    pub coul_d: f64,
    pub coul_f: f64,
    pub exp_s: f64,
    pub exp_p: f64,
    pub exp_d: f64,
    pub exp_f: f64,
    pub exp_d2: f64,
    pub exp_f2: f64,
    pub coeff_d1: f64,
    pub coeff_d2: f64,
    pub coeff_f1: f64,
    pub coeff_f2: f64,
}

/// Looks up the parameters of an element by its (any case) symbol; `None` when the table has no such element.
pub fn atom_parameters(symbol: &str) -> Option<Atom> {
    let sym = symbol.to_ascii_uppercase();
    let first = ROWS.iter().position(|r| r.sym == sym)?;
    let mut atom = Atom::default();
    for r in ROWS[first..].iter().take_while(|r| r.sym == sym) {
        atom.at_number = r.atnum;
        atom.num_valence = r.nval;
        let n = if r.hii != 0.0 { r.nquant } else { 0 };
        match r.ang.to_ascii_lowercase() {
            b's' => {
                atom.ns = n;
                atom.exp_s = r.exp1;
                atom.coul_s = r.hii;
            }
            b'p' => {
                atom.np = n;
                atom.exp_p = r.exp1;
                atom.coul_p = r.hii;
            }
            b'd' | b'f' => {
                let is_d = r.ang.eq_ignore_ascii_case(&b'd');
                let (mut c1, mut c2) = (r.c1, 0.0);
                let mut e2 = 0.0;
                if r.nzeta == 2 {
                    e2 = r.exp2;
                    c2 = r.c2;
                    // renormalise the double zeta contraction
                    let mut temp = 4.0 * (r.exp1 * r.exp2 / (r.exp1 + r.exp2).powf(2.0));
                    temp = temp.powf(f64::from(r.nquant) + 0.5);
                    temp = (r.c1 * r.c1 + r.c2 * r.c2 + 2.0 * temp * r.c1 * r.c2).sqrt();
                    temp = 1.0 / temp;
                    c1 *= temp;
                    c2 *= temp;
                }
                if is_d {
                    atom.nd = n;
                    atom.exp_d = r.exp1;
                    atom.coul_d = r.hii;
                    atom.coeff_d1 = c1;
                    if r.nzeta == 2 {
                        atom.exp_d2 = e2;
                        atom.coeff_d2 = c2;
                    }
                } else {
                    atom.nf = n;
                    atom.exp_f = r.exp1;
                    atom.coul_f = r.hii;
                    atom.coeff_f1 = c1;
                    if r.nzeta == 2 {
                        atom.exp_f2 = e2;
                        atom.coeff_f2 = c2;
                    }
                }
            }
            _ => {}
        }
    }
    Some(atom)
}
