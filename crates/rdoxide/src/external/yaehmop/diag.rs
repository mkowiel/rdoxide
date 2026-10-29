//! Symmetric eigenproblems: the Cholesky reduction of the generalized problem and the Householder /
//! implicit QL solver (EISPACK `tred2`/`tql2`) YAeHMOP's `cboris` is built from.

/// All eigenvalues (ascending) and eigenvectors of the real symmetric matrix `v` (n x n, row-major,
/// overwritten with the eigenvectors stored as columns). `None` when the QL iteration does not converge.
fn symmetric_eigen(n: usize, v: &mut [f64]) -> Option<Vec<f64>> {
    let mut d = vec![0.0f64; n];
    let mut e = vec![0.0f64; n];
    let at = |i: usize, j: usize| i * n + j;
    for j in 0..n {
        d[j] = v[at(n - 1, j)];
    }
    // Householder reduction to tridiagonal form
    for i in (1..n).rev() {
        let mut scale = 0.0;
        let mut h = 0.0;
        for k in 0..i {
            scale += d[k].abs();
        }
        if scale == 0.0 {
            e[i] = d[i - 1];
            for j in 0..i {
                d[j] = v[at(i - 1, j)];
                v[at(i, j)] = 0.0;
                v[at(j, i)] = 0.0;
            }
        } else {
            for k in 0..i {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            let mut f = d[i - 1];
            let mut g = h.sqrt();
            if f > 0.0 {
                g = -g;
            }
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for j in 0..i {
                e[j] = 0.0;
            }
            for j in 0..i {
                f = d[j];
                v[at(j, i)] = f;
                g = e[j] + v[at(j, j)] * f;
                for k in j + 1..i {
                    g += v[at(k, j)] * d[k];
                    e[k] += v[at(k, j)] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for j in 0..i {
                e[j] /= h;
                f += e[j] * d[j];
            }
            let hh = f / (h + h);
            for j in 0..i {
                e[j] -= hh * d[j];
            }
            for j in 0..i {
                f = d[j];
                g = e[j];
                for k in j..i {
                    v[at(k, j)] -= f * e[k] + g * d[k];
                }
                d[j] = v[at(i - 1, j)];
                v[at(i, j)] = 0.0;
            }
        }
        d[i] = h;
    }
    // accumulate the transformations
    for i in 0..n.saturating_sub(1) {
        v[at(n - 1, i)] = v[at(i, i)];
        v[at(i, i)] = 1.0;
        let h = d[i + 1];
        if h != 0.0 {
            for k in 0..=i {
                d[k] = v[at(k, i + 1)] / h;
            }
            for j in 0..=i {
                let mut g = 0.0;
                for k in 0..=i {
                    g += v[at(k, i + 1)] * v[at(k, j)];
                }
                for k in 0..=i {
                    v[at(k, j)] -= g * d[k];
                }
            }
        }
        for k in 0..=i {
            v[at(k, i + 1)] = 0.0;
        }
    }
    for j in 0..n {
        d[j] = v[at(n - 1, j)];
        v[at(n - 1, j)] = 0.0;
    }
    v[at(n - 1, n - 1)] = 1.0;
    e[0] = 0.0;
    // implicit QL
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    let mut f = 0.0;
    let mut tst1: f64 = 0.0;
    let eps = 2.0f64.powi(-52);
    for l in 0..n {
        tst1 = tst1.max(d[l].abs() + e[l].abs());
        let mut m = l;
        while m < n {
            if e[m].abs() <= eps * tst1 {
                break;
            }
            m += 1;
        }
        if m > l {
            let mut iter = 0;
            loop {
                iter += 1;
                if iter > 200 {
                    return None;
                }
                let mut g = d[l];
                let mut p = (d[l + 1] - g) / (2.0 * e[l]);
                let mut r = p.hypot(1.0);
                if p < 0.0 {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                let dl1 = d[l + 1];
                let mut h = g - d[l];
                for i in l + 2..n {
                    d[i] -= h;
                }
                f += h;
                p = d[m];
                let mut c = 1.0;
                let mut c2 = c;
                let mut c3 = c;
                let el1 = e[l + 1];
                let mut s = 0.0;
                let mut s2 = 0.0;
                for i in (l..m).rev() {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = p.hypot(e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for k in 0..n {
                        h = v[at(k, i + 1)];
                        v[at(k, i + 1)] = s * v[at(k, i)] + c * h;
                        v[at(k, i)] = c * v[at(k, i)] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                if e[l].abs() <= eps * tst1 {
                    break;
                }
            }
        }
        d[l] += f;
        e[l] = 0.0;
    }
    // sort ascending
    for i in 0..n.saturating_sub(1) {
        let mut k = i;
        let mut p = d[i];
        for j in i + 1..n {
            if d[j] < p {
                k = j;
                p = d[j];
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            for j in 0..n {
                v.swap(at(j, i), at(j, k));
            }
        }
    }
    Some(d)
}

/// Solves `H x = e S x` for symmetric `H` and positive definite `S`, both given by their upper triangles
/// (row-major, n x n). Returns the eigenvalues in ascending order and the eigenvectors as rows
/// (`vect[mo * n + ao]`), normalised to `x^T S x = 1`. `None` when `S` is not positive definite or the QL
/// iteration fails.
pub fn generalized_eigen(
    n: usize,
    h_upper: &[f64],
    s_upper: &[f64],
) -> Option<(Vec<f64>, Vec<f64>)> {
    // Cholesky factor S = L L^T
    let mut l = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = s_upper[j * n + i];
            for k in 0..j {
                sum -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if sum <= 0.0 {
                    return None;
                }
                l[i * n + i] = sum.sqrt();
            } else {
                l[i * n + j] = sum / l[j * n + j];
            }
        }
    }
    // C = L^-1 H L^-T
    let mut a = vec![0.0f64; n * n];
    for i in 0..n {
        for j in i..n {
            let x = h_upper[i * n + j];
            a[i * n + j] = x;
            a[j * n + i] = x;
        }
    }
    // a <- L^-1 a
    for j in 0..n {
        for i in 0..n {
            let mut sum = a[i * n + j];
            for k in 0..i {
                sum -= l[i * n + k] * a[k * n + j];
            }
            a[i * n + j] = sum / l[i * n + i];
        }
    }
    // a <- a L^-T
    for i in 0..n {
        for j in 0..n {
            let mut sum = a[i * n + j];
            for k in 0..j {
                sum -= a[i * n + k] * l[j * n + k];
            }
            a[i * n + j] = sum / l[j * n + j];
        }
    }
    for i in 0..n {
        for j in 0..i {
            let m = 0.5 * (a[i * n + j] + a[j * n + i]);
            a[i * n + j] = m;
            a[j * n + i] = m;
        }
    }
    let vals = symmetric_eigen(n, &mut a)?;
    // x = L^-T y, one eigenvector per row of the result
    let mut vect = vec![0.0f64; n * n];
    for mo in 0..n {
        for ao in (0..n).rev() {
            let mut sum = a[ao * n + mo];
            for k in ao + 1..n {
                sum -= l[k * n + ao] * vect[mo * n + k];
            }
            vect[mo * n + ao] = sum / l[ao * n + ao];
        }
    }
    Some((vals, vect))
}
