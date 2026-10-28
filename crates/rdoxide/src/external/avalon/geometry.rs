//! Port of Avalon Toolkit's `geometry.c`: 2D point set transformations and the helpers that
//! place substituents next to existing graphs. The C code's quirks (a second neighbour index
//! that always names the seed, a substituent score that uses the x coordinate twice) are kept.

/// `USE_INWARDS`, `H_GEOMETRY` and `NO_RAY_CROSSING` of `geometry.h`
pub const USE_INWARDS: i32 = 1;
pub const H_GEOMETRY: i32 = 2;
pub const NO_RAY_CROSSING: i32 = 4;

/// `TransformPoints`: rotates and moves `r` (without scaling) such that `p1`/`p2` come as close
/// as possible to `p1p`/`p2p`
pub fn transform_points(
    r: &mut [[f64; 2]],
    p1: [f64; 2],
    p2: [f64; 2],
    p1p: [f64; 2],
    p2p: [f64; 2],
) {
    let (mut ls, mut lsp) = (0.0, 0.0);
    let (mut dp, mut dpp, mut mp, mut mpp) = ([0.0; 2], [0.0; 2], [0.0; 2], [0.0; 2]);
    for i in 0..2 {
        dp[i] = p1[i] - p2[i];
        dpp[i] = p1p[i] - p2p[i];
        mp[i] = (p1[i] + p2[i]) / 2.0;
        mpp[i] = (p1p[i] + p2p[i]) / 2.0;
        ls += dp[i] * dp[i];
        lsp += dpp[i] * dpp[i];
    }
    let (a, b, c);
    if ls < 0.00001 || lsp < 0.00001 {
        a = 1.0;
        b = 0.0;
        c = 1.0;
    } else {
        a = (dp[0] * dpp[0] + dp[1] * dpp[1]) / ls;
        b = (dpp[0] * dp[1] - dp[0] * dpp[1]) / ls;
        c = (lsp / ls).sqrt();
    }
    for p in r.iter_mut() {
        let rp = [p[0] - mp[0], p[1] - mp[1]];
        p[0] = a * rp[0] + b * rp[1];
        p[1] = -b * rp[0] + a * rp[1];
        let rp = [p[0] / c, p[1] / c];
        p[0] = rp[0] + mpp[0];
        p[1] = rp[1] + mpp[1];
    }
}

/// `PointSetMatchTransformation`: rotates and translates (and optionally reflects) `points` such
/// that `from[0..nmatch]` comes as close as possible to `to[0..nmatch]`
pub fn point_set_match_transformation(
    points: &mut [[f64; 2]],
    from: &[[f64; 2]],
    to: &[[f64; 2]],
    nmatch: usize,
    reflection: bool,
) {
    let flipx = if reflection { -1.0 } else { 1.0 };
    let (mut dx1, mut dy1, mut dx2, mut dy2) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..nmatch {
        dx1 += flipx * from[i][0];
        dy1 += from[i][1];
        dx2 += to[i][0];
        dy2 += to[i][1];
    }
    let n = nmatch as f64;
    dx1 /= n;
    dy1 /= n;
    dx2 /= n;
    dy2 /= n;
    let (mut sxxp, mut syyp, mut sxyp, mut sxpy) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..nmatch {
        sxxp += (flipx * from[i][0] - dx1) * (to[i][0] - dx2);
        syyp += (from[i][1] - dy1) * (to[i][1] - dy2);
        sxyp += (flipx * from[i][0] - dx1) * (to[i][1] - dy2);
        sxpy += (from[i][1] - dy1) * (to[i][0] - dx2);
    }
    let denominator = ((sxxp + syyp) * (sxxp + syyp) + (sxpy - sxyp) * (sxpy - sxyp)).sqrt();
    let (mut calpha, mut salpha);
    if denominator < 1.0e-5 {
        calpha = 1.0;
        salpha = 0.0;
    } else {
        calpha = (sxxp + syyp) / denominator;
        salpha = (sxpy - sxyp) / denominator;
    }
    if salpha.abs() < 1.0e-10 {
        salpha = 0.0;
    }
    if calpha.abs() < 1.0e-10 {
        calpha = 0.0;
    }
    for p in points.iter_mut() {
        let xh = flipx * p[0] - dx1;
        let yh = p[1] - dy1;
        let xhh = calpha * xh + salpha * yh;
        let yhh = -salpha * xh + calpha * yh;
        p[0] = xhh + dx2;
        p[1] = yhh + dy2;
    }
}

fn len2(t: [f64; 2]) -> f64 {
    (t[0] * t[0] + t[1] * t[1]).sqrt()
}

fn dist(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)).sqrt()
}

fn dist2(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)
}

/// `NextGraphNeighbour`: the smallest neighbour of `node` that is larger than `last_neighbour`
pub fn next_graph_neighbour(node: i32, edges: &[[u32; 2]], last_neighbour: i32) -> i32 {
    let mut maxnode = 0u32;
    for e in edges {
        maxnode = maxnode.max(e[0]).max(e[1]);
    }
    let mut result = maxnode as i64 + 1;
    for e in edges {
        if i64::from(e[0]) == i64::from(node as u32) && i64::from(e[1]) < result {
            result = i64::from(e[1]);
        }
        if i64::from(e[1]) == i64::from(node as u32) && i64::from(e[0]) < result {
            result = i64::from(e[0]);
        }
    }
    if result == maxnode as i64 + 1 || result <= i64::from(last_neighbour) {
        -1
    } else {
        result as i32
    }
}

/// `HasCusp`: true if the neighbours form a cusp or umbrella around `m`
pub fn has_cusp(m: [f64; 2], neigh: &[[f64; 2]]) -> bool {
    let nneigh = neigh.len();
    if nneigh <= 1 {
        return false;
    }
    if nneigh == 2 {
        return true;
    }
    let mut mincosa = 2.0;
    let mut tb = [0.0f64; 2];
    for i in 0..nneigh {
        let mut ti = [neigh[i][0] - m[0], neigh[i][1] - m[1]];
        let len = len2(ti);
        ti[0] /= len;
        ti[1] /= len;
        for nj in neigh.iter().skip(i + 1) {
            let mut tj = [nj[0] - m[0], nj[1] - m[1]];
            let len = len2(tj);
            tj[0] /= len;
            tj[1] /= len;
            let cosa = ti[0] * tj[0] + ti[1] * tj[1];
            if cosa < -0.999 {
                return true;
            }
            if cosa < mincosa {
                mincosa = cosa;
                tb[0] = (ti[0] + tj[0]) / 2.0;
                tb[1] = (ti[1] + tj[1]) / 2.0;
                let len = len2(tb);
                tb[0] /= len;
                tb[1] /= len;
            }
        }
    }
    for n in neigh {
        let mut ti = [n[0] - m[0], n[1] - m[1]];
        let len = len2(ti);
        ti[0] /= len;
        ti[1] /= len;
        let cosa = ti[0] * tb[0] + ti[1] * tb[1];
        if cosa < 0.0 {
            return false;
        }
    }
    true
}

/// `PointsIntoRing`: true if the vector `coordinates[seed] -> p1` crosses a non-adjacent edge
pub fn points_into_ring(
    seed: usize,
    p1: [f64; 2],
    coordinates: &[[f64; 2]],
    edges: &[[u32; 2]],
) -> bool {
    let p0 = coordinates[seed];
    for e in edges {
        if e[0] as usize == seed || e[1] as usize == seed {
            continue;
        }
        let pp0 = coordinates[e[0] as usize];
        let pp1 = coordinates[e[1] as usize];
        let denominator = (pp1[0] - pp0[0]) * (p0[1] - p1[1]) - (pp1[1] - pp0[1]) * (p0[0] - p1[0]);
        if denominator.abs() < 1.0e-7 {
            continue;
        }
        let mut lambda = (p0[0] - pp0[0]) * (p0[1] - p1[1]) - (p0[1] - pp0[1]) * (p0[0] - p1[0]);
        lambda /= denominator;
        if !(0.0..=1.0).contains(&lambda) {
            continue;
        }
        let mut mu = (pp1[0] - pp0[0]) * (p0[1] - pp0[1]) - (p0[0] - pp0[0]) * (pp1[1] - pp0[1]);
        mu /= denominator;
        if mu <= 0.0 {
            continue;
        }
        return true;
    }
    false
}

/// `NextSubstituentPoint`: the best point to add a substituent to node `seed` of the graph; `point`
/// is left unchanged when there is no candidate
pub fn next_substituent_point(
    point: &mut [f64; 2],
    coordinates: &[[f64; 2]],
    edges: &[[u32; 2]],
    seed: usize,
    flags: i32,
) {
    const NPHI: usize = 2;
    const SINPHI: [f64; 4] = [0.866025, -0.866025, 0.707107, -0.707107];
    const COSPHI: [f64; 4] = [-0.5, -0.5, 0.707107, 0.707107];
    let m = coordinates[seed];
    let nneigh = edges
        .iter()
        .filter(|e| e[0] as usize == seed || e[1] as usize == seed)
        .count();
    let mut cand: Vec<[f64; 2]> = Vec::new();
    let mut neigh: Vec<[f64; 2]> = Vec::new();
    let mut neigh_degree: Vec<i32> = Vec::new();
    let mut nphi_used = false;
    if nneigh == 0 {
        let mut th = [1.514, 0.0];
        if flags & H_GEOMETRY != 0 {
            th[0] *= 0.666;
        }
        cand.push([th[0] + m[0], th[1] + m[1]]);
        cand.push([-th[0] + m[0], -th[1] + m[1]]);
        for i in 0..2 * NPHI {
            cand.push([
                m[0] + th[0] * COSPHI[i] - th[1] * SINPHI[i],
                m[1] + th[0] * SINPHI[i] + th[1] * COSPHI[i],
            ]);
        }
        for i in 0..2 * NPHI {
            cand.push([
                m[0] - th[0] * COSPHI[i] + th[1] * SINPHI[i],
                m[1] - th[0] * SINPHI[i] - th[1] * COSPHI[i],
            ]);
        }
        nphi_used = true;
    } else {
        let mut neigh_index: Vec<u32> = Vec::new();
        for e in edges {
            if e[0] as usize == seed {
                neigh.push(coordinates[e[1] as usize]);
                neigh_index.push(e[1]);
            }
            if e[1] as usize == seed {
                neigh.push(coordinates[e[0] as usize]);
                // the C code records edges[i][1], i.e. the seed itself
                neigh_index.push(e[1]);
            }
        }
        neigh_degree.resize(neigh.len(), 0);
        for e in edges {
            for (j, &ni) in neigh_index.iter().enumerate() {
                if e[0] == ni {
                    neigh_degree[j] += 1;
                }
                if e[1] == ni {
                    neigh_degree[j] += 1;
                }
            }
        }
        let nn = neigh.len();
        let mut len = 0.0;
        for n in &neigh {
            len += dist(m[0], m[1], n[0], n[1]);
        }
        len /= nn as f64;
        let d0 = dist(m[0], m[1], neigh[0][0], neigh[0][1]);
        if d0 > 0.001 && nn <= 2 {
            let th = [
                (neigh[0][0] - m[0]) * len / d0,
                (neigh[0][1] - m[1]) * len / d0,
            ];
            for i in 0..NPHI {
                cand.push([
                    th[0] * COSPHI[i] - th[1] * SINPHI[i] + m[0],
                    th[0] * SINPHI[i] + th[1] * COSPHI[i] + m[1],
                ]);
            }
            nphi_used = true;
        }
        let cusp = has_cusp(m, &neigh);
        for i in 0..nn {
            for j in i + 1..nn {
                let mut t = [neigh[i][0] - m[0], neigh[i][1] - m[1]];
                let d = (t[0] * t[0] + t[1] * t[1]).sqrt();
                if d > 0.001 {
                    t[0] /= d;
                    t[1] /= d;
                } else {
                    continue;
                }
                let mut th = [neigh[j][0] - m[0], neigh[j][1] - m[1]];
                let d = (th[0] * th[0] + th[1] * th[1]).sqrt();
                if d > 0.001 {
                    th[0] /= d;
                    th[1] /= d;
                } else {
                    continue;
                }
                let cos_alpha = t[0] * th[0] + t[1] * th[1];
                t[0] += th[0];
                t[1] += th[1];
                let d = (t[0] * t[0] + t[1] * t[1]).sqrt();
                if d > 0.001 {
                    if cusp {
                        cand.push([m[0] - t[0] * len / d, m[1] - t[1] * len / d]);
                    }
                    if (flags & USE_INWARDS) != 0 && nn > 2 && cos_alpha > -0.9 {
                        cand.push([m[0] + t[0] * len / d, m[1] + t[1] * len / d]);
                    } else {
                        // the C code computes this point but never stores it
                    }
                }
            }
        }
    }
    let mut value = 1.0e7;
    for (i, c) in cand.iter().enumerate() {
        let cusp = has_cusp(m, &neigh);
        let mut into_ring = false;
        if !cusp {
            into_ring = points_into_ring(seed, *c, coordinates, edges);
        }
        let mut vtmp_local = 0.0;
        for (j, n) in neigh.iter().enumerate() {
            vtmp_local += f64::from(neigh_degree[j]) * 1.0 / (0.1 + dist2(c[0], c[1], n[0], n[1]));
        }
        let mut vtmp = 0.0;
        for co in coordinates {
            vtmp += 1.0 / (0.1 + dist2(c[0], c[1], co[0], co[1]));
        }
        if i < NPHI && nphi_used {
            vtmp *= 1.85;
        }
        if cusp {
            vtmp += 100.0 * vtmp_local;
        } else {
            vtmp += 10.0 * vtmp_local;
        }
        if into_ring {
            vtmp *= 10.0;
        }
        if vtmp < value {
            value = vtmp;
            *point = *c;
        }
    }
}

/// `TransformGraphSubstituent`: rotates and moves the substituent graph `scoords` such that,
/// linked to the graph `mcoords` between `apo` and `seed`, it does not overlap too much
pub fn transform_graph_substituent(
    scoords: &mut [[f64; 2]],
    sedges: &[[u32; 2]],
    apo: usize,
    mcoords: &[[f64; 2]],
    medges: &[[u32; 2]],
    seed: usize,
) {
    const NPHI: usize = 4;
    const SINPHI: [f64; 4] = [0.866025, -0.866025, 0.866025, -0.866025];
    const COSPHI: [f64; 4] = [-0.5, -0.5, 0.5, 0.5];
    let nsnodes = scoords.len();
    let nmnodes = mcoords.len();
    if nsnodes <= 1 || nmnodes <= 1 {
        for i in 0..nsnodes {
            scoords[i][0] -= scoords[apo][0] - mcoords[seed][0];
            scoords[i][1] -= scoords[apo][1] - mcoords[seed][1];
        }
        return;
    }
    let m1 = mcoords[seed];
    let h = next_graph_neighbour(seed as i32, medges, -1);
    let m2 = if h == -1 {
        [m1[0] + 1.5, m1[1]]
    } else {
        mcoords[h as usize]
    };
    let s1 = scoords[apo];
    let h = next_graph_neighbour(apo as i32, sedges, -1);
    let s2 = if h == -1 {
        [s1[0] + 1.5, s1[1]]
    } else {
        scoords[h as usize]
    };
    let mut igood = 0;
    let mut value = 1.0e7;
    let d = dist2(s1[0], s1[1], s2[0], s2[1]);
    let mut r = vec![[0.0; 2]; nsnodes];
    for i in 0..NPHI {
        let th = [s2[0] - s1[0], s2[1] - s1[1]];
        let t = [
            th[0] * COSPHI[i] - th[1] * SINPHI[i] + s1[0],
            th[0] * SINPHI[i] + th[1] * COSPHI[i] + s1[1],
        ];
        if (0..nsnodes)
            .any(|j| j != apo && d * 0.01 > dist2(t[0], t[1], scoords[j][0], scoords[j][1]))
        {
            continue;
        }
        r.copy_from_slice(scoords);
        transform_points(&mut r, s1, t, m1, m2);
        let mut vtmp = 0.0;
        for rj in &r {
            for mk in mcoords {
                vtmp += 1.0 / (0.1 + dist2(rj[0], rj[1], mk[0], mk[0]));
            }
        }
        if vtmp < value {
            igood = i;
            value = vtmp;
        }
    }
    let th = [s2[0] - s1[0], s2[1] - s1[1]];
    let t = [
        th[0] * COSPHI[igood] - th[1] * SINPHI[igood] + s1[0],
        th[0] * SINPHI[igood] + th[1] * COSPHI[igood] + s1[1],
    ];
    r.copy_from_slice(scoords);
    transform_points(&mut r, s1, t, m1, m2);
    scoords.copy_from_slice(&r);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::scan::fmt_gp;

    struct Tok<'a>(std::str::SplitWhitespace<'a>);

    impl Tok<'_> {
        fn f(&mut self) -> f64 {
            self.0.next().expect("token").parse().expect("float")
        }

        fn i(&mut self) -> i64 {
            self.0.next().expect("token").parse().expect("int")
        }

        fn pts(&mut self, n: usize) -> Vec<[f64; 2]> {
            (0..n).map(|_| [self.f(), self.f()]).collect()
        }

        fn edges(&mut self, n: usize) -> Vec<[u32; 2]> {
            (0..n).map(|_| [self.i() as u32, self.i() as u32]).collect()
        }
    }

    fn show(p: &[[f64; 2]]) -> String {
        p.iter()
            .map(|q| format!("{} {}", fmt_gp(q[0], 17), fmt_gp(q[1], 17)))
            .collect::<Vec<_>>()
            .join("\n")
    }

    /// every function must give bit-identical results to the C code
    #[test]
    fn geometry_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_geometry.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let cmd = v["cmd"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let mut t = Tok(cmd.split_whitespace());
            let name = t.0.next().unwrap();
            let got = match name {
                "TP" => {
                    let n = t.i() as usize;
                    let (p1, p2, q1, q2) = (t.pts(1)[0], t.pts(1)[0], t.pts(1)[0], t.pts(1)[0]);
                    let mut p = t.pts(n);
                    transform_points(&mut p, p1, p2, q1, q2);
                    show(&p)
                }
                "PSM" => {
                    let (n, nm, refl) = (t.i() as usize, t.i() as usize, t.i() != 0);
                    let mut p = t.pts(n);
                    let (f, to) = (t.pts(nm), t.pts(nm));
                    point_set_match_transformation(&mut p, &f, &to, nm, refl);
                    show(&p)
                }
                "NGN" => {
                    let (node, last, ne) = (t.i() as i32, t.i() as i32, t.i() as usize);
                    let e = t.edges(ne);
                    next_graph_neighbour(node, &e, last).to_string()
                }
                "CUSP" => {
                    let m = t.pts(1)[0];
                    let n = t.i() as usize;
                    let p = t.pts(n);
                    i32::from(has_cusp(m, &p)).to_string()
                }
                "PIR" => {
                    let seed = t.i() as usize;
                    let p1 = t.pts(1)[0];
                    let nn = t.i() as usize;
                    let c = t.pts(nn);
                    let ne = t.i() as usize;
                    let e = t.edges(ne);
                    i32::from(points_into_ring(seed, p1, &c, &e)).to_string()
                }
                "NSP" => {
                    let (seed, flags) = (t.i() as usize, t.i() as i32);
                    let mut pt = t.pts(1)[0];
                    let nn = t.i() as usize;
                    let c = t.pts(nn);
                    let ne = t.i() as usize;
                    let e = t.edges(ne);
                    next_substituent_point(&mut pt, &c, &e, seed, flags);
                    show(&[pt])
                }
                "TGS" => {
                    let ns = t.i() as usize;
                    let mut s = t.pts(ns);
                    let nse = t.i() as usize;
                    let se = t.edges(nse);
                    let apo = t.i() as usize;
                    let nm = t.i() as usize;
                    let m = t.pts(nm);
                    let nme = t.i() as usize;
                    let me = t.edges(nme);
                    let seed = t.i() as usize;
                    transform_graph_substituent(&mut s, &se, apo, &m, &me, seed);
                    show(&s)
                }
                other => panic!("unknown command {other}"),
            };
            if got != want {
                bad.push(format!("{cmd}\nwant {want}\ngot  {got}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }
}
