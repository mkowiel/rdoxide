//! MaxMin and Leader pickers; expected picks from the reference RDKit.

use rdoxide::sim_div_pickers::{leader, max_min};

const POINTS: [f64; 10] = [0.0, 1.0, 2.5, 4.0, 4.2, 7.0, 9.0, 9.5, 10.0, 15.0];

fn line(i: u32, j: u32) -> rdoxide::Result<f64> {
    Ok((POINTS[i as usize] - POINTS[j as usize]).abs())
}

fn packed_matrix() -> Vec<f64> {
    (0..10)
        .flat_map(|i| (0..i).map(move |j| (POINTS[i] - POINTS[j]).abs()))
        .collect()
}

#[test]
fn max_min_from_a_distance_matrix() {
    let picks = max_min::pick(&packed_matrix(), 10, 4, &[0], 1).unwrap();
    assert_eq!(picks, [0, 9, 5, 3]);
    // seeded random first pick
    assert_eq!(
        max_min::pick(&packed_matrix(), 10, 4, &[], 7).unwrap(),
        [0, 9, 5, 3]
    );
}

#[test]
fn max_min_stops_at_the_threshold_and_reports_the_last_distance() {
    let res = max_min::lazy_pick(line, 10, 10, &[0], 1, 3.0).unwrap();
    assert_eq!(res.picks, [0, 9, 5]);
    assert_eq!(res.threshold, 7.0);
    // nothing to pick: no distance was evaluated
    assert_eq!(
        max_min::lazy_pick(line, 10, 1, &[0], 1, -1.0)
            .unwrap()
            .threshold,
        -1.0
    );
}

#[test]
fn argument_errors() {
    assert!(max_min::lazy_pick(line, 0, 0, &[], 1, -1.0).is_err());
    assert!(max_min::lazy_pick(line, 3, 5, &[], 1, -1.0).is_err());
    assert!(max_min::lazy_pick(line, 3, 2, &[7], 1, -1.0).is_err());
    assert!(
        max_min::pick(&[0.0], 3, 2, &[], 1).is_err(),
        "matrix too small"
    );
    assert!(leader::lazy_pick(line, 10, 0, &[11], 1.0, 1).is_err());
}

#[test]
fn leader_picks_items_outside_the_threshold() {
    assert_eq!(
        leader::lazy_pick(line, 10, 0, &[], 1.5, 1).unwrap(),
        [0, 2, 4, 5, 6, 9]
    );
    assert_eq!(
        leader::lazy_pick(line, 10, 2, &[9], 1.5, 1).unwrap(),
        [9, 0]
    );
}

#[test]
fn leader_gives_the_same_picks_on_any_number_of_threads() {
    // large enough to be split over threads
    let n = 40_000u32;
    let dist = |i: u32, j: u32| {
        let (a, b) = (
            f64::from(i * 7919 % 1000) / 1000.0,
            f64::from(j * 7919 % 1000) / 1000.0,
        );
        Ok((a - b).abs() + f64::from(i != j) * 1e-3)
    };
    let serial = leader::lazy_pick(dist, n, 0, &[], 0.05, 1).unwrap();
    assert!(serial.len() > 5);
    for threads in [2, 4, 0] {
        assert_eq!(
            serial,
            leader::lazy_pick(dist, n, 0, &[], 0.05, threads).unwrap()
        );
    }
}

mod hierarchical {
    use rdoxide::ml::cluster::murtagh::{Linkage, cluster_points, hierarchical_cluster};
    use rdoxide::sim_div_pickers::hierarchical;

    const LINE: [f64; 7] = [0.0, 1.0, 3.0, 7.0, 8.0, 15.0, 16.5];

    fn packed() -> Vec<f64> {
        (0..LINE.len())
            .flat_map(|i| (0..i).map(move |j| (LINE[i] - LINE[j]).abs()))
            .collect()
    }

    #[test]
    fn clusters_and_representatives() {
        let d = packed();
        for method in [Linkage::Ward, Linkage::Slink, Linkage::Upgma] {
            let clusters = hierarchical::cluster(&d, 7, 3, method).unwrap();
            assert_eq!(
                clusters,
                [vec![0, 1, 2], vec![3, 4], vec![5, 6]],
                "{method:?}"
            );
        }
        assert_eq!(
            hierarchical::pick(&d, 7, 3, Linkage::Ward).unwrap(),
            [1, 3, 5]
        );
        assert_eq!(
            hierarchical::pick(&d, 7, 3, Linkage::Slink).unwrap(),
            [0, 3, 5]
        );
    }

    #[test]
    fn argument_errors() {
        let d = packed();
        assert!(hierarchical::cluster(&d, 7, 8, Linkage::Ward).is_err());
        assert!(hierarchical::cluster(&d, 7, 0, Linkage::Ward).is_err());
        assert!(hierarchical::cluster(&d[..5], 7, 3, Linkage::Ward).is_err());
    }

    #[test]
    fn point_clustering_matches_the_reference() {
        let x = [0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 5.0, 5.0, 6.0, 5.0];
        let h = cluster_points(&x, 5, 2, Linkage::Slink, 1).unwrap();
        assert_eq!((h.ia, h.ib), (vec![1, 4, 1, 1, 0], vec![2, 5, 3, 4, 0]));
        assert_eq!(h.crit, [1.0, 1.0, 4.0, 34.0, 0.0]);
    }

    #[test]
    fn the_clustering_overwrites_the_matrix_like_upstream() {
        let mut d = [1.0, 4.0, 3.0, 9.0, 8.0, 2.0, 7.0, 6.0, 5.0, 10.0];
        let h = hierarchical_cluster(5, Linkage::Ward, &mut d).unwrap();
        assert_eq!(h.ia, [1, 3, 1, 1, 0]);
        assert_eq!(h.ib, [2, 4, 5, 3, 0]);
        assert_eq!(d[0], 0.5);
        assert_eq!(d[1], 5.333333333333333);
    }
}
