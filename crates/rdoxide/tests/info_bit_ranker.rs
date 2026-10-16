//! `InfoBitRanker` and `BitCorrMatGenerator`; the expectations are the
//! examples of RDKit 2026.09.1's `InfoBitRanker.h` and `testRanker.py`.

use rdoxide::datastructs::bit_vect::{BitVect, ExplicitBitVect, SparseBitVect};
use rdoxide::ml::info_bit_ranker::{BitCorrMatGenerator, InfoBitRanker, InfoType};
use rdoxide::rdgeneral::cxx_sort::LibcxxPriorityQueue;

/// 0001 0 / 0101 0 / 0010 1 / 1110 1
fn examples() -> Vec<(ExplicitBitVect, u32)> {
    [
        ([3].as_slice(), 0),
        ([1, 3].as_slice(), 0),
        ([2].as_slice(), 1),
        ([0, 1, 2].as_slice(), 1),
    ]
    .into_iter()
    .map(|(bits, act)| {
        let mut bv = ExplicitBitVect::new(4, false);
        for &b in bits {
            bv.set_bit(b).unwrap();
        }
        (bv, act)
    })
    .collect()
}

fn ranked(info: InfoType, bias: Option<&[i32]>) -> Vec<f64> {
    let mut r = InfoBitRanker::new(4, 2, info);
    if let Some(b) = bias {
        r.set_bias_list(b).unwrap();
    }
    for (bv, act) in examples() {
        r.accumulate_votes(&bv, act).unwrap();
    }
    r.get_top_n(3).unwrap().to_vec()
}

fn check(res: &[f64], expected: [(f64, f64, f64, f64); 3]) {
    assert_eq!(res.len(), 12);
    for (row, e) in res.chunks(4).zip(expected) {
        assert_eq!(row[0], e.0);
        assert!((row[1] - e.1).abs() < 5e-4, "{row:?} vs {e:?}");
        assert_eq!((row[2], row[3]), (e.2, e.3));
    }
}

#[test]
fn the_four_measures() {
    check(
        &ranked(InfoType::Entropy, None),
        [
            (3.0, 1.0, 2.0, 0.0),
            (2.0, 1.0, 0.0, 2.0),
            (0.0, 0.311, 0.0, 1.0),
        ],
    );
    check(
        &ranked(InfoType::BiasEntropy, Some(&[1])),
        [
            (2.0, 1.0, 0.0, 2.0),
            (0.0, 0.311, 0.0, 1.0),
            (1.0, 0.0, 1.0, 1.0),
        ],
    );
    check(
        &ranked(InfoType::ChiSquare, None),
        [
            (3.0, 4.0, 2.0, 0.0),
            (2.0, 4.0, 0.0, 2.0),
            (0.0, 1.333, 0.0, 1.0),
        ],
    );
    check(
        &ranked(InfoType::BiasChiSquare, Some(&[1])),
        [
            (2.0, 4.0, 0.0, 2.0),
            (0.0, 1.333, 0.0, 1.0),
            (1.0, 0.0, 1.0, 1.0),
        ],
    );
}

#[test]
fn sparse_and_explicit_vectors_vote_alike() {
    let mut a = InfoBitRanker::new(4, 2, InfoType::Entropy);
    let mut b = InfoBitRanker::new(4, 2, InfoType::Entropy);
    for (bv, act) in examples() {
        let mut sv = SparseBitVect::new(4);
        for bit in bv.on_bits() {
            sv.set_bit(bit as u32).unwrap();
        }
        a.accumulate_votes(&bv, act).unwrap();
        b.accumulate_votes(&sv, act).unwrap();
    }
    assert_eq!(a.num_instances(), 4);
    assert_eq!(a.get_top_n(4).unwrap(), b.get_top_n(4).unwrap());
}

#[test]
fn mask_bits_limit_the_ranking() {
    let mut r = InfoBitRanker::new(4, 2, InfoType::Entropy);
    r.set_mask_bits(&[0, 2]).unwrap();
    for (bv, act) in examples() {
        r.accumulate_votes(&bv, act).unwrap();
    }
    let mut ids: Vec<i32> = r
        .get_top_n(2)
        .unwrap()
        .chunks(4)
        .map(|row| row[0] as i32)
        .collect();
    ids.sort_unstable();
    assert_eq!(ids, [0, 2]);
    assert!(r.get_top_n(3).is_err());
}

#[test]
fn bad_arguments_are_errors() {
    let mut r = InfoBitRanker::new(4, 2, InfoType::Entropy);
    assert!(r.get_top_n(5).is_err());
    assert!(
        r.accumulate_votes(&ExplicitBitVect::new(4, false), 2)
            .is_err()
    );
    assert!(
        r.accumulate_votes(&ExplicitBitVect::new(5, false), 0)
            .is_err()
    );
    assert!(r.set_bias_list(&[0, 0]).is_err());
    assert!(r.set_bias_list(&[3]).is_err());
}

#[test]
fn nothing_ranked_gives_zero_info() {
    // test3Issue140: with every bit equally informative the ranking works
    let mut r = InfoBitRanker::new(2, 2, InfoType::Entropy);
    for (bits, act) in [([0u32, 1].as_slice(), 0), ([].as_slice(), 1)] {
        let mut bv = ExplicitBitVect::new(2, false);
        for &b in bits {
            bv.set_bit(b).unwrap();
        }
        r.accumulate_votes(&bv, act).unwrap();
    }
    assert_eq!(r.get_top_n(1).unwrap().len(), 4);
}

#[test]
fn table_of_the_top_bits() {
    let mut r = InfoBitRanker::new(4, 2, InfoType::Entropy);
    for (bv, act) in examples() {
        r.accumulate_votes(&bv, act).unwrap();
    }
    r.get_top_n(2).unwrap();
    assert_eq!(
        r.top_bits_table(),
        "         Bit InfoContent     class0     class1\n           3           1         2         0\n           2           1         0         2\n"
    );
}

#[test]
fn correlation_counts() {
    let mut g = BitCorrMatGenerator::new();
    g.set_bit_id_list(&[0, 1, 2]);
    for (bv, _) in examples() {
        g.collect_votes(&bv).unwrap();
    }
    assert_eq!(g.num_examples(), 4);
    // pairs (1,0), (2,0), (2,1): only the last example has 0 on
    assert_eq!(g.corr_mat(), [1.0, 1.0, 1.0]);
}

#[test]
fn libcxx_priority_queue_pops_in_order() {
    let mut q = LibcxxPriorityQueue::new(|a: &i32, b: &i32| a > b);
    for v in [5, 1, 4, 1, 3, 9, 2, 6] {
        q.push(v);
    }
    let mut out = Vec::new();
    while let Some(&t) = q.top() {
        out.push(t);
        q.pop();
    }
    assert_eq!(out, [1, 1, 2, 3, 4, 5, 6, 9]);
}
