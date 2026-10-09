//! FPB reader tests on files copied from RDKit; expectations from RDKit's
//! own `testFPB.py` / `testFPB.cpp`.

use std::path::PathBuf;
use std::sync::{Arc, RwLock};

use rdoxide::datastructs::BitVect;
use rdoxide::datastructs::fpb_reader::FpbReader;
use rdoxide::datastructs::multi_fpb_reader::MultiFpbReader;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data/fpb")
        .join(name)
}

fn reader(name: &str, lazy: bool) -> FpbReader {
    let mut r = FpbReader::open(data(name), lazy).unwrap();
    r.init().unwrap();
    r
}

#[test]
fn basics_are_the_same_in_memory_and_lazily() {
    for lazy in [false, true] {
        let r = reader("zim.head100.fpb", lazy);
        assert_eq!((r.len().unwrap(), r.n_bits().unwrap()), (100, 2048));
        assert_eq!(r.get_id(0).unwrap(), "ZINC00902219");
        assert_eq!(r.get_id(3).unwrap(), "ZINC04803506");
        let fp = r.get_fp(0).unwrap();
        assert_eq!(fp.num_on_bits(), 17);
        assert_eq!(
            fp.on_bits(),
            [
                1, 80, 183, 222, 227, 231, 482, 650, 807, 811, 831, 888, 1335, 1411, 1664, 1820,
                1917
            ]
        );
        assert!(r.get_id(100).is_err() && r.get_fp(100).is_err());
    }
}

#[test]
fn similarity_searches() {
    for lazy in [false, true] {
        let r = reader("zim.head100.fpb", lazy);
        let q = r.get_bytes(0).unwrap();
        assert!((r.get_tanimoto(1, &q).unwrap() - 0.3704).abs() < 1e-4);
        let hits = r.get_tanimoto_neighbors(&q, 0.3, true).unwrap();
        assert_eq!(hits.iter().map(|h| h.1).take(2).collect::<Vec<_>>(), [0, 1]);
        assert_eq!(hits.len(), 5);
        // the popcount screen only prunes, it never changes the answer
        assert_eq!(hits, r.get_tanimoto_neighbors(&q, 0.3, false).unwrap());
        assert_eq!(
            hits,
            r.get_tversky_neighbors(&q, 1.0, 1.0, 0.3, true).unwrap()
        );
        assert_eq!(
            r.get_containing_neighbors(&r.get_bytes(1).unwrap())
                .unwrap(),
            [1, 2, 3, 4]
        );
    }
}

#[test]
fn errors() {
    assert!(FpbReader::open(data("missing.fpb"), false).is_err());
    let r = FpbReader::open(data("zim.head100.fpb"), false).unwrap();
    assert!(r.len().is_err(), "not initialized");
    let r = reader("zim.head100.fpb", false);
    assert!(
        r.get_tanimoto_neighbors(&r.get_bytes(0).unwrap(), 1.5, true)
            .is_err()
    );
    assert!(
        r.get_tanimoto_neighbors(&[0u8; 4], 0.5, true).is_err(),
        "short query"
    );
}

#[test]
fn multi_reader_threads_do_not_change_the_result() {
    let mut multi = MultiFpbReader::new(false);
    for name in ["zinc_random200.1.patt.fpb", "zinc_random200.2.patt.fpb"] {
        let r = FpbReader::open(data(name), false).unwrap();
        multi.add_reader(Arc::new(RwLock::new(r))).unwrap();
    }
    multi.init().unwrap();
    assert_eq!((multi.len(), multi.n_bits().unwrap()), (2, 1024));
    let q = reader("zinc_random200.2.patt.fpb", false)
        .get_bytes(3)
        .unwrap();
    let serial = multi.get_tanimoto_neighbors(&q, 0.5, 1).unwrap();
    assert!(serial.iter().any(|h| h.1 == 3 && h.2 == 1 && h.0 == 1.0));
    for threads in [2, 4, 0] {
        assert_eq!(
            serial,
            multi.get_tanimoto_neighbors(&q, 0.5, threads).unwrap()
        );
    }
    let containing = multi.get_containing_neighbors(&q, 2).unwrap();
    assert!(containing.contains(&(3, 1)));
}

#[test]
fn init_on_search_defers_reading() {
    let mut multi = MultiFpbReader::new(true);
    let r = FpbReader::open(data("zim.head100.fpb"), true).unwrap();
    multi.add_reader(Arc::new(RwLock::new(r))).unwrap();
    let q = reader("zim.head100.fpb", false).get_bytes(0).unwrap();
    assert_eq!(multi.get_tanimoto_neighbors(&q, 0.99, 1).unwrap().len(), 1);
}
