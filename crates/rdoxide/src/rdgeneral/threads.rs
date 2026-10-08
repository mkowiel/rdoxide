//! Thread helpers: `getNumThreadsToUse` of `Code/RDGeneral/RDThreads.h` and an
//! order-preserving parallel map for the batch APIs (RDKit runs these loops
//! on `std::async` tasks; scoped threads give the same behaviour without a
//! dependency).

use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

/// `getNumThreadsToUse`: a positive target is used as is, zero or a negative
/// target means "all hardware threads minus `|target|`" (at least one).
pub fn num_threads_to_use(target: i32) -> usize {
    if target >= 1 {
        return target as usize;
    }
    let hardware = thread::available_parallelism().map_or(1, |n| n.get());
    let reserved = target.unsigned_abs() as usize;
    if hardware > reserved {
        hardware - reserved
    } else {
        1
    }
}

/// Maps `f` over `items` on up to `num_threads_to_use(target)` threads.
///
/// The result is in the order of `items` whatever the thread count. Items are
/// handed out one at a time, so uneven items balance across the threads. A
/// panic in `f` is propagated to the caller.
pub fn par_map<T, R, F>(items: &[T], target: i32, f: F) -> Vec<R>
where
    T: Sync,
    R: Send,
    F: Fn(&T) -> R + Sync,
{
    let workers = num_threads_to_use(target).min(items.len());
    if workers <= 1 {
        return items.iter().map(f).collect();
    }
    let next = AtomicUsize::new(0);
    let mut slots: Vec<Option<R>> = items.iter().map(|_| None).collect();
    thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                scope.spawn(|| {
                    let mut done = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(item) = items.get(i) else { break };
                        done.push((i, f(item)));
                    }
                    done
                })
            })
            .collect();
        for handle in handles {
            match handle.join() {
                Ok(done) => {
                    for (i, r) in done {
                        slots[i] = Some(r);
                    }
                }
                Err(payload) => std::panic::resume_unwind(payload),
            }
        }
    });
    slots
        .into_iter()
        .map(|r| r.expect("every index is claimed by exactly one worker"))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positive_targets_are_used_as_is() {
        assert_eq!(num_threads_to_use(1), 1);
        assert_eq!(num_threads_to_use(64), 64);
    }

    #[test]
    fn non_positive_targets_leave_at_least_one_thread() {
        let hardware = thread::available_parallelism().map_or(1, |n| n.get());
        assert_eq!(num_threads_to_use(0), hardware);
        assert_eq!(num_threads_to_use(-1), hardware.saturating_sub(1).max(1));
        assert_eq!(num_threads_to_use(i32::MIN), 1);
    }

    #[test]
    fn results_keep_the_input_order() {
        let items: Vec<u64> = (0..1000).collect();
        for target in [1, 2, 7, 64] {
            let out = par_map(&items, target, |&x| x * x);
            assert_eq!(out, items.iter().map(|&x| x * x).collect::<Vec<_>>());
        }
    }

    #[test]
    fn empty_and_single_inputs() {
        assert!(par_map::<u8, u8, _>(&[], 8, |&x| x).is_empty());
        assert_eq!(par_map(&[5u8], 8, |&x| x + 1), vec![6]);
    }

    #[test]
    fn every_item_is_processed_once() {
        let calls = AtomicUsize::new(0);
        let items = vec![(); 257];
        let out = par_map(&items, 8, |_| calls.fetch_add(1, Ordering::SeqCst));
        assert_eq!(calls.load(Ordering::SeqCst), 257);
        let mut seen = out;
        seen.sort_unstable();
        assert_eq!(seen, (0..257).collect::<Vec<_>>());
    }

    #[test]
    #[should_panic(expected = "boom")]
    fn panics_propagate() {
        par_map(&[1, 2, 3, 4], 4, |&x| {
            if x == 3 {
                panic!("boom");
            }
            x
        });
    }
}
