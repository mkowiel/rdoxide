//! Thread helpers: `getNumThreadsToUse` of `Code/RDGeneral/RDThreads.h` and an
//! order-preserving parallel map for the batch APIs (RDKit runs these loops
//! on `std::async` tasks; scoped threads give the same behaviour without a
//! dependency).

use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

/// Name of the environment variable that caps the threads rdoxide uses.
pub const NUM_THREADS_ENV: &str = "RDOXIDE_NUM_THREADS";

/// 0: not resolved yet; otherwise the ceiling on worker threads.
static MAX_THREADS: AtomicUsize = AtomicUsize::new(0);

fn hardware_threads() -> usize {
    thread::available_parallelism().map_or(1, |n| n.get())
}

/// Interprets a thread setting: a positive value is that many threads, zero
/// (or "all", or an empty string) all hardware threads and a negative value
/// all but that many. `None` for text that is not a number.
pub fn parse_thread_setting(text: &str) -> Option<usize> {
    let text = text.trim();
    if text.is_empty() || text.eq_ignore_ascii_case("all") {
        return Some(hardware_threads());
    }
    let n: i64 = text.parse().ok()?;
    Some(resolve_setting(n))
}

fn resolve_setting(n: i64) -> usize {
    let hardware = hardware_threads();
    if n > 0 {
        n as usize
    } else {
        hardware.saturating_sub(n.unsigned_abs() as usize).max(1)
    }
}

/// The ceiling on the threads any rdoxide operation uses: the value of the
/// `RDOXIDE_NUM_THREADS` environment variable, or all hardware threads when it
/// is unset. [`set_max_threads`] overrides it at run time.
pub fn max_threads() -> usize {
    let current = MAX_THREADS.load(Ordering::Relaxed);
    if current != 0 {
        return current;
    }
    let resolved = std::env::var(NUM_THREADS_ENV)
        .ok()
        .and_then(|v| parse_thread_setting(&v))
        .unwrap_or_else(hardware_threads);
    MAX_THREADS.store(resolved, Ordering::Relaxed);
    resolved
}

/// Sets the ceiling on worker threads: positive is that many threads, zero all
/// hardware threads and negative all but that many.
pub fn set_max_threads(n: i32) {
    MAX_THREADS.store(resolve_setting(i64::from(n)), Ordering::Relaxed);
}

/// `getNumThreadsToUse` with a global ceiling: a positive target asks for that
/// many threads, zero asks for the default (the ceiling) and a negative target
/// for the ceiling minus `|target|`. The result is at least one and never
/// above [`max_threads`].
pub fn num_threads_to_use(target: i32) -> usize {
    let ceiling = max_threads();
    let wanted = if target >= 1 {
        target as usize
    } else {
        ceiling
            .saturating_sub(target.unsigned_abs() as usize)
            .max(1)
    };
    wanted.min(ceiling)
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

/// Runs `f(first_row, rows)` on contiguous blocks of whole rows of the
/// row-major matrix `data` (rows of `row_len` items), one block per worker
/// thread. Fewer than `min_rows_per_worker` rows per worker are not worth a
/// thread, so small matrices run on the calling thread.
pub fn par_row_blocks_mut<T, F>(data: &mut [T], row_len: usize, min_rows_per_worker: usize, f: F)
where
    T: Send,
    F: Fn(usize, &mut [T]) + Sync,
{
    if row_len == 0 || data.is_empty() {
        return;
    }
    let rows = data.len() / row_len;
    let workers = num_threads_to_use(0).min(rows / min_rows_per_worker.max(1));
    if workers <= 1 {
        f(0, data);
        return;
    }
    let rows_per_block = rows.div_ceil(workers);
    thread::scope(|scope| {
        let f = &f;
        let handles: Vec<_> = data
            .chunks_mut(rows_per_block * row_len)
            .enumerate()
            .map(|(b, block)| scope.spawn(move || f(b * rows_per_block, block)))
            .collect();
        for handle in handles {
            if let Err(payload) = handle.join() {
                std::panic::resume_unwind(payload);
            }
        }
    });
}

/// [`par_row_blocks_mut`] over two matrices with the same shape.
pub fn par_row_blocks_mut2<A, B, F>(
    a: &mut [A],
    b: &mut [B],
    row_len: usize,
    min_rows_per_worker: usize,
    f: F,
) where
    A: Send,
    B: Send,
    F: Fn(usize, &mut [A], &mut [B]) + Sync,
{
    if row_len == 0 || a.is_empty() {
        return;
    }
    debug_assert_eq!(a.len(), b.len());
    let rows = a.len() / row_len;
    let workers = num_threads_to_use(0).min(rows / min_rows_per_worker.max(1));
    if workers <= 1 {
        f(0, a, b);
        return;
    }
    let rows_per_block = rows.div_ceil(workers);
    thread::scope(|scope| {
        let f = &f;
        let handles: Vec<_> = a
            .chunks_mut(rows_per_block * row_len)
            .zip(b.chunks_mut(rows_per_block * row_len))
            .enumerate()
            .map(|(i, (ba, bb))| scope.spawn(move || f(i * rows_per_block, ba, bb)))
            .collect();
        for handle in handles {
            if let Err(payload) = handle.join() {
                std::panic::resume_unwind(payload);
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    /// The tests that change the process-wide ceiling must not overlap.
    static CEILING: Mutex<()> = Mutex::new(());

    #[test]
    fn targets_are_capped_by_the_ceiling() {
        let _guard = CEILING.lock().unwrap_or_else(|p| p.into_inner());
        // the ceiling is process-global: set it explicitly for the checks
        set_max_threads(8);
        assert_eq!(num_threads_to_use(1), 1);
        assert_eq!(num_threads_to_use(4), 4);
        assert_eq!(num_threads_to_use(64), 8);
        assert_eq!(num_threads_to_use(0), 8);
        assert_eq!(num_threads_to_use(-3), 5);
        assert_eq!(num_threads_to_use(i32::MIN), 1);
        set_max_threads(2);
        assert_eq!(num_threads_to_use(0), 2);
        assert_eq!(num_threads_to_use(16), 2);
        set_max_threads(0);
        assert_eq!(max_threads(), hardware_threads());
    }

    #[test]
    fn thread_settings_parse() {
        let hw = hardware_threads();
        assert_eq!(parse_thread_setting("3"), Some(3));
        assert_eq!(parse_thread_setting(" all "), Some(hw));
        assert_eq!(parse_thread_setting(""), Some(hw));
        assert_eq!(parse_thread_setting("0"), Some(hw));
        assert_eq!(
            parse_thread_setting("-1"),
            Some(hw.saturating_sub(1).max(1))
        );
        assert_eq!(parse_thread_setting("lots"), None);
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
    fn row_blocks_cover_every_row_once() {
        let _guard = CEILING.lock().unwrap_or_else(|p| p.into_inner());
        for (rows, row_len) in [(1usize, 3usize), (10, 4), (257, 5), (1000, 1)] {
            for ceiling in [1, 3, 8] {
                set_max_threads(ceiling);
                let mut data = vec![0u32; rows * row_len];
                par_row_blocks_mut(&mut data, row_len, 1, |first, block| {
                    for (r, row) in block.chunks_mut(row_len).enumerate() {
                        row.fill((first + r) as u32 + 1);
                    }
                });
                let expected: Vec<u32> = (0..rows * row_len)
                    .map(|k| (k / row_len) as u32 + 1)
                    .collect();
                assert_eq!(data, expected, "rows={rows} ceiling={ceiling}");
                let (mut a, mut b) = (vec![0u32; rows * row_len], vec![0u8; rows * row_len]);
                par_row_blocks_mut2(&mut a, &mut b, row_len, 1, |first, ba, bb| {
                    for (r, (ra, rb)) in ba
                        .chunks_mut(row_len)
                        .zip(bb.chunks_mut(row_len))
                        .enumerate()
                    {
                        ra.fill((first + r) as u32);
                        rb.fill(1);
                    }
                });
                assert!(b.iter().all(|&x| x == 1));
                assert_eq!(a[(rows - 1) * row_len], (rows - 1) as u32);
            }
        }
        set_max_threads(0);
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
