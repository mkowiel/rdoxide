//! The limit on the number of recursive comparisons (upstream
//! `CIPLabeler.cpp`: `IterationBudget`, `ScopedIterationBudget`,
//! `ScopedPreliminaryBudget`, `decrementRemainingCallCountAndCheck`).

use crate::RDKitError;
use std::cell::Cell;

/// The comparisons allowed for the constitutional rules alone
const CONSTITUTIONAL_RULE_TIMEOUT: u32 = 2000;

/// `MaxIterationsExceeded`
pub(super) const MAX_ITERATIONS_MESSAGE: &str = "Max Iterations Exceeded in CIP label calculation";

pub(super) fn max_iterations_exceeded() -> RDKitError {
    RDKitError::Runtime(MAX_ITERATIONS_MESSAGE.into())
}

pub(super) fn is_max_iterations_exceeded(e: &RDKitError) -> bool {
    matches!(e, RDKitError::Runtime(m) if m == MAX_ITERATIONS_MESSAGE)
}

#[derive(Debug, Clone, Copy, Default)]
struct IterationBudget {
    has_global_limit: bool,
    remaining_global: u32,
    in_preliminary_pass: bool,
    remaining_preliminary: u32,
}

thread_local! {
    static BUDGET: Cell<IterationBudget> = Cell::new(IterationBudget::default());
}

/// `decrementRemainingCallCountAndCheck`: false when the budget is exhausted.
pub(super) fn consume() -> bool {
    BUDGET.with(|cell| {
        let mut b = cell.get();
        if (b.in_preliminary_pass && b.remaining_preliminary == 0)
            || (b.has_global_limit && b.remaining_global == 0)
        {
            return false;
        }
        if b.in_preliminary_pass {
            b.remaining_preliminary -= 1;
        }
        if b.has_global_limit {
            b.remaining_global -= 1;
        }
        cell.set(b);
        true
    })
}

/// Sets the global budget of a call and restores the previous one on drop.
pub(super) struct ScopedIterationBudget {
    previous: IterationBudget,
}

impl ScopedIterationBudget {
    pub(super) fn new(max_recursive_iterations: u32) -> Self {
        let previous = BUDGET.with(Cell::get);
        BUDGET.with(|c| {
            c.set(IterationBudget {
                has_global_limit: max_recursive_iterations != 0,
                remaining_global: max_recursive_iterations,
                in_preliminary_pass: false,
                remaining_preliminary: 0,
            })
        });
        ScopedIterationBudget { previous }
    }
}

impl Drop for ScopedIterationBudget {
    fn drop(&mut self) {
        BUDGET.with(|c| c.set(self.previous));
    }
}

/// Limits the comparisons of the preliminary (constitutional) pass.
pub(super) struct ScopedPreliminaryBudget;

impl ScopedPreliminaryBudget {
    pub(super) fn new() -> Self {
        BUDGET.with(|c| {
            let mut b = c.get();
            b.in_preliminary_pass = true;
            b.remaining_preliminary = CONSTITUTIONAL_RULE_TIMEOUT;
            c.set(b);
        });
        ScopedPreliminaryBudget
    }
}

impl Drop for ScopedPreliminaryBudget {
    fn drop(&mut self) {
        BUDGET.with(|c| {
            let mut b = c.get();
            b.in_preliminary_pass = false;
            b.remaining_preliminary = 0;
            c.set(b);
        });
    }
}
