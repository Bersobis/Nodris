//! Numerical tolerances for the solver.
//!
//! The constants live here so the factorisation and the checks around it agree
//! on what "singular" and "converged" mean, instead of each call site typing
//! its own literal.

use super::Real;

/// Smallest pivot the banded factorisation accepts. A pivot at or below this
/// magnitude makes the system singular.
pub const PIVOT: Real = 1.0e-12;

/// Largest relative residual a solve may leave before it is called suspect.
pub const RESIDUAL: Real = 1.0e-8;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tolerances_are_ordered_as_the_solver_needs() {
        // black_box keeps these runtime checks from being folded into
        // compile-time constants, which is what the test wants: it guards the
        // ordering, not the literals.
        let pivot = std::hint::black_box(PIVOT);
        let residual = std::hint::black_box(RESIDUAL);
        assert!(pivot > 0.0);
        assert!(residual > 0.0);
        // A pivot tolerance tighter than the residual tolerance would reject
        // systems the solver could otherwise solve accurately.
        assert!(pivot <= residual);
    }
}
