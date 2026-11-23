//! Named numerical tolerances for the solver.
//!
//! Each tolerance states what it is measured against. An absolute tolerance on
//! a length, a pivot and a determinant would mix three different units, so the
//! scale-dependent checks below are relative to a quantity of the same kind.
//! The constants live here so the factorisation and the checks around it agree
//! on what "singular" and "converged" mean.

use super::Real;

/// Absolute separation below which two points are treated as coincident.
/// Measured in the coordinate units of the model.
pub const LENGTH_EPS: Real = 1.0e-12;

/// Relative determinant threshold for `Mat3::inverse`. Compared with the cube
/// of the matrix's largest entry, so a uniformly scaled matrix has the same
/// verdict as the unscaled one.
pub const DET_RELATIVE: Real = 1.0e-12;

/// Relative pivot threshold for the banded factorisation. Compared with the
/// original diagonal of the equation being eliminated, so stiffness matrices
/// in different units behave the same way.
pub const PIVOT_RELATIVE: Real = 1.0e-12;

/// Largest relative residual accepted after a solve: `‖Ax − b‖∞ / ‖b‖∞`.
pub const RESIDUAL: Real = 1.0e-8;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tolerances_are_positive_and_ordered() {
        // black_box keeps these runtime checks from being folded into
        // compile-time constants, which is what the test wants: it guards the
        // ordering, not the literals.
        let pivot = std::hint::black_box(PIVOT_RELATIVE);
        let residual = std::hint::black_box(RESIDUAL);
        let length = std::hint::black_box(LENGTH_EPS);
        let det = std::hint::black_box(DET_RELATIVE);
        assert!(pivot > 0.0 && residual > 0.0 && length > 0.0 && det > 0.0);
        // A pivot tolerance tighter than the residual tolerance would reject
        // systems the solver could otherwise solve accurately.
        assert!(pivot <= residual);
    }
}
