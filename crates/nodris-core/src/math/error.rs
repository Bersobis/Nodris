//! Typed failures from the banded solver.
//!
//! A string error tells a reader what happened but not where it happened or
//! what to do next. [`SolveError`] carries the offending equation number and
//! the reason, so a caller can report the degree of freedom that failed and
//! decide whether the run can continue.

use std::fmt;

use super::tolerance::RESIDUAL;
use super::Real;

/// Why a banded system could not be solved.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SolveError {
    /// The factorisation reached a pivot at or below the relative threshold
    /// [`PIVOT_RELATIVE`](super::tolerance::PIVOT_RELATIVE). The system is
    /// singular or too ill-conditioned to solve.
    Singular { equation: usize },
    /// The factorisation reached a negative pivot: the matrix is symmetric but
    /// not positive definite, so the unpivoted `L·D·Lᵀ` factorisation does
    /// not apply.
    NotPositiveDefinite { equation: usize },
    /// A pivot was not a finite number, so the matrix already holds NaN or
    /// infinity.
    NotFinite { equation: usize },
    /// The right-hand side held a non-finite value. Loads enter the solver
    /// from outside, so this is checked separately from the matrix pivots.
    NonFiniteRightHandSide { equation: usize },
    /// The solved system left a relative residual above [`RESIDUAL`].
    ResidualTooLarge { equation: usize, residual: Real },
}

impl SolveError {
    /// The equation whose pivot, right-hand side or residual failed.
    pub fn equation(&self) -> usize {
        match self {
            SolveError::Singular { equation }
            | SolveError::NotPositiveDefinite { equation }
            | SolveError::NotFinite { equation }
            | SolveError::NonFiniteRightHandSide { equation }
            | SolveError::ResidualTooLarge { equation, .. } => *equation,
        }
    }
}

impl fmt::Display for SolveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolveError::Singular { equation } => {
                write!(f, "singular stiffness matrix at equation {equation}")
            }
            SolveError::NotPositiveDefinite { equation } => write!(
                f,
                "stiffness matrix is not positive definite at equation {equation}"
            ),
            SolveError::NotFinite { equation } => write!(
                f,
                "stiffness matrix has a non-finite pivot at equation {equation}"
            ),
            SolveError::NonFiniteRightHandSide { equation } => {
                write!(f, "right-hand side is not finite at equation {equation}")
            }
            SolveError::ResidualTooLarge { equation, residual } => write!(
                f,
                "relative residual {residual} exceeds {RESIDUAL} at equation {equation}"
            ),
        }
    }
}

impl std::error::Error for SolveError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_names_its_equation() {
        let errors = [
            SolveError::Singular { equation: 3 },
            SolveError::NotPositiveDefinite { equation: 4 },
            SolveError::NotFinite { equation: 5 },
            SolveError::NonFiniteRightHandSide { equation: 6 },
            SolveError::ResidualTooLarge {
                equation: 7,
                residual: 1.0,
            },
        ];
        for error in errors {
            assert!(error.to_string().contains(&error.equation().to_string()));
        }
    }
}
