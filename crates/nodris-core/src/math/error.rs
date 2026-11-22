//! Typed failures from the banded solver.
//!
//! A string error tells a reader what happened but not where it happened or
//! what to do next. [`SolveError`] carries the offending equation number and
//! the reason, so a caller can report the degree of freedom that failed and
//! decide whether the run can continue.

use std::fmt;

/// Why a banded system could not be solved.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SolveError {
    /// The factorisation reached a pivot at or below [`super::tolerance::PIVOT`].
    /// The system is singular or too ill-conditioned to solve.
    Singular { equation: usize },
    /// The factorisation reached a negative pivot: the matrix is symmetric but
    /// not positive definite, so the unpivoted `L·D·Lᵀ` factorisation does
    /// not apply.
    NotPositiveDefinite { equation: usize },
    /// A pivot was not a finite number, so the matrix already holds NaN or
    /// infinity.
    NotFinite { equation: usize },
}

impl SolveError {
    /// The equation whose pivot failed.
    pub fn equation(&self) -> usize {
        match self {
            SolveError::Singular { equation }
            | SolveError::NotPositiveDefinite { equation }
            | SolveError::NotFinite { equation } => *equation,
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
        ];
        for error in errors {
            assert!(error.to_string().contains(&error.equation().to_string()));
        }
    }
}
