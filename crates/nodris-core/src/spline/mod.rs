//! B-spline spaces and tensor-product patch geometry.
//!
//! The isogeometric element builds on these types: a [`BSplineSpace`] owns one
//! direction's knot vector and basis, and a [`BSplinePatch`] combines two of
//! them with a control net. Knot vectors are open (clamped), indices are
//! zero-based, and the parameter domain is the closed interval
//! `[knots[degree], knots[n]]`.

mod patch;
mod space;

pub use patch::{BSplinePatch, PatchBasis, SplineError};
pub use space::BSplineSpace;
