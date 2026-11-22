//! Fixed-size geometry types and the banded storage used by the solver.
//!
//! [`Mat3`] and [`Vec3`] are `Copy` and allocation-free. [`BandedMatrix`] owns
//! one `Vec` of coefficients and is the only type here that allocates.
//!
//! Element assembly runs these operations in inner loops, and the kernel is a
//! library with no dependencies, so the types are written out rather than
//! pulled from a linear-algebra crate.

pub mod error;
pub mod linalg;
mod mat3;
pub mod tolerance;
mod vec3;

pub use error::SolveError;
pub use linalg::BandedMatrix;
pub use mat3::Mat3;
pub use vec3::Vec3;

/// Analysis real. One name for the type, so switching precision later is a
/// single change.
pub type Real = f64;

/// Smallest magnitude the kernel treats as non-zero: a length when normalising
/// a direction, a determinant when inverting, and a point separation when
/// checking coincidence.
pub const EPS: Real = 1.0e-12;
