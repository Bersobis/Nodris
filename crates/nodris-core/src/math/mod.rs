//! Fixed-size geometry types and the banded storage used by the solver.
//!
//! [`Mat3`] and [`Vec3`] are `Copy` and allocation-free. [`BandedMatrix`] owns
//! one `Vec` of coefficients and is the only type here that allocates.
//!
//! Element assembly runs these operations in inner loops, and the kernel is a
//! library with no dependencies, so the types are written out rather than
//! pulled from a linear-algebra crate.

pub mod error;
pub mod index;
pub mod linalg;
mod mat3;
pub mod tolerance;
mod vec3;

pub use error::SolveError;
pub use index::{NewIndex, OldIndex, Permutation};
pub use linalg::{renumber, BandedMatrix, Factored};
pub use mat3::Mat3;
pub use vec3::Vec3;

/// Analysis real. One name for the type, but not a free switch: the tolerances
/// in [`tolerance`] and the seeded test generator in [`linalg`] are tuned for
/// `f64` and have to be reviewed if this changes.
pub type Real = f64;
