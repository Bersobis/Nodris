//! Fixed-size geometry types and, from Phase 1 onwards, the linear solver.
//!
//! Everything here is `Copy` and allocation-free. Element assembly runs these
//! operations in inner loops, and the kernel is a library with no dependencies,
//! so the types are written out rather than pulled from a linear-algebra crate.

pub mod mat3;
pub mod vec3;

pub use mat3::Mat3;
pub use vec3::Vec3;

/// Analysis real. One name for the type, so switching precision later is a
/// single change.
pub type Real = f64;

/// Smallest length treated as non-zero when normalising a direction.
pub const EPS: Real = 1.0e-12;
