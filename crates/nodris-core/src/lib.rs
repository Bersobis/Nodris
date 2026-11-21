//! Structural analysis kernel for Nodris.
//!
//! The crate holds the numerical kernel, the model schema, the element
//! formulations and the analysis pipeline. It depends on no other workspace
//! crate. The command line, the scene and the desktop application are written
//! against it.
//!
//! The first module is [`math`]: fixed-size vectors and matrices with the
//! operations element assembly needs.
//!
//! Nodris is a rebuild of Stodola, reconstructed from a plan that numbers its
//! commits, `P0-02` and so on. Those ids refer to the plan.

pub mod math;

pub use math::{Mat3, Vec3};
