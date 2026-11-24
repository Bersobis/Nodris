//! Structural analysis kernel for Nodris.
//!
//! The crate holds the numerical kernel, the model schema, the element
//! formulations and the analysis pipeline. It depends on no other workspace
//! crate. The command line, the scene and the desktop application are written
//! against it.
//!
//! The first module is [`math`]: vectors, matrices and banded storage.
//! [`units`] is the one place where section tables and material strengths
//! cross into analysis units. [`model`] holds the entities those values
//! describe, with stable ids.
//!
//! The solver's guarantees and limits are recorded in
//! `docs/SOLVER-GUARANTEES.md` at the workspace root.
//!
//! Nodris is a rebuild of Stodola, reconstructed from a plan that numbers its
//! commits, `P0-02` and so on. Those ids refer to the plan.

pub mod math;
pub mod model;
pub mod units;

pub use math::{
    renumber, BandedMatrix, Factored, Mat3, NewIndex, OldIndex, Permutation, SolveError, Vec3,
};
pub use model::{
    ElementKind, Material, MaterialId, Member, MemberId, Model, ModelError, ModelIndex, Node,
    NodeId, Section, SectionId,
};
