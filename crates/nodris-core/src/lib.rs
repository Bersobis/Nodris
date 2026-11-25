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
//! describe, with stable ids. [`element`] is the formulation contract built
//! on that model and maths, including the one six-DOF-per-node order.
//!
//! The solver's guarantees and limits are recorded in
//! `docs/SOLVER-GUARANTEES.md` at the workspace root.
//!
//! Nodris is a rebuild of Stodola, reconstructed from a plan that numbers its
//! commits, `P0-02` and so on. Those ids refer to the plan.

pub mod element;
pub mod math;
pub mod model;
pub mod units;

pub use element::{
    end_forces, local_axes, local_frame_stiffness, member_dofs, node_dofs, to_global, Element,
    ElementError, FrameElement, InitialStress, ShellElement, ShellForces, TrussElement, DOF_NAMES,
    DOF_PER_NODE,
};
pub use math::{
    renumber, BandedMatrix, Factored, Mat3, NewIndex, OldIndex, Permutation, SolveError, Vec3,
};
pub use model::{
    Area, AreaId, AreaKind, AreaLoad, BeamCoupling, BeamEnd, Combination, CombinationId,
    ElementKind, Load, LoadCase, LoadCaseId, Material, MaterialId, Member, MemberId, Model,
    ModelError, ModelIndex, NodalLoad, Node, NodeId, ObjectKind, ObjectRef, PatchCoupling,
    PatchDirection, PatchEdge, PatchId, PatchLoad, PatchTie, PatchTieId, Restraints, Section,
    SectionId, Severity, SplinePatch, Support, SupportDefinition, SupportDefinitionId,
    TieConstraint, ValidationIssue,
};
