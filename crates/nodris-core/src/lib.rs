//! Structural analysis kernel for Stodola.
//!
//! This crate is the product's centre of gravity: the numerical kernel, the
//! model schema, the element formulations and the analysis pipeline all live
//! here. Analysis code never depends on UI or scene code; the CLI and the
//! desktop application are consumers of this crate.
//!
//! The crate starts empty on purpose. Modules arrive with the phase that owns
//! them, each with the test that would fail without it.
