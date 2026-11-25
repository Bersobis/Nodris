//! Persistence and interchange for Nodris model documents.
//!
//! The on-disk format is JSON so a model can be diffed in git and read by any
//! tool. A document declares `schemaVersion`; the published contract is
//! `schema/model-v1.schema.json`. `interchange` is the one place that reads
//! and writes versioned documents, with strict and permissive import modes.
//!
//! Strict import is for fixtures and review: the document must declare v1 and
//! use only fields the schema names. Permissive import keeps root fields this
//! build does not know and records every deviation as an import note. Nested
//! unknown fields are reported but not yet kept, because that needs the
//! migration path planned for a later step.

mod model_json;

pub mod interchange;

pub use interchange::{
    document_to_value, read_model_document, write_model_document, ImportMode, ImportNote,
    ImportNoteKind, InterchangeError, ModelDocument, ModelRead, MODEL_SCHEMA_V1, ROOT_KEYS,
    ROOT_OPTIONAL, ROOT_REQUIRED, SCHEMA_KEY, SCHEMA_VERSION,
};
