//! Reads and writes the versioned model document.
//!
//! Strict import accepts only a declared v1 document with exactly the fields
//! the schema names. Permissive import accepts an unversioned document and
//! keeps fields this build does not know, recording each deviation as a note.
//! Type errors and unknown enum spellings fail in both modes because neither
//! can be repaired without guessing at the user's data.
//!
//! A v0 document stays v0 when written back. Changing the declared version is
//! a migration, not a side effect of reading or saving.
//!
//! Duplicate object keys are resolved by the JSON library, which keeps the
//! last value. Reporting them needs the canonical parser planned for a later
//! step.

use std::fmt;

use nodris_core::model::Model;
use serde_json::{Map, Value};

pub use crate::model_json::ROOT_KEYS;
use crate::model_json::{model_from_value, model_to_value};
#[cfg(test)]
use crate::model_json::{
    AREA_KEYS, AREA_LOAD_KEYS, COMBINATION_KEYS, MATERIAL_KEYS, MEMBER_KEYS, MEMBER_LOAD_KEYS,
    NODAL_LOAD_KEYS, NODE_KEYS, PATCH_KEYS, PATCH_LOAD_KEYS, PATCH_TIE_KEYS, RESTRAINT_KEYS,
    SECTION_KEYS, SUPPORT_DEFINITION_KEYS, SUPPORT_KEYS,
};

/// Version this build writes.
pub const SCHEMA_VERSION: u32 = 1;

/// Name of the version field in a model document.
pub const SCHEMA_KEY: &str = "schemaVersion";

/// Published machine-readable contract for a v1 document.
pub const MODEL_SCHEMA_V1: &str = include_str!("../schema/model-v1.schema.json");

/// Root fields a strict v1 document must carry.
pub const ROOT_REQUIRED: &[&str] = ROOT_KEYS;

/// Root fields that may be omitted. The v1 document has none: permissive
/// import treats a missing collection as empty, but strict import requires it.
pub const ROOT_OPTIONAL: &[&str] = &[];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ImportMode {
    Strict,
    Permissive,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ImportNoteKind {
    /// The document had no version field and was read as v0.
    LegacyDocument,
    /// A field absent from the document was read as an empty collection.
    MissingCollection,
    /// A field this build does not know was found and kept.
    UnknownField,
    /// Migration changed the declared version.
    Migrated,
    /// Migration was asked for on a document that is already current.
    AlreadyCurrent,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ImportNote {
    pub kind: ImportNoteKind,
    /// Path of the affected field, for example `$.nodes[2].extra`.
    pub path: String,
    pub message: String,
}

impl ImportNote {
    pub(crate) fn new(
        kind: ImportNoteKind,
        path: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            kind,
            path: path.into(),
            message: message.into(),
        }
    }
}

#[derive(Debug)]
pub enum InterchangeError {
    InvalidJson(serde_json::Error),
    ExpectedType {
        path: String,
        expected: &'static str,
    },
    MissingField {
        path: String,
    },
    UnknownField {
        path: String,
    },
    InvalidNumber {
        path: String,
    },
    InvalidEnum {
        path: String,
        value: String,
    },
    UnsupportedVersion(u32),
    InvalidEntity {
        path: String,
        message: String,
    },
    ModelValidation(String),
}

impl fmt::Display for InterchangeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidJson(error) => write!(f, "invalid JSON: {error}"),
            Self::ExpectedType { path, expected } => {
                write!(f, "{path} must be {expected}")
            }
            Self::MissingField { path } => write!(f, "{path} is required"),
            Self::UnknownField { path } => write!(f, "{path} is not in the v1 schema"),
            Self::InvalidNumber { path } => write!(f, "{path} must be a finite number"),
            Self::InvalidEnum { path, value } => {
                write!(f, "{path} has unknown spelling {value:?}")
            }
            Self::UnsupportedVersion(version) => {
                write!(
                    f,
                    "schema version {version} is newer than {}",
                    SCHEMA_VERSION
                )
            }
            Self::InvalidEntity { path, message } => write!(f, "{path}: {message}"),
            Self::ModelValidation(message) => write!(f, "invalid model: {message}"),
        }
    }
}

impl std::error::Error for InterchangeError {}

impl From<serde_json::Error> for InterchangeError {
    fn from(error: serde_json::Error) -> Self {
        Self::InvalidJson(error)
    }
}

/// A model with the fields of its source document that are not model data.
#[derive(Clone, Debug)]
pub struct ModelDocument {
    pub model: Model,
    /// Declared version. Zero means the source had no version field.
    pub schema_version: u32,
    /// Root fields unknown to this schema, kept for a later write.
    pub unknown_root: Map<String, Value>,
}

#[derive(Clone, Debug)]
pub struct ModelRead {
    pub document: ModelDocument,
    pub notes: Vec<ImportNote>,
}

/// The result of a migration. It is new text; the source is never edited.
#[derive(Clone, Debug)]
pub struct Migration {
    pub from_version: u32,
    pub to_version: u32,
    pub text: String,
    pub notes: Vec<ImportNote>,
}

/// Rejects fields the schema does not name, or records them in permissive mode.
pub(crate) fn check_keys(
    object: &Map<String, Value>,
    path: &str,
    allowed: &[&str],
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<(), InterchangeError> {
    for key in object.keys() {
        if allowed.contains(&key.as_str()) {
            continue;
        }
        let field_path = if path == "$" {
            format!("$.{key}")
        } else {
            format!("{path}.{key}")
        };
        match mode {
            ImportMode::Strict => {
                return Err(InterchangeError::UnknownField { path: field_path });
            }
            ImportMode::Permissive => notes.push(ImportNote::new(
                ImportNoteKind::UnknownField,
                field_path.clone(),
                format!("{field_path} is not in the v1 schema and was kept"),
            )),
        }
    }
    Ok(())
}

fn read_version(
    root: &Map<String, Value>,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<u32, InterchangeError> {
    let Some(value) = root.get(SCHEMA_KEY) else {
        return match mode {
            ImportMode::Strict => Err(InterchangeError::MissingField {
                path: format!("$.{SCHEMA_KEY}"),
            }),
            ImportMode::Permissive => {
                notes.push(ImportNote::new(
                    ImportNoteKind::LegacyDocument,
                    format!("$.{SCHEMA_KEY}"),
                    "no schema version, read as legacy v0",
                ));
                Ok(0)
            }
        };
    };
    let version = value
        .as_u64()
        .and_then(|number| u32::try_from(number).ok())
        .ok_or_else(|| InterchangeError::ExpectedType {
            path: format!("$.{SCHEMA_KEY}"),
            expected: "non-negative integer",
        })?;
    if version > SCHEMA_VERSION {
        return Err(InterchangeError::UnsupportedVersion(version));
    }
    if version != SCHEMA_VERSION && mode == ImportMode::Strict {
        return Err(InterchangeError::ExpectedType {
            path: format!("$.{SCHEMA_KEY}"),
            expected: "1",
        });
    }
    Ok(version)
}

/// Reads text as a v1 document, or as v0 when permissive.
pub fn read_model_document(text: &str, mode: ImportMode) -> Result<ModelRead, InterchangeError> {
    let value: Value = serde_json::from_str(text)?;
    let root = value
        .as_object()
        .ok_or_else(|| InterchangeError::ExpectedType {
            path: "$".to_string(),
            expected: "object",
        })?;
    let mut notes = Vec::new();
    check_keys(root, "$", ROOT_KEYS, mode, &mut notes)?;
    let schema_version = read_version(root, mode, &mut notes)?;
    let model = model_from_value(root, mode, &mut notes)?;
    let unknown_root = root
        .iter()
        .filter(|(key, _)| !ROOT_KEYS.contains(&key.as_str()))
        .map(|(key, value)| (key.clone(), value.clone()))
        .collect();
    Ok(ModelRead {
        document: ModelDocument {
            model,
            schema_version,
            unknown_root,
        },
        notes,
    })
}

/// Builds the JSON tree for a document, including any preserved root fields.
pub fn document_to_value(document: &ModelDocument) -> Result<Value, InterchangeError> {
    let mut root = match model_to_value(&document.model)? {
        Value::Object(root) => root,
        Value::Array(_) | Value::Null | Value::Bool(_) | Value::Number(_) | Value::String(_) => {
            return Err(InterchangeError::ExpectedType {
                path: "$".to_string(),
                expected: "object",
            })
        }
    };
    root.insert(SCHEMA_KEY.to_string(), Value::from(document.schema_version));
    for (key, value) in &document.unknown_root {
        root.insert(key.clone(), value.clone());
    }
    Ok(Value::Object(root))
}

/// Serialises a JSON tree with sorted object keys and a trailing newline.
fn canonical_text(value: &Value) -> Result<String, InterchangeError> {
    let mut text = serde_json::to_string_pretty(value)?;
    text.push('\n');
    Ok(text)
}

/// Writes a document with sorted object keys and a trailing newline.
pub fn write_model_document(document: &ModelDocument) -> Result<String, InterchangeError> {
    canonical_text(&document_to_value(document)?)
}

/// Rewrites a legacy document as v1 without discarding fields this build does
/// not know. The parsed tree is edited directly, so unknown fields survive at
/// every depth. The input text is never modified.
pub fn migrate_model_text(raw: &str) -> Result<Migration, InterchangeError> {
    let read = read_model_document(raw, ImportMode::Permissive)?;
    let from_version = read.document.schema_version;
    let mut tree: Value = serde_json::from_str(raw)?;
    let root = tree
        .as_object_mut()
        .ok_or_else(|| InterchangeError::ExpectedType {
            path: "$".to_string(),
            expected: "object",
        })?;
    root.insert(SCHEMA_KEY.to_string(), Value::from(SCHEMA_VERSION));
    let mut notes = read.notes;
    notes.push(ImportNote::new(
        if from_version == SCHEMA_VERSION {
            ImportNoteKind::AlreadyCurrent
        } else {
            ImportNoteKind::Migrated
        },
        format!("$.{SCHEMA_KEY}"),
        if from_version == SCHEMA_VERSION {
            "document already declares v1; canonical formatting applied".to_string()
        } else {
            format!("declared schema version changed from v{from_version} to v{SCHEMA_VERSION}")
        },
    ));
    Ok(Migration {
        from_version,
        to_version: SCHEMA_VERSION,
        text: canonical_text(&Value::Object(root.clone()))?,
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nodris_core::model::{
        Area, AreaId, AreaKind, Combination, CombinationId, ElementKind, Load, LoadCase,
        LoadCaseId, Material, MaterialId, Member, MemberId, Model, NodalLoad, Node, NodeId,
        PatchEdge, PatchId, PatchTie, PatchTieId, Restraints, Section, SectionId, SplinePatch,
        Support, SupportDefinition, SupportDefinitionId,
    };

    fn sample_document() -> ModelDocument {
        let mut model = Model::new();
        model
            .add_node(Node::new(NodeId::new(1), 1.0, 2.0, 3.0))
            .unwrap();
        model
            .add_node(Node::new(NodeId::new(2), 4.0, 5.0, 6.0))
            .unwrap();
        model
            .add_node(Node::new(NodeId::new(3), 7.0, 8.0, 9.0))
            .unwrap();
        model
            .add_node(Node::new(NodeId::new(4), 10.0, 11.0, 12.0))
            .unwrap();
        model
            .add_node(Node::new(NodeId::new(12), 13.0, 14.0, 15.0))
            .unwrap();
        model
            .add_node(Node::new(NodeId::new(13), 16.0, 17.0, 18.0))
            .unwrap();
        model
            .add_node(Node::new(NodeId::new(14), 19.0, 20.0, 21.0))
            .unwrap();
        model
            .add_node(Node::new(NodeId::new(15), 22.0, 23.0, 24.0))
            .unwrap();
        model
            .add_material(Material::new(
                MaterialId::new(3),
                "S355",
                7850.0,
                210_000_000.0,
                0.3,
            ))
            .unwrap();
        model
            .add_section(Section::new(
                SectionId::new(4),
                "IPE",
                0.01,
                1.0e-5,
                2.0e-5,
                3.0e-5,
                4.0e-7,
            ))
            .unwrap();
        model
            .add_member(Member::new(
                MemberId::new(5),
                NodeId::new(1),
                NodeId::new(2),
                SectionId::new(4),
                0.25,
                ElementKind::Frame,
            ))
            .unwrap();
        model
            .add_area(Area::new(
                AreaId::new(6),
                [
                    NodeId::new(1),
                    NodeId::new(2),
                    NodeId::new(3),
                    NodeId::new(4),
                ],
                AreaKind::Shell,
            ))
            .unwrap();
        model
            .add_patch(SplinePatch::unweighted(
                PatchId::new(7),
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
                vec![
                    vec![NodeId::new(1), NodeId::new(2)],
                    vec![NodeId::new(3), NodeId::new(4)],
                ],
            ))
            .unwrap();
        model
            .add_patch(SplinePatch::unweighted(
                PatchId::new(16),
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
                vec![
                    vec![NodeId::new(12), NodeId::new(13)],
                    vec![NodeId::new(14), NodeId::new(15)],
                ],
            ))
            .unwrap();
        model
            .add_load_case(LoadCase::new(LoadCaseId::new(8), "LC").with_self_weight())
            .unwrap();
        model
            .add_load(Load::Nodal {
                case: LoadCaseId::new(8),
                load: NodalLoad::new(NodeId::new(1), [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]),
            })
            .unwrap();
        model
            .add_combination(Combination::new(
                CombinationId::new(9),
                "ULS",
                vec![(LoadCaseId::new(8), 1.35)],
            ))
            .unwrap();
        model
            .add_support_definition(SupportDefinition::new(
                SupportDefinitionId::new(10),
                "PIN",
                Restraints::new(true, true, true, false, false, false),
            ))
            .unwrap();
        model
            .add_support(Support::new(NodeId::new(1), SupportDefinitionId::new(10)))
            .unwrap();
        model
            .add_patch_tie(PatchTie::new(
                PatchTieId::new(11),
                PatchId::new(7),
                PatchEdge::Xi0,
                PatchId::new(16),
                PatchEdge::Xi0,
            ))
            .unwrap();
        model.apply_supports().unwrap();
        ModelDocument {
            model,
            schema_version: SCHEMA_VERSION,
            unknown_root: Map::new(),
        }
    }

    #[test]
    fn a_document_round_trips_in_strict_mode() {
        let document = sample_document();
        let text = write_model_document(&document).unwrap();
        let read = read_model_document(&text, ImportMode::Strict).unwrap();
        assert!(read.notes.is_empty());
        assert_eq!(read.document.schema_version, SCHEMA_VERSION);
        assert_eq!(
            write_model_document(&read.document).unwrap(),
            text,
            "a strict round trip must be byte-stable"
        );
    }

    #[test]
    fn every_entity_survives_the_round_trip() {
        let document = sample_document();
        let text = write_model_document(&document).unwrap();
        let read = read_model_document(&text, ImportMode::Strict).unwrap();
        let original = document_to_value(&document).unwrap();
        let restored = document_to_value(&read.document).unwrap();
        assert_eq!(restored, original);
    }

    #[test]
    fn keys_are_sorted_in_the_written_document() {
        let text = write_model_document(&sample_document()).unwrap();
        assert!(text.lines().nth(1).unwrap().starts_with("  \"areas\""));
        assert!(text.ends_with('\n'));
    }

    #[test]
    fn strict_mode_rejects_an_unknown_root_field() {
        let text = write_model_document(&sample_document()).unwrap();
        let patched = text.replace("  \"areas\": [", "  \"futureField\": 7,\n  \"areas\": [");
        assert!(patched != text);
        let error = read_model_document(&patched, ImportMode::Strict).unwrap_err();
        assert!(matches!(error, InterchangeError::UnknownField { .. }));
        assert!(error.to_string().contains("$.futureField"));
    }

    #[test]
    fn permissive_mode_keeps_an_unknown_root_field() {
        let text = write_model_document(&sample_document()).unwrap();
        let patched = text.replace(
            "  \"areas\": [",
            "  \"futureField\": {\"deep\": [1, 2]},\n  \"areas\": [",
        );
        let read = read_model_document(&patched, ImportMode::Permissive).unwrap();
        assert!(read.notes.iter().any(|note| {
            note.kind == ImportNoteKind::UnknownField && note.path == "$.futureField"
        }));
        assert_eq!(
            read.document.unknown_root.get("futureField"),
            Some(&serde_json::json!({"deep": [1, 2]}))
        );
        let written = write_model_document(&read.document).unwrap();
        assert!(written.contains("\"futureField\""));
        assert_eq!(
            read_model_document(&written, ImportMode::Permissive)
                .unwrap()
                .document
                .unknown_root,
            read.document.unknown_root
        );
    }

    #[test]
    fn a_nested_unknown_field_is_noted_but_not_kept() {
        let text = write_model_document(&sample_document()).unwrap();
        let patched = text.replace("\"x\": 1.0,", "\"x\": 1.0,\n      \"localNote\": 9,");
        assert!(patched != text);
        let strict = read_model_document(&patched, ImportMode::Strict).unwrap_err();
        assert!(matches!(strict, InterchangeError::UnknownField { .. }));
        let permissive = read_model_document(&patched, ImportMode::Permissive).unwrap();
        assert!(permissive.notes.iter().any(|note| {
            note.kind == ImportNoteKind::UnknownField && note.path.contains("localNote")
        }));
        assert!(!write_model_document(&permissive.document)
            .unwrap()
            .contains("localNote"));
    }

    #[test]
    fn an_explicit_v0_document_is_rejected_only_in_strict_mode() {
        let mut document = sample_document();
        document.schema_version = 0;
        let text = write_model_document(&document).unwrap();
        let strict = read_model_document(&text, ImportMode::Strict).unwrap_err();
        assert!(matches!(strict, InterchangeError::ExpectedType { .. }));
        let permissive = read_model_document(&text, ImportMode::Permissive).unwrap();
        assert_eq!(permissive.document.schema_version, 0);
        assert!(permissive.notes.is_empty());
    }

    #[test]
    fn a_document_without_a_version_field_is_read_as_v0() {
        let text = write_model_document(&sample_document()).unwrap();
        let unversioned: String = text
            .lines()
            .filter(|line| !line.contains("\"schemaVersion\""))
            .collect::<Vec<_>>()
            .join("\n");
        let strict = read_model_document(&unversioned, ImportMode::Strict).unwrap_err();
        assert!(matches!(strict, InterchangeError::MissingField { .. }));
        let read = read_model_document(&unversioned, ImportMode::Permissive).unwrap();
        assert_eq!(read.document.schema_version, 0);
        assert!(read
            .notes
            .iter()
            .any(|note| note.kind == ImportNoteKind::LegacyDocument));
    }

    #[test]
    fn a_newer_version_is_rejected_in_both_modes() {
        let mut document = sample_document();
        document.schema_version = 2;
        let text = write_model_document(&document).unwrap();
        for mode in [ImportMode::Strict, ImportMode::Permissive] {
            let error = read_model_document(&text, mode).unwrap_err();
            assert!(matches!(error, InterchangeError::UnsupportedVersion(2)));
        }
    }

    #[test]
    fn a_missing_collection_is_empty_in_permissive_mode() {
        let text = write_model_document(&sample_document()).unwrap();
        let mut tree: Value = serde_json::from_str(&text).unwrap();
        tree.as_object_mut().unwrap().remove("areas");
        let trimmed = serde_json::to_string(&tree).unwrap();
        let strict = read_model_document(&trimmed, ImportMode::Strict).unwrap_err();
        assert!(matches!(strict, InterchangeError::MissingField { .. }));
        let read = read_model_document(&trimmed, ImportMode::Permissive).unwrap();
        assert!(read.document.model.areas().is_empty());
        assert!(read.notes.iter().any(|note| {
            note.kind == ImportNoteKind::MissingCollection && note.path == "$.areas"
        }));
    }

    #[test]
    fn a_wrong_type_is_fatal_in_both_modes() {
        let text = write_model_document(&sample_document()).unwrap();
        let patched = text.replace("\"x\": 1.0,", "\"x\": \"one\",");
        for mode in [ImportMode::Strict, ImportMode::Permissive] {
            let error = read_model_document(&patched, mode).unwrap_err();
            assert!(matches!(error, InterchangeError::ExpectedType { .. }));
        }
    }

    #[test]
    fn an_unknown_enum_is_fatal_in_both_modes() {
        let text = write_model_document(&sample_document()).unwrap();
        let patched = text.replace("\"kind\": \"frame\"", "\"kind\": \"beam\"");
        for mode in [ImportMode::Strict, ImportMode::Permissive] {
            let error = read_model_document(&patched, mode).unwrap_err();
            assert!(matches!(error, InterchangeError::InvalidEnum { .. }));
        }
    }

    #[test]
    fn malformed_json_is_reported_with_its_position() {
        let error = read_model_document("{ not json }", ImportMode::Strict).unwrap_err();
        assert!(matches!(error, InterchangeError::InvalidJson(_)));
        assert!(error.to_string().contains("line"));
    }

    #[test]
    fn the_bundled_schema_is_a_valid_v1_document() {
        let schema: Value = serde_json::from_str(MODEL_SCHEMA_V1).unwrap();
        assert_eq!(schema.get("type"), Some(&Value::from("object")));
        let version = schema
            .get("properties")
            .and_then(|properties| properties.get(SCHEMA_KEY))
            .and_then(|field| field.get("const"))
            .and_then(Value::as_u64);
        assert_eq!(version, Some(u64::from(SCHEMA_VERSION)));
        let required = schema
            .get("required")
            .and_then(Value::as_array)
            .unwrap()
            .iter()
            .map(|value| value.as_str().unwrap())
            .collect::<Vec<_>>();
        assert_eq!(required, ROOT_REQUIRED);
    }

    #[test]
    fn schema_fields_match_the_reader_fields() {
        let schema: Value = serde_json::from_str(MODEL_SCHEMA_V1).unwrap();
        let definitions = schema.get("$defs").unwrap().as_object().unwrap();
        let pairs = [
            ("node", NODE_KEYS),
            ("material", MATERIAL_KEYS),
            ("section", SECTION_KEYS),
            ("member", MEMBER_KEYS),
            ("area", AREA_KEYS),
            ("patch", PATCH_KEYS),
            ("loadCase", crate::model_json::LOAD_CASE_KEYS),
            ("nodalLoad", NODAL_LOAD_KEYS),
            ("memberLoad", MEMBER_LOAD_KEYS),
            ("areaLoad", AREA_LOAD_KEYS),
            ("patchLoad", PATCH_LOAD_KEYS),
            ("combination", COMBINATION_KEYS),
            ("supportDefinition", SUPPORT_DEFINITION_KEYS),
            ("restraints", RESTRAINT_KEYS),
            ("support", SUPPORT_KEYS),
            ("patchTie", PATCH_TIE_KEYS),
        ];
        for (name, keys) in pairs {
            let declared = definitions
                .get(name)
                .unwrap()
                .get("properties")
                .unwrap()
                .as_object()
                .unwrap()
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>();
            let mut expected = keys.to_vec();
            expected.sort_unstable();
            assert_eq!(
                declared, expected,
                "field drift in schema definition {name}"
            );
        }
    }

    #[test]
    fn an_empty_model_round_trips() {
        let document = ModelDocument {
            model: Model::new(),
            schema_version: SCHEMA_VERSION,
            unknown_root: Map::new(),
        };
        let text = write_model_document(&document).unwrap();
        let read = read_model_document(&text, ImportMode::Strict).unwrap();
        assert_eq!(write_model_document(&read.document).unwrap(), text);
    }
}
