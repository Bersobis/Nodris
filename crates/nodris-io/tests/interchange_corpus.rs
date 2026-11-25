use nodris_io::{
    migrate_model_text, read_model_document, write_model_document, ImportMode, ImportNoteKind,
    InterchangeError, SCHEMA_VERSION,
};
use serde_json::Value;

const V1: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/model-v1-minimal.json"
));
const V0: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/model-v0-legacy.json"
));

#[test]
fn a_v1_fixture_round_trips_byte_canonically() {
    let read = read_model_document(V1, ImportMode::Strict).unwrap();
    assert!(read.notes.is_empty());
    assert_eq!(write_model_document(&read.document).unwrap(), V1);
}

#[test]
fn a_v1_fixture_is_already_current_for_migration() {
    let migration = migrate_model_text(V1).unwrap();
    assert_eq!(migration.from_version, SCHEMA_VERSION);
    assert_eq!(migration.to_version, SCHEMA_VERSION);
    assert_eq!(migration.text, V1);
    assert!(migration
        .notes
        .iter()
        .any(|note| note.kind == ImportNoteKind::AlreadyCurrent));
}

#[test]
fn a_legacy_fixture_changes_only_the_declared_version() {
    let source: Value = serde_json::from_str(V0).unwrap();
    let migration = migrate_model_text(V0).unwrap();
    assert_eq!(migration.from_version, 0);
    assert_eq!(migration.to_version, SCHEMA_VERSION);
    assert!(migration
        .notes
        .iter()
        .any(|note| note.kind == ImportNoteKind::Migrated));

    let mut expected = source;
    expected
        .as_object_mut()
        .unwrap()
        .insert("schemaVersion".into(), Value::from(SCHEMA_VERSION));
    let migrated: Value = serde_json::from_str(&migration.text).unwrap();
    assert_eq!(migrated, expected);
}

#[test]
fn migration_keeps_unknown_fields_at_every_depth() {
    let migration = migrate_model_text(V0).unwrap();
    let tree: Value = serde_json::from_str(&migration.text).unwrap();
    assert_eq!(tree["futureRoot"]["keep"][1]["deep"], Value::Bool(true));
    assert_eq!(tree["nodes"][0]["legacyTag"], Value::from("kept"));
    assert_eq!(tree["materials"][0]["extra"]["cadLayer"], Value::from("A"));
    assert_eq!(
        tree["combinations"][0]["factors"][0]["note"],
        Value::from("legacy")
    );
    assert_eq!(tree["patches"][0]["sourceTable"], Value::from("legacy"));
}

#[test]
fn a_migrated_document_is_canonical_and_readable() {
    let migration = migrate_model_text(V0).unwrap();
    assert_eq!(
        migration.text,
        serde_json::to_string_pretty(&serde_json::from_str::<Value>(&migration.text).unwrap())
            .unwrap()
            + "\n"
    );
    let read = read_model_document(&migration.text, ImportMode::Permissive).unwrap();
    assert_eq!(read.document.schema_version, SCHEMA_VERSION);
    assert_eq!(
        read.document.unknown_root["futureRoot"]["keep"][1]["deep"],
        Value::Bool(true)
    );
    assert_eq!(read.document.model.nodes().len(), 8);
}

#[test]
fn migration_reports_a_bad_legacy_document() {
    let invalid = V0.replace("\"rotation\": 0.25", "\"rotation\": \"crooked\"");
    assert!(invalid != V0);
    let error = migrate_model_text(&invalid).unwrap_err();
    assert!(matches!(error, InterchangeError::ExpectedType { .. }));
    assert!(error.to_string().contains("rotation"));
}

#[test]
fn a_broken_reference_names_the_object_path() {
    let invalid = V1.replace("\"section\": 4", "\"section\": 999");
    assert!(invalid != V1);
    let error = read_model_document(&invalid, ImportMode::Permissive).unwrap_err();
    assert!(matches!(error, InterchangeError::InvalidEntity { .. }));
    let message = error.to_string();
    assert!(message.contains("$.members[]"), "{message}");
    assert!(message.contains("999"), "{message}");
}
