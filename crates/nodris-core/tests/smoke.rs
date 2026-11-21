//! Smoke test for the kernel crate.
//!
//! It runs only compiled code: if the crate fails to build or its manifest
//! loses the version, this is the first test to notice.

#[test]
fn version_has_three_numeric_parts() {
    let version = env!("CARGO_PKG_VERSION");
    let parts: Vec<&str> = version.split('.').collect();
    assert_eq!(
        parts.len(),
        3,
        "version {version} should read major.minor.patch"
    );
    for part in parts {
        assert!(
            part.parse::<u32>().is_ok(),
            "version {version} has a non-numeric part: {part}"
        );
    }
}
