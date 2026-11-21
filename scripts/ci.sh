#!/usr/bin/env bash
# The gate the CI workflow runs. Run it from anywhere; it moves to the
# repository root first, so a local run exercises exactly what the pipeline
# does.
#
#     bash scripts/ci.sh
#
# Steps: rustfmt, cargo check, clippy with warnings denied, tests.
set -euo pipefail

cd "$(dirname "$0")/.."

echo "== rustfmt"
cargo fmt --all -- --check

echo "== cargo check"
cargo check --workspace --all-targets

echo "== clippy"
cargo clippy --workspace --all-targets -- -D warnings

echo "== tests"
cargo test --workspace

echo "all gates passed"
