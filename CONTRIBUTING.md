# Contributing

## Commits

Message format:

```
<type>(<scope>): <imperative summary>

What changed, in two or three sentences. Name the public contract if one
changed.

Why: the problem this solves, or the invariant it protects.

Verified: the command that was run.
```

Types: `feat`, `fix`, `refactor`, `perf`, `test`, `docs`, `chore`, `ci`,
`release`.

Scopes name the layer or crate being changed: `core`, `math`, `model`,
`element`, `analysis`, `loads`, `tie`, `modeling`, `design`, `steel`, `io`,
`cli`, `desktop`, `scene`, `section`.

A commit should take about twenty minutes to review. Larger changes are split,
usually into a refactor and then the feature, or the feature and then its
tests.

Refactor commits keep behaviour identical. When a refactor would change
behaviour, make that a separate commit.

`target/` and other build output stay out of the repository. `Cargo.lock` is
committed because the build reads it.

## Review

- Every feature commit includes a test that fails without the change. If a test
  cannot express it, the commit body explains why.
- Public types and functions have doc comments. A new public contract gets its
  test in the same commit.
- Unit conversions happen at defined boundaries. A conversion factor buried in
  analysis code is hard to find and easy to get wrong.
- Nothing under the analysis modules imports a scene, a viewport or a
  command-line type.
- Tolerances in interaction code are named constants rather than literals
  scattered through event handlers.
- Solver failures return typed errors. A NaN reaching a result means one of
  those error paths was skipped.
- A new default is a field the reader can see, and its source is documented.

## Checks

Run before committing code:

```
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

These are the commands the CI workflow will run, with the same flags.
