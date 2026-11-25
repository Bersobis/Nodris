# Nodris

Nodris is a structural analysis program for frames, trusses, shells and
isogeometric patches. It is a rebuild of my first vibe-coded attempt called
Stodola. The work is split into phases, and ids such as `P0-02` label them.

Models are stored as JSON. The solver is a library crate; a command line and a
desktop window will call it.

## Status

The kernel and the model document reader exist so far. The kernel has no
elements and no user interface, so the project compiles but cannot analyse a
model yet.

## Layout

| Path | Contents |
|---|---|
| `crates/nodris-core` | numerical kernel, model schema, elements, analysis |
| `crates/nodris-io` | versioned model documents, strict and permissive import |
| [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | layering, crate boundaries and the reasons for them |
| [CONTRIBUTING.md](CONTRIBUTING.md) | commit and review conventions |

## Building

```
cargo check --workspace
cargo test --workspace
```

The compiler version is pinned in `rust-toolchain.toml`. The kernel uses no
external crates; the IO crate uses `serde_json` to parse documents.

## Rules

The command line and the viewport depend on `nodris-core`. Nothing in
`nodris-core` refers back to either of them.

Inside the kernel, coordinates and forces use one unit system. Section
properties are stored in the units of the table they came from and converted
when they enter the analysis.

[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) has the layer diagram and the crate
plan.

## Licence

MIT. See [LICENSE](LICENSE).
