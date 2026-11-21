# Architecture

This document records the layering, the crate boundaries and the reasons for
them. Each section is added by the phase that implements the layer it
describes; the ones marked unwritten have no content yet.

## 1. Layers

```text
        model  ──►  element  ──►  assembly  ──►  solve
          │            │             │            │
          └────────────┴─────────────┴────────────┴──►  results
                                                          │
                                                     report / io
```

The model holds what the user entered: nodes, members, areas, patches,
supports, load cases and combinations. Elements turn that into stiffness
contributions. Assembly numbers the equations and builds the system. The solver
produces displacements, and result recovery produces the forces, reactions and
stresses that reports, tables and the viewport display.

## 2. Crates

| Crate | State | Responsibility | Why it is a separate crate |
|---|---|---|---|
| `nodris-core` | exists | numerical kernel, model schema, elements, analysis, results | Every other crate depends on it, so it must depend on nothing but the standard library. |
| `nodris-io` | planned | canonical JSON, schema, migration, CSV, fingerprints | File formats and migration change on their own schedule, and the kernel should not carry them. |
| `nodris-cli` | planned | scriptable front end over the kernel | Tests can run the command line without a display. |
| `nodris-scene` | planned | renderer-neutral scene description | The viewport can be replaced or tested without depending on analysis types. |
| `nodris-app` | planned | commands and edits, without a UI toolkit | Both the keyboard loop and the window need the command state, so it is stored outside both. |
| `nodris-desktop` | planned | native window over `nodris-scene` and `nodris-app` | Window and platform code is kept away from the solver. |

The planned crates have no directories yet. When one is added, its row changes
from planned to exists in the same commit.

## 3. Dependency rule

Dependencies point one way only:

```text
desktop ─► app ─► core
    │        │       ▲
    └──► scene ──────┘

cli ─► io ─► core
```

`nodris-core` does not refer to a command line, a scene or a window. Both front
ends are written against it. A module in the kernel that imports a scene or UI
type breaks the rule and fails review.

## 4. Units and coordinates

*To be written in Phase 1, with the conversion helpers and their tests.*

The kernel will keep coordinates and forces in one system. Section properties
will keep the units of the table they come from, and conversion will happen
where a section enters the analysis. The table of units goes here once the code
enforces it.

## 5. Performance

*To be written in Phase 5, with the first measurements.*

No timings are recorded yet. The first benchmark will be written here together
with the command that produced it.

## 6. Document format

*To be written in Phase 2, with schema v1.*

The model format will be JSON with stable object ids and a published schema.
Phase 2 adds schema v1. Canonical output, strict and permissive readers, and
migration belong to Phase 10.
