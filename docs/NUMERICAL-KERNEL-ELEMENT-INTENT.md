# Element contract

This document records the contract between an element formulation and the
analysis pipeline in `nodris-core` as it stands at the end of Phase 3. It
covers the degrees of freedom, the local axes, the sign conventions and what a
`None` geometric stiffness means. The solver contract lives in
[Solver guarantees and limits](SOLVER-GUARANTEES.md).

The implementation is in
[`element/mod.rs`](../crates/nodris-core/src/element/mod.rs), with the frame in
[`element/frame3d.rs`](../crates/nodris-core/src/element/frame3d.rs) and the bar
in [`element/truss3d.rs`](../crates/nodris-core/src/element/truss3d.rs).

---

## Degrees of freedom

Every node has six degrees of freedom in one fixed order, listed once in
`DOF_NAMES`:

| Index | Name | Meaning |
|---|---|---|
| 0 | tx | translation along global X |
| 1 | ty | translation along global Y |
| 2 | tz | translation along global Z |
| 3 | rx | rotation about global X |
| 4 | ry | rotation about global Y |
| 5 | rz | rotation about global Z |

`node_dofs` and `member_dofs` are the only places that turn a node index into
global numbers. Node *n* owns global degrees of freedom `6n` through `6n + 5`.

An element reports its own degrees of freedom in `dofs()`, in local order. A
frame has twelve: the six of the start node, then the six of the end node. A
bar also has twelve, because a node keeps six degrees of freedom whether or not
the bar resists them.

## Local axes

`local_axes` returns the direction cosines of a member as the rows of a 3 × 3
matrix.

- Local x runs from the first node to the second.
- Local z stays as close to global Z as the member allows.
- Local y completes the right-handed triad.
- A vertical member has no unique near-Z direction, so the reference becomes
  global X.
- `rotation` spins the section about the member axis, which is how a rotated
  I-beam is modelled.

The member direction and the rotation must both be finite. `local_axes` returns
`None` otherwise, and `FrameElement::new` turns that into
`ElementError::DegenerateMember`.

`axes()` hands the direction cosines to the assembler, so that end forces can
be transformed for post-processing without every consumer re-deriving them.

## Stiffness

`stiffness()` returns the global stiffness, row-major, `n × n` with
`n = dofs().len()`.

A frame builds its matrix in local axes from four rigidities: `EA` in kN, `EIy`
and `EIz` in kN·m², and `GJ` in kN·m². Deflection in local y is resisted by
`EIz` with rotation about local z, and deflection in local z is resisted by
`EIy` with rotation about local y. The local matrix reaches global axes through
the congruence `Tᵀ K T`, where `T` carries the direction cosines in four 3 × 3
blocks, one per translation or rotation group.

A bar builds only the axial term, `(EA/L) c cᵀ`, in the translational blocks.
Its rotational diagonal stays zero. The assembler deactivates those degrees of
freedom rather than giving them a penalty stiffness, because a penalty hides
the problem instead of reporting it.

## Sign conventions

End forces come back in local axes, and tension is positive.

`end_forces` computes `f = k·d − p`. The `d` vector holds global-axis
displacements in local degree-of-freedom order, and `p` is the equivalent nodal
load of the distributed load. The subtraction is what makes the result correct
for a loaded member.

For a frame under a uniform load `q` in local axes, the equivalent nodal load
is `qL/2` at each end in the load direction, plus fixed-end couples of
`qL²/12`. The couple at the start node has the opposite sign to the couple at
the end node. The signs follow the rotation convention used by the geometric
stiffness: `rz` is `dv/dx` and `ry` is `-dw/dx`.

`equivalent_load` returns that same equivalent load in global axes, so an
assembler can add it straight to the global load vector.

A bar returns axial end forces only, `[-N, 0, ..., +N, 0, ...]` in local axes,
with `N` positive in tension.

## Units

| Quantity | Unit |
|---|---|
| Force | kN |
| Length | m |
| Mass | t |
| Stiffness | kN/m for translation, kN·m for rotation |
| Distributed load | kN/m |
| Pressure | kN/m² |

Mass is in tonnes so that `kN / (m·t) = s⁻²` for modal analysis.

## Lumped mass

`lumped_mass` returns mass in local degree-of-freedom order, with half the
mass at each end.

A frame includes rotary inertia on the rotational degrees of freedom even though
Euler–Bernoulli theory drops it: the eigenproblem needs a positive mass matrix,
and the true value is far too small to move the bending modes. A bar leaves
rotary inertia out, because its rotational degrees of freedom are deactivated.

## Geometric stiffness

`geometric_stiffness` takes an `InitialStress` and returns the global
geometric stiffness, or `None`.

`None` means the formulation has no initial-stress matrix. A second-order
analysis leaves those degrees of freedom untouched, and it must not substitute
a zero matrix, because a zero matrix changes the tangent stiffness without
saying so.

`InitialStress` carries an axial force in kN and the membrane resultants `nx`,
`ny` and `nxy` in kN/m. Tension is positive, so compression subtracts from the
elastic stiffness and softens the transverse degrees of freedom. That softening
is what produces P-Delta amplification and buckling.

A frame reads the axial force only. The membrane components belong to shells
and patches.

## What the contract does not cover

- Only the shell membrane is implemented. The Mindlin plate bending, the
  shell geometric stiffness and the isogeometric patch arrive in later steps
  of Phase 4, and they reuse this degree-of-freedom order, so that ties and
  loads can be shared between the formulations.
- The assembler does not exist yet. Numbering, restraints and the deactivation
  of zero-diagonal degrees of freedom are its job.
- Strong-axis bending and warping are not modelled. The section carries a
  warping constant for the viewer and later use, but no term in the frame
  matrix uses it.

## What the tests pin

| Rule | Where |
|---|---|
| Local matrix symmetric, rigid-body modes free | `element/frame3d.rs` |
| Cantilever deflection, fixed-fixed end moments, method of joints | `tests/arch_compare.rs` |
| Compressed cantilever reaches the Euler load | `element/frame3d.rs` |
| A bar carries axial force only | `element/truss3d.rs` |
| Mass split shared by frame and bar | `element/mod.rs` |
| Membrane constant strain and edge traction | `element/shell4.rs` |
