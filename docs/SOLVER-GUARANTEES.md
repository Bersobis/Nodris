# Solver guarantees and limits

This document records the contract of the `nodris-core` banded solver as it
exists at Phase 1. It is a statement about the numerical kernel, not about a
finished structural product: no element formulations, loads, reports or design
checks are guaranteed by it.

The implementation is in
[`math/linalg.rs`](../crates/nodris-core/src/math/linalg.rs), with typed errors
in [`math/error.rs`](../crates/nodris-core/src/math/error.rs), thresholds in
[`math/tolerance.rs`](../crates/nodris-core/src/math/tolerance.rs) and the
unit-boundary rules in [`units.rs`](../crates/nodris-core/src/units.rs).

---

## What the solver guarantees

### Input contract

The solver accepts a symmetric positive-definite (SPD) matrix in
[`BandedMatrix`](../crates/nodris-core/src/math/linalg.rs) storage.

- Symmetry is a caller obligation. Only the lower triangle is stored; use
  `add_symmetric` once per symmetric pair. Adding both `(i, j)` and `(j, i)`
  doubles the off-diagonal and no later check can detect it.
- Every contribution must be finite. `add_symmetric` panics in every build
  profile on a non-finite value or an entry outside the stored band. A dropped
  stiffness term would produce plausible but wrong displacements, so the write
  is refused rather than ignored.
- The matrix must be SPD. Restraints and zero-diagonal deactivation are the
  assembler's responsibility in a later phase; the solver does not repair a
  singular stiffness by adding penalties.

### Factorisation

`BandedMatrix::factor` consumes the matrix and performs an unpivoted
`L·D·Lᵀ` factorisation. On success it returns a `Factored` that owns both the
factor and an unfactored copy of the matrix; the unfactored copy exists only for
the residual check. A stiffness change therefore requires a new factor, and a
stale factor cannot be reused accidentally.

Every pivot is checked before it is used:

| Condition | Typed failure |
|---|---|
| Pivot is not finite | `SolveError::NotFinite` |
| Pivot is negative | `SolveError::NotPositiveDefinite` |
| Pivot is at or below `PIVOT_RELATIVE × original diagonal` | `SolveError::Singular` |

`PIVOT_RELATIVE` is `1e-12`. It is relative to the original diagonal of the
same equation, so uniformly rescaling an SPD system does not change the verdict.
All error variants carry the failing equation number; callers receive a typed
value, not a string or a NaN.

### Solve

Only `Factored::solve` can solve a system. It overwrites the right-hand side
with the solution and guarantees:

1. every right-hand-side value is checked for finiteness before substitution;
2. forward, diagonal and backward substitution run on the stored band;
3. the relative residual `‖A·x − b‖∞ / ‖b‖∞` is checked against `RESIDUAL`
   (`1e-8`) using the retained unfactored matrix (absolute residual when
   `b = 0`);
4. a residual above the threshold returns `SolveError::ResidualTooLarge` with
   the worst equation and the measured value.

A successful solve therefore returns finite values whose backward residual is
within the named threshold. A non-finite load cannot put NaNs into a result.

### Verification evidence

The guarantees above are exercised by the unit tests in
`crates/nodris-core/src/math/linalg.rs`:

- `factor_and_solve_reproduce_the_closed_form_3x3_system`: a hand-solved
  tridiagonal system;
- `solve_agrees_with_the_dense_oracle_on_seeded_spd_systems`: seeded SPD
  systems compared with an independent dense solver;
- `zero_diagonal_is_reported_as_singular`, `negative_pivot_is_reported_as_not_positive_definite`,
  `non_finite_pivot_is_reported_as_not_finite` and
  `pivot_near_the_relative_threshold_is_rejected`: typed failure paths;
- `pivot_rejection_does_not_depend_on_the_units`: scale invariance;
- `relative_residual_does_not_depend_on_the_scale_of_the_load` and
  `excessive_relative_residual_is_reported`: residual contract.

---

## What the solver does not guarantee

- No pivoting. There is no row or column permutation, so indefinite and
  semidefinite matrices are refused, not solved. A negative pivot means the
  input is outside the contract; it is not a signal to try a different method
  inside this solver.
- No iterative refinement. One triangular solve is followed by a residual
  check, not by correction iterations. A residual below `1e-8` is a backward
  check, not a forward-error bound and not a condition-number estimate.
- No condition estimate. An ill-conditioned but SPD matrix can pass the
  pivot and residual checks while still returning a solution with large forward
  error. The solver does not quantify that risk.
- No general sparse support. The matrix must fit the fixed half-bandwidth
  given at construction. Entries outside it are errors, not implicit zeros to
  be discovered later.
- No unsymmetric or complex systems. Storage and factorisation assume real
  symmetric coefficients.
- No model semantics. The solver knows nothing about nodes, elements,
  supports, load cases or units. Section-table values must already have crossed
  the conversion boundary in [`units.rs`](../crates/nodris-core/src/units.rs);
  analysis quantities stay in kN, m and kN/m².
- No product-level claims. Passing this kernel contract says nothing about
  element formulations, equilibrium of an assembled structure, design-check
  compliance or engineering correctness. Those claims belong to their own
  phases and their own verification.

---

## Contract changes

Changing a pivot rule, a residual threshold, the error variants or the band
layout changes this public contract. Update this document and the tests in the
same commit, and name the change in the commit message.
