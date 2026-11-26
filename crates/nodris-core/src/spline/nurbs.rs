//! Rational (NURBS) evaluation on top of the polynomial patch.
//!
//! A NURBS patch stores one weight per control point and evaluates
//! `R_A = N_A w_A / W`, where `W = Σ N_B w_B` is the weighted sum of the whole
//! polynomial basis. The derivatives follow from the quotient rule. A patch
//! without weights never enters this module, which keeps the polynomial path
//! free of the division.

use super::PatchBasis;
use crate::math::Real;

/// Converts polynomial basis values and parametric derivatives to rational ones
/// in place.
///
/// `weights` is row-major in the control net, the same order as [`PatchBasis`].
/// Every entry is finite and positive, which
/// [`BSplinePatch::from_grid_weighted`](super::BSplinePatch::from_grid_weighted)
/// checks when the patch is built, so `W` stays positive and the quotient is
/// defined over the whole patch.
pub(crate) fn rationalise(basis: &mut PatchBasis, weights: &[Real]) {
    let mut denominator = 0.0;
    let mut denominator_du = 0.0;
    let mut denominator_dv = 0.0;
    for (((&value, &du), &dv), &weight) in basis
        .values
        .iter()
        .zip(&basis.du)
        .zip(&basis.dv)
        .zip(weights)
    {
        denominator += value * weight;
        denominator_du += du * weight;
        denominator_dv += dv * weight;
    }
    let inverse = 1.0 / denominator;
    let inverse_squared = inverse * inverse;

    for (((value, du), dv), &weight) in basis
        .values
        .iter_mut()
        .zip(&mut basis.du)
        .zip(&mut basis.dv)
        .zip(weights)
    {
        let weighted = *value * weight;
        *value = weighted * inverse;
        *du = (*du * weight * denominator - weighted * denominator_du) * inverse_squared;
        *dv = (*dv * weight * denominator - weighted * denominator_dv) * inverse_squared;
    }
}

#[cfg(test)]
mod tests {
    use super::super::{BSplinePatch, SplineError};
    use crate::math::{Real, Vec3};

    /// A quarter cylinder of radius one: the middle control point sits at the
    /// corner with weight `1/√2`, which is what turns the two quadratic
    /// segments into an arc. The `v` direction is a straight extrusion.
    fn quarter_cylinder() -> BSplinePatch {
        let weight = std::f64::consts::FRAC_1_SQRT_2;
        let control = vec![
            vec![Vec3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 1.0)],
            vec![Vec3::new(1.0, 1.0, 0.0), Vec3::new(1.0, 1.0, 1.0)],
            vec![Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 1.0, 1.0)],
        ];
        let weights = vec![1.0, 1.0, weight, weight, 1.0, 1.0];
        BSplinePatch::from_grid_weighted(&control, &weights, 2, 1).expect("a valid rational net")
    }

    #[test]
    fn a_rational_arc_keeps_its_exact_radius() {
        let patch = quarter_cylinder();
        for step in 0..=20 {
            let u = step as Real / 20.0;
            for &v in &[0.0, 0.5, 1.0] {
                let point = patch.point(u, v);
                let radius = (point.x * point.x + point.y * point.y).sqrt();
                assert!(
                    (radius - 1.0).abs() < 1.0e-12,
                    "radius {radius} at ({u}, {v})"
                );
                assert!(
                    (point.z - v).abs() < 1.0e-12,
                    "height {} at v = {v}",
                    point.z
                );
            }
        }
    }

    #[test]
    fn rational_derivatives_match_a_finite_difference() {
        let patch = quarter_cylinder();
        let step = 1.0e-6;
        for &u in &[0.15, 0.4, 0.65, 0.9] {
            // Interior parameters only: a central difference at the domain
            // edge would evaluate the basis outside the patch.
            for &v in &[0.2, 0.5, 0.8] {
                let (du, dv) = patch.tangents(u, v);
                let difference_u =
                    (patch.point(u + step, v) - patch.point(u - step, v)) / (2.0 * step);
                let difference_v =
                    (patch.point(u, v + step) - patch.point(u, v - step)) / (2.0 * step);
                assert!(
                    (du - difference_u).norm() < 1.0e-7,
                    "du {du:?} against {difference_u:?} at ({u}, {v})"
                );
                assert!(
                    (dv - difference_v).norm() < 1.0e-7,
                    "dv {dv:?} against {difference_v:?} at ({u}, {v})"
                );
            }
        }
    }

    #[test]
    fn the_rational_basis_partitions_to_one() {
        let patch = quarter_cylinder();
        for step in 0..=10 {
            let u = step as Real / 10.0;
            for &v in &[0.0, 0.25, 1.0] {
                let basis = patch.basis(u, v);
                let values: Real = basis.values.iter().sum();
                let du: Real = basis.du.iter().sum();
                let dv: Real = basis.dv.iter().sum();
                assert!((values - 1.0).abs() < 1.0e-12, "sum {values} at ({u}, {v})");
                assert!(du.abs() < 1.0e-10, "du sum {du} at ({u}, {v})");
                assert!(dv.abs() < 1.0e-10, "dv sum {dv} at ({u}, {v})");
            }
        }
    }

    #[test]
    fn unit_weights_reproduce_the_polynomial_patch() {
        let net = vec![
            vec![
                Vec3::ZERO,
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 2.0, 0.0),
            ],
            vec![
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 1.0),
                Vec3::new(1.0, 2.0, 0.0),
            ],
            vec![
                Vec3::new(2.0, 0.0, 0.0),
                Vec3::new(2.0, 1.0, 0.0),
                Vec3::new(2.0, 2.0, 0.0),
            ],
        ];
        let polynomial = BSplinePatch::from_grid(&net, 2, 2).expect("a 3 by 3 net");
        let weighted =
            BSplinePatch::from_grid_weighted(&net, &[1.0; 9], 2, 2).expect("unit weights");

        for &u in &[0.0, 0.3, 0.5, 0.9, 1.0] {
            for &v in &[0.0, 0.6, 1.0] {
                assert!((polynomial.point(u, v) - weighted.point(u, v)).norm() < 1.0e-13);
                let (pu, pv) = polynomial.tangents(u, v);
                let (wu, wv) = weighted.tangents(u, v);
                assert!((pu - wu).norm() < 1.0e-13);
                assert!((pv - wv).norm() < 1.0e-13);
            }
        }
    }

    #[test]
    fn weights_must_match_the_net_and_be_positive() {
        let net = vec![vec![Vec3::ZERO, Vec3::X], vec![Vec3::Y, Vec3::Z]];
        assert_eq!(
            BSplinePatch::from_grid_weighted(&net, &[1.0, 1.0, 1.0], 1, 1),
            Err(SplineError::WeightCount {
                count: 3,
                expected: 4,
            })
        );
        for &bad in &[0.0, -1.0, f64::NAN] {
            assert!(
                matches!(
                    BSplinePatch::from_grid_weighted(&net, &[1.0, bad, 1.0, 1.0], 1, 1),
                    Err(SplineError::InvalidWeight { index: 1 })
                ),
                "weight {bad} was accepted"
            );
        }
    }
}
