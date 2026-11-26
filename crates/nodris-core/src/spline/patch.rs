//! Tensor-product B-spline patches over a rectangular control net.

use std::fmt;

use super::nurbs::rationalise;
use super::BSplineSpace;
use crate::math::{Real, Vec3};
use crate::model::PatchDirection;

/// Basis values and first parametric derivatives at one point of a patch.
///
/// Every vector is `n_u * n_v` long and row-major in the control net, so
/// `iu * n_v + iv` names control point `(iu, iv)`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PatchBasis {
    /// `N_A(u, v)`.
    pub values: Vec<Real>,
    /// `∂N_A/∂u`.
    pub du: Vec<Real>,
    /// `∂N_A/∂v`.
    pub dv: Vec<Real>,
}

/// Why a control net cannot define a patch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SplineError {
    /// The net is empty or its rows have different lengths.
    RaggedControlNet,
    /// A direction has no more control points than its degree.
    TooFewControlPoints {
        /// The direction that failed.
        direction: PatchDirection,
        /// Control points available.
        count: usize,
        /// Requested degree.
        degree: usize,
    },
    /// The weight count differs from the number of control points.
    WeightCount {
        /// Weights supplied.
        count: usize,
        /// Control points in the net.
        expected: usize,
    },
    /// A weight is zero, negative or not finite.
    InvalidWeight {
        /// Flat control-point index, row-major.
        index: usize,
    },
}

impl fmt::Display for SplineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RaggedControlNet => {
                write!(f, "control net is empty or its rows differ in length")
            }
            Self::TooFewControlPoints {
                direction,
                count,
                degree,
            } => write!(
                f,
                "{} direction has {count} control points for degree {degree}",
                direction.as_str()
            ),
            Self::WeightCount { count, expected } => {
                write!(f, "{count} weights for {expected} control points")
            }
            Self::InvalidWeight { index } => {
                write!(f, "weight {index} is not finite and positive")
            }
        }
    }
}

impl std::error::Error for SplineError {}

/// A clamped tensor-product B-spline patch.
#[derive(Clone, Debug, PartialEq)]
pub struct BSplinePatch {
    u: BSplineSpace,
    v: BSplineSpace,
    control: Vec<Vec<Vec3>>,
    /// One weight per control point, row-major, or `None` for a polynomial
    /// patch.
    weights: Option<Vec<Real>>,
}

impl BSplinePatch {
    /// Builds a patch from a rectangular control net and the two degrees.
    ///
    /// Each direction uses the clamped open-uniform knots from
    /// [`BSplineSpace::open_uniform`].
    pub fn from_grid(
        control: &[Vec<Vec3>],
        degree_u: usize,
        degree_v: usize,
    ) -> Result<Self, SplineError> {
        let columns = control.first().ok_or(SplineError::RaggedControlNet)?.len();
        if columns == 0 || control.iter().any(|row| row.len() != columns) {
            return Err(SplineError::RaggedControlNet);
        }
        if control.len() <= degree_u {
            return Err(SplineError::TooFewControlPoints {
                direction: PatchDirection::U,
                count: control.len(),
                degree: degree_u,
            });
        }
        if columns <= degree_v {
            return Err(SplineError::TooFewControlPoints {
                direction: PatchDirection::V,
                count: columns,
                degree: degree_v,
            });
        }
        Ok(Self {
            u: BSplineSpace::open_uniform(degree_u, control.len()),
            v: BSplineSpace::open_uniform(degree_v, columns),
            control: control.to_vec(),
            weights: None,
        })
    }

    /// Builds a rational patch: the same clamped spaces and control net, with
    /// one weight per control point.
    ///
    /// Weights are row-major like the net, `iu * n_v + iv`, and must be finite
    /// and positive. A zero or negative weight would let the weighted sum
    /// vanish inside the patch, and the rational point with it.
    pub fn from_grid_weighted(
        control: &[Vec<Vec3>],
        weights: &[Real],
        degree_u: usize,
        degree_v: usize,
    ) -> Result<Self, SplineError> {
        let mut patch = Self::from_grid(control, degree_u, degree_v)?;
        if weights.len() != patch.n_cp() {
            return Err(SplineError::WeightCount {
                count: weights.len(),
                expected: patch.n_cp(),
            });
        }
        if let Some(index) = weights.iter().position(|&w| !(w.is_finite() && w > 0.0)) {
            return Err(SplineError::InvalidWeight { index });
        }
        patch.weights = Some(weights.to_vec());
        Ok(patch)
    }

    /// Basis in the `u` direction.
    pub fn space_u(&self) -> &BSplineSpace {
        &self.u
    }

    /// Basis in the `v` direction.
    pub fn space_v(&self) -> &BSplineSpace {
        &self.v
    }

    /// Control points along `u`.
    pub fn n_u(&self) -> usize {
        self.u.n()
    }

    /// Control points along `v`.
    pub fn n_v(&self) -> usize {
        self.v.n()
    }

    /// Total control points, `n_u * n_v`.
    pub fn n_cp(&self) -> usize {
        self.n_u() * self.n_v()
    }

    /// The control net, one row of control points per `u` index.
    pub fn control_net(&self) -> &[Vec<Vec3>] {
        &self.control
    }

    /// Parametric domain `((u0, u1), (v0, v1))`.
    pub fn domain(&self) -> ((Real, Real), (Real, Real)) {
        (self.u.domain(), self.v.domain())
    }

    /// Basis values and first parametric derivatives at `(u, v)`.
    ///
    /// Rational when the patch carries weights, polynomial otherwise. The
    /// rational form is `R_A = N_A w_A / W` with `W = Σ N_B w_B`, and its
    /// derivative is the quotient rule applied to that ratio.
    pub fn basis(&self, u: Real, v: Real) -> PatchBasis {
        let n_u = self.n_u();
        let n_v = self.n_v();
        // Every row of the net shares the same `v` basis, so evaluate it once
        // rather than once per row.
        let basis_v: Vec<Real> = (0..n_v).map(|iv| self.v.basis(iv, v)).collect();
        let derivative_v: Vec<Real> = (0..n_v).map(|iv| self.v.basis_derivative(iv, v)).collect();
        let mut basis = PatchBasis {
            values: vec![0.0; n_u * n_v],
            du: vec![0.0; n_u * n_v],
            dv: vec![0.0; n_u * n_v],
        };
        for iu in 0..n_u {
            let nu = self.u.basis(iu, u);
            let dnu = self.u.basis_derivative(iu, u);
            let row = iu * n_v;
            for (offset, (&nv, &dnv)) in basis_v.iter().zip(&derivative_v).enumerate() {
                let index = row + offset;
                basis.values[index] = nu * nv;
                basis.du[index] = dnu * nv;
                basis.dv[index] = nu * dnv;
            }
        }
        if let Some(weights) = &self.weights {
            rationalise(&mut basis, weights);
        }
        basis
    }

    /// The patch map `x(u, v)`.
    pub fn point(&self, u: Real, v: Real) -> Vec3 {
        let basis = self.basis(u, v);
        self.combine(&basis.values)
    }

    /// The two parametric tangents `(∂x/∂u, ∂x/∂v)`.
    pub fn tangents(&self, u: Real, v: Real) -> (Vec3, Vec3) {
        let basis = self.basis(u, v);
        (self.combine(&basis.du), self.combine(&basis.dv))
    }

    fn combine(&self, weights: &[Real]) -> Vec3 {
        self.control
            .iter()
            .flatten()
            .zip(weights)
            .fold(Vec3::ZERO, |sum, (&point, &weight)| sum + point * weight)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bilinear_net() -> Vec<Vec<Vec3>> {
        vec![
            vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 2.0, 0.0)],
            vec![Vec3::new(3.0, 0.0, 0.0), Vec3::new(3.0, 2.0, 0.0)],
        ]
    }

    #[test]
    fn a_bilinear_patch_reproduces_its_control_grid() {
        let patch = BSplinePatch::from_grid(&bilinear_net(), 1, 1).expect("a 2 by 2 net");
        assert_eq!(patch.domain(), ((0.0, 1.0), (0.0, 1.0)));
        assert_eq!(patch.point(0.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
        assert_eq!(patch.point(1.0, 0.0), Vec3::new(3.0, 0.0, 0.0));
        assert_eq!(patch.point(0.0, 1.0), Vec3::new(0.0, 2.0, 0.0));
        assert_eq!(patch.point(1.0, 1.0), Vec3::new(3.0, 2.0, 0.0));

        let centre = patch.point(0.5, 0.5);
        assert!((centre - Vec3::new(1.5, 1.0, 0.0)).norm() < 1.0e-12);

        let (du, dv) = patch.tangents(0.3, 0.7);
        assert!((du - Vec3::new(3.0, 0.0, 0.0)).norm() < 1.0e-12);
        assert!((dv - Vec3::new(0.0, 2.0, 0.0)).norm() < 1.0e-12);
    }

    #[test]
    fn the_patch_basis_partitions_to_one() {
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
        let patch = BSplinePatch::from_grid(&net, 2, 2).expect("a 3 by 3 net");
        for &u in &[0.0, 0.2, 0.5, 0.8, 1.0] {
            for &v in &[0.0, 0.35, 1.0] {
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
    fn a_ragged_or_undersized_net_is_rejected() {
        let empty: Vec<Vec<Vec3>> = Vec::new();
        assert_eq!(
            BSplinePatch::from_grid(&empty, 1, 1),
            Err(SplineError::RaggedControlNet)
        );

        let ragged = vec![vec![Vec3::ZERO], vec![Vec3::ZERO, Vec3::X]];
        assert_eq!(
            BSplinePatch::from_grid(&ragged, 1, 1),
            Err(SplineError::RaggedControlNet)
        );

        let single = vec![vec![Vec3::ZERO]];
        assert!(matches!(
            BSplinePatch::from_grid(&single, 1, 1),
            Err(SplineError::TooFewControlPoints {
                direction: PatchDirection::U,
                ..
            })
        ));
    }
}
