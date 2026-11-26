//! Clamped B-spline spaces: knot vectors, basis values and derivatives.

use crate::math::Real;

/// A knot span narrower than this is treated as empty.
const KNOT_EPS: Real = 1.0e-12;

/// A univariate B-spline space: degree, basis count and open knot vector.
///
/// `knots` holds `n + degree + 1` entries, with the first and last repeated
/// `degree + 1` times, and the parameter domain is
/// `[knots[degree], knots[n]]`. Indices are zero-based throughout.
#[derive(Clone, Debug, PartialEq)]
pub struct BSplineSpace {
    degree: usize,
    n: usize,
    knots: Vec<Real>,
}

impl BSplineSpace {
    /// Clamped open-uniform space on `[0, 1]` with `n` basis functions.
    ///
    /// # Panics
    ///
    /// Panics when `n <= degree`, because such a space has no basis function
    /// with a non-empty support.
    pub fn open_uniform(degree: usize, n: usize) -> Self {
        assert!(
            n > degree,
            "a B-spline space needs at least p + 1 basis functions"
        );
        let mut knots = Vec::with_capacity(n + degree + 1);
        knots.extend(std::iter::repeat_n(0.0, degree + 1));
        let width = 1.0 / (n - degree) as Real;
        for interior in 1..(n - degree) {
            knots.push(interior as Real * width);
        }
        knots.extend(std::iter::repeat_n(1.0, degree + 1));
        Self { degree, n, knots }
    }

    /// Polynomial degree.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Number of basis functions.
    pub fn n(&self) -> usize {
        self.n
    }

    /// The knot vector.
    pub fn knots(&self) -> &[Real] {
        &self.knots
    }

    /// Parametric domain, as a closed interval.
    pub fn domain(&self) -> (Real, Real) {
        (self.knots[self.degree], self.knots[self.n])
    }

    /// Value of the `index`-th basis function at `u`.
    ///
    /// Outside the domain the value is zero. At the two clamped ends the
    /// closed-domain convention applies: the first basis is one at the left
    /// end and the last basis is one at the right end, so a patch corner lands
    /// exactly on its control point.
    pub fn basis(&self, index: usize, u: Real) -> Real {
        if index >= self.n {
            return 0.0;
        }
        let (first, last) = self.domain();
        if u == first {
            return if index == 0 { 1.0 } else { 0.0 };
        }
        if u == last {
            return if index == self.n - 1 { 1.0 } else { 0.0 };
        }
        if u < first || u > last {
            return 0.0;
        }
        basis_rec(&self.knots, self.degree, index, u)
    }

    /// First derivative of the `index`-th basis function at `u`.
    ///
    /// Uses the analytic Cox–de Boor recurrence
    /// `p/(t_{i+p} - t_i) N_{i,p-1} - p/(t_{i+p+1} - t_{i+1}) N_{i+1,p-1}`,
    /// with a zero-width span contributing nothing. This is the derivative the
    /// IGA assembly needs, not a finite difference.
    pub fn basis_derivative(&self, index: usize, u: Real) -> Real {
        let (first, last) = self.domain();
        if self.degree == 0 || index >= self.n || u < first || u > last {
            return 0.0;
        }
        let degree = self.degree;
        let left = if index + degree < self.knots.len() {
            let width = self.knots[index + degree] - self.knots[index];
            if width > KNOT_EPS {
                degree as Real / width * basis_rec(&self.knots, degree - 1, index, u)
            } else {
                0.0
            }
        } else {
            0.0
        };
        let right = if index + 1 < self.n && index + degree + 1 < self.knots.len() {
            let width = self.knots[index + degree + 1] - self.knots[index + 1];
            if width > KNOT_EPS {
                degree as Real / width * basis_rec(&self.knots, degree - 1, index + 1, u)
            } else {
                0.0
            }
        } else {
            0.0
        };
        left - right
    }
}

/// Cox–de Boor recurrence with `0/0 := 0`.
///
/// Degree-zero spans are half-open, `[t_i, t_{i+1})`, so at a repeated knot
/// exactly one base function is active instead of every zero-width span
/// counting. The final non-empty span is closed at the right clamp, which is
/// what makes the recurrence return the left-hand value (and therefore the
/// left-hand derivative) at the end of the domain.
fn basis_rec(knots: &[Real], degree: usize, index: usize, u: Real) -> Real {
    if degree == 0 {
        let Some(&right) = knots.get(index + 1) else {
            return 0.0;
        };
        let inside = knots[index] <= u && u < right;
        let closes_right = right == knots[knots.len() - 1] && (u - right).abs() <= KNOT_EPS;
        return if inside || closes_right { 1.0 } else { 0.0 };
    }
    let left = if index + degree < knots.len() {
        let width = knots[index + degree] - knots[index];
        if width > KNOT_EPS {
            (u - knots[index]) / width * basis_rec(knots, degree - 1, index, u)
        } else {
            0.0
        }
    } else {
        0.0
    };
    let right = if index + degree + 1 < knots.len() {
        let width = knots[index + degree + 1] - knots[index + 1];
        if width > KNOT_EPS {
            (knots[index + degree + 1] - u) / width * basis_rec(knots, degree - 1, index + 1, u)
        } else {
            0.0
        }
    } else {
        0.0
    };
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_quadratic_basis_matches_the_closed_form() {
        let space = BSplineSpace::open_uniform(2, 3);
        assert_eq!(space.domain(), (0.0, 1.0));
        for &u in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let expected = [(1.0 - u) * (1.0 - u), 2.0 * u * (1.0 - u), u * u];
            for (index, &value) in expected.iter().enumerate() {
                let actual = space.basis(index, u);
                assert!(
                    (actual - value).abs() < 1.0e-12,
                    "N_{index}({u}) = {actual}, expected {value}"
                );
            }
        }
    }

    #[test]
    fn the_basis_partitions_to_one() {
        for degree in 0..=3 {
            for n in (degree + 1)..=(degree + 5) {
                let space = BSplineSpace::open_uniform(degree, n);
                for step in 0..=20 {
                    let u = step as Real / 20.0;
                    let sum: Real = (0..n).map(|index| space.basis(index, u)).sum();
                    assert!(
                        (sum - 1.0).abs() < 1.0e-12,
                        "degree {degree}, n {n}: basis sum {sum} at u = {u}"
                    );
                }
            }
        }
    }

    #[test]
    fn the_derivatives_sum_to_zero() {
        for degree in 1..=3 {
            for n in (degree + 1)..=(degree + 5) {
                let space = BSplineSpace::open_uniform(degree, n);
                for step in 0..=20 {
                    let u = step as Real / 20.0;
                    let sum: Real = (0..n).map(|index| space.basis_derivative(index, u)).sum();
                    assert!(
                        sum.abs() < 1.0e-10,
                        "degree {degree}, n {n}: derivative sum {sum} at u = {u}"
                    );
                }
            }
        }
    }

    #[test]
    fn the_derivative_closes_the_right_clamp() {
        // The basis is continuous at the closed right end, so a backward
        // difference there is the reference for the one-sided derivative.
        let space = BSplineSpace::open_uniform(3, 6);
        let step = 1.0e-7;
        for index in 0..space.n() {
            let backward = (space.basis(index, 1.0) - space.basis(index, 1.0 - step)) / step;
            let derivative = space.basis_derivative(index, 1.0);
            assert!(
                (derivative - backward).abs() < 1.0e-5,
                "N'_{index}(1) = {derivative}, backward difference {backward}"
            );
        }
    }

    #[test]
    fn a_derivative_matches_a_finite_difference() {
        let space = BSplineSpace::open_uniform(3, 6);
        let u = 0.37;
        let step = 1.0e-6;
        for index in 0..space.n() {
            let difference =
                (space.basis(index, u + step) - space.basis(index, u - step)) / (2.0 * step);
            let derivative = space.basis_derivative(index, u);
            assert!(
                (derivative - difference).abs() < 1.0e-7,
                "N'_{index}({u}) = {derivative}, finite difference {difference}"
            );
        }
    }
}
