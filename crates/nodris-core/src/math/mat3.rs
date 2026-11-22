//! Three-by-three matrix used for local axes and element transformations.

use super::{Real, Vec3, EPS};

/// A row-major 3x3 matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mat3 {
    rows: [[Real; 3]; 3],
}

impl Mat3 {
    pub const ZERO: Self = Self {
        rows: [[0.0; 3]; 3],
    };

    pub const IDENTITY: Self = Self {
        rows: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    };

    pub const fn from_rows(rows: [[Real; 3]; 3]) -> Self {
        Self { rows }
    }

    /// Row `index`, copied out. Rows are the rows of the matrix as written.
    ///
    /// # Panics
    ///
    /// Panics on an index above 2.
    pub fn row(self, index: usize) -> [Real; 3] {
        self.rows[index]
    }

    pub fn col(self, index: usize) -> [Real; 3] {
        [
            self.rows[0][index],
            self.rows[1][index],
            self.rows[2][index],
        ]
    }

    pub fn transpose(self) -> Self {
        let r = self.rows;
        Self::from_rows([
            [r[0][0], r[1][0], r[2][0]],
            [r[0][1], r[1][1], r[2][1]],
            [r[0][2], r[1][2], r[2][2]],
        ])
    }

    pub fn mul_vec(self, v: Vec3) -> Vec3 {
        Vec3::new(
            self.rows[0][0] * v.x + self.rows[0][1] * v.y + self.rows[0][2] * v.z,
            self.rows[1][0] * v.x + self.rows[1][1] * v.y + self.rows[1][2] * v.z,
            self.rows[2][0] * v.x + self.rows[2][1] * v.y + self.rows[2][2] * v.z,
        )
    }

    pub fn mul(&self, other: &Self) -> Self {
        let mut rows = [[0.0; 3]; 3];
        for (i, row) in rows.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value = self.rows[i][0] * other.rows[0][j]
                    + self.rows[i][1] * other.rows[1][j]
                    + self.rows[i][2] * other.rows[2][j];
            }
        }
        Self::from_rows(rows)
    }

    /// `self^T * other * self`, with `self` as the transformation.
    ///
    /// The receiver is the transform, so element code reads
    /// `axes.congruence(&block)` exactly as the Stodola formulations do, and a
    /// ported call keeps its shape.
    pub fn congruence(&self, other: &Self) -> Self {
        self.transpose().mul(other).mul(self)
    }

    pub fn det(self) -> Real {
        let r = self.rows;
        r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
    }

    /// Inverse, or `None` when the matrix is singular or not finite.
    ///
    /// Uses the adjugate, which is exact for the 3x3 case and needs no
    /// pivoting. A determinant at or below [`EPS`] is rejected rather than
    /// divided by, because the resulting entries would be meaningless.
    pub fn inverse(self) -> Option<Self> {
        let det = self.det();
        if !det.is_finite() || det.abs() <= EPS {
            return None;
        }
        let r = self.rows;
        let adjugate = Self::from_rows([
            [
                r[1][1] * r[2][2] - r[1][2] * r[2][1],
                r[0][2] * r[2][1] - r[0][1] * r[2][2],
                r[0][1] * r[1][2] - r[0][2] * r[1][1],
            ],
            [
                r[1][2] * r[2][0] - r[1][0] * r[2][2],
                r[0][0] * r[2][2] - r[0][2] * r[2][0],
                r[0][2] * r[1][0] - r[0][0] * r[1][2],
            ],
            [
                r[1][0] * r[2][1] - r[1][1] * r[2][0],
                r[0][1] * r[2][0] - r[0][0] * r[2][1],
                r[0][0] * r[1][1] - r[0][1] * r[1][0],
            ],
        ]);
        let mut rows = [[0.0; 3]; 3];
        for (i, row) in rows.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value = adjugate.rows[i][j] / det;
            }
        }
        Some(Self::from_rows(rows))
    }
}

impl Default for Mat3 {
    fn default() -> Self {
        Self::ZERO
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rotation_z(angle: Real) -> Mat3 {
        let (s, c) = angle.sin_cos();
        Mat3::from_rows([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    }

    #[test]
    fn identity_leaves_a_vector_alone() {
        let v = Vec3::new(1.0, -2.0, 3.0);
        assert_eq!(Mat3::IDENTITY.mul_vec(v), v);
        assert_eq!(Mat3::IDENTITY.det(), 1.0);
    }

    #[test]
    fn rows_and_columns_are_read_in_their_own_order() {
        let m = Mat3::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
        assert_eq!(m.row(1), [4.0, 5.0, 6.0]);
        assert_eq!(m.col(1), [2.0, 5.0, 8.0]);
        assert_eq!(m.row(2), [7.0, 8.0, 9.0]);
    }

    #[test]
    fn transpose_of_a_product_reverses_the_factors() {
        let a = rotation_z(0.3);
        let b = Mat3::from_rows([[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]);
        assert_eq!(a.mul(&b).transpose(), b.transpose().mul(&a.transpose()));
    }

    #[test]
    fn inverse_undoes_the_transform() {
        let m = rotation_z(0.7).mul(&Mat3::from_rows([
            [2.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 4.0],
        ]));
        let inverse = m.inverse().expect("invertible");
        let product = m.mul(&inverse);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((product.row(i)[j] - expected).abs() < 1.0e-12);
            }
        }
    }

    #[test]
    fn singular_matrix_has_no_inverse() {
        let rank_deficient = Mat3::from_rows([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [0.0, 0.0, 0.0]]);
        assert_eq!(rank_deficient.det(), 0.0);
        assert!(rank_deficient.inverse().is_none());
    }

    #[test]
    fn empty_matrix_is_singular() {
        assert!(Mat3::ZERO.inverse().is_none());
    }

    #[test]
    fn congruence_applies_the_transform_on_both_sides() {
        // A transform that swaps x and y swaps rows and columns of the matrix.
        // Worked by hand, so this fails if the receiver and the argument are
        // ever exchanged.
        let t = Mat3::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let k = Mat3::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
        let expected = Mat3::from_rows([[5.0, 4.0, 6.0], [2.0, 1.0, 3.0], [8.0, 7.0, 9.0]]);
        assert_eq!(t.congruence(&k), expected);
    }

    #[test]
    fn congruence_matches_the_explicit_product() {
        let k = Mat3::from_rows([[2.0, 1.0, 0.0], [1.0, 3.0, 0.5], [0.0, 0.5, 4.0]]);
        let t = rotation_z(-0.4);
        let expected = t.transpose().mul(&k).mul(&t);
        assert_eq!(t.congruence(&k), expected);
        // A rotation preserves the determinant of the transformed matrix.
        assert!((t.congruence(&k).det() - k.det()).abs() < 1.0e-12);
    }
}
