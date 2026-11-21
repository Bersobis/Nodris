//! Three-dimensional vector used for coordinates, directions and lengths.

use core::ops::{Add, Div, Index, Mul, Neg, Sub};

use super::{Real, EPS};

/// A point or direction in three dimensions.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: Real,
    pub y: Real,
    pub z: Real,
}

impl Vec3 {
    pub const ZERO: Self = Self {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    pub const X: Self = Self::new(1.0, 0.0, 0.0);
    pub const Y: Self = Self::new(0.0, 1.0, 0.0);
    pub const Z: Self = Self::new(0.0, 0.0, 1.0);

    pub const fn new(x: Real, y: Real, z: Real) -> Self {
        Self { x, y, z }
    }

    pub const fn from_array(values: [Real; 3]) -> Self {
        Self::new(values[0], values[1], values[2])
    }

    pub const fn to_array(self) -> [Real; 3] {
        [self.x, self.y, self.z]
    }

    pub fn dot(self, other: Self) -> Real {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Right-handed cross product: `self x other`.
    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    pub fn norm_squared(self) -> Real {
        self.dot(self)
    }

    /// Euclidean length.
    pub fn norm(self) -> Real {
        self.norm_squared().sqrt()
    }

    /// Unit vector in the same direction.
    ///
    /// Returns `None` for a zero-length vector or one with a non-finite
    /// component: both have no direction to hand back, and a caller that
    /// silently receives a NaN axis produces wrong element matrices.
    pub fn normalized(self) -> Option<Self> {
        let norm = self.norm();
        if !norm.is_finite() || norm <= EPS {
            return None;
        }
        Some(self / norm)
    }

    pub fn distance(self, other: Self) -> Real {
        (self - other).norm()
    }

    /// True when the vector is shorter than [`EPS`], the length at which the
    /// kernel treats two points as coincident.
    pub fn is_zero(self) -> bool {
        self.norm() < EPS
    }

    /// Rotate about an arbitrary axis by `angle` radians (Rodrigues' formula).
    ///
    /// Returns `None` when the axis has no direction, the same rule
    /// [`Vec3::normalized`] applies.
    pub fn rotated(self, axis: Self, angle: Real) -> Option<Self> {
        let k = axis.normalized()?;
        let (s, c) = angle.sin_cos();
        Some(self * c + k.cross(self) * s + k * (k.dot(self) * (1.0 - c)))
    }
}

impl Add for Vec3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Neg for Vec3 {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl Mul<Real> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: Real) -> Self {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Div<Real> for Vec3 {
    type Output = Self;

    fn div(self, rhs: Real) -> Self {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl Index<usize> for Vec3 {
    type Output = Real;

    /// Component by axis index: 0 for x, 1 for y, 2 for z.
    ///
    /// # Panics
    ///
    /// Panics on an index above 2.
    fn index(&self, index: usize) -> &Real {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Vec3 index {index} is outside 0..3"),
        }
    }
}

impl From<[Real; 3]> for Vec3 {
    fn from(values: [Real; 3]) -> Self {
        Self::from_array(values)
    }
}

impl From<Vec3> for [Real; 3] {
    fn from(v: Vec3) -> Self {
        v.to_array()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dot_and_norm_agree_with_pythagoras() {
        let v = Vec3::new(3.0, 0.0, 4.0);
        assert_eq!(v.dot(Vec3::new(0.0, 1.0, 0.0)), 0.0);
        assert_eq!(v.norm_squared(), 25.0);
        assert_eq!(v.norm(), 5.0);
    }

    #[test]
    fn cross_product_follows_the_right_hand_rule() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        assert_eq!(x.cross(y), Vec3::new(0.0, 0.0, 1.0));
        assert_eq!(y.cross(x), Vec3::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn cross_product_is_orthogonal_to_both_inputs() {
        let a = Vec3::new(2.0, -1.0, 0.5);
        let b = Vec3::new(-3.0, 4.0, 1.0);
        let c = a.cross(b);
        assert!(c.dot(a).abs() < 1.0e-12);
        assert!(c.dot(b).abs() < 1.0e-12);
        // |a x b| = |a||b| sin(theta), checked through the scalar triple
        // product identity |c|^2 = |a|^2|b|^2 - (a.b)^2.
        let expected = a.norm_squared() * b.norm_squared() - a.dot(b).powi(2);
        assert!((c.norm_squared() - expected).abs() < 1.0e-12);
    }

    #[test]
    fn parallel_vectors_have_no_cross_product() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let c = a.cross(a * 2.0);
        assert_eq!(c.norm(), 0.0);
    }

    #[test]
    fn normalisation_gives_unit_length() {
        let unit = Vec3::new(3.0, 0.0, 4.0).normalized().expect("non-zero");
        assert!((unit.norm() - 1.0).abs() < 1.0e-12);
        assert_eq!(unit, Vec3::new(0.6, 0.0, 0.8));
    }

    #[test]
    fn zero_vector_does_not_normalise() {
        assert!(Vec3::ZERO.normalized().is_none());
        assert!(Vec3::new(f64::NAN, 0.0, 0.0).normalized().is_none());
    }

    #[test]
    fn distance_is_the_length_between_points() {
        let a = Vec3::new(1.0, 2.0, 2.0);
        let b = Vec3::new(1.0, 2.0, 5.0);
        assert_eq!(a.distance(b), 3.0);
        assert_eq!(a.distance(a), 0.0);
    }

    #[test]
    fn zero_length_is_recognised() {
        assert!(Vec3::ZERO.is_zero());
        assert!(Vec3::new(1.0e-13, 0.0, 0.0).is_zero());
        assert!(!Vec3::new(1.0e-6, 0.0, 0.0).is_zero());
    }

    #[test]
    fn rotation_about_z_moves_x_to_y() {
        let rotated = Vec3::X
            .rotated(Vec3::Z, std::f64::consts::FRAC_PI_2)
            .expect("unit axis");
        assert!((rotated - Vec3::Y).norm() < 1.0e-12);
        assert!((rotated.norm() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn rotation_about_the_vector_itself_changes_nothing() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let rotated = v.rotated(v, 0.9).expect("non-zero axis");
        assert!((rotated - v).norm() < 1.0e-12);
    }

    #[test]
    fn zero_axis_does_not_rotate() {
        assert!(Vec3::X.rotated(Vec3::ZERO, 1.0).is_none());
    }

    #[test]
    fn index_and_array_conversions_agree() {
        let v = Vec3::new(4.0, 5.0, 6.0);
        assert_eq!(v[0], 4.0);
        assert_eq!(v[2], 6.0);
        let array: [Real; 3] = v.into();
        assert_eq!(Vec3::from(array), v);
        assert_eq!(v.to_array(), array);
    }
}
