//! Four-node quadrilateral shell element (24 DOF).
//!
//! The element combines a bilinear isoparametric membrane field for `u, v`
//! with a Mindlin–Reissner plate field for `w` and the rotations `θx, θy`.
//! Both use 2 × 2 Gauss quadrature; the plate carries the MITC4 assumed
//! transverse-shear field so a thin plate does not lock.
//!
//! The joints are projected onto a best-fit plane, so a mildly warped quad
//! still behaves as a flat one. A joint order that gives a negative Jacobian
//! is rejected instead of producing a matrix with the wrong sign.
//!
//! Local DOF order at each of the four joints is `[u, v, w, θx, θy, θz]`, the
//! same order the frame uses, so shells and patches can share degrees of
//! freedom. The drilling rotation `θz` has no stiffness; the assembler
//! deactivates it as it does for a truss.
//!
//! Sign conventions:
//!
//! * membrane forces `n11`, `n22` and `n12` are positive in tension, per unit
//!   width;
//! * bending moments `m11`, `m22` and `m12` are positive when they put the
//!   positive 3 face in tension;
//! * a pressure acts on the positive 3 face.

use super::{Element, ElementError, InitialStress, ShellForces, DOF_PER_NODE};
use crate::math::{Mat3, Real, Vec3};
use crate::model::Material;

const GAUSS: [Real; 2] = [-0.577_350_269_189_625_7, 0.577_350_269_189_625_7];
const SHELL_DOF: usize = 4 * DOF_PER_NODE;

/// Local `w` degree of freedom at each joint.
const TRANSVERSE_DOF: [usize; 4] = [
    2,
    DOF_PER_NODE + 2,
    2 * DOF_PER_NODE + 2,
    3 * DOF_PER_NODE + 2,
];

/// How close the shell normal may sit to the horizontal before the joint
/// order stops being readable as clockwise or counter-clockwise.
const UPRIGHT_TOLERANCE: Real = 1.0e-6;

/// A four-node shell element carrying membrane stiffness.
#[derive(Clone, Debug)]
pub struct ShellElement {
    /// Area of the projected quadrilateral.
    area: Real,
    thickness: Real,
    axes: Mat3,
    k_local: [Real; SHELL_DOF * SHELL_DOF],
    k_global: [Real; SHELL_DOF * SHELL_DOF],
    /// Joint coordinates in the element plane, centred on the centroid.
    x: [Real; 4],
    y: [Real; 4],
    dm: [[Real; 3]; 3],
    db: [[Real; 3]; 3],
    ds: Real,
    dofs: [usize; SHELL_DOF],
    mass: [Real; SHELL_DOF],
}

impl ShellElement {
    /// Builds a shell from four joint positions.
    ///
    /// The thickness is a property of the element rather than of the joints,
    /// so the same four nodes can carry shells of different thickness.
    pub fn new(
        positions: [Vec3; 4],
        thickness: Real,
        material: &Material,
        node_indices: [usize; 4],
    ) -> Result<Self, ElementError> {
        if !(thickness.is_finite() && thickness > 0.0) {
            return Err(ElementError::InvalidThickness { value: thickness });
        }
        let modulus = material.elastic_modulus;
        if !(modulus.is_finite() && modulus > 0.0) {
            return Err(ElementError::InvalidRigidity {
                name: "E",
                value: modulus,
            });
        }
        let poisson = material.poisson_ratio;
        if !poisson.is_finite() || poisson.abs() >= 1.0 {
            return Err(ElementError::InvalidPoissonRatio { value: poisson });
        }

        let axes = shell_axes(positions)?;
        let centre = (positions[0] + positions[1] + positions[2] + positions[3]) * 0.25;
        let e1 = axes.row(0);
        let e2 = axes.row(1);
        let mut x = [0.0; 4];
        let mut y = [0.0; 4];
        for joint in 0..4 {
            let offset = positions[joint] - centre;
            x[joint] = offset.dot(Vec3::new(e1[0], e1[1], e1[2]));
            y[joint] = offset.dot(Vec3::new(e2[0], e2[1], e2[2]));
        }

        // A counter-clockwise quad has a positive Jacobian at every Gauss
        // point. This rejects a clockwise or crossed joint order.
        let mut projected_area = 0.0;
        for &xi in &GAUSS {
            for &eta in &GAUSS {
                let geometry = Geometry::at(&x, &y, xi, eta);
                if geometry.det < 0.0 {
                    return Err(ElementError::ClockwiseOrdering);
                }
                projected_area += geometry.det;
            }
        }
        if projected_area <= 1.0e-12 {
            return Err(ElementError::DegenerateArea);
        }

        let shear_modulus = modulus / (2.0 * (1.0 + poisson));
        let dm = plane_stress_matrix(modulus * thickness, poisson);
        let db = plane_stress_matrix(modulus * thickness * thickness * thickness / 12.0, poisson);
        let ds = 5.0 / 6.0 * shear_modulus * thickness;
        let mut k_local = [0.0; SHELL_DOF * SHELL_DOF];
        for &xi in &GAUSS {
            for &eta in &GAUSS {
                let geometry = Geometry::at(&x, &y, xi, eta);
                add_btb(&mut k_local, &membrane_b(&geometry), &dm, geometry.det);
                add_btb(&mut k_local, &bending_b(&geometry), &db, geometry.det);
                add_shear(
                    &mut k_local,
                    ds,
                    &shear_b(&x, &y, xi, eta, &geometry),
                    geometry.det,
                );
            }
        }

        let k_global = transform_to_global(&axes, &k_local);

        let mut dofs = [0usize; SHELL_DOF];
        for (joint, &node) in node_indices.iter().enumerate() {
            for component in 0..DOF_PER_NODE {
                dofs[joint * DOF_PER_NODE + component] = node * DOF_PER_NODE + component;
            }
        }

        let total_mass = material.density * thickness * projected_area / 1000.0;
        let rotary = material.density * thickness * thickness * thickness * projected_area
            / 12.0
            / 4.0
            / 1000.0;
        let mut mass = [0.0; SHELL_DOF];
        for joint in 0..4 {
            for component in 0..3 {
                mass[joint * DOF_PER_NODE + component] = total_mass * 0.25;
            }
            for component in 3..6 {
                mass[joint * DOF_PER_NODE + component] = rotary;
            }
        }

        Ok(Self {
            area: projected_area,
            thickness,
            axes,
            k_local,
            k_global,
            x,
            y,
            dm,
            db,
            ds,
            dofs,
            mass,
        })
    }

    /// Area of the projected quadrilateral.
    pub fn area(&self) -> Real {
        self.area
    }

    /// Thickness the element was built with.
    pub fn thickness(&self) -> Real {
        self.thickness
    }

    /// The local matrix, kept for post-processing.
    pub fn local_stiffness(&self) -> &[Real; SHELL_DOF * SHELL_DOF] {
        &self.k_local
    }

    /// Internal force resultants at the natural coordinate `(xi, eta)`.
    pub fn forces_at(&self, d_global: &[Real], xi: Real, eta: Real) -> ShellForces {
        let d = self.to_local(d_global);
        let geometry = Geometry::at(&self.x, &self.y, xi, eta);
        let strain = b_times_d(&membrane_b(&geometry), &d);
        let membrane = mat_vec(&self.dm, &strain);
        let curvature = b_times_d(&bending_b(&geometry), &d);
        let moments = mat_vec(&self.db, &curvature);
        let shear = shear_b(&self.x, &self.y, xi, eta, &geometry);
        ShellForces {
            n11: membrane[0],
            n22: membrane[1],
            n12: membrane[2],
            m11: moments[0],
            m22: moments[1],
            m12: moments[2],
            v13: self.ds * dot_row(&shear[0], &d),
            v23: self.ds * dot_row(&shear[1], &d),
        }
    }

    /// Initial-stress matrix from the membrane resultants [kN/m], tension
    /// positive.
    ///
    /// The transverse slopes are integrated against `nx`, `ny` and `nxy` with
    /// the same 2 × 2 rule as the elastic stiffness, so tension stiffens the
    /// plate and compression softens it. Only `w` carries the term, so a
    /// constant `w` produces no geometric force.
    pub fn geometric_stiffness_matrix(
        &self,
        stress: InitialStress,
    ) -> [Real; SHELL_DOF * SHELL_DOF] {
        let mut local = [0.0; SHELL_DOF * SHELL_DOF];
        if stress.nx == 0.0 && stress.ny == 0.0 && stress.nxy == 0.0 {
            return local;
        }
        for &xi in &GAUSS {
            for &eta in &GAUSS {
                let geometry = Geometry::at(&self.x, &self.y, xi, eta);
                for (row, &row_dof) in TRANSVERSE_DOF.iter().enumerate() {
                    for (column, &column_dof) in TRANSVERSE_DOF.iter().enumerate() {
                        let slope = stress.nx * geometry.dx[row] * geometry.dx[column]
                            + stress.ny * geometry.dy[row] * geometry.dy[column]
                            + stress.nxy
                                * (geometry.dx[row] * geometry.dy[column]
                                    + geometry.dy[row] * geometry.dx[column]);
                        local[row_dof * SHELL_DOF + column_dof] += slope * geometry.det;
                    }
                }
            }
        }
        transform_to_global(&self.axes, &local)
    }

    fn to_global(&self, local: &[Real; SHELL_DOF]) -> [Real; SHELL_DOF] {
        let transpose = self.axes.transpose();
        let mut out = [0.0; SHELL_DOF];
        for block in 0..8 {
            let vector = Vec3::new(local[3 * block], local[3 * block + 1], local[3 * block + 2]);
            let global = transpose.mul_vec(vector);
            out[3 * block] = global.x;
            out[3 * block + 1] = global.y;
            out[3 * block + 2] = global.z;
        }
        out
    }

    fn to_local(&self, global: &[Real]) -> [Real; SHELL_DOF] {
        assert_eq!(
            global.len(),
            SHELL_DOF,
            "a shell has twenty-four global degrees of freedom"
        );
        let mut out = [0.0; SHELL_DOF];
        for block in 0..8 {
            let vector = Vec3::new(
                global[3 * block],
                global[3 * block + 1],
                global[3 * block + 2],
            );
            let local = self.axes.mul_vec(vector);
            out[3 * block] = local.x;
            out[3 * block + 1] = local.y;
            out[3 * block + 2] = local.z;
        }
        out
    }

    /// Work-equivalent nodal loads of a uniform local pressure [kN/m²].
    ///
    /// The bilinear shape functions produce nodal forces only, never nodal
    /// moments.
    fn consistent_load(&self, q: Vec3) -> [Real; SHELL_DOF] {
        let mut weight = [0.0; 4];
        for &xi in &GAUSS {
            for &eta in &GAUSS {
                let geometry = Geometry::at(&self.x, &self.y, xi, eta);
                for (w, n) in weight.iter_mut().zip(geometry.n.iter()) {
                    *w += n * geometry.det;
                }
            }
        }
        let mut local = [0.0; SHELL_DOF];
        for joint in 0..4 {
            local[joint * DOF_PER_NODE] = q.x * weight[joint];
            local[joint * DOF_PER_NODE + 1] = q.y * weight[joint];
            local[joint * DOF_PER_NODE + 2] = q.z * weight[joint];
        }
        self.to_global(&local)
    }
}

impl Element for ShellElement {
    fn dofs(&self) -> &[usize] {
        &self.dofs
    }

    fn stiffness(&self) -> &[Real] {
        &self.k_global
    }

    fn axes(&self) -> &Mat3 {
        &self.axes
    }

    fn equivalent_load(&self, q: Vec3) -> Vec<Real> {
        self.consistent_load(q).to_vec()
    }

    fn end_forces(&self, d: &[Real], q: Vec3) -> Vec<Real> {
        let d_local = self.to_local(d);
        let load = self.consistent_load(q);
        let mut forces = [0.0; SHELL_DOF];
        for row in 0..SHELL_DOF {
            let stiffness_row = &self.k_local[row * SHELL_DOF..(row + 1) * SHELL_DOF];
            let stress: Real = stiffness_row
                .iter()
                .zip(d_local.iter())
                .map(|(k, displacement)| k * displacement)
                .sum();
            forces[row] = stress - load[row];
        }
        forces.to_vec()
    }

    fn length(&self) -> Real {
        self.area.sqrt()
    }

    fn lumped_mass(&self) -> &[Real] {
        &self.mass
    }

    fn geometric_stiffness(&self, stress: InitialStress) -> Option<Vec<Real>> {
        Some(self.geometric_stiffness_matrix(stress).to_vec())
    }
}

/// Local axes of a quadrilateral: 1 along the first usable edge, 3 normal.
///
/// The normal comes from Newell's polygon normal, which averages the corner
/// triangles and so tolerates a mildly warped quad. A listing that runs
/// clockwise about global Z is rejected rather than mirrored, because a
/// mirrored shell takes pressure on the wrong face.
fn shell_axes(positions: [Vec3; 4]) -> Result<Mat3, ElementError> {
    let mut normal = Vec3::ZERO;
    for joint in 0..4 {
        normal = normal + positions[joint].cross(positions[(joint + 1) % 4]);
    }
    let normal = normal.normalized().ok_or(ElementError::DegenerateArea)?;
    if normal.dot(Vec3::Z) < -UPRIGHT_TOLERANCE {
        return Err(ElementError::ClockwiseOrdering);
    }

    let mut e1 = Vec3::ZERO;
    for edge in [
        positions[1] - positions[0],
        positions[2] - positions[1],
        positions[3] - positions[2],
        positions[0] - positions[3],
    ] {
        let projected = edge - normal * edge.dot(normal);
        if projected.norm() > 1.0e-12 {
            e1 = projected.normalized().ok_or(ElementError::DegenerateArea)?;
            break;
        }
    }
    if e1.is_zero() {
        return Err(ElementError::DegenerateArea);
    }
    let e2 = normal.cross(e1);
    Ok(Mat3::from_rows([
        e1.to_array(),
        e2.to_array(),
        normal.to_array(),
    ]))
}

/// Bilinear shape-function values and geometry at one natural point.
struct Geometry {
    n: [Real; 4],
    dxi: [Real; 4],
    deta: [Real; 4],
    dx: [Real; 4],
    dy: [Real; 4],
    x_xi: Real,
    y_xi: Real,
    x_eta: Real,
    y_eta: Real,
    det: Real,
}

impl Geometry {
    fn at(x: &[Real; 4], y: &[Real; 4], xi: Real, eta: Real) -> Self {
        // Joint order: (-1,-1), (+1,-1), (+1,+1), (-1,+1).
        const XI: [Real; 4] = [-1.0, 1.0, 1.0, -1.0];
        const ETA: [Real; 4] = [-1.0, -1.0, 1.0, 1.0];

        let mut n = [0.0; 4];
        let mut dxi = [0.0; 4];
        let mut deta = [0.0; 4];
        for joint in 0..4 {
            n[joint] = 0.25 * (1.0 + xi * XI[joint]) * (1.0 + eta * ETA[joint]);
            dxi[joint] = 0.25 * XI[joint] * (1.0 + eta * ETA[joint]);
            deta[joint] = 0.25 * ETA[joint] * (1.0 + xi * XI[joint]);
        }

        let mut x_xi = 0.0;
        let mut y_xi = 0.0;
        let mut x_eta = 0.0;
        let mut y_eta = 0.0;
        for joint in 0..4 {
            x_xi += dxi[joint] * x[joint];
            y_xi += dxi[joint] * y[joint];
            x_eta += deta[joint] * x[joint];
            y_eta += deta[joint] * y[joint];
        }
        let det = x_xi * y_eta - y_xi * x_eta;
        let inverse = if det.abs() > 1.0e-30 { 1.0 / det } else { 0.0 };

        let mut dx = [0.0; 4];
        let mut dy = [0.0; 4];
        for joint in 0..4 {
            dx[joint] = (y_eta * dxi[joint] - y_xi * deta[joint]) * inverse;
            dy[joint] = (-x_eta * dxi[joint] + x_xi * deta[joint]) * inverse;
        }

        Self {
            n,
            dxi,
            deta,
            dx,
            dy,
            x_xi,
            y_xi,
            x_eta,
            y_eta,
            det,
        }
    }
}

/// Membrane strain–displacement matrix in local components.
fn membrane_b(geometry: &Geometry) -> [[Real; SHELL_DOF]; 3] {
    let mut b = [[0.0; SHELL_DOF]; 3];
    for joint in 0..4 {
        let u = joint * DOF_PER_NODE;
        let v = u + 1;
        b[0][u] = geometry.dx[joint];
        b[1][v] = geometry.dy[joint];
        b[2][u] = geometry.dy[joint];
        b[2][v] = geometry.dx[joint];
    }
    b
}

/// Bending curvature–displacement matrix: `κ1 = -∂θy/∂1`, `κ2 = ∂θx/∂2`,
/// `κ12 = ∂θx/∂1 - ∂θy/∂2`.
fn bending_b(geometry: &Geometry) -> [[Real; SHELL_DOF]; 3] {
    let mut b = [[0.0; SHELL_DOF]; 3];
    for joint in 0..4 {
        let rx = joint * DOF_PER_NODE + 3;
        let ry = joint * DOF_PER_NODE + 4;
        b[0][ry] = -geometry.dx[joint];
        b[1][rx] = geometry.dy[joint];
        b[2][rx] = geometry.dx[joint];
        b[2][ry] = -geometry.dy[joint];
    }
    b
}

/// Covariant transverse shear row `γ_ξζ` (`along_xi`) or `γ_ηζ`.
fn natural_shear_b(geometry: &Geometry, along_xi: bool) -> [Real; SHELL_DOF] {
    let mut row = [0.0; SHELL_DOF];
    for joint in 0..4 {
        let w = joint * DOF_PER_NODE + 2;
        let rx = joint * DOF_PER_NODE + 3;
        let ry = joint * DOF_PER_NODE + 4;
        if along_xi {
            row[w] = geometry.dxi[joint];
            row[ry] += geometry.x_xi * geometry.n[joint];
            row[rx] -= geometry.y_xi * geometry.n[joint];
        } else {
            row[w] = geometry.deta[joint];
            row[ry] += geometry.x_eta * geometry.n[joint];
            row[rx] -= geometry.y_eta * geometry.n[joint];
        }
    }
    row
}

/// MITC4 assumed transverse-shear field in physical local components.
///
/// The covariant shear is sampled at the four edge midpoints and interpolated
/// bilinearly. That is what stops a thin element from over-stiffening.
fn shear_b(
    x: &[Real; 4],
    y: &[Real; 4],
    xi: Real,
    eta: Real,
    geometry: &Geometry,
) -> [[Real; SHELL_DOF]; 2] {
    let bottom = natural_shear_b(&Geometry::at(x, y, 0.0, -1.0), true);
    let top = natural_shear_b(&Geometry::at(x, y, 0.0, 1.0), true);
    let left = natural_shear_b(&Geometry::at(x, y, -1.0, 0.0), false);
    let right = natural_shear_b(&Geometry::at(x, y, 1.0, 0.0), false);

    let inverse = 1.0 / geometry.det;
    let mut b = [[0.0; SHELL_DOF]; 2];
    for dof in 0..SHELL_DOF {
        let natural_xi = 0.5 * (1.0 - eta) * bottom[dof] + 0.5 * (1.0 + eta) * top[dof];
        let natural_eta = 0.5 * (1.0 - xi) * left[dof] + 0.5 * (1.0 + xi) * right[dof];
        b[0][dof] = (geometry.y_eta * natural_xi - geometry.y_xi * natural_eta) * inverse;
        b[1][dof] = (-geometry.x_eta * natural_xi + geometry.x_xi * natural_eta) * inverse;
    }
    b
}

/// Plane-stress constitutive matrix, with the shear term in engineering
/// strain. `rigidity` is `E·t` for membrane and `E·t³/12` for bending.
pub(super) fn plane_stress_matrix(rigidity: Real, poisson: Real) -> [[Real; 3]; 3] {
    let factor = rigidity / (1.0 - poisson * poisson);
    [
        [factor, factor * poisson, 0.0],
        [factor * poisson, factor, 0.0],
        [0.0, 0.0, factor * (1.0 - poisson) / 2.0],
    ]
}

/// Turns a local 24 × 24 matrix into global axes, block by block.
fn transform_to_global(
    axes: &Mat3,
    k_local: &[Real; SHELL_DOF * SHELL_DOF],
) -> [Real; SHELL_DOF * SHELL_DOF] {
    let mut global = [0.0; SHELL_DOF * SHELL_DOF];
    for row_block in 0..8 {
        for column_block in 0..8 {
            let rows = [0, 1, 2].map(|row| {
                [0, 1, 2].map(|column| {
                    k_local[(row_block * 3 + row) * SHELL_DOF + column_block * 3 + column]
                })
            });
            let transformed = axes.congruence(&Mat3::from_rows(rows));
            for row in 0..3 {
                for column in 0..3 {
                    global[(row_block * 3 + row) * SHELL_DOF + column_block * 3 + column] =
                        transformed.row(row)[column];
                }
            }
        }
    }
    global
}

fn b_times_d(b: &[[Real; SHELL_DOF]; 3], d: &[Real; SHELL_DOF]) -> [Real; 3] {
    let mut out = [0.0; 3];
    for row in 0..3 {
        for column in 0..SHELL_DOF {
            out[row] += b[row][column] * d[column];
        }
    }
    out
}

fn mat_vec(m: &[[Real; 3]; 3], v: &[Real; 3]) -> [Real; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

fn dot_row(row: &[Real; SHELL_DOF], d: &[Real; SHELL_DOF]) -> Real {
    row.iter().zip(d.iter()).map(|(a, b)| a * b).sum()
}

fn add_shear(
    k: &mut [Real; SHELL_DOF * SHELL_DOF],
    ds: Real,
    b: &[[Real; SHELL_DOF]; 2],
    weight: Real,
) {
    for row in 0..SHELL_DOF {
        for column in 0..SHELL_DOF {
            let value = ds * (b[0][row] * b[0][column] + b[1][row] * b[1][column]);
            k[row * SHELL_DOF + column] += value * weight;
        }
    }
}

fn add_btb(
    k: &mut [Real; SHELL_DOF * SHELL_DOF],
    b: &[[Real; SHELL_DOF]; 3],
    d: &[[Real; 3]; 3],
    weight: Real,
) {
    for row in 0..SHELL_DOF {
        for column in 0..SHELL_DOF {
            let mut value = 0.0;
            for a in 0..3 {
                for c in 0..3 {
                    value += b[a][row] * d[a][c] * b[c][column];
                }
            }
            k[row * SHELL_DOF + column] += value * weight;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Material, MaterialId};

    const MODULUS: Real = 210_000_000.0;
    const POISSON: Real = 0.3;
    const THICKNESS: Real = 0.1;

    fn material() -> Material {
        Material::new(MaterialId::new(1), "S355", 7850.0, MODULUS, POISSON)
    }

    /// A rectangle in the global xy plane, counter-clockwise from the origin.
    fn rectangle(width: Real, height: Real) -> [Vec3; 4] {
        [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(width, 0.0, 0.0),
            Vec3::new(width, height, 0.0),
            Vec3::new(0.0, height, 0.0),
        ]
    }

    fn shell(positions: [Vec3; 4]) -> ShellElement {
        ShellElement::new(positions, THICKNESS, &material(), [0, 1, 2, 3])
            .expect("a counter-clockwise rectangle with a positive thickness")
    }

    /// Displacements that put the membrane in uniaxial stress.
    ///
    /// The transverse strain is `-nu` times the axial one, so the transverse
    /// stress cancels and the state is pure uniaxial.
    fn uniaxial_field(positions: &[Vec3; 4], strain: Real) -> [Real; SHELL_DOF] {
        let mut d = [0.0; SHELL_DOF];
        for joint in 0..4 {
            d[joint * DOF_PER_NODE] = strain * positions[joint].x;
            d[joint * DOF_PER_NODE + 1] = -POISSON * strain * positions[joint].y;
        }
        d
    }

    fn largest_stiffness(element: &ShellElement) -> Real {
        element
            .local_stiffness()
            .iter()
            .fold(0.0_f64, |largest, value| largest.max(value.abs()))
    }

    #[test]
    fn a_constant_strain_field_gives_constant_stress() {
        let width = 3.0;
        let height = 2.0;
        let strain = 1.0e-4;
        let positions = rectangle(width, height);
        let element = shell(positions);
        let d = uniaxial_field(&positions, strain);

        let resultant = MODULUS * THICKNESS * strain;
        let centre = element.forces_at(&d, 0.0, 0.0);
        let corner = element.forces_at(&d, 0.9, -0.8);
        assert!(
            (centre.n11 - resultant).abs() < 1.0e-9 * resultant,
            "{centre:?}"
        );
        assert!(
            (corner.n11 - resultant).abs() < 1.0e-9 * resultant,
            "{corner:?}"
        );
        assert!(centre.n22.abs() < 1.0e-9 * resultant);
        assert!(centre.n12.abs() < 1.0e-9 * resultant);
    }

    #[test]
    fn the_nodal_forces_match_the_edge_traction() {
        let width = 3.0;
        let height = 2.0;
        let strain = 1.0e-4;
        let positions = rectangle(width, height);
        let element = shell(positions);
        let d = uniaxial_field(&positions, strain);
        let forces = element.end_forces(&d, Vec3::ZERO);

        let resultant = MODULUS * THICKNESS * strain * height / 2.0;
        assert!((forces[DOF_PER_NODE] - resultant).abs() < 1.0e-9 * resultant);
        assert!((forces[2 * DOF_PER_NODE] - resultant).abs() < 1.0e-9 * resultant);
        assert!((forces[0] + resultant).abs() < 1.0e-9 * resultant);
        assert!((forces[3 * DOF_PER_NODE] + resultant).abs() < 1.0e-9 * resultant);
        assert!(forces[1].abs() < 1.0e-9 * resultant);
    }

    #[test]
    fn a_rigid_translation_carries_no_force() {
        let element = shell(rectangle(3.0, 2.0));
        let mut d = [0.0; SHELL_DOF];
        for joint in 0..4 {
            d[joint * DOF_PER_NODE] = 1.5;
            d[joint * DOF_PER_NODE + 1] = -0.5;
        }
        let forces = element.end_forces(&d, Vec3::ZERO);
        let allowed = 1.0e-9 * largest_stiffness(&element);
        for force in forces {
            assert!(force.abs() < allowed, "rigid force {force}");
        }
    }

    #[test]
    fn a_clockwise_joint_order_is_rejected() {
        let positions = rectangle(3.0, 2.0);
        let clockwise = [positions[0], positions[3], positions[2], positions[1]];
        assert!(matches!(
            ShellElement::new(clockwise, THICKNESS, &material(), [0, 1, 2, 3]),
            Err(ElementError::ClockwiseOrdering)
        ));
    }

    #[test]
    fn a_degenerate_area_is_rejected() {
        let collapsed = [Vec3::ZERO, Vec3::ZERO, Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0)];
        assert!(matches!(
            ShellElement::new(collapsed, THICKNESS, &material(), [0, 1, 2, 3]),
            Err(ElementError::DegenerateArea)
        ));
    }

    #[test]
    fn bad_thickness_and_material_are_rejected() {
        let positions = rectangle(3.0, 2.0);
        assert!(matches!(
            ShellElement::new(positions, 0.0, &material(), [0, 1, 2, 3]),
            Err(ElementError::InvalidThickness { value: 0.0 })
        ));
        let mut stiff = material();
        stiff.poisson_ratio = 1.5;
        assert!(matches!(
            ShellElement::new(positions, THICKNESS, &stiff, [0, 1, 2, 3]),
            Err(ElementError::InvalidPoissonRatio { .. })
        ));
    }

    #[test]
    fn the_shell_uses_the_frame_dof_order() {
        let element = shell(rectangle(3.0, 2.0));
        let formulation: &dyn Element = &element;
        assert_eq!(formulation.dofs().len(), SHELL_DOF);
        assert_eq!(formulation.dofs()[0], 0);
        assert_eq!(formulation.dofs()[2 * DOF_PER_NODE], 2 * DOF_PER_NODE);
        assert_eq!(formulation.stiffness().len(), SHELL_DOF * SHELL_DOF);
        assert!((element.area() - 6.0).abs() < 1.0e-12);
        assert_eq!(element.thickness(), THICKNESS);
        // A shell reads the membrane resultants, so an axial-only state
        // contributes nothing.
        let axial_only = formulation
            .geometric_stiffness(InitialStress::axial(-1.0))
            .expect("a shell has an initial-stress matrix");
        assert!(axial_only.iter().all(|&value| value == 0.0));
    }

    #[test]
    fn the_stiffness_is_symmetric() {
        let element = shell(rectangle(3.0, 2.0));
        let stiffness = element.local_stiffness();
        let allowed = 1.0e-12 * largest_stiffness(&element);
        for row in 0..SHELL_DOF {
            for column in 0..SHELL_DOF {
                let a = stiffness[row * SHELL_DOF + column];
                let b = stiffness[column * SHELL_DOF + row];
                assert!((a - b).abs() <= allowed, "entries {a} and {b} differ");
            }
        }
    }

    fn plate_rigidity() -> Real {
        MODULUS * THICKNESS * THICKNESS * THICKNESS / (12.0 * (1.0 - POISSON * POISSON))
    }

    /// A constant-curvature field with Kirchhoff rotations, in local
    /// coordinates.
    fn curvature_field(positions: &[Vec3; 4], k11: Real, k22: Real) -> [Real; SHELL_DOF] {
        let centre = (positions[0] + positions[1] + positions[2] + positions[3]) * 0.25;
        let mut d = [0.0; SHELL_DOF];
        for joint in 0..4 {
            let x = positions[joint].x - centre.x;
            let y = positions[joint].y - centre.y;
            d[joint * DOF_PER_NODE + 2] = 0.5 * (k11 * x * x + k22 * y * y);
            d[joint * DOF_PER_NODE + 3] = k22 * y;
            d[joint * DOF_PER_NODE + 4] = -k11 * x;
        }
        d
    }

    #[test]
    fn a_constant_curvature_field_gives_the_plate_moments() {
        let positions = rectangle(3.0, 2.0);
        let element = shell(positions);
        let k11 = 1.0e-3;
        let k22 = -4.0e-4;
        let d = curvature_field(&positions, k11, k22);

        let rigidity = plate_rigidity();
        let allowed = 1.0e-9 * rigidity * (k11.abs() + k22.abs());
        let forces = element.forces_at(&d, 0.0, 0.0);
        let m11 = rigidity * (k11 + POISSON * k22);
        let m22 = rigidity * (k22 + POISSON * k11);
        assert!((forces.m11 - m11).abs() < allowed, "{forces:?}");
        assert!((forces.m22 - m22).abs() < allowed, "{forces:?}");
        assert!(forces.m12.abs() < allowed, "{forces:?}");
    }

    #[test]
    fn a_kirchhoff_field_carries_no_transverse_shear() {
        let positions = rectangle(3.0, 2.0);
        let element = shell(positions);
        let k11 = 1.0e-3;
        let k22 = -4.0e-4;
        let d = curvature_field(&positions, k11, k22);

        let allowed = 1.0e-9 * element.ds * (k11.abs() + k22.abs());
        let forces = element.forces_at(&d, 0.4, -0.3);
        assert!(forces.v13.abs() < allowed, "{forces:?}");
        assert!(forces.v23.abs() < allowed, "{forces:?}");
    }

    #[test]
    fn a_pressure_reaches_the_transverse_degrees_of_freedom() {
        let element = shell(rectangle(3.0, 2.0));
        let load = element.equivalent_load(Vec3::new(0.0, 0.0, -2.0));
        let total: Real = (0..4).map(|joint| load[joint * DOF_PER_NODE + 2]).sum();
        assert!((total + 12.0).abs() < 1.0e-9, "{total}");
        for joint in 0..4 {
            assert!(load[joint * DOF_PER_NODE].abs() < 1.0e-12);
            assert!(load[joint * DOF_PER_NODE + 1].abs() < 1.0e-12);
        }
    }

    /// A saddle-shaped `w` field, so the field is neither constant nor a rigid
    /// rotation.
    fn warped_field(positions: &[Vec3; 4], amplitude: Real) -> [Real; SHELL_DOF] {
        let centre = (positions[0] + positions[1] + positions[2] + positions[3]) * 0.25;
        let mut d = [0.0; SHELL_DOF];
        for joint in 0..4 {
            let x = positions[joint].x - centre.x;
            let y = positions[joint].y - centre.y;
            d[joint * DOF_PER_NODE + 2] = amplitude * x * y;
        }
        d
    }

    fn quadratic_form(matrix: &[Real], d: &[Real; SHELL_DOF]) -> Real {
        let mut sum = 0.0;
        for row in 0..SHELL_DOF {
            for column in 0..SHELL_DOF {
                sum += d[row] * matrix[row * SHELL_DOF + column] * d[column];
            }
        }
        sum
    }

    #[test]
    fn tension_stiffens_and_compression_softens_the_plate() {
        let positions = rectangle(3.0, 2.0);
        let element = shell(positions);
        let d = warped_field(&positions, 1.0e-3);

        let tension = element.geometric_stiffness_matrix(InitialStress::membrane(10.0, 10.0, 0.0));
        let compression =
            element.geometric_stiffness_matrix(InitialStress::membrane(-10.0, -10.0, 0.0));
        let tension_energy = quadratic_form(&tension, &d);
        let compression_energy = quadratic_form(&compression, &d);
        assert!(tension_energy > 0.0, "{tension_energy}");
        assert!(compression_energy < 0.0, "{compression_energy}");
        assert!(
            (tension_energy + compression_energy).abs() <= 1.0e-9 * tension_energy,
            "the sign must flip with the stress"
        );

        let zero = element.geometric_stiffness_matrix(InitialStress::default());
        assert_eq!(quadratic_form(&zero, &d), 0.0);

        let doubled = element.geometric_stiffness_matrix(InitialStress::membrane(20.0, 20.0, 0.0));
        assert!(
            (quadratic_form(&doubled, &d) - 2.0 * tension_energy).abs() <= 1.0e-9 * tension_energy,
            "the matrix must be linear in the resultants"
        );
    }

    #[test]
    fn a_rigid_translation_carries_no_geometric_force() {
        let element = shell(rectangle(3.0, 2.0));
        let stress = InitialStress::membrane(10.0, -4.0, 3.0);
        let matrix = element.geometric_stiffness_matrix(stress);
        let allowed = 1.0e-9
            * matrix
                .iter()
                .fold(0.0_f64, |peak, value| peak.max(value.abs()))
            * SHELL_DOF as Real;

        let mut translation = [0.0; SHELL_DOF];
        for joint in 0..4 {
            translation[joint * DOF_PER_NODE + 2] = 1.0;
        }
        for row in 0..SHELL_DOF {
            let force: Real = (0..SHELL_DOF)
                .map(|column| matrix[row * SHELL_DOF + column] * translation[column])
                .sum();
            assert!(force.abs() < allowed, "dof {row} force {force}");
        }
    }

    #[test]
    fn the_trait_reports_the_same_geometric_stiffness() {
        let element = shell(rectangle(3.0, 2.0));
        let stress = InitialStress::membrane(5.0, -2.0, 1.5);
        let formulation: &dyn Element = &element;
        let through_trait = formulation
            .geometric_stiffness(stress)
            .expect("a shell has an initial-stress matrix");
        let direct = element.geometric_stiffness_matrix(stress);
        assert_eq!(through_trait.len(), SHELL_DOF * SHELL_DOF);
        assert_eq!(through_trait.as_slice(), direct.as_slice());
    }
}
