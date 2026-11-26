//! Isogeometric patch element, six degrees of freedom per control point.
//!
//! The patch's own spline basis supplies the membrane, bending and transverse
//! shear fields, so the geometry and the displacement use the same functions.
//! Membrane and bending are integrated with the full Gauss rule for the spline
//! degree; shear with one order less. Under-integrating only the shear term is
//! what stops a thin patch from locking without introducing an assumed shear
//! field.
//!
//! The patch is treated as flat: control points are projected onto the plane
//! through the patch centre spanned by the parametric tangents. A mildly
//! curved patch therefore behaves as its flat projection, the same
//! approximation the four-node shell makes for a warped quad.
//!
//! Local DOF order at each control point is `[u, v, w, θx, θy, θz]`, the order
//! shell4 uses, so ties, loads and result post-processing can be shared. The
//! drilling rotation `θz` carries no stiffness.
//!
//! Sign conventions match shell4: membrane resultants are positive in tension,
//! bending moments are positive when they put the positive 3 face in tension,
//! and a pressure acts on the positive 3 face.

use super::shell4::plane_stress_matrix;
use super::{Element, ElementError, InitialStress, ShellForces, DOF_PER_NODE};
use crate::math::{Mat3, Real, Vec3};
use crate::model::Material;
use crate::spline::BSplinePatch;

/// A tensor-product patch carrying membrane, bending and shear stiffness.
#[derive(Clone, Debug)]
pub struct IgaPatchElement {
    patch: BSplinePatch,
    thickness: Real,
    axes: Mat3,
    area: Real,
    k_local: Vec<Real>,
    k_global: Vec<Real>,
    /// Control-point coordinates in the element plane, centred on the patch.
    x: Vec<Real>,
    y: Vec<Real>,
    dm: [[Real; 3]; 3],
    db: [[Real; 3]; 3],
    ds: Real,
    dofs: Vec<usize>,
    mass: Vec<Real>,
}

impl IgaPatchElement {
    /// Builds a patch element from its spline patch.
    ///
    /// `node_indices` gives the model node of every control point in the
    /// patch's row-major order, the same order [`BSplinePatch::basis`]
    /// returns.
    pub fn new(
        patch: BSplinePatch,
        thickness: Real,
        material: &Material,
        node_indices: &[usize],
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
        let n_cp = patch.n_cp();
        if node_indices.len() != n_cp {
            return Err(ElementError::NodeCount {
                expected: n_cp,
                count: node_indices.len(),
            });
        }

        let axes = patch_axes(&patch)?;
        let ((u0, u1), (v0, v1)) = patch.domain();
        let centre = patch.point(0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let e1 = Vec3::from_array(axes.row(0));
        let e2 = Vec3::from_array(axes.row(1));
        let mut x = Vec::with_capacity(n_cp);
        let mut y = Vec::with_capacity(n_cp);
        for point in patch.control_net().iter().flatten() {
            let offset = *point - centre;
            x.push(offset.dot(e1));
            y.push(offset.dot(e2));
        }

        let shear_modulus = modulus / (2.0 * (1.0 + poisson));
        let dm = plane_stress_matrix(modulus * thickness, poisson);
        let db = plane_stress_matrix(modulus * thickness * thickness * thickness / 12.0, poisson);
        let ds = 5.0 / 6.0 * shear_modulus * thickness;

        let size = n_cp * DOF_PER_NODE;
        let mut k_local = vec![0.0; size * size];
        let (nodes_u, weights_u) = gauss_rule(patch.space_u().degree() + 1);
        let (nodes_v, weights_v) = gauss_rule(patch.space_v().degree() + 1);
        // The rule is written on [-1, 1]; the patch parameter runs over
        // [u0, u1] × [v0, v1], so each node and weight is scaled onto it.
        let scale = 0.5 * (u1 - u0) * 0.5 * (v1 - v0);
        let mut area = 0.0;
        // Lumped mass share of each control point: `∫ N_A dA`. The shares sum
        // to the area because the basis partitions to one. They are equal only
        // when every basis function integrates alike, which is the case over a
        // single span but not across an interior knot.
        let mut lumped = vec![0.0; n_cp];
        for (iu, &node_u) in nodes_u.iter().enumerate() {
            for (iv, &node_v) in nodes_v.iter().enumerate() {
                let u = 0.5 * ((u1 - u0) * node_u + u0 + u1);
                let v = 0.5 * ((v1 - v0) * node_v + v0 + v1);
                let mapping = mapping(&patch, &x, &y, u, v);
                if !mapping.det.is_finite() || mapping.det <= 0.0 {
                    return Err(ElementError::DegenerateArea);
                }
                let weight = weights_u[iu] * weights_v[iv] * mapping.det * scale;
                area += weight;
                for (control_point, &value) in mapping.value.iter().enumerate() {
                    lumped[control_point] += value * weight;
                }
                add_btb(&mut k_local, &membrane_b(&mapping), &dm, weight, size);
                add_btb(&mut k_local, &bending_b(&mapping), &db, weight, size);
            }
        }

        // One order less for shear: the transverse shear field is lower degree
        // than the bending field, and integrating it fully over-stiffens a
        // thin patch.
        let (shear_u, shear_weights_u) = gauss_rule(patch.space_u().degree().max(1));
        let (shear_v, shear_weights_v) = gauss_rule(patch.space_v().degree().max(1));
        for (iu, &node_u) in shear_u.iter().enumerate() {
            for (iv, &node_v) in shear_v.iter().enumerate() {
                let u = 0.5 * ((u1 - u0) * node_u + u0 + u1);
                let v = 0.5 * ((v1 - v0) * node_v + v0 + v1);
                let mapping = mapping(&patch, &x, &y, u, v);
                let weight = shear_weights_u[iu] * shear_weights_v[iv] * mapping.det * scale;
                add_shear(&mut k_local, &shear_b(&mapping), ds, weight, size);
            }
        }

        let k_global = transform_to_global(&axes, &k_local, size);

        let mut dofs = Vec::with_capacity(size);
        for &node in node_indices {
            for component in 0..DOF_PER_NODE {
                dofs.push(node * DOF_PER_NODE + component);
            }
        }

        let translational = material.density * thickness / 1000.0;
        let rotary = material.density * thickness * thickness * thickness / 12.0 / 1000.0;
        let mut mass = vec![0.0; size];
        for (control_point, &share) in lumped.iter().enumerate() {
            for component in 0..3 {
                mass[control_point * DOF_PER_NODE + component] = translational * share;
            }
            for component in 3..DOF_PER_NODE {
                mass[control_point * DOF_PER_NODE + component] = rotary * share;
            }
        }

        Ok(Self {
            patch,
            thickness,
            axes,
            area,
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

    /// Area of the flat patch.
    pub fn area(&self) -> Real {
        self.area
    }

    /// Thickness the element was built with.
    pub fn thickness(&self) -> Real {
        self.thickness
    }

    /// The local matrix, kept for post-processing.
    pub fn local_stiffness(&self) -> &[Real] {
        &self.k_local
    }

    /// Internal force resultants at the parametric point `(u, v)`.
    pub fn forces_at(&self, d_global: &[Real], u: Real, v: Real) -> ShellForces {
        let d = self.to_local(d_global);
        let mapping = mapping(&self.patch, &self.x, &self.y, u, v);
        let strain = b_times_d(&membrane_b(&mapping), &d);
        let membrane = mat_vec(&self.dm, &strain);
        let curvature = b_times_d(&bending_b(&mapping), &d);
        let moments = mat_vec(&self.db, &curvature);
        let shear = shear_b(&mapping);
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

    fn to_global(&self, local: &[Real]) -> Vec<Real> {
        let transpose = self.axes.transpose();
        let mut out = vec![0.0; local.len()];
        for block in 0..local.len() / 3 {
            let vector = Vec3::new(local[3 * block], local[3 * block + 1], local[3 * block + 2]);
            let global = transpose.mul_vec(vector);
            out[3 * block] = global.x;
            out[3 * block + 1] = global.y;
            out[3 * block + 2] = global.z;
        }
        out
    }

    fn to_local(&self, global: &[Real]) -> Vec<Real> {
        assert_eq!(
            global.len(),
            self.dofs.len(),
            "a patch has six degrees of freedom per control point"
        );
        let mut out = vec![0.0; global.len()];
        for block in 0..global.len() / 3 {
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

    /// Work-equivalent control-point loads of a uniform local pressure.
    fn consistent_load(&self, q: Vec3) -> Vec<Real> {
        let mut local = vec![0.0; self.dofs.len()];
        let ((u0, u1), (v0, v1)) = self.patch.domain();
        let scale = 0.5 * (u1 - u0) * 0.5 * (v1 - v0);
        let (nodes_u, weights_u) = gauss_rule(self.patch.space_u().degree() + 1);
        let (nodes_v, weights_v) = gauss_rule(self.patch.space_v().degree() + 1);
        for (iu, &node_u) in nodes_u.iter().enumerate() {
            for (iv, &node_v) in nodes_v.iter().enumerate() {
                let u = 0.5 * ((u1 - u0) * node_u + u0 + u1);
                let v = 0.5 * ((v1 - v0) * node_v + v0 + v1);
                let mapping = mapping(&self.patch, &self.x, &self.y, u, v);
                let weight = weights_u[iu] * weights_v[iv] * mapping.det * scale;
                for (control_point, &value) in mapping.value.iter().enumerate() {
                    let dof = control_point * DOF_PER_NODE;
                    local[dof] += q.x * value * weight;
                    local[dof + 1] += q.y * value * weight;
                    local[dof + 2] += q.z * value * weight;
                }
            }
        }
        self.to_global(&local)
    }
}

impl Element for IgaPatchElement {
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
        self.consistent_load(q)
    }

    fn end_forces(&self, d: &[Real], q: Vec3) -> Vec<Real> {
        let d_local = self.to_local(d);
        let load = self.consistent_load(q);
        let size = self.dofs.len();
        let mut forces = vec![0.0; size];
        for row in 0..size {
            let stress: Real = (0..size)
                .map(|column| self.k_local[row * size + column] * d_local[column])
                .sum();
            forces[row] = stress - load[row];
        }
        forces
    }

    fn length(&self) -> Real {
        self.area.sqrt()
    }

    fn lumped_mass(&self) -> &[Real] {
        &self.mass
    }

    fn geometric_stiffness(&self, _stress: InitialStress) -> Option<Vec<Real>> {
        None
    }
}

/// Element axes: 1 along the parametric tangent, 3 along the surface normal.
fn patch_axes(patch: &BSplinePatch) -> Result<Mat3, ElementError> {
    let ((u0, u1), (v0, v1)) = patch.domain();
    let (du, dv) = patch.tangents(0.5 * (u0 + u1), 0.5 * (v0 + v1));
    let normal = du
        .cross(dv)
        .normalized()
        .ok_or(ElementError::DegenerateArea)?;
    let projected = du - normal * du.dot(normal);
    let e1 = projected.normalized().ok_or(ElementError::DegenerateArea)?;
    let e2 = normal.cross(e1);
    Ok(Mat3::from_rows([
        e1.to_array(),
        e2.to_array(),
        normal.to_array(),
    ]))
}

/// Basis values and physical gradients at one parametric point.
struct Mapping {
    value: Vec<Real>,
    dx: Vec<Real>,
    dy: Vec<Real>,
    det: Real,
}

fn mapping(patch: &BSplinePatch, x: &[Real], y: &[Real], u: Real, v: Real) -> Mapping {
    let basis = patch.basis(u, v);
    let mut x_xi = 0.0;
    let mut y_xi = 0.0;
    let mut x_eta = 0.0;
    let mut y_eta = 0.0;
    for (((&du, &dv), &x), &y) in basis.du.iter().zip(&basis.dv).zip(x).zip(y) {
        x_xi += du * x;
        y_xi += du * y;
        x_eta += dv * x;
        y_eta += dv * y;
    }
    let det = x_xi * y_eta - y_xi * x_eta;
    let inverse = if det.abs() > 1.0e-30 { 1.0 / det } else { 0.0 };
    let dx = basis
        .du
        .iter()
        .zip(&basis.dv)
        .map(|(&du, &dv)| (y_eta * du - y_xi * dv) * inverse)
        .collect();
    let dy = basis
        .du
        .iter()
        .zip(&basis.dv)
        .map(|(&du, &dv)| (-x_eta * du + x_xi * dv) * inverse)
        .collect();
    Mapping {
        value: basis.values,
        dx,
        dy,
        det,
    }
}

/// Membrane strain–displacement matrix in local components.
fn membrane_b(mapping: &Mapping) -> [Vec<Real>; 3] {
    let size = mapping.dx.len() * DOF_PER_NODE;
    let mut b = [vec![0.0; size], vec![0.0; size], vec![0.0; size]];
    for (control_point, (&dx, &dy)) in mapping.dx.iter().zip(&mapping.dy).enumerate() {
        let u = control_point * DOF_PER_NODE;
        let v = u + 1;
        b[0][u] = dx;
        b[1][v] = dy;
        b[2][u] = dy;
        b[2][v] = dx;
    }
    b
}

/// Bending curvature–displacement matrix: `κ1 = -∂θy/∂1`, `κ2 = ∂θx/∂2`,
/// `κ12 = ∂θx/∂1 - ∂θy/∂2`.
fn bending_b(mapping: &Mapping) -> [Vec<Real>; 3] {
    let size = mapping.dx.len() * DOF_PER_NODE;
    let mut b = [vec![0.0; size], vec![0.0; size], vec![0.0; size]];
    for (control_point, (&dx, &dy)) in mapping.dx.iter().zip(&mapping.dy).enumerate() {
        let rx = control_point * DOF_PER_NODE + 3;
        let ry = control_point * DOF_PER_NODE + 4;
        b[0][ry] = -dx;
        b[1][rx] = dy;
        b[2][rx] = dx;
        b[2][ry] = -dy;
    }
    b
}

/// Transverse shear rows `γ13 = ∂w/∂1 + θy` and `γ23 = ∂w/∂2 - θx`.
fn shear_b(mapping: &Mapping) -> [Vec<Real>; 2] {
    let size = mapping.dx.len() * DOF_PER_NODE;
    let mut b = [vec![0.0; size], vec![0.0; size]];
    for (control_point, ((&dx, &dy), &value)) in mapping
        .dx
        .iter()
        .zip(&mapping.dy)
        .zip(&mapping.value)
        .enumerate()
    {
        let w = control_point * DOF_PER_NODE + 2;
        let rx = control_point * DOF_PER_NODE + 3;
        let ry = control_point * DOF_PER_NODE + 4;
        b[0][w] = dx;
        b[0][ry] += value;
        b[1][w] = dy;
        b[1][rx] -= value;
    }
    b
}

fn add_btb(k: &mut [Real], b: &[Vec<Real>; 3], d: &[[Real; 3]; 3], weight: Real, size: usize) {
    for row in 0..size {
        for column in 0..size {
            let mut value = 0.0;
            for a in 0..3 {
                for c in 0..3 {
                    value += b[a][row] * d[a][c] * b[c][column];
                }
            }
            k[row * size + column] += value * weight;
        }
    }
}

fn add_shear(k: &mut [Real], b: &[Vec<Real>; 2], ds: Real, weight: Real, size: usize) {
    for row in 0..size {
        for column in 0..size {
            let value = ds * (b[0][row] * b[0][column] + b[1][row] * b[1][column]);
            k[row * size + column] += value * weight;
        }
    }
}

/// Turns a local matrix into global axes, block by block.
fn transform_to_global(axes: &Mat3, local: &[Real], size: usize) -> Vec<Real> {
    let mut global = vec![0.0; size * size];
    for row_block in 0..size / 3 {
        for column_block in 0..size / 3 {
            let rows = [0, 1, 2].map(|row| {
                [0, 1, 2]
                    .map(|column| local[(row_block * 3 + row) * size + column_block * 3 + column])
            });
            let transformed = axes.congruence(&Mat3::from_rows(rows));
            for row in 0..3 {
                for column in 0..3 {
                    global[(row_block * 3 + row) * size + column_block * 3 + column] =
                        transformed.row(row)[column];
                }
            }
        }
    }
    global
}

fn b_times_d(b: &[Vec<Real>; 3], d: &[Real]) -> [Real; 3] {
    let mut out = [0.0; 3];
    for (row, values) in b.iter().enumerate() {
        out[row] = values
            .iter()
            .zip(d)
            .map(|(value, displacement)| value * displacement)
            .sum();
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

fn dot_row(row: &[Real], d: &[Real]) -> Real {
    row.iter().zip(d).map(|(a, b)| a * b).sum()
}

/// Gauss–Legendre nodes and weights on `[-1, 1]`.
///
/// The spline degree sets the rule at run time, so the nodes are computed
/// rather than tabulated. Newton's method converges in a handful of steps for
/// the degrees a patch uses.
fn gauss_rule(count: usize) -> (Vec<Real>, Vec<Real>) {
    let count = count.max(1);
    let mut nodes = vec![0.0; count];
    let mut weights = vec![0.0; count];
    let half = count.div_ceil(2);
    for index in 0..half {
        let mut x = (std::f64::consts::PI * (index as Real + 0.75) / (count as Real + 0.5)).cos();
        for _ in 0..20 {
            let (polynomial, derivative) = legendre(count, x);
            let step = polynomial / derivative;
            x -= step;
            if step.abs() < 1.0e-15 {
                break;
            }
        }
        let (_, derivative) = legendre(count, x);
        let weight = 2.0 / ((1.0 - x * x) * derivative * derivative);
        nodes[index] = -x;
        weights[index] = weight;
        nodes[count - 1 - index] = x;
        weights[count - 1 - index] = weight;
    }
    (nodes, weights)
}

/// Legendre polynomial and its derivative at `x`.
fn legendre(n: usize, x: Real) -> (Real, Real) {
    let mut previous = 1.0;
    let mut current = x;
    for order in 2..=n {
        let next = (((2 * order - 1) as Real * x * current) - ((order - 1) as Real * previous))
            / order as Real;
        previous = current;
        current = next;
    }
    if n == 0 {
        return (1.0, 0.0);
    }
    let derivative = (n as Real) * (x * current - previous) / (x * x - 1.0);
    (current, derivative)
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

    /// A rectangular control net on the global xy plane. For degree two the
    /// evenly spaced net reproduces the rectangle exactly.
    fn rectangle(width: Real, height: Real, degree_u: usize, degree_v: usize) -> BSplinePatch {
        let rows = degree_u + 1;
        let columns = degree_v + 1;
        let mut net = Vec::with_capacity(rows);
        for row in 0..rows {
            let mut points = Vec::with_capacity(columns);
            for column in 0..columns {
                let x = width * row as Real / (rows - 1) as Real;
                let y = height * column as Real / (columns - 1) as Real;
                points.push(Vec3::new(x, y, 0.0));
            }
            net.push(points);
        }
        BSplinePatch::from_grid(&net, degree_u, degree_v).expect("a rectangular net")
    }

    /// A quadratic strip with one interior knot in `u`: four control points
    /// placed at the Greville abscissae so the map stays linear across the
    /// knot. The end basis functions integrate to a sixth of the line, the
    /// interior ones to a third.
    fn quadratic_strip(width: Real, height: Real) -> BSplinePatch {
        let mut net = Vec::with_capacity(4);
        for fraction in [0.0, 0.25, 0.75, 1.0] {
            let x = width * fraction;
            net.push(vec![Vec3::new(x, 0.0, 0.0), Vec3::new(x, height, 0.0)]);
        }
        BSplinePatch::from_grid(&net, 2, 1).expect("a quadratic strip")
    }

    /// A rectangular net rotated about global X, so the patch normal leaves
    /// the xy plane.
    fn rotated_rectangle(width: Real, height: Real, angle: Real) -> BSplinePatch {
        let (sin, cos) = angle.sin_cos();
        let point = |x: Real, y: Real| Vec3::new(x, y * cos, y * sin);
        let net = vec![
            vec![point(0.0, 0.0), point(0.0, height)],
            vec![point(width, 0.0), point(width, height)],
        ];
        BSplinePatch::from_grid(&net, 1, 1).expect("a rotated rectangular net")
    }

    fn element(patch: BSplinePatch) -> IgaPatchElement {
        let nodes: Vec<usize> = (0..patch.n_cp()).collect();
        IgaPatchElement::new(patch, THICKNESS, &material(), &nodes)
            .expect("a flat patch with a positive thickness")
    }

    fn control_points(element: &IgaPatchElement) -> Vec<Vec3> {
        element
            .patch
            .control_net()
            .iter()
            .flatten()
            .copied()
            .collect()
    }

    fn field(
        element: &IgaPatchElement,
        mut value: impl FnMut(Vec3) -> [Real; DOF_PER_NODE],
    ) -> Vec<Real> {
        let mut d = Vec::new();
        for point in control_points(element) {
            d.extend_from_slice(&value(point));
        }
        d
    }

    fn largest_stiffness(element: &IgaPatchElement) -> Real {
        element
            .local_stiffness()
            .iter()
            .fold(0.0_f64, |largest, value| largest.max(value.abs()))
    }

    #[test]
    fn the_gauss_rule_integrates_polynomials_exactly() {
        let (nodes, weights) = gauss_rule(3);
        let quadratic: Real = weights
            .iter()
            .zip(&nodes)
            .map(|(weight, node)| weight * node * node)
            .sum();
        let quartic: Real = weights
            .iter()
            .zip(&nodes)
            .map(|(weight, node)| weight * node * node * node * node)
            .sum();
        assert!((weights.iter().sum::<Real>() - 2.0).abs() < 1.0e-14);
        assert!((quadratic - 2.0 / 3.0).abs() < 1.0e-14);
        assert!((quartic - 2.0 / 5.0).abs() < 1.0e-14);
    }

    #[test]
    fn a_constant_strain_field_gives_constant_stress() {
        let element = element(rectangle(3.0, 2.0, 1, 1));
        assert!((element.area() - 6.0).abs() < 1.0e-12);

        let strain = 1.0e-4;
        // Uniaxial stress: the transverse strain removes the transverse
        // stress, so `n11` is the only resultant left.
        let d = field(&element, |point| {
            let mut values = [0.0; DOF_PER_NODE];
            values[0] = strain * point.x;
            values[1] = -POISSON * strain * point.y;
            values
        });
        let resultant = MODULUS * THICKNESS * strain;
        let centre = element.forces_at(&d, 0.5, 0.5);
        let corner = element.forces_at(&d, 0.9, 0.2);
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

    fn plate_rigidity() -> Real {
        MODULUS * THICKNESS * THICKNESS * THICKNESS / (12.0 * (1.0 - POISSON * POISSON))
    }

    /// Kirchhoff field on a quadratic patch: `w` is quadratic in the in-plane
    /// coordinates and the rotations are its slopes, so the transverse shear
    /// vanishes.
    ///
    /// A quadratic polynomial is not reproduced by sampling it at the control
    /// points, so `w` is written in Bernstein form: `[c/4, -c/4, c/4]` along
    /// each direction gives `c (ξ - 1/2)²`.
    fn curvature_field(
        element: &IgaPatchElement,
        width: Real,
        height: Real,
        k11: Real,
        k22: Real,
    ) -> Vec<Real> {
        let c11 = 0.5 * k11 * width * width;
        let c22 = 0.5 * k22 * height * height;
        let along_u = [0.25 * c11, -0.25 * c11, 0.25 * c11];
        let along_v = [0.25 * c22, -0.25 * c22, 0.25 * c22];
        let n_v = element.patch.n_v();
        let mut d = vec![0.0; element.dofs.len()];
        for (iu, &control_u) in along_u.iter().enumerate() {
            for (iv, &control_v) in along_v.iter().enumerate() {
                let dof = (iu * n_v + iv) * DOF_PER_NODE;
                let x = width * (iu as Real / 2.0 - 0.5);
                let y = height * (iv as Real / 2.0 - 0.5);
                d[dof + 2] = control_u + control_v;
                d[dof + 3] = k22 * y;
                d[dof + 4] = -k11 * x;
            }
        }
        d
    }

    #[test]
    fn a_constant_curvature_field_gives_the_plate_moments() {
        let element = element(rectangle(3.0, 2.0, 2, 2));
        let k11 = 1.0e-3;
        let k22 = -4.0e-4;
        let d = curvature_field(&element, 3.0, 2.0, k11, k22);

        let rigidity = plate_rigidity();
        let allowed = 1.0e-9 * rigidity * (k11.abs() + k22.abs());
        let forces = element.forces_at(&d, 0.5, 0.5);
        let m11 = rigidity * (k11 + POISSON * k22);
        let m22 = rigidity * (k22 + POISSON * k11);
        assert!((forces.m11 - m11).abs() < allowed, "{forces:?}");
        assert!((forces.m22 - m22).abs() < allowed, "{forces:?}");
        assert!(forces.m12.abs() < allowed, "{forces:?}");
    }

    #[test]
    fn a_kirchhoff_field_carries_no_transverse_shear() {
        let element = element(rectangle(3.0, 2.0, 2, 2));
        let k11 = 1.0e-3;
        let k22 = -4.0e-4;
        let d = curvature_field(&element, 3.0, 2.0, k11, k22);

        let allowed = 1.0e-9 * element.ds * (k11.abs() + k22.abs()) * 3.0;
        let forces = element.forces_at(&d, 0.4, 0.7);
        assert!(forces.v13.abs() < allowed, "{forces:?}");
        assert!(forces.v23.abs() < allowed, "{forces:?}");
    }

    #[test]
    fn a_rigid_translation_carries_no_force() {
        let element = element(rectangle(3.0, 2.0, 2, 2));
        let d = field(&element, |_| {
            let mut values = [0.0; DOF_PER_NODE];
            values[0] = 1.5;
            values[1] = -0.5;
            values[2] = 0.7;
            values
        });
        let forces = element.end_forces(&d, Vec3::ZERO);
        let allowed = 1.0e-9 * largest_stiffness(&element);
        for force in forces {
            assert!(force.abs() < allowed, "rigid force {force}");
        }
    }

    #[test]
    fn the_lumped_mass_integrates_the_basis() {
        let element = element(quadratic_strip(3.0, 2.0));
        let mass = element.lumped_mass();
        let expected = 7850.0 * THICKNESS * 6.0 / 1000.0;
        let total: Real = (0..8).map(|cp| mass[cp * DOF_PER_NODE]).sum();
        assert!((total - expected).abs() < 1.0e-9 * expected, "{total}");

        // Across one interior knot the end shares integrate to half of the
        // interior ones; an equal split would make all four the same.
        let end = mass[0];
        let interior = mass[2 * DOF_PER_NODE];
        assert!(
            (interior - 2.0 * end).abs() < 1.0e-9 * interior,
            "end {end} interior {interior}"
        );
        // Rotary inertia carries the same share as translation.
        assert!((mass[2 * DOF_PER_NODE + 3] / mass[3] - 2.0).abs() < 1.0e-9);
    }

    #[test]
    fn a_patch_out_of_plane_uses_its_own_axes() {
        let angle = 0.6;
        let element = element(rotated_rectangle(3.0, 2.0, angle));
        let axes = element.axes();
        let normal = Vec3::new(0.0, -angle.sin(), angle.cos());
        assert!((Vec3::from_array(axes.row(2)) - normal).norm() < 1.0e-12);
        for row in 0..3 {
            for column in 0..3 {
                let a = Vec3::from_array(axes.row(row));
                let b = Vec3::from_array(axes.row(column));
                let expected = if row == column { 1.0 } else { 0.0 };
                assert!((a.dot(b) - expected).abs() < 1.0e-12);
            }
        }

        // Uniaxial stress along the patch's own 1 axis, written as a global
        // displacement field. If the axes or the transform were wrong, the
        // transverse stress would not cancel.
        let e1 = Vec3::from_array(axes.row(0));
        let e2 = Vec3::from_array(axes.row(1));
        let strain = 1.0e-4;
        let points = control_points(&element);
        let mut d = vec![0.0; element.dofs.len()];
        for (cp, point) in points.iter().enumerate() {
            let displacement =
                e1 * (strain * point.dot(e1)) - e2 * (POISSON * strain * point.dot(e2));
            d[cp * DOF_PER_NODE] = displacement.x;
            d[cp * DOF_PER_NODE + 1] = displacement.y;
            d[cp * DOF_PER_NODE + 2] = displacement.z;
        }
        let resultant = MODULUS * THICKNESS * strain;
        let forces = element.forces_at(&d, 0.5, 0.5);
        assert!(
            (forces.n11 - resultant).abs() < 1.0e-9 * resultant,
            "{forces:?}"
        );
        assert!(forces.n22.abs() < 1.0e-9 * resultant, "{forces:?}");
        assert!(forces.n12.abs() < 1.0e-9 * resultant, "{forces:?}");

        let translation = field(&element, |_| {
            let mut values = [0.0; DOF_PER_NODE];
            values[0] = 1.5;
            values[1] = -0.5;
            values[2] = 0.7;
            values
        });
        let allowed = 1.0e-9 * largest_stiffness(&element);
        for force in element.end_forces(&translation, Vec3::ZERO) {
            assert!(force.abs() < allowed, "rigid force {force}");
        }
    }

    #[test]
    fn the_patch_uses_the_shell_dof_order() {
        let patch = rectangle(3.0, 2.0, 2, 2);
        let element = element(patch);
        let formulation: &dyn Element = &element;
        assert_eq!(formulation.dofs().len(), 9 * DOF_PER_NODE);
        assert_eq!(formulation.dofs()[0], 0);
        assert_eq!(formulation.dofs()[2 * DOF_PER_NODE], 2 * DOF_PER_NODE);
        assert_eq!(formulation.stiffness().len(), 54 * 54);
        assert!((element.area() - 6.0).abs() < 1.0e-12);
        assert_eq!(element.thickness(), THICKNESS);
    }

    #[test]
    fn the_stiffness_is_symmetric() {
        let element = element(rectangle(2.0, 1.0, 2, 1));
        let size = element.dofs.len();
        let stiffness = element.local_stiffness();
        let allowed = 1.0e-12 * largest_stiffness(&element);
        for row in 0..size {
            for column in 0..size {
                let a = stiffness[row * size + column];
                let b = stiffness[column * size + row];
                assert!((a - b).abs() <= allowed, "entries {a} and {b} differ");
            }
        }
    }

    #[test]
    fn a_pressure_reaches_every_control_point() {
        let element = element(rectangle(3.0, 2.0, 2, 2));
        let load = element.equivalent_load(Vec3::new(0.0, 0.0, -2.0));
        let total: Real = (0..element.dofs.len() / DOF_PER_NODE)
            .map(|control_point| load[control_point * DOF_PER_NODE + 2])
            .sum();
        assert!((total + 12.0).abs() < 1.0e-9, "{total}");
        for control_point in 0..element.dofs.len() / DOF_PER_NODE {
            assert!(load[control_point * DOF_PER_NODE].abs() < 1.0e-12);
            assert!(load[control_point * DOF_PER_NODE + 1].abs() < 1.0e-12);
        }
    }

    #[test]
    fn a_wrong_node_count_is_rejected() {
        let patch = rectangle(3.0, 2.0, 1, 1);
        assert!(matches!(
            IgaPatchElement::new(patch, THICKNESS, &material(), &[0, 1, 2]),
            Err(ElementError::NodeCount {
                expected: 4,
                count: 3
            })
        ));
    }
}
