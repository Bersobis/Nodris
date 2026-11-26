//! Multi-patch checks: a flat strip and a folded pair.
//!
//! Both meshes share the seam control points, so a seam node is one set of
//! global DOFs. The strip is solved under a uniaxial field that is exact in
//! the patch space; the fold is checked with a rigid rotation, which only
//! balances if each patch transforms its local matrix with its own axes.

use nodris_core::{
    BSplinePatch, BandedMatrix, Element, IgaPatchElement, Material, MaterialId, Vec3,
};

const MODULUS: f64 = 210_000_000.0;
const POISSON: f64 = 0.3;
const THICKNESS: f64 = 0.1;
const PATCH_LENGTH: f64 = 2.0;
const WIDTH: f64 = 2.0;
const STRAIN: f64 = 1.0e-4;
const DOF_PER_NODE: usize = 6;

fn material() -> Material {
    Material::new(MaterialId::new(1), "S355", 7850.0, MODULUS, POISSON)
}

/// A `3 × 3` control net on `x ∈ [x0, x0 + length]`, `y ∈ [0, width]`.
fn flat_net(x0: f64, length: f64) -> Vec<Vec<Vec3>> {
    [0.0, 0.5, 1.0]
        .into_iter()
        .map(|row| {
            [0.0, 0.5, 1.0]
                .into_iter()
                .map(|column| Vec3::new(x0 + length * row, WIDTH * column, 0.0))
                .collect()
        })
        .collect()
}

fn patch(net: &[Vec<Vec3>], degree_u: usize, degree_v: usize) -> BSplinePatch {
    BSplinePatch::from_grid(net, degree_u, degree_v).expect("a rectangular control net")
}

fn element(patch: BSplinePatch, nodes: &[usize]) -> IgaPatchElement {
    IgaPatchElement::new(patch, THICKNESS, &material(), nodes).expect("a flat patch")
}

/// Global stiffness from every element, dense.
fn assemble(elements: &[&IgaPatchElement], ndof: usize) -> Vec<f64> {
    let mut stiffness = vec![0.0; ndof * ndof];
    for element in elements {
        let size = element.dofs().len();
        for (row, &row_dof) in element.dofs().iter().enumerate() {
            for (column, &column_dof) in element.dofs().iter().enumerate() {
                stiffness[row_dof * ndof + column_dof] += element.stiffness()[row * size + column];
            }
        }
    }
    stiffness
}

fn largest(matrix: &[f64]) -> f64 {
    matrix
        .iter()
        .fold(0.0_f64, |peak, value| peak.max(value.abs()))
}

fn residual(matrix: &[f64], ndof: usize, d: &[f64]) -> f64 {
    let mut peak = 0.0_f64;
    for row in 0..ndof {
        let force: f64 = (0..ndof)
            .map(|column| matrix[row * ndof + column] * d[column])
            .sum();
        peak = peak.max(force.abs());
    }
    peak
}

/// Residual on the unconstrained rows only: a prescribed row holds the
/// reaction, not an imbalance.
fn free_residual(matrix: &[f64], ndof: usize, d: &[f64], equation: &[usize]) -> f64 {
    let mut peak = 0.0_f64;
    for (row, &free) in equation.iter().enumerate() {
        if free == usize::MAX {
            continue;
        }
        let force: f64 = (0..ndof)
            .map(|column| matrix[row * ndof + column] * d[column])
            .sum();
        peak = peak.max(force.abs());
    }
    peak
}

fn exact_field(position: Vec3) -> [f64; 2] {
    [STRAIN * position.x, -POISSON * STRAIN * position.y]
}

#[test]
fn a_constant_strain_crosses_the_strip_seam() {
    // Left patch 0..2, right patch 2..4. The right seam row reuses the left
    // seam nodes, as a tie would.
    let left_net = flat_net(0.0, PATCH_LENGTH);
    let right_net = flat_net(PATCH_LENGTH, PATCH_LENGTH);
    let left_nodes = [0, 1, 2, 3, 4, 5, 6, 7, 8];
    let right_nodes = [6, 7, 8, 9, 10, 11, 12, 13, 14];
    let left = element(patch(&left_net, 2, 2), &left_nodes);
    let right = element(patch(&right_net, 2, 2), &right_nodes);

    // A seam node is one seam displacement, not two.
    assert_eq!(&left_nodes[6..], &right_nodes[..3]);
    assert_eq!(left.dofs()[6 * DOF_PER_NODE], right.dofs()[0]);

    let positions: Vec<Vec3> = left_net
        .iter()
        .flatten()
        .chain(right_net.iter().skip(1).flatten())
        .copied()
        .collect();
    assert_eq!(positions.len(), 15);

    let ndof = positions.len() * DOF_PER_NODE;
    let stiffness = assemble(&[&left, &right], ndof);

    // Axial displacement at the two ends, transverse displacement at the
    // bottom edge. The rest is free; the uniaxial field satisfies all three.
    let mut value = vec![0.0; ndof];
    let mut equation = vec![usize::MAX; ndof];
    let mut equations = 0;
    for (node, position) in positions.iter().enumerate() {
        for component in 0..2 {
            let prescribed = if component == 0 {
                if position.x <= 1.0e-12 {
                    Some(0.0)
                } else if position.x >= 2.0 * PATCH_LENGTH - 1.0e-12 {
                    Some(STRAIN * 2.0 * PATCH_LENGTH)
                } else {
                    None
                }
            } else if position.y <= 1.0e-12 {
                Some(0.0)
            } else {
                None
            };
            let dof = node * DOF_PER_NODE + component;
            match prescribed {
                Some(displacement) => value[dof] = displacement,
                None => {
                    equation[dof] = equations;
                    equations += 1;
                }
            }
        }
    }

    let mut reduced = BandedMatrix::new(equations, equations);
    let mut load = vec![0.0; equations];
    for (row, &row_dof) in equation.iter().enumerate() {
        if row_dof == usize::MAX {
            continue;
        }
        for (column, &column_dof) in equation.iter().enumerate() {
            let entry = stiffness[row * ndof + column];
            if column_dof == usize::MAX {
                load[row_dof] -= entry * value[column];
            } else if column_dof >= row_dof {
                reduced.add_symmetric(row_dof, column_dof, entry);
            }
        }
    }

    let factored = reduced.factor().expect("the restrained strip is SPD");
    factored.solve(&mut load).expect("the strip solves");
    for (dof, &free) in equation.iter().enumerate() {
        if free != usize::MAX {
            value[dof] = load[free];
        }
    }

    let scale = STRAIN * 2.0 * PATCH_LENGTH;
    for (node, position) in positions.iter().enumerate() {
        let exact = exact_field(*position);
        assert!(
            (value[node * DOF_PER_NODE] - exact[0]).abs() <= 1.0e-9 * scale,
            "u at {position:?} is {}, expected {}",
            value[node * DOF_PER_NODE],
            exact[0]
        );
        assert!(
            (value[node * DOF_PER_NODE + 1] - exact[1]).abs() <= 1.0e-9 * scale,
            "v at {position:?} is {}, expected {}",
            value[node * DOF_PER_NODE + 1],
            exact[1]
        );
    }

    // Internal forces balance with the exact field on every free row, seam
    // included. The constrained rows carry the end reactions.
    let allowed = 1.0e-9 * largest(&stiffness) * scale * ndof as f64;
    let force = free_residual(&stiffness, ndof, &value, &equation);
    assert!(force < allowed, "strip residual {force}");
}

#[test]
fn a_rigid_rotation_balances_across_the_fold() {
    let left_net = flat_net(0.0, PATCH_LENGTH);
    let folded: Vec<Vec<Vec3>> = [0.0, 0.5, 1.0]
        .into_iter()
        .map(|row| {
            [0.0, 0.5, 1.0]
                .into_iter()
                .map(|column| Vec3::new(PATCH_LENGTH, WIDTH * column, PATCH_LENGTH * row))
                .collect()
        })
        .collect();
    let left_nodes = [0, 1, 2, 3, 4, 5, 6, 7, 8];
    let right_nodes = [6, 7, 8, 9, 10, 11, 12, 13, 14];
    let left = element(patch(&left_net, 2, 2), &left_nodes);
    let right = element(patch(&folded, 2, 2), &right_nodes);

    for column in 0..3 {
        assert!((left_net[2][column] - folded[0][column]).norm() < 1.0e-12);
    }
    assert_eq!(&left_nodes[6..], &right_nodes[..3]);
    assert_eq!(left.dofs()[6 * DOF_PER_NODE], right.dofs()[0]);
    let normal_left = Vec3::from_array(left.axes().row(2));
    let normal_right = Vec3::from_array(right.axes().row(2));
    assert!(normal_left.dot(normal_right).abs() < 1.0e-12);

    let positions: Vec<Vec3> = left_net
        .iter()
        .flatten()
        .chain(folded.iter().skip(1).flatten())
        .copied()
        .collect();
    let ndof = positions.len() * DOF_PER_NODE;
    let stiffness = assemble(&[&left, &right], ndof);

    // Rigid rotation about global Y, including the matching nodal rotation.
    let angle = 0.4;
    let mut displacement = vec![0.0; ndof];
    for (node, position) in positions.iter().enumerate() {
        displacement[node * DOF_PER_NODE] = angle * position.z;
        displacement[node * DOF_PER_NODE + 2] = -angle * position.x;
        displacement[node * DOF_PER_NODE + 4] = angle;
    }

    let allowed = 1.0e-9 * largest(&stiffness) * angle * 4.0 * ndof as f64;
    let force = residual(&stiffness, ndof, &displacement);
    assert!(force < allowed, "rigid rotation force {force}");
}
