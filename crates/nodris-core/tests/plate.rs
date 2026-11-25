//! Simply supported square plate checked against the Navier series.
//!
//! The element is a Mindlin plate, so its converged answer carries a small
//! transverse-shear contribution the Kirchhoff series does not have. At
//! `t/a = 1/40` that residual is a fraction of a percent, and the
//! discretisation error dominates: coarse meshes are too stiff and refinement
//! approaches the closed form from below.

use nodris_core::{BandedMatrix, Element, Material, MaterialId, ShellElement, Vec3};

const MODULUS: f64 = 210_000_000.0;
const POISSON: f64 = 0.3;
const THICKNESS: f64 = 0.1;
const SIDE: f64 = 4.0;
/// Uniform pressure [kN/m²], downward.
const PRESSURE: f64 = 1.0;
/// `w_centre = 0.00406 · q·a⁴ / D` for a simply supported square plate.
const NAVIER_COEFFICIENT: f64 = 0.00406;
const DOF_PER_NODE: usize = 6;
const SHELL_DOF: usize = 4 * DOF_PER_NODE;

/// Local transverse index. The mesh lies in the global xy plane, so the
/// element's local axes are the global ones and these indices address the
/// global stiffness directly.
const W: usize = 2;
/// Components kept: `w`, `θx`, `θy`. Membrane and drilling rotations do not
/// take part in transverse bending.
const KEPT: [usize; 3] = [W, 3, 4];

fn material() -> Material {
    Material::new(MaterialId::new(1), "S355", 7850.0, MODULUS, POISSON)
}

/// Solves an `n × n` mesh with `w` restrained on all four edges and returns
/// the centre deflection.
fn centre_deflection(n: usize) -> f64 {
    let nodes = n + 1;
    let step = SIDE / n as f64;
    let node = |i: usize, j: usize| j * nodes + i;

    // One equation per kept component, with the edge `w` DOFs left out.
    let mut position = vec![usize::MAX; nodes * nodes * KEPT.len()];
    let mut equations = 0;
    for j in 0..nodes {
        for i in 0..nodes {
            let on_edge = i == 0 || j == 0 || i == n || j == n;
            for (component, &local) in KEPT.iter().enumerate() {
                if on_edge && local == W {
                    continue;
                }
                position[node(i, j) * KEPT.len() + component] = equations;
                equations += 1;
            }
        }
    }

    let bandwidth = KEPT.len() * (nodes + 1) + 1;
    let mut stiffness = BandedMatrix::new(equations, bandwidth);
    let mut load = vec![0.0; equations];
    let material = material();

    for j in 0..n {
        for i in 0..n {
            let positions = [
                Vec3::new(i as f64 * step, j as f64 * step, 0.0),
                Vec3::new((i + 1) as f64 * step, j as f64 * step, 0.0),
                Vec3::new((i + 1) as f64 * step, (j + 1) as f64 * step, 0.0),
                Vec3::new(i as f64 * step, (j + 1) as f64 * step, 0.0),
            ];
            let joints = [
                node(i, j),
                node(i + 1, j),
                node(i + 1, j + 1),
                node(i, j + 1),
            ];
            let element = ShellElement::new(positions, THICKNESS, &material, joints)
                .expect("a counter-clockwise square");
            let element_load = element.equivalent_load(Vec3::new(0.0, 0.0, -PRESSURE));

            // Global equation and shell-local DOF of each kept component.
            let mut kept = [(0usize, 0usize); 4 * KEPT.len()];
            for (joint, &mesh_node) in joints.iter().enumerate() {
                for (component, &local) in KEPT.iter().enumerate() {
                    kept[joint * KEPT.len() + component] = (
                        mesh_node * KEPT.len() + component,
                        joint * DOF_PER_NODE + local,
                    );
                }
            }

            for (row, &(row_equation, row_local)) in kept.iter().enumerate() {
                let row_free = position[row_equation];
                if row_free == usize::MAX {
                    continue;
                }
                load[row_free] += element_load[row_local];
                for &(column_equation, column_local) in kept.iter().skip(row) {
                    let column_free = position[column_equation];
                    if column_free == usize::MAX {
                        continue;
                    }
                    let value = element.stiffness()[row_local * SHELL_DOF + column_local];
                    stiffness.add_symmetric(row_free, column_free, value);
                }
            }
        }
    }

    let factored = stiffness.factor().expect("the restrained plate is SPD");
    let mut deflection = load;
    factored
        .solve(&mut deflection)
        .expect("the plate system solves");
    deflection[position[node(n / 2, n / 2) * KEPT.len()]]
}

#[test]
fn a_refined_plate_mesh_converges_to_the_navier_value() {
    let rigidity = MODULUS * THICKNESS * THICKNESS * THICKNESS / (12.0 * (1.0 - POISSON * POISSON));
    let exact = NAVIER_COEFFICIENT * PRESSURE * SIDE * SIDE * SIDE * SIDE / rigidity;

    let coarse = centre_deflection(4).abs();
    let fine = centre_deflection(8).abs();
    let coarse_error = (coarse - exact).abs() / exact;
    let fine_error = (fine - exact).abs() / exact;

    assert!(
        fine_error < coarse_error,
        "8×8 error {fine_error} must beat 4×4 error {coarse_error}"
    );
    assert!(
        fine_error < 0.01,
        "8×8 still {fine_error} off the Navier value {exact}"
    );
}
