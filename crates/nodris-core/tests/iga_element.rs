//! The IGA plate checked against the Navier deflection.
//!
//! One quadratic patch covers the square, with `w` restrained on all four
//! edges and the rotations left free. The element is a Mindlin plate, so its
//! answer carries a small transverse-shear deflection the Kirchhoff series
//! does not have; at `t/a = 1/40` that is a fraction of a percent, and the
//! discretisation error dominates.

use nodris_core::{
    BSplinePatch, BandedMatrix, Element, IgaPatchElement, Material, MaterialId, Vec3,
};

const MODULUS: f64 = 210_000_000.0;
const POISSON: f64 = 0.3;
const THICKNESS: f64 = 0.1;
const SIDE: f64 = 4.0;
/// Uniform pressure [kN/m²], downward.
const PRESSURE: f64 = 1.0;
/// `w_centre = 0.00406 · q·a⁴ / D` for a simply supported square plate.
const NAVIER_COEFFICIENT: f64 = 0.00406;
const DOF_PER_NODE: usize = 6;
/// `w`, `θx`, `θy`; membrane and drilling rotations do not bend the plate.
const KEPT: [usize; 3] = [2, 3, 4];

fn material() -> Material {
    Material::new(MaterialId::new(1), "S355", 7850.0, MODULUS, POISSON)
}

/// A square patch of degree two with `spans` knot spans per direction. The
/// control net has `spans + 2` points per direction, and the returned value is
/// that count.
fn square_patch(spans: usize) -> (BSplinePatch, usize) {
    let control = spans + 2;
    let mut net = Vec::with_capacity(control);
    for row in 0..control {
        let mut points = Vec::with_capacity(control);
        for column in 0..control {
            points.push(Vec3::new(
                SIDE * row as f64 / (control - 1) as f64,
                SIDE * column as f64 / (control - 1) as f64,
                0.0,
            ));
        }
        net.push(points);
    }
    (
        BSplinePatch::from_grid(&net, 2, 2).expect("a square control net"),
        control,
    )
}

/// Solves the simply supported square and returns the centre deflection.
fn centre_deflection(spans: usize) -> f64 {
    let (patch, control) = square_patch(spans);
    let nodes: Vec<usize> = (0..patch.n_cp()).collect();
    let element =
        IgaPatchElement::new(patch, THICKNESS, &material(), &nodes).expect("a flat square patch");

    let mut position = vec![usize::MAX; control * control * KEPT.len()];
    let mut equations = 0;
    for row in 0..control {
        for column in 0..control {
            let on_edge = row == 0 || column == 0 || row + 1 == control || column + 1 == control;
            for (component, &local) in KEPT.iter().enumerate() {
                if on_edge && local == 2 {
                    continue;
                }
                position[(row * control + column) * KEPT.len() + component] = equations;
                equations += 1;
            }
        }
    }

    let element_load = element.equivalent_load(Vec3::new(0.0, 0.0, -PRESSURE));
    let local_size = element.dofs().len();
    let mut free = Vec::with_capacity(equations);
    let mut load = vec![0.0; equations];
    for row in 0..control {
        for column in 0..control {
            for (component, &local) in KEPT.iter().enumerate() {
                let equation = position[(row * control + column) * KEPT.len() + component];
                if equation == usize::MAX {
                    continue;
                }
                free.push((equation, (row * control + column) * DOF_PER_NODE + local));
                load[equation] += element_load[(row * control + column) * DOF_PER_NODE + local];
            }
        }
    }

    let mut stiffness = BandedMatrix::new(equations, equations);
    for (row, &(row_equation, row_local)) in free.iter().enumerate() {
        for &(column_equation, column_local) in free.iter().skip(row) {
            stiffness.add_symmetric(
                row_equation,
                column_equation,
                element.stiffness()[row_local * local_size + column_local],
            );
        }
    }

    let factored = stiffness.factor().expect("the restrained plate is SPD");
    let mut deflection = load;
    factored
        .solve(&mut deflection)
        .expect("the plate system solves");

    let centre = (control / 2) * control + control / 2;
    deflection[position[centre * KEPT.len()]]
}

#[test]
fn a_quadratic_patch_reproduces_the_navier_deflection() {
    let rigidity = MODULUS * THICKNESS * THICKNESS * THICKNESS / (12.0 * (1.0 - POISSON * POISSON));
    let exact = NAVIER_COEFFICIENT * PRESSURE * SIDE * SIDE * SIDE * SIDE / rigidity;

    let coarse = centre_deflection(4).abs();
    let fine = centre_deflection(6).abs();
    let coarse_error = (coarse - exact).abs() / exact;
    let fine_error = (fine - exact).abs() / exact;

    assert!(
        fine_error < coarse_error,
        "6-span error {fine_error} must beat 4-span error {coarse_error}"
    );
    assert!(
        fine_error < 0.015,
        "6-span still {fine_error} off the Navier value {exact}"
    );
}
