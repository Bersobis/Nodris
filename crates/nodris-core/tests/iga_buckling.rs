//! Euler buckling of a simply supported strip on the quadratic patch.
//!
//! The strip is narrow, with its rotations about the long axis and its
//! in-plane displacements restrained, so it bends in one plane only. Its
//! critical load per unit width is then `π²D/L²`, the Euler value for a beam
//! of unit width.

use nodris_core::{
    BSplinePatch, Element, IgaPatchElement, InitialStress, Material, MaterialId, Vec3,
};

const MODULUS: f64 = 210_000_000.0;
const POISSON: f64 = 0.3;
const THICKNESS: f64 = 0.05;
const LENGTH: f64 = 2.0;
const WIDTH: f64 = 0.2;
const SPANS: usize = 15;
const DOF_PER_NODE: usize = 6;

fn material() -> Material {
    Material::new(MaterialId::new(1), "S355", 7850.0, MODULUS, POISSON)
}

fn strip() -> IgaPatchElement {
    let control = SPANS + 2;
    let mut net = Vec::with_capacity(control);
    for row in 0..control {
        let x = LENGTH * row as f64 / (control - 1) as f64;
        net.push(vec![Vec3::new(x, 0.0, 0.0), Vec3::new(x, WIDTH, 0.0)]);
    }
    let patch = BSplinePatch::from_grid(&net, 2, 1).expect("a strip control net");
    let nodes: Vec<usize> = (0..patch.n_cp()).collect();
    IgaPatchElement::new(patch, THICKNESS, &material(), &nodes).expect("a flat strip")
}

fn quadratic_form(matrix: &[f64], size: usize, free: &[usize], d: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (row, &row_dof) in free.iter().enumerate() {
        let stress: f64 = free
            .iter()
            .enumerate()
            .map(|(column, &column_dof)| matrix[row_dof * size + column_dof] * d[column])
            .sum();
        sum += d[row] * stress;
    }
    sum
}

#[test]
fn a_compressed_strip_reaches_the_euler_load() {
    let element = strip();
    let control = SPANS + 2;

    // Simply supported ends: `w` fixed at x = 0 and x = L, rotations free.
    // Everything else is restrained, so the trial field has `w` and `θy`.
    let mut free = Vec::new();
    for row in 0..control {
        for column in 0..2 {
            let base = (row * 2 + column) * DOF_PER_NODE;
            if row != 0 && row + 1 != control {
                free.push(base + 2);
            }
            free.push(base + 4);
        }
    }
    let displacement: Vec<f64> = free
        .iter()
        .map(|&dof| {
            let x = LENGTH * (dof / DOF_PER_NODE / 2) as f64 / (control - 1) as f64;
            let wave = std::f64::consts::PI * x / LENGTH;
            if dof % DOF_PER_NODE == 2 {
                wave.sin()
            } else {
                -(std::f64::consts::PI / LENGTH) * wave.cos()
            }
        })
        .collect();

    let size = element.dofs().len();
    let geometric = element.geometric_stiffness_matrix(InitialStress::membrane(-1.0, 0.0, 0.0));
    let numerator = quadratic_form(element.stiffness(), size, &free, &displacement);
    let denominator = -quadratic_form(&geometric, size, &free, &displacement);
    assert!(denominator > 0.0, "compression must add a negative term");
    let critical = numerator / denominator;

    let rigidity = MODULUS * THICKNESS * THICKNESS * THICKNESS / (12.0 * (1.0 - POISSON * POISSON));
    let euler = std::f64::consts::PI * std::f64::consts::PI * rigidity / (LENGTH * LENGTH);
    assert!(
        (critical - euler).abs() <= 0.01 * euler,
        "patch critical load {critical}, Euler {euler}"
    );

    let tension = element.geometric_stiffness_matrix(InitialStress::membrane(1.0, 0.0, 0.0));
    for (compressed, stretched) in geometric.iter().zip(&tension) {
        assert_eq!(*compressed, -stretched);
    }
}

#[test]
fn a_zero_stress_field_has_no_geometric_stiffness() {
    let element = strip();
    let matrix = element.geometric_stiffness_matrix(InitialStress::default());
    assert_eq!(matrix.len(), element.dofs().len() * element.dofs().len());
    assert!(matrix.iter().all(|&value| value == 0.0));
    assert!(element
        .geometric_stiffness(InitialStress::default())
        .is_some());
}
