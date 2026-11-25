//! Closed-form checks for the frame and truss formulations.
//!
//! Equilibrium tests alone cannot catch a sign or transformation error that
//! happens to balance. These cases compare against formulas derived outside
//! the element code.

use nodris_core::{
    BandedMatrix, Element, FrameElement, Material, MaterialId, Section, SectionId, TrussElement,
    Vec3,
};

const E: f64 = 210_000_000.0;
const AREA: f64 = 0.01;
const INERTIA: f64 = 1.0e-4;
const TORSION: f64 = 1.0e-4;

fn material() -> Material {
    Material::new(MaterialId::new(1), "S355", 7850.0, E, 0.3)
}

fn section() -> Section {
    Section::new(
        SectionId::new(2),
        "reference",
        AREA,
        INERTIA,
        INERTIA,
        TORSION,
        0.0,
    )
}

/// Solves a reduced system assembled from the free DOFs of one element.
fn solve_reduced(stiffness: &[f64], dimension: usize, free: &[usize], load: &[f64]) -> Vec<f64> {
    let mut matrix = BandedMatrix::new(free.len(), free.len());
    for (row, &outer) in free.iter().enumerate() {
        for (column, &inner) in free.iter().enumerate().skip(row) {
            matrix.add_symmetric(row, column, stiffness[outer * dimension + inner]);
        }
    }
    let factored = matrix.factor().expect("reduced stiffness must be SPD");
    let mut rhs = load.to_vec();
    factored.solve(&mut rhs).expect("reduced system must solve");
    rhs
}

#[test]
fn cantilever_tip_deflection_matches_pl_cubed_over_three_ei() {
    let length = 2.0;
    let load = 10.0;
    let expected = load * length * length * length / (3.0 * E * INERTIA);
    let element = FrameElement::new(
        Vec3::ZERO,
        Vec3::new(length, 0.0, 0.0),
        0.0,
        &material(),
        &section(),
        0,
        1,
    )
    .unwrap();

    let free: Vec<usize> = (0..6).map(|offset| 6 + offset).collect();
    let mut rhs = vec![0.0; free.len()];
    rhs[1] = load;
    let displacement = solve_reduced(element.stiffness(), 12, &free, &rhs);

    assert!(
        (displacement[1] - expected).abs() <= 1.0e-9 * expected,
        "tip deflection {}, closed form {expected}",
        displacement[1]
    );
    assert!(displacement[0].abs() < 1.0e-12);
    assert!(displacement[2].abs() < 1.0e-12);
}

#[test]
fn fixed_fixed_end_moments_match_the_beam_formulas() {
    let length = 3.0;
    let pressure = 4.0;
    let element = FrameElement::new(
        Vec3::ZERO,
        Vec3::new(length, 0.0, 0.0),
        0.0,
        &material(),
        &section(),
        0,
        1,
    )
    .unwrap();

    let forces = element.end_forces(&[0.0; 12], Vec3::new(0.0, pressure, 0.0));
    let expected = pressure * length * length / 12.0;
    assert!(
        (forces[5] + expected).abs() <= 1.0e-9 * expected,
        "left end moment {}, expected -{expected}",
        forces[5]
    );
    assert!(
        (forces[11] - expected).abs() <= 1.0e-9 * expected,
        "right end moment {}, expected {expected}",
        forces[11]
    );
    assert!((forces[1] + pressure * length / 2.0).abs() < 1.0e-9);
    assert!((forces[7] + pressure * length / 2.0).abs() < 1.0e-9);
}

#[test]
fn a_two_bar_truss_matches_the_method_of_joints() {
    let load = 8.0;
    let length = 2.0_f64.sqrt();
    let apex_y = 1.0;
    let left = TrussElement::new(
        Vec3::new(-1.0, 0.0, 0.0),
        Vec3::new(0.0, apex_y, 0.0),
        0.0,
        &material(),
        &section(),
        0,
        2,
    )
    .unwrap();
    let right = TrussElement::new(
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, apex_y, 0.0),
        0.0,
        &material(),
        &section(),
        1,
        2,
    )
    .unwrap();

    // The apex is the end node of both bars. Out-of-plane motion is fixed for
    // this planar truss, leaving the apex x and y translations free.
    let free = [6, 7];
    let mut rhs = vec![0.0, -load];
    let mut matrix = BandedMatrix::new(2, 2);
    for bar in [&left, &right] {
        let stiffness = bar.stiffness();
        for (row, &outer) in free.iter().enumerate() {
            for (column, &inner) in free.iter().enumerate().skip(row) {
                matrix.add_symmetric(row, column, stiffness[outer * 12 + inner]);
            }
        }
    }
    let factored = matrix.factor().expect("truss stiffness must be SPD");
    factored.solve(&mut rhs).expect("truss system must solve");

    let mut d = [0.0; 12];
    d[6] = rhs[0];
    d[7] = rhs[1];
    let axial_expected = -load / (2.0 * (apex_y / length));
    for bar in [&left, &right] {
        let forces = bar.end_forces(&d, Vec3::ZERO);
        assert!(
            (forces[6] - axial_expected).abs() <= 1.0e-8 * axial_expected.abs(),
            "axial force {}, method of joints {axial_expected}",
            forces[6]
        );
        assert_eq!(forces[7], 0.0);
        assert_eq!(forces[8], 0.0);
    }
    assert!(rhs[0].abs() < 1.0e-12, "apex should not sway sideways");
}
