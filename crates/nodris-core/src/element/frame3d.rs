//! Local stiffness of a 12-DOF Euler–Bernoulli frame in its own axes.
//!
//! This module builds only the local matrix. Placing it into global axes is a
//! separate step, so a sign or ordering error in the beam formulas can be
//! tested without a coordinate transform.
//!
//! Local DOF order at each end is translation X, Y, Z, then rotation X, Y, Z.

use super::{member_dofs, Element, ElementError, DOF_PER_NODE};
use crate::math::{Mat3, Real, Vec3};
use crate::model::{Material, Section};

const FRAME_DOF: usize = 2 * DOF_PER_NODE;

/// Local 12 × 12 stiffness, row-major, in kN and kN·m.
///
/// `ea` is axial rigidity in kN, `eiy` and `eiz` are bending rigidities in
/// kN·m², `gj` is torsional rigidity in kN·m², and `length` is in metres.
/// Each must be finite and greater than zero.
pub fn local_frame_stiffness(
    ea: Real,
    eiy: Real,
    eiz: Real,
    gj: Real,
    length: Real,
) -> Result<[[Real; FRAME_DOF]; FRAME_DOF], ElementError> {
    for (name, value) in [("EA", ea), ("EIy", eiy), ("EIz", eiz), ("GJ", gj)] {
        if !(value.is_finite() && value > 0.0) {
            return Err(ElementError::InvalidRigidity { name, value });
        }
    }
    if !(length.is_finite() && length > 0.0) {
        return Err(ElementError::InvalidLength { length });
    }

    let axial = ea / length;
    let torsion = gj / length;

    // Deflection in local y is resisted by EIz, with rotation about local z.
    let shear_z = 12.0 * eiz / (length * length * length);
    let couple_z = 6.0 * eiz / (length * length);
    let rot_z = 4.0 * eiz / length;
    let cross_z = 2.0 * eiz / length;

    // Deflection in local z is resisted by EIy, with rotation about local y.
    let shear_y = 12.0 * eiy / (length * length * length);
    let couple_y = 6.0 * eiy / (length * length);
    let rot_y = 4.0 * eiy / length;
    let cross_y = 2.0 * eiy / length;

    let mut k = [[0.0; FRAME_DOF]; FRAME_DOF];

    k[0][0] = axial;
    k[0][6] = -axial;
    k[6][6] = axial;

    k[3][3] = torsion;
    k[3][9] = -torsion;
    k[9][9] = torsion;

    k[1][1] = shear_z;
    k[1][5] = couple_z;
    k[1][7] = -shear_z;
    k[1][11] = couple_z;
    k[5][5] = rot_z;
    k[5][7] = -couple_z;
    k[5][11] = cross_z;
    k[7][7] = shear_z;
    k[7][11] = -couple_z;
    k[11][11] = rot_z;

    k[2][2] = shear_y;
    k[2][4] = -couple_y;
    k[2][8] = -shear_y;
    k[2][10] = -couple_y;
    k[4][4] = rot_y;
    k[4][8] = couple_y;
    k[4][10] = cross_y;
    k[8][8] = shear_y;
    k[8][10] = couple_y;
    k[10][10] = rot_y;

    for row in 0..FRAME_DOF {
        for column in row + 1..FRAME_DOF {
            k[column][row] = k[row][column];
        }
    }

    Ok(k)
}

/// Direction cosines of a member, one local axis per row.
///
/// Local x runs from `p1` to `p2`. Local z stays near global Z, falling back
/// to global X for a vertical member. `rotation` spins the section about x.
pub fn local_axes(p1: Vec3, p2: Vec3, rotation: Real) -> Option<Mat3> {
    if !rotation.is_finite() {
        return None;
    }
    let ex = (p2 - p1).normalized()?;
    let mut ey = Vec3::Z.cross(ex);
    if ey.norm() < 1.0e-8 {
        ey = Vec3::X.cross(ex);
    }
    let mut ey = ey.normalized()?;
    let mut ez = ex.cross(ey);
    if rotation != 0.0 {
        let (sine, cosine) = rotation.sin_cos();
        let turned_y = ey * cosine + ez * sine;
        let turned_z = ez * cosine - ey * sine;
        ey = turned_y;
        ez = turned_z;
    }
    Some(Mat3::from_rows([
        ex.to_array(),
        ey.to_array(),
        ez.to_array(),
    ]))
}

/// Turns the local stiffness into global axes with `axesᵀ K axes`.
///
/// Translation and rotation blocks use the same direction cosines, so the
/// 12 × 12 transform is four copies of the 3 × 3 rotation along the diagonal.
pub fn to_global(
    axes: &Mat3,
    k_local: &[[Real; FRAME_DOF]; FRAME_DOF],
) -> [[Real; FRAME_DOF]; FRAME_DOF] {
    let mut out = [[0.0; FRAME_DOF]; FRAME_DOF];
    for block_row in 0..4 {
        for block_column in 0..4 {
            let rows = [0, 1, 2].map(|row| {
                [0, 1, 2].map(|column| k_local[3 * block_row + row][3 * block_column + column])
            });
            let transformed = axes.congruence(&Mat3::from_rows(rows));
            for row in 0..3 {
                for column in 0..3 {
                    out[3 * block_row + row][3 * block_column + column] =
                        transformed.row(row)[column];
                }
            }
        }
    }
    out
}

fn to_local(axes: &Mat3, d_global: &[Real; FRAME_DOF]) -> [Real; FRAME_DOF] {
    let mut out = [0.0; FRAME_DOF];
    for block in 0..4 {
        let global = Vec3::new(
            d_global[3 * block],
            d_global[3 * block + 1],
            d_global[3 * block + 2],
        );
        let local = axes.mul_vec(global);
        out[3 * block] = local.x;
        out[3 * block + 1] = local.y;
        out[3 * block + 2] = local.z;
    }
    out
}

fn force_to_global(axes: &Mat3, f_local: &[Real; FRAME_DOF]) -> [Real; FRAME_DOF] {
    let transpose = axes.transpose();
    let mut out = [0.0; FRAME_DOF];
    for block in 0..4 {
        let local = Vec3::new(
            f_local[3 * block],
            f_local[3 * block + 1],
            f_local[3 * block + 2],
        );
        let global = transpose.mul_vec(local);
        out[3 * block] = global.x;
        out[3 * block + 1] = global.y;
        out[3 * block + 2] = global.z;
    }
    out
}

/// Equivalent nodal loads of a uniform load `q` in local axes.
fn local_equivalent_load(q: Vec3, length: Real) -> [Real; FRAME_DOF] {
    let squared = length * length;
    let mut load = [0.0; FRAME_DOF];
    load[0] = q.x * length / 2.0;
    load[6] = load[0];
    load[1] = q.y * length / 2.0;
    load[7] = load[1];
    load[2] = q.z * length / 2.0;
    load[8] = load[2];
    load[5] = q.y * squared / 12.0;
    load[11] = -load[5];
    load[4] = -q.z * squared / 12.0;
    load[10] = -load[4];
    load
}

/// Local end forces for global-axis displacements, with `q` in local axes.
pub fn end_forces(
    k_local: &[[Real; FRAME_DOF]; FRAME_DOF],
    axes: &Mat3,
    d_global: &[Real; FRAME_DOF],
    q: Vec3,
    length: Real,
) -> [Real; FRAME_DOF] {
    let local = to_local(axes, d_global);
    let load = local_equivalent_load(q, length);
    let mut out = [0.0; FRAME_DOF];
    for row in 0..FRAME_DOF {
        let stress: Real = (0..FRAME_DOF)
            .map(|column| k_local[row][column] * local[column])
            .sum();
        out[row] = stress - load[row];
    }
    out
}

fn member_mass(density: Real, section: &Section, length: Real) -> [Real; FRAME_DOF] {
    let translational = density * section.area * length / 1000.0;
    let torsional = density * section.torsion_constant * length / 1000.0;
    let about_y = density * section.inertia_y * length / 1000.0;
    let about_z = density * section.inertia_z * length / 1000.0;
    let mut mass = [0.0; FRAME_DOF];
    for component in 0..3 {
        mass[component] = translational / 2.0;
        mass[component + DOF_PER_NODE] = translational / 2.0;
    }
    mass[3] = torsional / 2.0;
    mass[DOF_PER_NODE + 3] = torsional / 2.0;
    mass[4] = about_y / 2.0;
    mass[DOF_PER_NODE + 4] = about_y / 2.0;
    mass[5] = about_z / 2.0;
    mass[DOF_PER_NODE + 5] = about_z / 2.0;
    mass
}

/// A fully prepared two-node frame element.
#[derive(Clone, Debug)]
pub struct FrameElement {
    length: Real,
    axes: Mat3,
    k_local: [[Real; FRAME_DOF]; FRAME_DOF],
    k_global: [Real; FRAME_DOF * FRAME_DOF],
    dofs: [usize; FRAME_DOF],
    mass: [Real; FRAME_DOF],
    ea_over_l: Real,
}

impl FrameElement {
    /// Builds an element from two node positions and their node indices.
    pub fn new(
        p1: Vec3,
        p2: Vec3,
        rotation: Real,
        material: &Material,
        section: &Section,
        start: usize,
        end: usize,
    ) -> Result<Self, ElementError> {
        let axes = local_axes(p1, p2, rotation).ok_or(ElementError::DegenerateMember)?;
        let length = p1.distance(p2);
        if !(length.is_finite() && length > 0.0) {
            return Err(ElementError::InvalidLength { length });
        }
        if !(material.density.is_finite() && material.density >= 0.0) {
            return Err(ElementError::InvalidDensity {
                value: material.density,
            });
        }
        let denominator = 2.0 * (1.0 + material.poisson_ratio);
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(ElementError::InvalidPoissonRatio {
                value: material.poisson_ratio,
            });
        }
        let shear_modulus = material.elastic_modulus / denominator;
        if !(shear_modulus.is_finite() && shear_modulus > 0.0) {
            return Err(ElementError::InvalidPoissonRatio {
                value: material.poisson_ratio,
            });
        }

        let ea = material.elastic_modulus * section.area;
        let eiy = material.elastic_modulus * section.inertia_y;
        let eiz = material.elastic_modulus * section.inertia_z;
        let gj = shear_modulus * section.torsion_constant;
        let k_local = local_frame_stiffness(ea, eiy, eiz, gj, length)?;
        let transformed = to_global(&axes, &k_local);
        let k_global = transformed
            .into_iter()
            .flat_map(|row| row.into_iter())
            .collect::<Vec<_>>()
            .try_into()
            .expect("a 12 by 12 matrix has 144 entries");

        Ok(Self {
            length,
            axes,
            k_local,
            k_global,
            dofs: member_dofs(start, end),
            mass: member_mass(material.density, section, length),
            ea_over_l: ea / length,
        })
    }

    /// The local matrix, kept for post-processing and geometric stiffness.
    pub fn local_stiffness(&self) -> &[[Real; FRAME_DOF]; FRAME_DOF] {
        &self.k_local
    }

    /// Axial rigidity over length, used by axial force recovery.
    pub fn ea_over_l(&self) -> Real {
        self.ea_over_l
    }
}

impl Element for FrameElement {
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
        force_to_global(&self.axes, &local_equivalent_load(q, self.length)).to_vec()
    }

    fn end_forces(&self, d: &[Real], q: Vec3) -> Vec<Real> {
        assert_eq!(d.len(), FRAME_DOF, "frame end forces expect 12 global DOFs");
        let mut global = [0.0; FRAME_DOF];
        global.copy_from_slice(d);
        end_forces(&self.k_local, &self.axes, &global, q, self.length).to_vec()
    }

    fn length(&self) -> Real {
        self.length
    }

    fn lumped_mass(&self) -> &[Real] {
        &self.mass
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::ElementError;
    use crate::math::tolerance::RESIDUAL;

    fn stiffness() -> [[Real; FRAME_DOF]; FRAME_DOF] {
        local_frame_stiffness(100.0, 20.0, 10.0, 5.0, 2.0).unwrap()
    }

    fn force(k: &[[Real; FRAME_DOF]; FRAME_DOF], d: &[Real; FRAME_DOF]) -> [Real; FRAME_DOF] {
        let mut out = [0.0; FRAME_DOF];
        for row in 0..FRAME_DOF {
            for column in 0..FRAME_DOF {
                out[row] += k[row][column] * d[column];
            }
        }
        out
    }

    #[test]
    fn the_local_matrix_is_symmetric() {
        let k = stiffness();
        for row in 0..FRAME_DOF {
            for column in 0..FRAME_DOF {
                assert_eq!(k[row][column], k[column][row]);
            }
        }
    }

    #[test]
    fn entries_match_the_beam_formulas() {
        let k = stiffness();
        assert_eq!(k[0][0], 50.0);
        assert_eq!(k[0][6], -50.0);
        assert_eq!(k[3][3], 2.5);
        assert_eq!(k[1][1], 15.0);
        assert_eq!(k[5][5], 20.0);
        assert_eq!(k[1][11], 15.0);
        assert_eq!(k[2][2], 30.0);
        assert_eq!(k[4][4], 40.0);
        assert_eq!(k[2][4], -30.0);
    }

    #[test]
    fn rigid_translations_and_rotations_produce_no_force() {
        let k = stiffness();
        let length = 2.0;
        let mut translation_x = [0.0; FRAME_DOF];
        let mut translation_y = [0.0; FRAME_DOF];
        let mut translation_z = [0.0; FRAME_DOF];
        let mut rotation_x = [0.0; FRAME_DOF];
        let mut rotation_y = [0.0; FRAME_DOF];
        let mut rotation_z = [0.0; FRAME_DOF];
        translation_x[0] = 1.0;
        translation_x[DOF_PER_NODE] = 1.0;
        rotation_x[3] = 1.0;
        rotation_x[3 + DOF_PER_NODE] = 1.0;
        translation_y[1] = 1.0;
        translation_y[7] = 1.0;
        translation_z[2] = 1.0;
        translation_z[8] = 1.0;
        rotation_y[4] = 1.0;
        rotation_y[10] = 1.0;
        rotation_y[8] = -length;
        rotation_z[5] = 1.0;
        rotation_z[11] = 1.0;
        rotation_z[7] = length;

        let scale = k
            .iter()
            .flatten()
            .fold(0.0_f64, |largest, value| largest.max(value.abs()));
        let allowed = RESIDUAL * scale;
        for (name, mode) in [
            ("translation x", translation_x),
            ("translation y", translation_y),
            ("translation z", translation_z),
            ("rotation x", rotation_x),
            ("rotation y", rotation_y),
            ("rotation z", rotation_z),
        ] {
            for (dof, value) in force(&k, &mode).iter().enumerate() {
                assert!(
                    value.abs() <= allowed,
                    "{name} leaves force {value} at {dof}"
                );
            }
        }
    }

    #[test]
    fn the_matrix_carries_unit_strain_states() {
        let k = stiffness();
        let mut axial = [0.0; FRAME_DOF];
        axial[6] = 1.0;
        let forces = force(&k, &axial);
        assert_eq!(forces[0], -50.0);
        assert_eq!(forces[6], 50.0);
    }

    #[test]
    fn bad_lengths_and_rigidities_are_rejected() {
        assert_eq!(
            local_frame_stiffness(100.0, 20.0, 10.0, 5.0, 0.0),
            Err(ElementError::InvalidLength { length: 0.0 })
        );
        assert_eq!(
            local_frame_stiffness(100.0, -20.0, 10.0, 5.0, 2.0),
            Err(ElementError::InvalidRigidity {
                name: "EIy",
                value: -20.0,
            })
        );
        let nan = local_frame_stiffness(100.0, 20.0, 10.0, f64::NAN, 2.0);
        assert!(matches!(
            nan,
            Err(ElementError::InvalidRigidity { name: "GJ", .. })
        ));
    }

    #[test]
    fn an_error_names_the_bad_input() {
        let length_error = local_frame_stiffness(100.0, 20.0, 10.0, 5.0, -1.0).unwrap_err();
        assert!(length_error.to_string().contains("length"));
        let rigidity_error = local_frame_stiffness(100.0, -20.0, 10.0, 5.0, 2.0).unwrap_err();
        assert!(rigidity_error.to_string().contains("EIy"));
    }

    fn material() -> crate::model::Material {
        crate::model::Material::new(
            crate::model::MaterialId::new(1),
            "S355",
            7850.0,
            210_000_000.0,
            0.3,
        )
    }

    fn section() -> crate::model::Section {
        crate::model::Section::new(
            crate::model::SectionId::new(2),
            "IPE",
            0.01,
            1.0e-5,
            2.0e-5,
            3.0e-5,
            4.0e-7,
        )
    }

    #[test]
    fn local_axes_follow_the_member_and_its_rotation() {
        let horizontal = local_axes(Vec3::ZERO, Vec3::new(2.0, 0.0, 0.0), 0.0).unwrap();
        assert_eq!(horizontal.row(0), [1.0, 0.0, 0.0]);
        assert_eq!(horizontal.row(1), [0.0, 1.0, 0.0]);
        assert_eq!(horizontal.row(2), [0.0, 0.0, 1.0]);

        let vertical = local_axes(Vec3::ZERO, Vec3::new(0.0, 0.0, 2.0), 0.0).unwrap();
        assert_eq!(vertical.row(0), [0.0, 0.0, 1.0]);
        assert!((vertical.det() - 1.0).abs() < RESIDUAL);

        let turned = local_axes(
            Vec3::ZERO,
            Vec3::new(2.0, 0.0, 0.0),
            std::f64::consts::FRAC_PI_2,
        )
        .unwrap();
        for (got, want) in turned.row(1).iter().zip([0.0, 0.0, 1.0]) {
            assert!((got - want).abs() < 1.0e-12);
        }
        assert!(local_axes(Vec3::ZERO, Vec3::ZERO, 0.0).is_none());
        assert!(local_axes(Vec3::ZERO, Vec3::X, f64::NAN).is_none());
    }

    #[test]
    fn a_rigid_translation_is_force_free_for_a_rotated_member() {
        let axes = local_axes(Vec3::ZERO, Vec3::new(0.0, 3.0, 0.0), 0.4).unwrap();
        let mut d = [0.0; FRAME_DOF];
        for node in 0..2 {
            let base = node * DOF_PER_NODE;
            d[base] = 1.5;
            d[base + 1] = -0.5;
            d[base + 2] = 2.0;
        }
        let forces = end_forces(&stiffness(), &axes, &d, Vec3::ZERO, 2.0);
        let scale = stiffness()
            .iter()
            .flatten()
            .fold(0.0_f64, |largest, value| largest.max(value.abs()));
        for force in forces {
            assert!(force.abs() <= RESIDUAL * scale, "rigid force {force}");
        }
    }

    #[test]
    fn identity_axes_leave_the_local_matrix_unchanged() {
        let axes = local_axes(Vec3::ZERO, Vec3::new(2.0, 0.0, 0.0), 0.0).unwrap();
        assert_eq!(to_global(&axes, &stiffness()), stiffness());
    }

    #[test]
    fn a_cantilever_tip_load_satisfies_local_equilibrium() {
        let length = 2.0;
        let eiz = 10.0;
        let load = 3.0;
        let axes = local_axes(Vec3::ZERO, Vec3::new(length, 0.0, 0.0), 0.0).unwrap();
        let k = stiffness();
        let mut d = [0.0; FRAME_DOF];
        d[7] = load * length * length * length / (3.0 * eiz);
        d[11] = load * length * length / (2.0 * eiz);

        let forces = end_forces(&k, &axes, &d, Vec3::ZERO, length);
        for (got, want) in [
            (forces[7], load),
            (forces[1], -load),
            (forces[5], -load * length),
        ] {
            assert!((got - want).abs() < 1.0e-9, "got {got}, want {want}");
        }
    }

    #[test]
    fn a_uniform_load_balances_its_equivalent_nodal_loads() {
        let length = 2.0;
        let axes = local_axes(Vec3::ZERO, Vec3::new(length, 0.0, 0.0), 0.0).unwrap();
        let q = Vec3::new(0.0, 4.0, 0.0);
        let forces = end_forces(&stiffness(), &axes, &[0.0; FRAME_DOF], q, length);
        assert!((forces[1] + 4.0).abs() < 1.0e-12);
        assert!((forces[7] + 4.0).abs() < 1.0e-12);
        assert!((forces[1] + forces[7] + q.y * length).abs() < 1.0e-12);
    }

    #[test]
    fn a_frame_element_reaches_every_trait_hook() {
        let element = FrameElement::new(
            Vec3::ZERO,
            Vec3::new(2.0, 0.0, 0.0),
            0.0,
            &material(),
            &section(),
            0,
            1,
        )
        .unwrap();
        let formulation: &dyn Element = &element;
        assert_eq!(formulation.length(), 2.0);
        assert_eq!(formulation.dofs(), member_dofs(0, 1));
        assert_eq!(formulation.stiffness().len(), FRAME_DOF * FRAME_DOF);
        assert_eq!(formulation.lumped_mass().len(), FRAME_DOF);
        assert_eq!(formulation.axes(), &Mat3::IDENTITY);
        assert!((element.local_stiffness()[0][0] - 1.05e6).abs() < 1.0);
        assert!((element.ea_over_l() - 1.05e6).abs() < 1.0);
        assert_eq!(formulation.equivalent_load(Vec3::Y)[1], 1.0);
        let forces = formulation.end_forces(&[0.0; FRAME_DOF], Vec3::Y);
        assert!((forces[1] + 1.0).abs() < 1.0e-9);
    }

    #[test]
    fn unusable_material_inputs_are_rejected() {
        let section = section();
        let mut bad_poisson = material();
        bad_poisson.poisson_ratio = -1.0;
        assert!(matches!(
            FrameElement::new(
                Vec3::ZERO,
                Vec3::new(2.0, 0.0, 0.0),
                0.0,
                &bad_poisson,
                &section,
                0,
                1,
            ),
            Err(ElementError::InvalidPoissonRatio { .. })
        ));
        let mut bad_density = material();
        bad_density.density = -1.0;
        assert!(matches!(
            FrameElement::new(
                Vec3::ZERO,
                Vec3::new(2.0, 0.0, 0.0),
                0.0,
                &bad_density,
                &section,
                0,
                1,
            ),
            Err(ElementError::InvalidDensity { .. })
        ));
    }
}
