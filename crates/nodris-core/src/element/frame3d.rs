//! Local stiffness of a 12-DOF Euler–Bernoulli frame in its own axes.
//!
//! This module builds only the local matrix. Placing it into global axes is a
//! separate step, so a sign or ordering error in the beam formulas can be
//! tested without a coordinate transform.
//!
//! Local DOF order at each end is translation X, Y, Z, then rotation X, Y, Z.

use super::{ElementError, DOF_PER_NODE};
use crate::math::Real;

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
}
