//! Two-node axial bar.
//!
//! A truss is a modelling intent on a member, not a separate model entity.
//! The element contributes axial stiffness only, so the rotational DOFs at a
//! truss-only node keep a zero diagonal. The assembler deactivates those DOFs
//! rather than giving them a penalty stiffness.

use super::{member_dofs, Element, ElementError, InitialStress, DOF_PER_NODE};
use crate::element::frame3d::local_axes;
use crate::math::{Mat3, Real, Vec3};
use crate::model::{Material, Section};

const BAR_DOF: usize = 2 * DOF_PER_NODE;

/// A two-node bar carrying axial force only.
#[derive(Clone, Debug)]
pub struct TrussElement {
    length: Real,
    axes: Mat3,
    k_global: [Real; BAR_DOF * BAR_DOF],
    ea_over_l: Real,
    dofs: [usize; BAR_DOF],
    mass: [Real; BAR_DOF],
}

impl TrussElement {
    /// Builds a bar from two node positions and their node indices.
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

        let ea = material.elastic_modulus * section.area;
        if !(ea.is_finite() && ea > 0.0) {
            return Err(ElementError::InvalidRigidity {
                name: "EA",
                value: ea,
            });
        }
        let ea_over_l = ea / length;
        let axis = axes.row(0);
        let direction = [axis[0], axis[1], axis[2]];

        // K = (EA/L) · c cᵀ in the translational blocks; rotations stay zero.
        let mut k = [0.0; BAR_DOF * BAR_DOF];
        for row in 0..3 {
            for column in 0..3 {
                let entry = ea_over_l * direction[row] * direction[column];
                k[row * BAR_DOF + column] = entry;
                k[row * BAR_DOF + column + DOF_PER_NODE] = -entry;
                k[(row + DOF_PER_NODE) * BAR_DOF + column] = -entry;
                k[(row + DOF_PER_NODE) * BAR_DOF + column + DOF_PER_NODE] = entry;
            }
        }

        let mass_each = material.density * section.area * length / 2000.0;
        let mut mass = [0.0; BAR_DOF];
        for component in 0..3 {
            mass[component] = mass_each;
            mass[component + DOF_PER_NODE] = mass_each;
        }

        Ok(Self {
            length,
            axes,
            k_global: k,
            ea_over_l,
            dofs: member_dofs(start, end),
            mass,
        })
    }

    /// Axial rigidity over length, used by force recovery.
    pub fn ea_over_l(&self) -> Real {
        self.ea_over_l
    }
}

impl Element for TrussElement {
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
        let axis = self.axes.row(0);
        let at_each_end = q.x * self.length / 2.0;
        let mut load = vec![0.0; BAR_DOF];
        for component in 0..3 {
            let global = at_each_end * axis[component];
            load[component] = global;
            load[component + DOF_PER_NODE] = global;
        }
        load
    }

    fn end_forces(&self, d: &[Real], q: Vec3) -> Vec<Real> {
        assert_eq!(d.len(), BAR_DOF, "bar end forces expect 12 global DOFs");
        let axis = self.axes.row(0);
        let mut elongation = 0.0;
        for component in 0..3 {
            let direction = axis[component];
            elongation += direction * (d[component + DOF_PER_NODE] - d[component]);
        }
        let axial = self.ea_over_l * elongation - q.x * self.length / 2.0;
        let mut forces = vec![0.0; BAR_DOF];
        forces[0] = -axial;
        forces[DOF_PER_NODE] = axial;
        forces
    }

    fn length(&self) -> Real {
        self.length
    }

    fn lumped_mass(&self) -> &[Real] {
        &self.mass
    }

    fn geometric_stiffness(&self, stress: InitialStress) -> Option<Vec<Real>> {
        let axis = self.axes.row(0);
        let direction = [axis[0], axis[1], axis[2]];
        let factor = stress.axial / self.length;
        let mut kg = [0.0; BAR_DOF * BAR_DOF];
        for row in 0..3 {
            for column in 0..3 {
                let identity = if row == column { 1.0 } else { 0.0 };
                let entry = factor * (identity - direction[row] * direction[column]);
                kg[row * BAR_DOF + column] = entry;
                kg[row * BAR_DOF + column + DOF_PER_NODE] = -entry;
                kg[(row + DOF_PER_NODE) * BAR_DOF + column] = -entry;
                kg[(row + DOF_PER_NODE) * BAR_DOF + column + DOF_PER_NODE] = entry;
            }
        }
        Some(kg.to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::DOF_NAMES;

    fn material() -> Material {
        Material::new(
            crate::model::MaterialId::new(1),
            "S355",
            7850.0,
            210_000_000.0,
            0.3,
        )
    }

    fn section() -> Section {
        Section::new(
            crate::model::SectionId::new(2),
            "rod",
            0.01,
            0.0,
            0.0,
            0.0,
            0.0,
        )
    }

    #[test]
    fn a_bar_has_no_rotational_stiffness_or_mass() {
        let bar = TrussElement::new(
            Vec3::ZERO,
            Vec3::new(2.0, 0.0, 0.0),
            0.0,
            &material(),
            &section(),
            0,
            1,
        )
        .unwrap();
        let formulation: &dyn Element = &bar;
        for (index, name) in DOF_NAMES.iter().enumerate() {
            if name.starts_with('r') {
                assert_eq!(bar.k_global[(index) * BAR_DOF + index], 0.0);
                assert_eq!(
                    bar.k_global[(index + DOF_PER_NODE) * BAR_DOF + index + DOF_PER_NODE],
                    0.0
                );
                assert_eq!(bar.mass[index], 0.0);
            }
        }
        assert_eq!(formulation.stiffness().len(), BAR_DOF * BAR_DOF);
        assert_eq!(bar.k_global[0], bar.ea_over_l());
        assert_eq!(formulation.dofs(), member_dofs(0, 1));
        assert_eq!(formulation.length(), 2.0);
        assert_eq!(formulation.equivalent_load(Vec3::X)[0], 1.0);
    }

    #[test]
    fn a_two_bar_truss_carries_only_axial_force() {
        let e_a = material().elastic_modulus * 0.01;
        let length = 2.0_f64.sqrt();
        let apex = Vec3::new(0.0, 1.0, 0.0);
        let left = TrussElement::new(
            Vec3::new(-1.0, 0.0, 0.0),
            apex,
            0.0,
            &material(),
            &section(),
            0,
            2,
        )
        .unwrap();
        let right = TrussElement::new(
            Vec3::new(1.0, 0.0, 0.0),
            apex,
            0.0,
            &material(),
            &section(),
            1,
            2,
        )
        .unwrap();

        let mut d = [0.0; BAR_DOF];
        d[7] = 0.01;
        let mut vertical = 0.0;
        let mut horizontal = 0.0;
        for bar in [&left, &right] {
            let local = bar.end_forces(&d, Vec3::ZERO);
            assert_eq!(local[1], 0.0);
            assert_eq!(local[2], 0.0);
            assert_eq!(local[3], 0.0);
            let global = bar
                .axes()
                .transpose()
                .mul_vec(Vec3::new(local[6], local[7], local[8]));
            horizontal += global.x;
            vertical += global.y;
        }
        assert!(
            horizontal.abs() < 1.0e-6,
            "horizontal residual {horizontal}"
        );
        let expected = e_a * 0.01 / length;
        assert!(
            (vertical - expected).abs() < 1.0e-6,
            "vertical {vertical}, expected {expected}"
        );
    }

    #[test]
    fn compression_reduces_transverse_stiffness_but_not_axial() {
        let bar = TrussElement::new(
            Vec3::ZERO,
            Vec3::new(2.0, 0.0, 0.0),
            0.0,
            &material(),
            &section(),
            0,
            1,
        )
        .unwrap();
        let kg = bar
            .geometric_stiffness(InitialStress::axial(-10.0))
            .unwrap();
        assert_eq!(kg[0], 0.0);
        assert!((kg[7 * BAR_DOF + 7] + 5.0).abs() < 1.0e-9);
        assert_eq!(kg[3 * BAR_DOF + 3], 0.0);
    }

    #[test]
    fn bad_inputs_are_rejected() {
        assert!(matches!(
            TrussElement::new(Vec3::ZERO, Vec3::ZERO, 0.0, &material(), &section(), 0, 1,),
            Err(ElementError::DegenerateMember)
        ));
        let mut material = material();
        material.density = -1.0;
        assert!(matches!(
            TrussElement::new(Vec3::ZERO, Vec3::X, 0.0, &material, &section(), 0, 1,),
            Err(ElementError::InvalidDensity { .. })
        ));
    }
}
