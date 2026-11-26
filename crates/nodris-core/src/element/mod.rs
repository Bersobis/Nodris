//! Element formulations.
//!
//! Every formulation implements [`Element`], the contract the analysis
//! pipeline sees. An assembler reads the DOF numbers, the stiffness and the
//! load hook without branching on the formulation.
//!
//! There are six degrees of freedom per node, always in this order:
//! translation X, Y, Z, then rotation about X, Y, Z. [`node_dofs`] is the one
//! place that turns a node index into those global numbers.
//!
//! The full contract, including the local-axis convention, the sign
//! conventions for end forces and the meaning of a `None` geometric
//! stiffness, is written down in
//! [NUMERICAL-KERNEL-ELEMENT-INTENT.md](https://github.com/Bersobis/Nodris/blob/main/docs/NUMERICAL-KERNEL-ELEMENT-INTENT.md).

use std::fmt;

use crate::math::{Mat3, Real, Vec3};
use crate::model::Section;

mod frame3d;
mod iga_patch;
mod shell4;
mod truss3d;

pub use frame3d::{end_forces, local_axes, local_frame_stiffness, to_global, FrameElement};
pub use iga_patch::IgaPatchElement;
pub use shell4::ShellElement;
pub use truss3d::TrussElement;

/// Invalid input to an element formulation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ElementError {
    /// A member length must be finite and greater than zero.
    InvalidLength { length: Real },
    /// A rigidity must be finite and greater than zero.
    InvalidRigidity { name: &'static str, value: Real },
    /// The two nodes do not define a direction.
    DegenerateMember,
    /// Poisson ratio makes the shear modulus unusable.
    InvalidPoissonRatio { value: Real },
    /// Density must be finite and not negative.
    InvalidDensity { value: Real },
    /// Thickness must be finite and greater than zero.
    InvalidThickness { value: Real },
    /// The joints do not define a plane.
    DegenerateArea,
    /// The joints run clockwise or cross over.
    ClockwiseOrdering,
    /// The node list does not cover every control point.
    NodeCount {
        /// Control points in the patch.
        expected: usize,
        /// Node indices supplied.
        count: usize,
    },
}

impl fmt::Display for ElementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLength { length } => {
                write!(
                    f,
                    "length must be finite and greater than zero, got {length}"
                )
            }
            Self::InvalidRigidity { name, value } => {
                write!(
                    f,
                    "{name} must be finite and greater than zero, got {value}"
                )
            }
            Self::DegenerateMember => write!(f, "member nodes do not define a direction"),
            Self::InvalidPoissonRatio { value } => {
                write!(f, "Poisson ratio {value} gives an unusable shear modulus")
            }
            Self::InvalidDensity { value } => {
                write!(f, "density must be finite and not negative, got {value}")
            }
            Self::InvalidThickness { value } => {
                write!(
                    f,
                    "thickness must be finite and greater than zero, got {value}"
                )
            }
            Self::DegenerateArea => write!(f, "the joints do not define a plane"),
            Self::ClockwiseOrdering => {
                write!(f, "the joints run clockwise or cross over")
            }
            Self::NodeCount { expected, count } => {
                write!(
                    f,
                    "a patch has {expected} control points, {count} nodes were given"
                )
            }
        }
    }
}

impl std::error::Error for ElementError {}

/// The stress state a geometric stiffness is built from.
///
/// A member uses `axial` in kN, tension positive. A shell or patch uses the
/// membrane resultants in kN/m, also positive in tension. Keeping both in one
/// value lets every formulation use the same trait method.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct InitialStress {
    pub axial: Real,
    pub nx: Real,
    pub ny: Real,
    pub nxy: Real,
}

impl InitialStress {
    /// Axial force in a bar or frame, tension positive.
    pub fn axial(axial: Real) -> Self {
        Self {
            axial,
            ..Self::default()
        }
    }

    /// Membrane resultants in a shell or patch, tension positive.
    pub fn membrane(nx: Real, ny: Real, nxy: Real) -> Self {
        Self {
            nx,
            ny,
            nxy,
            ..Self::default()
        }
    }

    pub fn is_zero(&self) -> bool {
        self.axial == 0.0 && self.nx == 0.0 && self.ny == 0.0 && self.nxy == 0.0
    }
}

/// Internal force resultants at one point of a shell or patch, per unit width.
///
/// The four-node shell and the isogeometric patch both report this type, so
/// result post-processing does not depend on which area formulation produced
/// the section forces.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ShellForces {
    /// Membrane force in local 1, tension positive.
    pub n11: Real,
    /// Membrane force in local 2, tension positive.
    pub n22: Real,
    /// Membrane shear.
    pub n12: Real,
    /// Bending moment about local 2.
    pub m11: Real,
    /// Bending moment about local 1.
    pub m22: Real,
    /// Twisting moment.
    pub m12: Real,
    /// Transverse shear on the local 1 face.
    pub v13: Real,
    /// Transverse shear on the local 2 face.
    pub v23: Real,
}

/// Number of degrees of freedom at every node.
pub const DOF_PER_NODE: usize = 6;

/// The six nodal degrees of freedom in their fixed order.
pub const DOF_NAMES: [&str; DOF_PER_NODE] = ["tx", "ty", "tz", "rx", "ry", "rz"];

/// Global DOF numbers for one node, in the order given by [`DOF_NAMES`].
pub const fn node_dofs(node: usize) -> [usize; DOF_PER_NODE] {
    let mut dofs = [0; DOF_PER_NODE];
    let mut component = 0;
    while component < DOF_PER_NODE {
        dofs[component] = node * DOF_PER_NODE + component;
        component += 1;
    }
    dofs
}

/// Global DOF numbers for a two-node member, start then end.
pub const fn member_dofs(start: usize, end: usize) -> [usize; 2 * DOF_PER_NODE] {
    let mut dofs = [0; 2 * DOF_PER_NODE];
    let mut component = 0;
    while component < DOF_PER_NODE {
        dofs[component] = start * DOF_PER_NODE + component;
        dofs[component + DOF_PER_NODE] = end * DOF_PER_NODE + component;
        component += 1;
    }
    dofs
}

/// Lumped nodal mass of a member, in tonnes, in local DOF order.
///
/// Half the mass goes to each end. Rotary inertia is included even though
/// Euler–Bernoulli theory drops it: the eigenproblem needs a positive mass
/// matrix. A truss leaves it out, because its rotational degrees of freedom
/// are deactivated by the assembler.
pub(crate) fn member_lumped_mass(
    density: Real,
    section: &Section,
    length: Real,
    rotary: bool,
) -> [Real; 2 * DOF_PER_NODE] {
    let translational = density * section.area * length / 1000.0;
    let mut mass = [0.0; 2 * DOF_PER_NODE];
    for component in 0..3 {
        mass[component] = translational / 2.0;
        mass[component + DOF_PER_NODE] = translational / 2.0;
    }
    if rotary {
        let torsional = density * section.torsion_constant * length / 1000.0;
        let about_y = density * section.inertia_y * length / 1000.0;
        let about_z = density * section.inertia_z * length / 1000.0;
        mass[3] = torsional / 2.0;
        mass[DOF_PER_NODE + 3] = torsional / 2.0;
        mass[4] = about_y / 2.0;
        mass[DOF_PER_NODE + 4] = about_y / 2.0;
        mass[5] = about_z / 2.0;
        mass[DOF_PER_NODE + 5] = about_z / 2.0;
    }
    mass
}

/// The contract between a formulation and the analysis pipeline.
///
/// With `n = dofs().len()`, every vector below is in that local order. A
/// frame has `n = 12`, but nothing outside the formulation relies on it.
pub trait Element {
    /// Global DOF number for each local DOF, in local order.
    fn dofs(&self) -> &[usize];

    /// Global stiffness, row-major `n × n` with `n = dofs().len()`.
    fn stiffness(&self) -> &[Real];

    /// Direction cosines of the local axes, one axis per row.
    fn axes(&self) -> &Mat3;

    /// Consistent nodal load in global axes for `q` in local axes [kN/m].
    fn equivalent_load(&self, q: Vec3) -> Vec<Real>;

    /// Local end forces for global-axis displacements `d` in local DOF order,
    /// with load `q` in local axes [kN/m].
    fn end_forces(&self, d: &[Real], q: Vec3) -> Vec<Real>;

    /// Length along the member axis in metres.
    fn length(&self) -> Real;

    /// Lumped nodal mass in tonnes, in local DOF order.
    ///
    /// Mass is in tonnes so that `kN / (m·t) = s⁻²` for modal analysis.
    fn lumped_mass(&self) -> &[Real];

    /// Global geometric stiffness for the current initial-stress state.
    ///
    /// `None` means the formulation has no initial-stress matrix, so a
    /// second-order analysis leaves it unchanged.
    fn geometric_stiffness(&self, _stress: InitialStress) -> Option<Vec<Real>> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Section, SectionId};

    /// Relative comparison for values the two sides compute differently.
    fn close(left: Real, right: Real) -> bool {
        (left - right).abs() <= 1.0e-12 * left.abs().max(right.abs()).max(1.0)
    }

    #[test]
    fn the_member_mass_split_is_shared_between_frame_and_truss() {
        let section = Section::new(
            SectionId::new(1),
            "reference",
            0.01,
            2.0e-5,
            3.0e-5,
            4.0e-5,
            0.0,
        );
        let density = 7850.0;
        let length = 2.0;
        let frame = member_lumped_mass(density, &section, length, true);
        let truss = member_lumped_mass(density, &section, length, false);

        let translational = density * 0.01 * length / 2000.0;
        for component in 0..3 {
            assert!(close(frame[component], translational));
            assert!(close(frame[component + DOF_PER_NODE], translational));
            assert!(close(truss[component], translational));
            assert!(close(truss[component + DOF_PER_NODE], translational));
        }
        assert!(frame[3] > 0.0 && frame[4] > 0.0 && frame[5] > 0.0);
        assert!(close(truss[3], 0.0));
        assert!(close(truss[4], 0.0));
        assert!(close(truss[5], 0.0));
        assert!(close(truss[DOF_PER_NODE + 3], 0.0));
    }

    struct DummyElement {
        dofs: [usize; 2 * DOF_PER_NODE],
        stiffness: Vec<Real>,
        axes: Mat3,
        mass: Vec<Real>,
    }

    impl DummyElement {
        fn new() -> Self {
            let size = 2 * DOF_PER_NODE;
            Self {
                dofs: member_dofs(0, 1),
                stiffness: (0..size)
                    .flat_map(|row| {
                        (0..size).map(move |column| {
                            if row == column {
                                (row + 1) as Real
                            } else {
                                0.0
                            }
                        })
                    })
                    .collect(),
                axes: Mat3::from_rows([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
                mass: vec![0.5; size],
            }
        }
    }

    impl Element for DummyElement {
        fn dofs(&self) -> &[usize] {
            &self.dofs
        }

        fn stiffness(&self) -> &[Real] {
            &self.stiffness
        }

        fn axes(&self) -> &Mat3 {
            &self.axes
        }

        fn equivalent_load(&self, q: Vec3) -> Vec<Real> {
            let mut load = vec![0.0; self.dofs.len()];
            load[0] = q.x;
            load[1] = q.y;
            load[2] = q.z;
            load[DOF_PER_NODE] = q.x;
            load[DOF_PER_NODE + 1] = q.y;
            load[DOF_PER_NODE + 2] = q.z;
            load
        }

        fn end_forces(&self, d: &[Real], q: Vec3) -> Vec<Real> {
            let load = self.equivalent_load(q);
            (0..d.len())
                .map(|row| {
                    let stress: Real = (0..d.len())
                        .map(|column| self.stiffness[row * d.len() + column] * d[column])
                        .sum();
                    stress - load[row]
                })
                .collect()
        }

        fn length(&self) -> Real {
            2.0
        }

        fn lumped_mass(&self) -> &[Real] {
            &self.mass
        }
    }

    /// Places one element's coefficients into a global matrix of `size²`.
    fn assemble(element: &dyn Element, size: usize) -> Vec<Real> {
        let mut global = vec![0.0; size * size];
        let dofs = element.dofs();
        let local_size = dofs.len();
        for (local_row, &row) in dofs.iter().enumerate() {
            for (local_column, &column) in dofs.iter().enumerate() {
                global[row * size + column] +=
                    element.stiffness()[local_row * local_size + local_column];
            }
        }
        global
    }

    #[test]
    fn dof_helpers_use_one_order_everywhere() {
        assert_eq!(DOF_NAMES, ["tx", "ty", "tz", "rx", "ry", "rz"]);
        assert_eq!(node_dofs(0), [0, 1, 2, 3, 4, 5]);
        assert_eq!(node_dofs(2), [12, 13, 14, 15, 16, 17]);
        assert_eq!(member_dofs(0, 1), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
        assert_eq!(member_dofs(3, 4)[0], 3 * DOF_PER_NODE);
        assert_eq!(member_dofs(3, 4)[DOF_PER_NODE], 4 * DOF_PER_NODE);
    }

    #[test]
    fn a_dummy_element_assembles_through_the_trait() {
        let element = DummyElement::new();
        let formulation: &dyn Element = &element;
        let size = 2 * DOF_PER_NODE;
        let global = assemble(formulation, size);
        let local_size = formulation.dofs().len();
        for (local_row, &row) in formulation.dofs().iter().enumerate() {
            for (local_column, &column) in formulation.dofs().iter().enumerate() {
                let expected = formulation.stiffness()[local_row * local_size + local_column];
                assert_eq!(global[row * size + column], expected);
            }
        }
        for row in 0..size {
            for column in 0..size {
                assert_eq!(global[row * size + column], global[column * size + row]);
            }
        }
    }

    #[test]
    fn the_load_and_mass_hooks_are_reachable_through_the_trait() {
        let element = DummyElement::new();
        let formulation: &dyn Element = &element;
        let q = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(formulation.length(), 2.0);
        assert_eq!(
            formulation.axes(),
            &Mat3::from_rows([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        );
        assert_eq!(formulation.equivalent_load(q)[0], 1.0);
        assert_eq!(formulation.equivalent_load(q)[DOF_PER_NODE + 2], 3.0);
        assert_eq!(formulation.lumped_mass().len(), 2 * DOF_PER_NODE);
        assert!(formulation
            .geometric_stiffness(InitialStress::axial(-1.0))
            .is_none());
        let d = vec![0.0; 2 * DOF_PER_NODE];
        let forces = formulation.end_forces(&d, q);
        assert_eq!(forces[0], -1.0);
        assert_eq!(forces.len(), formulation.dofs().len());
    }
}
