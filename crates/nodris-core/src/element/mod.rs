//! Element formulations.
//!
//! Every formulation implements [`Element`], the contract the analysis
//! pipeline sees. An assembler reads the DOF numbers, the stiffness and the
//! load hook without branching on the formulation.
//!
//! There are six degrees of freedom per node, always in this order:
//! translation X, Y, Z, then rotation about X, Y, Z. [`node_dofs`] is the one
//! place that turns a node index into those global numbers.

use std::fmt;

use crate::math::{Mat3, Real, Vec3};

mod frame3d;

pub use frame3d::local_frame_stiffness;

/// Invalid input to an element formulation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ElementError {
    /// A member length must be finite and greater than zero.
    InvalidLength { length: Real },
    /// A rigidity must be finite and greater than zero.
    InvalidRigidity { name: &'static str, value: Real },
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
        }
    }
}

impl std::error::Error for ElementError {}

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
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let d = vec![0.0; 2 * DOF_PER_NODE];
        let forces = formulation.end_forces(&d, q);
        assert_eq!(forces[0], -1.0);
        assert_eq!(forces.len(), formulation.dofs().len());
    }
}
