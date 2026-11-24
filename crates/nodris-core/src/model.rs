//! Neutral model entities with stable ids.
//!
//! A model holds what an editor can change: nodes, materials and sections.
//! None of these is an analysis element yet; elements are derived later from
//! the same objects. Ids are `u32` values owned by the document, so a selection
//! or a reference survives a reorder of the vectors.
//!
//! Coordinates and section properties are stored in SI. A section table that
//! quotes mm, mm², mm⁴ or mm⁶ crosses [`crate::units`] once in `Section::from_tables`
//! and never again.

use std::collections::HashMap;
use std::fmt;

use crate::math::Real;
use crate::units::{area_from_mm2, inertia_from_mm4, stress_from_mpa, warping_from_mm6};

/// Stable identity of a node in the model document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeId(u32);

impl NodeId {
    /// Wraps a document id.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// The underlying document id.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// Stable identity of a material in the model document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MaterialId(u32);

impl MaterialId {
    /// Wraps a document id.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// The underlying document id.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// Stable identity of a section in the model document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SectionId(u32);

impl SectionId {
    /// Wraps a document id.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// The underlying document id.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// A point in model space, in metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Node {
    /// Stable document id.
    pub id: NodeId,
    /// X coordinate in metres.
    pub x: Real,
    /// Y coordinate in metres.
    pub y: Real,
    /// Z coordinate in metres.
    pub z: Real,
}

impl Node {
    /// Creates a node from its id and SI coordinates.
    pub const fn new(id: NodeId, x: Real, y: Real, z: Real) -> Self {
        Self { id, x, y, z }
    }
}

/// Material properties in analysis units.
#[derive(Clone, Debug, PartialEq)]
pub struct Material {
    /// Stable document id.
    pub id: MaterialId,
    /// Human-facing name from the material table.
    pub name: String,
    /// Mass density in kg/m³.
    pub density: Real,
    /// Elastic modulus in kN/m².
    pub elastic_modulus: Real,
    /// Poisson ratio, dimensionless.
    pub poisson_ratio: Real,
}

impl Material {
    /// Creates a material from SI values.
    pub fn new(
        id: MaterialId,
        name: impl Into<String>,
        density: Real,
        elastic_modulus: Real,
        poisson_ratio: Real,
    ) -> Self {
        Self {
            id,
            name: name.into(),
            density,
            elastic_modulus,
            poisson_ratio,
        }
    }

    /// Creates a material from a modulus quoted in MPa.
    ///
    /// This is the single conversion point for a table modulus. The stored
    /// value stays in kN/m² for the rest of the kernel.
    pub fn from_mpa(
        id: MaterialId,
        name: impl Into<String>,
        density: Real,
        elastic_modulus_mpa: Real,
        poisson_ratio: Real,
    ) -> Self {
        Self::new(
            id,
            name,
            density,
            stress_from_mpa(elastic_modulus_mpa),
            poisson_ratio,
        )
    }
}

/// Cross-section properties in SI.
#[derive(Clone, Debug, PartialEq)]
pub struct Section {
    /// Stable document id.
    pub id: SectionId,
    /// Human-facing name from the section table.
    pub name: String,
    /// Area in m².
    pub area: Real,
    /// Second moment about the local y axis, in m⁴.
    pub inertia_y: Real,
    /// Second moment about the local z axis, in m⁴.
    pub inertia_z: Real,
    /// Torsional constant in m⁴.
    pub torsion_constant: Real,
    /// Warping constant in m⁶.
    pub warping_constant: Real,
}

impl Section {
    /// Creates a section from SI values.
    pub fn new(
        id: SectionId,
        name: impl Into<String>,
        area: Real,
        inertia_y: Real,
        inertia_z: Real,
        torsion_constant: Real,
        warping_constant: Real,
    ) -> Self {
        Self {
            id,
            name: name.into(),
            area,
            inertia_y,
            inertia_z,
            torsion_constant,
            warping_constant,
        }
    }

    /// Creates a section from the units used by rolled-section tables.
    ///
    /// Areas arrive in mm², moments and torsional constants in mm⁴, and the
    /// warping constant in mm⁶. Each value is converted here, once, so an
    /// element never has to know which table it came from.
    pub fn from_tables(
        id: SectionId,
        name: impl Into<String>,
        area_mm2: Real,
        inertia_y_mm4: Real,
        inertia_z_mm4: Real,
        torsion_mm4: Real,
        warping_mm6: Real,
    ) -> Self {
        Self::new(
            id,
            name,
            area_from_mm2(area_mm2),
            inertia_from_mm4(inertia_y_mm4),
            inertia_from_mm4(inertia_z_mm4),
            inertia_from_mm4(torsion_mm4),
            warping_from_mm6(warping_mm6),
        )
    }
}

/// The first model entities, before members, loads or analysis.
#[derive(Clone, Debug, Default)]
pub struct Model {
    nodes: Vec<Node>,
    materials: Vec<Material>,
    sections: Vec<Section>,
}

impl Model {
    /// Creates an empty model.
    pub fn new() -> Self {
        Self::default()
    }

    /// All nodes in insertion order.
    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    /// All materials in insertion order.
    pub fn materials(&self) -> &[Material] {
        &self.materials
    }

    /// All sections in insertion order.
    pub fn sections(&self) -> &[Section] {
        &self.sections
    }

    /// Adds a node, refusing an id that is already in the model.
    pub fn add_node(&mut self, node: Node) -> Result<(), ModelError> {
        if self.nodes.iter().any(|existing| existing.id == node.id) {
            return Err(ModelError::DuplicateNode { id: node.id });
        }
        self.nodes.push(node);
        Ok(())
    }

    /// Adds a material, refusing an id that is already in the model.
    pub fn add_material(&mut self, material: Material) -> Result<(), ModelError> {
        if self
            .materials
            .iter()
            .any(|existing| existing.id == material.id)
        {
            return Err(ModelError::DuplicateMaterial { id: material.id });
        }
        self.materials.push(material);
        Ok(())
    }

    /// Adds a section, refusing an id that is already in the model.
    pub fn add_section(&mut self, section: Section) -> Result<(), ModelError> {
        if self
            .sections
            .iter()
            .any(|existing| existing.id == section.id)
        {
            return Err(ModelError::DuplicateSection { id: section.id });
        }
        self.sections.push(section);
        Ok(())
    }

    /// Builds the id-to-position lookup for the current contents.
    pub fn index(&self) -> ModelIndex {
        ModelIndex::build(self)
    }
}

/// A model that cannot contain two entities of one kind with the same id.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ModelError {
    /// A node id was used twice.
    DuplicateNode { id: NodeId },
    /// A material id was used twice.
    DuplicateMaterial { id: MaterialId },
    /// A section id was used twice.
    DuplicateSection { id: SectionId },
}

impl fmt::Display for ModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ModelError::DuplicateNode { id } => write!(f, "duplicate node id {}", id.get()),
            ModelError::DuplicateMaterial { id } => {
                write!(f, "duplicate material id {}", id.get())
            }
            ModelError::DuplicateSection { id } => {
                write!(f, "duplicate section id {}", id.get())
            }
        }
    }
}

impl std::error::Error for ModelError {}

/// Maps stable ids to positions in a [`Model`].
///
/// Lookups return vector positions, not references, so the caller decides
/// whether it needs the entity itself or a mutable handle into the model.
#[derive(Clone, Debug, Default)]
pub struct ModelIndex {
    nodes: HashMap<NodeId, usize>,
    materials: HashMap<MaterialId, usize>,
    sections: HashMap<SectionId, usize>,
}

impl ModelIndex {
    /// Builds an index from a model whose ids are already unique.
    pub fn build(model: &Model) -> Self {
        let mut index = Self::default();
        for (position, node) in model.nodes().iter().enumerate() {
            let inserted = index.nodes.insert(node.id, position);
            debug_assert!(inserted.is_none(), "duplicate node id in model");
        }
        for (position, material) in model.materials().iter().enumerate() {
            let inserted = index.materials.insert(material.id, position);
            debug_assert!(inserted.is_none(), "duplicate material id in model");
        }
        for (position, section) in model.sections().iter().enumerate() {
            let inserted = index.sections.insert(section.id, position);
            debug_assert!(inserted.is_none(), "duplicate section id in model");
        }
        index
    }

    /// Position of a node id in `Model::nodes`.
    pub fn node(&self, id: NodeId) -> Option<usize> {
        self.nodes.get(&id).copied()
    }

    /// Position of a material id in `Model::materials`.
    pub fn material(&self, id: MaterialId) -> Option<usize> {
        self.materials.get(&id).copied()
    }

    /// Position of a section id in `Model::sections`.
    pub fn section(&self, id: SectionId) -> Option<usize> {
        self.sections.get(&id).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(got: Real, want: Real) -> bool {
        (got - want).abs() <= 1.0e-12 * want.abs().max(1.0)
    }

    fn sample_model() -> Model {
        let mut model = Model::new();
        model
            .add_node(Node::new(NodeId::new(7), 1.0, 2.0, 3.0))
            .unwrap();
        model
            .add_node(Node::new(NodeId::new(3), 0.0, 0.0, 0.0))
            .unwrap();
        model
            .add_material(Material::from_mpa(
                MaterialId::new(2),
                "S355",
                7850.0,
                210_000.0,
                0.3,
            ))
            .unwrap();
        model
            .add_section(Section::from_tables(
                SectionId::new(9),
                "test",
                1.0e6,
                1.0e12,
                5.0e11,
                2.0e11,
                3.0e18,
            ))
            .unwrap();
        model
    }

    #[test]
    fn entity_ids_round_trip_through_their_wrappers() {
        assert_eq!(NodeId::new(11).get(), 11);
        assert_eq!(MaterialId::new(12).get(), 12);
        assert_eq!(SectionId::new(13).get(), 13);
    }

    #[test]
    fn index_finds_each_entity_by_id() {
        let model = sample_model();
        let index = model.index();

        assert_eq!(index.node(NodeId::new(7)), Some(0));
        assert_eq!(index.node(NodeId::new(3)), Some(1));
        assert_eq!(index.material(MaterialId::new(2)), Some(0));
        assert_eq!(index.section(SectionId::new(9)), Some(0));
    }

    #[test]
    fn index_returns_none_for_an_unknown_id() {
        let index = sample_model().index();
        assert_eq!(index.node(NodeId::new(99)), None);
        assert_eq!(index.material(MaterialId::new(99)), None);
        assert_eq!(index.section(SectionId::new(99)), None);
    }

    #[test]
    fn duplicate_ids_are_rejected_at_insertion() {
        let mut model = sample_model();

        assert_eq!(
            model.add_node(Node::new(NodeId::new(7), 4.0, 5.0, 6.0)),
            Err(ModelError::DuplicateNode { id: NodeId::new(7) })
        );
        assert_eq!(
            model.add_material(Material::new(
                MaterialId::new(2),
                "duplicate",
                1.0,
                1.0,
                0.0
            )),
            Err(ModelError::DuplicateMaterial {
                id: MaterialId::new(2)
            })
        );
        assert_eq!(
            model.add_section(Section::new(
                SectionId::new(9),
                "duplicate",
                1.0,
                1.0,
                1.0,
                1.0,
                1.0
            )),
            Err(ModelError::DuplicateSection {
                id: SectionId::new(9)
            })
        );
    }

    #[test]
    fn duplicate_errors_name_the_id_and_entity() {
        let error = ModelError::DuplicateNode {
            id: NodeId::new(42),
        };
        assert_eq!(error.to_string(), "duplicate node id 42");
    }

    #[test]
    fn section_table_values_convert_once() {
        let section = sample_model().sections()[0].clone();
        assert!(close(section.area, 1.0));
        assert!(close(section.inertia_y, 1.0));
        assert!(close(section.inertia_z, 5.0e-1));
        assert!(close(section.torsion_constant, 2.0e-1));
        assert!(close(section.warping_constant, 3.0));
    }

    #[test]
    fn material_modulus_converts_from_mpa() {
        let model = sample_model();
        let material = &model.materials()[0];
        assert!(close(material.elastic_modulus, 210_000_000.0));
        assert!(close(material.density, 7850.0));
        assert_eq!(material.name, "S355");
    }
}
