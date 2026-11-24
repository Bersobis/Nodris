//! Neutral model entities with stable ids.
//!
//! A model holds what an editor can change: nodes, materials, sections,
//! members, areas and spline patches. None of these is an analysis element yet;
//! elements are derived later from the same objects. Ids are `u32` values owned
//! by the document, so a selection or a reference survives a reorder of the
//! vectors.
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

/// Stable identity of a member in the model document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MemberId(u32);

impl MemberId {
    /// Wraps a document id.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// The underlying document id.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// What a member is meant to do in the structure.
///
/// This is user intent, not a formulation. Analysis chooses an element
/// implementation later from the same member and section.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElementKind {
    /// A member that carries axial force, shear, bending and torsion.
    Frame,
    /// A member that carries axial force only.
    Truss,
}

impl ElementKind {
    /// Stable spelling used by documents and reports.
    pub const fn as_str(self) -> &'static str {
        match self {
            ElementKind::Frame => "frame",
            ElementKind::Truss => "truss",
        }
    }

    /// Parses a document spelling, case-insensitively.
    pub fn parse(text: &str) -> Option<Self> {
        match text.trim().to_ascii_lowercase().as_str() {
            "frame" => Some(ElementKind::Frame),
            "truss" => Some(ElementKind::Truss),
            _ => None,
        }
    }
}

/// A line between two nodes, with the section and kind chosen by the user.
///
/// A member is intent. It holds no stiffness, no DOFs and no mesh; analysis
/// derives those from it without rewriting the document.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Member {
    /// Stable document id.
    pub id: MemberId,
    /// Node at the start of the member.
    pub start: NodeId,
    /// Node at the end of the member.
    pub end: NodeId,
    /// Section assigned to the member.
    pub section: SectionId,
    /// Rotation of the local axes about the member, in radians.
    pub rotation: Real,
    /// Frame or truss intent.
    pub kind: ElementKind,
}

impl Member {
    /// Creates a member from stable ids and analysis-unit rotation.
    pub const fn new(
        id: MemberId,
        start: NodeId,
        end: NodeId,
        section: SectionId,
        rotation: Real,
        kind: ElementKind,
    ) -> Self {
        Self {
            id,
            start,
            end,
            section,
            rotation,
            kind,
        }
    }
}

/// Stable identity of an area in the model document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AreaId(u32);

impl AreaId {
    /// Wraps a document id.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// The underlying document id.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// Stable identity of a spline patch in the model document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PatchId(u32);

impl PatchId {
    /// Wraps a document id.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// The underlying document id.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// What a four-node area is meant to model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AreaKind {
    /// A full shell with membrane and bending behaviour.
    Shell,
    /// In-plane stiffness only.
    Membrane,
    /// Out-of-plane bending only.
    Plate,
}

impl AreaKind {
    /// Stable spelling used by documents and reports.
    pub const fn as_str(self) -> &'static str {
        match self {
            AreaKind::Shell => "shell",
            AreaKind::Membrane => "membrane",
            AreaKind::Plate => "plate",
        }
    }

    /// Parses a document spelling, case-insensitively.
    pub fn parse(text: &str) -> Option<Self> {
        match text.trim().to_ascii_lowercase().as_str() {
            "shell" => Some(AreaKind::Shell),
            "membrane" => Some(AreaKind::Membrane),
            "plate" => Some(AreaKind::Plate),
            _ => None,
        }
    }
}

/// A four-node area prepared for shell meshing.
///
/// The nodes are references to model nodes, so moving a node rebuilds the
/// area geometry without changing the document object.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Area {
    /// Stable document id.
    pub id: AreaId,
    /// Corner nodes in circuit order.
    pub nodes: [NodeId; 4],
    /// Shell, membrane or plate intent.
    pub kind: AreaKind,
}

impl Area {
    /// Creates an area from four corner node ids.
    pub const fn new(id: AreaId, nodes: [NodeId; 4], kind: AreaKind) -> Self {
        Self { id, nodes, kind }
    }
}

/// The parametric direction a knot vector belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PatchDirection {
    /// The direction of the control-net rows.
    U,
    /// The direction of the control-net columns.
    V,
}

impl PatchDirection {
    /// Stable single-letter spelling used in error messages.
    pub const fn as_str(self) -> &'static str {
        match self {
            PatchDirection::U => "u",
            PatchDirection::V => "v",
        }
    }
}

/// An IGA patch defined by control-net node ids and open knot vectors.
///
/// The patch stores no evaluated points or basis values. Analysis rebuilds
/// them from the knots, the control net and the weights when it meshes.
#[derive(Clone, Debug, PartialEq)]
pub struct SplinePatch {
    /// Stable document id.
    pub id: PatchId,
    /// Polynomial degree in the U direction.
    pub degree_u: u8,
    /// Polynomial degree in the V direction.
    pub degree_v: u8,
    /// Open knot vector for the control-net rows.
    pub knots_u: Vec<Real>,
    /// Open knot vector for the control-net columns.
    pub knots_v: Vec<Real>,
    /// Control-net node ids, one row per U knot span.
    pub control_net: Vec<Vec<NodeId>>,
    /// NURBS weights in the same shape as the control net.
    pub weights: Vec<Vec<Real>>,
}

impl SplinePatch {
    /// Creates a patch and checks its knot vectors and control-net shape.
    pub fn new(
        id: PatchId,
        degree_u: u8,
        degree_v: u8,
        knots_u: Vec<Real>,
        knots_v: Vec<Real>,
        control_net: Vec<Vec<NodeId>>,
        weights: Vec<Vec<Real>>,
    ) -> Self {
        Self {
            id,
            degree_u,
            degree_v,
            knots_u,
            knots_v,
            control_net,
            weights,
        }
    }

    /// Creates a non-rational patch with every weight set to one.
    pub fn unweighted(
        id: PatchId,
        degree_u: u8,
        degree_v: u8,
        knots_u: Vec<Real>,
        knots_v: Vec<Real>,
        control_net: Vec<Vec<NodeId>>,
    ) -> Self {
        let weights = control_net.iter().map(|row| vec![1.0; row.len()]).collect();
        Self::new(
            id,
            degree_u,
            degree_v,
            knots_u,
            knots_v,
            control_net,
            weights,
        )
    }
}

fn check_knots(
    patch: PatchId,
    degree: u8,
    knots: &[Real],
    direction: PatchDirection,
) -> Result<(), ModelError> {
    if degree == 0 {
        return Err(ModelError::PatchDegree { patch, direction });
    }
    if knots.iter().any(|knot| !knot.is_finite())
        || knots.windows(2).any(|pair| pair[1] < pair[0])
        || knots.last() <= knots.first()
    {
        return Err(ModelError::PatchKnotsInvalid { patch, direction });
    }
    let p = degree as usize;
    let last = knots.len() - 1;
    let clamped_start = knots[..=p].iter().all(|knot| *knot == knots[0]);
    let clamped_end = knots[last - p..].iter().all(|knot| *knot == knots[last]);
    if !clamped_start || !clamped_end {
        return Err(ModelError::PatchKnotsNotClamped { patch, direction });
    }
    Ok(())
}

fn validate_patch_geometry(patch: &SplinePatch) -> Result<(), ModelError> {
    if patch.control_net.is_empty() || patch.control_net[0].is_empty() {
        return Err(ModelError::PatchNetShape { patch: patch.id });
    }
    let rows = patch.control_net.len();
    let columns = patch.control_net[0].len();
    if patch.control_net.iter().any(|row| row.len() != columns)
        || patch.weights.len() != rows
        || patch.weights.iter().any(|row| row.len() != columns)
    {
        return Err(ModelError::PatchNetShape { patch: patch.id });
    }
    if rows < patch.degree_u as usize + 1 || columns < patch.degree_v as usize + 1 {
        return Err(ModelError::PatchNetShape { patch: patch.id });
    }
    if patch.knots_u.len() != rows + patch.degree_u as usize + 1 {
        return Err(ModelError::PatchKnotCount {
            patch: patch.id,
            direction: PatchDirection::U,
        });
    }
    if patch.knots_v.len() != columns + patch.degree_v as usize + 1 {
        return Err(ModelError::PatchKnotCount {
            patch: patch.id,
            direction: PatchDirection::V,
        });
    }
    check_knots(patch.id, patch.degree_u, &patch.knots_u, PatchDirection::U)?;
    check_knots(patch.id, patch.degree_v, &patch.knots_v, PatchDirection::V)
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
    members: Vec<Member>,
    areas: Vec<Area>,
    patches: Vec<SplinePatch>,
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

    /// All members in insertion order.
    pub fn members(&self) -> &[Member] {
        &self.members
    }

    /// All areas in insertion order.
    pub fn areas(&self) -> &[Area] {
        &self.areas
    }

    /// All spline patches in insertion order.
    pub fn patches(&self) -> &[SplinePatch] {
        &self.patches
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

    /// Adds a member, refusing a duplicate id or a missing reference.
    pub fn add_member(&mut self, member: Member) -> Result<(), ModelError> {
        if self.members.iter().any(|existing| existing.id == member.id) {
            return Err(ModelError::DuplicateMember { id: member.id });
        }
        if !self.nodes.iter().any(|node| node.id == member.start) {
            return Err(ModelError::UnknownNode {
                member: member.id,
                node: member.start,
            });
        }
        if !self.nodes.iter().any(|node| node.id == member.end) {
            return Err(ModelError::UnknownNode {
                member: member.id,
                node: member.end,
            });
        }
        if !self
            .sections
            .iter()
            .any(|section| section.id == member.section)
        {
            return Err(ModelError::UnknownSection {
                member: member.id,
                section: member.section,
            });
        }
        self.members.push(member);
        Ok(())
    }

    /// Adds an area, refusing a duplicate id or a missing corner node.
    pub fn add_area(&mut self, area: Area) -> Result<(), ModelError> {
        if self.areas.iter().any(|existing| existing.id == area.id) {
            return Err(ModelError::DuplicateArea { id: area.id });
        }
        for node in area.nodes {
            if !self.nodes.iter().any(|existing| existing.id == node) {
                return Err(ModelError::UnknownAreaNode {
                    area: area.id,
                    node,
                });
            }
        }
        self.areas.push(area);
        Ok(())
    }

    /// Adds a patch, refusing a duplicate id, a bad net or a missing node.
    pub fn add_patch(&mut self, patch: SplinePatch) -> Result<(), ModelError> {
        if self.patches.iter().any(|existing| existing.id == patch.id) {
            return Err(ModelError::DuplicatePatch { id: patch.id });
        }
        validate_patch_geometry(&patch)?;
        for row in &patch.control_net {
            for node in row {
                if !self.nodes.iter().any(|existing| existing.id == *node) {
                    return Err(ModelError::UnknownPatchNode {
                        patch: patch.id,
                        node: *node,
                    });
                }
            }
        }
        self.patches.push(patch);
        Ok(())
    }

    /// Checks that every member references nodes and a section in the model.
    pub fn validate_references(&self) -> Result<(), ModelError> {
        let index = self.index();
        for member in &self.members {
            if index.node(member.start).is_none() {
                return Err(ModelError::UnknownNode {
                    member: member.id,
                    node: member.start,
                });
            }
            if index.node(member.end).is_none() {
                return Err(ModelError::UnknownNode {
                    member: member.id,
                    node: member.end,
                });
            }
            if index.section(member.section).is_none() {
                return Err(ModelError::UnknownSection {
                    member: member.id,
                    section: member.section,
                });
            }
        }
        for area in &self.areas {
            for node in area.nodes {
                if index.node(node).is_none() {
                    return Err(ModelError::UnknownAreaNode {
                        area: area.id,
                        node,
                    });
                }
            }
        }
        for patch in &self.patches {
            for row in &patch.control_net {
                for node in row {
                    if index.node(*node).is_none() {
                        return Err(ModelError::UnknownPatchNode {
                            patch: patch.id,
                            node: *node,
                        });
                    }
                }
            }
        }
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
    /// A member id was used twice.
    DuplicateMember { id: MemberId },
    /// A member references a node that is not in the model.
    UnknownNode { member: MemberId, node: NodeId },
    /// A member references a section that is not in the model.
    UnknownSection {
        member: MemberId,
        section: SectionId,
    },
    /// An area id was used twice.
    DuplicateArea { id: AreaId },
    /// A spline patch id was used twice.
    DuplicatePatch { id: PatchId },
    /// An area references a node that is not in the model.
    UnknownAreaNode { area: AreaId, node: NodeId },
    /// A patch references a node that is not in the model.
    UnknownPatchNode { patch: PatchId, node: NodeId },
    /// A patch control net or weight net is not rectangular.
    PatchNetShape { patch: PatchId },
    /// A knot vector does not match its control-net dimension.
    PatchKnotCount {
        patch: PatchId,
        direction: PatchDirection,
    },
    /// A knot vector holds a non-finite or decreasing value.
    PatchKnotsInvalid {
        patch: PatchId,
        direction: PatchDirection,
    },
    /// A knot vector is not open at both ends.
    PatchKnotsNotClamped {
        patch: PatchId,
        direction: PatchDirection,
    },
    /// A patch degree is zero.
    PatchDegree {
        patch: PatchId,
        direction: PatchDirection,
    },
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
            ModelError::DuplicateMember { id } => write!(f, "duplicate member id {}", id.get()),
            ModelError::UnknownNode { member, node } => write!(
                f,
                "member {} references unknown node {}",
                member.get(),
                node.get()
            ),
            ModelError::UnknownSection { member, section } => write!(
                f,
                "member {} references unknown section {}",
                member.get(),
                section.get()
            ),
            ModelError::DuplicateArea { id } => write!(f, "duplicate area id {}", id.get()),
            ModelError::DuplicatePatch { id } => write!(f, "duplicate patch id {}", id.get()),
            ModelError::UnknownAreaNode { area, node } => write!(
                f,
                "area {} references unknown node {}",
                area.get(),
                node.get()
            ),
            ModelError::UnknownPatchNode { patch, node } => write!(
                f,
                "patch {} references unknown node {}",
                patch.get(),
                node.get()
            ),
            ModelError::PatchNetShape { patch } => {
                write!(f, "patch {} has a non-rectangular net", patch.get())
            }
            ModelError::PatchKnotCount { patch, direction } => write!(
                f,
                "patch {} has a wrong {} knot count",
                patch.get(),
                direction.as_str()
            ),
            ModelError::PatchKnotsInvalid { patch, direction } => write!(
                f,
                "patch {} has an invalid {} knot vector",
                patch.get(),
                direction.as_str()
            ),
            ModelError::PatchKnotsNotClamped { patch, direction } => write!(
                f,
                "patch {} has an unclamped {} knot vector",
                patch.get(),
                direction.as_str()
            ),
            ModelError::PatchDegree { patch, direction } => write!(
                f,
                "patch {} has zero degree in {}",
                patch.get(),
                direction.as_str()
            ),
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
    members: HashMap<MemberId, usize>,
    areas: HashMap<AreaId, usize>,
    patches: HashMap<PatchId, usize>,
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
        for (position, member) in model.members().iter().enumerate() {
            let inserted = index.members.insert(member.id, position);
            debug_assert!(inserted.is_none(), "duplicate member id in model");
        }
        for (position, area) in model.areas().iter().enumerate() {
            let inserted = index.areas.insert(area.id, position);
            debug_assert!(inserted.is_none(), "duplicate area id in model");
        }
        for (position, patch) in model.patches().iter().enumerate() {
            let inserted = index.patches.insert(patch.id, position);
            debug_assert!(inserted.is_none(), "duplicate patch id in model");
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

    /// Position of a member id in `Model::members`.
    pub fn member(&self, id: MemberId) -> Option<usize> {
        self.members.get(&id).copied()
    }

    /// Position of an area id in `Model::areas`.
    pub fn area(&self, id: AreaId) -> Option<usize> {
        self.areas.get(&id).copied()
    }

    /// Position of a patch id in `Model::patches`.
    pub fn patch(&self, id: PatchId) -> Option<usize> {
        self.patches.get(&id).copied()
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
            .add_node(Node::new(NodeId::new(11), 2.0, 0.0, 0.0))
            .unwrap();
        model
            .add_node(Node::new(NodeId::new(12), 3.0, 0.0, 0.0))
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
            .add_member(Member::new(
                MemberId::new(4),
                NodeId::new(7),
                NodeId::new(3),
                SectionId::new(9),
                0.0,
                ElementKind::Frame,
            ))
            .unwrap();
        model
            .add_area(Area::new(
                AreaId::new(5),
                [
                    NodeId::new(7),
                    NodeId::new(3),
                    NodeId::new(11),
                    NodeId::new(12),
                ],
                AreaKind::Shell,
            ))
            .unwrap();
        model
            .add_patch(SplinePatch::unweighted(
                PatchId::new(6),
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
                vec![
                    vec![NodeId::new(7), NodeId::new(3)],
                    vec![NodeId::new(11), NodeId::new(12)],
                ],
            ))
            .unwrap();
        model
    }

    #[test]
    fn entity_ids_round_trip_through_their_wrappers() {
        assert_eq!(NodeId::new(11).get(), 11);
        assert_eq!(MaterialId::new(12).get(), 12);
        assert_eq!(SectionId::new(13).get(), 13);
        assert_eq!(MemberId::new(14).get(), 14);
        assert_eq!(AreaId::new(15).get(), 15);
        assert_eq!(PatchId::new(16).get(), 16);
    }

    #[test]
    fn index_finds_each_entity_by_id() {
        let model = sample_model();
        let index = model.index();

        assert_eq!(index.node(NodeId::new(7)), Some(0));
        assert_eq!(index.node(NodeId::new(3)), Some(1));
        assert_eq!(index.material(MaterialId::new(2)), Some(0));
        assert_eq!(index.section(SectionId::new(9)), Some(0));
        assert_eq!(index.member(MemberId::new(4)), Some(0));
        assert_eq!(index.area(AreaId::new(5)), Some(0));
        assert_eq!(index.patch(PatchId::new(6)), Some(0));
    }

    #[test]
    fn index_returns_none_for_an_unknown_id() {
        let index = sample_model().index();
        assert_eq!(index.node(NodeId::new(99)), None);
        assert_eq!(index.material(MaterialId::new(99)), None);
        assert_eq!(index.section(SectionId::new(99)), None);
        assert_eq!(index.member(MemberId::new(99)), None);
        assert_eq!(index.area(AreaId::new(99)), None);
        assert_eq!(index.patch(PatchId::new(99)), None);
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
        let duplicate = Member::new(
            MemberId::new(4),
            NodeId::new(7),
            NodeId::new(3),
            SectionId::new(9),
            0.0,
            ElementKind::Truss,
        );
        assert_eq!(
            model.add_member(duplicate),
            Err(ModelError::DuplicateMember {
                id: MemberId::new(4)
            })
        );
        let duplicate_area = Area::new(
            AreaId::new(5),
            [
                NodeId::new(7),
                NodeId::new(3),
                NodeId::new(11),
                NodeId::new(12),
            ],
            AreaKind::Plate,
        );
        assert_eq!(
            model.add_area(duplicate_area),
            Err(ModelError::DuplicateArea { id: AreaId::new(5) })
        );
        let duplicate_patch = SplinePatch::unweighted(
            PatchId::new(6),
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![NodeId::new(7), NodeId::new(3)],
                vec![NodeId::new(11), NodeId::new(12)],
            ],
        );
        assert_eq!(
            model.add_patch(duplicate_patch),
            Err(ModelError::DuplicatePatch {
                id: PatchId::new(6)
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

    #[test]
    fn element_kind_spellings_round_trip() {
        for kind in [ElementKind::Frame, ElementKind::Truss] {
            assert_eq!(ElementKind::parse(kind.as_str()), Some(kind));
        }
        assert_eq!(ElementKind::parse(" FRAME "), Some(ElementKind::Frame));
        assert_eq!(ElementKind::parse("cable"), None);
    }

    #[test]
    fn a_valid_member_passes_reference_validation() {
        assert_eq!(sample_model().validate_references(), Ok(()));
    }

    #[test]
    fn missing_node_reference_is_rejected() {
        let mut model = sample_model();
        let member = Member::new(
            MemberId::new(5),
            NodeId::new(7),
            NodeId::new(99),
            SectionId::new(9),
            0.0,
            ElementKind::Frame,
        );
        assert_eq!(
            model.add_member(member),
            Err(ModelError::UnknownNode {
                member: MemberId::new(5),
                node: NodeId::new(99),
            })
        );
    }

    #[test]
    fn missing_section_reference_is_rejected() {
        let mut model = sample_model();
        let member = Member::new(
            MemberId::new(5),
            NodeId::new(7),
            NodeId::new(3),
            SectionId::new(99),
            0.0,
            ElementKind::Truss,
        );
        assert_eq!(
            model.add_member(member),
            Err(ModelError::UnknownSection {
                member: MemberId::new(5),
                section: SectionId::new(99),
            })
        );
    }

    #[test]
    fn reference_errors_name_the_member_and_target() {
        let error = ModelError::UnknownNode {
            member: MemberId::new(5),
            node: NodeId::new(99),
        };
        assert_eq!(error.to_string(), "member 5 references unknown node 99");
    }

    #[test]
    fn area_kind_spellings_round_trip() {
        for kind in [AreaKind::Shell, AreaKind::Membrane, AreaKind::Plate] {
            assert_eq!(AreaKind::parse(kind.as_str()), Some(kind));
        }
        assert_eq!(AreaKind::parse(" SHELL "), Some(AreaKind::Shell));
        assert_eq!(AreaKind::parse("solid"), None);
    }

    #[test]
    fn area_reference_is_checked() {
        let mut model = sample_model();
        let area = Area::new(
            AreaId::new(8),
            [
                NodeId::new(7),
                NodeId::new(3),
                NodeId::new(11),
                NodeId::new(99),
            ],
            AreaKind::Membrane,
        );
        assert_eq!(
            model.add_area(area),
            Err(ModelError::UnknownAreaNode {
                area: AreaId::new(8),
                node: NodeId::new(99),
            })
        );
    }

    #[test]
    fn patch_reference_is_checked() {
        let mut model = sample_model();
        let patch = SplinePatch::unweighted(
            PatchId::new(8),
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![NodeId::new(7), NodeId::new(99)],
                vec![NodeId::new(11), NodeId::new(12)],
            ],
        );
        assert_eq!(
            model.add_patch(patch),
            Err(ModelError::UnknownPatchNode {
                patch: PatchId::new(8),
                node: NodeId::new(99),
            })
        );
    }

    #[test]
    fn patch_with_ragged_net_is_rejected() {
        let mut model = sample_model();
        let patch = SplinePatch::unweighted(
            PatchId::new(8),
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![vec![NodeId::new(7), NodeId::new(3)], vec![NodeId::new(11)]],
        );
        assert_eq!(
            model.add_patch(patch),
            Err(ModelError::PatchNetShape {
                patch: PatchId::new(8)
            })
        );
    }

    #[test]
    fn patch_with_wrong_knot_count_is_rejected() {
        let mut model = sample_model();
        let patch = SplinePatch::unweighted(
            PatchId::new(8),
            1,
            1,
            vec![0.0, 0.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![NodeId::new(7), NodeId::new(3)],
                vec![NodeId::new(11), NodeId::new(12)],
            ],
        );
        assert_eq!(
            model.add_patch(patch),
            Err(ModelError::PatchKnotCount {
                patch: PatchId::new(8),
                direction: PatchDirection::U,
            })
        );
    }

    #[test]
    fn patch_with_unclamped_knots_is_rejected() {
        let mut model = sample_model();
        let patch = SplinePatch::unweighted(
            PatchId::new(8),
            1,
            1,
            vec![0.0, 0.5, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![NodeId::new(7), NodeId::new(3)],
                vec![NodeId::new(11), NodeId::new(12)],
            ],
        );
        assert_eq!(
            model.add_patch(patch),
            Err(ModelError::PatchKnotsNotClamped {
                patch: PatchId::new(8),
                direction: PatchDirection::U,
            })
        );
    }

    #[test]
    fn patch_weight_shape_is_checked() {
        let patch = SplinePatch::new(
            PatchId::new(8),
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![NodeId::new(7), NodeId::new(3)],
                vec![NodeId::new(11), NodeId::new(12)],
            ],
            vec![vec![1.0, 1.0]],
        );
        assert_eq!(
            validate_patch_geometry(&patch),
            Err(ModelError::PatchNetShape {
                patch: PatchId::new(8)
            })
        );
    }

    #[test]
    fn patch_direction_has_a_stable_spelling() {
        assert_eq!(PatchDirection::U.as_str(), "u");
        assert_eq!(PatchDirection::V.as_str(), "v");
    }
}
