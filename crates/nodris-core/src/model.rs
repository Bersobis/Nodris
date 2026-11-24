//! Neutral model entities with stable ids.
//!
//! A model holds what an editor can change: nodes, materials, sections,
//! members, areas, patches, load cases, loads, combinations and supports.
//! [`Model::validate`] checks blocking structure; [`Model::validate_issues`]
//! adds object references, severity and orphan warnings. None of these is an
//! analysis element yet; elements are derived later from the same objects.
//! Ids are `u32` values owned by the document, so a selection or a reference
//! survives a reorder of the vectors.
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

    fn grid_shape(&self) -> Option<(usize, usize)> {
        let columns = self.control_net.first()?.len();
        if columns == 0 || self.control_net.iter().any(|row| row.len() != columns) {
            return None;
        }
        Some((self.control_net.len(), columns))
    }

    /// Control-point node ids along one patch edge.
    ///
    /// An empty vector means the control net is not rectangular.
    pub fn edge_nodes(&self, edge: PatchEdge) -> Vec<NodeId> {
        let Some((rows, columns)) = self.grid_shape() else {
            return Vec::new();
        };
        match edge {
            PatchEdge::Xi0 => self.control_net[0].clone(),
            PatchEdge::Xi1 => self.control_net[rows - 1].clone(),
            PatchEdge::Eta0 => self.control_net.iter().map(|row| row[0]).collect(),
            PatchEdge::Eta1 => self
                .control_net
                .iter()
                .map(|row| row[columns - 1])
                .collect(),
        }
    }

    /// Reports whether a node is one of the four control-net corners.
    ///
    /// A corner is where neighbouring patches in a grid meet without forming
    /// a seam.
    fn is_corner(&self, node: NodeId) -> bool {
        let Some((rows, columns)) = self.grid_shape() else {
            return false;
        };
        let last_row = rows - 1;
        let last_column = columns - 1;
        [
            (0, 0),
            (0, last_column),
            (last_row, 0),
            (last_row, last_column),
        ]
        .into_iter()
        .any(|(row, column)| self.control_net[row][column] == node)
    }

    /// Knot vector running along one patch edge.
    pub fn edge_knots(&self, edge: PatchEdge) -> &[Real] {
        match edge {
            PatchEdge::Xi0 | PatchEdge::Xi1 => &self.knots_v,
            PatchEdge::Eta0 | PatchEdge::Eta1 => &self.knots_u,
        }
    }
}

/// One parametric edge of a spline patch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PatchEdge {
    /// The edge at `xi = 0`.
    Xi0,
    /// The edge at `xi = 1`.
    Xi1,
    /// The edge at `eta = 0`.
    Eta0,
    /// The edge at `eta = 1`.
    Eta1,
}

impl PatchEdge {
    /// Every edge, useful for seam scans.
    pub const ALL: [PatchEdge; 4] = [
        PatchEdge::Xi0,
        PatchEdge::Xi1,
        PatchEdge::Eta0,
        PatchEdge::Eta1,
    ];

    /// Stable spelling used by reports and the CLI.
    pub const fn as_str(self) -> &'static str {
        match self {
            PatchEdge::Xi0 => "xi0",
            PatchEdge::Xi1 => "xi1",
            PatchEdge::Eta0 => "eta0",
            PatchEdge::Eta1 => "eta1",
        }
    }

    /// Parses a report spelling, case-insensitively.
    pub fn parse(text: &str) -> Option<Self> {
        match text.trim().to_ascii_lowercase().as_str() {
            "xi0" => Some(PatchEdge::Xi0),
            "xi1" => Some(PatchEdge::Xi1),
            "eta0" => Some(PatchEdge::Eta0),
            "eta1" => Some(PatchEdge::Eta1),
            _ => None,
        }
    }
}

/// Stable identity of a patch tie in the model document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PatchTieId(u32);

impl PatchTieId {
    /// Wraps a document id.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// The underlying document id.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// A declared master and slave edge pair that must move together.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PatchTie {
    /// Stable document id.
    pub id: PatchTieId,
    /// Patch supplying the constraint weights.
    pub master: PatchId,
    /// Edge of the master patch.
    pub master_edge: PatchEdge,
    /// Patch whose edge is eliminated.
    pub slave: PatchId,
    /// Edge of the slave patch.
    pub slave_edge: PatchEdge,
}

impl PatchTie {
    /// Creates a patch tie declaration.
    pub const fn new(
        id: PatchTieId,
        master: PatchId,
        master_edge: PatchEdge,
        slave: PatchId,
        slave_edge: PatchEdge,
    ) -> Self {
        Self {
            id,
            master,
            master_edge,
            slave,
            slave_edge,
        }
    }
}

/// Which end of a member was found on a patch edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BeamEnd {
    /// The member's start node.
    Start,
    /// The member's end node.
    End,
}

impl BeamEnd {
    /// Stable spelling used by reports and the CLI.
    pub const fn as_str(self) -> &'static str {
        match self {
            BeamEnd::Start => "start",
            BeamEnd::End => "end",
        }
    }
}

/// A member end detected on a patch edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BeamCoupling {
    /// Member whose end sits on the edge.
    pub member: MemberId,
    /// Which end was found.
    pub end: BeamEnd,
    /// Patch whose edge contains the node.
    pub patch: PatchId,
    /// Edge containing the node.
    pub edge: PatchEdge,
    /// Node at the member end.
    pub node: NodeId,
}

impl BeamCoupling {
    /// Creates a beam coupling record.
    pub const fn new(
        member: MemberId,
        end: BeamEnd,
        patch: PatchId,
        edge: PatchEdge,
        node: NodeId,
    ) -> Self {
        Self {
            member,
            end,
            patch,
            edge,
            node,
        }
    }
}

/// How two patch edges are coupled.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PatchCoupling {
    /// The edges already share every control-point node.
    Shared { node_count: usize },
    /// A declared `PatchTie` joins the two patches.
    Tied { tie: PatchTieId },
    /// The patches share control points that do not form one common edge.
    NeedTie { node: NodeId },
}

/// One exact control-point constraint row produced by a tie.
#[derive(Clone, Debug, PartialEq)]
pub struct TieConstraint {
    /// Slave node eliminated by this constraint.
    pub slave: NodeId,
    /// Master nodes carrying the displacement.
    pub masters: Vec<NodeId>,
    /// Interpolation weights, parallel to `masters`.
    pub weights: Vec<Real>,
}

fn patch_by_id(patches: &[SplinePatch], id: PatchId) -> Option<&SplinePatch> {
    patches.iter().find(|patch| patch.id == id)
}

/// Reports whether a tie joins the given patch pair, in either role.
///
/// The declared edges are not compared here. A tie is the statement that
/// these two patches are joined, and its edges were checked when it was
/// added.
fn tie_covers_patches(tie: &PatchTie, first: PatchId, second: PatchId) -> bool {
    (tie.master == first && tie.slave == second) || (tie.master == second && tie.slave == first)
}

/// The edge of `patch` that holds every node in `shared`.
fn edge_holding(patch: &SplinePatch, shared: &[NodeId]) -> Option<PatchEdge> {
    PatchEdge::ALL.into_iter().find(|edge| {
        let nodes = patch.edge_nodes(*edge);
        !nodes.is_empty() && shared.iter().all(|node| nodes.contains(node))
    })
}

/// Resolves how two patches meet.
///
/// A declared tie wins, because a tie is how two differently refined patches
/// are joined. Otherwise the patches are compared by the control points they
/// have in common: none means they are independent, a lone corner is a joint,
/// and anything else is a seam that must be fully shared.
fn coupling_between(
    patches: &[SplinePatch],
    ties: &[PatchTie],
    first: PatchId,
    second: PatchId,
) -> Option<PatchCoupling> {
    if first == second {
        return None;
    }
    let a = patch_by_id(patches, first)?;
    let b = patch_by_id(patches, second)?;
    if let Some(tie) = ties
        .iter()
        .find(|tie| tie_covers_patches(tie, first, second))
    {
        return Some(PatchCoupling::Tied { tie: tie.id });
    }
    let a_nodes: Vec<NodeId> = a.control_net.iter().flatten().copied().collect();
    let b_nodes: Vec<NodeId> = b.control_net.iter().flatten().copied().collect();
    let mut shared: Vec<NodeId> = Vec::new();
    for node in a_nodes {
        if b_nodes.contains(&node) && !shared.contains(&node) {
            shared.push(node);
        }
    }
    let first_shared = *shared.first()?;
    // Four patches around one corner share that corner and nothing else. That
    // is a joint, not a seam, so it needs no tie.
    if shared.len() == 1 && a.is_corner(first_shared) && b.is_corner(first_shared) {
        return None;
    }
    let (Some(edge_a), Some(edge_b)) = (edge_holding(a, &shared), edge_holding(b, &shared)) else {
        return Some(PatchCoupling::NeedTie { node: first_shared });
    };
    let count = shared.len();
    if count == a.edge_nodes(edge_a).len() && count == b.edge_nodes(edge_b).len() {
        return Some(PatchCoupling::Shared { node_count: count });
    }
    Some(PatchCoupling::NeedTie { node: first_shared })
}

/// Checks one tie against the other declarations.
///
/// `ties` may or may not already contain `tie`. Entries sharing its id are
/// skipped, so the caller can validate a candidate before pushing it.
fn check_tie(patches: &[SplinePatch], ties: &[PatchTie], tie: &PatchTie) -> Result<(), ModelError> {
    if tie.master == tie.slave {
        return Err(ModelError::PatchTieSamePatch {
            tie: tie.id,
            patch: tie.master,
        });
    }
    let master = patch_by_id(patches, tie.master).ok_or(ModelError::UnknownPatchTiePatch {
        tie: tie.id,
        patch: tie.master,
    })?;
    let slave = patch_by_id(patches, tie.slave).ok_or(ModelError::UnknownPatchTiePatch {
        tie: tie.id,
        patch: tie.slave,
    })?;
    let master_nodes = master.edge_nodes(tie.master_edge);
    let slave_nodes = slave.edge_nodes(tie.slave_edge);
    if master_nodes.is_empty()
        || master_nodes.len() != slave_nodes.len()
        || master.edge_knots(tie.master_edge) != slave.edge_knots(tie.slave_edge)
    {
        return Err(ModelError::PatchTieEdgeShape { tie: tie.id });
    }
    // One master may drive several slaves, but a slave node cannot follow two
    // different masters. Overlapping slave edges from different masters would
    // constrain the same node twice.
    for other in ties.iter().filter(|other| other.id != tie.id) {
        if other.master == tie.master {
            continue;
        }
        let other_slave = match patch_by_id(patches, other.slave) {
            Some(patch) => patch,
            None => continue,
        };
        if let Some(node) = slave_nodes
            .iter()
            .find(|node| other_slave.edge_nodes(other.slave_edge).contains(node))
        {
            return Err(ModelError::PatchTieSlaveTwoMasters {
                node: *node,
                first: other.master,
                second: tie.master,
            });
        }
    }
    Ok(())
}

/// Checks every declared tie plus the seams between all patch pairs.
fn check_ties(patches: &[SplinePatch], ties: &[PatchTie]) -> Result<(), ModelError> {
    for (position, tie) in ties.iter().enumerate() {
        if ties[..position]
            .iter()
            .any(|previous| previous.id == tie.id)
        {
            return Err(ModelError::DuplicatePatchTie { id: tie.id });
        }
        check_tie(patches, ties, tie)?;
    }
    for (position, first) in patches.iter().enumerate() {
        for second in &patches[position + 1..] {
            if let Some(PatchCoupling::NeedTie { node }) =
                coupling_between(patches, ties, first.id, second.id)
            {
                return Err(ModelError::TJunctionWithoutTie {
                    first: first.id,
                    second: second.id,
                    node,
                });
            }
        }
    }
    Ok(())
}

/// Stable identity of a load case in the model document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LoadCaseId(u32);

impl LoadCaseId {
    /// Wraps a document id.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// The underlying document id.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// Stable identity of a load combination in the model document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CombinationId(u32);

impl CombinationId {
    /// Wraps a document id.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// The underlying document id.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// A named group of loads that combinations can factor.
#[derive(Clone, Debug, PartialEq)]
pub struct LoadCase {
    /// Stable document id.
    pub id: LoadCaseId,
    /// Human-facing case name.
    pub name: String,
    /// Whether self-weight belongs to this case.
    pub self_weight: bool,
}

impl LoadCase {
    /// Creates a case without self-weight.
    pub fn new(id: LoadCaseId, name: impl Into<String>) -> Self {
        Self {
            id,
            name: name.into(),
            self_weight: false,
        }
    }

    /// Marks the case as carrying self-weight.
    pub fn with_self_weight(mut self) -> Self {
        self.self_weight = true;
        self
    }
}

/// A force and moment applied directly to a node.
///
/// Forces are in kN, moments in kN·m, in global axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NodalLoad {
    /// Node receiving the load.
    pub node: NodeId,
    /// Global force components X, Y, Z in kN.
    pub force: [Real; 3],
    /// Global moment components X, Y, Z in kN·m.
    pub moment: [Real; 3],
}

impl NodalLoad {
    /// Creates a nodal load from global components.
    pub const fn new(node: NodeId, force: [Real; 3], moment: [Real; 3]) -> Self {
        Self {
            node,
            force,
            moment,
        }
    }
}

/// A distributed load on a member, in global axes per unit length.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MemberLoad {
    /// Member receiving the load.
    pub member: MemberId,
    /// Global force per length, X, Y, Z in kN/m.
    pub force: [Real; 3],
    /// Global moment per length, X, Y, Z in kN.
    pub moment: [Real; 3],
}

impl MemberLoad {
    /// Creates a member load from global components.
    pub const fn new(member: MemberId, force: [Real; 3], moment: [Real; 3]) -> Self {
        Self {
            member,
            force,
            moment,
        }
    }
}

/// A pressure on an area, positive along global Z, in kN/m².
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AreaLoad {
    /// Area receiving the load.
    pub area: AreaId,
    /// Pressure in kN/m².
    pub pressure: Real,
}

impl AreaLoad {
    /// Creates an area pressure.
    pub const fn new(area: AreaId, pressure: Real) -> Self {
        Self { area, pressure }
    }
}

/// A pressure on a spline patch, positive along global Z, in kN/m².
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PatchLoad {
    /// Patch receiving the load.
    pub patch: PatchId,
    /// Pressure in kN/m².
    pub pressure: Real,
}

impl PatchLoad {
    /// Creates a patch pressure.
    pub const fn new(patch: PatchId, pressure: Real) -> Self {
        Self { patch, pressure }
    }
}

/// One applied load and the case it belongs to.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Load {
    /// A load applied to a node.
    Nodal {
        /// Owning case.
        case: LoadCaseId,
        /// Load value.
        load: NodalLoad,
    },
    /// A distributed load applied to a member.
    Member {
        /// Owning case.
        case: LoadCaseId,
        /// Load value.
        load: MemberLoad,
    },
    /// A pressure applied to an area.
    Area {
        /// Owning case.
        case: LoadCaseId,
        /// Load value.
        load: AreaLoad,
    },
    /// A pressure applied to a spline patch.
    Patch {
        /// Owning case.
        case: LoadCaseId,
        /// Load value.
        load: PatchLoad,
    },
}

impl Load {
    /// The case that owns this load.
    pub const fn case(&self) -> LoadCaseId {
        match self {
            Load::Nodal { case, .. }
            | Load::Member { case, .. }
            | Load::Area { case, .. }
            | Load::Patch { case, .. } => *case,
        }
    }

    fn is_finite(&self) -> bool {
        match self {
            Load::Nodal { load, .. } => load
                .force
                .iter()
                .chain(load.moment.iter())
                .all(|v| v.is_finite()),
            Load::Member { load, .. } => load
                .force
                .iter()
                .chain(load.moment.iter())
                .all(|v| v.is_finite()),
            Load::Area { load, .. } => load.pressure.is_finite(),
            Load::Patch { load, .. } => load.pressure.is_finite(),
        }
    }
}

/// A named linear combination of load cases.
///
/// A model with no combinations receives one implicit combination at solve
/// time; the document never has to store that default.
#[derive(Clone, Debug, PartialEq)]
pub struct Combination {
    /// Stable document id.
    pub id: CombinationId,
    /// Human-facing combination name.
    pub name: String,
    /// Case id and factor pairs.
    pub factors: Vec<(LoadCaseId, Real)>,
}

impl Combination {
    /// Creates a combination from case id and factor pairs.
    pub fn new(
        id: CombinationId,
        name: impl Into<String>,
        factors: Vec<(LoadCaseId, Real)>,
    ) -> Self {
        Self {
            id,
            name: name.into(),
            factors,
        }
    }
}

/// Six translational and rotational restraint flags at a node.
///
/// The order is fixed here so every consumer reads the same DOFs. `true`
/// means restrained.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Restraints {
    /// Translation along global X.
    pub tx: bool,
    /// Translation along global Y.
    pub ty: bool,
    /// Translation along global Z.
    pub tz: bool,
    /// Rotation about global X.
    pub rx: bool,
    /// Rotation about global Y.
    pub ry: bool,
    /// Rotation about global Z.
    pub rz: bool,
}

impl Restraints {
    /// Creates a restraint set from the six named flags.
    pub const fn new(tx: bool, ty: bool, tz: bool, rx: bool, ry: bool, rz: bool) -> Self {
        Self {
            tx,
            ty,
            tz,
            rx,
            ry,
            rz,
        }
    }

    /// A node with all six DOFs free.
    pub const fn free() -> Self {
        Self::new(false, false, false, false, false, false)
    }

    /// A node with all six DOFs restrained.
    pub const fn fixed() -> Self {
        Self::new(true, true, true, true, true, true)
    }
}

/// Stable identity of a named support definition.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SupportDefinitionId(u32);

impl SupportDefinitionId {
    /// Wraps a document id.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// The underlying document id.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// A named set of six-DOF restraints users apply by group.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportDefinition {
    /// Stable document id.
    pub id: SupportDefinitionId,
    /// Human-facing definition name.
    pub name: String,
    /// Restraints this definition applies.
    pub restraints: Restraints,
}

impl SupportDefinition {
    /// Creates a support definition.
    pub fn new(id: SupportDefinitionId, name: impl Into<String>, restraints: Restraints) -> Self {
        Self {
            id,
            name: name.into(),
            restraints,
        }
    }
}

/// One node claiming one support definition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Support {
    /// Node receiving the restraints.
    pub node: NodeId,
    /// Definition supplying the restraints.
    pub definition: SupportDefinitionId,
}

impl Support {
    /// Creates a support assignment.
    pub const fn new(node: NodeId, definition: SupportDefinitionId) -> Self {
        Self { node, definition }
    }
}

/// The kind of model object a validation issue points at.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ObjectKind {
    /// The model as a whole, with no single owning object.
    Model,
    /// A node.
    Node,
    /// A material.
    Material,
    /// A section.
    Section,
    /// A member.
    Member,
    /// An area.
    Area,
    /// A spline patch.
    Patch,
    /// A load case.
    LoadCase,
    /// A load combination.
    Combination,
    /// A named support definition.
    SupportDefinition,
    /// A declared patch tie.
    PatchTie,
}

impl ObjectKind {
    /// Stable spelling used by reports and the CLI.
    pub const fn as_str(self) -> &'static str {
        match self {
            ObjectKind::Model => "model",
            ObjectKind::Node => "node",
            ObjectKind::Material => "material",
            ObjectKind::Section => "section",
            ObjectKind::Member => "member",
            ObjectKind::Area => "area",
            ObjectKind::Patch => "patch",
            ObjectKind::LoadCase => "load_case",
            ObjectKind::Combination => "combination",
            ObjectKind::SupportDefinition => "support_definition",
            ObjectKind::PatchTie => "patch_tie",
        }
    }

    /// Parses a report spelling, case-insensitively.
    pub fn parse(text: &str) -> Option<Self> {
        match text.trim().to_ascii_lowercase().as_str() {
            "model" => Some(ObjectKind::Model),
            "node" => Some(ObjectKind::Node),
            "material" => Some(ObjectKind::Material),
            "section" => Some(ObjectKind::Section),
            "member" => Some(ObjectKind::Member),
            "area" => Some(ObjectKind::Area),
            "patch" => Some(ObjectKind::Patch),
            "load_case" => Some(ObjectKind::LoadCase),
            "combination" => Some(ObjectKind::Combination),
            "support_definition" => Some(ObjectKind::SupportDefinition),
            "patch_tie" => Some(ObjectKind::PatchTie),
            _ => None,
        }
    }
}

/// A stable reference to the object responsible for an issue.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ObjectRef {
    /// Kind of object.
    pub kind: ObjectKind,
    /// Document id of the object, or `None` for a model-wide issue.
    pub id: Option<u32>,
}

impl ObjectRef {
    /// Points at one object by kind and document id.
    pub const fn new(kind: ObjectKind, id: u32) -> Self {
        Self { kind, id: Some(id) }
    }

    /// Points at the model as a whole.
    pub const fn model() -> Self {
        Self {
            kind: ObjectKind::Model,
            id: None,
        }
    }
}

/// Whether a validation issue blocks analysis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Severity {
    /// The model cannot be analysed until this is fixed.
    Blocking,
    /// The model can be analysed, but the object deserves attention.
    Warning,
}

impl Severity {
    /// Stable spelling used by reports and the CLI.
    pub const fn as_str(self) -> &'static str {
        match self {
            Severity::Blocking => "blocking",
            Severity::Warning => "warning",
        }
    }

    /// True when analysis must stop.
    pub const fn is_blocking(self) -> bool {
        matches!(self, Severity::Blocking)
    }
}

/// One machine-readable validation result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ValidationIssue {
    /// Object the issue belongs to.
    pub object: ObjectRef,
    /// Whether the issue blocks analysis.
    pub severity: Severity,
    /// Human-readable explanation.
    pub message: String,
}

impl ValidationIssue {
    /// Creates an issue with an explicit severity.
    pub fn new(object: ObjectRef, severity: Severity, message: impl Into<String>) -> Self {
        Self {
            object,
            severity,
            message: message.into(),
        }
    }

    /// Creates a blocking issue.
    pub fn blocking(object: ObjectRef, message: impl Into<String>) -> Self {
        Self::new(object, Severity::Blocking, message)
    }

    /// Creates a warning issue.
    pub fn warning(object: ObjectRef, message: impl Into<String>) -> Self {
        Self::new(object, Severity::Warning, message)
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
    load_cases: Vec<LoadCase>,
    loads: Vec<Load>,
    combinations: Vec<Combination>,
    support_definitions: Vec<SupportDefinition>,
    supports: Vec<Support>,
    node_restraints: HashMap<NodeId, Restraints>,
    patch_ties: Vec<PatchTie>,
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

    /// All load cases in insertion order.
    pub fn load_cases(&self) -> &[LoadCase] {
        &self.load_cases
    }

    /// All applied loads in insertion order.
    pub fn loads(&self) -> &[Load] {
        &self.loads
    }

    /// All load combinations in insertion order.
    pub fn combinations(&self) -> &[Combination] {
        &self.combinations
    }

    /// All named support definitions in insertion order.
    pub fn support_definitions(&self) -> &[SupportDefinition] {
        &self.support_definitions
    }

    /// All node assignments to support definitions.
    pub fn supports(&self) -> &[Support] {
        &self.supports
    }

    /// All declared patch ties.
    pub fn patch_ties(&self) -> &[PatchTie] {
        &self.patch_ties
    }

    /// Reconciled restraints for a node, after `apply_supports`.
    pub fn node_restraints(&self, node: NodeId) -> Option<Restraints> {
        self.node_restraints.get(&node).copied()
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

    /// Adds a load case, refusing a duplicate id.
    pub fn add_load_case(&mut self, case: LoadCase) -> Result<(), ModelError> {
        if self
            .load_cases
            .iter()
            .any(|existing| existing.id == case.id)
        {
            return Err(ModelError::DuplicateLoadCase { id: case.id });
        }
        self.load_cases.push(case);
        Ok(())
    }

    /// Adds a load, refusing an unknown case, target or non-finite value.
    pub fn add_load(&mut self, load: Load) -> Result<(), ModelError> {
        let case = load.case();
        if !self.load_cases.iter().any(|existing| existing.id == case) {
            return Err(ModelError::UnknownLoadCase { case });
        }
        if !load.is_finite() {
            return Err(ModelError::NonFiniteLoad { case });
        }
        match load {
            Load::Nodal { load, .. } => {
                if !self.nodes.iter().any(|node| node.id == load.node) {
                    return Err(ModelError::UnknownLoadNode {
                        case,
                        node: load.node,
                    });
                }
            }
            Load::Member { load, .. } => {
                if !self.members.iter().any(|member| member.id == load.member) {
                    return Err(ModelError::UnknownLoadMember {
                        case,
                        member: load.member,
                    });
                }
            }
            Load::Area { load, .. } => {
                if !self.areas.iter().any(|area| area.id == load.area) {
                    return Err(ModelError::UnknownLoadArea {
                        case,
                        area: load.area,
                    });
                }
            }
            Load::Patch { load, .. } => {
                if !self.patches.iter().any(|patch| patch.id == load.patch) {
                    return Err(ModelError::UnknownLoadPatch {
                        case,
                        patch: load.patch,
                    });
                }
            }
        }
        self.loads.push(load);
        Ok(())
    }

    /// Adds a combination, refusing a duplicate id or a bad case factor.
    pub fn add_combination(&mut self, combination: Combination) -> Result<(), ModelError> {
        if self
            .combinations
            .iter()
            .any(|existing| existing.id == combination.id)
        {
            return Err(ModelError::DuplicateCombination { id: combination.id });
        }
        let mut seen: Vec<LoadCaseId> = Vec::new();
        for &(case, factor) in &combination.factors {
            if !factor.is_finite() {
                return Err(ModelError::NonFiniteCombinationFactor {
                    combination: combination.id,
                    case,
                });
            }
            if seen.contains(&case) {
                return Err(ModelError::DuplicateCombinationCase {
                    combination: combination.id,
                    case,
                });
            }
            seen.push(case);
            if !self.load_cases.iter().any(|existing| existing.id == case) {
                return Err(ModelError::UnknownCombinationCase {
                    combination: combination.id,
                    case,
                });
            }
        }
        self.combinations.push(combination);
        Ok(())
    }

    /// Adds a support definition, refusing a duplicate id.
    pub fn add_support_definition(
        &mut self,
        definition: SupportDefinition,
    ) -> Result<(), ModelError> {
        if self
            .support_definitions
            .iter()
            .any(|existing| existing.id == definition.id)
        {
            return Err(ModelError::DuplicateSupportDefinition { id: definition.id });
        }
        self.support_definitions.push(definition);
        Ok(())
    }

    /// Replaces a support definition by id so a group edit can be reconciled.
    pub fn set_support_definition(
        &mut self,
        definition: SupportDefinition,
    ) -> Result<SupportDefinition, ModelError> {
        let existing = self
            .support_definitions
            .iter_mut()
            .find(|current| current.id == definition.id)
            .ok_or(ModelError::UnknownSupportDefinition { id: definition.id })?;
        Ok(std::mem::replace(existing, definition))
    }

    /// Adds a support assignment, refusing a second claim on the same node.
    pub fn add_support(&mut self, support: Support) -> Result<(), ModelError> {
        if let Some(existing) = self.supports.iter().find(|s| s.node == support.node) {
            return Err(ModelError::MultiClaimedNode {
                node: support.node,
                first: existing.definition,
                second: support.definition,
            });
        }
        if !self.nodes.iter().any(|node| node.id == support.node) {
            return Err(ModelError::UnknownSupportNode {
                support: support.node,
            });
        }
        if !self
            .support_definitions
            .iter()
            .any(|definition| definition.id == support.definition)
        {
            return Err(ModelError::UnknownSupportDefinition {
                id: support.definition,
            });
        }
        self.supports.push(support);
        Ok(())
    }

    /// Removes a support assignment, leaving the mask stale until reconcile.
    pub fn remove_support(&mut self, node: NodeId) -> Option<Support> {
        let position = self
            .supports
            .iter()
            .position(|support| support.node == node)?;
        Some(self.supports.remove(position))
    }

    /// Rebuilds node restraints from the current assignments and definitions.
    ///
    /// The derived map is cleared first, so a removed or renamed assignment
    /// cannot leave a stale mask behind.
    pub fn apply_supports(&mut self) -> Result<(), ModelError> {
        self.node_restraints.clear();
        for support in &self.supports {
            let definition = self
                .support_definitions
                .iter()
                .find(|current| current.id == support.definition)
                .ok_or(ModelError::UnknownSupportDefinition {
                    id: support.definition,
                })?;
            self.node_restraints
                .insert(support.node, definition.restraints);
        }
        Ok(())
    }

    /// Adds a tie declaration, refusing a bad seam or a second master.
    pub fn add_patch_tie(&mut self, tie: PatchTie) -> Result<(), ModelError> {
        if self.patch_ties.iter().any(|existing| existing.id == tie.id) {
            return Err(ModelError::DuplicatePatchTie { id: tie.id });
        }
        check_tie(&self.patches, &self.patch_ties, &tie)?;
        self.patch_ties.push(tie);
        Ok(())
    }

    /// Resolves how two patch edges are coupled, if they touch at all.
    pub fn patch_coupling(&self, first: PatchId, second: PatchId) -> Option<PatchCoupling> {
        coupling_between(&self.patches, &self.patch_ties, first, second)
    }

    /// Finds member ends that coincide with patch edge control points.
    ///
    /// Only control-point coincidence is detected. A member end in the
    /// interior of an edge curve needs spline evaluation and is out of scope
    /// for the model layer.
    pub fn beam_couplings(&self) -> Vec<BeamCoupling> {
        let mut couplings = Vec::new();
        for member in &self.members {
            for (node, end) in [(member.start, BeamEnd::Start), (member.end, BeamEnd::End)] {
                for patch in &self.patches {
                    for edge in PatchEdge::ALL {
                        if patch.edge_nodes(edge).contains(&node) {
                            couplings.push(BeamCoupling::new(member.id, end, patch.id, edge, node));
                        }
                    }
                }
            }
        }
        couplings
    }

    /// Builds exact control-point constraint rows from declared ties.
    ///
    /// A node that sits on both the master and the slave edge is already one
    /// unknown, so it produces no row. Interior mortar constraints are an
    /// analysis concern and are not produced here.
    pub fn tie_constraints(&self) -> Result<Vec<TieConstraint>, ModelError> {
        self.validate_ties()?;
        let mut rows = Vec::new();
        for tie in &self.patch_ties {
            let master =
                patch_by_id(&self.patches, tie.master).ok_or(ModelError::UnknownPatchTiePatch {
                    tie: tie.id,
                    patch: tie.master,
                })?;
            let slave =
                patch_by_id(&self.patches, tie.slave).ok_or(ModelError::UnknownPatchTiePatch {
                    tie: tie.id,
                    patch: tie.slave,
                })?;
            let masters = master.edge_nodes(tie.master_edge);
            let slaves = slave.edge_nodes(tie.slave_edge);
            for (slave_node, master_node) in slaves.into_iter().zip(masters) {
                if slave_node == master_node {
                    continue;
                }
                rows.push(TieConstraint {
                    slave: slave_node,
                    masters: vec![master_node],
                    weights: vec![1.0],
                });
            }
        }
        Ok(rows)
    }

    fn validate_ties(&self) -> Result<(), ModelError> {
        check_ties(&self.patches, &self.patch_ties)
    }

    /// Checks the model's structural rules.
    ///
    /// This is the blocking contract. Non-blocking observations live in
    /// [`Model::validate_issues`].
    pub fn validate(&self) -> Result<(), ModelError> {
        self.validate_references()?;
        self.validate_ties()
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
        for load in &self.loads {
            let case = load.case();
            if index.load_case(case).is_none() {
                return Err(ModelError::UnknownLoadCase { case });
            }
            match load {
                Load::Nodal { load, .. } => {
                    if index.node(load.node).is_none() {
                        return Err(ModelError::UnknownLoadNode {
                            case,
                            node: load.node,
                        });
                    }
                }
                Load::Member { load, .. } => {
                    if index.member(load.member).is_none() {
                        return Err(ModelError::UnknownLoadMember {
                            case,
                            member: load.member,
                        });
                    }
                }
                Load::Area { load, .. } => {
                    if index.area(load.area).is_none() {
                        return Err(ModelError::UnknownLoadArea {
                            case,
                            area: load.area,
                        });
                    }
                }
                Load::Patch { load, .. } => {
                    if index.patch(load.patch).is_none() {
                        return Err(ModelError::UnknownLoadPatch {
                            case,
                            patch: load.patch,
                        });
                    }
                }
            }
        }
        for combination in &self.combinations {
            let mut seen: Vec<LoadCaseId> = Vec::new();
            for &(case, _) in &combination.factors {
                if seen.contains(&case) {
                    return Err(ModelError::DuplicateCombinationCase {
                        combination: combination.id,
                        case,
                    });
                }
                seen.push(case);
                if index.load_case(case).is_none() {
                    return Err(ModelError::UnknownCombinationCase {
                        combination: combination.id,
                        case,
                    });
                }
            }
        }
        let mut claimed: Vec<(NodeId, SupportDefinitionId)> = Vec::new();
        for support in &self.supports {
            if let Some((_, first)) = claimed.iter().find(|(node, _)| *node == support.node) {
                return Err(ModelError::MultiClaimedNode {
                    node: support.node,
                    first: *first,
                    second: support.definition,
                });
            }
            claimed.push((support.node, support.definition));
            if index.node(support.node).is_none() {
                return Err(ModelError::UnknownSupportNode {
                    support: support.node,
                });
            }
            if index.support_definition(support.definition).is_none() {
                return Err(ModelError::UnknownSupportDefinition {
                    id: support.definition,
                });
            }
        }
        Ok(())
    }

    /// Returns blocking errors and non-blocking observations as one list.
    ///
    /// Structural errors come first. If any exist, warnings are not computed
    /// because their reference sets are not trustworthy yet.
    pub fn validate_issues(&self) -> Vec<ValidationIssue> {
        if let Err(error) = self.validate() {
            return vec![error.to_issue()];
        }

        let mut referenced: Vec<NodeId> = Vec::new();
        for member in &self.members {
            referenced.push(member.start);
            referenced.push(member.end);
        }
        for area in &self.areas {
            referenced.extend(area.nodes);
        }
        for patch in &self.patches {
            for row in &patch.control_net {
                referenced.extend(row.iter().copied());
            }
        }
        for support in &self.supports {
            referenced.push(support.node);
        }
        for load in &self.loads {
            if let Load::Nodal { load, .. } = load {
                referenced.push(load.node);
            }
        }

        let mut issues = Vec::new();
        for node in &self.nodes {
            if !referenced.contains(&node.id) {
                issues.push(ValidationIssue::warning(
                    ObjectRef::new(ObjectKind::Node, node.id.get()),
                    "orphan node is not used by any object",
                ));
            }
        }
        issues
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
    /// A load case id was used twice.
    DuplicateLoadCase { id: LoadCaseId },
    /// A combination id was used twice.
    DuplicateCombination { id: CombinationId },
    /// A load references a case that is not in the model.
    UnknownLoadCase { case: LoadCaseId },
    /// A load references a node that is not in the model.
    UnknownLoadNode { case: LoadCaseId, node: NodeId },
    /// A load references a member that is not in the model.
    UnknownLoadMember { case: LoadCaseId, member: MemberId },
    /// A load references an area that is not in the model.
    UnknownLoadArea { case: LoadCaseId, area: AreaId },
    /// A load references a patch that is not in the model.
    UnknownLoadPatch { case: LoadCaseId, patch: PatchId },
    /// A load carries a non-finite value.
    NonFiniteLoad { case: LoadCaseId },
    /// A combination references a case that is not in the model.
    UnknownCombinationCase {
        combination: CombinationId,
        case: LoadCaseId,
    },
    /// A combination factor is not finite.
    NonFiniteCombinationFactor {
        combination: CombinationId,
        case: LoadCaseId,
    },
    /// A combination lists the same case twice.
    DuplicateCombinationCase {
        combination: CombinationId,
        case: LoadCaseId,
    },
    /// A support definition id was used twice.
    DuplicateSupportDefinition { id: SupportDefinitionId },
    /// A support assignment references an unknown definition.
    UnknownSupportDefinition { id: SupportDefinitionId },
    /// A support assignment references an unknown node.
    UnknownSupportNode { support: NodeId },
    /// Two support definitions claim the same node.
    MultiClaimedNode {
        node: NodeId,
        first: SupportDefinitionId,
        second: SupportDefinitionId,
    },
    /// A patch tie id was used twice.
    DuplicatePatchTie { id: PatchTieId },
    /// A patch tie references a patch that is not in the model.
    UnknownPatchTiePatch { tie: PatchTieId, patch: PatchId },
    /// A patch tie names the same patch as master and slave.
    PatchTieSamePatch { tie: PatchTieId, patch: PatchId },
    /// A tie joins edges whose control points cannot be matched one to one.
    PatchTieEdgeShape { tie: PatchTieId },
    /// A slave edge is shared by two different master patches.
    PatchTieSlaveTwoMasters {
        node: NodeId,
        first: PatchId,
        second: PatchId,
    },
    /// Two patches share part of an edge without a tie declaration.
    TJunctionWithoutTie {
        first: PatchId,
        second: PatchId,
        node: NodeId,
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
            ModelError::DuplicateLoadCase { id } => {
                write!(f, "duplicate load case id {}", id.get())
            }
            ModelError::DuplicateCombination { id } => {
                write!(f, "duplicate combination id {}", id.get())
            }
            ModelError::UnknownLoadCase { case } => {
                write!(f, "load references unknown case {}", case.get())
            }
            ModelError::UnknownLoadNode { case, node } => write!(
                f,
                "case {} references unknown load node {}",
                case.get(),
                node.get()
            ),
            ModelError::UnknownLoadMember { case, member } => write!(
                f,
                "case {} references unknown load member {}",
                case.get(),
                member.get()
            ),
            ModelError::UnknownLoadArea { case, area } => write!(
                f,
                "case {} references unknown load area {}",
                case.get(),
                area.get()
            ),
            ModelError::UnknownLoadPatch { case, patch } => write!(
                f,
                "case {} references unknown load patch {}",
                case.get(),
                patch.get()
            ),
            ModelError::NonFiniteLoad { case } => {
                write!(f, "case {} has a non-finite load value", case.get())
            }
            ModelError::UnknownCombinationCase { combination, case } => write!(
                f,
                "combination {} references unknown case {}",
                combination.get(),
                case.get()
            ),
            ModelError::NonFiniteCombinationFactor { combination, case } => write!(
                f,
                "combination {} has a non-finite factor for case {}",
                combination.get(),
                case.get()
            ),
            ModelError::DuplicateCombinationCase { combination, case } => write!(
                f,
                "combination {} lists case {} twice",
                combination.get(),
                case.get()
            ),
            ModelError::DuplicateSupportDefinition { id } => {
                write!(f, "duplicate support definition id {}", id.get())
            }
            ModelError::UnknownSupportDefinition { id } => {
                write!(f, "unknown support definition id {}", id.get())
            }
            ModelError::UnknownSupportNode { support } => {
                write!(f, "support references unknown node {}", support.get())
            }
            ModelError::MultiClaimedNode {
                node,
                first,
                second,
            } => write!(
                f,
                "node {} is claimed by definitions {} and {}",
                node.get(),
                first.get(),
                second.get()
            ),
            ModelError::DuplicatePatchTie { id } => {
                write!(f, "duplicate patch tie id {}", id.get())
            }
            ModelError::UnknownPatchTiePatch { tie, patch } => write!(
                f,
                "patch tie {} references unknown patch {}",
                tie.get(),
                patch.get()
            ),
            ModelError::PatchTieEdgeShape { tie } => write!(
                f,
                "patch tie {} joins edges with different control-point layouts",
                tie.get()
            ),
            ModelError::PatchTieSamePatch { tie, patch } => write!(
                f,
                "patch tie {} names patch {} as both master and slave",
                tie.get(),
                patch.get()
            ),
            ModelError::PatchTieSlaveTwoMasters {
                node,
                first,
                second,
            } => write!(
                f,
                "slave node {} is tied to masters {} and {}",
                node.get(),
                first.get(),
                second.get()
            ),
            ModelError::TJunctionWithoutTie {
                first,
                second,
                node,
            } => write!(
                f,
                "patches {} and {} share only part of an edge at node {}; add a tie or match the control points",
                first.get(),
                second.get(),
                node.get()
            ),
        }
    }
}

impl std::error::Error for ModelError {}

impl ModelError {
    /// Converts a blocking error into the shared issue contract.
    pub fn to_issue(&self) -> ValidationIssue {
        let object = match self {
            ModelError::DuplicateNode { id } => ObjectRef::new(ObjectKind::Node, id.get()),
            ModelError::DuplicateMaterial { id } => ObjectRef::new(ObjectKind::Material, id.get()),
            ModelError::DuplicateSection { id } => ObjectRef::new(ObjectKind::Section, id.get()),
            ModelError::DuplicateMember { id } => ObjectRef::new(ObjectKind::Member, id.get()),
            ModelError::UnknownNode { member, .. } => {
                ObjectRef::new(ObjectKind::Member, member.get())
            }
            ModelError::UnknownSection { member, .. } => {
                ObjectRef::new(ObjectKind::Member, member.get())
            }
            ModelError::DuplicateArea { id } => ObjectRef::new(ObjectKind::Area, id.get()),
            ModelError::DuplicatePatch { id } => ObjectRef::new(ObjectKind::Patch, id.get()),
            ModelError::UnknownAreaNode { area, .. } => {
                ObjectRef::new(ObjectKind::Area, area.get())
            }
            ModelError::UnknownPatchNode { patch, .. } => {
                ObjectRef::new(ObjectKind::Patch, patch.get())
            }
            ModelError::PatchNetShape { patch }
            | ModelError::PatchKnotCount { patch, .. }
            | ModelError::PatchKnotsInvalid { patch, .. }
            | ModelError::PatchKnotsNotClamped { patch, .. }
            | ModelError::PatchDegree { patch, .. } => {
                ObjectRef::new(ObjectKind::Patch, patch.get())
            }
            ModelError::DuplicateLoadCase { id } => ObjectRef::new(ObjectKind::LoadCase, id.get()),
            ModelError::DuplicateCombination { id } => {
                ObjectRef::new(ObjectKind::Combination, id.get())
            }
            ModelError::UnknownLoadCase { case }
            | ModelError::UnknownLoadNode { case, .. }
            | ModelError::UnknownLoadMember { case, .. }
            | ModelError::UnknownLoadArea { case, .. }
            | ModelError::UnknownLoadPatch { case, .. }
            | ModelError::NonFiniteLoad { case } => {
                ObjectRef::new(ObjectKind::LoadCase, case.get())
            }
            ModelError::UnknownCombinationCase { combination, .. }
            | ModelError::NonFiniteCombinationFactor { combination, .. }
            | ModelError::DuplicateCombinationCase { combination, .. } => {
                ObjectRef::new(ObjectKind::Combination, combination.get())
            }
            ModelError::DuplicateSupportDefinition { id }
            | ModelError::UnknownSupportDefinition { id } => {
                ObjectRef::new(ObjectKind::SupportDefinition, id.get())
            }
            ModelError::UnknownSupportNode { support } => {
                ObjectRef::new(ObjectKind::Node, support.get())
            }
            ModelError::MultiClaimedNode { node, .. } => {
                ObjectRef::new(ObjectKind::Node, node.get())
            }
            ModelError::DuplicatePatchTie { id }
            | ModelError::UnknownPatchTiePatch { tie: id, .. }
            | ModelError::PatchTieSamePatch { tie: id, .. }
            | ModelError::PatchTieEdgeShape { tie: id } => {
                ObjectRef::new(ObjectKind::PatchTie, id.get())
            }
            ModelError::PatchTieSlaveTwoMasters { first, .. }
            | ModelError::TJunctionWithoutTie { first, .. } => {
                ObjectRef::new(ObjectKind::Patch, first.get())
            }
        };
        ValidationIssue::blocking(object, self.to_string())
    }
}

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
    load_cases: HashMap<LoadCaseId, usize>,
    combinations: HashMap<CombinationId, usize>,
    support_definitions: HashMap<SupportDefinitionId, usize>,
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
        for (position, case) in model.load_cases().iter().enumerate() {
            let inserted = index.load_cases.insert(case.id, position);
            debug_assert!(inserted.is_none(), "duplicate load case id in model");
        }
        for (position, combination) in model.combinations().iter().enumerate() {
            let inserted = index.combinations.insert(combination.id, position);
            debug_assert!(inserted.is_none(), "duplicate combination id in model");
        }
        for (position, definition) in model.support_definitions().iter().enumerate() {
            let inserted = index.support_definitions.insert(definition.id, position);
            debug_assert!(
                inserted.is_none(),
                "duplicate support definition id in model"
            );
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

    /// Position of a load case id in `Model::load_cases`.
    pub fn load_case(&self, id: LoadCaseId) -> Option<usize> {
        self.load_cases.get(&id).copied()
    }

    /// Position of a combination id in `Model::combinations`.
    pub fn combination(&self, id: CombinationId) -> Option<usize> {
        self.combinations.get(&id).copied()
    }

    /// Position of a support definition id in `Model::support_definitions`.
    pub fn support_definition(&self, id: SupportDefinitionId) -> Option<usize> {
        self.support_definitions.get(&id).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(got: Real, want: Real) -> bool {
        (got - want).abs() <= 1.0e-12 * want.abs().max(1.0)
    }

    fn nodes_only() -> Model {
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
    }

    fn with_material_section() -> Model {
        let mut model = nodes_only();
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

    fn with_member() -> Model {
        let mut model = with_material_section();
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
    }

    fn full_model() -> Model {
        let mut model = with_member();
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
            .add_load_case(LoadCase::new(LoadCaseId::new(10), "LC1"))
            .unwrap();
        model
            .add_load(Load::Nodal {
                case: LoadCaseId::new(10),
                load: NodalLoad::new(NodeId::new(7), [1.0, 0.0, 0.0], [0.0; 3]),
            })
            .unwrap();
        model
            .add_combination(Combination::new(
                CombinationId::new(13),
                "ULS",
                vec![(LoadCaseId::new(10), 1.35)],
            ))
            .unwrap();
        model
            .add_support_definition(SupportDefinition::new(
                SupportDefinitionId::new(14),
                "PIN",
                Restraints::new(true, true, true, false, false, false),
            ))
            .unwrap();
        model
            .add_support(Support::new(NodeId::new(7), SupportDefinitionId::new(14)))
            .unwrap();
        model.apply_supports().unwrap();
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
        assert_eq!(LoadCaseId::new(17).get(), 17);
        assert_eq!(CombinationId::new(18).get(), 18);
        assert_eq!(SupportDefinitionId::new(19).get(), 19);
    }

    #[test]
    fn index_finds_each_entity_by_id() {
        let model = full_model();
        let index = model.index();

        assert_eq!(index.node(NodeId::new(7)), Some(0));
        assert_eq!(index.node(NodeId::new(3)), Some(1));
        assert_eq!(index.material(MaterialId::new(2)), Some(0));
        assert_eq!(index.section(SectionId::new(9)), Some(0));
        assert_eq!(index.member(MemberId::new(4)), Some(0));
        assert_eq!(index.area(AreaId::new(5)), Some(0));
        assert_eq!(index.patch(PatchId::new(6)), Some(0));
        assert_eq!(index.load_case(LoadCaseId::new(10)), Some(0));
        assert_eq!(index.combination(CombinationId::new(13)), Some(0));
        assert_eq!(
            index.support_definition(SupportDefinitionId::new(14)),
            Some(0)
        );
    }

    #[test]
    fn index_returns_none_for_an_unknown_id() {
        let index = full_model().index();
        assert_eq!(index.node(NodeId::new(99)), None);
        assert_eq!(index.material(MaterialId::new(99)), None);
        assert_eq!(index.section(SectionId::new(99)), None);
        assert_eq!(index.member(MemberId::new(99)), None);
        assert_eq!(index.area(AreaId::new(99)), None);
        assert_eq!(index.patch(PatchId::new(99)), None);
        assert_eq!(index.load_case(LoadCaseId::new(99)), None);
        assert_eq!(index.combination(CombinationId::new(99)), None);
        assert_eq!(index.support_definition(SupportDefinitionId::new(99)), None);
    }

    #[test]
    fn duplicate_ids_are_rejected_at_insertion() {
        let mut model = full_model();

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
        let section = Section::from_tables(
            SectionId::new(9),
            "test",
            1.0e6,
            1.0e12,
            5.0e11,
            2.0e11,
            3.0e18,
        );
        assert!(close(section.area, 1.0));
        assert!(close(section.inertia_y, 1.0));
        assert!(close(section.inertia_z, 5.0e-1));
        assert!(close(section.torsion_constant, 2.0e-1));
        assert!(close(section.warping_constant, 3.0));
    }

    #[test]
    fn material_modulus_converts_from_mpa() {
        let material = Material::from_mpa(MaterialId::new(2), "S355", 7850.0, 210_000.0, 0.3);
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
        assert_eq!(with_member().validate_references(), Ok(()));
    }

    #[test]
    fn missing_node_reference_is_rejected() {
        let mut model = with_material_section();
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
        let mut model = with_material_section();
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
        let mut model = nodes_only();
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
        let mut model = nodes_only();
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
        let mut model = nodes_only();
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
        let mut model = nodes_only();
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
        let mut model = nodes_only();
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

    fn add_nodes(model: &mut Model, ids: &[u32]) {
        for id in ids {
            model
                .add_node(Node::new(NodeId::new(*id), f64::from(*id), 0.0, 0.0))
                .unwrap();
        }
    }

    /// Two patches whose `xi0` rows overlap on two of three control points.
    fn patches_with_partial_edge() -> Model {
        let mut model = Model::new();
        add_nodes(&mut model, &[1, 2, 3, 4, 5, 6, 7, 8]);
        model
            .add_patch(SplinePatch::unweighted(
                PatchId::new(1),
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 0.5, 1.0, 1.0],
                vec![
                    vec![NodeId::new(1), NodeId::new(2), NodeId::new(3)],
                    vec![NodeId::new(4), NodeId::new(5), NodeId::new(6)],
                ],
            ))
            .unwrap();
        model
            .add_patch(SplinePatch::unweighted(
                PatchId::new(2),
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 0.5, 1.0, 1.0],
                vec![
                    vec![NodeId::new(1), NodeId::new(2), NodeId::new(7)],
                    vec![NodeId::new(4), NodeId::new(5), NodeId::new(8)],
                ],
            ))
            .unwrap();
        model
    }

    /// Two patches that already share a whole edge.
    fn patches_sharing_edge() -> Model {
        let mut model = Model::new();
        add_nodes(&mut model, &[1, 2, 3, 4, 5, 6]);
        model
            .add_patch(SplinePatch::unweighted(
                PatchId::new(1),
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
                vec![
                    vec![NodeId::new(1), NodeId::new(2)],
                    vec![NodeId::new(3), NodeId::new(4)],
                ],
            ))
            .unwrap();
        model
            .add_patch(SplinePatch::unweighted(
                PatchId::new(2),
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
                vec![
                    vec![NodeId::new(3), NodeId::new(4)],
                    vec![NodeId::new(5), NodeId::new(6)],
                ],
            ))
            .unwrap();
        model
    }

    #[test]
    fn patch_edges_follow_the_control_net() {
        let patch = SplinePatch::unweighted(
            PatchId::new(1),
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.5, 1.0, 1.0],
            vec![
                vec![NodeId::new(1), NodeId::new(2), NodeId::new(3)],
                vec![NodeId::new(4), NodeId::new(5), NodeId::new(6)],
            ],
        );
        assert_eq!(
            patch.edge_nodes(PatchEdge::Xi0),
            vec![NodeId::new(1), NodeId::new(2), NodeId::new(3)]
        );
        assert_eq!(
            patch.edge_nodes(PatchEdge::Xi1),
            vec![NodeId::new(4), NodeId::new(5), NodeId::new(6)]
        );
        assert_eq!(
            patch.edge_nodes(PatchEdge::Eta0),
            vec![NodeId::new(1), NodeId::new(4)]
        );
        assert_eq!(
            patch.edge_nodes(PatchEdge::Eta1),
            vec![NodeId::new(3), NodeId::new(6)]
        );
        assert_eq!(patch.edge_knots(PatchEdge::Xi0), &[0.0, 0.0, 0.5, 1.0, 1.0]);
        assert_eq!(patch.edge_knots(PatchEdge::Eta0), &[0.0, 0.0, 1.0, 1.0]);
    }

    #[test]
    fn patch_edge_spellings_round_trip() {
        for edge in PatchEdge::ALL {
            assert_eq!(PatchEdge::parse(edge.as_str()), Some(edge));
        }
        assert_eq!(PatchEdge::parse(" ETA1 "), Some(PatchEdge::Eta1));
        assert_eq!(PatchEdge::parse("xi2"), None);
    }

    #[test]
    fn a_fully_shared_edge_needs_no_tie() {
        let model = patches_sharing_edge();
        assert_eq!(
            model.patch_coupling(PatchId::new(1), PatchId::new(2)),
            Some(PatchCoupling::Shared { node_count: 2 })
        );
        assert_eq!(model.validate(), Ok(()));
    }

    /// Two patches that meet at one corner only, as in a grid.
    fn patches_touching_at_a_corner() -> Model {
        let mut model = Model::new();
        add_nodes(&mut model, &[1, 2, 3, 4, 5, 6, 7]);
        model
            .add_patch(SplinePatch::unweighted(
                PatchId::new(1),
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
                vec![
                    vec![NodeId::new(1), NodeId::new(2)],
                    vec![NodeId::new(3), NodeId::new(4)],
                ],
            ))
            .unwrap();
        model
            .add_patch(SplinePatch::unweighted(
                PatchId::new(2),
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
                vec![
                    vec![NodeId::new(4), NodeId::new(5)],
                    vec![NodeId::new(6), NodeId::new(7)],
                ],
            ))
            .unwrap();
        model
    }

    /// The tie that joins the overlapping `xi0` rows of patches 1 and 2.
    fn partial_edge_tie(id: u32) -> PatchTie {
        PatchTie::new(
            PatchTieId::new(id),
            PatchId::new(1),
            PatchEdge::Xi0,
            PatchId::new(2),
            PatchEdge::Xi0,
        )
    }

    /// Adds a patch whose `xi0` row has two control points where patch 1 has
    /// three.
    fn patches_with_uneven_edge() -> Model {
        let mut model = patches_with_partial_edge();
        add_nodes(&mut model, &[13]);
        model
            .add_patch(SplinePatch::unweighted(
                PatchId::new(4),
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
                vec![
                    vec![NodeId::new(1), NodeId::new(2)],
                    vec![NodeId::new(4), NodeId::new(13)],
                ],
            ))
            .unwrap();
        model
    }

    /// Adds a third patch that can reach the second patch's other slave edge.
    fn patches_with_three_patches() -> Model {
        let mut model = patches_with_partial_edge();
        add_nodes(&mut model, &[13, 14, 15, 16, 17, 18]);
        model
            .add_patch(SplinePatch::unweighted(
                PatchId::new(3),
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 0.5, 1.0, 1.0],
                vec![
                    vec![NodeId::new(13), NodeId::new(14), NodeId::new(17)],
                    vec![NodeId::new(15), NodeId::new(16), NodeId::new(18)],
                ],
            ))
            .unwrap();
        model
    }

    #[test]
    fn a_shared_corner_is_a_joint_not_a_seam() {
        let model = patches_touching_at_a_corner();
        assert_eq!(model.patch_coupling(PatchId::new(1), PatchId::new(2)), None);
        assert_eq!(model.validate(), Ok(()));
    }

    #[test]
    fn a_declared_tie_marks_the_pair_as_tied() {
        let mut model = patches_touching_at_a_corner();
        let tie = PatchTie::new(
            PatchTieId::new(1),
            PatchId::new(1),
            PatchEdge::Xi1,
            PatchId::new(2),
            PatchEdge::Xi0,
        );
        assert_eq!(model.add_patch_tie(tie), Ok(()));
        assert_eq!(
            model.patch_coupling(PatchId::new(1), PatchId::new(2)),
            Some(PatchCoupling::Tied {
                tie: PatchTieId::new(1)
            })
        );
    }

    #[test]
    fn a_tie_cannot_name_one_patch_twice() {
        let mut model = patches_with_partial_edge();
        let tie = PatchTie::new(
            PatchTieId::new(1),
            PatchId::new(1),
            PatchEdge::Xi0,
            PatchId::new(1),
            PatchEdge::Xi1,
        );
        assert_eq!(
            model.add_patch_tie(tie),
            Err(ModelError::PatchTieSamePatch {
                tie: PatchTieId::new(1),
                patch: PatchId::new(1),
            })
        );
    }

    #[test]
    fn coupling_with_an_unknown_patch_is_none() {
        let model = patches_with_partial_edge();
        assert_eq!(model.patch_coupling(PatchId::new(1), PatchId::new(9)), None);
        assert_eq!(model.patch_coupling(PatchId::new(9), PatchId::new(1)), None);
    }

    #[test]
    fn a_partial_edge_share_needs_a_tie() {
        let model = patches_with_partial_edge();
        assert_eq!(
            model.patch_coupling(PatchId::new(1), PatchId::new(2)),
            Some(PatchCoupling::NeedTie {
                node: NodeId::new(1)
            })
        );
        assert_eq!(
            model.validate(),
            Err(ModelError::TJunctionWithoutTie {
                first: PatchId::new(1),
                second: PatchId::new(2),
                node: NodeId::new(1),
            })
        );
    }

    #[test]
    fn a_tie_resolves_the_seam() {
        let mut model = patches_with_partial_edge();
        assert_eq!(model.add_patch_tie(partial_edge_tie(1)), Ok(()));
        assert_eq!(model.patch_ties().len(), 1);
        assert_eq!(
            model.patch_coupling(PatchId::new(1), PatchId::new(2)),
            Some(PatchCoupling::Tied {
                tie: PatchTieId::new(1)
            })
        );
        assert_eq!(model.validate(), Ok(()));
    }

    #[test]
    fn a_tie_needs_matching_control_point_counts() {
        let mut model = patches_with_uneven_edge();
        let tie = PatchTie::new(
            PatchTieId::new(1),
            PatchId::new(1),
            PatchEdge::Xi0,
            PatchId::new(4),
            PatchEdge::Xi0,
        );
        assert_eq!(
            model.add_patch_tie(tie),
            Err(ModelError::PatchTieEdgeShape {
                tie: PatchTieId::new(1)
            })
        );
    }

    #[test]
    fn a_tie_to_an_unknown_patch_is_rejected() {
        let mut model = patches_with_partial_edge();
        let tie = PatchTie::new(
            PatchTieId::new(1),
            PatchId::new(9),
            PatchEdge::Xi0,
            PatchId::new(2),
            PatchEdge::Xi0,
        );
        assert_eq!(
            model.add_patch_tie(tie),
            Err(ModelError::UnknownPatchTiePatch {
                tie: PatchTieId::new(1),
                patch: PatchId::new(9),
            })
        );
    }

    #[test]
    fn duplicate_tie_ids_are_rejected() {
        let mut model = patches_with_partial_edge();
        let tie = partial_edge_tie(1);
        model.add_patch_tie(tie).unwrap();
        assert_eq!(
            model.add_patch_tie(tie),
            Err(ModelError::DuplicatePatchTie {
                id: PatchTieId::new(1)
            })
        );
    }

    #[test]
    fn a_slave_edge_cannot_serve_two_masters() {
        let mut model = patches_with_three_patches();
        model.add_patch_tie(partial_edge_tie(1)).unwrap();
        let second = PatchTie::new(
            PatchTieId::new(2),
            PatchId::new(3),
            PatchEdge::Eta0,
            PatchId::new(2),
            PatchEdge::Eta0,
        );
        assert_eq!(
            model.add_patch_tie(second),
            Err(ModelError::PatchTieSlaveTwoMasters {
                node: NodeId::new(1),
                first: PatchId::new(1),
                second: PatchId::new(3),
            })
        );
    }

    #[test]
    fn a_tie_eliminates_only_the_unshared_slave_nodes() {
        let mut model = patches_with_partial_edge();
        model.add_patch_tie(partial_edge_tie(1)).unwrap();
        assert_eq!(
            model.tie_constraints().unwrap(),
            vec![TieConstraint {
                slave: NodeId::new(7),
                masters: vec![NodeId::new(3)],
                weights: vec![1.0],
            }]
        );
    }

    #[test]
    fn beam_ends_on_a_patch_edge_are_reported() {
        assert!(patches_with_partial_edge().beam_couplings().is_empty());
        assert_eq!(
            full_model().beam_couplings(),
            vec![
                BeamCoupling::new(
                    MemberId::new(4),
                    BeamEnd::Start,
                    PatchId::new(6),
                    PatchEdge::Xi0,
                    NodeId::new(7),
                ),
                BeamCoupling::new(
                    MemberId::new(4),
                    BeamEnd::Start,
                    PatchId::new(6),
                    PatchEdge::Eta0,
                    NodeId::new(7),
                ),
                BeamCoupling::new(
                    MemberId::new(4),
                    BeamEnd::End,
                    PatchId::new(6),
                    PatchEdge::Xi0,
                    NodeId::new(3),
                ),
                BeamCoupling::new(
                    MemberId::new(4),
                    BeamEnd::End,
                    PatchId::new(6),
                    PatchEdge::Eta1,
                    NodeId::new(3),
                ),
            ]
        );
    }

    #[test]
    fn tie_errors_name_the_tie_and_the_patch() {
        let shape = ModelError::PatchTieEdgeShape {
            tie: PatchTieId::new(1),
        };
        assert_eq!(
            shape.to_string(),
            "patch tie 1 joins edges with different control-point layouts"
        );
        assert_eq!(
            shape.to_issue().object,
            ObjectRef::new(ObjectKind::PatchTie, 1)
        );
        let seam = ModelError::TJunctionWithoutTie {
            first: PatchId::new(1),
            second: PatchId::new(2),
            node: NodeId::new(3),
        };
        assert!(seam.to_string().contains("patches 1 and 2"));
        assert_eq!(seam.to_issue().object, ObjectRef::new(ObjectKind::Patch, 1));
    }

    #[test]
    fn load_case_self_weight_is_explicit() {
        let plain = LoadCase::new(LoadCaseId::new(1), "D");
        let weighted = LoadCase::new(LoadCaseId::new(2), "SW").with_self_weight();
        assert!(!plain.self_weight);
        assert!(weighted.self_weight);
        assert_eq!(plain.name, "D");
    }

    #[test]
    fn full_model_load_references_validate() {
        assert_eq!(full_model().validate_references(), Ok(()));
    }

    #[test]
    fn load_without_case_is_rejected() {
        let mut model = nodes_only();
        let load = Load::Nodal {
            case: LoadCaseId::new(10),
            load: NodalLoad::new(NodeId::new(7), [1.0; 3], [0.0; 3]),
        };
        assert_eq!(
            model.add_load(load),
            Err(ModelError::UnknownLoadCase {
                case: LoadCaseId::new(10)
            })
        );
    }

    #[test]
    fn load_targets_are_checked() {
        let mut model = with_material_section();
        model
            .add_load_case(LoadCase::new(LoadCaseId::new(10), "LC"))
            .unwrap();
        let unknown_node = Load::Nodal {
            case: LoadCaseId::new(10),
            load: NodalLoad::new(NodeId::new(99), [1.0; 3], [0.0; 3]),
        };
        assert_eq!(
            model.add_load(unknown_node),
            Err(ModelError::UnknownLoadNode {
                case: LoadCaseId::new(10),
                node: NodeId::new(99),
            })
        );

        let unknown_member = Load::Member {
            case: LoadCaseId::new(10),
            load: MemberLoad::new(MemberId::new(99), [1.0; 3], [0.0; 3]),
        };
        assert_eq!(
            model.add_load(unknown_member),
            Err(ModelError::UnknownLoadMember {
                case: LoadCaseId::new(10),
                member: MemberId::new(99),
            })
        );
    }

    #[test]
    fn non_finite_load_is_rejected() {
        let mut model = nodes_only();
        model
            .add_load_case(LoadCase::new(LoadCaseId::new(10), "LC"))
            .unwrap();
        let load = Load::Nodal {
            case: LoadCaseId::new(10),
            load: NodalLoad::new(NodeId::new(7), [Real::NAN, 0.0, 0.0], [0.0; 3]),
        };
        assert_eq!(
            model.add_load(load),
            Err(ModelError::NonFiniteLoad {
                case: LoadCaseId::new(10)
            })
        );
    }

    #[test]
    fn combination_case_and_factor_are_checked() {
        let mut model = with_material_section();
        model
            .add_load_case(LoadCase::new(LoadCaseId::new(10), "LC"))
            .unwrap();
        let unknown = Combination::new(
            CombinationId::new(13),
            "ULS",
            vec![(LoadCaseId::new(99), 1.0)],
        );
        assert_eq!(
            model.add_combination(unknown),
            Err(ModelError::UnknownCombinationCase {
                combination: CombinationId::new(13),
                case: LoadCaseId::new(99),
            })
        );
        let non_finite = Combination::new(
            CombinationId::new(13),
            "ULS",
            vec![(LoadCaseId::new(10), Real::INFINITY)],
        );
        assert_eq!(
            model.add_combination(non_finite),
            Err(ModelError::NonFiniteCombinationFactor {
                combination: CombinationId::new(13),
                case: LoadCaseId::new(10),
            })
        );
        let repeated = Combination::new(
            CombinationId::new(13),
            "ULS",
            vec![(LoadCaseId::new(10), 1.0), (LoadCaseId::new(10), 1.0)],
        );
        assert_eq!(
            model.add_combination(repeated),
            Err(ModelError::DuplicateCombinationCase {
                combination: CombinationId::new(13),
                case: LoadCaseId::new(10),
            })
        );
    }

    #[test]
    fn duplicate_load_case_and_combination_ids_are_rejected() {
        let mut model = full_model();
        assert_eq!(
            model.add_load_case(LoadCase::new(LoadCaseId::new(10), "again")),
            Err(ModelError::DuplicateLoadCase {
                id: LoadCaseId::new(10)
            })
        );
        assert_eq!(
            model.add_combination(Combination::new(CombinationId::new(13), "again", vec![])),
            Err(ModelError::DuplicateCombination {
                id: CombinationId::new(13)
            })
        );
    }

    fn pin_definition() -> SupportDefinition {
        SupportDefinition::new(
            SupportDefinitionId::new(14),
            "PIN",
            Restraints::new(true, true, true, false, false, false),
        )
    }

    #[test]
    fn duplicate_support_definition_id_is_rejected() {
        let mut model = nodes_only();
        model.add_support_definition(pin_definition()).unwrap();
        assert_eq!(
            model.add_support_definition(pin_definition()),
            Err(ModelError::DuplicateSupportDefinition {
                id: SupportDefinitionId::new(14)
            })
        );
    }

    #[test]
    fn multi_claimed_node_is_rejected() {
        let mut model = nodes_only();
        model.add_support_definition(pin_definition()).unwrap();
        model
            .add_support_definition(SupportDefinition::new(
                SupportDefinitionId::new(15),
                "FIX",
                Restraints::fixed(),
            ))
            .unwrap();
        model
            .add_support(Support::new(NodeId::new(7), SupportDefinitionId::new(14)))
            .unwrap();
        assert_eq!(
            model.add_support(Support::new(NodeId::new(7), SupportDefinitionId::new(15))),
            Err(ModelError::MultiClaimedNode {
                node: NodeId::new(7),
                first: SupportDefinitionId::new(14),
                second: SupportDefinitionId::new(15),
            })
        );
    }

    #[test]
    fn support_references_are_checked() {
        let mut model = nodes_only();
        model.add_support_definition(pin_definition()).unwrap();
        assert_eq!(
            model.add_support(Support::new(NodeId::new(99), SupportDefinitionId::new(14))),
            Err(ModelError::UnknownSupportNode {
                support: NodeId::new(99)
            })
        );
        assert_eq!(
            model.add_support(Support::new(NodeId::new(7), SupportDefinitionId::new(99))),
            Err(ModelError::UnknownSupportDefinition {
                id: SupportDefinitionId::new(99)
            })
        );
    }

    #[test]
    fn apply_supports_builds_node_masks() {
        let mut model = nodes_only();
        model.add_support_definition(pin_definition()).unwrap();
        model
            .add_support(Support::new(NodeId::new(7), SupportDefinitionId::new(14)))
            .unwrap();
        model.apply_supports().unwrap();

        assert_eq!(
            model.node_restraints(NodeId::new(7)),
            Some(pin_definition().restraints)
        );
        assert_eq!(model.node_restraints(NodeId::new(3)), None);
    }

    #[test]
    fn editing_a_definition_reconciles_the_mask() {
        let mut model = nodes_only();
        model.add_support_definition(pin_definition()).unwrap();
        model
            .add_support(Support::new(NodeId::new(7), SupportDefinitionId::new(14)))
            .unwrap();
        model.apply_supports().unwrap();

        model
            .set_support_definition(SupportDefinition::new(
                SupportDefinitionId::new(14),
                "FIX",
                Restraints::fixed(),
            ))
            .unwrap();
        model.apply_supports().unwrap();
        assert_eq!(
            model.node_restraints(NodeId::new(7)),
            Some(Restraints::fixed())
        );
    }

    #[test]
    fn removing_a_support_clears_the_stale_mask() {
        let mut model = nodes_only();
        model.add_support_definition(pin_definition()).unwrap();
        model
            .add_support(Support::new(NodeId::new(7), SupportDefinitionId::new(14)))
            .unwrap();
        model.apply_supports().unwrap();
        assert!(model.node_restraints(NodeId::new(7)).is_some());

        assert_eq!(
            model.remove_support(NodeId::new(7)),
            Some(Support::new(NodeId::new(7), SupportDefinitionId::new(14)))
        );
        model.apply_supports().unwrap();
        assert_eq!(model.node_restraints(NodeId::new(7)), None);
    }

    #[test]
    fn set_support_definition_requires_an_existing_id() {
        let mut model = nodes_only();
        assert_eq!(
            model.set_support_definition(pin_definition()),
            Err(ModelError::UnknownSupportDefinition {
                id: SupportDefinitionId::new(14)
            })
        );
    }

    #[test]
    fn object_kind_spellings_round_trip() {
        for kind in [
            ObjectKind::Model,
            ObjectKind::Node,
            ObjectKind::Material,
            ObjectKind::Section,
            ObjectKind::Member,
            ObjectKind::Area,
            ObjectKind::Patch,
            ObjectKind::LoadCase,
            ObjectKind::Combination,
            ObjectKind::SupportDefinition,
            ObjectKind::PatchTie,
        ] {
            assert_eq!(ObjectKind::parse(kind.as_str()), Some(kind));
        }
        assert_eq!(ObjectKind::parse(" NODE "), Some(ObjectKind::Node));
        assert_eq!(ObjectKind::parse("load"), None);
    }

    #[test]
    fn severity_names_its_contract() {
        assert_eq!(Severity::Blocking.as_str(), "blocking");
        assert_eq!(Severity::Warning.as_str(), "warning");
        assert!(Severity::Blocking.is_blocking());
        assert!(!Severity::Warning.is_blocking());
    }

    #[test]
    fn model_ref_has_no_object_id() {
        assert_eq!(
            ObjectRef::model(),
            ObjectRef {
                kind: ObjectKind::Model,
                id: None,
            }
        );
    }

    #[test]
    fn blocking_error_maps_to_an_issue() {
        let issue = ModelError::UnknownNode {
            member: MemberId::new(4),
            node: NodeId::new(99),
        }
        .to_issue();
        assert_eq!(issue.object, ObjectRef::new(ObjectKind::Member, 4));
        assert_eq!(issue.severity, Severity::Blocking);
        assert!(issue.message.contains("member 4"));
    }

    #[test]
    fn full_model_validates_without_issues() {
        let model = full_model();
        assert_eq!(model.validate(), Ok(()));
        assert!(model.validate_issues().is_empty());
    }

    #[test]
    fn orphan_node_reports_a_warning_with_its_id() {
        let issues = with_member().validate_issues();
        assert_eq!(issues.len(), 2);
        assert!(issues.iter().any(|issue| {
            issue.object == ObjectRef::new(ObjectKind::Node, 11)
                && issue.severity == Severity::Warning
                && issue.message.contains("orphan")
        }));
        assert!(issues
            .iter()
            .any(|issue| issue.object == ObjectRef::new(ObjectKind::Node, 12)));
    }
}
