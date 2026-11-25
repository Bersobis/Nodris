//! Maps the v1 document onto the model entities.
//!
//! The field spellings here are the wire contract. The published schema under
//! `schema/` describes the same fields, and both are checked by the same tests.

use nodris_core::model::*;
use serde_json::{Map, Value};

use crate::interchange::{check_keys, ImportMode, ImportNote, ImportNoteKind, InterchangeError};

pub const ROOT_KEYS: &[&str] = &[
    "schemaVersion",
    "nodes",
    "materials",
    "sections",
    "members",
    "areas",
    "patches",
    "loadCases",
    "loads",
    "combinations",
    "supportDefinitions",
    "supports",
    "patchTies",
];

pub(crate) const NODE_KEYS: &[&str] = &["id", "x", "y", "z"];
pub(crate) const MATERIAL_KEYS: &[&str] =
    &["id", "name", "density", "elasticModulus", "poissonRatio"];
pub(crate) const SECTION_KEYS: &[&str] = &[
    "id",
    "name",
    "area",
    "inertiaY",
    "inertiaZ",
    "torsionConstant",
    "warpingConstant",
];
pub(crate) const MEMBER_KEYS: &[&str] = &["id", "start", "end", "section", "rotation", "kind"];
pub(crate) const AREA_KEYS: &[&str] = &["id", "nodes", "kind"];
pub(crate) const PATCH_KEYS: &[&str] = &[
    "id",
    "degreeU",
    "degreeV",
    "knotsU",
    "knotsV",
    "controlNet",
    "weights",
];
pub(crate) const LOAD_CASE_KEYS: &[&str] = &["id", "name", "selfWeight"];
pub(crate) const NODAL_LOAD_KEYS: &[&str] = &["case", "type", "node", "force", "moment"];
pub(crate) const MEMBER_LOAD_KEYS: &[&str] = &["case", "type", "member", "force", "moment"];
pub(crate) const AREA_LOAD_KEYS: &[&str] = &["case", "type", "area", "pressure"];
pub(crate) const PATCH_LOAD_KEYS: &[&str] = &["case", "type", "patch", "pressure"];
pub(crate) const COMBINATION_KEYS: &[&str] = &["id", "name", "factors"];
pub(crate) const FACTOR_KEYS: &[&str] = &["case", "factor"];
pub(crate) const SUPPORT_DEFINITION_KEYS: &[&str] = &["id", "name", "restraints"];
pub(crate) const RESTRAINT_KEYS: &[&str] = &["tx", "ty", "tz", "rx", "ry", "rz"];
pub(crate) const SUPPORT_KEYS: &[&str] = &["node", "definition"];
pub(crate) const PATCH_TIE_KEYS: &[&str] = &["id", "master", "masterEdge", "slave", "slaveEdge"];

fn number(path: &str, value: f64) -> Result<Value, InterchangeError> {
    serde_json::Number::from_f64(value)
        .map(Value::Number)
        .ok_or_else(|| InterchangeError::InvalidNumber {
            path: path.to_string(),
        })
}

fn numbers(path: &str, values: &[f64]) -> Result<Value, InterchangeError> {
    values
        .iter()
        .enumerate()
        .map(|(index, value)| number(&format!("{path}[{index}]"), *value))
        .collect::<Result<Vec<_>, _>>()
        .map(Value::Array)
}

fn vector(path: &str, values: [f64; 3]) -> Result<Value, InterchangeError> {
    numbers(path, &values)
}

fn matrix(path: &str, rows: &[Vec<f64>]) -> Result<Value, InterchangeError> {
    rows.iter()
        .enumerate()
        .map(|(index, row)| numbers(&format!("{path}[{index}]"), row))
        .collect::<Result<Vec<_>, _>>()
        .map(Value::Array)
}

fn ids(values: &[u32]) -> Value {
    Value::Array(values.iter().map(|id| Value::from(*id)).collect())
}

fn node_to_value(node: &Node) -> Result<Value, InterchangeError> {
    let mut object = Map::new();
    object.insert("id".into(), Value::from(node.id.get()));
    object.insert("x".into(), number("nodes[].x", node.x)?);
    object.insert("y".into(), number("nodes[].y", node.y)?);
    object.insert("z".into(), number("nodes[].z", node.z)?);
    Ok(Value::Object(object))
}

fn material_to_value(material: &Material) -> Result<Value, InterchangeError> {
    let mut object = Map::new();
    object.insert("id".into(), Value::from(material.id.get()));
    object.insert("name".into(), Value::from(material.name.clone()));
    object.insert(
        "density".into(),
        number("materials[].density", material.density)?,
    );
    object.insert(
        "elasticModulus".into(),
        number("materials[].elasticModulus", material.elastic_modulus)?,
    );
    object.insert(
        "poissonRatio".into(),
        number("materials[].poissonRatio", material.poisson_ratio)?,
    );
    Ok(Value::Object(object))
}

fn section_to_value(section: &Section) -> Result<Value, InterchangeError> {
    let mut object = Map::new();
    object.insert("id".into(), Value::from(section.id.get()));
    object.insert("name".into(), Value::from(section.name.clone()));
    object.insert("area".into(), number("sections[].area", section.area)?);
    object.insert(
        "inertiaY".into(),
        number("sections[].inertiaY", section.inertia_y)?,
    );
    object.insert(
        "inertiaZ".into(),
        number("sections[].inertiaZ", section.inertia_z)?,
    );
    object.insert(
        "torsionConstant".into(),
        number("sections[].torsionConstant", section.torsion_constant)?,
    );
    object.insert(
        "warpingConstant".into(),
        number("sections[].warpingConstant", section.warping_constant)?,
    );
    Ok(Value::Object(object))
}

fn member_to_value(member: &Member) -> Result<Value, InterchangeError> {
    let mut object = Map::new();
    object.insert("id".into(), Value::from(member.id.get()));
    object.insert("start".into(), Value::from(member.start.get()));
    object.insert("end".into(), Value::from(member.end.get()));
    object.insert("section".into(), Value::from(member.section.get()));
    object.insert(
        "rotation".into(),
        number("members[].rotation", member.rotation)?,
    );
    object.insert("kind".into(), Value::from(member.kind.as_str()));
    Ok(Value::Object(object))
}

fn area_to_value(area: &Area) -> Result<Value, InterchangeError> {
    let mut object = Map::new();
    object.insert("id".into(), Value::from(area.id.get()));
    object.insert("nodes".into(), ids(&area.nodes.map(|node| node.get())));
    object.insert("kind".into(), Value::from(area.kind.as_str()));
    Ok(Value::Object(object))
}

fn patch_to_value(patch: &SplinePatch) -> Result<Value, InterchangeError> {
    let mut object = Map::new();
    object.insert("id".into(), Value::from(patch.id.get()));
    object.insert("degreeU".into(), Value::from(patch.degree_u));
    object.insert("degreeV".into(), Value::from(patch.degree_v));
    object.insert(
        "knotsU".into(),
        numbers("patches[].knotsU", &patch.knots_u)?,
    );
    object.insert(
        "knotsV".into(),
        numbers("patches[].knotsV", &patch.knots_v)?,
    );
    object.insert(
        "controlNet".into(),
        Value::Array(
            patch
                .control_net
                .iter()
                .map(|row| ids(&row.iter().map(|id| id.get()).collect::<Vec<_>>()))
                .collect(),
        ),
    );
    object.insert(
        "weights".into(),
        matrix("patches[].weights", &patch.weights)?,
    );
    Ok(Value::Object(object))
}

fn load_case_to_value(case: &LoadCase) -> Result<Value, InterchangeError> {
    let mut object = Map::new();
    object.insert("id".into(), Value::from(case.id.get()));
    object.insert("name".into(), Value::from(case.name.clone()));
    object.insert("selfWeight".into(), Value::from(case.self_weight));
    Ok(Value::Object(object))
}

fn load_to_value(load: &Load) -> Result<Value, InterchangeError> {
    let mut object = Map::new();
    match load {
        Load::Nodal { case, load } => {
            object.insert("case".into(), Value::from(case.get()));
            object.insert("type".into(), Value::from("nodal"));
            object.insert("node".into(), Value::from(load.node.get()));
            object.insert("force".into(), vector("loads[].force", load.force)?);
            object.insert("moment".into(), vector("loads[].moment", load.moment)?);
        }
        Load::Member { case, load } => {
            object.insert("case".into(), Value::from(case.get()));
            object.insert("type".into(), Value::from("member"));
            object.insert("member".into(), Value::from(load.member.get()));
            object.insert("force".into(), vector("loads[].force", load.force)?);
            object.insert("moment".into(), vector("loads[].moment", load.moment)?);
        }
        Load::Area { case, load } => {
            object.insert("case".into(), Value::from(case.get()));
            object.insert("type".into(), Value::from("area"));
            object.insert("area".into(), Value::from(load.area.get()));
            object.insert(
                "pressure".into(),
                number("loads[].pressure", load.pressure)?,
            );
        }
        Load::Patch { case, load } => {
            object.insert("case".into(), Value::from(case.get()));
            object.insert("type".into(), Value::from("patch"));
            object.insert("patch".into(), Value::from(load.patch.get()));
            object.insert(
                "pressure".into(),
                number("loads[].pressure", load.pressure)?,
            );
        }
    }
    Ok(Value::Object(object))
}

fn combination_to_value(combination: &Combination) -> Result<Value, InterchangeError> {
    let factors = combination
        .factors
        .iter()
        .enumerate()
        .map(|(index, (case, factor))| {
            let mut object = Map::new();
            object.insert("case".into(), Value::from(case.get()));
            object.insert(
                "factor".into(),
                number(&format!("combinations[].factors[{index}].factor"), *factor)?,
            );
            Ok(Value::Object(object))
        })
        .collect::<Result<Vec<_>, InterchangeError>>()?;
    let mut object = Map::new();
    object.insert("id".into(), Value::from(combination.id.get()));
    object.insert("name".into(), Value::from(combination.name.clone()));
    object.insert("factors".into(), Value::Array(factors));
    Ok(Value::Object(object))
}

fn restraints_to_value(restraints: Restraints) -> Value {
    let mut object = Map::new();
    object.insert("tx".into(), Value::from(restraints.tx));
    object.insert("ty".into(), Value::from(restraints.ty));
    object.insert("tz".into(), Value::from(restraints.tz));
    object.insert("rx".into(), Value::from(restraints.rx));
    object.insert("ry".into(), Value::from(restraints.ry));
    object.insert("rz".into(), Value::from(restraints.rz));
    Value::Object(object)
}

fn support_definition_to_value(definition: &SupportDefinition) -> Value {
    let mut object = Map::new();
    object.insert("id".into(), Value::from(definition.id.get()));
    object.insert("name".into(), Value::from(definition.name.clone()));
    object.insert(
        "restraints".into(),
        restraints_to_value(definition.restraints),
    );
    Value::Object(object)
}

fn support_to_value(support: &Support) -> Value {
    let mut object = Map::new();
    object.insert("node".into(), Value::from(support.node.get()));
    object.insert("definition".into(), Value::from(support.definition.get()));
    Value::Object(object)
}

fn patch_tie_to_value(tie: &PatchTie) -> Value {
    let mut object = Map::new();
    object.insert("id".into(), Value::from(tie.id.get()));
    object.insert("master".into(), Value::from(tie.master.get()));
    object.insert("masterEdge".into(), Value::from(tie.master_edge.as_str()));
    object.insert("slave".into(), Value::from(tie.slave.get()));
    object.insert("slaveEdge".into(), Value::from(tie.slave_edge.as_str()));
    Value::Object(object)
}

fn array(values: &[Value]) -> Value {
    Value::Array(values.to_vec())
}

/// Builds the model part of a v1 document.
pub(crate) fn model_to_value(model: &Model) -> Result<Value, InterchangeError> {
    let mut root = Map::new();
    root.insert(
        "nodes".into(),
        array(
            &model
                .nodes()
                .iter()
                .map(node_to_value)
                .collect::<Result<Vec<_>, _>>()?,
        ),
    );
    root.insert(
        "materials".into(),
        array(
            &model
                .materials()
                .iter()
                .map(material_to_value)
                .collect::<Result<Vec<_>, _>>()?,
        ),
    );
    root.insert(
        "sections".into(),
        array(
            &model
                .sections()
                .iter()
                .map(section_to_value)
                .collect::<Result<Vec<_>, _>>()?,
        ),
    );
    root.insert(
        "members".into(),
        array(
            &model
                .members()
                .iter()
                .map(member_to_value)
                .collect::<Result<Vec<_>, _>>()?,
        ),
    );
    root.insert(
        "areas".into(),
        array(
            &model
                .areas()
                .iter()
                .map(area_to_value)
                .collect::<Result<Vec<_>, _>>()?,
        ),
    );
    root.insert(
        "patches".into(),
        array(
            &model
                .patches()
                .iter()
                .map(patch_to_value)
                .collect::<Result<Vec<_>, _>>()?,
        ),
    );
    root.insert(
        "loadCases".into(),
        array(
            &model
                .load_cases()
                .iter()
                .map(load_case_to_value)
                .collect::<Result<Vec<_>, _>>()?,
        ),
    );
    root.insert(
        "loads".into(),
        array(
            &model
                .loads()
                .iter()
                .map(load_to_value)
                .collect::<Result<Vec<_>, _>>()?,
        ),
    );
    root.insert(
        "combinations".into(),
        array(
            &model
                .combinations()
                .iter()
                .map(combination_to_value)
                .collect::<Result<Vec<_>, _>>()?,
        ),
    );
    root.insert(
        "supportDefinitions".into(),
        array(
            &model
                .support_definitions()
                .iter()
                .map(support_definition_to_value)
                .collect::<Vec<_>>(),
        ),
    );
    root.insert(
        "supports".into(),
        array(
            &model
                .supports()
                .iter()
                .map(support_to_value)
                .collect::<Vec<_>>(),
        ),
    );
    root.insert(
        "patchTies".into(),
        array(
            &model
                .patch_ties()
                .iter()
                .map(patch_tie_to_value)
                .collect::<Vec<_>>(),
        ),
    );
    Ok(Value::Object(root))
}

fn field<'a>(
    object: &'a Map<String, Value>,
    path: &str,
    key: &str,
) -> Result<&'a Value, InterchangeError> {
    object
        .get(key)
        .ok_or_else(|| InterchangeError::MissingField {
            path: format!("{path}.{key}"),
        })
}

fn object<'a>(value: &'a Value, path: &str) -> Result<&'a Map<String, Value>, InterchangeError> {
    value
        .as_object()
        .ok_or_else(|| InterchangeError::ExpectedType {
            path: path.to_string(),
            expected: "object",
        })
}

fn array_field<'a>(
    object: &'a Map<String, Value>,
    path: &str,
    key: &str,
) -> Result<&'a Vec<Value>, InterchangeError> {
    field(object, path, key)?
        .as_array()
        .ok_or_else(|| InterchangeError::ExpectedType {
            path: format!("{path}.{key}"),
            expected: "array",
        })
}

fn u32_field(object: &Map<String, Value>, path: &str, key: &str) -> Result<u32, InterchangeError> {
    field(object, path, key)?
        .as_u64()
        .and_then(|value| u32::try_from(value).ok())
        .ok_or_else(|| InterchangeError::ExpectedType {
            path: format!("{path}.{key}"),
            expected: "u32",
        })
}

fn f64_field(object: &Map<String, Value>, path: &str, key: &str) -> Result<f64, InterchangeError> {
    let value = field(object, path, key)?;
    let number = value
        .as_f64()
        .ok_or_else(|| InterchangeError::ExpectedType {
            path: format!("{path}.{key}"),
            expected: "number",
        })?;
    if number.is_finite() {
        Ok(number)
    } else {
        Err(InterchangeError::InvalidNumber {
            path: format!("{path}.{key}"),
        })
    }
}

fn str_field<'a>(
    object: &'a Map<String, Value>,
    path: &str,
    key: &str,
) -> Result<&'a str, InterchangeError> {
    field(object, path, key)?
        .as_str()
        .ok_or_else(|| InterchangeError::ExpectedType {
            path: format!("{path}.{key}"),
            expected: "string",
        })
}

fn bool_field(
    object: &Map<String, Value>,
    path: &str,
    key: &str,
) -> Result<bool, InterchangeError> {
    field(object, path, key)?
        .as_bool()
        .ok_or_else(|| InterchangeError::ExpectedType {
            path: format!("{path}.{key}"),
            expected: "boolean",
        })
}

fn f64_list(value: &Value, path: &str) -> Result<Vec<f64>, InterchangeError> {
    let values = value
        .as_array()
        .ok_or_else(|| InterchangeError::ExpectedType {
            path: path.to_string(),
            expected: "array",
        })?;
    values
        .iter()
        .enumerate()
        .map(|(index, value)| {
            let number = value
                .as_f64()
                .ok_or_else(|| InterchangeError::ExpectedType {
                    path: format!("{path}[{index}]"),
                    expected: "number",
                })?;
            if number.is_finite() {
                Ok(number)
            } else {
                Err(InterchangeError::InvalidNumber {
                    path: format!("{path}[{index}]"),
                })
            }
        })
        .collect()
}

fn u32_list(value: &Value, path: &str) -> Result<Vec<u32>, InterchangeError> {
    let values = value
        .as_array()
        .ok_or_else(|| InterchangeError::ExpectedType {
            path: path.to_string(),
            expected: "array",
        })?;
    values
        .iter()
        .enumerate()
        .map(|(index, value)| {
            value
                .as_u64()
                .and_then(|number| u32::try_from(number).ok())
                .ok_or_else(|| InterchangeError::ExpectedType {
                    path: format!("{path}[{index}]"),
                    expected: "u32",
                })
        })
        .collect()
}

fn matrix_field<T>(
    object: &Map<String, Value>,
    path: &str,
    key: &str,
    convert: fn(&Value, &str) -> Result<Vec<T>, InterchangeError>,
) -> Result<Vec<Vec<T>>, InterchangeError> {
    array_field(object, path, key)?
        .iter()
        .enumerate()
        .map(|(index, row)| convert(row, &format!("{path}.{key}[{index}]")))
        .collect()
}

fn add_entity(result: Result<(), ModelError>, path: &str) -> Result<(), InterchangeError> {
    result.map_err(|error| InterchangeError::InvalidEntity {
        path: path.to_string(),
        message: error.to_string(),
    })
}

fn read_node(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<Node, InterchangeError> {
    let object = object(value, path)?;
    check_keys(object, path, NODE_KEYS, mode, notes)?;
    Ok(Node::new(
        NodeId::new(u32_field(object, path, "id")?),
        f64_field(object, path, "x")?,
        f64_field(object, path, "y")?,
        f64_field(object, path, "z")?,
    ))
}

fn read_material(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<Material, InterchangeError> {
    let object = object(value, path)?;
    check_keys(object, path, MATERIAL_KEYS, mode, notes)?;
    Ok(Material::new(
        MaterialId::new(u32_field(object, path, "id")?),
        str_field(object, path, "name")?,
        f64_field(object, path, "density")?,
        f64_field(object, path, "elasticModulus")?,
        f64_field(object, path, "poissonRatio")?,
    ))
}

fn read_section(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<Section, InterchangeError> {
    let object = object(value, path)?;
    check_keys(object, path, SECTION_KEYS, mode, notes)?;
    Ok(Section::new(
        SectionId::new(u32_field(object, path, "id")?),
        str_field(object, path, "name")?,
        f64_field(object, path, "area")?,
        f64_field(object, path, "inertiaY")?,
        f64_field(object, path, "inertiaZ")?,
        f64_field(object, path, "torsionConstant")?,
        f64_field(object, path, "warpingConstant")?,
    ))
}

fn read_member(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<Member, InterchangeError> {
    let object = object(value, path)?;
    check_keys(object, path, MEMBER_KEYS, mode, notes)?;
    let spelling = str_field(object, path, "kind")?;
    let kind = ElementKind::parse(spelling).ok_or_else(|| InterchangeError::InvalidEnum {
        path: format!("{path}.kind"),
        value: spelling.to_string(),
    })?;
    Ok(Member::new(
        MemberId::new(u32_field(object, path, "id")?),
        NodeId::new(u32_field(object, path, "start")?),
        NodeId::new(u32_field(object, path, "end")?),
        SectionId::new(u32_field(object, path, "section")?),
        f64_field(object, path, "rotation")?,
        kind,
    ))
}

fn read_area(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<Area, InterchangeError> {
    let object = object(value, path)?;
    check_keys(object, path, AREA_KEYS, mode, notes)?;
    let spelling = str_field(object, path, "kind")?;
    let kind = AreaKind::parse(spelling).ok_or_else(|| InterchangeError::InvalidEnum {
        path: format!("{path}.kind"),
        value: spelling.to_string(),
    })?;
    let ids = u32_list(field(object, path, "nodes")?, &format!("{path}.nodes"))?;
    let nodes: [u32; 4] = ids.try_into().map_err(|_| InterchangeError::ExpectedType {
        path: format!("{path}.nodes"),
        expected: "four node ids",
    })?;
    Ok(Area::new(
        AreaId::new(u32_field(object, path, "id")?),
        nodes.map(NodeId::new),
        kind,
    ))
}

fn read_patch(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<SplinePatch, InterchangeError> {
    let object = object(value, path)?;
    check_keys(object, path, PATCH_KEYS, mode, notes)?;
    let net = matrix_field(object, path, "controlNet", u32_list)?;
    let control_net = net
        .into_iter()
        .map(|row| row.into_iter().map(NodeId::new).collect())
        .collect();
    Ok(SplinePatch::new(
        PatchId::new(u32_field(object, path, "id")?),
        u8_field(object, path, "degreeU")?,
        u8_field(object, path, "degreeV")?,
        f64_list(field(object, path, "knotsU")?, &format!("{path}.knotsU"))?,
        f64_list(field(object, path, "knotsV")?, &format!("{path}.knotsV"))?,
        control_net,
        matrix_field(object, path, "weights", f64_list)?,
    ))
}

fn u8_field(object: &Map<String, Value>, path: &str, key: &str) -> Result<u8, InterchangeError> {
    let value = u32_field(object, path, key)?;
    u8::try_from(value).map_err(|_| InterchangeError::ExpectedType {
        path: format!("{path}.{key}"),
        expected: "u8",
    })
}

fn read_load(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<Load, InterchangeError> {
    let object = object(value, path)?;
    let kind = str_field(object, path, "type")?;
    let allowed = match kind {
        "nodal" => NODAL_LOAD_KEYS,
        "member" => MEMBER_LOAD_KEYS,
        "area" => AREA_LOAD_KEYS,
        "patch" => PATCH_LOAD_KEYS,
        _ => {
            return Err(InterchangeError::InvalidEnum {
                path: format!("{path}.type"),
                value: kind.to_string(),
            })
        }
    };
    check_keys(object, path, allowed, mode, notes)?;
    let case = LoadCaseId::new(u32_field(object, path, "case")?);
    Ok(match kind {
        "nodal" => Load::Nodal {
            case,
            load: NodalLoad::new(
                NodeId::new(u32_field(object, path, "node")?),
                vector_field(object, path, "force")?,
                vector_field(object, path, "moment")?,
            ),
        },
        "member" => Load::Member {
            case,
            load: MemberLoad::new(
                MemberId::new(u32_field(object, path, "member")?),
                vector_field(object, path, "force")?,
                vector_field(object, path, "moment")?,
            ),
        },
        "area" => Load::Area {
            case,
            load: AreaLoad::new(
                AreaId::new(u32_field(object, path, "area")?),
                f64_field(object, path, "pressure")?,
            ),
        },
        _ => Load::Patch {
            case,
            load: PatchLoad::new(
                PatchId::new(u32_field(object, path, "patch")?),
                f64_field(object, path, "pressure")?,
            ),
        },
    })
}

fn vector_field(
    object: &Map<String, Value>,
    path: &str,
    key: &str,
) -> Result<[f64; 3], InterchangeError> {
    let values = f64_list(field(object, path, key)?, &format!("{path}.{key}"))?;
    values
        .try_into()
        .map_err(|_| InterchangeError::ExpectedType {
            path: format!("{path}.{key}"),
            expected: "three numbers",
        })
}

fn read_load_case(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<LoadCase, InterchangeError> {
    let object = object(value, path)?;
    check_keys(object, path, LOAD_CASE_KEYS, mode, notes)?;
    let mut case = LoadCase::new(
        LoadCaseId::new(u32_field(object, path, "id")?),
        str_field(object, path, "name")?,
    );
    if bool_field(object, path, "selfWeight")? {
        case = case.with_self_weight();
    }
    Ok(case)
}

fn read_combination(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<Combination, InterchangeError> {
    let fields = object(value, path)?;
    check_keys(fields, path, COMBINATION_KEYS, mode, notes)?;
    let factors = array_field(fields, path, "factors")?
        .iter()
        .enumerate()
        .map(|(index, value)| {
            let factor_path = format!("{path}.factors[{index}]");
            let factor_object = object(value, &factor_path)?;
            check_keys(factor_object, &factor_path, FACTOR_KEYS, mode, notes)?;
            Ok((
                LoadCaseId::new(u32_field(factor_object, &factor_path, "case")?),
                f64_field(factor_object, &factor_path, "factor")?,
            ))
        })
        .collect::<Result<Vec<_>, InterchangeError>>()?;
    Ok(Combination::new(
        CombinationId::new(u32_field(fields, path, "id")?),
        str_field(fields, path, "name")?,
        factors,
    ))
}

fn read_restraints(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<Restraints, InterchangeError> {
    let object = object(value, path)?;
    check_keys(object, path, RESTRAINT_KEYS, mode, notes)?;
    Ok(Restraints::new(
        bool_field(object, path, "tx")?,
        bool_field(object, path, "ty")?,
        bool_field(object, path, "tz")?,
        bool_field(object, path, "rx")?,
        bool_field(object, path, "ry")?,
        bool_field(object, path, "rz")?,
    ))
}

fn read_support_definition(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<SupportDefinition, InterchangeError> {
    let object = object(value, path)?;
    check_keys(object, path, SUPPORT_DEFINITION_KEYS, mode, notes)?;
    let restraints = read_restraints(
        field(object, path, "restraints")?,
        &format!("{path}.restraints"),
        mode,
        notes,
    )?;
    Ok(SupportDefinition::new(
        SupportDefinitionId::new(u32_field(object, path, "id")?),
        str_field(object, path, "name")?,
        restraints,
    ))
}

fn read_support(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<Support, InterchangeError> {
    let object = object(value, path)?;
    check_keys(object, path, SUPPORT_KEYS, mode, notes)?;
    Ok(Support::new(
        NodeId::new(u32_field(object, path, "node")?),
        SupportDefinitionId::new(u32_field(object, path, "definition")?),
    ))
}

fn read_patch_tie(
    value: &Value,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<PatchTie, InterchangeError> {
    let object = object(value, path)?;
    check_keys(object, path, PATCH_TIE_KEYS, mode, notes)?;
    let master = str_field(object, path, "masterEdge")?;
    let slave = str_field(object, path, "slaveEdge")?;
    Ok(PatchTie::new(
        PatchTieId::new(u32_field(object, path, "id")?),
        PatchId::new(u32_field(object, path, "master")?),
        PatchEdge::parse(master).ok_or_else(|| InterchangeError::InvalidEnum {
            path: format!("{path}.masterEdge"),
            value: master.to_string(),
        })?,
        PatchId::new(u32_field(object, path, "slave")?),
        PatchEdge::parse(slave).ok_or_else(|| InterchangeError::InvalidEnum {
            path: format!("{path}.slaveEdge"),
            value: slave.to_string(),
        })?,
    ))
}

fn each<'a>(
    root: &'a Map<String, Value>,
    key: &str,
    path: &str,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<Vec<&'a Value>, InterchangeError> {
    let field_path = format!("{path}.{key}");
    match root.get(key) {
        None => match mode {
            ImportMode::Strict => Err(InterchangeError::MissingField { path: field_path }),
            ImportMode::Permissive => {
                notes.push(ImportNote::new(
                    ImportNoteKind::MissingCollection,
                    field_path,
                    "absent collection read as empty",
                ));
                Ok(Vec::new())
            }
        },
        Some(value) => value
            .as_array()
            .ok_or(InterchangeError::ExpectedType {
                path: field_path,
                expected: "array",
            })
            .map(|values| values.iter().collect()),
    }
}

/// Builds a model from the model part of a v1 document.
pub(crate) fn model_from_value(
    root: &Map<String, Value>,
    mode: ImportMode,
    notes: &mut Vec<ImportNote>,
) -> Result<Model, InterchangeError> {
    let mut model = Model::new();
    for value in each(root, "nodes", "$", mode, notes)? {
        let node = read_node(value, "$.nodes[]", mode, notes)?;
        add_entity(model.add_node(node), "$.nodes[]")?;
    }
    for value in each(root, "materials", "$", mode, notes)? {
        let material = read_material(value, "$.materials[]", mode, notes)?;
        add_entity(model.add_material(material), "$.materials[]")?;
    }
    for value in each(root, "sections", "$", mode, notes)? {
        let section = read_section(value, "$.sections[]", mode, notes)?;
        add_entity(model.add_section(section), "$.sections[]")?;
    }
    for value in each(root, "members", "$", mode, notes)? {
        let member = read_member(value, "$.members[]", mode, notes)?;
        add_entity(model.add_member(member), "$.members[]")?;
    }
    for value in each(root, "areas", "$", mode, notes)? {
        let area = read_area(value, "$.areas[]", mode, notes)?;
        add_entity(model.add_area(area), "$.areas[]")?;
    }
    for value in each(root, "patches", "$", mode, notes)? {
        let patch = read_patch(value, "$.patches[]", mode, notes)?;
        add_entity(model.add_patch(patch), "$.patches[]")?;
    }
    for value in each(root, "loadCases", "$", mode, notes)? {
        let case = read_load_case(value, "$.loadCases[]", mode, notes)?;
        add_entity(model.add_load_case(case), "$.loadCases[]")?;
    }
    for value in each(root, "loads", "$", mode, notes)? {
        let load = read_load(value, "$.loads[]", mode, notes)?;
        add_entity(model.add_load(load), "$.loads[]")?;
    }
    for value in each(root, "combinations", "$", mode, notes)? {
        let combination = read_combination(value, "$.combinations[]", mode, notes)?;
        add_entity(model.add_combination(combination), "$.combinations[]")?;
    }
    for value in each(root, "supportDefinitions", "$", mode, notes)? {
        let definition = read_support_definition(value, "$.supportDefinitions[]", mode, notes)?;
        add_entity(
            model.add_support_definition(definition),
            "$.supportDefinitions[]",
        )?;
    }
    for value in each(root, "supports", "$", mode, notes)? {
        let support = read_support(value, "$.supports[]", mode, notes)?;
        add_entity(model.add_support(support), "$.supports[]")?;
    }
    for value in each(root, "patchTies", "$", mode, notes)? {
        let tie = read_patch_tie(value, "$.patchTies[]", mode, notes)?;
        add_entity(model.add_patch_tie(tie), "$.patchTies[]")?;
    }
    model
        .validate()
        .map_err(|error| InterchangeError::ModelValidation(error.to_string()))?;
    model
        .apply_supports()
        .map_err(|error| InterchangeError::ModelValidation(error.to_string()))?;
    Ok(model)
}
