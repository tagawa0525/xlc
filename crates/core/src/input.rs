//! 炉心の入力ファイル（TOML）。変わる頻度で4つのファイルに分ける（ADR-0105）。実行するのはケースのファイルで、
//! ケースがサイクルを、サイクルが炉心の形と台帳を指す。パスは、それを書いたファイルのあるディレクトリから辿る。
//!
//! 炉心の形（プラントでほとんど変わらない）:
//!
//! ```toml
//! pitch = 21.5                 # 集合体のピッチ (cm)
//! height = 406.3               # 炉心の高さ (cm)
//! map = [                      # 上の行が +y。`*` は集合体を装荷する位置、ほかは反射体の種類（`@90` で回す）
//!   "R  R  R",
//!   "R  *  R",
//!   "R  R  R",
//! ]
//! [reflector.R]
//! bounds = [0.0, 406.3]        # 軸方向の層の境界 (cm)。0 から炉心の高さまで
//! layers = ["reflector.msgpack"]
//! ```
//!
//! 台帳（サイクルをまたいで追記する）:
//!
//! ```toml
//! [design.A]
//! bounds = [0.0, 11.951, 377.711, 406.3]
//! layers = ["bottom.msgpack", "fuel-a.msgpack", "top.msgpack"]
//! [assembly.J01]
//! design = "A"
//! ```
//!
//! サイクル（`loading` は炉心の形の `map` と同じ行と列で、`*` の位置に ID、ほかに `.`。`@90` は反時計回りの回転）:
//!
//! ```toml
//! plant = "plant.toml"
//! inventory = "assemblies.toml"
//! loading = [".  .  .", ".  J01@90  .", ".  .  ."]
//! ```
//!
//! ケース（計算ごとに変わる）:
//!
//! ```toml
//! cycle = "cycle-1.toml"
//! axial-mesh = [0.0, 11.951, 200.0, 377.711, 406.3]   # 軸方向のノードの境界 (cm)
//! radial-nodes = 1                                     # 集合体1体の各方向のノードの数
//! # 各軸の（低い側, 高い側）の境界条件（reflective、zero-flux、vacuum、marshak のいずれか）。6面とも書く。
//! # vacuum は少数群定数の真空の辺の比 J/φ を使い、marshak は P1 の関係（J/φ = 1/2）を置く（ADR-0101）
//! boundary = { x = ["vacuum", "vacuum"], y = ["vacuum", "vacuum"], z = ["vacuum", "vacuum"] }
//! ```

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};

use serde::Deserialize;
use serde::de::DeserializeOwned;
use xlc_format::fewgroup::FewGroup;

use crate::nem::{Boundary, Mesh, Problem};

/// 入力ファイルの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum InputError {
    #[error(transparent)]
    Toml(#[from] toml::de::Error),
    #[error("{0}")]
    Invalid(String),
}

fn invalid<T>(message: impl Into<String>) -> Result<T, InputError> {
    Err(InputError::Invalid(message.into()))
}

fn parse_toml<T: DeserializeOwned>(text: &str) -> Result<T, InputError> {
    Ok(toml::from_str(text)?)
}

/// 軸方向の層: 境界の高さ (cm) と、層ごとの少数群定数のファイル。
#[derive(Debug, Clone, PartialEq)]
pub struct Column {
    pub bounds: Vec<f64>,
    pub layers: Vec<PathBuf>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ColumnEntry {
    bounds: Vec<f64>,
    layers: Vec<PathBuf>,
}

impl Column {
    fn new(what: &str, entry: ColumnEntry, base: &Path) -> Result<Self, InputError> {
        check_heights(what, &entry.bounds)?;
        if entry.layers.len() + 1 != entry.bounds.len() {
            return invalid(format!(
                "{what}: {} layers for {} bounds (need one fewer layer than bounds)",
                entry.layers.len(),
                entry.bounds.len()
            ));
        }
        Ok(Self {
            bounds: entry.bounds,
            layers: entry.layers.iter().map(|l| base.join(l)).collect(),
        })
    }

    fn height(&self) -> f64 {
        self.bounds.last().copied().unwrap_or(0.0)
    }
}

/// 高さの列が 0 から始まり、有限で、増えていくこと。
fn check_heights(what: &str, heights: &[f64]) -> Result<(), InputError> {
    if heights.len() < 2 || heights[0] != 0.0 {
        return invalid(format!(
            "{what}: heights must start at 0 and have at least two values, got {heights:?}"
        ));
    }
    if !heights.iter().all(|h| h.is_finite()) || heights.windows(2).any(|w| w[1] <= w[0]) {
        return invalid(format!("{what}: heights must increase, got {heights:?}"));
    }
    Ok(())
}

/// 名前や ID が、図の記号（`*`、`-`、`.`、空白、回転の `@`）と区別できること。
fn check_name(what: &str, name: &str) -> Result<(), InputError> {
    if name.is_empty()
        || matches!(name, "*" | "-" | ".")
        || name.contains('@')
        || name.chars().any(char::is_whitespace)
    {
        return invalid(format!(
            "{what} {name:?}: must not be empty, `*`, `-` or `.`, nor contain `@` or spaces"
        ));
    }
    Ok(())
}

/// 図の行を空白で区切り、行の長さがそろっていることを確かめる。
fn split_rows<'a>(what: &str, rows: &'a [String]) -> Result<Vec<Vec<&'a str>>, InputError> {
    let cells: Vec<Vec<&str>> = rows
        .iter()
        .map(|r| r.split_whitespace().collect())
        .collect();
    let Some(width) = cells.first().map(Vec::len) else {
        return invalid(format!("{what}: no rows"));
    };
    if width == 0 || cells.iter().any(|r| r.len() != width) {
        return invalid(format!(
            "{what}: every row needs the same positive number of entries"
        ));
    }
    Ok(cells)
}

/// 炉心の形の位置。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Site {
    /// 集合体を装荷する位置（`*`）
    Assembly,
    /// 反射体の種類と、反時計回りの 90° の回数（燃料に面する辺の向きを合わせる）
    Reflector { name: String, turns: u8 },
}

/// 炉心の形（ADR-0105）。
#[derive(Debug, Clone, PartialEq)]
pub struct Plant {
    pub pitch: f64,
    pub height: f64,
    pub reflectors: BTreeMap<String, Column>,
    /// 上の行（+y）から
    pub map: Vec<Vec<Site>>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct PlantFile {
    pitch: f64,
    height: f64,
    map: Vec<String>,
    #[serde(default)]
    reflector: BTreeMap<String, ColumnEntry>,
}

impl Plant {
    pub fn parse(text: &str, base: &Path) -> Result<Self, InputError> {
        let file: PlantFile = parse_toml(text)?;
        for (what, value) in [("pitch", file.pitch), ("height", file.height)] {
            if !(value > 0.0 && value.is_finite()) {
                return invalid(format!("{what} must be positive, got {value}"));
            }
        }
        let mut reflectors = BTreeMap::new();
        for (name, entry) in file.reflector {
            check_name("reflector", &name)?;
            let column = Column::new(&format!("reflector {name}"), entry, base)?;
            if column.height() != file.height {
                return invalid(format!(
                    "reflector {name}: the last bound {} differs from the height {}",
                    column.height(),
                    file.height
                ));
            }
            reflectors.insert(name, column);
        }
        let map = split_rows("map", &file.map)?
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|cell| match cell {
                        "*" => Ok(Site::Assembly),
                        "-" => invalid(
                            "map: `-` (no node) needs boundaries along the core outline (#0034)",
                        ),
                        _ => {
                            let (name, turns) = rotation("map", cell)?;
                            if !reflectors.contains_key(name) {
                                return invalid(format!("map: unknown reflector {name:?}"));
                            }
                            Ok(Site::Reflector {
                                name: name.to_string(),
                                turns,
                            })
                        }
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            pitch: file.pitch,
            height: file.height,
            reflectors,
            map,
        })
    }
}

/// 集合体の台帳（ADR-0105）。
#[derive(Debug, Clone, PartialEq)]
pub struct Inventory {
    pub designs: BTreeMap<String, Column>,
    /// ID → 設計の名前
    pub assemblies: BTreeMap<String, String>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InventoryFile {
    design: BTreeMap<String, ColumnEntry>,
    assembly: BTreeMap<String, AssemblyEntry>,
}

/// 集合体ごとの項目は、使う機能を作るときにここに足す
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct AssemblyEntry {
    design: String,
}

impl Inventory {
    pub fn parse(text: &str, base: &Path) -> Result<Self, InputError> {
        let file: InventoryFile = parse_toml(text)?;
        let mut designs = BTreeMap::new();
        for (name, entry) in file.design {
            let column = Column::new(&format!("design {name}"), entry, base)?;
            designs.insert(name, column);
        }
        let mut assemblies = BTreeMap::new();
        for (id, entry) in file.assembly {
            check_name("assembly ID", &id)?;
            if !designs.contains_key(&entry.design) {
                return invalid(format!("assembly {id}: unknown design {:?}", entry.design));
            }
            assemblies.insert(id, entry.design);
        }
        Ok(Self {
            designs,
            assemblies,
        })
    }
}

/// サイクルの装荷図の位置に置いた集合体。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Placed {
    pub id: String,
    /// 反時計回りの 90° の回数（0〜3）
    pub turns: u8,
}

/// サイクルの装荷（ADR-0105）。
#[derive(Debug, Clone, PartialEq)]
pub struct Cycle {
    pub plant: PathBuf,
    pub inventory: PathBuf,
    /// 上の行（+y）から。`.` の位置は `None`
    pub loading: Vec<Vec<Option<Placed>>>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CycleFile {
    plant: PathBuf,
    inventory: PathBuf,
    loading: Vec<String>,
}

/// `名前@角度` を名前と反時計回りの 90° の回数に分ける（角度がなければ 0）。
fn rotation<'a>(what: &str, cell: &'a str) -> Result<(&'a str, u8), InputError> {
    let Some((name, angle)) = cell.split_once('@') else {
        return Ok((cell, 0));
    };
    let turns = match angle {
        "0" => 0,
        "90" => 1,
        "180" => 2,
        "270" => 3,
        _ => {
            return invalid(format!(
                "{what} {cell:?}: the rotation must be 0, 90, 180 or 270"
            ));
        }
    };
    Ok((name, turns))
}

fn placed(cell: &str) -> Result<Option<Placed>, InputError> {
    if cell == "." {
        return Ok(None);
    }
    let (id, turns) = rotation("loading", cell)?;
    check_name("loading", id)?;
    Ok(Some(Placed {
        id: id.to_string(),
        turns,
    }))
}

impl Cycle {
    pub fn parse(text: &str, base: &Path) -> Result<Self, InputError> {
        let file: CycleFile = parse_toml(text)?;
        let loading = split_rows("loading", &file.loading)?
            .into_iter()
            .map(|row| row.into_iter().map(placed).collect::<Result<Vec<_>, _>>())
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            plant: base.join(file.plant),
            inventory: base.join(file.inventory),
            loading,
        })
    }
}

/// 計算のケース（ADR-0105）。
#[derive(Debug, Clone, PartialEq)]
pub struct Case {
    pub cycle: PathBuf,
    pub axial_mesh: Vec<f64>,
    pub radial_nodes: usize,
    /// [-x, +x, -y, +y, -z, +z]
    pub boundaries: [Boundary; 6],
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct CaseFile {
    cycle: PathBuf,
    axial_mesh: Vec<f64>,
    radial_nodes: usize,
    boundary: BoundaryEntry,
}

/// 各軸の（低い側, 高い側）
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BoundaryEntry {
    x: [Condition; 2],
    y: [Condition; 2],
    z: [Condition; 2],
}

#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "kebab-case")]
enum Condition {
    Reflective,
    ZeroFlux,
    Vacuum,
    Marshak,
}

impl Case {
    pub fn parse(text: &str, base: &Path) -> Result<Self, InputError> {
        let file: CaseFile = parse_toml(text)?;
        check_heights("axial-mesh", &file.axial_mesh)?;
        match file.radial_nodes {
            0 => return invalid("radial-nodes must be positive"),
            1 => {}
            n => {
                return invalid(format!(
                    "radial-nodes = {n}: ADFs of split nodes are not defined yet (#0033)"
                ));
            }
        }
        let b = file.boundary;
        let boundaries = [b.x[0], b.x[1], b.y[0], b.y[1], b.z[0], b.z[1]].map(|c| match c {
            Condition::Reflective => Boundary::Reflective,
            Condition::ZeroFlux => Boundary::ZeroFlux,
            Condition::Vacuum => Boundary::Vacuum,
            Condition::Marshak => Boundary::Marshak,
        });
        Ok(Self {
            cycle: base.join(file.cycle),
            axial_mesh: file.axial_mesh,
            radial_nodes: file.radial_nodes,
            boundaries,
        })
    }
}

/// ノードの物質: 少数群定数のファイルの番号と、反時計回りの 90° の回数。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LayoutMaterial {
    pub file: usize,
    pub turns: u8,
}

/// 4つの入力から組み立てた炉心のノードの並び。
#[derive(Debug, Clone, PartialEq)]
pub struct Layout {
    pub mesh: Mesh,
    /// 少数群定数のファイル（現れた順）
    pub files: Vec<PathBuf>,
    /// ファイルごとに、反射体の層か
    pub reflector: Vec<bool>,
    pub materials: Vec<LayoutMaterial>,
    /// ノードごとの物質の番号。並びは x が最も速く、次に y（下から）、z
    pub node_material: Vec<usize>,
    /// 上の行（+y）から、位置ごとの集合体の ID（反射体の位置は `None`）
    pub ids: Vec<Vec<Option<String>>>,
    pub boundaries: [Boundary; 6],
}

/// 位置ごとの層、回転、反射体か
type Placement<'a> = (&'a Column, u8, bool);

impl Layout {
    pub fn new(
        case: &Case,
        cycle: &Cycle,
        plant: &Plant,
        inventory: &Inventory,
    ) -> Result<Self, InputError> {
        let mesh_height = case.axial_mesh.last().copied().unwrap_or(0.0);
        if mesh_height != plant.height {
            return invalid(format!(
                "axial-mesh ends at {mesh_height}, but the core height is {}",
                plant.height
            ));
        }
        let rows = plant.map.len();
        let cols = plant.map.first().map_or(0, Vec::len);
        if cycle.loading.len() != rows || cycle.loading.iter().any(|r| r.len() != cols) {
            return invalid(format!(
                "loading must have the {rows} rows and {cols} columns of the plant map"
            ));
        }

        let mut seen = HashSet::new();
        let mut placements: Vec<Vec<Placement>> = Vec::with_capacity(rows);
        let mut ids = Vec::with_capacity(rows);
        for (r, (sites, loaded)) in plant.map.iter().zip(&cycle.loading).enumerate() {
            let mut row = Vec::with_capacity(cols);
            let mut id_row = Vec::with_capacity(cols);
            for (c, (site, placed)) in sites.iter().zip(loaded).enumerate() {
                let at = format!("loading row {} column {}", r + 1, c + 1);
                match (site, placed) {
                    (Site::Assembly, Some(p)) => {
                        let design = inventory.assemblies.get(&p.id).ok_or_else(|| {
                            InputError::Invalid(format!("{at}: {} is not in the inventory", p.id))
                        })?;
                        if !seen.insert(p.id.as_str()) {
                            return invalid(format!("{at}: {} is loaded twice", p.id));
                        }
                        let column = inventory.designs.get(design).ok_or_else(|| {
                            InputError::Invalid(format!("{at}: unknown design {design:?}"))
                        })?;
                        if column.height() != plant.height {
                            return invalid(format!(
                                "design {design}: the last bound {} differs from the height {}",
                                column.height(),
                                plant.height
                            ));
                        }
                        row.push((column, p.turns, false));
                        id_row.push(Some(p.id.clone()));
                    }
                    (Site::Assembly, None) => {
                        return invalid(format!("{at}: an assembly ID is needed at `*`"));
                    }
                    (Site::Reflector { name, turns }, None) => {
                        let column = plant.reflectors.get(name).ok_or_else(|| {
                            InputError::Invalid(format!("{at}: unknown reflector {name:?}"))
                        })?;
                        row.push((column, *turns, true));
                        id_row.push(None);
                    }
                    (Site::Reflector { name, .. }, Some(p)) => {
                        return invalid(format!(
                            "{at}: {} is placed on the reflector {name}; use `.`",
                            p.id
                        ));
                    }
                }
            }
            placements.push(row);
            ids.push(id_row);
        }

        // 層の境界がメッシュに乗っていること（同じ数で書く）
        let mesh = &case.axial_mesh;
        for (column, _, _) in placements.iter().flatten() {
            if let Some(b) = column.bounds.iter().find(|b| !mesh.contains(b)) {
                return invalid(format!("layer bound {b} is not on the axial-mesh"));
            }
        }

        let mut files: Vec<PathBuf> = Vec::new();
        let mut reflector: Vec<bool> = Vec::new();
        let mut file_index: HashMap<&Path, usize> = HashMap::new();
        let mut materials = Vec::new();
        let mut material_index: HashMap<LayoutMaterial, usize> = HashMap::new();
        let mut node_material = Vec::with_capacity(rows * cols * (mesh.len() - 1));
        for z in mesh.windows(2) {
            // 図の上の行が +y なので、y は下の行から
            for &(column, turns, is_reflector) in placements.iter().rev().flatten() {
                // 層の境界はメッシュに乗っているので、ノードはどれか1つの層に収まる
                let layer = column
                    .bounds
                    .windows(2)
                    .position(|b| b[0] <= z[0] && z[1] <= b[1]);
                let Some(path) = layer.map(|l| column.layers[l].as_path()) else {
                    return invalid(format!("no layer covers {}..{} cm", z[0], z[1]));
                };
                let file = *file_index.entry(path).or_insert_with(|| {
                    files.push(path.to_path_buf());
                    reflector.push(is_reflector);
                    files.len() - 1
                });
                if reflector[file] != is_reflector {
                    return invalid(format!(
                        "{} is used both as a reflector and in an assembly",
                        path.display()
                    ));
                }
                let material = LayoutMaterial { file, turns };
                let m = *material_index.entry(material).or_insert_with(|| {
                    materials.push(material);
                    materials.len() - 1
                });
                node_material.push(m);
            }
        }

        Ok(Self {
            mesh: Mesh {
                dx: vec![plant.pitch; cols],
                dy: vec![plant.pitch; rows],
                dz: mesh.windows(2).map(|w| w[1] - w[0]).collect(),
            },
            files,
            reflector,
            materials,
            node_material,
            ids,
            boundaries: case.boundaries,
        })
    }

    /// `files` の順の少数群定数から、NEM の問題を作る。反射体の層が核分裂すればエラー。
    pub fn problem(&self, constants: &[FewGroup]) -> Result<Problem, InputError> {
        if constants.len() != self.files.len() {
            return invalid(format!(
                "{} few-group constants for {} files",
                constants.len(),
                self.files.len()
            ));
        }
        for ((c, path), &is_reflector) in constants.iter().zip(&self.files).zip(&self.reflector) {
            if is_reflector && c.nu_fission.iter().any(|&nf| nf != 0.0) {
                return invalid(format!(
                    "reflector layer {} has fission; load fuel with an assembly ID",
                    path.display()
                ));
            }
        }
        Ok(Problem {
            mesh: self.mesh.clone(),
            materials: self
                .materials
                .iter()
                .map(|m| rotated(&constants[m.file], m.turns))
                .collect(),
            node_material: self.node_material.clone(),
            boundaries: self.boundaries,
        })
    }
}

/// 少数群定数の辺ごとの量（ADF、真空の辺の比、f* の印）を、反時計回りに 90° × `turns` 回す。
/// 辺の並びは（左、右、上、下）。90° 回すと、右の辺が上に、上が左に、左が下に、下が右に来る。
pub fn rotated(constants: &FewGroup, turns: u8) -> FewGroup {
    fn turn<T: Clone>(sides: &[T]) -> Vec<T> {
        // 新しい（左、右、上、下）は、元の（上、下、右、左）
        [2, 3, 1, 0]
            .iter()
            .filter_map(|&s| sides.get(s).cloned())
            .collect()
    }
    let mut out = constants.clone();
    for _ in 0..turns % 4 {
        out.adf = turn(&out.adf);
        out.vacuum_ratio = turn(&out.vacuum_ratio);
        out.normalized_adf = turn(&out.normalized_adf);
    }
    out
}

#[cfg(test)]
mod tests;
