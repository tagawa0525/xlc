//! 炉心の入力ファイル（TOML）。変わる頻度で4つのファイルに分ける（ADR-0105）。実行するのはケースのファイルで、
//! ケースがサイクルを、サイクルが炉心の形と台帳を指す。パスは、それを書いたファイルのあるディレクトリから辿る。
//!
//! 炉心の形（プラントでほとんど変わらない）:
//!
//! ```toml
//! pitch = 21.5                 # 集合体のピッチ (cm)
//! height = 406.3               # 炉心の高さ (cm)
//! map = ["*  R"]               # 上の行が +y。`*` は集合体を装荷する位置、ほかは反射体の種類（`@90` で回す）
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
//! loading = ["J01@90  ."]
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
//! # 炉心の状態（省略すると、軸のないテーブルだけを使う。ADR-0097）
//! [state]
//! boron = 1291        # ホウ素濃度 (ppm)
//! b10 = 0.199         # ホウ素の中の B-10 の原子分率（テーブルと同じ値）
//! moderator = 565     # 減速材の温度 (K)
//! fuel = 565          # 燃料の温度 (K)
//! ```

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Component, Path, PathBuf};

use serde::Deserialize;
use serde::de::DeserializeOwned;
use xlc_format::fewgroup::{FewGroup, Pins};

use crate::nem::{Boundary, Mesh, Problem, Solution};
use crate::state::State;

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
            layers: entry
                .layers
                .iter()
                .map(|l| normalized(&base.join(l)))
                .collect(),
        })
    }

    fn height(&self) -> f64 {
        self.bounds.last().copied().unwrap_or(0.0)
    }
}

/// `.` と、前の名前を打ち消す `..` を字面で取り除く（同じファイルを別の綴りで書いても1つにまとめるため）。
/// シンボリックリンクはたどらない。
fn normalized(path: &Path) -> PathBuf {
    let mut out = PathBuf::new();
    for component in path.components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir
                if matches!(out.components().next_back(), Some(Component::Normal(_))) =>
            {
                out.pop();
            }
            other => out.push(other),
        }
    }
    out
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
                        _ if cell.starts_with("*@") => invalid(format!(
                            "map: {cell:?}: rotate assemblies in the cycle's loading, not at `*`"
                        )),
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
    /// 炉心の状態。書けば分岐のテーブルをこの状態で内挿し、書かなければ軸のないテーブルだけを使う（ADR-0097）
    pub state: Option<State>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct CaseFile {
    cycle: PathBuf,
    axial_mesh: Vec<f64>,
    radial_nodes: usize,
    boundary: BoundaryEntry,
    state: Option<State>,
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
            state: file.state,
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
    /// ファイルごとに、反射体（の種類のどれか）に使う層か
    pub reflector: Vec<bool>,
    pub materials: Vec<LayoutMaterial>,
    /// ノードごとの物質の番号。並びは x が最も速く、次に y（下から）、z
    pub node_material: Vec<usize>,
    /// 上の行（+y）から、位置ごとの集合体の ID（反射体の位置は `None`）
    pub ids: Vec<Vec<Option<String>>>,
    pub boundaries: [Boundary; 6],
}

/// 装荷図の位置に置いたもの。
#[derive(Clone, Copy)]
struct Placement<'a> {
    column: &'a Column,
    /// 反時計回りの 90° の回数
    turns: u8,
    /// 反射体の種類か
    reflector: bool,
}

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
        // 置いた設計と反射体の種類（名前ごとに1回だけ検査する）
        let mut used: BTreeMap<String, &Column> = BTreeMap::new();
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
                        used.insert(format!("design {design}"), column);
                        row.push(Placement {
                            column,
                            turns: p.turns,
                            reflector: false,
                        });
                        id_row.push(Some(p.id.clone()));
                    }
                    (Site::Assembly, None) => {
                        return invalid(format!("{at}: an assembly ID is needed at `*`"));
                    }
                    (Site::Reflector { name, turns }, None) => {
                        let column = plant.reflectors.get(name).ok_or_else(|| {
                            InputError::Invalid(format!("{at}: unknown reflector {name:?}"))
                        })?;
                        used.insert(format!("reflector {name}"), column);
                        row.push(Placement {
                            column,
                            turns: *turns,
                            reflector: true,
                        });
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

        // 層の最後が炉心の高さで、層の境界がメッシュに乗っていること（同じ数で書く）
        let mesh = &case.axial_mesh;
        for (what, column) in &used {
            if column.height() != plant.height {
                return invalid(format!(
                    "{what}: the last bound {} differs from the height {}",
                    column.height(),
                    plant.height
                ));
            }
            if let Some(b) = column.bounds.iter().find(|b| !mesh.contains(b)) {
                return invalid(format!("{what}: layer bound {b} is not on the axial-mesh"));
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
            for placement in placements.iter().rev().flatten() {
                let Placement {
                    column,
                    turns,
                    reflector: is_reflector,
                } = *placement;
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
                    reflector.push(false);
                    files.len() - 1
                });
                // 軸方向の反射体の層などは、反射体と集合体の設計で同じファイルを使ってよい
                reflector[file] |= is_reflector;
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
                .map(|m| {
                    rotated(&constants[m.file], m.turns).map_err(|e| {
                        InputError::Invalid(format!("{}: {e}", self.files[m.file].display()))
                    })
                })
                .collect::<Result<_, _>>()?,
            node_material: self.node_material.clone(),
            boundaries: self.boundaries,
        })
    }

    /// 集合体の位置ごとの、軸方向に積分した核分裂率（集合体の位置の平均で割る）。上の行（+y）から、反射体の
    /// 位置は `None`。`problem` はこの並びから作った問題、`solution` はその解
    pub fn assembly_fission_rates(
        &self,
        problem: &Problem,
        solution: &Solution,
    ) -> Vec<Vec<Option<f64>>> {
        let (nx, ny) = (problem.mesh.dx.len(), problem.mesh.dy.len());
        let mut rates = vec![vec![0.0; nx]; ny];
        for (n, flux) in solution.flux.iter().enumerate() {
            let (x, y, z) = (n % nx, (n / nx) % ny, n / (nx * ny));
            let m = &problem.materials[problem.node_material[n]];
            let volume = problem.mesh.dx[x] * problem.mesh.dy[y] * problem.mesh.dz[z];
            let fission: f64 = m.fission.iter().zip(flux).map(|(f, phi)| f * phi).sum();
            // 図の上の行が +y
            rates[ny - 1 - y][x] += fission * volume;
        }
        let loaded: Vec<f64> = self
            .ids
            .iter()
            .zip(&rates)
            .flat_map(|(ids, r)| {
                ids.iter()
                    .zip(r)
                    .filter(|(id, _)| id.is_some())
                    .map(|(_, &r)| r)
            })
            .collect();
        let mean = loaded.iter().sum::<f64>() / loaded.len().max(1) as f64;
        self.ids
            .iter()
            .zip(&rates)
            .map(|(ids, r)| {
                ids.iter()
                    .zip(r)
                    .map(|(id, &r)| id.as_ref().map(|_| r / mean))
                    .collect()
            })
            .collect()
    }
}

/// 少数群定数の辺ごとの量（ADF、真空の辺の比、f* の印）と、ピンの形状（ADR-0096）を、反時計回りに 90° × `turns`
/// 回す。辺の並びは（左、右、上、下）。90° 回すと、右の辺が上に、上が左に、左が下に、下が右に来る。
/// 配列の長さは `FewGroup::validate` で先に確かめる。
pub fn rotated(constants: &FewGroup, turns: u8) -> Result<FewGroup, String> {
    fn turn<T: Clone>(sides: &[T]) -> Vec<T> {
        // 新しい（左、右、上、下）は、元の（上、下、右、左）
        [2, 3, 1, 0]
            .iter()
            .filter_map(|&s| sides.get(s).cloned())
            .collect()
    }
    constants.validate().map_err(|e| e.to_string())?;
    let mut out = constants.clone();
    for _ in 0..turns % 4 {
        out.adf = turn(&out.adf);
        out.vacuum_ratio = turn(&out.vacuum_ratio);
        out.normalized_adf = turn(&out.normalized_adf);
        out.pins = out.pins.as_ref().map(turn_pins);
    }
    Ok(out)
}

/// 90° 回したピンの形状。(x, y) は (H - y, x) に移る（H はノードの高さ）ので、新しい列の幅は元の行の高さを
/// 上から並べたもの、新しい行の高さは元の列の幅。新しい行 j、列 i のセルは、元の行（下から）rows - 1 - i、列 j の
/// セル。角（-x -y、+x -y、+x +y、-x +y）は1つずつ先に進む。長さは検査済みで、並べ替えだけなので丸めは入らない
fn turn_pins(pins: &Pins) -> Pins {
    let rows = pins.dy.len();
    Pins {
        dx: pins.dy.iter().rev().copied().collect(),
        dy: pins.dx.clone(),
        form: pins
            .form
            .iter()
            .map(|g| {
                (0..pins.dx.len())
                    .map(|j| {
                        (0..rows)
                            .filter_map(|i| g.get(rows - 1 - i).and_then(|row| row.get(j)).copied())
                            .collect()
                    })
                    .collect()
            })
            .collect(),
        corners: [3, 0, 1, 2]
            .iter()
            .filter_map(|&c| pins.corners.get(c).cloned())
            .collect(),
    }
}

#[cfg(test)]
mod tests;
