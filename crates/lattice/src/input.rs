//! 格子計算の入力ファイル（TOML）。
//!
//! ```toml
//! library = "m1.msgpack"          # 多群ライブラリ（cross-sections と排他）
//! self-shielding = false          # 燃料の自己遮蔽（library のとき。省略すると true）
//! # cross-sections = "c5g7.toml"  # 巨視的断面積のファイル（このときは material を書かない）
//! few-group = { output = "pincell.fg.msgpack", bounds = [1e7, 0.625, 1e-5] }   # 少数群定数（library のとき）
//! # pins = true を足すと、ピン出力の再構成の形状関数と角の不連続因子も書く（output では外周がすべて鏡面反射で、
//! # 角に余白のある格子。ADR-0096。nodes ではどの格子でも、ノードが格子の上の長方形の並びであること。ADR-0131）。
//! # split = n を足すと、格子全体を n × n の部分ノードに分けた定数を書く
//! # （外周がすべて鏡面反射の格子。nodes とは使えない。ADR-0108）
//! # 格子をノードに分けるときは output の代わりに nodes を書く（ADR-0095）。列と行は格子のセル（余白の薄い行と列を
//! # 含む）の 1 から数えた最初と最後。ノードは格子を重ならずに覆う。reflector は反射体のノード（燃料と接する辺の
//! # ADF と、燃料のノードに接する角の点の CDF を燃料の値で割る。燃料に角の点だけで接するノードにも付ける。
//! # 省略すると false）
//! # few-group = { bounds = [1e7, 0.625, 1e-5], nodes = [
//! #   { output = "fuel.msgpack", columns = [1, 19], rows = [1, 19] },
//! #   { output = "reflector.msgpack", columns = [20, 30], rows = [1, 19], reflector = true } ] }
//! pin-powers = "pins.msgpack"     # セルごとの核分裂率（平均が 1）
//! medium = "fuel.medium"          # 格子全体を細かい群のまま均質化した媒質（library のとき。外周がすべて鏡面反射の
//!                                 # 格子。軸方向の反射体の平板の燃料にする。ADR-0128）
//! pitch = 1.26                    # セルのピッチ (cm)
//! azimuthal = 32                  # (0, π) の方位角の数
//! spacing = 0.05                  # トラック間隔 (cm)
//! scattering-order = 1            # 散乱の非等方性のルジャンドルの次数（library のとき。省略するとライブラリにあるすべての次数）
//! symmetry = "quarter"            # 使う対称性を強制する（検証用。"none"、"quarter"（左右と上下の鏡映）、"octant"（8 つの
//!                                 # 操作すべて））。どれかの操作で格子が変わればエラー。省略すると、成り立つ対称操作を
//!                                 # すべて使う（ADR-0112）
//! boundary = { right = "vacuum", bottom = { albedo = 0.3 } }   # 外周の辺の境界条件（reflective, vacuum,
//!                                                              # { albedo = α }。書かない辺は鏡面反射）
//! rows = ["pin box", "box pin"]   # セルの配置を上の行から
//! gap = { width = 0.04, material = "water" }   # 格子の周りの余白（集合体の間の水のギャップ。省略すると無し）
//!
//! [material.fuel]
//! temperature = 900               # K。ライブラリにある温度であること
//! nuclides = { U-235 = 7.0e-4, U-238 = 2.2e-2, O-16 = 4.6e-2 }   # 個数密度 (1/(barn cm))
//!
//! [material.water]
//! temperature = 565
//! # nuclides の代わりに、ホウ酸水の密度 (g/cm³)、ホウ素の質量 ppm、B-10 の原子分率（ADR-0106）
//! water = { density = 0.743, boron = 1291, b10 = 0.199 }
//!
//! # 軸方向の反射体の平板の物質（ADR-0128。どちらも temperature を書かず、self-shielding = false の格子で使う）
//! [material.assembly]
//! medium = "fuel.medium"          # 単一集合体の格子計算が書いた細かい群の媒質
//! [material.plate]                # ほかの物質（核種かホウ酸水で書いたもの、または void）を体積の割合で混ぜる
//! mix = [{ material = "ss304", volume = 0.5 }, { material = "water", volume = 0.5 }]
//!
//! [cell.pin]                      # 帯を内側から。最後の帯は radius を省略してセルの外周まで
//! zones = [{ material = "fuel", radius = 0.4096, rings = 3 }, { material = "void", radius = 0.418, rings = 1 },
//!          { material = "water", rings = 2 }]   # void は真空（定義しない予約の名前）
//! sectors = 8                     # 各リングを 45° ごとの扇形に分ける（1, 4, 8。省略すると 1）
//! frame = { thickness = 0.0141, material = "zircaloy" }   # 外周の正方形の枠（支持格子。省略すると無し）
//!
//! [branches]                    # 分岐（省略すると1つの状態）。軸と書き方は lattice::branches（ADR-0097）
//! fuel = { temperature = [550, 565, 580], materials = ["fuel"] }
//!
//! [cell.box]
//! zones = [{ material = "water", rings = 1 }]
//! grid = 4                        # 均質なセルを 4 × 4 の正方形に分ける（省略すると 1）
//! ```
//!
//! セルの帯は外半径の内側を分割数の領域に分ける。外半径のある帯は等体積に（切り取る前の円の
//! 面積で）、最後の帯はセルの内接円までを半径方向に等間隔に分ける。外半径は外接円（ピッチ/√2）
//! まで書け、内接円より大きい円はセルの外周で切り取る。四隅の減速材を半径方向に分けるには、
//! 内接円より大きい外半径の帯を置き、最後の帯は分割数 1 にする。

use std::path::PathBuf;

use std::borrow::Cow;
use std::collections::BTreeMap;

use serde::Deserialize;
use xlc_format::library::is_nuclide_name;

use crate::branches::{Branches, BranchesEntry, Water};
use crate::fewgroup::{Node, node_grid};
use crate::geometry::{Cell, GeometryError, Lattice, PinCell};
use crate::symmetry::Symmetry;
use crate::tracking::{Boundaries, Boundary};
use crate::water::borated_water;

/// 解くときに使う格子の対称性の選び方（ADR-0112）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SymmetryChoice {
    /// 成り立つ対称操作をすべて使う（`symmetry` を省略したとき）
    Detect,
    /// この対称性を使う（検証用）。どれかの操作が成り立たなければエラー
    Force(Symmetry),
}

/// 真空の帯の物質の名前（予約語。この名前の物質は定義できない）
pub const VOID: &str = "void";

/// 断面積の出どころ。
#[derive(Debug, Clone, PartialEq)]
pub enum CrossSections {
    /// 多群ライブラリ（`material` で組成を与える）
    Library(PathBuf),
    /// 巨視的断面積のファイル（物質はファイルの名前で参照する）
    Macroscopic(PathBuf),
}

/// ライブラリから作る物質。
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialInput {
    pub name: String,
    /// 温度 (K)。ライブラリにある温度であること
    pub temperature: f64,
    /// 核種と個数密度 (1/(barn cm))
    pub composition: Vec<(String, f64)>,
}

/// 細かい群の媒質のファイルを読む物質（ADR-0128）。
#[derive(Debug, Clone, PartialEq)]
pub struct MediumInput {
    pub name: String,
    /// 媒質のファイル（入力のファイルのあるディレクトリから）
    pub file: PathBuf,
}

/// ほかの物質を体積の割合で混ぜた物質（ADR-0128）。
#[derive(Debug, Clone, PartialEq)]
pub struct Mixture {
    pub name: String,
    /// 成分の物質の名前（核種かホウ酸水で書いた物質、または `VOID`）と体積の割合。割合の和は 1
    pub parts: Vec<(String, f64)>,
}

/// セルの同心円の帯。
#[derive(Debug, Clone, PartialEq)]
pub struct Zone {
    pub material: String,
    /// 外半径 (cm)。`None` はセルの外周まで
    pub outer: Option<f64>,
    /// 帯を分ける領域の数
    pub rings: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CellInput {
    pub name: String,
    pub zones: Vec<Zone>,
    /// 均質なセルを分ける小さな正方形の1辺の数（`grid M`、書かなければ 1）
    pub grid: usize,
    /// 各リングを分ける扇形の数（`sectors S`、書かなければ 1）
    pub sectors: usize,
    /// 外周の正方形の枠（支持格子を体積を保って薄い枠にしたもの）
    pub frame: Option<Frame>,
}

/// セルの外周の枠。
#[derive(Debug, Clone, PartialEq)]
pub struct Frame {
    /// 厚さ (cm)
    pub thickness: f64,
    pub material: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LatticeInput {
    pub cross_sections: CrossSections,
    /// 燃料（核分裂する核種を含む物質）をピンごとに自己遮蔽する（ライブラリのとき）
    pub self_shielding: bool,
    /// 少数群定数の出力
    pub output: Option<FewGroupOutput>,
    pub pin_powers: Option<PathBuf>,
    /// 格子全体を細かい群のまま均質化した媒質の出力（ADR-0128）
    pub medium: Option<PathBuf>,
    pub pitch: f64,
    pub azimuthal: usize,
    pub spacing: f64,
    /// 散乱の非等方性のルジャンドルの次数（0 なら等方散乱。`None` ならライブラリにあるすべての次数）
    pub scattering_order: Option<usize>,
    pub boundaries: Boundaries,
    /// 解くときに使う格子の対称性の選び方（ADR-0112）
    pub symmetry: SymmetryChoice,
    /// 核種の組成かホウ酸水で書いた物質
    pub materials: Vec<MaterialInput>,
    /// 細かい群の媒質の物質（ADR-0128）
    pub media: Vec<MediumInput>,
    /// 体積の割合で混ぜた物質（ADR-0128）
    pub mixtures: Vec<Mixture>,
    pub cells: Vec<CellInput>,
    /// 上の行から順の、セルの名前の配置
    pub rows: Vec<Vec<String>>,
    /// 列の幅（左から, cm）。`None` ならどの列もピッチ
    pub widths: Option<Vec<f64>>,
    /// 行の高さ（上から, cm）。`None` ならどの行もピッチ
    pub heights: Option<Vec<f64>>,
    /// 格子の周りの余白
    pub gap: Option<Gap>,
    /// 分岐（ADR-0097）。`None` なら物質は書いたままの1つの状態
    pub branches: Option<Branches>,
}

/// 少数群定数の出力。
#[derive(Debug, Clone, PartialEq)]
pub struct FewGroupOutput {
    /// 群境界 (eV)
    pub bounds: Vec<f64>,
    pub files: FewGroupFiles,
    /// 定数に、ピン出力の再構成の形状関数と角の不連続因子を付ける（格子全体は ADR-0096、ノードごとは ADR-0131）
    pub pins: bool,
    /// 格子全体を各方向に分ける部分ノードの数（ADR-0108）
    pub split: usize,
}

/// 少数群定数のファイル。
#[derive(Debug, Clone, PartialEq)]
pub enum FewGroupFiles {
    /// 格子全体を1つのノードとして均質化する
    Lattice(PathBuf),
    /// ノードごとのファイル（ADR-0095）
    Nodes(Vec<(PathBuf, Node)>),
}

/// 格子の周りの余白（集合体の間の水のギャップなど）。
#[derive(Debug, Clone, PartialEq)]
pub struct Gap {
    /// 幅 (cm)
    pub width: f64,
    pub material: String,
}

/// 入力ファイルの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum InputError {
    #[error(transparent)]
    Toml(#[from] toml::de::Error),
    #[error("{0}")]
    Invalid(String),
    #[error("cell `{cell}`: material `{material}` is not defined")]
    UndefinedMaterial { cell: String, material: String },
    #[error("cell `{0}` is not defined")]
    UndefinedCell(String),
    #[error(
        "cell `{cell}`: the last zone starts at {start} cm, outside the inscribed circle, \
         so it cannot be split into rings; give the zone an outer radius instead"
    )]
    LastZoneOutside { cell: String, start: f64 },
    #[error("cell `{cell}`")]
    Cell { cell: String, source: GeometryError },
    #[error(transparent)]
    Geometry(#[from] GeometryError),
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct InputFile {
    library: Option<PathBuf>,
    /// 省略すると、ライブラリのときは true
    self_shielding: Option<bool>,
    cross_sections: Option<PathBuf>,
    /// 省略すると少数群定数を出さない
    few_group: Option<FewGroupEntry>,
    /// 省略するとピン出力を出さない
    pin_powers: Option<PathBuf>,
    /// 省略すると媒質を出さない
    medium: Option<PathBuf>,
    pitch: f64,
    azimuthal: usize,
    spacing: f64,
    /// 省略すると、ライブラリにあるすべての次数（巨視的断面積では 0）
    scattering_order: Option<usize>,
    /// 省略した辺は鏡面反射
    boundary: Option<BoundaryEntry>,
    /// 省略すると、成り立つ対称操作をすべて使う（ADR-0112）
    symmetry: Option<SymmetryEntry>,
    rows: Vec<String>,
    /// 省略すると、どの列もピッチの幅
    widths: Option<Vec<SizeEntry>>,
    /// 省略すると、どの行もピッチの高さ
    heights: Option<Vec<SizeEntry>>,
    /// 省略すると余白なし
    gap: Option<GapEntry>,
    /// 巨視的断面積のときは書かない
    material: Option<BTreeMap<String, MaterialEntry>>,
    cell: BTreeMap<String, CellEntry>,
    /// 省略すると分岐なし
    branches: Option<BranchesEntry>,
}

/// 解くときに使う格子の対称性
#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
enum SymmetryEntry {
    None,
    Quarter,
    Octant,
}

/// 同じ大きさの列（行）の並び
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SizeEntry {
    size: f64,
    /// 省略すると 1
    count: Option<usize>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GapEntry {
    width: f64,
    material: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FewGroupEntry {
    /// `nodes` と排他
    output: Option<PathBuf>,
    bounds: Vec<f64>,
    nodes: Option<Vec<NodeEntry>>,
    /// 省略すると false
    pins: Option<bool>,
    /// 省略すると 1
    split: Option<usize>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct NodeEntry {
    output: PathBuf,
    /// 1 から数えた最初と最後（TOML の配列は長さを確かめないので、読んでから確かめる）
    columns: Vec<usize>,
    rows: Vec<usize>,
    /// 省略すると false
    reflector: Option<bool>,
}

impl FewGroupEntry {
    fn into_output(self) -> Result<FewGroupOutput, String> {
        let pins = self.pins.unwrap_or(false);
        let split = self.split.unwrap_or(1);
        if split == 0 {
            return Err("`split` must be positive".into());
        }
        if split > 1 && self.nodes.is_some() {
            return Err(
                "`split` (sub-nodes) needs `output` (the whole lattice), not `nodes`".into(),
            );
        }
        let groups = self.bounds.len().saturating_sub(1);
        let files = match (self.output, self.nodes) {
            (Some(output), None) => FewGroupFiles::Lattice(output),
            (None, Some(nodes)) => {
                if nodes.is_empty() {
                    return Err("`few-group` needs some nodes".into());
                }
                let range = |values: &[usize], what: &str| match *values {
                    [first, last] if first >= 1 && first <= last => Ok(first - 1..last),
                    _ => Err(format!(
                        "node {what} must be the first and last, counted from 1, got {values:?}"
                    )),
                };
                let mut files: Vec<(PathBuf, Node)> = Vec::with_capacity(nodes.len());
                for entry in nodes {
                    if files.iter().any(|(output, _)| *output == entry.output) {
                        return Err(format!("two nodes write to {}", entry.output.display()));
                    }
                    let node = Node {
                        columns: range(&entry.columns, "columns")?,
                        rows: range(&entry.rows, "rows")?,
                        reflector: entry.reflector.unwrap_or(false),
                    };
                    files.push((entry.output, node));
                }
                // 形状関数と角の不連続因子は、2 群で、ノードが格子の上の長方形の並びのときに作る（ADR-0131）
                if pins {
                    if groups != 2 {
                        return Err(format!(
                            "`pins` with `nodes` needs 2 few groups, got {groups}"
                        ));
                    }
                    let nodes: Vec<Node> = files.iter().map(|(_, n)| n.clone()).collect();
                    node_grid(&nodes).map_err(|e| format!("`pins` with `nodes`: {e}"))?;
                }
                FewGroupFiles::Nodes(files)
            }
            _ => {
                return Err(
                    "`few-group` needs exactly one of `output` (the whole lattice) and `nodes`"
                        .into(),
                );
            }
        };
        Ok(FewGroupOutput {
            bounds: self.bounds,
            files,
            pins,
            split,
        })
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BoundaryEntry {
    /// "reflective"、"vacuum"、または { albedo = α }（`condition` で読む）
    left: Option<toml::Value>,
    right: Option<toml::Value>,
    top: Option<toml::Value>,
    bottom: Option<toml::Value>,
}

/// 辺の境界条件を読む（書ける条件を、誤りの文で示す）。
fn condition(value: toml::Value) -> Result<Boundary, String> {
    let expected = "expected \"reflective\", \"vacuum\" or { albedo = number }";
    match value {
        toml::Value::String(name) => match name.as_str() {
            "reflective" => Ok(Boundary::Reflective),
            "vacuum" => Ok(Boundary::Vacuum),
            _ => Err(format!("unknown side condition `{name}`; {expected}")),
        },
        toml::Value::Table(table) => match table.get("albedo") {
            Some(albedo) if table.len() == 1 => match albedo {
                toml::Value::Float(a) => Ok(Boundary::Albedo(*a)),
                toml::Value::Integer(a) => Ok(Boundary::Albedo(*a as f64)),
                _ => Err(format!("the albedo must be a number; {expected}")),
            },
            _ => Err(format!("unknown side condition {table}; {expected}")),
        },
        other => Err(format!("unknown side condition {other}; {expected}")),
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MaterialEntry {
    /// `nuclides` と `water` で書く物質だけ（媒質と混ぜた物質は持たない）
    temperature: Option<f64>,
    /// `nuclides`、`water`、`medium`、`mix` は排他
    nuclides: Option<BTreeMap<String, f64>>,
    water: Option<WaterEntry>,
    medium: Option<PathBuf>,
    mix: Option<Vec<MixPartEntry>>,
}

/// 混ぜた物質の成分
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MixPartEntry {
    material: String,
    /// 体積の割合
    volume: f64,
}

/// 書き方ごとに分けた物質
enum MaterialKind {
    Nuclides(MaterialInput),
    Medium(MediumInput),
    Mixture(Mixture),
}

/// 混ぜた物質の割合の和と 1 の差の許容値（入力の 10 進の割合の和の丸め）
const VOLUME_TOLERANCE: f64 = 1.0e-12;

/// ホウ酸水（ADR-0106）
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WaterEntry {
    /// 密度 (g/cm³)
    density: f64,
    /// ホウ素の質量 ppm
    boron: f64,
    /// ホウ素の中の B-10 の原子分率
    b10: f64,
}

impl WaterEntry {
    fn composition(&self, name: &str) -> Result<BTreeMap<String, f64>, String> {
        borated_water(self.density, self.boron, self.b10)
            .map_err(|message| format!("material `{name}`: {message}"))
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CellEntry {
    zones: Vec<ZoneEntry>,
    /// 省略すると 1（均質なセルだけ）
    grid: Option<usize>,
    /// 省略すると 1
    sectors: Option<usize>,
    frame: Option<FrameEntry>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FrameEntry {
    thickness: f64,
    material: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ZoneEntry {
    material: String,
    /// 最後の帯だけ省略する（セルの外周まで）
    radius: Option<f64>,
    rings: usize,
}

fn positive(v: f64) -> bool {
    v > 0.0 && v.is_finite()
}

impl LatticeInput {
    pub fn parse(text: &str) -> Result<Self, InputError> {
        let file: InputFile = toml::from_str(text)?;
        let invalid = |message: String| InputError::Invalid(message);
        let cross_sections = match (file.library, file.cross_sections) {
            (Some(path), None) => CrossSections::Library(path),
            (None, Some(path)) => CrossSections::Macroscopic(path),
            _ => {
                return Err(invalid(
                    "give exactly one of `library` and `cross-sections`".into(),
                ));
            }
        };
        let macroscopic = matches!(cross_sections, CrossSections::Macroscopic(_));
        if macroscopic && file.material.is_some() {
            return Err(invalid("materials come from the cross-section file".into()));
        }
        if macroscopic && file.self_shielding.is_some() {
            return Err(invalid(
                "`self-shielding` needs a library (nuclide data)".into(),
            ));
        }
        if macroscopic && file.few_group.is_some() {
            return Err(invalid(
                "`few-group` needs a library (group energies)".into(),
            ));
        }
        if macroscopic && file.medium.is_some() {
            return Err(invalid(
                "`medium` needs a library (group energies and scattering moments)".into(),
            ));
        }
        let scattering_order = file.scattering_order;
        if macroscopic && scattering_order.is_some_and(|order| order > 0) {
            return Err(invalid(
                "`scattering-order` needs a library (macroscopic cross sections are P0 only)"
                    .into(),
            ));
        }
        for (key, value) in [("pitch", file.pitch), ("spacing", file.spacing)] {
            if !positive(value) {
                return Err(invalid(format!("`{key}` must be positive, got {value}")));
            }
        }
        if file.azimuthal == 0 {
            return Err(invalid("`azimuthal` must be positive".into()));
        }
        let boundaries = match file.boundary {
            None => Boundaries::REFLECTIVE,
            Some(b) => {
                let side = |c: Option<toml::Value>| {
                    c.map_or(Ok(Boundary::Reflective), condition)
                        .map_err(invalid)
                };
                Boundaries {
                    left: side(b.left)?,
                    right: side(b.right)?,
                    top: side(b.top)?,
                    bottom: side(b.bottom)?,
                }
            }
        };
        if file.material.as_ref().is_some_and(|m| m.contains_key(VOID)) {
            return Err(invalid(format!(
                "`{VOID}` is reserved for void zones and cannot be defined"
            )));
        }
        let output = file
            .few_group
            .map(FewGroupEntry::into_output)
            .transpose()
            .map_err(invalid)?;
        // 格子全体の形状関数と部分ノードは、外周がすべて鏡面反射の格子でだけ作る（ADR-0096、ADR-0108。格子のノード
        // ごとの形状関数はどの格子でも作る。ADR-0131）。MOC を解く前に止める
        let whole = |o: &FewGroupOutput| matches!(o.files, FewGroupFiles::Lattice(_));
        if output
            .as_ref()
            .is_some_and(|o| (o.pins && whole(o)) || o.split > 1)
            && boundaries != Boundaries::REFLECTIVE
        {
            return Err(invalid(
                "`pins` and `split` need a lattice with every side reflective".into(),
            ));
        }
        // 媒質の無限媒質の解が格子の解に一致するのは、外周がすべて鏡面反射のとき（ADR-0128）
        if file.medium.is_some() && boundaries != Boundaries::REFLECTIVE {
            return Err(invalid(
                "`medium` needs a lattice with every side reflective".into(),
            ));
        }
        let entries = file.material.unwrap_or_default();
        let self_shielding = !macroscopic && file.self_shielding.unwrap_or(true);
        // 媒質は核種の組成を持たず、混ぜた物質の遮蔽は決めていない（ADR-0128）
        if self_shielding
            && entries
                .values()
                .any(|e| e.medium.is_some() || e.mix.is_some())
        {
            return Err(invalid(
                "materials with `medium` or `mix` need `self-shielding = false` (ADR-0128)".into(),
            ));
        }
        // ホウ酸水で書いた物質の量（分岐でホウ素濃度と密度を変える）
        let waters: BTreeMap<String, Water> = entries
            .iter()
            .filter_map(|(name, entry)| {
                entry.water.as_ref().map(|w| {
                    let water = Water {
                        density: w.density,
                        boron: w.boron,
                        b10: w.b10,
                    };
                    (name.clone(), water)
                })
            })
            .collect();
        let mut materials = Vec::new();
        let mut media = Vec::new();
        let mut mixtures = Vec::new();
        for (name, entry) in entries {
            match entry.into_input(name).map_err(invalid)? {
                MaterialKind::Nuclides(m) => materials.push(m),
                MaterialKind::Medium(m) => media.push(m),
                MaterialKind::Mixture(m) => mixtures.push(m),
            }
        }
        // 成分は核種かホウ酸水で書いた物質か真空（媒質と混ぜた物質は成分にしない）
        for mixture in &mixtures {
            for (part, _) in &mixture.parts {
                if part != VOID && !materials.iter().any(|m| &m.name == part) {
                    let kind = if media.iter().any(|m| &m.name == part) {
                        "a medium"
                    } else if mixtures.iter().any(|m| &m.name == part) {
                        "a mixture"
                    } else {
                        "not defined"
                    };
                    return Err(invalid(format!(
                        "material `{}`: component `{part}` is {kind}; mix materials written \
                         with nuclides or water, or `{VOID}`",
                        mixture.name
                    )));
                }
            }
        }
        if file.branches.is_some() && matches!(cross_sections, CrossSections::Macroscopic(_)) {
            return Err(invalid(
                "branches need a library (they change the compositions and temperatures)".into(),
            ));
        }
        let branches = file
            .branches
            .map(|b| b.build(&materials, &waters, &mixtures, &media))
            .transpose()
            .map_err(invalid)?;
        if branches.is_some() && file.pin_powers.is_some() {
            return Err(invalid(
                "pin-powers cannot be written with branches (form functions go into the table; \
                 ADR-0096)"
                    .into(),
            ));
        }
        let cells = file
            .cell
            .into_iter()
            .map(|(name, entry)| entry.into_input(name))
            .collect::<Result<Vec<_>, String>>()
            .map_err(invalid)?;
        let rows = file
            .rows
            .iter()
            .map(|row| row.split_whitespace().map(String::from).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        if rows.is_empty() || rows.iter().any(Vec::is_empty) {
            return Err(invalid("`rows` needs cell names in every row".into()));
        }
        if let Some(gap) = &file.gap {
            if !positive(gap.width) {
                return Err(invalid(format!(
                    "the gap width must be positive, got {}",
                    gap.width
                )));
            }
            if gap.material == VOID {
                return Err(invalid(format!(
                    "the gap must be a material, not `{VOID}`; leave out `gap` for no gap"
                )));
            }
        }
        // 列の幅と行の高さ（同じ大きさの並びを展開する）
        let sizes = |entries: Option<Vec<SizeEntry>>, what: &str| -> Result<_, InputError> {
            let Some(entries) = entries else {
                return Ok(None);
            };
            let mut sizes = Vec::new();
            for e in entries {
                let count = e.count.unwrap_or(1);
                if !positive(e.size) || count == 0 {
                    return Err(invalid(format!(
                        "{what} need positive sizes and counts, got {} x {count}",
                        e.size
                    )));
                }
                sizes.extend(std::iter::repeat_n(e.size, count));
            }
            Ok(Some(sizes))
        };
        let widths = sizes(file.widths, "widths")?;
        let heights = sizes(file.heights, "heights")?;
        if (widths.is_some() || heights.is_some()) && file.gap.is_some() {
            return Err(invalid(
                "`gap` and `widths`/`heights` exclude each other; give the gap as columns and rows"
                    .into(),
            ));
        }
        let gap = file.gap.map(|g| Gap {
            width: g.width,
            material: g.material,
        });
        Ok(Self {
            self_shielding,
            cross_sections,
            output,
            pin_powers: file.pin_powers,
            medium: file.medium,
            branches,
            pitch: file.pitch,
            azimuthal: file.azimuthal,
            spacing: file.spacing,
            scattering_order,
            boundaries,
            symmetry: match file.symmetry {
                None => SymmetryChoice::Detect,
                Some(SymmetryEntry::None) => SymmetryChoice::Force(Symmetry::NONE),
                Some(SymmetryEntry::Quarter) => SymmetryChoice::Force(Symmetry::QUARTER),
                Some(SymmetryEntry::Octant) => SymmetryChoice::Force(Symmetry::OCTANT),
            },
            materials,
            media,
            mixtures,
            cells,
            rows,
            widths,
            heights,
            gap,
        })
    }

    /// 入力で定義した物質の名前か（核種、ホウ酸水、媒質、混ぜた物質のどれか）
    fn defines(&self, name: &str) -> bool {
        self.materials.iter().any(|m| m.name == name)
            || self.media.iter().any(|m| m.name == name)
            || self.mixtures.iter().any(|m| m.name == name)
    }

    /// 形状と、領域ごとの物質の名前。
    pub fn geometry(&self) -> Result<(Lattice, Vec<String>), InputError> {
        let mut types = Vec::with_capacity(self.cells.len());
        let mut region_materials = Vec::with_capacity(self.cells.len());
        // [セルの種類] の、ピッチに合わないマスに置くときの長方形のセル（均質なセルだけ）
        let mut blocks: Vec<Option<Cell>> = Vec::with_capacity(self.cells.len());
        for cell in &self.cells {
            // 枠の厚さは、それを使って帯を分ける前に確かめる（帯の誤りとして報告しないため）
            if let Some(frame) = &cell.frame
                && !(frame.thickness > 0.0 && frame.thickness < 0.5 * self.pitch)
            {
                return Err(InputError::Cell {
                    cell: cell.name.clone(),
                    source: GeometryError::Frame(frame.thickness),
                });
            }
            // 最後の帯を等間隔に分ける内接円（枠があれば枠の内側の正方形の内接円）の半径
            let half = 0.5 * self.pitch - cell.frame.as_ref().map_or(0.0, |f| f.thickness);
            let materials = cell
                .zones
                .iter()
                .map(|z| &z.material)
                .chain(cell.frame.as_ref().map(|f| &f.material));
            let mut radii = Vec::new();
            let mut names = Vec::new();
            let mut inner = 0.0_f64;
            for material in materials {
                if let CrossSections::Library(_) = self.cross_sections
                    && material != VOID
                    && !self.defines(material)
                {
                    return Err(InputError::UndefinedMaterial {
                        cell: cell.name.clone(),
                        material: material.clone(),
                    });
                }
            }
            for zone in &cell.zones {
                let n = zone.rings as f64;
                match zone.outer {
                    // 等体積に分ける
                    Some(r) => {
                        radii.extend((1..=zone.rings).map(|k| {
                            (inner * inner + (r * r - inner * inner) * k as f64 / n).sqrt()
                        }))
                    }
                    // セルの内接円までを等間隔に分ける（最後の領域は内接円の外の角を含む）。
                    // 内接円の外から始まる帯は、分ける幅がない
                    None if zone.rings > 1 && inner >= half => {
                        return Err(InputError::LastZoneOutside {
                            cell: cell.name.clone(),
                            start: inner,
                        });
                    }
                    None => {
                        radii.extend((1..zone.rings).map(|k| inner + (half - inner) * k as f64 / n))
                    }
                }
                inner = zone.outer.unwrap_or(half);
                names.extend(std::iter::repeat_n(zone.material.clone(), zone.rings));
            }
            // 均質なセル（帯が1つで分割1、扇形と支持格子なし）は、ピッチの正方形のマスではピンセル（同心円なし）、
            // ほかのマスでは長方形のセルにする（種類はマスで決める。下の `blocks`）
            let homogeneous = cell.zones.len() == 1
                && cell.zones[0].rings == 1
                && cell.sectors == 1
                && cell.frame.is_none();
            blocks.push(homogeneous.then_some(Cell::Block {
                nx: cell.grid,
                ny: cell.grid,
            }));
            let pin = if cell.grid > 1 {
                names = vec![names[0].clone(); cell.grid * cell.grid];
                PinCell::grid(self.pitch, cell.grid).map(Cell::Pin)
            } else {
                names = names
                    .into_iter()
                    .flat_map(|n| std::iter::repeat_n(n, cell.sectors))
                    .collect();
                let pin = PinCell::with_sectors(self.pitch, radii, cell.sectors);
                match &cell.frame {
                    Some(frame) => {
                        names.push(frame.material.clone());
                        pin.and_then(|p| p.with_frame(frame.thickness))
                    }
                    None => pin,
                }
                .map(Cell::Pin)
            }
            .map_err(|source| InputError::Cell {
                cell: cell.name.clone(),
                source,
            })?;
            types.push(pin);
            region_materials.push(names);
        }
        let layout = self
            .rows
            .iter()
            .map(|row| {
                row.iter()
                    .map(|name| {
                        self.cells
                            .iter()
                            .position(|c| &c.name == name)
                            .ok_or_else(|| InputError::UndefinedCell(name.clone()))
                    })
                    .collect::<Result<Vec<usize>, InputError>>()
            })
            .collect::<Result<Vec<_>, InputError>>()?;
        let columns = layout.first().map_or(0, Vec::len);
        let pitches = |n: usize| vec![self.pitch; n];
        let (mut widths, mut heights) = (
            self.widths.clone().unwrap_or_else(|| pitches(columns)),
            self.heights
                .clone()
                .unwrap_or_else(|| pitches(layout.len())),
        );
        // ピッチの正方形でないマスの均質なセルは、長方形のセルの種類に置き換える
        let mut layout = layout;
        let mut block_type = vec![None; types.len()];
        for (r, row) in layout.iter_mut().enumerate() {
            for (c, t) in row.iter_mut().enumerate() {
                let slot = (widths.get(c), heights.get(r));
                if slot == (Some(&self.pitch), Some(&self.pitch)) {
                    continue;
                }
                if let Some(block) = &blocks[*t] {
                    let index = *block_type[*t].get_or_insert_with(|| {
                        types.push(block.clone());
                        region_materials.push(region_materials[*t].clone());
                        types.len() - 1
                    });
                    *t = index;
                }
            }
        }
        if let Some(gap) = &self.gap {
            if let CrossSections::Library(_) = self.cross_sections
                && !self.defines(&gap.material)
            {
                return Err(InputError::UndefinedMaterial {
                    cell: "gap".into(),
                    material: gap.material.clone(),
                });
            }
            // 余白は周りの薄い行と列の均質なセル
            let border = types.len();
            types.push(Cell::BLOCK);
            region_materials.push(vec![gap.material.clone()]);
            let edge = vec![border; columns + 2];
            let mut bordered = vec![edge.clone()];
            bordered.extend(layout.into_iter().map(|row| {
                std::iter::once(border)
                    .chain(row)
                    .chain(std::iter::once(border))
                    .collect()
            }));
            bordered.push(edge);
            layout = bordered;
            for sizes in [&mut widths, &mut heights] {
                sizes.insert(0, gap.width);
                sizes.push(gap.width);
            }
        }
        let materials: Vec<String> = layout
            .iter()
            .flatten()
            .flat_map(|&t| region_materials[t].iter().cloned())
            .collect();
        let lattice = Lattice::rectilinear(widths, heights, types, layout)?;
        Ok((lattice, materials))
    }
}

impl LatticeInput {
    /// 分岐の状態ごとの物質の組（軸の直積の順）。分岐がなければ、書いたままの物質の1つの組
    pub fn branch_states(&self) -> Cow<'_, [Vec<MaterialInput>]> {
        match &self.branches {
            Some(branches) => Cow::Borrowed(branches.states()),
            None => Cow::Owned(vec![self.materials.clone()]),
        }
    }
}

/// 核種の個数密度で書いた組成を確かめる。
fn check_nuclides(
    name: &str,
    nuclides: BTreeMap<String, f64>,
) -> Result<BTreeMap<String, f64>, String> {
    if nuclides.is_empty() {
        return Err(format!("material `{name}` has no nuclides"));
    }
    for (nuclide, &density) in &nuclides {
        if !is_nuclide_name(nuclide) {
            return Err(format!(
                "material `{name}`: nuclide name `{nuclide}` is not of the form U-235, Am-242m"
            ));
        }
        if !(density >= 0.0 && density.is_finite()) {
            return Err(format!(
                "material `{name}`: invalid density {density} of {nuclide}"
            ));
        }
    }
    if nuclides.values().all(|&density| density == 0.0) {
        return Err(format!(
            "material `{name}` has no nuclide with positive density; write `{VOID}` for void zones"
        ));
    }
    Ok(nuclides)
}

impl MaterialEntry {
    fn into_input(self, name: String) -> Result<MaterialKind, String> {
        let given = [
            self.nuclides.is_some(),
            self.water.is_some(),
            self.medium.is_some(),
            self.mix.is_some(),
        ];
        if given.iter().filter(|&&g| g).count() != 1 {
            return Err(format!(
                "material `{name}` needs exactly one of nuclides, water, medium and mix"
            ));
        }
        if let Some(file) = self.medium {
            if self.temperature.is_some() {
                return Err(format!(
                    "material `{name}`: a medium has no `temperature` (its states follow the branches)"
                ));
            }
            return Ok(MaterialKind::Medium(MediumInput { name, file }));
        }
        if let Some(parts) = self.mix {
            if self.temperature.is_some() {
                return Err(format!(
                    "material `{name}`: a mixture has no `temperature` (its components have theirs)"
                ));
            }
            return mixture(name, parts).map(MaterialKind::Mixture);
        }
        let temperature = match self.temperature {
            Some(t) if positive(t) => t,
            Some(t) => return Err(format!("material `{name}`: invalid temperature {t}")),
            None => return Err(format!("material `{name}` needs a `temperature`")),
        };
        let nuclides = match (self.nuclides, self.water) {
            (Some(nuclides), _) => check_nuclides(&name, nuclides)?,
            // 作り方から、名前は正しく、密度は正で、水素と酸素を必ず含む
            (None, Some(water)) => water.composition(&name)?,
            (None, None) => return Err(format!("material `{name}` needs nuclides or water")),
        };
        Ok(MaterialKind::Nuclides(MaterialInput {
            name,
            temperature,
            composition: nuclides.into_iter().collect(),
        }))
    }
}

/// 混ぜた物質の成分を確かめる: 1つ以上あり、重ならず、割合が正で、和が丸めの範囲で 1（ADR-0128）。成分の物質が
/// 定義されていることは、すべての物質を読んでから確かめる
fn mixture(name: String, parts: Vec<MixPartEntry>) -> Result<Mixture, String> {
    if parts.is_empty() {
        return Err(format!("material `{name}`: `mix` needs components"));
    }
    for (i, part) in parts.iter().enumerate() {
        if !positive(part.volume) {
            return Err(format!(
                "material `{name}`: component `{}` needs a positive volume fraction, got {}",
                part.material, part.volume
            ));
        }
        if parts[..i].iter().any(|p| p.material == part.material) {
            return Err(format!(
                "material `{name}`: component `{}` is listed twice",
                part.material
            ));
        }
    }
    let sum: f64 = parts.iter().map(|p| p.volume).sum();
    if (sum - 1.0).abs() > VOLUME_TOLERANCE {
        return Err(format!(
            "material `{name}`: the volume fractions sum to {sum}, not 1"
        ));
    }
    Ok(Mixture {
        name,
        parts: parts.into_iter().map(|p| (p.material, p.volume)).collect(),
    })
}

impl CellEntry {
    fn into_input(self, name: String) -> Result<CellInput, String> {
        let error = |message: &str| Err(format!("cell `{name}`: {message}"));
        let Some((last, inner)) = self.zones.split_last() else {
            return error("needs at least one zone");
        };
        if last.radius.is_some() {
            return error("the last zone extends to the cell boundary (no radius)");
        }
        let mut previous = 0.0;
        for zone in inner {
            match zone.radius {
                Some(r) if positive(r) && r > previous => previous = r,
                Some(r) => return error(&format!("radii must increase, got {r}")),
                None => return error("only the last zone may omit its radius"),
            }
        }
        if self.zones.iter().any(|z| z.rings == 0) {
            return error("every zone needs at least one ring");
        }
        let (grid, sectors) = match (self.grid, self.sectors) {
            (Some(_), Some(_)) => return error("`grid` and `sectors` exclude each other"),
            (Some(0), _) | (_, Some(0)) => return error("`grid` and `sectors` must be positive"),
            (Some(m), None) => {
                if self.zones.len() != 1 || self.zones[0].rings != 1 {
                    return error("`grid` needs a homogeneous cell (one zone with one ring)");
                }
                (m, 1)
            }
            (None, s) => (1, s.unwrap_or(1)),
        };
        if self.frame.is_some() && grid > 1 {
            return error("`frame` needs a cell with rings, not `grid`");
        }
        Ok(CellInput {
            name,
            zones: self
                .zones
                .into_iter()
                .map(|z| Zone {
                    material: z.material,
                    outer: z.radius,
                    rings: z.rings,
                })
                .collect(),
            grid,
            sectors,
            frame: self.frame.map(|f| Frame {
                thickness: f.thickness,
                material: f.material,
            }),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PIN_CELL: &str = r#"
# pin cell
library = "m1.msgpack"
few-group = { output = "pincell.fg.msgpack", bounds = [1e7, 0.625, 1e-5] }
pitch = 1.26
azimuthal = 16
spacing = 0.05
rows = ["pin"]

[material.fuel]   # UO2
temperature = 900
nuclides = { U-235 = 7.0e-4, U-238 = 2.2e-2, O-16 = 4.6e-2 }

[material.water]
temperature = 600
nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }

[cell.pin]
zones = [{ material = "fuel", radius = 0.4, rings = 2 }, { material = "water", rings = 2 }]
"#;

    const LATTICE: &str = r#"
cross-sections = "c5g7.toml"
pin-powers = "pins.msgpack"
pitch = 2.0
azimuthal = 8
spacing = 0.1
boundary = { right = "vacuum", bottom = "vacuum" }
rows = ["pin box", "box pin"]

[cell.pin]
zones = [
  { material = "UO2", radius = 0.5, rings = 1 },
  { material = "clad", radius = 0.6, rings = 1 },
  { material = "water", rings = 1 },
]

[cell.box]
zones = [{ material = "water", rings = 1 }]
"#;

    fn replaced(text: &str, from: &str, to: &str) -> Result<LatticeInput, InputError> {
        assert!(text.contains(from), "{from}");
        LatticeInput::parse(&text.replace(from, to))
    }

    #[test]
    fn symmetry_is_detected_unless_given() {
        let given = |value: &str| {
            replaced(
                PIN_CELL,
                "spacing = 0.05",
                &format!("spacing = 0.05\nsymmetry = \"{value}\""),
            )
        };
        // 省略すると成り立つ対称操作をすべて使い、名前を書くとその対称性を強制する（ADR-0112）
        assert_eq!(
            LatticeInput::parse(PIN_CELL).unwrap().symmetry,
            SymmetryChoice::Detect
        );
        assert_eq!(
            given("quarter").unwrap().symmetry,
            SymmetryChoice::Force(Symmetry::QUARTER)
        );
        assert_eq!(
            given("octant").unwrap().symmetry,
            SymmetryChoice::Force(Symmetry::OCTANT)
        );
        assert_eq!(
            given("none").unwrap().symmetry,
            SymmetryChoice::Force(Symmetry::NONE)
        );
        assert!(given("half").is_err());
        assert!(given("auto").is_err());
    }

    /// 媒質の燃料と、鋼と水を混ぜた板と水を並べた平板（ADR-0128）
    const SLAB: &str = r#"
library = "m1.msgpack"
self-shielding = false
pitch = 1.26
azimuthal = 16
spacing = 0.05
boundary = { left = "vacuum", right = "vacuum" }
rows = ["fuel plate water"]

[material.fuel]
medium = "fuel.medium"

[material.ss304]
temperature = 600
nuclides = { Fe-56 = 8.0e-2 }

[material.water]
temperature = 600
water = { density = 0.743, boron = 1300, b10 = 0.199 }

[material.plate]
mix = [{ material = "ss304", volume = 0.5 }, { material = "water", volume = 0.5 }]

[cell.fuel]
zones = [{ material = "fuel", rings = 1 }]

[cell.plate]
zones = [{ material = "plate", rings = 1 }]

[cell.water]
zones = [{ material = "water", rings = 1 }]
"#;

    #[test]
    fn materials_may_be_media_or_mixtures() {
        // 物質は、核種の組成、ホウ酸水、媒質のファイル、ほかの物質を体積の割合で混ぜたもの、のどれか（ADR-0128）
        let input = LatticeInput::parse(SLAB).unwrap();
        assert_eq!(
            input.media,
            vec![MediumInput {
                name: "fuel".into(),
                file: PathBuf::from("fuel.medium"),
            }]
        );
        assert_eq!(
            input.mixtures,
            vec![Mixture {
                name: "plate".into(),
                parts: vec![("ss304".into(), 0.5), ("water".into(), 0.5)],
            }]
        );
        let names: Vec<&str> = input.materials.iter().map(|m| m.name.as_str()).collect();
        assert_eq!(names, ["ss304", "water"]);
        let (_, regions) = input.geometry().unwrap();
        assert_eq!(regions, ["fuel", "plate", "water"]);
        // 真空も成分にできる
        let void = replaced(
            SLAB,
            "{ material = \"water\", volume = 0.5 }",
            "{ material = \"void\", volume = 0.5 }",
        )
        .unwrap();
        assert_eq!(void.mixtures[0].parts[1], ("void".to_string(), 0.5));
    }

    #[test]
    fn invalid_media_and_mixtures_are_errors() {
        let bad = |from: &str, to: &str, word: &str| {
            let err = replaced(SLAB, from, to).unwrap_err().to_string();
            assert!(err.contains(word), "{to}: {err}");
        };
        let medium = "medium = \"fuel.medium\"";
        let parts = "mix = [{ material = \"ss304\", volume = 0.5 }, { material = \"water\", volume = 0.5 }]";
        // 媒質と混ぜた物質は温度を持たず、ほかの書き方と排他
        bad(
            medium,
            &format!("temperature = 600\n{medium}"),
            "temperature",
        );
        bad(parts, &format!("temperature = 600\n{parts}"), "temperature");
        bad(
            medium,
            &format!("{medium}\nnuclides = {{ U-235 = 1e-3 }}"),
            "fuel",
        );
        bad(parts, &format!("{parts}\n{medium}"), "plate");
        // 割合の和が 1 でない、割合が正でない、成分が空か重なる
        bad("volume = 0.5 }]", "volume = 0.6 }]", "sum");
        bad("volume = 0.5 }, {", "volume = -0.5 }, {", "volume");
        bad(parts, "mix = []", "plate");
        bad(
            "\"water\", volume = 0.5 }]",
            "\"ss304\", volume = 0.5 }]",
            "twice",
        );
        // 成分は核種かホウ酸水で書いた物質か真空で、知らない物質、媒質、混ぜた物質はエラー
        bad(
            "\"water\", volume = 0.5 }]",
            "\"steel\", volume = 0.5 }]",
            "steel",
        );
        bad(
            "\"water\", volume = 0.5 }]",
            "\"fuel\", volume = 0.5 }]",
            "fuel",
        );
        bad(
            "\"water\", volume = 0.5 }]",
            "\"plate\", volume = 0.5 }]",
            "plate",
        );
        // 媒質と混ぜた物質は自己遮蔽の対象にならないので、平板の格子は自己遮蔽しない（ADR-0128）
        bad("self-shielding = false\n", "", "self-shielding");
    }

    #[test]
    fn branches_list_mixtures_with_their_components() {
        // 分岐の軸は、変える成分を持つ混ぜた物質を挙げ、媒質は燃料の温度の軸に挙げる（媒質の状態を選ぶ。ADR-0128）
        let branches = "\n[branches]\n\
            boron = { ppm = [1300], materials = [\"water\", \"plate\"] }\n\
            moderator = { temperature = [600], density = [0.743], materials = [\"water\", \"ss304\", \"plate\"] }\n\
            fuel = { temperature = [600], materials = [\"fuel\"] }\n";
        let with = |text: String| LatticeInput::parse(&format!("{SLAB}{text}"));
        assert!(with(branches.into()).is_ok());
        let bad = |from: &str, to: &str, word: &str| {
            assert!(branches.contains(from), "{from}");
            let err = with(branches.replace(from, to)).unwrap_err().to_string();
            assert!(err.contains(word), "{to}: {err}");
        };
        // 変える成分を持つ混ぜた物質が挙がっていない
        bad("[\"water\", \"plate\"] }", "[\"water\"] }", "plate");
        bad("\"ss304\", \"plate\"]", "\"ss304\"]", "plate");
        // 挙げた混ぜた物質に、その軸が変える成分がない
        bad("[\"fuel\"]", "[\"fuel\", \"plate\"]", "plate");
        // 媒質は燃料の温度の軸に挙げる（媒質の状態を選ぶ）
        let without = branches
            .replace("\"ss304\", \"plate\"]", "\"plate\"]")
            .replace("[\"fuel\"]", "[\"ss304\", \"plate\"]");
        let err = with(without).unwrap_err().to_string();
        assert!(err.contains("medium") && err.contains("fuel"), "{err}");
        // 媒質を燃料の温度の軸のほかに挙げない
        bad(
            "\"ss304\", \"plate\"]",
            "\"ss304\", \"plate\", \"fuel\"]",
            "fuel",
        );
    }

    #[test]
    fn single_lattices_may_write_a_medium() {
        // 外周がすべて鏡面反射の格子は、細かい群のまま均質化した媒質を書ける（ADR-0128）
        let with =
            |text: &str, extra: &str| replaced(text, "pitch = ", &format!("{extra}\npitch = "));
        let input = with(PIN_CELL, "medium = \"fuel.medium\"").unwrap();
        assert_eq!(input.medium, Some(PathBuf::from("fuel.medium")));
        assert_eq!(LatticeInput::parse(PIN_CELL).unwrap().medium, None);
        // 鏡面反射でない辺があると、媒質の無限媒質の解が格子の解と一致しない
        let open = with(
            PIN_CELL,
            "medium = \"fuel.medium\"\nboundary = { right = \"vacuum\" }",
        );
        let err = open.unwrap_err().to_string();
        assert!(
            err.contains("medium") && err.contains("reflective"),
            "{err}"
        );
        // 群のエネルギーと散乱の次数は、ライブラリから来る
        let macroscopic = replaced(
            LATTICE,
            "boundary = { right = \"vacuum\", bottom = \"vacuum\" }",
            "medium = \"fuel.medium\"",
        );
        let err = macroscopic.unwrap_err().to_string();
        assert!(err.contains("medium") && err.contains("library"), "{err}");
    }

    #[test]
    fn void_zones_need_no_material() {
        let input = replaced(
            PIN_CELL,
            "{ material = \"fuel\", radius = 0.4, rings = 2 }",
            "{ material = \"fuel\", radius = 0.4, rings = 2 }, { material = \"void\", radius = 0.418, rings = 1 }",
        )
        .unwrap();
        // 燃料 2、真空 1、水 2 の領域
        let (_, regions) = input.geometry().unwrap();
        assert_eq!(regions[2], VOID);
    }

    #[test]
    fn materials_need_some_nuclide() {
        let result = replaced(
            PIN_CELL,
            "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
            "nuclides = { H-1 = 0.0, O-16 = 0.0 }",
        );
        assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
    }

    #[test]
    fn void_cannot_be_defined() {
        let result = replaced(PIN_CELL, "[material.water]", "[material.void]");
        assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
    }

    #[test]
    fn self_shielding_is_on_unless_turned_off() {
        assert!(LatticeInput::parse(PIN_CELL).unwrap().self_shielding);
        let off = replaced(
            PIN_CELL,
            "library = \"m1.msgpack\"",
            "library = \"m1.msgpack\"\nself-shielding = false",
        )
        .unwrap();
        assert!(!off.self_shielding);
    }

    #[test]
    fn self_shielding_needs_a_library() {
        let result = replaced(
            LATTICE,
            "cross-sections = \"c5g7.toml\"",
            "cross-sections = \"c5g7.toml\"\nself-shielding = true",
        );
        assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
    }

    #[test]
    fn scattering_order_is_optional() {
        // 省略すると、ライブラリにある次数をすべて使う（決めるのは物質を作った後）
        assert_eq!(
            LatticeInput::parse(PIN_CELL).unwrap().scattering_order,
            None
        );
        let p1 = replaced(
            PIN_CELL,
            "pitch = 1.26",
            "pitch = 1.26\nscattering-order = 1",
        )
        .unwrap();
        assert_eq!(p1.scattering_order, Some(1));
        // 巨視的断面積は P0 だけ
        let result = replaced(
            LATTICE,
            "cross-sections = \"c5g7.toml\"",
            "cross-sections = \"c5g7.toml\"\nscattering-order = 1",
        );
        assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
    }

    #[test]
    fn parses_pin_cell() {
        let input = LatticeInput::parse(PIN_CELL).unwrap();
        assert_eq!(
            input.cross_sections,
            CrossSections::Library(PathBuf::from("m1.msgpack"))
        );
        assert_eq!(
            input.output,
            Some(FewGroupOutput {
                bounds: vec![1e7, 0.625, 1e-5],
                files: FewGroupFiles::Lattice(PathBuf::from("pincell.fg.msgpack")),
                pins: false,
                split: 1,
            })
        );
        assert_eq!(input.pin_powers, None);
        assert_eq!(
            (input.pitch, input.azimuthal, input.spacing),
            (1.26, 16, 0.05)
        );
        assert_eq!(input.boundaries, Boundaries::REFLECTIVE);
        // 物質と核種は名前の順
        assert_eq!(input.materials.len(), 2);
        assert_eq!(input.materials[0].name, "fuel");
        assert_eq!(input.materials[0].temperature, 900.0);
        assert_eq!(
            input.materials[0].composition[1],
            ("U-235".to_string(), 7.0e-4)
        );
        assert_eq!(
            input.cells[0].zones,
            vec![
                Zone {
                    material: "fuel".into(),
                    outer: Some(0.4),
                    rings: 2
                },
                Zone {
                    material: "water".into(),
                    outer: None,
                    rings: 2
                },
            ]
        );
        assert_eq!((input.cells[0].grid, input.cells[0].sectors), (1, 1));
        assert_eq!(input.rows, vec![vec!["pin".to_string()]]);
    }

    #[test]
    fn pin_cell_geometry_splits_fuel_by_volume_and_water_by_radius() {
        let input = LatticeInput::parse(PIN_CELL).unwrap();
        let (lattice, materials) = input.geometry().unwrap();
        assert_eq!(materials, vec!["fuel", "fuel", "water", "water"]);
        // 燃料は等体積の2つのリング、水はペレットとセルの内接円の中点で分ける
        let v = lattice.volumes();
        assert!((v[0] - v[1]).abs() < 1e-14);
        let r = (0.4 + 0.63) / 2.0;
        assert!((v[2] - std::f64::consts::PI * (r * r - 0.16)).abs() < 1e-14);
    }

    #[test]
    fn zones_may_reach_beyond_the_inscribed_circle() {
        let geometry = |zones: &str| {
            let text = PIN_CELL.replace(
                "zones = [{ material = \"fuel\", radius = 0.4, rings = 2 }, { material = \"water\", rings = 2 }]",
                zones,
            );
            LatticeInput::parse(&text).and_then(|i| i.geometry())
        };
        // 内接円（0.63）の外の 0.8 までの水の帯は、セルの外周で切り取った2つのリングになる
        let (lattice, materials) = geometry(
            r#"zones = [
  { material = "fuel", radius = 0.4, rings = 2 },
  { material = "water", radius = 0.8, rings = 2 },
  { material = "water", rings = 1 },
]"#,
        )
        .unwrap();
        assert_eq!(materials.len(), 5);
        let v = lattice.volumes();
        assert!((v.iter().sum::<f64>() - 1.26 * 1.26).abs() < 1e-14);
        // 最後の帯は内接円の外から始まるので、等間隔には分けられない
        assert!(matches!(
            geometry(
                r#"zones = [
  { material = "fuel", radius = 0.4, rings = 2 },
  { material = "water", radius = 0.8, rings = 2 },
  { material = "water", rings = 2 },
]"#
            ),
            Err(InputError::LastZoneOutside { .. })
        ));
    }

    #[test]
    fn parses_lattice() {
        let input = LatticeInput::parse(LATTICE).unwrap();
        assert_eq!(
            input.cross_sections,
            CrossSections::Macroscopic(PathBuf::from("c5g7.toml"))
        );
        assert_eq!(input.output, None);
        assert_eq!(input.pin_powers, Some(PathBuf::from("pins.msgpack")));
        assert_eq!(
            input.boundaries,
            Boundaries {
                right: Boundary::Vacuum,
                bottom: Boundary::Vacuum,
                ..Boundaries::REFLECTIVE
            }
        );
        let (lattice, materials) = input.geometry().unwrap();
        assert_eq!(lattice.width(), 4.0);
        // セルは配置の順（上の行から）
        assert_eq!(
            materials,
            vec![
                "UO2", "clad", "water", "water", "water", "UO2", "clad", "water"
            ]
        );
    }

    #[test]
    fn columns_and_rows_may_have_their_own_sizes() {
        // ピンの格子の右に、厚さ 0.5 のバッフルと 1.5 の水の列を足す（反射体）
        let reflector = |rows: &str, extra: &str| {
            let text = LATTICE
                .replace(
                    "rows = [\"pin box\", \"box pin\"]",
                    &format!("widths = [{{ size = 2.0, count = 2 }}, {{ size = 0.5 }}, {{ size = 1.5 }}]\n{extra}rows = {rows}"),
                )
                .replace(
                    "[cell.box]",
                    "[cell.baffle]\nzones = [{ material = \"clad\", rings = 1 }]\n\n[cell.box]",
                );
            LatticeInput::parse(&text).and_then(|input| input.geometry())
        };
        let (lattice, materials) =
            reflector("[\"pin box baffle box\", \"box pin baffle box\"]", "").unwrap();
        assert_eq!(lattice.widths(), &[2.0, 2.0, 0.5, 1.5]);
        assert_eq!(lattice.heights(), &[2.0, 2.0]);
        // ピン 3 + 水 1 + バッフル 1 + 水 1 の行が 2 つ
        assert_eq!(lattice.regions(), 12);
        assert_eq!(materials[4], "clad");
        // 均質なセルは、ピッチのマスではピンセル、ほかのマスでは長方形のセル
        assert!(lattice.cell_type(1).pin().is_some());
        assert!(lattice.cell_type(3).pin().is_none());
        assert_eq!(lattice.pin_rows_and_columns(), (vec![0, 1], vec![0, 1]));
        // ピンはピッチの幅の列にだけ置ける
        assert!(reflector("[\"pin box pin box\", \"box pin baffle box\"]", "").is_err());
        // 列の数がそろわない、余白と同時には使えない
        assert!(reflector("[\"pin box baffle\", \"box pin baffle\"]", "").is_err());
        assert!(
            reflector(
                "[\"pin box baffle box\", \"box pin baffle box\"]",
                "gap = { width = 0.1, material = \"water\" }\n"
            )
            .is_err()
        );
    }

    #[test]
    fn sides_may_have_an_albedo() {
        let input = replaced(
            LATTICE,
            "boundary = { right = \"vacuum\", bottom = \"vacuum\" }",
            "boundary = { right = \"vacuum\", top = { albedo = 0.3 } }",
        )
        .unwrap();
        assert_eq!(input.boundaries.top, Boundary::Albedo(0.3));
        assert_eq!(input.boundaries.right, Boundary::Vacuum);
        assert_eq!(input.boundaries.bottom, Boundary::Reflective);
        assert!(
            replaced(
                LATTICE,
                "boundary = { right = \"vacuum\", bottom = \"vacuum\" }",
                "boundary = { top = { albedo = 0.3, extra = 1 } }",
            )
            .is_err()
        );
        // 書き間違いには、書ける条件を示す
        let error = replaced(
            LATTICE,
            "boundary = { right = \"vacuum\", bottom = \"vacuum\" }",
            "boundary = { right = \"vacum\" }",
        )
        .unwrap_err()
        .to_string();
        assert!(
            error.contains("reflective") && error.contains("albedo"),
            "{error}"
        );
    }

    #[test]
    fn gap_surrounds_the_lattice() {
        let text = LATTICE.replace(
            "rows = [",
            "gap = { width = 0.25, material = \"water\" }\nrows = [",
        );
        let (lattice, materials) = LatticeInput::parse(&text).unwrap().geometry().unwrap();
        assert_eq!(lattice.width(), 4.5);
        assert_eq!(lattice.widths(), &[0.25, 2.0, 2.0, 0.25]);
        // 余白は周りの薄い行と列のセル（各 1 領域）で、領域はセルの順に並ぶ
        assert_eq!(materials.len(), lattice.regions());
        assert_eq!(materials.len(), 8 + 12);
        assert!(materials[..4].iter().all(|m| m == "water"));
        assert_eq!(materials[4..9], ["water", "UO2", "clad", "water", "water"]);
        assert_eq!(materials.iter().filter(|m| *m == "UO2").count(), 2);
        // 幅が正でない、物質がない、真空（巨視的断面積でも）
        for gap in [
            "{ width = 0.0, material = \"water\" }",
            "{ width = 0.1 }",
            "{ width = 0.1, material = \"void\" }",
        ] {
            let text = LATTICE.replace("rows = [", &format!("gap = {gap}\nrows = ["));
            let result = LatticeInput::parse(&text).and_then(|input| input.geometry().map(|_| ()));
            assert!(result.is_err(), "{gap}");
        }
    }

    #[test]
    fn frames_are_cell_edges() {
        let with_frame = |frame: &str| {
            replaced(
                LATTICE,
                "  { material = \"water\", rings = 1 },\n]",
                &format!("  {{ material = \"water\", rings = 2 }},\n]\nframe = {frame}"),
            )
        };
        let input = with_frame("{ thickness = 0.05, material = \"clad\" }").unwrap();
        let (lattice, regions) = input.geometry().unwrap();
        // ピンセル 4 + 枠 1 が 2 つ、均質なセル 1 が 2 つ
        assert_eq!(lattice.regions(), 2 * 5 + 2);
        assert_eq!(regions[4], "clad");
        // 厚すぎる枠は、帯の配置ではなく枠の誤りとして報告する
        for thickness in ["1.0", "-0.1", "nan"] {
            let frame = format!("{{ thickness = {thickness}, material = \"clad\" }}");
            let error = with_frame(&frame).and_then(|input| input.geometry().map(|_| ()));
            assert!(
                matches!(
                    error,
                    Err(InputError::Cell {
                        source: GeometryError::Frame(_),
                        ..
                    })
                ),
                "{thickness}: {error:?}"
            );
        }
        assert!(with_frame("{ thickness = 0.05 }").is_err());
    }

    #[test]
    fn grid_and_sectors_split_cells() {
        let text = LATTICE.replace(
            "[cell.box]\nzones = [{ material = \"water\", rings = 1 }]",
            "[cell.box]\nzones = [{ material = \"water\", rings = 1 }]\ngrid = 3",
        );
        let input = LatticeInput::parse(&text).unwrap();
        let (lattice, _) = input.geometry().unwrap();
        assert_eq!(lattice.regions(), 3 + 9 + 9 + 3);
        let text = LATTICE.replace(
            "  { material = \"water\", rings = 1 },\n]",
            "  { material = \"water\", rings = 1 },\n]\nsectors = 8",
        );
        let (lattice, _) = LatticeInput::parse(&text).unwrap().geometry().unwrap();
        assert_eq!(lattice.regions(), 2 * 24 + 2);
    }

    #[test]
    fn invalid_inputs_are_errors() {
        // 断面積の出どころは library と cross-sections のどちらか一方
        assert!(LatticeInput::parse(&format!("cross-sections = \"a.toml\"\n{PIN_CELL}")).is_err());
        assert!(replaced(PIN_CELL, "library = \"m1.msgpack\"\n", "").is_err());
        // 巨視的断面積では material を書かない、少数群定数を出さない
        assert!(
            LatticeInput::parse(&format!(
                "{LATTICE}\n[material.UO2]\ntemperature = 600\nnuclides = {{ U-235 = 1e-3 }}"
            ))
            .is_err()
        );
        assert!(
            replaced(
                LATTICE,
                "pitch = 2.0",
                "pitch = 2.0\nfew-group = { output = \"a\", bounds = [1e7, 1e-5] }"
            )
            .is_err()
        );
        // 知らないキー、知らない辺と条件
        assert!(replaced(PIN_CELL, "pitch = 1.26", "pitch = 1.26\nalbedo = 1").is_err());
        assert!(replaced(LATTICE, "right = \"vacuum\"", "east = \"vacuum\"").is_err());
        assert!(replaced(LATTICE, "right = \"vacuum\"", "right = \"white\"").is_err());
        // 正でない値、温度、個数密度、核種名
        assert!(replaced(PIN_CELL, "pitch = 1.26", "pitch = 0").is_err());
        assert!(replaced(PIN_CELL, "spacing = 0.05", "spacing = -0.05").is_err());
        assert!(replaced(PIN_CELL, "azimuthal = 16", "azimuthal = 0").is_err());
        assert!(replaced(PIN_CELL, "temperature = 900", "temperature = -900").is_err());
        assert!(replaced(PIN_CELL, "U-238 = 2.2e-2", "U-238 = -2.2e-2").is_err());
        assert!(replaced(PIN_CELL, "U-238 = 2.2e-2", "U238 = 2.2e-2").is_err());
        assert!(
            replaced(
                PIN_CELL,
                "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
                "nuclides = {}"
            )
            .is_err()
        );
        // 帯: 外半径が昇順でない、最後の帯に外半径がある、途中の帯に外半径がない、分割数が 0
        assert!(replaced(LATTICE, "radius = 0.6", "radius = 0.4").is_err());
        assert!(
            replaced(
                PIN_CELL,
                "{ material = \"water\", rings = 2 }",
                "{ material = \"water\", radius = 0.6, rings = 2 }"
            )
            .is_err()
        );
        assert!(
            replaced(
                PIN_CELL,
                "{ material = \"fuel\", radius = 0.4, rings = 2 }",
                "{ material = \"fuel\", rings = 2 }"
            )
            .is_err()
        );
        assert!(replaced(PIN_CELL, "rings = 2 }]", "rings = 0 }]").is_err());
        // grid は均質なセル（帯が1つで分割数 1）だけ、grid と sectors は同時に使わない
        assert!(replaced(PIN_CELL, "rings = 2 }]", "rings = 2 }]\ngrid = 2").is_err());
        assert!(
            replaced(
                LATTICE,
                "zones = [{ material = \"water\", rings = 1 }]",
                "zones = [{ material = \"water\", rings = 1 }]\ngrid = 2\nsectors = 4"
            )
            .is_err()
        );
        assert!(
            replaced(
                LATTICE,
                "zones = [{ material = \"water\", rings = 1 }]",
                "zones = [{ material = \"water\", rings = 1 }]\ngrid = 0"
            )
            .is_err()
        );
    }

    #[test]
    fn geometry_checks_names_and_layout() {
        let geometry = |text: String| LatticeInput::parse(&text).and_then(|i| i.geometry());
        // 配置に知らないセル、行の長さが行の数と違う
        assert!(geometry(LATTICE.replace("\"box pin\"", "\"box rod\"")).is_err());
        assert!(geometry(LATTICE.replace("\"box pin\"", "\"box pin pin\"")).is_err());
        // ライブラリのときは、セルの物質が material で定義されていること
        assert!(
            geometry(PIN_CELL.replace(
                "{ material = \"water\", rings = 2 }",
                "{ material = \"steam\", rings = 2 }"
            ))
            .is_err()
        );
        // 帯の外半径がセルの外接円（ピッチ 2 で √2）を超える、扇形の数が 1, 4, 8 でない
        assert!(geometry(LATTICE.replace("radius = 0.6", "radius = 1.5")).is_err());
        assert!(
            geometry(LATTICE.replace(
                "  { material = \"water\", rings = 1 },\n]",
                "  { material = \"water\", rings = 1 },\n]\nsectors = 3"
            ))
            .is_err()
        );
    }

    #[test]
    fn cell_errors_keep_the_cause_in_the_chain_only() {
        // 原因の文は連鎖（source）にだけ置く。自分の文にも埋め込むと、連鎖を表示するときに2回出る
        let text = LATTICE.replace("radius = 0.6", "radius = 1.5");
        let err = LatticeInput::parse(&text).unwrap().geometry().unwrap_err();
        let cause = std::error::Error::source(&err).expect("the geometry error is the cause");
        assert!(!err.to_string().contains(&cause.to_string()), "{err}");
    }

    #[test]
    fn toml_errors_are_not_repeated_in_the_chain() {
        // 包んだ TOML のエラーは、自分の文と原因の連鎖のどちらか一方にだけ出す
        let err = LatticeInput::parse("pitch = ").unwrap_err();
        let chain = format!("{:#}", anyhow::Error::from(err));
        assert_eq!(chain.matches("TOML parse error").count(), 1, "{chain}");
    }

    #[test]
    fn parses_few_group_nodes() {
        // 格子をノードに分け、ノードごとに少数群定数を書く（ADR-0095）。列と行は格子のセル（余白の薄い行と
        // 列を含む）の、1 から数えた最初と最後。LATTICE は巨視的断面積なので、ライブラリに置き換える
        let text = LATTICE.replace("cross-sections = \"c5g7.toml\"\n", "");
        let input = replaced(
            &text,
            "pin-powers = \"pins.msgpack\"",
            "library = \"m1.msgpack\"\nfew-group = { bounds = [1e7, 1e-5], nodes = [\n  \
             { output = \"fuel.msgpack\", columns = [1, 1], rows = [1, 2] },\n  \
             { output = \"reflector.msgpack\", columns = [2, 2], rows = [1, 2], reflector = true },\n] }",
        );
        let input = input.unwrap_or_else(|e| panic!("{e}"));
        let expected = FewGroupOutput {
            bounds: vec![1e7, 1e-5],
            files: FewGroupFiles::Nodes(vec![
                (
                    PathBuf::from("fuel.msgpack"),
                    Node {
                        columns: 0..1,
                        rows: 0..2,
                        reflector: false,
                    },
                ),
                (
                    PathBuf::from("reflector.msgpack"),
                    Node {
                        columns: 1..2,
                        rows: 0..2,
                        reflector: true,
                    },
                ),
            ]),
            pins: false,
            split: 1,
        };
        assert_eq!(input.output, Some(expected));
    }

    #[test]
    fn few_group_may_ask_for_pin_form_functions() {
        // ピン出力の再構成の形状関数と角の不連続因子を付ける（ADR-0096）
        let few_group = |entry: &str| {
            replaced(
                PIN_CELL,
                "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
                &format!("few-group = {{ bounds = [1e7, 1e-5]{entry} }}"),
            )
        };
        let input = few_group(", output = \"a\", pins = true").unwrap();
        assert!(input.output.is_some_and(|o| o.pins));
        let input = few_group(", output = \"a\", pins = false").unwrap();
        assert!(input.output.is_some_and(|o| !o.pins));
        // 格子のノードごとの定数にも形状関数と角の不連続因子を付けられる（ADR-0131）。2 群で、ノードが格子の上の
        // 長方形の並びであることを、MOC を解く前に確かめる
        let with_nodes = |bounds: &str, nodes: &str| {
            replaced(
                PIN_CELL,
                "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
                &format!("few-group = {{ bounds = {bounds}, nodes = [{nodes}], pins = true }}"),
            )
        };
        let one = "{ output = \"a\", columns = [1, 1], rows = [1, 1] }";
        let input = with_nodes("[1e7, 0.625, 1e-5]", one).unwrap();
        assert!(input.output.is_some_and(|o| o.pins));
        let err = with_nodes("[1e7, 1e-5]", one).unwrap_err().to_string();
        assert!(err.contains("2 few groups"), "{err}");
        // 左のノードが 2 行にまたがり、右は 1 行ずつ: 列と行の切れ目の長方形の並びにならない
        let uneven = "{ output = \"a\", columns = [1, 1], rows = [1, 2] }, \
                      { output = \"b\", columns = [2, 2], rows = [1, 1] }, \
                      { output = \"c\", columns = [2, 2], rows = [2, 2] }";
        let err = with_nodes("[1e7, 0.625, 1e-5]", uneven)
            .unwrap_err()
            .to_string();
        assert!(err.contains("grid"), "{err}");
    }

    #[test]
    fn few_group_may_split_the_lattice() {
        // 格子全体を n × n の部分ノードに分けた定数（ADR-0108）。省略すると 1。nodes とは使えない
        let few_group = |entry: &str| {
            replaced(
                PIN_CELL,
                "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
                &format!("few-group = {{ bounds = [1e7, 1e-5]{entry} }}"),
            )
        };
        let input = few_group(", output = \"a\", split = 2").unwrap();
        assert_eq!(input.output.map(|o| o.split), Some(2));
        let input = few_group(", output = \"a\"").unwrap();
        assert_eq!(input.output.map(|o| o.split), Some(1));
        assert!(few_group(", output = \"a\", split = 0").is_err());
        let nodes = ", nodes = [{ output = \"a\", columns = [1, 1], rows = [1, 1] }], split = 2";
        let err = few_group(nodes).unwrap_err().to_string();
        assert!(err.contains("split"), "{err}");
        // 外周がすべて鏡面反射でない格子は、MOC を解く前に入力の段階で止める（pins も同じ）
        for entry in [
            ", output = \"a\", split = 2",
            ", output = \"a\", pins = true",
        ] {
            let open = replaced(
                PIN_CELL,
                "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
                &format!(
                    "few-group = {{ bounds = [1e7, 1e-5]{entry} }}\nboundary = {{ right = \"vacuum\" }}"
                ),
            );
            let err = open.unwrap_err().to_string();
            assert!(err.contains("reflective"), "{entry}: {err}");
        }
    }

    #[test]
    fn invalid_few_group_nodes_are_errors() {
        let few_group = |entry: &str| {
            replaced(
                PIN_CELL,
                "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
                &format!("few-group = {{ bounds = [1e7, 1e-5]{entry} }}"),
            )
        };
        let node = |output: &str, columns: &str, rows: &str| {
            format!("{{ output = \"{output}\", columns = {columns}, rows = {rows} }}")
        };
        let good = node("a", "[1, 1]", "[1, 1]");
        assert!(few_group(&format!(", nodes = [{good}]")).is_ok());
        for (entry, what) in [
            (String::new(), "neither output nor nodes"),
            (
                format!(", output = \"x\", nodes = [{good}]"),
                "both output and nodes",
            ),
            (", nodes = []".to_string(), "no nodes"),
            (
                format!(", nodes = [{}]", node("a", "[0, 1]", "[1, 1]")),
                "column 0",
            ),
            (
                format!(", nodes = [{}]", node("a", "[2, 1]", "[1, 1]")),
                "reversed columns",
            ),
            (
                format!(", nodes = [{}]", node("a", "[1, 1]", "[3, 2]")),
                "reversed rows",
            ),
            (
                format!(", nodes = [{good}, {}]", node("a", "[2, 2]", "[1, 1]")),
                "same output",
            ),
            (
                format!(", nodes = [{}]", node("a", "[1, 2, 3]", "[1, 1]")),
                "three numbers",
            ),
        ] {
            assert!(
                matches!(
                    few_group(&entry),
                    Err(InputError::Invalid(_) | InputError::Toml(_))
                ),
                "{what}"
            );
        }
    }

    #[test]
    fn borated_water_matches_vera_number_densities() {
        // VERA の Problem 5 の減速材（565 K、0.743 g/cm³、B-10 19.9 at%）の個数密度の表（ADR-0106）
        let table = [
            (1285.0, [4.96231e-2, 1.05835e-5, 4.25999e-5, 2.48116e-2]),
            (1291.0, [4.96228e-2, 1.06329e-5, 4.27988e-5, 2.48114e-2]),
            (1170.0, [4.96288e-2, 9.63633e-6, 3.87874e-5, 2.48144e-2]),
            (1230.0, [4.96258e-2, 1.01305e-5, 4.07765e-5, 2.48129e-2]),
        ];
        for (ppm, expected) in table {
            let input = replaced(
                PIN_CELL,
                "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
                &format!("water = {{ density = 0.743, boron = {ppm}, b10 = 0.199 }}"),
            )
            .unwrap();
            let water = input.materials.iter().find(|m| m.name == "water").unwrap();
            let names: Vec<&str> = water.composition.iter().map(|(n, _)| n.as_str()).collect();
            assert_eq!(names, ["B-10", "B-11", "H-1", "O-16"]);
            let [h, b10, b11, o] = expected;
            for ((_, got), want) in water.composition.iter().zip([b10, b11, h, o]) {
                // 表は有効数字 6 桁
                assert!(
                    ((got - want) / want).abs() < 1e-5,
                    "{ppm} ppm: {got} vs {want}"
                );
            }
        }
    }

    #[test]
    fn water_without_boron_has_no_boron_nuclides() {
        let input = replaced(
            PIN_CELL,
            "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
            "water = { density = 0.743, boron = 0, b10 = 0.199 }",
        )
        .unwrap();
        let water = input.materials.iter().find(|m| m.name == "water").unwrap();
        let names: Vec<&str> = water.composition.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(names, ["H-1", "O-16"]);
        assert!(water.composition.iter().all(|&(_, n)| n > 0.0));
    }

    #[test]
    fn water_keeps_only_present_boron_isotopes() {
        for (b10, names) in [
            (1.0, ["B-10", "H-1", "O-16"]),
            (0.0, ["B-11", "H-1", "O-16"]),
        ] {
            let input = replaced(
                PIN_CELL,
                "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
                &format!("water = {{ density = 0.743, boron = 1291, b10 = {b10:?} }}"),
            )
            .unwrap();
            let water = input.materials.iter().find(|m| m.name == "water").unwrap();
            let got: Vec<&str> = water.composition.iter().map(|(n, _)| n.as_str()).collect();
            assert_eq!(got, names);
        }
    }

    #[test]
    fn invalid_water_is_an_error() {
        let water = "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }";
        let bad = |to: &str| replaced(PIN_CELL, water, to).is_err();
        assert!(bad("water = { density = 0.0, boron = 1291, b10 = 0.199 }"));
        assert!(bad("water = { density = 0.743, boron = -1, b10 = 0.199 }"));
        assert!(bad("water = { density = 0.743, boron = 1e6, b10 = 0.199 }"));
        assert!(bad("water = { density = 0.743, boron = 1291, b10 = 1.5 }"));
        // 3つとも必須、nuclides と同時には書けない
        assert!(bad("water = { density = 0.743, boron = 1291 }"));
        assert!(bad(&format!(
            "{water}\nwater = {{ density = 0.743, boron = 1291, b10 = 0.199 }}"
        )));
        assert!(bad(""));
    }
}
