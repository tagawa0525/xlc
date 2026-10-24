//! 格子計算の入力ファイル（TOML）。
//!
//! ```toml
//! library = "m1.msgpack"          # 多群ライブラリ（cross-sections と排他）
//! self-shielding = false          # 燃料の自己遮蔽（library のとき。省略すると true）
//! track-table = 16               # 自己遮蔽の減速計算の飛程の表の、ln Σ の 10 倍あたりの格子点の数（省略すると 16。
//!                                 # false で格子点ごとに求める。ADR-0156）
//! ufg-tolerance = 1e-4            # 自己遮蔽の減速計算の格子を、リングの物質の巨視的な断面積をこの相対差で再現できる点に
//!                                 # 間引く（省略すると間引かない。ADR-0185）
//! ufg-scattering = "free-gas"     # 自己遮蔽の減速計算の散乱核の標的（"static" か "free-gas"。省略すると static。ADR-0203）
//! ufg-free-gas = { upper = 400, cells = 4, sweeps = 2 }   # free-gas の条件（どれも省略できる。ADR-0203）: 自由ガスの核を
//!                                 # 使う入射エネルギーの上限（kT の倍数。省略すると上限なし）、核の差を積分する区間の
//!                                 # 細かさ、上方散乱のスイープの回数（省略すると収束するまで）
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
//! # axial = true を足すと、格子の x を炉心の z として定数を書く（軸方向の反射体の平板。上下の辺は鏡面反射で、
//! # ノードは格子のすべての行にわたる。pins と split とは使えない。ADR-0128）
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
//!
//! [depletion]                     # 燃焼（ADR-0173。分岐、pin-powers、媒質、ファイルの few-group とは使えない）
//! chain = "chain.msgpack"         # 燃焼チェーン
//! power = 38.41                   # 初めの重金属の出力密度 (W/g)
//! burnup = [0.1, 0.5, 1.0]        # 燃焼度の点 (MWd/kg)
//! materials = ["fuel"]            # 燃焼させる物質（核種で書いたもの）
//! output = "history.msgpack"      # 燃焼の履歴（燃焼点ごとの k∞、数密度、ピン出力。ADR-0181）
//! few-group = { bounds = [1e7, 0.625, 1e-5] }   # 燃焼点ごとの少数群定数（履歴に入る。省略すると作らない）
//! shielding = "every-transport"  # 自己遮蔽を解く輸送（"every-transport" か "step-start"。省略すると every-transport。ADR-0197）
//!
//! [history]                       # 燃焼の履歴の燃焼点からの分岐（ADR-0187。[branches] とファイルの few-group と組む）
//! file = "history.msgpack"        # 燃焼の履歴（`[depletion]` の output）
//! burnup = [0, 10, 20]            # 分岐を解く燃焼点 (MWd/kg)。履歴の燃焼度から選ぶ
//! off-nominal-burnup = [0, 20]    # 公称でない状態を解く燃焼点（省略すると burnup のすべて。ADR-0200）
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
use crate::depletion::{DepletionFewGroup, DepletionInput, ShieldingFrequency};
use crate::fewgroup::{Node, node_grid};
use crate::geometry::{Cell, GeometryError, Lattice, PinCell};
use crate::history::HistoryInput;
use crate::shielding::ShieldingOptions;
use crate::symmetry::Symmetry;
use crate::tracking::{Boundaries, Boundary};
use crate::ufg::{Collision, FreeGas, Scattering};
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
    /// 自己遮蔽の UFG の減速計算の飛程の求め方（ADR-0156）
    pub collision: Collision,
    /// 自己遮蔽の UFG の減速計算の格子の間引きの許容差（ADR-0185。`None` は間引かない）
    pub thinning: Option<f64>,
    /// 自己遮蔽の UFG の減速計算の散乱核の標的の扱い（ADR-0203）
    pub scattering: Scattering,
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
    /// 燃焼（ADR-0173）。`None` なら燃焼させない
    pub depletion: Option<DepletionInput>,
    /// 燃焼の履歴の燃焼点からの分岐（ADR-0187）。`None` なら物質は入力のまま
    pub history: Option<HistoryInput>,
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
    /// 格子の x を炉心の z とした軸方向の定数を書く（ADR-0128）
    pub axial: bool,
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

/// `ufg-free-gas` の値
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct FreeGasEntry {
    upper: Option<f64>,
    cells: Option<i64>,
    sweeps: Option<i64>,
}

/// `track-table` の値: 10 倍あたりの格子点の数か、false（格子点ごとに求める）
#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum TrackTableEntry {
    Density(i64),
    Exact(bool),
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct InputFile {
    library: Option<PathBuf>,
    /// 省略すると、ライブラリのときは true
    self_shielding: Option<bool>,
    /// 省略すると ADR-0156 の既定の細かさ
    track_table: Option<TrackTableEntry>,
    /// 省略すると間引かない（ADR-0185）
    ufg_tolerance: Option<f64>,
    /// 省略すると静止した標的（ADR-0203）
    ufg_scattering: Option<String>,
    /// `ufg-scattering = "free-gas"` の条件。省略すると ADR-0203 の既定
    ufg_free_gas: Option<FreeGasEntry>,
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
    /// 省略すると燃焼させない
    depletion: Option<DepletionEntry>,
    /// 省略すると燃焼の履歴を読まない
    history: Option<HistoryEntry>,
}

/// 燃焼の履歴の燃焼点からの分岐（ADR-0187、ADR-0200）
#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct HistoryEntry {
    file: PathBuf,
    burnup: Vec<f64>,
    off_nominal_burnup: Option<Vec<f64>>,
}

/// 燃焼の条件（ADR-0173）
#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct DepletionEntry {
    chain: PathBuf,
    power: f64,
    burnup: Vec<f64>,
    materials: Vec<String>,
    output: PathBuf,
    /// 省略すると少数群定数を作らない
    few_group: Option<DepletionFewGroupEntry>,
    /// 省略すると輸送のたびに解く（ADR-0197）
    shielding: Option<ShieldingEntry>,
}

/// 燃焼の自己遮蔽を解く輸送（ADR-0197）
#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
enum ShieldingEntry {
    EveryTransport,
    StepStart,
}

/// 燃焼点ごとの格子全体の少数群定数（燃焼の履歴に書く）
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct DepletionFewGroupEntry {
    bounds: Vec<f64>,
    /// 省略すると false
    pins: Option<bool>,
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
    /// 省略すると false
    axial: Option<bool>,
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
        let axial = self.axial.unwrap_or(false);
        // 形状関数と部分ノードは径方向の量（ADR-0128）
        if axial && (pins || split > 1) {
            return Err("`axial` constants take no `pins` or `split`".into());
        }
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
            axial,
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

/// `ufg-scattering` と `ufg-free-gas` から散乱核の標的の扱い（ADR-0203）
fn parse_scattering(kind: Option<&str>, entry: Option<FreeGasEntry>) -> Result<Scattering, String> {
    match kind {
        None | Some("static") => match entry {
            None => Ok(Scattering::Static),
            Some(_) => Err("`ufg-free-gas` needs `ufg-scattering = \"free-gas\"`".into()),
        },
        Some("free-gas") => {
            return Ok(Scattering::Static);
            #[allow(unreachable_code)]
            let defaults = FreeGas::default();
            let Some(entry) = entry else {
                return Ok(Scattering::FreeGas(defaults));
            };
            let upper = match entry.upper {
                Some(v) if !positive(v) => {
                    return Err(format!(
                        "`ufg-free-gas.upper` must be a positive multiple of kT, got {v}"
                    ));
                }
                v => v.or(defaults.upper),
            };
            let count = |name: &str, v: Option<i64>| -> Result<Option<usize>, String> {
                match v {
                    Some(n) if n < 1 => Err(format!(
                        "`ufg-free-gas.{name}` must be a positive integer, got {n}"
                    )),
                    Some(n) => Ok(Some(n as usize)),
                    None => Ok(None),
                }
            };
            Ok(Scattering::FreeGas(FreeGas {
                upper,
                cells: count("cells", entry.cells)?.unwrap_or(defaults.cells),
                sweeps: count("sweeps", entry.sweeps)?.or(defaults.sweeps),
            }))
        }
        Some(other) => Err(format!(
            "`ufg-scattering` must be \"static\" or \"free-gas\", not {other:?}"
        )),
    }
}

fn positive(v: f64) -> bool {
    v > 0.0 && v.is_finite()
}

impl LatticeInput {
    /// 自己遮蔽の条件。入力で選べるのは飛程の求め方（`track-table`）、格子の間引き（`ufg-tolerance`）、散乱核の標的の扱い
    /// （`ufg-scattering`、`ufg-free-gas`）だけで、ほかは既定の値
    pub fn shielding_options(&self) -> ShieldingOptions {
        ShieldingOptions {
            collision: self.collision,
            thinning: self.thinning,
            scattering: self.scattering,
            ..ShieldingOptions::default()
        }
    }

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
        // 平板の格子の y は炉心の径方向で、一様（ADR-0128）
        if output.as_ref().is_some_and(|o| o.axial)
            && (boundaries.top != Boundary::Reflective || boundaries.bottom != Boundary::Reflective)
        {
            return Err(invalid(
                "`axial` constants need the top and bottom sides reflective".into(),
            ));
        }
        let entries = file.material.unwrap_or_default();
        let self_shielding = !macroscopic && file.self_shielding.unwrap_or(true);
        // 飛程の表は自己遮蔽の減速計算でだけ使う
        if !self_shielding && file.track_table.is_some() {
            return Err(invalid(
                "`track-table` is used only by the self-shielding (a library without `self-shielding = false`)"
                    .into(),
            ));
        }
        if !self_shielding && file.ufg_tolerance.is_some() {
            return Err(invalid(
                "`ufg-tolerance` is used only by the self-shielding (a library without `self-shielding = false`)"
                    .into(),
            ));
        }
        let thinning = file.ufg_tolerance;
        if let Some(t) = thinning
            && !(t > 0.0 && t < 1.0)
        {
            return Err(invalid(format!(
                "`ufg-tolerance` must be between 0 and 1, got {t}"
            )));
        }
        if !self_shielding && (file.ufg_scattering.is_some() || file.ufg_free_gas.is_some()) {
            return Err(invalid(
                "`ufg-scattering` and `ufg-free-gas` are used only by the self-shielding (a library without `self-shielding = false`)"
                    .into(),
            ));
        }
        let scattering =
            parse_scattering(file.ufg_scattering.as_deref(), file.ufg_free_gas).map_err(invalid)?;
        let collision = match file.track_table {
            None => Collision::default(),
            Some(TrackTableEntry::Exact(false)) => Collision::Exact,
            Some(TrackTableEntry::Density(density)) if density > 0 => Collision::Table {
                density: density as usize,
            },
            Some(other) => {
                return Err(invalid(format!(
                    "`track-table` must be a positive number of points per decade or false, not {other:?}"
                )));
            }
        };
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
            .map(|b| {
                b.build(
                    &materials,
                    &waters,
                    &mixtures,
                    &media,
                    file.history.is_some(),
                )
            })
            .transpose()
            .map_err(invalid)?;
        // 細かい群の媒質のファイルは、軸の直積の状態を持つ（ADR-0128）
        if file.medium.is_some() && branches.as_ref().is_some_and(|b| !b.is_product()) {
            return Err(invalid(
                "medium needs the branches to solve the product of their axes (no subgrids. ADR-0128)".into(),
            ));
        }
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
        let depletion = file
            .depletion
            .map(|d| {
                // 燃焼点からの分岐は履歴を読む別の入力（`[history]`。ADR-0187）で解く。媒質は決めていない（ADR-0173）。燃焼点ごとのピン出力は
                // 燃焼の履歴に入る（ADR-0181）
                if macroscopic || branches.is_some() {
                    return Err("depletion needs a library and no branches (ADR-0173)".to_string());
                }
                if !self_shielding {
                    return Err("depletion needs the self-shielding (ADR-0173)".into());
                }
                if file.pin_powers.is_some() || file.medium.is_some() {
                    return Err(
                        "depletion writes neither `pin-powers` (the pin powers of each point go into the depletion \
                         history) nor `medium`"
                            .into(),
                    );
                }
                // 燃焼点ごとの定数は燃焼の履歴に入る（`[depletion]` の `few-group`）。ファイルの `few-group` は書かれないので受けない
                if output.is_some() {
                    return Err(
                        "with [depletion], give the few-group bounds in [depletion] (`few-group = { bounds = [...] }`); \
                         the constants of each point go into the depletion history"
                            .into(),
                    );
                }
                let few_group = match d.few_group {
                    Some(f) => {
                        if f.bounds.len() < 2 {
                            return Err("the depletion few-group bounds need at least two energies".into());
                        }
                        let pins = f.pins.unwrap_or(false);
                        // 格子全体の形状関数は外周がすべて鏡面反射の格子でだけ作る（ADR-0096）
                        if pins && boundaries != Boundaries::REFLECTIVE {
                            return Err("depletion `pins` need a lattice with every side reflective".into());
                        }
                        Some(DepletionFewGroup {
                            bounds: f.bounds,
                            pins,
                        })
                    }
                    None => None,
                };
                if !positive(d.power) {
                    return Err(format!("the depletion power must be positive, got {}", d.power));
                }
                if d.burnup.is_empty()
                    || !d.burnup.iter().all(|&b| positive(b))
                    || !d.burnup.windows(2).all(|w| w[0] < w[1])
                {
                    return Err(format!(
                        "the burnup points must be positive and increasing, got {:?}",
                        d.burnup
                    ));
                }
                if d.materials.is_empty() {
                    return Err("depletion needs the materials to deplete".into());
                }
                for (k, name) in d.materials.iter().enumerate() {
                    if d.materials[..k].contains(name) {
                        return Err(format!("depleted material `{name}` is listed twice"));
                    }
                    if !materials.iter().any(|m| &m.name == name) {
                        return Err(format!(
                            "depleted material `{name}` is not a material written with nuclides or water"
                        ));
                    }
                    // 燃焼の領域の物質は `<物質>@<セル>` の名前で解く
                    if name.contains('@') {
                        return Err(format!("depleted material `{name}` must not contain `@`"));
                    }
                }
                Ok(DepletionInput {
                    chain: d.chain,
                    power: d.power,
                    burnup: d.burnup,
                    materials: d.materials,
                    output: d.output,
                    few_group,
                    shielding: match d.shielding {
                        None | Some(ShieldingEntry::EveryTransport) => {
                            ShieldingFrequency::EveryTransport
                        }
                        Some(ShieldingEntry::StepStart) => ShieldingFrequency::StepStart,
                    },
                })
            })
            .transpose()
            .map_err(invalid)?;
        let history = file
            .history
            .map(|h| {
                // core は 3 つの分岐の軸と燃焼度で引く。燃焼の輸送と同じく自己遮蔽で解き、格子全体の 1 つのテーブルを書く
                if branches.is_none() {
                    return Err("[history] needs [branches] (ADR-0187)".to_string());
                }
                if !self_shielding {
                    return Err("[history] needs the self-shielding (ADR-0187)".into());
                }
                if depletion.is_some() || file.medium.is_some() || !media.is_empty() || !mixtures.is_empty() {
                    return Err(
                        "[history] excludes [depletion], `medium`, media and mixtures (ADR-0187)".into(),
                    );
                }
                if !matches!(
                    &output,
                    Some(FewGroupOutput {
                        files: FewGroupFiles::Lattice(_),
                        axial: false,
                        ..
                    })
                ) {
                    return Err(
                        "[history] writes one table of the lattice: give `few-group = { output = ..., bounds = [...] }`"
                            .into(),
                    );
                }
                if h.burnup.is_empty()
                    || !h.burnup.iter().all(|&b| b >= 0.0 && b.is_finite())
                    || !h.burnup.windows(2).all(|w| w[0] < w[1])
                {
                    return Err(format!(
                        "the history burnup points must be non-negative and increasing, got {:?}",
                        h.burnup
                    ));
                }
                // 公称でない状態を解く燃焼点。書かなければすべての燃焼点（ADR-0200 の案 B1）
                let off_nominal = h.off_nominal_burnup.unwrap_or_else(|| h.burnup.clone());
                // core は公称の状態との差を外挿しないので、最初と最後の燃焼点を含む
                if off_nominal.first() != h.burnup.first()
                    || off_nominal.last() != h.burnup.last()
                    || !off_nominal.iter().all(|b| h.burnup.contains(b))
                    || !off_nominal.windows(2).all(|w| w[0] < w[1])
                {
                    return Err(format!(
                        "the off-nominal burnup points must be increasing points of `burnup` from its first to its \
                         last, got {off_nominal:?}"
                    ));
                }
                Ok(HistoryInput {
                    file: h.file,
                    burnup: h.burnup,
                    off_nominal,
                })
            })
            .transpose()
            .map_err(invalid)?;
        Ok(Self {
            self_shielding,
            collision,
            thinning,
            scattering,
            cross_sections,
            output,
            pin_powers: file.pin_powers,
            medium: file.medium,
            branches,
            depletion,
            history,
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
    if parts.iter().all(|p| p.material == VOID) {
        return Err(format!(
            "material `{name}`: a mixture of `{VOID}` alone is not a material; write `{VOID}` in the zones"
        ));
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
mod tests;
