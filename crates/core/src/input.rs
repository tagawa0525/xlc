//! 炉心の入力ファイル（TOML）。変わる頻度で4つのファイルに分ける（ADR-0105）。実行するのはケースのファイルで、
//! ケースがサイクルを、サイクルが炉心の形と台帳を指す。パスは、それを書いたファイルのあるディレクトリから辿る。
//!
//! 炉心の形（プラントでほとんど変わらない）:
//!
//! ```toml
//! pitch = 21.5                 # 集合体のピッチ (cm)
//! height = 406.3               # 炉心の高さ (cm)
//! map = ["*  R"]               # 上の行が +y。`*` は集合体を装荷する位置、`-` はノードを置かない位置（炉心の外形の外。
//!                              # 図の端まで `-` だけをたどって出られること。ADR-0132）、ほかは反射体の種類（`@90` で回す）
//! [reflector.R]
//! bounds = [0.0, 50.0, 406.3]  # 軸方向の層の境界 (cm)。0 から炉心の高さまで
//! # 層ごとの少数群定数。径方向の真空の比 r を持たない定数（軸方向の反射体のスラブの定数）の層は、真空の面を
//! # Marshak の条件にすることを `vacuum = "marshak"` で明示する（ADR-0132）
//! layers = [{ file = "water.msgpack", vacuum = "marshak" }, "reflector.msgpack"]
//! # 制御棒（ADR-0127。棒のない炉心は書かない）。先端の高さは bottom + ステップ × step (cm)
//! [rods]
//! step = 1.5875                # 1 ステップの長さ (cm)
//! bottom = 17.031              # 0 ステップの先端の高さ (cm)
//! steps = 230                  # 最大のステップ
//! top = 448.711                # 棒の区間を切る高さ (cm)。上の層のテーブルは制御棒の軸を持たない（ADR-0128）
//! map = [".  ."]               # `map` と同じ行と列。`*` の位置にバンクの名前か `.`
//! [rods.bank.D]
//! type = "rcca"
//! [rods.type.rcca]             # 先端から上への区間: テーブルの制御棒の軸の値と長さ (cm)
//! segments = [{ value = "aic", length = 101.6 }, { value = "b4c", length = 259.08 }]
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
//! radial-nodes = 1   # 集合体1体の各方向のノードの数。テーブルの分割の数と同じにする（反射体のある炉心は 1。ADR-0108）
//! # 各軸の（低い側, 高い側）の境界条件（reflective、zero-flux、vacuum、marshak のいずれか）。6面とも書く。
//! # vacuum は少数群定数の真空の辺の比 J/φ を使い、marshak は P1 の関係（J/φ = 1/2）を置く（ADR-0128）。
//! # outline は `map` の `-` に面する面の条件で、`-` があるときだけ、同じ値のどれかを書く（ADR-0132）
//! boundary = { x = ["vacuum", "vacuum"], y = ["vacuum", "vacuum"], z = ["vacuum", "vacuum"] }
//! # map に `-` があれば、boundary に outline = "vacuum" のように足す
//! # 再構成したピン出力（省略すると書かない。どのノードも形状関数を持つこと。ADR-0096）。状態が1つ以下のときだけ
//! # 書ける（この例は状態が2つなので書かない）
//! # pin-powers = "pins.msgpack"
//! # 核分裂のある層ごとのピン出力と軸方向の出力（省略すると書かない）。pin-powers と同じく、状態が1つ以下のときだけ
//! # 書ける
//! # axial-pin-powers = "layers.msgpack"
//! # 名前つきの炉心の状態。状態ごとに分岐のテーブルを内挿して解き、名前の順に出す（ADR-0097、ADR-0146）。
//! # 省略すると、軸のないテーブルだけを使って1つ解く
//! [state.aro]
//! boron = 1291        # ホウ素濃度 (ppm)
//! b10 = 0.199         # ホウ素の中の B-10 の原子分率（テーブルと同じ値）
//! moderator = 565     # 減速材の温度 (K)
//! fuel = 565          # 燃料の温度 (K)
//! # burnup = 10.0     # 燃焼度 (MWd/kg)。炉心のすべての集合体に同じ値。燃焼度の軸を持つテーブルで要る（ADR-0187）
//! [state.critical]    # boron の代わりに、k が目標になるホウ素濃度を探す（ADR-0145）
//! b10 = 0.199
//! moderator = 565
//! fuel = 565
//! # 範囲 (ppm。両端の k が目標を挟む)、目標の k、|k − 目標| の許容、固有値計算の回数の上限（両端を含む）
//! search = { boron = [1100.0, 1400.0], k = 1.0, tolerance = 1.0e-6, max-solves = 20 }
//! # 炉心の形に制御棒があれば、どの状態にもバンクの位置（ステップ）をすべて書く（ADR-0127）:
//! # banks = { D = 167 }
//! # 部分挿入のノードの扱い（ADR-0127。省略すると中性子束の重みと ADR-0127 の細かいメッシュの幅）
//! [partial-rods]
//! weighting = "flux"   # flux（ADR-0127）か volume（採らなかった体積の重み。比べるときだけで、fine-mesh は書かない）
//! fine-mesh = 0.25     # 細かいメッシュの幅の上限 (cm)
//! ```
//!
//! 燃焼（ADR-0190）のケースは `[depletion]` を持ち、名前つきの状態（`burnup` を書かないもの）をちょうど 1 つ持つ。その状態を
//! すべてのステップに使い、ノードごとの燃焼度を 0 から進める。状態が `search` を持てば、点と中点の解のすべてで臨界ホウ素を
//! 探す（ADR-0191）。`radial-nodes` は 1 で、`pin-powers` と `axial-pin-powers` は書けない:
//!
//! ```toml
//! [depletion]
//! specific-power = 38.41                      # 定格の出力密度 (W/g、炉心の平均、初めの重金属)
//! steps = [{ efpd = 9.0 }, { efpd = 32.0 }]   # ステップの終わりの累積の EFPD。0 より大きく増える
//! history = "history.msgpack"                 # 燃焼の履歴の出力
//! ```

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Component, Path, PathBuf};

use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use xlc_format::fewgroup::{FewGroup, Pins};
use xlc_format::table::ROD_OUT;

use crate::nem::{Boundary, Mesh, Problem, Solution, fission_rates};
use crate::rods::{PartialNode, Segment, Weighting};
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
    /// 層ごとに、真空の面を Marshak の条件（r = 1/2）にするか。径方向の真空の比 r を持たない定数（軸方向の
    /// 反射体のスラブの定数）を、真空の面に接する位置に置くときに入力で明示する（ADR-0132）
    pub marshak: Vec<bool>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ColumnEntry {
    bounds: Vec<f64>,
    layers: Vec<LayerEntry>,
}

/// 層: ファイルのパスか、`{ file = "...", vacuum = "marshak" }`。`#[serde(untagged)]` は書き誤りの原因（どの値が
/// 悪いか）を捨てるので、文字列かテーブルかで分けて読む
#[derive(Deserialize)]
#[serde(try_from = "toml::Value")]
enum LayerEntry {
    File(PathBuf),
    Table(LayerTable),
}

impl TryFrom<toml::Value> for LayerEntry {
    type Error = String;

    fn try_from(value: toml::Value) -> Result<Self, String> {
        match value {
            toml::Value::String(path) => Ok(Self::File(PathBuf::from(path))),
            table @ toml::Value::Table(_) => table
                .try_into()
                .map(Self::Table)
                .map_err(|e| format!("layer {{ file, vacuum = \"marshak\" }}: {e}")),
            other => Err(format!(
                "a layer is a file path or {{ file = \"...\", vacuum = \"marshak\" }}, got {other}"
            )),
        }
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct LayerTable {
    file: PathBuf,
    vacuum: VacuumRule,
}

/// 層の真空の面の条件。書けるのは Marshak だけ（真空の比 r を使うのは書かないとき）
#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
enum VacuumRule {
    Marshak,
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
        let (layers, marshak) = entry
            .layers
            .into_iter()
            .map(|l| match l {
                LayerEntry::File(file) => (normalized(&base.join(file)), false),
                LayerEntry::Table(LayerTable {
                    file,
                    vacuum: VacuumRule::Marshak,
                }) => (normalized(&base.join(file)), true),
            })
            .unzip();
        Ok(Self {
            bounds: entry.bounds,
            layers,
            marshak,
        })
    }

    fn height(&self) -> f64 {
        self.bounds.last().copied().unwrap_or(0.0)
    }
}

/// `.` と、前の名前を打ち消す `..` を字面で取り除く（同じファイルを別の綴りで書いても1つにまとめるため）。
/// symlinkはたどらない。
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
    /// ノードを置かない位置（`-`。炉心の外形の外。ADR-0132）
    Empty,
}

/// 炉心の形（ADR-0105）。
#[derive(Debug, Clone, PartialEq)]
pub struct Plant {
    pub pitch: f64,
    pub height: f64,
    pub reflectors: BTreeMap<String, Column>,
    /// 上の行（+y）から
    pub map: Vec<Vec<Site>>,
    /// 制御棒の位置、バンク、型（ADR-0127）。棒のない炉心は `None`
    pub rods: Option<Rods>,
}

/// 炉心の形の制御棒（ADR-0127）。位置 (cm) は下の炉心の端からの高さ。
#[derive(Debug, Clone, PartialEq)]
pub struct Rods {
    /// 1 ステップの長さ (cm)
    pub step: f64,
    /// 0 ステップの先端の高さ (cm)
    pub bottom: f64,
    /// 最大のステップ（全引き抜き）
    pub steps: u32,
    /// 棒の区間を切る高さ (cm)。その上の層のテーブルは制御棒の軸を持たない（ADR-0128 の「制御棒の軸」）
    pub top: f64,
    /// 上の行（+y）から、位置ごとのバンクの名前（棒のない位置は `None`）
    pub map: Vec<Vec<Option<String>>>,
    /// バンク → 型の名前
    pub banks: BTreeMap<String, String>,
    /// 型 → 先端から上への区間（制御棒の軸の値の名前、長さ (cm)）
    pub types: BTreeMap<String, Vec<(String, f64)>>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct PlantFile {
    pitch: f64,
    height: f64,
    map: Vec<String>,
    #[serde(default)]
    reflector: BTreeMap<String, ColumnEntry>,
    rods: Option<RodsFile>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RodsFile {
    step: f64,
    bottom: f64,
    steps: u32,
    top: f64,
    map: Vec<String>,
    bank: BTreeMap<String, BankEntry>,
    #[serde(rename = "type")]
    types: BTreeMap<String, RodTypeEntry>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BankEntry {
    #[serde(rename = "type")]
    kind: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RodTypeEntry {
    segments: Vec<RodSegmentEntry>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RodSegmentEntry {
    value: String,
    length: f64,
}

impl Rods {
    /// 制御棒の入力を確かめて読む。`sites` は炉心の形の図
    fn new(file: RodsFile, sites: &[Vec<Site>]) -> Result<Self, InputError> {
        if !(file.step > 0.0 && file.step.is_finite()) || file.steps == 0 {
            return invalid(format!(
                "rods: step must be positive and steps at least 1, got {} cm and {}",
                file.step, file.steps
            ));
        }
        if !(file.bottom.is_finite() && file.top.is_finite() && file.top > file.bottom) {
            return invalid(format!(
                "rods: bottom and top must be finite with top above bottom, got {} and {}",
                file.bottom, file.top
            ));
        }
        let mut types = BTreeMap::new();
        for (name, entry) in file.types {
            check_name("rod type", &name)?;
            let segments: Vec<(String, f64)> = entry
                .segments
                .into_iter()
                .map(|s| (s.value, s.length))
                .collect();
            if segments.is_empty()
                || segments.iter().any(|(value, length)| {
                    !(*length > 0.0 && length.is_finite())
                        || value.is_empty()
                        || value == ROD_OUT
                        || value.chars().any(char::is_whitespace)
                })
            {
                return invalid(format!(
                    "rod type {name}: needs segments of positive lengths whose values are names \
                     other than {ROD_OUT:?}"
                ));
            }
            types.insert(name, segments);
        }
        let mut banks = BTreeMap::new();
        for (name, entry) in file.bank {
            check_name("bank", &name)?;
            if !types.contains_key(&entry.kind) {
                return invalid(format!("bank {name}: unknown rod type {:?}", entry.kind));
            }
            banks.insert(name, entry.kind);
        }
        let cells = split_rows("rods map", &file.map)?;
        if cells.len() != sites.len() || cells.iter().zip(sites).any(|(c, s)| c.len() != s.len()) {
            return invalid("rods map must have the rows and columns of the plant map");
        }
        let mut used = HashSet::new();
        let map = cells
            .iter()
            .zip(sites)
            .enumerate()
            .map(|(r, (row, sites))| {
                row.iter()
                    .zip(sites)
                    .enumerate()
                    .map(|(c, (&cell, site))| match (cell, site) {
                        (".", _) => Ok(None),
                        (bank, Site::Assembly) if banks.contains_key(bank) => {
                            used.insert(bank.to_string());
                            Ok(Some(bank.to_string()))
                        }
                        (bank, Site::Assembly) => invalid(format!(
                            "rods map row {} column {}: unknown bank {bank:?}",
                            r + 1,
                            c + 1
                        )),
                        (bank, _) => invalid(format!(
                            "rods map row {} column {}: bank {bank:?} is not at an assembly site (`*`)",
                            r + 1,
                            c + 1
                        )),
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        if let Some(bank) = banks.keys().find(|b| !used.contains(*b)) {
            return invalid(format!("bank {bank} has no position in the rods map"));
        }
        Ok(Self {
            step: file.step,
            bottom: file.bottom,
            steps: file.steps,
            top: file.top,
            map,
            banks,
            types,
        })
    }
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
                        "-" => Ok(Site::Empty),
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
        check_open(&map)?;
        let rods = file.rods.map(|r| Rods::new(r, &map)).transpose()?;
        Ok(Self {
            pitch: file.pitch,
            height: file.height,
            reflectors,
            map,
            rods,
        })
    }
}

/// どの `-` も、図の端まで `-` だけを上下左右の隣でたどって出られることと、ノードのある位置がひとつながりであることを
/// 確かめる。ノードに囲まれた `-` は、炉心の中に真空の穴を開けることになり、意図した入力と考えにくい（ADR-0132）
fn check_open(map: &[Vec<Site>]) -> Result<(), InputError> {
    let rows = map.len();
    let cols = map.first().map_or(0, Vec::len);
    let empty = |r: usize, c: usize| map[r][c] == Site::Empty;
    // 図の端の `-` から、`-` をたどって届く位置
    let mut reached = vec![vec![false; cols]; rows];
    let mut stack: Vec<(usize, usize)> = (0..rows)
        .flat_map(|r| (0..cols).map(move |c| (r, c)))
        .filter(|&(r, c)| (r == 0 || c == 0 || r + 1 == rows || c + 1 == cols) && empty(r, c))
        .collect();
    while let Some((r, c)) = stack.pop() {
        if reached[r][c] {
            continue;
        }
        reached[r][c] = true;
        let beside = [
            r.checked_sub(1).map(|r| (r, c)),
            (r + 1 < rows).then_some((r + 1, c)),
            c.checked_sub(1).map(|c| (r, c)),
            (c + 1 < cols).then_some((r, c + 1)),
        ];
        stack.extend(beside.into_iter().flatten().filter(|&(r, c)| empty(r, c)));
    }
    for (r, row) in map.iter().enumerate() {
        for (c, site) in row.iter().enumerate() {
            if *site == Site::Empty && !reached[r][c] {
                return invalid(format!(
                    "map row {} column {}: `-` is enclosed by nodes; it would open a vacuum hole \
                     in the core",
                    r + 1,
                    c + 1
                ));
            }
        }
    }
    // ノードのある位置が、上下左右の隣でひとつながりであること。`-` で切り離されたノードは別の問題になる
    let mut joined = vec![vec![false; cols]; rows];
    let start = (0..rows)
        .flat_map(|r| (0..cols).map(move |c| (r, c)))
        .find(|&(r, c)| !empty(r, c));
    let Some(start) = start else {
        return invalid("map: no position has a node");
    };
    let mut stack = vec![start];
    while let Some((r, c)) = stack.pop() {
        if joined[r][c] {
            continue;
        }
        joined[r][c] = true;
        let beside = [
            r.checked_sub(1).map(|r| (r, c)),
            (r + 1 < rows).then_some((r + 1, c)),
            c.checked_sub(1).map(|c| (r, c)),
            (c + 1 < cols).then_some((r, c + 1)),
        ];
        stack.extend(beside.into_iter().flatten().filter(|&(r, c)| !empty(r, c)));
    }
    for (r, row) in joined.iter().enumerate() {
        for (c, &joined) in row.iter().enumerate() {
            if !empty(r, c) && !joined {
                return invalid(format!(
                    "map row {} column {}: the node is cut off from the rest of the core by `-`",
                    r + 1,
                    c + 1
                ));
            }
        }
    }
    Ok(())
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
    /// 炉心の形の `-` に面する面の条件（ADR-0132）
    pub outline: Option<Boundary>,
    /// 名前つきの状態（`[state.<名前>]`。名前の順）。状態ごとに分岐のテーブルを内挿して解く（ADR-0097）。空なら
    /// 軸のないテーブルだけを使って1つ解く（ADR-0146）
    pub states: Vec<CaseState>,
    /// 再構成したピン出力を書くファイル（ADR-0096）。状態が1つ以下のときだけ
    pub pin_powers: Option<PathBuf>,
    /// 部分挿入のノードの扱い（`[partial-rods]`）。書かなければ ADR-0127 の値
    pub partial_rods: PartialRods,
    /// `[partial-rods]` を書いたか（棒のない炉心に書けばエラーにする）
    pub partial_rods_given: bool,
    /// 核分裂のある層ごとのピン出力と軸方向の出力を書くファイル。状態が1つ以下のときだけ
    pub axial_pin_powers: Option<PathBuf>,
    /// 燃焼のステップ（ADR-0190）
    pub depletion: Option<DepletionInput>,
}

/// ケースの `[depletion]`（ADR-0190）。
#[derive(Debug, Clone, PartialEq)]
pub struct DepletionInput {
    /// 定格の出力密度 (W/g、炉心の平均、初めの重金属)
    pub specific_power: f64,
    /// ステップ（0 より大きく増える累積の EFPD の順）
    pub steps: Vec<DepletionStep>,
    /// 燃焼の履歴を書くファイル
    pub history: PathBuf,
}

/// 燃焼のステップ。ステップの中の出力は定格で一定（ADR-0190）
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DepletionStep {
    /// ステップの終わりの累積の EFPD
    pub efpd: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct DepletionEntry {
    specific_power: f64,
    steps: Vec<StepEntry>,
    history: PathBuf,
}

/// 項目は、効くようになるときに足す（出力の割合、入口の温度など。ADR-0190）
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct StepEntry {
    efpd: f64,
}

impl DepletionInput {
    fn new(entry: DepletionEntry, base: &Path) -> Result<Self, InputError> {
        let p = entry.specific_power;
        if !(p > 0.0 && p.is_finite()) {
            return invalid(format!(
                "depletion: specific-power must be positive, got {p}"
            ));
        }
        let efpd: Vec<f64> = entry.steps.iter().map(|s| s.efpd).collect();
        if efpd.is_empty()
            || !efpd.iter().all(|e| e.is_finite())
            || efpd[0] <= 0.0
            || !efpd.windows(2).all(|w| w[0] < w[1])
        {
            return invalid(format!(
                "depletion: the step efpd must increase from above 0, got {efpd:?}"
            ));
        }
        Ok(Self {
            specific_power: p,
            steps: efpd
                .into_iter()
                .map(|efpd| DepletionStep { efpd })
                .collect(),
            history: normalized(&base.join(entry.history)),
        })
    }

    /// ケースの状態と出力の組み合わせを確かめる
    fn check(
        &self,
        states: &[CaseState],
        pin_powers: bool,
        radial_nodes: usize,
    ) -> Result<(), InputError> {
        let [state] = states else {
            return invalid(format!(
                "depletion needs exactly one state, got {}",
                states.len()
            ));
        };
        if state.burnup.is_some() {
            return invalid(format!(
                "state {}: give no burnup; the depletion advances the node burnups (ADR-0190)",
                state.name
            ));
        }
        if pin_powers {
            return invalid(
                "depletion: pin-powers and axial-pin-powers are not written yet (ADR-0190)",
            );
        }
        // テーブルの重金属の質量は格子全体の 1 つの値で、部分ノードごとの質量を持たない
        if radial_nodes != 1 {
            return invalid(format!(
                "depletion needs radial-nodes = 1, got {radial_nodes} (ADR-0190)"
            ));
        }
        Ok(())
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct CaseFile {
    cycle: PathBuf,
    axial_mesh: Vec<f64>,
    radial_nodes: usize,
    boundary: BoundaryEntry,
    state: Option<BTreeMap<String, StateEntry>>,
    pin_powers: Option<PathBuf>,
    partial_rods: Option<PartialRodsEntry>,
    axial_pin_powers: Option<PathBuf>,
    depletion: Option<DepletionEntry>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct StateEntry {
    boron: Option<f64>,
    b10: f64,
    moderator: f64,
    fuel: f64,
    burnup: Option<f64>,
    search: Option<SearchEntry>,
    /// バンク → 位置（ステップ）。棒のない炉心は書かない（`Layout::new` が確かめる）
    #[serde(default)]
    banks: BTreeMap<String, u32>,
}

/// 部分挿入のノードの細かいメッシュの幅 (cm)。ADR-0127 で収束を測って決めた値
pub const FINE_MESH: f64 = 0.5;

/// 部分挿入のノードの扱い（ADR-0127）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PartialRods {
    /// 区間の定数の混ぜ方。体積の重みは ADR-0127 で採らなかった案で、比べるときだけ書く
    pub weighting: Weighting,
    /// 細かいメッシュの幅の上限 (cm)
    pub fine_mesh: f64,
}

impl Default for PartialRods {
    fn default() -> Self {
        Self {
            weighting: Weighting::Flux,
            fine_mesh: FINE_MESH,
        }
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct PartialRodsEntry {
    weighting: Option<WeightingEntry>,
    fine_mesh: Option<f64>,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "kebab-case")]
enum WeightingEntry {
    Flux,
    Volume,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct SearchEntry {
    boron: [f64; 2],
    k: f64,
    tolerance: f64,
    max_solves: usize,
}

impl CaseState {
    fn new(name: String, entry: StateEntry) -> Result<Self, InputError> {
        check_name("state", &name)?;
        let what = |message: String| invalid(format!("state {name}: {message}"));
        let boron = match (entry.boron, entry.search) {
            (Some(ppm), None) => Boron::Given(ppm),
            (None, Some(search)) => {
                let [low, high] = search.boron;
                if !(low >= 0.0 && low < high && high.is_finite()) {
                    return what(format!(
                        "search boron must be a range from low to high ppm, got {:?}",
                        search.boron
                    ));
                }
                if !(search.k > 0.0 && search.k.is_finite()) {
                    return what(format!("search k must be positive, got {}", search.k));
                }
                if !(search.tolerance > 0.0 && search.tolerance.is_finite()) {
                    return what(format!(
                        "search tolerance must be positive, got {}",
                        search.tolerance
                    ));
                }
                // 範囲の両端だけでは探せない
                if search.max_solves < 3 {
                    return what(format!(
                        "search max-solves must be at least 3 (the two ends and a point between), got {}",
                        search.max_solves
                    ));
                }
                Boron::Search(BoronSearch {
                    bracket: search.boron,
                    k: search.k,
                    tolerance: search.tolerance,
                    max_solves: search.max_solves,
                })
            }
            (Some(_), Some(_)) => return what("give boron or search, not both".into()),
            (None, None) => return what("give boron or search".into()),
        };
        if let Some(b) = entry.burnup
            && !(b >= 0.0 && b.is_finite())
        {
            return what(format!("burnup must be non-negative, got {b}"));
        }
        Ok(Self {
            name,
            boron,
            b10: entry.b10,
            moderator: entry.moderator,
            fuel: entry.fuel,
            burnup: entry.burnup,
            banks: entry.banks,
        })
    }
}

/// 各軸の（低い側, 高い側）
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BoundaryEntry {
    x: [Condition; 2],
    y: [Condition; 2],
    z: [Condition; 2],
    outline: Option<Condition>,
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
        // 名前のない `[state]` に項目を書いた入力は、項目の名前が状態の名前に見え、型のエラーでは修正方法が分からないので、
        // 先に見つけて示す（状態の中身のエラーは、行を示す toml の読み込みに任せる）
        let raw: toml::Table = parse_toml(text)?;
        if let Some(toml::Value::Table(states)) = raw.get("state")
            && let Some((name, _)) = states.iter().find(|(_, v)| !v.is_table())
        {
            return invalid(format!(
                "[state] has a value {name:?}; name each state as [state.<name>] (ADR-0146)"
            ));
        }
        let file: CaseFile = parse_toml(text)?;
        check_heights("axial-mesh", &file.axial_mesh)?;
        if file.radial_nodes == 0 {
            return invalid("radial-nodes must be positive");
        }
        let b = file.boundary;
        let condition = |c: Condition| match c {
            Condition::Reflective => Boundary::Reflective,
            Condition::ZeroFlux => Boundary::ZeroFlux,
            Condition::Vacuum => Boundary::Vacuum,
            Condition::Marshak => Boundary::Marshak,
        };
        let boundaries = [b.x[0], b.x[1], b.y[0], b.y[1], b.z[0], b.z[1]].map(condition);
        // BTreeMap なので名前の順
        let states = match file.state {
            None => Vec::new(),
            Some(states) if states.is_empty() => {
                return invalid("[state] has no states; name each as [state.<name>]");
            }
            Some(states) => states
                .into_iter()
                .map(|(name, entry)| CaseState::new(name, entry))
                .collect::<Result<Vec<_>, InputError>>()?,
        };
        let mut partial_rods = PartialRods::default();
        let partial_rods_given = file.partial_rods.is_some();
        if let Some(entry) = file.partial_rods {
            if let Some(w) = entry.weighting {
                partial_rods.weighting = match w {
                    WeightingEntry::Flux => Weighting::Flux,
                    WeightingEntry::Volume => Weighting::Volume,
                };
            }
            if let Some(fine) = entry.fine_mesh {
                if !(fine > 0.0 && fine.is_finite()) {
                    return invalid(format!(
                        "partial-rods: fine-mesh must be positive, got {fine}"
                    ));
                }
                // 体積の重みは細かいメッシュを解かないので、幅を書いても結果は変わらない
                if partial_rods.weighting == Weighting::Volume {
                    return invalid(
                        "partial-rods: fine-mesh has no effect with weighting = \"volume\"",
                    );
                }
                partial_rods.fine_mesh = fine;
            }
        }
        let depletion = match file.depletion {
            Some(entry) => {
                let d = DepletionInput::new(entry, base)?;
                d.check(
                    &states,
                    file.pin_powers.is_some() || file.axial_pin_powers.is_some(),
                    file.radial_nodes,
                )?;
                Some(d)
            }
            None => None,
        };
        if file.pin_powers.is_some() && states.len() > 1 {
            return invalid(format!(
                "pin-powers needs at most one state, got {}",
                states.len()
            ));
        }
        if file.axial_pin_powers.is_some() && states.len() > 1 {
            return invalid(format!(
                "axial-pin-powers needs at most one state, got {}",
                states.len()
            ));
        }
        let pin_powers = file.pin_powers.map(|p| normalized(&base.join(p)));
        let axial_pin_powers = file.axial_pin_powers.map(|p| normalized(&base.join(p)));
        if let (Some(a), Some(b)) = (&pin_powers, &axial_pin_powers)
            && a == b
        {
            return invalid(format!(
                "pin-powers and axial-pin-powers name the same file {}",
                a.display()
            ));
        }
        Ok(Self {
            cycle: base.join(file.cycle),
            axial_mesh: file.axial_mesh,
            radial_nodes: file.radial_nodes,
            boundaries,
            outline: b.outline.map(condition),
            states,
            pin_powers,
            partial_rods,
            partial_rods_given,
            axial_pin_powers,
            depletion,
        })
    }
}

/// ケースの名前つきの状態（ADR-0146）。
#[derive(Debug, Clone, PartialEq)]
pub struct CaseState {
    pub name: String,
    pub boron: Boron,
    /// ホウ素の中の B-10 の原子分率（テーブルと同じ値）
    pub b10: f64,
    /// 減速材の温度 (K)
    pub moderator: f64,
    /// 燃料の温度 (K)
    pub fuel: f64,
    /// 燃焼度 (MWd/kg)。炉心のすべての集合体に同じ値（ADR-0187）。燃焼していない炉心は `None`
    pub burnup: Option<f64>,
    /// バンク → 位置（ステップ）。炉心の形に制御棒があれば、そのバンクをすべて書く（ADR-0127）
    pub banks: BTreeMap<String, u32>,
}

/// 状態のホウ素濃度: 与えた値か、探索（ADR-0145）。
#[derive(Debug, Clone, PartialEq)]
pub enum Boron {
    /// ppm
    Given(f64),
    Search(BoronSearch),
}

/// k が目標になるホウ素濃度の探索（ADR-0145）。燃焼の履歴には入力と同じ名前で書く（ADR-0191）
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BoronSearch {
    /// 探す範囲 (ppm)。下端 < 上端
    #[serde(rename = "boron")]
    pub bracket: [f64; 2],
    /// 目標の k
    pub k: f64,
    /// |k − 目標| の許容
    pub tolerance: f64,
    /// 固有値計算の回数の上限（範囲の両端を含む）
    pub max_solves: usize,
}

impl CaseState {
    /// ホウ素濃度を `boron` (ppm) にした、テーブルを内挿する状態。
    pub fn at(&self, boron: f64) -> State {
        State {
            boron,
            b10: self.b10,
            moderator: self.moderator,
            fuel: self.fuel,
            burnup: self.burnup,
        }
    }
}

/// ノードの物質: 少数群定数のファイルの番号と、反時計回りの 90° の回数。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LayoutMaterial {
    pub file: usize,
    /// テーブルの使い方（`Layout::tables`）の番号
    pub table: usize,
    pub turns: u8,
    /// テーブルの部分ノードの番号（回す前の並び。-y の行から、各行は -x から。ADR-0108）
    pub sub: usize,
    /// 炉心の形の反射体の位置に置いたか（同じファイルでも集合体の位置のものと分ける。ADR-0095）
    pub reflector: bool,
    /// 真空の面を Marshak の条件にする層か（`Column::marshak`）
    pub marshak: bool,
}

/// テーブルの使い方: ファイルと、制御棒の区間の値。
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TableUse {
    /// `Layout::files` の番号
    pub file: usize,
    /// 制御棒の区間の材料（制御棒の軸の値）。棒のない区間は `None`（軸があれば `out` を選ぶ）
    pub rod: Option<String>,
    /// 棒を切る高さ（`Rods::top`）より上の区間か。棒はあるが、層のテーブルは制御棒の軸を持たないこと
    pub cut: bool,
    /// 燃焼のケース（ADR-0190）の集合体の位置のノードは、ノードごとの燃焼度で引くので、その直方体の位置。ほかは `None`
    pub node: Option<usize>,
}

/// 部分挿入のノード（ADR-0127）。
#[derive(Debug, Clone, PartialEq)]
pub struct LayoutPartial {
    /// 直方体の位置
    pub site: usize,
    /// 下からの区間の（テーブルの使い方の番号, 長さ (cm)）
    pub segments: Vec<(usize, f64)>,
    pub turns: u8,
    pub sub: usize,
    pub marshak: bool,
}

/// 4つの入力から組み立てた炉心のノードの並び。
#[derive(Debug, Clone, PartialEq)]
pub struct Layout {
    pub mesh: Mesh,
    /// 少数群定数のファイル（現れた順）
    pub files: Vec<PathBuf>,
    /// テーブルの使い方（現れた順）。`problem` と `partial_nodes` の定数はこの並び
    pub tables: Vec<TableUse>,
    /// 部分挿入のノード。物質の番号は `materials` の後ろに、この順に続く
    pub partial: Vec<LayoutPartial>,
    /// ファイルごとに、反射体（の種類のどれか）に使う層か
    pub reflector: Vec<bool>,
    pub materials: Vec<LayoutMaterial>,
    /// 位置ごとの物質の番号。並びは x が最も速く、次に y（下から）、z。ノードのない位置（`-`）は `None`。
    /// 部分挿入のノードの番号は `materials.len()` 以上（`partial` の順）
    pub node_material: Vec<Option<usize>>,
    /// 集合体の位置を各方向に分けたノードの数（ケースの `radial-nodes`。ADR-0108）
    pub split: usize,
    /// 上の行（+y）から、位置ごとの集合体の ID（反射体の位置とノードのない位置は `None`）
    pub ids: Vec<Vec<Option<String>>>,
    /// 上の行（+y）から、ノードのない位置（`-`。ADR-0132）か
    pub empty: Vec<Vec<bool>>,
    pub boundaries: [Boundary; 6],
    /// `-` に面する面の条件（ADR-0132）
    pub outline: Option<Boundary>,
    /// 部分挿入のノードの扱い（ケースの `[partial-rods]`）
    pub partial_rods: PartialRods,
    inputs: Box<LayoutInputs>,
}

/// 装荷図の位置に置いたもの。
#[derive(Clone, Copy)]
struct Placement<'a> {
    column: &'a Column,
    /// 反時計回りの 90° の回数
    turns: u8,
    /// 反射体の種類か
    reflector: bool,
    /// 制御棒の先端の高さ (cm) と型の区間
    rod: Option<(f64, &'a [(String, f64)])>,
}

/// ノードの中の、制御棒の状態が一様な区間の値
#[derive(Debug, Clone, PartialEq)]
enum RodState {
    Out,
    In(String),
    Cut,
}

/// 高さ [z0, z1] のノードを、先端 `tip` から上に `segments` が並ぶ棒（`top` より上は切る）で、状態が一様な区間
/// （下から、長さ）に分ける。ノードの面から相対 1e-9 より近い境界は面に寄せる（先端の高さはステップの和なので、
/// 丸めで面からずれた境界が非常に細い区間を作らないように）
fn rod_intervals(
    tip: f64,
    segments: &[(String, f64)],
    top: f64,
    z0: f64,
    z1: f64,
) -> Vec<(RodState, f64)> {
    let eps = 1.0e-9 * (z1 - z0);
    let mut points = vec![tip];
    for (_, length) in segments {
        let last = points[points.len() - 1];
        points.push(last + length);
    }
    let state_at = |z: f64| {
        if z >= top {
            return RodState::Cut;
        }
        match points.windows(2).position(|w| w[0] <= z && z < w[1]) {
            Some(i) => RodState::In(segments[i].0.clone()),
            None => RodState::Out,
        }
    };
    let mut bounds = vec![z0];
    let mut inner: Vec<f64> = points
        .iter()
        .chain(std::iter::once(&top))
        .copied()
        .filter(|&p| p > z0 + eps && p < z1 - eps)
        .collect();
    inner.sort_by(f64::total_cmp);
    bounds.extend(inner);
    bounds.push(z1);
    let mut out: Vec<(RodState, f64)> = Vec::new();
    for w in bounds.windows(2) {
        let state = state_at(0.5 * (w[0] + w[1]));
        match out.last_mut() {
            Some((last, length)) if *last == state => *length += w[1] - w[0],
            _ => out.push((state, w[1] - w[0])),
        }
    }
    out
}

/// 炉心の形の制御棒とケースの状態のバンクの位置が合うことを確かめる: 棒のある炉心は状態を持ち、どの状態もバンクを
/// ちょうどすべて書き、ステップが範囲の中にある。棒のない炉心にバンクの位置や `[partial-rods]` を書けばエラー
fn check_banks(case: &Case, plant: &Plant) -> Result<(), InputError> {
    let Some(rods) = &plant.rods else {
        if let Some(state) = case.states.iter().find(|s| !s.banks.is_empty()) {
            return invalid(format!(
                "state {}: gives bank steps, but the plant has no rods",
                state.name
            ));
        }
        if case.partial_rods_given {
            return invalid("the case has [partial-rods], but the plant has no rods");
        }
        return Ok(());
    };
    if case.states.is_empty() {
        return invalid(
            "the plant has rods; give the core states with the bank steps \
             ([state.<name>] with banks = { name = steps })",
        );
    }
    for state in &case.states {
        check_state_banks(rods, &state.banks).map_err(|e| match e {
            InputError::Invalid(m) => InputError::Invalid(format!("state {}: {m}", state.name)),
            e => e,
        })?;
    }
    Ok(())
}

/// バンクの位置が炉心の形の制御棒のバンクをちょうどすべて書き、ステップが範囲の中にあること
fn check_state_banks(rods: &Rods, banks: &BTreeMap<String, u32>) -> Result<(), InputError> {
    if let Some(bank) = rods.banks.keys().find(|b| !banks.contains_key(*b)) {
        return invalid(format!("banks: no steps for bank {bank}"));
    }
    if let Some(bank) = banks.keys().find(|b| !rods.banks.contains_key(*b)) {
        return invalid(format!("banks: unknown bank {bank}"));
    }
    if let Some((bank, steps)) = banks.iter().find(|(_, s)| **s > rods.steps) {
        return invalid(format!(
            "banks: bank {bank} at {steps} steps is beyond the {} steps of the rods",
            rods.steps
        ));
    }
    Ok(())
}

/// 位置ごと（炉心の形の図の行と列）の、制御棒の先端の高さ (cm) と型の区間。`banks` が `None` なら棒を入れない
/// （`Layout::new` の並び。状態ごとの並びは `Layout::with_banks`）。棒の位置のバンクの位置か型がなければエラー
#[allow(clippy::type_complexity)]
fn rod_tips<'a>(
    plant: &'a Plant,
    banks: Option<&BTreeMap<String, u32>>,
) -> Result<Vec<Vec<Option<(f64, &'a [(String, f64)])>>>, InputError> {
    let rows = plant.map.len();
    let cols = plant.map.first().map_or(0, Vec::len);
    let (Some(rods), Some(banks)) = (&plant.rods, banks) else {
        return Ok(vec![vec![None; cols]; rows]);
    };
    rods.map
        .iter()
        .map(|row| {
            row.iter()
                .map(|bank| {
                    let Some(bank) = bank else {
                        return Ok(None);
                    };
                    let steps = banks
                        .get(bank)
                        .ok_or_else(|| InputError::Invalid(format!("bank {bank} has no step")))?;
                    let kind = rods
                        .banks
                        .get(bank)
                        .and_then(|t| rods.types.get(t))
                        .ok_or_else(|| {
                            InputError::Invalid(format!("bank {bank} has no rod type"))
                        })?;
                    Ok(Some((
                        rods.bottom + f64::from(*steps) * rods.step,
                        kind.as_slice(),
                    )))
                })
                .collect()
        })
        .collect()
}

/// `Layout::with_banks` が並びを作り直すための入力
#[derive(Debug, Clone, PartialEq)]
struct LayoutInputs {
    case: Case,
    cycle: Cycle,
    plant: Plant,
    inventory: Inventory,
}

impl Layout {
    /// 4つの入力からノードの並びを作る。制御棒は入れない（状態ごとのバンクの位置で `with_banks` が入れる）が、
    /// 炉心の形の制御棒とケースの状態のバンクの位置が合うことは、ここで確かめる
    pub fn new(
        case: &Case,
        cycle: &Cycle,
        plant: &Plant,
        inventory: &Inventory,
    ) -> Result<Self, InputError> {
        check_banks(case, plant)?;
        Self::build(case, cycle, plant, inventory, None)
    }

    /// バンクの位置 `banks`（ケースの状態の `banks`）で制御棒を入れた並び（ADR-0127）
    pub fn with_banks(&self, banks: &BTreeMap<String, u32>) -> Result<Self, InputError> {
        let LayoutInputs {
            case,
            cycle,
            plant,
            inventory,
        } = &*self.inputs;
        match &plant.rods {
            Some(rods) => check_state_banks(rods, banks)?,
            None if !banks.is_empty() => {
                return invalid("bank steps are given, but the plant has no rods");
            }
            // 棒のない炉心の並びは状態によらない
            None => return Ok(self.clone()),
        }
        Self::build(case, cycle, plant, inventory, Some(banks))
    }

    fn build(
        case: &Case,
        cycle: &Cycle,
        plant: &Plant,
        inventory: &Inventory,
        banks: Option<&BTreeMap<String, u32>>,
    ) -> Result<Self, InputError> {
        let mesh_height = case.axial_mesh.last().copied().unwrap_or(0.0);
        if mesh_height != plant.height {
            return invalid(format!(
                "axial-mesh ends at {mesh_height}, but the core height is {}",
                plant.height
            ));
        }
        // 制御棒の先端の高さ（位置ごと、上の行から）と型の区間
        let rod_tips = rod_tips(plant, banks)?;
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
        // ノードのない位置（`-`）は `None`
        let mut placements: Vec<Vec<Option<Placement>>> = Vec::with_capacity(rows);
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
                        row.push(Some(Placement {
                            column,
                            turns: p.turns,
                            reflector: false,
                            rod: rod_tips[r][c],
                        }));
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
                        row.push(Some(Placement {
                            column,
                            turns: *turns,
                            reflector: true,
                            rod: None,
                        }));
                        id_row.push(None);
                    }
                    (Site::Empty, None) => {
                        row.push(None);
                        id_row.push(None);
                    }
                    (Site::Empty, Some(p)) => {
                        return invalid(format!(
                            "{at}: {} is placed where the plant map has no node (`-`); use `.`",
                            p.id
                        ));
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

        let empty = placements.iter().flatten().any(Option::is_none);
        match (empty, case.outline) {
            (true, None) => {
                return invalid(
                    "the plant map has positions without nodes (`-`); give the condition of the faces \
                     next to them as `outline` in the case's boundary",
                );
            }
            (false, Some(_)) => {
                return invalid(
                    "the case's boundary has an outline condition, but the plant map has no `-`",
                );
            }
            _ => {}
        }

        // 反射体の定数は、分けた均質な問題の不連続因子をまだ作れない（ADR-0108 の次の段階）
        let split = case.radial_nodes;
        if split > 1 && placements.iter().flatten().flatten().any(|p| p.reflector) {
            return invalid(format!(
                "radial-nodes = {split}: cores with reflector sites cannot be split yet (ADR-0108)"
            ));
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

        // 棒を切る高さがノードの途中にあれば、1 つのノードに切った区間と切らない区間ができ、切った区間の層の
        // テーブルが制御棒の軸を持たないことと合わない
        if let Some(rods) = &plant.rods
            && rods.top > 0.0
            && rods.top < plant.height
            && !mesh.contains(&rods.top)
        {
            return invalid(format!(
                "rods: top {} cm is not on the axial-mesh",
                rods.top
            ));
        }

        let mut files: Vec<PathBuf> = Vec::new();
        let mut reflector: Vec<bool> = Vec::new();
        let mut file_index: HashMap<&Path, usize> = HashMap::new();
        let mut tables: Vec<TableUse> = Vec::new();
        let mut table_index: HashMap<TableUse, usize> = HashMap::new();
        let mut partial: Vec<LayoutPartial> = Vec::new();
        // 部分挿入のノードの、node_material の位置（物質の番号は物質を数え終えてから付ける）
        let mut pending: Vec<usize> = Vec::new();
        let mut materials = Vec::new();
        let mut material_index: HashMap<LayoutMaterial, usize> = HashMap::new();
        let mut node_material = Vec::with_capacity(rows * cols * split * split * (mesh.len() - 1));
        // 部分ノードの位置（炉心の向きで、-x と -y から）に来る、回す前の部分ノードの番号。90° 回すと元の (i, j) は
        // (n − 1 − j, i) に移るので、逆にたどる
        let sub_of = |mut i: usize, mut j: usize, turns: u8| {
            for _ in 0..turns % 4 {
                (i, j) = (j, split - 1 - i);
            }
            i + split * j
        };
        let sites = placements.iter().rev().flat_map(|row| {
            (0..split).flat_map(move |sj| {
                row.iter()
                    .flat_map(move |p| (0..split).map(move |si| (p, si, sj)))
            })
        });
        let sites: Vec<(&Option<Placement>, usize, usize)> = sites.collect();
        let last_plane = mesh.len() - 2;
        for (plane, z) in mesh.windows(2).enumerate() {
            // 図の上の行が +y なので、y は下の行から。集合体の位置の中は部分ノードの行と列
            for &(placement, si, sj) in &sites {
                let Some(Placement {
                    column,
                    turns,
                    reflector: is_reflector,
                    rod,
                }) = *placement
                else {
                    node_material.push(None);
                    continue;
                };
                // 層の境界はメッシュに乗っているので、ノードはどれか1つの層に収まる
                let layer = column
                    .bounds
                    .windows(2)
                    .position(|b| b[0] <= z[0] && z[1] <= b[1]);
                let Some((path, marshak)) =
                    layer.map(|l| (column.layers[l].as_path(), column.marshak[l]))
                else {
                    return invalid(format!("no layer covers {}..{} cm", z[0], z[1]));
                };
                let file = *file_index.entry(path).or_insert_with(|| {
                    files.push(path.to_path_buf());
                    reflector.push(false);
                    files.len() - 1
                });
                // 軸方向の反射体の層などは、反射体と集合体の設計で同じファイルを使ってよい
                reflector[file] |= is_reflector;
                let intervals = match rod {
                    Some((tip, segments)) => {
                        let top = plant.rods.as_ref().map_or(f64::INFINITY, |r| r.top);
                        rod_intervals(tip, segments, top, z[0], z[1])
                    }
                    None => vec![(RodState::Out, z[1] - z[0])],
                };
                let site = node_material.len();
                let mut use_of = |state: &RodState| {
                    let key = TableUse {
                        file,
                        rod: match state {
                            RodState::In(value) => Some(value.clone()),
                            RodState::Out | RodState::Cut => None,
                        },
                        cut: *state == RodState::Cut,
                        node: (case.depletion.is_some() && !is_reflector).then_some(site),
                    };
                    *table_index.entry(key.clone()).or_insert_with(|| {
                        tables.push(key);
                        tables.len() - 1
                    })
                };
                if let [_, _, ..] = intervals.as_slice() {
                    // 棒のない区間が真空側にある 1 ノードの問題は、z の面の不連続因子が負になって解けない（#0102）
                    let faces = [(0, 4), (last_plane, 5)];
                    if let Some(&(_, face)) = faces.iter().find(|&&(p, f)| {
                        p == plane
                            && matches!(case.boundaries[f], Boundary::Vacuum | Boundary::Marshak)
                    }) {
                        return invalid(format!(
                            "a partially rodded node ({}..{} cm) faces the {} z end of the core, \
                             which is a vacuum or Marshak boundary; put an axial reflector layer \
                             between the rods and that end (the solve stops on a negative z-face \
                             discontinuity factor there, #0102)",
                            z[0],
                            z[1],
                            if face == 4 { "lower" } else { "upper" }
                        ));
                    }
                    partial.push(LayoutPartial {
                        site: node_material.len(),
                        segments: intervals
                            .iter()
                            .map(|(state, length)| (use_of(state), *length))
                            .collect(),
                        turns,
                        sub: sub_of(si, sj, turns),
                        marshak,
                    });
                    pending.push(node_material.len());
                    node_material.push(None);
                    continue;
                }
                let table = use_of(&intervals[0].0);
                let material = LayoutMaterial {
                    file,
                    table,
                    turns,
                    sub: sub_of(si, sj, turns),
                    reflector: is_reflector,
                    marshak,
                };
                let m = *material_index.entry(material).or_insert_with(|| {
                    materials.push(material);
                    materials.len() - 1
                });
                node_material.push(Some(m));
            }
        }
        for (i, &site) in pending.iter().enumerate() {
            node_material[site] = Some(materials.len() + i);
        }

        Ok(Self {
            mesh: Mesh {
                dx: vec![plant.pitch / split as f64; cols * split],
                dy: vec![plant.pitch / split as f64; rows * split],
                dz: mesh.windows(2).map(|w| w[1] - w[0]).collect(),
            },
            files,
            tables,
            partial,
            reflector,
            materials,
            node_material,
            split,
            empty: placements
                .iter()
                .map(|row| row.iter().map(Option::is_none).collect())
                .collect(),
            ids,
            boundaries: case.boundaries,
            outline: case.outline,
            partial_rods: case.partial_rods,
            inputs: Box::new(LayoutInputs {
                case: case.clone(),
                cycle: cycle.clone(),
                plant: plant.clone(),
                inventory: inventory.clone(),
            }),
        })
    }

    /// `tables` の順（テーブルの使い方ごと）の、部分ノードごとの少数群定数から、NEM の問題を作る。反射体の層が
    /// 核分裂すればエラー。部分挿入のノードの物質は、最初の区間の定数を置いておく（`partial_nodes` の定数で
    /// 置き換えて解く。`rods::solve_rodded`）
    pub fn problem(&self, constants: &[Vec<FewGroup>]) -> Result<Problem, InputError> {
        if constants.len() != self.tables.len() {
            return invalid(format!(
                "{} few-group constants for {} table uses",
                constants.len(),
                self.tables.len()
            ));
        }
        let sub = self.split * self.split;
        for (subs, t) in constants.iter().zip(&self.tables) {
            let path = &self.files[t.file];
            if subs.len() != sub {
                return invalid(format!(
                    "{} has {} sub-nodes per state; radial-nodes = {} needs {sub} (ADR-0108)",
                    path.display(),
                    subs.len(),
                    self.split
                ));
            }
            if self.reflector[t.file]
                && subs
                    .iter()
                    .any(|c| c.nu_fission.iter().any(|&nf| nf != 0.0))
            {
                return invalid(format!(
                    "reflector layer {} has fission; load fuel with an assembly ID",
                    path.display()
                ));
            }
        }
        let mut materials = self
            .materials
            .iter()
            .map(|m| self.constants_of(constants, m.table, m.turns, m.sub, m.marshak))
            .collect::<Result<Vec<_>, _>>()?;
        for p in &self.partial {
            let &(first, _) = p
                .segments
                .first()
                .ok_or_else(|| InputError::Invalid("a partial node has no intervals".into()))?;
            materials.push(self.constants_of(constants, first, p.turns, p.sub, p.marshak)?);
        }
        Ok(Problem {
            mesh: self.mesh.clone(),
            materials,
            reflectors: (0..self.materials.len())
                .filter(|&m| self.materials[m].reflector)
                .collect(),
            node_material: self.node_material.clone(),
            boundaries: self.boundaries,
            outline: self.outline,
        })
    }

    /// 部分挿入のノードと、その区間の定数（`problem` と同じ `constants`）。物質の番号は `problem` の並び
    pub fn partial_nodes(
        &self,
        constants: &[Vec<FewGroup>],
    ) -> Result<Vec<PartialNode>, InputError> {
        self.partial
            .iter()
            .enumerate()
            .map(|(i, p)| {
                Ok(PartialNode {
                    site: p.site,
                    material: self.materials.len() + i,
                    segments: p
                        .segments
                        .iter()
                        .map(|&(table, length)| {
                            Ok(Segment {
                                constants: self
                                    .constants_of(constants, table, p.turns, p.sub, p.marshak)?,
                                length,
                            })
                        })
                        .collect::<Result<_, InputError>>()?,
                })
            })
            .collect()
    }

    /// テーブルの使い方 `table` の部分ノード `sub` の定数を、`turns` 回して、Marshak の層なら径方向の辺に
    /// r = 1/2 を置いたもの
    fn constants_of(
        &self,
        constants: &[Vec<FewGroup>],
        table: usize,
        turns: u8,
        sub: usize,
        marshak: bool,
    ) -> Result<FewGroup, InputError> {
        let file = || self.files[self.tables[table].file].display();
        let source = constants
            .get(table)
            .and_then(|subs| subs.get(sub))
            .ok_or_else(|| InputError::Invalid(format!("{}: no constants", file())))?;
        let mut c =
            rotated(source, turns).map_err(|e| InputError::Invalid(format!("{}: {e}", file())))?;
        if marshak {
            // Marshak の条件は r = 1/2 の真空の条件と同じ式（nem の `face_inputs`）なので、径方向の辺に
            // r = 1/2 を置く。格子計算の r を持つ辺があれば、どちらを使うか分からないのでエラー
            if c.vacuum_ratio.iter().take(4).any(Option::is_some) {
                return invalid(format!(
                    "{}: the layer is marked vacuum = \"marshak\", but its constants have a \
                     vacuum ratio J/phi on a radial side; drop the mark to use the ratio",
                    file()
                ));
            }
            let groups = c.groups();
            for r in c.vacuum_ratio.iter_mut().take(4) {
                *r = Some(vec![0.5; groups]);
            }
        }
        Ok(c)
    }

    /// ケースの軸方向のノードの境界 (cm。下から)
    pub fn axial_mesh(&self) -> &[f64] {
        &self.inputs.case.axial_mesh
    }

    /// 直方体の位置 `site` のノードのテーブルのファイル（`files` の番号）。部分挿入のノードは区間のファイル（区間は同じ
    /// 層のテーブルの制御棒の値なので、どれも同じファイル）。ノードのない位置は `None`
    pub fn site_file(&self, site: usize) -> Option<usize> {
        let m = self.node_material.get(site).copied().flatten()?;
        let table = match self.materials.get(m) {
            Some(material) => material.table,
            None => {
                self.partial
                    .get(m - self.materials.len())?
                    .segments
                    .first()?
                    .0
            }
        };
        Some(self.tables[table].file)
    }

    /// z の面ごとの、核分裂率の割合（下から。全体の和が 1）。メッシュの違う計算を、面を足し合わせて比べるのに使う
    pub fn axial_fission_fractions(&self, problem: &Problem, solution: &Solution) -> Vec<f64> {
        let mut rates = vec![0.0; problem.mesh.dz.len()];
        for (n, rate) in fission_rates(problem, solution).into_iter().enumerate() {
            rates[problem.mesh.position(n)[2]] += rate;
        }
        let total: f64 = rates.iter().sum();
        rates.iter().map(|r| r / total).collect()
    }

    /// 集合体の位置ごとの、軸方向に積分した核分裂率（集合体の位置の平均で割る）。上の行（+y）から、反射体の
    /// 位置とノードのない位置（`empty`）は `None`。`problem` はこの並びから作った問題、`solution` はその解
    pub fn assembly_fission_rates(
        &self,
        problem: &Problem,
        solution: &Solution,
    ) -> Vec<Vec<Option<f64>>> {
        let (nx, ny) = (problem.mesh.dx.len(), problem.mesh.dy.len());
        // 集合体の位置の数（部分ノードを足し合わせる）
        let (cols, rows) = (nx / self.split, ny / self.split);
        let mut rates = vec![vec![0.0; cols]; rows];
        for (n, rate) in fission_rates(problem, solution).into_iter().enumerate() {
            let [x, y, _] = problem.mesh.position(n);
            // 図の上の行が +y
            rates[rows - 1 - y / self.split][x / self.split] += rate;
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
/// 回す。辺の並びは（左、右、上、下、−z、+z）。90° 回すと、右の辺が上に、上が左に、左が下に、下が右に来て、
/// z の辺は変わらない（ADR-0128）。配列の長さは `FewGroup::validate` で先に確かめる。
pub fn rotated(constants: &FewGroup, turns: u8) -> Result<FewGroup, String> {
    fn turn<T: Clone>(sides: &[T]) -> Vec<T> {
        // 新しい（左、右、上、下、−z、+z）は、元の（上、下、右、左、−z、+z）
        [2, 3, 1, 0, 4, 5]
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
        normalized_corners: [3, 0, 1, 2]
            .iter()
            .filter_map(|&c| pins.normalized_corners.get(c).copied())
            .collect(),
        // 辺（左、右、上、下）の量: 新しい辺は元の（上、下、右、左）
        halves: [2, 3, 1, 0]
            .iter()
            .filter_map(|&s| pins.halves.get(s).copied())
            .collect(),
    }
}

#[cfg(test)]
mod tests;
