//! 少数群定数の分岐のテーブル（`lattice` から `core` に渡すファイル。ADR-0097、ADR-0200）。
//!
//! フィールド名付きの map 形式の MessagePack（ADR-0011、ADR-0013）:
//! `format` ("xlc-table"), `bounds`（群境界 (eV)、降順）, `axes`（軸ごとに `name`（`burnup`、`boron`、`moderator`、`fuel`、`rod`）。
//! 連続の軸は `points`（昇順。燃焼度は MWd/kg（ADR-0187）、ほかは ppm か K）を持ち、`moderator` は格子点ごとの `density` (g/cm³)、`boron` は `b10`（B-10 の
//! 原子分率）も持つ。`boron`、`moderator`、`fuel` は、内挿の変数 `transform`（`linear` か `sqrt`。ADR-0097）と、公称の値 `nominal`（格子点の
//! 1 つ。部分格子の組み立ての基準。ADR-0200。3 つとも持つか、どれも持たない）を持つ。離散の軸の `rod` は `points` の代わりに `values`
//! （制御棒の値の名前の列。`out` は棒のない状態）を持つ）,
//! `split`（集合体を各方向に分けたノードの数。ADR-0108）, `heavy_metal`（初めの重金属の質量 (g/cm、高さ 1 cm あたり)。燃焼の履歴から
//! 作ったテーブルだけが書き、ほかはキーを書かない。ADR-0190）,
//! `points`（解いた状態の点ごとに、軸の順の格子点の番号（離散の軸は値の番号）。辞書式の昇順。ADR-0200）,
//! `states`（`points` の順で、点ごとに split² の部分ノードの定数が続く（部分ノードは -y の行から、各行は -x から）。各状態は少数群定数の
//! `bounds` 以外の項目。`height` は数か nil、`pins` は `dx`、`dy`、`form`、`corners`、`normalized_corners`、`halves` を持つ map か nil で、
//! どちらもキーは省けない）。軸のないテーブルは点が1つ（番号の列は空）。
//!
//! 点の集合は、制御棒の値ごとに、`boron`、`moderator`、`fuel` の軸の直積か（公称の値がないとき）、公称の状態と、公称と違う軸の組
//! （成分）ごとの部分格子の和（cut-HDMR。成分の集合は部分集合について閉じる。ADR-0200）。点の成分は、公称と違う軸の組として点の
//! 番号から決まる。燃焼度の軸を持つテーブルは、点ごとに解いた燃焼度の点が違ってよく、すべての燃焼度の点で解いていない点があれば、
//! 公称の状態（制御棒の軸があれば `out`）がすべての燃焼度の点を持つ。

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use crate::fewgroup::{FewGroup, FewGroupError, Pins};
use crate::msgpack::{self, MsgpackError};

/// 連続の軸の内挿の変数（ADR-0097）。core は、この変数について格子点の間を線形に内挿する
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Transform {
    /// 値そのもの
    Linear,
    /// 値の平方根（燃料の温度。共鳴の Doppler 拡がりの効果が √T にほぼ線形）
    Sqrt,
}

impl Transform {
    /// 内挿の変数の値
    pub fn apply(self, x: f64) -> f64 {
        match self {
            Self::Linear => x,
            Self::Sqrt => x.sqrt(),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Linear => "linear",
            Self::Sqrt => "sqrt",
        }
    }
}

/// テーブルの軸。点は昇順。
#[derive(Debug, Clone, PartialEq)]
pub enum Axis {
    /// 燃焼度 (MWd/kg、初めの重金属)。燃焼の履歴の燃焼点から作る（ADR-0187）。内挿は ADR-0192 の 3 次 Hermite
    Burnup { points: Vec<f64> },
    /// ホウ素濃度 (ppm) と、ホウ素の中の B-10 の原子分率
    Boron {
        ppm: Vec<f64>,
        b10: f64,
        nominal: Option<f64>,
        transform: Transform,
    },
    /// 減速材の温度 (K) と、格子点ごとの密度 (g/cm³)
    Moderator {
        temperature: Vec<f64>,
        density: Vec<f64>,
        nominal: Option<f64>,
        transform: Transform,
    },
    /// 燃料の温度 (K)
    Fuel {
        temperature: Vec<f64>,
        nominal: Option<f64>,
        transform: Transform,
    },
    /// 制御棒（離散の軸。ADR-0097、ADR-0127）: 値の名前の列。棒のない状態は `out`、ほかは棒の区間の材料の名前
    /// （炉心の形の制御棒の型が使う名前）。core は内挿せず、ノードの区間の値で選ぶ
    Rod { values: Vec<String> },
}

/// 制御棒の軸の、棒のない状態の値の名前
pub const ROD_OUT: &str = "out";

impl Axis {
    /// ホウ素の軸。内挿の変数は ppm（ADR-0097）
    pub fn boron(ppm: Vec<f64>, b10: f64) -> Self {
        Self::Boron {
            ppm,
            b10,
            nominal: None,
            transform: Transform::Linear,
        }
    }

    /// 減速材の軸。内挿の変数は温度（ADR-0097）
    pub fn moderator(temperature: Vec<f64>, density: Vec<f64>) -> Self {
        Self::Moderator {
            temperature,
            density,
            nominal: None,
            transform: Transform::Linear,
        }
    }

    /// 燃料の温度の軸。内挿の変数は √T（ADR-0097）
    pub fn fuel(temperature: Vec<f64>) -> Self {
        Self::Fuel {
            temperature,
            nominal: None,
            transform: Transform::Sqrt,
        }
    }

    /// 公称の値を除いた軸
    pub fn without_nominal(mut self) -> Self {
        match &mut self {
            Self::Boron { nominal, .. }
            | Self::Moderator { nominal, .. }
            | Self::Fuel { nominal, .. } => *nominal = None,
            Self::Burnup { .. } | Self::Rod { .. } => {}
        }
        self
    }

    /// 公称の値を `value` にした軸（`boron`、`moderator`、`fuel` だけ。ほかの軸はそのまま）
    pub fn with_nominal(mut self, value: f64) -> Self {
        match &mut self {
            Self::Boron { nominal, .. }
            | Self::Moderator { nominal, .. }
            | Self::Fuel { nominal, .. } => *nominal = Some(value),
            Self::Burnup { .. } | Self::Rod { .. } => {}
        }
        self
    }

    pub fn name(&self) -> &'static str {
        match self {
            Self::Burnup { .. } => "burnup",
            Self::Boron { .. } => "boron",
            Self::Moderator { .. } => "moderator",
            Self::Fuel { .. } => "fuel",
            Self::Rod { .. } => "rod",
        }
    }

    /// 連続の軸の格子点。離散の軸（`Rod`）は空で、値の数は `len`
    pub fn points(&self) -> &[f64] {
        match self {
            Self::Burnup { points } => points,
            Self::Boron { ppm, .. } => ppm,
            Self::Moderator { temperature, .. } | Self::Fuel { temperature, .. } => temperature,
            Self::Rod { .. } => &[],
        }
    }

    /// 部分格子で組み立てる軸（ホウ素、減速材、燃料）の、公称の値と内挿の変数。ほかの軸は `None`
    pub fn state_axis(&self) -> Option<(Option<f64>, Transform)> {
        match self {
            Self::Boron {
                nominal, transform, ..
            }
            | Self::Moderator {
                nominal, transform, ..
            }
            | Self::Fuel {
                nominal, transform, ..
            } => Some((*nominal, *transform)),
            Self::Burnup { .. } | Self::Rod { .. } => None,
        }
    }

    /// 公称の値の格子点の番号（公称の値のない軸と、部分格子で組み立てない軸は `None`）
    pub fn nominal_index(&self) -> Option<usize> {
        let (nominal, _) = self.state_axis()?;
        let nominal = nominal?;
        self.points().iter().position(|&p| p == nominal)
    }

    /// 軸の格子点（離散の軸は値）の数
    pub fn len(&self) -> usize {
        match self {
            Self::Rod { values } => values.len(),
            _ => self.points().len(),
        }
    }

    /// 格子点も値もないか（`validate_axes` はエラーにする）
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// 軸の並びを確かめる: 同じ軸が2つなく、点が範囲の中で増え、減速材の密度が点ごとに正で、B-10 の割合が 0〜1。公称の値は
/// 格子点の 1 つで、部分格子で組み立てる軸（ホウ素、減速材、燃料）のすべてが持つか、どれも持たない。
pub fn validate_axes(axes: &[Axis]) -> Result<(), String> {
    let nominals: Vec<bool> = axes
        .iter()
        .filter_map(|a| a.state_axis().map(|(n, _)| n.is_some()))
        .collect();
    if nominals.iter().any(|&n| n) && !nominals.iter().all(|&n| n) {
        return Err(
            "give the nominal value to all of the boron, moderator and fuel axes or to none".into(),
        );
    }
    let mut names = Vec::new();
    for axis in axes {
        let name = axis.name();
        if names.contains(&name) {
            return Err(format!("axis {name} appears twice"));
        }
        names.push(name);
        if let Axis::Rod { values } = axis {
            // 値は名前で選ぶので、空でなく、区切りの空白を含まず、重ならない
            let distinct = values
                .iter()
                .enumerate()
                .all(|(i, v)| !values[..i].contains(v));
            if values.is_empty()
                || !distinct
                || values
                    .iter()
                    .any(|v| v.is_empty() || v.chars().any(char::is_whitespace))
            {
                return Err(format!(
                    "axis rod: values must be distinct names without spaces, got {values:?}"
                ));
            }
            continue;
        }
        let points = axis.points();
        // 燃焼度とホウ素濃度は 0 を含み、温度は正
        let valid = |p: f64| match axis {
            Axis::Burnup { .. } | Axis::Boron { .. } => p >= 0.0 && p.is_finite(),
            Axis::Moderator { .. } | Axis::Fuel { .. } => p > 0.0 && p.is_finite(),
            // 制御棒の軸は点を持たない（上で値の名前を確かめて次の軸に進む）
            Axis::Rod { .. } => false,
        };
        if points.is_empty()
            || !points.iter().all(|&p| valid(p))
            || !points.windows(2).all(|w| w[0] < w[1])
        {
            return Err(format!(
                "axis {name}: points must be valid and increasing, got {points:?}"
            ));
        }
        if let Some((Some(nominal), _)) = axis.state_axis()
            && !points.contains(&nominal)
        {
            return Err(format!(
                "axis {name}: the nominal value {nominal} is not a point of {points:?}"
            ));
        }
        match axis {
            Axis::Boron { b10, .. } if !(0.0..=1.0).contains(b10) => {
                return Err(format!("axis boron: invalid B-10 atom fraction {b10}"));
            }
            Axis::Moderator {
                temperature,
                density,
                ..
            } if density.len() != temperature.len()
                || !density.iter().all(|&d| d > 0.0 && d.is_finite()) =>
            {
                return Err(format!(
                    "axis moderator: {} positive densities are needed, got {density:?}",
                    temperature.len()
                ));
            }
            _ => {}
        }
    }
    Ok(())
}

/// 軸の直積のすべての点（最初の軸がいちばん遅く変わる。辞書式の昇順）
pub fn product_points(axes: &[Axis]) -> Vec<Vec<usize>> {
    let counts: Vec<usize> = axes.iter().map(Axis::len).collect();
    let total: usize = counts.iter().product();
    (0..total)
        .map(|index| {
            let mut at = vec![0; counts.len()];
            let mut rest = index;
            for (k, &n) in counts.iter().enumerate().rev() {
                at[k] = rest % n;
                rest /= n;
            }
            at
        })
        .collect()
}

/// 点 `at`（部分格子で組み立てる軸の番号の列）の、公称 `nominal` と違う軸の組（`at` の中の位置の bit）
pub fn component(at: &[usize], nominal: &[usize]) -> u32 {
    at.iter()
        .zip(nominal)
        .enumerate()
        .filter(|(_, (a, n))| a != n)
        .fold(0, |mask, (i, _)| mask | (1 << i))
}

/// 点の番号が軸の範囲の中で辞書式に増え、点の集合が制御棒の値ごとに直積か、公称の状態と成分の部分格子の和（成分の集合は部分
/// 集合について閉じ、各成分の部分格子の点をすべて持つ）で、すべての燃焼度の点で解いていない点があれば公称の状態がすべての燃焼度の
/// 点を持つこと（モジュールの説明）
fn validate_points(axes: &[Axis], points: &[Vec<usize>]) -> Result<(), String> {
    for p in points {
        if p.len() != axes.len() || p.iter().zip(axes).any(|(&i, a)| i >= a.len()) {
            return Err(format!(
                "point {p:?} does not index the axes {:?}",
                axes.iter().map(Axis::name).collect::<Vec<_>>()
            ));
        }
    }
    if let Some(w) = points.windows(2).find(|w| w[0] >= w[1]) {
        return Err(format!(
            "points must increase in lexicographic order, got {:?} before {:?}",
            w[0], w[1]
        ));
    }
    let burnup = axes.iter().position(|a| matches!(a, Axis::Burnup { .. }));
    let rod = axes.iter().position(|a| matches!(a, Axis::Rod { .. }));
    let state_axes: Vec<usize> = (0..axes.len())
        .filter(|&a| axes[a].state_axis().is_some())
        .collect();
    let nominal: Option<Vec<usize>> = state_axes
        .iter()
        .map(|&a| axes[a].nominal_index())
        .collect();
    // 燃焼度の軸を除いた点ごとに、解いた燃焼度の点
    let mut burnups: BTreeMap<Vec<usize>, Vec<usize>> = BTreeMap::new();
    for p in points {
        let key: Vec<usize> = (0..axes.len())
            .filter(|&a| Some(a) != burnup)
            .map(|a| p[a])
            .collect();
        burnups
            .entry(key)
            .or_default()
            .push(burnup.map_or(0, |b| p[b]));
    }
    // 燃焼度の軸を除いた番号の列の中の、軸 `a` の位置
    let in_key = |a: usize| a - usize::from(burnup.is_some_and(|b| b < a));
    let rod_key = rod.map(in_key);
    let rods = rod.map_or(1, |r| axes[r].len());
    for r in 0..rods {
        let at: BTreeSet<Vec<usize>> = burnups
            .keys()
            .filter(|key| rod_key.is_none_or(|k| key[k] == r))
            .map(|key| state_axes.iter().map(|&a| key[in_key(a)]).collect())
            .collect();
        let value = rod.map_or_else(String::new, |ra| match &axes[ra] {
            Axis::Rod { values } => format!(" of rod value {}", values[r]),
            _ => String::new(),
        });
        if at.is_empty() {
            return Err(format!("no points{value}"));
        }
        let expected = match &nominal {
            None => state_axes.iter().map(|&a| axes[a].len()).product(),
            Some(nominal) => {
                let masks: BTreeSet<u32> = at.iter().map(|p| component(p, nominal)).collect();
                if let Some(m) = masks.iter().find(|&&m| {
                    (0..state_axes.len())
                        .any(|b| m & (1 << b) != 0 && !masks.contains(&(m & !(1 << b))))
                }) {
                    return Err(format!(
                        "the points{value} hold the subgrid of the axes {:?} but not all of its subsets",
                        names_of(axes, &state_axes, *m)
                    ));
                }
                masks
                    .iter()
                    .map(|&m| {
                        (0..state_axes.len())
                            .filter(|b| m & (1 << b) != 0)
                            .map(|b| axes[state_axes[b]].len() - 1)
                            .product::<usize>()
                    })
                    .sum()
            }
        };
        if at.len() != expected {
            return Err(format!(
                "the {} points{value} are not {}",
                at.len(),
                if nominal.is_some() {
                    "the full subgrids of their components"
                } else {
                    "the product of the axes (tables without nominal values)"
                }
            ));
        }
    }
    if let Some(b) = burnup {
        let all = axes[b].len();
        // core は公称の状態との差を外挿しないので、どの点も燃焼度の軸の最初と最後の点で解く
        if burnups
            .values()
            .any(|solved| solved.first() != Some(&0) || solved.last() != Some(&(all - 1)))
        {
            return Err(
                "every point must be solved at the first and the last burnup (the differences from the nominal \
                 state are not extrapolated. ADR-0200)"
                    .into(),
            );
        }
        if burnups.values().any(|solved| solved.len() != all) {
            let nominal = nominal.ok_or(
                "points solved at only some burnups need the nominal values (the burnups of the other points are \
                 interpolated as differences from the nominal state. ADR-0200)",
            )?;
            let out = match rod.map(|ra| &axes[ra]) {
                Some(Axis::Rod { values }) => Some(
                    values
                        .iter()
                        .position(|v| v == ROD_OUT)
                        .ok_or("points solved at only some burnups need the rod value out")?,
                ),
                _ => None,
            };
            let mut key = vec![0; axes.len() - 1];
            for (&a, &n) in state_axes.iter().zip(&nominal) {
                key[in_key(a)] = n;
            }
            if let (Some(k), Some(out)) = (rod_key, out) {
                key[k] = out;
            }
            if burnups.get(&key).is_none_or(|solved| solved.len() != all) {
                return Err(
                    "points solved at only some burnups need the nominal state at every burnup"
                        .into(),
                );
            }
        }
    }
    Ok(())
}

/// 成分 `mask` の軸の名前
fn names_of(axes: &[Axis], state_axes: &[usize], mask: u32) -> Vec<&'static str> {
    state_axes
        .iter()
        .enumerate()
        .filter(|(b, _)| mask & (1 << b) != 0)
        .map(|(_, &a)| axes[a].name())
        .collect()
}

/// 分岐のテーブル。
#[derive(Debug, Clone, PartialEq)]
pub struct Table {
    pub axes: Vec<Axis>,
    /// 集合体を各方向に分けたノードの数（ADR-0108）。分けなければ 1
    pub split: usize,
    /// [点] の、軸の順の格子点の番号（離散の軸は値の番号）。辞書式の昇順（モジュールの説明の点の集合）
    pub points: Vec<Vec<usize>>,
    /// `points` の点ごとに、split × split の部分ノードの定数（部分ノードがいちばん速く変わる。部分ノードは -y の行から、各行は
    /// -x から）
    pub states: Vec<FewGroup>,
    /// 初めの重金属（Z ≥ 90）の質量 (g/cm、格子全体の高さ 1 cm あたり)。燃焼の履歴から作ったテーブルだけが持ち、core が
    /// ノードの燃焼度を進めるのに使う（ADR-0190）
    pub heavy_metal: Option<f64>,
}

/// テーブルのファイルの読み書きのエラー。
#[derive(Debug, thiserror::Error)]
pub enum TableError {
    #[error(transparent)]
    Msgpack(#[from] MsgpackError),
    #[error("expected format xlc-table, found {0}")]
    Format(String),
    #[error("state {state}")]
    State {
        state: usize,
        #[source]
        source: FewGroupError,
    },
    #[error("{0}")]
    Invalid(String),
}

const FORMAT: &str = "xlc-table";

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct TableFile {
    format: String,
    bounds: Vec<f64>,
    axes: Vec<AxisFile>,
    split: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    heavy_metal: Option<f64>,
    points: Vec<Vec<usize>>,
    states: Vec<StateFile>,
}

/// ファイルの軸（テーブルと媒質で共有する。ADR-0128）
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct AxisFile {
    name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    points: Option<Vec<f64>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    values: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    density: Option<Vec<f64>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    b10: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    nominal: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    transform: Option<String>,
}

impl From<&Axis> for AxisFile {
    fn from(axis: &Axis) -> Self {
        Self {
            name: axis.name().into(),
            points: match axis {
                Axis::Rod { .. } => None,
                _ => Some(axis.points().to_vec()),
            },
            values: match axis {
                Axis::Rod { values } => Some(values.clone()),
                _ => None,
            },
            density: match axis {
                Axis::Moderator { density, .. } => Some(density.clone()),
                _ => None,
            },
            b10: match axis {
                Axis::Boron { b10, .. } => Some(*b10),
                _ => None,
            },
            nominal: axis.state_axis().and_then(|(n, _)| n),
            transform: axis.state_axis().map(|(_, t)| t.name().into()),
        }
    }
}

impl TryFrom<AxisFile> for Axis {
    type Error = String;

    fn try_from(a: AxisFile) -> Result<Self, String> {
        let transform = match a.transform.as_deref() {
            None => None,
            Some("linear") => Some(Transform::Linear),
            Some("sqrt") => Some(Transform::Sqrt),
            Some(other) => {
                return Err(format!(
                    "axis {}: unknown transform {other} (linear or sqrt)",
                    a.name
                ));
            }
        };
        let nominal = a.nominal;
        match (
            a.name.as_str(),
            a.points,
            a.values,
            a.density,
            a.b10,
            transform,
        ) {
            ("burnup", Some(points), None, None, None, None) if nominal.is_none() => {
                Ok(Axis::Burnup { points })
            }
            ("boron", Some(ppm), None, None, Some(b10), Some(transform)) => Ok(Axis::Boron {
                ppm,
                b10,
                nominal,
                transform,
            }),
            ("moderator", Some(temperature), None, Some(density), None, Some(transform)) => {
                Ok(Axis::Moderator {
                    temperature,
                    density,
                    nominal,
                    transform,
                })
            }
            ("fuel", Some(temperature), None, None, None, Some(transform)) => Ok(Axis::Fuel {
                temperature,
                nominal,
                transform,
            }),
            ("rod", None, Some(values), None, None, None) if nominal.is_none() => {
                Ok(Axis::Rod { values })
            }
            (name, ..) => Err(format!(
                "axis {name}: unknown name or wrong fields (burnup, boron, moderator and fuel have points, \
                 boron has b10, moderator has density, boron, moderator and fuel have a transform and may have a \
                 nominal value, rod has values instead of points)"
            )),
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StateFile {
    k_inf: f64,
    total: Vec<f64>,
    absorption: Vec<f64>,
    nu_fission: Vec<f64>,
    fission: Vec<f64>,
    chi: Vec<f64>,
    diffusion: Vec<f64>,
    scatter: Vec<Vec<f64>>,
    adf: Vec<Vec<f64>>,
    vacuum_ratio: Vec<Option<Vec<f64>>>,
    normalized_adf: Vec<bool>,
    /// nil を書いても、キーは省けない（ADR-0129）
    #[serde(deserialize_with = "crate::msgpack::nil_or")]
    height: Option<f64>,
    /// nil を書いても、キーは省けない（`deserialize_with` を付けると、serde は省いた `Option` を `None` にしない）
    #[serde(deserialize_with = "crate::msgpack::nil_or")]
    pins: Option<PinsFile>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PinsFile {
    dx: Vec<f64>,
    dy: Vec<f64>,
    form: Vec<Vec<Vec<f64>>>,
    corners: Vec<Vec<f64>>,
    normalized_corners: Vec<bool>,
    halves: Vec<bool>,
}

impl Table {
    /// 軸のない、状態が1つのテーブル。
    pub fn single(state: FewGroup) -> Self {
        Self {
            axes: Vec::new(),
            split: 1,
            points: vec![Vec::new()],
            states: vec![state],
            heavy_metal: None,
        }
    }

    /// 点 `at`（軸の順の番号）の `points` の番号
    pub fn point_index(&self, at: &[usize]) -> Option<usize> {
        self.points.binary_search_by(|p| p.as_slice().cmp(at)).ok()
    }

    /// 状態ごとの部分ノードの数（split²）
    pub fn sub_nodes(&self) -> usize {
        self.split * self.split
    }

    /// 軸、点の集合、状態の数、状態の値、状態の間で共通であるべき量（群境界、f* の印、真空の辺、z の厚さ）を確かめる。
    pub fn validate(&self) -> Result<(), TableError> {
        let invalid = |message: String| Err(TableError::Invalid(message));
        validate_axes(&self.axes).map_err(TableError::Invalid)?;
        if self.split == 0 {
            return invalid("split must be positive".into());
        }
        if let Some(m) = self.heavy_metal
            && !(m > 0.0 && m.is_finite())
        {
            return invalid(format!("heavy_metal must be positive, got {m}"));
        }
        validate_points(&self.axes, &self.points).map_err(TableError::Invalid)?;
        let points = self.points.len();
        let sub = self.sub_nodes();
        if self.states.len() != points * sub {
            return invalid(format!(
                "{} states for {points} points and {sub} sub-nodes per point",
                self.states.len()
            ));
        }
        if self.states.is_empty() {
            return invalid("no states".into());
        }
        for (i, state) in self.states.iter().enumerate() {
            // 状態の間で比べる量は、状態 0 の同じ部分ノードと比べる
            let first = &self.states[i % sub];
            state
                .validate()
                .map_err(|source| TableError::State { state: i, source })?;
            let vacuum = |s: &FewGroup| {
                s.vacuum_ratio
                    .iter()
                    .map(Option::is_some)
                    .collect::<Vec<_>>()
            };
            // 群境界は表全体の量（ファイルには1つだけ書く）なので、状態 0 の最初の定数と比べる
            if state.bounds != self.states[0].bounds
                || state.normalized_adf != first.normalized_adf
                || vacuum(state) != vacuum(first)
                || state.height != first.height
            {
                return invalid(format!(
                    "state {i}: group bounds, f* sides, vacuum sides and the height must equal those of the same sub-node of state 0"
                ));
            }
            // 形状関数は状態ごとに内挿するので、どの状態も持つか持たないかで、セルの境界も同じ
            let cells = |s: &FewGroup| {
                s.pins.as_ref().map(|p| {
                    (
                        p.dx.clone(),
                        p.dy.clone(),
                        p.normalized_corners.clone(),
                        p.halves.clone(),
                    )
                })
            };
            if cells(state) != cells(first) {
                return invalid(format!(
                    "state {i}: pins (form functions) and their cells must equal those of the same sub-node of state 0"
                ));
            }
        }
        Ok(())
    }

    /// MessagePack（フィールド名付きの map 形式）で書く。書く前に `validate` で確かめる。
    pub fn to_msgpack(&self) -> Result<Vec<u8>, TableError> {
        self.validate()?;
        let bounds = self
            .states
            .first()
            .map(|s| s.bounds.clone())
            .unwrap_or_default();
        let axes = self.axes.iter().map(AxisFile::from).collect();
        let states = self
            .states
            .iter()
            .map(|s| StateFile {
                k_inf: s.k_inf,
                total: s.total.clone(),
                absorption: s.absorption.clone(),
                nu_fission: s.nu_fission.clone(),
                fission: s.fission.clone(),
                chi: s.chi.clone(),
                diffusion: s.diffusion.clone(),
                scatter: s.scatter.clone(),
                adf: s.adf.clone(),
                vacuum_ratio: s.vacuum_ratio.clone(),
                normalized_adf: s.normalized_adf.clone(),
                height: s.height,
                pins: s.pins.as_ref().map(|p| PinsFile {
                    dx: p.dx.clone(),
                    dy: p.dy.clone(),
                    form: p.form.clone(),
                    corners: p.corners.clone(),
                    normalized_corners: p.normalized_corners.clone(),
                    halves: p.halves.clone(),
                }),
            })
            .collect();
        Ok(msgpack::to_named(&TableFile {
            format: FORMAT.into(),
            bounds,
            axes,
            split: self.split,
            heavy_metal: self.heavy_metal,
            points: self.points.clone(),
            states,
        })?)
    }

    /// MessagePack から読み、`validate` で値を確かめる。
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self, TableError> {
        let file: TableFile = msgpack::from_value(msgpack::read_value(bytes)?)?;
        if file.format != FORMAT {
            return Err(TableError::Format(file.format));
        }
        let axes = file
            .axes
            .into_iter()
            .map(Axis::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map_err(TableError::Invalid)?;
        let states = file
            .states
            .into_iter()
            .map(|s| FewGroup {
                bounds: file.bounds.clone(),
                k_inf: s.k_inf,
                total: s.total,
                absorption: s.absorption,
                nu_fission: s.nu_fission,
                fission: s.fission,
                chi: s.chi,
                diffusion: s.diffusion,
                scatter: s.scatter,
                adf: s.adf,
                vacuum_ratio: s.vacuum_ratio,
                normalized_adf: s.normalized_adf,
                height: s.height,
                pins: s.pins.map(|p| Pins {
                    dx: p.dx,
                    dy: p.dy,
                    form: p.form,
                    corners: p.corners,
                    normalized_corners: p.normalized_corners,
                    halves: p.halves,
                }),
            })
            .collect();
        let table = Self {
            axes,
            split: file.split,
            points: file.points,
            states,
            heavy_metal: file.heavy_metal,
        };
        table.validate()?;
        Ok(table)
    }
}

#[cfg(test)]
mod tests;
