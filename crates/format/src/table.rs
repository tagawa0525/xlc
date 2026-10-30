//! 少数群定数の分岐のテーブル（`lattice` から `core` に渡すファイル。ADR-0097、ADR-0200）。
//!
//! フィールド名付きの map 形式の MessagePack（ADR-0011、ADR-0013）:
//! `format` ("xlc-table"), `bounds`（群境界 (eV)、降順）, `axes`（軸ごとに `name`（`burnup`、`boron`、`moderator`、`fuel`、`rod`）。
//! 連続の軸は `points`（昇順。燃焼度は MWd/kg（ADR-0187）、ほかは ppm か K）を持ち、`moderator` は格子点ごとの `density` (g/cm³)、`boron` は `b10`（B-10 の
//! 原子分率）も持つ。`boron`、`moderator`、`fuel` は、内挿の変数 `transform`（`linear` か `sqrt`。ADR-0097）と、公称の値 `nominal`（格子点の
//! 1 つ。部分格子の組み立ての基準。ADR-0200。3 つとも持つか、どれも持たない）を持つ。離散の軸の `rod` は `points` の代わりに `values`
//! （制御棒の値の名前の列。`out` は棒のない状態）を持つ）,
//! `split`（集合体を各方向に分けたノードの数。ADR-0108）, `heavy_metal`（初めの重金属の質量 (g/cm、高さ 1 cm あたり)。燃焼の履歴から
//! 作ったテーブルだけが書き、ほかはキーを書かない。ADR-0190）, `micro`（ミクロ燃焼の核種の `nuclides`、燃焼度の点ごとの均質化した
//! 数密度 `densities`、`points` の点ごとの核種ごとの微視的断面積 `states`（`total`、`nu_fission`、`fission`、`scatter`、`reactions`（`mt` と
//! 経路ごとの `channels`）、`thermal_fission`（nil か値。キーは省けない）、`derivative`（自分の数密度による格子全体の少数群定数の微分の `total`、
//! `nu_fission`、`fission`、`scatter`。数密度が 0 の点は σᵢ。微分を持つ核種だけが書き、ほかはキーを書かない（nil はエラー）。ADR-0221））。燃焼度の軸を持ち、分けない
//! テーブルだけが書き、ほかはキーを書かない。ADR-0216）,
//! `histories`（履歴の変数ごとの補正の残り R の列。変数の名前 `variable`（`moderator-density`）、公称の燃焼の値 `nominal`（H⁰）、履歴の
//! 分岐の値 `branch`（Hᵇ）、制御棒の値ごとの `rods`（制御棒の値の名前 `rod`、R を持つ燃焼度の点 `burnup`（燃焼度の軸の点で、最初と最後を
//! 含み、テーブルが履歴の分岐の条件の状態を解いた点）、点ごとの `residuals`（`total`、`absorption`、`nu_fission`、`fission`、`chi`、
//! `diffusion`、`scatter`、`adf`））。履歴の分岐の条件は、変数の軸（減速材の密度は `moderator`）だけが値 `branch` の格子点で、ほかの軸は
//! 公称の状態。ミクロ燃焼のデータを持つテーブルだけが書き、持たなければキーを書かない。ADR-0215）,
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
    /// ミクロ燃焼の核種の均質化した数密度と微視的断面積（ADR-0187、ADR-0216）。燃焼の履歴から作ったテーブルだけが持てる
    pub micro: Option<Micro>,
    /// 履歴の変数ごとの、ミクロ燃焼で追わない核種の分の補正の残り R（ADR-0215）。ミクロ燃焼のデータを持つテーブルだけが持てる
    pub histories: Vec<HistoryCorrection>,
}

/// ミクロ燃焼で追う核種の、燃焼度の点ごとの格子で均質化した数密度と、状態ごとの少数群の微視的断面積（ADR-0216）。格子全体の
/// 均質化なので、部分ノードに分けないテーブル（`split` が 1）だけが持つ
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Micro {
    /// 追う核種の名前
    pub nuclides: Vec<String>,
    /// [燃焼度の軸の点][核種] の均質化した数密度 (1/(barn cm))
    pub densities: Vec<Vec<f64>>,
    /// [`Table::points` の点][核種] の微視的断面積
    pub states: Vec<Vec<MicroNuclide>>,
}

/// 1 つの状態の、1 つの核種の少数群の微視的断面積 (barn)（ADR-0216）
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MicroNuclide {
    pub total: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub fission: Vec<f64>,
    /// P0 の散乱行列 [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
    /// ライブラリのチェーンの反応（ADR-0172）ごとの、経路ごとの断面積
    pub reactions: Vec<MicroReaction>,
    /// 核分裂率のうち、収率の境界（ADR-0183）より下の割合。核分裂しない核種は nil（キーは省けない）
    #[serde(deserialize_with = "crate::msgpack::nil_or")]
    pub thermal_fission: Option<f64>,
    /// 自分の数密度による格子全体の少数群定数の微分（ADR-0221）。微分を持たない核種はキーを書かない（nil はエラー）
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "crate::msgpack::present"
    )]
    pub derivative: Option<MicroDerivative>,
}

/// 1 つの状態の、1 つの核種の均質化した数密度 N̄ᵢ による格子全体の少数群定数の微分 ∂Σ/∂N̄ᵢ (barn)（ADR-0221）。すべての核種と
/// 燃焼させない物質の分を含む。ただし N̄ᵢ が 0 の点は、相対の摂動が解けないので lattice は σᵢ（自分の項だけ）を書く（#0206）
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MicroDerivative {
    pub total: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub fission: Vec<f64>,
    /// P0 の散乱行列 [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
}

/// チェーンの反応の代表の MT と、経路（ライブラリとチェーンの経路の並び）ごとの少数群の断面積 (barn)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MicroReaction {
    pub mt: u32,
    pub channels: Vec<Vec<f64>>,
}

/// 履歴の変数の 1 つの補正（ADR-0215）。core はノードの少数群定数に、w = (H − H⁰) / (Hᵇ − H⁰) を掛けた残り R を足す
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HistoryCorrection {
    /// 履歴の変数の名前（[`MODERATOR_DENSITY`]）
    pub variable: String,
    /// 公称の燃焼の履歴の変数の値 H⁰（減速材の密度は g/cm³）
    pub nominal: f64,
    /// 履歴の分岐の変数の値 Hᵇ
    pub branch: f64,
    /// 制御棒の値ごとの残り
    pub rods: Vec<HistoryRod>,
}

/// 1 つの制御棒の値の、燃焼度の点ごとの残り R（ADR-0215）
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HistoryRod {
    /// 制御棒の値の名前（制御棒の軸の値。軸のないテーブルは [`ROD_OUT`]）
    pub rod: String,
    /// R を持つ燃焼度の点 (MWd/kg)。燃焼度の軸の点で、最初と最後の点を含み、テーブルが履歴の分岐の条件の状態を解いた点
    pub burnup: Vec<f64>,
    /// `burnup` の点ごとの残り
    pub residuals: Vec<Residual>,
}

/// 残り R の少数群定数の項目（状態の少数群定数と同じ単位）。k∞、真空の辺の比、形状関数は持たない
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Residual {
    pub total: Vec<f64>,
    pub absorption: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub fission: Vec<f64>,
    pub chi: Vec<f64>,
    pub diffusion: Vec<f64>,
    /// [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
    /// [辺][群]
    pub adf: Vec<Vec<f64>>,
}

/// 減速材の密度の履歴の変数の名前（ADR-0215）
pub const MODERATOR_DENSITY: &str = "moderator-density";

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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    micro: Option<Micro>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    histories: Vec<HistoryCorrection>,
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
            micro: None,
            histories: Vec::new(),
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
        if let Some(micro) = &self.micro {
            self.validate_micro(micro).map_err(TableError::Invalid)?;
        }
        self.validate_histories().map_err(TableError::Invalid)?;
        Ok(())
    }

    /// 履歴の補正（ADR-0215）: ミクロ燃焼のデータを持つテーブルの、重ならない既知の変数で、H⁰ と Hᵇ が減速材の軸の公称と公称でない
    /// 格子点の密度。制御棒の値ごとの燃焼度の点は、その条件の状態（変数の軸だけが Hᵇ の格子点で、ほかの軸は公称）をテーブルが解いた
    /// 燃焼度の点で、増え、最初と最後の点を含み、点ごとに状態と同じ群と辺の数の有限の残りを持つ
    fn validate_histories(&self) -> Result<(), String> {
        if self.histories.is_empty() {
            return Ok(());
        }
        if self.micro.is_none() {
            return Err(
                "history residuals need micro data (they are left after the micro correction. ADR-0215)".into(),
            );
        }
        let position = |f: fn(&Axis) -> bool| self.axes.iter().position(f);
        let burnup = position(|a| matches!(a, Axis::Burnup { .. }))
            .ok_or("history residuals need a burnup axis")?;
        let moderator = position(|a| matches!(a, Axis::Moderator { .. }));
        let rod = position(|a| matches!(a, Axis::Rod { .. }));
        let burnups = self.axes[burnup].points();
        let groups = self.states[0].bounds.len() - 1;
        let sides = self.states[0].adf.len();
        for (h, history) in self.histories.iter().enumerate() {
            let variable = &history.variable;
            if self.histories[..h].iter().any(|o| &o.variable == variable) {
                return Err(format!("history variable {variable} appears twice"));
            }
            if variable != MODERATOR_DENSITY {
                return Err(format!(
                    "unknown history variable {variable:?} ({MODERATOR_DENSITY}. ADR-0215)"
                ));
            }
            // 減速材の密度: H⁰ は公称の格子点の密度、Hᵇ はほかの 1 つの格子点の密度で、その格子点が履歴の分岐の条件
            let (m, density) = match moderator.map(|m| (m, &self.axes[m])) {
                Some((m, Axis::Moderator { density, .. })) => (m, density),
                _ => {
                    return Err(format!(
                        "history variable {variable} needs a moderator axis"
                    ));
                }
            };
            let nominal = self.axes[m].nominal_index().ok_or_else(|| {
                format!(
                    "history variable {variable}: the branch condition is the nominal state but for the moderator, \
                     and the table has no nominal values"
                )
            })?;
            if history.nominal != density[nominal] {
                return Err(format!(
                    "history variable {variable}: nominal value {} is not the density {} of the nominal moderator point",
                    history.nominal, density[nominal]
                ));
            }
            // 変数の差 Hᵇ − H⁰ で core が重みを割る
            if history.branch == history.nominal {
                return Err(format!(
                    "history variable {variable}: branch value {} equals the nominal value",
                    history.branch
                ));
            }
            let condition = match density
                .iter()
                .enumerate()
                .filter(|&(i, &d)| i != nominal && d == history.branch)
                .map(|(i, _)| i)
                .collect::<Vec<_>>()[..]
            {
                [only] => only,
                _ => {
                    return Err(format!(
                        "history variable {variable}: branch value {} is not the density of one off-nominal moderator \
                         point of {density:?}",
                        history.branch
                    ));
                }
            };
            if history.rods.is_empty() {
                return Err(format!("history variable {variable}: no rod values"));
            }
            for (r, entry) in history.rods.iter().enumerate() {
                let at = || format!("history variable {variable}, rod value {}", entry.rod);
                if history.rods[..r].iter().any(|o| o.rod == entry.rod) {
                    return Err(format!("{}: the rod value appears twice", at()));
                }
                let rod_index = match rod.map(|k| &self.axes[k]) {
                    Some(Axis::Rod { values }) => {
                        Some(values.iter().position(|v| v == &entry.rod).ok_or_else(|| {
                            format!("{}: not a value of the rod axis {values:?}", at())
                        })?)
                    }
                    _ if entry.rod == ROD_OUT => None,
                    _ => {
                        return Err(format!(
                            "{}: tables without a rod axis hold the rod value {ROD_OUT} only",
                            at()
                        ));
                    }
                };
                // 条件の状態の点（燃焼度の番号は後で入れる）
                let mut point = vec![0; self.axes.len()];
                for (k, axis) in self.axes.iter().enumerate() {
                    point[k] = match axis {
                        Axis::Burnup { .. } => 0,
                        Axis::Rod { .. } => rod_index.unwrap_or(0),
                        _ if k == m => condition,
                        _ => axis.nominal_index().unwrap_or(0),
                    };
                }
                if !entry.burnup.windows(2).all(|w| w[0] < w[1]) {
                    return Err(format!(
                        "{}: burnup points must increase, got {:?}",
                        at(),
                        entry.burnup
                    ));
                }
                for &b in &entry.burnup {
                    let index = burnups.iter().position(|&p| p == b).ok_or_else(|| {
                        format!("{}: burnup {b} is not a point of the burnup axis", at())
                    })?;
                    point[burnup] = index;
                    if self.point_index(&point).is_none() {
                        return Err(format!(
                            "{}: the table has not solved the branch condition at burnup {b} (the residuals are \
                             taken where the table holds lattice values. ADR-0215)",
                            at()
                        ));
                    }
                }
                if entry.burnup.first() != burnups.first() || entry.burnup.last() != burnups.last()
                {
                    return Err(format!(
                        "{}: the burnup points must include the first and the last of the burnup axis (core does not \
                         extrapolate the residuals)",
                        at()
                    ));
                }
                if entry.residuals.len() != entry.burnup.len() {
                    return Err(format!(
                        "{}: {} residuals for {} burnup points",
                        at(),
                        entry.residuals.len(),
                        entry.burnup.len()
                    ));
                }
                let sized = |v: &[f64]| v.len() == groups && v.iter().all(|x| x.is_finite());
                if let Some(k) = entry.residuals.iter().position(|q| {
                    ![
                        &q.total,
                        &q.absorption,
                        &q.nu_fission,
                        &q.fission,
                        &q.chi,
                        &q.diffusion,
                    ]
                    .iter()
                    .all(|v| sized(v))
                        || q.scatter.len() != groups
                        || !q.scatter.iter().all(|row| sized(row))
                        || q.adf.len() != sides
                        || !q.adf.iter().all(|side| sized(side))
                }) {
                    return Err(format!(
                        "{}: residual {k} must have {groups} finite groups ({sides} ADF sides)",
                        at()
                    ));
                }
            }
        }
        Ok(())
    }

    /// ミクロ燃焼のデータ（ADR-0216）が、燃焼度の軸と点と群の数に合い、値が有限で、核種ごとの反応の並びと核分裂の有無と数密度の
    /// 微分（ADR-0221）の有無が状態の間で同じこと（状態の間で内挿するため）
    fn validate_micro(&self, micro: &Micro) -> Result<(), String> {
        let burnups = self
            .axes
            .iter()
            .find(|a| matches!(a, Axis::Burnup { .. }))
            .ok_or("micro data need a burnup axis (they are written from a depletion history)")?
            .len();
        if self.split != 1 {
            return Err(format!(
                "micro data are homogenized over the whole lattice, so the table must not be split (split {})",
                self.split
            ));
        }
        let names = &micro.nuclides;
        if names.is_empty() {
            return Err("micro data need nuclides".into());
        }
        if let Some(name) = names
            .iter()
            .enumerate()
            .find(|&(i, n)| !crate::library::is_nuclide_name(n) || names[..i].contains(n))
            .map(|(_, n)| n)
        {
            return Err(format!(
                "micro nuclide {name:?} is not a distinct nuclide name"
            ));
        }
        if micro.densities.len() != burnups
            || micro.densities.iter().any(|row| {
                row.len() != names.len() || !row.iter().all(|&n| n >= 0.0 && n.is_finite())
            })
        {
            return Err(format!(
                "micro densities must be {burnups} rows of {} nonnegative values",
                names.len()
            ));
        }
        if micro.states.len() != self.points.len() {
            return Err(format!(
                "{} micro states for {} points",
                micro.states.len(),
                self.points.len()
            ));
        }
        let groups = self.states[0].bounds.len() - 1;
        let sized = |v: &[f64]| v.len() == groups && v.iter().all(|x| x.is_finite());
        let nonnegative = |v: &[f64]| sized(v) && v.iter().all(|&x| x >= 0.0);
        // 状態の間で比べる、核種ごとの反応の MT と経路の数、核分裂の有無、微分の有無
        let shape = |n: &MicroNuclide| {
            (
                n.reactions
                    .iter()
                    .map(|r| (r.mt, r.channels.len()))
                    .collect::<Vec<_>>(),
                n.thermal_fission.is_some(),
                n.derivative.is_some(),
            )
        };
        for (p, state) in micro.states.iter().enumerate() {
            if state.len() != names.len() {
                return Err(format!(
                    "micro state {p}: {} nuclides for {}",
                    state.len(),
                    names.len()
                ));
            }
            for (i, n) in state.iter().enumerate() {
                let at = || format!("micro state {p}, {}", names[i]);
                let fissile = n.fission.iter().any(|&f| f > 0.0);
                if !nonnegative(&n.total)
                    || !nonnegative(&n.nu_fission)
                    || !nonnegative(&n.fission)
                    || n.scatter.len() != groups
                    || !n.scatter.iter().all(|row| sized(row))
                {
                    return Err(format!(
                        "{}: cross sections must have {groups} finite groups (total and fission nonnegative)",
                        at()
                    ));
                }
                // 微分は負になりうる（組成が増えると減る少数群定数）ので、群の数と有限だけを確かめる
                if let Some(d) = &n.derivative
                    && !(sized(&d.total)
                        && sized(&d.nu_fission)
                        && sized(&d.fission)
                        && d.scatter.len() == groups
                        && d.scatter.iter().all(|row| sized(row)))
                {
                    return Err(format!(
                        "{}: density derivatives must have {groups} finite groups",
                        at()
                    ));
                }
                let mts: Vec<u32> = n.reactions.iter().map(|r| r.mt).collect();
                if mts.iter().enumerate().any(|(k, mt)| mts[..k].contains(mt))
                    || !n
                        .reactions
                        .iter()
                        .all(|r| r.channels.iter().all(|c| nonnegative(c)))
                {
                    return Err(format!(
                        "{}: reactions must be distinct and their channels {groups} nonnegative groups",
                        at()
                    ));
                }
                match n.thermal_fission {
                    Some(f) if fissile && (0.0..=1.0).contains(&f) => {}
                    None if !fissile => {}
                    _ => {
                        return Err(format!(
                            "{}: the fraction of fissions below the yield cutoff must be in 0..1 for fissile nuclides \
                             and nil for the others, got {:?}",
                            at(),
                            n.thermal_fission
                        ));
                    }
                }
                if shape(n) != shape(&micro.states[0][i]) {
                    return Err(format!(
                        "{}: the reactions, the fission and the density derivative differ from those of state 0",
                        at()
                    ));
                }
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
            micro: self.micro.clone(),
            histories: self.histories.clone(),
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
            micro: file.micro,
            histories: file.histories,
        };
        table.validate()?;
        Ok(table)
    }
}

#[cfg(test)]
mod tests;
