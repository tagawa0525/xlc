//! 履歴の補正の残り R（ADR-0215）。公称の燃焼のテーブルと、履歴の分岐の燃焼の履歴を自分の条件で解いた分岐のテーブル（`[history]`）と、
//! 公称の燃焼の組成のうち追う核種だけを履歴の分岐の値にした組成（ハイブリッド）を同じ条件で解いたテーブルから、燃焼度の点 u ごとに
//!
//! U(u) = f(履歴の分岐の組成、条件 j) − f(ハイブリッドの組成、条件 j)（追わない核種の分）
//! L(u) = f(ハイブリッドの組成、条件 j) − f(公称の組成、条件 j) − Σᵢ cᵢ(条件 j、u) (N̄ᵢ^j(u) − N̄ᵢ^0(u))（線形化の誤差）
//!
//! を作り（和が残り R）、履歴の分岐の均質化した数密度 N̄ᵢ^j(u) と一緒に公称のテーブルの `histories` に書き足す。L の第 2 項と係数 cᵢ は
//! 公称のテーブルの条件 j の状態の値で、cᵢ は core の線形の補正と同じ関数（`xlc_format::table::add_linear_correction`。数密度の微分を
//! 持つ核種は ∂Σ/∂N̄ᵢ、ほかは σᵢ。ADR-0216、ADR-0221）。N̄ᵢ^j と N̄ᵢ^0 はそれぞれのテーブルの均質化した数密度（ハイブリッドの N̄ᵢ は
//! 分岐と同じ）。補正を持たない項目（χ、拡散係数、ADF）の L は 2 つの定数の差。
//!
//! 同じ点ごとに、追う核種のミクロ燃焼の反応率の微視的断面積の差 Δσᵢ(u) = σᵢ(履歴の分岐の組成、条件 j) − σᵢ(公称の組成、条件 j) も書く
//! （ADR-0215 の案 S2。核分裂、チェーンの反応の経路、収率の境界より下の核分裂の断面積）。
//!
//! 使い方: lattice histories <入力ファイル>（ファイルのパスは入力ファイルのあるディレクトリから）
//!
//! ```toml
//! nominal = "table-nominal.msgpack"            # 公称のテーブル（燃焼度の軸、公称の値、ミクロ燃焼のデータを持つ）
//! output = "table-nominal-histories.msgpack"   # R を書き足した公称のテーブル
//!
//! [[history]]
//! variable = "moderator-density"               # 履歴の変数
//! # 制御棒の値ごとの、履歴の分岐の自分の条件のテーブルと、同じ条件のハイブリッドのテーブル
//! branches = { out = { branch = "table-density-branch.msgpack", hybrid = "table-hybrid-density-branch.msgpack" } }
//! ```
//!
//! 履歴の分岐のテーブルは、燃焼度の軸と、点が 1 つの分岐の軸（制御棒の軸は持てば値が 1 つ）を持ち、その状態（条件 j）が公称の
//! テーブルの格子点で、変数の軸（減速材の密度は `moderator`）だけが公称と違うこと。燃焼度の点は、公称のテーブルが条件 j の状態を
//! 解いた燃焼度の点と同じこと。ハイブリッドのテーブルは、分岐のテーブルと同じ軸と点、追う核種、初めの重金属、群境界と、同じ均質化した
//! 数密度を持つこと。変数の値は、公称の燃焼が公称の格子点の値 H⁰、履歴の分岐が条件 j の格子点の値 Hᵇ。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, anyhow, bail};
use serde::Deserialize;
use xlc_format::fewgroup::FewGroup;
use xlc_format::table::{
    Axis, HistoryCondition, HistoryCorrection, HistoryRod, MicroNuclide, MicroReaction, ROD_OUT,
    RateResidual, Residual, Table, add_linear_correction,
};

use crate::checkpoint::write_synced;

/// `lattice histories` の入力。
#[derive(Debug, Clone, PartialEq)]
pub struct ResidualInput {
    /// 公称のテーブル
    pub nominal: PathBuf,
    /// R を書き足した公称のテーブルの出力
    pub output: PathBuf,
    /// 履歴の変数ごとの、制御棒の値と履歴の分岐の自分の条件のテーブル
    pub histories: Vec<HistoryEntry>,
}

/// 1 つの履歴の変数の入力。
#[derive(Debug, Clone, PartialEq)]
pub struct HistoryEntry {
    pub variable: String,
    /// (制御棒の値, その値のテーブル)。制御棒の値の名前の順
    pub branches: Vec<(String, BranchFiles)>,
}

/// 1 つの制御棒の値の、履歴の分岐の自分の条件のテーブルと、同じ条件のハイブリッドのテーブル
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BranchFiles {
    pub branch: PathBuf,
    pub hybrid: PathBuf,
}

/// 1 つの制御棒の値の、履歴の分岐の自分の条件のテーブルと、追う核種だけを履歴の分岐の値にした公称の組成を同じ条件で解いたテーブル
#[derive(Debug, Clone, PartialEq)]
pub struct BranchTables {
    pub branch: Table,
    pub hybrid: Table,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InputFile {
    nominal: PathBuf,
    output: PathBuf,
    history: Vec<HistoryFile>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct HistoryFile {
    variable: String,
    branches: BTreeMap<String, BranchFiles>,
}

impl ResidualInput {
    /// TOML の入力を読む。履歴と、変数ごとの分岐のテーブルが 1 つ以上あり、出力が公称のテーブルと別のファイルであること
    pub fn parse(text: &str) -> anyhow::Result<Self> {
        let file: InputFile = toml::from_str(text)?;
        if file.history.is_empty() {
            bail!("give a [[history]]");
        }
        // 読んだテーブルを上書きすると、R を作る前の公称のテーブルが残らない
        if file.output == file.nominal {
            bail!("output must differ from the nominal table");
        }
        let histories = file
            .history
            .into_iter()
            .map(|h| {
                if h.branches.is_empty() {
                    bail!(
                        "history {}: give the branch tables per rod value",
                        h.variable
                    );
                }
                Ok(HistoryEntry {
                    variable: h.variable,
                    branches: h.branches.into_iter().collect(),
                })
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        Ok(Self {
            nominal: file.nominal,
            output: file.output,
            histories,
        })
    }
}

/// 公称のテーブル `nominal` の、変数 `variable` の履歴の補正を、制御棒の値ごとの履歴の分岐の自分の条件のテーブルとハイブリッドの
/// テーブル `branches` から作る。変数の値 H⁰ と Hᵇ は制御棒の値によらず同じこと
pub fn history_correction(
    nominal: &Table,
    variable: &str,
    branches: &[(String, BranchTables)],
) -> anyhow::Result<HistoryCorrection> {
    let mut values: Option<(f64, f64)> = None;
    let mut rods = Vec::with_capacity(branches.len());
    for (rod, branch) in branches {
        let (h0, hb, entry) = residuals(nominal, variable, rod, branch)
            .with_context(|| format!("history {variable}, rod value {rod}"))?;
        match values {
            Some(v) if v != (h0, hb) => bail!(
                "history {variable}: rod value {rod} gives the values {:?}, another gives {v:?}",
                (h0, hb)
            ),
            _ => values = Some((h0, hb)),
        }
        rods.push(entry);
    }
    let (nominal_value, branch_value) =
        values.ok_or_else(|| anyhow!("history {variable}: no branch tables"))?;
    Ok(HistoryCorrection {
        variable: variable.into(),
        nominal: nominal_value,
        branch: branch_value,
        rods,
    })
}

/// 分岐の軸（ホウ素、減速材、燃料）の、点が 1 つのテーブルのその値（ホウ素は (ppm, B-10 の割合)、減速材は (温度, 密度)、燃料は温度）
#[derive(Debug, Clone, Copy, PartialEq)]
enum Value {
    Boron(f64, f64),
    Moderator(f64, f64),
    Fuel(f64),
}

/// 軸 `axis` の格子点 `i` の値（燃焼度と制御棒の軸は `None`）
fn value(axis: &Axis, i: usize) -> Option<Value> {
    match axis {
        Axis::Boron { ppm, b10, .. } => Some(Value::Boron(ppm[i], *b10)),
        Axis::Moderator {
            temperature,
            density,
            ..
        } => Some(Value::Moderator(temperature[i], density[i])),
        Axis::Fuel { temperature, .. } => Some(Value::Fuel(temperature[i])),
        Axis::Burnup { .. } | Axis::Rod { .. } => None,
    }
}

/// 公称のテーブル `nominal` の制御棒の値 `rod` の残りを、履歴の分岐の自分の条件のテーブルとハイブリッドのテーブル `tables` から作る。
/// 返すのは (H⁰, Hᵇ, 残り)
pub fn residuals(
    nominal: &Table,
    variable: &str,
    rod: &str,
    tables: &BranchTables,
) -> anyhow::Result<(f64, f64, HistoryRod)> {
    let BranchTables { branch, hybrid } = tables;
    let micro = nominal
        .micro
        .as_ref()
        .ok_or_else(|| anyhow!("the nominal table has no micro data (ADR-0216)"))?;
    let branch_micro = branch
        .micro
        .as_ref()
        .ok_or_else(|| anyhow!("the branch table has no micro data (ADR-0216)"))?;
    if branch_micro.nuclides != micro.nuclides {
        bail!(
            "the branch table tracks the micro nuclides {:?}, the nominal table {:?}",
            branch_micro.nuclides,
            micro.nuclides
        );
    }
    // 同じ初めの燃料の燃焼の履歴どうし（燃焼度の点が同じ組成の量を指す）
    if branch.heavy_metal.is_none() || branch.heavy_metal != nominal.heavy_metal {
        bail!(
            "the branch table has the initial heavy metal {:?} g/cm, the nominal table {:?}",
            branch.heavy_metal,
            nominal.heavy_metal
        );
    }
    if branch.states[0].bounds != nominal.states[0].bounds {
        bail!("the branch and the nominal tables have other group bounds");
    }
    check_hybrid(branch, hybrid)?;
    // 分岐のテーブルの条件 j: 燃焼度のほかの軸の点が 1 つ（制御棒の軸は持てば値が `rod`）
    let mut condition: Vec<Value> = Vec::new();
    for axis in &branch.axes {
        match axis {
            Axis::Burnup { .. } => {}
            Axis::Rod { values } => {
                if values.as_slice() != [rod.to_string()] {
                    bail!("the branch table holds the rod values {values:?}, not one value {rod}");
                }
            }
            _ if axis.len() != 1 => bail!(
                "the branch table must be one state per burnup; its {} axis has {} points",
                axis.name(),
                axis.len()
            ),
            _ => condition.extend(value(axis, 0)),
        }
    }
    let branch_rod = branch.axes.iter().any(|a| matches!(a, Axis::Rod { .. }));
    if !branch_rod && rod != ROD_OUT {
        bail!("the branch table has no rod axis, so its rod value is {ROD_OUT}, not {rod}");
    }
    // 変数の値 Hᵇ は、分岐のテーブルの条件の減速材の密度
    let hb = condition
        .iter()
        .find_map(|v| match v {
            Value::Moderator(_, density) => Some(*density),
            _ => None,
        })
        .ok_or_else(|| anyhow!("the branch table has no moderator axis for {variable}"))?;
    // 公称のテーブルの条件 j の格子点（テーブルの検査と同じ規則）
    let HistoryCondition {
        nominal: h0,
        mut point,
        burnup_axis,
    } = nominal
        .history_condition(variable, hb, rod)
        .map_err(|e| anyhow!(e))?;
    // 分岐のテーブルの条件の値（減速材は温度と密度の組）が、その格子点の値と同じこと。公称のテーブルに無い軸の値も、分岐の
    // テーブルに無い軸の値も推測しない（違う値の効きが R に混ざる）
    let at_point: Vec<Value> = nominal
        .axes
        .iter()
        .zip(&point)
        .filter_map(|(axis, &i)| value(axis, i))
        .collect();
    let name = |v: &Value| match v {
        Value::Boron(..) => "boron",
        Value::Moderator(..) => "moderator",
        Value::Fuel(..) => "fuel",
    };
    if let Some(lost) = condition.iter().find(|v| !at_point.contains(v)) {
        bail!(
            "the branch condition {lost:?} on the {} axis is not the value at the condition point {at_point:?} of \
             the nominal table (the moderator at the density {hb}, the other axes nominal)",
            name(lost)
        );
    }
    if let Some(missing) = at_point.iter().find(|v| !condition.contains(v)) {
        bail!(
            "the branch table has no {} axis for the condition value {missing:?} of the nominal table",
            name(missing)
        );
    }
    // 公称のテーブルが条件 j を解いた燃焼度の点が、分岐のテーブルの燃焼度の点
    let burnups = nominal.axes[burnup_axis].points();
    let solved: Vec<usize> = (0..burnups.len())
        .filter(|&u| {
            point[burnup_axis] = u;
            nominal.point_index(&point).is_some()
        })
        .collect();
    let solved_values: Vec<f64> = solved.iter().map(|&u| burnups[u]).collect();
    let branch_burnups = branch
        .axes
        .iter()
        .find(|a| matches!(a, Axis::Burnup { .. }))
        .ok_or_else(|| anyhow!("the branch table has no burnup axis"))?
        .points();
    if branch_burnups != solved_values.as_slice() {
        bail!(
            "the branch table has the burnup points {branch_burnups:?}, the nominal table solved the branch \
             condition at {solved_values:?}"
        );
    }
    // 分岐のテーブルは燃焼度の点ごとに状態 1 つで、点は燃焼度の順
    let rates = solved
        .iter()
        .enumerate()
        .map(|(k, &u)| {
            point[burnup_axis] = u;
            let p = nominal
                .point_index(&point)
                .ok_or_else(|| anyhow!("no nominal point {point:?}"))?;
            branch_micro.states[k]
                .iter()
                .zip(&micro.states[p])
                .zip(&micro.nuclides)
                .map(|((b, n), name)| {
                    rate_difference(b, n).with_context(|| format!("burnup {}, {name}", burnups[u]))
                })
                .collect::<anyhow::Result<Vec<_>>>()
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let linearization = solved
        .iter()
        .enumerate()
        .map(|(k, &u)| {
            point[burnup_axis] = u;
            let p = nominal
                .point_index(&point)
                .ok_or_else(|| anyhow!("no nominal point {point:?}"))?;
            let delta: Vec<f64> = branch_micro.densities[k]
                .iter()
                .zip(&micro.densities[u])
                .map(|(nj, n0)| nj - n0)
                .collect();
            let mut corrected = nominal.states[p].clone();
            add_linear_correction(&mut corrected, &micro.states[p], &delta);
            Ok(difference(&hybrid.states[k], &corrected))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let untracked = branch
        .states
        .iter()
        .zip(&hybrid.states)
        .map(|(b, h)| difference(b, h))
        .collect();
    Ok((
        h0,
        hb,
        HistoryRod {
            rod: rod.into(),
            burnup: solved_values,
            untracked,
            linearization,
            densities: branch_micro.densities.clone(),
            micro: rates,
        },
    ))
}

/// ハイブリッドのテーブル `hybrid` が、分岐のテーブル `branch` と同じ条件と燃焼度の点（軸と点）、追う核種、初めの重金属、群境界と、
/// 同じ均質化した数密度を持つこと（U と L は同じ状態の 2 つの定数の差で、ハイブリッドは追う核種が分岐と同じ組成。体積の平均の
/// N̄ᵢ が違えば、置き換えた核種か履歴が違う）
fn check_hybrid(branch: &Table, hybrid: &Table) -> anyhow::Result<()> {
    if hybrid.axes != branch.axes || hybrid.points != branch.points {
        bail!("the hybrid table solves other conditions or burnup points than the branch table");
    }
    let micro = hybrid
        .micro
        .as_ref()
        .ok_or_else(|| anyhow!("the hybrid table has no micro data (ADR-0216)"))?;
    let branch_micro = branch
        .micro
        .as_ref()
        .ok_or_else(|| anyhow!("the branch table has no micro data (ADR-0216)"))?;
    if micro.nuclides != branch_micro.nuclides {
        bail!(
            "the hybrid table tracks the micro nuclides {:?}, the branch table {:?}",
            micro.nuclides,
            branch_micro.nuclides
        );
    }
    if micro.densities != branch_micro.densities {
        bail!("the hybrid table holds other tracked densities than the branch table");
    }
    if hybrid.heavy_metal != branch.heavy_metal {
        bail!(
            "the hybrid table has the initial heavy metal {:?} g/cm, the branch table {:?}",
            hybrid.heavy_metal,
            branch.heavy_metal
        );
    }
    if hybrid.states[0].bounds != branch.states[0].bounds {
        bail!("the hybrid and the branch tables have other group bounds");
    }
    Ok(())
}

/// 反応率の微視的断面積の差 `a − b`（`a` は履歴の分岐、`b` は公称の条件 j の状態）。収率の境界より下の割合は、割合 × 群の和の σf の
/// 差にする（core が補正した σf の和との比を割合にする）。反応の並び（MT と経路の数）か核分裂の有無が違えばエラー
fn rate_difference(a: &MicroNuclide, b: &MicroNuclide) -> anyhow::Result<RateResidual> {
    let shape = |n: &MicroNuclide| {
        n.reactions
            .iter()
            .map(|r| (r.mt, r.channels.len()))
            .collect::<Vec<_>>()
    };
    if shape(a) != shape(b) {
        bail!(
            "the branch table has the reactions {:?}, the nominal table {:?}",
            shape(a),
            shape(b)
        );
    }
    let minus = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(p, q)| p - q).collect::<Vec<_>>();
    let below = |n: &MicroNuclide, f: f64| f * n.fission.iter().sum::<f64>();
    let thermal_fission = match (a.thermal_fission, b.thermal_fission) {
        (Some(fa), Some(fb)) => Some(below(a, fa) - below(b, fb)),
        (None, None) => None,
        (fa, fb) => bail!(
            "the branch table has the fraction of fissions below the yield cutoff {fa:?}, the nominal table {fb:?}"
        ),
    };
    Ok(RateResidual {
        fission: minus(&a.fission, &b.fission),
        reactions: a
            .reactions
            .iter()
            .zip(&b.reactions)
            .map(|(x, y)| MicroReaction {
                mt: x.mt,
                channels: x
                    .channels
                    .iter()
                    .zip(&y.channels)
                    .map(|(p, q)| minus(p, q))
                    .collect(),
            })
            .collect(),
        thermal_fission,
    })
}

/// 少数群定数の差 `a − b` の、残りの項目
fn difference(a: &FewGroup, b: &FewGroup) -> Residual {
    let minus = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(p, q)| p - q).collect::<Vec<_>>();
    let rows = |x: &[Vec<f64>], y: &[Vec<f64>]| {
        x.iter()
            .zip(y)
            .map(|(p, q)| minus(p, q))
            .collect::<Vec<_>>()
    };
    Residual {
        total: minus(&a.total, &b.total),
        absorption: minus(&a.absorption, &b.absorption),
        nu_fission: minus(&a.nu_fission, &b.nu_fission),
        fission: minus(&a.fission, &b.fission),
        chi: minus(&a.chi, &b.chi),
        diffusion: minus(&a.diffusion, &b.diffusion),
        scatter: rows(&a.scatter, &b.scatter),
        adf: rows(&a.adf, &b.adf),
    }
}

/// 入力ファイル `path` を読み、R を書き足した公称のテーブルを書く。公称のテーブルが既に R を持てばエラー（重ねて足さない）
pub fn run(path: &Path) -> anyhow::Result<()> {
    let shown = |p: &Path| p.display().to_string();
    let text = std::fs::read_to_string(path).with_context(|| shown(path))?;
    let input = ResidualInput::parse(&text).with_context(|| shown(path))?;
    let base = path.parent().unwrap_or(Path::new("."));
    let read = |file: &Path| -> anyhow::Result<Table> {
        let file = base.join(file);
        let bytes = std::fs::read(&file).with_context(|| shown(&file))?;
        Table::from_msgpack(&bytes).with_context(|| shown(&file))
    };
    let mut nominal = read(&input.nominal)?;
    if !nominal.histories.is_empty() {
        bail!(
            "{}: the nominal table already holds history residuals",
            shown(&input.nominal)
        );
    }
    let mut corrections = Vec::with_capacity(input.histories.len());
    for history in &input.histories {
        let branches = history
            .branches
            .iter()
            .map(|(rod, files)| {
                Ok((
                    rod.clone(),
                    BranchTables {
                        branch: read(&files.branch)?,
                        hybrid: read(&files.hybrid)?,
                    },
                ))
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        corrections.push(history_correction(&nominal, &history.variable, &branches)?);
    }
    nominal.histories = corrections;
    let output = base.join(&input.output);
    // 読んだテーブルを別の書き方のパスで上書きしない（作り直しに格子計算が要る）。まだ無い出力は読んだファイルでない
    if output.exists() {
        let canonical = |p: &Path| std::fs::canonicalize(p).with_context(|| shown(p));
        let target = canonical(&output)?;
        for file in std::iter::once(&input.nominal).chain(
            input
                .histories
                .iter()
                .flat_map(|h| h.branches.iter())
                .flat_map(|(_, f)| [&f.branch, &f.hybrid]),
        ) {
            if canonical(&base.join(file))? == target {
                bail!(
                    "{}: the output is the table {}, which it reads",
                    shown(&output),
                    shown(file)
                );
            }
        }
    }
    let bytes = nominal.to_msgpack().with_context(|| shown(&output))?;
    write_synced(&output, &bytes)
}

#[cfg(test)]
mod tests;
