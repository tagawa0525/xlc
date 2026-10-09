//! 炉心の状態からの分岐のテーブルの内挿（ADR-0097）。
//!
//! 軸ごとに隣の2点で線形に内挿する（ホウ素濃度は ppm、減速材は温度、燃料は √T）。格子点の範囲の外はエラーにし、
//! 外挿しない。点が1つの軸は、状態がその値に等しいときだけ使える。状態を与えるときはテーブルが3つの軸をすべて
//! 持ち、与えないときは軸のないテーブルだけを使う。

use serde::Deserialize;
use xlc_format::fewgroup::FewGroup;
use xlc_format::table::{Axis, Table};

/// 炉心の状態（ケースの `[state]`）。
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    /// ホウ素濃度 (ppm)
    pub boron: f64,
    /// ホウ素の中の B-10 の原子分率（テーブルの値と同じであること）
    pub b10: f64,
    /// 減速材の温度 (K)。密度はテーブルの格子点の組で決まる
    pub moderator: f64,
    /// 燃料の温度 (K)
    pub fuel: f64,
}

/// 内挿する軸の、隣の2点の番号と、上の点の重み。
struct Bracket {
    low: usize,
    high: usize,
    weight: f64,
}

/// `points` の中で `x` を挟む2点（`transform` をかけた値で線形）。
fn bracket(
    name: &str,
    points: &[f64],
    x: f64,
    transform: impl Fn(f64) -> f64,
) -> Result<Bracket, String> {
    let out = || format!("{name}: state {x} is outside the table points {points:?}");
    if !x.is_finite() {
        return Err(out());
    }
    if let [only] = points {
        return if x == *only {
            Ok(Bracket {
                low: 0,
                high: 0,
                weight: 0.0,
            })
        } else {
            Err(out())
        };
    }
    let i = points
        .windows(2)
        .position(|w| w[0] <= x && x <= w[1])
        .ok_or_else(out)?;
    let (a, b) = (transform(points[i]), transform(points[i + 1]));
    Ok(Bracket {
        low: i,
        high: i + 1,
        weight: (transform(x) - a) / (b - a),
    })
}

/// `acc += w * x`（同じ形の量どうし）
fn add(acc: &mut [f64], w: f64, x: &[f64]) {
    for (a, x) in acc.iter_mut().zip(x) {
        *a += w * x;
    }
}

/// テーブルを炉心の状態で内挿した少数群定数。`state` が `None` なら、軸のないテーブルの状態をそのまま返す。
pub fn interpolate(table: &Table, state: Option<&State>) -> Result<FewGroup, String> {
    let Some(state) = state else {
        if !table.axes.is_empty() {
            return Err("the table has branch axes; give the core [state]".into());
        }
        return table
            .states
            .first()
            .cloned()
            .ok_or_else(|| "the table has no states".to_string());
    };
    // テーブルの軸の順に、隣の2点
    let mut brackets = Vec::with_capacity(table.axes.len());
    for axis in &table.axes {
        brackets.push(match axis {
            Axis::Boron { ppm, b10 } => {
                if *b10 != state.b10 {
                    return Err(format!(
                        "boron: the table has b10 {b10}, the state {}",
                        state.b10
                    ));
                }
                bracket("boron", ppm, state.boron, |x| x)?
            }
            Axis::Moderator { temperature, .. } => {
                bracket("moderator", temperature, state.moderator, |x| x)?
            }
            Axis::Fuel { temperature } => bracket("fuel", temperature, state.fuel, f64::sqrt)?,
        });
    }
    for name in ["boron", "moderator", "fuel"] {
        if !table.axes.iter().any(|a| a.name() == name) {
            return Err(format!("the table has no {name} axis for the core state"));
        }
    }

    let first = table.states.first().ok_or("the table has no states")?;
    let zeros = |v: &[f64]| vec![0.0; v.len()];
    let mut out = FewGroup {
        bounds: first.bounds.clone(),
        k_inf: 0.0,
        total: zeros(&first.total),
        absorption: zeros(&first.absorption),
        nu_fission: zeros(&first.nu_fission),
        fission: zeros(&first.fission),
        chi: zeros(&first.chi),
        diffusion: zeros(&first.diffusion),
        scatter: first.scatter.iter().map(|r| zeros(r)).collect(),
        adf: first.adf.iter().map(|r| zeros(r)).collect(),
        vacuum_ratio: first
            .vacuum_ratio
            .iter()
            .map(|r| r.as_deref().map(zeros))
            .collect(),
        normalized_adf: first.normalized_adf.clone(),
    };
    // 隣の 2^d 個の格子点の重みつきの和（最初の軸がいちばん遅く変わる並び）
    for corner in 0..1usize << brackets.len() {
        let mut index = 0;
        let mut weight = 1.0;
        for (k, (axis, b)) in table.axes.iter().zip(&brackets).enumerate() {
            let upper = corner >> (brackets.len() - 1 - k) & 1 == 1;
            index = index * axis.points().len() + if upper { b.high } else { b.low };
            weight *= if upper { b.weight } else { 1.0 - b.weight };
        }
        if weight == 0.0 {
            continue;
        }
        let s = table
            .states
            .get(index)
            .ok_or_else(|| format!("the table has no state {index}"))?;
        out.k_inf += weight * s.k_inf;
        add(&mut out.total, weight, &s.total);
        add(&mut out.absorption, weight, &s.absorption);
        add(&mut out.nu_fission, weight, &s.nu_fission);
        add(&mut out.fission, weight, &s.fission);
        add(&mut out.chi, weight, &s.chi);
        add(&mut out.diffusion, weight, &s.diffusion);
        for (acc, x) in out.scatter.iter_mut().zip(&s.scatter) {
            add(acc, weight, x);
        }
        for (acc, x) in out.adf.iter_mut().zip(&s.adf) {
            add(acc, weight, x);
        }
        // 真空の辺はどの状態でも同じ（テーブルの検査）
        for (acc, x) in out.vacuum_ratio.iter_mut().zip(&s.vacuum_ratio) {
            if let (Some(acc), Some(x)) = (acc, x) {
                add(acc, weight, x);
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests;
