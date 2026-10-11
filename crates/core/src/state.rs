//! 炉心の状態からの分岐のテーブルの内挿（ADR-0097）。
//!
//! 軸ごとに隣の2点で線形に内挿する（ホウ素濃度は ppm、減速材は温度、燃料は √T）。格子点の範囲の外はエラーにし、
//! 外挿しない。点が1つの軸は、状態がその値に等しいときだけ使える。状態を与えるときはテーブルが3つの軸をすべて
//! 持ち、与えないときは軸のないテーブルだけを使う。

use serde::Deserialize;
use xlc_format::fewgroup::{FewGroup, Pins};
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

/// `w` 倍した少数群定数（群境界と f* の印はそのまま）。項目を足したらここで扱いを決める（`..` を使わない）
fn scaled(f: &FewGroup, w: f64) -> FewGroup {
    let FewGroup {
        bounds,
        k_inf,
        total,
        absorption,
        nu_fission,
        fission,
        chi,
        diffusion,
        scatter,
        adf,
        vacuum_ratio,
        normalized_adf,
        pins,
    } = f;
    let v = |x: &[f64]| x.iter().map(|x| w * x).collect::<Vec<_>>();
    FewGroup {
        bounds: bounds.clone(),
        k_inf: w * k_inf,
        total: v(total),
        absorption: v(absorption),
        nu_fission: v(nu_fission),
        fission: v(fission),
        chi: v(chi),
        diffusion: v(diffusion),
        scatter: scatter.iter().map(|r| v(r)).collect(),
        adf: adf.iter().map(|r| v(r)).collect(),
        vacuum_ratio: vacuum_ratio.iter().map(|r| r.as_deref().map(v)).collect(),
        normalized_adf: normalized_adf.clone(),
        // セルの境界はそのまま、形状関数と角の不連続因子は w 倍
        pins: pins.as_ref().map(|p| Pins {
            dx: p.dx.clone(),
            dy: p.dy.clone(),
            form: p
                .form
                .iter()
                .map(|g| g.iter().map(|row| v(row)).collect())
                .collect(),
            corners: p.corners.iter().map(|c| v(c)).collect(),
            normalized_corners: p.normalized_corners.clone(),
            halves: p.halves.clone(),
        }),
    }
}

/// `acc += w * f`。項目を足したらここで扱いを決める（`..` を使わない）。真空の辺、形状関数の有無とセルは
/// どの状態でも同じはず（テーブルの検査）で、違えばエラー
fn accumulate(acc: &mut FewGroup, w: f64, f: &FewGroup) -> Result<(), String> {
    let FewGroup {
        bounds: _,
        k_inf,
        total,
        absorption,
        nu_fission,
        fission,
        chi,
        diffusion,
        scatter,
        adf,
        vacuum_ratio,
        normalized_adf: _,
        pins,
    } = f;
    acc.k_inf += w * k_inf;
    add(&mut acc.total, w, total);
    add(&mut acc.absorption, w, absorption);
    add(&mut acc.nu_fission, w, nu_fission);
    add(&mut acc.fission, w, fission);
    add(&mut acc.chi, w, chi);
    add(&mut acc.diffusion, w, diffusion);
    for (a, x) in acc.scatter.iter_mut().zip(scatter) {
        add(a, w, x);
    }
    for (a, x) in acc.adf.iter_mut().zip(adf) {
        add(a, w, x);
    }
    for (a, x) in acc.vacuum_ratio.iter_mut().zip(vacuum_ratio) {
        match (a, x) {
            (Some(a), Some(x)) => add(a, w, x),
            (None, None) => {}
            _ => return Err("the table states differ in their vacuum sides".into()),
        }
    }
    match (&mut acc.pins, pins) {
        (Some(a), Some(x))
            if a.dx == x.dx
                && a.dy == x.dy
                && a.normalized_corners == x.normalized_corners
                && a.halves == x.halves =>
        {
            for (a, x) in a.form.iter_mut().flatten().zip(x.form.iter().flatten()) {
                add(a, w, x);
            }
            for (a, x) in a.corners.iter_mut().zip(&x.corners) {
                add(a, w, x);
            }
        }
        (None, None) => {}
        _ => return Err("the table states differ in their pins (form functions)".into()),
    }
    Ok(())
}

/// `acc += w * x`（同じ形の量どうし）
fn add(acc: &mut [f64], w: f64, x: &[f64]) {
    for (a, x) in acc.iter_mut().zip(x) {
        *a += w * x;
    }
}

/// テーブルを炉心の状態で内挿した、部分ノードごと（`Table::states` の部分ノードの順。分けていなければ1つ）の少数群
/// 定数。`state` が `None` なら、軸のないテーブルの状態をそのまま返す。
pub fn interpolate(table: &Table, state: Option<&State>) -> Result<Vec<FewGroup>, String> {
    let sub = table.sub_nodes();
    let Some(state) = state else {
        if !table.axes.is_empty() {
            return Err("the table has branch axes; give the core [state]".into());
        }
        if table.states.len() != sub || sub == 0 {
            return Err("the table has no states".into());
        }
        return Ok(table.states.clone());
    };
    // 軸がないことを、ほかの軸の範囲の外より先に示す（テーブルの検査で同じ軸は2つない）
    for name in ["boron", "moderator", "fuel"] {
        if !table.axes.iter().any(|a| a.name() == name) {
            return Err(format!("the table has no {name} axis for the core state"));
        }
    }
    // テーブルの軸の順に、隣の2点
    let mut brackets = Vec::with_capacity(table.axes.len());
    for axis in &table.axes {
        brackets.push(match axis {
            Axis::Boron { ppm, b10 } => {
                if *b10 != state.b10 {
                    return Err(format!(
                        "boron: the table has b10 {b10} and the state {}; they must be exactly equal",
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

    let mut out: Vec<FewGroup> = (0..sub)
        .map(|q| table.states.get(q).map(|s| scaled(s, 0.0)))
        .collect::<Option<_>>()
        .ok_or("the table has no states")?;
    // 隣の 2^d 個の格子点の重みつきの和（最初の軸がいちばん遅く変わる並び）。部分ノードごとに足す
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
        for (q, acc) in out.iter_mut().enumerate() {
            let s = table
                .states
                .get(index * sub + q)
                .ok_or_else(|| format!("the table has no state {index} sub-node {q}"))?;
            accumulate(acc, weight, s)?;
        }
    }
    // 凸結合なので値の検査は通るはずだが、ソルバに渡す前に確かめる
    for f in &out {
        f.validate().map_err(|e| e.to_string())?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests;
