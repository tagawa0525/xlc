//! 炉心の状態からの分岐のテーブルの内挿（ADR-0097）。
//!
//! 軸ごとに隣の2点で線形に内挿する（ホウ素濃度は ppm、減速材は温度、燃料は √T）。格子点の範囲の外はエラーにし、
//! 外挿しない。点が1つの軸は、状態がその値に等しいときだけ使える。状態を与えるときはテーブルが3つの軸をすべて
//! 持ち、与えないときは軸のないテーブルだけを使う。制御棒の離散の軸は内挿せず、値の名前で選ぶ（棒がなければ
//! `out`。ADR-0127）。燃焼度は、各点の傾きをその点と両隣の 3 点の放物線の微分にした 3 次 Hermite で内挿する
//! （ADR-0192）。燃焼度 0 の点は傾きに使わず、0 から次の点までは線形にする。傾きに使う点が 2 つの軸は線形になる。
//! 燃焼度を与える状態では、核分裂のあるテーブルは燃焼度の軸を持つこと（燃焼度 0 は軸のない新しい燃料のテーブルでも
//! よい）。核分裂のないテーブル（反射体）は燃焼しないので、軸なしにどの燃焼度でも使う（ADR-0187）。

use xlc_format::fewgroup::{FewGroup, Pins};
use xlc_format::table::{Axis, ROD_OUT, Table};

/// テーブルを内挿する炉心の状態（ケースの状態 `input::CaseState` のホウ素濃度を決めたもの。ADR-0145、ADR-0146）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct State {
    /// ホウ素濃度 (ppm)
    pub boron: f64,
    /// ホウ素の中の B-10 の原子分率（テーブルの値と同じであること）
    pub b10: f64,
    /// 減速材の温度 (K)。密度はテーブルの格子点の組で決まる
    pub moderator: f64,
    /// 燃料の温度 (K)
    pub fuel: f64,
    /// 燃焼度 (MWd/kg、初めの重金属)。燃焼していない炉心は `None`（ADR-0187）
    pub burnup: Option<f64>,
}

/// 内挿する軸の、格子点の番号と重み（番号の昇順。重みの和は 1）。
type Terms = Vec<(usize, f64)>;

/// `points` の中で `x` を挟む2点（`transform` をかけた値で線形）。`x` が格子点なら、その点の重みが厳密に 1
fn bracket(
    name: &str,
    points: &[f64],
    x: f64,
    transform: impl Fn(f64) -> f64,
) -> Result<Terms, String> {
    let out = || format!("{name}: state {x} is outside the table points {points:?}");
    if !x.is_finite() {
        return Err(out());
    }
    if let [only] = points {
        return if x == *only {
            Ok(vec![(0, 1.0)])
        } else {
            Err(out())
        };
    }
    let i = points
        .windows(2)
        .position(|w| w[0] <= x && x <= w[1])
        .ok_or_else(out)?;
    let (a, b) = (transform(points[i]), transform(points[i + 1]));
    let weight = (transform(x) - a) / (b - a);
    Ok(vec![(i, 1.0 - weight), (i + 1, weight)])
}

/// 燃焼度の軸の点の重み（ADR-0192）。点の間は 3 次 Hermite で、点 k の傾きは、k と両隣（端では内側の 2 点）の 3 点を
/// 通る放物線の微分。燃焼度 0 の点は Xe と Sm のない新しい燃料で、次の点（平衡に近い）との間の変化は燃焼度でなく時間で
/// 決まるので、傾きに使わず、0 から次の点までは線形にする。傾きに使う点が 2 つなら線形（2 点の放物線はない）。
/// 自然スプラインや PCHIP を使わない理由は ADR-0192（測った誤差）。値は 4 点の線形結合なので、少数群定数の項目ごとに
/// 同じ重みで足せる
fn burnup_terms(points: &[f64], x: f64) -> Result<Terms, String> {
    let linear = bracket("burnup", points, x, |x| x)?;
    let first = usize::from(points[0] == 0.0);
    let [(i, _), (_, t)] = linear[..] else {
        return Ok(linear);
    };
    if i < first || points.len() - first < 3 {
        return Ok(linear);
    }
    let smooth = &points[first..];
    let j = i - first;
    let h = smooth[j + 1] - smooth[j];
    // Hermite の基底（t = 0、1 で厳密に 0 か 1 になる積の形）
    let h00 = (1.0 + 2.0 * t) * (1.0 - t) * (1.0 - t);
    let h01 = t * t * (3.0 - 2.0 * t);
    let h10 = t * (1.0 - t) * (1.0 - t);
    let h11 = -t * t * (1.0 - t);
    let mut weights = [0.0; 4];
    // weights[m] は smooth[base + m] の重み。base は区間の左の外側の点（j = 0 なら区間の左の点で、m = 3 を使わない）。
    // slope は base..base + 4 の中の点しか返さない（窓は区間の両端の点を含む 3 点）
    let base = j.saturating_sub(1);
    weights[j - base] += h00;
    weights[j + 1 - base] += h01;
    for (k, hk) in [(j, h10), (j + 1, h11)] {
        for (p, c) in slope(smooth, k) {
            weights[p - base] += h * hk * c;
        }
    }
    Ok(weights
        .iter()
        .enumerate()
        .filter(|&(_, w)| *w != 0.0)
        .map(|(m, &w)| (first + base + m, w))
        .collect())
}

/// 点 `k` の傾きの、3 点の値に掛ける係数（点の番号と係数）。`k` と両隣、端では端の 3 点を通る放物線の微分
fn slope(points: &[f64], k: usize) -> [(usize, f64); 3] {
    let s = k.saturating_sub(1).min(points.len() - 3);
    let p = [points[s], points[s + 1], points[s + 2]];
    let x = points[k];
    // 3 点のラグランジュの基底の微分: l_0'(x) = (2x − p_1 − p_2) / ((p_0 − p_1)(p_0 − p_2))、ほかも添字を入れ替えた形
    let coefficient = |m: usize| {
        let (a, b) = match m {
            0 => (p[1], p[2]),
            1 => (p[0], p[2]),
            _ => (p[0], p[1]),
        };
        (2.0 * x - a - b) / ((p[m] - a) * (p[m] - b))
    };
    [
        (s, coefficient(0)),
        (s + 1, coefficient(1)),
        (s + 2, coefficient(2)),
    ]
}

/// `w` 倍した少数群定数（群境界、f* の印、z の厚さはそのまま）。項目を足したらここで扱いを決める（`..` を使わない）
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
        height,
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
        height: *height,
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

/// `acc += w * f`。項目を足したらここで扱いを決める（`..` を使わない）。真空の辺、z の厚さ、形状関数の有無と
/// セルはどの状態でも同じはず（テーブルの検査）で、違えばエラー
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
        height,
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
    if acc.height != *height {
        return Err("the table states differ in their height".into());
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
/// 定数。`state` が `None` なら、軸のないテーブルの状態をそのまま返す。`rod` は制御棒の区間の材料の名前で、`None`
/// は棒のない区間（制御棒の軸があれば `out` を選ぶ）。棒のある区間は、その値を持つ制御棒の軸が要る。
pub fn interpolate(
    table: &Table,
    state: Option<&State>,
    rod: Option<&str>,
) -> Result<Vec<FewGroup>, String> {
    let sub = table.sub_nodes();
    if let Some(rod) = rod
        && !table.axes.iter().any(|a| matches!(a, Axis::Rod { .. }))
    {
        return Err(format!(
            "the rod value {rod:?} is in a layer whose table has no rod axis"
        ));
    }
    let Some(state) = state else {
        if !table.axes.is_empty() {
            return Err(
                "the table has branch axes; give the case named states [state.<name>]".into(),
            );
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
    let has_burnup = table.axes.iter().any(|a| matches!(a, Axis::Burnup { .. }));
    // 燃焼度 0 は新しい燃料のテーブルの状態そのもの
    if let Some(burnup) = state.burnup
        && burnup != 0.0
        && !has_burnup
        && table
            .states
            .iter()
            .any(|s| s.nu_fission.iter().any(|&v| v > 0.0))
    {
        return Err(format!(
            "the table has fission but no burnup axis for the burnup {burnup} of the state (ADR-0187)"
        ));
    }
    // テーブルの軸の順に、内挿に使う格子点と重み
    let mut terms = Vec::with_capacity(table.axes.len());
    for axis in &table.axes {
        terms.push(match axis {
            Axis::Burnup { points } => {
                let burnup = state.burnup.ok_or(
                    "the table has a burnup axis; give the state a burnup (ADR-0187)",
                )?;
                burnup_terms(points, burnup)?
            }
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
            Axis::Rod { values } => {
                let wanted = rod.unwrap_or(ROD_OUT);
                let i = values.iter().position(|v| v == wanted).ok_or_else(|| {
                    format!("rod: the table has the values {values:?}, not {wanted:?}")
                })?;
                vec![(i, 1.0)]
            }
        });
    }

    let mut out: Vec<FewGroup> = (0..sub)
        .map(|q| table.states.get(q).map(|s| scaled(s, 0.0)))
        .collect::<Option<_>>()
        .ok_or("the table has no states")?;
    // 軸ごとの格子点の直積の重みつきの和（最初の軸がいちばん遅く変わる並び）。部分ノードごとに足す
    let mut pick = vec![0usize; terms.len()];
    'corners: loop {
        let mut index = 0;
        let mut weight = 1.0;
        for ((axis, t), &p) in table.axes.iter().zip(&terms).zip(&pick) {
            let (i, w) = t[p];
            index = index * axis.len() + i;
            weight *= w;
        }
        if weight != 0.0 {
            for (q, acc) in out.iter_mut().enumerate() {
                let s = table
                    .states
                    .get(index * sub + q)
                    .ok_or_else(|| format!("the table has no state {index} sub-node {q}"))?;
                accumulate(acc, weight, s)?;
            }
        }
        // 最後の軸から桁を進める
        let mut k = terms.len();
        loop {
            if k == 0 {
                break 'corners;
            }
            k -= 1;
            pick[k] += 1;
            if pick[k] < terms[k].len() {
                break;
            }
            pick[k] = 0;
        }
    }
    // 燃焼度の 3 次 Hermite は凸結合でない（重みに負がある）ので、ソルバに渡す前に値を確かめる。負の重みで内挿したときは、
    // 値の検査（validate）に加えて断面積と拡散係数が負でないことを確かめ、どちらのエラーにも燃焼度の点を足す示唆を添える
    // （折れ曲がりの近くの 0 に近い量が負になりうる。ADR-0192）
    let convex = terms.iter().flatten().all(|&(_, w)| w >= 0.0);
    for f in &out {
        let hint = |e: String| {
            if convex {
                e
            } else {
                format!(
                    "burnup {}: {e}; the cubic burnup interpolation overshot, add burnup points near the bend (ADR-0192)",
                    state.burnup.unwrap_or(0.0)
                )
            }
        };
        f.validate().map_err(|e| hint(e.to_string()))?;
        if !convex {
            non_negative(f).map_err(hint)?;
        }
    }
    Ok(out)
}

/// 断面積と拡散係数が負でないこと（`FewGroup::validate` は有限かどうかだけを見る）。項目を足したらここで扱いを決める
/// （`..` を使わない）
fn non_negative(f: &FewGroup) -> Result<(), String> {
    let FewGroup {
        bounds: _,
        k_inf: _,
        total,
        absorption,
        nu_fission,
        fission,
        chi,
        diffusion,
        scatter,
        // 不連続因子、真空の辺の比、形状関数は validate が範囲を確かめる
        adf: _,
        vacuum_ratio: _,
        normalized_adf: _,
        height: _,
        pins: _,
    } = f;
    let vectors = [
        ("total", total),
        ("absorption", absorption),
        ("nu_fission", nu_fission),
        ("fission", fission),
        ("chi", chi),
        ("diffusion", diffusion),
    ];
    for (name, v) in vectors
        .into_iter()
        .chain(scatter.iter().map(|row| ("scatter", row)))
    {
        if let Some(bad) = v.iter().find(|&&x| x < 0.0) {
            return Err(format!("the interpolation gave a negative {name} {bad}"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
