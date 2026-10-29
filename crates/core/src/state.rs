//! 炉心の状態からの分岐のテーブルの内挿（ADR-0097、ADR-0200）。
//!
//! 軸ごとに隣の2点で、軸の内挿の変数（テーブルの軸の `transform`。ADR-0097 の ppm、温度、√T）について線形に内挿する。格子点の範囲の外はエラーにし、
//! 外挿しない。点が1つの軸は、状態がその値に等しいときだけ使える。状態を与えるときはテーブルが3つの軸をすべて
//! 持ち、与えないときは軸のないテーブルだけを使う。制御棒の離散の軸は内挿せず、値の名前で選ぶ（棒がなければ
//! `out`。ADR-0127）。燃焼度は、各点の傾きをその点と両隣の 3 点の放物線の微分にした 3 次 Hermite で内挿する
//! （ADR-0192）。燃焼度 0 の点は傾きに使わず、0 から次の点までは線形にする。傾きに使う点が 2 つの軸は線形になる。
//! 燃焼度を与える状態では、核分裂のあるテーブルは燃焼度の軸を持つこと（燃焼度 0 は軸のない新しい燃料のテーブルでも
//! よい）。核分裂のないテーブル（反射体）は燃焼しないので、軸なしにどの燃焼度でも使う（ADR-0187）。
//!
//! テーブルが軸の直積でなく、公称の状態と成分（公称と違う軸の組）の部分格子の点を持つときは、制御棒の値の点の中で cut-HDMR の式
//! f(x) = Σ_{u∈U} Σ_{v⊆u} (−1)^{|u|−|v|} f(x_v, x₀) で組み立てる（U は成分の集合、f(x_v, x₀) は v の軸だけ x の値にし、ほかを公称に
//! した状態を、その部分格子の中で多重線形に内挿した値。ADR-0200）。f(x_v, x₀) の係数をまとめて、成分 v の部分格子に
//! c_v = Σ_{u∈U, u⊇v} (−1)^{|u|−|v|} を掛けて足す。軸の直積は c が 3 つの軸の成分だけ 1 で、ほかの成分は 0 なので、直積の多重線形
//! そのものになる。すべての燃焼度の点で解いていない点は、その点の燃焼度の点での公称の状態との差を燃焼度について内挿し、
//! 公称の状態をすべての燃焼度の点で内挿した値に足す（燃焼度の内挿はどちらも ADR-0192）。

use std::collections::BTreeSet;

use xlc_format::fewgroup::{FewGroup, Pins};
use xlc_format::table::{Axis, MicroNuclide, ROD_OUT, Table, component};

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
    check_rod(table, rod)?;
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
    let weights = state_weights(table, state, rod)?;

    let mut out: Vec<FewGroup> = (0..sub)
        .map(|q| table.states.get(q).map(|s| scaled(s, 0.0)))
        .collect::<Option<_>>()
        .ok_or("the table has no states")?;
    // 点の順（軸の直積なら、最初の軸がいちばん遅く変わる順）に、部分ノードごとに足す
    for &(point, weight) in &weights.points {
        if weight != 0.0 {
            for (q, acc) in out.iter_mut().enumerate() {
                let s = table
                    .states
                    .get(point * sub + q)
                    .ok_or_else(|| format!("the table has no state {point} sub-node {q}"))?;
                accumulate(acc, weight, s)?;
            }
        }
    }
    // 燃焼度の 3 次 Hermite と cut-HDMR の組み立ては凸結合でない（重みに負がある）ので、ソルバに渡す前に値を確かめる。負の重みで
    // 組み立てたときは、値の検査（validate）に加えて断面積と拡散係数が負でないことを確かめ、どちらのエラーにも負の重みの出どころと
    // 示唆を添える（折れ曲がりの近くや、解いていない格子点の近くの 0 に近い量が負になりうる。ADR-0192、ADR-0200）
    for f in &out {
        f.validate().map_err(|e| weights.hint(e.to_string()))?;
        if !weights.causes.is_empty() {
            non_negative(f).map_err(|e| weights.hint(e))?;
        }
    }
    Ok(out)
}

/// テーブルの微視的断面積と均質化した数密度を、炉心の状態で内挿した値（ADR-0216）
#[derive(Debug, Clone, PartialEq)]
pub struct MicroAt {
    /// 燃焼度で内挿した、公称の履歴の均質化した数密度 N̄ᵢ⁰ (1/(barn cm)。テーブルの核種の順)
    pub densities: Vec<f64>,
    /// 状態と燃焼度で内挿した微視的断面積（テーブルの核種の順）
    pub nuclides: Vec<MicroNuclide>,
}

/// テーブルの微視的断面積を、少数群定数と同じ重み（状態の組み立てと燃焼度の内挿）で、均質化した数密度を燃焼度の内挿だけで引く。
/// テーブルが `micro` を持たなければエラー
pub fn interpolate_micro(
    table: &Table,
    state: &State,
    rod: Option<&str>,
) -> Result<MicroAt, String> {
    let micro = table
        .micro
        .as_ref()
        .ok_or("the table has no micro data for the micro depletion (ADR-0216)")?;
    check_rod(table, rod)?;
    let weights = state_weights(table, state, rod)?;
    let Some(Axis::Burnup { points }) =
        table.axes.iter().find(|a| matches!(a, Axis::Burnup { .. }))
    else {
        return Err("micro data need a burnup axis (ADR-0216)".into());
    };
    let burnup = state
        .burnup
        .ok_or("the table has a burnup axis; give the state a burnup (ADR-0187)")?;
    let mut densities = vec![0.0; micro.nuclides.len()];
    for (k, h) in burnup_terms(points, burnup)? {
        let row = micro
            .densities
            .get(k)
            .ok_or_else(|| format!("the micro data have no densities at burnup point {k}"))?;
        add(&mut densities, h, row);
    }
    let first = micro
        .states
        .first()
        .ok_or("the micro data have no states")?;
    let mut nuclides: Vec<MicroNuclide> = first.iter().map(|n| scaled_micro(n, 0.0)).collect();
    // 少数群定数と同じ点の順に、同じ重みで足す（重みが 0 の点は使わない）
    for &(point, weight) in &weights.points {
        if weight != 0.0 {
            let state = micro
                .states
                .get(point)
                .ok_or_else(|| format!("the micro data have no state {point}"))?;
            for (acc, n) in nuclides.iter_mut().zip(state) {
                accumulate_micro(acc, weight, n);
            }
        }
    }
    for (name, n) in micro.nuclides.iter().zip(&nuclides) {
        // 収率の重みの割合は 0..1（負の重みの組み立てで外れうる。ADR-0183）
        if let Some(f) = n.thermal_fission
            && !(0.0..=1.0).contains(&f)
        {
            return Err(weights.hint(format!(
                "{name}: the interpolated fraction of fissions below the yield cutoff is {f}"
            )));
        }
        if !weights.causes.is_empty() {
            let negative = [&n.total, &n.nu_fission, &n.fission]
                .into_iter()
                .chain(n.reactions.iter().flat_map(|r| &r.channels))
                .flatten()
                .find(|&&x| x < 0.0);
            if let Some(x) = negative {
                return Err(weights.hint(format!(
                    "{name}: the interpolation gave a negative micro cross section {x}"
                )));
            }
        }
    }
    Ok(MicroAt {
        densities,
        nuclides,
    })
}

/// `w` 倍した微視的断面積。項目を足したらここで扱いを決める（`..` を使わない）
fn scaled_micro(n: &MicroNuclide, w: f64) -> MicroNuclide {
    let MicroNuclide {
        total,
        nu_fission,
        fission,
        scatter,
        reactions,
        thermal_fission,
    } = n;
    let v = |x: &[f64]| x.iter().map(|x| w * x).collect::<Vec<_>>();
    MicroNuclide {
        total: v(total),
        nu_fission: v(nu_fission),
        fission: v(fission),
        scatter: scatter.iter().map(|r| v(r)).collect(),
        reactions: reactions
            .iter()
            .map(|r| xlc_format::table::MicroReaction {
                mt: r.mt,
                channels: r.channels.iter().map(|c| v(c)).collect(),
            })
            .collect(),
        thermal_fission: thermal_fission.map(|f| w * f),
    }
}

/// `acc += w * n`（反応の並びと核分裂の有無はテーブルの検査で状態の間で同じ）。項目を足したらここで扱いを決める
fn accumulate_micro(acc: &mut MicroNuclide, w: f64, n: &MicroNuclide) {
    let MicroNuclide {
        total,
        nu_fission,
        fission,
        scatter,
        reactions,
        thermal_fission,
    } = n;
    add(&mut acc.total, w, total);
    add(&mut acc.nu_fission, w, nu_fission);
    add(&mut acc.fission, w, fission);
    for (a, x) in acc.scatter.iter_mut().zip(scatter) {
        add(a, w, x);
    }
    for (a, x) in acc.reactions.iter_mut().zip(reactions) {
        for (a, x) in a.channels.iter_mut().zip(&x.channels) {
            add(a, w, x);
        }
    }
    if let (Some(a), Some(x)) = (&mut acc.thermal_fission, thermal_fission) {
        *a += w * x;
    }
}

/// 棒のある区間（`rod`）は、制御棒の軸を持つテーブルで引く
fn check_rod(table: &Table, rod: Option<&str>) -> Result<(), String> {
    if let Some(rod) = rod
        && !table.axes.iter().any(|a| matches!(a, Axis::Rod { .. }))
    {
        return Err(format!(
            "the rod value {rod:?} is in a layer whose table has no rod axis"
        ));
    }
    Ok(())
}

/// 内挿に使うテーブルの点の重みと、負の重みの出どころ（負の重みで組み立てたときだけ。値の検査のエラーに添える）
struct Weights {
    /// [点] の重み（点の順。同じ点の重みはまとめる）
    points: Vec<(usize, f64)>,
    causes: Vec<String>,
}

impl Weights {
    /// エラーに、負の重みの出どころと示唆を添える
    fn hint(&self, e: String) -> String {
        if self.causes.is_empty() {
            e
        } else {
            format!("{e}; {}", self.causes.join("; "))
        }
    }
}

/// 炉心の状態 `state` と制御棒の値 `rod` での、テーブルの点の重み（モジュールの説明の組み立て）
fn state_weights(table: &Table, state: &State, rod: Option<&str>) -> Result<Weights, String> {
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
    // テーブルの軸の順に、内挿に使う格子点と重み。範囲の外は、組み立てに使わない軸でもエラーにする（燃焼度の軸は、すべての燃焼度の
    // 点で解いた点の重み）
    let mut terms = Vec::with_capacity(table.axes.len());
    let mut burnup = None;
    for (a, axis) in table.axes.iter().enumerate() {
        terms.push(match axis {
            Axis::Burnup { points } => {
                let value = state.burnup.ok_or(
                    "the table has a burnup axis; give the state a burnup (ADR-0187)",
                )?;
                burnup = Some((a, points.as_slice(), value));
                burnup_terms(points, value)?
            }
            Axis::Boron {
                ppm, b10, transform, ..
            } => {
                if *b10 != state.b10 {
                    return Err(format!(
                        "boron: the table has b10 {b10} and the state {}; they must be exactly equal",
                        state.b10
                    ));
                }
                bracket("boron", ppm, state.boron, |x| transform.apply(x))?
            }
            Axis::Moderator {
                temperature,
                transform,
                ..
            } => bracket("moderator", temperature, state.moderator, |x| {
                transform.apply(x)
            })?,
            Axis::Fuel {
                temperature,
                transform,
                ..
            } => bracket("fuel", temperature, state.fuel, |x| transform.apply(x))?,
            Axis::Rod { values } => {
                let wanted = rod.unwrap_or(ROD_OUT);
                let i = values.iter().position(|v| v == wanted).ok_or_else(|| {
                    format!("rod: the table has the values {values:?}, not {wanted:?}")
                })?;
                vec![(i, 1.0)]
            }
        });
    }
    let mut assembly = Assembly::new(table, &terms, burnup);
    let points = assembly.weights()?;

    // 負の重みの出どころと示唆（ADR-0192、ADR-0200）
    let mut causes = Vec::new();
    if assembly.negative_burnup {
        causes.push(format!(
            "the cubic burnup interpolation at burnup {} overshot, add burnup points near the bend (ADR-0192)",
            state.burnup.unwrap_or(0.0)
        ));
    }
    if assembly.negative_difference {
        causes.push(
            "a state solved at some burnups adds its difference from the nominal state, add its burnup points \
             (ADR-0200)"
                .to_string(),
        );
    }
    if assembly.negative_components {
        causes.push(
            "the cut-HDMR assembly of the subgrids went past the solved states, add the subgrid of more axes \
             (ADR-0200)"
                .to_string(),
        );
    }
    Ok(Weights { points, causes })
}

/// テーブルの点の、状態の値に掛ける重みの組み立て（モジュールの説明の cut-HDMR と、燃焼度の公称との差の内挿）
struct Assembly<'a> {
    table: &'a Table,
    /// [軸] の格子点と重み（燃焼度の軸は、すべての燃焼度の点で解いた点の重み）
    terms: &'a [Terms],
    /// 燃焼度の軸の番号、点、状態の燃焼度
    burnup: Option<(usize, &'a [f64], f64)>,
    /// 部分格子で組み立てる軸（ホウ素、減速材、燃料）の番号
    state_axes: Vec<usize>,
    /// 部分格子で組み立てる軸の公称の格子点（テーブルが公称の値を持たなければ `None`）
    nominal: Option<Vec<usize>>,
    /// 制御棒の軸の番号と、選んだ値の番号
    rod: Option<(usize, usize)>,
    /// 負の重みの出どころ
    negative_burnup: bool,
    negative_components: bool,
    /// 公称の状態との差で燃焼度を内挿した点がある（公称の状態を負の重みで引く）
    negative_difference: bool,
}

impl<'a> Assembly<'a> {
    fn new(table: &'a Table, terms: &'a [Terms], burnup: Option<(usize, &'a [f64], f64)>) -> Self {
        let state_axes: Vec<usize> = (0..table.axes.len())
            .filter(|&a| table.axes[a].state_axis().is_some())
            .collect();
        let nominal = state_axes
            .iter()
            .map(|&a| table.axes[a].nominal_index())
            .collect();
        let rod = table
            .axes
            .iter()
            .position(|a| matches!(a, Axis::Rod { .. }))
            .map(|a| (a, terms[a][0].0));
        Self {
            table,
            terms,
            burnup,
            state_axes,
            nominal,
            rod,
            negative_burnup: false,
            negative_components: false,
            negative_difference: false,
        }
    }

    /// 選んだ制御棒の値の点の成分と、その係数 c_v（係数が 0 の成分は除く）。公称の値のないテーブルは直積で、3 つの軸の成分だけ
    fn components(&self) -> Result<Vec<(u32, f64)>, String> {
        let all = (1u32 << self.state_axes.len()) - 1;
        let Some(nominal) = &self.nominal else {
            return Ok(vec![(all, 1.0)]);
        };
        let masks: BTreeSet<u32> = self
            .table
            .points
            .iter()
            .filter(|p| self.rod.is_none_or(|(a, r)| p[a] == r))
            .map(|p| {
                let at: Vec<usize> = self.state_axes.iter().map(|&a| p[a]).collect();
                component(&at, nominal)
            })
            .collect();
        if masks.is_empty() {
            return Err("the table has no states for the rod value".into());
        }
        // c_v = Σ_{u∈U, u⊇v} (−1)^{|u|−|v|}（整数なので丸めない）
        Ok(masks
            .iter()
            .map(|&v| {
                let c: i32 = masks
                    .iter()
                    .filter(|&&u| u & v == v)
                    .map(|&u| if (u ^ v).count_ones() % 2 == 0 { 1 } else { -1 })
                    .sum();
                (v, f64::from(c))
            })
            .filter(|&(_, c)| c != 0.0)
            .collect())
    }

    /// 点の番号（`at` の燃焼度の軸の番号は `k` にする）
    fn point(&self, at: &[usize], k: Option<usize>) -> Result<usize, String> {
        let mut at = at.to_vec();
        if let (Some((b, ..)), Some(k)) = (self.burnup, k) {
            at[b] = k;
        }
        self.table
            .point_index(&at)
            .ok_or_else(|| format!("the table has no state at the point {at:?}"))
    }

    /// [点] の重み（点の順。同じ点の重みはまとめる）
    fn weights(&mut self) -> Result<Vec<(usize, f64)>, String> {
        let mut weights: Vec<(usize, f64)> = Vec::new();
        self.negative_burnup = self
            .burnup
            .is_some_and(|(b, ..)| self.terms[b].iter().any(|&(_, w)| w < 0.0));
        for (v, c) in self.components()? {
            if c < 0.0 {
                self.negative_components = true;
            }
            // 軸ごとの格子点と重み: 成分の軸は状態の重み、成分でない状態の軸は公称の点、燃焼度の軸はあとで点ごとに
            let corner_terms: Vec<Terms> = (0..self.table.axes.len())
                .map(|a| match self.state_axes.iter().position(|&s| s == a) {
                    Some(i) if v & (1 << i) == 0 => {
                        let nominal = self.nominal.as_ref().map_or(0, |n| n[i]);
                        vec![(nominal, 1.0)]
                    }
                    _ if self.burnup.is_some_and(|(b, ..)| b == a) => vec![(0, 1.0)],
                    _ => self.terms[a].clone(),
                })
                .collect();
            let mut pick = vec![0usize; corner_terms.len()];
            'corners: loop {
                let at: Vec<usize> = corner_terms
                    .iter()
                    .zip(&pick)
                    .map(|(t, &p)| t[p].0)
                    .collect();
                self.corner(&at, &corner_terms, &pick, c, &mut weights)?;
                // 最後の軸から桁を進める
                let mut k = corner_terms.len();
                loop {
                    if k == 0 {
                        break 'corners;
                    }
                    k -= 1;
                    pick[k] += 1;
                    if pick[k] < corner_terms[k].len() {
                        break;
                    }
                    pick[k] = 0;
                }
            }
        }
        weights.sort_by_key(|&(p, _)| p);
        let mut merged: Vec<(usize, f64)> = Vec::with_capacity(weights.len());
        for (p, w) in weights {
            match merged.last_mut() {
                Some((last, sum)) if *last == p => *sum += w,
                _ => merged.push((p, w)),
            }
        }
        Ok(merged)
    }

    /// 成分の係数 `c` の部分格子の角 `at`（`corner_terms` の `pick` の格子点）の重みを `weights` に足す。燃焼度の軸があれば、
    /// その点を解いた燃焼度の点で
    fn corner(
        &mut self,
        at: &[usize],
        corner_terms: &[Terms],
        pick: &[usize],
        c: f64,
        weights: &mut Vec<(usize, f64)>,
    ) -> Result<(), String> {
        // 軸の順に掛けた重み（燃焼度の軸は `burnup_weight`）
        let product = |burnup_weight: f64| {
            let mut weight = 1.0;
            for (a, (t, &p)) in corner_terms.iter().zip(pick).enumerate() {
                weight *= if self.burnup.is_some_and(|(b, ..)| b == a) {
                    burnup_weight
                } else {
                    t[p].1
                };
            }
            c * weight
        };
        // 重みが 0 の角は使わない（解いた燃焼度の点の範囲の外でも、エラーにしない）
        if product(1.0) == 0.0 {
            return Ok(());
        }
        let Some((b, points, value)) = self.burnup else {
            weights.push((self.point(at, None)?, product(1.0)));
            return Ok(());
        };
        let solved: Vec<usize> = (0..points.len())
            .filter(|&k| self.point(at, Some(k)).is_ok())
            .collect();
        if solved.len() == points.len() {
            for &(k, h) in &self.terms[b] {
                weights.push((self.point(at, Some(k))?, product(h)));
            }
            return Ok(());
        }
        // 公称の状態との差を、解いた燃焼度の点で内挿する
        let nominal = self.nominal.as_ref().ok_or(
            "points solved at only some burnups need the nominal values of the table (ADR-0200)",
        )?;
        let mut origin = at.to_vec();
        for (&a, &n) in self.state_axes.iter().zip(nominal) {
            origin[a] = n;
        }
        if let (Some((a, _)), Some(Axis::Rod { values })) =
            (self.rod, self.rod.map(|(a, _)| &self.table.axes[a]))
        {
            origin[a] = values
                .iter()
                .position(|v| v == ROD_OUT)
                .ok_or("points solved at only some burnups need the rod value out (ADR-0200)")?;
        }
        let at_solved: Vec<f64> = solved.iter().map(|&k| points[k]).collect();
        let difference = burnup_terms(&at_solved, value).map_err(|e| {
            format!("the point {at:?} is solved at the burnups {at_solved:?} only: {e}")
        })?;
        let weight = product(1.0);
        self.negative_difference = true;
        for &(k, h) in &self.terms[b] {
            weights.push((self.point(&origin, Some(k))?, weight * h));
        }
        for (i, h) in difference {
            if h < 0.0 {
                self.negative_burnup = true;
            }
            let k = solved[i];
            weights.push((self.point(at, Some(k))?, weight * h));
            weights.push((self.point(&origin, Some(k))?, -weight * h));
        }
        Ok(())
    }
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
