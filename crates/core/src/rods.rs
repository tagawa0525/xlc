//! 制御棒の先端や棒の中の材料の境界がノードの途中にあるノード（部分挿入のノード）の定数（ADR-0127）。
//!
//! 部分挿入のノードは、制御棒の状態が一様な区間（下から）に分かれ、区間ごとに分岐のテーブルの定数を持つ。
//! 同じ列の、部分挿入のノードが続く範囲にその上と下の一様なノードを 1 つずつ加えた範囲（窓）で、横方向に積分した
//! 多群の拡散の式を細かいメッシュの差分で解き、区間の中性子束の積分を重みに定数を混ぜる。z の 2 つの面には、
//! 細かいメッシュの不均質な面の中性子束と、混ぜた定数で NEM の z の向きの式を解いた面の中性子束の比を、不連続因子として
//! 付ける。外側の反復のたびに、そのときの k、ノード平均、横方向の漏れで解き直す（`nem::Refresh`）。
//!
//! 細かいメッシュの差分の面では、区間の境界（ノードの中）は中性子束を連続に、ノードの面は GET のとおり
//! 不連続因子をかけた中性子束 f φ を連続にする（f は、その面に接する区間の定数の z の辺の ADF）。窓の端は、一様な
//! ノードを加えた側はそのノードの細かいメッシュの平均がノード平均（そのときの反復の値）に等しいことを条件にし、
//! 炉心の端までの側は core と同じ境界条件にする。

use std::ops::Range;

use xlc_format::fewgroup::FewGroup;
use xlc_numeric::lu::Lu;
use xlc_numeric::nem::{Inflow, Line, is_fissile, outflow};

use crate::nem::{
    Boundary, NemError, NodeSnapshot, Options, Problem, Refresh, Solution, solve, solve_with,
};
use crate::state::MicroAt;
use xlc_format::table::{MicroNuclide, MicroReaction};

/// 部分挿入のノードの、制御棒の状態が一様な区間。
#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    /// 区間の定数（層のテーブルを状態と区間の制御棒の値で内挿し、ノードの向きに回したもの）
    pub constants: FewGroup,
    /// z の長さ (cm)
    pub length: f64,
}

/// 部分挿入のノード。
#[derive(Debug, Clone, PartialEq)]
pub struct PartialNode {
    /// 直方体の位置（`Problem::node_material` の並び）
    pub site: usize,
    /// このノードだけが使う物質の番号。混ぜた定数で置き換える
    pub material: usize,
    /// 下（−z）からの区間。長さの和はノードの z の幅
    pub segments: Vec<Segment>,
}

/// 区間の定数の混ぜ方。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Weighting {
    /// 軸方向の 1D の細かいメッシュの解の中性子束で重みづけ、z の面に不連続因子を付ける（ADR-0127）
    Flux,
    /// 区間の長さで重みづける（ADR-0127 で採らなかった案。比べるためだけに使う）
    Volume,
}

/// 窓: 同じ列のノードの、z の位置の範囲と、[下、上] の端に一様なノードを加えたか（加えなければ炉心の端まで）。
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Window {
    pub nodes: Range<usize>,
    pub constrained: [bool; 2],
}

/// 列の z の位置ごとの部分挿入のノードの印から、窓を作る（ADR-0127 の「窓」）。部分挿入のノードが続く範囲に、
/// 上と下の一様なノードを 1 つずつ加え、同じ一様なノードを端に持つ 2 つの窓は 1 つにまとめる（まとめると両端が炉心の
/// 端になるときは、まとめない）。どちらの端にも一様なノードを加えられない窓は、細かいメッシュの解の大きさが決まらない
/// のでエラー。
pub(crate) fn windows(partial: &[bool]) -> Result<Vec<Window>, String> {
    let n = partial.len();
    let mut out: Vec<Window> = Vec::new();
    let mut z = 0;
    while z < n {
        if !partial[z] {
            z += 1;
            continue;
        }
        let start = z;
        while z < n && partial[z] {
            z += 1;
        }
        let window = Window {
            nodes: start.saturating_sub(1)..(z + 1).min(n),
            constrained: [start > 0, z < n],
        };
        match out.last_mut() {
            // 前の窓の上端の一様なノードが、この窓の下端の一様なノードと同じ。まとめると両端が炉心の端になるときは、
            // まとめずに、共有するノードでそれぞれを拘束する（同じノードが 2 つの窓に入る）
            Some(last)
                if window.nodes.start < last.nodes.end
                    && (last.constrained[0] || window.constrained[1]) =>
            {
                last.nodes.end = window.nodes.end;
                last.constrained[1] = window.constrained[1];
            }
            _ => out.push(window),
        }
    }
    if let Some(w) = out.iter().find(|w| w.constrained == [false, false]) {
        return Err(format!(
            "the partially rodded nodes at z positions {:?} have no uniform node above or below \
             in their column, so the axial fine-mesh flux has no scale",
            w.nodes
        ));
    }
    Ok(out)
}

/// 部分挿入のノードの [区間][群] の、混ぜた定数を作った重み（中性子束の重みでは区間の細かいメッシュの中性子束の積分、体積の
/// 重みでは区間の長さ）。`partial` の順
pub type SegmentWeights = Vec<Vec<Vec<f64>>>;

/// 部分挿入のノードを `weighting` で扱って炉心を解く。`fine` は細かいメッシュの幅の上限 (cm)。解と、混ぜた定数を
/// 入れた物質の並びと、その定数を混ぜた重みを返す。
pub fn solve_rodded(
    problem: &Problem,
    partial: Vec<PartialNode>,
    weighting: Weighting,
    fine: f64,
    options: &Options,
) -> Result<(Solution, Vec<FewGroup>, SegmentWeights), NemError> {
    let invalid = NemError::InvalidProblem;
    check_partial(problem, &partial).map_err(invalid)?;
    let mut initial = problem.clone();
    let lengths = lengths(&partial);
    for (p, lengths) in partial.iter().zip(&lengths) {
        let mut mixed = mix(&p.segments, lengths).map_err(invalid)?;
        // 体積の重みは z の辺の ADF も区間の長さで混ぜる。中性子束の重みでは、最初の反復の前の値
        for side in [4, 5] {
            for g in 0..mixed.groups() {
                mixed.adf[side][g] = weighted(&p.segments, lengths, g, |c| c.adf[side][g]);
            }
        }
        initial.materials[p.material] = mixed;
    }
    match weighting {
        Weighting::Volume => {
            let solution = solve(&initial, options)?;
            Ok((solution, initial.materials, lengths))
        }
        Weighting::Flux => {
            if !(fine > 0.0 && fine.is_finite()) {
                return Err(invalid(format!(
                    "the axial fine mesh must be positive, got {fine} cm"
                )));
            }
            let mut cusping = Cusping::new(problem, &partial, fine).map_err(invalid)?;
            let (solution, materials) = solve_with(&initial, options, &mut cusping)?;
            Ok((solution, materials, cusping.weights))
        }
    }
}

/// 部分挿入のノードの [区間][群] の長さの重み
fn lengths(partial: &[PartialNode]) -> SegmentWeights {
    partial
        .iter()
        .map(|p| {
            p.segments
                .iter()
                .map(|s| vec![s.length; s.constants.groups()])
                .collect()
        })
        .collect()
}

/// 部分挿入のノードの入力を確かめる: 位置にノードがあり、物質がそのノードだけのもので、区間が 2 つ以上、長さが正で
/// 和がノードの幅、区間の定数が核分裂を持ち（核分裂のない層の部分挿入は、軸方向の反射体の層の制御棒の軸（#0030）
/// とともに扱う）、z の辺に f* がなく、群、f* の辺、z の厚さ、真空の辺、ピンのセルが区間の間で同じ。
fn check_partial(problem: &Problem, partial: &[PartialNode]) -> Result<(), String> {
    // 区間の長さは層の境界と制御棒の位置の差の和なので、丸めを許す
    const TOLERANCE: f64 = 1.0e-9;
    let [nx, ny, _] = problem.mesh.counts();
    for (i, p) in partial.iter().enumerate() {
        if partial[..i].iter().any(|q| q.site == p.site) {
            return Err(format!(
                "position {} has two partially rodded nodes",
                p.site
            ));
        }
        let at = format!("partially rodded node at position {}", p.site);
        if problem.node_material.get(p.site).copied().flatten() != Some(p.material) {
            return Err(format!(
                "{at}: the node does not use material {}",
                p.material
            ));
        }
        if problem
            .node_material
            .iter()
            .enumerate()
            .any(|(site, m)| site != p.site && *m == Some(p.material))
        {
            return Err(format!(
                "{at}: material {} is used by other nodes too",
                p.material
            ));
        }
        let dz = problem.mesh.dz[p.site / (nx * ny)];
        let total: f64 = p.segments.iter().map(|s| s.length).sum();
        if p.segments.len() < 2
            || p.segments
                .iter()
                .any(|s| !(s.length > 0.0 && s.length.is_finite()))
            || (total - dz).abs() > TOLERANCE * dz
        {
            return Err(format!(
                "{at}: needs two or more intervals of positive lengths summing to the node height \
                 {dz} cm"
            ));
        }
        let first = &p.segments[0].constants;
        // 混ぜた定数の z の辺は細かいメッシュの不連続因子で置き換えるので、燃料の ADF で割った f* の印を残せない
        if first.normalized_adf[4] || first.normalized_adf[5] {
            return Err(format!(
                "{at}: the intervals have f* on a z side, which the mixed z-face factors replace"
            ));
        }
        let vacuum = |m: &FewGroup| {
            m.vacuum_ratio
                .iter()
                .map(Option::is_some)
                .collect::<Vec<_>>()
        };
        let cells = |m: &FewGroup| {
            m.pins.as_ref().map(|p| {
                (
                    p.dx.clone(),
                    p.dy.clone(),
                    p.normalized_corners.clone(),
                    p.halves.clone(),
                )
            })
        };
        for s in &p.segments {
            let c = &s.constants;
            c.validate().map_err(|e| format!("{at}: {e}"))?;
            if !is_fissile(&c.nu_fission).map_err(|e| format!("{at}: {e}"))? {
                return Err(format!(
                    "{at}: an interval has no fission; partial rods in non-fissile layers are not \
                     handled yet (#0030)"
                ));
            }
            if c.bounds != first.bounds
                || c.normalized_adf != first.normalized_adf
                || c.height != first.height
                || vacuum(c) != vacuum(first)
                || cells(c) != cells(first)
            {
                return Err(format!(
                    "{at}: the intervals differ in the groups, the f* sides, the height, the vacuum \
                     sides or the pin cells"
                ));
            }
        }
    }
    Ok(())
}

/// [区間][群] の重みが正で有限なこと（`mix` と `mix_micro` で同じ）
fn check_weights(weights: &[Vec<f64>]) -> Result<(), String> {
    if weights
        .iter()
        .flatten()
        .any(|&w| !(w > 0.0 && w.is_finite()))
    {
        return Err(format!(
            "the interval weights of a partially rodded node must be positive, got {weights:?}"
        ));
    }
    Ok(())
}

/// 区間ごとの量 `value` の、[区間][群] の重み `weights` の群 g の重みつき平均。少数群定数（`mix`）と微視的断面積
/// （`mix_micro`）の両方をこれで混ぜる（ADR-0217）
fn average<T>(items: &[T], weights: &[Vec<f64>], g: usize, value: impl Fn(&T) -> f64) -> f64 {
    let total: f64 = weights.iter().map(|w| w[g]).sum();
    items
        .iter()
        .zip(weights)
        .map(|(item, w)| w[g] * value(item))
        .sum::<f64>()
        / total
}

/// 区間の定数の量 `value` の、[区間][群] の重み `weights` の群 g の重みつき平均
fn weighted(
    segments: &[Segment],
    weights: &[Vec<f64>],
    g: usize,
    value: impl Fn(&FewGroup) -> f64,
) -> f64 {
    average(segments, weights, g, |s| value(&s.constants))
}

/// 区間の定数を、[区間][群] の中性子束の積分 `flux` を重みに混ぜる（ADR-0127 の「混ぜる量」）。反応の断面積、
/// 拡散係数、径方向の ADF と真空の辺の比、角の不連続因子は群の中性子束で、群の間の散乱は入射の群の中性子束で、χ は
/// 核分裂の中性子源で、形状関数は群の核分裂率で混ぜる。z の辺の ADF は呼ぶ側が決める（ここでは最初の区間の値）。
/// z の真空の辺の比は、その面に接する区間の値。k∞ は core の使わない量で、全群の中性子束で混ぜる。
fn mix(segments: &[Segment], flux: &[Vec<f64>]) -> Result<FewGroup, String> {
    let first = &segments
        .first()
        .ok_or("a partially rodded node has no intervals")?
        .constants;
    let last = &segments[segments.len() - 1].constants;
    let ng = first.groups();
    check_weights(flux)?;
    let by_flux = |value: &dyn Fn(&FewGroup, usize) -> f64| -> Vec<f64> {
        (0..ng)
            .map(|g| weighted(segments, flux, g, |c| value(c, g)))
            .collect()
    };
    // 区間の核分裂の中性子源 Σ νΣf Φ（どの群にも同じ値）と、群の核分裂率 Σf Φ
    let source: Vec<Vec<f64>> = segments
        .iter()
        .zip(flux)
        .map(|(s, f)| {
            let total: f64 = (0..ng).map(|g| s.constants.nu_fission[g] * f[g]).sum();
            vec![total; ng]
        })
        .collect();
    let fission: Vec<Vec<f64>> = segments
        .iter()
        .zip(flux)
        .map(|(s, f)| (0..ng).map(|g| s.constants.fission[g] * f[g]).collect())
        .collect();
    if source.iter().all(|s| s[0] <= 0.0) {
        return Err("a partially rodded node needs fission".into());
    }
    let all: Vec<Vec<f64>> = flux.iter().map(|f| vec![f.iter().sum(); ng]).collect();
    let mut mixed = first.clone();
    mixed.k_inf = weighted(segments, &all, 0, |c| c.k_inf);
    mixed.total = by_flux(&|c, g| c.total[g]);
    mixed.absorption = by_flux(&|c, g| c.absorption[g]);
    mixed.nu_fission = by_flux(&|c, g| c.nu_fission[g]);
    mixed.fission = by_flux(&|c, g| c.fission[g]);
    mixed.diffusion = by_flux(&|c, g| c.diffusion[g]);
    mixed.chi = (0..ng)
        .map(|g| weighted(segments, &source, g, |c| c.chi[g]))
        .collect();
    mixed.scatter = (0..ng)
        .map(|from| {
            (0..ng)
                .map(|to| weighted(segments, flux, from, |c| c.scatter[from][to]))
                .collect()
        })
        .collect();
    for side in 0..4 {
        mixed.adf[side] = by_flux(&|c, g| c.adf[side][g]);
        if mixed.vacuum_ratio[side].is_some() {
            mixed.vacuum_ratio[side] = Some(by_flux(&|c, g| {
                c.vacuum_ratio[side].as_ref().map_or(f64::NAN, |r| r[g])
            }));
        }
    }
    mixed.vacuum_ratio[5] = last.vacuum_ratio[5].clone();
    if let Some(pins) = &mut mixed.pins {
        // 群の核分裂率が 0 の群の形状関数は、どの区間でも 0 のはずなので、中性子束で混ぜる
        for (g, rows) in pins.form.iter_mut().enumerate() {
            let weights = if fission.iter().any(|f| f[g] > 0.0) {
                &fission
            } else {
                flux
            };
            for (j, row) in rows.iter_mut().enumerate() {
                for (i, v) in row.iter_mut().enumerate() {
                    *v = weighted(segments, weights, g, |c| {
                        c.pins.as_ref().map_or(f64::NAN, |p| p.form[g][j][i])
                    });
                }
            }
        }
        for (corner, values) in pins.corners.iter_mut().enumerate() {
            for (g, v) in values.iter_mut().enumerate() {
                *v = weighted(segments, flux, g, |c| {
                    c.pins.as_ref().map_or(f64::NAN, |p| p.corners[corner][g])
                });
            }
        }
    }
    mixed.validate().map_err(|e| e.to_string())?;
    Ok(mixed)
}

/// 区間ごとの追う核種の微視的断面積 `micro` を、[区間][群] の重み `weights`（そのノードの少数群定数を混ぜた重み。
/// `solve_rodded` が返す）で混ぜる（ADR-0217）。反応の断面積と経路の断面積は群の重みで、散乱行列は入射の群の重みで
/// （`mix` と同じ）、収率の境界より下の核分裂の割合は、ノード平均の中性子束 `flux` で燃焼が使う核分裂率のうちの区間の分
/// Σ_G σf,s,G (w_s,G / Σ_s w_s,G) φ̄_G で混ぜる。公称の数密度は区間で同じこと
pub(crate) fn mix_micro(
    micro: &[MicroAt],
    weights: &[Vec<f64>],
    flux: &[f64],
) -> Result<MicroAt, String> {
    let first = micro
        .first()
        .ok_or("a partially rodded node has no intervals")?;
    let ng = flux.len();
    if weights.len() != micro.len() || weights.iter().any(|w| w.len() != ng) {
        return Err(format!(
            "{} interval weights for {} intervals of {ng} groups",
            weights.len(),
            micro.len()
        ));
    }
    // 群の数の違う微視的断面積は、群ごとに引く前に止める
    let shaped = |n: &MicroNuclide| {
        [&n.total, &n.nu_fission, &n.fission]
            .into_iter()
            .chain(&n.scatter)
            .chain(n.reactions.iter().flat_map(|r| &r.channels))
            .all(|v| v.len() == ng)
            && n.scatter.len() == ng
    };
    if micro.iter().flat_map(|m| &m.nuclides).any(|n| !shaped(n)) {
        return Err(format!(
            "the interval micro cross sections must have the {ng} groups of the node flux"
        ));
    }
    check_weights(weights)?;
    // 公称の数密度は燃焼度だけの関数で、区間は同じ層のテーブルなので同じ値のはず
    if micro
        .iter()
        .any(|m| m.densities != first.densities || m.nuclides.len() != first.nuclides.len())
    {
        return Err(
            "the intervals of a partially rodded node differ in their nominal micro densities"
                .into(),
        );
    }
    let nuclides = (0..first.nuclides.len())
        .map(|i| {
            let segments: Vec<&MicroNuclide> = micro.iter().map(|m| &m.nuclides[i]).collect();
            let by_group = |value: &dyn Fn(&MicroNuclide) -> &[f64]| -> Vec<f64> {
                (0..ng)
                    .map(|g| average(&segments, weights, g, |n| value(n)[g]))
                    .collect()
            };
            let reference = segments[0];
            if segments.iter().any(|n| {
                n.thermal_fission.is_some() != reference.thermal_fission.is_some()
                    || n.reactions.len() != reference.reactions.len()
                    || n.reactions
                        .iter()
                        .zip(&reference.reactions)
                        .any(|(a, b)| a.mt != b.mt || a.channels.len() != b.channels.len())
            }) {
                return Err(format!(
                    "micro nuclide {i}: the intervals differ in fission or in the chain reactions"
                ));
            }
            let thermal_fission = match reference.thermal_fission {
                None => None,
                Some(_) => {
                    // 割合は率の比なので、燃焼が使う核分裂率 σ̄f,G φ̄_G のうちの区間の分で混ぜると熱の核分裂率の和を保つ。重み w は
                    // 群ごとにしか正規化されていない（細かいメッシュの解と NEM の群の比は一致しない。体積の重みは長さ）ので、群ごとの
                    // 割合 w / Σ_s w に φ̄ を掛ける。割合を持つ区間の σf の群の和は正（`interpolate_micro`）なので、率の和は正のはず
                    let totals: Vec<f64> = (0..ng)
                        .map(|g| weights.iter().map(|w| w[g]).sum())
                        .collect();
                    let (mut rate, mut thermal) = (0.0, 0.0);
                    for (n, w) in segments.iter().zip(weights) {
                        let r: f64 = (0..ng)
                            .map(|g| n.fission[g] * w[g] / totals[g] * flux[g])
                            .sum();
                        rate += r;
                        thermal += r * n.thermal_fission.unwrap_or(f64::NAN);
                    }
                    if rate.is_nan() || rate <= 0.0 {
                        return Err(format!(
                            "micro nuclide {i}: the interval fission rates sum to {rate}"
                        ));
                    }
                    Some(thermal / rate)
                }
            };
            Ok(MicroNuclide {
                total: by_group(&|n| &n.total),
                nu_fission: by_group(&|n| &n.nu_fission),
                fission: by_group(&|n| &n.fission),
                scatter: (0..ng)
                    .map(|from| {
                        (0..ng)
                            .map(|to| average(&segments, weights, from, |n| n.scatter[from][to]))
                            .collect()
                    })
                    .collect(),
                reactions: reference
                    .reactions
                    .iter()
                    .enumerate()
                    .map(|(r, reaction)| MicroReaction {
                        mt: reaction.mt,
                        channels: (0..reaction.channels.len())
                            .map(|c| by_group(&|n| &n.reactions[r].channels[c]))
                            .collect(),
                    })
                    .collect(),
                thermal_fission,
            })
        })
        .collect::<Result<_, String>>()?;
    Ok(MicroAt {
        densities: first.densities.clone(),
        nuclides,
    })
}

/// 窓の中のノード。
struct WindowNode {
    /// 直方体の位置
    site: usize,
    /// z の幅 (cm)
    height: f64,
    /// 部分挿入のノードなら、`Cusping::materials` の番号
    partial: Option<usize>,
    /// 下からの区間（一様なノードは自分の定数の 1 つ）
    segments: Vec<Segment>,
}

/// 細かいメッシュの 1 つのセル。
struct Cell {
    /// 窓の中のノードの番号
    node: usize,
    /// ノードの区間の番号
    segment: usize,
    /// 幅 (cm)
    width: f64,
    /// ノードの中の位置 ξ ∈ [−1/2, 1/2] の [下端, 上端]
    xi: [f64; 2],
}

/// 1 つの窓の 1D の問題。
struct WindowProblem {
    nodes: Vec<WindowNode>,
    cells: Vec<Cell>,
    /// [下、上] の端の炉心の境界条件。一様なノードを加えた端（その平均で拘束する）は `None`
    edges: [Option<Boundary>; 2],
}

/// 部分挿入のノードの定数を、外側の反復のたびに作り直す（`nem::Refresh`）。
struct Cusping {
    groups: usize,
    windows: Vec<WindowProblem>,
    /// [部分挿入のノード] の物質の番号
    materials: Vec<usize>,
    /// 窓のノードの直方体の位置（窓の順、各窓の中は下から）
    sites: Vec<usize>,
    /// [窓] の前の反復の細かいメッシュの [セル][群] の中性子束（核分裂の中性子源に使う。最初はノード平均）
    previous: Vec<Option<Vec<Vec<f64>>>>,
    /// [部分挿入のノード] の、最後に混ぜた定数の重み（[区間][群] の中性子束の積分。最初の更新の前は区間の長さ）
    weights: SegmentWeights,
}

impl Cusping {
    fn new(problem: &Problem, partial: &[PartialNode], fine: f64) -> Result<Self, String> {
        let [nx, ny, nz] = problem.mesh.counts();
        let layer = nx * ny;
        let groups = problem.materials.first().ok_or("no materials")?.groups();
        // 部分挿入のノードのある列ごと
        let mut columns: Vec<usize> = partial.iter().map(|p| p.site % layer).collect();
        columns.sort_unstable();
        columns.dedup();
        let mut windows_out = Vec::new();
        let mut sites = Vec::new();
        for column in columns {
            let at = |z: usize| column + layer * z;
            let index = |z: usize| partial.iter().position(|p| p.site == at(z));
            let marks: Vec<bool> = (0..nz).map(|z| index(z).is_some()).collect();
            for window in windows(&marks)? {
                let nodes = window
                    .nodes
                    .clone()
                    .map(|z| {
                        let site = at(z);
                        let height = problem.mesh.dz[z];
                        Ok(match index(z) {
                            Some(i) => WindowNode {
                                site,
                                height,
                                partial: Some(i),
                                segments: partial[i].segments.clone(),
                            },
                            None => {
                                let m = problem.node_material[site]
                                    .ok_or_else(|| format!("position {site} has no node"))?;
                                WindowNode {
                                    site,
                                    height,
                                    partial: None,
                                    segments: vec![Segment {
                                        constants: problem.materials[m].clone(),
                                        length: height,
                                    }],
                                }
                            }
                        })
                    })
                    .collect::<Result<Vec<_>, String>>()?;
                let edges = [
                    (!window.constrained[0]).then_some(problem.boundaries[4]),
                    (!window.constrained[1]).then_some(problem.boundaries[5]),
                ];
                let mut cells = Vec::new();
                for (n, node) in nodes.iter().enumerate() {
                    let mut bottom = 0.0;
                    for (s, segment) in node.segments.iter().enumerate() {
                        // 区間ごとに等しい幅のセル。長さが幅のちょうど整数倍のときに、丸めでセルが増えないようにする
                        let count = ((segment.length / fine) - 1.0e-9).ceil().max(1.0) as usize;
                        let width = segment.length / count as f64;
                        for c in 0..count {
                            let low = bottom + width * c as f64;
                            cells.push(Cell {
                                node: n,
                                segment: s,
                                width,
                                xi: [low / node.height - 0.5, (low + width) / node.height - 0.5],
                            });
                        }
                        bottom += segment.length;
                    }
                }
                sites.extend(nodes.iter().map(|n| n.site));
                windows_out.push(WindowProblem {
                    nodes,
                    cells,
                    edges,
                });
            }
        }
        Ok(Self {
            groups,
            previous: windows_out.iter().map(|_| None).collect(),
            windows: windows_out,
            materials: partial.iter().map(|p| p.material).collect(),
            sites,
            weights: lengths(partial),
        })
    }
}

impl Refresh for Cusping {
    fn sites(&self) -> &[usize] {
        &self.sites
    }

    fn refresh(
        &mut self,
        k: f64,
        nodes: &[NodeSnapshot],
    ) -> Result<Vec<(usize, FewGroup)>, String> {
        let mut out = Vec::new();
        let mut offset = 0;
        for (window, previous) in self.windows.iter().zip(&mut self.previous) {
            let snapshots = nodes
                .get(offset..offset + window.nodes.len())
                .ok_or("too few node snapshots for the windows")?;
            offset += window.nodes.len();
            let fission_flux = match previous.take() {
                Some(flux) => flux,
                None => window
                    .cells
                    .iter()
                    .map(|cell| snapshots[cell.node].flux.clone())
                    .collect(),
            };
            let solved = window.solve(self.groups, k, snapshots, &fission_flux)?;
            for (n, node) in window.nodes.iter().enumerate() {
                let Some(i) = node.partial else {
                    continue;
                };
                let (mixed, integrals) = solved.constants(window, n, k, &snapshots[n])?;
                out.push((self.materials[i], mixed));
                self.weights[i] = integrals;
            }
            *previous = Some(solved.flux);
        }
        Ok(out)
    }

    /// 部分挿入のノードの z の面で f* が掛けるのは、その面に接する区間（下の面は最初、上の面は最後の区間）の ADF
    fn f_star_partners(&self) -> Vec<(usize, [Vec<f64>; 2])> {
        self.windows
            .iter()
            .flat_map(|w| &w.nodes)
            .filter_map(|node| {
                let i = node.partial?;
                let first = &node.segments.first()?.constants;
                let last = &node.segments.last()?.constants;
                Some((
                    self.materials[i],
                    [first.adf[4].clone(), last.adf[5].clone()],
                ))
            })
            .collect()
    }
}

/// 窓の細かいメッシュの解。
struct FineSolution {
    /// [セル][群] の中性子束
    flux: Vec<Vec<f64>>,
    /// [ノード][下の面、上の面] の、[群] の不均質な面の中性子束 f φ と、+z の向きの正味の中性子流
    faces: Vec<[(Vec<f64>, Vec<f64>); 2]>,
}

/// 区間の定数の、ノードの面の z の辺の ADF（下の面は −z の辺、上の面は +z の辺）
fn face_adf(c: &FewGroup, upper: bool, g: usize) -> f64 {
    c.adf[if upper { 5 } else { 4 }][g]
}

impl WindowProblem {
    fn constants(&self, c: usize) -> &FewGroup {
        let cell = &self.cells[c];
        &self.nodes[cell.node].segments[cell.segment].constants
    }

    /// セル c の 2D/h（半分のセルの差分の係数）
    fn half(&self, c: usize, g: usize) -> f64 {
        2.0 * self.constants(c).diffusion[g] / self.cells[c].width
    }

    /// セル c と c + 1 の面の f。区間の境界（ノードの中）は 1、ノードの面は接する区間の z の辺の ADF。f* の辺
    /// （`normalized_adf`）は、core（`Core::adf`）と同じく f* × 向かいの区間の ADF（ADR-0095。両側が f* でない
    /// ことは core が確かめる）
    fn interface_adfs(&self, c: usize, g: usize) -> [f64; 2] {
        if self.cells[c].node == self.cells[c + 1].node {
            return [1.0, 1.0];
        }
        let (below, above) = (self.constants(c), self.constants(c + 1));
        let (fb, fa) = (face_adf(below, true, g), face_adf(above, false, g));
        [
            if below.normalized_adf[5] { fb * fa } else { fb },
            if above.normalized_adf[4] { fa * fb } else { fa },
        ]
    }

    /// セル c と c + 1 の面の J = α φ_c − β φ_{c+1}（+z の向き）の [α, β] と、不均質な面の中性子束
    /// ψ = p φ_c + q φ_{c+1} の [p, q]。半分のセルの差分 J = a_c(φ_c − ψ/f_c) = a_{c+1}(ψ/f_{c+1} − φ_{c+1})
    /// （a = 2D/h）から ψ = (a_c φ_c + a_{c+1} φ_{c+1}) / (a_c/f_c + a_{c+1}/f_{c+1})
    fn interface(&self, c: usize, g: usize) -> ([f64; 2], [f64; 2]) {
        let (a, b) = (self.half(c, g), self.half(c + 1, g));
        let [fa, fb] = self.interface_adfs(c, g);
        let (ra, rb) = (a / fa, b / fb);
        let sum = ra + rb;
        ([a * rb / sum, b * ra / sum], [a / sum, b / sum])
    }

    /// 炉心の端の面の、外向きの正味の中性子流 J_外 = γ φ_c の γ と、不均質な面の中性子束 ψ = s φ_c の s。
    /// セル c は端のセル、`upper` は上の端か
    fn edge(
        &self,
        boundary: Boundary,
        c: usize,
        upper: bool,
        g: usize,
    ) -> Result<(f64, f64), String> {
        let a = self.half(c, g);
        let constants = self.constants(c);
        let f = face_adf(constants, upper, g);
        Ok(match boundary {
            Boundary::Reflective => (0.0, f),
            Boundary::ZeroFlux => (a, 0.0),
            Boundary::Vacuum | Boundary::Marshak => {
                let r = match boundary {
                    Boundary::Marshak => 0.5,
                    _ => constants.vacuum_ratio[if upper { 5 } else { 4 }]
                        .as_ref()
                        .map(|r| r[g])
                        .ok_or(
                            "a partially rodded node at a vacuum z face needs the vacuum ratio",
                        )?,
                };
                // J_外 = r f φ_面 と J_外 = a(φ_c − φ_面) から φ_面 = a φ_c / (a + r f)
                let rf = r * f;
                (a * rf / (a + rf), f * a / (a + rf))
            }
        })
    }

    /// k、窓のノードの量（ノード平均と横方向の漏れ）で、細かいメッシュの式を解く。核分裂の中性子源は、
    /// [セル][群] の中性子束 `fission_flux`（前の反復の解）から作る固定した源とする。
    ///
    /// 源を固定すると、残る行列は除去と漏れと群の間の散乱だけで、正の解を返す（k を係数にして同時に解く形を採らない
    /// 理由は ADR-0127 の「検討した案」）。収束した解では源は解そのものから作った値に等しい
    fn solve(
        &self,
        ng: usize,
        k: f64,
        snapshots: &[NodeSnapshot],
        fission_flux: &[Vec<f64>],
    ) -> Result<FineSolution, String> {
        let nc = self.cells.len();
        // 拘束された端（一様なノードを加えた端）。その端の外側の面の中性子束が群ごとの未知数で、重ね合わせで決める
        let ends: Vec<usize> = (0..2).filter(|&e| self.edges[e].is_none()).collect();
        let params = ends.len() * ng;
        let columns = 1 + params;
        let param = |e: usize, g: usize| ends.iter().position(|&x| x == e).map(|i| i * ng + g);
        // ブロック三重対角: 対角 [セル] の ng × ng、上下の結合 [セル] の群ごと（対角の行列）、右辺 [セル] の
        // [群][右辺]（0 は横方向の漏れ、1 から端の面の中性子束 1 の応答）
        let mut diagonal = vec![vec![0.0; ng * ng]; nc];
        let mut upper = vec![vec![0.0; ng]; nc];
        let mut lower = vec![vec![0.0; ng]; nc];
        let mut rhs = vec![vec![0.0; ng * columns]; nc];
        for c in 0..nc {
            let cell = &self.cells[c];
            let m = self.constants(c);
            let h = cell.width;
            let transverse = &snapshots[cell.node].transverse;
            let source: f64 = (0..ng)
                .map(|from| m.nu_fission[from] * fission_flux[c][from])
                .sum::<f64>()
                / k;
            for g in 0..ng {
                for from in 0..ng {
                    diagonal[c][g * ng + from] += if from == g {
                        h * (m.total[g] - m.scatter[g][g])
                    } else {
                        -h * m.scatter[from][g]
                    };
                }
                rhs[c][g * columns] += h * m.chi[g] * source;
                // セルの横方向の漏れ: L̄ + l1 ξ + l2 (3ξ² − 1/4) のセルの平均
                let [xa, xb] = cell.xi;
                let [mean, l1, l2] = transverse[g];
                let leakage =
                    mean + l1 * 0.5 * (xa + xb) + l2 * (xa * xa + xa * xb + xb * xb - 0.25);
                rhs[c][g * columns] -= h * leakage;
            }
        }
        for c in 0..nc.saturating_sub(1) {
            for g in 0..ng {
                let ([alpha, beta], _) = self.interface(c, g);
                diagonal[c][g * ng + g] += alpha;
                upper[c][g] = -beta;
                lower[c + 1][g] = -alpha;
                diagonal[c + 1][g * ng + g] += beta;
            }
        }
        for (e, cell) in [0, nc - 1].into_iter().enumerate() {
            for g in 0..ng {
                match (self.edges[e], param(e, g)) {
                    (Some(boundary), _) => {
                        let (gamma, _) = self.edge(boundary, cell, e == 1, g)?;
                        diagonal[cell][g * ng + g] += gamma;
                    }
                    (None, Some(p)) => {
                        // 外側の面の中性子束 φ_面 で、外向きの流れ a(φ_c − φ_面)
                        let a = self.half(cell, g);
                        diagonal[cell][g * ng + g] += a;
                        rhs[cell][g * columns + 1 + p] += a;
                    }
                    (None, None) => return Err("a constrained end has no unknown".into()),
                }
            }
        }
        let x = block_tridiagonal(&diagonal, &upper, &lower, &rhs, ng, columns)?;
        // 拘束: 端の一様なノードの細かいメッシュの平均 = そのノード平均
        let mut theta = vec![0.0; params];
        if params > 0 {
            let mut a = vec![0.0; params * params];
            let mut b = vec![0.0; params];
            for (i, &e) in ends.iter().enumerate() {
                let node = if e == 0 { 0 } else { self.nodes.len() - 1 };
                let height = self.nodes[node].height;
                for g in 0..ng {
                    let row = i * ng + g;
                    let mut constant = 0.0;
                    for (c, cell) in self.cells.iter().enumerate() {
                        if cell.node != node {
                            continue;
                        }
                        let w = cell.width / height;
                        constant += w * x[c][g * columns];
                        for p in 0..params {
                            a[row * params + p] += w * x[c][g * columns + 1 + p];
                        }
                    }
                    b[row] = snapshots[node].flux[g] - constant;
                }
            }
            let lu = Lu::new(a, params)
                .ok_or("the constraints of the axial fine-mesh problem are singular")?;
            lu.solve(&mut b);
            theta = b;
        }
        let flux: Vec<Vec<f64>> = (0..nc)
            .map(|c| {
                (0..ng)
                    .map(|g| {
                        x[c][g * columns]
                            + (0..params)
                                .map(|p| theta[p] * x[c][g * columns + 1 + p])
                                .sum::<f64>()
                    })
                    .collect()
            })
            .collect();
        // ノードの面の不均質な中性子束と正味の中性子流
        let mut faces: Vec<[(Vec<f64>, Vec<f64>); 2]> = self
            .nodes
            .iter()
            .map(|_| std::array::from_fn(|_| (vec![0.0; ng], vec![0.0; ng])))
            .collect();
        for c in 0..nc.saturating_sub(1) {
            let (below, above) = (self.cells[c].node, self.cells[c + 1].node);
            if below == above {
                continue;
            }
            for (g, (&low, &high)) in flux[c].iter().zip(&flux[c + 1]).enumerate() {
                let ([alpha, beta], [p, q]) = self.interface(c, g);
                let current = alpha * low - beta * high;
                let psi = p * low + q * high;
                faces[below][1].0[g] = psi;
                faces[below][1].1[g] = current;
                faces[above][0].0[g] = psi;
                faces[above][0].1[g] = current;
            }
        }
        let last = self.nodes.len() - 1;
        for (e, (cell, node)) in [(0, 0), (nc - 1, last)].into_iter().enumerate() {
            for (g, &phi) in flux[cell].iter().enumerate() {
                let (psi, current) = match (self.edges[e], param(e, g)) {
                    (Some(boundary), _) => {
                        let (gamma, s) = self.edge(boundary, cell, e == 1, g)?;
                        let out = gamma * phi;
                        (s * phi, if e == 0 { -out } else { out })
                    }
                    (None, Some(p)) => {
                        let face = theta[p];
                        let a = self.half(cell, g);
                        let f = face_adf(self.constants(cell), e == 1, g);
                        let current = if e == 0 {
                            a * (face - phi)
                        } else {
                            a * (phi - face)
                        };
                        (f * face, current)
                    }
                    (None, None) => return Err("a constrained end has no unknown".into()),
                };
                faces[node][e].0[g] = psi;
                faces[node][e].1[g] = current;
            }
        }
        Ok(FineSolution { flux, faces })
    }
}

impl FineSolution {
    /// 窓のノード n（部分挿入のノード）の混ぜた定数と z の面の不連続因子、混ぜた重み（[区間][群] の中性子束の積分）
    fn constants(
        &self,
        window: &WindowProblem,
        n: usize,
        k: f64,
        snapshot: &NodeSnapshot,
    ) -> Result<(FewGroup, Vec<Vec<f64>>), String> {
        let node = &window.nodes[n];
        let ng = node.segments[0].constants.groups();
        let mut integrals = vec![vec![0.0; ng]; node.segments.len()];
        for (c, cell) in window.cells.iter().enumerate() {
            if cell.node == n {
                for (i, phi) in integrals[cell.segment].iter_mut().zip(&self.flux[c]) {
                    *i += cell.width * phi;
                }
            }
        }
        let mut mixed = mix(&node.segments, &integrals)?;
        let average: Vec<f64> = (0..ng)
            .map(|g| integrals.iter().map(|i| i[g]).sum::<f64>() / node.height)
            .collect();
        let [(psi_low, j_low), (psi_high, j_high)] = &self.faces[n];
        let homogeneous =
            homogeneous_faces(&mixed, k, node.height, &average, j_low, j_high, snapshot)?;
        let outer = [n == 0, n + 1 == window.nodes.len()];
        for (side, psi) in [psi_low, psi_high].into_iter().enumerate() {
            for g in 0..ng {
                mixed.adf[4 + side][g] = match window.edges[side] {
                    // 炉心の反射と中性子束ゼロの面は、core の面の条件が不連続因子によらない
                    Some(Boundary::Reflective | Boundary::ZeroFlux) if outer[side] => 1.0,
                    _ => {
                        let f = psi[g] / homogeneous[side][g];
                        if !(f > 0.0 && f.is_finite()) {
                            return Err(format!(
                                "node {}: the z-face discontinuity factor is {f} in group {g}",
                                node.site
                            ));
                        }
                        f
                    }
                };
            }
        }
        mixed.validate().map_err(|e| e.to_string())?;
        Ok((mixed, integrals))
    }
}

/// 1ノードの問題の反復の上限と、面の中性子束の（ノード平均に対する）変化の収束の判定（`xlc_numeric::nem` の
/// `surface_fluxes` と同じ）
const ONE_NODE_ITERATIONS: usize = 10_000;
const ONE_NODE_CONVERGED: f64 = 1.0e-14;

/// 混ぜた定数のノードを、core と同じ NEM の z の向きの式（横方向の漏れは `snapshot` の 2 次式）で、ノード平均を
/// `average`、下と上の面の +z の向きの正味の中性子流を `j_low`、`j_high` に固定して解いた、[下、上] の面の中性子束
/// （ADR-0103 の1ノードの問題を z の 1D にしたもの）。
fn homogeneous_faces(
    m: &FewGroup,
    k: f64,
    height: f64,
    average: &[f64],
    j_low: &[f64],
    j_high: &[f64],
    snapshot: &NodeSnapshot,
) -> Result<[Vec<f64>; 2], String> {
    let ng = m.groups();
    let mut coefficients = vec![[0.0; 4]; ng];
    let faces = |a: &[f64; 4], phi: f64| [phi + 0.5 * (a[1] - a[0]), phi + 0.5 * (a[1] + a[0])];
    for _ in 0..ONE_NODE_ITERATIONS {
        let mut change: f64 = 0.0;
        for g in 0..ng {
            let removal = m.total[g] - m.scatter[g][g];
            let coupling: Vec<f64> = (0..ng)
                .map(|from| {
                    let s = if from == g { 0.0 } else { m.scatter[from][g] };
                    s + m.chi[g] * m.nu_fission[from] / k
                })
                .collect();
            let [_, l1, l2] = snapshot.transverse[g];
            // 外向きの正味の中性子流を固定する: i = o − J_外（β = 1）。[+ 面, − 面]
            let line = outflow(&Line {
                diffusion: m.diffusion[g],
                width: height,
                removal,
                coupling: &coupling,
                coefficients: &coefficients,
                group: g,
                transverse: (l1, l2),
                inflow: [
                    Inflow {
                        beta: 1.0,
                        external: -j_high[g],
                    },
                    Inflow {
                        beta: 1.0,
                        external: j_low[g],
                    },
                ],
            });
            let updated = line.at(average[g]).coefficients;
            for (new, old) in faces(&updated, average[g])
                .iter()
                .zip(faces(&coefficients[g], average[g]))
            {
                let relative = (new - old).abs() / average[g];
                // NaN も収束していないとみなす（一度 NaN になれば、あとの比較で消さない）
                if relative.is_nan() || relative > change {
                    change = relative;
                }
            }
            coefficients[g] = updated;
        }
        if change <= ONE_NODE_CONVERGED {
            let mut out = [vec![0.0; ng], vec![0.0; ng]];
            for g in 0..ng {
                let [low, high] = faces(&coefficients[g], average[g]);
                out[0][g] = low;
                out[1][g] = high;
            }
            return Ok(out);
        }
    }
    Err(format!(
        "the one-node z problem of a partially rodded node did not converge in \
         {ONE_NODE_ITERATIONS} iterations"
    ))
}

/// ブロック三重対角の連立一次方程式を、複数の右辺について解く。`diagonal[c]` は行優先の ng × ng、`upper[c]` と
/// `lower[c]` は群ごとの対角の結合（セル c と c + 1、c と c − 1）、`rhs[c]` は [群][右辺] の行優先。
/// 前進: U_c = A_c − L_c M_{c−1}、Y_c = U_c⁻¹ (R_c − L_c Y_{c−1})、M_c = U_c⁻¹ B_c。後退: X_c = Y_c − M_c X_{c+1}
fn block_tridiagonal(
    diagonal: &[Vec<f64>],
    upper: &[Vec<f64>],
    lower: &[Vec<f64>],
    rhs: &[Vec<f64>],
    ng: usize,
    columns: usize,
) -> Result<Vec<Vec<f64>>, String> {
    let nc = diagonal.len();
    let mut m: Vec<Vec<f64>> = Vec::with_capacity(nc);
    let mut y: Vec<Vec<f64>> = Vec::with_capacity(nc);
    for c in 0..nc {
        let mut u = diagonal[c].clone();
        let mut r = rhs[c].clone();
        if let (Some(previous_m), Some(previous_y)) = (
            c.checked_sub(1).map(|p| &m[p]),
            c.checked_sub(1).map(|p| &y[p]),
        ) {
            for g in 0..ng {
                let l = lower[c][g];
                for j in 0..ng {
                    u[g * ng + j] -= l * previous_m[g * ng + j];
                }
                for j in 0..columns {
                    r[g * columns + j] -= l * previous_y[g * columns + j];
                }
            }
        }
        let lu = Lu::new(u, ng).ok_or("the axial fine-mesh problem is singular")?;
        // 列ごとに解く
        let solve_columns = |values: &[f64], width: usize| -> Vec<f64> {
            let mut out = vec![0.0; ng * width];
            for j in 0..width {
                let mut col: Vec<f64> = (0..ng).map(|g| values[g * width + j]).collect();
                lu.solve(&mut col);
                for g in 0..ng {
                    out[g * width + j] = col[g];
                }
            }
            out
        };
        let b: Vec<f64> = (0..ng * ng)
            .map(|i| {
                if i / ng == i % ng {
                    upper[c][i / ng]
                } else {
                    0.0
                }
            })
            .collect();
        m.push(solve_columns(&b, ng));
        y.push(solve_columns(&r, columns));
    }
    let mut x: Vec<Vec<f64>> = vec![Vec::new(); nc];
    for c in (0..nc).rev() {
        let mut v = y[c].clone();
        if c + 1 < nc {
            for g in 0..ng {
                for j in 0..columns {
                    v[g * columns + j] -= (0..ng)
                        .map(|h| m[c][g * ng + h] * x[c + 1][h * columns + j])
                        .sum::<f64>();
                }
            }
        }
        x[c] = v;
    }
    Ok(x)
}

#[cfg(test)]
mod tests;
