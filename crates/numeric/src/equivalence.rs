//! 核分裂のないノードの等価な拡散係数（Koebke の等価理論。ADR-0104）。
//!
//! 反射体のノードの断面積は反応率を保存して中性子束で重みづけするが、拡散係数は、均質な1ノードの解
//! （core と同じ核。核分裂のないノードは ANM）が、ノードの正味の中性子流の向きの両方の面で、参照の正味の
//! 中性子流と、同じ不連続因子をかけた参照の不均質な面の中性子束を再現するように決める（SPH 因子のように、
//! 等価性のための自由度として使う）。拡散係数は群ごとに1つなので、2次元のノードでは、正味の中性子流の大きい
//! 方向（x か y）の2つの面で不連続因子をそろえ、もう一方の方向の不連続因子は GET の1ノードの問題で決める。

use crate::nem::{OneNode, OneNodeError, homogeneous_faces, is_fissile};

/// 根の探索で D を 2 倍ずつ広げる回数の上限（D の 2⁻⁶⁰〜2⁶⁰ 倍）
const BRACKET_STEPS: usize = 60;

/// ln D の区間がこの幅より狭くなったら根とする（D の相対誤差。倍精度の丸めの程度）
const ROOT_WIDTH: f64 = 1.0e-14;

/// 群の間の結合（ほかの群の D を固定して1群ずつ解く Gauss–Seidel）の収束の判定と反復の上限
const CONVERGED: f64 = 1.0e-12;
const MAX_ITERATIONS: usize = 200;

/// 等価な拡散係数の結果。
#[derive(Debug, Clone, PartialEq)]
pub enum Equivalence {
    /// 2つの面の不連続因子をそろえ、面の均質な中性子束を正にする [群] の拡散係数
    Found(Vec<f64>),
    /// 群 `group` で、面の均質な中性子束が正の範囲に、均質な面の比を不均質な比 `target`（t = φ+/φ−）にする D がない。
    /// `range` は、根を探して評価した D のうち、その群の2つの面の均質な中性子束が正の点の、均質な比の最小と最大
    /// （そのような点がなければ `None`）
    Missing {
        group: usize,
        target: f64,
        range: Option<[f64; 2]>,
    },
}

/// 根を挟めないか、収束した根で群 `group` の面の均質な中性子束が正でないときの結果。正の面の点の均質な比の範囲 `range` に不均質な比
/// `target` が入っていれば、正の面の範囲にも根があるのに別の根に収束したので、エラーにする（黙って中性子流の D に
/// しない）。入っていなければ `Missing`
fn negative_root(
    group: usize,
    target: f64,
    range: Option<[f64; 2]>,
) -> Result<Equivalence, OneNodeError> {
    match range {
        Some([low, high]) if low <= target && target <= high => {
            Err(OneNodeError::Invalid(format!(
                "group {group}: the equivalent diffusion coefficient converged to a root with \
             non-positive faces, but positive-face ratios {low}..{high} reach the target {target}"
            )))
        }
        _ => Ok(Equivalence::Missing {
            group,
            target,
            range,
        }),
    }
}

/// 核分裂のないノードの拡散係数をどれにしたか（ADR-0134）。
#[derive(Debug, Clone, PartialEq)]
pub enum Choice {
    /// Koebke の等価な D
    Koebke,
    /// Koebke の D が正の面の範囲になく、中性子流で重みづけた D（`Equivalence::Missing` の量）
    CurrentWeighted {
        group: usize,
        target: f64,
        range: Option<[f64; 2]>,
    },
}

/// 核分裂のないノードの拡散係数（ADR-0134）: Koebke の等価な D（`equivalent_diffusion`）が面の均質な中性子束を正にする
/// 範囲にあればそれを、なければ `node` の拡散係数（中性子流で重みづけた D）を、すべての群で使う。中性子流で重みづけた D
/// で面が正になるかは、呼ぶ側の1ノードの問題（`surface_fluxes`）が確かめる
pub fn reflector_diffusion(
    node: &OneNode,
    heterogeneous: [&[f64]; 4],
) -> Result<(Vec<f64>, Choice), OneNodeError> {
    Ok(match equivalent_diffusion(node, heterogeneous)? {
        Equivalence::Found(d) => (d, Choice::Koebke),
        Equivalence::Missing {
            group,
            target,
            range,
        } => (
            node.diffusion.to_vec(),
            Choice::CurrentWeighted {
                group,
                target,
                range,
            },
        ),
    })
}

/// 核分裂のないノードの、等価な [群] の拡散係数。`node` の拡散係数は初めの値で、ほかは1ノードの問題の入力
/// （ノード平均、面の外向きの正味の流れ、境界）と同じ。`heterogeneous` は [面（-x, +x, -y, +y）][群] の参照の
/// 不均質な面の中性子束。
///
/// 正味の中性子流の大きい方向 d（Σ_g |J+ − J−|、J は面の外向きの正味の流れ）の2つの面で、均質な面の中性子束の比を
/// 不均質な面の中性子束の比 t_g = φ+/φ− にする D_g を、ほかの群の D を固定して ln D について二分法で探し、群の間で
/// 収束するまで繰り返す。源のない群では、cosh(κL) = c、c = (1 + r t)/(r + t)（r は向きをそろえた正味の流れの比）の
/// 閉じた形の解と同じで、c ≤ 1（t が流れの減る向きに増えるなど）なら解がない。根を挟めないか、収束した D で面の均質な
/// 中性子束が正でなければ `Equivalence::Missing` を返す（群の源の形と不均質な面の比が逆向きのときに起きる。ADR-0134）。
/// 正味の中性子流がない、核分裂のあるノード、群の間の反復が収束しないときはエラーにする。
pub fn equivalent_diffusion(
    node: &OneNode,
    heterogeneous: [&[f64]; 4],
) -> Result<Equivalence, OneNodeError> {
    let invalid = |what: &str| Err(OneNodeError::Invalid(what.into()));
    let fissile = is_fissile(node.nu_fission).map_err(|e| OneNodeError::Invalid(e.to_string()))?;
    if fissile {
        return invalid("equivalent diffusion coefficients are for nodes without fission");
    }
    let ng = node.flux.len();
    if heterogeneous
        .iter()
        .any(|face| face.len() != ng || !face.iter().all(|v| *v > 0.0 && v.is_finite()))
    {
        return invalid("heterogeneous face fluxes must be positive for every group");
    }
    let strength = |d: usize| -> f64 {
        (0..ng)
            .map(|g| (node.currents[2 * d + 1][g] - node.currents[2 * d][g]).abs())
            .sum()
    };
    let (sx, sy) = (strength(0), strength(1));
    if !(sx > 0.0 || sy > 0.0) {
        return Err(OneNodeError::NoNetCurrent);
    }
    let d = if sx >= sy { 0 } else { 1 };
    let ratio: Vec<f64> = (0..ng)
        .map(|g| heterogeneous[2 * d + 1][g] / heterogeneous[2 * d][g])
        .collect();
    let mut diffusion = node.diffusion.to_vec();
    // 群 g の残差 φ+ − t φ−（D が根の両側で符号が変わる）
    let residual = |diffusion: &[f64], g: usize| -> Result<f64, OneNodeError> {
        let faces = homogeneous_faces(&OneNode { diffusion, ..*node })?;
        Ok(faces[2 * d + 1][g] - ratio[g] * faces[2 * d][g])
    };
    // 群 g の D を初めの値の 2⁻ⁿ〜2ⁿ 倍（ほかの群は `diffusion`）にしたうち、2つの面の均質な中性子束が正の点の、均質な
    // 比の最小と最大（根がない根拠として返す。ADR-0134）
    let positive_range = |diffusion: &[f64], g: usize| -> Result<Option<[f64; 2]>, OneNodeError> {
        let mut trial = diffusion.to_vec();
        let mut range: Option<[f64; 2]> = None;
        for k in -(BRACKET_STEPS as i32)..=BRACKET_STEPS as i32 {
            trial[g] = node.diffusion[g] * 2f64.powi(k);
            let faces = homogeneous_faces(&OneNode {
                diffusion: &trial,
                ..*node
            })?;
            let (minus, plus) = (faces[2 * d][g], faces[2 * d + 1][g]);
            if minus > 0.0 && plus > 0.0 {
                let r = plus / minus;
                range = Some(range.map_or([r, r], |[lo, hi]| [lo.min(r), hi.max(r)]));
            }
        }
        Ok(range)
    };
    for _ in 0..MAX_ITERATIONS {
        let mut change: f64 = 0.0;
        for g in 0..ng {
            let start = diffusion[g];
            let at = |log_d: f64, diffusion: &mut Vec<f64>| -> Result<f64, OneNodeError> {
                diffusion[g] = log_d.exp();
                residual(diffusion, g)
            };
            let s0 = start.ln();
            let f0 = at(s0, &mut diffusion)?;
            // 根を挟む区間 [lo, hi]（f(lo) と f(hi) の符号が違う）を、D を 2 倍ずつ両側へ広げて探す
            let mut bracket = (f0 == 0.0).then_some((s0, s0, f0));
            for k in 1..=BRACKET_STEPS {
                if bracket.is_some() {
                    break;
                }
                for s in [
                    s0 + k as f64 * std::f64::consts::LN_2,
                    s0 - k as f64 * std::f64::consts::LN_2,
                ] {
                    let f = at(s, &mut diffusion)?;
                    if f == 0.0 || f.signum() != f0.signum() {
                        bracket = Some((s0.min(s), s0.max(s), if s < s0 { f } else { f0 }));
                        break;
                    }
                }
            }
            let Some((mut lo, mut hi, mut f_lo)) = bracket else {
                return negative_root(g, ratio[g], positive_range(&diffusion, g)?);
            };
            while hi - lo > ROOT_WIDTH {
                let mid = 0.5 * (lo + hi);
                let f = at(mid, &mut diffusion)?;
                if f == 0.0 {
                    lo = mid;
                    hi = mid;
                } else if f.signum() == f_lo.signum() {
                    lo = mid;
                    f_lo = f;
                } else {
                    hi = mid;
                }
            }
            diffusion[g] = (0.5 * (lo + hi)).exp();
            change = change.max((diffusion[g] - start).abs() / start);
        }
        if change <= CONVERGED {
            // 面の均質な中性子束が負になる根は、等価な D として使えない
            let faces = homogeneous_faces(&OneNode {
                diffusion: &diffusion,
                ..*node
            })?;
            return match (0..ng).find(|&g| !(faces[2 * d][g] > 0.0 && faces[2 * d + 1][g] > 0.0)) {
                Some(g) => negative_root(g, ratio[g], positive_range(&diffusion, g)?),
                None => Ok(Equivalence::Found(diffusion)),
            };
        }
    }
    Err(OneNodeError::EquivalenceNotConverged(MAX_ITERATIONS))
}

#[cfg(test)]
mod tests;
