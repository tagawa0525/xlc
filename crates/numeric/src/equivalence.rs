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

/// 核分裂のないノードの、等価な [群] の拡散係数。`node` の拡散係数は初めの値で、ほかは1ノードの問題の入力
/// （ノード平均、面の外向きの正味の流れ、境界）と同じ。`heterogeneous` は [面（-x, +x, -y, +y）][群] の参照の
/// 不均質な面の中性子束。
///
/// 正味の中性子流の大きい方向 d（Σ_g |J+ − J−|、J は面の外向きの正味の流れ）の2つの面で、均質な面の中性子束の比を
/// 不均質な面の中性子束の比 t_g = φ+/φ− にする D_g を、ほかの群の D を固定して ln D について二分法で探し、群の間で
/// 収束するまで繰り返す。源のない群では、cosh(κL) = c、c = (1 + r t)/(r + t)（r は向きをそろえた正味の流れの比）の
/// 閉じた形の解と同じで、c ≤ 1（t が流れの減る向きに増えるなど）なら解がない。解がない、正味の中性子流がない、
/// 核分裂のあるノードは、どれもエラーにする（黙って別の D にしない）。
pub fn equivalent_diffusion(
    node: &OneNode,
    heterogeneous: [&[f64]; 4],
) -> Result<Vec<f64>, OneNodeError> {
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
                return Err(OneNodeError::NoEquivalentDiffusion { group: g });
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
            return Ok(diffusion);
        }
    }
    Err(OneNodeError::EquivalenceNotConverged(MAX_ITERATIONS))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nem::{Edge, OneNode, OneNodeError, surface_fluxes};

    const WIDTH: f64 = 21.5;

    /// 群が結合しない2群の、源のない 1D の解 φ = a cosh κx + b sinh κx（x ∈ [0, L]）の
    /// （ノード平均、−x の面の外向きの正味の流れ、+x の面の外向きの正味の流れ、−x の面の中性子束、+x の面の中性子束）
    fn source_free(d: f64, removal: f64, b: f64) -> (f64, f64, f64, f64, f64) {
        let kappa = (removal / d).sqrt();
        let kl = kappa * WIDTH;
        let (ch, sh) = (kl.cosh(), kl.sinh());
        let a = 1.0;
        let left = -d * kappa * b;
        let right = -d * kappa * (a * sh + b * ch);
        let average = (a * sh + b * (ch - 1.0)) / kl;
        (average, -left, right, a, a * ch + b * sh)
    }

    fn node<'a>(
        diffusion: &'a [f64],
        total: &'a [f64],
        scatter: &'a [Vec<f64>],
        flux: &'a [f64],
        currents: [&'a [f64]; 4],
    ) -> OneNode<'a> {
        OneNode {
            diffusion,
            total,
            scatter,
            nu_fission: &[0.0, 0.0],
            chi: &[0.0, 0.0],
            k: 1.0,
            flux,
            widths: [WIDTH, WIDTH],
            currents,
            edges: [
                Edge::Reflective,
                Edge::Vacuum,
                Edge::Reflective,
                Edge::Reflective,
            ],
        }
    }

    #[test]
    fn source_free_groups_take_the_closed_form() {
        // 源のない群では、両方の面の正味の流れ J と面の中性子束の比 t = φ+/φ− を、同じ不連続因子で再現する
        // D は閉じた形 D = Σr (L / acosh c)²、c = (1 + r t)/(r + t)、r = J+/J−（向きをそろえた正味の流れの比）
        let truth = [1.3, 0.35];
        let total = [0.6, 1.7];
        let scatter = vec![vec![0.57, 0.0], vec![0.0, 1.66]];
        let removal = [0.03, 0.04];
        let solutions = [
            source_free(truth[0], removal[0], -0.95),
            source_free(truth[1], removal[1], -0.97),
        ];
        let flux: Vec<f64> = solutions.iter().map(|s| s.0).collect();
        let minus: Vec<f64> = solutions.iter().map(|s| s.1).collect();
        let plus: Vec<f64> = solutions.iter().map(|s| s.2).collect();
        let zero = [0.0, 0.0];
        // 不均質な面の中性子束は、均質な面の中性子束に群ごとの同じ因子を掛けたもの
        let factors = [1.3, 0.7];
        let het_minus: Vec<f64> = (0..2).map(|g| factors[g] * solutions[g].3).collect();
        let het_plus: Vec<f64> = (0..2).map(|g| factors[g] * solutions[g].4).collect();
        let guess = [0.9, 0.25];
        let found = equivalent_diffusion(
            &node(
                &guess,
                &total,
                &scatter,
                &flux,
                [&minus, &plus, &zero, &zero],
            ),
            [&het_minus, &het_plus, &flux, &flux],
        )
        .unwrap();
        for g in 0..2 {
            let r = plus[g] / -minus[g];
            let t = het_plus[g] / het_minus[g];
            let c = (1.0 + r * t) / (r + t);
            let closed = removal[g] * (WIDTH / c.acosh()).powi(2);
            assert!((closed - truth[g]).abs() < 1e-12 * truth[g], "{closed}");
            assert!((found[g] - truth[g]).abs() < 1e-10 * truth[g], "{found:?}");
        }
    }

    #[test]
    fn coupled_groups_recover_the_generating_d() {
        // 下方散乱と上方散乱のある2群。拡散係数 `truth` の均質な解（ANM）の面の中性子束に群ごとの因子を掛けて
        // 不均質な面の中性子束とすると、等価な拡散係数は `truth` に戻る
        let truth = [1.45, 0.35];
        let total = [0.6215, 1.6711];
        let scatter = vec![vec![0.5956, 0.02371], vec![0.000806, 1.6323]];
        let width = WIDTH;
        // −x から入り +x へ出る（外向きの正味の流れ）
        let minus = [-0.01136, -0.0005];
        let plus = [0.000927, 0.0002];
        let zero = [0.0, 0.0];
        // 1D のバランス (J+ + J−)/h + Σr φ̄_g − Σs(g'→g) φ̄_g' = 0 からノード平均
        let removal = [total[0] - scatter[0][0], total[1] - scatter[1][1]];
        let leak = [(plus[0] + minus[0]) / width, (plus[1] + minus[1]) / width];
        let det = removal[0] * removal[1] - scatter[1][0] * scatter[0][1];
        let flux = [
            (-leak[0] * removal[1] - scatter[1][0] * leak[1]) / det,
            (-leak[1] * removal[0] - scatter[0][1] * leak[0]) / det,
        ];
        let homogeneous = surface_fluxes(&node(
            &truth,
            &total,
            &scatter,
            &flux,
            [&minus, &plus, &zero, &zero],
        ))
        .unwrap();
        let factors = [0.92, 1.4];
        let het: Vec<Vec<f64>> = homogeneous
            .iter()
            .map(|face| (0..2).map(|g| factors[g] * face[g]).collect())
            .collect();
        let guess = [1.0, 0.3];
        let found = equivalent_diffusion(
            &node(
                &guess,
                &total,
                &scatter,
                &flux,
                [&minus, &plus, &zero, &zero],
            ),
            [&het[0], &het[1], &het[2], &het[3]],
        )
        .unwrap();
        for g in 0..2 {
            assert!((found[g] - truth[g]).abs() < 1e-9 * truth[g], "{found:?}");
        }
    }

    #[test]
    fn impossible_equivalences_are_errors() {
        let total = [0.6, 1.7];
        let scatter = vec![vec![0.57, 0.0], vec![0.0, 1.66]];
        let flux = [0.02, 0.01];
        let zero = [0.0, 0.0];
        let guess = [1.0, 0.3];
        // 流れが −x から +x へ減るのに、面の中性子束が +x で大きい（t ≥ 1）: どの D でも同じ因子にならない
        let minus = [-0.01, -0.002];
        let plus = [0.001, 0.0005];
        let result = equivalent_diffusion(
            &node(
                &guess,
                &total,
                &scatter,
                &flux,
                [&minus, &plus, &zero, &zero],
            ),
            [&[0.01, 0.005], &[0.012, 0.004], &flux, &flux],
        );
        assert!(
            matches!(
                result,
                Err(OneNodeError::NoEquivalentDiffusion { group: 0 })
            ),
            "{result:?}"
        );
        // 正味の流れがない（対称な位置）: 向きが定まらない
        let result = equivalent_diffusion(
            &node(&guess, &total, &scatter, &flux, [&zero; 4]),
            [&flux, &flux, &flux, &flux],
        );
        assert!(
            matches!(result, Err(OneNodeError::NoNetCurrent)),
            "{result:?}"
        );
        // 核分裂のあるノードは対象でない
        let mut fissile = node(
            &guess,
            &total,
            &scatter,
            &flux,
            [&minus, &plus, &zero, &zero],
        );
        fissile.nu_fission = &[0.0, 0.1];
        fissile.chi = &[1.0, 0.0];
        let result = equivalent_diffusion(&fissile, [&flux, &flux, &flux, &flux]);
        assert!(
            matches!(result, Err(OneNodeError::Invalid(_))),
            "{result:?}"
        );
    }
}
