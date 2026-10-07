//! 自由ガス模型による熱中性子の弾性散乱。
//!
//! 標的核を温度 T の理想気体（Maxwell 分布で運動する自由な原子）とみなし、実験室系の
//! 散乱を等方として、出射エネルギーの分布 S(E→E') を解析的に求める。
//!
//! S は自由原子の断面積を 1 としたときの値で、規格化しない。∫S dE' = N(E) は自由ガスの実効断面積と
//! 自由原子の断面積の比で、低エネルギーほど 1 より大きい。散乱行列では、温度に拡げた断面積
//! （一定の自由原子の断面積なら σ N(E)）に S/N(E) を掛ける。自由原子の断面積が一定なら σ S に等しく、
//! 詳細釣り合いが厳密に成り立つ。

use super::GroupStructure;
use super::integrate::GAUSS_LEGENDRE;
use crate::special::{erf, erfcx};

/// 熱領域の弾性散乱を自由ガス模型で扱う条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeGas {
    /// kT (eV)
    pub kt: f64,
    /// この入射エネルギー (eV) 未満で自由ガス模型を使う
    pub cutoff: f64,
}

/// 自由原子の散乱断面積を 1 としたときの、エネルギー移行の確率密度 S(E→E') (1/eV)。
/// `awr` は標的の質量比。
pub fn free_gas_kernel(e: f64, e_out: f64, awr: f64, kt: f64) -> f64 {
    if e <= 0.0 || e_out <= 0.0 {
        return 0.0;
    }
    let sqrt_a = awr.sqrt();
    let theta = (awr + 1.0) / (2.0 * sqrt_a);
    let zeta = (awr - 1.0) / (2.0 * sqrt_a);
    let (s, s_out) = ((e / kt).sqrt(), (e_out / kt).sqrt());
    // e^(-y²) erfcx(x) = e^(-y²+x²) erfc(x)。指数を先にまとめて桁あふれを避ける
    let scaled = |y: f64, x: f64| (-y * y).exp() * erfcx(x);
    let bracket = if e_out <= e {
        // erf(θ√ε' - ζ√ε) + erf(θ√ε' + ζ√ε) + e^(ε-ε') [erfc(θ√ε + ζ√ε') - erfc(θ√ε - ζ√ε')]
        let (x1, x2) = (theta * s + zeta * s_out, theta * s - zeta * s_out);
        let (u, v) = (zeta * s - theta * s_out, zeta * s + theta * s_out);
        if u > 0.0 {
            // 遠い下方散乱: O(1) の項が打ち消し合うので、同じ指数因子をもつ項どうしをまとめる
            (-u * u).exp() * (erfcx(u) - erfcx(x2)) + (-v * v).exp() * (erfcx(x1) - erfcx(v))
        } else {
            erf(-u) + erf(v) + scaled(v, x1) - scaled(u, x2)
        }
    } else {
        // erf(θ√ε' - ζ√ε) - erf(θ√ε' + ζ√ε) + e^(ε-ε') [erf(θ√ε - ζ√ε') + erf(θ√ε + ζ√ε')]
        let (p, q) = (theta * s_out - zeta * s, theta * s_out + zeta * s);
        let (r, t) = (theta * s - zeta * s_out, theta * s + zeta * s_out);
        if r < 0.0 {
            // 遠い上方散乱: 同様に、同じ指数因子をもつ項どうしをまとめる
            (-p * p).exp() * (erfcx(-r) - erfcx(p)) + (-q * q).exp() * (erfcx(q) - erfcx(t))
        } else {
            scaled(q, q) - scaled(p, p) + ((e - e_out) / kt).exp() * (erf(r) + erf(t))
        }
    };
    // 打ち消し合いの丸めで負になった裾は 0 にする
    (theta * theta / (2.0 * e) * bracket).max(0.0)
}

/// 入射エネルギー e での ∫_g' S(E→E') dE'（群ごと、和は N(E)）。
/// 群構造の外に散乱する分は端の群に入れる（ほかの反応と同じく、中性子を失わない）。
pub fn free_gas_group_integrals(groups: &GroupStructure, e: f64, awr: f64, kt: f64) -> Vec<f64> {
    let bounds = groups.bounds();
    let highest = bounds[0];
    // 核の幅の目安。重い核では E の近くに集中し、軽い核では kT 程度に広がる
    let width = 2.0 * (e * kt / awr).sqrt() + kt / awr;
    let reach = e + 30.0 * width + 40.0 * kt;
    // 下方散乱は静止標的なら αE まで届き、熱運動でさらに広がる
    let alpha = ((awr - 1.0) / (awr + 1.0)).powi(2);
    let lower = (alpha * e - 30.0 * width).max(0.0);
    let kernel = |x: f64| free_gas_kernel(e, x, awr, kt);

    // 区間の境界: 群境界、E' = E（核の折れ点）、αE、E の周りの核の幅ごとの点
    let mut points: Vec<f64> = bounds
        .iter()
        .copied()
        .chain((-60..=60).map(|k| e + 0.5 * f64::from(k) * width))
        .chain([lower, reach, e, alpha * e])
        .filter(|&x| x >= lower && x <= reach)
        .collect();
    points.sort_by(f64::total_cmp);
    points.dedup();

    let mut integrals = vec![0.0; groups.len()];
    for w in points.windows(2) {
        let (a, b) = (w[0], w[1]);
        let mid = 0.5 * (a + b);
        let g = groups
            .find(mid)
            .unwrap_or(if mid > highest { 0 } else { groups.len() - 1 });
        integrals[g] += if a > 0.0 {
            integrate_log(&kernel, a, b)
        } else {
            integrate_linear(&kernel, a, b)
        };
    }
    integrals
}

/// [a, b]（a > 0）をエネルギー比が一定以下になるよう対数的に分割して積分する。
pub(super) fn integrate_log(f: &impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let pieces = ((b / a).ln() / MAX_RATIO.ln()).ceil().max(1.0) as usize;
    let ratio = (b / a).powf(1.0 / pieces as f64);
    let mut lo = a;
    let mut sum = 0.0;
    for i in 0..pieces {
        let hi = if i + 1 == pieces { b } else { lo * ratio };
        sum += gauss_legendre(f, lo, hi);
        lo = hi;
    }
    sum
}

/// 0 を含む区間 [a, b] を等分して積分する。
pub(super) fn integrate_linear(f: &impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    const PIECES: usize = 16;
    let width = (b - a) / PIECES as f64;
    (0..PIECES)
        .map(|i| gauss_legendre(f, a + i as f64 * width, a + (i + 1) as f64 * width))
        .sum()
}

fn gauss_legendre(f: &impl Fn(f64) -> f64, lo: f64, hi: f64) -> f64 {
    let (mid, half) = (0.5 * (lo + hi), 0.5 * (hi - lo));
    GAUSS_LEGENDRE
        .iter()
        .map(|&(x, wk)| wk * f(mid + half * x))
        .sum::<f64>()
        * half
}

/// 1つの積分区間の上限と下限のエネルギー比の上限
const MAX_RATIO: f64 = 1.05;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::special::erf;

    const KT: f64 = 0.0253;

    fn assert_rel(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol * expected.abs(),
            "actual {actual}, expected {expected}"
        );
    }

    #[test]
    fn group_moments_keep_the_p0_integrals() {
        // P0 のモーメントは群積分と同じ。重い核の P1 は熱領域の上で静止した標的の 2/(3A) に近づく
        let groups = GroupStructure::wims69();
        for (awr, e) in [(1.0, 0.05), (15.86, 0.1), (236.0, 0.5)] {
            let moments = free_gas_group_moments(&groups, e, awr, KT, 3);
            let reference = free_gas_group_integrals(&groups, e, awr, KT);
            let total: f64 = reference.iter().sum();
            for (g, (m, r)) in moments[0].iter().zip(&reference).enumerate() {
                assert!((m - r).abs() <= 1e-6 * total, "A={awr} g={g}: {m} vs {r}");
            }
        }
        let moments = free_gas_group_moments(&groups, 2.0, 15.86, KT, 1);
        let mu: f64 = moments[1].iter().sum::<f64>() / moments[0].iter().sum::<f64>();
        assert!((mu - 2.0 / (3.0 * 15.86)).abs() < 2e-3, "{mu}");
    }

    #[test]
    fn free_gas_moments_match_the_s_alpha_beta_table() {
        // 自由ガスの S(α, β) の表（A = 1）から求めたモーメントと、表の補間の誤差の範囲で一致する
        use super::super::sab::ThermalSab;
        use super::super::sab::tests::{TEMPERATURE, free_gas_table};
        let data = free_gas_table();
        let sab = ThermalSab::new(&data, TEMPERATURE, 1.0).unwrap();
        let kt = crate::physics::BK * TEMPERATURE;
        let groups = GroupStructure::wims69();
        for e in [0.0253, 0.1, 0.5] {
            let ours = free_gas_group_moments(&groups, e, 1.0, kt, 3);
            // どちらも自由原子の断面積を 1 とした規格化（束縛原子の断面積 σb を含む）
            let table = sab.group_moments(&groups, e, 3);
            let total: f64 = ours[0].iter().sum();
            for l in 0..=3 {
                let (a, b): (f64, f64) = (ours[l].iter().sum(), table[l].iter().sum());
                assert!((a - b).abs() < 3e-3 * total, "E={e} P{l}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn hydrogen_has_closed_form() {
        // A=1: E' < E では erf(√(E'/kT))/E、E' > E では exp((E-E')/kT) erf(√(E/kT))/E
        let e = 0.05;
        for e_out in [0.001, 0.02, 0.0499] {
            assert_rel(
                free_gas_kernel(e, e_out, 1.0, KT),
                erf((e_out / KT).sqrt()) / e,
                1e-12,
            );
        }
        for e_out in [0.0501, 0.1, 0.3] {
            let expected = ((e - e_out) / KT).exp() * erf((e / KT).sqrt()) / e;
            assert_rel(free_gas_kernel(e, e_out, 1.0, KT), expected, 1e-12);
        }
    }

    #[test]
    fn detailed_balance() {
        // E exp(-E/kT) S(E→E') = E' exp(-E'/kT) S(E'→E)
        for awr in [1.0, 15.86, 236.0] {
            for (e, e_out) in [(0.01, 0.03), (0.05, 0.02), (1.0, 0.98), (2.0, 3.5)] {
                let forward = e * (-e / KT).exp() * free_gas_kernel(e, e_out, awr, KT);
                let backward = e_out * (-e_out / KT).exp() * free_gas_kernel(e_out, e, awr, KT);
                assert_rel(forward, backward, 1e-9);
            }
        }
    }

    #[test]
    fn integrates_to_effective_cross_section() {
        // ∫ S dE' = (1 + 1/(2Aε)) erf(√(Aε)) + exp(-Aε)/√(πAε)、ε = E/kT
        for awr in [1.0, 15.86, 236.0] {
            for e in [0.005, 0.0253, 0.5, 3.0] {
                let expected = effective_ratio(e, awr);
                // 台形公式で積分する（核の幅に比べて十分細かい格子）
                let width = (e * KT / awr).sqrt() + KT;
                let (lo, hi) = ((e - 60.0 * width).max(0.0), e + 60.0 * width);
                let n = 400_000;
                let h = (hi - lo) / n as f64;
                let integral: f64 = (0..=n)
                    .map(|i| {
                        let x = lo + i as f64 * h;
                        let f = free_gas_kernel(e, x, awr, KT);
                        if i == 0 || i == n { 0.5 * f } else { f }
                    })
                    .sum::<f64>()
                    * h;
                assert_rel(integral, expected, 1e-5);
            }
        }
    }

    #[test]
    fn integrates_down_to_alpha_e_when_e_is_far_above_kt() {
        // 下方散乱は αE（α = ((A-1)/(A+1))²）まで届く。重水素 20 K、E/kT ≈ 2300
        let (awr, kt, e) = (1.99, 20.0 * crate::physics::BK, 3.99);
        let groups = GroupStructure::wims69();
        let sum: f64 = free_gas_group_integrals(&groups, e, awr, kt).iter().sum();
        let a_eps = awr * e / kt;
        let expected = (1.0 + 1.0 / (2.0 * a_eps)) * erf(a_eps.sqrt())
            + (-a_eps).exp() / (std::f64::consts::PI * a_eps).sqrt();
        assert_rel(sum, expected, 1e-6);
    }

    #[test]
    fn upscatter_above_top_goes_to_first_group() {
        // 熱群だけの群構造でも、上限より上に散乱した分を失わない
        let groups = GroupStructure::new(vec![1.0, 0.01]).unwrap();
        let (awr, e) = (1.0, 0.9);
        let sum: f64 = free_gas_group_integrals(&groups, e, awr, KT).iter().sum();
        assert_rel(sum, effective_ratio(e, awr), 1e-6);
    }

    #[test]
    fn no_overflow_far_below_incident_energy() {
        // E/kT ≈ 160 で E' がずっと低くても有限
        let s = free_gas_kernel(4.0, 0.01, 15.86, KT);
        assert!(s.is_finite() && s >= 0.0, "{s}");
    }

    fn effective_ratio(e: f64, awr: f64) -> f64 {
        let a_eps = awr * e / KT;
        (1.0 + 1.0 / (2.0 * a_eps)) * erf(a_eps.sqrt())
            + (-a_eps).exp() / (std::f64::consts::PI * a_eps).sqrt()
    }

    #[test]
    fn group_integrals_sum_to_effective_ratio() {
        let groups = GroupStructure::wims69();
        for awr in [1.0, 15.86, 236.0] {
            for e in [1.0e-3, 0.0253, 0.9, 3.9] {
                let p = free_gas_group_integrals(&groups, e, awr, KT);
                assert_rel(p.iter().sum::<f64>(), effective_ratio(e, awr), 1e-6);
                assert!(p.iter().all(|&v| v >= 0.0));
            }
        }
        // 重い核ではほとんど同じ群に残る（群 [0.025, 0.03] eV の中央から入射）
        let e = 0.0275;
        let g = groups.find(e).unwrap();
        assert!(free_gas_group_integrals(&groups, e, 236.0, KT)[g] > 0.5);
    }
}
