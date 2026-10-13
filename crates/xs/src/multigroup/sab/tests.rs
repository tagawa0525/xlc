use super::*;
use crate::multigroup::{free_gas_group_integrals, free_gas_kernel};

pub(crate) const TEMPERATURE: f64 = 293.6;

/// 自由ガス（A = 1、σfree = 1 b）の対称な S(α, β) = e^(-(α²+β²)/(4α))/√(4πα) を細かい格子で表にしたもの。
/// 小さな α では S が β について幅 √α ほどで急に変わるので、実際の表と同じく β も原点の近くを細かくとる
pub(crate) fn free_gas_table() -> Mf7Incoherent {
    let alpha: Vec<f64> = (0..=800)
        .map(|i| 1.0e-5 * 10f64.powf(i as f64 / 100.0))
        .collect();
    let beta: Vec<f64> = std::iter::once(0.0)
        .chain((0..=750).map(|i| 1.0e-5 * 10f64.powf(i as f64 / 100.0)))
        .collect();
    let s = beta
        .iter()
        .flat_map(|&b| {
            alpha.iter().map(move |&a| {
                (-(a * a + b * b) / (4.0 * a)).exp() / (4.0 * std::f64::consts::PI * a).sqrt()
            })
        })
        .collect();
    Mf7Incoherent {
        lat: 0,
        lasym: 0,
        sigma_free: 1.0,
        awr: 1.0,
        e_max: 10.0,
        alpha,
        alpha_law: Interpolation::LogLin,
        beta,
        beta_law: Interpolation::LogLin,
        temperatures: vec![TEMPERATURE],
        s: vec![s],
        effective_temperature: crate::endf::Tab1::linear(vec![TEMPERATURE], vec![TEMPERATURE]),
    }
}

/// `free_gas_table` を1桁あたり 20 点に粗くした表。補間の誤差は大きいが、詳細バランスのように
/// 表の細かさによらない性質を、散乱行列の段階で安く確かめるのに使う。
pub(crate) fn coarse_free_gas_table() -> Mf7Incoherent {
    let fine = free_gas_table();
    let (na, nb) = (fine.alpha.len(), fine.beta.len());
    let keep_alpha: Vec<usize> = (0..na).step_by(5).collect();
    let keep_beta: Vec<usize> = std::iter::once(0).chain((1..nb).step_by(5)).collect();
    let s = keep_beta
        .iter()
        .flat_map(|&j| keep_alpha.iter().map(move |&i| (j, i)))
        .map(|(j, i)| fine.s[0][j * na + i])
        .collect();
    Mf7Incoherent {
        alpha: keep_alpha.iter().map(|&i| fine.alpha[i]).collect(),
        beta: keep_beta.iter().map(|&j| fine.beta[j]).collect(),
        s: vec![s],
        ..fine
    }
}

#[test]
fn free_gas_table_reproduces_free_gas_kernel() {
    // 表の補間（ln S を α, β について線形）の誤差の範囲で、自由ガスの解析式と一致する
    let data = free_gas_table();
    let sab = ThermalSab::new(&data, TEMPERATURE, 1.0, SabOutside::Zero).unwrap();
    let kt = BK * TEMPERATURE;
    let groups = GroupStructure::wims69();
    for e in [0.01, 0.0253, 0.1, 0.5] {
        let ours = sab.group_integrals(&groups, e);
        let reference = free_gas_group_integrals(&groups, e, 1.0, kt);
        let total: f64 = reference.iter().sum();
        let ours_total: f64 = ours.iter().sum();
        assert!(
            (ours_total / total - 1.0).abs() < 1e-3,
            "E={e}: {ours_total} vs {total}"
        );
        for (g, (o, r)) in ours.iter().zip(&reference).enumerate() {
            assert!(
                (o - r).abs() <= 2e-3 * r + 1e-6 * total,
                "E={e}, g={g}: {o} vs {r}"
            );
        }
        let e_out = 0.7 * e;
        let (k, k_ref) = (sab.kernel(e, e_out), free_gas_kernel(e, e_out, 1.0, kt));
        assert!((k / k_ref - 1.0).abs() < 2e-3, "E={e}: {k} vs {k_ref}");
    }
}

#[test]
fn kernel_moments_match_the_free_gas_double_differential() {
    // 自由ガス（A = 1、σfree = 1 b、σb = 4 b）の σ(E→E', μ) = σb/(2kT) √(E'/E) e^(-β/2) S(α, β) を
    // μ で細かく積分した P_l のモーメントと、表から求めたものが、表の補間の誤差の範囲で一致する
    let data = free_gas_table();
    let sab = ThermalSab::new(&data, TEMPERATURE, 1.0, SabOutside::Zero).unwrap();
    let kt = BK * TEMPERATURE;
    let reference = |e: f64, e_out: f64, l: usize| {
        let n = 200_000;
        let beta = (e_out - e) / kt;
        let f = |mu: f64| {
            let alpha = (e + e_out - 2.0 * mu * (e * e_out).sqrt()) / kt;
            let s = (-(alpha * alpha + beta * beta) / (4.0 * alpha)).exp()
                / (4.0 * std::f64::consts::PI * alpha).sqrt();
            4.0 / (2.0 * kt) * (e_out / e).sqrt() * (-0.5 * beta).exp() * s * legendre(mu, l)[l]
        };
        // 中点則（μ = 1 の近くでも α > 0）
        (0..n)
            .map(|i| f(-1.0 + (i as f64 + 0.5) * 2.0 / n as f64))
            .sum::<f64>()
            * 2.0
            / n as f64
    };
    for (e, e_out) in [(0.05, 0.03), (0.05, 0.07), (0.0253, 0.02), (0.2, 0.15)] {
        let moments = sab.kernel_moments(e, e_out, 3);
        let p0 = reference(e, e_out, 0);
        assert!(
            (moments[0] / p0 - 1.0).abs() < 2e-3,
            "{e} {e_out}: {moments:?} vs {p0}"
        );
        for (l, m) in moments.iter().enumerate().skip(1) {
            let r = reference(e, e_out, l);
            assert!((m - r).abs() < 2e-3 * p0, "{e} {e_out} P{l}: {} vs {r}", m);
        }
    }
}

#[test]
fn kernel_satisfies_detailed_balance() {
    // E σ0(E→E') e^(-E/kT) は E と E' について対称
    let data = free_gas_table();
    let sab = ThermalSab::new(&data, TEMPERATURE, 1.0, SabOutside::Zero).unwrap();
    let kt = BK * TEMPERATURE;
    for (e, e_out) in [(0.01, 0.05), (0.2, 0.03), (0.0253, 0.026)] {
        let forward = e * sab.kernel(e, e_out) * (-e / kt).exp();
        let backward = e_out * sab.kernel(e_out, e) * (-e_out / kt).exp();
        assert!(
            (forward / backward - 1.0).abs() < 1e-12,
            "{forward} vs {backward}"
        );
    }
}

/// `free_gas_table` の S を、条件 `zero`(α, β) を満たす格子点で 0 にしたもの
fn table_with_zeros(zero: impl Fn(f64, f64) -> bool) -> Mf7Incoherent {
    let mut data = free_gas_table();
    let na = data.alpha.len();
    for (j, &b) in data.beta.iter().enumerate() {
        for (i, &a) in data.alpha.iter().enumerate() {
            if zero(a, b) {
                data.s[0][j * na + i] = 0.0;
            }
        }
    }
    data
}

/// σb A/(4E) ∫ e^(-β/2) S(α, β) dα の短時間衝突近似の S（有効温度 `t_eff` K）での値を、α の中点則で求めたもの
/// （A = 1、σfree = 1 b、σb = 4 b）
fn short_collision_time_kernel(e: f64, e_out: f64, t_eff: f64) -> f64 {
    let kt = BK * TEMPERATURE;
    let r = TEMPERATURE / t_eff;
    let beta = (e_out - e) / kt;
    let b = beta.abs();
    let (lo, hi) = (
        (e.sqrt() - e_out.sqrt()).powi(2) / kt,
        (e.sqrt() + e_out.sqrt()).powi(2) / kt,
    );
    let n = 400_000;
    let h = (hi - lo) / n as f64;
    let integral: f64 = (0..n)
        .map(|i| {
            let a = lo + (i as f64 + 0.5) * h;
            (-(a - b).powi(2) * r / (4.0 * a) - 0.5 * b - 0.5 * beta).exp()
                / (4.0 * std::f64::consts::PI * a / r).sqrt()
        })
        .sum::<f64>()
        * h;
    4.0 / (4.0 * e) * integral
}

#[test]
fn short_collision_time_fills_alpha_above_the_table() {
    // α の上限が 1 の表では、0.5 eV からの散乱の大半が表の外にある。有効温度が物質の温度なら短時間衝突近似は
    // 自由ガスの S と同じなので、自由ガスの核を表の補間の誤差の範囲で再現する
    let mut data = free_gas_table();
    let na = data.alpha.len();
    let keep = data.alpha.partition_point(|&a| a <= 1.0);
    data.s[0] = (0..data.beta.len())
        .flat_map(|j| data.s[0][j * na..j * na + keep].to_vec())
        .collect();
    data.alpha.truncate(keep);
    // 0 とするなら、表が下方散乱を覆わないのでエラー
    assert!(matches!(
        ThermalSab::new(&data, TEMPERATURE, 1.0, SabOutside::Zero),
        Err(SabError::Coverage { .. })
    ));
    let sab = ThermalSab::new(&data, TEMPERATURE, 1.0, SabOutside::ShortCollisionTime).unwrap();
    let kt = BK * TEMPERATURE;
    let groups = GroupStructure::wims69();
    let e = 0.5;
    let ours: f64 = sab.group_integrals(&groups, e).iter().sum();
    let reference: f64 = free_gas_group_integrals(&groups, e, 1.0, kt).iter().sum();
    assert!(
        (ours / reference - 1.0).abs() < 2e-3,
        "{ours} vs {reference}"
    );
    for e_out in [0.1, 0.45, 0.55] {
        let (k, k_ref) = (sab.kernel(e, e_out), free_gas_kernel(e, e_out, 1.0, kt));
        assert!((k / k_ref - 1.0).abs() < 2e-3, "{e_out}: {k} vs {k_ref}");
    }
}

#[test]
fn short_collision_time_uses_the_effective_temperature() {
    // 表の α の範囲より大きい α だけに散乱する (E, E') では、核は有効温度の短時間衝突近似そのもの
    let mut data = free_gas_table();
    let t_eff = 2.0 * TEMPERATURE;
    data.effective_temperature = crate::endf::Tab1::linear(vec![TEMPERATURE], vec![t_eff]);
    let na = data.alpha.len();
    let keep = data.alpha.partition_point(|&a| a <= 1.0);
    data.s[0] = (0..data.beta.len())
        .flat_map(|j| data.s[0][j * na..j * na + keep].to_vec())
        .collect();
    data.alpha.truncate(keep);
    let sab = ThermalSab::new(&data, TEMPERATURE, 4.0, SabOutside::ShortCollisionTime).unwrap();
    for (e, e_out) in [(3.0, 1.0), (2.0, 0.5), (1.0, 2.0)] {
        let (k, k_ref) = (
            sab.kernel(e, e_out),
            short_collision_time_kernel(e, e_out, t_eff),
        );
        assert!(
            (k / k_ref - 1.0).abs() < 1e-6,
            "{e} → {e_out}: {k} vs {k_ref}"
        );
    }
}

#[test]
fn short_collision_time_fills_zeros_only_at_large_alpha_or_beta() {
    // NJOY THERMR と同じく、S = 0 の格子点に接する区間は α A ≥ 30 か β ≥ 30 のところだけ短時間衝突近似で補う。
    // 2 → 1.4 eV（β ≈ −24）の S は α ≈ 24 の近くに幅 7 ほどで集まる
    let kt = BK * TEMPERATURE;
    let (e, e_out) = (2.0, 1.4);
    let reference = free_gas_kernel(e, e_out, 1.0, kt);
    // α ≥ 32 の 0 は補う
    let data = table_with_zeros(|a, _| a >= 32.0);
    let sab = ThermalSab::new(&data, TEMPERATURE, 2.5, SabOutside::ShortCollisionTime).unwrap();
    let k = sab.kernel(e, e_out);
    assert!((k / reference - 1.0).abs() < 2e-3, "{k} vs {reference}");
    // α < 30 かつ β < 30 の 0 は補わない
    let data = table_with_zeros(|a, b| (1.0..=29.0).contains(&a) && b < 29.0);
    let sab = ThermalSab::new(&data, TEMPERATURE, 2.5, SabOutside::ShortCollisionTime).unwrap();
    let k = sab.kernel(e, e_out);
    assert!(k < 0.5 * reference, "{k} vs {reference}");
    // 0 とするなら、α ≥ 32 の 0 も補わない
    let data = table_with_zeros(|a, _| a >= 32.0);
    let sab = ThermalSab::new(&data, TEMPERATURE, 2.5, SabOutside::Zero).unwrap();
    let k = sab.kernel(e, e_out);
    assert!(k < 0.95 * reference, "{k} vs {reference}");
}

#[test]
fn far_downscatter_at_low_temperature_stays_finite() {
    // LAT=1 の表を kT の小さい温度で引くと、物質の温度の β は表の値よりずっと大きく、e^(-β/2) を別に計算すると
    // 桁あふれする。9.5 → 0.5 eV は表の β の上限の外で、短時間衝突近似だけの値
    let mut data = free_gas_table();
    data.lat = 1;
    let temperature = 20.0;
    data.temperatures = vec![temperature];
    data.effective_temperature = crate::endf::Tab1::linear(vec![temperature], vec![temperature]);
    let sab = ThermalSab::new(&data, temperature, 10.0, SabOutside::ShortCollisionTime).unwrap();
    let k = sab.kernel_moments(9.5, 0.5, 3);
    assert!(k.iter().all(|v| v.is_finite()), "{k:?}");
}

#[test]
fn tiny_s_counts_as_zero_like_thermr() {
    // THERMR は S ≤ e^-225 を 0 とみなす（thermr.f90 の sabflg）。そのような格子点も、α A ≥ 30 なら補う
    let kt = BK * TEMPERATURE;
    let (e, e_out) = (2.0, 1.4);
    let reference = free_gas_kernel(e, e_out, 1.0, kt);
    let mut data = table_with_zeros(|a, _| a >= 32.0);
    for v in data.s[0].iter_mut().filter(|v| **v == 0.0) {
        *v = 1.0e-120;
    }
    let sab = ThermalSab::new(&data, TEMPERATURE, 2.5, SabOutside::ShortCollisionTime).unwrap();
    let k = sab.kernel(e, e_out);
    assert!((k / reference - 1.0).abs() < 2e-3, "{k} vs {reference}");
}

#[test]
fn cutoff_above_the_evaluated_limit_is_error() {
    // 表の上限 B(4) より上では S(α, β) を使えない（表の外を補う扱いでも）
    let data = free_gas_table();
    for outside in [SabOutside::ShortCollisionTime, SabOutside::Zero] {
        assert!(matches!(
            ThermalSab::new(&data, TEMPERATURE, 1.5 * data.e_max, outside),
            Err(SabError::AboveEvaluated { .. })
        ));
    }
}

#[test]
fn empty_grid_is_error() {
    let mut data = free_gas_table();
    data.beta.clear();
    assert!(matches!(
        ThermalSab::new(&data, TEMPERATURE, 1.0, SabOutside::Zero),
        Err(SabError::EmptyGrid)
    ));
}

#[test]
fn rejects_missing_temperature_and_short_table() {
    let data = free_gas_table();
    assert!(ThermalSab::new(&data, 600.0, 1.0, SabOutside::Zero).is_err());
    // β の上限 316 は kT ≈ 0.0253 eV で 8 eV 分なので、10 eV までの下方散乱を覆わない
    assert!(ThermalSab::new(&data, TEMPERATURE, 10.0, SabOutside::Zero).is_err());
}
