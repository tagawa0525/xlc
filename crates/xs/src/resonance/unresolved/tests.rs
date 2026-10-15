use super::super::factors::{penetrability, phase_shift};
use super::*;
use crate::endf::{OrbitalL, Tab1, UrrLValue};
use crate::special::erfc;

const AWRI: f64 = 236.0;
const AP: f64 = 0.948;

/// テストのサンプルの数。ライブラリの数（`Sampling::LIBRARY`）は黒い極限のテストだけで使う
const TEST: Sampling = Sampling {
    ladders: 64,
    samples: 250,
    margin: 1.0,
};

fn point(e: f64, gno: f64) -> UrrPoint {
    UrrPoint {
        e,
        d: 20.0,
        gx: 0.0,
        gno,
        gg: 0.023,
        gf: 0.0,
    }
}

fn s_wave(points: Vec<UrrPoint>, int: Interpolation) -> Unresolved {
    Unresolved {
        lssf: 0,
        spi: 0.0,
        ap: AP,
        ap_energy: None,
        naps: 1,
        l_values: vec![UrrLValue {
            awri: AWRI,
            l: OrbitalL::S,
            j_values: vec![UrrJValue {
                aj: 0.5,
                int,
                amux: WidthDof::Fixed,
                amun: WidthDof::One,
                amuf: WidthDof::Fixed,
                points,
            }],
        }],
    }
}

/// 平均 1 の χ²_1 での E[1/(x+c)] = c^(-1/2) e^(c/2) √(π/2) erfc(√(c/2))。
fn inverse_moment(c: f64) -> f64 {
    (PI / (2.0 * c)).sqrt() * (0.5 * c).exp() * erfc((0.5 * c).sqrt())
}

#[test]
fn energy_outside_the_radius_table_is_error() {
    let mut urr = s_wave(
        vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
        Interpolation::LinLin,
    );
    urr.ap_energy = Some(Tab1::linear(vec![5.0e3, 1.0e4], vec![AP, AP]));
    let xs = UnresolvedXs::new(&urr).unwrap();
    assert!(xs.cross_sections(8.0e3).is_ok());
    assert!(matches!(
        xs.cross_sections(1.5e4),
        Err(UrrError::RadiusTable { .. })
    ));
}

#[test]
fn s_wave_without_fission_matches_closed_form() {
    // GNO は 5e3 eV で 1e-3、2e4 eV で 4e-3 なので、1e4 eV では lin-lin で 2e-3
    let urr = s_wave(
        vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
        Interpolation::LinLin,
    );
    let e = 1.0e4;
    let sigma = UnresolvedXs::new(&urr).unwrap().cross_sections(e).unwrap();

    let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
    let (d, gg) = (20.0, 0.023);
    let gn = 2.0e-3 * e.sqrt();
    // 中性子幅 gn x、分母 gn x + gg の期待値（c = gg/gn）
    let c = gg / gn;
    let e_inv = inverse_moment(c);
    let e_x = (1.0 - c * e_inv) / gn;
    let e_x2 = (1.0 - c + c * c * e_inv) / gn;
    let scale = 2.0 * PI * PI / (k * k) / d;
    let phase = (k * AP).sin().powi(2);
    let capture = scale * gn * gg * e_x;
    let elastic = scale * gn * gn * e_x2 - scale * 2.0 * gn * phase + 4.0 * PI * phase / (k * k);
    for (actual, expected) in [(sigma.capture, capture), (sigma.elastic, elastic)] {
        assert!(
            (actual - expected).abs() <= 1e-8 * expected,
            "{actual} vs {expected}"
        );
    }
    assert_eq!(sigma.fission, 0.0);
    assert!((sigma.total - (sigma.elastic + sigma.capture)).abs() <= 1e-12 * sigma.total);
}

#[test]
fn grid_fills_wide_steps_like_njoy() {
    let urr = s_wave(
        vec![
            point(1.5e3, 1.0e-3),
            point(3.0e3, 1.0e-3),
            point(1.0e4, 1.0e-3),
            point(1.0e5, 1.0e-3),
        ],
        Interpolation::LinLin,
    );
    let expected = vec![
        1.5e3, 1.7e3, 2.0e3, 2.5e3, 3.0e3, 3.5e3, 4.0e3, 5.0e3, 6.0e3, 7.2e3, 8.5e3, 1.0e4, 1.25e4,
        1.5e4, 1.7e4, 2.0e4, 2.5e4, 3.0e4, 3.5e4, 4.0e4, 5.0e4, 6.0e4, 7.2e4, 8.5e4, 1.0e5,
    ];
    assert_eq!(unresolved_grid(&urr, 1.5e3, 1.0e5), expected);
}

#[test]
fn fluctuating_fission_width_converges() {
    // U-234（ENDF/B-VIII.0）の l=1, J=1/2 の系列。中性子幅と核分裂幅の入れ子の積分で、
    // 内側の積分のばらつきのために外側が収束しなかった
    let point = |e, gno, gf| UrrPoint {
        e,
        d: 10.6,
        gx: 0.0,
        gno,
        gg: 0.02,
        gf,
    };
    let urr = Unresolved {
        lssf: 0,
        spi: 0.0,
        ap: 0.893,
        ap_energy: None,
        naps: 0,
        l_values: vec![UrrLValue {
            awri: 232.029,
            l: OrbitalL::P,
            j_values: vec![UrrJValue {
                aj: 0.5,
                int: Interpolation::LinLin,
                amux: WidthDof::Fixed,
                amun: WidthDof::One,
                amuf: WidthDof::One,
                points: vec![
                    point(4.0e4, 1.5872e-3, 6.7737e-4),
                    point(6.0e4, 1.4566e-3, 6.0163e-4),
                ],
            }],
        }],
    };
    let sigma = UnresolvedXs::new(&urr)
        .unwrap()
        .cross_sections(5.0e4)
        .unwrap();
    assert!(sigma.fission > 0.0 && sigma.fission.is_finite());
}

#[test]
fn three_fluctuating_widths_converge() {
    // U-236（ENDF/B-VIII.0）の l=2, J=3/2 の系列。中性子・核分裂・競合の幅がゆらぐ
    let point = |e, gx| UrrPoint {
        e,
        d: 9.18,
        gx,
        gno: 1.836e-3,
        gg: 0.019,
        gf: 3.4e-4,
    };
    let urr = Unresolved {
        lssf: 0,
        spi: 0.0,
        ap: 0.93543,
        ap_energy: None,
        naps: 0,
        l_values: vec![UrrLValue {
            awri: 234.0178,
            l: OrbitalL::D,
            j_values: vec![UrrJValue {
                aj: 1.5,
                int: Interpolation::LogLog,
                amux: WidthDof::One,
                amun: WidthDof::One,
                amuf: WidthDof::One,
                points: vec![
                    point(4.5e4, 0.0),
                    point(5.0e4, 2.03e-5),
                    point(6.0e4, 3.33e-4),
                ],
            }],
        }],
    };
    let sigma = UnresolvedXs::new(&urr)
        .unwrap()
        .cross_sections(5.0e4)
        .unwrap();
    assert!(sigma.fission > 0.0 && sigma.fission.is_finite());
}

/// 種 7 の、σ0 ごとのサンプリングした実効断面積
fn sampled(xs: &UnresolvedXs, e: f64, temperature: f64, sigma0: &[f64]) -> Bondarenko {
    xs.shielded(e, temperature, sigma0, Sigma::default(), TEST, 7)
        .unwrap()
}

/// 種 1〜16 の独立なサンプリングでの量 `of` の平均と、平均の標準誤差
fn over_seeds(
    xs: &UnresolvedXs,
    e: f64,
    temperature: f64,
    sigma0: &[f64],
    of: impl Fn(&Bondarenko) -> f64,
) -> (f64, f64) {
    let values: Vec<f64> = (1..=16)
        .map(|seed| {
            of(&xs
                .shielded(e, temperature, sigma0, Sigma::default(), TEST, seed)
                .unwrap())
        })
        .collect();
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
    (mean, (variance / n).sqrt())
}

/// サンプリングの平均が、期待値 `expected` に統計誤差の 4 倍（と相対 1e-3 の系統の余裕）の内で合う
fn assert_within_statistics((mean, error): (f64, f64), expected: f64, what: &str) {
    assert!(
        (mean - expected).abs() <= 4.0 * error + 1e-3 * expected.abs(),
        "{what}: {mean} ± {error} vs {expected}"
    );
}

#[test]
fn sampled_mean_matches_the_width_fluctuation_average() {
    // ラダーの平均の期待値は、幅のゆらぎを積分した無限希釈の平均断面積。捕獲は単準位と多準位で同じで、弾性散乱は
    // 多準位の干渉の分だけ違いうるので、捕獲で確かめる
    let urr = s_wave(
        vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
        Interpolation::LinLin,
    );
    let xs = UnresolvedXs::new(&urr).unwrap();
    let e = 1.0e4;
    let average = xs.cross_sections(e).unwrap();
    let capture = over_seeds(&xs, e, 300.0, &[1.0e10], |b| b.mean.capture);
    assert_within_statistics(capture, average.capture, "capture");
}

#[test]
fn bondarenko_average_at_infinite_dilution_is_the_sampled_mean() {
    let urr = s_wave(
        vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
        Interpolation::LinLin,
    );
    let xs = UnresolvedXs::new(&urr).unwrap();
    let b = sampled(&xs, 1.0e4, 300.0, &[1.0e10]);
    let [dilute] = b.shielded[..] else {
        panic!("one σ0");
    };
    for (shielded, mean) in [
        (dilute.capture, b.mean.capture),
        (dilute.elastic, b.mean.elastic),
        (dilute.total, b.mean.total),
        (dilute.current_total, b.mean.total),
    ] {
        assert!(
            (shielded - mean).abs() <= 1e-6 * mean,
            "{shielded} vs {mean}"
        );
    }
}

#[test]
fn isolated_resonances_match_the_direct_integral_over_one_spacing() {
    // 準位間隔 D に比べて共鳴が十分狭いと、ラダーの平均は 1 つの共鳴の区間 D の積分の、幅のゆらぎの期待値になる:
    // ⟨g(σ)⟩ = (1/D) E_x[∫ g(σ(E)) dE]（Γn = Γ̄n x、x は自由度 1 の χ²/1。x = t²、t は標準正規）。捕獲の実効断面積
    // ⟨σc/(σt + σ0)⟩/⟨1/(σt + σ0)⟩ を、単準位の ψ, χ を直接に数値積分したものと比べる
    let urr = s_wave(
        vec![point(5.0e3, 2.0e-6), point(2.0e4, 2.0e-6)],
        Interpolation::LinLin,
    );
    let xs = UnresolvedXs::new(&urr).unwrap();
    let (e, temperature, sigma0): (f64, f64, f64) = (1.0e4, 300.0, 1.0);
    let k = crate::physics::wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
    let rho = k * AP;
    let phi = phase_shift(OrbitalL::S, rho);
    let ab = 4.0 * PI / (k * k);
    let mean_gn = 2.0e-6 * penetrability(OrbitalL::S, rho) / rho * e.sqrt();
    let gg = 0.023;
    let delta = 2.0 * (e * crate::physics::BK * temperature / AWRI).sqrt();
    let potential = ab * phi.sin().powi(2);
    let d = 20.0;
    // t の Simpson 則（[0, 8]、偶関数なので片側）と、E の中点則（幅 Γ より十分細かい 1 meV。区間 D の全体）
    let (nt, ne) = (400, 20_000);
    let (mut num, mut den, mut weight) = (0.0, 0.0, 0.0);
    for i in 0..=nt {
        let t = 8.0 * i as f64 / nt as f64;
        let simpson = if i == 0 || i == nt {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let wt = simpson * (-0.5 * t * t).exp();
        let gn = mean_gn * t * t;
        let gamma = gn + gg;
        let y = 0.5 * gamma / delta;
        let (rn, rg) = (gn / gamma, gg / gamma);
        let (mut n_e, mut d_e) = (0.0, 0.0);
        for j in 0..ne {
            let de = -0.5 * d + d * (j as f64 + 0.5) / ne as f64;
            let w = crate::special::faddeeva(num_complex::Complex64::new(de / delta, y));
            let (psi, chi) = (PI.sqrt() * y * w.re, PI.sqrt() * y * w.im);
            let capture = ab * rn * rg * psi;
            let elastic = potential
                + ab * rn * ((rn - 2.0 * phi.sin().powi(2)) * psi + (2.0 * phi).sin() * chi);
            let total = capture + elastic;
            n_e += capture / (total + sigma0);
            d_e += 1.0 / (total + sigma0);
        }
        num += wt * n_e;
        den += wt * d_e;
        weight += wt;
    }
    let expected = (num / weight) / (den / weight);
    let capture = over_seeds(&xs, e, temperature, &[sigma0], |b| b.shielded[0].capture);
    assert_within_statistics(capture, expected, "capture");
}

#[test]
fn black_limit_stays_positive_below_the_mean() {
    // 準位間隔 3 eV に幅の大きい共鳴が並ぶ低いエネルギー（Pm-149 の 2.6 eV と同じ形。#0127）でも、σ0 = 1 b の実効断面積は
    // 正で、平均より小さい
    let mut low = point(2.6, 1.0e-3);
    low.d = 3.0;
    low.gg = 0.065;
    let mut high = low;
    high.e = 1.0e3;
    let mut urr = s_wave(vec![low, high], Interpolation::LinLin);
    urr.l_values[0].awri = 147.0;
    let xs = UnresolvedXs::new(&urr).unwrap();
    let b = xs
        .shielded(
            2.6,
            565.0,
            &[1.0e10, 10.0, 1.0],
            Sigma::default(),
            Sampling::LIBRARY,
            7,
        )
        .unwrap();
    for s in &b.shielded[1..] {
        assert!(
            s.capture > 0.0 && s.capture < b.mean.capture,
            "{} vs {}",
            s.capture,
            b.mean.capture
        );
        assert!(s.total > 0.0 && s.current_total > 0.0);
    }
    assert!(b.shielded[2].capture < b.shielded[1].capture);
}

#[test]
fn ladders_follow_the_seed_and_keep_it_across_temperatures() {
    let urr = s_wave(
        vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
        Interpolation::LinLin,
    );
    let xs = UnresolvedXs::new(&urr).unwrap();
    let sigma0 = [1.0e10, 10.0];
    let a = sampled(&xs, 1.0e4, 300.0, &sigma0);
    assert_eq!(a.shielded, sampled(&xs, 1.0e4, 300.0, &sigma0).shielded);
    let other = xs
        .shielded(1.0e4, 300.0, &sigma0, Sigma::default(), TEST, 8)
        .unwrap();
    assert_ne!(a.shielded, other.shielded);
    // 同じ種なら、温度を変えてもラダーとエネルギーは同じ: 実効断面積は温度の滑らかな関数になり、1 K ずつの差がほぼ等しい。
    // 温度ごとに別のラダーなら、差は種を変えた揺らぎで乱れる
    let at = |t: f64| sampled(&xs, 1.0e4, t, &sigma0).shielded[1].capture;
    let (d1, d2) = (at(301.0) - at(300.0), at(302.0) - at(301.0));
    let independent = (other.shielded[1].capture - a.shielded[1].capture).abs();
    assert!(0.05 * d1.abs() < independent, "{d1} vs {independent}");
    assert!((d2 - d1).abs() < 0.05 * d1.abs(), "{d1} vs {d2}");
    // 高い温度ほど遮蔽が弱い
    let hot = sampled(&xs, 1.0e4, 1200.0, &sigma0);
    assert!(hot.shielded[1].capture > a.shielded[1].capture);
}

#[test]
fn naps_two_needs_radius_table() {
    // NAPS=2 は、散乱半径の表（NRO=1）とは別に AP をチャネル半径に使う指定
    let mut urr = s_wave(vec![point(1.5e3, 1.0e-3)], Interpolation::LinLin);
    urr.naps = 2;
    assert!(UnresolvedXs::new(&urr).is_err());
}

#[test]
fn interpolation_must_agree_across_sequences() {
    let mut urr = s_wave(
        vec![point(1.5e3, 1.0e-3), point(1.0e5, 1.0e-3)],
        Interpolation::LogLog,
    );
    assert_eq!(interpolation(&urr).unwrap(), Interpolation::LogLog);
    let mut other = urr.l_values[0].j_values[0].clone();
    other.int = Interpolation::LinLin;
    urr.l_values[0].j_values.push(other);
    assert!(interpolation(&urr).is_err());
}

#[test]
fn negative_elastic_background_floors_the_elastic_at_one_microbarn() {
    // MF3 の弾性散乱の背景が負で、ラダーの谷で弾性散乱 + 背景が負になる点は、NJOY PURR と同じく 1 μb にする（ADR-0177）。
    // 背景が −10⁴ b なら、すべての点が 1 μb になる
    let urr = s_wave(
        vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
        Interpolation::LinLin,
    );
    let xs = UnresolvedXs::new(&urr).unwrap();
    let sigma0 = [1.0e10, 10.0];
    let background = Sigma {
        total: -1.0e4,
        elastic: -1.0e4,
        fission: 0.0,
        capture: 0.0,
    };
    let b = xs
        .shielded(1.0e4, 300.0, &sigma0, background, TEST, 7)
        .unwrap();
    let samples = TEST.ladders * TEST.samples;
    assert_eq!(b.floored, samples);
    assert!(
        (b.mean.elastic - 1.0e-6).abs() < 1.0e-15,
        "{}",
        b.mean.elastic
    );
    for s in &b.shielded {
        assert!((s.elastic - 1.0e-6).abs() < 1.0e-15, "{}", s.elastic);
        assert!((s.total - s.elastic - s.capture).abs() < 1.0e-9 * s.total);
    }
    // 背景がなければ、弾性散乱は |1 − U|² で負にならず、どの点も変えない
    let plain = sampled(&xs, 1.0e4, 300.0, &sigma0);
    assert_eq!(plain.floored, 0);
}

#[test]
fn empty_sampling_is_an_error() {
    // ラダーかエネルギーの数が 0 なら、平均が 0/0 になるので、値を返さずに止める
    let urr = s_wave(
        vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
        Interpolation::LinLin,
    );
    let xs = UnresolvedXs::new(&urr).unwrap();
    for sampling in [
        Sampling { ladders: 0, ..TEST },
        Sampling { samples: 0, ..TEST },
    ] {
        let result = xs.shielded(1.0e4, 300.0, &[1.0e10], Sigma::default(), sampling, 7);
        assert!(
            matches!(result, Err(UrrError::EmptySampling { .. })),
            "{result:?}"
        );
    }
}
