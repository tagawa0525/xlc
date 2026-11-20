use super::*;
use crate::endf::IncidentLaw;

fn assert_close(actual: &[f64], expected: &[f64]) {
    for (a, e) in actual.iter().zip(expected) {
        assert!((a - e).abs() <= 1e-14, "{actual:?} vs {expected:?}");
    }
}

#[test]
fn evaporation_spectrum_is_cut_at_e_minus_u() {
    // E' ≤ E - U = 2 MeV。群ごとの割合は P(2, x) の差（x = E'/θ）
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]).unwrap();
    let law = EnergyLaw::Evaporation {
        restriction: 1.0e6,
        theta: ClampedTab1::new(Tab1::linear(vec![1.0e6, 2.0e7], vec![1.0e5, 1.0e5])).unwrap(),
    };
    let grouped = GroupedLaw::new(&groups, &law);
    let mut out = vec![0.0; 3];
    grouped.accumulate(3.0e6, 2.0, &mut out);
    let p = |x: f64| gamma_p_2(x) / gamma_p_2(20.0);
    assert_close(
        &out,
        &[
            2.0 * (1.0 - p(10.0)),
            2.0 * (p(10.0) - p(1.0)),
            2.0 * p(1.0),
        ],
    );
    // E ≤ U では中性子を出さない
    let mut out = vec![0.0; 3];
    grouped.accumulate(1.0e6, 1.0, &mut out);
    assert_eq!(out, vec![0.0; 3]);
    // E - U が群境界を横切る点で、群ごとの割合の傾きが変わる
    let breakpoints = grouped.breakpoints();
    for e in [1.0e6, 1.1e6, 2.0e6] {
        assert!(breakpoints.contains(&e), "{e}: {breakpoints:?}");
    }
}

#[test]
fn maxwellian_interpolates_theta() {
    // θ(2 MeV) = 200 keV。E - U が群構造の上限を超える分は最上位の群に入る
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]).unwrap();
    let law = EnergyLaw::Maxwellian {
        restriction: -3.0e7,
        theta: ClampedTab1::new(Tab1::linear(vec![1.0e6, 3.0e6], vec![1.0e5, 3.0e5])).unwrap(),
    };
    let grouped = GroupedLaw::new(&groups, &law);
    let mut out = vec![0.0; 3];
    grouped.accumulate(2.0e6, 1.0, &mut out);
    let p = |x: f64| gamma_p_3_2(x) / gamma_p_3_2(160.0);
    assert_close(&out, &[1.0 - p(5.0), p(5.0) - p(0.5), p(0.5)]);
    assert!(grouped.breakpoints().contains(&3.0e6));
}

#[test]
fn integrates_piecewise_linear_density_by_group() {
    let groups = GroupStructure::new(vec![4.0, 2.0, 1.0]).unwrap();
    let density = Tab1::linear(vec![1.0, 3.0], vec![0.5, 0.5]);
    assert_eq!(outgoing_group_integrals(&groups, &density), vec![0.5, 0.5]);
    // 群構造の外に出る分は端の群に入れる（NJOY と同じ。中性子を失わない）
    let density = Tab1::linear(vec![0.0, 3.0, 5.0], vec![0.2, 0.2, 0.2]);
    let r = outgoing_group_integrals(&groups, &density);
    assert!(
        (r[0] - 0.6).abs() < 1e-15 && (r[1] - 0.4).abs() < 1e-15,
        "{r:?}"
    );
}

#[test]
fn linear_density_split_inside_a_segment() {
    let groups = GroupStructure::new(vec![2.0, 1.0, 0.5]).unwrap();
    let density = Tab1::linear(vec![0.0, 2.0], vec![0.0, 1.0]);
    let r = outgoing_group_integrals(&groups, &density);
    // ∫ E/2 dE。0.5 eV 未満の分（0.0625）は最下位の群に入る
    assert!(
        (r[0] - 0.75).abs() < 1e-15 && (r[1] - 0.25).abs() < 1e-15,
        "{r:?}"
    );
}

#[test]
fn center_of_mass_energy_spreads_in_laboratory() {
    // 重心系で E'cm にほぼ単色の中性子は、実験室系で
    // [(√E'cm - √Ec)², (√E'cm + √Ec)²] に一様に分布する（Ec = E/(A+1)²）
    let (e, awr) = (1.6e7, 15.0);
    let ec = e / (awr + 1.0_f64).powi(2);
    assert_eq!(ec, 62500.0);
    let distributions = OutgoingDistributions {
        breaks: vec![1],
        laws: vec![IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: false,
        }],
        incident: vec![e],
        spectra: vec![Tab1 {
            breaks: vec![2],
            laws: vec![Interpolation::Histogram],
            x: vec![1.0e6 - 1.0, 1.0e6 + 1.0],
            y: vec![0.5, 0.5],
        }],
    };
    let groups = GroupStructure::new(vec![2.0e6, 1.0e6, 1.0e5]).unwrap();
    let grouped = GroupedOutgoing::center_of_mass(&groups, &distributions, awr);
    let mut out = vec![0.0; 2];
    grouped.accumulate(e, 1.0, &mut out);
    // 実験室系の範囲は [562500, 1562500]
    assert!(
        (out[0] - 0.5625).abs() < 1e-4 && (out[1] - 0.4375).abs() < 1e-4,
        "{out:?}"
    );
}

/// 重心系の密度 `density`（線形の区間）を実験室系の群に振り分けた値の参照。√E' について細かい中点則で、E' ごとの実験室系の一様な
/// 範囲 [(√E' − √Ec)², (√E' + √Ec)²] と群の重なりを数える（群構造の外は端の群）
fn laboratory_reference(groups: &GroupStructure, density: &Tab1, ec: f64) -> Vec<f64> {
    const STEPS: usize = 400_000;
    let bounds = groups.bounds();
    let ng = groups.len();
    let c = ec.sqrt();
    let (u0, u1) = (density.x[0].sqrt(), density.x[density.x.len() - 1].sqrt());
    let du = (u1 - u0) / STEPS as f64;
    let mut out = vec![0.0; ng];
    for i in 0..STEPS {
        let u = u0 + (i as f64 + 0.5) * du;
        let f = density.eval(u * u).unwrap_or(0.0) * 2.0 * u * du;
        let (lo, hi) = ((u - c).powi(2), (u + c).powi(2));
        for g in 0..ng {
            let top = if g == 0 { f64::INFINITY } else { bounds[g] };
            let bottom = if g == ng - 1 { 0.0 } else { bounds[g + 1] };
            let overlap = hi.min(top) - lo.max(bottom);
            if overlap > 0.0 {
                out[g] += f * overlap / (hi - lo);
            }
        }
    }
    out
}

#[test]
fn center_of_mass_panel_wider_than_ec_is_split_exactly_into_low_groups() {
    // Pa-233 の MF6 MT18 の入射 1e-5 eV の形: E' = 0〜100 keV の 1 区間の線形の密度。Ec = 100 eV の程度の群（10〜1000 eV）にも
    // 実験室系の範囲が掛かる
    let density = Tab1::linear(vec![0.0, 1.0e5], vec![1.0e-6 / 0.55, 1.0e-5 / 0.55]);
    let groups =
        GroupStructure::new(vec![2.0e5, 1.0e4, 1.0e3, 300.0, 100.0, 30.0, 10.0, 1.0e-5]).unwrap();
    let ec = 100.0;
    let actual = laboratory_group_integrals(&groups, &density, ec);
    let expected = laboratory_reference(&groups, &density, ec);
    for (g, (a, e)) in actual.iter().zip(&expected).enumerate() {
        assert!(
            (a - e).abs() <= 1e-7 * e.abs() + 1e-15,
            "group {g}: {a:e} vs {e:e}"
        );
    }
    assert!(
        (actual.iter().sum::<f64>() - 1.0).abs() < 1e-12,
        "{actual:?}"
    );
}

#[test]
fn center_of_mass_groups_out_of_reach_stay_exactly_zero() {
    // 重心系で 1〜3 MeV の 30 区間の密度、Ec = 1 keV。実験室系の範囲 [(√E' − √Ec)², (√E' + √Ec)²] の届かない群は丸めの値でなく
    // ちょうど 0 で、どの群も負にならない（負の χ は媒質の検査で止まる）
    let x: Vec<f64> = (0..=30)
        .map(|i| 1.0e6 + 2.0e6 * f64::from(i) / 30.0)
        .collect();
    let y: Vec<f64> = x.iter().map(|e| 1.0 + (e / 1.0e6).sin()).collect();
    let density = Tab1::linear(x, y);
    let bounds: Vec<f64> = (0..200)
        .map(|i| 2.0e7 * (1.0e-12_f64).powf(f64::from(i) / 199.0))
        .collect();
    let groups = GroupStructure::new(bounds).unwrap();
    let ec = 1.0e3;
    let actual = laboratory_group_integrals(&groups, &density, ec);
    let (low, high) = (
        (1.0e6_f64.sqrt() - ec.sqrt()).powi(2),
        (3.0e6_f64.sqrt() + ec.sqrt()).powi(2),
    );
    for (g, value) in actual.iter().enumerate() {
        let (bottom, top) = groups.range(g);
        assert!(*value >= 0.0, "group {g}: {value:e}");
        if top <= low || bottom >= high {
            assert_eq!(*value, 0.0, "group {g} [{bottom:e}, {top:e}]");
        }
    }
    let total: f64 = density
        .x
        .windows(2)
        .zip(density.y.windows(2))
        .map(|(x, y)| 0.5 * (x[1] - x[0]) * (y[0] + y[1]))
        .sum();
    assert!((actual.iter().sum::<f64>() / total - 1.0).abs() < 1e-13);
}

#[test]
fn center_of_mass_motion_far_smaller_than_the_outgoing_energy_keeps_the_total() {
    // 入射 1e-5 eV の重い標的（A ≈ 232）の Ec ≈ 1.8e-10 eV。実験室系の範囲の幅 4√(E' Ec) は E' の丸めの桁に近く、範囲の両端の差から
    // 幅を求めると相対 1e-8 ほど狂う。群に入る量の和は重心系の密度の積分のまま
    let density = Tab1::linear(vec![0.0, 1.0e6, 2.0e7], vec![1.0e-7, 1.0e-7, 0.0]);
    let groups = GroupStructure::wims69();
    let ec = 1.0e-5 / 233.0_f64.powi(2);
    let actual = laboratory_group_integrals(&groups, &density, ec);
    let total = 0.5 * 1.0e6 * 2.0e-7 + 0.5 * 1.9e7 * 1.0e-7;
    assert!(
        (actual.iter().sum::<f64>() / total - 1.0).abs() < 1e-13,
        "{}",
        actual.iter().sum::<f64>() / total - 1.0
    );
}

/// 入射エネルギー 10 と 20 で、[0, 1] と [0, 3] の一様分布
fn widening(law: IncidentLaw) -> OutgoingDistributions {
    OutgoingDistributions {
        breaks: vec![2],
        laws: vec![law],
        incident: vec![10.0, 20.0],
        spectra: vec![
            Tab1::linear(vec![0.0, 1.0], vec![1.0, 1.0]),
            Tab1::linear(vec![0.0, 3.0], vec![1.0 / 3.0, 1.0 / 3.0]),
        ],
    }
}

#[test]
fn unit_base_keeps_the_shape_between_incident_energies() {
    // 単位基底では、中間の入射エネルギーで範囲が [0, 2] に伸びた一様分布になる
    let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
    let unit_base = IncidentLaw {
        law: Interpolation::LinLin,
        unit_base: true,
    };
    let grouped = GroupedOutgoing::new(&groups, &widening(unit_base));
    let mut out = vec![0.0; 3];
    grouped.accumulate(15.0, 1.0, &mut out);
    let expected = [0.0, 0.5, 0.5];
    for (o, e) in out.iter().zip(expected) {
        assert!((o - e).abs() < 1e-6, "{out:?}");
    }
    // ふつうの線形補間なら、2つの分布の平均（[0, 1] に 2/3、[1, 3] に 1/3）
    let cartesian = IncidentLaw {
        law: Interpolation::LinLin,
        unit_base: false,
    };
    let grouped = GroupedOutgoing::new(&groups, &widening(cartesian));
    let mut out = vec![0.0; 3];
    grouped.accumulate(15.0, 1.0, &mut out);
    let expected = [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0];
    for (o, e) in out.iter().zip(expected) {
        assert!((o - e).abs() < 1e-6, "{out:?}");
    }
}

#[test]
fn unit_base_keeps_histogram_steps() {
    // 出射分布がヒストグラム（段差がある）でも、確率の合計と形を保つ。
    // [0, 0.5] に密度 2 と [0, 1.5] に密度 2/3 の中間は、[0, 1] に密度 1
    let step = |width: f64| Tab1 {
        breaks: vec![3],
        laws: vec![Interpolation::Histogram],
        x: vec![0.0, 0.5 * width, width],
        y: vec![2.0 / width, 0.0, 0.0],
    };
    let distributions = OutgoingDistributions {
        breaks: vec![2],
        laws: vec![IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: true,
        }],
        incident: vec![10.0, 20.0],
        spectra: vec![step(1.0), step(3.0)],
    };
    let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
    let grouped = GroupedOutgoing::new(&groups, &distributions);
    let mut out = vec![0.0; 3];
    grouped.accumulate(15.0, 1.0, &mut out);
    let expected = [0.0, 0.0, 1.0];
    for (o, e) in out.iter().zip(expected) {
        assert!((o - e).abs() < 1e-9, "{out:?}");
    }
}

#[test]
fn unit_base_in_the_center_of_mass_matches_the_laboratory_for_a_heavy_target() {
    // 標的が十分重ければ、重心系と実験室系はほぼ同じ
    let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
    let law = IncidentLaw {
        law: Interpolation::LinLin,
        unit_base: true,
    };
    let lab = GroupedOutgoing::new(&groups, &widening(law));
    let cm = GroupedOutgoing::center_of_mass(&groups, &widening(law), 1.0e9);
    let (mut a, mut b) = (vec![0.0; 3], vec![0.0; 3]);
    lab.accumulate(15.0, 1.0, &mut a);
    cm.accumulate(15.0, 1.0, &mut b);
    for (x, y) in a.iter().zip(&b) {
        assert!((x - y).abs() < 1e-4, "{a:?} vs {b:?}");
    }
}

#[test]
fn incident_interpolation_follows_its_law() {
    // 線形・対数（LinLog）なら、入射エネルギーの対数の中点で2つの分布の平均
    let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
    let law = IncidentLaw {
        law: Interpolation::LinLog,
        unit_base: false,
    };
    let grouped = GroupedOutgoing::new(&groups, &widening(law));
    let mut out = vec![0.0; 3];
    grouped.accumulate((10.0_f64 * 20.0).sqrt(), 1.0, &mut out);
    let expected = [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0];
    for (o, e) in out.iter().zip(expected) {
        assert!((o - e).abs() < 1e-12, "{out:?}");
    }
}

#[test]
fn interpolates_between_incident_energies() {
    let groups = GroupStructure::new(vec![4.0, 2.0, 1.0]).unwrap();
    let distributions = OutgoingDistributions {
        breaks: vec![2],
        laws: vec![IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: false,
        }],
        incident: vec![10.0, 20.0],
        spectra: vec![
            Tab1::linear(vec![2.0, 4.0], vec![0.5, 0.5]),
            Tab1::linear(vec![1.0, 2.0], vec![1.0, 1.0]),
        ],
    };
    let grouped = GroupedOutgoing::new(&groups, &distributions);
    let mut out = vec![0.0; 2];
    grouped.accumulate(12.5, 2.0, &mut out);
    assert_eq!(out, vec![1.5, 0.5]);
    // 範囲外は端の分布
    let mut out = vec![0.0; 2];
    grouped.accumulate(30.0, 1.0, &mut out);
    assert_eq!(out, vec![0.0, 1.0]);
}

/// g(x) がヒストグラムで、x の 0〜1e5、1e5〜1e6、1e6〜2e6 に 0.1、0.45、0.45 の一般の蒸発スペクトル（LF=5）。
fn general_evaporation(restriction: f64, theta: f64) -> EnergyLaw {
    let shape = Tab1 {
        breaks: vec![4],
        laws: vec![crate::endf::Interpolation::Histogram],
        x: vec![0.0, 1.0e5, 1.0e6, 2.0e6],
        y: vec![1.0e-6, 5.0e-7, 4.5e-7, 0.0],
    };
    EnergyLaw::GeneralEvaporation {
        restriction,
        theta,
        shape,
    }
}

#[test]
fn general_evaporation_splits_the_tabulated_shape_by_group() {
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]).unwrap();
    let grouped = GroupedLaw::new(&groups, &general_evaporation(-3.0e7, 1.0));
    let mut out = vec![0.0; 3];
    grouped.accumulate(2.0e6, 2.0, &mut out);
    assert_close(&out, &[0.9, 0.9, 0.2]);
}

#[test]
fn general_evaporation_scales_the_outgoing_energy_by_theta() {
    // θ = 2 なので E' = 2x。群の境界 1e5、1e6 eV は x の 5e4、5e5
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]).unwrap();
    let grouped = GroupedLaw::new(&groups, &general_evaporation(-3.0e7, 2.0));
    let mut out = vec![0.0; 3];
    grouped.accumulate(2.0e6, 1.0, &mut out);
    assert_close(&out, &[0.7, 0.25, 0.05]);
}

#[test]
fn general_evaporation_is_cut_at_e_minus_u() {
    // E - U = 1 MeV で切り、残りの 0.55 で割り直す
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]).unwrap();
    let grouped = GroupedLaw::new(&groups, &general_evaporation(5.0e5, 1.0));
    let mut out = vec![0.0; 3];
    grouped.accumulate(1.5e6, 1.0, &mut out);
    assert_close(&out, &[0.0, 0.45 / 0.55, 0.1 / 0.55]);
    // E - U が g の表の点（θ 倍）を横切る点でも、群の割合の傾きが変わる
    let breakpoints = grouped.breakpoints();
    for e in [5.0e5, 6.0e5, 1.5e6, 2.5e6] {
        assert!(breakpoints.contains(&e), "{e}: {breakpoints:?}");
    }
}

/// ENDF-102 の 5.1.1.6 の g(E', E_F)（T_M = `temperature`）。
fn madland_nix_g(energy: f64, fragment: f64, temperature: f64) -> f64 {
    let u1 = (energy.sqrt() - fragment.sqrt()).powi(2) / temperature;
    let u2 = (energy.sqrt() + fragment.sqrt()).powi(2) / temperature;
    // u^(3/2) E1(u) は u → 0 で 0
    let term = |u: f64| {
        if u == 0.0 {
            0.0
        } else {
            u.powf(1.5) * crate::special::exp_e1(u).unwrap() * (-u).exp()
        }
    };
    let lower_gamma = |u: f64| 0.5 * std::f64::consts::PI.sqrt() * gamma_p_3_2(u);
    (term(u2) - term(u1) + lower_gamma(u2) - lower_gamma(u1))
        / (3.0 * (fragment * temperature).sqrt())
}

/// ENDF-102 の f(E → E') = (g(E', E_F(L)) + g(E', E_F(H)))/2 の [low, high] の積分。E' = s² と置き、s を √E_F で分けて
/// （g の傾きが変わる）細かい区間の Gauss-Legendre で求める。
fn madland_nix_integral(low: f64, high: f64, light: f64, heavy: f64, temperature: f64) -> f64 {
    let density = |e: f64| {
        0.5 * (madland_nix_g(e, light, temperature) + madland_nix_g(e, heavy, temperature))
    };
    let mut cuts = vec![low.sqrt(), high.sqrt()];
    cuts.extend(
        [light.sqrt(), heavy.sqrt()]
            .into_iter()
            .filter(|&s| s > low.sqrt() && s < high.sqrt()),
    );
    cuts.sort_by(f64::total_cmp);
    let mut sum = 0.0;
    for w in cuts.windows(2) {
        const PIECES: usize = 2000;
        let width = (w[1] - w[0]) / PIECES as f64;
        for i in 0..PIECES {
            let (mid, half) = (w[0] + (i as f64 + 0.5) * width, 0.5 * width);
            for &(t, wk) in &xlc_numeric::quadrature::GAUSS_LEGENDRE_8 {
                let s = mid + half * t;
                sum += wk * half * density(s * s) * 2.0 * s;
            }
        }
    }
    sum
}

/// ENDF/B-VII.1 の Am-241 の E_F と、入射エネルギーに線形の T_M。
fn madland_nix() -> EnergyLaw {
    EnergyLaw::MadlandNix {
        light: 1.029979e6,
        heavy: 5.467297e5,
        temperature: ClampedTab1::new(Tab1::linear(vec![1.0e-5, 2.0e7], vec![1.0e6, 1.2e6]))
            .unwrap(),
    }
}

#[test]
fn madland_nix_splits_the_endf_formula_by_group() {
    // 群境界に E_F(H) を含め、E_F の前後の両側（ENDF-102 の 5.4.10 の 3 つの領域）を通す。最上位の群は 2e7 eV より上の分を含む
    let bounds = vec![2.0e7, 1.0e7, 2.0e6, 1.0e6, 5.467297e5, 1.0e5, 1.0e3, 1.0e-5];
    let groups = GroupStructure::new(bounds.clone()).unwrap();
    let grouped = GroupedLaw::new(&groups, &madland_nix());
    for (e, temperature) in [(1.0e-5, 1.0e6), (1.0e7, 1.1e6)] {
        let mut out = vec![0.0; 7];
        grouped.accumulate(e, 2.0, &mut out);
        let (light, heavy) = (1.029979e6, 5.467297e5);
        // 2e8 eV より上は e^(-E'/T_M) で 1e-80 より小さい
        let mut edges = bounds.clone();
        edges[0] = 2.0e8;
        *edges.last_mut().unwrap() = 0.0;
        for g in 0..7 {
            let expected =
                2.0 * madland_nix_integral(edges[g + 1], edges[g], light, heavy, temperature);
            assert!(
                (out[g] - expected).abs() <= 1e-12,
                "E = {e}, group {g}: {} vs {expected}",
                out[g]
            );
        }
        // E' は 0〜∞（ENDF-102 の 5.1.1.6）で、入射エネルギーで切らない
        assert!((out.iter().sum::<f64>() - 2.0).abs() <= 1e-13, "{out:?}");
    }
}

#[test]
fn madland_nix_breaks_at_the_temperature_table() {
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    let breakpoints = GroupedLaw::new(&groups, &madland_nix()).breakpoints();
    assert_eq!(breakpoints, vec![1.0e-5, 2.0e7]);
}

#[test]
fn madland_nix_keeps_the_relative_accuracy_of_low_groups() {
    // E'/T_M が 1e-11 まで小さい群でも、群の割合の相対誤差が小さく、負にならない（累積は E'^(3/2) で 0 に近づく）。
    // 30〜80 eV と 80〜200 eV の群は、級数と閉じた式の切り替え（E_F(H) で 54.7 eV、E_F(L) で 103 eV）をまたぐ
    let bounds = vec![
        2.0e7, 1.0e6, 1.0e4, 1.0e3, 2.0e2, 80.0, 30.0, 1.0, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5,
    ];
    let n = bounds.len() - 1;
    let groups = GroupStructure::new(bounds.clone()).unwrap();
    let grouped = GroupedLaw::new(&groups, &madland_nix());
    let mut out = vec![0.0; n];
    grouped.accumulate(1.0e-5, 1.0, &mut out);
    let mut edges = bounds.clone();
    edges[0] = 2.0e8;
    *edges.last_mut().unwrap() = 0.0;
    for g in 1..n {
        let expected = madland_nix_integral(edges[g + 1], edges[g], 1.029979e6, 5.467297e5, 1.0e6);
        assert!(
            out[g] > 0.0 && (out[g] - expected).abs() <= 1e-9 * expected,
            "group {g}: {} vs {expected}",
            out[g]
        );
    }
}

/// ENDF-102 の 5.1.1.5 の Watt の分布の正規化しない密度 exp(-E'/a) sinh(√(b E')) の [low, high] の積分。E' = s² と置き、
/// 細かい区間の Gauss-Legendre で求める。
fn watt_integral(low: f64, high: f64, a: f64, b: f64) -> f64 {
    let (from, to) = (low.sqrt(), high.sqrt());
    const PIECES: usize = 4000;
    let width = (to - from) / PIECES as f64;
    let mut sum = 0.0;
    for i in 0..PIECES {
        let (mid, half) = (from + (i as f64 + 0.5) * width, 0.5 * width);
        for &(t, wk) in &xlc_numeric::quadrature::GAUSS_LEGENDRE_8 {
            let s = mid + half * t;
            sum += wk * half * (-s * s / a).exp() * (b.sqrt() * s).sinh() * 2.0 * s;
        }
    }
    sum
}

/// ENDF/B-VII.1 の U-233 の 1e-5 eV と 20 MeV の a と b を、入射エネルギーに線形にしたもの。b の表は 1.5 MeV にも点を持つ
fn watt(restriction: f64) -> EnergyLaw {
    EnergyLaw::Watt {
        restriction,
        a: ClampedTab1::new(Tab1::linear(vec![1.0e-5, 2.0e7], vec![9.77e5, 1.06e6])).unwrap(),
        b: ClampedTab1::new(Tab1::linear(
            vec![1.0e-5, 1.5e6, 2.0e7],
            vec![2.546e-6, 2.546e-6, 2.62e-6],
        ))
        .unwrap(),
    }
}

/// `watt` の入射エネルギー e の a と b
fn watt_parameters(e: f64) -> (f64, f64) {
    let a = 9.77e5 + (1.06e6 - 9.77e5) * (e - 1.0e-5) / (2.0e7 - 1.0e-5);
    let b = if e <= 1.5e6 {
        2.546e-6
    } else {
        2.546e-6 + (2.62e-6 - 2.546e-6) * (e - 1.5e6) / (2.0e7 - 1.5e6)
    };
    (a, b)
}

#[test]
fn watt_splits_the_endf_formula_by_group() {
    // E' は 0 ≤ E' ≤ E - U（ENDF-102 の 5.1.1.5）。U = -30 MeV（U-233 の評価）では最上位の群が 2e7 eV より上の分を含み、
    // U = 0 では E - U が群の中で分布を切る
    let bounds = vec![2.0e7, 1.0e7, 2.0e6, 1.0e6, 1.0e5, 1.0e3, 1.0e-5];
    let n = bounds.len() - 1;
    let groups = GroupStructure::new(bounds.clone()).unwrap();
    for (restriction, e) in [(-3.0e7, 1.0e-5), (-3.0e7, 1.0e7), (0.0, 5.0e6)] {
        let grouped = GroupedLaw::new(&groups, &watt(restriction));
        let mut out = vec![0.0; n];
        grouped.accumulate(e, 2.0, &mut out);
        let (a, b) = watt_parameters(e);
        let top = e - restriction;
        let mut edges: Vec<f64> = bounds.iter().map(|&edge| edge.min(top)).collect();
        edges[0] = top;
        *edges.last_mut().unwrap() = 0.0;
        let total = watt_integral(0.0, top, a, b);
        for g in 0..n {
            let expected = 2.0 * watt_integral(edges[g + 1], edges[g], a, b) / total;
            assert!(
                (out[g] - expected).abs() <= 1e-12,
                "U = {restriction}, E = {e}, group {g}: {} vs {expected}",
                out[g]
            );
        }
        assert!((out.iter().sum::<f64>() - 2.0).abs() <= 1e-13, "{out:?}");
    }
}

#[test]
fn watt_breaks_at_both_parameter_tables() {
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    let breakpoints = GroupedLaw::new(&groups, &watt(-3.0e7)).breakpoints();
    for e in [1.0e-5, 1.5e6, 2.0e7] {
        assert!(breakpoints.contains(&e), "{e}: {breakpoints:?}");
    }
}

#[test]
fn watt_keeps_the_relative_accuracy_of_low_groups() {
    // E'/a が 1e-11 まで小さい群でも、群の割合の相対誤差が小さく、負にならない（累積は E'^(3/2) で 0 に近づく）。
    // 5e5〜2e6 eV の群は、E'/a = 1 をまたぐ
    let bounds = vec![
        2.0e7, 2.0e6, 5.0e5, 1.0e4, 1.0e3, 30.0, 1.0, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5,
    ];
    let n = bounds.len() - 1;
    let groups = GroupStructure::new(bounds.clone()).unwrap();
    let grouped = GroupedLaw::new(&groups, &watt(-3.0e7));
    let mut out = vec![0.0; n];
    grouped.accumulate(1.0e-5, 1.0, &mut out);
    let (a, b) = watt_parameters(1.0e-5);
    let mut edges = bounds.clone();
    edges[0] = 3.0e7;
    *edges.last_mut().unwrap() = 0.0;
    let total = watt_integral(0.0, 3.0e7, a, b);
    for g in 1..n {
        let expected = watt_integral(edges[g + 1], edges[g], a, b) / total;
        assert!(
            out[g] > 0.0 && (out[g] - expected).abs() <= 1e-9 * expected,
            "group {g}: {} vs {expected}",
            out[g]
        );
    }
}

#[test]
fn watt_cumulative_keeps_its_accuracy_over_the_read_range_of_ab() {
    // 読む範囲の a b/4 = c（5 以下。ENDF の読み込みが上限を確かめる）で、閉じた式と級数の切り替え（x = E'/a = 1）の前後の累積が、
    // 数値積分と合う（a = 1 として b = 4c）
    for c in [1.0e-6, 1.0e-3, 0.05, 0.3, 0.62, 2.0, 5.0] {
        for x in [1.0e-6, 0.5, 0.999, 1.001, 2.0, 30.0] {
            let expected = watt_integral(0.0, x, 1.0, 4.0 * c) * (-c).exp();
            let actual = watt_cumulative(x, c);
            assert!(
                (actual - expected).abs() <= 1e-12 * expected,
                "c = {c}, x = {x}: {actual} vs {expected}"
            );
        }
    }
}

/// 位相空間分布の実験室系の群ごとの割合を、重心系のエネルギーについての細かい中点則で求めた参照。
/// 重心系の E' = E_max sin²θ ごとに、実験室系の [(√E' − √Ec)², (√E' + √Ec)²] の一様分布を振り分ける
fn phase_space_reference(
    groups: &GroupStructure,
    space: &PhaseSpace,
    emax: f64,
    ec: f64,
) -> Vec<f64> {
    let n = 200_000;
    let h = std::f64::consts::FRAC_PI_2 / n as f64;
    let m = space.exponent();
    let mut out = vec![0.0; groups.len()];
    for i in 0..n {
        let t = (i as f64 + 0.5) * h;
        let weight =
            space.normalization() * 2.0 * t.sin().powi(2) * t.cos().powf(2.0 * m + 1.0) * h;
        let v = emax.sqrt() * t.sin();
        let (lo, hi) = ((v - ec.sqrt()).powi(2), (v + ec.sqrt()).powi(2));
        uniform_fractions(groups, lo, hi, weight, &mut out);
    }
    out
}

#[test]
fn phase_space_spreads_by_the_centre_of_mass_motion() {
    let groups = GroupStructure::new(vec![
        2.0e7, 5.0e6, 3.0e6, 1.0e6, 5.0e5, 1.0e5, 1.0e4, 1.0e-5,
    ])
    .unwrap();
    // 重心の運動 Ec が E_max より小さいとき、大きいとき（実験室系の下限が 0 から離れる）
    for (emax, ec) in [(2.0e6, 4.0e5), (3.0e5, 1.5e6)] {
        for particles in [3, 4, 5] {
            let space = PhaseSpace {
                total_mass: 5.0,
                particles,
            };
            let mut out = vec![0.0; groups.len()];
            phase_space_fractions(&groups, &space, emax, ec, 2.0, &mut out);
            let reference = phase_space_reference(&groups, &space, emax, ec);
            for (g, (a, r)) in out.iter().zip(&reference).enumerate() {
                assert!(
                    (a - 2.0 * r).abs() < 2e-7,
                    "E_max {emax}, Ec {ec}, n {particles}, group {g}: {out:?} vs {reference:?}"
                );
            }
            let total: f64 = out.iter().sum();
            assert!((total - 2.0).abs() < 1e-12, "{total}");
        }
    }
}

#[test]
fn phase_space_below_threshold_moves_with_the_centre_of_mass() {
    // E_max ≤ 0 では、重心系で止まった中性子が重心とともに Ec で進む
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]).unwrap();
    let space = PhaseSpace {
        total_mass: 3.0,
        particles: 3,
    };
    let mut out = vec![0.0; 3];
    phase_space_fractions(&groups, &space, -1.0, 4.0e5, 1.5, &mut out);
    assert_eq!(out, vec![0.0, 1.5, 0.0]);
}
