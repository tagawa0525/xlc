use std::f64::consts::PI;

use super::*;

/// MatLib-1.2.2（PNNL-39258）の式 (2-1)、(2-2)、(2-4) の値（同じ式を別に計算したもの）
#[test]
fn uo2_conductivity_follows_matlib() {
    for (t, burnup, d, k) in [
        (1000.0, 0.0, 0.95, 3.434189679430292),
        (900.0, 30.0, 10.257 / 10.96, 2.597404435167642),
        (1500.0, 60.0, 0.95, 1.8701626340993773),
    ] {
        let value = uo2_conductivity(t, burnup, d);
        assert!(
            (value / k - 1.0).abs() < 1e-14,
            "{t} K, {burnup}: {value} vs {k}"
        );
    }
}

/// MatLib-1.2.2 の式 (3-1)、(3-2)
#[test]
fn zircaloy_conductivity_follows_matlib() {
    assert!((zircaloy_conductivity(600.0) - 16.475288).abs() < 1e-12);
    assert_eq!(zircaloy_conductivity(2098.0), 36.0);
}

/// ∫_{T₀}^{T} (a + bT′) dT′ = I の解（2 次方程式の正の根）
#[test]
fn inverse_kirchhoff_integrals_solve_linear_conductivities() {
    let (a, b) = (2.0, 3.0e-3);
    let k = |t: f64| a + b * t;
    let (from, integral) = (600.0, 4000.0);
    let t = inverse_kirchhoff(k, from, integral).unwrap();
    // b/2 (T² − T₀²) + a (T − T₀) = I
    let exact = (-a + (a * a + 2.0 * b * (integral + a * from + b * from * from / 2.0)).sqrt()) / b;
    assert!((t - exact).abs() < 1e-9, "{t} vs {exact}");
    assert_eq!(inverse_kirchhoff(k, from, 0.0).unwrap(), from);
    assert!(inverse_kirchhoff(k, from, -1.0).is_err());
    assert!(inverse_kirchhoff(|_| 0.0, from, 1.0).is_err());
}

/// 一定の熱伝導度の一様な発熱の円柱: 中心と表面の差は q′/(4πk)、体積平均は中心と表面の平均
#[test]
fn constant_conductivity_pellets_follow_the_parabola() {
    let (k, q, surface) = (3.0, 20_000.0, 700.0);
    let p = pellet(|_| k, surface, q).unwrap();
    let rise = q / (4.0 * PI * k);
    assert!((p.centre - surface - rise).abs() < 1e-9);
    assert!((p.average - surface - rise / 2.0).abs() < 1e-9);
}

/// Rowlands の実効温度は (4/9) T_c + (5/9) T_s、体積平均はペレットの平均
#[test]
fn effective_temperatures_weigh_the_centre_and_the_surface() {
    let p = Pellet {
        centre: 1300.0,
        surface: 800.0,
        average: 1000.0,
    };
    assert_eq!(Effective::VolumeAverage.of(&p), 1000.0);
    assert!((Effective::Rowlands.of(&p) - (4.0 * 1300.0 + 5.0 * 800.0) / 9.0).abs() < 1e-12);
}

/// 壁の温度は Dittus–Boelter と Jens–Lottes（FAST-1.2.2 の式 (4-15)）の小さい方
#[test]
fn walls_take_the_cooler_of_convection_and_nucleate_boiling() {
    let (p, t_sat) = (15.5, 617.9);
    let (h, q) = (40_000.0, 1.0e6);
    // 単相: T_b + q″/h = 590 + 25 は飽和温度より低い
    assert_eq!(wall_temperature(590.0, q, h, p, t_sat), 615.0);
    // 沸騰: T_sat + (q″/C)^(1/4) e^(−p/6.2) は 610 + 25 より低い
    let boiling = t_sat + (q / JENS_LOTTES).powf(0.25) * (-p / 6.2).exp();
    assert!(boiling < 635.0);
    assert_eq!(wall_temperature(610.0, q, h, p, t_sat), boiling);
}

/// Dittus–Boelter の h = 0.023 (k/D_e) Re^0.8 Pr^0.4 と、層流の Nu = 7.86 の大きい方（FAST-1.2.2 の式 (4-11)〜(4-13)）
#[test]
fn film_coefficients_follow_dittus_boelter() {
    let (p, t, g, de) = (15.5, 580.0, 3500.0, 0.0118);
    let w = water::region_one(p, t).unwrap();
    let mu = water::viscosity(t, w.density);
    let k = water::conductivity(p, t).unwrap();
    let re = g * de / mu;
    let pr = w.cp * 1e3 * mu / k;
    let expected = 0.023 * k / de * re.powf(0.8) * pr.powf(0.4);
    let h = film_coefficient(p, t, g, de).unwrap();
    assert!((h / expected - 1.0).abs() < 1e-14, "{h} vs {expected}");
    // 流れの止まった極限は層流
    let laminar = film_coefficient(p, t, 1e-6, de).unwrap();
    assert!((laminar / (7.86 * k / de) - 1.0).abs() < 1e-14);
}

/// 燃料棒の温度は外から順に上がり、熱流束の式を満たす
#[test]
fn rods_heat_up_from_the_coolant_inwards() {
    let rod = Rod {
        pellet_radius: 0.4096,
        clad_inner_radius: 0.418,
        clad_outer_radius: 0.475,
        density_fraction: 10.257 / 10.96,
        gap_conductance: 1.0e4,
    };
    let (p, bulk, g, de, q) = (15.513, 580.0, 3485.0, 0.0111795, 17_800.0);
    let t = rod
        .temperatures(q, bulk, p, g, de, 0.0, Effective::VolumeAverage)
        .unwrap();
    assert!(bulk < t.clad_outer && t.clad_outer < t.clad_inner);
    assert!(t.clad_inner < t.pellet.surface && t.pellet.surface < t.pellet.average);
    assert!(t.pellet.average < t.pellet.centre);
    assert_eq!(t.effective, t.pellet.average);
    // ギャップ: q′ = 2π r_p h_gap (T_s − T_ci)
    let gap = 2.0 * PI * rod.pellet_radius * 1e-2 * rod.gap_conductance;
    assert!((gap * (t.pellet.surface - t.clad_inner) / q - 1.0).abs() < 1e-12);
    // 被覆管: ∫ k dT = q′ ln(r_co / r_ci) / (2π)
    let clad = kirchhoff(zircaloy_conductivity, t.clad_outer, t.clad_inner);
    let expected = q * (rod.clad_outer_radius / rod.clad_inner_radius).ln() / (2.0 * PI);
    assert!((clad / expected - 1.0).abs() < 1e-10);
    // 線出力 0 は冷却材の温度
    let cold = rod
        .temperatures(0.0, bulk, p, g, de, 0.0, Effective::Rowlands)
        .unwrap();
    assert_eq!(cold.effective, bulk);
    assert!(
        rod.temperatures(-1.0, bulk, p, g, de, 0.0, Effective::Rowlands)
            .is_err()
    );
}
