use super::*;

/// `value` が `expected` と相対 `tolerance` 以内
fn close(value: f64, expected: f64, tolerance: f64) {
    assert!(
        ((value - expected) / expected).abs() <= tolerance,
        "{value} vs {expected} (relative {:e})",
        (value - expected) / expected
    );
}

/// IAPWS R7-97(2012) の表 5（領域 1 の式 (7) の検証の値。9 桁）
#[test]
fn region_one_matches_the_if97_verification_table() {
    // (p MPa, T K, v m³/kg, h kJ/kg, cp kJ/(kg K))
    let table = [
        (3.0, 300.0, 0.100215168e-2, 0.115331273e3, 0.417301218e1),
        (80.0, 300.0, 0.971180894e-3, 0.184142828e3, 0.401008987e1),
        (3.0, 500.0, 0.120241800e-2, 0.975542239e3, 0.465580682e1),
    ];
    for (p, t, v, h, cp) in table {
        let w = region_one(p, t).unwrap();
        close(1.0 / w.density, v, 1e-8);
        close(w.enthalpy, h, 1e-8);
        close(w.cp, cp, 1e-8);
    }
}

/// IAPWS R15-11 の表 7 の領域 1 の 2 点（IF97 の cv と (∂ρ/∂p)_T を含む。9 桁）
#[test]
fn region_one_derivatives_match_the_conductivity_verification_table() {
    // (p MPa, T K, ρ kg/m³, (∂ρ/∂p)_T kg/(m³ MPa), cp, cv kJ/(kg K))
    let table = [
        (
            20.0,
            620.0,
            0.613227777e3,
            0.520937820e1,
            0.763433705e1,
            0.303793441e1,
        ),
        (
            50.0,
            620.0,
            0.699226043e3,
            0.184869007e1,
            0.532047725e1,
            0.291692653e1,
        ),
    ];
    for (p, t, rho, drho, cp, cv) in table {
        let w = region_one(p, t).unwrap();
        close(w.density, rho, 1e-8);
        close(w.drho_dp, drho, 1e-8);
        close(w.cp, cp, 1e-8);
        close(w.cv, cv, 1e-8);
    }
}

/// IAPWS R7-97(2012) の表 35（式 (30)）と表 36（式 (31)）
#[test]
fn saturation_matches_the_if97_verification_tables() {
    for (t, p) in [
        (300.0, 0.353658941e-2),
        (500.0, 0.263889776e1),
        (600.0, 0.123443146e2),
    ] {
        close(saturation_pressure(t), p, 1e-8);
    }
    for (p, t) in [
        (0.1, 0.372755919e3),
        (1.0, 0.453035632e3),
        (10.0, 0.584149488e3),
    ] {
        close(saturation_temperature(p), t, 1e-8);
    }
}

/// 2250 psia の IF97 の密度は、VERA の仕様書（CASL-U-2012-0131-004）の付録 G の表 G-2 の 7 点の密度（g/cc、5 桁）に丸めで一致する
#[test]
fn densities_at_2250_psia_round_to_the_vera_table_g2() {
    let p = 2250.0 * PSI;
    for (t, rho) in [
        (550.0, 0.76972),
        (555.0, 0.76106),
        (560.0, 0.75207),
        (565.0, 0.74271),
        (570.0, 0.73294),
        (575.0, 0.72269),
        (580.0, 0.71190),
    ] {
        let d = region_one(p, t).unwrap().density / 1000.0;
        assert!((d - rho).abs() <= 5e-6, "{t} K: {d} vs {rho}");
    }
}

#[test]
fn temperatures_invert_the_enthalpy() {
    let p = 15.5;
    for t in [300.0, 450.0, 565.0, 600.0, 615.0] {
        let h = region_one(p, t).unwrap().enthalpy;
        let back = temperature(p, h).unwrap();
        assert!((back - t).abs() <= 1e-9, "{t}: {back}");
    }
}

#[test]
fn states_outside_the_liquid_region_are_errors() {
    // 15.5 MPa の飽和温度は約 617.9 K
    assert!(region_one(15.5, 620.0).is_err());
    assert!(region_one(15.5, 272.0).is_err());
    assert!(region_one(101.0, 400.0).is_err());
    assert!(region_one(f64::NAN, 400.0).is_err());
    let h = region_one(15.5, 617.0).unwrap().enthalpy;
    assert!(temperature(15.5, h + 50.0).is_err());
    assert!(temperature(15.5, f64::INFINITY).is_err());
}

/// IAPWS R12-08 の表 4（μ₂ = 1 の式 (10)。工業用の推奨の 3 節。小数 6 桁）
#[test]
fn viscosity_matches_the_r12_verification_table() {
    // (T K, ρ kg/m³, μ μPa s)
    let table = [
        (298.15, 998.0, 889.735100),
        (298.15, 1200.0, 1437.649467),
        (373.15, 1000.0, 307.883622),
        (433.15, 1.0, 14.538324),
        (433.15, 1000.0, 217.685358),
        (873.15, 1.0, 32.619287),
        (873.15, 100.0, 35.802262),
        (873.15, 600.0, 77.430195),
        (1173.15, 1.0, 44.217245),
        (1173.15, 100.0, 47.640433),
        (1173.15, 400.0, 64.154608),
    ];
    // 表は小数 6 桁なので、丸めの半分の 5e-7 μPa s まで
    for (t, rho, mu) in table {
        let value = viscosity(t, rho) * 1e6;
        assert!(
            (value - mu).abs() <= 5e-7,
            "{t} K, {rho} kg/m3: {value} vs {mu}"
        );
    }
}

/// IAPWS R15-11 の表 4（λ₂ = 0 の液の 2 点）と、工業用の推奨の表 7（領域 1 の 2 点。λ₀、λ₁、λ₂ と、ξ、Z(y)、μ）
#[test]
fn conductivity_matches_the_r15_verification_tables() {
    for (t, rho, lambda) in [(298.15, 998.0, 607.712868), (298.15, 1200.0, 799.038144)] {
        close(
            dilute_conductivity(t) * finite_density_conductivity(t, rho) * 1e3,
            lambda,
            1e-8,
        );
    }
    // (p MPa, T K, λ, λ₁, λ₂ mW/(m K), ξ nm, Z(y), μ μPa s)
    let table = [
        (
            20.0,
            620.0,
            0.481485195e3,
            0.966869008e1,
            0.126391714e2,
            0.377694973,
            0.166942638,
            0.709051068e2,
        ),
        (
            50.0,
            620.0,
            0.545038940e3,
            0.111212177e2,
            0.575816285e1,
            0.189692422,
            0.113592223,
            0.841527945e2,
        ),
    ];
    for (p, t, lambda, l1, l2, xi, z, mu) in table {
        let w = region_one(p, t).unwrap();
        close(dilute_conductivity(t) * 1e3, 0.484911627e2, 1e-8);
        close(finite_density_conductivity(t, w.density), l1, 1e-8);
        let critical = critical_enhancement(t, &w);
        close(critical.correlation_length, xi, 1e-7);
        close(critical.z, z, 1e-7);
        close(critical.conductivity * 1e3, l2, 1e-7);
        close(viscosity(t, w.density) * 1e6, mu, 1e-8);
        close(conductivity(p, t).unwrap() * 1e3, lambda, 1e-8);
    }
}
