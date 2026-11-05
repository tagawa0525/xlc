//! 圧縮水の物性（ADR-0229）: IAPWS-IF97 の領域 1 と飽和線（IAPWS R7-97(2012)）、IAPWS 2008 の粘性（R12-08）と IAPWS 2011 の
//! 熱伝導率（R15-11）の工業用の推奨。圧力は MPa、温度は K、密度は kg/m³、エンタルピーは kJ/kg、比熱は kJ/(kg K)。

/// 1 psi (MPa)
pub const PSI: f64 = 6.894_757_293_168e-3;

/// IF97 の比気体定数 (kJ/(kg K))、式 (1)
const R: f64 = 0.461_526;

/// 領域 1 の式 (7) の係数（表 2）: (I, J, n)
const REGION_ONE: [(i32, i32, f64); 34] = [
    (0, -2, 0.146_329_712_131_67),
    (0, -1, -0.845_481_871_691_14),
    (0, 0, -0.375_636_036_720_40e1),
    (0, 1, 0.338_551_691_683_85e1),
    (0, 2, -0.957_919_633_878_72),
    (0, 3, 0.157_720_385_132_28),
    (0, 4, -0.166_164_171_995_01e-1),
    (0, 5, 0.812_146_299_835_68e-3),
    (1, -9, 0.283_190_801_238_04e-3),
    (1, -7, -0.607_063_015_658_74e-3),
    (1, -1, -0.189_900_682_184_19e-1),
    (1, 0, -0.325_297_487_705_05e-1),
    (1, 1, -0.218_417_171_754_14e-1),
    (1, 3, -0.528_383_579_699_30e-4),
    (2, -3, -0.471_843_210_732_67e-3),
    (2, 0, -0.300_017_807_930_26e-3),
    (2, 1, 0.476_613_939_069_87e-4),
    (2, 3, -0.441_418_453_308_46e-5),
    (2, 17, -0.726_949_962_975_94e-15),
    (3, -4, -0.316_796_448_450_54e-4),
    (3, 0, -0.282_707_979_853_12e-5),
    (3, 6, -0.852_051_281_201_03e-9),
    (4, -5, -0.224_252_819_080_00e-5),
    (4, -2, -0.651_712_228_956_01e-6),
    (4, 10, -0.143_417_299_379_24e-12),
    (5, -8, -0.405_169_968_601_17e-6),
    (8, -11, -0.127_343_017_416_41e-8),
    (8, -6, -0.174_248_712_306_34e-9),
    (21, -29, -0.687_621_312_955_31e-18),
    (23, -31, 0.144_783_078_285_21e-19),
    (29, -38, 0.263_357_816_627_95e-22),
    (30, -39, -0.119_476_226_400_71e-22),
    (31, -40, 0.182_280_945_814_04e-23),
    (32, -41, -0.935_370_872_924_58e-25),
];

/// 飽和線の式 (29)〜(31) の係数 n₁〜n₁₀（表 34）
const SATURATION: [f64; 10] = [
    0.116_705_214_527_67e4,
    -0.724_213_167_032_06e6,
    -0.170_738_469_400_92e2,
    0.120_208_247_024_70e5,
    -0.323_255_503_223_33e7,
    0.149_151_086_135_30e2,
    -0.482_326_573_615_91e4,
    0.405_113_405_420_57e6,
    -0.238_555_575_678_49,
    0.650_175_348_447_98e3,
];

/// 領域 1 の温度の上限 (K)（領域 3 との境）と下限
const REGION_ONE_MAX_T: f64 = 623.15;
const REGION_ONE_MIN_T: f64 = 273.15;
/// 領域 1 の圧力の上限 (MPa)
const REGION_ONE_MAX_P: f64 = 100.0;

/// 領域 1 の状態の物性
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Liquid {
    /// kg/m³
    pub density: f64,
    /// kJ/kg
    pub enthalpy: f64,
    /// kJ/(kg K)
    pub cp: f64,
    /// kJ/(kg K)
    pub cv: f64,
    /// (∂ρ/∂p)_T (kg/(m³ MPa))
    pub drho_dp: f64,
}

/// 領域 1 の (p, T) の物性（式 (7) と表 3、4）。領域 1 の外（飽和圧力より低い圧力、温度と圧力の範囲の外、有限でない値）はエラー
pub fn region_one(p: f64, t: f64) -> Result<Liquid, String> {
    if !(p.is_finite() && t.is_finite()) {
        return Err(format!("water at {p} MPa and {t} K is not a finite state"));
    }
    if !(REGION_ONE_MIN_T..=REGION_ONE_MAX_T).contains(&t) || !(p > 0.0 && p <= REGION_ONE_MAX_P) {
        return Err(format!(
            "water at {p} MPa and {t} K is outside IAPWS-IF97 region 1 ({REGION_ONE_MIN_T}..{REGION_ONE_MAX_T} K, up to \
             {REGION_ONE_MAX_P} MPa)"
        ));
    }
    if p < saturation_pressure(t) {
        return Err(format!(
            "water at {p} MPa and {t} K is not liquid: the saturation pressure is {} MPa (boiling is not modelled. \
             ADR-0229)",
            saturation_pressure(t)
        ));
    }
    let pi = p / 16.53;
    let tau = 1386.0 / t;
    let (a, b) = (7.1 - pi, tau - 1.222);
    let (mut gp, mut gpp, mut gt, mut gtt, mut gpt) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for &(i, j, n) in &REGION_ONE {
        let (fi, fj) = (f64::from(i), f64::from(j));
        // a^(I−2) と b^(J−2) から上の冪を作る（I = 0 の a^(−1) は I を掛けて消えるので有限のまま）
        let ai2 = a.powi(i - 2);
        let bj2 = b.powi(j - 2);
        let (ai1, ai) = (ai2 * a, ai2 * a * a);
        let (bj1, bj) = (bj2 * b, bj2 * b * b);
        gp -= n * fi * ai1 * bj;
        gpp += n * fi * (fi - 1.0) * ai2 * bj;
        gt += n * ai * fj * bj1;
        gtt += n * ai * fj * (fj - 1.0) * bj2;
        gpt -= n * fi * ai1 * fj * bj1;
    }
    // v = π γ_π R T / p（R は kJ/(kg K)、p は kPa）
    let v = pi * gp * R * t / (p * 1e3);
    let cp = -tau * tau * gtt * R;
    let cv = R * (-tau * tau * gtt + (gp - tau * gpt).powi(2) / gpp);
    // (∂v/∂p)_T = R T γ_ππ / p*²（p* = 16.53 MPa を kPa にし、MPa あたりに直す）
    let dv_dp = R * t * gpp / (16.53e3 * 16.53);
    Ok(Liquid {
        density: 1.0 / v,
        enthalpy: tau * gt * R * t,
        cp,
        cv,
        drho_dp: -dv_dp / (v * v),
    })
}

/// 飽和圧力 (MPa)（式 (30)）
pub fn saturation_pressure(t: f64) -> f64 {
    let n = &SATURATION;
    let theta = t + n[8] / (t - n[9]);
    let a = theta * theta + n[0] * theta + n[1];
    let b = n[2] * theta * theta + n[3] * theta + n[4];
    let c = n[5] * theta * theta + n[6] * theta + n[7];
    (2.0 * c / (-b + (b * b - 4.0 * a * c).sqrt())).powi(4)
}

/// 飽和温度 (K)（式 (31)）
pub fn saturation_temperature(p: f64) -> f64 {
    let n = &SATURATION;
    let beta = p.powf(0.25);
    let e = beta * beta + n[2] * beta + n[5];
    let f = n[0] * beta * beta + n[3] * beta + n[6];
    let g = n[1] * beta * beta + n[4] * beta + n[7];
    let d = 2.0 * g / (-f - (f * f - 4.0 * e * g).sqrt());
    (n[9] + d - ((n[9] + d).powi(2) - 4.0 * (n[8] + n[9] * d)).sqrt()) / 2.0
}

/// エンタルピー `h` (kJ/kg) の領域 1 の温度 (K)。Newton 法で h(p, T) = h を解く（IF97 の逆の式 (11) は式 (7) と 25 mK まで
/// しか合わないので使わない）。飽和温度より上か、領域 1 の外ならエラー
pub fn temperature(p: f64, h: f64) -> Result<f64, String> {
    let top = saturation_temperature(p).min(REGION_ONE_MAX_T);
    let at_top = region_one(p, top)?.enthalpy;
    if !(h.is_finite() && h <= at_top) {
        return Err(format!(
            "water at {p} MPa with the enthalpy {h} kJ/kg is above the liquid at {top} K ({at_top} kJ/kg): boiling is \
             not modelled (ADR-0229)"
        ));
    }
    let mut t = 550.0_f64.min(top);
    for _ in 0..100 {
        let w = region_one(p, t)?;
        let step = (w.enthalpy - h) / w.cp;
        t = (t - step).clamp(REGION_ONE_MIN_T, top);
        if step.abs() <= 1e-11 * t {
            return Ok(t);
        }
    }
    Err(format!(
        "the temperature of water at {p} MPa and {h} kJ/kg did not converge"
    ))
}

/// R12-08 の基準の温度 (K)、密度 (kg/m³)。R15-11 も同じ
const T_STAR: f64 = 647.096;
const RHO_STAR: f64 = 322.0;
/// R15-11 の基準の圧力 (MPa) と比気体定数 (kJ/(kg K))（式 (2)、(6)）
const P_STAR: f64 = 22.064;
const R_CONDUCTIVITY: f64 = 0.461_518_05;

/// 粘性 (Pa s)。R12-08 の式 (10) で臨界の項 μ₂ = 1（3 節の工業用の推奨。ρ は IF97 の密度を渡す）
pub fn viscosity(t: f64, rho: f64) -> f64 {
    const H0: [f64; 4] = [1.677_52, 2.204_62, 0.636_656_4, -0.241_605];
    // 表 2 の H_ij（[i][j]。書いていない係数は 0）
    const H1: [[f64; 7]; 6] = [
        [
            5.200_94e-1,
            2.225_31e-1,
            -2.813_78e-1,
            1.619_13e-1,
            -3.253_72e-2,
            0.0,
            0.0,
        ],
        [
            8.508_95e-2,
            9.991_15e-1,
            -9.068_51e-1,
            2.573_99e-1,
            0.0,
            0.0,
            0.0,
        ],
        [-1.083_74, 1.887_97, -7.724_79e-1, 0.0, 0.0, 0.0, 0.0],
        [
            -2.895_55e-1,
            1.266_13,
            -4.898_37e-1,
            0.0,
            6.984_52e-2,
            0.0,
            -4.356_73e-3,
        ],
        [0.0, 0.0, -2.570_40e-1, 0.0, 0.0, 8.721_02e-3, 0.0],
        [0.0, 1.205_73e-1, 0.0, 0.0, 0.0, 0.0, -5.932_64e-4],
    ];
    let tr = t / T_STAR;
    let dr = rho / RHO_STAR;
    let mu0 = 100.0 * tr.sqrt() / H0.iter().zip(0..).map(|(h, i)| h / tr.powi(i)).sum::<f64>();
    let mu1 = (dr
        * H1.iter()
            .zip(0..)
            .map(|(row, i)| {
                (1.0 / tr - 1.0).powi(i)
                    * row
                        .iter()
                        .zip(0..)
                        .map(|(h, j)| h * (dr - 1.0).powi(j))
                        .sum::<f64>()
            })
            .sum::<f64>())
    .exp();
    mu0 * mu1 * 1e-6
}

/// 希薄気体の極限の熱伝導率 λ₀ (W/(m K))（R15-11 の式 (16)、表 1）
pub fn dilute_conductivity(t: f64) -> f64 {
    const L: [f64; 5] = [
        2.443_221e-3,
        1.323_095e-2,
        6.770_357e-3,
        -3.454_586e-3,
        4.096_266e-4,
    ];
    let tr = t / T_STAR;
    tr.sqrt() / L.iter().zip(0..).map(|(l, k)| l / tr.powi(k)).sum::<f64>() * 1e-3
}

/// 有限の密度の因子 λ₁（無次元。R15-11 の式 (17)、表 2）
pub fn finite_density_conductivity(t: f64, rho: f64) -> f64 {
    const L: [[f64; 6]; 5] = [
        [
            1.603_973_57,
            -0.646_013_523,
            0.111_443_906,
            0.102_997_357,
            -0.050_412_363_4,
            0.006_098_592_58,
        ],
        [
            2.337_718_42,
            -2.788_437_78,
            1.536_161_67,
            -0.463_045_512,
            0.083_282_701_9,
            -0.007_192_012_45,
        ],
        [
            2.196_505_29,
            -4.545_807_85,
            3.557_772_44,
            -1.409_449_78,
            0.275_418_278,
            -0.020_593_881_6,
        ],
        [
            -1.210_513_78,
            1.608_129_89,
            -0.621_178_141,
            0.071_637_322_4,
            0.0,
            0.0,
        ],
        [
            -2.720_337_0,
            4.575_863_31,
            -3.183_692_45,
            1.116_834_8,
            -0.192_683_05,
            0.012_913_842,
        ],
    ];
    let tr = t / T_STAR;
    let dr = rho / RHO_STAR;
    (dr * L
        .iter()
        .zip(0..)
        .map(|(row, i)| {
            (1.0 / tr - 1.0).powi(i)
                * row
                    .iter()
                    .zip(0..)
                    .map(|(l, j)| l * (dr - 1.0).powi(j))
                    .sum::<f64>()
        })
        .sum::<f64>())
    .exp()
}

/// 臨界の項 λ₂ とその途中の値
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Critical {
    /// λ₂ (W/(m K))
    pub conductivity: f64,
    /// 相関長 ξ (nm)
    pub correlation_length: f64,
    /// Z(y)
    pub z: f64,
}

/// 臨界の項 λ₂（R15-11 の式 (18)〜(26)。工業用の推奨の 3.1 節: cp、κ、ζ(T, ρ) は IF97、ζ(T_R, ρ) は式 (25)、μ は μ₂ = 1 の粘性）
pub fn critical_enhancement(t: f64, w: &Liquid) -> Critical {
    // 表 3
    const LAMBDA: f64 = 177.851_4;
    const Q_D_INVERSE: f64 = 0.40;
    const NU: f64 = 0.630;
    const GAMMA: f64 = 1.239;
    const XI0: f64 = 0.13;
    const GAMMA0: f64 = 0.06;
    const T_R: f64 = 1.5;
    // 表 6 の A_ij（[i][j]）と、式 (26) の j の範囲の上端
    const A: [[f64; 5]; 6] = [
        [
            6.537_868_071_995_16,
            6.527_177_592_817_99,
            5.355_005_298_961_24,
            1.552_259_599_066_81,
            1.119_999_264_199_94,
        ],
        [
            -5.611_499_549_233_48,
            -6.308_169_833_875_75,
            -3.964_156_899_254_46,
            0.464_621_290_821_181,
            0.595_748_562_571_649,
        ],
        [
            3.396_241_673_613_25,
            8.083_792_854_925_95,
            8.919_902_089_187_95,
            8.932_373_748_614_79,
            9.889_525_650_789_20,
        ],
        [
            -2.274_926_297_308_78,
            -9.822_405_101_976_03,
            -12.033_872_950_579_0,
            -11.032_196_006_112_6,
            -10.325_505_114_704_0,
        ],
        [
            10.263_185_466_270_9,
            12.135_841_379_139_5,
            9.194_948_651_943_02,
            6.167_809_999_333_60,
            4.668_612_944_574_14,
        ],
        [
            1.978_150_503_315_19,
            -5.543_496_645_712_95,
            -2.168_662_744_797_12,
            -0.965_458_722_086_812,
            -0.503_243_546_373_828,
        ],
    ];
    const BOUNDS: [f64; 4] = [0.310_559_006, 0.776_397_516, 1.242_236_025, 1.863_354_037];
    let tr = t / T_STAR;
    let dr = w.density / RHO_STAR;
    let j = BOUNDS.iter().take_while(|&&b| dr > b).count();
    let zeta_r = 1.0
        / A.iter()
            .zip(0..)
            .map(|(row, i)| row[j] * dr.powi(i))
            .sum::<f64>();
    // 3.4 節の注 2: IF97 の ζ と cp が負か 1e13 を超えたら 1e13 にする
    let cap = |x: f64| if (0.0..=1e13).contains(&x) { x } else { 1e13 };
    let zeta = cap(w.drho_dp * P_STAR / RHO_STAR);
    let cp = cap(w.cp);
    let dchi = (dr * (zeta - zeta_r * T_R / tr)).max(0.0);
    let xi = XI0 * (dchi / GAMMA0).powf(NU / GAMMA);
    let y = xi / Q_D_INVERSE;
    let z = if y < 1.2e-7 {
        0.0
    } else {
        let kappa_inv = w.cv / cp;
        2.0 / (std::f64::consts::PI * y)
            * (((1.0 - kappa_inv) * y.atan() + kappa_inv * y)
                - (1.0 - (-1.0 / (1.0 / y + y * y / (3.0 * dr * dr))).exp()))
    };
    let mu = viscosity(t, w.density) * 1e6;
    Critical {
        conductivity: LAMBDA * dr * (cp / R_CONDUCTIVITY) * tr / mu * z * 1e-3,
        correlation_length: xi,
        z,
    }
}

/// 熱伝導率 (W/(m K))（R15-11 の式 (15) の工業用の推奨。密度と熱力学の量は IF97 の領域 1）
pub fn conductivity(p: f64, t: f64) -> Result<f64, String> {
    let w = region_one(p, t)?;
    Ok(
        dilute_conductivity(t) * finite_density_conductivity(t, w.density)
            + critical_enhancement(t, &w).conductivity,
    )
}

#[cfg(test)]
mod tests;
