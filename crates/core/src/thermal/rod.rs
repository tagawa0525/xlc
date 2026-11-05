//! 燃料棒の径方向の 1 次元の定常の熱伝導（ADR-0229）。外から、冷却材の膜、被覆管、ギャップ、一様な発熱のペレット。
//! 被覆管とペレットは Kirchhoff の積分 ∫ k dT で解くので、径方向のメッシュを持たない。長さの入力は cm、線出力は W/m、
//! 熱伝導度は W/(m K)、熱伝達率は W/(m² K)。

use std::f64::consts::PI;
use std::sync::OnceLock;

use super::water;

/// Jens–Lottes の係数 (W/(m² K⁴))（FAST-1.2.2、PNNL-39254 の式 (4-15)）
pub const JENS_LOTTES: f64 = 2.553_509;
/// UO2 の理論密度 (g/cc)（MatLib-1.2.2、PNNL-39258 の 2.1.6 節）
pub const UO2_THEORETICAL_DENSITY: f64 = 10.96;

/// UO2 の熱伝導度 (W/(m K))。MatLib-1.2.2 の式 (2-1)、(2-2)（ガドリニアなし）と Lucuta の密度の補正の式 (2-4)。`burnup` は GWd/MTU
/// （= MWd/kg）、`density_fraction` は理論密度に対する比
pub fn uo2_conductivity(t: f64, burnup: f64, density_fraction: f64) -> f64 {
    let h = 1.0 / (1.0 + 396.0 * (-6380.0 / t).exp());
    let f = 0.00187 * burnup;
    let g = 0.038 * burnup.powf(0.28);
    let k95 = 1.0 / (0.0452 + 2.46e-4 * t + f + (1.0 - 0.9 * (-0.04 * burnup).exp()) * g * h)
        + 3.5e9 / (t * t) * (-16361.0 / t).exp();
    1.0789 * k95 * density_fraction / (1.0 + 0.5 * (1.0 - density_fraction))
}

/// Zircaloy などのジルコニウムの合金の被覆管の熱伝導度 (W/(m K))（MatLib-1.2.2 の式 (3-1)、2098 K 以上は式 (3-2)）
pub fn zircaloy_conductivity(t: f64) -> f64 {
    if t >= 2098.0 {
        return 36.0;
    }
    7.511 + 2.088e-2 * t - 1.45e-5 * t * t + 7.668e-9 * t * t * t
}

/// 5 点の Gauss–Legendre の [−1, 1] の点と重み
const GAUSS5: [(f64, f64); 5] = [
    (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
    (-0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
    (0.0, 0.568_888_888_888_888_9),
    (0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
    (0.906_179_845_938_664, 0.236_926_885_056_189_1),
];
/// 熱伝導度の積分の区間の幅の上限 (K)。熱伝導度は温度で滑らかなので、5 点の求積の誤差は 1e-12 より小さい
const PANEL: f64 = 50.0;

/// ∫_{a}^{b} k(T) dT（区間を `PANEL` 以下に分けた 5 点の Gauss–Legendre）
pub fn kirchhoff(k: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let panels = ((b - a).abs() / PANEL).ceil().max(1.0);
    let width = (b - a) / panels;
    (0..panels as usize)
        .map(|i| {
            let mid = a + (i as f64 + 0.5) * width;
            GAUSS5
                .iter()
                .map(|&(x, w)| w * k(mid + x * width / 2.0))
                .sum::<f64>()
                * width
                / 2.0
        })
        .sum()
}

/// ∫_{from}^{T} k dT = `integral`（≥ 0）の T（Newton 法）。熱伝導度が正でなければエラー
pub fn inverse_kirchhoff(k: impl Fn(f64) -> f64, from: f64, integral: f64) -> Result<f64, String> {
    if !(integral >= 0.0 && integral.is_finite()) {
        return Err(format!(
            "the conduction integral {integral} must be non-negative"
        ));
    }
    if integral == 0.0 {
        return Ok(from);
    }
    let conductivity = |t: f64| {
        let value = k(t);
        if value > 0.0 && value.is_finite() {
            Ok(value)
        } else {
            Err(format!("the conductivity at {t} K is {value}"))
        }
    };
    let mut t = from + integral / conductivity(from)?;
    for _ in 0..100 {
        let step = (kirchhoff(&k, from, t) - integral) / conductivity(t)?;
        t -= step;
        if step.abs() <= 1e-12 * t {
            return Ok(t);
        }
    }
    Err(format!(
        "the temperature for the conduction integral {integral} from {from} K did not converge"
    ))
}

/// ペレットの温度 (K)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pellet {
    pub centre: f64,
    pub surface: f64,
    /// 体積平均
    pub average: f64,
}

/// ペレットの体積平均の求積の点の数（u = 1 − r²/R² の [0, 1] の Gauss–Legendre）。T(u) は u で滑らかで、一定の熱伝導度では
/// 線形
const PELLET_POINTS: usize = 16;

/// u の [0, 1] の求積の点（昇順）と重み
fn pellet_points() -> &'static [(f64, f64)] {
    static POINTS: OnceLock<Vec<(f64, f64)>> = OnceLock::new();
    POINTS.get_or_init(|| {
        let mut points: Vec<(f64, f64)> = xlc_numeric::quadrature::gauss_legendre(PELLET_POINTS)
            .unwrap_or_default()
            .into_iter()
            .map(|(x, w)| ((x + 1.0) / 2.0, w / 2.0))
            .collect();
        points.sort_by(|a, b| a.0.total_cmp(&b.0));
        points
    })
}

/// 一様な発熱の線出力 `q` (W/m) のペレットの温度。表面の温度 `surface` から ∫_{T_s}^{T(u)} k dT = q u / (4π) を u の昇順に解く
pub fn pellet(k: impl Fn(f64) -> f64, surface: f64, q: f64) -> Result<Pellet, String> {
    let points = pellet_points();
    if points.len() != PELLET_POINTS {
        return Err("the Gauss-Legendre points for the pellet average are missing".into());
    }
    let scale = q / (4.0 * PI);
    let (mut t, mut u, mut average) = (surface, 0.0, 0.0);
    for &(x, w) in points {
        t = inverse_kirchhoff(&k, t, scale * (x - u))?;
        u = x;
        average += w * t;
    }
    let centre = inverse_kirchhoff(&k, t, scale * (1.0 - u))?;
    Ok(Pellet {
        centre,
        surface,
        average,
    })
}

/// テーブルの燃料の温度の軸に渡す実効温度（ADR-0229 の案 E1、E2。既定なし）
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Effective {
    /// ペレットの体積平均
    VolumeAverage,
    /// Rowlands: (4/9) T_c + (5/9) T_s
    Rowlands,
}

impl Effective {
    pub fn of(self, p: &Pellet) -> f64 {
        match self {
            Self::VolumeAverage => p.average,
            Self::Rowlands => (4.0 * p.centre + 5.0 * p.surface) / 9.0,
        }
    }
}

/// 表面の熱伝達率 (W/(m² K))。Dittus–Boelter と層流の Nu = 7.86 の大きい方（FAST-1.2.2 の式 (4-11)〜(4-13)）。物性は冷却材の
/// 温度 `t` (K) と圧力 `p` (MPa)、`mass_flux` は kg/(m² s)、`diameter` は水力直径 (m)
pub fn film_coefficient(p: f64, t: f64, mass_flux: f64, diameter: f64) -> Result<f64, String> {
    let w = water::region_one(p, t)?;
    let mu = water::viscosity(t, w.density);
    let k = water::conductivity(p, t)?;
    let re = mass_flux * diameter / mu;
    let pr = w.cp * 1e3 * mu / k;
    let turbulent = 0.023 * k / diameter * re.powf(0.8) * pr.powf(0.4);
    Ok(turbulent.max(7.86 * k / diameter))
}

/// 被覆管の外面の温度 (K): 単相の T_b + q″/h と、核沸騰の Jens–Lottes の T_sat + (q″/C)^{1/4} e^{−p/6.2 MPa}（FAST-1.2.2 の
/// 式 (4-15)）の小さい方（ADR-0229 の案 F1。FAST の表 2-19 の切り替えは壁の温度が跳ぶので使わない）
pub fn wall_temperature(bulk: f64, heat_flux: f64, h: f64, p: f64, saturation: f64) -> f64 {
    let convection = bulk + heat_flux / h;
    let boiling = saturation + (heat_flux / JENS_LOTTES).powf(0.25) * (-p / 6.2).exp();
    convection.min(boiling)
}

/// 燃料棒の寸法と物性の入力（`[thermal-hydraulics] rod`）
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rod {
    /// cm
    pub pellet_radius: f64,
    /// cm
    pub clad_inner_radius: f64,
    /// cm
    pub clad_outer_radius: f64,
    /// ペレットの密度の理論密度に対する比
    pub density_fraction: f64,
    /// W/(m² K)（ペレットの外面の熱流束に掛ける）
    pub gap_conductance: f64,
}

/// 燃料棒の温度 (K)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RodTemperatures {
    pub clad_outer: f64,
    pub clad_inner: f64,
    pub pellet: Pellet,
    /// テーブルの燃料の温度の軸に渡す温度
    pub effective: f64,
}

impl Rod {
    /// 線出力 `q` (W/m)、冷却材の温度 `bulk` (K)、圧力 `p` (MPa)、質量流束 `mass_flux` (kg/(m² s))、水力直径 `diameter` (m)、
    /// 燃焼度 `burnup` (MWd/kg) の燃料棒の温度
    #[allow(
        clippy::too_many_arguments,
        reason = "線出力、冷却材の状態、流路、燃焼度、実効温度の選び方は出どころが別で、まとめる型を作るより読みやすい"
    )]
    pub fn temperatures(
        &self,
        q: f64,
        bulk: f64,
        p: f64,
        mass_flux: f64,
        diameter: f64,
        burnup: f64,
        effective: Effective,
    ) -> Result<RodTemperatures, String> {
        if !(q >= 0.0 && q.is_finite()) {
            return Err(format!("the linear power {q} W/m must be non-negative"));
        }
        if q == 0.0 {
            let pellet = Pellet {
                centre: bulk,
                surface: bulk,
                average: bulk,
            };
            return Ok(RodTemperatures {
                clad_outer: bulk,
                clad_inner: bulk,
                pellet,
                effective: bulk,
            });
        }
        let (rp, rci, rco) = (
            self.pellet_radius * 1e-2,
            self.clad_inner_radius * 1e-2,
            self.clad_outer_radius * 1e-2,
        );
        let h = film_coefficient(p, bulk, mass_flux, diameter)?;
        let clad_outer = wall_temperature(
            bulk,
            q / (2.0 * PI * rco),
            h,
            p,
            water::saturation_temperature(p),
        );
        let clad_inner = inverse_kirchhoff(
            zircaloy_conductivity,
            clad_outer,
            q * (rco / rci).ln() / (2.0 * PI),
        )?;
        let surface = clad_inner + q / (2.0 * PI * rp * self.gap_conductance);
        let pellet = pellet(
            |t| uo2_conductivity(t, burnup, self.density_fraction),
            surface,
            q,
        )?;
        Ok(RodTemperatures {
            clad_outer,
            clad_inner,
            pellet,
            effective: effective.of(&pellet),
        })
    }
}

#[cfg(test)]
mod tests;
