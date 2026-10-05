//! 共鳴の計算に共通する、透過因子・シフト因子・位相シフト（NJOY facts, facphi）。

use crate::physics::AMASSN;

/// チャネル半径の経験式 a = 0.123 (m AWRI)^(1/3) + 0.08 (10^-12 cm)
pub(super) fn channel_radius(awri: f64) -> f64 {
    0.123 * (AMASSN * awri).cbrt() + 0.08
}

/// 透過因子 P_l(ρ)（NJOY facts）。
pub(super) fn penetrability(l: u32, rho: f64) -> f64 {
    let r2 = rho * rho;
    let r4 = r2 * r2;
    match l {
        0 => rho,
        1 => r2 * rho / (1.0 + r2),
        2 => r4 * rho / (9.0 + 3.0 * r2 + r4),
        3 => r4 * r2 * rho / (225.0 + 45.0 * r2 + 6.0 * r4 + r4 * r2),
        4 => r4 * r4 * rho / (11025.0 + 1575.0 * r2 + 135.0 * r4 + 10.0 * r4 * r2 + r4 * r4),
        _ => unreachable!("L={l} is rejected when reading MF2"),
    }
}

/// 硬球散乱の位相シフト φ_l(ρ)（NJOY facphi）。
pub(super) fn phase_shift(l: u32, rho: f64) -> f64 {
    let r2 = rho * rho;
    let phi = match l {
        0 => return rho,
        1 => return rho - rho.atan(),
        2 => rho - (3.0 * rho / (3.0 - r2)).atan(),
        3 => rho - ((15.0 * rho - rho * r2) / (15.0 - 6.0 * r2)).atan(),
        _ => {
            let top = 105.0 * rho - 10.0 * r2 * rho;
            let bot = 105.0 - 45.0 * r2 + r2 * r2;
            rho - (top / bot).atan()
        }
    };
    if phi / rho < 1.0e-6 { 0.0 } else { phi }
}

/// シフト因子 S_l(ρ)（NJOY facts）。
pub(super) fn shift(l: u32, rho: f64) -> f64 {
    let r2 = rho * rho;
    let r4 = r2 * r2;
    let r6 = r4 * r2;
    match l {
        0 => 0.0,
        1 => -1.0 / (1.0 + r2),
        2 => -(18.0 + 3.0 * r2) / (9.0 + 3.0 * r2 + r4),
        3 => -(675.0 + 90.0 * r2 + 6.0 * r4) / (225.0 + 45.0 * r2 + 6.0 * r4 + r6),
        4 => {
            -(44100.0 + 4725.0 * r2 + 270.0 * r4 + 10.0 * r6)
                / (11025.0 + 1575.0 * r2 + 135.0 * r4 + 10.0 * r6 + r4 * r4)
        }
        _ => unreachable!("L={l} is rejected when reading MF2"),
    }
}
