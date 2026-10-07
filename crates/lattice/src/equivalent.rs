//! Dancoff 係数を合わせた円筒の等価なセル。
//!
//! Dancoff 係数の定義は [`crate::moc::dancoff`] と同じ: 燃料を黒体にし、燃料の外の各リングに、無限に広がれば
//! スカラー中性子束が 1 になる一様な等方線源 Σt（散乱なし）を置く。燃料の吸収率 R と、燃料の表面に部分中性子流
//! 1/4 で入る値 R_iso = 2π r_f / 4 から C = 1 − R / R_iso とする。円筒のセルは外周を白色境界とし、衝突確率法
//! （[`Cylinder`]）で解く。線源が与えられ散乱もないので、R はリングの分け方によらない。

use std::f64::consts::PI;

use crate::cp::{CpError, Cylinder};

/// 黒体の燃料の光学的半径 Σ r。Σr = 1e4 との C の差は 1e-12 以下
const BLACK: f64 = 1.0e5;

/// 等価なセルの外径の、減速材の厚さに対する収束の判定
const TOLERANCE: f64 = 1.0e-12;

/// 合わせられる Dancoff 係数の下限。孤立したセルの [`cell_dancoff`] は丸めで 1e-11 程度ずれるので、それより
/// 十分大きくとる
const MIN_DANCOFF: f64 = 1.0e-9;

/// 減速材の厚さを倍にして Dancoff 係数を下回るまで探す回数の上限
const DOUBLINGS: usize = 64;

/// 等価なセルのエラー。
#[derive(Debug, thiserror::Error)]
pub enum EquivalentError {
    #[error("{totals} cross sections for {rings} rings")]
    Rings { totals: usize, rings: usize },
    #[error("ring {0}: total cross sections outside the fuel must be finite and not negative")]
    Total(usize),
    #[error("the moderator total cross section {0} must be positive and finite")]
    Moderator(f64),
    #[error("Dancoff factor {0} must be between {MIN_DANCOFF} and 1")]
    Range(f64),
    #[error("Dancoff factor {0} is not reached even without moderator")]
    Unreachable(f64),
    #[error("Dancoff factor {0} needs a cell larger than {1} cm")]
    TooSmall(f64, f64),
    #[error(transparent)]
    Cp(#[from] CpError),
}

/// 黒体の燃料（最も内側のリング）の円筒のセルの Dancoff 係数。`radii` はリングの外半径、`totals` はリングごとの
/// 全断面積 (1/cm)（燃料の値は使わない）。
pub fn cell_dancoff(radii: &[f64], totals: &[f64]) -> Result<f64, EquivalentError> {
    if totals.len() != radii.len() {
        return Err(EquivalentError::Rings {
            totals: totals.len(),
            rings: radii.len(),
        });
    }
    if let Some(i) = (1..totals.len()).find(|&i| !(totals[i].is_finite() && totals[i] >= 0.0)) {
        return Err(EquivalentError::Total(i));
    }
    let cell = Cylinder::new(radii.to_vec())?;
    let fuel = radii[0];
    let black = BLACK / fuel;
    let mut sigma = totals.to_vec();
    sigma[0] = black;
    let mut sources = totals.to_vec();
    sources[0] = 0.0;
    let flux = cell.fixed_source(&sigma, &sources)?;
    let absorbed = black * flux[0] * cell.volumes()[0];
    Ok(1.0 - absorbed / (PI * fuel / 2.0))
}

/// Dancoff 係数が `dancoff` になる等価なセルの外半径 (cm)。`radii` は減速材より内側のリング（燃料から被覆管まで）の
/// 外半径、`totals` はそのリングと外側の減速材のリングの全断面積 (1/cm)（燃料の値は使わない）。
pub fn equivalent_radius(
    radii: &[f64],
    totals: &[f64],
    dancoff: f64,
) -> Result<f64, EquivalentError> {
    if totals.len() != radii.len() + 1 {
        return Err(EquivalentError::Rings {
            totals: totals.len(),
            rings: radii.len() + 1,
        });
    }
    let moderator = totals[radii.len()];
    if !(moderator > 0.0 && moderator.is_finite()) {
        return Err(EquivalentError::Moderator(moderator));
    }
    if !(MIN_DANCOFF..1.0).contains(&dancoff) {
        return Err(EquivalentError::Range(dancoff));
    }
    let inner = radii.last().copied().ok_or(EquivalentError::Rings {
        totals: totals.len(),
        rings: 0,
    })?;
    // 減速材の厚さ t の Dancoff 係数（t について減る）
    let at = |t: f64| -> Result<f64, EquivalentError> {
        let mut cell = radii.to_vec();
        cell.push(inner + t);
        cell_dancoff(&cell, totals)
    };
    // 下回る厚さを探してから二分法
    let mut high = inner;
    let mut doublings = 0;
    while at(high)? >= dancoff {
        doublings += 1;
        if doublings > DOUBLINGS {
            return Err(EquivalentError::TooSmall(dancoff, inner + high));
        }
        high *= 2.0;
    }
    let mut low = 0.0;
    let mut reached = false;
    while high - low > TOLERANCE * (inner + high) {
        let mid = 0.5 * (low + high);
        if at(mid)? > dancoff {
            low = mid;
            reached = true;
        } else {
            high = mid;
        }
    }
    if !reached {
        return Err(EquivalentError::Unreachable(dancoff));
    }
    Ok(inner + 0.5 * (low + high))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 燃料・ギャップ（真空）・被覆管の外半径と、その外の減速材を含むリングの全断面積（燃料の値は使わない）
    const INNER: [f64; 3] = [0.4096, 0.418, 0.475];
    const TOTALS: [f64; 4] = [0.0, 0.0, 0.3, 1.2];

    fn cell(outer: f64) -> Vec<f64> {
        let mut radii = INNER.to_vec();
        radii.push(outer);
        radii
    }

    #[test]
    fn isolated_pin_is_not_shadowed() {
        // 減速材が十分に厚ければ、隣のピンの影はない
        let c = cell_dancoff(&cell(30.0), &TOTALS).unwrap();
        assert!(c.abs() < 1e-8, "{c}");
    }

    #[test]
    fn dancoff_falls_as_the_cell_grows() {
        let c = [0.55, 0.711, 0.9].map(|r| cell_dancoff(&cell(r), &TOTALS).unwrap());
        assert!(c[0] > c[1] && c[1] > c[2] && c[2] > 0.0, "{c:?}");
    }

    #[test]
    fn equivalent_radius_recovers_the_cell() {
        let c = cell_dancoff(&cell(0.711), &TOTALS).unwrap();
        let r = equivalent_radius(&INNER, &TOTALS, c).unwrap();
        assert!((r - 0.711).abs() < 1e-9, "{r}");
    }

    #[test]
    fn rejects_unreachable_dancoff() {
        // 減速材をなくしても届かない影、0 以下、1 以上、セルの Dancoff 係数の数値誤差より小さい値
        for c in [0.999, 0.0, -0.1, 1.0, f64::NAN, 1e-12] {
            assert!(equivalent_radius(&INNER, &TOTALS, c).is_err(), "{c}");
        }
        // 減速材の全断面積が 0、リングの数が合わない
        assert!(equivalent_radius(&INNER, &[0.0, 0.0, 0.3, 0.0], 0.3).is_err());
        assert!(cell_dancoff(&cell(0.711), &TOTALS[..3]).is_err());
    }
}
