//! 燃料を黒体にした固定源の MOC による、ピンごと・群ごとの Dancoff 係数。
//!
//! 燃料を黒体（全断面積が無限大）にし、燃料の外の各領域には、その群の全断面積 Σt と、無限に広がれば
//! スカラー中性子束が 1 になる一様な等方線源 Σt（角度中性子束の平衡値 ψ0 = 1/4π、散乱なし）を置く。
//! 孤立したピンでは、燃料に入る角度中性子束はどの方向でも ψ0 になる。格子では隣のピンの影で減るので、
//! ピンの燃料の吸収率 R と、同じトラックで燃料に入るたびに ψ0 を数えた値 R_iso から C = 1 − R / R_iso とする
//! （同じ求積で数えるので、角度の離散化の誤差は打ち消し合う）。散乱がないので、反復するのは境界の
//! 角度中性子束だけである。

use std::f64::consts::PI;

use super::{POLAR_SIN, Representatives, Sweeper, direction_weight, sum_tallies};
use crate::cmfd::{CmfdError, CoarseMesh};
use crate::material::Material;
use crate::tracking::Tracking;

/// 反復の打ち切り条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DancoffOptions {
    /// 境界の角度中性子束の変化の上限（ψ0 に対する比）。C の誤差もこの程度になる
    pub tolerance: f64,
    pub max_iterations: usize,
}

impl Default for DancoffOptions {
    fn default() -> Self {
        Self {
            tolerance: 1.0e-8,
            max_iterations: 10_000,
        }
    }
}

/// Dancoff 係数の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum DancoffError {
    #[error("{given} values for {regions} regions")]
    Regions { given: usize, regions: usize },
    #[error("no region outside the fuel")]
    NoModerator,
    #[error(
        "region {0}: total cross sections must be positive and finite, or zero in every group (void)"
    )]
    Total(usize),
    #[error("fuel in cell {0} touches the cell boundary; the pin must be surrounded by non-fuel")]
    FuelOnBoundary(usize),
    #[error("cell {cell}: negative Dancoff factor {value} beyond the convergence tolerance")]
    Negative { cell: usize, value: f64 },
    #[error("region {0} is in a lump but not black")]
    NotBlack(usize),
    #[error("lump {lump} is outside the {lumps} lumps")]
    Lump { lump: usize, lumps: usize },
    #[error("lump {0} has no face open to non-black regions")]
    Unexposed(usize),
    #[error("boundary fluxes not converged in {0} iterations")]
    NotConverged(usize),
    #[error(transparent)]
    Cmfd(#[from] CmfdError),
}

/// セルごとの群ごとの Dancoff 係数（燃料のないセルは `None`）。`totals` は [領域][群] の全断面積 (1/cm) で、
/// 燃料の領域（`fuel` が真）の値は使わない。
pub fn dancoff(
    tracking: &Tracking,
    totals: &[Vec<f64>],
    fuel: &[bool],
    options: &DancoffOptions,
) -> Result<Vec<Option<Vec<f64>>>, DancoffError> {
    let lattice = &tracking.lattice;
    let regions = lattice.volumes().len();
    for given in [totals.len(), fuel.len()] {
        if given != regions {
            return Err(DancoffError::Regions { given, regions });
        }
    }
    // 群の数は燃料の外の領域から取る（燃料の行は使わない）
    let ng = (0..regions)
        .find(|&r| !fuel[r])
        .map(|r| totals[r].len())
        .ok_or(DancoffError::NoModerator)?;
    for (region, t) in totals.iter().enumerate() {
        let valid = |v: &f64| *v > 0.0 && v.is_finite();
        let void = t.iter().all(|&v| v == 0.0);
        if !fuel[region] && (t.len() != ng || !(void || t.iter().all(valid))) {
            return Err(DancoffError::Total(region));
        }
    }
    let cells = lattice.cells();
    let cell_of: Vec<usize> = (0..regions).map(|r| lattice.region_cell(r)).collect();

    // 孤立したピンでの吸収率 R_iso / ψ0: 燃料に入るたびに方向の重みを足す
    let mut isolated = vec![0.0; cells];
    for track in &tracking.tracks {
        let weight = track_weight(tracking, track.azimuth);
        let n = track.segments.len();
        for forward in [true, false] {
            let mut previous: Option<usize> = None;
            for step in 0..n {
                let s = if forward { step } else { n - 1 - step };
                let region = track.segments[s].region;
                if fuel[region] {
                    let cell = cell_of[region];
                    // 燃料には同じセルの燃料の外の領域から入ること。格子の端やセルの境目から直接入ると、
                    // 孤立したピンへの入射として数えられない（逆向きのスイープで、出る側も確かめる）
                    match previous {
                        Some(p) if fuel[p] && cell_of[p] == cell => {}
                        Some(p) if !fuel[p] && cell_of[p] == cell => isolated[cell] += weight,
                        _ => return Err(DancoffError::FuelOnBoundary(cell)),
                    }
                }
                previous = Some(region);
            }
        }
    }

    let psi0 = 1.0 / (4.0 * PI);
    let delta = absorbed_by_black(tracking, totals, fuel, ng, options)?;
    let mut absorbed = vec![vec![0.0; ng]; cells];
    for region in (0..regions).filter(|&r| fuel[r]) {
        for (a, d) in absorbed[cell_of[region]]
            .iter_mut()
            .zip(&delta[region * ng..][..ng])
        {
            *a += d;
        }
    }
    absorbed
        .into_iter()
        .zip(&isolated)
        .enumerate()
        .map(|(cell, (rates, &iso))| {
            (iso > 0.0)
                .then(|| {
                    rates
                        .iter()
                        .map(|&r| factor(cell, r, psi0 * iso, options.tolerance))
                        .collect()
                })
                .transpose()
        })
        .collect()
}

/// 吸収率 `absorbed` と孤立したときの入射率 `isolated`（ψ0 を掛けたもの）からの Dancoff 係数 C = 1 − R/R_iso。
/// 孤立した塊では C = 0 が丸めや収束の判定の幅 `tolerance` で負になる。その幅の中なら 0 とする
fn factor(owner: usize, absorbed: f64, isolated: f64, tolerance: f64) -> Result<f64, DancoffError> {
    let c = 1.0 - absorbed / isolated;
    if c >= 0.0 {
        Ok(c)
    } else if c >= -tolerance {
        Ok(0.0)
    } else {
        Err(DancoffError::Negative {
            cell: owner,
            value: c,
        })
    }
}

/// 黒体（`black` が真の領域。全断面積が無限大）の外の各領域に、その群の全断面積 Σt と、無限に広がれば
/// スカラー中性子束が 1 になる一様な等方線源 Σt（ψ0 = 1/4π、散乱なし）を置いた固定源の MOC の、
/// [領域][群] の黒体の吸収率（減衰の重み付き和）。散乱がないので、反復するのは境界の角度中性子束だけ。
fn absorbed_by_black(
    tracking: &Tracking,
    totals: &[Vec<f64>],
    black: &[bool],
    ng: usize,
    options: &DancoffOptions,
) -> Result<Vec<f64>, DancoffError> {
    let regions = black.len();
    let materials: Vec<Material> = (0..regions)
        .map(|region| Material {
            total: if black[region] {
                vec![f64::INFINITY; ng]
            } else {
                totals[region].clone()
            },
            scatter: vec![vec![0.0; ng]; ng],
            anisotropic: Vec::new(),
            fission: vec![0.0; ng],
            nu_fission: vec![0.0; ng],
            fission_terms: Vec::new(),
        })
        .collect();
    let references: Vec<&Material> = materials.iter().collect();
    let mesh = CoarseMesh::new(tracking);
    let sweeper = Sweeper::new(
        tracking,
        &references,
        ng,
        &mesh,
        &[(0, 0)],
        Representatives::all(tracking, &mesh, 1),
    );
    let psi0 = 1.0 / (4.0 * PI);
    let equilibrium: Vec<Vec<f64>> = (0..regions)
        .map(|region| vec![if black[region] { 0.0 } else { psi0 }; ng])
        .collect();
    // 最初は境界から ψ0 が入ってくるとする
    let mut outgoing = sweeper.isotropic_outgoing(&vec![vec![1.0; ng]; regions]);
    let mut next = outgoing.clone();
    let mut tallies = sweeper.tallies(regions);
    for _ in 0..options.max_iterations {
        sweeper.sweep(&equilibrium, &outgoing, &mut next, &mut tallies)?;
        let change = outgoing
            .iter()
            .zip(&next)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        std::mem::swap(&mut outgoing, &mut next);
        if change <= options.tolerance * psi0 {
            return Ok(sum_tallies(&sweeper, &tallies).delta);
        }
    }
    Err(DancoffError::NotConverged(options.max_iterations))
}

/// トラックの1つの向きの、極角についての方向の重みの和 Σ_p 4π (w_a/2) ω_p sinθ_p δ_a
fn track_weight(tracking: &Tracking, azimuth: usize) -> f64 {
    (0..POLAR_SIN.len())
        .map(|p| direction_weight(tracking, azimuth, p))
        .sum()
}

/// 黒体にした塊（燃料でない共鳴の吸収体の領域など）の、群ごとの Dancoff 係数と平均弦長。
#[derive(Debug, Clone, PartialEq)]
pub struct Lump {
    pub dancoff: Vec<f64>,
    /// 平均弦長 l̄ = 4V/S (cm)。S は黒体でない領域に接する面（露出した面）
    pub chord: f64,
}

/// 塊ごとの Dancoff 係数と平均弦長。`black` は [領域] の黒体にするか、`lump` は [領域] の塊の番号（塊でない
/// 領域は `None`。塊は黒体の領域だけ）で、塊は 0..`lumps`。塊でない黒体の領域は、影を作るだけで結果を求めない。
/// 黒体の領域どうしの間を移る中性子は数えず、黒体でない領域から塊に入るたびに、孤立した塊への
/// 入射として数える（隣のセルとつながる支持格子の枠のように、露出した面が一部だけの塊にも使える）。
/// 等方な中性子束 1 での入射率 R_iso = S/4 から l̄ = 4V/S = V/R_iso とする。
pub fn lump_dancoff(
    tracking: &Tracking,
    totals: &[Vec<f64>],
    black: &[bool],
    lump: &[Option<usize>],
    lumps: usize,
    options: &DancoffOptions,
) -> Result<Vec<Lump>, DancoffError> {
    let regions = tracking.lattice.regions();
    for given in [totals.len(), black.len(), lump.len()] {
        if given != regions {
            return Err(DancoffError::Regions { given, regions });
        }
    }
    if let Some(&l) = lump.iter().flatten().find(|&&l| l >= lumps) {
        return Err(DancoffError::Lump { lump: l, lumps });
    }
    if let Some(region) = (0..regions).find(|&r| lump[r].is_some() && !black[r]) {
        return Err(DancoffError::NotBlack(region));
    }
    let ng = (0..regions)
        .find(|&r| !black[r])
        .map(|r| totals[r].len())
        .ok_or(DancoffError::NoModerator)?;
    for (region, t) in totals.iter().enumerate() {
        let valid = |v: &f64| *v > 0.0 && v.is_finite();
        let void = t.iter().all(|&v| v == 0.0);
        if !black[region] && (t.len() != ng || !(void || t.iter().all(valid))) {
            return Err(DancoffError::Total(region));
        }
    }
    // 孤立した塊への入射率 R_iso / ψ0: 黒体でない領域から塊に入るたびに方向の重みを足す
    // （外周から直接入るのは、鏡面反射の向こうの続きの塊からなので数えない）
    let mut isolated = vec![0.0; lumps];
    for track in &tracking.tracks {
        let weight = track_weight(tracking, track.azimuth);
        let n = track.segments.len();
        for forward in [true, false] {
            let mut previous: Option<usize> = None;
            for step in 0..n {
                let s = if forward { step } else { n - 1 - step };
                let region = track.segments[s].region;
                if let (Some(l), Some(p)) = (lump[region], previous)
                    && !black[p]
                {
                    isolated[l] += weight;
                }
                previous = Some(region);
            }
        }
    }
    let psi0 = 1.0 / (4.0 * PI);
    let delta = absorbed_by_black(tracking, totals, black, ng, options)?;
    let volumes = tracking.volumes(regions);
    let mut absorbed = vec![vec![0.0; ng]; lumps];
    let mut volume = vec![0.0; lumps];
    for (region, l) in lump.iter().enumerate() {
        if let Some(l) = *l {
            volume[l] += volumes[region];
            for (a, d) in absorbed[l].iter_mut().zip(&delta[region * ng..][..ng]) {
                *a += d;
            }
        }
    }
    (0..lumps)
        .map(|l| {
            let iso = psi0 * isolated[l];
            if iso <= 0.0 {
                return Err(DancoffError::Unexposed(l));
            }
            Ok(Lump {
                dancoff: absorbed[l]
                    .iter()
                    .map(|&r| factor(l, r, iso, options.tolerance))
                    .collect::<Result<_, _>>()?,
                chord: volume[l] / iso,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests;
