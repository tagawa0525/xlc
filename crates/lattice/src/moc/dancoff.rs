//! 燃料を黒体にした固定源の MOC による、ピンごと・群ごとの Dancoff 係数。
//!
//! 燃料を黒体（全断面積が無限大）にし、燃料の外の各領域には、その群の全断面積 Σt と、無限に広がれば
//! スカラー中性子束が 1 になる一様な等方線源 Σt（角度中性子束の平衡値 ψ0 = 1/4π、散乱なし）を置く。
//! 孤立したピンでは、燃料に入る角度中性子束はどの方向でも ψ0 になる。格子では隣のピンの影で減るので、
//! ピンの燃料の吸収率 R と、同じトラックで燃料に入るたびに ψ0 を数えた値 R_iso から C = 1 − R / R_iso とする
//! （同じ求積で数えるので、角度の離散化の誤差は打ち消し合う）。散乱がないので、反復するのは境界の
//! 角度中性子束だけである。

use std::f64::consts::PI;

use super::{POLAR_SIN, Sweeper, sum_tallies};
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
    #[error("region {0}: total cross sections must be positive and finite")]
    Total(usize),
    #[error("fuel in cell {0} touches the cell boundary; the pin must be surrounded by non-fuel")]
    FuelOnBoundary(usize),
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
    let ng = totals.first().map_or(0, Vec::len);
    for (region, t) in totals.iter().enumerate() {
        let valid = |v: &f64| *v > 0.0 && v.is_finite();
        if !fuel[region] && (t.len() != ng || !t.iter().all(valid)) {
            return Err(DancoffError::Total(region));
        }
    }
    let cells = lattice.cells_per_side().pow(2);
    let mut cell_of = vec![0; regions];
    for cell in 0..cells {
        for region in lattice.cell_regions(cell) {
            cell_of[region] = cell;
        }
    }

    // 黒体の燃料（減衰 1 − exp(−τ) が 1 になる）と、燃料の外の物質
    let materials: Vec<Material> = (0..regions)
        .map(|region| Material {
            total: if fuel[region] {
                vec![f64::INFINITY; ng]
            } else {
                totals[region].clone()
            },
            scatter: vec![vec![0.0; ng]; ng],
            fission: vec![0.0; ng],
            nu_fission: vec![0.0; ng],
            fission_terms: Vec::new(),
        })
        .collect();
    let references: Vec<&Material> = materials.iter().collect();
    let mesh = CoarseMesh::new(tracking);
    let sweeper = Sweeper::new(tracking, &references, ng, &mesh);

    // 孤立したピンでの吸収率 R_iso / ψ0: 燃料に入るたびに方向の重みを足す
    let mut isolated = vec![0.0; cells];
    for track in &tracking.tracks {
        let weight: f64 = (0..POLAR_SIN.len())
            .map(|p| sweeper.direction_weight(track.azimuth, p))
            .sum();
        let n = track.segments.len();
        for forward in [true, false] {
            let mut previous = None;
            for step in 0..n {
                let s = if forward { step } else { n - 1 - step };
                let region = track.segments[s].region;
                let current = fuel[region].then_some(cell_of[region]);
                if let Some(cell) = current {
                    // 燃料が格子の端やセルの境目に接していると、孤立したピンへの入射として数えられない
                    if step == 0 || step == n - 1 || previous.is_some_and(|c| c != cell) {
                        return Err(DancoffError::FuelOnBoundary(cell));
                    }
                    if previous.is_none() {
                        isolated[cell] += weight;
                    }
                }
                previous = current;
            }
        }
    }

    let psi0 = 1.0 / (4.0 * PI);
    let equilibrium: Vec<Vec<f64>> = (0..regions)
        .map(|region| vec![if fuel[region] { 0.0 } else { psi0 }; ng])
        .collect();
    // 最初は境界から ψ0 が入ってくるとする（燃料は境界に接しないので、端の領域はどれも燃料の外）
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
            let sums = sum_tallies(&tallies, ng);
            let mut absorbed = vec![vec![0.0; ng]; cells];
            for region in (0..regions).filter(|&r| fuel[r]) {
                for (a, d) in absorbed[cell_of[region]]
                    .iter_mut()
                    .zip(&sums.delta[region * ng..][..ng])
                {
                    *a += d;
                }
            }
            return Ok(absorbed
                .into_iter()
                .zip(&isolated)
                .map(|(rates, &iso)| {
                    (iso > 0.0).then(|| rates.iter().map(|r| 1.0 - r / (psi0 * iso)).collect())
                })
                .collect());
        }
    }
    Err(DancoffError::NotConverged(options.max_iterations))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};

    fn options() -> DancoffOptions {
        DancoffOptions::default()
    }

    /// ピンと水のセルを `layout`（0: ピン、1: 水）に並べ、燃料以外の群の全断面積を `total` とした Dancoff 係数
    fn factors(pitch: f64, layout: Vec<Vec<usize>>, total: &[f64]) -> Vec<Option<Vec<f64>>> {
        let pin = PinCell::new(pitch, vec![0.4096]).unwrap();
        let water = PinCell::new(pitch, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, water], layout).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let regions = lattice.volumes().len();
        // ピンのセルは (燃料, 水) の2領域、水のセルは1領域
        let mut fuel = Vec::with_capacity(regions);
        for cell in 0..lattice.cells_per_side().pow(2) {
            let range = lattice.cell_regions(cell);
            fuel.push(range.len() == 2);
            fuel.extend(std::iter::repeat_n(false, range.len() - 1));
        }
        let totals = vec![total.to_vec(); regions];
        dancoff(&tracking, &totals, &fuel, &options()).unwrap()
    }

    #[test]
    fn distant_pins_do_not_shadow_each_other() {
        // 隣のピンまで減速材の光学的厚さが 20 程度あるので、影はほぼない
        let c = factors(20.0, vec![vec![0]], &[1.0, 2.0]);
        for value in c[0].as_ref().unwrap() {
            assert!(value.abs() < 1e-6, "{value}");
        }
    }

    #[test]
    fn transparent_moderator_shadows_completely() {
        let c = factors(1.26, vec![vec![0]], &[1.0e-6]);
        let value = c[0].as_ref().unwrap()[0];
        assert!(value > 0.999, "{value}");
    }

    #[test]
    fn repeated_pins_match_the_single_cell_lattice() {
        let total = [0.7, 1.5];
        let single = factors(1.26, vec![vec![0]], &total);
        let four = factors(1.26, vec![vec![0, 0], vec![0, 0]], &total);
        let reference = single[0].as_ref().unwrap();
        assert!(reference[0] > reference[1], "{reference:?}");
        assert!(
            reference.iter().all(|&c| c > 0.0 && c < 1.0),
            "{reference:?}"
        );
        for cell in &four {
            for (c, r) in cell.as_ref().unwrap().iter().zip(reference) {
                assert!((c - r).abs() < 2e-3, "{c} vs {r}");
            }
        }
    }

    #[test]
    fn cells_without_fuel_have_no_factors() {
        let c = factors(1.26, vec![vec![0, 1], vec![1, 0]], &[1.0]);
        assert!(c[0].is_some() && c[3].is_some());
        assert!(c[1].is_none() && c[2].is_none());
    }
}
