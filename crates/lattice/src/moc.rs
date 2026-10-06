//! 平坦線源の特性曲線法（MOC）による固有値計算。

use std::f64::consts::PI;

use crate::material::Material;
use crate::tracking::Tracking;

/// 極角求積（Tabuchi-Yamamoto 3点、0 < θ < π/2 の半球で重みの和が 1）
const POLAR_SIN: [f64; 3] = [0.166_648, 0.537_707, 0.932_954];
const POLAR_WEIGHT: [f64; 3] = [0.046_233, 0.283_619, 0.670_148];

/// 反復の打ち切り条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Options {
    /// 実効増倍率の変化の上限
    pub k_tolerance: f64,
    /// 中性子束の相対変化の上限
    pub flux_tolerance: f64,
    pub max_iterations: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            k_tolerance: 1.0e-7,
            flux_tolerance: 1.0e-6,
            max_iterations: 20_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    pub k: f64,
    /// [領域][群] の領域平均スカラー中性子束（全体の生成率が 1 になるよう規格化）
    pub flux: Vec<Vec<f64>>,
    /// トラックから推定した領域の体積
    pub volumes: Vec<f64>,
    /// [群] の真空境界からの漏れ（中性子束と同じ規格化）
    pub leakage: Vec<f64>,
    pub iterations: usize,
}

/// `materials[i]` は領域 i の物質。
pub fn solve(
    tracking: &Tracking,
    materials: &[&Material],
    options: &Options,
) -> Result<Solution, String> {
    let regions = materials.len();
    let ng = materials.first().ok_or("no regions")?.groups();
    if materials.iter().any(|m| m.groups() != ng) {
        return Err("all materials must have the same number of groups".into());
    }
    for (i, m) in materials.iter().enumerate() {
        m.validate().map_err(|e| format!("region {i}: {e}"))?;
    }
    if materials
        .iter()
        .any(|m| m.total.iter().any(|&t| t <= 0.0 || !t.is_finite()))
    {
        return Err("total cross sections must be positive (void is not supported)".into());
    }
    let used = tracking
        .tracks
        .iter()
        .flat_map(|t| &t.segments)
        .map(|s| s.region)
        .max()
        .map_or(0, |r| r + 1);
    if used != regions {
        return Err(format!(
            "{regions} materials given, but the tracks cross {used} regions"
        ));
    }
    let volumes = tracking.volumes(regions);
    if let Some(i) = volumes.iter().position(|&v| v <= 0.0) {
        return Err(format!(
            "no track crosses region {i}; use a finer track spacing"
        ));
    }

    let np = POLAR_SIN.len();
    let width = np * ng;
    // 区間ごと・極角ごと・群ごとの 1 - exp(-τ)（スイープの前に一度だけ求める）
    let attenuation: Vec<Vec<f64>> = tracking
        .tracks
        .iter()
        .map(|track| {
            let mut values = Vec::with_capacity(track.segments.len() * width);
            for s in &track.segments {
                for sin in POLAR_SIN {
                    for t in &materials[s.region].total {
                        values.push(-(-t * s.length / sin).exp_m1());
                    }
                }
            }
            values
        })
        .collect();
    // 方向ごとの重み 4π (w_a/2) ω_p sinθ_p δ_a
    let direction_weight = |azimuth: usize, p: usize| {
        let a = &tracking.azimuths[azimuth];
        4.0 * PI * 0.5 * a.weight * POLAR_WEIGHT[p] * POLAR_SIN[p] * a.spacing
    };
    let production = |flux: &[Vec<f64>]| -> f64 {
        flux.iter()
            .zip(materials)
            .zip(&volumes)
            .map(|((phi, m), v)| {
                v * phi
                    .iter()
                    .zip(&m.nu_fission)
                    .map(|(f, n)| f * n)
                    .sum::<f64>()
            })
            .sum()
    };

    let mut flux = vec![vec![1.0; ng]; regions];
    let initial = production(&flux);
    if initial <= 0.0 {
        return Err("no fissile material".into());
    }
    for f in flux.iter_mut().flatten() {
        *f /= initial;
    }
    let mut k = 1.0;
    // 境界の角度中性子束 [トラック][向き (0: 順, 1: 逆)][極角][群]
    let slot = |track: usize, forward: bool| (2 * track + usize::from(!forward)) * width;
    let mut boundary = vec![0.0; 2 * tracking.tracks.len() * width];
    let mut next_boundary = boundary.clone();
    let mut psi = vec![0.0; width];
    let mut k_history = Remaining::default();
    let mut flux_history = Remaining::default();

    for iteration in 1..=options.max_iterations {
        // 等方線源を Σt で割ったもの（角度中性子束の平衡値）
        let equilibrium: Vec<Vec<f64>> = (0..regions)
            .map(|i| equilibrium_flux(materials[i], &flux[i], k))
            .collect();

        let mut delta = vec![vec![0.0; ng]; regions];
        for (index, track) in tracking.tracks.iter().enumerate() {
            let weights: Vec<f64> = (0..np)
                .map(|p| direction_weight(track.azimuth, p))
                .collect();
            for forward in [true, false] {
                psi.copy_from_slice(&boundary[slot(index, forward)..][..width]);
                let n = track.segments.len();
                for step in 0..n {
                    let s = if forward { step } else { n - 1 - step };
                    let region = track.segments[s].region;
                    let att = &attenuation[index][s * width..][..width];
                    let q = &equilibrium[region];
                    let d = &mut delta[region];
                    for (p, weight) in weights.iter().enumerate() {
                        for g in 0..ng {
                            let j = p * ng + g;
                            let change = (psi[j] - q[g]) * att[j];
                            psi[j] -= change;
                            d[g] += weight * change;
                        }
                    }
                }
                let link = if forward {
                    track.forward_link
                } else {
                    track.backward_link
                };
                match link {
                    Some(link) => next_boundary[slot(link.track, link.forward)..][..width]
                        .copy_from_slice(&psi),
                    None => todo!("leakage"),
                }
            }
        }
        std::mem::swap(&mut boundary, &mut next_boundary);

        let mut new_flux: Vec<Vec<f64>> = (0..regions)
            .map(|i| {
                let total = &materials[i].total;
                (0..ng)
                    .map(|g| 4.0 * PI * equilibrium[i][g] + delta[i][g] / (total[g] * volumes[i]))
                    .collect()
            })
            .collect();
        // べき乗法: 生成率の比で k を更新し、中性子束と境界の値を生成率 1 に規格化する
        let rate = production(&new_flux);
        let new_k = k * rate;
        for f in new_flux.iter_mut().flatten() {
            *f /= rate;
        }
        for b in &mut boundary {
            *b /= rate;
        }
        let change = new_flux
            .iter()
            .flatten()
            .zip(flux.iter().flatten())
            .map(|(new, old)| relative_change(*new, *old))
            .fold(0.0, f64::max);
        // 収束が遅いと1回の変化は小さくても真の解から遠いので、変化の縮小率から
        // 残りの誤差を見積もって判定する
        // 両方の履歴を毎回更新する（短絡評価で片方が飛ぶと、縮小率を離れた反復どうしで比べてしまう）
        let k_remaining = k_history.update((new_k - k).abs());
        let flux_remaining = flux_history.update(change);
        let converged =
            k_remaining < options.k_tolerance && flux_remaining < options.flux_tolerance;
        flux = new_flux;
        k = new_k;
        if converged {
            return Ok(Solution {
                k,
                flux,
                volumes,
                leakage: vec![],
                iterations: iteration,
            });
        }
    }
    Err(format!(
        "not converged in {} iterations (k = {k})",
        options.max_iterations
    ))
}

/// 中性子束の相対変化。両方 0 なら変化なし。
fn relative_change(new: f64, old: f64) -> f64 {
    let scale = new.abs().max(old.abs());
    if scale == 0.0 {
        0.0
    } else {
        (new - old).abs() / scale
    }
}

/// 反復の変化量の列から、収束までに残る誤差を見積もる。
///
/// 変化が縮小率 ρ で等比的に減るとき、残りの誤差は 最新の変化 × ρ/(1-ρ)。
#[derive(Debug, Default)]
struct Remaining {
    previous: Option<f64>,
}

impl Remaining {
    /// 最新の変化量を受け取り、残りの誤差の見積もり（変化量そのもの以上）を返す。
    fn update(&mut self, change: f64) -> f64 {
        let estimate = match self.previous {
            // 浮動小数点の不動点に達した
            _ if change == 0.0 => 0.0,
            Some(previous) if previous > 0.0 => {
                let ratio = change / previous;
                if ratio < 1.0 {
                    change * (1.0 + ratio / (1.0 - ratio))
                } else {
                    f64::INFINITY
                }
            }
            _ => f64::INFINITY,
        };
        self.previous = Some(change);
        estimate
    }
}

/// 散乱源と核分裂源（/k）の和を 4πΣt で割った、角度中性子束の平衡値。
fn equilibrium_flux(m: &Material, phi: &[f64], k: f64) -> Vec<f64> {
    let ng = phi.len();
    let mut source = vec![0.0; ng];
    for (from, row) in m.scatter.iter().enumerate() {
        for (to, s) in row.iter().enumerate() {
            source[to] += s * phi[from];
        }
    }
    for term in &m.fission_terms {
        let rate: f64 = term.nu_fission.iter().zip(phi).map(|(n, f)| n * f).sum();
        for (s, chi) in source.iter_mut().zip(&term.chi) {
            *s += chi * rate / k;
        }
    }
    source
        .iter()
        .zip(&m.total)
        .map(|(q, t)| q / (4.0 * PI * t))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};
    use crate::material::FissionTerm;
    use crate::tracking::{Boundaries, Boundary};

    fn two_group(total: [f64; 2], scatter: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Material {
        Material {
            total: total.to_vec(),
            scatter: scatter.iter().map(|r| r.to_vec()).collect(),
            fission: nu_fission.iter().map(|v| v / 2.5).collect(),
            nu_fission: nu_fission.to_vec(),
            fission_terms: if nu_fission.iter().any(|&v| v > 0.0) {
                vec![FissionTerm {
                    chi: vec![1.0, 0.0],
                    nu_fission: nu_fission.to_vec(),
                }]
            } else {
                vec![]
            },
        }
    }

    fn fuel() -> Material {
        two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25])
    }

    fn moderator() -> Material {
        two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0])
    }

    /// 無限媒質の k = νΣfᵀ (diag Σt - Σsᵀ)⁻¹ χ（χ = (1, 0)）
    fn infinite_medium_k(m: &Material) -> f64 {
        let a = [
            [m.total[0] - m.scatter[0][0], -m.scatter[1][0]],
            [-m.scatter[0][1], m.total[1] - m.scatter[1][1]],
        ];
        let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        let phi = [a[1][1] / det, -a[1][0] / det];
        m.nu_fission[0] * phi[0] + m.nu_fission[1] * phi[1]
    }

    #[test]
    fn homogeneous_cell_matches_infinite_medium() {
        let cell = PinCell::new(1.26, vec![0.3, 0.45]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
        let m = fuel();
        let solution = solve(&tracking, &[&m, &m, &m], &Options::default()).unwrap();
        let expected = infinite_medium_k(&m);
        assert!(
            (solution.k - expected).abs() < 1e-6,
            "{} vs {expected}",
            solution.k
        );
    }

    fn heterogeneous(n_azimuthal: usize, spacing: f64) -> Solution {
        let cell = PinCell::new(1.26, vec![0.25, 0.4096, 0.55]).unwrap();
        let tracking =
            Tracking::new(&Lattice::pin_cell(cell.clone()), n_azimuthal, spacing).unwrap();
        let (f, m) = (fuel(), moderator());
        solve(&tracking, &[&f, &f, &m, &m], &Options::default()).unwrap()
    }

    #[test]
    fn heterogeneous_cell_balances_neutrons() {
        let s = heterogeneous(16, 0.05);
        let materials = [fuel(), fuel(), moderator(), moderator()];
        let (mut absorption, mut production) = (0.0, 0.0);
        for (i, m) in materials.iter().enumerate() {
            for g in 0..2 {
                absorption += s.volumes[i] * m.absorption()[g] * s.flux[i][g];
                production += s.volumes[i] * m.nu_fission[g] * s.flux[i][g];
            }
        }
        assert!((absorption - production / s.k).abs() < 1e-6 * absorption);
        // 熱中性子束は減速材で高く、燃料の中心で低い
        assert!(s.flux[3][1] > s.flux[0][1]);
    }

    #[test]
    fn reflective_cell_has_no_leakage() {
        assert_eq!(heterogeneous(8, 0.1).leakage, vec![0.0, 0.0]);
    }

    #[test]
    fn vacuum_lattice_balances_production_with_absorption_and_leakage() {
        // 2×2 の格子（燃料のピン2本と水のセル2つ）の外周を真空にすると、漏れの分だけ k が下がる
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
        let (f, m) = (fuel(), moderator());
        let materials = [&f, &m, &m, &m, &f, &m];
        let solve_with = |boundaries| {
            let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
            solve(&tracking, &materials, &Options::default()).unwrap()
        };
        let all_vacuum = Boundaries {
            left: Boundary::Vacuum,
            right: Boundary::Vacuum,
            top: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
        };
        let s = solve_with(all_vacuum);
        let (mut absorption, mut production) = (0.0, 0.0);
        for (i, m) in materials.iter().enumerate() {
            for g in 0..2 {
                absorption += s.volumes[i] * m.absorption()[g] * s.flux[i][g];
                production += s.volumes[i] * m.nu_fission[g] * s.flux[i][g];
            }
        }
        let leakage: f64 = s.leakage.iter().sum();
        assert!(leakage > 0.0);
        assert!(
            (absorption + leakage - production / s.k).abs() < 1e-6 * absorption,
            "{absorption} + {leakage} vs {}",
            production / s.k
        );
        assert!(s.k < solve_with(Boundaries::REFLECTIVE).k);
    }

    #[test]
    fn refinement_converges() {
        let coarse = heterogeneous(16, 0.05);
        let fine = heterogeneous(32, 0.02);
        assert!(
            (coarse.k - fine.k).abs() < 1e-3,
            "{} vs {}",
            coarse.k,
            fine.k
        );
    }

    #[test]
    fn no_fission_is_error() {
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
        let m = moderator();
        assert!(solve(&tracking, &[&m, &m], &Options::default()).is_err());
    }

    #[test]
    fn region_without_tracks_is_error() {
        // 中心の極小の円には、中心からずれたトラックが1本も通らない
        let cell = PinCell::new(1.26, vec![1.0e-6, 0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 4, 0.3).unwrap();
        let m = fuel();
        assert!(tracking.volumes(3)[0] == 0.0);
        let err = solve(&tracking, &[&m, &m, &m], &Options::default()).unwrap_err();
        assert!(err.contains("region 0"), "{err}");
    }

    #[test]
    fn relative_change_handles_zero_flux() {
        assert_eq!(relative_change(0.0, 0.0), 0.0);
        assert_eq!(relative_change(1.0, 0.0), 1.0);
        assert_eq!(relative_change(0.0, 2.0), 1.0);
        assert!((relative_change(1.1, 1.0) - 0.1 / 1.1).abs() < 1e-15);
    }

    #[test]
    fn remaining_error_estimate() {
        let mut r = Remaining::default();
        assert_eq!(r.update(1.0), f64::INFINITY);
        // 縮小率 0.5 なら、残りは最新の変化と同じだけ
        assert!((r.update(0.5) - 1.0).abs() < 1e-15);
        // 変化が 0 になれば収束
        let mut r = Remaining::default();
        r.update(1.0e-16);
        assert_eq!(r.update(0.0), 0.0);
        let mut r = Remaining::default();
        r.update(0.0);
        assert_eq!(r.update(0.0), 0.0);
    }

    #[test]
    fn inconsistent_material_is_error() {
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
        let mut m = fuel();
        m.fission_terms.clear();
        let err = solve(&tracking, &[&m, &m], &Options::default()).unwrap_err();
        assert!(err.contains("fission terms"), "{err}");
    }

    #[test]
    fn region_count_must_match() {
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
        let m = fuel();
        assert!(solve(&tracking, &[&m], &Options::default()).is_err());
    }
}
