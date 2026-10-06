//! 平坦線源の特性曲線法（MOC）による固有値計算。

use std::f64::consts::PI;

use rayon::prelude::*;

use crate::cmfd::CoarseMesh;
use crate::material::{Material, MaterialError};
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

/// MOC の固有値計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum MocError {
    #[error("no regions")]
    NoRegions,
    #[error("all materials must have the same number of groups")]
    MixedGroups,
    #[error("region {region}")]
    Material {
        region: usize,
        source: MaterialError,
    },
    #[error("total cross sections must be positive (void is not supported)")]
    NonPositiveTotal,
    #[error("{materials} materials given, but the tracks cross {regions} regions")]
    RegionCount { materials: usize, regions: usize },
    #[error("no track crosses region {0}; use a finer track spacing")]
    Untracked(usize),
    #[error("no fissile material")]
    NoFission,
    #[error("not converged in {iterations} iterations (k = {k})")]
    NotConverged { iterations: usize, k: f64 },
}

/// `materials[i]` は領域 i の物質。
pub fn solve(
    tracking: &Tracking,
    materials: &[&Material],
    options: &Options,
) -> Result<Solution, MocError> {
    let regions = materials.len();
    let ng = materials.first().ok_or(MocError::NoRegions)?.groups();
    if materials.iter().any(|m| m.groups() != ng) {
        return Err(MocError::MixedGroups);
    }
    for (region, m) in materials.iter().enumerate() {
        m.validate()
            .map_err(|source| MocError::Material { region, source })?;
    }
    if materials
        .iter()
        .any(|m| m.total.iter().any(|&t| t <= 0.0 || !t.is_finite()))
    {
        return Err(MocError::NonPositiveTotal);
    }
    let used = tracking
        .tracks
        .iter()
        .flat_map(|t| &t.segments)
        .map(|s| s.region)
        .max()
        .map_or(0, |r| r + 1);
    if used != regions {
        return Err(MocError::RegionCount {
            materials: regions,
            regions: used,
        });
    }
    let volumes = tracking.volumes(regions);
    if let Some(i) = volumes.iter().position(|&v| v <= 0.0) {
        return Err(MocError::Untracked(i));
    }

    let mesh = CoarseMesh::new(tracking);
    let sweeper = Sweeper::new(tracking, materials, ng, &mesh);
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
        return Err(MocError::NoFission);
    }
    for f in flux.iter_mut().flatten() {
        *f /= initial;
    }
    let mut k = 1.0;
    // スイープの終わりに出ていく角度中性子束（`Sweeper::sweep` を参照）
    let mut outgoing = vec![0.0; sweeper.outgoing_len()];
    let mut next_outgoing = outgoing.clone();
    let mut tallies = sweeper.tallies(regions);
    let mut k_history = Remaining::default();
    let mut flux_history = Remaining::default();

    for iteration in 1..=options.max_iterations {
        // 等方線源を Σt で割ったもの（角度中性子束の平衡値）
        let equilibrium: Vec<Vec<f64>> = (0..regions)
            .map(|i| equilibrium_flux(materials[i], &flux[i], k))
            .collect();

        sweeper.sweep(&equilibrium, &outgoing, &mut next_outgoing, &mut tallies);
        std::mem::swap(&mut outgoing, &mut next_outgoing);
        let Sums {
            delta,
            mut leak,
            current: _,
        } = sum_tallies(&tallies, ng);

        let mut new_flux: Vec<Vec<f64>> = (0..regions)
            .map(|i| {
                let total = &materials[i].total;
                let delta = &delta[i * ng..][..ng];
                (0..ng)
                    .map(|g| 4.0 * PI * equilibrium[i][g] + delta[g] / (total[g] * volumes[i]))
                    .collect()
            })
            .collect();
        // べき乗法: 生成率の比で k を更新し、中性子束と境界の値を生成率 1 に規格化する
        let rate = production(&new_flux);
        let new_k = k * rate;
        for f in new_flux.iter_mut().flatten() {
            *f /= rate;
        }
        for l in &mut leak {
            *l /= rate;
        }
        for b in &mut outgoing {
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
                leakage: leak,
                iterations: iteration,
            });
        }
    }
    Err(MocError::NotConverged {
        iterations: options.max_iterations,
        k,
    })
}

/// トラックのスイープ。トラックと物質から決まる量（減衰の表など）を、反復の前に一度だけ求めて持つ。
struct Sweeper<'a> {
    tracking: &'a Tracking,
    /// 正味の中性子流を数える粗メッシュ
    mesh: &'a CoarseMesh,
    ng: usize,
    /// [極角][群] の数
    width: usize,
    /// [トラック] の区間ごと・極角ごと・群ごとの 1 - exp(-τ)
    attenuation: Vec<Vec<f64>>,
    /// [トラック][向き] のスイープの始めに入ってくる角度中性子束の出どころ（`incoming_sources`）
    sources: Vec<Option<usize>>,
    /// 方位角ごとのトラックの番号の範囲
    azimuths: Vec<std::ops::Range<usize>>,
}

impl<'a> Sweeper<'a> {
    fn new(
        tracking: &'a Tracking,
        materials: &[&Material],
        ng: usize,
        mesh: &'a CoarseMesh,
    ) -> Self {
        let width = POLAR_SIN.len() * ng;
        let attenuation = tracking
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
        Self {
            tracking,
            mesh,
            ng,
            width,
            attenuation,
            sources: incoming_sources(tracking),
            azimuths: azimuth_ranges(tracking),
        }
    }

    /// 出ていく角度中性子束の配列の長さ
    fn outgoing_len(&self) -> usize {
        2 * self.tracking.tracks.len() * self.width
    }

    /// 方位角ごとの部分和の入れ物。方位角ごとに順にスイープし、方位角の順に足すので、結果は
    /// スレッドの数や仕事の分け方によらない。方位角ごとのトラックの長さの和はほぼ等しく、
    /// 仕事の量も揃う
    fn tallies(&self, regions: usize) -> Vec<Tally> {
        self.azimuths
            .iter()
            .map(|_| Tally::new(regions, self.mesh.faces(), self.ng, self.width))
            .collect()
    }

    /// 方向ごとの重み 4π (w_a/2) ω_p sinθ_p δ_a
    fn direction_weight(&self, azimuth: usize, p: usize) -> f64 {
        let a = &self.tracking.azimuths[azimuth];
        4.0 * PI * 0.5 * a.weight * POLAR_WEIGHT[p] * POLAR_SIN[p] * a.spacing
    }

    /// 全トラックを両方向に1回ずつスイープする。
    ///
    /// `equilibrium` は [領域][群] の角度中性子束の平衡値（線源 / 4πΣt）。`incoming` は前回の
    /// スイープで出ていった角度中性子束、`outgoing` には今回出ていく分を書く。どちらも
    /// [トラック][向き (0: 順, 1: 逆)][極角][群] の並びで、各トラックは自分の場所にだけ書くので、
    /// トラックごとに並列にスイープできる。
    fn sweep(
        &self,
        equilibrium: &[Vec<f64>],
        incoming: &[f64],
        outgoing: &mut [f64],
        tallies: &mut [Tally],
    ) {
        let (ng, width) = (self.ng, self.width);
        let slot = |track: usize, forward: bool| 2 * track + usize::from(!forward);
        // 出ていく角度中性子束の配列を、方位角ごとの範囲に分ける
        let mut rest = outgoing;
        let mut outputs = Vec::with_capacity(self.azimuths.len());
        for range in &self.azimuths {
            let (head, tail) = rest.split_at_mut(2 * width * range.len());
            outputs.push(head);
            rest = tail;
        }
        tallies
            .par_iter_mut()
            .zip(outputs)
            .zip(&self.azimuths)
            .for_each(|((tally, outputs), range)| {
                tally.clear();
                for (index, out) in range.clone().zip(outputs.chunks_mut(2 * width)) {
                    let track = &self.tracking.tracks[index];
                    let weights: [f64; POLAR_SIN.len()] =
                        std::array::from_fn(|p| self.direction_weight(track.azimuth, p));
                    for forward in [true, false] {
                        let psi = &mut tally.psi;
                        match self.sources[slot(index, forward)] {
                            Some(source) => {
                                psi.copy_from_slice(&incoming[source * width..][..width])
                            }
                            // 真空境界から入ってくる角度中性子束は 0
                            None => psi.fill(0.0),
                        }
                        // 外周の面は順方向の向きで持つので、逆方向では向きを反対にする
                        let [start, end] = self.mesh.ends(index);
                        let (enter, exit, direction) = if forward {
                            (start, end, 1.0)
                        } else {
                            (end, start, -1.0)
                        };
                        add_current(
                            &mut tally.current,
                            &[(enter.0, direction * enter.1)],
                            psi,
                            &weights,
                        );
                        let n = track.segments.len();
                        let first = if forward { 0 } else { n - 1 };
                        let mut cell = self.mesh.cell(track.segments[first].region);
                        for step in 0..n {
                            let s = if forward { step } else { n - 1 - step };
                            let region = track.segments[s].region;
                            // 隣のセルに移るときに通る面
                            let next = self.mesh.cell(region);
                            if next != cell {
                                let path = self.mesh.path(cell, next);
                                add_current(&mut tally.current, path.faces(), psi, &weights);
                                cell = next;
                            }
                            let att = &self.attenuation[index][s * width..][..width];
                            let q = &equilibrium[region];
                            let d = &mut tally.delta[region * ng..][..ng];
                            for (p, weight) in weights.iter().enumerate() {
                                for g in 0..ng {
                                    let j = p * ng + g;
                                    let change = (psi[j] - q[g]) * att[j];
                                    psi[j] -= change;
                                    d[g] += weight * change;
                                }
                            }
                        }
                        add_current(
                            &mut tally.current,
                            &[(exit.0, direction * exit.1)],
                            psi,
                            &weights,
                        );
                        let link = if forward {
                            track.forward_link
                        } else {
                            track.backward_link
                        };
                        // 真空境界から出た分は漏れ
                        if link.is_none() {
                            for (p, weight) in weights.iter().enumerate() {
                                for (l, value) in tally.leak.iter_mut().zip(&psi[p * ng..][..ng]) {
                                    *l += weight * value;
                                }
                            }
                        }
                        out[usize::from(!forward) * width..][..width].copy_from_slice(psi);
                    }
                }
            });
    }
}

/// スイープで足し込む量（方位角ごとに持ち、スイープの後に足し合わせる）。
struct Tally {
    /// [領域][群] の角度中性子束の減衰の重み付き和
    delta: Vec<f64>,
    /// [群] の真空境界からの漏れ
    leak: Vec<f64>,
    /// [粗メッシュの面][群] の正味の中性子流（面の正の向きが正）
    current: Vec<f64>,
    /// スイープ中の角度中性子束 [極角][群]（作業用）
    psi: Vec<f64>,
}

impl Tally {
    fn new(regions: usize, faces: usize, groups: usize, width: usize) -> Self {
        Self {
            delta: vec![0.0; regions * groups],
            leak: vec![0.0; groups],
            current: vec![0.0; faces * groups],
            psi: vec![0.0; width],
        }
    }

    fn clear(&mut self) {
        self.delta.fill(0.0);
        self.leak.fill(0.0);
        self.current.fill(0.0);
    }
}

/// 方位角ごとの部分和を足し合わせたもの。
struct Sums {
    /// [領域][群] の角度中性子束の減衰の重み付き和
    delta: Vec<f64>,
    /// [群] の真空境界からの漏れ
    leak: Vec<f64>,
    /// [粗メッシュの面][群] の正味の中性子流
    current: Vec<f64>,
}

/// 方位角ごとの部分和を方位角の順に足す。
fn sum_tallies(tallies: &[Tally], ng: usize) -> Sums {
    // 要素ごとに並列に足す（1つの要素の中では方位角の順なので、結果は分け方によらない）
    let sum = |part: fn(&Tally) -> &[f64]| -> Vec<f64> {
        let mut sum = vec![0.0; part(&tallies[0]).len()];
        sum.par_chunks_mut(ng).enumerate().for_each(|(i, chunk)| {
            for tally in tallies {
                for (a, b) in chunk.iter_mut().zip(&part(tally)[i * ng..][..ng]) {
                    *a += b;
                }
            }
        });
        sum
    };
    Sums {
        delta: sum(|t| &t.delta),
        leak: sum(|t| &t.leak),
        current: sum(|t| &t.current),
    }
}

/// 面 `faces`（面の正の向きに通るなら重み正）を通る角度中性子束 `psi` [極角][群] の流れを、
/// `current` [面][群] に足す。
fn add_current(
    current: &mut [f64],
    faces: &[(usize, f64)],
    psi: &[f64],
    weights: &[f64; POLAR_SIN.len()],
) {
    let ng = psi.len() / weights.len();
    for &(face, sign) in faces {
        let c = &mut current[face * ng..][..ng];
        for (p, weight) in weights.iter().enumerate() {
            for (c, value) in c.iter_mut().zip(&psi[p * ng..][..ng]) {
                *c += sign * weight * value;
            }
        }
    }
}

/// 方位角ごとのトラックの番号の範囲（トラックは方位角の順に並ぶ）。
fn azimuth_ranges(tracking: &Tracking) -> Vec<std::ops::Range<usize>> {
    let mut ranges = Vec::with_capacity(tracking.azimuths.len());
    let mut start = 0;
    for chunk in tracking.tracks.chunk_by(|a, b| a.azimuth == b.azimuth) {
        assert_eq!(
            chunk[0].azimuth,
            ranges.len(),
            "tracks are not ordered by azimuth"
        );
        ranges.push(start..start + chunk.len());
        start += chunk.len();
    }
    ranges
}

/// [トラック][向き] のスイープの始めに入ってくる角度中性子束が、どのトラックのどの向きのスイープの
/// 終わりから来るか（同じ並びの番号。真空境界なら `None`）。
fn incoming_sources(tracking: &Tracking) -> Vec<Option<usize>> {
    let slot = |track: usize, forward: bool| 2 * track + usize::from(!forward);
    let mut sources = vec![None; 2 * tracking.tracks.len()];
    for (index, track) in tracking.tracks.iter().enumerate() {
        for (forward, link) in [(true, track.forward_link), (false, track.backward_link)] {
            if let Some(link) = link {
                let target = &mut sources[slot(link.track, link.forward)];
                assert!(
                    target.is_none(),
                    "two tracks reflect into track {} ({})",
                    link.track,
                    if link.forward { "forward" } else { "backward" }
                );
                *target = Some(slot(index, forward));
            }
        }
    }
    sources
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
        // 反復は相対変化で打ち切るので、収支の誤差は全体の量（生成/k）に比例する
        assert!(
            (absorption + leakage - production / s.k).abs() < 1e-6 * production / s.k,
            "{absorption} + {leakage} vs {}",
            production / s.k
        );
        assert!(s.k < solve_with(Boundaries::REFLECTIVE).k);
    }

    #[test]
    fn result_does_not_depend_on_the_thread_count() {
        // スイープの和を足す順がスレッドの分け方で変わると、結果が最下位のビットで揺れる
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let (f, m) = (fuel(), moderator());
        let materials = [&f, &m, &m, &m, &f, &m];
        let solve_on = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| solve(&tracking, &materials, &Options::default()).unwrap())
        };
        let serial = solve_on(1);
        for threads in [2, 5, 8] {
            assert_eq!(solve_on(threads), serial, "{threads} threads");
        }
    }

    #[test]
    fn sweep_currents_balance_each_coarse_cell() {
        // 1回のスイープで、セルの中の領域の減衰の和 Σ Δ（= 入る - 出る）と、セルの面から外向きの
        // 正味の中性子流の和は打ち消し合う。線源と入ってくる角度中性子束は任意の値でよい
        let pin = PinCell::new(1.26, vec![0.3, 0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
        let mesh = CoarseMesh::new(&tracking);
        // 2×2 の格子の中心を通るトラックがあり、セルの角を斜めに抜ける場合も確かめる
        let corner = |track: &crate::tracking::Track| {
            track.segments.windows(2).any(|w| {
                let (a, b) = (mesh.cell(w[0].region), mesh.cell(w[1].region));
                a / 2 != b / 2 && a % 2 != b % 2
            })
        };
        assert!(tracking.tracks.iter().any(corner));
        let (f, m) = (fuel(), moderator());
        let materials = [&f, &f, &m, &m, &m, &f, &f, &m];
        let ng = 2;
        let sweeper = Sweeper::new(&tracking, &materials, ng, &mesh);
        let equilibrium: Vec<Vec<f64>> = (0..materials.len())
            .map(|i| {
                (0..ng)
                    .map(|g| 0.1 + ((7 * i + 3 * g) % 11) as f64 * 0.01)
                    .collect()
            })
            .collect();
        let incoming: Vec<f64> = (0..sweeper.outgoing_len())
            .map(|j| 0.05 + ((37 * j) % 101) as f64 * 1.0e-3)
            .collect();
        let mut outgoing = vec![0.0; sweeper.outgoing_len()];
        let mut tallies = sweeper.tallies(materials.len());
        sweeper.sweep(&equilibrium, &incoming, &mut outgoing, &mut tallies);
        let sums = sum_tallies(&tallies, ng);
        for cell in 0..mesh.cells() {
            for g in 0..ng {
                let inside: f64 = lattice
                    .cell_regions(cell)
                    .map(|i| sums.delta[i * ng + g])
                    .sum();
                let outward: f64 = mesh
                    .cell_faces(cell)
                    .iter()
                    .map(|&(face, sign)| sign * sums.current[face * ng + g])
                    .sum();
                assert!(
                    (inside + outward).abs() < 1e-12 * inside.abs().max(outward.abs()),
                    "cell {cell} group {g}: {inside} + {outward}"
                );
            }
        }
    }

    /// CMFD あり・なしで解き、同じ解に、より少ない反復で収束することを確かめる。
    fn assert_cmfd_accelerates(tracking: &Tracking, materials: &[&Material]) {
        let plain = Options {
            acceleration: Acceleration::None,
            ..Options::default()
        };
        let slow = solve(tracking, materials, &plain).unwrap();
        let fast = solve(tracking, materials, &Options::default()).unwrap();
        assert!((fast.k - slow.k).abs() < 1e-6, "{} vs {}", fast.k, slow.k);
        for (a, b) in fast.flux.iter().flatten().zip(slow.flux.iter().flatten()) {
            assert!(relative_change(*a, *b) < 1e-5, "{a} vs {b}");
        }
        assert!(
            5 * fast.iterations <= slow.iterations,
            "{} vs {} iterations",
            fast.iterations,
            slow.iterations
        );
    }

    #[test]
    fn cmfd_accelerates_a_pin_cell() {
        let cell = PinCell::new(1.26, vec![0.25, 0.4096, 0.55]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.05).unwrap();
        let (f, m) = (fuel(), moderator());
        assert_cmfd_accelerates(&tracking, &[&f, &f, &m, &m]);
    }

    #[test]
    fn cmfd_accelerates_a_lattice_with_vacuum_sides() {
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let layout = vec![vec![0, 1, 1], vec![1, 0, 1], vec![1, 1, 1]];
        let lattice = Lattice::new(vec![pin, water], layout).unwrap();
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
        let (f, m) = (fuel(), moderator());
        // セルの順: ピン, 水, 水 / 水, ピン, 水 / 水, 水, 水（ピンは燃料と水の2領域）
        let materials = [&f, &m, &m, &m, &m, &f, &m, &m, &m, &m, &m];
        assert_cmfd_accelerates(&tracking, &materials);
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
        let err = solve(&tracking, &[&m, &m, &m], &Options::default())
            .unwrap_err()
            .to_string();
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
    fn material_errors_keep_the_cause_in_the_chain_only() {
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
        let mut m = fuel();
        m.nu_fission[1] *= 2.0;
        let err = solve(&tracking, &[&m, &m], &Options::default()).unwrap_err();
        let cause = std::error::Error::source(&err).expect("the material error is the cause");
        assert!(!err.to_string().contains(&cause.to_string()), "{err}");
    }

    #[test]
    fn inconsistent_material_is_error() {
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
        let mut m = fuel();
        m.fission_terms.clear();
        let err = solve(&tracking, &[&m, &m], &Options::default()).unwrap_err();
        assert!(
            matches!(
                err,
                MocError::Material {
                    region: 0,
                    source: MaterialError::NuFission
                }
            ),
            "{err}"
        );
    }

    #[test]
    fn region_count_must_match() {
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
        let m = fuel();
        assert!(solve(&tracking, &[&m], &Options::default()).is_err());
    }
}
