//! 平坦線源の特性曲線法（MOC）による固有値計算。

use std::f64::consts::PI;
use std::sync::{Mutex, PoisonError};

use rayon::prelude::*;

pub mod dancoff;
pub mod fixed_source;
mod harmonics;

use crate::cmfd::{CmfdError, CmfdOptions, Coarse, CoarseMesh};
use crate::material::{Material, MaterialError};
use crate::tracking::Tracking;
use harmonics::{planar_moments, real_harmonic};

/// 極角求積（Tabuchi-Yamamoto 3点、0 < θ < π/2 の半球で重みの和が 1）
const POLAR_SIN: [f64; 3] = [0.166_648, 0.537_707, 0.932_954];
const POLAR_WEIGHT: [f64; 3] = [0.046_233, 0.283_619, 0.670_148];

/// 反復の打ち切り条件と、散乱の非等方性の次数。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Options {
    /// 散乱の非等方性を展開するルジャンドルの次数 L（0 なら等方散乱）。L ≥ 1 では、真空でない物質に
    /// P1〜PL の散乱行列が要る
    pub scattering_order: usize,
    /// 実効増倍率の変化の上限
    pub k_tolerance: f64,
    /// 中性子束の相対変化の上限
    pub flux_tolerance: f64,
    pub max_iterations: usize,
    pub acceleration: Acceleration,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            scattering_order: 0,
            k_tolerance: 1.0e-7,
            flux_tolerance: 1.0e-6,
            max_iterations: 20_000,
            acceleration: Acceleration::Cmfd(CmfdOptions::default()),
        }
    }
}

/// べき乗法の加速。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Acceleration {
    /// 加速しない（スイープごとに生成率の比で k を更新する）
    None,
    /// スイープごとにピンセル単位の CMFD を解き、k と中性子束を更新する
    Cmfd(CmfdOptions),
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
    #[error("total cross sections must be positive in every group, or zero in every group (void)")]
    NonPositiveTotal,
    #[error("{materials} materials given, but the tracks cross {regions} regions")]
    RegionCount { materials: usize, regions: usize },
    #[error("no track crosses region {0}; use a finer track spacing")]
    Untracked(usize),
    #[error("cell {0} is all void")]
    VoidCell(usize),
    #[error("no fissile material")]
    NoFission,
    #[error("no neutron reaches group {0} (no fission spectrum or scattering into it)")]
    GroupWithoutSource(usize),
    #[error(
        "region {region}: scattering order {requested} needs P1 to P{requested} matrices, \
         but the material has them up to P{available}"
    )]
    ScatteringOrder {
        region: usize,
        available: usize,
        requested: usize,
    },
    #[error(transparent)]
    Cmfd(#[from] CmfdError),
    #[error("not converged in {iterations} iterations (k = {k})")]
    NotConverged { iterations: usize, k: f64 },
}

/// `materials[i]` は領域 i の物質。
#[tracing::instrument(name = "moc", skip_all, fields(regions = materials.len(), tracks = tracking.tracks.len()))]
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
    // どの群も正（ふつうの物質）か、どの群も 0（真空）であること
    if materials
        .iter()
        .any(|m| !m.is_void() && m.total.iter().any(|&t| t <= 0.0 || !t.is_finite()))
    {
        return Err(MocError::NonPositiveTotal);
    }
    let order = options.scattering_order;
    if let Some((region, m)) = materials
        .iter()
        .enumerate()
        .find(|(_, m)| !m.is_void() && m.anisotropic.len() < order)
    {
        return Err(MocError::ScatteringOrder {
            region,
            available: m.anisotropic.len(),
            requested: order,
        });
    }
    if let Some(g) = group_without_source(materials, ng) {
        return Err(MocError::GroupWithoutSource(g));
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
    // 全体が真空のセルは、粗メッシュ（CMFD）の拡散係数が定まらない
    let lattice = &tracking.lattice;
    if let Some(cell) = (0..lattice.cells_per_side().pow(2))
        .find(|&c| lattice.cell_regions(c).all(|r| materials[r].is_void()))
    {
        return Err(MocError::VoidCell(cell));
    }
    let volumes = tracking.volumes(regions);
    if let Some(i) = volumes.iter().position(|&v| v <= 0.0) {
        return Err(MocError::Untracked(i));
    }
    // 領域の平均はトラックから求めた体積で計算するので、細い領域でトラックが足りないと
    // 解析値からずれる。分割とトラック間隔が釣り合っているかを確かめられるように出す
    if tracing::enabled!(tracing::Level::DEBUG) {
        let lattice = &tracking.lattice;
        let exact = lattice.volumes();
        let (region, error) = worst_volume_error(&exact, &volumes);
        let cell = lattice.region_cell(region);
        tracing::debug!(
            regions,
            worst_region = region,
            cell,
            relative_error = error,
            "track volumes against the exact areas"
        );
        for (region, (exact, tracked)) in exact.iter().zip(&volumes).enumerate() {
            tracing::trace!(
                region,
                exact,
                tracked,
                relative_error = tracked / exact - 1.0
            );
        }
    }

    let mesh = CoarseMesh::new(tracking);
    // 角度中性子束のモーメント (l, m)。0 番目がスカラー中性子束
    let pairs = planar_moments(order);
    let sweeper = Sweeper::new(tracking, materials, ng, &mesh, &pairs);
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
    // [領域][モーメント (1 番目から)][群] の高次のモーメント
    let mut higher = vec![vec![0.0; (pairs.len() - 1) * ng]; regions];
    let initial = production(&flux);
    if initial <= 0.0 {
        return Err(MocError::NoFission);
    }
    for f in flux.iter_mut().flatten() {
        *f /= initial;
    }
    let mut k = 1.0;
    // スイープの終わりに出ていく角度中性子束（`Sweeper::sweep` を参照）。最初は端の領域の中性子束の
    // 等方な分布とする（0 から始めると、鏡面反射の辺でも外向きの流れがあるように見え、CMFD が
    // 漏れの大きな問題として解いてしまう）
    let mut outgoing = sweeper.isotropic_outgoing(&flux);
    let mut next_outgoing = outgoing.clone();
    let mut tallies = sweeper.tallies(regions);
    let mut k_history = Remaining::default();
    let mut flux_history = Remaining::default();

    for iteration in 1..=options.max_iterations {
        // 線源を Σt で割ったもの（角度中性子束の平衡値）の [モーメント][群] の展開係数
        let equilibrium: Vec<Vec<f64>> = (0..regions)
            .map(|i| equilibrium_moments(materials[i], &flux[i], &higher[i], k, &pairs))
            .collect();

        sweeper.sweep(&equilibrium, &outgoing, &mut next_outgoing, &mut tallies)?;
        std::mem::swap(&mut outgoing, &mut next_outgoing);
        let Sums {
            delta,
            void_flux,
            mut leak,
            mut current,
        } = sum_tallies(&tallies, ng);

        let mut new_flux: Vec<Vec<f64>> = (0..regions)
            .map(|i| {
                // 真空では、通り抜ける角度中性子束の経路長の積分から求める
                if let Some(v) = sweeper.void[i] {
                    return void_flux[v * ng..][..ng]
                        .iter()
                        .map(|f| f / volumes[i])
                        .collect();
                }
                let total = &materials[i].total;
                let delta = &delta[i * pairs.len() * ng..][..ng];
                (0..ng)
                    .map(|g| 4.0 * PI * equilibrium[i][g] + delta[g] / (total[g] * volumes[i]))
                    .collect()
            })
            .collect();
        // 高次のモーメント: φ_k = 4π e_k / (2l + 1) + Σ w Y_k Δψ / (Σt V)（真空では散乱しないので使わない）
        let mut new_higher: Vec<Vec<f64>> = (0..regions)
            .map(|i| {
                let total = &materials[i].total;
                let mut values = vec![0.0; (pairs.len() - 1) * ng];
                if sweeper.void[i].is_some() {
                    return values;
                }
                for (k, &(l, _)) in pairs.iter().enumerate().skip(1) {
                    let delta = &delta[(i * pairs.len() + k) * ng..][..ng];
                    let e = &equilibrium[i][k * ng..][..ng];
                    for (g, v) in values[(k - 1) * ng..][..ng].iter_mut().enumerate() {
                        *v = 4.0 * PI * e[g] / (2 * l + 1) as f64
                            + delta[g] / (total[g] * volumes[i]);
                    }
                }
                values
            })
            .collect();
        // べき乗法: 生成率の比で k を更新し、中性子束と境界の値を生成率 1 に規格化する
        let rate = production(&new_flux);
        let mut new_k = k * rate;
        for f in new_flux.iter_mut().chain(&mut new_higher).flatten() {
            *f /= rate;
        }
        for l in leak.iter_mut().chain(&mut current) {
            *l /= rate;
        }
        for b in &mut outgoing {
            *b /= rate;
        }
        // CMFD: 粗メッシュの固有値問題の解で、k と、セルごと・群ごとの中性子束の比を更新する。
        // 粗メッシュの解は生成率 1 に規格化されているので、更新後も生成率は 1 のまま
        if let Acceleration::Cmfd(cmfd) = &options.acceleration {
            let coarse = Coarse::homogenize(&mesh, materials, &volumes, &new_flux, &current)?;
            let (coarse_k, ratio) = coarse.solve(new_k, cmfd)?;
            new_k = coarse_k;
            // 高次のモーメントも同じ比で（角度分布の形は変えない）
            for (region, (f, h)) in new_flux.iter_mut().zip(&mut new_higher).enumerate() {
                let cell = mesh.cell(region);
                for (j, value) in f.iter_mut().chain(h.iter_mut()).enumerate() {
                    *value *= ratio[cell * ng + j % ng];
                }
            }
            // トラックの端の角度中性子束も、端のあるセルの比で更新する（領域の中性子束だけを
            // 更新すると、次のスイープで入ってくる流れと食い違う）
            sweeper.scale_outgoing(&mut outgoing, |region, g| ratio[mesh.cell(region) * ng + g]);
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
        tracing::debug!(
            iteration,
            k = new_k,
            k_change = new_k - k,
            flux_change = change
        );
        if tracing::enabled!(tracing::Level::TRACE) {
            let (index, largest) = new_flux
                .iter()
                .flatten()
                .zip(flux.iter().flatten())
                .map(|(new, old)| relative_change(*new, *old))
                .enumerate()
                .fold((0, 0.0), |a, (i, c)| if c > a.1 { (i, c) } else { a });
            let (region, group) = (index / ng, index % ng);
            tracing::trace!(
                region,
                cell = mesh.cell(region),
                group,
                old = flux[region][group],
                new = new_flux[region][group],
                change = largest,
                "largest flux change"
            );
        }
        let k_remaining = k_history.update((new_k - k).abs());
        let flux_remaining = flux_history.update(change);
        let converged =
            k_remaining < options.k_tolerance && flux_remaining < options.flux_tolerance;
        flux = new_flux;
        higher = new_higher;
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

/// トラックのスイープ。トラックと物質から決まる量（方向ごとの球面調和関数など）を、反復の前に一度だけ求めて持つ。
struct Sweeper<'a> {
    tracking: &'a Tracking,
    /// 正味の中性子流を数える粗メッシュ
    mesh: &'a CoarseMesh,
    ng: usize,
    /// [極角][群] の数
    width: usize,
    /// [領域][群] の全断面積。減衰 1 - exp(-τ) はスイープのたびにその場で計算する
    totals: Vec<f64>,
    /// 方位角ごとのトラックの番号の範囲
    azimuths: Vec<std::ops::Range<usize>>,
    /// [領域] の真空の領域の通し番号（真空でなければ `None`）。集計は真空の領域の分だけ持つ
    void: Vec<Option<usize>>,
    /// 角度中性子束のモーメントの数（1 なら等方散乱）
    moments: usize,
    /// [方位角][向き][極角][モーメント] の方向ごとの実球面調和関数の値
    harmonics: Vec<f64>,
    /// 非等方散乱のスイープの作業用の配列（再利用する。中身は方位角ごとに作り直すので、どれを使っても結果は同じ）
    pool: Mutex<Vec<Directional>>,
}

/// 非等方散乱のスイープで、方位角ごとに作り直す配列。どちらも [領域][向き][極角][群]。
struct Directional {
    /// 方向ごとの角度中性子束の平衡値 Σ_k Y_k(Ω) e_k
    source: Vec<f64>,
    /// 方向ごとの減衰の重み付き和
    delta: Vec<f64>,
}

impl<'a> Sweeper<'a> {
    fn new(
        tracking: &'a Tracking,
        materials: &[&Material],
        ng: usize,
        mesh: &'a CoarseMesh,
        pairs: &[(usize, isize)],
    ) -> Self {
        let width = POLAR_SIN.len() * ng;
        // 順方向は方位角 φ、逆方向は φ + π。z について対称なモーメントだけなので、極角は上半分だけでよい
        let mut harmonics = Vec::new();
        for azimuth in &tracking.azimuths {
            for side in 0..2 {
                let phi = azimuth.phi + side as f64 * PI;
                for sin in POLAR_SIN {
                    let mu = (1.0 - sin * sin).sqrt();
                    harmonics.extend(pairs.iter().map(|&(l, m)| real_harmonic(l, m, mu, phi)));
                }
            }
        }
        let totals = materials
            .iter()
            .flat_map(|m| m.total.iter().copied())
            .collect();
        Self {
            tracking,
            mesh,
            ng,
            width,
            totals,
            azimuths: (0..tracking.azimuths.len())
                .map(|a| tracking.azimuth_tracks(a))
                .collect(),
            void: materials
                .iter()
                .scan(0, |next, m| {
                    Some(m.is_void().then(|| {
                        *next += 1;
                        *next - 1
                    }))
                })
                .collect(),
            moments: pairs.len(),
            harmonics,
            pool: Mutex::new(Vec::new()),
        }
    }

    /// 方位角 `azimuth` の [向き][極角] ごとの実球面調和関数の値 [モーメント]
    fn harmonics(&self, azimuth: usize, side: usize, p: usize) -> &[f64] {
        let k = self.moments;
        &self.harmonics[((azimuth * 2 + side) * POLAR_SIN.len() + p) * k..][..k]
    }

    /// 方位角 `azimuth` の方向ごとの平衡値を展開係数 `equilibrium` から作り、減衰の和を 0 にした作業用の配列
    fn directional(&self, equilibrium: &[Vec<f64>], azimuth: usize) -> Directional {
        let (ng, width) = (self.ng, self.width);
        let size = equilibrium.len() * 2 * width;
        let reused = self
            .pool
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .pop();
        let mut d = reused.unwrap_or_else(|| Directional {
            source: vec![0.0; size],
            delta: vec![0.0; size],
        });
        d.delta.fill(0.0);
        for (region, e) in equilibrium.iter().enumerate() {
            // 真空の領域はスイープで線源を使わない
            if self.void[region].is_some() {
                continue;
            }
            for side in 0..2 {
                for p in 0..POLAR_SIN.len() {
                    let y = self.harmonics(azimuth, side, p);
                    let out =
                        &mut d.source[((2 * region + side) * POLAR_SIN.len() + p) * ng..][..ng];
                    for (g, o) in out.iter_mut().enumerate() {
                        *o = y.iter().enumerate().map(|(k, y)| y * e[k * ng + g]).sum();
                    }
                }
            }
        }
        d
    }

    /// 方向ごとの減衰の和を、モーメントごとの和 `delta` [領域][モーメント][群] に足し、作業用の配列を返す
    fn fold(&self, d: Directional, azimuth: usize, delta: &mut [f64]) {
        let (ng, k) = (self.ng, self.moments);
        for (region, sums) in delta.chunks_mut(k * ng).enumerate() {
            for side in 0..2 {
                for p in 0..POLAR_SIN.len() {
                    let y = self.harmonics(azimuth, side, p);
                    let values = &d.delta[((2 * region + side) * POLAR_SIN.len() + p) * ng..][..ng];
                    for (sum, y) in sums.chunks_mut(ng).zip(y) {
                        for (s, v) in sum.iter_mut().zip(values) {
                            *s += y * v;
                        }
                    }
                }
            }
        }
        self.pool
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .push(d);
    }

    /// 出ていく角度中性子束の配列の長さ
    fn outgoing_len(&self) -> usize {
        2 * self.tracking.tracks.len() * self.width
    }

    /// トラックのスイープの終わりの区間の領域
    fn exit_region(&self, track: usize, forward: bool) -> usize {
        let segments = &self.tracking.tracks[track].segments;
        if forward {
            segments[segments.len() - 1].region
        } else {
            segments[0].region
        }
    }

    /// スイープの終わりの領域の [領域][群] の中性子束 `flux` が等方に出ていくときの角度中性子束
    fn isotropic_outgoing(&self, flux: &[Vec<f64>]) -> Vec<f64> {
        let mut outgoing = vec![0.0; self.outgoing_len()];
        for (j, out) in outgoing.chunks_mut(self.width).enumerate() {
            let region = self.exit_region(j / 2, j % 2 == 0);
            for (o, f) in out.iter_mut().zip(flux[region].iter().cycle()) {
                *o = f / (4.0 * PI);
            }
        }
        outgoing
    }

    /// 出ていく角度中性子束に、スイープの終わりの領域と群ごとの係数 `factor` を掛ける
    fn scale_outgoing(&self, outgoing: &mut [f64], factor: impl Fn(usize, usize) -> f64) {
        for (j, out) in outgoing.chunks_mut(self.width).enumerate() {
            let region = self.exit_region(j / 2, j % 2 == 0);
            for (i, o) in out.iter_mut().enumerate() {
                *o *= factor(region, i % self.ng);
            }
        }
    }

    /// 方位角ごとの部分和の入れ物。方位角ごとに順にスイープし、方位角の順に足すので、結果は
    /// スレッドの数や仕事の分け方によらない。方位角ごとのトラックの長さの和はほぼ等しく、
    /// 仕事の量も揃う
    fn tallies(&self, regions: usize) -> Vec<Tally> {
        self.azimuths
            .iter()
            .map(|_| {
                let voids = self.void.iter().flatten().count();
                Tally::new(
                    regions * self.moments,
                    voids,
                    self.mesh.faces(),
                    self.ng,
                    self.width,
                )
            })
            .collect()
    }

    /// 方向ごとの重み 4π (w_a/2) ω_p sinθ_p δ_a
    fn direction_weight(&self, azimuth: usize, p: usize) -> f64 {
        let a = &self.tracking.azimuths[azimuth];
        4.0 * PI * 0.5 * a.weight * POLAR_WEIGHT[p] * POLAR_SIN[p] * a.spacing
    }

    /// 全トラックを両方向に1回ずつスイープする。
    ///
    /// `equilibrium` は [領域][モーメント][群] の角度中性子束の平衡値（線源 / Σt）の展開係数 e_k
    /// （平衡値は Σ_k Y_k(Ω) e_k。等方なら線源 / 4πΣt）。`tallies` の減衰の和は [領域][モーメント][群]。`incoming` は前回の
    /// スイープで出ていった角度中性子束、`outgoing` には今回出ていく分を書く。どちらも
    /// [トラック][向き (0: 順, 1: 逆)][極角][群] の並びで、各トラックは自分の場所にだけ書くので、
    /// トラックごとに並列にスイープできる。
    fn sweep(
        &self,
        equilibrium: &[Vec<f64>],
        incoming: &[f64],
        outgoing: &mut [f64],
        tallies: &mut [Tally],
    ) -> Result<(), CmfdError> {
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
            .enumerate()
            .try_for_each(|(azimuth, ((tally, outputs), range))| {
                tally.clear();
                // 非等方散乱では、線源と減衰の和を方向ごとに持つ
                let mut directional =
                    (self.moments > 1).then(|| self.directional(equilibrium, azimuth));
                for (index, out) in range.clone().zip(outputs.chunks_mut(2 * width)) {
                    let track = &self.tracking.tracks[index];
                    let weights: [f64; POLAR_SIN.len()] =
                        std::array::from_fn(|p| self.direction_weight(track.azimuth, p));
                    for forward in [true, false] {
                        let psi = &mut tally.psi;
                        match self.tracking.incoming_source(index, forward) {
                            Some(source) => psi.copy_from_slice(
                                &incoming[slot(source.track, source.forward) * width..][..width],
                            ),
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
                                add_current(&mut tally.current, path?.faces(), psi, &weights);
                                cell = next;
                            }
                            if let Some(v) = self.void[region] {
                                // 真空では角度中性子束は変わらない。3次元の経路長 (長さ / sinθ) を掛けて数える
                                let length = track.segments[s].length;
                                let v = &mut tally.void_flux[v * ng..][..ng];
                                for (p, weight) in weights.iter().enumerate() {
                                    let w = weight * length / POLAR_SIN[p];
                                    for (f, value) in v.iter_mut().zip(&psi[p * ng..][..ng]) {
                                        *f += w * value;
                                    }
                                }
                                continue;
                            }
                            let length = track.segments[s].length;
                            let total = &self.totals[region * ng..][..ng];
                            // 等方なら [群] をすべての極角で使い、非等方なら [極角][群]
                            let (q, d, stride) = match directional.as_mut() {
                                None => (
                                    &equilibrium[region][..ng],
                                    &mut tally.delta[region * ng..][..ng],
                                    0,
                                ),
                                Some(dir) => {
                                    let at = (2 * region + usize::from(!forward)) * width;
                                    (
                                        &dir.source[at..][..width],
                                        &mut dir.delta[at..][..width],
                                        ng,
                                    )
                                }
                            };
                            for (p, weight) in weights.iter().enumerate() {
                                let (q, d) = (&q[p * stride..][..ng], &mut d[p * stride..][..ng]);
                                let path = length / POLAR_SIN[p];
                                for g in 0..ng {
                                    let j = p * ng + g;
                                    let attenuation = -(-total[g] * path).exp_m1();
                                    let change = (psi[j] - q[g]) * attenuation;
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
                if let Some(d) = directional {
                    self.fold(d, azimuth, &mut tally.delta);
                }
                Ok(())
            })
    }
}

/// スイープで足し込む量（方位角ごとに持ち、スイープの後に足し合わせる）。
struct Tally {
    /// [領域][群] の角度中性子束の減衰の重み付き和
    delta: Vec<f64>,
    /// [真空の領域][群] の、真空の領域を通る角度中性子束の経路長の重み付き和
    void_flux: Vec<f64>,
    /// [群] の真空境界からの漏れ
    leak: Vec<f64>,
    /// [粗メッシュの面][群] の正味の中性子流（面の正の向きが正）
    current: Vec<f64>,
    /// スイープ中の角度中性子束 [極角][群]（作業用）
    psi: Vec<f64>,
}

impl Tally {
    fn new(regions: usize, voids: usize, faces: usize, groups: usize, width: usize) -> Self {
        Self {
            delta: vec![0.0; regions * groups],
            void_flux: vec![0.0; voids * groups],
            leak: vec![0.0; groups],
            current: vec![0.0; faces * groups],
            psi: vec![0.0; width],
        }
    }

    fn clear(&mut self) {
        self.delta.fill(0.0);
        self.void_flux.fill(0.0);
        self.leak.fill(0.0);
        self.current.fill(0.0);
    }
}

/// 方位角ごとの部分和を足し合わせたもの。
struct Sums {
    /// [領域][群] の角度中性子束の減衰の重み付き和
    delta: Vec<f64>,
    /// [真空の領域][群] の、真空の領域を通る角度中性子束の経路長の重み付き和
    void_flux: Vec<f64>,
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
        void_flux: sum(|t| &t.void_flux),
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

/// どの物質でも核分裂のスペクトルが 0 で、中性子の届く群からの散乱もない群（最初の1つ）。
/// その群の中性子束は常に 0 になり、断面積の誤りを疑うべき状態なので、解く前に止める。
fn group_without_source(materials: &[&Material], ng: usize) -> Option<usize> {
    let mut reached = vec![false; ng];
    for term in materials.iter().flat_map(|m| &m.fission_terms) {
        for (r, chi) in reached.iter_mut().zip(&term.chi) {
            *r |= *chi > 0.0;
        }
    }
    // 届く群からの散乱で届く群を、増えなくなるまで広げる
    let mut changed = true;
    while changed {
        changed = false;
        for m in materials {
            for (from, row) in m.scatter.iter().enumerate() {
                if !reached[from] {
                    continue;
                }
                for (to, s) in row.iter().enumerate() {
                    if *s > 0.0 && !reached[to] {
                        reached[to] = true;
                        changed = true;
                    }
                }
            }
        }
    }
    reached.iter().position(|r| !r)
}

/// トラックから求めた体積 `tracked` の、解析値 `exact` に対する相対誤差が最も大きい領域と、その誤差。
fn worst_volume_error(exact: &[f64], tracked: &[f64]) -> (usize, f64) {
    exact
        .iter()
        .zip(tracked)
        .map(|(e, t)| t / e - 1.0)
        .enumerate()
        .fold((0, 0.0), |worst, (region, error)| {
            if error.abs() > worst.1.abs() {
                (region, error)
            } else {
                worst
            }
        })
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
pub(crate) struct Remaining {
    previous: Option<f64>,
}

impl Remaining {
    /// 最新の変化量を受け取り、残りの誤差の見積もり（変化量そのもの以上）を返す。
    pub(crate) fn update(&mut self, change: f64) -> f64 {
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

/// 真空でないすべての物質が持つ散乱行列の最高次数（そういう物質がなければ 0）。
pub fn available_order(materials: &[&Material]) -> usize {
    materials
        .iter()
        .filter(|m| !m.is_void())
        .map(|m| m.anisotropic.len())
        .min()
        .unwrap_or(0)
}

/// 角度中性子束の平衡値（線源 / Σt）の [モーメント][群] の展開係数 e_k。線源は
/// q(Ω) = Σ_k (2l+1)/(4π) Y_k(Ω) Σ_g' σ_l(g'→g) φ_k,g'（と等方な核分裂源）なので、e_k = (2l+1) q_k / (4πΣt)。
/// `phi` はスカラー中性子束 [群]、`higher` は高次のモーメント [モーメント (1 番目から)][群]
fn equilibrium_moments(
    m: &Material,
    phi: &[f64],
    higher: &[f64],
    k: f64,
    pairs: &[(usize, isize)],
) -> Vec<f64> {
    let ng = phi.len();
    let mut e = equilibrium_flux(m, phi, k);
    for (index, &(l, _)) in pairs.iter().enumerate().skip(1) {
        let moment = &higher[(index - 1) * ng..][..ng];
        let mut source = vec![0.0; ng];
        // 真空には散乱行列がない（solve で、真空でない物質には次数の分があることを確かめている）
        if let Some(sigma) = m.anisotropic.get(l - 1) {
            for (from, row) in sigma.iter().enumerate() {
                for (to, s) in row.iter().enumerate() {
                    source[to] += s * moment[from];
                }
            }
        }
        let factor = (2 * l + 1) as f64 / (4.0 * PI);
        e.extend(source.iter().zip(&m.total).map(
            |(q, t)| {
                if *t == 0.0 { 0.0 } else { factor * q / t }
            },
        ));
    }
    e
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
    // 真空には線源がない（散乱も核分裂もないので source も 0）
    source
        .iter()
        .zip(&m.total)
        .map(|(q, t)| if *t == 0.0 { 0.0 } else { q / (4.0 * PI * t) })
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
            anisotropic: Vec::new(),
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

    /// 燃料・ギャップ・減速材のピンセルで、ギャップを `gap` にした解
    fn with_gap(gap: &Material) -> Solution {
        let cell = PinCell::new(1.26, vec![0.3, 0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.02).unwrap();
        let (f, m) = (fuel(), moderator());
        solve(&tracking, &[&f, gap, &m], &Options::default()).unwrap()
    }

    #[test]
    fn void_gap_matches_a_nearly_transparent_gap() {
        let void = with_gap(&Material::void(2));
        // 吸収も散乱もほとんどない物質（Σt だけがごく小さい）
        let thin = with_gap(&two_group(
            [1.0e-9, 1.0e-9],
            [[0.0, 0.0], [0.0, 0.0]],
            [0.0, 0.0],
        ));
        assert!((void.k - thin.k).abs() < 1e-7, "{} vs {}", void.k, thin.k);
        for g in 0..2 {
            let (a, b) = (void.flux[1][g], thin.flux[1][g]);
            assert!((a - b).abs() < 1e-6 * b, "group {g}: {a} vs {b}");
        }
    }

    #[test]
    fn rejects_partly_void_materials() {
        let cell = PinCell::new(1.26, vec![0.3, 0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.05).unwrap();
        let (f, m) = (fuel(), moderator());
        let partly = two_group([1.0, 0.0], [[0.5, 0.0], [0.0, 0.0]], [0.0, 0.0]);
        let result = solve(&tracking, &[&f, &partly, &m], &Options::default());
        assert!(
            matches!(result, Err(MocError::NonPositiveTotal)),
            "{result:?}"
        );
    }

    #[test]
    fn rejects_cells_that_are_all_void() {
        // ピンと、全体が真空のセル
        let pin = PinCell::new(1.26, vec![0.4]).unwrap();
        let empty = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, empty], vec![vec![0, 1], vec![1, 0]]).unwrap();
        let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
        let (f, m, v) = (fuel(), moderator(), Material::void(2));
        let result = solve(&tracking, &[&f, &m, &v, &v, &f, &m], &Options::default());
        assert!(matches!(result, Err(MocError::VoidCell(1))), "{result:?}");
    }

    #[test]
    fn scattering_without_total_is_not_void() {
        let cell = PinCell::new(1.26, vec![0.3, 0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.05).unwrap();
        let (f, m) = (fuel(), moderator());
        let broken = two_group([0.0, 0.0], [[0.1, 0.0], [0.0, 0.0]], [0.0, 0.0]);
        let result = solve(&tracking, &[&f, &broken, &m], &Options::default());
        assert!(
            matches!(result, Err(MocError::NonPositiveTotal)),
            "{result:?}"
        );
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

    /// P1〜P(L) の散乱行列を、P0 の行列に `factors` の係数を掛けて付けた物質
    fn with_moments(m: &Material, factors: &[f64]) -> Material {
        let mut m = m.clone();
        m.anisotropic = factors
            .iter()
            .map(|f| {
                m.scatter
                    .iter()
                    .map(|row| row.iter().map(|s| f * s).collect())
                    .collect()
            })
            .collect();
        m
    }

    fn order(scattering_order: usize) -> Options {
        Options {
            scattering_order,
            ..Options::default()
        }
    }

    /// 燃料のピンと水のセルを市松に並べた n×n の格子（外周は真空）
    fn vacuum_checkerboard(n: usize) -> Tracking {
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let layout = (0..n)
            .map(|i| (0..n).map(|j| (i + j) % 2).collect())
            .collect();
        let lattice = Lattice::new(vec![pin, water], layout).unwrap();
        let vacuum = Boundaries {
            left: Boundary::Vacuum,
            right: Boundary::Vacuum,
            top: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
        };
        Tracking::with_boundaries(&lattice, 16, 0.05, vacuum).unwrap()
    }

    /// 市松の格子の領域ごとの物質（ピンは燃料と水、水のセルは水）
    fn checkerboard_materials<'a>(n: usize, f: &'a Material, m: &'a Material) -> Vec<&'a Material> {
        (0..n * n)
            .flat_map(|c| {
                if (c / n + c % n).is_multiple_of(2) {
                    vec![f, m]
                } else {
                    vec![m]
                }
            })
            .collect()
    }

    #[test]
    fn isotropic_moments_reproduce_p0() {
        // 高次の散乱行列が 0 なら、高次のモーメントは線源に入らず、P0 と同じ解（丸め誤差まで）
        let tracking = vacuum_checkerboard(2);
        let (f, m) = (
            with_moments(&fuel(), &[0.0, 0.0]),
            with_moments(&moderator(), &[0.0, 0.0]),
        );
        let materials = checkerboard_materials(2, &f, &m);
        let p0 = solve(&tracking, &materials, &order(0)).unwrap();
        let p2 = solve(&tracking, &materials, &order(2)).unwrap();
        assert!((p2.k - p0.k).abs() < 1e-12 * p0.k, "{} vs {}", p2.k, p0.k);
        for (a, b) in p2.flux.iter().flatten().zip(p0.flux.iter().flatten()) {
            assert!((a - b).abs() < 1e-10 * b, "{a} vs {b}");
        }
    }

    #[test]
    fn anisotropy_does_not_change_an_infinite_medium() {
        // 角度中性子束が等方なので、高次の散乱は効かない
        let cell = PinCell::new(1.26, vec![0.3, 0.45]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
        let m = with_moments(&fuel(), &[0.4, 0.1, 0.02]);
        let solution = solve(&tracking, &[&m, &m, &m], &order(3)).unwrap();
        let expected = infinite_medium_k(&m);
        assert!(
            (solution.k - expected).abs() < 1e-6,
            "{} vs {expected}",
            solution.k
        );
    }

    #[test]
    fn forward_scattering_acts_like_the_transport_correction() {
        // 前方に偏った散乱では、外向きに流れる中性子が散乱しても外向きに進み続けるので、等方散乱より
        // 漏れが増えて k が下がる。その効果を拡散の近似で表すのが輸送補正（Σt と自群の散乱から
        // Σ_g' σ1(g→g') を引いた P0。拡散係数が大きくなる）なので、P1 の k は P0 より輸送補正の k に近い
        let n = 6;
        let tracking = vacuum_checkerboard(n);
        let (f, m) = (
            with_moments(&fuel(), &[0.05]),
            with_moments(&moderator(), &[0.5]),
        );
        let corrected = |m: &Material| {
            let mut c = m.clone();
            for g in 0..2 {
                let p1: f64 = m.anisotropic[0][g].iter().sum();
                c.total[g] -= p1;
                c.scatter[g][g] -= p1;
            }
            c.anisotropic.clear();
            c
        };
        let (fc, mc) = (corrected(&f), corrected(&m));
        let k = |f, m, options| {
            solve(&tracking, &checkerboard_materials(n, f, m), &options)
                .unwrap()
                .k
        };
        let (p0, p1, tc) = (
            k(&f, &m, order(0)),
            k(&f, &m, order(1)),
            k(&fc, &mc, order(0)),
        );
        assert!(p1 < p0, "P1 {p1} vs P0 {p0}");
        assert!(
            (p1 - tc).abs() < (p0 - tc).abs(),
            "P0 {p0}, P1 {p1}, transport-corrected {tc}"
        );
    }

    #[test]
    fn anisotropic_result_does_not_depend_on_the_thread_count() {
        let tracking = vacuum_checkerboard(2);
        let (f, m) = (
            with_moments(&fuel(), &[0.05, 0.01]),
            with_moments(&moderator(), &[0.5, 0.2]),
        );
        let materials = checkerboard_materials(2, &f, &m);
        let solve_on = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| solve(&tracking, &materials, &order(2)).unwrap())
        };
        let serial = solve_on(1);
        for threads in [2, 5] {
            assert_eq!(solve_on(threads), serial, "{threads} threads");
        }
    }

    #[test]
    fn available_order_is_the_lowest_among_materials() {
        let (f, m, void) = (
            with_moments(&fuel(), &[0.05, 0.01]),
            with_moments(&moderator(), &[0.5]),
            Material::void(2),
        );
        assert_eq!(available_order(&[&f, &void]), 2);
        assert_eq!(available_order(&[&f, &m, &void]), 1);
        assert_eq!(available_order(&[&f, &fuel()]), 0);
        assert_eq!(available_order(&[&void]), 0);
    }

    #[test]
    fn scattering_order_needs_the_moments() {
        let tracking = vacuum_checkerboard(2);
        let (f, m) = (with_moments(&fuel(), &[0.05]), moderator());
        let materials = checkerboard_materials(2, &f, &m);
        assert!(matches!(
            solve(&tracking, &materials, &order(1)),
            Err(MocError::ScatteringOrder {
                region: 1,
                available: 0,
                requested: 1
            })
        ));
        // 足りなければ使わない次数は無視する
        let m = with_moments(&moderator(), &[0.5]);
        let materials = checkerboard_materials(2, &f, &m);
        assert!(solve(&tracking, &materials, &order(1)).is_ok());
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
        let sweeper = Sweeper::new(&tracking, &materials, ng, &mesh, &[(0, 0)]);
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
        sweeper
            .sweep(&equilibrium, &incoming, &mut outgoing, &mut tallies)
            .unwrap();
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
    fn cmfd_accelerates_a_lattice_with_a_margin() {
        // 集合体の間の水のギャップのように、格子の周りに余白がある。縁のセルの粗メッシュは余白の分だけ広い
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let layout = vec![vec![0, 1, 0], vec![1, 0, 1], vec![0, 1, 0]];
        let lattice = Lattice::with_margin(vec![pin, water], layout, 0.3).unwrap();
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
        let (f, m) = (fuel(), moderator());
        // セルの順にピン（燃料と水）と水、その後に余白の 4 × 3 + 4 領域は水
        let mut materials = Vec::new();
        for row in [
            [true, false, true],
            [false, true, false],
            [true, false, true],
        ] {
            for pin in row {
                if pin {
                    materials.extend([&f, &m]);
                } else {
                    materials.push(&m);
                }
            }
        }
        materials.extend(std::iter::repeat_n(&m, 16));
        assert_eq!(materials.len(), lattice.regions());
        assert_cmfd_accelerates(&tracking, &materials);
    }

    #[test]
    fn cmfd_converges_with_an_optically_thick_reflector() {
        // C5G7 の反射体と同じく、熱群で光学的に厚い（Σt h ≈ 4）水のセルが真空境界まで続く。
        // 従来の CMFD はこのような粗メッシュで振動して収束しない
        let pin = PinCell::new(1.26, vec![0.54]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let mut layout = vec![vec![1; 6]; 6];
        for row in layout.iter_mut().take(3) {
            row[..3].fill(0);
        }
        let lattice = Lattice::new(vec![pin, water], layout).unwrap();
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&lattice, 8, 0.1, boundaries).unwrap();
        let f = fuel();
        let reflector = two_group([0.6, 3.2], [[0.5, 0.098], [0.002, 3.18]], [0.0, 0.0]);
        let mut materials = Vec::new();
        for cell in 0..36 {
            let (row, column) = (cell / 6, cell % 6);
            if row < 3 && column < 3 {
                materials.extend([&f, &reflector]);
            } else {
                materials.push(&reflector);
            }
        }
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
    fn group_without_neutrons_is_error() {
        // 3群目は χ = 0 で、どの群からも散乱で入らないので、中性子束が常に 0 になる
        let three_group = |nu_fission: [f64; 3]| {
            let chi = vec![1.0, 0.0, 0.0];
            Material {
                total: vec![0.5, 1.2, 0.8],
                scatter: vec![
                    vec![0.45, 0.02, 0.0],
                    vec![0.0, 0.9, 0.0],
                    vec![0.0, 0.0, 0.7],
                ],
                anisotropic: Vec::new(),
                fission: nu_fission.iter().map(|v| v / 2.5).collect(),
                nu_fission: nu_fission.to_vec(),
                fission_terms: if nu_fission.iter().any(|&v| v > 0.0) {
                    vec![FissionTerm {
                        chi,
                        nu_fission: nu_fission.to_vec(),
                    }]
                } else {
                    vec![]
                },
            }
        };
        let (f, m) = (three_group([0.008, 0.25, 0.0]), three_group([0.0; 3]));
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
        for acceleration in [Acceleration::None, Options::default().acceleration] {
            let options = Options {
                acceleration,
                ..Options::default()
            };
            let err = solve(&tracking, &[&f, &m], &options).unwrap_err();
            assert!(matches!(err, MocError::GroupWithoutSource(2)), "{err}");
        }
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
    fn worst_volume_error_finds_the_largest_relative_mismatch() {
        let exact = [1.0, 0.1, 2.0];
        let tracked = [1.001, 0.098, 2.002];
        let (region, error) = worst_volume_error(&exact, &tracked);
        assert_eq!(region, 1);
        assert!((error + 0.02).abs() < 1e-12, "{error}");
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
