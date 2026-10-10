//! 平坦線源の特性曲線法（MOC）による固有値計算。

use std::f64::consts::PI;
use std::sync::{Mutex, PoisonError};

use rayon::prelude::*;

pub mod dancoff;
pub mod fixed_source;
mod harmonics;
mod symmetric;

use crate::cmfd::{CmfdError, CmfdOptions, Coarse, CoarseMesh, FaceCells};
use crate::material::{Material, MaterialError};
use crate::symmetry::{Symmetry, SymmetryError};
use crate::tracking::{Boundary, Tracking};
use harmonics::{planar_moments, real_harmonic};
use symmetric::{RegionKind, Representatives};

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
    /// 解くときに使う格子の対称性（ADR-0112）。恒等だけなら全体を解く。どれかの操作で、形状、境界条件、
    /// トラック、移り合う領域の物質（`solve` に渡す同じ参照）が移った先と同じでなければエラー
    pub symmetry: Symmetry,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            scattering_order: 0,
            k_tolerance: 1.0e-7,
            flux_tolerance: 1.0e-6,
            max_iterations: 20_000,
            acceleration: Acceleration::Cmfd(CmfdOptions::default()),
            // 対称性は呼び出し側が決める。格子のプログラムは成り立つ対称操作をすべて使う（ADR-0112）
            symmetry: Symmetry::NONE,
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
    /// [群] の外周からの漏れ（真空の辺と、アルベドの辺で戻らない 1 − α の分。中性子束と同じ規格化）
    pub leakage: Vec<f64>,
    /// [外周の辺（左、右、上、下）][群] の辺の上の平均スカラー中性子束（中性子束と同じ規格化）。
    /// 集合体の不連続因子（ADF）に使う
    pub sides: Vec<Vec<f64>>,
    /// [外周の辺（左、右、上、下）][群] の辺の上で平均した、外向きの正味の中性子流（中性子束と同じ規格化）。
    /// CMFD のために数える粗メッシュの面の正味の流れを辺にわたって足し、辺の長さで割ったもの。
    /// ADF の1ノードの問題（ADR-0103）に使う
    pub currents: Vec<Vec<f64>>,
    /// [粗メッシュの面][群] の面の上の平均スカラー中性子束（中性子束と同じ規格化）。面の番号は
    /// `CoarseMesh`（格子のセルの面。縦の面を先に）。内側の面も含み、複数のノードの格子のノードの辺の
    /// ADF（ADR-0095）に使う
    pub face_fluxes: Vec<Vec<f64>>,
    /// [粗メッシュの面][群] の面の上で平均した、面の正の向き（縦の面は +x、横の面は -y）の正味の中性子流
    pub face_currents: Vec<Vec<f64>>,
    /// [領域][x, y][群] の領域平均の中性子流（角度中性子束の1次のモーメント。中性子束と同じ規格化）。散乱の次数が
    /// 1 以上なら1次のモーメントそのもの、0 ならスイープで数えた値。真空の領域は通り抜ける角度中性子束から数える。
    /// 非反射の格子の拡散係数の重み（ADR-0104）に使う
    pub region_currents: Vec<[Vec<f64>; 2]>,
    /// 格子の x、y の幅 (cm)
    pub size: [f64; 2],
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
    #[error("no track crosses cell face {0}; use a finer track spacing")]
    UntrackedFace(usize),
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
    #[error(transparent)]
    Symmetry(#[from] SymmetryError),
    #[error(
        "region {region} and its image {image} under the {operation} refer to different \
         materials (pass the same material for regions that the operation maps onto each \
         other; self-shielding does so for the operations it verifies, ADR-0114)"
    )]
    AsymmetricMaterial {
        operation: &'static str,
        region: usize,
        image: usize,
    },
    #[error("not converged in {iterations} iterations (k = {k})")]
    NotConverged { iterations: usize, k: f64 },
    #[error(
        "tolerances must be positive and finite, and CMFD ones no larger than MOC ones \
         (MOC k {k}, flux {flux}; CMFD k {cmfd_k}, flux {cmfd_flux})"
    )]
    Tolerances {
        k: f64,
        flux: f64,
        cmfd_k: f64,
        cmfd_flux: f64,
    },
}

/// `materials[i]` は領域 i の物質。
#[tracing::instrument(name = "moc", skip_all, fields(regions = materials.len(), tracks = tracking.tracks.len()))]
pub fn solve(
    tracking: &Tracking,
    materials: &[&Material],
    options: &Options,
) -> Result<Solution, MocError> {
    // CMFD を MOC より細かく解く前提（ADR-0110）。比が 1 以下なので、緩めた CMFD の精度はスイープによる変化を超えない
    let positive = |t: f64| t > 0.0 && t.is_finite();
    let (cmfd_k, cmfd_flux) = match &options.acceleration {
        Acceleration::Cmfd(c) => (c.k_tolerance, c.flux_tolerance),
        Acceleration::None => (options.k_tolerance, options.flux_tolerance),
    };
    if ![
        options.k_tolerance,
        options.flux_tolerance,
        cmfd_k,
        cmfd_flux,
    ]
    .into_iter()
    .all(positive)
        || cmfd_k > options.k_tolerance
        || cmfd_flux > options.flux_tolerance
    {
        return Err(MocError::Tolerances {
            k: options.k_tolerance,
            flux: options.flux_tolerance,
            cmfd_k,
            cmfd_flux,
        });
    }
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
    if let Some(cell) =
        (0..lattice.cells()).find(|&c| lattice.cell_regions(c).all(|r| materials[r].is_void()))
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
    let face_lengths = tracked_face_lengths(tracking, &mesh)?;
    if let Some(face) = face_lengths.iter().position(|&l| l <= 0.0) {
        return Err(MocError::UntrackedFace(face));
    }
    // 角度中性子束のモーメント (l, m)。0 番目がスカラー中性子束
    let pairs = planar_moments(order);
    let representatives = if options.symmetry == Symmetry::NONE {
        Representatives::all(tracking, &mesh, pairs.len())
    } else {
        Representatives::symmetric(tracking, materials, &mesh, &pairs, options.symmetry)?
    };
    let sweeper = Sweeper::new(tracking, materials, ng, &mesh, &pairs, representatives);
    // 状態（中性子束、モーメント、平衡値）は、領域の軌道の代表 `kept[i]` の分だけ持つ（対称性を使わなければ全ての領域）
    let reduced = &sweeper.representatives;
    let kept = &reduced.regions.representatives;
    let states = kept.len();
    let material = |i: usize| materials[kept[i]];
    let volume = |i: usize| volumes[kept[i]];
    // [代表] の軌道の領域の体積の和（生成率は全ての領域の和）
    let orbit_volumes = reduced.regions.sums(&volumes);
    let production = |flux: &[Vec<f64>]| -> f64 {
        flux.iter()
            .enumerate()
            .zip(&orbit_volumes)
            .map(|((i, phi), v)| {
                v * phi
                    .iter()
                    .zip(&material(i).nu_fission)
                    .map(|(f, n)| f * n)
                    .sum::<f64>()
            })
            .sum()
    };

    let mut flux = vec![vec![1.0; ng]; states];
    // [代表][モーメント (1 番目から)][群] の高次のモーメント
    let mut higher = vec![vec![0.0; (pairs.len() - 1) * ng]; states];
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
        let equilibrium: Vec<Vec<f64>> = (0..states)
            .map(|i| equilibrium_moments(material(i), &flux[i], &higher[i], k, &pairs))
            .collect();

        sweeper.sweep(&equilibrium, &outgoing, &mut next_outgoing, &mut tallies)?;
        std::mem::swap(&mut outgoing, &mut next_outgoing);
        let mut sums = sum_tallies(&sweeper, &tallies);
        sweeper
            .representatives
            .symmetrize(&mut sums, &sweeper.void, ng);
        let Sums {
            delta,
            void_flux,
            mut leak,
            mut current,
            face_flux,
            void_current,
            region_current,
        } = sums;

        let mut new_flux: Vec<Vec<f64>> = (0..states)
            .map(|i| {
                // 真空では、通り抜ける角度中性子束の経路長の積分から求める
                if let Some(v) = sweeper.void[kept[i]] {
                    return void_flux[v * ng..][..ng]
                        .iter()
                        .map(|f| f / volume(i))
                        .collect();
                }
                let total = &material(i).total;
                let delta = &delta[i * pairs.len() * ng..][..ng];
                (0..ng)
                    .map(|g| 4.0 * PI * equilibrium[i][g] + delta[g] / (total[g] * volume(i)))
                    .collect()
            })
            .collect();
        // 高次のモーメント: φ_k = 4π e_k / (2l + 1) + Σ w Y_k Δψ / (Σt V)（真空では散乱しないので使わない）
        let mut new_higher: Vec<Vec<f64>> = (0..states)
            .map(|i| {
                let total = &material(i).total;
                let mut values = vec![0.0; (pairs.len() - 1) * ng];
                if sweeper.void[kept[i]].is_some() {
                    return values;
                }
                for (k, &(l, _)) in pairs.iter().enumerate().skip(1) {
                    let delta = &delta[(i * pairs.len() + k) * ng..][..ng];
                    let e = &equilibrium[i][k * ng..][..ng];
                    for (g, v) in values[(k - 1) * ng..][..ng].iter_mut().enumerate() {
                        *v = 4.0 * PI * e[g] / (2 * l + 1) as f64
                            + delta[g] / (total[g] * volume(i));
                    }
                }
                values
            })
            .collect();
        // 領域の中性子流 [領域][x, y の順に群]: 真空は通り抜ける角度中性子束、散乱の次数が 1 以上なら1次の
        // モーメント（Y_11 = Ω_x、Y_1,-1 = Ω_y）、0 なら Σ w Ω Δψ / (Σt V)（等方な線源は中性子流を持たない）
        let moment_of = |pair: (usize, isize)| pairs.iter().position(|&q| q == pair);
        let (x_moment, y_moment) = (moment_of((1, 1)), moment_of((1, -1)));
        let mut new_currents: Vec<Vec<f64>> = (0..states)
            .map(|i| {
                if let Some(v) = sweeper.void[kept[i]] {
                    return void_current[v * 2 * ng..][..2 * ng]
                        .iter()
                        .map(|j| j / volume(i))
                        .collect();
                }
                match (x_moment, y_moment) {
                    (Some(x), Some(y)) => [x, y]
                        .iter()
                        .flat_map(|&k| new_higher[i][(k - 1) * ng..][..ng].to_vec())
                        .collect(),
                    _ => {
                        let total = &material(i).total;
                        region_current[i * 2 * ng..][..2 * ng]
                            .iter()
                            .enumerate()
                            .map(|(j, c)| c / (total[j % ng] * volume(i)))
                            .collect()
                    }
                }
            })
            .collect();
        // べき乗法: 生成率の比で k を更新し、中性子束と境界の値を生成率 1 に規格化する
        let rate = production(&new_flux);
        let mut new_k = k * rate;
        for f in new_flux.iter_mut().chain(&mut new_higher).flatten() {
            *f /= rate;
        }
        for c in new_currents.iter_mut().flatten() {
            *c /= rate;
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
            let cells = &reduced.cells;
            let coarse = Coarse::homogenize(
                &mesh,
                cells,
                materials,
                &volumes,
                &new_flux,
                &reduced.regions.index,
                &current,
            )?;
            // MOC がまだ収束から遠いうちは、CMFD を細かく解いても次のスイープで崩れる。スイープによる変化（CMFD の誤差を
            // 含まない）に合わせて精度を緩める
            let tuned = inexact_cmfd(
                cmfd,
                options,
                (new_k - k).abs(),
                largest_change(&new_flux, &flux),
            );
            // [代表のセル][群] の比
            let (coarse_k, ratio) = coarse.solve(new_k, &tuned)?;
            new_k = coarse_k;
            // 高次のモーメントも同じ比で（角度分布の形は変えない）
            for (i, ((f, h), c)) in new_flux
                .iter_mut()
                .zip(&mut new_higher)
                .zip(&mut new_currents)
                .enumerate()
            {
                let cell = cells.index[mesh.cell(kept[i])];
                for (j, value) in f
                    .iter_mut()
                    .chain(h.iter_mut())
                    .chain(c.iter_mut())
                    .enumerate()
                {
                    *value *= ratio[cell * ng + j % ng];
                }
            }
            // トラックの端の角度中性子束も、端のあるセルの比で更新する（領域の中性子束だけを
            // 更新すると、次のスイープで入ってくる流れと食い違う）
            sweeper.scale_outgoing(&mut outgoing, |region, g| {
                ratio[cells.index[mesh.cell(region)] * ng + g]
            });
        }
        let change = largest_change(&new_flux, &flux);
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
            let (i, group) = (index / ng, index % ng);
            tracing::trace!(
                region = kept[i],
                cell = mesh.cell(kept[i]),
                group,
                old = flux[i][group],
                new = new_flux[i][group],
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
            // 面と辺の上の平均: 長さで割る。入ってくる側は前の反復のスイープの値で、CMFD の比もかけていないが、
            // 収束した解ではどちらも反復の許容誤差の程度しか違わない
            let lattice = &tracking.lattice;
            let lengths = [
                lattice.height(),
                lattice.height(),
                lattice.width(),
                lattice.width(),
            ];
            // 面の中性子束は、トラックから推定した面の長さで割る（領域の体積と同じく、トラックの間隔の
            // 離散化で面を横切る本数が長さに比例しない分を打ち消す。外周の辺の全体では推定が厳密）
            let face_fluxes: Vec<Vec<f64>> = face_flux
                .chunks(ng)
                .zip(&face_lengths)
                .map(|(values, length)| values.iter().map(|v| v / (rate * length)).collect())
                .collect();
            let face_currents: Vec<Vec<f64>> = current
                .chunks(ng)
                .enumerate()
                .map(|(face, values)| {
                    let length = mesh.face_length(face);
                    values.iter().map(|v| v / length).collect()
                })
                .collect();
            return Ok(Solution {
                k,
                flux: reduced.expand(&flux),
                volumes,
                leakage: leak,
                sides: outer_sides(&mesh, &face_fluxes, 1.0, lengths),
                currents: outer_sides(&mesh, &face_currents, -1.0, lengths),
                face_fluxes,
                face_currents,
                region_currents: reduced.expand_vectors(&new_currents, ng),
                size: [lattice.width(), lattice.height()],
                iterations: iteration,
            });
        }
    }
    Err(MocError::NotConverged {
        iterations: options.max_iterations,
        k,
    })
}

/// 面の長さのトラックによる推定 [面]: トラックが面を横切るごとに、方位角の重み × トラックの間隔 / 方向が
/// 面の法線となす方位角の余弦（面を横切る割合を掛ける）を足す。領域の体積の推定（`Tracking::volumes`）と同じ
/// 形で、面の中性子束の推定と整合する。
fn tracked_face_lengths(tracking: &Tracking, mesh: &CoarseMesh) -> Result<Vec<f64>, MocError> {
    let mut lengths = vec![0.0; mesh.faces()];
    for (index, track) in tracking.tracks.iter().enumerate() {
        let azimuth = &tracking.azimuths[track.azimuth];
        let across = [azimuth.phi.cos().abs(), azimuth.phi.sin().abs()];
        let mut add = |faces: &[(usize, f64)]| {
            for &(face, share) in faces {
                let c = across[usize::from(!mesh.is_vertical(face))];
                lengths[face] += share.abs() * azimuth.weight * azimuth.spacing / c;
            }
        };
        let [start, end] = mesh.ends(index);
        add(&[start]);
        let mut cell = mesh.cell(track.segments[0].region);
        for segment in &track.segments {
            let next = mesh.cell(segment.region);
            if next != cell {
                add(mesh.path(cell, next)?.faces());
                cell = next;
            }
        }
        add(&[end]);
    }
    Ok(lengths)
}

/// [外周の辺（左、右、上、下）][群] の、外周の面の値 `values` [面][群]（面の上の平均）の辺の上の平均。
/// `lengths` は辺の長さ。`sign` が 1 なら値をそのまま（中性子束）、-1 なら外向きの向きにそろえて
/// （面の正の向きの正味の中性子流を、外向きの正味の中性子流にして）足す。
fn outer_sides(
    mesh: &CoarseMesh,
    values: &[Vec<f64>],
    sign: f64,
    lengths: [f64; 4],
) -> Vec<Vec<f64>> {
    let ng = values.first().map_or(0, Vec::len);
    let mut sums = vec![vec![0.0; ng]; 4];
    for (face, face_values) in values.iter().enumerate() {
        let FaceCells::Outer { cell, outward } = mesh.face_cells(face) else {
            continue;
        };
        // 外周の面は、その外向きの向きの辺の上にある
        let side = mesh
            .cell_faces(cell)
            .iter()
            .position(|&(f, s)| f == face && s == outward);
        if let Some(side) = side {
            let factor = if sign < 0.0 { outward } else { 1.0 };
            let length = mesh.face_length(face);
            for (sum, value) in sums[side].iter_mut().zip(face_values) {
                *sum += factor * value * length;
            }
        }
    }
    sums.iter()
        .zip(lengths)
        .map(|(values, length)| values.iter().map(|v| v / length).collect())
        .collect()
}

/// トラックのスイープ。方向や領域ごとに決まる量（球面調和関数、全断面積）を、反復の前に一度だけ求めて持つ。
/// 区間の数に比例する量（減衰 1 - exp(-τ) など）は持たない。集合体でトラックを細かくするとメモリに載らないため。
struct Sweeper<'a> {
    tracking: &'a Tracking,
    /// 正味の流れを数える粗メッシュ
    mesh: &'a CoarseMesh,
    ng: usize,
    /// [極角][群] の数
    width: usize,
    /// [領域][群] の全断面積。減衰 1 - exp(-τ) はスイープでトラックごとに計算する
    totals: Vec<f64>,
    /// [領域] の真空の領域の通し番号（真空でなければ `None`）。集計は真空の領域の分だけ持つ
    void: Vec<Option<usize>>,
    /// [方位角][縦の辺, 横の辺][極角] の、外周の辺の上の中性子束の重み
    side_weights: Vec<[[f64; POLAR_SIN.len()]; 2]>,
    /// 角度中性子束のモーメントの数（1 なら等方散乱）
    moments: usize,
    /// [方位角][向き][極角][モーメント] の方向ごとの実球面調和関数の値
    harmonics: Vec<f64>,
    /// 非等方散乱のスイープの作業用の配列（使い回す。中身は方位角ごとに作り直すので、どれを使っても結果は同じ）
    pool: Mutex<Vec<Directional>>,
    /// スイープするトラック（対称性を使わなければ全て）
    representatives: Representatives,
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
        representatives: Representatives,
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
            void: materials
                .iter()
                .scan(0, |next, m| {
                    Some(m.is_void().then(|| {
                        *next += 1;
                        *next - 1
                    }))
                })
                .collect(),
            side_weights: tracking
                .azimuths
                .iter()
                .map(|azimuth| {
                    // 辺の上の中性子束の重み: 流れの重み 4π (w_a/2) ω_p sinθ_p δ_a を sinθ_p と、方向が辺の法線と
                    // なす方位角の余弦（縦の辺（左右）では |cos φ|、横の辺（上下）では |sin φ|）で割る
                    [azimuth.phi.cos().abs(), azimuth.phi.sin().abs()].map(|c| {
                        std::array::from_fn(|p| {
                            4.0 * PI * 0.5 * azimuth.weight * POLAR_WEIGHT[p] * azimuth.spacing / c
                        })
                    })
                })
                .collect(),
            moments: pairs.len(),
            harmonics,
            pool: Mutex::new(Vec::new()),
            representatives,
        }
    }

    /// 面 `faces`（重みの絶対値が、その面を通る割合）を横切る角度中性子束 `psi` [極角][群] を、方位角の
    /// 重み `weights` [縦の面, 横の面][極角] で `face_flux` [面][群] に足す（面の中性子束 × 面の長さ）。
    fn add_face_flux(
        &self,
        face_flux: &mut [f64],
        faces: &[(usize, f64)],
        psi: &[f64],
        weights: &[[f64; POLAR_SIN.len()]; 2],
    ) {
        let ng = self.ng;
        for &(face, share) in faces {
            let weights = &weights[usize::from(!self.mesh.is_vertical(face))];
            let f = &mut face_flux[face * ng..][..ng];
            for (p, weight) in weights.iter().enumerate() {
                for (v, value) in f.iter_mut().zip(&psi[p * ng..][..ng]) {
                    *v += share.abs() * weight * value;
                }
            }
        }
    }

    /// 方位角 `azimuth` の [向き][極角] ごとの実球面調和関数の値 [モーメント]
    fn harmonics(&self, azimuth: usize, side: usize, p: usize) -> &[f64] {
        let k = self.moments;
        &self.harmonics[((azimuth * 2 + side) * POLAR_SIN.len() + p) * k..][..k]
    }

    /// 方位角 `azimuth` の方向ごとの平衡値を、代表の領域の展開係数 `equilibrium`（[代表][モーメント][群]）から作り、
    /// 減衰の和を 0 にした作業用の配列。組のスイープが通る領域 `regions`（組の中の番号の順）の分だけ持つ
    fn directional(
        &self,
        equilibrium: &[Vec<f64>],
        azimuth: usize,
        regions: &[usize],
    ) -> Directional {
        let (ng, width) = (self.ng, self.width);
        let size = regions.len() * 2 * width;
        let reused = self
            .pool
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .pop();
        let mut d = reused.unwrap_or_else(|| Directional {
            source: Vec::new(),
            delta: Vec::new(),
        });
        // 組ごとに領域の番号の付け方が違うので、前の組の値を残さない（真空の領域は 0 のまま）
        d.source.clear();
        d.source.resize(size, 0.0);
        d.delta.clear();
        d.delta.resize(size, 0.0);
        for (local, &region) in regions.iter().enumerate() {
            // 真空の領域はスイープで線源を使わない
            if self.void[region].is_some() {
                continue;
            }
            // 代表の値を、代表をこの領域へ移す操作で移す: e_k = s e'_k'
            let representatives = &self.representatives;
            let e = &equilibrium[representatives.regions.index[region]];
            let moments = representatives.moments(region);
            for side in 0..2 {
                for p in 0..POLAR_SIN.len() {
                    let y = self.harmonics(azimuth, side, p);
                    let out =
                        &mut d.source[((2 * local + side) * POLAR_SIN.len() + p) * ng..][..ng];
                    for (g, o) in out.iter_mut().enumerate() {
                        *o = y
                            .iter()
                            .zip(moments)
                            .map(|(y, &(k, sign))| y * (sign * e[k * ng + g]))
                            .sum();
                    }
                }
            }
        }
        d
    }

    /// 方向ごとの減衰の和を、モーメントごとの和 `delta` [組の中の領域][モーメント][群] に足し、作業用の配列を返す
    fn fold(&self, d: Directional, azimuth: usize, delta: &mut [f64]) {
        let (ng, k) = (self.ng, self.moments);
        for (local, sums) in delta.chunks_mut(k * ng).enumerate() {
            for side in 0..2 {
                for p in 0..POLAR_SIN.len() {
                    let y = self.harmonics(azimuth, side, p);
                    let values = &d.delta[((2 * local + side) * POLAR_SIN.len() + p) * ng..][..ng];
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

    /// 領域の中性子流を数えるか（等方散乱のとき。非等方では1次のモーメントを使う）
    fn counts_currents(&self) -> bool {
        self.moments == 1
    }

    /// 出ていく角度中性子束の配列（[代表][向き][極角][群]）の長さ
    fn outgoing_len(&self) -> usize {
        2 * self.representatives.tracks.len() * self.width
    }

    /// 代表の位置 `representative` のトラックのスイープの終わりの区間の領域
    fn exit_region(&self, representative: usize, forward: bool) -> usize {
        let track = self.representatives.tracks[representative];
        let segments = &self.tracking.tracks[track].segments;
        if forward {
            segments[segments.len() - 1].region
        } else {
            segments[0].region
        }
    }

    /// スイープの終わりの領域の中性子束（[代表][群] の `flux`）が等方に出ていくときの角度中性子束
    fn isotropic_outgoing(&self, flux: &[Vec<f64>]) -> Vec<f64> {
        let mut outgoing = vec![0.0; self.outgoing_len()];
        for (j, out) in outgoing.chunks_mut(self.width).enumerate() {
            let region = self.exit_region(j / 2, j % 2 == 0);
            let flux = &flux[self.representatives.regions.index[region]];
            for (o, f) in out.iter_mut().zip(flux.iter().cycle()) {
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

    /// 代表のトラックの組ごとの部分和の入れ物。組ごとに順にスイープし、組の順に足すので、結果は
    /// スレッドの数や仕事の分け方によらない。組は、対称性を使わなければ方位角ごと（方位角ごとのトラックの長さの和は
    /// ほぼ等しい）、使えば区間の数の和を揃えたもの
    fn tallies(&self, regions: usize) -> Vec<Tally> {
        let voids = self.void.iter().flatten().count();
        let representatives = &self.representatives;
        let tallies: Vec<Tally> = representatives
            .chunks
            .iter()
            .map(|range| {
                let tracks = &representatives.tracks[range.clone()];
                let segments = tracks
                    .iter()
                    .map(|&t| self.tracking.tracks[t].segments.len())
                    .max()
                    .unwrap_or(0);
                // 組のスイープが通る領域だけを持つ（対称性を使うと、組ごとに格子の一部しか通らない）
                let mut local = vec![NOT_CROSSED; regions];
                let mut crossed = Vec::new();
                for &t in tracks {
                    for s in &self.tracking.tracks[t].segments {
                        if local[s.region] == NOT_CROSSED {
                            local[s.region] = 0;
                            crossed.push(s.region);
                        }
                    }
                }
                crossed.sort_unstable();
                for (i, &r) in crossed.iter().enumerate() {
                    local[r] = i;
                }
                Tally::new(
                    crossed,
                    local,
                    self.moments,
                    voids,
                    self.mesh.faces(),
                    self.ng,
                    self.width,
                    segments,
                    self.counts_currents(),
                )
            })
            .collect();
        tracing::debug!(
            regions,
            crossed = tallies.iter().map(|t| t.regions.len()).sum::<usize>(),
            "regions held by the tallies of the track groups"
        );
        tallies
    }

    /// 全トラックを両方向に1回ずつスイープする。
    ///
    /// `equilibrium` は [代表][モーメント][群] の角度中性子束の平衡値（線源 / Σt）の展開係数 e_k
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
        let representatives = &self.representatives;
        // 出ていく角度中性子束の配列を、代表のトラックの組ごとの範囲に分ける
        let mut rest = outgoing;
        let mut outputs = Vec::with_capacity(representatives.chunks.len());
        for range in &representatives.chunks {
            let (head, tail) = rest.split_at_mut(2 * width * range.len());
            outputs.push(head);
            rest = tail;
        }
        tallies
            .par_iter_mut()
            .zip(outputs)
            .zip(&representatives.chunks)
            .try_for_each(|((tally, outputs), range)| {
                tally.clear();
                // 非等方散乱では、線源と減衰の和を方向ごとに持つ。組の中で方位角が変わるたびに作り直す
                let mut directional: Option<(usize, Directional)> = None;
                for (position, out) in range.clone().zip(outputs.chunks_mut(2 * width)) {
                    let index = representatives.tracks[position];
                    let track = &self.tracking.tracks[index];
                    if self.moments > 1
                        && directional
                            .as_ref()
                            .is_none_or(|(azimuth, _)| *azimuth != track.azimuth)
                    {
                        if let Some((azimuth, d)) = directional.take() {
                            self.fold(d, azimuth, &mut tally.delta);
                        }
                        directional = Some((
                            track.azimuth,
                            self.directional(equilibrium, track.azimuth, &tally.regions),
                        ));
                    }
                    // 代表のトラックは、移り合うトラックの数だけ重く数える（ADR-0113）
                    let multiplicity = representatives.multiplicity[position];
                    let weights: [f64; POLAR_SIN.len()] = std::array::from_fn(|p| {
                        multiplicity * direction_weight(self.tracking, track.azimuth, p)
                    });
                    let side_weights = &self.side_weights[track.azimuth]
                        .map(|weights| weights.map(|w| multiplicity * w));
                    // 減衰はトラックの両向きで同じなので、トラックごとに一度求める（真空の区間では使わない）
                    for (s, att) in track
                        .segments
                        .iter()
                        .zip(tally.attenuation.chunks_exact_mut(width))
                    {
                        if self.void[s.region].is_some() {
                            continue;
                        }
                        let total = &self.totals[s.region * ng..][..ng];
                        for (sin, att) in POLAR_SIN.iter().zip(att.chunks_exact_mut(ng)) {
                            let path = s.length / sin;
                            for (a, t) in att.iter_mut().zip(total) {
                                *a = -(-t * path).exp_m1();
                            }
                        }
                    }
                    for forward in [true, false] {
                        let psi = &mut tally.psi;
                        let exit_side = if forward {
                            track.end_side
                        } else {
                            track.start_side
                        };
                        match self.tracking.incoming_source(index, forward) {
                            // 出ていくときに境界から戻る分（アルベド）を掛けてある
                            Some(source) => psi.copy_from_slice(
                                &incoming
                                    [representatives.slot(source.track, source.forward) * width..]
                                    [..width],
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
                        // 進む向きの方位角の余弦と正弦（極角の sinθ は掛けていない）
                        let phi = self.tracking.azimuths[track.azimuth].phi;
                        let cosines = [direction * phi.cos(), direction * phi.sin()];
                        let entered = [(enter.0, direction * enter.1)];
                        add_current(&mut tally.current, &entered, psi, &weights);
                        self.add_face_flux(&mut tally.face_flux, &entered, psi, side_weights);
                        let n = track.segments.len();
                        let first = if forward { 0 } else { n - 1 };
                        let mut cell = self.mesh.cell(track.segments[first].region);
                        for step in 0..n {
                            let s = if forward { step } else { n - 1 - step };
                            let (region, length) =
                                (track.segments[s].region, track.segments[s].length);
                            // 隣のセルに移るときに通る面
                            let next = self.mesh.cell(region);
                            if next != cell {
                                let path = self.mesh.path(cell, next)?;
                                add_current(&mut tally.current, path.faces(), psi, &weights);
                                self.add_face_flux(
                                    &mut tally.face_flux,
                                    path.faces(),
                                    psi,
                                    side_weights,
                                );
                                cell = next;
                            }
                            if let Some(v) = self.void[region] {
                                // 真空では角度中性子束は変わらない。3次元の経路長 (長さ / sinθ) を掛けて数える
                                let flux = &mut tally.void_flux[v * ng..][..ng];
                                let current = &mut tally.void_current[v * 2 * ng..][..2 * ng];
                                for (p, weight) in weights.iter().enumerate() {
                                    let w = weight * length / POLAR_SIN[p];
                                    let cosines =
                                        [cosines[0] * POLAR_SIN[p], cosines[1] * POLAR_SIN[p]];
                                    for (g, value) in psi[p * ng..][..ng].iter().enumerate() {
                                        flux[g] += w * value;
                                        current[g] += w * cosines[0] * value;
                                        current[ng + g] += w * cosines[1] * value;
                                    }
                                }
                                continue;
                            }
                            let att = &tally.attenuation[s * width..][..width];
                            // 組の中の領域の番号
                            let local = tally.local[region];
                            // 等方なら [群] を全ての極角で使い、非等方なら [極角][群]
                            let (q, d, stride) = match directional.as_mut().map(|(_, d)| d) {
                                // スカラーの平衡値はどの操作でも変わらない
                                None => (
                                    &equilibrium[representatives.regions.index[region]][..ng],
                                    &mut tally.delta[local * ng..][..ng],
                                    0,
                                ),
                                Some(dir) => {
                                    let at = (2 * local + usize::from(!forward)) * width;
                                    (
                                        &dir.source[at..][..width],
                                        &mut dir.delta[at..][..width],
                                        ng,
                                    )
                                }
                            };
                            // 等方散乱では、中性子流も減衰の和から数える（非等方では1次のモーメントを使う）
                            let mut current = (stride == 0)
                                .then(|| &mut tally.region_current[local * 2 * ng..][..2 * ng]);
                            for (p, weight) in weights.iter().enumerate() {
                                let (q, d) = (&q[p * stride..][..ng], &mut d[p * stride..][..ng]);
                                let cosines =
                                    [cosines[0] * POLAR_SIN[p], cosines[1] * POLAR_SIN[p]];
                                for g in 0..ng {
                                    let j = p * ng + g;
                                    let change = (psi[j] - q[g]) * att[j];
                                    psi[j] -= change;
                                    d[g] += weight * change;
                                    if let Some(c) = current.as_deref_mut() {
                                        c[g] += weight * cosines[0] * change;
                                        c[ng + g] += weight * cosines[1] * change;
                                    }
                                }
                            }
                        }
                        let exited = [(exit.0, direction * exit.1)];
                        add_current(&mut tally.current, &exited, psi, &weights);
                        self.add_face_flux(&mut tally.face_flux, &exited, psi, side_weights);
                        // 境界から戻るのは割合 r（鏡面反射 1、アルベド α、真空 0）で、残りは漏れ。戻る分を
                        // 掛けて渡すので、つながった先のトラックはそのまま受け取る
                        let r = returned(self.tracking.boundaries().of(exit_side));
                        if r < 1.0 {
                            for (p, weight) in weights.iter().enumerate() {
                                for (l, value) in tally.leak.iter_mut().zip(&psi[p * ng..][..ng]) {
                                    *l += (1.0 - r) * weight * value;
                                }
                            }
                        }
                        let out = &mut out[usize::from(!forward) * width..][..width];
                        for (o, v) in out.iter_mut().zip(psi.iter()) {
                            *o = r * v;
                        }
                    }
                }
                if let Some((azimuth, d)) = directional {
                    self.fold(d, azimuth, &mut tally.delta);
                }
                Ok(())
            })
    }
}

/// 組のスイープが通らない領域の、組の中の番号
const NOT_CROSSED: usize = usize::MAX;

/// スイープで足し込む量（代表のトラックの組ごとに持ち、スイープの後に足し合わせる）。領域ごとの量は、組のスイープが通る領域
/// `regions` の分だけ持つ。
struct Tally {
    /// 組のスイープが通る領域（番号の順）
    regions: Vec<usize>,
    /// [領域] の組の中の番号（通らなければ `NOT_CROSSED`）
    local: Vec<usize>,
    /// [組の中の領域][モーメント][群] の角度中性子束の減衰の重み付き和
    delta: Vec<f64>,
    /// [真空の領域][群] の、真空の領域を通る角度中性子束の経路長の重み付き和
    void_flux: Vec<f64>,
    /// [群] の真空境界からの漏れ
    leak: Vec<f64>,
    /// [粗メッシュの面][群] の正味の流れ（面の正の向きが正）
    current: Vec<f64>,
    /// [粗メッシュの面][群] の、面を横切る角度中性子束の重み付き和（面の中性子束 × 面の長さ）
    face_flux: Vec<f64>,
    /// [真空の領域][x, y][群] の、真空の領域を通る角度中性子束 × 方向の余弦の経路長の重み付き和
    void_current: Vec<f64>,
    /// [組の中の領域][x, y][群] の、減衰の重み付き和 × 方向の余弦（等方散乱のときだけ。非等方では1次のモーメントを
    /// 使う）
    region_current: Vec<f64>,
    /// スイープ中の角度中性子束 [極角][群]（作業用）
    psi: Vec<f64>,
    /// スイープ中のトラックの [区間][極角][群] の減衰 1 - exp(-τ)（作業用。方位角で区間の最も多いトラックの分）
    attenuation: Vec<f64>,
}

impl Tally {
    /// `regions` は組のスイープが通る領域、`local` はその逆引き、`currents` は領域の中性子流を数えるか（等方散乱のとき）
    #[allow(clippy::too_many_arguments)]
    fn new(
        regions: Vec<usize>,
        local: Vec<usize>,
        moments: usize,
        voids: usize,
        faces: usize,
        groups: usize,
        width: usize,
        segments: usize,
        currents: bool,
    ) -> Self {
        let crossed = regions.len();
        Self {
            regions,
            local,
            delta: vec![0.0; crossed * moments * groups],
            void_flux: vec![0.0; voids * groups],
            void_current: vec![0.0; voids * 2 * groups],
            region_current: vec![0.0; if currents { crossed * 2 * groups } else { 0 }],
            leak: vec![0.0; groups],
            face_flux: vec![0.0; faces * groups],
            current: vec![0.0; faces * groups],
            psi: vec![0.0; width],
            attenuation: vec![0.0; segments * width],
        }
    }

    fn clear(&mut self) {
        self.delta.fill(0.0);
        self.void_flux.fill(0.0);
        self.void_current.fill(0.0);
        self.region_current.fill(0.0);
        self.leak.fill(0.0);
        self.face_flux.fill(0.0);
        self.current.fill(0.0);
    }
}

/// 方位角ごとの部分和を足し合わせたもの。
struct Sums {
    /// [代表][モーメント][群] の角度中性子束の減衰の重み付き和（全体をスイープしたときの、代表の領域の値）
    delta: Vec<f64>,
    /// [真空の領域][群] の、真空の領域を通る角度中性子束の経路長の重み付き和
    void_flux: Vec<f64>,
    /// [群] の真空境界からの漏れ
    leak: Vec<f64>,
    /// [粗メッシュの面][群] の正味の流れ
    current: Vec<f64>,
    /// [粗メッシュの面][群] の面の中性子束 × 面の長さ
    face_flux: Vec<f64>,
    /// [真空の領域][x, y][群] の中性子流の経路長の和
    void_current: Vec<f64>,
    /// [代表][x, y][群] の減衰 × 方向の余弦の和（等方散乱のときだけ）
    region_current: Vec<f64>,
}

/// 組ごとの部分和を組の順に足す。量の大きさはスイープ `sweeper` の代表の領域、真空の領域、面、モーメントの数で決まる
fn sum_tallies(sweeper: &Sweeper, tallies: &[Tally]) -> Sums {
    let ng = sweeper.ng;
    let representatives = &sweeper.representatives;
    let voids = sweeper.void.iter().flatten().count();
    let faces = sweeper.mesh.faces();
    // 領域ごとの量 [代表][幅] を、代表ごとに並列に足す（1つの代表の中では像の順、組の順なので、結果は分け方によらない）。
    // 対称性を使うときは、代表の領域に移る像の領域の値を集める
    let spread = |kind: RegionKind, width: usize, part: fn(&Tally) -> &[f64]| -> Vec<f64> {
        let mut sum = vec![0.0; representatives.regions.representatives.len() * width];
        if width == 0 {
            return sum;
        }
        let each = |from: usize, f: &mut dyn FnMut(&[f64])| {
            for tally in tallies {
                let local = tally.local[from];
                if local != NOT_CROSSED {
                    f(&part(tally)[local * width..][..width]);
                }
            }
        };
        sum.par_chunks_mut(width)
            .enumerate()
            .for_each(|(i, chunk)| representatives.gather(i, chunk, kind, &each));
        sum
    };
    // 要素ごとに並列に足す（どの組も長さ `len`）
    let sum = |len: usize, part: fn(&Tally) -> &[f64]| -> Vec<f64> {
        let mut sum = vec![0.0; len];
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
        delta: spread(RegionKind::Moments(ng), sweeper.moments * ng, |t| &t.delta),
        void_flux: sum(voids * ng, |t| &t.void_flux),
        leak: sum(ng, |t| &t.leak),
        current: sum(faces * ng, |t| &t.current),
        face_flux: sum(faces * ng, |t| &t.face_flux),
        void_current: sum(voids * 2 * ng, |t| &t.void_current),
        region_current: spread(
            RegionKind::Vector(ng),
            if sweeper.counts_currents() { 2 * ng } else { 0 },
            |t| &t.region_current,
        ),
    }
}

/// 境界から戻る角度中性子束の割合（鏡面反射 1、真空 0、アルベド α）
fn returned(boundary: Boundary) -> f64 {
    match boundary {
        Boundary::Reflective => 1.0,
        Boundary::Vacuum => 0.0,
        Boundary::Albedo(albedo) => albedo,
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

/// MOC の反復で CMFD を解く精度。スイープによる k と中性子束の変化（`k_change`, `flux_change`）が MOC の許容誤差の
/// 何倍あるかに比例して、設定の精度 `cmfd` を緩める（ADR-0110）。変化が許容誤差以下なら設定のまま。許容誤差は
/// `solve` の最初で正の有限と確かめてある
fn inexact_cmfd(cmfd: &CmfdOptions, moc: &Options, k_change: f64, flux_change: f64) -> CmfdOptions {
    CmfdOptions {
        k_tolerance: cmfd.k_tolerance * (k_change / moc.k_tolerance).max(1.0),
        flux_tolerance: cmfd.flux_tolerance * (flux_change / moc.flux_tolerance).max(1.0),
        ..*cmfd
    }
}

/// 領域と群のうちで最も大きい中性子束の相対変化
fn largest_change(new: &[Vec<f64>], old: &[Vec<f64>]) -> f64 {
    new.iter()
        .flatten()
        .zip(old.iter().flatten())
        .map(|(new, old)| relative_change(*new, *old))
        .fold(0.0, f64::max)
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

/// 真空でない全ての物質が持つ散乱行列の最高次数（そういう物質がなければ 0）。
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

/// トラックの1つの向きの、方向ごとの重み 4π (w_a/2) ω_p sinθ_p δ_a
fn direction_weight(tracking: &Tracking, azimuth: usize, p: usize) -> f64 {
    let a = &tracking.azimuths[azimuth];
    4.0 * PI * 0.5 * a.weight * POLAR_WEIGHT[p] * POLAR_SIN[p] * a.spacing
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Cell, Lattice, PinCell};
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
    fn sides_of_an_infinite_medium_see_the_mean_flux() {
        // 均質な無限媒質では中性子束は平らで等方なので、外周の辺の上の中性子束は領域の中性子束と同じ
        let cell = PinCell::new(1.26, vec![0.3, 0.45]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
        let m = fuel();
        let solution = solve(&tracking, &[&m, &m, &m], &Options::default()).unwrap();
        assert_eq!(solution.sides.len(), 4);
        for side in &solution.sides {
            for (s, f) in side.iter().zip(&solution.flux[0]) {
                assert!(
                    (s - f).abs() < 1e-8 * f,
                    "{side:?} vs {:?}",
                    solution.flux[0]
                );
            }
        }
    }

    #[test]
    fn faces_of_an_infinite_medium_see_the_mean_flux() {
        // 均質な無限媒質では、内側の面（セルの境）の上の中性子束も領域の中性子束と同じで、正味の流れは 0
        let cell = PinCell::new(1.26, vec![0.3, 0.45]).unwrap();
        let lattice = Lattice::new(vec![cell], vec![vec![0, 0], vec![0, 0]]).unwrap();
        let tracking = Tracking::new(&lattice, 8, 0.1).unwrap();
        let m = fuel();
        let solution = solve(&tracking, &[&m; 12], &Options::default()).unwrap();
        // 2×2 のセルの面は 縦 2×3 と横 3×2
        assert_eq!(solution.face_fluxes.len(), 12);
        assert_eq!(solution.face_currents.len(), 12);
        let mean = &solution.flux[0];
        for (face, values) in solution.face_fluxes.iter().enumerate() {
            for (g, (v, f)) in values.iter().zip(mean).enumerate() {
                assert!(
                    (v - f).abs() < 1e-8 * f,
                    "face {face} group {g}: {v} vs {f}"
                );
                let j = solution.face_currents[face][g];
                assert!(j.abs() < 1e-8 * f, "face {face} group {g}: current {j}");
            }
        }
    }

    #[test]
    fn faces_agree_with_sides_and_mirror_lattices() {
        // 外周の面を辺ごとに長さで重みづけして平均すると、辺の値（`sides`、`currents`）に一致する。
        // 真ん中が対称面の格子の、真ん中の内側の面の中性子束は、半分の格子の鏡面反射の辺の値に近い
        // （方位角がセルの比で変わる cyclic tracking の差だけ違う）
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let open = Boundaries {
            left: Boundary::Vacuum,
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let full = Lattice::new(
            vec![pin.clone(), water.clone()],
            vec![
                vec![1, 0, 0, 1],
                vec![1, 0, 0, 1],
                vec![1, 0, 0, 1],
                vec![1, 0, 0, 1],
            ],
        )
        .unwrap();
        let (f, m) = (fuel(), moderator());
        // 行ごとの領域: 水、ピン（燃料、減速材）、ピン、水
        let row = [&m, &f, &m, &f, &m, &m];
        let materials: Vec<&Material> = (0..4).flat_map(|_| row).collect();
        let tracking = Tracking::with_boundaries(&full, 32, 0.02, open).unwrap();
        let s = solve(&tracking, &materials, &Options::default()).unwrap();
        // 縦の面 (r, i) は r × 5 + i、横の面 (j, c) は 20 + j × 4 + c
        let vertical = |r: usize, i: usize| r * 5 + i;
        let horizontal = |j: usize, c: usize| 20 + j * 4 + c;
        let sides: [Vec<(usize, f64)>; 4] = [
            (0..4).map(|r| (vertical(r, 0), -1.0)).collect(),
            (0..4).map(|r| (vertical(r, 4), 1.0)).collect(),
            (0..4).map(|c| (horizontal(0, c), -1.0)).collect(),
            (0..4).map(|c| (horizontal(4, c), 1.0)).collect(),
        ];
        for (side, faces) in sides.iter().enumerate() {
            for g in 0..2 {
                let flux: f64 = faces.iter().map(|&(i, _)| s.face_fluxes[i][g]).sum::<f64>() / 4.0;
                let out: f64 = faces
                    .iter()
                    .map(|&(i, o)| o * s.face_currents[i][g])
                    .sum::<f64>()
                    / 4.0;
                let (expected, current) = (s.sides[side][g], s.currents[side][g]);
                assert!(
                    (flux - expected).abs() < 1e-12 * expected,
                    "{flux} vs {expected}"
                );
                assert!(
                    (out - current).abs() < 1e-12 * expected,
                    "{out} vs {current}"
                );
            }
        }
        let half = Lattice::rectilinear(
            vec![1.26; 2],
            vec![1.26; 4],
            vec![Cell::Pin(pin), Cell::Pin(water)],
            vec![vec![1, 0]; 4],
        )
        .unwrap();
        let tracking = Tracking::with_boundaries(
            &half,
            32,
            0.02,
            Boundaries {
                left: Boundary::Vacuum,
                ..Boundaries::REFLECTIVE
            },
        )
        .unwrap();
        let materials: Vec<&Material> = (0..4).flat_map(|_| [&m, &f, &m]).collect();
        let h = solve(&tracking, &materials, &Options::default()).unwrap();
        for g in 0..2 {
            let center: f64 = (0..4)
                .map(|r| s.face_fluxes[vertical(r, 2)][g])
                .sum::<f64>()
                / 4.0;
            // 全体の中性子束の規格化（生成率 1）は半分の格子の 2 倍の体積にわたるので、半分の値の 1/2
            let expected = 0.5 * h.sides[1][g];
            assert!(
                (center - expected).abs() < 2e-3 * expected,
                "group {g}: {center} vs {expected}"
            );
            let net: f64 = (0..4)
                .map(|r| s.face_currents[vertical(r, 2)][g])
                .sum::<f64>();
            assert!(net.abs() < 1e-6 * center, "{net}");
        }
    }

    #[test]
    fn side_currents_add_up_to_the_leakage() {
        // 真空の辺（右、下）の外向きの正味の流れ × 辺の長さの和は漏れに等しく、鏡面反射の辺（左、上）
        // では 0
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
        let (f, m) = (fuel(), moderator());
        let s = solve(&tracking, &[&f, &m, &m, &m, &f, &m], &Options::default()).unwrap();
        assert_eq!(s.size, [2.52, 2.52]);
        assert_eq!(s.currents.len(), 4);
        for g in 0..2 {
            let leakage = s.leakage[g];
            assert!(leakage > 0.0);
            for side in [0, 2] {
                assert!(
                    s.currents[side][g].abs() * 2.52 < 1e-6 * leakage,
                    "{:?}",
                    s.currents
                );
            }
            for side in [1, 3] {
                assert!(s.currents[side][g] > 0.0, "{:?}", s.currents);
            }
            let out = 2.52 * (s.currents[1][g] + s.currents[3][g]);
            assert!((out - leakage).abs() < 1e-6 * leakage, "{out} vs {leakage}");
        }
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
    fn albedo_sides_lie_between_vacuum_and_reflection() {
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
        let (f, m) = (fuel(), moderator());
        let materials = [&f, &m, &m, &m, &f, &m];
        let solve_with = |right: Boundary| {
            let boundaries = Boundaries {
                right,
                ..Boundaries::REFLECTIVE
            };
            let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
            solve(&tracking, &materials, &Options::default()).unwrap()
        };
        let reflective = solve_with(Boundary::Reflective);
        let vacuum = solve_with(Boundary::Vacuum);
        // α = 1 は鏡面反射、α = 0 は真空と同じ
        assert!((solve_with(Boundary::Albedo(1.0)).k - reflective.k).abs() < 1e-9);
        assert!((solve_with(Boundary::Albedo(0.0)).k - vacuum.k).abs() < 1e-9);
        // 間のアルベドでは k も間で、返らなかった分は漏れとして収支が合う
        let half = solve_with(Boundary::Albedo(0.5));
        assert!(vacuum.k < half.k && half.k < reflective.k, "{}", half.k);
        let (mut absorption, mut production) = (0.0, 0.0);
        for (i, m) in materials.iter().enumerate() {
            for g in 0..2 {
                absorption += half.volumes[i] * m.absorption()[g] * half.flux[i][g];
                production += half.volumes[i] * m.nu_fission[g] * half.flux[i][g];
            }
        }
        let leakage: f64 = half.leakage.iter().sum();
        assert!(
            (absorption + leakage - production / half.k).abs() < 1e-6 * production / half.k,
            "{absorption} + {leakage} vs {}",
            production / half.k
        );
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
        // 正味の流れの和は打ち消し合う。線源と入ってくる角度中性子束は任意の値でよい
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
        let sweeper = Sweeper::new(
            &tracking,
            &materials,
            ng,
            &mesh,
            &[(0, 0)],
            Representatives::all(&tracking, &mesh, 1),
        );
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
        let sums = sum_tallies(&sweeper, &tallies);
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
    fn cmfd_tolerance_follows_the_moc_change() {
        let cmfd = CmfdOptions::default();
        let moc = Options::default();
        // MOC の許容誤差以下の変化では設定のまま
        assert_eq!(
            inexact_cmfd(&cmfd, &moc, 0.5 * moc.k_tolerance, moc.flux_tolerance),
            cmfd
        );
        // 許容誤差の 100 倍と 30 倍の変化では、k と中性子束の精度をそれぞれ同じ倍率で緩める
        let tuned = inexact_cmfd(
            &cmfd,
            &moc,
            100.0 * moc.k_tolerance,
            30.0 * moc.flux_tolerance,
        );
        assert!((tuned.k_tolerance / cmfd.k_tolerance - 100.0).abs() < 1e-12);
        assert!((tuned.flux_tolerance / cmfd.flux_tolerance - 30.0).abs() < 1e-12);
        assert_eq!(tuned.max_iterations, cmfd.max_iterations);
    }

    #[test]
    fn rejects_tolerances_that_do_not_bound_cmfd() {
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
        let (f, m) = (fuel(), moderator());
        let cmfd = CmfdOptions::default();
        let with = |k: f64, flux: f64, cmfd: CmfdOptions| Options {
            k_tolerance: k,
            flux_tolerance: flux,
            acceleration: Acceleration::Cmfd(cmfd),
            ..Options::default()
        };
        // MOC の許容誤差が 0、負、NaN、無限大
        let bad = [
            with(0.0, 1e-6, cmfd),
            with(1e-7, -1e-6, cmfd),
            with(f64::NAN, 1e-6, cmfd),
            with(1e-7, f64::INFINITY, cmfd),
            // CMFD の許容誤差が 0 か、MOC より大きい（k と中性子束で別に）
            with(
                1e-7,
                1e-6,
                CmfdOptions {
                    k_tolerance: 0.0,
                    ..cmfd
                },
            ),
            with(
                1e-7,
                1e-6,
                CmfdOptions {
                    k_tolerance: 1e-6,
                    ..cmfd
                },
            ),
            with(
                1e-7,
                1e-6,
                CmfdOptions {
                    flux_tolerance: 1e-5,
                    ..cmfd
                },
            ),
        ];
        for options in bad {
            let result = solve(&tracking, &[&f, &m], &options);
            assert!(
                matches!(result, Err(MocError::Tolerances { .. })),
                "{options:?}: {result:?}"
            );
        }
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

    /// 燃料の均質なブロック（1 cm 角、4 × 4 の領域）を x に `columns` 個並べた帯の、[領域] の x 座標（中心が原点）
    fn strip(columns: usize, n: usize) -> (Lattice, Vec<f64>) {
        let lattice = Lattice::rectilinear(
            vec![1.0; columns],
            vec![1.0],
            vec![Cell::Block { nx: n, ny: n }],
            vec![vec![0; columns]],
        )
        .unwrap();
        let x = (0..columns)
            .flat_map(|c| {
                (0..n * n).map(move |r| {
                    -0.5 * columns as f64 + c as f64 + ((r % n) as f64 + 0.5) / n as f64
                })
            })
            .collect();
        (lattice, x)
    }

    /// P1 の散乱を P0 の 0.3 倍で持つ物質
    fn with_p1(mut m: Material) -> Material {
        m.anisotropic = vec![
            m.scatter
                .iter()
                .map(|row| row.iter().map(|s| 0.3 * s).collect())
                .collect(),
        ];
        m
    }

    #[test]
    fn region_currents_vanish_in_an_infinite_medium() {
        let (lattice, _) = strip(2, 4);
        let tracking = Tracking::new(&lattice, 16, 0.05).unwrap();
        let m = fuel();
        let s = solve(&tracking, &[&m; 32], &Options::default()).unwrap();
        assert_eq!(s.region_currents.len(), 32);
        for (r, current) in s.region_currents.iter().enumerate() {
            for (c, component) in current.iter().enumerate() {
                for (g, j) in component.iter().enumerate() {
                    assert!(
                        j.abs() < 1e-8 * s.flux[r][g],
                        "region {r} {c} group {g}: {j}"
                    );
                }
            }
        }
    }

    /// 燃料の帯（x ∈ [−2, 2]、右が真空、ほかは鏡面反射。1 cm のブロックを n × n の領域に分け、1つを真空にする）の、
    /// x の1次のモーメントの式
    ///   Σ_r V_r J_x,r = 2 × 1 × J_右 − Σ_r x_r V_r (q_r − Σt,r φ_r)
    /// （部分積分と領域ごとのバランス ∫∇·J = V (q − Σt φ)）の [群] の相対差。x_r は領域の中心なので、領域の中での
    /// 変化による誤差（領域の幅の2乗に比例）が残る
    fn first_moment_mismatch(n: usize, order: usize) -> Vec<f64> {
        let (lattice, x) = strip(4, n);
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&lattice, 32, 0.02, boundaries).unwrap();
        let (f, void) = (with_p1(fuel()), Material::void(2));
        let materials: Vec<&Material> = (0..4 * n * n)
            .map(|r| if r == n + 1 { &void } else { &f })
            .collect();
        let options = Options {
            scattering_order: order,
            ..Options::default()
        };
        let s = solve(&tracking, &materials, &options).unwrap();
        (0..2)
            .map(|g| {
                let mut moment = 0.0;
                let mut divergence = 0.0;
                for (r, m) in materials.iter().enumerate() {
                    let v = s.volumes[r];
                    moment += v * s.region_currents[r][0][g];
                    let scattered: f64 = (0..2)
                        .map(|from| m.scatter[from][g] * s.flux[r][from])
                        .sum();
                    let fission: f64 = m
                        .fission_terms
                        .iter()
                        .map(|t| {
                            t.chi[g]
                                * t.nu_fission
                                    .iter()
                                    .zip(&s.flux[r])
                                    .map(|(n, phi)| n * phi)
                                    .sum::<f64>()
                        })
                        .sum();
                    let q = scattered + fission / s.k;
                    divergence += x[r] * v * (q - m.total[g] * s.flux[r][g]);
                }
                // 真空の領域の中性子流も右向き
                assert!(s.region_currents[n + 1][0][g] > 0.0);
                let expected = 2.0 * s.currents[1][g] - divergence;
                (moment - expected).abs() / expected.abs()
            })
            .collect()
    }

    #[test]
    fn region_currents_satisfy_the_first_moment_balance() {
        // 等方散乱（スイープで数えた中性子流）と P1 の散乱（1次のモーメント）の両方で、領域を細かくすると差は
        // 幅の2乗で減る（0.25 cm で熱群 2%、0.125 cm で 0.6% 程度）
        for order in [0, 1] {
            let (coarse, fine) = (
                first_moment_mismatch(4, order),
                first_moment_mismatch(8, order),
            );
            for g in 0..2 {
                assert!(fine[g] < 1e-2, "order {order} group {g}: {fine:?}");
                assert!(
                    coarse[g] > 3.0 * fine[g],
                    "order {order} group {g}: {coarse:?} {fine:?}"
                );
            }
        }
    }

    #[test]
    fn region_currents_are_antisymmetric_in_a_mirror_lattice() {
        // 両側が真空の対称な帯では、x の中性子流は左右で符号が逆、y の中性子流は 0（領域ごとには推定の誤差がある）
        let (lattice, _) = strip(4, 4);
        let boundaries = Boundaries {
            left: Boundary::Vacuum,
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&lattice, 32, 0.02, boundaries).unwrap();
        let m = fuel();
        let s = solve(&tracking, &[&m; 64], &Options::default()).unwrap();
        let mirror = |r: usize| {
            let (cell, inner) = (r / 16, r % 16);
            (3 - cell) * 16 + (inner / 4) * 4 + 3 - inner % 4
        };
        let scale = s
            .region_currents
            .iter()
            .map(|c| c[0][0].abs())
            .fold(0.0, f64::max);
        for r in 0..64 {
            for g in 0..2 {
                let (a, b) = (
                    s.region_currents[r][0][g],
                    s.region_currents[mirror(r)][0][g],
                );
                assert!(
                    (a + b).abs() < 1e-6 * scale,
                    "region {r} group {g}: {a} vs {b}"
                );
                assert!(s.region_currents[r][1][g].abs() < 1e-2 * scale);
            }
        }
        // 1D の問題なので、セルの中の y の中性子流は（領域の推定の誤差で ±0.2% 程度あるが）セルで足すと打ち消す
        for cell in 0..4 {
            for g in 0..2 {
                let sum: f64 = (cell * 16..(cell + 1) * 16)
                    .map(|r| s.volumes[r] * s.region_currents[r][1][g])
                    .sum();
                assert!(sum.abs() < 1e-8 * scale, "cell {cell} group {g}: {sum}");
            }
        }
    }

    #[test]
    fn faces_of_a_mirror_symmetric_lattice_are_symmetric() {
        // 格子の中心（4 つのセルの角）や外周のセルの角をちょうど通るトラックがある、上下と左右で対称な格子。角を通る移動が
        // 丸めで片側の隣のセルを経由したことになると、鏡像の面の量が揃わない（#0051）
        let pin = PinCell::with_sectors(1.26, vec![0.35, 0.4096], 8).unwrap();
        let water = PinCell::grid(1.26, 2).unwrap();
        let layout = vec![
            vec![0, 1, 1, 0],
            vec![0, 0, 0, 0],
            vec![0, 0, 0, 0],
            vec![0, 1, 1, 0],
        ];
        let lattice = Lattice::with_margin(vec![pin, water], layout, 0.04).unwrap();
        let (f, w) = (fuel(), moderator());
        let materials: Vec<&Material> = (0..lattice.cells())
            .flat_map(|cell| match lattice.cell_type(cell) {
                Cell::Pin(pin) if !pin.is_grid() => (0..pin.regions())
                    .map(|r| if r < 8 { &f } else { &w })
                    .collect::<Vec<_>>(),
                cell => vec![&w; cell.regions()],
            })
            .collect();
        let vacuum = Boundaries {
            left: Boundary::Vacuum,
            right: Boundary::Vacuum,
            top: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
        };
        let tracking = Tracking::with_boundaries(&lattice, 8, 0.05, vacuum).unwrap();
        // 前提: 内側のセルの角をちょうど通るトラックがある（角から直線までの距離が丸めの程度）
        let edges = |sizes: &[f64]| -> Vec<f64> {
            let mut at = -0.5 * sizes.iter().sum::<f64>();
            sizes
                .iter()
                .map(|s| {
                    at += s;
                    at
                })
                .collect()
        };
        // 内側の線だけ（最後の要素は外周の辺）
        let (xs, ys) = (edges(lattice.widths()), edges(lattice.heights()));
        let (xs, ys) = (&xs[..xs.len() - 1], &ys[..ys.len() - 1]);
        let through_corner = tracking.tracks.iter().any(|t| {
            let d = (t.end.0 - t.start.0, t.end.1 - t.start.1);
            let norm = d.0.hypot(d.1);
            xs.iter().any(|&x| {
                ys.iter()
                    .any(|&y| ((x - t.start.0) * d.1 - (y - t.start.1) * d.0).abs() / norm < 1e-12)
            })
        });
        assert!(through_corner, "no track passes through a corner");
        let options = Options::default();
        let s = solve(&tracking, &materials, &options).unwrap();
        let (columns, rows) = (lattice.columns(), lattice.rows());
        // 面の量を、同じ面の同じ群の中性子束の大きさで比べる（正味の流れは 0 に近いことがある）
        let check = |face: usize, image: usize, sign: f64| {
            for g in 0..2 {
                let scale = s.face_fluxes[face][g];
                for (name, a, b) in [
                    ("flux", s.face_fluxes[face][g], s.face_fluxes[image][g]),
                    (
                        "current",
                        s.face_currents[face][g],
                        sign * s.face_currents[image][g],
                    ),
                ] {
                    assert!(
                        (a - b).abs() < options.flux_tolerance * scale,
                        "{name} of face {face} and {image}, group {g}: {a} vs {b}"
                    );
                }
            }
        };
        // 上下の鏡映: 縦の面は行が反転して流れ（x 向き）はそのまま、横の面は線が反転して流れ（y 向き）の符号が変わる
        for row in 0..rows {
            for line in 0..=columns {
                let face = crate::cmfd::vertical_face(columns, row, line);
                check(
                    face,
                    crate::cmfd::vertical_face(columns, rows - 1 - row, line),
                    1.0,
                );
            }
        }
        for line in 0..=rows {
            for column in 0..columns {
                let face = crate::cmfd::horizontal_face(rows, columns, line, column);
                let image = crate::cmfd::horizontal_face(rows, columns, rows - line, column);
                check(face, image, -1.0);
            }
        }
        // 左右の鏡映: 縦の面は線が反転して流れ（x 向き）の符号が変わり、横の面は列が反転して流れはそのまま
        for row in 0..rows {
            for line in 0..=columns {
                let face = crate::cmfd::vertical_face(columns, row, line);
                check(
                    face,
                    crate::cmfd::vertical_face(columns, row, columns - line),
                    -1.0,
                );
            }
        }
        for line in 0..=rows {
            for column in 0..columns {
                let face = crate::cmfd::horizontal_face(rows, columns, line, column);
                let image = crate::cmfd::horizontal_face(rows, columns, line, columns - 1 - column);
                check(face, image, 1.0);
            }
        }
    }

    /// 1/4 対称だが 1/8 対称でない格子と、その領域の物質（燃料 0、水 1、別の燃料 2、真空 3 の番号）。ピンは扇形 8
    /// （燃料、真空のギャップ、水）。`odd` なら 3 × 3 で、真ん中の列は幅 1.0 の 2 × 2 の水のブロック（左右の鏡映の軸が
    /// セルの中を通る長方形）、真ん中の行のピンは別の燃料。そうでなければ 4 × 4 で、ピンと 2 × 2 の水を左右と上下で
    /// 対称に（対角線では非対称に）置き、周りに水の余白を付けたもの（軸がセルの境を通る）
    fn quarter_lattice(odd: bool) -> (Lattice, Vec<usize>) {
        let pin = PinCell::with_sectors(1.26, vec![0.35, 0.4096], 8).unwrap();
        let lattice = if odd {
            Lattice::rectilinear(
                vec![1.26, 1.0, 1.26],
                vec![1.26; 3],
                vec![Cell::Pin(pin), Cell::Block { nx: 2, ny: 2 }],
                vec![vec![0, 1, 0]; 3],
            )
            .unwrap()
        } else {
            let water = PinCell::grid(1.26, 2).unwrap();
            let layout = vec![
                vec![0, 1, 1, 0],
                vec![0, 0, 0, 0],
                vec![0, 0, 0, 0],
                vec![0, 1, 1, 0],
            ];
            Lattice::with_margin(vec![pin, water], layout, 0.04).unwrap()
        };
        let columns = lattice.columns();
        let materials = (0..lattice.cells())
            .flat_map(|cell| match lattice.cell_type(cell) {
                Cell::Pin(pin) if !pin.is_grid() => {
                    let fuel = if odd && cell / columns == 1 { 2 } else { 0 };
                    (0..pin.regions())
                        .map(|r| [fuel, 3, 1][r / 8])
                        .collect::<Vec<_>>()
                }
                cell => vec![1; cell.regions()],
            })
            .collect();
        (lattice, materials)
    }

    fn symmetric(scattering_order: usize, symmetry: Symmetry) -> Options {
        Options {
            scattering_order,
            symmetry,
            ..Options::default()
        }
    }

    /// 全体を解いた解 `full` と、対称性を使って解いた解 `reduced` が、反復の判定の程度で一致すること
    fn assert_same_solution(full: &Solution, reduced: &Solution, options: &Options) {
        assert!(
            (full.k - reduced.k).abs() < options.k_tolerance,
            "k {} vs {}",
            full.k,
            reduced.k
        );
        // 中性子流は、同じ場所・同じ群の中性子束の大きさに対して比べる（鏡面反射の辺や対称線の上の中性子流は丸めの
        // 程度の 0 なので、自分の大きさに対しては比べられない）。`scale` は [位置][群] の尺度
        let close = |a: &[Vec<f64>], b: &[Vec<f64>], scale: &[Vec<f64>], what: &str| {
            for (i, ((x, y), s)) in a.iter().zip(b).zip(scale).enumerate() {
                for (g, ((x, y), s)) in x.iter().zip(y).zip(s).enumerate() {
                    assert!(
                        (x - y).abs() <= options.flux_tolerance * s.abs(),
                        "{what} {i} group {g}: {x} vs {y} (scale {s})"
                    );
                }
            }
        };
        for (r, (a, b)) in full.flux.iter().zip(&reduced.flux).enumerate() {
            for (g, (a, b)) in a.iter().zip(b).enumerate() {
                assert!(
                    relative_change(*a, *b) < options.flux_tolerance,
                    "region {r} group {g}: {a} vs {b}"
                );
            }
        }
        let faces = &full.face_fluxes;
        close(faces, &reduced.face_fluxes, faces, "face flux");
        close(
            &full.face_currents,
            &reduced.face_currents,
            faces,
            "face current",
        );
        close(&full.sides, &reduced.sides, &full.sides, "side flux");
        close(
            &full.currents,
            &reduced.currents,
            &full.sides,
            "side current",
        );
        // 漏れは外周にわたる和なので、辺の中性子束 × 辺の長さの和に対して（辺は左、右、上、下）
        let [width, height] = full.size;
        let boundary: Vec<f64> = (0..full.leakage.len())
            .map(|g| {
                [height, height, width, width]
                    .iter()
                    .zip(&full.sides)
                    .map(|(l, side)| l * side[g])
                    .sum()
            })
            .collect();
        close(
            std::slice::from_ref(&full.leakage),
            std::slice::from_ref(&reduced.leakage),
            std::slice::from_ref(&boundary),
            "leakage",
        );
        // 領域の中性子流は、領域の中性子束に対して（[領域][x, y の順に群]）
        let currents = |s: &Solution| -> Vec<Vec<f64>> {
            s.region_currents.iter().map(|c| c.concat()).collect()
        };
        let fluxes: Vec<Vec<f64>> = full
            .flux
            .iter()
            .map(|f| [&f[..], &f[..]].concat())
            .collect();
        close(
            &currents(full),
            &currents(reduced),
            &fluxes,
            "region current",
        );
        assert_eq!(full.volumes, reduced.volumes);
    }

    #[test]
    fn quarter_symmetry_reproduces_the_full_solution() {
        let (f, w, g, void) = (
            with_p1(fuel()),
            with_p1(moderator()),
            with_p1(two_group(
                [0.55, 1.4],
                [[0.5, 0.03], [0.0, 1.1]],
                [0.01, 0.3],
            )),
            Material::void(2),
        );
        let library = [&f, &w, &g, &void];
        let vacuum = Boundaries {
            left: Boundary::Vacuum,
            right: Boundary::Vacuum,
            top: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
        };
        // 左右だけ真空（1/8 対称ではない）
        let sides = Boundaries {
            left: Boundary::Vacuum,
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        // アルベドの辺（左右と上下で違う値）
        let albedo = Boundaries {
            left: Boundary::Albedo(0.3),
            right: Boundary::Albedo(0.3),
            top: Boundary::Albedo(0.6),
            bottom: Boundary::Albedo(0.6),
        };
        // 左右がアルベドで上下が鏡面反射
        let mixed = Boundaries {
            left: Boundary::Albedo(0.3),
            right: Boundary::Albedo(0.3),
            ..Boundaries::REFLECTIVE
        };
        // 軸がセルの中を通る長方形と、セルの境を通る正方形。鏡面反射、真空、アルベド、等方散乱と P1
        for odd in [true, false] {
            let (lattice, index) = quarter_lattice(odd);
            let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
            for boundaries in [Boundaries::REFLECTIVE, vacuum, sides, albedo, mixed] {
                let tracking = Tracking::with_boundaries(&lattice, 8, 0.05, boundaries).unwrap();
                for order in [0, 1] {
                    let full =
                        solve(&tracking, &materials, &symmetric(order, Symmetry::NONE)).unwrap();
                    let options = symmetric(order, Symmetry::QUARTER);
                    let quarter = solve(&tracking, &materials, &options).unwrap();
                    assert_same_solution(&full, &quarter, &options);
                }
            }
        }
    }

    /// 正方形の 10 通りの対称性（ADR-0112）
    fn all_symmetries() -> Vec<Symmetry> {
        [
            &[][..],
            &[2],
            &[4],
            &[5],
            &[6],
            &[7],
            &[2, 4, 5],
            &[2, 6, 7],
            &[1, 2, 3],
            &[1, 2, 3, 4, 5, 6, 7],
        ]
        .iter()
        .map(|ops| Symmetry::new(ops).unwrap())
        .collect()
    }

    /// 1/8 対称な格子と、その領域の物質（燃料 0、水 1、別の燃料 2、真空 3 の番号）。ピンは扇形 8（燃料、真空のギャップ、
    /// 水）。`odd` なら 3 × 3 で中央のピンは別の燃料（軸と対角線がセルの中を通る）。そうでなければ 4 × 4 で、辺の中央の
    /// 8 つが 2 × 2 の水、周りに水の余白（軸がセルの境を通る）
    fn octant_lattice(odd: bool) -> (Lattice, Vec<usize>) {
        let pin = PinCell::with_sectors(1.26, vec![0.35, 0.4096], 8).unwrap();
        let lattice = if odd {
            Lattice::new(vec![pin], vec![vec![0; 3]; 3]).unwrap()
        } else {
            let water = PinCell::grid(1.26, 2).unwrap();
            let layout = vec![
                vec![0, 1, 1, 0],
                vec![1, 0, 0, 1],
                vec![1, 0, 0, 1],
                vec![0, 1, 1, 0],
            ];
            Lattice::with_margin(vec![pin, water], layout, 0.04).unwrap()
        };
        let materials = (0..lattice.cells())
            .flat_map(|cell| match lattice.cell_type(cell) {
                Cell::Pin(pin) if !pin.is_grid() => {
                    let fuel = if odd && cell == 4 { 2 } else { 0 };
                    (0..pin.regions())
                        .map(|r| [fuel, 3, 1][r / 8])
                        .collect::<Vec<_>>()
                }
                cell => vec![1; cell.regions()],
            })
            .collect();
        (lattice, materials)
    }

    #[test]
    fn every_symmetry_reproduces_the_full_solution() {
        let (f, w, g, void) = (
            with_p1(fuel()),
            with_p1(moderator()),
            with_p1(two_group(
                [0.55, 1.4],
                [[0.5, 0.03], [0.0, 1.1]],
                [0.01, 0.3],
            )),
            Material::void(2),
        );
        let library = [&f, &w, &g, &void];
        let vacuum = Boundaries {
            left: Boundary::Vacuum,
            right: Boundary::Vacuum,
            top: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
        };
        let albedo = Boundaries {
            left: Boundary::Albedo(0.3),
            right: Boundary::Albedo(0.3),
            top: Boundary::Albedo(0.3),
            bottom: Boundary::Albedo(0.3),
        };
        // 1/8 対称な格子では 10 通りの対称性がどれも成り立つ。方位角 8（45° なし）は等方散乱で、10（45° あり）は P1 で。間隔
        // 0.05 cm では、方位角 10 の 45° のトラックが余白の角のセル（0.04 cm 角）の面を通らないので 0.03 cm にする
        for odd in [true, false] {
            let (lattice, index) = octant_lattice(odd);
            let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
            for boundaries in [Boundaries::REFLECTIVE, vacuum, albedo] {
                for (azimuths, order) in [(8, 0), (10, 1)] {
                    let tracking =
                        Tracking::with_boundaries(&lattice, azimuths, 0.03, boundaries).unwrap();
                    let full =
                        solve(&tracking, &materials, &symmetric(order, Symmetry::NONE)).unwrap();
                    for symmetry in all_symmetries() {
                        let options = symmetric(order, symmetry);
                        let reduced = solve(&tracking, &materials, &options)
                            .unwrap_or_else(|e| panic!("{}: {e}", symmetry.name()));
                        assert_same_solution(&full, &reduced, &options);
                    }
                }
            }
        }
    }

    #[test]
    fn axis_swapping_symmetries_need_a_square_lattice() {
        // 1/4 対称な長方形（真ん中の列が幅 1.0）では、軸を入れ替えない対称性は解けて全体の解と一致し、軸を入れ替える
        // 操作を含む対称性はエラー
        let (lattice, index) = quarter_lattice(true);
        let (f, w, void) = (fuel(), moderator(), Material::void(2));
        let library = [&f, &w, &f, &void];
        let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
        let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
        let full = solve(&tracking, &materials, &symmetric(0, Symmetry::NONE)).unwrap();
        for symmetry in all_symmetries() {
            let options = symmetric(0, symmetry);
            let result = solve(&tracking, &materials, &options);
            let swaps = symmetry
                .operations()
                .iter()
                .any(|&i| [1, 3, 6, 7].contains(&i));
            if swaps {
                assert!(
                    matches!(result, Err(MocError::Symmetry(_))),
                    "{}: {result:?}",
                    symmetry.name()
                );
            } else {
                assert_same_solution(&full, &result.unwrap(), &options);
            }
        }
    }

    #[test]
    fn octant_result_does_not_depend_on_the_thread_count() {
        let (lattice, index) = octant_lattice(true);
        let (f, w, void) = (with_p1(fuel()), with_p1(moderator()), Material::void(2));
        let library = [&f, &w, &f, &void];
        let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
        let tracking = Tracking::new(&lattice, 10, 0.05).unwrap();
        let options = symmetric(1, Symmetry::OCTANT);
        let solve_on = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| solve(&tracking, &materials, &options).unwrap())
        };
        let serial = solve_on(1);
        for threads in [2, 5] {
            assert_eq!(solve_on(threads), serial, "{threads} threads");
        }
    }

    #[test]
    fn quarter_symmetry_solves_only_representatives() {
        // 軸がセルの境を通る格子では、どの領域とセルも 4 つずつ移り合う
        let (lattice, index) = quarter_lattice(false);
        let (f, w, void) = (fuel(), moderator(), Material::void(2));
        let library = [&f, &w, &f, &void];
        let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
        let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
        let mesh = CoarseMesh::new(&tracking);
        let reduced =
            Representatives::symmetric(&tracking, &materials, &mesh, &[(0, 0)], Symmetry::QUARTER)
                .unwrap();
        assert_eq!(4 * reduced.regions.representatives.len(), lattice.regions());
        assert_eq!(4 * reduced.cells.representatives.len(), lattice.cells());
    }

    #[test]
    fn quarter_result_does_not_depend_on_the_thread_count() {
        // 代表のトラックの組と、代表の領域に集める和の順は、スレッドの分け方によらない
        let (lattice, index) = quarter_lattice(true);
        let (f, w, void) = (with_p1(fuel()), with_p1(moderator()), Material::void(2));
        let library = [&f, &w, &f, &void];
        let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
        let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
        let options = symmetric(1, Symmetry::QUARTER);
        let solve_on = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| solve(&tracking, &materials, &options).unwrap())
        };
        let serial = solve_on(1);
        for threads in [2, 5] {
            assert_eq!(solve_on(threads), serial, "{threads} threads");
        }
    }

    #[test]
    fn quarter_symmetry_rejects_asymmetric_materials() {
        // 左上のピンの1つの扇形だけ物質が違う
        let (lattice, mut index) = quarter_lattice(false);
        let sector = index.iter().position(|&i| i == 0).unwrap();
        index[sector] = 2;
        let (f, w, void) = (fuel(), moderator(), Material::void(2));
        let g = Material {
            total: vec![0.55, 1.4],
            ..fuel()
        };
        let library = [&f, &w, &g, &void];
        let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
        let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
        assert!(solve(&tracking, &materials, &symmetric(0, Symmetry::NONE)).is_ok());
        let result = solve(&tracking, &materials, &symmetric(0, Symmetry::QUARTER));
        assert!(
            matches!(result, Err(MocError::AsymmetricMaterial { .. })),
            "{result:?}"
        );
        // 操作は番号ではなく名前で示す
        let message = result.unwrap_err().to_string();
        assert!(
            message.contains("under the rotation by 180 degrees")
                || message.contains("under the mirror across"),
            "{message}"
        );
    }

    #[test]
    fn quarter_symmetry_needs_one_material_for_images() {
        // 値が同じでも別の物質（自己遮蔽で別々に解いたピンにあたる）は、移り合う領域に使えない。値の比較で決めると、
        // 別々に解いたピンが丸めの違いだけで対称かどうかが変わる
        let (lattice, index) = quarter_lattice(false);
        let (f, w, void) = (fuel(), moderator(), Material::void(2));
        let copy = fuel();
        let library = [&f, &w, &f, &void];
        let mut materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
        let last = index.iter().rposition(|&i| i == 0).unwrap();
        materials[last] = &copy;
        let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
        assert!(solve(&tracking, &materials, &symmetric(0, Symmetry::NONE)).is_ok());
        let result = solve(&tracking, &materials, &symmetric(0, Symmetry::QUARTER));
        assert!(
            matches!(result, Err(MocError::AsymmetricMaterial { .. })),
            "{result:?}"
        );
        // 事実（別の物質を指している）を示し、自己遮蔽で別々に解かれた可能性はヒントにとどめる
        let message = result.unwrap_err().to_string();
        assert!(
            message.contains("refer to different materials") && message.contains("ADR-0114"),
            "{message}"
        );
    }

    #[test]
    fn quarter_symmetry_rejects_an_asymmetric_boundary() {
        // 右だけ真空: 上下の鏡映しか残らない
        let (lattice, index) = quarter_lattice(false);
        let (f, w, void) = (fuel(), moderator(), Material::void(2));
        let library = [&f, &w, &f, &void];
        let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&lattice, 8, 0.05, boundaries).unwrap();
        let result = solve(&tracking, &materials, &symmetric(0, Symmetry::QUARTER));
        assert!(matches!(result, Err(MocError::Symmetry(_))), "{result:?}");
    }
}
