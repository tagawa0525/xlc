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
use crate::symmetry::{Detected, Symmetry, SymmetryError};
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
    /// [群] の外周からの漏れ（真空の辺と、アルベドの辺で戻らない 1 − α の分。中性子束と同じ規格化）
    pub leakage: Vec<f64>,
    /// [外周の辺（左、右、上、下）][群] の辺の上の平均スカラー中性子束（中性子束と同じ規格化）。
    /// 集合体の不連続因子（ADF）に使う
    pub sides: Vec<Vec<f64>>,
    /// [外周の辺（左、右、上、下）][群] の辺の上で平均した、外向きの正味の中性子流（中性子束と同じ規格化）。
    /// CMFD のために数える粗メッシュの面の正味の中性子流を辺にわたって足し、辺の長さで割ったもの。
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

/// 収束した反復の状態（ADR-0195）。同じトラッキングと対称操作で解く次の問題の反復を、ここから始める
/// （`solve_symmetric_from`）。端の角度中性子束の全体を持つので、複製しない（`Clone` を持たない）
#[derive(Debug)]
pub struct Start {
    /// 解いたトラッキングの領域とトラックの数
    sizes: [usize; 2],
    /// 代表の領域と代表のトラックの番号（解いたトラッキングと対称操作の目印）
    representatives: [Vec<usize>; 2],
    k: f64,
    /// [代表の領域][群] のスカラー中性子束
    flux: Vec<Vec<f64>>,
    /// [代表の領域][モーメント (1 番目から)][群] の高次のモーメント
    higher: Vec<Vec<f64>>,
    /// [代表のトラック][向き][極角][群] のトラックの端の角度中性子束
    outgoing: Vec<f64>,
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
    #[error("the symmetry was verified for another tracking")]
    ForeignSymmetry,
    #[error(
        "region {region} and its image {image} under the {operation} refer to different \
         materials (pass the same material for regions that the operation maps onto each \
         other; self-shielding does so for the operations it verifies, ADR-0112)"
    )]
    AsymmetricMaterial {
        operation: &'static str,
        region: usize,
        image: usize,
    },
    #[error(
        "the start was saved from another tracking, symmetry or scattering order (its regions, \
         tracks or sizes differ)"
    )]
    ForeignStart,
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
pub fn solve(
    tracking: &Tracking,
    materials: &[&Material],
    options: &Options,
) -> Result<Solution, MocError> {
    solve_with(tracking, materials, options, None, None).map(|(solution, _)| solution)
}

/// `solve` と同じ問題を、確かめた対称操作 `symmetry`（`symmetry::detect` か `symmetry::require` の結果。同じトラッキングの
/// もの）で移り合うトラックと領域の代表だけで解く（ADR-0112）。移り合う領域が同じ物質（同じ参照）でなければエラー
pub fn solve_symmetric(
    tracking: &Tracking,
    materials: &[&Material],
    options: &Options,
    symmetry: &Detected,
) -> Result<Solution, MocError> {
    solve_with(tracking, materials, options, Some(symmetry), None).map(|(solution, _)| solution)
}

/// `solve_symmetric` と同じく解き、反復を `start`（同じトラッキング、対称操作、散乱の次数で解いた前の問題の収束した状態）から
/// 始める（ADR-0195）。`None` なら平らな初期値から。収束した状態を解と一緒に返す。`start` の領域とトラックの数、代表の領域と
/// トラックの番号、配列の大きさのどれかが今の問題と違えばエラー
pub fn solve_symmetric_from(
    tracking: &Tracking,
    materials: &[&Material],
    options: &Options,
    symmetry: &Detected,
    start: Option<&Start>,
) -> Result<(Solution, Start), MocError> {
    solve_with(tracking, materials, options, Some(symmetry), start)
}

#[tracing::instrument(name = "moc", skip_all, fields(regions = materials.len(), tracks = tracking.tracks.len()))]
fn solve_with(
    tracking: &Tracking,
    materials: &[&Material],
    options: &Options,
    symmetry: Option<&Detected>,
    start: Option<&Start>,
) -> Result<(Solution, Start), MocError> {
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
    // 反復の前の用意（入力の確認、体積、粗メッシュ、代表、スイープの表、初期値）
    let setup = tracing::trace_span!("setup").entered();
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
    let representatives = match symmetry {
        Some(verified) if verified.symmetry() != Symmetry::NONE => {
            // 別の格子やトラッキングで確かめた像は使えない
            if !verified.fits(tracking) {
                return Err(MocError::ForeignSymmetry);
            }
            Representatives::by_images(tracking, materials, &mesh, &pairs, verified.images())?
        }
        _ => Representatives::all(tracking, &mesh, pairs.len()),
    };
    let sweeper = Sweeper::new(tracking, materials, ng, &mesh, &pairs, representatives);
    // 状態（中性子束、モーメント、平衡値）は、領域の軌道の代表 `kept[i]` の分だけ持つ（対称性を使わなければすべての領域）
    let reduced = &sweeper.representatives;
    let kept = &reduced.regions.representatives;
    let states = kept.len();
    let material = |i: usize| materials[kept[i]];
    let volume = |i: usize| volumes[kept[i]];
    // [代表] の軌道の領域の体積の和（生成率はすべての領域の和）
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

    let width = (pairs.len() - 1) * ng;
    let sizes = [regions, tracking.tracks.len()];
    // 前の解の中性子束、[代表][モーメント (1 番目から)][群] の高次のモーメント、k。なければ平らな中性子束と k = 1
    let (mut flux, mut higher, mut k) = match start {
        Some(s) => {
            // 代表が同じなら、配列の長さは代表の数で揃う。幅は群とモーメントの数で確かめる
            if s.sizes != sizes
                || s.representatives[0] != *kept
                || s.representatives[1] != reduced.tracks
                || s.flux.iter().any(|f| f.len() != ng)
                || s.higher.iter().any(|h| h.len() != width)
            {
                return Err(MocError::ForeignStart);
            }
            (s.flux.clone(), s.higher.clone(), s.k)
        }
        None => (
            vec![vec![1.0; ng]; states],
            vec![vec![0.0; width]; states],
            1.0,
        ),
    };
    // 前の解も、今の物質の生成率で 1 に規格化し直す
    let initial = production(&flux);
    if initial <= 0.0 {
        return Err(MocError::NoFission);
    }
    for f in flux.iter_mut().chain(&mut higher).flatten() {
        *f /= initial;
    }
    // スイープの終わりに出ていく角度中性子束（`Sweeper::sweep` を参照）。前の解がなければ端の領域の中性子束の
    // 等方な分布とする（0 から始めると、鏡面反射の辺でも外向きの流れがあるように見え、CMFD が
    // 漏れの大きな問題として解いてしまう）
    let mut outgoing = match start {
        Some(s) => s.outgoing.iter().map(|b| b / initial).collect(),
        None => sweeper.isotropic_outgoing(&flux),
    };
    let mut next_outgoing = outgoing.clone();
    let mut tallies = sweeper.tallies(regions);
    let mut k_history = Remaining::default();
    let mut flux_history = Remaining::default();
    drop(setup);

    for iteration in 1..=options.max_iterations {
        // 線源を Σt で割ったもの（角度中性子束の平衡値）の [モーメント][群] の展開係数
        let equilibrium: Vec<Vec<f64>> = tracing::trace_span!("equilibrium").in_scope(|| {
            (0..states)
                .into_par_iter()
                .map(|i| equilibrium_moments(material(i), &flux[i], &higher[i], k, &pairs))
                .collect()
        });

        // 段の時間は trace の区間で測る（main.rs の FmtSpan::CLOSE で、閉じるときに time.busy を出す）
        tracing::trace_span!("sweep").in_scope(|| {
            sweeper.sweep(&equilibrium, &outgoing, &mut next_outgoing, &mut tallies)
        })?;
        std::mem::swap(&mut outgoing, &mut next_outgoing);
        let sums = tracing::trace_span!("tallies").in_scope(|| {
            let mut sums = sum_tallies(&sweeper, &tallies);
            sweeper
                .representatives
                .symmetrize(&mut sums, &sweeper.void, ng);
            sums
        });
        let Sums {
            delta,
            void_flux,
            mut leak,
            mut current,
            face_flux,
            void_current,
            region_current,
        } = sums;

        // スイープの和から、代表ごとの中性子束、高次のモーメント、中性子流を求め、生成率 1 に規格化する。代表ごとの値は
        // ほかの代表によらないので、代表ごとに並列に求める（どの値も同じ式なので、結果は分け方によらない）
        let update = tracing::trace_span!("update").entered();
        let mut new_flux: Vec<Vec<f64>> = (0..states)
            .into_par_iter()
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
            .into_par_iter()
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
            .into_par_iter()
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
        // べき乗法: 生成率の比で k を更新し、中性子束と境界の値を生成率 1 に規格化する。生成率の和は代表の順に
        // 足す（並列に足すと、足す順が分け方で変わる）
        let rate = production(&new_flux);
        let mut new_k = k * rate;
        new_flux
            .par_iter_mut()
            .chain(new_higher.par_iter_mut())
            .chain(new_currents.par_iter_mut())
            .for_each(|values| values.iter_mut().for_each(|v| *v /= rate));
        for l in leak.iter_mut().chain(&mut current) {
            *l /= rate;
        }
        outgoing
            .par_chunks_mut(sweeper.width)
            .for_each(|values| values.iter_mut().for_each(|b| *b /= rate));
        drop(update);
        // CMFD: 粗メッシュの固有値問題の解で、k と、セルごと・群ごとの中性子束の比を更新する。
        // 粗メッシュの解は生成率 1 に規格化されているので、更新後も生成率は 1 のまま
        if let Acceleration::Cmfd(cmfd) = &options.acceleration {
            let cells = &reduced.cells;
            let coarse = tracing::trace_span!("homogenize").in_scope(|| {
                Coarse::homogenize(
                    &mesh,
                    cells,
                    materials,
                    &volumes,
                    &new_flux,
                    &reduced.regions.index,
                    &current,
                )
            })?;
            // MOC がまだ収束から遠いうちは、CMFD を細かく解いても次のスイープで崩れる。スイープによる変化（CMFD の誤差を
            // 含まない）に合わせて精度を緩める
            let tuned = inexact_cmfd(
                cmfd,
                options,
                (new_k - k).abs(),
                largest_change(&new_flux, &flux),
            );
            // [代表のセル][群] の比
            let (coarse_k, ratio) =
                tracing::trace_span!("coarse").in_scope(|| coarse.solve(new_k, &tuned))?;
            new_k = coarse_k;
            let _scale = tracing::trace_span!("scale").entered();
            // 高次のモーメントも同じ比で（角度分布の形は変えない）
            new_flux
                .par_iter_mut()
                .zip(new_higher.par_iter_mut())
                .zip(new_currents.par_iter_mut())
                .enumerate()
                .for_each(|(i, ((f, h), c))| {
                    let cell = cells.index[mesh.cell(kept[i])];
                    for (j, value) in f
                        .iter_mut()
                        .chain(h.iter_mut())
                        .chain(c.iter_mut())
                        .enumerate()
                    {
                        *value *= ratio[cell * ng + j % ng];
                    }
                });
            // トラックの端の角度中性子束も、端のあるセルの比で更新する（領域の中性子束だけを
            // 更新すると、次のスイープで入ってくる流れと食い違う）
            sweeper.scale_outgoing(&mut outgoing, |cell, g| ratio[cells.index[cell] * ng + g]);
        }
        let converge = tracing::trace_span!("converge").entered();
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
        drop(converge);
        if converged {
            let _finish = tracing::trace_span!("finish").entered();
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
            let solution = Solution {
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
            };
            let start = Start {
                sizes,
                representatives: [kept.clone(), reduced.tracks.clone()],
                k,
                flux,
                higher,
                outgoing,
            };
            return Ok((solution, start));
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
        add(start.faces());
        let mut cell = mesh.cell(track.segments[0].region);
        for segment in &track.segments {
            let next = mesh.cell(segment.region);
            if next != cell {
                add(mesh.path(cell, next)?.faces());
                cell = next;
            }
        }
        add(end.faces());
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
    /// 正味の中性子流を数える粗メッシュ
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
    /// 非等方散乱のスイープの作業用の配列（再利用する。中身は方位角ごとに作り直すので、どれを使っても結果は同じ）
    pool: Mutex<Vec<Directional>>,
    /// スイープするトラック（対称性を使わなければすべて）
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

    /// 出ていく角度中性子束に、トラックの出口の粗メッシュのセルと群ごとの係数 `factor` を掛ける。出口が外周の 2 つの
    /// 面の境にあれば、出ていく経路と同じく両側のセルの係数を 1/2 ずつ平均する（片側のセルの係数だけを掛けると、
    /// 両側の面に半分ずつ数えた流れと食い違う）
    fn scale_outgoing(&self, outgoing: &mut [f64], factor: impl Fn(usize, usize) -> f64 + Sync) {
        outgoing
            .par_chunks_mut(self.width)
            .enumerate()
            .for_each(|(j, out)| {
                let track = self.representatives.tracks[j / 2];
                let cells = self.mesh.end_cells(track, j % 2 == 0);
                for (i, o) in out.iter_mut().enumerate() {
                    let g = i % self.ng;
                    *o *= cells
                        .iter()
                        .map(|&(cell, w)| w * factor(cell, g))
                        .sum::<f64>();
                }
            });
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
                    // 代表のトラックは、移り合うトラックの数だけ重く数える（ADR-0112）
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
                        let entered = enter.scaled(direction);
                        add_current(&mut tally.current, entered.faces(), psi, &weights);
                        self.add_face_flux(
                            &mut tally.face_flux,
                            entered.faces(),
                            psi,
                            side_weights,
                        );
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
                            // 等方なら [群] をすべての極角で使い、非等方なら [極角][群]
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
                        let exited = exit.scaled(direction);
                        add_current(&mut tally.current, exited.faces(), psi, &weights);
                        self.add_face_flux(&mut tally.face_flux, exited.faces(), psi, side_weights);
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
    /// [粗メッシュの面][群] の正味の中性子流（面の正の向きが正）
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
    /// [粗メッシュの面][群] の正味の中性子流
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

/// 領域と群のうちで最も大きい中性子束の相対変化。領域ごとに並列に求める（最大は足す順によらず、NaN はどの順でも
/// `f64::max` が捨てるので、結果は分け方によらない）
fn largest_change(new: &[Vec<f64>], old: &[Vec<f64>]) -> f64 {
    new.par_iter()
        .zip(old)
        .map(|(new, old)| {
            new.iter()
                .zip(old)
                .map(|(new, old)| relative_change(*new, *old))
                .fold(0.0, f64::max)
        })
        .reduce(|| 0.0, f64::max)
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

/// トラックの1つの向きの、方向ごとの重み 4π (w_a/2) ω_p sinθ_p δ_a
fn direction_weight(tracking: &Tracking, azimuth: usize, p: usize) -> f64 {
    let a = &tracking.azimuths[azimuth];
    4.0 * PI * 0.5 * a.weight * POLAR_WEIGHT[p] * POLAR_SIN[p] * a.spacing
}

#[cfg(test)]
mod tests;
