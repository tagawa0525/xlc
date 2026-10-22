//! ピンごとの超微細群（UFG）の減速計算。
//!
//! 同心円のリングのセル（白色境界）で、全核種の弾性散乱（重心系で等方、標的は静止）による減速を、超微細な格子の上で
//! 高いエネルギーから順に解く。レサジー u = ln(E_high/E) と F_i(u) = E φ_i(E) で、リング i の散乱による発生密度は
//!
//! q_i(u) = Σ_k N_k/(1 − α_k) ∫_{u−Δk}^{u} σs,k(u') F_i(u') e^{−(u−u')} du'、α = ((A − 1)/(A + 1))²、Δ = ln(1/α)
//!
//! となり、衝突確率法の飛程 W（[`Cylinder::tracks`]）で F_j V_j = Σ_i q_i W_ij とつながる。σs F は格子点の間で線形とし、
//! 窓の積分は窓をずらしながら更新する（窓が範囲全体にわたる H でも、格子点の数に比例する手間で済む）。格子点での
//! 散乱（最後の区間の、その点の値の項）は未知の F_i を含むので、点ごとにリングの数の連立方程式を解く。
//!
//! 範囲の上（E_high より上）は、中性子束をどのリングでも 1/E（F = 1）とし、散乱断面積はライブラリの無限希釈の
//! 群定数とする（範囲の上限はアクチニドの分離共鳴の上端で、その上の遮蔽は σ0 の表で扱う）。群の上端より上に
//! 中性子はない。範囲の下限（熱中性子散乱の打ち切り）より下は計算しない。

use std::borrow::Cow;
use std::ops::Range;
use std::time::{Duration, Instant};

use rayon::prelude::*;

use xlc_format::library::{Library, LookupError, Pointwise, Shielded, TemperatureData};

/// 1つの窓（散乱で届くレサジーの幅 Δ）に入れる格子点の最少の数。xs のライブラリの重みの減速計算の同名の定数とは
/// 独立に選べる、この計算の格子の細かさの設定である。同じ減速の式を同じ格子の作り方で解くので、xs で収束を確かめた
/// 値（xs の `finer_grids_do_not_change_the_flux`）をそのまま使っている。この計算での収束はまだ確かめていない
const WINDOW_POINTS: usize = 16;

use crate::cp::{CpError, Cylinder};
use crate::input::MaterialInput;
use crate::track_table::{DEFAULT_DENSITY, Layout, MAX_DEVIATION, TableError, TrackTable};
use xlc_numeric::lu::solve;

/// ポイントワイズの断面積の列の数: [全, 弾性散乱, 核分裂, νσf, 捕獲]
const COLUMNS: usize = 5;
type Row = [f64; COLUMNS];

/// セルの UFG の範囲の上限の頭打ち (eV)。U-238 の分離共鳴の上端で、これより上は U-238 の非弾性散乱（45 keV から）と
/// 核分裂の線源が効き、弾性散乱だけの減速計算が成り立たない（ADR-0069）
pub const MAX_TOP: f64 = 2.0e4;

/// 飛程をまとめて並列に求める格子点の数
const CHUNK: usize = 4096;

/// 減速計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum UfgError {
    #[error("{given} ring materials for {rings} rings")]
    Rings { given: usize, rings: usize },
    #[error("range {range:?} must be two group bounds inside the pointwise range {pointwise:?}")]
    Range {
        range: [f64; 2],
        pointwise: Option<[f64; 2]>,
    },
    #[error(
        "the fuel has resolved resonances but the library has no pointwise data \
         (made without any nuclide with resolved resonances)"
    )]
    NoPointwise,
    #[error("{name}: atomic weight ratio {awr} must be positive")]
    Mass { name: String, awr: f64 },
    #[error("{name}: pointwise data are not defined at {energy} eV")]
    Undefined { name: String, energy: f64 },
    #[error("the slowing-down equations are singular at {0} eV")]
    Singular(f64),
    #[error("ring {ring}: non-finite group constants in group {group}")]
    NonFinite { ring: usize, group: usize },
    #[error(transparent)]
    Lookup(#[from] LookupError),
    #[error(transparent)]
    Cp(#[from] CpError),
    #[error(transparent)]
    Table(#[from] TableError),
    #[error(
        "ring {0}: the total cross section is 0 at some lethargy point, so the track table cannot cover it (use `track-table = false`)"
    )]
    ZeroTotal(usize),
}

/// セルの減速計算の結果。
#[derive(Debug, Clone, PartialEq)]
pub struct CellSpectrum {
    /// 減速計算した群（ライブラリの群の番号）
    pub groups: Range<usize>,
    /// [リング][群] の中性子束 ∫φ dE（相対値。群は `groups` の中の順）
    pub flux: Vec<Vec<f64>>,
    /// [リング][核種（物質の組成の順）] の、リングの中性子束で縮約した群定数 (barn)。真空のリングは空
    pub nuclides: Vec<Vec<Shielded>>,
}

/// セルの UFG の範囲 [下限, 上限] (eV)。下限はライブラリのポイントワイズの断面積の下限（熱中性子散乱の打ち切り以上の群の境界）、
/// 上限は燃料（核分裂する核種を含む物質）の核種の分離共鳴の上端のうち最も高いもの（[`MAX_TOP`] で頭打ち）以下で
/// 最も高い群の境界。`rings` はリングごとの物質（真空は `None`）。燃料に分離共鳴を持つ核種がなければ `None`
/// （減速計算で遮蔽する共鳴がない）。
pub fn cell_range(
    library: &Library,
    rings: &[Option<&MaterialInput>],
) -> Result<Option<[f64; 2]>, UfgError> {
    let mut top: Option<f64> = None;
    for m in rings.iter().flatten() {
        let data = m
            .composition
            .iter()
            .map(|(name, _)| library.lookup(name, m.temperature))
            .collect::<Result<Vec<_>, _>>()?;
        let fuel = data.iter().any(|(_, t)| {
            t.shielded
                .first()
                .is_some_and(|s| s.nu_fission.iter().any(|&v| v > 0.0))
        });
        if !fuel {
            continue;
        }
        for r in data.iter().filter_map(|(d, _)| d.resolved.as_ref()) {
            top = Some(top.map_or(r.top, |t| t.max(r.top)));
        }
    }
    top.map(|top| range_below(library, top)).transpose()
}

/// 物質 `material` の UFG の範囲 [下限, 上限] (eV)。[`cell_range`] と同じ決め方で、上限を燃料でなく `material` の核種
/// （核分裂するかによらない）の分離共鳴の上端から決める（燃料でない棒の遮蔽に使う）。分離共鳴を持つ核種がなければ `None`。
pub fn material_range(
    library: &Library,
    material: &MaterialInput,
) -> Result<Option<[f64; 2]>, UfgError> {
    let mut top: Option<f64> = None;
    for (name, _) in &material.composition {
        if let Some(r) = &library.lookup(name, material.temperature)?.0.resolved {
            top = Some(top.map_or(r.top, |t| t.max(r.top)));
        }
    }
    top.map(|top| range_below(library, top)).transpose()
}

/// 分離共鳴の上端 `top` (eV) からの UFG の範囲（[`cell_range`]）。
fn range_below(library: &Library, top: f64) -> Result<[f64; 2], UfgError> {
    let pointwise = library.pointwise_range.ok_or(UfgError::NoPointwise)?;
    let top = top.min(MAX_TOP).min(pointwise[1]);
    let low = pointwise[0];
    let high = library
        .groups
        .bounds()
        .iter()
        .copied()
        .find(|&b| b <= top)
        .filter(|&b| b > low)
        .ok_or(UfgError::Range {
            range: [low, top],
            pointwise: Some(pointwise),
        })?;
    Ok([low, high])
}

/// 1つの核種（名前と温度の組）の減速計算のデータ
struct Data<'a> {
    name: &'a str,
    /// 1/(1 − α)
    gain: f64,
    /// 窓の幅 Δ = ln(1/α)（α = 0 なら無限）
    width: f64,
    /// その温度の群定数とポイントワイズの断面積（S(α, β) の核種は内挿したもの）
    table: Cow<'a, TemperatureData>,
}

impl Data<'_> {
    fn pointwise(&self) -> &Pointwise {
        &self.table.pointwise
    }

    /// 範囲の上で使う無限希釈の弾性散乱の群定数
    fn elastic(&self) -> &[f64] {
        &self.table.shielded[0].elastic
    }
}

/// 減速計算の格子点ごとの飛程の求め方。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Collision {
    /// 格子点ごとに衝突確率法で求める（表の検証の基準）
    Exact,
    /// 物質ごとの ln Σ の格子（10 倍あたり `density` 点）で求めた値から内挿する（[`TrackTable`]。ADR-0156）
    Table { density: usize },
}

impl Default for Collision {
    /// ADR-0156 の既定
    fn default() -> Self {
        Self::Table {
            density: DEFAULT_DENSITY,
        }
    }
}

/// セル `cell` の、範囲 `range`（[`cell_range`]）の減速計算。`rings` はリングごとの物質（真空は `None`）。
/// 飛程は既定の求め方（[`Collision::default`]）で求める。
pub fn slowing_down(
    library: &Library,
    cell: &Cylinder,
    rings: &[Option<&MaterialInput>],
    range: [f64; 2],
) -> Result<CellSpectrum, UfgError> {
    slowing_down_with(library, cell, rings, range, Collision::default())
}

/// [`slowing_down`] の、飛程の求め方を `collision` で選ぶもの。
#[tracing::instrument(name = "ufg", level = "debug", skip_all)]
pub fn slowing_down_with(
    library: &Library,
    cell: &Cylinder,
    rings: &[Option<&MaterialInput>],
    range: [f64; 2],
    collision: Collision,
) -> Result<CellSpectrum, UfgError> {
    let volumes = cell.volumes();
    let n = volumes.len();
    if rings.len() != n {
        return Err(UfgError::Rings {
            given: rings.len(),
            rings: n,
        });
    }
    let bounds = library.groups.bounds();
    let [low, high] = range;
    let pointwise = library.pointwise_range.ok_or(UfgError::NoPointwise)?;
    let bad = || UfgError::Range {
        range,
        pointwise: Some(pointwise),
    };
    let top = bounds.iter().position(|&b| b == high).ok_or_else(bad)?;
    let bottom = bounds.iter().position(|&b| b == low).ok_or_else(bad)?;
    if !(low < high && low >= pointwise[0] && high <= pointwise[1]) {
        return Err(bad());
    }

    // 核種（名前と温度の組）ごとのデータと、リングごとの（核種の番号、個数密度）
    let mut data: Vec<Data> = Vec::new();
    let mut keys: Vec<(&str, f64)> = Vec::new();
    let mut entries: Vec<Vec<(usize, f64)>> = Vec::with_capacity(n);
    for m in rings {
        let mut ring = Vec::new();
        for (name, density) in m.iter().flat_map(|m| &m.composition) {
            let temperature = m.map_or(f64::NAN, |m| m.temperature);
            let found = keys
                .iter()
                .position(|&(k, t)| k == name && t == temperature);
            let k = match found {
                Some(k) => k,
                None => {
                    let (d, t) = library.lookup(name, temperature)?;
                    if !(d.awr > 0.0 && d.awr.is_finite()) {
                        return Err(UfgError::Mass {
                            name: name.clone(),
                            awr: d.awr,
                        });
                    }
                    let alpha = ((d.awr - 1.0) / (d.awr + 1.0)).powi(2);
                    data.push(Data {
                        name: &d.name,
                        gain: 1.0 / (1.0 - alpha),
                        width: -alpha.ln(),
                        table: t,
                    });
                    keys.push((name, temperature));
                    data.len() - 1
                }
            };
            ring.push((k, *density));
        }
        entries.push(ring);
    }

    // 格子（エネルギーの降順）: ポイントワイズの断面積の格子と群の境界に、どの核種の窓にも WINDOW_POINTS 点が入るよう点を足す
    let group_bounds = &bounds[top..=bottom];
    let mut energies: Vec<f64> = data
        .iter()
        .flat_map(|d| d.pointwise().energy.iter().copied())
        .filter(|&e| e > low && e < high)
        .chain(group_bounds.iter().copied())
        .collect();
    energies.sort_by(|a, b| b.total_cmp(a));
    energies.dedup();
    let step = data
        .iter()
        .map(|d| d.width / WINDOW_POINTS as f64)
        .fold(f64::INFINITY, f64::min);
    let (u, energies) = lethargy_grid(&energies, high, step, |e| group_bounds.contains(&e));
    // ポイントワイズの断面積の不連続点（同じエネルギーの2点）では、高い側と低い側の極限の2つの格子点を同じレサジーに置く。
    // 範囲の上端では低い側、下端では高い側の極限だけを使う
    let mut jumps: Vec<f64> = data
        .iter()
        .flat_map(|d| d.pointwise().energy.windows(2))
        .filter(|w| w[0] == w[1] && w[0] > low && w[0] < high)
        .map(|w| w[0])
        .collect();
    jumps.sort_by(f64::total_cmp);
    jumps.dedup();
    // （エネルギー、低い側の極限か）
    let mut points: Vec<(f64, bool)> = Vec::with_capacity(energies.len() + jumps.len());
    let mut lethargies: Vec<f64> = Vec::with_capacity(energies.len() + jumps.len());
    for (&ue, &e) in u.iter().zip(&energies) {
        let jump = jumps.binary_search_by(|x| x.total_cmp(&e)).is_ok();
        points.push((e, e == high));
        lethargies.push(ue);
        if jump {
            points.push((e, true));
            lethargies.push(ue);
        }
    }
    let u = lethargies;
    tracing::debug!(points = u.len(), rings = entries.len(), "UFG slowing down");

    let ng = bottom - top;
    let above = Above { bounds, top, high };
    let mut sweep = Sweep {
        data: &data,
        entries: &entries,
        volumes: &volumes,
        u: &u,
        elastic: vec![Vec::with_capacity(u.len()); data.len()],
        flux: vec![Vec::with_capacity(u.len()); n],
        above,
        lower: vec![0; data.len()],
        // 範囲の上端での窓の積分（範囲の上の中性子束は 1/E）
        windows: entries
            .iter()
            .map(|ring| {
                ring.iter()
                    .map(|&(k, _)| above.integral(&data[k], -data[k].width, 0.0, 0.0))
                    .collect()
            })
            .collect(),
        previous: Vec::new(),
        group: 0,
        group_flux: vec![vec![0.0; ng]; n],
        rates: entries
            .iter()
            .map(|ring| vec![vec![[0.0; COLUMNS]; ng]; ring.len()])
            .collect(),
        known: Vec::new(),
        own: entries.iter().map(|ring| vec![0.0; ring.len()]).collect(),
        slabs: vec![(0.0, Vec::new()); data.len()],
    };

    // 表で内挿するなら、軸ごとの全断面積の範囲を格子点の全体で先に求めて表を作る（表で衝突確率が減らないセルは作らない）
    let table = match collision {
        Collision::Exact => None,
        Collision::Table { density } => track_table(cell, &data, &entries, &points, density)?,
    };

    let mut cursors: Vec<usize> = data.iter().map(|d| d.pointwise().energy.len()).collect();
    // 段ごとの壁時計の時間（ログの内訳。衝突確率は格子点の間で並列に求める）。ピンの間も並列のとき、衝突確率を待つ
    // スレッドがほかのピンの仕事を盗むと、その時間も衝突確率に入るので、段の割合は 1 スレッドで測る
    let mut times = [Duration::ZERO; 3];
    for start in (0..u.len()).step_by(CHUNK) {
        let clock = Instant::now();
        let chunk = start..(start + CHUNK).min(u.len());
        // [点][核種] のポイントワイズの断面積 [全, 弾性散乱, 核分裂, νσf, 捕獲]
        let mut values: Vec<Vec<Row>> = Vec::with_capacity(chunk.len());
        for &(e, below) in &points[chunk.clone()] {
            let row = data
                .iter()
                .zip(cursors.iter_mut())
                .map(|(d, c)| evaluate(d, c, e, below))
                .collect::<Result<Vec<_>, _>>()?;
            values.push(row);
        }
        let evaluated = Instant::now();
        let tracks = values
            .par_iter()
            .map(|row| match &table {
                None => {
                    let totals: Vec<f64> = entries.iter().map(|ring| total(ring, row)).collect();
                    cell.tracks(&totals).map_err(UfgError::from)
                }
                Some((table, first)) => {
                    let sigma: Vec<f64> = first.iter().map(|&i| total(&entries[i], row)).collect();
                    match table.tracks(&sigma)? {
                        Some(tracks) => Ok(tracks),
                        // 軸を共有する物質の差が展開の上限を超える格子点（ADR-0182）
                        None => {
                            let totals: Vec<f64> =
                                entries.iter().map(|ring| total(ring, row)).collect();
                            cell.tracks(&totals).map_err(UfgError::from)
                        }
                    }
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        let collided = Instant::now();
        for ((m, row), w) in chunk.zip(values).zip(&tracks) {
            let (e, below) = points[m];
            sweep.point(m, row, w, e)?;
            // 群の境界が不連続点なら、高い側の格子点で次の群に移る（低い側との間は幅 0）
            if sweep.group < ng && !below && e == group_bounds[sweep.group + 1] {
                sweep.group += 1;
            }
        }
        times[0] += evaluated - clock;
        times[1] += collided - evaluated;
        times[2] += collided.elapsed();
    }
    let [evaluate, collision, slowing] = times.map(|t| t.as_secs_f64());
    tracing::debug!(evaluate, collision, slowing, "UFG time (s)");

    // 群定数
    let mut nuclides = Vec::with_capacity(n);
    for (i, ring) in sweep.rates.iter().enumerate() {
        let flux = &sweep.group_flux[i];
        let mut list = Vec::with_capacity(ring.len());
        for rate in ring {
            let column =
                |x: usize| -> Vec<f64> { rate.iter().zip(flux).map(|(r, f)| r[x] / f).collect() };
            let shielded = Shielded {
                total: column(0),
                elastic: column(1),
                fission: column(2),
                nu_fission: column(3),
                capture: column(4),
            };
            let finite = |v: &Vec<f64>| v.iter().position(|x| !x.is_finite());
            let bad = [
                &shielded.total,
                &shielded.elastic,
                &shielded.fission,
                &shielded.nu_fission,
                &shielded.capture,
            ]
            .into_iter()
            .find_map(finite);
            if let Some(g) = bad {
                return Err(UfgError::NonFinite {
                    ring: i,
                    group: top + g,
                });
            }
            list.push(shielded);
        }
        nuclides.push(list);
    }
    Ok(CellSpectrum {
        groups: top..bottom,
        flux: sweep.group_flux,
        nuclides,
    })
}

/// リングの全断面積 Σ_k N_k σ_k（`row` は [核種] のポイントワイズの断面積）
fn total(ring: &[(usize, f64)], row: &[Row]) -> f64 {
    ring.iter().map(|&(k, d)| d * row[k][0]).sum()
}

/// セルの飛程の表と、[物質] の最初のリングの番号。同じ組成（核種と温度、個数密度）のリングは全断面積が同じなので
/// 1 つの物質にまとめ、核種の並び（核種と温度）が同じ物質は 1 つの軸を共有する（ADR-0182）。範囲は格子点 `points` の
/// 全体での軸の値（中心）の最小と最大。表で求める衝突確率の数が格子点の数以上なら、表を作らずに `None`（格子点ごとに求める）
fn track_table(
    cell: &Cylinder,
    data: &[Data],
    entries: &[Vec<(usize, f64)>],
    points: &[(f64, bool)],
    density: usize,
) -> Result<Option<(TrackTable, Vec<usize>)>, UfgError> {
    let mut first: Vec<usize> = Vec::new();
    let classes: Vec<Option<usize>> = entries
        .iter()
        .enumerate()
        .map(|(i, ring)| {
            if ring.is_empty() {
                return None;
            }
            Some(match first.iter().position(|&r| entries[r] == *ring) {
                Some(m) => m,
                None => {
                    first.push(i);
                    first.len() - 1
                }
            })
        })
        .collect();
    // 物質ごとの軸: 核種の並びが同じ物質の最初のものの軸
    let nuclides = |i: usize| entries[i].iter().map(|&(k, _)| k);
    let mut leaders: Vec<usize> = Vec::new();
    let axes: Vec<usize> = first
        .iter()
        .map(
            |&i| match leaders.iter().position(|&r| nuclides(r).eq(nuclides(i))) {
                Some(a) => a,
                None => {
                    leaders.push(i);
                    leaders.len() - 1
                }
            },
        )
        .collect();
    let layout = Layout::new(cell, classes, axes)?;
    let mut ranges = vec![[f64::INFINITY, f64::NEG_INFINITY]; leaders.len()];
    // 展開の上限を超えて格子点ごとに求める格子点の数
    let mut outside = 0usize;
    let mut cursors: Vec<usize> = data.iter().map(|d| d.pointwise().energy.len()).collect();
    let mut sigma = vec![0.0; first.len()];
    for &(e, below) in points {
        let row = data
            .iter()
            .zip(cursors.iter_mut())
            .map(|(d, c)| evaluate(d, c, e, below))
            .collect::<Result<Vec<_>, _>>()?;
        for (s, &i) in sigma.iter_mut().zip(&first) {
            let t = total(&entries[i], &row);
            // 格子点ごとに求めるとき（[`Cylinder::tracks`]）と同じく、有限で負でない全断面積だけを受ける
            // （min と max は NaN を飛ばすので、ここで止める）
            if !(t.is_finite() && t >= 0.0) {
                return Err(CpError::Total(i).into());
            }
            // 表は ln Σ の格子なので、全断面積が 0 になる物質は表にできない（格子点ごとに求めれば解ける）
            if t == 0.0 {
                return Err(UfgError::ZeroTotal(i));
            }
            *s = t;
        }
        if !widen_ranges(&layout, &sigma, &mut ranges) {
            outside += 1;
        }
    }
    // どの格子点も展開の外なら、表を読む格子点がない
    if outside == points.len() {
        tracing::info!(
            points = points.len(),
            "UFG tracks point by point: every point is beyond the ring expansion"
        );
        return Ok(None);
    }
    // 表で求める衝突確率（格子点の組ごとの表の値と微分、展開の外の格子点）が格子点ごとより多いなら、表は速くならず、メモリも
    // 大きいので作らない。組の数が usize に収まらない表も同じ
    let evaluations = match TrackTable::count(&ranges, density) {
        Ok(nodes) => nodes
            .checked_mul(layout.evaluations())
            .and_then(|n| n.checked_add(outside)),
        Err(TableError::Size) => None,
        Err(e) => return Err(e.into()),
    };
    if evaluations.is_none_or(|n| n >= points.len()) {
        tracing::info!(
            points = points.len(),
            table = evaluations,
            "UFG tracks point by point: the track table needs as many collision probabilities"
        );
        return Ok(None);
    }
    // 軸を共有する物質のあるセル（リングごとの燃料）は、展開の外の格子点の数を残す（ADR-0182 の確認に使う）
    if layout.evaluations() > 1 {
        tracing::info!(
            points = points.len(),
            outside,
            evaluations,
            "UFG track table with ring expansion"
        );
    }
    Ok(Some((
        TrackTable::new(cell, layout, &ranges, density)?,
        first,
    )))
}

/// 物質ごとの全断面積 `sigma` の格子点が展開の中（[`MAX_DEVIATION`] 以下）なら、軸の範囲 `ranges` をその中心まで広げて
/// true。展開の外の格子点は表を読まない（格子点ごとに求める）ので、範囲を広げずに false
fn widen_ranges(layout: &Layout, sigma: &[f64], ranges: &mut [[f64; 2]]) -> bool {
    let (centers, deviation) = layout.centers(sigma);
    if deviation > MAX_DEVIATION {
        return false;
    }
    for (range, c) in ranges.iter_mut().zip(centers) {
        range[0] = range[0].min(c);
        range[1] = range[1].max(c);
    }
    true
}

/// 降順のエネルギー `energies` のレサジーの格子と、その点のエネルギー。間隔が `step` より広いところは点を足す。
/// 丸めで同じレサジーになる点は1つにまとめ、`keep` を満たす（群の境界の）エネルギーを残す。
fn lethargy_grid(
    energies: &[f64],
    high: f64,
    step: f64,
    keep: impl Fn(f64) -> bool,
) -> (Vec<f64>, Vec<f64>) {
    let mut u: Vec<f64> = Vec::with_capacity(energies.len());
    let mut points: Vec<f64> = Vec::with_capacity(energies.len());
    for &e in energies {
        let ue = (high / e).ln();
        if let Some(&last) = u.last() {
            if ue <= last {
                if keep(e)
                    && let Some(p) = points.last_mut()
                {
                    *p = e;
                }
                continue;
            }
            if step.is_finite() {
                let gap = ue - last;
                let extra = (gap / step).ceil() as usize;
                for k in 1..extra {
                    let fill = last + gap * k as f64 / extra as f64;
                    u.push(fill);
                    points.push(high * (-fill).exp());
                }
            }
        }
        u.push(ue);
        points.push(e);
    }
    (u, points)
}

/// 核種 `d` のポイントワイズの断面積 [全, 弾性散乱, 核分裂, νσf, 捕獲] をエネルギー `e` で lin-lin に補間する。`cursor` はエネルギーの
/// 降順に引くための位置。不連続点では、`below` なら低い側、そうでなければ高い側の極限の値になる。
fn evaluate(d: &Data, cursor: &mut usize, e: f64, below: bool) -> Result<Row, UfgError> {
    let x = &d.pointwise().energy;
    let defined =
        x.len() >= 2 && x.first().is_some_and(|&a| a <= e) && x.last().is_some_and(|&b| e <= b);
    if !defined {
        return Err(UfgError::Undefined {
            name: d.name.to_string(),
            energy: e,
        });
    }
    *cursor = (*cursor).min(x.len() - 1);
    while *cursor > 0 && x[*cursor] > e {
        *cursor -= 1;
    }
    if below {
        // 同じエネルギーの点のうち最初のもの（低い側の極限）
        while *cursor > 0 && x[*cursor - 1] == e {
            *cursor -= 1;
        }
    }
    let (j, t) = if below && *cursor > 0 && x[*cursor] == e {
        (*cursor - 1, 1.0)
    } else {
        let j = (*cursor).min(x.len() - 2);
        let width = x[j + 1] - x[j];
        (j, if width > 0.0 { (e - x[j]) / width } else { 0.0 })
    };
    let p = d.pointwise();
    let lerp = |y: &[f64]| y[j] + t * (y[j + 1] - y[j]);
    Ok([
        lerp(&p.total),
        lerp(&p.elastic),
        lerp(&p.fission),
        lerp(&p.nu_fission),
        lerp(&p.capture),
    ])
}

/// 範囲の上（u < 0）の群
#[derive(Clone, Copy)]
struct Above<'a> {
    bounds: &'a [f64],
    /// 範囲の上端の群の境界の番号
    top: usize,
    high: f64,
}

impl Above<'_> {
    /// 範囲の上（F = 1）での ∫_a^b σs(u') e^{−(U − u')} du'（a < b ≤ 0）。σs は無限希釈の群定数で、群の上端より
    /// 上は 0。
    fn integral(&self, d: &Data, a: f64, b: f64, big_u: f64) -> f64 {
        (0..self.top)
            .map(|g| {
                let ga = (self.high / self.bounds[g]).ln();
                let gb = (self.high / self.bounds[g + 1]).ln();
                let (x, y) = (a.max(ga), b.min(gb));
                if x < y {
                    d.elastic()[g] * ((-(big_u - y)).exp() - (-(big_u - x)).exp())
                } else {
                    0.0
                }
            })
            .sum()
    }
}

/// 区間 [a, b] で線形な g の ∫_a^b g(u') e^{−(U − u')} du' を、両端の値 g(a), g(b) の係数 (ca, cb) で表す。
fn segment(a: f64, b: f64, big_u: f64) -> (f64, f64) {
    let s = b - a;
    if s == 0.0 {
        // 不連続点の2つの格子点の間
        return (0.0, 0.0);
    }
    // 1 − e^{−s} と、(s − (1 − e^{−s})) / s
    let em = -(-s).exp_m1();
    let tail = (s - em) / s;
    let decay = (-(big_u - b)).exp();
    (decay * (em - tail), decay * tail)
}

/// 高いエネルギーから順に解く減速計算の状態
struct Sweep<'a> {
    data: &'a [Data<'a>],
    entries: &'a [Vec<(usize, f64)>],
    volumes: &'a [f64],
    u: &'a [f64],
    above: Above<'a>,
    /// [核種][格子点] の弾性散乱
    elastic: Vec<Vec<f64>>,
    /// [リング][格子点] の F
    flux: Vec<Vec<f64>>,
    /// 核種ごとの、窓の下端を含む区間の始まりの格子点
    lower: Vec<usize>,
    /// [リング][項] の窓の積分 ∫ σs F e^{−(u−u')} du'
    windows: Vec<Vec<f64>>,
    /// 1つ前の格子点のポイントワイズの断面積
    previous: Vec<Row>,
    /// 今の区間の群（`groups` の中の番号）
    group: usize,
    group_flux: Vec<Vec<f64>>,
    /// [リング][項][群] の反応率 ∫ F σ du [全, 弾性散乱, 核分裂, νσf, 捕獲]
    rates: Vec<Vec<Vec<Row>>>,
    /// 格子点ごとに再利用する作業用の配列: [リング][項] の窓の積分の既知の部分
    known: Vec<Vec<f64>>,
    /// 格子点ごとに再利用する作業用の配列: [リング][項] のその点の F_i の係数
    own: Vec<Vec<f64>>,
    /// [核種] の、窓の下端の区間の積分の範囲の上の部分と、格子点の値の重み
    slabs: Vec<(f64, Vec<(usize, f64)>)>,
}

impl Sweep<'_> {
    /// 格子点 m の中性子束を求め、窓の積分と群の積分を進める。`values` はその点のポイントワイズの断面積、`tracks` は飛程。
    fn point(
        &mut self,
        m: usize,
        values: Vec<Row>,
        tracks: &[Vec<f64>],
        energy: f64,
    ) -> Result<(), UfgError> {
        let n = self.volumes.len();
        for (column, v) in self.elastic.iter_mut().zip(&values) {
            column.push(v[1]);
        }
        // 窓の積分のうち既知の部分 [リング][項] と、その点の F_i の係数
        let known = &mut self.known;
        known.clone_from(&self.windows);
        let own = &mut self.own;
        if m > 0 {
            let (u0, u1) = (self.u[m - 1], self.u[m]);
            let decay = (-(u1 - u0)).exp();
            let (ca, cb) = segment(u0, u1, u1);
            // 核種ごとの窓の下端の区間 [u0 − Δ, u1 − Δ] の積分の、格子点の値の重み（範囲の上の部分は値で）
            let above = self.above;
            for ((d, lower), (outside, weights)) in self
                .data
                .iter()
                .zip(self.lower.iter_mut())
                .zip(self.slabs.iter_mut())
            {
                let (a, b) = (u0 - d.width, u1 - d.width);
                *outside = if a < 0.0 {
                    above.integral(d, a, b.min(0.0), u1)
                } else {
                    0.0
                };
                slab(self.u, a, b, u1, lower, weights);
            }
            for (i, ring) in self.entries.iter().enumerate() {
                let f = &self.flux[i];
                own[i].clear();
                own[i].extend(ring.iter().map(|&(k, _)| cb * self.elastic[k][m]));
                for (e, &(k, _)) in ring.iter().enumerate() {
                    let (outside, weights) = &self.slabs[k];
                    let leaving: f64 = outside
                        + weights
                            .iter()
                            .map(|&(p, w)| w * self.elastic[k][p] * f[p])
                            .sum::<f64>();
                    known[i][e] =
                        decay * known[i][e] + ca * self.elastic[k][m - 1] * f[m - 1] - leaving;
                }
            }
        } else {
            for row in own.iter_mut() {
                row.fill(0.0);
            }
        }
        // F_j V_j − Σ_i c_i F_i W_ij = Σ_i k_i W_ij
        let source = |i: usize, terms: &[f64]| -> f64 {
            self.entries[i]
                .iter()
                .zip(terms)
                .map(|(&(k, d), t)| d * self.data[k].gain * t)
                .sum()
        };
        let k: Vec<f64> = (0..n).map(|i| source(i, &known[i])).collect();
        let c: Vec<f64> = (0..n).map(|i| source(i, &own[i])).collect();
        // tracks[j][i] は、リング i の線源からリング j の中を飛ぶ長さ W_ij。表から内挿した飛程は対称でないので、
        // この向きで読む（格子点ごとに求めた W は対称で、向きによらない）
        let matrix: Vec<f64> = tracks
            .iter()
            .zip(self.volumes)
            .enumerate()
            .flat_map(|(j, (row, v))| {
                row.iter()
                    .zip(&c)
                    .enumerate()
                    .map(move |(i, (w, ci))| if i == j { v - w * ci } else { -w * ci })
            })
            .collect();
        let rhs: Vec<f64> = tracks
            .iter()
            .map(|row| row.iter().zip(&k).map(|(w, ki)| w * ki).sum())
            .collect();
        let f = solve(matrix, rhs).ok_or(UfgError::Singular(energy))?;
        for (i, ring) in self.entries.iter().enumerate() {
            for (e, _) in ring.iter().enumerate() {
                self.windows[i][e] = known[i][e] + own[i][e] * f[i];
            }
            self.flux[i].push(f[i]);
        }
        // 群の積分（レサジーの台形則）
        if m > 0 {
            let du = self.u[m] - self.u[m - 1];
            let g = self.group;
            for (i, ring) in self.entries.iter().enumerate() {
                let (f0, f1) = (self.flux[i][m - 1], f[i]);
                self.group_flux[i][g] += 0.5 * du * (f0 + f1);
                for (e, &(k, _)) in ring.iter().enumerate() {
                    let (v0, v1) = (self.previous[k], values[k]);
                    for x in 0..COLUMNS {
                        self.rates[i][e][g][x] += 0.5 * du * (f0 * v0[x] + f1 * v1[x]);
                    }
                }
            }
        }
        self.previous = values;
        Ok(())
    }
}

/// 窓の下端が [a, b] を動くときに窓から出る ∫ σs F e^{−(U − u')} du'（σs F は格子点の間で線形）のうち、範囲の中
/// （u ≥ 0）の部分を、格子点の σs F にかける重みで `weights` に置く。`lower` は下端を含む区間の始まりの格子点（進めていく）。
fn slab(u: &[f64], a: f64, b: f64, big_u: f64, lower: &mut usize, weights: &mut Vec<(usize, f64)>) {
    weights.clear();
    let mut x = a.max(0.0);
    if b <= x {
        return;
    }
    while *lower + 1 < u.len() && u[*lower + 1] <= x {
        *lower += 1;
    }
    let mut q = *lower;
    while x < b && q + 1 < u.len() {
        let (ua, ub) = (u[q], u[q + 1]);
        let y = b.min(ub);
        if y > x {
            let (c0, c1) = segment(x, y, big_u);
            // g(x), g(y) を格子点 q, q + 1 の値の線形補間で表す
            let (tx, ty) = ((x - ua) / (ub - ua), (y - ua) / (ub - ua));
            weights.push((q, c0 * (1.0 - tx) + c1 * (1.0 - ty)));
            weights.push((q + 1, c0 * tx + c1 * ty));
        }
        x = y;
        q += 1;
    }
}

#[cfg(test)]
mod tests;
