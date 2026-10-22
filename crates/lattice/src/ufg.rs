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
use std::collections::HashSet;
use std::ops::Range;
use std::time::{Duration, Instant};

use rayon::prelude::*;

use xlc_format::library::{Library, LookupError, Shielded, TemperatureData};

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
    #[error(
        "ring {0}: the background must give, for each ring, one entry per nuclide of its material (none for a void ring) with every group"
    )]
    Background(usize),
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
    /// エネルギーの格子点の数（ポイントワイズの断面積の不連続点は、高い側と低い側の2点を数える）
    pub points: usize,
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

/// 1つの核種（名前と温度の組）か、1つのリングの背景の減速計算のデータ
struct Data<'a> {
    name: &'a str,
    /// 1/(1 − α)
    gain: f64,
    /// 窓の幅 Δ = ln(1/α)（α = 0 なら無限）
    width: f64,
    source: Source<'a>,
}

/// 減速計算の断面積の出どころ
enum Source<'a> {
    /// その温度の群定数とポイントワイズの断面積（S(α, β) の核種は内挿したもの）
    Pointwise(Box<Cow<'a, TemperatureData>>),
    /// リングの背景の核種をまとめた、群の中で一定の巨視的な断面積 (1/cm)。[ライブラリの群] の [全, 弾性散乱, 0, 0, 0]
    /// （背景の核種の群定数は減速計算の入力のまま出すので、反応率の列は使わない）
    Groups { bounds: &'a [f64], rows: Vec<Row> },
}

impl Data<'_> {
    /// ポイントワイズの断面積のエネルギーの格子（昇順。背景は空）
    fn energies(&self) -> &[f64] {
        match &self.source {
            Source::Pointwise(t) => &t.pointwise.energy,
            Source::Groups { .. } => &[],
        }
    }

    /// 範囲の上での弾性散乱の群 `g` の値（核種は無限希釈の群定数、背景は巨視的な断面積）
    fn elastic(&self, g: usize) -> f64 {
        match &self.source {
            Source::Pointwise(t) => t.shielded[0].elastic[g],
            Source::Groups { rows, .. } => rows[g][1],
        }
    }

    /// 降順に引くときの `evaluate` の位置の初めの値
    fn first_cursor(&self) -> usize {
        match &self.source {
            Source::Pointwise(t) => t.pointwise.energy.len(),
            Source::Groups { .. } => 0,
        }
    }
}

/// 減速計算の格子点ごとの飛程の求め方。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Collision {
    /// 格子点ごとに衝突確率法で求める（表の検証の基準）
    Exact,
    /// 軸ごとの ln Σ の格子（10 倍あたり `density` 点）で求めた値から内挿する（[`TrackTable`]。ADR-0156）。核種の並びが同じ
    /// 物質は 1 つの軸を共有して 1 次の展開で内挿し、展開の外の格子点と、表で衝突確率が減らないセルは格子点ごとに求める
    /// （ADR-0182）
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
/// すべての核種をポイントワイズの断面積で解き、飛程は既定の求め方（[`Collision::default`]）で求める。
pub fn slowing_down(
    library: &Library,
    cell: &Cylinder,
    rings: &[Option<&MaterialInput>],
    range: [f64; 2],
) -> Result<CellSpectrum, UfgError> {
    slowing_down_with(
        library,
        cell,
        rings,
        &no_background(rings),
        range,
        Collision::default(),
        None,
    )
}

/// [`slowing_down_with`] の、背景の核種のない（すべての核種をポイントワイズの断面積で解く）`background`
pub fn no_background(rings: &[Option<&MaterialInput>]) -> Vec<Vec<Option<Shielded>>> {
    rings
        .iter()
        .map(|m| vec![None; m.map_or(0, |m| m.composition.len())])
        .collect()
}

/// [`slowing_down`] の、背景の核種と飛程の求め方を選ぶもの。`background` は [リング][核種（物質の組成の順）] の、群の中で
/// 一定とする群定数（ライブラリのすべての群）で、`None` の核種はポイントワイズの断面積で解く。リングの背景の核種は、減速計算では
/// 1 つの項にまとめ（[`lump`]）、群定数は核種のポイントワイズの断面積をリングの中性子束で縮約して返す（[`collapse`]。ADR-0185）。
#[tracing::instrument(name = "ufg", level = "debug", skip_all)]
pub fn slowing_down_with(
    library: &Library,
    cell: &Cylinder,
    rings: &[Option<&MaterialInput>],
    background: &[Vec<Option<Shielded>>],
    range: [f64; 2],
    collision: Collision,
    thinning: Option<f64>,
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
    if background.len() != n {
        return Err(UfgError::Background(background.len().min(n)));
    }
    let ng_all = library.groups.len();
    for (i, (m, b)) in rings.iter().zip(background).enumerate() {
        let shaped = b.len() == m.map_or(0, |m| m.composition.len())
            && b.iter().flatten().all(|s| {
                [&s.total, &s.elastic, &s.fission, &s.nu_fission, &s.capture]
                    .iter()
                    .all(|v| v.len() == ng_all)
            });
        if !shaped {
            return Err(UfgError::Background(i));
        }
    }

    let started = Instant::now();
    // 核種（名前と温度の組）とリングの背景ごとのデータと、リングごとの（データの番号、個数密度）。背景の項の個数密度は 1
    // （断面積が巨視的）
    let mut data: Vec<Data> = Vec::new();
    let mut keys: Vec<(&str, f64)> = Vec::new();
    let mut entries: Vec<Vec<(usize, f64)>> = Vec::with_capacity(n);
    // [リング][組成の核種] の群定数の出どころ
    let mut slots: Vec<Vec<Slot>> = Vec::with_capacity(n);
    // 背景の項を持つリングと、その項のデータの番号
    let mut lumps: Vec<(usize, usize)> = Vec::new();
    for (i, m) in rings.iter().enumerate() {
        let mut ring = Vec::new();
        let mut slot = Vec::new();
        let mut members: Vec<(&str, f64, &Shielded)> = Vec::new();
        for ((name, density), constant) in m.iter().flat_map(|m| &m.composition).zip(&background[i])
        {
            let temperature = m.map_or(f64::NAN, |m| m.temperature);
            if let Some(s) = constant {
                members.push((name, *density, s));
                slot.push(Slot::Background { name, temperature });
                continue;
            }
            let found = keys
                .iter()
                .position(|&(k, t)| k == name && t == temperature);
            let k = match found {
                Some(k) => k,
                None => {
                    let (d, t) = library.lookup(name, temperature)?;
                    let alpha = scattering_alpha(name, d.awr)?;
                    data.push(Data {
                        name: &d.name,
                        gain: 1.0 / (1.0 - alpha),
                        width: -alpha.ln(),
                        source: Source::Pointwise(Box::new(t)),
                    });
                    keys.push((name, temperature));
                    data.len() - 1
                }
            };
            slot.push(Slot::Entry(ring.len()));
            ring.push((k, *density));
        }
        if let (Some(m), false) = (m, members.is_empty()) {
            // 組成と背景が同じリングは、同じ項を共有する（同じ全断面積のリングを飛程の表で 1 つの物質にまとめるため）
            let same = |r: usize| {
                rings[r].is_some_and(|o| {
                    o.temperature == m.temperature && o.composition == m.composition
                }) && background[r] == background[i]
            };
            let k = match lumps.iter().find(|&&(r, _)| same(r)) {
                Some(&(_, k)) => k,
                None => {
                    data.push(lump(library, &members, bounds, top..bottom)?);
                    data.len() - 1
                }
            };
            lumps.push((i, k));
            ring.push((k, 1.0));
        }
        entries.push(ring);
        slots.push(slot);
    }

    let looked_up = Instant::now();
    // 格子（エネルギーの降順）: ポイントワイズの断面積の格子と群の境界に、どの核種の窓にも WINDOW_POINTS 点が入るよう点を足す
    let group_bounds = &bounds[top..=bottom];
    let inside = |e: f64| e > low && e < high;
    let mut energies: Vec<f64> = data
        .iter()
        .flat_map(|d| d.energies().iter().copied())
        .filter(|&e| inside(e))
        .chain(group_bounds.iter().copied())
        .collect();
    energies.sort_by(|a, b| b.total_cmp(a));
    energies.dedup();
    // 核種の格子の和（群の境界を含む）の点の数。ここから窓の点を足したものが格子
    let union = energies.len();
    if tracing::enabled!(tracing::Level::DEBUG) {
        log_nuclide_points(&data, group_bounds, inside);
    }
    // ポイントワイズの断面積の不連続点（同じエネルギーの2点）では、高い側と低い側の極限の2つの格子点を同じレサジーに置く。
    // 範囲の上端では低い側、下端では高い側の極限だけを使う
    // 背景の断面積は群の境界で跳ぶので、範囲の中の群の境界も不連続点にする
    let mut jumps: Vec<f64> = data
        .iter()
        .flat_map(|d| d.energies().windows(2))
        .filter(|w| w[0] == w[1] && w[0] > low && w[0] < high)
        .map(|w| w[0])
        .chain(
            group_bounds
                .iter()
                .copied()
                .filter(|&e| !lumps.is_empty() && inside(e)),
        )
        .collect();
    jumps.sort_by(f64::total_cmp);
    jumps.dedup();
    // 間引くなら、群の境界と不連続点を残し、リングの物質ごとの巨視的な全断面積と散乱断面積を許容差で再現できる点だけを残す（ADR-0185）
    let energies = match thinning {
        None => energies,
        Some(tolerance) => thin(&energies, &data, &entries, high, tolerance, |e| {
            group_bounds.contains(&e) || jumps.binary_search_by(|x| x.total_cmp(&e)).is_ok()
        })?,
    };
    let step = data
        .iter()
        .map(|d| d.width / WINDOW_POINTS as f64)
        .fold(f64::INFINITY, f64::min);
    let (u, energies) = lethargy_grid(&energies, high, step, |e| group_bounds.contains(&e));
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
    let grid_done = Instant::now();
    let table = match collision {
        Collision::Exact => None,
        Collision::Table { density } => track_table(cell, &data, &entries, &points, density)?,
    };
    let table_done = Instant::now();

    let mut cursors: Vec<usize> = data.iter().map(Data::first_cursor).collect();
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
        let exact = |row: &[Row]| {
            let totals: Vec<f64> = entries.iter().map(|ring| total(ring, row)).collect();
            cell.tracks(&totals).map_err(UfgError::from)
        };
        let tracks = values
            .par_iter()
            .map(|row| match &table {
                None => exact(row),
                Some((table, first)) => {
                    let sigma: Vec<f64> = first.iter().map(|&i| total(&entries[i], row)).collect();
                    match table.tracks(&sigma)? {
                        Some(tracks) => Ok(tracks),
                        // 軸を共有する物質の差が展開の上限を超える格子点（ADR-0182）
                        None => exact(row),
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
    // 点の数と、段ごとの時間 (s): lookup = 核種のデータの引き、grid = 格子点とスイープの用意、table = 飛程の表の作成
    tracing::debug!(
        points = u.len(),
        union,
        jumps = jumps.len(),
        nuclides = data.len(),
        rings = n,
        lookup = (looked_up - started).as_secs_f64(),
        grid = (grid_done - looked_up).as_secs_f64(),
        table = (table_done - grid_done).as_secs_f64(),
        evaluate,
        collision,
        slowing,
        "UFG slowing down"
    );

    // 群定数。背景の核種は、核種（名前と温度の組）ごとにポイントワイズの断面積を 1 度だけ引く
    let mut tables: Vec<((&str, f64), Cow<TemperatureData>)> = Vec::new();
    let mut nuclides = Vec::with_capacity(n);
    for (i, ring) in slots.iter().enumerate() {
        let flux = &sweep.group_flux[i];
        let mut list = Vec::with_capacity(ring.len());
        for slot in ring {
            let shielded = match *slot {
                Slot::Entry(e) if thinning.is_some() => {
                    // 間引いた格子では核種自身の共鳴の点が落ちうるので、核種の格子で縮約する
                    let Source::Pointwise(table) = &data[entries[i][e].0].source else {
                        return Err(UfgError::Background(i));
                    };
                    collapse(
                        data[entries[i][e].0].name,
                        table,
                        &u,
                        &points,
                        &sweep.flux[i],
                        high,
                        group_bounds,
                        flux,
                    )?
                }
                Slot::Entry(e) => {
                    let rate = &sweep.rates[i][e];
                    let column = |x: usize| -> Vec<f64> {
                        rate.iter().zip(flux).map(|(r, f)| r[x] / f).collect()
                    };
                    Shielded {
                        total: column(0),
                        elastic: column(1),
                        fission: column(2),
                        nu_fission: column(3),
                        capture: column(4),
                    }
                }
                Slot::Background { name, temperature } => {
                    let found = tables.iter().position(|(k, _)| *k == (name, temperature));
                    let k = match found {
                        Some(k) => k,
                        None => {
                            tables
                                .push(((name, temperature), library.lookup(name, temperature)?.1));
                            tables.len() - 1
                        }
                    };
                    collapse(
                        name,
                        &tables[k].1,
                        &u,
                        &points,
                        &sweep.flux[i],
                        high,
                        group_bounds,
                        flux,
                    )?
                }
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
        points: u.len(),
    })
}

/// 降順のエネルギー `energies` を、リングの物質（`entries` の同じ並びは 1 つ）ごとの巨視的な全断面積と弾性散乱断面積が、残した点の
/// 間のレサジーについての線形な補間で、どの点でも相対差 `tolerance` 以内に再現できる点に間引く。`keep` を満たす点（群の境界と
/// 不連続点）と両端は残す。残す点は、補間の傾きの許される範囲を点ごとに狭めて決める（swinging door）
fn thin(
    energies: &[f64],
    data: &[Data],
    entries: &[Vec<(usize, f64)>],
    high: f64,
    tolerance: f64,
    keep: impl Fn(f64) -> bool,
) -> Result<Vec<f64>, UfgError> {
    let _ = (data, entries, high, tolerance, keep);
    Ok(energies.to_vec())
}

/// リングの組成の核種の群定数の出どころ
enum Slot<'a> {
    /// ポイントワイズで解いた、リングの項の番号
    Entry(usize),
    /// 背景の核種（名前と温度）
    Background { name: &'a str, temperature: f64 },
}

/// 背景の核種の群定数: ポイントワイズの断面積を、リングの中性子束 `flux`（減速計算の格子点の値。点の間で u について線形）で
/// 縮約する。ポイントワイズで解いた核種と同じく σF をレサジーの台形則で積分し、点は減速計算の格子点と核種自身の格子の点の和
/// （格子点の間の核種の点では中性子束を内挿する）にする。中性子束は背景の核種を群の中で一定の断面積で入れて解いたもので、
/// 核種自身の共鳴によるへこみを持たない（ほかの核種の共鳴の影は入る。ADR-0185）。`group_flux` は群ごとの ∫F du
#[allow(clippy::too_many_arguments)]
fn collapse(
    name: &str,
    table: &TemperatureData,
    u: &[f64],
    points: &[(f64, bool)],
    flux: &[f64],
    high: f64,
    group_bounds: &[f64],
    group_flux: &[f64],
) -> Result<Shielded, UfgError> {
    let p = &table.pointwise;
    let d = Data {
        name,
        gain: 1.0,
        width: 0.0,
        source: Source::Pointwise(Box::new(Cow::Borrowed(table))),
    };
    let ng = group_flux.len();
    let mut rates = vec![[0.0; COLUMNS]; ng];
    let mut cursor = d.first_cursor();
    // 核種の格子の、まだ使っていない最も高いエネルギーの点の次（昇順の配列を降順に進める）
    let mut next = p.energy.len();
    let mut group = 0;
    // 1 つ前の点の (u, F, 断面積)
    let mut previous: Option<(f64, f64, Row)> = None;
    let mut add = |group: usize, a: &(f64, f64, Row), b: &(f64, f64, Row)| {
        let du = b.0 - a.0;
        for (r, (va, vb)) in rates[group].iter_mut().zip(a.2.iter().zip(&b.2)) {
            *r += 0.5 * du * (a.1 * va + b.1 * vb);
        }
    };
    for (m, &(e, below)) in points.iter().enumerate() {
        let here = (u[m], flux[m], evaluate(&d, &mut cursor, e, below)?);
        if let Some(mut last) = previous {
            let upper = points[m - 1].0;
            while next > 0 && p.energy[next - 1] >= upper {
                next -= 1;
            }
            while next > 0 && p.energy[next - 1] > e {
                next -= 1;
                let x = (high / p.energy[next]).ln();
                let f = last.1 + (here.1 - last.1) * (x - last.0) / (here.0 - last.0);
                let node = (
                    x,
                    f,
                    [
                        p.total[next],
                        p.elastic[next],
                        p.fission[next],
                        p.nu_fission[next],
                        p.capture[next],
                    ],
                );
                add(group, &last, &node);
                last = node;
            }
            add(group, &last, &here);
        }
        previous = Some(here);
        // 群の境界では、高い側の格子点で次の群に移る（減速計算と同じ）
        if group < ng && !below && e == group_bounds[group + 1] {
            group += 1;
        }
    }
    let column = |x: usize| -> Vec<f64> {
        rates
            .iter()
            .zip(group_flux)
            .map(|(r, f)| r[x] / f)
            .collect()
    };
    Ok(Shielded {
        total: column(0),
        elastic: column(1),
        fission: column(2),
        nu_fission: column(3),
        capture: column(4),
    })
}

/// 弾性散乱の α = ((A − 1)/(A + 1))²
fn scattering_alpha(name: &str, awr: f64) -> Result<f64, UfgError> {
    if !(awr > 0.0 && awr.is_finite()) {
        return Err(UfgError::Mass {
            name: name.to_string(),
            awr,
        });
    }
    Ok(((awr - 1.0) / (awr + 1.0)).powi(2))
}

/// 散乱 1 回あたりのレサジーの増加の平均 ξ = 1 + α ln α / (1 − α)（α = 0 で 1）
fn mean_lethargy_gain(alpha: f64) -> f64 {
    if alpha == 0.0 {
        1.0
    } else {
        1.0 + alpha * alpha.ln() / (1.0 - alpha)
    }
}

/// リングの背景の核種 `members`（名前、個数密度、群定数）を 1 つの項にまとめる。断面積は群ごとの巨視的な値 Σ N σ（群の中で一定）。
/// 散乱の α は、範囲 `groups` の減速能 Σ ξ_k N_k σs,k を保つ値（ξ は核種の α の ξ を、範囲の群の σs をレサジーの幅で平均した
/// N σs で重みづけした平均）。まとめた核種の散乱は窓の幅が 1 つになり、同じ α の核種だけなら厳密（ADR-0185）
fn lump<'a>(
    library: &Library,
    members: &[(&str, f64, &Shielded)],
    bounds: &'a [f64],
    groups: Range<usize>,
) -> Result<Data<'a>, UfgError> {
    let mut alphas = Vec::with_capacity(members.len());
    let (mut power, mut scattering) = (0.0, 0.0);
    for &(name, density, s) in members {
        let awr = library
            .nuclides
            .iter()
            .find(|n| n.name == name)
            .ok_or_else(|| LookupError::Nuclide(name.into()))?
            .awr;
        let alpha = scattering_alpha(name, awr)?;
        let (mut sum, mut width) = (0.0, 0.0);
        for g in groups.clone() {
            let du = (bounds[g] / bounds[g + 1]).ln();
            sum += s.elastic[g] * du;
            width += du;
        }
        let weight = density * sum / width;
        power += weight * mean_lethargy_gain(alpha);
        scattering += weight;
        alphas.push(alpha);
    }
    let lightest = alphas.iter().copied().fold(f64::INFINITY, f64::min);
    let heaviest = alphas.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    // 散乱のない背景（個数密度がすべて 0 など）は散乱の線源に効かないので、α は最も重い核種のものでよい
    let alpha = if scattering > 0.0 {
        let xi = power / scattering;
        // ξ は α について単調に減るので、まとめた核種の α の間で二分法で解く
        let (mut a, mut b) = (lightest, heaviest);
        for _ in 0..200 {
            let mid = 0.5 * (a + b);
            if mean_lethargy_gain(mid) > xi {
                a = mid;
            } else {
                b = mid;
            }
        }
        0.5 * (a + b)
    } else {
        heaviest
    };
    let ng = bounds.len() - 1;
    let rows = (0..ng)
        .map(|g| {
            let sum = |of: fn(&Shielded) -> &Vec<f64>| -> f64 {
                members.iter().map(|&(_, n, s)| n * of(s)[g]).sum()
            };
            [sum(|s| &s.total), sum(|s| &s.elastic), 0.0, 0.0, 0.0]
        })
        .collect();
    Ok(Data {
        name: "background",
        gain: 1.0 / (1.0 - alpha),
        width: -alpha.ln(),
        source: Source::Groups { bounds, rows },
    })
}

/// 核種ごとの格子の点の数を debug で出す。`own` は範囲の中の自分の点の数、`added` は、データの順に足したときに格子の和に
/// 加わる（それまでの核種と群の境界にない）エネルギーの数
fn log_nuclide_points(data: &[Data], group_bounds: &[f64], inside: impl Fn(f64) -> bool) {
    let mut seen: HashSet<u64> = group_bounds.iter().map(|e| e.to_bits()).collect();
    for d in data {
        let own: Vec<f64> = d
            .energies()
            .iter()
            .copied()
            .filter(|&e| inside(e))
            .collect();
        let added = own.iter().filter(|e| seen.insert(e.to_bits())).count();
        tracing::debug!(name = d.name, own = own.len(), added, "UFG nuclide points");
    }
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
    // 物質ごとの軸: 核種の並びが同じ物質の最初のものの軸。リングの背景の項は、どのリングのものも同じ核種とみなす（燃焼で
    // リングの組成が違っても、背景はリングごとの項になるため）
    let nuclides = |i: usize| {
        entries[i].iter().map(|&(k, _)| match data[k].source {
            Source::Pointwise(_) => Some(k),
            Source::Groups { .. } => None,
        })
    };
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
    let mut cursors: Vec<usize> = data.iter().map(Data::first_cursor).collect();
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

/// 核種 `d` のポイントワイズの断面積 [全, 弾性散乱, 核分裂, νσf, 捕獲] をエネルギー `e` で lin-lin に補間する（背景は `e` を含む群の値）。
/// `cursor` はエネルギーの降順に引くための位置（[`Data::first_cursor`] から始める）。不連続点（背景では群の境界）では、`below` なら
/// 低い側、そうでなければ高い側の極限の値になる。
fn evaluate(d: &Data, cursor: &mut usize, e: f64, below: bool) -> Result<Row, UfgError> {
    let table = match &d.source {
        Source::Pointwise(t) => t,
        Source::Groups { bounds, rows } => {
            // 群 g は (bounds[g + 1], bounds[g])。下端ちょうどの点は、低い側の極限なら次の群
            while *cursor + 1 < rows.len()
                && (e < bounds[*cursor + 1] || (below && e == bounds[*cursor + 1]))
            {
                *cursor += 1;
            }
            return Ok(rows[*cursor]);
        }
    };
    let x = &table.pointwise.energy;
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
    let p = &table.pointwise;
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
                    d.elastic(g) * ((-(big_u - y)).exp() - (-(big_u - x)).exp())
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
