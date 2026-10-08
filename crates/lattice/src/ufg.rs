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

use rayon::prelude::*;

use xlc_format::library::{Library, LookupError, Pointwise, Shielded, TemperatureData};

/// 1つの窓（散乱で届くレサジーの幅 Δ）に入れる格子点の最少の数。xs のライブラリの重みの減速計算の同名の定数とは
/// 独立に選べる、この計算の格子の細かさの設定である。同じ減速の式を同じ格子の作り方で解くので、xs で収束を確かめた
/// 値（xs の `finer_grids_do_not_change_the_flux`）をそのまま使っている。この計算での収束はまだ確かめていない
const WINDOW_POINTS: usize = 16;

use crate::cp::{CpError, Cylinder};
use crate::input::MaterialInput;
use crate::lu::solve;

/// セルの UFG の範囲の上限の頭打ち (eV)。U-238 の分離共鳴の上端で、これより上は U-238 の非弾性散乱（45 keV から）と
/// 核分裂の線源が効き、弾性散乱だけの減速計算が成り立たない（docs/plans/m2-assembly.md）
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
    let Some(top) = top else {
        return Ok(None);
    };
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
    Ok(Some([low, high]))
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

/// セル `cell` の、範囲 `range`（[`cell_range`]）の減速計算。`rings` はリングごとの物質（真空は `None`）。
pub fn slowing_down(
    library: &Library,
    cell: &Cylinder,
    rings: &[Option<&MaterialInput>],
    range: [f64; 2],
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
            .map(|ring| vec![vec![[0.0; 4]; ng]; ring.len()])
            .collect(),
    };

    let mut cursors: Vec<usize> = data.iter().map(|d| d.pointwise().energy.len()).collect();
    for start in (0..u.len()).step_by(CHUNK) {
        let chunk = start..(start + CHUNK).min(u.len());
        // [点][核種] のポイントワイズの断面積 [全, 弾性散乱, 核分裂, νσf]
        let mut values: Vec<Vec<[f64; 4]>> = Vec::with_capacity(chunk.len());
        for &(e, below) in &points[chunk.clone()] {
            let row = data
                .iter()
                .zip(cursors.iter_mut())
                .map(|(d, c)| evaluate(d, c, e, below))
                .collect::<Result<Vec<_>, _>>()?;
            values.push(row);
        }
        let tracks = values
            .par_iter()
            .map(|row| {
                let totals: Vec<f64> = entries
                    .iter()
                    .map(|ring| ring.iter().map(|&(k, d)| d * row[k][0]).sum())
                    .collect();
                cell.tracks(&totals)
            })
            .collect::<Result<Vec<_>, _>>()?;
        for ((m, row), w) in chunk.zip(values).zip(&tracks) {
            let (e, below) = points[m];
            sweep.point(m, row, w, e)?;
            // 群の境界が不連続点なら、高い側の格子点で次の群に移る（低い側との間は幅 0）
            if sweep.group < ng && !below && e == group_bounds[sweep.group + 1] {
                sweep.group += 1;
            }
        }
    }

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
            };
            let finite = |v: &Vec<f64>| v.iter().position(|x| !x.is_finite());
            let bad = [
                &shielded.total,
                &shielded.elastic,
                &shielded.fission,
                &shielded.nu_fission,
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

/// 核種 `d` のポイントワイズの断面積 [全, 弾性散乱, 核分裂, νσf] をエネルギー `e` で lin-lin に補間する。`cursor` はエネルギーの
/// 降順に引くための位置。不連続点では、`below` なら低い側、そうでなければ高い側の極限の値になる。
fn evaluate(d: &Data, cursor: &mut usize, e: f64, below: bool) -> Result<[f64; 4], UfgError> {
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
    previous: Vec<[f64; 4]>,
    /// 今の区間の群（`groups` の中の番号）
    group: usize,
    group_flux: Vec<Vec<f64>>,
    /// [リング][項][群] の反応率 ∫ F σ du [全, 弾性散乱, 核分裂, νσf]
    rates: Vec<Vec<Vec<[f64; 4]>>>,
}

impl Sweep<'_> {
    /// 格子点 m の中性子束を求め、窓の積分と群の積分を進める。`values` はその点のポイントワイズの断面積、`tracks` は飛程。
    fn point(
        &mut self,
        m: usize,
        values: Vec<[f64; 4]>,
        tracks: &[Vec<f64>],
        energy: f64,
    ) -> Result<(), UfgError> {
        let n = self.volumes.len();
        for (column, v) in self.elastic.iter_mut().zip(&values) {
            column.push(v[1]);
        }
        // 窓の積分のうち既知の部分 [リング][項] と、その点の F_i の係数
        let mut known = self.windows.clone();
        let mut own = vec![vec![0.0; 0]; n];
        if m > 0 {
            let (u0, u1) = (self.u[m - 1], self.u[m]);
            let decay = (-(u1 - u0)).exp();
            let (ca, cb) = segment(u0, u1, u1);
            // 核種ごとの窓の下端の区間 [u0 − Δ, u1 − Δ] の積分の、格子点の値の重み（範囲の上の部分は値で）
            let above = self.above;
            let slabs: Vec<(f64, Vec<(usize, f64)>)> = self
                .data
                .iter()
                .zip(self.lower.iter_mut())
                .map(|(d, lower)| {
                    let (a, b) = (u0 - d.width, u1 - d.width);
                    let outside = if a < 0.0 {
                        above.integral(d, a, b.min(0.0), u1)
                    } else {
                        0.0
                    };
                    (outside, slab(self.u, a, b, u1, lower))
                })
                .collect();
            for (i, ring) in self.entries.iter().enumerate() {
                let f = &self.flux[i];
                own[i] = ring.iter().map(|&(k, _)| cb * self.elastic[k][m]).collect();
                for (e, &(k, _)) in ring.iter().enumerate() {
                    let (outside, weights) = &slabs[k];
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
            own = self
                .entries
                .iter()
                .map(|ring| vec![0.0; ring.len()])
                .collect();
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
        // W は対称なので W_ij = W_ji
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
                    for x in 0..4 {
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
/// （u ≥ 0）の部分を、格子点の σs F にかける重みで返す。`lower` は下端を含む区間の始まりの格子点（進めていく）。
fn slab(u: &[f64], a: f64, b: f64, big_u: f64, lower: &mut usize) -> Vec<(usize, f64)> {
    let mut weights = Vec::new();
    let mut x = a.max(0.0);
    if b <= x {
        return weights;
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
    weights
}

#[cfg(test)]
mod tests {
    use super::*;
    use xlc_format::group::GroupStructure;
    use xlc_format::library::{NuclideData, Pointwise, Resolved, Shielded, TemperatureData};

    const TEMPERATURE: f64 = 600.0;
    const BOUNDS: [f64; 7] = [2.0e7, 2.0e4, 1.0e4, 1.0e3, 1.0e2, 4.0, 1.0e-5];

    /// ポイントワイズの断面積 (energy, total, elastic, fission) を持つ核種。群定数は 1/E の重みの平均
    struct Spec {
        name: &'static str,
        awr: f64,
        top: Option<f64>,
        points: Vec<[f64; 4]>,
    }

    fn constant(name: &'static str, awr: f64, total: f64, elastic: f64, fission: f64) -> Spec {
        let points = [4.0, 2.0e4].map(|e| [e, total, elastic, fission]).to_vec();
        Spec {
            name,
            awr,
            top: None,
            points,
        }
    }

    fn library(specs: Vec<Spec>) -> Library {
        let groups = GroupStructure::new(BOUNDS.to_vec()).unwrap();
        let ng = BOUNDS.len() - 1;
        let nuclides = specs
            .into_iter()
            .map(|spec| {
                let column = |k: usize| spec.points.iter().map(|p| p[k]).collect::<Vec<_>>();
                // 群定数はポイントワイズの断面積の両端の値を使う（範囲の外の減速源に使う）
                let first = spec.points[0];
                let group = |k: usize| vec![first[k]; ng];
                NuclideData {
                    name: spec.name.into(),
                    za: 0.0,
                    awr: spec.awr,
                    resolved: spec.top.map(|top| Resolved {
                        top,
                        potential: 10.0,
                    }),
                    thermal_scattering: false,
                    temperatures: vec![TemperatureData {
                        temperature: TEMPERATURE,
                        chi: vec![0.0; ng],
                        scatter: vec![vec![vec![0.0; ng]; ng]],
                        shielded: vec![Shielded {
                            total: group(1),
                            elastic: group(2),
                            fission: group(3),
                            nu_fission: group(3).iter().map(|f| 2.4 * f).collect(),
                        }],
                        pointwise: Pointwise {
                            energy: column(0),
                            total: column(1),
                            elastic: column(2),
                            fission: column(3),
                            nu_fission: column(3).iter().map(|f| 2.4 * f).collect(),
                        },
                    }],
                }
            })
            .collect();
        Library {
            groups,
            sigma0: vec![1.0e10],
            pointwise_range: Some([4.0, 2.0e4]),
            nuclides,
        }
    }

    fn material(composition: &[(&str, f64)]) -> MaterialInput {
        MaterialInput {
            name: "m".into(),
            temperature: TEMPERATURE,
            composition: composition
                .iter()
                .map(|&(n, d)| (n.to_string(), d))
                .collect(),
        }
    }

    /// 燃料・ギャップ・被覆管・水
    const RADII: [f64; 4] = [0.4096, 0.418, 0.475, 0.711];

    #[test]
    fn range_follows_the_resolved_top_of_the_fuel() {
        let fuel = material(&[("U-235", 0.01)]);
        let clad = material(&[("Zr-90", 0.04)]);
        let water = material(&[("H-1", 0.05)]);
        let range = |top: f64| {
            let mut fissile = constant("U-235", 233.0, 10.0, 8.0, 1.0);
            fissile.top = Some(top);
            // 燃料の外の核種の分離共鳴は範囲を決めない
            let mut zirconium = constant("Zr-90", 89.0, 6.0, 6.0, 0.0);
            zirconium.top = Some(1.0e6);
            let lib = library(vec![
                fissile,
                zirconium,
                constant("H-1", 0.99917, 20.0, 20.0, 0.0),
            ]);
            cell_range(&lib, &[Some(&fuel), None, Some(&clad), Some(&water)])
                .unwrap()
                .unwrap()
        };
        // 上端以下で最も高い群の境界
        assert_eq!(range(3.0e3), [4.0, 1.0e3]);
        // 20 keV で頭打ちにする
        assert_eq!(range(5.0e4), [4.0, 2.0e4]);
    }

    #[test]
    fn no_range_without_resolved_fuel() {
        // 燃料がない、または燃料に分離共鳴を持つ核種がなければ、減速計算する範囲はない
        let lib = library(vec![
            constant("H-1", 0.99917, 20.0, 20.0, 0.0),
            constant("U-235", 233.0, 10.0, 8.0, 1.0),
        ]);
        let water = material(&[("H-1", 0.05)]);
        let fuel = material(&[("U-235", 0.01)]);
        assert_eq!(cell_range(&lib, &[Some(&water)]).unwrap(), None);
        assert_eq!(
            cell_range(&lib, &[Some(&fuel), Some(&water)]).unwrap(),
            None
        );
    }

    #[test]
    fn resolved_fuel_needs_pointwise_data() {
        // 燃料に分離共鳴があるのにポイントワイズの断面積のないライブラリはエラー。燃料のないセルは範囲がないだけでよい
        let mut fissile = constant("U-235", 233.0, 10.0, 8.0, 1.0);
        fissile.top = Some(3.0e3);
        let mut lib = library(vec![fissile, constant("H-1", 0.99917, 20.0, 20.0, 0.0)]);
        lib.pointwise_range = None;
        let water = material(&[("H-1", 0.05)]);
        let fuel = material(&[("U-235", 0.01)]);
        assert_eq!(cell_range(&lib, &[Some(&water)]).unwrap(), None);
        assert!(matches!(
            cell_range(&lib, &[Some(&fuel), Some(&water)]),
            Err(UfgError::NoPointwise)
        ));
    }

    #[test]
    fn scatterers_keep_the_asymptotic_spectrum() {
        // 吸収がなく断面積が一定なら、どのリングでも中性子束はレサジーあたり一定（1/E）のまま。減速材は、窓が
        // 群の上端（中性子のない範囲）まで届かない H-2 にする（H-1 では範囲の上からの線源が e^{−7.6} だけ足りない）
        let lib = library(vec![
            constant("U-238", 236.0, 11.0, 11.0, 0.0),
            constant("Zr-90", 89.0, 6.5, 6.5, 0.0),
            constant("H-2", 1.9968, 3.4, 3.4, 0.0),
            constant("O-16", 15.86, 3.8, 3.8, 0.0),
        ]);
        let fuel = material(&[("U-238", 0.022), ("O-16", 0.045)]);
        let clad = material(&[("Zr-90", 0.043)]);
        let water = material(&[("H-2", 0.05), ("O-16", 0.025)]);
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        let spectrum = slowing_down(
            &lib,
            &cell,
            &[Some(&fuel), None, Some(&clad), Some(&water)],
            [4.0, 1.0e4],
        )
        .unwrap();
        assert_eq!(spectrum.groups, 2..5);
        for (ring, flux) in spectrum.flux.iter().enumerate() {
            for (k, f) in flux.iter().enumerate() {
                let g = spectrum.groups.start + k;
                let per_lethargy = f / (BOUNDS[g] / BOUNDS[g + 1]).ln();
                assert!(
                    (per_lethargy - 1.0).abs() < 1e-9,
                    "{ring} {g}: {per_lethargy}"
                );
            }
        }
        // 群定数はポイントワイズの断面積の値のまま
        assert!((spectrum.nuclides[0][0].total[1] - 11.0).abs() < 1e-12);
        assert!(spectrum.nuclides[1].is_empty());
    }

    #[test]
    fn discontinuities_keep_both_limits() {
        // 500 eV で 10 b から 20 b に跳ぶ微量の核種: 中性子束は 1/E のままなので、群平均は両側の値の
        // レサジーの重みの平均になる
        let step = Spec {
            name: "Pu-239",
            awr: 236.999,
            top: None,
            points: vec![
                [4.0, 10.0, 0.0, 0.0],
                [500.0, 10.0, 0.0, 0.0],
                [500.0, 20.0, 0.0, 0.0],
                [2.0e4, 20.0, 0.0, 0.0],
            ],
        };
        let lib = library(vec![
            constant("U-238", 236.0, 11.0, 11.0, 0.0),
            constant("H-2", 1.9968, 3.4, 3.4, 0.0),
            step,
        ]);
        let fuel = material(&[("U-238", 0.022), ("Pu-239", 1.0e-12)]);
        let water = material(&[("H-2", 0.05)]);
        let cell = Cylinder::new(vec![0.4, 0.7]).unwrap();
        let spectrum =
            slowing_down(&lib, &cell, &[Some(&fuel), Some(&water)], [4.0, 1.0e4]).unwrap();
        // 1 keV〜100 eV の群
        let k = 3 - spectrum.groups.start;
        let expected = (10.0 * 5.0f64.ln() + 20.0 * 2.0f64.ln()) / 10.0f64.ln();
        let total = spectrum.nuclides[0][1].total[k];
        assert!(
            (total - expected).abs() < 1e-9 * expected,
            "{total} vs {expected}"
        );
    }

    #[test]
    fn hydrogen_slowing_down_matches_the_analytic_solution() {
        // H と吸収だけの核種の均質な媒質: 衝突密度 q は dq/du = −r q で減る。
        // r = 1 − Σs/((1 − α)Σt)、q(0) は範囲の上（中性子束 1/E、群の上端より上は 0）からの散乱
        let awr: f64 = 0.99917;
        let lib = library(vec![
            constant("H-1", awr, 20.0, 20.0, 0.0),
            constant("Pu-240", 238.0, 2.0, 0.0, 0.0),
        ]);
        let mix = material(&[("H-1", 0.05), ("Pu-240", 0.01)]);
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        let spectrum = slowing_down(&lib, &cell, &[Some(&mix); 4], [4.0, 1.0e4]).unwrap();
        let alpha = ((awr - 1.0) / (awr + 1.0)).powi(2);
        let (sigma_s, sigma_t) = (0.05 * 20.0, 0.05 * 20.0 + 0.01 * 2.0);
        let r = 1.0 - sigma_s / ((1.0 - alpha) * sigma_t);
        let q0 = sigma_s * (1.0 - 1.0e4 / 2.0e7) / (1.0 - alpha);
        for flux in &spectrum.flux {
            for (k, f) in flux.iter().enumerate() {
                let g = spectrum.groups.start + k;
                let (ua, ub) = ((1.0e4 / BOUNDS[g]).ln(), (1.0e4 / BOUNDS[g + 1]).ln());
                let expected = q0 / (sigma_t * r) * ((-r * ua).exp() - (-r * ub).exp());
                assert!(
                    (f - expected).abs() < 1e-7 * expected,
                    "{g}: {f} vs {expected}"
                );
            }
        }
    }

    #[test]
    fn resonance_is_shielded_more_inside_the_fuel() {
        // 200 eV の強い吸収の共鳴: 燃料の内側ほど中性子束がへこみ、実効断面積が小さい
        let resonance = |e: f64| 1.0e4 / (1.0 + ((e - 200.0) / 0.5).powi(2));
        let mut energies: Vec<f64> = (0..=400).map(|k| 180.0 + 0.1 * k as f64).collect();
        energies.splice(0..0, [4.0, 100.0]);
        energies.extend([1.0e3, 2.0e4]);
        let points = energies
            .iter()
            .map(|&e| [e, 11.0 + resonance(e), 11.0, 0.0])
            .collect();
        let lib = library(vec![
            Spec {
                name: "U-238",
                awr: 236.0,
                top: Some(2.0e4),
                points,
            },
            constant("U-235", 233.0, 10.0, 8.0, 1.0),
            constant("H-1", 0.99917, 20.0, 20.0, 0.0),
            constant("O-16", 15.86, 3.8, 3.8, 0.0),
        ]);
        let fuel = material(&[("U-238", 0.022), ("U-235", 0.0005), ("O-16", 0.045)]);
        let water = material(&[("H-1", 0.05), ("O-16", 0.025)]);
        let cell = Cylinder::new(vec![0.2, 0.35, 0.4096, 0.418, 0.711]).unwrap();
        let rings = [Some(&fuel), Some(&fuel), Some(&fuel), None, Some(&water)];
        let spectrum = slowing_down(&lib, &cell, &rings, [4.0, 1.0e4]).unwrap();
        // 100 eV〜1 keV の群
        let k = 3 - spectrum.groups.start;
        let total = |ring: usize| spectrum.nuclides[ring][0].total[k];
        // 無限希釈（1/E の重み）の平均
        let dilute = {
            let (mut rate, mut flux) = (0.0, 0.0);
            for w in energies
                .windows(2)
                .filter(|w| w[0] >= 100.0 && w[1] <= 1.0e3)
            {
                let du = (w[1] / w[0]).ln();
                rate += 0.5 * du * (resonance(w[0]) + resonance(w[1]));
                flux += du;
            }
            11.0 + rate / flux
        };
        assert!(
            total(0) < total(1) && total(1) < total(2),
            "{:?}",
            [total(0), total(1), total(2)]
        );
        assert!(total(2) < dilute, "{} vs {dilute}", total(2));
    }

    #[test]
    fn rejects_a_range_off_the_group_bounds() {
        let lib = library(vec![constant("H-1", 0.99917, 20.0, 20.0, 0.0)]);
        let water = material(&[("H-1", 0.05)]);
        let cell = Cylinder::new(vec![0.5]).unwrap();
        for range in [[4.0, 5.0e3], [3.0, 1.0e4], [1.0e2, 1.0e2]] {
            assert!(
                slowing_down(&lib, &cell, &[Some(&water)], range).is_err(),
                "{range:?}"
            );
        }
        // リングの数が合わない
        assert!(slowing_down(&lib, &cell, &[Some(&water); 2], [4.0, 1.0e4]).is_err());
    }
}
