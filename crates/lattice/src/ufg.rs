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
//!
//! [`Scattering::FreeGas`] では、分離共鳴を持たない核種の散乱核を温度 T の自由ガスの核にする（[`free_gas`]、ADR-0203）。
//! [`FreeGas::heavy`] を選ぶと、分離共鳴を持つ核種も、その上限より下で同じ断面積一定の自由ガスの核にする（ADR-0204）。
//! [`FreeGas::sab`] を選ぶと、S(α, β) の核種（H in H2O の H）は、S(α, β) の上限より下の入射でライブラリの S(α, β) の移行で
//! 散乱する（[`free_gas::SabSource`]、ADR-0205）。
//! 上方散乱でまだ求めていない低いエネルギーの中性子束を使うので、自由ガスの核の届く範囲を掃き直す。範囲の下限より下からの
//! 上方散乱は、その中性子束を範囲の最後の格子点の値（1/E）とする（ADR-0203）。

use std::borrow::Cow;
use std::collections::HashSet;
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
use free_gas::{Gas, Thermal};

mod free_gas;
mod transfer;
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
/// 間引きで巨視的な断面積を並列に足す点の塊の大きさ
const THIN_POINTS: usize = 4096;

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
    #[error(
        "{name}: `sab` needs the S(α, β) transfer of the library, which xs leaves out when the S(α, β) limit does not reach the UFG range or the up-scattering reach goes below 0 eV (and libraries made before ADR-0205 lack)"
    )]
    NoSabTransfer { name: String },
    #[error(
        "{name}: `sab` needs the free gas of this S(α, β) nuclide with resolved resonances (`heavy`)"
    )]
    SabWithoutFreeGas { name: String },
    #[error("{name}: the S(α, β) limit {floor} eV must lie below the UFG range top {high} eV")]
    SabAboveRange { name: String, floor: f64, high: f64 },
    #[error(
        "{name}: the {what} of the S(α, β) transfer is {value} eV where the UFG needs {needed} eV"
    )]
    SabCoverage {
        name: String,
        what: &'static str,
        value: f64,
        needed: f64,
    },
    /// S(α, β) の散乱源の小区間が表の外か空（表の範囲は核種のデータを引くときに確かめるので、内部の不整合）
    #[error(
        "internal error: the S(α, β) source at {energy} eV has no table share at transfer {transfer} eV"
    )]
    SabTable { energy: f64, transfer: f64 },
    #[error(
        "`ufg-transfer = \"ring-flux\"` counts the transfers on the full grid only (not with `ufg-tolerance`)"
    )]
    ThinnedTransfers,
    #[error("the up-scattering of the free-gas kernel did not converge in {0} sweeps")]
    NotConverged(usize),
    #[error("ring {ring}: non-finite group constants in group {group}")]
    NonFinite { ring: usize, group: usize },
    /// 核種を含むリングの縮約の結果がない（内部の不整合）
    #[error("ring {ring}: internal error, nuclide {nuclide} was not collapsed")]
    Uncollapsed { ring: usize, nuclide: usize },
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

/// UFG の減速計算の散乱核の、標的の扱い（ADR-0203）。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Scattering {
    /// どの核種も静止した標的（ADR-0070）
    #[default]
    Static,
    /// 分離共鳴を持たない核種（ライブラリの `resolved` が nil）を温度 T の自由ガスにする。分離共鳴を持つ核種は、
    /// [`FreeGas::heavy`] を選べばその上限より下で自由ガス、選ばなければ静止した標的
    FreeGas(FreeGas),
}

/// 自由ガスの核の計算の条件（ADR-0203）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeGas {
    /// 自由ガスの核を使う入射エネルギーの上限（核種の kT の倍数）。それより上は静止した標的。`None` は範囲の全体
    pub upper: Option<f64>,
    /// 自由ガスと静止した標的の核の差を積分する区間の細かさ（上方散乱の裾は kT あたり、下端のぼけはその幅あたりの区間の数）
    pub cells: usize,
    /// 上方散乱のスイープの回数。`None` は、中性子束の相対変化が [`SWEEP_TOLERANCE`] 以下になるまで（[`MAX_SWEEPS`] 回まで）
    pub sweeps: Option<usize>,
    /// 分離共鳴を持つ核種の自由ガスの核を使う入射エネルギーの上限（核種の kT の倍数。無限大は範囲の全体）。`None` は分離共鳴を
    /// 持つ核種を静止した標的にする（ADR-0204）
    pub heavy: Option<f64>,
    /// S(α, β) を持つ核種（H in H2O の H）を、S(α, β) の上限より下の入射でライブラリの S(α, β) の移行で散乱させる。false は、その
    /// 核種も自由ガスの核にする（ADR-0205）
    pub sab: bool,
}

impl Default for FreeGas {
    /// ADR-0203 の既定
    fn default() -> Self {
        Self {
            upper: None,
            cells: 4,
            sweeps: None,
            heavy: None,
            sab: false,
        }
    }
}

/// 上方散乱のスイープを打ち切る、中性子束の相対変化の最大
pub const SWEEP_TOLERANCE: f64 = 1.0e-10;
/// 上方散乱のスイープの回数の上限
pub const MAX_SWEEPS: usize = 50;

/// 1つの核種（名前と温度の組）の減速計算のデータ
struct Data<'a> {
    name: &'a str,
    /// 1/(1 − α)
    gain: f64,
    /// 窓の幅 Δ = ln(1/α)（α = 0 なら無限）
    width: f64,
    /// その温度の群定数とポイントワイズの断面積（S(α, β) の核種は内挿したもの）
    table: Cow<'a, TemperatureData>,
    /// 自由ガスの核で散乱するときの核（[`Scattering::FreeGas`] の、分離共鳴を持たない核種と、[`FreeGas::heavy`] を選んだときの
    /// 分離共鳴を持つ核種）
    gas: Option<Gas>,
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
/// 飛程は既定の求め方（[`Collision::default`]）で求める。
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
        range,
        Collision::default(),
        None,
        Scattering::Static,
    )
}

/// [`slowing_down`] の、飛程の求め方を `collision` で選び、`thinning` が `Some(許容差)` なら格子を間引き（[`thin`]。
/// 核種の群定数は、核種自身のポイントワイズの断面積の格子でリングの中性子束を内挿して縮約する（[`collapse`]。ADR-0185））、
/// 散乱核の標的の扱いを `scattering` で選ぶもの。
pub fn slowing_down_with(
    library: &Library,
    cell: &Cylinder,
    rings: &[Option<&MaterialInput>],
    range: [f64; 2],
    collision: Collision,
    thinning: Option<f64>,
    scattering: Scattering,
) -> Result<CellSpectrum, UfgError> {
    solve_cell(
        library, cell, rings, range, collision, thinning, scattering, None, None,
    )
}

/// 1 つの核種の、UFG の範囲の群から出る弾性散乱の移行（ADR-0207。[`slowing_down_transfers`]）。
#[derive(Debug, Clone, PartialEq)]
pub struct NuclideTransfer {
    /// 入射の群（ライブラリの群の番号。減速計算した群）
    pub groups: Range<usize>,
    /// [入射の群 − `groups.start`][出射の群（ライブラリのすべての群）] の移行の微視的断面積 (barn)
    pub rows: Vec<Vec<f64>>,
}

/// [`slowing_down_with`] と同じ減速計算を解き、[リング][核種（物質の組成の順）] の弾性散乱の移行の行も返す（ADR-0207）。
#[allow(clippy::too_many_arguments)]
pub fn slowing_down_transfers(
    library: &Library,
    cell: &Cylinder,
    rings: &[Option<&MaterialInput>],
    range: [f64; 2],
    collision: Collision,
    thinning: Option<f64>,
    scattering: Scattering,
) -> Result<(CellSpectrum, Vec<Vec<NuclideTransfer>>), UfgError> {
    let mut transfers = Vec::new();
    let spectrum = solve_cell(
        library,
        cell,
        rings,
        range,
        collision,
        thinning,
        scattering,
        None,
        Some(&mut transfers),
    )?;
    Ok((spectrum, transfers))
}

/// 減速計算の格子点ごとのリングの中性子束（[`slowing_down_fine`]）。
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FineFlux {
    /// 格子点のエネルギー (eV)、降順（不連続点は同じエネルギーの 2 点）
    pub energies: Vec<f64>,
    /// [リング][格子点] のレサジーあたりの中性子束（[`CellSpectrum::flux`] と同じ規格化）
    pub flux: Vec<Vec<f64>>,
}

/// [`slowing_down_with`] と同じ減速計算を解き、格子点ごとのリングの中性子束も返す（群の中の重みを参照解の細かい
/// 中性子束と比べるため。格子点の数 × リングの数の中性子束を保つので、自己遮蔽の経路では使わない）。
#[allow(clippy::too_many_arguments)]
pub fn slowing_down_fine(
    library: &Library,
    cell: &Cylinder,
    rings: &[Option<&MaterialInput>],
    range: [f64; 2],
    collision: Collision,
    thinning: Option<f64>,
    scattering: Scattering,
) -> Result<(CellSpectrum, FineFlux), UfgError> {
    let mut fine = FineFlux::default();
    let spectrum = solve_cell(
        library,
        cell,
        rings,
        range,
        collision,
        thinning,
        scattering,
        Some(&mut fine),
        None,
    )?;
    Ok((spectrum, fine))
}

#[tracing::instrument(name = "ufg", level = "debug", skip_all)]
#[allow(clippy::too_many_arguments)]
fn solve_cell(
    library: &Library,
    cell: &Cylinder,
    rings: &[Option<&MaterialInput>],
    range: [f64; 2],
    collision: Collision,
    thinning: Option<f64>,
    scattering: Scattering,
    fine: Option<&mut FineFlux>,
    transfers: Option<&mut Vec<Vec<NuclideTransfer>>>,
) -> Result<CellSpectrum, UfgError> {
    // 間引いた格子では群定数を核種自身の格子で縮約するので、間引いた格子で数える移行と格子が揃わない
    if transfers.is_some() && thinning.is_some() {
        return Err(UfgError::ThinnedTransfers);
    }
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

    let started = Instant::now();
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
                    // 分離共鳴を持つ核種は `heavy` を選んだときだけ、その上限で自由ガスにする（ADR-0204）。`sab` を選べば、
                    // S(α, β) の核種はその上限より下で S(α, β) の移行で散乱する（ADR-0205）
                    let gas = match scattering {
                        Scattering::FreeGas(options) => {
                            let upper = match &d.resolved {
                                None => Some(options.upper),
                                Some(_) => options.heavy.map(Some),
                            };
                            let gas = upper.map(|upper| Gas::new(d.awr, t.temperature, upper));
                            // S(α, β) の移行は自由ガスの核の上限より上の入射を自由ガスに残すので、自由ガスの核を持たない核種
                            // （`heavy` のない分離共鳴を持つ核種）では使えない（黙って静止した標的にしない）
                            if gas.is_none() && options.sab && d.thermal_scattering {
                                return Err(UfgError::SabWithoutFreeGas { name: name.clone() });
                            }
                            match gas {
                                Some(gas) if options.sab && d.thermal_scattering => {
                                    let table = t.sab_transfer.as_ref().ok_or_else(|| {
                                        UfgError::NoSabTransfer { name: name.clone() }
                                    })?;
                                    check_sab_coverage(name, table, gas.kt(), [low, high])?;
                                    Some(gas.below_sab(table))
                                }
                                gas => gas,
                            }
                        }
                        Scattering::Static => None,
                    };
                    data.push(Data {
                        name: &d.name,
                        gain: 1.0 / (1.0 - alpha),
                        width: -alpha.ln(),
                        table: t,
                        gas,
                    });
                    keys.push((name, temperature));
                    data.len() - 1
                }
            };
            ring.push((k, *density));
        }
        entries.push(ring);
    }

    let looked_up = Instant::now();
    // 格子（エネルギーの降順）: ポイントワイズの断面積の格子と群の境界に、どの核種の窓にも WINDOW_POINTS 点が入るよう点を足す
    let group_bounds = &bounds[top..=bottom];
    let inside = |e: f64| e > low && e < high;
    // S(α, β) の上限では、静止した標的の窓の σs/N が上で値を持ち、下で 0 になる（その下は S(α, β) の移行で扱う）ので、両側の
    // 極限の 2 つの格子点を置く（ADR-0205）
    let floors: Vec<f64> = data
        .iter()
        .filter_map(|d| d.gas.and_then(|g| g.floor()))
        .filter(|&f| inside(f))
        .collect();
    let mut energies: Vec<f64> = data
        .iter()
        .flat_map(|d| d.pointwise().energy.iter().copied())
        .filter(|&e| inside(e))
        .chain(group_bounds.iter().copied())
        .chain(floors.iter().copied())
        .collect();
    // 全順序なので、並列の整列でも並びは同じ
    energies.par_sort_unstable_by(|a, b| b.total_cmp(a));
    energies.dedup();
    // 核種の格子の和（群の境界を含む）の点の数。ここから窓の点を足したものが格子
    let union = energies.len();
    if tracing::enabled!(tracing::Level::DEBUG) {
        log_nuclide_points(&data, group_bounds, inside);
    }
    // ポイントワイズの断面積の不連続点（同じエネルギーの2点）では、高い側と低い側の極限の2つの格子点を同じレサジーに置く。
    // 範囲の上端では低い側、下端では高い側の極限だけを使う
    let mut jumps: Vec<f64> = data
        .iter()
        .flat_map(|d| d.pointwise().energy.windows(2))
        .filter(|w| w[0] == w[1] && w[0] > low && w[0] < high)
        .map(|w| w[0])
        .collect();
    jumps.extend(&floors);
    jumps.sort_by(f64::total_cmp);
    jumps.dedup();
    // 間引くなら、群の境界と不連続点を残し、リングの物質ごとの巨視的な全断面積と弾性散乱断面積を許容差で再現できる点だけを
    // 残す（ADR-0185）
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
    // S(α, β) の上限は、丸めで同じレサジーになる点があっても残す（両側の極限の 2 点を置くため）
    let (u, energies) = lethargy_grid(&energies, high, step, |e| {
        group_bounds.contains(&e) || floors.contains(&e)
    });
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
        tally: thinning.is_none(),
        known: Vec::new(),
        own: entries.iter().map(|ring| vec![0.0; ring.len()]).collect(),
        slabs: vec![(0.0, Vec::new()); data.len()],
        thermal: None,
        extra: (vec![0.0; n], vec![0.0; n]),
    };

    // 表で内挿するなら、軸ごとの全断面積の範囲を格子点の全体で先に求めて表を作る（表で衝突確率が減らないセルは作らない）
    let grid_done = Instant::now();
    let table = match collision {
        Collision::Exact => None,
        Collision::Table { density } => track_table(cell, &data, &entries, &points, density)?,
    };
    let table_done = Instant::now();

    let mut cursors: Vec<usize> = data.iter().map(|d| d.pointwise().energy.len()).collect();
    // 段ごとの経過時間（ログの内訳。衝突確率は格子点の間で並列に求める）。ピンの間も並列のとき、衝突確率を待つ
    // スレッドがほかのピンの仕事を盗むと、その時間も衝突確率に入るので、段の割合は 1 スレッドで測る
    let mut times = [Duration::ZERO; 3];
    let pass = Pass {
        data: &data,
        entries: &entries,
        cell,
        table: &table,
        points: &points,
        group_bounds,
        ng,
    };
    // 自由ガスの核の届く範囲（start から）は、上方散乱が前のスイープの中性子束を使うので、収束するまで掃き直す
    let (thermal, sweeps) = match scattering {
        Scattering::Static => (None, None),
        Scattering::FreeGas(options) => {
            let gases: Vec<(usize, Gas)> = data
                .iter()
                .enumerate()
                .filter_map(|(k, d)| d.gas.map(|gas| (k, gas)))
                .collect();
            let energies: Vec<f64> = points.iter().map(|&(e, _)| e).collect();
            // D に掛ける σs は、上限より下の入射を補間する格子点からだけ引く（上限の低い分離共鳴を持つ核種で、範囲の全体を持たない）
            let sigma = gases
                .iter()
                .map(|&(k, gas)| {
                    let first = free_gas::first_below(&energies, gas.top());
                    let mut cursor = data[k].pointwise().energy.len();
                    let column = points[first..]
                        .iter()
                        .map(|&(e, below)| Ok(evaluate(&data[k], &mut cursor, e, below)?[1]))
                        .collect::<Result<Vec<f64>, UfgError>>()?;
                    Ok((first, column))
                })
                .collect::<Result<Vec<_>, UfgError>>()?;
            let rings = gases
                .iter()
                .map(|&(k, _)| {
                    entries
                        .iter()
                        .enumerate()
                        .flat_map(|(i, ring)| {
                            ring.iter()
                                .filter(move |&&(j, _)| j == k)
                                .map(move |&(_, d)| (i, d))
                        })
                        .collect()
                })
                .collect();
            // S(α, β) の移行で散乱する核種の表
            let tables = gases
                .iter()
                .map(|&(k, gas)| {
                    gas.floor()
                        .and(data[k].table.sab_transfer.as_ref())
                        .cloned()
                })
                .collect();
            let thermal = Thermal::new(
                gases.iter().map(|&(_, gas)| gas).collect(),
                sigma,
                rings,
                tables,
                options.cells,
                &energies,
                n,
            );
            (thermal, options.sweeps)
        }
    };
    sweep.thermal = thermal;
    let start = sweep.thermal.as_ref().map_or(u.len(), Thermal::start);
    pass.run(&mut sweep, &mut cursors, 0..start, &mut times)?;
    if sweep.thermal.is_some() {
        let saved = sweep.save();
        let saved_cursors = cursors.clone();
        let mut count = 0;
        loop {
            pass.run(&mut sweep, &mut cursors, start..u.len(), &mut times)?;
            count += 1;
            let change = sweep.thermal.as_mut().and_then(|t| t.finish(&sweep.flux));
            tracing::debug!(sweep = count, change, "UFG free-gas sweep");
            let done = match sweeps {
                Some(limit) => count >= limit,
                None => change.is_some_and(|c| c <= SWEEP_TOLERANCE),
            };
            if done {
                break;
            }
            if sweeps.is_none() && count >= MAX_SWEEPS {
                return Err(UfgError::NotConverged(MAX_SWEEPS));
            }
            sweep.restore(&saved);
            cursors.clone_from(&saved_cursors);
        }
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

    // 間引いた格子では核種自身の共鳴の点が落ちうるので、核種の格子で縮約する。核種ごとに、それを含むリングをまとめて 1 度に
    // 縮約する（[データ] の [(リング, 群定数)]）
    let mut collapsed: Vec<Vec<(usize, Shielded)>> = vec![Vec::new(); data.len()];
    if thinning.is_some() {
        collapsed = (0..data.len())
            .into_par_iter()
            .map(|k| {
                let rings: Vec<usize> = (0..n)
                    .filter(|&i| entries[i].iter().any(|&(j, _)| j == k))
                    .collect();
                let fluxes: Vec<&[f64]> = rings.iter().map(|&i| sweep.flux[i].as_slice()).collect();
                let group_fluxes: Vec<&[f64]> = rings
                    .iter()
                    .map(|&i| sweep.group_flux[i].as_slice())
                    .collect();
                let values = collapse(
                    &data[k],
                    &u,
                    &points,
                    &fluxes,
                    high,
                    group_bounds,
                    &group_fluxes,
                )?;
                Ok(rings.into_iter().zip(values).collect())
            })
            .collect::<Result<_, UfgError>>()?;
    }
    // 群定数
    let mut nuclides = Vec::with_capacity(n);
    for (i, ring) in sweep.rates.iter().enumerate() {
        let flux = &sweep.group_flux[i];
        let mut list = Vec::with_capacity(ring.len());
        for (rate, &(k, _)) in ring.iter().zip(&entries[i]) {
            let column =
                |x: usize| -> Vec<f64> { rate.iter().zip(flux).map(|(r, f)| r[x] / f).collect() };
            let shielded = if thinning.is_some() {
                collapsed[k]
                    .iter()
                    .find(|(r, _)| *r == i)
                    .map(|(_, s)| s.clone())
                    .ok_or(UfgError::Uncollapsed {
                        ring: i,
                        nuclide: k,
                    })?
            } else {
                Shielded {
                    total: column(0),
                    elastic: column(1),
                    fission: column(2),
                    nu_fission: column(3),
                    capture: column(4),
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
    if let Some(transfers) = transfers {
        // 自由ガスの核種の番号（`Thermal::new` に渡した順）
        let mut count = 0;
        let gas_of: Vec<Option<usize>> = data
            .iter()
            .map(|d| {
                d.gas.map(|_| {
                    count += 1;
                    count - 1
                })
            })
            .collect();
        let counted = Instant::now();
        *transfers = transfer::Counting {
            data: &data,
            entries: &entries,
            u: &u,
            points: &points,
            flux: &sweep.flux,
            elastic: &sweep.elastic,
            group_flux: &sweep.group_flux,
            nuclides: &nuclides,
            bounds,
            groups: top..bottom,
            high,
            thermal: sweep.thermal.as_ref(),
            gas_of: &gas_of,
        }
        .count()?;
        tracing::debug!(
            seconds = counted.elapsed().as_secs_f64(),
            "UFG transfers counted"
        );
    }
    // 頼まれたときだけ格子点ごとの中性子束を渡す（頼まれなければ、格子点の数 × リングの数の中性子束はここで捨てる）
    if let Some(fine) = fine {
        fine.energies = points.iter().map(|&(e, _)| e).collect();
        fine.flux = std::mem::take(&mut sweep.flux);
    }
    Ok(CellSpectrum {
        groups: top..bottom,
        flux: sweep.group_flux,
        nuclides,
        points: u.len(),
    })
}

/// 核種 `name` の S(α, β) の移行の表 `table` の上限が範囲 `range` の上端より下にあり（上端以上なら、範囲の上の入射を静止した
/// 標的と S(α, β) の移行の両方で数える）、表が範囲の下限の格子点への上方散乱（入射は下限 − TAIL kT から、移行は TAIL kT まで。
/// `kt` は物質の温度の kT）と、表の上限から下限への下方散乱を覆うこと（覆わない部分を黙って落とさない）
fn check_sab_coverage(
    name: &str,
    table: &xlc_format::library::SabTransfer,
    kt: f64,
    [low, high]: [f64; 2],
) -> Result<(), UfgError> {
    let reach = free_gas::TAIL * kt;
    let (first, floor) = (table.incident[0], table.incident[table.incident.len() - 1]);
    let (lowest, highest) = (table.transfer[0], table.transfer[table.transfer.len() - 1]);
    if floor >= high {
        return Err(UfgError::SabAboveRange {
            name: name.to_string(),
            floor,
            high,
        });
    }
    let fail = |what: &'static str, value: f64, needed: f64| UfgError::SabCoverage {
        name: name.to_string(),
        what,
        value,
        needed,
    };
    if first > low - reach {
        return Err(fail("lowest incident energy", first, low - reach));
    }
    if highest < reach {
        return Err(fail("largest transfer", highest, reach));
    }
    if lowest > low - floor {
        return Err(fail("smallest transfer", lowest, low - floor));
    }
    Ok(())
}

/// 格子点の範囲を高いエネルギーから順に解くスイープの、変わらない入力。
struct Pass<'a> {
    data: &'a [Data<'a>],
    entries: &'a [Vec<(usize, f64)>],
    cell: &'a Cylinder,
    table: &'a Option<(TrackTable, Vec<usize>)>,
    points: &'a [(f64, bool)],
    group_bounds: &'a [f64],
    ng: usize,
}

impl Pass<'_> {
    /// 格子点 `range` を順に解く（`cursors` は核種の格子のカーソル、`times` は段ごとの経過時間に足す）
    fn run(
        &self,
        sweep: &mut Sweep,
        cursors: &mut [usize],
        range: Range<usize>,
        times: &mut [Duration; 3],
    ) -> Result<(), UfgError> {
        let Self {
            data,
            entries,
            cell,
            table,
            points,
            group_bounds,
            ng,
        } = *self;
        for first in range.clone().step_by(CHUNK) {
            let clock = Instant::now();
            let chunk = first..(first + CHUNK).min(range.end);
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
                .map(|row| match table {
                    None => exact(row),
                    Some((table, first)) => {
                        let sigma: Vec<f64> =
                            first.iter().map(|&i| total(&entries[i], row)).collect();
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
                sweep.point(m, row, w, e, below)?;
                // 群の境界が不連続点なら、高い側の格子点で次の群に移る（低い側との間は幅 0）
                if sweep.group < ng && !below && e == group_bounds[sweep.group + 1] {
                    sweep.group += 1;
                }
            }
            times[0] += evaluated - clock;
            times[1] += collided - evaluated;
            times[2] += collided.elapsed();
        }
        Ok(())
    }
}

/// 降順のエネルギー `energies` を、リングの物質（`entries` の同じ並びは 1 つ）ごとの巨視的な全断面積と弾性散乱断面積が、残した点の
/// 間のレサジーについての線形な補間で、どの点でも相対差 `tolerance` 以内に再現できる点に間引く。`keep` を満たす点（群の境界と
/// 不連続点）と両端は残す。残す点は、補間の傾きの許される範囲を点ごとに狭めて決める（swinging door）。補間の起点の値は、範囲の
/// 上端と残した不連続点では減速計算と同じく低い側の極限、ほかの点ではその点の値
fn thin(
    energies: &[f64],
    data: &[Data],
    entries: &[Vec<(usize, f64)>],
    high: f64,
    tolerance: f64,
    keep: impl Fn(f64) -> bool,
) -> Result<Vec<f64>, UfgError> {
    let mut materials: Vec<&Vec<(usize, f64)>> = Vec::new();
    for ring in entries.iter().filter(|r| !r.is_empty()) {
        if !materials.contains(&ring) {
            materials.push(ring);
        }
    }
    // 物質のないセル（真空だけ）は減速計算そのものが止まるので、間引かずに返す
    if materials.is_empty() {
        return Ok(energies.to_vec());
    }
    let quantities = 2 * materials.len();
    // [核種][物質] の個数密度
    let densities: Vec<Vec<f64>> = (0..data.len())
        .map(|k| {
            materials
                .iter()
                .map(|m| m.iter().filter(|&&(j, _)| j == k).map(|&(_, n)| n).sum())
                .collect()
        })
        .collect();
    let present = |k: usize| densities[k].iter().any(|&n| n != 0.0);
    let values = macroscopic_totals(energies, data, &densities)?;
    // 減速計算と同じく、有限で負でない全断面積だけを受ける（NaN では傾きの範囲が決まらない）
    if let Some(q) = values.iter().position(|v| !(v.is_finite() && *v >= 0.0)) {
        let m = (q % quantities) / 2;
        let ring = entries.iter().position(|r| r == materials[m]).unwrap_or(0);
        return Err(CpError::Total(ring).into());
    }
    let values: Vec<&[f64]> = values.chunks_exact(quantities).collect();
    // 補間の起点（範囲の上端と、残した不連続点）の値は、減速計算と同じく低い側の極限
    let lower = |e: f64| -> Result<Vec<f64>, UfgError> {
        let mut sums = vec![0.0; quantities];
        for (k, d) in data.iter().enumerate().filter(|&(k, _)| present(k)) {
            let x = &d.pointwise().energy;
            let mut cursor = x
                .partition_point(|&v| v <= e)
                .min(x.len().saturating_sub(1));
            let row = evaluate(d, &mut cursor, e, true)?;
            for (pair, &n) in sums.as_chunks_mut::<2>().0.iter_mut().zip(&densities[k]) {
                pair[0] += n * row[0];
                pair[1] += n * row[1];
            }
        }
        Ok(sums)
    };
    let x: Vec<f64> = energies.iter().map(|&e| (high / e).ln()).collect();
    let last = energies.len().saturating_sub(1);
    let mut kept = Vec::with_capacity(energies.len() / 4);
    if energies.is_empty() {
        return Ok(kept);
    }
    kept.push(energies[0]);
    let mut a = 0;
    let mut origin = lower(energies[0])?;
    let fresh = vec![[f64::NEG_INFINITY, f64::INFINITY]; quantities];
    let mut bounds = fresh.clone();
    let mut i = 1;
    while i <= last {
        let dx = x[i] - x[a];
        if dx <= 0.0 && !(i == last || keep(energies[i])) {
            // 丸めで同じレサジーの点（格子の作り方と同じく 1 つにまとめる）
            i += 1;
            continue;
        }
        let fits = dx > 0.0
            && (0..quantities).all(|q| {
                let slope = (values[i][q] - origin[q]) / dx;
                bounds[q][0] <= slope && slope <= bounds[q][1]
            });
        if i == last || keep(energies[i]) {
            if !fits && i - 1 > a {
                kept.push(energies[i - 1]);
            }
            kept.push(energies[i]);
            a = i;
            origin = lower(energies[i])?;
            bounds.clone_from(&fresh);
            i += 1;
            continue;
        }
        if !fits {
            // 1 つ前の点まで補間でき、そこを残して始め直す
            kept.push(energies[i - 1]);
            a = i - 1;
            origin = values[a].to_vec();
            bounds.clone_from(&fresh);
            continue;
        }
        for (q, b) in bounds.iter_mut().enumerate() {
            let v = values[i][q];
            let band = tolerance * v.abs();
            b[0] = b[0].max((v - band - origin[q]) / dx);
            b[1] = b[1].min((v + band - origin[q]) / dx);
        }
        i += 1;
    }
    Ok(kept)
}

/// 降順のエネルギー `energies` の各点での、物質ごとの巨視的な全断面積と弾性散乱断面積（[点 × (物質 × 2)] の並び）。
/// `densities` は [核種][物質] の個数密度で、どの物質でも 0 の核種は足さない。点ごとに核種の順に足す
fn macroscopic_totals(
    energies: &[f64],
    data: &[Data],
    densities: &[Vec<f64>],
) -> Result<Vec<f64>, UfgError> {
    let quantities = 2 * densities.first().map_or(0, Vec::len);
    let present: Vec<(usize, &Data, &[f64])> = data
        .iter()
        .zip(densities)
        .enumerate()
        .filter(|(_, (_, n))| n.iter().any(|&n| n != 0.0))
        .map(|(k, (d, n))| (k, d, n.as_slice()))
        .collect();
    let mut values = vec![0.0; energies.len() * quantities];
    if quantities == 0 {
        return Ok(values);
    }
    // 点の塊ごとに並列に足し、和集合の格子の全体を核種ごとに流さない（16 本のピンが同時に流すとメモリの転送で詰まる。
    // ADR-0198）。塊の始めのカーソルは、その点以下の最後の格子点の次で、`bracket` はそこから下るので、上の点から下って
    // きたカーソルと同じ区間を返す。エラーは、核種の順、点の順に最初のもの（塊の並列の順によらない）
    let failures: Vec<(usize, usize, UfgError)> = values
        .par_chunks_mut(THIN_POINTS * quantities)
        .zip(energies.par_chunks(THIN_POINTS))
        .enumerate()
        .filter_map(|(b, (block, energies))| {
            for &(k, d, densities) in &present {
                let x = &d.pointwise().energy;
                let mut cursor = x.partition_point(|&v| v <= energies[0]);
                for (i, (row, &e)) in block.chunks_exact_mut(quantities).zip(energies).enumerate() {
                    let [total, elastic] = match evaluate_total_elastic(d, &mut cursor, e) {
                        Ok(v) => v,
                        Err(error) => return Some((k, b * THIN_POINTS + i, error)),
                    };
                    for (pair, &n) in row.as_chunks_mut::<2>().0.iter_mut().zip(densities) {
                        pair[0] += n * total;
                        pair[1] += n * elastic;
                    }
                }
            }
            None
        })
        .collect();
    match failures.into_iter().min_by_key(|&(k, i, _)| (k, i)) {
        Some((_, _, error)) => Err(error),
        None => Ok(values),
    }
}

/// 核種 `d` の群定数: ポイントワイズの断面積を、リングの中性子束 `flux`（減速計算の格子点の値。点の間で u について線形）で
/// 縮約する。ポイントワイズで解いた核種と同じく σF をレサジーの台形則で積分し、点は減速計算の格子点と核種自身の格子の点の和
/// （格子点の間の核種の点では中性子束を内挿する）にする。間引かない格子は核種の点を含むので、減速計算の中の積分と同じ値になる
/// （ADR-0185）。`fluxes` と `group_fluxes` は核種を含むリングごとの格子点の F と群ごとの ∫F du で、群定数はその順に返す
fn collapse(
    d: &Data,
    u: &[f64],
    points: &[(f64, bool)],
    fluxes: &[&[f64]],
    high: f64,
    group_bounds: &[f64],
    group_fluxes: &[&[f64]],
) -> Result<Vec<Shielded>, UfgError> {
    let p = d.pointwise();
    let ng = group_bounds.len() - 1;
    let rings = fluxes.len();
    // [リング][群] の反応率
    let mut rates = vec![vec![[0.0; COLUMNS]; ng]; rings];
    let mut cursor = p.energy.len();
    // 核種の格子の、まだ使っていない最も高いエネルギーの点の次（昇順の配列を降順に進める）
    let mut next = p.energy.len();
    let mut group = 0;
    // 1 つ前の点の (u, 断面積) と、リングごとの F
    let mut previous: Option<(f64, Row)> = None;
    let mut last_flux = vec![0.0; rings];
    let mut node_flux = vec![0.0; rings];
    let mut here_flux = vec![0.0; rings];
    let mut add =
        |group: usize, (ua, va): (f64, &Row), fa: &[f64], (ub, vb): (f64, &Row), fb: &[f64]| {
            let du = ub - ua;
            for (rate, (fa, fb)) in rates.iter_mut().zip(fa.iter().zip(fb)) {
                for (r, (a, b)) in rate[group].iter_mut().zip(va.iter().zip(vb)) {
                    *r += 0.5 * du * (fa * a + fb * b);
                }
            }
        };
    for (m, &(e, below)) in points.iter().enumerate() {
        let here = (u[m], evaluate(d, &mut cursor, e, below)?);
        for (f, flux) in here_flux.iter_mut().zip(fluxes) {
            *f = flux[m];
        }
        if let Some(mut last) = previous {
            let upper = points[m - 1].0;
            while next > 0 && p.energy[next - 1] >= upper {
                next -= 1;
            }
            while next > 0 && p.energy[next - 1] > e {
                next -= 1;
                let x = (high / p.energy[next]).ln();
                let t = (x - u[m - 1]) / (here.0 - u[m - 1]);
                for (f, flux) in node_flux.iter_mut().zip(fluxes) {
                    *f = flux[m - 1] + (flux[m] - flux[m - 1]) * t;
                }
                let row = [
                    p.total[next],
                    p.elastic[next],
                    p.fission[next],
                    p.nu_fission[next],
                    p.capture[next],
                ];
                add(group, (last.0, &last.1), &last_flux, (x, &row), &node_flux);
                last = (x, row);
                last_flux.clone_from(&node_flux);
            }
            add(
                group,
                (last.0, &last.1),
                &last_flux,
                (here.0, &here.1),
                &here_flux,
            );
        }
        previous = Some(here);
        last_flux.clone_from(&here_flux);
        // 群の境界では、高い側の格子点で次の群に移る（減速計算と同じ）
        if group < ng && !below && e == group_bounds[group + 1] {
            group += 1;
        }
    }
    Ok(rates
        .iter()
        .zip(group_fluxes)
        .map(|(rate, flux)| {
            let column =
                |x: usize| -> Vec<f64> { rate.iter().zip(*flux).map(|(r, f)| r[x] / f).collect() };
            Shielded {
                total: column(0),
                elastic: column(1),
                fission: column(2),
                nu_fission: column(3),
                capture: column(4),
            }
        })
        .collect())
}

/// 核種ごとの格子の点の数を debug で出す。`own` は範囲の中の自分の点の数、`added` は、データの順に足したときに格子の和に
/// 加わる（それまでの核種と群の境界にない）エネルギーの数
fn log_nuclide_points(data: &[Data], group_bounds: &[f64], inside: impl Fn(f64) -> bool) {
    let mut seen: HashSet<u64> = group_bounds.iter().map(|e| e.to_bits()).collect();
    for d in data {
        let own: Vec<f64> = d
            .pointwise()
            .energy
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
    let (j, t) = bracket(d, cursor, e, below)?;
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

/// 核種 `d` の全断面積と弾性散乱断面積の、エネルギー `e` での lin-lin の補間（[`evaluate`] の高い側の極限と同じ値）
fn evaluate_total_elastic(d: &Data, cursor: &mut usize, e: f64) -> Result<[f64; 2], UfgError> {
    let (j, t) = bracket(d, cursor, e, false)?;
    let p = d.pointwise();
    let lerp = |y: &[f64]| y[j] + t * (y[j + 1] - y[j]);
    Ok([lerp(&p.total), lerp(&p.elastic)])
}

/// [`evaluate`] の補間の区間の始まりの点と、区間の中の位置
fn bracket(d: &Data, cursor: &mut usize, e: f64, below: bool) -> Result<(usize, f64), UfgError> {
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
    Ok((j, t))
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
    /// 反応率を格子点で積分する（間引いた格子では核種の格子で縮約するので積分しない）
    tally: bool,
    /// 格子点ごとに再利用する作業用の配列: [リング][項] の窓の積分の既知の部分
    known: Vec<Vec<f64>>,
    /// 格子点ごとに再利用する作業用の配列: [リング][項] のその点の F_i の係数
    own: Vec<Vec<f64>>,
    /// [核種] の、窓の下端の区間の積分の範囲の上の部分と、格子点の値の重み
    slabs: Vec<(f64, Vec<(usize, f64)>)>,
    /// 自由ガスの核の散乱源（[`Scattering::FreeGas`] で、自由ガスの核種のあるとき）
    thermal: Option<Thermal>,
    /// 格子点ごとに再利用する作業用の配列: [リング] の自由ガスの核の散乱源の既知の部分と、その点の F_i の係数
    extra: (Vec<f64>, Vec<f64>),
}

/// 掃き直しの始めに戻す、スイープの状態
struct Saved {
    windows: Vec<Vec<f64>>,
    lower: Vec<usize>,
    previous: Vec<Row>,
    group: usize,
    group_flux: Vec<Vec<f64>>,
    rates: Vec<Vec<Vec<Row>>>,
    /// 解いた格子点の数
    len: usize,
}

impl Sweep<'_> {
    fn save(&self) -> Saved {
        Saved {
            windows: self.windows.clone(),
            lower: self.lower.clone(),
            previous: self.previous.clone(),
            group: self.group,
            group_flux: self.group_flux.clone(),
            rates: self.rates.clone(),
            len: self.flux.first().map_or(0, Vec::len),
        }
    }

    fn restore(&mut self, saved: &Saved) {
        self.windows.clone_from(&saved.windows);
        self.lower.clone_from(&saved.lower);
        self.previous.clone_from(&saved.previous);
        self.group = saved.group;
        self.group_flux.clone_from(&saved.group_flux);
        self.rates.clone_from(&saved.rates);
        for f in &mut self.flux {
            f.truncate(saved.len);
        }
        for e in &mut self.elastic {
            e.truncate(saved.len);
        }
    }

    /// 格子点 m の中性子束を求め、窓の積分と群の積分を進める。`values` はその点のポイントワイズの断面積、`tracks` は飛程、
    /// `below` は不連続点の低い側の極限か。
    fn point(
        &mut self,
        m: usize,
        values: Vec<Row>,
        tracks: &[Vec<f64>],
        energy: f64,
        below: bool,
    ) -> Result<(), UfgError> {
        let n = self.volumes.len();
        for ((column, v), d) in self.elastic.iter_mut().zip(&values).zip(self.data) {
            // 自由ガスの核の、静止した標的の窓の部分は σs/N（free_gas の R_static/N）。S(α, β) の上限より下は 0
            match d.gas {
                Some(gas) => column.push(v[1] * gas.window_scale(energy, below)),
                None => column.push(v[1]),
            }
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
        let mut k: Vec<f64> = (0..n).map(|i| source(i, &known[i])).collect();
        let mut c: Vec<f64> = (0..n).map(|i| source(i, &own[i])).collect();
        if let Some(thermal) = &mut self.thermal {
            let (extra_known, extra_own) = &mut self.extra;
            extra_known.fill(0.0);
            extra_own.fill(0.0);
            thermal.source(
                m,
                self.u,
                self.above.high,
                energy,
                &self.flux,
                extra_known,
                extra_own,
            )?;
            for i in 0..n {
                k[i] += extra_known[i];
                c[i] += extra_own[i];
            }
        }
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
        if let Some(thermal) = &mut self.thermal {
            thermal.record(m, self.u, &self.flux);
        }
        // 群の積分（レサジーの台形則）
        if m > 0 {
            let du = self.u[m] - self.u[m - 1];
            let g = self.group;
            for (i, ring) in self.entries.iter().enumerate() {
                let (f0, f1) = (self.flux[i][m - 1], f[i]);
                self.group_flux[i][g] += 0.5 * du * (f0 + f1);
                if !self.tally {
                    continue;
                }
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
