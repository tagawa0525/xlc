//! 格子計算の燃焼（ADR-0173）。
//!
//! 入力の `[depletion]` に挙げた燃料の物質を、自己遮蔽の対称操作で移り合うセルの組ごとに 1 つの燃焼の領域にし、全核種の
//! チェーン（ADR-0164）の数密度を領域ごとに持つ。時間は CE/CM で進める:
//!
//! ```text
//! n½ = exp(h/2 A(n₀)) n₀,  n₁ = exp(h A(n½)) n₀
//! ```
//!
//! A(n) は組成 n で自己遮蔽と MOC を解き、一定の出力密度に規格化した反応率の燃焼行列。行列指数は CRAM48（ADR-0161）。
//! 反応率は、領域ごと・核種ごとに、自己遮蔽した微視的断面積（SPH 因子を掛けない値）と物理的な中性子束の積を群で足し、
//! 領域の体積で平均したもの。出力は Σ V N σf φ Q（Q はチェーンの MF3 MT18 の QI）が、初めの重金属の質量 × 出力密度に
//! 等しくなるよう規格化する。2 次元の格子なので、体積と質量は高さ 1 cm あたり。
//!
//! 燃焼の履歴は、燃焼点ごとに、その点の組成で解いた輸送の k∞、ピン出力（`pin-powers` と同じ定義。ADR-0181）、少数群定数を持つ。

use std::collections::HashMap;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use anyhow::{Context, anyhow, bail};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use xlc_depletion::chain::{Chain, ChainError, Fission, NuclideRates, Reaction};
use xlc_depletion::{CramOrder, deplete};
use xlc_format::fewgroup::FewGroup;
use xlc_format::library::{Library, Shielded};
use xlc_format::msgpack;
use xlc_format::table::Table;

use crate::fewgroup::{homogenize, pin_shapes};
use crate::geometry::Lattice;
use crate::input::{LatticeInput, MaterialInput};
use crate::material::Material;
use crate::moc::{Options, Start, available_order, solve_symmetric_from};
use crate::power::{cell_fission_rates, pin_rows};
use crate::shielding::{FuelShielding, ShieldingReuse, shield_fuel_reusing, shielded_nuclide};
use crate::symmetry::Chosen;
use crate::tracking::Tracking;

/// 入力の `[depletion]`（ADR-0173）。
#[derive(Debug, Clone, PartialEq)]
pub struct DepletionInput {
    /// 燃焼チェーンのファイル（入力のファイルのあるディレクトリから）
    pub chain: PathBuf,
    /// 初めの重金属の質量あたりの出力密度 (W/g)
    pub power: f64,
    /// 燃焼度の点 (MWd/kg、初めの重金属)。0 より大きく増える。各点がステップの終わり
    pub burnup: Vec<f64>,
    /// 燃焼させる物質（核種で書いた物質）
    pub materials: Vec<String>,
    /// 燃焼の履歴のファイル（入力のファイルのあるディレクトリから）
    pub output: PathBuf,
    /// 燃焼点ごとに格子全体の少数群定数を作るなら、その条件
    pub few_group: Option<DepletionFewGroup>,
    /// 自己遮蔽を解く輸送（ADR-0197）
    pub shielding: ShieldingFrequency,
}

/// 燃焼の自己遮蔽を解く輸送（ADR-0197）
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShieldingFrequency {
    /// ステップの始まりと中点の輸送のどちらでも、その組成で解き直す（ADR-0173 の案 F1）
    EveryTransport,
    /// ステップの始まり（と最後の点）の輸送でだけ解く。中点の輸送は、始まりの自己遮蔽の微視的断面積と SPH 因子のまま、
    /// 燃焼の領域の物質を中点の数密度で作り直す（ADR-0173 の案 F2）
    StepStart,
}

/// 燃焼点ごとの格子全体の少数群定数の条件
#[derive(Debug, Clone, PartialEq)]
pub struct DepletionFewGroup {
    /// 群境界 (eV)
    pub bounds: Vec<f64>,
    /// ピン出力の再構成の形状関数と角の不連続因子も作る（ADR-0096）
    pub pins: bool,
}

/// 行列指数の次数（OpenMC の既定の CRAM48 と同じ。ADR-0173）
pub const CRAM_ORDER: CramOrder = CramOrder::Order48;

/// 中性子の質量 (u)。CODATA 2018
const NEUTRON_MASS: f64 = 1.008_664_915_95;
/// 原子質量単位 (g)。CODATA 2018
const ATOMIC_MASS_UNIT: f64 = 1.660_539_066_60e-24;
/// 1 eV (J)
const ELECTRON_VOLT: f64 = 1.602_176_634e-19;
/// 1 barn (cm²)。数密度 1/(barn cm) を 1/cm³ にするのにも使う
const BARN: f64 = 1.0e-24;
/// 1 日 (s)
const DAY: f64 = 86_400.0;
/// 重金属とする原子番号の下限
pub(crate) const HEAVY_METAL_Z: u32 = 90;
/// CRAM の近似の誤差で説明できる負の数密度の、領域の最大の数密度に対する大きさ（近似の誤差は初めの最大の 1e-15 倍の
/// 水準。ADR-0161）。これより大きな負の値はエラーにする（ADR-0173）
const NEGATIVE_TOLERANCE: f64 = 1.0e-12;

/// 燃焼度の点 `burnup` (MWd/kg) までの各ステップの長さ (s)。出力密度 `power` (W/g) で、1 MWd/kg = 1000 W d/g
pub fn step_times(power: f64, burnup: &[f64]) -> Vec<f64> {
    let mut previous = 0.0;
    burnup
        .iter()
        .map(|&b| {
            let step = (b - previous) * 1.0e3 / power * DAY;
            previous = b;
            step
        })
        .collect()
}

/// 核分裂の収率の熱の点 (eV)。境界 `yield_cutoff` より下の核分裂に、この点か、なければ境界より下で最も近い点を使う（ADR-0183）
const YIELD_THERMAL: f64 = 0.0253;
/// 核分裂の収率の高速の点 (eV)。境界 `yield_cutoff` より上の核分裂に、この点か、なければ境界以上で最も近い点を使う（ADR-0183）
const YIELD_FAST: f64 = 5.0e5;

/// 収率の熱と高速の点を分ける境界 (eV)。2 点の幾何平均 √(0.0253 · 5e5) ≈ 112.47 eV で、参照解の OpenMC に与える
/// `cutoff` と同じ（OpenMC 0.15.2 の既定の 112.0 eV ではない。ADR-0183）
fn yield_cutoff() -> f64 {
    (YIELD_THERMAL * YIELD_FAST).sqrt()
}

/// 核分裂 `fission` の収率のエネルギーの点ごとの重み（和は 1。ADR-0183）。群の核分裂率 `rates`（群境界 `group_bounds` は
/// 降順の eV）を `yield_cutoff` の下と上に分け、境界を含む群はレサジーの幅の比で分ける（OpenMC は連続エネルギーで数える）。
/// 下を熱の点、上を高速の点に置く。点の選び方は OpenMC 0.15.2 の `FissionYieldCutoffHelper.__init__` と同じ: 熱と高速の点が
/// ちょうどあればそれを使う。どちらかがなければ、境界が最初の点以下なら最初の点、最後の点以上なら最後の点だけを使い、そうで
/// なければ境界より下で熱に、境界以上で高速に最も近い点を使う。熱と高速の 2 点を使う核種で核分裂率がすべて 0 なら、最初の点に 1（`Chain::zero_rates` と同じ）
pub fn yield_weights(group_bounds: &[f64], rates: &[f64], fission: &Fission) -> Vec<f64> {
    let energies = &fission.energies;
    let cutoff = yield_cutoff();
    let mut weights = vec![0.0; energies.len()];
    let exact = |e: f64| energies.iter().position(|&p| p == e);
    let (thermal, fast) = match (exact(YIELD_THERMAL), exact(YIELD_FAST)) {
        (Some(t), Some(f)) => (t, f),
        (thermal, fast) => {
            if cutoff <= energies[0] {
                weights[0] = 1.0;
                return weights;
            }
            let last = energies.len() - 1;
            if cutoff >= energies[last] {
                weights[last] = 1.0;
                return weights;
            }
            let split = energies.partition_point(|&e| e < cutoff);
            // 境界は最初と最後の点の間にあるので、どちらの範囲も空でない。同じ近さなら前の点（Python の min と同じ）
            let nearest = |from: usize, to: usize, target: f64| {
                (from..to).fold(from, |best, k| {
                    if (energies[k] - target).abs() < (energies[best] - target).abs() {
                        k
                    } else {
                        best
                    }
                })
            };
            (
                thermal.unwrap_or_else(|| nearest(0, split, YIELD_THERMAL)),
                fast.unwrap_or_else(|| nearest(split, energies.len(), YIELD_FAST)),
            )
        }
    };
    let (below, above) = split_at(group_bounds, rates, cutoff);
    let total = below + above;
    if total <= 0.0 {
        weights[0] = 1.0;
        return weights;
    }
    weights[thermal] += below / total;
    weights[fast] += above / total;
    weights
}

/// 群の核分裂率 `rates`（群境界 `group_bounds` は降順の eV）の、`energy` より下と上の和。`energy` を含む群は
/// レサジーの幅の比で分ける
fn split_at(group_bounds: &[f64], rates: &[f64], energy: f64) -> (f64, f64) {
    let (mut below, mut above) = (0.0, 0.0);
    for (g, &rate) in rates.iter().enumerate() {
        let (upper, lower) = (group_bounds[g], group_bounds[g + 1]);
        if energy >= upper {
            below += rate;
        } else if energy <= lower {
            above += rate;
        } else {
            let part = (energy / lower).ln() / (upper / lower).ln();
            below += rate * part;
            above += rate * (1.0 - part);
        }
    }
    (below, above)
}

/// 反応ごとの群の断面積の出どころ（ADR-0173）
pub trait ReactionSource: Sync {
    /// 核種 `nuclide` の温度 `temperature` (K) のチェーンの反応 `reaction` の、経路（`reaction.channels`）ごとの群の微視的断面積
    /// (barn)。`shielded` はその領域で自己遮蔽した核種の群定数（SPH 因子を掛けない値）
    fn cross_sections(
        &self,
        nuclide: &str,
        temperature: f64,
        reaction: &Reaction,
        shielded: &Shielded,
    ) -> Result<Vec<Vec<f64>>, String>;
}

/// 多群ライブラリのチェーンの反応（ADR-0172）から引く。捕獲（MT102）は自己遮蔽した捕獲を、経路の無限希釈の断面積の比で群ごとに
/// 分ける（ADR-0172 の案 T2）。しきい値の反応は無限希釈の値のまま
pub struct LibraryReactions<'a>(pub &'a Library);

/// 捕獲の MT
const CAPTURE: u32 = 102;

impl ReactionSource for LibraryReactions<'_> {
    fn cross_sections(
        &self,
        nuclide: &str,
        temperature: f64,
        reaction: &Reaction,
        shielded: &Shielded,
    ) -> Result<Vec<Vec<f64>>, String> {
        let (_, data) = self
            .0
            .lookup(nuclide, temperature)
            .map_err(|e| e.to_string())?;
        let mt = reaction.mt;
        let group = data
            .reactions
            .iter()
            .find(|r| r.mt == mt)
            .ok_or_else(|| format!("{nuclide}: the library has no chain reaction MT{mt}"))?;
        // 経路の並びと状態はチェーンと同じ（ADR-0172）
        let same = group.channels.len() == reaction.channels.len()
            && group.channels.iter().zip(&reaction.channels).all(|(g, c)| {
                match (g.final_state, c.final_state) {
                    (None, None) => true,
                    (Some(g), Some(c)) => g.mf == c.mf && g.lfs == c.lfs,
                    _ => false,
                }
            });
        if !same {
            return Err(format!(
                "{nuclide}: the channels of MT{mt} in the library differ from the chain"
            ));
        }
        let dilute: Vec<Vec<f64>> = group
            .channels
            .iter()
            .map(|c| c.cross_section.clone())
            .collect();
        if mt != CAPTURE {
            return Ok(dilute);
        }
        let groups = shielded.capture.len();
        let mut channels = vec![vec![0.0; groups]; dilute.len()];
        for g in 0..groups {
            let sum: f64 = dilute.iter().map(|c| c[g]).sum();
            if sum > 0.0 {
                for (out, c) in channels.iter_mut().zip(&dilute) {
                    out[g] = shielded.capture[g] * c[g] / sum;
                }
            } else if shielded.capture[g] != 0.0 {
                return Err(format!(
                    "{nuclide}: shielded capture {} in group {g} without dilute channels",
                    shielded.capture[g]
                ));
            }
        }
        Ok(channels)
    }
}

/// 1 回の輸送で求めた、領域ごとの反応率（規格化したもの。[領域][チェーンの核種]）と、点に残す値
pub struct Evaluation<T> {
    pub rates: Vec<Vec<NuclideRates>>,
    pub extra: T,
}

/// 燃焼点: 時間 (s)、領域ごとのチェーンの数密度 (1/(barn cm))、その組成で解いた輸送の値
pub struct Point<T> {
    pub time: f64,
    pub densities: Vec<Vec<f64>>,
    pub extra: T,
}

/// 燃焼行列を組み立てる時間と、行列指数の時間
#[derive(Debug, Clone, Copy, Default)]
pub struct DepletionTimes {
    pub matrix: Duration,
    pub cram: Duration,
}

/// 領域ごとに、反応率 `rates` の燃焼行列で `n` を時間 `t` (s) だけ進める
fn advance(
    chain: &Chain,
    rates: &[Vec<NuclideRates>],
    n: &[Vec<f64>],
    t: f64,
    times: &mut DepletionTimes,
) -> Result<Vec<Vec<f64>>, ChainError> {
    let start = Instant::now();
    let matrices = rates
        .par_iter()
        .map(|r| chain.matrix(r))
        .collect::<Result<Vec<_>, _>>()?;
    times.matrix += start.elapsed();
    let start = Instant::now();
    let next = matrices
        .par_iter()
        .zip(n)
        .map(|(a, n)| deplete(a, n, t, CRAM_ORDER).map_err(ChainError::from))
        .collect::<Result<Vec<_>, _>>()?;
    times.cram += start.elapsed();
    Ok(next)
}

/// CE/CM（ADR-0173）で、初めの組成 `initial`（[領域][チェーンの核種]）からステップ `steps` (s) を進める。`solve` は組成から
/// 規格化した反応率を求める輸送で、ステップごとに始まりと中点の 2 回、最後の点で 1 回呼ぶ。2 つ目の引数は、その輸送の値
/// （`Evaluation::extra`）を点に残すか（中点の輸送は残さない）。点は始まりと各ステップの終わりで、点の値はその点の組成で解いた
/// 輸送のもの
pub fn ce_cm<T, E: From<ChainError>>(
    chain: &Chain,
    initial: Vec<Vec<f64>>,
    steps: &[f64],
    mut solve: impl FnMut(&[Vec<f64>], bool) -> Result<Evaluation<T>, E>,
) -> Result<(Vec<Point<T>>, DepletionTimes), E> {
    let mut times = DepletionTimes::default();
    let mut points = Vec::with_capacity(steps.len() + 1);
    let mut n = initial;
    let mut time = 0.0;
    for &h in steps {
        let start = solve(&n, true)?;
        let middle = advance(chain, &start.rates, &n, 0.5 * h, &mut times)?;
        let corrector = solve(&middle, false)?;
        let next = advance(chain, &corrector.rates, &n, h, &mut times)?;
        points.push(Point {
            time,
            densities: std::mem::replace(&mut n, next),
            extra: start.extra,
        });
        time += h;
    }
    let last = solve(&n, true)?;
    points.push(Point {
        time,
        densities: n,
        extra: last.extra,
    });
    Ok((points, times))
}

/// 燃焼の領域の 1 つの MOC の領域: 体積 (cm²)、物理的な中性子束 [群]、輸送する核種ごとの自己遮蔽した群定数
pub struct Member<'a> {
    pub volume: f64,
    pub flux: &'a [f64],
    pub cross_sections: &'a [Shielded],
}

/// 燃焼の領域の、中性子束の大きさを MOC の規格化のままとしたときの 1 群の反応率 (1/s、1 原子あたり)。
/// `composition` は輸送する核種のチェーンの添字（`Member::cross_sections` の順）。輸送しない核種は 0
pub fn region_rates(
    chain: &Chain,
    composition: &[usize],
    temperature: f64,
    members: &[Member],
    group_bounds: &[f64],
    reactions: &dyn ReactionSource,
) -> anyhow::Result<Vec<NuclideRates>> {
    let mut rates = chain.zero_rates();
    let volume: f64 = members.iter().map(|m| m.volume).sum();
    let groups = group_bounds.len() - 1;
    let collapse = |sigma: &[f64], flux: &[f64]| -> f64 {
        sigma.iter().zip(flux).map(|(s, f)| s * f).sum::<f64>()
    };
    for (k, &i) in composition.iter().enumerate() {
        let nuclide = &chain.nuclides[i];
        let mut group_fission = vec![0.0; groups];
        for m in members {
            let shielded = &m.cross_sections[k];
            for (g, f) in group_fission.iter_mut().enumerate() {
                *f += m.volume * shielded.fission[g] * m.flux[g];
            }
            for (r, reaction) in nuclide.reactions.iter().enumerate() {
                let channels = reactions
                    .cross_sections(&nuclide.name, temperature, reaction, shielded)
                    .map_err(|e| anyhow!(e))?;
                if channels.len() != reaction.channels.len() {
                    bail!(
                        "{}: {} channel cross sections for MT{} with {} channels",
                        nuclide.name,
                        channels.len(),
                        reaction.mt,
                        reaction.channels.len()
                    );
                }
                for (c, sigma) in channels.iter().enumerate() {
                    if sigma.len() != groups {
                        bail!(
                            "{}: {} groups in a channel of MT{} for {groups} groups",
                            nuclide.name,
                            sigma.len(),
                            reaction.mt
                        );
                    }
                    rates[i].reactions[r][c] += m.volume * collapse(sigma, m.flux) * BARN / volume;
                }
            }
        }
        let fission: f64 = group_fission.iter().sum();
        match &nuclide.fission {
            Some(f) => {
                rates[i].fission = fission * BARN / volume;
                rates[i].yield_weights = yield_weights(group_bounds, &group_fission, f);
            }
            None if fission > 0.0 => bail!(
                "{} fissions in the library but has no fission in the chain",
                nuclide.name
            ),
            None => {}
        }
    }
    Ok(rates)
}

/// ライブラリが核種 `name` を温度 `temperature` (K) そのもので持つか。熱中性子散乱の核種は温度を内挿して引けるが、輸送する
/// 核種は参照解（OpenMC の燃料の温度のデータのある核種）と揃えるので、内挿は数えない（ADR-0173）
pub fn has_data_at(library: &Library, name: &str, temperature: f64) -> bool {
    library
        .nuclides
        .iter()
        .find(|n| n.name == name)
        .is_some_and(|n| n.temperatures.iter().any(|t| t.temperature == temperature))
}

/// 燃焼の領域: 物質と、組の代表のセル、MOC の領域（組のすべてのセルの分）と、その体積の和 (cm²)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DepletionRegion {
    pub material: String,
    pub cell: usize,
    pub regions: Vec<usize>,
    pub volume: f64,
}

/// 燃焼の領域を自己遮蔽と MOC で解くときの物質の名前
pub fn transport_name(material: &str, cell: usize) -> String {
    format!("{material}@{cell}")
}

/// 解く物質の名前（`<物質>@<セル>`）が、入力の物質 `materials` の名前と重ならないこと
pub fn check_region_names(
    materials: &[MaterialInput],
    regions: &[DepletionRegion],
) -> anyhow::Result<()> {
    for d in regions {
        let name = transport_name(&d.material, d.cell);
        if materials.iter().any(|m| m.name == name) {
            bail!("material `{name}` clashes with the name of a depletion region");
        }
    }
    Ok(())
}

/// 領域ごとの解く物質の名前。燃焼させる物質 `depleted` の領域は、組の代表のセル `representatives`（[セル]）で
/// `<物質>@<セル>` にする
pub fn solved_names(
    lattice: &Lattice,
    region_materials: &[String],
    representatives: &[usize],
    depleted: &[String],
) -> Vec<String> {
    region_materials
        .iter()
        .enumerate()
        .map(|(r, name)| {
            if depleted.contains(name) {
                transport_name(name, representatives[lattice.region_cell(r)])
            } else {
                name.clone()
            }
        })
        .collect()
}

/// 燃焼の領域 `regions` の数密度 `densities`（[領域][核種]）の重金属の質量 (g、高さ 1 cm あたり)。`heavy(i)` は核種 i が重金属
/// （Z ≥ 90）ならその AWR、ほかは `None`。燃焼の出力の規格化と燃焼度の定義（ADR-0173）に使い、履歴からの分岐のテーブルにも
/// 同じ式で書く（ADR-0190）。重金属がなければエラー
pub(crate) fn heavy_metal_mass(
    regions: &[DepletionRegion],
    densities: &[Vec<f64>],
    heavy: impl Fn(usize) -> anyhow::Result<Option<f64>>,
) -> anyhow::Result<f64> {
    let mut mass = 0.0;
    for (d, n) in regions.iter().zip(densities) {
        for (i, &density) in n.iter().enumerate() {
            if density > 0.0
                && let Some(awr) = heavy(i)?
            {
                mass += d.volume * density / BARN * awr * NEUTRON_MASS * ATOMIC_MASS_UNIT;
            }
        }
    }
    if mass <= 0.0 {
        bail!("the depleted materials have no heavy metal (Z >= {HEAVY_METAL_Z})");
    }
    Ok(mass)
}

/// 燃焼の領域ごとの数密度 `n`（[領域][核種]）の負の値が、近似の誤差の水準であること（ADR-0173）。`name` は核種の名前
pub fn check_negative(
    regions: &[DepletionRegion],
    n: &[Vec<f64>],
    name: impl Fn(usize) -> String,
) -> anyhow::Result<()> {
    for (d, n) in regions.iter().zip(n) {
        let largest = n.iter().fold(0.0_f64, |a, &b| a.max(b));
        if let Some((i, v)) = n
            .iter()
            .enumerate()
            .find(|&(_, &v)| v < -NEGATIVE_TOLERANCE * largest)
        {
            bail!(
                "{} in material `{}` of cell {} is {v}, below the CRAM error",
                name(i),
                d.material,
                d.cell
            );
        }
    }
    Ok(())
}

/// 燃焼の領域 `d` を解く物質: 温度 `temperature`、輸送する核種 `transported`（番号）の数密度 `n`（負の値は 0。燃焼の
/// 輸送と履歴からの分岐で同じ。ADR-0173、ADR-0187）。`name` は核種の名前
pub fn region_material(
    d: &DepletionRegion,
    temperature: f64,
    n: &[f64],
    transported: &[usize],
    name: impl Fn(usize) -> String,
) -> MaterialInput {
    MaterialInput {
        name: transport_name(&d.material, d.cell),
        temperature,
        composition: transported
            .iter()
            .map(|&i| (name(i), n[i].max(0.0)))
            .collect(),
    }
}

/// 燃焼させる物質 `depleted` の領域を、セルの組の代表 `representatives`（[セル]）と物質ごとにまとめる（現れる順）
pub fn depletion_regions(
    lattice: &Lattice,
    region_materials: &[String],
    representatives: &[usize],
    depleted: &[String],
) -> Vec<DepletionRegion> {
    let volumes = lattice.volumes();
    let mut found: Vec<DepletionRegion> = Vec::new();
    for (r, name) in region_materials.iter().enumerate() {
        if !depleted.contains(name) {
            continue;
        }
        let cell = representatives[lattice.region_cell(r)];
        match found
            .iter_mut()
            .find(|d| d.cell == cell && &d.material == name)
        {
            Some(d) => {
                d.regions.push(r);
                d.volume += volumes[r];
            }
            None => found.push(DepletionRegion {
                material: name.clone(),
                cell,
                regions: vec![r],
                volume: volumes[r],
            }),
        }
    }
    found
}

/// 燃焼の履歴の 1 点
#[derive(Debug, Clone)]
pub struct HistoryPoint {
    /// 時間 (s)
    pub time: f64,
    /// 燃焼度 (MWd/kg、初めの重金属)
    pub burnup: f64,
    pub k: f64,
    /// [燃焼の領域][チェーンの核種] の数密度 (1/(barn cm))
    pub densities: Vec<Vec<f64>>,
    /// ピン出力（`pin-powers` と同じ。[行（上から）][列（左から）]、核分裂のあるセルの平均が 1。ADR-0181）
    pub pin_powers: Vec<Vec<f64>>,
    /// `[depletion]` の `few-group` があるときの、格子全体の少数群定数
    pub few_group: Option<FewGroup>,
    /// この点の組成で解いた輸送の核種の数
    pub transported: usize,
}

/// 輸送の各段の時間の和
#[derive(Debug, Clone, Copy, Default)]
pub struct TransportTimes {
    pub shielding: Duration,
    pub moc: Duration,
    /// 反応率と規格化（少数群定数を含む）
    pub rates: Duration,
}

/// 燃焼の履歴（ADR-0173）
#[derive(Debug, Clone)]
pub struct History {
    /// チェーンの核種の名前（数密度の列の順）
    pub nuclides: Vec<String>,
    pub regions: Vec<DepletionRegion>,
    /// 輸送できるチェーンの核種（ライブラリがすべての燃焼の領域の温度で持つもの。チェーンの順。ADR-0173）
    pub transportable: Vec<String>,
    pub points: Vec<HistoryPoint>,
    pub transport: TransportTimes,
    pub depletion: DepletionTimes,
    /// 輸送を解いた回数
    pub solves: usize,
    /// 自己遮蔽を解いた回数（ADR-0197）
    pub shieldings: usize,
}

const FORMAT: &str = "xlc-depletion-history";

#[derive(Serialize)]
struct HistoryFile<'a> {
    format: &'static str,
    nuclides: &'a [String],
    regions: &'a [DepletionRegion],
    transportable: &'a [String],
    points: Vec<PointFile<'a>>,
    seconds: SecondsFile,
    solves: usize,
    shieldings: usize,
}

/// 時間の内訳 (s)
#[derive(Serialize)]
struct SecondsFile {
    shielding: f64,
    moc: f64,
    rates: f64,
    matrix: f64,
    cram: f64,
}

#[derive(Serialize)]
struct PointFile<'a> {
    time: f64,
    burnup: f64,
    k_inf: f64,
    densities: &'a [Vec<f64>],
    pin_powers: &'a [Vec<f64>],
    /// `xlc-table` の状態が 1 つのテーブル（ADR-0097）の map か nil
    few_group: Option<rmpv::Value>,
}

/// 読む燃焼の履歴のファイル（書く形 `HistoryFile` と同じ項目）
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct HistoryFileIn {
    format: String,
    nuclides: Vec<String>,
    regions: Vec<DepletionRegion>,
    transportable: Vec<String>,
    points: Vec<PointFileIn>,
    seconds: SecondsFileIn,
    solves: usize,
    shieldings: usize,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SecondsFileIn {
    shielding: f64,
    moc: f64,
    rates: f64,
    matrix: f64,
    cram: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct PointFileIn {
    time: f64,
    burnup: f64,
    k_inf: f64,
    densities: Vec<Vec<f64>>,
    pin_powers: Vec<Vec<f64>>,
    /// nil を書いても、キーは省けない（`deserialize_with` を付けると、serde は省いた `Option` を `None` にしない）
    #[serde(deserialize_with = "xlc_format::msgpack::nil_or")]
    few_group: Option<rmpv::Value>,
}

/// 時間の内訳 (s) を `Duration` にする（負や有限でない値はエラー）
fn seconds(name: &str, s: f64) -> anyhow::Result<Duration> {
    Duration::try_from_secs_f64(s).map_err(|e| anyhow!("seconds.{name} {s}: {e}"))
}

impl History {
    /// `to_msgpack` の形を読む（ADR-0187）。核種の名前、燃焼度が 0 から増えること、各点の数密度が [燃焼の領域][核種] で
    /// あること、輸送できる核種がチェーンの核種であること、少数群定数が状態 1 つのテーブルであることを確かめる
    pub fn from_msgpack(bytes: &[u8]) -> anyhow::Result<Self> {
        let file: HistoryFileIn = msgpack::from_value(msgpack::read_value(bytes)?)?;
        if file.format != FORMAT {
            bail!("expected format {FORMAT}, found {}", file.format);
        }
        let (nuclides, regions, transportable) = (file.nuclides, file.regions, file.transportable);
        // 燃焼の輸送は燃焼の領域のない格子で止まる（`run`）ので、領域のない履歴は書かれない
        if regions.is_empty() {
            bail!("the history has no depletion regions");
        }
        // 数密度の列は核種の名前で引く
        for (k, name) in nuclides.iter().enumerate() {
            if nuclides[..k].contains(name) {
                bail!("the history lists {name} twice");
            }
        }
        if let Some(name) = transportable.iter().find(|n| !nuclides.contains(n)) {
            bail!("the transportable {name} is not a nuclide of the history");
        }
        let burnups: Vec<f64> = file.points.iter().map(|p| p.burnup).collect();
        if burnups.first() != Some(&0.0) || !burnups.windows(2).all(|w| w[0] < w[1]) {
            bail!("the history burnups must increase from 0, got {burnups:?}");
        }
        let points = file
            .points
            .into_iter()
            .map(|p| {
                if p.densities.len() != regions.len()
                    || p.densities.iter().any(|n| n.len() != nuclides.len())
                {
                    bail!(
                        "point at {} MWd/kg: the densities need {} regions of {} nuclides",
                        p.burnup,
                        regions.len(),
                        nuclides.len()
                    );
                }
                let few_group = p
                    .few_group
                    .map(|value| -> anyhow::Result<FewGroup> {
                        let mut bytes = Vec::new();
                        rmpv::encode::write_value(&mut bytes, &value)?;
                        let table = Table::from_msgpack(&bytes)?;
                        match (table.axes.is_empty(), table.states.as_slice()) {
                            (true, [one]) => Ok(one.clone()),
                            _ => bail!("the few-group constants of a point must be one state"),
                        }
                    })
                    .transpose()
                    .with_context(|| format!("point at {} MWd/kg", p.burnup))?;
                Ok(HistoryPoint {
                    time: p.time,
                    burnup: p.burnup,
                    k: p.k_inf,
                    densities: p.densities,
                    pin_powers: p.pin_powers,
                    few_group,
                    transported: transportable.len(),
                })
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        let s = &file.seconds;
        Ok(Self {
            nuclides,
            regions,
            transportable,
            points,
            transport: TransportTimes {
                shielding: seconds("shielding", s.shielding)?,
                moc: seconds("moc", s.moc)?,
                rates: seconds("rates", s.rates)?,
            },
            depletion: DepletionTimes {
                matrix: seconds("matrix", s.matrix)?,
                cram: seconds("cram", s.cram)?,
            },
            solves: file.solves,
            shieldings: file.shieldings,
        })
    }
    /// フィールド名付きの map 形式の MessagePack: `{format: "xlc-depletion-history", nuclides, regions: [{material, cell,
    /// regions, volume}], transportable, points: [{time, burnup, k_inf, densities, pin_powers, few_group}], seconds: {shielding, moc, rates,
    /// matrix, cram}（時間の内訳 (s)）, solves, shieldings}`
    pub fn to_msgpack(&self) -> anyhow::Result<Vec<u8>> {
        let points = self
            .points
            .iter()
            .map(|p| {
                let few_group = p
                    .few_group
                    .as_ref()
                    .map(|few| -> anyhow::Result<rmpv::Value> {
                        let bytes = Table::single(few.clone()).to_msgpack()?;
                        Ok(msgpack::read_value(&bytes)?)
                    })
                    .transpose()?;
                Ok(PointFile {
                    time: p.time,
                    burnup: p.burnup,
                    k_inf: p.k,
                    densities: &p.densities,
                    pin_powers: &p.pin_powers,
                    few_group,
                })
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        Ok(msgpack::to_named(&HistoryFile {
            format: FORMAT,
            nuclides: &self.nuclides,
            regions: &self.regions,
            transportable: &self.transportable,
            points,
            seconds: SecondsFile {
                shielding: self.transport.shielding.as_secs_f64(),
                moc: self.transport.moc.as_secs_f64(),
                rates: self.transport.rates.as_secs_f64(),
                matrix: self.depletion.matrix.as_secs_f64(),
                cram: self.depletion.cram.as_secs_f64(),
            },
            solves: self.solves,
            shieldings: self.shieldings,
        })?)
    }
}

/// 1 回の輸送で点に残す値
struct State {
    k: f64,
    /// 点に残す輸送だけ（中点の輸送は `None`）
    pin_powers: Option<Vec<Vec<f64>>>,
    few_group: Option<FewGroup>,
    transported: usize,
}

/// 自己遮蔽の結果 `start` の微視的断面積と SPH 因子のまま、燃焼の領域 `regions` の物質を組成 `materials`（[燃焼の領域]。
/// 核種の並びは `start` を解いたときと同じ）で作り直す（ADR-0197）。物質は、自己遮蔽と同じく、領域の核種の実効断面積
/// （`region_cross_sections`）から作った群定数を数密度で混ぜ、SPH 因子を掛ける。ほかの領域の物質はそのまま
pub(crate) fn remix_depleted(
    library: &Library,
    start: &FuelShielding,
    lattice: &Lattice,
    regions: &[DepletionRegion],
    materials: &[MaterialInput],
) -> anyhow::Result<FuelShielding> {
    if materials.len() != regions.len() {
        bail!(
            "{} materials for {} depletion regions",
            materials.len(),
            regions.len()
        );
    }
    // 領域ごとの SPH 因子 μ（[領域][群]）。物理的な中性子束は MOC の中性子束に μ を掛けたものなので、1 の中性子束から引く
    let ones = vec![vec![1.0; library.groups.len()]; lattice.regions()];
    let sph = start.physical_flux(lattice, &ones)?;
    // 物質の番号ごとの燃焼の領域。1 つの物質を、別の燃焼の領域か燃焼させない領域と共有していれば、作り直せない
    let mut depleted: Vec<Option<usize>> = vec![None; lattice.regions()];
    for (k, d) in regions.iter().enumerate() {
        for &r in &d.regions {
            depleted[r] = Some(k);
        }
    }
    let mut owner: Vec<Option<Option<usize>>> = vec![None; start.materials.len()];
    let mut first: Vec<Option<usize>> = vec![None; start.materials.len()];
    for (r, &k) in depleted.iter().enumerate() {
        let m = start.region_material[r];
        match owner[m] {
            None => {
                owner[m] = Some(k);
                first[m] = Some(r);
            }
            Some(o) if o == k => {
                let f = first[m].ok_or_else(|| anyhow!("material {m} has no first region"))?;
                if k.is_some() && sph[r] != sph[f] {
                    bail!("regions {f} and {r} share material {m} with different SPH factors");
                }
            }
            Some(_) => bail!(
                "material {m} of region {r} is shared across depletion regions or with a region not depleted"
            ),
        }
    }
    let rebuilt = (0..start.materials.len())
        .into_par_iter()
        .filter_map(|m| match (owner[m], first[m]) {
            (Some(Some(k)), Some(r)) => Some((m, k, r)),
            _ => None,
        })
        .map(|(m, k, r)| -> anyhow::Result<(usize, Material)> {
            let material = &materials[k];
            let cross_sections = start
                .region_cross_sections(r)
                .ok_or_else(|| anyhow!("region {r} has no shielded cross sections"))?;
            if cross_sections.len() != material.composition.len() {
                bail!(
                    "region {r}: {} shielded nuclides for the {} nuclides of `{}`",
                    cross_sections.len(),
                    material.composition.len(),
                    material.name
                );
            }
            let nuclides = material
                .composition
                .iter()
                .zip(cross_sections)
                .map(|((name, density), s)| {
                    Ok((
                        shielded_nuclide(library, name, material.temperature, s)?,
                        *density,
                    ))
                })
                .collect::<anyhow::Result<Vec<_>>>()?;
            let mixed =
                Material::from_nuclides(&nuclides, library.groups.len())?.scaled(&sph[r])?;
            Ok((m, mixed))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let mut remixed = start.clone();
    for (m, material) in rebuilt {
        remixed.materials[m] = material;
    }
    Ok(remixed)
}

/// 格子を `[depletion]` の条件で燃焼させる（ADR-0173）。`region_materials` は入力の領域ごとの物質の名前、`chosen` はその
/// 名前で選んだ対称操作。反応ごとの断面積は `reactions` から引く
#[allow(clippy::too_many_arguments)]
pub fn run(
    input: &LatticeInput,
    depletion: &DepletionInput,
    library: &Library,
    chain: &Chain,
    tracking: &Tracking,
    region_materials: &[String],
    chosen: &Chosen,
    reactions: &dyn ReactionSource,
) -> anyhow::Result<History> {
    run_with_moc(
        input,
        depletion,
        library,
        chain,
        tracking,
        region_materials,
        chosen,
        reactions,
        &Options::default(),
    )
}

/// `run` と同じく燃焼させ、MOC を `moc` の打ち切り条件と加速で解く（散乱の次数は `run` と同じく入力か物質から）
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_with_moc(
    input: &LatticeInput,
    depletion: &DepletionInput,
    library: &Library,
    chain: &Chain,
    tracking: &Tracking,
    region_materials: &[String],
    chosen: &Chosen,
    reactions: &dyn ReactionSource,
    moc: &Options,
) -> anyhow::Result<History> {
    chain.check()?;
    let lattice = &tracking.lattice;
    let volumes = lattice.volumes();
    let index: HashMap<&str, usize> = chain
        .nuclides
        .iter()
        .enumerate()
        .map(|(i, n)| (n.name.as_str(), i))
        .collect();
    let representatives = chosen.shielding().cell_representatives();
    let regions = depletion_regions(
        lattice,
        region_materials,
        &representatives,
        &depletion.materials,
    );
    for name in &depletion.materials {
        if !regions.iter().any(|d| &d.material == name) {
            bail!("depleted material `{name}` is not in the lattice");
        }
    }
    let base = |name: &str| -> anyhow::Result<&MaterialInput> {
        input
            .materials
            .iter()
            .find(|m| m.name == name)
            .ok_or_else(|| anyhow!("depleted material `{name}` is not written with nuclides"))
    };
    check_region_names(&input.materials, &regions)?;
    // 初めの数密度（チェーンの順）
    let initial = regions
        .iter()
        .map(|d| {
            let mut n = vec![0.0; chain.nuclides.len()];
            for (name, density) in &base(&d.material)?.composition {
                let i = index.get(name.as_str()).ok_or_else(|| {
                    anyhow!("material `{}`: {name} is not in the chain", d.material)
                })?;
                n[*i] = *density;
            }
            Ok(n)
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    // 輸送する核種の候補: チェーンの核種のうち、ライブラリがすべての燃焼の領域の温度で持つもの。ライブラリにあっても領域の温度の
    // データがない核種は、崩壊と生成だけを追い、輸送のたびに warn で示す（ADR-0173）
    let mut temperatures: Vec<f64> = Vec::new();
    for d in &regions {
        let t = base(&d.material)?.temperature;
        if !temperatures.contains(&t) {
            temperatures.push(t);
        }
    }
    // 輸送できる核種は領域の温度で決まる。燃焼させる物質の温度が違うときの集合（燃料と吸収体で別の集合にするか）は決めていない
    if temperatures.len() > 1 {
        bail!(
            "the depleted materials have several temperatures {temperatures:?}; the transported nuclides of several \
             depletion temperatures are not decided (#0139)"
        );
    }
    let in_library: Vec<bool> = chain
        .nuclides
        .iter()
        .map(|n| library.nuclides.iter().any(|d| d.name == n.name))
        .collect();
    let transportable: Vec<bool> = chain
        .nuclides
        .iter()
        .zip(&in_library)
        .map(|(n, &known)| {
            known
                && temperatures
                    .iter()
                    .all(|&t| has_data_at(library, &n.name, t))
        })
        .collect();
    // 初めの重金属の質量 (g、高さ 1 cm あたり)
    let mass = heavy_metal_mass(&regions, &initial, |i| {
        let nuclide = &chain.nuclides[i];
        if nuclide.z < HEAVY_METAL_Z {
            return Ok(None);
        }
        library
            .nuclides
            .iter()
            .find(|x| x.name == nuclide.name)
            .map(|data| Some(data.awr))
            .ok_or_else(|| anyhow!("{}: no mass in the library", nuclide.name))
    })?;
    let power = depletion.power * mass;
    let renamed = solved_names(
        lattice,
        region_materials,
        &representatives,
        &depletion.materials,
    );
    let groups = &library.groups;
    let mut transport = TransportTimes::default();
    let mut solves = 0;
    // 燃料のピンの Dancoff 係数は、燃料でない物質の全断面積が直前の輸送と同じ間だけ持ち越す。燃焼させる物質が燃料だけなら
    // 初めの輸送で 1 度だけ解き、燃料でない物質（吸収体）を燃焼させれば輸送ごとに解く（ADR-0194）
    let mut reuse = ShieldingReuse::new(tracking);
    // 直前の輸送の収束した MOC の状態。次の輸送の反復をここから始める（ADR-0195）
    let mut previous: Option<Start> = None;
    let mut shieldings = 0;
    // `StepStart` の、ステップの始まりの輸送の自己遮蔽。次の中点の輸送で使う（ADR-0197）
    let mut held: Option<FuelShielding> = None;
    let solve = |n: &[Vec<f64>], kept: bool| -> anyhow::Result<Evaluation<State>> {
        solves += 1;
        let start = Instant::now();
        // 負の数密度は近似の誤差の水準だけを許し、輸送には入れない（ADR-0173）
        check_negative(&regions, n, |i| chain.nuclides[i].name.clone())?;
        // 輸送できる核種はすべて、数密度が 0 でも輸送して反応率を求める（OpenMC も数密度を掛けない tally で、0 の核種に 1 原子
        // あたりの反応率を与える。始まりに 0 の Xe-135 も予測子で燃え尽きる）。どの領域でも同じ核種の並び（自己遮蔽は 1 つの
        // セルの燃料のリングに同じ並びを求める）
        let present = |i: usize| n.iter().any(|n| n[i] > 0.0);
        let composition: Vec<usize> = (0..chain.nuclides.len())
            .filter(|&i| transportable[i])
            .collect();
        // 崩壊だけを追う核種の数密度の和（領域の体積の重み）
        let total = |i: usize| -> f64 {
            regions
                .iter()
                .zip(n)
                .map(|(d, n)| d.volume * n[i].max(0.0))
                .sum::<f64>()
                / regions.iter().map(|d| d.volume).sum::<f64>()
        };
        let lacking: Vec<usize> = (0..chain.nuclides.len())
            .filter(|&i| in_library[i] && !transportable[i] && present(i))
            .collect();
        if !lacking.is_empty() {
            tracing::warn!(
                nuclides = lacking
                    .iter()
                    .map(|&i| chain.nuclides[i].name.as_str())
                    .collect::<Vec<_>>()
                    .join(" "),
                density = lacking.iter().map(|&i| total(i)).sum::<f64>(),
                "decay only: no library data at the temperature of the depletion regions"
            );
        }
        let absent: Vec<usize> = (0..chain.nuclides.len())
            .filter(|&i| !in_library[i] && present(i))
            .collect();
        tracing::info!(
            nuclides = absent.len(),
            density = absent.iter().map(|&i| total(i)).sum::<f64>(),
            "decay only: not in the library"
        );
        let mut materials: Vec<MaterialInput> = input
            .materials
            .iter()
            .filter(|m| !depletion.materials.contains(&m.name))
            .cloned()
            .collect();
        for (d, n) in regions.iter().zip(n) {
            materials.push(region_material(
                d,
                base(&d.material)?.temperature,
                n,
                &composition,
                |i| chain.nuclides[i].name.clone(),
            ));
        }
        // ce_cm が点に残さない輸送（`kept` でないもの）は、ステップの中点の輸送
        let shielded = match (depletion.shielding, kept) {
            (ShieldingFrequency::StepStart, false) => {
                let at_start = held.take().ok_or_else(|| {
                    anyhow!("the midpoint transport has no shielding from the step start")
                })?;
                let depleted = &materials[materials.len() - regions.len()..];
                remix_depleted(library, &at_start, lattice, &regions, depleted)
                    .context("remixing the depleted materials at the midpoint")?
            }
            _ => {
                shieldings += 1;
                let shielded = shield_fuel_reusing(
                    library,
                    &mut reuse,
                    &materials,
                    &renamed,
                    chosen.shielding(),
                    &input.shielding_options(),
                )
                .context("self-shielding")?;
                // 作り直しの前提（自己遮蔽の物質が、領域の実効断面積と数密度から作った群定数に SPH 因子を掛けたもの）を、
                // 最初の輸送で自己遮蔽の結果そのものと比べて確かめる。前提が崩れれば、中点の物質が始まりと食い違う
                if depletion.shielding == ShieldingFrequency::StepStart && solves == 1 {
                    let depleted = &materials[materials.len() - regions.len()..];
                    let remixed =
                        remix_depleted(library, &shielded, lattice, &regions, depleted)
                            .context("remixing the depleted materials at the first transport")?;
                    if remixed.materials != shielded.materials {
                        bail!(
                            "the depleted materials remixed from their shielded cross sections differ from the \
                             self-shielding"
                        );
                    }
                }
                shielded
            }
        };
        transport.shielding += start.elapsed();
        let start = Instant::now();
        let region_materials: Vec<&Material> = shielded
            .region_material
            .iter()
            .map(|&i| &shielded.materials[i])
            .collect();
        let options = Options {
            scattering_order: input
                .scattering_order
                .unwrap_or_else(|| available_order(&region_materials)),
            ..*moc
        };
        let (solution, converged) = solve_symmetric_from(
            tracking,
            &region_materials,
            &options,
            chosen.moc(),
            previous.as_ref(),
        )?;
        previous = Some(converged);
        transport.moc += start.elapsed();
        let start = Instant::now();
        let flux = shielded.physical_flux(lattice, &solution.flux)?;
        let mut rates = Vec::with_capacity(regions.len());
        let mut fission_power = 0.0;
        for (d, n) in regions.iter().zip(n) {
            let temperature = base(&d.material)?.temperature;
            let mut members = Vec::with_capacity(d.regions.len());
            for &r in &d.regions {
                let cross_sections = shielded.region_cross_sections(r).ok_or_else(|| {
                    anyhow!(
                        "region {r} of material `{}` has no cross sections",
                        d.material
                    )
                })?;
                if cross_sections.len() != composition.len() {
                    bail!(
                        "region {r} of material `{}`: {} shielded nuclides for {} transported",
                        d.material,
                        cross_sections.len(),
                        composition.len()
                    );
                }
                members.push(Member {
                    volume: volumes[r],
                    flux: &flux[r],
                    cross_sections,
                });
            }
            let region = region_rates(
                chain,
                &composition,
                temperature,
                &members,
                groups.bounds(),
                reactions,
            )?;
            for (i, r) in region.iter().enumerate() {
                if let Some(f) = &chain.nuclides[i].fission {
                    fission_power +=
                        d.volume * n[i].max(0.0) / BARN * r.fission * f.q * ELECTRON_VOLT;
                }
            }
            rates.push(region);
        }
        if !(fission_power > 0.0 && fission_power.is_finite()) {
            bail!("the depleted materials have no fission power ({fission_power})");
        }
        let scale = power / fission_power;
        for region in &mut rates {
            for r in region.iter_mut() {
                r.fission *= scale;
                for channels in &mut r.reactions {
                    for c in channels {
                        *c *= scale;
                    }
                }
            }
        }
        // 中点の輸送のピン出力と少数群定数は点に残らないので作らない
        let pin_powers = if kept {
            Some(pin_rows(
                &cell_fission_rates(&solution, lattice, &region_materials)?,
                lattice,
            ))
        } else {
            None
        };
        let few_group = match &depletion.few_group {
            Some(DepletionFewGroup { bounds, pins }) if kept => {
                let boundaries = tracking.boundaries();
                let mut few =
                    homogenize(&solution, &region_materials, groups, bounds, &boundaries)?;
                if *pins {
                    few.pins = Some(pin_shapes(
                        &solution,
                        lattice,
                        &region_materials,
                        groups,
                        bounds,
                        &boundaries,
                    )?);
                }
                Some(few)
            }
            _ => None,
        };
        transport.rates += start.elapsed();
        if depletion.shielding == ShieldingFrequency::StepStart && kept {
            held = Some(shielded);
        }
        tracing::info!(
            solve = solves,
            k = solution.k,
            iterations = solution.iterations,
            nuclides = composition.len(),
            "depletion transport"
        );
        Ok(Evaluation {
            rates,
            extra: State {
                k: solution.k,
                pin_powers,
                few_group,
                transported: composition.len(),
            },
        })
    };
    let steps = step_times(depletion.power, &depletion.burnup);
    let (points, depletion_times) = ce_cm(chain, initial, &steps, solve)?;
    let transportable = chain
        .nuclides
        .iter()
        .zip(&transportable)
        .filter(|(_, t)| **t)
        .map(|(n, _)| n.name.clone())
        .collect();
    let burnup = std::iter::once(0.0).chain(depletion.burnup.iter().copied());
    let points = points
        .into_iter()
        .zip(burnup)
        .map(|(p, burnup)| {
            Ok(HistoryPoint {
                time: p.time,
                burnup,
                k: p.extra.k,
                densities: p.densities,
                pin_powers: p.extra.pin_powers.ok_or_else(|| {
                    anyhow!("the transport of the point at {burnup} MWd/kg has no pin powers")
                })?,
                few_group: p.extra.few_group,
                transported: p.extra.transported,
            })
        })
        .collect::<anyhow::Result<_>>()?;
    Ok(History {
        nuclides: chain.nuclides.iter().map(|n| n.name.clone()).collect(),
        regions,
        transportable,
        points,
        transport,
        depletion: depletion_times,
        solves,
        shieldings,
    })
}

#[cfg(test)]
pub(crate) mod tests;
