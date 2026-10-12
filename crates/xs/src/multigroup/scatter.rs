//! 散乱行列（P0）。
//!
//! 弾性散乱と離散準位の非弾性散乱は、静止した標的との2体反応として扱う。重心系の角度分布は MF4 で、
//! なければ重心系で等方とする。
//! 熱領域の弾性散乱は、指定があれば自由ガス模型で扱う（標的の熱運動による上方散乱を含む）。
//! 連続的な中性子放出反応（MT91, (n,2n) など）は MF6 のエネルギー分布と多重度を使う。
//! MF6 がなく MF4/MF5 で与えられた反応は、MF5 のエネルギー分布と、反応で決まる多重度を使い、
//! MF4 の実験室系の角度分布をエネルギーと独立に掛ける。

use std::collections::BTreeMap;
use xlc_numeric::quadrature::{GAUSS_LEGENDRE_8, legendre};

use super::free_gas::free_gas_group_moments;
use super::outgoing::uniform_fractions;
use super::{
    FreeGas, GroupStructure, GroupedLaw, GroupedOutgoing, ThermalSab, Weight, group_integrals,
    group_integrals_multi,
};

use crate::endf::{Angular, Distribution, Frame, Mf4, Mf5, Mf6, Tab1};

/// 重心系等方の2体反応（Q 値 q、標的の質量比 awr）で、入射エネルギー e の中性子の
/// 出射エネルギーの範囲。出射エネルギーはこの範囲で一様に分布する。しきい値未満なら `None`。
pub fn two_body_range(e: f64, awr: f64, q: f64) -> Option<(f64, f64)> {
    // E' = E/(A+1)² (1 + A²β² + 2Aβμ)、β² = 1 + (A+1)Q/(A E)、μ は重心系の散乱角の余弦
    let beta2 = 1.0 + (awr + 1.0) * q / (awr * e);
    // しきい値 −Q (A + 1)/A 以上では β² ≥ 0 だが、ちょうどしきい値では丸めでわずかに負になることがある。
    // しきい値の断面積（xs はしきい値から表を始める。ADR-0100）の中性子を失わないよう、0 とする
    let beta2 = if beta2 < 0.0 && e >= -q * (awr + 1.0) / awr {
        0.0
    } else {
        beta2
    };
    if beta2 < 0.0 {
        return None;
    }
    let beta = beta2.sqrt();
    let c = e / (awr + 1.0).powi(2);
    Some((
        c * (1.0 - awr * beta).powi(2),
        c * (1.0 + awr * beta).powi(2),
    ))
}

/// 連続的な中性子放出反応の出射分布。
#[derive(Debug, Clone)]
pub enum Emission {
    /// MF6（多重度と出射分布）
    Mf6(Mf6),
    /// MF5（出射エネルギー分布。多重度は反応で決まる）
    Mf5(Mf5),
}

/// MF5 で与えられた反応の中性子の多重度（ENDF-6 の反応の定義による）。
/// 定義から決まらない反応（MT5 など）は `None`。
pub fn multiplicity(mt: u32) -> Option<f64> {
    match mt {
        4 | 22 | 23 | 28 | 29 | 32..=36 | 44 | 45 | 51..=91 => Some(1.0),
        11 | 16 | 24 | 30 | 41 => Some(2.0),
        17 | 25 | 42 => Some(3.0),
        37 => Some(4.0),
        _ => None,
    }
}

/// 熱領域の散乱の扱い。どちらも cutoff 未満の入射エネルギーで、弾性散乱の静止標的の2体反応を置き換える。
#[derive(Debug, Clone)]
pub enum Thermal<'a> {
    /// 自由ガス模型（温度に拡げた MT2 の断面積に、標的の熱運動による出射分布を掛ける）
    FreeGas(FreeGas),
    /// `sab.cutoff()` 未満は S(α, β)（MF3 の弾性散乱断面積は使わず、核そのものが断面積を与える）、その上
    /// `free_gas.cutoff` 未満は自由ガス模型（ADR-0133）
    Sab {
        sab: ThermalSab<'a>,
        free_gas: FreeGas,
    },
}

/// 散乱行列を作るときのエラー。
#[derive(Debug, thiserror::Error)]
pub enum ScatterError {
    #[error("free gas needs positive kT and a finite cutoff (kT {kt}, cutoff {cutoff})")]
    FreeGas { kt: f64, cutoff: f64 },
    #[error("free gas above S(α, β) must reach at least its limit {sab} eV, got {free_gas} eV")]
    FreeGasBelowSab { sab: f64, free_gas: f64 },
    #[error("MT{0}: unknown neutron multiplicity")]
    Multiplicity(u32),
    #[error("MT{0}: missing Q value")]
    QValue(u32),
    #[error("MT{0}: no outgoing neutron distribution")]
    NoDistribution(u32),
    #[error("MT{0}: a two-body (LAW=2) neutron distribution outside the discrete levels MT51-90")]
    TwoBody(u32),
    #[error("MT{0}: two-body angular distributions must be in the centre-of-mass frame")]
    LabTwoBody(u32),
    #[error("MT{0}: angular distributions of MF5 emissions must be in the laboratory frame")]
    CenterOfMassContinuum(u32),
}

/// 散乱行列を作る反応のデータ（どれもキーは MT）。
#[derive(Debug, Clone, Copy)]
pub struct Reactions<'a> {
    /// ポイントワイズの断面積
    pub tables: &'a BTreeMap<u32, Tab1>,
    /// 離散準位の非弾性散乱（MT51〜90）の Q 値 (eV)
    pub q_values: &'a BTreeMap<u32, f64>,
    /// 連続的な中性子放出反応の出射分布
    pub emissions: &'a BTreeMap<u32, Emission>,
    /// 角度分布（MF4）
    pub angular: &'a BTreeMap<u32, Mf4>,
}

/// 散乱行列 σ(g→g') (barn)。[入射の群][出射の群]。`thermal` は熱領域の弾性散乱の扱い
/// （`None` なら全域で静止標的）。
///
/// 中性子を放出するのに分布のない反応があればエラーにする（中性子を黙って失わないため）。
/// P0 の散乱行列 σ(g→g')（[`scattering_moments`] の 0 次）。
pub fn scattering_matrix(
    groups: &GroupStructure,
    weight: &Weight,
    awr: f64,
    reactions: &Reactions,
    thermal: Option<&Thermal>,
) -> Result<Vec<Vec<f64>>, ScatterError> {
    let moments = scattering_moments(groups, weight, awr, reactions, thermal, 0)?;
    Ok(moments.into_iter().next().unwrap_or_default())
}

/// 散乱行列のルジャンドルのモーメント σ_l(g→g') = ∫_g ∫_g' σ(E→E', μ) P_l(μ) dE' w(E) dE / ∫_g w dE
/// （[l][入射の群][出射の群]、l = 0..=order、μ は実験室系の散乱角の余弦）。
///
/// `angular` は MT ごとの角度分布（MF4）。2体反応（弾性散乱と離散準位の非弾性散乱）は重心系の分布と運動学から、
/// MF5 の連続的な放出は実験室系の分布から、熱領域（S(α, β) と自由ガス）はその角度分布から P1 以上を作る。
/// MF6 の連続的な放出の P1 以上はまだ角度分布を入れておらず 0 とする（ADR-0082）。
#[tracing::instrument(level = "debug", skip_all)]
pub fn scattering_moments(
    groups: &GroupStructure,
    weight: &Weight,
    awr: f64,
    reactions: &Reactions,
    thermal: Option<&Thermal>,
    order: usize,
) -> Result<Vec<Vec<Vec<f64>>>, ScatterError> {
    let Reactions {
        tables,
        q_values,
        emissions,
        angular,
    } = *reactions;
    // 2体反応の角度分布は重心系で使う（対象の評価済みデータはすべて重心系で与えている）
    for (&mt, mf4) in angular {
        if (mt == 2 || (51..=90).contains(&mt)) && mf4.frame != Frame::CenterOfMass {
            return Err(ScatterError::LabTwoBody(mt));
        }
        // MF5 の放出は出射エネルギーと角度を独立に掛けるので、角度分布は実験室系で使う。
        // 重心系で等方（O-18 など）は、これまでどおり P1 以上を 0 とする（重心系の運動の分を入れない近似）
        if matches!(emissions.get(&mt), Some(Emission::Mf5(_)))
            && mf4.frame != Frame::Laboratory
            && mf4.angular != Angular::Isotropic
        {
            return Err(ScatterError::CenterOfMassContinuum(mt));
        }
    }
    if let Some(Thermal::FreeGas(free_gas) | Thermal::Sab { free_gas, .. }) = thermal
        && !(free_gas.kt > 0.0 && free_gas.kt.is_finite() && free_gas.cutoff.is_finite())
    {
        return Err(ScatterError::FreeGas {
            kt: free_gas.kt,
            cutoff: free_gas.cutoff,
        });
    }
    if let Some(Thermal::Sab { sab, free_gas }) = thermal
        && free_gas.cutoff < sab.cutoff()
    {
        return Err(ScatterError::FreeGasBelowSab {
            sab: sab.cutoff(),
            free_gas: free_gas.cutoff,
        });
    }
    let ng = groups.len();
    let flux = group_integrals(groups, weight, &[], |_| 1.0);
    let has_n2n_partials = tables.keys().any(|mt| (875..=891).contains(mt));
    let has_inelastic_partials = tables.keys().any(|mt| (51..=91).contains(mt));
    let mut moments = vec![vec![vec![0.0; ng]; ng]; order + 1];
    for (&mt, table) in tables {
        // 部分反応のある MT4 と MT16 は和なので、部分反応の側で数える
        if !emits_neutrons(mt)
            || (mt == 4 && has_inelastic_partials)
            || (mt == 16 && has_n2n_partials)
        {
            continue;
        }
        // 弾性散乱は MF6 の有無にかかわらず2体反応（熱領域では自由ガス）で扱う
        let rows = if mt == 2 {
            elastic_rows(groups, weight, awr, table, thermal, angular.get(&2), order)
        } else if let Some(Emission::Mf6(data)) = emissions.get(&mt)
            && !data.two_body
        {
            isotropic(emission_rows(groups, weight, table, data), order)
        } else if let Some(Emission::Mf6(_)) = emissions.get(&mt)
            && !(51..=90).contains(&mt)
        {
            return Err(ScatterError::TwoBody(mt));
        } else if let Some(Emission::Mf5(data)) = emissions.get(&mt) {
            let count = multiplicity(mt).ok_or(ScatterError::Multiplicity(mt))?;
            spectrum_moments(groups, weight, table, data, count, angular.get(&mt), order)
        } else if (51..=90).contains(&mt) {
            let q = *q_values.get(&mt).ok_or(ScatterError::QValue(mt))?;
            two_body_rows(groups, weight, awr, q, table, angular.get(&mt), order)
        } else {
            return Err(ScatterError::NoDistribution(mt));
        };
        for (sum, part) in moments.iter_mut().zip(rows) {
            for (row, values) in sum.iter_mut().zip(part) {
                for (s, v) in row.iter_mut().zip(values) {
                    *s += v;
                }
            }
        }
    }
    for matrix in &mut moments {
        for (row, phi) in matrix.iter_mut().zip(&flux) {
            for v in row.iter_mut() {
                *v /= phi;
            }
        }
    }
    Ok(moments)
}

/// P0 の行だけを持ち、P1 以上を 0 とした [l][g][g']（角度分布をまだ入れていない反応）
fn isotropic(rows: Vec<Vec<f64>>, order: usize) -> Vec<Vec<Vec<f64>>> {
    let ng = rows.len();
    std::iter::once(rows)
        .chain((0..order).map(|_| vec![vec![0.0; ng]; ng]))
        .collect()
}

/// [g][l × ng + g'] の積分を [l][g][g'] に並べ替える
fn by_order(rows: Vec<Vec<f64>>, ng: usize, order: usize) -> Vec<Vec<Vec<f64>>> {
    (0..=order)
        .map(|l| {
            rows.iter()
                .map(|row| row[l * ng..(l + 1) * ng].to_vec())
                .collect()
        })
        .collect()
}

/// 角度分布（MF4）を散乱行列に使う反応: 弾性散乱、離散準位（2体反応）と、MF5 で与えられた連続的な放出。
/// MF6 の反応は MF6 の中の角度分布を使い、核分裂の中性子は等方とする
pub fn uses_angular(mt: u32, emission: Option<&Emission>) -> bool {
    mt == 2 || (51..=90).contains(&mt) || matches!(emission, Some(Emission::Mf5(_)))
}

/// 中性子を放出する反応（核分裂を除く）
pub fn emits_neutrons(mt: u32) -> bool {
    matches!(
        mt,
        2 | 4
            | 5
            | 11
            | 16
            | 17
            | 22..=25
            | 28..=30
            | 32..=37
            | 41
            | 42
            | 44
            | 45
            | 51..=91
            // ENDF-6 の MT152〜200 のうち、中性子を出さない 155, 182, 191〜193, 197 を除く
            | 152..=154
            | 156..=181
            | 183..=190
            | 194..=196
            | 198..=200
            | 875..=891
    )
}

/// 弾性散乱の ∫_g σ(E) w(E) P_l(E→g') dE（[l][g][g']）。`thermal` の S(α, β) と自由ガスの cutoff 未満では熱領域の
/// 扱い、それ以上では静止標的との2体反応とする。
fn elastic_rows(
    groups: &GroupStructure,
    weight: &Weight,
    awr: f64,
    table: &Tab1,
    thermal: Option<&Thermal>,
    angular: Option<&Mf4>,
    order: usize,
) -> Vec<Vec<Vec<f64>>> {
    let Some(thermal) = thermal else {
        return two_body_rows(groups, weight, awr, 0.0, table, angular, order);
    };
    let ng = groups.len();
    let (sab, free_gas) = match thermal {
        Thermal::FreeGas(free_gas) => (None, free_gas),
        Thermal::Sab { sab, free_gas } => (Some(sab), free_gas),
    };
    let sab_cutoff = sab.map_or(0.0, ThermalSab::cutoff);
    let mut breakpoints = table.x.clone();
    breakpoints.extend(two_body_kinks(groups, awr, 0.0));
    breakpoints.extend(angular.map(Mf4::energies).unwrap_or_default());
    breakpoints.extend([sab_cutoff, free_gas.cutoff]);
    let add = |out: &mut [f64], moments: Vec<Vec<f64>>, scale: f64| {
        for (l, row) in moments.into_iter().enumerate() {
            for (o, v) in out[l * ng..(l + 1) * ng].iter_mut().zip(row) {
                *o += scale * v;
            }
        }
    };
    let rows = group_integrals_multi(groups, weight, &breakpoints, ng * (order + 1), |e, out| {
        let sigma = table.eval(e).unwrap_or(0.0);
        match sab {
            // S(α, β) の核は断面積そのものなので、MF3 の弾性散乱断面積は使わない
            Some(sab) if e < sab_cutoff => add(out, sab.group_moments(groups, e, order), 1.0),
            // 温度に拡げた断面積はすでに熱運動による実効断面積の増加 N(E) を含むので、
            // 出射分布は和 N(E) で規格化する（NJOY THERMR と同じ）
            _ if e < free_gas.cutoff => {
                let moments = free_gas_group_moments(groups, e, awr, free_gas.kt, order);
                let total: f64 = moments[0].iter().sum();
                add(out, moments, sigma / total);
            }
            _ => add(
                out,
                two_body_moments(groups, e, awr, 0.0, sigma, order, angular),
                1.0,
            ),
        }
    });
    by_order(rows, ng, order)
}

/// 2体反応の ∫_g σ(E) w(E) P_l(E→g') dE（[l][g][g']）。
fn two_body_rows(
    groups: &GroupStructure,
    weight: &Weight,
    awr: f64,
    q: f64,
    table: &Tab1,
    angular: Option<&Mf4>,
    order: usize,
) -> Vec<Vec<Vec<f64>>> {
    let ng = groups.len();
    let mut breakpoints = table.x.clone();
    breakpoints.extend(two_body_kinks(groups, awr, q));
    breakpoints.extend(angular.map(Mf4::energies).unwrap_or_default());
    let rows = group_integrals_multi(groups, weight, &breakpoints, ng * (order + 1), |e, out| {
        let sigma = table.eval(e).unwrap_or(0.0);
        if sigma == 0.0 {
            return;
        }
        for (l, row) in two_body_moments(groups, e, awr, q, sigma, order, angular)
            .into_iter()
            .enumerate()
        {
            for (o, v) in out[l * ng..(l + 1) * ng].iter_mut().zip(row) {
                *o += v;
            }
        }
    });
    by_order(rows, ng, order)
}

/// 2体反応で、入射エネルギー `e` の中性子が群 g' に出る確率に、実験室系の散乱角の余弦 μ のルジャンドル多項式
/// P_l(μ) を掛けたもの [l][g']（`scale` 倍）。`angular` は重心系の角度分布（MF4。`None` か等方なら重心系で等方で、
/// l = 0 は一様な振り分け）。
///
/// 出射の速さの比 s = (A+1)√(E'/E) で積分する。重心系の余弦は μ_c = (s² − 1 − A²β²)/(2Aβ) で、
/// dμ_c = s ds/(Aβ)、実験室系の余弦は μ = (s² + 1 − A²β²)/(2s) と s の滑らかな式になる（A = 1 でも特異点がない）。
/// 角度分布があるときは f(μ_c) を掛けて P0 も同じく積分し、P0 の和が `scale` になるよう規格化する
/// （分布の規格化と数値積分の誤差で粒子数を変えないため）。
fn two_body_moments(
    groups: &GroupStructure,
    e: f64,
    awr: f64,
    q: f64,
    scale: f64,
    order: usize,
    angular: Option<&Mf4>,
) -> Vec<Vec<f64>> {
    // 入射エネルギーごとに一度だけ分布を求める
    let angular = angular
        .filter(|a| a.angular != Angular::Isotropic)
        .map(|a| a.at(e));
    let ng = groups.len();
    let mut out = vec![vec![0.0; ng]; order + 1];
    let Some((lo, hi)) = two_body_range(e, awr, q) else {
        return out;
    };
    let first = usize::from(angular.is_none());
    if angular.is_none() {
        uniform_fractions(groups, lo, hi, scale, &mut out[0]);
        if order == 0 {
            return out;
        }
    }
    let bounds = groups.bounds();
    let group = |x: f64| {
        groups
            .find(x)
            .unwrap_or(if x > bounds[0] { 0 } else { ng - 1 })
    };
    if hi <= lo {
        // しきい値ちょうど: 中性子は重心とともに前へ進む（μ = 1、P_l = 1）。角度分布があるときは
        // P0 もここで入れる（一様な振り分けをしていないため）
        let g = group(lo);
        for row in out.iter_mut().skip(first) {
            row[g] += scale;
        }
        return out;
    }
    let ab = awr * (1.0 + (awr + 1.0) * q / (awr * e)).sqrt();
    let c = 1.0 - ab * ab;
    let speed = |x: f64| (awr + 1.0) * (x / e).sqrt();
    // 出射エネルギーの区切り: 群境界と、角度分布の折れる μ_c（E' = E (1 + A²β² + 2Aβ μ_c)/(A+1)²）
    let kinks = angular
        .as_ref()
        .map(Distribution::kinks)
        .unwrap_or_default();
    let energy = |mu_c: f64| e * (1.0 + ab * ab + 2.0 * ab * mu_c) / (awr + 1.0).powi(2);
    let mut cuts: Vec<f64> = bounds
        .iter()
        .copied()
        .chain(kinks.into_iter().map(energy))
        .filter(|&b| b > lo && b < hi)
        .chain([lo, hi])
        .collect();
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    // 角度分布は μ_c の高次の多項式や折れ線になりうるので、範囲を細かく分ける
    // Legendre の分布は次数が高いほど細かく分ける（表の分布は折れ目で区切ってあるので少しでよい）
    let pieces = angular
        .as_ref()
        .map_or(1, |a| (a.order() / 2).max(MIN_ANGULAR_PIECES));
    for w in cuts.windows(2) {
        let g = group(0.5 * (w[0] + w[1]));
        let (sa, sb) = (speed(w[0]), speed(w[1]));
        let width = (sb - sa) / pieces as f64;
        let half = 0.5 * width;
        for piece in 0..pieces {
            let start = sa + width * piece as f64;
            for &(x, weight) in &GAUSS_LEGENDRE_8 {
                let s = start + half * (x + 1.0);
                let mu = (s * s + c) / (2.0 * s);
                let density = angular.as_ref().map_or(0.5, |a| {
                    let mu_c = ((s * s - 1.0 - ab * ab) / (2.0 * ab)).clamp(-1.0, 1.0);
                    a.density(mu_c)
                });
                let factor = scale * weight * half * s / ab * density;
                for (l, p) in legendre(mu, order).into_iter().enumerate().skip(first) {
                    out[l][g] += factor * p;
                }
            }
        }
    }
    if angular.is_some() {
        let total: f64 = out[0].iter().sum();
        if total > 0.0 {
            for v in out.iter_mut().flatten() {
                *v *= scale / total;
            }
        }
    }
    out
}

/// 角度分布のある2体反応で、区切った出射の速さの範囲をさらに等分する最少の数。Legendre の分布は次数の
/// 半分まで増やす（8 点の Gauss–Legendre は 15 次まで正確で、密度は μ_c の NL 次、μ_c は s の2次）
const MIN_ANGULAR_PIECES: usize = 4;

/// 出射エネルギーの範囲の端が群境界を横切る入射エネルギーと、しきい値。
/// 被積分関数はここで折れ曲がるので、積分区間の境界に入れる。
fn two_body_kinks(groups: &GroupStructure, awr: f64, q: f64) -> Vec<f64> {
    let bounds = groups.bounds();
    let (lowest, highest) = (bounds[groups.len()], bounds[0]);
    let threshold = if q < 0.0 { -q * (awr + 1.0) / awr } else { 0.0 };
    let start = threshold.max(lowest);
    if start >= highest {
        return vec![];
    }
    let mut kinks = vec![start];
    // 範囲の端は E について単調とは限らない（しきい値の近く）ので、標本点で符号の変化を探す
    const SAMPLES: usize = 400;
    let energy = |i: usize| start * (highest / start).powf(i as f64 / SAMPLES as f64);
    let ends = |e: f64| two_body_range(e, awr, q).unwrap_or((e, e));
    for &b in bounds {
        for end in [|r: (f64, f64)| r.0, |r: (f64, f64)| r.1] {
            let f = |e: f64| end(ends(e)) - b;
            for i in 0..SAMPLES {
                let (mut a, mut c) = (energy(i), energy(i + 1));
                if f(a).signum() == f(c).signum() {
                    continue;
                }
                for _ in 0..60 {
                    let m = 0.5 * (a + c);
                    if f(m).signum() == f(a).signum() {
                        a = m;
                    } else {
                        c = m;
                    }
                }
                kinks.push(0.5 * (a + c));
            }
        }
    }
    kinks
}

/// MF6 の中性子放出の ∫_g σ(E) y(E) w(E) P(E→g') dE。
fn emission_rows(
    groups: &GroupStructure,
    weight: &Weight,
    table: &Tab1,
    data: &Mf6,
) -> Vec<Vec<f64>> {
    let ng = groups.len();
    let mut rows = vec![vec![0.0; ng]; ng];
    for emission in &data.neutrons {
        let grouped = match emission.frame {
            Frame::Laboratory => GroupedOutgoing::new(groups, &emission.distributions),
            Frame::CenterOfMass => {
                GroupedOutgoing::center_of_mass(groups, &emission.distributions, data.awr)
            }
        };
        let mut breakpoints = table.x.clone();
        breakpoints.extend(&emission.multiplicity.x);
        breakpoints.extend(grouped.incident());
        let part = group_integrals_multi(groups, weight, &breakpoints, ng, |e, out| {
            let source =
                table.eval(e).unwrap_or(0.0) * emission.multiplicity.eval(e).unwrap_or(0.0);
            if source != 0.0 {
                grouped.accumulate(e, source, out);
            }
        });
        for (sum, row) in rows.iter_mut().zip(part) {
            for (s, v) in sum.iter_mut().zip(row) {
                *s += v;
            }
        }
    }
    rows
}

/// MF5 の出射エネルギー分布と多重度 `count`、実験室系の角度分布（MF4。エネルギーと独立）による
/// ∫_g count σ(E) w(E) a_l(E) Σ_k p_k(E) P_k(E→g') dE（[l][g][g']）。MF5 の分布は実験室系。
/// 角度分布がなければ P1 以上は 0。
fn spectrum_moments(
    groups: &GroupStructure,
    weight: &Weight,
    table: &Tab1,
    data: &Mf5,
    count: f64,
    angular: Option<&Mf4>,
    order: usize,
) -> Vec<Vec<Vec<f64>>> {
    let ng = groups.len();
    let angular = angular.filter(|a| a.angular != Angular::Isotropic && order > 0);
    let parts: Vec<(&Tab1, GroupedLaw)> = data
        .subsections
        .iter()
        .map(|(p, law)| (p, GroupedLaw::new(groups, law)))
        .collect();
    let mut breakpoints = table.x.clone();
    for (p, grouped) in &parts {
        breakpoints.extend(&p.x);
        breakpoints.extend(grouped.breakpoints());
    }
    breakpoints.extend(angular.map(Mf4::energies).unwrap_or_default());
    let rows = group_integrals_multi(groups, weight, &breakpoints, ng * (order + 1), |e, out| {
        let source = count * table.eval(e).unwrap_or(0.0);
        if source == 0.0 {
            return;
        }
        let (p0, higher) = out.split_at_mut(ng);
        for (p, grouped) in &parts {
            let pk = p.eval(e).unwrap_or(0.0);
            if pk != 0.0 {
                grouped.accumulate(e, source * pk, p0);
            }
        }
        if let Some(mf4) = angular {
            let a = mf4.moments(e, order);
            for (l, row) in higher.chunks_exact_mut(ng).enumerate() {
                for (o, v) in row.iter_mut().zip(p0.iter()) {
                    *o += a[l + 1] * v;
                }
            }
        }
    });
    by_order(rows, ng, order)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::{AngleTables, Interpolation, LegendreTable};

    /// ルジャンドルの次数ごとの、全群の和
    fn sums(moments: &[Vec<f64>]) -> Vec<f64> {
        moments.iter().map(|row| row.iter().sum()).collect()
    }

    #[test]
    fn hydrogen_moments_follow_the_kinematics() {
        // A = 1 の弾性散乱は μ_lab = √(E'/E) で E' は [0, E] に一様: ⟨P1⟩ = 2/3, ⟨P2⟩ = 1/4, ⟨P3⟩ = 0
        let groups = GroupStructure::new(vec![1.0e7, 1.0e4, 3.0e3, 1.0e3, 1.0, 1.0e-5]).unwrap();
        let m = two_body_moments(&groups, 5.0e3, 1.0, 0.0, 1.0, 3, None);
        let total = sums(&m);
        for (got, expected) in total.iter().zip([1.0, 2.0 / 3.0, 0.25, 0.0]) {
            assert!((got - expected).abs() < 1e-12, "{total:?}");
        }
        // P0 は一様な振り分けと同じ
        let mut uniform = vec![0.0; groups.len()];
        let (lo, hi) = two_body_range(5.0e3, 1.0, 0.0).unwrap();
        uniform_fractions(&groups, lo, hi, 1.0, &mut uniform);
        for (a, b) in m[0].iter().zip(&uniform) {
            assert!((a - b).abs() < 1e-14, "{:?} vs {uniform:?}", m[0]);
        }
    }

    #[test]
    fn heavy_target_moments_follow_the_kinematics() {
        // 重心系で等方な弾性散乱の実験室系の平均の余弦は 2/(3A)
        let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 9.0e2, 8.5e2, 1.0e-5]).unwrap();
        let m = two_body_moments(&groups, 1.0e3, 10.0, 0.0, 2.0, 1, None);
        let total = sums(&m);
        assert!((total[0] - 2.0).abs() < 1e-12);
        assert!(
            (total[1] / total[0] - 2.0 / 30.0).abs() < 1e-12,
            "{total:?}"
        );
    }

    /// 重心系の分布 f(μ_c) = 1/2 + 3/2 a μ_c（全エネルギーで同じ）
    fn linear_cm(a: f64) -> Mf4 {
        Mf4 {
            frame: Frame::CenterOfMass,
            angular: Angular::Legendre(LegendreTable {
                breaks: vec![2],
                laws: vec![Interpolation::LinLin],
                energies: vec![1.0e-5, 2.0e7],
                coefficients: vec![vec![a], vec![a]],
            }),
        }
    }

    #[test]
    fn isotropic_angular_data_changes_nothing() {
        let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 9.0e2, 8.5e2, 1.0e-5]).unwrap();
        let isotropic = Mf4 {
            frame: Frame::CenterOfMass,
            angular: Angular::Isotropic,
        };
        let plain = two_body_moments(&groups, 1.0e3, 10.0, 0.0, 2.0, 3, None);
        let with = two_body_moments(&groups, 1.0e3, 10.0, 0.0, 2.0, 3, Some(&isotropic));
        for (a, b) in plain.iter().flatten().zip(with.iter().flatten()) {
            assert!((a - b).abs() < 1e-12, "{plain:?} vs {with:?}");
        }
    }

    #[test]
    fn threshold_with_angles_keeps_the_particles() {
        // しきい値ちょうどでは出射エネルギーの幅が 0。角度分布があっても P0 の和は scale
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
        let (awr, q) = (10.0_f64, -1.0e6_f64);
        let threshold = -q * (awr + 1.0) / awr;
        let m = two_body_moments(&groups, threshold, awr, q, 2.0, 1, Some(&linear_cm(0.2)));
        let total = sums(&m);
        assert!((total[0] - 2.0).abs() < 1e-12, "{m:?}");
    }

    #[test]
    fn center_of_mass_distribution_shifts_the_outgoing_energies() {
        // A = 10、E = 1e3 eV の弾性散乱。E' = E (1 + A² + 2A μ_c)/(A+1)² は μ_c に線形なので、
        // 群 [E1, E2] に出る確率は ∫ f dμ_c（f = 1/2 + 3/2 a μ_c）
        let (awr, e, a) = (10.0_f64, 1.0e3_f64, 0.2);
        let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 9.0e2, 8.5e2, 1.0e-5]).unwrap();
        let mf4 = linear_cm(a);
        let m = two_body_moments(&groups, e, awr, 0.0, 1.0, 1, Some(&mf4));
        let mu_c = |ep: f64| ((awr + 1.0).powi(2) * ep / e - 1.0 - awr * awr) / (2.0 * awr);
        let cdf = |mu: f64| {
            let mu = mu.clamp(-1.0, 1.0);
            0.5 * (mu + 1.0) + 0.75 * a * (mu * mu - 1.0)
        };
        let b = groups.bounds();
        for g in 0..groups.len() {
            let expected = cdf(mu_c(b[g])) - cdf(mu_c(b[g + 1]));
            assert!((m[0][g] - expected).abs() < 1e-12, "group {g}: {:?}", m[0]);
        }
        // P1 の和は ∫ f(μ_c) μ_lab(μ_c) dμ_c、μ_lab = (1 + A μ_c)/√(1 + A² + 2A μ_c)
        let nodes = xlc_numeric::quadrature::gauss_legendre(64).unwrap();
        let expected: f64 = nodes
            .iter()
            .map(|&(x, w)| {
                let f = 0.5 + 1.5 * a * x;
                w * f * (1.0 + awr * x) / (1.0 + awr * awr + 2.0 * awr * x).sqrt()
            })
            .sum();
        let total = sums(&m);
        assert!(
            (total[1] - expected).abs() < 1e-12,
            "{total:?} vs {expected}"
        );
    }

    use crate::endf::{ClampedTab1, IncidentLaw};

    fn constant(value: f64) -> Tab1 {
        Tab1::linear(vec![1.0e-5, 2.0e7], vec![value, value])
    }

    #[test]
    fn elastic_ranges() {
        assert_eq!(two_body_range(10.0, 1.0, 0.0), Some((0.0, 10.0)));
        let (lo, hi) = two_body_range(10.0, 238.0, 0.0).unwrap();
        let alpha = (237.0_f64 / 239.0).powi(2);
        assert!((lo - alpha * 10.0).abs() < 1e-12 && (hi - 10.0).abs() < 1e-12);
    }

    #[test]
    fn two_body_range_holds_at_the_rounded_threshold() {
        // しきい値 −Q (A + 1)/A をそのまま計算した入射エネルギーでも、丸めで β² が負（−2.2e-16）になる組がある。
        // しきい値の断面積の中性子を失わないよう、しきい値以上では β² = 0 として出射エネルギーを返す
        let (awr, q): (f64, f64) = (23.885511445655023, -284446.41774354107);
        let threshold = -q * (awr + 1.0) / awr;
        assert!(1.0 + (awr + 1.0) * q / (awr * threshold) < 0.0);
        let expected = threshold / (awr + 1.0).powi(2);
        assert_eq!(
            two_body_range(threshold, awr, q),
            Some((expected, expected))
        );
    }

    #[test]
    fn inelastic_threshold() {
        // しきい値 -Q (A+1)/A
        let (awr, q) = (236.0, -4.49e4);
        let threshold = -q * (awr + 1.0) / awr;
        assert_eq!(two_body_range(threshold * 0.999, awr, q), None);
        let (lo, hi) = two_body_range(threshold * 1.5, awr, q).unwrap();
        assert!(lo < hi && hi < threshold * 1.5);
    }

    #[test]
    fn hydrogen_elastic_matrix() {
        // A=1 では出射エネルギーは [0, E] で一様。重み 1/E のもとで群内に残る割合は解析的に求まる
        let groups = GroupStructure::new(vec![1.0e4, 1.0e3, 1.0e2]).unwrap();
        let tables = BTreeMap::from([(2, constant(20.0))]);
        let m = scattering_matrix(
            &groups,
            &Weight::lwr(),
            1.0,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &BTreeMap::new(),
                angular: &BTreeMap::new(),
            },
            None,
        )
        .unwrap();
        let stay = 20.0 * (1.0 - 0.9 / 10.0_f64.ln());
        assert!((m[0][0] - stay).abs() < 1e-9 * stay, "{m:?}");
        // 群構造の下に出る分も最下位の群に入るので、行の和は σ
        for row in &m {
            assert!((row.iter().sum::<f64>() - 20.0).abs() < 1e-9);
        }
    }

    #[test]
    fn hydrogen_moments_keep_the_mean_cosine() {
        // 重心系で等方な A = 1 の弾性散乱は、入射エネルギーによらず実験室系の平均の余弦が 2/3 なので、
        // どの入射の群でも P1 の行の和は P0 の行の和の 2/3。P0 は scattering_matrix と同じ
        let groups = GroupStructure::new(vec![1.0e4, 1.0e3, 1.0e2]).unwrap();
        let tables = BTreeMap::from([(2, constant(20.0))]);
        let args = (&BTreeMap::new(), &BTreeMap::new());
        let moments = scattering_moments(
            &groups,
            &Weight::lwr(),
            1.0,
            &Reactions {
                tables: &tables,
                q_values: args.0,
                emissions: args.1,
                angular: &BTreeMap::new(),
            },
            None,
            3,
        )
        .unwrap();
        assert_eq!(moments.len(), 4);
        let p0 = scattering_matrix(
            &groups,
            &Weight::lwr(),
            1.0,
            &Reactions {
                tables: &tables,
                q_values: args.0,
                emissions: args.1,
                angular: &BTreeMap::new(),
            },
            None,
        )
        .unwrap();
        assert_eq!(moments[0], p0);
        for (row0, row1) in moments[0].iter().zip(&moments[1]) {
            let ratio = row1.iter().sum::<f64>() / row0.iter().sum::<f64>();
            assert!((ratio - 2.0 / 3.0).abs() < 1e-10, "{ratio}");
        }
    }

    #[test]
    fn discrete_inelastic_uses_q_value() {
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
        let (awr, q) = (236.0, -4.49e4);
        let tables = BTreeMap::from([(2, constant(10.0)), (51, constant(1.0))]);
        let q_values = BTreeMap::from([(51, q)]);
        let m = scattering_matrix(
            &groups,
            &Weight::lwr(),
            awr,
            &Reactions {
                tables: &tables,
                q_values: &q_values,
                emissions: &BTreeMap::new(),
                angular: &BTreeMap::new(),
            },
            None,
        )
        .unwrap();
        // 離散準位の分は高い群から低い群へも移る。行の和は σel + σ51（しきい値より上）
        assert!(m[0][1] > 0.0);
        assert!((m[0].iter().sum::<f64>() - 11.0).abs() < 1e-9);
    }

    #[test]
    fn two_body_mf6_level_uses_q_value() {
        // MF6 の LAW=2（2体反応）で与えられた準位は、MF6 がないときと同じく Q 値の2体反応
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
        let (awr, q) = (118.87, -1.17e6);
        let tables = BTreeMap::from([(2, constant(10.0)), (51, constant(1.0))]);
        let q_values = BTreeMap::from([(51, q)]);
        let two_body = Emission::Mf6(Mf6 {
            awr,
            neutrons: vec![],
            two_body: true,
        });
        let matrix = |emissions: &BTreeMap<u32, Emission>| {
            scattering_matrix(
                &groups,
                &Weight::lwr(),
                awr,
                &Reactions {
                    tables: &tables,
                    q_values: &q_values,
                    emissions,
                    angular: &BTreeMap::new(),
                },
                None,
            )
        };
        let with_mf6 = matrix(&BTreeMap::from([(51, two_body.clone())])).unwrap();
        assert_eq!(with_mf6, matrix(&BTreeMap::new()).unwrap());
        // 準位以外の反応が LAW=2 なら、Q 値の2体反応とは扱えない
        let tables = BTreeMap::from([(2, constant(10.0)), (16, constant(1.0))]);
        let err = scattering_matrix(
            &groups,
            &Weight::lwr(),
            awr,
            &Reactions {
                tables: &tables,
                q_values: &q_values,
                emissions: &BTreeMap::from([(16, two_body)]),
                angular: &BTreeMap::new(),
            },
            None,
        )
        .unwrap_err();
        assert!(matches!(err, ScatterError::TwoBody(16)), "{err}");
    }

    #[test]
    fn total_inelastic_without_partials_needs_a_distribution() {
        // MT51〜91 がなければ、MT4 自身を数える（分布がなければエラー）
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6]).unwrap();
        let tables = BTreeMap::from([(4, constant(1.0))]);
        let err = scattering_matrix(
            &groups,
            &Weight::lwr(),
            15.86,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &BTreeMap::new(),
                angular: &BTreeMap::new(),
            },
            None,
        )
        .unwrap_err();
        assert!(matches!(err, ScatterError::NoDistribution(4)), "{err}");
    }

    #[test]
    fn reactions_without_neutrons_need_no_distribution() {
        // MT155 (z,tα) などは MT152〜200 の範囲でも中性子を出さない
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6]).unwrap();
        for mt in [155, 182, 191, 192, 193, 197] {
            let tables = BTreeMap::from([(mt, constant(1.0))]);
            let m = scattering_matrix(
                &groups,
                &Weight::lwr(),
                15.86,
                &Reactions {
                    tables: &tables,
                    q_values: &BTreeMap::new(),
                    emissions: &BTreeMap::new(),
                    angular: &BTreeMap::new(),
                },
                None,
            )
            .unwrap();
            assert_eq!(m[0][0], 0.0, "MT{mt}");
        }
    }

    #[test]
    fn free_gas_matrix_satisfies_detailed_balance() {
        // 重み関数の熱領域（0.1 eV 以下）は kT = 0.0253 eV の Maxwell 分布なので、
        // 同じ温度の自由ガスなら群間でも φ_g σ(g→g') = φ_g' σ(g'→g) が成り立つ。
        // 自由原子の断面積 4 b を温度に拡げた 4 N(E) を与える（表の補間誤差は 1e-7 以下）
        // 熱領域だけを使うため WIMS 69群の 0.2 eV 以下に絞る。群構造の外へ散乱した分は端の群に入る（free_gas.rs）ので、最上位の群はバランスの対象にならず、対象の群（0.1 eV 以下）より上に1群以上を残す
        let groups = GroupStructure::new(
            GroupStructure::wims69()
                .bounds()
                .iter()
                .copied()
                .filter(|&e| e <= 0.2)
                .collect(),
        )
        .unwrap();
        let weight = Weight::lwr();
        let tables = BTreeMap::from([(2, broadened_constant(4.0, 15.86, 0.0253))]);
        let thermal = FreeGas {
            kt: 0.0253,
            cutoff: 4.0,
        };
        let m = scattering_matrix(
            &groups,
            &weight,
            15.86,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &BTreeMap::new(),
                angular: &BTreeMap::new(),
            },
            Some(&Thermal::FreeGas(thermal)),
        )
        .unwrap();
        let flux = group_integrals(&groups, &weight, &[], |_| 1.0);
        // 最下位の群は除く。1e-5 eV より下へ散乱した分はこの群に入れる（中性子を失わない）が、
        // 1e-5 eV より下から入射する中性子は積分しないので、この群だけはバランスが厳密でない
        let maxwellian: Vec<usize> = (0..groups.len() - 1)
            .filter(|&g| groups.range(g).1 <= 0.1)
            .collect();
        assert!(maxwellian.len() > 10);
        for &g in &maxwellian {
            // 標的の熱運動で、実効的な散乱断面積は自由原子の値より大きい
            let row: f64 = m[g].iter().sum();
            assert!(row > 4.0 && row < 6.0, "{row}");
            for &h in &maxwellian {
                let (forward, backward) = (flux[g] * m[g][h], flux[h] * m[h][g]);
                assert!(
                    (forward - backward).abs() <= 1e-6 * forward.max(backward) + 1e-12,
                    "{g}→{h}: {forward} vs {backward}"
                );
            }
        }
        // 熱群から上方散乱がある
        let g = groups.find(0.0253).unwrap();
        assert!(m[g][g - 1] > 0.0);
    }

    fn sab_matrix(elastic: f64) -> (GroupStructure, Vec<Vec<f64>>, f64) {
        use crate::multigroup::sab_tests::coarse_free_gas_table;
        let groups = GroupStructure::wims69();
        // 重み関数の Maxwell 分布（kT = 0.0253 eV）と同じ温度にする（自由ガスの S は α, β で書けば温度によらない）
        let temperature = 0.0253 / crate::physics::BK;
        let mut data = coarse_free_gas_table();
        data.temperatures = vec![temperature];
        let cutoff = 1.0;
        let sab = ThermalSab::new(
            &data,
            temperature,
            cutoff,
            crate::multigroup::SabOutside::Zero,
        )
        .unwrap();
        let tables = BTreeMap::from([(2, constant(elastic))]);
        let m = scattering_matrix(
            &groups,
            &Weight::lwr(),
            1.0,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &BTreeMap::new(),
                angular: &BTreeMap::new(),
            },
            Some(&Thermal::Sab {
                sab,
                free_gas: FreeGas { kt: 0.0253, cutoff },
            }),
        )
        .unwrap();
        (groups, m, cutoff)
    }

    /// A = 1 の S(α, β) を 1 eV 未満で使い、その上 `free_gas_upper` (eV) 未満を自由ガスで扱う、一定の弾性散乱断面積
    /// 20 b の散乱行列
    fn sab_and_free_gas(free_gas_upper: f64) -> Result<Vec<Vec<f64>>, ScatterError> {
        use crate::multigroup::sab_tests::coarse_free_gas_table;
        let kt = 0.0253;
        let temperature = kt / crate::physics::BK;
        let mut data = coarse_free_gas_table();
        data.temperatures = vec![temperature];
        data.effective_temperature = Tab1::linear(vec![temperature], vec![temperature]);
        let sab = ThermalSab::new(
            &data,
            temperature,
            1.0,
            crate::multigroup::SabOutside::ShortCollisionTime,
        )
        .unwrap();
        let tables = BTreeMap::from([(2, constant(20.0))]);
        scattering_matrix(
            &GroupStructure::wims69(),
            &Weight::lwr(),
            1.0,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &BTreeMap::new(),
                angular: &BTreeMap::new(),
            },
            Some(&Thermal::Sab {
                sab,
                free_gas: FreeGas {
                    kt,
                    cutoff: free_gas_upper,
                },
            }),
        )
    }

    #[test]
    fn free_gas_continues_above_the_sab_cutoff() {
        // S(α, β) の上限から自由ガスの上限までは、標的の熱運動で上方散乱があり、行の和は MF3 の弾性散乱断面積に
        // なる。その上は静止した標的なので上方散乱がない
        let groups = GroupStructure::wims69();
        let m = sab_and_free_gas(4.0).unwrap();
        let weight = Weight::lwr();
        let table = constant(20.0);
        let flux = group_integrals(&groups, &weight, &table.x, |_| 1.0);
        let sigma = group_integrals(&groups, &weight, &table.x, |e| table.eval(e).unwrap());
        let (mut free_gas, mut at_rest) = (0, 0);
        for g in 1..groups.len() {
            let (low, high) = groups.range(g);
            let up: f64 = m[g][..g].iter().sum();
            if low >= 1.0 && high <= 4.0 {
                free_gas += 1;
                assert!(up > 1e-4 * sigma[g] / flux[g], "group {g}: {up}");
                let row: f64 = m[g].iter().sum();
                let expected = sigma[g] / flux[g];
                assert!(
                    (row / expected - 1.0).abs() < 1e-9,
                    "group {g}: {row} vs {expected}"
                );
            } else if low >= 4.0 {
                at_rest += 1;
                assert_eq!(up, 0.0, "group {g}");
            }
        }
        assert!(free_gas >= 3 && at_rest >= 3);
        // 自由ガスの上限が S(α, β) の上限より低いのはエラー
        assert!(matches!(
            sab_and_free_gas(0.5),
            Err(ScatterError::FreeGasBelowSab { .. })
        ));
    }

    #[test]
    fn sab_matrix_satisfies_detailed_balance() {
        // 重み関数の熱領域（0.1 eV 以下）は kT = 0.0253 eV の Maxwell 分布で、S(α, β) も同じ温度なので、
        // 群間でも φ_g σ(g→g') = φ_g' σ(g'→g) が成り立つ（表の細かさによらない）。
        // 自由ガスとの物理的な一致は、核と群への積分の段階で確かめている（sab のテスト）
        let (groups, m, _) = sab_matrix(1.0);
        let flux = group_integrals(&groups, &Weight::lwr(), &[], |_| 1.0);
        let maxwellian: Vec<usize> = (0..groups.len() - 1)
            .filter(|&g| groups.range(g).1 <= 0.1)
            .collect();
        assert!(maxwellian.len() > 10);
        for &g in &maxwellian {
            for &h in &maxwellian {
                let (forward, backward) = (flux[g] * m[g][h], flux[h] * m[h][g]);
                assert!(
                    (forward - backward).abs() <= 1e-6 * forward.max(backward) + 1e-12,
                    "{g}→{h}: {forward} vs {backward}"
                );
            }
        }
    }

    #[test]
    fn sab_ignores_mf3_elastic_below_cutoff() {
        // S(α, β) の核は断面積そのものなので、cutoff 未満の行は MF3 の弾性散乱断面積によらない
        let (groups, one, cutoff) = sab_matrix(1.0);
        let (_, ten, _) = sab_matrix(10.0);
        for g in (0..groups.len()).filter(|&g| groups.range(g).1 <= cutoff) {
            assert_eq!(one[g], ten[g], "group {g}");
        }
    }

    #[test]
    fn elastic_uses_free_gas_even_with_mf6_data() {
        // MF6 に MT2 があっても、弾性散乱は2体反応と自由ガスで扱う
        let groups = GroupStructure::new(vec![1.0, 0.01, 1.0e-4]).unwrap();
        let tables = BTreeMap::from([(2, constant(4.0))]);
        let mf6 = BTreeMap::from([(
            2,
            Emission::Mf6(Mf6 {
                awr: 1.0,
                neutrons: vec![],
                two_body: false,
            }),
        )]);
        let thermal = FreeGas {
            kt: 0.0253,
            cutoff: 4.0,
        };
        let m = scattering_matrix(
            &groups,
            &Weight::lwr(),
            1.0,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &mf6,
                angular: &BTreeMap::new(),
            },
            Some(&Thermal::FreeGas(thermal)),
        )
        .unwrap();
        // 2体反応では上方散乱しない
        assert!(m[1][0] > 0.0, "{m:?}");
    }

    #[test]
    fn invalid_free_gas_is_error() {
        let groups = GroupStructure::new(vec![1.0, 0.01]).unwrap();
        let tables = BTreeMap::from([(2, constant(4.0))]);
        for (kt, cutoff) in [(0.0, 4.0), (-0.01, 4.0), (0.0253, f64::NAN)] {
            let thermal = FreeGas { kt, cutoff };
            let result = scattering_matrix(
                &groups,
                &Weight::lwr(),
                15.86,
                &Reactions {
                    tables: &tables,
                    q_values: &BTreeMap::new(),
                    emissions: &BTreeMap::new(),
                    angular: &BTreeMap::new(),
                },
                Some(&Thermal::FreeGas(thermal)),
            );
            assert!(result.is_err(), "kt {kt}, cutoff {cutoff}");
        }
    }

    #[test]
    fn neutron_emission_without_distribution_is_error() {
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6]).unwrap();
        let tables = BTreeMap::from([(16, constant(1.0))]);
        let err = scattering_matrix(
            &groups,
            &Weight::lwr(),
            15.86,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &BTreeMap::new(),
                angular: &BTreeMap::new(),
            },
            None,
        )
        .unwrap_err();
        assert!(matches!(err, ScatterError::NoDistribution(16)), "{err}");
    }

    #[test]
    fn mf6_multiplicity_scales_the_row() {
        use crate::endf::{Frame, Interpolation, NeutronEmission, OutgoingDistributions};
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
        let emission = NeutronEmission {
            multiplicity: constant(2.0),
            frame: Frame::Laboratory,
            distributions: OutgoingDistributions {
                breaks: vec![1],
                laws: vec![IncidentLaw {
                    law: Interpolation::LinLin,
                    unit_base: false,
                }],
                incident: vec![1.0e-5],
                spectra: vec![Tab1::linear(vec![1.0e-5, 1.0e6], vec![1.0e-6, 1.0e-6])],
            },
        };
        let mf6 = BTreeMap::from([(
            16,
            Emission::Mf6(Mf6 {
                awr: 15.86,
                neutrons: vec![emission],
                two_body: false,
            }),
        )]);
        let tables = BTreeMap::from([(16, constant(0.5))]);
        let m = scattering_matrix(
            &groups,
            &Weight::lwr(),
            15.86,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &mf6,
                angular: &BTreeMap::new(),
            },
            None,
        )
        .unwrap();
        // 多重度 2 × σ 0.5 がすべて低い群へ
        assert!(
            (m[0][1] - 1.0).abs() < 1e-6 && m[0][0].abs() < 1e-12,
            "{m:?}"
        );
    }

    /// θ = 100 keV、U = 1 MeV の蒸発スペクトルの MF5
    fn evaporation() -> Mf5 {
        use crate::endf::EnergyLaw;
        Mf5 {
            subsections: vec![(
                constant(1.0),
                EnergyLaw::Evaporation {
                    restriction: 1.0e6,
                    theta: ClampedTab1::new(constant(1.0e5)).unwrap(),
                },
            )],
        }
    }

    #[test]
    fn mf5_rows_use_the_multiplicity_of_the_reaction() {
        // (n,2n) は2個、(n,nα) は1個の中性子を出す。行の和は 多重度 × σ（しきい値より上）
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
        for (mt, count) in [(16, 2.0), (22, 1.0)] {
            let tables = BTreeMap::from([(mt, constant(0.5))]);
            let emissions = BTreeMap::from([(mt, Emission::Mf5(evaporation()))]);
            let m = scattering_matrix(
                &groups,
                &Weight::lwr(),
                15.86,
                &Reactions {
                    tables: &tables,
                    q_values: &BTreeMap::new(),
                    emissions: &emissions,
                    angular: &BTreeMap::new(),
                },
                None,
            )
            .unwrap();
            let row: f64 = m[0].iter().sum();
            assert!((row - count * 0.5).abs() < 1e-9, "MT{mt}: {m:?}");
            assert_eq!(m[1], vec![0.0, 0.0], "MT{mt}");
        }
    }

    #[test]
    fn multiplicities_follow_the_reaction_definitions() {
        for (mt, count) in [
            (11, 2.0),
            (16, 2.0),
            (17, 3.0),
            (22, 1.0),
            (24, 2.0),
            (25, 3.0),
            (28, 1.0),
            (30, 2.0),
            (37, 4.0),
            (41, 2.0),
            (42, 3.0),
            (91, 1.0),
        ] {
            assert_eq!(multiplicity(mt), Some(count), "MT{mt}");
        }
        assert_eq!(multiplicity(5), None);
    }

    /// 実験室系の分布 f(μ) = 1/2 + 3/2 a μ（全エネルギーで同じ）
    fn linear_lab(a: f64) -> Mf4 {
        Mf4 {
            frame: Frame::Laboratory,
            ..linear_cm(a)
        }
    }

    #[test]
    fn continuum_lab_angles_scale_the_spectrum() {
        // MF5 の出射エネルギーと MF4 の実験室系の角度は独立なので、P1 の行は a_1 × P0 の行
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
        let tables = BTreeMap::from([(16, constant(0.5))]);
        let emissions = BTreeMap::from([(16, Emission::Mf5(evaporation()))]);
        let angular = BTreeMap::from([(16, linear_lab(0.3))]);
        let m = scattering_moments(
            &groups,
            &Weight::lwr(),
            15.86,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &emissions,
                angular: &angular,
            },
            None,
            1,
        )
        .unwrap();
        for (p1, p0) in m[1].iter().flatten().zip(m[0].iter().flatten()) {
            assert!((p1 - 0.3 * p0).abs() < 1e-12, "{m:?}");
        }
        assert!(m[0][0].iter().sum::<f64>() > 0.0);
    }

    #[test]
    fn elastic_angles_move_the_outgoing_energies() {
        // 重心系で前方に偏る弾性散乱は、等方より高い出射エネルギーの群に多く出る（行の和は同じ）
        let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 9.0e2, 8.5e2, 1.0e-5]).unwrap();
        let tables = BTreeMap::from([(2, constant(1.0))]);
        let solve = |angular: &BTreeMap<u32, Mf4>| {
            scattering_matrix(
                &groups,
                &Weight::lwr(),
                10.0,
                &Reactions {
                    tables: &tables,
                    q_values: &BTreeMap::new(),
                    emissions: &BTreeMap::new(),
                    angular,
                },
                None,
            )
            .unwrap()
        };
        let isotropic = solve(&BTreeMap::new());
        let forward = solve(&BTreeMap::from([(2, linear_cm(0.2))]));
        let g = 1;
        assert!(
            forward[g][g] > isotropic[g][g] * 1.01,
            "{forward:?} vs {isotropic:?}"
        );
        let sum = |m: &Vec<Vec<f64>>| m[g].iter().sum::<f64>();
        assert!((sum(&forward) - sum(&isotropic)).abs() < 1e-12 * sum(&isotropic));
    }

    #[test]
    fn center_of_mass_continuum_angles_are_error() {
        // MF5 の放出の角度分布は実験室系で使う（重心系の分布を黙って実験室系として使わない）
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
        let tables = BTreeMap::from([(16, constant(0.5))]);
        let emissions = BTreeMap::from([(16, Emission::Mf5(evaporation()))]);
        let angular = BTreeMap::from([(16, linear_cm(0.3))]);
        let reactions = Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &emissions,
            angular: &angular,
        };
        let err =
            scattering_moments(&groups, &Weight::lwr(), 15.86, &reactions, None, 1).unwrap_err();
        assert!(
            matches!(err, ScatterError::CenterOfMassContinuum(16)),
            "{err}"
        );
    }

    #[test]
    fn isotropic_center_of_mass_continuum_is_accepted() {
        // O-18 の MF5 の反応は、重心系で等方（LTT=0, LCT=2）の MF4 を持つ。P1 以上は 0 として扱う
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
        let tables = BTreeMap::from([(16, constant(0.5))]);
        let emissions = BTreeMap::from([(16, Emission::Mf5(evaporation()))]);
        let isotropic = Mf4 {
            frame: Frame::CenterOfMass,
            angular: Angular::Isotropic,
        };
        let angular = BTreeMap::from([(16, isotropic)]);
        let reactions = Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &emissions,
            angular: &angular,
        };
        let m = scattering_moments(&groups, &Weight::lwr(), 15.86, &reactions, None, 1).unwrap();
        assert!(m[1].iter().flatten().all(|&v| v == 0.0));
        assert!(m[0][0].iter().sum::<f64>() > 0.0);
    }

    #[test]
    fn tabulated_center_of_mass_kinks_are_integrated_exactly() {
        // 重心系の f(μ_c) が μ_c = 0 で折れる表: (−1, 0.25), (0, 0.25), (1, 1.25)。E' は μ_c に線形なので
        // 群に出る確率は累積分布 F(μ_c) の差
        let (awr, e) = (10.0_f64, 1.0e3_f64);
        let table = Tab1::linear(vec![-1.0, 0.0, 1.0], vec![0.25, 0.25, 1.25]);
        let mf4 = Mf4 {
            frame: Frame::CenterOfMass,
            angular: Angular::Tabulated(AngleTables {
                breaks: vec![2],
                laws: vec![Interpolation::LinLin],
                energies: vec![1.0e-5, 2.0e7],
                tables: vec![table.clone(), table],
            }),
        };
        let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 9.0e2, 8.5e2, 1.0e-5]).unwrap();
        let m = two_body_moments(&groups, e, awr, 0.0, 1.0, 0, Some(&mf4));
        let mu_c = |ep: f64| ((awr + 1.0).powi(2) * ep / e - 1.0 - awr * awr) / (2.0 * awr);
        let cdf = |mu: f64| {
            let mu = mu.clamp(-1.0, 1.0);
            if mu < 0.0 {
                0.25 * (mu + 1.0)
            } else {
                0.25 + 0.25 * mu + 0.5 * mu * mu
            }
        };
        let b = groups.bounds();
        for g in 0..groups.len() {
            let expected = cdf(mu_c(b[g])) - cdf(mu_c(b[g + 1]));
            assert!((m[0][g] - expected).abs() < 1e-12, "group {g}: {:?}", m[0]);
        }
    }

    #[test]
    fn lab_frame_two_body_angles_are_error() {
        let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 1.0e-5]).unwrap();
        let tables = BTreeMap::from([(2, constant(1.0))]);
        let err = scattering_matrix(
            &groups,
            &Weight::lwr(),
            10.0,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &BTreeMap::new(),
                angular: &BTreeMap::from([(2, linear_lab(0.2))]),
            },
            None,
        )
        .unwrap_err();
        assert!(matches!(err, ScatterError::LabTwoBody(2)), "{err}");
    }

    #[test]
    fn mf5_without_known_multiplicity_is_error() {
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6]).unwrap();
        let tables = BTreeMap::from([(5, constant(1.0))]);
        let emissions = BTreeMap::from([(5, Emission::Mf5(evaporation()))]);
        let err = scattering_matrix(
            &groups,
            &Weight::lwr(),
            15.86,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &emissions,
                angular: &BTreeMap::new(),
            },
            None,
        )
        .unwrap_err();
        assert!(matches!(err, ScatterError::Multiplicity(5)), "{err}");
    }

    /// 一定の自由原子の断面積 `sigma` を温度に拡げた σ N(E)。N(E) = (1 + 1/(2a²)) erf(a) + e^(-a²)/(a√π)、
    /// a² = A E/kT。対数で等間隔の 20000 点の表（線形補間の相対誤差は 1e-7 以下）。
    fn broadened_constant(sigma: f64, awr: f64, kt: f64) -> Tab1 {
        let n = |e: f64| {
            let a = (awr * e / kt).sqrt();
            (1.0 + 0.5 / (a * a)) * crate::special::erf(a)
                + (-a * a).exp() / (a * std::f64::consts::PI.sqrt())
        };
        let points = 20000;
        let (lo, hi) = (1.0e-5_f64, 2.0e7_f64);
        let x: Vec<f64> = (0..points)
            .map(|i| lo * (hi / lo).powf(i as f64 / (points - 1) as f64))
            .collect();
        let y = x.iter().map(|&e| sigma * n(e)).collect();
        Tab1::linear(x, y)
    }

    #[test]
    fn free_gas_rows_sum_to_the_broadened_elastic_cross_section() {
        // 温度に拡げた断面積に、規格化した出射分布を掛ける（NJOY THERMR と同じ）。熱領域の共鳴も
        // 拡げた形のまま行の和に入る。散乱はすべて群構造の中に残すので、行の和は群平均の断面積
        let groups = GroupStructure::wims69();
        let weight = Weight::lwr();
        let table = Tab1::linear(
            vec![1.0e-5, 0.2, 0.3, 0.4, 2.0e7],
            vec![4.0, 4.0, 50.0, 4.0, 4.0],
        );
        let thermal = FreeGas {
            kt: 0.0253,
            cutoff: 4.0,
        };
        let m = scattering_matrix(
            &groups,
            &weight,
            235.0,
            &Reactions {
                tables: &BTreeMap::from([(2, table.clone())]),
                q_values: &BTreeMap::new(),
                emissions: &BTreeMap::new(),
                angular: &BTreeMap::new(),
            },
            Some(&Thermal::FreeGas(thermal)),
        )
        .unwrap();
        let flux = group_integrals(&groups, &weight, &table.x, |_| 1.0);
        let sigma = group_integrals(&groups, &weight, &table.x, |e| table.eval(e).unwrap());
        for g in 0..groups.len() {
            let expected = sigma[g] / flux[g];
            let row: f64 = m[g].iter().sum();
            assert!(
                (row - expected).abs() <= 1e-9 * expected,
                "group {g}: {row} vs {expected}"
            );
        }
    }
}
