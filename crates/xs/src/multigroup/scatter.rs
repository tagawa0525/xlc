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
mod tests;
