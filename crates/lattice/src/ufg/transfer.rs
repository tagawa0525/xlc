//! UFG の範囲の群から出る弾性散乱の移行を、減速計算と同じ格子点の中性子束と散乱核で数える（ADR-0207）。
//!
//! リング i の核種 k の、入射の群 h から出射のレサジー x より下（エネルギーの低い側）への移行の率
//!
//!   B(h, x) = ∫_h σs(u') F_i(u') P_k(u' → u > x) du'
//!
//! と、x より上への移行の率 A(h, x) = ∫_h σs F P_k(u' → u < x) du' を、ライブラリの群の境界 x ごとに求め、群 h から下の群 g への
//! 移行を B(h, x_g) − B(h, x_{g+1})、上の群 g への移行を A(h, x_{g+1}) − A(h, x_g) とする（x_g は群 g の上端のレサジー）。最後の群の
//! 下限より下と、最初の群の上限より上への出射は、その群に入れる。自分の群への移行は、リングの中性子束で縮約した弾性散乱断面積から
//! ほかの群への移行を引いたもの（行の和が弾性散乱断面積になる）。P は減速計算の散乱核で、入射ごとに次の和:
//!
//! - 静止した標的の窓: x − Δ < u' < x で (e^{−(x−u')} − α)/(1 − α)、u' ≤ x − Δ で 0。σs は窓に使う格子点の値（自由ガスの核の
//!   入射では σs/N、S(α, β) の移行の入射では 0。[`super::Sweep::point`]）で、σs F は格子点の間で線形
//! - 自由ガスの核の差 D（自由ガスの核を使う入射）: [`Gas::shares`]（出射の x の下と上の積分）を、入射について区間ごとに 2 次式で
//!   近似したもの（[`Gas::share_cells`]）
//! - S(α, β) の移行（S(α, β) の上限より下の入射）: 表の累積 C(E', E_x − E')（[`SabTransfer::cumulative_at`]）
//!
//! 範囲の下限より下への出射も同じ式で数える（核は入射だけで決まり、その下に UFG の格子は要らない）。
//!
//! 間引いた格子（ADR-0185）でも、減速計算の格子点で数える。弾性散乱断面積は核種自身の格子で縮約した群定数なので、自分の群への
//! 移行には、ほかの群への移行を間引いた格子で数えた誤差だけが入る（ADR-0208）。

use std::ops::Range;

use rayon::prelude::*;
use xlc_format::library::{SabTransfer, Shielded};
use xlc_numeric::quadrature::{GAUSS_LEGENDRE_2, GAUSS_LEGENDRE_4};

use super::free_gas::{SHARE_CELLS, ShareCell, Thermal, grades};
use super::{Data, NuclideTransfer, UfgError, segment};

/// 自由ガスの核の差の割合を 2 次式で近似する入射の区間の、D の山が境界にかからないところでの幅の上限（レサジー）。割合は下端の側の山の
/// 全体の積分で、入射のエネルギーの冪のように滑らかに変わる
const SHARE_WIDTH: f64 = 0.01;

/// 移行の行の要素に許す負の値の大きさ（弾性散乱断面積に対する比）。差で決める要素の丸め（倍精度の 1e−16 の数十倍）より十分大きく、
/// 自由ガスの核の差の近似の誤差（`correction_shares_split_the_kernel_at_an_outgoing_energy` で 1e−10）よりも大きい
const NEGATIVE: f64 = 1.0e-9;

/// 減速計算の結果から移行を数える入力。
pub(super) struct Counting<'a> {
    pub(super) data: &'a [Data<'a>],
    pub(super) entries: &'a [Vec<(usize, f64)>],
    /// 格子点のレサジー（昇順）と（エネルギー、低い側の極限か）
    pub(super) u: &'a [f64],
    pub(super) points: &'a [(f64, bool)],
    /// [リング][格子点] の F
    pub(super) flux: &'a [Vec<f64>],
    /// [核種][格子点] の静止した標的の窓の σs
    pub(super) elastic: &'a [Vec<f64>],
    /// [リング][群] の ∫ F du
    pub(super) group_flux: &'a [Vec<f64>],
    /// [リング][核種] の縮約した群定数
    pub(super) nuclides: &'a [Vec<Shielded>],
    /// ライブラリの群の境界
    pub(super) bounds: &'a [f64],
    /// 減速計算した群
    pub(super) groups: Range<usize>,
    /// 範囲の上端 (eV)
    pub(super) high: f64,
    pub(super) thermal: Option<&'a Thermal>,
    /// [核種] の自由ガスの核種の番号（[`Thermal::gas`]）
    pub(super) gas_of: &'a [Option<usize>],
}

/// 群 h の [境界] の、下への移行の率 B と上への移行の率 A
struct Rates {
    below: Vec<Vec<f64>>,
    above: Vec<Vec<f64>>,
}

impl Counting<'_> {
    /// [リング][核種（物質の組成の順）] の移行の行
    pub(super) fn count(&self) -> Result<Vec<Vec<NuclideTransfer>>, UfgError> {
        let ng = self.bounds.len() - 1;
        let x: Vec<f64> = self.bounds.iter().map(|&b| (self.high / b).ln()).collect();
        let segments = self.segment_groups();
        // 自由ガスの核種ごと、境界ごとの、割合の 2 次式の区間（リングによらない）
        let gases = self.thermal.map_or(0, |_| {
            self.gas_of.iter().flatten().max().map_or(0, |&l| l + 1)
        });
        let pairs: Vec<(usize, usize)> = (0..gases)
            .flat_map(|l| (1..ng).map(move |b| (l, b)))
            .collect();
        let grades = grades(SHARE_CELLS);
        let cells: Vec<[Vec<ShareCell>; 2]> = pairs
            .par_iter()
            .map(|&(l, b)| self.share_cells(l, b, x[b], &grades))
            .collect();
        let tasks: Vec<(usize, usize)> = self
            .entries
            .iter()
            .enumerate()
            .flat_map(|(i, ring)| (0..ring.len()).map(move |e| (i, e)))
            .collect();
        let rows = tasks
            .par_iter()
            .map(|&(i, e)| {
                let k = self.entries[i][e].0;
                let mut rates = self.static_rates(i, k, &x, &segments);
                if let Some(l) = self.gas_of[k] {
                    let mine = pairs.iter().zip(&cells).filter(|((m, _), _)| *m == l);
                    for (&(_, b), [below, above]) in mine {
                        for cell in below {
                            self.add_cell(i, l, cell, &segments, |h, v| rates.below[h][b] += v);
                        }
                        for cell in above {
                            self.add_cell(i, l, cell, &segments, |h, v| rates.above[h][b] += v);
                        }
                    }
                    self.add_sab(i, k, l, &segments, &mut rates)?;
                }
                self.rows(i, e, &rates)
            })
            .collect::<Result<Vec<_>, UfgError>>()?;
        let mut rows = rows.into_iter();
        Ok(self
            .entries
            .iter()
            .map(|ring| rows.by_ref().take(ring.len()).collect())
            .collect())
    }

    /// [区間] の群（`groups` の中の番号）。減速計算と同じく、群の境界の高い側の格子点で次の群に移る
    fn segment_groups(&self) -> Vec<usize> {
        let nh = self.groups.len();
        let mut out = Vec::with_capacity(self.u.len().saturating_sub(1));
        let mut group = 0;
        for (m, &(e, below)) in self.points.iter().enumerate() {
            if m > 0 {
                out.push(group);
            }
            if group < nh && !below && e == self.bounds[self.groups.start + group + 1] {
                group += 1;
            }
        }
        out
    }

    /// 静止した標的の窓による、リング `i` の核種 `k` の下への移行の率
    fn static_rates(&self, i: usize, k: usize, x: &[f64], segments: &[usize]) -> Rates {
        let (top, nh) = (self.groups.start, self.groups.len());
        let ng = self.bounds.len() - 1;
        let d = &self.data[k];
        let alpha = 1.0 - 1.0 / d.gain;
        let (u, flux, elastic) = (self.u, &self.flux[i], &self.elastic[k]);
        let g = |p: usize| elastic[p] * flux[p];
        // 群の区間の始まりと終わり
        let first = |h: usize| segments.partition_point(|&s| s < h);
        // 区間 [lo, 群の終わり] の ∫ g e^{−(x−u')} du' と ∫ g du'
        let window = |h: usize, lo: f64, big_x: f64| -> (f64, f64) {
            let (start, end) = (first(h), first(h + 1));
            let from = u[start..end]
                .partition_point(|&v| v <= lo)
                .saturating_sub(1);
            let (mut decayed, mut plain) = (0.0, 0.0);
            for p in start + from..end {
                let (ua, uz) = (u[p], u[p + 1]);
                let a = ua.max(lo);
                if uz <= a {
                    continue;
                }
                let ga = g(p) + (g(p + 1) - g(p)) * (a - ua) / (uz - ua);
                let (ca, cz) = segment(a, uz, big_x);
                decayed += ca * ga + cz * g(p + 1);
                plain += 0.5 * (uz - a) * (ga + g(p + 1));
            }
            (decayed, plain)
        };
        let mut below = vec![vec![0.0; ng + 1]; nh];
        for (h, row) in below.iter_mut().enumerate() {
            let (xa, xz) = (x[top + h], x[top + h + 1]);
            let mut whole: Option<(f64, f64)> = None;
            for b in top + h + 1..ng {
                let lo = x[b] - d.width;
                // 群のすべての入射が窓の下端より上（u' + Δ ≤ x）なら、x より下へは出ない。b とともに x は増える
                if xz <= lo {
                    break;
                }
                let value = if xa >= lo {
                    let (decayed, plain) = *whole.get_or_insert_with(|| window(h, xa, xz));
                    (-(x[b] - xz)).exp() * decayed - alpha * plain
                } else {
                    let (decayed, plain) = window(h, lo, x[b]);
                    decayed - alpha * plain
                };
                row[b] = d.gain * value;
            }
        }
        Rates {
            below,
            above: vec![vec![0.0; ng + 1]; nh],
        }
    }

    /// 自由ガスの核種 `l` の、境界 `b`（レサジー `x`）の割合の区間（[`Gas::share_cells`]）
    fn share_cells(&self, l: usize, b: usize, x: f64, grades: &[f64]) -> [Vec<ShareCell>; 2] {
        let Some(thermal) = self.thermal else {
            return [Vec::new(), Vec::new()];
        };
        let (gas, ..) = thermal.gas(l);
        let last = self.u[self.u.len() - 1];
        let lethargy = |e: f64| (self.high / e).ln();
        // 自由ガスの核を使う入射: S(α, β) の上限以上、自由ガスの上限未満（[`Gas::window_scale`]）
        let incident = [
            lethargy(gas.cap().min(self.high)).max(0.0),
            gas.floor().map_or(last, |f| lethargy(f).min(last)),
        ];
        gas.share_cells(
            self.high,
            x,
            self.bounds[b],
            incident,
            grades,
            SHARE_CELLS,
            SHARE_WIDTH,
        )
    }

    /// 割合の区間 `cell` の ∫ σs F (割合) du' を、リング `i` の格子の区間ごとに 4 点の Gauss-Legendre 則（σs と F は格子点の間で
    /// 線形、割合は 2 次式で、積は 4 次なので厳密）で求め、区間の群ごとに `add` に渡す
    fn add_cell(
        &self,
        i: usize,
        l: usize,
        cell: &ShareCell,
        segments: &[usize],
        mut add: impl FnMut(usize, f64),
    ) {
        let Some(thermal) = self.thermal else { return };
        let (_, offset, sigma, _) = thermal.gas(l);
        let (u, flux) = (self.u, &self.flux[i]);
        let start = u
            .partition_point(|&v| v <= cell.a)
            .saturating_sub(1)
            .max(offset);
        for p in start..u.len() - 1 {
            let (ua, uz) = (u[p], u[p + 1]);
            if ua >= cell.b {
                break;
            }
            let (a, z) = (ua.max(cell.a), uz.min(cell.b));
            if z <= a {
                continue;
            }
            let (mid, half) = (0.5 * (a + z), 0.5 * (z - a));
            let value: f64 = GAUSS_LEGENDRE_4
                .iter()
                .map(|&(t, w)| {
                    let s = mid + half * t;
                    let r = (s - ua) / (uz - ua);
                    let lerp = |v0: f64, v1: f64| v0 + (v1 - v0) * r;
                    let sig = lerp(sigma[p - offset], sigma[p + 1 - offset]);
                    w * half * sig * lerp(flux[p], flux[p + 1]) * cell.quadratic.at(s)
                })
                .sum();
            add(segments[p], value);
        }
    }

    /// S(α, β) の移行で散乱する核種 `k`（自由ガスの核種 `l`）の、S(α, β) の上限より下の入射の移行の率を `rates` に足す
    fn add_sab(
        &self,
        i: usize,
        k: usize,
        l: usize,
        segments: &[usize],
        rates: &mut Rates,
    ) -> Result<(), UfgError> {
        let Some(thermal) = self.thermal else {
            return Ok(());
        };
        let (gas, offset, sigma, table) = thermal.gas(l);
        let (Some(floor), Some(table)) = (gas.floor(), table) else {
            return Ok(());
        };
        let (u, flux) = (self.u, &self.flux[i]);
        let ng = self.bounds.len() - 1;
        let top = self.groups.start;
        let x_floor = (self.high / floor).ln();
        let start = u
            .partition_point(|&v| v <= x_floor)
            .saturating_sub(1)
            .max(offset);
        let (lowest, highest) = (table.transfer[0], table.transfer[table.transfer.len() - 1]);
        for p in start..u.len() - 1 {
            let (ua, uz) = (u[p], u[p + 1]);
            let a = ua.max(x_floor);
            if uz <= a {
                continue;
            }
            let h = segments[p];
            let (mid, half) = (0.5 * (a + uz), 0.5 * (uz - a));
            for &(t, w) in &GAUSS_LEGENDRE_2 {
                let s = mid + half * t;
                let r = (s - ua) / (uz - ua);
                let lerp = |v0: f64, v1: f64| v0 + (v1 - v0) * r;
                let weight = w
                    * half
                    * lerp(sigma[p - offset], sigma[p + 1 - offset])
                    * lerp(flux[p], flux[p + 1]);
                let e_in = self.high * (-s).exp();
                // 表の最小の移行より下への割合は、表の累積がそこで 0 のときだけ 0 と分かる（入射ごとに 1 度だけ確かめる）
                let mut checked = false;
                for b in 1..ng {
                    let shift = self.bounds[b] - e_in;
                    if b > top + h {
                        // 下へ: 出射が境界以下の割合 C(E', E_b − E')
                        if shift <= lowest {
                            if !checked && cumulative(table, e_in, lowest, self.data[k].name)? > 0.0
                            {
                                return Err(UfgError::SabCoverage {
                                    name: self.data[k].name.to_string(),
                                    what: "smallest transfer",
                                    value: lowest,
                                    needed: shift,
                                });
                            }
                            checked = true;
                            continue;
                        }
                        let c = cumulative(table, e_in, shift.min(highest), self.data[k].name)?;
                        rates.below[h][b] += weight * c;
                    } else {
                        // 上へ: 表は最大の移行で 1 になる
                        if shift >= highest {
                            continue;
                        }
                        let c = cumulative(table, e_in, shift.max(lowest), self.data[k].name)?;
                        rates.above[h][b] += weight * (1.0 - c);
                    }
                }
            }
        }
        Ok(())
    }

    /// リング `i` の核種 `e`（組成の番号）の移行の行。要素は有限で、負なら弾性散乱断面積の [`NEGATIVE`] 倍までの丸め（それより
    /// 負なら、数え方の近似が減速計算と食い違っているので止める）
    fn rows(&self, i: usize, e: usize, rates: &Rates) -> Result<NuclideTransfer, UfgError> {
        let ng = self.bounds.len() - 1;
        let top = self.groups.start;
        let elastic = &self.nuclides[i][e].elastic;
        let rows = (0..self.groups.len())
            .map(|h| {
                let big = top + h;
                let phi = self.group_flux[i][h];
                let (below, above) = (&rates.below[h], &rates.above[h]);
                let mut row = vec![0.0; ng];
                for (g, value) in row.iter_mut().enumerate() {
                    *value = if g > big {
                        let next = if g + 1 < ng { below[g + 1] } else { 0.0 };
                        (below[g] - next) / phi
                    } else if g < big {
                        let previous = if g >= 1 { above[g] } else { 0.0 };
                        (above[g + 1] - previous) / phi
                    } else {
                        0.0
                    };
                }
                // 間引いた格子でも、ほかの群への移行を 2 つの格子の弾性散乱の比で拡大縮小しない。共鳴の点が落ちた核種で、群の
                // 端から離れた共鳴の分まで出る割合に掛かる（ADR-0208 の案 G2）
                let others: f64 = row.iter().sum();
                row[big] = elastic[h] - others;
                if let Some((to, &value)) = row
                    .iter()
                    .enumerate()
                    .find(|&(_, &v)| !(v.is_finite() && v >= -NEGATIVE * elastic[h]))
                {
                    return Err(UfgError::Transfer {
                        ring: i,
                        nuclide: e,
                        from: big,
                        to,
                        value,
                    });
                }
                Ok(row)
            })
            .collect::<Result<_, UfgError>>()?;
        Ok(NuclideTransfer {
            groups: self.groups.clone(),
            rows,
        })
    }
}

/// 表の累積 C(`e_in`, `shift`)。入射と移行は表の中にある（入射は範囲の下限から S(α, β) の上限まで、移行は呼ぶ側で表の範囲に
/// 収める。表が範囲を覆うことは核種のデータを引くときに確かめてある）
fn cumulative(table: &SabTransfer, e_in: f64, shift: f64, name: &str) -> Result<f64, UfgError> {
    table
        .cumulative_at(e_in, shift)
        .ok_or(UfgError::SabCoverage {
            name: name.to_string(),
            what: "incident energy",
            value: e_in,
            needed: e_in,
        })
}
