//! CMFD（粗メッシュ有限差分）の粗メッシュの固有値問題（ADR-0003、ADR-0028、ADR-0232）。
//!
//! 結合係数の作り方（形状と、高次の解き方の正味の中性子流からの補正）と、細かい量への戻し方は、使う側のステージが持つ。
//! ここは、組み立てたセルごとの対角のブロック、隣との結合、核分裂のスペクトルの項から、固有値と中性子束の比を求める。
//! `xlc-numeric` に置かないのは、それに依存する xs と lattice の実行ファイルのバイトを変えないため（ADR-0232）。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use xlc_numeric::lu::Lu;

/// 粗メッシュのセルの、同じ核分裂スペクトルの項の和。核分裂源は χ (νΣf · φ)
#[derive(Debug, Clone, PartialEq)]
pub struct Spectrum {
    pub chi: Vec<f64>,
    /// χ の和（核分裂源の生成率に掛ける）
    pub chi_sum: f64,
    /// [群] の体積 × νΣf（セルの中性子束あたり）
    pub nu_fission: Vec<f64>,
}

/// 粗メッシュの固有値問題を解く打ち切り条件。k と中性子束の変化から見積もった残りの誤差で判定する
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tolerances {
    /// 実効増倍率の残りの誤差の上限
    pub k: f64,
    /// 中性子束の相対的な残りの誤差の上限
    pub flux: f64,
    pub max_iterations: usize,
}

/// 粗メッシュの固有値問題のエラー。`cell` は組み立てた側のセルの番号（`Coarse::new` の `labels`）
#[derive(Debug, thiserror::Error)]
pub enum CoarseError {
    #[error("the coarse-mesh flux of cell {cell} in group {group} became {flux}")]
    NonPositiveSolution {
        cell: usize,
        group: usize,
        flux: f64,
    },
    #[error("the coarse-mesh matrix of cell {0} is singular")]
    Singular(usize),
    #[error("the coarse-mesh problem has no fission")]
    NoFission,
    #[error("the coarse-mesh eigenvalue problem did not converge in {0} iterations")]
    NotConverged(usize),
}

/// 粗メッシュの固有値問題の解
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    pub k: f64,
    /// [セル × 群] の中性子束の比（解 / 組み立てた中性子束、どちらも生成率 1 に規格化）
    pub ratio: Vec<f64>,
    pub iterations: usize,
}

/// セルごとに組み立てた粗メッシュの問題。方程式は A_c φ_c + Σ 係数 φ_隣 = χ (V νΣf · φ_c) / k。
///
/// 反復ごとに全セルを回るので、セルごとの量は平らな配列に並べる（セルごとの `Vec` をたどると、メモリの読み出しが反復の
/// 時間の大半になる）。
#[derive(Debug, Clone)]
pub struct Coarse {
    groups: usize,
    /// [セル] のエラーの表示に使う番号
    labels: Vec<usize>,
    /// [セル × 群] の組み立てた中性子束。解の比はこれに対して返し、反復もここから始める
    flux: Vec<f64>,
    /// [セル] の核分裂のスペクトルの項の範囲（`spectra` の添字）。セル c の項は `spectrum_offsets[c]..spectrum_offsets[c + 1]`
    spectrum_offsets: Vec<usize>,
    /// [項] の χ の和、[項 × 群] の χ と体積 × νΣf
    chi_sums: Vec<f64>,
    chi: Vec<f64>,
    nu_fission: Vec<f64>,
    /// [セル × 群 × 群] の対角のブロックの逆行列（行優先）
    inverses: Vec<f64>,
    /// [セル] の隣との結合の範囲（`others` の添字）
    neighbor_offsets: Vec<usize>,
    /// [結合] の隣のセルと、[結合 × 群] の係数
    others: Vec<usize>,
    coefficients: Vec<f64>,
}

impl Coarse {
    /// `matrices` は [セル] の対角のブロック（群×群、行優先）、`neighbors` は [セル] の隣との結合 (隣のセル, [群] の係数)、
    /// `fission` は [セル] の核分裂のスペクトルの項、`flux` は [セル × 群] の組み立てた中性子束、`labels` は [セル] のエラーに
    /// 出す番号
    pub fn new(
        groups: usize,
        labels: Vec<usize>,
        flux: Vec<f64>,
        fission: Vec<Vec<Spectrum>>,
        matrices: Vec<Vec<f64>>,
        neighbors: Vec<Vec<(usize, Vec<f64>)>>,
    ) -> Result<Self, CoarseError> {
        let ng = groups;
        let mut inverses = Vec::with_capacity(matrices.len() * ng * ng);
        let mut column = vec![0.0; ng];
        for (a, &label) in matrices.into_iter().zip(&labels) {
            let lu = Lu::new(a, ng).ok_or(CoarseError::Singular(label))?;
            let start = inverses.len();
            inverses.resize(start + ng * ng, 0.0);
            for j in 0..ng {
                column.fill(0.0);
                column[j] = 1.0;
                lu.solve(&mut column);
                for (i, v) in column.iter().enumerate() {
                    inverses[start + i * ng + j] = *v;
                }
            }
        }
        let mut spectrum_offsets = vec![0];
        let (mut chi_sums, mut chi, mut nu_fission) = (Vec::new(), Vec::new(), Vec::new());
        for terms in fission {
            for s in terms {
                chi_sums.push(s.chi_sum);
                chi.extend(s.chi);
                nu_fission.extend(s.nu_fission);
            }
            spectrum_offsets.push(chi_sums.len());
        }
        let mut neighbor_offsets = vec![0];
        let (mut others, mut coefficients) = (Vec::new(), Vec::new());
        for list in neighbors {
            for (other, c) in list {
                others.push(other);
                coefficients.extend(c);
            }
            neighbor_offsets.push(others.len());
        }
        Ok(Self {
            groups,
            labels,
            flux,
            spectrum_offsets,
            chi_sums,
            chi,
            nu_fission,
            inverses,
            neighbor_offsets,
            others,
            coefficients,
        })
    }

    fn cells(&self) -> usize {
        self.labels.len()
    }

    /// セル c の核分裂のスペクトルの項の (χ の和, [群] の χ, [群] の体積 × νΣf)
    fn spectra(&self, c: usize) -> impl Iterator<Item = (f64, &[f64], &[f64])> {
        let ng = self.groups;
        (self.spectrum_offsets[c]..self.spectrum_offsets[c + 1]).map(move |t| {
            (
                self.chi_sums[t],
                &self.chi[t * ng..][..ng],
                &self.nu_fission[t * ng..][..ng],
            )
        })
    }

    /// [セル × 群] の中性子束 `phi` の、すべてのセルの生成率の和
    fn production(&self, phi: &[f64]) -> f64 {
        let ng = self.groups;
        (0..self.cells())
            .flat_map(|c| {
                let p = &phi[c * ng..][..ng];
                self.spectra(c)
                    .map(move |(sum, _, f)| sum * f.iter().zip(p).map(|(a, b)| a * b).sum::<f64>())
            })
            .sum()
    }

    /// 固有値問題を、`k` と組み立てた中性子束から解き始める。
    ///
    /// べき乗法の各反復で、ブロック Gauss-Seidel をセルの順に 1 回だけ行う。収束すると Gauss-Seidel でも解が変わらない
    /// ので、内側の反復を収束させなくても固有値問題の解になる。
    pub fn solve(&self, k: f64, tolerances: &Tolerances) -> Result<Solution, CoarseError> {
        let ng = self.groups;
        let start = self.production(&self.flux);
        if !(start > 0.0 && start.is_finite()) {
            return Err(CoarseError::NoFission);
        }
        let mut phi: Vec<f64> = self.flux.iter().map(|f| f / start).collect();
        let mut previous = phi.clone();
        let mut k = k;
        let mut k_history = Remaining::default();
        let mut flux_history = Remaining::default();
        let mut rhs = vec![0.0; ng];
        for iteration in 1..=tolerances.max_iterations {
            previous.copy_from_slice(&phi);
            let inverse_k = 1.0 / k;
            for c in 0..self.cells() {
                rhs.fill(0.0);
                for (_, chi, f) in self.spectra(c) {
                    let rate = f
                        .iter()
                        .zip(&previous[c * ng..][..ng])
                        .map(|(a, b)| a * b)
                        .sum::<f64>()
                        * inverse_k;
                    for (r, x) in rhs.iter_mut().zip(chi) {
                        *r += x * rate;
                    }
                }
                for link in self.neighbor_offsets[c]..self.neighbor_offsets[c + 1] {
                    let other = &phi[self.others[link] * ng..][..ng];
                    for ((r, a), p) in rhs
                        .iter_mut()
                        .zip(&self.coefficients[link * ng..][..ng])
                        .zip(other)
                    {
                        *r -= a * p;
                    }
                }
                let inverse = &self.inverses[c * ng * ng..][..ng * ng];
                for (i, out) in phi[c * ng..][..ng].iter_mut().enumerate() {
                    *out = inverse[i * ng..][..ng]
                        .iter()
                        .zip(&rhs)
                        .map(|(a, b)| a * b)
                        .sum();
                }
            }
            // 割る前に正であることを確かめる（負の生成率で割ると、負の解も正に見える）
            if let Some(j) = phi.iter().position(|&f| !(f > 0.0 && f.is_finite())) {
                return Err(CoarseError::NonPositiveSolution {
                    cell: self.labels[j / ng],
                    group: j % ng,
                    flux: phi[j],
                });
            }
            let rate = self.production(&phi);
            let new_k = k * rate;
            for f in &mut phi {
                *f /= rate;
            }
            let change = phi
                .iter()
                .zip(&previous)
                .map(|(a, b)| (a - b).abs() / a.max(*b))
                .fold(0.0, f64::max);
            let k_remaining = k_history.update((new_k - k).abs());
            let flux_remaining = flux_history.update(change);
            k = new_k;
            if k_remaining < tolerances.k && flux_remaining < tolerances.flux {
                let ratio = phi
                    .iter()
                    .zip(&self.flux)
                    .map(|(new, old)| new / (old / start))
                    .collect();
                return Ok(Solution {
                    k,
                    ratio,
                    iterations: iteration,
                });
            }
        }
        Err(CoarseError::NotConverged(tolerances.max_iterations))
    }
}

/// 反復の変化量の列から、収束までに残る誤差を見積もる。
///
/// 変化が縮小率 ρ で等比的に減るとき、最新の変化を含めた残りの誤差は 変化/(1-ρ)。丸め誤差の水準（相対 1e-14）以下の変化は、
/// 増減が揺れても収束とみなす（許容が丸め誤差の水準より小さいと、揺れる変化で縮小率が 1 を超え、止まらない）。
#[derive(Debug, Default)]
struct Remaining {
    previous: Option<f64>,
}

impl Remaining {
    fn update(&mut self, change: f64) -> f64 {
        const ROUNDING: f64 = 1.0e-14;
        let estimate = match self.previous {
            _ if change <= ROUNDING => 0.0,
            Some(previous) if previous > 0.0 && change < previous => {
                let ratio = change / previous;
                change * (1.0 + ratio / (1.0 - ratio))
            }
            _ => f64::INFINITY,
        };
        self.previous = Some(change);
        estimate
    }
}

#[cfg(test)]
mod tests;
