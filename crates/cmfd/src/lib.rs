//! CMFD（粗メッシュ有限差分）の粗メッシュの固有値問題（ADR-0003、ADR-0028、ADR-0232）。
//!
//! 結合係数の作り方（形状と、高次の解き方の正味の中性子流からの補正）と、細かい量への戻し方は、使う側のステージが持つ。
//! ここは、組み立てたセルごとの対角のブロック、隣との結合、核分裂のスペクトルの項から、固有値と中性子束の比を求める。
//! `xlc-numeric` に置かないのは、それに依存する xs の実行ファイルのバイトを、xs のコードを変えない変更で変えないため（ADR-0232）。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use std::ops::Range;

use rayon::prelude::*;

pub use blocks::BlockSolve;
use blocks::Blocks;

mod blocks;

/// 色ごとに並列に解くときの、1 つのスレッドがまとめて解くセルの数の下限（細かく分けると、仕事を配る時間が解く時間を超える）
const MIN_CELLS_PER_TASK: usize = 256;

/// 色ごとに解くときに、1 つのスレッドがまとめて解き、生成率の部分和を取るセルの数。和の順をスレッドの数によらず決めるため、
/// 固定の大きさで分ける
const CELLS_PER_SUM: usize = 256;

/// 粗メッシュのセルの、同じ核分裂スペクトルの項の和。核分裂源は χ (νΣf · φ)
#[derive(Debug, Clone, PartialEq)]
pub struct Spectrum {
    pub chi: Vec<f64>,
    /// χ の和（核分裂源の生成率に掛ける）
    pub chi_sum: f64,
    /// [群] の体積 × νΣf（セルの中性子束あたり）
    pub nu_fission: Vec<f64>,
}

impl Spectrum {
    /// セルの中性子束 `phi`（[群]）での体積 × νΣf · φ
    pub fn rate(&self, phi: &[f64]) -> f64 {
        self.nu_fission.iter().zip(phi).map(|(a, b)| a * b).sum()
    }
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
    /// 色ごとに解くとき、隣と同じ色のセル（`cell` と `other` は組み立てた側の番号）
    #[error("cell {cell} couples to cell {other} of its own color")]
    SameColor { cell: usize, other: usize },
    #[error("{given} colors are given for {cells} cells")]
    Colors { given: usize, cells: usize },
    /// 生成率の重み（`Coarse::weighted`）の数がセルの数と違う
    #[error("{given} production weights are given for {cells} cells")]
    Weights { given: usize, cells: usize },
    /// 生成率の重みが正の有限でない（`cell` は最初のそういう重みの、与えた順の番号）
    #[error("production weight {cell} is {value}, not a positive number")]
    Weight { cell: usize, value: f64 },
    /// 対角のブロックの成分の数が群の数の 2 乗でないか、群の数が 0（`cell` は組み立てた側の番号）
    #[error("the coarse-mesh block of cell {cell} has {len} entries for {groups} groups")]
    BlockSize {
        cell: usize,
        len: usize,
        groups: usize,
    },
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
    /// [セル] の対角のブロックの LU 分解（逆行列を掛ける形は後退安定でなく、条件数の悪いブロックで誤差が大きい。ADR-0275）
    blocks: Blocks,
    /// [セル] の隣との結合の範囲（`others` の添字）
    neighbor_offsets: Vec<usize>,
    /// [結合] の隣のセルと、[結合 × 群] の係数
    others: Vec<usize>,
    coefficients: Vec<f64>,
    /// 色ごとに解くときの色（`Coarse::colored`）。`None` はセルの順に 1 つずつ解く
    colors: Option<Colors>,
    /// [セル] の生成率の重み（`Coarse::weighted`。既定は 1）
    weights: Vec<f64>,
}

/// 色ごとに解くセルの並び。セルは色の順に並べ直して持ち、同じ色のセルは続いた範囲にある
#[derive(Debug, Clone)]
struct Colors {
    /// [色] の、並べ直したセルの範囲
    ranges: Vec<Range<usize>>,
    /// [与えたセル] の、並べ直した番号
    position: Vec<usize>,
}

impl Coarse {
    /// `matrices` は [セル] の対角のブロック（群×群、行優先）、`neighbors` は [セル] の隣との結合 (隣のセル, [群] の係数)、
    /// `fission` は [セル] の核分裂のスペクトルの項、`flux` は [セル × 群] の組み立てた中性子束、`labels` は [セル] のエラーに
    /// 出す番号、`solve` は対角のブロックの解き方
    pub fn new(
        solve: BlockSolve,
        groups: usize,
        labels: Vec<usize>,
        flux: Vec<f64>,
        fission: Vec<Vec<Spectrum>>,
        matrices: Vec<Vec<f64>>,
        neighbors: Vec<Vec<(usize, Vec<f64>)>>,
    ) -> Result<Self, CoarseError> {
        let ng = groups;
        let cells = labels.len();
        if let Some((c, a)) = matrices
            .iter()
            .enumerate()
            .find(|(_, a)| ng == 0 || a.len() != ng * ng)
        {
            return Err(CoarseError::BlockSize {
                cell: labels[c],
                len: a.len(),
                groups: ng,
            });
        }
        let blocks =
            Blocks::new(matrices, ng, solve).map_err(|c| CoarseError::Singular(labels[c]))?;
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
            blocks,
            neighbor_offsets,
            others,
            coefficients,
            colors: None,
            weights: vec![1.0; cells],
        })
    }

    /// `new` と同じ問題を、`colors`（[セル] の色の番号）の順に、同じ色のセルを並列に解く形で作る（ADR-0234）。同じ色のセルは
    /// 互いの値を読まないので、色の中の解く順とスレッドの数は解を変えない。隣と同じ色のセルがあればエラー
    pub fn colored(
        solve: BlockSolve,
        groups: usize,
        labels: Vec<usize>,
        flux: Vec<f64>,
        fission: Vec<Vec<Spectrum>>,
        matrices: Vec<Vec<f64>>,
        neighbors: Vec<Vec<(usize, Vec<f64>)>>,
        colors: &[usize],
    ) -> Result<Self, CoarseError> {
        let cells = labels.len();
        if colors.len() != cells {
            return Err(CoarseError::Colors {
                given: colors.len(),
                cells,
            });
        }
        for (c, list) in neighbors.iter().enumerate() {
            if let Some(&(other, _)) = list.iter().find(|(other, _)| colors[*other] == colors[c]) {
                return Err(CoarseError::SameColor {
                    cell: labels[c],
                    other: labels[other],
                });
            }
        }
        // 色の順（同じ色はセルの順）に並べ直す
        let mut order: Vec<usize> = (0..cells).collect();
        order.sort_by_key(|&c| colors[c]);
        let mut position = vec![0; cells];
        for (new, &old) in order.iter().enumerate() {
            position[old] = new;
        }
        let mut ranges: Vec<Range<usize>> = Vec::new();
        for (new, &old) in order.iter().enumerate() {
            match ranges.last_mut() {
                Some(last) if colors[order[last.start]] == colors[old] => last.end = new + 1,
                _ => ranges.push(new..new + 1),
            }
        }
        let ng = groups;
        let (mut fission, mut matrices, mut neighbors) = (fission, matrices, neighbors);
        let mut coarse = Self::new(
            solve,
            groups,
            order.iter().map(|&c| labels[c]).collect(),
            order
                .iter()
                .flat_map(|&c| flux[c * ng..][..ng].iter().copied())
                .collect(),
            order
                .iter()
                .map(|&c| std::mem::take(&mut fission[c]))
                .collect(),
            order
                .iter()
                .map(|&c| std::mem::take(&mut matrices[c]))
                .collect(),
            order
                .iter()
                .map(|&c| {
                    std::mem::take(&mut neighbors[c])
                        .into_iter()
                        .map(|(other, coefficients)| (position[other], coefficients))
                        .collect()
                })
                .collect(),
        )?;
        coarse.colors = Some(Colors { ranges, position });
        Ok(coarse)
    }

    /// 生成率を数えるときに、セルごとに `weights`（[セル]、与えたセルの順。正の有限）を掛ける問題にする。対称性の代表のセルを、
    /// その組のセルの数だけ数えるため。重みは解の比の大きさ（生成率 1 の規格化）を決め、固有値を変えない
    pub fn weighted(mut self, weights: Vec<f64>) -> Result<Self, CoarseError> {
        let cells = self.cells();
        if weights.len() != cells {
            return Err(CoarseError::Weights {
                given: weights.len(),
                cells,
            });
        }
        if let Some(cell) = weights.iter().position(|w| !(*w > 0.0 && w.is_finite())) {
            return Err(CoarseError::Weight {
                cell,
                value: weights[cell],
            });
        }
        self.weights = match &self.colors {
            None => weights,
            Some(colors) => {
                let mut ordered = vec![0.0; cells];
                for (c, &p) in colors.position.iter().enumerate() {
                    ordered[p] = weights[c];
                }
                ordered
            }
        };
        Ok(self)
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

    /// [セル × 群] の中性子束 `phi` の、すべてのセルの生成率の和（セルごとに重みを掛ける）
    fn production(&self, phi: &[f64]) -> f64 {
        let ng = self.groups;
        (0..self.cells())
            .flat_map(|c| {
                let p = &phi[c * ng..][..ng];
                self.cell_production(c, p)
            })
            .sum()
    }

    /// セル c の中性子束 `p`（[群]）の、核分裂のスペクトルの項ごとの重みを掛けた生成率
    fn cell_production<'a>(&'a self, c: usize, p: &'a [f64]) -> impl Iterator<Item = f64> + 'a {
        let weight = self.weights[c];
        self.spectra(c).map(move |(sum, _, f)| {
            weight * (sum * f.iter().zip(p).map(|(a, b)| a * b).sum::<f64>())
        })
    }

    /// セル c の式の右辺（`fixed`（[群]。無ければ 0）と核分裂源 χ (V νΣf · φ_前) / k の和から、隣との結合 × 隣の中性子束を
    /// 引いたもの）を `rhs` に書く。`neighbor` は隣のセルの [群] の中性子束
    fn source<'a>(
        &self,
        c: usize,
        rhs: &mut [f64],
        fixed: Option<&[f64]>,
        previous: &[f64],
        neighbor: impl Fn(usize) -> &'a [f64],
        inverse_k: f64,
    ) {
        let ng = self.groups;
        match fixed {
            Some(values) => rhs.copy_from_slice(values),
            None => rhs.fill(0.0),
        }
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
            let other = neighbor(self.others[link]);
            for ((r, a), p) in rhs
                .iter_mut()
                .zip(&self.coefficients[link * ng..][..ng])
                .zip(other)
            {
                *r -= a * p;
            }
        }
    }

    /// セル c の中性子束 `out` を、対角のブロックの式 A_c out = `rhs` の解にする
    fn invert(&self, c: usize, rhs: &[f64], out: &mut [f64]) {
        self.blocks.solve(c, rhs, out);
    }

    /// [セル × 群] の値 `j` が正でないことのエラー
    fn non_positive(&self, phi: &[f64], j: usize) -> CoarseError {
        let ng = self.groups;
        CoarseError::NonPositiveSolution {
            cell: self.labels[j / ng],
            group: j % ng,
            flux: phi[j],
        }
    }

    /// 色 `range` のセルを並列に解き、解いたセルの生成率の部分和を、決まった大きさ（`CELLS_PER_SUM`）の塊の順に返す（和の順を
    /// スレッドの数によらず決める）。隣は別の色なので、色の中で読む値は解く順によらない。塊に正でない中性子束があれば、その
    /// 塊の部分和の代わりに、最初の値の [セル × 群] の番号を返す（負の生成率で割ると、負の解も正に見える）
    fn sweep_color(
        &self,
        range: &Range<usize>,
        phi: &mut [f64],
        previous: &[f64],
        inverse_k: f64,
    ) -> Vec<Result<f64, usize>> {
        let ng = self.groups;
        let (before, rest) = phi.split_at_mut(range.start * ng);
        let (current, after) = rest.split_at_mut(range.len() * ng);
        let (before, after) = (&*before, &*after);
        // 隣と同じ色のセルが無いことは `colored` で確かめてある
        let neighbor = |other: usize| -> &[f64] {
            if other < range.start {
                &before[other * ng..][..ng]
            } else {
                &after[(other - range.end) * ng..][..ng]
            }
        };
        current
            .par_chunks_mut(CELLS_PER_SUM * ng)
            .enumerate()
            .map_init(
                || vec![0.0; ng],
                |rhs, (chunk, values)| {
                    let first = range.start + chunk * CELLS_PER_SUM;
                    for (i, out) in values.chunks_mut(ng).enumerate() {
                        self.source(first + i, rhs, None, previous, neighbor, inverse_k);
                        self.invert(first + i, rhs, out);
                    }
                    if let Some(j) = values.iter().position(|&f| !(f > 0.0 && f.is_finite())) {
                        return Err(first * ng + j);
                    }
                    Ok(values
                        .chunks(ng)
                        .enumerate()
                        .flat_map(|(i, p)| self.cell_production(first + i, p))
                        .sum())
                },
            )
            .collect()
    }

    /// ブロック Gauss-Seidel を 1 回行う。`x`（[セル × 群]、与えたセルの順）を前の値とし、セルの順（`colored` は色の順）に、
    /// A_c x_c = s_c + χ (V νΣf · x前_c) `inverse_k` − Σ 係数 x_隣 を解いて x_c を置き換える。隣はその時点の値（先に解いたセルは
    /// 新しい値）。s は `source`（[セル × 群]、与えたセルの順）で、無ければ 0
    pub fn sweep(&self, x: &mut [f64], source: Option<&[f64]>, inverse_k: f64) {
        let ng = self.groups;
        match &self.colors {
            None => {
                let previous = x.to_vec();
                self.pass(x, &previous, source, inverse_k);
            }
            Some(colors) => {
                // 色の順に並べ直す。同じ色のセルは互いに読まないので、色の順に 1 つずつ解くと、色ごとに並列に解くのと同じ
                let ordered = |values: &[f64]| -> Vec<f64> {
                    let mut out = vec![0.0; values.len()];
                    for (c, &p) in colors.position.iter().enumerate() {
                        out[p * ng..][..ng].copy_from_slice(&values[c * ng..][..ng]);
                    }
                    out
                };
                let previous = ordered(x);
                let fixed = source.map(ordered);
                let mut phi = previous.clone();
                self.pass(&mut phi, &previous, fixed.as_deref(), inverse_k);
                for (c, &p) in colors.position.iter().enumerate() {
                    x[c * ng..][..ng].copy_from_slice(&phi[p * ng..][..ng]);
                }
            }
        }
    }

    /// 並べ直したセルの順に、`phi` をその場で 1 回解く（`sweep`。`previous` と `fixed` も並べ直したセルの順）
    fn pass(&self, phi: &mut [f64], previous: &[f64], fixed: Option<&[f64]>, inverse_k: f64) {
        let ng = self.groups;
        let mut rhs = vec![0.0; ng];
        for c in 0..self.cells() {
            self.source(
                c,
                &mut rhs,
                fixed.map(|f| &f[c * ng..][..ng]),
                previous,
                |other| &phi[other * ng..][..ng],
                inverse_k,
            );
            self.invert(c, &rhs, &mut phi[c * ng..][..ng]);
        }
    }

    /// 固有値問題を、`k` と組み立てた中性子束から解き始める。
    ///
    /// べき乗法の各反復で、ブロック Gauss-Seidel を 1 回だけ行う（`new` はセルの順、`colored` は色の順に、同じ色のセルを並列に）。
    /// 収束すると Gauss-Seidel でも解が変わらないので、内側の反復を収束させなくても固有値問題の解になる。
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
            let inverse_k = 1.0 / k;
            let (rate, change) = match &self.colors {
                None => {
                    previous.copy_from_slice(&phi);
                    for c in 0..self.cells() {
                        self.source(
                            c,
                            &mut rhs,
                            None,
                            &previous,
                            |other| &phi[other * ng..][..ng],
                            inverse_k,
                        );
                        self.invert(c, &rhs, &mut phi[c * ng..][..ng]);
                    }
                    // 割る前に正であることを確かめる（負の生成率で割ると、負の解も正に見える）
                    if let Some(j) = phi.iter().position(|&f| !(f > 0.0 && f.is_finite())) {
                        return Err(self.non_positive(&phi, j));
                    }
                    let rate = self.production(&phi);
                    for f in &mut phi {
                        *f /= rate;
                    }
                    let change = phi
                        .iter()
                        .zip(&previous)
                        .map(|(a, b)| relative_change(*a, *b))
                        .fold(0.0, f64::max);
                    (rate, change)
                }
                Some(colors) => {
                    // `previous` は前の反復の終わりに `phi` と同じにしてある
                    let mut rate = 0.0;
                    for range in &colors.ranges {
                        for part in self.sweep_color(range, &mut phi, &previous, inverse_k) {
                            rate += part.map_err(|j| self.non_positive(&phi, j))?;
                        }
                    }
                    // 最大はどの順に取っても同じ
                    let change = phi
                        .par_iter_mut()
                        .zip(previous.par_iter_mut())
                        .with_min_len(MIN_CELLS_PER_TASK * ng)
                        .map(|(a, b)| {
                            *a /= rate;
                            let change = relative_change(*a, *b);
                            *b = *a;
                            change
                        })
                        .reduce(|| 0.0, f64::max);
                    (rate, change)
                }
            };
            let new_k = k * rate;
            let k_remaining = k_history.update((new_k - k).abs());
            let flux_remaining = flux_history.update(change);
            k = new_k;
            if k_remaining < tolerances.k && flux_remaining < tolerances.flux {
                let ratio: Vec<f64> = phi
                    .iter()
                    .zip(&self.flux)
                    .map(|(new, old)| new / (old / start))
                    .collect();
                // 色の順に並べ直したセルは、与えたセルの順に戻す
                let ratio = match &self.colors {
                    None => ratio,
                    Some(colors) => colors
                        .position
                        .iter()
                        .flat_map(|&p| ratio[p * ng..][..ng].iter().copied())
                        .collect(),
                };
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

/// 中性子束の相対変化（どちらも正）
fn relative_change(new: f64, old: f64) -> f64 {
    (new - old).abs() / new.max(old)
}

/// 反復の変化量の列から、収束までに残る誤差を見積もる（粗メッシュの固有値問題と、lattice の MOC と core の NEM の外側の反復）。
///
/// 変化が縮小率 ρ で等比的に減るとき、最新の変化を含めた残りの誤差は 変化/(1-ρ)。1e-14 以下の変化は、増減が揺れても収束とみなす
/// （許容が丸め誤差の水準より小さいと、揺れる変化で縮小率が 1 を超え、止まらない）。閾値は絶対値なので、変化は 1 の桁で測ったもの
/// （1 の桁の k の変化か、中性子束の相対変化）を渡す。1e-14 はその丸め誤差の水準。
#[derive(Debug, Default)]
pub struct Remaining {
    previous: Option<f64>,
}

impl Remaining {
    /// 最新の変化量 `change`（k の変化か、中性子束の相対変化の最大）を受け取り、残りの誤差の見積もりを返す
    pub fn update(&mut self, change: f64) -> f64 {
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
