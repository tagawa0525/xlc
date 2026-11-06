//! ノードを粗メッシュとし、面の部分中性子流で補正する CMFD（partial current CMFD。ADR-0232）。
//!
//! NEM のスイープの後のノードの平均中性子束 φ と面の部分中性子流から、面ごとに補正を決めた有限差分の固有値問題を作って
//! 解き、その解でノードの平均中性子束、展開係数、出ていく部分中性子流を更新する。
//!
//! 面の部分中性子流は、不均質な中性子束 ψ = f φ_s（面で連続。f はその側の ADF、φ_s は均質な面の中性子束）と正味の中性子流
//! J（面で連続）から、ĵ± = ψ/4 ± J/2 とする（面の両側で同じ値になる。NEM の均質な部分中性子流は ADF が 1 でない面で両側の
//! 値が違い、ADF の大きい側の入ってくる値は負になりうる）。方向 d の面でノード n とその隣 m が接するとき、n から出ていく
//! ĵ⁺ と n に入ってくる ĵ⁻ を
//!   ĵ⁺ = D̃/2 (f_n φ_n − f_m φ_m) + D̂⁺ φ_n
//!   ĵ⁻ = −D̃/2 (f_n φ_n − f_m φ_m) + D̂⁻ φ_m
//! と書く（正味の ĵ⁺ − ĵ⁻ は D̃ (f_n φ_n − f_m φ_m) + D̂⁺ φ_n − D̂⁻ φ_m）。D̃ は、ノードの中を平坦な中性子束とした有限差分で、
//! a = 2 D_n / h_n、b = 2 D_m / h_m として D̃ = a b / (a f_m + b f_n)（面で f φ_s が連続で、J = a (φ_n − φ_s,n) =
//! b (φ_s,m − φ_m) から）。D̂⁺ と D̂⁻ は、スイープの後の部分中性子流と φ に合うように決める。行列の n の対角に入る
//! A (D̃ f_n + D̂⁺) = A (ĵ⁺ + D̃ S / 2) / φ_n と、m との結合 −A (D̃ f_m + D̂⁻) = −A (ĵ⁻ + D̃ S / 2) / φ_m（S = f_n φ_n + f_m φ_m、
//! A は面の面積）は、部分中性子流と φ が正なら、それぞれ正と負になる（正味の中性子流で 1 つの D̂ を決める形は、結合が正に
//! なり、解が負になりうる。ADR-0232 の検討した案）。外周の面は、外向きの正味の中性子流をすべて D̂ φ_n で表す（D̂ = J / φ_n）。
//!
//! 解の φ' = r φ（r はノードと群ごとの比）で、面の部分中性子流を同じ式で求め直す:
//!   ĵ⁺' = r_n ĵ⁺ + D̃/2 f_m φ_m (r_n − r_m)、ĵ⁻' = r_m ĵ⁻ + D̃/2 f_n φ_n (r_m − r_n)
//! そこから n の均質な出ていく部分中性子流を o' = φ_s'/4 + J'/2 = (ĵ⁺' + ĵ⁻') / (2 f_n) + (ĵ⁺' − ĵ⁻') / 2 とする（隣の o' も
//! 同じ面の値から求めるので、NEM の面の条件で求める入ってくる部分中性子流は φ_s'/4 − J'/2 になる）。外周の面は r_n o。ノードの
//! 平均中性子束と展開係数の a3、a4 は r_n を掛け、a1、a2 は更新した面の中性子束とノードの平均から求め直す。

use rayon::prelude::*;
use xlc_cmfd::{Coarse, Spectrum, Tolerances};

use super::{CmfdError, CoarseRatios, Core, Face, MIN_NODES_PER_TASK, Options, Order, State};

/// 粗メッシュの固有値問題の反復の上限
const MAX_ITERATIONS: usize = 100_000;

/// CMFD を解く精度。NEM の許容と、スイープによる k と中性子束の変化（`k_change`、`flux_change`）の大きい方に、`ratios` の比を
/// 掛ける（NEM がまだ収束から遠いうちは、CMFD を細かく解いても次のスイープで崩れる。lattice と同じ考え方。ADR-0110）。
/// 比は 1 より十分小さくする。粗メッシュの反復は組み立てた中性子束から始まり、最初の数回の変化の縮小率は遅いモードの縮小率より
/// 小さいので、許容がスイープによる変化と同じ桁だと、遅いモードを直さないまま数回で止まり、外側の反復の変化が急に小さくなって
/// 許容の外で止まる（#0244）。CMFD の精度は収束した解を変えない（収束は、CMFD の更新を含めた外側の反復の変化で判定する）
pub(super) fn tolerances(
    nem: &Options,
    ratios: CoarseRatios,
    k_change: f64,
    flux_change: f64,
) -> Tolerances {
    Tolerances {
        k: ratios.k * nem.k_tolerance.max(k_change),
        flux: ratios.flux * nem.flux_tolerance.max(flux_change),
        max_iterations: MAX_ITERATIONS,
    }
}

/// ノードと隣のノードの面の、更新に使う量（`Core::coarse` が作り、`Core::accelerate` が使う）
#[derive(Debug, Clone, Copy, Default)]
struct Link {
    /// 隣のノード
    other: usize,
    /// D̃/2 f_m φ_m と D̃/2 f_n φ_n（ĵ⁺' と ĵ⁻' の比の差の係数）
    toward: f64,
    back: f64,
    /// ĵ⁺ と ĵ⁻
    outgoing: f64,
    incoming: f64,
    /// 両側の平均で割った自分の側の ADF
    adf: f64,
}

/// ノード 1 つの粗メッシュの式（`Core::assemble`）
struct Assembled {
    /// 対角のブロック（群×群、行優先）
    matrix: Vec<f64>,
    /// 隣との結合 (隣のノード, [群] の係数)
    neighbors: Vec<(usize, Vec<f64>)>,
    spectrum: Spectrum,
    /// [方向 × 面 × 群] の面の量
    links: Vec<Link>,
}

impl Core<'_> {
    /// スイープの後の `state` から粗メッシュの問題を作って `k` から解き、`state` を更新する。粗メッシュの問題の k を返す。粗メッシュの
    /// セルは `order` の順に解く（ADR-0234）。ノードごとの組み立てと更新は互いに独立なので並列に行う（どのスレッドでも同じ値）
    pub(super) fn accelerate(
        &self,
        state: &mut State,
        k: f64,
        tolerances: &Tolerances,
        order: Order,
    ) -> Result<f64, CmfdError> {
        let ng = self.groups;
        let (coarse, links) =
            tracing::trace_span!("homogenize").in_scope(|| self.coarse(state, order))?;
        let solution = tracing::trace_span!("coarse").in_scope(|| coarse.solve(k, tolerances))?;
        let _scale = tracing::trace_span!("scale").entered();
        tracing::trace!(
            iterations = solution.iterations,
            k = solution.k,
            "coarse-mesh eigenvalue"
        );
        let ratio = |n: usize, g: usize| solution.ratio[n * ng + g];
        let State {
            flux,
            coefficients: [c0, c1, c2],
            out_minus: [m0, m1, m2],
            out_plus: [p0, p1, p2],
        } = state;
        (
            flux.par_iter_mut(),
            c0.par_iter_mut(),
            c1.par_iter_mut(),
            c2.par_iter_mut(),
            m0.par_iter_mut(),
            m1.par_iter_mut(),
            m2.par_iter_mut(),
            p0.par_iter_mut(),
            p1.par_iter_mut(),
            p2.par_iter_mut(),
        )
            .into_par_iter()
            .enumerate()
            .with_min_len(MIN_NODES_PER_TASK)
            .for_each(|(n, (f, c0, c1, c2, m0, m1, m2, p0, p1, p2))| {
                for (d, (c, minus, plus_side)) in [(c0, m0, p0), (c1, m1, p1), (c2, m2, p2)]
                    .into_iter()
                    .enumerate()
                {
                    for (side, (plus, out)) in
                        [(false, minus), (true, plus_side)].into_iter().enumerate()
                    {
                        for (g, o) in out.iter_mut().enumerate() {
                            let r = ratio(n, g);
                            *o = match self.face(n, d, plus) {
                                Face::Node(_) => {
                                    let link = &links[((n * 3 + d) * 2 + side) * ng + g];
                                    let r_m = ratio(link.other, g);
                                    let outgoing = r * link.outgoing + link.toward * (r - r_m);
                                    let incoming = r_m * link.incoming + link.back * (r_m - r);
                                    0.5 * (outgoing + incoming) / link.adf
                                        + 0.5 * (outgoing - incoming)
                                }
                                Face::Edge(_) => r * *o,
                            };
                        }
                    }
                    for (g, coefficients) in c.iter_mut().enumerate() {
                        let r = ratio(n, g);
                        for x in coefficients {
                            *x *= r;
                        }
                    }
                }
                for (g, phi) in f.iter_mut().enumerate() {
                    *phi *= ratio(n, g);
                }
            });
        // a1 = φ_+ − φ_−、a2 = φ_+ + φ_− − 2 φ̄ は、更新した部分中性子流の面の中性子束 2 (o + i) とノードの平均から（比を掛けた
        // 値のままだと、面の量と食い違う）。入ってくる部分中性子流は隣の更新した値から求めるので、出ていく値をすべて更新してから。
        // 入ってくる部分中性子流は展開係数を読まないので、展開係数を外に出して並列に書く
        let mut coefficients = std::mem::take(&mut state.coefficients);
        let read = &*state;
        for (d, by_node) in coefficients.iter_mut().enumerate() {
            by_node
                .par_iter_mut()
                .enumerate()
                .with_min_len(MIN_NODES_PER_TASK)
                .for_each(|(n, row)| {
                    for (g, c) in row.iter_mut().enumerate() {
                        let (in_minus, in_plus) = self.incoming(read, n, d, g);
                        let plus = 2.0 * (read.out_plus[d][n][g] + in_plus);
                        let minus = 2.0 * (read.out_minus[d][n][g] + in_minus);
                        c[0] = plus - minus;
                        c[1] = plus + minus - 2.0 * read.flux[n][g];
                    }
                });
        }
        state.coefficients = coefficients;
        Ok(solution.k)
    }

    /// `state` のノードの平均中性子束と面の部分中性子流から、補正した有限差分の粗メッシュの問題と、[ノード × 方向 × 面 × 群] の
    /// 面の量（外周の面は使わない）を作る。`order` が `RedBlack` なら、セルを市松の色ごとに解く形にする
    fn coarse(&self, state: &State, order: Order) -> Result<(Coarse, Vec<Link>), CmfdError> {
        let ng = self.groups;
        let nodes = self.nodes();
        if let Some((n, g)) = (0..nodes)
            .flat_map(|n| (0..ng).map(move |g| (n, g)))
            .find(|&(n, g)| !(state.flux[n][g] > 0.0 && state.flux[n][g].is_finite()))
        {
            return Err(CmfdError::NonPositiveFlux {
                node: self.site(n),
                group: g,
                flux: state.flux[n][g],
            });
        }
        let assembled: Vec<Assembled> = (0..nodes)
            .into_par_iter()
            .with_min_len(MIN_NODES_PER_TASK)
            .map(|n| self.assemble(state, n))
            .collect();
        let mut links = Vec::with_capacity(nodes * 6 * ng);
        let mut matrices = Vec::with_capacity(nodes);
        let mut neighbors = Vec::with_capacity(nodes);
        let mut fission = Vec::with_capacity(nodes);
        let mut flux = Vec::with_capacity(nodes * ng);
        for (n, a) in assembled.into_iter().enumerate() {
            links.extend(a.links);
            matrices.push(a.matrix);
            neighbors.push(a.neighbors);
            fission.push(vec![a.spectrum]);
            flux.extend_from_slice(&state.flux[n]);
        }
        let labels = (0..nodes).map(|n| self.site(n)).collect();
        let coarse = match order {
            Order::Lexicographic => Coarse::new(ng, labels, flux, fission, matrices, neighbors)?,
            Order::RedBlack => {
                Coarse::colored(ng, labels, flux, fission, matrices, neighbors, &self.color)?
            }
        };
        Ok((coarse, links))
    }

    /// ノード n の粗メッシュの式と面の量。
    ///
    /// D̃ は有限差分の値か、部分中性子流の負の値を補って結合の符号を保つ値（ĵ + D̃ S / 2 ≥ −ĵ > 0 になる −4 ĵ / S）の大きい方。
    /// D̃ は収束した解を変えない（D̂ が部分中性子流に合わせる）。面の条件は両側の ADF の比だけで決まるので、ADF は両側の平均で
    /// 割った値を使う（両側に同じ数を掛けても反復が変わらない）。中性子束が正であることは `coarse` で確かめてある
    fn assemble(&self, state: &State, n: usize) -> Assembled {
        let ng = self.groups;
        let m = self.material(n);
        let volume = self.volume(n);
        let phi = &state.flux[n];
        let mut links = vec![Link::default(); 6 * ng];
        // 対角のブロック: 行が群 g の式、列が中性子束の群。除去と、ほかの群からの散乱
        let mut a = vec![0.0; ng * ng];
        for g in 0..ng {
            a[g * ng + g] = volume * (m.total[g] - m.scatter[g][g]);
            for from in (0..ng).filter(|&from| from != g) {
                a[g * ng + from] = -volume * m.scatter[from][g];
            }
        }
        let mut list = Vec::new();
        for d in 0..3 {
            let h = self.width(n, d);
            let area = volume / h;
            let mut couplings = [vec![0.0; ng], vec![0.0; ng]];
            for g in 0..ng {
                let (in_minus, in_plus) = self.incoming(state, n, d, g);
                for (side, plus) in [false, true].into_iter().enumerate() {
                    let (o, i) = if plus {
                        (state.out_plus[d][n][g], in_plus)
                    } else {
                        (state.out_minus[d][n][g], in_minus)
                    };
                    let Face::Node(other) = self.face(n, d, plus) else {
                        // 外周の面: 外向きの正味の中性子流を D̂ φ_n で表す
                        a[g * ng + g] += area * (o - i) / phi[g];
                        continue;
                    };
                    let (f_n, f_m) = (self.adf(n, d, plus, g), self.adf(other, d, !plus, g));
                    let mean = 0.5 * (f_n + f_m);
                    let (f_n, f_m) = (f_n / mean, f_m / mean);
                    // 不均質な部分中性子流 ψ/4 ± J/2（ψ = 2 f (o + i)、J = o − i）
                    let (quarter, half) = (0.5 * f_n * (o + i), 0.5 * (o - i));
                    let (outgoing, incoming) = (quarter + half, quarter - half);
                    let phi_m = state.flux[other][g];
                    let sum = f_n * phi[g] + f_m * phi_m;
                    let b_n = 2.0 * m.diffusion[g] / h;
                    let b_m = 2.0 * self.material(other).diffusion[g] / self.width(other, d);
                    let tilde = (b_n * b_m / (b_n * f_m + b_m * f_n))
                        .max(-4.0 * outgoing.min(incoming) / sum);
                    let common = 0.5 * tilde * sum;
                    a[g * ng + g] += area * (outgoing + common) / phi[g];
                    couplings[side][g] = -area * (incoming + common) / phi_m;
                    links[(d * 2 + side) * ng + g] = Link {
                        other,
                        toward: 0.5 * tilde * f_m * phi_m,
                        back: 0.5 * tilde * f_n * phi[g],
                        outgoing,
                        incoming,
                        adf: f_n,
                    };
                }
            }
            for (side, plus) in [false, true].into_iter().enumerate() {
                if let Face::Node(other) = self.face(n, d, plus) {
                    list.push((other, std::mem::take(&mut couplings[side])));
                }
            }
        }
        Assembled {
            matrix: a,
            neighbors: list,
            spectrum: Spectrum {
                chi: m.chi.clone(),
                chi_sum: m.chi.iter().sum(),
                nu_fission: m.nu_fission.iter().map(|f| volume * f).collect(),
            },
            links,
        }
    }
}
