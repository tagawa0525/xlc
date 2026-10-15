//! 非分離共鳴のラダーのサンプリングによる自己遮蔽（ADR-0177）。
//!
//! 系列（l, J）ごとに、平均のパラメータから共鳴の並び（ラダー）を乱数で作り、その上の単準位 Breit-Wigner の断面積を
//! 等間隔の点を乱数でずらしたエネルギーで求め、σ0 ごとの Bondarenko の平均 ⟨σx/(σt + σ0)⟩/⟨1/(σt + σ0)⟩ を直接に求める。
//! NJOY PURR（purr.f90 の unrest と ladr2）の「bondarenko cross sections by direct sampling」と同じ量。PURR と違い、弾性散乱は
//! 同じ系列の共鳴の干渉を多準位で入れ（`add_ladder`）、幅は χ² の分位点の表でなく正規乱数の 2 乗の和で引き、遠い共鳴の寄与も
//! 切らずに足す（窓の外の寄与は、窓の幅で誤差を抑える）。

use std::f64::consts::PI;

use num_complex::Complex64;

use super::{ShieldedSigma, Sigma};
use crate::endf::WidthDof;
use crate::special::faddeeva;

/// サンプルの数。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sampling {
    /// ラダーの組の数
    pub ladders: usize,
    /// ラダーごとの、断面積を求めるエネルギーの数
    pub samples: usize,
}

impl Sampling {
    /// ライブラリを作るときの数。統計誤差の測定（ADR-0177）で決める
    pub const LIBRARY: Self = Self {
        ladders: 512,
        samples: 250,
    };
}

/// エネルギーを引く窓の幅（系列の最大の準位間隔の倍数）
const SAMPLE_SPACINGS: f64 = 50.0;

/// 窓の外に並べる共鳴の数（準位間隔の倍数）の下限と上限。上限では、裾を落とす誤差は Γ/(π D × 上限) になる
const MIN_MARGIN_SPACINGS: f64 = 50.0;
const MAX_MARGIN_SPACINGS: f64 = 1000.0;

/// 窓の外の共鳴の裾を落とす誤差の目安（平均断面積に対する相対）。1 つの系列の、窓から L 離れた外の共鳴の
/// ローレンツの裾の寄与は平均の Γ/(π L) 倍なので、L = Γ/(π D TAIL) 個の準位間隔まで並べる
const TAIL: f64 = 1.0e-4;

/// 弾性散乱 + 負の背景が負になる点の弾性散乱 (barn)。NJOY PURR と同じ値
const ELASTIC_FLOOR: f64 = 1.0e-6;

/// |x| か y がこれより大きいと、w(z) を漸近形 i/(√π z) で求める（相対誤差 1/(2|z|²) ≈ 5e-5）
const ASYMPTOTIC: f64 = 100.0;

/// 1 つの系列の、エネルギー E での平均のパラメータ。幅 (eV) は平均で、自由度の χ² 分布でゆらぐ。
#[derive(Debug, Clone, Copy)]
pub(super) struct LadderSequence {
    /// 4π/k² × gJ (barn)
    pub(super) scale: f64,
    pub(super) gn: f64,
    pub(super) gf: f64,
    pub(super) gg: f64,
    pub(super) gx: f64,
    pub(super) dof_n: WidthDof,
    pub(super) dof_f: WidthDof,
    pub(super) dof_x: WidthDof,
    /// 準位間隔 D (eV)
    pub(super) d: f64,
    /// ドップラー幅 Δ (eV)
    pub(super) delta: f64,
    /// cos 2φ と sin 2φ（φ は位相シフト）
    pub(super) cos2: f64,
    pub(super) sin2p: f64,
}

/// σ0 ごとの実効断面積と、同じサンプルの平均（σ0 → ∞）。
#[derive(Debug, Clone, PartialEq)]
pub struct Bondarenko {
    pub shielded: Vec<ShieldedSigma>,
    pub mean: Sigma,
    /// 弾性散乱 + 負の背景が負になり、1 μb にした点の数
    pub floored: usize,
}

/// 系列 `sequences` のラダーで、σ0 ごとの実効断面積を求める。`potential` はポテンシャル散乱 (barn)、`background` は
/// MF3 の背景断面積。サンプルの数は `sampling`、乱数は種 `seed` から決まる。弾性散乱 + 負の背景が負の点は弾性散乱を 1 μb にし（数を `floored` に数える）、それでも全断面積が 0 以下になる点があれば `Err` でそのエネルギーの位置
/// （窓の中心からの eV）を返す。
pub(super) fn bondarenko(
    sequences: &[LadderSequence],
    potential: f64,
    background: Sigma,
    sigma0: &[f64],
    sampling: Sampling,
    seed: u64,
) -> Result<Bondarenko, f64> {
    let mut rng = Rng::new(seed);
    let d_max = sequences.iter().map(|q| q.d).fold(0.0, f64::max);
    let half = 0.5 * SAMPLE_SPACINGS * d_max;
    // σ0 ごとの Σ w σt, Σ w σel, Σ w σf, Σ w σc, Σ w, Σ w² σt, Σ w²（w = 1/(σt + σ0)）
    let mut sums = vec![[0.0; 7]; sigma0.len()];
    let mut mean = [0.0; 4];
    let mut count = 0usize;
    let mut floored = 0usize;
    // 弾性散乱、核分裂、捕獲の外の反応の背景（MF3 の全断面積との差）
    let other = background.total - background.elastic - background.fission - background.capture;
    let mut energies = vec![0.0; sampling.samples];
    let mut parts = vec![[0.0; 3]; sampling.samples];
    for _ in 0..sampling.ladders {
        // 等間隔の点を一様な乱数でずらす（系統抽出）。どの点も窓の中で一様なので偏りがなく、共鳴の間の谷を独立な一様の点より
        // 漏れなく拾う
        let offset = rng.uniform();
        let step = 2.0 * half / sampling.samples as f64;
        for (i, x) in energies.iter_mut().enumerate() {
            *x = -half + step * (i as f64 + offset);
        }
        parts.iter_mut().for_each(|p| *p = [0.0; 3]);
        for q in sequences {
            add_ladder(q, half, &mut rng, &energies, &mut parts);
        }
        for (p, &x) in parts.iter().zip(&energies) {
            let mut elastic = background.elastic + potential + p[0];
            // |1 − U|² は負にならないので、負になるのは MF3 の負の背景を足したときだけ。NJOY PURR（purr.f90 の unrest）と同じく
            // 1 μb にして、参照解の確率テーブルと扱いを揃える（ADR-0177）
            if elastic < 0.0 && background.elastic < 0.0 {
                elastic = ELASTIC_FLOOR;
                floored += 1;
            }
            let fission = background.fission + p[1];
            let capture = background.capture + p[2];
            let total = elastic + fission + capture + other;
            if total.is_nan() || total <= 0.0 {
                return Err(x);
            }
            for (sum, &s0) in sums.iter_mut().zip(sigma0) {
                let w = 1.0 / (total + s0);
                sum[0] += w * total;
                sum[1] += w * elastic;
                sum[2] += w * fission;
                sum[3] += w * capture;
                sum[4] += w;
                sum[5] += w * w * total;
                sum[6] += w * w;
            }
            mean[0] += total;
            mean[1] += elastic;
            mean[2] += fission;
            mean[3] += capture;
            count += 1;
        }
    }
    let m = mean.map(|v| v / count as f64);
    Ok(Bondarenko {
        shielded: sums
            .iter()
            .map(|s| ShieldedSigma {
                total: s[0] / s[4],
                elastic: s[1] / s[4],
                fission: s[2] / s[4],
                capture: s[3] / s[4],
                current_total: s[5] / s[6],
            })
            .collect(),
        mean: Sigma {
            total: m[0],
            elastic: m[1],
            fission: m[2],
            capture: m[3],
        },
        floored,
    })
}

/// 系列 `q` のラダーを 1 つ作り、エネルギー `energies`（窓 [−half, half] の中）での共鳴の断面積
/// [弾性散乱, 核分裂, 捕獲] を `parts` に足す。平均断面積と同じく、Γn か Γγ が 0 以下の系列は共鳴を持たない。
///
/// 捕獲と核分裂は単準位 Breit-Wigner で、ψ = √π y Re w（PURR の unrest と同じ）。弾性散乱は、同じ系列（同じ J）の共鳴どうしの
/// 干渉を多準位 Breit-Wigner の |1 − U|² で入れる。単準位の和（PURR）は干渉を共鳴ごとに数えるので、強い s 波の共鳴の干渉の谷で
/// 全断面積が負になる（ADR-0177）。|1 − U|² は 0 K で負にならず、ガウスの核で広げても負にならない。
///
/// 1 − U = A + B S、A = 1 − e^{−2iφ}、B = −i e^{−2iφ}、S = Σ_r Γn_r f_r、f_r = 1/(p_r − E)、p_r = E_r − iΓ_r/2 として、
/// |1 − U|² = |A|² + Σ_r 2 Re[a_r f_r]、a_r = Γn_r (C_r + i A)、C_r = Σ_s Γn_s/(p̄_s − p_r)（f_r f̄_s を部分分数に分けた）。
/// |A|² = 4 sin²φ はポテンシャル散乱で、系列の外で足す。ドップラーで広げた f_r は (i√π/Δ) w((E − E_r + iΓ_r/2)/Δ)
/// （w(z) = (i/π) ∫ e^{−t²}/(z − t) dt から）。σ = (π/k²) g |1 − U|² で、π/k² g = scale/4
fn add_ladder(
    q: &LadderSequence,
    half: f64,
    rng: &mut Rng,
    energies: &[f64],
    parts: &mut [[f64; 3]],
) {
    if q.gn <= 0.0 || q.gg <= 0.0 {
        return;
    }
    let mean_width = q.gn + q.gf + q.gg + q.gx;
    let margin =
        (mean_width / (PI * q.d * TAIL)).clamp(MIN_MARGIN_SPACINGS, MAX_MARGIN_SPACINGS) * q.d;
    let (low, high) = (-half - margin, half + margin);
    // 共鳴: (E_r, Γn, Γf, Γ)。最初の共鳴の位置は、最初の間隔の幅の中で一様（PURR の ladr2 と同じ）
    let mut resonances = Vec::new();
    let mut er = low + q.d * (4.0 / PI).sqrt() * rng.uniform();
    while er < high {
        let gn = q.gn * chi_square(q.dof_n, rng);
        let gf = q.gf * chi_square(q.dof_f, rng);
        let gx = q.gx * chi_square(q.dof_x, rng);
        resonances.push((er, gn, gf, gn + gf + q.gg + gx));
        er += wigner_spacing(q.d, rng);
    }
    let a = Complex64::new(1.0 - q.cos2, q.sin2p);
    let pole = |(e, _, _, g): (f64, f64, f64, f64)| Complex64::new(e, -0.5 * g);
    let coefficients: Vec<Complex64> = resonances
        .iter()
        .map(|&r| {
            let p = pole(r);
            let c: Complex64 = resonances.iter().map(|&s| s.1 / (pole(s).conj() - p)).sum();
            r.1 * (c + Complex64::i() * a)
        })
        .collect();
    let root_pi = PI.sqrt();
    // 弾性散乱 (π/k²) g 2 Re[a_r (i√π/Δ) w] = scale/2 · √π/Δ · Re[i a_r w]
    let elastic_scale = 0.5 * q.scale * root_pi / q.delta;
    for (&(er, gn, gf, gamma), &ar) in resonances.iter().zip(&coefficients) {
        let y = 0.5 * gamma / q.delta;
        // 捕獲 scale (Γn Γγ/Γ²) ψ、核分裂 scale (Γn Γf/Γ²) ψ
        let c = q.scale * gn / gamma * root_pi * y / gamma;
        let (fission_re, capture_re) = (c * gf, c * q.gg);
        let iar = Complex64::i() * ar * elastic_scale;
        for (p, &e) in parts.iter_mut().zip(energies) {
            let x = (e - er) / q.delta;
            let w = if x.abs() > ASYMPTOTIC || y > ASYMPTOTIC {
                // w(z) ≈ i/(√π z)
                Complex64::i() / (root_pi * Complex64::new(x, y))
            } else {
                faddeeva(Complex64::new(x, y))
            };
            p[0] += (iar * w).re;
            p[1] += fission_re * w.re;
            p[2] += capture_re * w.re;
        }
    }
}

/// 自由度 `dof` の χ² 分布を自由度で割った値（平均 1）。自由度なしは 1。
pub(super) fn chi_square(dof: WidthDof, rng: &mut Rng) -> f64 {
    let nu = dof.number();
    if nu == 0 {
        return 1.0;
    }
    let sum: f64 = (0..nu).map(|_| rng.normal().powi(2)).sum();
    sum / f64::from(nu)
}

/// 平均 `d` の Wigner 分布の準位間隔 s = d √(−4 ln u / π)（p(s) = (π s/2d²) exp(−π s²/(4d²)) の累積分布の逆関数）。
pub(super) fn wigner_spacing(d: f64, rng: &mut Rng) -> f64 {
    d * (-4.0 * rng.uniform().ln() / PI).sqrt()
}

/// 決定的な擬似乱数（xoshiro256**。種は SplitMix64 で広げる）。外部の依存を足さずに、種から同じ列を作る。
pub(super) struct Rng {
    state: [u64; 4],
    /// Box-Muller の 2 つ目の正規乱数
    spare: Option<f64>,
}

impl Rng {
    pub(super) fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            x ^ (x >> 31)
        };
        Self {
            state: [next(), next(), next(), next()],
            spare: None,
        }
    }

    fn next_u64(&mut self) -> u64 {
        let s = &mut self.state;
        let result = s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// (0, 1) の一様乱数（53 ビットの格子の中点で、0 と 1 を含まない）。
    pub(super) fn uniform(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64 + 0.5) * (1.0 / (1u64 << 53) as f64)
    }

    /// 標準正規乱数（Box-Muller）。
    pub(super) fn normal(&mut self) -> f64 {
        if let Some(v) = self.spare.take() {
            return v;
        }
        let (u, v) = (self.uniform(), self.uniform());
        let r = (-2.0 * u.ln()).sqrt();
        let a = 2.0 * PI * v;
        self.spare = Some(r * a.sin());
        r * a.cos()
    }
}

#[cfg(test)]
mod tests;
