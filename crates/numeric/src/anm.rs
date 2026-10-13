//! 解析的ノード法（ANM）の1ノードの計算（ADR-0102）。核分裂のないノードを、core と lattice の ADF の
//! 1ノードの問題（ADR-0103）が同じ式で解く。
//!
//! 横方向に積分した2群の1D の式 −D φ″ + A φ = −L(x)（A[g][g] = Σr_g、A[g][g'] = −Σs(g'→g)、L は横方向の漏れ）
//! を、ξ = x/h ∈ [−1/2, 1/2] で
//!
//!   φ(ξ) = φ̄ + p1 f1(ξ) + p2 f2(ξ) + Σ_k v_k (A_k C_k(ξ) + B_k S_k(ξ))
//!
//! と厳密に表す。f1 = ξ、f2 = 3ξ² − 1/4 は NEM と同じ横方向の漏れの1次・2次の項（L = L̄ + l1 f1 + l2 f2）で、
//! 特解は A p1 = −l1、A p2 = −l2。κ_k² と v_k は D⁻¹A の固有値と固有ベクトルで、C_k = (cosh(κ_k h ξ) − 平均) と
//! S_k = sinh(κ_k h ξ)（どちらも cosh(κ_k h/2) で割る）はノード平均が 0。ノード平均 φ̄ は NEM と同じく引数とし、
//! 0 次のモーメント（バランス）は 3D のバランス式で課す。4 つの係数 A_k, B_k を、2群の両方の面の入射の条件
//! i = β o + 外からの入射 から解く。

use crate::nem::{Inflow, Outflow};

/// 1つのノードの1つの方向の、核分裂のない2群の入力。
#[derive(Debug, Clone, Copy)]
pub struct AnmNode<'a> {
    /// [群] の拡散係数
    pub diffusion: &'a [f64],
    /// [群] の全断面積
    pub total: &'a [f64],
    /// [出る群][入る群] の散乱
    pub scatter: &'a [Vec<f64>],
    /// その方向のノードの幅 (cm)
    pub width: f64,
    /// [群] の横方向の漏れの f1, f2 の係数（`nem::transverse_moments`）
    pub transverse: &'a [(f64, f64)],
    /// [群][+ 面, - 面] の入射の条件
    pub inflow: &'a [[Inflow; 2]],
    /// [群] のノード平均（求める群の値は使わない）
    pub averages: &'a [f64],
}

/// ANM の核のエラー。
#[derive(Debug, thiserror::Error)]
pub enum AnmError {
    #[error("the analytic nodal method needs 2 groups, got {0}")]
    Groups(usize),
    #[error("{0}")]
    Modes(String),
    #[error("{0}")]
    Invalid(String),
    /// 面の条件の式（列をそろえた 4×4）の条件数が `MAX_CONDITION` を超える。群が結合して D⁻¹A の固有値が
    /// 重なる（近い）と、2つのモードの固有ベクトルがそろって起きる
    #[error(
        "the analytic nodal face conditions are ill-conditioned (condition number {condition:e}); \
         the eigenvalues of D^-1 A are (nearly) repeated with coupled groups"
    )]
    IllConditioned { condition: f64 },
}

/// 面の条件の式の条件数の上限。解の相対誤差はおよそ 条件数 × ε なので、ε の 10⁻⁸ 倍の逆数にして、8 桁の
/// 精度を保てないときはエラーにする（黙って精度を落とさない）
const MAX_CONDITION: f64 = 1.0e-8 / f64::EPSILON;

/// 1 − tanh(t)/t を級数で求める t の上限。直接の式は 1 から引く桁落ちで相対誤差がおよそ 3ε/t²、級数
/// t²/3 − 2t⁴/15 + 17t⁶/315 − 62t⁸/2835 の打ち切りの相対誤差は、次の項 1382t¹⁰/155925 から 0.0266 t⁸。
/// 両者が等しくなる t = (3ε/0.0266)^(1/10) ≈ 0.0436 で切り替え、どちらの側でも相対誤差を 4 × 10⁻¹³ 以下にする
const SERIES_BELOW: f64 = 0.0436;

/// cosh のモードの面の値（cosh(κh/2) で割り、ノード平均を引いたもの）1 − tanh(t)/t（t = κh/2）
fn even_face_value(t: f64) -> f64 {
    if t < SERIES_BELOW {
        let t2 = t * t;
        t2 * (1.0 / 3.0 - t2 * (2.0 / 15.0 - t2 * (17.0 / 315.0 - t2 * 62.0 / 2835.0)))
    } else {
        1.0 - t.tanh() / t
    }
}

/// 核分裂のない2群のノードの、1つの方向の定数で決まる量（固有値と固有ベクトル、面の値の表）。core はノードと
/// 方向ごとに1回作り、反復のたびに `outflow` で使う。
#[derive(Debug, Clone, PartialEq)]
pub struct AnmModes {
    diffusion: [f64; 2],
    width: f64,
    /// A[g][g']
    removal: [[f64; 2]; 2],
    /// [群][面（+, -）][未知数 (A_a, B_a, A_b, B_b)] の、面の中性子束と dφ/dξ への係数
    phi: [[[f64; 4]; 2]; 2],
    slope: [[[f64; 4]; 2]; 2],
}

impl AnmModes {
    /// 拡散係数、全断面積、散乱 [出る群][入る群]、幅から作る。2群で、値が有限で、除去が正で、D⁻¹A の固有値が
    /// 正の実数であること。
    pub fn new(
        diffusion: &[f64],
        total: &[f64],
        scatter: &[Vec<f64>],
        width: f64,
    ) -> Result<Self, AnmError> {
        let ng = diffusion.len();
        if ng != 2 {
            return Err(AnmError::Groups(ng));
        }
        if total.len() != 2 || scatter.len() != 2 || scatter.iter().any(|row| row.len() != 2) {
            return Err(AnmError::Invalid("ANM constants must have 2 groups".into()));
        }
        if !(width > 0.0 && width.is_finite()) {
            return Err(AnmError::Invalid(format!(
                "ANM width must be positive, got {width}"
            )));
        }
        if !diffusion.iter().all(|&d| d > 0.0 && d.is_finite()) {
            return Err(AnmError::Invalid(
                "ANM diffusion coefficients must be positive".into(),
            ));
        }
        if !(total.iter().all(|v| v.is_finite()) && scatter.iter().flatten().all(|v| v.is_finite()))
        {
            return Err(AnmError::Invalid(
                "ANM cross sections must be finite".into(),
            ));
        }
        let a: [[f64; 2]; 2] = std::array::from_fn(|g| {
            std::array::from_fn(|from| {
                if from == g {
                    total[g] - scatter[g][g]
                } else {
                    -scatter[from][g]
                }
            })
        });
        if !(a[0][0] > 0.0 && a[1][1] > 0.0) {
            return Err(AnmError::Modes(format!(
                "ANM removal cross sections must be positive, got {} and {}",
                a[0][0], a[1][1]
            )));
        }
        let m: [[f64; 2]; 2] =
            std::array::from_fn(|g| std::array::from_fn(|h| a[g][h] / diffusion[g]));
        let mut phi = [[[0.0; 4]; 2]; 2];
        let mut slope = [[[0.0; 4]; 2]; 2];
        for (k, (mu, v)) in modes(m)?.iter().enumerate() {
            let t = 0.5 * mu.sqrt() * width;
            let tanh = t.tanh();
            // C(±1/2) = 1 − tanh(t)/t、S(±1/2) = ±tanh t、dC/dξ(±1/2) = ±2t tanh t、dS/dξ(±1/2) = 2t
            let even = even_face_value(t);
            for g in 0..2 {
                for (f, sign) in [(0, 1.0), (1, -1.0)] {
                    phi[g][f][2 * k] = v[g] * even;
                    phi[g][f][2 * k + 1] = v[g] * sign * tanh;
                    slope[g][f][2 * k] = v[g] * sign * 2.0 * t * tanh;
                    slope[g][f][2 * k + 1] = v[g] * 2.0 * t;
                }
            }
        }
        Ok(Self {
            diffusion: [diffusion[0], diffusion[1]],
            width,
            removal: a,
            phi,
            slope,
        })
    }

    /// 群 `group` の、出ていく部分中性子流のノード平均についての線形の式（`nem::Outflow`。3次・4次の係数は 0）。
    /// `transverse` は [群] の横方向の漏れの f1, f2 の係数、`inflow` は [群][+ 面, - 面] の入射の条件、`averages` は
    /// [群] のノード平均（求める群の値は使わない）。
    pub fn outflow(
        &self,
        transverse: &[(f64, f64)],
        inflow: &[[Inflow; 2]],
        averages: &[f64],
        group: usize,
    ) -> Result<Outflow, AnmError> {
        if transverse.len() != 2 || inflow.len() != 2 || averages.len() != 2 || group >= 2 {
            return Err(AnmError::Invalid(
                "ANM leakages, inflows and averages must have 2 groups".into(),
            ));
        }
        let finite = transverse
            .iter()
            .all(|(a, b)| a.is_finite() && b.is_finite())
            && inflow
                .iter()
                .flatten()
                .all(|i| i.beta.is_finite() && i.external.is_finite())
            && averages.iter().all(|v| v.is_finite());
        if !finite {
            return Err(AnmError::Invalid(
                "ANM leakages, inflows and averages must be finite".into(),
            ));
        }
        let (d, h, a) = (self.diffusion, self.width, self.removal);
        let det_a = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        let inverse = |r: [f64; 2]| {
            [
                (a[1][1] * r[0] - a[0][1] * r[1]) / det_a,
                (a[0][0] * r[1] - a[1][0] * r[0]) / det_a,
            ]
        };
        let p1 = inverse([-transverse[0].0, -transverse[1].0]);
        let p2 = inverse([-transverse[0].1, -transverse[1].1]);
        let (phi, slope) = (&self.phi, &self.slope);
        // 面の条件の行列（[群][面] の行）: + 面 (1 − β) φ/4 − (1 + β) J/2 = e、− 面 (1 − β) φ/4 + (1 + β) J/2 = e、
        // J = −(D/h) dφ/dξ
        let row = |g: usize, f: usize| 2 * g + f;
        let weights = |g: usize, f: usize| {
            let Inflow { beta, .. } = inflow[g][f];
            let sign = if f == 0 { -1.0 } else { 1.0 };
            // φ と dφ/dξ にかける係数
            ((1.0 - beta) / 4.0, -sign * (1.0 + beta) / 2.0 * d[g] / h)
        };
        let mut matrix = [[0.0; 4]; 4];
        for g in 0..2 {
            for f in 0..2 {
                let (wp, ws) = weights(g, f);
                for j in 0..4 {
                    matrix[row(g, f)][j] = wp * phi[g][f][j] + ws * slope[g][f][j];
                }
            }
        }
        let lu = Lu4::conditioned(matrix)?;
        // 既知の部分（ノード平均 `averages`、特解、外からの入射）の面の値から右辺を作り、解いて、群 `group` の
        // 出ていく部分中性子流 o+ = φ/4 + J/2、o− = φ/4 − J/2 を返す
        let solve_out = |averages: [f64; 2], particular: bool| -> [f64; 2] {
            let known = |g: usize, f: usize| -> (f64, f64) {
                let sign = if f == 0 { 1.0 } else { -1.0 };
                let (q1, q2) = if particular {
                    (p1[g], p2[g])
                } else {
                    (0.0, 0.0)
                };
                // f1(±1/2) = ±1/2、f2(±1/2) = 1/2、f1' = 1、f2'(±1/2) = ±3
                (
                    averages[g] + sign * q1 / 2.0 + q2 / 2.0,
                    q1 + sign * 3.0 * q2,
                )
            };
            let mut rhs = [0.0; 4];
            for g in 0..2 {
                for f in 0..2 {
                    let (wp, ws) = weights(g, f);
                    let (p, s) = known(g, f);
                    let external = if particular {
                        inflow[g][f].external
                    } else {
                        0.0
                    };
                    rhs[row(g, f)] = external - wp * p - ws * s;
                }
            }
            let u = lu.solve(rhs);
            std::array::from_fn(|f| {
                let (p, s) = known(group, f);
                let value = p + (0..4).map(|j| phi[group][f][j] * u[j]).sum::<f64>();
                let derivative = s + (0..4).map(|j| slope[group][f][j] * u[j]).sum::<f64>();
                let current = -d[group] / h * derivative;
                let sign = if f == 0 { 1.0 } else { -1.0 };
                value / 4.0 + sign * current / 2.0
            })
        };
        let mut others = [averages[0], averages[1]];
        others[group] = 0.0;
        let mut unit = [0.0; 2];
        unit[group] = 1.0;
        Ok(Outflow {
            out: solve_out(others, true),
            slope: solve_out(unit, false),
            beta: [inflow[group][0].beta, inflow[group][1].beta],
            external: [inflow[group][0].external, inflow[group][1].external],
            higher: [0.0; 2],
            higher_slope: [0.0; 2],
        })
    }
}

/// 群 `group` の `AnmModes::outflow`（モードを毎回作る）。
pub fn outflow(node: &AnmNode, group: usize) -> Result<Outflow, AnmError> {
    AnmModes::new(node.diffusion, node.total, node.scatter, node.width)?.outflow(
        node.transverse,
        node.inflow,
        node.averages,
        group,
    )
}

/// D⁻¹A の2つの固有値 κ² と固有ベクトル。どちらも正の実数であること。群が結合しないときは単位ベクトルで、
/// 固有値が等しくてもよい。群が結合して固有値が重なる（近い）ときは、固有ベクトルがそろい、面の条件の式の
/// 条件数（`Lu4::conditioned`）でエラーになる。
fn modes(m: [[f64; 2]; 2]) -> Result<[(f64, [f64; 2]); 2], AnmError> {
    if m[0][1] == 0.0 && m[1][0] == 0.0 {
        return Ok([(m[0][0], [1.0, 0.0]), (m[1][1], [0.0, 1.0])]);
    }
    let half = 0.5 * (m[0][0] + m[1][1]);
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let discriminant = half * half - det;
    if !(discriminant.is_finite() && det.is_finite()) || discriminant < 0.0 {
        // 散乱の非対角の積が負でなければ判別式は (m00 − m11)²/4 + m01 m10 ≥ 0
        return Err(AnmError::Modes(format!(
            "D^-1 A has complex eigenvalues (matrix {m:?})"
        )));
    }
    let large = half + discriminant.sqrt();
    // 小さいほうは det / large（half − root の桁落ちを避ける）
    let values = [large, det / large];
    if values.iter().any(|&mu| mu.is_nan() || mu <= 0.0) {
        return Err(AnmError::Modes(format!(
            "D^-1 A must have positive eigenvalues, got {values:?}"
        )));
    }
    Ok(values.map(|mu| {
        // (M − μ) v = 0 の2つの行から、大きいほうの候補を使う
        let first = [m[0][1], mu - m[0][0]];
        let second = [mu - m[1][1], m[1][0]];
        let norm = |v: [f64; 2]| v[0].abs().max(v[1].abs());
        let v = if norm(first) >= norm(second) {
            first
        } else {
            second
        };
        let scale = norm(v);
        (mu, [v[0] / scale, v[1] / scale])
    }))
}

/// 4×4 の LU 分解（部分ピボット）。
struct Lu4 {
    lu: [[f64; 4]; 4],
    /// 列のスケール（列をそろえた行列を分解する）
    scale: [f64; 4],
    pivot: [usize; 4],
}

impl Lu4 {
    /// 列を最大の絶対値でそろえて分解し、1-ノルムの条件数が `MAX_CONDITION` 以下であることを確かめる
    /// （列のスケールは未知数の単位を変えるだけで、解の精度は列をそろえた行列の条件数で決まる）。
    fn conditioned(matrix: [[f64; 4]; 4]) -> Result<Self, AnmError> {
        let ill = |condition: f64| Err(AnmError::IllConditioned { condition });
        let mut scale = [0.0; 4];
        for (j, s) in scale.iter_mut().enumerate() {
            *s = (0..4).map(|i| matrix[i][j].abs()).fold(0.0, f64::max);
            if *s == 0.0 || !s.is_finite() {
                return ill(f64::INFINITY);
            }
        }
        let mut lu = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                lu[i][j] = matrix[i][j] / scale[j];
            }
        }
        let norm = (0..4)
            .map(|j| (0..4).map(|i| lu[i][j].abs()).sum::<f64>())
            .fold(0.0, f64::max);
        let mut pivot = [0, 1, 2, 3];
        for col in 0..4 {
            let best = (col..4)
                .max_by(|&i, &j| lu[i][col].abs().total_cmp(&lu[j][col].abs()))
                .unwrap_or(col);
            if lu[best][col] == 0.0 {
                return ill(f64::INFINITY);
            }
            lu.swap(col, best);
            pivot.swap(col, best);
            let pivot_row = lu[col];
            for row in lu.iter_mut().skip(col + 1) {
                let factor = row[col] / pivot_row[col];
                row[col] = factor;
                for (value, p) in row[col + 1..].iter_mut().zip(&pivot_row[col + 1..]) {
                    *value -= factor * p;
                }
            }
        }
        let factored = Self { lu, scale, pivot };
        // 逆行列の 1-ノルム（列ごとに単位ベクトルを解く。列をそろえた行列の逆）
        let inverse_norm = (0..4)
            .map(|j| {
                let mut e = [0.0; 4];
                e[j] = 1.0;
                factored
                    .solve_scaled(e)
                    .iter()
                    .map(|v| v.abs())
                    .sum::<f64>()
            })
            .fold(0.0, f64::max);
        let condition = norm * inverse_norm;
        if condition.is_nan() || condition > MAX_CONDITION {
            return ill(condition);
        }
        Ok(factored)
    }

    /// 列をそろえた行列の解
    fn solve_scaled(&self, b: [f64; 4]) -> [f64; 4] {
        let mut y: [f64; 4] = std::array::from_fn(|i| b[self.pivot[i]]);
        for i in 0..4 {
            for j in 0..i {
                y[i] -= self.lu[i][j] * y[j];
            }
        }
        for i in (0..4).rev() {
            for j in i + 1..4 {
                y[i] -= self.lu[i][j] * y[j];
            }
            y[i] /= self.lu[i][i];
        }
        y
    }

    /// 元の行列の解（列のスケールで割り戻す）
    fn solve(&self, b: [f64; 4]) -> [f64; 4] {
        let y = self.solve_scaled(b);
        std::array::from_fn(|j| y[j] / self.scale[j])
    }
}

#[cfg(test)]
mod tests;
