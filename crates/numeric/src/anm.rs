//! 解析的ノード法（ANM）の1ノードの計算（ADR-0102）。核分裂のないノードを、core と lattice の ADF の
//! 1ノードの問題（ADR-0094）が同じ式で解く。
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
    #[error("the analytic nodal face conditions are singular")]
    Singular,
}

/// 群 `group` の、出ていく部分中性子流のノード平均についての線形の式（`nem::Outflow`。3次・4次の係数は 0）。
/// ほかの群のノード平均と入射の条件は `node` の値を使う。
pub fn outflow(node: &AnmNode, group: usize) -> Result<Outflow, AnmError> {
    let ng = node.diffusion.len();
    if ng != 2 {
        return Err(AnmError::Groups(ng));
    }
    let lengths = [
        node.total.len(),
        node.scatter.len(),
        node.transverse.len(),
        node.inflow.len(),
        node.averages.len(),
    ];
    if lengths.iter().any(|&l| l != 2) || node.scatter.iter().any(|row| row.len() != 2) {
        return Err(AnmError::Invalid(
            "ANM constants, leakages, inflows and averages must have 2 groups".into(),
        ));
    }
    let h = node.width;
    if !(h > 0.0 && h.is_finite()) {
        return Err(AnmError::Invalid(format!(
            "ANM width must be positive, got {h}"
        )));
    }
    if !node.diffusion.iter().all(|&d| d > 0.0 && d.is_finite()) {
        return Err(AnmError::Invalid(
            "ANM diffusion coefficients must be positive".into(),
        ));
    }
    let d = node.diffusion;
    // A[g][g'] と M = D⁻¹ A
    let a: [[f64; 2]; 2] = std::array::from_fn(|g| {
        std::array::from_fn(|from| {
            if from == g {
                node.total[g] - node.scatter[g][g]
            } else {
                -node.scatter[from][g]
            }
        })
    });
    if !(a[0][0] > 0.0 && a[1][1] > 0.0) {
        return Err(AnmError::Modes(format!(
            "ANM removal cross sections must be positive, got {} and {}",
            a[0][0], a[1][1]
        )));
    }
    let m: [[f64; 2]; 2] = std::array::from_fn(|g| std::array::from_fn(|h| a[g][h] / d[g]));
    let modes = modes(m)?;
    let det_a = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let inverse = |r: [f64; 2]| {
        [
            (a[1][1] * r[0] - a[0][1] * r[1]) / det_a,
            (a[0][0] * r[1] - a[1][0] * r[0]) / det_a,
        ]
    };
    let p1 = inverse([-node.transverse[0].0, -node.transverse[1].0]);
    let p2 = inverse([-node.transverse[0].1, -node.transverse[1].1]);

    // 面（+, -）の中性子束と dφ/dξ の、係数 u = (A_a, B_a, A_b, B_b) についての係数
    // [群][面][未知数]
    let mut phi = [[[0.0; 4]; 2]; 2];
    let mut slope = [[[0.0; 4]; 2]; 2];
    for (k, (mu, v)) in modes.iter().enumerate() {
        let t = 0.5 * mu.sqrt() * h;
        let tanh = t.tanh();
        // C(±1/2) = 1 − tanh(t)/t、S(±1/2) = ±tanh t、dC/dξ(±1/2) = ±2t tanh t、dS/dξ(±1/2) = 2t
        let even = 1.0 - tanh / t;
        for g in 0..2 {
            for (f, sign) in [(0, 1.0), (1, -1.0)] {
                phi[g][f][2 * k] = v[g] * even;
                phi[g][f][2 * k + 1] = v[g] * sign * tanh;
                slope[g][f][2 * k] = v[g] * sign * 2.0 * t * tanh;
                slope[g][f][2 * k + 1] = v[g] * 2.0 * t;
            }
        }
    }
    // 面の条件の行列（[群][面] の行）: + 面 (1 − β) φ/4 − (1 + β) J/2 = e、− 面 (1 − β) φ/4 + (1 + β) J/2 = e、
    // J = −(D/h) dφ/dξ
    let mut matrix = [[0.0; 4]; 4];
    let row = |g: usize, f: usize| 2 * g + f;
    let weights = |g: usize, f: usize| {
        let Inflow { beta, .. } = node.inflow[g][f];
        let sign = if f == 0 { -1.0 } else { 1.0 };
        // φ と dφ/dξ にかける係数
        ((1.0 - beta) / 4.0, -sign * (1.0 + beta) / 2.0 * d[g] / h)
    };
    for g in 0..2 {
        for f in 0..2 {
            let (wp, ws) = weights(g, f);
            for j in 0..4 {
                matrix[row(g, f)][j] = wp * phi[g][f][j] + ws * slope[g][f][j];
            }
        }
    }
    // 既知の部分（ノード平均 `averages`、特解、外からの入射）の面の値から右辺を作り、解いて、群 `group` の
    // 出ていく部分中性子流 o+ = φ/4 + J/2、o− = φ/4 − J/2 を返す
    let solve_out = |averages: [f64; 2], particular: bool| -> Result<[f64; 2], AnmError> {
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
                    node.inflow[g][f].external
                } else {
                    0.0
                };
                rhs[row(g, f)] = external - wp * p - ws * s;
            }
        }
        let u = solve4(matrix, rhs).ok_or(AnmError::Singular)?;
        Ok(std::array::from_fn(|f| {
            let (p, s) = known(group, f);
            let value = p + (0..4).map(|j| phi[group][f][j] * u[j]).sum::<f64>();
            let derivative = s + (0..4).map(|j| slope[group][f][j] * u[j]).sum::<f64>();
            let current = -d[group] / h * derivative;
            let sign = if f == 0 { 1.0 } else { -1.0 };
            value / 4.0 + sign * current / 2.0
        }))
    };
    let mut others = [node.averages[0], node.averages[1]];
    others[group] = 0.0;
    let mut unit = [0.0; 2];
    unit[group] = 1.0;
    Ok(Outflow {
        out: solve_out(others, true)?,
        slope: solve_out(unit, false)?,
        beta: [node.inflow[group][0].beta, node.inflow[group][1].beta],
        external: [
            node.inflow[group][0].external,
            node.inflow[group][1].external,
        ],
        higher: [0.0; 2],
        higher_slope: [0.0; 2],
    })
}

/// D⁻¹A の2つの固有値 κ² と固有ベクトル。どちらも正の実数で、互いに異なること。
fn modes(m: [[f64; 2]; 2]) -> Result<[(f64, [f64; 2]); 2], AnmError> {
    let half = 0.5 * (m[0][0] + m[1][1]);
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let discriminant = half * half - det;
    let distinct = 1.0e-12 * half * half;
    if !(discriminant.is_finite() && det.is_finite()) || discriminant < -distinct {
        return Err(AnmError::Modes(format!(
            "D^-1 A has complex eigenvalues (matrix {m:?})"
        )));
    }
    let coupled = m[0][1] != 0.0 || m[1][0] != 0.0;
    if discriminant <= distinct && coupled {
        return Err(AnmError::Modes(format!(
            "D^-1 A has a repeated eigenvalue with coupled groups (matrix {m:?})"
        )));
    }
    let root = discriminant.max(0.0).sqrt();
    let values = [half + root, half - root];
    if values.iter().any(|&mu| mu.is_nan() || mu <= 0.0) {
        return Err(AnmError::Modes(format!(
            "D^-1 A must have positive eigenvalues, got {values:?}"
        )));
    }
    if !coupled {
        // 群が結合しないなら、固有ベクトルは単位ベクトル（固有値が等しくてもよい）
        return Ok([(m[0][0], [1.0, 0.0]), (m[1][1], [0.0, 1.0])]);
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

/// 4×4 の連立一次方程式（部分ピボットの Gauss の消去法）。特異なら `None`。
fn solve4(mut a: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    for col in 0..4 {
        let pivot = (col..4).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col] == 0.0 || a[pivot][col].is_nan() {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        let pivot_row = a[col];
        for r in col + 1..4 {
            let factor = a[r][col] / pivot_row[col];
            for (value, p) in a[r][col..].iter_mut().zip(&pivot_row[col..]) {
                *value -= factor * p;
            }
            b[r] -= factor * b[col];
        }
    }
    let mut x = [0.0; 4];
    for r in (0..4).rev() {
        let sum: f64 = (r + 1..4).map(|c| a[r][c] * x[c]).sum();
        x[r] = (b[r] - sum) / a[r][r];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nem::Inflow;

    /// 反射体（水と鋼）に近い2群の定数: 拡散係数、全断面積、散乱 [出る群][入る群]（上方散乱あり）
    fn reflector() -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
        (
            vec![1.25, 0.30],
            vec![0.60, 1.90],
            vec![vec![0.565, 0.032], vec![0.0015, 1.885]],
        )
    }

    /// 細かいメッシュ（セル中心の差分、`cells` 個）の参照解: -D φ'' + A φ = -L(x)
    /// （L = L̄ + l1 ξ + l2 (3ξ² - 1/4)、ξ = x/h - 1/2）。右の面は入射 `right`、左の面は入射 `left`
    /// （`None` なら反射）。[群] の（ノード平均、右の面の出ていく部分中性子流、左の面の出ていく部分中性子流）を返す
    fn finite_difference(
        width: f64,
        cells: usize,
        transverse: &[(f64, f64, f64)],
        right: &[f64],
        left: Option<&[f64]>,
    ) -> Vec<(f64, f64, f64)> {
        let (diffusion, total, scatter) = reflector();
        let dx = width / cells as f64;
        let mut flux = vec![vec![0.0; cells]; 2];
        // 入射 i の面: 面の中性子束 φf で φf/4 - J_in/2 = i、外向きの J = k(φ_端 - φf)、k = 2D/dx から
        // φf = (i + k φ_端/2) / (1/4 + k/2)
        let face = |d: f64, end: f64, incoming: f64| {
            let k = 2.0 * d / dx;
            (incoming + k / 2.0 * end) / (0.25 + k / 2.0)
        };
        for _ in 0..200 {
            let mut change: f64 = 0.0;
            for g in 0..2 {
                let other = 1 - g;
                let removal = total[g] - scatter[g][g];
                let d = diffusion[g];
                let k = 2.0 * d / dx;
                let f = 1.0 / (0.25 + k / 2.0);
                // 三重対角 a_i φ_{i-1} + b_i φ_i + c_i φ_{i+1} = r_i
                let mut a = vec![-d / (dx * dx); cells];
                let mut b = vec![2.0 * d / (dx * dx) + removal; cells];
                let mut c = vec![-d / (dx * dx); cells];
                let mut r: Vec<f64> = (0..cells)
                    .map(|i| {
                        let xi = (i as f64 + 0.5) / cells as f64 - 0.5;
                        let (mean, l1, l2) = transverse[g];
                        scatter[other][g] * flux[other][i]
                            - (mean + l1 * xi + l2 * (3.0 * xi * xi - 0.25))
                    })
                    .collect();
                // 入射の面: 外向きの J = k φ_端 (1 - f k/2) - k f i
                let open = k * (1.0 - f * k / 2.0) / dx;
                a[0] = 0.0;
                b[0] -= d / (dx * dx);
                if let Some(left) = left {
                    b[0] += open;
                    r[0] += k * f * left[g] / dx;
                }
                c[cells - 1] = 0.0;
                b[cells - 1] += -d / (dx * dx) + open;
                r[cells - 1] += k * f * right[g] / dx;
                // Thomas 法
                for i in 1..cells {
                    let m = a[i] / b[i - 1];
                    b[i] -= m * c[i - 1];
                    r[i] -= m * r[i - 1];
                }
                let mut new = vec![0.0; cells];
                new[cells - 1] = r[cells - 1] / b[cells - 1];
                for i in (0..cells - 1).rev() {
                    new[i] = (r[i] - c[i] * new[i + 1]) / b[i];
                }
                for (n, o) in new.iter().zip(&flux[g]) {
                    change = change.max((n - o).abs() / n.abs().max(1e-300));
                }
                flux[g] = new;
            }
            if change < 1e-15 {
                break;
            }
        }
        (0..2)
            .map(|g| {
                let d = diffusion[g];
                let k = 2.0 * d / dx;
                let average = flux[g].iter().sum::<f64>() / cells as f64;
                // 出ていく部分中性子流 φf/4 + J_out/2
                let out = |end: f64, incoming: f64| {
                    let phi = face(d, end, incoming);
                    phi / 4.0 + k * (end - phi) / 2.0
                };
                let right_out = out(flux[g][cells - 1], right[g]);
                let left_out = match left {
                    Some(left) => out(flux[g][0], left[g]),
                    // 反射: 正味の中性子流 0、面の中性子束は外挿（2次の精度）
                    None => (1.5 * flux[g][0] - 0.5 * flux[g][1]) / 4.0,
                };
                (average, right_out, left_out)
            })
            .collect()
    }

    #[test]
    fn non_fissile_slab_matches_fine_differences() {
        // 核分裂のない2群の平板（20 cm）。左は反射、右から入射。横方向の漏れの2次式を含む。ノード平均を
        // 参照の値にすると、出ていく部分中性子流は参照に一致し（ANM は1D の式を厳密に解く）、1D のバランス
        // （正味の漏れ + 除去 - 群の間の散乱 + 横方向の漏れの平均 = 0）も満たす
        let (diffusion, total, scatter) = reflector();
        let width = 20.0;
        let transverse = [(1.0e-4, -2.0e-4, 5.0e-5), (-3.0e-5, 4.0e-5, -1.0e-5)];
        let incoming = [0.02, 0.006];
        let reference = finite_difference(width, 20_000, &transverse, &incoming, None);
        let averages: Vec<f64> = reference.iter().map(|r| r.0).collect();
        let moments: Vec<(f64, f64)> = transverse.iter().map(|t| (t.1, t.2)).collect();
        let inflow: Vec<[Inflow; 2]> = (0..2)
            .map(|g| {
                [
                    Inflow {
                        beta: 0.0,
                        external: incoming[g],
                    },
                    Inflow {
                        beta: 1.0,
                        external: 0.0,
                    },
                ]
            })
            .collect();
        let node = AnmNode {
            diffusion: &diffusion,
            total: &total,
            scatter: &scatter,
            width,
            transverse: &moments,
            inflow: &inflow,
            averages: &averages,
        };
        for g in 0..2 {
            let faces = outflow(&node, g).unwrap().at(averages[g]);
            let (_, right, left) = reference[g];
            assert!(
                (faces.out[0] - right).abs() < 2e-6 * right,
                "group {g} right: {} vs {right}",
                faces.out[0]
            );
            assert!(
                (faces.out[1] - left).abs() < 2e-6 * left,
                "group {g} left: {} vs {left}",
                faces.out[1]
            );
            // 1D のバランス
            let net = (faces.out[0] - faces.incoming[0] + faces.out[1] - faces.incoming[1]) / width;
            let other = 1 - g;
            let balance = net + (total[g] - scatter[g][g]) * averages[g]
                - scatter[other][g] * averages[other]
                + transverse[g].0;
            assert!(balance.abs() < 1e-6 * net.abs(), "group {g}: {balance}");
        }
    }

    #[test]
    fn steep_reflector_keeps_positive_face_flux() {
        // 反射体の I 型（ADR-0102 の背景）のように、燃料側から入り真空側へ出るときの 21.5 cm のノード。
        // 面の中性子束は正（NEM の4次の多項式では負になる形）
        let (diffusion, total, scatter) = reflector();
        let width = 21.5;
        let incoming = [0.05, 0.004];
        let reference = finite_difference(
            width,
            20_000,
            &[(0.0, 0.0, 0.0); 2],
            &incoming,
            Some(&[0.0; 2]),
        );
        let averages: Vec<f64> = reference.iter().map(|r| r.0).collect();
        let inflow: Vec<[Inflow; 2]> = (0..2)
            .map(|g| {
                [
                    Inflow {
                        beta: 0.0,
                        external: incoming[g],
                    },
                    Inflow {
                        beta: 0.0,
                        external: 0.0,
                    },
                ]
            })
            .collect();
        let node = AnmNode {
            diffusion: &diffusion,
            total: &total,
            scatter: &scatter,
            width,
            transverse: &[(0.0, 0.0); 2],
            inflow: &inflow,
            averages: &averages,
        };
        for g in 0..2 {
            let faces = outflow(&node, g).unwrap().at(averages[g]);
            let (_, right, left) = reference[g];
            assert!((faces.out[0] - right).abs() < 2e-6 * right, "group {g}");
            assert!((faces.out[1] - left).abs() < 2e-5 * left, "group {g}");
            let phi = 2.0 * (faces.out[1] + faces.incoming[1]);
            assert!(phi > 0.0, "group {g}: {phi}");
        }
    }

    #[test]
    fn invalid_anm_nodes_are_errors() {
        let (diffusion, total, scatter) = reflector();
        let inflow = vec![[Inflow::default(); 2]; 2];
        let averages = [1.0, 1.0];
        let node = AnmNode {
            diffusion: &diffusion,
            total: &total,
            scatter: &scatter,
            width: 20.0,
            transverse: &[(0.0, 0.0); 2],
            inflow: &inflow,
            averages: &averages,
        };
        assert!(outflow(&node, 0).is_ok());
        // 2群でない
        let one = AnmNode {
            diffusion: &diffusion[..1],
            total: &total[..1],
            scatter: &[vec![0.5]],
            transverse: &[(0.0, 0.0)],
            inflow: &inflow[..1],
            averages: &averages[..1],
            ..node
        };
        assert!(matches!(outflow(&one, 0), Err(AnmError::Groups(1))));
        // 除去が負（固有値が正でない）
        let negative = [vec![0.7, 0.032], vec![0.0015, 1.885]];
        let bad = AnmNode {
            scatter: &negative,
            ..node
        };
        assert!(matches!(outflow(&bad, 0), Err(AnmError::Modes(_))));
        // 幅が正でない
        let zero = AnmNode { width: 0.0, ..node };
        assert!(matches!(outflow(&zero, 0), Err(AnmError::Invalid(_))));
    }
}
