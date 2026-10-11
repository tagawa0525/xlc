use super::*;

/// 集合体の大きさのノードの半幅 (cm)
const A: f64 = 10.75;

/// 固有値 μ の範囲: 0、その近く、級数と閉じた式の境（|μ| (a² + b²) = 1）の両側、熱群のモード程度まで
fn mus() -> Vec<f64> {
    let r2 = 2.0 * A * A;
    let mut v = vec![0.0];
    for s in [1.0, -1.0] {
        for m in [1e-14, 1e-8, 1e-3, 0.999 / r2, 1.001 / r2, 0.01, 0.3] {
            v.push(s * m);
        }
    }
    v
}

/// 面の平均の Gauss-Legendre の点と重み
fn rule() -> Vec<(f64, f64)> {
    crate::quadrature::gauss_legendre(FACE_POINTS).unwrap()
}

fn points() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0),
        (A, A),
        (-A, 0.3 * A),
        (0.7 * A, -A),
        (0.21 * A, -0.83 * A),
    ]
}

#[test]
fn functions_solve_the_helmholtz_equation() {
    // ∇²ψ = μψ を、中心差分（刻み 1e-3 cm）の 2 階微分と比べる
    let h = 1e-3;
    for mu in mus() {
        let basis = Basis::new(mu, A, A);
        for (x, y) in points() {
            let v = |x: f64, y: f64| basis.values(x, y);
            let (c, xp, xm, yp, ym) = (v(x, y), v(x + h, y), v(x - h, y), v(x, y + h), v(x, y - h));
            for k in 0..8 {
                let laplacian = (xp[k] + xm[k] + yp[k] + ym[k] - 4.0 * c[k]) / (h * h);
                let scale = c[k].abs().max(1.0) * mu.abs().max(1e-3);
                assert!(
                    (laplacian - mu * c[k]).abs() < 1e-3 * scale,
                    "μ {mu} function {k} at ({x}, {y}): {laplacian} vs {}",
                    mu * c[k]
                );
            }
        }
    }
}

#[test]
fn gradients_match_finite_differences() {
    let h = 1e-5;
    for mu in mus() {
        let basis = Basis::new(mu, A, A);
        for (x, y) in points() {
            let g = basis.gradients(x, y);
            let (xp, xm) = (basis.values(x + h, y), basis.values(x - h, y));
            let (yp, ym) = (basis.values(x, y + h), basis.values(x, y - h));
            for k in 0..8 {
                for (got, want) in [
                    (g[k][0], (xp[k] - xm[k]) / (2.0 * h)),
                    (g[k][1], (yp[k] - ym[k]) / (2.0 * h)),
                ] {
                    assert!(
                        (got - want).abs() < 1e-6 * want.abs().max(1.0),
                        "μ {mu} function {k}: {got} vs {want}"
                    );
                }
            }
        }
    }
}

#[test]
fn functions_are_continuous_in_mu() {
    // 級数と閉じた式は同じ関数で、μ → 0 で多項式に近づく（8 つの関数が独立なまま）
    for (x, y) in points() {
        let at = |mu: f64| Basis::new(mu, A, A).values(x, y);
        let zero = at(0.0);
        for mu in [1e-14, -1e-14] {
            for (k, (v, z)) in at(mu).iter().zip(&zero).enumerate() {
                assert!((v - z).abs() < 1e-9 * z.abs().max(1.0), "{k}");
            }
        }
        let r2 = 2.0 * A * A;
        for edge in [1.0 / r2, -1.0 / r2] {
            let (below, above) = (at(edge * (1.0 - 1e-12)), at(edge * (1.0 + 1e-12)));
            for k in 0..8 {
                assert!(
                    (below[k] - above[k]).abs() < 1e-10 * below[k].abs().max(1.0),
                    "function {k} at the series edge: {} vs {}",
                    below[k],
                    above[k]
                );
            }
        }
    }
}

#[test]
fn node_conditions_determine_the_coefficients() {
    // 4 つの角の値、ノード平均、面の平均の外向きの法線微分の 3 つの組み合わせ（x の両面の差、y の両面の差、
    // x と y の面の和の差）の 8 × 8 の行列が、どの μ でも正則。面の微分の重み付きの和 ∮∂ψ/∂n = μ∫ψ は平均と
    // 結びつくので条件に入れない（μ = 0 では、4次式 (−x⁴/2 + 3x²y² − y⁴/2)/24 の角の値と面の微分が定数と同じ
    // 形になり、平均を使わないと区別できない）
    for (a, b) in [(A, A), (A, 0.6 * A)] {
        for mu in mus() {
            let basis = Basis::new(mu, a, b);
            let mut rows = Vec::new();
            for (x, y) in [(-a, -b), (a, -b), (a, b), (-a, b)] {
                rows.push(basis.values(x, y).to_vec());
            }
            rows.push(basis.averages(-a, a, -b, b, &rule()).to_vec());
            let f: Vec<[f64; 8]> = (0..4)
                .map(|face| basis.face_derivatives(face, &rule()))
                .collect();
            rows.push((0..8).map(|k| f[1][k] - f[0][k]).collect());
            rows.push((0..8).map(|k| f[3][k] - f[2][k]).collect());
            rows.push(
                (0..8)
                    .map(|k| b * (f[0][k] + f[1][k]) - a * (f[2][k] + f[3][k]))
                    .collect(),
            );
            let flat: Vec<f64> = rows.concat();
            // 列ごとに大きさをそろえてから解く
            let mut scaled = flat.clone();
            for k in 0..8 {
                let norm = (0..8).map(|r| flat[r * 8 + k].abs()).fold(0.0, f64::max);
                assert!(norm > 0.0, "μ {mu} function {k} vanishes");
                for r in 0..8 {
                    scaled[r * 8 + k] /= norm;
                }
            }
            let mut x = vec![1.0; 8];
            let lu = crate::lu::Lu::new(scaled, 8).unwrap_or_else(|| panic!("μ {mu} singular"));
            lu.solve(&mut x);
            assert!(
                x.iter().all(|v| v.is_finite() && v.abs() < 1e6),
                "μ {mu}: {x:?}"
            );
        }
    }
}

#[test]
fn averages_match_a_fine_midpoint_sum() {
    for mu in [0.0, -0.004, 0.3] {
        let basis = Basis::new(mu, A, 0.6 * A);
        // ノードの一部の長方形（ピンのセルの大きさ）とノード全体
        for (x0, x1, y0, y1) in [(2.0, 3.26, -5.0, -3.74), (-A, A, -0.6 * A, 0.6 * A)] {
            // 中点則（誤差は刻みの2乗に比例）を点の数 n と 2n で求め、Richardson 外挿で 4 次にする
            let midpoint = |n: usize| {
                let mut sum = [0.0; 9];
                for i in 0..n {
                    for j in 0..n {
                        let x = x0 + (x1 - x0) * (i as f64 + 0.5) / n as f64;
                        let y = y0 + (y1 - y0) * (j as f64 + 0.5) / n as f64;
                        let v = basis.values(x, y);
                        for k in 0..8 {
                            sum[k] += v[k] / (n * n) as f64;
                        }
                        sum[8] += basis.particular(x, y) / (n * n) as f64;
                    }
                }
                sum
            };
            let (coarse, fine) = (midpoint(200), midpoint(400));
            let extrapolated: Vec<f64> =
                (0..9).map(|k| (4.0 * fine[k] - coarse[k]) / 3.0).collect();
            let (sum, particular) = (&extrapolated[..8], extrapolated[8]);
            let got = basis.averages(x0, x1, y0, y1, &rule());
            for k in 0..8 {
                assert!(
                    (got[k] - sum[k]).abs() < 1e-5 * sum[k].abs().max(1e-2),
                    "μ {mu} function {k}: {} vs {}",
                    got[k],
                    sum[k]
                );
            }
            let p = basis.particular_average(x0, x1, y0, y1, &rule());
            assert!(
                (p - particular).abs() < 1e-5 * particular.abs().max(1e-2),
                "{p} vs {particular}"
            );
        }
    }
}

#[test]
fn face_derivatives_average_the_outward_normal_derivative() {
    // 面（-x、+x、-y、+y）の平均を、細かい台形則と比べる
    for mu in [0.0, -0.004, 0.3] {
        let basis = Basis::new(mu, A, 0.6 * A);
        let n = 20000;
        for face in 0..4usize {
            let (axis, sign) = (face / 2, if face.is_multiple_of(2) { -1.0 } else { 1.0 });
            let along = if axis == 0 { 0.6 * A } else { A };
            let mut sum = [0.0; 8];
            for i in 0..=n {
                let t = -along + 2.0 * along * i as f64 / n as f64;
                let (x, y) = if axis == 0 {
                    (sign * A, t)
                } else {
                    (t, sign * 0.6 * A)
                };
                let w = if i == 0 || i == n { 0.5 } else { 1.0 };
                for (s, g) in sum.iter_mut().zip(basis.gradients(x, y)) {
                    *s += w * sign * g[axis] / n as f64;
                }
            }
            let got = basis.face_derivatives(face, &rule());
            for k in 0..8 {
                assert!(
                    (got[k] - sum[k]).abs() < 1e-6 * sum[k].abs().max(1e-3),
                    "μ {mu} face {face} function {k}: {} vs {}",
                    got[k],
                    sum[k]
                );
            }
        }
    }
}

#[test]
fn particular_solution_has_unit_source() {
    // ∇²P - μP = 1（軸方向の漏れなど、ノードの中で一様な線源の特解）
    let h = 1e-3;
    for mu in mus() {
        let basis = Basis::new(mu, A, A);
        for (x, y) in points() {
            let p = |x: f64, y: f64| basis.particular(x, y);
            let laplacian =
                (p(x + h, y) + p(x - h, y) + p(x, y + h) + p(x, y - h) - 4.0 * p(x, y)) / (h * h);
            let residual = laplacian - mu * p(x, y) - 1.0;
            assert!(
                residual.abs() < 1e-3 * p(x, y).abs().max(1.0) * mu.abs().max(1e-3),
                "μ {mu}: {residual}"
            );
        }
    }
}
