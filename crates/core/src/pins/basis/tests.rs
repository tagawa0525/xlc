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
        for k in 0..8 {
            for mu in [1e-14, -1e-14] {
                let v = at(mu)[k];
                assert!((v - zero[k]).abs() < 1e-9 * zero[k].abs().max(1.0), "{k}");
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
    // 4 つの角の値と、4 つの面の平均の外向きの法線微分の 8 × 8 の行列が、どの μ でも正則
    for mu in mus() {
        let basis = Basis::new(mu, A, A);
        let mut rows = Vec::new();
        for (x, y) in [(-A, -A), (A, -A), (A, A), (-A, A)] {
            rows.push(basis.values(x, y).to_vec());
        }
        for face in 0..4 {
            rows.push(basis.face_derivatives(face).to_vec());
        }
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
        let mut b = vec![1.0; 8];
        let lu = xlc_numeric::lu::Lu::new(scaled, 8).unwrap_or_else(|| panic!("μ {mu} singular"));
        lu.solve(&mut b);
        assert!(
            b.iter().all(|v| v.is_finite() && v.abs() < 1e8),
            "μ {mu}: {b:?}"
        );
    }
}

#[test]
fn face_derivatives_average_the_outward_normal_derivative() {
    // 面（-x、+x、-y、+y）の平均を、細かい台形則と比べる
    for mu in [0.0, -0.004, 0.3] {
        let basis = Basis::new(mu, A, 0.6 * A);
        let n = 20000;
        for face in 0..4 {
            let (axis, sign) = (face / 2, if face % 2 == 0 { -1.0 } else { 1.0 });
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
            let got = basis.face_derivatives(face);
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
