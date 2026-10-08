//! 解析的ノード法（ANM）の1ノードの計算（ADR-0102）。

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
                    // 反射: 正味の流れ 0、面の中性子束は外挿（2次の精度）
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
