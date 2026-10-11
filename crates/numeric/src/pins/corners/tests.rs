use super::*;
use crate::pins::FACE_POINTS;
use crate::pins::node::Constants;

fn rule() -> Vec<(f64, f64)> {
    crate::quadrature::gauss_legendre(FACE_POINTS).unwrap()
}

#[test]
fn known_corner_fluxes_scale_the_corner_estimates() {
    // core の形（未知数が角の点の不均質な中性子束 Φ、ノードの均質な角の値は Φ / C）で解いた Φ を既知とし、ノードの
    // 均質な角の値を s_p × 見込み（見込みのノードの間の比は 1/C と同じ、角の点ごとの大きさは任意の t_p）として解くと、
    // 外周の真空の辺の流れ r Φ は既知の量になり、s_p t_p が Φ_p に戻る（ADR-0125）
    let total = [1.0 / (3.0 * 1.4), 1.0 / (3.0 * 0.4)];
    let scatter = vec![
        vec![total[0] - 0.010 - 0.017, 0.017],
        vec![0.0, total[1] - 0.08],
    ];
    let (nu_fission, chi, diffusion) = ([0.006, 0.11], [1.0, 0.0], [1.4, 0.4]);
    let constants = Constants {
        total: &total,
        scatter: &scatter,
        nu_fission: &nu_fission,
        chi: &chi,
        diffusion: &diffusion,
    };
    let half = [10.75, 10.75];
    let model = NodeModel::new(&constants, 1.0, half, &rule()).unwrap();
    // 2 × 1 のノード。面の正味の中性子流は共通の面で連続、+x の外周は真空
    let conditions = [
        Conditions {
            half,
            average: [1.0, 0.3],
            corners: [[0.0; 2]; 4],
            currents: [[0.0; 2], [0.01, 0.002], [0.0; 2], [0.0; 2]],
            axial_leakage: [0.0; 2],
        },
        Conditions {
            half,
            average: [0.8, 0.25],
            corners: [[0.0; 2]; 4],
            currents: [[0.01, 0.002], [0.03, 0.004], [0.0; 2], [0.0; 2]],
            axial_leakage: [0.0; 2],
        },
    ];
    let cdf = |n: usize, c: usize, g: usize| 0.9 + 0.05 * (n + c + g) as f64;
    let outflow = |n: usize| {
        if n == 1 {
            [None, Some([0.4, 0.3]), None, None]
        } else {
            [None; 4]
        }
    };
    let nodes = |weights: &dyn Fn(usize, usize, usize) -> f64| -> Vec<CornerNode> {
        (0..2)
            .map(|n| CornerNode {
                model: &model,
                conditions: &conditions[n],
                diffusion,
                weights: std::array::from_fn(|c| std::array::from_fn(|g| weights(n, c, g))),
                outflow: outflow(n),
            })
            .collect()
    };
    let zero = vec![false; 6];
    let phi = corner_balance(2, 1, &nodes(&|n, c, g| 1.0 / cdf(n, c, g)), &zero, None).unwrap();
    // ノード n の角 c の角の点（x が速い 3 × 2）
    let point = |n: usize, c: usize| [n, n + 1, n + 4, n + 3][c];
    let t = |p: usize| 1.0 + 0.1 * p as f64;
    let s = corner_balance(
        2,
        1,
        &nodes(&|n, c, g| t(point(n, c)) / cdf(n, c, g)),
        &zero,
        Some(&phi),
    )
    .unwrap();
    for p in 0..6 {
        for g in 0..2 {
            let (got, want) = (s[p][g] * t(p), phi[p][g]);
            assert!(
                (got - want).abs() < 1e-10 * want.abs(),
                "point {p} group {g}: {got} vs {want}"
            );
        }
    }
}

/// 核分裂のない 2 群の正方形（中心が原点、半幅 a）の、細かいメッシュの拡散の解。一様な高速群の線源 1、-x と -y は
/// 反射、+x と +y は J = r φ。格子点の差分（境界の格子点は半分の体積、角は 1/4）を SOR で解く。[群] の
/// （ノード平均、+x の面の平均の正味の中性子流、+y の面の平均の正味の中性子流、[角] の値）
fn fine_mesh_diffusion(
    a: f64,
    diffusion: [f64; 2],
    removal: [f64; 2],
    transfer: f64,
    ratio: f64,
    cells: usize,
) -> [(f64, f64, f64, [f64; 4]); 2] {
    let n = cells + 1;
    let h = 2.0 * a / cells as f64;
    // 格子点 (i, j) の体積と台形の重み
    let weight = |i: usize| if i == 0 || i == cells { 0.5 } else { 1.0 };
    let mut out = [(0.0, 0.0, 0.0, [0.0; 4]); 2];
    let mut previous = vec![0.0; n * n];
    for g in 0..2 {
        let d = diffusion[g];
        let source: Vec<f64> = if g == 0 {
            vec![1.0; n * n]
        } else {
            previous.iter().map(|p| transfer * p).collect()
        };
        let mut phi = vec![0.0; n * n];
        let omega = 2.0 / (1.0 + (std::f64::consts::PI / n as f64).sin());
        for _ in 0..20000 {
            let mut change: f64 = 0.0;
            for j in 0..n {
                for i in 0..n {
                    let (wx, wy) = (weight(i), weight(j));
                    // 面の結合（隣の格子点との間の辺の長さ × D / h）
                    let mut diagonal = removal[g] * wx * wy * h * h;
                    let mut sum = source[i + n * j] * wx * wy * h * h;
                    for (di, dj) in [(-1i64, 0i64), (1, 0), (0, -1), (0, 1)] {
                        let (ii, jj) = (i as i64 + di, j as i64 + dj);
                        if ii < 0 || jj < 0 || ii >= n as i64 || jj >= n as i64 {
                            continue;
                        }
                        let length = if di != 0 { wy } else { wx } * h;
                        let c = d * length / h;
                        diagonal += c;
                        sum += c * phi[ii as usize + n * jj as usize];
                    }
                    // +x と +y の真空の辺
                    if i == cells {
                        diagonal += ratio * wy * h;
                    }
                    if j == cells {
                        diagonal += ratio * wx * h;
                    }
                    let new = (1.0 - omega) * phi[i + n * j] + omega * sum / diagonal;
                    change = change.max((new - phi[i + n * j]).abs());
                    phi[i + n * j] = new;
                }
            }
            if change < 1e-13 {
                break;
            }
        }
        let total: f64 = (0..n * n)
            .map(|k| weight(k % n) * weight(k / n) * phi[k])
            .sum::<f64>()
            * h
            * h;
        let face = |at: &dyn Fn(usize) -> usize| -> f64 {
            (0..n).map(|k| weight(k) * phi[at(k)]).sum::<f64>() * h * ratio / (2.0 * a)
        };
        out[g] = (
            total / (4.0 * a * a),
            face(&|k| cells + n * k),
            face(&|k| k + n * cells),
            [phi[0], phi[cells], phi[cells + n * cells], phi[n * cells]],
        );
        previous = phi;
    }
    out
}

#[test]
fn vacuum_corners_follow_the_fine_mesh_diffusion() {
    // 2 辺が真空の、核分裂のないノードで、面の平均の量から角の点のバランスで解いた均質な角の値は正で、細かい
    // メッシュの拡散の解の角の値に近い。真空の 2 辺の角は展開が表しきれず、拡散距離（熱群 4.3 cm）に比べてノードが
    // 大きいほど誤差が大きい（#0065）
    let (diffusion, removal, transfer, ratio) = ([1.4, 0.37], [0.03, 0.02], 0.025, 0.5);
    let total = [1.0 / (3.0 * diffusion[0]), 1.0 / (3.0 * diffusion[1])];
    let scatter = vec![
        vec![total[0] - removal[0], transfer],
        vec![0.0, total[1] - removal[1]],
    ];
    let constants = Constants {
        total: &total,
        scatter: &scatter,
        nu_fission: &[0.0, 0.0],
        chi: &[1.0, 0.0],
        diffusion: &diffusion,
    };
    // （半幅、2 辺が真空の角の相対誤差の上限）
    for (a, bound) in [(1.5, 0.02), (3.0, 0.06), (6.0, 0.20)] {
        let reference = fine_mesh_diffusion(a, diffusion, removal, transfer, ratio, 60);
        let half = [a, a];
        let model = NodeModel::new(&constants, 1.0, half, &rule()).unwrap();
        let conditions = Conditions {
            half,
            average: std::array::from_fn(|g| reference[g].0),
            corners: [[0.0; 2]; 4],
            currents: [
                [0.0; 2],
                std::array::from_fn(|g| reference[g].1),
                [0.0; 2],
                std::array::from_fn(|g| reference[g].2),
            ],
            axial_leakage: [-1.0, 0.0],
        };
        let node = CornerNode {
            model: &model,
            conditions: &conditions,
            diffusion,
            weights: [[1.0; 2]; 4],
            outflow: [None, Some([ratio; 2]), None, Some([ratio; 2])],
        };
        let corners = corner_balance(1, 1, &[node], &[false; 4], None).unwrap();
        // 角の点（x が速い）と、参照の角（-x -y、+x -y、+x +y、-x +y）
        for (p, c) in [0, 1, 3, 2].into_iter().enumerate() {
            let limit = if c == 2 { bound } else { 0.01 };
            for g in 0..2 {
                let (got, want) = (corners[p][g], reference[g].3[c]);
                assert!(
                    got > 0.0 && (got / want - 1.0).abs() < limit,
                    "half width {a}, corner {c}, group {g}: {got} vs {want}"
                );
            }
        }
    }
}
