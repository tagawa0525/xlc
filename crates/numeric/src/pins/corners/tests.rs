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
