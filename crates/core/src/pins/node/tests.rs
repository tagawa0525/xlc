use super::*;
use crate::pins::basis::FACE_POINTS;
use xlc_format::fewgroup::FewGroup;

fn rule() -> Vec<(f64, f64)> {
    xlc_numeric::quadrature::gauss_legendre(FACE_POINTS).unwrap()
}

/// 2群の燃料（χ = (1, 0)、上方散乱なし）
fn fuel() -> FewGroup {
    let (d1, d2) = (1.4, 0.4);
    let (a1, a2) = (0.010, 0.080);
    let s12 = 0.017;
    let (nf1, nf2) = (0.006, 0.11);
    let total = [1.0 / (3.0 * d1), 1.0 / (3.0 * d2)];
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: 1.0,
        total: total.to_vec(),
        absorption: vec![a1, a2],
        nu_fission: vec![nf1, nf2],
        fission: vec![nf1 / 2.4, nf2 / 2.4],
        chi: vec![1.0, 0.0],
        diffusion: vec![d1, d2],
        scatter: vec![vec![total[0] - a1 - s12, s12], vec![0.0, total[1] - a2]],
        adf: vec![vec![1.0, 1.0]; 4],
        vacuum_ratio: vec![None; 4],
        normalized_adf: vec![false; 4],
        pins: None,
    }
}

/// 拡散方程式 −D∇²φ + Aφ = q（q は一様）の厳密な解から、ノードの条件と、その解の値を求める関数を作る
struct Exact<'a> {
    flux: Box<dyn Fn(f64, f64) -> [f64; 2] + 'a>,
    gradient: Box<dyn Fn(f64, f64) -> [[f64; 2]; 2] + 'a>,
}

/// 厳密な解の、ノードの条件（角、平均、面の正味の中性子流、軸方向の漏れ）
fn conditions_of(exact: &Exact, m: &FewGroup, half: [f64; 2], axial: [f64; 2]) -> Conditions {
    let [a, b] = half;
    let rule = rule();
    let corners = [(-a, -b), (a, -b), (a, b), (-a, b)].map(|(x, y)| (exact.flux)(x, y));
    let mut average = [0.0; 2];
    for &(s, ws) in &rule {
        for &(t, wt) in &rule {
            let f = (exact.flux)(a * s, b * t);
            for g in 0..2 {
                average[g] += 0.25 * ws * wt * f[g];
            }
        }
    }
    // 面の正味の中性子流 J = −D ∂φ/∂(座標)（座標の正の向きを正）
    let mut currents = [[0.0; 2]; 4];
    for (face, current) in currents.iter_mut().enumerate() {
        let (axis, sign) = (face / 2, if face.is_multiple_of(2) { -1.0 } else { 1.0 });
        for &(t, w) in &rule {
            let (x, y) = if axis == 0 {
                (sign * a, t * b)
            } else {
                (t * a, sign * b)
            };
            let grad = (exact.gradient)(x, y);
            for (g, j) in current.iter_mut().enumerate() {
                *j += 0.5 * w * -m.diffusion[g] * grad[g][axis];
            }
        }
    }
    Conditions {
        half,
        average,
        corners,
        currents,
        axial_leakage: axial,
    }
}

#[test]
fn flat_conditions_give_a_flat_flux() {
    let m = fuel();
    // 無限の媒質の釣り合い: 熱群は φ2 Σa2 = φ1 Σs12
    let phi = [1.3, 1.3 * m.scatter[0][1] / m.absorption[1]];
    let conditions = Conditions {
        half: [10.75, 10.75],
        average: phi,
        corners: [phi; 4],
        currents: [[0.0; 2]; 4],
        axial_leakage: [0.0; 2],
    };
    // k は k∞（無限の媒質で平らな解が釣り合う）
    let k = (m.nu_fission[0] + m.nu_fission[1] * m.scatter[0][1] / m.absorption[1])
        / (m.absorption[0] + m.scatter[0][1]);
    let model = NodeModel::new(&m, k, conditions.half, &rule()).unwrap();
    let node = model.solve(&conditions);
    for (x, y) in [(0.0, 0.0), (5.0, -3.0), (10.75, 10.75)] {
        let f = node.flux(x, y);
        for g in 0..2 {
            assert!((f[g] - phi[g]).abs() < 1e-10 * phi[g], "{f:?}");
        }
    }
}

#[test]
fn solutions_in_the_span_are_reproduced() {
    // 2つのモードの平面波の和と、軸方向の漏れの一様な線源の特解からなる厳密な解。ノードの中の値、勾配、
    // 4 つの面の正味の中性子流（条件には 3 つの組み合わせしか入れない）を再現する
    let m = fuel();
    let k = 1.02;
    let half = [10.75, 8.0];
    let axial = [2e-4, 5e-5];
    // M = D⁻¹A の固有値と固有ベクトル（テストの中で独立に求める）
    let a = |g: usize, h: usize| {
        let d = if g == h { m.total[g] } else { 0.0 };
        d - m.scatter[h][g] - m.chi[g] * m.nu_fission[h] / k
    };
    let mm = [
        [a(0, 0) / m.diffusion[0], a(0, 1) / m.diffusion[0]],
        [a(1, 0) / m.diffusion[1], a(1, 1) / m.diffusion[1]],
    ];
    let (tr, det) = (
        mm[0][0] + mm[1][1],
        mm[0][0] * mm[1][1] - mm[0][1] * mm[1][0],
    );
    let disc = (tr * tr / 4.0 - det).sqrt();
    let mus = [tr / 2.0 - disc, tr / 2.0 + disc];
    assert!(mus[0] < 0.0 && mus[1] > 0.0, "{mus:?}");
    let v = mus.map(|mu| [mm[0][1], mu - mm[0][0]]);
    // 一様な線源 q = −(軸方向の漏れ) の定数の特解 φ_p = A⁻¹ q
    let q = [-axial[0], -axial[1]];
    let det_a = a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
    let particular = [
        (a(1, 1) * q[0] - a(0, 1) * q[1]) / det_a,
        (a(0, 0) * q[1] - a(1, 0) * q[0]) / det_a,
    ];
    let (k1, k2) = (mus[0].abs().sqrt(), mus[1].sqrt());
    let r = std::f64::consts::FRAC_1_SQRT_2;
    let (v0, v1) = (v[0], v[1]);
    // モード 0（μ < 0）: cos(k1 x) + 0.3 sin(k1 u)、モード 1（μ > 0）: 0.2 cosh(k2 y) + 0.05 sinh(k2 v)
    let psi = move |x: f64, y: f64| {
        let (u, w) = (r * (x + y), r * (x - y));
        [
            (k1 * x).cos() + 0.3 * (k1 * u).sin(),
            0.2 * (k2 * y).cosh() + 0.05 * (k2 * w).sinh(),
        ]
    };
    let dpsi = move |x: f64, y: f64| {
        let (u, w) = (r * (x + y), r * (x - y));
        let c0 = 0.3 * k1 * (k1 * u).cos() * r;
        let c1 = 0.05 * k2 * (k2 * w).cosh() * r;
        [
            [-k1 * (k1 * x).sin() + c0, c0],
            [c1, 0.2 * k2 * (k2 * y).sinh() - c1],
        ]
    };
    let exact = Exact {
        flux: Box::new(move |x, y| {
            let p = psi(x, y);
            std::array::from_fn(|g| particular[g] + v0[g] * p[0] + v1[g] * p[1])
        }),
        gradient: Box::new(move |x, y| {
            let d = dpsi(x, y);
            std::array::from_fn(|g| {
                std::array::from_fn(|axis| v0[g] * d[0][axis] + v1[g] * d[1][axis])
            })
        }),
    };
    let conditions = conditions_of(&exact, &m, half, axial);
    let model = NodeModel::new(&m, k, conditions.half, &rule()).unwrap();
    let node = model.solve(&conditions);
    for (x, y) in [(0.0, 0.0), (3.1, -7.2), (-10.0, 7.9), (10.75, 8.0)] {
        let (got, want) = (node.flux(x, y), (exact.flux)(x, y));
        let (gg, gw) = (node.gradient(x, y), (exact.gradient)(x, y));
        for g in 0..2 {
            assert!(
                (got[g] - want[g]).abs() < 1e-9 * want[g].abs().max(1e-3),
                "({x}, {y}) group {g}: {} vs {}",
                got[g],
                want[g]
            );
            for axis in 0..2 {
                assert!(
                    (gg[g][axis] - gw[g][axis]).abs() < 1e-9,
                    "gradient ({x}, {y}) group {g}"
                );
            }
        }
    }
    // 4 つの面の正味の中性子流
    let reproduced = conditions_of(
        &Exact {
            flux: Box::new(move |x, y| node.flux(x, y)),
            gradient: Box::new(move |x, y| node.gradient(x, y)),
        },
        &m,
        half,
        axial,
    );
    for face in 0..4usize {
        for g in 0..2 {
            let (got, want) = (reproduced.currents[face][g], conditions.currents[face][g]);
            assert!(
                (got - want).abs() < 1e-9 * want.abs().max(1e-6),
                "face {face} group {g}: {got} vs {want}"
            );
        }
    }
}

#[test]
fn invalid_nodes_are_errors() {
    let m = fuel();
    let conditions = Conditions {
        half: [10.75, 10.75],
        average: [1.0, 0.3],
        corners: [[1.0, 0.3]; 4],
        currents: [[0.0; 2]; 4],
        axial_leakage: [0.0; 2],
    };
    // 2群でない
    let mut three = m.clone();
    three.total.push(1.0);
    assert!(NodeModel::new(&three, 1.0, conditions.half, &rule()).is_err());
    // 固有値が複素数。非対角の項が負でない物理的な定数では起きない（M の非対角の積が負でない）ので、負の散乱で作る
    let mut complex = m.clone();
    complex.scatter[1][0] = -0.5;
    assert!(NodeModel::new(&complex, 1.0, conditions.half, &rule()).is_err());
}
