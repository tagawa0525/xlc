use super::*;
use crate::nem::Inflow;

/// 拡散係数、全断面積、散乱 [出る群][入る群]
type Constants = (Vec<f64>, Vec<f64>, Vec<Vec<f64>>);

/// 反射体（水と鋼）に近い2群の定数（上方散乱あり）
fn reflector() -> Constants {
    (
        vec![1.25, 0.30],
        vec![0.60, 1.90],
        vec![vec![0.565, 0.032], vec![0.0015, 1.885]],
    )
}

/// 細かいメッシュ（セル中心の差分、`cells` 個）の参照解: -D φ'' + A φ = -L(x)
/// （L = L̄ + l1 ξ + l2 (3ξ² - 1/4)、ξ = x/h - 1/2）。右の面は入射 `right`、左の面は入射 `left`
/// （`None` なら反射）。[群] の（ノード平均、右の面の出ていく部分中性子流、左の面の出ていく部分中性子流）を返す。
/// `faces` には [群] の（右の面の中性子束、右の外向きの正味の中性子流、左の面の中性子束、左の外向きの正味の中性子流）を書く
fn finite_difference(
    width: f64,
    cells: usize,
    transverse: &[(f64, f64, f64)],
    right: &[f64],
    left: Option<&[f64]>,
    faces: &mut Vec<[f64; 4]>,
) -> Vec<(f64, f64, f64)> {
    finite_difference_of(&reflector(), width, cells, transverse, right, left, faces)
}

/// 定数 `constants`（拡散係数、全断面積、散乱）の `finite_difference`
fn finite_difference_of(
    constants: &Constants,
    width: f64,
    cells: usize,
    transverse: &[(f64, f64, f64)],
    right: &[f64],
    left: Option<&[f64]>,
    faces: &mut Vec<[f64; 4]>,
) -> Vec<(f64, f64, f64)> {
    let (diffusion, total, scatter) = constants.clone();
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
            let right_phi = face(d, flux[g][cells - 1], right[g]);
            let (left_out, left_phi) = match left {
                Some(left) => (out(flux[g][0], left[g]), face(d, flux[g][0], left[g])),
                // 反射: 正味の中性子流 0、面の中性子束は外挿（2次の精度）
                None => {
                    let phi = 1.5 * flux[g][0] - 0.5 * flux[g][1];
                    (phi / 4.0, phi)
                }
            };
            faces.push([
                right_phi,
                2.0 * right_out - right_phi / 2.0,
                left_phi,
                2.0 * left_out - left_phi / 2.0,
            ]);
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
    let reference = finite_difference(width, 20_000, &transverse, &incoming, None, &mut Vec::new());
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
        &mut Vec::new(),
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

#[test]
fn one_node_problems_of_non_fissile_nodes_use_anm() {
    // ADF の1ノードの問題（ADR-0103）は、核分裂のないノードを ANM で解く（ADR-0102）。参照の 1D の解の
    // ノード平均と面の正味の中性子流を与えると、均質な面の中性子束は参照の面の中性子束になる（ADF は 1）。
    // 燃料側から入り真空側へ出る 21.5 cm の反射体でも正
    use crate::nem::{Edge, OneNode, surface_fluxes};
    let (diffusion, total, scatter) = reflector();
    let width = 21.5;
    let mut faces = Vec::new();
    let reference = finite_difference(
        width,
        20_000,
        &[(0.0, 0.0, 0.0); 2],
        &[0.0; 2],
        Some(&[0.05, 0.004]),
        &mut faces,
    );
    let averages: Vec<f64> = reference.iter().map(|r| r.0).collect();
    // 面（-x, +x, -y, +y）の外向きの正味の中性子流。左（-x）が入射の面
    let minus: Vec<f64> = faces.iter().map(|f| f[3]).collect();
    let plus: Vec<f64> = faces.iter().map(|f| f[1]).collect();
    let zero = [0.0; 2];
    let surface = surface_fluxes(&OneNode {
        diffusion: &diffusion,
        total: &total,
        scatter: &scatter,
        nu_fission: &[0.0; 2],
        chi: &[0.0; 2],
        k: 1.0,
        flux: &averages,
        widths: [width, width],
        currents: [&minus, &plus, &zero, &zero],
        edges: [
            Edge::Reflective,
            Edge::Vacuum,
            Edge::Reflective,
            Edge::Reflective,
        ],
    })
    .unwrap();
    for g in 0..2 {
        let [right, _, left, _] = faces[g];
        assert!(
            (surface[0][g] - left).abs() < 2e-5 * left,
            "{surface:?} {faces:?}"
        );
        assert!(
            (surface[1][g] - right).abs() < 2e-5 * right,
            "{surface:?} {faces:?}"
        );
        assert!(surface[1][g] > 0.0);
    }
}

/// 定数 `constants` の幅 `width` のノード（右から入射、左は真空、横方向の漏れなし）を、ANM と細かいメッシュの
/// 差分で解き、ANM の出ていく部分中性子流の、差分の解に対する最大の相対差（値か入射の大きいほうに対する）を返す
fn against_differences(constants: &Constants, width: f64, cells: usize) -> f64 {
    let (diffusion, total, scatter) = constants;
    let incoming = [0.05, 0.004];
    let reference = finite_difference_of(
        constants,
        width,
        cells,
        &[(0.0, 0.0, 0.0); 2],
        &incoming,
        Some(&[0.0; 2]),
        &mut Vec::new(),
    );
    let averages: Vec<f64> = reference.iter().map(|r| r.0).collect();
    let inflow: Vec<[Inflow; 2]> = (0..2)
        .map(|g| {
            [
                Inflow {
                    beta: 0.0,
                    external: incoming[g],
                },
                Inflow::default(),
            ]
        })
        .collect();
    let node = AnmNode {
        diffusion,
        total,
        scatter,
        width,
        transverse: &[(0.0, 0.0); 2],
        inflow: &inflow,
        averages: &averages,
    };
    let mut worst: f64 = 0.0;
    for g in 0..2 {
        let faces = outflow(&node, g).unwrap().at(averages[g]);
        let (_, right, left) = reference[g];
        // 厚いノードの真空側は入射より桁違いに小さいので、入射の大きさで割る
        worst = worst
            .max((faces.out[0] - right).abs() / right.max(incoming[g]))
            .max((faces.out[1] - left).abs() / left.max(incoming[g]));
    }
    worst
}

#[test]
fn thin_and_thick_nodes_match_analytic_and_fine_differences() {
    // 群が結合しない2群では、各群が φ = c1 cosh κx + c2 sinh κx（x ∈ [0, h]、左は真空、右から入射 i）と
    // 解析的に解ける。薄いノード（κh/2 ≈ 10⁻⁴）から厚いノード（κh/2 ≈ 80）まで、丸め誤差の水準で一致する
    let constants: Constants = (
        vec![1.25, 0.30],
        vec![0.60, 1.90],
        vec![vec![0.57, 0.0], vec![0.0, 1.88]],
    );
    let (diffusion, total, scatter) = &constants;
    let incoming = [0.05, 0.004];
    for width in [1.0e-3, 1.0, 200.0] {
        // [群] の（ノード平均、右の出ていく部分中性子流、左の出ていく部分中性子流）
        let exact: Vec<(f64, f64, f64)> = (0..2)
            .map(|g| {
                let d = diffusion[g];
                let kappa = ((total[g] - scatter[g][g]) / d).sqrt();
                let kh = kappa * width;
                let (c, sh) = (kh.cosh(), kh.sinh());
                // 左: 入射 φ/4 + J/2 = 0、J = −Dφ' から c1 = 2Dκ c2。右: φ/4 − J/2 = i
                let per_c2 = 2.0 * d * kappa;
                let right_per_c2 = (per_c2 * c + sh) / 4.0 + d * kappa * (per_c2 * sh + c) / 2.0;
                let c2 = incoming[g] / right_per_c2;
                let c1 = per_c2 * c2;
                let half = (0.5 * kh).sinh();
                let average = (c1 * sh + c2 * 2.0 * half * half) / kh;
                let phi_right = c1 * c + c2 * sh;
                let derivative_right = kappa * (c1 * sh + c2 * c);
                // 右の出ていく部分中性子流 φ/4 + J/2 = φ/4 − Dφ'/2、左は φ/4 − J/2
                let right = phi_right / 4.0 - d * derivative_right / 2.0;
                let left = c1 / 4.0 + d * kappa * c2 / 2.0;
                (average, right, left)
            })
            .collect();
        let averages: Vec<f64> = exact.iter().map(|e| e.0).collect();
        let inflow: Vec<[Inflow; 2]> = (0..2)
            .map(|g| {
                [
                    Inflow {
                        beta: 0.0,
                        external: incoming[g],
                    },
                    Inflow::default(),
                ]
            })
            .collect();
        let node = AnmNode {
            diffusion,
            total,
            scatter,
            width,
            transverse: &[(0.0, 0.0); 2],
            inflow: &inflow,
            averages: &averages,
        };
        for g in 0..2 {
            let faces = outflow(&node, g).unwrap().at(averages[g]);
            let (_, right, left) = exact[g];
            assert!(
                (faces.out[0] - right).abs() < 1e-12 * right,
                "{width} group {g}: {} vs {right}",
                faces.out[0]
            );
            // 厚いノードの真空側は入射より桁違いに小さいので、入射の大きさに対する丸め誤差で比べる
            assert!(
                (faces.out[1] - left).abs() < 1e-12 * left.max(incoming[g]),
                "{width} group {g}: {} vs {left}",
                faces.out[1]
            );
        }
    }
    // 群が結合する厚いノード（κh/2 ≈ 50）も、細かいメッシュの差分に一致する
    let thick = against_differences(&reflector(), 200.0, 40_000);
    assert!(thick < 1e-5, "{thick}");
}

#[test]
fn even_face_value_keeps_precision_for_small_arguments() {
    // 1 − tanh(t)/t は、t が小さいと直接の式では桁が落ちる。級数 t²/3 − 2t⁴/15 + 17t⁶/315 − 62t⁸/2835 と比べる
    for t in [1.0e-6_f64, 1.0e-4, 1.0e-2, 0.04] {
        let t2 = t * t;
        let series = t2 / 3.0 - 2.0 * t2 * t2 / 15.0 + 17.0 * t2.powi(3) / 315.0
            - 62.0 * t2.powi(4) / 2835.0;
        let value = even_face_value(t);
        assert!(
            (value - series).abs() < 1e-14 * series,
            "{t}: {value} vs {series}"
        );
    }
    // 大きい t では直接の式
    for t in [0.1_f64, 1.0, 50.0, 1.0e4] {
        let direct = 1.0 - t.tanh() / t;
        assert!((even_face_value(t) - direct).abs() < 1e-15, "{t}");
    }
}

#[test]
fn uncoupled_equal_modes_are_solved() {
    // 群の間の散乱がなく D⁻¹A の対角が等しい（固有値が重なる）ときも、固有ベクトルは単位ベクトルで解ける
    let constants: Constants = (
        vec![1.0, 0.5],
        vec![0.4, 0.3],
        vec![vec![0.36, 0.0], vec![0.0, 0.28]],
    );
    let error = against_differences(&constants, 20.0, 20_000);
    assert!(error < 2e-6, "{error}");
}

#[test]
fn ill_conditioned_modes_are_errors() {
    // 群が結合して D⁻¹A の固有値が重なる（または非常に近い）と、2つのモードの固有ベクトルがそろい、面の
    // 条件の式が悪条件になる。精度を黙って落とさず、エラーにする
    let inflow = vec![[Inflow::default(); 2]; 2];
    for gap in [0.0, 1.0e-14] {
        let diffusion = [1.0, 1.0];
        let total = [0.5, 0.48 + 0.02 + gap];
        // 除去は 0.12 と 0.12 + gap、下方散乱 0.02 で結合
        let scatter = [vec![0.38, 0.02], vec![0.0, 0.38]];
        let node = AnmNode {
            diffusion: &diffusion,
            total: &total,
            scatter: &scatter,
            width: 20.0,
            transverse: &[(0.0, 0.0); 2],
            inflow: &inflow,
            averages: &[1.0, 1.0],
        };
        let result = outflow(&node, 0);
        assert!(
            matches!(result, Err(AnmError::IllConditioned { .. })),
            "{gap}: {result:?}"
        );
    }
}

#[test]
fn non_finite_inputs_are_errors() {
    let (diffusion, total, scatter) = reflector();
    let inflow = vec![[Inflow::default(); 2]; 2];
    let base = AnmNode {
        diffusion: &diffusion,
        total: &total,
        scatter: &scatter,
        width: 20.0,
        transverse: &[(0.0, 0.0); 2],
        inflow: &inflow,
        averages: &[1.0, 1.0],
    };
    let nan_inflow = vec![
        [
            Inflow {
                beta: f64::NAN,
                external: 0.0,
            },
            Inflow::default(),
        ],
        [Inflow::default(); 2],
    ];
    let infinite_inflow = vec![
        [Inflow::default(); 2],
        [
            Inflow::default(),
            Inflow {
                beta: 0.0,
                external: f64::INFINITY,
            },
        ],
    ];
    for node in [
        AnmNode {
            transverse: &[(f64::NAN, 0.0), (0.0, 0.0)],
            ..base
        },
        AnmNode {
            inflow: &nan_inflow,
            ..base
        },
        AnmNode {
            inflow: &infinite_inflow,
            ..base
        },
        AnmNode {
            averages: &[1.0, f64::NAN],
            ..base
        },
        AnmNode {
            total: &[0.6, f64::INFINITY],
            ..base
        },
    ] {
        assert!(
            matches!(outflow(&node, 0), Err(AnmError::Invalid(_))),
            "{node:?}"
        );
    }
}
