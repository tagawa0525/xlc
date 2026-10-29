use super::*;
use crate::nem::{Boundary, Mesh, NodeSnapshot, Options, Problem, Solution, solve};

/// 2群の燃料（χ = (1, 0)、上方散乱なし）。`absorber` 倍の吸収を足すと、制御棒の入った定数の代わりになる
fn fuel(absorber: f64) -> FewGroup {
    let (d1, d2) = (1.4, 0.4);
    let (a1, a2) = (0.010 * absorber, 0.080 * absorber);
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
        adf: vec![vec![1.0, 1.0]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: None,
    }
}

const PITCH: f64 = 21.5;
const NODE: f64 = 20.0;
const NODES: usize = 10;

/// 集合体1体（径方向は反射）の高さ 200 cm の列を 20 cm のノードに分け、`tip` (cm) より上に制御棒（吸収を
/// `absorber` 倍）を入れる。z の両端は Marshak。`split` なら軸方向のメッシュを先端でも分けた参照（部分挿入の
/// ノードがない）、そうでなければ先端を含むノードを部分挿入のノードにした問題と、その区間
fn column(tip: f64, absorber: f64, split: bool) -> (Problem, Vec<PartialNode>) {
    let mut bounds: Vec<f64> = (0..=NODES).map(|i| NODE * i as f64).collect();
    if split && !bounds.contains(&tip) {
        bounds.push(tip);
        bounds.sort_by(f64::total_cmp);
    }
    let dz: Vec<f64> = bounds.windows(2).map(|w| w[1] - w[0]).collect();
    let mut materials = vec![fuel(1.0), fuel(absorber)];
    let mut node_material = Vec::new();
    let mut partial = Vec::new();
    for (z, w) in bounds.windows(2).enumerate() {
        if w[1] <= tip {
            node_material.push(Some(0));
        } else if w[0] >= tip {
            node_material.push(Some(1));
        } else {
            // 部分挿入のノードは自分だけの物質を持つ
            materials.push(fuel(1.0));
            node_material.push(Some(materials.len() - 1));
            partial.push(PartialNode {
                site: z,
                material: materials.len() - 1,
                segments: vec![
                    Segment {
                        constants: fuel(1.0),
                        length: tip - w[0],
                    },
                    Segment {
                        constants: fuel(absorber),
                        length: w[1] - tip,
                    },
                ],
            });
        }
    }
    let problem = Problem {
        mesh: Mesh {
            dx: vec![PITCH],
            dy: vec![PITCH],
            dz,
        },
        materials,
        reflectors: Vec::new(),
        node_material,
        boundaries: [
            Boundary::Reflective,
            Boundary::Reflective,
            Boundary::Reflective,
            Boundary::Reflective,
            Boundary::Marshak,
            Boundary::Marshak,
        ],
        outline: None,
    };
    (problem, partial)
}

/// 核分裂率を、高さ 20 cm の粗いノードごとに足したもの（先端で分けた参照の 2 つのノードは 1 つにまとめる）
fn coarse_fission(problem: &Problem, materials: &[FewGroup], solution: &Solution) -> Vec<f64> {
    let mut rates = vec![0.0; NODES];
    let mut z0 = 0.0;
    for (z, dz) in problem.mesh.dz.iter().enumerate() {
        let m = &materials[problem.node_material[z].unwrap()];
        let rate: f64 = m
            .fission
            .iter()
            .zip(&solution.flux[z])
            .map(|(f, phi)| f * phi)
            .sum();
        rates[((z0 + 0.5 * dz) / NODE) as usize] += rate * dz;
        z0 += dz;
    }
    rates
}

fn pcm(a: f64, b: f64) -> f64 {
    1.0e5 * (a - b)
}

/// 粗いノードの核分裂率の、参照に対する相対差の最大
fn max_relative(rates: &[f64], reference: &[f64]) -> f64 {
    rates
        .iter()
        .zip(reference)
        .map(|(r, q)| ((r - q) / q).abs())
        .fold(0.0, f64::max)
}

#[test]
fn identical_intervals_reproduce_the_uniform_node() {
    // 区間の定数が同じなら、混ぜた定数はその定数で、k も一様なノードの問題と同じ
    let (plain, _) = column(200.0, 1.0, false);
    let expected = solve(&plain, &Options::default()).unwrap();
    let (problem, partial) = column(128.0, 1.0, false);
    let (solution, materials, _) =
        solve_rodded(&problem, partial, Weighting::Flux, 0.5, &Options::default()).unwrap();
    // 不連続因子が 1 からずれる分だけ k も動く
    assert!(
        pcm(solution.k, expected.k).abs() < 0.1,
        "{}",
        pcm(solution.k, expected.k)
    );
    let mixed = &materials[problem.node_material[6].unwrap()];
    for g in 0..2 {
        assert!((mixed.absorption[g] - fuel(1.0).absorption[g]).abs() < 1.0e-12);
        // z の面の不連続因子は、細かいメッシュの差分と NEM の面の中性子束の差だけ 1 からずれる
        for side in [4, 5] {
            assert!((mixed.adf[side][g] - 1.0).abs() < 1.0e-3, "{:?}", mixed.adf);
        }
    }
}

#[test]
fn f_star_below_a_partial_node_keeps_the_split_mesh_k() {
    // 部分挿入のノードの下の一様なノードが +z の辺に f*（燃料の ADF で割った値。ADR-0095）を持つとき、その面は
    // f* × 接する区間の ADF でつながり、混ぜた定数の z の面の不連続因子を掛けない。先端をノードの境に置いた参照との
    // k の差は、f* のない列（flux_weighting_reproduces_the_mesh_split_at_the_tip）と同じく小さい
    let with_adfs = |(mut problem, mut partial): (Problem, Vec<PartialNode>)| {
        for m in &mut problem.materials {
            m.adf[4] = vec![0.8, 0.7];
            m.adf[5] = vec![0.8, 0.7];
        }
        for p in &mut partial {
            for s in &mut p.segments {
                s.constants.adf[4] = vec![0.8, 0.7];
                s.constants.adf[5] = vec![0.8, 0.7];
            }
        }
        let mut below = problem.materials[0].clone();
        below.adf[5] = vec![1.3, 1.5];
        below.normalized_adf[5] = true;
        problem.materials.push(below);
        problem.node_material[0] = Some(problem.materials.len() - 1);
        (problem, partial)
    };
    let (plain, _) = with_adfs(column(28.0, 1.6, true));
    let expected = solve(&plain, &Options::default()).unwrap();
    let (problem, partial) = with_adfs(column(28.0, 1.6, false));
    let (solution, _, _) =
        solve_rodded(&problem, partial, Weighting::Flux, 0.5, &Options::default()).unwrap();
    let dk = pcm(solution.k, expected.k);
    assert!(dk.abs() < 10.0, "{dk} pcm");
}

#[test]
fn volume_weighting_mixes_by_the_rodded_fraction() {
    let (problem, partial) = column(128.0, 1.6, false);
    let (_, materials, _) = solve_rodded(
        &problem,
        partial,
        Weighting::Volume,
        0.5,
        &Options::default(),
    )
    .unwrap();
    let mixed = &materials[problem.node_material[6].unwrap()];
    for g in 0..2 {
        let (out, rodded) = (fuel(1.0).absorption[g], fuel(1.6).absorption[g]);
        let expected = out + 0.6 * (rodded - out);
        assert!((mixed.absorption[g] - expected).abs() < 1.0e-14);
        assert_eq!(mixed.adf[4][g], 1.0);
    }
}

#[test]
fn flux_weighting_reproduces_the_mesh_split_at_the_tip() {
    // 先端をノードの境に置いた参照に対し、中性子束の重みと z の面の不連続因子（ADR-0127）の k と軸方向の
    // 核分裂率の差は、体積の重みより十分小さい
    for tip in [122.0, 128.0, 138.0] {
        let (reference, _) = column(tip, 1.6, true);
        let expected = solve(&reference, &Options::default()).unwrap();
        let expected_rates = coarse_fission(&reference, &reference.materials, &expected);
        let (problem, partial) = column(tip, 1.6, false);
        let run = |weighting: Weighting| {
            let (solution, materials, _) = solve_rodded(
                &problem,
                partial.clone(),
                weighting,
                0.5,
                &Options::default(),
            )
            .unwrap();
            let rates = coarse_fission(&problem, &materials, &solution);
            (
                pcm(solution.k, expected.k),
                max_relative(&rates, &expected_rates),
            )
        };
        let (dk_flux, dp_flux) = run(Weighting::Flux);
        let (dk_volume, dp_volume) = run(Weighting::Volume);
        assert!(dk_flux.abs() < 10.0, "tip {tip}: {dk_flux} pcm");
        assert!(
            dk_flux.abs() < dk_volume.abs() / 5.0,
            "tip {tip}: {dk_flux} vs {dk_volume}"
        );
        assert!(
            dp_flux < dp_volume / 5.0,
            "tip {tip}: {dp_flux} vs {dp_volume}"
        );
    }
}

#[test]
fn windows_add_one_uniform_node_on_each_side_and_merge_shared_ones() {
    let w = |partial: &[bool]| windows(partial).unwrap();
    assert_eq!(
        w(&[false, false, true, false, false]),
        vec![Window {
            nodes: 1..4,
            constrained: [true, true]
        }]
    );
    // 炉心の端の側は一様なノードを加えない
    assert_eq!(
        w(&[true, false, false]),
        vec![Window {
            nodes: 0..2,
            constrained: [false, true]
        }]
    );
    // 2 つの窓が同じ一様なノードを端に持てば 1 つにまとめる
    assert_eq!(
        w(&[false, true, false, true, false]),
        vec![Window {
            nodes: 0..5,
            constrained: [true, true]
        }]
    );
    // 部分挿入のノードが続けば 1 つの窓、間に一様なノードが 2 つあれば別の窓
    assert_eq!(
        w(&[false, true, true, false, false, true, false]),
        vec![
            Window {
                nodes: 0..4,
                constrained: [true, true]
            },
            Window {
                nodes: 4..7,
                constrained: [true, true]
            }
        ]
    );
    assert!(w(&[false, false]).is_empty());
}

#[test]
fn windows_reaching_both_core_edges_keep_the_shared_uniform_node_as_their_ends() {
    // まとめると両端が炉心の端になる 2 つの窓は、共有する一様なノードでそれぞれ拘束する
    assert_eq!(
        windows(&[true, false, true]).unwrap(),
        vec![
            Window {
                nodes: 0..2,
                constrained: [false, true]
            },
            Window {
                nodes: 1..3,
                constrained: [true, false]
            }
        ]
    );
}

#[test]
fn partial_nodes_at_both_ends_around_one_uniform_node_match_the_split_mesh() {
    // 3 ノードの列の下と上のノードが部分挿入（下は 10 cm、上は 5 cm より上に棒）。z の両端は反射（炉心の
    // Marshak の端の部分挿入のノードは #0102）
    let problem_of = |dz: Vec<f64>, node_material: Vec<Option<usize>>, materials| Problem {
        mesh: Mesh {
            dx: vec![PITCH],
            dy: vec![PITCH],
            dz,
        },
        materials,
        reflectors: Vec::new(),
        node_material,
        boundaries: [
            Boundary::Reflective,
            Boundary::Reflective,
            Boundary::Reflective,
            Boundary::Reflective,
            Boundary::Reflective,
            Boundary::Reflective,
        ],
        outline: None,
    };
    let segments = |below: f64| {
        vec![
            Segment {
                constants: fuel(1.0),
                length: below,
            },
            Segment {
                constants: fuel(1.6),
                length: NODE - below,
            },
        ]
    };
    let problem = problem_of(
        vec![NODE; 3],
        vec![Some(1), Some(0), Some(2)],
        vec![fuel(1.0), fuel(1.0), fuel(1.0)],
    );
    let partial = vec![
        PartialNode {
            site: 0,
            material: 1,
            segments: segments(10.0),
        },
        PartialNode {
            site: 2,
            material: 2,
            segments: segments(5.0),
        },
    ];
    let reference = problem_of(
        vec![10.0, 10.0, 20.0, 5.0, 15.0],
        vec![Some(0), Some(1), Some(0), Some(0), Some(1)],
        vec![fuel(1.0), fuel(1.6)],
    );
    let expected = solve(&reference, &Options::default()).unwrap();
    let (solution, _, _) =
        solve_rodded(&problem, partial, Weighting::Flux, 0.5, &Options::default()).unwrap();
    let dk = pcm(solution.k, expected.k);
    assert!(dk.abs() < 10.0, "{dk} pcm");
}

#[test]
fn windows_without_a_uniform_node_are_errors() {
    // どちらの側にも一様なノードがなければ、細かいメッシュの解の大きさが決まらない
    assert!(windows(&[true]).is_err());
    assert!(windows(&[true, true]).is_err());
}

#[test]
fn partial_nodes_without_fission_are_errors() {
    // 核分裂のない層の部分挿入（軸方向の反射体の層の制御棒の軸。#0030）はまだ扱わない
    let (mut problem, mut partial) = column(128.0, 1.6, false);
    for s in &mut partial[0].segments {
        s.constants.nu_fission = vec![0.0, 0.0];
        s.constants.fission = vec![0.0, 0.0];
    }
    let m = partial[0].material;
    problem.materials[m].nu_fission = vec![0.0, 0.0];
    problem.materials[m].fission = vec![0.0, 0.0];
    assert!(solve_rodded(&problem, partial, Weighting::Flux, 0.5, &Options::default()).is_err());
}

#[test]
fn partial_nodes_with_f_star_on_a_z_side_are_errors() {
    // 混ぜた定数の z の辺は細かいメッシュの不連続因子で置き換えるので、f*（燃料の ADF で割った値）を持てない
    let (mut problem, mut partial) = column(128.0, 1.6, false);
    for s in &mut partial[0].segments {
        s.constants.normalized_adf[4] = true;
    }
    let m = partial[0].material;
    problem.materials[m].normalized_adf[4] = true;
    assert!(solve_rodded(&problem, partial, Weighting::Flux, 0.5, &Options::default()).is_err());
}

#[test]
fn duplicated_partial_nodes_are_errors() {
    // 同じ位置の部分挿入のノードが 2 つあれば、どちらの区間を使うか分からない
    let (problem, mut partial) = column(128.0, 1.6, false);
    partial.push(partial[0].clone());
    assert!(
        solve_rodded(
            &problem,
            partial,
            Weighting::Volume,
            0.5,
            &Options::default()
        )
        .is_err()
    );
}

#[test]
fn one_node_z_problems_with_nan_do_not_converge() {
    // 群 0 の変化が NaN、群 1 が有限でも、収束したとしない（後の有限の値で NaN を消さない）
    let m = fuel(1.0);
    let snapshot = NodeSnapshot {
        flux: vec![f64::NAN, 1.0],
        transverse: vec![[0.0; 3]; 2],
    };
    let result = homogeneous_faces(
        &m,
        1.0,
        20.0,
        &[f64::NAN, 1.0],
        &[0.0, 0.0],
        &[0.0, 0.0],
        &snapshot,
    );
    assert!(result.is_err(), "{result:?}");
}

// 部分挿入のノードのミクロ燃焼（ADR-0217）

use xlc_format::table::{MicroNuclide, MicroReaction};

use crate::micro::correct;

/// 2 群の微視的断面積（barn）。`scale` 倍した値で区間ごとに変える
fn micro_nuclide(scale: f64, thermal: Option<f64>) -> MicroNuclide {
    let v = |x: [f64; 2]| x.iter().map(|x| scale * x).collect::<Vec<_>>();
    MicroNuclide {
        total: v([5.0, 60.0]),
        nu_fission: v([2.4, 100.0]),
        fission: v([1.0, 41.0]),
        scatter: vec![v([3.9, 0.1]), v([0.0, 10.0])],
        reactions: vec![MicroReaction {
            mt: 102,
            channels: vec![v([0.3, 8.0]), v([0.1, 2.0])],
        }],
        thermal_fission: thermal,
    }
}

/// 区間ごとの 2 核種の微視的断面積。棒の入った区間（`rodded`）は熱群を小さくし、収率の境界より下の割合も変える
fn segment_micro(rodded: bool) -> MicroAt {
    let (scale, thermal) = if rodded { (0.7, 0.8) } else { (1.0, 0.95) };
    MicroAt {
        densities: vec![7.0e-4, 1.0e-5],
        nuclides: vec![
            micro_nuclide(scale, Some(thermal)),
            micro_nuclide(0.5 * scale, None),
        ],
    }
}

/// 追う核種の数密度 (1/(barn cm))
const DENSITIES: [f64; 2] = [6.5e-4, 2.0e-5];

#[test]
fn mixed_micro_cross_sections_take_the_weights_of_the_mixed_constants() {
    // 区間の定数が、同じ下地に追う核種の Σᵢ Nᵢ σᵢ を足したものなら、同じ重みで混ぜた定数は、下地に Σᵢ Nᵢ（混ぜた σᵢ）を
    // 足したものと同じ（少数群定数と微視的断面積を同じ関数で混ぜる。ADR-0217）
    let zero = |m: &MicroAt| MicroAt {
        densities: vec![0.0; m.densities.len()],
        nuclides: m.nuclides.clone(),
    };
    let micro = [segment_micro(false), segment_micro(true)];
    let segments: Vec<Segment> = micro
        .iter()
        .map(|m| {
            let mut constants = fuel(1.0);
            correct(&mut constants, &zero(m), &DENSITIES);
            Segment {
                constants,
                length: 10.0,
            }
        })
        .collect();
    let weights = vec![vec![3.0, 1.0], vec![0.5, 0.2]];
    let mixed = mix(&segments, &weights).unwrap();
    let flux = [2.0, 0.5];
    let mixed_micro = mix_micro(&micro, &weights, &flux).unwrap();
    let mut expected = fuel(1.0);
    correct(&mut expected, &zero(&mixed_micro), &DENSITIES);
    let close = |a: &[f64], b: &[f64]| {
        a.iter()
            .zip(b)
            .all(|(x, y)| (x - y).abs() <= 1.0e-13 * y.abs().max(1.0e-30))
    };
    assert!(close(&mixed.total, &expected.total), "{mixed:?}");
    assert!(close(&mixed.absorption, &expected.absorption));
    assert!(close(&mixed.nu_fission, &expected.nu_fission));
    assert!(close(&mixed.fission, &expected.fission));
    for (a, b) in mixed.scatter.iter().zip(&expected.scatter) {
        assert!(close(a, b));
    }
    // 経路の断面積は群の中性子束で混ぜる
    let u235 = &mixed_micro.nuclides[0];
    for (c, channel) in u235.reactions[0].channels.iter().enumerate() {
        for g in 0..2 {
            let expected = (weights[0][g] * micro[0].nuclides[0].reactions[0].channels[c][g]
                + weights[1][g] * micro[1].nuclides[0].reactions[0].channels[c][g])
                / (weights[0][g] + weights[1][g]);
            assert!((channel[g] - expected).abs() <= 1.0e-13 * expected);
        }
    }
    // 収率の境界より下の割合は、燃焼が使う核分裂率 σ̄f,G φ̄_G のうちの区間の分で混ぜる（熱の核分裂率の和を保つ）
    let rate = |s: usize| -> f64 {
        (0..2)
            .map(|g| {
                micro[s].nuclides[0].fission[g] * weights[s][g] / (weights[0][g] + weights[1][g])
                    * flux[g]
            })
            .sum()
    };
    let rates = [rate(0), rate(1)];
    let solved_rate: f64 = (0..2).map(|g| u235.fission[g] * flux[g]).sum();
    assert!((rates[0] + rates[1] - solved_rate).abs() < 1.0e-13 * solved_rate);
    let thermal = (0.95 * rates[0] + 0.8 * rates[1]) / (rates[0] + rates[1]);
    assert!((u235.thermal_fission.unwrap() - thermal).abs() < 1.0e-14);
    assert_eq!(mixed_micro.nuclides[1].thermal_fission, None);
    // 公称の数密度は区間によらない（同じ層のテーブルの、同じ燃焼度の値）
    assert_eq!(mixed_micro.densities, micro[0].densities);
}

#[test]
fn mixing_micro_cross_sections_needs_matching_segments() {
    let micro = [segment_micro(false), segment_micro(true)];
    let flux = [1.0, 1.0];
    let even = vec![vec![1.0, 1.0]; 2];
    // 重みの区間と群の数
    assert!(mix_micro(&micro, &[vec![1.0, 1.0]], &flux).is_err());
    assert!(mix_micro(&micro, &[vec![1.0], vec![1.0]], &flux).is_err());
    // 中性子束と群の数の違う微視的断面積（群ごとに引く前に止める）
    let mut short = segment_micro(true);
    short.nuclides[1].fission.pop();
    assert!(mix_micro(&[segment_micro(false), short], &even, &flux).is_err());
    // 公称の数密度の違う区間（同じ層のテーブルの同じ燃焼度なら同じはず）
    let mut other = segment_micro(true);
    other.densities[0] *= 1.1;
    assert!(mix_micro(&[segment_micro(false), other], &even, &flux).is_err());
    // 核分裂する区間としない区間の混じった核種
    let mut fissionless = segment_micro(true);
    fissionless.nuclides[0].thermal_fission = None;
    assert!(mix_micro(&[segment_micro(false), fissionless], &even, &flux).is_err());
}

#[test]
fn rodded_solutions_return_the_weights_of_their_mixed_constants() {
    // 返す重みで区間の定数を混ぜると、解いた混ぜた定数（z の面の不連続因子のほか）になる。体積の重みは区間の長さ
    for weighting in [Weighting::Flux, Weighting::Volume] {
        let (problem, partial) = column(128.0, 1.6, false);
        let (_, materials, weights) = solve_rodded(
            &problem,
            partial.clone(),
            weighting,
            0.5,
            &Options::default(),
        )
        .unwrap();
        assert_eq!(weights.len(), 1);
        let solved = &materials[partial[0].material];
        let mixed = mix(&partial[0].segments, &weights[0]).unwrap();
        assert_eq!(mixed.total, solved.total, "{weighting:?}");
        assert_eq!(mixed.absorption, solved.absorption);
        assert_eq!(mixed.nu_fission, solved.nu_fission);
        assert_eq!(mixed.scatter, solved.scatter);
        assert_eq!(mixed.diffusion, solved.diffusion);
        assert_eq!(mixed.chi, solved.chi);
        if weighting == Weighting::Volume {
            assert_eq!(weights[0], vec![vec![8.0; 2], vec![12.0; 2]]);
        }
    }
}
