use super::*;
use crate::nem::{Boundary, Mesh, Options, Problem, Solution, solve};

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
    let (solution, materials) =
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
fn volume_weighting_mixes_by_the_rodded_fraction() {
    let (problem, partial) = column(128.0, 1.6, false);
    let (_, materials) = solve_rodded(
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
            let (solution, materials) = solve_rodded(
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
