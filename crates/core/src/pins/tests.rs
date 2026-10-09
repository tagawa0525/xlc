use super::*;
use crate::nem::{Boundary, Mesh, Options, Problem, solve};
use std::f64::consts::PI;
use xlc_format::fewgroup::{FewGroup, Pins};

/// 2群の燃料（χ = (1, 0)、上方散乱なし）。形状関数はセル `cells` × `cells` で 1、角の不連続因子は 1
fn fuel(width: f64, cells: usize) -> FewGroup {
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
        pins: Some(Pins {
            dx: vec![width / cells as f64; cells],
            dy: vec![width / cells as f64; cells],
            form: vec![vec![vec![1.0; cells]; cells]; 2],
            corners: vec![vec![1.0, 1.0]; 4],
        }),
    }
}

fn square(width: f64, nodes: usize, cells: usize, boundary: Boundary) -> Problem {
    let h = width / nodes as f64;
    let mut boundaries = [boundary; 6];
    boundaries[4] = Boundary::Reflective;
    boundaries[5] = Boundary::Reflective;
    Problem {
        mesh: Mesh {
            dx: vec![h; nodes],
            dy: vec![h; nodes],
            dz: vec![30.0],
        },
        materials: vec![fuel(h, cells)],
        node_material: vec![0; nodes * nodes],
        boundaries,
    }
}

#[test]
fn identical_reflective_nodes_give_flat_pins() {
    // 鏡面反射の中の同じノードは平らな解で、セルの出力はどれも Σ_g Σf,g φ_g
    let p = square(43.0, 2, 4, Boundary::Reflective);
    let s = solve(&p, &Options::default()).unwrap();
    let pins = reconstruct(&p, &s).unwrap();
    assert_eq!(pins.len(), 4);
    let m = &p.materials[0];
    for (n, map) in pins.iter().enumerate() {
        let expected: f64 = (0..2).map(|g| m.fission[g] * s.flux[n][g]).sum();
        assert_eq!((map.len(), map[0].len()), (4, 4));
        for v in map.iter().flatten() {
            // NEM の k と中性子束の収束の許容差（1e-8、1e-7）の程度
            assert!(
                (v - expected).abs() < 1e-6 * expected,
                "node {n}: {v} vs {expected}"
            );
        }
    }
}

#[test]
fn bare_square_matches_the_cosine_shape() {
    // 均質な裸の正方形（中性子束ゼロの境界）の基本モードは cos(πx/L) cos(πy/L)。ノードの中で再構成したセルの
    // 出力の形は、セルの平均の解析値に近づき、ノードを細かくすると誤差が減る（ADR-0096 の確認）
    let width = 172.0;
    let error = |nodes: usize| {
        let cells = 16 / nodes * 4;
        let p = square(width, nodes, cells, Boundary::ZeroFlux);
        let s = solve(&p, &Options::default()).unwrap();
        let pins = reconstruct(&p, &s).unwrap();
        // 炉心全体のセルの出力（行は -y から）と、セルの平均の cos × cos の解析値
        let total = nodes * cells;
        let cell = width / total as f64;
        let mut worst: f64 = 0.0;
        let (mut sum_got, mut sum_want) = (0.0, 0.0);
        let mut values = Vec::new();
        for j in 0..total {
            for i in 0..total {
                let (nj, ni) = (j / cells, i / cells);
                let got = pins[nj * nodes + ni][j % cells][i % cells];
                let mean = |k: usize| {
                    let (x0, x1) = (k as f64 * cell, (k + 1) as f64 * cell);
                    width / PI * ((PI * x0 / width).cos() - (PI * x1 / width).cos()) / cell
                };
                let want = mean(i) * mean(j);
                sum_got += got;
                sum_want += want;
                values.push((got, want));
            }
        }
        // 平均で規格化した出力の差の最大を、解析値の最大（中心）で割る
        let peak = values.iter().map(|v| v.1 / sum_want).fold(0.0, f64::max);
        for (got, want) in values {
            worst = worst.max((got / sum_got - want / sum_want).abs() / peak);
        }
        worst
    };
    // 2×2 の 86 cm のノードで 0.86%、4×4 の 43 cm で 0.076%（2026-10-09）
    let (coarse, fine) = (error(2), error(4));
    assert!(fine < 1e-3, "relative error {fine}");
    assert!(fine < coarse / 4.0, "{coarse} -> {fine}");
}

#[test]
fn nodes_without_form_functions_are_errors() {
    let mut p = square(43.0, 2, 4, Boundary::Reflective);
    let mut plain = p.materials[0].clone();
    plain.pins = None;
    p.materials.push(plain);
    p.node_material[3] = 1;
    let s = solve(&p, &Options::default()).unwrap();
    let err = reconstruct(&p, &s).unwrap_err().to_string();
    assert!(err.contains("form functions"), "{err}");
}

#[test]
fn cells_must_cover_the_node() {
    let mut p = square(43.0, 2, 4, Boundary::Reflective);
    if let Some(pins) = p.materials[0].pins.as_mut() {
        pins.dx[0] += 0.5;
    }
    let s = solve(&p, &Options::default()).unwrap();
    let err = reconstruct(&p, &s).unwrap_err().to_string();
    assert!(err.contains("width"), "{err}");
}
