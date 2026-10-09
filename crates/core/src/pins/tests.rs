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
            halves: vec![false; 4],
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

#[test]
fn core_rows_join_nodes_from_the_top() {
    // 2 × 1 のノード、2 層（厚さ 1 と 3 cm）。ノードのセルは [行（-y から）][列]。すべて 0 の行と列（余白）を除き、
    // 層の厚さで重みを付けて足し、核分裂のあるセルの平均で 1 に規格化する。行は上から
    let p = Problem {
        mesh: Mesh {
            dx: vec![2.0, 2.0],
            dy: vec![2.0],
            dz: vec![1.0, 3.0],
        },
        materials: vec![fuel(2.0, 2)],
        node_material: vec![0; 4],
        boundaries: [Boundary::Reflective; 6],
    };
    let left = vec![vec![1.0, 0.0], vec![2.0, 0.0]];
    let right = vec![vec![3.0, 0.0], vec![0.0, 0.0]];
    // 層 0 と層 1 で同じ形、層 1 は 2 倍
    let twice = |m: &Vec<Vec<f64>>| {
        m.iter()
            .map(|r| r.iter().map(|v| 2.0 * v).collect())
            .collect()
    };
    let maps = vec![left.clone(), right.clone(), twice(&left), twice(&right)];
    let rows = core_rows(&p, &maps).unwrap();
    // 層の和は 7 倍: 左 (7, 0 / 14, 0)、右 (21, 0 / 0, 0)。列 1 と 3 はすべて 0 なので除く。上の行は y の大きい行
    let (a, b, c) = (14.0, 7.0, 21.0);
    let mean = (a + b + c) / 3.0;
    assert_eq!(rows, vec![vec![a / mean, 0.0], vec![b / mean, c / mean]]);
    // ノードのセルの行の数が、同じ行のノードで違えばエラー
    let mut short = maps.clone();
    short[1].pop();
    assert!(core_rows(&p, &short).is_err());
}

#[test]
fn vacuum_corners_lose_r_times_the_corner_flux() {
    // 平板（x は Marshak、y と z は反射）の左下の角の点の箱は、ノード 0 だけを含み、箱の辺から出ていく流れ
    // J_x + J_y（J = −D∇φ）と、真空の辺から出ていく r Φ（r = 1/2、Φ は不均質な角の中性子束）の和が 0。面の ADF（1）と
    // 違う CDF（1.3）でも、真空の辺の流れは ADF でなく Φ で決まる
    let mut p = square(86.0, 4, 2, Boundary::Marshak);
    p.mesh.dy = vec![86.0 / 4.0];
    p.node_material = vec![0; 4];
    p.boundaries[2] = Boundary::Reflective;
    p.boundaries[3] = Boundary::Reflective;
    if let Some(pins) = p.materials[0].pins.as_mut() {
        pins.corners = vec![vec![1.3, 1.3]; 4];
    }
    let s = solve(&p, &Options::default()).unwrap();
    let rule = gauss_legendre(FACE_POINTS).unwrap();
    let (models, conditions) = layer(&p, &s, 0, &rule).unwrap();
    let corners = corner_fluxes(&p, 0, &models, &conditions).unwrap();
    // 左下の角の点（番号 0）と、ノード 0 の展開の角の点の勾配
    let m = &p.materials[0];
    let pins = m.pins.as_ref().unwrap();
    let mut c = conditions[0].clone();
    c.corners = corner_values(&corners, 4, 0, 0, pins);
    let node = models[0].solve(&c);
    let [a, b] = c.half;
    let grad = node.gradient(-a, -b);
    for g in 0..2 {
        let outward = m.diffusion[g] * (grad[g][0] + grad[g][1]);
        let expected = 0.5 * corners[0][g];
        assert!(
            (outward - expected).abs() < 1e-8 * expected.abs(),
            "group {g}: {outward} vs {expected}"
        );
    }
}

#[test]
fn layers_without_fission_are_skipped() {
    // 下の層が反射体（核分裂も形状関数もない）。その層の出力は空で、炉心全体の表は燃料の層だけから作る
    let mut p = square(43.0, 2, 4, Boundary::Reflective);
    p.mesh.dz = vec![20.0, 30.0];
    let mut water = p.materials[0].clone();
    water.nu_fission = vec![0.0; 2];
    water.fission = vec![0.0; 2];
    water.chi = vec![0.0; 2];
    water.pins = None;
    p.materials.push(water);
    p.node_material = vec![1, 1, 1, 1, 0, 0, 0, 0];
    let s = solve(&p, &Options::default()).unwrap();
    let maps = reconstruct(&p, &s).unwrap();
    assert!(maps[..4].iter().all(Vec::is_empty));
    assert!(maps[4..].iter().all(|m| m.len() == 4));
    let rows = core_rows(&p, &maps).unwrap();
    assert_eq!((rows.len(), rows[0].len()), (8, 8));
    let mean = rows.iter().flatten().sum::<f64>() / 64.0;
    assert!((mean - 1.0).abs() < 1e-12);
}

#[test]
fn core_rows_merge_halves_of_split_cells() {
    // 2 × 1 の部分ノード（ADR-0108）。左のノードの +x の辺と右のノードの -x の辺のセルは、1つのセルを半分ずつ分けた
    // もので、面積で重み付けして1つに戻す（幅 0.5 と 0.3 の半分のセル、出力の密度 2 と 4 → (2·0.5 + 4·0.3)/0.8 = 2.75）
    let mut left = fuel(2.0, 2);
    let mut right = fuel(2.0, 2);
    if let (Some(l), Some(r)) = (left.pins.as_mut(), right.pins.as_mut()) {
        l.dx = vec![1.5, 0.5];
        l.halves = vec![false, true, false, false];
        r.dx = vec![0.3, 1.7];
        r.halves = vec![true, false, false, false];
    }
    let p = Problem {
        mesh: Mesh {
            dx: vec![2.0, 2.0],
            dy: vec![2.0],
            dz: vec![1.0],
        },
        materials: vec![left, right],
        node_material: vec![0, 1],
        boundaries: [Boundary::Reflective; 6],
    };
    let maps = vec![
        vec![vec![1.0, 2.0], vec![1.0, 2.0]],
        vec![vec![4.0, 1.0], vec![4.0, 1.0]],
    ];
    let rows = core_rows(&p, &maps).unwrap();
    let mean = (1.0 + 2.75 + 1.0) / 3.0;
    for row in &rows {
        for (got, want) in row.iter().zip([1.0, 2.75, 1.0]) {
            assert!((got - want / mean).abs() < 1e-14, "{rows:?}");
        }
    }
    // 片方の辺だけに印があればエラー
    let mut one_sided = p.clone();
    if let Some(r) = one_sided.materials[1].pins.as_mut() {
        r.halves = vec![false; 4];
    }
    assert!(core_rows(&one_sided, &maps).is_err());
}
