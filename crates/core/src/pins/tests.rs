use super::*;
use crate::nem::{Boundary, Mesh, Options, Problem, solve};
use std::f64::consts::PI;
use xlc_format::fewgroup::{FewGroup, Pins};

/// どの位置にもノードを置く物質の並び
fn placed(materials: impl IntoIterator<Item = usize>) -> Vec<Option<usize>> {
    materials.into_iter().map(Some).collect()
}

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
        adf: vec![vec![1.0, 1.0]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: Some(Pins {
            dx: vec![width / cells as f64; cells],
            dy: vec![width / cells as f64; cells],
            form: vec![vec![vec![1.0; cells]; cells]; 2],
            corners: vec![vec![1.0, 1.0]; 4],
            normalized_corners: vec![false; 4],
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
        reflectors: Vec::new(),
        node_material: placed(vec![0; nodes * nodes]),
        boundaries,
        outline: None,
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
    p.node_material[3] = Some(1);
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
        reflectors: Vec::new(),
        node_material: placed(vec![0; 4]),
        boundaries: [Boundary::Reflective; 6],
        outline: None,
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
fn vacuum_corners_lose_r_times_the_homogeneous_corner_flux() {
    // 平板（x は Marshak、y と z は反射）の左下の角の点の箱は、ノード 0 だけを含み、箱の辺から出ていく流れ
    // J_x + J_y（J = −D∇φ）と、真空の辺から出ていく r f φ_c（r = 1/2、f は面の ADF、φ_c = Φ / CDF は均質な角の値）の
    // 和が 0。CDF（1.3）が ADF（1.2）と違っても、真空の辺の流れは面の条件と同じく ADF で決まる（ADR-0131）
    let mut p = square(86.0, 4, 2, Boundary::Marshak);
    p.mesh.dy = vec![86.0 / 4.0];
    p.node_material = placed(vec![0; 4]);
    p.boundaries[2] = Boundary::Reflective;
    p.boundaries[3] = Boundary::Reflective;
    p.materials[0].adf = vec![vec![1.2, 1.2]; 6];
    if let Some(pins) = p.materials[0].pins.as_mut() {
        pins.corners = vec![vec![1.3, 1.3]; 4];
    }
    let s = solve(&p, &Options::default()).unwrap();
    let rule = gauss_legendre(FACE_POINTS).unwrap();
    let constants = layer_constants(&p, 0).unwrap();
    let (models, conditions) = layer(&p, &s, 0, &constants, &rule).unwrap();
    let cdfs = corner_factors(&p, 0, &constants).unwrap();
    let corners = corner_fluxes(&p, 0, &constants, &models, &conditions, &cdfs).unwrap();
    // 左下の角の点（番号 0）と、ノード 0 の展開の角の点の勾配
    let m = &p.materials[0];
    let mut c = conditions[0].clone().unwrap();
    c.corners = corner_values(&corners, 4, 0, 0, cdfs[0].as_ref().unwrap());
    let node = models[0].as_ref().unwrap().solve(&c);
    let [a, b] = c.half;
    let grad = node.gradient(-a, -b);
    for g in 0..2 {
        let outward = m.diffusion[g] * (grad[g][0] + grad[g][1]);
        let expected = 0.5 * 1.2 * corners[0][g] / 1.3;
        assert!(
            (outward - expected).abs() < 1e-8 * expected.abs(),
            "group {g}: {outward} vs {expected}"
        );
    }
}

#[test]
fn vacuum_corners_take_the_ratio_and_adf_of_their_side() {
    // `vacuum` の外周の辺では、角の点の流れは、その辺（-y は下）の比 r と ADF f を均質な角の値に掛けた r f φ_c
    // （ADR-0131）。ほかの辺の r と ADF は効かない
    let mut p = square(86.0, 4, 2, Boundary::Vacuum);
    p.mesh.dx = vec![86.0 / 4.0];
    p.node_material = placed(vec![0; 4]);
    p.boundaries[0] = Boundary::Reflective;
    p.boundaries[1] = Boundary::Reflective;
    // 辺（左、右、上、下、−z、+z）
    p.materials[0].adf = vec![
        vec![1.0, 1.0],
        vec![1.0, 1.0],
        vec![0.9, 0.8],
        vec![1.2, 1.1],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    p.materials[0].vacuum_ratio = vec![
        None,
        None,
        Some(vec![0.45, 0.35]),
        Some(vec![0.3, 0.4]),
        None,
        None,
    ];
    if let Some(pins) = p.materials[0].pins.as_mut() {
        pins.corners = vec![vec![1.3, 1.25]; 4];
    }
    let s = solve(&p, &Options::default()).unwrap();
    let rule = gauss_legendre(FACE_POINTS).unwrap();
    let constants = layer_constants(&p, 0).unwrap();
    let (models, conditions) = layer(&p, &s, 0, &constants, &rule).unwrap();
    let cdfs = corner_factors(&p, 0, &constants).unwrap();
    let corners = corner_fluxes(&p, 0, &constants, &models, &conditions, &cdfs).unwrap();
    // 左下の角の点（番号 0）と、ノード 0 の展開の角の点の勾配
    let m = &p.materials[0];
    let mut c = conditions[0].clone().unwrap();
    c.corners = corner_values(&corners, 1, 0, 0, cdfs[0].as_ref().unwrap());
    let node = models[0].as_ref().unwrap().solve(&c);
    let [a, b] = c.half;
    let grad = node.gradient(-a, -b);
    for g in 0..2 {
        let outward = m.diffusion[g] * (grad[g][0] + grad[g][1]);
        let expected = [0.3, 0.4][g] * [1.2, 1.1][g] * corners[0][g] / [1.3, 1.25][g];
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
    p.node_material = placed(vec![1, 1, 1, 1, 0, 0, 0, 0]);
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
        reflectors: Vec::new(),
        node_material: placed(vec![0, 1]),
        boundaries: [Boundary::Reflective; 6],
        outline: None,
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
    // ノードのない位置（ADR-0132）に面する辺に印があっても、分け合う相手がないのでエラー
    let mut beside_empty = p.clone();
    beside_empty.mesh.dy = vec![2.0, 2.0];
    beside_empty.node_material = vec![Some(0), Some(1), Some(0), None];
    beside_empty.outline = Some(Boundary::Reflective);
    let mut maps_empty = maps.clone();
    maps_empty.push(maps[0].clone());
    maps_empty.push(Vec::new());
    let err = core_rows(&beside_empty, &maps_empty).unwrap_err();
    assert!(err.contains("split"), "{err}");
    // 印のない辺でも、同じ列の境でほかの行が分け合っていれば、ノードのない位置の 0 と平均することになるのでエラー
    let mut unmarked = beside_empty.clone();
    let mut plain = unmarked.materials[0].clone();
    if let Some(q) = plain.pins.as_mut() {
        q.halves = vec![false; 4];
    }
    unmarked.materials.push(plain);
    unmarked.node_material[2] = Some(2);
    let err = core_rows(&unmarked, &maps_empty).unwrap_err();
    assert!(err.contains("split"), "{err}");
}

#[test]
fn cell_sizes_match_to_rounding_between_nodes() {
    // 切ったセルの幅は、回した集合体と回していない集合体で最後の桁が違いうる（座標の引き算の丸め）。同じ列の
    // ノードのセルの幅は丸めの程度で比べ、それより大きく違えばエラー
    let with_dx = |dx: Vec<f64>| {
        let mut m = fuel(2.56, 2);
        if let Some(p) = m.pins.as_mut() {
            p.dx = dx;
            p.dy = vec![1.28, 1.28];
        }
        m
    };
    let mut p = Problem {
        mesh: Mesh {
            dx: vec![2.56],
            dy: vec![2.56, 2.56],
            dz: vec![1.0],
        },
        materials: vec![
            with_dx(vec![1.93, 0.6300000000000008]),
            with_dx(vec![1.93, 0.629999999999999]),
        ],
        reflectors: Vec::new(),
        node_material: placed(vec![0, 1]),
        boundaries: [Boundary::Reflective; 6],
        outline: None,
    };
    let maps = vec![vec![vec![1.0, 2.0]; 2]; 2];
    assert!(core_rows(&p, &maps).is_ok());
    if let Some(pins) = p.materials[1].pins.as_mut() {
        pins.dx = vec![1.929, 0.631];
    }
    assert!(core_rows(&p, &maps).is_err());
}

/// 核分裂のない2群の反射体。形状関数は 0（`cells` × `cells`）、角の不連続因子は `corners`
fn reflector(width: f64, cells: usize, corners: Vec<Vec<f64>>, normalized: Vec<bool>) -> FewGroup {
    let (d1, d2) = (1.2, 0.2);
    let (a1, a2) = (0.001, 0.02);
    let s12 = 0.03;
    let total = [1.0 / (3.0 * d1), 1.0 / (3.0 * d2)];
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: 1.0,
        total: total.to_vec(),
        absorption: vec![a1, a2],
        nu_fission: vec![0.0, 0.0],
        fission: vec![0.0, 0.0],
        chi: vec![1.0, 0.0],
        diffusion: vec![d1, d2],
        scatter: vec![vec![total[0] - a1 - s12, s12], vec![0.0, total[1] - a2]],
        adf: vec![vec![1.0, 1.0]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: Some(Pins {
            dx: vec![width / cells as f64; cells],
            dy: vec![width / cells as f64; cells],
            form: vec![vec![vec![0.0; cells]; cells]; 2],
            corners,
            normalized_corners: normalized,
            halves: vec![false; 4],
        }),
    }
}

/// 燃料（-x）と反射体（+x）の 1 × 2 のノード（外周は鏡面反射）。燃料の角の不連続因子は `fuel_corners`
fn fuel_and_reflector(fuel_corners: Vec<Vec<f64>>, reflector: FewGroup) -> Problem {
    let mut f = fuel(21.5, 4);
    if let Some(p) = f.pins.as_mut() {
        p.corners = fuel_corners;
    }
    Problem {
        mesh: Mesh {
            dx: vec![21.5; 2],
            dy: vec![21.5],
            dz: vec![30.0],
        },
        materials: vec![f, reflector],
        reflectors: Vec::new(),
        node_material: placed(vec![0, 1]),
        boundaries: [Boundary::Reflective; 6],
        outline: None,
    }
}

#[test]
fn normalized_reflector_corners_take_the_fuel_corner_factor() {
    // 燃料の CDF で割った反射体の角（ADR-0131）の CDF は、C* × 角の点に接する燃料のノードの CDF（調和平均）。
    // 同じ値を割らずに書いた反射体と同じ再構成になり、反射体のノードの出力は 0
    let fuel_corners = vec![
        vec![1.1, 0.9],
        vec![1.2, 0.8],
        vec![1.25, 0.85],
        vec![1.05, 0.95],
    ];
    let star = [vec![0.8, 1.3], vec![0.7, 1.4]];
    let outer = vec![1.1, 0.7];
    let normalized = reflector(
        21.5,
        2,
        vec![
            star[0].clone(),
            outer.clone(),
            outer.clone(),
            star[1].clone(),
        ],
        vec![true, false, false, true],
    );
    // 反射体の -x -y の角は燃料の +x -y の角、-x +y の角は燃料の +x +y の角に接する
    let times = |s: &[f64], f: &[f64]| vec![s[0] * f[0], s[1] * f[1]];
    let plain = reflector(
        21.5,
        2,
        vec![
            times(&star[0], &fuel_corners[1]),
            outer.clone(),
            outer,
            times(&star[1], &fuel_corners[2]),
        ],
        vec![false; 4],
    );
    let run = |r: FewGroup| {
        let p = fuel_and_reflector(fuel_corners.clone(), r);
        let s = solve(&p, &Options::default()).unwrap();
        reconstruct(&p, &s).unwrap()
    };
    let (a, b) = (run(normalized), run(plain));
    for (x, y) in a[0].iter().flatten().zip(b[0].iter().flatten()) {
        assert!((x - y).abs() < 1e-12 * y.abs(), "{x} vs {y}");
    }
    assert!(a[1].iter().flatten().all(|&v| v == 0.0));
}

#[test]
fn normalized_corners_need_fuel_at_the_corner_point() {
    // 燃料の CDF で割った角の点に、燃料（核分裂のある）のノードがなければ、割っていない反射体の角があってもエラー
    let a = reflector(
        21.5,
        2,
        vec![vec![1.0, 1.0]; 4],
        vec![false, true, true, false],
    );
    let b = reflector(21.5, 2, vec![vec![1.0, 1.0]; 4], vec![false; 4]);
    let mut p = fuel_and_reflector(vec![vec![1.0, 1.0]; 4], a);
    p.materials.push(b);
    p.mesh.dx = vec![21.5; 3];
    p.node_material = placed(vec![0, 1, 2]);
    let s = solve(&p, &Options::default()).unwrap();
    let err = reconstruct(&p, &s).unwrap_err();
    assert!(err.contains("fuel"), "{err}");
}

/// 行列（上の行から）を反時計回りに 90° 回す
fn turned<T: Clone>(grid: &[Vec<T>]) -> Vec<Vec<T>> {
    let (rows, cols) = (grid.len(), grid[0].len());
    (0..cols)
        .map(|i| (0..rows).map(|j| grid[j][cols - 1 - i].clone()).collect())
        .collect()
}

#[test]
fn reflector_sites_turned_with_the_core_keep_k_and_pins() {
    // L 型の反射体（燃料と左と下の辺で接し、角の 3 つの点で燃料に接する。f* と C* の印を持つ）と燃料 3 体の 2 × 2 を、
    // 装荷図の回転（ADR-0105）で東西南北の 4 つの向きに置く。燃料も同じだけ回せば、k は向きによらず、ピン出力は
    // 炉心ごと回したものになる（辺、角、形状関数、f* と C* の印を同じ規則で回す）
    use crate::input::{Case, Cycle, Inventory, Layout, Plant};
    use std::path::Path;
    let cells = 3;
    let mut f = fuel(21.5, cells);
    f.adf = vec![
        vec![1.02, 0.97],
        vec![0.99, 1.05],
        vec![1.04, 0.95],
        vec![0.98, 1.03],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    if let Some(p) = f.pins.as_mut() {
        p.form = (0..2)
            .map(|g| {
                (0..cells)
                    .map(|j| {
                        (0..cells)
                            .map(|i| 1.0 + 0.05 * ((g + 3 * i + 7 * j) % 5) as f64)
                            .collect()
                    })
                    .collect()
            })
            .collect();
        p.corners = vec![
            vec![1.1, 0.9],
            vec![1.2, 0.8],
            vec![1.25, 0.85],
            vec![1.05, 0.95],
        ];
    }
    let mut r = reflector(
        21.5,
        cells,
        vec![
            vec![0.8, 1.3],
            vec![0.75, 1.35],
            vec![1.1, 0.7],
            vec![0.7, 1.4],
        ],
        vec![true, true, false, true],
    );
    r.adf = vec![
        vec![1.2, 0.4],
        vec![1.1, 0.6],
        vec![0.9, 1.2],
        vec![1.15, 0.45],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    r.normalized_adf = vec![true, false, false, true, false, false];
    let base = Path::new("dir");
    let inventory = Inventory::parse(
        "[design.F]\nbounds = [0.0, 30.0]\nlayers = [\"f.msgpack\"]\n\
         [assembly.J01]\ndesign = \"F\"\n[assembly.J02]\ndesign = \"F\"\n\
         [assembly.J03]\ndesign = \"F\"\n",
        base,
    )
    .unwrap();
    let case = Case::parse(
        "cycle = \"cycle.toml\"\naxial-mesh = [0.0, 30.0]\nradial-nodes = 1\n\
         boundary = { x = [\"reflective\", \"reflective\"], y = [\"reflective\", \"reflective\"], \
         z = [\"reflective\", \"reflective\"] }\n",
        base,
    )
    .unwrap();
    // 向き 0: 反射体は右上（東の列の北）、燃料は左上、左下、右下
    let mut sites = vec![
        vec![("*", "J01"), ("R", ".")],
        vec![("*", "J02"), ("*", "J03")],
    ];
    let mut reference: Option<(f64, Vec<Vec<f64>>)> = None;
    for turns in 0..4 {
        let angle = 90 * turns;
        let row = |cells: &[(&str, &str)], site: bool| {
            cells
                .iter()
                .map(|&(s, id)| {
                    let name = if site { s } else { id };
                    match name {
                        "*" | "." => name.to_string(),
                        _ => format!("{name}@{angle}"),
                    }
                })
                .collect::<Vec<_>>()
                .join("  ")
        };
        let quoted = |site: bool| {
            sites
                .iter()
                .map(|cells| format!("\"{}\"", row(cells, site)))
                .collect::<Vec<_>>()
                .join(", ")
        };
        let plant = Plant::parse(
            &format!(
                "pitch = 21.5\nheight = 30.0\nmap = [{}]\n[reflector.R]\nbounds = [0.0, 30.0]\n\
                 layers = [\"r.msgpack\"]\n",
                quoted(true)
            ),
            base,
        )
        .unwrap();
        let cycle = Cycle::parse(
            &format!(
                "plant = \"plant.toml\"\ninventory = \"assemblies.toml\"\nloading = [{}]\n",
                quoted(false)
            ),
            base,
        )
        .unwrap();
        let layout = Layout::new(&case, &cycle, &plant, &inventory).unwrap();
        let constants: Vec<Vec<FewGroup>> = layout
            .files
            .iter()
            .map(|file| match file.to_str() {
                Some("dir/f.msgpack") => vec![f.clone()],
                Some("dir/r.msgpack") => vec![r.clone()],
                other => panic!("{other:?}"),
            })
            .collect();
        let problem = layout.problem(&constants).unwrap();
        // 向きで反復の道筋が変わるので、許容差を詰めて収束させた解で比べる
        let tight = Options {
            k_tolerance: 1e-14,
            flux_tolerance: 1e-13,
            ..Options::default()
        };
        let solution = solve(&problem, &tight).unwrap();
        let maps = reconstruct(&problem, &solution).unwrap();
        let rows = core_rows(&problem, &maps).unwrap();
        match &reference {
            None => reference = Some((solution.k, rows)),
            Some((k, expected)) => {
                assert!(
                    (solution.k - k).abs() < 1e-12 * k,
                    "{angle}°: {} vs {k}",
                    solution.k
                );
                let mut expected = expected.clone();
                for _ in 0..turns {
                    expected = turned(&expected);
                }
                let shape = |m: &[Vec<f64>]| m.iter().map(Vec::len).collect::<Vec<_>>();
                assert_eq!(shape(&rows), shape(&expected), "{angle}°");
                for (a, b) in rows.iter().flatten().zip(expected.iter().flatten()) {
                    assert!(
                        (a - b).abs() < 1e-10 * b.abs().max(1.0),
                        "{angle}°: {a} vs {b}"
                    );
                }
            }
        }
        sites = turned(&sites);
    }
}

#[test]
fn outline_faces_reconstruct_like_cuboid_faces() {
    // 外形（ノードのない位置）に面する真空の辺の流れは、直方体の外周の真空の辺と同じ r f φ_c（ADR-0131、ADR-0132）。
    // 燃料と反射体の 1 × 2 で、反射体の +x を外周にした問題と、+x にノードのない位置を足した問題は同じ出力になる
    let mut r = reflector(21.5, 4, vec![vec![1.1, 0.9]; 4], vec![false; 4]);
    r.adf[1] = vec![1.3, 0.7];
    r.vacuum_ratio[1] = Some(vec![0.5, 0.35]);
    let mut cuboid = fuel_and_reflector(vec![vec![1.0, 1.0]; 4], r);
    cuboid.boundaries[1] = Boundary::Vacuum;
    let mut stepped = cuboid.clone();
    stepped.mesh.dx.push(21.5);
    stepped.node_material.push(None);
    stepped.boundaries[1] = Boundary::ZeroFlux;
    stepped.outline = Some(Boundary::Vacuum);
    let a = solve(&cuboid, &Options::default()).unwrap();
    let b = solve(&stepped, &Options::default()).unwrap();
    let (pa, pb) = (
        reconstruct(&cuboid, &a).unwrap(),
        reconstruct(&stepped, &b).unwrap(),
    );
    assert!(pb[2].is_empty());
    // ノードのある位置の表は同じ大きさ
    for n in 0..2 {
        assert_eq!(pa[n].len(), pb[n].len(), "node {n}");
        assert!(
            pa[n].iter().zip(&pb[n]).all(|(x, y)| x.len() == y.len()),
            "node {n}"
        );
    }
    assert_eq!(pa[0].len(), 4);
    for (x, y) in pa[0].iter().flatten().zip(pb[0].iter().flatten()) {
        assert!((x - y).abs() < 1e-9 * x, "{x} vs {y}");
    }
    let (ra, rb) = (
        core_rows(&cuboid, &pa).unwrap(),
        core_rows(&stepped, &pb).unwrap(),
    );
    assert_eq!(ra.len(), rb.len());
    assert!(ra.iter().zip(&rb).all(|(x, y)| x.len() == y.len()));
    for (x, y) in ra.iter().flatten().zip(rb.iter().flatten()) {
        assert!((x - y).abs() < 1e-9 * x.abs().max(1.0), "{x} vs {y}");
    }
}

#[test]
fn concave_outline_corners_keep_the_symmetry() {
    // 2 × 2 の燃料の位置のうち +x +y にノードを置かない L 字の炉心。外形の凹んだ角の点は 3 つのノードに囲まれ、
    // 2 本の半分の辺が外形に面する。体系は x と y の入れ替えで対称なので、(1, 0) と (0, 1) のノードの出力は転置で
    // 一致する。ノードのない位置の出力は空で、炉心全体の表では 0
    let mut p = square(43.0, 2, 4, Boundary::Marshak);
    p.boundaries[0] = Boundary::Reflective;
    p.boundaries[2] = Boundary::Reflective;
    p.node_material[3] = None;
    p.outline = Some(Boundary::Marshak);
    let s = solve(&p, &Options::default()).unwrap();
    let pins = reconstruct(&p, &s).unwrap();
    assert!(pins[3].is_empty());
    for (j, row) in pins[1].iter().enumerate() {
        for (i, &a) in row.iter().enumerate() {
            let b = pins[2][i][j];
            assert!((a - b).abs() < 1e-6 * a, "({i}, {j}): {a} vs {b}");
        }
    }
    let rows = core_rows(&p, &pins).unwrap();
    assert_eq!((rows.len(), rows[0].len()), (8, 8));
    // 上の行が +y なので、ノードのない +x +y の位置は上の 4 行の右の 4 列
    assert!(rows[..4].iter().all(|r| r[4..].iter().all(|&v| v == 0.0)));
}

#[test]
fn nodes_without_fission_take_the_cells_of_the_fuel_beside_them() {
    // 階段の外形の炉心では、燃料の列（行）に、セルの幅の違う反射体のノード（向きの違う I 型、角型）が並ぶ。核分裂の
    // ないノードのセルは 0 なので、その列と行のセルは燃料のノードから取り、反射体のノードは 0 で埋める
    let mut p = Problem {
        mesh: Mesh {
            dx: vec![21.5, 21.5],
            dy: vec![21.5, 21.5],
            dz: vec![1.0],
        },
        materials: vec![
            fuel(21.5, 4),
            reflector(21.5, 3, vec![vec![1.0, 1.0]; 4], vec![false; 4]),
        ],
        reflectors: Vec::new(),
        // (0, 0) が燃料、ほかの 3 つが反射体
        node_material: placed(vec![0, 1, 1, 1]),
        boundaries: [Boundary::Reflective; 6],
        outline: None,
    };
    let fuel_map: Vec<Vec<f64>> = (0..4)
        .map(|r| (0..4).map(|c| 1.0 + (4 * r + c) as f64).collect())
        .collect();
    let zeros = vec![vec![0.0; 3]; 3];
    let maps = vec![fuel_map.clone(), zeros.clone(), zeros.clone(), zeros];
    let rows = core_rows(&p, &maps).unwrap();
    let mean = fuel_map.iter().flatten().sum::<f64>() / 16.0;
    let expected: Vec<Vec<f64>> = fuel_map
        .iter()
        .rev()
        .map(|r| r.iter().map(|v| v / mean).collect())
        .collect();
    assert_eq!(rows, expected);
    // 燃料のない列と行はすべて 0 で表に残らないので、そこに並ぶ反射体のノードのセルの大きさは違ってよい
    p.materials
        .push(reflector(21.5, 2, vec![vec![1.0, 1.0]; 4], vec![false; 4]));
    p.node_material[3] = Some(2);
    let mut maps = maps;
    maps[3] = vec![vec![0.0; 2]; 2];
    assert_eq!(core_rows(&p, &maps).unwrap(), expected);
    // 核分裂のないノードの出力に 0 でない値があれば、0 で埋めずにエラー
    maps[1][0][0] = 1.0e-3;
    let err = core_rows(&p, &maps).unwrap_err();
    assert!(err.contains("no fission but nonzero"), "{err}");
    // 燃料のない列と行にある反射体のノードでも同じ
    maps[1][0][0] = 0.0;
    maps[3][0][0] = 1.0e-3;
    let err = core_rows(&p, &maps).unwrap_err();
    assert!(err.contains("no fission but nonzero"), "{err}");
}

#[test]
fn split_cells_between_fuel_and_reflector_are_errors() {
    // 燃料と反射体の 2 × 2。列 1 は燃料（下）と反射体（上）。上の行で燃料の列 0 と反射体の列 1 が分け合うセルを
    // 持つと、燃料の半分のセルを反射体の 0 と平均することになるのでエラー
    let mut f = fuel(21.5, 4);
    let mut r = reflector(21.5, 4, vec![vec![1.0, 1.0]; 4], vec![false; 4]);
    if let (Some(a), Some(b)) = (f.pins.as_mut(), r.pins.as_mut()) {
        a.halves[1] = true;
        a.halves[0] = true;
        b.halves[0] = true;
    }
    let p = Problem {
        mesh: Mesh {
            dx: vec![21.5, 21.5],
            dy: vec![21.5, 21.5],
            dz: vec![1.0],
        },
        materials: vec![f, r],
        reflectors: Vec::new(),
        node_material: placed(vec![0, 0, 0, 1]),
        boundaries: [Boundary::Reflective; 6],
        outline: None,
    };
    let ones = vec![vec![1.0; 4]; 4];
    let maps = vec![ones.clone(), ones.clone(), ones, vec![vec![0.0; 4]; 4]];
    let err = core_rows(&p, &maps).unwrap_err();
    assert!(err.contains("split cell"), "{err}");
}

#[test]
fn split_cells_between_columns_without_fission_are_dropped() {
    // 燃料と反射体 2 つの 3 × 1。核分裂のないノードだけの列はセルを持たず、反射体どうしが分け合うセル（`halves`）は
    // どちらも 0 なので、戻すものがなく表は燃料だけ
    let mut left = reflector(21.5, 4, vec![vec![1.0, 1.0]; 4], vec![false; 4]);
    let mut right = left.clone();
    if let (Some(a), Some(b)) = (left.pins.as_mut(), right.pins.as_mut()) {
        a.halves[1] = true;
        b.halves[0] = true;
    }
    let p = Problem {
        mesh: Mesh {
            dx: vec![21.5; 3],
            dy: vec![21.5],
            dz: vec![1.0],
        },
        materials: vec![fuel(21.5, 4), left, right],
        reflectors: Vec::new(),
        node_material: placed(vec![0, 1, 2]),
        boundaries: [Boundary::Reflective; 6],
        outline: None,
    };
    let zeros = vec![vec![0.0; 4]; 4];
    let maps = vec![vec![vec![1.0; 4]; 4], zeros.clone(), zeros];
    assert_eq!(core_rows(&p, &maps).unwrap(), vec![vec![1.0; 4]; 4]);
}

#[test]
fn layers_take_the_cells_of_the_fuel_in_any_layer() {
    // 2 × 1 × 2 のノード。層 0 は燃料が 2 つ、層 1 は燃料と、セルの数の違う反射体。列 1 のセルは、層 1 でも層 0 の
    // 燃料のものを使い、層の表の形がそろう
    let p = Problem {
        mesh: Mesh {
            dx: vec![21.5, 21.5],
            dy: vec![21.5],
            dz: vec![1.0, 1.0],
        },
        materials: vec![
            fuel(21.5, 4),
            reflector(21.5, 3, vec![vec![1.0, 1.0]; 4], vec![false; 4]),
        ],
        reflectors: Vec::new(),
        node_material: placed(vec![0, 0, 0, 1]),
        boundaries: [Boundary::Reflective; 6],
        outline: None,
    };
    let ones = vec![vec![1.0; 4]; 4];
    let maps = vec![ones.clone(), ones.clone(), ones, vec![vec![0.0; 3]; 3]];
    let rows = core_rows(&p, &maps).unwrap();
    // 層の和は、列 0 のセルが 2、列 1 のセルが 1（層 1 の反射体は 0）。平均 1.5 で規格化する
    assert_eq!(
        rows,
        vec![[vec![2.0 / 1.5; 4], vec![1.0 / 1.5; 4]].concat(); 4]
    );
}

#[test]
fn positions_with_only_some_values_are_errors() {
    // 層の位置ごとの量（定数、モデル、条件、CDF）は、ノードのある位置ではすべて、ない位置ではどれも持たない。
    // 一部だけ持つのは配列の作り方の誤りなのでエラー
    let p = square(43.0, 1, 4, Boundary::Reflective);
    let s = solve(&p, &Options::default()).unwrap();
    let rule = gauss_legendre(FACE_POINTS).unwrap();
    let constants = layer_constants(&p, 0).unwrap();
    let (models, conditions) = layer(&p, &s, 0, &constants, &rule).unwrap();
    let cdfs = corner_factors(&p, 0, &constants).unwrap();
    assert!(
        parts(0, &constants, &models, &conditions, &cdfs)
            .unwrap()
            .is_some()
    );
    let none = vec![None];
    assert!(
        parts(0, &[None], &none, &[None], &[None])
            .unwrap()
            .is_none()
    );
    assert!(parts(0, &constants, &models, &conditions, &[None]).is_err());
}
