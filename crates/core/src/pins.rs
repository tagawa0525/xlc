//! ピン出力の再構成（ADR-0096）。
//!
//! 軸方向の層ごとに、各ノードの中の均質な中性子束を、NEM の解のノード平均、面の正味の中性子流、軸方向の漏れ（一様な
//! 線源）と、4 つの角の均質な中性子束から求める（`node`）。角の中性子束は、角の点を囲む小さな箱の中性子収支
//! （corner point balance）で決める。箱の辺を通って出ていく流れ Σ_ノード (s_x J_x + s_y J_y) と、外周の境界の辺を
//! 通る流れ（反射は 0、真空は J = r f φ）の和を 0 とする。s はノードが角の点の正の側にあれば +1、負の側なら −1 で、
//! J はそのノードの展開の角の点での値。不均質な角の中性子束 Φ は角を囲むノードで共通とし、各ノードの均質な値は
//! Φ / CDF とする。ノードの展開は角の値について線形なので、Φ の連立一次方程式を組んで解く。
//!
//! セルの出力は Σ_g Σf,g φ̄_g（セルの均質な中性子束の平均）f_g（形状関数）。

use std::collections::HashMap;
use std::rc::Rc;

use xlc_format::fewgroup::{FewGroup, Pins};
use xlc_numeric::lu::Lu;
use xlc_numeric::quadrature::gauss_legendre;

use crate::nem::{Boundary, Problem, Solution};
use basis::FACE_POINTS;
use node::{Conditions, NodeFlux, NodeModel};

mod basis;
mod node;

/// セルの平均の Gauss-Legendre の点の数（各方向）。ピンのセル（1.26 cm）の中で κh ≲ 1 なので、6 点の誤差は 1e-12 より小さい
const CELL_POINTS: usize = 6;

/// [ノード] のセルの出力 [行（-y から）][列（-x から）]。ノードの並びは NEM と同じ（x が最も速く、次に y、z）。
/// 核分裂のない軸方向の層（上下の反射体）は出力が 0 なので、空の表にする。核分裂のある層では、どのノードも形状関数を
/// 持つこと（径方向の反射体の角の CDF は決めていない。#0020）
pub fn reconstruct(problem: &Problem, solution: &Solution) -> Result<Vec<Vec<Vec<f64>>>, String> {
    let face_rule = gauss_legendre(FACE_POINTS).map_err(|e| e.to_string())?;
    let cell_rule = gauss_legendre(CELL_POINTS).map_err(|e| e.to_string())?;
    let (nx, ny, nz) = (
        problem.mesh.dx.len(),
        problem.mesh.dy.len(),
        problem.mesh.dz.len(),
    );
    let mut out = vec![Vec::new(); nx * ny * nz];
    for k in 0..nz {
        let fissile = (0..nx * ny).any(|local| {
            let n = local + nx * ny * k;
            problem
                .node_material
                .get(n)
                .and_then(|&m| problem.materials.get(m))
                .is_some_and(|m| m.nu_fission.iter().any(|&v| v != 0.0))
        });
        if !fissile {
            continue;
        }
        let (models, conditions) = layer(problem, solution, k, &face_rule)?;
        let corners = corner_fluxes(problem, k, &models, &conditions)?;
        for j in 0..ny {
            for i in 0..nx {
                let local = i + nx * j;
                let n = local + nx * ny * k;
                let (m, pins) = node_constants(problem, n)?;
                let mut c = conditions[local].clone();
                c.corners = corner_values(&corners, nx, i, j, pins);
                let node = models[local].solve(&c);
                let map = cell_powers(&node, m, pins, c.half, &cell_rule);
                if map.iter().flatten().any(|v| !v.is_finite()) {
                    return Err(format!("node {n}: a reconstructed pin power is not finite"));
                }
                out[n] = map;
            }
        }
    }
    Ok(out)
}

/// ノード n の少数群定数と形状関数
fn node_constants(problem: &Problem, n: usize) -> Result<(&FewGroup, &Pins), String> {
    let m = problem
        .node_material
        .get(n)
        .and_then(|&m| problem.materials.get(m))
        .ok_or_else(|| format!("node {n} has no material"))?;
    let pins = m.pins.as_ref().ok_or_else(|| {
        format!(
            "node {n} has no form functions (pins); reflector nodes next to fuel wait for #0020"
        )
    })?;
    Ok((m, pins))
}

/// 層 k の [ノード（層の中で i + nx j）] のモデルと、角の値を 0 にした条件。モデルは物質とノードの大きさで決まるので、
/// 同じ組のノードで共有する
fn layer(
    problem: &Problem,
    solution: &Solution,
    k: usize,
    rule: &[(f64, f64)],
) -> Result<(Vec<Rc<NodeModel>>, Vec<Conditions>), String> {
    let mesh = &problem.mesh;
    let (nx, ny) = (mesh.dx.len(), mesh.dy.len());
    let mut cache: HashMap<(usize, u64, u64), Rc<NodeModel>> = HashMap::new();
    let mut models = Vec::with_capacity(nx * ny);
    let mut conditions = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let n = i + nx * (j + ny * k);
            let (m, pins) = node_constants(problem, n)?;
            for (sizes, width, name) in [
                (&pins.dx, mesh.dx[i], "width"),
                (&pins.dy, mesh.dy[j], "height"),
            ] {
                let sum: f64 = sizes.iter().sum();
                if (sum - width).abs() > 1e-9 * width {
                    return Err(format!(
                        "node {n}: the cells of the form functions have a {name} of {sum} cm, \
                         the node {width} cm"
                    ));
                }
            }
            let half = [0.5 * mesh.dx[i], 0.5 * mesh.dy[j]];
            let key = (
                problem.node_material[n],
                half[0].to_bits(),
                half[1].to_bits(),
            );
            let model = match cache.get(&key) {
                Some(model) => Rc::clone(model),
                None => {
                    let model = Rc::new(
                        NodeModel::new(m, solution.k, half, rule)
                            .map_err(|e| format!("node {n}: {e}"))?,
                    );
                    cache.insert(key, Rc::clone(&model));
                    model
                }
            };
            models.push(model);
            let current = |d: usize, side: usize| -> [f64; 2] {
                std::array::from_fn(|g| solution.currents[d][n][g][side])
            };
            let (z_minus, z_plus) = (current(2, 0), current(2, 1));
            conditions.push(Conditions {
                half,
                average: [solution.flux[n][0], solution.flux[n][1]],
                corners: [[0.0; 2]; 4],
                currents: [current(0, 0), current(0, 1), current(1, 0), current(1, 1)],
                axial_leakage: std::array::from_fn(|g| (z_plus[g] - z_minus[g]) / mesh.dz[k]),
            });
        }
    }
    Ok((models, conditions))
}

/// 炉心全体のセルの出力 [行（上から）][列（左から）]。`maps` は `reconstruct` のノードごとの出力で、軸方向の層は
/// 厚さで重みを付けて足す。すべて 0 の行と列（集合体の間の余白など）を除き、核分裂のあるセルの平均で 1 に規格化する
/// （`lattice` の `pin-powers` と同じ規格化）
pub fn core_rows(problem: &Problem, maps: &[Vec<Vec<f64>>]) -> Result<Vec<Vec<f64>>, String> {
    let mesh = &problem.mesh;
    let (nx, ny, nz) = (mesh.dx.len(), mesh.dy.len(), mesh.dz.len());
    if maps.len() != nx * ny * nz {
        return Err(format!(
            "{} pin maps for {} nodes",
            maps.len(),
            nx * ny * nz
        ));
    }
    // 層の和（ノードは層の中で i + nx j）。核分裂のない層の空の表は足さない
    let mut radial: Vec<Option<Vec<Vec<f64>>>> = vec![None; nx * ny];
    for (n, map) in maps.iter().enumerate() {
        if map.is_empty() {
            continue;
        }
        let (local, k) = (n % (nx * ny), n / (nx * ny));
        let sum = radial[local]
            .get_or_insert_with(|| map.iter().map(|row| vec![0.0; row.len()]).collect());
        if map.len() != sum.len() || map.iter().zip(sum.iter()).any(|(a, b)| a.len() != b.len()) {
            return Err(format!(
                "node {n}: its pin map differs in shape from its other layers"
            ));
        }
        for (s, v) in sum.iter_mut().flatten().zip(map.iter().flatten()) {
            *s += mesh.dz[k] * v;
        }
    }
    let radial: Vec<Vec<Vec<f64>>> = radial
        .into_iter()
        .enumerate()
        .map(|(local, map)| {
            map.ok_or_else(|| format!("node column {local} has no pin map in any layer"))
        })
        .collect::<Result<_, _>>()?;
    // ノードの行ごとにセルの行の数が、ノードの列ごとにセルの列の数がそろうこと
    let rows_of = |j: usize| radial[nx * j].len();
    let columns_of = |i: usize| radial[i].first().map_or(0, Vec::len);
    for j in 0..ny {
        for i in 0..nx {
            let map = &radial[i + nx * j];
            if map.len() != rows_of(j) || map.iter().any(|r| r.len() != columns_of(i)) {
                return Err(format!(
                    "node ({i}, {j}) has a pin map of a different size from its row or column of nodes"
                ));
            }
        }
    }
    // 上の行から: ノードの行 j を上から、ノードの中の行（-y から）も上から
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for j in (0..ny).rev() {
        for r in (0..rows_of(j)).rev() {
            rows.push(
                (0..nx)
                    .flat_map(|i| radial[i + nx * j][r].iter().copied())
                    .collect(),
            );
        }
    }
    let columns = rows.first().map_or(0, Vec::len);
    let kept: Vec<usize> = (0..columns)
        .filter(|&c| rows.iter().any(|r| r[c] != 0.0))
        .collect();
    let rows: Vec<Vec<f64>> = rows
        .into_iter()
        .filter(|r| r.iter().any(|&v| v != 0.0))
        .map(|r| kept.iter().map(|&c| r[c]).collect())
        .collect();
    let fissile: Vec<f64> = rows
        .iter()
        .flatten()
        .copied()
        .filter(|&v| v > 0.0)
        .collect();
    if fissile.is_empty() {
        return Err("no cell has fission".into());
    }
    let mean = fissile.iter().sum::<f64>() / fissile.len() as f64;
    Ok(rows
        .into_iter()
        .map(|r| r.into_iter().map(|v| v / mean).collect())
        .collect())
}

/// ノード (i, j) の 4 つの角の均質な中性子束 Φ / CDF
fn corner_values(
    corners: &[[f64; 2]],
    nx: usize,
    i: usize,
    j: usize,
    pins: &Pins,
) -> [[f64; 2]; 4] {
    let at = |ci: usize, cj: usize| corners[ci + (nx + 1) * cj];
    let points = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)];
    std::array::from_fn(|c| std::array::from_fn(|g| points[c][g] / pins.corners[c][g]))
}

/// 層 k の [角の点（(nx + 1) × (ny + 1)、x が速い）][群] の不均質な中性子束 Φ
fn corner_fluxes(
    problem: &Problem,
    k: usize,
    models: &[Rc<NodeModel>],
    conditions: &[Conditions],
) -> Result<Vec<[f64; 2]>, String> {
    let (nx, ny) = (problem.mesh.dx.len(), problem.mesh.dy.len());
    let points = (nx + 1) * (ny + 1);
    let unknowns = 2 * points;
    let mut matrix = vec![0.0; unknowns * unknowns];
    let mut rhs = vec![0.0; unknowns];
    // ノードの角 c（-x -y、+x -y、+x +y、-x +y）の、角の点の番号と、ノードが角の点の +x 側か、+y 側か
    let point = |i: usize, j: usize, c: usize| -> (usize, f64, f64) {
        let (ci, cj) = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)][c];
        let sx = if ci == i { 1.0 } else { -1.0 };
        let sy = if cj == j { 1.0 } else { -1.0 };
        (ci + (nx + 1) * cj, sx, sy)
    };
    for j in 0..ny {
        for i in 0..nx {
            let local = i + nx * j;
            let n = local + nx * ny * k;
            let (m, pins) = node_constants(problem, n)?;
            let model = &models[local];
            let base = conditions[local].clone();
            // 角の値が 0 のときの勾配と、角 c の群 h の値を 1 にしたときの増分。[角の点の角][群][軸]
            let gradients = |node: &NodeFlux| -> [[[f64; 2]; 2]; 4] {
                let [a, b] = base.half;
                [(-a, -b), (a, -b), (a, b), (-a, b)].map(|(x, y)| node.gradient(x, y))
            };
            let zero = gradients(&model.solve(&base));
            let mut response = [[[[[0.0; 2]; 2]; 4]; 2]; 4];
            for (c, per_corner) in response.iter_mut().enumerate() {
                for (h, per_group) in per_corner.iter_mut().enumerate() {
                    let mut unit = base.clone();
                    unit.corners[c][h] = 1.0;
                    let moved = gradients(&model.solve(&unit));
                    for (r, (m, z)) in per_group.iter_mut().zip(moved.iter().zip(&zero)) {
                        for (rg, (mg, zg)) in r.iter_mut().zip(m.iter().zip(z)) {
                            for (ra, (ma, za)) in rg.iter_mut().zip(mg.iter().zip(zg)) {
                                *ra = ma - za;
                            }
                        }
                    }
                }
            }
            let d = [m.diffusion[0], m.diffusion[1]];
            // ノードの角 `at` の点の式に、このノードの寄与 s_x J_x + s_y J_y（J = −D ∇φ）を足す
            for at in 0..4 {
                let (p, sx, sy) = point(i, j, at);
                for g in 0..2 {
                    let row = 2 * p + g;
                    let flow = |grad: [f64; 2]| -d[g] * (sx * grad[0] + sy * grad[1]);
                    rhs[row] -= flow(zero[at][g]);
                    for (c, per_corner) in response.iter().enumerate() {
                        let (q, _, _) = point(i, j, c);
                        for (h, per_group) in per_corner.iter().enumerate() {
                            // 均質な角の値は Φ / CDF
                            matrix[row * unknowns + 2 * q + h] +=
                                flow(per_group[at][g]) / pins.corners[c][h];
                        }
                    }
                }
            }
            // 外周の真空の辺を通る流れ。真空の辺の比 r は不均質な中性子束に対する比（ADR-0101）で、角の点の不均質な
            // 中性子束は Φ なので r Φ。反射は 0
            for face in 0..4 {
                let on_boundary = match face {
                    0 => i == 0,
                    1 => i == nx - 1,
                    2 => j == 0,
                    _ => j == ny - 1,
                };
                let ratio = |g: usize| -> Result<Option<f64>, String> {
                    match problem.boundaries[face] {
                        Boundary::Reflective | Boundary::ZeroFlux => Ok(None),
                        Boundary::Marshak => Ok(Some(0.5)),
                        // 辺（左、右、上、下）の番号は、面（-x、+x、-y、+y）の順に 0、1、3、2
                        Boundary::Vacuum => m.vacuum_ratio[[0, 1, 3, 2][face]]
                            .as_ref()
                            .map(|r| Some(r[g]))
                            .ok_or_else(|| format!("node {n}: vacuum side without a ratio")),
                    }
                };
                if !on_boundary {
                    continue;
                }
                // この辺に接するノードの角
                let touching: [usize; 2] = match face {
                    0 => [0, 3],
                    1 => [1, 2],
                    2 => [0, 1],
                    _ => [3, 2],
                };
                for c in touching {
                    let (p, _, _) = point(i, j, c);
                    for g in 0..2 {
                        if let Some(r) = ratio(g)? {
                            let row = 2 * p + g;
                            matrix[row * unknowns + row] += r;
                        }
                    }
                }
            }
        }
    }
    // 中性子束ゼロの外周の角は Φ = 0
    let b = &problem.boundaries;
    for cj in 0..=ny {
        for ci in 0..=nx {
            let zero_flux = (ci == 0 && b[0] == Boundary::ZeroFlux)
                || (ci == nx && b[1] == Boundary::ZeroFlux)
                || (cj == 0 && b[2] == Boundary::ZeroFlux)
                || (cj == ny && b[3] == Boundary::ZeroFlux);
            if zero_flux {
                for g in 0..2 {
                    let row = 2 * (ci + (nx + 1) * cj) + g;
                    matrix[row * unknowns..(row + 1) * unknowns].fill(0.0);
                    matrix[row * unknowns + row] = 1.0;
                    rhs[row] = 0.0;
                }
            }
        }
    }
    let lu = Lu::new(matrix, unknowns).ok_or("the corner point balance is singular")?;
    lu.solve(&mut rhs);
    if rhs.iter().any(|v| !v.is_finite()) {
        return Err("the corner point balance gave non-finite corner fluxes".into());
    }
    Ok((0..points).map(|p| [rhs[2 * p], rhs[2 * p + 1]]).collect())
}

/// セルの出力 Σ_g Σf,g φ̄_g f_g（[行（-y から）][列]）
fn cell_powers(
    node: &NodeFlux,
    m: &FewGroup,
    pins: &Pins,
    half: [f64; 2],
    rule: &[(f64, f64)],
) -> Vec<Vec<f64>> {
    let edges = |sizes: &[f64], start: f64| -> Vec<f64> {
        let mut e = Vec::with_capacity(sizes.len() + 1);
        let mut x = start;
        e.push(x);
        for s in sizes {
            x += s;
            e.push(x);
        }
        e
    };
    let (xs, ys) = (edges(&pins.dx, -half[0]), edges(&pins.dy, -half[1]));
    (0..pins.dy.len())
        .map(|j| {
            (0..pins.dx.len())
                .map(|i| {
                    let phi = node.average(xs[i], xs[i + 1], ys[j], ys[j + 1], rule);
                    (0..2)
                        .map(|g| m.fission[g] * phi[g] * pins.form[g][j][i])
                        .sum()
                })
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests;
