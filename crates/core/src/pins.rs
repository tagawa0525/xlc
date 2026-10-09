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

use xlc_format::fewgroup::{FewGroup, Pins};
use xlc_numeric::lu::Lu;
use xlc_numeric::quadrature::gauss_legendre;

use crate::nem::{Boundary, Problem, Solution};
use basis::FACE_POINTS;
use node::{Conditions, NodeFlux, NodeModel};

mod basis;
mod node;

/// 辺（左、右、上、下）の番号。面（-x、+x、-y、+y）の順に
const SIDES: [usize; 4] = [0, 1, 3, 2];

/// [ノード] のセルの出力 [行（-y から）][列（-x から）]。ノードの並びは NEM と同じ（x が最も速く、次に y、z）
pub fn reconstruct(problem: &Problem, solution: &Solution) -> Result<Vec<Vec<Vec<f64>>>, String> {
    let rule = gauss_legendre(FACE_POINTS).map_err(|e| e.to_string())?;
    let mesh = &problem.mesh;
    let (nx, ny, nz) = (mesh.dx.len(), mesh.dy.len(), mesh.dz.len());
    let index = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
    let material = |n: usize| -> Result<(&FewGroup, &Pins), String> {
        let m = problem
            .node_material
            .get(n)
            .and_then(|&m| problem.materials.get(m))
            .ok_or_else(|| format!("node {n} has no material"))?;
        let pins = m
            .pins
            .as_ref()
            .ok_or_else(|| format!("node {n} has no form functions (pins)"))?;
        Ok((m, pins))
    };
    let mut out = vec![Vec::new(); nx * ny * nz];
    for k in 0..nz {
        // [ノード（層の中で i + nx j）] のモデルと条件
        let mut models = Vec::with_capacity(nx * ny);
        let mut conditions = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let n = index(i, j, k);
                let (m, pins) = material(n)?;
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
                models.push(
                    NodeModel::new(m, solution.k, half, &rule)
                        .map_err(|e| format!("node {n}: {e}"))?,
                );
                let current = |d: usize, side: usize| -> [f64; 2] {
                    std::array::from_fn(|g| solution.currents[d][n][g][side])
                };
                let z = current(2, 1);
                let z_minus = current(2, 0);
                conditions.push(Conditions {
                    half,
                    average: [solution.flux[n][0], solution.flux[n][1]],
                    corners: [[0.0; 2]; 4],
                    currents: [current(0, 0), current(0, 1), current(1, 0), current(1, 1)],
                    axial_leakage: std::array::from_fn(|g| (z[g] - z_minus[g]) / mesh.dz[k]),
                });
            }
        }
        let corners = corner_fluxes(problem, k, &models, &conditions)?;
        for j in 0..ny {
            for i in 0..nx {
                let local = i + nx * j;
                let n = index(i, j, k);
                let (m, pins) = material(n)?;
                let mut c = conditions[local].clone();
                c.corners = corner_values(&corners, nx, i, j, pins);
                let node = models[local].solve(&c);
                out[n] = cell_powers(&node, m, pins, c.half, &rule);
            }
        }
    }
    Ok(out)
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
    models: &[NodeModel],
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
    let zero_flux = |ci: usize, cj: usize| -> bool {
        let b = &problem.boundaries;
        (ci == 0 && b[0] == Boundary::ZeroFlux)
            || (ci == nx && b[1] == Boundary::ZeroFlux)
            || (cj == 0 && b[2] == Boundary::ZeroFlux)
            || (cj == ny && b[3] == Boundary::ZeroFlux)
    };
    for j in 0..ny {
        for i in 0..nx {
            let local = i + nx * j;
            let n = i + nx * (j + ny * k);
            let m = &problem.materials[problem.node_material[n]];
            let pins = m
                .pins
                .as_ref()
                .ok_or_else(|| format!("node {n} has no form functions (pins)"))?;
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
            // 外周の真空の辺を通る流れ r f φ（φ = Φ / CDF）。反射は 0
            for (face, side) in SIDES.iter().enumerate() {
                let on_boundary = match face {
                    0 => i == 0,
                    1 => i == nx - 1,
                    2 => j == 0,
                    _ => j == ny - 1,
                };
                if !on_boundary {
                    continue;
                }
                let ratio = |g: usize| -> Result<Option<f64>, String> {
                    match problem.boundaries[face] {
                        Boundary::Reflective | Boundary::ZeroFlux => Ok(None),
                        Boundary::Marshak => Ok(Some(0.5)),
                        Boundary::Vacuum => m.vacuum_ratio[*side]
                            .as_ref()
                            .map(|r| Some(r[g]))
                            .ok_or_else(|| format!("node {n}: vacuum side without a ratio")),
                    }
                };
                if m.normalized_adf[*side] {
                    return Err(format!(
                        "node {n}: an outer side has an f* (normalized) ADF"
                    ));
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
                            matrix[row * unknowns + row] +=
                                r * m.adf[*side][g] / pins.corners[c][g];
                        }
                    }
                }
            }
        }
    }
    // 中性子束ゼロの外周の角は Φ = 0
    for cj in 0..=ny {
        for ci in 0..=nx {
            if zero_flux(ci, cj) {
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
