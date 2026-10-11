//! 角の点のバランス（corner point balance、ADR-0096）。
//!
//! 角の点を囲む小さな箱の中性子収支を、角の点ごと・群ごとの式とする。箱の辺を通って出ていく流れ
//! Σ_ノード (s_x J_x + s_y J_y)（s はノードが角の点の正の側にあれば +1、負の側なら −1。J はそのノードの展開の角の点
//! での値）と、外周の境界の辺を通る流れ（比 × ノードの均質な角の値）の和を 0 とする。ノードの角の均質な値は
//! 角の点の未知数 × ノードの角の重みとし、ノードの展開は角の値について線形なので、未知数の連立一次方程式を組んで解く。

use super::node::{Conditions, NodeFlux, NodeModel};
use crate::lu::Lu;

/// 角の点のバランスに入るノード。ノードの並びは x が速く、角は（-x -y、+x -y、+x +y、-x +y）、面は（-x、+x、-y、+y）。
#[derive(Debug, Clone)]
pub struct CornerNode<'a> {
    pub model: &'a NodeModel,
    /// ノードの条件（角の値は使わない）
    pub conditions: &'a Conditions,
    /// [群] の拡散係数
    pub diffusion: [f64; 2],
    /// [角][群] の、ノードの均質な角の値 / 角の点の未知数
    pub weights: [[f64; 2]; 4],
    /// [面][群] の、外周の辺を通って出ていく流れ / ノードの均質な角の値（真空の辺の比 r × 面の ADF。ADR-0131）。
    /// 外周でない面と、流れのない外周の面（反射）は `None`
    pub outflow: [Option<[f64; 2]>; 4],
}

/// nx × ny のノードの層の、[角の点（(nx + 1) × (ny + 1)、x が速い）][群] の未知数。`zero` の角の点（中性子束ゼロの
/// 外周の角）は 0 とする。ノードの均質な角の値は未知数 × 重みで、外周の辺の流れもそれに比例するので、式は
/// 均質な角の値だけで決まる。core は未知数を角の点の不均質な中性子束（重み 1/CDF）、lattice は均質な角の値の見込みに
/// 掛ける大きさ（重みは見込み。ADR-0131）とする
pub fn corner_balance(
    nx: usize,
    ny: usize,
    nodes: &[CornerNode],
    zero: &[bool],
) -> Result<Vec<[f64; 2]>, String> {
    let points = (nx + 1) * (ny + 1);
    if nodes.len() != nx * ny || zero.len() != points {
        return Err(format!(
            "corner balance: {} nodes and {} zero marks for a {nx} × {ny} layer",
            nodes.len(),
            zero.len(),
        ));
    }
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
            let node = &nodes[i + nx * j];
            let (model, base) = (node.model, node.conditions);
            // 角の値が 0 のときの勾配と、角 c の群 h の値を 1 にしたときの増分。[角の点の角][群][軸]
            let gradients = |flux: &NodeFlux| -> [[[f64; 2]; 2]; 4] {
                let [a, b] = base.half;
                [(-a, -b), (a, -b), (a, b), (-a, b)].map(|(x, y)| flux.gradient(x, y))
            };
            let mut zeroed = base.clone();
            zeroed.corners = [[0.0; 2]; 4];
            let zero_gradients = gradients(&model.solve(&zeroed));
            let mut response = [[[[[0.0; 2]; 2]; 4]; 2]; 4];
            for (c, per_corner) in response.iter_mut().enumerate() {
                for (h, per_group) in per_corner.iter_mut().enumerate() {
                    let mut unit = zeroed.clone();
                    unit.corners[c][h] = 1.0;
                    let moved = gradients(&model.solve(&unit));
                    for (r, (m, z)) in per_group.iter_mut().zip(moved.iter().zip(&zero_gradients)) {
                        for (rg, (mg, zg)) in r.iter_mut().zip(m.iter().zip(z)) {
                            for (ra, (ma, za)) in rg.iter_mut().zip(mg.iter().zip(zg)) {
                                *ra = ma - za;
                            }
                        }
                    }
                }
            }
            let d = node.diffusion;
            // ノードの角 `at` の点の式に、このノードの寄与 s_x J_x + s_y J_y（J = −D ∇φ）を足す
            for at in 0..4 {
                let (p, sx, sy) = point(i, j, at);
                for g in 0..2 {
                    let row = 2 * p + g;
                    let flow = |grad: [f64; 2]| -d[g] * (sx * grad[0] + sy * grad[1]);
                    rhs[row] -= flow(zero_gradients[at][g]);
                    for (c, per_corner) in response.iter().enumerate() {
                        let (q, _, _) = point(i, j, c);
                        for (h, per_group) in per_corner.iter().enumerate() {
                            matrix[row * unknowns + 2 * q + h] +=
                                flow(per_group[at][g]) * node.weights[c][h];
                        }
                    }
                }
            }
            // 外周の辺を通る流れ（比 × ノードの均質な角の値）
            for (face, ratio) in node.outflow.iter().enumerate() {
                let Some(ratio) = ratio else { continue };
                // この辺に接するノードの角
                let touching: [usize; 2] = match face {
                    0 => [0, 3],
                    1 => [1, 2],
                    2 => [0, 1],
                    _ => [3, 2],
                };
                for c in touching {
                    let (p, _, _) = point(i, j, c);
                    for (g, r) in ratio.iter().enumerate() {
                        let row = 2 * p + g;
                        matrix[row * unknowns + row] += r * node.weights[c][g];
                    }
                }
            }
        }
    }
    for (p, _) in zero.iter().enumerate().filter(|(_, z)| **z) {
        for g in 0..2 {
            let row = 2 * p + g;
            matrix[row * unknowns..(row + 1) * unknowns].fill(0.0);
            matrix[row * unknowns + row] = 1.0;
            rhs[row] = 0.0;
        }
    }
    let lu = Lu::new(matrix, unknowns).ok_or("the corner point balance is singular")?;
    lu.solve(&mut rhs);
    if rhs.iter().any(|v| !v.is_finite()) {
        return Err("the corner point balance gave non-finite corner fluxes".into());
    }
    Ok((0..points).map(|p| [rhs[2 * p], rhs[2 * p + 1]]).collect())
}

/// nx × ny のノードの層で、角の点 (ci, cj) を囲むノード（[i + nx j]）と、そのノードの、その点の角（-x -y、+x -y、
/// +x +y、-x +y）。層の外は除く
pub fn nodes_at_point(nx: usize, ny: usize, ci: usize, cj: usize) -> Vec<(usize, usize)> {
    // 角の点の左下、右下、右上、左上のノードと、それぞれの +x +y、-x +y、-x -y、+x -y の角
    [(1, 1, 2), (0, 1, 3), (0, 0, 0), (1, 0, 1)]
        .into_iter()
        .filter_map(|(di, dj, corner)| {
            let (i, j) = (ci.checked_sub(di)?, cj.checked_sub(dj)?);
            (i < nx && j < ny).then_some((i + nx * j, corner))
        })
        .collect()
}

/// 燃料の CDF で割った角の大きさ H_p: 角の点に接する燃料のノードの CDF の調和平均（ADR-0131）。lattice は CDF を
/// これで割り、core はこれを掛ける。燃料がなければ `None`
pub fn fuel_harmonic_mean(fuel: &[f64]) -> Option<f64> {
    (!fuel.is_empty()).then(|| fuel.len() as f64 / fuel.iter().map(|f| 1.0 / f).sum::<f64>())
}

#[cfg(test)]
mod tests;
