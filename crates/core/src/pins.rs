//! ピン出力の再構成（ADR-0096）。
//!
//! 軸方向の層ごとに、各ノードの中の均質な中性子束を、NEM の解のノード平均、面の正味の中性子流、軸方向の漏れ（一様な
//! 線源）と、4 つの角の均質な中性子束から求める（`node`）。角の中性子束は、角の点を囲む小さな箱の中性子収支
//! （corner point balance）で決める。箱の辺を通って出ていく流れ Σ_ノード (s_x J_x + s_y J_y) と、外周の境界の辺を
//! 通る流れ（反射は 0、真空は J = r f φ）の和を 0 とする。s はノードが角の点の正の側にあれば +1、負の側なら −1 で、
//! J はそのノードの展開の角の点での値。不均質な角の中性子束 Φ は角を囲むノードで共通とし、各ノードの均質な値は
//! Φ / CDF とする（燃料の CDF で割った反射体の角は `corner_factors`。ADR-0125）。ノードの展開は角の値について線形なので、Φ の連立一次方程式を組んで解く。
//!
//! セルの出力は Σ_g Σf,g φ̄_g（セルの均質な中性子束の平均）f_g（形状関数）。

use std::collections::HashMap;
use std::rc::Rc;

use xlc_format::fewgroup::{FewGroup, Pins};
use xlc_numeric::quadrature::gauss_legendre;

use crate::nem::{Boundary, Problem, Solution};
use xlc_numeric::nem::is_fissile;
use xlc_numeric::pins::corners::{CornerNode, corner_balance, fuel_harmonic_mean, nodes_at_point};
use xlc_numeric::pins::node::{Conditions, Constants, NodeFlux, NodeModel};
use xlc_numeric::pins::{CELL_POINTS, FACE_POINTS};

/// [ノード] のセルの出力 [行（-y から）][列（-x から）]。ノードの並びは NEM と同じ（x が最も速く、次に y、z）。
/// 核分裂のない軸方向の層（上下の反射体）は出力が 0 なので、空の表にする。核分裂のある層では、どのノードも形状関数を
/// 持つこと（径方向の反射体のノードは、形状関数 0 と角の不連続因子を持つ。ADR-0125）
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
        let constants = layer_constants(problem, k)?;
        let (models, conditions) = layer(problem, solution, k, &constants, &face_rule)?;
        let cdfs = corner_factors(problem, k, &constants)?;
        let corners = corner_fluxes(problem, k, &constants, &models, &conditions, &cdfs)?;
        for j in 0..ny {
            for i in 0..nx {
                let local = i + nx * j;
                let n = local + nx * ny * k;
                let (m, pins) = constants[local];
                let mut c = conditions[local].clone();
                c.corners = corner_values(&corners, nx, i, j, &cdfs[local]);
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

/// 再構成のノードの展開の定数
fn expansion_constants(m: &FewGroup) -> Constants<'_> {
    Constants {
        total: &m.total,
        scatter: &m.scatter,
        nu_fission: &m.nu_fission,
        chi: &m.chi,
        diffusion: &m.diffusion,
    }
}

/// ノード n の少数群定数と形状関数
fn node_constants(problem: &Problem, n: usize) -> Result<(&FewGroup, &Pins), String> {
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
}

/// 層 k の [ノード（層の中で i + nx j）] の少数群定数と形状関数
fn layer_constants(problem: &Problem, k: usize) -> Result<Vec<(&FewGroup, &Pins)>, String> {
    let (nx, ny) = (problem.mesh.dx.len(), problem.mesh.dy.len());
    (0..nx * ny)
        .map(|local| node_constants(problem, local + nx * ny * k))
        .collect()
}

/// 層 k の [ノード（層の中で i + nx j）] のモデルと、角の値を 0 にした条件。モデルは物質とノードの大きさで決まるので、
/// 同じ組のノードで共有する
fn layer(
    problem: &Problem,
    solution: &Solution,
    k: usize,
    constants: &[(&FewGroup, &Pins)],
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
            let (m, pins) = constants[i + nx * j];
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
                        NodeModel::new(&expansion_constants(m), solution.k, half, rule)
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
/// 厚さで重みを付けて足す（核分裂のない層の空の表は足さない）。集合体を分けた部分ノードが分け合ったセル（形状関数の
/// `halves`。ADR-0108）は、面積で重み付けして1つに戻す。すべて 0 の行と列（集合体の間の余白など）を除き、核分裂の
/// あるセルの平均で 1 に規格化する（`lattice` の `pin-powers` と同じ規格化）
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
    // [行（-y から）][列] の層の和
    let mut total: Option<Vec<Vec<f64>>> = None;
    for k in 0..nz {
        let layer = &maps[nx * ny * k..nx * ny * (k + 1)];
        if layer.iter().all(Vec::is_empty) {
            continue;
        }
        let grid = layer_grid(problem, k, layer)?;
        let sum = total.get_or_insert_with(|| grid.iter().map(|r| vec![0.0; r.len()]).collect());
        if sum.len() != grid.len() || sum.iter().zip(&grid).any(|(a, b)| a.len() != b.len()) {
            return Err(format!(
                "layer {k}: its pin map differs in shape from the other layers"
            ));
        }
        for (s, v) in sum.iter_mut().flatten().zip(grid.iter().flatten()) {
            *s += mesh.dz[k] * v;
        }
    }
    let total = total.ok_or("no cell has fission")?;
    // 上の行から
    let rows: Vec<Vec<f64>> = total.into_iter().rev().collect();
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

/// セルの大きさの並びが、丸めの程度（全体の 1e-9）で等しいか。切ったセルの幅は、回した集合体と回していない集合体で
/// 最後の桁が違いうる（座標の引き算の丸め）
fn same_sizes(a: &[f64], b: &[f64]) -> bool {
    let scale: f64 = b.iter().sum();
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= 1e-9 * scale)
}

/// 層 k の炉心全体のセルの出力 [行（-y から）][列]。分け合ったセルは面積で重み付けして1つに戻す
fn layer_grid(
    problem: &Problem,
    k: usize,
    layer: &[Vec<Vec<f64>>],
) -> Result<Vec<Vec<f64>>, String> {
    let mesh = &problem.mesh;
    let (nx, ny) = (mesh.dx.len(), mesh.dy.len());
    let pins = |i: usize, j: usize| -> Result<&Pins, String> {
        node_constants(problem, i + nx * (j + ny * k)).map(|(_, p)| p)
    };
    // ノードの列ごとのセルの幅と、ノードの行ごとのセルの高さは、その列（行）のどのノードでも同じ
    let mut widths: Vec<&[f64]> = Vec::with_capacity(nx);
    for i in 0..nx {
        let first = &pins(i, 0)?.dx;
        for j in 0..ny {
            if !same_sizes(&pins(i, j)?.dx, first) {
                return Err(format!(
                    "layer {k}: the nodes of column {i} have different cell widths"
                ));
            }
        }
        widths.push(first);
    }
    let mut heights: Vec<&[f64]> = Vec::with_capacity(ny);
    for j in 0..ny {
        let first = &pins(0, j)?.dy;
        for i in 0..nx {
            if !same_sizes(&pins(i, j)?.dy, first) {
                return Err(format!(
                    "layer {k}: the nodes of row {j} have different cell heights"
                ));
            }
        }
        heights.push(first);
    }
    // ノードの列 i と i + 1（行 j と j + 1）の間で分け合ったセルがあるか。両側の印がそろうこと
    let shared = |a: &Pins, a_side: usize, b: &Pins, b_side: usize| -> Result<bool, String> {
        if a.halves[a_side] != b.halves[b_side] {
            return Err(format!(
                "layer {k}: a split cell is marked on only one side"
            ));
        }
        Ok(a.halves[a_side])
    };
    let mut merge_x = Vec::with_capacity(nx.saturating_sub(1));
    for i in 0..nx.saturating_sub(1) {
        let first = shared(pins(i, 0)?, 1, pins(i + 1, 0)?, 0)?;
        for j in 1..ny {
            if shared(pins(i, j)?, 1, pins(i + 1, j)?, 0)? != first {
                return Err(format!(
                    "layer {k}: split cells differ along the column boundary {i}"
                ));
            }
        }
        merge_x.push(first);
    }
    let mut merge_y = Vec::with_capacity(ny.saturating_sub(1));
    for j in 0..ny.saturating_sub(1) {
        let first = shared(pins(0, j)?, 2, pins(0, j + 1)?, 3)?;
        for i in 1..nx {
            if shared(pins(i, j)?, 2, pins(i, j + 1)?, 3)? != first {
                return Err(format!(
                    "layer {k}: split cells differ along the row boundary {j}"
                ));
            }
        }
        merge_y.push(first);
    }
    // 分ける前の全体の表と、列の幅、行の高さ
    let mut grid: Vec<Vec<f64>> = Vec::new();
    for j in 0..ny {
        for r in 0..heights[j].len() {
            let mut row = Vec::new();
            for i in 0..nx {
                let map = &layer[i + nx * j];
                let cells = map
                    .get(r)
                    .filter(|cells| cells.len() == widths[i].len())
                    .ok_or_else(|| {
                        format!("layer {k}: node ({i}, {j}) has a pin map of a wrong size")
                    })?;
                row.extend_from_slice(cells);
            }
            grid.push(row);
        }
    }
    let concat = |sizes: &[&[f64]]| {
        sizes
            .iter()
            .flat_map(|s| s.iter().copied())
            .collect::<Vec<f64>>()
    };
    // 列ごと、行ごとに、分け合ったセルの対（境の両側の番号）
    let pairs = |sizes: &[&[f64]], merge: &[bool]| -> Vec<usize> {
        let mut at = 0;
        let mut out = Vec::new();
        for (n, s) in sizes.iter().enumerate() {
            at += s.len();
            if merge.get(n).copied().unwrap_or(false) {
                // 境の左（下）のセルは at - 1、右（上）は at
                out.push(at - 1);
            }
        }
        out
    };
    let column_pairs = pairs(&widths, &merge_x);
    let row_pairs = pairs(&heights, &merge_y);
    let merge = |values: &[f64], sizes: &[f64], pairs: &[usize]| -> Vec<f64> {
        let mut out = Vec::with_capacity(values.len());
        let mut c = 0;
        while c < values.len() {
            if pairs.contains(&c) {
                let (w0, w1) = (sizes[c], sizes[c + 1]);
                out.push((values[c] * w0 + values[c + 1] * w1) / (w0 + w1));
                c += 2;
            } else {
                out.push(values[c]);
                c += 1;
            }
        }
        out
    };
    let (all_widths, all_heights) = (concat(&widths), concat(&heights));
    let merged_columns: Vec<Vec<f64>> = grid
        .iter()
        .map(|row| merge(row, &all_widths, &column_pairs))
        .collect();
    let columns = merged_columns.first().map_or(0, Vec::len);
    let transposed: Vec<Vec<f64>> = (0..columns)
        .map(|c| merged_columns.iter().map(|r| r[c]).collect())
        .collect();
    let merged: Vec<Vec<f64>> = transposed
        .iter()
        .map(|col| merge(col, &all_heights, &row_pairs))
        .collect();
    let rows = merged.first().map_or(0, Vec::len);
    Ok((0..rows)
        .map(|r| merged.iter().map(|col| col[r]).collect())
        .collect())
}

/// ノード (i, j) の 4 つの角の均質な中性子束 Φ / CDF
fn corner_values(
    corners: &[[f64; 2]],
    nx: usize,
    i: usize,
    j: usize,
    cdf: &[[f64; 2]; 4],
) -> [[f64; 2]; 4] {
    let at = |ci: usize, cj: usize| corners[ci + (nx + 1) * cj];
    let points = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)];
    std::array::from_fn(|c| std::array::from_fn(|g| points[c][g] / cdf[c][g]))
}

/// 層 k の [ノード（層の中で i + nx j）][角][群] の CDF。燃料の CDF で割った角（`normalized_corners`。ADR-0125）は、
/// C* × 角の点に接する燃料（核分裂のある）のノードの、その点の CDF の調和平均とする。印の付いた角の点に燃料の
/// ノードがなければエラー
fn corner_factors(
    problem: &Problem,
    k: usize,
    nodes: &[(&FewGroup, &Pins)],
) -> Result<Vec<[[f64; 2]; 4]>, String> {
    let (nx, ny) = (problem.mesh.dx.len(), problem.mesh.dy.len());
    // 燃料は核分裂のあるノード（lattice が C* を割るときと同じ判定。ADR-0125）
    let fissile: Vec<bool> = nodes
        .iter()
        .map(|(m, _)| is_fissile(&m.nu_fission))
        .collect::<Result<_, _>>()
        .map_err(|e| e.to_string())?;
    (0..nx * ny)
        .map(|local| {
            let (i, j) = (local % nx, local / nx);
            let pins = nodes[local].1;
            let mut cdf = [[0.0; 2]; 4];
            for (c, out) in cdf.iter_mut().enumerate() {
                let own = &pins.corners[c];
                if !pins.normalized_corners[c] {
                    *out = [own[0], own[1]];
                    continue;
                }
                let (ci, cj) = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)][c];
                let fuel: Vec<(usize, usize)> = nodes_at_point(nx, ny, ci, cj)
                    .into_iter()
                    .filter(|&(m, _)| fissile[m])
                    .collect();
                for (g, o) in out.iter_mut().enumerate() {
                    let factors: Vec<f64> =
                        fuel.iter().map(|&(m, d)| nodes[m].1.corners[d][g]).collect();
                    let harmonic = fuel_harmonic_mean(&factors).ok_or_else(|| {
                        format!(
                            "node {}: corner {c} is normalized by the fuel factor, but no fuel node \
                             touches the corner point",
                            local + nx * ny * k
                        )
                    })?;
                    *o = own[g] * harmonic;
                }
            }
            Ok(cdf)
        })
        .collect()
}

/// 層 k の [角の点（(nx + 1) × (ny + 1)、x が速い）][群] の不均質な中性子束 Φ。各ノードの均質な角の値は Φ / CDF
fn corner_fluxes(
    problem: &Problem,
    k: usize,
    constants: &[(&FewGroup, &Pins)],
    models: &[Rc<NodeModel>],
    conditions: &[Conditions],
    cdfs: &[[[f64; 2]; 4]],
) -> Result<Vec<[f64; 2]>, String> {
    let (nx, ny) = (problem.mesh.dx.len(), problem.mesh.dy.len());
    let mut nodes = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let local = i + nx * j;
            let n = local + nx * ny * k;
            let (m, _) = constants[local];
            // 外周の真空の辺を通る流れ r f φ_c。面の条件（ADR-0101 の J = r f φ）を均質な角の値 φ_c = Φ / CDF に
            // 課す（ADR-0131）。不均質な Φ に r を掛けないのは、均質な角の値を均質な量だけで決めるため。反射と
            // 中性子束ゼロは 0
            let on_boundary = [i == 0, i == nx - 1, j == 0, j == ny - 1];
            let mut outflow = [None; 4];
            for (face, out) in outflow.iter_mut().enumerate() {
                if !on_boundary[face] {
                    continue;
                }
                // 辺（左、右、上、下）の番号は、面（-x、+x、-y、+y）の順に 0、1、3、2
                let side = [0, 1, 3, 2][face];
                let ratio = match problem.boundaries[face] {
                    Boundary::Reflective | Boundary::ZeroFlux => continue,
                    Boundary::Marshak => [0.5; 2],
                    Boundary::Vacuum => {
                        let r = m.vacuum_ratio[side]
                            .as_ref()
                            .ok_or_else(|| format!("node {n}: vacuum side without a ratio"))?;
                        [r[0], r[1]]
                    }
                };
                // 外周の辺に f* がないことは、NEM の解で確かめてある（`check_normalized_adfs`）
                *out = Some(std::array::from_fn(|g| ratio[g] * m.adf[side][g]));
            }
            nodes.push(CornerNode {
                model: &models[local],
                conditions: &conditions[local],
                diffusion: [m.diffusion[0], m.diffusion[1]],
                // 均質な角の値は Φ / CDF
                weights: std::array::from_fn(|c| std::array::from_fn(|g| 1.0 / cdfs[local][c][g])),
                outflow,
            });
        }
    }
    // 中性子束ゼロの外周の角は Φ = 0
    let b = &problem.boundaries;
    let zero: Vec<bool> = (0..=ny)
        .flat_map(|cj| {
            (0..=nx).map(move |ci| {
                (ci == 0 && b[0] == Boundary::ZeroFlux)
                    || (ci == nx && b[1] == Boundary::ZeroFlux)
                    || (cj == 0 && b[2] == Boundary::ZeroFlux)
                    || (cj == ny && b[3] == Boundary::ZeroFlux)
            })
        })
        .collect();
    corner_balance(nx, ny, &nodes, &zero)
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
