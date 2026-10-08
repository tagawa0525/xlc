//! MOC の解から、少数群定数（`xlc_format::fewgroup::FewGroup`）を作る。

use std::ops::Range;

use xlc_format::fewgroup::FewGroup;
use xlc_format::group::GroupStructure;
use xlc_numeric::nem::{Edge, OneNode, OneNodeError, surface_fluxes};

use crate::b1::{B1Error, Medium};
use crate::cmfd::{horizontal_face, vertical_face};
use crate::geometry::Lattice;
use crate::material::Material;
use crate::moc::Solution;
use crate::tracking::{Boundaries, Boundary};

/// 少数群への縮約のエラー。
#[derive(Debug, thiserror::Error)]
pub enum HomogenizeError {
    #[error("{materials} materials for {regions} regions")]
    RegionCount { materials: usize, regions: usize },
    #[error("few-group bounds must be bounds of the library groups")]
    NotLibraryBounds,
    #[error("few-group bounds must be descending and span the library groups")]
    Span,
    #[error("a few group has no flux")]
    NoFlux,
    #[error("library group {0} has no flux")]
    NoFineFlux(usize),
    #[error("discontinuity factors cannot be made with an albedo side")]
    Albedo,
    #[error(
        "region {region} has no P1 scattering moments; diffusion coefficients need a library with them"
    )]
    NoP1 { region: usize },
    #[error("few group {group} has transport cross section {transport}")]
    Transport { group: usize, transport: f64 },
    #[error(
        "few group {group} has no net current in the lattice; current-weighted diffusion coefficients need one"
    )]
    NoNetCurrent { group: usize },
    #[error("B1 critical spectrum")]
    B1(#[from] B1Error),
    #[error("one-node problem for the discontinuity factors")]
    OneNode(#[from] OneNodeError),
    #[error(
        "vacuum side {side}, few group {group}: net current / surface flux is {ratio}, outside (0, 1]"
    )]
    VacuumRatio {
        side: usize,
        group: usize,
        ratio: f64,
    },
    #[error("nodes: {0}")]
    Nodes(String),
}

/// 格子のセルの長方形の範囲で、1つの少数群定数にまとめるもの（ADR-0095）。
#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    /// セルの列（左から、0 始まり）の範囲
    pub columns: Range<usize>,
    /// セルの行（上から、0 始まり）の範囲
    pub rows: Range<usize>,
    /// 反射体のノード: 反射体でないノード（燃料）と接する辺の ADF を、その燃料のノードの同じ面の ADF で割った
    /// f* = f_反射体 / f_燃料 で持つ
    pub reflector: bool,
}

/// 辺（左、右、上、下）と、1ノードの問題の面（-x, +x, -y, +y）の対応: 面 i は辺 `FACES[i]`
const FACES: [usize; 4] = [0, 1, 3, 2];

/// 辺の向かいの辺（左と右、上と下）
const OPPOSITE: [usize; 4] = [1, 0, 3, 2];

/// MOC の解を使い、領域・群ごとの反応率を保存するようにセル全体で均質化し、
/// `coarse_bounds`（`fine` の境界の一部、降順）の少数群に縮約する。
/// 拡散係数は `diffusion_coefficients` で作る（ADR-0104）。格子全体を1つのノード（集合体）として均質化する。
///
/// 不連続因子（ADF）は、辺の上の不均質な中性子束（MOC）を、同じ辺の均質な中性子束で割る（一般化等価理論、
/// ADR-0103）。均質な中性子束は、均質化した定数、MOC の k とノード平均の中性子束、辺ごとの正味の中性子流で
/// 解いた1ノードの問題（`xlc_numeric::nem::surface_fluxes`。core と同じ核で、核分裂のあるノードは NEM、ない
/// ノードは ANM。ADR-0102）の面の値とする。外周がすべて
/// 鏡面反射なら正味の中性子流は 0 で、均質な解は平らになり、体積平均の中性子束で割るのと同じになる。
/// アルベドの辺は core に対応する境界がないのでエラーにする。
///
/// 真空の辺には、MOC の辺の外向きの正味の中性子流 / 辺の上の不均質な中性子束（どちらも辺にわたる平均で、
/// 少数群の中の細かい群の和）を持たせる。core の真空の面の条件 J = r f φ に使う（ADR-0101）。入射がないので
/// 0 < r ≤ 1 で、外れればエラーにする。
pub fn homogenize(
    solution: &Solution,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
    boundaries: &Boundaries,
) -> Result<FewGroup, HomogenizeError> {
    let sides = sides_of(boundaries);
    if sides.iter().any(|b| matches!(b, Boundary::Albedo(_))) {
        return Err(HomogenizeError::Albedo);
    }
    check_regions(solution, materials)?;
    let group_of = coarse_groups(fine, coarse_bounds)?;
    let coarse = coarse_bounds.len() - 1;
    let regions: Vec<usize> = (0..materials.len()).collect();
    let (mut few, average) = constants(
        solution,
        materials,
        &regions,
        &group_of,
        coarse_bounds,
        reflective(boundaries),
    )?;
    let collapse = |values: &[Vec<f64>]| collapse(values, &group_of, coarse);
    let (heterogeneous, currents) = (collapse(&solution.sides), collapse(&solution.currents));
    few.vacuum_ratio = vacuum_ratios(&sides, &heterogeneous, &currents)?;
    let edges = FACES.map(|side| match sides[side] {
        Boundary::Vacuum => Edge::Vacuum,
        _ => Edge::Reflective,
    });
    few.adf = discontinuity_factors(
        &few,
        solution.k,
        &average,
        solution.size,
        &heterogeneous,
        &currents,
        edges,
    )?;
    Ok(few)
}

/// 格子をノード `nodes`（格子のセルを重ならずに覆う長方形）に分け、ノードごとに `homogenize` と同じ形で
/// 少数群定数を作る（ADR-0095）。返す定数はノードの順。
///
/// - 断面積は、ノードの領域だけで反応率を保存する。拡散係数は、ノードの4辺がすべて外周の鏡面反射の辺
///   （格子全体が1つのノード）なら B1、ほかは outflow（ADR-0098）。`k_inf` は格子の k
/// - 辺の不均質な中性子束と正味の中性子流は、辺の上の粗メッシュの面の値（MOC の `face_fluxes`、
///   `face_currents`）を、面の長さで重みづけして平均する。内側の辺も同じ
/// - ADF は1ノードの問題（ADR-0094）で作る。外周の辺は格子の境界（鏡面反射か真空）、内側の辺は隣の
///   ノード（その幅と、参照の横方向の漏れ。core と同じ2次式）とする。内側の辺は、向かいにちょうど1つの
///   ノードがあり、辺の長さがそろうこと
/// - 外周の真空の辺には J/φ_het を持たせる（ADR-0101）
/// - 反射体のノード（`Node::reflector`）の、反射体でないノードと接する辺の ADF は、その隣のノードの同じ面の
///   ADF で割った f* とし、`normalized_adf` に印を付ける。反射体のノードは、反射体でないノードと少なくとも
///   1辺で接すること
pub fn homogenize_nodes(
    solution: &Solution,
    lattice: &Lattice,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
    boundaries: &Boundaries,
    nodes: &[Node],
) -> Result<Vec<FewGroup>, HomogenizeError> {
    let sides = sides_of(boundaries);
    if sides.iter().any(|b| matches!(b, Boundary::Albedo(_))) {
        return Err(HomogenizeError::Albedo);
    }
    check_regions(solution, materials)?;
    let group_of = coarse_groups(fine, coarse_bounds)?;
    let coarse = coarse_bounds.len() - 1;
    let neighbors = node_neighbors(lattice, nodes)?;
    let (widths, heights) = (lattice.widths(), lattice.heights());
    let (columns, rows) = (lattice.columns(), lattice.rows());

    // ノードごとの定数、ノード平均の中性子束、辺（左、右、上、下）の不均質な中性子束と外向きの正味の中性子流
    struct Homogenized {
        few: FewGroup,
        average: Vec<f64>,
        size: [f64; 2],
        heterogeneous: Vec<Vec<f64>>,
        currents: Vec<Vec<f64>>,
    }
    let mut homogenized = Vec::with_capacity(nodes.len());
    for (n, node) in nodes.iter().enumerate() {
        let regions: Vec<usize> = node
            .rows
            .clone()
            .flat_map(|r| node.columns.clone().map(move |c| r * columns + c))
            .flat_map(|cell| lattice.cell_regions(cell))
            .collect();
        let whole = neighbors[n].iter().all(Option::is_none);
        let (mut few, average) = constants(
            solution,
            materials,
            &regions,
            &group_of,
            coarse_bounds,
            whole && reflective(boundaries),
        )?;
        // 辺の上の面（面の番号、長さ、面の正の向きが外向きなら 1、逆なら -1）
        let (c0, c1, r0, r1) = (
            node.columns.start,
            node.columns.end,
            node.rows.start,
            node.rows.end,
        );
        let side_faces: [Vec<(usize, f64, f64)>; 4] = [
            node.rows
                .clone()
                .map(|r| (vertical_face(columns, r, c0), heights[r], -1.0))
                .collect(),
            node.rows
                .clone()
                .map(|r| (vertical_face(columns, r, c1), heights[r], 1.0))
                .collect(),
            node.columns
                .clone()
                .map(|c| (horizontal_face(rows, columns, r0, c), widths[c], -1.0))
                .collect(),
            node.columns
                .clone()
                .map(|c| (horizontal_face(rows, columns, r1, c), widths[c], 1.0))
                .collect(),
        ];
        let mean = |faces: &[(usize, f64, f64)], values: &[Vec<f64>], signed: bool| {
            let length: f64 = faces.iter().map(|f| f.1).sum();
            let mut sum = vec![0.0; group_of.len()];
            for &(face, l, outward) in faces {
                let sign = if signed { outward } else { 1.0 };
                for (s, v) in sum.iter_mut().zip(&values[face]) {
                    *s += sign * l * v;
                }
            }
            sum.iter().map(|s| s / length).collect::<Vec<f64>>()
        };
        let fine_heterogeneous: Vec<Vec<f64>> = side_faces
            .iter()
            .map(|faces| mean(faces, &solution.face_fluxes, false))
            .collect();
        let fine_currents: Vec<Vec<f64>> = side_faces
            .iter()
            .map(|faces| mean(faces, &solution.face_currents, true))
            .collect();
        let heterogeneous = collapse(&fine_heterogeneous, &group_of, coarse);
        let currents = collapse(&fine_currents, &group_of, coarse);
        // 外周の辺だけが格子の境界条件を持つ
        let node_sides: [Boundary; 4] = std::array::from_fn(|side| match neighbors[n][side] {
            None => sides[side],
            Some(_) => Boundary::Reflective,
        });
        few.vacuum_ratio = vacuum_ratios(&node_sides, &heterogeneous, &currents)?;
        tracing::info!(node = n, average_flux = ?average, "node");
        homogenized.push(Homogenized {
            few,
            average,
            size: [widths[c0..c1].iter().sum(), heights[r0..r1].iter().sum()],
            heterogeneous,
            currents,
        });
    }
    // 横方向の漏れ [ノード][方向（x, y）][少数群]: x の方向の横方向の漏れは y の正味の漏れ（上と下の辺の
    // 外向きの正味の中性子流の和 / 高さ）、y の方向は x の正味の漏れ
    let transverse: Vec<[Vec<f64>; 2]> = homogenized
        .iter()
        .map(|h| {
            let net = |a: usize, b: usize, width: f64| -> Vec<f64> {
                (0..coarse)
                    .map(|g| (h.currents[a][g] + h.currents[b][g]) / width)
                    .collect()
            };
            [net(2, 3, h.size[1]), net(0, 1, h.size[0])]
        })
        .collect();
    let mut adfs = Vec::with_capacity(nodes.len());
    for (n, h) in homogenized.iter().enumerate() {
        let edges: [Edge; 4] = FACES.map(|side| {
            // 左右の辺は x の方向
            let d = usize::from(side >= 2);
            match neighbors[n][side] {
                Some(m) => Edge::Neighbor {
                    width: homogenized[m].size[d],
                    leakage: &transverse[m][d],
                },
                None if sides[side] == Boundary::Vacuum => Edge::Vacuum,
                None => Edge::Reflective,
            }
        });
        adfs.push(discontinuity_factors(
            &h.few,
            solution.k,
            &h.average,
            h.size,
            &h.heterogeneous,
            &h.currents,
            edges,
        )?);
    }
    let mut result = Vec::with_capacity(nodes.len());
    for (n, (node, h)) in nodes.iter().zip(homogenized).enumerate() {
        let mut few = h.few;
        few.adf = adfs[n].clone();
        if node.reflector {
            for side in 0..4 {
                if let Some(m) = neighbors[n][side]
                    && !nodes[m].reflector
                {
                    for (f, fuel) in few.adf[side].iter_mut().zip(&adfs[m][OPPOSITE[side]]) {
                        *f /= fuel;
                    }
                    few.normalized_adf[side] = true;
                }
            }
            if !few.normalized_adf.contains(&true) {
                return Err(HomogenizeError::Nodes(format!(
                    "reflector node {n} touches no fuel (non-reflector) node"
                )));
            }
        }
        result.push(few);
    }
    Ok(result)
}

/// [ノード][辺（左、右、上、下）] の、辺の向かいのノード（外周の辺は `None`）。ノードが格子のセルを重ならずに
/// 覆い、内側の辺の向かいがちょうど1つのノードで、長さがそろうことを確かめる。
fn node_neighbors(
    lattice: &Lattice,
    nodes: &[Node],
) -> Result<Vec<[Option<usize>; 4]>, HomogenizeError> {
    let error = |message: String| Err(HomogenizeError::Nodes(message));
    if nodes.is_empty() {
        return error("no nodes".into());
    }
    let (columns, rows) = (lattice.columns(), lattice.rows());
    let mut owner = vec![None; columns * rows];
    for (n, node) in nodes.iter().enumerate() {
        if node.columns.is_empty() || node.rows.is_empty() {
            return error(format!("node {n} has no cells"));
        }
        if node.columns.end > columns || node.rows.end > rows {
            return error(format!(
                "node {n} (columns {:?}, rows {:?}) is outside the {rows} x {columns} lattice",
                node.columns, node.rows
            ));
        }
        for r in node.rows.clone() {
            for c in node.columns.clone() {
                if let Some(other) = owner[r * columns + c].replace(n) {
                    return error(format!(
                        "nodes {other} and {n} overlap at row {r}, column {c}"
                    ));
                }
            }
        }
    }
    if let Some(cell) = owner.iter().position(Option::is_none) {
        return error(format!(
            "no node covers the cell at row {}, column {}",
            cell / columns,
            cell % columns
        ));
    }
    nodes
        .iter()
        .enumerate()
        .map(|(n, node)| {
            let (c0, c1, r0, r1) = (
                node.columns.start,
                node.columns.end,
                node.rows.start,
                node.rows.end,
            );
            // 辺の向かいのセル（外周なら空）
            let across: [Vec<usize>; 4] = [
                if c0 == 0 {
                    vec![]
                } else {
                    node.rows.clone().map(|r| r * columns + c0 - 1).collect()
                },
                if c1 == columns {
                    vec![]
                } else {
                    node.rows.clone().map(|r| r * columns + c1).collect()
                },
                if r0 == 0 {
                    vec![]
                } else {
                    node.columns
                        .clone()
                        .map(|c| (r0 - 1) * columns + c)
                        .collect()
                },
                if r1 == rows {
                    vec![]
                } else {
                    node.columns.clone().map(|c| r1 * columns + c).collect()
                },
            ];
            let mut neighbors = [None; 4];
            for (side, cells) in across.iter().enumerate() {
                let Some(&first) = cells.first() else {
                    continue;
                };
                let m = owner[first].unwrap_or(n);
                let same = if side < 2 {
                    nodes[m].rows == node.rows
                } else {
                    nodes[m].columns == node.columns
                };
                if !same || cells.iter().any(|&cell| owner[cell] != Some(m)) {
                    return Err(HomogenizeError::Nodes(format!(
                        "side {side} of node {n} does not face exactly one node with the same side"
                    )));
                }
                neighbors[side] = Some(m);
            }
            Ok(neighbors)
        })
        .collect()
}

/// 外周の辺（左、右、上、下）の境界条件
fn sides_of(boundaries: &Boundaries) -> [Boundary; 4] {
    [
        boundaries.left,
        boundaries.right,
        boundaries.top,
        boundaries.bottom,
    ]
}

/// 領域 `regions` を中性子束と体積で重みづけして均質化した少数群定数（ADF と真空の辺の比は空、
/// `normalized_adf` はすべて false）と、ノード平均の少数群の中性子束。`b1` なら拡散係数を B1 で、
/// でなければ outflow で作る。
fn constants(
    solution: &Solution,
    materials: &[&Material],
    regions: &[usize],
    group_of: &[usize],
    coarse_bounds: &[f64],
    b1: bool,
) -> Result<(FewGroup, Vec<f64>), HomogenizeError> {
    let coarse = coarse_bounds.len() - 1;
    let mut flux = vec![0.0; coarse];
    let mut total = vec![0.0; coarse];
    let mut nu_fission = vec![0.0; coarse];
    let mut fission = vec![0.0; coarse];
    let mut scatter = vec![vec![0.0; coarse]; coarse];
    let mut chi = vec![0.0; coarse];
    for &region in regions {
        let (m, phi, volume) = (
            materials[region],
            &solution.flux[region],
            solution.volumes[region],
        );
        for (g, &big) in group_of.iter().enumerate() {
            let w = volume * phi[g];
            flux[big] += w;
            total[big] += w * m.total[g];
            nu_fission[big] += w * m.nu_fission[g];
            fission[big] += w * m.fission[g];
            for (to, s) in m.scatter[g].iter().enumerate() {
                scatter[big][group_of[to]] += w * s;
            }
        }
        for term in &m.fission_terms {
            let rate: f64 = term
                .nu_fission
                .iter()
                .zip(phi)
                .map(|(n, f)| volume * n * f)
                .sum();
            for (g, c) in term.chi.iter().enumerate() {
                chi[group_of[g]] += c * rate;
            }
        }
    }
    if flux.iter().any(|&f| f <= 0.0) {
        return Err(HomogenizeError::NoFlux);
    }
    let per_flux =
        |values: Vec<f64>| -> Vec<f64> { values.iter().zip(&flux).map(|(v, f)| v / f).collect() };
    let total = per_flux(total);
    let scatter: Vec<Vec<f64>> = scatter
        .into_iter()
        .zip(&flux)
        .map(|(row, f)| row.iter().map(|v| v / f).collect())
        .collect();
    let produced: f64 = chi.iter().sum();
    let volume: f64 = regions.iter().map(|&r| solution.volumes[r]).sum();
    let average = flux.iter().map(|f| f / volume).collect();
    let few = FewGroup {
        bounds: coarse_bounds.to_vec(),
        k_inf: solution.k,
        absorption: total
            .iter()
            .zip(&scatter)
            .map(|(t, row)| t - row.iter().sum::<f64>())
            .collect(),
        diffusion: diffusion(solution, materials, regions, group_of, coarse, b1)?,
        total,
        nu_fission: per_flux(nu_fission),
        fission: per_flux(fission),
        chi: if produced > 0.0 {
            chi.iter().map(|c| c / produced).collect()
        } else {
            chi
        },
        scatter,
        adf: Vec::new(),
        vacuum_ratio: Vec::new(),
        normalized_adf: vec![false; 4],
    };
    Ok((few, average))
}

/// [辺] ごとの細かい群の値を少数群に縮約する（細かい群の和。辺の平均の中性子束と正味の中性子流は群について
/// 足せる）
fn collapse(values: &[Vec<f64>], group_of: &[usize], coarse: usize) -> Vec<Vec<f64>> {
    values
        .iter()
        .map(|side| {
            let mut sum = vec![0.0; coarse];
            for (g, v) in side.iter().enumerate() {
                sum[group_of[g]] += v;
            }
            sum
        })
        .collect()
}

/// [辺（左、右、上、下）][少数群] の ADF: 辺の不均質な中性子束 `heterogeneous` を、1ノードの問題（ADR-0094）
/// の均質な面の中性子束で割る。`currents` は辺の外向きの正味の中性子流、`edges` は1ノードの問題の面
/// （-x, +x, -y, +y）の外側の扱い。
fn discontinuity_factors(
    few: &FewGroup,
    k: f64,
    average: &[f64],
    size: [f64; 2],
    heterogeneous: &[Vec<f64>],
    currents: &[Vec<f64>],
    edges: [Edge; 4],
) -> Result<Vec<Vec<f64>>, HomogenizeError> {
    tracing::debug!(
        diffusion = ?few.diffusion,
        total = ?few.total,
        scatter = ?few.scatter,
        nu_fission = ?few.nu_fission,
        k,
        size = ?size,
        "one-node constants"
    );
    for (side, het) in heterogeneous.iter().enumerate() {
        tracing::debug!(
            side,
            heterogeneous = ?het,
            outward_current = ?currents[side],
            average = ?average,
            "one-node inputs"
        );
    }
    let surface = surface_fluxes(&OneNode {
        diffusion: &few.diffusion,
        total: &few.total,
        scatter: &few.scatter,
        nu_fission: &few.nu_fission,
        chi: &few.chi,
        k,
        flux: average,
        widths: size,
        currents: FACES.map(|side| currents[side].as_slice()),
        edges,
    })?;
    for (side, &face) in FACES.iter().enumerate() {
        tracing::debug!(side, homogeneous = ?surface[face], "one-node face flux");
    }
    Ok((0..4)
        .map(|side| {
            heterogeneous[side]
                .iter()
                .zip(&surface[FACES[side]])
                .map(|(h, s)| h / s)
                .collect()
        })
        .collect())
}

/// [辺] の真空の辺の比 J/φ_het（[少数群]）。真空でない辺は `None`。
fn vacuum_ratios(
    sides: &[Boundary; 4],
    heterogeneous: &[Vec<f64>],
    currents: &[Vec<f64>],
) -> Result<Vec<Option<Vec<f64>>>, HomogenizeError> {
    (0..4)
        .map(|side| {
            if sides[side] != Boundary::Vacuum {
                return Ok(None);
            }
            heterogeneous[side]
                .iter()
                .zip(&currents[side])
                .enumerate()
                .map(|(group, (phi, j))| {
                    let ratio = j / phi;
                    if ratio > 0.0 && ratio <= 1.0 {
                        Ok(ratio)
                    } else {
                        Err(HomogenizeError::VacuumRatio { side, group, ratio })
                    }
                })
                .collect::<Result<Vec<f64>, _>>()
                .map(Some)
        })
        .collect()
}

fn reflective(boundaries: &Boundaries) -> bool {
    [
        boundaries.left,
        boundaries.right,
        boundaries.top,
        boundaries.bottom,
    ]
    .iter()
    .all(|&b| b == Boundary::Reflective)
}

fn check_regions(solution: &Solution, materials: &[&Material]) -> Result<(), HomogenizeError> {
    if materials.len() != solution.flux.len() {
        return Err(HomogenizeError::RegionCount {
            materials: materials.len(),
            regions: solution.flux.len(),
        });
    }
    Ok(())
}

/// [細かい群] の少数群の番号。少数群の境界は細かい群の境界の一部で、全体を覆う降順であること。
fn coarse_groups(
    fine: &GroupStructure,
    coarse_bounds: &[f64],
) -> Result<Vec<usize>, HomogenizeError> {
    let ng = fine.len();
    let positions = coarse_bounds
        .iter()
        .map(|b| fine.bounds().iter().position(|f| f == b))
        .collect::<Option<Vec<usize>>>()
        .ok_or(HomogenizeError::NotLibraryBounds)?;
    let covers = positions.first() == Some(&0) && positions.last() == Some(&ng);
    if !covers || positions.windows(2).any(|w| w[0] >= w[1]) {
        return Err(HomogenizeError::Span);
    }
    let mut group_of = vec![0; ng];
    for (big, w) in positions.windows(2).enumerate() {
        group_of[w[0]..w[1]].fill(big);
    }
    Ok(group_of)
}

/// 少数群の拡散係数（ADR-0104）。どちらもライブラリの P1 の散乱行列を使い、真空でない物質に P1 の
/// モーメントがなければエラーにする。
///
/// - 外周がすべて鏡面反射の格子: 細かい群のまま中性子束と体積で均質化し、B1 の臨界のバックリングでの
///   基本モードから D_g = J_g / (B φ_g) とし、B1 の中性子束で縮約する（D_G = Σ J_g / Σ B φ_g）。
///   核分裂のない媒質では臨界のバックリングが定まらないので、B² = 0 の中性子流の式（P1）を MOC の
///   中性子束で解き、MOC の中性子束で縮約する
/// - それ以外（真空かアルベドの辺がある）: outflow の輸送断面積 Σtr,g = Σt,g − Σ_g' Σs1(g→g') を、
///   領域と細かい群について、MOC の領域の中性子流 × 体積を格子の正味の中性子流の向きへ射影した重みで少数群に
///   平均し、D_G = 1/(3⟨Σtr⟩_G) とする（ADR-0104。真空の領域は Σtr = 0 で、中性子流は重みに入る）
pub fn diffusion_coefficients(
    solution: &Solution,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
    boundaries: &Boundaries,
) -> Result<Vec<f64>, HomogenizeError> {
    check_regions(solution, materials)?;
    let group_of = coarse_groups(fine, coarse_bounds)?;
    let regions: Vec<usize> = (0..materials.len()).collect();
    diffusion(
        solution,
        materials,
        &regions,
        &group_of,
        coarse_bounds.len() - 1,
        reflective(boundaries),
    )
}

/// 検査済みの領域と群の対応で、領域 `regions` の拡散係数を作る（`diffusion_coefficients`）。
fn diffusion(
    solution: &Solution,
    materials: &[&Material],
    regions: &[usize],
    group_of: &[usize],
    coarse: usize,
    reflective: bool,
) -> Result<Vec<f64>, HomogenizeError> {
    if let Some(&region) = regions
        .iter()
        .find(|&&r| !materials[r].is_void() && materials[r].anisotropic.is_empty())
    {
        return Err(HomogenizeError::NoP1 { region });
    }
    if reflective {
        b1_diffusion(solution, materials, regions, group_of, coarse)
    } else {
        outflow_diffusion(solution, materials, regions, group_of, coarse)
    }
}

fn b1_diffusion(
    solution: &Solution,
    materials: &[&Material],
    regions: &[usize],
    group_of: &[usize],
    coarse: usize,
) -> Result<Vec<f64>, HomogenizeError> {
    let ng = group_of.len();
    let mut flux = vec![0.0; ng];
    let mut medium = Medium {
        total: vec![0.0; ng],
        scatter: vec![vec![0.0; ng]; ng],
        scatter1: vec![vec![0.0; ng]; ng],
        nu_fission: vec![0.0; ng],
        chi: vec![0.0; ng],
    };
    for &region in regions {
        let (m, phi, volume) = (
            materials[region],
            &solution.flux[region],
            solution.volumes[region],
        );
        for g in 0..ng {
            let w = volume * phi[g];
            flux[g] += w;
            medium.total[g] += w * m.total[g];
            medium.nu_fission[g] += w * m.nu_fission[g];
            for (to, s) in m.scatter[g].iter().enumerate() {
                medium.scatter[g][to] += w * s;
            }
            // 真空（P1 のモーメントを持たない）は散乱しない
            if let Some(p1) = m.anisotropic.first() {
                for (to, s) in p1[g].iter().enumerate() {
                    medium.scatter1[g][to] += w * s;
                }
            }
        }
        for term in &m.fission_terms {
            let rate: f64 = term
                .nu_fission
                .iter()
                .zip(phi)
                .map(|(n, f)| volume * n * f)
                .sum();
            for (c, chi) in term.chi.iter().zip(medium.chi.iter_mut()) {
                *chi += c * rate;
            }
        }
    }
    if let Some(g) = flux.iter().position(|&f| f <= 0.0) {
        return Err(HomogenizeError::NoFineFlux(g));
    }
    for (g, f) in flux.iter().enumerate() {
        medium.total[g] /= f;
        medium.nu_fission[g] /= f;
        for v in medium.scatter[g]
            .iter_mut()
            .chain(medium.scatter1[g].iter_mut())
        {
            *v /= f;
        }
    }
    let produced: f64 = medium.chi.iter().sum();
    let (current, weight) = if produced > 0.0 {
        for c in &mut medium.chi {
            *c /= produced;
        }
        let mode = medium.critical()?;
        (mode.current, mode.flux)
    } else {
        // 核分裂のない媒質: B² = 0 の P1 の中性子流を、MOC の中性子束（領域で積分した値）で解く
        (medium.zero_buckling_current(&flux)?, flux)
    };
    let mut sum = vec![0.0; coarse];
    let mut total = vec![0.0; coarse];
    for g in 0..ng {
        sum[group_of[g]] += current[g];
        total[group_of[g]] += weight[g];
    }
    if total.iter().any(|&f| f <= 0.0) {
        return Err(HomogenizeError::NoFlux);
    }
    Ok(sum.iter().zip(&total).map(|(j, f)| j / f).collect())
}

/// 中性子流で重みづけした outflow の輸送断面積の拡散係数（ADR-0104）: Σtr,g = Σt,g − Σ_g' Σs1(g→g') を、
/// 領域の中性子流 × 体積をノードの正味の中性子流の向き ê_G（Σ_{g∈G} Σ_r V_r J⃗_r,g の向き）へ射影した符号つきの重み
/// w_r,g で少数群に平均し、D_G = 1/(3⟨Σtr⟩_G) とする（Fick の法則を中性子流について足した形）。真空の領域は Σtr = 0 で
/// 重みに入る。正味の中性子流が 0 なら向きが定まらず、⟨Σtr⟩ が正でなければ D が定まらないので、どちらもエラーにする。
fn outflow_diffusion(
    solution: &Solution,
    materials: &[&Material],
    regions: &[usize],
    group_of: &[usize],
    coarse: usize,
) -> Result<Vec<f64>, HomogenizeError> {
    let mut net = vec![[0.0; 2]; coarse];
    for &region in regions {
        let (current, volume) = (&solution.region_currents[region], solution.volumes[region]);
        for (g, &big) in group_of.iter().enumerate() {
            for (n, c) in net[big].iter_mut().zip(current) {
                *n += volume * c[g];
            }
        }
    }
    let mut transport = vec![0.0; coarse];
    let mut weight = vec![0.0; coarse];
    let mut magnitude = vec![0.0; coarse];
    for &region in regions {
        let (m, current, volume) = (
            materials[region],
            &solution.region_currents[region],
            solution.volumes[region],
        );
        for (g, &big) in group_of.iter().enumerate() {
            let length = net[big][0].hypot(net[big][1]);
            let w = volume * (current[0][g] * net[big][0] + current[1][g] * net[big][1]) / length;
            let out: f64 = m.anisotropic.first().map_or(0.0, |p1| p1[g].iter().sum());
            transport[big] += w * (m.total[g] - out);
            weight[big] += w;
            magnitude[big] += w.abs();
        }
    }
    (0..coarse)
        .map(|group| {
            let length = net[group][0].hypot(net[group][1]);
            if !(length > 0.0 && length.is_finite()) {
                return Err(HomogenizeError::NoNetCurrent { group });
            }
            // 打ち消しの程度: 1 なら領域と細かい群の中性子流がすべてノードの正味の中性子流の向き
            tracing::debug!(
                group,
                cancellation = weight[group] / magnitude[group],
                "current weights"
            );
            let transport = transport[group] / weight[group];
            if transport > 0.0 && transport.is_finite() {
                Ok(1.0 / (3.0 * transport))
            } else {
                Err(HomogenizeError::Transport { group, transport })
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Cell, Lattice, PinCell};
    use crate::material::FissionTerm;
    use crate::moc::{Options, solve};
    use crate::tracking::{Boundary, Tracking};

    fn two_group(total: [f64; 2], scatter: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Material {
        // P1 の散乱は P0 の 0.3 倍（前方散乱）
        let p1 = scatter
            .iter()
            .map(|r| r.iter().map(|s| 0.3 * s).collect())
            .collect();
        Material {
            total: total.to_vec(),
            scatter: scatter.iter().map(|r| r.to_vec()).collect(),
            anisotropic: vec![p1],
            fission: nu_fission.iter().map(|v| v / 2.5).collect(),
            nu_fission: nu_fission.to_vec(),
            fission_terms: if nu_fission.iter().any(|&v| v > 0.0) {
                vec![FissionTerm {
                    chi: vec![1.0, 0.0],
                    nu_fission: nu_fission.to_vec(),
                }]
            } else {
                vec![]
            },
        }
    }

    fn solved() -> (Solution, Vec<Material>, GroupStructure) {
        let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let moderator = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let materials = vec![fuel.clone(), moderator.clone(), moderator];
        let cell = PinCell::new(1.26, vec![0.4096, 0.55]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 16, 0.05).unwrap();
        let refs: Vec<&Material> = materials.iter().collect();
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        (
            solution,
            materials,
            GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
        )
    }

    #[test]
    fn homogeneous_cell_has_unit_adf() {
        let m = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let cell = PinCell::new(1.26, vec![0.4096]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.05).unwrap();
        let solution = solve(&tracking, &[&m, &m], &Options::default()).unwrap();
        let fine = GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap();
        let fg = homogenize(
            &solution,
            &[&m, &m],
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        for value in fg.adf.iter().flatten() {
            assert!((value - 1.0).abs() < 1e-8, "{:?}", fg.adf);
        }
    }

    #[test]
    fn reflective_adfs_divide_by_the_average_flux() {
        // 外周がすべて鏡面反射なら、1ノードの問題の均質な解は平らで、ADF は辺の上の中性子束 / 体積平均
        // の中性子束（ADR-0093 までの定義）
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(
            &solution,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        let volume: f64 = solution.volumes.iter().sum();
        for g in 0..2 {
            let average: f64 = solution
                .flux
                .iter()
                .zip(&solution.volumes)
                .map(|(f, v)| v * f[g])
                .sum::<f64>()
                / volume;
            for (side, values) in fg.adf.iter().enumerate() {
                let expected = solution.sides[side][g] / average;
                assert!(
                    (values[g] - expected).abs() < 1e-6 * expected,
                    "side {side} group {g}: {} vs {expected}",
                    values[g]
                );
            }
        }
    }

    /// 燃料と水の 2×2 の格子（右が真空）の解
    fn open_lattice() -> (Solution, Vec<Material>, GroupStructure, Boundaries) {
        let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let cell = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, cell], vec![vec![0, 1], vec![1, 0]]).unwrap();
        let boundaries = vacuum();
        let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
        let materials = vec![
            fuel.clone(),
            water.clone(),
            water.clone(),
            water.clone(),
            fuel,
            water,
        ];
        let refs: Vec<&Material> = materials.iter().collect();
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        (solution, materials, fine(), boundaries)
    }

    #[test]
    fn open_lattice_adfs_divide_by_the_one_node_flux() {
        // 真空の辺がある格子の ADF は、参照の k、ノード平均の中性子束、辺の正味の中性子流で解いた1ノードの
        // NEM の均質な面の中性子束で割る（ADR-0103）。辺は（左、右、上、下）で、-x が左、+x が右、
        // -y が下、+y が上
        use xlc_numeric::nem::{Edge, OneNode, surface_fluxes};
        let (solution, materials, fine, boundaries) = open_lattice();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(&solution, &refs, &fine, fine.bounds(), &boundaries).unwrap();
        let volume: f64 = solution.volumes.iter().sum();
        let average: Vec<f64> = (0..2)
            .map(|g| {
                solution
                    .flux
                    .iter()
                    .zip(&solution.volumes)
                    .map(|(f, v)| v * f[g])
                    .sum::<f64>()
                    / volume
            })
            .collect();
        let c = &solution.currents;
        let surface = surface_fluxes(&OneNode {
            diffusion: &fg.diffusion,
            total: &fg.total,
            scatter: &fg.scatter,
            nu_fission: &fg.nu_fission,
            chi: &fg.chi,
            k: solution.k,
            flux: &average,
            widths: solution.size,
            currents: [&c[0], &c[1], &c[3], &c[2]],
            edges: [
                Edge::Reflective,
                Edge::Vacuum,
                Edge::Reflective,
                Edge::Reflective,
            ],
        })
        .unwrap();
        for (side, face) in [(0, 0), (1, 1), (2, 3), (3, 2)] {
            for (g, homogeneous) in surface[face].iter().enumerate() {
                let expected = solution.sides[side][g] / homogeneous;
                let adf = fg.adf[side][g];
                assert!(
                    (adf - expected).abs() < 1e-12 * expected,
                    "side {side} group {g}: {adf} vs {expected}"
                );
            }
        }
        // 真空の辺では、均質な解が平らでないので、体積平均で割った値と違う
        for (g, mean) in average.iter().enumerate() {
            let flat = solution.sides[1][g] / mean;
            assert!(
                (fg.adf[1][g] - flat).abs() > 1e-2 * flat,
                "{} vs {flat}",
                fg.adf[1][g]
            );
        }
    }

    #[test]
    fn vacuum_sides_carry_the_current_ratio() {
        // 真空の辺は、MOC の辺の外向きの正味の中性子流 / 辺の上の不均質な中性子束を、少数群ごとに細かい群の
        // 和の比で持つ（ADR-0101）。鏡面反射の辺は持たない
        let (solution, materials, fine, boundaries) = open_lattice();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(&solution, &refs, &fine, fine.bounds(), &boundaries).unwrap();
        for side in [0, 2, 3] {
            assert_eq!(fg.vacuum_ratio[side], None, "side {side}");
        }
        let ratio = fg.vacuum_ratio[1].clone().unwrap();
        for (g, r) in ratio.iter().enumerate() {
            let expected = solution.currents[1][g] / solution.sides[1][g];
            assert!((r - expected).abs() < 1e-14 * expected, "{r} vs {expected}");
            // 出ていく角度中性子束だけなので P1 の 1/2 より大きく、1 より小さい
            assert!(*r > 0.5 && *r < 1.0, "{r}");
        }
        let one = homogenize(&solution, &refs, &fine, &[1.0e7, 1.0e-5], &boundaries).unwrap();
        let collapsed = one.vacuum_ratio[1].clone().unwrap();
        let expected =
            solution.currents[1].iter().sum::<f64>() / solution.sides[1].iter().sum::<f64>();
        assert!((collapsed[0] - expected).abs() < 1e-14 * expected);
        // 外周がすべて鏡面反射の格子は、どの辺も持たない
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(
            &solution,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        assert_eq!(fg.vacuum_ratio, vec![None; 4]);
    }

    #[test]
    fn adf_rejects_albedo_sides() {
        // core にアルベドの境界がないので、アルベドの辺の格子では ADF を作らない（ADR-0103）
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let albedo = Boundaries {
            top: Boundary::Albedo(0.5),
            ..Boundaries::REFLECTIVE
        };
        assert!(matches!(
            homogenize(&solution, &refs, &fine, fine.bounds(), &albedo),
            Err(HomogenizeError::Albedo)
        ));
    }

    #[test]
    fn fuel_pin_raises_the_thermal_adf() {
        // 燃料が中心にあると、外周（減速材）の熱中性子束はセルの平均より高い。4 辺は対称で等しい
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(
            &solution,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        assert!(fg.adf[0][1] > 1.0, "{:?}", fg.adf);
        for side in &fg.adf {
            for (a, b) in side.iter().zip(&fg.adf[0]) {
                assert!((a - b).abs() < 1e-6 * b, "{:?}", fg.adf);
            }
        }
    }

    #[test]
    fn one_group_collapse_preserves_k() {
        // 反応率を保存するので、1群の νΣf/Σa はセルの k∞ に一致する
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let one = homogenize(
            &solution,
            &refs,
            &fine,
            &[1.0e7, 1.0e-5],
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        assert_eq!(one.groups(), 1);
        let k = one.nu_fission[0] / one.absorption[0];
        assert!(
            (k - solution.k).abs() < 1e-6 * solution.k,
            "{k} vs {}",
            solution.k
        );
        assert_eq!(one.k_inf, solution.k);
        assert!((one.chi[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn two_group_constants_reproduce_k_inf() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(
            &solution,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        // 均質な2群の無限媒質の k = νΣfᵀ (diag Σt - Σsᵀ)⁻¹ χ（χ = (1, 0)）
        let a = [
            [fg.total[0] - fg.scatter[0][0], -fg.scatter[1][0]],
            [-fg.scatter[0][1], fg.total[1] - fg.scatter[1][1]],
        ];
        let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        let phi = [a[1][1] / det, -a[1][0] / det];
        let k = fg.nu_fission[0] * phi[0] + fg.nu_fission[1] * phi[1];
        assert!((k - solution.k).abs() < 1e-6 * solution.k);
        for g in 0..2 {
            let removal = fg.total[g] - fg.scatter[g].iter().sum::<f64>();
            assert!((fg.absorption[g] - removal).abs() < 1e-14);
        }
    }

    #[test]
    fn coarse_bounds_must_be_fine_bounds() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        assert!(
            homogenize(
                &solution,
                &refs,
                &fine,
                &[1.0e7, 1.0, 1.0e-5],
                &Boundaries::REFLECTIVE
            )
            .is_err()
        );
        assert!(
            homogenize(
                &solution,
                &refs,
                &fine,
                &[1.0e7, 0.625],
                &Boundaries::REFLECTIVE
            )
            .is_err()
        );
    }

    #[test]
    fn reflective_lattice_takes_d_from_the_b1_mode() {
        // 細かい群 = 少数群のとき、D は均質化した断面積の B1 の臨界の基本モードの J/(Bφ) そのもの
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(
            &solution,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        let p1 = |m: &Material, g: usize, to: usize| m.anisotropic[0][g][to];
        let mut medium = Medium {
            total: vec![0.0; 2],
            scatter: vec![vec![0.0; 2]; 2],
            scatter1: vec![vec![0.0; 2]; 2],
            nu_fission: vec![0.0; 2],
            chi: vec![1.0, 0.0],
        };
        for g in 0..2 {
            let mut flux = 0.0;
            for ((m, phi), v) in materials.iter().zip(&solution.flux).zip(&solution.volumes) {
                let w = v * phi[g];
                flux += w;
                medium.total[g] += w * m.total[g];
                medium.nu_fission[g] += w * m.nu_fission[g];
                for to in 0..2 {
                    medium.scatter[g][to] += w * m.scatter[g][to];
                    medium.scatter1[g][to] += w * p1(m, g, to);
                }
            }
            medium.total[g] /= flux;
            medium.nu_fission[g] /= flux;
            for to in 0..2 {
                medium.scatter[g][to] /= flux;
                medium.scatter1[g][to] /= flux;
            }
        }
        let mode = medium.critical().unwrap();
        // B² の符号は k∞ − 1 の符号
        assert_eq!(mode.buckling > 0.0, solution.k > 1.0, "{}", mode.buckling);
        for (d, e) in fg.diffusion.iter().zip(mode.diffusion()) {
            assert!((d - e).abs() < 1e-9 * e, "{d} vs {e}");
        }
        // 1群に縮約した D は B1 の中性子流と中性子束の和の比
        let one = homogenize(
            &solution,
            &refs,
            &fine,
            &[1.0e7, 1.0e-5],
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        let expected = mode.current.iter().sum::<f64>() / mode.flux.iter().sum::<f64>();
        assert!((one.diffusion[0] - expected).abs() < 1e-9 * expected);
    }

    #[test]
    fn diffusion_needs_p1_moments() {
        let (solution, mut materials, fine) = solved();
        for m in &mut materials {
            m.anisotropic.clear();
        }
        let refs: Vec<&Material> = materials.iter().collect();
        for boundaries in [
            Boundaries::REFLECTIVE,
            Boundaries {
                left: Boundary::Vacuum,
                ..Boundaries::REFLECTIVE
            },
        ] {
            assert!(matches!(
                diffusion_coefficients(&solution, &refs, &fine, fine.bounds(), &boundaries),
                Err(HomogenizeError::NoP1 { region: 0 })
            ));
        }
    }

    /// 2つの領域の、中性子束と体積がわかっている解
    fn two_regions() -> Solution {
        Solution {
            k: 1.0,
            flux: vec![vec![2.0, 1.0], vec![1.0, 3.0]],
            volumes: vec![1.0, 2.0],
            leakage: vec![0.0; 2],
            sides: vec![vec![0.0; 2]; 4],
            currents: vec![vec![0.0; 2]; 4],
            face_fluxes: vec![vec![0.0; 2]; 4],
            face_currents: vec![vec![0.0; 2]; 4],
            // [領域][x, y][群]。体積 (1, 2) を掛けて ê_G へ射影した重みは、群 0 で (0.5, 2)、群 1 で (1.0, −0.2)
            // （群 1 は領域 1 で逆向き）。どちらの群も ê = (0.6, 0.8)
            region_currents: vec![
                [vec![0.3, 0.6], vec![0.4, 0.8]],
                [vec![0.6, -0.06], vec![0.8, -0.08]],
            ],
            size: [1.0, 3.0],
            iterations: 1,
        }
    }

    fn fine() -> GroupStructure {
        GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap()
    }

    fn vacuum() -> Boundaries {
        Boundaries {
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        }
    }

    #[test]
    fn open_lattice_without_p1_scattering_uses_one_third_total() {
        // P1 が 0 なら Σtr = Σt。重みは中性子流（`two_regions` のコメント）。D = 1/(3⟨Σt⟩)
        let mut a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
        let mut b = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        for m in [&mut a, &mut b] {
            m.anisotropic = vec![vec![vec![0.0; 2]; 2]];
        }
        let d = diffusion_coefficients(
            &two_regions(),
            &[&a, &b],
            &fine(),
            fine().bounds(),
            &vacuum(),
        )
        .unwrap();
        let expected = [
            1.0 / (3.0 * (0.5 * 0.5 + 2.0 * 0.6) / 2.5),
            1.0 / (3.0 * (1.0 * 1.2 - 0.2 * 2.0) / 0.8),
        ];
        for (d, e) in d.iter().zip(expected) {
            assert!((d - e).abs() < 1e-14 * e, "{d} vs {e}");
        }
    }

    #[test]
    fn open_lattice_weights_the_outflow_transport_cross_section_by_currents() {
        // P1 = 0.3 P0: Σtr,g = Σt,g − 0.3 Σ_g' Σs0(g→g')。⟨Σtr⟩ を、ノードの正味の中性子流の向きへ射影した
        // 領域の中性子流 × 体積で重みづけし（符号つき）、D = 1/(3⟨Σtr⟩)（ADR-0104）
        let a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
        let b = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let albedo = Boundaries {
            top: Boundary::Albedo(0.5),
            ..Boundaries::REFLECTIVE
        };
        let tr = |t: f64, out: f64| t - 0.3 * out;
        let sums = [
            0.5 * tr(0.5, 0.47) + 2.0 * tr(0.6, 0.598),
            1.0 * tr(1.2, 0.9) - 0.2 * tr(2.0, 1.972),
        ];
        let expected = [1.0 / (3.0 * sums[0] / 2.5), 1.0 / (3.0 * sums[1] / 0.8)];
        let d =
            diffusion_coefficients(&two_regions(), &[&a, &b], &fine(), fine().bounds(), &albedo)
                .unwrap();
        for (d, e) in d.iter().zip(expected) {
            assert!((d - e).abs() < 1e-14 * e, "{d} vs {e}");
        }
        // 1群への縮約も同じ向き（ê = (0.6, 0.8)）への射影で重みづけする
        let one = diffusion_coefficients(
            &two_regions(),
            &[&a, &b],
            &fine(),
            &[1.0e7, 1.0e-5],
            &albedo,
        )
        .unwrap();
        let expected = 1.0 / (3.0 * (sums[0] + sums[1]) / 3.3);
        assert!(
            (one[0] - expected).abs() < 1e-14 * expected,
            "{one:?} vs {expected}"
        );
    }

    #[test]
    fn current_weighting_counts_a_void_gap_with_its_current() {
        // 真空の隙間（ピンのギャップ）は Σtr = 0 で、中性子流は重みに入る
        let a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
        let void = Material::void(2);
        let d = diffusion_coefficients(
            &two_regions(),
            &[&a, &void],
            &fine(),
            fine().bounds(),
            &vacuum(),
        )
        .unwrap();
        let expected = [
            1.0 / (3.0 * 0.5 * tr_of(&a, 0) / 2.5),
            1.0 / (3.0 * 1.0 * tr_of(&a, 1) / 0.8),
        ];
        for (d, e) in d.iter().zip(expected) {
            assert!((d - e).abs() < 1e-14 * e, "{d} vs {e}");
        }
    }

    #[test]
    fn current_weighting_needs_a_net_current_and_a_positive_mean() {
        let a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
        let b = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        // ノードの正味の中性子流が 0（対称な位置）なら、向きが定まらない
        let mut zero = two_regions();
        zero.region_currents = vec![[vec![0.0; 2], vec![0.0; 2]]; 2];
        assert!(matches!(
            diffusion_coefficients(&zero, &[&a, &b], &fine(), fine().bounds(), &vacuum()),
            Err(HomogenizeError::NoNetCurrent { group: 0 })
        ));
        // 逆向きの重みが大きく、⟨Σtr⟩ が正でない
        let mut opposite = two_regions();
        // 群 0 の重みは (1.0, −0.9)（和 0.1）で、⟨Σtr⟩ = (0.359 − 0.9 × 0.4206)/0.1 < 0
        opposite.region_currents = vec![
            [vec![1.0, 0.6], vec![0.0, 0.8]],
            [vec![-0.45, -0.06], vec![0.0, -0.08]],
        ];
        assert!(matches!(
            diffusion_coefficients(&opposite, &[&a, &b], &fine(), fine().bounds(), &vacuum()),
            Err(HomogenizeError::Transport { group: 0, .. })
        ));
    }

    #[test]
    fn open_lattice_d_uses_the_moc_region_currents() {
        // 格子の解の領域の中性子流で、手で重みづけした値と一致する
        let (solution, materials, fine, boundaries) = open_lattice();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(&solution, &refs, &fine, fine.bounds(), &boundaries).unwrap();
        for g in 0..2 {
            let net: [f64; 2] = std::array::from_fn(|c| {
                solution
                    .region_currents
                    .iter()
                    .zip(&solution.volumes)
                    .map(|(j, v)| v * j[c][g])
                    .sum()
            });
            let length = net[0].hypot(net[1]);
            let (mut weight, mut sum) = (0.0, 0.0);
            for ((j, v), m) in solution
                .region_currents
                .iter()
                .zip(&solution.volumes)
                .zip(&materials)
            {
                let w = v * (j[0][g] * net[0] + j[1][g] * net[1]) / length;
                weight += w;
                sum += w * tr_of(m, g);
            }
            let expected = 1.0 / (3.0 * sum / weight);
            assert!(
                (fg.diffusion[g] - expected).abs() < 1e-12 * expected,
                "group {g}"
            );
        }
    }

    fn tr_of(m: &Material, g: usize) -> f64 {
        m.total[g] - m.anisotropic[0][g].iter().sum::<f64>()
    }

    #[test]
    fn reflective_lattice_without_fission_uses_the_p1_current_of_the_moc_flux() {
        // 核分裂がないと臨界のバックリングは定まらない。B² = 0 の中性子流の式
        // 3 ⟨Σt⟩_g d_g − 3 Σ_g' ⟨Σs1⟩(g'→g) d_g' = Φ_g を MOC の中性子束 Φ で解き、D_G = Σ d / Σ Φ
        let a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
        let b = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let flux = [4.0, 7.0];
        let hom = |g: usize, f: &dyn Fn(&Material) -> f64| {
            let w = [[2.0, 1.0], [2.0, 6.0]];
            (w[0][g] * f(&a) + w[1][g] * f(&b)) / flux[g]
        };
        let t0 = hom(0, &|m| m.total[0]);
        let t1 = hom(1, &|m| m.total[1]);
        let s = |from: usize, to: usize| hom(from, &|m| m.anisotropic[0][from][to]);
        // [3(t0 − s00), −3 s10; −3 s01, 3(t1 − s11)] d = Φ
        let m = [
            [3.0 * (t0 - s(0, 0)), -3.0 * s(1, 0)],
            [-3.0 * s(0, 1), 3.0 * (t1 - s(1, 1))],
        ];
        let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        let d = [
            (m[1][1] * flux[0] - m[0][1] * flux[1]) / det,
            (m[0][0] * flux[1] - m[1][0] * flux[0]) / det,
        ];
        let got = diffusion_coefficients(
            &two_regions(),
            &[&a, &b],
            &fine(),
            fine().bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        for g in 0..2 {
            let e = d[g] / flux[g];
            assert!((got[g] - e).abs() < 1e-13 * e, "{} vs {e}", got[g]);
        }
        let one = diffusion_coefficients(
            &two_regions(),
            &[&a, &b],
            &fine(),
            &[1.0e7, 1.0e-5],
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        let e = (d[0] + d[1]) / (flux[0] + flux[1]);
        assert!((one[0] - e).abs() < 1e-13 * e);
    }

    /// 2つの核分裂性の物質（χ が違い、熱群から高速群への上方散乱がある）
    fn fissile(
        total: [f64; 2],
        scatter: [[f64; 2]; 2],
        p1: [[f64; 2]; 2],
        nu_fission: [f64; 2],
        chi: [f64; 2],
    ) -> Material {
        Material {
            total: total.to_vec(),
            scatter: scatter.iter().map(|r| r.to_vec()).collect(),
            anisotropic: vec![p1.iter().map(|r| r.to_vec()).collect()],
            fission: nu_fission.iter().map(|v| v / 2.5).collect(),
            nu_fission: nu_fission.to_vec(),
            fission_terms: vec![FissionTerm {
                chi: chi.to_vec(),
                nu_fission: nu_fission.to_vec(),
            }],
        }
    }

    /// 2群の B1 の基本モード（手で書いた 2 × 2 の式）。(k, φ, d)
    fn oracle_mode(
        total: [f64; 2],
        s0: [[f64; 2]; 2],
        s1: [[f64; 2]; 2],
        nu_fission: [f64; 2],
        chi: [f64; 2],
        b2: f64,
    ) -> (f64, [f64; 2], [f64; 2]) {
        let gamma = |t: f64| {
            let x = b2.abs().sqrt() / t;
            if b2 == 0.0 {
                1.0
            } else if b2 > 0.0 {
                x * x * x.atan() / (3.0 * (x - x.atan()))
            } else {
                x * x * x.atanh() / (3.0 * (x.atanh() - x))
            }
        };
        // d = M⁻¹ φ、M = [3γ0 t0 − 3 s1(0→0), −3 s1(1→0); −3 s1(0→1), 3γ1 t1 − 3 s1(1→1)]
        let m = [
            [
                3.0 * gamma(total[0]) * total[0] - 3.0 * s1[0][0],
                -3.0 * s1[1][0],
            ],
            [
                -3.0 * s1[0][1],
                3.0 * gamma(total[1]) * total[1] - 3.0 * s1[1][1],
            ],
        ];
        let det_m = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        let inv = [
            [m[1][1] / det_m, -m[0][1] / det_m],
            [-m[1][0] / det_m, m[0][0] / det_m],
        ];
        // A = diag(t) − S0ᵀ + B² M⁻¹、A φ = χ
        let a = [
            [
                total[0] - s0[0][0] + b2 * inv[0][0],
                -s0[1][0] + b2 * inv[0][1],
            ],
            [
                -s0[0][1] + b2 * inv[1][0],
                total[1] - s0[1][1] + b2 * inv[1][1],
            ],
        ];
        let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        let phi = [
            (a[1][1] * chi[0] - a[0][1] * chi[1]) / det,
            (a[0][0] * chi[1] - a[1][0] * chi[0]) / det,
        ];
        let d = [
            inv[0][0] * phi[0] + inv[0][1] * phi[1],
            inv[1][0] * phi[0] + inv[1][1] * phi[1],
        ];
        (nu_fission[0] * phi[0] + nu_fission[1] * phi[1], phi, d)
    }

    #[test]
    fn b1_diffusion_matches_a_hand_computed_two_region_lattice() {
        let a = fissile(
            [0.5, 1.2],
            [[0.42, 0.05], [0.01, 0.95]],
            [[0.15, 0.01], [0.002, 0.3]],
            [0.016, 0.45],
            [1.0, 0.0],
        );
        let b = fissile(
            [0.6, 1.6],
            [[0.5, 0.08], [0.02, 1.4]],
            [[0.2, 0.02], [0.004, 0.5]],
            [0.008, 0.19],
            [0.6, 0.4],
        );
        let solution = Solution {
            k: 1.0,
            flux: vec![vec![2.0, 1.0], vec![1.5, 2.5]],
            volumes: vec![1.0, 0.5],
            leakage: vec![0.0; 2],
            sides: vec![vec![1.0; 2]; 4],
            currents: vec![vec![0.0; 2]; 4],
            face_fluxes: vec![vec![1.0; 2]; 4],
            face_currents: vec![vec![0.0; 2]; 4],
            region_currents: vec![[vec![0.0; 2], vec![0.0; 2]]; 2],
            size: [1.0, 0.5],
            iterations: 1,
        };
        // 重み V φ: A は (2, 1)、B は (0.75, 1.25)。中性子束 (2.75, 2.25)
        let w = [[2.0, 1.0], [0.75, 1.25]];
        let flux = [2.75, 2.25];
        let hom = |g: usize, xa: f64, xb: f64| (w[0][g] * xa + w[1][g] * xb) / flux[g];
        let total = [hom(0, 0.5, 0.6), hom(1, 1.2, 1.6)];
        let s0 = [
            [hom(0, 0.42, 0.5), hom(0, 0.05, 0.08)],
            [hom(1, 0.01, 0.02), hom(1, 0.95, 1.4)],
        ];
        let s1 = [
            [hom(0, 0.15, 0.2), hom(0, 0.01, 0.02)],
            [hom(1, 0.002, 0.004), hom(1, 0.3, 0.5)],
        ];
        let nu_fission = [hom(0, 0.016, 0.008), hom(1, 0.45, 0.19)];
        let rate_a = 2.0 * 0.016 + 1.0 * 0.45;
        let rate_b = 0.75 * 0.008 + 1.25 * 0.19;
        let chi = [
            (rate_a + 0.6 * rate_b) / (rate_a + rate_b),
            0.4 * rate_b / (rate_a + rate_b),
        ];
        let k = |b2: f64| oracle_mode(total, s0, s1, nu_fission, chi, b2);
        // 二分法（k は B² について単調に減る）
        let k_inf = k(0.0).0;
        assert!(k_inf > 1.0, "{k_inf}");
        let (mut low, mut high) = (0.0, 1.0);
        assert!(k(high).0 < 1.0);
        for _ in 0..200 {
            let middle = 0.5 * (low + high);
            if k(middle).0 > 1.0 {
                low = middle;
            } else {
                high = middle;
            }
        }
        let (_, phi, d) = k(0.5 * (low + high));
        let fg = homogenize(
            &solution,
            &[&a, &b],
            &fine(),
            fine().bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        for g in 0..2 {
            let e = d[g] / phi[g];
            assert!(
                (fg.diffusion[g] - e).abs() < 1e-8 * e,
                "{} vs {e}",
                fg.diffusion[g]
            );
        }
        // 1群: B1 の中性子束で縮約する。MOC の中性子束で重みづけした値とは違う
        let one = homogenize(
            &solution,
            &[&a, &b],
            &fine(),
            &[1.0e7, 1.0e-5],
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        let e = (d[0] + d[1]) / (phi[0] + phi[1]);
        assert!(
            (one.diffusion[0] - e).abs() < 1e-8 * e,
            "{} vs {e}",
            one.diffusion[0]
        );
        let moc = (d[0] / phi[0] * flux[0] + d[1] / phi[1] * flux[1]) / (flux[0] + flux[1]);
        assert!((moc - e).abs() > 1e-3 * e, "{moc} vs {e}");
    }

    /// 燃料のピンと水のセルを横に並べた帯（1行2列、右が真空）と、その解
    fn strip() -> (Lattice, Solution, Vec<Material>, Boundaries) {
        let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let cell = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::rectilinear(
            vec![1.26; 2],
            vec![1.26],
            vec![Cell::Pin(pin), Cell::Pin(cell)],
            vec![vec![0, 1]],
        )
        .unwrap();
        let boundaries = vacuum();
        let tracking = Tracking::with_boundaries(&lattice, 32, 0.02, boundaries).unwrap();
        let materials = vec![fuel, water.clone(), water];
        let refs: Vec<&Material> = materials.iter().collect();
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        (lattice, solution, materials, boundaries)
    }

    fn node(columns: std::ops::Range<usize>, rows: std::ops::Range<usize>) -> Node {
        Node {
            columns,
            rows,
            reflector: false,
        }
    }

    #[test]
    fn nodes_of_an_infinite_medium_are_the_medium() {
        // 1つの物質の格子（外周がすべて鏡面反射）を2つのノードに分けると、どちらも物質の断面積そのもので、
        // 面の正味の中性子流は 0、ADF は 1。複数のノードの格子の拡散係数は outflow（ADR-0098）で、
        // 1/(3(Σt − Σs1))
        let m = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let lattice = Lattice::new(vec![cell], vec![vec![0, 0], vec![0, 0]]).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.05).unwrap();
        let refs = vec![&m; 8];
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        let fine = fine();
        let nodes = [node(0..1, 0..2), node(1..2, 0..2)];
        let fgs = homogenize_nodes(
            &solution,
            &lattice,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
            &nodes,
        )
        .unwrap();
        assert_eq!(fgs.len(), 2);
        for fg in &fgs {
            for g in 0..2 {
                assert!((fg.total[g] - m.total[g]).abs() < 1e-12, "{:?}", fg.total);
                assert!((fg.nu_fission[g] - m.nu_fission[g]).abs() < 1e-12);
                let transport = m.total[g] - m.anisotropic[0][g].iter().sum::<f64>();
                let d = 1.0 / (3.0 * transport);
                assert!(
                    (fg.diffusion[g] - d).abs() < 1e-12 * d,
                    "{:?}",
                    fg.diffusion
                );
            }
            for value in fg.adf.iter().flatten() {
                assert!((value - 1.0).abs() < 1e-6, "{:?}", fg.adf);
            }
            assert_eq!(fg.vacuum_ratio, vec![None; 4]);
            assert_eq!(fg.normalized_adf, vec![false; 4]);
        }
    }

    #[test]
    fn strip_nodes_make_get_constants_and_normalized_adfs() {
        // 燃料と水の帯を、燃料のノードと反射体のノードに分ける（ADR-0095）
        use xlc_numeric::nem::{Edge, OneNode, surface_fluxes};
        let (lattice, solution, materials, boundaries) = strip();
        let refs: Vec<&Material> = materials.iter().collect();
        let fine = fine();
        let nodes = [
            node(0..1, 0..1),
            Node {
                reflector: true,
                ..node(1..2, 0..1)
            },
        ];
        let fgs = homogenize_nodes(
            &solution,
            &lattice,
            &refs,
            &fine,
            fine.bounds(),
            &boundaries,
            &nodes,
        )
        .unwrap();
        let (fuel, reflector) = (&fgs[0], &fgs[1]);
        // 定数は、ノードの領域だけの反応率の保存（燃料は領域 0, 1、反射体は領域 2）
        let average = |regions: std::ops::Range<usize>, g: usize, x: &dyn Fn(&Material) -> f64| {
            let (mut rate, mut flux) = (0.0, 0.0);
            for r in regions {
                let w = solution.volumes[r] * solution.flux[r][g];
                rate += w * x(&materials[r]);
                flux += w;
            }
            (rate / flux, flux / 1.26_f64.powi(2))
        };
        let mut fluxes = [[0.0; 2]; 2];
        for (node_fluxes, (fg, regions)) in fluxes.iter_mut().zip([(fuel, 0..2), (reflector, 2..3)])
        {
            for (g, node_flux) in node_fluxes.iter_mut().enumerate() {
                let (total, flux) = average(regions.clone(), g, &|m| m.total[g]);
                assert!((fg.total[g] - total).abs() < 1e-12 * total);
                // 燃料のノードの D は、領域の中性子流をノードの正味の中性子流の向きへ射影した重みの outflow の
                // Σtr から（ADR-0104）。反射体のノードの D は等価な D（下で x の面の不連続因子がそろうことを見る）
                if regions.start == 0 {
                    let net: [f64; 2] = std::array::from_fn(|c| {
                        regions
                            .clone()
                            .map(|r| solution.volumes[r] * solution.region_currents[r][c][g])
                            .sum()
                    });
                    let length = net[0].hypot(net[1]);
                    let (mut weight, mut sum) = (0.0, 0.0);
                    for r in regions.clone() {
                        let j = &solution.region_currents[r];
                        let w =
                            solution.volumes[r] * (j[0][g] * net[0] + j[1][g] * net[1]) / length;
                        let m = &materials[r];
                        weight += w;
                        sum += w * (m.total[g] - m.anisotropic[0][g].iter().sum::<f64>());
                    }
                    let d = 1.0 / (3.0 * sum / weight);
                    assert!(
                        (fg.diffusion[g] - d).abs() < 1e-12 * d,
                        "{} vs {d}",
                        fg.diffusion[g]
                    );
                }
                *node_flux = flux;
            }
            assert_eq!(fg.k_inf, solution.k);
        }
        // 面の値: 縦の面は 0（左）、1（燃料と反射体の境）、2（右）、横の面は 3, 4（上）と 5, 6（下）
        let het = |face: usize| solution.face_fluxes[face].clone();
        let out = |face: usize, sign: f64| -> Vec<f64> {
            solution.face_currents[face]
                .iter()
                .map(|j| sign * j)
                .collect()
        };
        let one_node =
            |fg: &FewGroup, flux: &[f64], faces: [usize; 4], signs: [f64; 4], edges: [Edge; 4]| {
                let c: Vec<Vec<f64>> = (0..4).map(|i| out(faces[i], signs[i])).collect();
                surface_fluxes(&OneNode {
                    diffusion: &fg.diffusion,
                    total: &fg.total,
                    scatter: &fg.scatter,
                    nu_fission: &fg.nu_fission,
                    chi: &fg.chi,
                    k: solution.k,
                    flux,
                    widths: [1.26, 1.26],
                    currents: [&c[0], &c[1], &c[2], &c[3]],
                    edges,
                })
                .unwrap()
            };
        // 隣の横方向（y）の漏れは、鏡面反射の上下の面の正味の中性子流（収束の許容差の程度で 0）
        let y_leakage = |top: usize, bottom: usize| -> Vec<f64> {
            (0..2)
                .map(|g| {
                    (-solution.face_currents[top][g] + solution.face_currents[bottom][g]) / 1.26
                })
                .collect()
        };
        let (fuel_y, reflector_y) = (y_leakage(3, 5), y_leakage(4, 6));
        let fuel_surface = one_node(
            fuel,
            &fluxes[0],
            [0, 1, 5, 3],
            [-1.0, 1.0, 1.0, -1.0],
            [
                Edge::Reflective,
                Edge::Neighbor {
                    width: 1.26,
                    leakage: &reflector_y,
                },
                Edge::Reflective,
                Edge::Reflective,
            ],
        );
        let reflector_surface = one_node(
            reflector,
            &fluxes[1],
            [1, 2, 6, 4],
            [-1.0, 1.0, 1.0, -1.0],
            [
                Edge::Neighbor {
                    width: 1.26,
                    leakage: &fuel_y,
                },
                Edge::Vacuum,
                Edge::Reflective,
                Edge::Reflective,
            ],
        );
        for g in 0..2 {
            // 燃料の右の辺の ADF は、境の面の不均質な中性子束 / 均質な面の中性子束
            let expected = het(1)[g] / fuel_surface[1][g];
            assert!((fuel.adf[1][g] - expected).abs() < 1e-10 * expected);
            // 反射体の左の辺は f* = f_反射体 / f_燃料
            let f_reflector = het(1)[g] / reflector_surface[0][g];
            let expected = f_reflector / fuel.adf[1][g];
            assert!(
                (reflector.adf[0][g] - expected).abs() < 1e-10 * expected,
                "{} vs {expected}",
                reflector.adf[0][g]
            );
            // 真空の辺の ADF と J/φ
            let expected = het(2)[g] / reflector_surface[1][g];
            assert!((reflector.adf[1][g] - expected).abs() < 1e-10 * expected);
            // 反射体のノードの x の2つの面の不連続因子は等しい（Koebke の等価な D。ADR-0104）
            let left = reflector.adf[0][g] * fuel.adf[1][g];
            assert!(
                (left - reflector.adf[1][g]).abs() < 1e-9 * left,
                "group {g}: {left} vs {}",
                reflector.adf[1][g]
            );
            let ratio = reflector.vacuum_ratio[1].as_ref().unwrap()[g];
            let expected = solution.face_currents[2][g] / het(2)[g];
            assert!((ratio - expected).abs() < 1e-12 * expected);
        }
        assert_eq!(fuel.normalized_adf, vec![false; 4]);
        assert_eq!(reflector.normalized_adf, vec![true, false, false, false]);
        assert_eq!(fuel.vacuum_ratio, vec![None; 4]);
        assert_eq!(
            reflector
                .vacuum_ratio
                .iter()
                .map(Option::is_some)
                .collect::<Vec<_>>(),
            vec![false, true, false, false]
        );
    }

    #[test]
    fn invalid_nodes_are_errors() {
        let (lattice, solution, materials, boundaries) = strip();
        let refs: Vec<&Material> = materials.iter().collect();
        let fine = fine();
        let run = |nodes: &[Node], boundaries: &Boundaries| {
            homogenize_nodes(
                &solution,
                &lattice,
                &refs,
                &fine,
                fine.bounds(),
                boundaries,
                nodes,
            )
        };
        let reflector = |n: Node| Node {
            reflector: true,
            ..n
        };
        let cases: Vec<(Vec<Node>, &str)> = vec![
            (vec![], "no nodes"),
            (vec![node(0..1, 0..1)], "not covered"),
            (vec![node(0..2, 0..1), node(1..2, 0..1)], "overlap"),
            (vec![node(0..1, 0..1), node(1..3, 0..1)], "outside"),
            (vec![node(0..1, 0..1), node(1..1, 0..1)], "empty"),
            (
                vec![reflector(node(0..1, 0..1)), reflector(node(1..2, 0..1))],
                "no fuel",
            ),
        ];
        for (nodes, what) in &cases {
            assert!(
                matches!(run(nodes, &boundaries), Err(HomogenizeError::Nodes(_))),
                "{what}"
            );
        }
        let albedo = Boundaries {
            top: Boundary::Albedo(0.5),
            ..boundaries
        };
        assert!(matches!(
            run(&[node(0..1, 0..1), node(1..2, 0..1)], &albedo),
            Err(HomogenizeError::Albedo)
        ));
        // 隣が1つのノードで面がそろわない（左のノードの右の辺に2つのノードが接する）
        let m = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let square = Lattice::new(vec![cell], vec![vec![0, 0], vec![0, 0]]).unwrap();
        let tracking = Tracking::new(&square, 16, 0.05).unwrap();
        let refs = vec![&m; 8];
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        let misaligned = homogenize_nodes(
            &solution,
            &square,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
            &[node(0..1, 0..2), node(1..2, 0..1), node(1..2, 1..2)],
        );
        assert!(
            matches!(misaligned, Err(HomogenizeError::Nodes(_))),
            "{misaligned:?}"
        );
    }
}
