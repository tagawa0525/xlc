//! MOC の解から、少数群定数（`xlc_format::fewgroup::FewGroup`）を作る。

use std::ops::Range;

use xlc_format::fewgroup::{FewGroup, Pins};
use xlc_format::group::GroupStructure;
use xlc_numeric::equivalence::equivalent_diffusion;
use xlc_numeric::nem::{Edge, OneNode, OneNodeError, is_fissile, surface_fluxes};
use xlc_numeric::pins::corners::{CornerNode, corner_balance, fuel_harmonic_mean, nodes_at_point};
use xlc_numeric::pins::node::{Conditions, Constants, NodeModel};
use xlc_numeric::pins::{CELL_POINTS, FACE_POINTS};
use xlc_numeric::quadrature::gauss_legendre;

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
    #[error("pins: {0}")]
    Pins(String),
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
///   （格子全体が1つのノード）なら B1、ほかは、核分裂のあるノードではノードの領域の中性子流で重みづけした outflow、
///   核分裂のないノードでは Koebke の等価理論の D（正味の中性子流の向きの2つの面の不連続因子をそろえる。
///   ADR-0104）。`k_inf` は格子の k
/// - 辺の不均質な中性子束と正味の中性子流は、辺の上の粗メッシュの面の値（MOC の `face_fluxes`、
///   `face_currents`）を、面の長さで重みづけして平均する。内側の辺も同じ
/// - ADF は1ノードの問題（ADR-0094）で作る。外周の辺は格子の境界（鏡面反射か真空）、内側の辺は隣の
///   ノード（その幅と、参照の横方向の漏れ。core と同じ2次式）とする。内側の辺は、向かいにちょうど1つの
///   ノードがあり、辺の長さがそろうこと
/// - 外周の真空の辺には J/φ_het を持たせる（ADR-0101）
/// - 反射体のノード（`Node::reflector`）の、反射体でないノードと接する辺の ADF は、その隣のノードの同じ面の
///   ADF で割った f* とし、`normalized_adf` に印を付ける。反射体のノードは、反射体でないノードと少なくとも
///   1辺で接すること
/// - `pins` なら、形状関数と角の不連続因子（`node_pins`。ADR-0125）も付ける
#[allow(clippy::too_many_arguments)]
pub fn homogenize_nodes(
    solution: &Solution,
    lattice: &Lattice,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
    boundaries: &Boundaries,
    nodes: &[Node],
    pins: bool,
) -> Result<Vec<FewGroup>, HomogenizeError> {
    let NodeReferences {
        references: homogenized,
        group_of,
        neighbors,
    } = node_references(
        solution,
        lattice,
        materials,
        fine,
        coarse_bounds,
        boundaries,
        nodes,
    )?;
    let sides = sides_of(boundaries);
    let coarse = coarse_bounds.len() - 1;
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
        // 核分裂のないノードの拡散係数は、Koebke の等価理論の D（ADR-0104）。中性子流の重みの D から探す
        let mut few = h.few.clone();
        let fissile = is_fissile(&few.nu_fission)
            .map_err(|e| HomogenizeError::OneNode(OneNodeError::Invalid(e.to_string())))?;
        if !fissile {
            few.diffusion = equivalent_diffusion(
                &OneNode {
                    diffusion: &h.few.diffusion,
                    total: &few.total,
                    scatter: &few.scatter,
                    nu_fission: &few.nu_fission,
                    chi: &few.chi,
                    k: solution.k,
                    flux: &h.average,
                    widths: h.size,
                    currents: FACES.map(|side| h.currents[side].as_slice()),
                    edges,
                },
                FACES.map(|side| h.heterogeneous[side].as_slice()),
            )?;
            tracing::info!(
                node = n,
                current_weighted = ?h.few.diffusion,
                equivalent = ?few.diffusion,
                "equivalent diffusion coefficients"
            );
        }
        let adf = discontinuity_factors(
            &few,
            solution.k,
            &h.average,
            h.size,
            &h.heterogeneous,
            &h.currents,
            edges,
        )?;
        few.adf = adf;
        adfs.push(few);
    }
    let shapes = if pins {
        Some(node_pins(
            solution,
            lattice,
            materials,
            &group_of,
            nodes,
            &neighbors,
            &homogenized,
            &adfs,
        )?)
    } else {
        None
    };
    let mut result = Vec::with_capacity(nodes.len());
    for (n, node) in nodes.iter().enumerate() {
        let mut few = adfs[n].clone();
        few.pins = shapes.as_ref().map(|p| p[n].clone());
        if node.reflector {
            for side in 0..4 {
                if let Some(m) = neighbors[n][side]
                    && !nodes[m].reflector
                {
                    for (f, fuel) in few.adf[side].iter_mut().zip(&adfs[m].adf[OPPOSITE[side]]) {
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

/// ノードの均質化の参照の量: ノードの少数群定数（ADF は空）、ノード平均の中性子束、ノードの幅と高さ、辺（左、右、上、
/// 下）の不均質な中性子束と外向きの正味の中性子流。
#[derive(Debug, Clone)]
pub(crate) struct Reference {
    pub(crate) few: FewGroup,
    pub(crate) average: Vec<f64>,
    pub(crate) size: [f64; 2],
    pub(crate) heterogeneous: Vec<Vec<f64>>,
    pub(crate) currents: Vec<Vec<f64>>,
}

/// ノードの均質化の参照の量と、確かめたノードの並び。
#[derive(Debug, Clone)]
pub(crate) struct NodeReferences {
    /// [ノード]
    pub(crate) references: Vec<Reference>,
    /// [ライブラリの群] の少数群
    pub(crate) group_of: Vec<usize>,
    /// [ノード][辺（左、右、上、下）] の、辺の向かいのノード（`node_neighbors`）
    pub(crate) neighbors: Vec<[Option<usize>; 4]>,
}

/// ノードごとの `Reference`（`homogenize_nodes` の前半）。ノードの並びと、アルベドの辺がないことを確かめる
pub(crate) fn node_references(
    solution: &Solution,
    lattice: &Lattice,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
    boundaries: &Boundaries,
    nodes: &[Node],
) -> Result<NodeReferences, HomogenizeError> {
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
    let mut references = Vec::with_capacity(nodes.len());
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
        references.push(Reference {
            few,
            average,
            size: [widths[c0..c1].iter().sum(), heights[r0..r1].iter().sum()],
            heterogeneous,
            currents,
        });
    }
    Ok(NodeReferences {
        references,
        group_of,
        neighbors,
    })
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

/// ピン出力の再構成の形状関数と角の不連続因子（ADR-0096）。格子全体を1つのノードとし、外周がすべて鏡面反射
/// の格子でだけ作る（均質な解が平らなので、角の均質な中性子束はノード平均。ADR-0094。ほかの格子は #0020）。
///
/// - セルの大きさは、格子の列の幅と行の高さ（行は下から）
/// - 形状関数 f_g（セル）= セルの少数群 g の核分裂率の平均 / ノードの少数群 g の核分裂率の平均。核分裂のないセルは 0。
///   core はセルの幾何の面積で再構成するので、面積で重みを付けた平均がちょうど 1 になるよう、トラックから推定した
///   体積でなく解析的な体積で平均する
/// - 角の不連続因子 = 余白の角のセルの中性子束の平均 / ノード平均の中性子束。面の ADF（`homogenize`）と同じ均質な
///   中性子束を基準にするよう、トラックから推定した体積で平均する。角のセルがピンのセル（余白がない）ならエラー
pub fn pin_shapes(
    solution: &Solution,
    lattice: &Lattice,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
    boundaries: &Boundaries,
) -> Result<Pins, HomogenizeError> {
    if !reflective(boundaries) {
        return Err(HomogenizeError::Pins(
            "form functions need a lattice with every side reflective (other lattices: #0020)"
                .into(),
        ));
    }
    check_regions(solution, materials)?;
    let volumes = lattice.volumes();
    if volumes.len() != solution.flux.len() {
        return Err(HomogenizeError::RegionCount {
            materials: volumes.len(),
            regions: solution.flux.len(),
        });
    }
    let group_of = coarse_groups(fine, coarse_bounds)?;
    let coarse = coarse_bounds.len() - 1;
    let (columns, rows) = (lattice.columns(), lattice.rows());
    // [セル][少数群] の、解析的な体積での核分裂率の積分と、トラックの体積での中性子束の積分。セルの面積も両方
    let mut fission = vec![vec![0.0; coarse]; lattice.cells()];
    let mut flux = vec![vec![0.0; coarse]; lattice.cells()];
    let (mut area, mut track_area) = (vec![0.0; lattice.cells()], vec![0.0; lattice.cells()]);
    for cell in 0..lattice.cells() {
        for region in lattice.cell_regions(cell) {
            let (m, phi) = (materials[region], &solution.flux[region]);
            let (volume, track) = (volumes[region], solution.volumes[region]);
            area[cell] += volume;
            track_area[cell] += track;
            for (g, &big) in group_of.iter().enumerate() {
                fission[cell][big] += volume * phi[g] * m.fission[g];
                flux[cell][big] += track * phi[g];
            }
        }
    }
    let node = |values: &[Vec<f64>], areas: &[f64], g: usize| -> f64 {
        values.iter().map(|v| v[g]).sum::<f64>() / areas.iter().sum::<f64>()
    };
    let mut form = Vec::with_capacity(coarse);
    for g in 0..coarse {
        let mean = node(&fission, &area, g);
        if mean.is_nan() || mean <= 0.0 {
            return Err(HomogenizeError::Pins(format!(
                "few group {g} has no fission; form functions need fission"
            )));
        }
        // 行は下から（格子の行は上から）
        form.push(
            (0..rows)
                .rev()
                .map(|r| {
                    (0..columns)
                        .map(|c| fission[r * columns + c][g] / area[r * columns + c] / mean)
                        .collect()
                })
                .collect(),
        );
    }
    // 角（-x -y、+x -y、+x +y、-x +y）のセル。格子の行は上から
    let corner_cells = [(rows - 1) * columns, rows * columns - 1, columns - 1, 0];
    let mut corners = Vec::with_capacity(4);
    for cell in corner_cells {
        if lattice.cell_type(cell).pin().is_some() {
            return Err(HomogenizeError::Pins(
                "corner discontinuity factors need a margin (gap) at the lattice corners".into(),
            ));
        }
        corners.push(
            (0..coarse)
                .map(|g| flux[cell][g] / track_area[cell] / node(&flux, &track_area, g))
                .collect(),
        );
    }
    Ok(Pins {
        dx: lattice.widths().to_vec(),
        dy: lattice.heights().iter().rev().copied().collect(),
        form,
        corners,
        normalized_corners: vec![false; 4],
        halves: vec![false; 4],
    })
}

/// 点 (x, y)（格子の -x、-y の辺から）の [少数群] の不均質な中性子束: 点を頂点とする 4 つの象限のうち格子の中にある
/// ものの領域の中性子束の平均（ADR-0108）
fn point_flux(
    solution: &Solution,
    lattice: &Lattice,
    group_of: &[usize],
    coarse: usize,
    x: f64,
    y: f64,
) -> Vec<f64> {
    let (width, height) = (
        lattice.widths().iter().sum::<f64>(),
        lattice.heights().iter().sum::<f64>(),
    );
    let e = 1e-9 * width.max(height);
    let mut sum = vec![0.0; coarse];
    let mut count = 0.0;
    for (sx, sy) in [(1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0)] {
        let (px, py) = (x + sx * e, y + sy * e);
        if px <= 0.0 || px >= width || py <= 0.0 || py >= height {
            continue;
        }
        let region = lattice.region_at(px - 0.5 * width, py - 0.5 * height);
        for (g, &big) in group_of.iter().enumerate() {
            sum[big] += solution.flux[region][g];
        }
        count += 1.0;
    }
    sum.iter().map(|s| s / count).collect()
}

/// 格子の上の nx × ny の長方形の並びにしたノード。
#[derive(Debug, Clone, PartialEq)]
pub struct NodeGrid {
    /// 列の切れ目（格子の列は左から）
    pub column_cuts: Vec<usize>,
    /// 行の切れ目（格子の行は上から）
    pub row_cuts: Vec<usize>,
    pub nx: usize,
    pub ny: usize,
    /// [i + nx j]（j は -y から）のノードの番号
    pub grid: Vec<usize>,
}

/// ノードの列と行の範囲の切れ目で、ノードが nx × ny の長方形の並びをちょうど 1 回ずつ埋めることを確かめる
/// （形状関数と角の不連続因子の角の点のバランス。ADR-0125）
pub fn node_grid(nodes: &[Node]) -> Result<NodeGrid, String> {
    let cuts = |ranges: &mut dyn Iterator<Item = &Range<usize>>| -> Vec<usize> {
        let mut c: Vec<usize> = ranges.flat_map(|r| [r.start, r.end]).collect();
        c.sort_unstable();
        c.dedup();
        c
    };
    let column_cuts = cuts(&mut nodes.iter().map(|n| &n.columns));
    let row_cuts = cuts(&mut nodes.iter().map(|n| &n.rows));
    if column_cuts.len() < 2 || row_cuts.len() < 2 {
        return Err("the nodes do not make a grid".into());
    }
    let (nx, ny) = (column_cuts.len() - 1, row_cuts.len() - 1);
    let mut grid = vec![None; nx * ny];
    for (n, node) in nodes.iter().enumerate() {
        let i = column_cuts.iter().position(|&c| c == node.columns.start);
        let r = row_cuts.iter().position(|&c| c == node.rows.start);
        let slot = match (i, r) {
            (Some(i), Some(r))
                if column_cuts.get(i + 1) == Some(&node.columns.end)
                    && row_cuts.get(r + 1) == Some(&node.rows.end) =>
            {
                &mut grid[i + nx * (ny - 1 - r)]
            }
            _ => {
                return Err(format!(
                    "node {n} is not one rectangle of the grid that the nodes make"
                ));
            }
        };
        if slot.replace(n).is_some() {
            return Err(format!("node {n} overlaps another node of the grid"));
        }
    }
    let grid = grid
        .into_iter()
        .collect::<Option<Vec<usize>>>()
        .ok_or_else(|| "the nodes do not fill the grid".to_string())?;
    Ok(NodeGrid {
        column_cuts,
        row_cuts,
        nx,
        ny,
        grid,
    })
}

/// 格子のノードごとの形状関数と角の不連続因子（ADR-0125）。`fews` は ADF（f* で割る前）と拡散係数を決めた
/// ノードの定数。ノードは格子の上の nx × ny の長方形の並びであること。
///
/// - 角の点 p の不均質な中性子束 Φ_p は `point_flux`
/// - ノード n の角 c の均質な値の見込みは h̃ = φ_x φ_y / φ̄（φ_x、φ_y は角 c を含む x と y の辺の、1ノードの問題の
///   均質な面の中性子束 = 辺の不均質な中性子束 / ADF、φ̄ はノード平均）
/// - 均質な角の値を s_p h̃ とし、core と同じ展開（参照のノード平均と辺の正味の中性子流、格子の k）と角の点の
///   バランス（`corner_balance`。真空の辺の流れは既知の r Φ_p）で s_p を解く。s_p が正でなければエラー
/// - CDF = Φ_p / (s_p h̃)。反射体のノードの、角の点が燃料（核分裂のある）のノードにも接する角は、それらのノードの
///   CDF の調和平均で割り、`normalized_corners` に印を付ける
/// - 形状関数 f_g（セル）= セルの少数群 g の核分裂率の平均（解析的な体積）/（Σf,g × 再構成した均質な中性子束の
///   セルの平均）。核分裂のない群とセルは 0
#[allow(clippy::too_many_arguments)]
fn node_pins(
    solution: &Solution,
    lattice: &Lattice,
    materials: &[&Material],
    group_of: &[usize],
    nodes: &[Node],
    neighbors: &[[Option<usize>; 4]],
    references: &[Reference],
    fews: &[FewGroup],
) -> Result<Vec<Pins>, HomogenizeError> {
    let error = |message: String| HomogenizeError::Pins(message);
    let coarse = references.first().map_or(0, |r| r.average.len());
    if coarse != 2 {
        return Err(error(format!(
            "corner discontinuity factors need 2 few groups, got {coarse}"
        )));
    }
    let (widths, heights) = (lattice.widths(), lattice.heights());
    let columns = lattice.columns();
    let NodeGrid {
        column_cuts,
        row_cuts,
        nx,
        ny,
        grid,
    } = node_grid(nodes).map_err(error)?;
    // 角の点の座標（-x、-y の辺から）
    let x_at: Vec<f64> = column_cuts
        .iter()
        .map(|&c| widths[..c].iter().sum())
        .collect();
    let height: f64 = heights.iter().sum();
    let y_at: Vec<f64> = row_cuts
        .iter()
        .rev()
        .map(|&r| height - heights[..r].iter().sum::<f64>())
        .collect();
    let points = (nx + 1) * (ny + 1);
    let phi: Vec<[f64; 2]> = (0..points)
        .map(|p| {
            let f = point_flux(
                solution,
                lattice,
                group_of,
                coarse,
                x_at[p % (nx + 1)],
                y_at[p / (nx + 1)],
            );
            [f[0], f[1]]
        })
        .collect();
    // [ノード] の並びの位置 (i, j)。どのノードも並びにちょうど 1 回ある（上で確かめた）
    let mut position = vec![(0, 0); nodes.len()];
    for (k, &n) in grid.iter().enumerate() {
        position[n] = (k % nx, k / nx);
    }
    // ノード n の角 c（-x -y、+x -y、+x +y、-x +y）の角の点
    let point = |n: usize, c: usize| {
        let (i, j) = position[n];
        let (ci, cj) = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)][c];
        ci + (nx + 1) * cj
    };
    // [ノード][角][群] の見込み h̃。辺（左、右、上、下）の均質な面の中性子束は不均質な中性子束 / ADF
    let estimates: Vec<[[f64; 2]; 4]> = references
        .iter()
        .zip(fews)
        .map(|(r, f)| {
            let face = |side: usize, g: usize| r.heterogeneous[side][g] / f.adf[side][g];
            // 角の x の辺と y の辺
            let corner_sides = [(0, 3), (1, 3), (1, 2), (0, 2)];
            std::array::from_fn(|c| {
                let (xs, ys) = corner_sides[c];
                std::array::from_fn(|g| face(xs, g) * face(ys, g) / r.average[g])
            })
        })
        .collect();
    let rule = gauss_legendre(FACE_POINTS).map_err(|e| error(e.to_string()))?;
    fn constants(f: &FewGroup) -> Constants<'_> {
        Constants {
            total: &f.total,
            scatter: &f.scatter,
            nu_fission: &f.nu_fission,
            chi: &f.chi,
            diffusion: &f.diffusion,
        }
    }
    let models: Vec<NodeModel> = references
        .iter()
        .zip(fews)
        .enumerate()
        .map(|(n, (r, f))| {
            NodeModel::new(
                &constants(f),
                solution.k,
                [0.5 * r.size[0], 0.5 * r.size[1]],
                &rule,
            )
            .map_err(|e| error(format!("node {n}: {e}")))
        })
        .collect::<Result<_, _>>()?;
    let conditions: Vec<Conditions> = references
        .iter()
        .map(|r| {
            // 辺（左、右、上、下）の外向きの流れを、面（-x、+x、-y、+y）の座標の正の向きの流れにする
            let current =
                |side: usize, sign: f64| [sign * r.currents[side][0], sign * r.currents[side][1]];
            Conditions {
                half: [0.5 * r.size[0], 0.5 * r.size[1]],
                average: [r.average[0], r.average[1]],
                corners: [[0.0; 2]; 4],
                currents: [
                    current(0, -1.0),
                    current(1, 1.0),
                    current(3, -1.0),
                    current(2, 1.0),
                ],
                axial_leakage: [0.0; 2],
            }
        })
        .collect();
    let corner_nodes: Vec<CornerNode> = grid
        .iter()
        .map(|&n| CornerNode {
            model: &models[n],
            conditions: &conditions[n],
            diffusion: [fews[n].diffusion[0], fews[n].diffusion[1]],
            weights: estimates[n],
            // 外周の真空の辺（面 -x、+x、-y、+y は辺 0、1、3、2）の流れ r Φ。鏡面反射は 0
            outflow: std::array::from_fn(|face| {
                let side = FACES[face];
                match neighbors[n][side] {
                    None => fews[n].vacuum_ratio[side].as_ref().map(|r| [r[0], r[1]]),
                    Some(_) => None,
                }
            }),
        })
        .collect();
    let scales =
        corner_balance(nx, ny, &corner_nodes, &vec![false; points], Some(&phi)).map_err(error)?;
    if let Some(p) = scales
        .iter()
        .position(|s| s.iter().any(|v| v.is_nan() || *v <= 0.0))
    {
        return Err(error(format!(
            "corner point {p}: the homogeneous corner flux scale {:?} is not positive",
            scales[p]
        )));
    }
    // [ノード][角][群] の CDF
    let cdf: Vec<[[f64; 2]; 4]> = (0..nodes.len())
        .map(|n| {
            std::array::from_fn(|c| {
                let p = point(n, c);
                std::array::from_fn(|g| phi[p][g] / (scales[p][g] * estimates[n][c][g]))
            })
        })
        .collect();
    let volumes = lattice.volumes();
    let cell_rule = gauss_legendre(CELL_POINTS).map_err(|e| error(e.to_string()))?;
    let mut out = Vec::with_capacity(nodes.len());
    for (n, node) in nodes.iter().enumerate() {
        let mut corners = Vec::with_capacity(4);
        let mut normalized_corners = vec![false; 4];
        for c in 0..4 {
            let p = point(n, c);
            // 角の点を囲むノードのうち燃料（核分裂のある）のノードの、その点の CDF。core は、どのノードが反射体
            // だったかを知らず、核分裂の有無で燃料を見分ける
            let (ci, cj) = (p % (nx + 1), p / (nx + 1));
            let mut fuel: Vec<(usize, usize)> = Vec::with_capacity(4);
            for (k, d) in nodes_at_point(nx, ny, ci, cj) {
                let m = grid[k];
                if is_fissile(&fews[m].nu_fission).map_err(|e| error(format!("node {m}: {e}")))? {
                    fuel.push((m, d));
                }
            }
            let harmonic = |g: usize| {
                let factors: Vec<f64> = fuel.iter().map(|&(m, d)| cdf[m][d][g]).collect();
                fuel_harmonic_mean(&factors)
            };
            match (node.reflector, harmonic(0), harmonic(1)) {
                (true, Some(h0), Some(h1)) => {
                    corners.push(vec![cdf[n][c][0] / h0, cdf[n][c][1] / h1]);
                    normalized_corners[c] = true;
                }
                _ => corners.push(cdf[n][c].to_vec()),
            }
        }
        // 形状関数: ノードの均質な解の角の値は s_p h̃
        let mut c = conditions[n].clone();
        c.corners = std::array::from_fn(|k| {
            let p = point(n, k);
            std::array::from_fn(|g| scales[p][g] * estimates[n][k][g])
        });
        let flux = models[n].solve(&c);
        let dx: Vec<f64> = widths[node.columns.clone()].to_vec();
        let dy: Vec<f64> = heights[node.rows.clone()].iter().rev().copied().collect();
        let [a, b] = c.half;
        let mut form = vec![vec![vec![0.0; dx.len()]; dy.len()]; 2];
        let mut y0 = -b;
        for (j, h) in dy.iter().enumerate() {
            let mut x0 = -a;
            for (i, w) in dx.iter().enumerate() {
                // 行は -y から（格子の行は上から）
                let cell = (node.rows.end - 1 - j) * columns + node.columns.start + i;
                let (mut rate, mut area) = ([0.0; 2], 0.0);
                for region in lattice.cell_regions(cell) {
                    area += volumes[region];
                    for (g, &big) in group_of.iter().enumerate() {
                        rate[big] += volumes[region]
                            * solution.flux[region][g]
                            * materials[region].fission[g];
                    }
                }
                let average = flux.average(x0, x0 + w, y0, y0 + h, &cell_rule);
                for g in 0..2 {
                    let rate = rate[g] / area;
                    if rate == 0.0 {
                        continue;
                    }
                    let denominator = fews[n].fission[g] * average[g];
                    if denominator.is_nan() || denominator <= 0.0 {
                        return Err(error(format!(
                            "node {n} cell ({j}, {i}) group {g}: Σf × homogeneous flux is {denominator}"
                        )));
                    }
                    form[g][j][i] = rate / denominator;
                }
                x0 += w;
            }
            y0 += h;
        }
        out.push(Pins {
            dx,
            dy,
            form,
            corners,
            normalized_corners,
            halves: vec![false; 4],
        });
    }
    Ok(out)
}

/// 集合体を n × n に分けた部分ノードの少数群定数（ADR-0108）。並びは -y の行から、各行は -x から。外周がすべて
/// 鏡面反射の単一の格子でだけ作る（均質な解がどう分けても平らなので、均質な面の中性子束はノード平均）。
///
/// - 断面積、拡散係数、k∞ は集合体の定数 `few` のもの
/// - 辺の ADF: 集合体の外周の辺では、その辺を n 個に分けた区間の不均質な中性子束 / ノード平均の中性子束（区間に
///   重なる格子のセルの面の中性子束を、重なる長さで平均する。区間の端がセルの面の途中にあれば、その面の中性子束は
///   面の上で一様とする）。内側の辺は 1
/// - 形状関数（`few.pins` があれば）: 部分ノードに重なるセルに分け、分ける線にかかるセルは半分ずつにして元の値を
///   持たせ、その辺に `halves` の印を付ける
/// - 角の不連続因子: 角の点に接する格子の中の領域（角の点を頂点とする象限のうち格子の中にあるもの）の中性子束の
///   平均 / ノード平均の中性子束（トラックの体積）
#[allow(clippy::too_many_arguments)]
pub fn split_constants(
    few: &FewGroup,
    solution: &Solution,
    lattice: &Lattice,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
    boundaries: &Boundaries,
    n: usize,
) -> Result<Vec<FewGroup>, HomogenizeError> {
    if n == 0 {
        return Err(HomogenizeError::Nodes("split must be positive".into()));
    }
    if !reflective(boundaries) {
        return Err(HomogenizeError::Nodes(
            "splits need a lattice with every side reflective (other lattices: ADR-0108)".into(),
        ));
    }
    check_regions(solution, materials)?;
    let group_of = coarse_groups(fine, coarse_bounds)?;
    let coarse = coarse_bounds.len() - 1;
    let volume: f64 = solution.volumes.iter().sum();
    let mut average = vec![0.0; coarse];
    for (phi, v) in solution.flux.iter().zip(&solution.volumes) {
        for (g, &big) in group_of.iter().enumerate() {
            average[big] += v * phi[g] / volume;
        }
    }
    if average.iter().any(|a| a.is_nan() || *a <= 0.0) {
        return Err(HomogenizeError::NoFlux);
    }
    let (widths, heights) = (lattice.widths(), lattice.heights());
    let (columns, rows) = (lattice.columns(), lattice.rows());
    let (width, height) = (widths.iter().sum::<f64>(), heights.iter().sum::<f64>());
    // 格子のセルの境界（-x、-y の辺から。行は下から）
    let running = |sizes: &mut dyn Iterator<Item = f64>| -> Vec<f64> {
        let mut e = vec![0.0];
        let mut at = 0.0;
        for s in sizes {
            at += s;
            e.push(at);
        }
        e
    };
    let x_edges = running(&mut widths.iter().copied());
    let y_edges = running(&mut heights.iter().rev().copied());
    let y_sizes: Vec<f64> = heights.iter().rev().copied().collect();
    let face_flux = |face: usize| -> Vec<f64> {
        let mut out = vec![0.0; coarse];
        for (g, &big) in group_of.iter().enumerate() {
            out[big] += solution.face_fluxes[face][g];
        }
        out
    };
    // 辺 side（左、右、上、下）の区間 q の不均質な中性子束
    let segment = |side: usize, q: usize| -> Vec<f64> {
        let along = if side < 2 { height } else { width };
        let (lo, hi) = (
            along * q as f64 / n as f64,
            along * (q + 1) as f64 / n as f64,
        );
        let mut sum = vec![0.0; coarse];
        let cells = if side < 2 { rows } else { columns };
        for k in 0..cells {
            // 辺に沿った位置（-y か -x から）の k 番目のセル
            let (a, b) = if side < 2 {
                (y_edges[k], y_edges[k + 1])
            } else {
                (x_edges[k], x_edges[k + 1])
            };
            let overlap = (b.min(hi) - a.max(lo)).max(0.0);
            if overlap == 0.0 {
                continue;
            }
            let face = match side {
                // 左と右の辺は縦の面。格子の行は上から
                0 => crate::cmfd::vertical_face(columns, rows - 1 - k, 0),
                1 => crate::cmfd::vertical_face(columns, rows - 1 - k, columns),
                2 => crate::cmfd::horizontal_face(rows, columns, 0, k),
                _ => crate::cmfd::horizontal_face(rows, columns, rows, k),
            };
            for (s, f) in sum.iter_mut().zip(face_flux(face)) {
                *s += overlap * f;
            }
        }
        sum.iter().map(|s| s / (hi - lo)).collect()
    };
    // 点 (x, y)（-x、-y の辺から）の不均質な中性子束
    let point = |x: f64, y: f64| point_flux(solution, lattice, &group_of, coarse, x, y);
    let tolerance = 1e-9 * width.max(height);
    // 区間 [lo, hi] に重なるセルの (番号, 重なる大きさ) と、両端がセルの途中か
    // 丸ごと入るセルは元の大きさをそのまま使う（境界の座標の引き算の丸めを入れない）
    let cut = |sizes: &[f64], edges: &[f64], lo: f64, hi: f64| -> (Vec<(usize, f64)>, bool, bool) {
        let cells: Vec<(usize, f64)> = (0..edges.len() - 1)
            .filter_map(|k| {
                let whole = edges[k] >= lo - tolerance && edges[k + 1] <= hi + tolerance;
                let overlap = if whole {
                    sizes[k]
                } else {
                    edges[k + 1].min(hi) - edges[k].max(lo)
                };
                (overlap > tolerance).then_some((k, overlap))
            })
            .collect();
        let inside = |at: f64| edges.iter().all(|e| (e - at).abs() > tolerance);
        let at_lo = lo > tolerance && inside(lo);
        let at_hi = hi < edges[edges.len() - 1] - tolerance && inside(hi);
        (cells, at_lo, at_hi)
    };
    let mut out = Vec::with_capacity(n * n);
    for j in 0..n {
        for i in 0..n {
            let mut sub = few.clone();
            // 辺（左、右、上、下）: 集合体の外周なら区間の値、内側なら 1
            let outer = [i == 0, i == n - 1, j == n - 1, j == 0];
            sub.adf = (0..4)
                .map(|side| {
                    if outer[side] {
                        let q = if side < 2 { j } else { i };
                        segment(side, q)
                            .iter()
                            .zip(&average)
                            .map(|(s, a)| s / a)
                            .collect()
                    } else {
                        vec![1.0; coarse]
                    }
                })
                .collect();
            if let Some(pins) = &few.pins {
                let (x0, x1) = (
                    width * i as f64 / n as f64,
                    width * (i + 1) as f64 / n as f64,
                );
                let (y0, y1) = (
                    height * j as f64 / n as f64,
                    height * (j + 1) as f64 / n as f64,
                );
                let (cols, cut_left, cut_right) = cut(widths, &x_edges, x0, x1);
                let (rws, cut_bottom, cut_top) = cut(&y_sizes, &y_edges, y0, y1);
                // 集合体の外周の角は集合体の CDF をそのまま使う（同じ量を計算し直さない）
                let corner = |x: f64, y: f64| -> Vec<f64> {
                    let at = |v: f64, end: f64| (v - end).abs() <= tolerance;
                    let outer = [(0.0, 0.0), (width, 0.0), (width, height), (0.0, height)]
                        .iter()
                        .position(|&(cx, cy)| at(x, cx) && at(y, cy));
                    match outer {
                        Some(c) => pins.corners[c].clone(),
                        None => point(x, y)
                            .iter()
                            .zip(&average)
                            .map(|(p, a)| p / a)
                            .collect(),
                    }
                };
                sub.pins = Some(Pins {
                    dx: cols.iter().map(|&(_, w)| w).collect(),
                    dy: rws.iter().map(|&(_, h)| h).collect(),
                    form: pins
                        .form
                        .iter()
                        .map(|g| {
                            rws.iter()
                                .map(|&(r, _)| cols.iter().map(|&(c, _)| g[r][c]).collect())
                                .collect()
                        })
                        .collect(),
                    corners: vec![
                        corner(x0, y0),
                        corner(x1, y0),
                        corner(x1, y1),
                        corner(x0, y1),
                    ],
                    normalized_corners: vec![false; 4],
                    halves: vec![cut_left, cut_right, cut_top, cut_bottom],
                });
            }
            out.push(sub);
        }
    }
    Ok(out)
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
        pins: None,
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
            false,
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
            false,
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
                false,
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
            false,
        );
        assert!(
            matches!(misaligned, Err(HomogenizeError::Nodes(_))),
            "{misaligned:?}"
        );
    }

    /// `nodes`（x が速く、-y の行から並べた nx × ny の格子）の定数と形状（ADR-0125）で、core の形の角の点の釣り合い
    /// （未知数が角の点の不均質な中性子束、均質な角の値は Φ / CDF）を解くと、格子の角の点の不均質な中性子束（角の点に
    /// 接する領域の中性子束の平均）に戻り、再構成したセルの出力（Σf × 均質な中性子束のセルの平均 × 形状関数）が
    /// MOC のセルの核分裂率に戻ることを確かめる。ノード平均と面の正味の中性子流は、core の NEM が GET で再現する参照の値
    fn assert_pins_reproduce(
        lattice: &Lattice,
        solution: &Solution,
        materials: &[Material],
        boundaries: &Boundaries,
        nodes: &[Node],
        nx: usize,
        ny: usize,
    ) {
        use xlc_numeric::pins::FACE_POINTS;
        use xlc_numeric::pins::corners::{CornerNode, corner_balance};
        use xlc_numeric::pins::node::{Conditions, Constants, NodeModel};
        let refs: Vec<&Material> = materials.iter().collect();
        let fine = fine();
        let fgs = homogenize_nodes(
            solution,
            lattice,
            &refs,
            &fine,
            fine.bounds(),
            boundaries,
            nodes,
            true,
        )
        .unwrap();
        let references = node_references(
            solution,
            lattice,
            &refs,
            &fine,
            fine.bounds(),
            boundaries,
            nodes,
        )
        .unwrap()
        .references;
        let rule = xlc_numeric::quadrature::gauss_legendre(FACE_POINTS).unwrap();
        let cells_rule = xlc_numeric::quadrature::gauss_legendre(6).unwrap();
        let pins: Vec<&Pins> = fgs.iter().map(|f| f.pins.as_ref().unwrap()).collect();
        // ノード n の角 c の角の点（x が速い (nx + 1) × (ny + 1)）
        let point = |n: usize, c: usize| {
            let (i, j) = (n % nx, n / nx);
            let (ci, cj) = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)][c];
            ci + (nx + 1) * cj
        };
        // 炉心の CDF: 燃料の CDF で割った角は、角の点に接する燃料（核分裂のある）のノードの CDF の調和平均を掛ける
        let refs_fissile = |f: &FewGroup| f.nu_fission.iter().any(|&v| v > 0.0);
        let cdf = |n: usize, c: usize, g: usize| -> f64 {
            if !pins[n].normalized_corners[c] {
                return pins[n].corners[c][g];
            }
            let p = point(n, c);
            let fuel: Vec<f64> = (0..nodes.len())
                .filter(|&m| refs_fissile(&fgs[m]))
                .flat_map(|m| {
                    (0..4)
                        .filter(move |&d| point(m, d) == p)
                        .map(move |d| (m, d))
                })
                .map(|(m, d)| pins[m].corners[d][g])
                .collect();
            assert!(!fuel.is_empty(), "node {n} corner {c}");
            let harmonic = fuel.len() as f64 / fuel.iter().map(|f| 1.0 / f).sum::<f64>();
            pins[n].corners[c][g] * harmonic
        };
        let models: Vec<NodeModel> = references
            .iter()
            .zip(&fgs)
            .map(|(r, f)| {
                let constants = Constants {
                    total: &f.total,
                    scatter: &f.scatter,
                    nu_fission: &f.nu_fission,
                    chi: &f.chi,
                    diffusion: &f.diffusion,
                };
                NodeModel::new(
                    &constants,
                    solution.k,
                    [0.5 * r.size[0], 0.5 * r.size[1]],
                    &rule,
                )
                .unwrap()
            })
            .collect();
        let conditions: Vec<Conditions> = references
            .iter()
            .map(|r| {
                // 辺（左、右、上、下）の外向きの流れを、面（-x、+x、-y、+y）の座標の正の向きの流れにする
                let current = |side: usize, sign: f64| {
                    [sign * r.currents[side][0], sign * r.currents[side][1]]
                };
                Conditions {
                    half: [0.5 * r.size[0], 0.5 * r.size[1]],
                    average: [r.average[0], r.average[1]],
                    corners: [[0.0; 2]; 4],
                    currents: [
                        current(0, -1.0),
                        current(1, 1.0),
                        current(3, -1.0),
                        current(2, 1.0),
                    ],
                    axial_leakage: [0.0; 2],
                }
            })
            .collect();
        let corner_nodes: Vec<CornerNode> = (0..nodes.len())
            .map(|n| CornerNode {
                model: &models[n],
                conditions: &conditions[n],
                diffusion: [fgs[n].diffusion[0], fgs[n].diffusion[1]],
                weights: std::array::from_fn(|c| std::array::from_fn(|g| 1.0 / cdf(n, c, g))),
                outflow: std::array::from_fn(|face| {
                    fgs[n].vacuum_ratio[[0, 1, 3, 2][face]]
                        .as_ref()
                        .map(|r| [r[0], r[1]])
                }),
            })
            .collect();
        let zero = vec![false; (nx + 1) * (ny + 1)];
        let phi = corner_balance(nx, ny, &corner_nodes, &zero, None).unwrap();
        // 角の点の不均質な中性子束: 角の点に接する格子の中の領域の中性子束の平均（少数群は 2 群のまま）
        let (width, height): (f64, f64) = (
            lattice.widths().iter().sum(),
            lattice.heights().iter().sum(),
        );
        let x_edges: Vec<f64> = (0..=nx)
            .map(|i| nodes[..i].iter().map(|m| r_width(lattice, m)).sum())
            .collect();
        let y_edges: Vec<f64> = (0..=ny)
            .map(|j| (0..j).map(|jj| r_height(lattice, &nodes[jj * nx])).sum())
            .collect();
        for cj in 0..=ny {
            for ci in 0..=nx {
                let (x, y) = (x_edges[ci], y_edges[cj]);
                let e = 1e-9 * width.max(height);
                let mut want = [0.0; 2];
                let mut count = 0.0;
                for (sx, sy) in [(1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0)] {
                    let (px, py) = (x + sx * e, y + sy * e);
                    if px <= 0.0 || px >= width || py <= 0.0 || py >= height {
                        continue;
                    }
                    let region = lattice.region_at(px - 0.5 * width, py - 0.5 * height);
                    for (w, f) in want.iter_mut().zip(&solution.flux[region]) {
                        *w += f;
                    }
                    count += 1.0;
                }
                let got = phi[ci + (nx + 1) * cj];
                for g in 0..2 {
                    let want = want[g] / count;
                    assert!(
                        (got[g] - want).abs() < 1e-9 * want,
                        "corner point ({ci}, {cj}) group {g}: {} vs {want}",
                        got[g]
                    );
                }
            }
        }
        // セルの出力が MOC のセルの核分裂率（解析的な体積で平均）に戻る
        let volumes = lattice.volumes();
        for (n, node) in nodes.iter().enumerate() {
            let mut c = conditions[n].clone();
            c.corners = std::array::from_fn(|k| {
                let p = phi[point(n, k)];
                std::array::from_fn(|g| p[g] / cdf(n, k, g))
            });
            let flux = models[n].solve(&c);
            let [a, b] = c.half;
            let mut y0 = -b;
            for (j, dy) in pins[n].dy.iter().enumerate() {
                let mut x0 = -a;
                for (i, dx) in pins[n].dx.iter().enumerate() {
                    let average = flux.average(x0, x0 + dx, y0, y0 + dy, &cells_rule);
                    // 行は -y から（格子の行は上から）
                    let row = node.rows.end - 1 - j;
                    let cell = row * lattice.columns() + node.columns.start + i;
                    let (mut rate, mut area) = ([0.0; 2], 0.0);
                    for region in lattice.cell_regions(cell) {
                        area += volumes[region];
                        for (g, r) in rate.iter_mut().enumerate() {
                            *r += volumes[region]
                                * solution.flux[region][g]
                                * refs[region].fission[g];
                        }
                    }
                    for g in 0..2 {
                        let want = rate[g] / area;
                        let got = fgs[n].fission[g] * average[g] * pins[n].form[g][j][i];
                        assert!(
                            (got - want).abs() < 1e-9 * want.max(1e-12),
                            "node {n} cell ({j}, {i}) group {g}: {got} vs {want}"
                        );
                    }
                    x0 += dx;
                }
                y0 += dy;
            }
        }
    }

    /// ノードの幅と高さ
    fn r_width(lattice: &Lattice, node: &Node) -> f64 {
        lattice.widths()[node.columns.clone()].iter().sum()
    }

    fn r_height(lattice: &Lattice, node: &Node) -> f64 {
        lattice.heights()[node.rows.clone()].iter().sum()
    }

    #[test]
    fn strip_pins_reproduce_the_lattice() {
        // 燃料と水の帯（右が真空）を燃料と反射体のノードに分けた定数と形状は、定数を作った体系の角の中性子束と
        // セルの核分裂率を再現する（ADR-0125）
        let (lattice, solution, materials, boundaries) = strip();
        let nodes = [
            node(0..1, 0..1),
            Node {
                reflector: true,
                ..node(1..2, 0..1)
            },
        ];
        assert_pins_reproduce(&lattice, &solution, &materials, &boundaries, &nodes, 2, 1);
    }

    /// 燃料のピン 2 × 2 の右と上に水のセルを 1 列ずつ並べた 3 × 3 と、外周の境界 `boundaries` の解。燃料は右と上の
    /// 辺に接しない
    fn l_shaped(boundaries: Boundaries) -> (Lattice, Solution, Vec<Material>, Boundaries) {
        let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let cell = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::rectilinear(
            vec![1.26; 3],
            vec![1.26; 3],
            vec![Cell::Pin(pin), Cell::Pin(cell)],
            vec![vec![1, 1, 1], vec![0, 0, 1], vec![0, 0, 1]],
        )
        .unwrap();
        let tracking = Tracking::with_boundaries(&lattice, 32, 0.02, boundaries).unwrap();
        // ピンの中心の領域だけ燃料、ほかは水
        let materials: Vec<Material> = (0..lattice.cells())
            .flat_map(|cell| {
                let fuelled = lattice
                    .cell_type(cell)
                    .pin()
                    .is_some_and(|p| !p.radii().is_empty());
                lattice
                    .cell_regions(cell)
                    .enumerate()
                    .map(|(k, _)| {
                        if fuelled && k == 0 {
                            fuel.clone()
                        } else {
                            water.clone()
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        let refs: Vec<&Material> = materials.iter().collect();
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        (lattice, solution, materials, boundaries)
    }

    #[test]
    fn l_shaped_pins_reproduce_the_lattice() {
        // 燃料のノード（2 × 2 のピン）の上と右に反射体のノード、右上に水のノード（反射体の印なし。燃料に角だけで
        // 接する）。反射体の角は燃料の面と、右上のノードに接する（ADR-0125）
        let (lattice, solution, materials, boundaries) = l_shaped(vacuum());
        let reflector = |n: Node| Node {
            reflector: true,
            ..n
        };
        // x が速く、-y の行（格子の行 1、2）から
        let nodes = [
            node(0..2, 1..3),
            reflector(node(2..3, 1..3)),
            reflector(node(0..2, 0..1)),
            node(2..3, 0..1),
        ];
        assert_pins_reproduce(&lattice, &solution, &materials, &boundaries, &nodes, 2, 2);
    }

    #[test]
    fn non_positive_corner_scales_are_errors() {
        // 右と上の 2 辺が真空の格子では、右上の角の点の均質な値の大きさ s_p の熱群が負になる（#0065）。正でない
        // s_p はエラー（ADR-0125）
        let (lattice, solution, materials, boundaries) = l_shaped(Boundaries {
            top: Boundary::Vacuum,
            ..vacuum()
        });
        let refs: Vec<&Material> = materials.iter().collect();
        let fine = fine();
        let reflector = |n: Node| Node {
            reflector: true,
            ..n
        };
        let nodes = [
            node(0..2, 1..3),
            reflector(node(2..3, 1..3)),
            reflector(node(0..2, 0..1)),
            node(2..3, 0..1),
        ];
        let err = homogenize_nodes(
            &solution,
            &lattice,
            &refs,
            &fine,
            fine.bounds(),
            &boundaries,
            &nodes,
            true,
        )
        .unwrap_err();
        assert!(
            matches!(&err, HomogenizeError::Pins(m) if m.contains("not positive")),
            "{err:?}"
        );
    }

    #[test]
    fn pins_need_reflector_nodes_without_fission() {
        // 燃料の CDF で割る角は、核分裂のあるノードを燃料とする（ADR-0125）。反射体の印のノードに核分裂があると、
        // 自分の CDF で割ることになり、core で掛け戻す値と食い違うので、形状関数を作るときはエラー
        let (lattice, solution, materials, boundaries) = strip();
        let refs: Vec<&Material> = materials.iter().collect();
        let fine = fine();
        let nodes = [
            Node {
                reflector: true,
                ..node(0..1, 0..1)
            },
            node(1..2, 0..1),
        ];
        let err = homogenize_nodes(
            &solution,
            &lattice,
            &refs,
            &fine,
            fine.bounds(),
            &boundaries,
            &nodes,
            true,
        )
        .unwrap_err();
        assert!(
            matches!(&err, HomogenizeError::Nodes(m) if m.contains("fission")),
            "{err:?}"
        );
    }

    #[test]
    fn strip_corner_factors_follow_the_face_discontinuity() {
        // 帯（上と下は鏡面反射）では均質な解が x にしか変わらず、境の角の均質な値の見込みは境の面の均質な中性子束に
        // なる。反射体の燃料に接する角の C* は、境の面の f* に一致する（ADR-0125）
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
            true,
        )
        .unwrap();
        let reflector = &fgs[1];
        let pins = reflector.pins.as_ref().unwrap();
        // 燃料に接する角は -x -y と -x +y
        assert_eq!(pins.normalized_corners, vec![true, false, false, true]);
        assert!(
            !fgs[0]
                .pins
                .as_ref()
                .unwrap()
                .normalized_corners
                .contains(&true)
        );
        for c in [0, 3] {
            for g in 0..2 {
                let (got, want) = (pins.corners[c][g], reflector.adf[0][g]);
                assert!(
                    (got - want).abs() < 1e-6 * want,
                    "corner {c} group {g}: {got} vs {want}"
                );
            }
        }
        // 反射体の形状関数は 0
        assert!(pins.form.iter().flatten().flatten().all(|&f| f == 0.0));
    }

    #[test]
    fn one_reflective_node_keeps_the_assembly_pins() {
        // 外周がすべて鏡面反射の格子を1つのノードにすると、均質な解は平らで、CDF は角の不均質な中性子束 /
        // ノード平均（`pin_shapes`）に一致する。形状関数の分母は Σf × 均質な中性子束（トラックの体積で平均した
        // ノード平均）で、`pin_shapes` の分母（解析的な体積で平均した核分裂率）との比だけ違う
        let (lattice, solution, materials) = margined(0.04);
        let refs: Vec<&Material> = materials.iter().collect();
        let fine = fine();
        let assembly = pin_shapes(
            &solution,
            &lattice,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        let fgs = homogenize_nodes(
            &solution,
            &lattice,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
            &[node(0..4, 0..4)],
            true,
        )
        .unwrap();
        let pins = fgs[0].pins.as_ref().unwrap();
        assert_eq!((&pins.dx, &pins.dy), (&assembly.dx, &assembly.dy));
        // 平らさは MOC の鏡面反射の辺の正味の中性子流（収束の許容差の程度で 0）の程度
        let flat = |got: f64, want: f64| (got - want).abs() < 1e-6 * want;
        for (got, want) in pins
            .corners
            .iter()
            .flatten()
            .zip(assembly.corners.iter().flatten())
        {
            assert!(flat(*got, *want), "{got} vs {want}");
        }
        let volumes = lattice.volumes();
        let area: f64 = volumes.iter().sum();
        let (track_volume, track) = (solution.volumes.iter().sum::<f64>(), &solution.volumes);
        for g in 0..2 {
            let analytic = (0..refs.len())
                .map(|r| volumes[r] * solution.flux[r][g] * refs[r].fission[g])
                .sum::<f64>()
                / area;
            let average = (0..refs.len())
                .map(|r| track[r] * solution.flux[r][g])
                .sum::<f64>()
                / track_volume;
            let ratio = analytic / (fgs[0].fission[g] * average);
            for (got, want) in pins.form[g]
                .iter()
                .flatten()
                .zip(assembly.form[g].iter().flatten())
            {
                assert!(
                    (got - want * ratio).abs() < 1e-6 * want.max(1e-3),
                    "group {g}: {got} vs {want} × {ratio}"
                );
            }
        }
        assert_eq!(pins.normalized_corners, vec![false; 4]);
    }

    /// 2×2 のピン（対角に燃料）の周りに幅 0.04 の余白を付けた格子、領域ごとの物質、鏡面反射の解
    fn margined(margin: f64) -> (Lattice, Solution, Vec<Material>) {
        let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let cell = PinCell::new(1.26, vec![]).unwrap();
        let lattice =
            Lattice::with_margin(vec![pin, cell], vec![vec![0, 1], vec![1, 0]], margin).unwrap();
        // ピンの中心の領域だけ燃料、ほかは水
        let materials: Vec<Material> = (0..lattice.cells())
            .flat_map(|cell| {
                let fuelled = lattice
                    .cell_type(cell)
                    .pin()
                    .is_some_and(|p| !p.radii().is_empty());
                lattice
                    .cell_regions(cell)
                    .enumerate()
                    .map(|(k, _)| {
                        if fuelled && k == 0 {
                            fuel.clone()
                        } else {
                            water.clone()
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let refs: Vec<&Material> = materials.iter().collect();
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        (lattice, solution, materials)
    }

    #[test]
    fn pin_shapes_average_to_one_and_take_corner_fluxes() {
        // 形状関数はセルの群の核分裂率の平均 / ノードの群の核分裂率の平均で、面積で重みを付けたノードの平均は 1。
        // 角の不連続因子は、余白の角のセルの中性子束 / ノード平均の中性子束（ADR-0094、ADR-0096）。形状関数は
        // 解析的な体積で、角の不連続因子はトラックの体積で平均する
        let (lattice, solution, materials) = margined(0.04);
        let refs: Vec<&Material> = materials.iter().collect();
        let pins = pin_shapes(
            &solution,
            &lattice,
            &refs,
            &fine(),
            fine().bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        let close = |a: f64, b: f64| (a - b).abs() < 1e-12 * b.abs().max(1.0);
        // 対角の対称は、MOC の収束の程度で成り立つ
        let symmetric = |a: f64, b: f64| (a - b).abs() < 1e-6 * b.abs();
        // セルの大きさは格子の列の幅と行の高さ（行は下から）
        assert_eq!(pins.dx, vec![0.04, 1.26, 1.26, 0.04]);
        assert_eq!(pins.dy, vec![0.04, 1.26, 1.26, 0.04]);
        let area = |j: usize, i: usize| pins.dx[i] * pins.dy[j];
        for g in 0..2 {
            let form = &pins.form[g];
            assert_eq!((form.len(), form[0].len()), (4, 4));
            let mean: f64 = (0..4)
                .flat_map(|j| (0..4).map(move |i| (j, i)))
                .map(|(j, i)| area(j, i) * form[j][i])
                .sum::<f64>()
                / (2.6 * 2.6);
            assert!(close(mean, 1.0), "group {g}: {mean}");
            // 燃料は上の行（下から 2 行目）の左と、下の行（下から 1 行目）の右。対称なので等しい
            assert!(
                form[2][1] > 1.0 && symmetric(form[2][1], form[1][2]),
                "{form:?}"
            );
            for (j, i) in [
                (2, 2),
                (1, 1),
                (0, 0),
                (0, 3),
                (3, 0),
                (3, 3),
                (0, 1),
                (2, 0),
            ] {
                assert_eq!(form[j][i], 0.0, "cell ({j}, {i})");
            }
        }
        // 角のセル（格子の行は上から）: -x -y は左下、+x -y は右下、+x +y は右上、-x +y は左上
        let corner_cells = [(3, 0), (3, 3), (0, 3), (0, 0)];
        // 面の ADF と同じく、トラックから推定した体積で平均する
        let volumes = &solution.volumes;
        let volume: f64 = volumes.iter().sum();
        for g in 0..2 {
            let average = solution
                .flux
                .iter()
                .zip(volumes)
                .map(|(f, v)| v * f[g])
                .sum::<f64>()
                / volume;
            for (c, (r, col)) in corner_cells.into_iter().enumerate() {
                let cell = r * 4 + col;
                let regions = lattice.cell_regions(cell);
                let (flux, v) = regions.fold((0.0, 0.0), |(f, v), i| {
                    (f + volumes[i] * solution.flux[i][g], v + volumes[i])
                });
                let expected = flux / v / average;
                assert!(close(pins.corners[c][g], expected), "corner {c} group {g}");
            }
            // 左下と右上は水のピンの角、左上と右下は燃料のピンの角で、対称
            assert!(symmetric(pins.corners[0][g], pins.corners[2][g]));
            assert!(symmetric(pins.corners[1][g], pins.corners[3][g]));
        }
    }

    #[test]
    fn pin_shapes_need_a_reflective_lattice_with_corner_margins() {
        let (lattice, solution, materials) = margined(0.04);
        let refs: Vec<&Material> = materials.iter().collect();
        let shapes = |boundaries: &Boundaries| {
            pin_shapes(
                &solution,
                &lattice,
                &refs,
                &fine(),
                fine().bounds(),
                boundaries,
            )
        };
        // 外周が鏡面反射でない格子の角の均質な中性子束は決めていない（#0020）
        let err = shapes(&vacuum()).unwrap_err().to_string();
        assert!(err.contains("reflective"), "{err}");
        // 角の中性子束は余白の角のセルで取る
        let (lattice, solution, materials) = margined(0.0);
        let refs: Vec<&Material> = materials.iter().collect();
        let err = pin_shapes(
            &solution,
            &lattice,
            &refs,
            &fine(),
            fine().bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("margin"), "{err}");
    }

    /// `layout` のピン（0 が燃料、1 が水）の周りに幅 0.04 の余白を付けた格子、物質、鏡面反射の解と、集合体の定数
    fn assembly(layout: Vec<Vec<usize>>) -> (Lattice, Solution, Vec<Material>, FewGroup) {
        let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let cell = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::with_margin(vec![pin, cell], layout, 0.04).unwrap();
        let materials: Vec<Material> = (0..lattice.cells())
            .flat_map(|cell| {
                let fuelled = lattice
                    .cell_type(cell)
                    .pin()
                    .is_some_and(|p| !p.radii().is_empty());
                lattice
                    .cell_regions(cell)
                    .enumerate()
                    .map(|(k, _)| {
                        if fuelled && k == 0 {
                            fuel.clone()
                        } else {
                            water.clone()
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let refs: Vec<&Material> = materials.iter().collect();
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        let mut few = homogenize(
            &solution,
            &refs,
            &fine(),
            fine().bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        few.pins = Some(
            pin_shapes(
                &solution,
                &lattice,
                &refs,
                &fine(),
                fine().bounds(),
                &Boundaries::REFLECTIVE,
            )
            .unwrap(),
        );
        (lattice, solution, materials, few)
    }

    fn split_of(
        (lattice, solution, materials, few): &(Lattice, Solution, Vec<Material>, FewGroup),
        n: usize,
    ) -> Result<Vec<FewGroup>, HomogenizeError> {
        let refs: Vec<&Material> = materials.iter().collect();
        split_constants(
            few,
            solution,
            lattice,
            &refs,
            &fine(),
            fine().bounds(),
            &Boundaries::REFLECTIVE,
            n,
        )
    }

    #[test]
    fn one_by_one_split_is_the_assembly() {
        // 分けなければ、集合体の定数そのもの（ADR-0108）
        let a = assembly(vec![vec![0, 1], vec![1, 0]]);
        let subs = split_of(&a, 1).unwrap();
        assert_eq!(subs.len(), 1);
        // 集合体の ADF は、MOC の外周の正味の中性子流（鏡面反射なので収束の程度で 0）を条件にした1ノードの NEM で作る。
        // 分けた定数は均質な解を平らとするので、MOC の収束の程度で一致する
        let close = |x: f64, y: f64| (x - y).abs() < 1e-8 * y.abs();
        for (s, w) in subs[0].adf.iter().flatten().zip(a.3.adf.iter().flatten()) {
            assert!(close(*s, *w), "{:?} vs {:?}", subs[0].adf, a.3.adf);
        }
        assert_eq!(subs[0].pins, a.3.pins);
        assert_eq!(subs[0].total, a.3.total);
    }

    #[test]
    fn two_by_two_split_takes_segment_fluxes() {
        // 2 × 2 のピン（左上と右下が燃料）を 2 × 2 に分ける。分ける線はピンのセルの境に乗る。外周の辺の ADF は、
        // その辺の半分の区間の面の中性子束 / ノード平均（トラックの体積）。内側の辺は 1。断面積は集合体のもの
        let a = assembly(vec![vec![0, 1], vec![1, 0]]);
        let (lattice, solution, _, few) = &a;
        let subs = split_of(&a, 2).unwrap();
        assert_eq!(subs.len(), 4);
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
        // 左の辺（縦の面の線 0）の下の区間は行 2、3（行は上から。行 3 は余白）
        let (rows, columns) = (lattice.rows(), lattice.columns());
        let heights = lattice.heights();
        let face = |row: usize| crate::cmfd::vertical_face(columns, row, 0);
        for (g, avg) in average.iter().enumerate() {
            let lower = (solution.face_fluxes[face(2)][g] * heights[2]
                + solution.face_fluxes[face(3)][g] * heights[3])
                / (heights[2] + heights[3]);
            let got = subs[0].adf[0][g];
            assert!((got - lower / avg).abs() < 1e-12 * got, "group {g}");
            // 内側の辺（部分ノード 0 の右と上）は 1
            assert_eq!((subs[0].adf[1][g], subs[0].adf[2][g]), (1.0, 1.0));
        }
        assert!(rows == 4 && columns == 4);
        for s in &subs {
            assert_eq!(
                (&s.total, &s.diffusion, s.k_inf),
                (&few.total, &few.diffusion, few.k_inf)
            );
        }
        // 部分ノード 0（左下）の形状関数は、集合体の左下の 2 × 2 のセル。分ける線はセルの境なので半分のセルはない
        let pins = subs[0].pins.as_ref().unwrap();
        let whole = few.pins.as_ref().unwrap();
        assert_eq!(pins.dx, vec![0.04, 1.26]);
        assert_eq!(pins.dy, vec![0.04, 1.26]);
        assert_eq!(pins.halves, vec![false; 4]);
        for g in 0..2 {
            assert_eq!(
                pins.form[g],
                vec![
                    whole.form[g][0][..2].to_vec(),
                    whole.form[g][1][..2].to_vec()
                ]
            );
        }
        // 集合体の外周の角（-x -y）は集合体の CDF（足し合わせの順の丸めまで）
        for (got, want) in pins.corners[0].iter().zip(&whole.corners[0]) {
            assert!((got - want).abs() < 1e-12 * want);
        }
        // 集合体の中心（部分ノード 0 の +x +y の角）は、4 つのピンのセルの角の領域の中性子束の平均
        let center: Vec<f64> = (0..2)
            .map(|g| {
                let e = 1e-7;
                [(e, e), (-e, e), (-e, -e), (e, -e)]
                    .iter()
                    .map(|&(x, y)| solution.flux[lattice.region_at(x, y)][g])
                    .sum::<f64>()
                    / 4.0
                    / average[g]
            })
            .collect();
        for (got, want) in pins.corners[2].iter().zip(&center) {
            assert!((got - want).abs() < 1e-12 * want);
        }
    }

    #[test]
    fn odd_lattices_split_through_the_middle_pins() {
        // 3 × 3 のピンを 2 × 2 に分けると、中央の行と列のピンのセルを半分ずつに分ける
        let a = assembly(vec![vec![0, 1, 0], vec![1, 0, 1], vec![0, 1, 0]]);
        let subs = split_of(&a, 2).unwrap();
        let pins = subs[0].pins.as_ref().unwrap();
        let close = |got: &[f64], want: &[f64]| {
            got.len() == want.len() && got.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-12)
        };
        assert!(close(&pins.dx, &[0.04, 1.26, 0.63]), "{:?}", pins.dx);
        assert!(close(&pins.dy, &[0.04, 1.26, 0.63]), "{:?}", pins.dy);
        assert_eq!(pins.halves, vec![false, true, true, false]);
        let whole = a.3.pins.as_ref().unwrap();
        // 半分のセルは元のセルの形状関数を持つ（中央のピンは行 2（下から）、列 2）
        assert_eq!(pins.form[0][2][2], whole.form[0][2][2]);
        // 右上の部分ノードは中央のピンの残りの半分から始まる
        let top_right = subs[3].pins.as_ref().unwrap();
        assert!(
            close(&top_right.dx, &[0.63, 1.26, 0.04]),
            "{:?}",
            top_right.dx
        );
        assert_eq!(top_right.halves, vec![true, false, false, true]);
    }

    #[test]
    fn splits_need_a_reflective_lattice() {
        let a = assembly(vec![vec![0, 1], vec![1, 0]]);
        let refs: Vec<&Material> = a.2.iter().collect();
        let err = split_constants(
            &a.3,
            &a.1,
            &a.0,
            &refs,
            &fine(),
            fine().bounds(),
            &vacuum(),
            2,
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("reflective"), "{err}");
        assert!(split_of(&a, 0).is_err());
    }
}
