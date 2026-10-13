//! MOC の解から、少数群定数（`xlc_format::fewgroup::FewGroup`）を作る。

use std::ops::Range;

use xlc_format::fewgroup::{FewGroup, Pins};
use xlc_format::group::GroupStructure;
use xlc_numeric::equivalence::{Choice, reflector_diffusion};
use xlc_numeric::nem::{Edge, OneNode, OneNodeError, is_fissile, surface_fluxes};
use xlc_numeric::pins::corners::{CornerNode, corner_balance, fuel_harmonic_mean, nodes_at_point};
use xlc_numeric::pins::node::{Conditions, Constants, NodeModel};
use xlc_numeric::pins::{CELL_POINTS, FACE_POINTS};
use xlc_numeric::quadrature::gauss_legendre;

use crate::b1::{B1Error, Medium};
use crate::cmfd::{horizontal_face, vertical_face};
use crate::geometry::Lattice;
use crate::material::Material;
use crate::medium::fine_medium;
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
    #[error("medium: {0}")]
    Medium(String),
}

/// 格子のセルの長方形の範囲で、1つの少数群定数にまとめるもの（ADR-0095）。
#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    /// セルの列（左から、0 始まり）の範囲
    pub columns: Range<usize>,
    /// セルの行（上から、0 始まり）の範囲
    pub rows: Range<usize>,
    /// 反射体のノード: 反射体でないノード（燃料）と接する辺の ADF を、その燃料のノードの同じ面の ADF で割った
    /// f* = f_反射体 / f_燃料 で持ち、燃料のノードに接する角の点の CDF を燃料の CDF で割って持つ（ADR-0131）。
    /// 燃料と辺で接しなくてよい
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
/// 鏡面反射なら正味の流れは 0 で、均質な解は平らになり、体積平均の中性子束で割るのと同じになる。
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
    Ok(with_z_sides(few))
}

/// 格子をノード `nodes`（格子のセルを重ならずに覆う長方形）に分け、ノードごとに `homogenize` と同じ形で
/// 少数群定数を作る（ADR-0095）。返す定数はノードの順。
///
/// - 断面積は、ノードの領域だけで反応率を保存する。拡散係数は、ノードの4辺がすべて外周の鏡面反射の辺
///   （格子全体が1つのノード）か、領域がすべて1つの核分裂のある媒質（`media` は領域ごとの媒質の番号、媒質でなければ `None`。ADR-0141）なら B1、ほかは、核分裂の
///   あるノードではノードの領域の中性子流で重みづけした outflow、
///   核分裂のないノードでは Koebke の等価理論の D（正味の中性子流の向きの2つの面の不連続因子をそろえる。
///   ADR-0104）。`k_inf` は格子の k
/// - 辺の不均質な中性子束と正味の中性子流は、辺の上の粗メッシュの面の値（MOC の `face_fluxes`、
///   `face_currents`）を、面の長さで重みづけして平均する。内側の辺も同じ
/// - ADF は1ノードの問題（ADR-0094）で作る。外周の辺は格子の境界（鏡面反射か真空）、内側の辺は隣の
///   ノード（その幅と、参照の横方向の漏れ。core と同じ2次式）とする。内側の辺は、向かいにちょうど1つの
///   ノードがあり、辺の長さがそろうこと
/// - 外周の真空の辺には J/φ_het を持たせる（ADR-0101）
/// - 反射体のノード（`Node::reflector`）の、反射体でないノードと接する辺の ADF は、その隣のノードの同じ面の
///   ADF で割った f* とし、`normalized_adf` に印を付ける。反射体どうしの辺と、燃料に角の点だけで接する反射体の
///   ノード（角型の格子の右上。ADR-0095）や燃料と接しない反射体のノード（軸方向のノズル、炉心板、水の層。
///   ADR-0128）の辺は f* でない GET の ADF を持つ
/// - `pins` なら、形状関数と角の不連続因子（`node_pins`。ADR-0131）も付ける
#[allow(clippy::too_many_arguments)]
pub fn homogenize_nodes(
    solution: &Solution,
    lattice: &Lattice,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
    boundaries: &Boundaries,
    nodes: &[Node],
    media: &[Option<usize>],
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
        media,
    )?;
    let sides = sides_of(boundaries);
    let coarse = coarse_bounds.len() - 1;
    // 横方向の漏れ [ノード][方向（x, y）][少数群]: x の方向の横方向の漏れは y の正味の漏れ（上と下の辺の
    // 外向きの正味の流れの和 / 高さ）、y の方向は x の正味の漏れ
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
        // 核分裂のないノードの拡散係数は、Koebke の等価理論の D（ADR-0104）。中性子流の重みの D から探し、正の面の範囲に
        // なければ中性子流の重みの D（ADR-0134）
        let mut few = h.few.clone();
        let fissile = is_fissile(&few.nu_fission)
            .map_err(|e| HomogenizeError::OneNode(OneNodeError::Invalid(e.to_string())))?;
        if !fissile {
            let (diffusion, choice) = reflector_diffusion(
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
            few.diffusion = diffusion;
            match choice {
                Choice::Koebke => tracing::info!(
                    node = n,
                    current_weighted = ?h.few.diffusion,
                    equivalent = ?few.diffusion,
                    "equivalent diffusion coefficients"
                ),
                Choice::CurrentWeighted {
                    group,
                    target,
                    range,
                } => tracing::info!(
                    node = n,
                    current_weighted = ?few.diffusion,
                    group,
                    heterogeneous_ratio = target,
                    homogeneous_ratio_range = ?range,
                    "no Koebke diffusion coefficient with positive faces; current-weighted (ADR-0134)"
                ),
            }
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
    // 反射体のノードの検査を、形状関数の角の点の釣り合いより先にする（誤りの入力を、釣り合いの別のエラーで隠さない）
    for (n, node) in nodes.iter().enumerate() {
        if !node.reflector {
            continue;
        }
        // 燃料の CDF で割る角は、核分裂のあるノードを燃料とする（ADR-0131）。核分裂のある反射体のノードは自分の CDF で
        // 割ることになり、core で掛け戻す値と食い違う
        let fissile = is_fissile(&adfs[n].nu_fission)
            .map_err(|e| HomogenizeError::Nodes(format!("node {n}: {e}")))?;
        if pins && fissile {
            return Err(HomogenizeError::Nodes(format!(
                "reflector node {n} has fission; form functions need reflector nodes without fission"
            )));
        }
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
        }
        result.push(with_z_sides(few));
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
#[allow(clippy::too_many_arguments)]
pub(crate) fn node_references(
    solution: &Solution,
    lattice: &Lattice,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
    boundaries: &Boundaries,
    nodes: &[Node],
    media: &[Option<usize>],
) -> Result<NodeReferences, HomogenizeError> {
    let sides = sides_of(boundaries);
    if sides.iter().any(|b| matches!(b, Boundary::Albedo(_))) {
        return Err(HomogenizeError::Albedo);
    }
    check_regions(solution, materials)?;
    if media.len() != materials.len() {
        return Err(HomogenizeError::Nodes(format!(
            "{} medium entries for {} regions",
            media.len(),
            materials.len()
        )));
    }
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
        // 領域がすべて1つの核分裂のある媒質のノードは、その媒質の B1 の D（単一集合体の格子の D と同じ。ADR-0141）。
        // 中性子流の重みは、反射体の隣でノードの中の向きが変わると打ち消して定まらない。核分裂のない媒質のノードは、
        // 媒質でないノードと同じく後で Koebke の D を中性子流の重みの D から探す（ADR-0104、ADR-0134）
        let one_medium = media[regions[0]].is_some()
            && regions.iter().all(|&r| media[r] == media[regions[0]])
            && is_fissile(&materials[regions[0]].nu_fission)
                .map_err(|e| HomogenizeError::Nodes(format!("node {n}: {e}")))?;
        let (mut few, average) = constants(
            solution,
            materials,
            &regions,
            &group_of,
            coarse_bounds,
            (whole && reflective(boundaries)) || one_medium,
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
/// の格子でだけ作る（均質な解が平らなので、角の均質な中性子束はノード平均。ADR-0094。ほかの格子は `nodes` の
/// 形状関数。ADR-0131）。
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
            "form functions of the whole lattice need every side reflective (other lattices: `nodes` with `pins`, ADR-0131)"
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
/// （形状関数と角の不連続因子の角の点の釣り合い。ADR-0131）
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

/// 角の点の均質な値の大きさ s_p がどれも正であること（ADR-0131）。正でなければ CDF を作れないので、別の値に
/// 替えずにエラーにする
fn positive_scales(scales: &[[f64; 2]]) -> Result<(), String> {
    match scales
        .iter()
        .position(|s| s.iter().any(|v| v.is_nan() || *v <= 0.0))
    {
        Some(p) => Err(format!(
            "corner point {p}: the homogeneous corner flux scale {:?} is not positive",
            scales[p]
        )),
        None => Ok(()),
    }
}

/// 格子のノードごとの形状関数と角の不連続因子（ADR-0131）。`fews` は ADF（f* で割る前）と拡散係数を決めた
/// ノードの定数。ノードは格子の上の nx × ny の長方形の並びであること。
///
/// - 角の点 p の不均質な中性子束 Φ_p は `point_flux`
/// - ノード n の角 c の均質な値の見込みは h̃ = φ_x φ_y / φ̄（φ_x、φ_y は角 c を含む x と y の辺の、1ノードの問題の
///   均質な面の中性子束 = 辺の不均質な中性子束 / ADF、φ̄ はノード平均）
/// - 均質な角の値を s_p h̃ とし、core と同じ展開（参照のノード平均と辺の正味の流れ、格子の k）と角の点の
///   釣り合い（`corner_balance`。真空の辺の流れは r × ADF × 均質な角の値。ADR-0131）で s_p を解く。s_p が正で
///   なければエラー
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
    let corner_nodes: Vec<Option<CornerNode>> = grid
        .iter()
        .map(|&n| {
            Some(CornerNode {
                model: &models[n],
                conditions: &conditions[n],
                diffusion: [fews[n].diffusion[0], fews[n].diffusion[1]],
                weights: estimates[n],
                // 外周の真空の辺（面 -x、+x、-y、+y は辺 0、1、3、2）の流れ r f × 均質な角の値（f は f* で割る前の
                // ADF。core と同じ。ADR-0131）。鏡面反射は 0
                outflow: std::array::from_fn(|face| {
                    let side = FACES[face];
                    match neighbors[n][side] {
                        None => fews[n].vacuum_ratio[side]
                            .as_ref()
                            .map(|r| std::array::from_fn(|g| r[g] * fews[n].adf[side][g])),
                        Some(_) => None,
                    }
                }),
            })
        })
        .collect();
    let scales = corner_balance(nx, ny, &corner_nodes, &vec![false; points]).map_err(error)?;
    for (p, (s, f)) in scales.iter().zip(&phi).enumerate() {
        tracing::debug!(
            point = p,
            column = p % (nx + 1),
            row_from_bottom = p / (nx + 1),
            scale = ?s,
            heterogeneous = ?f,
            "corner point"
        );
    }
    positive_scales(&scales).map_err(error)?;
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
            // 辺（左、右、上、下、−z、+z）: 集合体の外周なら区間の値、内側と z の辺（`with_z_sides`）なら 1
            let outer = [i == 0, i == n - 1, j == n - 1, j == 0, false, false];
            sub.adf = (0..outer.len())
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

/// 2次元の格子の定数の辺（左、右、上、下）に z の辺（−z、+z）を足す。軸方向に一様な格子では、z の面の不均質な
/// 中性子束も均質な解も z に平らでノード平均に等しいので、GET の ADF は厳密に 1 になる（ADR-0128）。z の辺は真空の比と
/// f* の印を持たない
fn with_z_sides(mut few: FewGroup) -> FewGroup {
    let coarse = few.groups();
    few.adf.extend([vec![1.0; coarse], vec![1.0; coarse]]);
    few.vacuum_ratio.extend([None, None]);
    few.normalized_adf.extend([false, false]);
    few
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
/// `normalized_adf` は径方向の4辺がすべて false、z の厚さはなし）と、ノード平均の少数群の中性子束。`b1` なら拡散係数を B1 で、
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
        height: None,
        pins: None,
    };
    Ok((few, average))
}

/// [辺] ごとの細かい群の値を少数群に縮約する（細かい群の和。辺の平均の中性子束と正味の流れは群について
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
    // 媒質の均質化（ADR-0128）と同じ重みの P0 と P1。P1 のモーメントがあることは `diffusion` で確かめてある
    let (fine, flux) = fine_medium(solution, materials, regions)?;
    let ng = group_of.len();
    let medium = Medium {
        total: fine.total,
        scatter1: fine
            .scatter
            .get(1)
            .cloned()
            .unwrap_or_else(|| vec![vec![0.0; ng]; ng]),
        scatter: fine.scatter.into_iter().next().unwrap_or_default(),
        nu_fission: fine.nu_fission,
        chi: fine.chi,
    };
    let produced: f64 = medium.chi.iter().sum();
    let (current, weight) = if produced > 0.0 {
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
            // 打ち消しの程度: 1 なら領域と細かい群の中性子流がすべてノードの正味の流れの向き
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
mod tests;
