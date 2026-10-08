//! 格子の対称性（正方形の回転と鏡映）による、ピンの組。
//!
//! 正方形の 8 つの対称操作のうち、セルの配置と大きさ、領域の物質、外周の境界条件を変えないものを調べ、それで
//! 移り合うセルを同じ組にする。組の中のピンは、自己遮蔽の結果が同じなので1回だけ解けばよい（VERA の集合体の
//! ように 1/8 対称なら約 1/8 になる）。対称性は前提とせず、真空境界や位置で違う物質のある検証の体系では、成り立つ
//! 操作だけを使う。

use crate::geometry::Cell;
use crate::tracking::{Boundary, Tracking};

/// 対称性のエラー。
#[derive(Debug, thiserror::Error)]
pub enum SymmetryError {
    #[error("{given} region materials for {regions} regions")]
    Regions { given: usize, regions: usize },
}

/// 正方形の 8 つの対称操作。格子の中心を原点とするセルの座標 (u, v)（u は右、v は下向き）に掛ける行列
/// [a, b, c, d]: (u, v) → (a u + b v, c u + d v)。恒等、回転（3つ）、左右・上下・対角線・反対角線の鏡映
const OPERATIONS: [[i64; 4]; 8] = [
    [1, 0, 0, 1],
    [0, -1, 1, 0],
    [-1, 0, 0, -1],
    [0, 1, -1, 0],
    [-1, 0, 0, 1],
    [1, 0, 0, -1],
    [0, 1, 1, 0],
    [0, -1, -1, 0],
];

/// セルごとの、同じ組の代表のセル（組の中で最も小さい番号）。`region_material` は領域ごとの物質の番号。
pub fn cell_orbits(
    tracking: &Tracking,
    region_material: &[usize],
) -> Result<Vec<usize>, SymmetryError> {
    let lattice = &tracking.lattice;
    let regions = lattice.regions();
    if region_material.len() != regions {
        return Err(SymmetryError::Regions {
            given: region_material.len(),
            regions,
        });
    }
    let (nc, nr) = (lattice.columns(), lattice.rows());
    let cells = nc * nr;
    // セルの座標は2倍して整数にする: U = 2c − (nc − 1)、V = 2r − (nr − 1)
    let (oc, or) = (nc as i64 - 1, nr as i64 - 1);
    let transform =
        |op: &[i64; 4], (u, v): (i64, i64)| (op[0] * u + op[1] * v, op[2] * u + op[3] * v);
    // セルの像（格子の外に出れば None）
    let map = |op: &[i64; 4], cell: usize| -> Option<usize> {
        let (u, v) = (2 * (cell % nc) as i64 - oc, 2 * (cell / nc) as i64 - or);
        let (u, v) = transform(op, (u, v));
        let (c, r) = (u + oc, v + or);
        let inside = c >= 0 && r >= 0 && c % 2 == 0 && r % 2 == 0;
        let (c, r) = (c / 2, r / 2);
        (inside && (c as usize) < nc && (r as usize) < nr).then(|| r as usize * nc + c as usize)
    };
    // セルのリングごとの物質（扇形や小さな正方形で物質が違うセルは None: 回すと領域の対応が変わる）
    let rings = |cell: usize| -> Option<Vec<usize>> {
        let materials = &region_material[lattice.cell_regions(cell)];
        let width = match lattice.cell_type(cell) {
            Cell::Pin(pin) if !pin.is_grid() => pin.sectors(),
            _ => materials.len(),
        };
        materials
            .chunks(width)
            .map(|ring| ring.iter().all(|&m| m == ring[0]).then_some(ring[0]))
            .collect()
    };
    let ring_materials: Vec<Option<Vec<usize>>> = (0..cells).map(rings).collect();
    // 外周の辺の向き (u, v) と境界条件
    let b = tracking.boundaries();
    let sides = [
        ((-1, 0), b.left),
        ((1, 0), b.right),
        ((0, -1), b.top),
        ((0, 1), b.bottom),
    ];
    let side = |direction: (i64, i64)| -> Option<Boundary> {
        sides.iter().find(|(d, _)| *d == direction).map(|(_, b)| *b)
    };
    let valid: Vec<&[i64; 4]> = OPERATIONS
        .iter()
        .filter(|op| {
            let identity = **op == OPERATIONS[0];
            let boundaries = sides
                .iter()
                .all(|&(d, condition)| side(transform(op, d)) == Some(condition));
            // 軸を入れ替える操作（x 軸の像が y 軸の上）では、セルの幅と高さも入れ替わる
            let swap = transform(op, (1, 0)).0 == 0;
            let cells_match = (0..cells).all(|cell| {
                let Some(image) = map(op, cell) else {
                    return false;
                };
                let (w, h) = lattice.cell_size(cell);
                let size = if swap { (h, w) } else { (w, h) };
                size == lattice.cell_size(image)
                    && lattice.cell_type(cell) == lattice.cell_type(image)
                    && ring_materials[cell].is_some()
                    && ring_materials[cell] == ring_materials[image]
            });
            identity || (boundaries && cells_match)
        })
        .collect();
    Ok((0..cells)
        .map(|cell| {
            valid
                .iter()
                .filter_map(|op| map(op, cell))
                .min()
                .unwrap_or(cell)
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};
    use crate::tracking::{Boundaries, Boundary, Tracking};

    /// 3×3 の格子。`water` のセルは水、ほかはピン（燃料と水の 2 領域）。物質の番号: 燃料 0、水 1
    fn lattice(water: &[usize]) -> (Lattice, Vec<usize>) {
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let box_ = PinCell::new(1.26, vec![]).unwrap();
        let layout: Vec<Vec<usize>> = (0..3)
            .map(|r| {
                (0..3)
                    .map(|c| usize::from(water.contains(&(3 * r + c))))
                    .collect()
            })
            .collect();
        let materials = (0..9)
            .flat_map(|cell| {
                if water.contains(&cell) {
                    vec![1]
                } else {
                    vec![0, 1]
                }
            })
            .collect();
        (Lattice::new(vec![pin, box_], layout).unwrap(), materials)
    }

    fn orbits(boundaries: Boundaries, materials: Option<Vec<usize>>) -> Vec<usize> {
        let (lattice, default) = lattice(&[4]);
        let tracking = Tracking::with_boundaries(&lattice, 8, 0.2, boundaries).unwrap();
        cell_orbits(&tracking, &materials.unwrap_or(default)).unwrap()
    }

    #[test]
    fn reflective_lattice_has_all_symmetries() {
        // 角の 4 本と辺の 4 本がそれぞれ同じ組、中央の水は単独
        assert_eq!(
            orbits(Boundaries::REFLECTIVE, None),
            vec![0, 1, 0, 1, 4, 1, 0, 1, 0]
        );
    }

    #[test]
    fn vacuum_sides_keep_the_diagonal_mirror() {
        // 右と下が真空: 左上から右下への対角線の鏡映だけが残る
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        assert_eq!(orbits(boundaries, None), vec![0, 1, 2, 1, 4, 5, 2, 5, 8]);
    }

    #[test]
    fn materials_break_the_symmetry() {
        // 右下の角のピンだけ燃料の物質が違う（番号 2）: その角を動かさない対角線の鏡映だけが残る
        let (_, mut materials) = lattice(&[4]);
        let last = materials.len() - 2;
        materials[last] = 2;
        assert_eq!(
            orbits(Boundaries::REFLECTIVE, Some(materials)),
            vec![0, 1, 2, 1, 4, 5, 2, 5, 8]
        );
    }

    #[test]
    fn region_materials_must_cover_the_lattice() {
        let (lattice, materials) = lattice(&[4]);
        let tracking = Tracking::new(&lattice, 8, 0.2).unwrap();
        assert!(cell_orbits(&tracking, &materials[1..]).is_err());
    }
}
