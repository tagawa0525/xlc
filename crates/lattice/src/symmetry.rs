//! 格子の対称性（正方形の回転と鏡映）による、ピンの組。

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
