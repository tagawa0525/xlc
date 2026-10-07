//! 燃料を黒体にした固定源の MOC による、ピンごと・群ごとの Dancoff 係数。

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};

    fn options() -> DancoffOptions {
        DancoffOptions::default()
    }

    /// ピンと水のセルを `layout`（0: ピン、1: 水）に並べ、燃料以外の群の全断面積を `total` とした Dancoff 係数
    fn factors(pitch: f64, layout: Vec<Vec<usize>>, total: &[f64]) -> Vec<Option<Vec<f64>>> {
        let pin = PinCell::new(pitch, vec![0.4096]).unwrap();
        let water = PinCell::new(pitch, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, water], layout).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let regions = lattice.volumes().len();
        // ピンのセルは (燃料, 水) の2領域、水のセルは1領域
        let mut fuel = Vec::with_capacity(regions);
        for cell in 0..lattice.cells_per_side().pow(2) {
            let range = lattice.cell_regions(cell);
            fuel.push(range.len() == 2);
            fuel.extend(std::iter::repeat_n(false, range.len() - 1));
        }
        let totals = vec![total.to_vec(); regions];
        dancoff(&tracking, &totals, &fuel, &options()).unwrap()
    }

    #[test]
    fn distant_pins_do_not_shadow_each_other() {
        // 隣のピンまで減速材の光学的厚さが 20 程度あるので、影はほぼない
        let c = factors(20.0, vec![vec![0]], &[1.0, 2.0]);
        for value in c[0].as_ref().unwrap() {
            assert!(value.abs() < 1e-6, "{value}");
        }
    }

    #[test]
    fn transparent_moderator_shadows_completely() {
        let c = factors(1.26, vec![vec![0]], &[1.0e-6]);
        let value = c[0].as_ref().unwrap()[0];
        assert!(value > 0.999, "{value}");
    }

    #[test]
    fn repeated_pins_match_the_single_cell_lattice() {
        let total = [0.7, 1.5];
        let single = factors(1.26, vec![vec![0]], &total);
        let four = factors(1.26, vec![vec![0, 0], vec![0, 0]], &total);
        let reference = single[0].as_ref().unwrap();
        assert!(reference[0] > reference[1], "{reference:?}");
        assert!(
            reference.iter().all(|&c| c > 0.0 && c < 1.0),
            "{reference:?}"
        );
        for cell in &four {
            for (c, r) in cell.as_ref().unwrap().iter().zip(reference) {
                assert!((c - r).abs() < 2e-3, "{c} vs {r}");
            }
        }
    }

    #[test]
    fn cells_without_fuel_have_no_factors() {
        let c = factors(1.26, vec![vec![0, 1], vec![1, 0]], &[1.0]);
        assert!(c[0].is_some() && c[3].is_some());
        assert!(c[1].is_none() && c[2].is_none());
    }
}
