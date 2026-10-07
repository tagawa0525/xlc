//! 散乱のない固定源の MOC（ESSM の背景断面積を求めるのに使う）。

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};

    fn pin(radii: Vec<f64>) -> Tracking {
        let cell = PinCell::new(1.26, radii).unwrap();
        Tracking::new(&Lattice::pin_cell(cell), 16, 0.02).unwrap()
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-9 * b.abs().max(1.0)
    }

    #[test]
    fn uniform_source_ratio_gives_that_flux() {
        let tracking = pin(vec![0.4]);
        let flux = fixed_source(
            &tracking,
            &[vec![2.0, 0.3], vec![0.5, 1.2]],
            &[vec![6.0, 0.9], vec![1.5, 3.6]],
            &FixedSourceOptions::default(),
        )
        .unwrap();
        for f in flux.iter().flatten() {
            assert!(close(*f, 3.0), "{f}");
        }
    }

    #[test]
    fn void_regions_carry_the_surrounding_flux() {
        let tracking = pin(vec![0.4, 0.45]);
        let flux = fixed_source(
            &tracking,
            &[vec![1.0], vec![0.0], vec![1.0]],
            &[vec![2.0], vec![0.0], vec![2.0]],
            &FixedSourceOptions::default(),
        )
        .unwrap();
        for f in flux.iter().flatten() {
            assert!(close(*f, 2.0), "{f}");
        }
    }

    #[test]
    fn black_regions_absorb_everything() {
        let tracking = pin(vec![0.4]);
        let flux = fixed_source(
            &tracking,
            &[vec![f64::INFINITY], vec![1.0]],
            &[vec![0.0], vec![1.0]],
            &FixedSourceOptions::default(),
        )
        .unwrap();
        assert_eq!(flux[0][0], 0.0);
        assert!(flux[1][0] > 0.0 && flux[1][0] < 1.0, "{}", flux[1][0]);
    }
}
