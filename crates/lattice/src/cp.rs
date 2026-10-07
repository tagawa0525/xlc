//! 円筒のセル（同心円のリング、外周は白色境界）の衝突確率法。

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * b.abs().max(1.0)
    }

    #[test]
    fn bickley_function_has_known_values() {
        // Ki3(0) = π/4、Ki3 は単調に減る
        assert!(close(ki3(0.0), PI / 4.0, 1e-12));
        // Ki3(1) = ∫_0^{π/2} cos²θ exp(−1/cosθ) dθ（別の求積で求めた値）
        assert!(close(ki3(1.0), 0.146_276_052_3, 1e-8), "{}", ki3(1.0));
        assert!(ki3(2.0) < ki3(1.0) && ki3(20.0) > 0.0);
    }

    /// 燃料・ギャップ・被覆管・水のセル
    const RADII: [f64; 4] = [0.4096, 0.418, 0.475, 0.711];

    #[test]
    fn probabilities_are_conserved_and_reciprocal() {
        let totals = [0.6, 0.0, 0.3, 1.2];
        let p = Cylinder::new(RADII.to_vec())
            .unwrap()
            .probabilities(&totals)
            .unwrap();
        let volumes = Cylinder::new(RADII.to_vec()).unwrap().volumes();
        for i in [0, 2, 3] {
            // 白色境界では、どこかで必ず衝突する
            let sum: f64 = (0..4).map(|j| p[i][j]).sum();
            assert!(close(sum, 1.0, 1e-8), "row {i}: {sum}");
            for j in [0, 2, 3] {
                let (a, b) = (
                    volumes[i] * totals[i] * p[i][j],
                    volumes[j] * totals[j] * p[j][i],
                );
                assert!(close(a, b, 1e-8), "{i},{j}: {a} vs {b}");
            }
            // 真空では衝突しない
            assert_eq!(p[i][1], 0.0);
        }
    }

    #[test]
    fn homogeneous_cell_is_an_infinite_medium() {
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        let flux = cell.fixed_source(&[0.8; 4], &[2.4; 4]).unwrap();
        for f in &flux {
            assert!(close(*f, 3.0, 1e-8), "{f}");
        }
    }

    #[test]
    fn void_ring_carries_the_surrounding_flux() {
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        let flux = cell
            .fixed_source(&[0.8, 0.0, 0.8, 0.8], &[2.4, 0.0, 2.4, 2.4])
            .unwrap();
        for i in [0, 2, 3] {
            assert!(close(flux[i], 3.0, 1e-8), "{i}: {}", flux[i]);
        }
    }

    #[test]
    fn rejects_bad_radii() {
        assert!(Cylinder::new(vec![]).is_err());
        assert!(Cylinder::new(vec![0.5, 0.4]).is_err());
        assert!(Cylinder::new(vec![0.0, 0.4]).is_err());
    }
}
