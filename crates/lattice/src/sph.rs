//! 等価なセルの SPH 因子。

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::Material;

    /// 群 0 が範囲の上、群 1〜3 が範囲の中
    const GROUPS: std::ops::Range<usize> = 1..4;
    const RADII: [f64; 4] = [0.3, 0.4096, 0.418, 0.7];

    fn material(total: [f64; 4], scattering: f64) -> Material {
        // 散乱は自群と次の群へ（下方散乱だけ）
        let scatter = (0..4)
            .map(|g| {
                let mut row = vec![0.0; 4];
                row[g] = 0.6 * scattering * total[g];
                if g + 1 < 4 {
                    row[g + 1] = 0.4 * scattering * total[g];
                }
                row
            })
            .collect();
        Material {
            total: total.to_vec(),
            scatter,
            fission: vec![0.0; 4],
            nu_fission: vec![0.0; 4],
            fission_terms: Vec::new(),
        }
    }

    fn materials(fuel_total: f64) -> Vec<Material> {
        let fuel = material([0.5, fuel_total, 1.0, 0.8], 0.3);
        let water = material([0.6, 1.0, 1.1, 1.2], 0.99);
        vec![fuel.clone(), fuel, Material::void(4), water]
    }

    #[test]
    fn consistent_reference_needs_no_correction() {
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        let reference = multigroup_flux(&cell, &materials(2.0), GROUPS, &[1.0]).unwrap();
        let factors = sph_factors(&cell, &materials(2.0), GROUPS, &[1.0], &reference).unwrap();
        for row in &factors {
            for mu in row {
                assert!((mu - 1.0).abs() < 1e-10, "{mu}");
            }
        }
    }

    #[test]
    fn factors_reproduce_the_reference() {
        // 燃料の全断面積を変えた解を参照にする
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        let reference = multigroup_flux(&cell, &materials(3.0), GROUPS, &[1.0]).unwrap();
        let original = materials(2.0);
        let factors = sph_factors(&cell, &original, GROUPS, &[1.0], &reference).unwrap();
        let scaled: Vec<Material> = original
            .iter()
            .zip(&factors)
            .map(|(m, mu)| {
                let mut all = vec![1.0; 4];
                all[GROUPS].copy_from_slice(mu);
                m.scaled(&all)
            })
            .collect();
        let flux = multigroup_flux(&cell, &scaled, GROUPS, &[1.0]).unwrap();
        let volumes = cell.volumes();
        for k in 0..GROUPS.len() {
            // 真空でないリングでは μ φ̃ が参照に一致する
            for r in [0, 1, 3] {
                let (a, b) = (factors[r][k] * flux[r][k], reference[r][k]);
                assert!((a - b).abs() < 1e-8 * b, "{r} {k}: {a} vs {b}");
            }
            // 真空のリングは 1、真空でないリングの中性子束と体積の積の和を保つ
            assert_eq!(factors[2][k], 1.0);
            let sum =
                |f: &Vec<Vec<f64>>| [0, 1, 3].iter().map(|&r| volumes[r] * f[r][k]).sum::<f64>();
            assert!((sum(&flux) - sum(&reference)).abs() < 1e-8 * sum(&reference));
        }
        // 燃料の全断面積を大きくした参照なので、燃料の μ は 1 と違う
        assert!((factors[0][0] - 1.0).abs() > 1e-3, "{:?}", factors[0]);
    }

    #[test]
    fn rejects_mismatched_sizes() {
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        assert!(multigroup_flux(&cell, &materials(2.0)[..3], GROUPS, &[1.0]).is_err());
        assert!(multigroup_flux(&cell, &materials(2.0), GROUPS, &[]).is_err());
        let reference = vec![vec![1.0; 2]; 4];
        assert!(sph_factors(&cell, &materials(2.0), GROUPS, &[1.0], &reference).is_err());
    }
}
