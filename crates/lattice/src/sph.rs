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
    fn factors_reproduce_the_reference_shape() {
        // 燃料の全断面積を変えた解を参照にする。群の釣り合い（反応率の和 = 上の群からの線源）は μ では変えられない
        // ので、群ごとのリング間の形 μ φ̃ ∝ φ_ref を再現し、大きさは Σ V φ_ref / μ = Σ V φ_ref で決める
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
            // 真空でないリングで μ φ̃ / φ_ref が等しい
            let ratio = |r: usize| factors[r][k] * flux[r][k] / reference[r][k];
            for r in [1, 3] {
                assert!(
                    (ratio(r) - ratio(0)).abs() < 1e-8 * ratio(0),
                    "{r} {k}: {} vs {}",
                    ratio(r),
                    ratio(0)
                );
            }
            let preserved: f64 = [0, 1, 3]
                .iter()
                .map(|&r| volumes[r] * reference[r][k])
                .sum();
            let normalized: f64 = [0, 1, 3]
                .iter()
                .map(|&r| volumes[r] * reference[r][k] / factors[r][k])
                .sum();
            assert!((normalized - preserved).abs() < 1e-10 * preserved);
            // 真空のリングには掛ける断面積がないので 1
            assert_eq!(factors[2][k], 1.0);
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
