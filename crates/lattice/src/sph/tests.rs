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
        anisotropic: Vec::new(),
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
    let factors = sph_factors(
        &cell,
        &materials(2.0),
        &[false; 4],
        GROUPS,
        &[1.0],
        &reference,
    )
    .unwrap();
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
    let factors = sph_factors(&cell, &original, &[false; 4], GROUPS, &[1.0], &reference).unwrap();
    let scaled: Vec<Material> = original
        .iter()
        .zip(&factors)
        .map(|(m, mu)| {
            let mut all = vec![1.0; 4];
            all[GROUPS].copy_from_slice(mu);
            m.scaled(&all).unwrap()
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
fn fixed_rings_keep_one() {
    // 掛ける先のないリング（等価なセルの枠）は μ = 1 のまま、残りのリングで形を合わせて規格化する
    let cell = Cylinder::new(RADII.to_vec()).unwrap();
    let reference = multigroup_flux(&cell, &materials(3.0), GROUPS, &[1.0]).unwrap();
    let original = materials(2.0);
    let fixed = [false, false, false, true];
    let factors = sph_factors(&cell, &original, &fixed, GROUPS, &[1.0], &reference).unwrap();
    let scaled: Vec<Material> = original
        .iter()
        .zip(&factors)
        .map(|(m, mu)| {
            let mut all = vec![1.0; 4];
            all[GROUPS].copy_from_slice(mu);
            m.scaled(&all).unwrap()
        })
        .collect();
    let flux = multigroup_flux(&cell, &scaled, GROUPS, &[1.0]).unwrap();
    let volumes = cell.volumes();
    for k in 0..GROUPS.len() {
        assert_eq!(factors[3][k], 1.0);
        let ratio = |r: usize| factors[r][k] * flux[r][k] / reference[r][k];
        assert!((ratio(1) - ratio(0)).abs() < 1e-8 * ratio(0));
        let preserved: f64 = [0, 1].iter().map(|&r| volumes[r] * reference[r][k]).sum();
        let normalized: f64 = [0, 1]
            .iter()
            .map(|&r| volumes[r] * reference[r][k] / factors[r][k])
            .sum();
        assert!((normalized - preserved).abs() < 1e-10 * preserved);
    }
    // 固定するリングの数が合わなければエラー
    assert!(sph_factors(&cell, &original, &[false; 3], GROUPS, &[1.0], &reference).is_err());
}

#[test]
fn rejects_rings_without_flux() {
    // 範囲の上から散乱してくる中性子がなければ、多群の中性子束が 0 になり、μ が決まらない
    let cell = Cylinder::new(RADII.to_vec()).unwrap();
    let reference = multigroup_flux(&cell, &materials(2.0), GROUPS, &[1.0]).unwrap();
    let mut isolated = materials(2.0);
    for m in &mut isolated {
        m.scatter[0] = vec![0.0; 4];
    }
    assert!(matches!(
        sph_factors(&cell, &isolated, &[false; 4], GROUPS, &[1.0], &reference),
        Err(SphError::Flux { .. })
    ));
}

#[test]
fn rejects_references_without_flux() {
    // 真空でないリングの参照の中性子束が 0 や NaN なら、μ が決まらない
    let cell = Cylinder::new(RADII.to_vec()).unwrap();
    for bad in [0.0, f64::NAN, -1.0] {
        let mut reference = multigroup_flux(&cell, &materials(2.0), GROUPS, &[1.0]).unwrap();
        reference[0][1] = bad;
        assert!(
            matches!(
                sph_factors(
                    &cell,
                    &materials(2.0),
                    &[false; 4],
                    GROUPS,
                    &[1.0],
                    &reference
                ),
                Err(SphError::ReferenceValue { ring: 0, group: 2 })
            ),
            "{bad}"
        );
    }
}

#[test]
fn rejects_mismatched_sizes() {
    let cell = Cylinder::new(RADII.to_vec()).unwrap();
    assert!(multigroup_flux(&cell, &materials(2.0)[..3], GROUPS, &[1.0]).is_err());
    assert!(multigroup_flux(&cell, &materials(2.0), GROUPS, &[]).is_err());
    let reference = vec![vec![1.0; 2]; 4];
    assert!(
        sph_factors(
            &cell,
            &materials(2.0),
            &[false; 4],
            GROUPS,
            &[1.0],
            &reference
        )
        .is_err()
    );
    // 物質の群より先まで求める群の範囲
    let reference = vec![vec![1.0; 5]; 4];
    assert!(
        sph_factors(
            &cell,
            &materials(2.0),
            &[false; 4],
            1..6,
            &[1.0],
            &reference
        )
        .is_err()
    );
}
