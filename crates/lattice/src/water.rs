//! ホウ酸水の組成（ADR-0106）。

use std::collections::BTreeMap;

/// AME2020 の原子質量 (u)（M. Wang et al., Chinese Phys. C 45, 030003 (2021)）
const MASS_H1: f64 = 1.007_825_032_23;
const MASS_O16: f64 = 15.994_914_619_57;
const MASS_B10: f64 = 10.012_936_95;
const MASS_B11: f64 = 11.009_305_36;
/// アボガドロ定数 × 10⁻²⁴（個数密度を 1/(barn cm) にする）
const AVOGADRO: f64 = 0.602_214_076;

/// 密度 (g/cm³)、ホウ素の質量 ppm、ホウ素の中の B-10 の原子分率から、核種の名前の順の個数密度 (1/(barn cm)) を
/// 作る。水は H-1 2 個と O-16 1 個の分子とし、原子分率が 0 の同位体は入れない。
pub fn borated_water(density: f64, boron: f64, b10: f64) -> Result<BTreeMap<String, f64>, String> {
    if !(density > 0.0 && density.is_finite()) {
        return Err(format!("invalid water density {density}"));
    }
    if !(0.0..1.0e6).contains(&boron) {
        return Err(format!("invalid boron {boron} ppm"));
    }
    if !(0.0..=1.0).contains(&b10) {
        return Err(format!("invalid B-10 atom fraction {b10}"));
    }
    let w = boron * 1.0e-6;
    let water = density * (1.0 - w) * AVOGADRO / (2.0 * MASS_H1 + MASS_O16);
    let mut composition = BTreeMap::from([
        ("H-1".to_string(), 2.0 * water),
        ("O-16".to_string(), water),
    ]);
    if boron > 0.0 {
        let mass = b10 * MASS_B10 + (1.0 - b10) * MASS_B11;
        let atoms = density * w * AVOGADRO / mass;
        // 原子分率が 0 の同位体は入れない（ライブラリに要らない核種を求めない）
        for (nuclide, fraction) in [("B-10", b10), ("B-11", 1.0 - b10)] {
            if fraction > 0.0 {
                composition.insert(nuclide.to_string(), fraction * atoms);
            }
        }
    }
    Ok(composition)
}
