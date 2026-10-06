//! セルごとの核分裂率（ピン出力）。

use crate::geometry::Lattice;
use crate::material::Material;
use crate::moc::Solution;

/// セルごとの核分裂率 Σ_領域 Σ_g Σf φ V を、核分裂のあるセルの平均が 1 になるよう規格化する。
/// セルの順は `Lattice` と同じ（上の行から、各行は左から）。
pub fn cell_fission_rates(
    solution: &Solution,
    lattice: &Lattice,
    materials: &[&Material],
) -> Result<Vec<f64>, String> {
    let cells = lattice.cells_per_side().pow(2);
    let rates: Vec<f64> = (0..cells)
        .map(|cell| {
            lattice
                .cell_regions(cell)
                .map(|i| {
                    let phi = &solution.flux[i];
                    solution.volumes[i]
                        * materials[i]
                            .fission
                            .iter()
                            .zip(phi)
                            .map(|(f, p)| f * p)
                            .sum::<f64>()
                })
                .sum()
        })
        .collect();
    let fissile: Vec<f64> = rates.iter().copied().filter(|&r| r > 0.0).collect();
    if fissile.is_empty() {
        return Err("no cell has fission".into());
    }
    let mean = fissile.iter().sum::<f64>() / fissile.len() as f64;
    Ok(rates.iter().map(|r| r / mean).collect())
}

/// 核分裂率を、格子の行ごとに1行ずつの表にする。
pub fn to_text(rates: &[f64], n: usize) -> String {
    rates
        .chunks(n)
        .map(|row| {
            let values: Vec<String> = row.iter().map(|r| format!("{r:.5}")).collect();
            values.join(" ") + "\n"
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::PinCell;
    use crate::material::FissionTerm;
    use crate::moc::{Options, solve};
    use crate::tracking::Tracking;

    fn material(total: f64, scatter: f64, nu_fission: f64) -> Material {
        Material {
            total: vec![total],
            scatter: vec![vec![scatter]],
            fission: vec![nu_fission / 2.5],
            nu_fission: vec![nu_fission],
            fission_terms: if nu_fission > 0.0 {
                vec![FissionTerm {
                    chi: vec![1.0],
                    nu_fission: vec![nu_fission],
                }]
            } else {
                vec![]
            },
        }
    }

    #[test]
    fn symmetric_pins_have_equal_rates() {
        // 対角にピンを置いた 2×2 の格子は、鏡面反射で2本のピンが等価
        let pin = PinCell::new(1.26, vec![0.4]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
        let (fuel, moderator) = (material(0.5, 0.3, 0.3), material(1.0, 0.99, 0.0));
        let materials = [&fuel, &moderator, &moderator, &moderator, &fuel, &moderator];
        let tracking = Tracking::new(&lattice, 16, 0.05).unwrap();
        let solution = solve(&tracking, &materials, &Options::default()).unwrap();
        let rates = cell_fission_rates(&solution, &lattice, &materials).unwrap();
        assert_eq!(rates.len(), 4);
        assert_eq!((rates[1], rates[2]), (0.0, 0.0));
        assert!(
            (rates[0] - 1.0).abs() < 1e-6 && (rates[3] - 1.0).abs() < 1e-6,
            "{rates:?}"
        );
    }

    #[test]
    fn text_has_one_line_per_row() {
        assert_eq!(
            to_text(&[1.0, 0.0, 0.5, 1.5], 2),
            "1.00000 0.00000\n0.50000 1.50000\n"
        );
    }
}
