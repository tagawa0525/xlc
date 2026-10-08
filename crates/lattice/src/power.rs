//! セルごとの核分裂率（ピン出力）。

use crate::geometry::Lattice;
use crate::material::Material;
use crate::moc::Solution;
use xlc_format::msgpack::MsgpackError;

/// 核分裂のあるセルがない。
#[derive(Debug, thiserror::Error)]
#[error("no cell has fission")]
pub struct NoFission;

/// セルごとの核分裂率 Σ_領域 Σ_g Σf φ V を、核分裂のあるセルの平均が 1 になるよう規格化する。
/// セルの順は `Lattice` と同じ（上の行から、各行は左から）。
pub fn cell_fission_rates(
    solution: &Solution,
    lattice: &Lattice,
    materials: &[&Material],
) -> Result<Vec<f64>, NoFission> {
    let cells = lattice.cells();
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
        return Err(NoFission);
    }
    let mean = fissile.iter().sum::<f64>() / fissile.len() as f64;
    Ok(rates.iter().map(|r| r / mean).collect())
}

/// 核分裂率（`cell_fission_rates`）を、フィールド名付きの map 形式の MessagePack にする: `format`
/// ("xlc-pin-powers") と、ピンセルのある行と列の核分裂率 `rows`（上の行から、各行は左から）。余白や反射体の
/// ようにピンセルのない行と列は書かない。
pub fn to_msgpack(rates: &[f64], lattice: &Lattice) -> Result<Vec<u8>, MsgpackError> {
    #[derive(serde::Serialize)]
    struct PinPowers {
        format: &'static str,
        rows: Vec<Vec<f64>>,
    }
    xlc_format::msgpack::to_named(&PinPowers {
        format: "xlc-pin-powers",
        rows: pin_rows(rates, lattice),
    })
}

/// ピンセルのある行と列だけを取り出した、行ごとの値
fn pin_rows(rates: &[f64], lattice: &Lattice) -> Vec<Vec<f64>> {
    let columns = lattice.columns();
    let pin = |cell: usize| lattice.cell_type(cell).pin().is_some();
    let rows: Vec<usize> = (0..lattice.rows())
        .filter(|&r| (0..columns).any(|c| pin(r * columns + c)))
        .collect();
    let kept: Vec<usize> = (0..columns)
        .filter(|&c| (0..lattice.rows()).any(|r| pin(r * columns + c)))
        .collect();
    rows.iter()
        .map(|&r| kept.iter().map(|&c| rates[r * columns + c]).collect())
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
            anisotropic: Vec::new(),
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
    fn pin_powers_are_written_by_rows() {
        #[derive(serde::Deserialize, Debug, PartialEq)]
        #[serde(deny_unknown_fields)]
        struct PinPowers {
            format: String,
            rows: Vec<Vec<f64>>,
        }
        let pin = PinCell::new(1.26, vec![0.4]).unwrap();
        let lattice = Lattice::new(vec![pin], vec![vec![0, 0], vec![0, 0]]).unwrap();
        let bytes = to_msgpack(&[1.0, 0.0, 0.5, 1.5], &lattice).unwrap();
        let value: PinPowers = rmp_serde::from_slice(&bytes).unwrap();
        assert_eq!(
            value,
            PinPowers {
                format: "xlc-pin-powers".into(),
                rows: vec![vec![1.0, 0.0], vec![0.5, 1.5]],
            }
        );
        // フィールド名付きの map 形式
        assert!(xlc_format::msgpack::read_value(&bytes).is_ok());
        // 余白の行と列は書かない
        let pin = PinCell::new(1.26, vec![0.4]).unwrap();
        let margin = Lattice::with_margin(vec![pin], vec![vec![0]], 0.1).unwrap();
        let mut rates = vec![0.0; 9];
        rates[4] = 1.0;
        let value: PinPowers =
            rmp_serde::from_slice(&to_msgpack(&rates, &margin).unwrap()).unwrap();
        assert_eq!(value.rows, vec![vec![1.0]]);
    }
}
