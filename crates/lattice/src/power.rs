//! セルごとの核分裂率（ピン出力）。

use crate::geometry::Lattice;
use crate::material::Material;
use crate::moc::Solution;
use xlc_format::msgpack::MsgpackError;
use xlc_format::pins::PinPowers;

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

/// 核分裂率（`cell_fission_rates`）の、ピンセルのある行と列の値（`xlc_format::pins`）。余白や反射体のように
/// ピンセルのない行と列は書かない。
pub fn to_msgpack(rates: &[f64], lattice: &Lattice) -> Result<Vec<u8>, MsgpackError> {
    PinPowers {
        rows: pin_rows(rates, lattice),
    }
    .to_msgpack()
}

/// ピンセルのある行と列だけを取り出した、行ごとの値
fn pin_rows(rates: &[f64], lattice: &Lattice) -> Vec<Vec<f64>> {
    let columns = lattice.columns();
    let (rows, kept) = lattice.pin_rows_and_columns();
    rows.iter()
        .map(|&r| kept.iter().map(|&c| rates[r * columns + c]).collect())
        .collect()
}

#[cfg(test)]
mod tests;
