//! 円筒のセルの衝突確率の、弦の求積の点の数ごとの時間と差を書く（#0286 の衝突確率の近似の数値）。
//!
//! 使い方: cp_quadrature
//!
//! セルは 2A の燃料のピンの等価なセル（燃料 5 リング（等体積）、真空のギャップ、被覆管、水）。燃料の全断面積を
//! 0.3〜3000 /cm に振り、R_ij = W_ij Σ_j（真空の列は W_ij）を、小区間あたり 64 点の求積と比べる。差は ADR-0156 の測り方
//! （真空でない列は |ΔR_ij| / V_i、真空の列はその列の最大の値との比）。時間は 1 回の `Cylinder::tracks` の平均（1 スレッド）。

use std::time::Instant;

use anyhow::Context;
use xlc_lattice::cp::Cylinder;

const FUEL: usize = 5;
const GAP: usize = FUEL;

fn radii() -> Vec<f64> {
    let outer: f64 = 0.4096;
    let mut r: Vec<f64> = (1..=FUEL)
        .map(|k| outer * (k as f64 / FUEL as f64).sqrt())
        .collect();
    r.extend([0.418, 0.475, 0.711]);
    r
}

fn collisions(cell: &Cylinder, totals: &[f64]) -> anyhow::Result<Vec<f64>> {
    let w = cell.tracks(totals)?;
    let n = totals.len();
    let mut r = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            r[i * n + j] = if totals[j] > 0.0 {
                w[i][j] * totals[j]
            } else {
                w[i][j]
            };
        }
    }
    Ok(r)
}

fn difference(got: &[f64], exact: &[f64], volumes: &[f64]) -> f64 {
    let n = volumes.len();
    let mut worst: f64 = 0.0;
    for i in 0..n {
        for j in 0..n {
            let d = (got[i * n + j] - exact[i * n + j]).abs();
            let scale = if j == GAP {
                (0..n).map(|k| exact[k * n + j]).fold(0.0, f64::max)
            } else {
                volumes[i]
            };
            worst = worst.max(d / scale);
        }
    }
    worst
}

fn main() -> anyhow::Result<()> {
    let reference = Cylinder::with_points(radii(), 64).context("reference cell")?;
    let volumes = reference.volumes();
    let fuels = [0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0];
    let totals = |fuel: f64| {
        let mut t = vec![fuel; FUEL];
        t.extend([0.0, 0.3, 1.0]);
        t
    };
    let exact: Vec<Vec<f64>> = fuels
        .iter()
        .map(|&f| collisions(&reference, &totals(f)))
        .collect::<anyhow::Result<_>>()?;
    println!("points microseconds worst {fuels:?}");
    for points in [16, 12, 10, 8, 6, 4, 3, 2] {
        let cell = Cylinder::with_points(radii(), points).context("cell")?;
        let repeat = 2000;
        let clock = Instant::now();
        for k in 0..repeat {
            std::hint::black_box(cell.tracks(&totals(fuels[k % fuels.len()]))?);
        }
        let micro = clock.elapsed().as_secs_f64() / repeat as f64 * 1e6;
        let diffs: Vec<String> = fuels
            .iter()
            .zip(&exact)
            .map(|(&f, e)| {
                collisions(&cell, &totals(f))
                    .map(|r| format!("{:.1e}", difference(&r, e, &volumes)))
            })
            .collect::<anyhow::Result<_>>()?;
        println!("{points} {micro:.1} {}", diffs.join(" "));
    }
    Ok(())
}
