//! 燃料のリングの全断面積がリングごとに違うセルで、燃料の平均（中心）のまわりの展開の誤差を書く（ADR-0182 の
//! 検討した案の数値）。
//!
//! 使い方: ring_expansion
//!
//! セルは 2O の Gd のピン（燃料 8 リング、真空のギャップ、被覆管、水 6 リング）。燃料の i 番目のリングの全断面積を
//! Σ̂ e^{δ_i}（δ は体積の重みの平均が 0）とし、R_ij = W_ij Σ_j（真空の列は W_ij）を次の近似と比べる:
//!
//! - 0 次: 燃料のリングをすべて Σ̂ とした R
//! - 1 次: 0 次に Σ δ_i ∂R/∂δ_i を足したもの（中心差分、幅 1e-3）
//! - 2 次: 1 次に ½ Σ δ_i δ_j ∂²R/∂δ_i∂δ_j を足したもの
//!
//! 微分は被覆管の全断面積 0.9 /cm で求め、被覆管を 0.25 と 2.0 /cm にした行は、0 次の値だけをその被覆管で求める
//! （微分を被覆管の 1 点だけで持つ案の誤差）。差は ADR-0156 の測り方: 真空でない列は |ΔR_ij| / V_i、真空の列はその列の
//! 最大の値との比。δ の形は、Gd の割合 g の燃料で Gd の数密度が内側から 0.95〜0.08 倍（外側から燃えた形）のときの
//! ln(1 − g + g p_i / 0.66)（0.66 は割合の平均の目安）を中心に合わせたもの。

use xlc_lattice::cp::Cylinder;

/// 燃料のリングの数
const FUEL: usize = 8;
/// 中心差分の幅（ln Σ）
const STEP: f64 = 1e-3;
const WATER: f64 = 1.4;

struct Cell {
    cylinder: Cylinder,
    volumes: Vec<f64>,
}

impl Cell {
    fn rings(&self) -> usize {
        self.volumes.len()
    }

    /// 燃料のリングの ln Σ の中心からの差 `d`、中心 `center`、被覆管 `clad` での R（[リング × リング]）
    fn collisions(&self, center: f64, d: &[f64], clad: f64) -> Vec<f64> {
        let mut totals: Vec<f64> = d.iter().map(|x| center * x.exp()).collect();
        totals.push(0.0);
        totals.push(clad);
        totals.resize(self.rings(), WATER);
        let w = self.cylinder.tracks(&totals).expect("finite totals");
        let n = self.rings();
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
        r
    }

    /// `got` と `exact` の差（モジュールの説明の測り方）
    fn difference(&self, got: &[f64], exact: &[f64]) -> f64 {
        let n = self.rings();
        let mut worst: f64 = 0.0;
        for i in 0..n {
            for j in 0..n {
                let d = (got[i * n + j] - exact[i * n + j]).abs();
                // 燃料の後ろの列 FUEL は真空のギャップ
                let scale = if j == FUEL {
                    (0..n).map(|k| exact[k * n + j]).fold(0.0, f64::max)
                } else {
                    self.volumes[i]
                };
                worst = worst.max(d / scale);
            }
        }
        worst
    }
}

fn main() {
    let mut radii = vec![
        0.1448154688,
        0.2048,
        0.2508277497,
        0.2896309376,
        0.3238172324,
        0.3547240054,
        0.3831457164,
        0.4096,
    ];
    radii.extend([0.418, 0.475, 0.5317, 0.5822, 0.63, 0.6916, 0.75, 0.80]);
    let cylinder = Cylinder::new(radii).expect("increasing radii");
    let volumes = cylinder.volumes();
    let cell = Cell { cylinder, volumes };
    let n = cell.rings();
    let fuel_volume: f64 = cell.volumes[..FUEL].iter().sum();
    let left = [0.95, 0.93, 0.9, 0.85, 0.75, 0.55, 0.3, 0.08];
    let clad0 = 0.9;
    println!("center gd amplitude clad zero first second");
    for center in [0.4, 3.0, 30.0, 300.0] {
        let zero = vec![0.0; FUEL];
        let r0 = cell.collisions(center, &zero, clad0);
        let shifted = |pairs: &[(usize, f64)]| {
            let mut d = zero.clone();
            for &(k, x) in pairs {
                d[k] = x;
            }
            cell.collisions(center, &d, clad0)
        };
        let mut gradient = vec![vec![0.0; n * n]; FUEL];
        let mut hessian = vec![vec![vec![0.0; n * n]; FUEL]; FUEL];
        for a in 0..FUEL {
            let (p, q) = (shifted(&[(a, STEP)]), shifted(&[(a, -STEP)]));
            for x in 0..n * n {
                gradient[a][x] = (p[x] - q[x]) / (2.0 * STEP);
                hessian[a][a][x] = (p[x] - 2.0 * r0[x] + q[x]) / (STEP * STEP);
            }
            for b in 0..a {
                let pp = shifted(&[(a, STEP), (b, STEP)]);
                let pm = shifted(&[(a, STEP), (b, -STEP)]);
                let mp = shifted(&[(a, -STEP), (b, STEP)]);
                let mm = shifted(&[(a, -STEP), (b, -STEP)]);
                for x in 0..n * n {
                    let v = (pp[x] - pm[x] - mp[x] + mm[x]) / (4.0 * STEP * STEP);
                    hessian[a][b][x] = v;
                    hessian[b][a][x] = v;
                }
            }
        }
        for g in [0.01, 0.03, 0.1, 0.3, 0.6, 0.9] {
            let raw: Vec<f64> = left
                .iter()
                .map(|p: &f64| (1.0 - g + g * p / 0.66).ln())
                .collect();
            let mean: f64 = raw
                .iter()
                .zip(&cell.volumes)
                .map(|(x, v)| x * v)
                .sum::<f64>()
                / fuel_volume;
            let d: Vec<f64> = raw.iter().map(|x| x - mean).collect();
            let amplitude = d.iter().fold(0.0f64, |a, b| a.max(b.abs()));
            for clad in [clad0, 0.25, 2.0] {
                let exact = cell.collisions(center, &d, clad);
                let base = cell.collisions(center, &zero, clad);
                let first: Vec<f64> = (0..n * n)
                    .map(|x| base[x] + (0..FUEL).map(|a| d[a] * gradient[a][x]).sum::<f64>())
                    .collect();
                let second: Vec<f64> = (0..n * n)
                    .map(|x| {
                        first[x]
                            + 0.5
                                * (0..FUEL)
                                    .flat_map(|a| (0..FUEL).map(move |b| (a, b)))
                                    .map(|(a, b)| d[a] * d[b] * hessian[a][b][x])
                                    .sum::<f64>()
                    })
                    .collect();
                println!(
                    "{center} {g} {amplitude:.3} {clad} {:.2e} {:.2e} {:.2e}",
                    cell.difference(&base, &exact),
                    cell.difference(&first, &exact),
                    cell.difference(&second, &exact),
                );
            }
        }
    }
}
