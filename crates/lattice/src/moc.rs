//! 平坦線源の特性曲線法（MOC）による固有値計算。

use std::f64::consts::PI;

use crate::material::Material;
use crate::tracking::Tracking;

/// 極角求積（Tabuchi-Yamamoto 3点、0 < θ < π/2 の半球で重みの和が 1）
const POLAR_SIN: [f64; 3] = [0.166_648, 0.537_707, 0.932_954];
const POLAR_WEIGHT: [f64; 3] = [0.046_233, 0.283_619, 0.670_148];

/// 反復の打ち切り条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Options {
    /// 実効増倍率の変化の上限
    pub k_tolerance: f64,
    /// 中性子束の相対変化の上限
    pub flux_tolerance: f64,
    pub max_iterations: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            k_tolerance: 1.0e-7,
            flux_tolerance: 1.0e-6,
            max_iterations: 20_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    pub k: f64,
    /// [領域][群] の領域平均スカラー中性子束（全体の生成率が 1 になるよう規格化）
    pub flux: Vec<Vec<f64>>,
    /// トラックから推定した領域の体積
    pub volumes: Vec<f64>,
    pub iterations: usize,
}

/// `materials[i]` は領域 i の物質。
pub fn solve(
    tracking: &Tracking,
    materials: &[&Material],
    options: &Options,
) -> Result<Solution, String> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::PinCell;
    use crate::material::FissionTerm;

    fn two_group(total: [f64; 2], scatter: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Material {
        Material {
            total: total.to_vec(),
            scatter: scatter.iter().map(|r| r.to_vec()).collect(),
            fission: nu_fission.iter().map(|v| v / 2.5).collect(),
            nu_fission: nu_fission.to_vec(),
            fission_terms: if nu_fission.iter().any(|&v| v > 0.0) {
                vec![FissionTerm {
                    chi: vec![1.0, 0.0],
                    nu_fission: nu_fission.to_vec(),
                }]
            } else {
                vec![]
            },
        }
    }

    fn fuel() -> Material {
        two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25])
    }

    fn moderator() -> Material {
        two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0])
    }

    /// 無限媒質の k = νΣfᵀ (diag Σt - Σsᵀ)⁻¹ χ（χ = (1, 0)）
    fn infinite_medium_k(m: &Material) -> f64 {
        let a = [
            [m.total[0] - m.scatter[0][0], -m.scatter[1][0]],
            [-m.scatter[0][1], m.total[1] - m.scatter[1][1]],
        ];
        let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        let phi = [a[1][1] / det, -a[1][0] / det];
        m.nu_fission[0] * phi[0] + m.nu_fission[1] * phi[1]
    }

    #[test]
    fn homogeneous_cell_matches_infinite_medium() {
        let cell = PinCell::new(1.26, vec![0.3, 0.45]).unwrap();
        let tracking = Tracking::new(&cell, 8, 0.1).unwrap();
        let m = fuel();
        let solution = solve(&tracking, &[&m, &m, &m], &Options::default()).unwrap();
        let expected = infinite_medium_k(&m);
        assert!(
            (solution.k - expected).abs() < 1e-6,
            "{} vs {expected}",
            solution.k
        );
    }

    fn heterogeneous(n_azimuthal: usize, spacing: f64) -> Solution {
        let cell = PinCell::new(1.26, vec![0.25, 0.4096, 0.55]).unwrap();
        let tracking = Tracking::new(&cell, n_azimuthal, spacing).unwrap();
        let (f, m) = (fuel(), moderator());
        solve(&tracking, &[&f, &f, &m, &m], &Options::default()).unwrap()
    }

    #[test]
    fn heterogeneous_cell_balances_neutrons() {
        let s = heterogeneous(16, 0.05);
        let materials = [fuel(), fuel(), moderator(), moderator()];
        let (mut absorption, mut production) = (0.0, 0.0);
        for (i, m) in materials.iter().enumerate() {
            for g in 0..2 {
                absorption += s.volumes[i] * m.absorption()[g] * s.flux[i][g];
                production += s.volumes[i] * m.nu_fission[g] * s.flux[i][g];
            }
        }
        assert!((absorption - production / s.k).abs() < 1e-6 * absorption);
        // 熱中性子束は減速材で高く、燃料の中心で低い
        assert!(s.flux[3][1] > s.flux[0][1]);
    }

    #[test]
    fn refinement_converges() {
        let coarse = heterogeneous(16, 0.05);
        let fine = heterogeneous(32, 0.02);
        assert!(
            (coarse.k - fine.k).abs() < 1e-3,
            "{} vs {}",
            coarse.k,
            fine.k
        );
    }

    #[test]
    fn no_fission_is_error() {
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let tracking = Tracking::new(&cell, 8, 0.1).unwrap();
        let m = moderator();
        assert!(solve(&tracking, &[&m, &m], &Options::default()).is_err());
    }

    #[test]
    fn region_count_must_match() {
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let tracking = Tracking::new(&cell, 8, 0.1).unwrap();
        let m = fuel();
        assert!(solve(&tracking, &[&m], &Options::default()).is_err());
    }
}
