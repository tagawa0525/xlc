//! 3次元直交座標の多群ノード展開法（NEM）。
//!
//! 各ノード・各方向で横方向に積分した1次元の拡散方程式を、ノード平均と4次までの多項式
//! f1 = ξ, f2 = 3ξ² - 1/4, f3 = ξ³ - ξ/4, f4 = (ξ² - 1/20)(ξ² - 1/4)（ξ ∈ [-1/2, 1/2]）で解く。
//! 3次・4次の係数は f1, f2 で重み付けした残差から、横方向の漏れは隣のノードを含む
//! 2次式で近似する。ノード間は部分中性子流でつなぐ。

use lattice::fewgroup::FewGroup;

/// 外周の境界条件。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// 正味の中性子流 0（対称面）
    Reflective,
    /// 境界での中性子束 0
    ZeroFlux,
}

/// 直方体の炉心のノード分割。幅は cm。
#[derive(Debug, Clone, PartialEq)]
pub struct Mesh {
    pub dx: Vec<f64>,
    pub dy: Vec<f64>,
    pub dz: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Problem {
    pub mesh: Mesh,
    pub materials: Vec<FewGroup>,
    /// ノードごとの物質番号。並びは x が最も速く、次に y、z
    pub node_material: Vec<usize>,
    /// [-x, +x, -y, +y, -z, +z]
    pub boundaries: [Boundary; 6],
}

/// 反復の打ち切り条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Options {
    pub k_tolerance: f64,
    pub flux_tolerance: f64,
    pub max_iterations: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            k_tolerance: 1.0e-8,
            flux_tolerance: 1.0e-7,
            max_iterations: 10_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    pub k: f64,
    /// [ノード][群] のノード平均中性子束（全体の生成率が 1 になるよう規格化）
    pub flux: Vec<Vec<f64>>,
    /// 外周から外へ出る正味の中性子数（群ごと）
    pub leakage: Vec<f64>,
    pub iterations: usize,
}

pub fn solve(problem: &Problem, options: &Options) -> Result<Solution, String> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// 2群の均質な物質（χ = (1, 0)、上方散乱なし）
    fn material() -> FewGroup {
        let (d1, d2) = (1.4, 0.4);
        let (a1, a2) = (0.010, 0.080);
        let s12 = 0.017;
        let (nf1, nf2) = (0.006, 0.11);
        let total = [1.0 / (3.0 * d1), 1.0 / (3.0 * d2)];
        FewGroup {
            bounds: vec![1.0e7, 0.625, 1.0e-5],
            k_inf: 0.0,
            total: total.to_vec(),
            absorption: vec![a1, a2],
            nu_fission: vec![nf1, nf2],
            fission: vec![nf1 / 2.4, nf2 / 2.4],
            chi: vec![1.0, 0.0],
            diffusion: vec![d1, d2],
            scatter: vec![vec![total[0] - a1 - s12, s12], vec![0.0, total[1] - a2]],
        }
    }

    /// 2群・χ = (1, 0)・上方散乱なしの均質な体系の k（バックリング b2）
    fn analytic_k(m: &FewGroup, b2: f64) -> f64 {
        let s12 = m.scatter[0][1];
        let r1 = m.absorption[0] + s12 + m.diffusion[0] * b2;
        let r2 = m.absorption[1] + m.diffusion[1] * b2;
        (m.nu_fission[0] + m.nu_fission[1] * s12 / r2) / r1
    }

    fn cuboid(size: [f64; 3], nodes: [usize; 3], boundaries: [Boundary; 6]) -> Problem {
        let split = |l: f64, n: usize| vec![l / n as f64; n];
        Problem {
            mesh: Mesh {
                dx: split(size[0], nodes[0]),
                dy: split(size[1], nodes[1]),
                dz: split(size[2], nodes[2]),
            },
            materials: vec![material()],
            node_material: vec![0; nodes.iter().product()],
            boundaries,
        }
    }

    #[test]
    fn reflective_box_gives_k_infinity() {
        let p = cuboid([30.0, 30.0, 40.0], [2, 3, 2], [Boundary::Reflective; 6]);
        let s = solve(&p, &Options::default()).unwrap();
        let expected = analytic_k(&material(), 0.0);
        assert!((s.k - expected).abs() < 1e-8, "{} vs {expected}", s.k);
    }

    #[test]
    fn bare_cuboid_matches_analytic_buckling() {
        let size = [60.0, 80.0, 100.0];
        let b2: f64 = size.iter().map(|l| (PI / l).powi(2)).sum();
        let expected = analytic_k(&material(), b2);
        let p = cuboid(size, [4, 4, 4], [Boundary::ZeroFlux; 6]);
        let s = solve(&p, &Options::default()).unwrap();
        assert!(
            (s.k - expected).abs() < 2e-4 * expected,
            "{} vs {expected}",
            s.k
        );
        // 細かくすると解析解に近づく
        let fine = solve(
            &cuboid(size, [8, 8, 8], [Boundary::ZeroFlux; 6]),
            &Options::default(),
        )
        .unwrap();
        assert!((fine.k - expected).abs() < (s.k - expected).abs());
    }

    #[test]
    fn symmetric_half_core_matches_full_core() {
        let full = cuboid([60.0, 60.0, 60.0], [6, 4, 4], [Boundary::ZeroFlux; 6]);
        let mut half_boundaries = [Boundary::ZeroFlux; 6];
        half_boundaries[0] = Boundary::Reflective;
        let half = cuboid([30.0, 60.0, 60.0], [3, 4, 4], half_boundaries);
        let (f, h) = (
            solve(&full, &Options::default()).unwrap(),
            solve(&half, &Options::default()).unwrap(),
        );
        assert!((f.k - h.k).abs() < 1e-7, "{} vs {}", f.k, h.k);
    }

    #[test]
    fn neutron_balance_holds() {
        let p = cuboid([60.0, 80.0, 100.0], [3, 4, 5], [Boundary::ZeroFlux; 6]);
        let s = solve(&p, &Options::default()).unwrap();
        let m = material();
        let volume = 60.0 * 80.0 * 100.0 / 60.0;
        let (mut absorption, mut production) = (0.0, 0.0);
        for phi in &s.flux {
            for g in 0..2 {
                absorption += volume * m.absorption[g] * phi[g];
                production += volume * m.nu_fission[g] * phi[g];
            }
        }
        let leakage: f64 = s.leakage.iter().sum();
        assert!(leakage > 0.0);
        assert!(
            (production / s.k - absorption - leakage).abs() < 1e-6 * production,
            "{production} / {} vs {absorption} + {leakage}",
            s.k
        );
    }

    #[test]
    fn inconsistent_problem_is_error() {
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
        p.node_material.pop();
        assert!(solve(&p, &Options::default()).is_err());
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
        p.node_material[0] = 1;
        assert!(solve(&p, &Options::default()).is_err());
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
        p.mesh.dx[0] = 0.0;
        assert!(solve(&p, &Options::default()).is_err());
    }
}
