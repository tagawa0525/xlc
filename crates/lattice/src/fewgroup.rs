//! 均質化した少数群定数（`lattice` から `core` に渡すファイル）。
//!
//! ```text
//! xlc-fewgroup 1
//! groups <G>
//! bounds <G+1 個の群境界 (eV)、降順>
//! k-inf <k∞>
//! total <G 個>
//! absorption <G 個>
//! nu-fission <G 個>
//! fission <G 個>
//! chi <G 個>
//! diffusion <G 個>
//! scatter
//! <G 行 × G 個: Σ(g→g')、行が入射の群>
//! ```

use xs::multigroup::GroupStructure;

use crate::material::Material;
use crate::moc::Solution;

/// セル全体で均質化した少数群定数。断面積は 1/cm、拡散係数は cm。
#[derive(Debug, Clone, PartialEq)]
pub struct FewGroup {
    /// 群境界 (eV)、降順
    pub bounds: Vec<f64>,
    pub k_inf: f64,
    pub total: Vec<f64>,
    pub absorption: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub fission: Vec<f64>,
    pub chi: Vec<f64>,
    pub diffusion: Vec<f64>,
    /// [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
}

impl FewGroup {
    pub fn groups(&self) -> usize {
        self.total.len()
    }

    pub fn to_text(&self) -> String {
        todo!()
    }

    /// エラーは「line N: ...」の形の文字列で返す。
    pub fn parse(text: &str) -> Result<Self, String> {
        todo!()
    }
}

/// MOC の解を使い、領域・群ごとの反応率を保存するようにセル全体で均質化し、
/// `coarse_bounds`（`fine` の境界の一部、降順）の少数群に縮約する。
/// 拡散係数は D = 1/(3Σt) とする（M0 の近似）。
pub fn homogenize(
    solution: &Solution,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
) -> Result<FewGroup, String> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::PinCell;
    use crate::material::FissionTerm;
    use crate::moc::{Options, solve};
    use crate::tracking::Tracking;

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

    fn solved() -> (Solution, Vec<Material>, GroupStructure) {
        let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let moderator = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let materials = vec![fuel.clone(), moderator.clone(), moderator];
        let cell = PinCell::new(1.26, vec![0.4096, 0.55]).unwrap();
        let tracking = Tracking::new(&cell, 16, 0.05).unwrap();
        let refs: Vec<&Material> = materials.iter().collect();
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        (
            solution,
            materials,
            GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]),
        )
    }

    #[test]
    fn one_group_collapse_preserves_k() {
        // 反応率を保存するので、1群の νΣf/Σa はセルの k∞ に一致する
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let one = homogenize(&solution, &refs, &fine, &[1.0e7, 1.0e-5]).unwrap();
        assert_eq!(one.groups(), 1);
        let k = one.nu_fission[0] / one.absorption[0];
        assert!(
            (k - solution.k).abs() < 1e-6 * solution.k,
            "{k} vs {}",
            solution.k
        );
        assert_eq!(one.k_inf, solution.k);
        assert!((one.chi[0] - 1.0).abs() < 1e-12);
        assert!((one.diffusion[0] - 1.0 / (3.0 * one.total[0])).abs() < 1e-15);
    }

    #[test]
    fn two_group_constants_reproduce_k_inf() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(&solution, &refs, &fine, fine.bounds()).unwrap();
        // 均質な2群の無限媒質の k = νΣfᵀ (diag Σt - Σsᵀ)⁻¹ χ（χ = (1, 0)）
        let a = [
            [fg.total[0] - fg.scatter[0][0], -fg.scatter[1][0]],
            [-fg.scatter[0][1], fg.total[1] - fg.scatter[1][1]],
        ];
        let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        let phi = [a[1][1] / det, -a[1][0] / det];
        let k = fg.nu_fission[0] * phi[0] + fg.nu_fission[1] * phi[1];
        assert!((k - solution.k).abs() < 1e-6 * solution.k);
        for g in 0..2 {
            let removal = fg.total[g] - fg.scatter[g].iter().sum::<f64>();
            assert!((fg.absorption[g] - removal).abs() < 1e-14);
        }
    }

    #[test]
    fn coarse_bounds_must_be_fine_bounds() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        assert!(homogenize(&solution, &refs, &fine, &[1.0e7, 1.0, 1.0e-5]).is_err());
        assert!(homogenize(&solution, &refs, &fine, &[1.0e7, 0.625]).is_err());
    }

    #[test]
    fn text_round_trips() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(&solution, &refs, &fine, fine.bounds()).unwrap();
        let text = fg.to_text();
        assert!(text.starts_with("xlc-fewgroup 1\n"));
        assert_eq!(FewGroup::parse(&text).unwrap(), fg);
    }

    #[test]
    fn parse_errors_have_line_numbers() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let text = homogenize(&solution, &refs, &fine, fine.bounds())
            .unwrap()
            .to_text();
        let broken: String = text
            .lines()
            .map(|l| {
                if l.starts_with("diffusion") {
                    "diffusion 1.0"
                } else {
                    l
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let err = FewGroup::parse(&broken).unwrap_err();
        let line = text
            .lines()
            .position(|l| l.starts_with("diffusion"))
            .unwrap()
            + 1;
        assert!(err.starts_with(&format!("line {line}:")), "{err}");
        assert!(FewGroup::parse("xlc-fewgroup 9\n").is_err());
    }
}
