//! 物質の巨視的断面積。

use xs::library::Library;

/// 1つの核分裂性核種の寄与。核分裂源は核種ごとの χ × νΣf の和として正確に持つ。
#[derive(Debug, Clone, PartialEq)]
pub struct FissionTerm {
    pub chi: Vec<f64>,
    /// νΣf (1/cm)
    pub nu_fission: Vec<f64>,
}

/// 物質の多群巨視的断面積 (1/cm)。
#[derive(Debug, Clone, PartialEq)]
pub struct Material {
    pub total: Vec<f64>,
    /// [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
    pub fission: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub fission_terms: Vec<FissionTerm>,
}

impl Material {
    /// 核種名と個数密度 (1/(barn·cm)) の組から作る。
    pub fn mix(library: &Library, composition: &[(String, f64)]) -> Result<Self, String> {
        todo!()
    }

    pub fn groups(&self) -> usize {
        self.total.len()
    }

    /// 吸収断面積 Σt - Σ_g' Σs(g→g')（多重度のある反応では負になりうる）。
    pub fn absorption(&self) -> Vec<f64> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use xs::library::Nuclide;
    use xs::multigroup::GroupStructure;

    fn library() -> Library {
        let nuclide = |name: &str, fissile: bool| Nuclide {
            name: name.into(),
            za: 0.0,
            awr: 1.0,
            temperature: 293.6,
            total: vec![10.0, 20.0],
            fission: if fissile {
                vec![1.0, 2.0]
            } else {
                vec![0.0, 0.0]
            },
            nu_fission: if fissile {
                vec![2.5, 5.0]
            } else {
                vec![0.0, 0.0]
            },
            chi: if fissile {
                vec![1.0, 0.0]
            } else {
                vec![0.0, 0.0]
            },
            scatter: vec![vec![6.0, 2.0], vec![0.5, 15.0]],
        };
        Library {
            groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]),
            nuclides: vec![nuclide("F", true), nuclide("M", false)],
        }
    }

    #[test]
    fn mixes_number_densities() {
        let m = Material::mix(&library(), &[("F".into(), 0.01), ("M".into(), 0.03)]).unwrap();
        assert_eq!(m.groups(), 2);
        assert!((m.total[1] - 0.8).abs() < 1e-15);
        assert!((m.scatter[0][1] - 0.08).abs() < 1e-15);
        assert!((m.nu_fission[1] - 0.05).abs() < 1e-15);
        assert_eq!(m.fission_terms.len(), 1);
        assert!((m.fission_terms[0].nu_fission[0] - 0.025).abs() < 1e-15);
        let a = m.absorption();
        assert!((a[0] - (0.4 - 0.32)).abs() < 1e-15);
    }

    #[test]
    fn unknown_nuclide_is_error() {
        assert!(Material::mix(&library(), &[("X".into(), 0.01)]).is_err());
        assert!(Material::mix(&library(), &[("F".into(), -0.01)]).is_err());
        assert!(Material::mix(&library(), &[]).is_err());
    }
}
