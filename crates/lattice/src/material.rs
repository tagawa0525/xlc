//! 物質の巨視的断面積。

use xs::library::Library;

/// χ の和と 1 の差の許容値
const CHI_TOLERANCE: f64 = 1.0e-6;

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
        if composition.is_empty() {
            return Err("a material needs at least one nuclide".into());
        }
        let ng = library.groups.len();
        let mut m = Self {
            total: vec![0.0; ng],
            scatter: vec![vec![0.0; ng]; ng],
            fission: vec![0.0; ng],
            nu_fission: vec![0.0; ng],
            fission_terms: Vec::new(),
        };
        let add = |sum: &mut [f64], values: &[f64], n: f64| {
            for (s, v) in sum.iter_mut().zip(values) {
                *s += n * v;
            }
        };
        for (name, density) in composition {
            if !(*density >= 0.0 && density.is_finite()) {
                return Err(format!("{name}: invalid number density {density}"));
            }
            let nuclide = library
                .nuclide(name)
                .ok_or_else(|| format!("nuclide {name} is not in the library"))?;
            add(&mut m.total, &nuclide.total, *density);
            add(&mut m.fission, &nuclide.fission, *density);
            add(&mut m.nu_fission, &nuclide.nu_fission, *density);
            for (row, from) in m.scatter.iter_mut().zip(&nuclide.scatter) {
                add(row, from, *density);
            }
            if nuclide.nu_fission.iter().any(|&v| v > 0.0) {
                m.fission_terms.push(FissionTerm {
                    chi: nuclide.chi.clone(),
                    nu_fission: nuclide.nu_fission.iter().map(|v| density * v).collect(),
                });
            }
        }
        m.validate()?;
        Ok(m)
    }

    /// 配列の大きさ、χ の規格化、νΣf と核分裂の項の和の一致を検査する。
    pub fn validate(&self) -> Result<(), String> {
        let ng = self.total.len();
        let sized = |v: &[f64]| v.len() == ng;
        if ng == 0
            || !sized(&self.fission)
            || !sized(&self.nu_fission)
            || self.scatter.len() != ng
            || !self.scatter.iter().all(|row| sized(row))
        {
            return Err(format!("cross section arrays must all have {ng} groups"));
        }
        let mut nu_fission = vec![0.0; ng];
        for (t, term) in self.fission_terms.iter().enumerate() {
            if !sized(&term.chi) || !sized(&term.nu_fission) {
                return Err(format!("fission term {t} must have {ng} groups"));
            }
            let sum: f64 = term.chi.iter().sum();
            if (sum - 1.0).abs() > CHI_TOLERANCE {
                return Err(format!("chi of fission term {t} sums to {sum}"));
            }
            for (n, v) in nu_fission.iter_mut().zip(&term.nu_fission) {
                *n += v;
            }
        }
        let consistent = nu_fission
            .iter()
            .zip(&self.nu_fission)
            .all(|(a, b)| (a - b).abs() <= 1.0e-12 * b.abs());
        if !consistent {
            return Err("nu-fission does not match the sum of fission terms".into());
        }
        Ok(())
    }

    pub fn groups(&self) -> usize {
        self.total.len()
    }

    /// 吸収断面積 Σt - Σ_g' Σs(g→g')（多重度のある反応では負になりうる）。
    pub fn absorption(&self) -> Vec<f64> {
        self.total
            .iter()
            .zip(&self.scatter)
            .map(|(t, row)| t - row.iter().sum::<f64>())
            .collect()
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
    fn fissile_nuclide_needs_normalized_chi() {
        let mut lib = library();
        lib.nuclides[0].chi = vec![0.0, 0.0];
        assert!(Material::mix(&lib, &[("F".into(), 0.01)]).is_err());
    }

    #[test]
    fn validate_detects_inconsistent_data() {
        let good = Material::mix(&library(), &[("F".into(), 0.01), ("M".into(), 0.03)]).unwrap();
        assert!(good.validate().is_ok());
        let mut short_row = good.clone();
        short_row.scatter[0].pop();
        assert!(short_row.validate().is_err());
        let mut missing_term = good.clone();
        missing_term.fission_terms.clear();
        assert!(missing_term.validate().is_err());
        let mut bad_chi = good.clone();
        bad_chi.fission_terms[0].chi = vec![0.5, 0.0];
        assert!(bad_chi.validate().is_err());
    }

    #[test]
    fn unknown_nuclide_is_error() {
        assert!(Material::mix(&library(), &[("X".into(), 0.01)]).is_err());
        assert!(Material::mix(&library(), &[("F".into(), -0.01)]).is_err());
        assert!(Material::mix(&library(), &[]).is_err());
    }
}
