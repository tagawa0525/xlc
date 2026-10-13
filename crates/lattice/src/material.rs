//! 物質の巨視的断面積。

use std::collections::BTreeMap;

use xlc_format::library::{Library, LookupError, Nuclide};
use xlc_format::medium::MediumState;

use crate::input::MaterialInput;

/// χ の和と 1 の差の許容値
const CHI_TOLERANCE: f64 = 1.0e-6;

/// 1つの核分裂性核種の寄与。核分裂源は核種ごとの χ × νΣf の和として正確に持つ。
#[derive(Debug, Clone, PartialEq)]
pub struct FissionTerm {
    pub chi: Vec<f64>,
    /// νΣf (1/cm)
    pub nu_fission: Vec<f64>,
}

/// 物質 `parts`（成分と体積の割合）を混ぜた（核種、温度 (K)、個数密度 (1/(barn·cm))）の並び（核種の名前、温度の順）。
/// 個数密度は成分の個数密度 × 割合の和で、同じ温度の同じ核種は足し、違う温度は別の項にする（ADR-0128）
pub fn mixed_nuclides(parts: &[(&MaterialInput, f64)]) -> Vec<(String, f64, f64)> {
    // 温度は正なので、ビット列の順は値の順と同じ
    let mut sum: BTreeMap<(String, u64), f64> = BTreeMap::new();
    for (part, fraction) in parts {
        for (name, density) in &part.composition {
            *sum.entry((name.clone(), part.temperature.to_bits()))
                .or_insert(0.0) += density * fraction;
        }
    }
    sum.into_iter()
        .map(|((name, t), density)| (name, f64::from_bits(t), density))
        .collect()
}

/// 物質の多群巨視的断面積 (1/cm)。
#[derive(Debug, Clone, PartialEq)]
pub struct Material {
    pub total: Vec<f64>,
    /// P0 の散乱行列 [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
    /// P1 以上の散乱行列 [l - 1][入射の群][出射の群]（l = 1..=L）。空なら等方散乱
    pub anisotropic: Vec<Vec<Vec<f64>>>,
    pub fission: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub fission_terms: Vec<FissionTerm>,
}

/// 物質の断面積を作る・確かめるときのエラー。
#[derive(Debug, thiserror::Error)]
pub enum MaterialError {
    #[error("a material needs at least one nuclide")]
    NoNuclides,
    #[error("{nuclide}: invalid number density {density}")]
    Density { nuclide: String, density: f64 },
    #[error(transparent)]
    Lookup(#[from] LookupError),
    #[error("cross section arrays must all have {0} groups")]
    Groups(usize),
    #[error("fission term {term} must have {groups} groups")]
    FissionTermGroups { term: usize, groups: usize },
    #[error("chi of fission term {term} sums to {sum}")]
    Chi { term: usize, sum: f64 },
    #[error("nu-fission does not match the sum of fission terms")]
    NuFission,
    #[error("{nuclide} has scattering moments up to P{order}, the others up to P{expected}")]
    ScatteringOrder {
        nuclide: String,
        order: usize,
        expected: usize,
    },
}

impl Material {
    /// 核種名と個数密度 (1/(barn·cm)) の組から、温度 `temperature` (K) の無限希釈の群定数で作る。
    pub fn mix(
        library: &Library,
        composition: &[(String, f64)],
        temperature: f64,
    ) -> Result<Self, MaterialError> {
        let nuclides = composition
            .iter()
            .map(|(name, density)| Ok((library.nuclide(name, temperature)?, *density)))
            .collect::<Result<Vec<_>, MaterialError>>()?;
        Self::from_nuclides(&nuclides, library.groups.len())
    }

    /// 物質 `parts` を体積の割合で混ぜた物質（ADR-0128）。成分ごとの温度のライブラリの値を使う（`mixed_nuclides`）。
    pub fn mixture(
        library: &Library,
        parts: &[(&MaterialInput, f64)],
    ) -> Result<Self, MaterialError> {
        let nuclides = mixed_nuclides(parts)
            .into_iter()
            .map(|(name, temperature, density)| Ok((library.nuclide(&name, temperature)?, density)))
            .collect::<Result<Vec<_>, MaterialError>>()?;
        Self::from_nuclides(&nuclides, library.groups.len())
    }

    /// 細かい群の均質な媒質（ADR-0128）から作る。P0 の散乱を `scatter`、P1 以上を `anisotropic` に、核分裂を
    /// 媒質の χ の1つの項に持つ（核分裂がなければ項はない）。
    pub fn from_medium(medium: &MediumState) -> Result<Self, MaterialError> {
        let ng = medium.groups();
        let Some((p0, higher)) = medium.scatter.split_first() else {
            return Err(MaterialError::Groups(ng));
        };
        let fissile = medium.nu_fission.iter().any(|&v| v > 0.0);
        let m = Self {
            total: medium.total.clone(),
            scatter: p0.clone(),
            anisotropic: higher.to_vec(),
            fission: medium.fission.clone(),
            nu_fission: medium.nu_fission.clone(),
            fission_terms: if fissile {
                vec![FissionTerm {
                    chi: medium.chi.clone(),
                    nu_fission: medium.nu_fission.clone(),
                }]
            } else {
                Vec::new()
            },
        };
        m.validate()?;
        Ok(m)
    }

    /// 核種の群定数 (barn) と個数密度 (1/(barn·cm)) の組から作る。
    pub fn from_nuclides(nuclides: &[(Nuclide, f64)], ng: usize) -> Result<Self, MaterialError> {
        if nuclides.is_empty() {
            return Err(MaterialError::NoNuclides);
        }
        let order = nuclides[0].0.anisotropic.len();
        let mut m = Self {
            total: vec![0.0; ng],
            scatter: vec![vec![0.0; ng]; ng],
            anisotropic: vec![vec![vec![0.0; ng]; ng]; order],
            fission: vec![0.0; ng],
            nu_fission: vec![0.0; ng],
            fission_terms: Vec::new(),
        };
        let add = |sum: &mut [f64], values: &[f64], n: f64| {
            for (s, v) in sum.iter_mut().zip(values) {
                *s += n * v;
            }
        };
        for (nuclide, density) in nuclides {
            if !(*density >= 0.0 && density.is_finite()) {
                return Err(MaterialError::Density {
                    nuclide: nuclide.name.clone(),
                    density: *density,
                });
            }
            if nuclide.anisotropic.len() != order {
                return Err(MaterialError::ScatteringOrder {
                    nuclide: nuclide.name.clone(),
                    order: nuclide.anisotropic.len(),
                    expected: order,
                });
            }
            add(&mut m.total, &nuclide.total, *density);
            add(&mut m.fission, &nuclide.fission, *density);
            add(&mut m.nu_fission, &nuclide.nu_fission, *density);
            let sums = std::iter::once(&mut m.scatter).chain(&mut m.anisotropic);
            let moments = std::iter::once(&nuclide.scatter).chain(&nuclide.anisotropic);
            for (sum, moment) in sums.zip(moments) {
                for (row, from) in sum.iter_mut().zip(moment) {
                    add(row, from, *density);
                }
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
    pub fn validate(&self) -> Result<(), MaterialError> {
        let ng = self.total.len();
        let sized = |v: &[f64]| v.len() == ng;
        if ng == 0
            || !sized(&self.fission)
            || !sized(&self.nu_fission)
            || !std::iter::once(&self.scatter)
                .chain(&self.anisotropic)
                .all(|m| m.len() == ng && m.iter().all(|row| sized(row)))
        {
            return Err(MaterialError::Groups(ng));
        }
        let mut nu_fission = vec![0.0; ng];
        for (t, term) in self.fission_terms.iter().enumerate() {
            if !sized(&term.chi) || !sized(&term.nu_fission) {
                return Err(MaterialError::FissionTermGroups {
                    term: t,
                    groups: ng,
                });
            }
            let sum: f64 = term.chi.iter().sum();
            if (sum - 1.0).abs() > CHI_TOLERANCE {
                return Err(MaterialError::Chi { term: t, sum });
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
            return Err(MaterialError::NuFission);
        }
        Ok(())
    }

    pub fn groups(&self) -> usize {
        self.total.len()
    }

    /// 真空（どの群の断面積も 0）。
    pub fn void(ng: usize) -> Self {
        Self {
            total: vec![0.0; ng],
            scatter: vec![vec![0.0; ng]; ng],
            anisotropic: Vec::new(),
            fission: vec![0.0; ng],
            nu_fission: vec![0.0; ng],
            fission_terms: Vec::new(),
        }
    }

    /// 真空か（どの反応もどの群でも 0）
    pub fn is_void(&self) -> bool {
        let zero = |v: &[f64]| v.iter().all(|&x| x == 0.0);
        zero(&self.total)
            && zero(&self.fission)
            && zero(&self.nu_fission)
            && std::iter::once(&self.scatter)
                .chain(&self.anisotropic)
                .flatten()
                .all(|row| zero(row))
            && self.fission_terms.is_empty()
    }

    /// 群 g の反応（全断面積、群 g からの各次数の散乱、核分裂、νΣf）に `factors[g]` を掛けた物質（SPH 因子を入れる）。
    pub fn scaled(&self, factors: &[f64]) -> Result<Self, MaterialError> {
        if factors.len() != self.groups() {
            return Err(MaterialError::Groups(self.groups()));
        }
        let scale = |v: &[f64]| {
            v.iter()
                .zip(factors)
                .map(|(x, f)| x * f)
                .collect::<Vec<_>>()
        };
        let rows = |matrix: &[Vec<f64>]| {
            matrix
                .iter()
                .zip(factors)
                .map(|(row, f)| row.iter().map(|x| x * f).collect())
                .collect()
        };
        Ok(Self {
            total: scale(&self.total),
            scatter: rows(&self.scatter),
            anisotropic: self.anisotropic.iter().map(|m| rows(m)).collect(),
            fission: scale(&self.fission),
            nu_fission: scale(&self.nu_fission),
            fission_terms: self
                .fission_terms
                .iter()
                .map(|t| FissionTerm {
                    chi: t.chi.clone(),
                    nu_fission: scale(&t.nu_fission),
                })
                .collect(),
        })
    }

    /// 核分裂で生まれる中性子の行列 Σ_核種 νΣf(g) χ(g') [入射の群][出射の群]
    pub fn fission_matrix(&self) -> Vec<Vec<f64>> {
        let ng = self.groups();
        let mut matrix = vec![vec![0.0; ng]; ng];
        for term in &self.fission_terms {
            for (row, n) in matrix.iter_mut().zip(&term.nu_fission) {
                for (value, chi) in row.iter_mut().zip(&term.chi) {
                    *value += n * chi;
                }
            }
        }
        matrix
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
mod tests;
