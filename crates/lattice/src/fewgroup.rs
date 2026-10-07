//! 均質化した少数群定数（`lattice` から `core` に渡すファイル）。
//!
//! フィールド名付きの map 形式の MessagePack（全体計画の「ファイルの約束事」）:
//! `format` ("xlc-fewgroup"), `bounds`（群境界 (eV)、降順）, `k_inf`, `total`, `absorption`, `nu_fission`,
//! `fission`, `chi`, `diffusion`（G 個ずつ）, `scatter`（[入射の群][出射の群]）。

use serde::{Deserialize, Serialize};
use xlc_xs::msgpack::{self, MsgpackError};
use xlc_xs::multigroup::GroupStructure;

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

/// 少数群定数のファイルの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum FewGroupError {
    #[error(transparent)]
    Msgpack(#[from] MsgpackError),
    #[error("expected format xlc-fewgroup, found {0}")]
    Format(String),
    #[error("{0}")]
    Invalid(String),
}

/// 少数群への縮約のエラー。
#[derive(Debug, thiserror::Error)]
pub enum HomogenizeError {
    #[error("{materials} materials for {regions} regions")]
    RegionCount { materials: usize, regions: usize },
    #[error("few-group bounds must be bounds of the library groups")]
    NotLibraryBounds,
    #[error("few-group bounds must be descending and span the library groups")]
    Span,
    #[error("a few group has no flux")]
    NoFlux,
}

/// ファイルの中身（フィールド名付きの map 形式で書く）。
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct FewGroupFile {
    format: String,
    bounds: Vec<f64>,
    k_inf: f64,
    total: Vec<f64>,
    absorption: Vec<f64>,
    nu_fission: Vec<f64>,
    fission: Vec<f64>,
    chi: Vec<f64>,
    diffusion: Vec<f64>,
    scatter: Vec<Vec<f64>>,
}

impl FewGroup {
    pub fn groups(&self) -> usize {
        self.total.len()
    }

    /// MessagePack（フィールド名付きの map 形式）で書く。
    pub fn to_msgpack(&self) -> Result<Vec<u8>, MsgpackError> {
        msgpack::to_named(&FewGroupFile {
            format: FORMAT.into(),
            bounds: self.bounds.clone(),
            k_inf: self.k_inf,
            total: self.total.clone(),
            absorption: self.absorption.clone(),
            nu_fission: self.nu_fission.clone(),
            fission: self.fission.clone(),
            chi: self.chi.clone(),
            diffusion: self.diffusion.clone(),
            scatter: self.scatter.clone(),
        })
    }

    /// MessagePack から読み、群境界の順序、配列の長さ、値が有限であることを確かめる。
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self, FewGroupError> {
        let file: FewGroupFile = msgpack::from_value(msgpack::read_value(bytes)?)?;
        if file.format != FORMAT {
            return Err(FewGroupError::Format(file.format));
        }
        let invalid = |message: String| Err(FewGroupError::Invalid(message));
        let b = &file.bounds;
        if b.len() < 2
            || !b.windows(2).all(|w| w[0] > w[1])
            || !b.iter().all(|&x| x > 0.0 && x.is_finite())
        {
            return invalid(format!("bounds must be positive and descending, got {b:?}"));
        }
        let ng = b.len() - 1;
        let vectors = [
            ("total", &file.total),
            ("absorption", &file.absorption),
            ("nu_fission", &file.nu_fission),
            ("fission", &file.fission),
            ("chi", &file.chi),
            ("diffusion", &file.diffusion),
        ];
        for (name, v) in vectors
            .into_iter()
            .chain(file.scatter.iter().map(|row| ("scatter", row)))
        {
            if v.len() != ng {
                return invalid(format!("`{name}` has {} values for {ng} groups", v.len()));
            }
            if let Some(bad) = v.iter().find(|x| !x.is_finite()) {
                return invalid(format!("`{name}` has a non-finite value {bad}"));
            }
        }
        if file.scatter.len() != ng {
            return invalid(format!(
                "`scatter` has {} rows for {ng} groups",
                file.scatter.len()
            ));
        }
        if !(file.k_inf > 0.0 && file.k_inf.is_finite()) {
            return invalid(format!("invalid k-inf {}", file.k_inf));
        }
        Ok(Self {
            bounds: file.bounds,
            k_inf: file.k_inf,
            total: file.total,
            absorption: file.absorption,
            nu_fission: file.nu_fission,
            fission: file.fission,
            chi: file.chi,
            diffusion: file.diffusion,
            scatter: file.scatter,
        })
    }
}

const FORMAT: &str = "xlc-fewgroup";

/// MOC の解を使い、領域・群ごとの反応率を保存するようにセル全体で均質化し、
/// `coarse_bounds`（`fine` の境界の一部、降順）の少数群に縮約する。
/// 拡散係数は D = 1/(3Σt) とする（M0 の近似）。
pub fn homogenize(
    solution: &Solution,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
) -> Result<FewGroup, HomogenizeError> {
    let ng = fine.len();
    if materials.len() != solution.flux.len() {
        return Err(HomogenizeError::RegionCount {
            materials: materials.len(),
            regions: solution.flux.len(),
        });
    }
    // 少数群の境界を細かい群の境界の位置に対応させる
    let positions = coarse_bounds
        .iter()
        .map(|b| fine.bounds().iter().position(|f| f == b))
        .collect::<Option<Vec<usize>>>()
        .ok_or(HomogenizeError::NotLibraryBounds)?;
    let covers = positions.first() == Some(&0) && positions.last() == Some(&ng);
    if !covers || positions.windows(2).any(|w| w[0] >= w[1]) {
        return Err(HomogenizeError::Span);
    }
    let coarse = positions.len() - 1;
    let mut group_of = vec![0; ng];
    for (big, w) in positions.windows(2).enumerate() {
        group_of[w[0]..w[1]].fill(big);
    }

    let mut flux = vec![0.0; coarse];
    let mut total = vec![0.0; coarse];
    let mut nu_fission = vec![0.0; coarse];
    let mut fission = vec![0.0; coarse];
    let mut scatter = vec![vec![0.0; coarse]; coarse];
    let mut chi = vec![0.0; coarse];
    for ((m, phi), volume) in materials.iter().zip(&solution.flux).zip(&solution.volumes) {
        for g in 0..ng {
            let w = volume * phi[g];
            let big = group_of[g];
            flux[big] += w;
            total[big] += w * m.total[g];
            nu_fission[big] += w * m.nu_fission[g];
            fission[big] += w * m.fission[g];
            for (to, s) in m.scatter[g].iter().enumerate() {
                scatter[big][group_of[to]] += w * s;
            }
        }
        for term in &m.fission_terms {
            let rate: f64 = term
                .nu_fission
                .iter()
                .zip(phi)
                .map(|(n, f)| volume * n * f)
                .sum();
            for (g, c) in term.chi.iter().enumerate() {
                chi[group_of[g]] += c * rate;
            }
        }
    }
    if flux.iter().any(|&f| f <= 0.0) {
        return Err(HomogenizeError::NoFlux);
    }
    let per_flux =
        |values: Vec<f64>| -> Vec<f64> { values.iter().zip(&flux).map(|(v, f)| v / f).collect() };
    let total = per_flux(total);
    let scatter: Vec<Vec<f64>> = scatter
        .into_iter()
        .zip(&flux)
        .map(|(row, f)| row.iter().map(|v| v / f).collect())
        .collect();
    let produced: f64 = chi.iter().sum();
    Ok(FewGroup {
        bounds: coarse_bounds.to_vec(),
        k_inf: solution.k,
        absorption: total
            .iter()
            .zip(&scatter)
            .map(|(t, row)| t - row.iter().sum::<f64>())
            .collect(),
        diffusion: total.iter().map(|t| 1.0 / (3.0 * t)).collect(),
        total,
        nu_fission: per_flux(nu_fission),
        fission: per_flux(fission),
        chi: if produced > 0.0 {
            chi.iter().map(|c| c / produced).collect()
        } else {
            chi
        },
        scatter,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};
    use crate::material::FissionTerm;
    use crate::moc::{Options, solve};
    use crate::tracking::Tracking;

    fn two_group(total: [f64; 2], scatter: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Material {
        Material {
            total: total.to_vec(),
            scatter: scatter.iter().map(|r| r.to_vec()).collect(),
            anisotropic: Vec::new(),
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
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 16, 0.05).unwrap();
        let refs: Vec<&Material> = materials.iter().collect();
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        (
            solution,
            materials,
            GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
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

    fn sample() -> FewGroup {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        homogenize(&solution, &refs, &fine, fine.bounds()).unwrap()
    }

    fn sample_file() -> FewGroupFile {
        let fg = sample();
        FewGroupFile {
            format: "xlc-fewgroup".into(),
            bounds: fg.bounds,
            k_inf: fg.k_inf,
            total: fg.total,
            absorption: fg.absorption,
            nu_fission: fg.nu_fission,
            fission: fg.fission,
            chi: fg.chi,
            diffusion: fg.diffusion,
            scatter: fg.scatter,
        }
    }

    fn read_modified(change: impl FnOnce(&mut FewGroupFile)) -> Result<FewGroup, FewGroupError> {
        let mut file = sample_file();
        change(&mut file);
        FewGroup::from_msgpack(&rmp_serde::to_vec_named(&file).unwrap())
    }

    #[test]
    fn round_trips_bit_exactly() {
        let fg = sample();
        assert_eq!(
            FewGroup::from_msgpack(&fg.to_msgpack().unwrap()).unwrap(),
            fg
        );
    }

    #[test]
    fn invalid_files_are_rejected() {
        assert!(read_modified(|f| f.format = "xlc-multigroup".into()).is_err());
        assert!(
            read_modified(|f| {
                f.diffusion.pop();
            })
            .is_err()
        );
        assert!(
            read_modified(|f| {
                f.scatter[0].pop();
            })
            .is_err()
        );
        assert!(read_modified(|f| f.bounds.reverse()).is_err());
        assert!(read_modified(|f| f.total[0] = f64::NAN).is_err());
        assert!(read_modified(|f| f.k_inf = f64::INFINITY).is_err());
        // 配列形式、末尾の余り
        let array = rmp_serde::to_vec(&sample_file()).unwrap();
        assert!(matches!(
            FewGroup::from_msgpack(&array),
            Err(FewGroupError::Msgpack(MsgpackError::NotMap))
        ));
        let mut bytes = sample().to_msgpack().unwrap();
        bytes.push(0xc0);
        assert!(FewGroup::from_msgpack(&bytes).is_err());
    }
}
