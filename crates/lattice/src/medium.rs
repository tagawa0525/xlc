//! 格子を細かい群のまま中性子束と体積で均質化した媒質（ADR-0128）。

use xlc_format::medium::{Medium, MediumState};
use xlc_format::table::Axis;

use crate::fewgroup::HomogenizeError;
use crate::material::Material;
use crate::moc::{Solution, available_order};
use crate::tracking::{Boundaries, Boundary};

/// 外周がすべて鏡面反射の格子全体を、細かい群（ライブラリの群）のまま中性子束と体積で均質化した媒質。
///
/// 全断面積、散乱行列のすべてのモーメント（P0 から、真空でない物質がみな持つ次数まで）、核分裂、νΣf は領域の
/// 中性子束 × 体積で重みづけし、χ は核種の項の χ を生成率で重みづけした1つにする。反応率と生成率を保存するので、
/// この媒質の無限媒質の解は、格子の k と平均の細かい群の中性子束に一致する。P1 以上のモーメントもスカラー中性子束で
/// 重みづけするのは、鏡面反射の格子では正味の中性子流が 0 で、中性子流の重みが定まらないため（B1 の拡散係数の
/// 均質化と同じ）
pub fn homogenize_medium(
    solution: &Solution,
    materials: &[&Material],
    boundaries: &Boundaries,
) -> Result<MediumState, HomogenizeError> {
    let sides = [
        boundaries.left,
        boundaries.right,
        boundaries.top,
        boundaries.bottom,
    ];
    if sides.iter().any(|&b| b != Boundary::Reflective) {
        return Err(HomogenizeError::Medium(
            "a medium needs a lattice with every side reflective".into(),
        ));
    }
    if materials.len() != solution.flux.len() {
        return Err(HomogenizeError::RegionCount {
            materials: materials.len(),
            regions: solution.flux.len(),
        });
    }
    let regions: Vec<usize> = (0..materials.len()).collect();
    fine_medium(solution, materials, &regions).map(|(medium, _)| medium)
}

/// 領域 `regions` を細かい群のまま中性子束と体積で均質化した媒質と、[細かい群] の領域で積分した中性子束。
/// χ は生成率で規格化する（核分裂がなければ 0）
pub(crate) fn fine_medium(
    solution: &Solution,
    materials: &[&Material],
    regions: &[usize],
) -> Result<(MediumState, Vec<f64>), HomogenizeError> {
    let ng = solution.flux.first().map_or(0, Vec::len);
    let inside: Vec<&Material> = regions.iter().map(|&r| materials[r]).collect();
    let moments = 1 + available_order(&inside);
    let mut flux = vec![0.0; ng];
    let mut medium = MediumState {
        total: vec![0.0; ng],
        scatter: vec![vec![vec![0.0; ng]; ng]; moments],
        fission: vec![0.0; ng],
        nu_fission: vec![0.0; ng],
        chi: vec![0.0; ng],
    };
    for &region in regions {
        let (m, phi, volume) = (
            materials[region],
            &solution.flux[region],
            solution.volumes[region],
        );
        for g in 0..ng {
            let w = volume * phi[g];
            flux[g] += w;
            medium.total[g] += w * m.total[g];
            medium.fission[g] += w * m.fission[g];
            medium.nu_fission[g] += w * m.nu_fission[g];
            // 真空は高次のモーメントを持たない（散乱しない）
            let own = std::iter::once(&m.scatter).chain(&m.anisotropic);
            for (sum, moment) in medium.scatter.iter_mut().zip(own) {
                for (s, x) in sum[g].iter_mut().zip(&moment[g]) {
                    *s += w * x;
                }
            }
        }
        for term in &m.fission_terms {
            let rate: f64 = term
                .nu_fission
                .iter()
                .zip(phi)
                .map(|(n, f)| volume * n * f)
                .sum();
            for (c, chi) in term.chi.iter().zip(medium.chi.iter_mut()) {
                *chi += c * rate;
            }
        }
    }
    if let Some(g) = flux.iter().position(|&f| f <= 0.0) {
        return Err(HomogenizeError::NoFineFlux(g));
    }
    for (g, f) in flux.iter().enumerate() {
        medium.total[g] /= f;
        medium.fission[g] /= f;
        medium.nu_fission[g] /= f;
        for moment in &mut medium.scatter {
            for v in &mut moment[g] {
                *v /= f;
            }
        }
    }
    let produced: f64 = medium.chi.iter().sum();
    if produced > 0.0 {
        for c in &mut medium.chi {
            *c /= produced;
        }
    }
    Ok((medium, flux))
}

/// 平板の格子計算が読む媒質が、ライブラリの群境界 `bounds` と、平板の分岐の軸 `axes`（名前、点、密度、B-10 の
/// 割合。分岐がなければ空）に一致することを確かめる。分岐の状態ごとに、媒質の同じ状態を使う（ADR-0128）
pub fn check_medium(medium: &Medium, bounds: &[f64], axes: &[Axis]) -> Result<(), String> {
    if medium.bounds != bounds {
        return Err(format!(
            "the medium has the group bounds {:?}, the library {bounds:?}",
            medium.bounds
        ));
    }
    if medium.axes != axes {
        return Err(format!(
            "the medium has the branch axes {:?}, the lattice {axes:?}; make the medium with the \
             same branches",
            medium.axes
        ));
    }
    Ok(())
}

/// 媒質 `media`（名前）の散乱のモーメントの次数が、ほかの真空でない物質と同じことを確かめる（ADR-0128）。`named` は
/// 状態の物質の（名前、物質）
pub fn check_orders(named: &[(String, Material)], media: &[String]) -> Result<(), String> {
    for (name, medium) in named.iter().filter(|(n, _)| media.contains(n)) {
        let order = medium.anisotropic.len();
        if let Some((other, m)) = named
            .iter()
            .find(|(_, m)| !m.is_void() && m.anisotropic.len() != order)
        {
            return Err(format!(
                "the medium `{name}` has scattering moments up to P{order}, the material `{other}` \
                 up to P{}",
                m.anisotropic.len()
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
