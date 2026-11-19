//! 無限希釈の吸収体（ADR-0261）。
//!
//! 格子計算の燃焼が輸送しない核種のうち、ライブラリの放射化の断面積（[`Library::activation`]。無限希釈の群の、チェーンの反応と
//! 核分裂）を持つものを、燃焼の領域の吸収体として数える。自己遮蔽は輸送する核種だけで解き、その後で:
//!
//! - 輸送: MOC に渡す燃焼の領域の物質の全断面積に Σᵢ Nᵢ (σ(n,γ) + σ(n,p) + σ(n,α)) を足す（[`absorption`]、[`overlay`]）。
//!   散乱と核分裂は足さない（放射化の断面積は散乱、ν、χ を持たない）。中性子を増やす (n,2n)、(n,3n)、(n,4n) も足さない
//! - 燃焼の反応率: 無限希釈の断面積と物理的な中性子束から、チェーンの反応と核分裂を数える（[`add_rates`]）

use anyhow::{anyhow, bail};
use xlc_depletion::chain::{Chain, NuclideRates};
use xlc_format::library::Library;

use crate::depletion::{DepletionRegion, same_channels, yield_weights};
use crate::material::Material;

/// 輸送で吸収に数える反応の MT: (n,γ)、(n,p)、(n,α)（中性子を 1 つ失い、出さない反応）
const ABSORPTION: [u32; 3] = [102, 103, 107];
/// 1 barn (cm²)
const BARN: f64 = 1.0e-24;

/// ライブラリが核種 `name` の放射化の断面積を温度 `temperature` (K) そのもので持つか
pub fn has_data_at(library: &Library, name: &str, temperature: f64) -> bool {
    library.activation(name, temperature).is_ok()
}

/// 無限希釈の吸収体 `nuclides`（名前と数密度 (1/(barn cm))）の、温度 `temperature` (K) の群ごとの巨視的な吸収断面積 (1/cm)。
/// 数密度が 0 の核種は引かない
pub fn absorption<'a>(
    library: &Library,
    nuclides: impl IntoIterator<Item = (&'a str, f64)>,
    temperature: f64,
) -> anyhow::Result<Vec<f64>> {
    let mut sum = vec![0.0; library.groups.len()];
    for (name, density) in nuclides {
        if !(density >= 0.0 && density.is_finite()) {
            bail!("{name}: invalid number density {density}");
        }
        if density == 0.0 {
            continue;
        }
        let data = library.activation(name, temperature)?;
        for reaction in data.reactions.iter().filter(|r| ABSORPTION.contains(&r.mt)) {
            for channel in &reaction.channels {
                if channel.cross_section.len() != sum.len() {
                    bail!(
                        "{name}: {} groups in MT{} for {} groups",
                        channel.cross_section.len(),
                        reaction.mt,
                        sum.len()
                    );
                }
                for (s, sigma) in sum.iter_mut().zip(&channel.cross_section) {
                    *s += density * sigma;
                }
            }
        }
    }
    Ok(sum)
}

/// 自己遮蔽した物質 `materials`（[領域] の番号は `region_material`）の、燃焼の領域 `regions` の物質の全断面積に、領域ごとの吸収
/// `absorption`（[燃焼の領域][群]、1/cm）を足した物質と番号。MOC に渡す物質は SPH 因子 μ を掛けた形なので、吸収にも領域の因子
/// `sph`（[領域][群]）を掛ける（物理的な中性子束 μφ との積が反応率になる）。燃焼の領域でない領域の物質と番号は変えない。同じ物質と
/// 同じ因子の領域は 1 つの物質を共有する
pub fn overlay(
    materials: &[Material],
    region_material: &[usize],
    regions: &[DepletionRegion],
    absorption: &[Vec<f64>],
    sph: &[Vec<f64>],
) -> anyhow::Result<(Vec<Material>, Vec<usize>)> {
    if absorption.len() != regions.len() || sph.len() != region_material.len() {
        bail!(
            "{} absorptions for {} depletion regions, {} SPH factors for {} regions",
            absorption.len(),
            regions.len(),
            sph.len(),
            region_material.len()
        );
    }
    let mut unique = materials.to_vec();
    let mut index = region_material.to_vec();
    // (元の物質, 燃焼の領域) ごとの、重ねた物質の（最初の領域, 番号）。因子の違う領域は別の物質にする
    let mut made: Vec<(usize, usize, usize, usize)> = Vec::new();
    for (k, (d, sigma)) in regions.iter().zip(absorption).enumerate() {
        for &r in &d.regions {
            let m = *region_material.get(r).ok_or_else(|| {
                anyhow!("depletion region {k} has region {r} outside the lattice")
            })?;
            let found = made
                .iter()
                .find(|&&(material, region, first, _)| {
                    material == m && region == k && sph[first] == sph[r]
                })
                .map(|&(_, _, _, at)| at);
            index[r] = match found {
                Some(at) => at,
                None => {
                    let scaled: Vec<f64> =
                        sigma.iter().zip(&sph[r]).map(|(s, mu)| s * mu).collect();
                    let base = materials
                        .get(m)
                        .ok_or_else(|| anyhow!("region {r} has no material {m}"))?;
                    unique.push(base.absorbing(&scaled)?);
                    made.push((m, k, r, unique.len() - 1));
                    unique.len() - 1
                }
            };
        }
    }
    Ok((unique, index))
}

/// 燃焼の領域の無限希釈の吸収体 `dilute`（チェーンの添字）の、中性子束の大きさを MOC の規格化のままとしたときの 1 群の反応率
/// (1/s、1 原子あたり) を `rates`（[チェーンの核種]）に入れる。`members` は領域の MOC の領域の（体積 (cm²)、物理的な中性子束 [群]）。
/// 反応の経路はチェーンと同じ並びであること（輸送する核種の `LibraryReactions` と同じ）
pub fn add_rates(
    chain: &Chain,
    dilute: &[usize],
    library: &Library,
    temperature: f64,
    members: &[(f64, &[f64])],
    group_bounds: &[f64],
    rates: &mut [NuclideRates],
) -> anyhow::Result<()> {
    let groups = group_bounds.len() - 1;
    let volume: f64 = members.iter().map(|(v, _)| v).sum();
    // 体積で重みを付けた領域の中性子束 [群]
    let mut flux = vec![0.0; groups];
    for (v, phi) in members {
        if phi.len() != groups {
            bail!("{} flux groups for {groups} groups", phi.len());
        }
        for (f, p) in flux.iter_mut().zip(*phi) {
            *f += v * p;
        }
    }
    for &i in dilute {
        let nuclide = &chain.nuclides[i];
        let name = &nuclide.name;
        let data = library.activation(name, temperature)?;
        if data.fission.len() != groups {
            bail!(
                "{name}: {} fission groups for {groups} groups",
                data.fission.len()
            );
        }
        for (r, reaction) in nuclide.reactions.iter().enumerate() {
            let mt = reaction.mt;
            let group = data
                .reactions
                .iter()
                .find(|g| g.mt == mt)
                .ok_or_else(|| anyhow!("{name}: the library has no chain reaction MT{mt}"))?;
            if !same_channels(group, reaction) {
                bail!("{name}: the channels of MT{mt} in the library differ from the chain");
            }
            for (c, channel) in group.channels.iter().enumerate() {
                if channel.cross_section.len() != groups {
                    bail!(
                        "{name}: {} groups in a channel of MT{mt} for {groups} groups",
                        channel.cross_section.len()
                    );
                }
                let rate: f64 = channel
                    .cross_section
                    .iter()
                    .zip(&flux)
                    .map(|(s, f)| s * f)
                    .sum();
                rates[i].reactions[r][c] = rate * BARN / volume;
            }
        }
        let group_fission: Vec<f64> = data.fission.iter().zip(&flux).map(|(s, f)| s * f).collect();
        let fission: f64 = group_fission.iter().sum();
        match &nuclide.fission {
            Some(f) => {
                rates[i].fission = fission * BARN / volume;
                rates[i].yield_weights = yield_weights(group_bounds, &group_fission, f);
            }
            None if fission > 0.0 => {
                bail!("{name} fissions in the library but has no fission in the chain")
            }
            None => {}
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
