//! 冷却材（ホウ素の軸のホウ酸水）の核種の、格子かノードで均質化した個数密度と少数群の微視的断面積（ADR-0236）。
//!
//! テーブルは B-10 の値を持ち（`xlc_format::table::Boron10`）、core は状態の B-10 の割合がテーブルと違うとき、個数密度の差に掛けて少数群定数に
//! 足す。定義はミクロ燃焼の微視的断面積（ADR-0216 の「微視的断面積の項目と定義」）の、燃焼の領域をホウ酸水の領域にしたもの:
//! N̄ = Σ_r V_r N_r / V、σ_G = Σ_r Σ_{g∈G} V_r μ_r,g φ_r,g N_r s_r,g / (N̄ Φ_G)、Φ_G = Σ_r Σ_{g∈G} V_r φ_r,g（r はノードの領域のうち
//! ホウ酸水のもの、Φ と V はノードのすべての領域の和、φ は MOC の中性子束、μ は SPH 因子、s はその領域の物質を作った核種の群定数）。
//! N̄ が 0 なら N_r の代わりに 1 で重みを付け、ホウ酸水の体積の割合で割る（ADR-0216 の案 Z1）。

use std::collections::HashMap;

use anyhow::{Context, anyhow, bail};
use xlc_format::library::{Library, LookupError, Nuclide, Shielded};
use xlc_format::table::Boron10;

use crate::fewgroup::coarse_groups;
use crate::geometry::Lattice;
use crate::input::MaterialInput;
use crate::moc::Solution;
use crate::shielding::{FuelShielding, shielded_nuclide_with};
use crate::ufg::NuclideTransfer;

/// 状態の自己遮蔽と、それを解いた物質（ホウ素の軸の物質を状態の組成で作り直したときは、作り直す前の参照の物質。ADR-0201）
#[derive(Clone, Copy)]
pub struct Shielding<'a> {
    pub fuel: &'a FuelShielding,
    pub materials: &'a [MaterialInput],
}

/// 解いた状態の、ホウ酸水の核種の量を作るのに要るもの
#[derive(Clone, Copy)]
pub struct Coolant<'a> {
    pub library: &'a Library,
    /// 自己遮蔽した状態なら、その自己遮蔽（なければ、物質はライブラリの無限希釈の値で作ってある）
    pub shielding: Option<Shielding<'a>>,
    /// 状態の物質（ホウ酸水の組成はここから取る）
    pub state: &'a [MaterialInput],
    /// [領域] の物質の名前
    pub region_materials: &'a [String],
    /// ホウ素の軸のホウ酸水の物質の名前
    pub waters: &'a [String],
    pub lattice: &'a Lattice,
    pub solution: &'a Solution,
    /// 少数群の境界（降順）
    pub bounds: &'a [f64],
}

/// 実効断面積、UFG の範囲の移行、物質の名前が同じ領域の組の鍵（同じ配列を指す領域は、核種の群定数を 1 度だけ作る）
type Key = (usize, usize, usize);

impl Coolant<'_> {
    /// 領域 `regions`（格子全体か 1 つのノード）で均質化した、ホウ酸水の核種 `name` の個数密度と微視的断面積
    pub fn micro(&self, name: &str, regions: &[usize]) -> anyhow::Result<Boron10> {
        let library = self.library;
        let solution = self.solution;
        let fine = library.groups.len();
        let group_of = coarse_groups(&library.groups, self.bounds)?;
        let coarse = self.bounds.len() - 1;
        if solution.flux.len() != self.lattice.regions()
            || solution.volumes.len() != self.lattice.regions()
            || self.region_materials.len() != self.lattice.regions()
        {
            bail!(
                "the solution and the region materials do not match the {} regions of the lattice",
                self.lattice.regions()
            );
        }
        // 領域ごとの SPH 因子 μ（物理的な中性子束は MOC の中性子束に μ を掛けたもの）。自己遮蔽しなければ 1
        let sph = match self.shielding {
            Some(s) => Some(
                s.fuel
                    .physical_flux(self.lattice, &vec![vec![1.0; fine]; self.lattice.regions()])?,
            ),
            None => None,
        };
        let mut flux = vec![0.0; coarse];
        let mut volume = 0.0;
        let mut water_volume = 0.0;
        let mut atoms = 0.0;
        // ホウ酸水の領域: (領域, 個数密度, 物質)
        let mut waters = Vec::new();
        for &r in regions {
            let v = solution.volumes[r];
            volume += v;
            for (g, f) in solution.flux[r].iter().enumerate() {
                flux[group_of[g]] += v * f;
            }
            let material = &self.region_materials[r];
            if !self.waters.contains(material) {
                continue;
            }
            let m = self
                .state
                .iter()
                .find(|m| &m.name == material)
                .ok_or_else(|| anyhow!("material `{material}` is not in the state"))?;
            let n = m
                .composition
                .iter()
                .find(|(n, _)| n == name)
                .map_or(0.0, |(_, n)| *n);
            water_volume += v;
            atoms += v * n;
            waters.push((r, n, m));
        }
        if waters.is_empty() {
            bail!("the regions have no borated water of the boron axis");
        }
        if flux.iter().any(|&f| !(f > 0.0 && f.is_finite())) {
            bail!("the few-group fluxes {flux:?} must be positive");
        }
        let mean = atoms / volume;
        // 同じ群定数を作る領域の組ごとに、重みを付けた細かい群の中性子束を足す
        let mut groups: Vec<(Key, Vec<f64>)> = Vec::new();
        let mut members: HashMap<Key, (usize, &MaterialInput)> = HashMap::new();
        for &(r, n, m) in &waters {
            let weight = if mean > 0.0 { n } else { 1.0 };
            if weight == 0.0 {
                continue;
            }
            let (effective, transfers) = self.region_data(r);
            let key = (
                effective.map_or(0, |e| e.as_ptr() as usize),
                transfers.map_or(0, |t| t.as_ptr() as usize),
                m as *const MaterialInput as usize,
            );
            let weighted = (0..fine).map(|g| {
                weight
                    * solution.volumes[r]
                    * solution.flux[r][g]
                    * sph.as_ref().map_or(1.0, |s| s[r][g])
            });
            match groups.iter_mut().find(|(k, _)| *k == key) {
                Some((_, sum)) => {
                    for (s, w) in sum.iter_mut().zip(weighted) {
                        *s += w;
                    }
                }
                None => {
                    groups.push((key, weighted.collect()));
                    members.insert(key, (r, m));
                }
            }
        }
        let mut total = vec![0.0; coarse];
        let mut scatter = vec![vec![0.0; coarse]; coarse];
        for (key, weighted) in &groups {
            let (r, m) = members[key];
            let nuclide = self
                .nuclide(name, r, m)
                .with_context(|| format!("{name} in region {r} of `{}`", m.name))?;
            for (g, &w) in weighted.iter().enumerate() {
                let big = group_of[g];
                total[big] += w * nuclide.total[g];
                for (to, x) in nuclide.scatter[g].iter().enumerate() {
                    scatter[big][group_of[to]] += w * x;
                }
            }
        }
        // 均質化した個数密度 × 少数群の中性子束で割る。個数密度が 0 なら、ホウ酸水の体積の割合（重み 1 の個数密度）
        let density = if mean > 0.0 {
            mean
        } else {
            water_volume / volume
        };
        Ok(Boron10 {
            density: mean,
            total: total
                .iter()
                .zip(&flux)
                .map(|(t, f)| t / (density * f))
                .collect(),
            scatter: scatter
                .iter()
                .zip(&flux)
                .map(|(row, f)| row.iter().map(|x| x / (density * f)).collect())
                .collect(),
        })
    }

    /// 領域 `r` の自己遮蔽の実効断面積と UFG の範囲の移行（自己遮蔽しなければ、どちらもない）
    fn region_data(&self, r: usize) -> (Option<&[Shielded]>, Option<&[NuclideTransfer]>) {
        match self.shielding {
            Some(s) => (s.fuel.region_cross_sections(r), s.fuel.region_transfers(r)),
            None => (None, None),
        }
    }

    /// 領域 `r`（状態の物質 `m`）の核種 `name` の群定数。物質を作ったものと同じく、自己遮蔽したなら、自己遮蔽を解いた物質の組成の中の
    /// 位置の実効断面積と移行から作り、その組成にない核種と自己遮蔽しない状態は無限希釈の値（`branches::remix_varied`、`Material::mix`）
    fn nuclide(&self, name: &str, r: usize, m: &MaterialInput) -> anyhow::Result<Nuclide> {
        let library = self.library;
        let Some(s) = self.shielding else {
            return Ok(library.nuclide(name, m.temperature)?);
        };
        let solved = s
            .materials
            .iter()
            .find(|x| x.name == m.name)
            .ok_or_else(|| anyhow!("material `{}` is not in the self-shielding", m.name))?;
        if solved.temperature != m.temperature {
            bail!(
                "material `{}` is at {} K in the state and {} K in the self-shielding",
                m.name,
                m.temperature,
                solved.temperature
            );
        }
        let index = solved.composition.iter().position(|(n, _)| n == name);
        let (effective, transfers) = self.region_data(r);
        let shielded = match (index, effective) {
            (Some(i), Some(e)) => {
                if e.len() != solved.composition.len() {
                    bail!(
                        "region {r}: {} shielded nuclides for the {} nuclides of `{}`",
                        e.len(),
                        solved.composition.len(),
                        m.name
                    );
                }
                e[i].clone()
            }
            (Some(_), None) => bail!("region {r} has no shielded cross sections"),
            (None, _) => {
                let infinite = *library.sigma0.first().ok_or(LookupError::NoSigma0)?;
                library.shielded(name, m.temperature, &vec![infinite; library.groups.len()])?
            }
        };
        // 自己遮蔽の組成にない核種は無限希釈で、UFG の範囲の移行も持たない
        let transfer = index.and_then(|i| transfers.map(|t| &t[i]));
        Ok(shielded_nuclide_with(
            library,
            name,
            m.temperature,
            &shielded,
            transfer,
        )?)
    }
}
