//! 影の計算（黒体にした吸収体の Dancoff 係数）に使う、黒体でない物質の全断面積。
//!
//! Dancoff 係数は、影を受ける吸収体 A が吸収するエネルギーでの、黒体でない領域の透過である。そこでポイントワイズの断面積で
//! 解ける群（ポイントワイズの断面積の範囲の下限から [`MAX_TOP`] まで、範囲に収まる群）では、物質 m の全断面積を A の吸収率で
//! 重みづけして平均する:
//!
//! Σ_m,g = ∫_g Σt,m w dE / ∫_g w dE、w = Σa,A φ_A、φ_A = 1/(E (Σt,A + 1/l̄_A))
//!
//! φ_A は NR 近似の中性子束で、l̄_A は孤立した A の平均弦長（Dancoff 係数を入れると、係数とこの重みが循環する
//! ので入れない）。Σt,A は A のすべての核種の全断面積で、A 自身のポテンシャル散乱を含む（NR の中性子束の分子
//! (Σp + 1/l̄) は定数で、重みの比では消える）。Σa = Σt − Σe（弾性散乱以外）で、範囲の上限より下では捕獲と
//! 核分裂に、U-235 の低い準位（77 eV から）などのわずかな非弾性散乱が加わる。
//!
//! 範囲の外の群と、A に分離共鳴を持つ核種がないとき（重みに解くべき共鳴の構造がない）は、無限希釈の全断面積とする
//! （docs/plans/m3-core-hzp.md）。

use std::borrow::Cow;

use xlc_format::library::{Library, LookupError, TemperatureData};

use crate::input::MaterialInput;
use crate::ufg::{MAX_TOP, lethargy_grid};

/// 重みの積分のレサジーの刻みの上限。ポイントワイズの断面積の点の間をこれ以下に分けて台形則で積分する（誤差は刻みの2乗に比例する）
const STEP: f64 = 1.0e-3;

/// 影を受ける吸収体 A（黒体にする物質）。
#[derive(Debug, Clone, PartialEq)]
pub struct Absorber<'a> {
    /// 吸収体の物質と、その体積 (cm²)。複数なら巨視的断面積を体積で平均する
    pub materials: Vec<(&'a MaterialInput, f64)>,
    /// 孤立した吸収体の平均弦長 (cm)
    pub chord: f64,
}

/// 影の計算の全断面積のエラー。
#[derive(Debug, thiserror::Error)]
pub enum ShadowError {
    #[error("material `{material}`: {source}")]
    Lookup {
        material: String,
        source: LookupError,
    },
    #[error("{name}: pointwise data are not defined at {energy} eV")]
    Undefined { name: String, energy: f64 },
    #[error("absorber chord {0} cm must be positive and finite")]
    Chord(f64),
    #[error("absorber volume {0} cm² must be positive and finite")]
    Volume(f64),
    #[error("group {0}: the absorber does not absorb")]
    NoAbsorption(usize),
    #[error("group {0}: non-finite shadow cross section")]
    NonFinite(usize),
    #[error("the absorber has resolved resonances but the library has no pointwise range")]
    NoPointwise,
}

/// 吸収体 `absorber` の影の計算に使う、`materials` の物質ごと・群ごとの全断面積 (1/cm)（モジュールの説明を参照）。
pub fn shadow_totals(
    library: &Library,
    absorber: &Absorber,
    materials: &[&MaterialInput],
) -> Result<Vec<Vec<f64>>, ShadowError> {
    if !(absorber.chord > 0.0 && absorber.chord.is_finite()) {
        return Err(ShadowError::Chord(absorber.chord));
    }
    if let Some(&(_, v)) = absorber
        .materials
        .iter()
        .find(|&&(_, v)| !(v > 0.0 && v.is_finite()))
    {
        return Err(ShadowError::Volume(v));
    }
    let mut tables = Tables::default();
    let mut totals = materials
        .iter()
        .map(|m| dilute(library, &mut tables, m))
        .collect::<Result<Vec<_>, _>>()?;
    // 吸収体の核種（表の番号、体積で平均した個数密度）と、物質ごとの核種（表の番号、個数密度）。個数密度が 0 の核種は
    // 何にも寄与しないので除く
    let volume: f64 = absorber.materials.iter().map(|&(_, v)| v).sum();
    let mut absorbing = Vec::new();
    for &(m, v) in &absorber.materials {
        for (name, density) in m.composition.iter().filter(|(_, n)| *n > 0.0) {
            absorbing.push((tables.index(library, m, name)?, density * v / volume));
        }
    }
    let targets = materials
        .iter()
        .map(|m| {
            m.composition
                .iter()
                .filter(|(_, n)| *n > 0.0)
                .map(|(name, density)| Ok((tables.index(library, m, name)?, *density)))
                .collect::<Result<Vec<_>, ShadowError>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    if !absorbing.iter().any(|&(k, _)| tables.resolved[k]) {
        return Ok(totals);
    }
    let range = library.pointwise_range.ok_or(ShadowError::NoPointwise)?;
    let (low, high) = (range[0], range[1].min(MAX_TOP));
    let mut involved: Vec<usize> = absorbing
        .iter()
        .chain(targets.iter().flatten())
        .map(|&(k, _)| k)
        .collect();
    involved.sort_unstable();
    involved.dedup();
    let escape = 1.0 / absorber.chord;
    for g in 0..library.groups.len() {
        let (lower, upper) = library.groups.range(g);
        if !(lower >= low && upper <= high) {
            continue;
        }
        // 格子（降順）: 群の境界と、核種のポイントワイズの断面積の点。点の間は STEP 以下に分ける
        let mut energies: Vec<f64> = involved
            .iter()
            .flat_map(|&k| tables.tables[k].pointwise.energy.iter().copied())
            .filter(|&e| e > lower && e < upper)
            .chain([upper, lower])
            .collect();
        energies.sort_by(|a, b| b.total_cmp(a));
        energies.dedup();
        let (u, points) = lethargy_grid(&energies, upper, STEP, |e| e == upper || e == lower);
        // 点 e での [核種] の [全, 弾性散乱]。不連続点では、`below` なら低い側の極限
        let sample = |e: f64, below: bool| -> Result<Vec<[f64; 2]>, ShadowError> {
            let mut values = vec![[0.0; 2]; tables.tables.len()];
            for &k in &involved {
                values[k] = tables.at(k, e, below)?;
            }
            Ok(values)
        };
        // 点での重み Σa,A/(Σt,A + 1/l̄) と、物質ごとの Σt
        let weigh = |values: &[[f64; 2]]| -> (f64, Vec<f64>) {
            let (mut total, mut absorption) = (0.0, 0.0);
            for &(k, n) in &absorbing {
                total += n * values[k][0];
                absorption += n * (values[k][0] - values[k][1]);
            }
            let sigma = targets
                .iter()
                .map(|t| t.iter().map(|&(k, n)| n * values[k][0]).sum())
                .collect();
            (absorption / (total + escape), sigma)
        };
        let mut numerator = vec![0.0; materials.len()];
        let mut denominator = 0.0;
        for i in 0..points.len().saturating_sub(1) {
            let du = u[i + 1] - u[i];
            // 区間は points[i] から下。上端は低い側、下端は高い側の極限
            let (wa, sa) = weigh(&sample(points[i], true)?);
            let (wb, sb) = weigh(&sample(points[i + 1], false)?);
            denominator += 0.5 * du * (wa + wb);
            for (n, (a, b)) in numerator.iter_mut().zip(sa.iter().zip(&sb)) {
                *n += 0.5 * du * (a * wa + b * wb);
            }
        }
        if !(denominator > 0.0 && denominator.is_finite()) {
            return Err(ShadowError::NoAbsorption(g));
        }
        for (total, n) in totals.iter_mut().zip(&numerator) {
            let value = n / denominator;
            if !value.is_finite() {
                return Err(ShadowError::NonFinite(g));
            }
            total[g] = value;
        }
    }
    Ok(totals)
}

/// 核種（名前と温度の組）ごとの表。同じ核種を何度も引かない
#[derive(Default)]
struct Tables<'a> {
    keys: Vec<(&'a str, f64)>,
    tables: Vec<Cow<'a, TemperatureData>>,
    /// 分離共鳴を持つか
    resolved: Vec<bool>,
}

impl<'a> Tables<'a> {
    /// 物質 `material` の核種 `name` の表の番号
    fn index(
        &mut self,
        library: &'a Library,
        material: &MaterialInput,
        name: &'a str,
    ) -> Result<usize, ShadowError> {
        let t = material.temperature;
        if let Some(k) = self.keys.iter().position(|&(n, tk)| n == name && tk == t) {
            return Ok(k);
        }
        let (data, table) = library
            .lookup(name, t)
            .map_err(|source| ShadowError::Lookup {
                material: material.name.clone(),
                source,
            })?;
        self.keys.push((name, t));
        self.tables.push(table);
        self.resolved.push(data.resolved.is_some());
        Ok(self.tables.len() - 1)
    }

    /// 表 `k` のポイントワイズの断面積 [全, 弾性散乱] の、エネルギー `e` での lin-lin の補間。不連続点（同じエネルギーの2点。
    /// 低い側の極限が先）では、`below` なら低い側、そうでなければ高い側の極限
    fn at(&self, k: usize, e: f64, below: bool) -> Result<[f64; 2], ShadowError> {
        let p = &self.tables[k].pointwise;
        let x = &p.energy;
        let undefined = || ShadowError::Undefined {
            name: self.keys[k].0.to_string(),
            energy: e,
        };
        let (Some(&first), Some(&last)) = (x.first(), x.last()) else {
            return Err(undefined());
        };
        if !(x.len() >= 2 && first <= e && e <= last) {
            return Err(undefined());
        }
        // j は e 以上の最初の点（`below`）または e 以下の最後の点
        let j = if below {
            x.partition_point(|&v| v < e)
        } else {
            x.partition_point(|&v| v <= e) - 1
        };
        if x[j] == e {
            return Ok([p.total[j], p.elastic[j]]);
        }
        // e は x[a] と x[a + 1] の間
        let a = if below { j - 1 } else { j };
        let t = (e - x[a]) / (x[a + 1] - x[a]);
        let lerp = |y: &[f64]| y[a] + t * (y[a + 1] - y[a]);
        Ok([lerp(&p.total), lerp(&p.elastic)])
    }
}

/// 物質の無限希釈の全断面積 (1/cm)
fn dilute<'a>(
    library: &'a Library,
    tables: &mut Tables<'a>,
    material: &'a MaterialInput,
) -> Result<Vec<f64>, ShadowError> {
    let mut total = vec![0.0; library.groups.len()];
    for (name, density) in &material.composition {
        let k = tables.index(library, material, name)?;
        for (sum, sigma) in total.iter_mut().zip(&tables.tables[k].shielded[0].total) {
            *sum += density * sigma;
        }
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::Material;
    use xlc_format::group::GroupStructure;
    use xlc_format::library::{NuclideData, Pointwise, Resolved, Shielded};

    const T: f64 = 600.0;
    /// 群 1（100 eV〜10 keV）と群 2（4〜100 eV）がポイントワイズの断面積の範囲
    const BOUNDS: [f64; 6] = [1.0e7, 1.0e4, 1.0e2, 4.0, 0.625, 1.0e-5];

    /// 30 eV または 60 eV の、半幅 0.5 eV、高さ 1e4 barn の三角形の共鳴（裾がなく、離れた共鳴とは重ならない）
    fn resonance(e: f64, at: f64) -> f64 {
        1.0e4 * (1.0 - (e - at).abs() / 0.5).max(0.0)
    }

    /// ポイントワイズの断面積 [エネルギー, 全, 弾性散乱] の核種。群定数の全断面積は `dilute`
    fn nuclide(name: &str, points: Vec<[f64; 3]>, dilute: f64, resolved: bool) -> NuclideData {
        let ng = BOUNDS.len() - 1;
        let column = |k: usize| points.iter().map(|p| p[k]).collect::<Vec<_>>();
        let table = Shielded {
            total: vec![dilute; ng],
            elastic: vec![dilute; ng],
            fission: vec![0.0; ng],
            nu_fission: vec![0.0; ng],
        };
        NuclideData {
            name: name.into(),
            za: 0.0,
            awr: 50.0,
            resolved: resolved.then_some(Resolved {
                top: 1.0e4,
                potential: 5.0,
            }),
            thermal_scattering: false,
            temperatures: vec![TemperatureData {
                temperature: T,
                chi: vec![0.0; ng],
                scatter: vec![vec![vec![0.0; ng]; ng]],
                shielded: vec![table.clone(), table],
                pointwise: Pointwise {
                    energy: column(0),
                    total: column(1),
                    elastic: column(2),
                    fission: vec![0.0; points.len()],
                    nu_fission: vec![0.0; points.len()],
                },
            }],
        }
    }

    /// 共鳴を細かく表す点の並び（昇順）
    fn energies() -> Vec<f64> {
        let mut e = vec![4.0, 10.0];
        for at in [30.0, 60.0] {
            e.extend((0..=400).map(|k| at - 2.0 + 0.01 * k as f64));
        }
        e.extend([100.0, 1.0e3]);
        e.extend((0..=400).map(|k| 2998.0 + 0.01 * k as f64));
        e.push(1.0e4);
        e
    }

    fn library() -> Library {
        let along = |f: &dyn Fn(f64) -> [f64; 2]| -> Vec<[f64; 3]> {
            energies()
                .into_iter()
                .map(|e| {
                    let [t, s] = f(e);
                    [e, t, s]
                })
                .collect()
        };
        // 30 eV と 3 keV で捕獲する吸収体
        let absorber = along(&|e| [10.0 + resonance(e, 30.0) + resonance(e, 3.0e3), 10.0]);
        // E について線形な全断面積（4 eV で 2、10 keV で 12）
        let linear = |e: f64| 2.0 + 10.0 * (e - 4.0) / (1.0e4 - 4.0);
        Library {
            pointwise_range: Some([4.0, 1.0e4]),
            groups: GroupStructure::new(BOUNDS.to_vec()).unwrap(),
            sigma0: vec![1.0e10, 1.0],
            nuclides: vec![
                nuclide("A", absorber, 12.0, true),
                // 吸収が一定の吸収体（分離共鳴を持つとする）
                nuclide("uniform", along(&|_| [10.0, 5.0]), 10.0, true),
                // 分離共鳴を持たない吸収体
                nuclide("plain", along(&|_| [10.0, 5.0]), 10.0, false),
                nuclide("flat", along(&|_| [6.5, 6.5]), 7.0, false),
                nuclide("linear", along(&|e| [linear(e), linear(e)]), 7.0, false),
                nuclide(
                    "apart",
                    along(&|e| [5.0 + resonance(e, 60.0), 5.0]),
                    9.0,
                    true,
                ),
                nuclide(
                    "overlap",
                    along(&|e| [5.0 + resonance(e, 30.0), 5.0]),
                    9.0,
                    true,
                ),
            ],
        }
    }

    fn material(name: &str, nuclide: &str, density: f64) -> MaterialInput {
        MaterialInput {
            name: name.into(),
            temperature: T,
            composition: vec![(nuclide.into(), density)],
        }
    }

    fn totals(absorber: &str, target: &str) -> Vec<f64> {
        let library = library();
        let a = material("absorber", absorber, 0.02);
        let m = material("target", target, 0.04);
        let absorber = Absorber {
            materials: vec![(&a, 1.0)],
            chord: 1.0,
        };
        let mut result = shadow_totals(&library, &absorber, &[&m]).unwrap();
        result.remove(0)
    }

    fn close(a: f64, b: f64, tolerance: f64) -> bool {
        (a - b).abs() <= tolerance * b.abs()
    }

    #[test]
    fn flat_cross_sections_are_kept_inside_the_range() {
        // ポイントワイズの断面積が一定なら、重みによらずその値。範囲の外の群は無限希釈の値
        let t = totals("A", "flat");
        let expected = [7.0, 6.5, 6.5, 7.0, 7.0].map(|s| 0.04 * s);
        for (a, b) in t.iter().zip(expected) {
            assert!(close(*a, b, 1e-12), "{t:?}");
        }
    }

    #[test]
    fn uniform_absorption_gives_the_lethargy_average() {
        // 吸収体の断面積が一定なら重みはレサジーで一様で、E について線形な σ の平均は
        // 2 + 10/(1e4 − 4) ((E₂ − E₁)/ln(E₂/E₁) − 4)
        let t = totals("uniform", "linear");
        for g in [1, 2] {
            let (lo, hi) = (BOUNDS[g + 1], BOUNDS[g]);
            let average = 2.0 + 10.0 / (1.0e4 - 4.0) * ((hi - lo) / (hi / lo).ln() - 4.0);
            assert!(close(t[g], 0.04 * average, 1e-6), "group {g}: {t:?}");
        }
    }

    #[test]
    fn resonances_apart_from_the_absorber_are_not_seen() {
        // 吸収体が 30 eV でだけ吸収するとき、60 eV の共鳴は見えず、共鳴の外の 5 barn
        let t = totals("A", "apart");
        assert!(close(t[2], 0.04 * 5.0, 1e-9), "{t:?}");
    }

    #[test]
    fn overlapping_resonances_shadow() {
        // 同じ 30 eV の共鳴は、吸収体が吸収するエネルギーで不透明
        let t = totals("A", "overlap");
        assert!(t[2] > 2.0 * 0.04 * 5.0, "{t:?}");
    }

    #[test]
    fn absorbers_without_resolved_resonances_keep_infinite_dilution() {
        let library = library();
        let t = totals("plain", "apart");
        let m = material("target", "apart", 0.04);
        assert_eq!(t, Material::mix(&library, &m.composition, T).unwrap().total);
    }

    #[test]
    fn rejects_non_positive_chords() {
        let library = library();
        let a = material("absorber", "A", 0.02);
        let absorber = Absorber {
            materials: vec![(&a, 1.0)],
            chord: 0.0,
        };
        let result = shadow_totals(&library, &absorber, &[&a]);
        assert!(matches!(result, Err(ShadowError::Chord(_))), "{result:?}");
    }
}
