//! 燃料の自己遮蔽。背景断面積は ESSM（埋め込み自己遮蔽法）で求める（docs/plans/m2-assembly.md）。
//!
//! 共鳴の群ごとの固定源の輸送方程式 Ω·∇ψ + Σ_i (Σa,i + Σp,i) ψ = Σ_i Σp,i / 4π（NR 近似）を格子全体の MOC で
//! 解き、ピンの燃料の平均の中性子束 φ から Σb = Σa φ / (1 − φ) とする。ライブラリの σ0 は吸収体自身の散乱を
//! 含まないので、核種 R の σ0 は (Σb − N_R σp,R + Σ_{j≠R} N_j σa,j) / N_R（均質な媒質では厳密な σ0 に戻る）。
//! 吸収は σ0 によるので、σ0 と実効断面積を収束するまで交互に求める。
//!
//! 等価原理（Dancoff 係数と Carlvik の2項有理近似）の関数は、置き換える前の方式として残っている。

use xlc_xs::library::{Library, LookupError, Nuclide, Shielded};

use crate::input::{MaterialInput, VOID};
use crate::material::{Material, MaterialError};
use crate::moc::fixed_source::{FixedSourceError, FixedSourceOptions, fixed_source};
use crate::tracking::Tracking;

/// 格子の燃料の自己遮蔽の反復の条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShieldingOptions {
    /// 燃料の各核種の全断面積の相対変化の上限
    pub tolerance: f64,
    pub max_iterations: usize,
    pub fixed_source: FixedSourceOptions,
}

impl Default for ShieldingOptions {
    fn default() -> Self {
        Self {
            tolerance: 1.0e-6,
            max_iterations: 100,
            fixed_source: FixedSourceOptions::default(),
        }
    }
}

/// 格子の燃料の自己遮蔽のエラー。
#[derive(Debug, thiserror::Error)]
pub enum SelfShieldingError {
    #[error("material `{0}` is not defined")]
    UnknownMaterial(String),
    #[error("{given} region materials for {regions} regions")]
    Regions { given: usize, regions: usize },
    #[error("cell {0}: more than one fuel material")]
    MixedFuel(usize),
    #[error("cell {0}: non-finite shielded cross sections")]
    NonFinite(usize),
    #[error("shielding not converged in {0} iterations")]
    NotConverged(usize),
    #[error("material `{name}`: {source}")]
    Material { name: String, source: MaterialError },
    #[error(transparent)]
    FixedSource(#[from] FixedSourceError),
    #[error(transparent)]
    Shielding(#[from] ShieldingError),
}

impl From<LookupError> for SelfShieldingError {
    fn from(error: LookupError) -> Self {
        Self::Shielding(ShieldingError::Lookup(error))
    }
}

/// 自己遮蔽した物質。燃料はセルごとに、それ以外は物質ごとに1つずつ持つ。
#[derive(Debug, Clone, PartialEq)]
pub struct FuelShielding {
    pub materials: Vec<Material>,
    /// [領域] の `materials` の番号
    pub region_material: Vec<usize>,
    /// 燃料のあるセルごとの自己遮蔽の結果（セルの順）
    pub pins: Vec<PinShielding>,
}

/// 1本のピンの燃料の自己遮蔽の結果。
#[derive(Debug, Clone, PartialEq)]
pub struct PinShielding {
    pub cell: usize,
    /// [核種][群] の背景断面積 σ0 (barn)（`nuclides` の順）
    pub sigma0: Vec<Vec<f64>>,
    /// 燃料の核種の自己遮蔽した群定数 (barn) と個数密度 (1/(barn cm))
    pub nuclides: Vec<(Nuclide, f64)>,
    /// 核種ごとの実効断面積（`nuclides` の順。弾性散乱を含む）
    pub effective: Vec<Shielded>,
}

/// 核種の、ESSM の方程式に使う吸収 σa（全 − 弾性）とポテンシャル散乱 σp (barn) の群ごとの値。
struct Parts {
    absorption: Vec<f64>,
    potential: Vec<f64>,
}

/// 1本のピンの燃料の反復の状態。
struct Pin {
    cell: usize,
    regions: Vec<usize>,
    material: usize,
    /// 核種ごとの今の実効断面積と σ0
    tables: Vec<Shielded>,
    sigma0: Vec<Vec<f64>>,
}

/// 格子の燃料（核分裂する核種を含む物質）を、ESSM の背景断面積で自己遮蔽する。
/// 被覆管と減速材は無限希釈のまま。`region_materials` は領域ごとの物質の名前。
pub fn shield_fuel(
    library: &Library,
    tracking: &Tracking,
    materials: &[MaterialInput],
    region_materials: &[String],
    options: &ShieldingOptions,
) -> Result<FuelShielding, SelfShieldingError> {
    let lattice = &tracking.lattice;
    let regions = lattice.volumes().len();
    if region_materials.len() != regions {
        return Err(SelfShieldingError::Regions {
            given: region_materials.len(),
            regions,
        });
    }
    let ng = library.groups.len();
    let infinite = *library.sigma0.first().ok_or(LookupError::NoSigma0)?;
    let mut dilute = materials
        .iter()
        .map(|m| {
            Material::mix(library, &m.composition, m.temperature).map_err(|source| {
                SelfShieldingError::Material {
                    name: m.name.clone(),
                    source,
                }
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    // 真空の帯は、物質の一覧の最後に置く真空の物質にする
    let void = dilute.len();
    dilute.push(Material::void(ng));
    let region_index = region_materials
        .iter()
        .map(|name| {
            if name == VOID {
                return Ok(void);
            }
            materials
                .iter()
                .position(|m| &m.name == name)
                .ok_or_else(|| SelfShieldingError::UnknownMaterial(name.clone()))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let is_fuel: Vec<bool> = dilute
        .iter()
        .map(|m| m.nu_fission.iter().any(|&v| v > 0.0))
        .collect();

    // 物質ごと・核種ごとの σa と σp（無限希釈）
    let parts = |m: &MaterialInput| -> Result<Vec<Parts>, SelfShieldingError> {
        m.composition
            .iter()
            .map(|(name, _)| {
                let table = library.shielded(name, m.temperature, &vec![infinite; ng])?;
                let potential = library
                    .nuclides
                    .iter()
                    .find(|n| &n.name == name)
                    .and_then(|n| n.potential);
                Ok(Parts {
                    absorption: absorption(&table),
                    potential: potential.map_or(table.elastic.clone(), |p| vec![p; ng]),
                })
            })
            .collect()
    };
    let material_parts = materials.iter().map(parts).collect::<Result<Vec<_>, _>>()?;
    // 燃料以外の領域の、ESSM の全断面積と線源（反復で変わらない）
    let fixed_terms = |index: usize| -> (Vec<f64>, Vec<f64>) {
        let Some(m) = materials.get(index) else {
            return (vec![0.0; ng], vec![0.0; ng]);
        };
        let mut total = vec![0.0; ng];
        let mut source = vec![0.0; ng];
        for ((_, density), p) in m.composition.iter().zip(&material_parts[index]) {
            for g in 0..ng {
                total[g] += density * (p.absorption[g] + p.potential[g]);
                source[g] += density * p.potential[g];
            }
        }
        (total, source)
    };

    // ピンごとの燃料の領域と物質
    let mut pins = Vec::new();
    for cell in 0..lattice.cells_per_side().pow(2) {
        let fuel_regions: Vec<usize> = lattice
            .cell_regions(cell)
            .filter(|&r| is_fuel[region_index[r]])
            .collect();
        let Some(&first) = fuel_regions.first() else {
            continue;
        };
        let material = region_index[first];
        if fuel_regions.iter().any(|&r| region_index[r] != material) {
            return Err(SelfShieldingError::MixedFuel(cell));
        }
        let m = &materials[material];
        let tables = m
            .composition
            .iter()
            .map(|(name, _)| library.shielded(name, m.temperature, &vec![infinite; ng]))
            .collect::<Result<Vec<_>, _>>()?;
        pins.push(Pin {
            cell,
            regions: fuel_regions,
            material,
            sigma0: vec![vec![infinite; ng]; tables.len()],
            tables,
        });
    }
    let volumes = tracking.volumes(regions);

    for _ in 0..options.max_iterations {
        // 領域ごとの ESSM の全断面積と線源
        let (mut totals, mut sources): (Vec<Vec<f64>>, Vec<Vec<f64>>) =
            region_index.iter().map(|&i| fixed_terms(i)).unzip();
        for pin in &pins {
            let m = &materials[pin.material];
            let (mut total, mut source) = (vec![0.0; ng], vec![0.0; ng]);
            for (((_, density), table), p) in m
                .composition
                .iter()
                .zip(&pin.tables)
                .zip(&material_parts[pin.material])
            {
                let a = absorption(table);
                for g in 0..ng {
                    total[g] += density * (a[g] + p.potential[g]);
                    source[g] += density * p.potential[g];
                }
            }
            for &r in &pin.regions {
                totals[r].clone_from(&total);
                sources[r].clone_from(&source);
            }
        }
        let flux = fixed_source(tracking, &totals, &sources, &options.fixed_source)?;

        let mut change = 0.0_f64;
        for pin in &mut pins {
            let m = &materials[pin.material];
            let parts = &material_parts[pin.material];
            let volume: f64 = pin.regions.iter().map(|&r| volumes[r]).sum();
            let absorptions: Vec<Vec<f64>> = pin.tables.iter().map(absorption).collect();
            let mut new_tables = Vec::with_capacity(pin.tables.len());
            for (k, (name, density)) in m.composition.iter().enumerate() {
                if *density == 0.0 {
                    // 個数密度が 0 の核種は何にも寄与しないので、無限希釈のまま
                    new_tables.push(pin.tables[k].clone());
                    continue;
                }
                let sigma0: Vec<f64> = (0..ng)
                    .map(|g| {
                        let phi = pin
                            .regions
                            .iter()
                            .map(|&r| volumes[r] * flux[r][g])
                            .sum::<f64>()
                            / volume;
                        let sigma_a: f64 = m
                            .composition
                            .iter()
                            .zip(&absorptions)
                            .map(|((_, n), a)| n * a[g])
                            .sum();
                        let others: f64 = m
                            .composition
                            .iter()
                            .zip(&absorptions)
                            .enumerate()
                            .filter(|&(j, _)| j != k)
                            .map(|(_, ((_, n), a))| n * a[g])
                            .sum();
                        // 吸収がない（中性子束が減らない）群は無限希釈
                        if sigma_a <= 0.0 || phi >= 1.0 {
                            return infinite;
                        }
                        let background = sigma_a * phi / (1.0 - phi);
                        let s = (background - density * parts[k].potential[g] + others) / density;
                        // 表の最初の σ0 は無限希釈そのもの。それより大きい σ0 は無限希釈とする
                        s.min(infinite)
                    })
                    .collect();
                new_tables.push(library.shielded(name, m.temperature, &sigma0)?);
                pin.sigma0[k] = sigma0;
            }
            let finite = |s: &Shielded| {
                [&s.total, &s.elastic, &s.fission, &s.nu_fission]
                    .iter()
                    .all(|v| v.iter().all(|x| x.is_finite()))
            };
            if !new_tables.iter().all(finite) {
                return Err(SelfShieldingError::NonFinite(pin.cell));
            }
            for (new, old) in new_tables.iter().zip(&pin.tables) {
                for (a, b) in new.total.iter().zip(&old.total) {
                    change = change.max(((a - b) / b).abs());
                }
            }
            pin.tables = new_tables;
        }
        if change <= options.tolerance {
            return finish(library, materials, dilute, region_index, pins);
        }
    }
    Err(SelfShieldingError::NotConverged(options.max_iterations))
}

/// 吸収（全 − 弾性。共鳴の群では捕獲と核分裂の和）。
fn absorption(table: &Shielded) -> Vec<f64> {
    table
        .total
        .iter()
        .zip(&table.elastic)
        .map(|(t, e)| t - e)
        .collect()
}

/// 収束した実効断面積から、ピンごとの物質と結果を作る。
fn finish(
    library: &Library,
    materials: &[MaterialInput],
    dilute: Vec<Material>,
    region_index: Vec<usize>,
    pins: Vec<Pin>,
) -> Result<FuelShielding, SelfShieldingError> {
    let mut shielded = FuelShielding {
        materials: dilute,
        region_material: region_index,
        pins: Vec::with_capacity(pins.len()),
    };
    for pin in pins {
        let material = &materials[pin.material];
        let nuclides = material
            .composition
            .iter()
            .zip(&pin.tables)
            .map(|((name, density), s)| {
                Ok((
                    shielded_nuclide(library, name, material.temperature, s)?,
                    *density,
                ))
            })
            .collect::<Result<Vec<_>, SelfShieldingError>>()?;
        let mixed = Material::from_nuclides(&nuclides, library.groups.len()).map_err(|source| {
            SelfShieldingError::Material {
                name: material.name.clone(),
                source,
            }
        })?;
        let index = shielded.materials.len();
        shielded.materials.push(mixed);
        for &r in &pin.regions {
            shielded.region_material[r] = index;
        }
        shielded.pins.push(PinShielding {
            cell: pin.cell,
            sigma0: pin.sigma0,
            nuclides,
            effective: pin.tables,
        });
    }
    Ok(shielded)
}

/// 無限希釈の群定数を、自己遮蔽した値に置き換える。散乱行列は無限希釈の形のまま、行の和を
/// 弾性散乱の変化の分だけ変える（M1 で決めたとおり）。
fn shielded_nuclide(
    library: &Library,
    name: &str,
    temperature: f64,
    shielded: &Shielded,
) -> Result<Nuclide, SelfShieldingError> {
    let mut nuclide = library.nuclide(name, temperature)?;
    let ng = library.groups.len();
    let infinite = library.sigma0.first().ok_or(LookupError::NoSigma0)?;
    let dilute = library.shielded(name, temperature, &vec![*infinite; ng])?;
    for (g, row) in nuclide.scatter.iter_mut().enumerate() {
        let sum: f64 = row.iter().sum();
        if sum > 0.0 {
            let factor = (sum + shielded.elastic[g] - dilute.elastic[g]) / sum;
            row.iter_mut().for_each(|v| *v *= factor);
        }
    }
    nuclide.total.clone_from(&shielded.total);
    nuclide.fission.clone_from(&shielded.fission);
    nuclide.nu_fission.clone_from(&shielded.nu_fission);
    Ok(nuclide)
}

/// 自己遮蔽の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum ShieldingError {
    #[error("Dancoff factor {0} must be from 0 to 1")]
    Dancoff(f64),
    #[error("mean chord length {0} cm must be positive")]
    Chord(f64),
    #[error("number density {0} must be positive")]
    Density(f64),
    #[error("{given} rational terms for {groups} groups")]
    Groups { given: usize, groups: usize },
    #[error(transparent)]
    Lookup(#[from] LookupError),
}

/// 有理近似の1つの項。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RationalTerm {
    pub alpha: f64,
    pub beta: f64,
}

/// 1本のピンの燃料の等価原理のパラメータ: ペレットの平均弦長 l̄ (cm) と、群ごとの有理近似の項。
#[derive(Debug, Clone, PartialEq)]
pub struct Equivalence {
    pub chord: f64,
    pub terms: Vec<[RationalTerm; 2]>,
}

/// Dancoff 係数 `dancoff` での Carlvik の2項有理近似の係数。
///
/// A = (1 − C)/C として α₁,₂ = [(5A + 6) ∓ √(A² + 36A + 36)] / (2(A + 1))、
/// β₁ = ((4A + 6)/(A + 1) − α₁)/(α₂ − α₁)、β₂ = 1 − β₁。分子と分母に C を掛けた形
/// α₁,₂ = [(5 + C) ∓ √(1 + 34C + C²)] / 2、(4A + 6)/(A + 1) = 4 + 2C で計算し、C = 0（孤立した棒）でも使えるようにする。
pub fn carlvik(dancoff: f64) -> Result<[RationalTerm; 2], ShieldingError> {
    if !(0.0..=1.0).contains(&dancoff) {
        return Err(ShieldingError::Dancoff(dancoff));
    }
    let c = dancoff;
    let root = (1.0 + 34.0 * c + c * c).sqrt();
    let (alpha1, alpha2) = (0.5 * (5.0 + c - root), 0.5 * (5.0 + c + root));
    let beta1 = (4.0 + 2.0 * c - alpha1) / (alpha2 - alpha1);
    Ok([
        RationalTerm {
            alpha: alpha1,
            beta: beta1,
        },
        RationalTerm {
            alpha: alpha2,
            beta: 1.0 - beta1,
        },
    ])
}

/// 燃料の核種 `name`（個数密度 `density`、1/(barn cm)）の実効断面積。`background` は群ごとの、
/// 同じ燃料の他の核種の全断面積の和 (1/cm)。
pub fn effective(
    library: &Library,
    name: &str,
    temperature: f64,
    density: f64,
    background: &[f64],
    equivalence: &Equivalence,
) -> Result<Shielded, ShieldingError> {
    if !(equivalence.chord > 0.0 && equivalence.chord.is_finite()) {
        return Err(ShieldingError::Chord(equivalence.chord));
    }
    if !(density > 0.0 && density.is_finite()) {
        return Err(ShieldingError::Density(density));
    }
    if equivalence.terms.len() != background.len() {
        return Err(ShieldingError::Groups {
            given: equivalence.terms.len(),
            groups: background.len(),
        });
    }
    let sigma0 = |i: usize| -> Vec<f64> {
        background
            .iter()
            .zip(&equivalence.terms)
            .map(|(b, terms)| (b + terms[i].alpha / equivalence.chord) / density)
            .collect()
    };
    let (sigma0_1, sigma0_2) = (sigma0(0), sigma0(1));
    let tables = [
        library.shielded(name, temperature, &sigma0_1)?,
        library.shielded(name, temperature, &sigma0_2)?,
    ];
    let sigma0 = [sigma0_1, sigma0_2];
    let groups = background.len();
    // 群ごとの重み β_i φ_i
    let weights: Vec<[f64; 2]> = (0..groups)
        .map(|g| {
            std::array::from_fn(|i| {
                let s = sigma0[i][g];
                equivalence.terms[g][i].beta * s / (tables[i].total[g] + s)
            })
        })
        .collect();
    let combine = |of: fn(&Shielded) -> &Vec<f64>| -> Vec<f64> {
        (0..groups)
            .map(|g| {
                let [w1, w2] = weights[g];
                (w1 * of(&tables[0])[g] + w2 * of(&tables[1])[g]) / (w1 + w2)
            })
            .collect()
    };
    Ok(Shielded {
        total: combine(|x| &x.total),
        elastic: combine(|x| &x.elastic),
        fission: combine(|x| &x.fission),
        nu_fission: combine(|x| &x.nu_fission),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use xlc_xs::library::{Library, NuclideData, Shielded, TemperatureData};
    use xlc_xs::multigroup::GroupStructure;

    const T: f64 = 900.0;

    /// 群 0 だけが σ0 に依存する核種
    fn library() -> Library {
        let shielded = |total: f64, elastic: f64| Shielded {
            total: vec![total, 20.0],
            elastic: vec![elastic, 15.0],
            fission: vec![0.0, 0.0],
            nu_fission: vec![0.0, 0.0],
        };
        Library {
            groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
            sigma0: vec![1.0e10, 100.0, 1.0],
            nuclides: vec![NuclideData {
                name: "U-238".into(),
                za: 92238.0,
                awr: 236.0,
                potential: None,
                temperatures: vec![TemperatureData {
                    temperature: T,
                    chi: vec![0.0, 0.0],
                    scatter: vec![vec![6.0, 2.0], vec![0.5, 15.0]],
                    shielded: vec![shielded(10.0, 8.0), shielded(8.0, 7.0), shielded(5.0, 5.0)],
                }],
            }],
        }
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * b.abs().max(1.0)
    }

    #[test]
    fn isolated_rod_gives_carlvik_coefficients() {
        let [a, b] = carlvik(0.0).unwrap();
        assert!(close(a.alpha, 2.0) && close(a.beta, 2.0), "{a:?}");
        assert!(close(b.alpha, 3.0) && close(b.beta, -1.0), "{b:?}");
    }

    #[test]
    fn black_lattice_has_no_escape() {
        // C = 1 では燃料から出た中性子はすべて隣の燃料に入るので、脱出の項がなくなる
        let [a, b] = carlvik(1.0).unwrap();
        assert!(close(a.alpha, 0.0) && close(a.beta, 1.0), "{a:?}");
        assert!(close(b.alpha, 6.0) && close(b.beta, 0.0), "{b:?}");
    }

    #[test]
    fn coefficients_are_ordered_and_betas_sum_to_one() {
        for c in [0.1, 0.5, 0.9] {
            let [a, b] = carlvik(c).unwrap();
            assert!(a.alpha < b.alpha, "{c}");
            assert!(close(a.beta + b.beta, 1.0), "{c}");
        }
    }

    #[test]
    fn rejects_dancoff_factors_outside_zero_to_one() {
        for c in [-0.1, 1.1, f64::NAN] {
            assert!(matches!(carlvik(c), Err(ShieldingError::Dancoff(_))), "{c}");
        }
    }

    #[test]
    fn rejects_non_positive_chord_and_density() {
        let library = library();
        let equivalence = |chord| Equivalence {
            chord,
            terms: vec![carlvik(0.0).unwrap(); 2],
        };
        let shield = |density, chord| {
            effective(
                &library,
                "U-238",
                T,
                density,
                &[0.2; 2],
                &equivalence(chord),
            )
        };
        assert!(matches!(shield(0.02, 0.0), Err(ShieldingError::Chord(_))));
        assert!(matches!(shield(0.0, 0.8), Err(ShieldingError::Density(_))));
        assert!(matches!(
            shield(f64::NAN, 0.8),
            Err(ShieldingError::Density(_))
        ));
    }

    #[test]
    fn without_escape_the_background_alone_sets_sigma0() {
        let library = library();
        let equivalence = Equivalence {
            chord: 0.8192,
            terms: vec![carlvik(1.0).unwrap(); 2],
        };
        // σ0 = 0.2 / 0.02 = 10 b
        let result = effective(&library, "U-238", T, 0.02, &[0.2, 0.2], &equivalence).unwrap();
        let expected = library.shielded("U-238", T, &[10.0, 10.0]).unwrap();
        for (a, b) in result.total.iter().zip(&expected.total) {
            assert!(close(*a, *b), "{a} vs {b}");
        }
    }

    #[test]
    fn two_terms_are_weighted_by_the_narrow_resonance_flux() {
        let library = library();
        let (chord, density, background) = (0.8192, 0.02, 0.05);
        let terms = carlvik(0.0).unwrap();
        let equivalence = Equivalence {
            chord,
            terms: vec![terms; 2],
        };
        let result = effective(
            &library,
            "U-238",
            T,
            density,
            &[background; 2],
            &equivalence,
        )
        .unwrap();
        let (mut numerator, mut denominator) = (0.0, 0.0);
        for term in terms {
            let sigma0 = (background + term.alpha / chord) / density;
            let table = library.shielded("U-238", T, &[sigma0; 2]).unwrap();
            let flux = sigma0 / (table.total[0] + sigma0);
            numerator += term.beta * flux * table.total[0];
            denominator += term.beta * flux;
        }
        assert!(close(result.total[0], numerator / denominator));
        // σ0 に依存しない群はそのまま
        assert!(close(result.total[1], 20.0));
    }

    mod lattice {
        use super::super::*;
        use crate::geometry::{Lattice, PinCell};
        use crate::input::MaterialInput;
        use crate::material::Material;
        use crate::tracking::Tracking;
        use xlc_xs::library::{Library, NuclideData, Shielded, TemperatureData};
        use xlc_xs::multigroup::GroupStructure;

        const T: f64 = 600.0;

        /// 共鳴のある U-238（群 0 が σ0 に依存）、核分裂する U-235、水素
        fn library() -> Library {
            let data = |name: &str, tables: Vec<Shielded>, fissile: bool| NuclideData {
                name: name.into(),
                za: 0.0,
                awr: 1.0,
                // 共鳴のある U-238 だけがポテンシャル散乱を持つ（ほかは弾性散乱を使う）
                potential: (name == "U-238").then_some(11.29),
                temperatures: vec![TemperatureData {
                    temperature: T,
                    chi: if fissile {
                        vec![1.0, 0.0]
                    } else {
                        vec![0.0, 0.0]
                    },
                    scatter: vec![vec![6.0, 1.0], vec![0.0, 4.0]],
                    shielded: tables,
                }],
            };
            let table = |total: [f64; 2], elastic: [f64; 2], fission: [f64; 2]| Shielded {
                total: total.to_vec(),
                elastic: elastic.to_vec(),
                fission: fission.to_vec(),
                nu_fission: fission.iter().map(|f| 2.4 * f).collect(),
            };
            let u238 = |t: f64, e: f64| table([t, 10.0], [e, 4.0], [0.0, 0.0]);
            let constant = |total, elastic, fission| vec![table(total, elastic, fission); 3];
            Library {
                groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
                sigma0: vec![1.0e10, 100.0, 1.0],
                nuclides: vec![
                    data(
                        "U-238",
                        vec![u238(30.0, 8.0), u238(20.0, 7.5), u238(9.0, 7.0)],
                        false,
                    ),
                    data(
                        "U-235",
                        constant([10.0, 50.0], [6.0, 4.0], [2.0, 40.0]),
                        true,
                    ),
                    // 減速材は共鳴の群（群 0）でほとんど吸収しない（ESSM は遠くの中性子束を 1 とみなす）
                    data(
                        "H-1",
                        constant([20.0, 30.0], [19.99, 29.5], [0.0, 0.0]),
                        false,
                    ),
                ],
            }
        }

        fn materials() -> Vec<MaterialInput> {
            vec![
                MaterialInput {
                    name: "fuel".into(),
                    temperature: T,
                    composition: vec![("U-238".into(), 0.022), ("U-235".into(), 7.0e-4)],
                },
                MaterialInput {
                    name: "water".into(),
                    temperature: T,
                    composition: vec![("H-1".into(), 0.05)],
                },
            ]
        }

        /// ピン（燃料と水）と水のセルを `layout` に並べた格子と、領域ごとの物質の名前
        fn geometry(layout: Vec<Vec<usize>>) -> (Tracking, Vec<String>) {
            let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
            let water = PinCell::new(1.26, vec![]).unwrap();
            let names = |t: usize| -> Vec<String> {
                if t == 0 {
                    vec!["fuel".into(), "water".into()]
                } else {
                    vec!["water".into()]
                }
            };
            let regions = layout.iter().flatten().flat_map(|&t| names(t)).collect();
            let lattice = Lattice::new(vec![pin, water], layout).unwrap();
            (Tracking::new(&lattice, 16, 0.02).unwrap(), regions)
        }

        fn fuel_total(shielded: &FuelShielding, region: usize) -> f64 {
            shielded.materials[shielded.region_material[region]].total[0]
        }

        #[test]
        fn shields_fuel_and_keeps_moderator_dilute() {
            let library = library();
            let (tracking, regions) = geometry(vec![vec![0]]);
            let shielded = shield_fuel(
                &library,
                &tracking,
                &materials(),
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            let dilute = |name: &str| {
                let m = materials().into_iter().find(|m| m.name == name).unwrap();
                Material::mix(&library, &m.composition, T).unwrap()
            };
            assert!(fuel_total(&shielded, 0) < dilute("fuel").total[0]);
            assert_eq!(
                shielded.materials[shielded.region_material[1]],
                dilute("water")
            );
        }

        #[test]
        fn reports_each_pin_with_its_nuclides() {
            let library = library();
            let (tracking, regions) = geometry(vec![vec![0, 0], vec![0, 1]]);
            let shielded = shield_fuel(
                &library,
                &tracking,
                &materials(),
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            let cells: Vec<usize> = shielded.pins.iter().map(|p| p.cell).collect();
            assert_eq!(cells, vec![0, 1, 2]);
            let pin = &shielded.pins[0];
            // 核種ごと・群ごとの背景断面積
            assert_eq!(pin.sigma0.len(), 2);
            assert!(pin.sigma0.iter().all(|s| s.len() == 2));
            let names: Vec<&str> = pin.nuclides.iter().map(|(n, _)| n.name.as_str()).collect();
            assert_eq!(names, vec!["U-238", "U-235"]);
            // 実効断面積（弾性散乱を含む）も核種ごとに持ち、全断面積は核種の群定数と同じ
            assert_eq!(pin.effective.len(), 2);
            assert_eq!(pin.effective[0].total, pin.nuclides[0].0.total);
            assert!(pin.effective[0].elastic[0] > 0.0);
            // 物質の全断面積は、核種の全断面積と個数密度の積の和
            let material = &shielded.materials[shielded.region_material[0]];
            let sum: f64 = pin
                .nuclides
                .iter()
                .map(|(n, density)| density * n.total[0])
                .sum();
            assert!((material.total[0] - sum).abs() < 1e-12 * sum);
        }

        #[test]
        fn pins_next_to_water_are_shielded_less() {
            let library = library();
            // 領域: セル0 (燃料, 水), セル1 (燃料, 水), セル2 (燃料, 水), セル3 (水)
            let (tracking, regions) = geometry(vec![vec![0, 0], vec![0, 1]]);
            let shielded = shield_fuel(
                &library,
                &tracking,
                &materials(),
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            let (far, beside, below) = (
                fuel_total(&shielded, 0),
                fuel_total(&shielded, 2),
                fuel_total(&shielded, 4),
            );
            // 水のセルの隣のピンは影が少なく（Dancoff 係数が小さく）、遮蔽が弱い
            assert!(beside > far, "{beside} vs {far}");
            assert!(
                (beside - below).abs() < 1e-3 * beside,
                "{beside} vs {below}"
            );
        }

        #[test]
        fn void_gap_is_kept_void() {
            let library = library();
            let pin = PinCell::new(1.26, vec![0.4096, 0.418]).unwrap();
            let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.02).unwrap();
            let regions: Vec<String> = ["fuel", crate::input::VOID, "water"]
                .map(String::from)
                .to_vec();
            let shielded = shield_fuel(
                &library,
                &tracking,
                &materials(),
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            assert!(shielded.materials[shielded.region_material[1]].is_void());
            assert!(fuel_total(&shielded, 0) > 0.0);
        }

        #[test]
        fn nuclides_without_density_stay_dilute() {
            let library = library();
            let (tracking, regions) = geometry(vec![vec![0]]);
            let shield = |materials: &[MaterialInput]| {
                let s = shield_fuel(
                    &library,
                    &tracking,
                    materials,
                    &regions,
                    &ShieldingOptions::default(),
                )
                .unwrap();
                s.materials[s.region_material[0]].clone()
            };
            let mut with_zero = materials();
            with_zero[0].composition.push(("H-1".into(), 0.0));
            assert_eq!(shield(&with_zero), shield(&materials()));
        }

        #[test]
        fn rejects_non_finite_shielded_cross_sections() {
            let mut library = library();
            library.nuclides[0].temperatures[0].shielded[1].total[0] = f64::NAN;
            let (tracking, regions) = geometry(vec![vec![0]]);
            let result = shield_fuel(
                &library,
                &tracking,
                &materials(),
                &regions,
                &ShieldingOptions::default(),
            );
            assert!(
                matches!(result, Err(SelfShieldingError::NonFinite(0))),
                "{result:?}"
            );
        }

        #[test]
        fn fuel_may_fill_any_rings() {
            // 燃料の形の制約はない（ESSM は燃料の領域の中性子束だけを使う）
            let library = library();
            let (tracking, mut regions) = geometry(vec![vec![0]]);
            regions.swap(0, 1);
            let shielded = shield_fuel(
                &library,
                &tracking,
                &materials(),
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            assert_eq!(shielded.pins.len(), 1);
        }

        #[test]
        fn homogeneous_fuel_gets_the_exact_background() {
            // 格子全体が燃料（均質な無限媒質）なら、U-238 の σ0 は Σ_{j≠R} N_j σt,j / N_R
            let library = library();
            let cell = PinCell::new(1.26, vec![]).unwrap();
            let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.05).unwrap();
            let (n238, n235) = (0.022, 0.01);
            let materials = vec![MaterialInput {
                name: "fuel".into(),
                temperature: T,
                composition: vec![("U-238".into(), n238), ("U-235".into(), n235)],
            }];
            let shielded = shield_fuel(
                &library,
                &tracking,
                &materials,
                &["fuel".into()],
                &ShieldingOptions::default(),
            )
            .unwrap();
            let pin = &shielded.pins[0];
            // U-235 の全断面積は σ0 によらない
            let expected: Vec<f64> = [10.0, 50.0].iter().map(|t| n235 * t / n238).collect();
            for (a, b) in pin.sigma0[0].iter().zip(&expected) {
                assert!((a - b).abs() < 1e-8 * b, "{a} vs {b}");
            }
            let table = library.shielded("U-238", T, &expected).unwrap();
            for (a, b) in pin.effective[0].total.iter().zip(&table.total) {
                assert!((a - b).abs() < 1e-8 * b, "{a} vs {b}");
            }
        }
    }
}
