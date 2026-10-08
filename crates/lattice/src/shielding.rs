//! 燃料のピンの自己遮蔽。
//!
//! 燃料の分離共鳴の範囲（[`crate::ufg::cell_range`]）は、ピンごとに Dancoff 係数を合わせた円筒の等価なセル
//! （[`crate::equivalent`]）の超微細群の減速計算（[`crate::ufg`]）で、燃料のリング、被覆管、減速材の群定数を
//! それぞれのリングの中性子束で縮約して求める。その範囲の外の群は、燃料を等価原理（Dancoff 係数を入れた
//! Carlvik の2項有理近似）で遮蔽し、被覆管と減速材は無限希釈とする。
//!
//! 燃料のピンの等価なセルに入らない構造材（支持格子の枠、案内管など）のうち分離共鳴を持つ核種を含むものは、
//! 塊ごとに Dancoff 係数と平均弦長を求め、燃料と同じ等価原理で遮蔽する（[`shield_lumps`]）。
//!
//! 等価原理では、燃料の脱出確率を P(x) = Σ_i β_i x / (x + α_i)（x = Σt l̄）で近似し、核種 k の背景断面積を
//! σ0,i = (Σ_{j≠k} N_j σt,j + α_i / l̄) / N_k とする。実効断面積は、2つの σ0 でのライブラリの値を、
//! NR 近似の群平均の中性子束 φ_i = σ0,i / (σt(σ0,i) + σ0,i) と β_i で重みづけして合わせる
//! （ライブラリの重み w/(σt + σ0) での群平均では、この式が厳密に成り立つ）。

use xlc_format::library::{Library, LookupError, Nuclide, Shielded};

use crate::cp::{CpError, Cylinder};
use crate::equivalent::{EquivalentError, MIN_DANCOFF, equivalent_radius};
use crate::geometry::{Lattice, PinCell};
use crate::input::{MaterialInput, VOID};
use crate::material::{Material, MaterialError};
use crate::moc::dancoff::{DancoffError, DancoffOptions, dancoff, lump_dancoff};
use crate::sph::{SphError, sph_factors};
use crate::symmetry::{SymmetryError, cell_orbits};
use crate::tracking::Tracking;
use crate::ufg::{UfgError, cell_range, slowing_down};

/// 格子の燃料の自己遮蔽の反復の条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShieldingOptions {
    /// 燃料の各核種の全断面積の相対変化の上限
    pub tolerance: f64,
    pub max_iterations: usize,
    pub dancoff: DancoffOptions,
}

impl Default for ShieldingOptions {
    fn default() -> Self {
        Self {
            tolerance: 1.0e-6,
            max_iterations: 100,
            dancoff: DancoffOptions::default(),
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
    #[error("cell {0}: fuel must fill the innermost rings of a ring cell")]
    FuelShape(usize),
    #[error("cell {0}: more than one fuel material")]
    MixedFuel(usize),
    #[error(
        "cell {0}: each ring outside the fuel must be one material, and the outer rings the moderator"
    )]
    CellShape(usize),
    #[error("cell {0}: non-finite shielded cross sections")]
    NonFinite(usize),
    #[error("cell {cell}: shielding not converged in {iterations} iterations")]
    NotConverged { cell: usize, iterations: usize },
    #[error("material `{name}`: {source}")]
    Material { name: String, source: MaterialError },
    #[error(transparent)]
    Dancoff(#[from] DancoffError),
    #[error(transparent)]
    Equivalent(#[from] EquivalentError),
    #[error(transparent)]
    Ufg(#[from] UfgError),
    #[error(transparent)]
    Cp(#[from] CpError),
    #[error(transparent)]
    Sph(#[from] SphError),
    #[error(transparent)]
    Symmetry(#[from] SymmetryError),
    #[error(transparent)]
    Shielding(#[from] ShieldingError),
}

impl From<LookupError> for SelfShieldingError {
    fn from(error: LookupError) -> Self {
        Self::Shielding(ShieldingError::Lookup(error))
    }
}

/// 自己遮蔽した物質。燃料はセルごとに、構造材の塊は対称なセルの組と物質ごとに、それ以外は物質ごとに1つずつ持つ。
#[derive(Debug, Clone, PartialEq)]
pub struct FuelShielding {
    pub materials: Vec<Material>,
    /// [領域] の `materials` の番号
    pub region_material: Vec<usize>,
    /// 燃料のあるセルごとの自己遮蔽の結果（セルの順）
    pub pins: Vec<PinShielding>,
    /// 構造材の塊ごとの自己遮蔽の結果
    pub lumps: Vec<LumpShielding>,
}

/// 構造材の1つの塊（[`shield_lumps`]）の自己遮蔽の結果。
#[derive(Debug, Clone, PartialEq)]
pub struct LumpShielding {
    /// 塊のあるセルの組の代表（最も小さいセル番号）
    pub cell: usize,
    pub material: String,
    /// 塊の領域（組のすべてのセルの分）
    pub regions: Vec<usize>,
    /// 群ごとの Dancoff 係数
    pub dancoff: Vec<f64>,
    /// 平均弦長 (cm)
    pub chord: f64,
    /// 物質の核種ごとの実効断面積（物質の組成の順）
    pub effective: Vec<Shielded>,
}

/// 1本のピンの燃料の自己遮蔽の結果。格子の対称性で移り合うピンは、組の代表（最も小さいセル番号）の結果を
/// 共有する（`cell` 以外は代表の値）。
#[derive(Debug, Clone, PartialEq)]
pub struct PinShielding {
    pub cell: usize,
    /// 群ごとの Dancoff 係数
    pub dancoff: Vec<f64>,
    /// 燃料の核種の自己遮蔽した群定数 (barn) と個数密度 (1/(barn cm))
    pub nuclides: Vec<(Nuclide, f64)>,
    /// 核種ごとの実効断面積（`nuclides` の順。弾性散乱を含む）。UFG の範囲の群は、燃料のリングの値を中性子束と
    /// 体積で平均したもの。`nuclides` とともに SPH 因子を掛けていない物理的な値（参照解と比べる用）で、MOC の
    /// 中性子束と組み合わせて反応率を求めるには、領域の物質（SPH 因子を掛けたもの）を使う
    pub effective: Vec<Shielded>,
    /// UFG の減速計算に使った等価なセルの外半径 (cm)。燃料に分離共鳴を持つ核種がなければ `None`
    pub equivalent: Option<f64>,
    /// 等価なセルのリングごと・UFG の範囲の群ごとの SPH 因子（[`crate::sph`]）。その群の領域の物質に掛けてある
    pub sph: Option<Vec<Vec<f64>>>,
}

/// 格子の燃料（核分裂する核種を含む物質）のピンを自己遮蔽する（モジュールの説明を参照）。燃料のないセルは
/// 無限希釈のまま。`region_materials` は領域ごとの物質の名前。
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
    let dilute = materials
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
    let mut dilute = dilute;
    dilute.push(Material::void(library.groups.len()));
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
    let fuel: Vec<bool> = region_index
        .iter()
        .map(|&i| dilute[i].nu_fission.iter().any(|&v| v > 0.0))
        .collect();
    // 燃料の行は使わない
    let totals: Vec<Vec<f64>> = region_index
        .iter()
        .zip(&fuel)
        .map(|(&i, &f)| {
            if f {
                Vec::new()
            } else {
                dilute[i].total.clone()
            }
        })
        .collect();
    // 燃料のあるセルごとに、燃料の物質と、燃料が埋める内側のリングの外半径。Dancoff 係数を
    // 求める前に確かめる（燃料が外側にあると、Dancoff 係数のほうが先に失敗して原因が分かりにくい）
    let mut pins = Vec::new();
    for cell in 0..lattice.cells() {
        let range = lattice.cell_regions(cell);
        let local: Vec<bool> = range.clone().map(|r| fuel[r]).collect();
        if !local.contains(&true) {
            continue;
        }
        // 燃料はピンセルの内側のリングにだけ置ける
        let Some(pin) = lattice.cell_type(cell).pin() else {
            return Err(SelfShieldingError::FuelShape(cell));
        };
        let n_fuel = local.iter().take_while(|&&f| f).count();
        let rings = n_fuel / pin.sectors();
        // 燃料は内側のリングをちょうど埋め、その外側に燃料がないこと
        if pin.is_grid()
            || n_fuel % pin.sectors() != 0
            || rings == 0
            || rings > pin.radii().len()
            || local[n_fuel..].contains(&true)
        {
            return Err(SelfShieldingError::FuelShape(cell));
        }
        let fuel_regions = range.start..range.start + n_fuel;
        let index = region_index[range.start];
        if region_index[fuel_regions.clone()]
            .iter()
            .any(|&i| i != index)
        {
            return Err(SelfShieldingError::MixedFuel(cell));
        }
        pins.push((cell, fuel_regions, index, rings));
    }
    let factors = dancoff(tracking, &totals, &fuel, &options.dancoff)?;

    let mut shielded = FuelShielding {
        materials: dilute,
        region_material: region_index.clone(),
        pins: Vec::new(),
        lumps: Vec::new(),
    };
    // 格子の対称性で移り合うピンは同じ結果になるので、組の代表（最も小さいセル番号。セルの順で最初に現れる）
    // だけを解き、ほかは代表の物質を共有する
    let orbits = cell_orbits(tracking, &region_index)?;
    let mut solved: Vec<(usize, Solved)> = Vec::new();
    for (cell, fuel_regions, index, rings) in pins {
        let factors = factors[cell]
            .as_ref()
            .ok_or(SelfShieldingError::FuelShape(cell))?;
        if let Some((_, s)) = solved.iter().find(|(c, _)| *c == orbits[cell]) {
            assign(
                lattice,
                cell,
                &s.assigned,
                &fuel_regions,
                &mut shielded.region_material,
            )?;
            shielded.pins.push(PinShielding {
                cell,
                ..s.pin.clone()
            });
            continue;
        }
        let radius = lattice
            .cell_type(cell)
            .pin()
            .ok_or(SelfShieldingError::FuelShape(cell))?
            .radii()[rings - 1];
        let equivalence = Equivalence {
            chord: 2.0 * radius,
            terms: factors
                .iter()
                .map(|&c| carlvik(c))
                .collect::<Result<_, _>>()?,
        };
        let material = &materials[index];
        let carlvik = shield_material(library, material, &equivalence, options, cell)?;
        let pin_cell = PinRings {
            lattice,
            cell,
            rings,
            region_index: &region_index,
            void,
        };
        let ufg = shield_cell(
            library,
            &pin_cell,
            materials,
            &shielded.materials,
            factors,
            &carlvik,
        )?;
        let effective = match &ufg {
            Some(c) => c.effective.clone(),
            None => carlvik,
        };
        let nuclides = material
            .composition
            .iter()
            .zip(&effective)
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
        let pin = PinShielding {
            cell,
            dancoff: factors.clone(),
            nuclides,
            effective,
            equivalent: ufg.as_ref().map(|c| c.radius),
            sph: ufg.as_ref().map(|c| c.sph.clone()),
        };
        let assigned = match ufg {
            // リングごとの物質（真空のリングは真空のまま）
            Some(c) => Assigned::Rings {
                materials: c
                    .rings
                    .into_iter()
                    .map(|m| {
                        m.map(|m| {
                            shielded.materials.push(m);
                            shielded.materials.len() - 1
                        })
                    })
                    .collect(),
                start: c.start,
            },
            None => {
                shielded.materials.push(mixed);
                Assigned::Fuel(shielded.materials.len() - 1)
            }
        };
        assign(
            lattice,
            cell,
            &assigned,
            &fuel_regions,
            &mut shielded.region_material,
        )?;
        shielded.pins.push(pin.clone());
        solved.push((cell, Solved { pin, assigned }));
    }
    shield_lumps(
        library,
        tracking,
        materials,
        &region_index,
        &orbits,
        options,
        &mut shielded,
    )?;
    Ok(shielded)
}

/// 燃料でなく、燃料のピンの等価なセルでも置き換えていない領域のうち、分離共鳴を持つ核種を含む物質の領域
/// （支持格子の枠、案内管など）を塊にまとめ、等価原理で遮蔽する。塊は、セルの中で接する（トラックで隣り合う）
/// 同じ物質の領域のひとまとまりで、対称に移り合うセルの組の対応するものは1つにする。分離共鳴を持つ核種の組が
/// 同じ物質ごとに、その物質の領域（ピンの等価なセルで遮蔽した被覆管なども含む）を同時に黒体にして、塊ごとの
/// Dancoff 係数と平均弦長を求める（[`lump_dancoff`]。共鳴を持つ核種の違う物質は、共鳴が別のエネルギーにあるので
/// 影を作らない）。
fn shield_lumps(
    library: &Library,
    tracking: &Tracking,
    materials: &[MaterialInput],
    region_index: &[usize],
    orbits: &[usize],
    options: &ShieldingOptions,
    shielded: &mut FuelShielding,
) -> Result<(), SelfShieldingError> {
    let lattice = &tracking.lattice;
    let regions = region_index.len();
    let is_fuel = |m: &Material| m.nu_fission.iter().any(|&v| v > 0.0);
    // 物質ごとの、分離共鳴を持つ核種（個数密度が正のもの）の名前の組。燃料と、持たない物質は空
    let resonant = materials
        .iter()
        .zip(&shielded.materials)
        .map(|(m, dilute)| {
            let mut names = Vec::new();
            if is_fuel(dilute) {
                return Ok(names);
            }
            for (name, density) in &m.composition {
                if *density > 0.0 && library.lookup(name, m.temperature)?.0.resolved.is_some() {
                    names.push(name.as_str());
                }
            }
            names.sort_unstable();
            names.dedup();
            Ok(names)
        })
        .collect::<Result<Vec<Vec<&str>>, LookupError>>()?;
    // 真空は物質の一覧の外
    let resonant_of = |index: usize| resonant.get(index).map_or(&[][..], Vec::as_slice);
    let candidate: Vec<bool> = (0..regions)
        .map(|r| {
            shielded.region_material[r] == region_index[r]
                && !resonant_of(region_index[r]).is_empty()
        })
        .collect();
    if !candidate.contains(&true) {
        return Ok(());
    }
    // セルの中でトラックが続けて通る、同じ物質の候補の領域をまとめる。まとまりの代表は最も小さい領域の番号
    let mut parent: Vec<usize> = (0..regions).collect();
    let root = |parent: &mut Vec<usize>, mut r: usize| {
        while parent[r] != r {
            parent[r] = parent[parent[r]];
            r = parent[r];
        }
        r
    };
    for track in &tracking.tracks {
        for pair in track.segments.windows(2) {
            let (a, b) = (pair[0].region, pair[1].region);
            if a != b
                && candidate[a]
                && candidate[b]
                && region_index[a] == region_index[b]
                && lattice.region_cell(a) == lattice.region_cell(b)
            {
                let (ra, rb) = (root(&mut parent, a), root(&mut parent, b));
                parent[ra.max(rb)] = ra.min(rb);
            }
        }
    }
    // 塊は (セルの組の代表, 物質の番号, まとまりの代表のセルの中での番号) ごと。対称な操作はリングごとの物質が
    // 同じセルどうしでしか成り立たない（[`cell_orbits`]）ので、組のセルの対応するまとまりは同じ番号になる
    let mut keys: Vec<(usize, usize, usize)> = Vec::new();
    let mut lump = vec![None; regions];
    let mut members: Vec<Vec<usize>> = Vec::new();
    for region in (0..regions).filter(|&r| candidate[r]) {
        let cell = lattice.region_cell(region);
        let first = lattice.cell_regions(cell).start;
        let key = (
            orbits[cell],
            region_index[region],
            root(&mut parent, region) - first,
        );
        let k = match keys.iter().position(|&k| k == key) {
            Some(k) => k,
            None => {
                keys.push(key);
                members.push(Vec::new());
                keys.len() - 1
            }
        };
        lump[region] = Some(k);
        members[k].push(region);
    }
    // 影の計算の全断面積は、燃料のピンの Dancoff 係数と同じく燃料でない領域は無限希釈とし、燃料はピンの
    // SPH 因子を掛けない実効断面積（物理的な値。SPH 因子は MOC の反応率を合わせるためのもの）とする。
    // 対称に移り合うピンは同じ核種を共有するので、組ごとに1回求める
    let ng = library.groups.len();
    let mut totals: Vec<Vec<f64>> = region_index
        .iter()
        .map(|&i| shielded.materials[i].total.clone())
        .collect();
    let mut fuel_totals: Vec<(usize, Vec<f64>)> = Vec::new();
    for pin in &shielded.pins {
        let orbit = orbits[pin.cell];
        let total = match fuel_totals.iter().find(|(o, _)| *o == orbit) {
            Some((_, total)) => total.clone(),
            None => {
                let fuel = Material::from_nuclides(&pin.nuclides, ng).map_err(|source| {
                    SelfShieldingError::Material {
                        name: format!("fuel in cell {}", pin.cell),
                        source,
                    }
                })?;
                fuel_totals.push((orbit, fuel.total.clone()));
                fuel.total
            }
        };
        for region in lattice.cell_regions(pin.cell) {
            if is_fuel(&shielded.materials[region_index[region]]) {
                totals[region].clone_from(&total);
            }
        }
    }
    // 分離共鳴を持つ核種の組ごとに、その組の物質の領域を黒体にし、その組の塊の結果を求める。塊の番号は組の中での
    // 番号に付け直す
    let mut found = Vec::with_capacity(keys.len());
    let mut sets: Vec<&[&str]> = keys.iter().map(|&(_, i, _)| resonant_of(i)).collect();
    sets.sort_unstable();
    sets.dedup();
    for set in sets {
        let group: Vec<usize> = (0..keys.len())
            .filter(|&k| resonant_of(keys[k].1) == set)
            .collect();
        let black: Vec<bool> = region_index
            .iter()
            .map(|&i| resonant_of(i) == set)
            .collect();
        let local: Vec<Option<usize>> = lump
            .iter()
            .map(|l| l.and_then(|l| group.iter().position(|&k| k == l)))
            .collect();
        let lumps = lump_dancoff(
            tracking,
            &totals,
            &black,
            &local,
            group.len(),
            &options.dancoff,
        )?;
        found.extend(group.into_iter().zip(lumps));
    }
    found.sort_by_key(|&(k, _)| k);
    let mut lump_material = Vec::with_capacity(keys.len());
    for ((&(cell, index, _), (_, found)), regions) in keys.iter().zip(&found).zip(members) {
        let equivalence = Equivalence {
            chord: found.chord,
            terms: found
                .dancoff
                .iter()
                .map(|&c| carlvik(c))
                .collect::<Result<_, _>>()?,
        };
        let material = &materials[index];
        let tables = shield_material(library, material, &equivalence, options, cell)?;
        let nuclides = material
            .composition
            .iter()
            .zip(&tables)
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
        shielded.materials.push(mixed);
        lump_material.push(shielded.materials.len() - 1);
        shielded.lumps.push(LumpShielding {
            cell,
            material: material.name.clone(),
            regions,
            dancoff: found.dancoff.clone(),
            chord: found.chord,
            effective: tables,
        });
    }
    for (region, l) in lump.iter().enumerate() {
        if let Some(l) = l {
            shielded.region_material[region] = lump_material[*l];
        }
    }
    Ok(())
}

/// 組の代表のピンの結果
struct Solved {
    pin: PinShielding,
    assigned: Assigned,
}

/// ピンのセルの領域に割り当てる物質の番号
enum Assigned {
    /// 等価なセルのリングごとの物質（真空は `None`）と、減速材のリングの始まり
    Rings {
        materials: Vec<Option<usize>>,
        start: usize,
    },
    /// 燃料の領域だけ（UFG の範囲がないとき）
    Fuel(usize),
}

/// ピンのセル `cell` の領域に、`assigned` の物質の番号を割り当てる。
fn assign(
    lattice: &Lattice,
    cell: usize,
    assigned: &Assigned,
    fuel_regions: &std::ops::Range<usize>,
    region_material: &mut [usize],
) -> Result<(), SelfShieldingError> {
    match assigned {
        Assigned::Rings { materials, start } => {
            let pin = lattice
                .cell_type(cell)
                .pin()
                .ok_or(SelfShieldingError::FuelShape(cell))?;
            let first = lattice.cell_regions(cell).start;
            for (targets, material) in ring_targets(pin, first, *start).iter().zip(materials) {
                if let Some(m) = material {
                    for &region in targets {
                        region_material[region] = *m;
                    }
                }
            }
        }
        Assigned::Fuel(m) => {
            for region in fuel_regions.clone() {
                region_material[region] = *m;
            }
        }
    }
    Ok(())
}

/// 等価なセルのリングごとの、格子の領域（ピンセル `geometry` の最初の領域が `first`）。リング `start` は減速材で、
/// セルのそれより外のリングをすべて含む（外周の枠は含めず、無限希釈のままにする）。
fn ring_targets(geometry: &PinCell, first: usize, start: usize) -> Vec<Vec<usize>> {
    let (rings, sectors) = (geometry.radii().len(), geometry.sectors());
    let regions = |ring: usize| first + ring * sectors..first + (ring + 1) * sectors;
    (0..start)
        .map(|j| regions(j).collect())
        .chain([(start..=rings).flat_map(regions).collect()])
        .collect()
}

/// 燃料のあるピンのセルと、その領域の物質
struct PinRings<'a> {
    lattice: &'a Lattice,
    cell: usize,
    /// 燃料のリングの数
    rings: usize,
    /// [領域] の物質の番号
    region_index: &'a [usize],
    /// 真空の物質の番号
    void: usize,
}

/// ピンのセルの UFG の減速計算の結果
struct CellShielding {
    /// 等価なセルの外半径 (cm)
    radius: f64,
    /// 等価なセルのリングごとの物質（真空は `None`）。格子の領域は [`ring_targets`]
    rings: Vec<Option<Material>>,
    /// 減速材のリングの始まり
    start: usize,
    /// 燃料の核種の実効断面積。UFG の範囲の群は燃料のリングの値を中性子束と体積で平均し、範囲の外は等価原理の値
    effective: Vec<Shielded>,
    /// [リング][UFG の範囲の群] の SPH 因子
    sph: Vec<Vec<f64>>,
}

/// ピンのセルを、Dancoff 係数を合わせた円筒の等価なセルの UFG の減速計算で自己遮蔽する。燃料、燃料と減速材の間の
/// リング（ギャップ、被覆管）、減速材（セルの外側の物質。等価なセルでは1つのリング）のそれぞれの核種に、そのリングの
/// 中性子束で縮約した群定数を与える。UFG の範囲の外の群は、燃料は等価原理の値（`carlvik`）、それ以外は無限希釈。
/// 燃料に分離共鳴を持つ核種がなければ `None`。
fn shield_cell(
    library: &Library,
    pin: &PinRings,
    materials: &[MaterialInput],
    dilute: &[Material],
    dancoff: &[f64],
    carlvik: &[Shielded],
) -> Result<Option<CellShielding>, SelfShieldingError> {
    let cell = pin.cell;
    let geometry = pin
        .lattice
        .cell_type(cell)
        .pin()
        .ok_or(SelfShieldingError::CellShape(cell))?;
    let (radii, sectors) = (geometry.radii(), geometry.sectors());
    let first = pin.lattice.cell_regions(cell).start;
    let regions = |ring: usize| first + ring * sectors..first + (ring + 1) * sectors;
    // リングごとの物質の番号（最後はいちばん外の円の外）。扇形はどれも同じ物質であること
    let ring_index = (0..=radii.len())
        .map(|r| {
            let index = pin.region_index[regions(r).start];
            if pin.region_index[regions(r)].iter().all(|&i| i == index) {
                Ok(index)
            } else {
                Err(SelfShieldingError::CellShape(cell))
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let moderator = ring_index[radii.len()];
    // 減速材のリングの始まり: これより外はすべて減速材
    let start = (pin.rings..=radii.len())
        .find(|&r| ring_index[r..].iter().all(|&i| i == moderator))
        .unwrap_or(radii.len());
    if moderator == pin.void {
        return Err(SelfShieldingError::CellShape(cell));
    }
    let input = |index: usize| (index != pin.void).then(|| &materials[index]);
    let fuel = &materials[ring_index[0]];
    let inputs: Vec<Option<&MaterialInput>> = (0..pin.rings)
        .map(|_| Some(fuel))
        .chain(ring_index[pin.rings..start].iter().map(|&i| input(i)))
        .chain([input(moderator)])
        .collect();
    let Some(range) = cell_range(library, &inputs)? else {
        return Ok(None);
    };
    let bounds = library.groups.bounds();
    let position = |e: f64| {
        bounds.iter().position(|&b| b == e).ok_or(UfgError::Range {
            range,
            pointwise: library.pointwise_range,
        })
    };
    let (top, bottom) = (position(range[1])?, position(range[0])?);

    // 等価なセルの外半径: UFG の範囲の群ごとに Dancoff 係数を合わせた外径の平均
    let inner: Vec<f64> = [radii[pin.rings - 1]]
        .into_iter()
        .chain(radii[pin.rings..start].iter().copied())
        .collect();
    let mut sum = 0.0;
    for (g, &c) in dancoff.iter().enumerate().take(bottom).skip(top) {
        let totals: Vec<f64> = [0.0]
            .into_iter()
            .chain(
                ring_index[pin.rings..start]
                    .iter()
                    .map(|&i| dilute[i].total[g]),
            )
            .chain([dilute[moderator].total[g]])
            .collect();
        // Dancoff 係数がセルの数値誤差より小さいピンは、孤立したピンと区別できない（C の違いは MIN_DANCOFF 以下）
        sum += equivalent_radius(&inner, &totals, c.max(MIN_DANCOFF))?;
    }
    let radius = sum / (bottom - top) as f64;
    let cylinder = Cylinder::new(radii[..start].iter().copied().chain([radius]).collect())?;
    let spectrum = slowing_down(library, &cylinder, &inputs, range)?;
    let volumes = cylinder.volumes();

    let ng = library.groups.len();
    let infinite = *library.sigma0.first().ok_or(LookupError::NoSigma0)?;
    let mut rings = Vec::with_capacity(inputs.len());
    for (j, input) in inputs.iter().enumerate() {
        let Some(m) = input else {
            rings.push(None);
            continue;
        };
        let nuclides = m
            .composition
            .iter()
            .enumerate()
            .map(|(k, (name, density))| {
                let mut full = if j < pin.rings {
                    carlvik[k].clone()
                } else {
                    library.shielded(name, m.temperature, &vec![infinite; ng])?
                };
                let ufg = &spectrum.nuclides[j][k];
                for g in top..bottom {
                    full.total[g] = ufg.total[g - top];
                    full.elastic[g] = ufg.elastic[g - top];
                    full.fission[g] = ufg.fission[g - top];
                    full.nu_fission[g] = ufg.nu_fission[g - top];
                }
                Ok((
                    shielded_nuclide(library, name, m.temperature, &full)?,
                    *density,
                ))
            })
            .collect::<Result<Vec<_>, SelfShieldingError>>()?;
        let material = Material::from_nuclides(&nuclides, ng).map_err(|source| {
            SelfShieldingError::Material {
                name: m.name.clone(),
                source,
            }
        })?;
        rings.push(Some(material));
    }

    // SPH 因子: 同じセルを多群で解いた中性子束のリング間の形が UFG と一致するよう、リングの物質に掛ける。
    // 範囲の上の群の中性子束は UFG と同じく F = 1（群のレサジー幅）
    let multigroup: Vec<Material> = rings
        .iter()
        .map(|m| m.clone().unwrap_or_else(|| Material::void(ng)))
        .collect();
    let above: Vec<f64> = (0..top).map(|g| (bounds[g] / bounds[g + 1]).ln()).collect();
    let sph = sph_factors(&cylinder, &multigroup, top..bottom, &above, &spectrum.flux)?;
    for (ring, (material, mu)) in rings.iter_mut().zip(&sph).enumerate() {
        if let Some(m) = material {
            let mut all = vec![1.0; ng];
            all[top..bottom].copy_from_slice(mu);
            *m = m
                .scaled(&all)
                .map_err(|source| SelfShieldingError::Material {
                    name: format!("ring {ring} of cell {cell}"),
                    source,
                })?;
        }
    }

    // 燃料の平均: 燃料のリングの中性子束と体積の重み
    let effective = carlvik
        .iter()
        .enumerate()
        .map(|(k, c)| {
            let mut s = c.clone();
            for g in top..bottom {
                let weight = |j: usize| spectrum.flux[j][g - top] * volumes[j];
                let total_weight: f64 = (0..pin.rings).map(weight).sum();
                let average = |x: fn(&Shielded) -> &Vec<f64>| -> f64 {
                    (0..pin.rings)
                        .map(|j| weight(j) * x(&spectrum.nuclides[j][k])[g - top])
                        .sum::<f64>()
                        / total_weight
                };
                s.total[g] = average(|v| &v.total);
                s.elastic[g] = average(|v| &v.elastic);
                s.fission[g] = average(|v| &v.fission);
                s.nu_fission[g] = average(|v| &v.nu_fission);
            }
            s
        })
        .collect();
    Ok(Some(CellShielding {
        radius,
        rings,
        start,
        effective,
        sph,
    }))
}

/// 1本のピンの燃料の核種の実効断面積を、全断面積が収束するまで反復して求める。
fn shield_material(
    library: &Library,
    material: &MaterialInput,
    equivalence: &Equivalence,
    options: &ShieldingOptions,
    cell: usize,
) -> Result<Vec<Shielded>, SelfShieldingError> {
    let t = material.temperature;
    let ng = library.groups.len();
    let infinite = *library.sigma0.first().ok_or(LookupError::NoSigma0)?;
    let mut totals = material
        .composition
        .iter()
        .map(|(name, _)| Ok(library.nuclide(name, t)?.total))
        .collect::<Result<Vec<_>, SelfShieldingError>>()?;
    for _ in 0..options.max_iterations {
        let tables = material
            .composition
            .iter()
            .enumerate()
            .map(|(k, (name, density))| {
                // 個数密度が 0 の核種は何にも寄与しないので、無限希釈のままにする
                if *density == 0.0 {
                    return Ok(library.shielded(name, t, &vec![infinite; ng])?);
                }
                // 同じ燃料の他の核種の全断面積の和 (1/cm)
                let background: Vec<f64> = (0..ng)
                    .map(|g| {
                        (material.composition.iter().zip(&totals))
                            .enumerate()
                            .filter(|&(j, _)| j != k)
                            .map(|(_, ((_, n), sigma))| n * sigma[g])
                            .sum()
                    })
                    .collect();
                effective(library, name, t, *density, &background, equivalence)
            })
            .collect::<Result<Vec<_>, ShieldingError>>()?;
        let finite = |s: &Shielded| {
            [&s.total, &s.elastic, &s.fission, &s.nu_fission]
                .iter()
                .all(|v| v.iter().all(|x| x.is_finite()))
        };
        if !tables.iter().all(finite) {
            return Err(SelfShieldingError::NonFinite(cell));
        }
        let change = tables
            .iter()
            .zip(&totals)
            .flat_map(|(s, old)| s.total.iter().zip(old).map(|(a, b)| ((a - b) / b).abs()))
            .fold(0.0, f64::max);
        totals = tables.iter().map(|s| s.total.clone()).collect();
        if change <= options.tolerance {
            return Ok(tables);
        }
    }
    Err(SelfShieldingError::NotConverged {
        cell,
        iterations: options.max_iterations,
    })
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
    // 弾性散乱の変化を P0 の行の和に入れ、どの次数の行も同じ比で変える（角度分布の形は変えない）
    for g in 0..ng {
        let sum: f64 = nuclide.scatter[g].iter().sum();
        if sum > 0.0 {
            let factor = (sum + shielded.elastic[g] - dilute.elastic[g]) / sum;
            std::iter::once(&mut nuclide.scatter)
                .chain(&mut nuclide.anisotropic)
                .for_each(|m| m[g].iter_mut().for_each(|v| *v *= factor));
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
    use xlc_format::group::GroupStructure;
    use xlc_format::library::{Library, NuclideData, Shielded, TemperatureData};

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
            pointwise_range: Some([4.0, 1.0e4]),
            groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
            sigma0: vec![1.0e10, 100.0, 1.0],
            nuclides: vec![NuclideData {
                name: "U-238".into(),
                za: 92238.0,
                awr: 236.0,
                resolved: None,
                thermal_scattering: false,
                temperatures: vec![TemperatureData {
                    pointwise: Default::default(),
                    temperature: T,
                    chi: vec![0.0, 0.0],
                    scatter: vec![vec![vec![6.0, 2.0], vec![0.5, 15.0]]],
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
        use xlc_format::group::GroupStructure;
        use xlc_format::library::{Library, NuclideData, Resolved, Shielded, TemperatureData};

        const T: f64 = 600.0;

        /// 共鳴のある U-238（群 0 が σ0 に依存）、核分裂する U-235、水素
        fn library() -> Library {
            let data = |name: &str, tables: Vec<Shielded>, fissile: bool| NuclideData {
                name: name.into(),
                za: 0.0,
                awr: 1.0,
                resolved: None,
                thermal_scattering: false,
                temperatures: vec![TemperatureData {
                    pointwise: Default::default(),
                    temperature: T,
                    chi: if fissile {
                        vec![1.0, 0.0]
                    } else {
                        vec![0.0, 0.0]
                    },
                    scatter: vec![vec![vec![6.0, 1.0], vec![0.0, 4.0]]],
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
            let constant = |total, fission| vec![table(total, [6.0, 4.0], fission); 3];
            Library {
                pointwise_range: Some([4.0, 1.0e4]),
                groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
                sigma0: vec![1.0e10, 100.0, 1.0],
                nuclides: vec![
                    data(
                        "U-238",
                        vec![u238(30.0, 8.0), u238(20.0, 7.5), u238(9.0, 7.0)],
                        false,
                    ),
                    data("U-235", constant([10.0, 50.0], [2.0, 40.0]), true),
                    data("H-1", constant([20.0, 30.0], [0.0, 0.0]), false),
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
        fn shielding_scales_every_scattering_moment() {
            // 群 0 の弾性散乱が無限希釈の 8 から 7 に減ると、群 0 の散乱の行はどの次数も (7 - 1) / 7 倍
            let mut library = library();
            library.nuclides[0].temperatures[0]
                .scatter
                .push(vec![vec![2.0, 0.5], vec![0.0, 1.0]]);
            let table = library.shielded("U-238", T, &[1.0; 2]).unwrap();
            let nuclide = shielded_nuclide(&library, "U-238", T, &table).unwrap();
            let factor = 6.0 / 7.0;
            assert_eq!(nuclide.scatter[0], vec![6.0 * factor, factor]);
            assert_eq!(nuclide.anisotropic[0][0], vec![2.0 * factor, 0.5 * factor]);
            assert_eq!(nuclide.anisotropic[0][1], vec![0.0, 1.0]);
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
            assert_eq!(pin.dancoff.len(), 2);
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

        /// 分離共鳴を持つ（群 0 が σ0 に依存する）Zr を加えた `library()`
        fn library_with_zirconium() -> Library {
            let mut library = library();
            let mut zr = library.nuclides[2].clone();
            zr.name = "Zr-90".into();
            zr.resolved = Some(Resolved {
                top: 1.0e4,
                potential: 6.0,
            });
            zr.temperatures[0].shielded = [15.0, 10.0, 5.0]
                .iter()
                .map(|&t| Shielded {
                    total: vec![t, 8.0],
                    elastic: vec![6.0, 4.0],
                    fission: vec![0.0, 0.0],
                    nu_fission: vec![0.0, 0.0],
                })
                .collect();
            library.nuclides.push(zr);
            library
        }

        #[test]
        fn structural_frames_are_shielded_as_lumps() {
            let library = library_with_zirconium();
            let mut materials = materials();
            materials.push(MaterialInput {
                name: "zirc".into(),
                temperature: T,
                composition: vec![("Zr-90".into(), 0.04)],
            });
            // 領域: セルごとに (燃料, 水, 枠)
            let pin = PinCell::new(1.26, vec![0.4096])
                .unwrap()
                .with_frame(0.05)
                .unwrap();
            let lattice = Lattice::new(vec![pin], vec![vec![0, 0], vec![0, 0]]).unwrap();
            let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
            let regions: Vec<String> = (0..4)
                .flat_map(|_| ["fuel", "water", "zirc"].map(String::from))
                .collect();
            let shielded = shield_fuel(
                &library,
                &tracking,
                &materials,
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            let dilute = |name: &str| {
                let m = materials.iter().find(|m| m.name == name).unwrap();
                Material::mix(&library, &m.composition, T).unwrap()
            };
            let frame = shielded.region_material[2];
            assert!(shielded.materials[frame].total[0] < dilute("zirc").total[0]);
            // 対称に移り合う枠は同じ物質を共有する
            assert!((0..4).all(|c| shielded.region_material[3 * c + 2] == frame));
            // 分離共鳴を持たない水は無限希釈のまま
            assert_eq!(
                shielded.materials[shielded.region_material[1]],
                dilute("water")
            );
        }

        /// `library_with_zirconium()` に、Zr-90 と同じ群定数で分離共鳴を持つ核種 `resonant` と持たない核種
        /// `plain` を加えたもの
        fn library_with_copies(resonant: &str, plain: &str) -> Library {
            let mut library = library_with_zirconium();
            let zr = library.nuclides[3].clone();
            let mut copy = |name: &str, resolved: bool| {
                let mut n = zr.clone();
                n.name = name.into();
                if !resolved {
                    n.resolved = None;
                }
                library.nuclides.push(n);
            };
            copy(resonant, true);
            copy(plain, false);
            library
        }

        /// 枠（物質 zirc、厚さ 0.05）付きのセルの種類（半径と、内側からの物質の名前）を `layout` に並べた格子と、
        /// 領域ごとの物質の名前
        fn framed(
            types: &[(Vec<f64>, Vec<&str>)],
            layout: Vec<Vec<usize>>,
        ) -> (Tracking, Vec<String>) {
            let cells = types
                .iter()
                .map(|(radii, _)| {
                    PinCell::new(1.26, radii.clone())
                        .unwrap()
                        .with_frame(0.05)
                        .unwrap()
                })
                .collect();
            let regions = layout
                .iter()
                .flatten()
                .flat_map(|&t| types[t].1.iter().chain(&["zirc"]).map(|n| n.to_string()))
                .collect();
            let lattice = Lattice::new(cells, layout).unwrap();
            (Tracking::new(&lattice, 16, 0.02).unwrap(), regions)
        }

        fn with_material(name: &str, nuclide: &str) -> Vec<MaterialInput> {
            let mut materials = materials();
            for (name, nuclide) in [("zirc", "Zr-90"), (name, nuclide)] {
                materials.push(MaterialInput {
                    name: name.into(),
                    temperature: T,
                    composition: vec![(nuclide.into(), 0.04)],
                });
            }
            materials
        }

        #[test]
        fn reports_each_lump() {
            let library = library_with_zirconium();
            let types = [(vec![0.4096], vec!["fuel", "water"])];
            let (tracking, regions) = framed(&types, vec![vec![0, 0], vec![0, 0]]);
            let shielded = shield_fuel(
                &library,
                &tracking,
                &with_material("unused", "H-1"),
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            // 対称に移り合う4つの枠は1つの塊。露出した面は枠の内側だけ
            assert_eq!(shielded.lumps.len(), 1, "{:?}", shielded.lumps);
            let lump = &shielded.lumps[0];
            assert_eq!((lump.cell, lump.material.as_str()), (0, "zirc"));
            assert_eq!(lump.regions, vec![2, 5, 8, 11]);
            let (p, t): (f64, f64) = (1.26, 0.05);
            let chord = (p * p - (p - 2.0 * t).powi(2)) / (p - 2.0 * t);
            assert!((lump.chord - chord).abs() < 1e-2 * chord, "{}", lump.chord);
            assert_eq!(lump.dancoff.len(), 2);
            // 実効断面積は、塊の物質の核種の全断面積
            let frame = &shielded.materials[shielded.region_material[2]];
            assert!((frame.total[0] - 0.04 * lump.effective[0].total[0]).abs() < 1e-12);
        }

        #[test]
        fn touching_pieces_are_one_lump() {
            // 管（外径 0.6）は枠の内側（半辺 0.58）に届き、同じ物質の枠と1つの塊になる
            let library = library_with_zirconium();
            let types = [
                (vec![0.4096], vec!["fuel", "water"]),
                (vec![0.5, 0.6], vec!["water", "zirc", "water"]),
            ];
            let (tracking, regions) = framed(&types, vec![vec![0, 1], vec![1, 0]]);
            let shielded = shield_fuel(
                &library,
                &tracking,
                &with_material("unused", "H-1"),
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            // セル1の領域: 水, 管, 水, 枠
            assert_eq!(shielded.region_material[4], shielded.region_material[6]);
        }

        #[test]
        fn touching_materials_with_the_same_resonant_nuclides_are_one_body() {
            // 同じ核種で密度の違う物質 zirc2 を zirc が囲む。2つは1つの黒体で、Dancoff 係数と平均弦長を共有する
            let library = library_with_zirconium();
            let types = [
                (vec![0.4096], vec!["fuel", "water"]),
                (
                    vec![0.3, 0.35, 0.4, 0.45],
                    vec!["water", "zirc", "zirc2", "zirc", "water"],
                ),
            ];
            let (tracking, regions) = framed(&types, vec![vec![0, 1], vec![1, 0]]);
            let shielded = shield_fuel(
                &library,
                &tracking,
                &with_material("zirc2", "Zr-90"),
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            let lump = |material: &str| {
                shielded
                    .lumps
                    .iter()
                    .find(|l| l.material == material && l.cell == 1)
                    .unwrap()
            };
            let (outer, inner) = (lump("zirc"), lump("zirc2"));
            assert_eq!(outer.dancoff, inner.dancoff);
            assert_eq!(outer.chord, inner.chord);
            assert_eq!(outer.regions, vec![4, 6, 10, 12]);
            assert_eq!(inner.regions, vec![5, 11]);
        }

        #[test]
        fn other_resonant_materials_do_not_shadow_a_lump() {
            // 管の物質が分離共鳴を持つかどうかで、同じセルの枠（別の物質）の遮蔽は変わらない（群定数は同じ）
            let library = library_with_copies("Fe-56", "Fe-56x");
            let types = [
                (vec![0.4096], vec!["fuel", "water"]),
                (vec![0.45, 0.5], vec!["water", "tube", "water"]),
            ];
            let (tracking, regions) = framed(&types, vec![vec![0, 1], vec![1, 0]]);
            let frame = |nuclide: &str| {
                let shielded = shield_fuel(
                    &library,
                    &tracking,
                    &with_material("tube", nuclide),
                    &regions,
                    &ShieldingOptions::default(),
                )
                .unwrap();
                // セル1の領域: 水, 管, 水, 枠
                shielded.materials[shielded.region_material[6]].total[0]
            };
            let (resonant, plain) = (frame("Fe-56"), frame("Fe-56x"));
            assert!(
                (resonant - plain).abs() < 1e-12 * plain,
                "{resonant} vs {plain}"
            );
        }

        #[test]
        fn separate_pieces_of_a_material_are_separate_lumps() {
            // 管と枠が同じ物質でも、離れた領域は別の塊として形に合った平均弦長で遮蔽する
            let library = library_with_zirconium();
            let types = [
                (vec![0.4096], vec!["fuel", "water"]),
                (vec![0.45, 0.5], vec!["water", "zirc", "water"]),
            ];
            let (tracking, regions) = framed(&types, vec![vec![0, 1], vec![1, 0]]);
            let shielded = shield_fuel(
                &library,
                &tracking,
                &with_material("unused", "H-1"),
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            // セル1の領域: 水, 管, 水, 枠
            let (tube, frame) = (shielded.region_material[4], shielded.region_material[6]);
            assert_ne!(tube, frame);
            let total = |m: usize| shielded.materials[m].total[0];
            assert!(total(tube) != total(frame), "{}", total(tube));
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

        /// 分離共鳴（200 eV の共鳴）を持つ U-238 と、ポイントワイズの断面積を持つ核種のライブラリ。群 1, 2 が UFG の範囲
        fn resonant_library() -> Library {
            let bounds = vec![1.0e7, 1.0e4, 1.0e2, 4.0, 0.625, 1.0e-5];
            let ng = bounds.len() - 1;
            let resonance = |e: f64| 1.0e4 / (1.0 + ((e - 200.0) / 0.5).powi(2));
            let mut energies: Vec<f64> = (0..=400).map(|k| 180.0 + 0.1 * k as f64).collect();
            energies.splice(0..0, [4.0, 100.0]);
            energies.extend([1.0e3, 1.0e4]);
            // (名前, awr, ポイントワイズの断面積 [全, 弾性散乱, 核分裂], 分離共鳴の上端)
            let nuclide = |name: &str, awr: f64, points: Vec<[f64; 4]>, top: Option<f64>| {
                let column = |k: usize| points.iter().map(|p| p[k]).collect::<Vec<_>>();
                let (total, elastic, fission) = (points[0][1], points[0][2], points[0][3]);
                // σ0 の表: 無限希釈から下がる（Carlvik で引く範囲の外の群に使う）
                let table = |factor: f64| Shielded {
                    total: vec![total * factor; ng],
                    elastic: vec![elastic; ng],
                    fission: vec![fission * factor; ng],
                    nu_fission: vec![2.4 * fission * factor; ng],
                };
                NuclideData {
                    name: name.into(),
                    za: 0.0,
                    awr,
                    resolved: top.map(|top| xlc_format::library::Resolved {
                        top,
                        potential: 10.0,
                    }),
                    thermal_scattering: false,
                    temperatures: vec![TemperatureData {
                        temperature: T,
                        chi: (0..ng)
                            .map(|g| if fission > 0.0 && g == 0 { 1.0 } else { 0.0 })
                            .collect(),
                        // 弾性散乱は自群に 9 割、次の群に 1 割（行の和が弾性散乱。最後の群は自群だけ）
                        scatter: vec![
                            (0..ng)
                                .map(|g| {
                                    (0..ng)
                                        .map(|h| match (h == g, h == g + 1, g + 1 == ng) {
                                            (true, _, true) => elastic,
                                            (true, _, false) => 0.9 * elastic,
                                            (_, true, _) => 0.1 * elastic,
                                            _ => 0.0,
                                        })
                                        .collect()
                                })
                                .collect(),
                        ],
                        shielded: vec![table(1.0), table(0.9), table(0.8)],
                        pointwise: xlc_format::library::Pointwise {
                            energy: column(0),
                            total: column(1),
                            elastic: column(2),
                            fission: column(3),
                            nu_fission: column(3).iter().map(|f| 2.4 * f).collect(),
                        },
                    }],
                }
            };
            let flat = |t: f64, e: f64, f: f64| vec![[4.0, t, e, f], [1.0e4, t, e, f]];
            let u238 = energies
                .iter()
                .map(|&e| [e, 11.0 + resonance(e), 11.0, 0.0])
                .collect();
            Library {
                pointwise_range: Some([4.0, 1.0e4]),
                groups: GroupStructure::new(bounds).unwrap(),
                sigma0: vec![1.0e10, 100.0, 1.0],
                nuclides: vec![
                    nuclide("U-238", 236.0, u238, Some(1.0e4)),
                    nuclide("U-235", 233.0, flat(10.0, 8.0, 1.0), Some(2.0e3)),
                    nuclide("O-16", 15.86, flat(3.8, 3.8, 0.0), None),
                    nuclide("Zr-90", 89.0, flat(6.5, 6.5, 0.0), None),
                    nuclide("H-1", 0.99917, flat(20.0, 20.0, 0.0), None),
                ],
            }
        }

        /// 被覆管が物質 `clad`、枠（厚さ 0.05）が分離共鳴を持つ Zr-91 の物質 grid の、UFG で遮蔽するピンセル。
        /// 枠の塊の Dancoff 係数と、`black` の領域を黒体にして求め直した Dancoff 係数（燃料はピンの SPH 因子を
        /// 掛けない実効断面積、ほかは無限希釈の全断面積）
        fn frame_dancoff(clad: &str, black: &[usize]) -> (Vec<f64>, Vec<f64>) {
            let mut library = resonant_library();
            let mut zr91 = library.nuclides[3].clone();
            zr91.name = "Zr-91".into();
            zr91.resolved = Some(xlc_format::library::Resolved {
                top: 1.0e3,
                potential: 6.5,
            });
            library.nuclides.push(zr91);
            let material = |name: &str, composition: &[(&str, f64)]| MaterialInput {
                name: name.into(),
                temperature: T,
                composition: composition.iter().map(|&(n, d)| (n.into(), d)).collect(),
            };
            let materials = vec![
                material(
                    "fuel",
                    &[("U-238", 0.022), ("U-235", 5.0e-4), ("O-16", 0.045)],
                ),
                material("clad", &[("Zr-90", 0.043)]),
                material("grid", &[("Zr-91", 0.043)]),
                material("grid2", &[("Zr-91", 0.03)]),
                material("water", &[("H-1", 0.05), ("O-16", 0.025)]),
            ];
            let pin = PinCell::new(1.26, vec![0.2, 0.35, 0.4096, 0.418, 0.475])
                .unwrap()
                .with_frame(0.05)
                .unwrap();
            let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.02).unwrap();
            let regions: Vec<String> = [
                "fuel",
                "fuel",
                "fuel",
                crate::input::VOID,
                clad,
                "water",
                "grid",
            ]
            .map(String::from)
            .to_vec();
            let shielded = shield_fuel(
                &library,
                &tracking,
                &materials,
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            let pin = &shielded.pins[0];
            assert!(pin.sph.is_some());
            let fuel = Material::from_nuclides(&pin.nuclides, library.groups.len()).unwrap();
            let totals: Vec<Vec<f64>> = regions
                .iter()
                .map(|name| match name.as_str() {
                    "fuel" => fuel.total.clone(),
                    crate::input::VOID => vec![0.0; library.groups.len()],
                    name => {
                        let m = materials.iter().find(|m| m.name == name).unwrap();
                        Material::mix(&library, &m.composition, T).unwrap().total
                    }
                })
                .collect();
            let lump: Vec<Option<usize>> = (0..regions.len())
                .map(|r| black.contains(&r).then_some(usize::from(r != 6)))
                .collect();
            let lumps = 1 + usize::from(black.len() > 1);
            let black: Vec<bool> = lump.iter().map(Option::is_some).collect();
            let expected = lump_dancoff(
                &tracking,
                &totals,
                &black,
                &lump,
                lumps,
                &DancoffOptions::default(),
            )
            .unwrap();
            assert_eq!(shielded.lumps.len(), 1, "{:?}", shielded.lumps);
            (
                shielded.lumps[0].dancoff.clone(),
                expected[0].dancoff.clone(),
            )
        }

        #[test]
        fn lumps_see_the_physical_fuel() {
            // 燃料は SPH 因子を掛けない実効断面積で、塊の影の計算に入る
            let (dancoff, expected) = frame_dancoff("clad", &[6]);
            for (c, e) in dancoff.iter().zip(&expected) {
                assert!((c - e).abs() < 1e-9, "{dancoff:?} vs {expected:?}");
            }
        }

        #[test]
        fn same_material_in_the_pin_shadows_a_lump() {
            // 被覆管が枠と同じ物質なら、等価なセルで遮蔽した被覆管も、枠と同じ共鳴の影を作る（黒体）
            let (dancoff, expected) = frame_dancoff("grid", &[4, 6]);
            for (c, e) in dancoff.iter().zip(&expected) {
                assert!((c - e).abs() < 1e-9, "{dancoff:?} vs {expected:?}");
            }
        }

        #[test]
        fn materials_with_the_same_resonant_nuclides_shadow_a_lump() {
            // 被覆管が枠と別の物質でも、分離共鳴を持つ核種が同じなら同じ共鳴の影を作る（黒体）
            let (dancoff, expected) = frame_dancoff("grid2", &[4, 6]);
            for (c, e) in dancoff.iter().zip(&expected) {
                assert!((c - e).abs() < 1e-9, "{dancoff:?} vs {expected:?}");
            }
        }

        #[test]
        fn pieces_of_a_lump_follow_contact() {
            // 小さな正方形のセルの縦の列と、扇形の角度の 0 と 2π をまたぐ扇形は、それぞれ1つの塊
            let library = library_with_zirconium();
            let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
            let grid = PinCell::grid(1.26, 3).unwrap();
            let sectored = PinCell::with_sectors(1.26, vec![0.5, 0.6], 4).unwrap();
            let layout = vec![vec![0, 1, 0], vec![2, 0, 0], vec![0, 0, 0]];
            let names = |t: usize| -> Vec<String> {
                let names = match t {
                    0 => "fuel water",
                    1 => "water zirc water water zirc water water zirc water",
                    _ => "water water water water zirc water water zirc water water water water",
                };
                names.split(' ').map(String::from).collect()
            };
            let regions: Vec<String> = layout.iter().flatten().flat_map(|&t| names(t)).collect();
            let lattice = Lattice::new(vec![pin, grid, sectored], layout).unwrap();
            let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
            let shielded = shield_fuel(
                &library,
                &tracking,
                &with_material("unused", "H-1"),
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            let m = &shielded.region_material;
            // 正方形のセル（セル1）は領域 2〜10、扇形のセル（セル3）は 13〜24（リング 1 の扇形は 17〜20）
            assert!(m[3] == m[6] && m[6] == m[9], "{m:?}");
            assert_eq!(m[17], m[20], "{m:?}");
            assert_eq!(shielded.lumps.len(), 2, "{:?}", shielded.lumps);
        }

        #[test]
        fn fuel_rings_are_shielded_by_the_cell() {
            let library = resonant_library();
            let materials = vec![
                MaterialInput {
                    name: "fuel".into(),
                    temperature: T,
                    composition: vec![
                        ("U-238".into(), 0.022),
                        ("U-235".into(), 5.0e-4),
                        ("O-16".into(), 0.045),
                    ],
                },
                MaterialInput {
                    name: "clad".into(),
                    temperature: T,
                    composition: vec![("Zr-90".into(), 0.043)],
                },
                MaterialInput {
                    name: "water".into(),
                    temperature: T,
                    composition: vec![("H-1".into(), 0.05), ("O-16".into(), 0.025)],
                },
            ];
            let pin = PinCell::new(1.26, vec![0.2, 0.35, 0.4096, 0.418, 0.475]).unwrap();
            let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.02).unwrap();
            let regions: Vec<String> =
                ["fuel", "fuel", "fuel", crate::input::VOID, "clad", "water"]
                    .map(String::from)
                    .to_vec();
            let shielded = shield_fuel(
                &library,
                &tracking,
                &materials,
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            let total = |region: usize, g: usize| {
                shielded.materials[shielded.region_material[region]].total[g]
            };
            // UFG の範囲の群（100 eV〜10 keV）では、燃料の内側ほど共鳴が遮蔽される
            assert!(
                total(0, 1) < total(1, 1) && total(1, 1) < total(2, 1),
                "{:?}",
                [total(0, 1), total(1, 1), total(2, 1)]
            );
            // 範囲の外の群は、燃料のリングによらず等価原理の値
            assert_eq!(total(0, 0), total(2, 0));
            // 等価なセルの外径は、正方形のセルと面積の等しい円に近い
            let radius = shielded.pins[0].equivalent.unwrap();
            assert!((0.6..0.8).contains(&radius), "{radius}");
            assert!(shielded.materials[shielded.region_material[3]].is_void());
            // SPH 因子は等価なセルのリング（燃料3、ギャップ、被覆管、減速材）と UFG の範囲の群（1, 2）ごと。
            // 真空のギャップは 1、燃料では 1 と違う
            let sph = shielded.pins[0].sph.as_ref().unwrap();
            assert_eq!(sph.len(), 6);
            assert!(sph.iter().all(|row| row.len() == 2));
            assert_eq!(sph[3], vec![1.0, 1.0]);
            assert!(sph[0].iter().any(|mu| (mu - 1.0).abs() > 1e-6), "{sph:?}");
        }

        #[test]
        fn symmetric_pins_are_solved_once() {
            // 3×3 の中央が水のセル: 角の 4 本と辺の 4 本がそれぞれ同じ組になり、等価なセルの減速計算は組ごとに1回
            let library = resonant_library();
            let materials = vec![
                MaterialInput {
                    name: "fuel".into(),
                    temperature: T,
                    composition: vec![
                        ("U-238".into(), 0.022),
                        ("U-235".into(), 5.0e-4),
                        ("O-16".into(), 0.045),
                    ],
                },
                MaterialInput {
                    name: "clad".into(),
                    temperature: T,
                    composition: vec![("Zr-90".into(), 0.043)],
                },
                MaterialInput {
                    name: "water".into(),
                    temperature: T,
                    composition: vec![("H-1".into(), 0.05), ("O-16".into(), 0.025)],
                },
            ];
            let pin = PinCell::new(1.26, vec![0.2, 0.35, 0.4096, 0.418, 0.475]).unwrap();
            let water = PinCell::new(1.26, vec![]).unwrap();
            let layout = vec![vec![0, 0, 0], vec![0, 1, 0], vec![0, 0, 0]];
            let lattice = Lattice::new(vec![pin, water], layout).unwrap();
            let tracking = Tracking::new(&lattice, 16, 0.05).unwrap();
            let pin_regions = ["fuel", "fuel", "fuel", crate::input::VOID, "clad", "water"];
            let regions: Vec<String> = (0..9)
                .flat_map(|cell| {
                    if cell == 4 {
                        vec!["water"]
                    } else {
                        pin_regions.to_vec()
                    }
                })
                .map(String::from)
                .collect();
            let shielded = shield_fuel(
                &library,
                &tracking,
                &materials,
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            // 希釈した物質（3 と真空）に、2 組 × 真空でない 5 リングの物質
            assert_eq!(shielded.materials.len(), 4 + 2 * 5);
            // 角どうし、辺どうしは同じ物質を共有する。燃料のリングの最初の領域で比べる
            let first = |cell: usize| shielded.region_material[lattice.cell_regions(cell).start];
            assert_eq!(first(0), first(8));
            assert_eq!(first(1), first(7));
            assert_ne!(first(0), first(1));
            // ピンごとの結果は 8 本とも残り、組の中では同じ
            assert_eq!(shielded.pins.len(), 8);
            let radius = |cell: usize| {
                shielded
                    .pins
                    .iter()
                    .find(|p| p.cell == cell)
                    .unwrap()
                    .equivalent
                    .unwrap()
            };
            assert_eq!(radius(0), radius(8));
            assert_eq!(radius(1), radius(3));
            // 報告する Dancoff 係数も、外径を決めた代表の値
            let dancoff = |cell: usize| {
                &shielded
                    .pins
                    .iter()
                    .find(|p| p.cell == cell)
                    .unwrap()
                    .dancoff
            };
            assert_eq!(dancoff(0), dancoff(8));
        }

        #[test]
        fn isolated_pins_are_shielded() {
            // 水が厚く、Dancoff 係数がセルの数値誤差より小さいピンも、孤立したピンとして解ける
            let library = resonant_library();
            let materials = vec![
                MaterialInput {
                    name: "fuel".into(),
                    temperature: T,
                    composition: vec![
                        ("U-238".into(), 0.022),
                        ("U-235".into(), 5.0e-4),
                        ("O-16".into(), 0.045),
                    ],
                },
                MaterialInput {
                    name: "water".into(),
                    temperature: T,
                    composition: vec![("H-1".into(), 0.05), ("O-16".into(), 0.025)],
                },
            ];
            let pin = PinCell::new(30.0, vec![0.4096]).unwrap();
            let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.05).unwrap();
            let regions: Vec<String> = ["fuel", "water"].map(String::from).to_vec();
            let shielded = shield_fuel(
                &library,
                &tracking,
                &materials,
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            let radius = shielded.pins[0].equivalent.unwrap();
            assert!(radius > 3.0, "{radius}");
        }

        #[test]
        fn rejects_fuel_outside_the_inner_rings() {
            let library = library();
            let (tracking, mut regions) = geometry(vec![vec![0]]);
            regions.swap(0, 1);
            let result = shield_fuel(
                &library,
                &tracking,
                &materials(),
                &regions,
                &ShieldingOptions::default(),
            );
            assert!(
                matches!(result, Err(SelfShieldingError::FuelShape(0))),
                "{result:?}"
            );
        }
    }
}
