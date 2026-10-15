//! 燃料のピンの自己遮蔽。
//!
//! 燃料の分離共鳴の範囲（[`crate::ufg::cell_range`]）は、ピンごとに Dancoff 係数を合わせた円筒の等価なセル
//! （[`crate::equivalent`]）の超微細群の減速計算（[`crate::ufg`]）で、燃料のリング、被覆管、減速材の群定数を
//! それぞれのリングの中性子束で縮約して求める。その範囲の外の群は、燃料を等価原理（Dancoff 係数を入れた
//! Carlvik の2項有理近似）で遮蔽し、被覆管と減速材は無限希釈とする。
//!
//! 燃料のピンの等価なセルに入らない構造材（支持格子の枠、案内管など）のうち分離共鳴を持つ核種を含むものは、
//! 塊ごとに Dancoff 係数と平均弦長を求め、燃料と同じ等価原理で遮蔽する（[`shield_lumps`]）。ピンセルの内側のリングを
//! 埋める塊（制御棒の吸収体などの棒）は、UFG の範囲を燃料のピンと同じく等価なセルの減速計算で遮蔽する（ADR-0142）。
//!
//! 等価原理では、燃料の脱出確率を P(x) = Σ_i β_i x / (x + α_i)（x = Σt l̄）で近似し、核種 k の背景断面積を
//! σ0,i = (Σ_{j≠k} N_j σt,j + α_i / l̄) / N_k とする。実効断面積は、2つの σ0 でのライブラリの値を、
//! NR 近似の群平均の中性子束 φ_i = (σp + σ0,i) / (σt(σ0,i) + σ0,i) と β_i で重みづけして合わせる。σp は核種の
//! ポテンシャル散乱で、NR 近似の中性子束 (σp + σ0) w / (σt + σ0) の線源に入る（共鳴の間で中性子束が w になる）。
//! ライブラリの重みが w/(σt + σ0) の群（Bondarenko）では、この式は NR 近似の範囲で厳密に成り立つ。分離共鳴の
//! 範囲の群は flux calculator の中性子束で重みづけしてあり、この式は近似になる（ADR-0089）。
//! σp は分離共鳴の範囲の s 波の散乱半径の値をすべての群に使う（非分離共鳴の範囲の p 波などは入れない）。

use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::hash_map::Entry;

use xlc_format::library::{Library, LookupError, Nuclide, Shielded};

use crate::cp::{CpError, Cylinder};
use crate::equivalent::{EquivalentError, MIN_DANCOFF, OuterRing, equivalent_radius, outer_radius};
use crate::geometry::{Lattice, PinCell};
use crate::input::{MaterialInput, VOID};
use crate::material::{Material, MaterialError};
use crate::moc::dancoff::{DancoffError, DancoffOptions, Lump, dancoff, lump_dancoff};
use crate::sph::{SphError, sph_factors};
use crate::symmetry::{Detected, Image, SymmetryError, detect};
use crate::tracking::Tracking;
use crate::ufg::{
    CellSpectrum, Collision, UfgError, cell_range, material_range, slowing_down_with,
};

/// 格子の燃料の自己遮蔽の反復の条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShieldingOptions {
    /// 燃料の各核種の全断面積の相対変化の上限
    pub tolerance: f64,
    pub max_iterations: usize,
    pub dancoff: DancoffOptions,
    /// 燃料のピンと棒の UFG の減速計算の飛程の求め方
    pub collision: Collision,
}

impl Default for ShieldingOptions {
    fn default() -> Self {
        Self {
            tolerance: 1.0e-6,
            max_iterations: 100,
            dancoff: DancoffOptions::default(),
            collision: Collision::default(),
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
    /// 解くと決めたピンの結果がない（内部の不整合）
    #[error("cell {0}: internal error, the pin was not shielded")]
    Unsolved(usize),
    #[error(
        "cell {0}: the fuel rings differ in nuclides, their order or temperature, or a ring mixes materials"
    )]
    MixedFuel(usize),
    #[error(
        "cell {0}: each ring outside the fuel must be one material, and the outer rings the moderator"
    )]
    CellShape(usize),
    #[error(
        "region {0}: a symmetry operation maps its structural lump onto several lumps (internal error)"
    )]
    LumpImage(usize),
    #[error("structural body {0} has no Dancoff factor")]
    Body(usize),
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
    #[error("the symmetry was verified for another lattice or tracking")]
    ForeignSymmetry,
    #[error("a flux of {given} groups for {groups} groups")]
    FluxGroups { given: usize, groups: usize },
    #[error("cell {0}: the SPH factors do not match the rings and groups of the cell")]
    SphLayout(usize),
    #[error("the lattice is not the one that was self-shielded")]
    ForeignLattice,
    #[error(transparent)]
    Shielding(#[from] ShieldingError),
}

impl From<LookupError> for SelfShieldingError {
    fn from(error: LookupError) -> Self {
        Self::Shielding(ShieldingError::Lookup(error))
    }
}

/// 自己遮蔽した物質。燃料はピンの組ごとに、構造材の塊は対称操作で移り合う塊の組と物質ごとに、それ以外は物質ごとに
/// 1つずつ持つ。
#[derive(Debug, Clone, PartialEq)]
pub struct FuelShielding {
    pub materials: Vec<Material>,
    /// [領域] の `materials` の番号
    pub region_material: Vec<usize>,
    /// 燃料のあるセルごとの自己遮蔽の結果（セルの順）
    pub pins: Vec<PinShielding>,
    /// 構造材の塊ごとの自己遮蔽の結果
    pub lumps: Vec<LumpShielding>,
    /// 自己遮蔽した格子（リングと領域の対応を、同じ格子でだけ使う）
    lattice: Lattice,
}

/// 構造材の1つの塊（[`shield_lumps`]）の自己遮蔽の結果。
#[derive(Debug, Clone, PartialEq)]
pub struct LumpShielding {
    /// 塊のあるセル（対称操作で移り合う塊をまとめた組では、組の代表の領域（最も小さい番号）のセル。ADR-0112）
    pub cell: usize,
    pub material: String,
    /// 塊の領域（組のすべてのセルの分）
    pub regions: Vec<usize>,
    /// 群ごとの Dancoff 係数
    pub dancoff: Vec<f64>,
    /// 平均弦長 (cm)
    pub chord: f64,
    /// 棒（ピンセルの内側のリングを埋める塊）を UFG の範囲で遮蔽した等価なセルの減速材のリングの外半径 (cm)。
    /// 棒でない塊と、分離共鳴を持つ核種のない棒は `None`（ADR-0142）
    pub equivalent: Option<f64>,
    /// 棒のリングの SPH 因子（[棒のリング][群 − `groups.start`]）と、その群。棒でない塊は `None`
    pub sph: Option<(Vec<Vec<f64>>, std::ops::Range<usize>)>,
    /// 物質の核種ごとの実効断面積（物質の組成の順）。棒の UFG の範囲の群は、棒のリングの減速計算の平均
    pub effective: Vec<Shielded>,
}

impl LumpShielding {
    /// 塊の領域 `region` の、セルの中のリングの番号
    fn ring_of(lattice: &Lattice, region: usize) -> Option<usize> {
        let home = lattice.region_cell(region);
        let sectors = lattice.cell_type(home).pin()?.sectors();
        Some((region - lattice.cell_regions(home).start) / sectors)
    }
}

/// 1本のピンの燃料の自己遮蔽の結果。格子の対称性で移り合うピンは、組の代表（最も小さいセル番号）の結果を
/// 共有する（`cell` 以外は代表の値）。
#[derive(Debug, Clone, PartialEq)]
pub struct PinShielding {
    pub cell: usize,
    /// 群ごとの Dancoff 係数
    pub dancoff: Vec<f64>,
    /// 燃料の核種の自己遮蔽した群定数 (barn) と個数密度 (1/(barn cm))。燃料のリングの物質が違うときは、個数密度は
    /// リングの体積で平均したもの（等価原理はその組成で解く。ADR-0173）
    pub nuclides: Vec<(Nuclide, f64)>,
    /// 核種ごとの実効断面積（`nuclides` の順。弾性散乱を含む）。UFG の範囲の群は、燃料のリングの値を中性子束と
    /// 体積で平均したもの。`nuclides` とともに SPH 因子を掛けていない物理的な値（参照解と比べる用）で、MOC の
    /// 中性子束と組み合わせて反応率を求めるには、領域の物質（SPH 因子を掛けたもの）を使う
    pub effective: Vec<Shielded>,
    /// [燃料のリング（内側から）][核種（`nuclides` の順）] の実効断面積。UFG の範囲の群はそのリングの減速計算の値、範囲の外は
    /// 等価原理の値。`effective` と同じく SPH 因子を掛けない物理的な値で、物理的な中性子束（`FuelShielding::physical_flux`）と
    /// 掛けてリングの反応率にする（ADR-0173）
    pub rings: Vec<Vec<Shielded>>,
    /// UFG の減速計算に使った等価なセルの減速材のリングの外半径 (cm)。枠のあるセルでは、その外に枠と同じ面積の
    /// 枠のリングがある。燃料に分離共鳴を持つ核種がなければ `None`
    pub equivalent: Option<f64>,
    /// SPH 因子（[`crate::sph`]）。燃料に分離共鳴を持つ核種がなく UFG の範囲がなければ `None`
    pub sph: Option<Sph>,
}

/// 等価なセルのリングごと・UFG の範囲の群ごとの SPH 因子。その群の、リングの物質を置いた格子の領域の物質に掛けてある。
#[derive(Debug, Clone, PartialEq)]
pub struct Sph {
    /// [リング][群 − `groups.start`]。枠のあるセルの最後は枠のリングの値で、格子には使わないので 1
    pub factors: Vec<Vec<f64>>,
    /// `factors` の列の群（ライブラリの群の番号）
    pub groups: std::ops::Range<usize>,
    /// 等価なセルの減速材のリングの番号（リングを格子の領域に対応させる）
    moderator_ring: usize,
}

impl FuelShielding {
    /// MOC の中性子束 `flux`（[領域][群]）を物理的な中性子束にする。SPH 因子 μ を掛けた領域と群では、MOC の中性子束は
    /// 物理的な中性子束の 1/μ になる（反応率は μσ × φ/μ で物理的な値のまま）ので、μ を掛けて戻す（燃料のピンの等価なセルのリングと、棒の塊のリング）
    pub fn physical_flux(
        &self,
        lattice: &Lattice,
        flux: &[Vec<f64>],
    ) -> Result<Vec<Vec<f64>>, SelfShieldingError> {
        if lattice != &self.lattice {
            return Err(SelfShieldingError::ForeignLattice);
        }
        if flux.len() != lattice.regions() {
            return Err(SelfShieldingError::Regions {
                given: flux.len(),
                regions: lattice.regions(),
            });
        }
        let groups = self.materials.first().map_or(0, Material::groups);
        if let Some(f) = flux.iter().find(|f| f.len() != groups) {
            return Err(SelfShieldingError::FluxGroups {
                given: f.len(),
                groups,
            });
        }
        let mut physical = flux.to_vec();
        for pin in &self.pins {
            let Some(sph) = &pin.sph else {
                continue;
            };
            let layout = || SelfShieldingError::SphLayout(pin.cell);
            let rings = pin.ring_regions(lattice).ok_or_else(layout)?;
            // 因子はリングごと。枠のあるセルでは最後に、格子に使わない枠のリングの行（因子は 1）が 1 つ多い
            let framed = lattice
                .cell_type(pin.cell)
                .pin()
                .is_some_and(|p| p.frame().is_some());
            if sph.factors.len() != rings.len() + usize::from(framed) {
                return Err(layout());
            }
            let (used, frame) = sph.factors.split_at(rings.len());
            let width = sph.groups.len();
            if sph.groups.end > groups
                || sph.factors.iter().any(|row| row.len() != width)
                || frame.iter().flatten().any(|&mu| mu != 1.0)
            {
                return Err(layout());
            }
            for (factors, regions) in used.iter().zip(&rings) {
                for &region in regions {
                    for (g, mu) in sph.groups.clone().zip(factors) {
                        physical[region][g] *= mu;
                    }
                }
            }
        }
        // 棒のリングの因子（ADR-0142）
        for lump in &self.lumps {
            let Some((factors, range)) = &lump.sph else {
                continue;
            };
            let layout = || SelfShieldingError::SphLayout(lump.cell);
            if range.end > groups {
                return Err(layout());
            }
            for &region in &lump.regions {
                let ring = LumpShielding::ring_of(lattice, region).ok_or_else(layout)?;
                let row = factors.get(ring).ok_or_else(layout)?;
                if row.len() != range.len() {
                    return Err(layout());
                }
                for (g, mu) in range.clone().zip(row) {
                    physical[region][g] *= mu;
                }
            }
        }
        Ok(physical)
    }
}

impl PinShielding {
    /// SPH 因子のリングごとに、そのリングの物質を置いた格子 `lattice` の領域（枠のリングは格子に使わないので含めない）。
    /// UFG の範囲がないか、セルが同心円のセルでなければ `None`
    pub fn ring_regions(&self, lattice: &Lattice) -> Option<Vec<Vec<usize>>> {
        let sph = self.sph.as_ref()?;
        let pin = lattice.cell_type(self.cell).pin()?;
        Some(ring_targets(
            pin,
            lattice.cell_regions(self.cell).start,
            sph.moderator_ring,
        ))
    }
}

/// 格子の燃料（核分裂する核種を含む物質）のピンを自己遮蔽する（モジュールの説明を参照）。燃料のないセルは
/// 無限希釈のまま。`region_materials` は領域ごとの物質の名前。移り合うピンと塊の組は、領域ごとの物質の名前で成り立つ
/// 対称操作（`symmetry::detect`）で作る。成り立つ対称操作をすでに求めてあれば `shield_fuel_with` を使う
pub fn shield_fuel(
    library: &Library,
    tracking: &Tracking,
    materials: &[MaterialInput],
    region_materials: &[String],
    options: &ShieldingOptions,
) -> Result<FuelShielding, SelfShieldingError> {
    let symmetry = detect(tracking, region_materials)?;
    shield_fuel_with(
        library,
        tracking,
        materials,
        region_materials,
        &symmetry,
        options,
    )
}

/// `shield_fuel` と同じ。移り合うピンと塊の組は、確かめた対称操作 `symmetry`（同じトラッキングと領域ごとの物質の名前で
/// `symmetry::detect` が求めたもの）で作る（ADR-0112）。`symmetry` の操作で領域ごとの物質が移り合わなければエラー
#[tracing::instrument(name = "shielding", skip_all)]
pub fn shield_fuel_with(
    library: &Library,
    tracking: &Tracking,
    materials: &[MaterialInput],
    region_materials: &[String],
    symmetry: &Detected,
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
        // 燃料のリングごとの物質。リングの扇形は同じ物質で、リングどうしは別の物質でもよいが、核種の並びと温度は同じ
        // （等価原理を体積平均の組成で、UFG の範囲をリングごとの組成で解くため。ADR-0173）
        let sectors = pin.sectors();
        let ring_materials: Vec<usize> = (0..rings)
            .map(|j| region_index[range.start + j * sectors])
            .collect();
        let first = &materials[ring_materials[0]];
        let same_nuclides = |m: &MaterialInput| {
            m.temperature == first.temperature
                && m.composition
                    .iter()
                    .map(|(n, _)| n)
                    .eq(first.composition.iter().map(|(n, _)| n))
        };
        if (0..rings).any(|j| {
            let ring = range.start + j * sectors..range.start + (j + 1) * sectors;
            region_index[ring].iter().any(|&i| i != ring_materials[j])
        }) || !ring_materials.iter().all(|&i| same_nuclides(&materials[i]))
        {
            return Err(SelfShieldingError::MixedFuel(cell));
        }
        pins.push((cell, fuel_regions, ring_materials, rings));
    }
    let factors = tracing::debug_span!("dancoff")
        .in_scope(|| dancoff(tracking, &totals, &fuel, &options.dancoff))?;

    let mut shielded = FuelShielding {
        materials: dilute,
        region_material: region_index.clone(),
        pins: Vec::new(),
        lumps: Vec::new(),
        lattice: lattice.clone(),
    };
    // 格子の対称性で移り合うピンは同じ結果になるので、組の代表（最も小さいセル番号。セルの順で最初に現れる）
    // だけを解き、ほかは代表の物質を共有する。使う操作は、MOC と同じ確かめ方（トラックの像まで）で、領域ごとの入力の
    // 物質が移り合うもの（ADR-0112）。渡された操作が、この格子の物質で成り立つことを確かめる
    if !symmetry.fits(tracking) {
        return Err(SelfShieldingError::ForeignSymmetry);
    }
    for image in symmetry.images() {
        if let Some((r, &i)) = image
            .regions
            .iter()
            .enumerate()
            .find(|&(r, &i)| region_index[r] != region_index[i])
        {
            return Err(SelfShieldingError::Symmetry(SymmetryError::Broken {
                operation: image.name,
                reason: format!(
                    "region {r} has `{}`, but its image {i} has `{}`",
                    region_materials[r], region_materials[i]
                ),
            }));
        }
    }
    let orbits = symmetry.cell_representatives();
    // 解くピンを、下のループと同じ判定（同じ組の代表がもう解かれていれば使い回す）で先に決める。[ピン] の、解く
    // なら `solutions` の番号
    let mut planned: Vec<usize> = Vec::new();
    let plan: Vec<Option<usize>> = pins
        .iter()
        .map(|&(cell, ..)| {
            if planned.contains(&orbits[cell]) {
                None
            } else {
                planned.push(cell);
                Some(planned.len() - 1)
            }
        })
        .collect();
    // 解くピンは互いに独立なので並列に解く（各ピンの中でも衝突確率を並列に求めるが、ピンの減速計算は逐次なので、
    // ピンの間でも並列にしてコアを埋める）。ピンの計算は元の物質（`dilute`。領域の番号はこの範囲）だけを読む。
    // 結果の物質の追加と領域への割り当て、エラーは、下のループでピンの順に逐次に扱うので、逐次に解いたときと同じ
    let dilute_materials: &[Material] = &shielded.materials;
    let solving: Vec<(usize, &[usize], usize)> = pins
        .iter()
        .zip(&plan)
        .filter(|(_, p)| p.is_some())
        .map(|((cell, _, ring_materials, rings), _)| (*cell, ring_materials.as_slice(), *rings))
        .collect();
    // ピンの区間の親は自己遮蔽の区間にする（rayon のスレッドの今の区間は、盗んだ仕事の元のピンの区間でありうる）
    let shielding = tracing::Span::current();
    let mut solutions: Vec<Option<Result<PinSolution, SelfShieldingError>>> = solving
        .par_iter()
        .map(|&(cell, ring_materials, rings)| {
            let _pin = tracing::debug_span!(parent: &shielding, "pin", cell).entered();
            let factors = factors[cell]
                .as_ref()
                .ok_or(SelfShieldingError::FuelShape(cell))?;
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
            let material = pin_fuel(lattice, cell, materials, ring_materials);
            let carlvik = tracing::debug_span!("carlvik")
                .in_scope(|| shield_material(library, &material, &equivalence, options, cell))?;
            let pin_cell = PinRings {
                lattice,
                cell,
                rings,
                region_index: &region_index,
                void,
                collision: options.collision,
            };
            let ufg = shield_cell(
                library,
                &pin_cell,
                materials,
                dilute_materials,
                factors,
                &carlvik,
                |inputs| cell_range(library, inputs),
            )?;
            Ok(PinSolution { carlvik, ufg })
        })
        .map(Some)
        .collect();
    let mut solved: Vec<(usize, Solved)> = Vec::new();
    for ((cell, fuel_regions, ring_materials, rings), planned) in pins.into_iter().zip(plan) {
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
        let material = pin_fuel(lattice, cell, materials, &ring_materials);
        // 先に並列で解いた結果（同じ判定で解くと決めたピン）。エラーはここでピンの順に返す
        let PinSolution { carlvik, ufg } = planned
            .and_then(|k| solutions.get_mut(k).and_then(Option::take))
            .ok_or(SelfShieldingError::Unsolved(cell))??;
        let (effective, ring_effective) = match &ufg {
            Some(c) => (c.effective.clone(), c.center.clone()),
            None => (carlvik.clone(), vec![carlvik.clone(); rings]),
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
        let pin = PinShielding {
            cell,
            dancoff: factors.clone(),
            nuclides,
            effective,
            rings: ring_effective,
            equivalent: ufg.as_ref().map(|c| c.radius),
            sph: ufg.as_ref().map(|c| Sph {
                factors: c.sph.clone(),
                groups: c.groups.clone(),
                moderator_ring: c.start,
            }),
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
            // リングの物質ごとに、そのリングの組成と等価原理の値で作る（同じ物質のリングは 1 つを共有する）
            None => {
                let mut made: Vec<(usize, usize)> = Vec::new();
                let mut per_ring = Vec::with_capacity(rings);
                for &index in &ring_materials {
                    if let Some(&(_, m)) = made.iter().find(|(i, _)| *i == index) {
                        per_ring.push(m);
                        continue;
                    }
                    let ring = &materials[index];
                    let nuclides = ring
                        .composition
                        .iter()
                        .zip(&carlvik)
                        .map(|((name, density), s)| {
                            Ok((
                                shielded_nuclide(library, name, ring.temperature, s)?,
                                *density,
                            ))
                        })
                        .collect::<Result<Vec<_>, SelfShieldingError>>()?;
                    let mixed = Material::from_nuclides(&nuclides, library.groups.len()).map_err(
                        |source| SelfShieldingError::Material {
                            name: ring.name.clone(),
                            source,
                        },
                    )?;
                    shielded.materials.push(mixed);
                    let m = shielded.materials.len() - 1;
                    made.push((index, m));
                    per_ring.push(m);
                }
                Assigned::Fuel(per_ring)
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
    tracing::debug_span!("lumps").in_scope(|| {
        shield_lumps(
            library,
            tracking,
            materials,
            &region_index,
            (symmetry.images(), &orbits),
            options,
            &mut shielded,
        )
    })?;
    tracing::info!(
        pins = shielded.pins.len(),
        solved = solved.len(),
        "fuel pins and their symmetry representatives"
    );
    Ok(shielded)
}

/// 燃料でなく、燃料のピンの等価なセルでも置き換えていない領域のうち、分離共鳴を持つ核種を含む物質の領域
/// （支持格子の枠、案内管など）を、等価原理で遮蔽する。
///
/// 分離共鳴を持つ核種の組が同じ物質は、共鳴が同じエネルギーにあり互いに影を作る。そこで、セルの中で接する
/// （[`Cell::adjacent`](crate::geometry::Cell::adjacent)）組の同じ領域のひとまとまりを1つの黒体とし、組ごとに
/// その組の物質の領域（ピンの等価なセルで遮蔽した被覆管なども含む）を同時に黒体にして、黒体ごとの Dancoff 係数と
/// 平均弦長を求める（[`lump_dancoff`]。組の違う物質は影を作らない）。黒体の中の物質ごとの領域を塊とし、
/// 黒体の Dancoff 係数と平均弦長で遮蔽する。対称に移り合うセルの組の対応する黒体は1つにする。塊が棒（[`rod_rings`]）
/// なら、UFG の範囲の群は、塊の Dancoff 係数を合わせた等価なセルの減速計算の、棒のリングの平均にする（範囲は
/// [`material_range`]。棒の物質の核種で決める）。
fn shield_lumps(
    library: &Library,
    tracking: &Tracking,
    materials: &[MaterialInput],
    region_index: &[usize],
    (images, orbits): (&[Image], &[usize]),
    options: &ShieldingOptions,
    shielded: &mut FuelShielding,
) -> Result<(), SelfShieldingError> {
    let lattice = &tracking.lattice;
    let n_regions = region_index.len();
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
    let candidate: Vec<bool> = (0..n_regions)
        .map(|r| {
            shielded.region_material[r] == region_index[r]
                && !resonant_of(region_index[r]).is_empty()
        })
        .collect();
    if !candidate.contains(&true) {
        return Ok(());
    }
    // セルの中で接する、組の同じ候補の領域をまとめる。まとまりの代表は最も小さい領域の番号
    let mut parent: Vec<usize> = (0..n_regions).collect();
    let root = |parent: &mut Vec<usize>, mut r: usize| {
        while parent[r] != r {
            parent[r] = parent[parent[r]];
            r = parent[r];
        }
        r
    };
    for cell in 0..lattice.cells() {
        let first = lattice.cell_regions(cell).start;
        for (a, b) in lattice.cell_type(cell).adjacent() {
            let (a, b) = (first + a, first + b);
            if candidate[a]
                && candidate[b]
                && resonant_of(region_index[a]) == resonant_of(region_index[b])
            {
                let (ra, rb) = (root(&mut parent, a), root(&mut parent, b));
                parent[ra.max(rb)] = ra.min(rb);
            }
        }
    }
    // 使う対称操作で移したまとまりを同じ組にし、組を1つの黒体にする（ADR-0112）。まとまりの像は、まとまりの領域の
    // 像を含むまとまり（セルの中の接し方は対称操作で保たれるので、まとまりの領域はどれも同じまとまりに移る）。組の
    // 代表は、像のまとまりの代表の領域の最も小さいもの
    let mut lump_root = vec![usize::MAX; n_regions];
    for region in (0..n_regions).filter(|&r| candidate[r]) {
        let representative = root(&mut parent, region);
        let mut least = representative;
        for image in images {
            let map = &image.regions;
            let moved = root(&mut parent, map[region]);
            // 領域の像のまとまりが、代表の領域の像のまとまりと違えば、対称操作が接し方を保っていない
            if moved != root(&mut parent, map[representative]) || !candidate[map[region]] {
                return Err(SelfShieldingError::LumpImage(region));
            }
            least = least.min(moved);
        }
        lump_root[region] = least;
    }
    // 黒体は組の代表ごと（セルは代表の領域のセル）、塊は (黒体, 物質の番号) ごと
    let mut bodies: HashMap<usize, usize> = HashMap::new();
    // [黒体] の (セル, 核種の組)
    let mut body_of: Vec<(usize, &[&str])> = Vec::new();
    let mut body = vec![None; n_regions];
    let mut lumps: HashMap<(usize, usize), usize> = HashMap::new();
    // [塊] の (黒体, 物質の番号, 領域)
    let mut lump_of: Vec<(usize, usize, Vec<usize>)> = Vec::new();
    for region in (0..n_regions).filter(|&r| candidate[r]) {
        let key = lump_root[region];
        let b = *bodies.entry(key).or_insert_with(|| {
            body_of.push((lattice.region_cell(key), resonant_of(region_index[region])));
            body_of.len() - 1
        });
        body[region] = Some(b);
        let index = region_index[region];
        let l = *lumps.entry((b, index)).or_insert_with(|| {
            lump_of.push((b, index, Vec::new()));
            lump_of.len() - 1
        });
        lump_of[l].2.push(region);
    }
    // 影の計算の全断面積は、燃料のピンの Dancoff 係数と同じく燃料でない領域は無限希釈とし、燃料はピンの
    // SPH 因子を掛けない実効断面積（物理的な値。SPH 因子は MOC の反応率を合わせるためのもの）とする。
    // 対称に移り合うピンは同じ核種を共有するので、組ごとに1回求める
    let ng = library.groups.len();
    let mut totals: Vec<Vec<f64>> = region_index
        .iter()
        .map(|&i| shielded.materials[i].total.clone())
        .collect();
    let mut fuel_totals: HashMap<usize, Vec<f64>> = HashMap::new();
    for pin in &shielded.pins {
        let orbit = orbits[pin.cell];
        let total = match fuel_totals.entry(orbit) {
            Entry::Occupied(e) => e.into_mut(),
            Entry::Vacant(e) => {
                let fuel = Material::from_nuclides(&pin.nuclides, ng).map_err(|source| {
                    SelfShieldingError::Material {
                        name: format!("fuel in cell {}", pin.cell),
                        source,
                    }
                })?;
                e.insert(fuel.total)
            }
        };
        for region in lattice.cell_regions(pin.cell) {
            if is_fuel(&shielded.materials[region_index[region]]) {
                totals[region].clone_from(total);
            }
        }
    }
    // 核種の組ごとに、その組の物質の領域を黒体にし、その組の黒体の結果を求める。黒体の番号は組の中での番号に付け直す
    let mut found: Vec<Option<Lump>> = vec![None; body_of.len()];
    let mut sets: Vec<&[&str]> = body_of.iter().map(|&(_, set)| set).collect();
    sets.sort_unstable();
    sets.dedup();
    for set in sets {
        let group: Vec<usize> = (0..body_of.len())
            .filter(|&b| body_of[b].1 == set)
            .collect();
        let mut number = vec![None; body_of.len()];
        for (k, &b) in group.iter().enumerate() {
            number[b] = Some(k);
        }
        let black: Vec<bool> = region_index
            .iter()
            .map(|&i| resonant_of(i) == set)
            .collect();
        let local: Vec<Option<usize>> = body.iter().map(|b| b.and_then(|b| number[b])).collect();
        let results = lump_dancoff(
            tracking,
            &totals,
            &black,
            &local,
            group.len(),
            &options.dancoff,
        )?;
        for (b, result) in group.into_iter().zip(results) {
            found[b] = Some(result);
        }
    }
    // 塊の遮蔽は互いに独立なので並列に解く（ピンと同じ）。解くときに読む物質は、ここまでに足した物質（ピンの結果を
    // 含む）で、塊の結果の物質は下のループで塊の順に足すので、逐次に解いたときと同じ番号と値になる。エラーも塊の順に返す
    let parent = tracing::Span::current();
    let dilute: &[Material] = &shielded.materials;
    let solutions: Vec<Result<LumpSolution, SelfShieldingError>> = lump_of
        .par_iter()
        .map(|(b, index, lump_regions)| {
            let _lump = tracing::debug_span!(parent: &parent, "lump", body = b).entered();
            let cell = body_of[*b].0;
            let found = found[*b].as_ref().ok_or(SelfShieldingError::Body(*b))?;
            let equivalence = Equivalence {
                chord: found.chord,
                terms: found
                    .dancoff
                    .iter()
                    .map(|&c| carlvik(c))
                    .collect::<Result<_, _>>()?,
            };
            let material = &materials[*index];
            let tables = shield_material(library, material, &equivalence, options, cell)?;
            // 棒は、UFG の範囲を塊の Dancoff 係数を合わせた等価なセルの減速計算で遮蔽し、リングごとに SPH 因子を掛けた
            // 群定数を置く（燃料のピンと同じ。ADR-0142）
            let rod = match rod_rings(lattice, cell, lump_regions) {
                Some(rings) => {
                    let pin = PinRings {
                        lattice,
                        cell,
                        rings,
                        region_index,
                        void: materials.len(),
                        collision: options.collision,
                    };
                    shield_cell(
                        library,
                        &pin,
                        materials,
                        dilute,
                        &found.dancoff,
                        &tables,
                        |_| material_range(library, material),
                    )?
                    .map(|c| (c, rings))
                }
                None => None,
            };
            Ok(LumpSolution { tables, rod })
        })
        .collect();
    for ((b, index, lump_regions), solution) in lump_of.into_iter().zip(solutions) {
        let cell = body_of[b].0;
        let found = found[b].as_ref().ok_or(SelfShieldingError::Body(b))?;
        let material = &materials[index];
        let LumpSolution { tables, rod } = solution?;
        let (tables, equivalent, sph) = match rod {
            Some((c, rings)) => {
                // 棒のリングの物質。塊の領域のリングは、その領域のセル（対称操作の像）の中の番号から求める
                let mut ring_material = Vec::with_capacity(rings);
                for m in c.rings.into_iter().take(rings) {
                    let m = m.ok_or(SelfShieldingError::CellShape(cell))?;
                    shielded.materials.push(m);
                    ring_material.push(shielded.materials.len() - 1);
                }
                for &region in &lump_regions {
                    let shape = || SelfShieldingError::CellShape(lattice.region_cell(region));
                    let ring = LumpShielding::ring_of(lattice, region).ok_or_else(shape)?;
                    shielded.region_material[region] =
                        *ring_material.get(ring).ok_or_else(shape)?;
                }
                let factors = c.sph.into_iter().take(rings).collect();
                (c.effective, Some(c.radius), Some((factors, c.groups)))
            }
            None => {
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
                let mixed = Material::from_nuclides(&nuclides, ng).map_err(|source| {
                    SelfShieldingError::Material {
                        name: material.name.clone(),
                        source,
                    }
                })?;
                shielded.materials.push(mixed);
                let m = shielded.materials.len() - 1;
                for &region in &lump_regions {
                    shielded.region_material[region] = m;
                }
                (tables, None, None)
            }
        };
        shielded.lumps.push(LumpShielding {
            cell,
            material: material.name.clone(),
            regions: lump_regions,
            dancoff: found.dancoff.clone(),
            chord: found.chord,
            equivalent,
            sph,
            effective: tables,
        });
    }
    Ok(())
}

/// ピンの燃料の物質。リングがすべて同じ物質ならその物質、違えばリングの体積で個数密度を平均した物質（等価原理に使う。
/// 核種の並びと温度が同じことは確かめてある。ADR-0173）
fn pin_fuel<'a>(
    lattice: &Lattice,
    cell: usize,
    materials: &'a [MaterialInput],
    ring_materials: &[usize],
) -> std::borrow::Cow<'a, MaterialInput> {
    let first = &materials[ring_materials[0]];
    if ring_materials.iter().all(|&i| i == ring_materials[0]) {
        return std::borrow::Cow::Borrowed(first);
    }
    let range = lattice.cell_regions(cell);
    let volumes = lattice.volumes();
    // 燃料はピンセルのリングにだけ置ける（呼ぶ前に確かめてある）
    let sectors = lattice.cell_type(cell).pin().map_or(1, PinCell::sectors);
    let ring_volume = |j: usize| -> f64 {
        volumes[range.start + j * sectors..range.start + (j + 1) * sectors]
            .iter()
            .sum()
    };
    let weights: Vec<f64> = (0..ring_materials.len()).map(ring_volume).collect();
    let total: f64 = weights.iter().sum();
    let composition = first
        .composition
        .iter()
        .enumerate()
        .map(|(k, (name, _))| {
            let sum: f64 = ring_materials
                .iter()
                .zip(&weights)
                .map(|(&i, w)| w * materials[i].composition[k].1)
                .sum();
            (name.clone(), sum / total)
        })
        .collect();
    std::borrow::Cow::Owned(MaterialInput {
        name: format!("fuel of cell {cell}"),
        temperature: first.temperature,
        composition,
    })
}

/// 組の代表のピンの結果
struct Solved {
    pin: PinShielding,
    assigned: Assigned,
}

/// 代表のピンの等価原理と UFG の減速計算の結果
struct PinSolution {
    carlvik: Vec<Shielded>,
    ufg: Option<CellShielding>,
}

/// ピンのセルの領域に割り当てる物質の番号
enum Assigned {
    /// 等価なセルのリングごとの物質（真空は `None`）と、減速材のリングの始まり
    Rings {
        materials: Vec<Option<usize>>,
        start: usize,
    },
    /// 燃料のリングごとの物質（UFG の範囲がないとき。燃料の領域だけ）
    Fuel(Vec<usize>),
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
        Assigned::Fuel(rings) => {
            let sectors = fuel_regions.len() / rings.len().max(1);
            for (k, region) in fuel_regions.clone().enumerate() {
                region_material[region] = rings[k / sectors];
            }
        }
    }
    Ok(())
}

/// 等価なセルのリングごとの、格子の領域（ピンセル `geometry` の最初の領域が `first`）。リング `start` は減速材で、
/// セルのそれより外のリングをすべて含む（外周の枠は含めない。分離共鳴を持つ核種を含む枠は、塊として
/// [`shield_lumps`] で遮蔽する）。
fn ring_targets(geometry: &PinCell, first: usize, start: usize) -> Vec<Vec<usize>> {
    let (rings, sectors) = (geometry.radii().len(), geometry.sectors());
    let regions = |ring: usize| first + ring * sectors..first + (ring + 1) * sectors;
    (0..start)
        .map(|j| regions(j).collect())
        .chain([(start..=rings).flat_map(regions).collect()])
        .collect()
}

/// 中心の物質（燃料か、燃料でない棒の物質）のリングのあるピンのセルと、その領域の物質
struct PinRings<'a> {
    lattice: &'a Lattice,
    cell: usize,
    /// 中心の物質のリングの数
    rings: usize,
    /// [領域] の物質の番号
    region_index: &'a [usize],
    /// 真空の物質の番号
    void: usize,
    /// 減速計算の飛程の求め方
    collision: Collision,
}

/// 塊の遮蔽の結果: 等価原理の実効断面積と、棒なら等価なセルの減速計算の結果と棒のリングの数
struct LumpSolution {
    tables: Vec<Shielded>,
    rod: Option<(CellShielding, usize)>,
}

/// ピンのセルの UFG の減速計算の結果
struct CellShielding {
    /// 等価なセルの減速材のリングの外半径 (cm)
    radius: f64,
    /// 等価なセルの減速材までのリングごとの物質（真空は `None`）。格子の領域は [`ring_targets`]。枠のリングの
    /// 物質は格子に使わないので含めない（枠は塊として遮蔽する）
    rings: Vec<Option<Material>>,
    /// 減速材のリングの始まり
    start: usize,
    /// 燃料の核種の実効断面積。UFG の範囲の群は燃料のリングの値を中性子束と体積で平均し、範囲の外は等価原理の値
    effective: Vec<Shielded>,
    /// [中心のリング][核種] の実効断面積（UFG の範囲の群はリングの値、範囲の外は等価原理の値）
    center: Vec<Vec<Shielded>>,
    /// [リング][UFG の範囲の群] の SPH 因子
    sph: Vec<Vec<f64>>,
    /// `sph` の列の群
    groups: std::ops::Range<usize>,
}

/// 等価なセルの減速計算（[`solve_equivalent_cell`]）
struct EquivalentCell<'a> {
    /// 減速材のリングの外半径 (cm)
    radius: f64,
    cylinder: Cylinder,
    /// リングごとの物質（真空は `None`）。中心の物質のリング、減速材の前までの帯、減速材、枠があれば枠
    inputs: Vec<Option<&'a MaterialInput>>,
    /// 減速材のリングの始まり
    start: usize,
    /// 最後のリングが枠か
    frame: bool,
    /// 減速計算した群
    groups: std::ops::Range<usize>,
    spectrum: CellSpectrum,
}

/// セル `pin.cell` の内側の `pin.rings` のリング（燃料は核種の並びの同じ物質、棒は同じ物質）を中心に、群ごとの Dancoff 係数 `dancoff` を合わせた
/// 円筒の等価なセルを作り、減速計算を解く。範囲は `range_of`（等価なセルのリングごとの物質から決める）で、`None` なら
/// 解かずに `None`。リングの扇形はどれも同じ物質で、減速材（セルのいちばん外の物質）は真空でないこと。
fn solve_equivalent_cell<'a>(
    library: &Library,
    pin: &PinRings,
    materials: &'a [MaterialInput],
    dilute: &[Material],
    dancoff: &[f64],
    range_of: impl FnOnce(&[Option<&MaterialInput>]) -> Result<Option<[f64; 2]>, UfgError>,
) -> Result<Option<EquivalentCell<'a>>, SelfShieldingError> {
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
    // 外周の枠: (物質の番号, 面積)。等価なセルでは減速材のリングの外に同じ面積のリングとして置く（枠の影と、
    // 枠の物質が減速材の中性子束を変える効果を入れるため）
    let frame = geometry
        .frame()
        .map(|(region, area)| (pin.region_index[first + region], area));
    // 中心のリングは、それぞれのリングの物質（燃料のリングは核種の並びが同じ。ADR-0173）
    let inputs: Vec<Option<&MaterialInput>> = ring_index[..pin.rings]
        .iter()
        .map(|&i| input(i))
        .chain(ring_index[pin.rings..start].iter().map(|&i| input(i)))
        .chain([input(moderator)])
        .chain(frame.map(|(i, _)| input(i)))
        .collect();
    let Some(range) = range_of(&inputs)? else {
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
        let outer = frame.map(|(i, area)| OuterRing {
            area,
            total: dilute[i].total[g],
        });
        sum += equivalent_radius(&inner, &totals, outer, c.max(MIN_DANCOFF))?;
    }
    let radius = sum / (bottom - top) as f64;
    let cylinder = Cylinder::new(
        radii[..start]
            .iter()
            .copied()
            .chain([radius])
            .chain(frame.map(|(_, area)| outer_radius(radius, area)))
            .collect(),
    )?;
    let spectrum = slowing_down_with(library, &cylinder, &inputs, range, pin.collision)?;
    Ok(Some(EquivalentCell {
        radius,
        cylinder,
        inputs,
        start,
        frame: frame.is_some(),
        groups: top..bottom,
        spectrum,
    }))
}

/// 中心のリングの群定数の平均: `base` の核種ごとの値の、UFG の範囲の群を、中心の `rings` のリングの値を中性子束と体積で
/// 平均したものに置き換える（範囲の外の群は `base` のまま）。
fn center_average(base: &[Shielded], equivalent: &EquivalentCell, rings: usize) -> Vec<Shielded> {
    let EquivalentCell {
        cylinder,
        groups,
        spectrum,
        ..
    } = equivalent;
    let volumes = cylinder.volumes();
    base.iter()
        .enumerate()
        .map(|(k, c)| {
            let mut s = c.clone();
            for g in groups.clone() {
                let weight = |j: usize| spectrum.flux[j][g - groups.start] * volumes[j];
                let total_weight: f64 = (0..rings).map(weight).sum();
                let average = |x: fn(&Shielded) -> &Vec<f64>| -> f64 {
                    (0..rings)
                        .map(|j| weight(j) * x(&spectrum.nuclides[j][k])[g - groups.start])
                        .sum::<f64>()
                        / total_weight
                };
                s.total[g] = average(|v| &v.total);
                s.elastic[g] = average(|v| &v.elastic);
                s.fission[g] = average(|v| &v.fission);
                s.nu_fission[g] = average(|v| &v.nu_fission);
                s.capture[g] = average(|v| &v.capture);
            }
            s
        })
        .collect()
}

/// 塊の領域 `regions` が、セル `cell` のピンセルの内側のリングを（扇形すべてで）ちょうど埋めていれば、そのリングの数
/// （塊は棒）。
fn rod_rings(lattice: &Lattice, cell: usize, regions: &[usize]) -> Option<usize> {
    let geometry = lattice.cell_type(cell).pin()?;
    if geometry.is_grid() {
        return None;
    }
    let first = lattice.cell_regions(cell).start;
    let sectors = geometry.sectors();
    let mut here: Vec<usize> = regions
        .iter()
        .copied()
        .filter(|&r| lattice.region_cell(r) == cell)
        .collect();
    here.sort_unstable();
    let rings = here.len() / sectors;
    (rings > 0
        && rings <= geometry.radii().len()
        && here.len().is_multiple_of(sectors)
        && here.iter().copied().eq(first..first + rings * sectors))
    .then_some(rings)
}

/// ピンのセルを、Dancoff 係数を合わせた円筒の等価なセルの UFG の減速計算で自己遮蔽する。燃料、燃料と減速材の間の
/// リング（ギャップ、被覆管）、減速材（セルの外側の物質。等価なセルでは1つのリング）のそれぞれの核種に、そのリングの
/// 中性子束で縮約した群定数を与える。UFG の範囲の外の群は、中心の物質（燃料か棒）は等価原理の値（`carlvik`）、それ以外は無限希釈。
/// `range_of` が中心の物質の UFG の範囲を返さなければ（分離共鳴を持つ核種がなければ）`None`。
fn shield_cell(
    library: &Library,
    pin: &PinRings,
    materials: &[MaterialInput],
    dilute: &[Material],
    dancoff: &[f64],
    carlvik: &[Shielded],
    range_of: impl FnOnce(&[Option<&MaterialInput>]) -> Result<Option<[f64; 2]>, UfgError>,
) -> Result<Option<CellShielding>, SelfShieldingError> {
    let cell = pin.cell;
    let Some(equivalent) =
        solve_equivalent_cell(library, pin, materials, dilute, dancoff, range_of)?
    else {
        return Ok(None);
    };
    let EquivalentCell {
        radius,
        cylinder,
        inputs,
        start,
        frame,
        groups,
        spectrum,
    } = &equivalent;
    let (top, bottom) = (groups.start, groups.end);
    let bounds = library.groups.bounds();

    let ng = library.groups.len();
    let infinite = *library.sigma0.first().ok_or(LookupError::NoSigma0)?;
    let mut rings = Vec::with_capacity(inputs.len());
    let mut center = Vec::with_capacity(pin.rings);
    for (j, input) in inputs.iter().enumerate() {
        let Some(m) = input else {
            rings.push(None);
            continue;
        };
        let mut fulls = Vec::with_capacity(m.composition.len());
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
                    full.capture[g] = ufg.capture[g - top];
                }
                let nuclide = shielded_nuclide(library, name, m.temperature, &full)?;
                fulls.push(full);
                Ok((nuclide, *density))
            })
            .collect::<Result<Vec<_>, SelfShieldingError>>()?;
        if j < pin.rings {
            center.push(fulls);
        }
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
    // 枠のリングの断面積は格子に使わないので、因子を掛けない
    let mut fixed = vec![false; multigroup.len()];
    if *frame && let Some(last) = fixed.last_mut() {
        *last = true;
    }
    let sph = tracing::debug_span!("sph").in_scope(|| {
        sph_factors(
            cylinder,
            &multigroup,
            &fixed,
            top..bottom,
            &above,
            &spectrum.flux,
        )
    })?;
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

    // 枠のリングの物質は格子に使わない（枠は塊として遮蔽する）
    if *frame {
        rings.pop();
    }

    // 中心の物質（燃料か棒）の平均: そのリングの中性子束と体積の重み
    let effective = center_average(carlvik, &equivalent, pin.rings);
    Ok(Some(CellShielding {
        radius: *radius,
        rings,
        start: *start,
        effective,
        center,
        sph,
        groups: top..bottom,
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
    // 個数密度が 0 の核種は何にも寄与しないので、無限希釈のまま（反復で変わらないので 1 度だけ引く。燃焼では輸送できる核種の
    // 多くが 0 で入る）
    let dilute = material
        .composition
        .iter()
        .map(|(name, density)| {
            (*density == 0.0)
                .then(|| library.shielded(name, t, &vec![infinite; ng]))
                .transpose()
        })
        .collect::<Result<Vec<_>, _>>()?;
    for _ in 0..options.max_iterations {
        let tables = material
            .composition
            .iter()
            .enumerate()
            .map(|(k, (name, density))| {
                if let Some(table) = &dilute[k] {
                    return Ok(table.clone());
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
            [&s.total, &s.elastic, &s.fission, &s.nu_fission, &s.capture]
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
/// 弾性散乱の変化の分だけ変える（ADR-0076）。
pub fn shielded_nuclide(
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
    // ライブラリの最初の σ0 は無限希釈とみなす値（xlc_format::library）。燃焼でできたばかりの微量の核種はそれより大きな σ0 を
    // 持つ（数密度が非正規数なら無限大）が、σt ≪ σ0 で群定数も重みの中性子束ももう σ0 によらないので、無限希釈の値に置く
    // （ADR-0173。重みの中性子束の差は σt / σ0 の水準）。NaN は min で消えるので、比べて置き換える（NaN は表の検査で止める）
    let infinite = *library.sigma0.first().ok_or(LookupError::NoSigma0)?;
    let dilute = |s: Vec<f64>| -> Vec<f64> {
        s.into_iter()
            .map(|x| if x > infinite { infinite } else { x })
            .collect()
    };
    let (sigma0_1, sigma0_2) = (dilute(sigma0(0)), dilute(sigma0(1)));
    let tables = [
        library.shielded(name, temperature, &sigma0_1)?,
        library.shielded(name, temperature, &sigma0_2)?,
    ];
    let sigma0 = [sigma0_1, sigma0_2];
    let groups = background.len();
    // 分離共鳴を持たない核種の表は、滑らかな断面積の群の中の形でしか σ0 によらず、2つの表がほぼ等しいので、
    // 重みの σp は効かない
    let potential = library
        .nuclides
        .iter()
        .find(|n| n.name == name)
        .ok_or_else(|| LookupError::Nuclide(name.into()))?
        .resolved
        .as_ref()
        .map_or(0.0, |r| r.potential);
    // 群ごとの重み β_i φ_i
    let weights: Vec<[f64; 2]> = (0..groups)
        .map(|g| {
            std::array::from_fn(|i| {
                let s = sigma0[i][g];
                equivalence.terms[g][i].beta * (potential + s) / (tables[i].total[g] + s)
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
        capture: combine(|x| &x.capture),
    })
}

#[cfg(test)]
mod tests;
