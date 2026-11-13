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
//! 等価原理では、燃料の脱出確率を P(x) = Σ_i β_i x / (x + α_i)（x = Σt l̄）で近似し、核種 k のバックグラウンド断面積を
//! σ0,i = (Σ_{j≠k} N_j σt,j + α_i / l̄) / N_k とする。実効断面積は、2つの σ0 でのライブラリの値を、
//! NR 近似の群平均の中性子束 φ_i = (σp + σ0,i) / (σt(σ0,i) + σ0,i) と β_i で重みづけして合わせる。σp は核種の
//! ポテンシャル散乱で、NR 近似の中性子束 (σp + σ0) w / (σt + σ0) の線源に入る（共鳴の間で中性子束が w になる）。
//! ライブラリの重みが w/(σt + σ0) の群（Bondarenko）では、この式は NR 近似の範囲で厳密に成り立つ。分離共鳴の
//! 範囲の群は flux calculator の中性子束で重みづけしてあり、この式は近似になる（ADR-0089）。
//! σp は分離共鳴の範囲の s 波の散乱半径の値をすべての群に使う（非分離共鳴の範囲の p 波などは入れない）。

use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, OnceLock, PoisonError};
use std::time::{Duration, Instant};

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
    CellSpectrum, CellTangent, CensusCall, Collision, MemoryError, NuclideTransfer, Scattering,
    ShareCellCache, TangentCost, UfgError, cell_range, held_cells, hold_total, kernel_cap,
    material_range, slowing_down_tangent, slowing_down_transfers, slowing_down_with,
};

/// 物質の核種ごとの、UFG の範囲の群の弾性散乱の移行（[`UfgTransfer::RingFlux`]）
type Transfers = Arc<Vec<NuclideTransfer>>;

/// ライブラリの行の和が無限希釈の弾性散乱より小さくてよい丸めの大きさ（行の和に対する比）。depletion.msgpack の UFG の範囲では、その差は
/// 7e−16 以下
const NON_ELASTIC_ROUNDING: f64 = 1.0e-12;

/// リングごとの UFG の範囲の 1 群の捕獲断面積をログに書く核種（自己遮蔽の条件を比べる。#0286）。Gd は核分裂しないので、捕獲が
/// 吸収
const RING_CAPTURES: [&str; 5] = ["U-238", "Pu-239", "Pu-240", "Gd-155", "Gd-157"];

/// UFG の減速計算の等価なセルの、弦の求積の小区間あたりの点の数の既定（ADR-0252）。[`crate::cp::POINTS`]（16）は等価原理の
/// セルなど UFG の外の衝突確率も使う既定で、UFG の厳密な経路として入力（`ufg-cp-points = 16`）で選べる
pub const UFG_CP_POINTS: usize = 6;

/// 格子の燃料の自己遮蔽の反復の条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShieldingOptions {
    /// 燃料の各核種の全断面積の相対変化の上限
    pub tolerance: f64,
    pub max_iterations: usize,
    pub dancoff: DancoffOptions,
    /// 燃料のピンと棒の UFG の減速計算の飛程の求め方
    pub collision: Collision,
    /// 燃料のピンと棒の UFG の減速計算の等価なセルの、弦の求積の小区間あたりの点の数（[`Cylinder::with_points`]。既定は
    /// [`UFG_CP_POINTS`]。減らすと衝突確率を求める時間が点の数に比例して減り、求積の誤差が増える。#0286、ADR-0252）
    pub cp_points: usize,
    /// 燃料のピンと棒の UFG の減速計算の格子を、リングの物質の巨視的な断面積をこの相対差で再現できる点に間引く（`None` は
    /// 間引かない。ADR-0185）
    pub thinning: Option<f64>,
    /// 自由ガスの散乱源の区間をスイープの間で持つ上限を決める全体の量 (B)。減速計算の数えた作業の配列を除いた残りで持つ（`None` は
    /// 自己遮蔽のたびに使用可能メモリから求める。ADR-0213）
    pub hold_memory: Option<u64>,
    /// 燃料のピンと棒の UFG の減速計算の散乱核の標的の扱い（ADR-0203）
    pub scattering: Scattering,
    /// 等価なセルのリングの物質の、UFG の範囲の群の散乱行列の移行の形（ADR-0207）
    pub transfer: UfgTransfer,
}

impl Default for ShieldingOptions {
    fn default() -> Self {
        Self {
            tolerance: 1.0e-6,
            max_iterations: 100,
            dancoff: DancoffOptions::default(),
            collision: Collision::default(),
            cp_points: UFG_CP_POINTS,
            thinning: None,
            hold_memory: None,
            scattering: Scattering::default(),
            transfer: UfgTransfer::default(),
        }
    }
}

/// 等価なセルのリングの物質の、UFG の範囲の群の散乱行列の群の間の移行の形（ADR-0207。既定は ADR-0209）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UfgTransfer {
    /// ライブラリの無限希釈の形のまま、行の和を自己遮蔽した弾性散乱に合わせる（ADR-0076）
    InfiniteDilution,
    /// 弾性散乱の移行を、UFG の減速計算でリングの中性子束と同じ散乱核で数える
    #[default]
    RingFlux,
}

/// 比で求める状態（ADR-0250）の、比が決めない量の作り方。どれも段 1 で `every-state` と比べて決める候補で、`Default` は ADR-0250 が
/// 決めるまで使う基準の値
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct RatioOptions {
    /// 等価なセルのリングの UFG の範囲の移行
    pub transfer: RatioTransfer,
    /// 等価なセルの燃料でないリング（被覆管、減速材）の群定数
    pub rings: RatioRings,
    /// リングの物質に掛ける SPH 因子
    pub sph: RatioSph,
    /// 燃料のリングの UFG の範囲の群定数
    pub fuel: RatioFuel,
}

/// 比で求める状態の、燃料のリングの UFG の範囲の群定数。`Equivalence` のほかは、比の残る差をピンの平均の変化とリングの間の分布に
/// 分ける診断（#0285）で、状態の UFG の減速計算を解く（時間は状態ごとに解くのと同じ）
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RatioFuel {
    /// 基準のリングの実効断面積に、状態と基準の等価原理の値の比を掛ける（ADR-0250）
    #[default]
    Equivalence,
    /// 基準のリングの実効断面積に、状態と基準の UFG のピンの平均の比を掛ける（リングの間の分布は基準のまま）
    StateAverage,
    /// 状態の UFG のリングの実効断面積と移行
    State,
    /// 基準のリングの実効断面積に、状態と基準のそれぞれを、同じ等価なセルの減速計算を静止した標的（[`Scattering::Static`]。移行は
    /// 無限希釈）で解いた値の比を掛ける。比はリングごとなので、リングの間の分布の変化も入る。静止した標的と自由ガスの差（絶対値の
    /// 差。ADR-0209）は比で打ち消す
    StaticRatio,
}

/// 比で求める状態の、等価なセルのリングの UFG の範囲の移行（ADR-0250 の「比が決めない量」）
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RatioTransfer {
    /// 基準の移行の行をそのまま使う
    #[default]
    Reference,
    /// 基準の移行の行を、入射の群ごとに状態と基準の弾性散乱の比で縮める（行の和が状態の弾性散乱になる）
    ElasticRatio,
}

/// 比で求める状態の、等価なセルの燃料でないリング（被覆管、減速材）の群定数（ADR-0250 の「比が決めない量」）。UFG の範囲の外は
/// どれも状態の温度の無限希釈
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RatioRings {
    /// UFG の範囲は基準のリングの実効断面積（状態の数密度で混ぜる）
    #[default]
    Reference,
    /// UFG の範囲は基準のリングの実効断面積に、状態と基準の温度の無限希釈の値の比を掛けたもの（ADR-0201 の案 M の形）
    ShieldingFactor,
    /// UFG の範囲も状態の温度の無限希釈（移行の行を持たない）
    InfiniteDilution,
    /// UFG の範囲は基準のリングの実効断面積に、状態と基準を静止した標的で解いた減速計算の値の比を掛けたもの（[`RatioFuel::StaticRatio`]）
    StaticRatio,
    /// 状態の UFG のリングの実効断面積と移行（診断。状態の UFG の減速計算を解く）
    State,
    /// UFG の範囲は状態の UFG のリングの実効断面積で、移行は基準の行（[`RatioTransfer`] で縮める。診断）
    StateSections,
}

/// 比で求める状態の SPH 因子（ADR-0250 の「SPH 因子は基準のものを使う」の近似を分けて測る）
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RatioSph {
    /// 基準の SPH 因子
    #[default]
    Reference,
    /// 状態の UFG の減速計算を解いて求めた SPH 因子（診断。時間は状態ごとに解き直すのと同じ）
    State,
    /// 基準の SPH 因子に、状態と基準を静止した標的で解いた減速計算の SPH 因子の比を掛けたもの（[`RatioFuel::StaticRatio`]）
    StaticRatio,
}

/// 比で求める状態（ADR-0250）の基準: 今の UFG で解いた自己遮蔽 `shielding` と、それに渡した物質 `materials`（状態の物質と同じ名前を
/// 同じ順に持つ）。`light` はこの基準の静止した標的の減速計算を持ち回る置き場で、同じ基準の状態の間で共有し、別の基準とは共有しない
#[derive(Debug, Clone, Copy)]
pub struct RatioAnchor<'a> {
    pub shielding: &'a FuelShielding,
    pub materials: &'a [MaterialInput],
    pub options: RatioOptions,
    pub light: &'a AnchorLight,
}

/// 比の基準 1 つの、静止した標的の減速計算（[`RatioFuel::StaticRatio`] などの比の分母）の結果のセルごとの置き場。基準の物質と
/// Dancoff 係数は状態によらないので、最初の状態で解いた結果を同じ基準のほかの状態で使う（解き直しても結果はビットで同じ）
#[derive(Default)]
pub struct AnchorLight {
    /// 基準の物質の無限希釈の群定数（[`mixed`]）
    dilute: OnceLock<Vec<Material>>,
    /// （セル、結果）
    cells: Mutex<Vec<(usize, Arc<CellShielding>)>>,
}

impl std::fmt::Debug for AnchorLight {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let cells: Vec<usize> = self.lock().iter().map(|(c, _)| *c).collect();
        f.debug_struct("AnchorLight")
            .field("cells", &cells)
            .finish()
    }
}

impl AnchorLight {
    /// 基準の物質 `materials` の無限希釈の群定数（最初に求めたものを持つ）
    fn dilute(
        &self,
        library: &Library,
        materials: &[MaterialInput],
    ) -> Result<&[Material], SelfShieldingError> {
        if let Some(d) = self.dilute.get() {
            return Ok(d);
        }
        let d = mixed(library, materials)?;
        Ok(self.dilute.get_or_init(|| d))
    }

    /// セル `cell` の結果。持っていなければ `solve` で解いて持つ（ピンは並列に解くが、1 回の自己遮蔽で同じセルは 1 度だけ解く
    /// ので、解く間は置き場を占めない）
    fn cell(
        &self,
        cell: usize,
        solve: impl FnOnce() -> Result<CellShielding, SelfShieldingError>,
    ) -> Result<Arc<CellShielding>, SelfShieldingError> {
        let held = |cells: &[(usize, Arc<CellShielding>)]| {
            cells
                .iter()
                .find(|(c, _)| *c == cell)
                .map(|(_, s)| s.clone())
        };
        if let Some(s) = held(&self.lock()) {
            return Ok(s);
        }
        let solved = Arc::new(solve()?);
        let mut cells = self.lock();
        if let Some(s) = held(&cells) {
            return Ok(s);
        }
        cells.push((cell, solved.clone()));
        Ok(solved)
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, Vec<(usize, Arc<CellShielding>)>> {
        // 置き場の中身は push だけで変わり、途中の状態を残さないので、ほかのスレッドが panic した後も使える
        self.cells.lock().unwrap_or_else(PoisonError::into_inner)
    }
}

/// 燃料のピンの UFG の範囲の解き方
#[derive(Debug, Clone, Copy)]
enum PinMethod<'a> {
    /// 等価なセルの減速計算（`Some` なら、その核種の個数密度の向きの 1 次の変化も。[`shield_fuel_tangent`]）
    Ufg(Option<&'a [String]>),
    /// 基準の UFG の実効断面積に、状態と基準の等価原理の比を掛ける（[`shield_fuel_ratio`]）
    Ratio(RatioAnchor<'a>),
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
    /// 解くと決めたピンの結果がないか、2 度解いた（内部の不整合）
    #[error("cell {0}: internal error, the pin was not shielded exactly once")]
    Unsolved(usize),
    /// 1 次の変化と一緒に解いた減速計算が、向きのない減速計算と違う（内部の不整合。[`shield_fuel_tangent`]）
    #[error(
        "cell {0}: internal error, the slowing down solved with its tangent differs from the plain one"
    )]
    TangentSpectrum(usize),
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
    Memory(#[from] MemoryError),
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
    #[error(
        "{name}: transfer rows of {given} (groups, rows or row entries) for {groups} (internal error)"
    )]
    TransferGroups {
        name: String,
        given: usize,
        groups: usize,
    },
    #[error(
        "{name}: the scattering row of group {group} sums to {row} b, below the infinitely dilute elastic {elastic} b"
    )]
    ElasticAboveRow {
        name: String,
        group: usize,
        row: f64,
        elastic: f64,
    },
    #[error(
        "{name}: a scattering moment of group {group} to group {to} exceeds its P0 element (no element-wise P_l/P0)"
    )]
    MomentAboveP0 {
        name: String,
        group: usize,
        to: usize,
    },
    #[error("the lattice is not the one that was self-shielded")]
    ForeignLattice,
    /// 比の基準の物質の名前と順が状態と違う（ADR-0250）
    #[error("the ratio anchor has other materials than the state")]
    AnchorMaterials,
    /// 比の基準に、このセルのピンか、等価なセルのリングの実効断面積が無いか、リングの核種が状態と違う（ADR-0250）
    #[error(
        "cell {0}: the ratio anchor does not shield the pin with the rings and nuclides of the state"
    )]
    AnchorPin(usize),
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
    /// [領域] の、物質の核種ごとの物理的な実効断面積（[`FuelShielding::region_cross_sections`]）
    region_effective: Vec<Option<Arc<Vec<Shielded>>>>,
    /// [領域] の、物質の核種ごとの UFG の範囲の群の弾性散乱の移行（[`FuelShielding::region_transfers`]）
    region_transfers: Vec<Option<Transfers>>,
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
/// 共有する（`cell` 以外は代表の値。核種ごとの値は複製せず、代表と同じものを指す）。
#[derive(Debug, Clone, PartialEq)]
pub struct PinShielding {
    pub cell: usize,
    /// 群ごとの Dancoff 係数
    pub dancoff: Vec<f64>,
    /// 燃料の核種の自己遮蔽した群定数 (barn) と個数密度 (1/(barn cm))。散乱行列は ADR-0076 の形（[`UfgTransfer::RingFlux`] でも、
    /// リングの移行の行は平均せず、格子の物質だけがその行を使う）。燃料のリングの物質が違うときは、個数密度は
    /// リングの体積で平均したもの（等価原理はその組成で解く。ADR-0173）
    pub nuclides: Arc<Vec<(Nuclide, f64)>>,
    /// 核種ごとの実効断面積（`nuclides` の順。弾性散乱を含む）。UFG の範囲の群は、燃料のリングの値を中性子束と
    /// 体積で平均したもの。`nuclides` とともに SPH 因子を掛けていない物理的な値（参照解と比べる用）で、MOC の
    /// 中性子束と組み合わせて反応率を求めるには、領域の物質（SPH 因子を掛けたもの）を使う
    pub effective: Arc<Vec<Shielded>>,
    /// [燃料のリング（内側から）][核種（`nuclides` の順）] の実効断面積。UFG の範囲の群はそのリングの減速計算の値、範囲の外は
    /// 等価原理の値。`effective` と同じく SPH 因子を掛けない物理的な値で、物理的な中性子束（`FuelShielding::physical_flux`）と
    /// 掛けてリングの反応率にする（ADR-0173）。領域の値（[`FuelShielding::region_cross_sections`]）と同じものを共有する
    pub rings: Vec<Arc<Vec<Shielded>>>,
    /// UFG の減速計算に使った等価なセルの減速材のリングの外半径 (cm)。支持格子のあるセルでは、その外に枠と同じ面積の
    /// 支持格子のリングがある。燃料に分離共鳴を持つ核種がなければ `None`
    pub equivalent: Option<f64>,
    /// SPH 因子（[`crate::sph`]）。燃料に分離共鳴を持つ核種がなく UFG の範囲がなければ `None`
    pub sph: Option<Sph>,
    /// UFG の減速計算の、核種の個数密度の向きごとの 1 次の変化（[`shield_fuel_tangent`] で頼んだときだけ）
    pub tangent: Option<Arc<PinTangent>>,
}

/// 燃料のピンの UFG の減速計算の、核種の個数密度の向き（等価なセルのどのリングでも N に比例）ごとの 1 次の変化
/// （[`shield_fuel_tangent`]。ADR-0247 の段 I2-2 の測定）。中心のリング（燃料）だけを持つ
#[derive(Debug, Clone, PartialEq)]
pub struct PinTangent {
    /// 向きの核種（[向き] の順）
    pub nuclides: Vec<String>,
    /// UFG の範囲の群（ライブラリの群の番号）と、減速計算の格子点の数
    pub groups: std::ops::Range<usize>,
    pub points: usize,
    /// 向きのない減速計算と、1 次の変化と一緒に解いた減速計算の経過時間 (s)
    pub plain_seconds: f64,
    pub tangent_seconds: f64,
    pub cost: TangentCost,
    /// [中心のリング][群] の中性子束と、[向き][中心のリング][群] のその 1 次の変化
    pub flux: Vec<Vec<f64>>,
    pub flux_change: Vec<Vec<Vec<f64>>>,
    /// [中心のリング][反応][群] の巨視的断面積 Σ_k N_k σ_k (1/cm) と、[向き][中心のリング][反応][群] の暗黙の分
    /// Σ_k N_k δσ_k（個数密度を (1 + δ) 倍したときの δ についての微分。陽の分 Σ_k δN_k σ_k を含まない）。反応は
    /// [`PinTangent::REACTIONS`] の順
    pub macroscopic: Vec<[Vec<f64>; 4]>,
    pub implicit: Vec<Vec<[Vec<f64>; 4]>>,
}

impl PinTangent {
    /// [`PinTangent::macroscopic`] の反応
    pub const REACTIONS: [&'static str; 4] = ["total", "capture", "nu-fission", "elastic"];

    fn new(
        names: &[String],
        inputs: &[Option<&MaterialInput>],
        rings: usize,
        spectrum: &CellSpectrum,
        tangent: CellTangent,
        plain_seconds: f64,
        tangent_seconds: f64,
    ) -> Self {
        fn reactions(s: &Shielded) -> [&Vec<f64>; 4] {
            [&s.total, &s.capture, &s.nu_fission, &s.elastic]
        }
        // [リング][反応][群] の Σ_k N_k x_k（x は核種ごとの群定数）
        let macroscopic = |nuclides: &[Vec<Shielded>]| -> Vec<[Vec<f64>; 4]> {
            (0..rings)
                .map(|j| {
                    let composition = inputs[j].map_or(&[][..], |m| m.composition.as_slice());
                    std::array::from_fn(|x| {
                        let ng = spectrum.groups.len();
                        (0..ng)
                            .map(|g| {
                                composition
                                    .iter()
                                    .zip(&nuclides[j])
                                    .map(|((_, n), s)| n * reactions(s)[x][g])
                                    .sum()
                            })
                            .collect()
                    })
                })
                .collect()
        };
        Self {
            nuclides: names.to_vec(),
            groups: spectrum.groups.clone(),
            points: spectrum.points,
            plain_seconds,
            tangent_seconds,
            cost: tangent.cost.clone(),
            flux: spectrum.flux[..rings].to_vec(),
            flux_change: tangent.flux.iter().map(|f| f[..rings].to_vec()).collect(),
            macroscopic: macroscopic(&spectrum.nuclides),
            implicit: tangent.nuclides.iter().map(|n| macroscopic(n)).collect(),
        }
    }
}

/// 等価なセルのリングごと・UFG の範囲の群ごとの SPH 因子。その群の、リングの物質を置いた格子の領域の物質に掛けてある。
#[derive(Debug, Clone, PartialEq)]
pub struct Sph {
    /// [リング][群 − `groups.start`]。支持格子のあるセルの最後は支持格子のリングの値で、格子には使わないので 1
    pub factors: Vec<Vec<f64>>,
    /// `factors` の列の群（ライブラリの群の番号）
    pub groups: std::ops::Range<usize>,
    /// 等価なセルの減速材のリングの番号（リングを格子の領域に対応させる）
    moderator_ring: usize,
}

impl FuelShielding {
    /// 領域 `region` の物質を遮蔽した、核種ごとの実効断面積（物質の組成の順。SPH 因子を掛けない物理的な値で、物理的な中性子束
    /// [`FuelShielding::physical_flux`] と掛けて反応率にする）。燃料のピンの等価なセルのリング（燃料、IFBA の被覆などの燃料と減速材の
    /// 間のリング、減速材）と遮蔽した塊の領域は遮蔽した値、ほかの物質の領域は MOC と同じ無限希釈の値。真空は `None`
    pub fn region_cross_sections(&self, region: usize) -> Option<&[Shielded]> {
        self.region_effective
            .get(region)?
            .as_deref()
            .map(Vec::as_slice)
    }

    /// 領域 `region` の物質の核種ごとの、UFG の範囲の群の弾性散乱の移行（[`UfgTransfer::RingFlux`] で遮蔽した等価なセルのリングの
    /// 領域だけ。核種の順は [`FuelShielding::region_cross_sections`] と同じ）。領域の物質は、これと実効断面積から
    /// [`shielded_nuclide_with`] で作った群定数を混ぜたもの
    pub fn region_transfers(&self, region: usize) -> Option<&[NuclideTransfer]> {
        self.region_transfers
            .get(region)?
            .as_deref()
            .map(Vec::as_slice)
    }

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
            // 因子はリングごと。支持格子のあるセルでは最後に、格子に使わない支持格子のリングの行（因子は 1）が 1 つ多い
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
    /// SPH 因子のリングごとに、そのリングの物質を置いた格子 `lattice` の領域（支持格子のリングは格子に使わないので含めない）。
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

impl FuelShielding {
    /// MOC の中性子束 `flux`（[領域][群]。生成率で規格化したもの）から、燃料の領域の [`RING_CAPTURES`] の核種の捕獲と核分裂の反応率
    /// Σ_g N σ_g φ_g V（物理的な中性子束と実効断面積）を、領域の物質の名前の `@` の前（燃焼の領域の元の物質。リングの位置）ごとに足して、
    /// target `ring_rates` の DEBUG のログに出す（比で求める状態のリングの反応率を状態ごとに解いたものと比べる。ADR-0250）。
    /// `region_materials` は [領域] の物質の名前、`materials` はその物質
    pub fn log_ring_rates(
        &self,
        lattice: &Lattice,
        flux: &[Vec<f64>],
        region_materials: &[String],
        materials: &[MaterialInput],
    ) -> Result<(), SelfShieldingError> {
        if !tracing::enabled!(target: "ring_rates", tracing::Level::DEBUG) {
            return Ok(());
        }
        let physical = self.physical_flux(lattice, flux)?;
        let volumes = lattice.volumes();
        // (物質, 核種) ごとの [捕獲, 核分裂]（出す順は最初に現れた順）
        let mut rates: Vec<(&str, &str, [f64; 2])> = Vec::new();
        for (region, name) in region_materials.iter().enumerate() {
            let is_fuel = self.materials[self.region_material[region]]
                .nu_fission
                .iter()
                .any(|&v| v > 0.0);
            let (Some(sections), Some(material)) = (
                self.region_cross_sections(region),
                materials.iter().find(|m| &m.name == name),
            ) else {
                continue;
            };
            if !is_fuel {
                continue;
            }
            let base = name.split('@').next().unwrap_or(name);
            for ((nuclide, density), s) in material.composition.iter().zip(sections) {
                if !RING_CAPTURES.contains(&nuclide.as_str()) {
                    continue;
                }
                let rate = |x: &[f64]| -> f64 {
                    density
                        * volumes[region]
                        * x.iter()
                            .zip(&physical[region])
                            .map(|(a, b)| a * b)
                            .sum::<f64>()
                };
                let value = [rate(&s.capture), rate(&s.fission)];
                match rates
                    .iter_mut()
                    .find(|(m, n, _)| *m == base && *n == nuclide)
                {
                    Some((_, _, sum)) => {
                        sum[0] += value[0];
                        sum[1] += value[1];
                    }
                    None => rates.push((base, nuclide, value)),
                }
            }
        }
        for (material, nuclide, [capture, fission]) in rates {
            tracing::debug!(
                target: "ring_rates",
                material,
                nuclide,
                capture,
                fission,
                "ring reaction rates per neutron produced"
            );
        }
        Ok(())
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
pub fn shield_fuel_with(
    library: &Library,
    tracking: &Tracking,
    materials: &[MaterialInput],
    region_materials: &[String],
    symmetry: &Detected,
    options: &ShieldingOptions,
) -> Result<FuelShielding, SelfShieldingError> {
    shield_fuel_reusing(
        library,
        &mut ShieldingReuse::new(tracking),
        materials,
        region_materials,
        symmetry,
        options,
    )
}

/// 1 つのトラッキングで自己遮蔽を繰り返して解くとき（燃焼の輸送ごと、分岐の状態ごと）に持ち越す、燃料の組成によらない前処理。
/// 燃料のピンの Dancoff 係数は、燃料を黒体にし、燃料でない物質の全断面積だけから求めるので、直前に解いたときと入力（領域ごとの
/// 物質、燃料でない物質の全断面積、`DancoffOptions`）が同じなら、その結果を使う（ADR-0194）。呼び出し側が複数（状態と数密度の
/// 微分の摂動の組成。[`ShieldingReuse::sharing`]）なら、その数だけ、直前に解いた違う入力の結果を持つ
#[derive(Debug)]
pub struct ShieldingReuse<'a> {
    tracking: &'a Tracking,
    /// 解いた入力と結果（最後に使った順。`callers` 個まで）
    dancoff: Vec<DancoffInputs>,
    callers: usize,
    /// 燃料のピンの Dancoff 係数を解いた回数
    dancoff_solves: usize,
    /// 等価なセルの減速計算（燃料のピン、棒の塊、比の静止した標的）を解いた回数
    slowing_downs: AtomicUsize,
    /// 移行の計数の自由ガスの核の割合の区間の表。格子、中性子束、個数密度に依らず、輸送とピンの間で使える（ADR-0213）
    shares: ShareCellCache,
    /// 呼び出し側が最後に渡した燃焼点（[`ShieldingReuse::at_point`]）
    point: Option<usize>,
}

/// 燃料のピンの Dancoff 係数を求めた入力と結果
#[derive(Debug)]
struct DancoffInputs {
    /// [領域] の物質の番号
    region_index: Vec<usize>,
    /// [物質] の全断面積（領域に置いた燃料でない物質だけ。ほかは `None`）
    totals: Vec<Option<Vec<f64>>>,
    options: DancoffOptions,
    /// [セル] の群ごとの Dancoff 係数
    factors: Arc<Vec<Option<Vec<f64>>>>,
}

impl<'a> ShieldingReuse<'a> {
    /// トラッキング `tracking` で解く自己遮蔽の前処理（まだ何も持たない）
    pub fn new(tracking: &'a Tracking) -> Self {
        Self::sharing(tracking, NonZeroUsize::MIN)
    }

    /// トラッキング `tracking` で解く自己遮蔽の、`callers` 人の呼び出し側（状態と、数密度の微分の摂動の組成）が共有する前処理。
    /// 摂動が燃料でない物質を変えても（燃料でない燃焼の領域の核種の微分）、ほかの呼び出し側の直前の結果を追い出さない
    pub fn sharing(tracking: &'a Tracking, callers: NonZeroUsize) -> Self {
        Self {
            tracking,
            dancoff: Vec::new(),
            callers: callers.get(),
            dancoff_solves: 0,
            slowing_downs: AtomicUsize::new(0),
            shares: ShareCellCache::default(),
            point: None,
        }
    }

    /// これからの自己遮蔽が燃焼点 `point`（呼び出し側の番号）のものとする。直前に渡した燃焼点と違えば、自由ガスの核の区間と D の値の
    /// 表（ADR-0248）を次の燃焼点に進め、直前の燃焼点で使わなかった格子点を捨てる。表は同じ燃焼点の中では上限に届かない限り格子点を
    /// 捨てないので、呼ばなければ、すべての自己遮蔽を 1 つの燃焼点として扱う
    pub fn at_point(&mut self, point: usize) {
        if self.point.is_some_and(|p| p != point) {
            self.shares.kernels().next_point();
            self.shares.next_point_census();
        }
        self.point = Some(point);
    }

    /// 前処理のトラッキング
    pub fn tracking(&self) -> &'a Tracking {
        self.tracking
    }

    /// 燃料のピンの Dancoff 係数を解いた回数
    pub fn dancoff_solves(&self) -> usize {
        self.dancoff_solves
    }

    /// 等価なセルの減速計算（燃料のピン、棒の塊、比の静止した標的）を解いた回数
    pub fn slowing_downs(&self) -> usize {
        self.slowing_downs.load(Ordering::Relaxed)
    }

    /// 物質 `dilute`（無限希釈）を領域に `region_index` で置いた格子の、燃料のピンの Dancoff 係数（[セル]）。`fuel` は
    /// [領域] の燃料か
    fn pin_dancoff(
        &mut self,
        dilute: &[Material],
        region_index: &[usize],
        fuel: &[bool],
        options: &DancoffOptions,
    ) -> Result<Arc<Vec<Option<Vec<f64>>>>, SelfShieldingError> {
        // 領域に置いた燃料でない物質の全断面積だけを比べる（燃料は黒体で、どの領域にも置かない物質とともに値を読まない）
        let mut totals: Vec<Option<Vec<f64>>> = vec![None; dilute.len()];
        for (&i, &f) in region_index.iter().zip(fuel) {
            if !f && totals[i].is_none() {
                totals[i] = Some(dilute[i].total.clone());
            }
        }
        if let Some(k) = self.dancoff.iter().position(|d| {
            d.region_index == region_index && d.totals == totals && d.options == *options
        }) {
            tracing::debug!("pin Dancoff factors reused");
            // 使った入力を前に移す（追い出すのは最も長く使っていない入力）
            let hit = self.dancoff.remove(k);
            let factors = hit.factors.clone();
            self.dancoff.insert(0, hit);
            return Ok(factors);
        }
        let region_totals: Vec<Vec<f64>> = region_index
            .iter()
            .map(|&i| totals[i].clone().unwrap_or_default())
            .collect();
        let factors = Arc::new(
            tracing::info_span!("dancoff")
                .in_scope(|| dancoff(self.tracking, &region_totals, fuel, options))?,
        );
        self.dancoff_solves += 1;
        self.dancoff.truncate(self.callers - 1);
        self.dancoff.insert(
            0,
            DancoffInputs {
                region_index: region_index.to_vec(),
                totals,
                options: *options,
                factors: factors.clone(),
            },
        );
        Ok(factors)
    }
}

/// `shield_fuel_with` と同じ。トラッキングは `reuse` のもので、燃料の組成によらない前処理を `reuse` から持ち越し、`reuse` に
/// 残す（[`ShieldingReuse`]）。結果は `shield_fuel_with` とビットで同じ
pub fn shield_fuel_reusing(
    library: &Library,
    reuse: &mut ShieldingReuse,
    materials: &[MaterialInput],
    region_materials: &[String],
    symmetry: &Detected,
    options: &ShieldingOptions,
) -> Result<FuelShielding, SelfShieldingError> {
    shield_fuel_tangent(
        library,
        reuse,
        materials,
        region_materials,
        symmetry,
        options,
        None,
    )
}

/// [`shield_fuel_reusing`] の自己遮蔽で、`tangent` を与えれば、燃料のピンの UFG の減速計算を、その核種の個数密度の向き（等価なセルの
/// どのリングでも N に比例）ごとの 1 次の変化と一緒にも解き（[`slowing_down_tangent`]。ADR-0247 の案 I2 の段 I2-2）、
/// [`PinShielding::tangent`] に書く。自己遮蔽の結果は `tangent` によらずビットで同じ（向きのない減速計算も解いて比べ、違えば止める）
#[tracing::instrument(name = "shielding", skip_all)]
pub fn shield_fuel_tangent(
    library: &Library,
    reuse: &mut ShieldingReuse,
    materials: &[MaterialInput],
    region_materials: &[String],
    symmetry: &Detected,
    options: &ShieldingOptions,
    tangent: Option<&[String]>,
) -> Result<FuelShielding, SelfShieldingError> {
    shield_fuel_by(
        library,
        reuse,
        materials,
        region_materials,
        symmetry,
        options,
        PinMethod::Ufg(tangent),
    )
}

/// [`shield_fuel_reusing`] の自己遮蔽で、燃料のピンの UFG の範囲の群を減速計算で解かず、基準 `anchor`（同じ格子を今の UFG で解いた
/// 自己遮蔽）の等価なセルのリングの実効断面積に、状態と基準の等価原理の値の比 σ_eq(状態) / σ_eq(基準) を掛けて求める（ADR-0250）。
/// σ_eq は範囲の外の群と同じ等価原理（状態と基準のそれぞれの Dancoff 係数と組成と温度）で、反応ごと・群ごとに比を取る（基準の値が
/// 0 なら比は 1）。範囲の外の群、燃料でない物質、構造材の塊は状態の物質で解く（[`shield_fuel_reusing`] と同じ）。比が決めない量
/// （燃料でないリングの群定数、移行、SPH 因子）は `anchor.options` で選ぶ。基準の物質が状態と同じなら、結果は基準とビットで同じ
#[tracing::instrument(name = "ratio-shielding", skip_all)]
pub fn shield_fuel_ratio(
    library: &Library,
    reuse: &mut ShieldingReuse,
    materials: &[MaterialInput],
    region_materials: &[String],
    symmetry: &Detected,
    options: &ShieldingOptions,
    anchor: RatioAnchor,
) -> Result<FuelShielding, SelfShieldingError> {
    if anchor.materials.len() != materials.len()
        || anchor
            .materials
            .iter()
            .zip(materials)
            .any(|(a, m)| a.name != m.name)
    {
        return Err(SelfShieldingError::AnchorMaterials);
    }
    if anchor.shielding.lattice != reuse.tracking.lattice {
        return Err(SelfShieldingError::ForeignLattice);
    }
    shield_fuel_by(
        library,
        reuse,
        materials,
        region_materials,
        symmetry,
        options,
        PinMethod::Ratio(anchor),
    )
}

/// [`shield_fuel_tangent`] と [`shield_fuel_ratio`] の自己遮蔽。燃料のピンの UFG の範囲を `method` で解く
fn shield_fuel_by(
    library: &Library,
    reuse: &mut ShieldingReuse,
    materials: &[MaterialInput],
    region_materials: &[String],
    symmetry: &Detected,
    options: &ShieldingOptions,
    method: PinMethod,
) -> Result<FuelShielding, SelfShieldingError> {
    let tangent = match method {
        PinMethod::Ufg(tangent) => tangent,
        PinMethod::Ratio(_) => None,
    };
    let tracking = reuse.tracking;
    let lattice = &tracking.lattice;
    let regions = lattice.volumes().len();
    if region_materials.len() != regions {
        return Err(SelfShieldingError::Regions {
            given: region_materials.len(),
            regions,
        });
    }
    let dilute = tracing::info_span!("mix").in_scope(|| {
        materials
            .iter()
            .map(|m| {
                Material::mix(library, &m.composition, m.temperature).map_err(|source| {
                    SelfShieldingError::Material {
                        name: m.name.clone(),
                        source,
                    }
                })
            })
            .collect::<Result<Vec<_>, _>>()
    })?;
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
    let factors = reuse.pin_dancoff(&dilute, &region_index, &fuel, &options.dancoff)?;
    let shares = &reuse.shares;
    let solves = &reuse.slowing_downs;
    // 自己遮蔽のたびに求める（使用可能メモリは燃焼の間に変わる）。`total` は前の自己遮蔽から残った自由ガスの核の表（ADR-0248）を
    // 除いて新たに確保できる量
    let held = shares.kernels().start();
    let total = hold_total(options.hold_memory, held)?;
    tracing::debug!(bytes = total, "UFG held memory");

    let mut shielded = FuelShielding {
        materials: dilute,
        region_material: region_index.clone(),
        pins: Vec::new(),
        lumps: Vec::new(),
        lattice: lattice.clone(),
        region_effective: vec![None; regions],
        region_transfers: vec![None; regions],
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
    // 解くピンを、下のループと同じ判定（同じ組の代表がもう解かれていれば再利用する）で先に決める。[ピン] の、解く
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
    // 解くピンは互いに独立なので並列に解く（各ピンの中でも飛程の表、衝突確率、縮約を並列に求めるが、ピンの減速計算は
    // 逐次なので、ピンの間でも並列にしてコアを埋める）。ピンの計算は元の物質（`dilute`。領域の番号はこの範囲）だけを読む。
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
    // `broadcast` のピンはスレッドごとに 1 本なので、同時に解くピンはスレッドの数と解くピンの数の小さい方。自由ガスの核の区間と
    // D の値の表（リングの中性子束の移行で、ピンと塊、自己遮蔽の間で共有する）は、同時に解く減速計算の 1 つと同じ量を持ち、
    // ピンと塊はその残りを分ける（ADR-0248）。上限は残った表を含む全体（`total` + 残り）から決め、ピンと塊には表がまだ増えうる量を
    // 除いて渡す（新たに確保するのは `total` まで）
    let threads = rayon::current_num_threads();
    let at_once = threads.min(solving.len());
    let growth = shares.kernels().limit(match options.transfer {
        UfgTransfer::RingFlux if at_once > 0 => {
            kernel_cap(total, held, at_once, shares.kernels().demand())
        }
        // 表を使わない移行と、解くピンの無い自己遮蔽では表を増やさない（残りは塊が引く）
        _ => held,
    });
    let hold = held_cells(total - growth, at_once);
    let solve = |pin: &(usize, &[usize], usize)| -> Result<PinSolution, SelfShieldingError> {
        let &(cell, ring_materials, rings) = pin;
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
            cp_points: options.cp_points,
            thinning: options.thinning,
            scattering: options.scattering,
            transfer: options.transfer,
            shares,
            hold,
            tangent,
            solves,
        };
        let ufg = match method {
            PinMethod::Ufg(_) => shield_cell(
                library,
                &pin_cell,
                materials,
                dilute_materials,
                factors,
                &carlvik,
                |inputs| cell_range(library, inputs),
            )?,
            PinMethod::Ratio(anchor) => tracing::debug_span!("ratio").in_scope(|| {
                ratio_cell(
                    library,
                    &pin_cell,
                    materials,
                    dilute_materials,
                    factors,
                    &carlvik,
                    ring_materials,
                    options,
                    &anchor,
                )
            })?,
        };
        Ok(PinSolution { carlvik, ufg })
    };
    // ピンは `par_iter` の仕事にせず、プールの各スレッドが 1 本ずつ番号の順に取る。`broadcast` の仕事はそのスレッドだけが
    // 実行するので、ピンの中の並列を待つスレッドが別のピンを丸ごと盗んで、待っていたピンの完了を遅らせることがない。
    // ピンを取り終えたスレッドは、残りのピンの中の並列を手伝う（ADR-0198）。`broadcast` はプールのすべてのスレッドが
    // 仕事を終えるまで返らないので、プールの外（ほかの並列の仕事が走っていないところ）から呼ぶ
    let next = AtomicUsize::new(0);
    // 2 度解いたピンのセル（なければ usize::MAX）
    let twice = AtomicUsize::new(usize::MAX);
    let slots: Vec<OnceLock<Result<PinSolution, SelfShieldingError>>> =
        solving.iter().map(|_| OnceLock::new()).collect();
    rayon::broadcast(|_| {
        loop {
            let k = next.fetch_add(1, Ordering::Relaxed);
            let Some(pin) = solving.get(k) else { break };
            // 番号 k を取ったのはこのスレッドだけなので、`set` は失敗しない。失敗すれば番号の取り方の誤りで、解を捨てずに止める
            if slots[k].set(solve(pin)).is_err() {
                twice.store(pin.0, Ordering::Relaxed);
            }
        }
    });
    if let cell @ 0..usize::MAX = twice.into_inner() {
        return Err(SelfShieldingError::Unsolved(cell));
    }
    let mut solutions: Vec<Option<Result<PinSolution, SelfShieldingError>>> =
        slots.into_iter().map(OnceLock::into_inner).collect();
    let mut solved: Vec<(usize, Solved)> = Vec::new();
    // 代表のピンのエネルギーの格子点の数（解いた順）と、移り合うピンへの複製の時間
    let mut points: Vec<usize> = Vec::new();
    // 解いたピンの 1 回の量の最大 (B)。次の自己遮蔽の自由ガスの核の表の上限に使う（ADR-0248）
    let mut demand = 0;
    let mut replicate = Duration::ZERO;
    let serial = Instant::now();
    for ((cell, fuel_regions, ring_materials, rings), planned) in pins.into_iter().zip(plan) {
        let factors = factors[cell]
            .as_ref()
            .ok_or(SelfShieldingError::FuelShape(cell))?;
        if let Some((_, s)) = solved.iter().find(|(c, _)| *c == orbits[cell]) {
            let copying = Instant::now();
            assign(
                lattice,
                cell,
                &s.assigned,
                &fuel_regions,
                &mut shielded.region_material,
                &mut shielded.region_effective,
                &mut shielded.region_transfers,
            )?;
            shielded.pins.push(PinShielding {
                cell,
                ..s.pin.clone()
            });
            replicate += copying.elapsed();
            continue;
        }
        let material = pin_fuel(lattice, cell, materials, &ring_materials);
        // 先に並列で解いた結果（同じ判定で解くと決めたピン）。エラーはここでピンの順に返す
        let PinSolution { carlvik, mut ufg } = planned
            .and_then(|k| solutions.get_mut(k).and_then(Option::take))
            .ok_or(SelfShieldingError::Unsolved(cell))??;
        points.extend(ufg.as_ref().map(|c| c.points));
        demand = demand.max(ufg.as_ref().map_or(0, |c| c.demand));
        // 範囲の外の燃料のリングは、どれも等価原理の同じ微視的断面積（組成の並びはリングで同じ）
        let micro = Arc::new(carlvik);
        let (effective, ring_effective) = match &mut ufg {
            Some(c) => (
                Arc::new(std::mem::take(&mut c.effective)),
                c.ring_nuclides[..rings]
                    .iter()
                    .map(|n| n.clone().ok_or(SelfShieldingError::CellShape(cell)))
                    .collect::<Result<Vec<_>, _>>()?,
            ),
            None => (micro.clone(), vec![micro.clone(); rings]),
        };
        let nuclides = material
            .composition
            .iter()
            .zip(effective.iter())
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
            nuclides: Arc::new(nuclides),
            effective,
            rings: ring_effective,
            equivalent: ufg.as_ref().map(|c| c.radius),
            sph: ufg.as_ref().map(|c| Sph {
                factors: c.sph.clone(),
                groups: c.groups.clone(),
                moderator_ring: c.start,
            }),
            tangent: ufg.as_ref().and_then(|c| c.tangent.clone()),
        };
        let assigned = match ufg {
            // リングごとの物質（真空のリングは真空のまま）
            Some(c) => Assigned::Rings {
                materials: c
                    .rings
                    .into_iter()
                    .zip(c.ring_nuclides)
                    .zip(c.ring_transfers)
                    .map(|((m, nuclides), transfers)| {
                        m.zip(nuclides).map(|(m, nuclides)| {
                            shielded.materials.push(m);
                            (shielded.materials.len() - 1, nuclides, transfers)
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
                        per_ring.push((m, micro.clone()));
                        continue;
                    }
                    let ring = &materials[index];
                    let nuclides = ring
                        .composition
                        .iter()
                        .zip(micro.iter())
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
                    per_ring.push((m, micro.clone()));
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
            &mut shielded.region_effective,
            &mut shielded.region_transfers,
        )?;
        shielded.pins.push(pin.clone());
        solved.push((cell, Solved { pin, assigned }));
    }
    let serial = serial.elapsed().as_secs_f64();
    tracing::info!(
        replicate = replicate.as_secs_f64(),
        solved = serial - replicate.as_secs_f64(),
        "serial loop over the pins (s): replicate = the copies to the symmetric pins, solved = the rest (the materials of the solved pins)"
    );
    tracing::info_span!("lumps").in_scope(|| {
        shield_lumps(
            library,
            tracking,
            materials,
            &region_index,
            (symmetry.images(), &orbits),
            options,
            (shares, solves, total - growth),
            &mut shielded,
        )
    })?;
    // 遮蔽しなかった物質の領域は、MOC の物質（`Material::mix`）と同じ無限希釈の値（σ0 の表の最初）。物質ごとに 1 度だけ作る
    let dilute_span = tracing::info_span!("dilute").entered();
    let mut dilute_nuclides: Vec<Option<Arc<Vec<Shielded>>>> = vec![None; materials.len()];
    for (region, &index) in region_index.iter().enumerate() {
        if index == void || shielded.region_effective[region].is_some() {
            continue;
        }
        if dilute_nuclides[index].is_none() {
            let m = &materials[index];
            let tables = m
                .composition
                .iter()
                .map(|(name, _)| Ok(library.lookup(name, m.temperature)?.1.shielded[0].clone()))
                .collect::<Result<Vec<_>, SelfShieldingError>>()?;
            dilute_nuclides[index] = Some(Arc::new(tables));
        }
        shielded.region_effective[region] = dilute_nuclides[index].clone();
    }
    drop(dilute_span);
    // 自由ガスの核の評価のうち、格子点と核が同じでピンの間（distinct）、直前の自己遮蔽の間（new）、燃焼点の中のそれより前の
    // 自己遮蔽の間（earlier）で重なる分と、核種ごとの核で数えた区間の数。数えるのは ufg の DEBUG のときで、ここの DEBUG と別に
    // 選べるので、出さないときも閉じて鍵を減らす
    let CensusCall {
        total,
        distinct,
        new,
        earlier,
        bytes,
        keys,
        nuclides,
    } = shares.census();
    if tracing::enabled!(tracing::Level::DEBUG) {
        tracing::debug!(
            total,
            distinct,
            new,
            earlier,
            distinct_bytes = bytes[0],
            earlier_bytes = bytes[1],
            point_bytes = bytes[2],
            keys,
            "UFG free-gas evaluations: total = sum over the pins, distinct = without the points shared between pins, new = also without the points of the previous shielding, earlier = also without the points of the earlier shieldings of the burnup point; bytes = the table size of those points"
        );
        for (kelvin, [total, distinct, earlier, point]) in nuclides {
            tracing::debug!(
                kelvin,
                total,
                distinct,
                earlier,
                point,
                "UFG free-gas cells by nuclide kernel: total = sum over the pins, distinct = without the cells shared between pins, earlier = also without the cells of the earlier shieldings of the burnup point, point = all the cells of the burnup point so far"
            );
        }
    }
    // 自由ガスの核の D の値の表（ADR-0248）: 値を使った核種と格子点の組、作った組、上限で入れなかった組の数と、持つ大きさ、上限、
    // 解いたピンの 1 回の量の最大（次の自己遮蔽の上限に使う）
    shares.kernels().set_demand(demand);
    let [reused, built, full] = shares.kernels().counts();
    tracing::info!(
        reused,
        built,
        full,
        bytes = shares.kernels().bytes(),
        cap = shares.kernels().cap(),
        demand,
        "UFG free-gas kernel points"
    );
    tracing::info!(
        pins = shielded.pins.len(),
        solved = solved.len(),
        points_min = ?points.iter().min(),
        points_max = ?points.iter().max(),
        points_sum = points.iter().sum::<usize>(),
        "fuel pins and their symmetry representatives"
    );
    Ok(shielded)
}

/// 同じ格子の 2 つの自己遮蔽 `nominal` と `perturbed`（数密度の微分の摂動の組成。ADR-0221）の、燃料のピンの核種ごとの実効断面積
/// （ピンの平均）の相対差 |σ(摂動) − σ(公称)|/σ(公称) を、核種ごと（ピンと群の最大と、その群）と群ごと（核種とピンの最大）にログに
/// 出す（摂動で変わる量とエネルギーの範囲。`bounds` は群の境界 (eV)）。群ごとの最大は、摂動した核種 `changed` とほかの核種に分ける。
/// 断面積は全断面積と捕獲で、公称が 0 の値は比べない
pub fn log_effective_change(
    bounds: &[f64],
    changed: &str,
    nominal: &FuelShielding,
    perturbed: &FuelShielding,
) {
    if !tracing::enabled!(tracing::Level::DEBUG) {
        return;
    }
    let relative = |a: f64, b: f64| if a > 0.0 { (b - a).abs() / a } else { 0.0 };
    // [核種] の（名前、[全, 捕獲] の最大、その群）と [群] の [摂動した核種, ほかの核種][全, 捕獲] の最大
    let mut by_nuclide: Vec<(String, [f64; 2], [usize; 2])> = Vec::new();
    let mut by_group: Vec<[[f64; 2]; 2]> = Vec::new();
    // 比べられない組（ピンの数か並び、核種の並びの違う摂動）は、名前と群を取り違えたログを出さずに、比べないことを出す
    if nominal.pins.len() != perturbed.pins.len() {
        tracing::warn!(
            nominal = nominal.pins.len(),
            perturbed = perturbed.pins.len(),
            "the perturbed self-shielding has another pin count; not compared"
        );
        return;
    }
    for (p, q) in nominal.pins.iter().zip(&perturbed.pins) {
        let names = |s: &PinShielding| -> Vec<String> {
            s.nuclides.iter().map(|(n, _)| n.name.clone()).collect()
        };
        if p.cell != q.cell
            || names(p) != names(q)
            || p.effective.len() != q.effective.len()
            || p.effective.len() != p.nuclides.len()
            || p.effective
                .iter()
                .zip(q.effective.iter())
                .any(|(a, b)| a.total.len() != b.total.len() || a.capture.len() != b.capture.len())
        {
            tracing::warn!(
                nominal = p.cell,
                perturbed = q.cell,
                "the perturbed self-shielding has another pin or nuclide order; not compared"
            );
            return;
        }
        for ((nuclide, _), (s, t)) in p
            .nuclides
            .iter()
            .zip(p.effective.iter().zip(q.effective.iter()))
        {
            let k = match by_nuclide.iter().position(|(n, ..)| *n == nuclide.name) {
                Some(k) => k,
                None => {
                    by_nuclide.push((nuclide.name.clone(), [0.0; 2], [0; 2]));
                    by_nuclide.len() - 1
                }
            };
            let pairs = [(&s.total, &t.total), (&s.capture, &t.capture)];
            for (x, (a, b)) in pairs.into_iter().enumerate() {
                for (g, (&a, &b)) in a.iter().zip(b).enumerate() {
                    let d = relative(a, b);
                    if by_group.len() <= g {
                        by_group.resize(g + 1, [[0.0; 2]; 2]);
                    }
                    let side = usize::from(nuclide.name != changed);
                    by_group[g][side][x] = by_group[g][side][x].max(d);
                    if d > by_nuclide[k].1[x] {
                        by_nuclide[k].1[x] = d;
                        by_nuclide[k].2[x] = g;
                    }
                }
            }
        }
    }
    for (name, worst, at) in &by_nuclide {
        tracing::debug!(
            name = name.as_str(),
            total = worst[0],
            total_group = at[0],
            capture = worst[1],
            capture_group = at[1],
            "effective cross section change by the perturbation (largest relative over the pins and groups)"
        );
    }
    for (g, worst) in by_group.iter().enumerate() {
        tracing::debug!(
            group = g,
            upper = bounds.get(g).copied(),
            lower = bounds.get(g + 1).copied(),
            total = worst[0][0],
            capture = worst[0][1],
            others_total = worst[1][0],
            others_capture = worst[1][1],
            "effective cross section change by the perturbation in a group (largest relative over the pins; the perturbed nuclide and the others)"
        );
    }
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
#[allow(clippy::too_many_arguments)]
fn shield_lumps(
    library: &Library,
    tracking: &Tracking,
    materials: &[MaterialInput],
    region_index: &[usize],
    (images, orbits): (&[Image], &[usize]),
    options: &ShieldingOptions,
    (shares, solves, total): (&ShareCellCache, &AtomicUsize, u64),
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
    // 塊の遮蔽は互いに独立なので並列に解く（塊は 6 個ほどで、ピンのようにスレッドごとに取る形にはしていない。#0171）。解くときに読む物質は、ここまでに足した物質（ピンの結果を
    // 含む）で、塊の結果の物質は下のループで塊の順に足すので、逐次に解いたときと同じ番号と値になる。エラーも塊の順に返す
    let parent = tracing::Span::current();
    let dilute: &[Material] = &shielded.materials;
    // 塊は `par_iter` なので、減速計算を待つスレッドがほかの塊を盗んで入れ子に解きうる。同時に解きうる減速計算は棒の塊の数まで
    let rods = lump_of
        .iter()
        .filter(|(b, _, lump_regions)| rod_rings(lattice, body_of[*b].0, lump_regions).is_some())
        .count();
    // `total` は、自己遮蔽の全体の量から自由ガスの核の区間と D の値の表がまだ増えうる量（ピンの前に決めたもの。ADR-0248）を除いた量
    let hold = held_cells(total, rods);
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
                        cp_points: options.cp_points,
                        thinning: options.thinning,
                        scattering: options.scattering,
                        transfer: options.transfer,
                        shares,
                        hold,
                        tangent: None,
                        solves,
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
                let ring_iter = c
                    .rings
                    .into_iter()
                    .zip(c.ring_nuclides)
                    .zip(c.ring_transfers);
                for ((m, nuclides), transfers) in ring_iter.take(rings) {
                    let shape = || SelfShieldingError::CellShape(cell);
                    let (m, nuclides) = (m.ok_or_else(shape)?, nuclides.ok_or_else(shape)?);
                    shielded.materials.push(m);
                    ring_material.push((shielded.materials.len() - 1, nuclides, transfers));
                }
                for &region in &lump_regions {
                    let shape = || SelfShieldingError::CellShape(lattice.region_cell(region));
                    let ring = LumpShielding::ring_of(lattice, region).ok_or_else(shape)?;
                    let (m, nuclides, transfers) = ring_material.get(ring).ok_or_else(shape)?;
                    shielded.region_material[region] = *m;
                    shielded.region_effective[region] = Some(nuclides.clone());
                    shielded.region_transfers[region] = transfers.clone();
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
                let nuclides = Arc::new(tables.clone());
                for &region in &lump_regions {
                    shielded.region_material[region] = m;
                    shielded.region_effective[region] = Some(nuclides.clone());
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

/// 等価なセルのリングの物質の番号と、その物質の核種ごとの実効断面積と移行
type RingAssignment = (usize, Arc<Vec<Shielded>>, Option<Transfers>);

/// ピンのセルの領域に割り当てる物質の番号と、その物質の核種ごとの実効断面積
enum Assigned {
    /// 等価なセルのリングごとの物質（真空は `None`）と、減速材のリングの始まり
    Rings {
        materials: Vec<Option<RingAssignment>>,
        start: usize,
    },
    /// 燃料のリングごとの物質（UFG の範囲がないとき。燃料の領域だけ）
    Fuel(Vec<(usize, Arc<Vec<Shielded>>)>),
}

/// ピンのセル `cell` の領域に、`assigned` の物質の番号を割り当てる。
fn assign(
    lattice: &Lattice,
    cell: usize,
    assigned: &Assigned,
    fuel_regions: &std::ops::Range<usize>,
    region_material: &mut [usize],
    region_effective: &mut [Option<Arc<Vec<Shielded>>>],
    region_transfers: &mut [Option<Transfers>],
) -> Result<(), SelfShieldingError> {
    match assigned {
        Assigned::Rings { materials, start } => {
            let pin = lattice
                .cell_type(cell)
                .pin()
                .ok_or(SelfShieldingError::FuelShape(cell))?;
            let first = lattice.cell_regions(cell).start;
            for (targets, material) in ring_targets(pin, first, *start).iter().zip(materials) {
                if let Some((m, nuclides, transfers)) = material {
                    for &region in targets {
                        region_material[region] = *m;
                        region_effective[region] = Some(nuclides.clone());
                        region_transfers[region] = transfers.clone();
                    }
                }
            }
        }
        Assigned::Fuel(rings) => {
            let sectors = fuel_regions.len() / rings.len().max(1);
            for (k, region) in fuel_regions.clone().enumerate() {
                let (m, nuclides) = &rings[k / sectors];
                region_material[region] = *m;
                region_effective[region] = Some(nuclides.clone());
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
    /// 等価なセルの弦の求積の小区間あたりの点の数
    cp_points: usize,
    /// 減速計算の格子の間引きの許容差
    thinning: Option<f64>,
    /// 減速計算の散乱核の標的の扱い
    scattering: Scattering,
    /// リングの物質の UFG の範囲の群の散乱行列の移行の形
    transfer: UfgTransfer,
    /// 移行の計数の自由ガスの核の割合の区間の表と、1 回の減速計算の作業の配列とスイープの間で持つ区間を合わせた量（[`held_cells`]。
    /// ADR-0213）
    shares: &'a ShareCellCache,
    hold: usize,
    /// 減速計算の 1 次の変化を求める核種（[`shield_fuel_tangent`]）
    tangent: Option<&'a [String]>,
    /// 等価なセルの減速計算を解いた回数（[`ShieldingReuse::slowing_downs`]）
    solves: &'a AtomicUsize,
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
    /// 等価なセルの減速材までのリングごとの物質（真空は `None`）。格子の領域は [`ring_targets`]。支持格子のリングの
    /// 物質は格子に使わないので含めない（枠は塊として遮蔽する）
    rings: Vec<Option<Material>>,
    /// 減速材のリングの始まり
    start: usize,
    /// 燃料の核種の実効断面積。UFG の範囲の群は燃料のリングの値を中性子束と体積で平均し、範囲の外は等価原理の値
    effective: Vec<Shielded>,
    /// [リング][核種] の実効断面積（`rings` と同じリング。真空は `None`）。中心のリングは UFG の範囲の群がリングの値、範囲の外が
    /// 等価原理の値で、ほかのリングは範囲の外が無限希釈
    ring_nuclides: Vec<Option<Arc<Vec<Shielded>>>>,
    /// [リング][核種] の UFG の範囲の群の弾性散乱の移行（`rings` と同じリング。[`UfgTransfer::RingFlux`] のときだけ）
    ring_transfers: Vec<Option<Transfers>>,
    /// [リング][UFG の範囲の群] の SPH 因子
    sph: Vec<Vec<f64>>,
    /// `sph` の列の群
    groups: std::ops::Range<usize>,
    /// 減速計算のエネルギーの格子点の数
    points: usize,
    /// 減速計算が上限なしで持つ量 (B。[`CellSpectrum::demand`])
    demand: u64,
    /// 減速計算の 1 次の変化（[`shield_fuel_tangent`] で頼んだときだけ）
    tangent: Option<Arc<PinTangent>>,
    /// [中心のリング][UFG の範囲の群] の、等価なセルの減速計算の中性子束と体積の積（比の残る差のログの重み。#0285）
    center_flux: Vec<Vec<f64>>,
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
    /// [リング][核種] の弾性散乱の移行（[`UfgTransfer::RingFlux`] のときだけ）
    transfers: Option<Vec<Vec<NuclideTransfer>>>,
    /// 減速計算の 1 次の変化（[`shield_fuel_tangent`] で頼んだときだけ）
    tangent: Option<PinTangent>,
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
    // 支持格子の物質が減速材の中性子束を変える効果を入れるため）
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
    let sum = tracing::debug_span!("equivalent-radius").in_scope(
        || -> Result<f64, SelfShieldingError> {
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
            Ok(sum)
        },
    )?;
    let radius = sum / (bottom - top) as f64;
    tracing::debug!(radius, "equivalent cell outer radius (cm)");
    let cylinder = Cylinder::with_points(
        radii[..start]
            .iter()
            .copied()
            .chain([radius])
            .chain(frame.map(|(_, area)| outer_radius(radius, area)))
            .collect(),
        pin.cp_points,
    )?;
    let solve = || -> Result<_, UfgError> {
        Ok(match pin.transfer {
            UfgTransfer::InfiniteDilution => (
                slowing_down_with(
                    library,
                    &cylinder,
                    &inputs,
                    range,
                    pin.collision,
                    pin.thinning,
                    pin.scattering,
                    pin.hold,
                )?,
                None,
            ),
            UfgTransfer::RingFlux => {
                let (spectrum, transfers) = slowing_down_transfers(
                    library,
                    &cylinder,
                    &inputs,
                    range,
                    pin.collision,
                    pin.thinning,
                    pin.scattering,
                    pin.shares,
                    pin.hold,
                )?;
                (spectrum, Some(transfers))
            }
        })
    };
    let (spectrum, transfers) = solve()?;
    pin.solves.fetch_add(1, Ordering::Relaxed);
    // 1 次の変化は、同じ設定の減速計算と一緒に解き直す（時間を向きのない減速計算と比べる測定。ADR-0247 の段 I2-2）。比べる
    // 向きのない減速計算も解き直す（最初の減速計算が自由ガスの核の D の値の表を埋めるので、同じ表を引く 2 度目と比べる）
    let tangent = match pin.tangent {
        None => None,
        Some(names) => {
            let clock = Instant::now();
            solve()?;
            let plain = clock.elapsed().as_secs_f64();
            let directions: Vec<Vec<Vec<f64>>> = names
                .iter()
                .map(|name| {
                    inputs
                        .iter()
                        .map(|m| {
                            m.iter()
                                .flat_map(|m| &m.composition)
                                .map(|(n, d)| if n == name { *d } else { 0.0 })
                                .collect()
                        })
                        .collect()
                })
                .collect();
            let clock = Instant::now();
            let shares = (pin.transfer == UfgTransfer::RingFlux).then_some(pin.shares);
            let (again, again_transfers, tangent) = slowing_down_tangent(
                library,
                &cylinder,
                &inputs,
                range,
                pin.collision,
                pin.thinning,
                pin.scattering,
                shares,
                pin.hold,
                &directions,
            )?;
            let seconds = clock.elapsed().as_secs_f64();
            // 1 次の変化を求めても減速計算は同じ（違えば 1 次の変化は別の問題の微分）
            if again != spectrum || again_transfers != transfers {
                return Err(SelfShieldingError::TangentSpectrum(cell));
            }
            Some(PinTangent::new(
                names, &inputs, pin.rings, &spectrum, tangent, plain, seconds,
            ))
        }
    };
    Ok(Some(EquivalentCell {
        radius,
        cylinder,
        inputs,
        start,
        frame: frame.is_some(),
        groups: top..bottom,
        spectrum,
        transfers,
        tangent,
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
        transfers,
        tangent,
    } = &equivalent;
    let (top, bottom) = (groups.start, groups.end);
    let bounds = library.groups.bounds();

    let ng = library.groups.len();
    let infinite = *library.sigma0.first().ok_or(LookupError::NoSigma0)?;
    let mut rings = Vec::with_capacity(inputs.len());
    let mut ring_nuclides = Vec::with_capacity(inputs.len());
    let mut ring_transfers = Vec::with_capacity(inputs.len());
    tracing::debug_span!("ring-materials").in_scope(|| -> Result<(), SelfShieldingError> {
        for (j, input) in inputs.iter().enumerate() {
            let Some(m) = input else {
                rings.push(None);
                ring_nuclides.push(None);
                ring_transfers.push(None);
                continue;
            };
            let ring_transfer = transfers.as_ref().map(|t| &t[j]);
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
                    let transfer = ring_transfer.map(|t| &t[k]);
                    let nuclide =
                        shielded_nuclide_with(library, name, m.temperature, &full, transfer)?;
                    fulls.push(full);
                    Ok((nuclide, *density))
                })
                .collect::<Result<Vec<_>, SelfShieldingError>>()?;
            // 燃料のリングの [`RING_CAPTURES`] の核種の、UFG の範囲の群をリングの中性子束で縮約した 1 群の捕獲断面積（自己遮蔽の
            // 条件を比べるログ。#0286）
            if j < pin.rings && tracing::enabled!(tracing::Level::DEBUG) {
                let flux = &spectrum.flux[j];
                for ((name, _), s) in m.composition.iter().zip(&spectrum.nuclides[j]) {
                    if RING_CAPTURES.contains(&name.as_str()) {
                        let rate: f64 = s.capture.iter().zip(flux).map(|(c, f)| c * f).sum();
                        let capture = rate / flux.iter().sum::<f64>();
                        tracing::debug!(
                            cell,
                            ring = j,
                            nuclide = name.as_str(),
                            capture,
                            "UFG ring capture"
                        );
                    }
                }
            }
            ring_nuclides.push(Some(Arc::new(fulls)));
            ring_transfers.push(ring_transfer.map(|t| Arc::new(t.clone())));
            let material = Material::from_nuclides(&nuclides, ng).map_err(|source| {
                SelfShieldingError::Material {
                    name: m.name.clone(),
                    source,
                }
            })?;
            rings.push(Some(material));
        }
        Ok(())
    })?;

    // SPH 因子: 同じセルを多群で解いた中性子束のリング間の形が UFG と一致するよう、リングの物質に掛ける。
    // 範囲の上の群の中性子束は UFG と同じく F = 1（群のレサジー幅）
    let multigroup: Vec<Material> = rings
        .iter()
        .map(|m| m.clone().unwrap_or_else(|| Material::void(ng)))
        .collect();
    let above: Vec<f64> = (0..top).map(|g| (bounds[g] / bounds[g + 1]).ln()).collect();
    // 支持格子のリングの断面積は格子に使わないので、因子を掛けない
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

    // 支持格子のリングの物質は格子に使わない（枠は塊として遮蔽する）
    if *frame {
        rings.pop();
        ring_nuclides.pop();
        ring_transfers.pop();
    }

    // 中心の物質（燃料か棒）の平均: そのリングの中性子束と体積の重み
    let effective = center_average(carlvik, &equivalent, pin.rings);
    let volumes = cylinder.volumes();
    let center_flux = (0..pin.rings)
        .map(|j| spectrum.flux[j].iter().map(|f| f * volumes[j]).collect())
        .collect();
    Ok(Some(CellShielding {
        radius: *radius,
        rings,
        start: *start,
        effective,
        ring_nuclides,
        ring_transfers,
        sph,
        groups: top..bottom,
        points: spectrum.points,
        demand: spectrum.demand,
        tangent: tangent.clone().map(Arc::new),
        center_flux,
    }))
}

/// 状態と基準の値の比（基準が 0 なら 1。比で変えない）
fn ratio(state: f64, anchor: f64) -> f64 {
    if anchor != 0.0 { state / anchor } else { 1.0 }
}

/// `target` の UFG の範囲の群 `groups` の反応ごとの値を、`base` の値に `state` / `anchor` の比を掛けたものにする
fn scale_by_ratio(
    target: &mut Shielded,
    base: &Shielded,
    state: &Shielded,
    anchor: &Shielded,
    groups: std::ops::Range<usize>,
) {
    fn reactions(s: &Shielded) -> [&Vec<f64>; 5] {
        [&s.total, &s.elastic, &s.fission, &s.nu_fission, &s.capture]
    }
    let (b, s, a) = (reactions(base), reactions(state), reactions(anchor));
    let out = [
        &mut target.total,
        &mut target.elastic,
        &mut target.fission,
        &mut target.nu_fission,
        &mut target.capture,
    ];
    for (x, values) in out.into_iter().enumerate() {
        for g in groups.clone() {
            values[g] = b[x][g] * ratio(s[x][g], a[x][g]);
        }
    }
}

/// ピンのセル `pin.cell` の等価なセルのリングの物質を、基準 `anchor` の同じセルの UFG の結果から比で作る（[`shield_fuel_ratio`]）。
/// 燃料のリングは、基準のリングの実効断面積に、状態の等価原理の値 `carlvik`（状態の Dancoff 係数 `dancoff` で解いたもの）と基準の
/// 等価原理の値（基準の Dancoff 係数と組成で解き直す）の比を掛ける（診断の [`RatioFuel`] では状態の UFG から作る）。target
/// `ratio_check` の DEBUG のときは、燃料のリングの [`RING_CAPTURES`] の核種の捕獲の比を群ごとにログに出す（[`log_ratio_check`]）。
/// 燃料でないリング、移行、SPH 因子は `anchor.options` で作り、
/// リングの物質は [`shield_cell`] と同じ作り方（核種の群定数を数密度で混ぜて SPH 因子を掛ける）。基準のピンに UFG の範囲が無ければ
/// `None`（範囲の外の等価原理だけ）
#[allow(clippy::too_many_arguments)]
fn ratio_cell(
    library: &Library,
    pin: &PinRings,
    materials: &[MaterialInput],
    dilute: &[Material],
    dancoff: &[f64],
    carlvik: &[Shielded],
    ring_materials: &[usize],
    options: &ShieldingOptions,
    anchor: &RatioAnchor,
) -> Result<Option<CellShielding>, SelfShieldingError> {
    let cell = pin.cell;
    let mismatch = || SelfShieldingError::AnchorPin(cell);
    let anchor_pin = anchor
        .shielding
        .pins
        .iter()
        .find(|p| p.cell == cell)
        .ok_or_else(mismatch)?;
    let (Some(anchor_sph), Some(radius)) = (&anchor_pin.sph, anchor_pin.equivalent) else {
        return Ok(None);
    };
    let groups = anchor_sph.groups.clone();
    let start = anchor_sph.moderator_ring;
    // 基準の等価原理の値（範囲の外の群と同じ式。基準の Dancoff 係数と、基準の物質の組成と温度）
    let radius_fuel = pin
        .lattice
        .cell_type(cell)
        .pin()
        .ok_or(SelfShieldingError::FuelShape(cell))?
        .radii()[pin.rings - 1];
    let anchor_equivalence = Equivalence {
        chord: 2.0 * radius_fuel,
        terms: anchor_pin
            .dancoff
            .iter()
            .map(|&c| self::carlvik(c))
            .collect::<Result<_, _>>()?,
    };
    let anchor_fuel = pin_fuel(pin.lattice, cell, anchor.materials, ring_materials);
    let anchor_carlvik =
        shield_material(library, &anchor_fuel, &anchor_equivalence, options, cell)?;
    let geometry = pin
        .lattice
        .cell_type(cell)
        .pin()
        .ok_or(SelfShieldingError::FuelShape(cell))?;
    let targets = ring_targets(geometry, pin.lattice.cell_regions(cell).start, start);
    // 状態の UFG の減速計算（診断の燃料のリングと状態の SPH 因子とログのときだけ。時間は状態ごとに解くのと同じ）
    let check = tracing::enabled!(target: "ratio_check", tracing::Level::DEBUG);
    let state_ufg = if matches!(
        anchor.options.fuel,
        RatioFuel::StateAverage | RatioFuel::State
    ) || anchor.options.sph == RatioSph::State
        || matches!(
            anchor.options.rings,
            RatioRings::State | RatioRings::StateSections
        )
        || check
    {
        Some(
            shield_cell(
                library,
                pin,
                materials,
                dilute,
                dancoff,
                carlvik,
                |inputs| cell_range(library, inputs),
            )?
            .ok_or_else(mismatch)?,
        )
    } else {
        None
    };
    let state_cell = || state_ufg.as_ref().ok_or_else(mismatch);
    // 静止した標的の減速計算の、状態と基準の結果（比のどれかが使うときだけ）。基準の結果は同じ基準の状態の間で持ち回る
    let light = if anchor.options.fuel == RatioFuel::StaticRatio
        || anchor.options.rings == RatioRings::StaticRatio
        || anchor.options.sph == RatioSph::StaticRatio
    {
        let static_pin = PinRings {
            scattering: Scattering::Static,
            transfer: UfgTransfer::InfiniteDilution,
            tangent: None,
            ..*pin
        };
        let range = |inputs: &[Option<&MaterialInput>]| cell_range(library, inputs);
        let at_state = tracing::debug_span!("static-state").in_scope(|| {
            shield_cell(
                library,
                &static_pin,
                materials,
                dilute,
                dancoff,
                carlvik,
                range,
            )
        })?;
        let at_anchor = anchor.light.cell(cell, || {
            let anchor_dilute = anchor.light.dilute(library, anchor.materials)?;
            tracing::debug_span!("static-anchor")
                .in_scope(|| {
                    shield_cell(
                        library,
                        &static_pin,
                        anchor.materials,
                        anchor_dilute,
                        &anchor_pin.dancoff,
                        &anchor_carlvik,
                        range,
                    )
                })?
                .ok_or_else(mismatch)
        })?;
        Some((at_state.ok_or_else(mismatch)?, at_anchor))
    } else {
        None
    };
    let light_rings = |j: usize| -> Result<(&[Shielded], &[Shielded]), SelfShieldingError> {
        let (s, a) = light.as_ref().ok_or_else(mismatch)?;
        let (Some(s), Some(a)) = (&s.ring_nuclides[j], &a.ring_nuclides[j]) else {
            return Err(mismatch());
        };
        Ok((s.as_slice(), a.as_slice()))
    };
    let ng = library.groups.len();
    let infinite = *library.sigma0.first().ok_or(LookupError::NoSigma0)?;
    let rows = |transfer: &NuclideTransfer, full: &Shielded, base: &Shielded| -> NuclideTransfer {
        match anchor.options.transfer {
            RatioTransfer::Reference => transfer.clone(),
            RatioTransfer::ElasticRatio => NuclideTransfer {
                groups: transfer.groups.clone(),
                rows: transfer
                    .groups
                    .clone()
                    .zip(&transfer.rows)
                    .map(|(g, row)| {
                        let r = ratio(full.elastic[g], base.elastic[g]);
                        row.iter().map(|v| v * r).collect()
                    })
                    .collect(),
            },
        }
    };
    let mut rings = Vec::with_capacity(targets.len());
    let mut ring_nuclides = Vec::with_capacity(targets.len());
    let mut ring_transfers = Vec::with_capacity(targets.len());
    for (j, regions) in targets.iter().enumerate() {
        let region = *regions.first().ok_or_else(mismatch)?;
        let index = pin.region_index[region];
        if index == pin.void {
            rings.push(None);
            ring_nuclides.push(None);
            ring_transfers.push(None);
            continue;
        }
        let (m, a) = (&materials[index], &anchor.materials[index]);
        let base = anchor
            .shielding
            .region_cross_sections(region)
            .ok_or_else(mismatch)?;
        let base_transfers = anchor.shielding.region_transfers(region);
        if base.len() != a.composition.len()
            || base_transfers.is_some_and(|t| t.len() != base.len())
        {
            return Err(mismatch());
        }
        // 基準と状態は同じ燃焼点の、ホウ素が公称の状態なので、リングの核種の並びは同じ（違えば別の組成の基準）
        if !m
            .composition
            .iter()
            .map(|(n, _)| n)
            .eq(a.composition.iter().map(|(n, _)| n))
        {
            return Err(mismatch());
        }
        let fuel = j < pin.rings;
        let mut fulls = Vec::with_capacity(m.composition.len());
        let mut transfers = Vec::with_capacity(m.composition.len());
        let mut nuclides = Vec::with_capacity(m.composition.len());
        for (k, (name, density)) in m.composition.iter().enumerate() {
            let (full, transfer) = if fuel {
                let by_ratio = |state: &Shielded, at_anchor: &Shielded| {
                    let mut full = carlvik[k].clone();
                    scale_by_ratio(&mut full, &base[k], state, at_anchor, groups.clone());
                    let transfer = base_transfers.map(|t| rows(&t[k], &full, &base[k]));
                    (full, transfer)
                };
                match anchor.options.fuel {
                    RatioFuel::Equivalence => by_ratio(&carlvik[k], &anchor_carlvik[k]),
                    RatioFuel::StateAverage => {
                        by_ratio(&state_cell()?.effective[k], &anchor_pin.effective[k])
                    }
                    RatioFuel::StaticRatio => {
                        let (s, a) = light_rings(j)?;
                        by_ratio(&s[k], &a[k])
                    }
                    RatioFuel::State => {
                        let s = state_cell()?;
                        let rings = s.ring_nuclides[j].as_ref().ok_or_else(mismatch)?;
                        (
                            rings[k].clone(),
                            s.ring_transfers[j].as_ref().map(|t| t[k].clone()),
                        )
                    }
                }
            } else {
                let mut full = library.shielded(name, m.temperature, &vec![infinite; ng])?;
                match anchor.options.rings {
                    RatioRings::InfiniteDilution => (full, None),
                    RatioRings::Reference => {
                        for g in groups.clone() {
                            full.total[g] = base[k].total[g];
                            full.elastic[g] = base[k].elastic[g];
                            full.fission[g] = base[k].fission[g];
                            full.nu_fission[g] = base[k].nu_fission[g];
                            full.capture[g] = base[k].capture[g];
                        }
                        (
                            full,
                            base_transfers.map(|t| rows(&t[k], &base[k], &base[k])),
                        )
                    }
                    RatioRings::ShieldingFactor => {
                        let state = full.clone();
                        let at_anchor =
                            library.shielded(name, a.temperature, &vec![infinite; ng])?;
                        scale_by_ratio(&mut full, &base[k], &state, &at_anchor, groups.clone());
                        let transfer = base_transfers.map(|t| rows(&t[k], &full, &base[k]));
                        (full, transfer)
                    }
                    RatioRings::StaticRatio => {
                        let (s, a) = light_rings(j)?;
                        scale_by_ratio(&mut full, &base[k], &s[k], &a[k], groups.clone());
                        let transfer = base_transfers.map(|t| rows(&t[k], &full, &base[k]));
                        (full, transfer)
                    }
                    RatioRings::State | RatioRings::StateSections => {
                        let s = state_cell()?;
                        let ring = s.ring_nuclides[j].as_ref().ok_or_else(mismatch)?;
                        let full = ring[k].clone();
                        let transfer = match anchor.options.rings {
                            RatioRings::State => s.ring_transfers[j].as_ref().map(|t| t[k].clone()),
                            _ => base_transfers.map(|t| rows(&t[k], &full, &base[k])),
                        };
                        (full, transfer)
                    }
                }
            };
            let nuclide =
                shielded_nuclide_with(library, name, m.temperature, &full, transfer.as_ref())?;
            nuclides.push((nuclide, *density));
            fulls.push(full);
            transfers.push(transfer);
        }
        let material = Material::from_nuclides(&nuclides, ng).map_err(|source| {
            SelfShieldingError::Material {
                name: m.name.clone(),
                source,
            }
        })?;
        rings.push(Some(material));
        ring_nuclides.push(Some(Arc::new(fulls)));
        // 移行の行は、基準が持てばどの核種にもあり、無限希釈の移行で解いた基準と、燃料でないリングを無限希釈にしたときは無い
        let transfers: Option<Vec<NuclideTransfer>> = transfers.into_iter().collect();
        ring_transfers.push(base_transfers.and(transfers).map(Arc::new));
    }
    let sph = match anchor.options.sph {
        RatioSph::Reference => anchor_sph.factors.clone(),
        RatioSph::State => state_cell()?.sph.clone(),
        RatioSph::StaticRatio => {
            let (s, a) = light.as_ref().ok_or_else(mismatch)?;
            if s.sph.len() != anchor_sph.factors.len() || a.sph.len() != s.sph.len() {
                return Err(mismatch());
            }
            anchor_sph
                .factors
                .iter()
                .zip(&s.sph)
                .zip(&a.sph)
                .map(|((base, state), at_anchor)| {
                    base.iter()
                        .zip(state)
                        .zip(at_anchor)
                        .map(|((b, s), a)| b * ratio(*s, *a))
                        .collect()
                })
                .collect()
        }
    };
    if sph.len() < rings.len() || sph.iter().any(|row| row.len() != groups.len()) {
        return Err(mismatch());
    }
    for (ring, (material, mu)) in rings.iter_mut().zip(&sph).enumerate() {
        if let Some(m) = material {
            let mut all = vec![1.0; ng];
            all[groups.clone()].copy_from_slice(mu);
            *m = m
                .scaled(&all)
                .map_err(|source| SelfShieldingError::Material {
                    name: format!("ring {ring} of cell {cell}"),
                    source,
                })?;
        }
    }
    // ピンの平均（範囲の外は状態の等価原理の値）も、基準の平均に同じ比を掛ける（状態の UFG の平均の比なら、状態の UFG の平均）
    let effective = match anchor.options.fuel {
        RatioFuel::Equivalence => carlvik
            .iter()
            .zip(anchor_pin.effective.iter())
            .zip(&anchor_carlvik)
            .map(|((state, base), at_anchor)| {
                let mut s = state.clone();
                scale_by_ratio(&mut s, base, state, at_anchor, groups.clone());
                s
            })
            .collect(),
        RatioFuel::StateAverage | RatioFuel::State => state_cell()?.effective.clone(),
        RatioFuel::StaticRatio => {
            let (s, a) = light.as_ref().ok_or_else(mismatch)?;
            carlvik
                .iter()
                .zip(anchor_pin.effective.iter())
                .zip(s.effective.iter().zip(&a.effective))
                .map(|((state, base), (light_state, light_anchor))| {
                    let mut out = state.clone();
                    scale_by_ratio(&mut out, base, light_state, light_anchor, groups.clone());
                    out
                })
                .collect()
        }
    };
    if check {
        log_ratio_check(
            library,
            pin,
            dancoff,
            &RatioCheck {
                materials,
                carlvik,
                anchor_carlvik: &anchor_carlvik,
                anchor_pin,
                anchor,
                state: state_cell()?,
                targets: &targets,
            },
        )?;
    }
    // この状態で解いた減速計算の格子点の数と量（持ち回った基準の静止した標的の減速計算は数えない）
    let solved = || state_ufg.iter().chain(light.as_ref().map(|(s, _)| s));
    let points = solved().map(|c| c.points).max().unwrap_or(0);
    let demand = solved().map(|c| c.demand).max().unwrap_or(0);
    let radius = match anchor.options.fuel {
        RatioFuel::State => state_cell()?.radius,
        RatioFuel::Equivalence | RatioFuel::StateAverage | RatioFuel::StaticRatio => radius,
    };
    Ok(Some(CellShielding {
        radius,
        rings,
        start,
        effective,
        ring_nuclides,
        ring_transfers,
        sph,
        groups,
        points,
        demand,
        tangent: None,
        center_flux: Vec::new(),
    }))
}

/// 物質 `materials` の無限希釈の群定数に、最後に真空の物質を足したもの（[`shield_fuel_by`] の物質の番号と同じ並び）
fn mixed(
    library: &Library,
    materials: &[MaterialInput],
) -> Result<Vec<Material>, SelfShieldingError> {
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
    dilute.push(Material::void(library.groups.len()));
    Ok(dilute)
}

/// [`log_ratio_check`] が比べる値
struct RatioCheck<'a> {
    /// 状態の物質
    materials: &'a [MaterialInput],
    /// 状態と基準の等価原理の値
    carlvik: &'a [Shielded],
    anchor_carlvik: &'a [Shielded],
    /// 基準の同じセルのピン
    anchor_pin: &'a PinShielding,
    anchor: &'a RatioAnchor<'a>,
    /// 状態の UFG の減速計算
    state: &'a CellShielding,
    /// 等価なセルのリングごとの格子の領域（[`ring_targets`]）
    targets: &'a [Vec<usize>],
}

/// 比で求める状態のピン `pin.cell` の燃料のリングの [`RING_CAPTURES`] の核種の、UFG の範囲の群ごとの捕獲の、状態と基準の比を target
/// `ratio_check` の DEBUG のログに出す（比の残る差の原因を分ける。#0285）: `equivalence` は等価原理の比（[`RatioFuel::Equivalence`]
/// が掛ける比）、`state` はリングの UFG の比、`mean` は UFG のピンの平均の比、`dancoff` は基準の物質のまま状態の Dancoff 係数
/// `dancoff` で解き直した UFG のリングの比（減速材の物質の変化を除き、等価なセルの外径だけを変えた分）。`rate` は状態の UFG の
/// そのリングと群の捕獲の反応率（数密度、実効断面積、等価なセルの中性子束と体積の積。群とリングの重み）
fn log_ratio_check(
    library: &Library,
    pin: &PinRings,
    dancoff: &[f64],
    check: &RatioCheck,
) -> Result<(), SelfShieldingError> {
    let cell = pin.cell;
    let mismatch = || SelfShieldingError::AnchorPin(cell);
    let anchor = check.anchor;
    // 基準の物質（真空は最後）で、状態の Dancoff 係数の等価なセルを解く
    let anchor_dilute = mixed(library, anchor.materials)?;
    let moved = shield_cell(
        library,
        pin,
        anchor.materials,
        &anchor_dilute,
        dancoff,
        check.anchor_carlvik,
        |inputs| cell_range(library, inputs),
    )?
    .ok_or_else(mismatch)?;
    let state = check.state;
    let groups = state.groups.clone();
    for (j, regions) in check.targets.iter().enumerate().take(pin.rings) {
        let region = *regions.first().ok_or_else(mismatch)?;
        let base = anchor
            .shielding
            .region_cross_sections(region)
            .ok_or_else(mismatch)?;
        let rings = state.ring_nuclides[j].as_ref().ok_or_else(mismatch)?;
        let moved_rings = moved.ring_nuclides[j].as_ref().ok_or_else(mismatch)?;
        let composition = &check.materials[pin.region_index[region]].composition;
        for (k, (name, density)) in composition.iter().enumerate() {
            if !RING_CAPTURES.contains(&name.as_str()) {
                continue;
            }
            for g in groups.clone() {
                let capture = |s: &Shielded| s.capture[g];
                tracing::debug!(
                    target: "ratio_check",
                    cell,
                    ring = j,
                    nuclide = name.as_str(),
                    group = g,
                    equivalence = ratio(capture(&check.carlvik[k]), capture(&check.anchor_carlvik[k])),
                    state = ratio(capture(&rings[k]), capture(&base[k])),
                    mean = ratio(capture(&state.effective[k]), capture(&check.anchor_pin.effective[k])),
                    dancoff = ratio(capture(&moved_rings[k]), capture(&base[k])),
                    rate = density * capture(&rings[k]) * state.center_flux[j][g - groups.start],
                    "ratio check of the ring capture"
                );
            }
        }
    }
    Ok(())
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
    for iteration in 1..=options.max_iterations {
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
            tracing::debug!(
                iterations = iteration,
                nuclides = tables.len(),
                dilute = dilute.iter().flatten().count(),
                "equivalence iterations"
            );
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
    shielded_nuclide_with(library, name, temperature, shielded, None)
}

/// [`shielded_nuclide`] の、UFG の範囲の群の弾性散乱の移行 `transfer` のある群の散乱行列を、その移行で作るもの（ADR-0207）。
/// その群の P0 の行は、移行の行（弾性散乱。行の和は `shielded` の弾性散乱）に、ライブラリの行の和のうち無限希釈の弾性散乱を超える
/// 分（非弾性散乱など。UFG の減速計算は扱わない）をライブラリの行の形で配ったものを足す（ライブラリの行の和が正なら、行の和は
/// [`shielded_nuclide`] と同じ。行の和が 0 の群は、[`shielded_nuclide`] が 0 のまま残す弾性散乱を移行の行で入れる）。P1 以上の行は、
/// 要素ごとにライブラリの P_l/P0 を保つ（ライブラリの P0 が 0 の要素は 0）。ほかの群と、`transfer` がなければすべての群は
/// [`shielded_nuclide`] のまま
pub fn shielded_nuclide_with(
    library: &Library,
    name: &str,
    temperature: f64,
    shielded: &Shielded,
    transfer: Option<&NuclideTransfer>,
) -> Result<Nuclide, SelfShieldingError> {
    let mut nuclide = library.nuclide(name, temperature)?;
    let ng = library.groups.len();
    let infinite = library.sigma0.first().ok_or(LookupError::NoSigma0)?;
    let dilute = library.shielded(name, temperature, &vec![*infinite; ng])?;
    // 移行の行で作る群の、ライブラリの P0〜P_L の行（下で ADR-0076 の比を掛ける前のもの）
    let library_rows: Vec<(usize, Vec<Vec<f64>>)> = match transfer {
        None => Vec::new(),
        Some(t) => {
            if t.groups.end > ng {
                return Err(SelfShieldingError::TransferGroups {
                    name: name.into(),
                    given: t.groups.end,
                    groups: ng,
                });
            }
            if t.rows.len() != t.groups.len() {
                return Err(SelfShieldingError::TransferGroups {
                    name: name.into(),
                    given: t.rows.len(),
                    groups: t.groups.len(),
                });
            }
            t.groups
                .clone()
                .map(|g| {
                    let rows = std::iter::once(&nuclide.scatter)
                        .chain(&nuclide.anisotropic)
                        .map(|m| m[g].clone())
                        .collect();
                    (g, rows)
                })
                .collect()
        }
    };
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
    if let Some(t) = transfer {
        for ((g, rows), row) in library_rows.into_iter().zip(&t.rows) {
            if row.len() != ng {
                return Err(SelfShieldingError::TransferGroups {
                    name: name.into(),
                    given: row.len(),
                    groups: ng,
                });
            }
            let library_row = &rows[0];
            let sum: f64 = library_row.iter().sum();
            // 行の和が無限希釈の弾性散乱より小さい行（丸めを超えて。和が 0 の行も）は、弾性散乱でない分を配れない（ライブラリの
            // 不整合）
            let other = if sum > 0.0 {
                (sum - dilute.elastic[g]) / sum
            } else if dilute.elastic[g] > 0.0 {
                f64::NEG_INFINITY
            } else {
                0.0
            };
            if other < -NON_ELASTIC_ROUNDING {
                return Err(SelfShieldingError::ElasticAboveRow {
                    name: name.into(),
                    group: g,
                    row: sum,
                    elastic: dilute.elastic[g],
                });
            }
            let p0: Vec<f64> = row
                .iter()
                .zip(library_row)
                .map(|(t, r)| t + other * r)
                .collect();
            // 分布の Legendre のモーメントは |P_l| ≤ P0 なので、要素ごとの P_l/P0 は 1 以下（超えればライブラリの不整合）
            for library_moment in &rows[1..] {
                if let Some(to) = (0..ng).find(|&h| {
                    library_moment[h].abs() > library_row[h] * (1.0 + NON_ELASTIC_ROUNDING)
                }) {
                    return Err(SelfShieldingError::MomentAboveP0 {
                        name: name.into(),
                        group: g,
                        to,
                    });
                }
            }
            for (moment, library_moment) in nuclide.anisotropic.iter_mut().zip(&rows[1..]) {
                moment[g] = library_moment
                    .iter()
                    .zip(library_row)
                    .zip(&p0)
                    .map(|((a, r), p)| if *r != 0.0 { a * p / r } else { 0.0 })
                    .collect();
            }
            nuclide.scatter[g] = p0;
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
pub(crate) mod tests;
