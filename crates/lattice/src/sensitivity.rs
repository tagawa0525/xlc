//! k の核種の数密度の感度（ADR-0247 の段 2。#0279）。
//!
//! 核種 i の数密度を、それを持つ領域ごとに N_i,r に比例して変える向き（ADR-0221 の摂動と同じ向き）の k の微分
//! S_i = Σ_r N_i,r ∂k/∂N_i,r を、前向きと随伴の解（`moc::solve_importance`）から 1 次の摂動で求める。自己遮蔽は状態の解のまま
//! 固定し（実効断面積と SPH 因子を変えない）、領域の物質の微分は、その核種の自己遮蔽した群定数（数密度 1）に SPH 因子を掛けたもの。
//! 均質化した数密度 N̄_i で割った dk/dN̄_i = S_i / N̄_i も書く。
//!
//! `direct` を与えれば、同じ向きに数密度を (1 ± ε) 倍と (1 ± ε/2) 倍した固有値計算（自己遮蔽は固定）の中心差分も書く（検証の
//! 参照。2 つの刻みから Richardson の外挿で ε² の誤差を除く）。

use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, anyhow, bail};
use serde::Serialize;
use xlc_format::library::Library;

use crate::cmfd::CmfdOptions;
use crate::fewgroup::{coarse_groups, homogenize};
use crate::geometry::Lattice;
use crate::input::{MaterialInput, VOID};
use crate::material::Material;
use crate::moc::{
    Acceleration, Generalized, GeneralizedMethod, GeneralizedOptions, MocError, Options,
    ResponseImportance, Solution, Start, response_importance, solve_generalized, solve_importance,
    solve_symmetric_from,
};
use crate::shielding::{FuelShielding, shielded_nuclide_with};
use crate::symmetry::Detected;
use crate::tracking::Tracking;

/// 物質の和を物質の微分の和と比べる相対の許容値（同じ群定数を足す順の丸めの違いだけを許す）
const RECONSTRUCTION_TOLERANCE: f64 = 1.0e-10;

/// 応答の値をテーブルの少数群定数（`fewgroup::homogenize`）と比べる相対の許容値。同じ正の項の和を、掛ける順と足す順を変えて
/// 割った丸めの違いだけを許す（n 項の和の丸めは n ε まで。2A の本番の形の格子は領域と細かい群の項が約 10⁵ で 1e-11 ほど。
/// 3×3 の格子で 4e-13 だった）
const TABLE_TOLERANCE: f64 = 1.0e-10;

/// `[sensitivity]` の入力（ADR-0247）
#[derive(Debug, Clone, PartialEq)]
pub struct SensitivityInput {
    /// 感度を求める核種
    pub nuclides: Vec<String>,
    /// 結果の TOML（入力ファイルのディレクトリからの相対パス）
    pub output: PathBuf,
    /// 直接の摂動の割合 ε（`None` なら直接の摂動を解かない）
    pub direct: Option<f64>,
    /// この入力のすべての MOC の判定（`None` なら `Options` の既定）
    pub k_tolerance: Option<f64>,
    pub flux_tolerance: Option<f64>,
    /// 一般化随伴の固定源の問題の比較（ADR-0247 の段 3。`None` なら解かない）
    pub generalized: Option<GeneralizedInput>,
    /// 少数群定数の数密度の微分（ADR-0247 の段 4。`None` なら求めない）
    pub few_group: Option<FewGroupInput>,
}

/// `[sensitivity.few-group]` の入力: 格子全体の少数群定数（テーブルの全断面積、νΣf、Σf、P0 の散乱行列）の核種の数密度の微分を、
/// 応答ごとの一般化随伴から求める（ADR-0247 の段 4）
#[derive(Debug, Clone, PartialEq)]
pub struct FewGroupInput {
    /// 少数群の境界 (eV、降順。ライブラリの群の境界から選ぶ)
    pub bounds: Vec<f64>,
    /// 一般化随伴の GMRES の判定と再始動の回数（`GeneralizedOptions`）
    pub tolerance: f64,
    pub restart: usize,
    /// 自己遮蔽を摂動の組成で解き直す核種（`nuclides` から選ぶ）。暗黙の効きを、自己遮蔽した物質の差を一般化随伴に掛ける 1 次の
    /// 変化（ADR-0247 の案 I1）と、自己遮蔽を解き直した直接の摂動（案 I3 の参照）の両方で、直接の摂動と同じ刻みの中心差分で求める
    pub shielded: Vec<String>,
}

/// `[sensitivity.generalized]` の入力: 少数群定数の比の応答の一般化随伴を、収束の方法ごとに解いて比べる（ADR-0247 の段 3）
#[derive(Debug, Clone, PartialEq)]
pub struct GeneralizedInput {
    /// 少数群の境界 (eV、降順。ライブラリの群の境界から選ぶ)
    pub bounds: Vec<f64>,
    pub responses: Vec<Response>,
    pub methods: Vec<Method>,
    /// 一般化随伴の反復の判定（`GeneralizedOptions::tolerance`）。方法ごとに、それぞれの判定で解く
    pub tolerances: Vec<f64>,
    /// `methods` のうち、すべての源を一緒にスイープしても解く方法（`GeneralizedOptions::batch`）
    pub batched: Vec<Method>,
    /// 差を比べる参照の方法と判定（`None` なら比べない）
    pub reference: Option<(Method, f64)>,
}

/// 格子全体で均質化した少数群定数の比の応答 R = ⟨h₁, φ⟩ / ⟨h₂, φ⟩（h₂ は少数群 G の細かい群で 1。番号は 0 始まり）
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Response {
    /// Σt,G（h₁ = Σt）
    Total(usize),
    /// νΣf,G
    NuFission(usize),
    /// Σf,G
    Fission(usize),
    /// P0 の散乱 Σs,G→H（h₁ は細かい群 g ∈ G の、H の細かい群への散乱の和）
    Scatter(usize, usize),
}

impl std::fmt::Display for Response {
    /// 入力と結果の名前（少数群は 1 始まり）
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Response::Total(g) => write!(f, "total-{}", g + 1),
            Response::NuFission(g) => write!(f, "nu-fission-{}", g + 1),
            Response::Fission(g) => write!(f, "fission-{}", g + 1),
            Response::Scatter(g, h) => write!(f, "scatter-{}-{}", g + 1, h + 1),
        }
    }
}

/// 一般化随伴の収束の方法（`GeneralizedMethod`。CMFD の判定は状態を解いた MOC のもの）
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Method {
    SourceIteration,
    Cmfd,
    Gmres { restart: usize },
}

impl std::fmt::Display for Method {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Method::SourceIteration => write!(f, "source-iteration"),
            Method::Cmfd => write!(f, "cmfd"),
            Method::Gmres { .. } => write!(f, "gmres"),
        }
    }
}

impl SensitivityInput {
    /// `base` の判定を、入力の判定に替えたもの。CMFD の判定は既定と入力の細かい方（ADR-0110 の、CMFD を MOC より細かく解く前提）
    pub fn options(&self, base: Options) -> Options {
        let k_tolerance = self.k_tolerance.unwrap_or(base.k_tolerance);
        let flux_tolerance = self.flux_tolerance.unwrap_or(base.flux_tolerance);
        let acceleration = match base.acceleration {
            Acceleration::Cmfd(c) => Acceleration::Cmfd(CmfdOptions {
                k_tolerance: c.k_tolerance.min(k_tolerance),
                flux_tolerance: c.flux_tolerance.min(flux_tolerance),
                ..c
            }),
            Acceleration::None => Acceleration::None,
        };
        Options {
            k_tolerance,
            flux_tolerance,
            acceleration,
            ..base
        }
    }
}

/// 感度を求める状態（解いた物質と、その作り方）
pub struct State<'a> {
    pub library: &'a Library,
    /// 状態の物質の入力（領域の物質の名前で引く）
    pub materials: &'a [MaterialInput],
    /// [領域] の物質の名前
    pub region_names: &'a [String],
    /// 解いた物質と、[領域] のその番号
    pub unique: &'a [Material],
    pub region_index: &'a [usize],
    /// 状態の燃料の自己遮蔽（自己遮蔽しなければ `None`）
    pub shielding: Option<&'a FuelShielding>,
    pub lattice: &'a Lattice,
    pub tracking: &'a Tracking,
    pub symmetry: &'a Detected,
}

/// 結果の TOML
#[derive(Debug, Serialize)]
pub struct Report {
    pub k: f64,
    pub adjoint_k: f64,
    pub forward_iterations: usize,
    pub adjoint_iterations: usize,
    pub k_tolerance: f64,
    pub flux_tolerance: f64,
    /// 直接の摂動の割合 ε（解かなければ書かない）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub epsilon: Option<f64>,
    pub seconds: Seconds,
    pub nuclides: Vec<NuclideSensitivity>,
    /// 一般化随伴の比較（`[sensitivity.generalized]` を与えたときだけ）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub generalized: Option<GeneralizedReport>,
    /// 少数群定数の微分（`[sensitivity.few-group]` を与えたときだけ）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub few_group: Option<FewGroupReport>,
}

/// 少数群定数の微分の結果（ADR-0247 の段 4）
#[derive(Debug, Serialize)]
pub struct FewGroupReport {
    pub bounds: Vec<f64>,
    /// 応答（テーブルの少数群定数）の名前（`values` と核種の値の [応答] の順）
    pub responses: Vec<String>,
    /// [応答] の格子全体の少数群定数
    pub values: Vec<f64>,
    /// [応答] の値と、テーブルに書く `fewgroup::homogenize` の値の相対の差の最大
    pub table_difference: f64,
    /// テーブルの少数群の核分裂のスペクトル（反応度に直す 2 群の k∞ に使う）
    pub chi: Vec<f64>,
    /// 一般化随伴の GMRES の判定と再始動の回数
    pub tolerance: f64,
    pub restart: usize,
    /// [応答] の一般化随伴のスイープの回数
    pub sweeps: Vec<usize>,
    pub seconds: FewGroupSeconds,
    pub nuclides: Vec<FewGroupNuclide>,
}

/// 少数群定数の微分の段ごとの時間 (s)
#[derive(Debug, Serialize)]
pub struct FewGroupSeconds {
    /// 一般化随伴と、その角度中性子束の積のスイープ
    pub generalized: f64,
    /// 核種ごとの 1 次の変化
    pub derivatives: f64,
    /// 自己遮蔽を解き直した摂動（自己遮蔽、1 次の変化、固有値計算）
    pub shielded: f64,
}

/// 核種の、少数群定数の数密度の微分
#[derive(Debug, Serialize)]
pub struct FewGroupNuclide {
    pub name: String,
    /// [応答] の dR/dδ（数密度を (1 + δ) 倍。自己遮蔽を固定した陽の分で、一般化随伴から。ADR-0247 の案 I0）。∂R/∂N̄ は
    /// `[[nuclides]]` の density で割る
    pub relative: Vec<f64>,
    /// 自己遮蔽を固定した直接の摂動（`direct` を与えたときだけ）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub direct: Option<Steps>,
    /// 自己遮蔽を解き直した摂動（`shielded` に挙げた核種だけ）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub shielded: Option<Shielded>,
}

/// [応答] の、刻み ε と ε/2 の中心差分と Richardson の外挿
#[derive(Debug, Serialize)]
pub struct Steps {
    /// [刻み][応答]
    pub central: [Vec<f64>; 2],
    pub richardson: Vec<f64>,
}

/// 刻み ε と ε/2 の中心差分と Richardson の外挿
#[derive(Debug, Serialize)]
pub struct ScalarSteps {
    pub central: [f64; 2],
    pub richardson: f64,
}

/// 自己遮蔽を摂動の組成で解き直した物質による変化
#[derive(Debug, Serialize)]
pub struct Shielded {
    /// [応答] の、状態の一般化随伴に掛けた 1 次の変化の中心差分（ADR-0247 の案 I1。陽の分を含む）
    pub first_order: Steps,
    /// [応答] の直接の摂動（案 I3 の参照）
    pub direct: Steps,
    /// k の、随伴に掛けた 1 次の変化と直接の摂動
    pub k_first_order: ScalarSteps,
    pub k_direct: ScalarSteps,
}

/// 一般化随伴の比較の結果（ADR-0247 の段 3）
#[derive(Debug, Serialize)]
pub struct GeneralizedReport {
    /// 応答の名前（`runs` の [応答] の順）
    pub responses: Vec<String>,
    pub runs: Vec<GeneralizedRun>,
    /// 差を比べた参照（`reference` を与えたときだけ）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reference: Option<GeneralizedRun>,
}

/// 1 つの収束の方法で、すべての応答の一般化随伴を解いた結果
#[derive(Debug, Serialize)]
pub struct GeneralizedRun {
    pub method: String,
    pub tolerance: f64,
    /// すべての源を一緒にスイープしたか
    pub batch: bool,
    /// すべての応答を解いた時間 (s)。粗メッシュの加速の問題を作る時間を含む
    pub seconds: f64,
    /// [応答] のスイープの回数
    pub sweeps: Vec<usize>,
    /// [応答] の ⟨Γ, Fφ⟩ / ⟨|Γ|, Fφ⟩
    pub orthogonality: Vec<f64>,
    /// [応答] の、参照との差の最大を参照の |Γ| の最大で割ったもの（参照がないか、解けなかったときは書かない）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub difference: Option<Vec<f64>>,
    /// 解けなかった理由（発散、打ち切り、粗メッシュの問題のエラー。そのときスイープの回数と ⟨Γ, Fφ⟩ は空）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub failure: Option<String>,
}

/// 段ごとの時間 (s)
#[derive(Debug, Serialize)]
pub struct Seconds {
    pub forward: f64,
    pub adjoint: f64,
    /// 前向きと随伴の角度中性子束の積のスイープと、内積の係数
    pub products: f64,
    /// 核種ごとの物質の微分と 1 次の摂動の内積
    pub derivatives: f64,
    /// 直接の摂動の固有値計算の和（解かなければ 0）
    pub direct: f64,
}

#[derive(Debug, Serialize)]
pub struct NuclideSensitivity {
    pub name: String,
    /// 格子全体で体積の重みを付けて均質化した数密度 N̄ (1/(barn cm))
    pub density: f64,
    /// S = Σ_r N_r ∂k/∂N_r（数密度を (1 + δ) 倍したときの dk/dδ）
    pub relative: f64,
    /// dk/dN̄ = S / N̄ (barn cm)。N̄ が 0 なら書かない（向きが定まらない。#0206）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub derivative: Option<f64>,
    /// 直接の摂動（`direct` を与えたときだけ）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub direct: Option<Direct>,
}

/// 数密度を (1 + δ) 倍した固有値計算の k と、その中心差分
#[derive(Debug, Serialize)]
pub struct Direct {
    /// δ = ε, −ε, ε/2, −ε/2 の k
    pub k: [f64; 4],
    /// 刻み ε と ε/2 の中心差分 (k(δ) − k(−δ)) / 2δ
    pub central: [f64; 2],
    /// Richardson の外挿 D(ε/2) + (D(ε/2) − D(ε)) / 3
    pub richardson: f64,
}

/// 自己遮蔽を解き直す口: 物質の入力から、自己遮蔽した物質と [領域] のその番号を作る（状態と同じ自己遮蔽の作り方）
pub type Reshield<'a> =
    dyn FnMut(&[MaterialInput]) -> anyhow::Result<(Vec<Material>, Vec<usize>)> + 'a;

/// 直接の摂動の刻み δ = ε, −ε, ε/2, −ε/2
fn steps_of(epsilon: f64) -> [f64; 4] {
    [epsilon, -epsilon, epsilon / 2.0, -epsilon / 2.0]
}

/// 刻み `steps_of(epsilon)` の値 `values` の中心差分と Richardson の外挿
fn scalar_steps(values: [f64; 4], epsilon: f64) -> ScalarSteps {
    let central = [
        (values[0] - values[1]) / (2.0 * epsilon),
        (values[2] - values[3]) / epsilon,
    ];
    ScalarSteps {
        central,
        richardson: central[1] + (central[1] - central[0]) / 3.0,
    }
}

/// [刻み][応答] の値 `values` の、応答ごとの中心差分と Richardson の外挿
fn response_steps(values: &[Vec<f64>; 4], epsilon: f64) -> Steps {
    let each: Vec<ScalarSteps> = (0..values[0].len())
        .map(|i| scalar_steps(values.each_ref().map(|v| v[i]), epsilon))
        .collect();
    Steps {
        central: [0, 1].map(|j| each.iter().map(|s| s.central[j]).collect()),
        richardson: each.iter().map(|s| s.richardson).collect(),
    }
}

/// [応答] の R（`response_values` の比）
fn ratios(values: &[(f64, f64)]) -> Vec<f64> {
    values.iter().map(|(a, b)| a / b).collect()
}

/// 少数群定数の微分の、状態の解から作る部分（`[sensitivity.few-group]`）
struct FewGroupState<'a> {
    input: &'a FewGroupInput,
    responses: Vec<Response>,
    group_of: Vec<usize>,
    solve: FewGroupSolve,
    /// 応答の値とテーブルの値の相対の差の最大
    table_difference: f64,
    /// テーブルの少数群の核分裂のスペクトル
    chi: Vec<f64>,
    /// 一般化随伴と係数の時間 (s)
    seconds: f64,
}

/// 状態 `state` を、`start` から（`None` なら平らな初期値から）`options` で解き、感度を求める。前向きの解と収束した状態も返す。
/// `[sensitivity.few-group]` の `shielded` は、自己遮蔽を `reshield` で解き直す（なければエラー）
pub fn evaluate(
    input: &SensitivityInput,
    state: &State,
    options: &Options,
    start: Option<&Start>,
    mut reshield: Option<&mut Reshield>,
) -> anyhow::Result<(Solution, Start, Report)> {
    // 解く前に分かる入力の誤りは、前向きと随伴の固有値計算の前に止める
    if let Some(f) = &input.few_group {
        coarse_groups(&state.library.groups, &f.bounds)
            .context("[sensitivity.few-group] bounds")?;
        if !f.shielded.is_empty() && reshield.is_none() {
            bail!(
                "[sensitivity.few-group] shielded solves the self-shielding of the state again, which needs a \
                 self-shielded state solved from its own composition (not from a reference state)"
            );
        }
    }
    let materials: Vec<&Material> = state
        .region_index
        .iter()
        .map(|&i| &state.unique[i])
        .collect();
    let clock = Instant::now();
    let (solution, converged, importance) = solve_importance(
        state.tracking,
        &materials,
        options,
        Some(state.symmetry),
        start,
    )?;
    let solve_seconds = clock.elapsed().as_secs_f64();
    let clock = Instant::now();
    let derivatives = Derivatives::new(state)?;
    let total_volume: f64 = solution.volumes.iter().sum();
    let mut nuclides = Vec::with_capacity(input.nuclides.len());
    // [核種] の（解いた物質の番号、物質の微分 N_r × 核種の群定数）
    let mut weighted_changes: Vec<Vec<(usize, Material)>> =
        Vec::with_capacity(input.nuclides.len());
    for name in &input.nuclides {
        // 組成にない核種は感度 0 と区別できないので止める（名前の誤りや、輸送しない核種の候補）
        if !derivatives.lists(name) {
            bail!("{name} is in no material of the solved state");
        }
        let changes = derivatives.changes(name);
        let mut density = 0.0;
        let mut weighted: Vec<(usize, Material)> = Vec::new();
        for (u, n, d) in &changes {
            for (region, &i) in state.region_index.iter().enumerate() {
                if i == *u {
                    density += solution.volumes[region] * n;
                }
            }
            weighted.push((*u, d.scaled(&vec![*n; d.groups()])?));
        }
        density /= total_volume;
        let relative = importance
            .k_derivative(&region_changes(&weighted, state.region_index))
            .with_context(|| format!("{name}: k derivative"))?;
        nuclides.push(NuclideSensitivity {
            name: name.clone(),
            density,
            relative,
            derivative: (density > 0.0).then(|| relative / density),
            direct: None,
        });
        // 少数群定数の微分だけが後で使う（核種と物質の数だけの群定数の行列を持たない）
        if input.few_group.is_some() {
            weighted_changes.push(weighted);
        }
    }
    let derivative_seconds = clock.elapsed().as_secs_f64();
    let few_group = match &input.few_group {
        Some(f) => Some(few_group_state(
            f,
            state,
            &materials,
            options,
            (&solution, &converged, importance.adjoint_start()),
        )?),
        None => None,
    };
    let clock = Instant::now();
    let mut few_group_nuclides = Vec::new();
    if let Some(fg) = &few_group {
        for (nuclide, weighted) in input.nuclides.iter().zip(&weighted_changes) {
            let relative = fg
                .solve
                .first_order(
                    &fg.responses,
                    &fg.group_of,
                    &solution,
                    &region_changes(weighted, state.region_index),
                )
                .with_context(|| format!("{nuclide}: few-group derivatives"))?;
            few_group_nuclides.push(FewGroupNuclide {
                name: nuclide.clone(),
                relative,
                direct: None,
                shielded: None,
            });
        }
    }
    let few_group_derivative_seconds = clock.elapsed().as_secs_f64();
    let clock = Instant::now();
    if let Some(epsilon) = input.direct {
        for (position, nuclide) in nuclides.iter_mut().enumerate() {
            let changes = derivatives.changes(&nuclide.name);
            let mut k = [0.0; 4];
            let mut responses: [Vec<f64>; 4] = Default::default();
            for (j, step) in steps_of(epsilon).into_iter().enumerate() {
                let mut unique = state.unique.to_vec();
                for (u, n, d) in &changes {
                    unique[*u] = state.unique[*u].added(d, step * n)?;
                }
                let perturbed: Vec<&Material> =
                    state.region_index.iter().map(|&i| &unique[i]).collect();
                let (s, _) = solve_symmetric_from(
                    state.tracking,
                    &perturbed,
                    options,
                    state.symmetry,
                    Some(&converged),
                )
                .with_context(|| format!("{}: direct perturbation {step}", nuclide.name))?;
                k[j] = s.k;
                if let Some(fg) = &few_group {
                    responses[j] = ratios(&response_values(
                        &fg.responses,
                        &fg.group_of,
                        &perturbed,
                        &s,
                    ));
                }
            }
            let steps = scalar_steps(k, epsilon);
            nuclide.direct = Some(Direct {
                k,
                central: steps.central,
                richardson: steps.richardson,
            });
            if let Some(f) = few_group_nuclides.get_mut(position) {
                f.direct = Some(response_steps(&responses, epsilon));
            }
        }
    }
    let direct_seconds = clock.elapsed().as_secs_f64();
    let clock = Instant::now();
    if let (Some(fg), Some(epsilon)) = (&few_group, input.direct) {
        for name in &fg.input.shielded {
            let Some(target) = few_group_nuclides.iter_mut().find(|f| &f.name == name) else {
                bail!("[sensitivity.few-group] shielded {name} is not a [sensitivity] nuclide");
            };
            let reshield = reshield.as_deref_mut().ok_or_else(|| {
                anyhow!(
                    "[sensitivity.few-group] shielded solves the self-shielding of the state again, which needs a \
                     self-shielded state solved from its own composition (not from a reference state)"
                )
            })?;
            target.shielded = Some(
                shielded_steps(
                    name,
                    epsilon,
                    state,
                    &materials,
                    options,
                    (&converged, &importance),
                    fg,
                    &solution,
                    reshield,
                )
                .with_context(|| {
                    format!("{name}: perturbation with the self-shielding solved again")
                })?,
            );
        }
    }
    let shielded_seconds = clock.elapsed().as_secs_f64();
    let generalized = match &input.generalized {
        Some(g) => Some(compare_generalized(
            g,
            state,
            &materials,
            options,
            (&solution, &converged),
            importance.adjoint_start(),
        )?),
        None => None,
    };
    let seconds = importance.seconds;
    tracing::info!(
        solve_seconds,
        derivative_seconds,
        direct_seconds,
        shielded_seconds,
        "k sensitivity"
    );
    let report = Report {
        k: importance.k(),
        adjoint_k: importance.adjoint_k(),
        forward_iterations: solution.iterations,
        adjoint_iterations: importance.adjoint_iterations,
        k_tolerance: options.k_tolerance,
        flux_tolerance: options.flux_tolerance,
        epsilon: input.direct,
        seconds: Seconds {
            forward: seconds[0],
            adjoint: seconds[1],
            products: seconds[2],
            derivatives: derivative_seconds,
            direct: direct_seconds,
        },
        nuclides,
        generalized,
        few_group: few_group.map(|fg| FewGroupReport {
            bounds: fg.input.bounds.clone(),
            responses: fg.responses.iter().map(ToString::to_string).collect(),
            values: ratios(&fg.solve.values),
            table_difference: fg.table_difference,
            chi: fg.chi.clone(),
            tolerance: fg.input.tolerance,
            restart: fg.input.restart,
            sweeps: fg.solve.sweeps.clone(),
            seconds: FewGroupSeconds {
                generalized: fg.seconds,
                derivatives: few_group_derivative_seconds,
                shielded: shielded_seconds,
            },
            nuclides: few_group_nuclides,
        }),
    };
    Ok((solution, converged, report))
}

/// 解いた物質ごとの物質の微分 `weighted`（解いた物質の番号、微分）を、それを指す [領域] の（領域、微分）にしたもの
fn region_changes<'a>(
    weighted: &'a [(usize, Material)],
    region_index: &[usize],
) -> Vec<(usize, &'a Material)> {
    region_index
        .iter()
        .enumerate()
        .filter_map(|(region, i)| {
            weighted
                .iter()
                .find(|(u, _)| u == i)
                .map(|(_, d)| (region, d))
        })
        .collect()
}

/// 少数群定数の応答の値を確かめ、一般化随伴を解く。応答の値が `fewgroup::homogenize` のテーブルの値と丸めの範囲で一致しなければ、
/// 微分がテーブルの定数のものでないので止める
fn few_group_state<'a>(
    input: &'a FewGroupInput,
    state: &State,
    materials: &[&Material],
    options: &Options,
    (solution, forward, adjoint): (&Solution, &Start, &Start),
) -> anyhow::Result<FewGroupState<'a>> {
    let groups = &state.library.groups;
    let group_of =
        coarse_groups(groups, &input.bounds).context("[sensitivity.few-group] bounds")?;
    let responses = table_responses(input.bounds.len() - 1);
    let table = homogenize(
        solution,
        materials,
        groups,
        &input.bounds,
        &state.tracking.boundaries(),
    )
    .context("[sensitivity.few-group] table constants")?;
    let values = ratios(&response_values(&responses, &group_of, materials, solution));
    let mut table_difference = 0.0_f64;
    for (response, value) in responses.iter().zip(&values) {
        let constant = match *response {
            Response::Total(g) => table.total[g],
            Response::NuFission(g) => table.nu_fission[g],
            Response::Fission(g) => table.fission[g],
            Response::Scatter(g, h) => table.scatter[g][h],
        };
        let difference = if constant == 0.0 && *value == 0.0 {
            0.0
        } else {
            ((value - constant) / constant).abs()
        };
        if difference.is_nan() || difference > TABLE_TOLERANCE {
            bail!(
                "response {response} is {value}, but the table constant is {constant}, so its derivative \
                 would not be that of the table"
            );
        }
        table_difference = table_difference.max(difference);
    }
    let settings = GeneralizedOptions {
        method: GeneralizedMethod::Gmres {
            restart: input.restart,
        },
        tolerance: input.tolerance,
        max_iterations: options.max_iterations,
        batch: false,
    };
    let clock = Instant::now();
    let solve = FewGroupSolve::new(
        state.tracking,
        materials,
        options,
        Some(state.symmetry),
        (solution, forward, adjoint),
        &responses,
        &group_of,
        &settings,
    )?;
    let seconds = clock.elapsed().as_secs_f64();
    tracing::info!(seconds, sweeps = ?solve.sweeps, "few-group generalized adjoints");
    Ok(FewGroupState {
        input,
        responses,
        group_of,
        solve,
        table_difference,
        chi: table.chi,
        seconds,
    })
}

/// 核種 `name` の数密度を刻み `steps_of(epsilon)` で変えた組成の自己遮蔽を `reshield` で解き直し、その物質で、少数群定数と k の
/// 1 次の変化（状態の一般化随伴と随伴に掛ける。案 I1）と、直接の摂動（案 I3 の参照）の中心差分を求める
#[allow(clippy::too_many_arguments)]
fn shielded_steps(
    name: &str,
    epsilon: f64,
    state: &State,
    materials: &[&Material],
    options: &Options,
    (converged, importance): (&Start, &crate::moc::Importance),
    fg: &FewGroupState,
    solution: &Solution,
    reshield: &mut Reshield,
) -> anyhow::Result<Shielded> {
    let mut first: [Vec<f64>; 4] = Default::default();
    let mut direct: [Vec<f64>; 4] = Default::default();
    let mut k_first = [0.0; 4];
    let mut k_direct = [0.0; 4];
    for (j, step) in steps_of(epsilon).into_iter().enumerate() {
        let (unique, index) = reshield(&scaled_inputs(state.materials, name, step))?;
        if index.len() != materials.len() || index.iter().any(|&i| i >= unique.len()) {
            bail!(
                "the self-shielding solved again gives {} regions for {}",
                index.len(),
                materials.len()
            );
        }
        let perturbed: Vec<&Material> = index.iter().map(|&i| &unique[i]).collect();
        // 物質の変わった領域（真空は変わらない）
        let changed: Vec<(usize, &Material)> = perturbed
            .iter()
            .enumerate()
            .filter(|&(r, m)| *m != materials[r])
            .map(|(r, &m)| (r, m))
            .collect();
        if changed.iter().any(|&(r, _)| materials[r].is_void()) {
            bail!("the self-shielding solved again changes a void region");
        }
        // 1 次の変化は物質について線形なので、変わった領域の (変えた物質の項) − (状態の物質の項) にする
        let base: Vec<(usize, &Material)> =
            changed.iter().map(|&(r, _)| (r, materials[r])).collect();
        let order = |list: &[(usize, &Material)]| {
            fg.solve
                .first_order(&fg.responses, &fg.group_of, solution, list)
        };
        first[j] = order(&changed)?
            .iter()
            .zip(order(&base)?)
            .map(|(a, b)| a - b)
            .collect();
        k_first[j] = importance.k_derivative(&changed)? - importance.k_derivative(&base)?;
        let (s, _) = solve_symmetric_from(
            state.tracking,
            &perturbed,
            options,
            state.symmetry,
            Some(converged),
        )
        .with_context(|| format!("direct perturbation {step}"))?;
        k_direct[j] = s.k;
        direct[j] = ratios(&response_values(
            &fg.responses,
            &fg.group_of,
            &perturbed,
            &s,
        ));
    }
    Ok(Shielded {
        first_order: response_steps(&first, epsilon),
        direct: response_steps(&direct, epsilon),
        k_first_order: scalar_steps(k_first, epsilon),
        k_direct: scalar_steps(k_direct, epsilon),
    })
}

/// 少数群 `coarse` 個のテーブルの少数群定数の応答（全断面積、νΣf、Σf、P0 の散乱行列 G→H の順）
fn table_responses(coarse: usize) -> Vec<Response> {
    (0..coarse)
        .map(Response::Total)
        .chain((0..coarse).map(Response::NuFission))
        .chain((0..coarse).map(Response::Fission))
        .chain((0..coarse).flat_map(|g| (0..coarse).map(move |h| Response::Scatter(g, h))))
        .collect()
}

/// 応答の少数群（散乱は入射の少数群）
fn few_group_of(response: Response) -> usize {
    match response {
        Response::Total(g)
        | Response::NuFission(g)
        | Response::Fission(g)
        | Response::Scatter(g, _) => g,
    }
}

/// 応答の分子 h₁ の、物質 `m` の細かい群 `g` の値（`group_of[g]` が応答の少数群でなければ 0）
fn numerator(response: Response, m: &Material, g: usize, group_of: &[usize]) -> f64 {
    if group_of[g] != few_group_of(response) {
        return 0.0;
    }
    match response {
        Response::Total(_) => m.total[g],
        Response::NuFission(_) => m.nu_fission[g],
        Response::Fission(_) => m.fission[g],
        Response::Scatter(_, to) => m.scatter[g]
            .iter()
            .zip(group_of)
            .filter(|&(_, &h)| h == to)
            .map(|(s, _)| s)
            .sum(),
    }
}

/// [応答] の（⟨h₁, φ⟩、⟨h₂, φ⟩）。h₂ は応答の少数群の細かい群で 1（真空を含むすべての領域）なので、R = ⟨h₁, φ⟩ / ⟨h₂, φ⟩ は
/// `fewgroup::homogenize` の格子全体の少数群定数（MOC の中性子束と体積の重み。物質は SPH 因子を掛けたもの）
fn response_values(
    responses: &[Response],
    group_of: &[usize],
    materials: &[&Material],
    solution: &Solution,
) -> Vec<(f64, f64)> {
    responses
        .iter()
        .map(|&response| {
            let few = few_group_of(response);
            let (mut a, mut b) = (0.0, 0.0);
            for ((m, phi), v) in materials.iter().zip(&solution.flux).zip(&solution.volumes) {
                for (g, f) in phi.iter().enumerate() {
                    a += v * numerator(response, m, g, group_of) * f;
                    if group_of[g] == few {
                        b += v * f;
                    }
                }
            }
            (a, b)
        })
        .collect()
}

/// 少数群定数の応答の一般化随伴と、それから求める 1 次の変化の係数
struct FewGroupSolve {
    /// [応答] の（⟨h₁, φ⟩、⟨h₂, φ⟩）
    values: Vec<(f64, f64)>,
    /// [応答] の中性子束を通じた項の係数。どの物質でも分子が 0 の応答（3 群より多い少数群の上方散乱など）と、分子が分母に比例する
    /// 応答は、値が中性子束の形によらないので一般化随伴を解かず `None`
    importances: Vec<Option<ResponseImportance>>,
    /// [応答] の一般化随伴のスイープの回数（解かない応答は 0）
    sweeps: Vec<usize>,
}

impl FewGroupSolve {
    /// 応答 `responses` の一般化随伴を `settings` で解き、係数を作る。`forward` と `adjoint` は前向きと随伴の固有値計算の収束した状態
    #[allow(clippy::too_many_arguments)]
    fn new(
        tracking: &Tracking,
        materials: &[&Material],
        options: &Options,
        symmetry: Option<&Detected>,
        (solution, forward, adjoint): (&Solution, &Start, &Start),
        responses: &[Response],
        group_of: &[usize],
        settings: &GeneralizedOptions,
    ) -> anyhow::Result<Self> {
        let values = response_values(responses, group_of, materials, solution);
        let mut importances = Vec::with_capacity(responses.len());
        let mut sweeps = Vec::with_capacity(responses.len());
        for (&response, &(a, b)) in responses.iter().zip(&values) {
            if a == 0.0 {
                importances.push(None);
                sweeps.push(0);
                continue;
            }
            // 源ごとに解いて係数にし、一般化随伴の状態（トラックの端の角度中性子束の全体）を応答の数だけ持たない
            let sources = response_sources(&[response], group_of, materials, solution)?;
            // 分子が分母に比例する応答（どの領域でも同じ断面積）は中性子束の形によらず、源は丸めだけ。中性子束を通じた項は 0
            if flat_source(&sources[0], b) {
                importances.push(None);
                sweeps.push(0);
                continue;
            }
            let gamma = solve_generalized(
                tracking, materials, options, symmetry, forward, adjoint, &sources, settings,
            )
            .with_context(|| format!("generalized adjoint of the few-group constant {response}"))?;
            importances.push(Some(response_importance(
                tracking,
                materials,
                options,
                symmetry,
                forward,
                &gamma[0],
                &sources[0],
            )?));
            sweeps.push(gamma[0].iterations);
        }
        Ok(Self {
            values,
            importances,
            sweeps,
        })
    }

    /// 領域ごとの物質の変化 `changes`（領域、δM）による [応答] の 1 次の変化 δR = R (⟨δh₁, φ⟩ / ⟨h₁, φ⟩ + ⟨Γ, (δF / k − δL) φ⟩)
    /// （分母の重みは物質によらない）。R が 0 の応答は ⟨δh₁, φ⟩ / ⟨h₂, φ⟩
    fn first_order(
        &self,
        responses: &[Response],
        group_of: &[usize],
        solution: &Solution,
        changes: &[(usize, &Material)],
    ) -> anyhow::Result<Vec<f64>> {
        responses
            .iter()
            .zip(&self.values)
            .zip(&self.importances)
            .map(|((&response, &(a, b)), importance)| {
                let mut direct = 0.0;
                for &(region, d) in changes {
                    let (phi, v) = (&solution.flux[region], solution.volumes[region]);
                    for (g, f) in phi.iter().enumerate() {
                        direct += v * numerator(response, d, g, group_of) * f;
                    }
                }
                let Some(importance) = importance else {
                    return Ok(direct / b);
                };
                let flux = importance
                    .flux_term(changes)
                    .with_context(|| format!("response {response}"))?;
                Ok(a / b * (direct / a + flux))
            })
            .collect()
    }
}

/// 比の応答の源 `source`（S = h₁ / a − h₂ / b）が丸めだけか。h₂ / b は 1 / b で、h₁ / a は (h₁ / R) / b なので、どちらの項も
/// 1 / b の桁（R は h₁ の重みの平均）。その 1e-12 倍以下なら、差は足す順の丸め
fn flat_source(source: &[Vec<f64>], b: f64) -> bool {
    let largest = source.iter().flatten().fold(0.0_f64, |m, x| m.max(x.abs()));
    largest <= 1e-12 / b
}

/// 物質の入力 `materials` のうち、核種 `name` を正の数密度で持つ組成の数密度を (1 + `step`) 倍したもの（`Derivatives::changes` と
/// 同じ向き）
fn scaled_inputs(materials: &[MaterialInput], name: &str, step: f64) -> Vec<MaterialInput> {
    materials
        .iter()
        .map(|m| MaterialInput {
            composition: m
                .composition
                .iter()
                .map(|(n, d)| {
                    let scaled = if n == name && *d > 0.0 {
                        d * (1.0 + step)
                    } else {
                        *d
                    };
                    (n.clone(), scaled)
                })
                .collect(),
            ..m.clone()
        })
        .collect()
}

/// 応答ごとの一般化随伴の源 S = h₁ / ⟨h₁, φ⟩ − h₂ / ⟨h₂, φ⟩（[応答][領域][群]。⟨h, φ⟩ = Σ V h φ）。`group_of` は細かい群の少数群
fn response_sources(
    responses: &[Response],
    group_of: &[usize],
    materials: &[&Material],
    solution: &Solution,
) -> anyhow::Result<Vec<Vec<Vec<f64>>>> {
    let ng = group_of.len();
    responses
        .iter()
        .zip(response_values(responses, group_of, materials, solution))
        .map(|(&response, (a, b))| {
            if a == 0.0 || !a.is_finite() || b <= 0.0 {
                bail!("response {response}: <h, phi> is {a} / {b}, so the ratio is not defined");
            }
            let few = few_group_of(response);
            Ok(materials
                .iter()
                .map(|m| {
                    (0..ng)
                        .map(|g| {
                            let denominator = f64::from(u8::from(group_of[g] == few));
                            numerator(response, m, g, group_of) / a - denominator / b
                        })
                        .collect()
                })
                .collect())
        })
        .collect()
}

/// `[sensitivity.generalized]` の比較: 応答の源を、方法ごと（と参照）に `solve_generalized` で解き、時間、スイープの回数、
/// ⟨Γ, Fφ⟩ の残り、参照との差を書く
fn compare_generalized(
    input: &GeneralizedInput,
    state: &State,
    materials: &[&Material],
    options: &Options,
    (solution, forward): (&Solution, &Start),
    adjoint: &Start,
) -> anyhow::Result<GeneralizedReport> {
    let group_of = coarse_groups(&state.library.groups, &input.bounds)
        .context("[sensitivity.generalized] bounds")?;
    let sources = response_sources(&input.responses, &group_of, materials, solution)?;
    let run = |method: Method,
               tolerance: f64,
               batch: bool|
     -> anyhow::Result<(GeneralizedRun, Vec<Generalized>)> {
        let settings = GeneralizedOptions {
            method: match method {
                Method::SourceIteration => GeneralizedMethod::SourceIteration,
                Method::Cmfd => GeneralizedMethod::Cmfd(generalized_cmfd(&options.acceleration)?),
                Method::Gmres { restart } => GeneralizedMethod::Gmres { restart },
            },
            tolerance,
            max_iterations: options.max_iterations,
            batch,
        };
        let clock = Instant::now();
        let solved = solve_generalized(
            state.tracking,
            materials,
            options,
            Some(state.symmetry),
            forward,
            adjoint,
            &sources,
            &settings,
        );
        let seconds = clock.elapsed().as_secs_f64();
        // 発散と打ち切りは方法の結果として書き、ほかの方法の比較を続ける（参照の失敗は下で止める）
        let gamma = match solved {
            Ok(gamma) => gamma,
            Err(
                e @ (MocError::GeneralizedDiverged { .. }
                | MocError::GeneralizedNotConverged { .. }
                | MocError::Cmfd(_)),
            ) => {
                tracing::warn!(%method, tolerance, batch, seconds, error = %e, "generalized adjoints failed");
                let failed = GeneralizedRun {
                    method: method.to_string(),
                    tolerance,
                    batch,
                    seconds,
                    sweeps: Vec::new(),
                    orthogonality: Vec::new(),
                    difference: None,
                    failure: Some(e.to_string()),
                };
                return Ok((failed, Vec::new()));
            }
            Err(e) => {
                return Err(
                    anyhow::Error::new(e).context(format!("generalized adjoint by {method}"))
                );
            }
        };
        let orthogonality = gamma
            .iter()
            .map(|g| fission_overlap(materials, solution, &g.flux))
            .collect();
        tracing::info!(%method, tolerance, batch, seconds, "generalized adjoints");
        Ok((
            GeneralizedRun {
                method: method.to_string(),
                tolerance,
                batch,
                seconds,
                sweeps: gamma.iter().map(|g| g.iterations).collect(),
                orthogonality,
                difference: None,
                failure: None,
            },
            gamma,
        ))
    };
    let reference = input
        .reference
        .map(|(method, tolerance)| run(method, tolerance, false))
        .transpose()?;
    if let Some((
        GeneralizedRun {
            failure: Some(f),
            method,
            ..
        },
        _,
    )) = &reference
    {
        bail!("the reference generalized adjoint by {method} failed: {f}");
    }
    let mut runs = Vec::new();
    let cases = input.methods.iter().flat_map(|&method| {
        let batches = if input.batched.contains(&method) {
            &[false, true][..]
        } else {
            &[false][..]
        };
        input
            .tolerances
            .iter()
            .flat_map(move |&t| batches.iter().map(move |&b| (method, t, b)))
    });
    for (method, tolerance, batch) in cases {
        let (mut result, gamma) = run(method, tolerance, batch)?;
        if let Some((_, reference)) = &reference
            && result.failure.is_none()
        {
            result.difference = Some(
                gamma
                    .iter()
                    .zip(reference)
                    .map(|(a, b)| largest_difference(&a.flux, &b.flux))
                    .collect(),
            );
        }
        runs.push(result);
    }
    Ok(GeneralizedReport {
        responses: input.responses.iter().map(ToString::to_string).collect(),
        runs,
        reference: reference.map(|(r, _)| r),
    })
}

/// 一般化随伴の CMFD の判定: 状態を解いた MOC の CMFD のもの。状態を CMFD なしで解いたなら、入力が粗メッシュの判定を決めていないので
/// エラー
fn generalized_cmfd(acceleration: &Acceleration) -> anyhow::Result<CmfdOptions> {
    match acceleration {
        Acceleration::Cmfd(c) => Ok(*c),
        Acceleration::None => bail!(
            "[sensitivity.generalized] cmfd takes the CMFD tolerances of the state, which was solved without CMFD"
        ),
    }
}

/// ⟨Γ, Fφ⟩ / ⟨|Γ|, Fφ⟩（`flux` は [領域][群] の Γ）
fn fission_overlap(materials: &[&Material], solution: &Solution, flux: &[Vec<f64>]) -> f64 {
    let (mut sum, mut size) = (0.0, 0.0);
    for (r, m) in materials.iter().enumerate() {
        for t in &m.fission_terms {
            let rate: f64 = t
                .nu_fission
                .iter()
                .zip(&solution.flux[r])
                .map(|(n, x)| n * x)
                .sum();
            let v = solution.volumes[r] * rate;
            sum += v * t.chi.iter().zip(&flux[r]).map(|(c, x)| c * x).sum::<f64>();
            size += v * t
                .chi
                .iter()
                .zip(&flux[r])
                .map(|(c, x)| c * x.abs())
                .sum::<f64>();
        }
    }
    if size == 0.0 { 0.0 } else { sum / size }
}

/// 領域と群の |a − b| の最大を、|b| の最大で割ったもの（有限でない値があれば NaN。`f64::max` は NaN を捨てるので先に見る）
fn largest_difference(a: &[Vec<f64>], b: &[Vec<f64>]) -> f64 {
    let values = || a.iter().flatten().zip(b.iter().flatten());
    if values().any(|(x, y)| !(x.is_finite() && y.is_finite())) {
        return f64::NAN;
    }
    let (mut change, mut scale) = (0.0_f64, 0.0_f64);
    for (x, y) in values() {
        change = change.max((x - y).abs());
        scale = scale.max(y.abs());
    }
    change / scale
}

/// 解いた物質ごとの、核種の数密度の微分（数密度 1 の群定数。自己遮蔽した物質では、その領域の実効断面積と移行の群定数に
/// SPH 因子を掛けたもの）
struct Derivatives {
    /// [解いた物質] の（核種の名前、数密度、微分）。真空と、組成に核種のない物質は空
    materials: Vec<Vec<(String, f64, Material)>>,
}

impl Derivatives {
    /// 解いた物質ごとに、それを指す最初の領域の物質の入力と自己遮蔽から微分を作り、組成の数密度で足した和が解いた物質に
    /// （丸めの範囲で）一致することを確かめる。一致しなければ、物質の作り方がこの微分の前提（核種の群定数の数密度の重みの和に
    /// SPH 因子を掛けたもの）と違うので止める
    fn new(state: &State) -> anyhow::Result<Self> {
        let ng = state.library.groups.len();
        let sph = match state.shielding {
            Some(s) => {
                Some(s.physical_flux(state.lattice, &vec![vec![1.0; ng]; state.lattice.regions()])?)
            }
            None => None,
        };
        let mut materials = Vec::with_capacity(state.unique.len());
        for (u, solved) in state.unique.iter().enumerate() {
            // どの領域も指さない物質（自己遮蔽の前の入力の物質など）は解いた問題に入らない
            let Some(region) = state.region_index.iter().position(|&i| i == u) else {
                materials.push(Vec::new());
                continue;
            };
            let name = &state.region_names[region];
            if name == VOID {
                materials.push(Vec::new());
                continue;
            }
            let input = state
                .materials
                .iter()
                .find(|m| &m.name == name)
                .ok_or_else(|| anyhow!("material `{name}` is not a material of nuclides"))?;
            let effective = state
                .shielding
                .map(|s| (s.region_cross_sections(region), s.region_transfers(region)));
            let mut list = Vec::with_capacity(input.composition.len());
            for (k, (nuclide, density)) in input.composition.iter().enumerate() {
                let d = match effective {
                    Some((cross_sections, transfers)) => {
                        let shielded = cross_sections
                            .and_then(|c| c.get(k))
                            .ok_or_else(|| anyhow!("material `{name}`: no shielded {nuclide}"))?;
                        let data = shielded_nuclide_with(
                            state.library,
                            nuclide,
                            input.temperature,
                            shielded,
                            transfers.and_then(|t| t.get(k)),
                        )?;
                        let mu = &sph
                            .as_ref()
                            .map_or_else(|| vec![1.0; ng], |s| s[region].clone());
                        Material::from_nuclides(&[(data, 1.0)], ng)?.scaled(mu)?
                    }
                    None => {
                        Material::mix(state.library, &[(nuclide.clone(), 1.0)], input.temperature)?
                    }
                };
                list.push((nuclide.clone(), *density, d));
            }
            check_sum(name, solved, &list)?;
            materials.push(list);
        }
        Ok(Self { materials })
    }

    /// 核種 `name` をどれかの物質の組成に持つか（数密度 0 を含む）
    fn lists(&self, name: &str) -> bool {
        self.materials.iter().flatten().any(|(n, _, _)| n == name)
    }

    /// 核種 `name` を正の数密度で持つ解いた物質の（番号、数密度、微分）
    fn changes(&self, name: &str) -> Vec<(usize, f64, &Material)> {
        self.materials
            .iter()
            .enumerate()
            .flat_map(|(u, list)| {
                list.iter()
                    .filter(|(n, density, _)| n == name && *density > 0.0)
                    .map(move |(_, density, d)| (u, *density, d))
            })
            .collect()
    }
}

/// 物質の値をすべて絶対値にする（`check_sum` の項の大きさ。核分裂の項は νΣf だけを絶対値にし、χ は非負なので核分裂の行列の
/// 要素も項の絶対値の和になる）
fn absolute(m: &mut Material) {
    for v in std::iter::once(&mut m.total)
        .chain(std::iter::once(&mut m.fission))
        .chain(std::iter::once(&mut m.nu_fission))
        .chain(m.fission_terms.iter_mut().map(|t| &mut t.nu_fission))
        .chain(m.scatter.iter_mut())
        .chain(m.anisotropic.iter_mut().flatten())
    {
        for x in v.iter_mut() {
            *x = x.abs();
        }
    }
}

/// 微分 `list` を数密度で足した物質が、解いた物質 `solved` と、項の絶対値の和に比べて丸めの範囲で一致するか
fn check_sum(
    name: &str,
    solved: &Material,
    list: &[(String, f64, Material)],
) -> anyhow::Result<()> {
    let Some((first, rest)) = list.split_first() else {
        bail!("material `{name}` has no nuclides");
    };
    let mut sum = first.2.scaled(&vec![first.1; first.2.groups()])?;
    // 項の絶対値の和（核種の寄与が打ち消し合う要素では、足す順の丸めの差は和でなく項の大きさで決まる）
    let mut size = first.2.scaled(&vec![first.1; first.2.groups()])?;
    absolute(&mut size);
    for (_, density, d) in rest {
        sum = sum.added(d, *density)?;
        let mut term = d.scaled(&vec![*density; d.groups()])?;
        absolute(&mut term);
        size = size.added(&term, 1.0)?;
    }
    let close = |a: &[f64], b: &[f64], scale: &[f64]| {
        a.iter()
            .zip(b)
            .zip(scale)
            .all(|((x, y), s)| (x - y).abs() <= RECONSTRUCTION_TOLERANCE * s.max(f64::MIN_POSITIVE))
    };
    let rows = |m: &Material| -> Vec<f64> {
        std::iter::once(&m.scatter)
            .chain(&m.anisotropic)
            .flatten()
            .flatten()
            .copied()
            .collect()
    };
    let terms = |m: &Material| -> Vec<f64> { m.fission_matrix().into_iter().flatten().collect() };
    if sum.anisotropic.len() != solved.anisotropic.len()
        || !close(&sum.total, &solved.total, &size.total)
        || !close(&rows(&sum), &rows(solved), &rows(&size))
        || !close(&sum.nu_fission, &solved.nu_fission, &size.nu_fission)
        || !close(&terms(&sum), &terms(solved), &terms(&size))
    {
        bail!(
            "material `{name}`: the nuclide cross sections weighted by their densities do not add up to the solved \
             material, so its density derivatives would not be those of the solved problem"
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests;
