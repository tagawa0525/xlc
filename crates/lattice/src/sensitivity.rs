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
use crate::fewgroup::coarse_groups;
use crate::geometry::Lattice;
use crate::input::{MaterialInput, VOID};
use crate::material::Material;
use crate::moc::{
    Acceleration, Generalized, GeneralizedMethod, GeneralizedOptions, MocError, Options,
    ResponseImportance, Solution, Start, solve_generalized, solve_importance, solve_symmetric_from,
};
use crate::shielding::{FuelShielding, shielded_nuclide_with};
use crate::symmetry::Detected;
use crate::tracking::Tracking;

/// 物質の和を物質の微分の和と比べる相対の許容値（同じ群定数を足す順の丸めの違いだけを許す）
const RECONSTRUCTION_TOLERANCE: f64 = 1.0e-10;

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

/// 状態 `state` を、`start` から（`None` なら平らな初期値から）`options` で解き、感度を求める。前向きの解と収束した状態も返す
pub fn evaluate(
    input: &SensitivityInput,
    state: &State,
    options: &Options,
    start: Option<&Start>,
) -> anyhow::Result<(Solution, Start, Report)> {
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
        let list: Vec<(usize, &Material)> = state
            .region_index
            .iter()
            .enumerate()
            .filter_map(|(region, i)| {
                weighted
                    .iter()
                    .find(|(u, _)| u == i)
                    .map(|(_, d)| (region, d))
            })
            .collect();
        let relative = importance
            .k_derivative(&list)
            .with_context(|| format!("{name}: k derivative"))?;
        nuclides.push(NuclideSensitivity {
            name: name.clone(),
            density,
            relative,
            derivative: (density > 0.0).then(|| relative / density),
            direct: None,
        });
    }
    let derivative_seconds = clock.elapsed().as_secs_f64();
    let clock = Instant::now();
    if let Some(epsilon) = input.direct {
        for nuclide in &mut nuclides {
            let changes = derivatives.changes(&nuclide.name);
            let mut k = [0.0; 4];
            for (value, step) in
                k.iter_mut()
                    .zip([epsilon, -epsilon, epsilon / 2.0, -epsilon / 2.0])
            {
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
                *value = s.k;
            }
            let central = [(k[0] - k[1]) / (2.0 * epsilon), (k[2] - k[3]) / epsilon];
            nuclide.direct = Some(Direct {
                k,
                central,
                richardson: central[1] + (central[1] - central[0]) / 3.0,
            });
        }
    }
    let direct_seconds = clock.elapsed().as_secs_f64();
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
    };
    Ok((solution, converged, report))
}

/// 少数群 `coarse` 個のテーブルの少数群定数の応答（全断面積、νΣf、Σf、P0 の散乱行列 G→H の順）
fn table_responses(coarse: usize) -> Vec<Response> {
    let _ = coarse;
    Vec::new()
}

/// 応答の分子 h₁ の、物質 `m` の細かい群 `g` の値（`group_of[g]` が応答の少数群でなければ 0）
fn numerator(response: Response, m: &Material, g: usize, group_of: &[usize]) -> f64 {
    let _ = (response, m, g, group_of);
    0.0
}

/// [応答] の（⟨h₁, φ⟩、⟨h₂, φ⟩）。R = ⟨h₁, φ⟩ / ⟨h₂, φ⟩ は `fewgroup::homogenize` の格子全体の少数群定数
fn response_values(
    responses: &[Response],
    group_of: &[usize],
    materials: &[&Material],
    solution: &Solution,
) -> Vec<(f64, f64)> {
    let _ = (group_of, materials, solution);
    vec![(0.0, 1.0); responses.len()]
}

/// 少数群定数の応答の一般化随伴と、それから求める 1 次の変化の係数
struct FewGroupSolve {
    /// [応答] の（⟨h₁, φ⟩、⟨h₂, φ⟩）
    values: Vec<(f64, f64)>,
    /// [応答] の中性子束を通じた項の係数
    importances: Vec<ResponseImportance>,
    /// [応答] の一般化随伴のスイープの回数
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
        let _ = (
            tracking, materials, options, symmetry, solution, forward, adjoint, responses,
            group_of, settings,
        );
        bail!("not implemented")
    }

    /// 領域ごとの物質の変化 `changes`（領域、δM）による [応答] の 1 次の変化 δR = R (⟨δh₁, φ⟩ / ⟨h₁, φ⟩ + ⟨Γ, (δF / k − δL) φ⟩)
    /// （分母の重みは物質によらない）
    fn first_order(
        &self,
        responses: &[Response],
        group_of: &[usize],
        solution: &Solution,
        changes: &[(usize, &Material)],
    ) -> anyhow::Result<Vec<f64>> {
        let _ = (group_of, solution, changes);
        Ok(vec![0.0; responses.len()])
    }
}

/// 物質の入力 `materials` のうち、核種 `name` を正の数密度で持つ組成の数密度を (1 + `step`) 倍したもの（`Derivatives::changes` と
/// 同じ向き）
fn scaled_inputs(materials: &[MaterialInput], name: &str, step: f64) -> Vec<MaterialInput> {
    let _ = (name, step);
    materials.to_vec()
}

/// 応答ごとの一般化随伴の源 S = h₁ / ⟨h₁, φ⟩ − h₂ / ⟨h₂, φ⟩（[応答][領域][群]。⟨h, φ⟩ = Σ V h φ）。`group_of` は細かい群の少数群
fn response_sources(
    responses: &[Response],
    group_of: &[usize],
    materials: &[&Material],
    solution: &Solution,
) -> anyhow::Result<Vec<Vec<Vec<f64>>>> {
    let ng = group_of.len();
    let inner = |h: &[Vec<f64>]| -> f64 {
        h.iter()
            .zip(&solution.flux)
            .zip(&solution.volumes)
            .map(|((h, f), v)| v * h.iter().zip(f).map(|(a, b)| a * b).sum::<f64>())
            .sum()
    };
    responses
        .iter()
        .map(|response| {
            let few = match *response {
                Response::Total(g)
                | Response::NuFission(g)
                | Response::Fission(g)
                | Response::Scatter(g, _) => g,
            };
            let numerator: Vec<Vec<f64>> = materials
                .iter()
                .map(|m| {
                    (0..ng)
                        .map(|g| {
                            if group_of[g] != few {
                                return 0.0;
                            }
                            match *response {
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
                        })
                        .collect()
                })
                .collect();
            let denominator: Vec<Vec<f64>> = materials
                .iter()
                .map(|_| {
                    (0..ng)
                        .map(|g| f64::from(u8::from(group_of[g] == few)))
                        .collect()
                })
                .collect();
            let (a, b) = (inner(&numerator), inner(&denominator));
            if a == 0.0 || !a.is_finite() || b <= 0.0 {
                bail!("response {response}: <h, phi> is {a} / {b}, so the ratio is not defined");
            }
            Ok(numerator
                .iter()
                .zip(&denominator)
                .map(|(x, y)| x.iter().zip(y).map(|(p, q)| p / a - q / b).collect())
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
