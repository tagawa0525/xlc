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
use crate::geometry::Lattice;
use crate::input::{MaterialInput, VOID};
use crate::material::Material;
use crate::moc::{Acceleration, Options, Solution, Start, solve_importance, solve_symmetric_from};
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
    };
    Ok((solution, converged, report))
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

/// 微分 `list` を数密度で足した物質が、解いた物質 `solved` と一致するか
fn check_sum(
    name: &str,
    solved: &Material,
    list: &[(String, f64, Material)],
) -> anyhow::Result<()> {
    let Some((first, rest)) = list.split_first() else {
        bail!("material `{name}` has no nuclides");
    };
    let mut sum = first.2.scaled(&vec![first.1; first.2.groups()])?;
    for (_, density, d) in rest {
        sum = sum.added(d, *density)?;
    }
    let close = |a: &[f64], b: &[f64]| {
        a.iter().zip(b).all(|(x, y)| {
            (x - y).abs() <= RECONSTRUCTION_TOLERANCE * x.abs().max(y.abs()).max(f64::MIN_POSITIVE)
        })
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
        || !close(&sum.total, &solved.total)
        || !close(&rows(&sum), &rows(solved))
        || !close(&sum.nu_fission, &solved.nu_fission)
        || !close(&terms(&sum), &terms(solved))
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
