//! 3次元直交座標の多群ノード展開法（NEM）。
//!
//! 各ノード・各方向で横方向に積分した1次元の拡散方程式を、ノード平均と4次までの多項式
//! f1 = ξ, f2 = 3ξ² - 1/4, f3 = ξ³ - ξ/4, f4 = (ξ² - 1/20)(ξ² - 1/4)（ξ ∈ [-1/2, 1/2]）で解く。
//! 3次・4次の係数は f1, f2 で重み付けした残差から、横方向の漏れは隣のノードを含む
//! 2次式で近似する。ノード間は部分中性子流でつなぐ。
//!
//! 面では、一般化等価理論（GET。K. S. Smith, "Assembly homogenization techniques for light water
//! reactor analysis", Prog. Nucl. Energy 17, 303 (1986)）の条件で、不連続因子（ADF）f をかけた
//! 不均質な中性子束 f φ と正味の中性子流を連続にする。径方向の ADF は少数群定数の辺ごとの値で、
//! 2次元の格子計算からは軸方向の ADF は出ない（軸方向に一様な格子では 1）ので、z の面では 1 とする。

use xlc_format::fewgroup::FewGroup;
use xlc_numeric::nem::{Beside, Inflow, Line, Outflow, outflow, transverse_moments};

/// 外周の境界条件。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// 正味の中性子流 0（対称面）
    Reflective,
    /// 境界での中性子束 0
    ZeroFlux,
    /// 外から入る部分中性子流 0（Marshak の真空境界）
    Vacuum,
}

/// 直方体の炉心のノード分割。幅は cm。
#[derive(Debug, Clone, PartialEq)]
pub struct Mesh {
    pub dx: Vec<f64>,
    pub dy: Vec<f64>,
    pub dz: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Problem {
    pub mesh: Mesh,
    pub materials: Vec<FewGroup>,
    /// ノードごとの物質番号。並びは x が最も速く、次に y、z
    pub node_material: Vec<usize>,
    /// [-x, +x, -y, +y, -z, +z]
    pub boundaries: [Boundary; 6],
}

/// 反復の打ち切り条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Options {
    pub k_tolerance: f64,
    pub flux_tolerance: f64,
    pub max_iterations: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            k_tolerance: 1.0e-8,
            flux_tolerance: 1.0e-7,
            max_iterations: 10_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    pub k: f64,
    /// [ノード][群] のノード平均中性子束（全体の生成率が 1 になるよう規格化）
    pub flux: Vec<Vec<f64>>,
    /// 外周から外へ出る正味の中性子数（群ごと）
    pub leakage: Vec<f64>,
    pub iterations: usize,
}

/// NEM の固有値計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum NemError {
    /// 問題の設定（メッシュ、物質の配置、断面積）の誤り
    #[error("{0}")]
    InvalidProblem(String),
    #[error("no fissile material")]
    NoFission,
    #[error(
        "node {node} has a non-positive flux; the iteration converged to a \
         non-fundamental mode (use smaller nodes)"
    )]
    NonFundamentalMode { node: usize },
    /// 真空の面の横方向の漏れに、面の中性子束 / ノード平均を使えない
    #[error("node {node} next to a vacuum face has a flux {flux}; it must be positive")]
    VacuumFaceFlux { node: usize, flux: f64 },
    #[error("not converged in {iterations} iterations (k = {k})")]
    NotConverged { iterations: usize, k: f64 },
}

#[tracing::instrument(name = "nem", skip_all)]
pub fn solve(problem: &Problem, options: &Options) -> Result<Solution, NemError> {
    let core = Core::new(problem).map_err(NemError::InvalidProblem)?;
    iterate(&core, options).map(|(solution, _)| solution)
}

/// べき乗法の反復。収束した解と、そのときの部分中性子流などの状態を返す。
fn iterate(core: &Core, options: &Options) -> Result<(Solution, State), NemError> {
    let ng = core.groups;
    let nodes = core.nodes();
    let mut state = State::new(nodes, ng);
    let mut k = 1.0;
    let production = |state: &State| -> f64 {
        (0..nodes)
            .map(|n| {
                let m = core.material(n);
                core.volume(n)
                    * (0..ng)
                        .map(|g| m.nu_fission[g] * state.flux[n][g])
                        .sum::<f64>()
            })
            .sum()
    };
    let initial = production(&state);
    if initial <= 0.0 {
        return Err(NemError::NoFission);
    }
    state.scale(1.0 / initial);
    let mut k_history = Remaining::default();
    let mut flux_history = Remaining::default();

    for iteration in 1..=options.max_iterations {
        let previous = state.flux.clone();
        for g in 0..ng {
            core.sweep_group(&mut state, g, k)?;
        }
        // べき乗法: 生成率の比で k を更新し、全体を生成率 1 に規格化する
        let rate = production(&state);
        let new_k = k * rate;
        state.scale(1.0 / rate);
        let change = state
            .flux
            .iter()
            .flatten()
            .zip(previous.iter().flatten())
            .map(|(new, old)| relative_change(*new, *old))
            .fold(0.0, f64::max);
        tracing::debug!(
            iteration,
            k = new_k,
            k_change = new_k - k,
            flux_change = change
        );
        let k_remaining = k_history.update((new_k - k).abs());
        let flux_remaining = flux_history.update(change);
        k = new_k;
        if k_remaining < options.k_tolerance && flux_remaining < options.flux_tolerance {
            // 基本モードの中性子束はどこでも正。正でないノードがあれば、別のモードに
            // 収束している（横方向の漏れが除去に比べて非常に大きい粗いノードで起きる）
            if let Some(n) = state
                .flux
                .iter()
                .position(|phi| phi.iter().any(|&f| f <= 0.0))
            {
                return Err(NemError::NonFundamentalMode { node: n });
            }
            let solution = Solution {
                k,
                leakage: core.leakage(&state),
                flux: state.flux.clone(),
                iterations: iteration,
            };
            return Ok((solution, state));
        }
    }
    Err(NemError::NotConverged {
        iterations: options.max_iterations,
        k,
    })
}

/// 検査済みの問題と、ノードの位置の計算。
struct Core<'a> {
    problem: &'a Problem,
    /// 方向ごとのノード幅
    widths: [&'a [f64]; 3],
    counts: [usize; 3],
    groups: usize,
}

/// 反復で更新する量。
struct State {
    /// [ノード][群]
    flux: Vec<Vec<f64>>,
    /// [方向][ノード][群] の展開係数 a1〜a4
    coefficients: [Vec<Vec<[f64; 4]>>; 3],
    /// [方向][ノード][群] の、ノードの - 面・+ 面から出ていく部分中性子流
    out_minus: [Vec<Vec<f64>>; 3],
    out_plus: [Vec<Vec<f64>>; 3],
}

impl State {
    fn new(nodes: usize, ng: usize) -> Self {
        let field = |value: f64| vec![vec![value; ng]; nodes];
        Self {
            flux: field(1.0),
            coefficients: std::array::from_fn(|_| vec![vec![[0.0; 4]; ng]; nodes]),
            // 平坦な中性子束 φ では、部分中性子流は φ/4
            out_minus: std::array::from_fn(|_| field(0.25)),
            out_plus: std::array::from_fn(|_| field(0.25)),
        }
    }

    /// 方程式は線形なので、すべての量を同じ倍率で変えてよい
    fn scale(&mut self, factor: f64) {
        for v in self.flux.iter_mut().flatten() {
            *v *= factor;
        }
        for d in 0..3 {
            for c in self.coefficients[d].iter_mut().flatten().flatten() {
                *c *= factor;
            }
            for v in self.out_minus[d].iter_mut().flatten() {
                *v *= factor;
            }
            for v in self.out_plus[d].iter_mut().flatten() {
                *v *= factor;
            }
        }
    }
}

impl<'a> Core<'a> {
    fn new(problem: &'a Problem) -> Result<Self, String> {
        let mesh = &problem.mesh;
        let widths = [mesh.dx.as_slice(), mesh.dy.as_slice(), mesh.dz.as_slice()];
        if widths
            .iter()
            .any(|w| w.is_empty() || w.iter().any(|&h| !(h > 0.0 && h.is_finite())))
        {
            return Err("node widths must be positive".into());
        }
        let counts = widths.map(<[f64]>::len);
        let nodes: usize = counts.iter().product();
        if problem.node_material.len() != nodes {
            return Err(format!(
                "{} node materials given for {nodes} nodes",
                problem.node_material.len()
            ));
        }
        if let Some(&m) = problem
            .node_material
            .iter()
            .find(|&&m| m >= problem.materials.len())
        {
            return Err(format!("material {m} does not exist"));
        }
        let groups = problem.materials.first().ok_or("no materials")?.groups();
        for (i, m) in problem.materials.iter().enumerate() {
            // 長さ、有限であること、不連続因子は少数群定数の検査（ファイルの読み込みと同じ）
            m.validate().map_err(|e| format!("material {i}: {e}"))?;
            if m.groups() != groups {
                return Err(format!("material {i} must have {groups} groups"));
            }
            if m.diffusion.iter().any(|&d| d.is_nan() || d <= 0.0) {
                return Err(format!(
                    "material {i}: diffusion coefficients must be positive"
                ));
            }
            if (0..groups).any(|g| {
                let removal = m.total[g] - m.scatter[g][g];
                removal.is_nan() || removal <= 0.0
            }) {
                return Err(format!(
                    "material {i}: removal cross sections must be positive"
                ));
            }
        }
        Ok(Self {
            problem,
            widths,
            counts,
            groups,
        })
    }

    fn nodes(&self) -> usize {
        self.counts.iter().product()
    }

    fn position(&self, n: usize) -> [usize; 3] {
        let [nx, ny, _] = self.counts;
        [n % nx, (n / nx) % ny, n / (nx * ny)]
    }

    fn index(&self, p: [usize; 3]) -> usize {
        let [nx, ny, _] = self.counts;
        p[0] + nx * (p[1] + ny * p[2])
    }

    /// 方向 d の隣のノード（`step` が -1 なら - 側、+1 なら + 側）。外周なら `None`。
    fn neighbor(&self, n: usize, d: usize, step: isize) -> Option<usize> {
        let mut p = self.position(n);
        let moved = p[d]
            .checked_add_signed(step)
            .filter(|&i| i < self.counts[d])?;
        p[d] = moved;
        Some(self.index(p))
    }

    fn width(&self, n: usize, d: usize) -> f64 {
        self.widths[d][self.position(n)[d]]
    }

    fn volume(&self, n: usize) -> f64 {
        (0..3).map(|d| self.width(n, d)).product()
    }

    fn material(&self, n: usize) -> &FewGroup {
        &self.problem.materials[self.problem.node_material[n]]
    }

    /// ノード n の方向 d の - 面（`plus` が false）・+ 面（true）の不連続因子。
    /// 少数群定数の辺は（左、右、上、下）で、-x が左、+x が右、-y が下、+y が上。z の面は 1
    fn adf(&self, n: usize, d: usize, plus: bool, g: usize) -> f64 {
        let side = match (d, plus) {
            (0, false) => 0,
            (0, true) => 1,
            (1, true) => 2,
            (1, false) => 3,
            _ => return 1.0,
        };
        self.material(n).adf[side][g]
    }

    /// ノード n の - 面・+ 面に入ってくる部分中性子流。
    fn incoming(&self, state: &State, n: usize, d: usize, g: usize) -> (f64, f64) {
        let ((beta_l, ext_l), (beta_r, ext_r)) = self.face_inputs(state, n, d, g);
        (
            beta_l * state.out_minus[d][n][g] + ext_l,
            beta_r * state.out_plus[d][n][g] + ext_r,
        )
    }

    /// ノード n の - 面・+ 面の入射の条件 (β, 外からの入射)。入射は β × 出ていく部分中性子流
    /// + 外からの入射。
    ///
    /// 隣のノードとの面では、自分の側の ADF を f、隣の側を f'、出ていく部分中性子流を o、o'
    /// とすると、f(o + i) = f'(o' + i')（不均質な中性子束）と o - i = i' - o'（正味の中性子流）から
    ///   i = ((f' - f) o + 2 f' o') / (f + f')
    /// なので、β = (f' - f)/(f + f')、外からの入射は 2 f' o'/(f + f')（ADF が等しければ i = o'）。
    ///
    /// 外周の境界条件も GET のとおり不均質な中性子束 f φ = 2f(o + i) と正味の中性子流 o - i に課し、
    /// 外からの入射は 0:
    /// - 反射: 正味の中性子流 0 で、β = 1（f によらない）
    /// - 中性子束ゼロ: f φ = 0 は φ = 0 と同じで、β = -1（f によらない）
    /// - 真空（Marshak）: 不均質な入射の部分中性子流 f φ/4 - (o - i)/2 = (f(o + i) - (o - i))/2 = 0
    ///   から β = (1 - f)/(1 + f)。隣の面の式で f' = 1、o' = 0 としたものに等しく、f = 1 なら 0
    fn face_inputs(&self, state: &State, n: usize, d: usize, g: usize) -> ((f64, f64), (f64, f64)) {
        let face = |step: isize, index: usize| match self.neighbor(n, d, step) {
            Some(m) => {
                let plus = step > 0;
                let (own, other) = (self.adf(n, d, plus, g), self.adf(m, d, !plus, g));
                let out = if plus {
                    state.out_minus[d][m][g]
                } else {
                    state.out_plus[d][m][g]
                };
                (
                    (other - own) / (own + other),
                    2.0 * other * out / (own + other),
                )
            }
            None => match self.problem.boundaries[2 * d + index] {
                Boundary::Reflective => (1.0, 0.0),
                Boundary::ZeroFlux => (-1.0, 0.0),
                Boundary::Vacuum => {
                    let f = self.adf(n, d, step > 0, g);
                    ((1.0 - f) / (1.0 + f), 0.0)
                }
            },
        };
        (face(-1, 0), face(1, 1))
    }

    /// ノード n の方向 d の正味の漏れ (J+ - J-)/h（単位体積あたり）
    fn net_leakage(&self, state: &State, n: usize, d: usize, g: usize) -> f64 {
        let (in_minus, in_plus) = self.incoming(state, n, d, g);
        let j_plus = state.out_plus[d][n][g] - in_plus;
        let j_minus = in_minus - state.out_minus[d][n][g];
        (j_plus - j_minus) / self.width(n, d)
    }

    fn sweep_group(&self, state: &mut State, g: usize, k: f64) -> Result<(), NemError> {
        let nodes = self.nodes();
        // 横方向の漏れは群のスイープの始めの値を使う
        let leakage: Vec<[f64; 3]> = (0..nodes)
            .map(|n| std::array::from_fn(|d| self.net_leakage(state, n, d, g)))
            .collect();
        for n in 0..nodes {
            let m = self.material(n);
            let removal = m.total[g] - m.scatter[g][g];
            // 群 g への源の結合係数 c(g'→g)（自群の核分裂を含み、自群の散乱を除く）
            let coupling: Vec<f64> = (0..self.groups)
                .map(|from| {
                    let scatter = if from == g { 0.0 } else { m.scatter[from][g] };
                    scatter + m.chi[g] * m.nu_fission[from] / k
                })
                .collect();
            let source: f64 = coupling
                .iter()
                .zip(&state.flux[n])
                .map(|(c, f)| c * f)
                .sum();

            // 出ていく部分中性子流はノード平均について線形（o = p + q φ̄）なので、
            // ノードのバランス式と同時に解く（古い値を使うと高速群で振動する）
            let mut parts = [Outflow::default(); 3];
            let mut balance_rhs = source;
            let mut balance_diag = removal;
            for (d, part) in parts.iter_mut().enumerate() {
                let h = self.width(n, d);
                // ノード n の群 g の値はまだ更新していないので、スイープの始めの値と同じ
                let transverse = self.transverse_shape(state, &leakage, n, d, g)?;
                // 入射は i = β o + 外からの入射（face_inputs）
                let ((beta_l, ext_l), (beta_r, ext_r)) = self.face_inputs(state, n, d, g);
                *part = outflow(&Line {
                    diffusion: m.diffusion[g],
                    width: h,
                    removal,
                    coupling: &coupling,
                    coefficients: &state.coefficients[d][n],
                    group: g,
                    transverse,
                    inflow: [
                        Inflow {
                            beta: beta_r,
                            external: ext_r,
                        },
                        Inflow {
                            beta: beta_l,
                            external: ext_l,
                        },
                    ],
                });
                // 正味の漏れ (o_R + o_L - i_R - i_L)/h（i = β o + 外からの入射）をバランス式に入れる
                for f in 0..2 {
                    balance_rhs -= ((1.0 - part.beta[f]) * part.out[f] - part.external[f]) / h;
                    balance_diag += (1.0 - part.beta[f]) * part.slope[f] / h;
                }
            }
            let average = balance_rhs / balance_diag;
            state.flux[n][g] = average;
            for (d, part) in parts.iter().enumerate() {
                let faces = part.at(average);
                state.out_plus[d][n][g] = faces.out[0];
                state.out_minus[d][n][g] = faces.out[1];
                state.coefficients[d][n][g] = faces.coefficients;
            }
        }
        Ok(())
    }

    /// 方向 d の横方向の漏れ（他の2方向の漏れの和）を、両隣のノードを含む2次式で近似し、
    /// ノード n での f1, f2 の係数を返す。
    ///
    /// 外周では、反射境界なら外に鏡像のノードがあるとみなす。中性子束ゼロと真空の境界では、
    /// 漏れは中性子束に比例するとみなし、面の漏れをノードの漏れ × 面の中性子束 / ノード平均
    /// とする（中性子束ゼロなら 0）。真空の面のノード平均が正でなければエラー。
    fn transverse_shape(
        &self,
        state: &State,
        leakage: &[[f64; 3]],
        n: usize,
        d: usize,
        g: usize,
    ) -> Result<(f64, f64), NemError> {
        let transverse =
            |m: usize| -> f64 { (0..3).filter(|&e| e != d).map(|e| leakage[m][e]).sum() };
        let h = self.width(n, d);
        let side = |step: isize, face: usize| {
            Ok(match self.neighbor(n, d, step) {
                Some(m) => Beside::Node {
                    width: self.width(m, d),
                    leakage: transverse(m),
                },
                None => match self.problem.boundaries[2 * d + face] {
                    Boundary::Reflective => Beside::Node {
                        width: h,
                        leakage: transverse(n),
                    },
                    Boundary::ZeroFlux => Beside::Face(0.0),
                    Boundary::Vacuum => {
                        let flux = state.flux[n][g];
                        if !(flux > 0.0 && flux.is_finite()) {
                            return Err(NemError::VacuumFaceFlux { node: n, flux });
                        }
                        let (in_minus, in_plus) = self.incoming(state, n, d, g);
                        let face_flux = if step > 0 {
                            2.0 * (state.out_plus[d][n][g] + in_plus)
                        } else {
                            2.0 * (state.out_minus[d][n][g] + in_minus)
                        };
                        Beside::Face(transverse(n) * face_flux / flux)
                    }
                },
            })
        };
        Ok(transverse_moments(
            h,
            side(-1, 0)?,
            transverse(n),
            side(1, 1)?,
        ))
    }

    /// 外周の面から外へ出る正味の中性子数（群ごと）
    fn leakage(&self, state: &State) -> Vec<f64> {
        let mut total = vec![0.0; self.groups];
        for n in 0..self.nodes() {
            for d in 0..3 {
                let area = self.volume(n) / self.width(n, d);
                for (g, t) in total.iter_mut().enumerate() {
                    let (in_minus, in_plus) = self.incoming(state, n, d, g);
                    if self.neighbor(n, d, -1).is_none() {
                        *t += area * (state.out_minus[d][n][g] - in_minus);
                    }
                    if self.neighbor(n, d, 1).is_none() {
                        *t += area * (state.out_plus[d][n][g] - in_plus);
                    }
                }
            }
        }
        total
    }
}

/// 中性子束の相対変化。両方 0 なら変化なし。
fn relative_change(new: f64, old: f64) -> f64 {
    let scale = new.abs().max(old.abs());
    if scale == 0.0 {
        0.0
    } else {
        (new - old).abs() / scale
    }
}

/// 反復の変化量の列から、収束までに残る誤差を見積もる。
///
/// 変化が縮小率 ρ で等比的に減るとき、最新の変化を含めた残りの誤差は 変化/(1-ρ)。
/// 丸め誤差の水準（相対 1e-14）以下の変化は、増減が揺れても収束とみなす。
#[derive(Debug, Default)]
struct Remaining {
    previous: Option<f64>,
}

impl Remaining {
    fn update(&mut self, change: f64) -> f64 {
        const ROUNDING: f64 = 1.0e-14;
        let estimate = match self.previous {
            _ if change <= ROUNDING => 0.0,
            Some(previous) if previous > 0.0 && change < previous => {
                let ratio = change / previous;
                change * (1.0 + ratio / (1.0 - ratio))
            }
            _ => f64::INFINITY,
        };
        self.previous = Some(change);
        estimate
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// 2群の均質な物質（χ = (1, 0)、上方散乱なし）
    fn material() -> FewGroup {
        let (d1, d2) = (1.4, 0.4);
        let (a1, a2) = (0.010, 0.080);
        let s12 = 0.017;
        let (nf1, nf2) = (0.006, 0.11);
        let total = [1.0 / (3.0 * d1), 1.0 / (3.0 * d2)];
        FewGroup {
            bounds: vec![1.0e7, 0.625, 1.0e-5],
            // NEM は使わないが、少数群定数として正しい値にする
            k_inf: 1.0,
            total: total.to_vec(),
            absorption: vec![a1, a2],
            nu_fission: vec![nf1, nf2],
            fission: vec![nf1 / 2.4, nf2 / 2.4],
            chi: vec![1.0, 0.0],
            diffusion: vec![d1, d2],
            scatter: vec![vec![total[0] - a1 - s12, s12], vec![0.0, total[1] - a2]],
            adf: vec![vec![1.0, 1.0]; 4],
        }
    }

    /// 2群・χ = (1, 0)・上方散乱なしの均質な体系の k（バックリング b2）
    fn analytic_k(m: &FewGroup, b2: f64) -> f64 {
        let s12 = m.scatter[0][1];
        let r1 = m.absorption[0] + s12 + m.diffusion[0] * b2;
        let r2 = m.absorption[1] + m.diffusion[1] * b2;
        (m.nu_fission[0] + m.nu_fission[1] * s12 / r2) / r1
    }

    fn cuboid(size: [f64; 3], nodes: [usize; 3], boundaries: [Boundary; 6]) -> Problem {
        let split = |l: f64, n: usize| vec![l / n as f64; n];
        Problem {
            mesh: Mesh {
                dx: split(size[0], nodes[0]),
                dy: split(size[1], nodes[1]),
                dz: split(size[2], nodes[2]),
            },
            materials: vec![material()],
            node_material: vec![0; nodes.iter().product()],
            boundaries,
        }
    }

    #[test]
    fn reflective_box_gives_k_infinity() {
        let p = cuboid([30.0, 30.0, 40.0], [2, 3, 2], [Boundary::Reflective; 6]);
        let s = solve(&p, &Options::default()).unwrap();
        let expected = analytic_k(&material(), 0.0);
        assert!((s.k - expected).abs() < 1e-8, "{} vs {expected}", s.k);
    }

    #[test]
    fn bare_cuboid_matches_analytic_buckling() {
        // 誤差は主に横方向の漏れの2次近似によるもので、ノードを細かくすると速く減る
        let size = [60.0, 80.0, 100.0];
        let b2: f64 = size.iter().map(|l| (PI / l).powi(2)).sum();
        let expected = analytic_k(&material(), b2);
        let error = |n: usize| {
            let p = cuboid(size, [n; 3], [Boundary::ZeroFlux; 6]);
            (solve(&p, &Options::default()).unwrap().k - expected).abs() / expected
        };
        let (coarse, fine) = (error(4), error(8));
        assert!(fine < 2e-4, "relative error {fine}");
        assert!(fine < coarse / 8.0, "{coarse} -> {fine}");
    }

    #[test]
    fn symmetric_half_core_matches_full_core() {
        let full = cuboid([60.0, 60.0, 60.0], [6, 4, 4], [Boundary::ZeroFlux; 6]);
        let mut half_boundaries = [Boundary::ZeroFlux; 6];
        half_boundaries[0] = Boundary::Reflective;
        let half = cuboid([30.0, 60.0, 60.0], [3, 4, 4], half_boundaries);
        let (f, h) = (
            solve(&full, &Options::default()).unwrap(),
            solve(&half, &Options::default()).unwrap(),
        );
        assert!((f.k - h.k).abs() < 1e-7, "{} vs {}", f.k, h.k);
    }

    #[test]
    fn neutron_balance_holds() {
        let p = cuboid([60.0, 80.0, 100.0], [3, 4, 5], [Boundary::ZeroFlux; 6]);
        let s = solve(&p, &Options::default()).unwrap();
        let m = material();
        let volume = 60.0 * 80.0 * 100.0 / 60.0;
        let (mut absorption, mut production) = (0.0, 0.0);
        for phi in &s.flux {
            for ((f, a), nf) in phi.iter().zip(&m.absorption).zip(&m.nu_fission) {
                absorption += volume * a * f;
                production += volume * nf * f;
            }
        }
        let leakage: f64 = s.leakage.iter().sum();
        assert!(leakage > 0.0);
        assert!(
            (production / s.k - absorption - leakage).abs() < 1e-6 * production,
            "{production} / {} vs {absorption} + {leakage}",
            s.k
        );
    }

    #[test]
    fn assembly_sized_nodes_converge() {
        // 集合体程度のノード（16.7 cm）の裸の立方体。外周の境界条件の遅れで発散していた
        let size = [200.0; 3];
        let b2: f64 = size.iter().map(|l| (PI / l).powi(2)).sum();
        let expected = analytic_k(&material(), b2);
        let p = cuboid(size, [12; 3], [Boundary::ZeroFlux; 6]);
        let s = solve(&p, &Options::default()).unwrap();
        assert!(
            (s.k - expected).abs() < 1e-4 * expected,
            "{} vs {expected}",
            s.k
        );
    }

    #[test]
    fn single_reflective_node_direction_stays_bounded() {
        // z に 300 cm のノード1つで両側が反射。その方向の解は一様で、展開係数は 0 のまま
        // でなければならない。前の反復の a1, a2 から a3, a4 を求めると、丸め誤差が反復ごとに
        // 約 -5 倍と -2.3 倍になり、400 回ほどで NaN になっていた（#0025）
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[0] = Boundary::ZeroFlux;
        boundaries[1] = Boundary::ZeroFlux;
        let p = cuboid([60.0, 20.0, 300.0], [6, 1, 1], boundaries);
        let tight = Options {
            k_tolerance: 1e-14,
            flux_tolerance: 1e-13,
            max_iterations: 10_000,
        };
        let converged = solve(&p, &tight).unwrap();
        let core = Core::new(&p).unwrap();
        // 許容差 0 では収束しないので、上限の回数まで反復を続ける
        let long = Options {
            k_tolerance: 0.0,
            flux_tolerance: 0.0,
            max_iterations: 2000,
        };
        match iterate(&core, &long) {
            Err(NemError::NotConverged { k, .. }) => assert!(
                (k - converged.k).abs() < 1e-12 * converged.k,
                "{k} vs {}",
                converged.k
            ),
            Err(e) => panic!("expected NotConverged, got {e}"),
            Ok((s, _)) => panic!("converged with zero tolerances (k = {})", s.k),
        }
    }

    #[test]
    fn non_fundamental_mode_is_error() {
        // 断面 10 cm の裸の細い柱を1ノードで表すと、横方向の漏れが除去の約5倍になり、
        // 符号が交互に変わる解に収束する。基本モードではないので誤りとして返す
        let p = cuboid([10.0, 10.0, 100.0], [1, 1, 10], [Boundary::ZeroFlux; 6]);
        let err = solve(&p, &Options::default()).unwrap_err().to_string();
        assert!(err.contains("positive"), "{err}");
    }

    #[test]
    fn changes_at_rounding_level_count_as_converged() {
        let mut r = Remaining::default();
        r.update(1.0e-15);
        assert_eq!(r.update(2.0e-15), 0.0);
        let mut r = Remaining::default();
        assert_eq!(r.update(1.0), f64::INFINITY);
        // 縮小率 0.5 なら、最新の変化を含めて残りは 2 倍
        assert!((r.update(0.5) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn inconsistent_problem_is_error() {
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
        p.node_material.pop();
        assert!(solve(&p, &Options::default()).is_err());
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
        p.node_material[0] = 1;
        assert!(solve(&p, &Options::default()).is_err());
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
        p.mesh.dx[0] = 0.0;
        assert!(solve(&p, &Options::default()).is_err());
        // 有限でない断面積は、反復に入る前にエラーにする
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
        p.materials[0].nu_fission[1] = f64::NAN;
        let err = solve(&p, &Options::default()).unwrap_err().to_string();
        assert!(err.contains("finite"), "{err}");
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
        p.materials[0].diffusion[0] = f64::INFINITY;
        assert!(
            solve(&p, &Options::default())
                .unwrap_err()
                .to_string()
                .contains("finite")
        );
    }

    /// 1群の均質な物質（拡散係数 1 cm）
    fn one_group() -> FewGroup {
        let (d, absorption, nu_fission) = (1.0, 0.07, 0.08);
        let total = 1.0 / (3.0 * d);
        FewGroup {
            bounds: vec![1.0e7, 1.0e-5],
            // NEM は使わないが、少数群定数として正しい値にする
            k_inf: 1.0,
            total: vec![total],
            absorption: vec![absorption],
            nu_fission: vec![nu_fission],
            fission: vec![nu_fission / 2.4],
            chi: vec![1.0],
            diffusion: vec![d],
            scatter: vec![vec![total - absorption]],
            adf: vec![vec![1.0]; 4],
        }
    }

    /// 1群の均質な体系の k（バックリング b2）
    fn one_group_k(m: &FewGroup, b2: f64) -> f64 {
        m.nu_fission[0] / (m.absorption[0] + m.diffusion[0] * b2)
    }

    /// (0, π/x_max) にある f の最小の根（符号の変わる区間を探して二分法）
    fn smallest_root(f: impl Fn(f64) -> f64, upper: f64) -> f64 {
        let steps = 10_000;
        let mut a = upper * 1e-6;
        for i in 1..=steps {
            let b = upper * i as f64 / steps as f64;
            if f(a).signum() != f(b).signum() {
                let (mut lo, mut hi) = (a, b);
                for _ in 0..200 {
                    let mid = 0.5 * (lo + hi);
                    if f(lo).signum() == f(mid).signum() {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                return 0.5 * (lo + hi);
            }
            a = b;
        }
        panic!("no root below {upper}");
    }

    /// x 方向の平板（y, z は反射）
    fn slab(
        dx: Vec<f64>,
        materials: Vec<FewGroup>,
        node_material: Vec<usize>,
        x: Boundary,
    ) -> Problem {
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[0] = x;
        boundaries[1] = x;
        Problem {
            mesh: Mesh {
                dx,
                dy: vec![20.0],
                dz: vec![20.0],
            },
            materials,
            node_material,
            boundaries,
        }
    }

    /// x 方向に 6 ノードの平板。燃料と、断面積の違う燃料が交互に並ぶ
    fn two_material_slab(adf_a: Vec<Vec<f64>>, adf_b: Vec<Vec<f64>>) -> Problem {
        let mut a = material();
        a.adf = adf_a;
        let mut b = material();
        b.diffusion = vec![1.2, 0.5];
        b.absorption = vec![0.012, 0.10];
        b.nu_fission = vec![0.004, 0.09];
        b.adf = adf_b;
        slab(
            vec![10.0; 6],
            vec![a, b],
            vec![0, 1, 0, 1, 0, 1],
            Boundary::ZeroFlux,
        )
    }

    /// 収束した状態も返す
    fn solve_state(p: &Problem) -> (Core<'_>, Solution, State) {
        let core = Core::new(p).unwrap();
        let (s, state) = iterate(&core, &Options::default()).unwrap();
        (core, s, state)
    }

    /// 内側のすべての面で、不均質な中性子束 f φ（φ = 2(J+ + J-)）と正味の中性子流が連続であること。
    /// 辺の番号は（左、右、上、下）で、-x は左、+x は右、-y は下、+y は上、z の面の ADF は 1
    fn assert_interface_conditions(p: &Problem) {
        let (core, _, state) = solve_state(p);
        let adf = |n: usize, d: usize, plus: bool, g: usize| {
            let side = match (d, plus) {
                (0, false) => 0,
                (0, true) => 1,
                (1, true) => 2,
                (1, false) => 3,
                _ => return 1.0,
            };
            p.materials[p.node_material[n]].adf[side][g]
        };
        let mut checked = 0;
        for n in 0..core.nodes() {
            for d in 0..3 {
                let Some(m) = core.neighbor(n, d, 1) else {
                    continue;
                };
                for g in 0..core.groups {
                    // n の + 面と m の - 面
                    let (_, in_n) = core.incoming(&state, n, d, g);
                    let (in_m, _) = core.incoming(&state, m, d, g);
                    let (out_n, out_m) = (state.out_plus[d][n][g], state.out_minus[d][m][g]);
                    let (phi_n, phi_m) = (2.0 * (out_n + in_n), 2.0 * (out_m + in_m));
                    let (f_n, f_m) = (adf(n, d, true, g), adf(m, d, false, g));
                    assert!(
                        (phi_n / phi_m - f_m / f_n).abs() < 1e-6 * f_m / f_n,
                        "node {n} dir {d} group {g}: {} vs {}",
                        phi_n / phi_m,
                        f_m / f_n
                    );
                    let (j_n, j_m) = (out_n - in_n, in_m - out_m);
                    assert!(
                        (j_n - j_m).abs() < 1e-6 * phi_n,
                        "node {n} dir {d} group {g}: {j_n} vs {j_m}"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 0);
    }

    #[test]
    fn unit_adfs_leave_results_unchanged() {
        // ADF を入れる前（e3b627c）の結果
        let p = two_material_slab(vec![vec![1.0, 1.0]; 4], vec![vec![1.0, 1.0]; 4]);
        let s = solve(&p, &Options::default()).unwrap();
        let k = 8.510692818894939e-1;
        assert!((s.k - k).abs() < 1e-12, "{} vs {k}", s.k);
        let flux = [
            [0.0005911585345818675, 0.00012519241066775465],
            [0.0016614051334569624, 0.0003510047054023206],
        ];
        for (phi, expected) in s.flux.iter().zip(flux) {
            for (a, b) in phi.iter().zip(expected) {
                assert!((a - b).abs() < 1e-10 * b, "{a} vs {b}");
            }
        }
        assert_interface_conditions(&p);
    }

    #[test]
    fn adfs_make_homogeneous_flux_jump_at_faces() {
        // 群と辺で違う値にして、辺の対応も確かめる
        let a = vec![
            vec![1.05, 1.3],
            vec![1.1, 1.25],
            vec![0.9, 1.4],
            vec![0.95, 1.2],
        ];
        let b = vec![
            vec![0.97, 0.8],
            vec![0.93, 0.85],
            vec![1.02, 0.75],
            vec![0.99, 0.9],
        ];
        let p = two_material_slab(a.clone(), b.clone());
        assert_interface_conditions(&p);
        let unit = two_material_slab(vec![vec![1.0, 1.0]; 4], vec![vec![1.0, 1.0]; 4]);
        let (k, k_unit) = (
            solve(&p, &Options::default()).unwrap().k,
            solve(&unit, &Options::default()).unwrap().k,
        );
        assert!((k - k_unit).abs() > 1e-4, "{k} vs {k_unit}");

        // x と y に 3 × 2 の市松模様、z に 2 層
        let mut q = cuboid([30.0, 20.0, 40.0], [3, 2, 2], [Boundary::ZeroFlux; 6]);
        let mut m = material();
        m.adf = b;
        m.absorption[1] = 0.09;
        q.materials[0].adf = a;
        q.materials.push(m);
        q.node_material = (0..12).map(|n| (n + n / 3) % 2).collect();
        assert_interface_conditions(&q);
    }

    #[test]
    fn two_region_slab_matches_analytic_jump() {
        // 同じ断面積の領域 A (0, a) と B (a, a + b) の1群の平板、両端は中性子束 0。
        // φ_A = sin(Bx)、φ_B = c sin(B(a + b - x)) を f_A φ_A = f_B φ_B と中性子流の連続に入れると
        // f_A sin(Ba) cos(Bb) + f_B cos(Ba) sin(Bb) = 0
        let (a, b) = (30.0, 50.0);
        let (f_a, f_b) = (1.0, 1.5);
        let buckling = smallest_root(
            |x| f_a * (x * a).sin() * (x * b).cos() + f_b * (x * a).cos() * (x * b).sin(),
            PI / a,
        );
        let mut m_a = one_group();
        m_a.adf = vec![vec![f_a]; 4];
        let mut m_b = one_group();
        m_b.adf = vec![vec![f_b]; 4];
        let expected = one_group_k(&m_a, buckling * buckling);
        let p = slab(
            vec![5.0; 16],
            vec![m_a, m_b],
            (0..16).map(|i| usize::from(i >= 6)).collect(),
            Boundary::ZeroFlux,
        );
        let s = solve(&p, &Options::default()).unwrap();
        assert!(
            (s.k - expected).abs() < 1e-8 * expected,
            "{} vs {expected}",
            s.k
        );
        // ADF がなければ k は B = π/(a + b) のもの
        let plain = one_group_k(&one_group(), (PI / (a + b)).powi(2));
        assert!((expected - plain).abs() > 1e-3 * plain);
    }

    #[test]
    fn uniform_adf_scaling_keeps_k_without_vacuum() {
        // 外周が反射と中性子束 0 だけなら、面の条件は ADF の比だけで決まる（真空の面では変わる。
        // outer_adf_at_vacuum_matches_analytic）
        let a = vec![
            vec![1.05, 1.3],
            vec![1.1, 1.25],
            vec![0.9, 1.4],
            vec![0.95, 1.2],
        ];
        let b = vec![
            vec![0.97, 0.8],
            vec![0.93, 0.85],
            vec![1.02, 0.75],
            vec![0.99, 0.9],
        ];
        // 群ごとに違う倍率で、すべての物質とすべての辺の ADF を変える
        let scale = |adf: &[Vec<f64>]| -> Vec<Vec<f64>> {
            adf.iter().map(|s| vec![2.0 * s[0], 0.5 * s[1]]).collect()
        };
        let base = solve(
            &two_material_slab(a.clone(), b.clone()),
            &Options::default(),
        )
        .unwrap();
        let scaled = solve(
            &two_material_slab(scale(&a), scale(&b)),
            &Options::default(),
        )
        .unwrap();
        assert!(
            (base.k - scaled.k).abs() < 1e-9,
            "{} vs {}",
            base.k,
            scaled.k
        );
        for (x, y) in base.flux.iter().flatten().zip(scaled.flux.iter().flatten()) {
            assert!((x - y).abs() < 1e-6 * x, "{x} vs {y}");
        }
    }

    #[test]
    fn invalid_adfs_are_errors() {
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
        p.materials[0].adf.pop();
        assert!(solve(&p, &Options::default()).is_err());
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
        p.materials[0].adf[1].pop();
        assert!(solve(&p, &Options::default()).is_err());
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
            p.materials[0].adf[2][1] = bad;
            let err = solve(&p, &Options::default()).unwrap_err().to_string();
            assert!(err.contains("discontinuity"), "{err}");
        }
    }

    #[test]
    fn vacuum_face_has_no_incoming_heterogeneous_current() {
        // 不均質な中性子束 f φ = 2f(o + i) と正味の中性子流 o - i から、入射の部分中性子流は
        // f φ/4 - (o - i)/2 = (f(o + i) - (o - i))/2 で、真空の面ではこれが 0
        let mut m = material();
        m.adf[0] = vec![1.3, 0.8];
        m.adf[1] = vec![0.7, 1.6];
        let p = slab(vec![5.0; 6], vec![m.clone()], vec![0; 6], Boundary::Vacuum);
        let (core, s, state) = solve_state(&p);
        for g in 0..2 {
            let (in_minus, _) = core.incoming(&state, 0, 0, g);
            let (_, in_plus) = core.incoming(&state, 5, 0, g);
            let (out_minus, out_plus) = (state.out_minus[0][0][g], state.out_plus[0][5][g]);
            for (f, o, i) in [
                (m.adf[0][g], out_minus, in_minus),
                (m.adf[1][g], out_plus, in_plus),
            ] {
                assert!(o > 0.0);
                let incoming = (f * (o + i) - (o - i)) / 2.0;
                assert!(incoming.abs() < 1e-14 * o, "{incoming} for {o}");
            }
            // 外へ出る正味の中性子数は、外周の面の正味の中性子流の和
            let area = 20.0 * 20.0;
            let out = area * (out_minus - in_minus + out_plus - in_plus);
            assert!((s.leakage[g] - out).abs() < 1e-12 * out);
        }
    }

    #[test]
    fn vacuum_k_lies_between_zero_flux_and_reflective() {
        let k = |x: Boundary| {
            let p = slab(vec![5.0; 6], vec![material()], vec![0; 6], x);
            solve(&p, &Options::default()).unwrap().k
        };
        let (zero, vacuum, reflective) = (
            k(Boundary::ZeroFlux),
            k(Boundary::Vacuum),
            k(Boundary::Reflective),
        );
        assert!(
            zero < vacuum && vacuum < reflective,
            "{zero} {vacuum} {reflective}"
        );
    }

    /// 1群の裸の平板（幅 l、外周の面の ADF f）に、不均質な中性子束 f φ で Marshak の条件
    /// J- = f φ/4 + (D/2) dφ/dx = 0 を課した解 φ = cos(Bx) は tan(Bl/2) = f/(2DB) を満たす
    /// （外挿距離は 2D/f）
    fn marshak_buckling(m: &FewGroup, l: f64, f: f64) -> f64 {
        let d = m.diffusion[0];
        smallest_root(
            |b| (b * l / 2.0).sin() * 2.0 * d * b - f * (b * l / 2.0).cos(),
            PI / l,
        )
    }

    #[test]
    fn outer_adf_at_vacuum_matches_analytic() {
        // 1種類の物質なので内側の面には跳びがなく、外周の真空の面の ADF だけが効く。
        // ADF をすべて同じ倍率で変えると、k は解析解のとおり変わる
        let l = 40.0;
        // f が大きいほど外挿距離 2D/f が短く、漏れが増えて k は下がる
        let mut previous = f64::INFINITY;
        for f in [0.7, 1.0, 1.5, 2.0] {
            let mut m = one_group();
            m.adf = vec![vec![f]; 4];
            let b = marshak_buckling(&m, l, f);
            let expected = one_group_k(&m, b * b);
            let p = slab(vec![l / 8.0; 8], vec![m], vec![0; 8], Boundary::Vacuum);
            let s = solve(&p, &Options::default()).unwrap();
            assert!(
                (s.k - expected).abs() < 1e-8 * expected,
                "f = {f}: {} vs {expected}",
                s.k
            );
            assert!(expected < previous, "f = {f}");
            previous = expected;
        }
    }

    #[test]
    fn bare_slab_matches_marshak_analytic() {
        let m = one_group();
        let l = 40.0;
        let b = marshak_buckling(&m, l, 1.0);
        let expected = one_group_k(&m, b * b);
        let p = slab(vec![l / 8.0; 8], vec![m], vec![0; 8], Boundary::Vacuum);
        let s = solve(&p, &Options::default()).unwrap();
        assert!(
            (s.k - expected).abs() < 1e-8 * expected,
            "{} vs {expected}",
            s.k
        );
        // 外挿距離 2D を無視した中性子束 0 の解とは区別できる
        let zero = one_group_k(&one_group(), (PI / l).powi(2));
        assert!((expected - zero).abs() > 1e-3 * zero);
    }

    #[test]
    fn bare_cuboid_matches_marshak_analytic() {
        // 各方向の解は分離できるので B² は各方向の Marshak の解の和。横方向の漏れの近似を通る
        let m = one_group();
        let size = [40.0, 50.0, 60.0];
        let b2: f64 = size
            .iter()
            .map(|&l| marshak_buckling(&m, l, 1.0).powi(2))
            .sum();
        let expected = one_group_k(&m, b2);
        let error = |n: usize| {
            let mut p = cuboid(size, [n; 3], [Boundary::Vacuum; 6]);
            p.materials = vec![one_group()];
            (solve(&p, &Options::default()).unwrap().k - expected).abs() / expected
        };
        let (coarse, fine) = (error(4), error(8));
        assert!(fine < 1e-4, "relative error {fine}");
        assert!(fine < coarse / 8.0, "{coarse} -> {fine}");
    }

    #[test]
    fn side_mapping_matches_two_region_analytic() {
        // 2領域の平板を x、y、z の各方向に置き、ADF の跳びを解析解（two_region_slab_matches_analytic_jump）
        // と比べる。跳びを境の2つのノードだけに置く: 境の A のノードは - 面で 1、+ 面で f_A、境の B の
        // ノードは - 面で f_B、+ 面で 1、ほかのノードと辺は 1 か関係のない値（その方向に隣がなく、
        // 反射の面では ADF は効かない）。- 面と + 面の辺を取り違えれば、跳びが境からずれて解析解から外れる
        let (a, b) = (30.0, 50.0);
        let (f_a, f_b) = (1.0, 1.5);
        let analytic = |f_a: f64, f_b: f64| {
            let buckling = smallest_root(
                |x| f_a * (x * a).sin() * (x * b).cos() + f_b * (x * a).cos() * (x * b).sin(),
                PI / a,
            );
            one_group_k(&one_group(), buckling * buckling)
        };
        // 方向ごとの（- 面, + 面）の辺（左、右、上、下の番号）。z の面の ADF は 1 なので跳びはない
        let cases = [
            (Some((0, 1)), analytic(f_a, f_b)),
            (Some((3, 2)), analytic(f_a, f_b)),
            (None, analytic(1.0, 1.0)),
        ];
        for (d, (sides, expected)) in cases.into_iter().enumerate() {
            let bulk = one_group();
            let mut edge_a = one_group();
            let mut edge_b = one_group();
            match sides {
                Some((minus, plus)) => {
                    edge_a.adf = vec![vec![3.0]; 4];
                    edge_a.adf[minus] = vec![1.0];
                    edge_a.adf[plus] = vec![f_a];
                    edge_b.adf = vec![vec![0.5]; 4];
                    edge_b.adf[minus] = vec![f_b];
                    edge_b.adf[plus] = vec![1.0];
                }
                None => {
                    edge_a.adf = vec![vec![f_a]; 4];
                    edge_b.adf = vec![vec![f_b]; 4];
                }
            }
            let mut widths = [vec![20.0], vec![20.0], vec![20.0]];
            widths[d] = vec![5.0; 16];
            let mut boundaries = [Boundary::Reflective; 6];
            boundaries[2 * d] = Boundary::ZeroFlux;
            boundaries[2 * d + 1] = Boundary::ZeroFlux;
            let [dx, dy, dz] = widths;
            let p = Problem {
                mesh: Mesh { dx, dy, dz },
                materials: vec![bulk, edge_a, edge_b],
                node_material: (0..16)
                    .map(|i| match i {
                        5 => 1,
                        6 => 2,
                        _ => 0,
                    })
                    .collect(),
                boundaries,
            };
            let k = solve(&p, &Options::default()).unwrap().k;
            assert!(
                (k - expected).abs() < 1e-8 * expected,
                "direction {d}: {k} vs {expected}"
            );
        }
    }

    #[test]
    fn vacuum_corner_converges_with_refinement() {
        // 40 cm 角の 2D（-x, -y は反射、+x, +y は真空）。20 cm 角の4区画に2種類の物質を市松に置き、
        // 真空の角にも物質の境がある分離できない問題で、細かいノードの解に近づくこと。
        // 横方向の漏れの真空の面での扱いと、物質の境の ADF を通る
        let mut a = material();
        a.adf = vec![vec![1.1, 1.3]; 4];
        let mut b = material();
        b.diffusion = vec![1.2, 0.5];
        b.absorption = vec![0.012, 0.10];
        b.nu_fission = vec![0.004, 0.09];
        b.adf = vec![vec![0.95, 0.8]; 4];
        let k = |n: usize| {
            let mut boundaries = [Boundary::Reflective; 6];
            boundaries[1] = Boundary::Vacuum;
            boundaries[3] = Boundary::Vacuum;
            let p = Problem {
                mesh: Mesh {
                    dx: vec![40.0 / n as f64; n],
                    dy: vec![40.0 / n as f64; n],
                    dz: vec![20.0],
                },
                materials: vec![a.clone(), b.clone()],
                node_material: (0..n * n)
                    .map(|i| usize::from((i % n < n / 2) != (i / n < n / 2)))
                    .collect(),
                boundaries,
            };
            solve(&p, &Options::default()).unwrap().k
        };
        let reference = k(32);
        let error = |n: usize| (k(n) - reference).abs() / reference;
        let (coarse, fine) = (error(4), error(8));
        // 測った値は 9.6e-5 -> 1.1e-5
        assert!(coarse < 2e-4, "relative error {coarse}");
        assert!(fine < coarse / 6.0, "{coarse} -> {fine}");
    }

    #[test]
    fn no_flux_next_to_a_vacuum_face_is_error() {
        // 真空の面の横方向の漏れは面の中性子束 / ノード平均を使うので、ノード平均が正でなければ止める
        let p = slab(vec![5.0; 6], vec![material()], vec![0; 6], Boundary::Vacuum);
        let core = Core::new(&p).unwrap();
        for bad in [0.0, -1.0, f64::NAN] {
            let mut state = State::new(6, 2);
            state.flux[0][1] = bad;
            let err = core
                .sweep_group(&mut state, 1, 1.0)
                .unwrap_err()
                .to_string();
            assert!(err.contains("vacuum"), "{err}");
        }
        // 真空の面がなければ問題にしない
        let p = slab(
            vec![5.0; 6],
            vec![material()],
            vec![0; 6],
            Boundary::ZeroFlux,
        );
        let core = Core::new(&p).unwrap();
        let mut state = State::new(6, 2);
        state.flux[0][1] = 0.0;
        assert!(core.sweep_group(&mut state, 1, 1.0).is_ok());
    }

    /// 細かいノードの解を参照（不均質な解）とし、粗いノード（x, y に `blocks` 個ずつの細かいノードの
    /// まとまり）ごとに反応率を保存して定数を均質化し、GET の1ノードの問題（ADR-0094）で ADF を作った
    /// 粗いノードの問題と、細かいノードの解を返す。z は1層。
    ///
    /// 1ノードの問題の外周の扱い（`Edge`）は、粗いノードの面が外周なら細かいノードの問題の境界とする。
    /// 内側の面（隣の粗いノード）の扱いは1ノードの問題にないので、内側の面は、その方向の横方向の漏れが
    /// 0 の問題（y と z に1つの細かいノードで、どちらも反射）でだけ使う
    fn coarsen_with_get(fine: &Problem, blocks: [usize; 2]) -> (Problem, Solution) {
        use xlc_numeric::nem::{Edge, OneNode, surface_fluxes};
        let (core, solution, state) = solve_state(fine);
        let core = &core;
        let [nx, ny, nz] = core.counts;
        assert_eq!(nz, 1);
        let [bx, by] = blocks;
        let (cx, cy) = (nx / bx, ny / by);
        assert_eq!((cx * bx, cy * by), (nx, ny));
        assert!(
            cx * cy == 1 || ny == 1,
            "inner faces need zero transverse leakage"
        );
        let ng = core.groups;
        let mut materials = Vec::new();
        for (cj, ci) in (0..cy).flat_map(|j| (0..cx).map(move |i| (j, i))) {
            let members: Vec<usize> = (cj * by..(cj + 1) * by)
                .flat_map(|j| (ci * bx..(ci + 1) * bx).map(move |i| core.index([i, j, 0])))
                .collect();
            // 反応率を保存する均質化（D は中性子束の重みの平均。GET の ADF が差を吸収する）
            let mut volume = 0.0;
            let mut flux = vec![0.0; ng];
            let mut m = FewGroup {
                bounds: fine.materials[0].bounds.clone(),
                k_inf: solution.k,
                total: vec![0.0; ng],
                absorption: vec![0.0; ng],
                nu_fission: vec![0.0; ng],
                fission: vec![0.0; ng],
                chi: vec![0.0; ng],
                diffusion: vec![0.0; ng],
                scatter: vec![vec![0.0; ng]; ng],
                adf: vec![vec![1.0; ng]; 4],
            };
            for &n in &members {
                let (v, f, x) = (core.volume(n), &solution.flux[n], core.material(n));
                volume += v;
                let production: f64 = (0..ng).map(|g| x.nu_fission[g] * f[g]).sum();
                for g in 0..ng {
                    let w = v * f[g];
                    flux[g] += w;
                    m.total[g] += w * x.total[g];
                    m.nu_fission[g] += w * x.nu_fission[g];
                    m.fission[g] += w * x.fission[g];
                    m.diffusion[g] += w * x.diffusion[g];
                    m.chi[g] += v * x.chi[g] * production;
                    for to in 0..ng {
                        m.scatter[g][to] += w * x.scatter[g][to];
                    }
                }
            }
            for (g, f) in flux.iter().enumerate() {
                for value in [
                    &mut m.total[g],
                    &mut m.nu_fission[g],
                    &mut m.fission[g],
                    &mut m.diffusion[g],
                ] {
                    *value /= f;
                }
                for value in &mut m.scatter[g] {
                    *value /= f;
                }
                m.absorption[g] = m.total[g] - m.scatter[g].iter().sum::<f64>();
            }
            let produced: f64 = m.chi.iter().sum();
            for c in &mut m.chi {
                *c /= produced;
            }
            let average: Vec<f64> = flux.iter().map(|f| f / volume).collect();

            // 面（-x, +x, -y, +y）の不均質な中性子束と外向きの正味の中性子流の、面にわたる平均
            let mut het = [vec![0.0; ng], vec![0.0; ng], vec![0.0; ng], vec![0.0; ng]];
            let mut current = het.clone();
            for (face, d, plus) in [(0, 0, false), (1, 0, true), (2, 1, false), (3, 1, true)] {
                // 1つのノードで両側が反射の方向では、解はその方向に一様で、面の中性子束はノード平均。
                // その方向の展開係数は使わない（両側の正味の中性子流が 0 に決まり、前の反復の a1, a2 から
                // a3, a4 を求める反復が丸め誤差を増やす。値は解に効かない）
                if core.counts[d] == 1
                    && fine.boundaries[2 * d] == Boundary::Reflective
                    && fine.boundaries[2 * d + 1] == Boundary::Reflective
                {
                    het[face].clone_from(&average);
                    continue;
                }
                let on_face: Vec<usize> = members
                    .iter()
                    .copied()
                    .filter(|&n| {
                        let p = core.position(n);
                        let lowest = if d == 0 { ci * bx } else { cj * by };
                        let count = if d == 0 { bx } else { by };
                        p[d] == if plus { lowest + count - 1 } else { lowest }
                    })
                    .collect();
                let length: f64 = on_face
                    .iter()
                    .map(|&n| core.volume(n) / core.width(n, d))
                    .sum();
                for &n in &on_face {
                    let w = core.volume(n) / core.width(n, d) / length;
                    for g in 0..ng {
                        let (in_minus, in_plus) = core.incoming(&state, n, d, g);
                        let (o, i) = if plus {
                            (state.out_plus[d][n][g], in_plus)
                        } else {
                            (state.out_minus[d][n][g], in_minus)
                        };
                        het[face][g] += w * 2.0 * (o + i);
                        current[face][g] += w * (o - i);
                    }
                }
            }
            let edge = |d: usize, plus: bool| {
                let outer = if d == 0 {
                    if plus { ci == cx - 1 } else { ci == 0 }
                } else if plus {
                    cj == cy - 1
                } else {
                    cj == 0
                };
                match (outer, fine.boundaries[2 * d + usize::from(plus)]) {
                    (true, Boundary::Vacuum) => Edge::Vacuum,
                    (true, Boundary::Reflective) | (false, _) => Edge::Reflective,
                    (true, Boundary::ZeroFlux) => panic!("no zero-flux edge in GET"),
                }
            };
            let widths = [
                core.widths[0][ci * bx..(ci + 1) * bx].iter().sum(),
                core.widths[1][cj * by..(cj + 1) * by].iter().sum(),
            ];
            let surface = surface_fluxes(&OneNode {
                diffusion: &m.diffusion,
                total: &m.total,
                scatter: &m.scatter,
                nu_fission: &m.nu_fission,
                chi: &m.chi,
                k: solution.k,
                flux: &average,
                widths,
                currents: [&current[0], &current[1], &current[2], &current[3]],
                edges: [edge(0, false), edge(0, true), edge(1, false), edge(1, true)],
            })
            .unwrap();
            // 少数群定数の辺は（左、右、上、下）で、-x が左、+x が右、-y が下、+y が上
            for (side, face) in [(0, 0), (1, 1), (2, 3), (3, 2)] {
                m.adf[side] = (0..ng).map(|g| het[face][g] / surface[face][g]).collect();
            }
            materials.push(m);
        }
        let merge = |widths: &[f64], b: usize| -> Vec<f64> {
            widths.chunks(b).map(|c| c.iter().sum()).collect()
        };
        let coarse = Problem {
            mesh: Mesh {
                dx: merge(&fine.mesh.dx, bx),
                dy: merge(&fine.mesh.dy, by),
                dz: fine.mesh.dz.clone(),
            },
            node_material: (0..materials.len()).collect(),
            materials,
            boundaries: fine.boundaries,
        };
        (coarse, solution)
    }

    /// 粗いノードの問題の k とノード平均の中性子束が、細かいノードの解（同じ規格化: 全体の生成率 1）を
    /// 再現することを確かめ、ADF をすべて 1 にしたときの k の相対差を返す
    fn assert_get_equivalence(coarse: &Problem, fine: &Solution, blocks: [usize; 2], nx: usize) {
        let s = solve(coarse, &Options::default()).unwrap();
        assert!(
            (s.k - fine.k).abs() < 1e-7 * fine.k,
            "{} vs {}",
            s.k,
            fine.k
        );
        let [bx, by] = blocks;
        let cx = nx / bx;
        for (c, phi) in s.flux.iter().enumerate() {
            let (ci, cj) = (c % cx, c / cx);
            let members: Vec<usize> = (cj * by..(cj + 1) * by)
                .flat_map(|j| (ci * bx..(ci + 1) * bx).map(move |i| i + nx * j))
                .collect();
            for (g, value) in phi.iter().enumerate() {
                // 細かいノードはどれも同じ大きさ
                let expected =
                    members.iter().map(|&n| fine.flux[n][g]).sum::<f64>() / members.len() as f64;
                assert!(
                    (value - expected).abs() < 1e-6 * expected,
                    "node {c} group {g}: {value} vs {expected}"
                );
            }
        }
        let mut unit = coarse.clone();
        for m in &mut unit.materials {
            m.adf = vec![vec![1.0; m.groups()]; 4];
        }
        let k_unit = solve(&unit, &Options::default()).unwrap().k;
        // ADF がなければ、粗いノードの NEM の誤差の分だけ k が外れる
        assert!(
            (k_unit - fine.k).abs() > 1e-4 * fine.k,
            "{k_unit} vs {}",
            fine.k
        );
    }

    /// 細かいノードの物質: 燃料、断面積の違う燃料、核分裂の少ない燃料
    fn three_fuels() -> Vec<FewGroup> {
        let a = material();
        let mut b = material();
        b.diffusion = vec![1.2, 0.5];
        b.absorption = vec![0.012, 0.10];
        b.nu_fission = vec![0.004, 0.09];
        let mut c = material();
        c.absorption = vec![0.009, 0.06];
        c.nu_fission = vec![0.002, 0.04];
        c.fission = vec![0.002 / 2.4, 0.04 / 2.4];
        vec![a, b, c]
    }

    #[test]
    fn get_adfs_reproduce_a_two_node_slab() {
        // x 方向に 20 cm の粗いノード2つ（細かいノード 2.5 cm × 8 ずつ、中に2種類の物質）。-x は反射、
        // +x は真空。粗いノードの間の面で、ADF をかけた不均質な中性子束と正味の中性子流がつながる
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[1] = Boundary::Vacuum;
        let fine = Problem {
            mesh: Mesh {
                dx: vec![2.5; 16],
                dy: vec![20.0],
                dz: vec![20.0],
            },
            materials: three_fuels(),
            node_material: vec![0, 1, 0, 1, 1, 0, 1, 0, 2, 2, 1, 1, 2, 2, 0, 0],
            boundaries,
        };
        let (coarse, solution) = coarsen_with_get(&fine, [8, 1]);
        assert_get_equivalence(&coarse, &solution, [8, 1], 16);
    }

    #[test]
    fn get_adfs_reproduce_a_node_with_vacuum_sides() {
        // 20 cm 角の粗いノード1つ（細かいノード 5 cm 角 × 16、3種類の物質）。-x と -y は反射、+x と +y は
        // 真空。真空の面での横方向の漏れの扱いと、2方向の横方向の漏れを通る
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[1] = Boundary::Vacuum;
        boundaries[3] = Boundary::Vacuum;
        let fine = Problem {
            mesh: Mesh {
                dx: vec![5.0; 4],
                dy: vec![5.0; 4],
                dz: vec![20.0],
            },
            materials: three_fuels(),
            node_material: vec![0, 1, 0, 2, 1, 0, 2, 0, 0, 2, 1, 1, 2, 0, 1, 0],
            boundaries,
        };
        let (coarse, solution) = coarsen_with_get(&fine, [4, 4]);
        assert_get_equivalence(&coarse, &solution, [4, 4], 4);
    }
}
