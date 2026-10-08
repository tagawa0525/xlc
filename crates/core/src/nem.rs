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
            let sized = |v: &[f64]| v.len() == groups;
            let consistent = m.groups() == groups
                && sized(&m.diffusion)
                && sized(&m.nu_fission)
                && sized(&m.chi)
                && sized(&m.absorption)
                && m.scatter.len() == groups
                && m.scatter.iter().all(|row| sized(row));
            if !consistent {
                return Err(format!("material {i} must have {groups} groups"));
            }
            let all_finite = [&m.total, &m.diffusion, &m.nu_fission, &m.chi, &m.absorption]
                .into_iter()
                .chain(&m.scatter)
                .flatten()
                .all(|v| v.is_finite());
            if !all_finite {
                return Err(format!("material {i}: cross sections must be finite"));
            }
            let adf_valid = m.adf.len() == 4
                && m.adf.iter().all(|side| {
                    side.len() == groups && side.iter().all(|&f| f > 0.0 && f.is_finite())
                });
            if !adf_valid {
                return Err(format!(
                    "material {i}: discontinuity factors must be 4 sides of {groups} positive values"
                ));
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
            let mut parts = [Direction::default(); 3];
            let mut balance_rhs = source;
            let mut balance_diag = removal;
            for (d, part) in parts.iter_mut().enumerate() {
                let h = self.width(n, d);
                let alpha = m.diffusion[g] / h;
                // 源のモーメント（他の群・核分裂の展開係数から）
                let mut s = [0.0; 4];
                for (from, c) in coupling.iter().enumerate() {
                    for (si, a) in s.iter_mut().zip(&state.coefficients[d][n][from]) {
                        *si += c * a;
                    }
                }
                // ノード n の群 g の値はまだ更新していないので、スイープの始めの値と同じ
                let (l1, l2) = self.transverse_shape(state, &leakage, n, d, g)?;
                s[0] -= l1;
                s[1] -= l2;

                let [a1, a2, _, _] = state.coefficients[d][n][g];
                let dh2 = m.diffusion[g] / (h * h);
                // f1, f2 で重み付けした残差から 3次・4次の係数
                let a3 = (removal * a1 / 12.0 - s[0] / 12.0 + s[2] / 120.0)
                    / (dh2 / 2.0 + removal / 120.0);
                let a4 = (removal * a2 / 20.0 - s[1] / 20.0 + s[3] / 700.0)
                    / (dh2 / 5.0 + removal / 700.0);

                // 面の中性子束と正味の中性子流の連続から、出ていく部分中性子流 o = (o_R, o_L) は
                //   M o = A i + c,  M = [[1+8α, 4α], [4α, 1+8α]],  A = [[1-8α, -4α], [-4α, 1-8α]]
                // を満たす（c はノード平均と a3, a4 の項）。入射は i = β o + 外からの入射
                // （face_inputs）なので、それも同時に解く
                let ((beta_l, ext_l), (beta_r, ext_r)) = self.face_inputs(state, n, d, g);
                let (p, q) = (1.0 + 8.0 * alpha, 1.0 - 8.0 * alpha);
                let w = 4.0 * alpha;
                // K = M - A B（B = diag(β_R, β_L)）
                let k = [
                    [p - q * beta_r, w + w * beta_l],
                    [w + w * beta_r, p - q * beta_l],
                ];
                let det = k[0][0] * k[1][1] - k[0][1] * k[1][0];
                let solve2 = |r: [f64; 2]| {
                    [
                        (k[1][1] * r[0] - k[0][1] * r[1]) / det,
                        (k[0][0] * r[1] - k[1][0] * r[0]) / det,
                    ]
                };
                // φ̄ = 0 のときの値と、φ̄ に対する傾き
                let base = solve2([
                    q * ext_r - w * ext_l - alpha * a3 / 2.0 - alpha * a4 / 5.0,
                    q * ext_l - w * ext_r + alpha * a3 / 2.0 - alpha * a4 / 5.0,
                ]);
                let slope = solve2([6.0 * alpha, 6.0 * alpha]);
                *part = Direction {
                    out: base,
                    slope,
                    beta: [beta_r, beta_l],
                    external: [ext_r, ext_l],
                    a3,
                    a4,
                };
                // 正味の漏れ (o_R + o_L - i_R - i_L)/h（i = β o + 外からの入射）をバランス式に入れる
                for f in 0..2 {
                    balance_rhs -= ((1.0 - part.beta[f]) * part.out[f] - part.external[f]) / h;
                    balance_diag += (1.0 - part.beta[f]) * part.slope[f] / h;
                }
            }
            let average = balance_rhs / balance_diag;
            state.flux[n][g] = average;
            for (d, part) in parts.iter().enumerate() {
                let out: [f64; 2] = std::array::from_fn(|f| part.out[f] + part.slope[f] * average);
                let incoming: [f64; 2] =
                    std::array::from_fn(|f| part.beta[f] * out[f] + part.external[f]);
                state.out_plus[d][n][g] = out[0];
                state.out_minus[d][n][g] = out[1];
                let (phi_r, phi_l) = (2.0 * (out[0] + incoming[0]), 2.0 * (out[1] + incoming[1]));
                state.coefficients[d][n][g] = [
                    phi_r - phi_l,
                    phi_r + phi_l - 2.0 * average,
                    part.a3,
                    part.a4,
                ];
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
        // ノードの中心を原点とした q(x) = c0 + c1 x + c2 x² への条件 (行, 右辺)
        let average = |a: f64, b: f64, value: f64| {
            ([1.0, 0.5 * (a + b), (a * a + a * b + b * b) / 3.0], value)
        };
        let side = |step: isize, face: usize| {
            let sign = step as f64;
            let near = sign * 0.5 * h;
            Ok(match self.neighbor(n, d, step) {
                Some(m) => {
                    let far = near + sign * self.width(m, d);
                    average(near.min(far), near.max(far), transverse(m))
                }
                None => match self.problem.boundaries[2 * d + face] {
                    Boundary::Reflective => {
                        let far = near + sign * h;
                        average(near.min(far), near.max(far), transverse(n))
                    }
                    Boundary::ZeroFlux => ([1.0, near, near * near], 0.0),
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
                        ([1.0, near, near * near], transverse(n) * face_flux / flux)
                    }
                },
            })
        };
        let rows = [
            side(-1, 0)?,
            average(-0.5 * h, 0.5 * h, transverse(n)),
            side(1, 1)?,
        ];
        let c = solve3(
            [rows[0].0, rows[1].0, rows[2].0],
            [rows[0].1, rows[1].1, rows[2].1],
        );
        Ok((c[1] * h, c[2] * h * h / 3.0))
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

/// 3×3 の連立一次方程式（Cramer の公式）。
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(a);
    std::array::from_fn(|j| {
        let mut m = a;
        for (row, value) in m.iter_mut().zip(b) {
            row[j] = value;
        }
        det(m) / d
    })
}

/// ノードの1方向の途中結果。`out` は φ̄ = 0 のときの出ていく部分中性子流で、
/// 実際の値は `slope` × φ̄ を足したもの。
#[derive(Debug, Clone, Copy, Default)]
struct Direction {
    /// [+ 面, - 面]
    out: [f64; 2],
    slope: [f64; 2],
    /// 外周の面での入射 i = β o の β（隣のノードがある面では 0）
    beta: [f64; 2],
    /// 隣のノードから入ってくる部分中性子流（外周の面では 0）
    external: [f64; 2],
    a3: f64,
    a4: f64,
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
            k_inf: 0.0,
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
            k_inf: 0.0,
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
        // ADF を入れる前（47c63a7）の結果
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
}
