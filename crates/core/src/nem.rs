//! 3次元直交座標の多群ノード展開法（NEM）。
//!
//! 各ノード・各方向で横方向に積分した1次元の拡散方程式を、ノード平均と4次までの多項式
//! f1 = ξ, f2 = 3ξ² - 1/4, f3 = ξ³ - ξ/4, f4 = (ξ² - 1/20)(ξ² - 1/4)（ξ ∈ [-1/2, 1/2]）で解く。
//! 3次・4次の係数は f1, f2 で重み付けした残差から、横方向の漏れは隣のノードを含む
//! 2次式で近似する。ノード間は部分中性子流でつなぐ。

use lattice::fewgroup::FewGroup;

/// 外周の境界条件。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// 正味の中性子流 0（対称面）
    Reflective,
    /// 境界での中性子束 0
    ZeroFlux,
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

pub fn solve(problem: &Problem, options: &Options) -> Result<Solution, String> {
    let core = Core::new(problem)?;
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
        return Err("no fissile material".into());
    }
    state.scale(1.0 / initial);
    let mut k_history = Remaining::default();
    let mut flux_history = Remaining::default();

    for iteration in 1..=options.max_iterations {
        let previous = state.flux.clone();
        for g in 0..ng {
            core.sweep_group(&mut state, g, k);
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
        if std::env::var("NEM_TRACE").is_ok() && iteration < 40 {
            eprintln!(
                "it {iteration} k {new_k} flux0 {:?} c {:?}",
                state.flux[0], state.coefficients[0][0]
            );
        }
        let k_remaining = k_history.update((new_k - k).abs());
        let flux_remaining = flux_history.update(change);
        k = new_k;
        if k_remaining < options.k_tolerance && flux_remaining < options.flux_tolerance {
            return Ok(Solution {
                k,
                leakage: core.leakage(&state),
                flux: state.flux,
                iterations: iteration,
            });
        }
    }
    Err(format!(
        "not converged in {} iterations (k = {k})",
        options.max_iterations
    ))
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

    /// ノード n の - 面・+ 面に入ってくる部分中性子流。外周では境界条件から求める。
    fn incoming(&self, state: &State, n: usize, d: usize, g: usize) -> (f64, f64) {
        let from = |step: isize, own_out: f64, face: usize| match self.neighbor(n, d, step) {
            Some(m) if step < 0 => state.out_plus[d][m][g],
            Some(m) => state.out_minus[d][m][g],
            None => match self.problem.boundaries[2 * d + face] {
                Boundary::Reflective => own_out,
                // 面の中性子束 2(J+ + J-) = 0
                Boundary::ZeroFlux => -own_out,
            },
        };
        (
            from(-1, state.out_minus[d][n][g], 0),
            from(1, state.out_plus[d][n][g], 1),
        )
    }

    /// ノード n の方向 d の正味の漏れ (J+ - J-)/h（単位体積あたり）
    fn net_leakage(&self, state: &State, n: usize, d: usize, g: usize) -> f64 {
        let (in_minus, in_plus) = self.incoming(state, n, d, g);
        let j_plus = state.out_plus[d][n][g] - in_plus;
        let j_minus = in_minus - state.out_minus[d][n][g];
        (j_plus - j_minus) / self.width(n, d)
    }

    fn sweep_group(&self, state: &mut State, g: usize, k: f64) {
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
                let (l1, l2) = self.transverse_shape(&leakage, n, d);
                s[0] -= l1;
                s[1] -= l2;

                let [a1, a2, _, _] = state.coefficients[d][n][g];
                let dh2 = m.diffusion[g] / (h * h);
                // f1, f2 で重み付けした残差から 3次・4次の係数
                let a3 = (removal * a1 / 12.0 - s[0] / 12.0 + s[2] / 120.0)
                    / (dh2 / 2.0 + removal / 120.0);
                let a4 = (removal * a2 / 20.0 - s[1] / 20.0 + s[3] / 700.0)
                    / (dh2 / 5.0 + removal / 700.0);

                // 面の中性子束と正味の中性子流の連続から求めた、出ていく部分中性子流の
                // φ̄ = 0 のときの値。φ̄ の係数は両面とも 6α/(1 + 12α)
                let (in_l, in_r) = self.incoming(state, n, d, g);
                let b_r = in_r * (1.0 - 8.0 * alpha)
                    - 4.0 * alpha * in_l
                    - alpha * a3 / 2.0
                    - alpha * a4 / 5.0;
                let b_l = in_l * (1.0 - 8.0 * alpha) - 4.0 * alpha * in_r + alpha * a3 / 2.0
                    - alpha * a4 / 5.0;
                let diag = 1.0 + 8.0 * alpha;
                let det = diag * diag - 16.0 * alpha * alpha;
                *part = Direction {
                    out_r: (diag * b_r - 4.0 * alpha * b_l) / det,
                    out_l: (diag * b_l - 4.0 * alpha * b_r) / det,
                    slope: 6.0 * alpha / (1.0 + 12.0 * alpha),
                    in_r,
                    in_l,
                    a3,
                    a4,
                };
                // 正味の漏れ (o_R + o_L - i_R - i_L)/h をバランス式に入れる
                balance_rhs -= (part.out_r + part.out_l - in_r - in_l) / h;
                balance_diag += 2.0 * part.slope / h;
            }
            let average = balance_rhs / balance_diag;
            state.flux[n][g] = average;
            for (d, part) in parts.iter().enumerate() {
                let out_r = part.out_r + part.slope * average;
                let out_l = part.out_l + part.slope * average;
                state.out_plus[d][n][g] = out_r;
                state.out_minus[d][n][g] = out_l;
                let (phi_r, phi_l) = (2.0 * (out_r + part.in_r), 2.0 * (out_l + part.in_l));
                state.coefficients[d][n][g] = [
                    phi_r - phi_l,
                    phi_r + phi_l - 2.0 * average,
                    part.a3,
                    part.a4,
                ];
            }
        }
    }

    /// 方向 d の横方向の漏れ（他の2方向の漏れの和）を、両隣のノードを含む2次式で近似し、
    /// ノード n での f1, f2 の係数を返す。
    ///
    /// 外周では、反射境界なら外に鏡像のノードがあるとみなし、中性子束ゼロの境界なら
    /// 面で漏れが 0 になる（漏れは中性子束に比例する）という条件を使う。
    fn transverse_shape(&self, leakage: &[[f64; 3]], n: usize, d: usize) -> (f64, f64) {
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
            match self.neighbor(n, d, step) {
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
                },
            }
        };
        let rows = [
            side(-1, 0),
            average(-0.5 * h, 0.5 * h, transverse(n)),
            side(1, 1),
        ];
        let c = solve3(
            [rows[0].0, rows[1].0, rows[2].0],
            [rows[0].1, rows[1].1, rows[2].1],
        );
        (c[1] * h, c[2] * h * h / 3.0)
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

/// ノードの1方向の途中結果。`out_r`, `out_l` は φ̄ = 0 のときの出ていく部分中性子流で、
/// 実際の値は `slope` × φ̄ を足したもの。
#[derive(Debug, Clone, Copy, Default)]
struct Direction {
    out_r: f64,
    out_l: f64,
    slope: f64,
    in_r: f64,
    in_l: f64,
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

/// 反復の変化量の列から、収束までに残る誤差を見積もる（変化が縮小率 ρ で等比的に減るとき、
/// 残りの誤差は 最新の変化 × ρ/(1-ρ)）。
#[derive(Debug, Default)]
struct Remaining {
    previous: Option<f64>,
}

impl Remaining {
    fn update(&mut self, change: f64) -> f64 {
        let estimate = match self.previous {
            _ if change == 0.0 => 0.0,
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
    }
}
