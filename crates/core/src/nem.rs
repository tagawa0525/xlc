//! 3次元直交座標の多群ノード展開法（NEM）。
//!
//! 各ノード・各方向で横方向に積分した1次元の拡散方程式を、ノード平均と4次までの多項式
//! f1 = ξ, f2 = 3ξ² - 1/4, f3 = ξ³ - ξ/4, f4 = (ξ² - 1/20)(ξ² - 1/4)（ξ ∈ [-1/2, 1/2]）で解く。
//! 3次・4次の係数は f1, f2 で重み付けした残差から、横方向の漏れは隣のノードを含む
//! 2次式で近似する。ノード間は部分中性子流でつなぐ。
//!
//! 面では、一般化等価理論（GET。K. S. Smith, "Assembly homogenization techniques for light water
//! reactor analysis", Prog. Nucl. Energy 17, 303 (1986)）の条件で、不連続因子（ADF）f をかけた
//! 不均質な中性子束 f φ と正味の中性子流を連続にする。ADF は少数群定数の辺（左、右、上、下、−z、+z）ごとの値で、
//! x、y、z のどの面も同じ規則で使う（2次元の格子計算の定数の z の辺は 1。軸方向の反射体の定数は平板の格子計算の
//! 値。ADR-0128）。炉心の形の反射体の位置どうしの面も同じ（ADR-0132）。外周の真空の面は、格子計算の真空の辺の比
//! r = J/φ_het を不均質な中性子束に課す（ADR-0128）。
//!
//! ノードは直方体の格子の上に並べ、ノードのない位置（炉心の外形の外）を持てる。ノードのない位置に面する面は、
//! `Problem::outline` の条件の外周の面として扱う（ADR-0132）。

use xlc_format::fewgroup::FewGroup;
use xlc_numeric::anm::AnmModes;
use xlc_numeric::nem::{Beside, Inflow, Line, Outflow, is_fissile, outflow, transverse_moments};

/// 外周の境界条件。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// 正味の中性子流 0（対称面）
    Reflective,
    /// 境界での中性子束 0
    ZeroFlux,
    /// 真空: 不均質な中性子束 f φ と正味の中性子流に J = r f φ を課す。r は、その面のノードの物質の
    /// 少数群定数の、その面に当たる辺の `vacuum_ratio`（格子計算の真空の辺の J/φ_het。ADR-0128）
    Vacuum,
    /// Marshak の条件（不均質な入射の部分中性子流 0）。r = 1/2 の `Vacuum` で、P1 の関係を置くことを
    /// 入力で明示したときだけ使う
    Marshak,
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
    /// 炉心の形の反射体の位置に置いた物質の番号。f* の辺がほかの反射体の位置に向く（回転を取り違えた定数）ことを
    /// 見つけるのに使う（ADR-0095）。面の ADF には効かない（反射体どうしの面も各側の ADF を使う。ADR-0132）
    pub reflectors: Vec<usize>,
    /// 直方体の位置ごとの物質番号。並びは x が最も速く、次に y、z。ノードのない位置は `None`
    pub node_material: Vec<Option<usize>>,
    /// 直方体の外周の [-x, +x, -y, +y, -z, +z] の条件
    pub boundaries: [Boundary; 6],
    /// ノードのない位置に面する面の条件（ADR-0132）。ノードのない位置があるときだけ持つ
    pub outline: Option<Boundary>,
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
    /// [位置][群] のノード平均中性子束（全体の生成率が 1 になるよう規格化）。ノードのない位置は空
    pub flux: Vec<Vec<f64>>,
    /// 外周から外へ出る正味の中性子数（群ごと）
    pub leakage: Vec<f64>,
    /// [方向][位置][群] の [- 面, + 面] の正味の中性子流。各面の上の平均で、単位面積あたり、座標の正の向きを
    /// 正とし、`flux` と同じ規格化。ノードのない位置は空
    pub currents: [Vec<Vec<[f64; 2]>>; 3],
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
                return Err(NemError::NonFundamentalMode { node: core.site(n) });
            }
            let currents: [_; 3] = std::array::from_fn(|d| core.currents(&state, d));
            let solution = Solution {
                k,
                leakage: core.leakage(&currents),
                currents: currents.map(|c| core.on_sites(c)),
                flux: core.on_sites(state.flux.clone()),
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

/// 検査済みの問題と、ノードの位置の計算。ノードは、直方体の位置のうちノードのあるものだけを、直方体の並び
/// （x が最も速く、次に y、z）の順に詰めた番号で数える。
struct Core<'a> {
    problem: &'a Problem,
    /// 方向ごとのノード幅
    widths: [&'a [f64]; 3],
    /// 直方体の方向ごとの位置の数
    counts: [usize; 3],
    /// [ノード] の直方体の位置
    cells: Vec<[usize; 3]>,
    /// [直方体の位置] のノードの番号。ノードのない位置は `None`
    slots: Vec<Option<usize>>,
    /// [ノード] の物質番号
    material_of: Vec<usize>,
    /// ノードのない位置に面する面の条件（`Problem::outline`）
    outline: Boundary,
    groups: usize,
    /// [物質] の核分裂のある物質か（ノードの核を選ぶ。ADR-0102）
    fissile: Vec<bool>,
    /// [物質] の、炉心の形の反射体の位置に置いた物質か（`Problem::reflectors`）
    reflector: Vec<bool>,
    /// [ノード] の、核分裂のないノードの方向ごとの ANM のモード（定数と幅で決まるので1回だけ作る）
    anm: Vec<Option<[AnmModes; 3]>>,
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
        let sites: usize = counts.iter().product();
        if problem.node_material.len() != sites {
            return Err(format!(
                "{} node materials given for {sites} positions",
                problem.node_material.len()
            ));
        }
        if let Some(&m) = problem
            .reflectors
            .iter()
            .find(|&&m| m >= problem.materials.len())
        {
            return Err(format!(
                "reflector material {m} is out of the {} materials",
                problem.materials.len()
            ));
        }
        if let Some(m) = problem
            .node_material
            .iter()
            .flatten()
            .find(|&&m| m >= problem.materials.len())
        {
            return Err(format!("material {m} does not exist"));
        }
        let absent = problem.node_material.iter().any(Option::is_none);
        match (absent, problem.outline) {
            (true, None) => {
                return Err(
                    "positions without nodes need the condition of the faces next to them \
                     (the outline)"
                        .into(),
                );
            }
            (false, Some(_)) => {
                return Err("an outline condition is given, but every position has a node".into());
            }
            _ => {}
        }
        let mut cells = Vec::new();
        let mut slots = Vec::with_capacity(sites);
        let mut material_of = Vec::new();
        for (site, m) in problem.node_material.iter().enumerate() {
            slots.push(m.map(|m| {
                cells.push([
                    site % counts[0],
                    (site / counts[0]) % counts[1],
                    site / (counts[0] * counts[1]),
                ]);
                material_of.push(m);
                cells.len() - 1
            }));
        }
        if cells.is_empty() {
            return Err("no position has a node".into());
        }
        let groups = problem.materials.first().ok_or("no materials")?.groups();
        for (i, m) in problem.materials.iter().enumerate() {
            // 長さ、有限であること、不連続因子は少数群定数の検査（ファイルの読み込みと同じ）
            m.validate().map_err(|e| format!("material {i}: {e}"))?;
            if m.groups() != groups {
                return Err(format!("material {i} must have {groups} groups"));
            }
            // 核分裂のないノードは2群の ANM で解く（ADR-0102）
            if !is_fissile(&m.nu_fission).map_err(|e| format!("material {i}: {e}"))? && groups != 2
            {
                return Err(format!(
                    "material {i} has no fission and is solved with the 2-group analytic nodal \
                     method, but has {groups} groups"
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
        let fissile = problem
            .materials
            .iter()
            .map(|m| is_fissile(&m.nu_fission))
            .collect::<Result<Vec<bool>, _>>()
            .map_err(|e| e.to_string())?;
        let mut core = Self {
            problem,
            widths,
            counts,
            cells,
            slots,
            material_of,
            // ノードのない位置がなければ使わない（ある問題が外形の条件を持つことは上で確かめた）
            outline: problem.outline.unwrap_or(Boundary::Reflective),
            groups,
            fissile,
            reflector: (0..problem.materials.len())
                .map(|m| problem.reflectors.contains(&m))
                .collect(),
            anm: Vec::new(),
        };
        core.anm = (0..core.nodes())
            .map(|n| {
                let material = core.material_of[n];
                if core.fissile[material] {
                    return Ok(None);
                }
                let m = &problem.materials[material];
                let modes = |d: usize| {
                    AnmModes::new(&m.diffusion, &m.total, &m.scatter, core.width(n, d))
                        .map_err(|e| format!("material {material} in node {}: {e}", core.site(n)))
                };
                Ok(Some([modes(0)?, modes(1)?, modes(2)?]))
            })
            .collect::<Result<Vec<_>, String>>()?;
        core.check_heights()?;
        core.check_vacuum_ratios()?;
        core.check_normalized_adfs()?;
        Ok(core)
    }

    /// z の厚さ（`height`）を持つ定数を置いたノードの z の幅が、丸めの範囲でその厚さに等しいことを確かめる。
    /// 軸方向の定数の ADF と等価な D は、平板のノードの幅で解いた1ノードの問題の量なので、幅が違えば等価性が
    /// 崩れる（ADR-0128）
    fn check_heights(&self) -> Result<(), String> {
        // 層の境界の差（例えば 61.053 − 55.0）と格子の列の幅の和の丸め（相対 1e-15 程度）を許す
        const TOLERANCE: f64 = 1.0e-9;
        for n in 0..self.nodes() {
            let Some(height) = self.material(n).height else {
                continue;
            };
            let dz = self.width(n, 2);
            if (dz - height).abs() > TOLERANCE * dz.max(height) {
                return Err(format!(
                    "material {} is made for a height of {height} cm, but node {} is {dz} cm \
                     high; give the axial reflector layer one node of its height",
                    self.material_of[n],
                    self.site(n)
                ));
            }
        }
        Ok(())
    }

    /// f* の辺（`normalized_adf`）の面が、核分裂のある隣のノードの f* でない辺と接することを確かめる（ADR-0095）。
    /// f* は燃料の ADF で割った値なので、核分裂のない隣は、上下や左右を取り違えた定数である（ADR-0128）
    fn check_normalized_adfs(&self) -> Result<(), String> {
        for n in 0..self.nodes() {
            for d in 0..3 {
                for plus in [false, true] {
                    let side = side_of(d, plus);
                    if !self.material(n).normalized_adf[side] {
                        continue;
                    }
                    let material = self.material_of[n];
                    let node = self.site(n);
                    let Some(m) = self.neighbor(n, d, if plus { 1 } else { -1 }) else {
                        return Err(format!(
                            "material {material} has f* (an ADF over the fuel's) on side {side}, \
                             but node {node} has no neighbor there"
                        ));
                    };
                    let other = self.site(m);
                    if self.material(m).normalized_adf[side_of(d, !plus)] {
                        return Err(format!(
                            "nodes {node} and {other} both have f* on their shared face"
                        ));
                    }
                    // f* は燃料の ADF で割った値なので、径方向で反射体の位置に向く f* は、燃料に面する辺を取り違えた
                    // （回転を取り違えた）定数である
                    if d < 2 && self.reflector_face(n, m) {
                        return Err(format!(
                            "material {material} has f* on side {side}, but node {node} faces the \
                             reflector site of node {other} there; check the rotation of the reflector"
                        ));
                    }
                    if !self.fissile[self.material_of[m]] {
                        return Err(format!(
                            "material {material} has f* (an ADF over the fuel's) on side {side}, \
                             but its neighbor node {other} has no fission"
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    /// `Vacuum` の面（直方体の外周と外形。z の面も同じ）のノードの物質が、その面の辺の比を持つことを確かめる
    /// （ADR-0128）
    fn check_vacuum_ratios(&self) -> Result<(), String> {
        const NAMES: [&str; 6] = ["-x", "+x", "-y", "+y", "-z", "+z"];
        for n in 0..self.nodes() {
            for (face, name) in NAMES.iter().enumerate() {
                let (d, plus) = (face / 2, face % 2 == 1);
                if self.boundary(n, d, plus) != Some(Boundary::Vacuum) {
                    continue;
                }
                let side = side_of(d, plus);
                let material = self.material_of[n];
                if self.problem.materials[material].vacuum_ratio[side].is_none() {
                    return Err(format!(
                        "material {material} has no vacuum ratio J/phi on side {side} for the \
                         vacuum face {name} of node {}; make it in a lattice with that side \
                         vacuum, or use \"marshak\" for the P1 (Marshak) condition",
                        self.site(n)
                    ));
                }
            }
        }
        Ok(())
    }

    fn nodes(&self) -> usize {
        self.cells.len()
    }

    fn position(&self, n: usize) -> [usize; 3] {
        self.cells[n]
    }

    /// ノード n の直方体の位置の番号（`Problem::node_material` と `Solution` の並び）
    fn site(&self, n: usize) -> usize {
        let [nx, ny, _] = self.counts;
        let p = self.cells[n];
        p[0] + nx * (p[1] + ny * p[2])
    }

    /// [ノード] の値を [直方体の位置] の並びに置き直す。ノードのない位置は空
    fn on_sites<T>(&self, values: Vec<Vec<T>>) -> Vec<Vec<T>> {
        let mut out: Vec<Vec<T>> = (0..self.slots.len()).map(|_| Vec::new()).collect();
        for (n, v) in values.into_iter().enumerate() {
            out[self.site(n)] = v;
        }
        out
    }

    /// 方向 d に `step`（-1 か +1）だけ動いた直方体の位置。直方体の外なら `None`
    fn moved(&self, n: usize, d: usize, step: isize) -> Option<usize> {
        let mut p = self.position(n);
        p[d] = p[d]
            .checked_add_signed(step)
            .filter(|&i| i < self.counts[d])?;
        let [nx, ny, _] = self.counts;
        Some(p[0] + nx * (p[1] + ny * p[2]))
    }

    /// 方向 d の隣のノード（`step` が -1 なら - 側、+1 なら + 側）。直方体の外周か、ノードのない位置なら `None`。
    fn neighbor(&self, n: usize, d: usize, step: isize) -> Option<usize> {
        self.slots[self.moved(n, d, step)?]
    }

    /// ノード n の方向 d の - 面（`plus` が false）・+ 面（true）の境界条件。隣にノードがあれば `None`
    fn boundary(&self, n: usize, d: usize, plus: bool) -> Option<Boundary> {
        self.neighbor(n, d, if plus { 1 } else { -1 })
            .is_none()
            .then(|| self.edge(n, d, plus))
    }

    /// 隣にノードのない面の条件: 直方体の外周ならその面の条件、ノードのない位置に面すれば外形の条件
    fn edge(&self, n: usize, d: usize, plus: bool) -> Boundary {
        match self.moved(n, d, if plus { 1 } else { -1 }) {
            None => self.problem.boundaries[2 * d + usize::from(plus)],
            Some(_) => self.outline,
        }
    }

    fn width(&self, n: usize, d: usize) -> f64 {
        self.widths[d][self.position(n)[d]]
    }

    fn volume(&self, n: usize) -> f64 {
        (0..3).map(|d| self.width(n, d)).product()
    }

    fn material(&self, n: usize) -> &FewGroup {
        &self.problem.materials[self.material_of[n]]
    }

    /// ノード n の方向 d の - 面（`plus` が false）・+ 面（true）の不連続因子（辺の対応は `side_of`）。
    /// f* の辺（`normalized_adf`）は、f* × 隣のノードの同じ面の ADF（ADR-0095。隣があり、隣の辺が f* で
    /// ないことは `check_normalized_adfs` で確かめてある）。反射体の位置どうしの面も、燃料どうしの面と同じく
    /// 各側の ADF を使う（ADR-0132）
    fn adf(&self, n: usize, d: usize, plus: bool, g: usize) -> f64 {
        let side = side_of(d, plus);
        let neighbor = self.neighbor(n, d, if plus { 1 } else { -1 });
        let material = self.material(n);
        let own = material.adf[side][g];
        if !material.normalized_adf[side] {
            return own;
        }
        neighbor.map_or(f64::NAN, |m| {
            own * self.material(m).adf[side_of(d, !plus)][g]
        })
    }

    /// 隣り合うノード n と m が、どちらも炉心の形の反射体の位置か
    fn reflector_face(&self, n: usize, m: usize) -> bool {
        let at = |n: usize| self.reflector[self.material_of[n]];
        at(n) && at(m)
    }

    /// ノード n の方向 d の外周の面の、真空の面の比 r = J/φ_het。`Marshak` は 1/2
    fn vacuum_ratio(&self, n: usize, d: usize, plus: bool, g: usize) -> f64 {
        match self.edge(n, d, plus) {
            Boundary::Marshak => 0.5,
            // `Vacuum` の面の物質が辺の比を持つことは `check_vacuum_ratios` で確かめてある
            _ => self.material(n).vacuum_ratio[side_of(d, plus)]
                .as_deref()
                .map_or(f64::NAN, |r| r[g]),
        }
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
    /// - 真空: 格子計算の真空の辺の比 r で J = r f φ、つまり o - i = 2 r f (o + i) から
    ///   β = (1 - 2rf)/(1 + 2rf)（ADR-0128）
    /// - Marshak: r = 1/2 で、不均質な入射の部分中性子流 f φ/4 - (o - i)/2 = 0 と同じ。β = (1 - f)/(1 + f) は
    ///   隣の面の式で f' = 1、o' = 0 としたものに等しく、f = 1 なら 0
    fn face_inputs(&self, state: &State, n: usize, d: usize, g: usize) -> ((f64, f64), (f64, f64)) {
        let face = |step: isize| match self.neighbor(n, d, step) {
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
            None => match self.edge(n, d, step > 0) {
                Boundary::Reflective => (1.0, 0.0),
                Boundary::ZeroFlux => (-1.0, 0.0),
                Boundary::Vacuum | Boundary::Marshak => {
                    let plus = step > 0;
                    let rf = self.vacuum_ratio(n, d, plus, g) * self.adf(n, d, plus, g);
                    ((1.0 - 2.0 * rf) / (1.0 + 2.0 * rf), 0.0)
                }
            },
        };
        (face(-1), face(1))
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
        // 横方向の漏れは群のスイープの始めの値を使う。核分裂のないノードの ANM が使うほかの群の漏れは、その群の
        // 値がこのスイープで変わらないので、使うノードの分だけその場で求める
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
                // 入射は i = β o + 外からの入射（face_inputs）。[+ 面, - 面]
                let inflow_of = |group: usize| {
                    let ((beta_l, ext_l), (beta_r, ext_r)) = self.face_inputs(state, n, d, group);
                    [
                        Inflow {
                            beta: beta_r,
                            external: ext_r,
                        },
                        Inflow {
                            beta: beta_l,
                            external: ext_l,
                        },
                    ]
                };
                *part = match &self.anm[n] {
                    None => outflow(&Line {
                        diffusion: m.diffusion[g],
                        width: h,
                        removal,
                        coupling: &coupling,
                        coefficients: &state.coefficients[d][n],
                        group: g,
                        transverse: self.transverse_shape(state, &|i| leakage[i], n, d, g)?,
                        inflow: inflow_of(g),
                    }),
                    // 核分裂のないノードは2群の ANM（ADR-0102）。ほかの群はスイープのそのときの値を使う
                    Some(modes) => {
                        let transverse = (0..self.groups)
                            .map(|group| {
                                if group == g {
                                    self.transverse_shape(state, &|i| leakage[i], n, d, g)
                                } else {
                                    let other = |i: usize| -> [f64; 3] {
                                        std::array::from_fn(|e| {
                                            self.net_leakage(state, i, e, group)
                                        })
                                    };
                                    self.transverse_shape(state, &other, n, d, group)
                                }
                            })
                            .collect::<Result<Vec<_>, _>>()?;
                        let inflow: Vec<[Inflow; 2]> = (0..self.groups).map(inflow_of).collect();
                        modes[d]
                            .outflow(&transverse, &inflow, &state.flux[n], g)
                            .map_err(|e| {
                                NemError::InvalidProblem(format!("node {}: {e}", self.site(n)))
                            })?
                    }
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
        leakage: &dyn Fn(usize) -> [f64; 3],
        n: usize,
        d: usize,
        g: usize,
    ) -> Result<(f64, f64), NemError> {
        let transverse = |m: usize| -> f64 {
            let l = leakage(m);
            (0..3).filter(|&e| e != d).map(|e| l[e]).sum()
        };
        let h = self.width(n, d);
        let side = |step: isize| {
            Ok(match self.neighbor(n, d, step) {
                Some(m) => Beside::Node {
                    width: self.width(m, d),
                    leakage: transverse(m),
                },
                None => match self.edge(n, d, step > 0) {
                    Boundary::Reflective => Beside::Node {
                        width: h,
                        leakage: transverse(n),
                    },
                    Boundary::ZeroFlux => Beside::Face(0.0),
                    Boundary::Vacuum | Boundary::Marshak => {
                        let flux = state.flux[n][g];
                        if !(flux > 0.0 && flux.is_finite()) {
                            return Err(NemError::VacuumFaceFlux {
                                node: self.site(n),
                                flux,
                            });
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
        Ok(transverse_moments(h, side(-1)?, transverse(n), side(1)?))
    }

    /// 方向 d の [ノード][群] の、- 面と + 面の正味の中性子流（座標の正の向きを正）
    fn currents(&self, state: &State, d: usize) -> Vec<Vec<[f64; 2]>> {
        (0..self.nodes())
            .map(|n| {
                (0..self.groups)
                    .map(|g| {
                        let (in_minus, in_plus) = self.incoming(state, n, d, g);
                        [
                            in_minus - state.out_minus[d][n][g],
                            state.out_plus[d][n][g] - in_plus,
                        ]
                    })
                    .collect()
            })
            .collect()
    }

    /// 外周の面から外へ出る正味の中性子数（群ごと）。面の正味の中性子流（`currents`）の外周の面の和
    fn leakage(&self, currents: &[Vec<Vec<[f64; 2]>>; 3]) -> Vec<f64> {
        let mut total = vec![0.0; self.groups];
        for n in 0..self.nodes() {
            for (d, by_node) in currents.iter().enumerate() {
                let area = self.volume(n) / self.width(n, d);
                for (t, [minus, plus]) in total.iter_mut().zip(&by_node[n]) {
                    if self.neighbor(n, d, -1).is_none() {
                        *t -= area * minus;
                    }
                    if self.neighbor(n, d, 1).is_none() {
                        *t += area * plus;
                    }
                }
            }
        }
        total
    }
}

/// 方向 d の - 面（`plus` が false）・+ 面（true）に当たる少数群定数の辺（左、右、上、下、−z、+z）の番号。
/// -x が左、+x が右、-y が下、+y が上
fn side_of(d: usize, plus: bool) -> usize {
    match (d, plus) {
        (0, false) => 0,
        (0, true) => 1,
        (1, true) => 2,
        (1, false) => 3,
        (_, false) => 4,
        (_, true) => 5,
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

    /// どの位置にもノードを置く物質の並び
    fn placed(materials: impl IntoIterator<Item = usize>) -> Vec<Option<usize>> {
        materials.into_iter().map(Some).collect()
    }

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
            adf: vec![vec![1.0, 1.0]; 6],
            vacuum_ratio: vec![None; 6],
            normalized_adf: vec![false; 6],
            height: None,
            pins: None,
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
            reflectors: Vec::new(),
            node_material: placed(vec![0; nodes.iter().product()]),
            boundaries,
            outline: None,
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

    /// 解の面の正味の中性子流が、各ノードのバランスを（反復の収束の程度で）満たし、隣のノードと面で一致し、
    /// 反射の面で 0、外周の面の和が漏れに等しく、Marshak の面で J = f φ/4·2 = f (o + i) であること
    fn assert_face_currents(p: &Problem) {
        let (core, s, state) = solve_state(p);
        let ng = core.groups;
        let scale = s
            .currents
            .iter()
            .flatten()
            .flatten()
            .flatten()
            .fold(0.0, |a: f64, j| a.max(j.abs()));
        assert!(scale > 0.0);
        for n in 0..core.nodes() {
            let m = core.material(n);
            let phi = &s.flux[n];
            let fission: f64 = (0..ng).map(|g| m.nu_fission[g] * phi[g]).sum::<f64>() / s.k;
            for g in 0..ng {
                let leakage: f64 = (0..3)
                    .map(|d| {
                        let [minus, plus] = s.currents[d][n][g];
                        (plus - minus) / core.width(n, d)
                    })
                    .sum();
                let scattering_in: f64 = (0..ng)
                    .filter(|&h| h != g)
                    .map(|h| m.scatter[h][g] * phi[h])
                    .sum();
                let removal = (m.total[g] - m.scatter[g][g]) * phi[g];
                let residual = leakage + removal - scattering_in - m.chi[g] * fission;
                assert!(
                    residual.abs() < 1e-6 * removal,
                    "node {n} group {g}: {residual}"
                );
            }
        }
        let mut outer = vec![0.0; ng];
        for n in 0..core.nodes() {
            for d in 0..3 {
                let area = core.volume(n) / core.width(n, d);
                for (g, out) in outer.iter_mut().enumerate() {
                    let [minus, plus] = s.currents[d][n][g];
                    for (step, index, current) in [(-1, 0, -minus), (1, 1, plus)] {
                        if let Some(next) = core.neighbor(n, d, step) {
                            let other = s.currents[d][next][g][1 - index];
                            let own = if step > 0 { plus } else { minus };
                            assert!(
                                (own - other).abs() < 1e-12 * scale,
                                "node {n} direction {d}: {own} vs {other}"
                            );
                            continue;
                        }
                        *out += area * current;
                        match p.boundaries[2 * d + index] {
                            Boundary::Reflective => {
                                assert!(current.abs() < 1e-12 * scale, "reflective face of {n}");
                            }
                            Boundary::Marshak => {
                                let (o, i) = if step > 0 {
                                    let (_, i) = core.incoming(&state, n, d, g);
                                    (state.out_plus[d][n][g], i)
                                } else {
                                    let (i, _) = core.incoming(&state, n, d, g);
                                    (state.out_minus[d][n][g], i)
                                };
                                let f = core.adf(n, d, step > 0, g);
                                assert!(
                                    (current - f * (o + i)).abs() < 1e-12 * scale,
                                    "Marshak face of {n}"
                                );
                            }
                            _ => assert!(current > 0.0, "outer face of {n} must leak"),
                        }
                    }
                }
            }
        }
        for (g, (out, leak)) in outer.iter().zip(&s.leakage).enumerate() {
            assert!(
                (out - leak).abs() < 1e-12 * leak.abs().max(1e-300),
                "group {g}: {out} vs {leak}"
            );
        }
    }

    #[test]
    fn face_currents_balance_each_node() {
        // ピン出力の再構成（ADR-0096）は面の正味の中性子流を条件にする
        let mut boundaries = [Boundary::ZeroFlux; 6];
        boundaries[0] = Boundary::Reflective;
        boundaries[3] = Boundary::Reflective;
        assert_face_currents(&cuboid([60.0, 80.0, 100.0], [3, 4, 5], boundaries));
    }

    #[test]
    fn face_currents_are_continuous_across_adf_jumps() {
        // 燃料と反射体（核分裂がないので ANM。ADR-0102）の 2D で、ADF が面で違い、+x は Marshak。
        // 正味の中性子流は ADF の跳びによらず面で連続
        let mut fuel = material();
        fuel.adf = vec![
            vec![1.1, 0.9],
            vec![0.95, 1.2],
            vec![1.05, 0.85],
            vec![1.0, 1.1],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let mut water = reflector();
        water.adf = vec![
            vec![0.8, 1.3],
            vec![1.2, 0.7],
            vec![0.9, 1.4],
            vec![1.1, 0.75],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[1] = Boundary::Marshak;
        boundaries[2] = Boundary::ZeroFlux;
        let p = Problem {
            mesh: Mesh {
                dx: vec![10.0, 10.0, 10.0, 10.0],
                dy: vec![10.0, 10.0],
                dz: vec![20.0],
            },
            materials: vec![fuel, water],
            reflectors: Vec::new(),
            node_material: placed(vec![0, 0, 1, 1, 0, 1, 1, 1]),
            boundaries,
            outline: None,
        };
        assert_face_currents(&p);
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
        p.node_material[0] = Some(1);
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
            adf: vec![vec![1.0]; 6],
            vacuum_ratio: vec![None; 6],
            normalized_adf: vec![false; 6],
            height: None,
            pins: None,
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
            reflectors: Vec::new(),
            node_material: placed(node_material),
            boundaries,
            outline: None,
        }
    }

    /// x 方向に1列の問題（y と z は1ノード）を z 方向に並べ直した問題。x と z の幅、境界条件、少数群定数の
    /// 辺（左、右）と（−z、+z）を入れ替える。物質の配置の並びはそのまま（ノードの番号は x でも z でも列の順）
    fn to_z(p: &Problem) -> Problem {
        assert_eq!((p.mesh.dy.len(), p.mesh.dz.len()), (1, 1));
        fn swap<T>(sides: &mut [T]) {
            sides.swap(0, 4);
            sides.swap(1, 5);
        }
        let mut q = p.clone();
        q.mesh.dx = p.mesh.dz.clone();
        q.mesh.dz = p.mesh.dx.clone();
        swap(&mut q.boundaries);
        for m in &mut q.materials {
            swap(&mut m.adf);
            swap(&mut m.vacuum_ratio);
            swap(&mut m.normalized_adf);
        }
        q
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
            p.materials[p.node_material[n].unwrap()].adf[side][g]
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
        // ADF を入れる前（e3b627c）と同じ反復（a3, a4 を前の反復の a1, a2 から求める。#0025 の修正の前で、
        // 既定の許容差では e3b627c の結果に 1e-12 で一致していた）を、許容差を詰めて収束させた値。
        // 既定の許容差の値は反復の道筋で 1e-9 程度変わるので、収束させた値で比べる
        let p = two_material_slab(vec![vec![1.0, 1.0]; 6], vec![vec![1.0, 1.0]; 6]);
        let tight = Options {
            k_tolerance: 1e-14,
            flux_tolerance: 1e-13,
            max_iterations: 10_000,
        };
        let s = solve(&p, &tight).unwrap();
        let k = 8.510692798259939e-1;
        assert!((s.k - k).abs() < 1e-12, "{} vs {k}", s.k);
        let flux = [
            [0.0005911585762519943, 0.00012519241937033706],
            [0.001661405206417996, 0.0003510047208700676],
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
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let b = vec![
            vec![0.97, 0.8],
            vec![0.93, 0.85],
            vec![1.02, 0.75],
            vec![0.99, 0.9],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let p = two_material_slab(a.clone(), b.clone());
        assert_interface_conditions(&p);
        let unit = two_material_slab(vec![vec![1.0, 1.0]; 6], vec![vec![1.0, 1.0]; 6]);
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
        q.node_material = placed((0..12).map(|n| (n + n / 3) % 2));
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
        m_a.adf = vec![vec![f_a]; 6];
        let mut m_b = one_group();
        m_b.adf = vec![vec![f_b]; 6];
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
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let b = vec![
            vec![0.97, 0.8],
            vec![0.93, 0.85],
            vec![1.02, 0.75],
            vec![0.99, 0.9],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
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
        let p = slab(vec![5.0; 6], vec![m.clone()], vec![0; 6], Boundary::Marshak);
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
            k(Boundary::Marshak),
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
            m.adf = vec![vec![f]; 6];
            let b = marshak_buckling(&m, l, f);
            let expected = one_group_k(&m, b * b);
            let p = slab(vec![l / 8.0; 8], vec![m], vec![0; 8], Boundary::Marshak);
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
        let p = slab(vec![l / 8.0; 8], vec![m], vec![0; 8], Boundary::Marshak);
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
            let mut p = cuboid(size, [n; 3], [Boundary::Marshak; 6]);
            p.materials = vec![one_group()];
            (solve(&p, &Options::default()).unwrap().k - expected).abs() / expected
        };
        let (coarse, fine) = (error(4), error(8));
        assert!(fine < 1e-4, "relative error {fine}");
        assert!(fine < coarse / 8.0, "{coarse} -> {fine}");
    }

    /// 外周の x の2つの面（辺の左と右）に、真空の面の比 r を群ごとに持たせる
    fn with_x_ratios(mut m: FewGroup, r: &[f64]) -> FewGroup {
        m.vacuum_ratio[0] = Some(r.to_vec());
        m.vacuum_ratio[1] = Some(r.to_vec());
        m
    }

    #[test]
    fn half_ratio_vacuum_equals_marshak() {
        // r = 1/2 の真空の面の条件 J = r f φ は、Marshak の条件 f φ/4 = J/2 と同じ（ADR-0128）
        let mut m = material();
        m.adf[0] = vec![1.3, 0.8];
        m.adf[1] = vec![0.7, 1.6];
        let marshak = solve(
            &slab(vec![5.0; 6], vec![m.clone()], vec![0; 6], Boundary::Marshak),
            &Options::default(),
        )
        .unwrap();
        let p = slab(
            vec![5.0; 6],
            vec![with_x_ratios(m, &[0.5, 0.5])],
            vec![0; 6],
            Boundary::Vacuum,
        );
        let vacuum = solve(&p, &Options::default()).unwrap();
        assert!(
            (vacuum.k - marshak.k).abs() < 1e-14,
            "{} vs {}",
            vacuum.k,
            marshak.k
        );
        for (a, b) in vacuum
            .flux
            .iter()
            .flatten()
            .zip(marshak.flux.iter().flatten())
        {
            assert!((a - b).abs() < 1e-12 * b, "{a} vs {b}");
        }
    }

    #[test]
    fn vacuum_face_imposes_the_ratio() {
        // 真空の面では、不均質な中性子束 f φ = 2f(o + i) と正味の中性子流 o - i が J = r f φ を満たす
        let mut m = material();
        m.adf[0] = vec![1.3, 0.8];
        m.adf[1] = vec![0.7, 1.6];
        let r = [0.62, 0.55];
        let p = slab(
            vec![5.0; 6],
            vec![with_x_ratios(m.clone(), &r)],
            vec![0; 6],
            Boundary::Vacuum,
        );
        let (core, _, state) = solve_state(&p);
        for (g, ratio) in r.into_iter().enumerate() {
            let (in_minus, _) = core.incoming(&state, 0, 0, g);
            let (_, in_plus) = core.incoming(&state, 5, 0, g);
            let (out_minus, out_plus) = (state.out_minus[0][0][g], state.out_plus[0][5][g]);
            for (f, o, i) in [
                (m.adf[0][g], out_minus, in_minus),
                (m.adf[1][g], out_plus, in_plus),
            ] {
                let (current, heterogeneous) = (o - i, 2.0 * f * (o + i));
                assert!(current > 0.0);
                assert!(
                    (current - ratio * heterogeneous).abs() < 1e-14 * current,
                    "{current} vs {}",
                    ratio * heterogeneous
                );
            }
        }
    }

    #[test]
    fn vacuum_ratio_matches_analytic() {
        // 1群の裸の平板に J = r f φ を課した解 φ = cos(Bx) は D B tan(Bl/2) = r f を満たす。
        // r が大きいほど漏れが増えて k は下がる
        let l = 40.0;
        for f in [0.8, 1.3] {
            let mut previous = f64::INFINITY;
            for r in [0.5, 0.55, 0.6, 0.9] {
                let mut m = one_group();
                m.adf = vec![vec![f]; 6];
                let d = m.diffusion[0];
                let b = smallest_root(
                    |b| (b * l / 2.0).sin() * d * b - r * f * (b * l / 2.0).cos(),
                    PI / l,
                );
                let expected = one_group_k(&m, b * b);
                let p = slab(
                    vec![l / 8.0; 8],
                    vec![with_x_ratios(m, &[r])],
                    vec![0; 8],
                    Boundary::Vacuum,
                );
                // z の面も同じ規則（ADR-0128）
                for p in [p.clone(), to_z(&p)] {
                    let s = solve(&p, &Options::default()).unwrap();
                    assert!(
                        (s.k - expected).abs() < 1e-8 * expected,
                        "f = {f}, r = {r}: {} vs {expected}",
                        s.k
                    );
                }
                assert!(expected < previous, "f = {f}, r = {r}");
                previous = expected;
            }
        }
    }

    #[test]
    fn vacuum_needs_ratios() {
        // 物質がその面の辺の比を持たなければエラー。黙って Marshak の条件にしない
        let p = slab(vec![5.0; 6], vec![material()], vec![0; 6], Boundary::Vacuum);
        let err = solve(&p, &Options::default()).unwrap_err().to_string();
        assert!(err.contains("ratio") && err.contains("marshak"), "{err}");
        // 片側だけ持つ
        let mut m = material();
        m.vacuum_ratio[0] = Some(vec![0.5, 0.5]);
        let p = slab(vec![5.0; 6], vec![m], vec![0; 6], Boundary::Vacuum);
        assert!(solve(&p, &Options::default()).is_err());
        // z の面も、その面の辺（+z）の比を使う（ADR-0128）。径方向の辺がすべて比を持っていても、+z の辺が
        // 持たなければエラー
        let mut m = material();
        m.vacuum_ratio = (0..6).map(|s| (s < 4).then(|| vec![0.5, 0.5])).collect();
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[5] = Boundary::Vacuum;
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], boundaries);
        p.materials = vec![m];
        let err = solve(&p, &Options::default()).unwrap_err().to_string();
        assert!(err.contains("side 5") && err.contains("marshak"), "{err}");
        p.materials[0].vacuum_ratio[5] = Some(vec![0.5, 0.5]);
        assert!(solve(&p, &Options::default()).is_ok());
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
        // 方向ごとの（- 面, + 面）の辺（左、右、上、下、−z、+z の番号）。z の面も少数群定数の辺の ADF を使う
        // （ADR-0128）
        let expected = analytic(f_a, f_b);
        for (d, (minus, plus)) in [(0, 1), (3, 2), (4, 5)].into_iter().enumerate() {
            let bulk = one_group();
            let mut edge_a = one_group();
            let mut edge_b = one_group();
            edge_a.adf = vec![vec![3.0]; 6];
            edge_a.adf[minus] = vec![1.0];
            edge_a.adf[plus] = vec![f_a];
            edge_b.adf = vec![vec![0.5]; 6];
            edge_b.adf[minus] = vec![f_b];
            edge_b.adf[plus] = vec![1.0];
            let mut widths = [vec![20.0], vec![20.0], vec![20.0]];
            widths[d] = vec![5.0; 16];
            let mut boundaries = [Boundary::Reflective; 6];
            boundaries[2 * d] = Boundary::ZeroFlux;
            boundaries[2 * d + 1] = Boundary::ZeroFlux;
            let [dx, dy, dz] = widths;
            let p = Problem {
                mesh: Mesh { dx, dy, dz },
                materials: vec![bulk, edge_a, edge_b],
                reflectors: Vec::new(),
                node_material: placed((0..16).map(|i| match i {
                    5 => 1,
                    6 => 2,
                    _ => 0,
                })),
                boundaries,
                outline: None,
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
        a.adf = vec![vec![1.1, 1.3]; 6];
        let mut b = material();
        b.diffusion = vec![1.2, 0.5];
        b.absorption = vec![0.012, 0.10];
        b.nu_fission = vec![0.004, 0.09];
        b.adf = vec![vec![0.95, 0.8]; 6];
        let k = |n: usize| {
            let mut boundaries = [Boundary::Reflective; 6];
            boundaries[1] = Boundary::Marshak;
            boundaries[3] = Boundary::Marshak;
            let p = Problem {
                mesh: Mesh {
                    dx: vec![40.0 / n as f64; n],
                    dy: vec![40.0 / n as f64; n],
                    dz: vec![20.0],
                },
                materials: vec![a.clone(), b.clone()],
                reflectors: Vec::new(),
                node_material: placed(
                    (0..n * n).map(|i| usize::from((i % n < n / 2) != (i / n < n / 2))),
                ),
                boundaries,
                outline: None,
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
        let p = slab(
            vec![5.0; 6],
            vec![material()],
            vec![0; 6],
            Boundary::Marshak,
        );
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
    /// まとまり）ごとに反応率を保存して定数を均質化し、GET の1ノードの問題（ADR-0103）で ADF を作った
    /// 粗いノードの問題と、細かいノードの解を返す。z は1層。
    ///
    /// 1ノードの問題の外周の扱い（`Edge`）は、粗いノードの面が外周なら細かいノードの問題の境界、内側の面
    /// なら隣の粗いノード（幅と参照の横方向の漏れ。ADR-0095）とする
    fn coarsen_with_get(fine: &Problem, blocks: [usize; 2]) -> (Problem, Solution) {
        use xlc_numeric::nem::{Edge, OneNode, surface_fluxes};
        let (core, solution, state) = solve_state(fine);
        let core = &core;
        let [nx, ny, nz] = core.counts;
        assert_eq!(nz, 1);
        let [bx, by] = blocks;
        let (cx, cy) = (nx / bx, ny / by);
        assert_eq!((cx * bx, cy * by), (nx, ny));
        let ng = core.groups;
        // [粗いノード] の定数、ノード平均、面（-x, +x, -y, +y）の不均質な中性子束と外向きの正味の中性子流、幅
        let mut nodes = Vec::new();
        for (cj, ci) in (0..cy).flat_map(|j| (0..cx).map(move |i| (j, i))) {
            let members: Vec<usize> = (cj * by..(cj + 1) * by)
                .flat_map(|j| {
                    (ci * bx..(ci + 1) * bx).map(move |i| core.slots[i + nx * j].unwrap())
                })
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
                adf: vec![vec![1.0; ng]; 6],
                vacuum_ratio: vec![None; 6],
                normalized_adf: vec![false; 6],
                height: None,
                pins: None,
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
            if produced > 0.0 {
                for c in &mut m.chi {
                    *c /= produced;
                }
            }
            let average: Vec<f64> = flux.iter().map(|f| f / volume).collect();

            // 面（-x, +x, -y, +y）の不均質な中性子束と外向きの正味の中性子流の、面にわたる平均
            let mut het = [vec![0.0; ng], vec![0.0; ng], vec![0.0; ng], vec![0.0; ng]];
            let mut current = het.clone();
            for (face, d, plus) in [(0, 0, false), (1, 0, true), (2, 1, false), (3, 1, true)] {
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
            let widths: [f64; 2] = [
                core.widths[0][ci * bx..(ci + 1) * bx].iter().sum(),
                core.widths[1][cj * by..(cj + 1) * by].iter().sum(),
            ];
            nodes.push((m, average, het, current, widths));
        }
        // [粗いノード][方向] の横方向の漏れ（もう一方の方向の正味の漏れ）
        let transverse: Vec<[Vec<f64>; 2]> = nodes
            .iter()
            .map(|(_, _, _, current, widths)| {
                let net = |d: usize| -> Vec<f64> {
                    (0..ng)
                        .map(|g| (current[2 * d][g] + current[2 * d + 1][g]) / widths[d])
                        .collect()
                };
                [net(1), net(0)]
            })
            .collect();
        let mut materials = Vec::new();
        for (c, (mut m, average, het, current, widths)) in nodes.clone().into_iter().enumerate() {
            let (ci, cj) = (c % cx, c / cx);
            let edge = |d: usize, plus: bool| {
                let neighbor = match (d, plus) {
                    (0, false) => (ci > 0).then(|| c - 1),
                    (0, true) => (ci + 1 < cx).then_some(c + 1),
                    (_, false) => (cj > 0).then(|| c - cx),
                    (_, true) => (cj + 1 < cy).then_some(c + cx),
                };
                match (neighbor, fine.boundaries[2 * d + usize::from(plus)]) {
                    (Some(other), _) => Edge::Neighbor {
                        width: nodes[other].4[d],
                        leakage: &transverse[other][d],
                    },
                    (None, Boundary::Marshak) => Edge::Vacuum,
                    (None, Boundary::Reflective) => Edge::Reflective,
                    (None, Boundary::ZeroFlux) => panic!("no zero-flux edge in GET"),
                    (None, Boundary::Vacuum) => panic!("the fine reference uses marshak"),
                }
            };
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
                // 参照の真空の面の比 J/φ_het（ADR-0128。参照は Marshak の条件なので 1/2）
                let (d, plus) = (face / 2, face % 2 == 1);
                if edge(d, plus) == Edge::Vacuum {
                    m.vacuum_ratio[side] =
                        Some((0..ng).map(|g| current[face][g] / het[face][g]).collect());
                }
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
            reflectors: Vec::new(),
            node_material: placed(0..materials.len()),
            materials,
            // 粗いノードの真空の面は、参照から作った比 J/φ_het を使う
            boundaries: fine.boundaries.map(|b| match b {
                Boundary::Marshak => Boundary::Vacuum,
                other => other,
            }),
            outline: None,
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
            m.adf = vec![vec![1.0; m.groups()]; 6];
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
        boundaries[1] = Boundary::Marshak;
        let fine = Problem {
            mesh: Mesh {
                dx: vec![2.5; 16],
                dy: vec![20.0],
                dz: vec![20.0],
            },
            materials: three_fuels(),
            reflectors: Vec::new(),
            node_material: placed(vec![0, 1, 0, 1, 1, 0, 1, 0, 2, 2, 1, 1, 2, 2, 0, 0]),
            boundaries,
            outline: None,
        };
        let (coarse, solution) = coarsen_with_get(&fine, [8, 1]);
        assert_get_equivalence(&coarse, &solution, [8, 1], 16);
    }

    #[test]
    fn get_adfs_reproduce_a_node_with_vacuum_sides() {
        // 20 cm 角の粗いノード1つ（細かいノード 5 cm 角 × 16、3種類の物質）。-x と -y は反射、+x と +y は
        // 真空。真空の面での横方向の漏れの扱いと、2方向の横方向の漏れを通る
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[1] = Boundary::Marshak;
        boundaries[3] = Boundary::Marshak;
        let fine = Problem {
            mesh: Mesh {
                dx: vec![5.0; 4],
                dy: vec![5.0; 4],
                dz: vec![20.0],
            },
            materials: three_fuels(),
            reflectors: Vec::new(),
            node_material: placed(vec![0, 1, 0, 2, 1, 0, 2, 0, 0, 2, 1, 1, 2, 0, 1, 0]),
            boundaries,
            outline: None,
        };
        let (coarse, solution) = coarsen_with_get(&fine, [4, 4]);
        assert_get_equivalence(&coarse, &solution, [4, 4], 4);
    }

    /// 核分裂のない2群の物質（水に近い。上方散乱あり）
    fn reflector() -> FewGroup {
        let mut m = material();
        m.diffusion = vec![1.25, 0.30];
        m.total = vec![0.60, 1.90];
        m.scatter = vec![vec![0.565, 0.032], vec![0.0015, 1.885]];
        m.absorption = vec![0.003, 0.0135];
        m.nu_fission = vec![0.0; 2];
        m.fission = vec![0.0; 2];
        m.chi = vec![0.0; 2];
        m
    }

    #[test]
    fn get_adfs_reproduce_fuel_and_a_wide_reflector() {
        // 燃料（20 cm）と反射体（21.5 cm。細かいノード 2.5 cm × 8 と 2.6875 cm × 8）の 1D、+x は真空。
        // 粗い反射体のノードは core と1ノードの問題の両方で ANM で解き（ADR-0102）、GET の定数で細かいノードの
        // 解を再現する
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[1] = Boundary::Marshak;
        let mut materials = three_fuels();
        materials.push(reflector());
        let mut dx = vec![2.5; 8];
        dx.extend(vec![21.5 / 8.0; 8]);
        let fine = Problem {
            mesh: Mesh {
                dx,
                dy: vec![20.0],
                dz: vec![20.0],
            },
            materials,
            reflectors: Vec::new(),
            node_material: placed(vec![0, 1, 0, 1, 1, 0, 2, 0, 3, 3, 3, 3, 3, 3, 3, 3]),
            boundaries,
            outline: None,
        };
        let (coarse, solution) = coarsen_with_get(&fine, [8, 1]);
        assert!(coarse.materials[1].nu_fission.iter().all(|&v| v == 0.0));
        assert_get_equivalence(&coarse, &solution, [8, 1], 16);
    }

    #[test]
    fn non_fissile_nodes_need_two_groups() {
        // ANM の核は2群（ADR-0102）
        let mut m = reflector();
        m.bounds = vec![1.0e7, 1.0e-5];
        m.total = vec![0.6];
        m.absorption = vec![0.01];
        m.diffusion = vec![1.0];
        m.scatter = vec![vec![0.59]];
        m.nu_fission = vec![0.0];
        m.fission = vec![0.0];
        m.chi = vec![0.0];
        m.adf = vec![vec![1.0]; 6];
        let mut fuel = m.clone();
        fuel.nu_fission = vec![0.02];
        fuel.fission = vec![0.008];
        fuel.chi = vec![1.0];
        let p = slab(
            vec![10.0; 2],
            vec![fuel, m],
            vec![0, 1],
            Boundary::Reflective,
        );
        let err = solve(&p, &Options::default());
        assert!(matches!(err, Err(NemError::InvalidProblem(_))), "{err:?}");
    }

    #[test]
    fn negative_nu_fission_is_an_error() {
        // νΣf で核の種類（NEM か ANM）を選ぶので、負の νΣf は推測せずにエラー（ADR-0102）
        let mut m = material();
        m.nu_fission[0] = -1.0e-3;
        let p = slab(vec![10.0; 2], vec![m], vec![0, 0], Boundary::Reflective);
        let err = solve(&p, &Options::default());
        assert!(matches!(err, Err(NemError::InvalidProblem(_))), "{err:?}");
    }

    #[test]
    fn get_adfs_reproduce_inner_faces_with_transverse_leakage() {
        // 2×2 の粗いノード（細かいノード 2.5 cm × 4 × 4 ずつ）。+x と +y が真空なので、内側の面にも横方向の
        // 漏れがあり、1ノードの問題の内側の辺は隣の粗いノードの幅と漏れで2次式を作る（ADR-0095）
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[1] = Boundary::Marshak;
        boundaries[3] = Boundary::Marshak;
        let pattern = [0, 1, 2, 0, 1, 1, 0, 2];
        let fine = Problem {
            mesh: Mesh {
                dx: vec![2.5; 8],
                dy: vec![2.5; 8],
                dz: vec![20.0],
            },
            materials: three_fuels(),
            reflectors: Vec::new(),
            node_material: placed((0..64).map(|n| pattern[(n * 5 + n / 8) % 8])),
            boundaries,
            outline: None,
        };
        let (coarse, solution) = coarsen_with_get(&fine, [4, 4]);
        assert_get_equivalence(&coarse, &solution, [4, 4], 8);
    }

    #[test]
    fn normalized_adfs_take_the_neighbor_adf() {
        // 反射体の定数のように、辺の ADF を f* = f_自分 / f_隣 で持ち、印を付けた辺では、core はその面の
        // ADF を f* × 隣のノードの面の ADF にする（ADR-0095）。隣の ADF が作ったときと同じなら、元の ADF の
        // 解と同じ
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[1] = Boundary::Marshak;
        let fine = Problem {
            mesh: Mesh {
                dx: vec![2.5; 16],
                dy: vec![20.0],
                dz: vec![20.0],
            },
            materials: three_fuels(),
            reflectors: Vec::new(),
            node_material: placed(vec![0, 1, 0, 1, 1, 0, 1, 0, 2, 2, 1, 1, 2, 2, 0, 0]),
            boundaries,
            outline: None,
        };
        let (coarse, _) = coarsen_with_get(&fine, [8, 1]);
        let raw = solve(&coarse, &Options::default()).unwrap();
        let mut normalized = coarse.clone();
        let fuel = normalized.materials[0].adf[1].clone();
        let reflector = &mut normalized.materials[1];
        reflector.adf[0] = reflector.adf[0]
            .iter()
            .zip(&fuel)
            .map(|(f, g)| f / g)
            .collect();
        reflector.normalized_adf[0] = true;
        let s = solve(&normalized, &Options::default()).unwrap();
        assert!((s.k - raw.k).abs() < 1e-12 * raw.k, "{} vs {}", s.k, raw.k);
        for (a, b) in s.flux.iter().flatten().zip(raw.flux.iter().flatten()) {
            assert!((a - b).abs() < 1e-10 * b, "{a} vs {b}");
        }
        // 隣の ADF を変えると f* はそれに掛かるが、面の条件は両側の ADF の比だけで決まるので、解は変わらない
        // （f* の方法が比を反射体の問題から取る理由）
        let mut other = normalized.clone();
        other.materials[0].adf[1] = fuel.iter().map(|f| 1.1 * f).collect();
        let k = solve(&other, &Options::default()).unwrap().k;
        assert!((k - raw.k).abs() < 1e-12 * raw.k, "{k} vs {}", raw.k);
        // z の面の f* も同じ（ADR-0128）
        let s = solve(&to_z(&normalized), &Options::default()).unwrap();
        assert!((s.k - raw.k).abs() < 1e-12 * raw.k, "{} vs {}", s.k, raw.k);
        for (a, b) in s.flux.iter().flatten().zip(raw.flux.iter().flatten()) {
            assert!((a - b).abs() < 1e-10 * b, "{a} vs {b}");
        }
    }

    #[test]
    fn normalized_adfs_need_a_fissile_neighbor() {
        // f* は燃料の ADF で割った値なので、隣が核分裂のない物質ならエラー（上下を取り違えた軸方向の反射体の
        // 定数を見つける。ADR-0128）。x でも z でも同じ
        let mut star = reflector();
        star.normalized_adf[0] = true;
        let p = slab(
            vec![10.0; 3],
            vec![material(), star, reflector()],
            vec![0, 1, 1],
            Boundary::ZeroFlux,
        );
        for p in [p.clone(), to_z(&p)] {
            let err = solve(&p, &Options::default()).unwrap_err().to_string();
            assert!(err.contains("fission"), "{err}");
        }
        // 燃料に向いていれば解ける
        let mut q = p.clone();
        q.node_material = placed(vec![0, 1, 2]);
        for q in [q.clone(), to_z(&q)] {
            assert!(solve(&q, &Options::default()).is_ok());
        }
    }

    #[test]
    fn normalized_adfs_toward_assembly_layers_take_unit_adfs() {
        // 反射体の位置の f* の辺が、x と y の面で集合体の位置の核分裂のない層（端栓、プレナム）に向くときは、両側の
        // ADF を 1 にする（ADR-0132。I 型の定数は燃料と組で作ったもので、核分裂のない層と組で作った定数はない）。
        // 下の行は燃料と反射体、上の行は核分裂のない集合体の層と反射体
        let mut star = reflector();
        star.normalized_adf[0] = true;
        star.adf[0] = vec![1.3, 0.6];
        let mut layer = reflector();
        layer.adf[1] = vec![0.8, 1.4];
        let problem = |star: FewGroup, layer: FewGroup, reflectors: Vec<usize>| {
            let mut boundaries = [Boundary::Reflective; 6];
            boundaries[1] = Boundary::Marshak;
            plane(
                vec![20.0, 21.5],
                vec![20.0, 20.0],
                vec![material(), star, layer],
                vec![Some(0), Some(1), Some(2), Some(1)],
                boundaries,
                None,
            )
            .with_reflectors(reflectors)
        };
        let p = problem(star.clone(), layer.clone(), vec![1]);
        // 上の行の面だけ両側を 1 にした定数（下の行の反射体は f* のまま）
        let mut unit_layer = layer.clone();
        unit_layer.adf[1] = vec![1.0, 1.0];
        let mut unit_star = star.clone();
        unit_star.normalized_adf[0] = false;
        unit_star.adf[0] = vec![1.0, 1.0];
        let mut q = problem(star.clone(), unit_layer, vec![1]);
        q.materials.push(unit_star);
        q.node_material[3] = Some(3);
        let (s, t) = (
            solve(&p, &Options::default()).unwrap(),
            solve(&q, &Options::default()).unwrap(),
        );
        assert!((s.k - t.k).abs() < 1e-12 * t.k, "{} vs {}", s.k, t.k);
        for (a, b) in s.flux.iter().flatten().zip(t.flux.iter().flatten()) {
            assert!((a - b).abs() < 1e-10 * b, "{a} vs {b}");
        }
        // 核分裂のない層が反射体の位置なら、回転を取り違えた定数としてエラー
        let err = solve(
            &problem(star.clone(), layer.clone(), vec![1, 2]),
            &Options::default(),
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("rotation"), "{err}");
        // f* の側が反射体の位置でなければ（燃料の位置に置いた反射体の定数）、今のとおりエラー
        let err = solve(&problem(star, layer, vec![]), &Options::default())
            .unwrap_err()
            .to_string();
        assert!(err.contains("fission"), "{err}");
    }

    #[test]
    fn heights_must_match_the_node() {
        // z の厚さを持つ（軸方向の）定数は、z の幅が丸めの範囲で同じノードにだけ置ける（ADR-0128）
        let mut p = cuboid([30.0, 30.0, 30.0], [1, 1, 3], [Boundary::Reflective; 6]);
        p.materials[0].height = Some(10.0);
        assert!(solve(&p, &Options::default()).is_ok());
        p.materials[0].height = Some(10.0 * (1.0 + 1e-14));
        assert!(solve(&p, &Options::default()).is_ok());
        for h in [9.9, 20.0] {
            p.materials[0].height = Some(h);
            let err = solve(&p, &Options::default()).unwrap_err().to_string();
            assert!(err.contains("height"), "{err}");
        }
    }

    /// x 方向の問題を y 方向に並べ直す（x と y の幅、境界条件、少数群定数の辺（左、右）と（下、上）を入れ替える）
    fn to_y(p: &Problem) -> Problem {
        let mut q = p.clone();
        q.mesh.dx = p.mesh.dy.clone();
        q.mesh.dy = p.mesh.dx.clone();
        q.boundaries.swap(0, 2);
        q.boundaries.swap(1, 3);
        // 辺は（左、右、上、下、−z、+z）: 新しい下は元の左、新しい上は元の右、新しい左右は元の下と上。z はそのまま
        fn turn<T: Clone>(sides: &[T]) -> Vec<T> {
            [3, 2, 1, 0, 4, 5]
                .iter()
                .map(|&s| sides[s].clone())
                .collect()
        }
        for m in &mut q.materials {
            m.adf = turn(&m.adf);
            m.vacuum_ratio = turn(&m.vacuum_ratio);
            m.normalized_adf = turn(&m.normalized_adf);
        }
        q
    }

    #[test]
    fn faces_between_reflector_sites_keep_their_adfs() {
        // 炉心の形の反射体の位置に置いた物質（`reflectors`）どうしの x と y の面も、燃料どうしの面と同じく各側の
        // 少数群定数の ADF を使う（ADR-0132）
        let mut a = reflector();
        a.adf = vec![
            vec![1.3, 0.8],
            vec![0.7, 1.6],
            vec![1.2, 0.9],
            vec![0.85, 1.4],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let mut b = reflector();
        b.adf = vec![
            vec![1.1, 0.6],
            vec![0.9, 1.2],
            vec![1.05, 0.75],
            vec![0.95, 1.3],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let x = |reflectors: Vec<usize>, a: &FewGroup, b: &FewGroup| {
            let mut p = slab(
                vec![20.0, 21.5, 21.5],
                vec![material(), a.clone(), b.clone()],
                vec![0, 1, 2],
                Boundary::Reflective,
            );
            p.boundaries[1] = Boundary::Marshak;
            p.reflectors = reflectors;
            p
        };
        let unit = |m: &FewGroup, side: usize| {
            let mut m = m.clone();
            m.adf[side] = vec![1.0; 2];
            m
        };
        let same = |p: &Problem, q: &Problem| {
            let (s, t) = (
                solve(p, &Options::default()).unwrap(),
                solve(q, &Options::default()).unwrap(),
            );
            let close = (s.k - t.k).abs() < 1e-12 * t.k
                && s.flux
                    .iter()
                    .flatten()
                    .zip(t.flux.iter().flatten())
                    .all(|(u, v)| (u - v).abs() < 1e-10 * v);
            (close, s.k, t.k)
        };
        for turn in [false, true] {
            let along = |p: Problem| if turn { to_y(&p) } else { p };
            // a の右（+x）と b の左（−x）の面
            let marked = along(x(vec![1, 2], &a, &b));
            let plain = along(x(vec![], &a, &b));
            let (close, k, k_plain) = same(&marked, &plain);
            assert!(close, "turn {turn}: {k} vs {k_plain}");
            // ADF は効いている（1 にすれば解が変わる）
            let unit = along(x(vec![1, 2], &unit(&a, 1), &unit(&b, 0)));
            assert!(!same(&marked, &unit).0, "turn {turn}");
        }
        // 範囲の外の番号はエラー
        let err = solve(&x(vec![3], &a, &b), &Options::default());
        assert!(matches!(err, Err(NemError::InvalidProblem(_))), "{err:?}");
        // f* の辺は燃料の ADF で割った値なので、反射体の位置に向いていれば（回転を取り違えた定数）エラー
        let mut turned = a.clone();
        turned.normalized_adf[1] = true;
        let err = solve(&x(vec![1, 2], &turned, &b), &Options::default());
        assert!(matches!(err, Err(NemError::InvalidProblem(_))), "{err:?}");
    }

    #[test]
    fn normalized_adfs_need_a_plain_neighbor() {
        // f* の辺は隣のノードの ADF を掛けるので、外周の面と、隣も f* の面はエラー
        let mut p = two_material_slab(vec![vec![1.0, 1.0]; 6], vec![vec![1.0, 1.0]; 6]);
        p.materials[0].normalized_adf[0] = true;
        assert!(matches!(
            solve(&p, &Options::default()),
            Err(NemError::InvalidProblem(_))
        ));
        let mut p = two_material_slab(vec![vec![1.0, 1.0]; 6], vec![vec![1.0, 1.0]; 6]);
        p.materials[0].normalized_adf[1] = true;
        p.materials[1].normalized_adf[0] = true;
        let err = solve(&p, &Options::default());
        assert!(matches!(err, Err(NemError::InvalidProblem(_))), "{err:?}");
        // 片側だけなら解ける
        p.materials[0].normalized_adf[1] = false;
        assert!(solve(&p, &Options::default()).is_ok());
    }

    /// z が1層（20 cm、上下は反射）の nx × ny の格子の問題。`sites` は y の下の行から、各行は x の小さい方から
    fn plane(
        dx: Vec<f64>,
        dy: Vec<f64>,
        materials: Vec<FewGroup>,
        sites: Vec<Option<usize>>,
        boundaries: [Boundary; 6],
        outline: Option<Boundary>,
    ) -> Problem {
        Problem {
            mesh: Mesh {
                dx,
                dy,
                dz: vec![20.0],
            },
            materials,
            reflectors: Vec::new(),
            node_material: sites,
            boundaries,
            outline,
        }
    }

    /// どの辺にも真空の比と、辺ごとに違う ADF を持つ反射体
    fn sided_reflector() -> FewGroup {
        let mut r = reflector();
        // 辺（左、右、上、下、−z、+z）
        r.adf = vec![
            vec![1.1, 0.8],
            vec![1.3, 0.7],
            vec![0.9, 1.2],
            vec![1.2, 0.6],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        r.vacuum_ratio = vec![
            Some(vec![0.45, 0.40]),
            Some(vec![0.50, 0.35]),
            Some(vec![0.42, 0.38]),
            Some(vec![0.47, 0.33]),
            None,
            None,
        ];
        r
    }

    #[test]
    fn outline_faces_match_the_cuboid_faces() {
        // ノードのない位置に面する面は、同じノードを直方体の外周に置いて同じ条件にした面と同じ解になる（ADR-0132）。
        // ノードのない位置の側の直方体の外周の条件は効かない
        let materials = vec![material(), sided_reflector()];
        let h = 21.5;
        for condition in [
            Boundary::Vacuum,
            Boundary::Marshak,
            Boundary::ZeroFlux,
            Boundary::Reflective,
        ] {
            // （直方体だけの問題、ノードのない位置を足した問題、比べる位置の対応）
            let mut x_plus = [Boundary::Reflective; 6];
            x_plus[1] = condition;
            let mut y_plus = [Boundary::Reflective; 6];
            y_plus[3] = condition;
            let mut x_minus = [Boundary::Reflective; 6];
            x_minus[0] = condition;
            // ノードのない位置の外の、直方体の外周の面だけを変える
            let far = |face: usize| {
                let mut b = [Boundary::Reflective; 6];
                b[face] = Boundary::ZeroFlux;
                b
            };
            let cases = [
                (
                    plane(
                        vec![h; 2],
                        vec![h],
                        materials.clone(),
                        placed([0, 1]),
                        x_plus,
                        None,
                    ),
                    plane(
                        vec![h; 3],
                        vec![h],
                        materials.clone(),
                        vec![Some(0), Some(1), None],
                        far(1),
                        Some(condition),
                    ),
                    vec![0, 1],
                ),
                (
                    plane(
                        vec![h],
                        vec![h; 2],
                        materials.clone(),
                        placed([0, 1]),
                        y_plus,
                        None,
                    ),
                    plane(
                        vec![h],
                        vec![h; 3],
                        materials.clone(),
                        vec![Some(0), Some(1), None],
                        far(3),
                        Some(condition),
                    ),
                    vec![0, 1],
                ),
                (
                    plane(
                        vec![h; 2],
                        vec![h],
                        materials.clone(),
                        placed([1, 0]),
                        x_minus,
                        None,
                    ),
                    plane(
                        vec![h; 3],
                        vec![h],
                        materials.clone(),
                        vec![None, Some(1), Some(0)],
                        far(0),
                        Some(condition),
                    ),
                    vec![1, 2],
                ),
            ];
            for (cuboid, stepped, sites) in cases {
                let a = solve(&cuboid, &Options::default()).unwrap();
                let b = solve(&stepped, &Options::default()).unwrap();
                assert!(
                    (a.k - b.k).abs() < 1e-12 * a.k,
                    "{condition:?}: {} vs {}",
                    a.k,
                    b.k
                );
                for (n, &site) in sites.iter().enumerate() {
                    for (x, y) in a.flux[n].iter().zip(&b.flux[site]) {
                        assert!((x - y).abs() < 1e-10 * x, "{condition:?}: {x} vs {y}");
                    }
                }
                let empty = (0..b.flux.len()).find(|s| !sites.contains(s)).unwrap();
                assert!(b.flux[empty].is_empty());
                assert!(b.currents.iter().all(|c| c[empty].is_empty()));
                for (x, y) in a.leakage.iter().zip(&b.leakage) {
                    assert!((x - y).abs() < 1e-10 * x.abs().max(1e-12), "{x} vs {y}");
                }
            }
        }
    }

    /// 4 × 4 の集合体の位置の四隅にノードを置かない、燃料だけの階段の炉心。各位置を `split` × `split` に分ける
    fn stepped_core(split: usize, outline: Boundary) -> Problem {
        let n = 4 * split;
        let corner = |i: usize| i / split == 0 || i / split == 3;
        let sites = (0..n * n)
            .map(|s| (!(corner(s % n) && corner(s / n))).then_some(0))
            .collect();
        let mut boundaries = [outline; 6];
        boundaries[4] = Boundary::Reflective;
        boundaries[5] = Boundary::Reflective;
        plane(
            vec![20.0 / split as f64; n],
            vec![20.0 / split as f64; n],
            vec![material()],
            sites,
            boundaries,
            Some(outline),
        )
    }

    #[test]
    fn stepped_core_keeps_its_symmetry_and_balance() {
        // 4 回の回転の対称性を持つ階段の炉心の中性子束は同じ対称性を持ち、生成 = 吸収 + 外周と外形からの漏れ
        let p = stepped_core(1, Boundary::Marshak);
        let s = solve(&p, &Options::default()).unwrap();
        for (site, phi) in s.flux.iter().enumerate() {
            let (i, j) = (site % 4, site / 4);
            let turned = &s.flux[(3 - j) + 4 * i];
            assert_eq!(phi.len(), turned.len());
            for (a, b) in phi.iter().zip(turned) {
                assert!((a - b).abs() < 1e-6 * a, "site {site}: {a} vs {b}");
            }
        }
        let m = material();
        let volume = 20.0 * 20.0 * 20.0;
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
    fn stepped_core_converges_with_refinement() {
        // 外形の中性子束ゼロの面は、ノードを細かくすると k が収束する（外形の凹んだ角を含む）。凹んだ角の特異性
        // のため、直方体（20 cm のノードで差分法と 3 × 10⁻⁴ で一致）より遅く、20 cm のノードでは k が差分法の
        // 細かいメッシュの解（0.8602。セル中心の差分、0.25 cm）より 2.5% 大きい。四隅を強い吸収体のノードで埋めた
        // 直方体も、同じ 20 cm のノードで 0.8838 と同じ程度に外れるので、外形の面の扱いでなく NEM の凹んだ角の
        // 近似の誤差である（ADR-0132）
        let k = |split: usize| {
            solve(
                &stepped_core(split, Boundary::ZeroFlux),
                &Options::default(),
            )
            .unwrap()
            .k
        };
        let ks = [k(1), k(2), k(4), k(8)];
        for w in ks.windows(3) {
            assert!((w[2] - w[1]).abs() < 0.7 * (w[1] - w[0]).abs(), "{ks:?}");
        }
        // 細かくすると、差分法の解に上から近づく
        assert!(ks.windows(2).all(|w| w[1] < w[0]), "{ks:?}");
        assert!((ks[3] - 0.8602).abs() < 3e-3 * 0.8602, "{ks:?}");
    }

    #[test]
    fn outline_needs_a_condition_and_ratios() {
        let h = 21.5;
        // ノードのない位置があるのに外形の条件がない
        let p = plane(
            vec![h; 3],
            vec![h],
            vec![material(), sided_reflector()],
            vec![Some(0), Some(1), None],
            [Boundary::Reflective; 6],
            None,
        );
        let err = solve(&p, &Options::default()).unwrap_err().to_string();
        assert!(err.contains("outline"), "{err}");
        // ノードのない位置がないのに外形の条件がある
        let p = plane(
            vec![h; 2],
            vec![h],
            vec![material(), sided_reflector()],
            placed([0, 1]),
            [Boundary::Reflective; 6],
            Some(Boundary::Vacuum),
        );
        let err = solve(&p, &Options::default()).unwrap_err().to_string();
        assert!(err.contains("outline"), "{err}");
        // 外形の真空の面のノードが、その辺の真空の比を持たない
        let mut r = sided_reflector();
        r.vacuum_ratio[1] = None;
        let p = plane(
            vec![h; 3],
            vec![h],
            vec![material(), r],
            vec![Some(0), Some(1), None],
            [Boundary::Reflective; 6],
            Some(Boundary::Vacuum),
        );
        let err = solve(&p, &Options::default()).unwrap_err().to_string();
        assert!(err.contains("vacuum ratio"), "{err}");
        // ノードが1つもない
        let p = plane(
            vec![h],
            vec![h],
            vec![material()],
            vec![None],
            [Boundary::Reflective; 6],
            Some(Boundary::Vacuum),
        );
        assert!(solve(&p, &Options::default()).is_err());
    }

    impl Problem {
        fn with_reflectors(mut self, reflectors: Vec<usize>) -> Self {
            self.reflectors = reflectors;
            self
        }
    }
}
