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
//! 値。ADR-0128）。炉心の形の反射体の位置どうしの面も同じ（ADR-0132）。ただし、反射体の位置の f* の辺が集合体の
//! 位置の核分裂のない軸方向の層に向く x と y の面は、両側の ADF を 1 とする（ADR-0132）。外周の真空の面は、格子計算の
//! 真空の辺の比 r = J/φ_het を不均質な中性子束に課す（ADR-0128）。
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
    /// 見つけるのと（ADR-0095）、f* の辺が集合体の位置の核分裂のない軸方向の層に向く面の ADF を両側 1 にするのに使う
    /// （ADR-0132）。反射体どうしの面は各側の ADF を使う（ADR-0132）
    pub reflectors: Vec<usize>,
    /// 直方体の位置ごとの物質番号。並びは x が最も速く、次に y、z。ノードのない位置は `None`
    pub node_material: Vec<Option<usize>>,
    /// 直方体の外周の [-x, +x, -y, +y, -z, +z] の条件
    pub boundaries: [Boundary; 6],
    /// ノードのない位置に面する面の条件（ADR-0132）。ノードのない位置があるときだけ持つ
    pub outline: Option<Boundary>,
}

impl Mesh {
    /// 方向ごとの位置の数
    pub fn counts(&self) -> [usize; 3] {
        [self.dx.len(), self.dy.len(), self.dz.len()]
    }

    /// 位置 p の番号（x が最も速く、次に y、z。`Problem::node_material` と `Solution` の並び）
    pub fn site(&self, p: [usize; 3]) -> usize {
        let [nx, ny, _] = self.counts();
        p[0] + nx * (p[1] + ny * p[2])
    }

    /// 番号 `site` の位置
    pub fn position(&self, site: usize) -> [usize; 3] {
        let [nx, ny, _] = self.counts();
        [site % nx, (site / nx) % ny, site / (nx * ny)]
    }

    /// 位置 p から方向 d に `step`（-1 か +1）だけ動いた位置。直方体の外なら `None`
    pub fn moved(&self, mut p: [usize; 3], d: usize, step: isize) -> Option<[usize; 3]> {
        p[d] = p[d]
            .checked_add_signed(step)
            .filter(|&i| i < self.counts()[d])?;
        Some(p)
    }
}

impl Problem {
    /// 位置 p の方向 d の - 面（`plus` が false）・+ 面（true）の境界条件。隣の位置にノードがあれば `None`。
    /// 直方体の外周ならその面の条件、ノードのない位置に面すれば外形の条件（ADR-0132）。外形の条件がなければエラー
    pub fn face_condition(
        &self,
        p: [usize; 3],
        d: usize,
        plus: bool,
    ) -> Result<Option<Boundary>, String> {
        match self.mesh.moved(p, d, if plus { 1 } else { -1 }) {
            None => Ok(Some(self.boundaries[2 * d + usize::from(plus)])),
            Some(q) if self.node_material[self.mesh.site(q)].is_some() => Ok(None),
            Some(_) => self
                .outline
                .map(Some)
                .ok_or_else(|| "positions without nodes need the outline condition".into()),
        }
    }
}

/// 直方体の位置ごとの核分裂率 Σ_g Σf,g φ_g V（`solution` の規格化。ノードのない位置は 0）。物質の番号は NEM の解で確かめてある
pub fn fission_rates(problem: &Problem, solution: &Solution) -> Vec<f64> {
    problem
        .node_material
        .iter()
        .zip(&solution.flux)
        .enumerate()
        .map(|(site, (material, flux))| {
            let Some(m) = material.map(|m| &problem.materials[m]) else {
                return 0.0;
            };
            let [x, y, z] = problem.mesh.position(site);
            let volume = problem.mesh.dx[x] * problem.mesh.dy[y] * problem.mesh.dz[z];
            m.fission
                .iter()
                .zip(flux)
                .map(|(f, phi)| f * phi)
                .sum::<f64>()
                * volume
        })
        .collect()
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
    let mut core = Core::new(problem).map_err(NemError::InvalidProblem)?;
    iterate(&mut core, options, None).map(|(solution, _)| solution)
}

/// 外側の反復のたびに、いくつかの物質の定数を作り直すもの（制御棒の部分挿入のノード。ADR-0127）。
pub trait Refresh {
    /// 定数を作るのに要るノードの直方体の位置（`Problem::node_material` の並び）
    fn sites(&self) -> &[usize];
    /// そのときの k と、`sites` の順のノードの量から、置き換える（物質の番号, 定数）の組。置き換えた物質は、
    /// 元の物質と群の数、核分裂の有無、f* の辺、真空の辺、z の厚さが同じであること
    fn refresh(&mut self, k: f64, nodes: &[NodeSnapshot])
    -> Result<Vec<(usize, FewGroup)>, String>;
    /// 置き換える物質の z の面で、隣の f* の辺（`normalized_adf`）が掛ける ADF の（物質の番号, [−z、+z][群]）。
    /// 置き換えた定数の z の面の不連続因子は、ノードの中を混ぜたことによる不連続を表し、f*（燃料の格子計算の ADF
    /// で割った値。ADR-0095）の相手の ADF ではないので、`Core::adf` はこの値を掛ける
    fn f_star_partners(&self) -> Vec<(usize, [Vec<f64>; 2])>;
}

/// `Refresh` に渡す、外側の反復の始めのノードの量。
#[derive(Debug, Clone, PartialEq)]
pub struct NodeSnapshot {
    /// [群] のノード平均の中性子束
    pub flux: Vec<f64>,
    /// [群] の、z の向きの式の横方向の漏れ（x と y の正味の漏れの和、単位体積あたり）の、ノード平均と f1、f2 の
    /// 係数（NEM の z の向きの式と同じ2次式）
    pub transverse: Vec<[f64; 3]>,
}

/// `solve` に、外側の反復（全群のスイープ）の前ごとに `refresh` で物質を置き換えることを加えたもの。解と、最後に
/// 置き換えたあとの物質の定数（`Problem::materials` の並び）を返す。
#[tracing::instrument(name = "nem", skip_all)]
pub fn solve_with(
    problem: &Problem,
    options: &Options,
    refresh: &mut dyn Refresh,
) -> Result<(Solution, Vec<FewGroup>), NemError> {
    let mut core = Core::new(problem).map_err(NemError::InvalidProblem)?;
    for (material, adfs) in refresh.f_star_partners() {
        core.set_f_star_partner(material, adfs)
            .map_err(NemError::InvalidProblem)?;
    }
    let (solution, _) = iterate(&mut core, options, Some(refresh))?;
    Ok((solution, core.materials))
}

/// べき乗法の反復。収束した解と、そのときの部分中性子流などの状態を返す。`refresh` があれば、反復の始めごとに
/// 物質を置き換える。
fn iterate(
    core: &mut Core,
    options: &Options,
    mut refresh: Option<&mut dyn Refresh>,
) -> Result<(Solution, State), NemError> {
    let ng = core.groups;
    let nodes = core.nodes();
    let mut state = State::new(nodes, ng);
    let mut k = 1.0;
    let production = |core: &Core, state: &State| -> f64 {
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
    let initial = production(core, &state);
    if initial <= 0.0 {
        return Err(NemError::NoFission);
    }
    state.scale(1.0 / initial);
    let mut k_history = Remaining::default();
    let mut flux_history = Remaining::default();

    // `Refresh::sites` のノードの番号
    let refreshed: Vec<usize> = match &refresh {
        Some(r) => {
            let node_of: std::collections::HashMap<usize, usize> =
                (0..nodes).map(|n| (core.site(n), n)).collect();
            r.sites()
                .iter()
                .map(|site| {
                    node_of.get(site).copied().ok_or_else(|| {
                        NemError::InvalidProblem(format!("position {site} has no node"))
                    })
                })
                .collect::<Result<_, _>>()?
        }
        None => Vec::new(),
    };

    for iteration in 1..=options.max_iterations {
        let previous = state.flux.clone();
        if let Some(r) = refresh.as_deref_mut() {
            let snapshots = refreshed
                .iter()
                .map(|&n| core.snapshot(&state, n))
                .collect::<Result<Vec<_>, _>>()?;
            for (material, constants) in
                r.refresh(k, &snapshots).map_err(NemError::InvalidProblem)?
            {
                core.replace_material(material, constants)
                    .map_err(NemError::InvalidProblem)?;
            }
        }
        for g in 0..ng {
            core.sweep_group(&mut state, g, k)?;
        }
        // べき乗法: 生成率の比で k を更新し、全体を生成率 1 に規格化する
        let rate = production(core, &state);
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
    /// [物質] の定数。部分挿入のノードの物質は、外側の反復のたびに `Refresh` が置き換える（ADR-0127）
    materials: Vec<FewGroup>,
    /// 方向ごとのノード幅
    widths: [&'a [f64]; 3],
    /// [ノード] の直方体の位置
    cells: Vec<[usize; 3]>,
    /// [ノード] の物質番号
    material_of: Vec<usize>,
    /// [ノード][方向][- 面、+ 面] の隣（隣のノードか、境界の条件。`Problem::face_condition`）
    faces: Vec<[[Face; 2]; 3]>,
    /// [ノード][方向][- 面、+ 面] の、ADF を 1 にする面か（`unpaired_face`）
    unit_adf: Vec<[[bool; 2]; 3]>,
    /// [物質] の、隣の f* の辺が z の面で掛ける [−z、+z][群] の ADF（`Refresh::f_star_partners`）。`None` は
    /// その物質の z の辺の ADF
    f_star_partner: Vec<Option<[Vec<f64>; 2]>>,
    groups: usize,
    /// [物質] の核分裂のある物質か（ノードの核を選ぶ。ADR-0102）
    fissile: Vec<bool>,
    /// [物質] の、炉心の形の反射体の位置に置いた物質か（`Problem::reflectors`）
    reflector: Vec<bool>,
    /// [ノード] の、核分裂のないノードの方向ごとの ANM のモード（定数と幅で決まるので1回だけ作る）
    anm: Vec<Option<[AnmModes; 3]>>,
}

/// ノードの面の隣。
#[derive(Debug, Clone, Copy, PartialEq)]
enum Face {
    /// 隣のノード
    Node(usize),
    /// 直方体の外周か、ノードのない位置に面する面の境界条件
    Edge(Boundary),
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
        let sites: usize = mesh.counts().iter().product();
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
        // ノードのない位置は径方向の図の位置（ADR-0132）なので、z の面には現れない
        let layer = mesh.dx.len() * mesh.dy.len();
        if let Some(site) = (layer..sites).find(|&site| {
            problem.node_material[site].is_none() != problem.node_material[site % layer].is_none()
        }) {
            return Err(format!(
                "position {site} has a node at some heights only; positions without nodes must be \
                 the same at every height"
            ));
        }
        let mut cells = Vec::new();
        let mut slots = Vec::with_capacity(sites);
        let mut material_of = Vec::new();
        for (site, m) in problem.node_material.iter().enumerate() {
            slots.push(m.map(|m| {
                cells.push(mesh.position(site));
                material_of.push(m);
                cells.len() - 1
            }));
        }
        if cells.is_empty() {
            return Err("no position has a node".into());
        }
        let faces = cells
            .iter()
            .map(|&p| {
                let face = |d: usize, plus: bool| -> Result<Face, String> {
                    Ok(match problem.face_condition(p, d, plus)? {
                        Some(b) => Face::Edge(b),
                        // 隣にノードがある（`face_condition` が確かめた）
                        None => mesh
                            .moved(p, d, if plus { 1 } else { -1 })
                            .and_then(|q| slots[mesh.site(q)])
                            .map(Face::Node)
                            .ok_or("a face with a neighbor has no neighbor node")?,
                    })
                };
                let mut out = [[Face::Edge(Boundary::Reflective); 2]; 3];
                for (d, sides) in out.iter_mut().enumerate() {
                    *sides = [face(d, false)?, face(d, true)?];
                }
                Ok(out)
            })
            .collect::<Result<Vec<_>, String>>()?;
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
            materials: problem.materials.clone(),
            widths,
            cells,
            material_of,
            faces,
            unit_adf: Vec::new(),
            f_star_partner: vec![None; problem.materials.len()],
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
                let m = &core.materials[material];
                let modes = |d: usize| {
                    AnmModes::new(&m.diffusion, &m.total, &m.scatter, core.width(n, d))
                        .map_err(|e| format!("material {material} in node {}: {e}", core.site(n)))
                };
                Ok(Some([modes(0)?, modes(1)?, modes(2)?]))
            })
            .collect::<Result<Vec<_>, String>>()?;
        core.check_heights()?;
        core.unit_adf = (0..core.nodes())
            .map(|n| {
                std::array::from_fn(|d| {
                    std::array::from_fn(|side| match core.faces[n][d][side] {
                        Face::Node(m) => core.unpaired_face(n, m, d),
                        Face::Edge(_) => false,
                    })
                })
            })
            .collect();
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
                    if !self.fissile[self.material_of[m]] && !self.unit_adf[n][d][usize::from(plus)]
                    {
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
                if self.materials[material].vacuum_ratio[side].is_none() {
                    return Err(format!(
                        "material {material} has no vacuum ratio J/phi on side {side} for the \
                         vacuum face {name} of node {}; make it in a lattice with that side \
                         vacuum, mark its layer (plant or inventory) with vacuum = \"marshak\", or use \
                         \"marshak\" for the face for the P1 (Marshak) condition",
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

    /// ノード n の `NodeSnapshot`
    fn snapshot(&self, state: &State, n: usize) -> Result<NodeSnapshot, NemError> {
        let transverse = (0..self.groups)
            .map(|g| {
                let leakage = |i: usize| -> [f64; 3] {
                    std::array::from_fn(|d| self.net_leakage(state, i, d, g))
                };
                let (l1, l2) = self.transverse_shape(state, &leakage, n, 2, g)?;
                let own = leakage(n);
                Ok([own[0] + own[1], l1, l2])
            })
            .collect::<Result<_, NemError>>()?;
        Ok(NodeSnapshot {
            flux: state.flux[n].clone(),
            transverse,
        })
    }

    /// 物質 `material` の定数を置き換える（`Refresh`）。`Core::new` の検査の前提（群の数、核分裂の有無、f* の辺、
    /// 真空の辺、z の厚さ）を変える定数はエラー。核分裂のないノードの ANM のモードは作り直す
    fn replace_material(&mut self, material: usize, constants: FewGroup) -> Result<(), String> {
        let old = self
            .materials
            .get(material)
            .ok_or_else(|| format!("material {material} does not exist"))?;
        constants
            .validate()
            .map_err(|e| format!("material {material}: {e}"))?;
        let vacuum = |m: &FewGroup| {
            m.vacuum_ratio
                .iter()
                .map(Option::is_some)
                .collect::<Vec<_>>()
        };
        let fissile = is_fissile(&constants.nu_fission).map_err(|e| e.to_string())?;
        if constants.groups() != self.groups
            || fissile != self.fissile[material]
            || constants.normalized_adf != old.normalized_adf
            || vacuum(&constants) != vacuum(old)
            || constants.height != old.height
        {
            return Err(format!(
                "material {material}: the refreshed constants change the groups, the fission, the f* \
                 sides, the vacuum sides or the height"
            ));
        }
        if constants.diffusion.iter().any(|&d| d.is_nan() || d <= 0.0)
            || (0..self.groups).any(|g| {
                let removal = constants.total[g] - constants.scatter[g][g];
                removal.is_nan() || removal <= 0.0
            })
        {
            return Err(format!(
                "material {material}: diffusion and removal cross sections must be positive"
            ));
        }
        self.materials[material] = constants;
        for n in 0..self.nodes() {
            if self.material_of[n] != material || self.anm[n].is_none() {
                continue;
            }
            let m = &self.materials[material];
            let modes = |d: usize| {
                AnmModes::new(&m.diffusion, &m.total, &m.scatter, self.width(n, d))
                    .map_err(|e| format!("material {material} in node {}: {e}", self.site(n)))
            };
            self.anm[n] = Some([modes(0)?, modes(1)?, modes(2)?]);
        }
        Ok(())
    }

    fn position(&self, n: usize) -> [usize; 3] {
        self.cells[n]
    }

    /// ノード n の直方体の位置の番号（`Problem::node_material` と `Solution` の並び）
    fn site(&self, n: usize) -> usize {
        self.problem.mesh.site(self.cells[n])
    }

    /// [ノード] の値を [直方体の位置] の並びに置き直す。ノードのない位置は空
    fn on_sites<T>(&self, values: Vec<Vec<T>>) -> Vec<Vec<T>> {
        let mut out: Vec<Vec<T>> = (0..self.problem.node_material.len())
            .map(|_| Vec::new())
            .collect();
        for (n, v) in values.into_iter().enumerate() {
            out[self.site(n)] = v;
        }
        out
    }

    /// ノード n の方向 d の - 面（`plus` が false）・+ 面（true）の隣
    fn face(&self, n: usize, d: usize, plus: bool) -> Face {
        self.faces[n][d][usize::from(plus)]
    }

    /// 方向 d の隣のノード（`step` が -1 なら - 側、+1 なら + 側）。直方体の外周か、ノードのない位置なら `None`。
    fn neighbor(&self, n: usize, d: usize, step: isize) -> Option<usize> {
        match self.face(n, d, step > 0) {
            Face::Node(m) => Some(m),
            Face::Edge(_) => None,
        }
    }

    /// ノード n の方向 d の - 面・+ 面の境界条件。隣にノードがあれば `None`
    fn boundary(&self, n: usize, d: usize, plus: bool) -> Option<Boundary> {
        match self.face(n, d, plus) {
            Face::Node(_) => None,
            Face::Edge(b) => Some(b),
        }
    }

    fn width(&self, n: usize, d: usize) -> f64 {
        self.widths[d][self.position(n)[d]]
    }

    fn volume(&self, n: usize) -> f64 {
        (0..3).map(|d| self.width(n, d)).product()
    }

    fn material(&self, n: usize) -> &FewGroup {
        &self.materials[self.material_of[n]]
    }

    /// ノード n の方向 d の - 面（`plus` が false）・+ 面（true）の不連続因子（辺の対応は `side_of`）。
    /// f* の辺（`normalized_adf`）は、f* × 隣のノードの同じ面の ADF（ADR-0095。隣があり、隣の辺が f* で
    /// ないことは `check_normalized_adfs` で確かめてある）。z の面で隣の物質に `f_star_partner` があれば、隣の
    /// ADF の代わりにその値を掛ける（ADR-0127）。反射体の位置どうしの面も、燃料どうしの面と同じく
    /// 各側の ADF を使う（ADR-0132）。反射体の位置の f* の辺が集合体の位置の核分裂のない層に向く面は 1
    /// （`unpaired_face`）
    fn adf(&self, n: usize, d: usize, plus: bool, g: usize) -> f64 {
        let side = side_of(d, plus);
        if self.unit_adf[n][d][usize::from(plus)] {
            return 1.0;
        }
        let neighbor = self.neighbor(n, d, if plus { 1 } else { -1 });
        let material = self.material(n);
        let own = material.adf[side][g];
        if !material.normalized_adf[side] {
            return own;
        }
        neighbor.map_or(f64::NAN, |m| {
            let facing = side_of(d, !plus);
            let partner = match &self.f_star_partner[self.material_of[m]] {
                Some(z) if d == 2 => z[usize::from(!plus)][g],
                _ => self.material(m).adf[facing][g],
            };
            own * partner
        })
    }

    /// 物質 `material` の z の面で、隣の f* の辺が掛ける [−z、+z][群] の ADF を `adfs` にする
    fn set_f_star_partner(&mut self, material: usize, adfs: [Vec<f64>; 2]) -> Result<(), String> {
        if material >= self.materials.len()
            || adfs
                .iter()
                .any(|a| a.len() != self.groups || a.iter().any(|&f| !(f > 0.0 && f.is_finite())))
        {
            return Err(format!(
                "material {material}: the ADFs for f* neighbors must exist and be positive in every \
                 group, got {adfs:?}"
            ));
        }
        self.f_star_partner[material] = Some(adfs);
        Ok(())
    }

    /// 方向 d で隣り合うノード n と m の面が、反射体の位置の f* の辺と、集合体の位置の核分裂のない軸方向の層（端栓、
    /// プレナム。z の厚さを持つスラブの定数。ADR-0128）の x か y の面か。f* は燃料の ADF で割った値で、核分裂のない
    /// 層と組で作った定数はないので、この面は両側の ADF を 1 とする（ADR-0132。近似）
    fn unpaired_face(&self, n: usize, m: usize, d: usize) -> bool {
        if d >= 2 {
            return false;
        }
        let toward = |a: usize, b: usize| {
            let plus = self.position(b)[d] > self.position(a)[d];
            self.reflector[self.material_of[a]]
                && !self.reflector[self.material_of[b]]
                && !self.fissile[self.material_of[b]]
                && self.material(b).height.is_some()
                && self.material(a).normalized_adf[side_of(d, plus)]
        };
        toward(n, m) || toward(m, n)
    }

    /// 隣り合うノード n と m が、どちらも炉心の形の反射体の位置か
    fn reflector_face(&self, n: usize, m: usize) -> bool {
        let at = |n: usize| self.reflector[self.material_of[n]];
        at(n) && at(m)
    }

    /// ノード n の方向 d の、条件 `boundary` の外周の面の、真空の面の比 r = J/φ_het。`Marshak` は 1/2
    fn vacuum_ratio(&self, n: usize, d: usize, plus: bool, boundary: Boundary, g: usize) -> f64 {
        match boundary {
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
        let face = |step: isize| match self.face(n, d, step > 0) {
            Face::Node(m) => {
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
            Face::Edge(boundary) => match boundary {
                Boundary::Reflective => (1.0, 0.0),
                Boundary::ZeroFlux => (-1.0, 0.0),
                Boundary::Vacuum | Boundary::Marshak => {
                    let plus = step > 0;
                    let rf = self.vacuum_ratio(n, d, plus, boundary, g) * self.adf(n, d, plus, g);
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
            Ok(match self.face(n, d, step > 0) {
                Face::Node(m) => Beside::Node {
                    width: self.width(m, d),
                    leakage: transverse(m),
                },
                Face::Edge(boundary) => match boundary {
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
mod tests;
