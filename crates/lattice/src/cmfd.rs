//! CMFD（粗メッシュ有限差分）による MOC の加速。
//!
//! 粗メッシュは格子のセルとする（格子の周りの余白も、薄い行と列のセルになる）。
//! スイープでセルの面を通る正味の中性子の流れを数え、
//! 拡散の結合係数を流れに合うよう補正した粗メッシュの固有値問題を解いて、
//! 細かい領域の中性子束をセルごとの比で更新する。

use crate::material::Material;
use crate::moc::Remaining;
use crate::tracking::Tracking;
use xlc_numeric::lu::Lu;

/// 格子のセルを粗メッシュとし、トラックが外周で通る面を前もって求めたもの。
///
/// 面は、縦の面（x 方向に通る）を先に、横の面（y 方向に通る）を後に番号を付ける。
/// 縦の面 (r, i) はセルの行 r の、左から i 番目の格子線の上（i = 0..=列の数）、横の面 (j, c) は
/// セルの列 c の、上から j 番目の格子線の上（j = 0..=行の数）にある。面の正の向きは、
/// 縦の面では +x（右）、横の面では -y（下、行の番号が増える向き）とする。
///
/// 内側の面の通過は前もって持たず、スイープで隣り合う区間のセルを比べて求める。前もって持つと、
/// 区間の数に比例するメモリと、スイープで読む量が増えるため。
#[derive(Debug, Clone, PartialEq)]
pub struct CoarseMesh {
    /// 列（左から）ごとのセルの幅 (cm)
    widths: Vec<f64>,
    /// 行（上から）ごとのセルの高さ (cm)
    heights: Vec<f64>,
    /// [領域] の属するセル
    cell_of_region: Vec<usize>,
    /// [トラック] の始点で入る面と終点で出る面（順方向に通るときの向き付き）
    ends: Vec<[(usize, f64); 2]>,
}

/// 面の両側のセル。面の正の向きは `low` から `high` へ。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FaceCells {
    /// 内側の面
    Inner { low: usize, high: usize },
    /// 外周の面。`outward` は外向きが面の正の向きなら 1、逆なら -1
    Outer { cell: usize, outward: f64 },
}

/// 隣のセルへ移るときに通る面と、その向き（面の正の向きなら 1、逆なら -1）を掛けた重み。
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FacePath {
    faces: [(usize, f64); MAX_FACES],
    count: usize,
}

/// 1回の移動で通る面の数の上限。セルの角をちょうど通るときは、隣り合う2つのセルを経由する
/// 2通りの経路（面を2つずつ）に 1/2 ずつ分ける
const MAX_FACES: usize = 4;

impl FacePath {
    fn new(faces: &[(usize, f64)]) -> Self {
        let mut array = [(0, 0.0); MAX_FACES];
        array[..faces.len()].copy_from_slice(faces);
        Self {
            faces: array,
            count: faces.len(),
        }
    }

    pub fn faces(&self) -> &[(usize, f64)] {
        &self.faces[..self.count]
    }
}

impl CoarseMesh {
    pub fn new(tracking: &Tracking) -> Self {
        let lattice = &tracking.lattice;
        let cell_of_region = (0..lattice.regions())
            .map(|r| lattice.region_cell(r))
            .collect();
        let half = (0.5 * lattice.width(), 0.5 * lattice.height());
        let mut mesh = Self {
            widths: lattice.widths().to_vec(),
            heights: lattice.heights().to_vec(),
            cell_of_region,
            ends: Vec::with_capacity(tracking.tracks.len()),
        };
        for track in &tracking.tracks {
            let first = mesh.cell(track.segments[0].region);
            let last = mesh.cell(track.segments[track.segments.len() - 1].region);
            let start = mesh.outer_face(first, track.start, half);
            let end = mesh.outer_face(last, track.end, half);
            // 始点では入る（外向きの逆）、終点では出る
            mesh.ends.push([(start.0, -start.1), end]);
        }
        mesh
    }

    /// 列の数
    fn columns(&self) -> usize {
        self.widths.len()
    }

    /// 行の数
    fn rows(&self) -> usize {
        self.heights.len()
    }

    pub fn cells(&self) -> usize {
        self.columns() * self.rows()
    }

    pub fn faces(&self) -> usize {
        let (nc, nr) = (self.columns(), self.rows());
        nr * (nc + 1) + (nr + 1) * nc
    }

    /// 領域の属するセル
    pub fn cell(&self, region: usize) -> usize {
        self.cell_of_region[region]
    }

    /// セルの4つの面と、外向きが面の正の向きなら 1、逆なら -1（左、右、上、下の順）。
    pub fn cell_faces(&self, cell: usize) -> [(usize, f64); 4] {
        let (r, c) = (cell / self.columns(), cell % self.columns());
        [
            (self.vertical(r, c), -1.0),
            (self.vertical(r, c + 1), 1.0),
            (self.horizontal(r, c), -1.0),
            (self.horizontal(r + 1, c), 1.0),
        ]
    }

    /// [トラック] の始点で入る面と終点で出る面（順方向に通るときの向き付き）
    pub(crate) fn ends(&self, track: usize) -> [(usize, f64); 2] {
        self.ends[track]
    }

    /// 面の両側のセル。
    pub fn face_cells(&self, face: usize) -> FaceCells {
        let (nc, nr) = (self.columns(), self.rows());
        let vertical = nr * (nc + 1);
        // 面が格子線 `line`（0..=last）の上にあり、線に沿った位置で決まるセルの番号を `cell(k)`
        // （k は線の手前側のセルの行または列）で求める
        let side = |line: usize, last: usize, cell: &dyn Fn(usize) -> usize| {
            if line == 0 {
                FaceCells::Outer {
                    cell: cell(0),
                    outward: -1.0,
                }
            } else if line == last {
                FaceCells::Outer {
                    cell: cell(last - 1),
                    outward: 1.0,
                }
            } else {
                FaceCells::Inner {
                    low: cell(line - 1),
                    high: cell(line),
                }
            }
        };
        if face < vertical {
            let (row, line) = (face / (nc + 1), face % (nc + 1));
            side(line, nc, &|column| row * nc + column)
        } else {
            let (line, column) = ((face - vertical) / nc, (face - vertical) % nc);
            side(line, nr, &|row| row * nc + column)
        }
    }

    /// 縦の面（x 方向に通る）なら true
    pub fn is_vertical(&self, face: usize) -> bool {
        face < self.rows() * (self.columns() + 1)
    }

    /// 面の長さ (cm)
    pub fn face_length(&self, face: usize) -> f64 {
        let nc = self.columns();
        if self.is_vertical(face) {
            self.heights[face / (nc + 1)]
        } else {
            self.widths[(face - self.rows() * (nc + 1)) % nc]
        }
    }

    /// 内側の面の長さと、面の低い側と高い側のセルの、面に垂直な方向の幅
    fn face_geometry(&self, face: usize) -> (f64, f64, f64) {
        let (nc, nr) = (self.columns(), self.rows());
        let vertical = nr * (nc + 1);
        if face < vertical {
            let (row, line) = (face / (nc + 1), face % (nc + 1));
            (self.heights[row], self.widths[line - 1], self.widths[line])
        } else {
            let (line, column) = ((face - vertical) / nc, (face - vertical) % nc);
            (
                self.widths[column],
                self.heights[line - 1],
                self.heights[line],
            )
        }
    }

    fn vertical(&self, row: usize, line: usize) -> usize {
        vertical_face(self.columns(), row, line)
    }

    fn horizontal(&self, line: usize, column: usize) -> usize {
        horizontal_face(self.rows(), self.columns(), line, column)
    }

    /// セル `from` から隣のセル `to` へ移るときに通る面と向き。
    pub(crate) fn path(&self, from: usize, to: usize) -> Result<FacePath, CmfdError> {
        let nc = self.columns();
        let (r1, c1) = (from / nc, from % nc);
        let (r2, c2) = (to / nc, to % nc);
        let sign = |a: usize, b: usize| if b > a { 1.0 } else { -1.0 };
        match (r1.abs_diff(r2), c1.abs_diff(c2)) {
            (0, 1) => Ok(FacePath::new(&[(
                self.vertical(r1, c1.max(c2)),
                sign(c1, c2),
            )])),
            (1, 0) => Ok(FacePath::new(&[(
                self.horizontal(r1.max(r2), c1),
                sign(r1, r2),
            )])),
            // 角を通る: (r1, c2) を経由する経路と (r2, c1) を経由する経路に半分ずつ
            (1, 1) => Ok(FacePath::new(&[
                (self.vertical(r1, c1.max(c2)), 0.5 * sign(c1, c2)),
                (self.horizontal(r1.max(r2), c2), 0.5 * sign(r1, r2)),
                (self.horizontal(r1.max(r2), c1), 0.5 * sign(r1, r2)),
                (self.vertical(r2, c1.max(c2)), 0.5 * sign(c1, c2)),
            ])),
            _ => Err(CmfdError::NonAdjacent { from, to }),
        }
    }

    /// トラックの端の点 `point`（外周の上）にある、セル `cell` の外周の面と、
    /// 外向きが面の正の向きなら 1、逆なら -1。
    fn outer_face(&self, cell: usize, point: (f64, f64), half: (f64, f64)) -> (usize, f64) {
        let (nc, nr) = (self.columns(), self.rows());
        let (r, c) = (cell / nc, cell % nc);
        // 外向きの向き。トラックは外周の角を通らないので、点に最も近い辺が点のある辺
        let sides = [
            ((point.0 + half.0).abs(), self.vertical(r, 0), -1.0),
            ((point.0 - half.0).abs(), self.vertical(r, nc), 1.0),
            ((point.1 - half.1).abs(), self.horizontal(0, c), -1.0),
            ((point.1 + half.1).abs(), self.horizontal(nr, c), 1.0),
        ];
        let [first, rest @ ..] = sides;
        let (_, face, outward) =
            rest.into_iter().fold(
                first,
                |nearest, side| {
                    if side.0 < nearest.0 { side } else { nearest }
                },
            );
        (face, outward)
    }
}

/// `columns` 列の格子の、セルの行 `row` の左から `line` 番目の格子線（0..=列の数）の上の縦の面の番号
/// （`CoarseMesh` の番号）
pub fn vertical_face(columns: usize, row: usize, line: usize) -> usize {
    row * (columns + 1) + line
}

/// `rows` 行 `columns` 列の格子の、セルの列 `column` の上から `line` 番目の格子線（0..=行の数）の上の
/// 横の面の番号（`CoarseMesh` の番号。縦の面のあと）
pub fn horizontal_face(rows: usize, columns: usize, line: usize, column: usize) -> usize {
    rows * (columns + 1) + line * columns + column
}

/// CMFD の固有値問題の打ち切り条件。
///
/// MOC の反復の判定（`moc::Options`）より十分小さくする。CMFD の打ち切りの誤差が、MOC の
/// 反復ごとの変化に紛れて収束の判定を誤らせないため。MOC の反復では、スイープによる変化が MOC の判定の
/// 何倍あるかに比例して緩めて解く（ADR-0110）。ここの値は、スイープによる変化が MOC の判定以下になったときの精度で、
/// 収束と判定した反復でも、スイープによる変化が判定より大きければ、その比だけ緩めたままになる。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CmfdOptions {
    /// 実効増倍率の変化の上限
    pub k_tolerance: f64,
    /// 中性子束の相対変化の上限
    pub flux_tolerance: f64,
    pub max_iterations: usize,
}

impl Default for CmfdOptions {
    fn default() -> Self {
        Self {
            k_tolerance: 1.0e-10,
            flux_tolerance: 1.0e-8,
            max_iterations: 100_000,
        }
    }
}

/// CMFD の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum CmfdError {
    #[error("the MOC flux averaged over cell {cell} in group {group} is not positive ({flux})")]
    NonPositiveFlux {
        cell: usize,
        group: usize,
        flux: f64,
    },
    #[error("the coarse-mesh flux of cell {cell} in group {group} became {flux}")]
    NonPositiveSolution {
        cell: usize,
        group: usize,
        flux: f64,
    },
    #[error("a track moves from cell {from} to cell {to}, which do not touch")]
    NonAdjacent { from: usize, to: usize },
    #[error("the coarse-mesh matrix of cell {0} is singular")]
    Singular(usize),
    #[error("the coarse-mesh eigenvalue problem did not converge in {0} iterations")]
    NotConverged(usize),
}

/// 粗メッシュのセルの、同じ核分裂スペクトルの項の和。核分裂源は χ (νΣf · φ)
#[derive(Debug, Clone)]
struct Spectrum {
    chi: Vec<f64>,
    /// χ の和（核分裂源の生成率に掛ける）
    chi_sum: f64,
    /// [群] の体積 × νΣf（均質化した中性子束あたり）
    nu_fission: Vec<f64>,
}

impl Spectrum {
    /// 中性子束 `phi`（[群]）での νΣf · φ
    fn rate(&self, phi: &[f64]) -> f64 {
        self.nu_fission.iter().zip(phi).map(|(a, b)| a * b).sum()
    }
}

/// セルごとに均質化した粗メッシュの問題（流れに合わせて補正した拡散方程式）。
#[derive(Debug, Clone)]
pub(crate) struct Coarse {
    groups: usize,
    /// [セル][群] の均質化した中性子束
    flux: Vec<f64>,
    /// [セル] の核分裂源 V χ νΣf を、核分裂のスペクトルごとの項の和として持つ（密な群×群の行列にすると、
    /// 右辺の計算が群の数の2乗になる）
    fission: Vec<Vec<Spectrum>>,
    /// [セル] の対角のブロック（群×群）の LU 分解
    blocks: Vec<Lu>,
    /// [セル] の隣のセルとの結合: (隣のセル, [群] の係数)。方程式は
    /// A_c φ_c + Σ 係数 φ_隣 = V F φ_c / k
    neighbors: Vec<Vec<(usize, Vec<f64>)>>,
}

impl Coarse {
    /// 細かい領域の [領域][群] の中性子束 `flux` と体積から、セルごとに断面積を均質化し、
    /// スイープで数えた [面][群] の正味の流れ `current` に合うよう結合係数を補正する。
    pub fn homogenize(
        mesh: &CoarseMesh,
        materials: &[&Material],
        volumes: &[f64],
        flux: &[Vec<f64>],
        current: &[f64],
    ) -> Result<Self, CmfdError> {
        let ng = materials[0].groups();
        let cells = mesh.cells();
        let mut volume = vec![0.0; cells];
        let mut phi = vec![0.0; cells * ng];
        // 反応率の和（後で中性子束で割る）
        let mut total = vec![0.0; cells * ng];
        let mut scatter = vec![0.0; cells * ng * ng];
        let mut fission: Vec<Vec<Spectrum>> = vec![Vec::new(); cells];
        for (region, m) in materials.iter().enumerate() {
            let (c, v, f) = (mesh.cell(region), volumes[region], &flux[region]);
            volume[c] += v;
            for g in 0..ng {
                phi[c * ng + g] += v * f[g];
                total[c * ng + g] += v * m.total[g] * f[g];
                for to in 0..ng {
                    scatter[(c * ng + to) * ng + g] += v * m.scatter[g][to] * f[g];
                }
            }
            for term in &m.fission_terms {
                // 同じスペクトルの項はまとめる
                let terms = &mut fission[c];
                let index = match terms.iter().position(|s| s.chi == term.chi) {
                    Some(i) => i,
                    None => {
                        terms.push(Spectrum {
                            chi: term.chi.clone(),
                            chi_sum: term.chi.iter().sum(),
                            nu_fission: vec![0.0; ng],
                        });
                        terms.len() - 1
                    }
                };
                for (r, (n, phi)) in terms[index]
                    .nu_fission
                    .iter_mut()
                    .zip(term.nu_fission.iter().zip(f))
                {
                    *r += v * n * phi;
                }
            }
        }
        for c in 0..cells {
            for g in 0..ng {
                let reaction = phi[c * ng + g];
                if reaction <= 0.0 || !reaction.is_finite() {
                    return Err(CmfdError::NonPositiveFlux {
                        cell: c,
                        group: g,
                        flux: reaction,
                    });
                }
                phi[c * ng + g] = reaction / volume[c];
            }
        }
        // 断面積 = 反応率 / (体積 × 中性子束)。行列は φ に掛けて反応率になる形で持つので、
        // 列（入射の群）ごとに中性子束で割る
        let sigma_t = |c: usize, g: usize| total[c * ng + g] / (volume[c] * phi[c * ng + g]);
        // 光学的に厚いセルでは従来の CMFD が不安定になるので、拡散係数に θ(Σt h) h を足す
        // （optimally diffusive CMFD）。h は面に垂直な方向のセルの幅
        let diffusion = |c: usize, g: usize, h: f64| {
            let sigma = sigma_t(c, g);
            1.0 / (3.0 * sigma) + artificial_diffusion(sigma * h) * h
        };
        for c in 0..cells {
            for to in 0..ng {
                for g in 0..ng {
                    let j = (c * ng + to) * ng + g;
                    scatter[j] /= phi[c * ng + g];
                }
            }
        }

        for (c, terms) in fission.iter_mut().enumerate() {
            for s in terms.iter_mut() {
                for (g, n) in s.nu_fission.iter_mut().enumerate() {
                    *n /= phi[c * ng + g];
                }
            }
        }

        // 対角のブロックの対角項: 全断面積 + 面からの漏れ。隣との結合
        let mut diagonal: Vec<f64> = (0..cells * ng).map(|j| total[j] / phi[j]).collect();
        let mut neighbors: Vec<Vec<(usize, Vec<f64>)>> = vec![Vec::new(); cells];
        for face in 0..mesh.faces() {
            let net = &current[face * ng..][..ng];
            match mesh.face_cells(face) {
                FaceCells::Inner { low, high } => {
                    // 面の正の向きの流れ J = -D̃ (φ_h - φ_l) - D̂ (φ_h + φ_l)（面の長さあたり）。
                    // D̃ = 2 D_l D_h / (w_l D_h + w_h D_l)（w は面に垂直な方向のセルの幅）
                    let (length, wl, wh) = mesh.face_geometry(face);
                    let mut on_high = vec![0.0; ng];
                    let mut on_low = vec![0.0; ng];
                    for g in 0..ng {
                        let (dl, dh) = (diffusion(low, g, wl), diffusion(high, g, wh));
                        let coupling = 2.0 * dl * dh / (wl * dh + wh * dl);
                        let (fl, fh) = (phi[low * ng + g], phi[high * ng + g]);
                        let correction = -(net[g] / length + coupling * (fh - fl)) / (fh + fl);
                        // 低い側の式では J が出ていく流れ、高い側の式では入ってくる流れ
                        diagonal[low * ng + g] += length * (coupling - correction);
                        on_high[g] = length * (-coupling - correction);
                        diagonal[high * ng + g] += length * (coupling + correction);
                        on_low[g] = length * (-coupling + correction);
                    }
                    neighbors[low].push((high, on_high));
                    neighbors[high].push((low, on_low));
                }
                // 外周の面: 外向きの流れを、すべて補正の項 D̂ φ で表す
                FaceCells::Outer { cell, outward } => {
                    for g in 0..ng {
                        diagonal[cell * ng + g] += outward * net[g] / phi[cell * ng + g];
                    }
                }
            }
        }
        let blocks = (0..cells)
            .map(|c| {
                let mut a: Vec<f64> = scatter[c * ng * ng..][..ng * ng]
                    .iter()
                    .map(|s| -s)
                    .collect();
                for g in 0..ng {
                    a[g * ng + g] += diagonal[c * ng + g];
                }
                Lu::new(a, ng).ok_or(CmfdError::Singular(c))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            groups: ng,
            flux: phi,
            fission,
            blocks,
            neighbors,
        })
    }

    /// 固有値問題を `k` から解き始め、k と [セル][群] の中性子束の比（解 / 均質化した中性子束、
    /// どちらも生成率 1 に規格化）を返す。
    ///
    /// べき乗法の各反復で、ブロック Gauss-Seidel をセルの順に1回だけ行う。収束すると
    /// Gauss-Seidel でも解が変わらないので、内側の反復を収束させなくても固有値問題の解になる。
    pub fn solve(&self, k: f64, options: &CmfdOptions) -> Result<(f64, Vec<f64>), CmfdError> {
        let ng = self.groups;
        let cells = self.blocks.len();
        let production = |phi: &[f64]| -> f64 {
            (0..cells)
                .flat_map(|c| {
                    let p = &phi[c * ng..][..ng];
                    self.fission[c].iter().map(move |s| s.chi_sum * s.rate(p))
                })
                .sum()
        };
        let start = production(&self.flux);
        let mut phi: Vec<f64> = self.flux.iter().map(|f| f / start).collect();
        let mut k = k;
        let mut k_history = Remaining::default();
        let mut flux_history = Remaining::default();
        let mut rhs = vec![0.0; ng];
        for iteration in 1..=options.max_iterations {
            let previous = phi.clone();
            for c in 0..cells {
                rhs.fill(0.0);
                for s in &self.fission[c] {
                    let rate = s.rate(&previous[c * ng..][..ng]) / k;
                    for (r, x) in rhs.iter_mut().zip(&s.chi) {
                        *r += x * rate;
                    }
                }
                for (other, coefficients) in &self.neighbors[c] {
                    for g in 0..ng {
                        rhs[g] -= coefficients[g] * phi[other * ng + g];
                    }
                }
                self.blocks[c].solve(&mut rhs);
                phi[c * ng..][..ng].copy_from_slice(&rhs);
            }
            let rate = production(&phi);
            let new_k = k * rate;
            for f in &mut phi {
                *f /= rate;
            }
            if let Some(j) = phi.iter().position(|&f| f <= 0.0 || !f.is_finite()) {
                return Err(CmfdError::NonPositiveSolution {
                    cell: j / ng,
                    group: j % ng,
                    flux: phi[j],
                });
            }
            let change = phi
                .iter()
                .zip(&previous)
                .map(|(a, b)| (a - b).abs() / a.abs().max(b.abs()))
                .fold(0.0, f64::max);
            let k_remaining = k_history.update((new_k - k).abs());
            let flux_remaining = flux_history.update(change);
            k = new_k;
            if k_remaining < options.k_tolerance && flux_remaining < options.flux_tolerance {
                let ratio = phi
                    .iter()
                    .zip(&self.flux)
                    .map(|(new, old)| new / (old / start))
                    .collect();
                tracing::debug!(iterations = iteration, k, "coarse-mesh eigenvalue");
                return Ok((k, ratio));
            }
        }
        Err(CmfdError::NotConverged(options.max_iterations))
    }
}

/// optimally diffusive CMFD の θ（拡散係数に足す θ h の係数）を、粗メッシュのセルの光学的厚さ
/// τ = Σt h から求める。
///
/// A. Zhu et al., "An optimally diffusive Coarse Mesh Finite Difference method to accelerate
/// neutron transport calculations", Ann. Nucl. Energy 95 (2016) 116–124 の式 (23) と Table 1。
/// フーリエ解析で求めた最適値の当てはめで、τ < 1 では 0（従来の CMFD）、τ ≥ 14 では飽和値。
/// C5G7 の反射体（熱群で τ ≈ 4）では従来の CMFD が振動して収束せず、partial current CMFD に
/// 当たる θ = 1/4 より、細かい領域の多い場合に反復が少ない（ADR-0062 を参照）
pub(crate) fn artificial_diffusion(tau: f64) -> f64 {
    /// 多項式を使う範囲の下限と上限
    const THIN: f64 = 1.0;
    const THICK: f64 = 14.0;
    /// τ ≥ 14 での値
    const SATURATED: f64 = 0.127;
    /// θ = Σ a_i τ^i の係数 a_0..a_6
    const COEFFICIENTS: [f64; 7] = [
        -5.542780e-02,
        8.740501e-02,
        -2.152599e-02,
        3.145553e-03,
        -2.683648e-04,
        1.222516e-05,
        -2.284879e-07,
    ];
    if tau < THIN {
        0.0
    } else if tau >= THICK {
        SATURATED
    } else {
        COEFFICIENTS.iter().rev().fold(0.0, |sum, a| sum * tau + a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};
    use crate::material::{FissionTerm, Material};

    fn two_by_two() -> Tracking {
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
        Tracking::new(&lattice, 16, 0.05).unwrap()
    }

    #[test]
    fn face_numbering() {
        let mesh = CoarseMesh::new(&two_by_two());
        assert_eq!((mesh.cells(), mesh.faces()), (4, 12));
        // 左上のセルの右の面は、右上のセルの左の面
        let (right, out) = mesh.cell_faces(0)[1];
        let (left, into) = mesh.cell_faces(1)[0];
        assert_eq!((right, out, into), (left, 1.0, -1.0));
        // 左上のセルの下の面は、左下のセルの上の面
        assert_eq!(mesh.cell_faces(0)[3].0, mesh.cell_faces(2)[2].0);
        // 全セルの面を数えると、内側の面は2回、外周の面は1回
        let mut count = vec![0; mesh.faces()];
        for cell in 0..mesh.cells() {
            for (face, _) in mesh.cell_faces(cell) {
                count[face] += 1;
            }
        }
        assert_eq!(count.iter().filter(|&&c| c == 2).count(), 4);
        assert_eq!(count.iter().filter(|&&c| c == 1).count(), 8);
    }

    #[test]
    fn homogeneous_cell_gives_the_infinite_medium_k() {
        // 4辺が鏡面反射の均質なセルでは面の流れが 0 で、粗メッシュの k は無限媒質の k
        // k∞ = νΣfᵀ (diag Σt - Σsᵀ)⁻¹ χ（χ = (1, 0)）
        let m = Material {
            total: vec![0.5, 1.2],
            scatter: vec![vec![0.45, 0.02], vec![0.0, 0.9]],
            anisotropic: Vec::new(),
            fission: vec![0.008 / 2.5, 0.25 / 2.5],
            nu_fission: vec![0.008, 0.25],
            fission_terms: vec![FissionTerm {
                chi: vec![1.0, 0.0],
                nu_fission: vec![0.008, 0.25],
            }],
        };
        let a = [[0.5 - 0.45, 0.0], [-0.02, 1.2 - 0.9]];
        let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        let expected = (0.008 * a[1][1] - 0.25 * a[1][0]) / det;
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
        let mesh = CoarseMesh::new(&tracking);
        let volumes = tracking.volumes(2);
        // 均質なので、領域の中性子束の形によらない
        let flux = vec![vec![1.0, 0.3], vec![0.8, 0.5]];
        let current = vec![0.0; mesh.faces() * 2];
        let coarse = Coarse::homogenize(&mesh, &[&m, &m], &volumes, &flux, &current).unwrap();
        let (k, ratio) = coarse.solve(1.0, &CmfdOptions::default()).unwrap();
        assert!((k - expected).abs() < 1e-9, "{k} vs {expected}");
        assert_eq!(ratio.len(), 2);
    }

    #[test]
    fn artificial_diffusion_follows_the_optically_thick_fit() {
        // Zhu ら (2016): 光学的に薄いと 0（従来の CMFD）、厚いと 0.127 に飽和し、Σt h = 4 で約 0.1
        assert_eq!(artificial_diffusion(0.5), 0.0);
        assert!((artificial_diffusion(4.0) - 0.1).abs() < 0.01);
        assert_eq!(artificial_diffusion(20.0), 0.127);
        // 多項式は飽和値に滑らかにつながる
        assert!((artificial_diffusion(14.0 - 1e-9) - 0.127).abs() < 0.002);
        let samples: Vec<f64> = (2..28)
            .map(|i| artificial_diffusion(0.5 * i as f64))
            .collect();
        assert!(samples.windows(2).all(|w| w[1] >= w[0]), "{samples:?}");
    }

    /// 面の外向きの重みの和（セル `cell` から出る分が正）
    fn outflow(mesh: &CoarseMesh, cell: usize, faces: &[(usize, f64)]) -> f64 {
        mesh.cell_faces(cell)
            .iter()
            .map(|&(face, out)| {
                faces
                    .iter()
                    .filter(|&&(f, _)| f == face)
                    .map(|&(_, w)| out * w)
                    .sum::<f64>()
            })
            .sum()
    }

    #[test]
    fn paths_leave_one_cell_and_enter_the_other() {
        let mesh = CoarseMesh::new(&two_by_two());
        for from in 0..4 {
            for to in (0..4).filter(|&to| to != from) {
                let path = mesh.path(from, to).unwrap();
                let faces = path.faces();
                assert_eq!(outflow(&mesh, from, faces), 1.0);
                assert_eq!(outflow(&mesh, to, faces), -1.0);
                // 角を通るときの経由するセルは、入って出るので正味 0
                for other in (0..4).filter(|&c| c != from && c != to) {
                    assert_eq!(outflow(&mesh, other, faces), 0.0);
                }
            }
        }
    }

    #[test]
    fn cells_that_do_not_touch_have_no_path() {
        let water = PinCell::new(1.26, vec![]).unwrap();
        let lattice = Lattice::new(vec![water], vec![vec![0; 3]; 3]).unwrap();
        let mesh = CoarseMesh::new(&Tracking::new(&lattice, 8, 0.2).unwrap());
        assert!(matches!(
            mesh.path(0, 2),
            Err(CmfdError::NonAdjacent { from: 0, to: 2 })
        ));
    }

    #[test]
    fn faces_know_their_cells() {
        let mesh = CoarseMesh::new(&two_by_two());
        // 左上のセルの左の面は外周（外向きは面の負の向き）、右の面は右上のセルとの間
        let [(left, _), (right, _), (top, _), (bottom, _)] = mesh.cell_faces(0);
        assert_eq!(
            mesh.face_cells(left),
            FaceCells::Outer {
                cell: 0,
                outward: -1.0
            }
        );
        assert_eq!(mesh.face_cells(right), FaceCells::Inner { low: 0, high: 1 });
        assert_eq!(
            mesh.face_cells(top),
            FaceCells::Outer {
                cell: 0,
                outward: -1.0
            }
        );
        assert_eq!(
            mesh.face_cells(bottom),
            FaceCells::Inner { low: 0, high: 2 }
        );
    }

    #[test]
    fn track_ends_are_on_the_outer_faces_of_their_cells() {
        let tracking = two_by_two();
        let mesh = CoarseMesh::new(&tracking);
        for (index, track) in tracking.tracks.iter().enumerate() {
            let [start, end] = mesh.ends(index);
            let first = mesh.cell(track.segments[0].region);
            let last = mesh.cell(track.segments[track.segments.len() - 1].region);
            assert_eq!(outflow(&mesh, first, &[start]), -1.0, "track {index}");
            assert_eq!(outflow(&mesh, last, &[end]), 1.0, "track {index}");
        }
    }
}
