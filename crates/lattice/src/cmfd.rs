//! CMFD（粗メッシュ有限差分）による MOC の加速。
//!
//! 粗メッシュは格子のセルとする（格子の周りの余白も、薄い行と列のセルになる）。
//! スイープでセルの面を通る正味の中性子の流れを数え、
//! 拡散の結合係数を流れに合うよう補正した粗メッシュの固有値問題を解いて、
//! 細かい領域の中性子束をセルごとの比で更新する。

use crate::material::Material;
use crate::moc::Remaining;
use crate::symmetry::Orbits;
use crate::tracking::Tracking;
use rayon::prelude::*;
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
    /// [トラック] の始点で入る経路と終点で出る経路（順方向に通るときの向き付き）。ふつうは端のセルの外周の面 1 つ
    /// （重み ±1）。端が外周の 2 つの面の境にあれば、自分の外周の面、隣のセルとの間の面、隣の外周の面の 3 つ（重み ±1/2）
    ends: Vec<[FacePath; 2]>,
    /// [トラック] の始点と終点の角度中性子束が属するセルと重み（`ends` の経路と同じ分け方）
    end_cells: Vec<[EndCells; 2]>,
}

/// トラックの端の角度中性子束が属するセルと重み（1 つで重み 1 か、2 つで 1/2 ずつ）
#[derive(Debug, Clone, Copy, PartialEq)]
struct EndCells {
    cells: [(usize, f64); 2],
    count: usize,
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

    /// 重みを `factor` 倍した経路（逆向きに通るときは -1 倍）
    pub fn scaled(&self, factor: f64) -> Self {
        let mut path = *self;
        for (_, weight) in &mut path.faces[..path.count] {
            *weight *= factor;
        }
        path
    }
}

impl CoarseMesh {
    pub fn new(tracking: &Tracking) -> Self {
        let lattice = &tracking.lattice;
        let cell_of_region = (0..lattice.regions())
            .map(|r| lattice.region_cell(r))
            .collect();
        let half = (0.5 * lattice.width(), 0.5 * lattice.height());
        let coincident = lattice.coincident();
        let lines = (lattice.column_lines(), lattice.row_lines());
        let mut mesh = Self {
            widths: lattice.widths().to_vec(),
            heights: lattice.heights().to_vec(),
            cell_of_region,
            ends: Vec::with_capacity(tracking.tracks.len()),
            end_cells: Vec::with_capacity(tracking.tracks.len()),
        };
        for track in &tracking.tracks {
            let first = mesh.cell(track.segments[0].region);
            let last = mesh.cell(track.segments[track.segments.len() - 1].region);
            let (start, start_cells) =
                mesh.end_path(first, track.start, half, lines, coincident, false);
            let (end, end_cells) = mesh.end_path(last, track.end, half, lines, coincident, true);
            mesh.ends.push([start, end]);
            mesh.end_cells.push([start_cells, end_cells]);
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

    /// トラック `track` の始点で入る経路と終点で出る経路（順方向に通るときの向き付き。`ends` のフィールドを参照）
    pub(crate) fn ends(&self, track: usize) -> [FacePath; 2] {
        self.ends[track]
    }

    /// トラック `track` の終点（`exit` が真）か始点の角度中性子束が属するセルと重み
    pub(crate) fn end_cells(&self, track: usize, exit: bool) -> &[(usize, f64)] {
        let end = &self.end_cells[track][usize::from(exit)];
        &end.cells[..end.count]
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

    /// トラックの端の点 `point`（外周の上）で、セル `cell` から出る（`exiting`）か、セルに入る経路。端が外周の 2 つの面の境
    /// （格子線と外周の交わる点。同じ点とみなす距離 `coincident` 以内）にあるときは、トラックが代表する帯の半分ずつが
    /// 境の両側の面を通るとして、2 通りの経路に 1/2 ずつ分ける: 自分のセルの外周の面と、隣のセルを経由する経路（セルの間の
    /// 面と隣のセルの外周の面）。片側の面だけに数えると、出ていくトラックと、反射してつながる入るトラックが境の別の側に
    /// 数えられ、鏡面反射の面に正味の中性子流が残る
    fn end_path(
        &self,
        cell: usize,
        point: (f64, f64),
        half: (f64, f64),
        (xs, ys): (&[f64], &[f64]),
        coincident: f64,
        exiting: bool,
    ) -> (FacePath, EndCells) {
        let (nc, nr) = (self.columns(), self.rows());
        let (r, c) = (cell / nc, cell % nc);
        let faces = self.cell_faces(cell);
        // 点に最も近い辺が点のある辺（左、右、上、下。トラックは外周の角を通らない）
        let distances = [
            (point.0 + half.0).abs(),
            (point.0 - half.0).abs(),
            (point.1 - half.1).abs(),
            (point.1 + half.1).abs(),
        ];
        let side = (1..4).fold(0, |nearest, s| {
            if distances[s] < distances[nearest] {
                s
            } else {
                nearest
            }
        });
        let own = faces[side];
        // 辺に沿った格子線のうち、点がセルの両端の線にあれば、その向こうのセルと、そちらへ移る面（cell_faces の番号）
        let neighbour = if side >= 2 {
            if c > 0 && (point.0 - xs[c]).abs() <= coincident {
                Some((cell - 1, 0))
            } else if c + 1 < nc && (point.0 - xs[c + 1]).abs() <= coincident {
                Some((cell + 1, 1))
            } else {
                None
            }
        } else if r > 0 && (point.1 - ys[r]).abs() <= coincident {
            Some((cell - nc, 2))
        } else if r + 1 < nr && (point.1 - ys[r + 1]).abs() <= coincident {
            Some((cell + nc, 3))
        } else {
            None
        };
        let sign = if exiting { 1.0 } else { -1.0 };
        match neighbour {
            None => (
                FacePath::new(&[(own.0, sign * own.1)]),
                EndCells {
                    cells: [(cell, 1.0), (cell, 0.0)],
                    count: 1,
                },
            ),
            Some((other, toward)) => {
                let theirs = self.cell_faces(other)[side];
                // 自分のセルから隣へ移る向きの重み。入るときは隣から自分へ移るので逆
                let between = faces[toward];
                (
                    FacePath::new(&[
                        (own.0, 0.5 * sign * own.1),
                        (between.0, 0.5 * sign * between.1),
                        (theirs.0, 0.5 * sign * theirs.1),
                    ]),
                    EndCells {
                        cells: [(cell, 0.5), (other, 0.5)],
                        count: 2,
                    },
                )
            }
        }
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
/// 収束と判定した反復でも、スイープによる変化が判定より大きければ、その比だけ緩めたままになる。正の有限で、MOC の
/// 判定以下であること（`moc::solve` が確かめる）。
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
    #[error(
        "the coarse-mesh block of cell {0} after eliminating the cells before it in the line is singular"
    )]
    SingularElimination(usize),
    #[error(
        "cell {cell} of a one-row or one-column lattice couples to cell {other}, which is not next to it"
    )]
    Unchained { cell: usize, other: usize },
    /// シフトした逆反復の分解が特異か、最初の解が正でない（`Coarse::solve_chain` がシフトなしで解き直す）
    #[error("the Wielandt shift of the coarse-mesh eigenvalue problem is not above the eigenvalue")]
    ShiftBelowEigenvalue,
    #[error("the coarse-mesh eigenvalue problem did not converge in {0} iterations")]
    NotConverged(usize),
}

/// 1 行か 1 列の格子の粗メッシュの逆反復のシフト k_s = k (1 + WIELANDT) の幅（ADR-0158）。収束の比は
/// (1/k₀ − 1/k_s) / (1/k₁ − 1/k_s)（k₀ > k₁ は 1 番目と 2 番目の固有値）。k が k₀ / (1 + WIELANDT) 以下なら
/// k_s ≤ k₀ になり、シフトなしで解き直す（`Coarse::solve_chain`）
const WIELANDT: f64 = 0.01;

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
    /// [代表] の解くセル（エラーの表示用）
    cells: Vec<usize>,
    /// [代表] の組のセルの数（生成率の和の重み）
    sizes: Vec<f64>,
    /// [セル][群] の均質化した中性子束
    flux: Vec<f64>,
    /// [セル] の核分裂源 V χ νΣf を、核分裂のスペクトルごとの項の和として持つ（密な群×群の行列にすると、
    /// 右辺の計算が群の数の2乗になる）
    fission: Vec<Vec<Spectrum>>,
    solver: Solver,
}

/// 粗メッシュの行列と、その解き方。方程式は A_c φ_c + Σ 係数 φ_隣 = V F φ_c / k
#[derive(Debug, Clone)]
enum Solver {
    /// べき乗法の反復ごとにブロック Gauss-Seidel を 1 回
    GaussSeidel {
        /// [セル] の対角のブロック（群×群）の LU 分解
        blocks: Vec<Lu>,
        /// [セル] の隣のセルとの結合: (隣のセル, [群] の係数)
        neighbors: Vec<Vec<(usize, Vec<f64>)>>,
    },
    /// 1 行か 1 列の格子で、ブロック三重対角の直接解法の逆反復
    Chain(Chain),
}

/// 代表のセル c の隣が c − 1 と c + 1 だけの結合（1 行か 1 列の格子）。行列はブロック三重対角になる。
#[derive(Debug, Clone)]
struct Chain {
    /// 群の数
    groups: usize,
    /// [代表] の対角のブロック A_c（群×群、行優先）
    matrices: Vec<Vec<f64>>,
    /// [代表][群] の c − 1 との結合の係数（c = 0 では 0）
    lower: Vec<f64>,
    /// [代表][群] の c + 1 との結合の係数（最後の代表では 0）
    upper: Vec<f64>,
    /// 対角の外の成分（結合とブロックの対角の外）がどれも 0 以下（Z 行列）で、隣との結合がどの群でも負（空間で切れて
    /// いない。群の間が散乱と核分裂でつながることは確かめない）。
    /// シフトの条件（ADR-0158）
    z_matrix: bool,
}

impl Chain {
    /// [代表] の対角のブロック `matrices` と隣との結合 `neighbors`（`Coarse` と同じ形）から作る
    fn new(
        groups: usize,
        cells_of: &[usize],
        matrices: Vec<Vec<f64>>,
        neighbors: &[Vec<(usize, Vec<f64>)>],
    ) -> Result<Self, CmfdError> {
        let cells = matrices.len();
        let mut lower = vec![0.0; cells * groups];
        let mut upper = vec![0.0; cells * groups];
        // 代表は行か列の順に並び、対称線をはさむ隣は自分の像として対角に入っている
        for (c, list) in neighbors.iter().enumerate() {
            for (o, coefficients) in list {
                let side = match *o {
                    o if o + 1 == c => &mut lower,
                    o if o == c + 1 => &mut upper,
                    o => {
                        return Err(CmfdError::Unchained {
                            cell: cells_of[c],
                            other: cells_of[o],
                        });
                    }
                };
                for (s, a) in side[c * groups..][..groups].iter_mut().zip(coefficients) {
                    *s += a;
                }
            }
        }
        // 隣との結合がどの群でも負（0 ならチェーンが空間で切れて既約でなくなり、最初の解の符号で判定できない）
        let z_matrix = lower[groups..].iter().all(|&a| a < 0.0)
            && upper[..(cells - 1) * groups].iter().all(|&a| a < 0.0)
            && matrices.iter().all(|a| {
                a.iter()
                    .enumerate()
                    .all(|(index, &v)| index / groups == index % groups || v <= 0.0)
            });
        Ok(Self {
            groups,
            matrices,
            lower,
            upper,
            z_matrix,
        })
    }

    /// A − μ F（F は [代表] の核分裂源 χ νΣf）をブロック LU に分解する。対角のブロック S_c = B_c − L_c X_(c−1)
    /// （B_c = A_c − μ F_c、L_c は c − 1 との結合の対角行列）の LU 分解と、X_c = S_c⁻¹ U_c（U_c は c + 1 との結合の
    /// 対角行列）を持つ
    fn factor(
        &self,
        fission: &[Vec<Spectrum>],
        shift: f64,
        cells: &[usize],
    ) -> Result<ChainFactors, CmfdError> {
        let ng = self.groups;
        let n = self.matrices.len();
        let mut pivots = Vec::with_capacity(n);
        let mut eliminations: Vec<Vec<f64>> = Vec::with_capacity(n - 1);
        let mut column = vec![0.0; ng];
        for (c, (matrix, spectra)) in self.matrices.iter().zip(fission).enumerate() {
            let mut s = matrix.clone();
            if shift != 0.0 {
                for spectrum in spectra {
                    for (i, x) in spectrum.chi.iter().enumerate() {
                        for (j, f) in spectrum.nu_fission.iter().enumerate() {
                            s[i * ng + j] -= shift * x * f;
                        }
                    }
                }
            }
            if let Some(previous) = eliminations.last() {
                let lower = &self.lower[c * ng..][..ng];
                for (i, l) in lower.iter().enumerate() {
                    for j in 0..ng {
                        s[i * ng + j] -= l * previous[i * ng + j];
                    }
                }
            }
            let lu = Lu::new(s, ng).ok_or(CmfdError::SingularElimination(cells[c]))?;
            if c + 1 < n {
                let mut x = vec![0.0; ng * ng];
                for (j, u) in self.upper[c * ng..][..ng].iter().enumerate() {
                    column.fill(0.0);
                    column[j] = *u;
                    lu.solve(&mut column);
                    for (i, v) in column.iter().enumerate() {
                        x[i * ng + j] = *v;
                    }
                }
                eliminations.push(x);
            }
            pivots.push(lu);
        }
        Ok(ChainFactors {
            pivots,
            eliminations,
        })
    }

    /// `b`（[代表][群]）を (A − μ F) φ = b の解 φ で置き換える
    fn substitute(&self, factors: &ChainFactors, b: &mut [f64]) {
        let ng = self.groups;
        let n = factors.pivots.len();
        for c in 0..n {
            if c > 0 {
                let (done, rest) = b.split_at_mut(c * ng);
                let previous = &done[(c - 1) * ng..];
                for ((r, l), y) in rest[..ng]
                    .iter_mut()
                    .zip(&self.lower[c * ng..][..ng])
                    .zip(previous)
                {
                    *r -= l * y;
                }
            }
            factors.pivots[c].solve(&mut b[c * ng..][..ng]);
        }
        for c in (0..n - 1).rev() {
            let (head, tail) = b.split_at_mut((c + 1) * ng);
            let next = &tail[..ng];
            let x = &factors.eliminations[c];
            for (i, r) in head[c * ng..].iter_mut().enumerate() {
                let sum: f64 = x[i * ng..][..ng].iter().zip(next).map(|(a, y)| a * y).sum();
                *r -= sum;
            }
        }
    }
}

/// `Chain::factor` の分解
struct ChainFactors {
    /// [代表] の S_c の LU 分解
    pivots: Vec<Lu>,
    /// [代表（最後を除く）] の X_c = S_c⁻¹ U_c（群×群、行優先）
    eliminations: Vec<Vec<f64>>,
}

impl Coarse {
    /// 細かい領域の中性子束と体積から、代表のセル（`orbits`）ごとに断面積を均質化し、スイープで数えた [面][群] の正味の
    /// 流れ `current` に合うよう結合係数を補正する。領域 r の中性子束は `flux[state[r]]`（[群]）。代表でない隣のセルとの
    /// 結合は、その代表との結合にする（対称な解では中性子束が同じ）。
    pub fn homogenize(
        mesh: &CoarseMesh,
        orbits: &Orbits,
        materials: &[&Material],
        volumes: &[f64],
        flux: &[Vec<f64>],
        state: &[usize],
        current: &[f64],
    ) -> Result<Self, CmfdError> {
        let ng = materials[0].groups();
        let cells = orbits.representatives.len();
        // [代表] のセルの領域（番号の順）
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); cells];
        for region in 0..materials.len() {
            if let Some(c) = orbits.position(mesh.cell(region)) {
                members[c].push(region);
            }
        }
        let mut volume = vec![0.0; cells];
        let mut phi = vec![0.0; cells * ng];
        // 反応率の和（後で中性子束で割る）
        let mut total = vec![0.0; cells * ng];
        let mut scatter = vec![0.0; cells * ng * ng];
        let mut fission: Vec<Vec<Spectrum>> = vec![Vec::new(); cells];
        // セルごとに並列に足す。1 つのセルの中では領域の番号の順に足すので、結果は分け方によらない
        volume
            .par_iter_mut()
            .zip(phi.par_chunks_mut(ng))
            .zip(total.par_chunks_mut(ng))
            .zip(scatter.par_chunks_mut(ng * ng))
            .zip(fission.par_iter_mut())
            .zip(&members)
            .for_each(|(((((volume, phi), total), scatter), terms), regions)| {
                for &region in regions {
                    let m = materials[region];
                    let (v, f) = (volumes[region], &flux[state[region]]);
                    *volume += v;
                    for g in 0..ng {
                        phi[g] += v * f[g];
                        total[g] += v * m.total[g] * f[g];
                        for to in 0..ng {
                            scatter[to * ng + g] += v * m.scatter[g][to] * f[g];
                        }
                    }
                    for term in &m.fission_terms {
                        // 同じスペクトルの項はまとめる
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
            });
        for c in 0..cells {
            for g in 0..ng {
                let reaction = phi[c * ng + g];
                if reaction <= 0.0 || !reaction.is_finite() {
                    return Err(CmfdError::NonPositiveFlux {
                        cell: orbits.representatives[c],
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
        scatter
            .par_chunks_mut(ng * ng)
            .zip(phi.par_chunks(ng))
            .for_each(|(scatter, phi)| {
                for row in scatter.chunks_mut(ng) {
                    for (s, f) in row.iter_mut().zip(phi) {
                        *s /= f;
                    }
                }
            });

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
        for (c, &cell) in orbits.representatives.iter().enumerate() {
            for (face, _) in mesh.cell_faces(cell) {
                let net = &current[face * ng..][..ng];
                match mesh.face_cells(face) {
                    FaceCells::Inner { low, high } => {
                        // 面の正の向きの流れ J = -D̃ (φ_h - φ_l) - D̂ (φ_h + φ_l)（面の長さあたり）。
                        // D̃ = 2 D_l D_h / (w_l D_h + w_h D_l)（w は面に垂直な方向のセルの幅）。このセルから見て、
                        // 外向きの流れ s J（s は外向きが面の正の向きなら 1）を -D̃ (φ_隣 - φ_自) - s D̂ (φ_自 + φ_隣) で表す
                        let (length, wl, wh) = mesh.face_geometry(face);
                        let (outward, other, own_width, other_width) = if low == cell {
                            (1.0, high, wl, wh)
                        } else {
                            (-1.0, low, wh, wl)
                        };
                        let o = orbits.index[other];
                        let mut on_other = vec![0.0; ng];
                        for g in 0..ng {
                            let (dc, dn) =
                                (diffusion(c, g, own_width), diffusion(o, g, other_width));
                            let coupling = 2.0 * dc * dn / (own_width * dn + other_width * dc);
                            let (fc, fo) = (phi[c * ng + g], phi[o * ng + g]);
                            // s D̂
                            let correction =
                                -(outward * net[g] / length + coupling * (fo - fc)) / (fc + fo);
                            diagonal[c * ng + g] += length * (coupling - correction);
                            on_other[g] = length * (-coupling - correction);
                        }
                        // 隣が自分の像（対称線をはさむ）なら、隣の中性子束は自分と同じ
                        if o == c {
                            for (d, a) in diagonal[c * ng..][..ng].iter_mut().zip(&on_other) {
                                *d += a;
                            }
                        } else {
                            neighbors[c].push((o, on_other));
                        }
                    }
                    // 外周の面: 外向きの流れを、すべて補正の項 D̂ φ で表す
                    FaceCells::Outer { outward, .. } => {
                        for g in 0..ng {
                            diagonal[c * ng + g] += outward * net[g] / phi[c * ng + g];
                        }
                    }
                }
            }
        }
        let matrices = (0..cells)
            .map(|c| {
                let mut a: Vec<f64> = scatter[c * ng * ng..][..ng * ng]
                    .iter()
                    .map(|s| -s)
                    .collect();
                for g in 0..ng {
                    a[g * ng + g] += diagonal[c * ng + g];
                }
                a
            })
            .collect();
        Self::assemble(
            mesh.rows() == 1 || mesh.columns() == 1,
            ng,
            orbits.representatives.clone(),
            orbits.sums(&vec![1.0; orbits.index.len()]),
            phi,
            fission,
            matrices,
            neighbors,
        )
    }

    /// 均質化した値から問題を組み立てる。`matrices` は [代表] の対角のブロック（群×群、行優先）、`line` は格子が
    /// 1 行か 1 列であること
    #[allow(clippy::too_many_arguments)]
    fn assemble(
        line: bool,
        groups: usize,
        cells: Vec<usize>,
        sizes: Vec<f64>,
        flux: Vec<f64>,
        fission: Vec<Vec<Spectrum>>,
        matrices: Vec<Vec<f64>>,
        neighbors: Vec<Vec<(usize, Vec<f64>)>>,
    ) -> Result<Self, CmfdError> {
        // 代表が 1 つなら Gauss-Seidel の 1 回が直接解法と同じ。2 次元の格子は、対称性で代表がチェーンになっても
        // Gauss-Seidel のまま（ADR-0158）
        let solver = if line && cells.len() >= 2 {
            Solver::Chain(Chain::new(groups, &cells, matrices, &neighbors)?)
        } else {
            // セルごとに並列に分解し、特異なものはセルの順で最初のものを返す
            let blocks = matrices
                .into_par_iter()
                .map(|a| Lu::new(a, groups))
                .collect::<Vec<_>>()
                .into_iter()
                .zip(&cells)
                .map(|(lu, &cell)| lu.ok_or(CmfdError::Singular(cell)))
                .collect::<Result<Vec<_>, _>>()?;
            Solver::GaussSeidel { blocks, neighbors }
        };
        Ok(Self {
            groups,
            cells,
            sizes,
            flux,
            fission,
            solver,
        })
    }

    /// [代表のセル][群] の中性子束 `phi` の、すべてのセルの生成率の和（代表のセルは組のセルの数だけ数える）
    fn production(&self, phi: &[f64]) -> f64 {
        let ng = self.groups;
        (0..self.cells.len())
            .flat_map(|c| {
                let p = &phi[c * ng..][..ng];
                let size = self.sizes[c];
                self.fission[c]
                    .iter()
                    .map(move |s| size * (s.chi_sum * s.rate(p)))
            })
            .sum()
    }

    /// 固有値問題を `k` から解き始め、k と [代表のセル][群] の中性子束の比（解 / 均質化した中性子束、
    /// どちらも生成率 1 に規格化）を返す。
    pub fn solve(&self, k: f64, options: &CmfdOptions) -> Result<(f64, Vec<f64>), CmfdError> {
        match &self.solver {
            Solver::Chain(chain) => self.solve_chain(chain, k, options),
            Solver::GaussSeidel { blocks, neighbors } => {
                self.solve_gauss_seidel(blocks, neighbors, k, options)
            }
        }
    }

    /// 結合がチェーンのとき、`k` の少し上の k_s = k (1 + `WIELANDT`) にシフトした逆反復で解く（ADR-0158）。
    ///
    /// 行列が Z 行列なら、A − F/k_s も Z 行列で、k_s が固有値 k₀ より上のときだけ正則な M 行列になる。そのとき
    /// (A − F/k_s)⁻¹ は非負で、正の右辺 F φ の解は正。k_s が k₀ 以下なら、正の解 ψ があると (A − F/k_s) ψ ≥ 0 から
    /// M 行列になって矛盾するので、最初の反復の解が正にならないか、行列が特異になる。そのときはシフトなしで解き
    /// 直す。Z 行列でなければ（補正の D̂ で結合が正になる）、この判定ができないので、シフトしない
    fn solve_chain(
        &self,
        chain: &Chain,
        k: f64,
        options: &CmfdOptions,
    ) -> Result<(f64, Vec<f64>), CmfdError> {
        if !chain.z_matrix {
            tracing::debug!("coarse-mesh matrix is not an irreducible Z-matrix; solving unshifted");
            return self.inverse_iteration(chain, k, 0.0, options);
        }
        let shift = 1.0 / (k * (1.0 + WIELANDT));
        match self.inverse_iteration(chain, k, shift, options) {
            Err(CmfdError::ShiftBelowEigenvalue) => {
                tracing::debug!("coarse-mesh shift below the eigenvalue; solving unshifted");
                self.inverse_iteration(chain, k, 0.0, options)
            }
            result => result,
        }
    }

    /// べき乗法の各反復で (A − μ F) φ = F φ_前 をブロック三重対角の直接解法で解く（`shift` = μ = 1/k_s の逆反復）。
    /// 解は φ / (1/k − μ) に近づくので、生成率 p から 1/k = μ + 1/p
    fn inverse_iteration(
        &self,
        chain: &Chain,
        k: f64,
        shift: f64,
        options: &CmfdOptions,
    ) -> Result<(f64, Vec<f64>), CmfdError> {
        let ng = self.groups;
        let cells = self.cells.len();
        let start = self.production(&self.flux);
        let mut phi: Vec<f64> = self.flux.iter().map(|f| f / start).collect();
        let mut k = k;
        let factors = match chain.factor(&self.fission, shift, &self.cells) {
            Err(CmfdError::SingularElimination(_)) if shift != 0.0 => {
                return Err(CmfdError::ShiftBelowEigenvalue);
            }
            factors => factors?,
        };
        let mut k_history = Remaining::default();
        let mut flux_history = Remaining::default();
        let mut next = vec![0.0; cells * ng];
        for iteration in 1..=options.max_iterations {
            for c in 0..cells {
                let rhs = &mut next[c * ng..][..ng];
                rhs.fill(0.0);
                for s in &self.fission[c] {
                    let rate = s.rate(&phi[c * ng..][..ng]);
                    for (r, x) in rhs.iter_mut().zip(&s.chi) {
                        *r += x * rate;
                    }
                }
            }
            chain.substitute(&factors, &mut next);
            // 割る前の解が正であることを確かめる（負の生成率で割ると、負の解も正に見える）
            let rate = self.production(&next);
            let new_k = 1.0 / (shift + 1.0 / rate);
            if let Some(j) = next.iter().position(|&f| f <= 0.0 || !f.is_finite()) {
                // Z 行列なら、k_s が固有値より上のときの解は正なので、最初の解の符号でシフトを判定する
                if shift != 0.0 && iteration == 1 {
                    return Err(CmfdError::ShiftBelowEigenvalue);
                }
                return Err(CmfdError::NonPositiveSolution {
                    cell: self.cells[j / ng],
                    group: j % ng,
                    flux: next[j],
                });
            }
            for f in &mut next {
                *f /= rate;
            }
            let change = next
                .iter()
                .zip(&phi)
                .map(|(a, b)| (a - b).abs() / a.abs().max(b.abs()))
                .fold(0.0, f64::max);
            std::mem::swap(&mut phi, &mut next);
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

    /// べき乗法の各反復で、ブロック Gauss-Seidel をセルの順に1回だけ行う。収束すると
    /// Gauss-Seidel でも解が変わらないので、内側の反復を収束させなくても固有値問題の解になる。
    fn solve_gauss_seidel(
        &self,
        blocks: &[Lu],
        neighbors: &[Vec<(usize, Vec<f64>)>],
        k: f64,
        options: &CmfdOptions,
    ) -> Result<(f64, Vec<f64>), CmfdError> {
        let ng = self.groups;
        let cells = blocks.len();
        let start = self.production(&self.flux);
        let mut phi: Vec<f64> = self.flux.iter().map(|f| f / start).collect();
        let mut k = k;
        let mut k_history = Remaining::default();
        let mut flux_history = Remaining::default();
        let mut rhs = vec![0.0; ng];
        for iteration in 1..=options.max_iterations {
            let previous = phi.clone();
            let inverse_k = 1.0 / k;
            for c in 0..cells {
                rhs.fill(0.0);
                for s in &self.fission[c] {
                    let rate = s.rate(&previous[c * ng..][..ng]) * inverse_k;
                    for (r, x) in rhs.iter_mut().zip(&s.chi) {
                        *r += x * rate;
                    }
                }
                for (other, coefficients) in &neighbors[c] {
                    for g in 0..ng {
                        rhs[g] -= coefficients[g] * phi[other * ng + g];
                    }
                }
                blocks[c].solve(&mut rhs);
                phi[c * ng..][..ng].copy_from_slice(&rhs);
            }
            let rate = self.production(&phi);
            let new_k = k * rate;
            for f in &mut phi {
                *f /= rate;
            }
            if let Some(j) = phi.iter().position(|&f| f <= 0.0 || !f.is_finite()) {
                return Err(CmfdError::NonPositiveSolution {
                    cell: self.cells[j / ng],
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
mod tests;
