//! CMFD（粗メッシュ有限差分）による MOC の加速。
//!
//! 粗メッシュは格子のセルとする。スイープでセルの面を通る正味の中性子の流れを数え、
//! 拡散の結合係数を流れに合うよう補正した粗メッシュの固有値問題を解いて、
//! 細かい領域の中性子束をセルごとの比で更新する。

use crate::tracking::Tracking;

/// 格子のセルを粗メッシュとし、トラックが外周で通る面を前もって求めたもの。
///
/// 面は、縦の面（x 方向に通る）を先に、横の面（y 方向に通る）を後に番号を付ける。
/// 縦の面 (r, i) はセルの行 r の、左から i 番目の格子線の上（i = 0..=n）、横の面 (j, c) は
/// セルの列 c の、上から j 番目の格子線の上（j = 0..=n）にある。面の正の向きは、
/// 縦の面では +x（右）、横の面では -y（下、行の番号が増える向き）とする。
///
/// 内側の面の通過は前もって持たず、スイープで隣り合う区間のセルを比べて求める。スイープは
/// 減衰の表の読み込みでメモリの転送速度に頭打ちになるので、読む量を増やさないため。
#[derive(Debug, Clone, PartialEq)]
pub struct CoarseMesh {
    /// 1辺のセルの数
    n: usize,
    /// [領域] の属するセル
    cell_of_region: Vec<usize>,
    /// [トラック] の始点で入る面と終点で出る面（順方向に通るときの向き付き）
    ends: Vec<[(usize, f64); 2]>,
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
        let n = lattice.cells_per_side();
        let mut cell_of_region = vec![0; lattice.regions()];
        for cell in 0..n * n {
            for region in lattice.cell_regions(cell) {
                cell_of_region[region] = cell;
            }
        }
        let half = 0.5 * lattice.width();
        let mut mesh = Self {
            n,
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

    pub fn cells(&self) -> usize {
        self.n * self.n
    }

    pub fn faces(&self) -> usize {
        2 * self.n * (self.n + 1)
    }

    /// 領域の属するセル
    pub fn cell(&self, region: usize) -> usize {
        self.cell_of_region[region]
    }

    /// セルの4つの面と、外向きが面の正の向きなら 1、逆なら -1（左、右、上、下の順）。
    pub fn cell_faces(&self, cell: usize) -> [(usize, f64); 4] {
        let (r, c) = (cell / self.n, cell % self.n);
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

    fn vertical(&self, row: usize, line: usize) -> usize {
        row * (self.n + 1) + line
    }

    fn horizontal(&self, line: usize, column: usize) -> usize {
        self.n * (self.n + 1) + line * self.n + column
    }

    /// セル `from` から隣のセル `to` へ移るときに通る面と向き。
    pub(crate) fn path(&self, from: usize, to: usize) -> FacePath {
        let (r1, c1) = (from / self.n, from % self.n);
        let (r2, c2) = (to / self.n, to % self.n);
        let sign = |a: usize, b: usize| if b > a { 1.0 } else { -1.0 };
        match (r1.abs_diff(r2), c1.abs_diff(c2)) {
            (0, 1) => FacePath::new(&[(self.vertical(r1, c1.max(c2)), sign(c1, c2))]),
            (1, 0) => FacePath::new(&[(self.horizontal(r1.max(r2), c1), sign(r1, r2))]),
            // 角を通る: (r1, c2) を経由する経路と (r2, c1) を経由する経路に半分ずつ
            (1, 1) => FacePath::new(&[
                (self.vertical(r1, c1.max(c2)), 0.5 * sign(c1, c2)),
                (self.horizontal(r1.max(r2), c2), 0.5 * sign(r1, r2)),
                (self.horizontal(r1.max(r2), c1), 0.5 * sign(r1, r2)),
                (self.vertical(r2, c1.max(c2)), 0.5 * sign(c1, c2)),
            ]),
            _ => panic!("consecutive segments in cells {from} and {to} that do not touch"),
        }
    }

    /// トラックの端の点 `point`（外周の上）にある、セル `cell` の外周の面と、
    /// 外向きが面の正の向きなら 1、逆なら -1。
    fn outer_face(&self, cell: usize, point: (f64, f64), half: f64) -> (usize, f64) {
        let (r, c) = (cell / self.n, cell % self.n);
        // 外向きの向き。トラックは外周の角を通らないので、点に最も近い辺が点のある辺
        let sides = [
            ((point.0 + half).abs(), self.vertical(r, 0), -1.0),
            ((point.0 - half).abs(), self.vertical(r, self.n), 1.0),
            ((point.1 - half).abs(), self.horizontal(0, c), -1.0),
            ((point.1 + half).abs(), self.horizontal(self.n, c), 1.0),
        ];
        let (_, face, outward) = sides
            .into_iter()
            .min_by(|a, b| a.0.total_cmp(&b.0))
            .expect("four sides");
        (face, outward)
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
                let path = mesh.path(from, to);
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
