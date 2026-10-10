//! 対称性を使う MOC（ADR-0113）: スイープする代表のトラックと、状態を持つ代表の領域と粗メッシュのセル。
//!
//! 対称操作の群 G で移り合うトラックの組（軌道 O）から代表を 1 本ずつ選び、軌道の大きさ |O| を掛けた重みでスイープする
//! （tally の T）。解が対称なら、全体をスイープしたときの量は F(r) = (1/|G|) Σ_g S_g T(g⁻¹ r) になる。S_g はスカラーでは 1、
//! ベクトルでは g の行列、モーメントでは符号つきの置換、面の正味の中性子流では面の向きの符号。領域ごとの量（中性子束、
//! モーメント、平衡値）と CMFD のセルは、軌道の代表の分だけ持つ。代表でない領域の値は、代表をその領域へ移す操作で
//! 代表の値から求める。対称性を使わないときは、群を恒等だけにした同じ形で扱う。

use std::ops::Range;

use super::harmonics::real_harmonic;
use super::{MocError, Sums};
use crate::cmfd::CoarseMesh;
use crate::material::Material;
use crate::symmetry::{Group, Image, Orbits, SymmetryError, images};
use crate::tracking::Tracking;

/// スイープする代表のトラック、状態を持つ代表の領域とセル、tally の集め方。
pub(super) struct Representatives {
    /// スイープするトラック（番号の順）
    pub tracks: Vec<usize>,
    /// [代表] の、移り合うトラックの数（スイープの重みに掛ける）
    pub multiplicity: Vec<f64>,
    /// 並列にスイープする組（`tracks` の位置の範囲）。組はスレッドの数によらない
    pub chunks: Vec<Range<usize>>,
    /// [トラック][向き (0: 順, 1: 逆)] の終わりの角度中性子束を持つ、代表のスイープの番号（代表の位置 × 2 + 向き）
    slots: Vec<usize>,
    /// 群の操作の像（恒等を含む。対称性を使わなければ恒等だけ）
    images: Vec<TallyImage>,
    /// 領域の軌道。状態は代表の領域の分だけ持つ
    pub regions: Orbits,
    /// [領域] の、代表をその領域へ移す操作（`images` の番号）
    placement: Vec<usize>,
    /// 粗メッシュのセルの軌道（CMFD は代表のセルだけで解く）
    pub cells: Orbits,
}

/// 1つの対称操作による tally の像。
struct TallyImage {
    /// [領域] の像の領域
    regions: Vec<usize>,
    /// [領域] の、その領域に移る領域（`regions` の逆）
    sources: Vec<usize>,
    /// [面] の像の面と、面の正の向きの像が像の面の正の向きなら 1、逆なら -1
    faces: Vec<(usize, f64)>,
    /// [モーメント k] の Y_k(gΩ) = s Y_k'(Ω) の (k', s)
    moments: Vec<(usize, f64)>,
    /// (x, y) に掛ける行列 [a, b, c, d]
    matrix: [f64; 4],
}

impl Representatives {
    /// すべてのトラックをスイープし、すべての領域とセルの状態を持つ（対称性を使わない）。組は方位角ごと。`moments` は
    /// 角度中性子束のモーメントの数
    pub fn all(tracking: &Tracking, mesh: &CoarseMesh, moments: usize) -> Self {
        let n = tracking.tracks.len();
        let regions = tracking.lattice.regions();
        let identity: Vec<usize> = (0..regions).collect();
        Self {
            tracks: (0..n).collect(),
            multiplicity: vec![1.0; n],
            chunks: (0..tracking.azimuths.len())
                .map(|a| tracking.azimuth_tracks(a))
                .collect(),
            slots: (0..2 * n).collect(),
            images: vec![TallyImage {
                regions: identity.clone(),
                sources: identity,
                faces: (0..mesh.faces()).map(|f| (f, 1.0)).collect(),
                moments: (0..moments).map(|k| (k, 1.0)).collect(),
                matrix: [1.0, 0.0, 0.0, 1.0],
            }],
            regions: Orbits::identity(regions),
            placement: vec![0; regions],
            cells: Orbits::identity(mesh.cells()),
        }
    }

    /// 群 `group` で移り合うトラックの代表だけをスイープし、領域とセルの軌道の代表の状態だけを持つ（ADR-0112）。形状と
    /// トラック（`symmetry::images`）が群のどれかの操作で移り合わないか、移り合う領域の物質 `materials` が同じ物質
    /// （同じ参照）でなければエラー
    pub fn group(
        tracking: &Tracking,
        materials: &[&Material],
        mesh: &CoarseMesh,
        pairs: &[(usize, isize)],
        group: Group,
    ) -> Result<Self, MocError> {
        let images = images(tracking, &group.operations())?;
        Self::by_images(tracking, materials, mesh, pairs, &images)
    }

    /// 群の操作の像 `images`（恒等を含む）から作る
    fn by_images(
        tracking: &Tracking,
        materials: &[&Material],
        mesh: &CoarseMesh,
        pairs: &[(usize, isize)],
        images: &[Image],
    ) -> Result<Self, MocError> {
        for image in images {
            for (region, &other) in image.regions.iter().enumerate() {
                // 同じ物質（同じ参照）であること。値が等しいだけの別の物質は認めない: 自己遮蔽で別々に解いたピンは
                // 丸めの程度で違い、値で比べると対称かどうかが丸めで変わる（ADR-0113）
                if !std::ptr::eq(materials[region], materials[other]) {
                    return Err(MocError::AsymmetricMaterial {
                        operation: image.name,
                        region,
                        image: other,
                    });
                }
            }
        }
        // 軌道の代表は、軌道の中で最も小さい番号のもの（トラック、領域、セル）。群は恒等を含むので、軌道は空でない
        let smallest =
            |of: &dyn Fn(&Image) -> usize, own: usize| images.iter().map(of).fold(own, usize::min);
        let n = tracking.tracks.len();
        let representative: Vec<usize> = (0..n).map(|t| smallest(&|g| g.tracks[t].0, t)).collect();
        let tracks: Vec<usize> = (0..n).filter(|&t| representative[t] == t).collect();
        let mut position = vec![usize::MAX; n];
        for (i, &t) in tracks.iter().enumerate() {
            position[t] = i;
        }
        let multiplicity = tracks
            .iter()
            .map(|&t| {
                let mut orbit: Vec<usize> = images.iter().map(|g| g.tracks[t].0).collect();
                orbit.sort_unstable();
                orbit.dedup();
                orbit.len() as f64
            })
            .collect();
        // トラック t の向き d の終わりは、t に移る代表 r の、像で同じ向きになる向きの終わりと同じ値
        let mut slots = vec![0; 2 * n];
        for t in 0..n {
            let r = representative[t];
            let (_, same) = images
                .iter()
                .map(|g| g.tracks[r])
                .find(|&(image, _)| image == t)
                .ok_or_else(|| {
                    broken(format!(
                        "track {t} is not an image of its representative {r}"
                    ))
                })?;
            for forward in [true, false] {
                let source = if same { forward } else { !forward };
                slots[2 * t + usize::from(!forward)] = 2 * position[r] + usize::from(!source);
            }
        }
        let chunks = balanced_chunks(tracking, &tracks, tracking.azimuths.len());

        let count = tracking.lattice.regions();
        let region_representative: Vec<usize> =
            (0..count).map(|r| smallest(&|g| g.regions[r], r)).collect();
        let regions = Orbits::new(&region_representative);
        // 代表 ρ を領域 r へ移す操作（どれでもよい: 代表を動かさない操作で値が変わる成分は、対称な解では 0）
        let placement = region_representative
            .iter()
            .enumerate()
            .map(|(r, &rho)| {
                images
                    .iter()
                    .position(|g| g.regions[rho] == r)
                    .ok_or_else(|| {
                        broken(format!("region {r} is not an image of its representative"))
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let cell_representative: Vec<usize> = (0..mesh.cells())
            .map(|c| smallest(&|g| g.cells[c], c))
            .collect();
        let cells = Orbits::new(&cell_representative);
        tracing::info!(
            tracks = n,
            representative_tracks = tracks.len(),
            regions = count,
            representative_regions = regions.representatives.len(),
            cells = mesh.cells(),
            representative_cells = cells.representatives.len(),
            "symmetry"
        );
        let images = images
            .iter()
            .map(|image| tally_image(image, mesh, pairs))
            .collect::<Result<_, _>>()?;
        Ok(Self {
            tracks,
            multiplicity,
            chunks,
            slots,
            images,
            regions,
            placement,
            cells,
        })
    }

    /// トラック `track` の向き `forward` のスイープの終わりの角度中性子束を持つ、代表のスイープの番号
    pub fn slot(&self, track: usize, forward: bool) -> usize {
        self.slots[2 * track + usize::from(!forward)]
    }

    /// 領域 `region` の [モーメント k] の値を、その代表の値から求める (k', s)（値は s × 代表の k' 番目）
    pub fn moments(&self, region: usize) -> &[(usize, f64)] {
        &self.images[self.placement[region]].moments
    }

    /// 代表の位置 `representative` の領域の、全体をスイープしたときの量を `out` に足す。`each(from, f)` は、領域 `from` の
    /// 代表のトラックの組ごとの値（軌道の大きさを掛けたもの）を組の順に `f` に渡す。像の領域の値を `kind` に従って移して
    /// 平均する（対称性を使わなければ、その領域の値をそのまま足す）
    pub fn gather(
        &self,
        representative: usize,
        out: &mut [f64],
        kind: RegionKind,
        each: &RegionValues,
    ) {
        let to = self.regions.representatives[representative];
        let share = 1.0 / self.images.len() as f64;
        for image in &self.images {
            each(image.sources[to], &mut |values| match kind {
                // [モーメント][群]: Y_k(gΩ) = s Y_k'(Ω)
                RegionKind::Moments(ng) => {
                    for (k, &(source, sign)) in image.moments.iter().enumerate() {
                        for (o, v) in out[k * ng..][..ng]
                            .iter_mut()
                            .zip(&values[source * ng..][..ng])
                        {
                            *o += share * sign * v;
                        }
                    }
                }
                // [x, y][群] のベクトル
                RegionKind::Vector(ng) => turn(out, values, &image.matrix, share, ng),
            });
        }
    }

    /// 代表のトラックの、真空の領域と面の tally `sums`（軌道の大きさを掛けたもの）を、全体をスイープしたときの値にする
    /// （対称性を使わなければそのまま）。領域ごとの量（`delta` と `region_current`）は `gather` で集めてある。`void` は
    /// [領域] の真空の領域の通し番号
    pub fn symmetrize(&self, sums: &mut Sums, void: &[Option<usize>], ng: usize) {
        if self.images.len() == 1 {
            return;
        }
        let share = 1.0 / self.images.len() as f64;
        let mut void_flux = vec![0.0; sums.void_flux.len()];
        let mut void_current = vec![0.0; sums.void_current.len()];
        let mut current = vec![0.0; sums.current.len()];
        let mut face_flux = vec![0.0; sums.face_flux.len()];
        for image in &self.images {
            for (from, &to) in image.regions.iter().enumerate() {
                // 物質が同じことを確かめてあるので、真空の領域の像は真空の領域
                if let (Some(v), Some(w)) = (void[from], void[to]) {
                    for (d, s) in void_flux[w * ng..][..ng]
                        .iter_mut()
                        .zip(&sums.void_flux[v * ng..][..ng])
                    {
                        *d += share * s;
                    }
                    turn(
                        &mut void_current[w * 2 * ng..][..2 * ng],
                        &sums.void_current[v * 2 * ng..][..2 * ng],
                        &image.matrix,
                        share,
                        ng,
                    );
                }
            }
            for (from, &(to, sign)) in image.faces.iter().enumerate() {
                for g in 0..ng {
                    current[to * ng + g] += share * sign * sums.current[from * ng + g];
                    face_flux[to * ng + g] += share * sums.face_flux[from * ng + g];
                }
            }
        }
        sums.void_flux = void_flux;
        sums.void_current = void_current;
        sums.current = current;
        sums.face_flux = face_flux;
    }

    /// [代表][群] のスカラーの値 `values` を、すべての領域の [領域][群] の値にする
    pub fn expand(&self, values: &[Vec<f64>]) -> Vec<Vec<f64>> {
        self.regions
            .index
            .iter()
            .map(|&i| values[i].clone())
            .collect()
    }

    /// [代表][x, y の順に群] のベクトル `values` を、すべての領域の [領域][x, y][群] にする
    pub fn expand_vectors(&self, values: &[Vec<f64>], ng: usize) -> Vec<[Vec<f64>; 2]> {
        self.regions
            .index
            .iter()
            .zip(&self.placement)
            .map(|(&i, &operation)| {
                let mut x = vec![0.0; 2 * ng];
                turn(&mut x, &values[i], &self.images[operation].matrix, 1.0, ng);
                let y = x.split_off(ng);
                [x, y]
            })
            .collect()
    }
}

/// 領域の番号と、その領域の組ごとの値を受け取る関数を受け取り、値を組の順に渡す関数（`Representatives::gather`）
pub(super) type RegionValues<'a> = dyn Fn(usize, &mut dyn FnMut(&[f64])) + 'a;

/// 領域ごとの量の並び（`Representatives::gather`）。
#[derive(Clone, Copy)]
pub(super) enum RegionKind {
    /// [モーメント][群]（群の数）
    Moments(usize),
    /// [x, y][群]（群の数）
    Vector(usize),
}

/// 群の性質が崩れたときのエラー
fn broken(reason: String) -> SymmetryError {
    SymmetryError::Broken {
        operation: "group",
        reason,
    }
}

/// [x, y][群] のベクトル `from` を行列 `m` で移し、`share` を掛けて `to` に足す
fn turn(to: &mut [f64], from: &[f64], m: &[f64; 4], share: f64, ng: usize) {
    for g in 0..ng {
        let (x, y) = (from[g], from[ng + g]);
        to[g] += share * (m[0] * x + m[1] * y);
        to[ng + g] += share * (m[2] * x + m[3] * y);
    }
}

/// 代表のトラック `tracks` を、区間の数の和が揃うように順に `count` 個以下の組に分ける
fn balanced_chunks(tracking: &Tracking, tracks: &[usize], count: usize) -> Vec<Range<usize>> {
    let total: usize = tracks
        .iter()
        .map(|&t| tracking.tracks[t].segments.len())
        .sum();
    let mut chunks = Vec::with_capacity(count);
    let (mut start, mut sum) = (0, 0);
    for (i, &t) in tracks.iter().enumerate() {
        sum += tracking.tracks[t].segments.len();
        // 組 k は、区間の数の累積が total (k + 1) / count に達したところで閉じる
        if sum * count >= total * (chunks.len() + 1) {
            chunks.push(start..i + 1);
            start = i + 1;
        }
    }
    if start < tracks.len() {
        chunks.push(start..tracks.len());
    }
    chunks
}

/// 対称操作の像 `image` を、tally の領域、面、モーメントの像にする
fn tally_image(
    image: &Image,
    mesh: &CoarseMesh,
    pairs: &[(usize, isize)],
) -> Result<TallyImage, MocError> {
    let m = image.matrix;
    // セルの辺（左、右、上、下）の外向きの向き (x, y)（y は上向き）
    const OUTWARD: [(f64, f64); 4] = [(-1.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.0, -1.0)];
    let mut faces = vec![(usize::MAX, 0.0); mesh.faces()];
    for (cell, &other) in image.cells.iter().enumerate() {
        let (own, theirs) = (mesh.cell_faces(cell), mesh.cell_faces(other));
        for (side, (x, y)) in OUTWARD.iter().enumerate() {
            let turned = (m[0] * x + m[1] * y, m[2] * x + m[3] * y);
            // 行列の成分は 0 と ±1 なので、像の向きはちょうど辺の向きのどれか
            let Some(target) = OUTWARD.iter().position(|&o| o == turned) else {
                return Err(broken(format!("side {side} of cell {cell} has no image side")).into());
            };
            // 面の正の向きは (外向きの符号) × 外向き。像では (外向きの符号) × 像の辺の外向き
            let (face, sign) = own[side];
            let (to, their_sign) = theirs[target];
            faces[face] = (to, sign * their_sign);
        }
    }
    // Y_k(gΩ) = s Y_k'(Ω) を、いくつかの方向で値を比べて求める（l が同じものの中で）
    let samples = [(0.3, 0.4), (0.7, 1.3), (0.5, 2.9), (0.9, 4.4)];
    let mut moments = Vec::with_capacity(pairs.len());
    for &(l, mk) in pairs {
        let found = pairs
            .iter()
            .enumerate()
            .filter(|(_, (l2, _))| *l2 == l)
            .flat_map(|(k, &(_, m2))| [(k, m2, 1.0), (k, m2, -1.0)])
            .find(|&(_, m2, sign)| {
                samples.iter().all(|&(mu, phi): &(f64, f64)| {
                    let (x, y) = (phi.cos(), phi.sin());
                    let turned = (m[2] * x + m[3] * y).atan2(m[0] * x + m[1] * y);
                    let lhs = real_harmonic(l, mk, mu, turned);
                    let rhs = sign * real_harmonic(l, m2, mu, phi);
                    (lhs - rhs).abs() < 1.0e-10
                })
            });
        let Some((k, _, sign)) = found else {
            return Err(broken(format!(
                "moment ({l}, {mk}) is not a signed permutation of the others"
            ))
            .into());
        };
        moments.push((k, sign));
    }
    // 像の領域は1対1（物質と体積の確認の前提。`symmetry::images` のテスト）
    let mut sources = vec![usize::MAX; image.regions.len()];
    for (from, &to) in image.regions.iter().enumerate() {
        sources[to] = from;
    }
    if let Some(to) = sources.iter().position(|&from| from == usize::MAX) {
        return Err(broken(format!("no region maps to region {to}")).into());
    }
    Ok(TallyImage {
        regions: image.regions.clone(),
        sources,
        faces,
        moments,
        matrix: m,
    })
}
