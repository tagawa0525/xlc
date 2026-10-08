//! 燃料を黒体にした固定源の MOC による、ピンごと・群ごとの Dancoff 係数。
//!
//! 燃料を黒体（全断面積が無限大）にし、燃料の外の各領域には、その群の全断面積 Σt と、無限に広がれば
//! スカラー中性子束が 1 になる一様な等方線源 Σt（角度中性子束の平衡値 ψ0 = 1/4π、散乱なし）を置く。
//! 孤立したピンでは、燃料に入る角度中性子束はどの方向でも ψ0 になる。格子では隣のピンの影で減るので、
//! ピンの燃料の吸収率 R と、同じトラックで燃料に入るたびに ψ0 を数えた値 R_iso から C = 1 − R / R_iso とする
//! （同じ求積で数えるので、角度の離散化の誤差は打ち消し合う）。散乱がないので、反復するのは境界の
//! 角度中性子束だけである。

use std::f64::consts::PI;

use super::{POLAR_SIN, Sweeper, direction_weight, sum_tallies};
use crate::cmfd::{CmfdError, CoarseMesh};
use crate::material::Material;
use crate::tracking::Tracking;

/// 反復の打ち切り条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DancoffOptions {
    /// 境界の角度中性子束の変化の上限（ψ0 に対する比）。C の誤差もこの程度になる
    pub tolerance: f64,
    pub max_iterations: usize,
}

impl Default for DancoffOptions {
    fn default() -> Self {
        Self {
            tolerance: 1.0e-8,
            max_iterations: 10_000,
        }
    }
}

/// Dancoff 係数の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum DancoffError {
    #[error("{given} values for {regions} regions")]
    Regions { given: usize, regions: usize },
    #[error("no region outside the fuel")]
    NoModerator,
    #[error(
        "region {0}: total cross sections must be positive and finite, or zero in every group (void)"
    )]
    Total(usize),
    #[error("fuel in cell {0} touches the cell boundary; the pin must be surrounded by non-fuel")]
    FuelOnBoundary(usize),
    #[error("cell {cell}: negative Dancoff factor {value} beyond the convergence tolerance")]
    Negative { cell: usize, value: f64 },
    #[error("region {0} is in a lump but not black")]
    NotBlack(usize),
    #[error("lump {lump} is outside the {lumps} lumps")]
    Lump { lump: usize, lumps: usize },
    #[error("lump {0} has no face open to non-black regions")]
    Unexposed(usize),
    #[error("boundary fluxes not converged in {0} iterations")]
    NotConverged(usize),
    #[error(transparent)]
    Cmfd(#[from] CmfdError),
}

/// セルごとの群ごとの Dancoff 係数（燃料のないセルは `None`）。`totals` は [領域][群] の全断面積 (1/cm) で、
/// 燃料の領域（`fuel` が真）の値は使わない。
pub fn dancoff(
    tracking: &Tracking,
    totals: &[Vec<f64>],
    fuel: &[bool],
    options: &DancoffOptions,
) -> Result<Vec<Option<Vec<f64>>>, DancoffError> {
    let lattice = &tracking.lattice;
    let regions = lattice.volumes().len();
    for given in [totals.len(), fuel.len()] {
        if given != regions {
            return Err(DancoffError::Regions { given, regions });
        }
    }
    // 群の数は燃料の外の領域から取る（燃料の行は使わない）
    let ng = (0..regions)
        .find(|&r| !fuel[r])
        .map(|r| totals[r].len())
        .ok_or(DancoffError::NoModerator)?;
    for (region, t) in totals.iter().enumerate() {
        let valid = |v: &f64| *v > 0.0 && v.is_finite();
        let void = t.iter().all(|&v| v == 0.0);
        if !fuel[region] && (t.len() != ng || !(void || t.iter().all(valid))) {
            return Err(DancoffError::Total(region));
        }
    }
    let cells = lattice.cells();
    let cell_of: Vec<usize> = (0..regions).map(|r| lattice.region_cell(r)).collect();

    // 孤立したピンでの吸収率 R_iso / ψ0: 燃料に入るたびに方向の重みを足す
    let mut isolated = vec![0.0; cells];
    for track in &tracking.tracks {
        let weight = track_weight(tracking, track.azimuth);
        let n = track.segments.len();
        for forward in [true, false] {
            let mut previous: Option<usize> = None;
            for step in 0..n {
                let s = if forward { step } else { n - 1 - step };
                let region = track.segments[s].region;
                if fuel[region] {
                    let cell = cell_of[region];
                    // 燃料には同じセルの燃料の外の領域から入ること。格子の端やセルの境目から直接入ると、
                    // 孤立したピンへの入射として数えられない（逆向きのスイープで、出る側も確かめる）
                    match previous {
                        Some(p) if fuel[p] && cell_of[p] == cell => {}
                        Some(p) if !fuel[p] && cell_of[p] == cell => isolated[cell] += weight,
                        _ => return Err(DancoffError::FuelOnBoundary(cell)),
                    }
                }
                previous = Some(region);
            }
        }
    }

    let psi0 = 1.0 / (4.0 * PI);
    let delta = absorbed_by_black(tracking, totals, fuel, ng, options)?;
    let mut absorbed = vec![vec![0.0; ng]; cells];
    for region in (0..regions).filter(|&r| fuel[r]) {
        for (a, d) in absorbed[cell_of[region]]
            .iter_mut()
            .zip(&delta[region * ng..][..ng])
        {
            *a += d;
        }
    }
    absorbed
        .into_iter()
        .zip(&isolated)
        .enumerate()
        .map(|(cell, (rates, &iso))| {
            (iso > 0.0)
                .then(|| {
                    rates
                        .iter()
                        .map(|&r| factor(cell, r, psi0 * iso, options.tolerance))
                        .collect()
                })
                .transpose()
        })
        .collect()
}

/// 吸収率 `absorbed` と孤立したときの入射率 `isolated`（ψ0 を掛けたもの）からの Dancoff 係数 C = 1 − R/R_iso。
/// 孤立した塊では C = 0 が丸めや収束の判定の幅 `tolerance` で負になる。その幅の中なら 0 とする
fn factor(owner: usize, absorbed: f64, isolated: f64, tolerance: f64) -> Result<f64, DancoffError> {
    let c = 1.0 - absorbed / isolated;
    if c >= 0.0 {
        Ok(c)
    } else if c >= -tolerance {
        Ok(0.0)
    } else {
        Err(DancoffError::Negative {
            cell: owner,
            value: c,
        })
    }
}

/// 黒体（`black` が真の領域。全断面積が無限大）の外の各領域に、その群の全断面積 Σt と、無限に広がれば
/// スカラー中性子束が 1 になる一様な等方線源 Σt（ψ0 = 1/4π、散乱なし）を置いた固定源の MOC の、
/// [領域][群] の黒体の吸収率（減衰の重み付き和）。散乱がないので、反復するのは境界の角度中性子束だけ。
fn absorbed_by_black(
    tracking: &Tracking,
    totals: &[Vec<f64>],
    black: &[bool],
    ng: usize,
    options: &DancoffOptions,
) -> Result<Vec<f64>, DancoffError> {
    let regions = black.len();
    let materials: Vec<Material> = (0..regions)
        .map(|region| Material {
            total: if black[region] {
                vec![f64::INFINITY; ng]
            } else {
                totals[region].clone()
            },
            scatter: vec![vec![0.0; ng]; ng],
            anisotropic: Vec::new(),
            fission: vec![0.0; ng],
            nu_fission: vec![0.0; ng],
            fission_terms: Vec::new(),
        })
        .collect();
    let references: Vec<&Material> = materials.iter().collect();
    let mesh = CoarseMesh::new(tracking);
    let sweeper = Sweeper::new(tracking, &references, ng, &mesh, &[(0, 0)]);
    let psi0 = 1.0 / (4.0 * PI);
    let equilibrium: Vec<Vec<f64>> = (0..regions)
        .map(|region| vec![if black[region] { 0.0 } else { psi0 }; ng])
        .collect();
    // 最初は境界から ψ0 が入ってくるとする
    let mut outgoing = sweeper.isotropic_outgoing(&vec![vec![1.0; ng]; regions]);
    let mut next = outgoing.clone();
    let mut tallies = sweeper.tallies(regions);
    for _ in 0..options.max_iterations {
        sweeper.sweep(&equilibrium, &outgoing, &mut next, &mut tallies)?;
        let change = outgoing
            .iter()
            .zip(&next)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        std::mem::swap(&mut outgoing, &mut next);
        if change <= options.tolerance * psi0 {
            return Ok(sum_tallies(&tallies, ng).delta);
        }
    }
    Err(DancoffError::NotConverged(options.max_iterations))
}

/// トラックの1つの向きの、極角についての方向の重みの和 Σ_p 4π (w_a/2) ω_p sinθ_p δ_a
fn track_weight(tracking: &Tracking, azimuth: usize) -> f64 {
    (0..POLAR_SIN.len())
        .map(|p| direction_weight(tracking, azimuth, p))
        .sum()
}

/// 黒体にした塊（燃料でない共鳴の吸収体の領域など）の、群ごとの Dancoff 係数と平均弦長。
#[derive(Debug, Clone, PartialEq)]
pub struct Lump {
    pub dancoff: Vec<f64>,
    /// 平均弦長 l̄ = 4V/S (cm)。S は黒体でない領域に接する面（露出した面）
    pub chord: f64,
}

/// 塊ごとの Dancoff 係数と平均弦長。`black` は [領域] の黒体にするか、`lump` は [領域] の塊の番号（塊でない
/// 領域は `None`。塊は黒体の領域だけ）で、塊は 0..`lumps`。塊でない黒体の領域は、影を作るだけで結果を求めない。
/// 黒体の領域どうしの間を移る中性子は数えず、黒体でない領域から塊に入るたびに、孤立した塊への
/// 入射として数える（隣のセルとつながる支持格子の枠のように、露出した面が一部だけの塊にも使える）。
/// 等方な中性子束 1 での入射率 R_iso = S/4 から l̄ = 4V/S = V/R_iso とする。
pub fn lump_dancoff(
    tracking: &Tracking,
    totals: &[Vec<f64>],
    black: &[bool],
    lump: &[Option<usize>],
    lumps: usize,
    options: &DancoffOptions,
) -> Result<Vec<Lump>, DancoffError> {
    let regions = tracking.lattice.regions();
    if totals.len() != regions {
        return Err(DancoffError::Regions {
            given: totals.len(),
            regions,
        });
    }
    let exposure = exposure(tracking, black, lump, lumps)?;
    let ng = (0..regions)
        .find(|&r| !black[r])
        .map(|r| totals[r].len())
        .ok_or(DancoffError::NoModerator)?;
    for (region, t) in totals.iter().enumerate() {
        let valid = |v: &f64| *v > 0.0 && v.is_finite();
        let void = t.iter().all(|&v| v == 0.0);
        if !black[region] && (t.len() != ng || !(void || t.iter().all(valid))) {
            return Err(DancoffError::Total(region));
        }
    }
    let delta = absorbed_by_black(tracking, totals, black, ng, options)?;
    let mut absorbed = vec![vec![0.0; ng]; lumps];
    for (region, l) in lump.iter().enumerate() {
        if let Some(l) = *l {
            for (a, d) in absorbed[l].iter_mut().zip(&delta[region * ng..][..ng]) {
                *a += d;
            }
        }
    }
    absorbed
        .iter()
        .zip(&exposure)
        .enumerate()
        .map(|(l, (rates, &(iso, volume)))| {
            Ok(Lump {
                dancoff: rates
                    .iter()
                    .map(|&r| factor(l, r, iso, options.tolerance))
                    .collect::<Result<_, _>>()?,
                chord: volume / iso,
            })
        })
        .collect()
}

/// 塊ごとの平均弦長 l̄ = 4V/S (cm)（[`lump_dancoff`] の `chord`）。トラックの幾何だけで決まり、全断面積も MOC も
/// 要らない。引数は [`lump_dancoff`] と同じ。
pub fn lump_chords(
    tracking: &Tracking,
    black: &[bool],
    lump: &[Option<usize>],
    lumps: usize,
) -> Result<Vec<f64>, DancoffError> {
    Ok(exposure(tracking, black, lump, lumps)?
        .into_iter()
        .map(|(iso, volume)| volume / iso)
        .collect())
}

/// 塊ごとの、孤立した塊への入射率 R_iso（ψ0 を掛けたもの。正であること）と体積。引数を確かめる。
fn exposure(
    tracking: &Tracking,
    black: &[bool],
    lump: &[Option<usize>],
    lumps: usize,
) -> Result<Vec<(f64, f64)>, DancoffError> {
    let regions = tracking.lattice.regions();
    for given in [black.len(), lump.len()] {
        if given != regions {
            return Err(DancoffError::Regions { given, regions });
        }
    }
    if let Some(&l) = lump.iter().flatten().find(|&&l| l >= lumps) {
        return Err(DancoffError::Lump { lump: l, lumps });
    }
    if let Some(region) = (0..regions).find(|&r| lump[r].is_some() && !black[r]) {
        return Err(DancoffError::NotBlack(region));
    }
    // 孤立した塊への入射率 R_iso / ψ0: 黒体でない領域から塊に入るたびに方向の重みを足す
    // （外周から直接入るのは、鏡面反射の向こうの続きの塊からなので数えない）
    let mut isolated = vec![0.0; lumps];
    for track in &tracking.tracks {
        let weight = track_weight(tracking, track.azimuth);
        let n = track.segments.len();
        for forward in [true, false] {
            let mut previous: Option<usize> = None;
            for step in 0..n {
                let s = if forward { step } else { n - 1 - step };
                let region = track.segments[s].region;
                if let (Some(l), Some(p)) = (lump[region], previous)
                    && !black[p]
                {
                    isolated[l] += weight;
                }
                previous = Some(region);
            }
        }
    }
    let psi0 = 1.0 / (4.0 * PI);
    let volumes = tracking.volumes(regions);
    let mut volume = vec![0.0; lumps];
    for (region, l) in lump.iter().enumerate() {
        if let Some(l) = *l {
            volume[l] += volumes[region];
        }
    }
    isolated
        .iter()
        .zip(volume)
        .enumerate()
        .map(|(l, (&iso, v))| {
            let iso = psi0 * iso;
            if iso <= 0.0 {
                return Err(DancoffError::Unexposed(l));
            }
            Ok((iso, v))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};

    fn options() -> DancoffOptions {
        DancoffOptions::default()
    }

    /// ピンと水のセルを `layout`（0: ピン、1: 水）に並べ、燃料以外の群の全断面積を `total` とした Dancoff 係数
    fn factors(pitch: f64, layout: Vec<Vec<usize>>, total: &[f64]) -> Vec<Option<Vec<f64>>> {
        let pin = PinCell::new(pitch, vec![0.4096]).unwrap();
        let water = PinCell::new(pitch, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, water], layout).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let regions = lattice.volumes().len();
        // ピンのセルは (燃料, 水) の2領域、水のセルは1領域
        let mut fuel = Vec::with_capacity(regions);
        for cell in 0..lattice.cells() {
            let range = lattice.cell_regions(cell);
            fuel.push(range.len() == 2);
            fuel.extend(std::iter::repeat_n(false, range.len() - 1));
        }
        let totals = vec![total.to_vec(); regions];
        dancoff(&tracking, &totals, &fuel, &options()).unwrap()
    }

    #[test]
    fn distant_pins_do_not_shadow_each_other() {
        // 隣のピンまで減速材の光学的厚さが 20 程度あるので、影はほぼない
        let c = factors(20.0, vec![vec![0]], &[1.0, 2.0]);
        for value in c[0].as_ref().unwrap() {
            assert!(value.abs() < 1e-6, "{value}");
        }
    }

    #[test]
    fn transparent_moderator_shadows_completely() {
        let c = factors(1.26, vec![vec![0]], &[1.0e-6]);
        let value = c[0].as_ref().unwrap()[0];
        assert!(value > 0.999, "{value}");
    }

    #[test]
    fn repeated_pins_match_the_single_cell_lattice() {
        let total = [0.7, 1.5];
        let single = factors(1.26, vec![vec![0]], &total);
        let four = factors(1.26, vec![vec![0, 0], vec![0, 0]], &total);
        let reference = single[0].as_ref().unwrap();
        assert!(reference[0] > reference[1], "{reference:?}");
        assert!(
            reference.iter().all(|&c| c > 0.0 && c < 1.0),
            "{reference:?}"
        );
        for cell in &four {
            for (c, r) in cell.as_ref().unwrap().iter().zip(reference) {
                assert!((c - r).abs() < 2e-3, "{c} vs {r}");
            }
        }
    }

    #[test]
    fn void_gap_matches_a_nearly_transparent_gap() {
        // 燃料・ギャップ・被覆管・水のピン
        let cell = PinCell::new(1.26, vec![0.4096, 0.418, 0.475]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.02).unwrap();
        let fuel = [true, false, false, false];
        let factor = |gap: f64| {
            let totals = vec![vec![], vec![gap], vec![0.3], vec![1.2]];
            dancoff(&tracking, &totals, &fuel, &options()).unwrap()[0]
                .clone()
                .unwrap()[0]
        };
        let (void, thin) = (factor(0.0), factor(1.0e-9));
        assert!((void - thin).abs() < 1e-7, "{void} vs {thin}");
    }

    #[test]
    fn fuel_rows_are_not_read() {
        // 燃料の行は使わないので、空でもよい
        let cell = PinCell::new(1.26, vec![0.4096]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.05).unwrap();
        let c = dancoff(&tracking, &[vec![], vec![1.2]], &[true, false], &options()).unwrap();
        assert_eq!(c[0].as_ref().map(Vec::len), Some(1));
    }

    #[test]
    fn isolated_annulus_has_its_own_chord() {
        // 光学的に厚い水の中の孤立した円環（内径 0.5、外径 0.55）: C ≈ 0、l̄ = 4V/S = 2(r₂ − r₁)。
        // 水が薄いと、内側の水は中性子が入れず空になり、外側も遠くの円環の影を受ける
        let cell = PinCell::new(4.0, vec![0.5, 0.55]).unwrap();
        let lattice = Lattice::pin_cell(cell);
        let tracking = Tracking::new(&lattice, 32, 0.005).unwrap();
        let totals = vec![vec![5.0]; 3];
        let lumps = lump_dancoff(
            &tracking,
            &totals,
            &[false, true, false],
            &[None, Some(0), None],
            1,
            &options(),
        )
        .unwrap();
        assert!(lumps[0].dancoff[0] < 0.05, "{lumps:?}");
        assert!((lumps[0].chord - 0.1).abs() < 1e-3 * 0.1, "{lumps:?}");
    }

    #[test]
    fn connected_frames_see_only_their_open_faces() {
        // 2×2 のセルの外周の枠（厚さ 0.05）は隣の枠とつながる。セルごとの塊の l̄ = 4V/S は、露出した内側の面
        // S = 4(p − 2t) だけで決まる
        let (p, t) = (1.26, 0.05);
        let cell = PinCell::new(p, vec![]).unwrap().with_frame(t).unwrap();
        let lattice = Lattice::new(vec![cell], vec![vec![0, 0], vec![0, 0]]).unwrap();
        let tracking = Tracking::new(&lattice, 32, 0.005).unwrap();
        let totals = vec![vec![0.3]; lattice.regions()];
        // セルごとに（水, 枠）の2領域
        let lump: Vec<Option<usize>> = (0..8).map(|r| (r % 2 == 1).then_some(r / 2)).collect();
        let black: Vec<bool> = lump.iter().map(Option::is_some).collect();
        let lumps = lump_dancoff(&tracking, &totals, &black, &lump, 4, &options()).unwrap();
        let expected = (p * p - (p - 2.0 * t).powi(2)) / (p - 2.0 * t);
        for l in &lumps {
            assert!((l.chord - expected).abs() < 1e-3 * expected, "{lumps:?}");
            assert!((0.0..=1.0).contains(&l.dancoff[0]), "{lumps:?}");
            assert!(
                (l.dancoff[0] - lumps[0].dancoff[0]).abs() < 1e-6,
                "{lumps:?}"
            );
        }
    }

    #[test]
    fn black_regions_outside_lumps_shadow_them() {
        // 4つの枠をすべて黒体にし、1つだけを塊として数えると、4つとも塊にしたときと同じ結果
        let (p, t) = (1.26, 0.05);
        let cell = PinCell::new(p, vec![]).unwrap().with_frame(t).unwrap();
        let lattice = Lattice::new(vec![cell], vec![vec![0, 0], vec![0, 0]]).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.01).unwrap();
        let totals = vec![vec![0.3]; lattice.regions()];
        let black: Vec<bool> = (0..8).map(|r| r % 2 == 1).collect();
        let all: Vec<Option<usize>> = (0..8).map(|r| (r % 2 == 1).then_some(r / 2)).collect();
        let one: Vec<Option<usize>> = (0..8).map(|r| (r == 1).then_some(0)).collect();
        let expected = lump_dancoff(&tracking, &totals, &black, &all, 4, &options()).unwrap();
        let lumps = lump_dancoff(&tracking, &totals, &black, &one, 1, &options()).unwrap();
        assert_eq!(lumps.len(), 1);
        assert!(
            (lumps[0].chord - expected[0].chord).abs() < 1e-12,
            "{lumps:?}"
        );
        assert!(
            (lumps[0].dancoff[0] - expected[0].dancoff[0]).abs() < 1e-12,
            "{lumps:?} vs {expected:?}"
        );
    }

    #[test]
    fn chords_come_from_the_tracks_alone() {
        // 塊の平均弦長は、Dancoff 係数を求めたときの値と同じ。複数の塊を1つにまとめると V と R_iso の和の比
        let (p, t) = (1.26, 0.05);
        let cell = PinCell::new(p, vec![]).unwrap().with_frame(t).unwrap();
        let lattice = Lattice::new(vec![cell], vec![vec![0, 0], vec![0, 0]]).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.01).unwrap();
        let totals = vec![vec![0.3]; lattice.regions()];
        let lump: Vec<Option<usize>> = (0..8).map(|r| (r % 2 == 1).then_some(r / 2)).collect();
        let black: Vec<bool> = lump.iter().map(Option::is_some).collect();
        let expected = lump_dancoff(&tracking, &totals, &black, &lump, 4, &options()).unwrap();
        let chords = lump_chords(&tracking, &black, &lump, 4).unwrap();
        assert_eq!(chords.len(), 4);
        for (c, e) in chords.iter().zip(&expected) {
            assert!((c - e.chord).abs() < 1e-12, "{chords:?} vs {expected:?}");
        }
        let merged: Vec<Option<usize>> = lump.iter().map(|l| l.map(|_| 0)).collect();
        let whole = lump_chords(&tracking, &black, &merged, 1).unwrap();
        let volumes = tracking.volumes(lattice.regions());
        let volume: f64 = (0..8).filter(|r| r % 2 == 1).map(|r| volumes[r]).sum();
        let entries: f64 = (0..4).map(|l| volumes[2 * l + 1] / expected[l].chord).sum();
        assert!((whole[0] - volume / entries).abs() < 1e-12, "{whole:?}");
    }

    #[test]
    fn rejects_lumps_that_are_not_black() {
        let cell = PinCell::new(4.0, vec![0.5, 0.55]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.05).unwrap();
        let totals = vec![vec![5.0]; 3];
        let result = lump_dancoff(
            &tracking,
            &totals,
            &[false, false, false],
            &[None, Some(0), None],
            1,
            &options(),
        );
        assert!(
            matches!(result, Err(DancoffError::NotBlack(1))),
            "{result:?}"
        );
    }

    #[test]
    fn rejects_lump_numbers_out_of_range() {
        let cell = PinCell::new(4.0, vec![0.5, 0.55]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.05).unwrap();
        let totals = vec![vec![5.0]; 3];
        let result = lump_dancoff(
            &tracking,
            &totals,
            &[false, true, false],
            &[None, Some(1), None],
            1,
            &options(),
        );
        assert!(
            matches!(result, Err(DancoffError::Lump { lump: 1, lumps: 1 })),
            "{result:?}"
        );
    }

    #[test]
    fn rejects_fuel_clipped_by_the_cell() {
        // 燃料の半径がセルの半ピッチより大きいと、隣のセルから直接燃料に入る
        let pin = PinCell::new(1.0, vec![0.6]).unwrap();
        let water = PinCell::new(1.0, vec![]).unwrap();
        let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 1]]).unwrap();
        let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
        let mut fuel = vec![false; lattice.volumes().len()];
        fuel[0] = true;
        let totals = vec![vec![1.0]; fuel.len()];
        let result = dancoff(&tracking, &totals, &fuel, &options());
        assert!(
            matches!(result, Err(DancoffError::FuelOnBoundary(0))),
            "{result:?}"
        );
    }

    #[test]
    fn cells_without_fuel_have_no_factors() {
        let c = factors(1.26, vec![vec![0, 1], vec![1, 0]], &[1.0]);
        assert!(c[0].is_some() && c[3].is_some());
        assert!(c[1].is_none() && c[2].is_none());
    }
}
