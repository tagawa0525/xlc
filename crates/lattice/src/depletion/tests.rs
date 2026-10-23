use xlc_depletion::chain::FinalState;
use xlc_depletion::chain::{
    Chain, Channel, Decay, DecayEnergy, DecayMode, Fission, Nuclide, Product, Reaction, Yield,
};
use xlc_depletion::{CramOrder, deplete};
use xlc_format::group::GroupStructure;
use xlc_format::library::{
    ChannelState, GroupChannel, GroupReaction, Library, NuclideData, Pointwise, Shielded,
    TemperatureData,
};

use std::f64::consts::LN_2;

use super::*;
use crate::input::LatticeInput;
use crate::shielding::{FuelShielding, shield_fuel_with};
use crate::symmetry::choose;
use crate::tracking::Tracking;

const T: f64 = 600.0;
const XE135: usize = 2;
const XE136: usize = 3;
const CS135: usize = 4;
const BA140: usize = 5;
pub(crate) const U235: usize = 6;
const U236: usize = 7;
const U238: usize = 8;
/// U-235 の核分裂の Q (eV)
const Q: f64 = 1.93e8;

fn stable(name: &str, z: u32, a: u32) -> Nuclide {
    Nuclide {
        name: name.into(),
        z,
        a,
        state: 0,
        excitation: 0.0,
        decay: None,
        reactions: vec![],
        fission: None,
    }
}

fn capture(target: usize) -> Reaction {
    Reaction {
        mt: 102,
        channels: vec![Channel {
            final_state: None,
            products: vec![Product {
                nuclide: target,
                count: 1,
            }],
        }],
    }
}

/// H-1、O-16、Xe-135（β− で Cs-135、(n,γ) で Xe-136）、Xe-136、Cs-135、Ba-140、U-235（(n,γ) で U-236、核分裂で Xe-135 と
/// Ba-140 を合わせて 2 個）、U-236、U-238 のチェーン。核分裂生成物の原子の数は核分裂の 2 倍になる
pub(crate) fn chain() -> Chain {
    let mut xe135 = stable("Xe-135", 54, 135);
    xe135.decay = Some(Decay {
        half_life: 32_904.0,
        energy: DecayEnergy {
            light: 0.0,
            electromagnetic: 0.0,
            heavy: 0.0,
        },
        modes: vec![DecayMode {
            mode: "beta-".into(),
            ratio: 1.0,
            products: vec![Product {
                nuclide: CS135,
                count: 1,
            }],
        }],
    });
    xe135.reactions.push(capture(XE136));
    let mut u235 = stable("U-235", 92, 235);
    u235.reactions.push(capture(U236));
    let table = vec![
        Yield {
            nuclide: XE135,
            value: 0.06,
        },
        Yield {
            nuclide: BA140,
            value: 1.94,
        },
    ];
    u235.fission = Some(Fission {
        q: Q,
        yield_source: "U-235".into(),
        energies: vec![0.0253, 5e5],
        yields: vec![table.clone(), table],
    });
    Chain {
        nuclides: vec![
            stable("H-1", 1, 1),
            stable("O-16", 8, 16),
            xe135,
            stable("Xe-136", 54, 136),
            stable("Cs-135", 55, 135),
            stable("Ba-140", 56, 140),
            u235,
            stable("U-236", 92, 236),
            stable("U-238", 92, 238),
        ],
    }
}

/// 2 群（0.625 eV で分ける）の、分離共鳴のないライブラリ。全断面積は散乱、核分裂、捕獲の和。散乱は P1 まで
pub(crate) fn library() -> Library {
    // (名前, awr, 散乱 [[g→g']], 核分裂, 捕獲)
    let nuclide =
        |name: &str, awr: f64, scatter: [[f64; 2]; 2], fission: [f64; 2], capture: [f64; 2]| {
            let total: Vec<f64> = (0..2)
                .map(|g| scatter[g][0] + scatter[g][1] + fission[g] + capture[g])
                .collect();
            let table = Shielded {
                total: total.clone(),
                elastic: (0..2).map(|g| scatter[g][0] + scatter[g][1]).collect(),
                fission: fission.to_vec(),
                nu_fission: fission.iter().map(|f| 2.43 * f).collect(),
                capture: capture.to_vec(),
            };
            NuclideData {
                name: name.into(),
                // ZA = 1000 Z + A（履歴からの分岐が重金属を Z で選ぶ。ADR-0190）
                za: match name {
                    "H-1" => 1001.0,
                    "O-16" => 8016.0,
                    "Xe-135" => 54135.0,
                    "U-235" => 92235.0,
                    "U-238" => 92238.0,
                    _ => panic!("no ZA for {name}"),
                },
                awr,
                resolved: None,
                thermal_scattering: false,
                temperatures: vec![TemperatureData {
                    temperature: T,
                    chi: if fission[0] > 0.0 {
                        vec![1.0, 0.0]
                    } else {
                        vec![0.0, 0.0]
                    },
                    // P1 は P0 の 1 割（拡散係数の重みに要る）
                    scatter: [1.0, 0.1]
                        .iter()
                        .map(|f| {
                            scatter
                                .iter()
                                .map(|row| row.iter().map(|s| f * s).collect())
                                .collect()
                        })
                        .collect(),
                    shielded: vec![table.clone(), table.clone(), table],
                    pointwise: Pointwise::default(),
                    // 捕獲は状態に分けない 1 つの経路
                    reactions: vec![GroupReaction {
                        mt: 102,
                        channels: vec![GroupChannel {
                            final_state: None,
                            cross_section: capture.to_vec(),
                        }],
                    }],
                }],
            }
        };
    Library {
        groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
        sigma0: vec![1.0e10, 100.0, 1.0],
        pointwise_range: None,
        nuclides: vec![
            nuclide(
                "H-1",
                0.99917,
                [[3.9, 0.3], [0.0, 20.0]],
                [0.0; 2],
                [0.0, 0.3],
            ),
            nuclide(
                "O-16",
                15.858,
                [[3.7, 0.05], [0.0, 3.8]],
                [0.0; 2],
                [0.0; 2],
            ),
            nuclide(
                "Xe-135",
                133.75,
                [[5.0, 0.0], [0.0, 5.0]],
                [0.0; 2],
                [1.0, 2.0e6],
            ),
            nuclide(
                "U-235",
                233.025,
                [[10.0, 0.0], [0.0, 10.0]],
                [1.2, 500.0],
                [0.8, 90.0],
            ),
            nuclide(
                "U-238",
                236.006,
                [[9.5, 0.0], [0.0, 9.0]],
                [0.0; 2],
                [0.45, 2.7],
            ),
        ],
    }
}

pub(crate) const PIN: &str = r#"
library = "unused.msgpack"
pitch = 1.26
azimuthal = 8
spacing = 0.05
rows = ["pin"]

[material.fuel]
temperature = 600
nuclides = { U-235 = 7.0e-4, U-238 = 2.2e-2, O-16 = 4.6e-2 }

[material.water]
temperature = 600
nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }

[cell.pin]
zones = [{ material = "fuel", radius = 0.4096, rings = 1 }, { material = "water", rings = 1 }]

[depletion]
chain = "unused.msgpack"
power = 40.0
burnup = [0.05, 0.5, 2.0]
materials = ["fuel"]
output = "unused.msgpack"
"#;

/// 入力 `text` を、テストのチェーンとライブラリで燃焼させる
fn deplete_input(text: &str) -> History {
    deplete_with(text, &library())
}

/// 入力 `text` を、テストのチェーンとライブラリ `library` で燃焼させる
pub(crate) fn deplete_with(text: &str, library: &Library) -> History {
    deplete_with_moc(text, library, &Options::default())
}

/// `deplete_with` の MOC を `moc` の打ち切り条件と加速で解く
pub(crate) fn deplete_with_moc(text: &str, library: &Library, moc: &Options) -> History {
    let input = LatticeInput::parse(text).unwrap();
    let depletion = input.depletion.clone().unwrap();
    let (lattice, region_materials) = input.geometry().unwrap();
    let tracking = Tracking::new(&lattice, input.azimuthal, input.spacing).unwrap();
    let forced = match input.symmetry {
        crate::input::SymmetryChoice::Detect => None,
        crate::input::SymmetryChoice::Force(s) => Some(s),
    };
    let chosen = choose(&tracking, &region_materials, forced).unwrap();
    run_with_moc(
        &input,
        &depletion,
        library,
        &chain(),
        &tracking,
        &region_materials,
        &chosen,
        &LibraryReactions(library),
        moc,
    )
    .unwrap()
}

#[test]
fn steps_last_the_burnup_over_the_power() {
    // 1 MWd/kg = 1000 W d/g なので、40 W/g で 0.5 MWd/kg は 12.5 日
    let steps = step_times(40.0, &[0.5, 1.5]);
    assert_eq!(steps.len(), 2);
    assert!((steps[0] - 12.5 * 86_400.0).abs() <= 1e-9 * steps[0]);
    assert!((steps[1] - 25.0 * 86_400.0).abs() <= 1e-9 * steps[1]);
}

#[test]
fn yield_cutoff_is_the_geometric_mean_of_the_thermal_and_fast_points() {
    // ADR-0183: √(0.0253 · 5e5) ≈ 112.47 eV（参照解の OpenMC に与える cutoff と同じ）
    assert_eq!(YIELD_THERMAL, 0.0253);
    assert_eq!(YIELD_FAST, 5e5);
    let expected = (0.0253f64 * 5e5).sqrt();
    assert!((yield_cutoff() - expected).abs() <= 1e-15 * expected);
}

/// 収率の点だけを持つ核分裂（重みの計算は Q と生成物を見ない）
fn yield_points(energies: &[f64]) -> Fission {
    Fission {
        q: Q,
        yield_source: "U-235".into(),
        energies: energies.to_vec(),
        yields: vec![vec![]; energies.len()],
    }
}

fn assert_weights(weights: &[f64], expected: &[f64]) {
    assert_eq!(weights.len(), expected.len(), "{weights:?} vs {expected:?}");
    for (w, e) in weights.iter().zip(expected) {
        assert!((w - e).abs() <= 1e-14, "{weights:?} vs {expected:?}");
    }
}

/// 群は 1e7〜1e3、1e3〜1、1〜1e-5 eV、核分裂率は 1、2、3
const YIELD_GROUPS: [f64; 4] = [1e7, 1e3, 1.0, 1e-5];
const YIELD_FISSION: [f64; 3] = [1.0, 2.0, 3.0];

/// cutoff より下の核分裂の割合。cutoff は 2 つ目の群をレサジーの幅の比で分ける
fn below_cutoff() -> f64 {
    (3.0 + 2.0 * yield_cutoff().ln() / 1e3f64.ln()) / 6.0
}

#[test]
fn yield_weights_split_the_group_at_the_cutoff_by_lethargy() {
    let weights = yield_weights(&YIELD_GROUPS, &YIELD_FISSION, &yield_points(&[0.0253, 5e5]));
    let below = below_cutoff();
    assert_weights(&weights, &[below, 1.0 - below]);
}

#[test]
fn yield_weights_take_the_thermal_and_fast_points_like_openmc() {
    // OpenMC 0.15.2 の FissionYieldCutoffHelper と同じ点の選び方（ADR-0183）
    let below = below_cutoff();
    let weights =
        |energies: &[f64]| yield_weights(&YIELD_GROUPS, &YIELD_FISSION, &yield_points(energies));
    // U-235 の 3 点、Pu-239 の 4 点: 熱と 500 keV の点だけ。14 MeV と 2 MeV の点は使わない
    assert_weights(&weights(&[0.0253, 5e5, 1.4e7]), &[below, 1.0 - below, 0.0]);
    assert_weights(
        &weights(&[0.0253, 5e5, 2e6, 1.4e7]),
        &[below, 1.0 - below, 0.0, 0.0],
    );
    // U-238 の 2 点（熱の点がなく、cutoff が最初の点より下）: すべての核分裂に最初の点
    assert_weights(&weights(&[5e5, 1.4e7]), &[1.0, 0.0]);
    // 1 点（cutoff が最後の点より上）: その点
    assert_weights(&weights(&[0.0253]), &[1.0]);
    // 熱の点がなければ cutoff より下で熱に最も近い点、速の点がなければ cutoff より上で速に最も近い点
    assert_weights(&weights(&[0.001, 0.03, 5e5]), &[0.0, below, 1.0 - below]);
    assert_weights(&weights(&[0.0253, 1e6, 1.4e7]), &[below, 1.0 - below, 0.0]);
}

#[test]
fn fissionless_regions_put_the_yield_weight_on_the_first_point() {
    // Chain::zero_rates と同じ（核分裂率が 0 なら燃焼行列に効かない）
    assert_weights(
        &yield_weights(&YIELD_GROUPS, &[0.0; 3], &yield_points(&[0.001, 0.03, 5e5])),
        &[1.0, 0.0, 0.0],
    );
    assert_weights(
        &yield_weights(
            &YIELD_GROUPS,
            &[0.0; 3],
            &yield_points(&[0.0253, 5e5, 1.4e7]),
        ),
        &[1.0, 0.0, 0.0],
    );
}

/// A (Xe-135 の位置) を (n,γ) で B に変える 2 核種のチェーン。反応率は A の数密度に比例し、dA/dt = −c A² の解析解を持つ
#[allow(
    clippy::type_complexity,
    reason = "the closure that counts the calls has no nameable type"
)]
fn quadratic() -> (
    Chain,
    impl FnMut(&[Vec<f64>], TransportAt) -> Result<Evaluation<(usize, TransportAt)>, ChainError>,
) {
    let mut a = stable("Xe-135", 54, 135);
    a.reactions.push(capture(1));
    let chain = Chain {
        nuclides: vec![a, stable("Xe-136", 54, 136)],
    };
    let zero = chain.zero_rates();
    let mut calls = 0;
    let solve = move |n: &[Vec<f64>], at: TransportAt| {
        calls += 1;
        let mut rates = zero.clone();
        rates[0].reactions[0][0] = 0.3 * n[0][0];
        Ok(Evaluation {
            rates: vec![rates],
            extra: (calls, at),
        })
    };
    (chain, solve)
}

#[test]
fn ce_cm_solves_the_transport_twice_per_step_and_at_the_end() {
    let (chain, mut solve) = quadratic();
    let mut midpoints = 0;
    let counting = |n: &[Vec<f64>], at: TransportAt| {
        midpoints += usize::from(at == TransportAt::Midpoint);
        solve(n, at)
    };
    let (points, _) = ce_cm(&chain, vec![vec![1.0, 0.0]], &[1.0, 1.0, 1.0], counting).unwrap();
    assert_eq!(points.len(), 4);
    assert_eq!(midpoints, 3);
    // 点の値は、その点の組成で解いた輸送（始まりの点は 1、3、5 回目、最後の点は 7 回目）で、その輸送に点の組成と伝える
    // （中点の輸送は少数群定数などを作らなくてよい）
    let calls: Vec<(usize, TransportAt)> = points.iter().map(|p| p.extra).collect();
    let point = TransportAt::Point;
    assert_eq!(calls, vec![(1, point), (3, point), (5, point), (7, point)]);
    let times: Vec<f64> = points.iter().map(|p| p.time).collect();
    assert_eq!(times, vec![0.0, 1.0, 2.0, 3.0]);
}

#[test]
fn ce_cm_with_constant_rates_is_one_matrix_exponential() {
    let chain = chain();
    let mut rates = chain.zero_rates();
    rates[U235].fission = 1e-9;
    rates[U235].reactions[0][0] = 3e-10;
    rates[XE135].reactions[0][0] = 2e-5;
    let mut n0 = vec![0.0; chain.nuclides.len()];
    n0[U235] = 7e-4;
    let solve = |_: &[Vec<f64>], _: TransportAt| -> Result<Evaluation<()>, ChainError> {
        Ok(Evaluation {
            rates: vec![rates.clone()],
            extra: (),
        })
    };
    let (points, _) = ce_cm(&chain, vec![n0.clone()], &[1e6], solve).unwrap();
    let expected = deplete(&chain.matrix(&rates).unwrap(), &n0, 1e6, CramOrder::Order48).unwrap();
    for (a, b) in points[1].densities[0].iter().zip(&expected) {
        assert!((a - b).abs() <= 1e-14 * 7e-4, "{a} vs {b}");
    }
}

#[test]
fn ce_cm_converges_at_second_order() {
    // dA/dt = −c A² の解は A(t) = 1 / (1 + c t)
    let error = |steps: usize| {
        let (chain, solve) = quadratic();
        let h = 2.0 / steps as f64;
        let (points, _) = ce_cm(&chain, vec![vec![1.0, 0.0]], &vec![h; steps], solve).unwrap();
        let a = points.last().unwrap().densities[0][0];
        (a - 1.0 / (1.0 + 0.3 * 2.0)).abs()
    };
    let (coarse, fine) = (error(4), error(8));
    let order = (coarse / fine).log2();
    assert!(
        (1.8..2.2).contains(&order),
        "order {order} ({coarse} {fine})"
    );
}

#[test]
fn region_rates_average_sigma_phi_over_the_region_volume() {
    let chain = chain();
    let library = library();
    // 2 つの領域（体積 1 と 3）の U-235 と Xe-135 の微視的断面積と中性子束
    let shielded = |name: &str| library.shielded(name, T, &[1e10, 1e10]).unwrap();
    let (u, xe) = (shielded("U-235"), shielded("Xe-135"));
    let members = [
        Member {
            volume: 1.0,
            flux: &[2.0, 1.0],
            cross_sections: &[u.clone(), xe.clone()][..],
        },
        Member {
            volume: 3.0,
            flux: &[1.0, 3.0],
            cross_sections: &[u.clone(), xe.clone()][..],
        },
    ];
    let rates = region_rates(
        &chain,
        &[U235, XE135],
        T,
        &members,
        library.groups.bounds(),
        &LibraryReactions(&library),
    )
    .unwrap();
    // 体積で平均した中性子束は (1·[2, 1] + 3·[1, 3]) / 4 = [1.25, 2.5]、反応率は σ (cm²) × φ
    let phi = [1.25, 2.5];
    let fission = (1.2 * phi[0] + 500.0 * phi[1]) * 1e-24;
    assert!((rates[U235].fission - fission).abs() <= 1e-12 * fission);
    let capture = (1.0 * phi[0] + 2.0e6 * phi[1]) * 1e-24;
    assert!((rates[XE135].reactions[0][0] - capture).abs() <= 1e-12 * capture);
    // 輸送しない核種は 0
    assert_eq!(rates[U236], chain.zero_rates()[U236]);
    // 収率の重みは、熱の点に熱の群（0.625 eV 以下）と速い群の cutoff より下の部分、500 keV の点に残り。群の核分裂率
    // （体積 × σf × φ）は速い群 1·1.2·2 + 3·1.2·1 = 6、熱の群 1·500·1 + 3·500·3 = 5000
    let below = (yield_cutoff() / 0.625).ln() / (1.0e7f64 / 0.625).ln();
    let thermal = (5000.0 + 6.0 * below) / 5006.0;
    assert_weights(&rates[U235].yield_weights, &[thermal, 1.0 - thermal]);
}

/// 2 つの状態（MF9 の LFS 0 と 2）に分ける捕獲と、(n,2n) を持つ U-235 のライブラリと、その反応のチェーンの形
fn split_capture() -> (Library, Reaction, Reaction) {
    let mut library = library();
    let u235 = library
        .nuclides
        .iter_mut()
        .find(|n| n.name == "U-235")
        .unwrap();
    let state = |lfs| {
        Some(ChannelState {
            mf: 9,
            lfs: Some(lfs),
        })
    };
    u235.temperatures[0].reactions = vec![
        GroupReaction {
            mt: 16,
            channels: vec![GroupChannel {
                final_state: None,
                cross_section: vec![0.1, 0.0],
            }],
        },
        GroupReaction {
            mt: 102,
            channels: vec![
                GroupChannel {
                    final_state: state(0),
                    cross_section: vec![1.0, 3.0],
                },
                GroupChannel {
                    final_state: state(2),
                    cross_section: vec![3.0, 1.0],
                },
            ],
        },
    ];
    let channel = |lfs| Channel {
        final_state: Some(FinalState {
            mf: 9,
            lfs: Some(lfs),
        }),
        products: vec![Product {
            nuclide: U236,
            count: 1,
        }],
    };
    let capture = Reaction {
        mt: 102,
        channels: vec![channel(0), channel(2)],
    };
    let n2n = Reaction {
        mt: 16,
        channels: vec![Channel {
            final_state: None,
            products: vec![],
        }],
    };
    (library, capture, n2n)
}

/// 自己遮蔽した捕獲 [2, 8] の群定数
fn shielded_capture() -> Shielded {
    Shielded {
        total: vec![20.0, 600.0],
        elastic: vec![10.0, 10.0],
        fission: vec![1.2, 500.0],
        nu_fission: vec![2.9, 1215.0],
        capture: vec![2.0, 8.0],
    }
}

#[test]
fn capture_channels_split_the_shielded_capture_by_the_dilute_ratios() {
    let (library, capture, _) = split_capture();
    let channels = LibraryReactions(&library)
        .cross_sections("U-235", T, &capture, &shielded_capture())
        .unwrap();
    // 群 1 は 1:3、群 2 は 3:1 で分ける
    assert_eq!(channels, vec![vec![0.5, 6.0], vec![1.5, 2.0]]);
}

#[test]
fn threshold_reactions_take_the_dilute_cross_sections() {
    let (library, _, n2n) = split_capture();
    let channels = LibraryReactions(&library)
        .cross_sections("U-235", T, &n2n, &shielded_capture())
        .unwrap();
    assert_eq!(channels, vec![vec![0.1, 0.0]]);
}

#[test]
fn library_channels_must_match_the_chain() {
    let (library, capture, _) = split_capture();
    let source = LibraryReactions(&library);
    let shielded = shielded_capture();
    // 経路の数が違う
    let one = Reaction {
        channels: capture.channels[..1].to_vec(),
        ..capture.clone()
    };
    assert!(source.cross_sections("U-235", T, &one, &shielded).is_err());
    // 状態が違う
    let mut other = capture.clone();
    other.channels[1].final_state = Some(FinalState {
        mf: 10,
        lfs: Some(2),
    });
    assert!(
        source
            .cross_sections("U-235", T, &other, &shielded)
            .is_err()
    );
    // ライブラリにない反応
    let n3n = Reaction { mt: 17, ..capture };
    let error = source
        .cross_sections("U-235", T, &n3n, &shielded)
        .unwrap_err();
    assert!(error.contains("U-235") && error.contains("MT17"), "{error}");
}

#[test]
fn reports_each_burnup_point_from_the_initial_composition() {
    let history = deplete_input(PIN);
    assert_eq!(history.points.len(), 4);
    let burnup: Vec<f64> = history.points.iter().map(|p| p.burnup).collect();
    assert_eq!(burnup, vec![0.0, 0.05, 0.5, 2.0]);
    let times: Vec<f64> = history.points.iter().map(|p| p.time).collect();
    assert!(
        (times[3] - 50.0 * 86_400.0).abs() <= 1e-9 * times[3],
        "{times:?}"
    );
    assert_eq!(history.regions.len(), 1);
    assert_eq!(history.regions[0].material, "fuel");
    let initial = &history.points[0].densities[0];
    assert_eq!(initial[U235], 7.0e-4);
    assert_eq!(initial[U238], 2.2e-2);
    assert_eq!(initial[XE135], 0.0);
    // Xe-135 が生まれ、k は始まりより下がる
    assert!(history.points[1].densities[0][XE135] > 0.0);
    assert!(history.points[3].k < history.points[0].k);
    assert!(history.points.iter().all(|p| p.few_group.is_none()));
}

#[test]
fn every_transportable_nuclide_gets_rates_from_the_first_transport() {
    // OpenMC は数密度が 0 の核種にも 1 原子あたりの反応率を与える（燃焼の tally は数密度を掛けない）。始まりに 0 の Xe-135 も
    // 最初の輸送から反応率を持ち、予測子で燃え尽きる（ADR-0173）
    let history = deplete_input(PIN);
    assert!(history.points.iter().all(|p| p.transported == 5));
}

#[test]
fn transport_needs_data_at_exactly_the_region_temperature() {
    // 温度の内挿ができる核種（熱中性子散乱）でも、領域の温度のデータがなければ輸送しない（OpenMC の参照解と同じ集合）
    let mut library = library();
    let mut h1 = library
        .nuclides
        .iter()
        .find(|n| n.name == "H-1")
        .unwrap()
        .clone();
    h1.name = "H-2".into();
    h1.thermal_scattering = true;
    let mut other = h1.temperatures[0].clone();
    h1.temperatures[0].temperature = 550.0;
    other.temperature = 650.0;
    h1.temperatures.push(other);
    library.nuclides.push(h1);
    assert!(library.lookup("H-2", T).is_ok());
    assert!(!has_data_at(&library, "H-2", T));
    assert!(has_data_at(&library, "H-2", 650.0));
    assert!(has_data_at(&library, "U-235", T));
    assert!(!has_data_at(&library, "Ba-140", T));
}

#[test]
fn fission_products_count_the_energy_released() {
    // 核分裂生成物の原子は核分裂の 2 倍。核分裂の数 × Q が、出力密度 × 重金属の質量 × 時間に等しい
    let history = deplete_input(PIN);
    let region = &history.regions[0];
    let last = history.points.last().unwrap();
    let products: f64 = [XE135, XE136, CS135, BA140]
        .iter()
        .map(|&k| last.densities[0][k])
        .sum();
    let fissions = 0.5 * products * 1e24 * region.volume;
    let energy = fissions * Q * 1.602_176_634e-19;
    let mass = (7.0e-4 * 233.025 + 2.2e-2 * 236.006)
        * 1.008_664_915_95
        * 1.660_539_066_60e-24
        * 1e24
        * region.volume;
    let expected = 40.0 * mass * last.time;
    assert!(
        (energy - expected).abs() <= 1e-3 * expected,
        "{energy} vs {expected}"
    );
}

#[test]
fn few_group_constants_are_reported_at_each_point() {
    let text = format!("{PIN}few-group = {{ bounds = [1e7, 0.625, 1e-5] }}\n");
    let history = deplete_input(&text);
    for point in &history.points {
        let few = point.few_group.as_ref().unwrap();
        assert!(
            (few.k_inf - point.k).abs() <= 1e-6,
            "{} vs {}",
            few.k_inf,
            point.k
        );
    }
}

#[test]
fn pins_of_one_symmetry_orbit_deplete_as_one_region() {
    let lattice = PIN.replace("rows = [\"pin\"]", "rows = [\"pin pin\", \"pin pin\"]");
    let detected = deplete_input(&lattice);
    assert_eq!(detected.regions.len(), 1);
    assert_eq!(detected.regions[0].regions.len(), 4);
    let none =
        deplete_input(&lattice.replace("spacing = 0.05", "spacing = 0.05\nsymmetry = \"none\""));
    assert_eq!(none.regions.len(), 4);
    // 移り合うピンは同じように燃える
    let xe = |h: &History, r: usize| h.points[3].densities[r][XE135];
    for r in 0..4 {
        assert!((xe(&none, r) - xe(&detected, 0)).abs() <= 1e-6 * xe(&detected, 0));
    }
}

#[test]
fn fuel_rings_deplete_as_their_own_regions() {
    let text = PIN
        .replace(
            "zones = [{ material = \"fuel\", radius = 0.4096, rings = 1 }",
            "zones = [{ material = \"inner\", radius = 0.3, rings = 1 }, { material = \"outer\", radius = 0.4096, rings = 1 }",
        )
        .replace("[material.fuel]", "[material.inner]")
        .replace(
            "[material.water]",
            "[material.outer]\ntemperature = 600\nnuclides = { U-235 = 7.0e-4, U-238 = 2.2e-2, O-16 = 4.6e-2 }\n\n[material.water]",
        )
        .replace("materials = [\"fuel\"]", "materials = [\"inner\", \"outer\"]");
    let history = deplete_input(&text);
    let names: Vec<&str> = history
        .regions
        .iter()
        .map(|r| r.material.as_str())
        .collect();
    assert_eq!(names, vec!["inner", "outer"]);
    // 外側のリングは熱中性子束が高く、U-235 が速く減る
    let u235 = |r: usize| history.points[3].densities[r][U235];
    assert!(u235(1) < u235(0), "{} {}", u235(0), u235(1));
}

#[test]
fn history_is_written_as_a_map_with_the_format_name() {
    let history = deplete_input(PIN);
    let bytes = history.to_msgpack().unwrap();
    let value: rmpv::Value = rmp_serde::from_slice(&bytes).unwrap();
    let map = value.as_map().unwrap();
    let get = |key: &str| {
        map.iter()
            .find(|(k, _)| k.as_str() == Some(key))
            .map(|(_, v)| v)
    };
    assert_eq!(
        get("format").and_then(|v| v.as_str()),
        Some("xlc-depletion-history")
    );
    assert_eq!(
        get("nuclides").and_then(|v| v.as_array()).map(Vec::len),
        Some(9)
    );
    assert_eq!(
        get("points").and_then(|v| v.as_array()).map(Vec::len),
        Some(4)
    );
    // 時間の内訳（秒）と輸送の回数（記録の時間の表に使う）
    let seconds = get("seconds").and_then(|v| v.as_map()).unwrap();
    let keys: Vec<&str> = seconds.iter().filter_map(|(k, _)| k.as_str()).collect();
    assert_eq!(keys, ["shielding", "moc", "rates", "matrix", "cram"]);
    assert_eq!(get("solves").and_then(rmpv::Value::as_u64), Some(7));
    assert_eq!(get("shieldings").and_then(rmpv::Value::as_u64), Some(7));
}

/// 履歴のファイルの、各点の `pin_powers`（[点][行][列]）
fn written_pin_powers(history: &History) -> Vec<Vec<Vec<f64>>> {
    let bytes = history.to_msgpack().unwrap();
    let value: rmpv::Value = rmp_serde::from_slice(&bytes).unwrap();
    let field = |map: &rmpv::Value, key: &str| -> rmpv::Value {
        map.as_map()
            .unwrap()
            .iter()
            .find(|(k, _)| k.as_str() == Some(key))
            .map(|(_, v)| v.clone())
            .unwrap_or_else(|| panic!("no `{key}`"))
    };
    field(&value, "points")
        .as_array()
        .unwrap()
        .iter()
        .map(|p| {
            field(p, "pin_powers")
                .as_array()
                .unwrap()
                .iter()
                .map(|row| {
                    row.as_array()
                        .unwrap()
                        .iter()
                        .map(|x| x.as_f64().unwrap())
                        .collect()
                })
                .collect()
        })
        .collect()
}

/// 2×2 の格子の右下を水のセルにした入力（対角線の鏡映で、右上と左下のピンが移り合う）
fn with_water_cell(text: &str) -> String {
    format!(
        "{}\n[cell.box]\nzones = [{{ material = \"water\", rings = 1 }}]\n",
        text.replace("rows = [\"pin\"]", "rows = [\"pin pin\", \"pin box\"]")
    )
}

#[test]
fn pin_powers_are_reported_at_each_point() {
    // ピン出力は lattice の `pin-powers` と同じ: セルごとの核分裂率を、核分裂のあるセルの平均で 1 にしたもの（ADR-0181）
    let history = deplete_input(&with_water_cell(PIN));
    let powers = written_pin_powers(&history);
    assert_eq!(powers.len(), history.points.len());
    for rows in &powers {
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().all(|r| r.len() == 2));
        assert_eq!(rows[1][1], 0.0);
        let mean = (rows[0][0] + rows[0][1] + rows[1][0]) / 3.0;
        assert!((mean - 1.0).abs() <= 1e-12, "{rows:?}");
        assert!((rows[0][1] - rows[1][0]).abs() <= 1e-9, "{rows:?}");
        // 水に面したピンは熱中性子束が高く、出力が高い
        assert!(rows[0][1] > rows[0][0], "{rows:?}");
    }
}

#[test]
fn the_pin_powers_of_a_point_come_from_its_own_composition() {
    // 点の値は、その点の組成で解いた輸送のもの（予測子の始まり）。中点の輸送の値なら、初めの点の値が最初のステップの長さで変わる
    let short = deplete_input(&with_water_cell(PIN));
    let long = deplete_input(&with_water_cell(
        &PIN.replace("burnup = [0.05, 0.5, 2.0]", "burnup = [5.0]"),
    ));
    let (short, long) = (written_pin_powers(&short), written_pin_powers(&long));
    // MOC は並列に解くので、ビットの一致でなく丸めの水準で比べる
    for (a, b) in short[0].iter().flatten().zip(long[0].iter().flatten()) {
        assert!((a - b).abs() <= 1e-12, "{a} vs {b}");
    }
    // 燃焼で組成が変われば、ピン出力も変わる
    assert_ne!(short[0], short[3]);
}

#[test]
fn a_half_life_without_reactions_halves_the_density() {
    let chain = chain();
    let rates = chain.zero_rates();
    let mut n0 = vec![0.0; chain.nuclides.len()];
    n0[XE135] = 1.0;
    let solve = |_: &[Vec<f64>], _: TransportAt| -> Result<Evaluation<()>, ChainError> {
        Ok(Evaluation {
            rates: vec![rates.clone()],
            extra: (),
        })
    };
    let (points, _) = ce_cm(&chain, vec![n0], &[32_904.0], solve).unwrap();
    assert!((points[1].densities[0][XE135] - 0.5).abs() <= 1e-12);
}

#[test]
fn nuclides_without_data_at_the_fuel_temperature_only_decay() {
    // 燃料の温度（600 K）のデータのない核種は輸送せず、崩壊と生成だけを追う（ADR-0173）
    let mut library = library();
    let mut ba140 = library
        .nuclides
        .iter()
        .find(|n| n.name == "O-16")
        .unwrap()
        .clone();
    ba140.name = "Ba-140".into();
    ba140.temperatures[0].temperature = 300.0;
    library.nuclides.push(ba140);
    let history = deplete_with(PIN, &library);
    let last = history.points.last().unwrap();
    assert!(last.densities[0][BA140] > 0.0);
    // 輸送できる核種の集合（記録で OpenMC と比べる）
    assert_eq!(
        history.transportable,
        ["H-1", "O-16", "Xe-135", "U-235", "U-238"].map(String::from)
    );
    // 輸送しないので、Ba-140 があっても k は同じ
    let without = deplete_input(PIN);
    assert_eq!(last.k, without.points.last().unwrap().k);
}

/// 経路の断面積を 1 群だけ返す受け口（群の数の誤り）
struct ShortChannels;

impl ReactionSource for ShortChannels {
    fn cross_sections(
        &self,
        _nuclide: &str,
        _temperature: f64,
        reaction: &Reaction,
        _shielded: &Shielded,
    ) -> Result<Vec<Vec<f64>>, String> {
        Ok(vec![vec![1.0]; reaction.channels.len()])
    }
}

#[test]
fn channel_cross_sections_need_every_group() {
    let chain = chain();
    let library = library();
    let u = library.shielded("U-235", T, &[1e10, 1e10]).unwrap();
    let members = [Member {
        volume: 1.0,
        flux: &[1.0, 1.0],
        cross_sections: &[u][..],
    }];
    let result = region_rates(
        &chain,
        &[U235],
        T,
        &members,
        library.groups.bounds(),
        &ShortChannels,
    );
    assert!(result.is_err());
}

/// 燃料のピンの中の薄い被覆 `coat`（IFBA の型）と、別のセルの吸収体の棒 `poison`（Pyrex の型）も燃焼させる PIN の格子
fn with_absorbers() -> String {
    PIN
        .replace("rows = [\"pin\"]", "rows = [\"pin rod\"]")
        .replace(
            "zones = [{ material = \"fuel\", radius = 0.4096, rings = 1 }",
            "zones = [{ material = \"fuel\", radius = 0.4, rings = 1 }, { material = \"coat\", radius = 0.4096, rings = 1 }",
        )
        .replace(
            "[material.water]",
            "[material.coat]\ntemperature = 600\nnuclides = { Xe-135 = 1.0e-4, O-16 = 4.6e-2 }\n\n\
             [material.poison]\ntemperature = 600\nnuclides = { Xe-135 = 1.0e-4, O-16 = 4.6e-2 }\n\n[material.water]",
        )
        .replace(
            "[depletion]",
            "[cell.rod]\nzones = [{ material = \"poison\", radius = 0.4, rings = 1 }, { material = \"water\", rings = 1 }]\n\n[depletion]",
        )
        .replace("materials = [\"fuel\"]", "materials = [\"fuel\", \"coat\", \"poison\"]")
}

#[test]
fn absorbers_outside_the_fuel_deplete_by_their_own_rates() {
    // どちらの吸収体も Xe-135 で、崩壊だけのときより速く減る（捕獲の反応率を持つ）。燃焼の領域は、対称操作の組と物質の対（ADR-0173）
    let history = deplete_input(&with_absorbers());
    let names: Vec<&str> = history
        .regions
        .iter()
        .map(|r| r.material.as_str())
        .collect();
    assert_eq!(names, vec!["fuel", "coat", "poison"]);
    let last = history.points.last().unwrap();
    let decay_only = 1.0e-4 * (-LN_2 * last.time / 32_904.0).exp();
    for (name, densities) in names.iter().zip(&last.densities).skip(1) {
        let xe = densities[XE135];
        assert!(xe < 0.5 * decay_only, "{name}: {xe} vs {decay_only}");
    }
}

#[test]
fn depleted_materials_share_one_temperature() {
    // 輸送できる核種は燃焼の領域の温度で決まる。吸収体を燃料と違う温度で燃焼させるときの集合は決めていないので止める
    let text = PIN
        .replace(
            "zones = [{ material = \"fuel\", radius = 0.4096, rings = 1 }",
            "zones = [{ material = \"fuel\", radius = 0.4, rings = 1 }, { material = \"coat\", radius = 0.4096, rings = 1 }",
        )
        .replace(
            "[material.water]",
            "[material.coat]\ntemperature = 300\nnuclides = { O-16 = 4.6e-2 }\n\n[material.water]",
        )
        .replace("materials = [\"fuel\"]", "materials = [\"fuel\", \"coat\"]");
    let input = LatticeInput::parse(&text).unwrap();
    let depletion = input.depletion.clone().unwrap();
    let (lattice, region_materials) = input.geometry().unwrap();
    let tracking = Tracking::new(&lattice, input.azimuthal, input.spacing).unwrap();
    let chosen = choose(&tracking, &region_materials, None).unwrap();
    let library = library();
    let error = run(
        &input,
        &depletion,
        &library,
        &chain(),
        &tracking,
        &region_materials,
        &chosen,
        &LibraryReactions(&library),
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("several temperatures"), "{error}");
}

#[test]
fn history_round_trips_through_its_file() {
    let text = format!("{PIN}few-group = {{ bounds = [1e7, 0.625, 1e-5] }}\n");
    let history = deplete_input(&text);
    let read = History::from_msgpack(&history.to_msgpack().unwrap()).unwrap();
    assert_eq!(read.nuclides, history.nuclides);
    assert_eq!(read.regions, history.regions);
    assert_eq!(read.transportable, history.transportable);
    assert_eq!(read.solves, history.solves);
    assert_eq!(read.shieldings, history.shieldings);
    assert_eq!(read.points.len(), history.points.len());
    for (r, h) in read.points.iter().zip(&history.points) {
        assert_eq!(r.time.to_bits(), h.time.to_bits());
        assert_eq!(r.burnup.to_bits(), h.burnup.to_bits());
        assert_eq!(r.k.to_bits(), h.k.to_bits());
        assert_eq!(r.densities, h.densities);
        assert_eq!(r.pin_powers, h.pin_powers);
        assert_eq!(r.few_group, h.few_group);
    }
}

#[test]
fn reading_a_history_checks_the_format_name() {
    let mut bytes = deplete_input(PIN).to_msgpack().unwrap();
    let at = bytes
        .windows(21)
        .position(|w| w == b"xlc-depletion-history")
        .unwrap();
    bytes[at + 4] = b'D';
    assert!(History::from_msgpack(&bytes).is_err());
}

#[test]
fn histories_without_depletion_regions_are_errors() {
    let bytes = deplete_input(PIN).to_msgpack().unwrap();
    let value: rmpv::Value = rmp_serde::from_slice(&bytes).unwrap();
    // 燃焼の領域と、各点の数密度（[領域][核種]）を空にする
    let emptied = |v: rmpv::Value, key: &str| match v {
        rmpv::Value::Map(entries) => rmpv::Value::Map(
            entries
                .into_iter()
                .map(|(k, v)| {
                    if k.as_str() == Some(key) {
                        (k, rmpv::Value::Array(Vec::new()))
                    } else {
                        (k, v)
                    }
                })
                .collect(),
        ),
        other => other,
    };
    let value = match emptied(value, "regions") {
        rmpv::Value::Map(entries) => rmpv::Value::Map(
            entries
                .into_iter()
                .map(|(k, v)| match (k.as_str(), v) {
                    (Some("points"), rmpv::Value::Array(points)) => (
                        k,
                        rmpv::Value::Array(
                            points
                                .into_iter()
                                .map(|p| emptied(p, "densities"))
                                .collect(),
                        ),
                    ),
                    (_, v) => (k, v),
                })
                .collect(),
        ),
        other => other,
    };
    let mut empty = Vec::new();
    rmpv::encode::write_value(&mut empty, &value).unwrap();
    let err = History::from_msgpack(&empty).unwrap_err();
    assert!(format!("{err:#}").contains("regions"), "{err:#}");
}

#[test]
fn histories_with_repeated_nuclides_are_errors() {
    let mut history = deplete_input(PIN);
    // 核種の名前で数密度の列を引くので、同じ名前が 2 つあれば列を選べない
    history.nuclides[1] = history.nuclides[0].clone();
    let err = History::from_msgpack(&history.to_msgpack().unwrap()).unwrap_err();
    assert!(format!("{err:#}").contains("twice"), "{err:#}");
}

#[test]
fn every_transport_solves_the_shielding_by_default() {
    // 自己遮蔽の頻度を書かない入力は、ステップの始まりと中点の両方で解き直す（ADR-0173 の案 F1）
    let history = deplete_input(PIN);
    assert_eq!(history.solves, 7);
    assert_eq!(history.shieldings, 7);
}

#[test]
fn step_start_shielding_skips_the_midpoint_transports() {
    let every = deplete_input(PIN);
    let start = deplete_input(&format!("{PIN}shielding = \"step-start\"\n"));
    // 輸送はステップに 2 回と最後の点で 1 回のまま、自己遮蔽は始まりと最後の点だけ
    assert_eq!(start.solves, 7);
    assert_eq!(start.shieldings, 4);
    // テストのライブラリの群定数は σ0 によらないので、中点の微視的断面積は組成を変えて解き直しても同じになり、中点の
    // 物質は始まりの微視的断面積と中点の数密度で作り直したものと、丸めの誤差の範囲で等しい
    let close = |a: f64, b: f64| (a - b).abs() <= 1e-12 * a.abs().max(b.abs());
    for (s, e) in start.points.iter().zip(&every.points) {
        assert!(close(s.k, e.k), "{} vs {}", s.k, e.k);
        for (ns, ne) in s.densities[0].iter().zip(&e.densities[0]) {
            assert!(close(*ns, *ne), "{ns} vs {ne}");
        }
    }
}

/// PIN の格子を、入力の物質のまま自己遮蔽した結果と、燃料の燃焼の領域と、入力の物質（燃料と水）
fn shielded_pin() -> (
    Library,
    Lattice,
    FuelShielding,
    Vec<DepletionRegion>,
    Vec<MaterialInput>,
) {
    let library = library();
    let input = LatticeInput::parse(PIN).unwrap();
    let (lattice, region_materials) = input.geometry().unwrap();
    let tracking = Tracking::new(&lattice, input.azimuthal, input.spacing).unwrap();
    let chosen = choose(&tracking, &region_materials, None).unwrap();
    let representatives = chosen.shielding().cell_representatives();
    let regions = depletion_regions(
        &lattice,
        &region_materials,
        &representatives,
        &["fuel".to_string()],
    );
    let shielded = shield_fuel_with(
        &library,
        &tracking,
        &input.materials,
        &region_materials,
        chosen.shielding(),
        &input.shielding_options(),
    )
    .unwrap();
    (library, lattice, shielded, regions, input.materials)
}

/// `shielded` の燃料の燃焼の領域の物質を、入力の物質 `materials` の燃料の組成で作り直したもの
fn remixed_pin(
    library: &Library,
    lattice: &Lattice,
    shielded: &FuelShielding,
    regions: &[DepletionRegion],
    materials: &[MaterialInput],
) -> anyhow::Result<FuelShielding> {
    let fuel = materials.iter().position(|m| m.name == "fuel").unwrap();
    let mut remixed = shielded.clone();
    for (m, material) in remix_depleted(library, shielded, lattice, regions, materials, &[fuel])? {
        remixed.materials[m] = material;
    }
    Ok(remixed)
}

#[test]
fn remixing_the_shielded_composition_reproduces_the_shielded_materials() {
    let (library, lattice, shielded, regions, materials) = shielded_pin();
    let remixed = remixed_pin(&library, &lattice, &shielded, &regions, &materials).unwrap();
    assert_eq!(remixed, shielded);
}

#[test]
fn remixed_materials_take_the_new_densities_and_keep_the_microscopic_cross_sections() {
    let (library, lattice, shielded, regions, mut materials) = shielded_pin();
    let fuel = materials.iter_mut().find(|m| m.name == "fuel").unwrap();
    for (_, density) in &mut fuel.composition {
        *density *= 2.0;
    }
    let remixed = remixed_pin(&library, &lattice, &shielded, &regions, &materials).unwrap();
    for r in 0..lattice.regions() {
        let before = &shielded.materials[shielded.region_material[r]];
        let after = &remixed.materials[remixed.region_material[r]];
        assert_eq!(
            remixed.region_cross_sections(r),
            shielded.region_cross_sections(r)
        );
        if regions[0].regions.contains(&r) {
            for (a, b) in after.total.iter().zip(&before.total) {
                assert!((a - 2.0 * b).abs() <= 1e-14 * a, "{a} vs 2 x {b}");
            }
        } else {
            assert_eq!(after, before);
        }
    }
}

#[test]
fn remixing_needs_the_shielded_nuclides_of_each_region() {
    let (library, lattice, shielded, regions, mut materials) = shielded_pin();
    let fuel = materials.iter_mut().find(|m| m.name == "fuel").unwrap();
    fuel.composition.push(("Xe-135".into(), 1e-9));
    let error = remixed_pin(&library, &lattice, &shielded, &regions, &materials).unwrap_err();
    assert!(format!("{error:#}").contains("nuclides"), "{error:#}");
}

#[test]
fn step_start_shielding_remixes_depleted_absorbers_too() {
    // 自己遮蔽しない吸収体（無限希釈の物質）も、中点で始まりの微視的断面積のまま、自己遮蔽と同じ作り方で作り直す。O-16 の弾性散乱を
    // 散乱行列の行の和と違う値にし、自己遮蔽した群定数の散乱行列の補正 (Σs + σe − σe) / Σs が丸めで 1 にならないライブラリでも、
    // 無限希釈の物質と同じになる。群定数は σ0 によらないので、輸送のたびに解いたときと丸めの誤差の範囲で等しい
    let mut library = library();
    let oxygen = library
        .nuclides
        .iter_mut()
        .find(|n| n.name == "O-16")
        .unwrap();
    for table in &mut oxygen.temperatures[0].shielded {
        table.elastic = vec![1.1, 1.1];
    }
    let every = deplete_with(&with_absorbers(), &library);
    let start = deplete_with(
        &format!("{}shielding = \"step-start\"\n", with_absorbers()),
        &library,
    );
    assert_eq!(start.shieldings, 4);
    let close = |a: f64, b: f64| (a - b).abs() <= 1e-12 * a.abs().max(b.abs());
    for (s, e) in start.points.iter().zip(&every.points) {
        assert!(close(s.k, e.k), "{} vs {}", s.k, e.k);
        for (ns, ne) in s
            .densities
            .iter()
            .flatten()
            .zip(e.densities.iter().flatten())
        {
            assert!(close(*ns, *ne), "{ns} vs {ne}");
        }
    }
}

#[test]
fn step_start_shielding_keeps_the_start_cross_sections_at_the_midpoint() {
    // U-238 の捕獲を σ0 で変え、自己遮蔽が組成に依るライブラリでは、中点で解き直す F1 と始まりの微視的断面積を使う F2 は、最初の点
    // （同じ組成の同じ輸送）だけが同じで、ステップを進めた点は違う
    let mut library = library();
    let uranium = library
        .nuclides
        .iter_mut()
        .find(|n| n.name == "U-238")
        .unwrap();
    for (table, factor) in uranium.temperatures[0]
        .shielded
        .iter_mut()
        .zip([1.0, 0.8, 0.5])
    {
        let capture = table.capture[1] * factor;
        table.total[1] += capture - table.capture[1];
        table.capture[1] = capture;
    }
    let every = deplete_with(PIN, &library);
    let start = deplete_with(&format!("{PIN}shielding = \"step-start\"\n"), &library);
    assert_eq!(start.points[0].k.to_bits(), every.points[0].k.to_bits());
    for (s, e) in start.points.iter().zip(&every.points).skip(1) {
        assert_ne!(s.k, e.k, "at {} MWd/kg", s.burnup);
    }
}
