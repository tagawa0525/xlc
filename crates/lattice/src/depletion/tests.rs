use xlc_depletion::chain::{
    Chain, Channel, Decay, DecayEnergy, DecayMode, Fission, Nuclide, Product, Reaction, Yield,
};
use xlc_depletion::{CramOrder, deplete};
use xlc_format::group::GroupStructure;
use xlc_format::library::{Library, NuclideData, Pointwise, Shielded, TemperatureData};

use super::*;
use crate::geometry::{Lattice, PinCell};
use crate::input::LatticeInput;
use crate::symmetry::choose;
use crate::tracking::Tracking;

const T: f64 = 600.0;
const XE135: usize = 2;
const XE136: usize = 3;
const CS135: usize = 4;
const BA140: usize = 5;
const U235: usize = 6;
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
fn chain() -> Chain {
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

/// 2 群（0.625 eV で分ける）の、分離共鳴のないライブラリ。全断面積は散乱、核分裂、捕獲の和
fn library() -> Library {
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
            };
            NuclideData {
                name: name.into(),
                za: 0.0,
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
                    scatter: vec![scatter.iter().map(|row| row.to_vec()).collect()],
                    shielded: vec![table.clone(), table.clone(), table],
                    pointwise: Pointwise::default(),
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
                [0.05, 0.0],
                [0.4, 2.7],
            ),
        ],
    }
}

/// 捕獲を、ライブラリの全断面積から散乱と核分裂を引いて与える受け口（テストのライブラリは全断面積がその和）
struct Captures<'a>(&'a Library);

impl ReactionSource for Captures<'_> {
    fn cross_sections(
        &self,
        nuclide: &str,
        temperature: f64,
        reaction: &Reaction,
        shielded: &Shielded,
    ) -> Result<Vec<Vec<f64>>, String> {
        if reaction.mt != 102 {
            return Err(format!("{nuclide}: no MT{}", reaction.mt));
        }
        let data = self
            .0
            .nuclide(nuclide, temperature)
            .map_err(|e| e.to_string())?;
        let capture = (0..shielded.total.len())
            .map(|g| shielded.total[g] - data.scatter[g].iter().sum::<f64>() - shielded.fission[g])
            .collect();
        Ok(vec![capture])
    }
}

const PIN: &str = r#"
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
    let input = LatticeInput::parse(text).unwrap();
    let depletion = input.depletion.clone().unwrap();
    let (lattice, region_materials) = input.geometry().unwrap();
    let tracking = Tracking::new(&lattice, input.azimuthal, input.spacing).unwrap();
    let forced = match input.symmetry {
        crate::input::SymmetryChoice::Detect => None,
        crate::input::SymmetryChoice::Force(s) => Some(s),
    };
    let chosen = choose(&tracking, &region_materials, forced).unwrap();
    let library = library();
    run(
        &input,
        &depletion,
        &library,
        &chain(),
        &tracking,
        &region_materials,
        &chosen,
        &Captures(&library),
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
fn yield_weights_split_groups_at_the_bounds_by_lethargy() {
    // 群は 1e7〜1e3、1e3〜1、1〜1e-5 eV。境界 112.5 eV は 2 つ目の群をレサジーで分ける
    let weights = yield_weights(&[1e7, 1e3, 1.0, 1e-5], &[1.0, 2.0, 3.0], &[112.5]);
    let below = 2.0 * 112.5f64.ln() / 1e3f64.ln();
    let expected = [(3.0 + below) / 6.0, (1.0 + 2.0 - below) / 6.0];
    assert_eq!(weights.len(), 2);
    for (w, e) in weights.iter().zip(expected) {
        assert!((w - e).abs() <= 1e-14, "{weights:?} vs {expected:?}");
    }
    // 1 点の収率は重み 1。核分裂のない領域は最初の点に 1 を置く（Chain::zero_rates と同じ）
    assert_eq!(
        yield_weights(&[1e7, 1e3, 1e-5], &[1.0, 2.0], &[]),
        vec![1.0]
    );
    assert_eq!(
        yield_weights(&[1e7, 1e3, 1e-5], &[0.0, 0.0], &[112.5]),
        vec![1.0, 0.0]
    );
}

/// A (Xe-135 の位置) を (n,γ) で B に変える 2 核種のチェーン。反応率は A の数密度に比例し、dA/dt = −c A² の解析解を持つ
fn quadratic() -> (
    Chain,
    impl FnMut(&[Vec<f64>]) -> Result<Evaluation<usize>, ChainError>,
) {
    let mut a = stable("Xe-135", 54, 135);
    a.reactions.push(capture(1));
    let chain = Chain {
        nuclides: vec![a, stable("Xe-136", 54, 136)],
    };
    let zero = chain.zero_rates();
    let mut calls = 0;
    let solve = move |n: &[Vec<f64>]| {
        calls += 1;
        let mut rates = zero.clone();
        rates[0].reactions[0][0] = 0.3 * n[0][0];
        Ok(Evaluation {
            rates: vec![rates],
            extra: calls,
        })
    };
    (chain, solve)
}

#[test]
fn ce_cm_solves_the_transport_twice_per_step_and_at_the_end() {
    let (chain, solve) = quadratic();
    let (points, _) = ce_cm(&chain, vec![vec![1.0, 0.0]], &[1.0, 1.0, 1.0], solve).unwrap();
    assert_eq!(points.len(), 4);
    // 点の値は、その点の組成で解いた輸送（始まりの点は 1、3、5 回目、最後の点は 7 回目）
    let calls: Vec<usize> = points.iter().map(|p| p.extra).collect();
    assert_eq!(calls, vec![1, 3, 5, 7]);
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
    let solve = |_: &[Vec<f64>]| -> Result<Evaluation<()>, ChainError> {
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
        &Captures(&library),
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
    // 収率の重みは、熱の群（0.625 eV 以下）と、境界 112.5 eV より下の速い群の部分
    let weights = &rates[U235].yield_weights;
    assert!(
        (weights.iter().sum::<f64>() - 1.0).abs() <= 1e-14,
        "{weights:?}"
    );
}

#[test]
fn the_library_has_no_reaction_cross_sections_yet() {
    let library = library();
    let shielded = library.shielded("U-235", T, &[1e10, 1e10]).unwrap();
    let error = LibraryReactions
        .cross_sections("U-235", T, &capture(U236), &shielded)
        .unwrap_err();
    assert!(
        error.contains("U-235") && error.contains("MT102"),
        "{error}"
    );
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
    let text = PIN.replace(
        "rows = [\"pin\"]",
        "rows = [\"pin\"]\nfew-group = { output = \"unused.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
    );
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
}

#[test]
fn a_half_life_without_reactions_halves_the_density() {
    let chain = chain();
    let rates = chain.zero_rates();
    let mut n0 = vec![0.0; chain.nuclides.len()];
    n0[XE135] = 1.0;
    let solve = |_: &[Vec<f64>]| -> Result<Evaluation<()>, ChainError> {
        Ok(Evaluation {
            rates: vec![rates.clone()],
            extra: (),
        })
    };
    let (points, _) = ce_cm(&chain, vec![n0], &[32_904.0], solve).unwrap();
    assert!((points[1].densities[0][XE135] - 0.5).abs() <= 1e-12);
}
