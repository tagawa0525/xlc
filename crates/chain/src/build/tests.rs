use xlc_depletion::chain::{FinalState, Nuclide, Product};

use super::*;
use crate::endf::{DecayData, DecayMode, FinalStateEntry, Material, RawYield, Step};

fn material(z: u32, a: u32, state: u32, excitation: f64) -> Material {
    Material {
        mat: 0,
        z,
        a,
        state,
        excitation,
    }
}

fn stable(z: u32, a: u32) -> DecayEvaluation {
    DecayEvaluation {
        material: material(z, a, 0, 0.0),
        decay: DecayData {
            stable: true,
            half_life: 0.0,
            energy: [0.0; 3],
            modes: vec![],
        },
    }
}

fn unstable(z: u32, a: u32, state: u32, half_life: f64, modes: Vec<DecayMode>) -> DecayEvaluation {
    DecayEvaluation {
        material: material(z, a, state, if state > 0 { 1e5 } else { 0.0 }),
        decay: DecayData {
            stable: false,
            half_life,
            energy: [1.0, 2.0, 3.0],
            modes,
        },
    }
}

fn mode(steps: &[Step], final_state: u32, ratio: f64) -> DecayMode {
    DecayMode {
        steps: steps.to_vec(),
        final_state,
        ratio,
    }
}

fn neutron(z: u32, a: u32, state: u32, mts: &[u32]) -> NeutronEvaluation {
    NeutronEvaluation {
        material: material(z, a, state, 0.0),
        mts: mts.to_vec(),
        final_states: vec![],
    }
}

fn find<'a>(built: &'a Built, name: &str) -> &'a Nuclide {
    built
        .chain
        .nuclides
        .iter()
        .find(|n| n.name == name)
        .unwrap_or_else(|| panic!("{name} is not in the chain"))
}

fn index(built: &Built, name: &str) -> usize {
    built
        .chain
        .nuclides
        .iter()
        .position(|n| n.name == name)
        .unwrap_or_else(|| panic!("{name} is not in the chain"))
}

fn one(built: &Built, name: &str) -> Product {
    Product {
        nuclide: index(built, name),
        count: 1,
    }
}

/// H-1、He-4 と、名前で与えた崩壊の評価
fn with_light(mut decays: Vec<DecayEvaluation>) -> Vec<DecayEvaluation> {
    decays.push(stable(1, 1));
    decays.push(stable(2, 4));
    decays
}

#[test]
fn orders_nuclides_by_z_a_and_state_and_names_them() {
    let neutron_itself = stable(0, 1);
    let decays = with_light(vec![
        unstable(
            95,
            242,
            1,
            4.4e9,
            vec![mode(&[Step::IsomericTransition], 0, 1.0)],
        ),
        stable(95, 242),
        neutron_itself,
    ]);
    let built = build(&decays, &[], &[]).unwrap();
    let names: Vec<_> = built
        .chain
        .nuclides
        .iter()
        .map(|n| n.name.as_str())
        .collect();
    assert_eq!(names, ["H-1", "He-4", "Am-242", "Am-242m"]);
    assert_eq!(find(&built, "Am-242m").excitation, 1e5);
    built.chain.check().unwrap();
}

#[test]
fn applies_each_decay_step_and_emits_light_nuclides() {
    let decays = with_light(vec![
        unstable(
            55,
            140,
            0,
            60.0,
            vec![
                mode(&[Step::BetaMinus], 0, 0.5),
                mode(&[Step::BetaMinus, Step::Neutron], 1, 0.25),
                mode(&[Step::EcBetaPlus, Step::Proton], 0, 0.125),
                mode(&[Step::Proton, Step::Proton], 0, 0.0625),
                mode(&[Step::Alpha], 0, 0.0625),
            ],
        ),
        stable(56, 140),
        unstable(
            56,
            139,
            1,
            1.0,
            vec![mode(&[Step::IsomericTransition], 0, 1.0)],
        ),
        stable(56, 139),
        stable(53, 139),
        stable(53, 138),
        stable(53, 136),
    ]);
    let built = build(&decays, &[], &[]).unwrap();
    let modes = &find(&built, "Cs-140").decay.as_ref().unwrap().modes;
    let products: Vec<_> = modes.iter().map(|m| m.products.clone()).collect();
    assert_eq!(products[0], vec![one(&built, "Ba-140")]);
    assert_eq!(products[1], vec![one(&built, "Ba-139m")]);
    assert_eq!(products[2], vec![one(&built, "I-139"), one(&built, "H-1")]);
    assert_eq!(
        products[3],
        vec![
            one(&built, "I-138"),
            Product {
                nuclide: index(&built, "H-1"),
                count: 2
            }
        ]
    );
    assert_eq!(products[4], vec![one(&built, "I-136"), one(&built, "He-4")]);
    let names: Vec<_> = modes.iter().map(|m| m.mode.as_str()).collect();
    assert_eq!(names, ["beta-", "beta-,n", "ec/beta+,p", "p,p", "alpha"]);
    let decay = find(&built, "Cs-140").decay.as_ref().unwrap();
    assert_eq!(decay.half_life, 60.0);
    assert_eq!(
        (
            decay.energy.light,
            decay.energy.electromagnetic,
            decay.energy.heavy
        ),
        (1.0, 2.0, 3.0)
    );
}

#[test]
fn removes_the_parent_without_products_by_spontaneous_fission() {
    let decays = with_light(vec![
        unstable(
            98,
            252,
            0,
            8.3e7,
            vec![
                mode(&[Step::Alpha], 0, 0.97),
                mode(&[Step::SpontaneousFission], 0, 0.03),
            ],
        ),
        stable(96, 248),
    ]);
    let built = build(&decays, &[], &[]).unwrap();
    let modes = &find(&built, "Cf-252").decay.as_ref().unwrap().modes;
    assert_eq!(modes[1].mode, "sf");
    assert!(modes[1].products.is_empty());
}

#[test]
fn rescales_decay_ratios_that_sum_close_to_one() {
    let decays = with_light(vec![
        unstable(
            92,
            238,
            0,
            1.4e17,
            vec![
                mode(&[Step::Alpha], 0, 1.0),
                mode(&[Step::SpontaneousFission], 0, 5.46e-7),
            ],
        ),
        stable(90, 234),
    ]);
    let built = build(&decays, &[], &[]).unwrap();
    let modes = &find(&built, "U-238").decay.as_ref().unwrap().modes;
    let sum = 1.0 + 5.46e-7;
    assert_eq!(modes[0].ratio, 1.0 / sum);
    assert_eq!(modes[1].ratio, 5.46e-7 / sum);
}

#[test]
fn rejects_decay_ratios_far_from_one() {
    let decays = with_light(vec![
        unstable(92, 238, 0, 1.4e17, vec![mode(&[Step::Alpha], 0, 0.9)]),
        stable(90, 234),
    ]);
    assert!(build(&decays, &[], &[]).is_err());
}

#[test]
fn drops_modes_with_zero_ratio() {
    let decays = with_light(vec![
        unstable(
            53,
            135,
            0,
            2.4e4,
            vec![
                mode(&[Step::BetaMinus], 0, 1.0),
                mode(&[Step::Alpha], 0, 0.0),
            ],
        ),
        stable(54, 135),
        stable(51, 131),
    ]);
    let built = build(&decays, &[], &[]).unwrap();
    assert_eq!(find(&built, "I-135").decay.as_ref().unwrap().modes.len(), 1);
}

#[test]
fn treats_an_unstable_nuclide_without_half_life_as_stable() {
    let decays = with_light(vec![
        unstable(
            54,
            136,
            0,
            0.0,
            vec![mode(&[Step::BetaMinus, Step::BetaMinus], 0, 1.0)],
        ),
        stable(56, 136),
    ]);
    let built = build(&decays, &[], &[]).unwrap();
    assert!(find(&built, "Xe-136").decay.is_none());
    assert!(
        built.warnings.iter().any(|w| w.contains("Xe-136")),
        "{:?}",
        built.warnings
    );
}

#[test]
fn keeps_the_loss_but_no_product_for_a_missing_daughter() {
    let decays = with_light(vec![unstable(
        28,
        48,
        0,
        1e-3,
        vec![mode(&[Step::EcBetaPlus], 0, 1.0)],
    )]);
    let built = build(&decays, &[], &[]).unwrap();
    let modes = &find(&built, "Ni-48").decay.as_ref().unwrap().modes;
    assert!(modes[0].products.is_empty());
    assert!(
        built.warnings.iter().any(|w| w.contains("Co-48")),
        "{:?}",
        built.warnings
    );
}

#[test]
fn requires_h1_and_he4() {
    assert!(build(&[stable(2, 4)], &[], &[]).is_err());
}

#[test]
fn maps_mf3_reactions_to_ground_state_products() {
    let decays = with_light(vec![
        stable(26, 56),
        stable(26, 57),
        stable(26, 55),
        stable(26, 54),
        stable(26, 53),
        stable(25, 56),
        stable(24, 53),
    ]);
    // (n,2n) は 875〜891 の部分だけ、(n,p) は 600〜649 の部分だけ
    let neutrons = [neutron(26, 56, 0, &[1, 2, 17, 37, 102, 107, 600, 601, 875])];
    let built = build(&decays, &neutrons, &[]).unwrap();
    let fe56 = find(&built, "Fe-56");
    let got: Vec<_> = fe56
        .reactions
        .iter()
        .map(|r| (r.mt, r.channels[0].products.clone()))
        .collect();
    assert_eq!(
        got,
        vec![
            (16, vec![one(&built, "Fe-55")]),
            (17, vec![one(&built, "Fe-54")]),
            (37, vec![one(&built, "Fe-53")]),
            (102, vec![one(&built, "Fe-57")]),
            (103, vec![one(&built, "Mn-56"), one(&built, "H-1")]),
            (107, vec![one(&built, "Cr-53"), one(&built, "He-4")]),
        ]
    );
    assert!(
        fe56.reactions
            .iter()
            .all(|r| r.channels.len() == 1 && r.channels[0].final_state.is_none())
    );
    assert!(fe56.fission.is_none());
}

#[test]
fn splits_reactions_into_the_final_states_of_mf8() {
    let decays = with_light(vec![
        unstable(95, 241, 0, 1.4e10, vec![mode(&[Step::Alpha], 0, 1.0)]),
        unstable(95, 242, 0, 5.8e4, vec![mode(&[Step::BetaMinus], 0, 1.0)]),
        unstable(
            95,
            242,
            1,
            4.4e9,
            vec![mode(&[Step::IsomericTransition], 0, 1.0)],
        ),
        stable(93, 237),
        stable(96, 242),
    ]);
    let mut am241 = neutron(95, 241, 0, &[1, 102]);
    am241.final_states = vec![(
        102,
        vec![
            FinalStateEntry {
                za: 95242,
                excitation: 0.0,
                mf: 9,
                lfs: 0,
            },
            FinalStateEntry {
                za: 95242,
                excitation: 1e5,
                mf: 9,
                lfs: 2,
            },
        ],
    )];
    let built = build(&decays, &[am241], &[]).unwrap();
    let channels = &find(&built, "Am-241").reactions[0].channels;
    assert_eq!(channels.len(), 2);
    assert_eq!(
        channels[0].final_state,
        Some(FinalState {
            mf: 9,
            lfs: Some(0)
        })
    );
    assert_eq!(channels[0].products, vec![one(&built, "Am-242")]);
    assert_eq!(
        channels[1].final_state,
        Some(FinalState {
            mf: 9,
            lfs: Some(2)
        })
    );
    assert_eq!(channels[1].products, vec![one(&built, "Am-242m")]);
}

fn entry(za: u32, excitation: f64, mf: u32, lfs: u32) -> FinalStateEntry {
    FinalStateEntry {
        za,
        excitation,
        mf,
        lfs,
    }
}

#[test]
fn sends_a_final_state_to_the_ground_state_when_the_product_has_no_isomer() {
    // decay に準安定状態のない生成物の準位（短寿命で、即発の γ で基底状態に落ちる）
    let decays = with_light(vec![stable(95, 241), stable(95, 242)]);
    let mut am241 = neutron(95, 241, 0, &[102]);
    am241.final_states = vec![(
        102,
        vec![entry(95242, 0.0, 10, 0), entry(95242, 3e5, 10, 3)],
    )];
    let built = build(&decays, &[am241], &[]).unwrap();
    let channels = &find(&built, "Am-241").reactions[0].channels;
    assert_eq!(
        channels[1].final_state,
        Some(FinalState {
            mf: 10,
            lfs: Some(3)
        })
    );
    assert_eq!(channels[1].products, vec![one(&built, "Am-242")]);
    assert!(
        built.warnings.iter().any(|w| w.contains("Am-241")),
        "{:?}",
        built.warnings
    );
}

/// Ir-190 の基底状態と 2 つの準安定状態（26.1 keV と 376.4 keV）
fn ir190_decays() -> Vec<DecayEvaluation> {
    let isomer = |state, excitation| DecayEvaluation {
        material: material(77, 190, state, excitation),
        ..unstable(
            77,
            190,
            state,
            1.0,
            vec![mode(&[Step::IsomericTransition], 0, 1.0)],
        )
    };
    with_light(vec![
        stable(77, 191),
        stable(77, 190),
        isomer(1, 26_100.0),
        isomer(2, 376_400.0),
    ])
}

#[test]
fn matches_excited_states_to_isomers_by_the_nearest_energy() {
    // ENDF/B-VIII.0 の Ir-191 (n,2n) の形: 26.0 keV と 376.4 keV は準安定状態、36.15 keV は decay にない準位
    let mut ir191 = neutron(77, 191, 0, &[16]);
    ir191.final_states = vec![(
        16,
        vec![
            entry(77190, 0.0, 10, 0),
            entry(77190, 26_000.0, 10, 1),
            entry(77190, 36_150.0, 10, 2),
            entry(77190, 376_400.0, 10, 3),
        ],
    )];
    let built = build(&ir190_decays(), &[ir191], &[]).unwrap();
    let products: Vec<_> = find(&built, "Ir-191").reactions[0]
        .channels
        .iter()
        .map(|c| c.products.clone())
        .collect();
    assert_eq!(
        products,
        vec![
            vec![one(&built, "Ir-190")],
            vec![one(&built, "Ir-190m")],
            vec![one(&built, "Ir-190")],
            vec![one(&built, "Ir-190m2")],
        ]
    );
    // 26.0 keV と 26.1 keV の違い（1e-3 を超える）と、36.15 keV の準位を基底状態に落としたことを示す
    let warnings = built.warnings.join("\n");
    assert!(warnings.contains("26000"), "{warnings}");
    assert!(warnings.contains("36150"), "{warnings}");
}

#[test]
fn keeps_the_remainder_when_mf8_lists_only_isomers() {
    // ENDF/B-VIII.0 の In-115 (n,γ) の形: MF8 は準安定状態だけを挙げ、残りが基底状態
    let decays = with_light(vec![
        stable(49, 115),
        unstable(49, 116, 0, 14.0, vec![mode(&[Step::BetaMinus], 0, 1.0)]),
        DecayEvaluation {
            material: material(49, 116, 1, 127_267.0),
            ..unstable(49, 116, 1, 3.3e3, vec![mode(&[Step::BetaMinus], 0, 1.0)])
        },
        stable(50, 116),
    ]);
    let mut in115 = neutron(49, 115, 0, &[102]);
    in115.final_states = vec![(102, vec![entry(49116, 127_269.7, 9, 1)])];
    let built = build(&decays, &[in115], &[]).unwrap();
    let channels = &find(&built, "In-115").reactions[0].channels;
    assert_eq!(channels.len(), 2);
    assert_eq!(
        channels[0].final_state,
        Some(FinalState {
            mf: 9,
            lfs: Some(1)
        })
    );
    assert_eq!(channels[0].products, vec![one(&built, "In-116m")]);
    assert_eq!(
        channels[1].final_state,
        Some(FinalState { mf: 9, lfs: None })
    );
    assert_eq!(channels[1].products, vec![one(&built, "In-116")]);
}

#[test]
fn rejects_a_final_state_of_another_nuclide() {
    let decays = with_light(vec![stable(95, 241), stable(95, 242), stable(95, 240)]);
    let mut am241 = neutron(95, 241, 0, &[102]);
    am241.final_states = vec![(
        102,
        vec![FinalStateEntry {
            za: 95240,
            excitation: 0.0,
            mf: 9,
            lfs: 0,
        }],
    )];
    assert!(build(&decays, &[am241], &[]).is_err());
}

#[test]
fn keeps_the_loss_but_no_product_for_a_missing_reaction_product() {
    let decays = with_light(vec![stable(94, 245), stable(94, 246)]);
    let built = build(&decays, &[neutron(94, 245, 0, &[102, 103])], &[]).unwrap();
    let reactions = &find(&built, "Pu-245").reactions;
    assert_eq!(reactions[1].mt, 103);
    // 生成物の Np-245 がないので、H-1 だけを生む
    assert_eq!(reactions[1].channels[0].products, vec![one(&built, "H-1")]);
    assert!(
        built.warnings.iter().any(|w| w.contains("Np-245")),
        "{:?}",
        built.warnings
    );
}

#[test]
fn rejects_neutron_data_for_a_nuclide_outside_decay() {
    let decays = with_light(vec![]);
    assert!(build(&decays, &[neutron(26, 56, 0, &[102])], &[]).is_err());
}

fn raw(za: u32, state: u32, value: f64) -> RawYield {
    RawYield { za, state, value }
}

fn fp_decays() -> Vec<DecayEvaluation> {
    with_light(vec![
        unstable(92, 235, 0, 2.2e16, vec![mode(&[Step::Alpha], 0, 1.0)]),
        unstable(92, 238, 0, 1.4e17, vec![mode(&[Step::Alpha], 0, 1.0)]),
        unstable(92, 239, 0, 1.4e3, vec![mode(&[Step::BetaMinus], 0, 1.0)]),
        stable(90, 231),
        stable(90, 234),
        stable(93, 239),
        unstable(53, 135, 0, 2.4e4, vec![mode(&[Step::BetaMinus], 0, 1.0)]),
        stable(54, 135),
        stable(52, 135),
        stable(55, 140),
    ])
}

fn u235_yields() -> YieldEvaluation {
    YieldEvaluation {
        material: material(92, 235, 0, 0.0),
        tables: vec![
            (5e5, vec![raw(53135, 0, 0.05), raw(54135, 0, 0.004)]),
            (
                0.0253,
                vec![
                    raw(53135, 0, 0.06),
                    // decay にない準安定状態は基底状態へ
                    raw(54135, 1, 0.001),
                    raw(54135, 0, 0.002),
                    // decay にない基底状態は β− の向きへ（Sb-140 → Te-140 → I-140 → Xe-140 → Cs-140）
                    raw(51140, 0, 1e-8),
                    raw(52135, 0, 0.0),
                ],
            ),
        ],
    }
}

#[test]
fn keeps_all_yield_energies_in_ascending_order() {
    let built = build(
        &fp_decays(),
        &[neutron(92, 235, 0, &[18, 102])],
        &[u235_yields()],
    )
    .unwrap();
    let fission = find(&built, "U-235").fission.as_ref().unwrap();
    assert_eq!(fission.yield_source, "U-235");
    assert_eq!(fission.energies, vec![0.0253, 5e5]);
    let thermal: Vec<_> = fission.yields[0]
        .iter()
        .map(|y| (built.chain.nuclides[y.nuclide].name.as_str(), y.value))
        .collect();
    assert_eq!(
        thermal,
        vec![("I-135", 0.06), ("Xe-135", 0.003), ("Cs-140", 1e-8)]
    );
    assert!(
        built.warnings.iter().any(|w| w.contains("Xe-135m")),
        "{:?}",
        built.warnings
    );
}

#[test]
fn borrows_yields_from_the_nearest_mass_of_the_same_element() {
    let u238_yields = YieldEvaluation {
        material: material(92, 238, 0, 0.0),
        tables: vec![(5e5, vec![raw(53135, 0, 0.07)])],
    };
    let neutrons = [
        neutron(92, 235, 0, &[18]),
        neutron(92, 238, 0, &[18]),
        neutron(92, 239, 0, &[18]),
    ];
    let built = build(&fp_decays(), &neutrons, &[u235_yields(), u238_yields]).unwrap();
    let fission = find(&built, "U-239").fission.as_ref().unwrap();
    assert_eq!(fission.yield_source, "U-238");
    assert_eq!(fission.energies, vec![5e5]);
    assert!(
        built.warnings.iter().any(|w| w.contains("U-239")),
        "{:?}",
        built.warnings
    );
}

#[test]
fn prefers_the_other_state_of_the_same_mass() {
    let decays = with_light(vec![
        unstable(95, 242, 0, 5.8e4, vec![mode(&[Step::BetaMinus], 0, 1.0)]),
        unstable(
            95,
            242,
            1,
            4.4e9,
            vec![mode(&[Step::IsomericTransition], 0, 1.0)],
        ),
        unstable(95, 241, 0, 1.4e10, vec![mode(&[Step::Alpha], 0, 1.0)]),
        stable(96, 242),
        stable(93, 237),
        stable(53, 135),
    ]);
    let yields = |state, a| YieldEvaluation {
        material: material(95, a, state, 0.0),
        tables: vec![(0.0253, vec![raw(53135, 0, 0.06)])],
    };
    let neutrons = [
        neutron(95, 241, 0, &[18]),
        neutron(95, 242, 0, &[18]),
        neutron(95, 242, 1, &[18]),
    ];
    let built = build(&decays, &neutrons, &[yields(0, 241), yields(1, 242)]).unwrap();
    assert_eq!(
        find(&built, "Am-242")
            .fission
            .as_ref()
            .unwrap()
            .yield_source,
        "Am-242m"
    );
}

#[test]
fn borrows_from_the_nearest_element_when_the_element_has_no_yields() {
    let decays = with_light(vec![stable(91, 233), stable(92, 235), stable(53, 135)]);
    let yields = YieldEvaluation {
        material: material(92, 235, 0, 0.0),
        tables: vec![(0.0253, vec![raw(53135, 0, 0.06)])],
    };
    let neutrons = [neutron(91, 233, 0, &[19]), neutron(92, 235, 0, &[18])];
    let built = build(&decays, &neutrons, &[yields]).unwrap();
    assert_eq!(
        find(&built, "Pa-233")
            .fission
            .as_ref()
            .unwrap()
            .yield_source,
        "U-235"
    );
}

#[test]
fn ignores_neutron_data_of_the_neutron_itself() {
    let built = build(&with_light(vec![]), &[neutron(0, 1, 0, &[2])], &[]).unwrap();
    assert_eq!(built.chain.nuclides.len(), 2);
}

#[test]
fn merges_a_daughter_that_is_also_the_emitted_light_nuclide() {
    // Li-8 の β−α は娘も He-4、Be-7 の (n,α) は残りの核も He-4
    let decays = with_light(vec![
        unstable(
            3,
            8,
            0,
            0.84,
            vec![mode(&[Step::BetaMinus, Step::Alpha], 0, 1.0)],
        ),
        unstable(4, 7, 0, 4.6e6, vec![mode(&[Step::EcBetaPlus], 0, 1.0)]),
        stable(3, 7),
    ]);
    let built = build(&decays, &[neutron(4, 7, 0, &[107])], &[]).unwrap();
    let two_alphas = vec![Product {
        nuclide: index(&built, "He-4"),
        count: 2,
    }];
    let li8 = find(&built, "Li-8").decay.as_ref().unwrap();
    assert_eq!(li8.modes[0].products, two_alphas);
    assert_eq!(
        find(&built, "Be-7").reactions[0].channels[0].products,
        two_alphas
    );
}
