use std::collections::BTreeMap;
use std::f64::consts::LN_2;

use serde::Serialize;

use super::*;

const I135: usize = 2;
const XE135: usize = 3;
const XE136: usize = 4;
const CS135: usize = 5;
const U235: usize = 6;
const U236: usize = 7;

fn product(nuclide: usize) -> Product {
    Product { nuclide, count: 1 }
}

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

fn decaying(name: &str, z: u32, a: u32, half_life: f64, modes: Vec<DecayMode>) -> Nuclide {
    Nuclide {
        decay: Some(Decay {
            half_life,
            energy: DecayEnergy {
                light: 1e5,
                electromagnetic: 2e5,
                heavy: 0.0,
            },
            modes,
        }),
        ..stable(name, z, a)
    }
}

fn mode(mode: &str, ratio: f64, products: Vec<Product>) -> DecayMode {
    DecayMode {
        mode: mode.into(),
        ratio,
        products,
    }
}

fn capture(target: usize) -> Reaction {
    Reaction {
        mt: 102,
        channels: vec![Channel {
            final_state: None,
            products: vec![product(target)],
        }],
    }
}

/// H-1、He-4、I-135 → Xe-135 → Cs-135、Xe-135 (n,γ) Xe-136、U-235 の α（娘はチェーンにない）と (n,γ) と核分裂
fn small_chain() -> Chain {
    let mut xe135 = decaying(
        "Xe-135",
        54,
        135,
        32_904.0,
        vec![mode("beta-", 1.0, vec![product(CS135)])],
    );
    xe135.reactions.push(capture(XE136));
    let mut u235 = decaying(
        "U-235",
        92,
        235,
        2.221e16,
        vec![mode("alpha", 1.0, vec![product(1)])],
    );
    u235.reactions.push(capture(U236));
    u235.fission = Some(Fission {
        q: 1.93e8,
        yield_source: "U-235".into(),
        energies: vec![0.0253, 5e5],
        yields: vec![
            vec![
                Yield {
                    nuclide: I135,
                    value: 0.06,
                },
                Yield {
                    nuclide: XE135,
                    value: 0.002,
                },
            ],
            vec![
                Yield {
                    nuclide: I135,
                    value: 0.05,
                },
                Yield {
                    nuclide: XE135,
                    value: 0.004,
                },
            ],
        ],
    });
    Chain {
        nuclides: vec![
            stable("H-1", 1, 1),
            stable("He-4", 2, 4),
            decaying(
                "I-135",
                53,
                135,
                23_652.0,
                vec![mode("beta-", 1.0, vec![product(XE135)])],
            ),
            xe135,
            stable("Xe-136", 54, 136),
            stable("Cs-135", 55, 135),
            u235,
            stable("U-236", 92, 236),
        ],
    }
}

fn invalid(chain: &Chain) -> String {
    match chain.check() {
        Err(ChainError::Invalid(message)) => message,
        other => panic!("expected an invalid chain, got {other:?}"),
    }
}

#[test]
fn round_trips_bit_exactly() {
    let mut chain = small_chain();
    chain.nuclides[XE135].decay.as_mut().unwrap().half_life = 0.1 + 0.2;
    let bytes = chain.to_msgpack().unwrap();
    assert_eq!(Chain::from_msgpack(&bytes).unwrap(), chain);
}

#[test]
fn written_as_a_map_with_the_format_name() {
    let bytes = small_chain().to_msgpack().unwrap();
    let value = xlc_format::msgpack::read_value(&bytes).unwrap();
    let format = xlc_format::msgpack::map(&value)
        .unwrap()
        .iter()
        .find(|(k, _)| k.as_str() == Some("format"))
        .map(|(_, v)| v.as_str().map(str::to_string));
    assert_eq!(format, Some(Some("xlc-chain".to_string())));
}

#[test]
fn rejects_another_format() {
    #[derive(Serialize)]
    struct File {
        format: &'static str,
        nuclides: Vec<Nuclide>,
    }
    let bytes = xlc_format::msgpack::to_named(&File {
        format: "xlc-multigroup",
        nuclides: small_chain().nuclides,
    })
    .unwrap();
    assert!(matches!(
        Chain::from_msgpack(&bytes),
        Err(ChainError::Format(f)) if f == "xlc-multigroup"
    ));
}

#[test]
fn rejects_nested_structs_written_as_arrays() {
    #[derive(Serialize)]
    struct File {
        format: &'static str,
        nuclides: Vec<Nuclide>,
    }
    // 構造体を配列で書く形（フィールド名なし）
    let bytes = rmp_serde::to_vec(&File {
        format: "xlc-chain",
        nuclides: small_chain().nuclides,
    })
    .unwrap();
    assert!(Chain::from_msgpack(&bytes).is_err());
}

#[test]
fn rejects_an_omitted_decay() {
    #[derive(Serialize)]
    struct Partial {
        name: &'static str,
        z: u32,
        a: u32,
        state: u32,
        excitation: f64,
        reactions: Vec<Reaction>,
        fission: Option<Fission>,
    }
    #[derive(Serialize)]
    struct File {
        format: &'static str,
        nuclides: Vec<Partial>,
    }
    let bytes = xlc_format::msgpack::to_named(&File {
        format: "xlc-chain",
        nuclides: vec![Partial {
            name: "H-1",
            z: 1,
            a: 1,
            state: 0,
            excitation: 0.0,
            reactions: vec![],
            fission: None,
        }],
    })
    .unwrap();
    assert!(Chain::from_msgpack(&bytes).is_err());
}

#[test]
fn checks_the_written_chain() {
    let mut chain = small_chain();
    chain.nuclides[0].name = "H1".into();
    assert!(chain.to_msgpack().is_err());
}

#[test]
fn rejects_names_out_of_the_convention() {
    let mut chain = small_chain();
    chain.nuclides[U236].name = "U236".into();
    assert!(invalid(&chain).contains("U236"));
}

#[test]
fn rejects_nuclides_out_of_order() {
    let mut chain = small_chain();
    chain.nuclides.swap(0, 1);
    invalid(&chain);
}

#[test]
fn rejects_products_outside_the_chain() {
    let mut chain = small_chain();
    chain.nuclides[I135].decay.as_mut().unwrap().modes[0].products[0].nuclide = 8;
    invalid(&chain);
    let mut chain = small_chain();
    chain.nuclides[U235].fission.as_mut().unwrap().yields[1][0].nuclide = 99;
    invalid(&chain);
}

#[test]
fn rejects_decay_ratios_that_do_not_sum_to_one() {
    let mut chain = small_chain();
    chain.nuclides[I135].decay.as_mut().unwrap().modes[0].ratio = 0.9;
    invalid(&chain);
}

#[test]
fn rejects_non_positive_half_lives() {
    let mut chain = small_chain();
    chain.nuclides[I135].decay.as_mut().unwrap().half_life = 0.0;
    invalid(&chain);
}

#[test]
fn rejects_yields_that_do_not_match_the_energies() {
    let mut chain = small_chain();
    chain.nuclides[U235].fission.as_mut().unwrap().energies = vec![0.0253];
    invalid(&chain);
    let mut chain = small_chain();
    chain.nuclides[U235].fission.as_mut().unwrap().energies = vec![5e5, 0.0253];
    invalid(&chain);
}

#[test]
fn rejects_non_positive_or_non_finite_fission_q() {
    for q in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut chain = small_chain();
        chain.nuclides[U235].fission.as_mut().unwrap().q = q;
        invalid(&chain);
    }
}

#[test]
fn rejects_unlabelled_channels_of_a_split_reaction() {
    let mut chain = small_chain();
    let reaction = &mut chain.nuclides[U235].reactions[0];
    reaction.channels.push(reaction.channels[0].clone());
    invalid(&chain);
}

#[test]
fn bounds_the_yield_energies_at_geometric_means() {
    let fission = small_chain().nuclides[U235].fission.clone().unwrap();
    let bounds = fission.yield_bounds();
    assert_eq!(bounds.len(), 1);
    assert!((bounds[0] - (0.0253f64 * 5e5).sqrt()).abs() <= 1e-12 * bounds[0]);
    let one = Fission {
        energies: vec![0.0253],
        yields: vec![vec![]],
        ..fission
    };
    assert!(one.yield_bounds().is_empty());
}

#[test]
fn zero_rates_follow_the_shape_of_the_chain() {
    let chain = small_chain();
    let rates = chain.zero_rates();
    assert_eq!(rates.len(), chain.nuclides.len());
    assert_eq!(rates[XE135].reactions, vec![vec![0.0]]);
    assert!(rates[I135].reactions.is_empty());
    assert_eq!(rates[U235].fission, 0.0);
    assert_eq!(rates[U235].yield_weights.len(), 2);
    assert!(rates[I135].yield_weights.is_empty());
}

/// 同じ位置の項を足した (行, 列) → 値
fn summed(entries: &[(usize, usize, f64)]) -> BTreeMap<(usize, usize), f64> {
    let mut map = BTreeMap::new();
    for &(row, col, value) in entries {
        *map.entry((row, col)).or_insert(0.0) += value;
    }
    map
}

#[test]
fn assembles_decay_reaction_and_fission_terms() {
    let chain = small_chain();
    let mut rates = chain.zero_rates();
    rates[XE135].reactions[0][0] = 1e-4;
    rates[U235].reactions[0][0] = 3e-9;
    rates[U235].fission = 1e-8;
    rates[U235].yield_weights = vec![0.75, 0.25];
    let got = summed(&chain.matrix_entries(&rates).unwrap());

    let lambda = |t: f64| LN_2 / t;
    let (l_i, l_xe, l_u) = (lambda(23_652.0), lambda(32_904.0), lambda(2.221e16));
    let mut expected = BTreeMap::new();
    expected.insert((I135, I135), -l_i);
    expected.insert((XE135, I135), l_i);
    expected.insert((XE135, XE135), -l_xe - 1e-4);
    expected.insert((CS135, XE135), l_xe);
    expected.insert((XE136, XE135), 1e-4);
    expected.insert((U235, U235), -l_u - 3e-9 - 1e-8);
    expected.insert((1, U235), l_u);
    expected.insert((U236, U235), 3e-9);
    expected.insert((I135, U235), 1e-8 * (0.75 * 0.06 + 0.25 * 0.05));
    expected.insert((XE135, U235), 1e-8 * (0.75 * 0.002 + 0.25 * 0.004));

    assert_eq!(
        got.keys().collect::<Vec<_>>(),
        expected.keys().collect::<Vec<_>>()
    );
    for (key, value) in &expected {
        let g = got[key];
        assert!(
            (g - value).abs() <= 1e-14 * value.abs(),
            "{key:?}: {g} vs {value}"
        );
    }
}

#[test]
fn counts_light_products_by_multiplicity() {
    let mut chain = small_chain();
    chain.nuclides[U235].decay.as_mut().unwrap().modes[0].products = vec![Product {
        nuclide: 1,
        count: 2,
    }];
    let got = summed(&chain.matrix_entries(&chain.zero_rates()).unwrap());
    let l_u = LN_2 / 2.221e16;
    assert!((got[&(1, U235)] - 2.0 * l_u).abs() <= 1e-14 * l_u);
}

#[test]
fn rejects_rates_of_another_shape() {
    let chain = small_chain();
    let mut rates = chain.zero_rates();
    rates[XE135].reactions[0].push(0.0);
    assert!(matches!(
        chain.matrix_entries(&rates),
        Err(ChainError::Rates(_))
    ));
    let mut rates = chain.zero_rates();
    rates.pop();
    assert!(matches!(
        chain.matrix_entries(&rates),
        Err(ChainError::Rates(_))
    ));
}

#[test]
fn rejects_negative_or_non_finite_rates() {
    let chain = small_chain();
    for bad in [-1e-9, f64::NAN, f64::INFINITY] {
        let mut rates = chain.zero_rates();
        rates[XE135].reactions[0][0] = bad;
        assert!(
            matches!(chain.matrix_entries(&rates), Err(ChainError::Rates(_))),
            "{bad}"
        );
    }
}

#[test]
fn rejects_yield_weights_that_do_not_sum_to_one() {
    let chain = small_chain();
    let mut rates = chain.zero_rates();
    rates[U235].fission = 1e-8;
    rates[U235].yield_weights = vec![0.5, 0.25];
    assert!(matches!(
        chain.matrix_entries(&rates),
        Err(ChainError::Rates(_))
    ));
}

#[test]
fn rejects_fission_rates_of_nuclides_without_fission() {
    let chain = small_chain();
    let mut rates = chain.zero_rates();
    rates[XE135].fission = 1e-8;
    assert!(matches!(
        chain.matrix_entries(&rates),
        Err(ChainError::Rates(_))
    ));
}

#[test]
fn the_matrix_depletes_like_the_entries() {
    let chain = small_chain();
    let mut rates = chain.zero_rates();
    rates[U235].fission = 1e-8;
    rates[U235].yield_weights = vec![1.0, 0.0];
    let a = chain.matrix(&rates).unwrap();
    let mut n0 = vec![0.0; chain.nuclides.len()];
    n0[U235] = 1.0;
    let n = crate::deplete(&a, &n0, 1e5, crate::CramOrder::Order48).unwrap();
    assert!(n[I135] > 0.0 && n[XE135] > 0.0, "{n:?}");
}

#[test]
fn rejects_names_that_disagree_with_z_a_and_state() {
    for name in ["U-238", "Np-235", "U-235m", "Xe-136"] {
        let mut chain = small_chain();
        chain.nuclides[U235].name = name.into();
        assert!(invalid(&chain).contains(name), "{name}");
    }
}

#[test]
fn rejects_a_yield_source_outside_the_chain() {
    let mut chain = small_chain();
    chain.nuclides[U235].fission.as_mut().unwrap().yield_source = "Pu-239".into();
    invalid(&chain);
}

#[test]
fn rejects_repeated_or_unknown_reactions() {
    let mut chain = small_chain();
    let capture = chain.nuclides[U235].reactions[0].clone();
    chain.nuclides[U235].reactions.push(capture);
    invalid(&chain);
    let mut chain = small_chain();
    chain.nuclides[U235].reactions[0].mt = 18;
    invalid(&chain);
}

fn split_capture(states: &[(u32, Option<u32>)]) -> Chain {
    let mut chain = small_chain();
    chain.nuclides[U235].reactions[0].channels = states
        .iter()
        .map(|&(mf, lfs)| Channel {
            final_state: Some(FinalState { mf, lfs }),
            products: vec![product(U236)],
        })
        .collect();
    chain
}

#[test]
fn accepts_a_remainder_channel_of_a_split_reaction() {
    split_capture(&[(9, Some(1)), (9, None)]).check().unwrap();
}

#[test]
fn rejects_repeated_final_states_or_other_files() {
    invalid(&split_capture(&[(9, Some(1)), (9, Some(1))]));
    invalid(&split_capture(&[(6, Some(0)), (6, Some(1))]));
}

#[test]
fn rejects_the_same_state_in_two_channels() {
    // 同じ状態を別のファイルから、または基底状態を LFS 0 と残りの両方で数えない
    invalid(&split_capture(&[(9, Some(1)), (10, Some(1))]));
    invalid(&split_capture(&[(9, Some(0)), (9, None)]));
    invalid(&split_capture(&[(9, Some(1)), (9, None), (10, None)]));
}

#[test]
fn rejects_products_that_are_the_parent() {
    let mut chain = small_chain();
    chain.nuclides[U235].reactions[0].channels[0].products[0].nuclide = U235;
    invalid(&chain);
    let mut chain = small_chain();
    chain.nuclides[I135].decay.as_mut().unwrap().modes[0].products[0].nuclide = I135;
    invalid(&chain);
    let mut chain = small_chain();
    chain.nuclides[U235].fission.as_mut().unwrap().yields[0][0].nuclide = U235;
    invalid(&chain);
}

#[test]
fn rejects_empty_yield_tables_and_sources_without_fission() {
    let mut chain = small_chain();
    chain.nuclides[U235].fission.as_mut().unwrap().yields[1].clear();
    invalid(&chain);
    let mut chain = small_chain();
    chain.nuclides[U235].fission.as_mut().unwrap().yield_source = "U-236".into();
    invalid(&chain);
}

#[test]
fn rejects_a_mass_number_below_one() {
    let mut chain = small_chain();
    chain.nuclides[0].a = 0;
    chain.nuclides[0].name = "H-0".into();
    invalid(&chain);
}

#[test]
fn rejects_a_nested_struct_written_as_an_array() {
    // 最上位は map のまま、崩壊のエネルギー（入れ子の構造体）だけを配列にする
    let bytes = small_chain().to_msgpack().unwrap();
    let mut value = xlc_format::msgpack::read_value(&bytes).unwrap();
    let rmpv::Value::Map(top) = &mut value else {
        panic!("not a map")
    };
    let nuclides = top
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("nuclides"))
        .unwrap();
    let rmpv::Value::Array(nuclides) = &mut nuclides.1 else {
        panic!("not an array")
    };
    let rmpv::Value::Map(i135) = &mut nuclides[I135] else {
        panic!("not a map")
    };
    let decay = i135
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("decay"))
        .unwrap();
    let rmpv::Value::Map(decay) = &mut decay.1 else {
        panic!("not a map")
    };
    let energy = decay
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("energy"))
        .unwrap();
    energy.1 = rmpv::Value::Array(vec![1e5.into(), 2e5.into(), 0.0.into()]);
    let mut changed = Vec::new();
    rmpv::encode::write_value(&mut changed, &value).unwrap();
    assert!(matches!(
        Chain::from_msgpack(&changed),
        Err(ChainError::Msgpack(MsgpackError::NotMap))
    ));
}

#[test]
fn rejects_unknown_decay_modes_and_products_of_spontaneous_fission() {
    let mut chain = small_chain();
    chain.nuclides[I135].decay.as_mut().unwrap().modes[0].mode = "beta".into();
    invalid(&chain);
    let mut chain = small_chain();
    chain.nuclides[U235].decay.as_mut().unwrap().modes[0].mode = "sf".into();
    invalid(&chain);
}

#[test]
fn checks_the_chain_before_assembling() {
    let mut chain = small_chain();
    chain.nuclides[I135].decay.as_mut().unwrap().modes[0].ratio = 0.9;
    let rates = chain.zero_rates();
    assert!(matches!(
        chain.matrix_entries(&rates),
        Err(ChainError::Invalid(_))
    ));
}

#[test]
fn leaves_out_zero_rate_terms() {
    let chain = small_chain();
    let entries = chain.matrix_entries(&chain.zero_rates()).unwrap();
    assert!(entries.iter().all(|&(_, _, v)| v != 0.0), "{entries:?}");
}

#[test]
fn rejects_a_product_listed_twice() {
    let mut chain = small_chain();
    chain.nuclides[U235].reactions[0].channels[0]
        .products
        .push(product(U236));
    invalid(&chain);
    let mut chain = small_chain();
    let first = chain.nuclides[U235].fission.as_ref().unwrap().yields[0][0];
    chain.nuclides[U235].fission.as_mut().unwrap().yields[0].push(first);
    invalid(&chain);
}

#[test]
fn rejects_an_excited_ground_state() {
    let mut chain = small_chain();
    chain.nuclides[U236].excitation = 1e3;
    invalid(&chain);
}
