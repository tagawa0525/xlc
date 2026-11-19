use xlc_depletion::chain::{Channel, FinalState, Fission, Nuclide, Product, Reaction, Yield};
use xlc_format::group::GroupStructure;
use xlc_format::library::{
    ActivationData, ActivationTemperature, ChannelState, GroupChannel, GroupReaction,
};

use super::*;

const T: f64 = 600.0;
/// (n,2n)、捕獲の 2 つの状態、(n,p)、(n,α)、核分裂の群の断面積 (barn)
const N2N: [f64; 2] = [0.3, 0.0];
const CAPTURE_GROUND: [f64; 2] = [2.0, 40.0];
const CAPTURE_ISOMER: [f64; 2] = [1.0, 10.0];
const NP: [f64; 2] = [0.02, 0.5];
const NA: [f64; 2] = [0.01, 0.25];
const FISSION: [f64; 2] = [0.5, 7.0];

fn reaction(mt: u32, channels: &[(Option<ChannelState>, [f64; 2])]) -> GroupReaction {
    GroupReaction {
        mt,
        channels: channels
            .iter()
            .map(|(final_state, sigma)| GroupChannel {
                final_state: *final_state,
                cross_section: sigma.to_vec(),
            })
            .collect(),
    }
}

/// 2 群の、輸送する核種のないライブラリ。放射化の核種 Am-242 は (n,2n)、2 つの状態への捕獲、(n,p)、(n,α)、核分裂を持つ
fn library() -> Library {
    let [ground, isomer] = [0, 1].map(|lfs| {
        Some(ChannelState {
            mf: 9,
            lfs: Some(lfs),
        })
    });
    Library {
        groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
        sigma0: vec![1.0e10],
        pointwise_range: None,
        inverse_velocity: vec![5.0e-8, 2.0e-6],
        nuclides: Vec::new(),
        activation: vec![ActivationData {
            name: "Am-242".into(),
            za: 95242.0,
            awr: 239.98,
            temperatures: vec![ActivationTemperature {
                temperature: T,
                fission: FISSION.to_vec(),
                reactions: vec![
                    reaction(16, &[(None, N2N)]),
                    reaction(102, &[(ground, CAPTURE_GROUND), (isomer, CAPTURE_ISOMER)]),
                    reaction(103, &[(None, NP)]),
                    reaction(107, &[(None, NA)]),
                ],
            }],
        }],
    }
}

fn chain_reaction(mt: u32, states: &[Option<FinalState>]) -> Reaction {
    Reaction {
        mt,
        channels: states
            .iter()
            .map(|final_state| Channel {
                final_state: *final_state,
                products: vec![Product {
                    nuclide: 1,
                    count: 1,
                }],
            })
            .collect(),
    }
}

/// Am-242（ライブラリと同じ反応と核分裂）と、その生成物のチェーン
fn chain() -> Chain {
    let nuclide = |name: &str, a: u32| Nuclide {
        name: name.into(),
        z: 95,
        a,
        state: 0,
        excitation: 0.0,
        decay: None,
        reactions: vec![],
        fission: None,
    };
    let [ground, isomer] = [0, 1].map(|lfs| {
        Some(FinalState {
            mf: 9,
            lfs: Some(lfs),
        })
    });
    let mut am242 = nuclide("Am-242", 242);
    am242.reactions = vec![
        chain_reaction(16, &[None]),
        chain_reaction(102, &[ground, isomer]),
        chain_reaction(103, &[None]),
        chain_reaction(107, &[None]),
    ];
    let table = vec![Yield {
        nuclide: 1,
        value: 2.0,
    }];
    am242.fission = Some(Fission {
        q: 2.0e8,
        yield_source: "Am-242".into(),
        energies: vec![0.0253, 5e5],
        yields: vec![table.clone(), table],
    });
    Chain {
        nuclides: vec![am242, nuclide("Am-243", 243)],
    }
}

fn close(got: f64, want: f64) {
    assert!((got - want).abs() <= 1e-13 * want.abs(), "{got} vs {want}");
}

#[test]
fn absorption_sums_capture_proton_and_alpha_reactions_only() {
    let library = library();
    let density = 3.0e-5;
    let got = absorption(&library, [("Am-242", density)], T).unwrap();
    // 中性子を増やす (n,2n) と、生成を数えられない核分裂は入らない
    for g in 0..2 {
        close(
            got[g],
            density * (CAPTURE_GROUND[g] + CAPTURE_ISOMER[g] + NP[g] + NA[g]),
        );
    }
    // 数密度が 0 の核種は引かない（データのない温度でも止まらない）
    assert_eq!(
        absorption(&library, [("Am-242", 0.0)], 900.0).unwrap(),
        vec![0.0, 0.0]
    );
    assert!(absorption(&library, [("Am-242", density)], 900.0).is_err());
    assert!(absorption(&library, [("Am-243", density)], T).is_err());
    assert!(absorption(&library, [("Am-242", -1.0)], T).is_err());
}

#[test]
fn activation_data_are_matched_at_the_exact_temperature() {
    let library = library();
    assert!(has_data_at(&library, "Am-242", T));
    assert!(!has_data_at(&library, "Am-242", 900.0));
    assert!(!has_data_at(&library, "Am-243", T));
}

fn region(regions: &[usize]) -> DepletionRegion {
    DepletionRegion {
        material: "fuel".into(),
        cell: 0,
        cells: vec![[0, 0]],
        regions: regions.to_vec(),
        volume: 1.0,
        sector: None,
    }
}

/// 全断面積が `total`、群の中の散乱が `scatter` の 2 群の物質
fn material(total: [f64; 2], scatter: [f64; 2]) -> Material {
    let mut m = Material::void(2);
    m.total = total.to_vec();
    m.scatter = vec![vec![scatter[0], 0.0], vec![0.0, scatter[1]]];
    m
}

#[test]
fn overlay_adds_the_absorption_times_the_sph_factor_to_the_depletion_regions() {
    let fuel = material([1.0, 2.0], [0.4, 0.5]);
    let water = material([3.0, 4.0], [2.9, 3.9]);
    let materials = [fuel.clone(), water];
    // 領域 0〜2 は燃料（0 と 1 が 1 つ目、2 が 2 つ目の燃焼の領域）、3 は水
    let region_material = [0, 0, 0, 1];
    let regions = [region(&[0, 1]), region(&[2])];
    let absorption = [vec![0.1, 0.2], vec![0.3, 0.4]];
    let sph = [
        vec![1.0, 0.9],
        vec![1.0, 0.9],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    let (unique, index) =
        overlay(&materials, &region_material, &regions, &absorption, &sph).unwrap();
    // 同じ物質、同じ燃焼の領域、同じ因子の領域は 1 つの物質を共有し、燃焼の領域でない領域は元の物質のまま
    assert_eq!(index, [2, 2, 3, 1]);
    assert_eq!(unique.len(), 4);
    assert_eq!(unique[..2], materials);
    assert_eq!(unique[2].total, [1.0 + 0.1, 2.0 + 0.9 * 0.2]);
    assert_eq!(unique[3].total, [1.0 + 0.3, 2.0 + 0.4]);
    for m in &unique[2..] {
        // 散乱と核分裂は変わらず、足した分はすべて吸収になる
        assert_eq!(m.scatter, fuel.scatter);
        assert_eq!(m.nu_fission, fuel.nu_fission);
        assert_eq!(m.fission, fuel.fission);
    }
    close(unique[2].absorption()[1], 2.0 - 0.5 + 0.9 * 0.2);
}

#[test]
fn overlay_splits_a_material_whose_regions_have_other_sph_factors() {
    let materials = [material([1.0, 2.0], [0.4, 0.5])];
    let regions = [region(&[0, 1])];
    let absorption = [vec![0.1, 0.2]];
    let sph = [vec![1.0, 0.9], vec![1.0, 0.8]];
    let (unique, index) = overlay(&materials, &[0, 0], &regions, &absorption, &sph).unwrap();
    assert_eq!(index, [1, 2]);
    assert_eq!(unique[1].total[1], 2.0 + 0.9 * 0.2);
    assert_eq!(unique[2].total[1], 2.0 + 0.8 * 0.2);
}

#[test]
fn overlay_needs_an_absorption_per_depletion_region_and_a_factor_per_region() {
    let materials = [material([1.0, 2.0], [0.4, 0.5])];
    let regions = [region(&[0])];
    let sph = [vec![1.0, 1.0]];
    assert!(overlay(&materials, &[0], &regions, &[], &sph).is_err());
    assert!(overlay(&materials, &[0], &regions, &[vec![0.1, 0.2]], &[]).is_err());
    assert!(overlay(&materials, &[0], &regions, &[vec![0.1]], &sph).is_err());
    assert!(overlay(&materials, &[0], &regions, &[vec![-0.1, 0.2]], &sph).is_err());
}

#[test]
fn rates_fold_the_dilute_cross_sections_with_the_volume_weighted_flux() {
    let (library, chain) = (library(), chain());
    let bounds = library.groups.bounds().to_vec();
    let mut rates = chain.zero_rates();
    // 体積 1 と 3 の領域の中性子束
    let (first, second) = ([2.0, 1.0], [4.0, 0.5]);
    let members: [(f64, &[f64]); 2] = [(1.0, &first), (3.0, &second)];
    add_rates(&chain, &[0], &library, T, &members, &bounds, &mut rates).unwrap();
    let flux = [(1.0 * 2.0 + 3.0 * 4.0) / 4.0, (1.0 * 1.0 + 3.0 * 0.5) / 4.0];
    let rate = |sigma: [f64; 2]| (sigma[0] * flux[0] + sigma[1] * flux[1]) * 1.0e-24;
    // (n,2n) も核分裂も、燃焼の反応率には数える
    close(rates[0].reactions[0][0], rate(N2N));
    close(rates[0].reactions[1][0], rate(CAPTURE_GROUND));
    close(rates[0].reactions[1][1], rate(CAPTURE_ISOMER));
    close(rates[0].reactions[2][0], rate(NP));
    close(rates[0].reactions[3][0], rate(NA));
    close(rates[0].fission, rate(FISSION));
    // 収率の重みは、収率の境界 √(0.0253 · 5e5) eV より下と上の核分裂率の割合。境界を含む上の群はレサジーの幅の比で分ける（ADR-0183）
    let (fast, thermal) = (FISSION[0] * flux[0], FISSION[1] * flux[1]);
    let below = ((0.0253_f64 * 5.0e5).sqrt() / 0.625).ln() / (1.0e7_f64 / 0.625).ln();
    close(
        rates[0].yield_weights[0],
        (thermal + below * fast) / (fast + thermal),
    );
    close(
        rates[0].yield_weights[1],
        (1.0 - below) * fast / (fast + thermal),
    );
    // 対象でない核種の反応率は 0 のまま
    assert_eq!(rates[1], chain.zero_rates()[1]);
}

#[test]
fn rates_need_the_chain_reactions_and_fission_of_the_library() {
    let library = library();
    let bounds = library.groups.bounds().to_vec();
    let flux = [1.0, 1.0];
    let members: [(f64, &[f64]); 1] = [(1.0, &flux)];
    let run = |chain: &Chain| {
        let mut rates = chain.zero_rates();
        add_rates(chain, &[0], &library, T, &members, &bounds, &mut rates)
            .map_err(|e| e.to_string())
    };
    assert!(run(&chain()).is_ok());
    // チェーンの反応がライブラリにない
    let mut missing = chain();
    missing.nuclides[0]
        .reactions
        .push(chain_reaction(17, &[None]));
    assert!(run(&missing).unwrap_err().contains("MT17"));
    // 経路の状態がチェーンと違う
    let mut other = chain();
    other.nuclides[0].reactions[1].channels.pop();
    assert!(run(&other).unwrap_err().contains("channels of MT102"));
    // ライブラリで核分裂する核種が、チェーンで核分裂を持たない
    let mut stable = chain();
    stable.nuclides[0].fission = None;
    assert!(
        run(&stable)
            .unwrap_err()
            .contains("no fission in the chain")
    );
    // 対象の核種がライブラリの温度のデータを持たない
    let mut rates = chain().zero_rates();
    assert!(
        add_rates(
            &chain(),
            &[0],
            &library,
            900.0,
            &members,
            &bounds,
            &mut rates
        )
        .is_err()
    );
}
