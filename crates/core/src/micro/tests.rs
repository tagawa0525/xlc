use xlc_depletion::chain::{
    Channel, Decay, DecayEnergy, DecayMode, Fission, Nuclide, NuclideRates, Product, Reaction,
    Yield,
};
use xlc_depletion::{CramOrder, deplete};
use xlc_format::table::MicroReaction;

use super::*;

/// 1 eV (J)
const EV: f64 = 1.602_176_634e-19;
/// 1 barn (cm²)
const BARN: f64 = 1.0e-24;

fn nuclide(name: &str, z: u32, a: u32) -> Nuclide {
    Nuclide {
        name: name.into(),
        z,
        a,
        state: 0,
        excitation: 0.0,
        decay: None,
        reactions: Vec::new(),
        fission: None,
    }
}

fn product(nuclide: usize) -> Product {
    Product { nuclide, count: 1 }
}

/// β- 崩壊で `daughter` になる
fn beta(half_life: f64, daughter: usize) -> Option<Decay> {
    Some(Decay {
        half_life,
        energy: DecayEnergy {
            light: 0.0,
            electromagnetic: 0.0,
            heavy: 0.0,
        },
        modes: vec![DecayMode {
            mode: "beta-".into(),
            ratio: 1.0,
            products: vec![product(daughter)],
        }],
    })
}

/// 捕獲で `to` になる（状態に分けない 1 つの経路）
fn capture(to: usize) -> Reaction {
    Reaction {
        mt: 102,
        channels: vec![Channel {
            final_state: None,
            products: vec![product(to)],
        }],
    }
}

// チェーンの添字
const SB135: usize = 0;
const TE135: usize = 1;
const I135: usize = 2;
const XE135: usize = 3;
const XE136: usize = 4;
const CS135: usize = 5;
const BA140: usize = 6;
const U235: usize = 7;
const U236: usize = 8;
const U237: usize = 9;
const NP237: usize = 10;

/// 小さなチェーン: U-235 の核分裂は Te-135、I-135、Xe-135、Ba-140 を生み（熱と高速の 2 点の収率）、Te-135 → I-135 → Xe-135 → Cs-135 と
/// 崩壊する。Sb-135 は Te-135 へ崩壊するが、どこからも生まれない。U-235 の捕獲は U-236、U-236 の捕獲は U-237 で、U-237 は Np-237 へ
/// 崩壊する。Xe-135 の捕獲は Xe-136。I-135 と Np-237 の捕獲の生成物はチェーンにない
fn chain() -> Chain {
    let mut nuclides = vec![
        nuclide("Sb-135", 51, 135),
        nuclide("Te-135", 52, 135),
        nuclide("I-135", 53, 135),
        nuclide("Xe-135", 54, 135),
        nuclide("Xe-136", 54, 136),
        nuclide("Cs-135", 55, 135),
        nuclide("Ba-140", 56, 140),
        nuclide("U-235", 92, 235),
        nuclide("U-236", 92, 236),
        nuclide("U-237", 92, 237),
        nuclide("Np-237", 93, 237),
    ];
    nuclides[SB135].decay = beta(1.7, TE135);
    nuclides[TE135].decay = beta(19.0, I135);
    nuclides[I135].decay = beta(6.57 * 3600.0, XE135);
    nuclides[XE135].decay = beta(9.14 * 3600.0, CS135);
    nuclides[XE135].reactions = vec![capture(XE136)];
    nuclides[U235].reactions = vec![capture(U236)];
    nuclides[U236].reactions = vec![capture(U237)];
    nuclides[U237].decay = beta(6.75 * 86400.0, NP237);
    // 生成物がチェーンにない捕獲
    let lost = Reaction {
        mt: 102,
        channels: vec![Channel {
            final_state: None,
            products: Vec::new(),
        }],
    };
    nuclides[I135].reactions = vec![lost.clone()];
    nuclides[NP237].reactions = vec![lost];
    let yields = |scale: f64| {
        vec![
            Yield {
                nuclide: TE135,
                value: 0.03 * scale,
            },
            Yield {
                nuclide: I135,
                value: 0.03,
            },
            Yield {
                nuclide: XE135,
                value: 0.002,
            },
            Yield {
                nuclide: BA140,
                value: 0.06,
            },
        ]
    };
    nuclides[U235].fission = Some(Fission {
        q: 193.0e6,
        yield_source: "U-235".into(),
        energies: vec![0.0253, 5.0e5],
        yields: vec![yields(1.0), yields(0.8)],
    });
    let chain = Chain { nuclides };
    chain.check().unwrap();
    chain
}

/// 2 群の微視的断面積（捕獲の 1 経路。`fission` が 0 なら核分裂しない）
fn sigma(capture: [f64; 2], fission: [f64; 2], thermal: f64) -> MicroNuclide {
    let fissile = fission.iter().any(|&f| f > 0.0);
    MicroNuclide {
        total: vec![
            capture[0] + fission[0] + 4.0,
            capture[1] + fission[1] + 10.0,
        ],
        nu_fission: fission.iter().map(|f| 2.43 * f).collect(),
        fission: fission.to_vec(),
        scatter: vec![vec![3.9, 0.1], vec![0.0, 10.0]],
        reactions: vec![MicroReaction {
            mt: 102,
            channels: vec![capture.to_vec()],
        }],
        thermal_fission: fissile.then_some(thermal),
    }
}

const TRACKED: [&str; 5] = ["U-235", "U-236", "Np-237", "Xe-135", "I-135"];

/// 追う核種（`TRACKED` の順）の微視的断面積
fn micro_at(densities: Vec<f64>) -> MicroAt {
    MicroAt {
        densities,
        nuclides: vec![
            sigma([10.0, 100.0], [1.5, 500.0], 0.9),
            sigma([0.5, 5.0], [0.0, 0.0], 0.0),
            sigma([1.0, 150.0], [0.0, 0.0], 0.0),
            sigma([0.0, 2.0e6], [0.0, 0.0], 0.0),
            sigma([0.0, 7.0], [0.0, 0.0], 0.0),
        ],
    }
}

fn names(n: &[&str]) -> Vec<String> {
    n.iter().map(|s| s.to_string()).collect()
}

fn tracked_chain() -> MicroChain {
    MicroChain::new(&chain(), &names(&TRACKED), &micro_at(vec![0.0; 5]).nuclides).unwrap()
}

#[test]
fn auxiliary_nuclides_are_the_reached_untracked_precursors_of_tracked_ones() {
    let micro = tracked_chain();
    assert_eq!(micro.nuclides, names(&TRACKED));
    // Te-135（収率で生まれ I-135 へ）と U-237（捕獲で生まれ Np-237 へ）。Sb-135 は生まれず、Xe-136、Cs-135、Ba-140 は追う核種へ流れ込まない
    assert_eq!(
        micro.solved,
        names(&[
            "Te-135", "I-135", "Xe-135", "U-235", "U-236", "U-237", "Np-237"
        ])
    );
    assert_eq!(micro.auxiliary(), 2);
    let n = micro.initial(&[1.0, 2.0, 3.0, 4.0, 5.0]);
    assert_eq!(n, [0.0, 5.0, 4.0, 1.0, 2.0, 0.0, 3.0]);
    assert_eq!(micro.tracked(&n), [1.0, 2.0, 3.0, 4.0, 5.0]);
}

#[test]
fn tracked_nuclides_must_be_in_the_chain_with_its_reactions_and_fission() {
    let reference = micro_at(vec![0.0; 5]).nuclides;
    // チェーンにない核種
    let mut missing = names(&TRACKED);
    missing[2] = "Pu-239".into();
    assert!(MicroChain::new(&chain(), &missing, &reference).is_err());
    // テーブルの反応の MT がチェーンと違う
    let mut other = reference.clone();
    other[1].reactions[0].mt = 16;
    assert!(MicroChain::new(&chain(), &names(&TRACKED), &other).is_err());
    // 経路の数が違う
    let mut other = reference.clone();
    other[0].reactions[0].channels.push(vec![0.0, 0.0]);
    assert!(MicroChain::new(&chain(), &names(&TRACKED), &other).is_err());
    // テーブルは核分裂するが、チェーンに収率がない
    let mut other = reference.clone();
    other[2] = sigma([1.0, 150.0], [0.1, 0.0], 0.0);
    assert!(MicroChain::new(&chain(), &names(&TRACKED), &other).is_err());
    // 同じ核種を 2 度挙げる
    let mut twice = names(&TRACKED);
    twice[4] = "U-235".into();
    assert!(MicroChain::new(&chain(), &twice, &reference).is_err());
}

#[test]
fn rates_make_the_tracked_fission_power_the_node_power_density() {
    let micro = tracked_chain();
    let at = micro_at(vec![0.0; 5]);
    let n = micro.initial(&[7.0e-4, 1.0e-5, 1.0e-6, 1.0e-9, 1.0e-8]);
    let flux = [3.0, 1.0];
    let power = 104.5;
    let rates = micro.rates(&at, &flux, &n, power).unwrap();
    // 規格化の倍率: 追う核分裂性の核種は U-235 だけ
    let fission_shape = 1.5 * 3.0 + 500.0 * 1.0;
    let scale = power / (7.0e-4 * fission_shape * 193.0e6 * EV);
    let close = |a: f64, b: f64| (a - b).abs() <= 1e-14 * b.abs();
    assert!(close(rates[0].fission, fission_shape * scale * BARN));
    assert!(close(
        rates[0].reactions[0][0],
        (10.0 * 3.0 + 100.0 * 1.0) * scale * BARN
    ));
    assert!(close(rates[3].reactions[0][0], 2.0e6 * 1.0 * scale * BARN));
    assert_eq!(rates[1].fission, 0.0);
    assert!(rates[1].yield_weights.is_empty());
    // 熱の点に 0.9、高速の点に 0.1
    assert_eq!(rates[0].yield_weights, [0.9, 1.0 - 0.9]);
    // 出力: Σ N σf φ Q = 出力密度（数密度 1/(barn cm) と BARN で 1/cm³）
    let fission_power = 7.0e-4 / BARN * rates[0].fission * 193.0e6 * EV;
    assert!(close(fission_power, power), "{fission_power}");
    // 追う核分裂性の核種の数密度がなければ、出力を出せない
    assert!(
        micro
            .rates(&at, &flux, &micro.initial(&[0.0; 5]), power)
            .is_err()
    );
}

#[test]
fn solved_nuclides_deplete_as_in_the_whole_chain() {
    let micro = tracked_chain();
    let at = micro_at(vec![0.0; 5]);
    let tracked = [7.0e-4, 1.0e-5, 1.0e-6, 1.0e-9, 1.0e-8];
    let n = micro.initial(&tracked);
    let rates = micro.rates(&at, &[3.0e13, 1.0e13], &n, 104.5).unwrap();
    let a = micro.matrix(&rates).unwrap();
    let t = 30.0 * 86400.0;
    let got = deplete(&a, &n, t, CramOrder::Order48).unwrap();
    // チェーンの全核種で、同じ反応率を追う核種に与えて解く
    let whole = chain();
    let mut all = whole.zero_rates();
    for (name, r) in TRACKED.iter().zip(&rates) {
        let i = whole.nuclides.iter().position(|x| x.name == *name).unwrap();
        all[i] = NuclideRates {
            reactions: r.reactions.clone(),
            fission: r.fission,
            yield_weights: if r.yield_weights.is_empty() {
                all[i].yield_weights.clone()
            } else {
                r.yield_weights.clone()
            },
        };
    }
    let mut n0 = vec![0.0; whole.nuclides.len()];
    for (name, &v) in TRACKED.iter().zip(&tracked) {
        n0[whole.nuclides.iter().position(|x| x.name == *name).unwrap()] = v;
    }
    let want = deplete(&whole.matrix(&all).unwrap(), &n0, t, CramOrder::Order48).unwrap();
    for (name, g) in micro.solved.iter().zip(&got) {
        let w = want[whole.nuclides.iter().position(|x| &x.name == name).unwrap()];
        assert!(
            (g - w).abs() <= 1e-12 * w.abs().max(1e-30),
            "{name}: {g} vs {w}"
        );
    }
}

#[test]
fn yield_weights_put_the_thermal_fraction_on_the_cutoff_points() {
    // 熱と高速の点があれば、その 2 点
    assert_eq!(
        yield_weights(&[0.0253, 5.0e5, 1.4e7], 0.75),
        [0.75, 0.25, 0.0]
    );
    // 熱の点がなければ、境界（112.47 eV）より下で 0.0253 eV に最も近い点
    assert_eq!(yield_weights(&[0.1, 5.0e5, 1.4e7], 0.75), [0.75, 0.25, 0.0]);
    // 高速の点がなければ、境界以上で 5e5 eV に最も近い点
    assert_eq!(
        yield_weights(&[0.0253, 4.0e5, 2.0e6], 0.75),
        [0.75, 0.25, 0.0]
    );
    // 境界が最初の点以下なら最初の点、最後の点以上なら最後の点だけ
    assert_eq!(yield_weights(&[4.0e5, 2.0e6], 0.75), [1.0, 0.0]);
    assert_eq!(yield_weights(&[0.0253, 1.0], 0.75), [0.0, 1.0]);
}

/// 2 群の少数群定数
fn constants() -> FewGroup {
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: 1.2,
        total: vec![0.52, 1.3],
        absorption: vec![0.01, 0.08],
        nu_fission: vec![0.006, 0.12],
        fission: vec![0.0025, 0.05],
        chi: vec![1.0, 0.0],
        diffusion: vec![1.4, 0.4],
        scatter: vec![vec![0.49, 0.02], vec![0.001, 1.22]],
        adf: vec![vec![1.0, 1.02]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: None,
    }
}

#[test]
fn corrections_add_the_tracked_density_differences_to_the_reaction_items() {
    let nominal = vec![7.0e-4, 1.0e-5, 1.0e-6, 1.0e-9, 1.0e-8];
    let at = micro_at(nominal.clone());
    // 公称の数密度なら、補正は 0 で定数はそのまま
    let mut f = constants();
    correct(&mut f, &at, &nominal);
    assert_eq!(f, constants());
    let n = [6.0e-4, 2.0e-5, 1.0e-6, 3.0e-9, 1.0e-8];
    let mut f = constants();
    correct(&mut f, &at, &n);
    let base = constants();
    let delta: Vec<f64> = n.iter().zip(&nominal).map(|(a, b)| a - b).collect();
    let sum = |pick: &dyn Fn(&MicroNuclide) -> f64| -> f64 {
        at.nuclides
            .iter()
            .zip(&delta)
            .map(|(s, d)| pick(s) * d)
            .sum()
    };
    let close = |a: f64, b: f64| (a - b).abs() <= 1e-15 * b.abs().max(1.0);
    for g in 0..2 {
        assert!(close(f.total[g], base.total[g] + sum(&|s| s.total[g])));
        assert!(close(
            f.nu_fission[g],
            base.nu_fission[g] + sum(&|s| s.nu_fission[g])
        ));
        assert!(close(
            f.fission[g],
            base.fission[g] + sum(&|s| s.fission[g])
        ));
        // 吸収は全断面積から散乱の行の和を引いたもの
        let absorption = sum(&|s| s.total[g] - s.scatter[g].iter().sum::<f64>());
        assert!(close(f.absorption[g], base.absorption[g] + absorption));
        for h in 0..2 {
            assert!(close(
                f.scatter[g][h],
                base.scatter[g][h] + sum(&|s| s.scatter[g][h])
            ));
        }
    }
    // χ、拡散係数、ADF、k∞ は補正しない
    assert_eq!(
        (f.chi, f.diffusion, f.adf, f.k_inf),
        (base.chi, base.diffusion, base.adf, base.k_inf)
    );
}

#[test]
fn invalid_rates_are_errors_in_the_matrix() {
    let micro = tracked_chain();
    let at = micro_at(vec![0.0; 5]);
    let n = micro.initial(&[7.0e-4, 1.0e-5, 1.0e-6, 1.0e-9, 1.0e-8]);
    let rates = micro.rates(&at, &[3.0, 1.0], &n, 104.5).unwrap();
    for bad in [f64::NAN, -1.0e-9, f64::INFINITY] {
        let mut r = rates.clone();
        r[1].reactions[0][0] = bad;
        assert!(micro.matrix(&r).is_err(), "{bad}");
        let mut r = rates.clone();
        r[0].fission = bad;
        assert!(micro.matrix(&r).is_err(), "{bad}");
    }
}
