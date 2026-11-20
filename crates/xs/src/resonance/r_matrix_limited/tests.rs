use super::*;
use crate::endf::{OrbitalL, ParticlePair, RMatrixLimited, RmlChannel, RmlResonance, SpinGroup};
use crate::physics::wave_number_factor;
use crate::resonance::reich_moore_tests::{AP, AWRI, single_level_exact};

fn pair(ma: f64, mb: f64, q: f64, mt: u32) -> ParticlePair {
    ParticlePair {
        ma,
        mb,
        za: 0.0,
        zb: 92.0,
        ia: if ma == 0.0 { 1.0 } else { 0.5 },
        ib: 0.0,
        q,
        pnt: 0,
        shf: 0,
        mt,
    }
}

fn channel(pair: usize, l: OrbitalL, sch: f64) -> RmlChannel {
    RmlChannel {
        pair,
        l,
        sch,
        bnd: 0.0,
        ape: AP,
        apt: AP,
        background: None,
    }
}

/// 捕獲と弾性散乱の粒子対に、`levels` の Q 値の非弾性散乱の粒子対（MT51 から）を加える。
fn pairs(levels: &[f64]) -> Vec<ParticlePair> {
    let mut pairs = vec![pair(0.0, AWRI + 1.0, 0.0, 102), pair(1.0, AWRI, 0.0, 2)];
    pairs.extend(
        levels
            .iter()
            .zip(51..)
            .map(|(&q, mt)| pair(1.0, AWRI, q, mt)),
    );
    pairs
}

/// Reich-Moore の `single_level` と同じ、l=0 の単一準位。
fn single_level(gn: f64, gg: f64) -> RMatrixLimited {
    RMatrixLimited {
        pairs: pairs(&[]),
        spin_groups: vec![SpinGroup {
            aj: 0.5,
            pj: 0.0,
            channels: vec![channel(0, OrbitalL::S, 0.0), channel(1, OrbitalL::S, 0.5)],
            resonances: vec![RmlResonance {
                er: 6.674,
                widths: vec![gg, gn],
            }],
        }],
    }
}

fn assert_close(actual: f64, expected: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1e-6),
        "{what}: actual {actual}, expected {expected}"
    );
}

#[test]
fn single_level_matches_breit_wigner() {
    let xs = RMatrixLimitedXs::new(&single_level(1.5e-3, 2.3e-2)).unwrap();
    assert_eq!(xs.reactions(), vec![1, 2, 102]);
    for e in [1.0e-3, 1.0, 6.0, 6.674, 6.69, 7.0, 100.0, 1.0e4] {
        let expected = single_level_exact(e, 1.5e-3, 2.3e-2, 0.0);
        let [total, elastic, capture] = xs.cross_sections(e).unwrap()[..] else {
            panic!("expected total, elastic and capture");
        };
        assert_close(total, expected.total, &format!("total at {e}"));
        assert_close(elastic, expected.elastic, &format!("elastic at {e}"));
        assert_close(capture, expected.capture, &format!("capture at {e}"));
    }
}

/// 弾性散乱と、Q = -1 keV の準位への非弾性散乱のチャネルを持つ共鳴（捕獲なし）。
pub(crate) fn inelastic_level(gg: f64) -> RMatrixLimited {
    RMatrixLimited {
        pairs: pairs(&[-1.0e3]),
        spin_groups: vec![SpinGroup {
            aj: 0.5,
            pj: 0.0,
            channels: vec![
                channel(0, OrbitalL::S, 0.0),
                channel(1, OrbitalL::S, 0.5),
                channel(2, OrbitalL::S, 0.5),
            ],
            resonances: vec![RmlResonance {
                er: 2.0e3,
                widths: vec![gg, 3.0, -2.0],
            }],
        }],
    }
}

#[test]
fn inelastic_channel_opens_at_threshold() {
    let xs = RMatrixLimitedXs::new(&inelastic_level(0.0)).unwrap();
    assert_eq!(xs.reactions(), vec![1, 2, 102, 51]);
    // 閾値（実験室系）は |Q| (A+1)/A
    let threshold = 1.0e3 * (AWRI + 1.0) / AWRI;
    let below = xs.cross_sections(threshold * 0.999).unwrap();
    assert_eq!(below[3], 0.0);
    for e in [threshold * 1.001, 1.5e3, 1.999e3, 2.0e3, 2.001e3, 3.0e3] {
        let [total, elastic, capture, level] = xs.cross_sections(e).unwrap()[..] else {
            panic!("expected four reactions");
        };
        assert!(level > 0.0, "level at {e}");
        // 捕獲がなければ、全断面積は弾性散乱と非弾性散乱の和（ユニタリ性）
        assert!(capture.abs() <= 1e-9 * total, "capture {capture} at {e}");
        assert_close(elastic + level, total, &format!("unitarity at {e}"));
    }
}

#[test]
fn capture_takes_the_rest_of_the_absorption() {
    let xs = RMatrixLimitedXs::new(&inelastic_level(0.5)).unwrap();
    let [total, elastic, capture, level] = xs.cross_sections(2.0e3).unwrap()[..] else {
        panic!("expected four reactions");
    };
    assert!(capture > 0.0);
    assert_close(elastic + capture + level, total, "sum of parts");
}

#[test]
fn pole_of_a_resonance_without_capture_is_the_limit() {
    // 捕獲がないと共鳴エネルギーは R 行列の極になるが、断面積は連続。両側の値の平均（傾きの分が
    // 打ち消し合う）と比べる
    let xs = RMatrixLimitedXs::new(&inelastic_level(0.0)).unwrap();
    let at = xs.cross_sections(2.0e3).unwrap();
    let below = xs.cross_sections(2.0e3 * (1.0 - 1.0e-9)).unwrap();
    let above = xs.cross_sections(2.0e3 * (1.0 + 1.0e-9)).unwrap();
    for ((a, b), c) in at.iter().zip(&below).zip(&above) {
        assert_close(*a, 0.5 * (b + c), "pole");
    }
}

#[test]
fn resonance_without_widths_does_not_count_as_a_pole() {
    // 幅がすべて 0 の共鳴は何も寄与しないので、同じエネルギーに極があっても問題ない
    let mut data = inelastic_level(0.0);
    data.spin_groups[0].resonances.push(RmlResonance {
        er: 2.0e3,
        widths: vec![0.0, 0.0, 0.0],
    });
    let xs = RMatrixLimitedXs::new(&data).unwrap();
    let single = RMatrixLimitedXs::new(&inelastic_level(0.0)).unwrap();
    assert_eq!(
        xs.cross_sections(2.0e3).unwrap(),
        single.cross_sections(2.0e3).unwrap()
    );
}

#[test]
fn rejects_width_of_a_channel_closed_at_the_resonance() {
    let mut data = inelastic_level(0.5);
    // 共鳴が閾値より下にあると、非弾性散乱の幅は定義できない
    data.spin_groups[0].resonances[0].er = 500.0;
    let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
    assert!(error.contains("closed"), "{error}");
}

/// 捕獲と中性子の 1 チャネル（軌道角運動量 `l`）の単一準位。中性子の粒子対のシフト因子の指定は `shf`、
/// 境界条件は `bnd`。
fn one_channel_level(l: OrbitalL, shf: i64, bnd: f64, gn: f64, gg: f64) -> RMatrixLimited {
    let mut data = single_level(gn, gg);
    data.pairs[1].shf = shf;
    let group = &mut data.spin_groups[0];
    group.aj = 0.5 + f64::from(l.number());
    group.channels[1].l = l;
    group.channels[1].bnd = bnd;
    data
}

/// 1 チャネルの R 行列 R(E) = γ²/(E_λ − E − iΓγ/2) + R_bkg(E) の、ENDF-102 の式 (D.42)〜(D.45) による
/// 全断面積、弾性散乱、捕獲。X = P R / (1 − R L)、L = (S − B) + iP、U = e^{−2iφ}(1 + 2iX)、
/// σ_t = 2πg/k² Re(1 − U)、σ_el = πg/k² |1 − U|²。
/// `shift` がなければ S − B = 0。
fn one_channel_exact(
    data: &RMatrixLimited,
    e: f64,
    shift_on: bool,
    background: impl Fn(f64) -> f64,
) -> [f64; 3] {
    use crate::resonance::factors::{penetrability, phase_shift, shift};
    let group = &data.spin_groups[0];
    let channel = &group.channels[1];
    let resonance = &group.resonances[0];
    let (gg, gn) = (resonance.widths[0], resonance.widths[1]);
    let k = |e: f64| wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
    let p = |e: f64| penetrability(channel.l, k(e) * channel.apt);
    let gamma2 = gn.abs() / (2.0 * p(resonance.er.abs()));
    let r =
        Complex64::new(gamma2, 0.0) / Complex64::new(resonance.er - e, -0.5 * gg) + background(e);
    let s_minus_b = if shift_on {
        shift(channel.l, k(e) * channel.apt) - channel.bnd
    } else {
        0.0
    };
    let l = Complex64::new(s_minus_b, p(e));
    let x = p(e) * r / (1.0 - r * l);
    let phi = phase_shift(channel.l, k(e) * channel.ape);
    // 1 − U = 2i sin φ e^{−iφ} − 2iX e^{−2iφ}（低いエネルギーで 1 − Re U の桁落ちを避ける）
    let one_minus_u = 2.0 * Complex64::I * phi.sin() * Complex64::from_polar(1.0, -phi)
        - 2.0 * Complex64::I * x * Complex64::from_polar(1.0, -2.0 * phi);
    let g = (2.0 * group.aj + 1.0) / 2.0;
    let pifac = PI * g / (k(e) * k(e));
    let total = 2.0 * pifac * one_minus_u.re;
    let elastic = pifac * one_minus_u.norm_sqr();
    // 1 − |U|² = 4(Im X − |X|²)
    let capture = 4.0 * pifac * (x.im - x.norm_sqr());
    [total, elastic, capture]
}

fn assert_matches_one_channel(data: &RMatrixLimited, shift_on: bool, bkg: impl Fn(f64) -> f64) {
    let xs = RMatrixLimitedXs::new(data).unwrap();
    for e in [
        1.0e-5, 1.0e-3, 1.0, 6.0, 6.674, 6.69, 7.0, 100.0, 1.0e4, 1.0e5,
    ] {
        let expected = one_channel_exact(data, e, shift_on, &bkg);
        let actual = xs.cross_sections(e).unwrap();
        assert_close(actual[0], expected[0], &format!("total at {e}"));
        assert_close(actual[1], expected[1], &format!("elastic at {e}"));
        // xs は捕獲を全断面積から引いて出すので、全断面積の丸め（数 ulp）の分だけずれうる
        let capture = (actual[2] - expected[2]).abs();
        assert!(
            capture <= 1e-9 * expected[2] + 8.0 * f64::EPSILON * expected[0],
            "capture at {e}: actual {}, expected {}",
            actual[2],
            expected[2]
        );
    }
}

#[test]
fn shift_factor_enters_the_level_matrix() {
    // l=1、B = −1（S_1(0) と同じ。ENDF/B-VIII.1 の Cu-63 の p 波の形）
    let data = one_channel_level(OrbitalL::P, 1, -1.0, 1.5e-3, 2.3e-2);
    assert_matches_one_channel(&data, true, |_| 0.0);
    // S − B が 0 でない所で、シフト因子を使わないときと違う
    let shifted = RMatrixLimitedXs::new(&data).unwrap();
    let plain =
        RMatrixLimitedXs::new(&one_channel_level(OrbitalL::P, 0, -1.0, 1.5e-3, 2.3e-2)).unwrap();
    let (a, b) = (
        shifted.cross_sections(1.0e4).unwrap(),
        plain.cross_sections(1.0e4).unwrap(),
    );
    assert!(
        (a[2] - b[2]).abs() > 1e-6 * b[2],
        "capture {} and {}",
        a[2],
        b[2]
    );
}

#[test]
fn shift_factor_of_an_s_wave_with_zero_boundary_changes_nothing() {
    let shifted = one_channel_level(OrbitalL::S, 1, 0.0, 1.5e-3, 2.3e-2);
    assert_matches_one_channel(&shifted, true, |_| 0.0);
}

#[test]
fn shift_factor_flag_minus_one_is_not_calculated() {
    // ENDF-102（2005 年版）は SHF=-1 を計算しない指定とする（ENDF/B-VII.1 の Cl-35）
    let data = one_channel_level(OrbitalL::P, -1, -1.0, 1.5e-3, 2.3e-2);
    assert_matches_one_channel(&data, false, |_| 0.0);
}

#[test]
fn pole_with_shift_factor_is_the_limit() {
    let data = one_channel_level(OrbitalL::P, 1, -1.0, 1.5e-3, 0.0);
    let xs = RMatrixLimitedXs::new(&data).unwrap();
    let er = 6.674;
    let at = xs.cross_sections(er).unwrap();
    let below = xs.cross_sections(er * (1.0 - 1.0e-9)).unwrap();
    let above = xs.cross_sections(er * (1.0 + 1.0e-9)).unwrap();
    for ((a, b), c) in at.iter().zip(&below).zip(&above) {
        assert_close(*a, 0.5 * (b + c), "pole");
    }
}

#[test]
fn rejects_unknown_shift_factor_flag() {
    let data = one_channel_level(OrbitalL::P, 2, -1.0, 1.5e-3, 2.3e-2);
    let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
    assert!(error.contains("SHF=2"), "{error}");
}

#[test]
fn rejects_shift_factor_for_capture() {
    let mut data = single_level(1.5e-3, 2.3e-2);
    data.pairs[0].shf = 1;
    let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
    assert!(error.contains("SHF=1"), "{error}");
}

/// ENDF/B-VIII.1 の Sr-88 の形の、SAMMY の対数のバックグラウンド（LBK=2）。
fn sammy_background() -> crate::endf::LogBackground {
    crate::endf::LogBackground {
        ed: 0.0,
        eu: 9.55e5,
        r0: -4.3e-2,
        r1: 2.8e-8,
        r2: 1.0e-14,
        s0: 1.0e-2,
        s1: 3.0e-9,
    }
}

/// ENDF-102 の式 (D.77)。
fn sammy_background_at(e: f64) -> f64 {
    let b = sammy_background();
    b.r0 + b.r1 * e + b.r2 * e * e
        - b.s1 * (b.eu - b.ed)
        - (b.s0 + b.s1 * e) * ((b.eu - e) / (e - b.ed)).ln()
}

#[test]
fn background_r_matrix_adds_to_the_channel() {
    let mut data = one_channel_level(OrbitalL::S, 0, 0.0, 1.5e-3, 2.3e-2);
    data.spin_groups[0].channels[1].background = Some(sammy_background());
    assert_matches_one_channel(&data, false, sammy_background_at);
}

#[test]
fn background_r_matrix_with_shift_factor() {
    let mut data = one_channel_level(OrbitalL::P, 1, -1.0, 1.5e-3, 2.3e-2);
    data.spin_groups[0].channels[1].background = Some(sammy_background());
    assert_matches_one_channel(&data, true, sammy_background_at);
}

#[test]
fn rejects_energies_outside_the_background_range() {
    let mut data = one_channel_level(OrbitalL::S, 0, 0.0, 1.5e-3, 2.3e-2);
    data.spin_groups[0].channels[1].background = Some(crate::endf::LogBackground {
        eu: 1.0e3,
        ..sammy_background()
    });
    let xs = RMatrixLimitedXs::new(&data).unwrap();
    assert!(xs.cross_sections(10.0).is_ok());
    let error = xs.cross_sections(2.0e3).unwrap_err().to_string();
    assert!(error.contains("background"), "{error}");
}

#[test]
fn rejects_background_on_the_capture_channel() {
    let mut data = single_level(1.5e-3, 2.3e-2);
    data.spin_groups[0].channels[0].background = Some(sammy_background());
    let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
    assert!(error.contains("background"), "{error}");
}

#[test]
fn rejects_charged_pairs() {
    let mut data = single_level(1.5e-3, 2.3e-2);
    data.pairs[1].za = 1.0;
    assert!(RMatrixLimitedXs::new(&data).is_err());
}

/// 核分裂の粒子対（MT=18、透過因子を計算しない PNT=-1）。
fn fission_pair() -> ParticlePair {
    ParticlePair {
        pnt: -1,
        ..pair(1.0, AWRI, 0.0, 18)
    }
}

/// Reich-Moore の `single_level` と同じ l=0 の単一準位に、核分裂のチャネルを 1 つ加えたもの。
pub(crate) fn fission_level(gn: f64, gg: f64, gf: f64) -> RMatrixLimited {
    let mut data = single_level(gn, gg);
    data.pairs.push(fission_pair());
    let group = &mut data.spin_groups[0];
    group.channels.push(channel(2, OrbitalL::S, 0.0));
    group.resonances[0].widths.push(gf);
    data
}

#[test]
fn fission_channel_matches_single_level_breit_wigner() {
    let (gn, gg, gf) = (1.5e-3, 2.3e-2, 4.0e-2);
    let xs = RMatrixLimitedXs::new(&fission_level(gn, gg, gf)).unwrap();
    assert_eq!(xs.reactions(), vec![1, 2, 102, 18]);
    for e in [1.0e-3, 1.0, 6.0, 6.674, 6.69, 7.0, 100.0, 1.0e4] {
        let expected = single_level_exact(e, gn, gg, gf);
        let [total, elastic, capture, fission] = xs.cross_sections(e).unwrap()[..] else {
            panic!("expected total, elastic, capture and fission");
        };
        assert_close(total, expected.total, &format!("total at {e}"));
        assert_close(elastic, expected.elastic, &format!("elastic at {e}"));
        assert_close(capture, expected.capture, &format!("capture at {e}"));
        assert_close(fission, expected.fission, &format!("fission at {e}"));
    }
}

#[test]
fn two_fission_channels_match_reich_moore() {
    use crate::endf::{ReichMoore, RmLValue, RmResonance};
    use crate::resonance::ReichMooreXs;
    // 負の共鳴、核分裂の幅の符号が違う共鳴、核分裂の幅のない共鳴を含め、準位とチャネルの間の干渉を比べる
    // （gn, gg, gfa, gfb）
    let levels = [
        (-3.0, 2.0e-3, 4.0e-2, 0.1, -5.0e-2),
        (0.3, 1.0e-4, 4.0e-2, 6.0e-2, 2.0e-2),
        (7.8, 1.5e-3, 3.5e-2, -0.3, 0.1),
        (11.0, 2.0e-3, 4.0e-2, 0.0, 0.0),
    ];
    let rm = ReichMoore {
        spi: 0.0,
        ap: AP,
        naps: 1,
        l_values: vec![RmLValue {
            awri: AWRI,
            apl: 0.0,
            l: OrbitalL::S,
            resonances: levels
                .iter()
                .map(|&(er, gn, gg, gfa, gfb)| RmResonance {
                    er,
                    aj: 0.5,
                    gn,
                    gg,
                    gfa,
                    gfb,
                })
                .collect(),
        }],
    };
    let mut pairs = pairs(&[]);
    pairs.push(fission_pair());
    let rml = RMatrixLimited {
        pairs,
        spin_groups: vec![SpinGroup {
            aj: 0.5,
            pj: 0.0,
            channels: vec![
                channel(0, OrbitalL::S, 0.0),
                channel(1, OrbitalL::S, 0.5),
                channel(2, OrbitalL::S, 0.0),
                channel(2, OrbitalL::S, 0.0),
            ],
            resonances: levels
                .iter()
                .map(|&(er, gn, gg, gfa, gfb)| RmlResonance {
                    er,
                    widths: vec![gg, gn, gfa, gfb],
                })
                .collect(),
        }],
    };
    let reich_moore = ReichMooreXs::new(&rm);
    let xs = RMatrixLimitedXs::new(&rml).unwrap();
    assert_eq!(xs.reactions(), vec![1, 2, 102, 18]);
    for e in [1.0e-5, 0.0253, 0.3, 1.0, 5.0, 7.8, 9.0, 11.0, 100.0] {
        let expected = reich_moore.cross_sections(e).unwrap();
        let [total, elastic, capture, fission] = xs.cross_sections(e).unwrap()[..] else {
            panic!("expected total, elastic, capture and fission");
        };
        assert_close(total, expected.total, &format!("total at {e}"));
        assert_close(elastic, expected.elastic, &format!("elastic at {e}"));
        assert_close(capture, expected.capture, &format!("capture at {e}"));
        assert_close(fission, expected.fission, &format!("fission at {e}"));
    }
}

#[test]
fn rejects_penetrability_for_fission() {
    let mut data = fission_level(1.5e-3, 2.3e-2, 4.0e-2);
    data.pairs[2].pnt = 1;
    let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
    assert!(error.contains("PNT=1"), "{error}");
}

#[test]
fn rejects_two_fission_pairs() {
    let mut data = fission_level(1.5e-3, 2.3e-2, 4.0e-2);
    data.pairs.push(fission_pair());
    let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
    assert!(error.contains("MT18"), "{error}");
}

#[test]
fn peaks_of_shifted_levels_sit_at_the_observed_energy() {
    // p 波の幅の広い準位（ρ ≈ 0.55）は、シフト因子で −γ²(S − B) ≈ −4.5 eV だけ動き、半値幅（約 3 eV）を超える
    let mut data = one_channel_level(OrbitalL::P, 1, -1.0, 5.0, 1.0);
    data.spin_groups[0].resonances[0].er = 7.0e4;
    let xs = RMatrixLimitedXs::new(&data).unwrap();
    let [(peak, half_width)] = xs.peaks().collect::<Vec<_>>()[..] else {
        panic!("expected one peak");
    };
    let capture = |e: f64| xs.cross_sections(e).unwrap()[2];
    let maximum = (-4000..=4000)
        .map(|i| 7.0e4 + f64::from(i) * 5.0e-3)
        .max_by(|a, b| capture(*a).total_cmp(&capture(*b)))
        .unwrap();
    assert!(
        (maximum - 7.0e4).abs() > half_width,
        "the level moves by {} eV",
        maximum - 7.0e4
    );
    assert!(
        (peak - maximum).abs() < 0.1 * half_width,
        "peak {peak} eV, capture maximum {maximum} eV, half width {half_width} eV"
    );
}

/// 陽子の質量（中性子の質量の単位）
const PROTON: f64 = 0.998_623_5;

/// `single_level` の準位に、陽子の放出（MT600、陽子と電荷 16 の残留核、Q 値 `q`）の l = 0 のチャネル（幅 `gp`）を加える。
fn proton_level(q: f64, gp: f64) -> RMatrixLimited {
    let mut data = single_level(1.5e-3, 2.3e-2);
    data.pairs.push(ParticlePair {
        za: 1.0,
        zb: 16.0,
        pnt: 1,
        ..pair(PROTON, AWRI, q, 600)
    });
    let group = &mut data.spin_groups[0];
    group.channels.push(channel(2, OrbitalL::S, 0.5));
    group.resonances[0].widths.push(gp);
    data
}

/// 陽子のチャネルの Coulomb の透過因子。η と ρ は重心系のエネルギー E·A/(A+1) + Q から、
/// η = zZα √(μc²/2E_cm)、ρ = √(2μc² E_cm) a / ħc で求める（xs の波数 k からの ηk/k と別の道筋）。
fn proton_penetrability(e: f64, q: f64) -> f64 {
    const ALPHA: f64 = 7.297_352_569_3e-3;
    const HBAR_C: f64 = 1.973_269_804e-5; // eV cm
    let mu_c2 = PROTON * AWRI / (PROTON + AWRI) * crate::physics::AMASSN * 931.494_102_42e6;
    let e_cm = e * AWRI / (AWRI + 1.0) + q;
    let eta = 16.0 * ALPHA * (mu_c2 / (2.0 * e_cm)).sqrt();
    let rho = (2.0 * mu_c2 * e_cm).sqrt() / HBAR_C * AP * 1.0e-12;
    crate::resonance::coulomb::penetrability(0, eta, rho)
}

#[test]
fn proton_channel_matches_single_level_breit_wigner_with_coulomb_penetrability() {
    let (gn, gg, gp, q) = (1.5e-3, 2.3e-2, 1.0e-2, 6.0e5);
    let xs = RMatrixLimitedXs::new(&proton_level(q, gp)).unwrap();
    assert_eq!(xs.reactions(), vec![1, 2, 102, 600]);
    let er: f64 = 6.674;
    for e in [1.0e-3_f64, 1.0, 6.0, 6.674, 6.69, 7.0, 100.0, 1.0e4] {
        // 幅は透過因子に比例する。中性子の l = 0 の透過因子は ρ ∝ √E
        let gn_e = gn * (e / er).sqrt();
        let gp_e = gp * proton_penetrability(e, q) / proton_penetrability(er, q);
        let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
        let width = gn_e + gg + gp_e;
        let expected = PI / (k * k) * gn_e * gp_e / ((e - er).powi(2) + width * width / 4.0);
        let proton = xs.cross_sections(e).unwrap()[3];
        assert_close(proton, expected, &format!("(n,p) at {e}"));
    }
}

#[test]
fn endothermic_proton_channel_is_closed_below_its_threshold() {
    // Q = −1 MeV の閾値（実験室系）は 1 MeV × (A+1)/A。共鳴は閾値の上の 2 MeV に置く
    let mut data = proton_level(-1.0e6, 1.0e-2);
    data.spin_groups[0].resonances[0].er = 2.0e6;
    let xs = RMatrixLimitedXs::new(&data).unwrap();
    let threshold = 1.0e6 * (AWRI + 1.0) / AWRI;
    assert_eq!(xs.cross_sections(threshold * 0.999).unwrap()[3], 0.0);
    assert!(xs.cross_sections(2.0e6).unwrap()[3] > 0.0);
}

#[test]
fn rejects_widths_whose_coulomb_penetrability_underflows() {
    // 共鳴で陽子の重心系のエネルギーが 1 meV になり、η ≈ 8e4 の障壁の下で透過因子が 0 になる
    let q = 1.0e-3 - 6.674 * AWRI / (AWRI + 1.0);
    let error = RMatrixLimitedXs::new(&proton_level(q, 1.0e-2))
        .unwrap_err()
        .to_string();
    assert!(error.contains("underflows"), "{error}");
}

#[test]
fn rejects_shift_factor_for_charged_pairs() {
    let mut data = proton_level(6.0e5, 1.0e-2);
    data.pairs[2].shf = 1;
    let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
    assert!(error.contains("SHF=1"), "{error}");
}

#[test]
fn rejects_charged_pairs_outside_charged_particle_levels() {
    let mut data = proton_level(6.0e5, 1.0e-2);
    data.pairs[2].mt = 103;
    let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
    assert!(error.contains("MT=103"), "{error}");
}

#[test]
fn rejects_charged_pairs_without_penetrability() {
    let mut data = proton_level(6.0e5, 1.0e-2);
    data.pairs[2].pnt = -1;
    let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
    assert!(error.contains("PNT=-1"), "{error}");
}

#[test]
fn lists_widths_whose_coulomb_penetrability_njoy_approximates() {
    // Q = 0.6 MeV の陽子は η ≈ 3.3 で、NJOY の近似の式の外。Q = 50 keV では η ≈ 11 で、障壁の下（転回点 2η ≫ ρ ≈ 0.5）
    let above = RMatrixLimitedXs::new(&proton_level(6.0e5, 1.0e-2)).unwrap();
    assert!(above.njoy_approximated_widths().is_empty());
    let below = RMatrixLimitedXs::new(&proton_level(5.0e4, 1.0e-2)).unwrap();
    assert_eq!(below.njoy_approximated_widths(), vec![600]);
}
