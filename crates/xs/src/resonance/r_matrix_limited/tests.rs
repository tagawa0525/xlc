use super::*;
use crate::endf::{OrbitalL, ParticlePair, RMatrixLimited, RmlChannel, RmlResonance, SpinGroup};
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

#[test]
fn rejects_shift_factor() {
    let mut data = single_level(1.5e-3, 2.3e-2);
    data.pairs[1].shf = 1;
    let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
    assert!(error.contains("SHF=1"), "{error}");
}

#[test]
fn rejects_charged_pairs() {
    let mut data = single_level(1.5e-3, 2.3e-2);
    data.pairs[1].za = 1.0;
    assert!(RMatrixLimitedXs::new(&data).is_err());
}
