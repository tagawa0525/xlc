use super::*;
use crate::endf::{RmLValue, RmResonance};
use crate::physics::wave_number_factor;
use std::f64::consts::PI;

pub(crate) const AWRI: f64 = 236.006;
pub(crate) const AP: f64 = 0.948;

pub(crate) fn single_level(gn: f64, gg: f64, gf: f64) -> ReichMoore {
    ReichMoore {
        spi: 0.0,
        ap: AP,
        naps: 1,
        l_values: vec![RmLValue {
            awri: AWRI,
            apl: 0.0,
            l: OrbitalL::S,
            resonances: vec![RmResonance {
                er: 6.674,
                aj: 0.5,
                gn,
                gg,
                gfa: gf,
                gfb: 0.0,
            }],
        }],
    }
}

/// スピン 0 の標的、l=0 の単一準位での厳密解（Breit-Wigner 型）。
pub(crate) fn single_level_exact(e: f64, gn: f64, gg: f64, gf: f64) -> Sigma {
    let er: f64 = 6.674;
    let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
    let kr = wave_number_factor() * AWRI / (AWRI + 1.0) * er.sqrt();
    let pifac = PI / (k * k);
    let phi = k * AP;
    // l=0 の透過因子は ρ なので、中性子幅は ρ に比例する
    let gn_e = gn * (k * AP) / (kr * AP);
    let width = gn_e + gg + gf;
    let d = er - e;
    // U = e^{-2iφ} (1 + iΓn / (d - iΓ/2))
    let denom = d * d + width * width / 4.0;
    let (qr, qi) = (-gn_e * width / 2.0 / denom, gn_e * d / denom);
    let (wr, wi) = (1.0 + qr, qi);
    let (c, s) = ((2.0 * phi).cos(), -(2.0 * phi).sin());
    let ui = c * wi + s * wr;
    // 1 - Re U を差を取らずに計算する（低エネルギーでの桁落ちを避ける）
    let one_minus_ur = 2.0 * phi.sin().powi(2) - (2.0 * phi).cos() * qr - (2.0 * phi).sin() * qi;
    let total = 2.0 * pifac * one_minus_ur;
    let elastic = pifac * (one_minus_ur.powi(2) + ui * ui);
    let fission = pifac * gn_e * gf / denom;
    Sigma {
        total,
        elastic,
        fission,
        capture: total - elastic - fission,
    }
}

fn assert_close(actual: f64, expected: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1e-6),
        "{what}: actual {actual}, expected {expected}"
    );
}

fn assert_sigma(actual: Sigma, expected: Sigma, e: f64) {
    assert_close(actual.total, expected.total, &format!("total at {e}"));
    assert_close(actual.elastic, expected.elastic, &format!("elastic at {e}"));
    assert_close(actual.fission, expected.fission, &format!("fission at {e}"));
    assert_close(actual.capture, expected.capture, &format!("capture at {e}"));
}

#[test]
fn singular_level_matrix_is_an_error() {
    // M = I + K が特異になる K。物理的な共鳴パラメータ（K の実部が半正定値）では起きない
    let mut group = ChannelGroup::default();
    group.r[0][0] = -1.0;
    assert!(group.terms(0.5, 0.1).is_err());
}

#[test]
fn single_level_without_fission_matches_breit_wigner() {
    let xs = ReichMooreXs::new(&single_level(1.5e-3, 2.3e-2, 0.0));
    for e in [1.0e-3, 1.0, 6.0, 6.674, 6.69, 7.0, 100.0, 1.0e4] {
        assert_sigma(
            xs.cross_sections(e).unwrap(),
            single_level_exact(e, 1.5e-3, 2.3e-2, 0.0),
            e,
        );
    }
}

#[test]
fn single_level_with_fission_matches_breit_wigner() {
    let xs = ReichMooreXs::new(&single_level(1.5e-3, 2.3e-2, 1.0e-1));
    for e in [1.0e-3, 1.0, 6.674, 6.8, 100.0] {
        assert_sigma(
            xs.cross_sections(e).unwrap(),
            single_level_exact(e, 1.5e-3, 2.3e-2, 1.0e-1),
            e,
        );
    }
}

#[test]
fn resonance_without_capture_width_is_the_limit_at_its_energy() {
    // Γγ = 0 では E = Er で R 行列が発散するが、断面積は Breit-Wigner の値に連続につながる
    let e = 6.674;
    for gf in [0.0, 1.0e-1] {
        let sigma = ReichMooreXs::new(&single_level(1.5e-3, 0.0, gf))
            .cross_sections(e)
            .unwrap();
        let exact = single_level_exact(e, 1.5e-3, 0.0, gf);
        assert_close(sigma.total, exact.total, "total");
        assert_close(sigma.elastic, exact.elastic, "elastic");
        assert_close(sigma.fission, exact.fission, "fission");
        // 捕獲の幅がないので、捕獲は差の丸めの誤差も残さず厳密に 0
        assert_eq!(sigma.capture, 0.0);
    }
}

#[test]
fn capture_is_exactly_zero_without_capture_width() {
    // 極でないエネルギーでも、捕獲の幅がなければ捕獲は差の丸めの誤差を残さない
    let xs = ReichMooreXs::new(&single_level(1.5e-3, 0.0, 1.0e-1));
    for e in [1.0e-3, 1.0, 6.7, 100.0] {
        assert_eq!(xs.cross_sections(e).unwrap().capture, 0.0, "at {e}");
    }
}

#[test]
fn two_poles_in_one_group_are_an_error() {
    // 同じ組に、捕獲幅のない共鳴が同じエネルギーに2つある場合（階数 2 の極限）には対応しない
    let mut data = single_level(1.5e-3, 0.0, 0.0);
    let mut second = data.l_values[0].resonances[0];
    second.gn = 2.0e-3;
    data.l_values[0].resonances.push(second);
    let error = ReichMooreXs::new(&data).cross_sections(6.674).unwrap_err();
    assert!(matches!(error, ReichMooreError::Poles { .. }), "{error}");
}

#[test]
fn negative_fission_amplitude_gives_same_single_level_result() {
    // 単一準位では核分裂幅の符号（振幅の符号）は結果に影響しない
    let xs = ReichMooreXs::new(&single_level(1.5e-3, 2.3e-2, -1.0e-1));
    let e = 6.7;
    assert_sigma(
        xs.cross_sections(e).unwrap(),
        single_level_exact(e, 1.5e-3, 2.3e-2, 1.0e-1),
        e,
    );
}

#[test]
fn awri_of_another_l_does_not_change_result() {
    // NAPS=0（チャネル半径を AWRI から求める）で、共鳴のない l の AWRI が違っても、
    // 共鳴のある l の断面積は変わらない
    let mut alone = single_level(1.5e-3, 2.3e-2, 1.0e-2);
    alone.naps = 0;
    let mut with_other = alone.clone();
    with_other.l_values.insert(
        0,
        RmLValue {
            awri: 1.0,
            apl: 0.0,
            l: OrbitalL::P,
            resonances: vec![],
        },
    );
    for e in [1.0, 6.674, 100.0] {
        assert_eq!(
            ReichMooreXs::new(&with_other).cross_sections(e).unwrap(),
            ReichMooreXs::new(&alone).cross_sections(e).unwrap(),
            "at {e}"
        );
    }
}

#[test]
fn p_wave_potential_scattering_sums_over_all_j() {
    // スピン 0 の標的の l=1 では J=1/2, 3/2 があり、gJ の和は 2l+1 = 3。
    // 共鳴から遠いエネルギーでは、どちらの J も硬球散乱 (4π/k²) gJ sin²φ1 だけになる
    let data = ReichMoore {
        spi: 0.0,
        ap: AP,
        naps: 1,
        l_values: vec![RmLValue {
            awri: AWRI,
            apl: 0.0,
            l: OrbitalL::P,
            resonances: vec![RmResonance {
                er: 1.0e4,
                aj: 1.5,
                gn: 1.0e-9,
                gg: 1.0e-9,
                gfa: 0.0,
                gfb: 0.0,
            }],
        }],
    };
    let e: f64 = 100.0;
    let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
    let rho = k * AP;
    let phi = rho - rho.atan();
    let expected = 4.0 * PI / (k * k) * 3.0 * phi.sin().powi(2);
    let sigma = ReichMooreXs::new(&data).cross_sections(e).unwrap();
    assert!((sigma.elastic - expected).abs() <= 1e-6 * expected);
    assert!((sigma.total - expected).abs() <= 1e-6 * expected);
}
