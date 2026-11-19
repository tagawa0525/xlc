use super::*;
use crate::endf::BreitWigner;
use crate::endf::{BwLValue, BwResonance, OrbitalL};
use crate::physics::{AMASSN, wave_number_factor};
use crate::resonance::reich_moore_tests::{AP, AWRI, single_level_exact};
use std::f64::consts::PI;

fn levels(spi: f64, naps: i64, l_values: Vec<BwLValue>) -> BreitWigner {
    BreitWigner {
        spi,
        ap: AP,
        ap_energy: None,
        naps,
        l_values,
    }
}

fn resonance(er: f64, aj: f64, (gn, gg, gf): (f64, f64, f64)) -> BwResonance {
    BwResonance { er, aj, gn, gg, gf }
}

fn assert_close(actual: f64, expected: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1e-6),
        "{what}: actual {actual}, expected {expected}"
    );
}

/// 剛体球の s 波と p 波の透過因子 P、シフト因子 S、位相シフト φ（ENDF-102 の表 D.1）
fn hard_sphere(l: OrbitalL, rho: f64) -> (f64, f64, f64) {
    match l {
        OrbitalL::S => (rho, 0.0, rho),
        OrbitalL::P => (
            rho.powi(3) / (1.0 + rho * rho),
            -1.0 / (1.0 + rho * rho),
            rho - rho.atan(),
        ),
        _ => panic!("only s and p waves"),
    }
}

/// ENDF-102 の付録 D の SLBW の式 (D.1)〜(D.7) を直接書いたもの（NAPS=0 のチャネル半径）。全断面積は部分反応の和
fn endf_formula(data: &BreitWigner, e: f64) -> Sigma {
    let mut sigma = Sigma::default();
    for lv in &data.l_values {
        let wave_number =
            |e: f64| wave_number_factor() * lv.awri / (lv.awri + 1.0) * e.abs().sqrt();
        let a = 0.123 * (AMASSN * lv.awri).cbrt() + 0.08;
        let k = wave_number(e);
        let pifac = PI / (k * k);
        let (p, s, _) = hard_sphere(lv.l, k * a);
        let (_, _, phi) = hard_sphere(lv.l, k * data.ap);
        let fl = f64::from(lv.l.number());
        sigma.elastic += 4.0 * pifac * (2.0 * fl + 1.0) * phi.sin().powi(2);
        for r in &lv.resonances {
            let (pr, sr, _) = hard_sphere(lv.l, wave_number(r.er) * a);
            let gn = p * r.gn / pr;
            let width = gn + r.gg + r.gf;
            let shifted = r.er + (sr - s) / (2.0 * pr) * r.gn;
            let g = (2.0 * r.aj + 1.0) / (2.0 * (2.0 * data.spi + 1.0));
            let d = e - shifted;
            let denom = d * d + width * width / 4.0;
            sigma.elastic += pifac
                * g
                * (gn * gn - 2.0 * gn * width * phi.sin().powi(2)
                    + 2.0 * d * gn * (2.0 * phi).sin())
                / denom;
            sigma.capture += pifac * g * gn * r.gg / denom;
            sigma.fission += pifac * g * gn * r.gf / denom;
        }
    }
    sigma.total = sigma.elastic + sigma.fission + sigma.capture;
    sigma
}

#[test]
fn single_s_wave_level_matches_exact_breit_wigner() {
    // 1 準位では準位の間の干渉が無いので、l=0 の SLBW は厳密な Breit-Wigner 形
    let widths = (1.5e-3, 2.3e-2, 1.0e-2);
    let data = levels(
        0.0,
        1,
        vec![BwLValue {
            awri: AWRI,
            l: OrbitalL::S,
            resonances: vec![resonance(6.674, 0.5, widths)],
        }],
    );
    let xs = SlbwXs::new(&data);
    for e in [1.0e-3, 1.0, 6.674, 6.7, 100.0] {
        let (a, b) = (
            xs.cross_sections(e).unwrap(),
            single_level_exact(e, widths.0, widths.1, widths.2),
        );
        assert_close(a.elastic, b.elastic, &format!("elastic at {e}"));
        assert_close(a.capture, b.capture, &format!("capture at {e}"));
        assert_close(a.fission, b.fission, &format!("fission at {e}"));
        assert_close(a.total, b.total, &format!("total at {e}"));
    }
}

/// 標的スピン 3/2、s 波の 3 準位（束縛準位を含む）と p 波の 1 準位（NAPS=0）
fn several_levels() -> BreitWigner {
    levels(
        1.5,
        0,
        vec![
            BwLValue {
                awri: AWRI,
                l: OrbitalL::S,
                resonances: vec![
                    resonance(-3.0, 1.0, (2.0e-3, 4.0e-2, 0.0)),
                    resonance(5.0, 2.0, (3.0e-3, 3.0e-2, 1.0e-2)),
                    resonance(12.0, 1.0, (8.0e-3, 3.5e-2, 2.0e-2)),
                ],
            },
            BwLValue {
                awri: AWRI,
                l: OrbitalL::P,
                resonances: vec![resonance(30.0, 2.5, (1.0e-4, 3.0e-2, 0.0))],
            },
        ],
    )
}

#[test]
fn several_levels_sum_the_endf_single_level_terms() {
    // 準位ごとの項の和（ENDF-102 の (D.2)、(D.4)、(D.6)）。p 波はシフト因子で共鳴エネルギーがずれる
    let data = several_levels();
    let xs = SlbwXs::new(&data);
    for e in [1.0e-3, 0.0253, 1.0, 4.9, 5.0, 8.0, 12.0, 29.9, 30.0, 100.0] {
        let (a, b) = (xs.cross_sections(e).unwrap(), endf_formula(&data, e));
        assert_close(a.elastic, b.elastic, &format!("elastic at {e}"));
        assert_close(a.capture, b.capture, &format!("capture at {e}"));
        assert_close(a.fission, b.fission, &format!("fission at {e}"));
        assert_close(a.total, b.total, &format!("total at {e}"));
    }
}

/// s 波の 2 準位。9 eV で Γn / |E - Er| が位相シフトの約 2 倍で、準位ごとの項がそれぞれ硬球散乱と同じ大きさで負になり、
/// 弾性散乱が負になる
pub(crate) fn cooperating_levels() -> BreitWigner {
    levels(
        0.0,
        0,
        vec![BwLValue {
            awri: AWRI,
            l: OrbitalL::S,
            resonances: vec![
                resonance(10.0, 0.5, (1.3e-2, 1.0e-3, 0.0)),
                resonance(10.5, 0.5, (2.0e-2, 1.0e-3, 0.0)),
            ],
        }],
    )
}

#[test]
fn cooperating_levels_give_the_negative_elastic_of_the_formula() {
    // 2 つの準位の、共鳴の下の負の干渉の裾が重なると、SLBW の弾性散乱は負になる（ENDF-102 の 2.4.9.1）。
    // 式の値をそのまま返し、負の扱いは MF3 のバックグラウンドと足したあとに決める
    let data = cooperating_levels();
    let xs = SlbwXs::new(&data);
    let elastic = xs.cross_sections(9.0).unwrap().elastic;
    assert!(elastic < 0.0, "elastic {elastic}");
    assert_close(elastic, endf_formula(&data, 9.0).elastic, "elastic at 9 eV");
}
