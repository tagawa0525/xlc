use super::*;
use crate::endf::{BwLValue, BwResonance};
use crate::physics::wave_number_factor;
use crate::resonance::reich_moore_tests::{AP, AWRI, single_level_exact};
use std::f64::consts::PI;

fn single(l: OrbitalL, spi: f64, aj: f64, er: f64, (gn, gg, gf): (f64, f64, f64)) -> Mlbw {
    Mlbw {
        spi,
        ap: AP,
        naps: 1,
        l_values: vec![BwLValue {
            awri: AWRI,
            l,
            resonances: vec![BwResonance { er, aj, gn, gg, gf }],
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
fn single_s_wave_level_matches_exact_breit_wigner() {
    // l=0 ではレベルシフトがなく、単一準位の MLBW は厳密な Breit-Wigner 形
    let widths = (1.5e-3, 2.3e-2, 1.0e-2);
    let xs = MlbwXs::new(&single(OrbitalL::S, 0.0, 0.5, 6.674, widths));
    for e in [1.0e-3, 1.0, 6.674, 6.7, 100.0] {
        let (a, b) = (
            xs.cross_sections(e),
            single_level_exact(e, widths.0, widths.1, widths.2),
        );
        assert_close(a.elastic, b.elastic, &format!("elastic at {e}"));
        assert_close(a.capture, b.capture, &format!("capture at {e}"));
        assert_close(a.fission, b.fission, &format!("fission at {e}"));
        assert_close(a.total, b.total, &format!("total at {e}"));
    }
}

#[test]
fn missing_spin_states_add_hard_sphere_scattering() {
    // 標的スピン 1/2・l=1 では J=1 が2つのチャネルスピンから作れるが、gJ は1回しか数えないので、
    // 足りない分 (2l+1) - ΣgJ の硬球散乱を足す。共鳴から遠ければ全体で (4π/k²)(2l+1) sin²φ1
    let xs = MlbwXs::new(&single(OrbitalL::P, 0.5, 2.0, 1.0e4, (1.0e-9, 1.0e-9, 0.0)));
    let e: f64 = 100.0;
    let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
    let rho = k * AP;
    let phi = rho - rho.atan();
    let expected = 4.0 * PI / (k * k) * 3.0 * phi.sin().powi(2);
    let sigma = xs.cross_sections(e);
    assert!(
        (sigma.elastic - expected).abs() <= 1e-6 * expected,
        "{} vs {expected}",
        sigma.elastic
    );
}

#[test]
fn level_shift_vanishes_at_resonance_energy() {
    // E = Er ではシフトが 0 なので、捕獲のピークは π g Γn Γγ / (k² Γ²/4)
    let (er, gn, gg) = (2.0e3, 0.3, 0.05);
    let xs = MlbwXs::new(&single(OrbitalL::P, 0.0, 1.5, er, (gn, gg, 0.0)));
    let k = wave_number_factor() * AWRI / (AWRI + 1.0) * er.sqrt();
    let g = 2.0;
    let width = gn + gg;
    let expected = PI * g * gn * gg / (k * k * width * width / 4.0);
    assert_close(xs.cross_sections(er).capture, expected, "capture peak");
}
