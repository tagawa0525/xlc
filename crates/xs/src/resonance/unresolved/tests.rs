use super::*;
use crate::endf::{OrbitalL, Tab1, UrrLValue};
use crate::special::erfc;

const AWRI: f64 = 236.0;
const AP: f64 = 0.948;

fn point(e: f64, gno: f64) -> UrrPoint {
    UrrPoint {
        e,
        d: 20.0,
        gx: 0.0,
        gno,
        gg: 0.023,
        gf: 0.0,
    }
}

fn s_wave(points: Vec<UrrPoint>, int: Interpolation) -> Unresolved {
    Unresolved {
        lssf: 0,
        spi: 0.0,
        ap: AP,
        ap_energy: None,
        naps: 1,
        l_values: vec![UrrLValue {
            awri: AWRI,
            l: OrbitalL::S,
            j_values: vec![UrrJValue {
                aj: 0.5,
                int,
                amux: WidthDof::Fixed,
                amun: WidthDof::One,
                amuf: WidthDof::Fixed,
                points,
            }],
        }],
    }
}

/// 平均 1 の χ²_1 での E[1/(x+c)] = c^(-1/2) e^(c/2) √(π/2) erfc(√(c/2))。
fn inverse_moment(c: f64) -> f64 {
    (PI / (2.0 * c)).sqrt() * (0.5 * c).exp() * erfc((0.5 * c).sqrt())
}

#[test]
fn energy_outside_the_radius_table_is_error() {
    let mut urr = s_wave(
        vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
        Interpolation::LinLin,
    );
    urr.ap_energy = Some(Tab1::linear(vec![5.0e3, 1.0e4], vec![AP, AP]));
    let xs = UnresolvedXs::new(&urr).unwrap();
    assert!(xs.cross_sections(8.0e3).is_ok());
    assert!(matches!(
        xs.cross_sections(1.5e4),
        Err(UrrError::RadiusTable { .. })
    ));
}

#[test]
fn s_wave_without_fission_matches_closed_form() {
    // GNO は 5e3 eV で 1e-3、2e4 eV で 4e-3 なので、1e4 eV では lin-lin で 2e-3
    let urr = s_wave(
        vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
        Interpolation::LinLin,
    );
    let e = 1.0e4;
    let sigma = UnresolvedXs::new(&urr).unwrap().cross_sections(e).unwrap();

    let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
    let (d, gg) = (20.0, 0.023);
    let gn = 2.0e-3 * e.sqrt();
    // 中性子幅 gn x、分母 gn x + gg の期待値（c = gg/gn）
    let c = gg / gn;
    let e_inv = inverse_moment(c);
    let e_x = (1.0 - c * e_inv) / gn;
    let e_x2 = (1.0 - c + c * c * e_inv) / gn;
    let scale = 2.0 * PI * PI / (k * k) / d;
    let phase = (k * AP).sin().powi(2);
    let capture = scale * gn * gg * e_x;
    let elastic = scale * gn * gn * e_x2 - scale * 2.0 * gn * phase + 4.0 * PI * phase / (k * k);
    for (actual, expected) in [(sigma.capture, capture), (sigma.elastic, elastic)] {
        assert!(
            (actual - expected).abs() <= 1e-8 * expected,
            "{actual} vs {expected}"
        );
    }
    assert_eq!(sigma.fission, 0.0);
    assert!((sigma.total - (sigma.elastic + sigma.capture)).abs() <= 1e-12 * sigma.total);
}

#[test]
fn grid_fills_wide_steps_like_njoy() {
    let urr = s_wave(
        vec![
            point(1.5e3, 1.0e-3),
            point(3.0e3, 1.0e-3),
            point(1.0e4, 1.0e-3),
            point(1.0e5, 1.0e-3),
        ],
        Interpolation::LinLin,
    );
    let expected = vec![
        1.5e3, 1.7e3, 2.0e3, 2.5e3, 3.0e3, 3.5e3, 4.0e3, 5.0e3, 6.0e3, 7.2e3, 8.5e3, 1.0e4, 1.25e4,
        1.5e4, 1.7e4, 2.0e4, 2.5e4, 3.0e4, 3.5e4, 4.0e4, 5.0e4, 6.0e4, 7.2e4, 8.5e4, 1.0e5,
    ];
    assert_eq!(unresolved_grid(&urr, 1.5e3, 1.0e5), expected);
}

#[test]
fn fluctuating_fission_width_converges() {
    // U-234（ENDF/B-VIII.0）の l=1, J=1/2 の系列。中性子幅と核分裂幅の入れ子の積分で、
    // 内側の積分のばらつきのために外側が収束しなかった
    let point = |e, gno, gf| UrrPoint {
        e,
        d: 10.6,
        gx: 0.0,
        gno,
        gg: 0.02,
        gf,
    };
    let urr = Unresolved {
        lssf: 0,
        spi: 0.0,
        ap: 0.893,
        ap_energy: None,
        naps: 0,
        l_values: vec![UrrLValue {
            awri: 232.029,
            l: OrbitalL::P,
            j_values: vec![UrrJValue {
                aj: 0.5,
                int: Interpolation::LinLin,
                amux: WidthDof::Fixed,
                amun: WidthDof::One,
                amuf: WidthDof::One,
                points: vec![
                    point(4.0e4, 1.5872e-3, 6.7737e-4),
                    point(6.0e4, 1.4566e-3, 6.0163e-4),
                ],
            }],
        }],
    };
    let sigma = UnresolvedXs::new(&urr)
        .unwrap()
        .cross_sections(5.0e4)
        .unwrap();
    assert!(sigma.fission > 0.0 && sigma.fission.is_finite());
}

#[test]
fn three_fluctuating_widths_converge() {
    // U-236（ENDF/B-VIII.0）の l=2, J=3/2 の系列。中性子・核分裂・競合の幅がゆらぐ
    let point = |e, gx| UrrPoint {
        e,
        d: 9.18,
        gx,
        gno: 1.836e-3,
        gg: 0.019,
        gf: 3.4e-4,
    };
    let urr = Unresolved {
        lssf: 0,
        spi: 0.0,
        ap: 0.93543,
        ap_energy: None,
        naps: 0,
        l_values: vec![UrrLValue {
            awri: 234.0178,
            l: OrbitalL::D,
            j_values: vec![UrrJValue {
                aj: 1.5,
                int: Interpolation::LogLog,
                amux: WidthDof::One,
                amun: WidthDof::One,
                amuf: WidthDof::One,
                points: vec![
                    point(4.5e4, 0.0),
                    point(5.0e4, 2.03e-5),
                    point(6.0e4, 3.33e-4),
                ],
            }],
        }],
    };
    let sigma = UnresolvedXs::new(&urr)
        .unwrap()
        .cross_sections(5.0e4)
        .unwrap();
    assert!(sigma.fission > 0.0 && sigma.fission.is_finite());
}

#[test]
fn shielded_cross_sections_at_infinite_dilution_match_averages() {
    let urr = s_wave(
        vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
        Interpolation::LinLin,
    );
    let xs = UnresolvedXs::new(&urr).unwrap();
    let e = 1.0e4;
    let average = xs.cross_sections(e).unwrap();
    let [dilute] = xs.shielded(e, 300.0, &[1.0e10], Sigma::default()).unwrap()[..] else {
        panic!("one σ0");
    };
    for (shielded, infinite) in [
        (dilute.capture, average.capture),
        (dilute.elastic, average.elastic),
        (dilute.total, average.total),
        (dilute.current_total, average.total),
    ] {
        assert!(
            (shielded - infinite).abs() <= 1e-6 * infinite,
            "{shielded} vs {infinite}"
        );
    }
}

#[test]
fn sequence_without_capture_width_is_skipped_as_in_averages() {
    // Γγ = 0 の系列は、平均断面積と同じく 0 として扱う（x → 0 で θ = Γ/Δ が 0 に近づくため）
    let mut zero_capture = point(5.0e3, 1.0e-3);
    zero_capture.gg = 0.0;
    let urr = s_wave(vec![zero_capture], Interpolation::LinLin);
    let xs = UnresolvedXs::new(&urr).unwrap();
    let average = xs.cross_sections(1.0e4).unwrap();
    let [dilute] = xs
        .shielded(1.0e4, 300.0, &[1.0e10], Sigma::default())
        .unwrap()[..]
    else {
        panic!("one σ0");
    };
    assert!((dilute.elastic - average.elastic).abs() <= 1e-12 * average.elastic);
    assert_eq!(dilute.capture, 0.0);
}

#[test]
fn capture_is_shielded_more_at_small_background_and_low_temperature() {
    let urr = s_wave(
        vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
        Interpolation::LinLin,
    );
    let xs = UnresolvedXs::new(&urr).unwrap();
    let sigma0 = [1.0e10, 1.0e3, 10.0, 1.0];
    let cold = xs
        .shielded(1.0e4, 300.0, &sigma0, Sigma::default())
        .unwrap();
    let hot = xs
        .shielded(1.0e4, 1200.0, &sigma0, Sigma::default())
        .unwrap();
    for w in cold.windows(2) {
        assert!(
            w[1].capture < w[0].capture,
            "{} then {}",
            w[0].capture,
            w[1].capture
        );
    }
    assert!(hot[2].capture > cold[2].capture);
}

#[test]
fn naps_two_needs_radius_table() {
    // NAPS=2 は、散乱半径の表（NRO=1）とは別に AP をチャネル半径に使う指定
    let mut urr = s_wave(vec![point(1.5e3, 1.0e-3)], Interpolation::LinLin);
    urr.naps = 2;
    assert!(UnresolvedXs::new(&urr).is_err());
}

#[test]
fn interpolation_must_agree_across_sequences() {
    let mut urr = s_wave(
        vec![point(1.5e3, 1.0e-3), point(1.0e5, 1.0e-3)],
        Interpolation::LogLog,
    );
    assert_eq!(interpolation(&urr).unwrap(), Interpolation::LogLog);
    let mut other = urr.l_values[0].j_values[0].clone();
    other.int = Interpolation::LinLin;
    urr.l_values[0].j_values.push(other);
    assert!(interpolation(&urr).is_err());
}
