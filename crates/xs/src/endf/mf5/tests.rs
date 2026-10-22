use super::*;
use crate::endf::{Tape, rec};

/// 入射エネルギー 1e-5 と 2e7 eV で、出射分布が [lo, hi] の一様分布の MF5。
pub(crate) fn two_point_mf5(
    lf: &str,
    first: (&str, &str, &str),
    second: (&str, &str, &str),
) -> String {
    two_point_mf5_for(18, lf, first, second)
}

/// `two_point_mf5` の分布を、MF5 の反応 `mt` として書く。
pub(crate) fn two_point_mf5_for(
    mt: u32,
    lf: &str,
    first: (&str, &str, &str),
    second: (&str, &str, &str),
) -> String {
    let r = |f: [&str; 6]| rec(f, 1, 5, mt);
    [
        r(["9.223500+4", "2.330248+2", "0", "0", "1", "0"]),
        r(["0.0", "0.0", "0", lf, "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["1.0-5", "1.0", "2.0+7", "1.0", "", ""]),
        r(["0.0", "0.0", "0", "0", "1", "2"]),
        // 入射エネルギーの間はふつうの線形補間
        r(["2", "2", "", "", "", ""]),
        r(["0.0", "1.0-5", "0", "0", "1", "2"]),
        r(["2", "1", "", "", "", ""]),
        r([first.0, first.2, first.1, "0.0", "", ""]),
        r(["0.0", "2.0+7", "0", "0", "1", "2"]),
        r(["2", "1", "", "", "", ""]),
        r([second.0, second.2, second.1, "0.0", "", ""]),
    ]
    .concat()
}

#[test]
fn reads_tabulated_distribution() {
    let text = two_point_mf5(
        "1",
        ("1.0+6", "2.0+7", "5.263158-8"),
        ("1.0-5", "1.0+6", "1.0-6"),
    );
    let tape = Tape::parse(&text).unwrap();
    let mf5 = Mf5::parse(&tape.section(1, 5, 18).unwrap()).unwrap();
    assert_eq!(mf5.subsections.len(), 1);
    let (p, EnergyLaw::Tabulated(d)) = &mf5.subsections[0] else {
        panic!("{:?}", mf5.subsections[0].1);
    };
    assert_eq!(p.eval(1.0), Some(1.0));
    assert_eq!(d.incident, vec![1.0e-5, 2.0e7]);
    assert_eq!(
        d.laws,
        vec![IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: false
        }]
    );
    assert_eq!(d.spectra[0].x, vec![1.0e6, 2.0e7]);
    assert_eq!(d.spectra[1].y, vec![1.0e-6, 0.0]);
}

/// p(E) = 1 で、θ(E) が 2 点の解析式の MF5。U は p(E) の TAB1 の C1。
fn analytic_mf5(lf: &str, restriction: &str) -> String {
    let r = |f: [&str; 6]| rec(f, 1, 5, 16);
    [
        r(["8.017000+3", "1.683600+1", "0", "0", "1", "0"]),
        r([restriction, "0.0", "0", lf, "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["4.0+6", "1.0", "2.0+7", "1.0", "", ""]),
        r(["0.0", "0.0", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["4.0+6", "2.0+5", "2.0+7", "6.0+5", "", ""]),
    ]
    .concat()
}

#[test]
fn reads_analytic_spectra() {
    for (lf, restriction) in [("7", "-3.0+7"), ("9", "1.0+6")] {
        let tape = Tape::parse(&analytic_mf5(lf, restriction)).unwrap();
        let mf5 = Mf5::parse(&tape.section(1, 5, 16).unwrap()).unwrap();
        let [(p, law)] = &mf5.subsections[..] else {
            panic!("{:?}", mf5.subsections);
        };
        assert_eq!(p.eval(1.0e7), Some(1.0));
        let (EnergyLaw::Maxwellian { restriction, theta }
        | EnergyLaw::Evaporation { restriction, theta }) = law
        else {
            panic!("{law:?}");
        };
        assert_eq!(theta.eval(1.2e7), 4.0e5);
        match lf {
            "7" => {
                assert!(matches!(law, EnergyLaw::Maxwellian { .. }) && *restriction == -3.0e7)
            }
            _ => assert!(matches!(law, EnergyLaw::Evaporation { .. }) && *restriction == 1.0e6),
        }
    }
}

#[test]
fn theta_must_be_positive() {
    let text = analytic_mf5("9", "1.0+6").replace("2.0+5", "0.0+0");
    let tape = Tape::parse(&text).unwrap();
    assert!(Mf5::parse(&tape.section(1, 5, 16).unwrap()).is_err());
}

#[test]
fn other_lf_is_error() {
    let text = two_point_mf5("11", ("1.0+6", "2.0+7", "1.0"), ("1.0-5", "1.0+6", "1.0"));
    let tape = Tape::parse(&text).unwrap();
    let err = Mf5::parse(&tape.section(1, 5, 18).unwrap()).unwrap_err();
    assert_eq!(err.line, 2);
}
