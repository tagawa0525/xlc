use super::*;
use crate::endf::Tape;
use crate::endf::mf5_tests::two_point_mf5;
use crate::multigroup::group_integrals;

fn parse_mf5(text: &str) -> Mf5 {
    let tape = Tape::parse(text).unwrap();
    Mf5::parse(&tape.section(1, 5, 18).unwrap()).unwrap()
}

#[test]
fn spectrum_follows_incident_interpolation() {
    // 入射 1e-5 eV では高い群だけ、2e7 eV では低い群だけに出る。その間は E に線形
    let mf5 = parse_mf5(&two_point_mf5(
        "1",
        ("1.0+6", "2.0+7", "5.263158-8"),
        ("1.0-5", "1.0+6", "1.0-6"),
    ));
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    let weight = Weight::lwr();
    let constant = Tab1::linear(vec![1.0e-5, 2.0e7], vec![1.0, 1.0]);
    let chi = fission_spectrum(
        &groups,
        &weight,
        &constant,
        &[FissionNeutrons {
            nubar: &constant,
            spectrum: &mf5,
        }],
    );

    let t = |e: f64| (e - 1.0e-5) / (2.0e7 - 1.0e-5);
    let high: f64 = group_integrals(&groups, &weight, &[], |e| 1.0 - t(e))
        .iter()
        .sum();
    let all: f64 = group_integrals(&groups, &weight, &[], |_| 1.0).iter().sum();
    assert!((chi[0] - high / all).abs() < 1e-6, "{chi:?}");
    assert!((chi.iter().sum::<f64>() - 1.0).abs() < 1e-12);
}

#[test]
fn each_part_of_the_neutrons_enters_chi_by_its_nubar() {
    // 即発は ν = 2 ですべて高い群に、遅発は ν = 0.5 ですべて低い群に生まれる
    let high = parse_mf5(&two_point_mf5(
        "1",
        ("1.0+6", "2.0+7", "5.263158-8"),
        ("1.0+6", "2.0+7", "5.263158-8"),
    ));
    let low = parse_mf5(&two_point_mf5(
        "1",
        ("1.0-5", "1.0+6", "1.0-6"),
        ("1.0-5", "1.0+6", "1.0-6"),
    ));
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    let fission = Tab1::linear(vec![1.0e-5, 2.0e7], vec![1.0, 1.0]);
    let prompt = Tab1::linear(vec![1.0e-5, 2.0e7], vec![2.0, 2.0]);
    let delayed = Tab1::linear(vec![1.0e-5, 2.0e7], vec![0.5, 0.5]);
    let chi = fission_spectrum(
        &groups,
        &Weight::lwr(),
        &fission,
        &[
            FissionNeutrons {
                nubar: &prompt,
                spectrum: &high,
            },
            FissionNeutrons {
                nubar: &delayed,
                spectrum: &low,
            },
        ],
    );
    assert!((chi[0] - 0.8).abs() < 1e-9, "{chi:?}");
    assert!((chi[1] - 0.2).abs() < 1e-9, "{chi:?}");
}
