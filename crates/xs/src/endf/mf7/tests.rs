use super::*;
use crate::endf::Tape;
use crate::endf::tape::tests::rec;

fn mf7(lln: &str, second_alpha: &str, s: [&str; 3]) -> Result<Mf7Incoherent, Error> {
    mf7_with(lln, second_alpha, s, "0.0", "1")
}

/// `list_beta` は β = 0 の2つ目の温度の LIST の β、`lt` は β = 0 の TAB1 の LT。
fn mf7_with(
    lln: &str,
    second_alpha: &str,
    s: [&str; 3],
    list_beta: &str,
    lt: &str,
) -> Result<Mf7Incoherent, Error> {
    let lines: Vec<[&str; 6]> = vec![
        ["1.001000+3", "9.991673-1", "0", "1", "0", "0"],
        // LIST: 0 0 LLN 0 NI NS、B = σfree M0, εmax, A0, Emax, 0, M0
        ["0.0", "0.0", lln, "0", "6", "0"],
        ["4.0", "10.0", "9.9917-1", "5.0", "0.0", "2.0"],
        // TAB2: β は 2 点、ln S が β について線形（INT=4）
        ["0.0", "0.0", "0", "0", "1", "2"],
        ["2", "4", "", "", "", ""],
        // β = 0: TAB1（T0 β LT 0 NR NP）と、2つ目の温度の LIST（T β LI 0 NP 0）
        ["2.936000+2", "0.0", lt, "0", "1", "3"],
        ["3", "4", "", "", "", ""],
        ["0.1", s[0], "1.0", s[1], "10.0", s[2]],
        ["6.000000+2", list_beta, "4", "0", "3", "0"],
        ["0.9", "0.6", "0.02", "", "", ""],
        // β = 2
        ["2.936000+2", "2.0", "1", "0", "1", "3"],
        ["3", "4", "", "", "", ""],
        ["0.1", "0.2", second_alpha, "0.3", "10.0", "0.005"],
        ["6.000000+2", "2.0", "4", "0", "3", "0"],
        ["0.25", "0.35", "0.006", "", "", ""],
        // 有効温度の TAB1
        ["0.0", "0.0", "0", "0", "1", "2"],
        ["2", "2", "", "", "", ""],
        ["2.936000+2", "1400.0", "6.000000+2", "1500.0", "", ""],
    ];
    let text: String = lines.iter().map(|f| rec(*f, 1, 7, 4)).collect();
    let tape = Tape::parse(&text).unwrap();
    Mf7Incoherent::parse(&tape.section(1, 7, 4).unwrap())
}

#[test]
fn reads_principal_scatterer() {
    let data = mf7("0", "1.0", ["1.0", "0.5", "0.01"]).unwrap();
    assert_eq!((data.lat, data.lasym), (1, 0));
    assert_eq!((data.sigma_free, data.awr, data.e_max), (2.0, 0.99917, 5.0));
    assert_eq!(data.alpha, vec![0.1, 1.0, 10.0]);
    assert_eq!(data.beta, vec![0.0, 2.0]);
    assert_eq!(
        (data.alpha_law, data.beta_law),
        (Interpolation::LogLin, Interpolation::LogLin)
    );
    assert_eq!(data.temperatures, vec![293.6, 600.0]);
    assert_eq!(data.s[0], vec![1.0, 0.5, 0.01, 0.2, 0.3, 0.005]);
    assert_eq!(data.s[1], vec![0.9, 0.6, 0.02, 0.25, 0.35, 0.006]);
    assert_eq!(data.effective_temperature.eval(600.0), Some(1500.0));
}

#[test]
fn rejects_logarithm_of_s() {
    // ENDF/B-VIII.0 の熱中性子散乱のデータはすべて LLN=0 なので、ln S での格納は扱わない
    assert!(mf7("1", "1.0", ["1.0", "0.5", "0.01"]).is_err());
}

#[test]
fn alpha_grid_must_be_shared_by_all_beta() {
    assert!(mf7("0", "1.5", ["1.0", "0.5", "0.01"]).is_err());
}

#[test]
fn temperature_records_must_have_the_same_beta() {
    assert!(mf7_with("0", "1.0", ["1.0", "0.5", "0.01"], "0.5", "1").is_err());
}

#[test]
fn number_of_temperatures_must_not_be_negative() {
    assert!(mf7_with("0", "1.0", ["1.0", "0.5", "0.01"], "0.0", "-1").is_err());
}
