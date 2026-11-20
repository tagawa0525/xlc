use super::*;
use crate::endf::rec;
use crate::multigroup::ScatterError;

/// 物質 1 の MF の節（TAB1 の 2 点）。`q` は QM と QI、`x` と `y` は 2 点のエネルギー (eV) と値
fn section(mf: u32, mt: u32, q: &str, x: [&str; 2], y: [&str; 2]) -> String {
    section_of(1, mf, mt, q, x, y)
}

fn section_of(mat: i32, mf: u32, mt: u32, q: &str, x: [&str; 2], y: [&str; 2]) -> String {
    [
        rec(
            ["9.223800+4", "2.360060+2", "0", "0", "0", "0"],
            mat,
            mf,
            mt,
        ),
        rec([q, q, "0", "0", "1", "2"], mat, mf, mt),
        rec(["2", "2", "", "", "", ""], mat, mf, mt),
        rec([x[0], y[0], x[1], y[1], "", ""], mat, mf, mt),
        rec(["", "", "", "", "", ""], mat, mf, 0),
    ]
    .concat()
}

const FULL: [&str; 2] = ["1.0-5", "2.0+7"];

/// 全断面積、弾性散乱、捕獲だけを MF3 に持つ評価の節（xs がそのまま読める）
fn plain() -> Vec<String> {
    vec![
        section(3, 1, "0.0", FULL, ["11.0", "11.0"]),
        section(3, 2, "0.0", FULL, ["10.0", "10.0"]),
        section(3, 102, "0.0", FULL, ["1.0", "1.0"]),
    ]
}

fn tape(sections: &[String]) -> Tape {
    Tape::parse(&sections.concat()).unwrap()
}

#[test]
fn an_evaluation_that_xs_reads_passes() {
    assert!(check_evaluation(&tape(&plain())).is_ok());
}

#[test]
fn a_cross_section_below_its_threshold_stops_the_check() {
    // Q = −1 MeV の反応が 1e-5 eV で 0 でない
    let mut sections = plain();
    sections.insert(2, section(3, 16, "-1.0+6", FULL, ["0.2", "0.2"]));
    let error = check_evaluation(&tape(&sections)).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("below even the centre-of-mass threshold"),
        "{error}"
    );
}

#[test]
fn a_neutron_emission_without_a_distribution_stops_the_check() {
    // (n,2n) の断面積があり、MF4〜6 に出射の分布が無い
    let mut sections = plain();
    let threshold = ["1.004237+6", "2.0+7"];
    sections.insert(2, section(3, 16, "-1.0+6", threshold, ["0.2", "0.2"]));
    let error = check_evaluation(&tape(&sections)).unwrap_err();
    assert!(
        matches!(
            error,
            ProcessError::Scatter(ScatterError::NoDistribution(16))
        ),
        "{error}"
    );
}

#[test]
fn fission_without_nubar_stops_the_check() {
    let mut sections = plain();
    sections.insert(2, section(3, 18, "1.9+8", FULL, ["1.0", "1.0"]));
    let error = check_evaluation(&tape(&sections)).unwrap_err();
    assert!(matches!(error, ProcessError::FissionData), "{error}");
}

#[test]
fn a_final_state_without_its_table_stops_the_check() {
    // MF8 が捕獲の基底状態を MF9 で挙げ、MF9 が無い
    let mut sections = plain();
    sections.push(
        [
            rec(["9.223800+4", "2.360060+2", "0", "0", "1", "1"], 1, 8, 102),
            rec(["9.223900+4", "0.0", "9", "0", "0", "0"], 1, 8, 102),
            rec(["", "", "", "", "", ""], 1, 8, 0),
        ]
        .concat(),
    );
    let error = check_evaluation(&tape(&sections)).unwrap_err();
    assert!(matches!(error, ProcessError::Endf(_)), "{error}");
}

#[test]
fn a_tape_of_two_materials_stops_the_check() {
    let mut text = plain().concat();
    text += &section_of(2, 3, 1, "0.0", FULL, ["11.0", "11.0"]);
    let tape = Tape::parse(&text).unwrap();
    assert_eq!(tape.materials().len(), 2);
    let error = check_evaluation(&tape).unwrap_err();
    assert!(matches!(error, ProcessError::Materials(_)), "{error}");
}

#[test]
fn the_report_lists_each_file_and_the_reason_it_stopped() {
    let checked = [
        Checked {
            endf: "n-001_H_001.endf".into(),
            error: None,
        },
        Checked {
            endf: "n-023_V_000.endf".into(),
            error: Some("line 3: NIS=2 isotopes are not supported".into()),
        },
    ];
    assert_eq!(
        report_toml(&checked),
        "[[evaluation]]\n\
         endf = \"n-001_H_001.endf\"\n\
         \n\
         [[evaluation]]\n\
         endf = \"n-023_V_000.endf\"\n\
         error = \"line 3: NIS=2 isotopes are not supported\"\n"
    );
}
