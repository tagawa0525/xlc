use super::*;
use crate::endf::{Tape, rec};

#[test]
fn reads_tabulated_nubar() {
    let text = [
        rec(
            ["9.223500+4", "2.330250+2", "0", "2", "0", "0"],
            9228,
            1,
            452,
        ),
        rec(["0.0", "0.0", "0", "0", "1", "2"], 9228, 1, 452),
        rec(["2", "2", "", "", "", ""], 9228, 1, 452),
        rec(["1.0-5", "2.4367", "2.0+7", "5.0", "", ""], 9228, 1, 452),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let nu = parse_nubar(&tape.section(9228, 1, 452).unwrap()).unwrap();
    assert_eq!(nu.x, vec![1.0e-5, 2.0e7]);
    assert_eq!(nu.y, vec![2.4367, 5.0]);
}

#[test]
fn polynomial_nubar_is_error() {
    let text = [
        rec(
            ["9.223500+4", "2.330250+2", "0", "1", "0", "0"],
            9228,
            1,
            452,
        ),
        rec(["0.0", "0.0", "0", "0", "2", "0"], 9228, 1, 452),
        rec(["2.4", "1.0-7", "", "", "", ""], 9228, 1, 452),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let err = parse_nubar(&tape.section(9228, 1, 452).unwrap()).unwrap_err();
    assert_eq!(err.line, 1);
}

/// MF1 MT452 が多項式（LNU=1）の係数 `coefficients` で、MF3 MT18 が 1e-5〜2e7 eV の表のテープ。
fn polynomial_nubar_tape(coefficients: &[&str], fission: bool) -> Tape {
    let mut fields = ["", "", "", "", "", ""];
    for (f, c) in fields.iter_mut().zip(coefficients) {
        *f = c;
    }
    let count = coefficients.len().to_string();
    let mut text = [
        rec(
            ["9.323900+4", "2.369990+2", "0", "1", "0", "0"],
            9352,
            1,
            452,
        ),
        rec(["0.0", "0.0", "0", "0", &count, "0"], 9352, 1, 452),
        rec(fields, 9352, 1, 452),
    ]
    .concat();
    if fission {
        text += &[
            rec(
                ["9.323900+4", "2.369990+2", "0", "0", "0", "0"],
                9352,
                3,
                18,
            ),
            rec(["1.9+8", "1.9+8", "0", "0", "1", "2"], 9352, 3, 18),
            rec(["2", "2", "", "", "", ""], 9352, 3, 18),
            rec(["1.0-5", "0.0", "2.0+7", "2.42", "", ""], 9352, 3, 18),
        ]
        .concat();
    }
    Tape::parse(&text).unwrap()
}

#[test]
fn polynomial_total_nubar_is_the_line_over_the_fission_range() {
    // ENDF/B-VII.0 の Np-239 の形。ν(E) = C1 + C2 E は、核分裂断面積のある範囲（ENDF-102 の 1.2.2）で有効
    let tape = polynomial_nubar_tape(&["2.74", "1.3-7"], true);
    let nu = total_nubar(&tape, 9352).unwrap().unwrap();
    assert_eq!(nu.x, vec![1.0e-5, 2.0e7]);
    for e in [1.0e-5, 0.0253, 1.0e6, 1.4e7, 2.0e7] {
        let expected = 2.74 + 1.3e-7 * e;
        let actual = nu.eval(e).unwrap();
        assert!(
            (actual - expected).abs() <= 1e-14 * expected,
            "{e}: {actual} vs {expected}"
        );
    }
}

#[test]
fn constant_polynomial_total_nubar_holds_over_the_fission_range() {
    let tape = polynomial_nubar_tape(&["2.5"], true);
    let nu = total_nubar(&tape, 9352).unwrap().unwrap();
    assert_eq!(nu.x, vec![1.0e-5, 2.0e7]);
    assert_eq!(nu.y, vec![2.5, 2.5]);
}

#[test]
fn curved_polynomial_total_nubar_is_error() {
    // 2 次以上の項は区分線形の表で厳密に表せない（対象の評価には無い）
    let tape = polynomial_nubar_tape(&["2.74", "1.3-7", "1.0-15"], true);
    assert!(total_nubar(&tape, 9352).is_err());
}

#[test]
fn polynomial_total_nubar_without_fission_is_error() {
    let tape = polynomial_nubar_tape(&["2.74", "1.3-7"], false);
    assert!(total_nubar(&tape, 9352).is_err());
}

#[test]
fn tabulated_total_nubar_is_the_table() {
    let text = [
        rec(
            ["9.223500+4", "2.330250+2", "0", "2", "0", "0"],
            9228,
            1,
            452,
        ),
        rec(["0.0", "0.0", "0", "0", "1", "2"], 9228, 1, 452),
        rec(["2", "2", "", "", "", ""], 9228, 1, 452),
        rec(["1.0-5", "2.4367", "2.0+7", "5.0", "", ""], 9228, 1, 452),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let nu = total_nubar(&tape, 9228).unwrap().unwrap();
    assert_eq!(nu.y, vec![2.4367, 5.0]);
    assert_eq!(total_nubar(&tape, 9999).unwrap(), None);
}

fn delayed_nubar_tape(ldg: &str) -> Tape {
    let text = [
        rec(
            ["9.223500+4", "2.330250+2", ldg, "2", "0", "0"],
            9228,
            1,
            455,
        ),
        // 族ごとの崩壊定数（2 族）
        rec(["0.0", "0.0", "0", "0", "2", "0"], 9228, 1, 455),
        rec(["1.3336-2", "3.2739-2", "", "", "", ""], 9228, 1, 455),
        rec(["0.0", "0.0", "0", "0", "1", "2"], 9228, 1, 455),
        rec(["2", "2", "", "", "", ""], 9228, 1, 455),
        rec(["1.0-5", "1.585-2", "2.0+7", "9.0-3", "", ""], 9228, 1, 455),
    ]
    .concat();
    Tape::parse(&text).unwrap()
}

#[test]
fn reads_the_decay_constants_and_the_delayed_nubar() {
    let tape = delayed_nubar_tape("0");
    let delayed = parse_delayed(&tape.section(9228, 1, 455).unwrap()).unwrap();
    assert_eq!(delayed.decay_constants, vec![1.3336e-2, 3.2739e-2]);
    assert_eq!(delayed.nubar.x, vec![1.0e-5, 2.0e7]);
    assert_eq!(delayed.nubar.y, vec![1.585e-2, 9.0e-3]);
}

#[test]
fn energy_dependent_decay_constants_are_error() {
    let tape = delayed_nubar_tape("1");
    let err = parse_delayed(&tape.section(9228, 1, 455).unwrap()).unwrap_err();
    assert_eq!(err.line, 1);
}
