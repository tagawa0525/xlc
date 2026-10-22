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
fn reads_delayed_nubar_after_the_decay_constants() {
    let tape = delayed_nubar_tape("0");
    let nu = parse_delayed_nubar(&tape.section(9228, 1, 455).unwrap()).unwrap();
    assert_eq!(nu.x, vec![1.0e-5, 2.0e7]);
    assert_eq!(nu.y, vec![1.585e-2, 9.0e-3]);
}

#[test]
fn energy_dependent_decay_constants_are_error() {
    let tape = delayed_nubar_tape("1");
    let err = parse_delayed_nubar(&tape.section(9228, 1, 455).unwrap()).unwrap_err();
    assert_eq!(err.line, 1);
}
