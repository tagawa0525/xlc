//! MF1 のうち、核分裂で放出される中性子数 ν（MT452, MT456）。

use super::Error;
use super::tab1::Tab1;
use super::tape::Section;

/// ν(E) の表を読む。M0 では表形式（LNU=2）だけを扱う。
pub fn parse_nubar(section: &Section<'_>) -> Result<Tab1, Error> {
    todo!()
}

#[cfg(test)]
mod tests {
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
}
