//! MF3（反応断面積）セクション。

use super::Error;
use super::tab1::Tab1;
use super::tape::Section;

/// MF3 の1反応の断面積。エネルギーは eV、断面積は barn。
#[derive(Debug, Clone, PartialEq)]
pub struct Mf3 {
    pub za: f64,
    /// 中性子質量を単位とする原子質量
    pub awr: f64,
    /// 質量差 Q 値 (eV)
    pub qm: f64,
    /// 反応 Q 値 (eV)
    pub qi: f64,
    /// 崩壊の仕方を表すフラグ（LR）
    pub lr: i64,
    pub sigma: Tab1,
}

impl Mf3 {
    pub fn parse(section: &Section<'_>) -> Result<Self, Error> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::Tape;
    use crate::endf::tape::tests::rec;

    #[test]
    fn parses_h1_capture() {
        // ENDF/B-VIII.0 H-1 MT102 の冒頭6点
        let text = [
            rec(
                ["1.001000+3", "9.991673-1", "0", "0", "0", "0"],
                125,
                3,
                102,
            ),
            rec(
                ["2.224648+6", "2.224648+6", "0", "0", "1", "6"],
                125,
                3,
                102,
            ),
            rec(["6", "5", "", "", "", ""], 125, 3, 102),
            rec(
                [
                    "1.000000-5",
                    "1.672869+1",
                    "2.000000-5",
                    "1.182897+1",
                    "5.000000-5",
                    "7.481299+0",
                ],
                125,
                3,
                102,
            ),
            rec(
                [
                    "1.000000-4",
                    "5.290077+0",
                    "2.000000-4",
                    "3.740649+0",
                    "5.000000-4",
                    "2.365794+0",
                ],
                125,
                3,
                102,
            ),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let mf3 = Mf3::parse(&tape.section(125, 3, 102).unwrap()).unwrap();
        assert_eq!(mf3.za, 1001.0);
        assert_eq!(mf3.awr, 0.9991673);
        assert_eq!(mf3.qi, 2.224648e6);
        assert_eq!(mf3.lr, 0);
        assert_eq!(mf3.sigma.x.len(), 6);
        assert_eq!(mf3.sigma.eval(2.0e-4), Some(3.740649));
    }
}
