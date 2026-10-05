//! テープ全体の行の管理と、セクション内のレコードの読み込み。

use super::Error;
use super::tab1::Tab1;

/// 1本の ENDF テープ。
#[derive(Debug, Clone)]
pub struct Tape {
    lines: Vec<Line>,
}

#[derive(Debug, Clone)]
struct Line {
    /// 1始まりの行番号
    number: usize,
    text: String,
    mat: i32,
    mf: u32,
    mt: u32,
}

/// MAT/MF/MT が同じ連続した行（末尾の SEND レコードは含まない）。
#[derive(Debug, Clone, Copy)]
pub struct Section<'a> {
    lines: &'a [Line],
}

/// セクションの先頭からレコードを順に読む。
#[derive(Debug, Clone)]
pub struct Reader<'a> {
    lines: &'a [Line],
    pos: usize,
}

/// CONT レコード（HEAD レコードも同じ形）。
#[derive(Debug, Clone, PartialEq)]
pub struct Cont {
    pub c1: f64,
    pub c2: f64,
    pub l1: i64,
    pub l2: i64,
    pub n1: i64,
    pub n2: i64,
}

/// LIST レコード。NPL は `values.len()` で表す。
#[derive(Debug, Clone, PartialEq)]
pub struct List {
    pub c1: f64,
    pub c2: f64,
    pub l1: i64,
    pub l2: i64,
    pub n2: i64,
    pub values: Vec<f64>,
}

/// TAB1 レコード。NR と NP は `table` の長さで表す。
#[derive(Debug, Clone, PartialEq)]
pub struct Tab1Record {
    pub c1: f64,
    pub c2: f64,
    pub l1: i64,
    pub l2: i64,
    pub table: Tab1,
}

impl Tape {
    /// テープ全体の文字列を読み込む。各行の制御欄（MAT/MF/MT）を検査する。
    pub fn parse(text: &str) -> Result<Self, Error> {
        todo!()
    }

    /// 指定した MAT/MF/MT のセクションを返す。
    pub fn section(&self, mat: i32, mf: u32, mt: u32) -> Option<Section<'_>> {
        todo!()
    }

    /// 指定した MAT/MF に含まれる MT を、テープ上の順に返す。
    pub fn mts(&self, mat: i32, mf: u32) -> Vec<u32> {
        todo!()
    }
}

impl<'a> Section<'a> {
    pub fn reader(&self) -> Reader<'a> {
        todo!()
    }
}

impl Reader<'_> {
    pub fn cont(&mut self) -> Result<Cont, Error> {
        todo!()
    }

    pub fn list(&mut self) -> Result<List, Error> {
        todo!()
    }

    pub fn tab1(&mut self) -> Result<Tab1Record, Error> {
        todo!()
    }
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    use crate::endf::Interpolation;

    /// 6個のデータ欄と制御欄から1行を作る。
    pub(crate) fn rec(fields: [&str; 6], mat: i32, mf: u32, mt: u32) -> String {
        let data: String = fields.iter().map(|f| format!("{f:>11}")).collect();
        format!("{data}{mat:4}{mf:2}{mt:3}\n")
    }

    const BLANK: [&str; 6] = ["", "", "", "", "", ""];

    fn sample_tape() -> String {
        [
            " $Rev:: 1451     $  $Date:: 2018-01-17#$                             1 0  0\n"
                .to_string(),
            rec(["1.001000+3", "9.991673-1", "0", "0", "0", "0"], 125, 3, 1),
            rec(["0.0", "0.0", "0", "0", "1", "2"], 125, 3, 1),
            rec(["2", "2", "", "", "", ""], 125, 3, 1),
            rec(
                ["1.000000-5", "3.0", "2.000000+7", "1.0", "", ""],
                125,
                3,
                1,
            ),
            rec(BLANK, 125, 3, 0),
            rec(
                ["1.001000+3", "9.991673-1", "0", "0", "0", "0"],
                125,
                3,
                102,
            ),
            rec(BLANK, 125, 3, 0),
            rec(BLANK, 125, 0, 0),
            rec(BLANK, 0, 0, 0),
            rec(BLANK, -1, 0, 0),
        ]
        .concat()
    }

    #[test]
    fn finds_sections_and_mts() {
        let tape = Tape::parse(&sample_tape()).unwrap();
        assert_eq!(tape.mts(125, 3), vec![1, 102]);
        assert!(tape.section(125, 3, 2).is_none());
        let mut r = tape.section(125, 3, 102).unwrap().reader();
        assert_eq!(r.cont().unwrap().c1, 1001.0);
        // SEND はセクションに含めない
        assert!(r.cont().is_err());
    }

    #[test]
    fn reads_records_in_order() {
        let tape = Tape::parse(&sample_tape()).unwrap();
        let mut r = tape.section(125, 3, 1).unwrap().reader();
        let head = r.cont().unwrap();
        assert_eq!(
            head,
            Cont {
                c1: 1001.0,
                c2: 0.9991673,
                l1: 0,
                l2: 0,
                n1: 0,
                n2: 0
            }
        );
        let t = r.tab1().unwrap().table;
        assert_eq!(t.breaks, vec![2]);
        assert_eq!(t.laws, vec![Interpolation::LinLin]);
        assert_eq!(t.x, vec![1.0e-5, 2.0e7]);
        assert_eq!(t.y, vec![3.0, 1.0]);
    }

    #[test]
    fn reads_tab1_spanning_lines_with_two_ranges() {
        let text = [
            rec(["1.0", "2.0", "3", "4", "2", "4"], 1, 3, 1),
            rec(["2", "1", "4", "5", "", ""], 1, 3, 1),
            rec(["1.0", "1.0", "2.0", "2.0", "3.0", "3.0"], 1, 3, 1),
            rec(["4.0", "4.0", "", "", "", ""], 1, 3, 1),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let rec = tape.section(1, 3, 1).unwrap().reader().tab1().unwrap();
        assert_eq!((rec.c1, rec.c2, rec.l1, rec.l2), (1.0, 2.0, 3, 4));
        assert_eq!(rec.table.breaks, vec![2, 4]);
        assert_eq!(
            rec.table.laws,
            vec![Interpolation::Histogram, Interpolation::LogLog]
        );
        assert_eq!(rec.table.x, vec![1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn reads_list_spanning_lines() {
        let text = [
            rec(["1.0", "2.0", "3", "4", "7", "5"], 1, 2, 151),
            rec(["1", "2", "3", "4", "5", "6"], 1, 2, 151),
            rec(["7", "", "", "", "", ""], 1, 2, 151),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let list = tape.section(1, 2, 151).unwrap().reader().list().unwrap();
        assert_eq!(list.n2, 5);
        assert_eq!(list.values, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
    }

    #[test]
    fn accepts_lines_without_sequence_number() {
        // 76〜80桁目の行番号がないファイルがある
        let line = rec(["1.0", "", "", "", "", ""], 1, 3, 1);
        assert_eq!(line.trim_end_matches('\n').len(), 75);
        let tape = Tape::parse(&line).unwrap();
        assert!(tape.section(1, 3, 1).is_some());
    }

    #[test]
    fn truncated_record_reports_line() {
        let text = [
            rec(["1.0", "2.0", "0", "0", "1", "4"], 1, 3, 1),
            rec(["4", "2", "", "", "", ""], 1, 3, 1),
            rec(["1.0", "1.0", "2.0", "2.0", "3.0", "3.0"], 1, 3, 1),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let err = tape.section(1, 3, 1).unwrap().reader().tab1().unwrap_err();
        assert_eq!(err.line, 4);
    }

    #[test]
    fn invalid_number_reports_line() {
        let text = [
            rec(["1.0", "", "", "", "", ""], 1, 3, 1),
            rec(["x.y", "", "", "", "", ""], 1, 3, 1),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let mut r = tape.section(1, 3, 1).unwrap().reader();
        r.cont().unwrap();
        assert_eq!(r.cont().unwrap_err().line, 2);
    }

    #[test]
    fn unknown_interpolation_is_error() {
        let text = [
            rec(["0.0", "0.0", "0", "0", "1", "2"], 1, 3, 1),
            rec(["2", "6", "", "", "", ""], 1, 3, 1),
            rec(["1.0", "1.0", "2.0", "2.0", "", ""], 1, 3, 1),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let err = tape.section(1, 3, 1).unwrap().reader().tab1().unwrap_err();
        assert_eq!(err.line, 2);
    }

    #[test]
    fn invalid_control_field_is_error() {
        let text = [
            rec(["1.0", "", "", "", "", ""], 1, 3, 1),
            format!("{:66}  1 3 xx\n", ""),
        ]
        .concat();
        assert_eq!(Tape::parse(&text).unwrap_err().line, 2);
    }

    #[test]
    fn short_line_is_error() {
        assert_eq!(Tape::parse(" 1.0 125 3  1\n").unwrap_err().line, 1);
    }
}
