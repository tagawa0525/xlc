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
        " $Rev:: 1451     $  $Date:: 2018-01-17#$                             1 0  0\n".to_string(),
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
fn lists_sections_with_their_mts() {
    let tape = Tape::parse(&sample_tape()).unwrap();
    let sections = tape.sections(125, 3);
    let mts: Vec<u32> = sections.iter().map(|(mt, _)| *mt).collect();
    assert_eq!(mts, vec![1, 102]);
    let mut r = sections[1].1.reader();
    assert_eq!(r.cont().unwrap().c1, 1001.0);
    assert!(r.cont().is_err());
    assert!(tape.sections(125, 2).is_empty());
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
fn lists_materials() {
    let tape = Tape::parse(&sample_tape()).unwrap();
    assert_eq!(tape.materials(), vec![125]);
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
fn at_end_after_last_record() {
    let tape = Tape::parse(&sample_tape()).unwrap();
    let mut r = tape.section(125, 3, 1).unwrap().reader();
    assert!(!r.at_end());
    r.cont().unwrap();
    r.tab1().unwrap();
    assert!(r.at_end());
}

#[test]
fn reads_tab2() {
    let text = [
        rec(["0.0", "1.5", "3", "4", "2", "5"], 1, 5, 18),
        rec(["2", "1", "5", "22", "", ""], 1, 5, 18),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let t = tape.section(1, 5, 18).unwrap().reader().tab2().unwrap();
    assert_eq!((t.c2, t.l1, t.l2, t.nz), (1.5, 3, 4, 5));
    assert_eq!(t.breaks, vec![2, 5]);
    assert_eq!(t.codes, vec![1, 22]);
}

#[test]
fn tab2_last_break_must_equal_nz() {
    let text = [
        rec(["0.0", "0.0", "0", "0", "1", "5"], 1, 5, 18),
        rec(["4", "2", "", "", "", ""], 1, 5, 18),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    assert!(tape.section(1, 5, 18).unwrap().reader().tab2().is_err());
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

#[test]
fn huge_count_is_error_not_abort() {
    // 壊れたヘッダの件数で巨大な確保をせず、セクションの終わりとして報告する
    let tab1 = [
        rec(["0.0", "0.0", "0", "0", "1", "99999999999"], 1, 3, 1),
        rec(["2", "2", "", "", "", ""], 1, 3, 1),
    ]
    .concat();
    let tape = Tape::parse(&tab1).unwrap();
    let err = tape.section(1, 3, 1).unwrap().reader().tab1().unwrap_err();
    assert_eq!(err.line, 3);

    let list = rec(["0.0", "0.0", "0", "0", "99999999999", "0"], 1, 2, 151);
    let tape = Tape::parse(&list).unwrap();
    let err = tape
        .section(1, 2, 151)
        .unwrap()
        .reader()
        .list()
        .unwrap_err();
    assert_eq!(err.line, 2);
}

#[test]
fn non_ascii_text_counts_as_one_column_per_character() {
    // MF1/MT451 の説明文に非 ASCII 文字が入っていても、制御欄の桁はずれない
    let text = [
        format!("{:66}{:4}{:2}{:3}\n", " 10 µs, café", 125, 1, 451),
        rec(["1.0", "", "", "", "", ""], 125, 3, 1),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    assert_eq!(tape.mts(125, 1), vec![451]);
    assert!(tape.section(125, 3, 1).is_some());
}

#[test]
fn non_ascii_in_numeric_field_is_error() {
    let text = format!("{:>11}{:55}{:4}{:2}{:3}\n", "1.0µ", "", 1, 3, 1);
    let tape = Tape::parse(&text).unwrap();
    let err = tape.section(1, 3, 1).unwrap().reader().cont().unwrap_err();
    assert_eq!(err.line, 1);
}

fn log_table(law: i64, x: [&str; 2], y: [&str; 2]) -> String {
    let law = law.to_string();
    [
        rec(["0.0", "0.0", "0", "0", "1", "2"], 1, 3, 1),
        rec(["2", &law, "", "", "", ""], 1, 3, 1),
        rec([x[0], y[0], x[1], y[1], "", ""], 1, 3, 1),
    ]
    .concat()
}

#[test]
fn negative_y_with_log_y_law_is_error() {
    for law in [4, 5] {
        let tape = Tape::parse(&log_table(law, ["1.0", "2.0"], ["-1.0", "1.0"])).unwrap();
        let err = tape.section(1, 3, 1).unwrap().reader().tab1().unwrap_err();
        assert_eq!(err.line, 1, "INT={law}");
    }
}

#[test]
fn non_positive_x_with_log_x_law_is_error() {
    for law in [3, 5] {
        let tape = Tape::parse(&log_table(law, ["0.0", "2.0"], ["1.0", "1.0"])).unwrap();
        let err = tape.section(1, 3, 1).unwrap().reader().tab1().unwrap_err();
        assert_eq!(err.line, 1, "INT={law}");
    }
}

#[test]
fn zero_y_with_log_y_law_is_accepted() {
    let tape = Tape::parse(&log_table(5, ["1.0", "2.0"], ["0.0", "1.0"])).unwrap();
    assert!(tape.section(1, 3, 1).unwrap().reader().tab1().is_ok());
}
