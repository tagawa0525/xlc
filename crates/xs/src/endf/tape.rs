//! テープ全体の行の管理と、セクション内のレコードの読み込み。

use super::Error;
use super::number::{parse_int, parse_real};
use super::tab1::{Interpolation, Tab1};

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
        let lines = text
            .lines()
            .enumerate()
            .map(|(i, line)| Line::parse(i + 1, line))
            .collect::<Result<_, _>>()?;
        Ok(Self { lines })
    }

    /// 指定した MAT/MF/MT のセクションを返す。
    pub fn section(&self, mat: i32, mf: u32, mt: u32) -> Option<Section<'_>> {
        let id = (mat, mf, mt);
        let start = self.lines.iter().position(|l| l.id() == id)?;
        let len = self.lines[start..]
            .iter()
            .take_while(|l| l.id() == id)
            .count();
        Some(Section {
            lines: &self.lines[start..start + len],
        })
    }

    /// 指定した MAT/MF に含まれる MT を、テープ上の順に返す。
    pub fn mts(&self, mat: i32, mf: u32) -> Vec<u32> {
        let mut mts = Vec::new();
        for line in &self.lines {
            if line.mat == mat && line.mf == mf && line.mt != 0 && mts.last() != Some(&line.mt) {
                mts.push(line.mt);
            }
        }
        mts
    }
}

const FIELD_WIDTH: usize = 11;

impl Line {
    fn parse(number: usize, text: &str) -> Result<Self, Error> {
        let control = |range: std::ops::Range<usize>, name: &str| {
            let field = text
                .get(range)
                .ok_or_else(|| Error::new(number, "line is shorter than 75 columns"))?;
            parse_int(field)
                .ok_or_else(|| Error::new(number, format!("invalid {name} field `{field}`")))
        };
        let mat = control(66..70, "MAT")?;
        let mf = control(70..72, "MF")?;
        let mt = control(72..75, "MT")?;
        let invalid = |name: &str| Error::new(number, format!("{name} out of range"));
        Ok(Self {
            number,
            text: text.to_string(),
            mat: mat.try_into().map_err(|_| invalid("MAT"))?,
            mf: mf.try_into().map_err(|_| invalid("MF"))?,
            mt: mt.try_into().map_err(|_| invalid("MT"))?,
        })
    }

    fn id(&self) -> (i32, u32, u32) {
        (self.mat, self.mf, self.mt)
    }

    fn field(&self, index: usize) -> Result<&str, Error> {
        self.text
            .get(index * FIELD_WIDTH..(index + 1) * FIELD_WIDTH)
            .ok_or_else(|| Error::new(self.number, "data field is not ASCII"))
    }

    fn real(&self, index: usize) -> Result<f64, Error> {
        let field = self.field(index)?;
        parse_real(field)
            .ok_or_else(|| Error::new(self.number, format!("invalid real field `{field}`")))
    }

    fn int(&self, index: usize) -> Result<i64, Error> {
        let field = self.field(index)?;
        parse_int(field)
            .ok_or_else(|| Error::new(self.number, format!("invalid integer field `{field}`")))
    }
}

impl<'a> Section<'a> {
    pub fn reader(&self) -> Reader<'a> {
        Reader {
            lines: self.lines,
            pos: 0,
        }
    }
}

const FIELDS_PER_LINE: usize = 6;

impl<'a> Reader<'a> {
    pub fn cont(&mut self) -> Result<Cont, Error> {
        let line = self.next_line()?;
        Ok(Cont {
            c1: line.real(0)?,
            c2: line.real(1)?,
            l1: line.int(2)?,
            l2: line.int(3)?,
            n1: line.int(4)?,
            n2: line.int(5)?,
        })
    }

    pub fn list(&mut self) -> Result<List, Error> {
        let head = self.cont()?;
        let npl = self.count(head.n1)?;
        Ok(List {
            c1: head.c1,
            c2: head.c2,
            l1: head.l1,
            l2: head.l2,
            n2: head.n2,
            values: self.values(npl, Line::real)?,
        })
    }

    pub fn tab1(&mut self) -> Result<Tab1Record, Error> {
        let head = self.cont()?;
        let head_line = self.lines[self.pos - 1].number;
        let nr = self.count(head.n1)?;
        let np = self.count(head.n2)?;
        let ranges = self.values(2 * nr, |line, i| Ok((line.int(i)?, line.number)))?;
        let points = self.values(2 * np, Line::real)?;

        let mut breaks = Vec::with_capacity(nr);
        let mut laws = Vec::with_capacity(nr);
        for &[(nbt, nbt_line), (code, code_line)] in ranges.as_chunks::<2>().0 {
            let previous = breaks.last().copied().unwrap_or(0);
            let nbt = usize::try_from(nbt)
                .ok()
                .filter(|&b| b > previous && b <= np)
                .ok_or_else(|| Error::new(nbt_line, format!("invalid NBT {nbt}")))?;
            let law = Interpolation::from_code(code).ok_or_else(|| {
                Error::new(code_line, format!("unsupported interpolation law {code}"))
            })?;
            breaks.push(nbt);
            laws.push(law);
        }
        if breaks.last() != Some(&np) {
            return Err(Error::new(head_line, "last NBT must equal NP"));
        }

        let x: Vec<f64> = points.iter().step_by(2).copied().collect();
        let y: Vec<f64> = points.iter().skip(1).step_by(2).copied().collect();
        if x.windows(2)
            .any(|w| w[0].partial_cmp(&w[1]).is_none_or(|o| o.is_gt()))
        {
            return Err(Error::new(head_line, "x values are not in ascending order"));
        }
        Ok(Tab1Record {
            c1: head.c1,
            c2: head.c2,
            l1: head.l1,
            l2: head.l2,
            table: Tab1 { breaks, laws, x, y },
        })
    }

    fn next_line(&mut self) -> Result<&'a Line, Error> {
        let lines = self.lines;
        let line = lines.get(self.pos).ok_or_else(|| {
            let end = lines.last().map_or(0, |l| l.number + 1);
            Error::new(end, "unexpected end of section")
        })?;
        self.pos += 1;
        Ok(line)
    }

    fn count(&self, n: i64) -> Result<usize, Error> {
        usize::try_from(n).map_err(|_| {
            Error::new(
                self.lines[self.pos - 1].number,
                format!("negative count {n}"),
            )
        })
    }

    /// 1行6個ずつ並んだ n 個の値を読む。最終行の余りの欄は読まない。
    fn values<T>(
        &mut self,
        n: usize,
        read: impl Fn(&Line, usize) -> Result<T, Error>,
    ) -> Result<Vec<T>, Error> {
        let mut values = Vec::with_capacity(n);
        while values.len() < n {
            let line = self.next_line()?;
            let count = (n - values.len()).min(FIELDS_PER_LINE);
            for i in 0..count {
                values.push(read(line, i)?);
            }
        }
        Ok(values)
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
}
