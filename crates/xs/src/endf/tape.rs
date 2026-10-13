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

/// TAB2 レコード（続く NZ 個のレコードの補間則を与える）。
/// 補間則は ENDF の値のまま持つ（unit-base の 21〜25 などを使う側で解釈する）。
#[derive(Debug, Clone, PartialEq)]
pub struct Tab2Record {
    pub c1: f64,
    pub c2: f64,
    pub l1: i64,
    pub l2: i64,
    pub breaks: Vec<usize>,
    pub codes: Vec<i64>,
    pub nz: usize,
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

    /// テープに含まれる物質を、テープ上の順に返す。
    pub fn materials(&self) -> Vec<i32> {
        let mut materials = Vec::new();
        for line in &self.lines {
            // 先頭のテープ識別行（MF=0）は物質ではない
            if line.mat > 0 && line.mf > 0 && !materials.contains(&line.mat) {
                materials.push(line.mat);
            }
        }
        materials
    }

    /// 指定した MAT/MF のセクションを、MT とともにテープ上の順に返す。
    pub fn sections(&self, mat: i32, mf: u32) -> Vec<(u32, Section<'_>)> {
        self.lines
            .chunk_by(|a, b| a.id() == b.id())
            .filter(|lines| {
                let line = &lines[0];
                line.mat == mat && line.mf == mf && line.mt != 0
            })
            .map(|lines| (lines[0].mt, Section { lines }))
            .collect()
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
/// MT 欄の終わりの桁
const CONTROL_END: usize = 75;

impl Line {
    fn parse(number: usize, text: &str) -> Result<Self, Error> {
        // ENDF の桁は文字単位なので、非 ASCII 文字を1桁の `?` に置き換えて
        // バイト位置と桁位置をそろえる（数値欄にあれば数値の解析でエラーになる）
        let text = if text.is_ascii() {
            text.to_string()
        } else {
            text.chars()
                .map(|c| if c.is_ascii() { c } else { '?' })
                .collect()
        };
        if text.len() < CONTROL_END {
            return Err(Error::new(number, "line is shorter than 75 columns"));
        }
        let control = |range: std::ops::Range<usize>, name: &str| {
            let field = &text[range];
            parse_int(field)
                .ok_or_else(|| Error::new(number, format!("invalid {name} field `{field}`")))
        };
        let mat = control(66..70, "MAT")?;
        let mf = control(70..72, "MF")?;
        let mt = control(72..CONTROL_END, "MT")?;
        let invalid = |name: &str| Error::new(number, format!("{name} out of range"));
        Ok(Self {
            number,
            mat: mat.try_into().map_err(|_| invalid("MAT"))?,
            mf: mf.try_into().map_err(|_| invalid("MF"))?,
            mt: mt.try_into().map_err(|_| invalid("MT"))?,
            text,
        })
    }

    fn id(&self) -> (i32, u32, u32) {
        (self.mat, self.mf, self.mt)
    }

    /// 行は75桁以上の ASCII であることを `parse` で保証している。
    fn field(&self, index: usize) -> &str {
        &self.text[index * FIELD_WIDTH..(index + 1) * FIELD_WIDTH]
    }

    fn real(&self, index: usize) -> Result<f64, Error> {
        let field = self.field(index);
        parse_real(field)
            .ok_or_else(|| Error::new(self.number, format!("invalid real field `{field}`")))
    }

    fn int(&self, index: usize) -> Result<i64, Error> {
        let field = self.field(index);
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

    pub fn tab2(&mut self) -> Result<Tab2Record, Error> {
        let head = self.cont()?;
        let head_line = self.line();
        let nr = self.count(head.n1)?;
        let nz = self.count(head.n2)?;
        let ranges = self.values(self.pairs(nr)?, |line, i| line.int(i))?;
        let mut breaks = Vec::with_capacity(nr);
        let mut codes = Vec::with_capacity(nr);
        for &[nbt, code] in ranges.as_chunks::<2>().0 {
            let previous = breaks.last().copied().unwrap_or(0);
            let nbt = usize::try_from(nbt)
                .ok()
                .filter(|&b| b > previous && b <= nz)
                .ok_or_else(|| Error::new(head_line, format!("invalid NBT {nbt}")))?;
            breaks.push(nbt);
            codes.push(code);
        }
        if breaks.last() != Some(&nz) {
            return Err(Error::new(head_line, "last NBT must equal NZ"));
        }
        Ok(Tab2Record {
            c1: head.c1,
            c2: head.c2,
            l1: head.l1,
            l2: head.l2,
            breaks,
            codes,
            nz,
        })
    }

    pub fn tab1(&mut self) -> Result<Tab1Record, Error> {
        let head = self.cont()?;
        let head_line = self.lines[self.pos - 1].number;
        let nr = self.count(head.n1)?;
        let np = self.count(head.n2)?;
        let ranges = self.values(self.pairs(nr)?, |line, i| Ok((line.int(i)?, line.number)))?;
        let points = self.values(self.pairs(np)?, Line::real)?;

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

        let x = points.iter().step_by(2).copied().collect();
        let y = points.iter().skip(1).step_by(2).copied().collect();
        let table = Tab1 { breaks, laws, x, y };
        table
            .check()
            .map_err(|message| Error::new(head_line, message))?;
        Ok(Tab1Record {
            c1: head.c1,
            c2: head.c2,
            l1: head.l1,
            l2: head.l2,
            table,
        })
    }

    /// セクションのレコードをすべて読み終えたかどうか。
    pub fn at_end(&self) -> bool {
        self.pos >= self.lines.len()
    }

    /// 最後に読んだ行の行番号。まだ何も読んでいなければ、セクションの先頭行の番号。
    pub fn line(&self) -> usize {
        let index = self.pos.saturating_sub(1);
        self.lines.get(index).map_or(0, |l| l.number)
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

    /// n 組の対に含まれる値の個数。
    fn pairs(&self, n: usize) -> Result<usize, Error> {
        n.checked_mul(2).ok_or_else(|| {
            Error::new(
                self.lines[self.pos - 1].number,
                format!("count {n} is too large"),
            )
        })
    }

    /// 1行6個ずつ並んだ n 個の値を読む。最終行の余りの欄は読まない。
    fn values<T>(
        &mut self,
        n: usize,
        read: impl Fn(&Line, usize) -> Result<T, Error>,
    ) -> Result<Vec<T>, Error> {
        // 件数はテープの値なので、セクションに残っている欄の数を超えて確保しない
        let available = (self.lines.len() - self.pos) * FIELDS_PER_LINE;
        let mut values = Vec::with_capacity(n.min(available));
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
pub(crate) mod tests;
