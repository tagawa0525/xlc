//! MF8 の反応の生成物の状態（LMF = 9 か 10 で、割合を MF9 の多重度か MF10 の生成断面積で与えるもの）。
//!
//! 燃焼鎖（xlc-chain。ADR-0164）と、鎖の反応の経路ごとの群の断面積（ADR-0172）が同じ読み込みを使う。

use super::tape::Tape;
use super::{Error, Tab1};

/// 中性子の反応の生成物の 1 状態（MF8 の MT = 反応の、LMF = 9 か 10 の項目）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FinalStateEntry {
    /// 生成物の ZA（1000 Z + A）
    pub za: u32,
    /// 励起エネルギー (eV)（ELFS）
    pub excitation: f64,
    /// 割合を与えるファイル（LMF）
    pub mf: u32,
    pub lfs: u32,
}

/// ZA（と ZAP）の実数を整数にするときの許容。ENDF/B-VIII.0 の Am-241 の MF8 は生成物を 95242.01 と書くので、それを
/// 受け、半端な値（書き誤り以上のずれ）は止める
const ZA_TOLERANCE: f64 = 0.05;

/// ZA の実数を 1000 Z + A にする（鎖の MF1 と収率の読み込みも使う）。
pub fn za(value: f64, line: usize) -> Result<u32, Error> {
    let rounded = value.round();
    if !(rounded >= 1.0 && (value - rounded).abs() <= ZA_TOLERANCE && rounded < 1e6) {
        return Err(Error::new(line, format!("ZA {value} is not an integer")));
    }
    // 範囲は上で確かめた
    Ok(rounded as u32)
}

/// 整数欄の非負の値。
fn count(value: i64, what: &str, line: usize) -> Result<u32, Error> {
    u32::try_from(value).map_err(|_| Error::new(line, format!("invalid {what} {value}")))
}

/// 物質 `mat` の MF8/MT`mt` が挙げる、MF9 か MF10 で割合を与える状態（テープの順）。MF8 がないか、どの状態も
/// MF9、MF10 でなければ `None`。一部の状態だけが MF9、MF10 ならエラー。
pub fn final_states(tape: &Tape, mat: i32, mt: u32) -> Result<Option<Vec<FinalStateEntry>>, Error> {
    let Some(section) = tape.section(mat, 8, mt) else {
        return Ok(None);
    };
    let mut r = section.reader();
    let head = r.cont()?;
    let ns = count(head.n1, "NS", r.line())?;
    let mut entries = Vec::new();
    for _ in 0..ns {
        let (c1, c2, l1, l2) = match head.n2 {
            0 => {
                let list = r.list()?;
                (list.c1, list.c2, list.l1, list.l2)
            }
            1 => {
                let cont = r.cont()?;
                (cont.c1, cont.c2, cont.l1, cont.l2)
            }
            no => return Err(Error::new(r.line(), format!("invalid NO {no}"))),
        };
        let line = r.line();
        let lmf = count(l1, "LMF", line)?;
        if lmf == 9 || lmf == 10 {
            entries.push(FinalStateEntry {
                za: za(c1, line)?,
                excitation: c2,
                mf: lmf,
                lfs: count(l2, "LFS", line)?,
            });
        }
    }
    if entries.is_empty() {
        return Ok(None);
    }
    if entries.len() != ns as usize {
        return Err(Error::new(
            section.reader().line(),
            format!("MAT {mat} MF8/MT{mt}: only some final states give MF9 or MF10"),
        ));
    }
    Ok(Some(entries))
}

/// MF9（多重度）か MF10（生成断面積）の、MT`mt` の状態 LFS = `lfs` の表。
pub fn final_state_table(tape: &Tape, mat: i32, mf: u32, mt: u32, lfs: u32) -> Result<Tab1, Error> {
    let section = tape.section(mat, mf, mt).ok_or_else(|| {
        Error::new(
            0,
            format!("MAT {mat} has no MF{mf}/MT{mt} for its MF8 final states"),
        )
    })?;
    let mut r = section.reader();
    let head = r.cont()?;
    let ns = count(head.n1, "NS", r.line())?;
    for _ in 0..ns {
        let record = r.tab1()?;
        if i64::from(lfs) == record.l2 {
            return Ok(record.table);
        }
    }
    Err(Error::new(
        r.line(),
        format!("MAT {mat} MF{mf}/MT{mt} has no final state LFS={lfs}"),
    ))
}
