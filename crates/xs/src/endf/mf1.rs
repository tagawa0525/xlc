//! MF1 のうち、核分裂で放出される中性子数 ν（MT452、MT455、MT456）。

use super::Error;
use super::mf3::Mf3;
use super::tab1::Tab1;
use super::tape::{Section, Tape};

/// ν(E) の表を読む。表形式（LNU=2）だけを扱う（MT456 の ν_p。多項式の MT452 は [`total_nubar`]）。
pub fn parse_nubar(section: &Section<'_>) -> Result<Tab1, Error> {
    let mut r = section.reader();
    let head = r.cont()?;
    let lnu = head.l2;
    if lnu != 2 {
        return Err(Error::new(r.line(), format!("LNU={lnu} is not supported")));
    }
    Ok(r.tab1()?.table)
}

/// 物質 `mat` の MF1 MT452 の ν(E)。MT452 がなければ `None`。
///
/// 多項式（LNU=1）の ν(E) = Σ C_n E^(n-1) は、核分裂断面積のある範囲で有効（ENDF-102 の 1.2.2）なので、MF3 MT18 の表の
/// 両端を結ぶ lin-lin の表にする。定数と 1 次の項だけを扱う（2 次以上の項は区分線形の表で厳密に表せない）。
pub fn total_nubar(tape: &Tape, mat: i32) -> Result<Option<Tab1>, Error> {
    let Some(section) = tape.section(mat, 1, 452) else {
        return Ok(None);
    };
    let mut r = section.reader();
    let lnu = r.cont()?.l2;
    if lnu != 1 {
        return parse_nubar(&section).map(Some);
    }
    let coefficients = r.list()?.values;
    let (c1, c2) = match coefficients[..] {
        [c1] => (c1, 0.0),
        [c1, c2] => (c1, c2),
        _ => {
            return Err(Error::new(
                r.line(),
                format!(
                    "LNU=1 with NC={} is not supported (only constant and linear terms)",
                    coefficients.len()
                ),
            ));
        }
    };
    let fission = tape.section(mat, 3, 18).ok_or_else(|| {
        Error::new(
            r.line(),
            "LNU=1 needs MF3 MT18 for the energy range of the polynomial",
        )
    })?;
    let x = Mf3::parse(&fission)?.sigma.x;
    let (Some(&low), Some(&high)) = (x.first(), x.last()) else {
        return Err(Error::new(r.line(), "MF3 MT18 has no points"));
    };
    let ends = vec![low, high];
    let y = ends.iter().map(|e| c1 + c2 * e).collect();
    Ok(Some(Tab1::linear(ends, y)))
}

/// 遅発中性子の族ごとの崩壊定数と ν_d(E)（MF1 MT455）。
#[derive(Debug, Clone, PartialEq)]
pub struct DelayedNubar {
    /// 族ごとの崩壊定数 (1/s)。LIST の順
    pub decay_constants: Vec<f64>,
    pub nubar: Tab1,
}

/// 遅発中性子の族ごとの崩壊定数（LDG=0 の LIST）と ν_d(E)（MT455）の表を読む。
/// 崩壊定数がエネルギーに依存する形（LDG=1）と、表形式でない ν（LNU=1）は扱わない。
pub fn parse_delayed(section: &Section<'_>) -> Result<DelayedNubar, Error> {
    let mut r = section.reader();
    let head = r.cont()?;
    let (ldg, lnu) = (head.l1, head.l2);
    if ldg != 0 {
        return Err(Error::new(r.line(), format!("LDG={ldg} is not supported")));
    }
    if lnu != 2 {
        return Err(Error::new(r.line(), format!("LNU={lnu} is not supported")));
    }
    let decay_constants = r.list()?.values;
    Ok(DelayedNubar {
        decay_constants,
        nubar: r.tab1()?.table,
    })
}

#[cfg(test)]
mod tests;
