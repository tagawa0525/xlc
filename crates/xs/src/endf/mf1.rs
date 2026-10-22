//! MF1 のうち、核分裂で放出される中性子数 ν（MT452、MT455、MT456）。

use super::Error;
use super::tab1::Tab1;
use super::tape::Section;

/// ν(E) の表を読む。M0 では表形式（LNU=2）だけを扱う。
pub fn parse_nubar(section: &Section<'_>) -> Result<Tab1, Error> {
    let mut r = section.reader();
    let head = r.cont()?;
    let lnu = head.l2;
    if lnu != 2 {
        return Err(Error::new(r.line(), format!("LNU={lnu} is not supported")));
    }
    Ok(r.tab1()?.table)
}

/// 遅発中性子の ν_d(E)（MT455）の表を読む。族ごとの崩壊定数（LDG=0 の LIST）は読み飛ばす。
/// 崩壊定数がエネルギーに依存する形（LDG=1）と、表形式でない ν（LNU=1）は扱わない。
pub fn parse_delayed_nubar(section: &Section<'_>) -> Result<Tab1, Error> {
    let mut r = section.reader();
    let head = r.cont()?;
    let (ldg, lnu) = (head.l1, head.l2);
    if ldg != 0 {
        return Err(Error::new(r.line(), format!("LDG={ldg} is not supported")));
    }
    if lnu != 2 {
        return Err(Error::new(r.line(), format!("LNU={lnu} is not supported")));
    }
    r.list()?;
    Ok(r.tab1()?.table)
}

#[cfg(test)]
mod tests;
