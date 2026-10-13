//! MF1 のうち、核分裂で放出される中性子数 ν（MT452, MT456）。

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

#[cfg(test)]
mod tests;
