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
        let mut reader = section.reader();
        let head = reader.cont()?;
        let tab = reader.tab1()?;
        Ok(Self {
            za: head.c1,
            awr: head.c2,
            qm: tab.c1,
            qi: tab.c2,
            lr: tab.l2,
            sigma: tab.table,
        })
    }
}

#[cfg(test)]
mod tests;
