//! ENDF-6 形式のテープの読み込み。
//!
//! 1行80桁の固定長書式で、66桁目までが11桁×6個のデータ欄、
//! 67〜75桁目が MAT/MF/MT の制御欄である（76〜80桁目の行番号は読まない）。

mod mf1;
mod mf2;
mod mf3;
mod mf4;
mod mf5;
mod mf6;
mod mf7;
mod mf8;
mod number;
mod tab1;
mod tape;

use std::fmt;

pub use mf1::parse_nubar;
pub use mf2::{
    BwLValue, BwResonance, Mf2, Mlbw, OrbitalL, ParticlePair, RMatrixLimited, RangeKind,
    ReichMoore, ResonanceRange, RmLValue, RmResonance, RmlChannel, RmlResonance, SpinGroup,
    Unresolved, UrrJValue, UrrLValue, UrrPoint, WidthDof, scattering_radius_at,
};
pub use mf3::Mf3;
pub use mf4::{AngleTables, Angular, Distribution, LegendreTable, Mf4};
pub use mf5::{EnergyLaw, IncidentLaw, Mf5, OutgoingDistributions, parse_outgoing};
pub use mf6::{Frame, Mf6, NeutronEmission};
pub use mf7::Mf7Incoherent;
pub use mf8::{FinalStateEntry, final_state_table, final_states};
pub use tab1::{ClampedTab1, Interpolation, Tab1};
pub use tape::{Cont, List, Reader, Section, Tab1Record, Tab2Record, Tape};

/// ENDF の読み込みで起きたエラー。`line` はテープ上の1始まりの行番号。
#[derive(Debug, Clone, PartialEq)]
pub struct Error {
    pub line: usize,
    pub message: String,
}

impl Error {
    pub(crate) fn new(line: usize, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for Error {}

#[cfg(test)]
pub(crate) use mf5::tests as mf5_tests;
#[cfg(test)]
pub(crate) use tape::tests::rec;
