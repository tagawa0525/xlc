//! 物理定数。NJOY2016 (phys.f90) と同じ値を使い、処理結果を比較できるようにする。

/// 1 eV のエネルギー (erg)
pub const EV: f64 = 1.602_176_634e-12;
/// 光速 (cm/s)
pub const CLIGHT: f64 = 2.997_924_58e10;
/// 原子質量単位 (g)
pub const AMU: f64 = 931.494_102_42e6 * EV / (CLIGHT * CLIGHT);
/// ħ (erg s)
pub const HBAR: f64 = 6.582_119_569e-16 * EV;
/// 中性子質量 (amu)
pub const AMASSN: f64 = 1.008_664_915_95;
/// Boltzmann 定数 (eV/K)
pub const BK: f64 = 8.617_333_262e-5;

/// E (eV) の中性子の波数を 10^-12 cm 単位で表すための係数 sqrt(2 m_n eV) × 10^-12 / ħ。
/// 実際の波数は、これに AWRI/(AWRI+1) × sqrt(E) を掛けて求める。
pub fn wave_number_factor() -> f64 {
    (2.0 * AMASSN * AMU * EV).sqrt() * 1.0e-12 / HBAR
}
