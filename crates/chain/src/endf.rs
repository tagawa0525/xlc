//! 燃焼チェーンに使う ENDF の評価の読み込み（decay の MF8/MT457、中性子の MF3 の MT と MF8 の生成物の状態、nfy の MT454）。
//!
//! 行とレコードは xlc-xs の `endf::Tape` で読む（ADR-0166）。

use anyhow::{Result, anyhow, bail};
use xlc_xs::endf::{Tape, final_states};

/// 評価の対象（MF1/MT451）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    pub mat: i32,
    pub z: u32,
    pub a: u32,
    /// 準安定状態の番号（LISO）
    pub state: u32,
    /// 励起エネルギー (eV)（ELIS）
    pub excitation: f64,
}

/// 崩壊の 1 段。RTYP の 1 桁にあたる。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Step {
    BetaMinus,
    EcBetaPlus,
    IsomericTransition,
    Alpha,
    Neutron,
    SpontaneousFission,
    Proton,
}

/// 崩壊の型の 1 つ（MF8/MT457 の RTYP、RFS、BR）。
#[derive(Debug, Clone, PartialEq)]
pub struct DecayMode {
    pub steps: Vec<Step>,
    /// 娘の準安定状態の番号（RFS）
    pub final_state: u32,
    pub ratio: f64,
}

/// MF8/MT457 の崩壊のデータ。
#[derive(Debug, Clone, PartialEq)]
pub struct DecayData {
    /// NST = 1
    pub stable: bool,
    /// 半減期 (s)。安定なら 0
    pub half_life: f64,
    /// 平均の崩壊エネルギー (eV): 軽い粒子、電磁、重い粒子
    pub energy: [f64; 3],
    pub modes: Vec<DecayMode>,
}

/// decay のサブライブラリの 1 評価。
#[derive(Debug, Clone, PartialEq)]
pub struct DecayEvaluation {
    pub material: Material,
    pub decay: DecayData,
}

pub use xlc_xs::endf::FinalStateEntry;

/// 中性子のサブライブラリの 1 評価のうち、チェーンに使うもの。
#[derive(Debug, Clone, PartialEq)]
pub struct NeutronEvaluation {
    pub material: Material,
    /// MF3 の MT（昇順）
    pub mts: Vec<u32>,
    /// MF8 が生成物の状態を MF9 か MF10 で挙げる反応の MT と、その状態（MT の昇順）
    pub final_states: Vec<(u32, Vec<FinalStateEntry>)>,
}

/// 核分裂の 1 つの生成物の独立収率。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RawYield {
    pub za: u32,
    /// 準安定状態の番号（FPS）
    pub state: u32,
    pub value: f64,
}

/// nfy の 1 評価（MT454 の独立収率）。
#[derive(Debug, Clone, PartialEq)]
pub struct YieldEvaluation {
    pub material: Material,
    /// 入射エネルギー (eV) ごとの収率（テープの順）
    pub tables: Vec<(f64, Vec<RawYield>)>,
}

/// ZA（と ZAP）の実数を整数にするときの許容。ENDF/B-VIII.0 の Am-241 の MF8 は生成物を 95242.01 と書くので、それを
/// 受け、半端な値（書き誤り以上のずれ）は止める
const ZA_TOLERANCE: f64 = 0.05;

/// ZA の実数を (Z, A) にする。
fn za(value: f64, line: usize) -> Result<(u32, u32)> {
    let rounded = value.round();
    if !(rounded >= 1.0 && (value - rounded).abs() <= ZA_TOLERANCE && rounded < 1e6) {
        bail!("line {line}: ZA {value} is not an integer");
    }
    // 範囲は上で確かめた
    let za = rounded as u32;
    Ok((za / 1000, za % 1000))
}

/// 状態の番号などの、整数で書かれるべき実数。
fn whole(value: f64, what: &str, line: usize) -> Result<u32> {
    if !(value >= 0.0 && value.fract() == 0.0 && value < 1e3) {
        bail!("line {line}: {what} {value} is not a small non-negative integer");
    }
    Ok(value as u32)
}

/// 整数欄の非負の値。
fn count(value: i64, what: &str, line: usize) -> Result<u32> {
    u32::try_from(value).map_err(|_| anyhow!("line {line}: invalid {what} {value}"))
}

/// MF1/MT451 から評価の対象を読む。テープは物質を 1 つだけ持つ。
pub fn read_material(tape: &Tape) -> Result<Material> {
    let mat = match tape.materials()[..] {
        [mat] => mat,
        ref other => bail!("the tape must hold one material, found {other:?}"),
    };
    let section = tape
        .section(mat, 1, 451)
        .ok_or_else(|| anyhow!("MAT {mat} has no MF1/MT451"))?;
    let mut r = section.reader();
    let head = r.cont()?;
    let (z, a) = za(head.c1, r.line())?;
    let second = r.cont()?;
    Ok(Material {
        mat,
        z,
        a,
        state: count(second.l2, "LISO", r.line())?,
        excitation: second.c1,
    })
}

/// RTYP の 1 桁の崩壊の段。
fn step(digit: char, line: usize) -> Result<Step> {
    Ok(match digit {
        '1' => Step::BetaMinus,
        '2' => Step::EcBetaPlus,
        '3' => Step::IsomericTransition,
        '4' => Step::Alpha,
        '5' => Step::Neutron,
        '6' => Step::SpontaneousFission,
        '7' => Step::Proton,
        _ => bail!("line {line}: decay type digit {digit} is not a decay of the nuclide"),
    })
}

/// RTYP（各桁が順に起きる段。1.5 は β− のあと中性子の放出）を段の並びにする。10（不明）は止める。
fn steps(rtyp: f64, line: usize) -> Result<Vec<Step>> {
    if !(1.0..10.0).contains(&rtyp) {
        bail!("line {line}: unsupported decay type RTYP {rtyp}");
    }
    // 小数の桁は 6 桁まで（ENDF の 11 桁の欄に書ける桁）。末尾の 0 は段ではない
    let text = format!("{rtyp:.6}");
    text.chars()
        .filter(|&c| c != '.')
        .collect::<String>()
        .trim_end_matches('0')
        .chars()
        .map(|c| step(c, line))
        .collect()
}

/// decay の評価のテープを読む。
pub fn read_decay(text: &str) -> Result<DecayEvaluation> {
    let tape = Tape::parse(text)?;
    let material = read_material(&tape)?;
    let mat = material.mat;
    let section = tape
        .section(mat, 8, 457)
        .ok_or_else(|| anyhow!("MAT {mat} has no MF8/MT457"))?;
    let mut r = section.reader();
    let head = r.cont()?;
    let stable = match head.n1 {
        0 => false,
        1 => true,
        nst => bail!("line {}: invalid NST {nst}", r.line()),
    };
    let times = r.list()?;
    let line = r.line();
    let v = &times.values;
    if v.len() < 6 {
        bail!(
            "line {line}: {} average decay energies, expected at least 6",
            v.len()
        );
    }
    let energy = [v[0], v[2], v[4]];
    let spins = r.list()?;
    let line = r.line();
    let ndk =
        usize::try_from(spins.n2).map_err(|_| anyhow!("line {line}: invalid NDK {}", spins.n2))?;
    if spins.values.len() < 6 * ndk {
        bail!(
            "line {line}: {} values for {ndk} decay modes",
            spins.values.len()
        );
    }
    let modes = spins.values[..6 * ndk]
        .as_chunks::<6>()
        .0
        .iter()
        .map(|m| {
            Ok(DecayMode {
                steps: steps(m[0], line)?,
                final_state: whole(m[1], "RFS", line)?,
                ratio: m[4],
            })
        })
        .collect::<Result<_>>()?;
    Ok(DecayEvaluation {
        material,
        decay: DecayData {
            stable,
            half_life: times.c1,
            energy,
            modes,
        },
    })
}

/// 中性子の評価のテープを読む。`reaction_mts` の MF8 だけ生成物の状態を読む。
pub fn read_neutron(text: &str, reaction_mts: &[u32]) -> Result<NeutronEvaluation> {
    let tape = Tape::parse(text)?;
    let material = read_material(&tape)?;
    let mat = material.mat;
    let mut mts = tape.mts(mat, 3);
    mts.sort_unstable();
    let mut final_states_of = Vec::new();
    let mut wanted = reaction_mts.to_vec();
    wanted.sort_unstable();
    for mt in wanted {
        if let Some(entries) = final_states(&tape, mat, mt)? {
            final_states_of.push((mt, entries));
        }
    }
    Ok(NeutronEvaluation {
        material,
        mts,
        final_states: final_states_of,
    })
}

/// nfy の評価のテープを読む。
pub fn read_yields(text: &str) -> Result<YieldEvaluation> {
    let tape = Tape::parse(text)?;
    let material = read_material(&tape)?;
    let mat = material.mat;
    let section = tape
        .section(mat, 8, 454)
        .ok_or_else(|| anyhow!("MAT {mat} has no MF8/MT454"))?;
    let mut r = section.reader();
    let head = r.cont()?;
    let energies = count(head.l1, "LE + 1", r.line())?;
    let mut tables = Vec::new();
    for _ in 0..energies {
        let list = r.list()?;
        let line = r.line();
        let nfp = usize::try_from(list.n2)
            .map_err(|_| anyhow!("line {line}: invalid NFP {}", list.n2))?;
        if list.values.len() != 4 * nfp {
            bail!(
                "line {line}: {} values for {nfp} fission products",
                list.values.len()
            );
        }
        let yields = list
            .values
            .as_chunks::<4>()
            .0
            .iter()
            .map(|y| {
                let (z, a) = za(y[0], line)?;
                Ok(RawYield {
                    za: 1000 * z + a,
                    state: whole(y[1], "FPS", line)?,
                    value: y[2],
                })
            })
            .collect::<Result<_>>()?;
        tables.push((list.c1, yields));
    }
    if !r.at_end() {
        bail!("MAT {mat} MF8/MT454 has records after its {energies} energies");
    }
    Ok(YieldEvaluation { material, tables })
}

#[cfg(test)]
mod tests;
