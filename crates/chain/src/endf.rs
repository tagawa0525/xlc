//! 燃焼チェーンに使う ENDF の評価の読み込み（decay の MF8/MT457、中性子の MF3 の MT と MF8 の生成物の状態、nfy の MT454）。
//!
//! 行とレコードは xlc-xs の `endf::Tape` で読む（ADR-0166）。

use std::borrow::Cow;
use std::collections::HashSet;

use anyhow::{Context, Result, anyhow, bail};
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
    /// MF3 MT18 の QI (eV)。MT18 がなければ `None`
    pub fission_q: Option<f64>,
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

/// ZA の実数を (Z, A) にする（許容は MF8 の読み込みと同じ。`xlc_xs::endf::za`）。
fn za(value: f64, line: usize) -> Result<(u32, u32)> {
    let za = xlc_xs::endf::za(value, line)?;
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

/// 行の制御欄（MAT、MF、MT）。欄は `Tape::parse` が検査した。非 ASCII の文字は 1 桁に数える（`Tape::parse` と同じ）。
/// `Tape` は行を見せないので読み直す。重なりの検査を `Tape` に移せば消せる（#0304）
fn control(line: &str) -> (i64, i64, i64) {
    let line: Cow<str> = if line.is_ascii() {
        Cow::Borrowed(line)
    } else {
        Cow::Owned(
            line.chars()
                .map(|c| if c.is_ascii() { c } else { '?' })
                .collect(),
        )
    };
    let field = |range: std::ops::Range<usize>| {
        line.get(range)
            .map_or(0, |f| f.trim().parse::<i64>().unwrap_or(0))
    };
    (field(66..70), field(70..72), field(72..75))
}

/// テープの物質をテープ上の順に返す（ENDF/B-VII.0 の decay は全核種を 1 つのテープに並べる）。セクションは MAT、MF、MT で
/// 引くと最初の連続した行だけが返るので、同じ (MAT, MF, MT) が別の行を挟んで再び現れるテープ（同じ MAT の物質が 2 つ）は止める
fn materials(text: &str, tape: &Tape) -> Result<Vec<i32>> {
    let materials = tape.materials();
    if materials.is_empty() {
        bail!("the tape holds no material");
    }
    let mut seen = HashSet::new();
    let mut previous = None;
    for line in text.lines() {
        let id = control(line);
        // MT = 0 は SEND、FEND、MEND、TEND とテープの識別行で、セクションではない
        if previous == Some(id) || id.2 == 0 {
            previous = Some(id);
            continue;
        }
        previous = Some(id);
        if !seen.insert(id) {
            bail!(
                "MAT {} MF{} MT{} appears in two places on the tape; each material needs its own MAT",
                id.0,
                id.1,
                id.2
            );
        }
    }
    Ok(materials)
}

/// テープの物質ごとに `read` で評価を読む。エラーには物質の MAT を添える
fn each_material<T>(text: &str, read: impl Fn(&Tape, Material) -> Result<T>) -> Result<Vec<T>> {
    let tape = Tape::parse(text)?;
    materials(text, &tape)?
        .into_iter()
        .map(|mat| {
            read_material(&tape, mat)
                .and_then(|material| read(&tape, material))
                .with_context(|| format!("MAT {mat}"))
        })
        .collect()
}

/// MF1/MT451 から物質 `mat` の評価の対象を読む。
pub fn read_material(tape: &Tape, mat: i32) -> Result<Material> {
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

/// decay の評価のテープを読む（物質ごとに 1 評価、テープ上の順）。
pub fn read_decay(text: &str) -> Result<Vec<DecayEvaluation>> {
    each_material(text, decay_of)
}

fn decay_of(tape: &Tape, material: Material) -> Result<DecayEvaluation> {
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

/// 中性子の評価のテープを読む（物質ごとに 1 評価、テープ上の順）。`reaction_mts` の MF8 だけ生成物の状態を読む。
pub fn read_neutron(text: &str, reaction_mts: &[u32]) -> Result<Vec<NeutronEvaluation>> {
    each_material(text, |tape, material| {
        neutron_of(tape, material, reaction_mts)
    })
}

fn neutron_of(tape: &Tape, material: Material, reaction_mts: &[u32]) -> Result<NeutronEvaluation> {
    let mat = material.mat;
    let mut mts = tape.mts(mat, 3);
    mts.sort_unstable();
    // MF3 の 2 つ目のレコード（CONT）の C2 が QI
    let fission_q = match tape.section(mat, 3, 18) {
        Some(section) => {
            let mut r = section.reader();
            r.cont()?;
            Some(r.cont()?.c2)
        }
        None => None,
    };
    let mut final_states_of = Vec::new();
    let mut wanted = reaction_mts.to_vec();
    wanted.sort_unstable();
    for mt in wanted {
        if let Some(entries) = final_states(tape, mat, mt)? {
            final_states_of.push((mt, entries));
        }
    }
    Ok(NeutronEvaluation {
        material,
        mts,
        final_states: final_states_of,
        fission_q,
    })
}

/// nfy の評価のテープを読む（物質ごとに 1 評価、テープ上の順）。
pub fn read_yields(text: &str) -> Result<Vec<YieldEvaluation>> {
    each_material(text, yields_of)
}

fn yields_of(tape: &Tape, material: Material) -> Result<YieldEvaluation> {
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
