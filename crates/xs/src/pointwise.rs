//! 0 K のポイントワイズの断面積（lin-lin で補間できる表）の作成。
//!
//! MF3 の各反応を許容誤差まで線形化し、分離共鳴領域では共鳴パラメータから計算した断面積を
//! 背景断面積に足す（NJOY の RECONR に相当）。

use std::collections::BTreeMap;

use crate::doppler::{broaden, broadening_limit};
use crate::endf::{Error, Interpolation, Mf2, Mf3, RangeKind, Tab1, Tape, Unresolved};
use crate::resonance::{
    CrossSections, MlbwXs, ReichMooreXs, Sigma, UnresolvedXs, unresolved_grid,
    unresolved_interpolation,
};

/// 線形化・共鳴計算の打ち切り判定: 区間の中点で、線形補間と真値の差が
/// `tol × |真値|` 以下になるまで区間を二分する。
pub const DEFAULT_TOLERANCE: f64 = 1.0e-3;

/// 共鳴の寄与を足す反応（全断面積・弾性散乱・核分裂・捕獲）。
pub const RESONANCE_MTS: [u32; 4] = [1, 2, 18, 102];

/// これより狭い区間は二分しない（x の相対幅）。
const MIN_RELATIVE_WIDTH: f64 = 1.0e-9;

/// 表を lin-lin 補間で `tol` の精度になるまで点を加えて線形化する。
/// 不連続点（x の重複）と、ヒストグラム補間の区間端の段差は、x の重複で表す。
pub fn linearize(table: &Tab1, tol: f64) -> Tab1 {
    let mut x = Vec::with_capacity(table.x.len());
    let mut y = Vec::with_capacity(table.y.len());
    for k in 0..table.x.len() {
        let (xk, yk) = (table.x[k], table.y[k]);
        if x.last() != Some(&xk) || y.last() != Some(&yk) {
            x.push(xk);
            y.push(yk);
        }
        let Some(&x1) = table.x.get(k + 1) else {
            break;
        };
        if x1 == xk {
            continue;
        }
        let law = table.law_of_interval(k);
        let f = |e: f64| [law.interpolate(e, (xk, yk), (x1, table.y[k + 1]))];
        // 右端は表の値をそのまま使う（補間し直すと丸めで値がずれ、偽の不連続点になる）。
        // ヒストグラムだけは区間内の値で終わり、次の点で段差になる
        let right = if law == Interpolation::Histogram {
            f(x1)
        } else {
            [table.y[k + 1]]
        };
        refine(&f, (xk, [yk]), (x1, right), tol, &mut |e, v| {
            x.push(e);
            y.push(v[0]);
        });
    }
    Tab1::linear(x, y)
}

/// 区間 (a, b) を、中点での線形補間の誤差が `tol` 以下になるまで二分し、
/// 加えた点と右端の点を左から順に `push` に渡す（左端は渡さない）。
fn refine<const N: usize>(
    f: &impl Fn(f64) -> [f64; N],
    a: (f64, [f64; N]),
    b: (f64, [f64; N]),
    tol: f64,
    push: &mut impl FnMut(f64, [f64; N]),
) {
    // 右側の区間から積み、左側から取り出す
    let mut stack = vec![b];
    let mut left = a;
    while let Some(right) = stack.pop() {
        let mid = 0.5 * (left.0 + right.0);
        let narrow = right.0 - left.0 <= MIN_RELATIVE_WIDTH * mid.abs();
        if !narrow {
            let exact = f(mid);
            let converged = (0..N).all(|i| {
                let linear = 0.5 * (left.1[i] + right.1[i]);
                (linear - exact[i]).abs() <= tol * exact[i].abs()
            });
            if !converged {
                stack.push(right);
                stack.push((mid, exact));
                continue;
            }
        }
        push(right.0, right.1);
        left = right;
    }
}

/// [el, eh] での共鳴の寄与を、反応ごとに lin-lin の表にする。
/// `nodes` は必ず格子に入れるエネルギー（共鳴のピークと半値幅の位置など）。
pub fn resonance_tables(
    xs: &impl CrossSections,
    nodes: &[f64],
    el: f64,
    eh: f64,
    tol: f64,
) -> [Tab1; 4] {
    let f = |e: f64| {
        let s = xs.cross_sections(e);
        RESONANCE_MTS.map(|mt| s.component(mt))
    };
    let mut seeds: Vec<f64> = nodes
        .iter()
        .copied()
        .filter(|&e| e > el && e < eh)
        .collect();
    seeds.push(el);
    seeds.push(eh);
    seeds.sort_by(f64::total_cmp);
    seeds.dedup();

    let mut x = vec![el];
    let mut values = vec![f(el)];
    for w in seeds.windows(2) {
        let left = (w[0], *values.last().expect("starts with el"));
        refine(&f, left, (w[1], f(w[1])), tol, &mut |e, v| {
            x.push(e);
            values.push(v);
        });
    }
    std::array::from_fn(|i| Tab1::linear(x.clone(), values.iter().map(|v| v[i]).collect()))
}

/// 2つの lin-lin の表の和。どちらの表も定義域の外では 0 とし、
/// 和の不連続点は x の重複で表す。
pub fn add(a: &Tab1, b: &Tab1) -> Tab1 {
    let mut breakpoints: Vec<f64> = a.x.iter().chain(&b.x).copied().collect();
    breakpoints.sort_by(f64::total_cmp);
    breakpoints.dedup();
    let last = breakpoints.len() - 1;
    let mut x = Vec::with_capacity(breakpoints.len());
    let mut y = Vec::with_capacity(breakpoints.len());
    for (i, &e) in breakpoints.iter().enumerate() {
        let left = left_limit(a, e) + left_limit(b, e);
        let right = right_limit(a, e) + right_limit(b, e);
        if i > 0 {
            x.push(e);
            y.push(left);
        }
        if i == 0 || (i < last && right != left) {
            x.push(e);
            y.push(right);
        }
    }
    Tab1::linear(x, y)
}

/// lin-lin の表の、e の左側からの極限値（定義域の外は 0）。
pub(crate) fn left_limit(t: &Tab1, e: f64) -> f64 {
    let (first, last) = (t.x[0], t.x[t.x.len() - 1]);
    if e <= first || e > last {
        return 0.0;
    }
    let k = t.x.partition_point(|&xi| xi < e);
    if t.x[k] == e {
        // 表の点の上では補間せず値をそのまま返す（丸めによる偽の不連続点を防ぐ）
        return t.y[k];
    }
    let (x0, x1) = (t.x[k - 1], t.x[k]);
    let (y0, y1) = (t.y[k - 1], t.y[k]);
    y0 + (y1 - y0) * (e - x0) / (x1 - x0)
}

/// lin-lin の表の、e の右側からの極限値（定義域の外は 0）。
pub(crate) fn right_limit(t: &Tab1, e: f64) -> f64 {
    let (first, last) = (t.x[0], t.x[t.x.len() - 1]);
    if e < first || e >= last {
        return 0.0;
    }
    t.eval(e).expect("inside the table range")
}

/// テープ中の1核種の、全反応の 0 K のポイントワイズの断面積。キーは MT。
///
/// NJOY の RECONR と同じく、全断面積などの冗長な反応は MF3 の値を使わず、
/// 部分反応の和として作り直す。全断面積と部分反応の和を一致させ、中性子の収支を保つためである。
pub fn reconstruct(tape: &Tape, mat: i32, tol: f64) -> Result<BTreeMap<u32, Tab1>, Error> {
    reconstruct_at(tape, mat, tol, 0.0)
}

/// テープ中の1核種の、温度 `temperature` (K) のポイントワイズの断面積。キーは MT。
///
/// 0 K の部分反応をドップラー拡がりで広げてから、冗長な反応を部分反応の和で作り直す。
/// `temperature` が 0 なら `reconstruct` と同じ。
pub fn reconstruct_at(
    tape: &Tape,
    mat: i32,
    tol: f64,
    temperature: f64,
) -> Result<BTreeMap<u32, Tab1>, Error> {
    Ok(ColdParts::read(tape, mat, tol)?.at(temperature))
}

/// 1核種の 0 K の部分反応（冗長な反応を和で作り直す前）。共鳴の再構成は重いので、
/// 複数の温度を扱うときは一度だけ作り、温度ごとに `at` で広げる。
#[derive(Debug, Clone)]
pub struct ColdParts {
    tables: BTreeMap<u32, Tab1>,
    mf2: Option<Mf2>,
    awr: Option<f64>,
    available: Vec<u32>,
    tol: f64,
}

impl ColdParts {
    /// MF3 と共鳴パラメータから、許容誤差 `tol` で 0 K の部分反応を作る。
    #[tracing::instrument(name = "reconstruct", level = "debug", skip_all)]
    pub fn read(tape: &Tape, mat: i32, tol: f64) -> Result<Self, Error> {
        let available = tape.mts(mat, 3);
        let mut tables = BTreeMap::new();
        let mut awr = None;
        for &mt in available.iter().filter(|&&mt| !is_dropped(mt, &available)) {
            let section = tape.section(mat, 3, mt).expect("listed by mts");
            let mf3 = Mf3::parse(&section)?;
            awr.get_or_insert(mf3.awr);
            tables.insert(mt, linearize(&mf3.sigma, tol));
        }
        let mut mf2 = None;
        if let Some(section) = tape.section(mat, 2, 151) {
            // MT18 を部分反応の和で作り直す場合は、共鳴の核分裂を第1チャンス核分裂に入れる
            let fission_mt = if is_rebuilt(18, &available) { 19 } else { 18 };
            let parsed = Mf2::parse(&section)?;
            add_resonances(&mut tables, &parsed, fission_mt, tol)
                .map_err(|message| Error::new(section.reader().line(), message))?;
            mf2 = Some(parsed);
        }
        Ok(Self {
            tables,
            mf2,
            awr,
            available,
            tol,
        })
    }

    /// 温度 `temperature` (K) のポイントワイズの断面積。キーは MT。
    pub fn at(&self, temperature: f64) -> BTreeMap<u32, Tab1> {
        let mut tables = self.tables.clone();
        if let (true, Some(awr)) = (temperature > 0.0, self.awr) {
            broaden_parts(&mut tables, self.mf2.as_ref(), awr, temperature, self.tol);
        }
        let total = add_all(
            tables
                .iter()
                .filter(|(mt, _)| is_part_of_total(**mt))
                .map(|(_, t)| t),
        );
        for sum in SUMS
            .iter()
            .filter(|sum| is_rebuilt(sum.mt, &self.available))
        {
            let parts = tables.iter().filter(|(mt, _)| (sum.is_part)(**mt));
            if let Some(t) = add_all(parts.map(|(_, t)| t)) {
                tables.insert(sum.mt, t);
            }
        }
        if let Some(total) = total {
            tables.insert(1, total);
        }
        tables
    }
}

/// 上限エネルギーより下から始まる反応を、まとめてドップラー拡がりで広げる。
fn broaden_parts(
    tables: &mut BTreeMap<u32, Tab1>,
    mf2: Option<&Mf2>,
    awr: f64,
    temperature: f64,
    tol: f64,
) {
    let candidates = || tables.iter().filter(|(mt, _)| is_broadened(**mt));
    let e_max = broadening_limit(mf2, candidates().map(|(_, t)| t));
    let mts: Vec<u32> = candidates()
        .filter(|(_, t)| t.x[0] < e_max)
        .map(|(mt, _)| *mt)
        .collect();
    if mts.is_empty() {
        return;
    }
    let cold: Vec<Tab1> = mts.iter().map(|mt| tables[mt].clone()).collect();
    for (mt, hot) in mts
        .into_iter()
        .zip(broaden(&cold, awr, temperature, e_max, tol))
    {
        tables.insert(mt, hot);
    }
}

/// ドップラー拡がりで広げる反応か（NJOY BROADR と同じ規則）。広げない反応は、拡げる上限の判定にも使わない。
///
/// 非弾性散乱の和（MT3, MT4）、MT46〜49、生成断面積など（MT201〜599）は広げない。
/// NJOY は MT18 を広げて MT19 を広げないが、ここでは MT19 がある場合に MT18 を部分反応の和で
/// 作り直すので、MT19 を広げる（MT18 を広げるのと同じ結果になる）。
pub fn is_broadened(mt: u32) -> bool {
    !matches!(mt, 3 | 4 | 46..=49 | 201..=599)
}

/// 部分反応の和で作り直す反応。
struct Sum {
    mt: u32,
    /// MT がこの和の部分反応かどうか
    is_part: fn(u32) -> bool,
    /// MF3 にこの MT があれば、和で作り直す
    trigger: fn(u32) -> bool,
}

const SUMS: [Sum; 8] = [
    Sum {
        mt: 4,
        is_part: |mt| (51..=91).contains(&mt),
        trigger: |mt| (51..=91).contains(&mt),
    },
    // NJOY と同じく、MT18 は MT19 がある場合だけ作り直す
    Sum {
        mt: 18,
        is_part: |mt| matches!(mt, 19..=21 | 38),
        trigger: |mt| mt == 19,
    },
    Sum {
        mt: 16,
        is_part: |mt| (875..=891).contains(&mt),
        trigger: |mt| (875..=891).contains(&mt),
    },
    Sum {
        mt: 103,
        is_part: |mt| (600..=649).contains(&mt),
        trigger: |mt| (600..=649).contains(&mt),
    },
    Sum {
        mt: 104,
        is_part: |mt| (650..=699).contains(&mt),
        trigger: |mt| (650..=699).contains(&mt),
    },
    Sum {
        mt: 105,
        is_part: |mt| (700..=749).contains(&mt),
        trigger: |mt| (700..=749).contains(&mt),
    },
    Sum {
        mt: 106,
        is_part: |mt| (750..=799).contains(&mt),
        trigger: |mt| (750..=799).contains(&mt),
    },
    Sum {
        mt: 107,
        is_part: |mt| (800..=849).contains(&mt),
        trigger: |mt| (800..=849).contains(&mt),
    },
];

/// `mt` が、MF3 の値を使わず部分反応の和で作り直す反応かどうか。
fn is_rebuilt(mt: u32, available: &[u32]) -> bool {
    SUMS.iter()
        .any(|sum| sum.mt == mt && available.iter().any(|&m| (sum.trigger)(m)))
}

/// MF3 から読まない反応: 和で作り直す反応と、断面積でない量（NJOY reconr lunion と同じ規則）。
fn is_dropped(mt: u32, available: &[u32]) -> bool {
    is_rebuilt(mt, available)
        || matches!(mt, 1 | 3 | 101 | 120 | 151 | 451..=599 | 851..=870)
        || ((251..=300).contains(&mt) && mt != 261)
        || mt > 891
}

/// 全断面積の和に含める反応。生成断面積（MT201〜）などは含めない。
/// 和で作り直す反応（MT4 など）は、部分反応があれば `is_dropped` で外れているので、
/// ここに残っているものは部分反応を持たない独立した反応である。
fn is_part_of_total(mt: u32) -> bool {
    !(mt == 10 || (46..=49).contains(&mt) || (201..=599).contains(&mt))
}

fn add_all<'a>(tables: impl Iterator<Item = &'a Tab1>) -> Option<Tab1> {
    tables.fold(None, |acc, t| match acc {
        None => Some(t.clone()),
        Some(acc) => Some(add(&acc, t)),
    })
}

/// 分離共鳴の寄与を弾性散乱・核分裂・捕獲に足す。
fn add_resonances(
    tables: &mut BTreeMap<u32, Tab1>,
    mf2: &Mf2,
    fission_mt: u32,
    tol: f64,
) -> Result<(), String> {
    for range in &mf2.ranges {
        let [_, elastic, fission, capture] = match &range.kind {
            RangeKind::ScatteringRadius { .. } => continue,
            RangeKind::Unresolved(urr) if urr.lssf == 1 => continue,
            RangeKind::UnresolvedLrf1 { lssf: 1 } => continue,
            RangeKind::Unresolved(urr) if urr.lssf == 0 => {
                // 平均で置き換えると、重なった部分の分離共鳴が消える（NJOY は分離共鳴の上端まで分離共鳴を使う）
                let resolved_top = mf2
                    .ranges
                    .iter()
                    .filter(|r| matches!(r.kind, RangeKind::ReichMoore(_) | RangeKind::Mlbw(_)))
                    .map(|r| r.eh)
                    .fold(f64::NEG_INFINITY, f64::max);
                if range.el < resolved_top {
                    return Err(format!(
                        "unresolved range from {} eV overlaps the resolved range up to {resolved_top} eV",
                        range.el
                    ));
                }
                replace_with_unresolved(tables, urr, range.el, range.eh, fission_mt, tol)?;
                continue;
            }
            RangeKind::Unresolved(Unresolved { lssf, .. }) | RangeKind::UnresolvedLrf1 { lssf } => {
                return Err(format!(
                    "unresolved range with LSSF={lssf} is not supported"
                ));
            }
            RangeKind::ReichMoore(rm) => {
                let peaks = rm
                    .l_values
                    .iter()
                    .flat_map(|lv| &lv.resonances)
                    .map(|r| (r.er, 0.5 * (r.gn + r.gg + r.gfa.abs() + r.gfb.abs())));
                let nodes = resonance_nodes(peaks);
                resonance_tables(&ReichMooreXs::new(rm), &nodes, range.el, range.eh, tol)
            }
            RangeKind::Mlbw(bw) => {
                let peaks = bw
                    .l_values
                    .iter()
                    .flat_map(|lv| &lv.resonances)
                    .map(|r| (r.er, 0.5 * (r.gn + r.gg + r.gf)));
                let nodes = resonance_nodes(peaks);
                resonance_tables(&MlbwXs::new(bw), &nodes, range.el, range.eh, tol)
            }
        };
        for (mt, contribution) in [(2, elastic), (fission_mt, fission), (102, capture)] {
            if contribution.y.iter().all(|&v| v == 0.0) {
                continue;
            }
            let combined = match tables.get(&mt) {
                Some(background) => add(background, &contribution),
                None => contribution,
            };
            tables.insert(mt, combined);
        }
    }
    Ok(())
}

/// LSSF=0 の非分離共鳴の範囲 [el, eh] で、弾性散乱・核分裂・捕獲を、平均断面積と MF3 の背景断面積の和に
/// 置き換える（NJOY RECONR と同じ）。和は格子点で求め、格子点の間は MF2 の補間則で補間してから線形化する。
fn replace_with_unresolved(
    tables: &mut BTreeMap<u32, Tab1>,
    urr: &Unresolved,
    el: f64,
    eh: f64,
    fission_mt: u32,
    tol: f64,
) -> Result<(), String> {
    // ENDF の読み込みのエラー（行番号つき）にまとめるので、ここでは文にする
    let xs = UnresolvedXs::new(urr).map_err(|e| e.to_string())?;
    let law = unresolved_interpolation(urr).map_err(|e| e.to_string())?;
    let grid = unresolved_grid(urr, el, eh);
    let averages: Vec<Sigma> = grid.iter().map(|&e| xs.cross_sections(e)).collect();
    for (mt, part) in [(2, 0), (fission_mt, 1), (102, 2)] {
        let background = tables.get(&mt);
        let values: Vec<f64> = grid
            .iter()
            .zip(&averages)
            .map(|(&e, s)| {
                // 上端では範囲の内側（左側）の値。MF3 は上端で範囲の外の値に段差で変わることがある
                let inside = background.map_or(0.0, |t| {
                    if e < eh {
                        right_limit(t, e)
                    } else {
                        left_limit(t, e)
                    }
                });
                [s.elastic, s.fission, s.capture][part] + inside
            })
            .collect();
        if background.is_none() && values.iter().all(|&v| v == 0.0) {
            continue;
        }
        let unresolved = linearize(
            &Tab1 {
                breaks: vec![grid.len()],
                laws: vec![law],
                x: grid.clone(),
                y: values,
            },
            tol,
        );
        let combined = match background {
            Some(t) => add(&outside(t, el, eh), &unresolved),
            None => unresolved,
        };
        tables.insert(mt, combined);
    }
    Ok(())
}

/// lin-lin の表の [el, eh] の中を 0 にした表（範囲の端は不連続点になる）。
fn outside(t: &Tab1, el: f64, eh: f64) -> Tab1 {
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for (&e, &v) in t.x.iter().zip(&t.y).filter(|(e, _)| **e < el) {
        x.push(e);
        y.push(v);
    }
    if el > t.x[0] {
        x.extend([el, el]);
        y.extend([left_limit(t, el), 0.0]);
    }
    if eh < t.x[t.x.len() - 1] {
        x.extend([eh, eh]);
        y.extend([0.0, right_limit(t, eh)]);
    }
    for (&e, &v) in t.x.iter().zip(&t.y).filter(|(e, _)| **e > eh) {
        x.push(e);
        y.push(v);
    }
    Tab1::linear(x, y)
}

/// 共鳴のピークと、その半値幅の位置。`peaks` は (共鳴エネルギー, 半値幅)。
fn resonance_nodes(peaks: impl Iterator<Item = (f64, f64)>) -> Vec<f64> {
    peaks
        .flat_map(|(er, half_width)| [er - half_width, er, er + half_width])
        .collect()
}

impl Sigma {
    fn component(&self, mt: u32) -> f64 {
        match mt {
            1 => self.total,
            2 => self.elastic,
            18 => self.fission,
            102 => self.capture,
            _ => unreachable!("not a resonance reaction: MT{mt}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::resonance::reich_moore_tests::{single_level, single_level_exact};

    fn linear(x: &[f64], y: &[f64]) -> Tab1 {
        Tab1 {
            breaks: vec![x.len()],
            laws: vec![Interpolation::LinLin],
            x: x.to_vec(),
            y: y.to_vec(),
        }
    }

    fn assert_linear(t: &Tab1) {
        assert_eq!(t.laws, vec![Interpolation::LinLin]);
        assert_eq!(t.breaks, vec![t.x.len()]);
    }

    #[test]
    fn linearizes_one_over_v_to_tolerance() {
        let source = Tab1 {
            breaks: vec![2],
            laws: vec![Interpolation::LogLog],
            x: vec![1.0e-5, 2.0e7],
            y: vec![1.0e3, 1.0e3 * (1.0e-5_f64 / 2.0e7).sqrt()],
        };
        let t = linearize(&source, 1.0e-3);
        assert_linear(&t);
        assert!(t.x.len() > 100);
        for w in t.x.windows(2) {
            let mid = 0.5 * (w[0] + w[1]);
            let exact = source.eval(mid).unwrap();
            let approx = t.eval(mid).unwrap();
            assert!((approx - exact).abs() <= 1.0e-3 * exact, "at {mid}");
        }
    }

    #[test]
    fn linearize_does_not_duplicate_continuous_points() {
        // 0.2 + (0.9 - 0.2) は 0.9 と 1ulp 違う。端点の値を補間し直すと偽の不連続点ができる
        let source = linear(&[1.0, 2.0, 3.0], &[0.2, 0.9, 0.9]);
        let t = linearize(&source, 1.0e-3);
        assert_eq!(t.x, vec![1.0, 2.0, 3.0]);
        assert_eq!(t.y, vec![0.2, 0.9, 0.9]);
    }

    #[test]
    fn linearize_keeps_discontinuities() {
        let source = linear(&[1.0, 2.0, 2.0, 3.0], &[1.0, 1.0, 5.0, 5.0]);
        let t = linearize(&source, 1.0e-3);
        assert_eq!(t.x, vec![1.0, 2.0, 2.0, 3.0]);
        assert_eq!(t.eval(2.0), Some(5.0));
    }

    #[test]
    fn resonance_tables_converge_at_midpoints() {
        let (gn, gg, gf) = (1.5e-3, 2.3e-2, 1.0e-2);
        let xs = ReichMooreXs::new(&single_level(gn, gg, gf));
        let tables = resonance_tables(&xs, &[6.674], 1.0e-5, 1.0e3, 1.0e-3);
        for (t, mt) in tables.iter().zip(RESONANCE_MTS) {
            assert_linear(t);
            assert_eq!((t.x[0], *t.x.last().unwrap()), (1.0e-5, 1.0e3));
            for w in t.x.windows(2) {
                let mid = 0.5 * (w[0] + w[1]);
                let exact = single_level_exact(mid, gn, gg, gf).component(mt);
                let approx = t.eval(mid).unwrap();
                assert!(
                    (approx - exact).abs() <= 1.0e-3 * exact.abs() + 1.0e-12,
                    "MT{mt} at {mid}: {approx} vs {exact}"
                );
            }
        }
        // 共鳴のピークが格子に入っている
        assert!(tables[0].x.contains(&6.674));
    }

    #[test]
    fn add_inserts_discontinuities_at_range_ends() {
        let a = linear(&[1.0, 10.0], &[1.0, 1.0]);
        let b = linear(&[2.0, 5.0], &[3.0, 3.0]);
        let sum = add(&a, &b);
        assert_linear(&sum);
        assert_eq!(sum.x, vec![1.0, 2.0, 2.0, 5.0, 5.0, 10.0]);
        assert_eq!(sum.y, vec![1.0, 1.0, 4.0, 4.0, 1.0, 1.0]);
    }

    #[test]
    fn add_keeps_existing_discontinuity_of_a() {
        let a = linear(&[1.0, 4.0, 4.0, 10.0], &[1.0, 1.0, 2.0, 2.0]);
        let b = linear(&[1.0, 10.0], &[1.0, 1.0]);
        let sum = add(&a, &b);
        assert_eq!(sum.x, vec![1.0, 4.0, 4.0, 10.0]);
        assert_eq!(sum.y, vec![2.0, 2.0, 3.0, 3.0]);
    }

    /// MF3 の1反応（2点の lin-lin の表）
    fn mf3_section(mt: u32, x: [&str; 2], y: [&str; 2]) -> String {
        use crate::endf::rec;
        [
            rec(["1.001000+3", "1.0", "0", "0", "0", "0"], 1, 3, mt),
            rec(["0.0", "0.0", "0", "0", "1", "2"], 1, 3, mt),
            rec(["2", "2", "", "", "", ""], 1, 3, mt),
            rec([x[0], y[0], x[1], y[1], "", ""], 1, 3, mt),
            rec(["", "", "", "", "", ""], 1, 3, 0),
        ]
        .concat()
    }

    #[test]
    fn redundant_reactions_are_sums_of_their_parts() {
        let full = ["1.0-5", "2.0+7"];
        let text = [
            // MT1, MT4, MT18 にはわざと部分反応の和と違う値を入れておく
            mf3_section(1, full, ["100.0", "100.0"]),
            mf3_section(2, full, ["10.0", "10.0"]),
            mf3_section(4, ["1.0+6", "2.0+7"], ["50.0", "50.0"]),
            mf3_section(18, full, ["70.0", "70.0"]),
            mf3_section(19, full, ["0.25", "0.25"]),
            mf3_section(20, ["5.0+6", "2.0+7"], ["0.125", "0.125"]),
            mf3_section(51, ["1.0+6", "2.0+7"], ["1.0", "1.0"]),
            mf3_section(52, ["2.0+6", "2.0+7"], ["2.0", "2.0"]),
            mf3_section(102, full, ["0.5", "0.5"]),
            // 生成断面積は全断面積に含めない
            mf3_section(203, ["1.0+6", "2.0+7"], ["3.0", "3.0"]),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let tables = reconstruct(&tape, 1, 1.0e-3).unwrap();
        let at = |mt: u32, e: f64| tables[&mt].eval(e).unwrap();

        assert_eq!(at(1, 1.0), 10.0 + 0.25 + 0.5);
        assert_eq!(at(1, 1.5e6), 10.0 + 0.25 + 1.0 + 0.5);
        assert_eq!(at(1, 3.0e6), 10.0 + 0.25 + 1.0 + 2.0 + 0.5);
        assert_eq!(at(1, 6.0e6), 10.0 + 0.25 + 0.125 + 1.0 + 2.0 + 0.5);
        assert_eq!(at(4, 1.5e6), 1.0);
        assert_eq!(at(4, 3.0e6), 3.0);
        assert_eq!(at(18, 1.0), 0.25);
        assert_eq!(at(18, 6.0e6), 0.375);
        assert_eq!(at(203, 3.0e6), 3.0);
    }

    #[test]
    fn reconstruct_at_temperature_broadens_parts_and_sums_total() {
        let full = ["1.0-5", "2.0+7"];
        let text = [
            mf3_section(1, full, ["100.0", "100.0"]),
            mf3_section(2, full, ["10.0", "10.0"]),
            mf3_section(102, full, ["0.5", "0.5"]),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let cold = reconstruct(&tape, 1, 1.0e-3).unwrap();
        assert_eq!(reconstruct_at(&tape, 1, 1.0e-3, 0.0).unwrap(), cold);

        let hot = reconstruct_at(&tape, 1, 1.0e-3, 600.0).unwrap();
        // 一定の断面積は、拡げると低エネルギーで大きくなる
        let e = 1.0e-3;
        assert!(hot[&2].eval(e).unwrap() > 1.5 * cold[&2].eval(e).unwrap());
        for &e in &hot[&1].x {
            let parts = hot[&2].eval(e).unwrap() + hot[&102].eval(e).unwrap();
            assert!(
                (hot[&1].eval(e).unwrap() - parts).abs() <= 1.0e-12 * parts,
                "at {e}"
            );
        }
    }

    #[test]
    fn broadening_limit_ignores_reactions_njoy_does_not_broaden() {
        // 共鳴パラメータがない核種で、NJOY BROADR が見ない MT4（部分反応なし）と生成断面積は
        // 上限を下げない。上限は格子の点として出力に現れる
        let full = ["1.0-5", "2.0+7"];
        let text = [
            mf3_section(2, full, ["10.0", "10.0"]),
            mf3_section(4, ["1.0+6", "2.0+7"], ["1.0", "1.0"]),
            mf3_section(102, full, ["0.5", "0.5"]),
            mf3_section(203, ["1.0+5", "2.0+7"], ["3.0", "3.0"]),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let hot = reconstruct_at(&tape, 1, 1.0e-3, 600.0).unwrap();
        assert!(hot[&2].x.contains(&crate::doppler::MAX_LIMIT));
    }

    /// 非分離共鳴だけの MF2（LSSF=0、s 波の1系列、パラメータは一定）
    fn unresolved_mf2() -> String {
        use crate::endf::rec;
        [
            rec(["1.001000+3", "1.0", "0", "0", "1", "0"], 1, 2, 151),
            rec(["1.001000+3", "1.0", "0", "0", "1", "0"], 1, 2, 151),
            rec(["1.0+3", "1.0+5", "2", "2", "0", "1"], 1, 2, 151),
            rec(["0.0", "0.948", "0", "0", "1", "0"], 1, 2, 151),
            rec(["2.360060+2", "0.0", "0", "0", "1", "0"], 1, 2, 151),
            rec(["0.5", "0.0", "2", "0", "18", "2"], 1, 2, 151),
            rec(["0.0", "0.0", "0.0", "1.0", "0.0", "0.0"], 1, 2, 151),
            rec(["1.0+3", "20.0", "0.0", "1.0-3", "2.3-2", "0.0"], 1, 2, 151),
            rec(["1.0+5", "20.0", "0.0", "1.0-3", "2.3-2", "0.0"], 1, 2, 151),
            rec(["", "", "", "", "", ""], 1, 2, 0),
        ]
        .concat()
    }

    #[test]
    fn unresolved_averages_replace_background_with_lssf0() {
        let full = ["1.0-5", "2.0+7"];
        let text = [
            unresolved_mf2(),
            mf3_section(2, full, ["10.0", "10.0"]),
            mf3_section(102, full, ["0.5", "0.5"]),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let tables = reconstruct(&tape, 1, 1.0e-3).unwrap();
        let mf2 = Mf2::parse(&tape.section(1, 2, 151).unwrap()).unwrap();
        let RangeKind::Unresolved(urr) = &mf2.ranges[0].kind else {
            panic!("expected unresolved");
        };
        let xs = UnresolvedXs::new(urr).unwrap();
        let at = |mt: u32, e: f64| tables[&mt].eval(e).unwrap();
        // 格子点（2 keV は補助格子の点）では、平均断面積と MF3 の背景断面積の和
        let average = xs.cross_sections(2.0e3);
        for (mt, expected) in [(2, average.elastic + 10.0), (102, average.capture + 0.5)] {
            assert!(
                (at(mt, 2.0e3) - expected).abs() <= 1e-12 * expected,
                "MT{mt}"
            );
        }
        // 範囲の外は MF3 のまま
        assert_eq!((at(2, 10.0), at(2, 2.0e5)), (10.0, 10.0));
        assert_eq!(at(1, 2.0e3), at(2, 2.0e3) + at(102, 2.0e3));
    }

    #[test]
    fn unresolved_background_at_top_comes_from_inside() {
        // MF3 は非分離共鳴の範囲の中では背景だけで、上端 1e5 eV で範囲の外の値に段差で変わる。
        // 上端の格子点の背景は、範囲の内側の値を使う
        use crate::endf::rec;
        let capture = [
            rec(["1.001000+3", "1.0", "0", "0", "0", "0"], 1, 3, 102),
            rec(["0.0", "0.0", "0", "0", "1", "4"], 1, 3, 102),
            rec(["4", "2", "", "", "", ""], 1, 3, 102),
            rec(["1.0-5", "0.5", "1.0+5", "0.5", "1.0+5", "7.0"], 1, 3, 102),
            rec(["2.0+7", "7.0", "", "", "", ""], 1, 3, 102),
            rec(["", "", "", "", "", ""], 1, 3, 0),
        ]
        .concat();
        let text = [
            unresolved_mf2(),
            mf3_section(2, ["1.0-5", "2.0+7"], ["10.0", "10.0"]),
            capture,
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let tables = reconstruct(&tape, 1, 1.0e-3).unwrap();
        let mf2 = Mf2::parse(&tape.section(1, 2, 151).unwrap()).unwrap();
        let RangeKind::Unresolved(urr) = &mf2.ranges[0].kind else {
            panic!("expected unresolved");
        };
        let average = UnresolvedXs::new(urr)
            .unwrap()
            .cross_sections(1.0e5)
            .capture;
        let t = &tables[&102];
        assert!((left_limit(t, 1.0e5) - (average + 0.5)).abs() <= 1e-12 * average);
        assert_eq!(right_limit(t, 1.0e5), 7.0);
    }

    #[test]
    fn rejects_lssf0_range_overlapping_resolved_range() {
        // 分離共鳴 [1e-5, 2e3] と、下端がその中にある LSSF=0 の非分離共鳴 [1e3, 1e5]。
        // 重なりで分離共鳴を平均で上書きしないよう、エラーにする（NJOY は分離共鳴の上端まで分離共鳴を使う）
        use crate::endf::rec;
        let mf2 = [
            rec(["9.223500+4", "2.360060+2", "0", "0", "1", "0"], 1, 2, 151),
            rec(["9.223500+4", "1.0", "0", "0", "2", "0"], 1, 2, 151),
            rec(["1.0-5", "2.0+3", "1", "3", "0", "1"], 1, 2, 151),
            rec(["0.0", "9.48-1", "0", "0", "1", "0"], 1, 2, 151),
            rec(["2.360060+2", "0.0", "0", "0", "6", "1"], 1, 2, 151),
            rec(
                ["6.674", "0.5", "1.5-3", "2.3-2", "1.0-2", "0.0"],
                1,
                2,
                151,
            ),
            rec(["1.0+3", "1.0+5", "2", "2", "0", "1"], 1, 2, 151),
            rec(["0.0", "0.948", "0", "0", "1", "0"], 1, 2, 151),
            rec(["2.360060+2", "0.0", "0", "0", "1", "0"], 1, 2, 151),
            rec(["0.5", "0.0", "2", "0", "18", "2"], 1, 2, 151),
            rec(["0.0", "0.0", "0.0", "1.0", "0.0", "0.0"], 1, 2, 151),
            rec(["1.0+3", "20.0", "0.0", "1.0-3", "2.3-2", "0.0"], 1, 2, 151),
            rec(["1.0+5", "20.0", "0.0", "1.0-3", "2.3-2", "0.0"], 1, 2, 151),
            rec(["", "", "", "", "", ""], 1, 2, 0),
        ]
        .concat();
        let full = ["1.0-5", "2.0+7"];
        let text = [
            mf2,
            mf3_section(2, full, ["10.0", "10.0"]),
            mf3_section(102, full, ["0.5", "0.5"]),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        assert!(reconstruct(&tape, 1, 1.0e-3).is_err());
    }

    fn single_level_mf2() -> String {
        use crate::endf::rec;
        [
            rec(["9.223500+4", "2.360060+2", "0", "0", "1", "0"], 1, 2, 151),
            rec(["9.223500+4", "1.0", "0", "0", "1", "0"], 1, 2, 151),
            rec(["1.0-5", "1.0+3", "1", "3", "0", "1"], 1, 2, 151),
            rec(["0.0", "9.48-1", "0", "0", "1", "0"], 1, 2, 151),
            rec(["2.360060+2", "0.0", "0", "0", "6", "1"], 1, 2, 151),
            rec(
                ["6.674", "0.5", "1.5-3", "2.3-2", "1.0-2", "0.0"],
                1,
                2,
                151,
            ),
            rec(["", "", "", "", "", ""], 1, 2, 0),
        ]
        .concat()
    }

    #[test]
    fn mt18_is_kept_without_mt19() {
        // NJOY と同じく、MT18 を冗長とみなすのは MT19 がある場合だけ
        let full = ["1.0-5", "2.0+7"];
        let text = [
            single_level_mf2(),
            mf3_section(2, full, ["0.0", "0.0"]),
            mf3_section(18, full, ["1.0", "1.0"]),
            mf3_section(20, ["5.0+6", "2.0+7"], ["0.125", "0.125"]),
            mf3_section(102, full, ["0.0", "0.0"]),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let tables = reconstruct(&tape, 1, 1.0e-3).unwrap();
        assert!(!tables.contains_key(&19));
        // 分離共鳴の外では MF3 の MT18 がそのまま残る
        assert_eq!(tables[&18].eval(1.0e4), Some(1.0));
        // 共鳴の核分裂は MT18 に足される
        assert!(tables[&18].eval(6.674).unwrap() > 2.0);
    }

    #[test]
    fn resonance_fission_goes_to_mt19_when_mt18_is_rebuilt() {
        let full = ["1.0-5", "2.0+7"];
        let text = [
            single_level_mf2(),
            mf3_section(2, full, ["0.0", "0.0"]),
            mf3_section(18, full, ["70.0", "70.0"]),
            mf3_section(19, full, ["1.0", "1.0"]),
            mf3_section(102, full, ["0.0", "0.0"]),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let tables = reconstruct(&tape, 1, 1.0e-3).unwrap();
        assert_eq!(tables[&18].eval(1.0e4), Some(1.0));
        // MT18 は MF3 の値（70 b）ではなく、共鳴を足した MT19 の和
        let fission = tables[&18].eval(6.674).unwrap();
        assert!(fission > 2.0);
        assert_eq!(Some(fission), tables[&19].eval(6.674));
        let parts: f64 = [2, 19, 102]
            .iter()
            .map(|mt| tables[mt].eval(6.674).unwrap())
            .sum();
        assert!((tables[&1].eval(6.674).unwrap() - parts).abs() <= 1e-9 * parts);
    }

    #[test]
    fn add_does_not_duplicate_continuous_points() {
        // 0.2 + (0.9 - 0.2) は 0.9 と 1ulp 違う。表の点の上で補間し直すと偽の不連続点ができる
        let a = linear(&[1.0, 2.0, 3.0], &[0.2, 0.9, 0.9]);
        let b = linear(&[1.0, 3.0], &[0.0, 0.0]);
        let sum = add(&a, &b);
        assert_eq!(sum.x, vec![1.0, 2.0, 3.0]);
        assert_eq!(sum.y, vec![0.2, 0.9, 0.9]);
    }

    #[test]
    fn n2n_is_rebuilt_from_its_partials() {
        // ENDF-6 では MT16 は MT875〜891 の和
        let full = ["1.0-5", "2.0+7"];
        let text = [
            mf3_section(2, full, ["10.0", "10.0"]),
            mf3_section(16, ["1.0+7", "2.0+7"], ["5.0", "5.0"]),
            mf3_section(875, ["1.0+7", "2.0+7"], ["1.0", "1.0"]),
            mf3_section(876, ["1.0+7", "2.0+7"], ["2.0", "2.0"]),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        let tables = reconstruct(&tape, 1, 1.0e-3).unwrap();
        assert_eq!(tables[&16].eval(1.5e7), Some(3.0));
        assert_eq!(tables[&1].eval(1.5e7), Some(13.0));
    }

    #[test]
    fn add_without_jump_at_shared_end_has_no_duplicate() {
        let a = linear(&[1.0, 10.0], &[1.0, 1.0]);
        let b = linear(&[1.0, 5.0], &[0.0, 2.0]);
        let sum = add(&a, &b);
        assert_eq!(sum.x, vec![1.0, 5.0, 5.0, 10.0]);
        assert_eq!(sum.y, vec![1.0, 3.0, 1.0, 1.0]);
    }
}
