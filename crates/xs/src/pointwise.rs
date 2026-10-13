//! 0 K のポイントワイズの断面積（lin-lin で補間できる表）の作成。
//!
//! MF3 の各反応を許容誤差まで線形化し、分離共鳴領域では共鳴パラメータから計算した断面積を
//! 背景断面積に足す（NJOY の RECONR に相当）。

use std::collections::BTreeMap;
use std::convert::Infallible;

use crate::doppler::{broaden, broadening_limit};
use crate::endf::{Error, Interpolation, Mf2, Mf3, RangeKind, Tab1, Tape, Unresolved};
use crate::resonance::{
    MlbwXs, RMatrixLimitedXs, ReichMooreXs, Sigma, UnresolvedXs, unresolved_grid,
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
        let value = |e: f64| [law.interpolate(e, (xk, yk), (x1, table.y[k + 1]))];
        // 右端は表の値をそのまま使う（補間し直すと丸めで値がずれ、偽の不連続点になる）。
        // ヒストグラムだけは区間内の値で終わり、次の点で段差になる
        let right = if law == Interpolation::Histogram {
            value(x1)
        } else {
            [table.y[k + 1]]
        };
        let f = |e: f64| Ok::<_, Infallible>(value(e));
        let Ok(()) = refine(&f, (xk, [yk]), (x1, right), tol, &mut |e, v| {
            x.push(e);
            y.push(v[0]);
        });
    }
    Tab1::linear(x, y)
}

/// 区間 (a, b) を、中点での線形補間の誤差が `tol` 以下になるまで二分し、
/// 加えた点と右端の点を左から順に `push` に渡す（左端は渡さない）。
fn refine<V: AsRef<[f64]> + Clone, E>(
    f: &impl Fn(f64) -> Result<V, E>,
    a: (f64, V),
    b: (f64, V),
    tol: f64,
    push: &mut impl FnMut(f64, V),
) -> Result<(), E> {
    // 右側の区間から積み、左側から取り出す
    let mut stack = vec![b];
    let mut left = a;
    while let Some(right) = stack.pop() {
        let mid = 0.5 * (left.0 + right.0);
        let narrow = right.0 - left.0 <= MIN_RELATIVE_WIDTH * mid.abs();
        if !narrow {
            let exact = f(mid)?;
            let converged = (left.1.as_ref().iter())
                .zip(right.1.as_ref())
                .zip(exact.as_ref())
                .all(|((&l, &r), &exact)| {
                    let linear = 0.5 * (l + r);
                    (linear - exact).abs() <= tol * exact.abs()
                });
            if !converged {
                stack.push(right);
                stack.push((mid, exact));
                continue;
            }
        }
        push(right.0, right.1.clone());
        left = right;
    }
    Ok(())
}

/// [el, eh] での共鳴の寄与 `f`（反応ごとの値）を、反応ごとに lin-lin の表にする。
/// `nodes` は必ず格子に入れるエネルギー（共鳴のピークと半値幅の位置など）。
pub fn resonance_tables<V: AsRef<[f64]> + Clone, E>(
    f: impl Fn(f64) -> Result<V, E>,
    nodes: &[f64],
    el: f64,
    eh: f64,
    tol: f64,
) -> Result<Vec<Tab1>, E> {
    let mut seeds: Vec<f64> = nodes
        .iter()
        .copied()
        .filter(|&e| e > el && e < eh)
        .collect();
    seeds.push(el);
    seeds.push(eh);
    seeds.sort_by(f64::total_cmp);
    seeds.dedup();

    let first = f(el)?;
    let mut x = vec![el];
    let mut values = vec![first.clone()];
    // 区間の左端の値（refine は右端を最後に積むので、次の区間の左端は今の区間の右端）
    let mut left_value = first;
    for w in seeds.windows(2) {
        let right = (w[1], f(w[1])?);
        let next_left = right.1.clone();
        refine(&f, (w[0], left_value), right, tol, &mut |e, v| {
            x.push(e);
            values.push(v);
        })?;
        left_value = next_left;
    }
    let count = values.first().map_or(0, |v| v.as_ref().len());
    Ok((0..count)
        .map(|i| Tab1::linear(x.clone(), values.iter().map(|v| v.as_ref()[i]).collect()))
        .collect())
}

/// `Sigma` の成分の表（`RESONANCE_MTS` の順）を、反応ごとの寄与にする。全断面積は部分反応の和で
/// 作り直すので除き、核分裂は `fission_mt` に入れる。
fn sigma_contributions(tables: Vec<Tab1>, fission_mt: u32) -> Vec<(u32, Tab1)> {
    RESONANCE_MTS
        .into_iter()
        .map(|mt| if mt == 18 { fission_mt } else { mt })
        .zip(tables)
        .filter(|&(mt, _)| mt != 1)
        .collect()
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
    // 上端は右側からの極限なので 0。下端より下は eval が None を返す
    if t.x.last().is_none_or(|&last| e >= last) {
        return 0.0;
    }
    t.eval(e).unwrap_or(0.0)
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

/// 運動学のしきい値 E_th = −Q (A + 1)/A より下から始まる反応の表を、しきい値から始める（[ADR-0100]）。
///
/// しきい値より下の断面積はありえない。評価済みデータには、重心系のしきい値 |Q| と E_th の間（Q 値と AWR の丸めや、
/// しきい値の書き方の違い）から始まる表がある。E_th より下の点は捨て、E_th に 0 を置き（NJOY RECONR と同じく
/// 最初の点をしきい値に移すが、NJOY は最初の点の値を残して段差にする。こちらは断面積がしきい値から連続に立ち上がる
/// よう 0 にする。評価済みデータでは 0 でない値もその反応の最大値の 5e-6 以下）、E_th から次の点までは、しきい値を
/// またいでいた区間の補間則を引き継ぐ（NJOY と同じ）。ただし ln y をとる補間則（INT=4, 5）は 0 の端で成り立たないので lin-lin にする。
/// その先はファイルのまま。ちょうど E_th にある点は残す。`awr` は運動学（2体の散乱行列）と同じ
/// AWR（最初の MF3 の見出し）。
///
/// 次は推測で直さずにエラーにする: AWR が正でない、|Q| より下に値が 0 でない点がある、E_th より上に点がない。
///
/// [ADR-0100]: ../../../docs/adr/0100-move-tables-to-kinematic-thresholds.md
fn start_at_threshold(
    mt: u32,
    mf3: &Mf3,
    awr: f64,
) -> Result<(Tab1, Option<ThresholdShift>), String> {
    let sigma = &mf3.sigma;
    if mf3.qi >= 0.0 {
        return Ok((sigma.clone(), None));
    }
    if awr <= 0.0 || awr.is_nan() {
        return Err(format!(
            "MT{mt}: a threshold reaction needs a positive AWR, got {awr}"
        ));
    }
    let q = -mf3.qi;
    let threshold = q * (awr + 1.0) / awr;
    let Some(&first) = sigma.x.first() else {
        return Ok((sigma.clone(), None));
    };
    if first >= threshold {
        return Ok((sigma.clone(), None));
    }
    if let Some((e, y)) = sigma
        .x
        .iter()
        .zip(&sigma.y)
        .find(|&(&e, &y)| e < q && y != 0.0)
    {
        return Err(format!(
            "MT{mt}: the cross section is {y} b at {e} eV, below even the centre-of-mass threshold \
             |Q| = {q} eV"
        ));
    }
    // E_th より下の点の数
    let below = sigma.x.partition_point(|&e| e < threshold);
    if !sigma.x[below..].iter().any(|&e| e > threshold) {
        return Err(format!(
            "MT{mt}: the table ends at or below the kinematic threshold {threshold} eV"
        ));
    }
    let at_threshold = sigma.x[below] == threshold;
    let (x, y, breaks, laws, dropped) = if at_threshold {
        // ちょうど E_th の点から始める。区間はどれもファイルのまま
        let (breaks, laws) = sigma
            .breaks
            .iter()
            .zip(&sigma.laws)
            .filter_map(|(&b, &law)| {
                let shifted = b.checked_sub(below)?;
                (shifted >= 2).then_some((shifted, law))
            })
            .unzip();
        (
            sigma.x[below..].to_vec(),
            sigma.y[below..].to_vec(),
            breaks,
            laws,
            below,
        )
    } else {
        // E_th に 0 を置く（しきい値から連続に立ち上がる）。区切りの番号は、捨てた点の数から置いた1点を引いた分だけ
        // 詰める
        let mut x = vec![threshold];
        x.extend_from_slice(&sigma.x[below..]);
        let mut y = vec![0.0];
        y.extend_from_slice(&sigma.y[below..]);
        // しきい値をまたいでいた区間（点 below − 1 と below の間）の補間則を引き継ぐ。ln y をとる補間則は 0 の端で
        // 成り立たないので lin-lin にする
        let crossing = sigma
            .breaks
            .iter()
            .position(|&b| b > below)
            .map_or(Interpolation::LinLin, |k| sigma.laws[k]);
        let takes_log_of_y = matches!(crossing, Interpolation::LogLin | Interpolation::LogLog);
        let first_law = if takes_log_of_y && (y[0] <= 0.0 || y[1] <= 0.0) {
            Interpolation::LinLin
        } else {
            crossing
        };
        let (mut breaks, mut laws) = (vec![2], vec![first_law]);
        for (&b, &law) in sigma.breaks.iter().zip(&sigma.laws) {
            if let Some(shifted) = (b + 1).checked_sub(below)
                && shifted >= 3
            {
                breaks.push(shifted);
                laws.push(law);
            }
        }
        (x, y, breaks, laws, below - 1)
    };
    // E_th より下の 0 でない値（捨てるか 0 にした）のうち最も大きいものの、反応の最大値に対する比
    let maximum = sigma.y.iter().copied().fold(0.0, f64::max);
    let below_max = sigma.y[..below].iter().map(|v| v.abs()).fold(0.0, f64::max);
    let zeroed = if maximum > 0.0 {
        below_max / maximum
    } else {
        0.0
    };
    tracing::debug!(
        "MT{mt}: the table starts at {first} eV, below the kinematic threshold {threshold} eV; \
         starting it at the threshold and dropping {dropped} more points below it; the largest \
         value below it, now 0, was {zeroed:.2e} of the maximum"
    );
    let shift = ThresholdShift {
        mt,
        first,
        threshold,
        dropped,
        zeroed,
    };
    Ok((Tab1 { breaks, laws, x, y }, Some(shift)))
}

/// 運動学のしきい値に合わせた1反応（`start_at_threshold`）。
#[derive(Debug, Clone, Copy, PartialEq)]
struct ThresholdShift {
    mt: u32,
    /// 表の最初の点 (eV)
    first: f64,
    /// 運動学のしきい値 (eV)
    threshold: f64,
    /// 捨てた点の数（しきい値に移した最初の点を除く）
    dropped: usize,
    /// しきい値より下の 0 でない値（捨てるか 0 にした）のうち最も大きいものの、反応の最大値に対する比
    zeroed: f64,
}

impl ThresholdShift {
    /// 最初の点を動かした幅 (eV)
    fn shift(&self) -> f64 {
        self.threshold - self.first
    }

    /// しきい値に対する、動かした幅の比
    fn relative(&self) -> f64 {
        self.shift() / self.threshold
    }
}

/// 1核種でしきい値に合わせた反応のまとめ（核種ごとに1回だけ警告するため）。
#[derive(Debug, Clone, Copy, PartialEq)]
struct ThresholdSummary {
    reactions: usize,
    /// 動かした幅 (eV) が最も大きい反応
    largest: ThresholdShift,
    /// しきい値に対する動かした幅の比が最も大きい反応
    largest_relative: ThresholdShift,
    /// 捨てた点の数の和
    dropped: usize,
    /// しきい値より下の 0 でない値を 0 にした反応の数
    zeroed: usize,
    /// そのうち、反応の最大値に対する比が最も大きかった反応
    largest_zeroed: Option<ThresholdShift>,
}

/// しきい値に合わせた反応をまとめる。なければ `None`。
fn summarize_thresholds(shifts: &[ThresholdShift]) -> Option<ThresholdSummary> {
    let largest = *shifts
        .iter()
        .max_by(|a, b| a.shift().total_cmp(&b.shift()))?;
    let largest_relative = *shifts
        .iter()
        .max_by(|a, b| a.relative().total_cmp(&b.relative()))?;
    Some(ThresholdSummary {
        reactions: shifts.len(),
        largest,
        largest_relative,
        dropped: shifts.iter().map(|s| s.dropped).sum(),
        zeroed: shifts.iter().filter(|s| s.zeroed > 0.0).count(),
        largest_zeroed: shifts
            .iter()
            .filter(|s| s.zeroed > 0.0)
            .max_by(|a, b| a.zeroed.total_cmp(&b.zeroed))
            .copied(),
    })
}

impl ColdParts {
    /// MF3 と共鳴パラメータから、許容誤差 `tol` で 0 K の部分反応を作る。
    #[tracing::instrument(name = "reconstruct", level = "debug", skip_all)]
    pub fn read(tape: &Tape, mat: i32, tol: f64) -> Result<Self, Error> {
        let sections = tape.sections(mat, 3);
        let available: Vec<u32> = sections.iter().map(|(mt, _)| *mt).collect();
        let mut tables = BTreeMap::new();
        let mut awr = None;
        let mut shifts = Vec::new();
        // しきい値は、運動学（2体の散乱行列）と同じく最初の MF3 の見出しの AWR で求める
        let head_awr = match sections.first() {
            Some((_, section)) => Mf3::parse(section)?.awr,
            None => 0.0,
        };
        for (mt, section) in sections
            .iter()
            .filter(|(mt, _)| !is_dropped(*mt, &available))
        {
            let mt = *mt;
            let mf3 = Mf3::parse(section)?;
            awr.get_or_insert(mf3.awr);
            let (sigma, shift) = start_at_threshold(mt, &mf3, head_awr)
                .map_err(|message| Error::new(section.reader().line(), message))?;
            shifts.extend(shift);
            tables.insert(mt, linearize(&sigma, tol));
        }
        if let Some(summary) = summarize_thresholds(&shifts) {
            let (largest, relative) = (summary.largest, summary.largest_relative);
            let zeroed = summary.largest_zeroed.map_or(String::new(), |z| {
                format!(
                    "; {} reactions had nonzero values below the threshold, now 0 (the largest \
                     {:.2e} of the reaction maximum in MT{})",
                    summary.zeroed, z.zeroed, z.mt
                )
            });
            tracing::warn!(
                "{} reactions start below their kinematic threshold -Q(A+1)/A and now start at it \
                 with a zero cross section; the largest shift is {:.4e} eV in MT{} and the largest \
                 relative shift is {:.2e} of the threshold in MT{}; {} points below the thresholds \
                 were dropped{zeroed}",
                summary.reactions,
                largest.shift(),
                largest.mt,
                relative.relative(),
                relative.mt,
                summary.dropped
            );
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
        let contributions = match &range.kind {
            RangeKind::ScatteringRadius { .. } => continue,
            RangeKind::Unresolved(urr) if urr.lssf == 1 => continue,
            RangeKind::Unresolved(urr) if urr.lssf == 0 => {
                // 平均で置き換えると、重なった部分の分離共鳴が消える（NJOY は分離共鳴の上端まで分離共鳴を使う）
                let resolved_top = mf2
                    .ranges
                    .iter()
                    .filter(|r| r.kind.is_resolved())
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
            RangeKind::Unresolved(Unresolved { lssf, .. }) => {
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
                let xs = ReichMooreXs::new(rm);
                let f = |e: f64| xs.cross_sections(e).map(|sigma| sigma.components());
                // ENDF の読み込みのエラー（行番号つき）にまとめるので、ここでは文にする
                let tables = resonance_tables(f, &nodes, range.el, range.eh, tol)
                    .map_err(|e| e.to_string())?;
                sigma_contributions(tables, fission_mt)
            }
            RangeKind::RMatrixLimited(rml) => {
                // ENDF の読み込みのエラー（行番号つき）にまとめるので、ここでは文にする
                let xs = RMatrixLimitedXs::new(rml).map_err(|e| e.to_string())?;
                let nodes = resonance_nodes(xs.peaks());
                let f = |e: f64| xs.cross_sections(e);
                let tables = resonance_tables(f, &nodes, range.el, range.eh, tol)
                    .map_err(|e| e.to_string())?;
                // 全断面積は部分反応の和で作り直す
                xs.reactions()
                    .into_iter()
                    .zip(tables)
                    .filter(|&(mt, _)| mt != 1)
                    .collect()
            }
            RangeKind::Mlbw(bw) => {
                let peaks = bw
                    .l_values
                    .iter()
                    .flat_map(|lv| &lv.resonances)
                    .map(|r| (r.er, 0.5 * (r.gn + r.gg + r.gf)));
                // 散乱半径の表の点で位相シフトが折れるので、その点も格子に入れる
                let breaks = bw.ap_energy.iter().flat_map(|t| t.x.iter().copied());
                let nodes: Vec<f64> = resonance_nodes(peaks).into_iter().chain(breaks).collect();
                let xs = MlbwXs::new(bw);
                let f = |e: f64| xs.cross_sections(e).map(|sigma| sigma.components());
                // ENDF の読み込みのエラー（行番号つき）にまとめるので、ここでは文にする
                let tables = resonance_tables(f, &nodes, range.el, range.eh, tol)
                    .map_err(|e| e.to_string())?;
                sigma_contributions(tables, fission_mt)
            }
        };
        for (mt, contribution) in contributions {
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
    let averages: Vec<Sigma> = grid
        .iter()
        .map(|&e| xs.cross_sections(e))
        .collect::<Result<_, _>>()
        .map_err(|e| e.to_string())?;
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
    /// `RESONANCE_MTS` の順の成分
    fn components(&self) -> [f64; 4] {
        [self.total, self.elastic, self.fission, self.capture]
    }
}

#[cfg(test)]
mod tests;
