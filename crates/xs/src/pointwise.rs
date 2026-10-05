//! 0 K のポイントワイズの断面積（lin-lin で補間できる表）の作成。
//!
//! MF3 の各反応を許容誤差まで線形化し、分離共鳴領域では共鳴パラメータから計算した断面積を
//! 背景断面積に足す（NJOY の RECONR に相当）。

use std::collections::BTreeMap;

use crate::endf::{Error, Mf2, Mf3, RangeKind, ReichMoore, Tab1, Tape};
use crate::resonance::{ReichMooreXs, Sigma};

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
        let right = f(x1);
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
pub fn resonance_tables(xs: &ReichMooreXs, nodes: &[f64], el: f64, eh: f64, tol: f64) -> [Tab1; 4] {
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
fn left_limit(t: &Tab1, e: f64) -> f64 {
    let (first, last) = (t.x[0], t.x[t.x.len() - 1]);
    if e <= first || e > last {
        return 0.0;
    }
    let k = t.x.partition_point(|&xi| xi < e);
    let (x0, x1) = (t.x[k - 1], t.x[k]);
    let (y0, y1) = (t.y[k - 1], t.y[k]);
    y0 + (y1 - y0) * (e - x0) / (x1 - x0)
}

/// lin-lin の表の、e の右側からの極限値（定義域の外は 0）。
fn right_limit(t: &Tab1, e: f64) -> f64 {
    let (first, last) = (t.x[0], t.x[t.x.len() - 1]);
    if e < first || e >= last {
        return 0.0;
    }
    t.eval(e).expect("inside the table range")
}

/// テープ中の1核種の、全反応の 0 K のポイントワイズの断面積。キーは MT。
pub fn reconstruct(tape: &Tape, mat: i32, tol: f64) -> Result<BTreeMap<u32, Tab1>, Error> {
    let mut tables = BTreeMap::new();
    for mt in tape.mts(mat, 3) {
        let section = tape.section(mat, 3, mt).expect("listed by mts");
        tables.insert(mt, linearize(&Mf3::parse(&section)?.sigma, tol));
    }
    let Some(section) = tape.section(mat, 2, 151) else {
        return Ok(tables);
    };
    let mf2 = Mf2::parse(&section)?;
    for range in &mf2.ranges {
        match &range.kind {
            RangeKind::ScatteringRadius { .. } => {}
            RangeKind::Unresolved { lssf: 1 } => {}
            RangeKind::Unresolved { lssf } => {
                return Err(Error::new(
                    section.reader().line(),
                    format!("unresolved range with LSSF={lssf} is not supported"),
                ));
            }
            RangeKind::ReichMoore(rm) => {
                let nodes = resonance_nodes(rm);
                let xs = ReichMooreXs::new(rm);
                let contributions = resonance_tables(&xs, &nodes, range.el, range.eh, tol);
                for (mt, contribution) in RESONANCE_MTS.into_iter().zip(contributions) {
                    if mt == 18 && contribution.y.iter().all(|&v| v == 0.0) {
                        continue;
                    }
                    let combined = match tables.get(&mt) {
                        Some(background) => add(background, &contribution),
                        None => contribution,
                    };
                    tables.insert(mt, combined);
                }
            }
        }
    }
    Ok(tables)
}

/// 共鳴のピークと、その半値幅の位置。
fn resonance_nodes(rm: &ReichMoore) -> Vec<f64> {
    rm.l_values
        .iter()
        .flat_map(|lv| &lv.resonances)
        .flat_map(|r| {
            let half_width = 0.5 * (r.gn + r.gg + r.gfa.abs() + r.gfb.abs());
            [r.er - half_width, r.er, r.er + half_width]
        })
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
    use crate::endf::Interpolation;
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

    #[test]
    fn add_without_jump_at_shared_end_has_no_duplicate() {
        let a = linear(&[1.0, 10.0], &[1.0, 1.0]);
        let b = linear(&[1.0, 5.0], &[0.0, 2.0]);
        let sum = add(&a, &b);
        assert_eq!(sum.x, vec![1.0, 5.0, 5.0, 10.0]);
        assert_eq!(sum.y, vec![1.0, 3.0, 1.0, 1.0]);
    }
}
