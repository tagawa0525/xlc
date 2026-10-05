//! 0 K のポイントワイズの断面積（lin-lin で補間できる表）の作成。
//!
//! MF3 の各反応を許容誤差まで線形化し、分離共鳴領域では共鳴パラメータから計算した断面積を
//! 背景断面積に足す（NJOY の RECONR に相当）。

use std::collections::BTreeMap;

use crate::endf::{Error, Tab1, Tape};
use crate::resonance::{ReichMooreXs, Sigma};

/// 線形化・共鳴計算の打ち切り判定: 区間の中点で、線形補間と真値の差が
/// `tol × |真値|` 以下になるまで区間を二分する。
pub const DEFAULT_TOLERANCE: f64 = 1.0e-3;

/// 共鳴の寄与を足す反応（全断面積・弾性散乱・核分裂・捕獲）。
pub const RESONANCE_MTS: [u32; 4] = [1, 2, 18, 102];

/// 表を lin-lin 補間で `tol` の精度になるまで点を加えて線形化する。
pub fn linearize(table: &Tab1, tol: f64) -> Tab1 {
    todo!()
}

/// [el, eh] での共鳴の寄与を、反応ごとに lin-lin の表にする。
pub fn resonance_tables(xs: &ReichMooreXs, nodes: &[f64], el: f64, eh: f64, tol: f64) -> [Tab1; 4] {
    todo!()
}

/// 2つの lin-lin の表の和。`b` は定義域の外で 0 とし、その端は不連続点（x の重複）で表す。
pub fn add(a: &Tab1, b: &Tab1) -> Tab1 {
    todo!()
}

/// テープ中の1核種の、全反応の 0 K のポイントワイズの断面積。キーは MT。
pub fn reconstruct(tape: &Tape, mat: i32, tol: f64) -> Result<BTreeMap<u32, Tab1>, Error> {
    todo!()
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
