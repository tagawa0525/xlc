use super::*;
use crate::endf::ResonanceRange;
use crate::resonance::r_matrix_limited_tests::inelastic_level;
use crate::resonance::reich_moore_tests::{AWRI, single_level, single_level_exact};

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
    let f = |e: f64| xs.cross_sections(e).map(|sigma| sigma.components());
    let tables = resonance_tables(f, &[6.674], 1.0e-5, 1.0e3, 1.0e-3).unwrap();
    for (i, (t, mt)) in tables.iter().zip(RESONANCE_MTS).enumerate() {
        assert_linear(t);
        assert_eq!((t.x[0], *t.x.last().unwrap()), (1.0e-5, 1.0e3));
        for w in t.x.windows(2) {
            let mid = 0.5 * (w[0] + w[1]);
            let exact = single_level_exact(mid, gn, gg, gf).components()[i];
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

/// 弾性散乱の MF3（しきい反応の `threshold_section` と同じ AWR 9。しきい値は最初の MF3 の AWR で求める）
fn elastic_with_awr_9() -> String {
    use crate::endf::rec;
    [
        rec(["4.009000+3", "9.0", "0", "0", "0", "0"], 1, 3, 2),
        rec(["0.0", "0.0", "0", "0", "1", "2"], 1, 3, 2),
        rec(["2", "2", "", "", "", ""], 1, 3, 2),
        rec(["1.0-5", "10.0", "2.0+7", "10.0", "", ""], 1, 3, 2),
        rec(["", "", "", "", "", ""], 1, 3, 0),
    ]
    .concat()
}

/// しきい反応の MF3（AWR 9、Q = −1000 eV なので運動学のしきい値は 1000 × 10 / 9 eV、lin-lin）
fn threshold_section(points: &[(&str, &str)]) -> String {
    threshold_section_with("9.0", "-1.0+3", "2", points)
}

/// しきい反応の MF3（AWR `awr`、Q 値 `q`、補間則 `law`）
fn threshold_section_with(awr: &str, q: &str, law: &str, points: &[(&str, &str)]) -> String {
    use crate::endf::rec;
    let np = points.len().to_string();
    let mut lines = vec![
        rec(["4.009000+3", awr, "0", "0", "0", "0"], 1, 3, 51),
        rec([q, q, "0", "0", "1", &np], 1, 3, 51),
        rec([&np, law, "", "", "", ""], 1, 3, 51),
    ];
    for chunk in points.chunks(3) {
        let mut fields = ["", "", "", "", "", ""];
        for (k, (x, y)) in chunk.iter().enumerate() {
            fields[2 * k] = x;
            fields[2 * k + 1] = y;
        }
        lines.push(rec(fields, 1, 3, 51));
    }
    lines.push(rec(["", "", "", "", "", ""], 1, 3, 0));
    lines.concat()
}

/// `threshold_section_with` の表を `start_at_threshold` に通す
fn adjusted(
    awr: &str,
    q: &str,
    law: &str,
    points: &[(&str, &str)],
) -> Result<(Tab1, Option<ThresholdShift>), String> {
    let tape = Tape::parse(&threshold_section_with(awr, q, law, points)).unwrap();
    let mf3 = Mf3::parse(&tape.section(1, 3, 51).unwrap()).unwrap();
    start_at_threshold(51, &mf3, mf3.awr)
}

#[test]
fn interval_from_the_threshold_keeps_laws_valid_at_zero() {
    // しきい値から次の点までは、ファイルの区間の補間則を引き継ぐ（NJOY と同じ）。y に線形な補間則
    // （lin-log の INT=3）は 0 の端でも成り立つので、区間が長くても形を変えない
    let threshold: f64 = 1.0e3 * 10.0 / 9.0;
    let points = [("1.0+3", "0.0"), ("2.0+3", "1.0"), ("2.0+7", "1.0")];
    let (table, _) = adjusted("9.0", "-1.0+3", "3", &points).unwrap();
    assert_eq!(table.x[0], threshold);
    let e = (threshold * 2.0e3).sqrt();
    let expected = (e / threshold).ln() / (2.0e3 / threshold).ln();
    assert!((table.eval(e).unwrap() - expected).abs() < 1e-12);
    // ln y の補間則（INT=4, 5）は 0 の端で成り立たないので、その区間だけ lin-lin にする
    let (table, _) = adjusted("9.0", "-1.0+3", "5", &points).unwrap();
    assert_eq!(table.laws[0], Interpolation::LinLin);
    let mid = 0.5 * (threshold + 2.0e3);
    assert!((table.eval(mid).unwrap() - 0.5).abs() < 1e-12);
    assert_eq!(table.laws.last(), Some(&Interpolation::LogLog));
}

#[test]
fn point_at_the_threshold_is_kept() {
    // Q = −900 eV、AWR 9 でしきい値はちょうど 1000 eV。しきい値の点は捨てず、その値も変えない
    let points = [("9.5+2", "0.0"), ("1.0+3", "0.2"), ("2.0+3", "1.0")];
    let (table, shift) = adjusted("9.0", "-9.0+2", "2", &points).unwrap();
    assert_eq!(table.x, vec![1.0e3, 2.0e3]);
    assert_eq!(table.y, vec![0.2, 1.0]);
    assert_eq!(shift.map(|s| s.dropped), Some(1));
}

#[test]
fn nonzero_values_below_q_are_an_error() {
    // 表が重心系のしきい値 |Q| より下から値を持つのは、しきい値の書き方の違いでは説明できない
    let points = [("5.0+2", "0.1"), ("2.0+3", "1.0")];
    assert!(adjusted("9.0", "-1.0+3", "2", &points).is_err());
}

#[test]
fn tables_without_points_above_the_threshold_are_an_error() {
    let points = [("1.0+3", "0.0"), ("1.1+3", "0.1")];
    assert!(adjusted("9.0", "-1.0+3", "2", &points).is_err());
}

#[test]
fn thresholds_need_a_positive_awr() {
    let points = [("1.0+3", "0.0"), ("2.0+3", "1.0")];
    assert!(adjusted("0.0", "-1.0+3", "2", &points).is_err());
}

#[test]
fn reactions_start_at_the_kinematic_threshold() {
    // 運動学のしきい値 −Q (A + 1)/A より下の断面積はありえないので、最初の点をしきい値に上げる
    // （NJOY RECONR と同じ）
    let threshold = 1.0e3 * 10.0 / 9.0;
    let text = [
        elastic_with_awr_9(),
        threshold_section(&[("1.0+3", "0.0"), ("2.0+3", "1.0"), ("2.0+7", "1.0")]),
    ]
    .concat();
    let tables = reconstruct(&Tape::parse(&text).unwrap(), 1, 1.0e-3).unwrap();
    let level = &tables[&51];
    assert_eq!(level.x[0], threshold);
    assert_eq!(level.eval(1.05e3), None);
    let mid = 0.5 * (threshold + 2.0e3);
    let expected = 0.5;
    assert!((level.eval(mid).unwrap() - expected).abs() < 1e-12);
}

#[test]
fn points_below_the_threshold_are_dropped() {
    // しきい値より下の点が2つ以上あれば、2つ目からは捨てる（NJOY はしきい値のすぐ上に詰めて並べる）
    let threshold = 1.0e3 * 10.0 / 9.0;
    let text = [
        elastic_with_awr_9(),
        threshold_section(&[
            ("1.0+3", "0.0"),
            ("1.1+3", "0.01"),
            ("2.0+3", "1.0"),
            ("2.0+7", "1.0"),
        ]),
    ]
    .concat();
    let tables = reconstruct(&Tape::parse(&text).unwrap(), 1, 1.0e-3).unwrap();
    let level = &tables[&51];
    assert_eq!(level.x[0], threshold);
    assert_eq!(level.y[0], 0.0);
    let mid = 0.5 * (threshold + 2.0e3);
    assert!((level.eval(mid).unwrap() - 0.5).abs() < 1e-12);
}

#[test]
fn nonzero_values_below_the_threshold_become_zero() {
    // しきい値より下の断面積はなく、しきい値から連続に立ち上がるので、しきい値の値は 0 にする（NJOY は
    // 最初の点の値を残して段差にする）。しきい値より上のデータは変えない
    let threshold: f64 = 1.0e3 * 10.0 / 9.0;
    let points = [
        ("1.0+3", "2.0-3"),
        ("1.05+3", "3.0-3"),
        ("2.0+3", "1.0"),
        ("2.0+7", "1.0"),
    ];
    let (table, shift) = adjusted("9.0", "-1.0+3", "2", &points).unwrap();
    assert_eq!((table.x[0], table.y[0]), (threshold, 0.0));
    assert_eq!(&table.y[1..], &[1.0, 1.0]);
    let mid = 0.5 * (threshold + 2.0e3);
    assert!((table.eval(mid).unwrap() - 0.5).abs() < 1e-12);
    // 0 にした値のうち最も大きいもの（反応の最大値に対する比）
    assert_eq!(shift.map(|s| s.zeroed), Some(3.0e-3));
}

#[test]
fn threshold_shift_is_reported_per_reaction() {
    let text = threshold_section(&[
        ("1.0+3", "0.0"),
        ("1.1+3", "0.01"),
        ("2.0+3", "1.0"),
        ("2.0+7", "1.0"),
    ]);
    let tape = Tape::parse(&text).unwrap();
    let mf3 = Mf3::parse(&tape.section(1, 3, 51).unwrap()).unwrap();
    let (_, shift) = start_at_threshold(51, &mf3, mf3.awr).unwrap();
    let threshold = 1.0e3 * 10.0 / 9.0;
    assert_eq!(
        shift,
        Some(ThresholdShift {
            mt: 51,
            first: 1.0e3,
            threshold,
            dropped: 1,
            // しきい値より下の 0.01 b を 0 にした（反応の最大値 1 b に対する比）
            zeroed: 0.01,
        })
    );
    // しきい値より上から始まる表は動かさない
    let text = threshold_section(&[("2.0+3", "0.0"), ("2.0+7", "1.0")]);
    let tape = Tape::parse(&text).unwrap();
    let mf3 = Mf3::parse(&tape.section(1, 3, 51).unwrap()).unwrap();
    assert_eq!(start_at_threshold(51, &mf3, mf3.awr).unwrap().1, None);
}

#[test]
fn threshold_shifts_are_summarized_per_nuclide() {
    // 核種ごとに1回だけ警告する: 合わせた反応の数、最も大きい移動（eV）とその MT、捨てた点の数の和
    assert_eq!(summarize_thresholds(&[]), None);
    let shift = |mt: u32, first: f64, threshold: f64, dropped: usize, zeroed: f64| ThresholdShift {
        mt,
        first,
        threshold,
        dropped,
        zeroed,
    };
    let shifts = [
        shift(51, 1.0e3, 1.1e3, 0, 0.0),
        shift(16, 9.0e6, 9.0005e6, 2, 2.0e-7),
        shift(52, 2.0e3, 2.05e3, 1, 1.0e-6),
    ];
    let summary = summarize_thresholds(&shifts).unwrap();
    assert_eq!(summary.reactions, 3);
    assert_eq!(summary.dropped, 3);
    // eV で最も大きいのは MT16 の 500 eV（相対では MT51 の 1/11 のほうが大きい）
    assert_eq!(summary.largest, shifts[1]);
    assert!((summary.largest.relative() - 500.0 / 9.0005e6).abs() < 1e-15);
    // 相対で最も大きいのは MT51 の 100 / 1100
    assert_eq!(summary.largest_relative, shifts[0]);
    // しきい値より下の 0 でない値を 0 にした反応の数と、反応の最大値に対して最も大きかった値
    assert_eq!(summary.zeroed, 2);
    assert_eq!(summary.largest_zeroed, Some(shifts[2]));
    let none = summarize_thresholds(&shifts[..1]).unwrap();
    assert_eq!((none.zeroed, none.largest_zeroed), (0, None));
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
    let average = xs.cross_sections(2.0e3).unwrap();
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
        .unwrap()
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

#[test]
fn r_matrix_limited_adds_level_inelastic_scattering() {
    let mf2 = Mf2 {
        za: 0.0,
        awr: AWRI,
        ranges: vec![ResonanceRange {
            el: 1.0,
            eh: 5.0e3,
            kind: RangeKind::RMatrixLimited(inelastic_level(0.5)),
        }],
    };
    let mut tables = BTreeMap::from([(2, linear(&[1.0, 5.0e3], &[1.0, 1.0]))]);
    add_resonances(&mut tables, &mf2, 18, 1.0e-3).unwrap();
    // 捕獲のない共鳴でも、共鳴エネルギー（格子点）で止まらない
    let mut capture_free = mf2.clone();
    let RangeKind::RMatrixLimited(rml) = &mut capture_free.ranges[0].kind else {
        panic!("expected R-matrix limited");
    };
    *rml = inelastic_level(0.0);
    add_resonances(&mut BTreeMap::new(), &capture_free, 18, 1.0e-3).unwrap();
    // 非弾性散乱は閾値より上でだけ値を持つ
    let threshold = 1.0e3 * (AWRI + 1.0) / AWRI;
    let level = &tables[&51];
    assert_eq!(level.eval(0.9 * threshold), Some(0.0));
    assert!(level.eval(2.0e3).unwrap() > 0.0);
    assert!(tables[&102].eval(2.0e3).unwrap() > 0.0);
    assert!(tables[&2].eval(2.0e3).unwrap() > 1.0);
}

#[test]
fn mlbw_tables_have_points_at_the_breaks_of_the_radius_table() {
    use crate::endf::rec;
    // AP(E) は 25 eV で折れる。位相シフトも折れるので、中点の判定だけでは折れ目の点が入らないことがある
    let full = ["1.0-5", "2.0+7"];
    let mf2: String = [
        ["1.001000+3", "1.0", "0", "0", "1", "0"],
        ["1.001000+3", "1.0", "0", "0", "1", "0"],
        ["1.0-5", "5.0+1", "1", "2", "1", "0"],
        ["0.0", "0.0", "0", "0", "1", "3"],
        ["3", "2", "", "", "", ""],
        ["1.0-5", "0.70", "2.5+1", "0.60", "5.0+1", "0.70"],
        ["0.0", "0.0", "0", "0", "1", "0"],
        ["1.0", "0.0", "0", "0", "6", "1"],
        ["1.0+4", "0.5", "2.0-9", "1.0-9", "1.0-9", "0.0"],
        ["", "", "", "", "", ""],
    ]
    .iter()
    .enumerate()
    .map(|(i, f)| {
        if i == 9 {
            rec(*f, 1, 2, 0)
        } else {
            rec(*f, 1, 2, 151)
        }
    })
    .collect();
    let text = [
        mf2,
        mf3_section(1, full, ["1.0", "1.0"]),
        mf3_section(2, full, ["1.0", "1.0"]),
        mf3_section(102, full, ["0.0", "0.0"]),
    ]
    .concat();
    let tables = reconstruct(&Tape::parse(&text).unwrap(), 1, 1.0e-3).unwrap();
    assert!(tables[&2].x.contains(&25.0), "{:?}", tables[&2].x);
}

#[test]
fn rebuilt_sums_follow_the_chain_reaction_ranges() {
    // 鎖の反応の部分の MT は、再構成が代表の MT に作り直す和の部分と同じ（xs は代表の MT だけを読む。ADR-0172）
    use xlc_format::library::CHAIN_REACTIONS;
    for reaction in &CHAIN_REACTIONS {
        let partials: Vec<u32> = (1..1000)
            .filter(|&mt| mt != reaction.mt && reaction.includes(mt))
            .collect();
        if partials.is_empty() {
            assert!(
                SUMS.iter().all(|s| s.mt != reaction.mt),
                "MT{}",
                reaction.mt
            );
            continue;
        }
        let sum = SUMS.iter().find(|s| s.mt == reaction.mt).unwrap();
        let parts: Vec<u32> = (1..1000).filter(|&mt| (sum.is_part)(mt)).collect();
        assert_eq!(parts, partials, "MT{}", reaction.mt);
        assert!(
            partials.iter().all(|&mt| (sum.trigger)(mt)),
            "MT{}",
            reaction.mt
        );
    }
}
