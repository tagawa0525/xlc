use super::*;
use crate::endf::ResonanceRange;
use crate::resonance::r_matrix_limited_tests::{fission_level, inelastic_level};
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

/// MF3 の1反応（Q = 0、lin-lin、点は `points`）
fn mf3_points(mt: u32, points: &[(&str, &str)]) -> String {
    use crate::endf::rec;
    let np = points.len().to_string();
    let mut lines = vec![
        rec(["1.001000+3", "1.0", "0", "0", "0", "0"], 1, 3, mt),
        rec(["0.0", "0.0", "0", "0", "1", &np], 1, 3, mt),
        rec([&np, "2", "", "", "", ""], 1, 3, mt),
    ];
    for chunk in points.chunks(3) {
        let mut fields = ["", "", "", "", "", ""];
        for (k, (x, y)) in chunk.iter().enumerate() {
            fields[2 * k] = x;
            fields[2 * k + 1] = y;
        }
        lines.push(rec(fields, 1, 3, mt));
    }
    lines.push(rec(["", "", "", "", "", ""], 1, 3, 0));
    lines.concat()
}

/// 弾性散乱が一定の 10 b で、捕獲が `capture` の点を持つテープ
fn capture_tape(capture: &[(&str, &str)]) -> Tape {
    let text = [
        mf3_section(2, ["1.0-5", "2.0+7"], ["10.0", "10.0"]),
        mf3_points(102, capture),
    ]
    .concat();
    Tape::parse(&text).unwrap()
}

/// 1 eV と 2 eV の 2 点が負の捕獲（零点は 0.5 eV と 2 + 2/3 eV）
const NEGATIVE_CAPTURE: [(&str, &str); 5] = [
    ("1.0-5", "1.0"),
    ("1.0", "-1.0"),
    ("2.0", "-2.0"),
    ("3.0", "1.0"),
    ("2.0+7", "1.0"),
];

fn declared(mt: u32, from: f64, to: f64, take_from: Option<u32>) -> NegativeToZero {
    NegativeToZero {
        mt,
        from,
        to,
        take_from,
    }
}

/// `ColdParts::read` を呼んだスレッドの warn の数
fn read_with_warnings(
    tape: &Tape,
    negatives: &[NegativeToZero],
) -> (Result<ColdParts, Error>, usize) {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use tracing_subscriber::Registry;
    use tracing_subscriber::layer::{Context, Layer, SubscriberExt};

    struct Count(Arc<AtomicUsize>);
    impl<S: tracing::Subscriber> Layer<S> for Count {
        fn on_event(&self, event: &tracing::Event<'_>, _: Context<'_, S>) {
            if *event.metadata().level() == tracing::Level::WARN {
                self.0.fetch_add(1, Ordering::Relaxed);
            }
        }
    }
    let count = Arc::new(AtomicUsize::new(0));
    let subscriber = Registry::default().with(Count(count.clone()));
    let result = tracing::subscriber::with_default(subscriber, || {
        ColdParts::read(tape, 1, 1.0e-3, negatives)
    });
    (result, count.load(Ordering::Relaxed))
}

#[test]
fn negative_cross_sections_are_an_error_at_any_energy() {
    // ドップラー拡がりの格子の判定は正の値を前提にし、負の区間を際限なく分ける。広げる上限（分離共鳴の無い核種では 6.5 MeV）より
    // 上の値も、拡がりの積分の窓と群の平均に入るので、エネルギーによらず 0 K の再構成で止める（ADR-0267）
    for capture in [
        [("1.0-5", "-1.0"), ("2.0+7", "1.0")],
        [("1.0-5", "3.0"), ("7.0+6", "-1.0")],
        [("1.0-5", "1.0"), ("2.0+7", "-1.0")],
    ] {
        let error = ColdParts::read(&capture_tape(&capture), 1, 1.0e-3, &[]).unwrap_err();
        assert!(error.message.contains("MT102"), "{error}");
        assert!(error.message.contains("negative-to-zero"), "{error}");
    }
    // 負の反応がいくつあっても、1 つの文にまとめて示す
    let text = [
        mf3_section(2, ["1.0-5", "2.0+7"], ["-10.0", "10.0"]),
        mf3_points(102, &NEGATIVE_CAPTURE),
    ]
    .concat();
    let error = ColdParts::read(&Tape::parse(&text).unwrap(), 1, 1.0e-3, &[]).unwrap_err();
    assert!(
        error.message.contains("MT2 ") && error.message.contains("MT102 "),
        "{error}"
    );
}

#[test]
fn declared_negative_values_are_raised_to_zero_with_a_warning() {
    // 宣言した範囲の負の点を 0 にし、符号の変わる区間に零点を足す（元の lin-lin の表を 0 で切った関数）。範囲の外の点はそのまま
    let tape = capture_tape(&NEGATIVE_CAPTURE);
    let (cold, warnings) = read_with_warnings(&tape, &[declared(102, 1.0e-5, 10.0, None)]);
    let tables = cold.unwrap().at(0.0);
    assert_eq!(warnings, 1);
    let capture = &tables[&102];
    assert!(capture.y.iter().all(|&v| v >= 0.0));
    let original = linear(
        &[1.0e-5, 1.0, 2.0, 3.0, 2.0e7],
        &[1.0, -1.0, -2.0, 1.0, 1.0],
    );
    let midpoints: Vec<f64> = capture.x.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect();
    for &e in capture.x.iter().chain(&midpoints) {
        let expected = original.eval(e).unwrap().max(0.0);
        let actual = capture.eval(e).unwrap();
        assert!(
            (actual - expected).abs() <= 1.0e-12,
            "at {e} eV: {actual} vs {expected}"
        );
    }
    assert_eq!(capture.eval(3.0), Some(1.0));
    // 全断面積は 0 にした捕獲との和
    assert_eq!(tables[&1].eval(2.0), Some(10.0));
    // 宣言しない反応は変えない
    assert_eq!(tables[&2], linear(&[1.0e-5, 2.0e7], &[10.0, 10.0]));
}

/// 宣言で `ColdParts::read` が止まったときの文
fn declaration_error(tape: &Tape, negatives: &[NegativeToZero]) -> String {
    ColdParts::read(tape, 1, 1.0e-3, negatives)
        .map(drop)
        .unwrap_err()
        .message
}

#[test]
fn declarations_that_do_not_match_the_tables_are_an_error() {
    let tape = capture_tape(&NEGATIVE_CAPTURE);
    let error = |negatives: &[NegativeToZero]| declaration_error(&tape, negatives);
    // 範囲に負の値が無い宣言は止める（評価が替わって宣言が古くなったことに気付けるように）
    let message = error(&[
        declared(102, 1.0e-5, 10.0, None),
        declared(102, 5.0, 10.0, None),
    ]);
    assert!(
        message.contains("MT102 at 5..10 eV: no value is negative"),
        "{message}"
    );
    // 宣言の範囲の外に残る負の値は、宣言しないときと同じく止める
    let message = error(&[declared(102, 1.5, 10.0, None)]);
    assert!(
        message.contains("the 0 K cross sections are negative: MT102"),
        "{message}"
    );
    // 表の無い反応、有限でないか逆転した範囲、引く先の反応が無いか同じ反応の宣言も、それぞれの理由で止める
    for (negatives, reason) in [
        (
            declared(103, 1.0e-5, 10.0, None),
            "the evaluation has no MT103",
        ),
        (
            declared(102, 10.0, 1.0e-5, None),
            "the range must be finite",
        ),
        (
            declared(102, f64::NAN, 10.0, None),
            "the range must be finite",
        ),
        (declared(102, 1.0e-5, 10.0, Some(4)), "no MT4 to take from"),
        (
            declared(102, 1.0e-5, 10.0, Some(102)),
            "must name another reaction",
        ),
    ] {
        let message = error(std::slice::from_ref(&negatives));
        assert!(message.contains(reason), "{negatives:?}: {message}");
    }
}

#[test]
fn taking_the_raised_amount_must_keep_the_total_and_the_source_non_negative() {
    // 全断面積に入らない反応（生成断面積など）との間で引くと、全断面積が変わる
    let text = [
        mf3_section(2, ["1.0-5", "2.0+7"], ["10.0", "10.0"]),
        mf3_points(102, &NEGATIVE_CAPTURE),
        mf3_section(203, ["1.0-5", "2.0+7"], ["10.0", "10.0"]),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let message = declaration_error(&tape, &[declared(102, 1.0e-5, 10.0, Some(203))]);
    assert!(message.contains("part of the total"), "{message}");
    // 引いた先が負になる宣言は、引く先の誤りとして止める（引いた先を別に 0 にすると、全断面積が戻らない）
    let text = [
        mf3_section(2, ["1.0-5", "2.0+7"], ["1.0", "1.0"]),
        mf3_points(102, &NEGATIVE_CAPTURE),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let message = declaration_error(&tape, &[declared(102, 1.0e-5, 10.0, Some(2))]);
    assert!(message.contains("makes MT2 negative"), "{message}");
}

#[test]
fn declarations_cannot_reach_the_background_of_an_lssf0_unresolved_range() {
    // LSSF=0 の非分離共鳴の自己遮蔽は、テープの MF3 のバックグラウンドを読むので、その範囲（1 keV〜100 keV）に宣言は届かない
    let text = [
        unresolved_mf2(),
        mf3_section(2, ["1.0-5", "2.0+7"], ["10.0", "10.0"]),
        mf3_points(102, &NEGATIVE_CAPTURE),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let message = declaration_error(&tape, &[declared(102, 1.0e-5, 2.0e3, None)]);
    assert!(message.contains("unresolved range"), "{message}");
    // 宣言の範囲は両端を含むので、下端 1 keV で終わる宣言も重なる
    let message = declaration_error(&tape, &[declared(102, 1.0e-5, 1.0e3, None)]);
    assert!(message.contains("unresolved range"), "{message}");
    // 共鳴の反応でなくても、和（全断面積など）が自己遮蔽のバックグラウンドと食い違うので止める
    let text = [
        unresolved_mf2(),
        mf3_section(2, ["1.0-5", "2.0+7"], ["10.0", "10.0"]),
        mf3_section(102, ["1.0-5", "2.0+7"], ["1.0", "1.0"]),
        mf3_points(16, &[("1.0-5", "1.0"), ("2.0+3", "-1.0"), ("2.0+7", "1.0")]),
    ]
    .concat();
    let message = declaration_error(
        &Tape::parse(&text).unwrap(),
        &[declared(16, 1.0e3, 1.0e4, None)],
    );
    assert!(message.contains("unresolved range"), "{message}");
    // 範囲の外の宣言は当てる
    assert!(ColdParts::read(&tape, 1, 1.0e-3, &[declared(102, 1.0e-5, 10.0, Some(2))]).is_ok());
}

#[test]
fn declarations_apply_in_the_same_way_in_any_order() {
    // 引く先（MT2）の負の点を 0 にする宣言と、MT2 から引く宣言は、並べた順によらず同じ結果になる（引かない宣言を先に当てる）。
    // 2.5 eV の MT2 は 0 にした後に捕獲の 0.5 b を引かれて負になるので、どちらの順でも引いたせいで負になると示して止まる
    let text = [
        mf3_points(
            2,
            &[
                ("1.0-5", "10.0"),
                ("2.4", "10.0"),
                ("2.5", "-0.1"),
                ("2.6", "10.0"),
                ("2.0+7", "10.0"),
            ],
        ),
        mf3_points(102, &NEGATIVE_CAPTURE),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let (raise, take) = (
        declared(2, 2.45, 2.55, None),
        declared(102, 1.0e-5, 10.0, Some(2)),
    );
    let first = declaration_error(&tape, &[raise.clone(), take.clone()]);
    let second = declaration_error(&tape, &[take.clone(), raise]);
    assert_eq!(first, second);
    assert!(first.contains("makes MT2 negative"), "{first}");
    // 引く先の負の値を宣言しなければ、引いたせいでなく、引く先がもともと負で引く分を出せないことを示して止める
    let message = declaration_error(&tape, &[take]);
    assert!(message.contains("MT2 is negative where"), "{message}");
}

#[test]
fn a_take_from_and_a_declaration_of_its_source_elsewhere_give_the_same_tables_in_any_order() {
    // MT2 の 100 eV の負の点を 0 にする宣言と、10 eV までの捕獲の分を MT2 から引く宣言は、範囲が重ならないので両方とも当たり、
    // 並べた順によらず同じ表になる
    let text = [
        mf3_points(
            2,
            &[
                ("1.0-5", "10.0"),
                ("90.0", "10.0"),
                ("100.0", "-1.0"),
                ("110.0", "10.0"),
                ("2.0+7", "10.0"),
            ],
        ),
        mf3_points(102, &NEGATIVE_CAPTURE),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let (raise, take) = (
        declared(2, 50.0, 200.0, None),
        declared(102, 1.0e-5, 10.0, Some(2)),
    );
    let first = ColdParts::read(&tape, 1, 1.0e-3, &[raise.clone(), take.clone()]).unwrap();
    let second = ColdParts::read(&tape, 1, 1.0e-3, &[take, raise]).unwrap();
    assert_eq!(first.at(0.0), second.at(0.0));
    assert_eq!(first.at(0.0)[&2].eval(100.0), Some(0.0));
}

#[test]
fn a_take_from_cannot_take_from_a_reaction_that_takes_from_another() {
    // 引く先が別の `take-from` の宣言の反応なら、当てる順で引く先の値が変わるので止める
    let text = [
        mf3_points(
            2,
            &[
                ("1.0-5", "10.0"),
                ("90.0", "10.0"),
                ("100.0", "-1.0"),
                ("110.0", "10.0"),
                ("2.0+7", "10.0"),
            ],
        ),
        mf3_points(102, &NEGATIVE_CAPTURE),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let (a, b) = (
        declared(102, 1.0e-5, 10.0, Some(2)),
        declared(2, 50.0, 200.0, Some(102)),
    );
    let first = declaration_error(&tape, &[a.clone(), b.clone()]);
    let second = declaration_error(&tape, &[b, a]);
    assert!(first.contains("also declared with `take-from`"), "{first}");
    assert!(
        second.contains("also declared with `take-from`"),
        "{second}"
    );
}

#[test]
fn the_raised_amount_is_taken_from_the_declared_reaction() {
    // `take-from` の反応から、0 にするために足した分を引くので、全断面積（部分反応の和）は宣言の前と同じ。ENDF/B-VII.1 の Y-89 は、
    // VII.0 の負の MT91 を 0 にし、全断面積を変えずに弾性散乱を減らした
    let tape = capture_tape(&NEGATIVE_CAPTURE);
    let tables = ColdParts::read(&tape, 1, 1.0e-3, &[declared(102, 1.0e-5, 10.0, Some(2))])
        .unwrap()
        .at(0.0);
    let original_total = linear(
        &[1.0e-5, 1.0, 2.0, 3.0, 2.0e7],
        &[11.0, 9.0, 8.0, 11.0, 11.0],
    );
    let total = &tables[&1];
    for e in total.x.iter().copied().chain([0.25, 1.5, 2.5]) {
        let (actual, expected) = (total.eval(e).unwrap(), original_total.eval(e).unwrap());
        assert!(
            (actual - expected).abs() <= 1.0e-12,
            "at {e} eV: {actual} vs {expected}"
        );
    }
    assert_eq!(tables[&2].eval(2.0), Some(8.0));
    assert!(tables[&102].y.iter().all(|&v| v >= 0.0));
    // 範囲の外の弾性散乱は変えない
    assert_eq!(tables[&2].eval(100.0), Some(10.0));
    assert_eq!(tables[&2].y.last(), Some(&10.0));
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
    let RangeKind::Unresolved(urr) = &mf2.isotopes[0].ranges[0].kind else {
        panic!("expected unresolved");
    };
    let xs = UnresolvedXs::new(urr).unwrap();
    let at = |mt: u32, e: f64| tables[&mt].eval(e).unwrap();
    // 格子点（2 keV は補助格子の点）では、平均断面積と MF3 のバックグラウンド断面積の和
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
    // MF3 は非分離共鳴の範囲の中ではバックグラウンドだけで、上端 1e5 eV で範囲の外の値に段差で変わる。
    // 上端の格子点のバックグラウンドは、範囲の内側の値を使う
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
    let RangeKind::Unresolved(urr) = &mf2.isotopes[0].ranges[0].kind else {
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
    let mf2 = one_isotope(vec![ResonanceRange {
        el: 1.0,
        eh: 5.0e3,
        kind: RangeKind::RMatrixLimited(inelastic_level(0.5)),
    }]);
    let mut tables = BTreeMap::from([(2, linear(&[1.0, 5.0e3], &[1.0, 1.0]))]);
    add_resonances(&mut tables, &mf2, 18, 1.0e-3).unwrap();
    // 捕獲のない共鳴でも、共鳴エネルギー（格子点）で止まらない
    let mut capture_free = mf2.clone();
    let RangeKind::RMatrixLimited(rml) = &mut capture_free.isotopes[0].ranges[0].kind else {
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
fn r_matrix_limited_fission_goes_to_the_fission_mt() {
    let mf2 = one_isotope(vec![ResonanceRange {
        el: 1.0,
        eh: 50.0,
        kind: RangeKind::RMatrixLimited(fission_level(1.5e-3, 2.3e-2, 4.0e-2)),
    }]);
    // MF3 が MT19 を持つと、共鳴の核分裂は MT19 に足す（MT18 は部分反応の和で作り直す）
    for fission_mt in [18, 19] {
        let mut tables = BTreeMap::new();
        add_resonances(&mut tables, &mf2, fission_mt, 1.0e-3).unwrap();
        let expected = single_level_exact(6.674, 1.5e-3, 2.3e-2, 4.0e-2).fission;
        let fission = tables[&fission_mt].eval(6.674).unwrap();
        assert!(
            (fission - expected).abs() <= 1e-3 * expected,
            "MT{fission_mt}: {fission} vs {expected}"
        );
        let other = if fission_mt == 18 { 19 } else { 18 };
        assert!(
            !tables.contains_key(&other),
            "MT{other} with MT{fission_mt}"
        );
    }
}

#[test]
fn negative_resolved_elastic_with_its_background_is_raised_to_zero() {
    // R 行列の弾性散乱は負にならないが、MF3 の負のバックグラウンドとの和は負になりうる（ENDF/B-VIII.1 の Fe-54 の 3.5 keV）。
    // 和が負になるところは、和の lin-lin の表を 0 で切った値にする（零点を格子に入れる）。範囲の外と、負でないところはそのまま
    let (el, eh) = (1.0, 5.0e3);
    let mf2 = one_isotope(vec![ResonanceRange {
        el,
        eh,
        kind: RangeKind::RMatrixLimited(inelastic_level(0.5)),
    }]);
    let mut resonance = BTreeMap::new();
    add_resonances(&mut resonance, &mf2, 18, 1.0e-3).unwrap();
    let inside = |e: &f64| el < *e && *e < eh;
    let values = || {
        resonance[&2]
            .x
            .iter()
            .zip(&resonance[&2].y)
            .filter(|(e, _)| inside(e))
    };
    let lowest = values().map(|(_, v)| *v).fold(f64::INFINITY, f64::min);
    let highest = values().map(|(_, v)| *v).fold(0.0, f64::max);
    assert!(lowest >= 0.0 && lowest < highest);
    let background = linear(
        &[1.0e-5, 2.0e7],
        &[-(lowest + 0.01 * (highest - lowest)); 2],
    );
    let unclamped = add(&background, &resonance[&2]);

    let mut tables = BTreeMap::from([(2, background.clone())]);
    add_resonances(&mut tables, &mf2, 18, 1.0e-3).unwrap();
    let elastic = &tables[&2];
    let negative = |t: &Tab1| {
        t.x.iter()
            .zip(&t.y)
            .filter(|(e, v)| inside(e) && **v < 0.0)
            .count()
    };
    assert!(negative(&unclamped) > 0);
    assert_eq!(negative(elastic), 0);
    // 表の点と、点の間の中点のどちらでも、和を 0 で切った値と同じ
    let midpoints = elastic.x.windows(2).map(|w| 0.5 * (w[0] + w[1]));
    for e in elastic.x.iter().copied().chain(midpoints).filter(inside) {
        let expected = unclamped.eval(e).unwrap().max(0.0);
        let actual = elastic.eval(e).unwrap();
        assert!(
            (actual - expected).abs() <= 1.0e-9 * highest,
            "at {e} eV: {actual} vs {expected}"
        );
    }
    // 範囲の外はバックグラウンドのまま
    assert_eq!(elastic.eval(0.5), background.eval(0.5));
}

/// 同位体が 1 つの MF2（AWR は `AWRI`）
fn one_isotope(ranges: Vec<ResonanceRange>) -> Mf2 {
    Mf2 {
        za: 0.0,
        awr: AWRI,
        isotopes: vec![crate::endf::Isotope {
            zai: 0.0,
            abn: 1.0,
            ranges,
        }],
    }
}

/// 分離共鳴の SLBW の範囲 [1, 50] eV だけの MF2
fn slbw_mf2(data: crate::endf::BreitWigner) -> Mf2 {
    one_isotope(vec![ResonanceRange {
        el: 1.0,
        eh: 50.0,
        kind: RangeKind::Slbw(data),
    }])
}

#[test]
fn slbw_range_adds_the_single_level_formula_to_the_background() {
    use crate::resonance::SlbwXs;
    let bw = crate::resonance::slbw_tests::cooperating_levels();
    let mut tables = BTreeMap::from([(2, linear(&[1.0e-5, 2.0e7], &[30.0, 30.0]))]);
    add_resonances(&mut tables, &slbw_mf2(bw.clone()), 18, 1.0e-3).unwrap();
    let xs = SlbwXs::new(&bw);
    for e in [2.0, 9.0, 10.0, 10.5, 30.0] {
        let sigma = xs.cross_sections(e).unwrap();
        let (elastic, capture) = (tables[&2].eval(e).unwrap(), tables[&102].eval(e).unwrap());
        assert!(
            (elastic - (30.0 + sigma.elastic)).abs() <= 2.0e-3 * elastic,
            "elastic at {e}: {elastic} vs 30 + {}",
            sigma.elastic
        );
        assert!(
            (capture - sigma.capture).abs() <= 2.0e-3 * capture,
            "capture at {e}: {capture} vs {}",
            sigma.capture
        );
    }
    // バックグラウンドが打ち消すので、SLBW の負の弾性散乱（9 eV）はそのまま残る
    assert!(xs.cross_sections(9.0).unwrap().elastic < 0.0);
}

#[test]
fn negative_slbw_elastic_with_its_background_is_raised_to_zero() {
    // SLBW の弾性散乱と MF3 のバックグラウンドの和が負になるところは 0 にする（ENDF-102 の 2.4.9.1）。
    // 全断面積は部分反応の和で作り直すので、0 にした分だけ全断面積も変わる
    use crate::resonance::SlbwXs;
    let bw = crate::resonance::slbw_tests::cooperating_levels();
    let xs = SlbwXs::new(&bw);
    let background = 0.5 * -xs.cross_sections(9.0).unwrap().elastic;
    let mut tables = BTreeMap::from([(2, linear(&[1.0e-5, 2.0e7], &[background, background]))]);
    add_resonances(&mut tables, &slbw_mf2(bw), 18, 1.0e-3).unwrap();
    let elastic = &tables[&2];
    assert_eq!(elastic.eval(9.0), Some(0.0));
    assert!(elastic.y.iter().all(|&v| v >= 0.0), "{:?}", elastic.y);
    // 負になる区間の外では、和はそのまま
    let e = 10.0;
    let expected = background + xs.cross_sections(e).unwrap().elastic;
    let actual = elastic.eval(e).unwrap();
    assert!(
        (actual - expected).abs() <= 2.0e-3 * expected,
        "{actual} vs {expected}"
    );
}

#[test]
fn negative_slbw_elastic_is_raised_against_both_sides_of_a_background_step() {
    // MF3 の弾性散乱が段差を持つ点では、和の表は段差の左右の両方の値を持つ。SLBW の弾性散乱が負になる 9 eV で、範囲の上端に
    // 段差がある範囲（[1, 9] eV）と、範囲の中に段差がある範囲（[1, 50] eV）のどちらでも、和は両側で負にならない
    let bw = crate::resonance::slbw_tests::cooperating_levels();
    let step = linear(&[1.0e-5, 9.0, 9.0, 2.0e7], &[0.0, 0.0, 30.0, 30.0]);
    for eh in [9.0, 50.0] {
        let mf2 = one_isotope(vec![ResonanceRange {
            el: 1.0,
            eh,
            kind: RangeKind::Slbw(bw.clone()),
        }]);
        let mut tables = BTreeMap::from([(2, step.clone())]);
        add_resonances(&mut tables, &mf2, 18, 1.0e-3).unwrap();
        let elastic = &tables[&2];
        assert!(
            elastic.y.iter().all(|&v| v >= 0.0),
            "range up to {eh} eV: {:?}",
            elastic
                .x
                .iter()
                .zip(&elastic.y)
                .filter(|(_, v)| **v < 0.0)
                .collect::<Vec<_>>()
        );
    }
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
    // チェーンの反応の部分の MT は、再構成が代表の MT に作り直す和の部分と同じ（xs は代表の MT だけを読む。ADR-0172）
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

/// MLBW の 1 つの範囲 [1e-5, 1e3] eV を持つ同位体の、MF2 の isotope の CONT からの行。値は ZAI、ABN、AP、共鳴エネルギー
fn mlbw_isotope(
    zai: &'static str,
    abn: &'static str,
    ap: &'static str,
    er: &'static str,
) -> Vec<[&'static str; 6]> {
    vec![
        [zai, abn, "0", "0", "1", "0"],
        ["1.0-5", "1.0+3", "1", "2", "0", "0"],
        ["0.0", ap, "0", "0", "1", "0"],
        // l=0: AWRI QX L LRX 6*NRS NRS、値は ER AJ GT GN GG GF
        ["6.338000+1", "0.0", "0", "0", "6", "1"],
        [er, "0.5", "2.3", "2.0", "0.3", "0.0"],
    ]
}

/// 同位体を並べた MF2 と、MF3 の全断面積・弾性散乱・捕獲（どれも一定）の評価
fn evaluation_with_isotopes(isotopes: &[Vec<[&'static str; 6]>]) -> String {
    use crate::endf::rec;
    let nis = isotopes.len().to_string();
    let head = ["3.000000+4", "6.483400+1", "0", "0", nis.as_str(), "0"];
    let full = ["1.0-5", "2.0+7"];
    std::iter::once(head)
        .chain(isotopes.iter().flatten().copied())
        .map(|f| rec(f, 1, 2, 151))
        .chain([
            rec(["", "", "", "", "", ""], 1, 2, 0),
            mf3_section(1, full, ["4.0", "4.0"]),
            mf3_section(2, full, ["3.0", "3.0"]),
            mf3_section(102, full, ["1.0", "1.0"]),
        ])
        .collect()
}

#[test]
fn resonances_of_isotopes_are_summed_with_their_abundances() {
    // ENDF-102 の 2.1: 天然の元素の MF2 は同位体ごとのパラメータを持ち、断面積は同位体の断面積を存在比 ABN（個数の割合）で
    // 重み付けした和。ここでは MF3 のバックグラウンドに 0.6 σ(Zn-64) + 0.4 σ(Zn-66) を足したもの
    use crate::resonance::MlbwXs;
    let a = mlbw_isotope("3.006400+4", "6.0-1", "5.6-1", "1.0+2");
    let b = mlbw_isotope("3.006600+4", "4.0-1", "6.3-1", "3.0+2");
    let alone = |isotope: &[[&'static str; 6]]| {
        let mut lines = isotope.to_vec();
        lines[0][1] = "1.0";
        let tape = Tape::parse(&evaluation_with_isotopes(&[lines])).unwrap();
        let mf2 = crate::endf::Mf2::parse(&tape.section(1, 2, 151).unwrap()).unwrap();
        let RangeKind::Mlbw(bw) = mf2.ranges().next().unwrap().kind.clone() else {
            panic!("expected MLBW");
        };
        MlbwXs::new(&bw)
    };
    let (xs_a, xs_b) = (alone(&a), alone(&b));
    let tape = Tape::parse(&evaluation_with_isotopes(&[a, b])).unwrap();
    let tables = reconstruct(&tape, 1, 1.0e-3).unwrap();
    // 同位体ごとの表は自分の格子点で厳密で、和の表はもう一方の同位体を補間するので、差は線形化の許容誤差の中にある。
    // 表の点と、点の間の中点で確かめる
    for (mt, background, part) in [
        (
            2,
            3.0,
            (|s: &crate::resonance::Sigma| s.elastic) as fn(&_) -> f64,
        ),
        (102, 1.0, |s: &crate::resonance::Sigma| s.capture),
    ] {
        let table = &tables[&mt];
        let midpoints = table.x.windows(2).map(|w| 0.5 * (w[0] + w[1]));
        let inside: Vec<f64> = table
            .x
            .iter()
            .copied()
            .chain(midpoints)
            .filter(|&e| 1.0e-5 < e && e < 1.0e3)
            .collect();
        assert!(inside.len() > 40, "MT{mt}: {inside:?}");
        for e in inside {
            let (sa, sb) = (
                xs_a.cross_sections(e).unwrap(),
                xs_b.cross_sections(e).unwrap(),
            );
            let expected = background + 0.6 * part(&sa) + 0.4 * part(&sb);
            let actual = table.eval(e).unwrap();
            assert!(
                (actual - expected).abs() <= 1.0e-3 * expected,
                "MT{mt} at {e} eV: {actual} vs {expected}"
            );
        }
    }
}
