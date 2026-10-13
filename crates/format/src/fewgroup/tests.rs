use super::*;

/// テストの少数群定数（2群、径方向と z の真空の辺を1つずつと、f* の辺を1つ持つ。z の厚さを持つ）
pub(crate) fn sample() -> FewGroup {
    // 10 進で短く書けない値にして、ビット単位の往復を確かめる
    let x = |v: f64| v / 3.0;
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: x(3.5),
        total: vec![x(1.6), x(2.2)],
        absorption: vec![x(0.03), x(0.3)],
        nu_fission: vec![x(0.02), x(0.4)],
        fission: vec![x(0.008), x(0.16)],
        chi: vec![1.0, 0.0],
        diffusion: vec![x(4.1), x(1.2)],
        scatter: vec![vec![x(1.5), x(0.05)], vec![x(0.002), x(1.9)]],
        adf: vec![vec![x(3.0), x(3.3)]; 6],
        vacuum_ratio: vec![
            None,
            Some(vec![x(1.6), x(1.7)]),
            None,
            None,
            None,
            Some(vec![x(1.4), x(1.1)]),
        ],
        normalized_adf: vec![true, false, false, false, false, false],
        height: Some(x(17.0)),
        pins: None,
    }
}

/// テストのピンの形状（余白の列 1 と、ピンの列 2 の 3×2 のセル。2群）
pub(crate) fn sample_pins() -> Pins {
    let x = |v: f64| v / 3.0;
    Pins {
        dx: vec![0.04, 1.26, 1.26],
        dy: vec![1.26, 1.26],
        form: vec![
            vec![vec![0.0, x(3.1), x(2.9)], vec![0.0, x(3.2), x(2.8)]],
            vec![vec![0.0, x(3.3), x(2.7)], vec![0.0, x(3.0), x(3.0)]],
        ],
        corners: vec![vec![x(2.5), x(3.5)]; 4],
        normalized_corners: vec![true, false, false, false],
        halves: vec![false, true, false, false],
    }
}

#[test]
fn pins_are_checked() {
    // ピン出力の再構成の形状関数と角の不連続因子（ADR-0096）。列の幅と行の高さは正で、形状関数は
    // [群][行][列] でセルの数だけあり、負でない。角の不連続因子は 4 つの角に群の数だけあり、正
    let with = |pins: Pins| FewGroup {
        pins: Some(pins),
        ..sample()
    };
    assert!(with(sample_pins()).validate().is_ok());
    let bad = |change: fn(&mut Pins), word: &str| {
        let mut pins = sample_pins();
        change(&mut pins);
        let err = with(pins).validate().unwrap_err().to_string();
        assert!(err.contains(word), "{err}");
    };
    bad(|p| p.dx[0] = 0.0, "pins");
    bad(|p| p.dx[1] = -1.26, "pins");
    bad(|p| p.dy.clear(), "pins");
    bad(|p| p.dy[1] = f64::INFINITY, "pins");
    bad(|p| p.form.pop().map(drop).unwrap_or(()), "form");
    bad(|p| p.form[0].pop().map(drop).unwrap_or(()), "form");
    bad(|p| p.form[1][0].pop().map(drop).unwrap_or(()), "form");
    bad(|p| p.form[0][1][2] = -0.1, "form");
    bad(|p| p.form[0][1][2] = f64::NAN, "form");
    bad(|p| p.corners.pop().map(drop).unwrap_or(()), "corners");
    bad(|p| p.corners[2] = vec![1.0], "corners");
    bad(|p| p.corners[3][1] = 0.0, "corners");
    // 辺ごとの、隣の部分ノードと分け合ったセルの印は 4 辺ぶん（ADR-0108）
    bad(|p| p.halves.truncate(3), "halves");
    // 燃料の CDF で割った角の印は 4 つの角ぶん（ADR-0131）
    bad(|p| p.normalized_corners.truncate(3), "normalized_corners");
}

#[test]
fn validate_checks_values_in_memory() {
    // ファイルを経ずに作った値（core の Problem など）も同じ検査を通す
    assert!(sample().validate().is_ok());
    let mut fg = sample();
    fg.adf.pop();
    assert!(fg.validate().is_err());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut fg = sample();
        fg.adf[1][0] = bad;
        let err = fg.validate().unwrap_err().to_string();
        assert!(err.contains("discontinuity"), "{err}");
    }
    let mut fg = sample();
    fg.diffusion.pop();
    assert!(fg.validate().is_err());
    let mut fg = sample();
    fg.k_inf = 0.0;
    assert!(fg.validate().is_err());
    // 群境界が降順でない、散乱の行が短い
    let mut fg = sample();
    fg.bounds.reverse();
    assert!(fg.validate().is_err());
    let mut fg = sample();
    fg.scatter[0].pop();
    assert!(fg.validate().is_err());
}

#[test]
fn vacuum_ratios_are_checked() {
    // 真空の辺の J/φ は、出ていく角度中性子束の μ の平均なので 0 < r ≤ 1（ADR-0101）
    let mut fg = sample();
    fg.vacuum_ratio[1] = Some(vec![1.0, 0.5]);
    assert!(fg.validate().is_ok());
    for bad in [0.0, -0.1, 1.0 + 1e-12, f64::NAN, f64::INFINITY] {
        let mut fg = sample();
        fg.vacuum_ratio[1] = Some(vec![bad, 0.5]);
        let err = fg.validate().unwrap_err().to_string();
        assert!(err.contains("vacuum_ratio"), "{err}");
    }
    let mut fg = sample();
    fg.vacuum_ratio.pop();
    assert!(fg.validate().is_err());
    let mut fg = sample();
    fg.vacuum_ratio[0] = Some(vec![0.5]);
    assert!(fg.validate().is_err());
    // 真空の辺がなければ、どの辺も持たない
    let mut fg = sample();
    fg.vacuum_ratio = vec![None; 6];
    assert!(fg.validate().is_ok());
}

#[test]
fn sides_include_the_z_faces() {
    // 辺は（左、右、上、下、−z、+z）の 6 つ（ADR-0128）。ADF、真空の比、f* の印のどれも 6 辺ぶん
    let bad = |change: fn(&mut FewGroup), word: &str| {
        let mut fg = sample();
        change(&mut fg);
        let err = fg.validate().unwrap_err().to_string();
        assert!(err.contains(word), "{err}");
    };
    bad(|f| f.adf.truncate(4), "adf");
    bad(|f| f.vacuum_ratio.truncate(4), "vacuum_ratio");
    bad(|f| f.normalized_adf.truncate(4), "normalized_adf");
    bad(|f| f.adf[5][1] = 0.0, "discontinuity");
    bad(|f| f.vacuum_ratio[5] = Some(vec![1.5, 0.5]), "vacuum_ratio");
    // z の真空の辺にも f* の印は付けない
    bad(|f| f.normalized_adf[5] = true, "normalized_adf");
    let mut fg = sample();
    fg.normalized_adf[4] = true;
    assert!(fg.validate().is_ok());
}

#[test]
fn height_is_positive_when_given() {
    // z の厚さは軸方向の定数だけが持ち（ADR-0128）、持つなら正
    let mut fg = sample();
    fg.height = None;
    assert!(fg.validate().is_ok());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut fg = sample();
        fg.height = Some(bad);
        let err = fg.validate().unwrap_err().to_string();
        assert!(err.contains("height"), "{err}");
    }
}

#[test]
fn normalized_sides_are_checked() {
    // 反射体の定数は、燃料と接する辺の ADF を燃料の ADF で割った f* で持ち、その辺に印を付ける
    // （ADR-0095）。印は4辺ぶんで、外周の真空の辺（燃料と接しない）には付けない
    assert!(sample().validate().is_ok());
    let mut fg = sample();
    fg.normalized_adf.pop();
    let err = fg.validate().unwrap_err().to_string();
    assert!(err.contains("normalized_adf"), "{err}");
    let mut fg = sample();
    fg.normalized_adf[1] = true;
    let err = fg.validate().unwrap_err().to_string();
    assert!(err.contains("normalized_adf"), "{err}");
    let mut fg = sample();
    fg.normalized_adf.push(false);
    assert!(fg.validate().is_err());
}
