use super::*;

fn state(k: f64) -> FewGroup {
    FewGroup {
        k_inf: k / 3.0,
        ..crate::fewgroup::tests::sample()
    }
}

/// ホウ素 2 点 × 減速材 2 点 × 燃料 1 点
fn sample() -> Table {
    let axes = vec![
        Axis::boron(vec![900.0, 1700.0], 0.199),
        Axis::moderator(vec![550.0, 580.0], vec![0.76972, 0.71190]),
        Axis::fuel(vec![565.0]),
    ];
    Table {
        points: product_points(&axes),
        axes,
        split: 1,
        states: (0..4).map(|i| state(3.0 + f64::from(i))).collect(),
        heavy_metal: None,
        micro: None,
        histories: Vec::new(),
    }
}

#[test]
fn round_trips_bit_exactly() {
    let table = sample();
    assert_eq!(
        Table::from_msgpack(&table.to_msgpack().unwrap()).unwrap(),
        table
    );
    let single = Table::single(state(3.5));
    assert!(single.axes.is_empty());
    assert_eq!(
        Table::from_msgpack(&single.to_msgpack().unwrap()).unwrap(),
        single
    );
    // z の厚さのない（2次元の格子計算の）定数
    let flat = Table::single(FewGroup {
        height: None,
        ..state(3.5)
    });
    assert_eq!(
        Table::from_msgpack(&flat.to_msgpack().unwrap()).unwrap(),
        flat
    );
}

#[test]
fn invalid_tables_are_errors() {
    let bad = |change: fn(&mut Table)| {
        let mut table = sample();
        change(&mut table);
        // 書くときも読むときも確かめる
        assert!(table.validate().is_err());
        assert!(table.to_msgpack().is_err());
    };
    // 状態の数が軸の点の数の積と違う
    bad(|t| {
        t.states.pop();
    });
    // 点が増えない、空、有限でない
    bad(|t| t.axes[0] = Axis::boron(vec![1700.0, 900.0], 0.199));
    bad(|t| t.axes[2] = Axis::fuel(vec![]));
    bad(|t| t.axes[2] = Axis::fuel(vec![f64::NAN]));
    bad(|t| t.axes[2] = Axis::fuel(vec![-1.0]));
    bad(|t| t.axes[0] = Axis::boron(vec![-1.0, 900.0], 0.199));
    // 密度の数が点の数と違う、密度が正でない、B-10 の割合が範囲の外
    bad(|t| t.axes[1] = Axis::moderator(vec![550.0, 580.0], vec![0.76972]));
    bad(|t| t.axes[1] = Axis::moderator(vec![550.0, 580.0], vec![0.76972, 0.0]));
    bad(|t| t.axes[0] = Axis::boron(vec![900.0, 1700.0], 1.5));
    // 同じ軸が2つ
    bad(|t| t.axes[2] = Axis::boron(vec![1300.0], 0.199));
    // 群境界、f* の印、真空の辺が状態の間で違う
    bad(|t| t.states[3].bounds = vec![2.0e7, 0.625, 1.0e-5]);
    bad(|t| t.states[3].normalized_adf = vec![false; 6]);
    bad(|t| t.states[3].vacuum_ratio[1] = None);
    bad(|t| t.states[3].vacuum_ratio[5] = None);
    // z の厚さも状態によらない（ADR-0128）
    bad(|t| t.states[2].height = None);
    bad(|t| t.states[2].height = Some(12.0));
    // 状態そのものの誤り
    bad(|t| t.states[1].total[0] = f64::NAN);
}

#[test]
fn reading_checks_the_format_name() {
    let mut bytes = Table::single(state(3.5)).to_msgpack().unwrap();
    let at = bytes.windows(9).position(|w| w == b"xlc-table").unwrap();
    bytes[at + 4] = b'T';
    assert!(Table::from_msgpack(&bytes).is_err());
    assert!(Table::from_msgpack(&[0x90]).is_err());
    // 末尾の余り
    let mut bytes = Table::single(state(3.5)).to_msgpack().unwrap();
    bytes.push(0xc0);
    assert!(Table::from_msgpack(&bytes).is_err());
}

#[test]
fn reading_rejects_wrong_axis_fields() {
    // 軸の欄を書き換えたファイルを読む（ADR-0013 の厳密な読み込み）
    let edit = |change: fn(&mut rmpv::Value)| {
        let bytes = sample().to_msgpack().unwrap();
        let mut value = msgpack::read_value(&bytes).unwrap();
        change(&mut value);
        let mut out = Vec::new();
        rmpv::encode::write_value(&mut out, &value).unwrap();
        Table::from_msgpack(&out)
    };
    // axes の n 番目の軸の map
    fn axis(value: &mut rmpv::Value, n: usize) -> &mut Vec<(rmpv::Value, rmpv::Value)> {
        let rmpv::Value::Map(entries) = value else {
            panic!()
        };
        let (_, axes) = entries
            .iter_mut()
            .find(|(k, _)| k.as_str() == Some("axes"))
            .unwrap();
        let rmpv::Value::Array(axes) = axes else {
            panic!()
        };
        let rmpv::Value::Map(fields) = &mut axes[n] else {
            panic!()
        };
        fields
    }
    assert!(edit(|_| {}).is_ok());
    // 知らない名前
    assert!(edit(|v| axis(v, 2)[0].1 = "clad".into()).is_err());
    // ホウ素に b10 がない、減速材に density がない
    assert!(edit(|v| axis(v, 0).retain(|(k, _)| k.as_str() != Some("b10"))).is_err());
    assert!(edit(|v| axis(v, 1).retain(|(k, _)| k.as_str() != Some("density"))).is_err());
    // 燃料に density や b10 がある
    assert!(edit(|v| axis(v, 2).push(("b10".into(), 0.2.into()))).is_err());
    assert!(
        edit(|v| axis(v, 2).push(("density".into(), vec![rmpv::Value::from(0.7)].into()))).is_err()
    );
}

#[test]
fn pins_round_trip_and_share_cells() {
    // 形状関数と角の不連続因子も状態ごとに持ち、セルの境界はすべての状態で同じ（ADR-0096）
    let with_pins = |mut table: Table| {
        for (i, s) in table.states.iter_mut().enumerate() {
            let mut pins = crate::fewgroup::tests::sample_pins();
            pins.form[0][1][1] += f64::from(u8::try_from(i).unwrap()) / 7.0;
            s.pins = Some(pins);
        }
        table
    };
    let table = with_pins(sample());
    assert_eq!(
        Table::from_msgpack(&table.to_msgpack().unwrap()).unwrap(),
        table
    );
    let mut table = with_pins(sample());
    if let Some(p) = table.states[2].pins.as_mut() {
        p.dx[1] = 1.27;
    }
    let err = table.validate().unwrap_err().to_string();
    assert!(err.contains("pins"), "{err}");
    // 分け合ったセルの印も状態の間で同じ
    let mut table = with_pins(sample());
    if let Some(p) = table.states[1].pins.as_mut() {
        p.halves[0] = !p.halves[0];
    }
    let err = table.validate().unwrap_err().to_string();
    assert!(err.contains("pins"), "{err}");
    // 燃料の CDF で割った角の印も状態の間で同じ（ADR-0131）
    let mut table = with_pins(sample());
    if let Some(p) = table.states[1].pins.as_mut() {
        p.normalized_corners[2] = !p.normalized_corners[2];
    }
    let err = table.validate().unwrap_err().to_string();
    assert!(err.contains("pins"), "{err}");
    let mut table = with_pins(sample());
    table.states[3].pins = None;
    let err = table.validate().unwrap_err().to_string();
    assert!(err.contains("pins"), "{err}");
}

#[test]
fn states_need_the_pins_key() {
    // 形状関数のない状態も `pins` に nil を書く。キーのないファイル（古い形式）は読まない
    let bytes = sample().to_msgpack().unwrap();
    let mut value = msgpack::read_value(&bytes).unwrap();
    let rmpv::Value::Map(entries) = &mut value else {
        panic!()
    };
    let (_, states) = entries
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("states"))
        .unwrap();
    let rmpv::Value::Array(states) = states else {
        panic!()
    };
    let rmpv::Value::Map(fields) = &mut states[0] else {
        panic!()
    };
    assert!(
        fields
            .iter()
            .any(|(k, v)| k.as_str() == Some("pins") && v.is_nil())
    );
    fields.retain(|(k, _)| k.as_str() != Some("pins"));
    let mut out = Vec::new();
    rmpv::encode::write_value(&mut out, &value).unwrap();
    let err = Table::from_msgpack(&out).unwrap_err().to_string();
    assert!(err.contains("pins") || err.contains("missing"), "{err}");
}

#[test]
fn split_tables_hold_sub_nodes_per_state() {
    // 集合体を 2 × 2 に分けたテーブル（ADR-0108）は、状態ごとに部分ノードの定数を 4 つ持つ（部分ノードがいちばん
    // 速く変わる）
    let mut table = sample();
    table.split = 2;
    table.states = (0..16).map(|i| state(3.0 + f64::from(i))).collect();
    assert_eq!(
        Table::from_msgpack(&table.to_msgpack().unwrap()).unwrap(),
        table
    );
    assert_eq!(table.sub_nodes(), 4);
    // 状態の数 × 部分ノードの数でないとエラー
    let mut short = table.clone();
    short.states.pop();
    assert!(short.validate().is_err());
    let mut zero = table.clone();
    zero.split = 0;
    assert!(zero.validate().is_err());
    // 状態の間で比べる量（真空の辺など）は、同じ部分ノードどうしで比べる
    // （部分ノード 1 だけ印を変えるのはよいが、ある状態の部分ノード 1 だけ違えばエラー）
    let mut mixed = table.clone();
    for state in 0..4 {
        mixed.states[4 * state + 1].normalized_adf[0] = false;
    }
    assert!(mixed.validate().is_ok());
    mixed.states[5].normalized_adf[0] = true;
    let err = mixed.validate().unwrap_err().to_string();
    assert!(err.contains("state"), "{err}"); // 群境界は表全体の量なので、部分ノードで違えばエラー（ファイルには1つしか書かない）
    let mut bounds = table.clone();
    for state in 0..4 {
        bounds.states[4 * state + 1].bounds = vec![2.0e7, 0.625, 1.0e-5];
    }
    let err = bounds.validate().unwrap_err().to_string();
    assert!(err.contains("bounds"), "{err}");
}

#[test]
fn states_need_the_height_key() {
    // z の厚さのない定数も `height` に nil を書く（ADR-0128、ADR-0129）。キーのないファイルは読まない
    let table = Table::single(FewGroup {
        height: None,
        ..state(3.5)
    });
    let mut value = msgpack::read_value(&table.to_msgpack().unwrap()).unwrap();
    let rmpv::Value::Map(entries) = &mut value else {
        panic!()
    };
    let (_, states) = entries
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("states"))
        .unwrap();
    let rmpv::Value::Array(states) = states else {
        panic!()
    };
    let rmpv::Value::Map(fields) = &mut states[0] else {
        panic!()
    };
    assert!(
        fields
            .iter()
            .any(|(k, v)| k.as_str() == Some("height") && v.is_nil())
    );
    fields.retain(|(k, _)| k.as_str() != Some("height"));
    let mut out = Vec::new();
    rmpv::encode::write_value(&mut out, &value).unwrap();
    let err = Table::from_msgpack(&out).unwrap_err().to_string();
    assert!(err.contains("height") || err.contains("missing"), "{err}");
}

/// sample の 3 つの連続の軸に、制御棒の離散の軸（なし、AIC、B4C）を足したもの（ADR-0097、ADR-0127）
fn rodded() -> Table {
    let mut table = sample();
    table.axes.push(Axis::Rod {
        values: vec!["out".into(), "aic".into(), "b4c".into()],
    });
    table.points = product_points(&table.axes);
    table.states = (0..12).map(|i| state(3.0 + f64::from(i))).collect();
    table
}

#[test]
fn rod_axes_round_trip_with_named_values() {
    let table = rodded();
    assert_eq!(table.axes[3].name(), "rod");
    assert_eq!(table.axes[3].len(), 3);
    assert_eq!(
        Table::from_msgpack(&table.to_msgpack().unwrap()).unwrap(),
        table
    );
}

#[test]
fn rod_axes_need_distinct_nonempty_values() {
    let bad = |values: Vec<&str>| {
        let mut table = rodded();
        table.axes[3] = Axis::Rod {
            values: values.into_iter().map(String::from).collect(),
        };
        table.points = product_points(&table.axes);
        table.states.truncate(4 * table.axes[3].len().max(1));
        assert!(table.validate().is_err());
    };
    bad(vec![]);
    bad(vec!["out", "out"]);
    bad(vec!["out", ""]);
    bad(vec!["out", "a b"]);
}

/// ファイルの値の、axes の n 番目の軸の map
fn axis_fields(value: &mut rmpv::Value, n: usize) -> &mut Vec<(rmpv::Value, rmpv::Value)> {
    let rmpv::Value::Map(entries) = value else {
        panic!()
    };
    let (_, axes) = entries
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("axes"))
        .unwrap();
    let rmpv::Value::Array(axes) = axes else {
        panic!()
    };
    let rmpv::Value::Map(fields) = &mut axes[n] else {
        panic!()
    };
    fields
}

#[test]
fn reading_rejects_points_on_rod_axes_and_values_on_continuous_axes() {
    let edit = |change: fn(&mut rmpv::Value)| {
        let bytes = rodded().to_msgpack().unwrap();
        let mut value = msgpack::read_value(&bytes).unwrap();
        change(&mut value);
        let mut out = Vec::new();
        rmpv::encode::write_value(&mut out, &value).unwrap();
        Table::from_msgpack(&out)
    };
    assert!(edit(|_| {}).is_ok());
    // 制御棒の軸に点がある、値がない
    assert!(
        edit(|v| axis_fields(v, 3).push(("points".into(), vec![rmpv::Value::from(1.0)].into())))
            .is_err()
    );
    assert!(edit(|v| axis_fields(v, 3).retain(|(k, _)| k.as_str() != Some("values"))).is_err());
    // 連続の軸に値の名前がある、点がない
    assert!(
        edit(|v| axis_fields(v, 2).push(("values".into(), vec![rmpv::Value::from("out")].into())))
            .is_err()
    );
    assert!(edit(|v| axis_fields(v, 2).retain(|(k, _)| k.as_str() != Some("points"))).is_err());
}

/// 燃焼度 2 点 × ホウ素 2 点 × 減速材 2 点 × 燃料 1 点（燃焼度がいちばん遅く変わる。ADR-0187）
fn burnup_sample() -> Table {
    let mut table = sample();
    table.axes.insert(
        0,
        Axis::Burnup {
            points: vec![0.0, 10.0],
        },
    );
    table.points = product_points(&table.axes);
    table.states = (0..8).map(|i| state(3.0 + f64::from(i))).collect();
    table
}

#[test]
fn burnup_axes_round_trip_bit_exactly() {
    let table = burnup_sample();
    assert_eq!(table.axes[0].name(), "burnup");
    assert_eq!(table.axes[0].points(), &[0.0, 10.0]);
    assert_eq!(
        Table::from_msgpack(&table.to_msgpack().unwrap()).unwrap(),
        table
    );
}

#[test]
fn burnup_points_must_be_nonnegative_and_increasing() {
    let bad = |points: Vec<f64>| {
        let mut table = burnup_sample();
        table.axes[0] = Axis::Burnup { points };
        assert!(table.validate().is_err());
    };
    bad(vec![-1.0, 10.0]);
    bad(vec![10.0, 0.0]);
    bad(vec![0.0, f64::NAN]);
    bad(vec![0.0, f64::INFINITY]);
    // 状態の数は燃焼度の点も数える
    let mut table = burnup_sample();
    table.states.truncate(4);
    assert!(table.validate().is_err());
    let mut table = burnup_sample();
    table.states.truncate(4);
    table.points.truncate(4);
    assert!(table.validate().is_err());
}

#[test]
fn heavy_metal_round_trips_and_is_written_only_when_held() {
    let mut table = burnup_sample();
    table.heavy_metal = Some(4321.0 / 7.0);
    let bytes = table.to_msgpack().unwrap();
    assert_eq!(Table::from_msgpack(&bytes).unwrap(), table);
    // 持たないテーブルはキーを書かない（ADR-0190）
    let fresh = sample();
    assert_eq!(fresh.heavy_metal, None);
    let value = msgpack::read_value(&fresh.to_msgpack().unwrap()).unwrap();
    let keys: Vec<&str> = value
        .as_map()
        .unwrap()
        .iter()
        .filter_map(|(k, _)| k.as_str())
        .collect();
    assert!(!keys.contains(&"heavy_metal"), "{keys:?}");
    assert_eq!(
        Table::from_msgpack(&fresh.to_msgpack().unwrap()).unwrap(),
        fresh
    );
}

#[test]
fn heavy_metal_must_be_positive() {
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut table = burnup_sample();
        table.heavy_metal = Some(bad);
        assert!(table.validate().is_err(), "{bad}");
    }
}

#[test]
fn points_round_trip_and_index_the_states() {
    let table = rodded();
    assert_eq!(table.points.len(), 12);
    assert_eq!(table.points[5], vec![0, 1, 0, 2]);
    assert_eq!(table.point_index(&[0, 1, 0, 2]), Some(5));
    assert_eq!(table.point_index(&[0, 2, 0, 0]), None);
    assert_eq!(
        Table::from_msgpack(&table.to_msgpack().unwrap()).unwrap(),
        table
    );
    // 軸のないテーブルの点は、番号のない 1 点
    assert_eq!(Table::single(state(3.5)).points, vec![Vec::<usize>::new()]);
}

#[test]
fn points_must_index_the_axes_in_increasing_order() {
    let bad = |change: fn(&mut Table)| {
        let mut table = sample();
        change(&mut table);
        assert!(table.validate().is_err());
        assert!(table.to_msgpack().is_err());
    };
    // 軸の範囲の外、番号の数が軸の数と違う、同じ点が 2 つ、辞書式に増えない
    bad(|t| t.points[3] = vec![1, 2, 0]);
    bad(|t| t.points[3] = vec![1, 1]);
    bad(|t| t.points[3] = vec![1, 0, 0]);
    bad(|t| t.points.swap(1, 2));
    // 状態の数が点の数 × 部分ノードの数と違う
    bad(|t| {
        t.points.pop();
    });
}

#[test]
fn reading_needs_the_points() {
    let bytes = sample().to_msgpack().unwrap();
    let mut value = msgpack::read_value(&bytes).unwrap();
    let rmpv::Value::Map(entries) = &mut value else {
        panic!()
    };
    entries.retain(|(k, _)| k.as_str() != Some("points"));
    let mut out = Vec::new();
    rmpv::encode::write_value(&mut out, &value).unwrap();
    let err = Table::from_msgpack(&out).unwrap_err().to_string();
    assert!(err.contains("points") || err.contains("missing"), "{err}");
}

#[test]
fn state_axes_carry_their_interpolation_variable_and_nominal_value() {
    // ADR-0097 の内挿の変数: ホウ素は ppm、減速材は温度、燃料は √T
    let table = sample();
    let transforms: Vec<_> = table
        .axes
        .iter()
        .map(|a| a.state_axis().unwrap().1)
        .collect();
    assert_eq!(
        transforms,
        [Transform::Linear, Transform::Linear, Transform::Sqrt]
    );
    assert_eq!(Transform::Sqrt.apply(16.0), 4.0);
    assert_eq!(Axis::Burnup { points: vec![0.0] }.state_axis(), None);
    let mut nominal = sample();
    nominal.axes = vec![
        Axis::boron(vec![900.0, 1700.0], 0.199).with_nominal(1700.0),
        Axis::moderator(vec![550.0, 580.0], vec![0.76972, 0.71190]).with_nominal(550.0),
        Axis::fuel(vec![565.0]).with_nominal(565.0),
    ];
    assert_eq!(nominal.axes[0].nominal_index(), Some(1));
    assert_eq!(nominal.axes[1].nominal_index(), Some(0));
    assert_eq!(
        Table::from_msgpack(&nominal.to_msgpack().unwrap()).unwrap(),
        nominal
    );
}

#[test]
fn reading_checks_the_transform_and_the_nominal_fields() {
    let edit = |change: fn(&mut rmpv::Value)| {
        let mut table = rodded();
        table.axes.insert(0, Axis::Burnup { points: vec![0.0] });
        table.points = product_points(&table.axes);
        let bytes = table.to_msgpack().unwrap();
        let mut value = msgpack::read_value(&bytes).unwrap();
        change(&mut value);
        let mut out = Vec::new();
        rmpv::encode::write_value(&mut out, &value).unwrap();
        Table::from_msgpack(&out)
    };
    assert!(edit(|_| {}).is_ok());
    // ホウ素、減速材、燃料は内挿の変数を書く。知らない変数はエラー
    assert!(edit(|v| axis_fields(v, 3).retain(|(k, _)| k.as_str() != Some("transform"))).is_err());
    assert!(
        edit(|v| {
            for (k, x) in axis_fields(v, 1).iter_mut() {
                if k.as_str() == Some("transform") {
                    *x = "log".into();
                }
            }
        })
        .is_err()
    );
    // 燃焼度と制御棒の軸は内挿の変数も公称の値も持たない
    assert!(edit(|v| axis_fields(v, 0).push(("transform".into(), "linear".into()))).is_err());
    assert!(edit(|v| axis_fields(v, 4).push(("nominal".into(), 0.0.into()))).is_err());
    assert!(edit(|v| axis_fields(v, 0).push(("nominal".into(), 0.0.into()))).is_err());
}

#[test]
fn nominal_values_are_points_given_to_all_state_axes_or_none() {
    let bad = |axes: Vec<Axis>| {
        let mut table = sample();
        table.axes = axes;
        assert!(table.validate().is_err());
    };
    // 格子点でない
    bad(vec![
        Axis::boron(vec![900.0, 1700.0], 0.199).with_nominal(1300.0),
        Axis::moderator(vec![550.0, 580.0], vec![0.76972, 0.71190]).with_nominal(550.0),
        Axis::fuel(vec![565.0]).with_nominal(565.0),
    ]);
    // 一部の軸だけ
    bad(vec![
        Axis::boron(vec![900.0, 1700.0], 0.199).with_nominal(900.0),
        Axis::moderator(vec![550.0, 580.0], vec![0.76972, 0.71190]),
        Axis::fuel(vec![565.0]),
    ]);
}

/// ホウ素、減速材、燃料が 3 点ずつで公称が中央の軸と、2 軸の組すべての部分格子の点（cut-HDMR の 2 次。ADR-0200）
fn pairs() -> Table {
    let axes = vec![
        Axis::boron(vec![900.0, 1300.0, 1700.0], 0.199).with_nominal(1300.0),
        Axis::moderator(vec![550.0, 565.0, 580.0], vec![0.76972, 0.74271, 0.71190])
            .with_nominal(565.0),
        Axis::fuel(vec![550.0, 565.0, 580.0]).with_nominal(565.0),
    ];
    // 3 つの軸がどれも公称でない点（3 次の成分）だけを除く
    let points: Vec<Vec<usize>> = product_points(&axes)
        .into_iter()
        .filter(|p| p.contains(&1))
        .collect();
    Table {
        states: (0..points.len()).map(|i| state(3.0 + i as f64)).collect(),
        points,
        axes,
        split: 1,
        heavy_metal: None,
        micro: None,
        histories: Vec::new(),
    }
}

#[test]
fn cut_hdmr_points_hold_full_subgrids_closed_under_subsets() {
    // 1 + Σ(nᵢ − 1) + Σᵢ<ⱼ(nᵢ − 1)(nⱼ − 1) = 1 + 6 + 12 状態（ADR-0200）
    let table = pairs();
    assert_eq!(table.points.len(), 19);
    assert_eq!(
        Table::from_msgpack(&table.to_msgpack().unwrap()).unwrap(),
        table
    );
    assert_eq!(component(&[0, 1, 2], &[1, 1, 1]), 0b101);
    let without = |keep: fn(&[usize]) -> bool| {
        let mut t = pairs();
        t.points.retain(|p| keep(p));
        t.states.truncate(t.points.len());
        t
    };
    // 部分格子の点が欠ける
    let err = without(|p| p != [0, 0, 1])
        .validate()
        .unwrap_err()
        .to_string();
    assert!(err.contains("subgrids"), "{err}");
    // ホウ素と減速材の部分格子を持ち、ホウ素だけの分岐を持たない（部分集合について閉じない）
    let err = without(|p| !(p[1] == 1 && p[2] == 1 && p[0] != 1))
        .validate()
        .unwrap_err()
        .to_string();
    assert!(err.contains("subsets"), "{err}");
    // 公称の値のないテーブルは直積
    let mut product = pairs();
    for axis in &mut product.axes {
        if let Axis::Boron { nominal, .. }
        | Axis::Moderator { nominal, .. }
        | Axis::Fuel { nominal, .. } = axis
        {
            *nominal = None;
        }
    }
    let err = product.validate().unwrap_err().to_string();
    assert!(err.contains("product"), "{err}");
    product.points = product_points(&product.axes);
    product.states = (0..27).map(|i| state(3.0 + f64::from(i))).collect();
    assert!(product.validate().is_ok());
    // 公称の値を持つ直積（3 軸の成分を含む、すべての成分の部分格子）
    let mut full = pairs();
    full.points = product_points(&full.axes);
    full.states = (0..27).map(|i| state(3.0 + f64::from(i))).collect();
    assert!(full.validate().is_ok());
}

#[test]
fn rod_values_hold_their_own_subgrids() {
    // 制御棒を入れた値は、公称の状態を制御棒を入れた状態にした別の組（ADR-0200）。組ごとに成分が違ってよい
    let mut table = pairs();
    table.axes.push(Axis::Rod {
        values: vec!["out".into(), "aic".into()],
    });
    let out: Vec<Vec<usize>> = pairs().points;
    // 制御棒を入れた組は 1 軸ずつの分岐だけ
    let rodded: Vec<Vec<usize>> = out
        .iter()
        .filter(|p| p.iter().filter(|&&i| i != 1).count() <= 1)
        .cloned()
        .collect();
    assert_eq!(rodded.len(), 7);
    let mut points: Vec<Vec<usize>> = out
        .iter()
        .map(|p| [p.as_slice(), &[0]].concat())
        .chain(rodded.iter().map(|p| [p.as_slice(), &[1]].concat()))
        .collect();
    points.sort();
    table.states = (0..points.len()).map(|i| state(3.0 + i as f64)).collect();
    table.points = points;
    assert!(table.validate().is_ok());
    // 値の点がない
    let mut empty = table.clone();
    empty.points.retain(|p| p[3] == 0);
    empty.states.truncate(empty.points.len());
    let err = empty.validate().unwrap_err().to_string();
    assert!(err.contains("aic"), "{err}");
}

/// pairs の点を燃焼度 0、10、20 MWd/kg で、公称の状態はすべての点で、ほかの点は 0 と 20 だけで解いたテーブル（ADR-0200 の
/// 案 B2）
fn thinned() -> Table {
    let mut table = pairs();
    table.axes.insert(
        0,
        Axis::Burnup {
            points: vec![0.0, 10.0, 20.0],
        },
    );
    table.points = (0..3)
        .flat_map(|b| {
            pairs()
                .points
                .into_iter()
                .filter(move |p| b != 1 || p == &[1, 1, 1])
                .map(move |p| [&[b][..], p.as_slice()].concat())
        })
        .collect();
    table.states = (0..table.points.len())
        .map(|i| state(3.0 + i as f64))
        .collect();
    table
}

#[test]
fn points_solved_at_some_burnups_need_the_nominal_state_at_every_burnup() {
    let table = thinned();
    assert_eq!(table.points.len(), 19 + 1 + 19);
    assert_eq!(
        Table::from_msgpack(&table.to_msgpack().unwrap()).unwrap(),
        table
    );
    // 公称の状態が燃焼度 10 を持たない
    let mut missing = thinned();
    missing.points.retain(|p| p != &[1, 1, 1, 1]);
    missing.states.pop();
    let err = missing.validate().unwrap_err().to_string();
    assert!(err.contains("every burnup"), "{err}");
    // 公称の値がなければ、燃焼度の点の差を作れない
    let mut product = thinned();
    product.axes = vec![
        product.axes[0].clone(),
        Axis::boron(vec![900.0, 1300.0, 1700.0], 0.199),
        Axis::moderator(vec![550.0, 565.0, 580.0], vec![0.76972, 0.74271, 0.71190]),
        Axis::fuel(vec![550.0, 565.0, 580.0]),
    ];
    product.points = (0..3)
        .flat_map(|b| {
            product_points(&product.axes[1..])
                .into_iter()
                .filter(move |p| b != 1 || p == &[1, 1, 1])
                .map(move |p| [&[b][..], p.as_slice()].concat())
        })
        .collect();
    product.states = (0..product.points.len())
        .map(|i| state(3.0 + i as f64))
        .collect();
    let err = product.validate().unwrap_err().to_string();
    assert!(err.contains("nominal"), "{err}");
    // 制御棒の値があれば、公称の状態は out の組のもの
    let mut rodded = thinned();
    rodded.axes.push(Axis::Rod {
        values: vec!["aic".into(), "out".into()],
    });
    let mut points: Vec<Vec<usize>> = rodded
        .points
        .iter()
        .flat_map(|p| {
            let in_ = [p.as_slice(), &[0]].concat();
            let out = [p.as_slice(), &[1]].concat();
            // 制御棒を入れた組は燃焼度 10 で解かない
            if p[0] == 1 { vec![out] } else { vec![in_, out] }
        })
        .collect();
    points.sort();
    rodded.states = (0..points.len()).map(|i| state(3.0 + i as f64)).collect();
    rodded.points = points;
    assert!(rodded.validate().is_ok());
    let mut no_out = rodded.clone();
    no_out.axes[4] = Axis::Rod {
        values: vec!["aic".into(), "b4c".into()],
    };
    let err = no_out.validate().unwrap_err().to_string();
    assert!(err.contains("out"), "{err}");
}

#[test]
fn points_solved_at_some_burnups_span_the_burnup_axis() {
    // core は燃焼度の差を外挿しないので、どの点も燃焼度の軸の最初と最後の点で解く（ADR-0200）
    let mut short = thinned();
    let last = short.points.len();
    short.points.retain(|p| !(p[0] == 2 && p[1..] != [1, 1, 1]));
    short.states.truncate(short.points.len());
    assert!(short.points.len() < last);
    let err = short.validate().unwrap_err().to_string();
    assert!(err.contains("first and the last"), "{err}");
}

/// 核種 `scale` 倍の値の微視的断面積（2 群。核分裂するなら収率の境界より下の割合を持つ）
fn micro_nuclide(scale: f64, fissile: bool) -> MicroNuclide {
    let x = |v: f64| scale * v / 3.0;
    MicroNuclide {
        total: vec![x(11.0), x(13.0)],
        nu_fission: if fissile {
            vec![x(1.1), x(70.0)]
        } else {
            vec![0.0; 2]
        },
        fission: if fissile {
            vec![x(0.4), x(29.0)]
        } else {
            vec![0.0; 2]
        },
        scatter: vec![vec![x(9.0), x(0.1)], vec![x(0.01), x(8.0)]],
        reactions: vec![
            MicroReaction {
                mt: 16,
                channels: vec![vec![x(0.002), 0.0]],
            },
            MicroReaction {
                mt: 102,
                channels: vec![vec![x(0.3), x(2.0)], vec![x(0.03), x(0.2)]],
            },
        ],
        thermal_fission: fissile.then_some(0.8 + scale / 1000.0),
        derivative: None,
    }
}

/// 燃焼度 2 点の `burnup_sample` に、核種 2 つ（核分裂するものとしないもの）のミクロ燃焼のデータを持たせたもの
fn micro_sample() -> Table {
    let mut table = burnup_sample();
    table.micro = Some(Micro {
        nuclides: vec!["U-235".into(), "Xe-135".into()],
        densities: vec![vec![7.0e-4 / 3.0, 0.0], vec![5.0e-4 / 3.0, 1.0e-8 / 3.0]],
        states: (0..table.points.len())
            .map(|p| {
                let p = 1.0 + p as f64;
                vec![micro_nuclide(p, true), micro_nuclide(2.0 * p, false)]
            })
            .collect(),
    });
    table
}

#[test]
fn micro_data_round_trip_bit_exactly_and_are_written_only_when_held() {
    let table = micro_sample();
    let bytes = table.to_msgpack().unwrap();
    assert_eq!(Table::from_msgpack(&bytes).unwrap(), table);
    // 持たないテーブルはキーを書かない（ADR-0216）
    let value = msgpack::read_value(&burnup_sample().to_msgpack().unwrap()).unwrap();
    let keys: Vec<&str> = value
        .as_map()
        .unwrap()
        .iter()
        .filter_map(|(k, _)| k.as_str())
        .collect();
    assert!(!keys.contains(&"micro"), "{keys:?}");
}

#[test]
fn micro_data_need_a_burnup_axis_one_sub_node_and_matching_sizes() {
    let bad = |change: &dyn Fn(&mut Table)| {
        let mut table = micro_sample();
        change(&mut table);
        assert!(table.validate().is_err());
        // 書く前にも確かめる
        assert!(table.to_msgpack().is_err());
    };
    // 燃焼度の軸のないテーブル
    bad(&|t| {
        t.axes.remove(0);
        t.points = product_points(&t.axes);
        t.states.truncate(t.points.len());
        let micro = t.micro.as_mut().unwrap();
        micro.states.truncate(t.points.len());
        micro.densities.truncate(1);
    });
    // 部分ノードに分けたテーブル
    bad(&|t| {
        t.split = 2;
        t.states = t.states.iter().flat_map(|s| vec![s.clone(); 4]).collect();
    });
    // 核種の名前
    bad(&|t| t.micro.as_mut().unwrap().nuclides[1] = "U-235".into());
    bad(&|t| t.micro.as_mut().unwrap().nuclides[1] = "xenon".into());
    bad(&|t| {
        let micro = t.micro.as_mut().unwrap();
        micro.nuclides.clear();
        micro.densities.iter_mut().for_each(Vec::clear);
        micro.states.iter_mut().for_each(Vec::clear);
    });
    // 数密度は燃焼度の点ごとに核種の数だけ、0 以上の有限の値
    bad(&|t| {
        t.micro.as_mut().unwrap().densities.pop();
    });
    bad(&|t| {
        t.micro.as_mut().unwrap().densities[0].pop();
    });
    bad(&|t| t.micro.as_mut().unwrap().densities[1][0] = -1.0e-30);
    bad(&|t| t.micro.as_mut().unwrap().densities[1][0] = f64::NAN);
    // 断面積は点ごとに核種の数だけ、群の数の有限の値
    bad(&|t| {
        t.micro.as_mut().unwrap().states.pop();
    });
    bad(&|t| {
        t.micro.as_mut().unwrap().states[2].pop();
    });
    bad(&|t| {
        t.micro.as_mut().unwrap().states[0][0].total.pop();
    });
    bad(&|t| {
        t.micro.as_mut().unwrap().states[0][0].scatter[1].pop();
    });
    bad(&|t| t.micro.as_mut().unwrap().states[0][1].total[0] = f64::INFINITY);
    bad(&|t| t.micro.as_mut().unwrap().states[0][0].fission[1] = -1.0);
    bad(&|t| {
        t.micro.as_mut().unwrap().states[0][0].reactions[1].channels[0].pop();
    });
    bad(&|t| t.micro.as_mut().unwrap().states[0][0].reactions[1].channels[0][0] = -0.5);
}

#[test]
fn micro_reactions_and_yield_fractions_agree_across_states() {
    let bad = |change: &dyn Fn(&mut Table)| {
        let mut table = micro_sample();
        change(&mut table);
        assert!(table.validate().is_err());
    };
    // 状態の間で内挿するので、核種ごとの反応の MT と経路の数はどの状態でも同じ
    bad(&|t| t.micro.as_mut().unwrap().states[3][0].reactions[0].mt = 17);
    bad(&|t| {
        t.micro.as_mut().unwrap().states[3][0].reactions.pop();
    });
    bad(&|t| {
        t.micro.as_mut().unwrap().states[3][0].reactions[1]
            .channels
            .pop();
    });
    // 1 つの核種に同じ反応が 2 つはない
    bad(&|t| {
        for state in &mut t.micro.as_mut().unwrap().states {
            state[0].reactions[1].mt = 16;
        }
    });
    // 収率の境界より下の割合は、核分裂する核種だけが 0〜1 で持つ
    bad(&|t| t.micro.as_mut().unwrap().states[0][0].thermal_fission = None);
    bad(&|t| t.micro.as_mut().unwrap().states[0][0].thermal_fission = Some(1.5));
    bad(&|t| t.micro.as_mut().unwrap().states[0][1].thermal_fission = Some(0.5));
}

/// 核種 `scale` 倍の値の数密度の微分（2 群。散乱行列の微分は負になりうる）
fn micro_derivative(scale: f64) -> MicroDerivative {
    let x = |v: f64| scale * v / 7.0;
    MicroDerivative {
        total: vec![x(12.0), x(15.0)],
        nu_fission: vec![x(1.3), x(-0.5)],
        fission: vec![x(0.5), x(-0.2)],
        scatter: vec![vec![x(9.5), x(-0.3)], vec![x(0.02), x(8.5)]],
    }
}

/// `micro_sample` の核分裂する核種（U-235）に、すべての状態で数密度の微分を持たせたもの
fn derivative_sample() -> Table {
    let mut table = micro_sample();
    for (p, state) in table.micro.as_mut().unwrap().states.iter_mut().enumerate() {
        state[0].derivative = Some(micro_derivative(1.0 + p as f64));
    }
    table
}

#[test]
fn micro_density_derivatives_round_trip_and_are_written_only_for_the_nuclides_that_hold_them() {
    let table = derivative_sample();
    let bytes = table.to_msgpack().unwrap();
    assert_eq!(Table::from_msgpack(&bytes).unwrap(), table);
    // 微分を持たない核種（Xe-135）は `derivative` のキーを書かない（ADR-0221）
    let value = msgpack::read_value(&bytes).unwrap();
    let keys_of = |nuclide: usize| -> Vec<String> {
        let micro = value
            .as_map()
            .unwrap()
            .iter()
            .find(|(k, _)| k.as_str() == Some("micro"))
            .unwrap()
            .1
            .clone();
        let states = micro
            .as_map()
            .unwrap()
            .iter()
            .find(|(k, _)| k.as_str() == Some("states"))
            .unwrap()
            .1
            .clone();
        states.as_array().unwrap()[0].as_array().unwrap()[nuclide]
            .as_map()
            .unwrap()
            .iter()
            .filter_map(|(k, _)| k.as_str().map(String::from))
            .collect()
    };
    assert!(
        keys_of(0).contains(&"derivative".to_string()),
        "{:?}",
        keys_of(0)
    );
    assert!(
        !keys_of(1).contains(&"derivative".to_string()),
        "{:?}",
        keys_of(1)
    );
    // 微分を持たないテーブルは、微分の項目を足す前と同じバイト
    let plain = micro_sample();
    assert_eq!(
        Table::from_msgpack(&plain.to_msgpack().unwrap()).unwrap(),
        plain
    );
}

#[test]
fn an_explicit_nil_density_derivative_is_an_error() {
    // 微分を持たない核種はキーを書かない（ADR-0221）。nil は書き手の誤りなので、持たないと読まない
    let bytes = micro_sample().to_msgpack().unwrap();
    let mut value = msgpack::read_value(&bytes).unwrap();
    let rmpv::Value::Map(top) = &mut value else {
        panic!("a table is a map")
    };
    let (_, micro) = top
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("micro"))
        .unwrap();
    let rmpv::Value::Map(micro) = micro else {
        panic!("micro is a map")
    };
    let (_, states) = micro
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("states"))
        .unwrap();
    let rmpv::Value::Array(points) = states else {
        panic!("micro states are an array")
    };
    let rmpv::Value::Array(nuclides) = &mut points[0] else {
        panic!("a micro state is an array")
    };
    let rmpv::Value::Map(nuclide) = &mut nuclides[1] else {
        panic!("a micro nuclide is a map")
    };
    nuclide.push((rmpv::Value::from("derivative"), rmpv::Value::Nil));
    let mut written = Vec::new();
    rmpv::encode::write_value(&mut written, &value).unwrap();
    assert!(Table::from_msgpack(&written).is_err());
}

#[test]
fn micro_density_derivatives_have_the_groups_and_are_held_in_every_state() {
    let bad = |change: &dyn Fn(&mut Table)| {
        let mut table = derivative_sample();
        change(&mut table);
        assert!(table.validate().is_err());
        assert!(table.to_msgpack().is_err());
    };
    fn derivative(t: &mut Table, p: usize) -> &mut MicroDerivative {
        t.micro.as_mut().unwrap().states[p][0]
            .derivative
            .as_mut()
            .unwrap()
    }
    bad(&|t| {
        derivative(t, 0).total.pop();
    });
    bad(&|t| {
        derivative(t, 1).nu_fission.push(0.0);
    });
    bad(&|t| {
        derivative(t, 2).fission.pop();
    });
    bad(&|t| {
        derivative(t, 0).scatter.pop();
    });
    bad(&|t| {
        derivative(t, 0).scatter[1].pop();
    });
    bad(&|t| derivative(t, 3).total[1] = f64::NAN);
    bad(&|t| derivative(t, 3).scatter[0][1] = f64::INFINITY);
    // 状態の間で内挿するので、微分を持つ核種はどの状態でも持つ
    bad(&|t| t.micro.as_mut().unwrap().states[2][0].derivative = None);
    bad(&|t| t.micro.as_mut().unwrap().states[1][1].derivative = Some(micro_derivative(1.0)));
    // 負の微分（組成が増えると減る少数群定数）は持てる
    let table = derivative_sample();
    assert!(
        table.micro.as_ref().unwrap().states[0][0]
            .derivative
            .as_ref()
            .unwrap()
            .nu_fission[1]
            < 0.0
    );
    table.validate().unwrap();
}

/// 値を `scale` 倍した残り（2 群、6 つの辺）
fn residual(scale: f64) -> Residual {
    let x = |v: f64| scale * v / 3.0;
    Residual {
        total: vec![x(0.01), x(-0.02)],
        absorption: vec![x(0.001), x(-0.003)],
        nu_fission: vec![x(-0.002), x(0.004)],
        fission: vec![x(-0.0008), x(0.0016)],
        chi: vec![0.0, 0.0],
        diffusion: vec![x(0.05), x(-0.01)],
        scatter: vec![vec![x(0.009), x(-0.0005)], vec![x(0.00002), x(-0.017)]],
        adf: vec![vec![x(0.001), x(-0.002)]; 6],
    }
}

/// `micro_nuclide` と同じ反応の並びの、値を `scale` 倍した反応率の微視的断面積の差（2 群。負の値を含む）
fn rate_residual(scale: f64, fissile: bool) -> RateResidual {
    let x = |v: f64| scale * v / 3.0;
    RateResidual {
        fission: if fissile {
            vec![x(-0.004), x(0.29)]
        } else {
            vec![0.0; 2]
        },
        reactions: vec![
            MicroReaction {
                mt: 16,
                channels: vec![vec![x(0.00002), 0.0]],
            },
            MicroReaction {
                mt: 102,
                channels: vec![vec![x(-0.003), x(0.02)], vec![x(0.0003), x(-0.002)]],
            },
        ],
        thermal_fission: fissile.then_some(x(-0.01)),
    }
}

/// 公称のテーブル: 燃焼度 0、10、20 MWd/kg、ホウ素 1 点、減速材 2 点（公称の 565 K と 0.743 g/cm³、履歴の分岐の条件の 600 K と
/// 0.661 g/cm³）、燃料 1 点。公称の状態はすべての燃焼度で、減速材の公称でない状態は 0 と 20 で解き、追う核種 1 つと、制御棒のない
/// 減速材の密度の履歴の補正（0 と 20 の点。追わない核種の分、線形化の誤差、履歴の分岐の均質化した数密度）を持つ
fn history_sample() -> Table {
    let axes = vec![
        Axis::Burnup {
            points: vec![0.0, 10.0, 20.0],
        },
        Axis::boron(vec![1300.0], 0.199).with_nominal(1300.0),
        Axis::moderator(vec![565.0, 600.0], vec![0.743, 0.661]).with_nominal(565.0),
        Axis::fuel(vec![565.0]).with_nominal(565.0),
    ];
    let points: Vec<Vec<usize>> = product_points(&axes)
        .into_iter()
        .filter(|p| p[0] != 1 || p[2] == 0)
        .collect();
    let n = points.len();
    Table {
        axes,
        split: 1,
        states: (0..n).map(|i| state(3.0 + i as f64)).collect(),
        heavy_metal: Some(400.0 / 3.0),
        micro: Some(Micro {
            nuclides: vec!["U-235".into()],
            densities: vec![vec![7.0e-4 / 3.0], vec![6.0e-4 / 3.0], vec![5.0e-4 / 3.0]],
            states: (0..n)
                .map(|p| vec![micro_nuclide(1.0 + p as f64, true)])
                .collect(),
        }),
        histories: vec![HistoryCorrection {
            variable: MODERATOR_DENSITY.into(),
            nominal: 0.743,
            branch: 0.661,
            rods: vec![HistoryRod {
                rod: ROD_OUT.into(),
                burnup: vec![0.0, 20.0],
                untracked: vec![residual(1.0), residual(2.0)],
                // 燃焼度 0 は履歴の分岐の数密度が公称と同じで、線形化の誤差は 2 つの格子計算の解の差の大きさ
                linearization: vec![residual(3.0e-9), residual(4.0)],
                densities: vec![vec![7.0e-4 / 3.0], vec![5.5e-4 / 3.0]],
                micro: RateResiduals {
                    untracked: vec![
                        vec![rate_residual(1.0, true)],
                        vec![rate_residual(2.0, true)],
                    ],
                    linearization: vec![
                        vec![rate_residual(-0.5e-9, true)],
                        vec![rate_residual(0.75, true)],
                    ],
                },
            }],
        }],
        points,
    }
}

/// 変えると検査で落ちる（書く前にも確かめる）。エラーの文は `expected` を含む
fn bad_history(change: &dyn Fn(&mut Table), expected: &str) {
    let mut table = history_sample();
    change(&mut table);
    let err = table.validate().unwrap_err().to_string();
    assert!(err.contains(expected), "{expected:?} not in {err}");
    // 折り返した文の空白が文に入らない
    assert!(!err.contains("  "), "{err}");
    assert!(table.to_msgpack().is_err());
}

#[test]
fn history_corrections_round_trip_bit_exactly_and_are_written_only_when_held() {
    let table = history_sample();
    table.validate().unwrap();
    let bytes = table.to_msgpack().unwrap();
    assert_eq!(Table::from_msgpack(&bytes).unwrap(), table);
    // 持たないテーブルはキーを書かない（ADR-0215）
    let mut none = history_sample();
    none.histories.clear();
    let value = msgpack::read_value(&none.to_msgpack().unwrap()).unwrap();
    let keys: Vec<&str> = value
        .as_map()
        .unwrap()
        .iter()
        .filter_map(|(k, _)| k.as_str())
        .collect();
    assert!(!keys.contains(&"histories"), "{keys:?}");
}

#[test]
fn history_corrections_need_micro_data_and_a_known_distinct_variable() {
    // R はミクロ燃焼の線形の補正を引いた残りなので、補正のデータを持つテーブルだけが持てる
    bad_history(&|t| t.micro = None, "micro");
    // 履歴の変数は減速材の密度だけ（制御棒と Pyrex は、それぞれの履歴の分岐を作るときに足す）
    bad_history(&|t| t.histories[0].variable = "boron".into(), "variable");
    bad_history(
        &|t| {
            let again = t.histories[0].clone();
            t.histories.push(again);
        },
        "twice",
    );
}

#[test]
fn history_values_are_the_moderator_densities_of_the_nominal_and_another_point() {
    // H⁰ は公称の格子点の密度
    bad_history(&|t| t.histories[0].nominal = 0.75, "nominal");
    // Hᵇ は公称でない格子点の密度（その格子点が履歴の分岐の条件）
    bad_history(&|t| t.histories[0].branch = 0.70, "branch");
    bad_history(&|t| t.histories[0].branch = 0.743, "branch");
    // 公称と同じ密度の格子点（温度だけが違う）は変数の差を作らない（core の重み (H − H⁰)/(Hᵇ − H⁰) が 0 で割る）
    bad_history(
        &|t| {
            t.axes[2] = Axis::moderator(vec![565.0, 600.0], vec![0.743, 0.743]).with_nominal(565.0);
            t.histories[0].branch = 0.743;
        },
        "branch",
    );
    // 公称の値のない（直積の）テーブルには、履歴の分岐の条件を決める公称の状態がない
    bad_history(
        &|t| {
            for axis in &mut t.axes {
                *axis = axis.clone().without_nominal();
            }
            t.points = product_points(&t.axes);
            t.states = (0..t.points.len()).map(|i| state(3.0 + i as f64)).collect();
            let micro = t.micro.as_mut().unwrap();
            micro.states = (0..t.points.len())
                .map(|p| vec![micro_nuclide(1.0 + p as f64, true)])
                .collect();
        },
        "nominal",
    );
}

#[test]
fn history_rods_are_distinct_values_of_the_rod_axis() {
    bad_history(&|t| t.histories[0].rods.clear(), "rod");
    bad_history(
        &|t| {
            let again = t.histories[0].rods[0].clone();
            t.histories[0].rods.push(again);
        },
        "rod",
    );
    // 制御棒の軸のないテーブルの値は out だけ
    bad_history(&|t| t.histories[0].rods[0].rod = "aic".into(), "rod");
}

#[test]
fn history_burnups_are_solved_points_of_the_condition_from_first_to_last() {
    // 公称でない状態を解いていない燃焼度の点（10 MWd/kg）の定数は内挿で、格子計算の値でない
    bad_history(
        &|t| {
            let rod = &mut t.histories[0].rods[0];
            rod.burnup = vec![0.0, 10.0, 20.0];
            rod.untracked.push(residual(3.0));
            rod.linearization.push(residual(3.0));
            rod.densities.push(vec![5.0e-4 / 3.0]);
        },
        "solved",
    );
    // 燃焼度の軸の点でない
    bad_history(&|t| t.histories[0].rods[0].burnup[1] = 15.0, "burnup");
    // core は外挿しないので、最初と最後の点を含む
    bad_history(
        &|t| {
            let rod = &mut t.histories[0].rods[0];
            rod.burnup.remove(0);
            rod.untracked.remove(0);
            rod.linearization.remove(0);
            rod.densities.remove(0);
        },
        "first and the last",
    );
    bad_history(
        &|t| {
            let rod = &mut t.histories[0].rods[0];
            rod.burnup.reverse();
        },
        "burnup",
    );
    // 点ごとに 1 つの追わない核種の分、線形化の誤差、数密度の列
    bad_history(
        &|t| {
            t.histories[0].rods[0].untracked.pop();
        },
        "untracked",
    );
    bad_history(
        &|t| {
            t.histories[0].rods[0].linearization.pop();
        },
        "linearization",
    );
    bad_history(
        &|t| {
            t.histories[0].rods[0].densities.pop();
        },
        "densities",
    );
}

#[test]
fn residuals_have_the_groups_and_sides_of_the_states_and_are_finite() {
    let cases: Vec<fn(&mut Residual)> = vec![
        |r| {
            r.total.pop();
        },
        |r| r.absorption.push(0.0),
        |r| {
            r.nu_fission.pop();
        },
        |r| {
            r.fission.pop();
        },
        |r| {
            r.chi.pop();
        },
        |r| {
            r.diffusion.pop();
        },
        |r| {
            r.scatter.pop();
        },
        |r| {
            r.scatter[1].pop();
        },
        |r| {
            r.adf.pop();
        },
        |r| {
            r.adf[2].pop();
        },
        |r| r.total[0] = f64::NAN,
        |r| r.scatter[0][1] = f64::INFINITY,
        |r| r.adf[5][1] = f64::NAN,
    ];
    for change in &cases {
        bad_history(
            &|t| change(&mut t.histories[0].rods[0].untracked[1]),
            "untracked",
        );
        bad_history(
            &|t| change(&mut t.histories[0].rods[0].linearization[1]),
            "linearization",
        );
    }
    // 残りは負になりうる
    history_sample().validate().unwrap();
}

#[test]
fn history_densities_are_finite_non_negative_rows_of_the_tracked_nuclides() {
    // core が線形化の誤差の指標を作る、履歴の分岐の均質化した数密度 N̄ᵢ^j（追う核種の順）
    bad_history(
        &|t| t.histories[0].rods[0].densities[1].push(0.0),
        "densities",
    );
    bad_history(
        &|t| t.histories[0].rods[0].densities[1].clear(),
        "densities",
    );
    bad_history(
        &|t| t.histories[0].rods[0].densities[1][0] = f64::NAN,
        "densities",
    );
    bad_history(
        &|t| t.histories[0].rods[0].densities[1][0] = -1.0e-9,
        "densities",
    );
}

/// 反応率の微視的断面積の差の 1 つの部分（追わない核種の分か線形化の誤差）を変える関数と、検査のエラーの文に入る語
type RatePartChange = (fn(&mut Vec<Vec<RateResidual>>), &'static str);

/// 制御棒の値の残りの、反応率の微視的断面積の差の部分を選ぶ関数と、その部分の名前
type RatePart = (
    fn(&mut HistoryRod) -> &mut Vec<Vec<RateResidual>>,
    &'static str,
);

const RATE_PARTS: [RatePart; 2] = [
    (|r| &mut r.micro.untracked, "untracked"),
    (|r| &mut r.micro.linearization, "linearization"),
];

#[test]
fn rate_residuals_match_the_micro_nuclides_and_their_reactions() {
    // 追わない核種の分と線形化の誤差のどちらも、同じ形の検査を受け、エラーの文は部分の名前を含む
    let cases: Vec<RatePartChange> = vec![
        // 点ごとに 1 つの列
        (
            |m| {
                m.pop();
            },
            "rate",
        ),
        // 列は追う核種の数
        (|m| m[1].push(rate_residual(1.0, true)), "rate"),
        (
            |m| {
                m[0].pop();
            },
            "rate",
        ),
        // 反応の並び（MT と経路の数）と核分裂の有無はミクロ燃焼のデータと同じ
        (|m| m[1][0].reactions[0].mt = 17, "reactions"),
        (
            |m| {
                m[1][0].reactions[1].channels.pop();
            },
            "reactions",
        ),
        (
            |m| {
                m[1][0].reactions.pop();
            },
            "reactions",
        ),
        (|m| m[0][0].thermal_fission = None, "fission"),
        // 群の数の有限の値
        (
            |m| {
                m[1][0].fission.pop();
            },
            "groups",
        ),
        (
            |m| {
                m[1][0].reactions[1].channels[0].push(0.0);
            },
            "groups",
        ),
        (|m| m[1][0].fission[0] = f64::NAN, "groups"),
        (
            |m| m[0][0].reactions[0].channels[0][1] = f64::INFINITY,
            "groups",
        ),
        (|m| m[1][0].thermal_fission = Some(f64::NAN), "fission"),
    ];
    for (part, name) in RATE_PARTS {
        for (change, expected) in &cases {
            bad_history(&|t| change(part(&mut t.histories[0].rods[0])), expected);
            bad_history(&|t| change(part(&mut t.histories[0].rods[0])), name);
        }
    }
    // 核分裂しない核種の差は nil
    let mut table = history_sample();
    let micro = table.micro.as_mut().unwrap();
    micro.nuclides.push("Xe-135".into());
    for (p, state) in micro.states.iter_mut().enumerate() {
        state.push(micro_nuclide(2.0 + p as f64, false));
    }
    for row in &mut micro.densities {
        row.push(1.0e-9);
    }
    for (part, name) in RATE_PARTS {
        for (k, row) in part(&mut table.histories[0].rods[0]).iter_mut().enumerate() {
            // 線形化の誤差の燃焼度 0 は、数密度が公称と同じ点の 0 の大きさ
            let tiny = name == "linearization" && k == 0;
            row.push(rate_residual(if tiny { 3.0e-9 } else { 3.0 }, false));
        }
    }
    for row in &mut table.histories[0].rods[0].densities {
        row.push(1.0e-9);
    }
    table.validate().unwrap();
    for (part, _) in RATE_PARTS {
        let mut changed = table.clone();
        part(&mut changed.histories[0].rods[0])[1][1].thermal_fission = Some(0.0);
        let err = changed.validate().unwrap_err().to_string();
        assert!(err.contains("fission"), "{err}");
    }
}

#[test]
fn branch_points_at_the_nominal_densities_hold_no_linearization() {
    // 履歴の分岐の追う核種の数密度が公称と同じ点では、ハイブリッドは公称の組成そのもので、線形化の誤差 L と Δσ_L は 2 つの格子計算の
    // 解の差の大きさ（core はその点の L と Δσ_L を足さないので、ほかの組成から作った値を黙って捨てない）
    history_sample().validate().unwrap();
    bad_history(
        &|t| t.histories[0].rods[0].linearization[0] = residual(3.0),
        "nominal densities",
    );
    bad_history(
        &|t| t.histories[0].rods[0].micro.linearization[0] = vec![rate_residual(1.0, true)],
        "nominal densities",
    );
    // 収率の境界より下の核分裂の断面積の差も
    bad_history(
        &|t| {
            t.histories[0].rods[0].micro.linearization[0][0].thermal_fission = Some(0.01);
        },
        "nominal densities",
    );
    // 追わない核種の分は、同じ組成でも 0 とは限らない（燃焼度 0 の residual(1.0)）。数密度が公称と違う点の L は大きくてよい
    let mut table = history_sample();
    table.histories[0].rods[0].linearization[1] = residual(400.0);
    table.validate().unwrap();
}

#[test]
fn each_linearization_item_at_the_nominal_densities_is_checked() {
    // L と Δσ_L の項目のどれか 1 つだけが 2 つの格子計算の解の差より大きくても止める（燃焼度 0 の点。ほかの項目は 0 の大きさのまま）
    type Change = fn(&mut HistoryRod);
    let items: [(&str, Change); 12] = [
        ("total", |r| r.linearization[0].total[1] = 1.0e-3),
        ("absorption", |r| r.linearization[0].absorption[0] = 1.0e-3),
        ("nu_fission", |r| r.linearization[0].nu_fission[1] = 1.0e-3),
        ("fission", |r| r.linearization[0].fission[0] = 1.0e-3),
        ("chi", |r| r.linearization[0].chi[0] = 1.0e-3),
        ("diffusion", |r| r.linearization[0].diffusion[1] = 1.0e-3),
        ("scatter", |r| r.linearization[0].scatter[1][0] = 1.0e-3),
        ("adf", |r| r.linearization[0].adf[0][0] = 1.0e-3),
        ("micro fission", |r| {
            r.micro.linearization[0][0].fission[1] = 1.0e-3
        }),
        ("micro reactions", |r| {
            r.micro.linearization[0][0].reactions[0].channels[0][0] = 1.0e-3;
        }),
        ("micro capture channel", |r| {
            r.micro.linearization[0][0].reactions[1].channels[1][1] = 1.0e-3;
        }),
        ("micro thermal fission", |r| {
            r.micro.linearization[0][0].thermal_fission = Some(1.0e-3);
        }),
    ];
    for (name, change) in items {
        let mut table = history_sample();
        change(&mut table.histories[0].rods[0]);
        assert!(table.validate().is_err(), "{name}: not checked");
        bad_history(
            &|t| change(&mut t.histories[0].rods[0]),
            "nominal densities",
        );
    }
}

#[test]
fn linearizations_at_the_nominal_densities_are_bounded_by_the_largest_value_of_each_item() {
    // 許容は項目ごとに、条件 j の状態の同じ項目の絶対値の最大の `NOMINAL_POINT_TOLERANCE` 倍（境界の内は通り、外は止める）
    let condition = |t: &Table| {
        let p = t.point_index(&[0, 0, 1, 0]).unwrap();
        (
            t.states[p].clone(),
            t.micro.as_ref().unwrap().states[p][0].clone(),
        )
    };
    let largest = |x: &[f64]| x.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let (f, sigma) = condition(&history_sample());
    let total = largest(&f.total);
    type Set = fn(&mut HistoryRod, f64);
    let items: [(&str, f64, Set); 4] = [
        ("total", total, |r, x| r.linearization[0].total[0] = x),
        // 散乱は行列の全体の最大（行の最大でない）。0 行目の最大より大きい値で、1 行目の最大の許容の内は通る
        ("scatter", largest(&f.scatter.concat()), |r, x| {
            r.linearization[0].scatter[0][1] = x
        }),
        // 捕獲は経路ごと
        (
            "capture",
            largest(&sigma.reactions[1].channels[0]),
            |r, x| r.micro.linearization[0][0].reactions[1].channels[0][1] = -x,
        ),
        // 収率の境界より下の核分裂の差は割合 × 群の和の σf なので、基準は群の和の σf（群の最大でない）
        ("thermal fission", sigma.fission.iter().sum(), |r, x| {
            r.micro.linearization[0][0].thermal_fission = Some(x)
        }),
    ];
    assert!(largest(&f.scatter[0]) < 0.99 * largest(&f.scatter.concat()));
    assert!(largest(&sigma.fission) < 0.99 * sigma.fission.iter().sum::<f64>());
    for (name, scale, set) in items {
        for (factor, fine) in [(0.99, true), (1.01, false)] {
            let mut table = history_sample();
            set(
                &mut table.histories[0].rods[0],
                factor * NOMINAL_POINT_TOLERANCE * scale,
            );
            assert_eq!(table.validate().is_ok(), fine, "{name} {factor}");
        }
    }
    // 基準は項目ごとの最大: 吸収の最大の許容を超える吸収は、全断面積の最大の許容に入っても止める
    let absorption = largest(&f.absorption);
    let value = 2.0 * NOMINAL_POINT_TOLERANCE * absorption;
    assert!(value < NOMINAL_POINT_TOLERANCE * total, "{value} {total}");
    bad_history(
        &|t| t.histories[0].rods[0].linearization[0].absorption[1] = value,
        "nominal densities",
    );
}

#[test]
fn linearizations_at_the_nominal_densities_of_items_that_are_zero_must_be_zero() {
    // 基準が 0 の項目（核分裂しない核種の σf）は、0 でない差を許さない
    let mut table = history_sample();
    let micro = table.micro.as_mut().unwrap();
    micro.nuclides.push("Xe-135".into());
    for (p, state) in micro.states.iter_mut().enumerate() {
        state.push(micro_nuclide(2.0 + p as f64, false));
    }
    for row in &mut micro.densities {
        row.push(1.0e-9);
    }
    let rod = &mut table.histories[0].rods[0];
    for part in [&mut rod.micro.untracked, &mut rod.micro.linearization] {
        for row in part.iter_mut() {
            row.push(rate_residual(0.0, false));
        }
    }
    for row in &mut rod.densities {
        row.push(1.0e-9);
    }
    table.validate().unwrap();
    table.histories[0].rods[0].micro.linearization[0][1].fission[0] = 1.0e-30;
    let err = table.validate().unwrap_err().to_string();
    assert!(err.contains("nominal densities"), "{err}");
}

#[test]
fn rate_residuals_of_one_unsplit_list_are_not_read() {
    // 反応率の微視的断面積の差は、追わない核種の分 `untracked` と線形化の誤差 `linearization` の map。分けない 1 つの列（前の形）は読まない
    let bytes = history_sample().to_msgpack().unwrap();
    let mut value = msgpack::read_value(&bytes).unwrap();
    let rmpv::Value::Map(entries) = &mut value else {
        panic!()
    };
    let (_, histories) = entries
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("histories"))
        .unwrap();
    let rmpv::Value::Array(histories) = histories else {
        panic!()
    };
    let rmpv::Value::Map(history) = &mut histories[0] else {
        panic!()
    };
    let (_, rods) = history
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("rods"))
        .unwrap();
    let rmpv::Value::Array(rods) = rods else {
        panic!()
    };
    let rmpv::Value::Map(rod) = &mut rods[0] else {
        panic!()
    };
    let (_, micro) = rod
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("micro"))
        .unwrap();
    let rmpv::Value::Map(parts) = micro else {
        panic!("{micro:?}")
    };
    let mut keys: Vec<&str> = parts.iter().filter_map(|(k, _)| k.as_str()).collect();
    keys.sort_unstable();
    assert_eq!(keys, ["linearization", "untracked"]);
    let untracked = parts
        .iter()
        .find(|(k, _)| k.as_str() == Some("untracked"))
        .unwrap()
        .1
        .clone();
    *micro = untracked;
    let mut out = Vec::new();
    rmpv::encode::write_value(&mut out, &value).unwrap();
    assert!(Table::from_msgpack(&out).is_err());
}

#[test]
fn linear_corrections_add_the_coefficients_times_the_density_differences() {
    // 係数は、数密度の微分を持つ核種はその微分、ほかは σᵢ。吸収の係数は全断面積の係数から散乱の行の和を引いたもの
    let plain = micro_nuclide(1.0, true);
    let mut derived = micro_nuclide(2.0, true);
    derived.derivative = Some(MicroDerivative {
        total: vec![5.0 / 3.0, 7.0 / 3.0],
        nu_fission: vec![0.5 / 3.0, -4.0 / 3.0],
        fission: vec![0.2 / 3.0, -1.5 / 3.0],
        scatter: vec![vec![4.0 / 3.0, 0.2 / 3.0], vec![0.01 / 3.0, 6.0 / 3.0]],
    });
    let base = state(3.0);
    let delta = [1.0e-5 / 3.0, -2.0e-6 / 7.0];
    let mut f = base.clone();
    add_linear_correction(&mut f, &[plain.clone(), derived.clone()], &delta);
    let c = derived.derivative.as_ref().unwrap();
    for g in 0..2 {
        let expected = |base: f64, a: f64, b: f64| base + delta[0] * a + delta[1] * b;
        assert_eq!(
            f.total[g],
            expected(base.total[g], plain.total[g], c.total[g])
        );
        let removal = |t: &[f64], s: &[Vec<f64>]| t[g] - s[g].iter().sum::<f64>();
        assert_eq!(
            f.absorption[g],
            expected(
                base.absorption[g],
                removal(&plain.total, &plain.scatter),
                removal(&c.total, &c.scatter)
            )
        );
        assert_eq!(
            f.nu_fission[g],
            expected(base.nu_fission[g], plain.nu_fission[g], c.nu_fission[g])
        );
        assert_eq!(
            f.fission[g],
            expected(base.fission[g], plain.fission[g], c.fission[g])
        );
        for h in 0..2 {
            assert_eq!(
                f.scatter[g][h],
                expected(base.scatter[g][h], plain.scatter[g][h], c.scatter[g][h])
            );
        }
    }
    // 補正を持たない項目は変えない
    assert_eq!(
        (&f.chi, &f.diffusion, &f.adf, f.k_inf),
        (&base.chi, &base.diffusion, &base.adf, base.k_inf)
    );
    // 差が 0 なら、定数はビットで同じ
    let mut same = base.clone();
    add_linear_correction(&mut same, &[plain, derived], &[0.0, 0.0]);
    assert_eq!(same, base);
}

#[test]
fn history_conditions_are_the_branch_point_of_the_variable_axis_and_nominal_elsewhere() {
    let table = history_sample();
    let condition = table
        .history_condition(MODERATOR_DENSITY, 0.661, ROD_OUT)
        .unwrap();
    assert_eq!(
        condition,
        HistoryCondition {
            nominal: 0.743,
            point: vec![0, 0, 1, 0],
            burnup_axis: 0,
        }
    );
    let fails = |table: &Table, variable: &str, branch: f64, rod: &str, expected: &str| {
        let err = table.history_condition(variable, branch, rod).unwrap_err();
        assert!(err.contains(expected), "{expected:?} not in {err}");
        assert!(!err.contains("  "), "{err}");
    };
    fails(&table, "boron", 0.661, ROD_OUT, "variable");
    // Hᵇ は公称でない 1 つの格子点の密度
    fails(&table, MODERATOR_DENSITY, 0.70, ROD_OUT, "branch");
    fails(&table, MODERATOR_DENSITY, 0.743, ROD_OUT, "branch");
    // 制御棒の軸のないテーブルの値は out だけ
    fails(&table, MODERATOR_DENSITY, 0.661, "aic", "rod");
    // 公称の値のないテーブルには、条件を決める公称の状態がない
    let mut direct = table.clone();
    for axis in &mut direct.axes {
        *axis = axis.clone().without_nominal();
    }
    fails(&direct, MODERATOR_DENSITY, 0.661, ROD_OUT, "nominal");
    // 制御棒の軸を持てば、その値の番号
    let mut rodded = table.clone();
    rodded.axes.push(Axis::Rod {
        values: vec![ROD_OUT.into(), "aic".into()],
    });
    let condition = rodded
        .history_condition(MODERATOR_DENSITY, 0.661, "aic")
        .unwrap();
    assert_eq!(condition.point, vec![0, 0, 1, 0, 1]);
    fails(&rodded, MODERATOR_DENSITY, 0.661, "b4c", "rod");
    // 燃焼度の軸が要る
    let mut fresh = table;
    fresh.axes.remove(0);
    fails(&fresh, MODERATOR_DENSITY, 0.661, ROD_OUT, "burnup");
}
