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
