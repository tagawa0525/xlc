use super::*;

fn state(k: f64) -> FewGroup {
    FewGroup {
        k_inf: k / 3.0,
        ..crate::fewgroup::sample()
    }
}

/// ホウ素 2 点 × 減速材 2 点 × 燃料 1 点
fn sample() -> Table {
    Table {
        axes: vec![
            Axis::Boron {
                ppm: vec![900.0, 1700.0],
                b10: 0.199,
            },
            Axis::Moderator {
                temperature: vec![550.0, 580.0],
                density: vec![0.76972, 0.71190],
            },
            Axis::Fuel {
                temperature: vec![565.0],
            },
        ],
        split: 1,
        states: (0..4).map(|i| state(3.0 + f64::from(i))).collect(),
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
    bad(|t| {
        t.axes[0] = Axis::Boron {
            ppm: vec![1700.0, 900.0],
            b10: 0.199,
        }
    });
    bad(|t| {
        t.axes[2] = Axis::Fuel {
            temperature: vec![],
        }
    });
    bad(|t| {
        t.axes[2] = Axis::Fuel {
            temperature: vec![f64::NAN],
        }
    });
    bad(|t| {
        t.axes[2] = Axis::Fuel {
            temperature: vec![-1.0],
        }
    });
    bad(|t| {
        t.axes[0] = Axis::Boron {
            ppm: vec![-1.0, 900.0],
            b10: 0.199,
        }
    });
    // 密度の数が点の数と違う、密度が正でない、B-10 の割合が範囲の外
    bad(|t| {
        t.axes[1] = Axis::Moderator {
            temperature: vec![550.0, 580.0],
            density: vec![0.76972],
        }
    });
    bad(|t| {
        t.axes[1] = Axis::Moderator {
            temperature: vec![550.0, 580.0],
            density: vec![0.76972, 0.0],
        }
    });
    bad(|t| {
        t.axes[0] = Axis::Boron {
            ppm: vec![900.0, 1700.0],
            b10: 1.5,
        }
    });
    // 同じ軸が2つ
    bad(|t| {
        t.axes[2] = Axis::Boron {
            ppm: vec![1300.0],
            b10: 0.199,
        }
    });
    // 群境界、f* の印、真空の辺が状態の間で違う
    bad(|t| t.states[3].bounds = vec![2.0e7, 0.625, 1.0e-5]);
    bad(|t| t.states[3].normalized_adf = vec![false; 4]);
    bad(|t| t.states[3].vacuum_ratio[1] = None);
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
            let mut pins = crate::fewgroup::sample_pins();
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
    assert!(err.contains("state"), "{err}");
}
