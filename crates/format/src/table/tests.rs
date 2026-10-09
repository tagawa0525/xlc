use super::*;

fn state(k: f64) -> FewGroup {
    // 10 進で短く書けない値にして、ビット単位の往復を確かめる
    let x = |v: f64| v / 3.0;
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: x(k),
        total: vec![x(1.6), x(2.2)],
        absorption: vec![x(0.03), x(0.3)],
        nu_fission: vec![x(0.02), x(0.4)],
        fission: vec![x(0.008), x(0.16)],
        chi: vec![1.0, 0.0],
        diffusion: vec![x(4.1), x(1.2)],
        scatter: vec![vec![x(1.5), x(0.05)], vec![x(0.002), x(1.9)]],
        adf: vec![vec![x(3.0), x(3.3)]; 4],
        vacuum_ratio: vec![None, Some(vec![x(1.6), x(1.7)]), None, None],
        normalized_adf: vec![true, false, false, false],
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
fn states_follow_the_axes_with_the_first_slowest() {
    let table = sample();
    // (ホウ素, 減速材, 燃料) の番号から状態の番号
    assert_eq!(table.index(&[0, 0, 0]), Some(0));
    assert_eq!(table.index(&[0, 1, 0]), Some(1));
    assert_eq!(table.index(&[1, 0, 0]), Some(2));
    assert_eq!(table.index(&[1, 1, 0]), Some(3));
    assert_eq!(table.index(&[2, 0, 0]), None);
    assert_eq!(table.index(&[0, 0]), None);
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
    // 群境界、f* の印が状態の間で違う
    bad(|t| t.states[3].bounds = vec![2.0e7, 0.625, 1.0e-5]);
    bad(|t| t.states[3].normalized_adf = vec![false; 4]);
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
