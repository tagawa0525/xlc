use super::*;

/// 3群、P0 と P1 の散乱、核分裂のある状態（10 進で短く書けない値にして、ビット単位の往復を確かめる）
fn state(scale: f64) -> MediumState {
    let x = |v: f64| scale * v / 3.0;
    MediumState {
        total: vec![x(0.6), x(1.1), x(2.0)],
        scatter: vec![
            vec![
                vec![x(0.5), x(0.08), x(0.001)],
                vec![0.0, x(1.0), x(0.05)],
                vec![0.0, x(0.002), x(1.9)],
            ],
            vec![
                vec![x(0.2), x(0.01), 0.0],
                vec![0.0, x(0.3), x(0.004)],
                vec![0.0, 0.0, x(0.4)],
            ],
        ],
        fission: vec![x(0.002), x(0.02), x(0.08)],
        nu_fission: vec![x(0.005), x(0.05), x(0.2)],
        chi: vec![0.75, 0.25, 0.0],
    }
}

/// ホウ素 2 点 × 減速材 1 点 × 燃料 2 点
fn sample() -> Medium {
    Medium {
        bounds: vec![2.0e7, 1.0e5, 0.625, 1.0e-5],
        axes: vec![
            Axis::boron(vec![900.0, 1700.0], 0.199),
            Axis::moderator(vec![565.0], vec![0.743]),
            Axis::fuel(vec![565.0, 900.0]),
        ],
        states: (1..=4).map(|i| state(f64::from(i))).collect(),
    }
}

#[test]
fn round_trips_bit_exactly() {
    let medium = sample();
    assert_eq!(
        Medium::from_msgpack(&medium.to_msgpack().unwrap()).unwrap(),
        medium
    );
    // 分岐のない媒質は状態が1つ
    let single = Medium {
        axes: vec![],
        states: vec![state(1.0)],
        ..sample()
    };
    assert_eq!(
        Medium::from_msgpack(&single.to_msgpack().unwrap()).unwrap(),
        single
    );
}

#[test]
fn invalid_media_are_errors() {
    let bad = |change: fn(&mut Medium), word: &str| {
        let mut medium = sample();
        change(&mut medium);
        // 書くときも読むときも確かめる
        let err = medium.validate().unwrap_err().to_string();
        assert!(err.contains(word), "{err}");
        assert!(medium.to_msgpack().is_err());
    };
    bad(|m| m.bounds.reverse(), "bounds");
    bad(|m| m.bounds.truncate(1), "bounds");
    // 状態の数が軸の点の数の積と違う、軸の誤り（テーブルと同じ検査）
    bad(
        |m| {
            m.states.pop();
        },
        "states",
    );
    bad(|m| m.axes[2] = Axis::fuel(vec![900.0, 565.0]), "axis");
    // 群の数と違う長さ、有限でない値、負の断面積
    bad(
        |m| {
            m.states[1].total.pop();
        },
        "total",
    );
    bad(
        |m| {
            m.states[1].scatter[1][2].pop();
        },
        "scatter",
    );
    bad(
        |m| {
            m.states[1].scatter[0].pop();
        },
        "scatter",
    );
    bad(|m| m.states[2].fission[0] = f64::NAN, "fission");
    bad(|m| m.states[2].nu_fission[1] = -0.1, "nu_fission");
    bad(|m| m.states[0].total[2] = 0.0, "total");
    // 散乱のモーメントがない、状態で次数が違う
    bad(|m| m.states[0].scatter.clear(), "scatter");
    bad(
        |m| {
            m.states[3].scatter.pop();
        },
        "moments",
    );
    // χ は核分裂があれば和が 1、なければ 0
    bad(|m| m.states[0].chi = vec![0.5, 0.25, 0.0], "chi");
    bad(|m| m.states[0].chi = vec![1.25, -0.25, 0.0], "chi");
    bad(
        |m| {
            m.states[0].nu_fission = vec![0.0; 3];
            m.states[0].fission = vec![0.0; 3];
        },
        "chi",
    );
}

#[test]
fn media_without_fission_have_zero_chi() {
    let mut medium = sample();
    for s in &mut medium.states {
        s.nu_fission = vec![0.0; 3];
        s.fission = vec![0.0; 3];
        s.chi = vec![0.0; 3];
    }
    assert!(medium.validate().is_ok());
}

#[test]
fn reading_is_strict() {
    let bytes = sample().to_msgpack().unwrap();
    // 形式の名前
    let mut renamed = bytes.clone();
    let at = renamed
        .windows(10)
        .position(|w| w == b"xlc-medium")
        .unwrap();
    renamed[at + 4] = b'M';
    assert!(Medium::from_msgpack(&renamed).is_err());
    // 末尾の余り
    let mut longer = bytes.clone();
    longer.push(0xc0);
    assert!(Medium::from_msgpack(&longer).is_err());
    // 状態の項目が欠けている、知らない項目がある
    let edit = |change: fn(&mut Vec<(rmpv::Value, rmpv::Value)>)| {
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
        change(fields);
        let mut out = Vec::new();
        rmpv::encode::write_value(&mut out, &value).unwrap();
        Medium::from_msgpack(&out)
    };
    assert!(edit(|_| {}).is_ok());
    assert!(edit(|f| f.retain(|(k, _)| k.as_str() != Some("chi"))).is_err());
    assert!(edit(|f| f.push(("adf".into(), rmpv::Value::Nil))).is_err());
}

#[test]
fn media_have_no_rod_axis() {
    // 媒質は制御棒の軸を持たない（ADR-0128）
    let mut medium = sample();
    medium.axes.push(Axis::Rod {
        values: vec!["out".into(), "b4c".into()],
    });
    medium.states = (1..=8).map(|i| state(f64::from(i))).collect();
    assert!(medium.validate().is_err());
}

#[test]
fn media_hold_no_nominal_values() {
    // 媒質は軸の直積の状態を持ち、部分格子の公称の状態を使わない（ADR-0200）
    let mut medium = sample();
    medium.axes = medium
        .axes
        .into_iter()
        .map(|a| {
            let value = a.points()[0];
            a.with_nominal(value)
        })
        .collect();
    assert!(medium.validate().is_err());
}
