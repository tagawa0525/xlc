//! 読むだけの入口は、出力の形式の名前で本番の読み手を選び、その読み手が読めない出力をエラーにする（ADR-0258）。
#![allow(
    clippy::unwrap_used,
    clippy::disallowed_macros,
    clippy::panic,
    reason = "a test: a failed read should stop it with the message"
)]

use rmpv::Value;
use xlc_format::medium::{Medium, MediumState};
use xlc_outputs::{Kind, read};

/// 2 群、軸のない 1 状態の媒質を、本番の書き手で書いたもの
fn medium() -> Vec<u8> {
    Medium {
        bounds: vec![2.0e7, 0.625, 1.0e-5],
        axes: vec![],
        states: vec![MediumState {
            total: vec![0.5, 1.0],
            scatter: vec![vec![vec![0.4, 0.05], vec![0.0, 0.9]]],
            fission: vec![0.0, 0.02],
            nu_fission: vec![0.0, 0.05],
            chi: vec![1.0, 0.0],
        }],
    }
    .to_msgpack()
    .unwrap()
}

/// 最上位の map を `edit` で書き換えて書き直す
fn edited(bytes: &[u8], edit: impl FnOnce(&mut Vec<(Value, Value)>)) -> Vec<u8> {
    let mut value = rmpv::decode::read_value(&mut &bytes[..]).unwrap();
    let Value::Map(entries) = &mut value else {
        panic!("not a map")
    };
    edit(entries);
    let mut out = Vec::new();
    rmpv::encode::write_value(&mut out, &value).unwrap();
    out
}

fn named(format: &str) -> Vec<u8> {
    let value = Value::Map(vec![(Value::from("format"), Value::from(format))]);
    let mut out = Vec::new();
    rmpv::encode::write_value(&mut out, &value).unwrap();
    out
}

#[test]
fn reads_an_output_the_writer_wrote() {
    read(Kind::Msgpack, &medium()).unwrap();
}

#[test]
fn an_output_of_a_renamed_format_is_unreadable() {
    let bytes = edited(&medium(), |entries| {
        for (k, v) in entries.iter_mut() {
            if k.as_str() == Some("format") {
                *v = Value::from("xlc-medium-0");
            }
        }
    });
    let err = read(Kind::Msgpack, &bytes).unwrap_err();
    assert!(format!("{err:#}").contains("xlc-medium-0"), "{err:#}");
}

#[test]
fn an_output_without_a_required_field_is_unreadable() {
    let bytes = edited(&medium(), |entries| {
        entries.retain(|(k, _)| k.as_str() != Some("states"));
    });
    let err = read(Kind::Msgpack, &bytes).unwrap_err();
    assert!(format!("{err:#}").contains("states"), "{err:#}");
}

#[test]
fn an_output_with_an_unknown_field_is_unreadable() {
    let bytes = edited(&medium(), |entries| {
        entries.push((Value::from("extra"), Value::from(1)));
    });
    read(Kind::Msgpack, &bytes).unwrap_err();
}

#[test]
fn a_format_without_a_rust_reader_is_unreadable() {
    let err = read(Kind::Msgpack, &named("xlc-pin-powers")).unwrap_err();
    assert!(format!("{err:#}").contains("xlc-pin-powers"), "{err:#}");
}

#[test]
fn a_file_without_a_format_name_is_unreadable() {
    read(Kind::Msgpack, b"pointwise-range = []\n").unwrap_err();
    let mut map = Vec::new();
    rmpv::encode::write_value(&mut map, &Value::Map(vec![])).unwrap();
    read(Kind::Msgpack, &map).unwrap_err();
}

/// 名前だけのファイルは、どの形式でも項目が足りずに止まる。止めたのがその形式の読み手であること（読み手が無い、でない）を見る
#[test]
fn every_format_with_a_rust_reader_goes_to_its_reader() {
    let formats = [
        xlc_format::library::FORMAT,
        xlc_format::table::FORMAT,
        xlc_format::medium::FORMAT,
        xlc_depletion::chain::FORMAT,
        xlc_lattice::depletion::FORMAT,
        xlc_core::depletion::FORMAT,
    ];
    let without_reader = format!("{:#}", read(Kind::Msgpack, &named("xlc-none")).unwrap_err());
    for format in formats {
        let err = format!("{:#}", read(Kind::Msgpack, &named(format)).unwrap_err());
        assert_ne!(err, without_reader.replace("xlc-none", format), "{format}");
    }
}

#[test]
fn reads_the_range_file_with_the_xs_reader() {
    read(
        Kind::XsRange,
        xlc_xs::range::to_toml(Some((4.0, 9118.0))).as_bytes(),
    )
    .unwrap();
    read(Kind::XsRange, xlc_xs::range::to_toml(None).as_bytes()).unwrap();
    read(Kind::XsRange, b"pointwise-range = [5.0, 4.0]\n").unwrap_err();
    read(Kind::XsRange, b"range = []\n").unwrap_err();
}
