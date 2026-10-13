//! VERA 1B のピンセルで、飛程を表から内挿した格子計算（ADR-0156 の既定の細かさ）が、格子点ごとに求めた計算と
//! k と吸収断面積で ADR-0156 の表の差の中に収まること。
//!
//! examples/vera/vera.msgpack（DVC のステージの出力で git で追跡しない。ADR-0130）を読むので、既定では実行しない。
//! 実行方法: `cargo test --release -p xlc-lattice --test track_table -- --ignored`
#![allow(
    clippy::unwrap_used,
    clippy::disallowed_macros,
    reason = "a test: a failed file or run should stop it with the message"
)]

use std::path::{Path, PathBuf};
use std::process::Command;

use xlc_format::medium::Medium;
use xlc_format::table::Table;

/// ADR-0156 の表の 1B、10 倍あたり 16 点の差の上限
const K_LIMIT: f64 = 1.0e-6;
const ABSORPTION_LIMIT: f64 = 2.3e-6;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// 1B の入力に、ライブラリの絶対パス、少数群定数と媒質の出力、飛程の求め方 `tracks` を足して実行し、
/// (k, 群ごとの吸収断面積) を返す
fn run(dir: &Path, tracks: &str) -> (f64, Vec<f64>) {
    let base = std::fs::read_to_string(root().join("examples/vera/p1/pincell-1B.toml")).unwrap();
    let library = root().join("examples/vera/vera.msgpack");
    let text = base.replace(
        "library = \"../vera.msgpack\"",
        &format!(
            "library = \"{}\"\n{tracks}\nmedium = \"medium.msgpack\"\n\
             few-group = {{ output = \"fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }}",
            library.display()
        ),
    );
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(dir.join("in.toml"), text).unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_lattice"))
        .arg(dir.join("in.toml"))
        .status()
        .unwrap();
    assert!(status.success());
    let few = Table::from_msgpack(&std::fs::read(dir.join("fg.msgpack")).unwrap()).unwrap();
    let medium = Medium::from_msgpack(&std::fs::read(dir.join("medium.msgpack")).unwrap()).unwrap();
    let state = &medium.states[0];
    let absorption = state
        .total
        .iter()
        .zip(&state.scatter[0])
        .map(|(t, row)| t - row.iter().sum::<f64>())
        .collect();
    (few.states[0].k_inf, absorption)
}

#[test]
#[ignore = "needs examples/vera/vera.msgpack (a DVC output, ADR-0130)"]
fn table_tracks_stay_within_the_adr_differences_in_1b() {
    let dir = std::env::temp_dir().join(format!("xlc-track-table-{}", std::process::id()));
    let (k, absorption) = run(&dir.join("table"), "track-table = 16");
    let (k_exact, exact) = run(&dir.join("exact"), "track-table = false");
    std::fs::remove_dir_all(&dir).unwrap();
    assert!((k - k_exact).abs() < K_LIMIT, "k {k} vs {k_exact}");
    let worst = absorption
        .iter()
        .zip(&exact)
        .filter(|(_, e)| **e != 0.0)
        .map(|(a, e)| (a / e - 1.0).abs())
        .fold(0.0, f64::max);
    assert!(worst <= ABSORPTION_LIMIT, "absorption {worst}");
}
