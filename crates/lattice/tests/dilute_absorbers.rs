//! 輸送しない核種を無限希釈の吸収体として数える燃焼（`[depletion]` の `dilute-absorbers`。ADR-0261）: ライブラリの放射化の核種の
//! 吸収を、燃焼の領域の巨視的断面積、燃焼の反応率、履歴からのブランチ計算の組成、少数群定数に数える。同じ断面積を持ち散乱を
//! 持たない輸送する核種として解いた格子計算と、k、数密度、少数群定数が丸めの程度で一致する（分離共鳴のないライブラリでは、
//! 輸送する核種も無限希釈で、SPH 因子は 1）。
#![allow(
    clippy::unwrap_used,
    clippy::disallowed_macros,
    clippy::disallowed_methods,
    clippy::panic,
    reason = "a test: a failed file or run should stop it with the message"
)]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use xlc_depletion::chain::{Chain, Channel, Fission, Nuclide, Product, Reaction, Yield};
use xlc_format::group::GroupStructure;
use xlc_format::library::{
    ActivationData, ActivationTemperature, DelayedGroups, GroupChannel, GroupReaction, Library,
    NuclideData, Pointwise, Shielded, TemperatureData,
};
use xlc_format::table::Table;
use xlc_lattice::depletion::History;

const T: f64 = 600.0;
/// Cs-133 の捕獲（群ごと、barn）
const CS133_CAPTURE: [f64; 2] = [5.0, 300.0];
/// k、数密度、少数群定数を比べる相対の許容（MOC の判定を 1e-11 にした 2 つの計算の、和の順の違いと止めた解の誤差）
const TOLERANCE: f64 = 1.0e-8;

/// 無限希釈の群定数 `(散乱 [[g→g']], 核分裂, 捕獲)` の、分離共鳴のない輸送する核種。全断面積は散乱、核分裂、捕獲の和。散乱は P1 まで
fn transported(
    name: &str,
    za: f64,
    awr: f64,
    scatter: [[f64; 2]; 2],
    fission: [f64; 2],
    capture: [f64; 2],
) -> NuclideData {
    let total: Vec<f64> = (0..2)
        .map(|g| scatter[g][0] + scatter[g][1] + fission[g] + capture[g])
        .collect();
    let table = Shielded {
        total,
        elastic: (0..2).map(|g| scatter[g][0] + scatter[g][1]).collect(),
        fission: fission.to_vec(),
        nu_fission: fission.iter().map(|f| 2.43 * f).collect(),
        capture: capture.to_vec(),
    };
    NuclideData {
        name: name.into(),
        za,
        awr,
        resolved: None,
        thermal_scattering: false,
        decay_constants: None,
        temperatures: vec![TemperatureData {
            sab_transfer: None,
            delayed: None,
            temperature: T,
            chi: if fission[0] > 0.0 {
                vec![1.0, 0.0]
            } else {
                vec![0.0, 0.0]
            },
            scatter: [1.0, 0.1]
                .iter()
                .map(|f| {
                    scatter
                        .iter()
                        .map(|row| row.iter().map(|s| f * s).collect())
                        .collect()
                })
                .collect(),
            shielded: vec![table.clone(), table.clone(), table],
            pointwise: Pointwise::default(),
            reactions: vec![GroupReaction {
                mt: 102,
                channels: vec![GroupChannel {
                    final_state: None,
                    cross_section: capture.to_vec(),
                }],
            }],
        }],
    }
}

/// 2 群（0.625 eV で分ける）の、分離共鳴のないライブラリ。Cs-133 と Cs-134 は持たない。U-235 は遅発中性子の 2 族を持つ
/// （履歴からのブランチ計算の動特性のデータに要る）
fn base_library() -> Library {
    let mut library = plain_library();
    let u235 = library
        .nuclides
        .iter_mut()
        .find(|n| n.name == "U-235")
        .unwrap();
    u235.decay_constants = Some(vec![0.0133, 0.0327]);
    for t in &mut u235.temperatures {
        t.delayed = Some(DelayedGroups {
            nu: vec![vec![0.004, 0.0041], vec![0.012, 0.0118]],
            chi: vec![vec![1.0, 0.0], vec![0.9, 0.1]],
        });
    }
    library
}

fn plain_library() -> Library {
    Library {
        inverse_velocity: vec![5.0e-8, 2.0e-6],
        activation: Vec::new(),
        groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
        sigma0: vec![1.0e10, 100.0, 1.0],
        pointwise_range: None,
        nuclides: vec![
            transported(
                "H-1",
                1001.0,
                0.99917,
                [[3.9, 0.3], [0.0, 20.0]],
                [0.0; 2],
                [0.0, 0.3],
            ),
            transported(
                "B-10",
                5010.0,
                9.9269,
                [[2.0, 0.0], [0.0, 2.0]],
                [0.0; 2],
                [1.0, 3000.0],
            ),
            transported(
                "B-11",
                5011.0,
                10.9147,
                [[4.0, 0.0], [0.0, 4.0]],
                [0.0; 2],
                [0.0, 0.005],
            ),
            transported(
                "O-16",
                8016.0,
                15.858,
                [[3.7, 0.05], [0.0, 3.8]],
                [0.0; 2],
                [0.0; 2],
            ),
            transported(
                "Xe-135",
                54135.0,
                133.75,
                [[5.0, 0.0], [0.0, 5.0]],
                [0.0; 2],
                [1.0, 2.0e6],
            ),
            transported(
                "U-235",
                92235.0,
                233.025,
                [[4.0, 0.0], [0.0, 15.0]],
                [1.2, 580.0],
                [0.3, 99.0],
            ),
            transported(
                "U-238",
                92238.0,
                236.006,
                [[8.0, 0.0], [0.0, 9.0]],
                [0.0; 2],
                [0.4, 2.7],
            ),
        ],
    }
}

/// Cs-133 を放射化の核種（捕獲だけ）として持つライブラリ
fn activation_library() -> Library {
    let mut library = base_library();
    library.activation.push(ActivationData {
        name: "Cs-133".into(),
        za: 55133.0,
        awr: 131.764,
        temperatures: vec![ActivationTemperature {
            temperature: T,
            fission: vec![0.0, 0.0],
            reactions: vec![GroupReaction {
                mt: 102,
                channels: vec![GroupChannel {
                    final_state: None,
                    cross_section: CS133_CAPTURE.to_vec(),
                }],
            }],
        }],
    });
    library
}

/// Cs-133 を、同じ捕獲を持ち散乱を持たない輸送する核種として持つライブラリ（核種はチェーンの順でなくてよい）
fn transporting_library() -> Library {
    let mut library = base_library();
    library.nuclides.push(transported(
        "Cs-133",
        55133.0,
        131.764,
        [[0.0; 2]; 2],
        [0.0; 2],
        CS133_CAPTURE,
    ));
    library
}

/// 燃料と減速材の核種と、U-235 の核分裂でできる Xe-135（0.1 個）と Cs-133（1.9 個）、Cs-133 の捕獲でできる Cs-134 のチェーン
/// （核種は (z, a) の順）
fn chain() -> Chain {
    let stable = |name: &str, z: u32, a: u32| Nuclide {
        name: name.into(),
        z,
        a,
        state: 0,
        excitation: 0.0,
        decay: None,
        reactions: vec![],
        fission: None,
    };
    let (xe135, cs133, cs134) = (4, 5, 6);
    let mut absorber = stable("Cs-133", 55, 133);
    absorber.reactions.push(Reaction {
        mt: 102,
        channels: vec![Channel {
            final_state: None,
            products: vec![Product {
                nuclide: cs134,
                count: 1,
            }],
        }],
    });
    let mut u235 = stable("U-235", 92, 235);
    let table = vec![
        Yield {
            nuclide: xe135,
            value: 0.1,
        },
        Yield {
            nuclide: cs133,
            value: 1.9,
        },
    ];
    u235.fission = Some(Fission {
        q: 1.93e8,
        yield_source: "U-235".into(),
        energies: vec![0.0253, 5e5],
        yields: vec![table.clone(), table],
    });
    Chain {
        nuclides: vec![
            stable("H-1", 1, 1),
            stable("B-10", 5, 10),
            stable("B-11", 5, 11),
            stable("O-16", 8, 16),
            stable("Xe-135", 54, 135),
            absorber,
            stable("Cs-134", 55, 134),
            u235,
            stable("U-238", 92, 238),
        ],
    }
}

/// ピンセル。MOC の判定は既定より 4 桁細かい（2 つの計算の k を止めた解の誤差より細かく比べる）
const LATTICE: &str = r#"
library = "library.msgpack"
pitch = 1.26
azimuthal = 8
spacing = 0.1
rows = ["pin"]

[moc]
k-tolerance = 1e-11
flux-tolerance = 1e-10
cmfd-k-tolerance = 1e-14
cmfd-flux-tolerance = 1e-12

[material.fuel]
temperature = 600
nuclides = { U-235 = 7.0e-4, U-238 = 2.2e-2, O-16 = 4.6e-2 }

[material.water]
temperature = 600
water = { density = 0.743, boron = 1300, b10 = 0.199 }

[cell.pin]
zones = [{ material = "fuel", radius = 0.4096, rings = 2 }, { material = "water", rings = 2 }]
"#;

const DEPLETION: &str = r#"
[depletion]
chain = "chain.msgpack"
power = 40.0
burnup = [0.5, 2.0]
materials = ["fuel"]
output = "history.msgpack"
few-group = { bounds = [1e7, 0.625, 1e-5] }
"#;

const DILUTE: &str = "dilute-absorbers = true\n";

/// 履歴からのブランチ計算（燃焼点 0、0.5、2.0。ホウ素のブランチ条件と、Xe-135 の数密度の微分）
const BRANCHES: &str = r#"
[branches]
boron = { ppm = [600, 1300], materials = ["water"] }
moderator = { temperature = [600], density = [0.743], materials = ["water"] }
fuel = { temperature = [600], materials = ["fuel"] }

[history]
file = "history.msgpack"
chain = "chain.msgpack"
burnup = [0, 0.5, 2.0]
micro = ["U-235", "Xe-135"]
derivative = ["Xe-135"]
"#;

fn branched() -> String {
    let library = "library = \"library.msgpack\"\n";
    assert!(LATTICE.contains(library));
    LATTICE.replace(
        library,
        &format!(
            "{library}few-group = {{ output = \"fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }}\n"
        ),
    )
}

fn scratch(name: &str) -> PathBuf {
    // cargo がテストのために用意する target の下の一時の置き場（失敗して残った置き場も cargo clean で消える）
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("dilute-absorbers-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// `dir` に `text` の入力、ライブラリ `library`、チェーンを置き、1 つのスレッドで格子計算を実行する
fn execute(dir: &Path, library: &Library, text: &str) -> Output {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(dir.join("library.msgpack"), library.to_msgpack().unwrap()).unwrap();
    std::fs::write(dir.join("chain.msgpack"), chain().to_msgpack().unwrap()).unwrap();
    std::fs::write(dir.join("in.toml"), text).unwrap();
    Command::new(env!("CARGO_BIN_EXE_lattice"))
        .arg(dir.join("in.toml"))
        .env("RAYON_NUM_THREADS", "1")
        .output()
        .unwrap()
}

fn run(dir: &Path, library: &Library, text: &str) {
    let output = execute(dir, library, text);
    assert!(
        output.status.success(),
        "lattice failed in {}:\n{}",
        dir.display(),
        String::from_utf8_lossy(&output.stderr)
    );
}

/// 失敗する格子計算の標準エラー
fn failure(dir: &Path, library: &Library, text: &str) -> String {
    let output = execute(dir, library, text);
    assert!(
        !output.status.success(),
        "lattice succeeded in {}",
        dir.display()
    );
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// `dir` で燃焼させた履歴（ファイルの中身と読んだもの）
fn deplete(dir: &Path, library: &Library, extra: &str) -> (Vec<u8>, History) {
    run(dir, library, &format!("{LATTICE}{DEPLETION}{extra}"));
    let bytes = std::fs::read(dir.join("history.msgpack")).unwrap();
    let history = History::from_msgpack(&bytes).unwrap();
    (bytes, history)
}

/// 履歴 `history` からのブランチ計算のテーブル
fn branch(dir: &Path, library: &Library, history: &[u8]) -> Table {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(dir.join("history.msgpack"), history).unwrap();
    run(dir, library, &format!("{}{BRANCHES}", branched()));
    Table::from_msgpack(&std::fs::read(dir.join("fg.msgpack")).unwrap()).unwrap()
}

fn assert_close(got: f64, want: f64, at: &str) {
    assert!(
        (got - want).abs() <= TOLERANCE * want.abs().max(1.0e-30),
        "{at}: {got} vs {want}"
    );
}

fn assert_all_close(got: &[f64], want: &[f64], at: &str) {
    assert_eq!(got.len(), want.len(), "{at}");
    for (g, w) in got.iter().zip(want) {
        assert_close(*g, *w, at);
    }
}

fn density(history: &History, point: usize, name: &str) -> f64 {
    let i = history.nuclides.iter().position(|n| n == name).unwrap();
    history.points[point].densities[0][i]
}

/// 履歴のファイルの最上位の map のキー
fn keys(history: &[u8]) -> Vec<String> {
    let value = xlc_format::msgpack::read_value(history).unwrap();
    value
        .as_map()
        .unwrap()
        .iter()
        .map(|(k, _)| k.as_str().unwrap().to_owned())
        .collect()
}

#[test]
fn dilute_absorbers_deplete_like_transported_nuclides_without_scattering() {
    let dir = scratch("depletion");
    let (bytes, dilute) = deplete(&dir.join("dilute"), &activation_library(), DILUTE);
    let (_, transported) = deplete(&dir.join("transported"), &transporting_library(), "");
    // 輸送する核種は増えず、無限希釈の吸収体は履歴に名前で残る
    assert!(!dilute.transportable.contains(&"Cs-133".to_owned()));
    assert!(transported.transportable.contains(&"Cs-133".to_owned()));
    assert!(keys(&bytes).contains(&"dilute_absorbers".to_owned()));
    assert_eq!(dilute.points.len(), 3);
    for (p, (d, t)) in dilute.points.iter().zip(&transported.points).enumerate() {
        let at = format!("point {p}");
        assert_close(d.k, t.k, &at);
        // 燃焼の反応率: Cs-133 は捕獲で減り、Cs-134 ができる
        assert_all_close(&d.densities[0], &t.densities[0], &at);
        // 少数群定数: 吸収は全断面積と吸収に入り、散乱と核分裂には入らない
        let (d, t) = (d.few_group.as_ref().unwrap(), t.few_group.as_ref().unwrap());
        assert_close(d.k_inf, t.k_inf, &at);
        assert_all_close(&d.total, &t.total, &at);
        assert_all_close(&d.absorption, &t.absorption, &at);
        assert_all_close(&d.nu_fission, &t.nu_fission, &at);
        assert_all_close(&d.diffusion, &t.diffusion, &at);
        for (d, t) in d.scatter.iter().zip(&t.scatter) {
            assert_all_close(d, t, &at);
        }
    }
    assert!(density(&dilute, 2, "Cs-134") > 0.0);
}

#[test]
fn without_the_option_activation_nuclides_only_decay() {
    let dir = scratch("option");
    let (bytes, dilute) = deplete(&dir.join("dilute"), &activation_library(), DILUTE);
    let (plain_bytes, plain) = deplete(&dir.join("plain"), &activation_library(), "");
    // 選ばない燃焼の履歴はキーを持たず、今の形のまま
    assert!(!keys(&plain_bytes).contains(&"dilute_absorbers".to_owned()));
    assert!(keys(&bytes).contains(&"dilute_absorbers".to_owned()));
    // 燃焼度 0 は吸収体が無く、同じ輸送
    assert_eq!(dilute.points[0].k, plain.points[0].k);
    for p in 1..3 {
        assert!(
            dilute.points[p].k < plain.points[p].k,
            "point {p}: {} vs {}",
            dilute.points[p].k,
            plain.points[p].k
        );
        assert!(density(&dilute, p, "Cs-133") < density(&plain, p, "Cs-133"));
    }
    assert_eq!(density(&plain, 2, "Cs-134"), 0.0);
}

#[test]
fn branches_of_a_history_with_dilute_absorbers_count_them_in_every_state() {
    let dir = scratch("branches");
    let (bytes, dilute) = deplete(&dir.join("dilute"), &activation_library(), DILUTE);
    let (transported_bytes, _) = deplete(&dir.join("transported"), &transporting_library(), "");
    let table = branch(&dir.join("dilute-branches"), &activation_library(), &bytes);
    let reference = branch(
        &dir.join("transported-branches"),
        &transporting_library(),
        &transported_bytes,
    );
    assert_eq!(table.points, reference.points);
    assert_eq!(table.states.len(), 3 * 2);
    let micro = table.micro.as_ref().unwrap();
    let reference_micro = reference.micro.as_ref().unwrap();
    for (s, (state, want)) in table.states.iter().zip(&reference.states).enumerate() {
        let at = format!("state {s} at {:?}", table.points[s]);
        assert_close(state.k_inf, want.k_inf, &at);
        assert_all_close(&state.total, &want.total, &at);
        assert_all_close(&state.absorption, &want.absorption, &at);
        // 基準条件（1300 ppm）の状態は、履歴の点の輸送
        if table.points[s][1] == 1 {
            assert_close(state.k_inf, dilute.points[table.points[s][0]].k, &at);
        }
        // 追う核種の微視的断面積と、摂動の組成を解き直す数密度の微分
        for (m, w) in micro.states[s].iter().zip(&reference_micro.states[s]) {
            assert_all_close(&m.total, &w.total, &at);
            match (&m.derivative, &w.derivative) {
                (Some(m), Some(w)) => {
                    for (m, w) in m.total.iter().zip(&w.total) {
                        // 前進差分は k の差を数密度の変化で割るので、止めた解の誤差が大きく入る
                        assert!((m - w).abs() <= 1.0e-4 * w.abs(), "{at}: {m} vs {w}");
                    }
                }
                (None, None) => {}
                _ => panic!("{at}: the derivatives differ"),
            }
        }
    }
}

#[test]
fn dilute_absorbers_need_activation_data_at_the_depletion_temperature() {
    let dir = scratch("errors");
    // 放射化の核種のないライブラリ
    let err = failure(
        &dir.join("none"),
        &base_library(),
        &format!("{LATTICE}{DEPLETION}{DILUTE}"),
    );
    assert!(err.contains("dilute-absorbers"), "{err}");
    // 履歴が数えた吸収体を、ブランチ計算のライブラリが持たない
    let (bytes, _) = deplete(&dir.join("dilute"), &activation_library(), DILUTE);
    let at = dir.join("branches");
    std::fs::create_dir_all(&at).unwrap();
    std::fs::write(at.join("history.msgpack"), &bytes).unwrap();
    let err = failure(&at, &base_library(), &format!("{}{BRANCHES}", branched()));
    assert!(err.contains("dilute absorbers"), "{err}");
    // 一般化摂動の微分は、重ねた物質で確かめていない（ADR-0261）
    let generalized = format!(
        "{}{BRANCHES}generalized = {{ tolerance = 1e-7, restart = 100, shielded = [\"U-235\"] }}\n",
        branched()
    );
    let err = failure(&at, &activation_library(), &generalized);
    assert!(err.contains("generalized derivatives"), "{err}");
}
