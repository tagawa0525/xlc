//! 制御棒の分岐（`[branches.rod]`。ADR-0224）: 1 つの入力で解いた制御棒の軸を持つテーブルの各状態が、基準の形状のセルを
//! 置き換えた別の入力の格子計算の状態とビットで一致する（値ごとの別の格子計算をまとめた形と同じ問題を解く）。細かい群の媒質は
//! 棒のない値の状態だけを書き、媒質を読む格子は、どの値でも媒質の同じ連続の軸の状態を使う（ADR-0128）。
#![allow(
    clippy::unwrap_used,
    clippy::disallowed_macros,
    clippy::disallowed_methods,
    clippy::panic,
    reason = "a test: a failed file or run should stop it with the message"
)]

use std::path::{Path, PathBuf};
use std::process::Command;

use xlc_format::group::GroupStructure;
use xlc_format::library::{
    GroupChannel, GroupReaction, Library, NuclideData, Pointwise, Shielded, TemperatureData,
};
use xlc_format::medium::Medium;
use xlc_format::table::{Axis, Table};

const T: f64 = 600.0;

/// 2 群（0.625 eV で分ける）の、分離共鳴のないライブラリ。全断面積は散乱、核分裂、捕獲の和。散乱は P1 まで
fn library() -> Library {
    // (名前, ZA, awr, 散乱 [[g→g']], 核分裂, 捕獲)
    let nuclide = |name: &str,
                   za: f64,
                   awr: f64,
                   scatter: [[f64; 2]; 2],
                   fission: [f64; 2],
                   capture: [f64; 2]| {
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
    };
    Library {
        inverse_velocity: vec![1.0e-6; 2],
        activation: Vec::new(),
        groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
        sigma0: vec![1.0e10, 100.0, 1.0],
        pointwise_range: None,
        nuclides: vec![
            nuclide(
                "H-1",
                1001.0,
                0.99917,
                [[3.9, 0.3], [0.0, 20.0]],
                [0.0; 2],
                [0.0, 0.3],
            ),
            nuclide(
                "O-16",
                8016.0,
                15.858,
                [[3.7, 0.05], [0.0, 3.8]],
                [0.0; 2],
                [0.0; 2],
            ),
            nuclide(
                "B-10",
                5010.0,
                9.9269,
                [[2.0, 0.0], [0.0, 2.0]],
                [0.0; 2],
                [1.0, 3000.0],
            ),
            nuclide(
                "B-11",
                5011.0,
                10.9147,
                [[4.0, 0.0], [0.0, 4.0]],
                [0.0; 2],
                [0.0, 0.005],
            ),
            nuclide(
                "Ag-109",
                47109.0,
                107.969,
                [[5.0, 0.01], [0.0, 5.0]],
                [0.0; 2],
                [2.0, 90.0],
            ),
            nuclide(
                "U-235",
                92235.0,
                233.025,
                [[4.0, 0.0], [0.0, 15.0]],
                [1.2, 580.0],
                [0.3, 99.0],
            ),
            nuclide(
                "U-238",
                92238.0,
                236.006,
                [[8.0, 0.0], [0.0, 9.0]],
                [0.3, 0.0],
                [0.4, 2.7],
            ),
        ],
    }
}

/// 2 × 2 の格子（燃料のピンと案内管）。案内管の中は水で、制御棒の値の分岐では吸収体の棒を入れたセルに置き換える
const ASSEMBLY: &str = r#"
library = "library.msgpack"
self-shielding = false
few-group = { output = "fg.msgpack", bounds = [1e7, 0.625, 1e-5] }
medium = "medium.msgpack"
pitch = 1.26
azimuthal = 8
spacing = 0.1
rows = ["pin gt", "gt pin"]

[material.fuel]
temperature = 600
nuclides = { U-235 = 7.0e-4, U-238 = 2.2e-2, O-16 = 4.6e-2 }

[material.absorber]
temperature = 600
nuclides = { Ag-109 = 4.0e-2 }

[material.strong]
temperature = 600
nuclides = { Ag-109 = 8.0e-2, B-10 = 1.0e-2 }

[material.water]
temperature = 600
water = { density = 0.743, boron = 1300, b10 = 0.199 }

[cell.pin]
zones = [{ material = "fuel", radius = 0.4096, rings = 2 }, { material = "water", rings = 2 }]

[cell.gt]
zones = [{ material = "water", radius = 0.56, rings = 2 }, { material = "water", rings = 1 }]

[cell.gt-a]
zones = [{ material = "absorber", radius = 0.38, rings = 2 }, { material = "void", radius = 0.39, rings = 1 },
         { material = "water", radius = 0.56, rings = 1 }, { material = "water", rings = 1 }]

[cell.gt-b]
zones = [{ material = "strong", radius = 0.37, rings = 3 }, { material = "water", radius = 0.56, rings = 1 },
         { material = "water", rings = 1 }]

[branches]
boron = { ppm = [600, 1300], materials = ["water"] }
moderator = { temperature = [600], density = [0.743], materials = ["water"] }
fuel = { temperature = [600], materials = ["fuel"] }
"#;

const RODS: &str = "\n[branches.rod]\nb = { gt = \"gt-b\" }\na = { gt = \"gt-a\" }\n";

fn scratch(name: &str) -> PathBuf {
    // cargo がテストのために用意する target の下の一時の置き場（失敗して残った置き場も cargo clean で消える）
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("rod-branch-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// `dir` に `text` の入力とライブラリを置き、1 つのスレッドで格子計算を実行する（MOC の和の順をそろえる）
fn run(dir: &Path, text: &str) {
    std::fs::write(dir.join("library.msgpack"), library().to_msgpack().unwrap()).unwrap();
    std::fs::write(dir.join("in.toml"), text).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_lattice"))
        .arg(dir.join("in.toml"))
        .env("RAYON_NUM_THREADS", "1")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "lattice failed in {}:\n{}",
        dir.display(),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn table(dir: &Path, file: &str) -> Table {
    Table::from_msgpack(&std::fs::read(dir.join(file)).unwrap()).unwrap()
}

/// 基準の形状の `rows` の gt を `cell` に置き換えた、制御棒の分岐のない入力
fn replaced(text: &str, cell: &str) -> String {
    let from = r#"rows = ["pin gt", "gt pin"]"#;
    assert!(text.contains(from));
    text.replace(from, &format!(r#"rows = ["pin {cell}", "{cell} pin"]"#))
}

#[test]
fn rod_branches_match_separate_lattices() {
    let dir = scratch("assembly");
    let combined = dir.join("combined");
    std::fs::create_dir_all(&combined).unwrap();
    run(&combined, &format!("{ASSEMBLY}{RODS}"));
    let rods = table(&combined, "fg.msgpack");
    let values = ["out", "a", "b"];
    assert_eq!(
        rods.axes.last(),
        Some(&Axis::Rod {
            values: values.iter().map(|v| v.to_string()).collect()
        })
    );
    for (r, cell) in ["gt", "gt-a", "gt-b"].iter().enumerate() {
        let separate = dir.join(values[r]);
        std::fs::create_dir_all(&separate).unwrap();
        run(&separate, &replaced(ASSEMBLY, cell));
        let one = table(&separate, "fg.msgpack");
        assert_eq!(rods.axes[..3], one.axes[..], "the continuous axes");
        for (point, state) in one.points.iter().zip(&one.states) {
            let mut at = point.clone();
            at.push(r);
            let n = rods.points.iter().position(|p| *p == at).unwrap();
            assert_eq!(
                &rods.states[n], state,
                "rod value {} at {point:?}",
                values[r]
            );
        }
    }
    assert_eq!(rods.points.len(), 2 * values.len());
    // 媒質は棒のない値の状態だけ（連続の軸の状態）
    let medium = std::fs::read(combined.join("medium.msgpack")).unwrap();
    let out = std::fs::read(dir.join("out").join("medium.msgpack")).unwrap();
    assert_eq!(
        Medium::from_msgpack(&medium).unwrap(),
        Medium::from_msgpack(&out).unwrap()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

/// 媒質を読む 1 行の平板（ADR-0128 の軸方向の平板の形）。右の層のセルを、制御棒の値で吸収体を混ぜた層に置き換える
const SLAB: &str = r#"
library = "library.msgpack"
self-shielding = false
few-group = { output = "fg.msgpack", bounds = [1e7, 0.625, 1e-5] }
pitch = 1.26
azimuthal = 8
spacing = 0.1
rows = ["fuel fuel layer layer"]

[material.medium]
medium = "assembly-medium.msgpack"

[material.absorber]
temperature = 600
nuclides = { Ag-109 = 4.0e-2 }

[material.water]
temperature = 600
water = { density = 0.743, boron = 1300, b10 = 0.199 }

[material.layer-a]
mix = [{ material = "water", volume = 0.9 }, { material = "absorber", volume = 0.1 }]

[cell.fuel]
zones = [{ material = "medium", rings = 1 }]

[cell.layer]
zones = [{ material = "water", rings = 1 }]

[cell.layer-a]
zones = [{ material = "layer-a", rings = 1 }]

[branches]
boron = { ppm = [600, 1300], materials = ["water", "layer-a"] }
moderator = { temperature = [600], density = [0.743], materials = ["water", "layer-a"] }
fuel = { temperature = [600], materials = ["medium"] }
"#;

#[test]
fn rod_branches_of_a_slab_use_the_same_medium_state_for_every_value() {
    let dir = scratch("slab");
    // 媒質は分岐の 2 つの状態を持つ（制御棒の分岐のない格子で作る）
    let assembly = dir.join("assembly");
    std::fs::create_dir_all(&assembly).unwrap();
    run(&assembly, ASSEMBLY);
    let medium = std::fs::read(assembly.join("medium.msgpack")).unwrap();
    let combined = dir.join("combined");
    let separate = dir.join("a");
    for d in [&combined, &separate] {
        std::fs::create_dir_all(d).unwrap();
        std::fs::write(d.join("assembly-medium.msgpack"), &medium).unwrap();
    }
    run(
        &combined,
        &format!("{SLAB}\n[branches.rod]\na = {{ layer = \"layer-a\" }}\n"),
    );
    let rods = table(&combined, "fg.msgpack");
    run(
        &separate,
        &SLAB.replace(
            r#"rows = ["fuel fuel layer layer"]"#,
            r#"rows = ["fuel fuel layer-a layer-a"]"#,
        ),
    );
    let one = table(&separate, "fg.msgpack");
    for (point, state) in one.points.iter().zip(&one.states) {
        let mut at = point.clone();
        at.push(1);
        let n = rods.points.iter().position(|p| *p == at).unwrap();
        assert_eq!(&rods.states[n], state, "rod value a at {point:?}");
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// 径方向の反射体の格子の形（ADR-0095）: 2 × 2 の燃料のピンと案内管の右に水の反射体の列を置き、右の辺は真空。燃料のノードと
/// 反射体のノードの定数を書く
fn reflector_lattice(text: &str) -> String {
    let from = r#"rows = ["pin gt", "gt pin"]"#;
    let output = r#"few-group = { output = "fg.msgpack", bounds = [1e7, 0.625, 1e-5] }"#;
    assert!(text.contains(from) && text.contains(output));
    text.replace(from, r#"rows = ["pin gt w", "gt pin w"]"#)
        .replace(
            output,
            "few-group = { bounds = [1e7, 0.625, 1e-5], nodes = [\n  \
             { output = \"fuel.msgpack\", columns = [1, 2], rows = [1, 2] },\n  \
             { output = \"reflector.msgpack\", columns = [3, 3], rows = [1, 2], reflector = true }] }\n\
             boundary = { right = \"vacuum\" }",
        )
        .replace("medium = \"medium.msgpack\"\n", "")
        .replace(
            "\n[branches]\n",
            "\n[cell.w]\nzones = [{ material = \"water\", rings = 1 }]\n\n[branches]\n",
        )
}

#[test]
fn nodes_without_replaced_cells_hold_the_rodless_states() {
    // 案内管は燃料のノードにだけあるので、燃料のノードのテーブルは制御棒の軸を持ち、反射体のノードのテーブルは棒のない値の
    // 状態だけを持つ（ADR-0227）
    let dir = scratch("reflector");
    let combined = dir.join("combined");
    std::fs::create_dir_all(&combined).unwrap();
    let base = reflector_lattice(ASSEMBLY);
    run(&combined, &format!("{base}{RODS}"));
    let fuel = table(&combined, "fuel.msgpack");
    let reflector = table(&combined, "reflector.msgpack");
    let values = ["out", "a", "b"];
    assert_eq!(
        fuel.axes.last(),
        Some(&Axis::Rod {
            values: values.iter().map(|v| v.to_string()).collect()
        })
    );
    assert_eq!(fuel.points.len(), 2 * values.len());
    for (r, cell) in ["gt", "gt-a", "gt-b"].iter().enumerate() {
        let separate = dir.join(values[r]);
        std::fs::create_dir_all(&separate).unwrap();
        let from = r#"rows = ["pin gt w", "gt pin w"]"#;
        assert!(base.contains(from));
        run(
            &separate,
            &base.replace(from, &format!(r#"rows = ["pin {cell} w", "{cell} pin w"]"#)),
        );
        let one = table(&separate, "fuel.msgpack");
        assert_eq!(fuel.axes[..3], one.axes[..], "the continuous axes");
        for (point, state) in one.points.iter().zip(&one.states) {
            let mut at = point.clone();
            at.push(r);
            let n = fuel.points.iter().position(|p| *p == at).unwrap();
            assert_eq!(
                &fuel.states[n], state,
                "rod value {} at {point:?}",
                values[r]
            );
        }
    }
    // 反射体のノードのテーブルは、基準の形状の別の入力のものと同じ
    assert_eq!(reflector, table(&dir.join("out"), "reflector.msgpack"));
    assert!(!reflector.axes.iter().any(|a| matches!(a, Axis::Rod { .. })));
    std::fs::remove_dir_all(&dir).unwrap();
}
