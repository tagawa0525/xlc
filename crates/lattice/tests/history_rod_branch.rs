//! 燃焼の履歴からの分岐（`[history]`。ADR-0187）と制御棒の分岐（`[branches.rod]`。ADR-0224）を組んだ 1 つの入力（ADR-0225）:
//! テーブルの各状態（少数群定数、追う核種の微視的断面積と数密度の微分）が、基準の形状のセルを置き換えた別の入力で同じ履歴を
//! 読んだ格子計算の状態とビットで一致する（値ごとの別の格子計算をまとめた形と同じ問題を解く）。制御棒の値の状態は、その値の
//! 燃焼点（`rod-burnup`）だけで解く。
#![allow(
    clippy::unwrap_used,
    clippy::disallowed_macros,
    clippy::disallowed_methods,
    clippy::panic,
    reason = "a test: a failed file or run should stop it with the message"
)]

use std::path::{Path, PathBuf};
use std::process::Command;

use xlc_depletion::chain::{Chain, Fission, Nuclide, Yield};
use xlc_format::group::GroupStructure;
use xlc_format::library::{
    GroupChannel, GroupReaction, Library, NuclideData, Pointwise, Shielded, TemperatureData,
};
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
            temperatures: vec![TemperatureData {
                sab_transfer: None,
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
                "O-16",
                8016.0,
                15.858,
                [[3.7, 0.05], [0.0, 3.8]],
                [0.0; 2],
                [0.0; 2],
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
                "Xe-135",
                54135.0,
                133.75,
                [[5.0, 0.0], [0.0, 5.0]],
                [0.0; 2],
                [1.0, 2.0e6],
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
                [0.0; 2],
                [0.4, 2.7],
            ),
        ],
    }
}

/// 燃料と減速材の核種と、U-235 の核分裂で 2 個できる Xe-135 のチェーン（核種は (z, a) の順）
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
    let xe135 = 4;
    let mut u235 = stable("U-235", 92, 235);
    let table = vec![Yield {
        nuclide: xe135,
        value: 2.0,
    }];
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
            u235,
            stable("U-238", 92, 238),
        ],
    }
}

/// 2 × 2 の格子（燃料のピンと案内管）。案内管の中は水で、制御棒の値では吸収体の棒を入れたセルに置き換える
const LATTICE: &str = r#"
library = "library.msgpack"
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
"#;

const DEPLETION: &str = r#"
[depletion]
chain = "chain.msgpack"
power = 40.0
burnup = [0.5, 2.0]
materials = ["fuel"]
output = "history.msgpack"
"#;

/// 履歴からの分岐の、制御棒の値のない部分（燃焼点 0、0.5、2.0。公称でない状態は 0 と 2.0）。格子は `branched` で少数群定数の
/// ファイルを足したもの
const BRANCHES: &str = r#"
[branches]
boron = { ppm = [600, 1300], materials = ["water"] }
moderator = { temperature = [600], density = [0.743], materials = ["water"] }
fuel = { temperature = [600], materials = ["fuel"] }

[history]
file = "history.msgpack"
burnup = [0, 0.5, 2.0]
off-nominal-burnup = [0, 2.0]
micro = ["U-235", "Xe-135"]
derivative = ["Xe-135"]
"#;

/// 格子 `LATTICE` に、格子全体の少数群定数のファイルを足したもの
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
        .join(format!("history-rod-branch-{}-{name}", std::process::id()));
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

fn table(dir: &Path) -> Table {
    Table::from_msgpack(&std::fs::read(dir.join("fg.msgpack")).unwrap()).unwrap()
}

/// 履歴の燃焼の入力と同じ格子の、`dir` の格子計算（`history` の燃焼の履歴を置く）
fn branch_run(dir: &Path, history: &[u8], text: &str) -> Table {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(dir.join("history.msgpack"), history).unwrap();
    run(dir, text);
    table(dir)
}

fn burnups(table: &Table) -> Vec<f64> {
    match &table.axes[0] {
        Axis::Burnup { points } => points.clone(),
        other => panic!("the first axis is not the burnup: {other:?}"),
    }
}

#[test]
fn history_rod_branches_match_separate_lattices() {
    let dir = scratch("assembly");
    let depletion = dir.join("depletion");
    std::fs::create_dir_all(&depletion).unwrap();
    std::fs::write(
        depletion.join("chain.msgpack"),
        chain().to_msgpack().unwrap(),
    )
    .unwrap();
    run(&depletion, &format!("{LATTICE}{DEPLETION}"));
    let history = std::fs::read(depletion.join("history.msgpack")).unwrap();

    let rods = branch_run(
        &dir.join("combined"),
        &history,
        &format!(
            "{}{BRANCHES}rod-burnup = {{ a = [0, 2.0] }}\n\n[branches.rod]\na = {{ gt = \"gt-a\" }}\n",
            branched()
        ),
    );
    assert_eq!(
        rods.axes.last(),
        Some(&Axis::Rod {
            values: vec!["out".into(), "a".into()]
        })
    );
    let out = branch_run(
        &dir.join("out"),
        &history,
        &format!("{}{BRANCHES}", branched()),
    );
    // 値 a は、案内管を置き換えた格子で、その値の燃焼点だけを読む別の入力
    let separate = format!("{}{BRANCHES}", branched())
        .replace(
            r#"rows = ["pin gt", "gt pin"]"#,
            r#"rows = ["pin gt-a", "gt-a pin"]"#,
        )
        .replace(
            "burnup = [0, 0.5, 2.0]\noff-nominal-burnup = [0, 2.0]",
            "burnup = [0, 2.0]",
        );
    assert!(separate.contains("gt-a pin") && separate.contains("burnup = [0, 2.0]\nmicro"));
    let a = branch_run(&dir.join("a"), &history, &separate);

    let all = burnups(&out);
    assert_eq!(burnups(&rods), all);
    assert_eq!(rods.axes[..rods.axes.len() - 1], out.axes[..]);
    assert_eq!(rods.heavy_metal, out.heavy_metal);
    let micro = rods.micro.as_ref().unwrap();
    let out_micro = out.micro.as_ref().unwrap();
    // 追う核種と燃焼度の点ごとの均質化した数密度は、棒のない値（履歴の組成で決まる）
    assert_eq!(micro.nuclides, out_micro.nuclides);
    assert_eq!(micro.densities, out_micro.densities);
    for (r, one) in [&out, &a].into_iter().enumerate() {
        let own = burnups(one);
        let one_micro = one.micro.as_ref().unwrap();
        for (n, point) in one.points.iter().enumerate() {
            let mut at = point.clone();
            at[0] = all.iter().position(|&b| b == own[point[0]]).unwrap();
            at.push(r);
            let m = rods.points.iter().position(|p| *p == at).unwrap();
            assert_eq!(rods.states[m], one.states[n], "rod value {r} at {point:?}");
            assert_eq!(
                micro.states[m], one_micro.states[n],
                "micro cross sections of rod value {r} at {point:?}"
            );
            assert!(
                micro.states[m]
                    .iter()
                    .any(|nuclide| nuclide.derivative.is_some()),
                "the density derivative of rod value {r} at {point:?}"
            );
        }
    }
    // 棒のない値は公称を 3 点とほかの状態を 2 点、値 a はどちらの状態も 2 点
    assert_eq!(out.points.len(), 5);
    assert_eq!(a.points.len(), 4);
    assert_eq!(rods.points.len(), 9);
    std::fs::remove_dir_all(&dir).unwrap();
}
