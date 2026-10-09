use super::*;
use crate::input::LatticeInput;

const INPUT: &str = r#"
library = "vera.msgpack"
few-group = { output = "fg.msgpack", bounds = [1e7, 0.625, 1e-5] }
pitch = 1.26
azimuthal = 16
spacing = 0.05
rows = ["pin"]

[material.fuel]
temperature = 565
nuclides = { U-235 = 7.0e-4, U-238 = 2.2e-2, O-16 = 4.6e-2 }

[material.clad]
temperature = 565
nuclides = { Zr-90 = 4.0e-2 }

[material.water]
temperature = 565
water = { density = 0.74271, boron = 1300, b10 = 0.199 }

[cell.pin]
zones = [{ material = "fuel", radius = 0.4, rings = 2 }, { material = "clad", radius = 0.475, rings = 1 },
         { material = "water", rings = 2 }]

[branches]
boron = { ppm = [900, 1700], materials = ["water"] }
moderator = { temperature = [550, 565, 580], density = [0.76972, 0.74271, 0.71190], materials = ["water", "clad"] }
fuel = { temperature = [550, 580], materials = ["fuel"] }
"#;

fn edited(from: &str, to: &str) -> Result<LatticeInput, crate::input::InputError> {
    assert!(INPUT.contains(from), "{from}");
    LatticeInput::parse(&INPUT.replace(from, to))
}

fn material<'a>(materials: &'a [MaterialInput], name: &str) -> &'a MaterialInput {
    materials.iter().find(|m| m.name == name).unwrap()
}

#[test]
fn axes_follow_boron_moderator_fuel() {
    let input = LatticeInput::parse(INPUT).unwrap();
    let branches = input.branches.as_ref().unwrap();
    assert_eq!(
        branches.axes(),
        vec![
            Axis::Boron {
                ppm: vec![900.0, 1700.0],
                b10: 0.199
            },
            Axis::Moderator {
                temperature: vec![550.0, 565.0, 580.0],
                density: vec![0.76972, 0.74271, 0.71190]
            },
            Axis::Fuel {
                temperature: vec![550.0, 580.0]
            },
        ]
    );
}

#[test]
fn states_change_the_listed_materials_with_the_first_axis_slowest() {
    let input = LatticeInput::parse(INPUT).unwrap();
    let states = input.branch_states();
    assert_eq!(states.len(), 2 * 3 * 2);
    // 状態 7 = ホウ素 1（1700 ppm）、減速材 0（550 K）、燃料 1（580 K）
    let s = &states[7];
    let water = material(s, "water");
    assert_eq!(water.temperature, 550.0);
    let expected: Vec<(String, f64)> = borated_water(0.76972, 1700.0, 0.199)
        .unwrap()
        .into_iter()
        .collect();
    assert_eq!(water.composition, expected);
    assert_eq!(material(s, "clad").temperature, 550.0);
    assert_eq!(
        material(s, "clad").composition,
        material(&input.materials, "clad").composition
    );
    assert_eq!(material(s, "fuel").temperature, 580.0);
    assert_eq!(
        material(s, "fuel").composition,
        material(&input.materials, "fuel").composition
    );
    // 最後の状態
    let s = &states[11];
    assert_eq!(material(s, "water").temperature, 580.0);
    assert_eq!(material(s, "fuel").temperature, 580.0);
}

#[test]
fn without_branches_there_is_one_state() {
    let i = INPUT.find("[branches]").unwrap();
    let input = LatticeInput::parse(&INPUT[..i]).unwrap();
    assert!(input.branches.is_none());
    assert_eq!(input.branch_states(), vec![input.materials.clone()]);
}

#[test]
fn branches_need_all_three_axes() {
    // core は状態を与えるとき3つの軸をすべて求める（ADR-0097）ので、一部の軸のテーブルは作らない
    for axis in ["boron = {", "moderator = {", "fuel = {"] {
        let i = INPUT.find(axis).unwrap();
        let j = i + INPUT[i..].find('\n').unwrap() + 1;
        let text = format!("{}{}", &INPUT[..i], &INPUT[j..]);
        let err = LatticeInput::parse(&text).unwrap_err().to_string();
        assert!(err.contains("all three axes"), "{err}");
    }
}

/// エラーになり、その文に `needle` を含むこと
fn rejects(from: &str, to: &str, needle: &str) {
    let err = edited(from, to).unwrap_err().to_string();
    assert!(err.contains(needle), "{err}");
}

#[test]
fn invalid_branches_are_errors() {
    let boron = "ppm = [900, 1700], materials = [\"water\"]";
    rejects(
        boron,
        "ppm = [900, 1700], materials = [\"oil\"]",
        "unknown material `oil`",
    );
    rejects(
        boron,
        "ppm = [900, 1700], materials = [\"water\", \"clad\"]",
        "`clad` is not written as water",
    );
    rejects(
        boron,
        "ppm = [900, 1700], materials = [\"water\", \"water\"]",
        "listed twice",
    );
    rejects(
        boron,
        "ppm = [1700, 900], materials = [\"water\"]",
        "increasing",
    );
    rejects(
        "materials = [\"water\", \"clad\"]",
        "materials = [\"clad\"]",
        "`water`",
    );
    rejects(
        "materials = [\"fuel\"]",
        "materials = [\"fuel\", \"clad\"]",
        "both the moderator and the fuel",
    );
    rejects(
        "materials = [\"fuel\"]",
        "materials = []",
        "list the materials",
    );
    rejects(
        "temperature = [550, 580]",
        "temperature = [580, 550]",
        "increasing",
    );
    rejects(
        "density = [0.76972, 0.74271, 0.71190]",
        "density = [0.76972, 0.74271]",
        "densities",
    );
    rejects(
        "density = [0.76972, 0.74271, 0.71190]",
        "density = [0.76972, nan, 0.71190]",
        "densities",
    );
    rejects(
        "density = [0.76972, 0.74271, 0.71190]",
        "density = [0.76972, 0.0, 0.71190]",
        "densities",
    );
    // ピン出力は分岐と一緒に出さない（状態ごとの形状関数は ADR-0096 でテーブルに入れる）
    rejects(
        "rows = [\"pin\"]",
        "rows = [\"pin\"]\npin-powers = \"pins.msgpack\"",
        "pin-powers",
    );
    let i = INPUT.find("boron = {").unwrap();
    let err = LatticeInput::parse(&INPUT[..i]).unwrap_err().to_string();
    assert!(err.contains("no axes"), "{err}");
}

#[test]
fn every_borated_water_follows_the_boron_and_moderator_axes() {
    // テーブルの状態は炉心の状態なので、ホウ酸水の一部だけを変えるとテーブルと計算が食い違う
    let two = |boron: &str, moderator: &str| {
        let text = INPUT
            .replace(
                "[cell.pin]",
                "[material.gapwater]\ntemperature = 565\nwater = { density = 0.74271, boron = 1300, b10 = 0.199 }\n\n[cell.pin]",
            )
            .replace("ppm = [900, 1700], materials = [\"water\"]", boron)
            .replace("materials = [\"water\", \"clad\"]", moderator);
        LatticeInput::parse(&text)
    };
    let both = "materials = [\"water\", \"gapwater\", \"clad\"]";
    assert!(
        two(
            "ppm = [900, 1700], materials = [\"water\", \"gapwater\"]",
            both
        )
        .is_ok()
    );
    let err = two("ppm = [900, 1700], materials = [\"water\"]", both)
        .unwrap_err()
        .to_string();
    assert!(err.contains("`gapwater`"), "{err}");
    let err = two(
        "ppm = [900, 1700], materials = [\"water\", \"gapwater\"]",
        "materials = [\"water\", \"clad\"]",
    )
    .unwrap_err()
    .to_string();
    assert!(err.contains("`gapwater`"), "{err}");
    // ホウ素の B-10 の割合が物質で違う
    let text = INPUT.replace(
        "[cell.pin]",
        "[material.gapwater]\ntemperature = 565\nwater = { density = 0.74271, boron = 1300, b10 = 0.2 }\n\n[cell.pin]",
    );
    let text = text
        .replace(
            "ppm = [900, 1700], materials = [\"water\"]",
            "ppm = [900, 1700], materials = [\"water\", \"gapwater\"]",
        )
        .replace("materials = [\"water\", \"clad\"]", both);
    let err = LatticeInput::parse(&text).unwrap_err().to_string();
    assert!(err.contains("b10"), "{err}");
}

#[test]
fn branches_need_a_library() {
    let text = INPUT.replace(
        "library = \"vera.msgpack\"",
        "cross-sections = \"c5g7.toml\"",
    );
    let i = text.find("[material.fuel]").unwrap();
    let j = text.find("[cell.pin]").unwrap();
    let text = format!("{}{}", &text[..i], &text[j..]).replace(
        "few-group = { output = \"fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }\n",
        "",
    );
    let err = LatticeInput::parse(&text).unwrap_err().to_string();
    assert!(err.contains("library"), "{err}");
}

#[test]
fn tables_collect_each_output_over_the_states() {
    // 状態ごとの出力の少数群定数（[状態][出力]）を、出力ごとのテーブル（状態は軸の直積の順）にまとめる
    let input = LatticeInput::parse(INPUT).unwrap();
    let axes = input.branches.as_ref().unwrap().axes().to_vec();
    let few = |k: f64| xlc_format::fewgroup::FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: k,
        total: vec![0.5, 1.0],
        absorption: vec![0.01, 0.08],
        nu_fission: vec![0.006, 0.1],
        fission: vec![0.0025, 0.04],
        chi: vec![1.0, 0.0],
        diffusion: vec![1.4, 0.4],
        scatter: vec![vec![0.47, 0.02], vec![0.0, 0.92]],
        adf: vec![vec![1.0, 1.0]; 4],
        vacuum_ratio: vec![None; 4],
        normalized_adf: vec![false; 4],
        pins: None,
    };
    let per_state: Vec<Vec<_>> = (0..12)
        .map(|s| vec![vec![few(1.0 + s as f64)], vec![few(2.0 + s as f64)]])
        .collect();
    let tables = tables(&axes, 1, per_state).unwrap();
    assert_eq!(tables.len(), 2);
    for (o, table) in tables.iter().enumerate() {
        assert_eq!(table.axes, axes);
        let ks: Vec<f64> = table.states.iter().map(|s| s.k_inf).collect();
        assert_eq!(
            ks,
            (0..12)
                .map(|s| 1.0 + o as f64 + s as f64)
                .collect::<Vec<_>>()
        );
    }
    // 状態ごとの出力の数が違う、状態の数が軸と違う
    let ragged: Vec<Vec<_>> = (0..12)
        .map(|s| vec![vec![few(1.0)]; if s == 3 { 1 } else { 2 }])
        .collect();
    assert!(super::tables(&axes, 1, ragged).is_err());
    assert!(super::tables(&axes, 1, vec![vec![vec![few(1.0)]]; 11]).is_err());
    // 集合体を 2 × 2 に分けた出力は、状態ごとに部分ノードの定数を 4 つ持つ（ADR-0108）。数が違えばエラー
    let split: Vec<Vec<_>> = (0..12)
        .map(|s| {
            vec![
                (0..4)
                    .map(|q| few(1.0 + s as f64 + q as f64 / 10.0))
                    .collect(),
            ]
        })
        .collect();
    let tables = super::tables(&axes, 2, split.clone()).unwrap();
    assert_eq!(tables[0].split, 2);
    assert_eq!(tables[0].states.len(), 48);
    assert_eq!(tables[0].states[5].k_inf, 2.1);
    assert!(super::tables(&axes, 3, split).is_err());
}
