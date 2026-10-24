use super::*;
use crate::input::LatticeInput;
use xlc_format::table::product_points;

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
            Axis::boron(vec![900.0, 1700.0], 0.199),
            Axis::moderator(vec![550.0, 565.0, 580.0], vec![0.76972, 0.74271, 0.71190]),
            Axis::fuel(vec![550.0, 580.0]),
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
fn references_keep_the_nominal_boron_at_the_moderator_and_fuel_points_of_each_state() {
    // 自己遮蔽を公称のホウ素で解くとき（ADR-0201）の参照: 減速材と燃料の温度の軸は状態の格子点、ホウ素は入力のまま
    let input = LatticeInput::parse(INPUT).unwrap();
    let branches = input.branches.as_ref().unwrap();
    let states = input.branch_states();
    // 減速材 3 点 × 燃料 2 点
    assert_eq!(branches.references().len(), 6);
    for (state, s) in states.iter().enumerate() {
        let reference = &branches.references()[branches.reference(state).unwrap()];
        assert_eq!(material(reference, "fuel"), material(s, "fuel"), "{state}");
        assert_eq!(material(reference, "clad"), material(s, "clad"), "{state}");
        let (water, nominal) = (material(reference, "water"), material(s, "water"));
        assert_eq!(water.temperature, nominal.temperature, "{state}");
        // ホウ酸水は状態の密度で、ホウ素だけが入力の 1300 ppm
        let density = [0.76972, 0.74271, 0.71190][(state / 2) % 3];
        let expected: Vec<(String, f64)> = borated_water(density, 1300.0, 0.199)
            .unwrap()
            .into_iter()
            .collect();
        assert_eq!(water.composition, expected, "{state}");
    }
    // ホウ素の点ごとに同じ参照を使う（ホウ素は最初の軸）
    assert_eq!(branches.reference(7), branches.reference(1));
    // 参照から作り直す物質は、ホウ素の軸が変える物質
    assert_eq!(branches.varied(), ["water".to_string()]);
    // 状態ごとに解き直すなら参照を使わない
    let every = edited(
        "fuel = { temperature = [550, 580], materials = [\"fuel\"] }",
        "fuel = { temperature = [550, 580], materials = [\"fuel\"] }\nshielding = \"every-state\"",
    )
    .unwrap();
    let every = every.branches.as_ref().unwrap();
    assert!((0..12).all(|state| every.reference(state).is_none()));
    // ホウ素が 1 点なら、状態ごとに参照が違い、再利用するものがないので状態ごとに解く
    let single = edited(
        "boron = { ppm = [900, 1700], materials = [\"water\"] }",
        "boron = { ppm = [900], materials = [\"water\"] }",
    )
    .unwrap();
    let single = single.branches.as_ref().unwrap();
    assert!((0..6).all(|state| single.reference(state).is_none()));
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

/// k∞ が `k` の 2 群の少数群定数
fn few(k: f64) -> xlc_format::fewgroup::FewGroup {
    xlc_format::fewgroup::FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: k,
        total: vec![0.5, 1.0],
        absorption: vec![0.01, 0.08],
        nu_fission: vec![0.006, 0.1],
        fission: vec![0.0025, 0.04],
        chi: vec![1.0, 0.0],
        diffusion: vec![1.4, 0.4],
        scatter: vec![vec![0.47, 0.02], vec![0.0, 0.92]],
        adf: vec![vec![1.0, 1.0]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: None,
    }
}

#[test]
fn tables_collect_each_output_over_the_states() {
    // 状態ごとの出力の少数群定数（[状態][出力]）を、出力ごとのテーブル（状態は軸の直積の順）にまとめる
    let input = LatticeInput::parse(INPUT).unwrap();
    let axes = input.branches.as_ref().unwrap().axes().to_vec();
    let per_state: Vec<Vec<_>> = (0..12)
        .map(|s| vec![vec![few(1.0 + s as f64)], vec![few(2.0 + s as f64)]])
        .collect();
    let points = product_points(&axes);
    let tables = tables(&axes, 1, &points, per_state).unwrap();
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
    assert!(super::tables(&axes, 1, &points, ragged).is_err());
    assert!(super::tables(&axes, 1, &points, vec![vec![vec![few(1.0)]]; 11]).is_err());
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
    let tables = super::tables(&axes, 2, &points, split.clone()).unwrap();
    assert_eq!(tables[0].split, 2);
    assert_eq!(tables[0].states.len(), 48);
    assert_eq!(tables[0].states[5].k_inf, 2.1);
    assert!(super::tables(&axes, 3, &points, split).is_err());
}

#[test]
fn each_state_starts_from_one_step_back_on_its_last_moved_axis() {
    // 状態は軸の直積の順（最初の軸がいちばん遅く変わる）。状態 0 は始める解がない。ほかの状態は、添字が 0 でない最後の軸を
    // 1 つ戻し、それより後の軸を 0 にした状態（その 1 つの軸だけが違う、前に解いた状態）から始める（ADR-0195）
    let cases: [(&[usize], &[Option<usize>]); 3] = [
        (
            &[2, 2, 2],
            &[
                None,
                Some(0),
                Some(0),
                Some(2),
                Some(0),
                Some(4),
                Some(4),
                Some(6),
            ],
        ),
        // 格子点が 1 つの軸は動かない
        (&[2, 1, 2], &[None, Some(0), Some(0), Some(2)]),
        // 軸のない分岐（状態 1 つ）
        (&[], &[None]),
    ];
    for (counts, expected) in cases {
        let mut starts = Starts::new(counts);
        for (state, source) in expected.iter().enumerate() {
            assert_eq!(
                starts.start().copied(),
                *source,
                "{counts:?}: state {state}"
            );
            starts.solved(state);
        }
    }
}

#[test]
fn restored_starts_need_every_axis_after_the_first_state() {
    // 状態 0 を解くと、すべての軸が解を持つ。書いた解の欠けた軸は、後の状態を平らな初期値から黙って始めさせる
    assert!(Starts::<f64>::restored(&[2, 3], (0..6).collect(), 0, vec![], &[None, None]).is_ok());
    for slots in [[None, Some(0)], [Some(0), None]] {
        let error = Starts::restored(&[2, 3], (0..6).collect(), 1, vec![0.5], &slots).unwrap_err();
        assert!(error.contains("no start"), "{error}");
    }
    let error =
        Starts::restored(&[2, 3], (0..6).collect(), 0, vec![0.5], &[Some(0), Some(0)]).unwrap_err();
    assert!(error.contains("state 0"), "{error}");
}

/// INPUT の軸を 3 点ずつ（入力の物質は中央の格子点）にし、`subgrids` を足した入力
fn three_by_three(subgrids: &str) -> String {
    let text = INPUT
        .replace(
            "boron = { ppm = [900, 1700], materials = [\"water\"] }",
            "boron = { ppm = [900, 1300, 1700], materials = [\"water\"] }",
        )
        .replace(
            "fuel = { temperature = [550, 580], materials = [\"fuel\"] }",
            "fuel = { temperature = [550, 565, 580], materials = [\"fuel\"] }",
        );
    format!("{text}{subgrids}\n")
}

const PAIRS: &str =
    r#"subgrids = [["boron", "moderator"], ["boron", "fuel"], ["moderator", "fuel"]]"#;

#[test]
fn the_nominal_state_is_the_grid_point_of_the_input_materials() {
    let input = LatticeInput::parse(&three_by_three("")).unwrap();
    let branches = input.branches.as_ref().unwrap();
    assert_eq!(branches.nominal(), Some(&[1, 1, 1][..]));
    // 状態の直積の番号（ホウ素 1、減速材 1、燃料 1）
    assert_eq!(branches.nominal_state(), Some(13));
    assert_eq!(
        branches.states()[13],
        input.materials,
        "the nominal state is the input"
    );
    // テーブルの軸は公称の値と、ADR-0097 の内挿の変数を持つ
    assert_eq!(
        branches.axes(),
        vec![
            Axis::boron(vec![900.0, 1300.0, 1700.0], 0.199).with_nominal(1300.0),
            Axis::moderator(vec![550.0, 565.0, 580.0], vec![0.76972, 0.74271, 0.71190])
                .with_nominal(565.0),
            Axis::fuel(vec![550.0, 565.0, 580.0]).with_nominal(565.0),
        ]
    );
    // 入力の物質が格子点にない軸があれば、公称の状態はない（INPUT のホウ素は 900 と 1700 ppm で、入力は 1300 ppm）
    let off = LatticeInput::parse(INPUT).unwrap();
    let off = off.branches.as_ref().unwrap();
    assert_eq!(off.nominal(), None);
    assert!(off.axes().iter().all(|a| a.nominal_index().is_none()));
    // 格子点が 1 つの軸は、その点が公称（入力の物質と違ってよい）
    let single = edited(
        "fuel = { temperature = [550, 580], materials = [\"fuel\"] }",
        "fuel = { temperature = [600], materials = [\"fuel\"] }",
    )
    .unwrap();
    let single = single.branches.as_ref().unwrap();
    assert_eq!(single.nominal(), None, "boron is still off the grid");
    let both = LatticeInput::parse(&three_by_three("").replace(
        "fuel = { temperature = [550, 565, 580], materials = [\"fuel\"] }",
        "fuel = { temperature = [600], materials = [\"fuel\"] }",
    ))
    .unwrap();
    assert_eq!(
        both.branches.as_ref().unwrap().nominal(),
        Some(&[1, 1, 0][..])
    );
}

#[test]
fn branches_without_a_history_solve_the_product_unless_told_the_subgrids() {
    // 燃焼度の軸のない分岐のテーブル（M3）は直積（ADR-0097）
    let input = LatticeInput::parse(&three_by_three("")).unwrap();
    let branches = input.branches.as_ref().unwrap();
    assert!((0..27).all(|s| branches.solved(s)));
    assert!(branches.is_product());
    // 公称の状態のない入力も直積なら解ける
    let off = LatticeInput::parse(INPUT).unwrap();
    assert!((0..12).all(|s| off.branches.as_ref().unwrap().solved(s)));
}

#[test]
fn pair_subgrids_solve_the_nominal_state_one_axis_branches_and_pairs() {
    // 1 + Σ(nᵢ − 1) + Σᵢ<ⱼ(nᵢ − 1)(nⱼ − 1) = 1 + 6 + 12 状態（ADR-0200）。3 つの軸がどれも公称でない状態だけを解かない
    let input = LatticeInput::parse(&three_by_three(PAIRS)).unwrap();
    let branches = input.branches.as_ref().unwrap();
    assert!(!branches.is_product());
    let solved: Vec<usize> = (0..27).filter(|&s| branches.solved(s)).collect();
    assert_eq!(solved.len(), 19);
    for s in 0..27 {
        let at = state_indices(&[3, 3, 3], s);
        assert_eq!(branches.solved(s), at.contains(&1), "{at:?}");
    }
    // 3 軸の部分格子を持てば直積、1 軸ずつなら 1 + 6
    let triple = LatticeInput::parse(&three_by_three(
        r#"subgrids = [["boron", "moderator", "fuel"]]"#,
    ))
    .unwrap();
    assert!(triple.branches.as_ref().unwrap().is_product());
    let one = LatticeInput::parse(&three_by_three(
        r#"subgrids = [["boron"], ["moderator"], ["fuel"]]"#,
    ))
    .unwrap();
    let one = one.branches.as_ref().unwrap();
    assert_eq!((0..27).filter(|&s| one.solved(s)).count(), 7);
    // 公称が端の格子点で、軸の点の数が違う（ADR-0200 の 2A の 40 MWd/kg の形: 3 × 2 × 2 で 10 状態）
    let uneven = LatticeInput::parse(
        &(INPUT
            .replace(
                "boron = { ppm = [900, 1700], materials = [\"water\"] }",
                "boron = { ppm = [900, 1100, 1300], materials = [\"water\"] }",
            )
            .replace(
                "moderator = { temperature = [550, 565, 580], density = [0.76972, 0.74271, 0.71190], materials = [\"water\", \"clad\"] }",
                "moderator = { temperature = [565, 580], density = [0.74271, 0.71190], materials = [\"water\", \"clad\"] }",
            )
            .replace(
                "fuel = { temperature = [550, 580], materials = [\"fuel\"] }",
                "fuel = { temperature = [565, 580], materials = [\"fuel\"] }",
            )
            + PAIRS
            + "\n"),
    )
    .unwrap();
    let uneven = uneven.branches.as_ref().unwrap();
    assert_eq!(uneven.nominal(), Some(&[2, 0, 0][..]));
    assert_eq!((0..12).filter(|&s| uneven.solved(s)).count(), 10);
}

#[test]
fn subgrids_need_the_nominal_state_and_name_distinct_axes() {
    // 入力の物質が格子点にない（INPUT のホウ素）と、部分格子を組み立てる公称の状態がない
    let text = format!("{INPUT}{PAIRS}\n");
    let err = LatticeInput::parse(&text).unwrap_err().to_string();
    assert!(err.contains("nominal") && err.contains("boron"), "{err}");
    for bad in [
        r#"subgrids = []"#,
        r#"subgrids = [[]]"#,
        r#"subgrids = [["boron", "boron"]]"#,
        r#"subgrids = [["boron", "rod"]]"#,
        // 点が 2 つ以上の燃料の軸を、どの組も持たない（燃料の分岐を解かず、core が黙って公称の燃料の値を使う）
        r#"subgrids = [["boron", "moderator"]]"#,
    ] {
        assert!(LatticeInput::parse(&three_by_three(bad)).is_err(), "{bad}");
    }
}

#[test]
fn references_of_subgrid_states_are_their_moderator_and_fuel_points() {
    // 公称のホウ素の自己遮蔽（ADR-0201）の参照は、部分格子の状態でも減速材と燃料の温度の格子点で選び、その参照の状態（ホウ素を
    // 公称にした状態）も部分格子の状態にある
    let input = LatticeInput::parse(&three_by_three(PAIRS)).unwrap();
    let branches = input.branches.as_ref().unwrap();
    let mut used = Vec::new();
    for s in (0..27).filter(|&s| branches.solved(s)) {
        let at = state_indices(&[3, 3, 3], s);
        let r = branches.reference(s).unwrap();
        let reference = &branches.references()[r];
        assert_eq!(state_indices(&[3, 3], r), at[1..], "{at:?}");
        let nominal_boron = 9 + r;
        assert!(branches.solved(nominal_boron), "{at:?}");
        assert_eq!(reference, &branches.states()[nominal_boron], "{at:?}");
        if !used.contains(&r) {
            used.push(r);
        }
    }
    // 自己遮蔽は、解く状態の減速材と燃料の点の数（9）だけ解く
    assert_eq!(used.len(), 9);
}

#[test]
fn media_need_the_product_of_the_branches() {
    // 細かい群の媒質のファイルは軸の直積の状態を持つ（ADR-0128）
    let text = format!("medium = \"fuel.medium\"\n{}", three_by_three(PAIRS));
    let err = LatticeInput::parse(&text).unwrap_err().to_string();
    assert!(err.contains("medium"), "{err}");
    let text = format!("medium = \"fuel.medium\"\n{}", three_by_three(""));
    assert!(LatticeInput::parse(&text).is_ok());
}

#[test]
fn tables_hold_the_solved_points() {
    let input = LatticeInput::parse(&three_by_three(PAIRS)).unwrap();
    let branches = input.branches.as_ref().unwrap();
    let axes = branches.axes().to_vec();
    let solved: Vec<usize> = (0..27).filter(|&s| branches.solved(s)).collect();
    let points: Vec<Vec<usize>> = solved
        .iter()
        .map(|&s| state_indices(&[3, 3, 3], s))
        .collect();
    let per_state: Vec<Vec<Vec<_>>> = (0..19).map(|s| vec![vec![few(1.0 + s as f64)]]).collect();
    let tables = tables(&axes, 1, &points, per_state.clone()).unwrap();
    assert_eq!(tables[0].points, points);
    tables[0].validate().unwrap();
    // 点の数と状態の数が違う
    assert!(super::tables(&axes, 1, &points[1..], per_state).is_err());
}

#[test]
fn states_outside_the_subgrids_start_from_the_last_solved_state() {
    // 部分格子の状態だけを解くとき、ADR-0195 の状態（添字が 0 でない最後の軸を 1 つ戻し、後の軸を 0 にした状態）を解いていれば
    // その解から、解いていなければ直前に解いた状態の解から始める。3 × 3 × 3 の 2 軸の組（公称は中央）
    let counts = [3, 3, 3];
    let present: Vec<usize> = (0..27)
        .filter(|&s| state_indices(&counts, s).contains(&1))
        .collect();
    let mut starts = Starts::over(&counts, present.clone());
    let mut previous = None;
    for &s in &present {
        let at = state_indices(&counts, s);
        let source = at.iter().rposition(|&i| i > 0).map(|a| {
            let mut q = at.clone();
            q[a] -= 1;
            q[a + 1..].iter_mut().for_each(|i| *i = 0);
            q.iter().fold(0, |n, &i| n * 3 + i)
        });
        let want = match source {
            Some(q) if present.contains(&q) => Some(q),
            Some(_) => previous,
            None => None,
        };
        assert_eq!(starts.start().copied(), want, "{at:?}");
        starts.solved(s);
        previous = Some(s);
    }
    // (0, 1, 0) の ADR-0195 の状態 (0, 0, 0) は解かないので、直前の (0, 0, 1) から
    assert_eq!(present[..2], [1, 3]);
}

#[test]
fn restored_starts_follow_the_solved_states() {
    // 書いた解の番号は、解く状態の並びから決まる軸ごとの状態と合うこと
    let counts = [3, 3, 3];
    let present: Vec<usize> = (0..27)
        .filter(|&s| state_indices(&counts, s).contains(&1))
        .collect();
    let mut straight = Starts::over(&counts, present.clone());
    for &s in &present[..5] {
        straight.solved(s as f64);
    }
    let (next, solutions, slots) = straight.saved();
    assert_eq!(next, 5);
    let solutions: Vec<f64> = solutions.into_iter().copied().collect();
    let restored =
        Starts::restored(&counts, present.clone(), next, solutions.clone(), &slots).unwrap();
    assert_eq!(restored.start(), straight.start());
    // 欠けた解
    let mut missing = slots.clone();
    let k = missing.iter().position(Option::is_some).unwrap();
    missing[k] = None;
    assert!(Starts::restored(&counts, present, next, solutions, &missing).is_err());
}

#[test]
fn state_index_inverts_state_indices() {
    let counts = [3, 1, 2];
    for state in 0..6 {
        assert_eq!(state_index(&counts, &state_indices(&counts, state)), state);
    }
}

#[test]
#[should_panic(expected = "solved")]
fn solving_past_the_last_state_is_a_bug() {
    let mut starts = Starts::over(&[2], vec![1]);
    starts.solved(0.5);
    starts.solved(1.5);
}
