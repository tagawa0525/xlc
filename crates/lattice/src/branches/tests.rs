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
        kappa_fission: None,
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
            starts.solved(state).unwrap();
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
        starts.solved(s).unwrap();
        previous = Some(s);
    }
    // (0, 1, 0) の ADR-0195 の状態 (0, 0, 0) は解かないので、直前の (0, 0, 1) から
    assert_eq!(present[..2], [1, 3]);
}

#[test]
fn states_step_back_to_the_previous_point_their_value_solves() {
    // 制御棒の値 × 燃焼度 3 点 × ホウ素 2 点で、棒のない値はすべて、値 1 は燃焼度 0 と 2 だけを解く（ADR-0225）。値 1 の燃焼度 2 の
    // 状態は、その値の 1 つ前の燃焼点（0）の同じ状態から始め、燃焼度だけを間引いた別の入力の格子計算と同じ始める解になる
    let counts = [2, 3, 2];
    let present: Vec<usize> = (0..12)
        .filter(|&s| {
            let at = state_indices(&counts, s);
            at[0] == 0 || at[1] != 1
        })
        .collect();
    let mut starts = Starts::over(&counts, present.clone());
    let mut seen = Vec::new();
    for &s in &present {
        seen.push((state_indices(&counts, s), starts.start().copied()));
        starts.solved(s).unwrap();
    }
    let index = |at: [usize; 3]| state_index(&counts, &at);
    let start_of = |at: [usize; 3]| {
        seen.iter()
            .find(|(a, _)| a[..] == at[..])
            .map(|(_, start)| *start)
            .unwrap()
    };
    assert_eq!(start_of([1, 2, 0]), Some(index([1, 0, 0])));
    assert_eq!(start_of([1, 2, 1]), Some(index([1, 2, 0])));
    // 棒のない値はすべての燃焼点を解くので、1 つ前の燃焼点から
    assert_eq!(start_of([0, 2, 0]), Some(index([0, 1, 0])));
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
        straight.solved(s as f64).unwrap();
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
fn solving_past_the_last_state_is_an_error() {
    let mut starts = Starts::over(&[2], vec![1]);
    starts.solved(0.5).unwrap();
    let err = starts.solved(1.5).unwrap_err();
    assert!(err.contains("solved"), "{err}");
}

/// INPUT の格子を、燃料のピンと案内管の 1 行にし、案内管に入れる棒のセルと制御棒の分岐を足したもの
fn with_rods(rods: &str) -> String {
    let text = INPUT.replace("rows = [\"pin\"]", "rows = [\"pin gt\"]").replace(
        "[branches]",
        "[material.absorber]\ntemperature = 565\nnuclides = { U-238 = 1.0e-2 }\n\n\
         [cell.gt]\nzones = [{ material = \"water\", radius = 0.56, rings = 2 }, { material = \"water\", rings = 1 }]\n\n\
         [cell.rod]\nzones = [{ material = \"absorber\", radius = 0.38, rings = 2 }, { material = \"clad\", radius = 0.56, rings = 1 },\n\
                  { material = \"water\", rings = 1 }]\n\n\
         [cell.other]\nzones = [{ material = \"clad\", rings = 1 }]\n\n[branches]",
    );
    format!("{text}{rods}")
}

#[test]
fn rod_branches_replace_cells_of_the_base_rows() {
    let input = LatticeInput::parse(&with_rods(
        "\n[branches.rod]\nb4c = { gt = \"other\" }\naic = { gt = \"rod\", pin = \"other\" }\n",
    ))
    .unwrap();
    let branches = input.branches.as_ref().unwrap();
    // 連続の軸の状態は今のまま、テーブルの軸は制御棒の軸（out を先に、値の名前の順）を最後に足す
    assert_eq!(branches.axes().len(), 3);
    let values: Vec<String> = ["out", "aic", "b4c"].map(String::from).to_vec();
    assert_eq!(
        branches.table_axes().last(),
        Some(&Axis::Rod {
            values: values.clone()
        })
    );
    assert_eq!(input.rod_values(), values);
    // 値 0 は基準の形状、ほかの値は置き換えたセルを並べた形状
    let base = input.geometry().unwrap();
    assert_eq!(input.rod_geometry(0).unwrap(), base);
    let rows = |cells: &str| {
        let text = with_rods("").replace("rows = [\"pin gt\"]", &format!("rows = [\"{cells}\"]"));
        LatticeInput::parse(&text).unwrap().geometry().unwrap()
    };
    assert_eq!(input.rod_geometry(1).unwrap(), rows("other rod"));
    assert_eq!(input.rod_geometry(2).unwrap(), rows("pin other"));
    // 制御棒の分岐のない入力は値 out だけ
    let plain = LatticeInput::parse(&with_rods("")).unwrap();
    assert_eq!(plain.rod_values(), vec!["out".to_string()]);
    assert_eq!(plain.branches.as_ref().unwrap().table_axes().len(), 3);
}

#[test]
fn invalid_rod_branches_are_errors() {
    let rejects = |rods: &str, needle: &str| {
        let err = LatticeInput::parse(&with_rods(&format!("\n[branches.rod]\n{rods}\n")))
            .unwrap_err()
            .to_string();
        assert!(err.contains(needle), "{rods}: {err}");
    };
    // out は基準の形状の値
    rejects("out = { gt = \"rod\" }", "`out`");
    rejects("aic = {}", "replaces no cell");
    rejects("", "give the rod values");
    rejects("\"a b\" = { gt = \"rod\" }", "without spaces");
    rejects("aic = { box = \"rod\" }", "`box` is not in the rows");
    rejects("aic = { gt = \"bar\" }", "`bar` is not defined");
    // 基準の形状の行にあるセルに置き換えると、棒のない状態に棒が残る
    rejects("aic = { gt = \"pin\" }", "`pin` is in the rows");
}

/// 径方向の反射体の格子の形（2 × 2 の燃料のピンと案内管の右に反射体の列）に、制御棒の分岐（値 aic は案内管 gt を棒に置き換える）と
/// ノードごとの出力 `nodes`（`few-group` の `nodes` の中身）を組んだ入力。`head` は先頭に足す行
fn reflector_with_rods(head: &str, nodes: &str) -> String {
    let text = with_rods("\n[branches.rod]\naic = { gt = \"rod\" }\n")
        .replace(
            "rows = [\"pin gt\"]",
            "rows = [\"pin gt box\", \"gt pin box\"]",
        )
        .replace(
            "few-group = { output = \"fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
            &format!("few-group = {{ bounds = [1e7, 0.625, 1e-5], nodes = [{nodes}] }}"),
        )
        .replace(
            "\n[branches]\n",
            "\n[cell.box]\nzones = [{ material = \"water\", rings = 1 }]\n\n[branches]\n",
        );
    assert!(text.contains("nodes = [") && text.contains("[cell.box]"));
    format!("{head}{text}")
}

const FUEL_AND_REFLECTOR: &str = "{ output = \"fuel.msgpack\", columns = [1, 2], rows = [1, 2] }, \
     { output = \"reflector.msgpack\", columns = [3, 3], rows = [1, 2], reflector = true }";

#[test]
fn rod_axes_go_to_nodes_with_replaced_cells() {
    // 案内管は燃料のノードにだけあるので、燃料のノードのテーブルだけが制御棒の軸を持つ（ADR-0227）
    let input = LatticeInput::parse(&reflector_with_rods("", FUEL_AND_REFLECTOR)).unwrap();
    assert_eq!(input.rod_axes(), vec![true, false]);
    // 余白の薄い行と列もノードの範囲に数える（余白の列 1、燃料の列 2〜3、反射体の列 4 と余白の列 5）
    let gap = reflector_with_rods(
        "gap = { width = 0.04, material = \"water\" }\n",
        "{ output = \"fuel.msgpack\", columns = [1, 3], rows = [1, 4] }, \
         { output = \"reflector.msgpack\", columns = [4, 5], rows = [1, 4], reflector = true }",
    );
    assert_eq!(
        LatticeInput::parse(&gap).unwrap().rod_axes(),
        vec![true, false]
    );
    // 制御棒の分岐がなければ、どのノードも軸を持たない
    let plain = reflector_with_rods("", FUEL_AND_REFLECTOR)
        .replace("[branches.rod]\naic = { gt = \"rod\" }\n", "");
    assert!(!plain.contains("[branches.rod]"));
    assert_eq!(
        LatticeInput::parse(&plain).unwrap().rod_axes(),
        vec![false, false]
    );
    // 格子全体の 1 つのテーブルは、制御棒の分岐があれば軸を持つ
    let lattice =
        LatticeInput::parse(&with_rods("\n[branches.rod]\naic = { gt = \"rod\" }\n")).unwrap();
    assert_eq!(lattice.rod_axes(), vec![true]);
    assert_eq!(
        LatticeInput::parse(&with_rods("")).unwrap().rod_axes(),
        vec![false]
    );
}

#[test]
fn node_rod_branches_outside_the_radial_reflector_lattice_are_errors() {
    let rejects = |text: String, needle: &str| {
        let err = LatticeInput::parse(&text).unwrap_err().to_string();
        assert!(err.contains(needle), "{needle}: {err}");
    };
    // 平板は上の層の値の定数とともに #0030 で通す
    rejects(
        reflector_with_rods("", FUEL_AND_REFLECTOR).replace(
            "bounds = [1e7, 0.625, 1e-5], nodes",
            "bounds = [1e7, 0.625, 1e-5], axial = true, nodes",
        ),
        "#0030",
    );
    // 反射体でないノードが 2 つ（カラーセット）
    rejects(
        reflector_with_rods(
            "",
            "{ output = \"a.msgpack\", columns = [1, 1], rows = [1, 2] }, \
             { output = \"b.msgpack\", columns = [2, 2], rows = [1, 2] }, \
             { output = \"reflector.msgpack\", columns = [3, 3], rows = [1, 2], reflector = true }",
        ),
        "#0214",
    );
    // 反射体のノードが置き換えるセルを含む（余白の列を数えないと、反射体のノードの列 3 は box の列に見える）
    rejects(
        reflector_with_rods(
            "gap = { width = 0.04, material = \"water\" }\n",
            "{ output = \"fuel.msgpack\", columns = [1, 2], rows = [1, 4] }, \
             { output = \"reflector.msgpack\", columns = [3, 5], rows = [1, 4], reflector = true }",
        ),
        "reflector node",
    );
    // 燃料のノードが置き換えるセルを含まない
    rejects(
        reflector_with_rods(
            "",
            "{ output = \"fuel.msgpack\", columns = [1, 1], rows = [1, 1] }, \
             { output = \"reflector.msgpack\", columns = [3, 3], rows = [1, 2], reflector = true }",
        ),
        "no cell",
    );
}

/// 制御棒の分岐（値 aic と b4c）と、燃焼の履歴からの分岐（燃焼点 0、10、20、40）の入力。連続の軸は入力の物質を公称の格子点に
/// 持つ（ADR-0200 の部分格子）。`history` は `[history]` に足す行
fn history_with_rods(history: &str) -> String {
    let text = with_rods(
        "\n[branches.rod]\naic = { gt = \"rod\" }\nb4c = { gt = \"other\" }\n\n[history]\nfile = \"h.msgpack\"\n\
         burnup = [0, 10, 20, 40]\n",
    )
    .replace("ppm = [900, 1700]", "ppm = [900, 1300]")
    .replace("temperature = [550, 580]", "temperature = [565, 580]");
    format!("{text}{history}")
}

#[test]
fn history_branches_take_the_rod_axis() {
    // 燃焼の履歴からの分岐も制御棒の軸を持つ（ADR-0225）。値ごとの燃焼点を書かなければ、制御棒の値の状態は公称でない状態の
    // 燃焼点で解く
    let input =
        LatticeInput::parse(&history_with_rods("off-nominal-burnup = [0, 20, 40]\n")).unwrap();
    assert_eq!(
        input.rod_values(),
        ["out", "aic", "b4c"].map(String::from).to_vec()
    );
    let history = input.history.as_ref().unwrap();
    assert_eq!(history.off_nominal, vec![0.0, 20.0, 40.0]);
    assert_eq!(history.rod_burnup("aic"), &[0.0, 20.0, 40.0][..]);
    assert_eq!(history.rod_burnup("b4c"), &[0.0, 20.0, 40.0][..]);
    // 値ごとに燃焼点を選べる（ADR-0200 の「組ごとに燃焼点を入力で選ぶ」）
    let input = LatticeInput::parse(&history_with_rods(
        "off-nominal-burnup = [0, 40]\nrod-burnup = { aic = [0, 10, 20, 40] }\n",
    ))
    .unwrap();
    let history = input.history.as_ref().unwrap();
    assert_eq!(history.rod_burnup("aic"), &[0.0, 10.0, 20.0, 40.0][..]);
    assert_eq!(history.rod_burnup("b4c"), &[0.0, 40.0][..]);
}

#[test]
fn rod_burnups_are_history_burnups_of_rod_values() {
    for (bad, why) in [
        (
            "rod-burnup = { pyrex = [0, 40] }\n",
            "not a value of [branches.rod]",
        ),
        (
            "rod-burnup = { out = [0, 40] }\n",
            "not a value of [branches.rod]",
        ),
        ("rod-burnup = { aic = [0, 20] }\n", "the last point"),
        ("rod-burnup = { aic = [10, 40] }\n", "the first point"),
        ("rod-burnup = { aic = [0, 15, 40] }\n", "not a point"),
        ("rod-burnup = { aic = [0, 20, 10, 40] }\n", "not increasing"),
        ("rod-burnup = { aic = [] }\n", "empty"),
    ] {
        assert!(
            LatticeInput::parse(&history_with_rods(bad)).is_err(),
            "{why}: {bad}"
        );
    }
    // 制御棒の分岐のない入力の値ごとの燃焼点は、解く値がないので書けない
    let text = history_with_rods("rod-burnup = { aic = [0, 40] }\n").replace(
        "[branches.rod]\naic = { gt = \"rod\" }\nb4c = { gt = \"other\" }\n",
        "",
    );
    assert!(!text.contains("[branches.rod]"));
    assert!(LatticeInput::parse(&text).is_err());
}

#[test]
fn rod_states_are_ordered_lexicographically_with_the_rod_axis_last() {
    // 値ごとに解いた状態（値 0 の 2 点、値 1 の 2 点）を、点の辞書式の順に並べる
    let points = vec![vec![0, 0], vec![1, 0], vec![0, 1], vec![1, 1]];
    assert_eq!(lexicographic(&points), vec![0, 2, 1, 3]);
}
