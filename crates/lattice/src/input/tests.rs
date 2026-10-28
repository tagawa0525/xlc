use super::*;
use crate::branches::BranchShielding;

const PIN_CELL: &str = r#"
# pin cell
library = "m1.msgpack"
few-group = { output = "pincell.fg.msgpack", bounds = [1e7, 0.625, 1e-5] }
pitch = 1.26
azimuthal = 16
spacing = 0.05
rows = ["pin"]

[material.fuel]   # UO2
temperature = 900
nuclides = { U-235 = 7.0e-4, U-238 = 2.2e-2, O-16 = 4.6e-2 }

[material.water]
temperature = 600
nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }

[cell.pin]
zones = [{ material = "fuel", radius = 0.4, rings = 2 }, { material = "water", rings = 2 }]
"#;

const LATTICE: &str = r#"
cross-sections = "c5g7.toml"
pin-powers = "pins.msgpack"
pitch = 2.0
azimuthal = 8
spacing = 0.1
boundary = { right = "vacuum", bottom = "vacuum" }
rows = ["pin box", "box pin"]

[cell.pin]
zones = [
  { material = "UO2", radius = 0.5, rings = 1 },
  { material = "clad", radius = 0.6, rings = 1 },
  { material = "water", rings = 1 },
]

[cell.box]
zones = [{ material = "water", rings = 1 }]
"#;

fn replaced(text: &str, from: &str, to: &str) -> Result<LatticeInput, InputError> {
    assert!(text.contains(from), "{from}");
    LatticeInput::parse(&text.replace(from, to))
}

#[test]
fn symmetry_is_detected_unless_given() {
    let given = |value: &str| {
        replaced(
            PIN_CELL,
            "spacing = 0.05",
            &format!("spacing = 0.05\nsymmetry = \"{value}\""),
        )
    };
    // 省略すると成り立つ対称操作をすべて使い、名前を書くとその対称性を強制する（ADR-0112）
    assert_eq!(
        LatticeInput::parse(PIN_CELL).unwrap().symmetry,
        SymmetryChoice::Detect
    );
    assert_eq!(
        given("quarter").unwrap().symmetry,
        SymmetryChoice::Force(Symmetry::QUARTER)
    );
    assert_eq!(
        given("octant").unwrap().symmetry,
        SymmetryChoice::Force(Symmetry::OCTANT)
    );
    assert_eq!(
        given("none").unwrap().symmetry,
        SymmetryChoice::Force(Symmetry::NONE)
    );
    assert!(given("half").is_err());
    assert!(given("auto").is_err());
}

/// 媒質の燃料と、鋼と水を混ぜた板と水を並べた平板（ADR-0128）
const SLAB: &str = r#"
library = "m1.msgpack"
self-shielding = false
pitch = 1.26
azimuthal = 16
spacing = 0.05
boundary = { left = "vacuum", right = "vacuum" }
rows = ["fuel plate water"]

[material.fuel]
medium = "fuel.medium"

[material.ss304]
temperature = 600
nuclides = { Fe-56 = 8.0e-2 }

[material.water]
temperature = 600
water = { density = 0.743, boron = 1300, b10 = 0.199 }

[material.plate]
mix = [{ material = "ss304", volume = 0.5 }, { material = "water", volume = 0.5 }]

[cell.fuel]
zones = [{ material = "fuel", rings = 1 }]

[cell.plate]
zones = [{ material = "plate", rings = 1 }]

[cell.water]
zones = [{ material = "water", rings = 1 }]
"#;

#[test]
fn materials_may_be_media_or_mixtures() {
    // 物質は、核種の組成、ホウ酸水、媒質のファイル、ほかの物質を体積の割合で混ぜたもの、のどれか（ADR-0128）
    let input = LatticeInput::parse(SLAB).unwrap();
    assert_eq!(
        input.media,
        vec![MediumInput {
            name: "fuel".into(),
            file: PathBuf::from("fuel.medium"),
        }]
    );
    assert_eq!(
        input.mixtures,
        vec![Mixture {
            name: "plate".into(),
            parts: vec![("ss304".into(), 0.5), ("water".into(), 0.5)],
        }]
    );
    let names: Vec<&str> = input.materials.iter().map(|m| m.name.as_str()).collect();
    assert_eq!(names, ["ss304", "water"]);
    let (_, regions) = input.geometry().unwrap();
    assert_eq!(regions, ["fuel", "plate", "water"]);
    // 真空も成分にできる
    let void = replaced(
        SLAB,
        "{ material = \"water\", volume = 0.5 }",
        "{ material = \"void\", volume = 0.5 }",
    )
    .unwrap();
    assert_eq!(void.mixtures[0].parts[1], ("void".to_string(), 0.5));
}

#[test]
fn invalid_media_and_mixtures_are_errors() {
    let bad = |from: &str, to: &str, word: &str| {
        let err = replaced(SLAB, from, to).unwrap_err().to_string();
        assert!(err.contains(word), "{to}: {err}");
    };
    let medium = "medium = \"fuel.medium\"";
    let parts =
        "mix = [{ material = \"ss304\", volume = 0.5 }, { material = \"water\", volume = 0.5 }]";
    // 媒質と混ぜた物質は温度を持たず、ほかの書き方と排他
    bad(
        medium,
        &format!("temperature = 600\n{medium}"),
        "temperature",
    );
    bad(parts, &format!("temperature = 600\n{parts}"), "temperature");
    bad(
        medium,
        &format!("{medium}\nnuclides = {{ U-235 = 1e-3 }}"),
        "fuel",
    );
    bad(parts, &format!("{parts}\n{medium}"), "plate");
    // 割合の和が 1 でない、割合が正でない、成分が空か重なる
    bad("volume = 0.5 }]", "volume = 0.6 }]", "sum");
    bad("volume = 0.5 }, {", "volume = -0.5 }, {", "volume");
    bad(parts, "mix = []", "plate");
    bad(
        "\"water\", volume = 0.5 }]",
        "\"ss304\", volume = 0.5 }]",
        "twice",
    );
    // 成分は核種かホウ酸水で書いた物質か真空で、知らない物質、媒質、混ぜた物質はエラー
    bad(
        "\"water\", volume = 0.5 }]",
        "\"steel\", volume = 0.5 }]",
        "steel",
    );
    bad(
        "\"water\", volume = 0.5 }]",
        "\"fuel\", volume = 0.5 }]",
        "fuel",
    );
    bad(
        "\"water\", volume = 0.5 }]",
        "\"plate\", volume = 0.5 }]",
        "plate",
    );
    // 真空だけの混ぜた物質は物質でない（帯に `void` と書く）
    bad(
        parts,
        "mix = [{ material = \"void\", volume = 1.0 }]",
        "void",
    );
    // 媒質と混ぜた物質は自己遮蔽の対象にならないので、平板の格子は自己遮蔽しない（ADR-0128）
    bad("self-shielding = false\n", "", "self-shielding");
}

#[test]
fn branches_list_mixtures_with_their_components() {
    // 分岐の軸は、変える成分を持つ混ぜた物質を挙げ、媒質は燃料の温度の軸に挙げる（媒質の状態を選ぶ。ADR-0128）
    let branches = "\n[branches]\n\
            boron = { ppm = [1300], materials = [\"water\", \"plate\"] }\n\
            moderator = { temperature = [600], density = [0.743], materials = [\"water\", \"ss304\", \"plate\"] }\n\
            fuel = { temperature = [600], materials = [\"fuel\"] }\n";
    let with = |text: String| LatticeInput::parse(&format!("{SLAB}{text}"));
    assert!(with(branches.into()).is_ok());
    let bad = |from: &str, to: &str, word: &str| {
        assert!(branches.contains(from), "{from}");
        let err = with(branches.replace(from, to)).unwrap_err().to_string();
        assert!(err.contains(word), "{to}: {err}");
    };
    // 変える成分を持つ混ぜた物質が挙がっていない
    bad("[\"water\", \"plate\"] }", "[\"water\"] }", "plate");
    bad("\"ss304\", \"plate\"]", "\"ss304\"]", "plate");
    // 挙げた混ぜた物質に、その軸が変える成分がない
    bad("[\"fuel\"]", "[\"fuel\", \"plate\"]", "plate");
    // 媒質は燃料の温度の軸に挙げる（媒質の状態を選ぶ）
    let without = branches
        .replace("\"ss304\", \"plate\"]", "\"plate\"]")
        .replace("[\"fuel\"]", "[\"ss304\", \"plate\"]");
    let err = with(without).unwrap_err().to_string();
    assert!(err.contains("medium") && err.contains("fuel"), "{err}");
    // 媒質を燃料の温度の軸のほかに挙げない
    bad(
        "\"ss304\", \"plate\"]",
        "\"ss304\", \"plate\", \"fuel\"]",
        "fuel",
    );
}

#[test]
fn mixtures_may_follow_both_temperature_axes() {
    // 減速材の温度の軸が変える成分と、燃料の温度の軸が変える成分を両方持つ混ぜた物質は、両方の軸に挙げる
    // （ADR-0128。混ぜた物質は成分から作るので、同じ物質を2つの温度の軸に挙げない規則は成分にだけ当てはめる）
    let text = SLAB.replace(
        "[cell.fuel]",
        "[material.pellet]\ntemperature = 600\nnuclides = { U-238 = 2.2e-2 }\n\n\
             [material.both]\nmix = [{ material = \"pellet\", volume = 0.5 }, { material = \"water\", volume = 0.5 }]\n\n\
             [cell.fuel]",
    );
    let branches = "\n[branches]\n\
            boron = { ppm = [1300], materials = [\"water\", \"plate\", \"both\"] }\n\
            moderator = { temperature = [600], density = [0.743], materials = [\"water\", \"ss304\", \"plate\", \"both\"] }\n\
            fuel = { temperature = [600], materials = [\"fuel\", \"pellet\", \"both\"] }\n";
    let input = LatticeInput::parse(&format!("{text}{branches}"));
    assert!(input.is_ok(), "{input:?}");
}

#[test]
fn single_lattices_may_write_a_medium() {
    // 外周がすべて鏡面反射の格子は、細かい群のまま均質化した媒質を書ける（ADR-0128）
    let with = |text: &str, extra: &str| replaced(text, "pitch = ", &format!("{extra}\npitch = "));
    let input = with(PIN_CELL, "medium = \"fuel.medium\"").unwrap();
    assert_eq!(input.medium, Some(PathBuf::from("fuel.medium")));
    assert_eq!(LatticeInput::parse(PIN_CELL).unwrap().medium, None);
    // 鏡面反射でない辺があると、媒質の無限媒質の解が格子の解と一致しない
    let open = with(
        PIN_CELL,
        "medium = \"fuel.medium\"\nboundary = { right = \"vacuum\" }",
    );
    let err = open.unwrap_err().to_string();
    assert!(
        err.contains("medium") && err.contains("reflective"),
        "{err}"
    );
    // 群のエネルギーと散乱の次数は、ライブラリから来る
    let macroscopic = replaced(
        LATTICE,
        "boundary = { right = \"vacuum\", bottom = \"vacuum\" }",
        "medium = \"fuel.medium\"",
    );
    let err = macroscopic.unwrap_err().to_string();
    assert!(err.contains("medium") && err.contains("library"), "{err}");
}

#[test]
fn void_zones_need_no_material() {
    let input = replaced(
        PIN_CELL,
        "{ material = \"fuel\", radius = 0.4, rings = 2 }",
        "{ material = \"fuel\", radius = 0.4, rings = 2 }, { material = \"void\", radius = 0.418, rings = 1 }",
    )
    .unwrap();
    // 燃料 2、真空 1、水 2 の領域
    let (_, regions) = input.geometry().unwrap();
    assert_eq!(regions[2], VOID);
}

#[test]
fn materials_need_some_nuclide() {
    let result = replaced(
        PIN_CELL,
        "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
        "nuclides = { H-1 = 0.0, O-16 = 0.0 }",
    );
    assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
}

#[test]
fn void_cannot_be_defined() {
    let result = replaced(PIN_CELL, "[material.water]", "[material.void]");
    assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
}

#[test]
fn self_shielding_is_on_unless_turned_off() {
    assert!(LatticeInput::parse(PIN_CELL).unwrap().self_shielding);
    let off = replaced(
        PIN_CELL,
        "library = \"m1.msgpack\"",
        "library = \"m1.msgpack\"\nself-shielding = false",
    )
    .unwrap();
    assert!(!off.self_shielding);
}

#[test]
fn track_table_defaults_to_the_adr_density() {
    // ADR-0156: 既定は 10 倍あたり 16 点、整数で細かさを選び、false で格子点ごとに求める
    let density = |line: &str| {
        replaced(
            PIN_CELL,
            "library = \"m1.msgpack\"",
            &format!("library = \"m1.msgpack\"\n{line}"),
        )
        .map(|input| input.collision)
    };
    assert_eq!(
        LatticeInput::parse(PIN_CELL).unwrap().collision,
        Collision::Table { density: 16 }
    );
    assert_eq!(
        density("track-table = 32").unwrap(),
        Collision::Table { density: 32 }
    );
    assert_eq!(density("track-table = false").unwrap(), Collision::Exact);
    for bad in ["track-table = true", "track-table = 0", "track-table = -4"] {
        assert!(density(bad).is_err(), "{bad}");
    }
}

#[test]
fn track_table_needs_self_shielding() {
    // 飛程の表は自己遮蔽の減速計算でだけ使うので、自己遮蔽しない入力に書くと止める
    let off = replaced(
        PIN_CELL,
        "library = \"m1.msgpack\"",
        "library = \"m1.msgpack\"\nself-shielding = false\ntrack-table = 32",
    );
    assert!(matches!(off, Err(InputError::Invalid(_))), "{off:?}");
    let macroscopic = replaced(
        LATTICE,
        "cross-sections = \"c5g7.toml\"",
        "cross-sections = \"c5g7.toml\"\ntrack-table = false",
    );
    assert!(
        matches!(macroscopic, Err(InputError::Invalid(_))),
        "{macroscopic:?}"
    );
}

#[test]
fn ufg_tolerance_thins_only_when_given() {
    // ADR-0185: 省略すると間引かず、0 と 1 の間の値で間引く
    let thinning = |line: &str| {
        replaced(
            PIN_CELL,
            "library = \"m1.msgpack\"",
            &format!("library = \"m1.msgpack\"\n{line}"),
        )
        .map(|input| input.shielding_options().thinning)
    };
    assert_eq!(LatticeInput::parse(PIN_CELL).unwrap().thinning, None);
    assert_eq!(thinning("ufg-tolerance = 1e-4").unwrap(), Some(1.0e-4));
    for bad in [
        "ufg-tolerance = 0",
        "ufg-tolerance = 1.5",
        "ufg-tolerance = -1e-4",
    ] {
        assert!(
            matches!(thinning(bad), Err(InputError::Invalid(_))),
            "{bad}"
        );
    }
    let off = replaced(
        PIN_CELL,
        "library = \"m1.msgpack\"",
        "library = \"m1.msgpack\"\nself-shielding = false\nufg-tolerance = 1e-4",
    );
    assert!(matches!(off, Err(InputError::Invalid(_))), "{off:?}");
}

#[test]
fn ufg_hold_memory_sets_the_held_memory_of_the_self_shielding() {
    // ADR-0213: 省略すると使用可能メモリから求め（`None`）、GiB で与えるとその量。負と有限でない値、自己遮蔽をしない入力ではエラー
    let hold = |line: &str| {
        replaced(
            PIN_CELL,
            "library = \"m1.msgpack\"",
            &format!("library = \"m1.msgpack\"\n{line}"),
        )
        .map(|input| input.shielding_options().hold_memory)
    };
    assert_eq!(
        LatticeInput::parse(PIN_CELL)
            .unwrap()
            .shielding_options()
            .hold_memory,
        None
    );
    assert_eq!(hold("ufg-hold-memory = 80").unwrap(), Some(80 << 30));
    assert_eq!(hold("ufg-hold-memory = 0.5").unwrap(), Some(1 << 29));
    assert_eq!(hold("ufg-hold-memory = 0").unwrap(), Some(0));
    for bad in [
        "ufg-hold-memory = -1",
        "ufg-hold-memory = inf",
        "ufg-hold-memory = nan",
    ] {
        assert!(matches!(hold(bad), Err(InputError::Invalid(_))), "{bad}");
    }
    let off = replaced(
        PIN_CELL,
        "library = \"m1.msgpack\"",
        "library = \"m1.msgpack\"\nself-shielding = false\nufg-hold-memory = 80",
    );
    assert!(matches!(off, Err(InputError::Invalid(_))), "{off:?}");
}

#[test]
fn ufg_scattering_defaults_to_the_fast_free_gas_with_sab() {
    // ADR-0209: 省略すると、ADR-0205 の速い方法（軽い核種と重い核種の 400 kT までの自由ガス、細かさ 2、スイープ 3 回、H in H2O の
    // S(α, β)）。`ufg-free-gas` の書かない項目もその値
    let fast = FreeGas {
        upper: Some(400.0),
        cells: 2,
        sweeps: Some(3),
        heavy: Some(400.0),
        sab: true,
    };
    assert_eq!(FreeGas::default(), fast);
    let scattering = |lines: &str| {
        replaced(
            PIN_CELL,
            "library = \"m1.msgpack\"",
            &format!("library = \"m1.msgpack\"\n{lines}"),
        )
        .map(|input| input.shielding_options().scattering)
    };
    assert_eq!(
        LatticeInput::parse(PIN_CELL).unwrap().scattering,
        Scattering::FreeGas(fast)
    );
    assert_eq!(
        scattering("ufg-scattering = \"free-gas\"").unwrap(),
        Scattering::FreeGas(fast)
    );
    assert_eq!(
        scattering("ufg-scattering = \"static\"").unwrap(),
        Scattering::Static
    );
    // `ufg-scattering` を省略しても条件を書ける
    assert_eq!(
        scattering("ufg-free-gas = { cells = 4 }").unwrap(),
        Scattering::FreeGas(FreeGas { cells: 4, ..fast })
    );
    assert_eq!(
        scattering(
            "ufg-scattering = \"free-gas\"\nufg-free-gas = { upper = 300, cells = 8, sweeps = 2 }"
        )
        .unwrap(),
        Scattering::FreeGas(FreeGas {
            upper: Some(300.0),
            cells: 8,
            sweeps: Some(2),
            ..fast
        })
    );
    // ADR-0203 の参照の値（上限なし、収束まで、重い核種は静止した標的、自由ガスの H）も書ける: inf、sweeps = inf、heavy = false、sab = false
    assert_eq!(
        scattering(
            "ufg-free-gas = { upper = inf, cells = 4, sweeps = inf, heavy = false, sab = false }"
        )
        .unwrap(),
        Scattering::FreeGas(FreeGas {
            upper: None,
            cells: 4,
            sweeps: None,
            heavy: None,
            sab: false,
        })
    );
    // ADR-0204: `heavy = inf` は重い核種の上限なしの自由ガス
    assert_eq!(
        scattering("ufg-free-gas = { heavy = inf }").unwrap(),
        Scattering::FreeGas(FreeGas {
            heavy: Some(f64::INFINITY),
            ..fast
        })
    );
    for bad in [
        "ufg-scattering = \"moving\"",
        "ufg-scattering = \"static\"\nufg-free-gas = { cells = 2 }",
        "ufg-free-gas = { cells = 0 }",
        "ufg-free-gas = { sweeps = 0 }",
        "ufg-free-gas = { sweeps = 2.5 }",
        "ufg-free-gas = { sweeps = -inf }",
        "ufg-free-gas = { sweeps = nan }",
        "ufg-free-gas = { upper = 0 }",
        "ufg-free-gas = { upper = -400 }",
        "ufg-free-gas = { upper = nan }",
        "ufg-free-gas = { lower = 4 }",
        "ufg-free-gas = { heavy = 0 }",
        "ufg-free-gas = { heavy = -400 }",
        "ufg-free-gas = { heavy = nan }",
        "ufg-free-gas = { heavy = true }",
        "ufg-free-gas = { sab = 1 }",
    ] {
        assert!(scattering(bad).is_err(), "{bad}");
    }
    for off in [
        "ufg-scattering = \"free-gas\"",
        "ufg-scattering = \"static\"",
        "ufg-free-gas = { cells = 4 }",
    ] {
        let off = replaced(
            PIN_CELL,
            "library = \"m1.msgpack\"",
            &format!("library = \"m1.msgpack\"\nself-shielding = false\n{off}"),
        );
        assert!(matches!(off, Err(InputError::Invalid(_))), "{off:?}");
    }
}

#[test]
fn ufg_transfer_defaults_to_the_ring_flux_rows() {
    // ADR-0209: 省略すると、等価なセルのリングの UFG の中性子束で移行を数える（ADR-0207）。"infinite-dilution" でライブラリの
    // 無限希釈の形（ADR-0076）
    let transfer = |lines: &str| {
        replaced(
            PIN_CELL,
            "library = \"m1.msgpack\"",
            &format!("library = \"m1.msgpack\"\n{lines}"),
        )
        .map(|input| input.shielding_options().transfer)
    };
    assert_eq!(
        LatticeInput::parse(PIN_CELL)
            .unwrap()
            .shielding_options()
            .transfer,
        UfgTransfer::RingFlux
    );
    assert_eq!(
        transfer("ufg-transfer = \"infinite-dilution\"").unwrap(),
        UfgTransfer::InfiniteDilution
    );
    assert_eq!(
        transfer("ufg-transfer = \"ring-flux\"").unwrap(),
        UfgTransfer::RingFlux
    );
    // 間引いた格子とも組める（ADR-0208）
    assert_eq!(
        transfer("ufg-tolerance = 1e-4").unwrap(),
        UfgTransfer::RingFlux
    );
    for bad in ["ufg-transfer = \"ring\"", "ufg-transfer = true"] {
        assert!(transfer(bad).is_err(), "{bad}");
    }
    let off = replaced(
        PIN_CELL,
        "library = \"m1.msgpack\"",
        "library = \"m1.msgpack\"\nself-shielding = false\nufg-transfer = \"ring-flux\"",
    );
    assert!(matches!(off, Err(InputError::Invalid(_))), "{off:?}");
}

#[test]
fn static_targets_and_infinite_dilution_reproduce_the_former_defaults() {
    // ADR-0209: 静止した標的と無限希釈の形を書いた入力の自己遮蔽の条件は、既定を切り替える前の既定（ADR-0070、ADR-0076）の条件と同じ。
    // 解く側は条件だけで決まる
    let input = replaced(
        PIN_CELL,
        "library = \"m1.msgpack\"",
        "library = \"m1.msgpack\"\nufg-scattering = \"static\"\nufg-transfer = \"infinite-dilution\"",
    )
    .unwrap();
    assert_eq!(
        input.shielding_options(),
        ShieldingOptions {
            tolerance: 1.0e-6,
            max_iterations: 100,
            dancoff: crate::moc::dancoff::DancoffOptions::default(),
            collision: Collision::Table { density: 16 },
            thinning: None,
            hold_memory: None,
            scattering: Scattering::Static,
            transfer: UfgTransfer::InfiniteDilution,
        }
    );
}

#[test]
fn self_shielding_needs_a_library() {
    let result = replaced(
        LATTICE,
        "cross-sections = \"c5g7.toml\"",
        "cross-sections = \"c5g7.toml\"\nself-shielding = true",
    );
    assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
}

#[test]
fn scattering_order_is_optional() {
    // 省略すると、ライブラリにある次数をすべて使う（決めるのは物質を作った後）
    assert_eq!(
        LatticeInput::parse(PIN_CELL).unwrap().scattering_order,
        None
    );
    let p1 = replaced(
        PIN_CELL,
        "pitch = 1.26",
        "pitch = 1.26\nscattering-order = 1",
    )
    .unwrap();
    assert_eq!(p1.scattering_order, Some(1));
    // 巨視的断面積は P0 だけ
    let result = replaced(
        LATTICE,
        "cross-sections = \"c5g7.toml\"",
        "cross-sections = \"c5g7.toml\"\nscattering-order = 1",
    );
    assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
}

#[test]
fn parses_pin_cell() {
    let input = LatticeInput::parse(PIN_CELL).unwrap();
    assert_eq!(
        input.cross_sections,
        CrossSections::Library(PathBuf::from("m1.msgpack"))
    );
    assert_eq!(
        input.output,
        Some(FewGroupOutput {
            bounds: vec![1e7, 0.625, 1e-5],
            files: FewGroupFiles::Lattice(PathBuf::from("pincell.fg.msgpack")),
            pins: false,
            split: 1,
            axial: false,
        })
    );
    assert_eq!(input.pin_powers, None);
    assert_eq!(
        (input.pitch, input.azimuthal, input.spacing),
        (1.26, 16, 0.05)
    );
    assert_eq!(input.boundaries, Boundaries::REFLECTIVE);
    // 物質と核種は名前の順
    assert_eq!(input.materials.len(), 2);
    assert_eq!(input.materials[0].name, "fuel");
    assert_eq!(input.materials[0].temperature, 900.0);
    assert_eq!(
        input.materials[0].composition[1],
        ("U-235".to_string(), 7.0e-4)
    );
    assert_eq!(
        input.cells[0].zones,
        vec![
            Zone {
                material: "fuel".into(),
                outer: Some(0.4),
                rings: 2
            },
            Zone {
                material: "water".into(),
                outer: None,
                rings: 2
            },
        ]
    );
    assert_eq!((input.cells[0].grid, input.cells[0].sectors), (1, 1));
    assert_eq!(input.rows, vec![vec!["pin".to_string()]]);
}

#[test]
fn pin_cell_geometry_splits_fuel_by_volume_and_water_by_radius() {
    let input = LatticeInput::parse(PIN_CELL).unwrap();
    let (lattice, materials) = input.geometry().unwrap();
    assert_eq!(materials, vec!["fuel", "fuel", "water", "water"]);
    // 燃料は等体積の2つのリング、水はペレットとセルの内接円の中点で分ける
    let v = lattice.volumes();
    assert!((v[0] - v[1]).abs() < 1e-14);
    let r = (0.4 + 0.63) / 2.0;
    assert!((v[2] - std::f64::consts::PI * (r * r - 0.16)).abs() < 1e-14);
}

#[test]
fn zones_may_reach_beyond_the_inscribed_circle() {
    let geometry = |zones: &str| {
        let text = PIN_CELL.replace(
            "zones = [{ material = \"fuel\", radius = 0.4, rings = 2 }, { material = \"water\", rings = 2 }]",
            zones,
        );
        LatticeInput::parse(&text).and_then(|i| i.geometry())
    };
    // 内接円（0.63）の外の 0.8 までの水の帯は、セルの外周で切り取った2つのリングになる
    let (lattice, materials) = geometry(
        r#"zones = [
  { material = "fuel", radius = 0.4, rings = 2 },
  { material = "water", radius = 0.8, rings = 2 },
  { material = "water", rings = 1 },
]"#,
    )
    .unwrap();
    assert_eq!(materials.len(), 5);
    let v = lattice.volumes();
    assert!((v.iter().sum::<f64>() - 1.26 * 1.26).abs() < 1e-14);
    // 最後の帯は内接円の外から始まるので、等間隔には分けられない
    assert!(matches!(
        geometry(
            r#"zones = [
  { material = "fuel", radius = 0.4, rings = 2 },
  { material = "water", radius = 0.8, rings = 2 },
  { material = "water", rings = 2 },
]"#
        ),
        Err(InputError::LastZoneOutside { .. })
    ));
}

#[test]
fn parses_lattice() {
    let input = LatticeInput::parse(LATTICE).unwrap();
    assert_eq!(
        input.cross_sections,
        CrossSections::Macroscopic(PathBuf::from("c5g7.toml"))
    );
    assert_eq!(input.output, None);
    assert_eq!(input.pin_powers, Some(PathBuf::from("pins.msgpack")));
    assert_eq!(
        input.boundaries,
        Boundaries {
            right: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        }
    );
    let (lattice, materials) = input.geometry().unwrap();
    assert_eq!(lattice.width(), 4.0);
    // セルは配置の順（上の行から）
    assert_eq!(
        materials,
        vec![
            "UO2", "clad", "water", "water", "water", "UO2", "clad", "water"
        ]
    );
}

#[test]
fn columns_and_rows_may_have_their_own_sizes() {
    // ピンの格子の右に、厚さ 0.5 のバッフルと 1.5 の水の列を足す（反射体）
    let reflector = |rows: &str, extra: &str| {
        let text = LATTICE
            .replace(
                "rows = [\"pin box\", \"box pin\"]",
                &format!("widths = [{{ size = 2.0, count = 2 }}, {{ size = 0.5 }}, {{ size = 1.5 }}]\n{extra}rows = {rows}"),
            )
            .replace(
                "[cell.box]",
                "[cell.baffle]\nzones = [{ material = \"clad\", rings = 1 }]\n\n[cell.box]",
            );
        LatticeInput::parse(&text).and_then(|input| input.geometry())
    };
    let (lattice, materials) =
        reflector("[\"pin box baffle box\", \"box pin baffle box\"]", "").unwrap();
    assert_eq!(lattice.widths(), &[2.0, 2.0, 0.5, 1.5]);
    assert_eq!(lattice.heights(), &[2.0, 2.0]);
    // ピン 3 + 水 1 + バッフル 1 + 水 1 の行が 2 つ
    assert_eq!(lattice.regions(), 12);
    assert_eq!(materials[4], "clad");
    // 均質なセルは、ピッチのマスではピンセル、ほかのマスでは長方形のセル
    assert!(lattice.cell_type(1).pin().is_some());
    assert!(lattice.cell_type(3).pin().is_none());
    assert_eq!(lattice.pin_rows_and_columns(), (vec![0, 1], vec![0, 1]));
    // ピンはピッチの幅の列にだけ置ける
    assert!(reflector("[\"pin box pin box\", \"box pin baffle box\"]", "").is_err());
    // 列の数がそろわない、余白と同時には使えない
    assert!(reflector("[\"pin box baffle\", \"box pin baffle\"]", "").is_err());
    assert!(
        reflector(
            "[\"pin box baffle box\", \"box pin baffle box\"]",
            "gap = { width = 0.1, material = \"water\" }\n"
        )
        .is_err()
    );
}

#[test]
fn sides_may_have_an_albedo() {
    let input = replaced(
        LATTICE,
        "boundary = { right = \"vacuum\", bottom = \"vacuum\" }",
        "boundary = { right = \"vacuum\", top = { albedo = 0.3 } }",
    )
    .unwrap();
    assert_eq!(input.boundaries.top, Boundary::Albedo(0.3));
    assert_eq!(input.boundaries.right, Boundary::Vacuum);
    assert_eq!(input.boundaries.bottom, Boundary::Reflective);
    assert!(
        replaced(
            LATTICE,
            "boundary = { right = \"vacuum\", bottom = \"vacuum\" }",
            "boundary = { top = { albedo = 0.3, extra = 1 } }",
        )
        .is_err()
    );
    // 書き間違いには、書ける条件を示す
    let error = replaced(
        LATTICE,
        "boundary = { right = \"vacuum\", bottom = \"vacuum\" }",
        "boundary = { right = \"vacum\" }",
    )
    .unwrap_err()
    .to_string();
    assert!(
        error.contains("reflective") && error.contains("albedo"),
        "{error}"
    );
}

#[test]
fn gap_surrounds_the_lattice() {
    let text = LATTICE.replace(
        "rows = [",
        "gap = { width = 0.25, material = \"water\" }\nrows = [",
    );
    let (lattice, materials) = LatticeInput::parse(&text).unwrap().geometry().unwrap();
    assert_eq!(lattice.width(), 4.5);
    assert_eq!(lattice.widths(), &[0.25, 2.0, 2.0, 0.25]);
    // 余白は周りの薄い行と列のセル（各 1 領域）で、領域はセルの順に並ぶ
    assert_eq!(materials.len(), lattice.regions());
    assert_eq!(materials.len(), 8 + 12);
    assert!(materials[..4].iter().all(|m| m == "water"));
    assert_eq!(materials[4..9], ["water", "UO2", "clad", "water", "water"]);
    assert_eq!(materials.iter().filter(|m| *m == "UO2").count(), 2);
    // 幅が正でない、物質がない、真空（巨視的断面積でも）
    for gap in [
        "{ width = 0.0, material = \"water\" }",
        "{ width = 0.1 }",
        "{ width = 0.1, material = \"void\" }",
    ] {
        let text = LATTICE.replace("rows = [", &format!("gap = {gap}\nrows = ["));
        let result = LatticeInput::parse(&text).and_then(|input| input.geometry().map(|_| ()));
        assert!(result.is_err(), "{gap}");
    }
}

#[test]
fn frames_are_cell_edges() {
    let with_frame = |frame: &str| {
        replaced(
            LATTICE,
            "  { material = \"water\", rings = 1 },\n]",
            &format!("  {{ material = \"water\", rings = 2 }},\n]\nframe = {frame}"),
        )
    };
    let input = with_frame("{ thickness = 0.05, material = \"clad\" }").unwrap();
    let (lattice, regions) = input.geometry().unwrap();
    // ピンセル 4 + 枠 1 が 2 つ、均質なセル 1 が 2 つ
    assert_eq!(lattice.regions(), 2 * 5 + 2);
    assert_eq!(regions[4], "clad");
    // 厚すぎる枠は、帯の配置ではなく枠の誤りとして報告する
    for thickness in ["1.0", "-0.1", "nan"] {
        let frame = format!("{{ thickness = {thickness}, material = \"clad\" }}");
        let error = with_frame(&frame).and_then(|input| input.geometry().map(|_| ()));
        assert!(
            matches!(
                error,
                Err(InputError::Cell {
                    source: GeometryError::Frame(_),
                    ..
                })
            ),
            "{thickness}: {error:?}"
        );
    }
    assert!(with_frame("{ thickness = 0.05 }").is_err());
}

#[test]
fn grid_and_sectors_split_cells() {
    let text = LATTICE.replace(
        "[cell.box]\nzones = [{ material = \"water\", rings = 1 }]",
        "[cell.box]\nzones = [{ material = \"water\", rings = 1 }]\ngrid = 3",
    );
    let input = LatticeInput::parse(&text).unwrap();
    let (lattice, _) = input.geometry().unwrap();
    assert_eq!(lattice.regions(), 3 + 9 + 9 + 3);
    let text = LATTICE.replace(
        "  { material = \"water\", rings = 1 },\n]",
        "  { material = \"water\", rings = 1 },\n]\nsectors = 8",
    );
    let (lattice, _) = LatticeInput::parse(&text).unwrap().geometry().unwrap();
    assert_eq!(lattice.regions(), 2 * 24 + 2);
}

#[test]
fn invalid_inputs_are_errors() {
    // 断面積の出どころは library と cross-sections のどちらか一方
    assert!(LatticeInput::parse(&format!("cross-sections = \"a.toml\"\n{PIN_CELL}")).is_err());
    assert!(replaced(PIN_CELL, "library = \"m1.msgpack\"\n", "").is_err());
    // 巨視的断面積では material を書かない、少数群定数を出さない
    assert!(
        LatticeInput::parse(&format!(
            "{LATTICE}\n[material.UO2]\ntemperature = 600\nnuclides = {{ U-235 = 1e-3 }}"
        ))
        .is_err()
    );
    assert!(
        replaced(
            LATTICE,
            "pitch = 2.0",
            "pitch = 2.0\nfew-group = { output = \"a\", bounds = [1e7, 1e-5] }"
        )
        .is_err()
    );
    // 知らないキー、知らない辺と条件
    assert!(replaced(PIN_CELL, "pitch = 1.26", "pitch = 1.26\nalbedo = 1").is_err());
    assert!(replaced(LATTICE, "right = \"vacuum\"", "east = \"vacuum\"").is_err());
    assert!(replaced(LATTICE, "right = \"vacuum\"", "right = \"white\"").is_err());
    // 正でない値、温度、個数密度、核種名
    assert!(replaced(PIN_CELL, "pitch = 1.26", "pitch = 0").is_err());
    assert!(replaced(PIN_CELL, "spacing = 0.05", "spacing = -0.05").is_err());
    assert!(replaced(PIN_CELL, "azimuthal = 16", "azimuthal = 0").is_err());
    assert!(replaced(PIN_CELL, "temperature = 900", "temperature = -900").is_err());
    assert!(replaced(PIN_CELL, "U-238 = 2.2e-2", "U-238 = -2.2e-2").is_err());
    assert!(replaced(PIN_CELL, "U-238 = 2.2e-2", "U238 = 2.2e-2").is_err());
    assert!(
        replaced(
            PIN_CELL,
            "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
            "nuclides = {}"
        )
        .is_err()
    );
    // 帯: 外半径が昇順でない、最後の帯に外半径がある、途中の帯に外半径がない、分割数が 0
    assert!(replaced(LATTICE, "radius = 0.6", "radius = 0.4").is_err());
    assert!(
        replaced(
            PIN_CELL,
            "{ material = \"water\", rings = 2 }",
            "{ material = \"water\", radius = 0.6, rings = 2 }"
        )
        .is_err()
    );
    assert!(
        replaced(
            PIN_CELL,
            "{ material = \"fuel\", radius = 0.4, rings = 2 }",
            "{ material = \"fuel\", rings = 2 }"
        )
        .is_err()
    );
    assert!(replaced(PIN_CELL, "rings = 2 }]", "rings = 0 }]").is_err());
    // grid は均質なセル（帯が1つで分割数 1）だけ、grid と sectors は同時に使わない
    assert!(replaced(PIN_CELL, "rings = 2 }]", "rings = 2 }]\ngrid = 2").is_err());
    assert!(
        replaced(
            LATTICE,
            "zones = [{ material = \"water\", rings = 1 }]",
            "zones = [{ material = \"water\", rings = 1 }]\ngrid = 2\nsectors = 4"
        )
        .is_err()
    );
    assert!(
        replaced(
            LATTICE,
            "zones = [{ material = \"water\", rings = 1 }]",
            "zones = [{ material = \"water\", rings = 1 }]\ngrid = 0"
        )
        .is_err()
    );
}

#[test]
fn geometry_checks_names_and_layout() {
    let geometry = |text: String| LatticeInput::parse(&text).and_then(|i| i.geometry());
    // 配置に知らないセル、行の長さが行の数と違う
    assert!(geometry(LATTICE.replace("\"box pin\"", "\"box rod\"")).is_err());
    assert!(geometry(LATTICE.replace("\"box pin\"", "\"box pin pin\"")).is_err());
    // ライブラリのときは、セルの物質が material で定義されていること
    assert!(
        geometry(PIN_CELL.replace(
            "{ material = \"water\", rings = 2 }",
            "{ material = \"steam\", rings = 2 }"
        ))
        .is_err()
    );
    // 帯の外半径がセルの外接円（ピッチ 2 で √2）を超える、扇形の数が 1, 4, 8 でない
    assert!(geometry(LATTICE.replace("radius = 0.6", "radius = 1.5")).is_err());
    assert!(
        geometry(LATTICE.replace(
            "  { material = \"water\", rings = 1 },\n]",
            "  { material = \"water\", rings = 1 },\n]\nsectors = 3"
        ))
        .is_err()
    );
}

#[test]
fn cell_errors_keep_the_cause_in_the_chain_only() {
    // 原因の文は連鎖（source）にだけ置く。自分の文にも埋め込むと、連鎖を表示するときに2回出る
    let text = LATTICE.replace("radius = 0.6", "radius = 1.5");
    let err = LatticeInput::parse(&text).unwrap().geometry().unwrap_err();
    let cause = std::error::Error::source(&err).expect("the geometry error is the cause");
    assert!(!err.to_string().contains(&cause.to_string()), "{err}");
}

#[test]
fn toml_errors_are_not_repeated_in_the_chain() {
    // 包んだ TOML のエラーは、自分の文と原因の連鎖のどちらか一方にだけ出す
    let err = LatticeInput::parse("pitch = ").unwrap_err();
    let chain = format!("{:#}", anyhow::Error::from(err));
    assert_eq!(chain.matches("TOML parse error").count(), 1, "{chain}");
}

#[test]
fn parses_few_group_nodes() {
    // 格子をノードに分け、ノードごとに少数群定数を書く（ADR-0095）。列と行は格子のセル（余白の薄い行と
    // 列を含む）の、1 から数えた最初と最後。LATTICE は巨視的断面積なので、ライブラリに置き換える
    let text = LATTICE.replace("cross-sections = \"c5g7.toml\"\n", "");
    let input = replaced(
        &text,
        "pin-powers = \"pins.msgpack\"",
        "library = \"m1.msgpack\"\nfew-group = { bounds = [1e7, 1e-5], nodes = [\n  \
             { output = \"fuel.msgpack\", columns = [1, 1], rows = [1, 2] },\n  \
             { output = \"reflector.msgpack\", columns = [2, 2], rows = [1, 2], reflector = true },\n] }",
    );
    let input = input.unwrap_or_else(|e| panic!("{e}"));
    let expected = FewGroupOutput {
        bounds: vec![1e7, 1e-5],
        files: FewGroupFiles::Nodes(vec![
            (
                PathBuf::from("fuel.msgpack"),
                Node {
                    columns: 0..1,
                    rows: 0..2,
                    reflector: false,
                },
            ),
            (
                PathBuf::from("reflector.msgpack"),
                Node {
                    columns: 1..2,
                    rows: 0..2,
                    reflector: true,
                },
            ),
        ]),
        pins: false,
        split: 1,
        axial: false,
    };
    assert_eq!(input.output, Some(expected));
}

#[test]
fn few_group_may_ask_for_pin_form_functions() {
    // ピン出力の再構成の形状関数と角の不連続因子を付ける（ADR-0096）
    let few_group = |entry: &str| {
        replaced(
            PIN_CELL,
            "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
            &format!("few-group = {{ bounds = [1e7, 1e-5]{entry} }}"),
        )
    };
    let input = few_group(", output = \"a\", pins = true").unwrap();
    assert!(input.output.is_some_and(|o| o.pins));
    let input = few_group(", output = \"a\", pins = false").unwrap();
    assert!(input.output.is_some_and(|o| !o.pins));
    // 格子のノードごとの定数にも形状関数と角の不連続因子を付けられる（ADR-0131）。2 群で、ノードが格子の上の
    // 長方形の並びであることを、MOC を解く前に確かめる
    let with_nodes = |bounds: &str, nodes: &str| {
        replaced(
            PIN_CELL,
            "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
            &format!("few-group = {{ bounds = {bounds}, nodes = [{nodes}], pins = true }}"),
        )
    };
    let one = "{ output = \"a\", columns = [1, 1], rows = [1, 1] }";
    let input = with_nodes("[1e7, 0.625, 1e-5]", one).unwrap();
    assert!(input.output.is_some_and(|o| o.pins));
    let err = with_nodes("[1e7, 1e-5]", one).unwrap_err().to_string();
    assert!(err.contains("2 few groups"), "{err}");
    // 左のノードが 2 行にまたがり、右は 1 行ずつ: 列と行の切れ目の長方形の並びにならない
    let uneven = "{ output = \"a\", columns = [1, 1], rows = [1, 2] }, \
                      { output = \"b\", columns = [2, 2], rows = [1, 1] }, \
                      { output = \"c\", columns = [2, 2], rows = [2, 2] }";
    let err = with_nodes("[1e7, 0.625, 1e-5]", uneven)
        .unwrap_err()
        .to_string();
    assert!(err.contains("grid"), "{err}");
}

#[test]
fn few_group_may_be_axial() {
    // 格子の x を炉心の z として定数を書く（ADR-0128）。上下の辺は鏡面反射で、形状関数と部分ノードは作らない
    let slab = |entry: &str, boundary: &str| {
        replaced(
            PIN_CELL,
            "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
            &format!(
                "few-group = {{ bounds = [1e7, 1e-5], nodes = [{{ output = \"a\", columns = [1, 1], rows = [1, 1] }}]{entry} }}\nboundary = {boundary}"
            ),
        )
    };
    let open = "{ left = \"vacuum\", right = \"vacuum\" }";
    let input = slab(", axial = true", open).unwrap();
    assert!(input.output.unwrap().axial);
    assert!(!slab("", open).unwrap().output.unwrap().axial);
    let err = slab(", axial = true", "{ top = \"vacuum\" }")
        .unwrap_err()
        .to_string();
    assert!(err.contains("axial") && err.contains("reflective"), "{err}");
    let err = slab(", axial = true, pins = true", open)
        .unwrap_err()
        .to_string();
    assert!(err.contains("axial") && err.contains("pins"), "{err}");
}

#[test]
fn few_group_may_split_the_lattice() {
    // 格子全体を n × n の部分ノードに分けた定数（ADR-0108）。省略すると 1。nodes とは使えない
    let few_group = |entry: &str| {
        replaced(
            PIN_CELL,
            "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
            &format!("few-group = {{ bounds = [1e7, 1e-5]{entry} }}"),
        )
    };
    let input = few_group(", output = \"a\", split = 2").unwrap();
    assert_eq!(input.output.map(|o| o.split), Some(2));
    let input = few_group(", output = \"a\"").unwrap();
    assert_eq!(input.output.map(|o| o.split), Some(1));
    assert!(few_group(", output = \"a\", split = 0").is_err());
    let nodes = ", nodes = [{ output = \"a\", columns = [1, 1], rows = [1, 1] }], split = 2";
    let err = few_group(nodes).unwrap_err().to_string();
    assert!(err.contains("split"), "{err}");
    // 外周がすべて鏡面反射でない格子は、MOC を解く前に入力の段階で止める（pins も同じ）
    for entry in [
        ", output = \"a\", split = 2",
        ", output = \"a\", pins = true",
    ] {
        let open = replaced(
            PIN_CELL,
            "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
            &format!(
                "few-group = {{ bounds = [1e7, 1e-5]{entry} }}\nboundary = {{ right = \"vacuum\" }}"
            ),
        );
        let err = open.unwrap_err().to_string();
        assert!(err.contains("reflective"), "{entry}: {err}");
    }
}

#[test]
fn invalid_few_group_nodes_are_errors() {
    let few_group = |entry: &str| {
        replaced(
            PIN_CELL,
            "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
            &format!("few-group = {{ bounds = [1e7, 1e-5]{entry} }}"),
        )
    };
    let node = |output: &str, columns: &str, rows: &str| {
        format!("{{ output = \"{output}\", columns = {columns}, rows = {rows} }}")
    };
    let good = node("a", "[1, 1]", "[1, 1]");
    assert!(few_group(&format!(", nodes = [{good}]")).is_ok());
    for (entry, what) in [
        (String::new(), "neither output nor nodes"),
        (
            format!(", output = \"x\", nodes = [{good}]"),
            "both output and nodes",
        ),
        (", nodes = []".to_string(), "no nodes"),
        (
            format!(", nodes = [{}]", node("a", "[0, 1]", "[1, 1]")),
            "column 0",
        ),
        (
            format!(", nodes = [{}]", node("a", "[2, 1]", "[1, 1]")),
            "reversed columns",
        ),
        (
            format!(", nodes = [{}]", node("a", "[1, 1]", "[3, 2]")),
            "reversed rows",
        ),
        (
            format!(", nodes = [{good}, {}]", node("a", "[2, 2]", "[1, 1]")),
            "same output",
        ),
        (
            format!(", nodes = [{}]", node("a", "[1, 2, 3]", "[1, 1]")),
            "three numbers",
        ),
    ] {
        assert!(
            matches!(
                few_group(&entry),
                Err(InputError::Invalid(_) | InputError::Toml(_))
            ),
            "{what}"
        );
    }
}

#[test]
fn borated_water_matches_vera_number_densities() {
    // VERA の Problem 5 の減速材（565 K、0.743 g/cm³、B-10 19.9 at%）の個数密度の表（ADR-0106）
    let table = [
        (1285.0, [4.96231e-2, 1.05835e-5, 4.25999e-5, 2.48116e-2]),
        (1291.0, [4.96228e-2, 1.06329e-5, 4.27988e-5, 2.48114e-2]),
        (1170.0, [4.96288e-2, 9.63633e-6, 3.87874e-5, 2.48144e-2]),
        (1230.0, [4.96258e-2, 1.01305e-5, 4.07765e-5, 2.48129e-2]),
    ];
    for (ppm, expected) in table {
        let input = replaced(
            PIN_CELL,
            "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
            &format!("water = {{ density = 0.743, boron = {ppm}, b10 = 0.199 }}"),
        )
        .unwrap();
        let water = input.materials.iter().find(|m| m.name == "water").unwrap();
        let names: Vec<&str> = water.composition.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(names, ["B-10", "B-11", "H-1", "O-16"]);
        let [h, b10, b11, o] = expected;
        for ((_, got), want) in water.composition.iter().zip([b10, b11, h, o]) {
            // 表は有効数字 6 桁
            assert!(
                ((got - want) / want).abs() < 1e-5,
                "{ppm} ppm: {got} vs {want}"
            );
        }
    }
}

#[test]
fn water_without_boron_has_no_boron_nuclides() {
    let input = replaced(
        PIN_CELL,
        "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
        "water = { density = 0.743, boron = 0, b10 = 0.199 }",
    )
    .unwrap();
    let water = input.materials.iter().find(|m| m.name == "water").unwrap();
    let names: Vec<&str> = water.composition.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(names, ["H-1", "O-16"]);
    assert!(water.composition.iter().all(|&(_, n)| n > 0.0));
}

#[test]
fn water_keeps_only_present_boron_isotopes() {
    for (b10, names) in [
        (1.0, ["B-10", "H-1", "O-16"]),
        (0.0, ["B-11", "H-1", "O-16"]),
    ] {
        let input = replaced(
            PIN_CELL,
            "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
            &format!("water = {{ density = 0.743, boron = 1291, b10 = {b10:?} }}"),
        )
        .unwrap();
        let water = input.materials.iter().find(|m| m.name == "water").unwrap();
        let got: Vec<&str> = water.composition.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(got, names);
    }
}

#[test]
fn invalid_water_is_an_error() {
    let water = "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }";
    let bad = |to: &str| replaced(PIN_CELL, water, to).is_err();
    assert!(bad("water = { density = 0.0, boron = 1291, b10 = 0.199 }"));
    assert!(bad("water = { density = 0.743, boron = -1, b10 = 0.199 }"));
    assert!(bad("water = { density = 0.743, boron = 1e6, b10 = 0.199 }"));
    assert!(bad("water = { density = 0.743, boron = 1291, b10 = 1.5 }"));
    // 3つとも必須、nuclides と同時には書けない
    assert!(bad("water = { density = 0.743, boron = 1291 }"));
    assert!(bad(&format!(
        "{water}\nwater = {{ density = 0.743, boron = 1291, b10 = 0.199 }}"
    )));
    assert!(bad(""));
}

/// 燃焼（ADR-0173）の節を足したピンセル
fn depleting(depletion: &str) -> Result<LatticeInput, InputError> {
    LatticeInput::parse(&format!(
        "{}\n[depletion]\n{depletion}",
        without_few_group()
    ))
}

/// 格子全体の少数群定数を書かないピンセル（燃焼は少数群定数を燃焼の履歴に書く）
fn without_few_group() -> String {
    PIN_CELL.replace(
        "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }\n",
        "",
    )
}

const DEPLETION: &str = r#"chain = "chain.msgpack"
power = 40.0
burnup = [0.1, 0.5, 2]
materials = ["fuel"]
output = "history.msgpack"
"#;

#[test]
fn reads_the_depletion_steps_and_materials() {
    let input = depleting(DEPLETION).unwrap();
    let depletion = input.depletion.unwrap();
    assert_eq!(depletion.chain, PathBuf::from("chain.msgpack"));
    assert_eq!(depletion.power, 40.0);
    assert_eq!(depletion.burnup, vec![0.1, 0.5, 2.0]);
    assert_eq!(depletion.materials, vec!["fuel".to_string()]);
    assert_eq!(depletion.output, PathBuf::from("history.msgpack"));
    assert!(LatticeInput::parse(PIN_CELL).unwrap().depletion.is_none());
}

#[test]
fn depletion_shields_at_every_transport_unless_told_to_shield_at_the_step_start() {
    let every = depleting(DEPLETION).unwrap().depletion.unwrap();
    assert_eq!(every.shielding, ShieldingFrequency::EveryTransport);
    let text = format!("{DEPLETION}shielding = \"every-transport\"\n");
    let written = depleting(&text).unwrap().depletion.unwrap();
    assert_eq!(written.shielding, ShieldingFrequency::EveryTransport);
    let text = format!("{DEPLETION}shielding = \"step-start\"\n");
    let start = depleting(&text).unwrap().depletion.unwrap();
    assert_eq!(start.shielding, ShieldingFrequency::StepStart);
    let text = format!("{DEPLETION}shielding = \"midpoint\"\n");
    assert!(depleting(&text).is_err());
}

#[test]
fn rejects_burnup_points_that_do_not_increase_from_zero() {
    for burnup in ["[]", "[0, 1]", "[1, 1]", "[2, 1]", "[-1]"] {
        let text = DEPLETION.replace("[0.1, 0.5, 2]", burnup);
        assert!(depleting(&text).is_err(), "{burnup}");
    }
}

#[test]
fn rejects_a_non_positive_power() {
    for power in ["0.0", "-40.0", "nan"] {
        let text = DEPLETION.replace("40.0", power);
        assert!(depleting(&text).is_err(), "{power}");
    }
}

#[test]
fn depleted_materials_are_defined_nuclide_materials() {
    for materials in [r#"[]"#, r#"["clad"]"#, r#"["fuel", "fuel"]"#] {
        let text = DEPLETION.replace(r#"["fuel"]"#, materials);
        assert!(depleting(&text).is_err(), "{materials}");
    }
    // ホウ酸水で書いた物質は燃焼させない（核分裂する核種がない）ので、名前を挙げてもよいかは燃焼の領域で決まる。ここでは
    // 核種で書いた物質であることだけを見る
    let text = without_few_group().replace(
        "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
        "water = { density = 0.7, boron = 1000, b10 = 0.2 }",
    );
    let text = format!(
        "{text}\n[depletion]\n{}",
        DEPLETION.replace(r#"["fuel"]"#, r#"["fuel", "water"]"#)
    );
    assert!(LatticeInput::parse(&text).is_ok());
}

#[test]
fn depletion_excludes_branches_pin_powers_media_and_file_constants() {
    let base = without_few_group();
    let branches = format!(
        "{base}\n[branches]\nfuel = {{ temperature = [600, 900], materials = [\"fuel\"] }}\n[depletion]\n{DEPLETION}"
    );
    assert!(LatticeInput::parse(&branches).is_err());
    let pins = format!("pin-powers = \"pins.msgpack\"\n{base}\n[depletion]\n{DEPLETION}");
    assert!(LatticeInput::parse(&pins).is_err());
    let medium = format!("medium = \"fuel.medium\"\n{base}\n[depletion]\n{DEPLETION}");
    assert!(LatticeInput::parse(&medium).is_err());
    // 少数群定数のファイルは書かない（燃焼点ごとの定数は燃焼の履歴に入る。`[depletion]` の `few-group`）
    assert!(LatticeInput::parse(&format!("{PIN_CELL}\n[depletion]\n{DEPLETION}")).is_err());
    let unshielded = format!("self-shielding = false\n{base}\n[depletion]\n{DEPLETION}");
    assert!(LatticeInput::parse(&unshielded).is_err());
    // 燃焼点ごとの形状関数は、外周がすべて鏡面反射の格子でだけ
    let vacuum = base.replace(
        "spacing = 0.05",
        "spacing = 0.05\nboundary = { right = \"vacuum\" }",
    );
    let pins = format!("{DEPLETION}few-group = {{ bounds = [1e7, 0.625, 1e-5], pins = true }}\n");
    assert!(LatticeInput::parse(&format!("{vacuum}\n[depletion]\n{pins}")).is_err());
    assert!(LatticeInput::parse(&format!("{base}\n[depletion]\n{pins}")).is_ok());
}

#[test]
fn depletion_takes_the_few_group_bounds_in_its_own_table() {
    let text = format!("{DEPLETION}few-group = {{ bounds = [1e7, 0.625, 1e-5], pins = true }}\n");
    let depletion = depleting(&text).unwrap().depletion.unwrap();
    let few = depletion.few_group.unwrap();
    assert_eq!(few.bounds, vec![1e7, 0.625, 1e-5]);
    assert!(few.pins);
    assert!(
        depleting(DEPLETION)
            .unwrap()
            .depletion
            .unwrap()
            .few_group
            .is_none()
    );
    let bad = format!("{DEPLETION}few-group = {{ bounds = [1e7] }}\n");
    assert!(depleting(&bad).is_err());
}

#[test]
fn depletion_rejects_a_yield_cutoff_entry() {
    // 収率の重みの境界と点は ADR-0183 で固定し、入力では選ばない
    let text =
        format!("{DEPLETION}yield-cutoff = {{ energy = 112.5, thermal = 0.0253, fast = 5e5 }}\n");
    let err = depleting(&text).unwrap_err().to_string();
    assert!(err.contains("yield-cutoff"), "{err}");
}

#[test]
fn depleted_material_names_leave_room_for_the_region_names() {
    // 燃焼の領域の物質は `<物質>@<セル>` の名前で解く（ADR-0173）ので、`@` を含む名前の物質は燃焼させられない
    let text = without_few_group()
        .replace("[material.fuel]", "[material.\"fuel@0\"]")
        .replace("material = \"fuel\"", "material = \"fuel@0\"");
    let text = format!(
        "{text}\n[depletion]\n{}",
        DEPLETION.replace(r#"["fuel"]"#, r#"["fuel@0"]"#)
    );
    assert!(LatticeInput::parse(&text).is_err());
}

/// 燃焼の履歴の燃焼点からの分岐（ADR-0187）のピンセル: 減速材をホウ酸水で書き、3 つの軸の分岐と `[history]` を足す
fn history_branches(history: &str) -> String {
    let water = PIN_CELL.replace(
        "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
        "water = { density = 0.7, boron = 1000, b10 = 0.2 }",
    );
    format!(
        "{water}\n[branches]\nboron = {{ ppm = [500, 1000], materials = [\"water\"] }}\n\
         moderator = {{ temperature = [600], density = [0.7], materials = [\"water\"] }}\n\
         fuel = {{ temperature = [900], materials = [\"fuel\"] }}\n[history]\n{history}"
    )
}

const HISTORY: &str = r#"file = "history.msgpack"
burnup = [0, 10, 20]
"#;

#[test]
fn reads_the_history_file_and_its_burnup_points() {
    let input = LatticeInput::parse(&history_branches(HISTORY)).unwrap();
    let history = input.history.unwrap();
    assert_eq!(history.file, PathBuf::from("history.msgpack"));
    assert_eq!(history.burnup, vec![0.0, 10.0, 20.0]);
    assert!(LatticeInput::parse(PIN_CELL).unwrap().history.is_none());
}

#[test]
fn history_burnup_points_are_nonnegative_and_increase() {
    for burnup in ["[]", "[-1]", "[1, 1]", "[2, 1]", "[0, nan]"] {
        let text = history_branches(&HISTORY.replace("[0, 10, 20]", burnup));
        assert!(LatticeInput::parse(&text).is_err(), "{burnup}");
    }
    let text = history_branches(&HISTORY.replace("[0, 10, 20]", "[5]"));
    assert!(LatticeInput::parse(&text).is_ok());
}

#[test]
fn history_branches_need_branches_shielding_and_a_lattice_table() {
    let ok = history_branches(HISTORY);
    assert!(LatticeInput::parse(&ok).is_ok());
    // 分岐の軸がない（core は 3 つの軸を求める）
    let start = ok.find("[branches]").unwrap();
    let end = ok.find("[history]").unwrap();
    let unbranched = format!("{}{}", &ok[..start], &ok[end..]);
    assert!(LatticeInput::parse(&unbranched).is_err());
    // テーブルのファイルがない、ノードごとのファイル
    let no_table = ok.replace(
        "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }\n",
        "",
    );
    assert!(LatticeInput::parse(&no_table).is_err());
    let nodes = ok.replace(
        "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
        "few-group = { bounds = [1e7, 0.625, 1e-5], nodes = [{ output = \"a.msgpack\", columns = [1, 1], rows = [1, 1] }] }",
    );
    assert!(LatticeInput::parse(&nodes).is_err());
    // 自己遮蔽なし、媒質の出力、燃焼
    assert!(LatticeInput::parse(&format!("self-shielding = false\n{ok}")).is_err());
    assert!(LatticeInput::parse(&format!("medium = \"fuel.medium\"\n{ok}")).is_err());
    let depleting = format!("{ok}\n[depletion]\n{DEPLETION}");
    assert!(LatticeInput::parse(&depleting).is_err());
}

#[test]
fn history_micro_nuclides_are_distinct_names_of_an_unsplit_table() {
    // 書かなければミクロ燃焼の微視的断面積を作らない（ADR-0216）
    let input = LatticeInput::parse(&history_branches(HISTORY)).unwrap();
    assert!(input.history.unwrap().micro.is_empty());
    let micro = |list: &str| history_branches(&format!("{HISTORY}micro = {list}\n"));
    let input = LatticeInput::parse(&micro(r#"["U-235", "Xe-135"]"#)).unwrap();
    assert_eq!(input.history.unwrap().micro, vec!["U-235", "Xe-135"]);
    for bad in ["[]", r#"["U-235", "U-235"]"#, r#"["uranium"]"#] {
        assert!(LatticeInput::parse(&micro(bad)).is_err(), "{bad}");
    }
    // 格子全体で均質化するので、部分ノードに分けたテーブルとは組めない
    let split = micro(r#"["U-235"]"#).replace(
        "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
        "few-group = { output = \"pincell.fg.msgpack\", bounds = [1e7, 0.625, 1e-5], split = 2 }",
    );
    let err = LatticeInput::parse(&split).unwrap_err();
    assert!(err.to_string().contains("split"), "{err}");
}

#[test]
fn branches_shield_at_the_nominal_boron_unless_told_to_shield_every_state() {
    // 分岐の自己遮蔽（ADR-0201）。書かなければ、ホウ素だけを変えた状態は公称のホウ素の自己遮蔽を使う。履歴のない分岐も同じ
    let with = |value: &str| {
        history_branches(HISTORY).replace(
            "[history]\n",
            &format!("shielding = \"{value}\"\n[history]\n"),
        )
    };
    let plain = |text: String| text[..text.find("[history]").unwrap()].to_string();
    let shielding = |text: &str| {
        LatticeInput::parse(text)
            .map(|input| input.branches.unwrap().shielding())
            .map_err(|e| e.to_string())
    };
    for text in [history_branches(HISTORY), plain(history_branches(HISTORY))] {
        assert_eq!(shielding(&text), Ok(BranchShielding::NominalBoron));
    }
    assert_eq!(
        shielding(&with("nominal-boron")),
        Ok(BranchShielding::NominalBoron)
    );
    for text in [with("every-state"), plain(with("every-state"))] {
        assert_eq!(shielding(&text), Ok(BranchShielding::EveryState));
    }
    // 値は nominal-boron と every-state だけ
    assert!(shielding(&with("nominal-moderator")).is_err());
    assert!(shielding(&with("nominal")).is_err());
}

#[test]
fn off_nominal_burnups_are_history_burnups_and_default_to_all() {
    // 公称でない状態を解く燃焼点（ADR-0200 の案 B2）。書かなければすべての燃焼点（案 B1）
    let input = LatticeInput::parse(&history_branches(HISTORY)).unwrap();
    assert_eq!(input.history.unwrap().off_nominal, vec![0.0, 10.0, 20.0]);
    let text = history_branches(&format!("{HISTORY}off-nominal-burnup = [0, 20]\n"));
    let history = LatticeInput::parse(&text).unwrap().history.unwrap();
    assert_eq!(history.off_nominal, vec![0.0, 20.0]);
    // 公称でない状態の差は外挿しないので、最初と最後の燃焼点を含む
    for bad in ["[]", "[5]", "[20, 0]", "[0, 0]", "[0, 10]", "[10, 20]"] {
        let text = history_branches(&format!("{HISTORY}off-nominal-burnup = {bad}\n"));
        assert!(LatticeInput::parse(&text).is_err(), "{bad}");
    }
}

/// `history_branches` の減速材と燃料を 2 点ずつにしたもの（入力の物質は 600 K と 900 K）
fn three_axes(history: &str) -> String {
    history_branches(history)
        .replace(
            "moderator = { temperature = [600], density = [0.7], materials = [\"water\"] }",
            "moderator = { temperature = [600, 650], density = [0.7, 0.66], materials = [\"water\"] }",
        )
        .replace(
            "fuel = { temperature = [900], materials = [\"fuel\"] }",
            "fuel = { temperature = [900, 1200], materials = [\"fuel\"] }",
        )
}

#[test]
fn history_branches_default_to_all_pair_subgrids() {
    // 燃焼度の軸を持つ分岐のテーブルの既定は、2 軸の組すべての部分格子（ADR-0200）。2 × 2 × 2 で 3 つの軸がどれも公称でない
    // 状態だけを解かない
    let input = LatticeInput::parse(&three_axes(HISTORY)).unwrap();
    let branches = input.branches.as_ref().unwrap();
    assert_eq!(branches.nominal(), Some(&[1, 0, 0][..]));
    assert_eq!((0..8).filter(|&s| branches.solved(s)).count(), 7);
    assert!(!branches.solved(3));
    // 3 軸の部分格子を書けば直積
    let product = three_axes(HISTORY).replace(
        "[history]",
        "subgrids = [[\"boron\", \"moderator\", \"fuel\"]]\n[history]",
    );
    let input = LatticeInput::parse(&product).unwrap();
    assert!(input.branches.as_ref().unwrap().is_product());
}
