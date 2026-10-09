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
    let input = edited(
        "[branches]\nboron = { ppm = [900, 1700], materials = [\"water\"] }\n",
        "[branches]\n",
    )
    .unwrap();
    // 減速材と燃料の軸だけ
    assert_eq!(input.branch_states().len(), 6);
    let i = INPUT.find("[branches]").unwrap();
    let input = LatticeInput::parse(&INPUT[..i]).unwrap();
    assert!(input.branches.is_none());
    assert_eq!(input.branch_states(), vec![input.materials.clone()]);
}

#[test]
fn invalid_branches_are_errors() {
    // 知らない物質、ホウ酸水でない物質のホウ素、ホウ酸水のない減速材の軸
    assert!(
        edited(
            "ppm = [900, 1700], materials = [\"water\"]",
            "ppm = [900, 1700], materials = [\"oil\"]"
        )
        .is_err()
    );
    assert!(
        edited(
            "ppm = [900, 1700], materials = [\"water\"]",
            "ppm = [900, 1700], materials = [\"clad\"]"
        )
        .is_err()
    );
    assert!(
        edited(
            "materials = [\"water\", \"clad\"]",
            "materials = [\"clad\"]"
        )
        .is_err()
    );
    // 同じ物質が2つの温度の軸にある
    assert!(edited("materials = [\"fuel\"]", "materials = [\"fuel\", \"clad\"]").is_err());
    // 点が増えない、密度の数が違う、空の物質の一覧
    assert!(edited("temperature = [550, 580]", "temperature = [580, 550]").is_err());
    assert!(
        edited(
            "density = [0.76972, 0.74271, 0.71190]",
            "density = [0.76972, 0.74271]"
        )
        .is_err()
    );
    assert!(edited("materials = [\"fuel\"]", "materials = []").is_err());
    // ピン出力は分岐と一緒に出さない（状態ごとの形状関数は ADR-0096 でテーブルに入れる）
    assert!(
        edited(
            "rows = [\"pin\"]",
            "rows = [\"pin\"]\npin-powers = \"pins.msgpack\""
        )
        .is_err()
    );
}
