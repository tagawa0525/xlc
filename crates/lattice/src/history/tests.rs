use xlc_format::library::{Library, Resolved, Shielded};
use xlc_format::table::Axis;

use super::*;
use crate::branches::{BranchShielder, remix_varied};
use crate::cmfd::CmfdOptions;
use crate::depletion::tests::{PIN, U235, deplete_with, deplete_with_moc, library};
use crate::depletion::{History, transport_name};
use crate::input::{LatticeInput, MaterialInput, SymmetryChoice};
use crate::material::Material;
use crate::moc::{Acceleration, Options, available_order, solve_symmetric};
use crate::shielding::shield_fuel_with;
use crate::symmetry::{Chosen, choose};
use crate::tracking::Tracking;

/// テストのライブラリに、ホウ酸水のホウ素を足したもの（B-10 は Xe-135、B-11 は O-16 のデータを借りる）
fn library_with_boron() -> Library {
    let mut library = library();
    let borrow = |from: &str, name: &str, za: f64, awr: f64| {
        let mut data = library
            .nuclides
            .iter()
            .find(|n| n.name == from)
            .unwrap()
            .clone();
        data.name = name.into();
        data.za = za;
        data.awr = awr;
        data
    };
    let (b10, b11) = (
        borrow("Xe-135", "B-10", 5010.0, 9.9269),
        borrow("O-16", "B-11", 5011.0, 10.9147),
    );
    library.nuclides.push(b10);
    library.nuclides.push(b11);
    library
}

/// 減速材をホウ酸水（1000 ppm）にしたテストのピンセル
fn borated(text: &str) -> String {
    text.replace(
        "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
        "water = { density = 0.7, boron = 1000, b10 = 0.2 }",
    )
}

/// 燃焼の入力 `PIN` の燃焼を外し、少数群定数のファイル、3 つの軸の分岐、`[history]` を足した入力
fn branch_input(burnup: &str, fuel: &str) -> String {
    let start = PIN.find("[depletion]").unwrap();
    let lattice = borated(&PIN[..start]).replace(
        "library = \"unused.msgpack\"",
        "library = \"unused.msgpack\"\nfew-group = { output = \"unused.msgpack\", bounds = [1e7, 0.625, 1e-5] }",
    );
    format!(
        "{lattice}\n[branches]\nboron = {{ ppm = [500, 1000], materials = [\"water\"] }}\n\
         moderator = {{ temperature = [600], density = [0.7], materials = [\"water\"] }}\n\
         fuel = {{ temperature = {fuel}, materials = [\"fuel\"] }}\n\
         [history]\nfile = \"unused.msgpack\"\nburnup = {burnup}\n"
    )
}

/// 解く準備（入力、トラック、領域ごとの物質の名前、対称操作）
struct Setup {
    input: LatticeInput,
    tracking: Tracking,
    region_materials: Vec<String>,
    chosen: Chosen,
}

fn setup(text: &str) -> Setup {
    let input = LatticeInput::parse(text).unwrap();
    let (lattice, region_materials) = input.geometry().unwrap();
    let tracking = Tracking::new(&lattice, input.azimuthal, input.spacing).unwrap();
    let forced = match input.symmetry {
        SymmetryChoice::Detect => None,
        SymmetryChoice::Force(s) => Some(s),
    };
    let chosen = choose(&tracking, &region_materials, forced).unwrap();
    Setup {
        input,
        tracking,
        region_materials,
        chosen,
    }
}

fn branches_of(
    text: &str,
    history: &History,
    library: &Library,
) -> anyhow::Result<HistoryBranches> {
    let s = setup(text);
    history_branches(
        &s.input,
        history,
        library,
        &s.tracking.lattice,
        &s.region_materials,
        &s.chosen,
    )
}

/// 分岐の状態の物質を、自己遮蔽と MOC（`moc` の打ち切り条件と加速）で解いた k
fn k_of(
    text: &str,
    library: &Library,
    materials: &[MaterialInput],
    region_materials: &[String],
    moc: &Options,
) -> f64 {
    let s = setup(text);
    let shielded = shield_fuel_with(
        library,
        &s.tracking,
        materials,
        region_materials,
        s.chosen.shielding(),
        &s.input.shielding_options(),
    )
    .unwrap();
    let regions: Vec<&Material> = shielded
        .region_material
        .iter()
        .map(|&i| &shielded.materials[i])
        .collect();
    let options = Options {
        scattering_order: available_order(&regions),
        ..*moc
    };
    solve_symmetric(&s.tracking, &regions, &options, s.chosen.moc())
        .unwrap()
        .k
}

fn history() -> (History, Library) {
    let library = library_with_boron();
    (deplete_with(&borated(PIN), &library), library)
}

/// MOC の判定を既定（k 1e-7、中性子束 1e-6）より 4 桁細かくした条件。CMFD は既定と同じ比で MOC より細かく解く（ADR-0110）
fn tight() -> Options {
    Options {
        k_tolerance: 1e-11,
        flux_tolerance: 1e-10,
        acceleration: Acceleration::Cmfd(CmfdOptions {
            k_tolerance: 1e-14,
            flux_tolerance: 1e-12,
            ..CmfdOptions::default()
        }),
        ..Options::default()
    }
}

#[test]
fn nominal_branches_reproduce_the_history_points() {
    // 履歴の燃焼度 0 より後の輸送は前の輸送の解から、分岐は平らな初期値から反復を始める（ADR-0195）ので、k は収束の誤差の
    // 程度で違う。このピンセル（k 0.08 程度）は止めた解の誤差が判定より大きく（#0160）、燃焼度 0.5 と 2.0 の点の k の相対差は
    // 判定の k の約 100〜300 倍で、判定に比例して縮む（判定 1e-8 で 1.2e-6、1e-9 で 2.1e-7、1e-11 で 1.5e-9。2026-10-23 の手元の
    // 測定）。履歴も分岐も判定を 1e-11 にして解き、相対 1e-8（測った差の約 7 倍）で比べる
    let library = library_with_boron();
    let history = deplete_with_moc(&borated(PIN), &library, &tight());
    // 履歴の燃焼度は 0、0.05、0.5、2.0。分岐はそのうち 3 点で、燃焼度がいちばん遅く変わる
    let text = branch_input("[0, 0.5, 2.0]", "[600]");
    let branches = branches_of(&text, &history, &library).unwrap();
    assert_eq!(
        branches.axes[0],
        Axis::Burnup {
            points: vec![0.0, 0.5, 2.0]
        }
    );
    assert_eq!(branches.axes.len(), 4);
    assert_eq!(branches.len(), 3 * 2);
    for (i, point) in [0, 2, 3].into_iter().enumerate() {
        // 公称（1000 ppm）は輸送の組成も状態も履歴の点と同じ。燃焼度 0 の点はどちらも平らな初期値から解くので同じ反復をたどる
        let nominal = k_of(
            &text,
            &library,
            &branches.state(2 * i + 1),
            &branches.region_materials,
            &tight(),
        );
        let want = history.points[point].k;
        assert!(
            (nominal - want).abs() <= if point == 0 { 1e-12 } else { 1e-8 } * want,
            "burnup {}: {nominal} vs {want}",
            history.points[point].burnup
        );
        // ホウ素を減らした分岐は組成が同じで、吸収が減る
        let diluted = k_of(
            &text,
            &library,
            &branches.state(2 * i),
            &branches.region_materials,
            &tight(),
        );
        assert!(diluted > nominal, "{diluted} vs {nominal}");
    }
}

#[test]
fn branch_regions_take_the_history_densities() {
    let (history, library) = history();
    let text = branch_input("[2.0]", "[600]");
    let branches = branches_of(&text, &history, &library).unwrap();
    let region = &history.regions[0];
    let name = transport_name(&region.material, region.cell);
    assert!(branches.region_materials.contains(&name));
    assert!(!branches.region_materials.contains(&"fuel".to_string()));
    for state in (0..branches.len()).map(|i| branches.state(i)) {
        let fuel = state.iter().find(|m| m.name == name).unwrap();
        assert_eq!(fuel.temperature, 600.0);
        // 輸送できる核種だけを、履歴の最後の点の数密度で（負は 0）
        let names: Vec<&str> = fuel.composition.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(names, history.transportable);
        for (nuclide, density) in &fuel.composition {
            let i = history.nuclides.iter().position(|n| n == nuclide).unwrap();
            assert_eq!(*density, history.points[3].densities[0][i].max(0.0));
        }
        // 燃焼させた物質そのものは解かない
        assert!(state.iter().all(|m| m.name != "fuel"));
    }
}

#[test]
fn branch_burnups_are_points_of_the_history() {
    let (history, library) = history();
    let err = branches_of(&branch_input("[0, 1.0]", "[600]"), &history, &library).unwrap_err();
    assert!(
        format!("{err:#}").contains("burnup 1 MWd/kg is not a point of the history"),
        "{err:#}"
    );
}

#[test]
fn the_branch_lattice_must_be_the_lattice_of_the_history() {
    let (history, library) = history();
    let text = branch_input("[0.5]", "[600]");
    // 初めの組成が違う
    let fuel = text.replace("U-235 = 7.0e-4", "U-235 = 7.1e-4");
    assert!(branches_of(&fuel, &history, &library).is_err());
    // 燃料のピンの大きさが違う（燃焼の領域の体積が違う）
    let radius = text.replace("radius = 0.4096", "radius = 0.41");
    assert!(branches_of(&radius, &history, &library).is_err());
}

#[test]
fn branch_temperatures_need_data_for_every_transported_nuclide() {
    // テストのライブラリは 600 K だけを持つ（近い温度で代用しない）
    let (history, library) = history();
    let err = branches_of(&branch_input("[0.5]", "[600, 900]"), &history, &library).unwrap_err();
    assert!(format!("{err:#}").contains("900"), "{err:#}");
}

#[test]
fn input_materials_named_like_depletion_regions_are_errors() {
    let (history, library) = history();
    let cell = history.regions[0].cell;
    let text = branch_input("[0.5]", "[600]").replace(
        "[material.water]",
        &format!(
            "[material.\"fuel@{cell}\"]\ntemperature = 600\nnuclides = {{ O-16 = 2.5e-2 }}\n\n[material.water]"
        ),
    );
    let err = branches_of(&text, &history, &library).unwrap_err();
    assert!(format!("{err:#}").contains("clashes"), "{err:#}");
}

#[test]
fn negative_history_densities_beyond_the_cram_error_are_errors() {
    let (mut history, library) = history();
    history.points[2].densities[0][U235] = -1.0e-3;
    let err = branches_of(&branch_input("[0.5]", "[600]"), &history, &library).unwrap_err();
    assert!(format!("{err:#}").contains("CRAM"), "{err:#}");
}

#[test]
fn the_branch_library_must_transport_the_nuclides_of_the_history() {
    // 燃焼のライブラリになかったチェーンの核種（U-236）を持つライブラリは、燃焼の輸送と違う核種を輸送するので使えない
    let (history, library) = history();
    let mut other = library.clone();
    let mut u236 = other
        .nuclides
        .iter()
        .find(|n| n.name == "U-238")
        .unwrap()
        .clone();
    u236.name = "U-236".into();
    other.nuclides.push(u236);
    let err = branches_of(&branch_input("[0.5]", "[600]"), &history, &other).unwrap_err();
    assert!(format!("{err:#}").contains("library"), "{err:#}");
}

#[test]
fn history_branches_carry_the_initial_heavy_metal_mass() {
    let (history, library) = history();
    let branches = branches_of(&branch_input("[0, 2.0]", "[600]"), &history, &library).unwrap();
    // 燃焼度 0 の組成（U-235 7.0e-4、U-238 2.2e-2 /(barn cm)）の U の質量 (g/cm)。O-16 は数えない（ADR-0190、ADR-0173）
    let awr = |name: &str| {
        library
            .nuclides
            .iter()
            .find(|n| n.name == name)
            .unwrap()
            .awr
    };
    let volume: f64 = history.regions.iter().map(|d| d.volume).sum();
    let grams_per_atom = |name: &str| awr(name) * 1.008_664_915_95 * 1.660_539_066_60e-24;
    let want =
        volume * 1.0e24 * (7.0e-4 * grams_per_atom("U-235") + 2.2e-2 * grams_per_atom("U-238"));
    assert!(
        (branches.heavy_metal - want).abs() <= 1e-13 * want,
        "{} vs {want}",
        branches.heavy_metal
    );
}

/// `library_with_boron()` の U-238 の熱群の捕獲を σ0 で変え（燃料の自己遮蔽が、Dancoff 係数を通して減速材の全断面積に依る）、
/// 分離共鳴を持ち群 0 が σ0 に依る Zr-90 を足し、すべての核種に 650 K のデータ（U-238 の群 0 の捕獲と、Zr-90 と H-1 の群 0 の
/// 全断面積を増やしたもの）を足したもの
fn library_for_reuse() -> Library {
    let mut library = library_with_boron();
    let mut zr = library
        .nuclides
        .iter()
        .find(|n| n.name == "O-16")
        .unwrap()
        .clone();
    zr.name = "Zr-90".into();
    zr.za = 40090.0;
    zr.awr = 89.132;
    zr.resolved = Some(Resolved {
        top: 1.0e4,
        potential: 6.0,
    });
    zr.temperatures[0].shielded = [15.0, 10.0, 5.0]
        .iter()
        .map(|&t| Shielded {
            total: vec![t, 8.0],
            elastic: vec![6.0, 4.0],
            fission: vec![0.0, 0.0],
            nu_fission: vec![0.0, 0.0],
            capture: vec![t - 6.0, 4.0],
        })
        .collect();
    library.nuclides.push(zr);
    for nuclide in &mut library.nuclides {
        if nuclide.name == "U-238" {
            for (table, factor) in nuclide.temperatures[0]
                .shielded
                .iter_mut()
                .zip([1.0, 0.8, 0.5])
            {
                let capture = table.capture[1] * factor;
                table.total[1] += capture - table.capture[1];
                table.capture[1] = capture;
            }
        }
        let mut hot = nuclide.temperatures[0].clone();
        hot.temperature = 650.0;
        let added = match nuclide.name.as_str() {
            "U-238" => 0.05,
            "Zr-90" => 0.5,
            "H-1" => 0.2,
            _ => 0.0,
        };
        for table in &mut hot.shielded {
            table.total[0] += added;
            if nuclide.name == "H-1" {
                table.elastic[0] += added;
            } else {
                table.capture[0] += added;
            }
        }
        nuclide.temperatures.push(hot);
    }
    library
}

/// 入力 `text` の `[branches]` に、状態ごとに自己遮蔽を解く指定（ADR-0201）を足したもの
fn every_state(text: &str) -> String {
    assert!(text.contains("[history]\n"));
    text.replace("[history]\n", "shielding = \"every-state\"\n[history]\n")
}

/// 入力 `text` のピンセルに、Zr-90 の支持格子の枠（物質 zirc。燃料のピンの等価なセルの外の塊として遮蔽する）を足したもの
fn framed(text: &str) -> String {
    assert!(text.contains("[cell.pin]\n"));
    text.replace(
        "[cell.pin]\n",
        "[material.zirc]\ntemperature = 600\nnuclides = { Zr-90 = 4.0e-2 }\n\n\
         [cell.pin]\nframe = { thickness = 0.05, material = \"zirc\" }\n",
    )
}

/// 減速材の軸を 2 点（600 K と 0.7 g/cc、650 K と 0.66 g/cc）にした `branch_input`
fn two_moderators(burnup: &str) -> String {
    let text = branch_input(burnup, "[600]");
    let axis = "moderator = { temperature = [600], density = [0.7], materials = [\"water\"] }";
    assert!(text.contains(axis));
    text.replace(
        axis,
        "moderator = { temperature = [600, 650], density = [0.7, 0.66], materials = [\"water\"] }",
    )
}

/// 状態の物質の一覧と、領域ごとのその番号
type Shielding = (Vec<Material>, Vec<usize>);

/// 入力 `text` の履歴からの分岐の状態を順に `BranchShielder` で自己遮蔽した物質と、自己遮蔽を解いた回数
fn shield_states(
    text: &str,
    history: &History,
    library: &Library,
) -> (HistoryBranches, Vec<Shielding>, usize) {
    let s = setup(text);
    let branches = history_branches(
        &s.input,
        history,
        library,
        &s.tracking.lattice,
        &s.region_materials,
        &s.chosen,
    )
    .unwrap();
    let mut shielder = BranchShielder::new(
        &s.tracking,
        &branches.region_materials,
        s.input.branches.as_ref().unwrap().varied(),
        s.chosen.shielding(),
        s.input.shielding_options(),
    );
    let shielded = (0..branches.len())
        .filter(|&i| branches.solved(i))
        .map(|i| {
            let reference = branches.reference(i);
            owned(
                shielder
                    .shield_state(library, &branches.state(i), reference, || {
                        let (point, r) = reference.unwrap();
                        branches.reference_state(point, r)
                    })
                    .unwrap(),
            )
        })
        .collect();
    let solved = shielder.shieldings();
    (branches, shielded, solved)
}

/// 領域 `region` の物質
/// 分岐の自己遮蔽の物質の一覧と [領域] のその番号
fn owned(state: crate::branches::StateShielding) -> Shielding {
    (state.materials().to_vec(), state.region_material().to_vec())
}

fn region_material(shielding: &Shielding, region: usize) -> &Material {
    &shielding.0[shielding.1[region]]
}

/// 物質の名前が `name` に合う領域
fn regions_of(names: &[String], name: impl Fn(&str) -> bool) -> Vec<usize> {
    names
        .iter()
        .enumerate()
        .filter(|(_, n)| name(n))
        .map(|(r, _)| r)
        .collect()
}

#[test]
fn nominal_boron_shielding_solves_once_per_burnup_moderator_and_fuel_point() {
    let library = library_for_reuse();
    let history = deplete_with(&borated(PIN), &library);
    let text = two_moderators("[0, 0.5]");
    let (branches, solved, solves) = shield_states(&every_state(&text), &history, &library);
    let (_, reused, reuses) = shield_states(&text, &history, &library);
    // 燃焼度 2 点 × ホウ素 2 点 × 減速材 2 点 × 燃料 1 点を、状態ごとに解く
    assert_eq!(branches.len(), 8);
    assert_eq!(solves, 8);
    // 既定は燃焼度、減速材、燃料の点ごとに 1 回（公称の 1000 ppm は軸の点）
    assert_eq!(reuses, 4);
    let fuel = regions_of(&branches.region_materials, |n| n.contains('@'));
    let water = regions_of(&branches.region_materials, |n| n == "water");
    assert!(!fuel.is_empty() && !water.is_empty());
    for point in 0..2 {
        for moderator in 0..2 {
            let at = |boron: usize| (point * 2 + boron) * 2 + moderator;
            // 公称のホウ素の状態は、その減速材で解いた自己遮蔽そのもの
            assert_eq!(reused[at(1)], solved[at(1)]);
            // ホウ素を変えた状態の燃料は公称のホウ素の状態の物質のまま。解き直すと Dancoff 係数を通して変わる
            for &r in &fuel {
                assert_eq!(
                    region_material(&reused[at(0)], r),
                    region_material(&solved[at(1)], r)
                );
                assert_ne!(
                    region_material(&reused[at(0)], r),
                    region_material(&solved[at(0)], r)
                );
            }
            // 減速材（無限希釈）は状態の組成で作る
            for &r in &water {
                assert_eq!(
                    region_material(&reused[at(0)], r),
                    region_material(&solved[at(0)], r)
                );
            }
        }
        // 減速材を変えた状態は解き直す（燃料の物質が違う）
        for &r in &fuel {
            assert_ne!(
                region_material(&reused[point * 4], r),
                region_material(&reused[point * 4 + 1], r)
            );
        }
    }
}

#[test]
fn branches_without_a_history_shield_at_the_nominal_boron_too() {
    let library = library_for_reuse();
    let history_text = two_moderators("[0]");
    let text = &history_text[..history_text.find("[history]").unwrap()];
    let s = setup(text);
    let branches = s.input.branches.as_ref().unwrap();
    let states = s.input.branch_states();
    let mut shielder = BranchShielder::new(
        &s.tracking,
        &s.region_materials,
        branches.varied(),
        s.chosen.shielding(),
        s.input.shielding_options(),
    );
    let shielded: Vec<Shielding> = states
        .iter()
        .enumerate()
        .map(|(i, state)| {
            let r = branches.reference(i).unwrap();
            owned(
                shielder
                    .shield_state(&library, state, Some((0, r)), || {
                        branches.references()[r].clone()
                    })
                    .unwrap(),
            )
        })
        .collect();
    // ホウ素 2 点 × 減速材 2 点で、減速材の点ごとに 1 回
    assert_eq!(states.len(), 4);
    assert_eq!(shielder.shieldings(), 2);
    for &r in &regions_of(&s.region_materials, |n| n == "fuel") {
        assert_eq!(
            region_material(&shielded[0], r),
            region_material(&shielded[2], r)
        );
        assert_ne!(
            region_material(&shielded[0], r),
            region_material(&shielded[1], r)
        );
    }
}

/// 枠（Zr-90 の塊）のある格子の、燃焼度 0.5 MWd/kg の分岐（減速材 2 点）と、その公称の状態（状態 2）で解いた自己遮蔽
fn framed_nominal(
    library: &Library,
) -> (
    Setup,
    HistoryBranches,
    Vec<MaterialInput>,
    crate::shielding::FuelShielding,
) {
    let history = deplete_with(&framed(&borated(PIN)), library);
    let s = setup(&framed(&two_moderators("[0.5]")));
    let branches = history_branches(
        &s.input,
        &history,
        library,
        &s.tracking.lattice,
        &s.region_materials,
        &s.chosen,
    )
    .unwrap();
    let nominal = branches.state(2);
    let shielded = shield_fuel_with(
        library,
        &s.tracking,
        &nominal,
        &branches.region_materials,
        s.chosen.shielding(),
        &s.input.shielding_options(),
    )
    .unwrap();
    (s, branches, nominal, shielded)
}

#[test]
fn remixing_the_nominal_composition_reproduces_the_shielded_materials() {
    // 作り直しの前提（自己遮蔽の物質が、領域の実効断面積と組成から作った群定数に SPH 因子を掛けたもの）の確かめ方: 公称の組成で
    // 作り直すと、自己遮蔽の物質とビットで同じ。無限希釈の水と、塊（枠）の両方で
    let library = library_for_reuse();
    let (s, branches, nominal, shielded) = framed_nominal(&library);
    let rebuilt = remix_varied(
        &library,
        &shielded,
        &s.tracking.lattice,
        &branches.region_materials,
        &["water".to_string(), "zirc".to_string()],
        &nominal,
        &nominal,
    )
    .unwrap();
    assert_eq!(rebuilt.len(), 2);
    for (m, material) in rebuilt {
        assert_eq!(material, shielded.materials[m], "material {m}");
    }
}

#[test]
fn remixed_materials_take_the_state_densities_at_the_nominal_temperature() {
    let library = library_for_reuse();
    let (s, branches, nominal, shielded) = framed_nominal(&library);
    let varied = ["zirc".to_string()];
    let frame = regions_of(&branches.region_materials, |n| n == "zirc")[0];
    // 枠の数密度を 2 倍にした状態は、同じ実効断面積で、巨視的断面積が 2 倍
    let mut state = nominal.clone();
    let zirc = state.iter_mut().find(|m| m.name == "zirc").unwrap();
    zirc.composition[0].1 *= 2.0;
    let rebuilt = remix_varied(
        &library,
        &shielded,
        &s.tracking.lattice,
        &branches.region_materials,
        &varied,
        &nominal,
        &state,
    )
    .unwrap();
    let m = shielded.region_material[frame];
    let (_, material) = rebuilt.iter().find(|(i, _)| *i == m).unwrap();
    for (a, b) in material.total.iter().zip(&shielded.materials[m].total) {
        assert!((a - 2.0 * b).abs() <= 1e-12 * a, "{a} vs {b}");
    }
    // 温度を変えた状態は、公称の温度の実効断面積から作れない
    let mut hot = nominal.clone();
    hot.iter_mut()
        .find(|m| m.name == "zirc")
        .unwrap()
        .temperature = 650.0;
    assert!(
        remix_varied(
            &library,
            &shielded,
            &s.tracking.lattice,
            &branches.region_materials,
            &varied,
            &nominal,
            &hot,
        )
        .is_err()
    );
}

#[test]
fn a_reference_that_differs_beyond_the_boron_axis_is_an_error() {
    // 参照から作るのはホウ素の軸の物質だけなので、ほかの物質が参照と違う状態（別の減速材の点の参照を渡した）は止める
    let library = library_for_reuse();
    let history = deplete_with(&borated(PIN), &library);
    let s = setup(&branch_input("[0.5]", "[600, 650]"));
    let branches = history_branches(
        &s.input,
        &history,
        &library,
        &s.tracking.lattice,
        &s.region_materials,
        &s.chosen,
    )
    .unwrap();
    let varied = s.input.branches.as_ref().unwrap().varied();
    let mut shielder = BranchShielder::new(
        &s.tracking,
        &branches.region_materials,
        varied,
        s.chosen.shielding(),
        s.input.shielding_options(),
    );
    // 状態 1 は燃料 650 K、参照 0 は燃料 600 K
    let err = shielder
        .shield_state(&library, &branches.state(1), Some((0, 0)), || {
            branches.reference_state(0, 0)
        })
        .unwrap_err();
    assert!(format!("{err:#}").contains("fuel@"), "{err:#}");
}

#[test]
fn off_nominal_states_are_solved_at_their_burnups_only() {
    // 燃焼度 0、0.5、2.0 の点で、公称の状態（ホウ素 1000 ppm、減速材 600 K、燃料 600 K）はすべての点、ほかの状態は 0 と 2.0 だけ
    // （ADR-0200 の案 B2）。ホウ素 2 点 × 減速材 2 点 × 燃料 1 点は 2 軸の組で直積
    let (history, library) = history();
    let text = format!(
        "{}off-nominal-burnup = [0, 2.0]\n",
        two_moderators("[0, 0.5, 2.0]")
    );
    let branches = branches_of(&text, &history, &library).unwrap();
    assert_eq!(branches.len(), 3 * 4);
    let solved: Vec<usize> = (0..12).filter(|&i| branches.solved(i)).collect();
    // 状態の番号は燃焼度 × 4 + (ホウ素 × 2 + 減速材)。公称は (1, 0)
    assert_eq!(solved, vec![0, 1, 2, 3, 6, 8, 9, 10, 11]);
    // 公称の状態が格子点にないと、公称でない状態の燃焼点を間引けない
    let off = text.replace(
        "boron = { ppm = [500, 1000], materials = [\"water\"] }",
        "boron = { ppm = [500, 700], materials = [\"water\"] }",
    );
    let err = format!("{:#}", branches_of(&off, &history, &library).unwrap_err());
    assert!(err.contains("nominal"), "{err}");
}

#[test]
fn pair_subgrid_states_shield_like_the_product_states() {
    // 2 軸の組の部分格子の状態を公称のホウ素の自己遮蔽（ADR-0201）で解くと、直積の同じ状態とビットで同じ。参照は、解く状態の
    // 燃焼度、減速材、燃料の点ごとに 1 度
    let library = library_for_reuse();
    let history = deplete_with(&borated(PIN), &library);
    let text = two_moderators("[0, 0.5]").replace(
        "fuel = { temperature = [600], materials = [\"fuel\"] }",
        "fuel = { temperature = [600, 650], materials = [\"fuel\"] }",
    );
    let product = text.replace(
        "[history]",
        "subgrids = [[\"boron\", \"moderator\", \"fuel\"]]\n[history]",
    );
    let (branches, pairs, pair_solves) = shield_states(&text, &history, &library);
    let (_, all, all_solves) = shield_states(&product, &history, &library);
    let solved: Vec<usize> = (0..branches.len())
        .filter(|&i| branches.solved(i))
        .collect();
    // 2 × 2 × 2 の 8 状態から、3 つの軸がどれも公称でない状態（ホウ素 500 ppm、650 K、650 K）を除く
    assert_eq!(solved.len(), 2 * 7);
    for (n, &i) in solved.iter().enumerate() {
        assert_eq!(pairs[n], all[i], "state {i}");
    }
    assert_eq!(pair_solves, 2 * 4);
    assert_eq!(all_solves, 2 * 4);
}

/// 分岐の状態 `index` を自己遮蔽と MOC（`tight`）で解いた、自己遮蔽の結果と MOC の解
fn solved(
    text: &str,
    library: &Library,
    branches: &HistoryBranches,
    index: usize,
) -> (Setup, crate::shielding::FuelShielding, crate::moc::Solution) {
    let s = setup(text);
    let shielded = shield_fuel_with(
        library,
        &s.tracking,
        &branches.state(index),
        &branches.region_materials,
        s.chosen.shielding(),
        &s.input.shielding_options(),
    )
    .unwrap();
    let regions: Vec<&Material> = shielded
        .region_material
        .iter()
        .map(|&i| &shielded.materials[i])
        .collect();
    let options = Options {
        scattering_order: available_order(&regions),
        ..tight()
    };
    let solution = solve_symmetric(&s.tracking, &regions, &options, s.chosen.moc()).unwrap();
    (s, shielded, solution)
}

/// 少数群の境界（テストのライブラリの群と同じ 2 群）
const BOUNDS: [f64; 3] = [1.0e7, 0.625, 1.0e-5];

#[test]
fn micro_densities_are_volume_averages_of_the_history_over_the_lattice() {
    let (history, library) = history();
    let text = format!(
        "{}micro = [\"U-235\", \"Xe-135\", \"U-238\"]\n",
        branch_input("[0, 0.5, 2.0]", "[600]")
    );
    let branches = branches_of(&text, &history, &library).unwrap();
    assert_eq!(branches.micro_nuclides(), ["U-235", "Xe-135", "U-238"]);
    // 格子の体積はピッチ 1.26 cm の正方形
    let lattice = 1.26 * 1.26;
    for (k, point) in [0, 2, 3].into_iter().enumerate() {
        for (j, name) in ["U-235", "Xe-135", "U-238"].into_iter().enumerate() {
            let i = history.nuclides.iter().position(|n| n == name).unwrap();
            let want: f64 = history
                .regions
                .iter()
                .zip(&history.points[point].densities)
                .map(|(d, n)| d.volume * n[i].max(0.0))
                .sum::<f64>()
                / lattice;
            let got = branches.micro_densities()[k][j];
            assert!(
                (got - want).abs() <= 1e-14 * want.abs(),
                "{name} at point {point}: {got} vs {want}"
            );
        }
    }
    // 燃焼度 0 の Xe-135 は 0
    assert_eq!(branches.micro_densities()[0][1], 0.0);
}

/// 追う核種の微視的断面積 `micro` と均質化した数密度 `densities` の Σᵢ σᵢ N̄ᵢ が、少数群定数のうち燃焼させた物質（名前に `@` を含む
/// 領域）の分（MOC の中性子束と体積で、解いた物質 `materials`（[領域] の番号 `region_material`）を均質化したもの）と丸めの範囲で
/// 同じこと（ADR-0216）
fn assert_micro_closes(
    micro: &[xlc_format::table::MicroNuclide],
    densities: &[f64],
    region_names: &[String],
    materials: &[Material],
    region_material: &[usize],
    solution: &crate::moc::Solution,
    what: &str,
) {
    let fuel: Vec<usize> = regions_of(region_names, |n| n.contains('@'));
    let flux: Vec<f64> = (0..2)
        .map(|g| {
            (0..solution.volumes.len())
                .map(|r| solution.volumes[r] * solution.flux[r][g])
                .sum()
        })
        .collect();
    let part = |value: &dyn Fn(&Material) -> f64, g: usize| -> f64 {
        fuel.iter()
            .map(|&r| {
                solution.volumes[r] * solution.flux[r][g] * value(&materials[region_material[r]])
            })
            .sum::<f64>()
            / flux[g]
    };
    let sum = |value: &dyn Fn(&xlc_format::table::MicroNuclide) -> f64| -> f64 {
        micro.iter().zip(densities).map(|(m, n)| n * value(m)).sum()
    };
    let close = |got: f64, want: f64, item: &str| {
        assert!(
            (got - want).abs() <= 1e-12 * want.abs().max(1e-300),
            "{what} {item}: {got} vs {want}"
        );
    };
    for g in 0..2 {
        close(sum(&|m| m.total[g]), part(&|m| m.total[g], g), "total");
        close(
            sum(&|m| m.nu_fission[g]),
            part(&|m| m.nu_fission[g], g),
            "nu-fission",
        );
        close(
            sum(&|m| m.fission[g]),
            part(&|m| m.fission[g], g),
            "fission",
        );
        for to in 0..2 {
            close(
                sum(&|m| m.scatter[g][to]),
                part(&|m| m.scatter[g][to], g),
                "scatter",
            );
        }
    }
}

/// 輸送するすべての核種を `micro` に挙げた分岐の入力
fn all_transported(history: &History, burnup: &str) -> String {
    let names: Vec<String> = history
        .transportable
        .iter()
        .map(|n| format!("\"{n}\""))
        .collect();
    format!(
        "{}micro = [{}]\n",
        branch_input(burnup, "[600]"),
        names.join(", ")
    )
}

#[test]
fn micro_cross_sections_close_the_homogenized_rates_of_the_depleted_materials() {
    // 輸送するすべての核種を追えば、Σᵢ σᵢ N̄ᵢ は少数群定数のうち燃焼させた物質の分と同じ
    let (history, library) = history();
    let text = all_transported(&history, "[0, 2.0]");
    let branches = branches_of(&text, &history, &library).unwrap();
    for index in 0..branches.len() {
        let (s, shielded, solution) = solved(&text, &library, &branches, index);
        let micro = branches
            .micro(
                index,
                &library,
                &shielded,
                &s.tracking.lattice,
                &solution,
                &BOUNDS,
            )
            .unwrap();
        let point = index / (branches.len() / 2);
        assert_micro_closes(
            &micro,
            &branches.micro_densities()[point],
            &branches.region_materials,
            &shielded.materials,
            &shielded.region_material,
            &solution,
            &format!("state {index}"),
        );
    }
}

#[test]
fn micro_cross_sections_take_the_fuel_shielding_of_the_branch_shielder() {
    // 分岐の自己遮蔽（ADR-0201）が返す燃料の自己遮蔽は、状態ごとに解いたもの（every-state）でも、ホウ素の軸の物質だけを作り直す
    // 参照のもの（nominal-boron）でも、返した物質の燃料と同じ実効断面積と SPH 因子を持つ
    let (history, library) = history();
    for text in [
        all_transported(&history, "[0, 2.0]"),
        every_state(&all_transported(&history, "[0, 2.0]")),
    ] {
        let s = setup(&text);
        let branches = branches_of(&text, &history, &library).unwrap();
        let mut shielder = BranchShielder::new(
            &s.tracking,
            &branches.region_materials,
            s.input.branches.as_ref().unwrap().varied(),
            s.chosen.shielding(),
            s.input.shielding_options(),
        );
        for index in 0..branches.len() {
            let reference = branches.reference(index);
            let shielded = shielder
                .shield_state(&library, &branches.state(index), reference, || {
                    let (point, r) = reference.unwrap();
                    branches.reference_state(point, r)
                })
                .unwrap();
            let (materials, region_material, fuel) = (
                shielded.materials(),
                shielded.region_material(),
                shielded.fuel(),
            );
            let regions: Vec<&Material> = region_material.iter().map(|&i| &materials[i]).collect();
            let options = Options {
                scattering_order: available_order(&regions),
                ..tight()
            };
            let solution =
                solve_symmetric(&s.tracking, &regions, &options, s.chosen.moc()).unwrap();
            let micro = branches
                .micro(
                    index,
                    &library,
                    fuel,
                    &s.tracking.lattice,
                    &solution,
                    &BOUNDS,
                )
                .unwrap();
            let point = index / (branches.len() / 2);
            assert_micro_closes(
                &micro,
                &branches.micro_densities()[point],
                &branches.region_materials,
                materials,
                region_material,
                &solution,
                &format!("state {index} (reference {reference:?})"),
            );
        }
    }
}

#[test]
fn micro_reactions_follow_the_library_chain_reactions_and_split_fissions_at_the_yield_cutoff() {
    let (history, library) = history();
    let text = format!(
        "{}micro = [\"U-235\", \"Xe-135\"]\n",
        branch_input("[0, 2.0]", "[600]")
    );
    let branches = branches_of(&text, &history, &library).unwrap();
    // 燃焼度 0 の公称（1000 ppm）の状態。Xe-135 の数密度は 0 で、燃料の中で一様に重みを付ける
    let index = 1;
    let (s, shielded, solution) = solved(&text, &library, &branches, index);
    let micro = branches
        .micro(
            index,
            &library,
            &shielded,
            &s.tracking.lattice,
            &solution,
            &BOUNDS,
        )
        .unwrap();
    let (u235, xe135) = (&micro[0], &micro[1]);
    // テストのライブラリのチェーンの反応は、状態に分けない捕獲だけ
    for m in [u235, xe135] {
        assert_eq!(m.reactions.len(), 1);
        assert_eq!(m.reactions[0].mt, 102);
        assert_eq!(m.reactions[0].channels.len(), 1);
        // 経路が 1 つなら、捕獲の経路の断面積は吸収から核分裂を引いたもの（テストのライブラリは散乱、核分裂、捕獲だけ）
        for g in 0..2 {
            let absorption = m.total[g] - m.scatter[g].iter().sum::<f64>();
            let capture = m.reactions[0].channels[0][g];
            assert!(
                (capture - (absorption - m.fission[g])).abs() <= 1e-12 * capture,
                "{capture} vs {absorption} - {}",
                m.fission[g]
            );
        }
    }
    assert!(xe135.total.iter().all(|s| s.is_finite() && *s > 0.0));
    assert_eq!(xe135.thermal_fission, None);
    // 収率の境界 √(0.0253 · 5e5) eV は速い群（0.625 eV〜10 MeV）の中にあり、レサジーの幅の比で分ける（ADR-0183）。少数群の
    // 核分裂率は σf × 中性子束で、細かい群と少数群が同じなので、割合はその比で決まる
    let cutoff = (0.0253_f64 * 5.0e5).sqrt();
    let below = (cutoff / 0.625).ln() / (1.0e7 / 0.625_f64).ln();
    let flux: Vec<f64> = (0..2)
        .map(|g| {
            (0..solution.volumes.len())
                .map(|r| solution.volumes[r] * solution.flux[r][g])
                .sum()
        })
        .collect();
    let rates = [u235.fission[0] * flux[0], u235.fission[1] * flux[1]];
    let want = (rates[1] + below * rates[0]) / (rates[0] + rates[1]);
    let got = u235.thermal_fission.unwrap();
    assert!((got - want).abs() <= 1e-12, "{got} vs {want}");
}

#[test]
fn micro_nuclides_must_be_transported_in_the_history() {
    let (history, library) = history();
    // U-236 はチェーンにあるがライブラリにない（輸送しない）。Pu-239 は履歴の核種でない
    for name in ["U-236", "Pu-239"] {
        let text = format!(
            "{}micro = [\"U-235\", \"{name}\"]\n",
            branch_input("[0.5]", "[600]")
        );
        let err = branches_of(&text, &history, &library).unwrap_err();
        assert!(format!("{err:#}").contains(name), "{err:#}");
    }
}
