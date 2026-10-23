use xlc_format::library::Library;
use xlc_format::table::Axis;

use super::*;
use crate::depletion::tests::{PIN, U235, deplete_with, library};
use crate::depletion::{History, transport_name};
use crate::input::{LatticeInput, MaterialInput, SymmetryChoice};
use crate::material::Material;
use crate::moc::{Options, available_order, solve_symmetric};
use crate::shielding::shield_fuel_with;
use crate::symmetry::{Chosen, choose};
use crate::tracking::Tracking;

/// テストのライブラリに、ホウ酸水のホウ素を足したもの（B-10 は Xe-135、B-11 は O-16 のデータを借りる）
fn library_with_boron() -> Library {
    let mut library = library();
    let borrow = |from: &str, name: &str, awr: f64| {
        let mut data = library
            .nuclides
            .iter()
            .find(|n| n.name == from)
            .unwrap()
            .clone();
        data.name = name.into();
        data.awr = awr;
        data
    };
    let (b10, b11) = (
        borrow("Xe-135", "B-10", 9.9269),
        borrow("O-16", "B-11", 10.9147),
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

/// 分岐の状態の物質を、自己遮蔽と MOC で解いた k
fn k_of(
    text: &str,
    library: &Library,
    materials: &[MaterialInput],
    region_materials: &[String],
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
        ..Options::default()
    };
    solve_symmetric(&s.tracking, &regions, &options, s.chosen.moc())
        .unwrap()
        .k
}

fn history() -> (History, Library) {
    let library = library_with_boron();
    (deplete_with(&borated(PIN), &library), library)
}

#[test]
fn nominal_branches_reproduce_the_history_points() {
    let (history, library) = history();
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
        // 公称（1000 ppm）は輸送の組成も状態も履歴の点と同じ
        let nominal = k_of(
            &text,
            &library,
            &branches.state(2 * i + 1),
            &branches.region_materials,
        );
        let want = history.points[point].k;
        assert!(
            (nominal - want).abs() <= 1e-12 * want,
            "burnup {}: {nominal} vs {want}",
            history.points[point].burnup
        );
        // ホウ素を減らした分岐は組成が同じで、吸収が減る
        let diluted = k_of(
            &text,
            &library,
            &branches.state(2 * i),
            &branches.region_materials,
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
